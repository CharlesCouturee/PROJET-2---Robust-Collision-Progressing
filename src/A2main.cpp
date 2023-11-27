/**
 * Provided code for particle system simulator.
 * This code provides the mouse interface for clicking and dragging particles, and the
 * code to draw the system.  When the simulator is running system.advanceTime is called
 * to numerically integrate the system forward.
 * @author kry
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"
#include "MatrixStack.h"
#include "Program.h"

#include "Text.hpp"

#include "ParticleSystem.hpp"
#include "TestSystems.hpp"
#include "SanityCheck.hpp"
#include "RecordVideo.hpp"

using namespace std;

ParticleSystem particleSystem;
SanityCheck sanityCheck;
TestSystems* testSystems;
bool run = false;
bool sanity = true;
float stepsize = 1.0f / 60.0f; // 60 Hz screen refresh
int substeps = 1;

// parameters for interacting with particles
float maxDist = 150;
float minDist = 50;
float grabThresh = 10;

// for grabbing particles
Particle* p1 = NULL;
Particle* p2 = NULL;
float d1 = 0;
float d2 = 0;

int xdown = 0;
int ydown = 0;
int xcurrent = 0;
int ycurrent = 0;
bool mouseDown = false;
bool mouseInWindow = false;
bool grabbed = false;
bool wasPinned = false;
bool closeToParticlePairLine = false;
bool canEdit = true; // when simulation hasn't been run/stepped yet

string scene = "";
int testSceneID = 0;
int showKeyBindings = 0;
// for openGL
GLFWwindow* window; // Main application window
string RES_DIR = ""; // Where data files live
shared_ptr<Program> progIM; // immediate mode

ImageRecorder imageRecorder;
bool recordFrames = false; // toggled by keyboard
bool recordFrame = false; // set to record when stepped in display

/** Finds the two closest particles for showing potential spring connections */
void findCloseParticles(int x, int y) {
	p1 = NULL;
	p2 = NULL;
	d1 = 0;
	d2 = 0;
	if (particleSystem.particles.size() > 0) {
		for (Particle* p : particleSystem.particles) {
			double d = p->distance(x, y);
			if (p1 == NULL || d < d1) {
				p2 = p1; d2 = d1; p1 = p; d1 = d;
			} else if (p2 == NULL || d < d2) {
				p2 = p; d2 = d;
			}
		}
	}
}

static void error_callback(int error, const char* description) {
	cerr << description << endl;
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if ((action != GLFW_PRESS) && (action != GLFW_REPEAT)) return;
	if (key == GLFW_KEY_ESCAPE) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	} else if (key == GLFW_KEY_SPACE) {
		run = !run;
		canEdit = false;
	} else if (key == GLFW_KEY_S) {
		for (int i = 0; i < substeps; i++) {
			testSystems->step(stepsize / substeps);
			particleSystem.advanceTime(stepsize / substeps);
		}
		canEdit = false;
	} else if (key == GLFW_KEY_R && !(mods & GLFW_MOD_SHIFT)) {
		if (particleSystem.name.compare("Alphabet Factory") == 0) {
			particleSystem.clearParticles();
			testSystems->createBox();
			testSystems->resetFactory();
		} else {
			particleSystem.resetParticles();
		}
		particleSystem.ccd.highestIterations = 0;
		sanity = true;
		if (!run) canEdit = true;
	} else if (key == GLFW_KEY_R && (mods & GLFW_MOD_SHIFT)) {
		recordFrames = !recordFrames;
	} else if (key == GLFW_KEY_C) {
		particleSystem.clearParticles();
		p1 = NULL;
		p2 = NULL;
	} else if (key == GLFW_KEY_D) {
		showKeyBindings = (showKeyBindings+1)%3;
	} else if (key == GLFW_KEY_T) {
		particleSystem.useGravity = !particleSystem.useGravity;
	} else if (key == GLFW_KEY_A) {
		sanity = true;
		particleSystem.ccd.highestIterations = 0;
		particleSystem.name = "Alphabet Factory";
		particleSystem.clearParticles();
		testSystems->resetFactory();
		particleSystem.useGravity = 1;
		testSystems->createBox();
		testSystems->runAlphabetFactory = true;
	} else if (key == GLFW_KEY_RIGHT) {
		sanity = true;
		testSceneID = testSystems->createSystem(++testSceneID);
		particleSystem.ccd.highestIterations = 0;
	} else if (key == GLFW_KEY_LEFT) {
		sanity = true;
		testSceneID = max(0, testSceneID - 1);
		testSceneID = testSystems->createSystem(testSceneID);
		particleSystem.ccd.highestIterations = 0;
	}
	if (key == GLFW_KEY_DELETE) {
		findCloseParticles(xcurrent, ycurrent);
		if (p1 != NULL && d1 < grabThresh) {
			particleSystem.remove(p1);
		}
	} else if (key == GLFW_KEY_Z) {
		for (Particle* p : particleSystem.particles) {
			p->v = glm::vec2(0, 0);
		}
	} else if (key == GLFW_KEY_UP) {
		substeps++;
	} else if (key == GLFW_KEY_DOWN) {
		substeps = max(substeps - 1, 1);
	}

	if (key == GLFW_KEY_F) {
		if ((mods & GLFW_MOD_SHIFT)) {
			int oldsubsteps = substeps;
			substeps = max(1, substeps / 2);
			stepsize *= (float)substeps / (float)oldsubsteps;
		} else {
			int oldsubsteps = substeps;
			substeps = min(128, substeps * 2);
			stepsize *= (float)substeps / (float)oldsubsteps;
		}
	}

	float scale = (mods & GLFW_MOD_SHIFT) ? 1.01f : 1 / 1.01f;
	if (key == GLFW_KEY_H) {
		stepsize *= scale;
	} else if (key == GLFW_KEY_V) {
		particleSystem.viscousDamping *= scale;
		if (particleSystem.viscousDamping == 0 && scale > 1) {
			particleSystem.viscousDamping = 0.01;
		} if (particleSystem.viscousDamping < 0.00999) {
			particleSystem.viscousDamping = 0;
		}
	} else if (key == GLFW_KEY_G) {
		particleSystem.gravity *= scale;
	} else if (key == GLFW_KEY_K) {
		particleSystem.springStiffness *= scale;
	} else if (key == GLFW_KEY_B) {
		particleSystem.springDamping *= scale;
		if (particleSystem.springDamping == 0 && scale > 1) {
			particleSystem.springDamping = 1;
		} else if (particleSystem.springDamping < 1) {
			particleSystem.springDamping = 0;
		}
	} else if (key == GLFW_KEY_O) {
		particleSystem.ccd.restitutionValue *= scale;
		if (particleSystem.ccd.restitutionValue > 1) {
			particleSystem.ccd.restitutionValue = 1;
		}
		if (particleSystem.ccd.restitutionValue == 0 && scale > 1) {
			particleSystem.ccd.restitutionValue = 0.01;
		} else if (particleSystem.ccd.restitutionValue < 0.01) {
			particleSystem.ccd.restitutionValue = 0;
		}
	} else if (key == GLFW_KEY_M) {
		particleSystem.ccd.thresholdDistance *= scale;
	} else if (key == GLFW_KEY_N) {
		particleSystem.ccd.repulsionStiffness *= scale;
	} else if (key == GLFW_KEY_I) {
		int increment = (mods & GLFW_MOD_SHIFT) ? -1 : 1;
		particleSystem.ccd.maxIterations = (int)(particleSystem.ccd.maxIterations + increment);
	} else if (key == GLFW_KEY_P) {
		particleSystem.ccd.useRepulsion = !particleSystem.ccd.useRepulsion;
	}
}

void mouse_pos_callback(GLFWwindow* window, double x, double y) {
	xcurrent = floor(x);
	ycurrent = floor(y);
	if (mouseDown) { // dragged
		particleSystem.mouseParticle.p.x = xcurrent;
		particleSystem.mouseParticle.p.y = ycurrent;
		if (grabbed) {
			p1->p = glm::vec2(xcurrent, ycurrent);
			p1->v = glm::vec2(0, 0);
			if (!run) {
				p1->p0 = p1->p;
				p1->v0 = p1->v;
				for (Spring* s : p1->springs) {
					s->recomputeRestLength();
				}
			}
		} else {
			findCloseParticles(xcurrent, ycurrent);
		}
	}
}

void mouse_callback(GLFWwindow* window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		double x;
		double y;
		glfwGetCursorPos(window, &x, &y);
		xcurrent = floor(x);
		ycurrent = floor(y);

		bool pressed = false;
		bool released = false;
		if (button == GLFW_MOUSE_BUTTON_LEFT) {
			if (GLFW_PRESS == action) {
				pressed = !mouseDown;
				mouseDown = true;
			} else if (GLFW_RELEASE == action) {
				released = mouseDown;
				mouseDown = false;
			}
		}

		if (!run) {
			if (pressed) {
				xdown = xcurrent;
				ydown = ycurrent;
				findCloseParticles(xcurrent, ycurrent);
				if (p1 != NULL && d1 < grabThresh) {
					wasPinned = p1->pinned;
					p1->pinned = true;
					grabbed = true;
					p1->p = glm::vec2(xcurrent, ycurrent);
					p1->v = glm::vec2(0, 0);
				}
			}

			if (released) {
				if (!grabbed && !run) {
					// were we within the threshold of a spring?
					if (closeToParticlePairLine) {
						if (!particleSystem.removeSpring(p1, p2) && canEdit) {
							particleSystem.createSpring(p1, p2);
						}
					} else if (canEdit) {
						Particle* p = particleSystem.createParticle(x, y, 0, 0);
						if (p1 != NULL && d1 < maxDist) {
							particleSystem.createSpring(p, p1);
						}
						if (p2 != NULL && d2 < maxDist) {
							particleSystem.createSpring(p, p2);
						}
					}
				} else if (grabbed && p1 != NULL) {
					p1->pinned = !wasPinned;
				}
				grabbed = false;
			}
		} else {
			if (pressed) {
				xdown = xcurrent;
				ydown = ycurrent;
				findCloseParticles(xcurrent, ycurrent);
				if (p1 != NULL && d1 < grabThresh) {
					particleSystem.mouseSpring.p1 = p1;
					particleSystem.useMouseSpring = true;
				}
				particleSystem.mouseParticle.p.x = xcurrent;
				particleSystem.mouseParticle.p.y = ycurrent;
			}
			if (released) {
				particleSystem.mouseSpring.p1 = NULL;
				particleSystem.useMouseSpring = false;
			}
		}
	}
}

static void init() {
	GLSL::checkVersion();

	// Check how many texture units are supported in the vertex shader
	int tmp;
	glGetIntegerv(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS, &tmp);
	cout << "GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS = " << tmp << endl;
	// Check how many uniforms are supported in the vertex shader
	glGetIntegerv(GL_MAX_VERTEX_UNIFORM_COMPONENTS, &tmp);
	cout << "GL_MAX_VERTEX_UNIFORM_COMPONENTS = " << tmp << endl;
	glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &tmp);
	cout << "GL_MAX_VERTEX_ATTRIBS = " << tmp << endl;

	// Set background color.
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	// Enable z-buffer test.
	glEnable(GL_DEPTH_TEST);

	// Initialize the GLSL programs.
	progIM = make_shared<Program>();
	progIM->setVerbose(true);
	progIM->setShaderNames(RES_DIR + "simple_vert.glsl", RES_DIR + "simple_frag.glsl");
	progIM->init();
	progIM->addUniform("P");
	progIM->addUniform("MV");
	progIM->setVerbose(false);

	initTextRender(RES_DIR);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glDisable(GL_DEPTH_TEST);
	particleSystem.init();
	testSystems = new TestSystems(&particleSystem, RES_DIR);
	testSystems->RES_DIR = RES_DIR;

	if (scene.empty()) {
		testSceneID = testSystems->createSystem(0);
	} else {
		particleSystem.name = scene;
		testSystems->parser.load(&particleSystem, RES_DIR + scene);
	}
	
	// If there were any OpenGL errors, this will print something.
	// You can intersperse this line in your code to find the exact location
	// of your OpenGL error.
	GLSL::checkError(GET_FILE_LINE);
}

/** draws a line from the given point to the given particle */
void drawLineToParticle(double x, double y, Particle* p, double d) {
	if (p == NULL) return;
	if (d > maxDist) return;
	double col = d < minDist ? 1 : (maxDist - d) / (maxDist - minDist);
	glColor4d(1 - col, 0, col, 0.75f);
	glBegin(GL_LINES);
	glVertex2d(x, y);
	glVertex2d(p->p.x, p->p.y);
	glEnd();
}

void display() {
	// set up projection for drawing in pixel units...

	if (run) {
		for (int i = 0; i < substeps; i++) {
			testSystems->step(stepsize / substeps);
			particleSystem.advanceTime(stepsize / substeps);
		}
		if ( recordFrames ) recordFrame = true;
	}

	if (!sanityCheck.sanityCheck(&particleSystem)) {
		run = false;
		sanity = false;
	}

	particleSystem.display();

	if (mouseDown) {
		if (!grabbed) {
			if (!run) {
				// check particle pair line
				if (p1 != NULL && p2 != NULL) {
					glm::vec2 v = p1->p - p2->p;
					v /= sqrt(v.x * v.x + v.y * v.y);
					double d = abs(v.x * (p1->p.y - ycurrent) - v.y * (p1->p.x - xcurrent));
					closeToParticlePairLine = d < grabThresh;
				}
				if (closeToParticlePairLine) {
					glColor4d(0, 1, 1, .5);
					glLineWidth(3.0f);
					glBegin(GL_LINES);
					glVertex2d(p1->p.x, p1->p.y);
					glVertex2d(p2->p.x, p2->p.y);
					glEnd();
				} else {
					glPointSize(5.0f);
					glLineWidth(2.0f);
					if (!run) {
						drawLineToParticle(xcurrent, ycurrent, p1, d1);
						drawLineToParticle(xcurrent, ycurrent, p2, d2);
					}
				}
			}
		} else {
			glPointSize(15.0f);
			glColor4d(0, 1, 0, 0.95);
			glBegin(GL_POINTS);
			glVertex2d(p1->p.x, p1->p.y);
			glEnd();
		}
	} else {
		//if ( mouseInWindow ) {
		findCloseParticles(xcurrent, ycurrent);
		if (p1 != NULL && d1 < grabThresh) {
			glPointSize(15.0f);
			glColor4d(0, 1, 0, 0.95);
			glBegin(GL_POINTS);
			glVertex2d(p1->p.x, p1->p.y);
			glEnd();
		}
		/*else if (p1 != NULL && p2 != NULL) {
			glm::vec2 v = p1->p - p2->p;
			v /= sqrt(v.x * v.x + v.y * v.y);
			double d = abs(v.x * (p1->p.y - ycurrent) - v.y * (p1->p.x - xcurrent));
			closeToParticlePairLine = d < grabThresh;
			if (closeToParticlePairLine) {
				glColor4d(0, 1, 1, .5);
				glLineWidth(3.0f);
				glBegin(GL_LINES);
				glVertex2d(p1->p.x, p1->p.y);
				glVertex2d(p2->p.x, p2->p.y);
				glEnd();
			}
		}*/
	}
}

static void render() {
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	float aspect = width / (float)height;
	glViewport(0, 0, width, height);
	particleSystem.width = width;
	particleSystem.height = height;

	// Clear framebuffer.
	if (sanity) {
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	} else {
		// Set background color.
		glClearColor(1.0f, 0.2f, 0.2f, 1.0f);
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Create matrix stacks.
	auto P = make_shared<MatrixStack>();
	auto MV = make_shared<MatrixStack>();

	// Setup projection for drawing in pixel coordinates
	float scale = 1;
	glm::mat4 projection = glm::ortho(0.0f, static_cast<float>(width * scale), static_cast<float>(height * scale), 0.0f);
	glm::mat4 modelview = glm::identity<glm::mat4>();

	progIM->bind();
	glUniformMatrix4fv(progIM->getUniform("P"), 1, GL_FALSE, &projection[0][0]);
	glUniformMatrix4fv(progIM->getUniform("MV"), 1, GL_FALSE, &modelview[0][0]);
	display();
	progIM->unbind();

	stringstream ss;
	ss << "\n";
	if (showKeyBindings == 0) {
		ss << "[left,right] - " << particleSystem.name << "\n";
		ss << "[h,H] h = " << stepsize << "\n";
		ss << "[v,V] c = " << particleSystem.viscousDamping << "\n";
		ss << "[b,B] b = " << particleSystem.springDamping << "\n";
		ss << "[k,K] k = " << particleSystem.springStiffness << "\n";
		ss << "[g,G] g = " << particleSystem.gravity << "\n";
		ss << "[t] use gravity = " << particleSystem.useGravity << "\n";
		ss << "[up,down] substeps = " << substeps << "\n";
		ss << "[i,I] ccd max iter = " << particleSystem.ccd.maxIterations << "\n";
		ss << "[p] use repulsion = " << particleSystem.ccd.useRepulsion << "\n";
		ss << "[m,M] repulsion threshold distance = " << particleSystem.ccd.thresholdDistance << "\n";
		ss << "[n,N] repulsion stiffness = " << particleSystem.ccd.repulsionStiffness << "\n";
		ss << "[o,O] restitution = " << particleSystem.ccd.restitutionValue << "\n";
		ss << "\n";
		ss << "computeTime = " << particleSystem.computeTime << "\n";
		ss << "CCD iterations = " << particleSystem.ccd.iterations << "\n";
		ss << "CCD highest iterations = " << particleSystem.ccd.highestIterations << "\n";
		ss << "\n";
		ss << "[R] recording" << (recordFrames?" ON ":" OFF ") << imageRecorder.filename << "\n";
		ss << "[d] show(hide) other keybindings" << "\n";
	} else if (showKeyBindings == 1) {
		ss << "[space] : start stop" << "\n";
		ss << "[s] step" << "\n";
		ss << "[r] reset" << "\n";
		ss << "[c] clear" << "\n";
		ss << "[z] zero velocities" << "\n";
		ss << "[a] alphabet factory scene" << "\n";
		ss << "[f,F] Faster(slower) h*=2 substeps*=2" << "\n";
	}
	string text = ss.str();
	RenderString(projection, modelview, 600, 50, 0.5, text);

	if (recordFrame) {
		recordFrame = false;
		imageRecorder.writeCurrentFrameToFile(window);
	}

	GLSL::checkError(GET_FILE_LINE);
}

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "Please specify the resource directory." << endl;
		return 0;
	} else if (argc < 3) {
		cout << "Unspecified scene to run, using default value " << scene << endl;
	} else {
		scene = string(argv[2]);
	}
	RES_DIR = argv[1] + string("/");

	// Set error callback.
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if (!glfwInit()) {
		return -1;
	}
	// https://en.wikipedia.org/wiki/OpenGL
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	// glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	// glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(1280, 720, "COMP 559 W23 - A2 CCD - CHARLES COUTURE", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}
	// Make the window's context current.
	glfwMakeContextCurrent(window);
	// Initialize GLEW.
	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}
	glGetError(); // A bug in glewInit() causes an error that we can safely ignore.
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	// Set vsync.
	glfwSwapInterval(1);
	// Set interaction callbacks.
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_callback);
	glfwSetCursorPosCallback(window, mouse_pos_callback);

	// Initialize scene.
	init();
	// Loop until the user closes the window.
	while (!glfwWindowShouldClose(window)) {
		// Render scene.
		render();
		// Swap front and back buffers.
		glfwSwapBuffers(window);
		// Poll for and process events.
		glfwPollEvents();
	}
	// Quit program.
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

