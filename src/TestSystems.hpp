#pragma once
#include "vector"
#include "Particle.hpp"
#include "Particlesystem.hpp"
#include "RobustCCD.hpp"
#include "Spring.hpp"
#include "AlphabetSoupFactory.hpp"
#include "XMLParser.hpp"

typedef glm::vec<2, double> vec2;

using namespace std;
# define M_PI           3.14159265358979323846

/**
 * Ugly class to procedurally help in creating a variety of interesting test scenes
 *
 * @author kry & Mercier-Aubin
 */
class TestSystems {

public:
	ParticleSystem* system;
	AlphabetSoupFactory alphabetSoupFactory;

	XMLParser parser = XMLParser();

	string RES_DIR;

	bool runAlphabetFactory = false;
	bool runPastaFactory = false;

	double stringParticleInterval = 3.0 / 60.0;
	double letterInterval = 500.0 / 60.0;
	double initialVelocity = 50;

	double elapsed = 0;
	double elapsedSince = letterInterval;
	int letterIndex = 0;

	TestSystems(ParticleSystem* system, string RES_DIR) : 
		alphabetSoupFactory(system, RES_DIR), 
		system(system) 
	{
		particles = &(system->particles);
		springs = &(system->springs);
		bendingSprings = &(system->bendingSprings);
	}

/**
 * reset the factory so that it always starts with the first letter
 */
	void resetFactory() {
		elapsed = 0;
		elapsedSince = letterInterval;
		letterIndex = 0;
		prevParticle0 = NULL;
		prevParticle1 = NULL;
	}

/**
 * Creates one of a number of simple test systems.
 *
 * Small systems are more useful for debugging!
 *
 * @param which
 * returns the ID of the system created (returns max defined if which is too large)
 */
	int createSystem(int which) {
		if (clearFirst) {
			system->clearParticles();
		}

		runAlphabetFactory = false;
		system->useGravity = true; // wanted for most of these examples
		system->ccd.useRepulsion = true;

		int ID = 0;
		if (which == ID++) {
			string scene = "interwovenHairsEasy.xml";
			system->name = std::to_string(which) + " - " + scene;
			parser.load(system, RES_DIR + scene);
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Simple";
			createSimple();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Pendulum Pair";
			createPendulumPair();
		} else if (which == ID++) {
			string scene = "pendulumPair.xml";
			system->name = std::to_string(which) + " - " + scene;
			parser.load(system, RES_DIR + scene);
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Bullet And Layers";
			createBulletAndLayers();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Chains";
			createChains();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - ZigZag";
			createZigZag();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Hairs Hard";
			createHairsHard();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Hairs Easy";
			createHairsEasy();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Triangle Truss";
			createTriangleTruss();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Chain Test";
			createChain(320, 100, 10, 200);
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Jacobi Test";
			createJacobiTest();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Rod Collision";
			createElasticRodCollision();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Bend Test";
			createSimpleBendTest();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Circles";
			createCircles();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Bent Hairs";
			createBentHairs();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - Sandwich";
			createSandwich();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - CHARLES COUTURE -260923463";
			// TODO: create your own here!
			createMyTest();
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - empty";
			// TODO: create your own here!
		} else if (which == ID++) {
			system->name = std::to_string(which) + " - empty";
			// TODO: create your own here!
		} else {
			system->name = std::to_string(which) + " - empty";
		}
		return ID - 1;
	}

	 /**
	  * Generates letters at regular intervals if the factory is turned on.
	  */
	void step(double stepSize) {
		if (!runAlphabetFactory && !runPastaFactory) return;

		elapsed += stepSize;
		elapsedSince += stepSize;

		if (runAlphabetFactory) {
			if (elapsedSince > letterInterval) {
				elapsedSince -= letterInterval;
				system->gravity = 9.8; // the fall is pretty big
				system->ccd.repulsionStiffness = 100;
				system->ccd.thresholdDistance = 2;
				system->ccd.restitutionValue = 0.1;
				int pos = letterIndex++ % alphabet.length();
				char letter = alphabet[pos];
				alphabetSoupFactory.createLetter(system, letter, 50, 200);
			}
		}

		if (runPastaFactory) {
			if (elapsedSince > stringParticleInterval) {
				elapsedSince -= stringParticleInterval;
				Particle* p = new Particle(100, 30, initialVelocity, sin(elapsed) * initialVelocity);
				particles->push_back(p);
				if (prevParticle1 != NULL) {
					Spring* s = new Spring(prevParticle1, p);
					s->recomputeRestLength(); // set rest length with current positions.
					springs->push_back(s);
				}
				if (prevParticle0 != NULL && prevParticle1 != NULL) {
					BendingSpring* bs = new BendingSpring(prevParticle0, prevParticle1, p);
					bs->theta0 = 0;
					//bs->kappa0 = 0;
					bendingSprings->push_back(bs);
				}
				prevParticle0 = prevParticle1;
				prevParticle1 = p;
			}
		}

	}
	
	void createMyTest()
	{
		Particle* p1 = new Particle(520, 100, 0, 0);
		Particle* p2 = new Particle(620, 100, 0, 0);
		p2->mass = 2;
		particles->push_back(p1);
		particles->push_back(p2);
		p1->pinned = true;
		springs->push_back(new Spring(p1, p2));

		Particle* p3 = new Particle(500, 100, 0, 0);
		Particle* p4 = new Particle(500, 200, 0, 0);
		particles->push_back(p3);
		particles->push_back(p4);
		p3->pinned = true;
		springs->push_back(new Spring(p3, p4));

		Particle* p5 = new Particle(480, 100, 0, 0);
		Particle* p6 = new Particle(480, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(460, 100, 0, 0);
		p6 = new Particle(460, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(440, 100, 0, 0);
		p6 = new Particle(440, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(420, 100, 0, 0);
		p6 = new Particle(420, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(400, 100, 0, 0);
		p6 = new Particle(400, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(380, 100, 0, 0);
		p6 = new Particle(380, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(360, 100, 0, 0);
		p6 = new Particle(360, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(340, 100, 0, 0);
		p6 = new Particle(340, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(320, 100, 0, 0);
		p6 = new Particle(320, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(300, 100, 0, 0);
		p6 = new Particle(300, 200, 0, 0);
		particles->push_back(p5);
		particles->push_back(p6);
		p5->pinned = true;
		springs->push_back(new Spring(p5, p6));

		p5 = new Particle(160, 81, 0, 0);
		p6 = new Particle(170, 150, 0, 0);
		p5->mass = 0.8;
		particles->push_back(p5);
		particles->push_back(p6);
		p6->pinned = true;
		springs->push_back(new Spring(p5, p6));

		Particle* p = new Particle(170, 40, 0, 0);
		particles->push_back(p);
		springs->push_back(new Spring(p, p6));

	}

/**
 * Creates a pinned box with an open top that matches the current display window
 */
	void createBox() {
		double h = system->height;
		double w = system->width;
		Particle* p1 = system->createParticle(5, 5, 0, 0); p1->pinned = true;
		Particle* p2 = system->createParticle(5, h - 5, 0, 0); p2->pinned = true;
		Particle* p3 = system->createParticle(w - 5, h - 5, 0, 0); p3->pinned = true;
		Particle* p4 = system->createParticle(w - 5, 5, 0, 0); p4->pinned = true;
		system->createSpring(p1, p2);
		system->createSpring(p2, p3);
		system->createSpring(p3, p4);
	}

private:
	vector<Particle*>* particles;
	vector<Spring*>* springs;
	vector<BendingSpring*>* bendingSprings;
	double mass = 1;
	bool clearFirst = true;

	/**
	 * Creates a free particle above a pinned line segment
	 */
	void createSimple() {
		Particle* p1 = new Particle(100, 400, 0, 0);
		Particle* p2 = new Particle(500, 400, 0, 0);
		particles->push_back(p1);
		particles->push_back(p2);
		p1->pinned = true;
		p2->pinned = true;
		springs->push_back(new Spring(p1, p2));
		Particle* p3 = new Particle(300, 50, 0, 0);
		particles->push_back(p3);
	}

	void createPendulumPair() {
		Particle* p1 = new Particle(320, 100, 0, 0);
		Particle* p2 = new Particle(520, 100, 0, 0);
		particles->push_back(p1);
		particles->push_back(p2);
		p1->pinned = true;
		springs->push_back(new Spring(p1, p2));
		Particle* p3 = new Particle(300, 150, 0, 0);
		Particle* p4 = new Particle(260, 150, 0, 0);
		particles->push_back(p3);
		particles->push_back(p4);
		p3->pinned = true;
		springs->push_back(new Spring(p3, p4));
	}

	/**
	 * Creates a vertical chain, haning from the specified location, with given nodes and length.
	 * Masses and set according to spring lengths.
	 * Does not include bend springs, but perhaps it should!
	 * @param xpos
	 * @param ypos
	 * @param N
	 * @param L
	 */
	void createChain(double xpos, double ypos, int N, double L) {
		double segLength = L / N;
		Particle* p1, * p2;
		p1 = new Particle(xpos, ypos, 0, 0);
		p1->mass = 0;
		p1->pinned = true;
		particles->push_back(p1);
		for (int i = 0; i <= N; i++) {
			ypos += i == 0 ? 20 : segLength;
			p2 = new Particle(xpos, ypos, 0, 0);
			p2->mass = 0;
			particles->push_back(p2);
			Spring* s = new Spring(p1, p2);
			springs->push_back(s);
			p1->mass += s->l0 / 20 / 2;
			p2->mass += s->l0 / 20 / 2;
			p1 = p2;
		}
	}

	/**
	 * Creates a row of chains with different numbers of elements.
	 * Would desire the same behaviour!
	 */
	void createChains() {
		int n = 5;
		for (int i = 0; i < 10; i++) {
			createChain(320 + i * 10, 100, n++, 200);
		}
	}

	void createTriangleTruss() {
		vec2 p = vec2(100, 100);
		vec2 d = vec2(20, 0);
		Particle* p1, * p2, * p3, * p4;
		p1 = new Particle(p.x - d.y, p.y + d.x, 0, 0);
		particles->push_back(p1);
		p2 = new Particle(p.x + d.y, p.y - d.x, 0, 0);
		particles->push_back(p2);
		springs->push_back(new Spring(p1, p2));
		p1->pinned = true;
		p2->pinned = true;
		p += d;
		p += d;
		int N = 10;
		for (int i = 1; i < N; i++) {
			// d.set( 20*Math.cos(i*Math.PI/N), 20*Math.sin(i*Math.PI/N) );
			p3 = new Particle(p.x - d.y, p.y + d.x, 0, 0);
			p4 = new Particle(p.x + d.y, p.y - d.x, 0, 0);
			particles->push_back(p3);
			particles->push_back(p4);
			springs->push_back(new Spring(p3, p1));
			//springs->push_back(new Spring(p3, p2));
			springs->push_back(new Spring(p4, p1));
			springs->push_back(new Spring(p4, p2));
			springs->push_back(new Spring(p4, p3));
			p1 = p3;
			p2 = p4;

			p += d;
			p += d;
		}
	}

	/**
	 * Creates a zig zag of particles with bending springs in between,
	 * a good test of non zero rest angle or curvature.
	 */
	void createZigZag() {
		int N = 20;
		int xpos = 100;
		Particle* p0, * p1, * p2;

		p0 = NULL;
		p1 = NULL;
		p2 = NULL;
		for (int i = 0; i < N; i++) {
			p2 = new Particle(xpos, 100 + 20 * (i % 2), 0, 0);
			particles->push_back(p2);
			if (i < 2) p2->pinned = true;
			if (p1 != NULL) springs->push_back(new Spring(p1, p2));
			if (p0 != NULL) bendingSprings->push_back(new BendingSpring(p0, p1, p2));
			p0 = p1;
			p1 = p2;
			xpos += 20;
		}
	}

	/**
	 * Like a chain, but creates with bending springs too!
	 * @param xpos
	 * @param ypos
	 * @param dx
	 * @param dy
	 */
	void createHair(double xpos, double ypos, double dx, double dy, int N) {
		Particle* p0 = NULL;
		Particle* p1 = NULL;
		Particle* p2 = NULL;
		for (int i = 0; i < N; i++) {
			p2 = new Particle(xpos + dx * i, ypos + dy * i, 0, 0);
			particles->push_back(p2);
			if (i < 2) p2->pinned = true;
			if (p1 != NULL) springs->push_back(new Spring(p1, p2));
			if (p0 != NULL) bendingSprings->push_back(new BendingSpring(p0, p1, p2));
			p0 = p1;
			p1 = p2;
		}
	}

	void createHairsHard() {
		int N = 20;
		int M = 5;
		int offset = 10;
		for (int j = 0; j < M; j++) {
			createHair(100, 100 + offset * 2 * j, 20, 0, 20);
			createHair(100 + (N - 2) * 20, 100 + offset * 2 * j + offset, -20, 0, 20);
		}
	}

	void createHairsEasy() {
		int N = 20;
		int M = 5;
		int offset = 15;
		for (int j = 0; j < M; j++) {
			createHair(100, 100 + offset * 2 * j, 20, 0, 20);
			createHair(100 + (N - 2) * 20, 100 + offset * 2 * j + offset, -20, 0, 20);
		}
	}

	void createBulletAndLayers() {
		int numHairs = 10;
		double spacing = 8;
		int N = 15;
		for (int i = 0; i < numHairs; i++) {
			createHair(300 + i * spacing, 100, 0, 10, N);
		}
		Particle* p = new Particle(200, 150, 20, 0);
		p->mass = 300;
		particles->push_back(p);
		system->useGravity = false;
		system->bendingStiffness = 1e3;
		// should change the base stiffness of the hairs instead!
	}

	/**
	 * Creates a test in which only symmetry preserved by Jacobi will lead to a
	 * balanced result.
	 */
	void createJacobiTest() {
		Particle* p1, * p2;
		// left end
		p1 = new Particle(250, 400, 0, 0);
		p1->pinned = true;
		particles->push_back(p1);
		p2 = new Particle(250, 250, 0, 0);
		particles->push_back(p2);
		springs->push_back(new Spring(p1, p2));

		// right end
		p1 = new Particle(500, 400, 0, 0);
		p1->pinned = true;
		particles->push_back(p1);
		p2 = new Particle(500, 250, 0, 0);
		particles->push_back(p2);
		springs->push_back(new Spring(p1, p2));

		// horizontal
		p1 = new Particle(250, 100, 0, 0);
		particles->push_back(p1);
		p2 = new Particle(500, 100, 0, 0);
		particles->push_back(p2);
		springs->push_back(new Spring(p1, p2));
	}

	void createElasticRodCollision() {
		// create a rod and set its velocity
		int N = 20;
		int xpos = 30;
		int ypos = 300;
		double dx = 100.0f / N;
		Particle* p1 = new Particle(xpos, ypos, 10, 0);
		particles->push_back(p1);
		xpos += dx;
		for (int i = 0; i < N; i++) {
			Particle* p2 = new Particle(xpos, ypos, 10, 0);
			particles->push_back(p2);
			springs->push_back(new Spring(p1, p2));
			p1 = p2;
			xpos += dx;
		}

		xpos = 300;
		p1 = new Particle(xpos, ypos, 0, 0);
		particles->push_back(p1);
		xpos += dx;
		for (int i = 0; i < N; i++) {
			Particle* p2 = new Particle(xpos, ypos, 0, 0);
			particles->push_back(p2);
			springs->push_back(new Spring(p1, p2));
			p1 = p2;
			xpos += dx;
		}
		system->useGravity = false;
	}

	void createSimpleBendTest() {
		int xpos = 300;
		int ypos = 300;
		double dx = 100;
		Particle* p0 = new Particle(xpos, ypos, 0, 0);
		p0->pinned = true;
		Particle* p1 = new Particle(xpos + dx, ypos, 0, 0);
		p1->pinned = true;
		Particle* p2 = new Particle(xpos, ypos + dx, 0, 0);
		particles->push_back(p0);
		particles->push_back(p1);
		particles->push_back(p2);
		Spring* s1 = new Spring(p0, p1);
		springs->push_back(s1);
		Spring* s2 = new Spring(p1, p2);
		springs->push_back(s2);
		BendingSpring* bs = new BendingSpring(p0, p1, p2);
		bendingSprings->push_back(bs);
	}

	/**
	 * Creates a particle with velocity between two
	 * springs, with one pinned, and the other pair being potentially
	 * quite massive.
	 * How many iterations of resolution will be necessary?
	 * Note the coefficient of restitution is set to 1 here!
	 */
	void createSandwich() {
		Particle* p1 = new Particle(100, 200, 0, 0);
		Particle* p2 = new Particle(100, 400, 0, 0);
		Particle* p3 = new Particle(100.1, 300, 200, 0);
		Particle* p4 = new Particle(100.2, 200, 0, 0);
		Particle* p5 = new Particle(100.2, 400, 0, 0);
		particles->push_back(p1);
		particles->push_back(p2);
		particles->push_back(p3);
		particles->push_back(p4);
		particles->push_back(p5);
		p1->pinned = true;
		p2->pinned = true;
		p3->mass = 2;
		p4->mass = mass;
		p5->mass = mass;
		springs->push_back(new Spring(p1, p2));
		springs->push_back(new Spring(p4, p5));
		system->ccd.restitutionValue;
		system->ccd.useRepulsion = false;
		system->useGravity = false;
	}
	/**
	 *
	 * @param cx
	 * @param cy
	 * @param r
	 * @param N
	 * @param k  create k of N points... close if k = N
	 * @param thetaOffset if zero starts at top
	 * @param pinned if true will do the first two
	 */
	void createCircle(double cx, double cy, double r, int N, int k, double thetaOffset, bool pinned) {
		Particle* p0 = NULL, * p1 = NULL, * p2 = NULL, * first = NULL, * second = NULL;
		for (int i = 0; i < k; i++) {
			double theta = i * 2 * M_PI / N;
			p2 = new Particle(cx - sin(theta + thetaOffset) * r, cy - cos(theta + thetaOffset) * r, 0, 0);
			p2->mass = 0;
			if (i == 0) first = p2;
			if (i == 1) second = p2;
			particles->push_back(p2);
			if (i < 2 && pinned) p2->pinned = true;
			if (p1 != NULL) {
				Spring* s = new Spring(p1, p2);
				springs->push_back(s);
				p1->mass += s->l0 / 20 / 2;
				p2->mass += s->l0 / 20 / 2;

			}
			if (p0 != NULL) bendingSprings->push_back(new BendingSpring(p0, p1, p2));
			p0 = p1;
			p1 = p2;
		}
		if (k == N) {
			Spring* s = new Spring(p1, first);
			springs->push_back(s);
			p1->mass += s->l0 / 20 / 2;
			first->mass += s->l0 / 20 / 2;
			bendingSprings->push_back(new BendingSpring(p0, p1, first));
			bendingSprings->push_back(new BendingSpring(p1, first, second));
		}
	}

	void createCircles() {
		createCircle(150, 150, 50, 20, 20, 0, false);
		createCircle(300, 150, 50, 30, 30, 0, false);
		createCircle(450, 150, 50, 40, 40, 0, false);
		createBox();
	}

	void createBentHairs() {
		int N = 20;
		double cx = 300;
		double cy = 200;
		double r = 50;
		vector<Particle*> circlePoints;
		for (int i = 0; i < N; i++) {
			double theta = i * 2 * M_PI / N;
			circlePoints.push_back( new Particle(cx - sin(theta) * r, cy - cos(theta) * r, 0, 0) );
			particles->push_back(circlePoints[i]);
		}

		int k = 10;
		for (int i = 0; i < N; i++) {
			springs->push_back(new Spring(circlePoints[i], circlePoints[(i + 1) % N]));
			bendingSprings->push_back(new BendingSpring(circlePoints[i], circlePoints[(i + 1) % N], circlePoints[(i + 2) % N]));
			double theta = i * 2 * M_PI / N;
			Particle* p0 = circlePoints[(i + N - 1) % N];
			Particle* p1 = circlePoints[i];
			for (int j = 1; j < k; j++) {
				double px = cx - sin(theta) * (r + j * 12);
				double py = cy - cos(theta) * (r + j * 12);
				Particle* p2 = new Particle(px, py, 0, 0);
				particles->push_back(p2);
				springs->push_back(new Spring(p1, p2));
				bendingSprings->push_back(new BendingSpring(p0, p1, p2));
				p0 = p1;
				p1 = p2;
			}
		}
		createBox();
	}

	Particle* prevParticle0 = NULL;
	Particle* prevParticle1 = NULL;

	const string alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";

};