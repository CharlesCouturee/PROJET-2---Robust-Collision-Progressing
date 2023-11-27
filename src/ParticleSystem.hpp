#pragma once
#include <vector>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"

#include "Particle.hpp"
#include "Spring.hpp"
#include "BendingSpring.hpp"
#include "RobustCCD.hpp"

typedef glm::vec<2, double> vec2;

using namespace std;

/**
 * Implementation of a simple particle system
 * @author kry & Mercier-Aubin
 */
class ParticleSystem  {
    
public:
    string name = "";
    std::vector<Particle*> particles;
    std::vector<Spring*> springs;
    std::vector<BendingSpring*> bendingSprings;
    RobustCCD ccd;

    /** and extra spring we'll use to pull on selected particles */
    bool useMouseSpring = false;
    Particle mouseParticle;
    Spring mouseSpring = Spring(&mouseParticle,&mouseParticle);

    /**
     * Creates an empty particle system
     */
    ParticleSystem(){
        ccd = RobustCCD();
        mouseSpring.k = 10;
        mouseSpring.l0 = 0;
    }

    ~ParticleSystem() {
        for (int i = 1; i < springs.size(); ++i) {
            delete springs[i];
        }
        for (int i = 1; i < particles.size(); ++i) {
            delete particles[i];
        }
    }
    
    /**
     * Resets the positions of all particles
     */
    void resetParticles() {
        for ( Particle* p : particles ) {
            p->reset();
        }
        time = 0;
    }
    
    /**
     * Deletes all particles
     */
    void clearParticles() {
        for (Particle* p : particles) { delete p; }
        particles.clear();
        for (Spring* s : springs) { delete s; }
        springs.clear();
        for (BendingSpring* s : bendingSprings) { delete s; }
        bendingSprings.clear();
    }
    
    /** Elapsed simulation time */
    double time = 0;
        
    /**
     *  Compute nodal forces using gravity and spring contributions
     */
    void computeForce() {

        // reset forces and set gravity
        for (Particle* p : particles) {
            vec2 viscousDampingForce = -viscousDamping * p->mass * p->v;
            p->f = vec2(0, 1) * p->mass * (useGravity ? gravity : 0)  + viscousDampingForce;
        }

        //applies spring forces
        for (Spring* s : springs) {
            s->apply();
        }

        if (useMouseSpring) {
            mouseSpring.apply();
        }

        for (BendingSpring* s : bendingSprings) {
            s->apply();
        }
    }
    
    void stepVelocities(float h) {
        for (Particle* p : particles) {
            if (p->pinned) {
                p->f = vec2(0, 0); // just to make sure!
                p->v = vec2(0, 0);
            }
            else {
                p->v += (h / p->mass) * p->f;
            }
        }
    }

    /**
     * Update the positions with current velocities given the time step
     * @param h time step
     */
    void stepPositions(double h) {
        for (Particle* p : particles) {
            if (p->pinned) continue;
            p->p += h * p->v;
            p->f = vec2(0, 0);
        }
    }

    /** Time in seconds that was necessary to advance the system */
    float computeTime = 0;
    
    /**
     * Advances the state of the system
     * @param elapsed
     */
    void advanceTime( double elapsed ) {
        for (Spring* s : springs) {
            s->k = springStiffness;
            s->c = springDamping;
        }
        
        double now = glfwGetTime();
        
        // Simple symplectic Euler integration
        computeForce();
        stepVelocities(elapsed);
        ccd.applyRepulsion(elapsed, &particles, &springs);
        bool resultOK = ccd.checkCollision(elapsed, &particles, &springs);
        stepPositions(elapsed);
        
        time += elapsed;
        computeTime = (glfwGetTime() - now);
    }
    
    void filter(VectorXf& v) {
        for ( Particle* p : particles ) {
            if ( !p->pinned ) continue;
            v[ p->index*2+0] = 0;
            v[ p->index*2+1] = 0;
        }
    }

    /**
     * Creates a new particle and adds it to the system
     * @param x
     * @param y
     * @param vx
     * @param vy
     * @return the new particle
     */
    Particle* createParticle( float x, float y, float vx, float vy ) {
        Particle* p = new Particle( x, y, vx, vy );
        p->index = particles.size();
        particles.push_back( p );
        return p;
    }
    
    void remove( Particle* p ) {
    	for ( Spring* s : p->springs ) {
    		Particle* other = s->p1 == p ? s->p2 : s->p1; 
            // remove from the other particle's spring list
            vector<Spring*>::iterator osi = std::remove(other->springs.begin(), other->springs.end(), s);
            other->springs.erase(osi);
            // and remove from the master list
            vector<Spring*>::iterator si = std::remove(springs.begin(), springs.end(), s);
            springs.erase(si);
            delete s;
    	}

    	particles.erase( std::remove(particles.begin(),particles.end(), p ) );
    	// reset indices of each particle :(
    	for ( int i = 0 ; i < particles.size(); i++ ) {
    		particles[i]->index = i;
    	}
    }
    
    /**
     * Creates a new spring between two particles and adds it to the system.
     * @param p1
     * @param p2
     * @return the new spring
     */
    Spring* createSpring( Particle* p1, Particle* p2 ) {
        Spring* s = new Spring( p1, p2 ); 
        springs.push_back( s );         
        return s;
    }
    
    /**
     * Removes a spring between p1 and p2 if it exists, does nothing otherwise
     * @param p1
     * @param p2
     * @return true if the spring was found and removed
     */
    bool removeSpring( Particle* p1 = new Particle, Particle* p2 = new Particle ) {
    	Spring* found = NULL;
    	for ( Spring* s : springs ) {
    		if ( ( s->p1 == p1 && s->p2 == p2 ) || ( s->p1 == p2 && s->p2 == p1 ) ) {
    			found = s;
    			break;
    		}
    	}
    	if ( found != NULL ) {
            found->p1->springs.erase(std::remove(found->p1->springs.begin(), found->p1->springs.end(), found));
            found->p2->springs.erase(std::remove(found->p2->springs.begin(), found->p2->springs.end(), found));
            springs.erase(std::remove(springs.begin(),springs.end(), found));
			return true;
    	}
    	return false;
    }
    
    void init() {
        // do nothing
    }

    int height = 0;
    int width = 0;

    void display() {

        glPointSize( 4 );
        glBegin( GL_POINTS );
        for ( Particle* p : particles ) {
            double alpha = 0.75;
            if ( p->pinned ) {
                glColor4d( 1, 0, 0, alpha );
            } else {
                glColor4d( p->color.x, p->color.y, p->color.z, alpha );
            }
            glVertex2d( p->p.x, p->p.y );
        }
        glEnd();
        
        glColor4d(0,.5,.5,1);
        glLineWidth(2.0f);
        glBegin( GL_LINES );
        for (Spring* s : springs) {
            glVertex2d( s->p1->p.x, s->p1->p.y );
            glVertex2d( s->p2->p.x, s->p2->p.y );
        }
        if (useMouseSpring) {
            glColor4d(1, 0, 0, 1);
            Spring* s = &mouseSpring;
            glVertex2d(s->p1->p.x, s->p1->p.y);
            glVertex2d(s->p2->p.x, s->p2->p.y);
        }
        glEnd();
    }
    
    bool useGravity = true;
    double gravity = 9.8;
    double springStiffness = 100;
    double springDamping = 0;
    double viscousDamping = 0;
    double bendingStiffness = 1e3;
};