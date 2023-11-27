#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

typedef glm::vec<2, double> vec2;

class Spring;
class BendingSpring;

/**
 * Particle class that contains particle properties (e.g., mass), 
 * initial positions and velocities, current position and velocities 
 * and a force accumulator for computing the total force acting on this particle.
 * @author kry
 */
class Particle {
public:
    /** Identifies this particles position in the particle list */
    int index;

    bool pinned = false;

    glm::vec3 color{ 0.0f, 0.95f, 0.0f };

    double mass = 1;

    vec2 p;
    vec2 v;
    vec2 vTemp;
    vec2 p0;
    vec2 v0;
    vec2 f;

    /**
     * A list of springs that use this particle.  This list is only needed
     * to adjust rest lengths when dragging particles around.
     * This is only used for UI... it is probably not needed for simulation
     */
    std::vector<Spring*> springs {};
    std::vector<BendingSpring*> bendingSprings{};
    /** Default constructor */
    Particle() {}

    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    Particle( double x, double y, double vx, double vy): vTemp(0,0), v(0,0) {
        p0 = vec2(x, y);
        v0 = vec2(vx, vy);
        reset();
    }

    /**
     * Resets the position of this particle
     */
    void reset() {
        p = p0;
        v = v0;
        f = vec2(0, 0);
    }

    /**
     * Clears all forces acting on this particle
     */
    void clearForce() {
        f = vec2(0, 0);
    }

    /**
     * Adds the given force to this particle
     * @param force
     */
    void addForce(vec2 force) {
        f += force;
    }

    /**
     * Computes the distance of a point to this particle
     * @param x
     * @param y
     * @return the distance
     */
    float distance(double x, double y) {
        vec2 diff = p - vec2(x, y);
        return (float) sqrt( diff.x*diff.x + diff.y*diff.y);
    }
};
