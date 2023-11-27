#pragma once
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Particle.hpp"

typedef glm::vec<2, double> vec2;

/**
 * Spring class for 599 assignment 1
 * @author kry 
 */
class Spring {

public:

    Particle* p1;
    Particle* p2;

    /** Spring stiffness, sometimes written k_s in equations */
    float k = 1;
    /** Spring damping (along spring direction), sometimes written k_d in equations */
    float c = 1;
    /** Rest length of this spring */
    double l0 = 0;

    /**
     * Creates a spring between two particles
     * @param p1
     * @param p2
     */
    Spring(Particle* p1, Particle* p2): p1(p1), p2(p2) {
        recomputeRestLength();
        p1->springs.push_back(this);
        p2->springs.push_back(this);
    }

    /**
     * Computes and sets the rest length based on the original position of the two particles
     */
    void recomputeRestLength() {
        vec2 diff = p1->p0 - p2->p0;
        l0 = sqrt(diff.x * diff.x + diff.y * diff.y);
    }

    /**
     * Applies the spring force by adding a force to each particle
     */
    void apply() {
        vec2 springDirection = p1->p - p2->p;
        double len = sqrt(springDirection.x * springDirection.x + springDirection.y * springDirection.y);
        double strain = (len - l0) / len;
        vec2 f = springDirection * (double) (-k * strain);
        p1->f += f;
        p2->f -= f;
    }
};