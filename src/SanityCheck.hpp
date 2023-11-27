#pragma once
#include "ParticleSystem.hpp"

/**
 * Implementation of a forgiving edge-edge intersection test.  Guaranteed to return false
 * only if there are definitely two edges intersecting.
 * @author kry
 */
class SanityCheck {

private:
    double eps = 1e-4;

    /**
     * Performs an approximate sanity check.  This code errs on the side of forgiveness
     * so as not to stop the simulation inappropriately.
     * @param system
     * @return true if no intersection, false if things "appear" to be OK
     */
public:
    bool sanityCheck(ParticleSystem* system) {
        for (Particle* p : system->particles) {
            if (isnan(p->p.x)) return false;
            if (isnan(p->p.y)) return false;
            if (isinf(p->p.x)) return false;
            if (isinf(p->p.y)) return false;
        }
        for (Spring* s1 : system->springs) {
            for (Spring* s2 : system->springs) {
                if (s1 == s2) continue;

                // discard the case where the two springs share a particle.
                // if they do happen to be superimposed, then there is likely to be
                // trouble at the next time step.  The one problem not caught will be
                // two segments passing through one another, but at least this leaves the 
                // simulation in a sane non-interpenetrating configuration.
                if (s1->p1 == s2->p1 || s1->p1 == s2->p2 || s1->p2 == s2->p1 || s1->p2 == s2->p2) continue;

                // we'll solve a little linear system for the intersection of two lines, 
                // but give up if the determinant is too close to zero

                double a = s1->p2->p.x - s1->p1->p.x;
                double b = -s2->p2->p.x + s2->p1->p.x;
                double c = s1->p2->p.y - s1->p1->p.y;
                double d = -s2->p2->p.y + s2->p1->p.y;

                double det = a * d - b * c;
                if (abs(det) < 1e-4) continue;// give up, but likely problems will be detected soon!
                double e = s2->p1->p.x - s1->p1->p.x;
                double f = s2->p1->p.y - s1->p1->p.y;
                double alpha1 = 1 / det * (d * e - b * f);
                double alpha2 = 1 / det * (-c * e + a * f);

                if (alpha1 > eps && alpha1 < 1 - eps && alpha2 > eps && alpha2 < 1 - eps) {
                    // bad news!
                    return false;
                }
            }
        }
        return true;
    }
};