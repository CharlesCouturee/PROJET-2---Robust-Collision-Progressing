#pragma once
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <Eigen/Dense>
using Eigen::MatrixXf;
using Eigen::VectorXf;

#include "Particle.hpp"

/**
 * Leaf Spring class for COMP 599 assignment 2
 * (do not distribute)
 * @author kry & Mercier-Aubin
 */
class BendingSpring {

public:
    Particle* p0;
    Particle* p1;
    Particle* p2;

    /** Modulation for bending stiffness */
    double kbMod = 1e4;

    /** Spring bending base stiffness */
    double kbBase = 1;

    /** Spring bending stiffness */
    double kb = kbBase * kbMod;

    /** Rest angle */
    double theta0 = 0;

    /**
     * Creates a spring connecting two particles.
     * The rest length should be set
     * @param p0
     * @param p1
     * @param p2
     */
    BendingSpring(Particle* p0, Particle* p1, Particle* p2):p0(p0),p1(p1),p2(p2) {
        computeRestAngle();
    }

    /**
     * Computes the rest length of the connected particles
     */
    void computeRestAngle() {
        double Ax = p0->p.x;
        double Ay = p0->p.y;
        double Bx = p1->p.x;
        double By = p1->p.y;
        double Cx = p2->p.x;
        double Cy = p2->p.y;
        double t1 = Ay - By;
        double t2 = Cx - Bx;
        double t4 = Bx - Ax;
        double t5 = Cy - By;
        theta0 = atan2(t1 * t2 + t4 * t5, t4 * t2 - t1 * t5);
    }

    /** Apply forces to particles */
    void apply() {
        kb = kbBase * kbMod;
        double Ax = p0->p.x;
        double Ay = p0->p.y;
        double Bx = p1->p.x;
        double By = p1->p.y;
        double Cx = p2->p.x;
        double Cy = p2->p.y;
        double t1 = Ay - By;
        double t2 = Cx - Bx;
        double t3 = Bx - Ax;
        double t4 = Cy - By;
        double t5 = t1 * t2 + t3 * t4;
        double t6 = -t1;
        double t7 = t3 * t2 + t6 * t4;
        t4 = -t4;
        double t8 = 1 / t7;
        double t9 = t8 * t5;
        double t10 = 1 + t5 * t5 * t8 * t8;
        t10 = 0.1e1 / t10;
        t5 = kb * (atan2(t5, t7) - theta0);
        if (abs(t5) < 1e-5) return;
        // note t8 is probably infinity in this case?
        p0->f.x += -t5 * t8 * (t4 + t9 * t2) * t10;
        p0->f.y += -t5 * t8 * (t2 - t9 * t4) * t10;
        p1->f.x += -t5 * t8 * (Cy - Ay - t9 * (Cx - 2 * Bx + Ax)) * t10;
        p1->f.y += -t5 * t8 * (-Cx + Ax - t9 * (Cy - 2 * By + Ay)) * t10;
        p2->f.x += -t5 * t8 * (t1 - t9 * t3) * t10;
        p2->f.y += -t5 * t8 * (t3 - t9 * t6) * t10;
    }

    /** Stiffness matrix */
    double K[6][6];

    void computeAngleBasedStiffnessKdx() {
        double Ax = p0->p.x;
        double Ay = p0->p.y;
        double Bx = p1->p.x;
        double By = p1->p.y;
        double Cx = p2->p.x;
        double Cy = p2->p.y;
        double t1 = -Cy + By;
        double t2 = Bx - Ax;
        double t3 = Cx - Bx;
        double t4 = By - Ay;
        double t5 = -t1;
        double t6 = t2 * t3 + t4 * t5;
        t5 = -t3 * t4 + t2 * t5;
        double t7 = -t3;
        double t8 = 0.1e1 / t6;
        double t9 = t5 * t8;
        double t10 = t9 * t7;
        double t11 = t1 - t10;
        double t12 = t8 * t11;
        double t13 = t8*t8;
        double t14 = 1 + t5*t5 * t13;
        t6 = atan2(t5, t6) - theta0;
        double t15 = 2 * t13;
        double t16 = t15 * t5;
        t14 = 1.0 / t14;
        double t17 = t14 * kb;
        double t18 = t9 * t1;
        double t19 = t3 - t18;
        double t20 = t8 * t19;
        t18 = -t3 + 2 * t18;
        double t21 = t6 * t16 * t19;
        double t22 = t17 * ((-t20 + t21) * t14 * t12 - t6 * t13 * (t18 * t7 - t1*t1));
        double t23 = Cy - Ay;
        double t24 = Cx - 2 * Bx + Ax;
        double t25 = t9 * t24;
        double t26 = t23 - t25;
        double t27 = t8 * t26;
        t10 = -t1 + 2 * t10;
        double t28 = t6 * t16 * t26;
        double t29 = (-t27 + t28) * t14;
        double t30 = t17 * (t29 * t12 - t6 * t13 * (t10 * t24 - t23 * t7 - t5));
        double t31 = -Cx + Ax;
        double t32 = Cy - 2 * By + Ay;
        double t33 = t31 - t9 * t32;
        double t34 = t8 * t33;
        double t35 = t6 * t16 * t33;
        double t36 = (-t34 + t35) * t14;
        t10 = t17 * (t36 * t12 - t6 * t8 * (t8 * (t10 * t32 - t31 * t7) + 1));
        double t37 = t9 * t2;
        double t38 = -t4 - t37;
        double t39 = t8 * t38;
        t37 = t4 + 2 * t37;
        double t40 = t1 * t2;
        double t41 = t6 * t16 * t38;
        double t42 = (-t39 + t41) * t14;
        double t43 = t17 * (t42 * t12 - t6 * t13 * (t37 * t7 - t40 + t5));
        t9 = t9 * t4;
        double t44 = t2 - t9;
        double t45 = t8 * t44;
        t9 = -t2 + 2 * t9;
        double t46 = t1 * t4;
        double t47 = t6 * t16 * t44;
        double t48 = (-t45 + t47) * t14;
        double t49 = t17 * (t48 * t12 - t6 * t8 * (t8 * (t9 * t7 - t46) - 1));
        t29 = t17 * (t29 * t20 - t6 * t8 * (t8 * (t18 * t24 - t23 * t1) - 1));
        double t50 = t17 * (t36 * t20 - t6 * t13 * (t18 * t32 - t31 * t1 - t5));
        t18 = t17 * (t42 * t20 - t6 * t8 * (t8 * (t18 * t2 + t1 * t4) + 1));
        t3 = t17 * (t48 * t20 - t6 * t13 * ((2 * t46 * t8 + 1) * t5 - t3 * t4 - t40));
        t25 = t17 * (t36 * t27 - t6 * t13 * ((-t23 + 2 * t25) * t32 - t31 * t24));
        t36 = t17 * (t42 * t27 - t6 * t13 * (t37 * t24 - t2 * t23 - t5));
        t23 = t17 * (t48 * t27 - t6 * t8 * (t8 * (t9 * t24 - t23 * t4) + 1));
        t8 = t17 * (t42 * t34 - t6 * t8 * (t8 * (t37 * t32 - t2 * t31) - 1));
        t9 = t17 * (t48 * t34 - t6 * t13 * (t9 * t32 - t31 * t4 - t5));
        t13 = t17 * (t48 * t39 - t6 * t13 * (t37 * t4 - t2 * t2));
        K[0][0] = t17 * (-t12*t12 * t14 + t6 * t15 * t7 * t11 + t6 * t12 * t14 * t16 * t11);
        K[0][1] = t22;
        K[0][2] = t30;
        K[0][3] = t10;
        K[0][4] = t43;
        K[0][5] = t49;
        K[1][0] = t22;
        K[1][1] = t17 * (-t20*t20 * t14 + t6 * t15 * t1 * t19 + t21 * t20 * t14);
        K[1][2] = t29;
        K[1][3] = t50;
        K[1][4] = t18;
        K[1][5] = t3;
        K[2][0] = t30;
        K[2][1] = t29;
        K[2][2] = t17 * (-t27*t27 * t14 - t6 * t15 * (t5 - t26 * t24) + t28 * t27 * t14);
        K[2][3] = t25;
        K[2][4] = t36;
        K[2][5] = t23;
        K[3][0] = t10;
        K[3][1] = t50;
        K[3][2] = t25;
        K[3][3] = t17 * (-t34*t34 * t14 - t6 * t15 * (t5 - t33 * t32) + t35 * t34 * t14);
        K[3][4] = t8;
        K[3][5] = t9;
        K[4][0] = t43;
        K[4][1] = t18;
        K[4][2] = t36;
        K[4][3] = t8;
        K[4][4] = t17 * (-t39*t39 * t14 + t6 * t15 * t2 * t38 + t41 * t39 * t14);
        K[4][5] = t13;
        K[5][0] = t49;
        K[5][1] = t3;
        K[5][2] = t23;
        K[5][3] = t9;
        K[5][4] = t13;
        K[5][5] = t17 * (-t45*t45 * t14 + t6 * t15 * t4 * t44 + t47 * t45 * t14);
    }

};