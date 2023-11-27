
#pragma once
#include "Particle.hpp"
#include "Spring.hpp"


class RobustCCD {
public:

	bool collision;
	bool useRepulsion;
	int maxIterations;
	double restitutionValue;
	double thresholdDistance;
	double eps;
	double repulsionStiffness;

	/** TODO: number of iterations in the last CCD processing loop, to keep an eye on how tricky the problem is */
	int iterations;
	/** TODO: keep track of highest count needed so far */
	int highestIterations = 0;

	// Variables for computing the roots for finding times of collisions
	Particle* particle_A;
	Particle* particle_B;
	Particle* particle_C;

	double a;
	double b;
	double c;

	double Adot1;
	double Adot2;
	double A1;
	double A2;

	double Bdot1;
	double Bdot2;
	double B1;
	double B2;

	double Cdot1;
	double Cdot2;
	double C1;
	double C2;

	double t;
	double t1;
	double t2;

	double alpha;

	int flag = 0;

	RobustCCD::RobustCCD() :
		eps(0.000001),
		collision(true),
		useRepulsion(true),
		maxIterations(100),
		restitutionValue(1.0),
		thresholdDistance(2),
		iterations(0),
		repulsionStiffness(100) {}

	/**
	 * Try to deal with contacts before they happen
	 * @param h
	 * @param system
	 */
	void applyRepulsion(double h, std::vector<Particle*>* particles, std::vector<Spring*>* springs)
	{
		//TODO: apply repulsion on all particles
		// - check for point-edge proximity between all particle-edge pairs
		// - find the normal
		// - take care to deal with segment end points carefully
		// - compute an appropriate  impulse if the distance is less than H
		// - make sure your impulse is in the correct direction!
		// - don't apply the impulse if the relative velocity is separating fast enough (see Bridson et al. 2002)
		// - distribute impulse to the three particles involved in the appropriate manner
		if (useRepulsion)
		{
			std::vector<Particle*> part = *particles;
			std::vector<Spring*> spr = *springs;
			for (Spring* s : spr)
			{
				// Setup A particle
				Particle* A = s->p1;

				// Setup B particle
				Particle* B = s->p2;

				for (Particle* p : part)
				{
					// Skip if Particle p is A or B
					if (p == particle_A || p == particle_B)
					{
						continue;
					}

					// Set particle C
					Particle* C = p;

					// Find normal
					vec2 AB = B->p - A->p;
					vec2 n;
					n[0] = -AB[1];
					n[1] = AB[0];
					n = glm::normalize(n);

					// Vector pointing to C
					vec2 AC = C->p - A->p;

					// Find normal direction
					if (glm::dot(n, AC) < 0)
					{
						//n *= -1;
						n = -n;
					}

					// Find relative velocity in the normal direction
					vec2 avg_rel;
					double avg_x = (A->v[0] + B->v[0]) / 2;
					double avg_y = (A->v[1] + B->v[1]) / 2;
					avg_rel[0] = avg_x;
					avg_rel[1] = avg_y;

					vec2 v_rel = C->v - avg_rel;
					double v_n = glm::dot(v_rel, n);

					// Skip if the particles are not getting closer
					if (v_n <= 0)
					{
						continue;
					}

					// Find the point closest to C on the line segment AB;
					double m1 = (B->p[1] - A->p[1]) / (B->p[0] - A->p[0]);
					double b1 = B->p[1] - (m1 * B->p[0]);

					if (m1 == 0) { continue; }

					double m2 = -1 / m1;
					double b2 = C->p[1] - (m2 * C->p[0]);

					double x = (b2 - b1) / (m1 - m2);
					double y = (m2 * x) + b2;
					vec2 C_closest;
					C_closest[0] = x;
					C_closest[1] = y;

					// Check alpha value
					// p2 = B		p1 = A		p0 = C_closest

					// Calculate alpha
					AB = A->p - B->p;
					double length_AB = glm::length(AB);

					AC = A->p - C_closest;
					double dot_product = glm::dot(AC, AB);
					double alpha_temp = dot_product / (length_AB * length_AB);

					// Check if there is an actual collision
					if (alpha_temp >= 0 && alpha_temp <= 1)
					{
						alpha = alpha_temp;
						checkRepulsion(A, B, C, n, v_n, h, C_closest);
						continue;
					}

					// Double check using epsillon parameter
					else if ((alpha_temp + eps >= 0 && alpha_temp + eps <= 1) || (alpha_temp - eps >= 0 && alpha_temp - eps <= 1))
					{
						alpha = alpha_temp;
						checkRepulsion(A, B, C, n, v_n, h, C_closest);
					}
				}
			}
		}
	}

	bool checkRepulsion(Particle* A, Particle* B, Particle* C, vec2 n, double v_n, double h, vec2 C_close)
	{
		// Calculate "d" value
		double d = thresholdDistance - (glm::dot((C->p - ((alpha)*B->p) - ((1-alpha) * A->p)), n));
		//double d = thresholdDistance - ( glm::dot((C->p - ((alpha)*B->p) - ((1-alpha)*A->p)), n) );
		if (v_n < ((0.1 * d) / h))
		{
			// Apply repulsion
			double j = -min(h*100*d, repulsionStiffness*( ((0.1*d) / h) - v_n) );
			j *= 0.25;

			// Setup masses of A,B and C
			double mA;
			double mB;
			double mC;

			// Check if A is pinned
			if (A->pinned)
			{
				mA = INFINITY;
			}
			else
			{
				mA = A->mass;
			}

			// Check if B is pinned
			if (B->pinned)
			{
				mB = INFINITY;
			}
			else
			{
				mB = B->mass;
			}

			// Check if C is pinned
			if (C->pinned)
			{
				mC = INFINITY;
			}
			else
			{
				mC = C->mass;
			}

			// Update velocities of A, B and C
			A->v += ((-1 * (1-alpha) * n * j) / mA);
			B->v += ((-1 * (alpha)*n * j) / mB);
			C->v += ((n * j) / mC);
			return true;
		}
		return false;
	}

	// Update velocities given the impulses
	void applyImpulses(double restitution, Particle* A, Particle* B, Particle* C)
	{
		cout << "a" << endl;
		vec2 At;
		vec2 Bt;
		vec2 Ct;
		//At = A->p + A->v * t;
		//Bt = B->p + B->v * t;
		//Ct = C->p + C->v * t;
	
		// Find the normal in either direction
		vec2 n;
		double dx = B->p[0] - A->p[0];
		//double dx = Bt[0] - At[0];
		double dy = B->p[1] - A->p[1];
		//double dy = Bt[1] - At[1];
		n[0] = (-1) * dy;
		n[1] = dx;
		n = glm::normalize(n);
	
		// Find parameter "V-"
		vec2 V_minus;
		//V_minus = (alpha*A->v) + ( (1-alpha)*B->v ) - (C->v);
		V_minus = ((1-alpha) * A->v) + (alpha * B->v) - (C->v);

		// Setup masses of A,B and C
		double mA;
		double mB;
		double mC;

		// Check if A is pinned
		if (A->pinned)
		{
			mA = INFINITY;
		}
		else
		{
			mA = A->mass;
		}

		// Check if B is pinned
		if (B->pinned)
		{
			mB = INFINITY;
		}
		else
		{
			mB = B->mass;
		}

		// Check if C is pinned
		if (C->pinned)
		{
			mC = INFINITY;
		}
		else
		{
			mC = C->mass;
		}

		// Find parameter "j" (impulse)
		double j;
		//j = ( (1 + restitution) * glm::dot(n, V_minus) ) / ( (alpha*alpha)/mA + ((1-alpha)*(1-alpha)) / mB + 1 / mC);
		j = ((1 + restitution) * glm::dot(n, V_minus)) / ((alpha * alpha) / mB + ((1 - alpha) * (1 - alpha)) / mA + 1 / mC);

		// Update velocities of A, B and C
		A->v += ((-1 * (1-alpha) * n * j) / mA);
		B->v += ((-1 * (alpha) * n * j) / mB);
		C->v += ((n * j) / mC);
	}

	/**
	* Checks all collisions in interval t to t+h
	* @param h
	* @param system
	* @return true if all collisions resolved
	*/
	bool checkCollision(double h, std::vector<Particle*>* particles, std::vector<Spring*>* springs) {
		// TODO: Most of your assignment code will go here
		// - put your iterative solving code here!
		// - you can use the nextCollision function below to deal with an individual particle-edge pair
		// For each particle-edge pair, find the roots for when the three particles are
		// co-linear, and then pick the first root on (0,h] which corresponds to an actually collision.
		// compute a collision response.  That is, compute an appropriate collision normal, compute the 
		// impulse, and then apply the impulse to the associated particles.  Be sure to deal with pinning 
		// constraints!
		std::vector<Particle*> part = *particles;
		std::vector<Spring*> spr = *springs;
		flag = 0;
		for (int i = 0; i < maxIterations; i++)
		{
			for (Spring* s : spr)
			{
				// Setup A particle
				particle_A = s->p1;

				// Setup B particle
				particle_B = s->p2;

				// Look for co-linearity with the other particles
				for (Particle* p : part)
				{
					// Skip if Particle p is A or B
					if (p == particle_A || p == particle_B)
					{
						continue;
					}

					//Setup C particle
					particle_C = p;

					// Check if there is a collision
					if (nextCollision(h, restitutionValue, particle_A, particle_B, particle_C))
					{
						flag == 1;
					}
				}
			}
		
			if (flag == 0)
			{
				break;
			}
			flag = 0;
		}
		
		return true;
	}

	/**
	* Processes next collision on (0,h] if it exists, between p0 and segment connecting p1 and p2.
	*
	* Watch out for
	* - normal direction
	* - epsilon tests
	* - segment connecting p1 and p2 at zero length! (unlikely?)
	*
	* @param h				timestep
	* @param restitution	bouncyness parameter
	* @param p0				particle
	* @param p1				line segment end point
	* @param p2				line segment end point
	* @return true if collision processed
	*/
	bool RobustCCD::nextCollision(double h, double restitution, Particle* A, Particle* B, Particle* C) 
	{
		//TODO: 
		// * - finds roots
		// * - checks that p0 actually falls on segment
		// * - processes collision by computing and applying impulses
		// * - returns true if collision processed

		// Setup position and velocity variables for each particle for clarity
		A1 = A->p[0];
		A2 = A->p[1];
		Adot1 = A->v[0];
		Adot2 = A->v[1];

		B1 = B->p[0];
		B2 = B->p[1];
		Bdot1 = B->v[0];
		Bdot2 = B->v[1];

		C1 = C->p[0];
		C2 = C->p[1];
		Cdot1 = C->v[0];
		Cdot2 = C->v[1];

		// Calculate the coeficcients
		a = Adot1 * Bdot2 - Adot2 * Bdot1 - Adot1 * Cdot2 + Adot2 * Cdot1 + Bdot1 * Cdot2 - Bdot2 * Cdot1;
		b = A1 * Bdot2 - A2 * Bdot1 + Adot1 * B2 - Adot2 * B1 - A1 * Cdot2 + A2 * Cdot1 - Adot1 * C2 + Adot2 * C1 + B1 * Cdot2 - B2 * Cdot1 + Bdot1 * C2 - Bdot2 * C1;
		c = A1 * B2 - A2 * B1 - A1 * C2 + A2 * C1 + B1 * C2 - B2 * C1;

		// Check special case
		if (a == 0 && b == 0 && c == 0)
		{
			// Co-linear
			if (checkAlpha(A, B, C))
			{
				// There is a collision
				//A->pinned = true;
				//B->pinned = true;
				//C->pinned = true;
				applyImpulses(restitution, A, B, C);
				return true;
			}

			return false;
		}

		// Check if we have to deal with a linear equation
		if (a == 0)
		{
			if (b == 0)
			{
				// Not co-linear
				return false;
			}

			t = (-1 * c) / b;
			if (t > 0 && t <= h + eps)
			{
				// Co-linear
				if (checkAlpha(A, B, C))
				{
					// There is a collision
					//A->pinned = true;
					//B->pinned = true;
					//C->pinned = true;
					applyImpulses(restitution, A, B, C);
					return true;
				}
			}

			return false;
		}

		/* If not, need to check the roots */
		
		// Take the discriminant
		double discriminant = (b * b) - (4 * a * c);

		// If discriminant is 0 there is no collision, just return true
		if (discriminant < 0)
		{
			// Not co-linear
			return false;
		}

		// If not, need to ckeck all cases of the roots
		t1 = ((-1 * b) - sqrt(discriminant)) / (2 * a);
		t2 = ((-1 * b) + sqrt(discriminant)) / (2 * a);

		// Both time are not in the interval
		if ( (t1 <= 0 && t2 <= 0) || (t1 > h + eps && t2 > h + eps))
		{
			// Not co-linear
			return false;
		}

		// t1 is valid ONLY
		else if ( (t1 > 0 && t1 <= h + eps) && (t2 <= 0 || t2 > h + eps))
		{
			// Co-linear
			t = t1;
			if (checkAlpha(A, B, C))
			{
				// There is a collision
				//A->pinned = true;
				//B->pinned = true;
				//C->pinned = true;
				applyImpulses(restitution, A, B, C);
				return true;
			}

			return false;
		}

		// t2 is valid ONLY
		else if ( (t2 > 0 && t2 <= h + eps) && (t1 <= 0 || t1 > h + eps))
		{
			// Co-linear
			t = t2;
			if (checkAlpha(A, B, C))
			{
				// There is a collision
				//A->pinned = true;
				//B->pinned = true;
				//C->pinned = true;
				applyImpulses(restitution, A, B, C);
				return true;
			}

			return false;
		}

		// Both are valid, need to take the smallest time
		else if ( (t1 > 0 && t1 <= h + eps) && (t2 > 0 && t2 <= h + eps))
		{
			if (t1 <= t2)
			{
				// Co-linear
				t = t1;
				if (checkAlpha(A, B, C))
				{
					// There is a collision
					//A->pinned = true;
					//B->pinned = true;
					//C->pinned = true;
					applyImpulses(restitution, A, B, C);
					return true;
				}

				return false;
			}
			
			else
			{
				// Co-linear
				t = t2;
				if (checkAlpha(A, B, C))
				{
					// There is a collision
					//A->pinned = true;
					//B->pinned = true;
					//C->pinned = true;
					applyImpulses(restitution, A, B, C);;
					return true;
				}
			}
		}
		//cout << "a: " << a << " b: " << b << " c: " << c << endl;
		return false;
	}

	// Helper function added to calculate alpha parameter
	bool RobustCCD::checkAlpha(Particle* A, Particle* B, Particle* C)
	{
		// p2 = B		p1 = A		p0 = C
		vec2 AB; 
		vec2 AC;
		
		vec2 At;
		vec2 Bt;
		vec2 Ct;
		At = A->p + A->v * t;
		Bt = B->p + B->v * t;
		Ct = C->p + C->v * t;

		// Calculate alpha
		//AB = A->p - B->p;
		AB = At - Bt;
		double length_AB = glm::length(AB);

		//AC = A->p - C->p;
		AC = At - Ct;

		double dot_product = glm::dot(AC, AB);
		double alpha_temp = dot_product / (length_AB * length_AB);
		// Check if there is an actual collision
		if (alpha_temp >= 0 && alpha_temp <= 1)
		{
			alpha = alpha_temp;
			return true;
		}
		
		// Double check using epsillon parameter
		else if ( (alpha_temp + eps >= 0 && alpha_temp + eps <= 1) || (alpha_temp - eps >= 0 && alpha_temp - eps <= 1) )
		{
			alpha = alpha_temp;
			return true;
		}
		
		// Return false if there was not a collision
		return false;
	}

};

