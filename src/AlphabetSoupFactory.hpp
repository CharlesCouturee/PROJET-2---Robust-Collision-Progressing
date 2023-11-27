#pragma once
#include "vector"
#include "Particle.hpp"
#include "ParticleSystem.hpp"
#include "Spring.hpp"
#include <freetype/ftglyph.h>
#include <freetype/freetype.h>
#include <freetype/ftimage.h>
#include <freetype/ftoutln.h>

using namespace std;

int MoveToFunc(const FT_Vector* to, void* user);
int LineToFunc(const FT_Vector* to, void* user);
int ConicToFunc(const FT_Vector* control, const FT_Vector* to, void* user);
int CubicToFunc(const FT_Vector* control1, const FT_Vector* control2, const FT_Vector* to, void* user);

/**
 * Letter generation helper class.  Deals with ugly issues such as trying to get an
 * approximately constant segment length piece-wise representation of letters on both
 * straight and curved segements. 
 *
 * May not work so well with different fonts
 *
 * @author kry & Mercier-Aubin
 */
class AlphabetSoupFactory {

public:

	double offsetx = 50;
	double offsety = 250;
	
	int segmentMax = 30 / 2;
	string text = "test";

	ParticleSystem* system;

	double fontsize = 64;
	string fontName = "consola.ttf";
	FT_Face face;
	FT_Error error;
	FT_Outline_Funcs funcs;

	int segmentMin = 15;

	AlphabetSoupFactory( ParticleSystem* system, string RES_DIR ) {
		this->system = system;
		int style = 0;
		FT_Library library;
		error = FT_Init_FreeType(&library);
		string fontLocation = RES_DIR + fontName;
		error = FT_New_Face(library, fontLocation.c_str(), 0, &face);
		if (error) {
			cout << "could not find " << fontName << endl;
			return;
		}
		error = FT_Set_Char_Size(face, fontsize, fontsize, 300, 300);

		funcs.move_to = MoveToFunc;
		funcs.line_to = LineToFunc;
		funcs.conic_to = ConicToFunc;
		funcs.cubic_to = CubicToFunc;
		funcs.shift = 0;
		funcs.delta = 0;
	}

	vec2 evaluateQuadraticBezier(double t, const vec2& P0, const vec2& P1, const vec2& P2) {
		double tinv = 1 - t;
		return tinv * tinv * P0 + 2 * tinv * t * P1 + t * t * P2;
	}

	vec2 evaluateCubicBezier(double t, const vec2& P0, const vec2& P1, const vec2& P2, const vec2& P3) {
		return (1 - t) * evaluateQuadraticBezier(t, P0, P1, P2) + t * evaluateQuadraticBezier(t, P1, P2, P3);
	}


	/**
	 * Creates the specified letter at the specified offset location
	 */
	void createLetter(ParticleSystem* system, char c, double offx, double offy ) {
		if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
			std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
			return;
		}
		FT_Int32  load_flags = FT_LOAD_DEFAULT;
		FT_Glyph glyph;
		error = FT_Get_Glyph(face->glyph, &glyph);
		FT_Outline Outg = face->glyph->outline;
		ox = offx;
		oy = offy;		
		FT_Outline_Decompose(&Outg, &funcs, (void*)this);
	}

	double ox; 	// offset of character to make
	double oy;
	vec2 start; // starting point of contour
	vector<vec2> curve; // part of a contour
	bool prevWasCurve; // flag for keeping track of curved regions
	vec2 prevTangent; // previous tangent to match with next in curved regions
	vector<vec2> particlesToMake; // coarse discretization to make

	/** Make a coarse discretization of the current piecewise curve built so far */
	void discretizeCurve() {
		if (curve.size() <= 1) return;
		double len = 0;
		for (int i = 1; i < curve.size(); i++) {
			vec2 d = curve[i] - curve[i - 1];
			len += sqrt(d.x * d.x + d.y * d.y);
		}
		if (len < 1e-4) { // skip this? probably shouldn't happen!
			prevWasCurve = false;
			//cout << "Empty curve, ignore";
			vec2 p = curve.back();
			curve.clear();
			curve.push_back(p);
			return;
		}
		// choose how many pieces
		int n = (int)ceil(len / segmentMin);
		double targetSegLen = len / n;
		double seglen = 0;
		int count = 0;
		// if partilesToMake already has points, then we don't want the first
		if (particlesToMake.size() == 0) {
			particlesToMake.push_back(curve[0]);
		}
		for (int i = 1; i < curve.size(); i++) {
			vec2 d = curve[i] - curve[i - 1];
			double l = sqrt(d.x * d.x + d.y * d.y);
			if (l < 1e-6) continue;
			seglen += l;
			while (seglen > targetSegLen) {
				seglen -= targetSegLen;
				double alpha = (l-seglen) / l;
				particlesToMake.push_back( (1-alpha) * curve[i - 1] + alpha * curve[i] );
				count++;
				if (count >= n - 1) break;
			}
			if (count >= n - 1) break;
		}
		particlesToMake.push_back(curve.back());
		prevWasCurve = false;
		curve.clear();
		curve.push_back( particlesToMake.back() ); // last point for starting next part of path
	}

	/** Create a the particles for the coarse discretization of the closed contour */
	void createClosedContour() {
		if (particlesToMake.size() < 3) {
			particlesToMake.clear();
			curve.clear();
			//cout << "skipping bad closed contour of less than 3 points" << endl;
			return; // ignore empty contours
		}
		Particle* p0 = NULL;
		Particle* p1 = NULL;
		Particle* p2 = NULL;
		Particle* first = NULL;
		Particle* second = NULL;
		// particles to make should include repeated particle to close loop, so igore it
		for (int i = 0; i < particlesToMake.size()-1; i++ ) {
			vec2 pos = particlesToMake[i];
			p2 = system->createParticle(pos.x + ox, -pos.y + oy, 0, 0);
			if (i == 0) first = p2;
			if (i == 1) second = p2;
			if (i > 0) system->springs.push_back(new Spring(p1, p2));
			if (i > 1)  system->bendingSprings.push_back(new BendingSpring(p0, p1, p2));
			p0 = p1;
			p1 = p2;
			p2 = NULL;
		}
		// close the loop
		system->springs.push_back(new Spring(p1, first));
		system->bendingSprings.push_back(new BendingSpring(p0, p1, first));
		system->bendingSprings.push_back(new BendingSpring(p1, first, second));
		particlesToMake.clear();
		curve.clear();
	}

	int MoveTo(const FT_Vector* to ) {
		//cout << "M " << to->x << " " << to->y << endl;
		start = vec2(to->x, to->y);
		curve.push_back(start);
		prevWasCurve = false;
		return 0;
	}
	int LineTo(const FT_Vector* to ) {
		//cout << "L " << to->x << " " << to->y << endl;
		if (prevWasCurve) {
			discretizeCurve();
		}
		vec2 p = vec2(to->x, to->y);
		curve.push_back(p);
		discretizeCurve();
		if (start.x == p.x && start.y == p.y ) { // close contour (expect move to next, or no more calls)
			createClosedContour();
		}
		prevWasCurve = false;
		return 0;
	}

	int ConicTo(const FT_Vector* control, const FT_Vector* to ) {
		//cout << "Q " << control->x << " " << control->y << " " << to->x << " " << to->y << endl;
		vec2 c0 = curve.back();
		vec2 c1 = vec2(control->x, control->y);
		vec2 c2 = vec2(to->x, to->y);
		vec2 startTangent = c1 - c0;
		double continuity = abs(prevTangent.x * startTangent.y - prevTangent.y * startTangent.x);
		continuity /= sqrt(prevTangent.x * prevTangent.x + prevTangent.y * prevTangent.y);
		continuity /= sqrt(startTangent.x * startTangent.x + startTangent.y * startTangent.y);
		if (!prevWasCurve || continuity > 0.1) {
			// no previous curve, or no continuity, so do coarse discretization
			discretizeCurve();
		}
		for (int i = 0; i <= 10; i++) { // arbitrary resolution for fine discretization
			double t = i / 10.0;
			curve.push_back(evaluateQuadraticBezier(t, c0, c1, c2));
		}
		if (start.x == c2.x && start.y == c2.y) { // close contour (expect move to next, or no more calls)
			discretizeCurve();
			createClosedContour();
		} else {
			prevTangent = c2 - c1;
			prevWasCurve = true;
		}
		return 0;
	}
	int CubicTo(const FT_Vector* control1, const FT_Vector* control2, const FT_Vector* to ) {
		//cout << "C " << control1->x << " " << control1->y << " " << control2->x << " " << control2->y << " " << to->x << " " << to->y << endl;
		vec2 c0 = curve.back();
		vec2 c1 = vec2(control1->x, control1->y);
		vec2 c2 = vec2(control2->x, control2->y);
		vec2 c3 = vec2(to->x, to->y);
		vec2 startTangent = c1 - c0;
		if (!prevWasCurve || abs(prevTangent.x * startTangent.y - prevTangent.y * startTangent.x) > 1e-4) {
			// no previous curve, or no continuity, so do coarse discretization
			discretizeCurve();
		}
		for (int i = 0; i <= 10; i++) { // arbitrary resolution for fine discretization
			double t = i / 10.0;
			curve.push_back(evaluateQuadraticBezier(t, c0, c1, c2));
		}
		if (start.x == c3.x && start.y == c3.y) { // close contour (expect move to next, or no more calls)
			discretizeCurve();
			createClosedContour();
		} else {
			prevTangent = c3 - c2;
			prevWasCurve = true;
		}
		return 0;
	}
};

int MoveToFunc(const FT_Vector* to, void* user) {
	return ((AlphabetSoupFactory*)user)->MoveTo(to);
}
int LineToFunc(const FT_Vector* to, void* user) {
	return ((AlphabetSoupFactory*)user)->LineTo(to);
}
int ConicToFunc(const FT_Vector* control, const FT_Vector* to, void* user) {
	return ((AlphabetSoupFactory*)user)->ConicTo(control, to);
}
int CubicToFunc(const FT_Vector* control1, const FT_Vector* control2, const FT_Vector* to, void* user) {
	return ((AlphabetSoupFactory*)user)->CubicTo(control1, control2, to);
}
