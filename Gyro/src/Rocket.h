#pragma once

// definition of classes for model rockets
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "Simbody.h"

using namespace std;
using namespace SimTK;

class Fin
{
public: 
	Fin(double density, double length, double chordlen, double thickness);

	Vec3 calcForces(const State& state, Vector_<SpatialVec>& bodyForces, const MobilizedBody& mobod,
		int printLevel, char ch, bool bApplyForces = 0);

	void getCoeffs(double airspeed, double aoa, double& Cn);

	double m_density;
	double m_length;
	double m_chordlen;
	double m_thickness;
	double m_area;
	double m_mass;
};
