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

	void getForces(const State& state, const MobilizedBody& mobod, int printLevel);

	double m_density;
	double m_length;
	double m_chordlen;
	double m_thickness;
	double m_mass;
};
