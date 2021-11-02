// Rocket.cpp
// Implementation of classes for model rockets

#include <vector>
#include "Simbody.h"
#include "Rotor.h"
#include "Rocket.h"

using namespace std;
using namespace SimTK;

// Constructor for a rocket fin
// Creates a simple rectangular fin with a flat-plate airfoil
// Inputs:
//   density = material density (kg/m^3)
//   length = spanwise dimension (X-axis)
//   chordlen = fore-aft dimension (Z-axis)
//   thickness = flat-plate thickness
//
Fin::Fin(double density, double length, double chordlen, double thickness)
{
	m_density = density;
	m_length = length;
	m_chordlen = chordlen;
	m_thickness = thickness;
	m_mass = density * length * chordlen * thickness;  
};

// Get aerodynamic forces on a rectangular flate-plate fin
// Assumes local right-handed reference frame centered at the center of mass,
// where:
//   Z-axis is from trailing edge to leading edge
//   Y-axis is normal to the plane of the fin
//   X-axis is from root to tip
//  Currently only calculates the normal force along Y, 
//    using the Y and Z components of the relative wind
//  Assumes that this force is centered at the center of mass
//  Inputs:
//    state = Simbody state
//    mobod = Simbody mobilized body representinmg a Fin
//    printlevel = 0 for none, 1 for some, 2 for all
//
void Fin::getForces(const State& state, const MobilizedBody& mobod, int printLevel)
{
	// in still air the relative wind velocity is the reciprocal of the fin velocity
	Vec3 velFin = mobod.getBodyOriginVelocity(state);	// in Ground frame
	// now get air velocity in fin frame
	Vec3 velRelW = mobod.expressGroundVectorInBodyFrame(state, -velFin);
	// reverse Vz for calculating AOA
	double Vy = velRelW[1];
	double Vz = -velRelW[2];
	double Vmag = sqrt(Vz * Vz + Vy * Vy);
	double AOA = atan2(Vy, Vz);
	if (printLevel)
	{
		printf("fin: Vmag %6.3f AOA %6.3f\r\n", Vmag, AOA);
	}
};
