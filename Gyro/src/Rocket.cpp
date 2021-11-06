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
	m_area = length * chordlen;
	m_mass = density * length * chordlen * thickness;  
};

// Get aerodynamic forces on a rectangular flate-plate fin
// Assumes local right-handed reference frame centered at the center of mass,
// where:
//   Z-axis is from trailing edge to leading edge (chord axis, -axial force axis)
//   Y-axis is normal to the plane of the fin (normal force axis)
//   X-axis is from root to tip (spanwise axis)
// Currently only calculates the normal force along Y, 
//    using the Y and Z components of the relative wind
// Assumes that this force is centered at the center of mass
// Inputs:
//   state = Simbody state
//   mobod = Simbody mobilized body representinmg a Fin
//   printlevel = 0 for none, 1 for some, 2 for all
//   ch = character to print after "Fin", usually the fin number
//   bApplyForces = true to apply forces to the Fin mobilized body
// Returns:
//   calculated force vector in the Grouind frame as a Vec3
//
Vec3 Fin::calcForces(const State& state, Vector_<SpatialVec>& bodyForces, 
	const MobilizedBody& mobod, int printLevel, char ch, bool bApplyForces)
{
	// In still air the relative wind velocity is the inverse of the Fin velocity vector
	// however, the angle of attack is the angle between the Fin +Z axis
	// and the direction that the relative wind is coming FROM so we use the inverse 
	// of the rel wind vector which is just the Fin velocity vector
	// for information, get unit vector along Fin Z-axis in Ground frame
	Vec3 UZG = mobod.expressVectorInGroundFrame(state, Vec3(0, 0, 1));
	Vec3 velFinG = mobod.getBodyOriginVelocity(state);	// Fin vel in Ground frame
	// inverse of relative wind velocity in Fin frame
	Vec3 invVelRelW = mobod.expressGroundVectorInBodyFrame(state, velFinG);
	// get magnitude and AOA of relative wind (chord is along Z axis)
	double Vy = invVelRelW[1];
	double Vz = invVelRelW[2];
	double Vmag = sqrt(Vz * Vz + Vy * Vy);
	double AOA = atan2(-Vy, Vz) * 180/Pi;
	double Cn;	// coefficient of normal force
	getCoeffs(Vmag, AOA, Cn);
	// calculate force
	const double rho = 1.204;     // density of dry air in kg / m ^ 3
	double Fy = 0.5*rho*Vmag*Vmag*m_area*Cn;
	Vec3 FyG = mobod.expressVectorInGroundFrame(state, Vec3(0, Fy, 0));
	if (bApplyForces)
	{
		mobod.applyForceToBodyPoint(state, Vec3(0, 0, 0), FyG, bodyForces);
	}
	if (printLevel)
	{
		printf("Fin%c: UZxyz  %6.3f, %6.3f, %6.3f, Vxyz %6.3f, %6.3f, %6.3f, RW mag %6.3f AOA %6.3f Fy %6.3f\r\n",
			ch, UZG[0], UZG[1], UZG[2], velFinG[0], velFinG[1], velFinG[2], Vmag, AOA, Fy);
	}
	return FyG;
};

// Calculate coefficient of normal force for a Fin from the relative wind
// Only considers Y and Z components of wind, ignores X (spanwise) component  
// Uses simplified version of technique from AeroVECTOR
// Inputs:
//   airspeed = magnitude of relative wind velocity
//   aoa = angle of attack in degrees
// Returns with:
//   Cn set to normal force coefficient
void Fin::getCoeffs(const double airspeed, const double aoa, double& Cn)
{
	double aoa_sign;
	if (aoa < 0)
		aoa_sign = -1.0;
	else
		aoa_sign = 1.0;
	double a = abs(aoa);
	if (a > 180.0)
		printf("Fin AOA = %f, out of range\r\n"); 
	if (a > 90.0)
		a = 180.0 - a;
	if (a > 50)
		Cn = aoa_sign * 1.15;
	else
		Cn = aoa_sign * a * (1.9 / 50.0);
}