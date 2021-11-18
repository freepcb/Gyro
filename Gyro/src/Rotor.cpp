// Implementation of classes for handling rotor aerodynamics

#include <vector>
#include "Simbody.h"
#include "Rotor.h"
using namespace std;
using namespace SimTK;

#define PI 3.14159265358979323846

// note that pitch is in degrees
RotorBlade::RotorBlade(Airfoil* af, double rootR, double tipR, double chordLen,
	double pitch, int nsegs)
{
	m_af = af;
	m_chordLen = chordLen;
	m_pitch = pitch;
	m_nsegs = nsegs;
	m_bladeLenX = tipR - rootR;	  // length of blade long X axis
	double segLen = m_bladeLenX / nsegs;  // length of each segment
	m_segArea = m_chordLen * segLen;	  // area of segment			
	// get position of segment centers from hub axis and blade center
	for (int iseg = 0; iseg < nsegs; iseg++)
	{
		m_seg_x.push_back(iseg * segLen + 0.5*segLen - (m_bladeLenX / 2));	 // x-coord of segment center
	}
	m_rey_coef = 70698 * m_chordLen;	// Reynolds number for v = 1 m/sec
};

// Implementation of RotorBlade::getForces() 
// Calculates aerodynamic forces on a rotor blade
// Assumes that the blade is a rectngular structure attached to a hub where 
//   the axis of rotation of the hub determines the forward direction of the blade
// The end of the blade closest to the hub is called the "root" and the opposite end
//   is called the "tip"
// Uses a local frame of reference where:
//   the origin is at the center of mass
//   the X axis is along the long axis of the blade
//   the Y axis is along the forward direction of rotation
//	 the Z axis may be considered "up" relative to the X-Y plane of rotation
// Assumes that the cross-section of the blade is an airfoil, where the angle 
//   between the chord of the iirfoil and the Y axis is called the "pitch"
// Divides the blade into sections to approximate continuous variation along the blade
// Calculates net forces and torques on/around the Z and Y local axes from the sections
//   and coverts these to equivalent forces at the root and tip ends of the blade
// Inputs:
//   state = SimBody state of the rigid-body simulation 
//   mobod = mobilized body for the blade as defined in SimBody 
//   windVelG = external wind velocity in the Ground frame (if present)
//   printLevel = 0 for none, 1 to print outputs, 2 to print details
// Outputs (these are arguments defined as references which are modified):
//   FZroot = equivalent force along local Z axis at the root
//   FZtip = equivalent force along local Z axis at the tip
//   FYroot = equivalent force along local Y axis at the root
//   FYtip = equivalent force along local Y axis at the tip
//
void RotorBlade::getForces(const State& state, const MobilizedBody& mobod, Vec3 windVelG, 
		int printLevel, double& FZroot, double& FZtip, double& FYroot, double& FYtip)
{
	static int mbi_of_first_blade = -1;

	const double rho = 1.204;     // density of dry air in kg / m ^ 3
	double netLift = 0.0;
	double netLiftMoment = 0.0;
	double netThrust = 0.0;
	double netThrustMoment = 0.0;
	// create vectors for printing
	vector<double> x_seg(m_nsegs);
	vector<double> spd(m_nsegs);
	vector<double> aoa(m_nsegs);
	vector<double> CL(m_nsegs);
	vector<double> CD(m_nsegs);
	vector<double> L(m_nsegs);
	vector<double> Lup(m_nsegs);
	vector<double> Lfwd(m_nsegs);
	vector<double> D(m_nsegs);
	vector<double> Dup(m_nsegs);
	vector<double> Dfwd(m_nsegs);
	// identify blade
	MobilizedBodyIndex mbi = mobod.getMobilizedBodyIndex();
	if (mbi_of_first_blade == -1)
		mbi_of_first_blade = mbi;
	int blade_num = mbi - mbi_of_first_blade + 1;
	// calculate relative wind at blade tip and root in G
	// there MUST be an easier way to do this!
//	double bladeLen = m_tipR - m_rootR;
	double segLen = m_bladeLenX / m_nsegs;
	// get velocities of centers of root and tip sections in G
	Vec3 tipStationVelG = mobod.findStationVelocityInGround(state, Vec3(m_bladeLenX /2-segLen/2, 0, 0));
	Vec3 rootStationVelG = mobod.findStationVelocityInGround(state, Vec3(-m_bladeLenX /2+segLen/2, 0, 0));
	Vec3 angVelG = mobod.getBodyAngularVelocity(state);
	double angVel = angVelG[2];
	// now make the velocities relative to the wind
	tipStationVelG = tipStationVelG - windVelG;
	rootStationVelG = rootStationVelG - windVelG;
	// now measure velocities along the blade Y and Z axes
	// first transform unit vectors along axes to G
	Vec3 bodyUnitYG = mobod.expressVectorInGroundFrame(state, Vec3(0, 1, 0)); // blade Y unit vector
	Vec3 bodyUnitZG = mobod.expressVectorInGroundFrame(state, Vec3(0, 0, 1)); // blade Z unit vector
	// take dot products of unit vectors and velocities
	double rootVelY = ~rootStationVelG * bodyUnitYG;	// airspeed along Y at root
	double rootVelZ = ~rootStationVelG * bodyUnitZG;    // airspeed along Z at root
	double tipVelY = ~tipStationVelG * bodyUnitYG;		// airspeed along Y at tip
	double tipVelZ = ~tipStationVelG * bodyUnitZG;		// airspeed along Z at tip
	for (int iseg = 0; iseg < m_nsegs; iseg++)
	{
		double x = m_seg_x[iseg];
		double fwdSpeed = rootVelY + iseg * (tipVelY-rootVelY) / (m_nsegs - 1);
		double upSpeed = rootVelZ + iseg * (tipVelZ - rootVelZ) / (m_nsegs - 1);
		double v_squared = fwdSpeed*fwdSpeed + upSpeed*upSpeed;
		double v = sqrt(v_squared);	// speed of relative wind
		if (v < 0.00001)
			v = 0.00001;
		double v_angle = atan2(-upSpeed, fwdSpeed) * 180 / PI; // airflow angle relative to Y
		double rey = m_rey_coef * fwdSpeed;
		double alpha = v_angle + m_pitch;		// airfoil angle of attack
		double cl = m_af->getClForAlpha(alpha);
		double cd = m_af->getCdForAlpha(alpha);
		double lift = 0.5*rho*v_squared*m_segArea*cl;
		double drag = 0.5*rho*v_squared*m_segArea*cd;
		double lift_up = lift * fwdSpeed / v;	// upward component of airfoil lift
		double lift_fwd = lift * -upSpeed / v;	// forward component of airfoil lift
		double drag_up = drag * -upSpeed / v;	// upward component of airfoil drag
		double drag_fwd = -drag * fwdSpeed / v; // forward component of airfoil drag
		netLift += lift_up + drag_up;
		netLiftMoment += (lift_up + drag_up) * x;
		netThrust += lift_fwd + drag_fwd;
		netThrustMoment += (lift_fwd + drag_fwd) * x;
		// save values for printing
		if (printLevel > 1)
		{
			x_seg[iseg] = x;
			spd[iseg] = fwdSpeed;
			aoa[iseg] = alpha;
			CL[iseg] = cl;
			CD[iseg] = cd;
			L[iseg] = lift;
			Lup[iseg] = lift_up;
			Lfwd[iseg] = lift_fwd;
			D[iseg] = drag;
			Dup[iseg] = drag_up;
			Dfwd[iseg] = drag_fwd;
		}
	}
	// return values
	double FZdiff = netLiftMoment / m_bladeLenX;
	FZroot = netLift / 2 - FZdiff;
	FZtip = netLift / 2 + FZdiff;
	double FYdiff = netThrustMoment / m_bladeLenX; //
	FYroot = netThrust / 2 - FYdiff;
	FYtip = netThrust / 2 + FYdiff;
	// now get forces in G
	Vec3 vecLiftG = mobod.expressVectorInGroundFrame(state, Vec3(0, 0, netLift));
	Vec3 vecThrustG = mobod.expressVectorInGroundFrame(state, Vec3(0, 0, netThrust));
	// apply to blade body at correct stations
//	mobod.applyForceToBodyPoint(state, Vec3(m_))
	// print
	Real t = state.getTime();
	if (printLevel)
	{ 
		if (blade_num == 1)
			printf("\r\nblades: t %.4f\r\n", t);
		printf("blade %d: UZ (%6.3f %6.3f %6.3f) VELOCITY(Y,Z) root(%6.3f,%6.3f) tip(%6.3f,%6.3f) FORCE(Y,Z) root(%6.3f,%6.3f) tip(%6.3f,%6.3f)\r\n",
			blade_num, bodyUnitZG[0], bodyUnitZG[1], bodyUnitZG[2], rootVelY, rootVelZ, tipVelY, tipVelZ, FYroot, FZroot,
			FYtip, FZtip);
	}
	if (printLevel > 1)
	{
		std::printf("rotor blade data: angVel %.5f, rpm %.5f\r\n", angVel, angVel * 60 / (2 * PI));
		std::printf("fspd: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", spd[iseg]);
		std::printf("\r\naoa:  ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", aoa[iseg]);
		std::printf("\r\nCL:   ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", CL[iseg]);
		std::printf("\r\nCD:   ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", CD[iseg]);
		std::printf("\r\nL:    ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", L[iseg]);
		std::printf("\r\nLup:  ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Lup[iseg]);
		std::printf("\r\nLfwd: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Lfwd[iseg]);
		std::printf("\r\nD:    ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", D[iseg]);
		std::printf("\r\nDup:  ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Dup[iseg]);
		std::printf("\r\nDfwd: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Dfwd[iseg]);
		std::printf("\r\nLift: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Lup[iseg] + Dup[iseg]);
		std::printf("\r\nLiftM:");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", (Lup[iseg] + Dup[iseg])*x_seg[iseg]);
		std::printf("\r\nThr:  ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Dfwd[iseg] + Lfwd[iseg]);
		std::printf("\r\nThrM: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", (Dfwd[iseg] + Lfwd[iseg])*x_seg[iseg]);
		std::printf("\r\n");
	}
};


