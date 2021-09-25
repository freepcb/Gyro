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
	m_rootR = rootR;
	m_tipR = tipR;
	m_chordLen = chordLen;
	m_pitch = pitch;
	m_nsegs = nsegs;
	double bladeLen = m_tipR - m_rootR;
	double segLen = bladeLen / nsegs;
	m_segArea = m_chordLen * segLen;
	for (int iseg = 0; iseg < nsegs; iseg++)	
		m_seg_r.push_back(iseg*segLen + m_rootR + 0.5*segLen); // radii of centers of segments	}
};

// Note: windVelG is the wind velocity in the Ground frame
// However, for convenience I sometimes don't allow motion of the rotor on the Z axis, and instead
// simulate the descent by setting windVelG to represent an upward airflow at the calculated descent
// velocity, so that the aerodynamic forces are the same 
void RotorBlade::getForces(const State& state, const MobilizedBody& mobod, Vec3 windVelG, 
		int printLevel, double& lift, double& torque)
{
	const double rho = 1.204;     // density of dry air in kg / m ^ 3
	double netLift = 0.0;
	double netTorque = 0.0;
	// create vectors for printing
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
	// calculate relative wind at blade tip and root in local frame
	// there MUST be an easier way to do this!
	double bladeLen = m_tipR - m_rootR;
	double segLen = bladeLen / m_nsegs;
	// get velocities of root and tip sections in G
	Vec3 tipStationVelG = mobod.findStationVelocityInGround(state, Vec3(bladeLen/2-segLen/2, 0, 0));
	Vec3 rootStationVelG = mobod.findStationVelocityInGround(state, Vec3(-bladeLen/2+segLen/2, 0, 0));
	Vec3 angVelG = mobod.getBodyAngularVelocity(state);
	double angVel = angVelG[2];
	// now make the velocities relative to the wind, not G
	tipStationVelG = tipStationVelG - windVelG;
	rootStationVelG = rootStationVelG - windVelG;
	// now measure components along the blade Y and Z axes, first transform the axes to G
	Vec3 bodyUnitYG = mobod.expressVectorInGroundFrame(state, Vec3(0, 1, 0)); // blade Y unit vector
	Vec3 bodyUnitZG = mobod.expressVectorInGroundFrame(state, Vec3(0, 0, 1)); // blade Z unit vector
	// these are dot products
	double tipVelY = ~tipStationVelG * bodyUnitYG;		// airspeed in Y
	double tipVelZ = ~tipStationVelG * bodyUnitZG;		// airspeed in Z
	double rootVelY = ~rootStationVelG * bodyUnitYG;	
	double rootVelZ = ~rootStationVelG * bodyUnitZG;	
	for (int iseg = 0; iseg < m_nsegs; iseg++)
	{
		double vertSpeed = rootVelZ + iseg * (tipVelZ - rootVelZ) / (m_nsegs - 1);
		double fwdSpeed = rootVelY + iseg * (tipVelY-rootVelY) / (m_nsegs - 1);
		double v_squared = fwdSpeed*fwdSpeed + vertSpeed*vertSpeed;
		double v = sqrt(v_squared);	// speed of relative wind
		if (v < 0.00001)
			v = 0.00001;
		double v_angle = atan2(vertSpeed, fwdSpeed) * 180 / PI; // airflow angle relative to Y
		double alpha = v_angle + m_pitch;		// airfoil angle of attack
		double cl = m_af->getClForAlpha(alpha);
		double cd = m_af->getCdForAlpha(alpha);
		double lift = 0.5*rho*v_squared*m_segArea*cl;
		double drag = 0.5*rho*v_squared*m_segArea*cd;
		double lift_up = lift * fwdSpeed / v;	// upward component of airfoil lift
		double lift_fwd = lift * -vertSpeed / v; // forward component of airfoil lift
		double drag_up = drag * -vertSpeed / v; // upward component of airfoil drag
		double drag_fwd = -drag * fwdSpeed / v; // forward component of airfoil drag
		netLift += lift_up + drag_up;
		netTorque += (lift_fwd + drag_fwd) * m_seg_r[iseg];
//		if (netLift > 10 || vertSpeed > 0.0)
//			std::printf("********** Error ************\r\n");
		// save values for printing
		if (printLevel > 0)
		{
			spd.push_back(fwdSpeed);
			aoa.push_back(alpha);
			CL.push_back(cl);
			CD.push_back(cd);
			L.push_back(lift);
			Lup.push_back(lift_up);
			Lfwd.push_back(lift_fwd);
			D.push_back(drag);
			Dup.push_back(drag_up);
			Dfwd.push_back(drag_fwd);
		}
	}
	// print
	if (printLevel > 0)
	{
		std::printf("rotor blade data: angVel %.5f, rpm %.5f, vertSpeed %.5f\r\n", angVel, angVel * 60 / (2 * PI), windVelG[2]);
		std::printf("fspd:");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", spd[iseg]);
		std::printf("\r\naoa: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", aoa[iseg]);
		std::printf("\r\nCL:  ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", CL[iseg]);
		std::printf("\r\nCD:  ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", CD[iseg]);
		std::printf("\r\nL:   ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", L[iseg]);
		std::printf("\r\nLup: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Lup[iseg]);
		std::printf("\r\nLfwd:");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Lfwd[iseg]);
		std::printf("\r\nD:   ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", D[iseg]);
		std::printf("\r\nDup: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Dup[iseg]);
		std::printf("\r\nDfwd:");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			std::printf("%8.4f ", Dfwd[iseg]);
		std::printf("\r\n");
	}
	// return values
	lift = netLift;
	torque = netTorque;
	if (printLevel > 0)
		std::printf("net lift: %.6f, net torque: %.6f\r\n", netLift, netTorque);
};


