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
	// get position of segment centers from hub axis and blade center
	for (int iseg = 0; iseg < nsegs; iseg++)
	{
		double x_from_root = iseg * segLen + 0.5*segLen;	 // x distance of segment center from root
		m_seg_x.push_back(x_from_root - (bladeLen / 2));	 // x distance from blade center
		m_seg_x_from_hub_Z.push_back(m_rootR + x_from_root); // x distance from hub axis
	}
	m_rey_coef = 70698 * m_chordLen;	// Reynolds number for v = 1 m/sec
};

// function getForces() 
// calculates aerodynamic forces on a rotor blade rotating around hub axis (Z-axis)
// divides the blade into m_nsegs sections to approximate continuous variation along blade
// then outputs net forces and net moments
// inputs:
//   windVelG = the wind velocity in the Ground frame
//   printLevel = 0 for none, 1 to print outputs, >1 to print details
// outputs:
//   lift = net force along Z-axis of blade (ie. "up")
//   rLift = radius of center of lift
//   thrust = net force acting along Y-axis of blade (ie. "fwd")
//   rThrust = radius of center of thrust
//
void RotorBlade::getForces(const State& state, const MobilizedBody& mobod, Vec3 windVelG, 
		int printLevel, double& lift, double& rLift, double& thrust, double& rThrust)
{
	const double rho = 1.204;     // density of dry air in kg / m ^ 3
	double netLift = 0.0;
	double netLiftMoment = 0.0;
	double netThrust = 0.0;
	double netThrustMoment = 0.0;
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
	// calculate relative wind at blade tip and root in G
	// there MUST be an easier way to do this!
	double bladeLen = m_tipR - m_rootR;
	double segLen = bladeLen / m_nsegs;
	// get velocities of centers of root and tip sections in G
	Vec3 tipStationVelG = mobod.findStationVelocityInGround(state, Vec3(bladeLen/2-segLen/2, 0, 0));
	Vec3 rootStationVelG = mobod.findStationVelocityInGround(state, Vec3(-bladeLen/2+segLen/2, 0, 0));
	Vec3 angVelG = mobod.getBodyAngularVelocity(state);
	double angVel = angVelG[2];
	// now make the velocities relative to the wind
	tipStationVelG = tipStationVelG - windVelG;
	rootStationVelG = rootStationVelG - windVelG;
	// now measure velocities along the blade Y and Z axes
	// first transform the axes to G
	Vec3 bodyUnitYG = mobod.expressVectorInGroundFrame(state, Vec3(0, 1, 0)); // blade Y unit vector
	Vec3 bodyUnitZG = mobod.expressVectorInGroundFrame(state, Vec3(0, 0, 1)); // blade Z unit vector
	// these are dot products
	double rootVelY = ~rootStationVelG * bodyUnitYG;	// airspeed along Y at root
	double rootVelZ = ~rootStationVelG * bodyUnitZG;    // airspeed along Z at root
	double tipVelY = ~tipStationVelG * bodyUnitYG;		// airspeed along Y at tip
	double tipVelZ = ~tipStationVelG * bodyUnitZG;		// airspeed along Z at tip
	for (int iseg = 0; iseg < m_nsegs; iseg++)
	{
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
		double lift_up = lift * fwdSpeed / v;	 // upward component of airfoil lift
		double lift_fwd = lift * -upSpeed / v; // forward component of airfoil lift
		double drag_up = drag * -upSpeed / v;  // upward component of airfoil drag
		double drag_fwd = -drag * fwdSpeed / v;  // forward component of airfoil drag
		netLift += lift_up + drag_up;
		netLiftMoment += (lift_up + drag_up) * m_seg_x_from_hub_Z[iseg];
		netThrust += lift_fwd + drag_fwd;
		netThrustMoment += (lift_fwd + drag_fwd) * m_seg_x_from_hub_Z[iseg];
		// save values for printing
		if (printLevel > 1)
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
	// return values
	lift = netLift;
	if (fabs(netLift) > 0.00001)
		rLift = netLiftMoment / netLift;
	else
		rLift = bladeLen / 2;
	thrust = netThrust;
	if (fabs(netThrust) > 0.00001)
		rThrust = netThrustMoment / netThrust;
	else
		rThrust = bladeLen / 2;
	// now get forces in G
	Vec3 vecLiftG = mobod.expressVectorInGroundFrame(state, Vec3(0, 0, netLift));
	Vec3 vecThrustG = mobod.expressVectorInGroundFrame(state, Vec3(0, 0, netThrust));
	// apply to blade body at correct stations
//	mobod.applyForceToBodyPoint(state, Vec3(m_))
	// print
	if (printLevel)
	{
		printf("blade:UZ %5.3f %5.3f %5.3f root:vY %6.2f vZ %6.2f tip:vY %6.2f vZ %6.2f net:L %6.3f rL %5.2f T %6.3f rT %6.3f\r\n",
			bodyUnitZG[0], bodyUnitZG[1], bodyUnitZG[2], rootVelY, rootVelZ, tipVelY, tipVelZ,
			lift, rLift, thrust, rThrust);
	}
	if (printLevel > 1)
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
};


