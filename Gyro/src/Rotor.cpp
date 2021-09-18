// Implementation of classes for handling rotor aerodynamics

#include <vector>
#include "Rotor.h"
using namespace std;

#define PI 3.14159265358979323846

// note that pitch is in degrees, not radians
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
	{
		m_seg_r.push_back(iseg*segLen + m_rootR + 0.5*segLen); // center of segment
	}
};

void RotorBlade::getForces(double angVel, double vertSpeed, int printLevel, 
	double& lift, double& torque)
{
	const double rho = 1.204;     // density of dry air in kg / m ^ 3
	double netLift = 0.0;
	double netTorque = 0.0;
	// create vectors for printing
	vector<double> spd;		spd.reserve(m_nsegs);
	vector<double> aoa;		aoa.reserve(m_nsegs);
	vector<double> CL;		CL.reserve(m_nsegs);
	vector<double> CD;		CD.reserve(m_nsegs);
	vector<double> L;		L.reserve(m_nsegs);
	vector<double> Lup;		Lup.reserve(m_nsegs);
	vector<double> Lfwd;	Lfwd.reserve(m_nsegs);
	vector<double> D;		D.reserve(m_nsegs);
	vector<double> Dup;		Dup.reserve(m_nsegs);
	vector<double> Dfwd;	Dfwd.reserve(m_nsegs);
	for (int iseg = 0; iseg < m_nsegs; iseg++)
	{
		double fwdSpeed = m_seg_r[iseg] * angVel;  // forward speed of this segment
		double v_squared = fwdSpeed*fwdSpeed + vertSpeed*vertSpeed;
		double v = sqrt(v_squared);	// speed of relative wind
		double v_angle = atan2(-vertSpeed, fwdSpeed) * 180 / PI; // degrees
		double alpha = v_angle + m_pitch;		// airfoil angle of attack
		double cl = m_af->getClForAlpha(alpha);
		double cd = m_af->getCdForAlpha(alpha);
		double lift = 0.5*rho*v_squared*m_segArea*cl;
		double drag = 0.5*rho*v_squared*m_segArea*cd;
		double lift_up = lift * fwdSpeed / v; // upward component of airfoil lift
		double lift_fwd = lift * -vertSpeed / v; // forward component of airfoil lift
		double drag_up = drag * -vertSpeed / v; // upward component of airfoil drag
		double drag_fwd = -drag * fwdSpeed / v; // forward component of airfoil drag
		netLift += lift_up + drag_up;
		netTorque += (lift_fwd + drag_fwd) * m_seg_r[iseg];
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
		printf("fspd:");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.2f ", spd[iseg]);
		printf("\r\naoa: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.2f ", aoa[iseg]);
		printf("\r\nCL:  ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.2f ", CL[iseg]);
		printf("\r\nCD:  ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.2f ", CD[iseg]);
		printf("\r\nL:   ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.3f ", L[iseg]);
		printf("\r\nLup: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.3f ", Lup[iseg]);
		printf("\r\nLfwd:");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.3f ", Lfwd[iseg]);
		printf("\r\nD:   ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.3f ", D[iseg]);
		printf("\r\nDup: ");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.3f ", Dup[iseg]);
		printf("\r\nDfwd:");
		for (int iseg = 0; iseg < m_nsegs; iseg++)
			printf("%6.3f ", Dfwd[iseg]);
		printf("\r\n");
	}
	// return values
	lift = netLift;
	torque = netTorque;
	if (printLevel > 0)
		printf("net lift: %.6f, net torque: %.6f\r\n", netLift, netTorque);
};


