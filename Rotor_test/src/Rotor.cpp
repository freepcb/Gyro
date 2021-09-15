// Implementation of classes for handling rotor aerodynamics

#include <vector>
#include "Rotor.h"
using namespace std;

#define PI 3.14159265358979323846

RotorBlade::RotorBlade(Airfoil* af, double root_r, double tip_r, double chord_len,
	double pitch, int nsegs)
{
	m_af = af;
	m_root_r = root_r;
	m_tip_r = tip_r;
	m_chord_len = chord_len;
	m_pitch = pitch;
	m_nsegs = nsegs;
	double bladeLen = m_tip_r - m_root_r;
	double segLen = bladeLen / nsegs;
	m_seg_area = m_chord_len * segLen;
	for (int iseg = 0; iseg < nsegs; iseg++)
	{
		m_seg_r.push_back(iseg*segLen + m_root_r + 0.5*segLen); // center of segment to axis 
	}
};

void RotorBlade::getForces(double angVel, double vertSpeed, 
	int printLevel, double& lift, double& torque)
{
	const double RHO = 1.204;     // density of dry air in kg / m ^ 3
	double netLift = 0.0;
	double netTorque = 0.0;
	// create vectors for printing
	vector<double> spd;
	vector<double> aoa;
	vector<double> CL;
	vector<double> CD;
	vector<double> L;
	vector<double> Lup;
	vector<double> Lfwd;
	vector<double> D;
	vector<double> Dup;
	vector<double> Dfwd;
	for (int iseg = 0; iseg < m_nsegs; iseg++)
	{
		double fwdSpeed = m_seg_r[iseg] * angVel;  // forward speed of this segment
		double v_squared = fwdSpeed*fwdSpeed + vertSpeed*vertSpeed;
		double v = sqrt(v_squared);	// speed of relative wind
		double v_angle = atan2(-vertSpeed, fwdSpeed) * 180 / PI;
		double alpha = v_angle + m_pitch;		// airfoil angle of attack
		double cl = m_af->getClForAlpha(alpha);
		double cd = m_af->getCdForAlpha(alpha);
		double lift = 0.5*RHO*v_squared*m_seg_area*cl;
		double drag = 0.5*RHO*v_squared*m_seg_area*cd;
		double lift_up = lift * fwdSpeed / v; // upward component of airfoil lift
		double lift_fwd = lift * -vertSpeed / v; // forward component of airfoil lift
		double drag_up = drag * -vertSpeed / v; // upward component of airfoil drag
		double drag_fwd = -drag * fwdSpeed / v; // forward component of airfoil drag
		netLift += lift_up + drag_up;
		netTorque += (lift_fwd + drag_fwd) * m_seg_r[iseg];
		// save values for printing
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


