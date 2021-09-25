// definition of classes for rotor aerodynamics
#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "Simbody.h"

using namespace std;
using namespace SimTK;

// Airfoil polar data, currently just for SG6042  
class Airfoil
{
public :
	Airfoil() {};
	~Airfoil() {};
private:
	string m_name_str = "SG6042";
	int m_arrayStart = -7;	// arrays are from -7 to 180 degrees
	double m_cl[188] = {
		-0.385, -0.432, -0.442, -0.346, -0.250, -0.134, -0.0346, 0.0455, 0.127, 0.228,
		0.402, 0.557, 0.709, 0.903, 1.123, 1.208, 1.289, 1.361, 1.3891, 1.2867,
		1.115, 0.943, 0.771, 0.737, 0.703, 0.668, 0.634, 0.6, 0.614, 0.628,
		0.642, 0.656, 0.67, 0.698, 0.726, 0.754, 0.782, 0.81, 0.844, 0.878,
		0.912, 0.946, 0.98, 0.988, 0.996, 1.004, 1.012, 1.02, 1.024, 1.028,
		1.032, 1.036, 1.04, 1.032, 1.024, 1.016, 1.008, 1, 0.99, 0.98,
		0.97, 0.96, 0.95, 0.934, 0.918, 0.902, 0.886, 0.87, 0.846, 0.822,
		0.798, 0.774,0.75, 0.724, 0.698, 0.672, 0.646, 0.62, 0.596, 0.572,
		0.548, 0.524, 0.5, 0.472, 0.444, 0.416, 0.388, 0.36, 0.33, 0.3,
		0.27, 0.24, 0.21, 0.186, 0.162, 0.138, 0.114, 0.09, 0.064, 0.038,
		0.012, -0.014, -0.04, -0.066, -0.092, -0.118, -0.144, -0.17, -0.194, -0.218,
		-0.242, -0.266, -0.29, -0.316, -0.342, -0.368, -0.394, -0.42, -0.45, -0.48,
		-0.51, -0.54, -0.57, -0.592, -0.614, -0.636, -0.658, -0.68, -0.7, -0.72,
		-0.74, -0.76, -0.78, -0.798, -0.816, -0.834, -0.852, -0.87, -0.888, -0.906,
		-0.924, -0.942, -0.96, -0.964, -0.968, -0.972, -0.976, -0.98, -0.964, -0.948,
		-0.932, -0.916, -0.9, -0.874, -0.848, -0.822, -0.796, -0.77, -0.748, -0.726,
		-0.704, -0.682, -0.66, -0.654, -0.648, -0.642, -0.636, -0.63, -0.64, -0.65,
		-0.66, -0.67, -0.68, -0.72, -0.76, -0.8, -0.84, -0.88000, -0.838, -0.796,
		-0.754, -0.712, -0.67000, -0.536, -0.402, -0.268, -0.134, 0.0
	};
	double m_cd[188] = {
		0.100, 0.093, 0.083, 0.063, 0.042, 0.033, 0.029, 0.026, 0.029, 0.033,
		0.039, 0.042, 0.044, 0.038, 0.029, 0.033, 0.040, 0.050, 0.062, 0.079,
		0.111, 0.144, 0.177, 0.210, 0.222, 0.245, 0.267, 0.290, 0.312, 0.334,
		0.356, 0.378, 0.400, 0.430, 0.460, 0.490, 0.520, 0.550, 0.582, 0.614,
		0.646, 0.678, 0.710, 0.748, 0.786, 0.824, 0.862, 0.900, 0.940, 0.980,
		1.020, 1.060, 1.100, 1.130, 1.160, 1.190, 1.220, 1.250, 1.280, 1.310,
		1.340, 1.370, 1.400, 1.420, 1.440, 1.460, 1.480, 1.500, 1.520, 1.540,
		1.560, 1.580, 1.600, 1.610, 1.620, 1.630, 1.640, 1.650, 1.680, 1.710,
		1.740, 1.770, 1.800, 1.790, 1.780, 1.770, 1.760, 1.750, 1.754, 1.758,
		1.762, 1.766, 1.770, 1.776, 1.782, 1.788, 1.794, 1.800, 1.794, 1.788,
		1.782, 1.776, 1.770, 1.766, 1.762, 1.758, 1.754, 1.750, 1.740, 1.730,
		1.720, 1.710, 1.700, 1.690, 1.680, 1.670, 1.660, 1.650, 1.640, 1.630,
		1.620, 1.610, 1.600, 1.590, 1.580, 1.570, 1.560, 1.550, 1.524, 1.498,
		1.472, 1.446, 1.420, 1.392, 1.364, 1.336, 1.308, 1.280, 1.244, 1.208,
		1.172, 1.136, 1.100, 1.052, 1.004, 0.956, 0.908, 0.860, 0.840, 0.820,
		0.800, 0.780, 0.760, 0.728, 0.696, 0.664, 0.632, 0.600, 0.560, 0.520,
		0.480, 0.440, 0.400, 0.380, 0.360, 0.340, 0.320, 0.300, 0.286, 0.272,
		0.258, 0.244, 0.230, 0.208, 0.186, 0.164, 0.142, 0.120, 0.104, 0.088,
		0.072, 0.056, 0.040, 0.036, 0.032, 0.028, 0.024, 0.020
	};
public:
	double getClForAlpha(double& alpha) 
	{
		int index = int(round(alpha)) - m_arrayStart;
		double cl = m_cl[index];
		return cl; 
	}
	double getCdForAlpha(double& alpha) 
	{
		int index = int(round(alpha)) - m_arrayStart;
		double cd = m_cd[index];
		return cd;
	}

};

// class to represent a single rotor blade attached to a hub
// assumes that the hub rotates around a vertical axis
// frame of reference:
//   origin = center of hub
//   z-axis = axis of rotation ("Up")
//   y-axis = direction of blade motion ("Forward")
//   x-axis = long axis of blade ("Radial")
// constructor arguments:
//   af = airfoil
//   rootR = radius of root of blade from axis
//   tipR = radius of tip of blade from axis
//   chordLen = chord length
//   pitch = pitch of the chord relative to the plane of the hub (degrees)
//   nsegs = blade is segmented to approximate variation of alpha along length
// if a flap hinge is added, then the Z-axis will no longer be the hub Z-axis

class RotorBlade
{
public:
	RotorBlade(Airfoil* af, double rootR, double tipR, double chordLen,
		double pitch, int nsegs);
	~RotorBlade() {};

	// for getForces():
	//   state = Simbody state
	//	 mobod = mobilized body for blade
	//   windVel is the vector of the airflow in the Ground frame
	//   angVel is the angular velocity of the blade
	void getForces(const State& state, const MobilizedBody& mobod, Vec3 windVelG, 
			int printLevel, double& lift, double& torque);
private:
	Airfoil* m_af;
	double m_rootR;
	double m_tipR;
	double m_chordLen;
	double m_pitch;
	int m_nsegs;
	double m_segArea;
	double m_angVel;
	double m_vertSpeed;
	int m_printLevel;
	vector<double> m_seg_r;	
};