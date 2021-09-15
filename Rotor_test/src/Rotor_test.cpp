#include "Simbody.h"
#include "Rotor.h"
#include <iostream>

//using namespace SimTK;

#define PI 3.14159265358979323846

int main()
{
	// world parameters
	const double gravity = -9.8;	// acceleration due to gravity
	const double airDensity = 1.204;	

	// rotor parameters
	Airfoil af;	  // default SG6042
	const double rotorMOI = 0.2;
	const int numBlades = 4;
	const int bladeNsegs = 18;	// number of segments to approximate spanwise variation
	const double bladeRootR = 0.1;  // distance of blade root from hub axis
	const double bladeTipR = 1.0;    // distance of blade tip from hub axis
	const double bladeChord = 0.055; // blade chord
	const double bladePitch = -6.0;  // blade pitch (degrees)
	double bladeLift;				// net upward lift from blade
	double bladeTorque;				// net torque from blade
	RotorBlade bl(&af, bladeRootR, bladeTipR, bladeChord, bladePitch, bladeNsegs);

	// airframe
	double airframeMass = 0.5;

	// loop variables
	const double timeStep = 0.01;
	double time = 0.0;
	double vertSpeed = -0.0001;
	double angVel = 0.0;
	double altitude = 0.0;
	int istep = 0;
	int printLevel = 0;
	while (true)
	{
		if((istep % 10) == 0)
			printLevel = 1;
		else
			printLevel = 0;;
		bl.getForces(angVel, vertSpeed, printLevel, bladeLift, bladeTorque);
		// multiply by number of blades
		bladeLift = bladeLift * numBlades;
		bladeTorque = bladeTorque * numBlades;
		angVel = angVel + timeStep * bladeTorque/rotorMOI;
		vertSpeed = vertSpeed + timeStep * (bladeLift/airframeMass + gravity);
		time += timeStep;
		if (printLevel)
		{
			double vertDrop = -0.5*gravity*pow((vertSpeed/gravity),2);
			printf("t %.3f, rpm %6.2f, vertSpeed %6.3f, vertDrop %6.3f\r\n\r\n",
				time, angVel * 60 / (2 * PI), vertSpeed, vertDrop);
		}
		istep++;
	}
}
