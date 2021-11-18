#include "Simbody.h"
#include <iostream>
#include "Rotor.h"
#include "Rocket.h"

using namespace SimTK; 

// defines
#define USE_ROTOR

// global variables for convenience and laziness
// indices for mobilized bodies so they can be found by name
MobilizedBodyIndex g_bi_motor;
MobilizedBodyIndex g_bi_mast;
MobilizedBodyIndex g_bi_hub;
MobilizedBodyIndex g_bi_blade1;
MobilizedBodyIndex g_bi_blade2;
MobilizedBodyIndex g_bi_blade3;
MobilizedBodyIndex g_bi_blade4;
MobilizedBodyIndex g_bi_fin1;
MobilizedBodyIndex g_bi_fin2;
MobilizedBodyIndex g_bi_fin3;
MobilizedBodyIndex g_bi_fin4;

// dimensions of rotor hub
const Real g_hubDensity = 1000.0;	// density of rotor blades in kg/m^3
const Real g_hubR = 0.1;			// hub radius, where blades are attached
const Real g_hubThick = 0.01;		// hub thickness
const Real g_hubM = g_hubDensity * g_hubThick * 2 * Pi * g_hubR * g_hubR; // mass

// dimensions of rotor blades
const Real g_bladeDensity = 1000.0;	// density of rotor blades in kg/m^3
const Real g_bladeLen = 1.8;		// length of blade
const Real g_bladeRootR = g_hubR;	// distance of blade root from hub axis
const Real g_bladeTipR = g_bladeLen + g_hubR;  // distance of blade tip from hub axis
const Real g_bladeChord = 0.2;		// chord length of airfoil
const Real g_bladeThick = 0.005;	// average thickness of blade with airfoil
const Real g_bladeM = g_bladeDensity* g_bladeLen * g_bladeChord * g_bladeThick; // blade mass
const int g_numBlades = 4;			// number of blades in rotor
const int g_bladeNsegs = 18;		// number of airfoil segments per blades 
const int g_bladeFlapAngleDeg = 30;	// fixed blade flap angle (degrees) 

// dimensions of fins
const Real finDensity = 1000.0;
const Real finThick = 0.003;		// thickness = 3 mm
const Real finLen = 0.5;		    // length (spanwise)
const Real finChord = 0.5;		    // chord length
const Real finM = finDensity * finLen * finChord * finThick;	// kg

// dimensions of motor (ie. fuselage)
const Real motorLen = 2;	// 2 m
const Real motorR = 0.1;	// 10 cm
const Real motorM = 5;	// 0.5 kg

// dimensions of cylindrical rotor mast
const Real mastDensity = 1000;		// 2 m
const Real mastLen = 3.0;		// 3 m
const Real mastR = 0.01;	// 10 cm
const Real mastM = mastDensity * mastLen * 2 * Pi * mastR * mastR;

// distance from tail to center of mass of motor and fins
Real g_motor_cm_Z = (motorM*motorLen / 2 + 4 * finM*finChord / 2) / (motorM + 4 * finM);

// physical and simulation constants
const Real g_gravity = -9.8;	// acceleration due to gravity in kg/sec^2
const Real g_mass = 1.5;	// airframe mass
Real g_last_time = 0.0;		// time of last update
Real g_lift = 0.0;			// rotor lift
Real g_torque = 0.0;		// rotor torque
Real g_angVel = 0.0;		// rotor ang vel
Real g_rpm = 0.0;			// rotor RPM
Real g_altitude = 0.0;		// rotor altitude;

// hub data
Vec3 g_hubCenter(0);	// hub center position
Vec3 g_hubLift;			// hub lift vector
Vec3 g_hubUZ;			// unit vector along hub Z axis
Rotation g_hubRotation;


// rotor objects and parameters
Airfoil af;						 // airfoil for blades
// distance of blade tip from hub axis
const double bladePitch = -6.0;  // blade pitch (degrees)
double bladeLift;				 // net upward lift from blade
double bladeTorque;				 // net torque from blade
// single blade object, since all blades are the same
RotorBlade bl(&af, g_bladeRootR, g_bladeTipR, g_bladeChord, bladePitch, g_bladeNsegs);
// single fin object
Fin fin(finDensity, finLen, finChord, finThick);


//==============================================================================
//                              SHOW DATA
//==============================================================================
// Generate text in the scene    
class ShowData : public DecorationGenerator {
public:
    explicit ShowData(const MultibodySystem& mbs) : m_mbs(mbs) {}
    void generateDecorations(const State&                state, 
                             Array_<DecorativeGeometry>& geometry) override;
private:
    const MultibodySystem& m_mbs;
};

//==============================================================================
//                         CUSTOM FORCE
//==============================================================================
class rotorForces : public Force::Custom::Implementation
{
public:
	rotorForces(SimbodyMatterSubsystem& matter) : matter(matter) {}

	void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
			Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
	{
		// calculate and apply normal forces on fins
		const MobilizedBody& fin1_mobody = matter.getMobilizedBody(g_bi_fin1);
		const MobilizedBody& fin2_mobody = matter.getMobilizedBody(g_bi_fin2);
		const MobilizedBody& fin3_mobody = matter.getMobilizedBody(g_bi_fin3);
		const MobilizedBody& fin4_mobody = matter.getMobilizedBody(g_bi_fin4);
		int finPrintLevel = 1;
		fin.calcForces(state, bodyForces, fin1_mobody, finPrintLevel, '1', true);
		fin.calcForces(state, bodyForces, fin2_mobody, finPrintLevel, '2', true);
		fin.calcForces(state, bodyForces, fin3_mobody, finPrintLevel, '3', true);
		fin.calcForces(state, bodyForces, fin4_mobody, finPrintLevel, '4', true);
#ifdef USE_ROTOR
		// get mobilized bodies for hub and blades
		const MobilizedBody& hub_mobody = matter.getMobilizedBody(g_bi_hub);
		g_hubRotation = hub_mobody.getBodyRotation(state);
		Vec3 uz = hub_mobody.expressVectorInGroundFrame(state, Vec3(0, 0, 1));
		g_hubUZ = uz;
		const MobilizedBody& blade1_mobody = matter.getMobilizedBody(g_bi_blade1);
		const MobilizedBody& blade2_mobody = matter.getMobilizedBody(g_bi_blade2);
		const MobilizedBody& blade3_mobody = matter.getMobilizedBody(g_bi_blade3);
		const MobilizedBody& blade4_mobody = matter.getMobilizedBody(g_bi_blade4);
		// print hub motion, same format as in event reporter
		SpatialVec vel = hub_mobody.getBodyVelocity(state);
//		printf("force: t %0.5f, hub: v (%.2f,%.2f,%.2f) RPM %.5f, alt %.2f\r\n",
//			state.getTime(), vel[1][0], vel[1][1], vel[1][2], vel[0][2] * 60 / (2 * Pi), g_altitude);
		g_angVel = vel[0][2];
		Vec3 angVel = vel[0];		// get ang velocity of hub
		Vec3 blade1XInG = blade1_mobody.expressVectorInGroundFrame(state, Vec3(1, 0, 0));
		double blade1Angle = (180 / Pi)*atan2(blade1XInG[1], blade1XInG[0]);
		Real t = state.getTime();
		// update aerodynamic forces
		// these variables will actually be set by bl.getForces()
		Real FZroot1 = 0; Real FZtip1 = 0; Real FYroot1 = 0; Real FYtip1 = 0;
		Real FZroot2 = 0; Real FZtip2 = 0; Real FYroot2 = 0; Real FYtip2 = 0;
		Real FZroot3 = 0; Real FZtip3 = 0; Real FYroot3 = 0; Real FYtip3 = 0;
		Real FZroot4 = 0; Real FZtip4 = 0; Real FYroot4 = 0; Real FYtip4 = 0;
		double xoff = bl.m_bladeLenX/2;
		int printLevel = 1;
		// set wind velocity to small updraft to avoid divide by zero or atan2() errors
		bl.getForces(state, blade1_mobody, Vec3(0, 0, 0.0001), printLevel, FZroot1, FZtip1, FYroot1, FYtip1);
		printLevel = 1;
		bl.getForces(state, blade2_mobody, Vec3(0, 0, 0.0001), printLevel, FZroot2, FZtip2, FYroot2, FYtip2);
		bl.getForces(state, blade3_mobody, Vec3(0, 0, 0.0001), printLevel, FZroot3, FZtip3, FYroot3, FYtip3);
		bl.getForces(state, blade4_mobody, Vec3(0, 0, 0.0001), printLevel, FZroot4, FZtip4, FYroot4, FYtip4);
		// calculate rotor thrust torque on the hub, just for information
		double torque1 = FYtip1 * (g_bladeLen + g_hubR) + FYroot1 * g_hubR;
		double torque2 = FYtip2 * (g_bladeLen + g_hubR) + FYroot2 * g_hubR;
		double torque3 = FYtip3 * (g_bladeLen + g_hubR) + FYroot3 * g_hubR;
		double torque4 = FYtip4 * (g_bladeLen + g_hubR) + FYroot4 * g_hubR;
		// calculate lift on the hub, just for information
		double lift1 = FZroot1 + FZtip1;
		double lift2 = FZroot2 + FZtip2;
		double lift3 = FZroot3 + FZtip3;
		double lift4 = FZroot4 + FZtip4;
		g_torque = torque1 + torque2 + torque3 + torque4;
		g_lift = lift1 + lift2 + lift3 + lift4;
		// get hub velocity and position
		Vec3 hubVel = vel[1];
		Vec3 hubPos = hub_mobody.getBodyOriginLocation(state);
		g_altitude = hubPos[2];
		if (printLevel)
			printf("forces: t %.4f, angle %5.1f, RPM %6.2f, alt %3.0f, ROD %5.2f, L: %.3f (%.3f+%.3f+%.3f+%.3f) T: %.3f (%.3f+%.3f+%.3f+%.3f)\r\n", 
				t, blade1Angle, angVel[2]*60/(2*Pi), hubPos[2], hubVel[2], 
				g_lift, lift1, lift2, lift3, lift4, g_torque, torque1, torque2, torque3, torque4);
		g_last_time = t;
		// apply lift and thrust to blades as external forces at the blade roots and tips
		// start by creating force vectors in blade frames
		Vec3 vecFZroot1(0, 0, FZroot1);
		Vec3 vecFZroot2(0, 0, FZroot2);
		Vec3 vecFZroot3(0, 0, FZroot3);
		Vec3 vecFZroot4(0, 0, FZroot4);
		Vec3 vecFZtip1(0, 0, FZtip1);
		Vec3 vecFZtip2(0, 0, FZtip2);
		Vec3 vecFZtip3(0, 0, FZtip3);
		Vec3 vecFZtip4(0, 0, FZtip4);
		Vec3 vecFYroot1(0, FYroot1, 0); 
		Vec3 vecFYroot2(0, FYroot2, 0); 
		Vec3 vecFYroot3(0, FYroot3, 0); 
		Vec3 vecFYroot4(0, FYroot4, 0);
		Vec3 vecFYtip1(0, FYtip1, 0); 
		Vec3 vecFYtip2(0, FYtip2, 0); 
		Vec3 vecFYtip3(0, FYtip3, 0); 
		Vec3 vecFYtip4(0, FYtip4, 0);
		// transform forces to G
		Vec3 vecFZrootInG1 = blade1_mobody.expressVectorInGroundFrame(state, vecFZroot1);
		Vec3 vecFZrootInG2 = blade2_mobody.expressVectorInGroundFrame(state, vecFZroot2);
		Vec3 vecFZrootInG3 = blade3_mobody.expressVectorInGroundFrame(state, vecFZroot3);
		Vec3 vecFZrootInG4 = blade4_mobody.expressVectorInGroundFrame(state, vecFZroot4);
		Vec3 vecFZtipInG1 = blade1_mobody.expressVectorInGroundFrame(state, vecFZtip1);
		Vec3 vecFZtipInG2 = blade2_mobody.expressVectorInGroundFrame(state, vecFZtip2);
		Vec3 vecFZtipInG3 = blade3_mobody.expressVectorInGroundFrame(state, vecFZtip3);
		Vec3 vecFZtipInG4 = blade4_mobody.expressVectorInGroundFrame(state, vecFZtip4);
		Vec3 vecFYrootInG1 = blade1_mobody.expressVectorInGroundFrame(state, vecFYroot1);
		Vec3 vecFYrootInG2 = blade2_mobody.expressVectorInGroundFrame(state, vecFYroot2);
		Vec3 vecFYrootInG3 = blade3_mobody.expressVectorInGroundFrame(state, vecFYroot3);
		Vec3 vecFYrootInG4 = blade4_mobody.expressVectorInGroundFrame(state, vecFYroot4);
		Vec3 vecFYtipInG1 = blade1_mobody.expressVectorInGroundFrame(state, vecFYtip1);
		Vec3 vecFYtipInG2 = blade2_mobody.expressVectorInGroundFrame(state, vecFYtip2);
		Vec3 vecFYtipInG3 = blade3_mobody.expressVectorInGroundFrame(state, vecFYtip3);
		Vec3 vecFYtipInG4 = blade4_mobody.expressVectorInGroundFrame(state, vecFYtip4);
		// apply forces to blade roots and tips
		blade1_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFZrootInG1, bodyForces);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFZrootInG2, bodyForces);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFZrootInG3, bodyForces);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFZrootInG4, bodyForces);
		blade1_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFZtipInG1, bodyForces);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFZtipInG2, bodyForces);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFZtipInG3, bodyForces);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFZtipInG4, bodyForces);
//#if 0
		blade1_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG1, bodyForces);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG2, bodyForces);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG3, bodyForces);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG4, bodyForces);
		blade1_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG1, bodyForces);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG2, bodyForces);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG3, bodyForces);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG4, bodyForces);
//#endif
		g_hubCenter = hub_mobody.findStationLocationInGround(state, Vec3(0));
		g_hubLift = vecFZrootInG1 + vecFZrootInG2 + vecFZrootInG3 + vecFZrootInG4
			+ vecFZtipInG1 + vecFZtipInG2 + vecFZtipInG3 + vecFZtipInG4;
#endif
	}
	Real calcPotentialEnergy(const State& state) const 
	{
		double energy = 0.0;
		return energy;
	}
	bool dependsOnlyOnPositions() const 
	{
		return false;
	}
	void SetWindInG(double windVel, double windAngle)
	{
		m_windVel = windVel;
		m_windAngle = windAngle;
	}
private:
	SimbodyMatterSubsystem& matter;
	double m_windVel;
	double m_windAngle;
};

//==============================================================================
//                           EVENT REPORTER
//==============================================================================
class MyEventReporter : public PeriodicEventReporter {
public:
	MyEventReporter(const MultibodySystem& system, const MobilizedBody& mobod,
		Real reportInterval)
		: PeriodicEventReporter(reportInterval), system(system), mobod(mobod) {}
	void handleEvent(const State& state) const {
		system.realize(state, Stage::Position);   // so that variables can be accessed
		SpatialVec vel = mobod.getBodyVelocity(state);
		SpatialVec velMobilizer = mobod.getMobilizerVelocity(state);
		Vec3 UZInG = mobod.expressVectorInGroundFrame(state, Vec3(0, 0, 1));
		g_rpm = velMobilizer[0][2] * 60 / (2 * Pi);
		printf("event:t %0.3f, hub:z (%5.2f,%5.2f,%5.2f) v (%5.2f,%5.2f,%5.2f) RPM %6.2f, alt %.2f, lift %5.2f, torque %5.3f\r\n",
			state.getTime(), UZInG[0], UZInG[1], UZInG[2], 
			vel[1][0], vel[1][1], vel[1][2], g_rpm, g_altitude,
			g_lift, g_torque);
	}
private:
	const MultibodySystem& system;
	const MobilizedBody& mobod;
};


//==============================================================================
//                                  MAIN
//==============================================================================
int main() {
  try 
  {   
    // Create the system.   
    MultibodySystem system; system.setUpDirection(ZAxis);
    SimbodyMatterSubsystem matter(system);
	GeneralForceSubsystem forces(system);
	Force::UniformGravity gravity(forces, matter, Vec3(0, 0, -9.8));

	// turn off drawing reference frames 
	matter.setShowDefaultGeometry(false); 

	// create cylindrical motor body with long axis = Z
	Inertia motorI = motorM * UnitInertia::cylinderAlongZ(motorR, motorLen/2); // MOI of motor 
	Body::Rigid motor_body(MassProperties(motorM, Vec3(0), motorI));
	// create decoration for the motor
	Rotation R;
	R.setRotationFromAngleAboutX(Pi / 2);
	Transform X(R);
	// note: DecorativeCylinder is oriented along Y so needs to be rotated
	// use half/length
	motor_body.addDecoration(X, DecorativeCylinder(motorR, motorLen/2).setColor(Green));
	// create MobilizedBody for the motor
	MobilizedBody::Free motor_mobody(matter.Ground(), Transform(),
		motor_body, Transform(Vec3(0,0, -motorLen/2)));	// joint at bottom of motor
	g_bi_motor = motor_mobody.getMobilizedBodyIndex();

	// create rotor mast
	Inertia mastI = mastM * UnitInertia::cylinderAlongZ(mastR, mastLen / 2); // MOI of mast 
	Body::Rigid mast_body(MassProperties(mastM, Vec3(0), mastI));
	// create decoration for the mast
	mast_body.addDecoration(X, DecorativeCylinder(mastR, mastLen / 2).setColor(Red));
	// F mobilizer frame
	Rotation R_mast_F;	// no rotation for pin joint with axis along Z
	R_mast_F.setRotationFromAngleAboutX(Pi / 2);
	Transform X_mast_F(R_mast_F, Vec3(0,0,0)); // at center of motor 
	// M mobilizer frame
	Rotation R_mast_M;
	R_mast_M.setRotationFromAngleAboutY(0); 
	// Sety angle of mast from motor, 0 = lateral, Pi/2 = behind, -Pi/2 = ahead 
	Rotation R_mast_angle;  
	R_mast_angle.setRotationFromAngleAboutZ(Pi/2);
	R_mast_M = R_mast_M * R_mast_angle;
	Transform X_mast_M(R_mast_M, Vec3(0, 0, -mastLen / 2));
	// create MobilizedBody for the mast
	MobilizedBody::Pin mast_mobody(motor_mobody, X_mast_F,
			mast_body, X_mast_M);	// joint at bottom of mast
	g_bi_mast = mast_mobody.getMobilizedBodyIndex();

	// create hub body
	Inertia hubI = g_hubM * UnitInertia::cylinderAlongZ(g_hubR, g_hubThick); // MOI of hub
	Body::Rigid hub_body(MassProperties(g_hubM, Vec3(0), hubI));
	printf("Hub: Ixx %.4f, Iyy %.4f, Izz %.4f\r\n",
		hubI.getMoments()[0],
		hubI.getMoments()[1],
		hubI.getMoments()[2]
		);
	// create decoration for the hub as a sphere, not a cylinder
	Rotation R_hub;
	R_hub.setRotationFromAngleAboutX(Pi/2);
	Transform X_hub(R_hub);
	hub_body.addDecoration(X_hub, DecorativeSphere(g_hubR).setColor(Red));
	// create MobilizedBody for the hub
	// with pin joint at top of mast and center of hub
	MobilizedBody::Pin hub_mobody(mast_mobody, Transform(Vec3(0, 0, mastLen / 2)),
			hub_body, Transform());
	g_bi_hub = hub_mobody.getMobilizedBodyIndex();

#ifdef USE_ROTOR
	// now create body for the blades
	const Vec3 bladeHalfLengths(g_bladeLen/2, g_bladeChord/2, g_bladeThick/2);
	Inertia bladeI = g_bladeM * UnitInertia::brick(bladeHalfLengths); // MOI of blade
	Body::Rigid blade_body(MassProperties(g_bladeM, Vec3(0), bladeI));
	printf("Blade: Ixx %.4f, Iyy %.4f, Izz %.4f\r\n",
		bladeI.getMoments()[0],
		bladeI.getMoments()[1],
		bladeI.getMoments()[2]
	);
	// create decoration for the blade
	blade_body.addDecoration(Transform(), DecorativeBrick(bladeHalfLengths).setColor(Blue));

	// create MobilizedBody for blade 1
	double flapAng = g_bladeFlapAngleDeg * Pi / 180;
	Transform X_bladeF1(Vec3(g_hubR, 0, 0));	// fixed mob frame at X edge of hub
	Rotation R_bladeF1;
	R_bladeF1.setRotationFromAngleAboutY(flapAng);	// set flap angle
	Transform X_bladeM1(R_bladeF1, Vec3(-g_bladeLen / 2, 0, 0));	// moving mob frame
	MobilizedBody::Weld blade1_mobody(hub_mobody, X_bladeF1,		// weld blade to hub
		blade_body, X_bladeM1);
	g_bi_blade1 = blade1_mobody.getMobilizedBodyIndex();

	// create MobilizedBody for blade 2
	Rotation R_bladeF2;
	R_bladeF2.setRotationFromAngleAboutZ(Pi/2);	// rotate 90 deg around Z
	Transform X_bladeF2(R_bladeF2, Vec3(0, g_hubR, 0)); // fixed mob frame at Y edge of hub 
	Rotation R_bladeM2;
	R_bladeM2.setRotationFromAngleAboutY(flapAng);	
	Transform X_bladeM2(R_bladeM2, Vec3(-g_bladeLen / 2, 0, 0));
	MobilizedBody::Weld blade2_mobody(hub_mobody, X_bladeF2,
		blade_body, X_bladeM2);
	g_bi_blade2 = blade2_mobody.getMobilizedBodyIndex();

	// create MobilizedBody for blade 3
	Rotation R_bladeF3;
	R_bladeF3.setRotationFromAngleAboutZ(Pi);	// rotate 180 deg around Z
	Transform X_bladeF3(R_bladeF3, Vec3(-g_hubR, 0, 0)); // fixed mob frame at -X edge of hub
	Rotation R_bladeM3;
	R_bladeM3.setRotationFromAngleAboutY(flapAng);	// flap angle
	Transform X_bladeM3(R_bladeM3, Vec3(-g_bladeLen / 2, 0, 0));
	MobilizedBody::Weld blade3_mobody(hub_mobody, X_bladeF3,
		blade_body, X_bladeM3);
	g_bi_blade3 = blade3_mobody.getMobilizedBodyIndex();

	// create MobilizedBody for blade 4
	Rotation R_bladeF4;
	R_bladeF4.setRotationFromAngleAboutZ(3*Pi/2); // rotate 270 deg  around Z
	Transform X_bladeF4(R_bladeF4, Vec3(0, -g_hubR, 0)); // fixed mob frame at -Y edge of hub
	Rotation R_bladeM4;
	R_bladeM4.setRotationFromAngleAboutY(flapAng);	// flap angle
	Transform X_bladeM4(R_bladeM4, Vec3(-g_bladeLen / 2, 0, 0));
	MobilizedBody::Weld blade4_mobody(hub_mobody, X_bladeF4,
		blade_body, X_bladeM4);
	g_bi_blade4 = blade4_mobody.getMobilizedBodyIndex();
#endif

	// now create body for fins, axes are X = span, Y = normal, Z toward forward edge
	const Vec3 finHalfLengths(finChord / 2, finThick / 2, finLen / 2);
	Inertia finI = finM * UnitInertia::brick(finHalfLengths); // MOI of fin
	Body::Rigid fin_body(MassProperties(finM, Vec3(0), finI));
	printf("Fin: Ixx %.4f, Iyy %.4f, Izz %.4f\r\n",
		finI.getMoments()[0],
		finI.getMoments()[1],
		finI.getMoments()[2]
	);
	// create decoration for the fin
	fin_body.addDecoration(Transform(), DecorativeBrick(finHalfLengths).setColor(Green));

	// now create different body for fin 1
	Body::Rigid fin1_body(MassProperties(finM, Vec3(0), finI));
	printf("Fin1: Ixx %.4f, Iyy %.4f, Izz %.4f\r\n",
		finI.getMoments()[0],
		finI.getMoments()[1],
		finI.getMoments()[2]
	);
	// create decoration for the fin
	fin1_body.addDecoration(Transform(), DecorativeBrick(finHalfLengths).setColor(Red));

	// create MobilizedBody for fin 1
	Rotation R_finF1;
	R_finF1.setRotationFromAngleAboutZ(0);
	Transform X_fin_F1(R_finF1, Vec3(motorR, 0, -motorLen / 2 + finChord / 2));	// fixed mob frame at X edge of motor at rear
	Transform X_fin_M1(Vec3(-finChord / 2, 0, 0));		// moving mob frame
	MobilizedBody::Weld fin1_mobody(motor_mobody, X_fin_F1,		// weld blade to hub
		fin1_body, X_fin_M1);
	g_bi_fin1 = fin1_mobody.getMobilizedBodyIndex();
	// create MobilizedBody for fin 2
	Rotation R_finF2;
	R_finF2.setRotationFromAngleAboutZ(Pi / 2);
	Transform X_fin_F2(R_finF2, Vec3(0, motorR, -motorLen / 2 + finChord / 2));	// fixed mob frame at X edge of motor at rear
	Transform X_fin_M2(Vec3(-finChord / 2, 0, 0));		// moving mob frame
	MobilizedBody::Weld fin2_mobody(motor_mobody, X_fin_F2,		// weld blade to hub
		fin_body, X_fin_M2);
	g_bi_fin2 = fin2_mobody.getMobilizedBodyIndex();

	// create MobilizedBody for fin 3
	Rotation R_finF3;
	R_finF3.setRotationFromAngleAboutZ(Pi);
	Transform X_fin_F3(R_finF3, Vec3(-motorR, 0, -motorLen / 2 + finChord / 2));	// fixed mob frame at X edge of motor at rear
	Transform X_fin_M3(Vec3(-finChord / 2, 0, 0));		// moving mob frame
	MobilizedBody::Weld fin3_mobody(motor_mobody, X_fin_F3,		// weld blade to hub
		fin_body, X_fin_M3);
	g_bi_fin3 = fin3_mobody.getMobilizedBodyIndex();

	// create MobilizedBody for fin 2
	Rotation R_finF4;
	R_finF4.setRotationFromAngleAboutZ(3*Pi / 2);
	Transform X_fin_F4(R_finF4, Vec3(0, -motorR, -motorLen / 2 + finChord / 2));	// fixed mob frame at X edge of motor at rear
	Transform X_fin_M4(Vec3(-finChord / 2, 0, 0));		// moving mob frame
	MobilizedBody::Weld fin4_mobody(motor_mobody, X_fin_F4,		// weld blade to hub
		fin_body, X_fin_M4);
	g_bi_fin4 = fin4_mobody.getMobilizedBodyIndex();


	// add event reporter, call every 0.1 second
	system.addEventReporter(new MyEventReporter(system, hub_mobody, 0.1));

	// Visualize a frame every 0.1 second
    Visualizer viz(system); 
	viz.setDesiredFrameRate(100);
    viz.addDecorationGenerator(new ShowData(system));
    system.addEventReporter(new Visualizer::Reporter(viz, 1./10.));
	viz.addFrameController(new Visualizer::BodyFollower(mast_mobody, Vec3(0), 
	         Vec3(0, 5, 1), UnitVec3(0, 0, 1)));
	viz.zoomCameraToShowAllGeometry();
	viz.setCameraFieldOfView(Pi / 3);

    // Initialize the state, including the custom Force 
	Force::Custom(forces, new rotorForces(matter));
	State state = system.realizeTopology();
	RungeKuttaMersonIntegrator integ(system);
	integ.setAccuracy(1e-4);

	// set initial conditions
	const double initial_ROD = 0;
	const double initial_RPM = 0;
	const double initial_ang_vel = initial_RPM * 2.0 * Pi/60.0;
	const double intial_motor_ang = 0.7*Pi; // radians, relative to Z
//#if 0
//	rotor.setQToFitTranslation(state, Vec3(0, 0, .5));
	// turn rocket on its side
	Rotation rMotor;
	rMotor.setRotationFromAngleAboutY(intial_motor_ang);
	motor_mobody.setQToFitRotation(state, rMotor);
	// set mast angle, Pi/2 = along rocket Z, 0 = 90 degrees
	// only works if Pin instead of Weld mobilizer
	Rotation rMast;
	rMast.setRotationFromAngleAboutZ(0);
	mast_mobody.setQToFitRotation(state, rMast);
	// set initial velocities
	hub_mobody.setUToFitAngularVelocity(state, Vec3(0, 0, initial_ang_vel));
#ifdef USE_ROTOR
	blade1_mobody.setUToFitAngularVelocity(state, Vec3(0, 0, initial_ang_vel));
	blade2_mobody.setUToFitAngularVelocity(state, Vec3(0, 0, initial_ang_vel));
	blade3_mobody.setUToFitAngularVelocity(state, Vec3(0, 0, initial_ang_vel));
	blade4_mobody.setUToFitAngularVelocity(state, Vec3(0, 0, initial_ang_vel));
	motor_mobody.setUToFitLinearVelocity(state, Vec3(0, 0, -initial_ROD));
	hub_mobody.setUToFitLinearVelocity(state, Vec3(0, 0, -initial_ROD));
	blade1_mobody.setUToFitLinearVelocity(state, Vec3(0, 0, -initial_ROD));
	blade2_mobody.setUToFitLinearVelocity(state, Vec3(0, 0, -initial_ROD));
	blade3_mobody.setUToFitLinearVelocity(state, Vec3(0, 0, -initial_ROD));
	blade4_mobody.setUToFitLinearVelocity(state, Vec3(0, 0, -initial_ROD));
#endif
//#endif
	// Simulate it.
    TimeStepper ts(system, integ);
    ts.initialize(state);
	ts.stepTo(100);

	printf("done\r\n");
  } 
  catch (const std::exception& e)
  {
    std::cout << "EXCEPTION: " << e.what() << std::endl;
	  return 1;
  }
  return 0;
}


//==============================================================================
//                      ShowData::generateDecorations()
//==============================================================================
void ShowData::generateDecorations(const State&                state,
                                     Array_<DecorativeGeometry>& geometry)
{
	char s[80];
	static int i = 0;
    m_mbs.realize(state, Stage::Dynamics);
    const Real E=m_mbs.calcEnergy(state);
    DecorativeText altitude;
	Vec3 textCenter = g_hubCenter - Vec3(0, 0, motorLen / 2);	
	altitude.setTransform(Transform(Rotation(), textCenter));
#if 0
    energy.setTransform(Vec3(-.2,0,.5))
            .setText("Energy: " + String(E, "%.6f"))
            .setScale(.09)
            .setColor(Black);
#endif
	sprintf(s, "               ALT %.0f RPM %.0f LIFT %.0f", 
		g_altitude, g_rpm, g_lift);
	altitude.setText(s);
	altitude.setScale(.12);
	altitude.setColor(Black);
	geometry.push_back(altitude);
#if 0
	// draw line from hub showing lift
	double liftLineHalfLen = g_lift / 400;
	DecorativeCylinder liftLine(0.05, liftLineHalfLen);
	Rotation R90Y;
	R90Y.setRotationFromAngleAboutX(Pi / 2);
	liftLine.setTransform(Transform(g_hubRotation*R90Y,g_hubCenter+g_hubUZ* liftLineHalfLen));
	liftLine.setColor(Red);
	geometry.push_back(liftLine);
#endif
	i++;
}
