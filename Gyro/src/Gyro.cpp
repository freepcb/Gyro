#include "Simbody.h"
#include <iostream>
#include "Rotor.h"

using namespace SimTK; 

// global variables for convenience
// indices for mobilized bodies so they can be called by name
int i_base = 1;
int i_motor = i_base++;
int i_hub = i_base++;
int i_bl1 = i_base++;
int i_bl2 = i_base++;
int i_bl3 = i_base++;
int i_bl4 = i_base++;

const Real gravity = -9.8;	// acceleration due to gravity
const Real g_mass = 0.5;	// airframe mass
Real g_last_time = 0.0;		// time of last update
Real g_lift = 0.0;			// rotor lift
Real g_torque = 0.0;		// rotor torque
Real g_angVel = 0.0;		// rotor ang vel
Real g_rpm = 0.0;			// rotor RPM
// Real g_velZ = -1.0;		// initial rotor vert speed to avoid divide-by-zero errors
Real g_altitude = 0.0;		// rotor altitude;
Vec3 g_hubCenter(0);

// rotor objects and parameters
Airfoil af;						 // airfoil for blades, uses SG6042 by defailt
const int numBlades = 2;		 // num blades in rotor
const int bladeNsegs = 18;		 // number of segments to approximate spanwise variation
const double bladeRootR = 0.1;   // distance of blade root from hub axis
const double bladeTipR = 1.0;    // distance of blade tip from hub axis
const double bladeChord = 0.055; // blade airfoil chord
const double bladePitch = -6.0;  // blade pitch (degrees)
double bladeLift;				 // net upward lift from blade
double bladeTorque;				 // net torque from blade
// single blade object, since all blades are the same
RotorBlade bl(&af, bladeRootR, bladeTipR, bladeChord, bladePitch, bladeNsegs);


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
		// get mobilized bodies for hub and blades
		MobilizedBodyIndex bi_hub(i_hub);	
		const MobilizedBody& hub_mobody = matter.getMobilizedBody(bi_hub);
		MobilizedBodyIndex bi_bl1(i_bl1), bi_bl2(i_bl2), bi_bl3(i_bl3), bi_bl4(i_bl4);
		const MobilizedBody& blade1_mobody = matter.getMobilizedBody(bi_bl1);
		const MobilizedBody& blade2_mobody = matter.getMobilizedBody(bi_bl2);
		const MobilizedBody& blade3_mobody = matter.getMobilizedBody(bi_bl3);
		const MobilizedBody& blade4_mobody = matter.getMobilizedBody(bi_bl4);
		// print hub motion, same format as in event reporter
		SpatialVec vel = hub_mobody.getBodyVelocity(state);
//		printf("force: t %0.5f, hub: v (%.2f,%.2f,%.2f) RPM %.5f, alt %.2f\r\n",
//			state.getTime(), vel[1][0], vel[1][1], vel[1][2], vel[0][2] * 60 / (2 * Pi), g_altitude);
		g_angVel = vel[0][2];

		Vec3 angVel = vel[0];		// get angles of both blades (should always be 360/num_blades apart)
		Vec3 blade1XInG = blade1_mobody.expressVectorInGroundFrame(state, Vec3(1, 0, 0));
		double blade1Angle = (180 / Pi)*atan2(blade1XInG[1], blade1XInG[0]);
		Real t = state.getTime();
		// update aerodynamic forces
		// these values will actually be set by bl.getForces()
		Real lift1 = 0;	Real rLift1 = 0; Real thrust1 = 0; Real rThrust1 = 0;
		Real lift2 = 0;	Real rLift2 = 0; Real thrust2 = 0; Real rThrust2 = 0;
		Real lift3 = 0; Real rLift3 = 0; Real thrust3 = 0; Real rThrust3 = 0;
		Real lift4 = 0;	Real rLift4 = 0; Real thrust4 = 0; Real rThrust4 = 0;
		const int printLevel = 0;
		bl.getForces(state, blade1_mobody, Vec3(0, 0, 0), printLevel, lift1, rLift1, thrust1, rThrust1);
		bl.getForces(state, blade2_mobody, Vec3(0, 0, 0), printLevel, lift2, rLift2, thrust2, rThrust2);
		bl.getForces(state, blade3_mobody, Vec3(0, 0, 0), printLevel, lift3, rLift3, thrust3, rThrust3);
		bl.getForces(state, blade4_mobody, Vec3(0, 0, 0), printLevel, lift4, rLift4, thrust4, rThrust4);
//		g_torque = torque1 + torque2 + torque3 + torque4;
//		g_lift = lift1 + lift2 + lift3 + lift4;
		Vec3 hubVel = vel[1];
		Vec3 hubPos = hub_mobody.getBodyOriginLocation(state);
		g_altitude = hubPos[2];
		if (printLevel)
			printf("t: %.5f, angle %5.1f, RPM %6.2f, alt %5.0f, ROD %5.2f, L: %.5f %.5f %.5f %.5f T: %.5f %.5f %.5f %.5f\r\n", 
				t, blade1Angle, angVel[2]*60/(2*Pi), hubPos[2], hubVel[2], 
				lift1, lift2, lift3, lift4, thrust1, thrust2, thrust3, thrust4);
		g_last_time = t;
		// apply lift and thrust to blades as external forces
		// start by creating vectors in blade frames
		Vec3 vecThrust1(0, thrust1, 0);	
		Vec3 vecThrust2(0, thrust2, 0);
		Vec3 vecThrust3(0, thrust3, 0);
		Vec3 vecThrust4(0, thrust4, 0);
		Vec3 vecLift1(0, 0, lift1);
		Vec3 vecLift2(0, 0, lift2);
		Vec3 vecLift3(0, 0, lift3);
		Vec3 vecLift4(0, 0, lift4);
		Vec3 LiftInG1 = blade1_mobody.expressVectorInGroundFrame(state, vecLift1);
		Vec3 ThrustInG1 = blade1_mobody.expressVectorInGroundFrame(state, vecThrust1);
		blade1_mobody.applyForceToBodyPoint(state, Vec3(rLift1, 0, 0), LiftInG1, bodyForces);
		blade1_mobody.applyForceToBodyPoint(state, Vec3(rThrust1, 0, 0), ThrustInG1, bodyForces);
		Vec3 LiftInG2 = blade2_mobody.expressVectorInGroundFrame(state, vecLift2);
		Vec3 ThrustInG2 = blade2_mobody.expressVectorInGroundFrame(state, vecThrust2);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(rLift2, 0, 0), LiftInG2, bodyForces);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(rThrust2, 0, 0), ThrustInG2, bodyForces);
		Vec3 LiftInG3 = blade3_mobody.expressVectorInGroundFrame(state, vecLift3);
		Vec3 ThrustInG3 = blade3_mobody.expressVectorInGroundFrame(state, vecThrust3);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(rLift3, 0, 0), LiftInG3, bodyForces);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(rThrust3, 0, 0), ThrustInG3, bodyForces);
		Vec3 LiftInG4 = blade4_mobody.expressVectorInGroundFrame(state, vecLift4);
		Vec3 ThrustInG4 = blade4_mobody.expressVectorInGroundFrame(state, vecThrust4);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(rLift4, 0, 0), LiftInG4, bodyForces);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(rThrust4, 0, 0), ThrustInG4, bodyForces);
		g_hubCenter = hub_mobody.findStationLocationInGround(state, Vec3(0));
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
		printf("event:t %0.5f, hub:z (%5.2f,%5.2f,%5.2f) v (%5.2f,%5.2f,%5.2f) RPM %6.2f, alt %.2f\r\n",
			state.getTime(), UZInG[0], UZInG[1], UZInG[2], 
			vel[1][0], vel[1][1], vel[1][2], velMobilizer[0][2] * 60 / (2 * Pi), g_altitude);
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
//	matter.setShowDefaultGeometry(false); // turn off frames 

	// keep track of MobilizedBody indices
	int bi = 0;

	// Constuct the rotor from a cylindrical hub with a pin joint
	// and multiple blades attached to the circumference of the hub
	// with welds or hinges
	const Real bladeDensity = 1000.0;	// density of rotor blades in kg/m^3
	const Real hubR = 0.05;
	const Real hubThick = 0.01;
	const Real bladeLen = 0.9;
	const Real bladeChord = 0.055;
	const Real bladeThick = 0.005;
	const Real bladeM = 0.25;

	// create motor body
	const Real motorLen = 2;	
	const Real motorR = 0.05;	// 5 cm
	const Real motorM = 0.5;	
	Inertia motorI = motorM * UnitInertia::cylinderAlongZ(motorR, motorLen/2); // MOI of motor 
	Body::Rigid motor_body(MassProperties(motorM, Vec3(0), motorI));
	// create decoration for the motor
	Rotation R;
	R.setRotationFromAngleAboutX(Pi / 2);
	Transform X(R);
	// note: DecorativeCylinder is oriented along Y so needs to be rotated
	// use half/length
	motor_body.addDecoration(X, DecorativeCylinder(motorR, motorLen/2).setColor(Red));
	// create MobilizedBody for the motor
	MobilizedBody::Free motor_mobody(matter.Ground(), Transform(),
		motor_body, Transform(Vec3(0,0, -motorLen/2)));	// joint at bottom of motor
	MobilizedBodyIndex bi_motor(i_motor);


	// create hub body
	Real hubM = 0;	// ignore mass properties of hub
	Inertia hubI = hubM * UnitInertia::cylinderAlongZ(hubR, hubThick); // MOI of hub
	Body::Rigid hub_body(MassProperties(hubM, Vec3(0), hubI));
	printf("Hub: Ixx %.4f, Iyy %.4f, Izz %.4f\r\n",
		hubI.getMoments()[0],
		hubI.getMoments()[1],
		hubI.getMoments()[2]
		);
	// create decoration for the hub
	Rotation hub_R;
	hub_R.setRotationFromAngleAboutX(Pi/2);
	Transform hubX(hub_R);
	hub_body.addDecoration(hubX, DecorativeSphere(hubR).setColor(Red));
	// create MobilizedBody for the hub
	// joint at top of motor and center of hub
	MobilizedBody::Pin hub_mobody(motor_mobody, Transform(Vec3(0, 0, motorLen/2)), 
		hub_body, Transform());
	MobilizedBodyIndex bi_hub(i_hub);

	// now create body for the blades
	const Vec3 bladeHalfLengths(bladeLen/2, bladeChord/2, bladeThick/2);
	Inertia bladeI = bladeM * UnitInertia::brick(bladeHalfLengths); // MOI of blade
	Body::Rigid blade_body(MassProperties(bladeM, Vec3(0), bladeI));
	printf("Blade: Ixx %.4f, Iyy %.4f, Izz %.4f\r\n",
		bladeI.getMoments()[0],
		bladeI.getMoments()[1],
		bladeI.getMoments()[2]
	);
	// create decoration for the blade
	blade_body.addDecoration(Transform(), DecorativeBrick(bladeHalfLengths).setColor(Blue));

	// create MobilizedBody for blade 1
	const double flapAng = Pi/8; 
	Transform X_F1(Vec3(hubR, 0, 0));
	Rotation blade_RF1;
	blade_RF1.setRotationFromAngleAboutY(flapAng);	// flap angle
	Transform X_M1(blade_RF1, Vec3(-bladeLen / 2, 0, 0));
	MobilizedBody::Weld blade1_mobody(hub_mobody, X_F1,
		blade_body, X_M1);
	MobilizedBodyIndex bi_blade1(i_bl1);

	// create MobilizedBody for blade 2
	Rotation blade_RF2;
	blade_RF2.setRotationFromAngleAboutZ(Pi/2);	// rotate 90 deg around Z
	Transform X_F2(blade_RF2, Vec3(0, hubR, 0)); // 
	Rotation blade_RM2;
	blade_RM2.setRotationFromAngleAboutY(flapAng);	// flap angle
	Transform X_M2(blade_RM2, Vec3(-bladeLen / 2, 0, 0));
	MobilizedBody::Weld blade2_mobody(hub_mobody, X_F2,
		blade_body, X_M2);
	MobilizedBodyIndex bi_blade2(i_bl2);

	// create MobilizedBody for blade 3
	Rotation blade_RF3;
	blade_RF3.setRotationFromAngleAboutZ(Pi);	// rotate 180 deg around Z, to reverse Y
	Transform X_F3(blade_RF3, Vec3(-hubR, 0, 0)); // shift opposite to blade 1
	Rotation blade_RM3;
	blade_RM3.setRotationFromAngleAboutY(flapAng);	// flap angle
	Transform X_M3(blade_RM3, Vec3(-bladeLen / 2, 0, 0));
	MobilizedBody::Weld blade3_mobody(hub_mobody, X_F3,
		blade_body, X_M3);
	MobilizedBodyIndex bi_blade3(i_bl3);

	// create MobilizedBody for blade 4
	Rotation blade_RF4;
	blade_RF4.setRotationFromAngleAboutZ(3*Pi/2);	// rotate 270 deg  around Z, to reverse Y
	Transform X_F4(blade_RF4, Vec3(0, -hubR, 0)); // shift opposite to blade 3
	Rotation blade_RM4;
	blade_RM4.setRotationFromAngleAboutY(flapAng);	// flap angle
	Transform X_M4(blade_RM4, Vec3(-bladeLen / 2, 0, 0));
	MobilizedBody::Weld blade4_mobody(hub_mobody, X_F4,
		blade_body, X_M4);
	MobilizedBodyIndex bi_blade4(i_bl4);


	// add event reporter, call every 0.01 sec
	system.addEventReporter(new MyEventReporter(system, hub_mobody, 0.1));

	// Visualize a frame every 1/10 s
    Visualizer viz(system); 
	viz.setDesiredFrameRate(100);
    viz.addDecorationGenerator(new ShowData(system));
    system.addEventReporter(new Visualizer::Reporter(viz, 1./10.));
	viz.addFrameController(new Visualizer::BodyFollower(hub_mobody, Vec3(0), 
	         Vec3(0, 5, 1), UnitVec3(0, 0, 1)));


    // Initialize the state, including the custom Force 
	Force::Custom(forces, new rotorForces(matter));
	State state = system.realizeTopology();
	RungeKuttaMersonIntegrator integ(system);
	integ.setAccuracy(1e-4);

	// set initial conditions
//	rotor.setQToFitTranslation(state, Vec3(0, 0, .5));
//	hub_mobody.setUToFitAngularVelocity(state, Vec3(0, 0, 5));
//	blade1_mobody.setUToFitAngularVelocity(state, Vec3(0, 0, 5));
//	blade2_mobody.setUToFitAngularVelocity(state, Vec3(0, 0, 5));
	motor_mobody.setUToFitLinearVelocity(state, Vec3(0, .1, -1));
//	hub_mobody.setUToFitLinearVelocity(state, Vec3(0, 0, -1));
//	blade1_mobody.setUToFitLinearVelocity(state, Vec3(0, 0, -1));
//	blade2_mobody.setUToFitLinearVelocity(state, Vec3(0, 0, -1));

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
	char s[40];
	static int i = 0;
    m_mbs.realize(state, Stage::Dynamics);
    const Real E=m_mbs.calcEnergy(state);
    DecorativeText altitude;
	altitude.setTransform(Transform(Rotation(), g_hubCenter));
#if 0
    energy.setTransform(Vec3(-.2,0,.5))
            .setText("Energy: " + String(E, "%.6f"))
            .setScale(.09)
            .setColor(Black);
#endif
	sprintf(s, "               ALT %4.0f", g_altitude);
	altitude.setText(s);
	altitude.setScale(.09);
	altitude.setColor(Black);
	geometry.push_back(altitude);
	i++;
}
