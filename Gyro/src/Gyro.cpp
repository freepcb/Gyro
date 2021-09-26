#include "Simbody.h"
#include <iostream>
#include "Rotor.h"

using namespace SimTK; 

// global variables for convenience
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
		MobilizedBodyIndex bi_hub(1);	// assumes hub is first body connected to Ground
		const MobilizedBody& hub_mobody = matter.getMobilizedBody(bi_hub);
		MobilizedBodyIndex bi_bl1(2), bi_bl2(3);	// assume blade body indices follow hub
		const MobilizedBody& blade1_mobody = matter.getMobilizedBody(bi_bl1);
		const MobilizedBody& blade2_mobody = matter.getMobilizedBody(bi_bl2);
		// get angular velocity of hub for reporting purposes
		Vec3 angVel = hub_mobody.getBodyAngularVelocity(state);
		g_angVel = angVel[2];
		// get angles of both blades (should always be 360/num_blades apart)
		Vec3 blade1XInG = blade1_mobody.expressVectorInGroundFrame(state, Vec3(1, 0, 0));
		Vec3 blade2XInG = blade2_mobody.expressVectorInGroundFrame(state, Vec3(1, 0, 0));
		double blade1Angle = (180 / Pi)*atan2(blade1XInG[1], blade1XInG[0]);
		double blade2Angle = (180 / Pi)*atan2(blade2XInG[1], blade2XInG[0]);
		// update aerodynamic forces
		Real t = state.getTime();
		Real lift1 = 0;		// values will actually be set by bl.getForces()	
		Real torque1 = 0;	//            "
		Real lift2 = 0;		//            "
		Real torque2 = 0;	//            "
		const int printLevel = 0;
		// added slight crosswind along X
//		bl.getForces(state, blade1_mobody, Vec3(.01, 0, -g_velZ), printLevel, lift1, torque1);
//		bl.getForces(state, blade2_mobody, Vec3(.01, 0, -g_velZ), printLevel, lift2, torque2);
		bl.getForces(state, blade1_mobody, Vec3(0, 0, 0), printLevel, lift1, torque1);
		bl.getForces(state, blade2_mobody, Vec3(0, 0, 0), printLevel, lift2, torque2);
		g_torque = torque1 + torque2;
		g_lift = lift1 + lift2;
		Vec3 hubVel = hub_mobody.getBodyOriginVelocity(state);
		Vec3 hubPos = hub_mobody.getBodyOriginLocation(state);
		g_altitude = hubPos[2];
//		if (printLevel)
			printf("t: %.5f, angle %5.1f, RPM %6.2f, alt %5.0f, ROD %5.2f, L1 %.5f, T1 %.5f, L2 %.5f, T2 %.5f\r\n", 
				t, blade1Angle, angVel[2]*60/(2*Pi), hubPos[2], hubVel[2], lift1, torque1, lift2, torque2);
		g_last_time = t;
		// apply torque to blades as external forces
		double thrust1 = torque1 / ((bladeRootR + bladeTipR) / 2);	// force to be applied at CM 
		double thrust2 = torque2 / ((bladeRootR + bladeTipR) / 2);	// force to be applied at CM 
//		double diff_lift = (lift2 - lift1)/2;
//		Vec3 vecForce1(0, force1, 0);	// force along local Y  to be transformed to Ground frame
//		Vec3 vecForce2(0, force2, 0);	// for now, same for both blades
		Vec3 vecForce1(0, thrust1, lift1);	// force along local Y  to be transformed to Ground frame
		Vec3 vecForce2(0, thrust2, lift2);
		Vec3 forceInG1 = blade1_mobody.expressVectorInGroundFrame(state, vecForce1);
		blade1_mobody.applyBodyForce(state, SpatialVec(Vec3(0), forceInG1), bodyForces);
		Vec3 forceInG2 = blade2_mobody.expressVectorInGroundFrame(state, vecForce2);
		blade2_mobody.applyBodyForce(state, SpatialVec(Vec3(0), forceInG2), bodyForces);
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
		Vec3 angVel = mobod.getBodyAngularVelocity(state);
		printf("event: t %0.5f, RPM %.5f, alt %.2f, lift %.5f, torque %.5f, \r\n",
			state.getTime(), angVel[2] * 60 / (2 * Pi), g_altitude, g_lift, g_torque);
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
	const Real hubR = 0.1;
	const Real hubThick = 0.01;
	const Real bladeLen = 0.9;
	const Real bladeChord = 0.055;
	const Real bladeThick = 0.005;
	const Real bladeM = 0.2;

	// create motor body
	const Real motorLen = 1;	// 20 cm
	const Real motorR = 0.05;	// 5 cm
	const Real motorM = 1;		// 5 kg
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
	bi++;
	MobilizedBodyIndex bi_motor(bi);


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
	MobilizedBody::Weld hub_mobody(motor_mobody, Transform(Vec3(0, 0, motorLen/2)), 
		hub_body, Transform());
	bi++;
	MobilizedBodyIndex bi_hub(bi);

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
	Real blade_cm_offset = hubR + bladeLen/2;  
	Transform X_F1(Vec3(hubR, 0, 0));
	Transform X_M1(Vec3(-bladeLen / 2, 0, 0));
	MobilizedBody::Weld blade_mobody1(hub_mobody, X_F1,
		blade_body, X_M1);
	bi++;
	MobilizedBodyIndex bi_blade1(bi);

	// create MobilizedBody for blade 2
	Rotation blade_R;
	blade_R.setRotationFromAngleAboutZ(Pi);	// rotate 180 around Z, to reverse Y
	Transform X_F2(blade_R,Vec3(-hubR, 0, 0)); // shift opposite to blade 1
	Transform X_M2(Vec3(-bladeLen / 2, 0, 0));
	MobilizedBody::Weld blade_mobody2(hub_mobody, X_F2,
		blade_body, X_M2);
	bi++;
	MobilizedBodyIndex bi_blade2(bi);

	// add event reporter, call every 0.01 sec
	system.addEventReporter(new MyEventReporter(system, hub_mobody, 0.01));

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
	integ.setAccuracy(1e-5);

	// set initial conditions
//	rotor.setQToFitTranslation(state, Vec3(0, 0, .5));
//	hub_mobody.setUToFitLinearVelocity(state, Vec3(0.3, 0, 0)); 
	
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
