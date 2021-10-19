#include "Simbody.h"
#include <iostream>
#include "Rotor.h"

using namespace SimTK; 

// global variables for convenience and laziness
// indices for mobilized bodies so they can be called by name
int i_base = 1;				// ground frame
int i_motor = i_base++;		// motor mount
int i_hub = i_base++;		// rotor hub
int i_bl1 = i_base++;		// rotor blades
int i_bl2 = i_base++;
int i_bl3 = i_base++;
int i_bl4 = i_base++;

// dimensions of rotor hub and blades
const Real bladeDensity = 1000.0;	// density of rotor blades in kg/m^3
const Real hubR = 0.05;				// hub radius, where blades are attached
const Real hubThick = 0.01;			// hub thickness
const Real bladeLen = 0.95;			// length of blade with airfoil
const double bladeRootR = hubR;		// distance of blade root from hub axis
const double bladeTipR = bladeLen + hubR;  // distance of blade tip from hub axis
const Real bladeChord = 0.055;		// chord length of airfoil
const Real bladeThick = 0.005;		// average thickness of blade with airfoil
const Real bladeM = bladeDensity* bladeLen*bladeChord*bladeThick; // blade mass
const int numBlades = 4;			// number of blades in rotor
const int bladeNsegs = 18;			// number of airfoil segments per blades 

// physical and simulation constants
const Real gravity = -9.8;	// acceleration due to gravity in kg/sec^2
const Real g_mass = 1.5;	// airframe mass
Real g_last_time = 0.0;		// time of last update
Real g_lift = 0.0;			// rotor lift
Real g_torque = 0.0;		// rotor torque
Real g_angVel = 0.0;		// rotor ang vel
Real g_rpm = 0.0;			// rotor RPM
Real g_altitude = 0.0;		// rotor altitude;
Vec3 g_hubCenter(0);

// rotor objects and parameters
Airfoil af;						 // airfoil for blades
// distance of blade tip from hub axis
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

		Vec3 angVel = vel[0];		// get angleof blade 1
		Vec3 blade1XInG = blade1_mobody.expressVectorInGroundFrame(state, Vec3(1, 0, 0));
		double blade1Angle = (180 / Pi)*atan2(blade1XInG[1], blade1XInG[0]);
		Real t = state.getTime();
		// update aerodynamic forces
		// these values will actually be set by bl.getForces()
		Real FZroot1 = 0; Real FZtip1 = 0; Real FYroot1 = 0; Real FYtip1 = 0;
		Real FZroot2 = 0; Real FZtip2 = 0; Real FYroot2 = 0; Real FYtip2 = 0;
		Real FZroot3 = 0; Real FZtip3 = 0; Real FYroot3 = 0; Real FYtip3 = 0;
		Real FZroot4 = 0; Real FZtip4 = 0; Real FYroot4 = 0; Real FYtip4 = 0;
		double xoff = bl.m_bladeLenX/2;
		const int printLevel = 2;
		bl.getForces(state, blade1_mobody, Vec3(0, 0, 0), printLevel, FZroot1, FZtip1, FYroot1, FYtip1);
		bl.getForces(state, blade2_mobody, Vec3(0, 0, 0), printLevel, FZroot2, FZtip2, FYroot2, FYtip2);
		bl.getForces(state, blade3_mobody, Vec3(0, 0, 0), printLevel, FZroot3, FZtip3, FYroot3, FYtip3);
		bl.getForces(state, blade4_mobody, Vec3(0, 0, 0), printLevel, FZroot4, FZtip4, FYroot4, FYtip4);
		double torque1 = FYtip1*bladeLen;
		double torque2 = FYtip2*bladeLen;
		double torque3 = FYtip3*bladeLen;
		double torque4 = FYtip4*bladeLen;
		double lift1 = FZroot1 + FZtip1;
		double lift2 = FZroot2 + FZtip2;
		double lift3 = FZroot3 + FZtip3;
		double lift4 = FZroot4 + FZtip4;
		g_torque = torque1 + torque2 + torque3 + torque4;
		g_lift = lift1 + lift2 + lift3 + lift4;
		Vec3 hubVel = vel[1];
		Vec3 hubPos = hub_mobody.getBodyOriginLocation(state);
		g_altitude = hubPos[2];
		if (printLevel)
			printf("forces: t %.4f, angle %5.1f, RPM %6.2f, alt %3.0f, ROD %5.2f, L: %.5f (%.5f+%.5f+%.5f+%.5f) T: %.5f (%.5f+%.5f+%.5f+%.5f)\r\n", 
				t, blade1Angle, angVel[2]*60/(2*Pi), hubPos[2], hubVel[2], 
				g_lift, lift1, lift2, lift3, lift4, g_torque, torque1, torque2, torque3, torque4);
		g_last_time = t;
		// apply lift and thrust to blades as external forces
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
		blade1_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG1, bodyForces);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG2, bodyForces);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG3, bodyForces);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG4, bodyForces);
		blade1_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG1, bodyForces);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG2, bodyForces);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG3, bodyForces);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG4, bodyForces);
		blade1_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG1, bodyForces);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG2, bodyForces);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG3, bodyForces);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(-xoff, 0, 0), vecFYrootInG4, bodyForces);
		blade1_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG1, bodyForces);
		blade2_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG2, bodyForces);
		blade3_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG3, bodyForces);
		blade4_mobody.applyForceToBodyPoint(state, Vec3(xoff, 0, 0), vecFYtipInG4, bodyForces);

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
		printf("event:t %0.3f, hub:z (%5.2f,%5.2f,%5.2f) v (%5.2f,%5.2f,%5.2f) RPM %6.2f, alt %.2f, lift %5.2f, torque %5.3f\r\n",
			state.getTime(), UZInG[0], UZInG[1], UZInG[2], 
			vel[1][0], vel[1][1], vel[1][2], velMobilizer[0][2] * 60 / (2 * Pi), g_altitude,
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
//	matter.setShowDefaultGeometry(false); // turn off frames 

	// keep track of MobilizedBody indices
	int bi = 0;

	// Constuct the rotor from a cylindrical hub with a pin joint
	// and multiple blades attached to the circumference of the hub
	// with welds or hinges
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
	const double flapAng = 0; 
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
	const double initial_ROD = 7.38;
	const double initial_RPM = 11.5;
	const double initial_ang_vel = initial_RPM * 2.0 * Pi/60.0;
//#if 0
//	rotor.setQToFitTranslation(state, Vec3(0, 0, .5));
	hub_mobody.setUToFitAngularVelocity(state, Vec3(0, 0, initial_ang_vel));
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
