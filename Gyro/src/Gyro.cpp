#include "Simbody.h"
#include <iostream>
#include "..\..\..\Gyro\Rotor_test\src\Rotor.h"

using namespace SimTK; 

// global variables for convenience
const Real gravity = -9.8;	// acceleration due to gravity
const Real g_mass = 0.5;	// airframe mass
Real g_last_time = 0.0;		// time of last update
Real g_lift = 0.0;			// rotor lift
Real g_torque = 0.0;		// rotor torque
Real g_rpm = 0.0;			// rotor RPM
Real g_velZ = -0.001;		// initial rotor vert speed to avoid divide-by-zero errors
Real g_altitude = 0.0;		// rotor altitude;
// rotor parameters
Airfoil af;						// use default SG6042 airfoil
const double rotorMOI = 0.2;
const int numBlades = 2;
const int bladeNsegs = 18;		// number of segments to approximate spanwise variation
const double bladeRootR = 0.1;  // distance of blade root from hub axis
const double bladeTipR = 1.0;    // distance of blade tip from hub axis
const double bladeChord = 0.055; // blade chord
const double bladePitch = -6.0;  // blade pitch (degrees)
double bladeLift;				// net upward lift from blade
double bladeTorque;				// net torque from blade
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
		MobilizedBodyIndex bi(1);
		MobilizerUIndex ui(0);
		const MobilizedBody& hub_mobody = matter.getMobilizedBody(bi);
		Vec3 angVel = hub_mobody.getBodyAngularVelocity(state);
		// update aerodynamic forces
		Real t = state.getTime();
		Real lift;		// will be set by bl.getForces()	
		Real torque;	// will be set by bl.getForces()
		const int printLevel = 1;
		bl.getForces(angVel[2], g_velZ, 0, lift, torque);	// don't print
		if(printLevel)
			printf("t: %.3f, lift %.3f, torque %.3f\r\n", t, lift, torque);
		// now apply them for time step
		Real dt = t - g_last_time;
		g_velZ = g_velZ + dt * (gravity + lift / g_mass);
		g_altitude = g_altitude + dt * g_velZ;
		g_torque = torque*2;	// for 2 blades
		g_lift = lift*2;		//     "
		g_last_time = t;
		// apply torque to the hub joint
		mobilityForces[0] = torque;
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
private:
	SimbodyMatterSubsystem& matter;
};

//==============================================================================
//                           EVENT REPORTER
//==============================================================================
class PositionReporter : public PeriodicEventReporter {
public:
	PositionReporter(const MultibodySystem& system, const MobilizedBody& mobod,
		Real reportInterval)
		: PeriodicEventReporter(reportInterval), system(system), mobod(mobod) {}
	void handleEvent(const State& state) const {
		system.realize(state, Stage::Position);   // so that variables can be accessed
//		Vec3 pos = mobod.getBodyOriginLocation(state);
		Vec3 angVel = mobod.getBodyAngularVelocity(state);
		printf("t %0.3f, RPM %.3f, ROD %.2f, alt %.2f, lift %.3f, torque %.3f, \r\n",
			g_last_time, angVel[2] * 60 / (2 * Pi), g_velZ, g_altitude, g_lift,
			g_torque);
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

	// Constuct the rotor from a cylindrical hub with a pin joint
	// and multiple blades attached to the circumference of the hub
	// with welds or hinges
	const Real density = 1000.0;	// density of rotor parts in kg/m^3
	const Real hubR = 0.05;
	const Real hubThick = 0.01;
	const Real bladeLen = 0.95;
	const Real bladeChord = 0.055;
	const Real bladeThick = 0.005;
	const Real bladeMass = density * bladeLen * bladeChord * bladeThick;
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
	Rotation R;
	R.setRotationFromAngleAboutX(Pi/2);
	Transform X(R);
	hub_body.addDecoration(X, DecorativeCylinder(hubR, hubThick).setColor(Blue));
	// create MobilizedBody for the hub
	MobilizedBody::Pin hub_mobody(matter.Ground(), Transform(),
		hub_body, Transform());

	// now create body for the blades
	const Vec3 bladeHalfLengths(bladeLen/2, bladeChord/2, bladeThick/2);
	Inertia bladeI = bladeMass * UnitInertia::brick(bladeHalfLengths); // MOI of blade
	Body::Rigid blade_body(MassProperties(bladeMass, Vec3(0), bladeI));
	printf("Blade: Ixx %.4f, Iyy %.4f, Izz %.4f\r\n",
		bladeI.getMoments()[0],
		bladeI.getMoments()[1],
		bladeI.getMoments()[2]
	);
	// create decoration for the blade
	blade_body.addDecoration(Transform(), DecorativeBrick(bladeHalfLengths).setColor(Blue));

	// create MobilizedBody for blade 1
	Real blade_cm_offset = hubR + bladeLen/2;
	Transform X_PF1(Vec3(blade_cm_offset, 0, 0));
	MobilizedBody::Weld blade_mobody1(hub_mobody, X_PF1,
		blade_body, Transform());
	// create MobilizedBody for blade 2
	R.setRotationFromAngleAboutZ(Pi);
	Transform X_PF2(R,Vec3(-blade_cm_offset, 0, 0));
	MobilizedBody::Weld blade_mobody2(hub_mobody, X_PF2,
		blade_body, Transform());

	// add event reporter, call every 0.1 sec
	system.addEventReporter(new PositionReporter(system, hub_mobody, 0.1));

	// Visualize a frame every 1/10 s
    Visualizer viz(system); 
	viz.setDesiredFrameRate(10);
    viz.addDecorationGenerator(new ShowData(system));
    system.addEventReporter(new Visualizer::Reporter(viz, 1./10));

    // Initialize the state, including the custom Force 
	Force::Custom(forces, new rotorForces(matter));
	State state = system.realizeTopology();
	RungeKuttaMersonIntegrator integ(system);
	integ.setAccuracy(1e-5);

	// set initial conditions
//	rotor.setQToFitTranslation(state, Vec3(0, 0, .5));
//	rotor.setUToFitAngularVelocity(state, Vec3(0, 0, 10)); // 10 rad/s
	MobilizedBodyIndex bi(1);
	MobilizerUIndex ui(0);
	const MobilizedBody& mobody1 = hub_mobody;
	
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
	char s[20];
	static int i = 0;
    m_mbs.realize(state, Stage::Dynamics);
    const Real E=m_mbs.calcEnergy(state);
    DecorativeText energy;
#if 0
    energy.setTransform(Vec3(-.2,0,.5))
            .setText("Energy: " + String(E, "%.6f"))
            .setScale(.09)
            .setColor(Black);
#endif
	sprintf(s, "   %3d", i);
	energy.setText(s);
	energy.setScale(.09);
	energy.setColor(Black);
	geometry.push_back(energy);
	i++;
}
