#include "Simbody.h"
#include <iostream>

using namespace SimTK;

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
class rotorTorque : public Force::Custom::Implementation 
{
public:
	rotorTorque(SimbodyMatterSubsystem& matter) : matter(matter) {}

	void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
			Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
	{
		Vec3 angVel(0);
		MobilizedBodyIndex bi(1);
		MobilizerUIndex ui(0);
		const MobilizedBody& hub_mobody = matter.getMobilizedBody(bi);
		mobilityForces[0] = 0.01;
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
		std::cout << state.getTime() << "\t" << angVel[2] << std::endl;
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
	// for now, assign moments of inertia for entire rotor to hub
	Real density = 1000.0;	// density of rotor in kg/m^3
	const Vec3 halfLengths(1.0, 0.055/2, 0.005/2);	// dimensions of hub + 2 blades
	Real rotor_mass = density * 2.* 0.055 * 0.0055; // mass of rotor
	Inertia hubI = rotor_mass * UnitInertia::brick(halfLengths); // MOI of rotor
	Body::Rigid hub(MassProperties(0.5, Vec3(0),
		UnitInertia::cylinderAlongZ(0.05, 0.01))); 
//	UnitInertia::brick(halfLengths)));  //TODO
	printf("Hub: Ixx %.4f, Iyy %.4f, Izz %.4f\r\n",
		hubI.getMoments()[0],
		hubI.getMoments()[1],
		hubI.getMoments()[2]
		);
	// create decoration for the hub
	Rotation R;
	R.setRotationFromAngleAboutX(Pi/2);
	Transform X(R);
	hub.addDecoration(X, DecorativeCylinder(0.05, 0.01).setColor(Blue));
	// create MobilizedBody for the hub
	MobilizedBody::Pin hub_mobody(matter.Ground(), Transform(),
		hub, Transform());

	// now create mobilized bodies for the blades
	const Vec3 bladeHalfLengths(0.95/2, 0.055/2, 0.005/2);	// dimensions of blade
	Real bladeMass = density * 0.95*0.055*0.005;
	Inertia bladeI = bladeMass * UnitInertia::brick(bladeHalfLengths); // MOI of blade
	Body::Rigid blade(MassProperties(bladeMass, Vec3(0),
		UnitInertia::brick(bladeHalfLengths)));

	// add event reporter, call every 0.1 sec
	system.addEventReporter(new PositionReporter(system, hub_mobody, 0.1));

	// Visualize a frame every 1/60 s, and include the energy.
    Visualizer viz(system); 
	viz.setDesiredFrameRate(10);
    viz.addDecorationGenerator(new ShowData(system));
    system.addEventReporter(new Visualizer::Reporter(viz, 1./10));

    // Initialize the state, including the custom Force 
	Force::Custom(forces, new rotorTorque(matter));
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
	ts.stepTo(10);

#if 0
	// set up loop
	double timeStepSize = 0.01;
	double timeTarget = timeStepSize;
	double timeEnd = 10.0;
	double timeCurrent = 0.0;
	double torque = 0.0;

	Vec3 angVel(0);
	SpatialVec Vel6(0);
	while (timeCurrent < timeEnd)
	{
		while ((timeTarget - timeCurrent) > 0.0001)
		{
			// step until we reach timeTarget, in case stepTo() returns early
			ts.stepTo(timeTarget);
			timeCurrent = ts.getTime();
		}
		// update rotor torque
		timeTarget += timeStepSize;
	}
#endif
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
