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
		const MobilizedBody& mobody1 = matter.getMobilizedBody(bi);
		angVel = mobody1.getBodyAngularVelocity(state);	
		printf("rpm = %.3f\r\n", 60*angVel[2]/(2*Pi));
		mobilityForces[0] = 0.05;
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
//                                  MAIN
//==============================================================================
int main() {
  try 
  {   
    // Create the system.   
    MultibodySystem system; system.setUpDirection(ZAxis);
    SimbodyMatterSubsystem matter(system);
	GeneralForceSubsystem forces(system);

	// Constuct the rotor as a single rectangular shape
	Real density = 1000.0;	// density of rotor in kg/m^3
	const Vec3 halfLengths(1.0, 0.055/2, 0.005/2);
	Real rotor_mass = density * 2.* 0.055 * 0.0055;
	Inertia rotor_MOI = rotor_mass * UnitInertia::brick(halfLengths);
	Body::Rigid rotorBody(MassProperties(0.5, Vec3(0),
							UnitInertia::brick(halfLengths)));
	printf("Moments of inertia: %.4f, %.4f, %.4f\r\n", rotor_MOI.getMoments()[0],
		rotor_MOI.getMoments()[1],
		rotor_MOI.getMoments()[2]);
	rotorBody.addDecoration(Transform(),
		DecorativeBrick(halfLengths).setColor(Blue));

	MobilizedBody::Pin rotor(matter.Ground(), Transform(),
		rotorBody, Transform());
	printf("Num mobodies %d\r\n", matter.getNumBodies());

	// Visualize a frame every 1/60 s, and include the energy.
    Visualizer viz(system); 
	viz.setDesiredFrameRate(10);
    viz.addDecorationGenerator(new ShowData(system));
    system.addEventReporter(new Visualizer::Reporter(viz, 1./10));

    // Initialize the system and state. 
	Force::Custom(forces, new rotorTorque(matter));
//	Force::DiscreteForces rotorDiscreteForce(forces, matter);
	State state = system.realizeTopology();

	// set initial conditions
	rotor.setQToFitTranslation(state, Vec3(0, 0, .5));
//	rotor.setUToFitAngularVelocity(state, Vec3(0, 0, 10)); // 10 rad/s
	MobilizedBodyIndex bi(1);
	MobilizerUIndex ui(0);
	const MobilizedBody& mobody1 = matter.getMobilizedBody(bi);
//	rotorDiscreteForce.setOneMobilityForce(state, mobody1, ui, 0.1);

    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-5);
    TimeStepper ts(system, integ);
    ts.initialize(state);

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
			const MobilizedBody& mobody1 = matter.getMobilizedBody(bi);
//			torque = rotorDiscreteForce.getOneMobilityForce(state, mobody1, ui);
			//			angVel = body1.getBodyAngularVelocity(state);
			//			Vel6 = mobody1.getMobilizerVelocity(state);
//			angVel = rotor.getBodyAngularVelocity(state);
			printf("%.3f %.3f ", timeCurrent, torque);
			printf("\r\n");
		}
		// update rotor torque
//		rotorDiscreteForce.clearAllMobilityForces(state);
		timeTarget += timeStepSize;
	}
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
//                              SHOW DATA
//==============================================================================

void ShowData::generateDecorations(const State&                state,
                                     Array_<DecorativeGeometry>& geometry)
{
	char s[10];
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
	sprintf(s, "%3d", i);
	energy.setText(s);
	energy.setScale(.09);
	energy.setColor(Black);
	geometry.push_back(energy);
	i++;
}
