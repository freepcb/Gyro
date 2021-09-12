#include "Simbody.h"
#include <iostream>

using namespace SimTK;

int main()
{
	MultibodySystem system;
	SimbodyMatterSubsystem matter(system);
	GeneralForceSubsystem forces(system);
	Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
	Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
	pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));
	
	// There are 4 frames of reference for a MobilizedBody: 
	//   P = parent body frame (grey in visualizer)
	//	 F = joint frame in parent body (blue)
	//   M = joint frame in child body (red)
	//   B = child body frame (grey)
	// Transforms: (include both translation and rotation)
	//   X_PF from P to F 
	//   X_BM from B to M 
	// the pin joint fixes the origin of M at the origin of F,
	// with the Z-axes of both frames aligned
	
	// pin joint between Ground and pendulum1: 
	//   F: origin = 1 unit above P origin, axes rotated around X-axis by 0.5 rad
	//   M: 1 unit along Y-axis from B origin
	Transform X_PF(Rotation(0.5, XAxis), Vec3(0, 1, 0));
	Transform X_BM(Vec3(0, 1, 0));
	MobilizedBody::Pin pendulum1(matter.Ground(), X_PF,
			pendulumBody, X_BM);
	MobilizedBody::Pin pendulum2(pendulum1, Transform(Vec3(0)),
		pendulumBody, Transform(Vec3(0, 1, 0)));
	// Set up visualization.
	Visualizer viz(system);
	system.addEventReporter(new Visualizer::Reporter(viz, 0.01));

	// Initialize the system and state.
	system.realizeTopology();
	State state = system.getDefaultState();
	pendulum2.setRate(state, 5.0);
	// Simulate it.
	RungeKuttaMersonIntegrator integ(system);
	TimeStepper ts(system, integ);
	ts.initialize(state);
	ts.stepTo(50.0);
}


#if 0
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
//                                  MAIN
//==============================================================================
int main() {
  try 
  {   
    // Create the system.   
    MultibodySystem system; system.setUpDirection(ZAxis);
    SimbodyMatterSubsystem matter(system);
    // No gravity or other forces

    matter.setShowDefaultGeometry(false); // turn off frames and other junk

    // Construct a single rigid body by welding together a cylindrical shaft
    // and a rectangular bar.
    Rotation YtoX(-Pi/2, ZAxis);
    Body::Rigid shaftBody(MassProperties(1, Vec3(0),
                            UnitInertia::cylinderAlongX(.02, .05)));
    shaftBody.addDecoration(YtoX, 
                            DecorativeCylinder(.02, .05).setColor(Red));

    const Vec3 halfLengths(.02,.04,.3);
    Body::Rigid barBody(MassProperties(2, Vec3(0),
                    UnitInertia::brick(halfLengths)));
    barBody.addDecoration(Transform(),
                          DecorativeBrick(halfLengths).setColor(Blue));

    MobilizedBody::Free shaft(matter.Ground(), Transform(),
                              shaftBody, Transform());
    MobilizedBody::Weld bar(shaft, Vec3(-.05,0,0),
                            barBody, Vec3(halfLengths[0],0,0));

    // Visualize a frame every 1/60 s, and include the energy.
    Visualizer viz(system); 
	viz.setDesiredFrameRate(10);
    viz.addDecorationGenerator(new ShowData(system));
    system.addEventReporter(new Visualizer::Reporter(viz, 1./10));

    // Initialize the system and state. 
    State state = system.realizeTopology();

    // Set initial conditions. Need a slight perturbation of angular velocity
    // to trigger the instability.
    shaft.setQToFitTranslation(state, Vec3(0,0,0.5));
    shaft.setUToFitAngularVelocity(state, Vec3(10,0,1e-10)); // 10 rad/s
    
    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-5);
    TimeStepper ts(system, integ);
    ts.initialize(state);
	// set up loop
	double timeStepSize = 0.1;
	double timeTarget = timeStepSize;
	double timeEnd = 10.0;
	double timeCurrent = 0.0;
	while (timeCurrent < timeEnd)
	{
		while((timeTarget-timeCurrent)>0.001)
		{
			// step until we reach timeTarget, in case stepTo() returns early
			ts.stepTo(timeTarget);
			timeCurrent = ts.getTime();
			printf("%0.3f\r\n", timeCurrent);
		}
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
#endif