# Gyro

This is a project to simulate a passive, free-wheeling rotor as a recovery method for a model rocket. In the initial implementation, it would just slow down the descent similar to a parachute, but eventually I would like to add control mechanisms so that the rocket would descend as a gyro-glider under autonomous or radio control.

The rigid body Physics engine is SimBody. Aerodynamic forces on the rotor are calculated using standard airfoil theory. The physics engine is needed to handle gyroscopic and Coriolis forces, which are signigicant for a rotorcraft.
