# Gyro

This is a C++ project to simulate a passive, free-wheeling rotor as a recovery device for a model rocket. In the initial implementation, it would just slow down the descent similar to a parachute, but eventually I would like to add control mechanisms so that the rocket would descend as a gyro-glider under autonomous or radio control.

Aerodynamic forces on the rotor are calculated using standard airfoil theory, currently for a NACA0015 airfoil. A Rigid Body Physics engine is used for the dynamics, mainly to handle gyroscopic and Coriolis forces which are significant for a rotorcraft, particularly at model scale. For this I am using Simbody, from Stanford and the NIH, which seems to provide better accuracy than the gaming engines, and the documentation is excellent.

I am using Microsoft Visual Studio 2017 on Windows 10
