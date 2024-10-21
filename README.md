# dymp (DYnamics-based Motion Planner)

## Abstract

dymp is a C++ library for trajectory optimization (TO) of legged robots.
While dymp provides core functionality for TO, another library [dymp_mpc](https://github.com/ytazz/dymp_mpc)
provides implementation of model predictive controllers that run on Choreonoid.

## Building

Use CMake.

### Dependent libraries

- Eigen3
- Intel MKL
- glfw3 (required for running samples)

## Running

There is a sample code of trajectory generation based on the centroidal dynamics model.

## Related paper

The centroidal dynamics model implemented in dymp is described in the following paper:
[Y. Tazaki, "Trajectory Generation for Legged Robots Based on a Closed-Form Solution of Centroidal Dynamics," in IEEE Robotics and Automation Letters, vol. 9, no. 11, pp. 9239-9246, Nov. 2024]
(https://ieeexplore.ieee.org/abstract/document/10669176).

The whole-body dynamics model is described in the following one:
[Humanoid Dance Simulation Using Choreonoid and Whole-Body Model Predictive Control]
(https://atompc-workshop.github.io/assets/pdf/paper2.pdf).

