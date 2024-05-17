# Kaisei Wieczorek-Numao (s5507922) - Particle-based fluid simulation - Computing for Graphics and Animation 23/24

The video of my running program can be found in my 'Images+video' folder.

For this project, I have computed a particle-based fluid simulation in C++ utilising the NGL library, offering an efficient simulation framework. This process involved exploring the underlying principles of fluid dynamics, including the use of the Navier-Stokes equation to devise an algorithm to model interactions between particles and external forces, such as gravity in a 3D environment. Furthermore, I have taken the Langrangian approach where each particle holds its own properties such as their position, which consistently gets updated per frame. I decided on taking this route as it simplified the simulation by focusing on individual particles as opposed to solving complex differential equations. 

* Assign variables to each particle

![alt tag](http://nccastaff.bournemouth.ac.uk/jmacey/GraphicsLib/Demos/BlankNGL.png)

## Background Research

## Techniques Used

## UML Diagram

![UMLdiagram](https://github.com/NCCA/cfgaa24programingassignment-KaiseiWN10/assets/160144511/3d1b7fa3-899d-42d7-9640-f9700ea25877)

## Conclusion

## References

Documentation - Particle Simulation using CUDA, S.Green, 2012
https://web.archive.org/web/20140725014123/https://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf

Documentation - Fluid Simulation Using Smoothed Particle Hydrodynamics, A.Strantzi, 2016
https://nccastaff.bournemouth.ac.uk/jmacey/MastersProject/MSc16/15/thesis.pdf

Documentation - Real-Time Fluid Dynamics for Games, J.Stam, Date unspecified 
http://graphics.cs.cmu.edu/nsp/course/15-464/Fall09/papers/StamFluidforGames.pdf

Video - Coding Adventure: Simulating Fluids, S.Lague, 2023
https://www.youtube.com/watch?v=rSKMYc1CQHE&t=2229s

Video - How to write an Eulerian fluid simulator with 200 lines of code, 2022
https://www.youtube.com/watch?v=iKAVRgIrUOU
