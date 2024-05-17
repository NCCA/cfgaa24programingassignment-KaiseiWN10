# Kaisei Wieczorek-Numao (s5507922) - Particle-based fluid simulation - Computing for Graphics and Animation 23/24

The video of my running program can be found in my 'Images+video' folder.

## Introduction
For this project, I have computed a particle-based fluid simulation in C++ utilising the NGL library, offering an efficient simulation framework. This process involved exploring the underlying principles of fluid dynamics, including the use of the Navier-Stokes equation to devise an algorithm to model interactions between particles and external forces, such as gravity in a 3D environment. Furthermore, I have taken the Langrangian approach where each particle holds its own properties such as their position, which consistently gets updated per frame. I decided on taking this route as it simplified the simulation by focusing on individual particles as opposed to solving more complex differential equations. 

## Implementation
Particle system properties:
* Position
* Velocity
* Size
* Colour
* Life
* Mass-Density
* Pressure
* Acceleration

This particle system heavily employs object-oriented programming principles with each particle possesing the attributes listed above. The program revolves around a single central class, 'Emitter', where calculations for these attributes are performed and applied, forming our fluidic motion. Moreover, the render() function stands out as it utilises provided shaders to visually represent the particles onto the screen. Equally important is the update() function which plays an essential role in dynamically refreshing the OpenGL scene. This function manages the lifespan and behaviour of particles within the system, ensuring the particles behave realistically and efficiently, leveraging spatial hashing and updating only alive particles.

![alt tag](http://nccastaff.bournemouth.ac.uk/jmacey/GraphicsLib/Demos/BlankNGL.png)

## Background Research



## Techniques Used

## UML Diagram

![UMLdiagram](https://github.com/NCCA/cfgaa24programingassignment-KaiseiWN10/assets/160144511/3d1b7fa3-899d-42d7-9640-f9700ea25877)

## Conclusion

In conclusion, through the process of completing this assignment I have gained a stronger understanding in the advantages of using object-oriented programming and applying it to creating real world simulations like this particle based fluid system, and I have come to understand that it is a powerful tool to use within computer graphics. To further develop this program, I would focus on improving the accuracy of the fluid motion, and this may involve implementing more particle attributes such as surface-tension. Further, I would also look at ways to optimise this program so that it can run slightly faster. I would try to do that mostly by utilising multithreading to distribute computational tasks across multiple CPU cores, where time consuming procedures such as the 'neighbour' search can be executed in parallel with OpenMP. This would be extremely beneficial as the main issue has been the simulation performance.

## References

Documentation - Particle Simulation using CUDA, S.Green, 2012
https://web.archive.org/web/20140725014123/https://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf

Documentation - Fluid Simulation Using Smoothed Particle Hydrodynamics, A.Strantzi, 2016
https://nccastaff.bournemouth.ac.uk/jmacey/MastersProject/MSc16/15/thesis.pdf

Documentation - Real-Time Fluid Dynamics for Games, J.Stam, Date unspecified 
http://graphics.cs.cmu.edu/nsp/course/15-464/Fall09/papers/StamFluidforGames.pdf

Documentation - Python Dynamic Mode Decomposition, N.Demo, M.Tezzele, G.Rozza, 2018
https://joss.theoj.org/papers/10.21105/joss.00530.pdf

Video - Coding Adventure: Simulating Fluids, S.Lague, 2023
https://www.youtube.com/watch?v=rSKMYc1CQHE&t=2229s

Video - How to write an Eulerian fluid simulator with 200 lines of code, 2022
https://www.youtube.com/watch?v=iKAVRgIrUOU
