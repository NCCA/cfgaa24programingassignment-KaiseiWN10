# Kaisei Wieczorek-Numao (s5507922) - Computing for Graphics and Animation 23/24
# Particle-based fluid simulation

The video of my running program can be found in my 'Images+video' folder.

## Introduction
For this project, I have computed a particle-based fluid simulation in C++ utilising the NGL and OpenGL libraries, offering an efficient simulation framework. This process involved exploring the underlying principles of fluid dynamics, including the use of the Navier-Stokes equation to devise an algorithm to model interactions between particles and external forces, such as gravity in a 3D environment. Furthermore, I have taken the Langrangian approach where each particle holds its own properties such as their position, which consistently gets updated per frame. I decided on taking this route as it simplified the simulation by focusing on individual particles as opposed to solving more complex differential equations, and also the code provided within the labs was a perfect foundation to build from.

![Screenshot from 2024-05-18 20-04-30](https://github.com/NCCA/cfgaa24programingassignment-KaiseiWN10/assets/160144511/2ba840f9-f69e-414b-9a6f-a0af7679a1fa)  ![Screenshot from 2024-05-18 20-04-35](https://github.com/NCCA/cfgaa24programingassignment-KaiseiWN10/assets/160144511/380f6323-e25c-402f-8220-885afc74cb95)

## Background Research
I have implemented one of the common particle-based methods used to simulate fluids, called 'Smoothed Particle Hydrodynamics', which is an example of a Langrangian method as stated above. Particles in this system have attributes like mass and velocity much like most particle simulations, however, further attributes, including density and pressure, are required to alter the positions to display fluid-like motions. Over a period of time, a set of particles are generated and move according to its assigned calculations until it reaches its life-time, where it is then destroyed.

Another important factor of the SPH method is its use of kernels. Kernels model a delta function which are relative to the particles' positions, and are used for calculating density, pressure and viscosity. For this program I have implemented the spiky, Poly6 (smoothing kernel) and the viscosity kernel and are all used for SPH approxmation but all play a different a different role as briefly stated below:

* Spiky kernel - used for pressure calculation. Produces sharper, stronger gradients than the Poly6 kernel.
* Poly6 kernel - Used for density calculation. Determines the influence of a particle over a certain distance. This decreases with distance meaning the particle that are closer have more of an effect.
* Viscosity kernel - Used for viscosity calculation. Simulates the viscous drag forces that act between neighbouring particles, smoothing out viscosity differences.

## Implementation
Particle system properties:
* Position
* Velocity
* Size
* Colour
* Life
* Mass-Density
* Pressure
* Viscosity
* Acceleration

This particle system heavily employs object-oriented programming principles with each particle possesing the attributes listed above. The program revolves around a single central class, 'Emitter', where calculations for these attributes are performed and applied, forming our fluidic motion. Moreover, the render() function stands out as it utilises provided shaders to visually represent the particles onto the screen. Equally important is the update() function which plays an essential role in dynamically refreshing the OpenGL scene. This function manages the lifespan and behaviour of particles within the system, ensuring the particles behave realistically and efficiently, leveraging spatial hashing and updating only alive particles.

One of the most performance-hindering procedures in this program is very much the neighbour search for each particle which requires updating at the start of every frame due to their position changes. I initially approached this functionality via a simple brute-force method. This is where a list of neighbours can be used as an attribute of the 'Particle' struct and depending on its distance with other particles they would be added to that particle's list. However, I very quickly came to the conclusion that it was an inefficient approach to the problem. As a result, I implemented a method called spatial-hashing where the simulation space is divided into a grid of cells and each particle is assigned a cell based on its position. This prevents the need of checking the distance of a particle with all other particles, only ones that are within distance, going from a time complexity of O(n^2) to O(n).  

                      <img width="214" alt="Screenshot 2024-05-17 at 22 56 57" src="https://github.com/NCCA/cfgaa24programingassignment-KaiseiWN10/assets/160144511/062abf93-dde4-429b-a74d-07da2a6202a6">

# Testing
Through this process, I found that it was very important to find the correct combination of values for my constants and kernel sizes to produce sensible results. Therefore, I began experimenting with different values and explored the effect they had on my overall simulation which was interesting to see. Results of some of my tests are shown below:


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
