# Kaisei Wieczorek-Numao (s5507922) CFGAA Assignment

Assignment Idea - 

**PARTICLE FLUID SIMULATION** 

Planning on recreating the dynamic-motion of a fluid through the movement of thousands of particles with individual values in 3D space. For this, I will need to calculate all the particle values for each frame which will be performance heavy so I will also need to consider ways of optimising this algorithm. 

Implement using spacial hashing and Eulerian's equation 

Furthermore, another factor to consider will be finding neighbouring particles and how the forces of them affect the particles we look at in our loop.


DATA STRUCTURES:
- [ ] I will use a 1D array to store data about my particles, each index in the array corresponding to a single particle in my array.
- [ ] A seperate 1D array will need to be created assigned to each particle which will hold the list of neighbouring particles - at the start of each frame this will need clearing and updated. This array will also affect the particle values that are calculated.


CLASSES:
- [ ] Particle - each particle will have it's own values like position, velocity, density etc.
- [ ] Bounding box - Checking the positions of the particles are found within a specified a range, if not implement collision with the bounding walls. 
- [ ] Simulation - Update the values of each particle for each frame