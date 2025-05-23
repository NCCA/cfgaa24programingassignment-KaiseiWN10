#include "Emitter.h"
#include <ngl/Random.h>
#include <fmt/format.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <ngl/Util.h>
#include <ngl/Transformation.h>
#include <ngl/ShaderLib.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/VAOFactory.h>
#include <omp.h>
#include <random>
#include <limits>
Emitter::Emitter(int _numParticles, int _maxAlive)
{
    m_numParticles=_numParticles;
    m_maxAlive=_maxAlive;
    pos.resize(m_numParticles);
    dir.resize(m_numParticles);
    colour.resize(m_numParticles);
    life.resize(m_numParticles);
    isAlive.resize(m_numParticles);
    densities.resize(m_numParticles);
    pressure.resize(m_numParticles);
    viscosity.resize(m_numParticles);
    calculated_properties.resize(m_numParticles);
    aver_density.resize(m_numParticles);
    acceleration.resize(m_numParticles);
    buoyancy.resize(m_numParticles);

    for(size_t i=0; i<m_numParticles; ++i)
    {
        createZeroParticle(i);
    }

    m_vao = ngl::vaoFactoryCast<ngl::MultiBufferVAO>(ngl::VAOFactory::createVAO(ngl::multiBufferVAO,GL_POINTS));
    m_vao->bind();
    m_vao->setData(ngl::MultiBufferVAO::VertexData(0,0)); // positions
    m_vao->setData(ngl::MultiBufferVAO::VertexData(1,0)); // colours
    m_vao->unbind();
}

ngl::Vec3 getRandomDir() {
    // Seed the random number generator
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(-1.0f, 1.0f);

    // Generate random values for x and y
    float randomX = dis(gen);
    float randomY = dis(gen);

    // Normalize the vector
    float length = sqrt(randomX * randomX + randomY * randomY);
    ngl::Vec3 dir;
    dir.m_x = randomX / length;
    dir.m_y = randomY / length;

    return dir;
}

void Emitter::createZeroParticle(size_t _p)
{
//    pos[_p].set(0.0f,0.0f,0.0f,0.01f);
//    dir[_p].set(0.0f,0.0f,0.0f);
    isAlive[_p]=false;
}

void Emitter::createDefaultParticle(size_t _p)
{
    pos[_p]=m_position;
    dir[_p]=m_emitDir * ngl::Random::randomPositiveNumber() +randomVectorOnSphere() * m_spread;
    dir[_p].m_y = std::abs(dir[_p].m_y);
    colour[_p] = {0.0f, 0.0f, 1.0f, 1.0f}; //setting initial colour to blue
    life[_p] = static_cast<int>(2.0f+ngl::Random::randomPositiveNumber(150));
    pos[_p].m_w= 0.01f;
    isAlive[_p] = true;
}

ngl::Vec3 Emitter::randomVectorOnSphere()
{
    auto phi = ngl::Random::randomPositiveNumber(M_PI * 2.0f);
    auto costheta = ngl::Random::randomNumber();
    auto theta = acosf(costheta);
    auto u = ngl::Random::randomPositiveNumber();
    auto r = 2.0f *std::cbrt(u);
    return ngl::Vec3(r*sinf(theta) *cosf(phi),
                     r* sinf(theta) * sinf(phi),
                     r*cosf(theta));
}

void Emitter::render() const
{
    glEnable(GL_PROGRAM_POINT_SIZE);
    //glPointSize(4.0f);
    m_vao->bind();
    m_vao->setData(0,ngl::AbstractVAO::VertexData(m_numParticles*sizeof(ngl::Vec4),pos[0].m_x));
    m_vao->setVertexAttributePointer(0,4,GL_FLOAT,0,0);

    m_vao->setData(1,ngl::AbstractVAO::VertexData(m_numParticles*sizeof(ngl::Vec3),colour[0].m_x));
    m_vao->setVertexAttributePointer(1,3,GL_FLOAT,0,0);
    m_vao->setNumIndices(m_numParticles);
    m_vao->draw();
    m_vao->unbind();
    glDisable(GL_PROGRAM_POINT_SIZE);
}

float SmoothingKernel(float radius, float dist)//Poly6 kernel
{
    if(dist < radius)
    {
        float scale = 315 / (64*M_PI * std::pow(abs(radius), 9)); //Poly6 kernel
        float v = radius*radius - dist*dist;
        return v*v*v*scale;
    }
    return 0;
}

float Emitter::magnitude(ngl::Vec4 currentParticle) // finding the magnitude of a passed in position
{
    return (float)std::sqrt(std::pow(currentParticle.m_x,2) +
                  std::pow(currentParticle.m_y,2) +
                  std::pow(currentParticle.m_z,2));
}

float Emitter::calcDensity(size_t currentParticle) //density calculation
{
    float density = 0;

    for(int p=0; p<pos.size(); p++) //going through each particle
    {
        float dist = magnitude(pos[p] - pos[currentParticle]);
        float influence = SmoothingKernel(1.0f, dist);
        density += 1.0f * influence; //1.0f SHOULD BE REPLACED BY MASS
    }
    return density;
}

float Emitter::averageNeighDensity(size_t currentParticle)//finding the average of a particle's neighbouring densities
{
    float totalDensity = 0.0f;
    int count = 0;

    ngl::Vec4 currentCell = positionToCellCoord(pos[currentParticle]); // Gets the cell coordinate of the current particle

    for (int dx = -1; dx <= 1; ++dx) // Iterates over the neighboring cells
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {
                ngl::Vec4 neighborCell(currentCell.m_x + dx, currentCell.m_y + dy, currentCell.m_z + dz, 0.0f);// Calculates the neighboring cell coordinate
                uint key = getKeyfromHash(hashCell(static_cast<int>(neighborCell.m_x), static_cast<int>(neighborCell.m_y), static_cast<int>(neighborCell.m_z)));// Retrieves the key for the neighbouring cell
                int cellStartIndex = startIndices[key];// Retrieve the start index of particles in the neighboring cell

                for (int i = cellStartIndex; i < pos.size(); ++i) // Iterate over particles in the neighboring cell
                {
                    if (spatialLookUp[i] != key) break; // End of neighboring cell

                    float dist = magnitude(pos[i] - pos[currentParticle]);// Calculate distance between current particle and neighboring particle

                    if (dist <= 1.0f)// Check if neighboring particle is within a certain radius
                    {
                        totalDensity += densities[i];
                        ++count;
                    }
                }
            }
        }
    }
    float averageDensity = (count > 0) ? totalDensity / count : 0.0f;// Calculate average density

    return averageDensity;
}

void Emitter::debugPressure(ngl::Vec4 dist)
{
    ngl::Vec3 normalised = {dist.m_x/ magnitude(dist),dist.m_y/ magnitude(dist),dist.m_z/ magnitude(dist)};
    ngl::Vec3 denominator = (M_PI * std::pow(1.0f, 6)) * normalised * std::pow((1.0f - magnitude(dist)), 2);
    //std::cout<<"denominator is: "<<"("<<denominator.m_x<<","<<denominator.m_y<<","<<denominator.m_z<<")"<<std::endl;
    //std::cout<<"magnitude is: "<<magnitude(dist)<<std::endl;
    //std::cout<<"dist is: "<<"("<<dist.m_x<<","<<dist.m_y<<","<<dist.m_z<<")"<<std::endl;
    //std::cout<<"normalized is: "<<"("<<normalised.m_x<<","<<normalised.m_y<<","<<normalised.m_z<<")"<<std::endl;    ;
}

ngl::Vec3 Emitter::getPressureKernel(ngl::Vec4 dist) // finding the pressure kernel used to calculate the particle pressure
{
    ngl::Vec3 normalised = {dist.m_x/ magnitude(dist),dist.m_y/ magnitude(dist),dist.m_z/ magnitude(dist)};
    ngl::Vec3 denominator = (M_PI * std::pow(1.0f, 6)) * normalised * std::pow((1.0f - magnitude(dist)), 2);

    if(denominator.m_x==0.0f||denominator.m_y==0.0f||denominator.m_z==0.0f)
    {
        ngl::Vec3 pressKernel = {0.0f, 0.0f, 0.0f};
        return pressKernel;
    }
    if(std::isnan(denominator.m_x)||std::isnan(denominator.m_y)||std::isnan(denominator.m_z))
    {
        ngl::Vec3 pressKernel = {0.0f, 0.0f, 0.0f};
        return pressKernel;
    }
    else
    {
        ngl::Vec3 pressKernel = {-45.0f / denominator.m_x, -45.0f / denominator.m_y, -45.0f / denominator.m_z};
        return pressKernel;
    }
}

ngl::Vec3 Emitter::calcPressure(size_t currentParticle) //pressure calculation
{
    ngl::Vec4 currentCell = positionToCellCoord(pos[currentParticle]); // Gets the cell coordinate of the current particle
    pressure[currentParticle] = ngl::Vec3(0.0f, 0.0f, 0.0f);

    for (int dx = -1; dx <= 1; ++dx) // Iterates over the neighboring cells
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {
                ngl::Vec4 neighborCell(currentCell.m_x + dx, currentCell.m_y + dy, currentCell.m_z + dz, 0.0f);// Calculates the neighbouring cell coordinate
                uint key = getKeyfromHash(hashCell(static_cast<int>(neighborCell.m_x), static_cast<int>(neighborCell.m_y),static_cast<int>(neighborCell.m_z)));// Retrieves the key for the neighbouring cell
                int cellStartIndex = startIndices[key];// Retrieves the start index of particles in the neighbouring cell

                for (int i = cellStartIndex; i < pos.size(); ++i) // Iterate over particles in the neighboring cell
                {
                    if (spatialLookUp[i] != key)
                    {
                        break; // End of neighboring cell
                    }

                    ngl::Vec4 dist = pos[currentParticle] - pos[i];// Calculate distance between current particle and neighboring particle

                    if(std::isnan(densities[i])||densities[i]==0.0f)
                    {
                        pressure[currentParticle].m_x += 0.0f;
                    }
                    else
                    {
                        pressure[currentParticle].m_x += ((pressure[i].m_x + pressure[currentParticle].m_x)/2) * (1.0f/densities[i]) * getPressureKernel(dist).m_x;//REPLACE 1.0F BY MASS
                        pressure[currentParticle].m_y += ((pressure[i].m_y + pressure[currentParticle].m_y)/2) * (1.0f/densities[i]) * getPressureKernel(dist).m_y;//REPLACE 1.0F BY MASS
                        pressure[currentParticle].m_z += ((pressure[i].m_z + pressure[currentParticle].m_z)/2) * (1.0f/densities[i]) * getPressureKernel(dist).m_z;//REPLACE 1.0F BY MASS
                    }
                }
            }//end of dz loop
        }//end of dy loop
    }//end of dx loop

    pressure[currentParticle] *= -1;
    return pressure[currentParticle];
}

float getViscKernel(float dist) // finding the viscosity kernel used to calculate the particle viscosity
{
    float viscKernel = 45.0f / (M_PI * std::pow(1.0f, 6)) * (1.0f - dist);
    return viscKernel;
}

ngl::Vec3 Emitter::calcViscosity(size_t currentParticle) // calculating viscosity of the passed in particle
{
    float const vConst = 3.5f;
    viscosity[currentParticle] = ngl::Vec3(0.0f, 0.0f, 0.0f);

    ngl::Vec4 currentCell = positionToCellCoord(pos[currentParticle]); // Gets the cell coordinate of the current particle

    for (int dx = -1; dx <= 1; ++dx) // Iterates over the neighbouring cells
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {
                ngl::Vec4 neighborCell(currentCell.m_x + dx, currentCell.m_y + dy, currentCell.m_z + dz, 0.0f);// Calculates the neighbouring cell coordinate
                uint key = getKeyfromHash(hashCell(static_cast<int>(neighborCell.m_x), static_cast<int>(neighborCell.m_y),static_cast<int>(neighborCell.m_z)));// Retrieves the key for the neighbouring cell
                int cellStartIndex = startIndices[key];// Retrieves the start index of particles in the neighbouring cell

                for (int i = cellStartIndex; i < pos.size(); ++i) // Iterate over particles in the neighbouring cell
                {
                    if (spatialLookUp[i] != key) break; // End of neighboring cell

                    float dist = magnitude(pos[i]-pos[currentParticle]);// Calculate distance between current particle and neighbouring particle
                    ngl::Vec3 dist2 = dir[i]-dir[currentParticle];// Calculate difference in velocity between current particle and neighbouring particle
                    float densityDenom = (1.0f / densities[i]);

                    if(std::isnan(densityDenom)||densityDenom==0)
                    {
                        densityDenom = 0;
                    }

                    viscosity[currentParticle].m_x += dist2.m_x*densityDenom*getViscKernel(dist);//REPLACE 1.0F BY MASS
                    viscosity[currentParticle].m_y += dist2.m_y*densityDenom*getViscKernel(dist);
                    viscosity[currentParticle].m_z += dist2.m_z*densityDenom*getViscKernel(dist);
                }
            }//end of dz loop
        }//end of dy loop
    }//end of dx loop
    viscosity[currentParticle] *= vConst;
    return viscosity[currentParticle];
}

ngl::Vec3 Emitter::calcBuoyancy(size_t currentParticle, ngl::Vec3 gravity) // calculating buoyancy for the passed in particle
{
    float const buoyConst = 0.0f;
    float const restDensity = 998.2f;

    return buoyConst * (densities[currentParticle] - restDensity) * gravity;
}

ngl::Vec4 Emitter::positionToCellCoord(ngl::Vec4 point) //converting particle position to it's cell coordinate
{
    ngl::Vec4 result;
    int cellX = static_cast<int>(point.m_x / 1.0f);
    int cellY = static_cast<int>(point.m_y / 1.0f);
    int cellZ = static_cast<int>(point.m_z / 1.0f);
    result.m_x = static_cast<float>(cellX);
    result.m_y = static_cast<float>(cellY);
    result.m_z = static_cast<float>(cellZ);
    result.m_w = 0.01f;
    return result;
}

uint Emitter::hashCell(int cellX, int cellY, int cellZ)//converts cell-coordinate into a single number
{
    uint a = (uint)cellX * 15823;
    uint b = (uint)cellY * 9737333;
    uint c = (uint)cellZ * 2345;
    return a+b+c;
}

uint Emitter::getKeyfromHash(uint hash)//wraps the hash value around the length of the array, so it can be as an index
{
    return hash % pos.size();
}

void Emitter::updateStartIndices() //updating startIndices array
{
    int currentStartIndex = -1;
    uint prevKey = std::numeric_limits<uint>::max();

    for (size_t p = 0; p < pos.size(); ++p)
    {
        uint currentKey = spatialLookUp[p];
        if (currentKey != prevKey)
        {
            startIndices[currentKey] = currentStartIndex + 1;
            currentStartIndex = p;
        }
        prevKey = currentKey;
    }
}

void Emitter::updateSpatialLookup() //updating spatialLookup array used for spatial-hashing
{
    spatialLookUp.resize(pos.size());
    startIndices.resize(pos.size());

    std::fill(startIndices.begin(), startIndices.end(), -1);

    int currentStartIndex = -1;

    for (size_t p = 0; p < pos.size(); ++p)
    {
        uint currentKey = spatialLookUp[p];
        spatialLookUp[p] = getKeyfromHash(hashCell(pos[p].m_x, pos[p].m_y, pos[p].m_z));

        if (currentKey != spatialLookUp[p])
        {
            startIndices[currentKey] = currentStartIndex + 1;
            currentStartIndex = p;
        }
    }
    updateStartIndices();
}

ngl::Vec3 Emitter:: updateVelocity(size_t currentParticle) //updating velocity values using calculated density values
{
    ngl::Vec4 currentCell = positionToCellCoord(pos[currentParticle]); // Gets the cell coordinate of the current particle

    float const h = 1.0f;
    float nei = 0.0f;
    float const velConst = 0.1f;

    for (int dx = -1; dx <= 1; ++dx) // Iterates over the neighboring cells
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {
                ngl::Vec4 neighborCell(currentCell.m_x + dx, currentCell.m_y + dy, currentCell.m_z + dz, 0.0f);// Calculates the neighboring cell coordinate
                uint key = getKeyfromHash(hashCell(static_cast<int>(neighborCell.m_x), static_cast<int>(neighborCell.m_y),static_cast<int>(neighborCell.m_z)));// Retrieves the key for the neighbouring cell
                int cellStartIndex = startIndices[key];// Retrieves the start index of particles in the neighboring cell
                for (int i = cellStartIndex; i < pos.size(); ++i) // Iterate over particles in the neighboring cell
                {
                    if (spatialLookUp[i] != key) break; // End of neighboring cell

                    float dist = magnitude(pos[i]-pos[currentParticle]);// Calculate distance between current particle and neighboring particle

                    float visc_Kernel = 45.0f / (M_PI * std::pow(h, 6)) * (0.3f - dist);
                    nei += ((2*1.0f)/(densities[currentParticle] + densities[i])) * SmoothingKernel(dist,h); //REPLACE 1.0 with MASS
                }
            }//end of dz loop
        }//end of dy loop
    }//end of dx loop
    dir[currentParticle].m_x += velConst + nei;
    dir[currentParticle].m_y += velConst + nei;
    dir[currentParticle].m_z += velConst + nei;

    return dir[currentParticle];
}

void Emitter:: checkBoundary(size_t currentParticle) //Ensuring the particles remain within the specified boundary
{
    float const collisionDamping = 0.6f; // used to reduce the velocity of the particle when colliding with a boundary
    ngl::Vec3 const boundary = {6.0f,6.0f,6.0f};//specified boundary

    if (abs(pos[currentParticle].m_x) > boundary.m_x || abs(pos[currentParticle].m_x) < -boundary.m_x)//checking x-component of particle
    {
        pos[currentParticle].m_x *= -1;
        dir[currentParticle].m_x *= -1 * collisionDamping;
    }
    if (abs(pos[currentParticle].m_y) > boundary.m_y || abs(pos[currentParticle].m_y) < -boundary.m_y)//checking y-component of particle
    {
        pos[currentParticle].m_y *= -1;
        dir[currentParticle].m_y *= -1 * collisionDamping;
    }
    if (abs(pos[currentParticle].m_z) > boundary.m_z || abs(pos[currentParticle].m_z) < -boundary.m_z)//checking z-component of particle
    {
        pos[currentParticle].m_z *= -1;
        dir[currentParticle].m_z *= -1 * collisionDamping;
    }
}

void Emitter::updateColour(size_t currentParticle)
{
    if(dir[currentParticle].m_x + dir[currentParticle].m_y + dir[currentParticle].m_z) // Increasing the r,g values of those with higher overall velocity
    {
        colour[currentParticle].m_r += 0.4f;
        colour[currentParticle].m_b += 0.3f;
    }
}

void Emitter::update()
{
    float _dt = 0.02f;

    ngl::Vec3 gravity(0.0f, -1.0f, 0.0f);

    // Choose number to birth
    int numberToBirth = 1000 + ngl::Random::randomPositiveNumber(1500);


    for (int i = 0; i < numberToBirth; ++i)
    {
        size_t index = 0;
        for (auto a : isAlive)
        {
            if (!a)
                break;
            else
                ++index;
        }
        createDefaultParticle(index);
    }

    float const accConst = 5.3f;

    for (size_t p = 0; p < m_numParticles; ++p)// Calculate dir[] for all particles
    {
        if (isAlive[p])
        {
            dir[p] += gravity * _dt * 0.5f; // Calculating velocity
        }
    }

    updateSpatialLookup();//spatial hashing - finding the neighbours of each particle

    // Calculate densities[] for all particles
    for (size_t p = 0; p < m_numParticles; ++p)
    {
        if (isAlive[p])
        {
            densities[p] = calcDensity(p); // Calculate average density for the current particle
        }
    }

    for (size_t p = 0; p < m_numParticles; ++p)
    {
        if (isAlive[p])
        {
            aver_density[p] = averageNeighDensity(p); // Calculate average density for the current particle
            pressure[p] = calcPressure(p);// calculate pressure for each particle
        }
    }

    for (size_t p = 0; p < m_numParticles; ++p)
    {
        viscosity[p] = calcViscosity(p);// calculate viscosity for each particle
        buoyancy[p] = calcBuoyancy(p, gravity);// calculate buoyancy for each particle
    }

    for (size_t p = 0; p < m_numParticles; ++p)
    {
        if (isAlive[p])
        {
            //calculates acceleration for each particle - adds all forces as stated in Navier-Stokes eq.
            acceleration[p].m_x = pressure[p].m_x + viscosity[p].m_x + accConst + buoyancy[p].m_x + gravity.m_x;
            acceleration[p].m_y = pressure[p].m_y + viscosity[p].m_y + accConst + buoyancy[p].m_y + gravity.m_y;
            acceleration[p].m_z = pressure[p].m_z + viscosity[p].m_z + accConst + buoyancy[p].m_z + gravity.m_z;
            dir[p] = updateVelocity(p);
        }
    }

    // Update pos[] for all particles
    for (size_t p = 0; p < m_numParticles; ++p)
    {
        if (isAlive[p])
        {
            // Update position
            pos[p] += dir[p] * _dt * acceleration[p];
            pos[p].m_w = 0.3f;
            updateColour(p);
            checkBoundary(p);

            if (--life[p] == 0 || pos[p].m_y <= 0.0f) //checking if the position of our particles are within our bounding box
            {
                createZeroParticle(p);
            }
        }
    }
}
