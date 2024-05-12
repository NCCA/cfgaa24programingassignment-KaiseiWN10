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
    calculated_properties.resize(m_numParticles);
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
    colour[_p] = ngl::Random::getRandomColour3();
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
    auto r = 1.0f *std::cbrt(u);
    return ngl::Vec3(r*sinf(theta) *cosf(phi),
                     r* sinf(theta) * sinf(phi),
                     r*cosf(theta));
}

void Emitter::render() const
{
    glEnable(GL_PROGRAM_POINT_SIZE);
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

float SmoothingKernel(float radius, float dist)
{
    if(dist < radius)
    {
        float scale = 315 / (64*M_PI * std::pow(abs(radius), 9));
        float v = radius*radius - dist*dist;
        return v*v*v*scale;
    }
    return 0;
}

float SpikyKernelPow3(float radius, float dist)
{
    if(dist < radius)
    {
        float scale = 15 / (M_PI * std::pow(radius, 6));
        float v = radius - dist;
        return v*v*v*v*scale;
    }
    return 0;
}

float SpikyKernelPow2(float radius, float dist)
{
    if(dist < radius)
    {
        float scale = 15 / (2 * M_PI * std::pow(radius, 5));
        float v = radius - dist;
        return v*v*v*v*scale;
    }
    return 0;
}

float SpikyPow3_Deriv(float radius, float dist)
{
    if(dist <= radius)
    {
        float scale = 45 / (pow(radius,6) * M_PI);
        float v = radius - dist;
        return  -v*v*scale;
    }
    return 0;
}

float SpikyPow2_Deriv(float radius, float dist)
{
    if(dist <= radius)
    {
        float scale = 15 / (pow(radius,5) * M_PI);
        float v = radius - dist;
        return  -v*scale;
    }
    return 0;
}

float densityKernel(float radius, float dist)
{
    return SpikyKernelPow2(radius, dist);
}

float nearDensityKernel(float radius, float dist)
{
    return SpikyKernelPow3(radius, dist);
}

float density_Deriv(float radius, float dist)
{
    return SpikyPow2_Deriv(radius, dist);
}

float nearDensity_Deriv(float radius, float dist)
{
    return SpikyPow3_Deriv(radius, dist);
}

static float SmoothingKernel_Deriv(float radius, float dist)
{
    if(dist>=radius) return 0;
    float scale = 12 / (std::pow(radius,4) * M_PI);
    return (dist-radius) * scale;
}

float Emitter::magnitude(ngl::Vec4 currentParticle)
{
    return std::sqrt(currentParticle.m_x*currentParticle.m_x +
                             currentParticle.m_y*currentParticle.m_y +
                             currentParticle.m_z*currentParticle.m_z);
}

float Emitter::calcDensity(size_t currentParticle)
{
    float density = 0;

    for(int p=0; p<pos.size(); p++) //going through each particle
    {
        //std::cout<<density << " ";
        //std::cout<<"size of particle array: " << sizeof(pos) << std::endl;
        //std::cout<<"currentpos " << pos[currentParticle].m_x << std::endl;
        float dist = magnitude(pos[p] - pos[currentParticle]);
        //std::cout<<dist << " ";
        float influence = SmoothingKernel(1.0f, dist);
        //std::cout<<influence<<" ";
        density += 1.0f * influence; //1.0f SHOULD BE REPLACED BY MASS
    }
    //std::cout<<density<<" ";
    return density;
}

void Emitter::updateDensities() //for optimisation - so it doesn't have to loop through all the particles to calc. density every time - just use cache values
//reduce computation time
{
#pragma omp parallel for
    for(int p=0; p<pos.size(); p++) //going through each particle
    {
        densities[p] = calcDensity(p);
    }
}

float Emitter::calcProperty(size_t currentParticle)
{
    float property = 0;

    for(int p=0; p<pos.size(); p++) //going through each particle
    {
        float dist = magnitude(pos[p] - pos[currentParticle]);
        float influence = SmoothingKernel(1.0f, dist);
        float density = calcDensity(p);
        property += -calculated_properties[p] * influence * 1.0f / density; //1.0f needs replacing with mass
    }
    return property;
}

std::vector<ngl::Vec3> Emitter::calcPropertyGradient(size_t currentParticle)
{
    for(int p=0; p<pos.size(); p++) //going through each particle
    {
        //gradient[p];
        float dist = magnitude(pos[p] - pos[currentParticle]);
        ngl::Vec4 dir = (pos[p] - pos[currentParticle]) / dist;
        float slope = SmoothingKernel_Deriv(1.0f, dist);//REPLACE 1.0f WITH SMOOTHING-RADIUS
        float density = densities[p];
        gradient[p] += -calculated_properties[p] * ngl::Vec3(dir.m_x, dir.m_y, dir.m_z) * slope * 1.0f / density;//REPLACE 1.0f WITH MASS
    }
    return gradient;
}

float Emitter::convertDensitytoPressure(float density)
{
    float densityError = density - targetDensity;//how far density is from what we want it to be
    float pressure = densityError * pressureMultiplier;
    return pressure;
}

float Emitter::calcSharedPressure(float densityA, float densityB)
{
    float pressureA = convertDensitytoPressure(densityA);
    float pressureB = convertDensitytoPressure(densityB);
    return (pressureA + pressureB) / 2;
}

ngl::Vec3 Emitter::calcPressure(size_t currentPos)
{
    ngl::Vec3 pressureForce;

    for(int p=0; p<pos.size(); p++) //going through each particle
    {
        if(currentPos == p) continue;

        ngl::Vec4 offset = pos[p] - pos[currentPos];
        float dist = magnitude(offset);
        ngl::Vec4 dir = dist == 0 ? getRandomDir() : offset / dist;
        float slope = SmoothingKernel_Deriv(1.0f, dist);//REPLACE 1.0f WITH SMOOTHING-RADIUS
        float density = densities[p];
        float sharedPressure = calcSharedPressure(density, densities[p]);
        pressureForce += sharedPressure * ngl::Vec3(dir.m_x, dir.m_y, dir.m_z) * slope * 1.0f /density;//REPLACE 1.0f WITH MASS
    }
    return pressureForce;
}

ngl::Vec4 Emitter::positionToCellCoord(ngl::Vec4 point)
{
    ngl::Vec4 result;
    int cellX = static_cast<int>(point.m_x / 1.0f);
    int cellY = static_cast<int>(point.m_y / 1.0f);
    int cellZ = static_cast<int>(point.m_z / 1.0f);
    result.m_x = static_cast<float>(cellX);
    result.m_y = static_cast<float>(cellY);
    result.m_z = static_cast<float>(cellZ);
    result.m_w = 0.0f;
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

void Emitter:: updateSpatialLookup(float radius)
{
    spatialLookUp.resize(pos.size());
    startIndices.resize(pos.size());
   for(size_t p=0; p<pos.size(); p++)//creates an unordered spatial-lookup
   {
       ngl::Vec4 cellCoord = positionToCellCoord(pos[p]);
       uint cellHash = hashCell(cellCoord.m_x,cellCoord.m_y,cellCoord.m_z);
       uint cellKey = getKeyfromHash(cellHash);
       spatialLookUp[p] = cellKey;
       startIndices[p] = std::numeric_limits<int>::max();
   }
//    {

//        startIndices[p] = std::numeric_limits<int>::max(); //reset start index
//    }
//
//    std::sort(spatialLookUp.begin(), spatialLookUp.end()); //sort by cell key
//
//    for(int p=0; p<points.size(); p++)//calculate start indices of each unique cell key in the spatial-lookup
//    {
//        uint key = spatialLookUp[p];
//        uint keyPrev = (p == 0) ? std::numeric_limits<uint>::max() : spatialLookUp[p - 1];
//        if(key != keyPrev)
//        {
//            startIndices[key] = p;
//        }
//    }
//
//    startIndices[0] = 0;

}

//void Emitter:: forEachPointWithinRadius(ngl::Vec3 currentPoint)
//{
//    ngl::Vec3 Centre = positionToCellCoord(currentPoint);
//    float sqrRadius = 1.0f * 1.0f; //1.0f MUST BE REPLACED BY RADIUS
//
//    //foreach((3,3,3) in cellOffsets)
//
//    uint key = getKeyfromHash(hashCell(Centre.m_x+Centre.m_x,Centre.m_y+Centre.m_y,Centre.m_z+Centre.m_z));
//    int cellStartIndex = startIndices[key];
//
//    for(int i = cellStartIndex; i < spatialLookUp.size(); i++)
//    {
//        if(spatialLookUp[i] != key) break;
//        float sqrDist = magnitude(points[i] - currentPoint);
//        if(sqrDist <= sqrRadius)
//        {
//
//        }
//    }
//}

void Emitter::update()
{
    float _dt=0.1f;
    ngl::Vec3 gravity(0,-9.87, 0);
//// choose number to birth
// find first not alive and set as new particle
    int numberToBirth=1000+ngl::Random::randomPositiveNumber(1500);

    for(int i=0; i<numberToBirth; ++i)
    {
        size_t index=0;
        for(auto a : isAlive )
        {
            if (!a)
                break;
            else
                ++index;
        }
        createDefaultParticle(index);
    }

    float collisionDamping = 0.35f;

    for(size_t p=0; p<m_numParticles; ++p)//loop through all particles
    {
        if(isAlive[p])
        {
            //for(size_t p=0; p<m_numParticles; ++p)
            //{
            updateSpatialLookup(1.0f);
                dir[p] += gravity * _dt *0.5f;//calculating velocity
                //std::cout << p << " ";
                //predictedPos[p] = pos[p] + dir[p] * _dt;
            //}

            for(size_t p=0; p<m_numParticles; ++p)
            {
                //densities[p] = calcDensity(p);//calculate densities
                //std::cout << densities[p] << " ";
            }

            for(size_t p=0; p<m_numParticles; ++p)
            {
                //pressure[p] = calcPressure(p);//calculating pressure-force
                //std::cout << pressure[p].m_x << " ";
                //pressure_Acc[p] = pressure[p] / densities[p];//calculating pressure-acceleration
                //dir[p] = pressure_Acc[p] * _dt;//updating velocity
            }

            for(size_t p=0; p<m_numParticles; ++p)
            {
                //std::cout << dir[p].m_y << " ";
                pos[p] += dir[p] * _dt;//updating position

                pos[p].m_w += 0.1f;

                //std::cout<<"size of particle array: " << sizeof(pos);
                if(--life[p] ==0 || pos[p].m_y <=0.0f)
                {
                    createZeroParticle(p);
                }
                if (abs(pos[p].m_x) > 20)
                {
                    pos[p].m_x *= -1;
                    dir[p].m_x *= -1 * collisionDamping;
                }
                if (abs(pos[p].m_y) > 20)
                {
                    pos[p].m_y *= -1;
                    dir[p].m_y *= -1 * collisionDamping;
                }

                if (abs(pos[p].m_z) > 20)
                {
                    pos[p].m_z *= -1;
                    dir[p].m_z *= -1 * collisionDamping;
                }
            }


        }

    }
}