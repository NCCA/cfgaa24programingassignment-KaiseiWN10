#ifndef EMITTER_H_
#define EMITTER_H_

#include <vector>
#include <string_view>
#include <ngl/Vec3.h>
#include <ngl/MultiBufferVAO.h>
#include <memory>
class Emitter
{
public :
    Emitter(int _numParticles, int _maxAlive);
    float calcDensity(size_t currentPos);
    float averageNeighDensity(size_t currentParticle);
    float magnitude(ngl::Vec4 currentParticle);
    ngl::Vec3 getPressureKernel(ngl::Vec4 dist);
    ngl::Vec3 calcPressure(size_t currentPos);
    ngl::Vec3 calcViscosity(size_t currentParticle);
    ngl::Vec3 calcBuoyancy(size_t currentParticle, ngl::Vec3 gravity);
    void updateSpatialLookup();
    ngl::Vec4 positionToCellCoord(ngl::Vec4 point);
    uint hashCell(int cellX, int cellY, int cellZ);
    uint getKeyfromHash(uint hash);
    void updateStartIndices();
    ngl::Vec3 updateVelocity(size_t currentParticle);
    void checkBoundary(size_t currentParticle);
    void updateColour(size_t currentParticle);
    void update();
    void render() const;

private :
    ngl::Vec3 randomVectorOnSphere();
    void createDefaultParticle(size_t _p);
    void createZeroParticle(size_t _p);
    size_t m_numParticles;
    ngl::Vec3 m_position={0,0,0};
    float m_spread = 15.0f;
    ngl::Vec3 m_emitDir = {0,20.0f,0};
    int m_maxAlive;
    std::unique_ptr<ngl::MultiBufferVAO> m_vao;

    std::vector<ngl::Vec4> pos;
    std::vector<ngl::Vec3> dir;
    std::vector<ngl::Vec3> colour;
    std::vector<int> life;
    std::vector<float> calculated_properties;
    std::vector<float> densities;
    std::vector<float> aver_density;
    std::vector<ngl::Vec3> buoyancy;
    std::vector<bool> isAlive;
    std::vector<ngl::Vec3> gradient;
    std::vector<ngl::Vec3> pressure;
    std::vector<ngl::Vec3> pressure_Acc;
    std::vector<uint> spatialLookUp;
    std::vector<uint> startIndices;
    std::vector<ngl::Vec4> points;
    std::vector<ngl::Vec3> predictedPos;
    std::vector<ngl::Vec3> viscosity;
    std::vector<ngl::Vec3> acceleration;
};

#endif
