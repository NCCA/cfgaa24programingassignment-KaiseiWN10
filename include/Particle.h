#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <ngl/Vec3.h>

struct Particle
{
  Particle()=default;

  Particle(ngl::Vec3 _pos, ngl::Vec3 _dir, int _life, float _size=0.1f, ngl::Vec3 _colour={1,1,1}, float _calculated_properties=0.0f, ngl::Vec3 _gradient={1,1,1},
           ngl::Vec3 _pressure={0,0,0}, float _aver_density=1.0f, ngl::Vec3 _viscosity={1,1,1}, ngl::Vec3 _acceleration={1,1,1}, ngl::Vec3 _buoyancy={0,0,0}) :
      pos{_pos},dir{_dir},life{_life},size{_size},colour{_colour},pressure{_pressure},aver_density{_aver_density},
      viscosity{_viscosity},acceleration{_acceleration},buoyancy{_buoyancy}
      {

      }
  ngl::Vec3 pos;
  ngl::Vec3 dir;
  ngl::Vec3 colour;
  int life=100;
  float size=0.01f;
  bool isAlive = false;
  ngl::Vec3 pressure;
  ngl::Vec3 viscosity;
  ngl::Vec3 acceleration;
  ngl::Vec3 buoyancy;
  float density = 0.0f;
  float aver_density = 1.0f;
};

#endif
