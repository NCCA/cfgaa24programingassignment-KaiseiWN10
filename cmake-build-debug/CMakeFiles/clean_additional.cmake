# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Debug")
  file(REMOVE_RECURSE
  "CMakeFiles/ParticleNGL_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/ParticleNGL_autogen.dir/ParseCache.txt"
  "ParticleNGL_autogen"
  )
endif()
