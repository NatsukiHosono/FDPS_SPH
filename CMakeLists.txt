
MESSAGE(STATUS "====================================================")
MESSAGE(STATUS "============ Configuring FDPS_SPH ==================")
MESSAGE(STATUS "====================================================")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)

PROJECT(FDPS_SPH)

#Bring the headers, such as Student.h into the project
INCLUDE_DIRECTORIES(src FDPS/src)

SET(DEFAULT_BUILD_TYPE "Release")

IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	MESSAGE(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
	SET(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
ENDIF()

SET(CMAKE_CXX_FLAGS "-DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION -DFAST_COMM_FOR_2ND_EXCHANGE")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops")


MESSAGE(STATUS "CMAKE_BUILD_TYPE:        ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "CMAKE_CXX_FLAGS:         ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "CMAKE_CXX_FLAGS_DEBUG:   ${CMAKE_CXX_FLAGS_DEBUG}")
MESSAGE(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")

ADD_EXECUTABLE(sph.out src/main.cpp)
