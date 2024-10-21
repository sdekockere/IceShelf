
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ice_shelf_project)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Find ROOT
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake/)
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Find GSL
find_package(GSL REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup ROOT include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include(${ROOT_USE_FILE})
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
set(sources ${PROJECT_SOURCE_DIR}/src/AirToIceRayTracing.cc ${PROJECT_SOURCE_DIR}/src/IceRayTracing.cc ${PROJECT_SOURCE_DIR}/src/IceShelfActionInitialization.cpp ${PROJECT_SOURCE_DIR}/src/IceShelfDetectorConstruction.cpp ${PROJECT_SOURCE_DIR}/src/IceShelfPrimaryGeneratorAction.cpp ${PROJECT_SOURCE_DIR}/src/IceShelfRunAction.cpp ${PROJECT_SOURCE_DIR}/src/IceShelfEventAction.cpp ${PROJECT_SOURCE_DIR}/src/IceShelfSteppingAction.cpp ${PROJECT_SOURCE_DIR}/src/SimplePhysicsList.cpp ${PROJECT_SOURCE_DIR}/src/IceShelfAntennas.cpp)
set(headers ${PROJECT_SOURCE_DIR}/include/AirToIceRayTracing.h ${PROJECT_SOURCE_DIR}/include/IceRayTracing.hh ${PROJECT_SOURCE_DIR}/include/IceShelfActionInitialization.hh ${PROJECT_SOURCE_DIR}/include/IceShelfDetectorConstruction.hh ${PROJECT_SOURCE_DIR}/include/IceShelfPrimaryGeneratorAction.hh ${PROJECT_SOURCE_DIR}/include/IceShelfRunAction.hh ${PROJECT_SOURCE_DIR}/include/IceShelfEventAction.hh ${PROJECT_SOURCE_DIR}/include/IceShelfSteppingAction.hh ${PROJECT_SOURCE_DIR}/include/IceShelfAnalysis.hh ${PROJECT_SOURCE_DIR}/include/SimplePhysicsList.hh ${PROJECT_SOURCE_DIR}/include/IceShelfGlobalVariables.hh ${PROJECT_SOURCE_DIR}/include/IceDensityModels.hh ${PROJECT_SOURCE_DIR}/include/IceShelfAntennas.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ice_shelf ice_shelf.cpp ${sources} ${headers})
target_link_libraries(ice_shelf ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} GSL::gsl GSL::gslcblas)

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(ice_shelf_project DEPENDS ice_shelf)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS ice_shelf DESTINATION bin)
