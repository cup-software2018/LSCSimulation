cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(LSCSimulation)
set(CMAKE_CXX_STANDARD 17)
#
# Locate the ROOT package
#
find_package(ROOT)
include(${ROOT_USE_FILE})

#
# Locate the Geant4 package
#
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
include_directories(${Geant4_INCLUDE_DIRS})
link_directories(${Geant4_LIBRARIES})

# -----------------------------------------------------------------------------

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set(PROJECT_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)

set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_BINDIR ${CMAKE_INSTALL_PREFIX}/bin)

include_directories(${PROJECT_INCLUDE_DIR})

add_subdirectory(MCObjs)
add_subdirectory(GLG4Sim)
add_subdirectory(LSCSim)
add_subdirectory(TrgSim)
add_subdirectory(RawObjs)

# install
install(DIRECTORY ${CMAKE_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})


