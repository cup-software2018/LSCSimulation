# ------------------------------------------------------------
# LSCSimulation top-level CMakeLists.txt
# ------------------------------------------------------------
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(LSCSimulation VERSION 1.0.0 LANGUAGES CXX)

# ---- Compilation Database ----
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---- C++ Standard ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Use standard GNU installation directories
include(GNUInstallDirs)

# ---- Relocatable RPATH Settings (Linux) ----
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")

# ---- Global Compiler Settings ----
add_compile_options(-O2)

# Allow unresolved symbols from Geant4's optional components (Qt3D, XercesC/GDML)
# that may not be installed on this system but are linked into Geant4 shared libs.
add_link_options(-Wl,--allow-shlib-undefined)

# ---- ROOT ----
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
include(${ROOT_USE_FILE})

# ---- Geant4 ----
# cmake/ contains shims for Qt5 3D and XercesC that may be missing on this system
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# Fall back to $G4BASE env var if Geant4_DIR is not set or was cached as NOTFOUND
if((NOT Geant4_DIR OR Geant4_DIR MATCHES "NOTFOUND") AND DEFINED ENV{G4BASE})
  set(Geant4_DIR "$ENV{G4BASE}/lib64/cmake/Geant4" CACHE PATH "Geant4 cmake config dir" FORCE)
endif()
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# ---- Submodules ----
add_subdirectory(MCObjs)
add_subdirectory(GLG4Sim)
add_subdirectory(LSCSim)
add_subdirectory(TrgSim)

# ---- Environment Setup Script ----
file(WRITE ${CMAKE_BINARY_DIR}/setup_lscsim.sh
  "#!/bin/bash
# LSCSimulation Environment Setup Script
# Generated automatically by CMake.

export LSCSIM_DIR=\"\$(builtin cd \"\$(dirname \"\${BASH_SOURCE[0]}\")\" && pwd)\"

export PATH=\"\$LSCSIM_DIR/bin:\$PATH\"
export LD_LIBRARY_PATH=\"\$LSCSIM_DIR/${CMAKE_INSTALL_LIBDIR}:\$LD_LIBRARY_PATH\"
export ROOT_INCLUDE_PATH=\"\$LSCSIM_DIR/include:\$ROOT_INCLUDE_PATH\"

echo \"[LSCSimulation] Environment loaded from: \$LSCSIM_DIR\"
")

install(FILES ${CMAKE_BINARY_DIR}/setup_lscsim.sh DESTINATION .)
