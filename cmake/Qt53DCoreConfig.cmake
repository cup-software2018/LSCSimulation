# Shim for Qt53DCore - satisfies Geant4's dependency check when Qt5 3D is not installed.
# Geant4 built with Qt3D support, but at runtime only uses it if available.
if(NOT TARGET Qt5::3DCore)
  add_library(Qt5::3DCore INTERFACE IMPORTED)
endif()
set(Qt53DCore_FOUND TRUE)
set(Qt53DCore_VERSION 5.15.9)
