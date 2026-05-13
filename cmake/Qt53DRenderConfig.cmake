# Shim for Qt53DRender - satisfies Geant4's dependency check when Qt5 3D is not installed.
if(NOT TARGET Qt5::3DRender)
  add_library(Qt5::3DRender INTERFACE IMPORTED)
endif()
set(Qt53DRender_FOUND TRUE)
set(Qt53DRender_VERSION 5.15.9)
