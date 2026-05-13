# Shim for Qt53DExtras - satisfies Geant4's dependency check when Qt5 3D is not installed.
if(NOT TARGET Qt5::3DExtras)
  add_library(Qt5::3DExtras INTERFACE IMPORTED)
endif()
set(Qt53DExtras_FOUND TRUE)
set(Qt53DExtras_VERSION 5.15.9)
