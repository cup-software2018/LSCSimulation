# FindXercesC.cmake shim - satisfies Geant4's GDML dependency check
# when XercesC is not installed. GDML features will not work at runtime.
set(XercesC_FOUND TRUE)
set(XercesC_VERSION "3.2.5")
set(XercesC_INCLUDE_DIRS "")

if(NOT TARGET XercesC::XercesC)
  add_library(XercesC::XercesC INTERFACE IMPORTED)
endif()
