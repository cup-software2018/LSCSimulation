#ifndef LSCDetectorConstruction_hh
#define LSCDetectorConstruction_hh

#include "G4OpticalSurface.hh"
#include "G4UImessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"

class G4VPhysicalVolume;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class LSCPMTSD;

class LSCDetectorConstruction : public G4VUserDetectorConstruction,
                                public G4UImessenger {
public:
  LSCDetectorConstruction();
  virtual ~LSCDetectorConstruction();

  virtual G4VPhysicalVolume * Construct();
  virtual void SetNewValue(G4UIcommand *, G4String);

  void SetGeometryDataFile(const char * file) { fGeometryDataFile = file; }
  void SetPMTPositionDataFile(const char * file)
  {
    fPMTPositionDataFile = file;
  }
  void SetMaterialDataFile(const char * file) { fMaterialDataFile = file; }

private:
  // Sensitive Detectors
  static LSCPMTSD * fPmt_SD;

  void ConstructMaterials();
  G4VPhysicalVolume * ConstructDetector();
  G4LogicalVolume * BuildCylindricalPMT(G4double pmtradius, G4double pmtlength,
                                        LSCPMTSD * detector);

  // Optical surface
  G4OpticalSurface * Photocathode_opsurf;
  G4OpticalSurface * Stainless_opsurf;
  G4OpticalSurface * Polyethylene_opsurf;
  G4OpticalSurface * Tyvek_opsurf;
  G4OpticalSurface * Teflon_opsurf;

  G4String fMaterialDataFile;
  G4String fGeometryDataFile;
  G4String fPMTPositionDataFile;

  G4UIdirectory * fDetectorDir;
  G4UIcmdWithABool * fGeomCheckOptCmd;
  G4UIcmdWithAString * fMaterialDataFileCmd;
  G4UIcmdWithAString * fGeometryDataFileCmd;
  G4UIcmdWithAString * fPMTPositionDataFileCmd;
  G4int fGeomCheck;
};

#endif
