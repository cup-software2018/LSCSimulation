#ifndef LSCDetectorConstruction_hh
#define LSCDetectorConstruction_hh

#include "G4UImessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"

class G4Box;
class G4Tubs;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UIdirectory;
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
  void SetMaterialDataFile(const char * file) { fMaterialDataFile = file; }

private:
  void ConstructMaterials();
  G4VPhysicalVolume * ConstructDetector();

  // Sensitive Detectors
  static LSCPMTSD * fPmt_SD;

  // Optical surface
  G4OpticalSurface * Photocathode_opsurf;
  G4OpticalSurface * Stainless_opsurf;
  G4OpticalSurface * Polyethylene_opsurf;
  G4OpticalSurface * Tyvek_opsurf;
  G4OpticalSurface * Teflon_opsurf;

  G4String fGeometryDataFile;
  G4String fMaterialDataFile;

  G4UIdirectory * fDetectorDir;
  G4UIcmdWithAString * fGeomCheckOptCmd;
  G4UIcmdWithAString * fSourcePositionCmd;
  G4UIcmdWithAString * fCalibrationOptCmd;
  G4int fGeomCheck;
};

#endif
