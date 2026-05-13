#include "GLG4Sim/GLG4param.hh"
#include "LSCSim/LSCDetectorConstruction.hh"
#include "LSCSim/LSCPMTSD.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "G4VisAttributes.hh"

LSCPMTSD * LSCDetectorConstruction::fPmt_SD = nullptr;

LSCDetectorConstruction::LSCDetectorConstruction()
  : G4VUserDetectorConstruction(),
    G4UImessenger(),
    fGeomCheck(0)
{
  fDetectorDir = new G4UIdirectory("/LSC/det/");

  fGeomCheckOptCmd = new G4UIcmdWithAnInteger("/LSC/det/geomcheck", this);
  fMaterialDataFileCmd = new G4UIcmdWithAString("/LSC/det/materialdata", this);
  fGeometryDataFileCmd = new G4UIcmdWithAString("/LSC/det/geometrydata", this);
  fPMTPositionDataFileCmd = new G4UIcmdWithAString("/LSC/det/pmtposdata", this);
}

LSCDetectorConstruction::~LSCDetectorConstruction()
{
  delete fDetectorDir;
  delete fGeomCheckOptCmd;
  delete fMaterialDataFileCmd;
  delete fGeometryDataFileCmd;
  delete fPMTPositionDataFileCmd;
}

void LSCDetectorConstruction::SetNewValue(G4UIcommand * command, G4String newValues)
{
  if (command == fGeomCheckOptCmd) { fGeomCheck = fGeomCheckOptCmd->GetNewIntValue(newValues); }
  if (command == fGeometryDataFileCmd) {
    if (fGeometryDataFile.empty()) fGeometryDataFile = newValues;
  }
  if (command == fMaterialDataFileCmd) {
    if (fMaterialDataFile.empty()) fMaterialDataFile = newValues;
    G4cout << fMaterialDataFile << G4endl;
  }
  if (command == fPMTPositionDataFileCmd) {
    if (fPMTPositionDataFile.empty()) fPMTPositionDataFile = newValues;
  }
}

G4VPhysicalVolume * LSCDetectorConstruction::Construct()
{
  ConstructMaterials();
  return ConstructDetector();
}

G4VPhysicalVolume * LSCDetectorConstruction::ConstructDetector()
{
  GLG4param & geom_db(GLG4param::GetDB());

  if (fGeometryDataFile.empty()) {
    G4String msg = "Error, geometry data file could not be opened.\n";
    G4cerr << msg << G4endl;
    G4Exception("LSCDetectorConstruction::ConstructDetector", "", FatalException, msg);
  }
  geom_db.ReadFile(fGeometryDataFile.c_str());

  // World (Rock)
  G4double worldX = cm * geom_db["worldx"];
  G4double worldY = cm * geom_db["worldy"];
  G4double worldZ = cm * geom_db["worldz"];
  auto WorldBox = new G4Box("WorldBox", worldX / 2, worldY / 2, worldZ / 2);
  auto WorldLog = new G4LogicalVolume(WorldBox, G4Material::GetMaterial("Rock"), "WorldLog");
  WorldLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto WorldPhys =
      new G4PVPlacement(0, G4ThreeVector(), WorldLog, "WorldPhys", 0, false, 0, fGeomCheck);

  G4SDManager * fSDman = G4SDManager::GetSDMpointer();
  LSCPMTSD * pmtSDInner = new LSCPMTSD("/LSC/PMT/inner");
  fSDman->AddNewDetector(pmtSDInner);

  ConstructDetector_Prototype(WorldPhys, pmtSDInner, geom_db);

  return WorldPhys;
}
