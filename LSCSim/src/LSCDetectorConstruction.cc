#include "LSCSim/LSCDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4String.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"
#include "G4UnionSolid.hh"
#include "GLG4Sim/GLG4param.hh"
#include "LSCSim/LSCPMTSD.hh"
#include "LSCSim/LSC_PMT_LogicalVolume.hh"

using namespace std;

LSCPMTSD * LSCDetectorConstruction::fPmt_SD = NULL;

LSCDetectorConstruction::LSCDetectorConstruction()
  : G4VUserDetectorConstruction(),
    G4UImessenger()
{
  fDetectorDir = new G4UIdirectory("/LSC/det/");

  fGeomCheck = 0;

  fGeomCheckOptCmd = new G4UIcmdWithAnInteger("/LSC/det/geomcheck", this);
  fMaterialDataFileCmd = new G4UIcmdWithAString("/LSC/det/materialdata", this);
  fGeometryDataFileCmd = new G4UIcmdWithAString("/LSC/det/geometrydata", this);
  fPMTPositionDataFileCmd = new G4UIcmdWithAString("/LSC/det/pmtposdata", this);
  fWhichDetectorCmd = new G4UIcmdWithAString("/LSC/det/detector", this);

  fGeometryDataFile = "";
  fPMTPositionDataFile = "";
  fMaterialDataFile = "";
  fWhichDetector = "LS";
}

LSCDetectorConstruction::~LSCDetectorConstruction()
{
  delete fDetectorDir;
  delete fGeomCheckOptCmd;
  delete fMaterialDataFileCmd;
  delete fGeometryDataFileCmd;
  delete fPMTPositionDataFileCmd;
}

void LSCDetectorConstruction::SetNewValue(G4UIcommand * command,
                                          G4String newValues)
{
  if (command == fGeomCheckOptCmd) {
    istringstream is(newValues);
    is >> fGeomCheck;
  }
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
  if (command == fWhichDetectorCmd) { fWhichDetector = newValues; }
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
    G4String msg = "Error, geometery data file could not be opened.\n";
    G4cerr << msg << G4endl;
    G4Exception("LSCDetectorConstruction::LSCDetectorConstruction", "",
                FatalException, msg);
  }
  else {
    geom_db.ReadFile(fGeometryDataFile.c_str());
  }

  fWhichDetector = "PROTO";
  G4cout << fWhichDetector << " will be constructed." << G4endl;

  // World (Rock)
  G4double worldX = cm * geom_db["worldx"];
  G4double worldY = cm * geom_db["worldy"];
  G4double worldZ = cm * geom_db["worldz"];
  auto WorldBox = new G4Box("WorldBox", worldX / 2, worldY / 2, worldZ / 2);
  auto WorldLog = new G4LogicalVolume(WorldBox, G4Material::GetMaterial("Rock"),
                                      "WorldLog", 0, 0, 0);
  WorldLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto WorldPhys = new G4PVPlacement(0, G4ThreeVector(), WorldLog, "WorldPhys",
                                     0, false, 0, fGeomCheck);


  // Veto
  G4VPhysicalVolume * VetoLiquidPhys = nullptr;
  if (fWhichDetector !=  "PROTO") {
    G4double vetoR = cm * geom_db["veto_radius"];
    G4double vetoH = cm * geom_db["veto_height"];
    G4double vetoT = cm * geom_db["veto_thickness"];
    auto VetoTankTubs =
      new G4Tubs("VetoTankTubs", 0, vetoR, vetoH / 2, 0, 360 * deg);
    auto VetoTankLog = new G4LogicalVolume(
					   VetoTankTubs, G4Material::GetMaterial("Steel"), "VetoTankLog", 0, 0, 0);
    VetoTankLog->SetVisAttributes(G4Colour::White());
    auto VetoTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), VetoTankLog, "VetoTankPhys",
                        WorldLog, false, 0, fGeomCheck);

    auto VetoLiquidTubs = new G4Tubs("VetoLiquidTubs", 0, vetoR - vetoT,
				     vetoH / 2 - vetoT, 0, 360 * deg);
    auto VetoLiquidLog =
      new G4LogicalVolume(VetoLiquidTubs, G4Material::GetMaterial("Water"),
                          "VetoLiquidLog", 0, 0, 0);
    VetoLiquidLog->SetVisAttributes(G4Colour(0, 0, 1, 0.1));
    auto VetoLiquidPhys =
      new G4PVPlacement(0, G4ThreeVector(), VetoLiquidLog, "VetoLiquidPhys",
                        VetoTankLog, false, 0, fGeomCheck);

    new G4LogicalBorderSurface("veto_logsurf1", VetoTankPhys, VetoLiquidPhys,
			       Stainless_opsurf);
    new G4LogicalBorderSurface("veto_logsurf2", VetoLiquidPhys, VetoTankPhys,
			       Stainless_opsurf);
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // --- PMT sensitive detector
  ///////////////////////////////////////////////////////////////////////////
  G4SDManager * fSDman = G4SDManager::GetSDMpointer();
  LSCPMTSD * pmtSDInner = new LSCPMTSD("/LSC/PMT/inner");
  fSDman->AddNewDetector(pmtSDInner);

  if (fWhichDetector == "LSCC")
    ConstructDetector_LSC_Cylinder(VetoLiquidPhys, pmtSDInner, geom_db);
  else if (fWhichDetector == "LSCS")
    ConstructDetector_LSC_Sphere(VetoLiquidPhys, pmtSDInner, geom_db);    
  else if (fWhichDetector == "PROTO")
    ConstructDetector_Prototype(WorldPhys, pmtSDInner, geom_db);

  return WorldPhys;
}
