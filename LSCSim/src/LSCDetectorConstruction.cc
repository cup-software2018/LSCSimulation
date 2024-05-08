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
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"

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

  fGeomCheck = false;

  fGeomCheckOptCmd = new G4UIcmdWithABool("/LSC/det/geomcheck", this);
  fMaterialDataFileCmd = new G4UIcmdWithAString("/LSC/det/materialdata", this);
  fGeometryDataFileCmd = new G4UIcmdWithAString("/LSC/det/geometrydata", this);
  fPMTPositionDataFileCmd = new G4UIcmdWithAString("/LSC/det/pmtposdata", this);

  fGeometryDataFile = "";
  fPMTPositionDataFile = "";
  fMaterialDataFile = "";
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
    istringstream is((const char *)newValues);
    is >> fGeomCheck;
  }
  if (command == fGeometryDataFileCmd) {
    if (fGeometryDataFile.empty()) fGeometryDataFile = newValues;
  }
  if (command == fMaterialDataFileCmd) {
    if (fMaterialDataFile.empty()) fMaterialDataFile = newValues;
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
    G4String msg = "Error, geometery data file could not be opened.\n";
    G4cerr << msg << G4endl;
    G4Exception("LSCDetectorConstruction::LSCDetectorConstruction", "",
                FatalException, msg);
  }
  else {
    geom_db.ReadFile(fGeometryDataFile.c_str());
  }

  // Experimental Hall
  G4double expHall_x = cm * geom_db["worldx"];
  G4double expHall_y = cm * geom_db["worldy"];
  G4double expHall_z = cm * geom_db["worldz"];

  G4Box * ExpHallBox =
      new G4Box("ExpHallBox", expHall_x / 2, expHall_y / 2, expHall_z / 2);

  G4LogicalVolume * ExpHallLog = new G4LogicalVolume(
      ExpHallBox, G4Material::GetMaterial("Air"), "ExpHallLog", 0, 0, 0);

  G4VPhysicalVolume * ExpHallPhys = new G4PVPlacement(
      0, G4ThreeVector(), ExpHallLog, "ExpHallPhys", 0, false, fGeomCheck);

  G4double Test_x = cm * geom_db["testx"];
  G4double Test_y = cm * geom_db["testy"];
  G4double Test_z = cm * geom_db["testz"];

  G4Box * TestBox = new G4Box("TestBox", Test_x / 2, Test_y / 2, Test_z / 2);

  G4LogicalVolume * TestLog = new G4LogicalVolume(
      TestBox, G4Material::GetMaterial("LS_LAB"), "TestLog", 0, 0, 0);

  G4VPhysicalVolume * TestPhys =
      new G4PVPlacement(0, G4ThreeVector(), TestLog, "TargetLSPhys", ExpHallLog,
                        false, fGeomCheck);

  ///////////////////////////////////////////////////////////////////////////
  // --- PMT sensitive detector
  ///////////////////////////////////////////////////////////////////////////
  G4SDManager * fSDman = G4SDManager::GetSDMpointer();
  LSCPMTSD * pmtSDInner = new LSCPMTSD("/LSC/PMT/inner");
  fSDman->AddNewDetector(pmtSDInner);

  // simple 2" pmt with PMT sensitive detector, no using fast simulation
  G4LogicalVolume * apmt = BuildCylindricalPMT(51.8, 200, pmtSDInner);

  new G4PVPlacement(0, G4ThreeVector(), "PMTPhys_0", apmt,
                    TestPhys, // physical parent
                    false, 0, fGeomCheck);
  
  return ExpHallPhys;
}

G4LogicalVolume * LSCDetectorConstruction::BuildCylindricalPMT(
    G4double pmtradius, G4double pmtlength, LSCPMTSD * detector)
{
  //****************** Build PMTs
  G4double ShieldThickness = 1. * mm;

  G4Tubs * PmtTubs = new G4Tubs("PmtTubs", 0.0, pmtradius, pmtlength / 2,
                                0. * deg, 360. * deg);

  pmtradius = pmtradius - ShieldThickness;

  G4LogicalVolume * PmtLog = new G4LogicalVolume(
      PmtTubs, G4Material::GetMaterial("Steel"), "PmtLog", 0, 0, 0);

  G4double pmtglasslength = pmtlength / 2.;
  G4Tubs * PmtGlassTubs = new G4Tubs("PmtGlassTubs", 0.0, pmtradius,
                                     pmtglasslength / 2, 0. * deg, 360. * deg);

  G4LogicalVolume * PmtGlassLog = new G4LogicalVolume(
      PmtGlassTubs, G4Material::GetMaterial("Glass"), "PmtGlassLog", 0, 0, 0);

  // PmtGlassLog->SetVisAttributes(fRedVisSolid);

  G4double pmtglassthickness = 3. * mm;
  G4double pmtvacuumradius = pmtradius - pmtglassthickness;
  G4double pmtvacuumlength = pmtglasslength - 2 * pmtglassthickness;

  G4Tubs * PmtVacuumTubs =
      new G4Tubs("PmtVacuumTubs", 0.0, pmtvacuumradius, pmtvacuumlength / 2,
                 0. * deg, 360. * deg);

  G4LogicalVolume * PmtVacuumLog =
      new G4LogicalVolume(PmtVacuumTubs, G4Material::GetMaterial("PMT_Vac"),
                          "PmtVacuumLog", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), PmtVacuumLog, "PmtVacuumPhys",
                    PmtGlassLog, false, 0, fGeomCheck);

  new G4PVPlacement(0, G4ThreeVector(0, 0, (pmtlength - pmtglasslength) / 2.),
                    PmtGlassLog, "PmtGlassPhys", PmtLog, false, 0, fGeomCheck);

  G4double frontThickness = pmtglassthickness;
  G4double cathodThickness = 10. * mm;
  G4double cathodRadius = pmtvacuumradius;

  G4Tubs * PhotocathodeTubs =
      new G4Tubs("PhotocathodeTubs", 0.0, cathodRadius, cathodThickness / 2,
                 0. * deg, 360. * deg);

  G4LogicalVolume * PhotocathodeLog = new G4LogicalVolume(
      PhotocathodeTubs, G4Material::GetMaterial("photocathode"),
      "PhotocathodeLog");
  auto RedVisSolid = new G4VisAttributes(G4Colour::Red());
  RedVisSolid->SetForceSolid(true);
  PhotocathodeLog->SetVisAttributes(RedVisSolid);
  PhotocathodeLog->SetSensitiveDetector(detector);

  G4VPhysicalVolume * PhotocathodePhys = new G4PVPlacement(
      0, G4ThreeVector(0, 0, (pmtvacuumlength - cathodThickness) / 2.),
      PhotocathodeLog, "PhotocathodePhys", PmtVacuumLog, false, 0, fGeomCheck);

  //**Create logical skin surfaces
  new G4LogicalSkinSurface("PhotocathodLogSkin", PhotocathodeLog,
                           Photocathode_opsurf);

  return PmtLog;
}
