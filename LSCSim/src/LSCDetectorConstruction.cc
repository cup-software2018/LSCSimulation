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
#include "G4UIcmdWithAnInteger.hh"
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

  fGeomCheck = 0;

  fGeomCheckOptCmd = new G4UIcmdWithAnInteger("/LSC/det/geomcheck", this);
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

  // World (Rock)
  G4double worldX = cm * geom_db["worldx"];
  G4double worldY = cm * geom_db["worldy"];
  G4double worldZ = cm * geom_db["worldz"];

  auto WorldBox = new G4Box("WorldBox", worldX / 2, worldY / 2, worldZ / 2);
  auto WorldLog = new G4LogicalVolume(WorldBox, G4Material::GetMaterial("Rock"),
                                      "WorldLog", 0, 0, 0);
  WorldLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto WorldPhys = new G4PVPlacement(0, G4ThreeVector(), WorldLog, "WorldPhys",
                                     0, false, fGeomCheck);

  // Buffer
  G4double bufferR = cm * geom_db["buffer_radius"];
  G4double bufferH = cm * geom_db["buffer_height"];
  G4double bufferT = cm * geom_db["buffer_tichkness"];

  auto BufferTankTubs =
      new G4Tubs("BufferTankTubs", 0, bufferR, bufferH / 2, 0, 360 * deg);
  auto BufferTankLog =
      new G4LogicalVolume(BufferTankTubs, G4Material::GetMaterial("Steel"),
                          "BufferTankLog", 0, 0, 0);
  BufferTankLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto BufferTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), BufferTankLog, "BufferTankPhys",
                        WorldLog, false, fGeomCheck);

  auto BufferLiquidTubs = new G4Tubs("BufferLiquidTubs", 0, bufferR - bufferT,
                               bufferH / 2 - bufferT, 0, 360 * deg);
  auto BufferLiquidLog = new G4LogicalVolume(
      BufferLiquidTubs, G4Material::GetMaterial("Water"), "BufferLog", 0, 0, 0);
  BufferLiquidLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto BufferLiquidPhys =
      new G4PVPlacement(0, G4ThreeVector(), BufferLiquidLog, "BufferLiquidPhys",
                        BufferTankLog, false, fGeomCheck);

  new G4LogicalBorderSurface("buffer_logsurf1", BufferTankPhys, BufferLiquidPhys,
                             Stainless_opsurf);
  new G4LogicalBorderSurface("buffer_logsurf2", BufferLiquidPhys, BufferTankPhys,
                             Stainless_opsurf);

  // Target
  G4double targetR = cm * geom_db["target_radius"];
  G4double targetH = cm * geom_db["target_height"];
  G4double targetT = cm * geom_db["target_tichkness"];

  auto TargetTankTubs =
      new G4Tubs("TargetTankTubs", 0, targetR, targetH / 2, 0, 360 * deg);
  auto TargetTankLog =
      new G4LogicalVolume(TargetTankTubs, G4Material::GetMaterial("Acrylic"),
                          "TargetTankLog", 0, 0, 0);
  // TargetTankLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto TargetTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), TargetTankLog, "TargetTankPhys",
                        BufferLiquidLog, false, fGeomCheck);

  auto TargetLSTubs = new G4Tubs("TargetLSTubs", 0, targetR - targetT,
                               targetH / 2 - targetT, 0, 360 * deg);
  auto TargetLSLog = new G4LogicalVolume(
      TargetLSTubs, G4Material::GetMaterial("LS_LAB"), "TargetLSLog", 0, 0, 0);
  TargetLSLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto TargetLSPhys =
      new G4PVPlacement(0, G4ThreeVector(), TargetLSLog, "TargetLSPhys",
                        TargetTankLog, false, fGeomCheck);

  ///////////////////////////////////////////////////////////////////////////
  // --- PMT sensitive detector
  ///////////////////////////////////////////////////////////////////////////
  G4SDManager * fSDman = G4SDManager::GetSDMpointer();
  LSCPMTSD * pmtSDInner = new LSCPMTSD("/LSC/PMT/inner");
  fSDman->AddNewDetector(pmtSDInner);

  ///////////////////////////////////////////////////////////////////////////
  // --- make the fundamental inner  PMT assembly
  ///////////////////////////////////////////////////////////////////////////
  auto _logiInnerPMT = new LSC_10inch_LogicalVolume(
      "InnerPMT", G4Material::GetMaterial("Water"),
      G4Material::GetMaterial("Glass"), Photocathode_opsurf,
      G4Material::GetMaterial("PMT_Vac"), G4Material::GetMaterial("Steel"),
      nullptr,
      pmtSDInner); // sensitive detector hook

  if (fPMTPositionDataFile.empty()) {
    G4String msg = "Error, pmt position data file could not be opened.\n";
    G4cerr << msg << G4endl;
    G4Exception("LSCDetectorConstruction::LSCDetectorConstruction", "",
                FatalException, msg);
  }

  char PMTname[64];

  double coord_x, coord_y, coord_z;
  int pmtno, nring, region;

  string line;
  ifstream pmtposfile(fPMTPositionDataFile.c_str());
  while (getline(pmtposfile, line)) {
    if (line.empty() || line[0] == '#') continue;

    istringstream sline(line);
    sline >> pmtno >> coord_x >> coord_y >> coord_z >> nring >> region;

    // coord_x *= cm;
    // coord_y *= cm;
    // coord_z *= cm;

    sprintf(PMTname, "InnerPMTPhys%d", pmtno);

    G4double r =
        sqrt(coord_x * coord_x + coord_y * coord_y + coord_z * coord_z);
    G4double dx = -coord_x / r;
    G4double dy = -coord_y / r;
    G4double dz = -coord_z / r;

    double angle_z = atan2(dx, dy);
    double angle_x = atan2(dz, sqrt(dx * dx + dy * dy));

    if (region != 0) {
      // top or bottom PMTs
      double normal_angle = (region > 0 ? -M_PI / 2 : M_PI / 2);
      angle_x = normal_angle;
    }
    else {
      angle_x = 0;
    }

    auto PMT_rotation = new G4RotationMatrix();
    PMT_rotation->rotateZ(angle_z);
    PMT_rotation->rotateX(M_PI / 2.0 - angle_x);

    G4ThreeVector pmtpos(coord_x, coord_y, coord_z);

    new G4PVPlacement(PMT_rotation, pmtpos, PMTname, _logiInnerPMT,
                      BufferLiquidPhys, false, pmtno - 1, fGeomCheck);
  }

  return WorldPhys;
}
