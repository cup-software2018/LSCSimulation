#include <fstream>
#include <sstream>
#include <string>

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "GLG4param.hh"
#include "LSCDetectorConstruction.hh"
#include "LSCPMTSD.hh"
#include "LSC_PMT_LogicalVolume.hh"

using namespace std;
using namespace CLHEP;

void LSCDetectorConstruction::ConstructDetector_Prototype(G4VPhysicalVolume * worldphys,
                                                          LSCPMTSD * pmtsd, GLG4param & geom_db)
{
  auto WorldLog = worldphys->GetLogicalVolume();

  // Buffer
  G4double bufferR = cm * geom_db["buffer_radius"];
  G4double bufferH = cm * geom_db["buffer_height"];
  G4double bufferT = cm * geom_db["buffer_thickness"];

  auto BufferTankTubs = new G4Tubs("BufferTankTubs", 0, bufferR, bufferH / 2, 0, 360 * deg);
  auto BufferTankLog = new G4LogicalVolume(
      BufferTankTubs, G4Material::GetMaterial("Steel_unpolished"), "BufferTankLog");
  auto visatt = new G4VisAttributes(G4Color(1, 1, 1, 0.5));
  visatt->SetForceSolid(true);
  BufferTankLog->SetVisAttributes(visatt);
  auto BufferTankPhys = new G4PVPlacement(0, G4ThreeVector(), BufferTankLog, "BufferTankPhys",
                                          WorldLog, false, 0, fGeomCheck);

  auto BufferLiquidTubs =
      new G4Tubs("BufferLiquidTubs", 0, bufferR - bufferT, bufferH / 2 - bufferT, 0, 360 * deg);
  auto BufferLiquidLog =
      new G4LogicalVolume(BufferLiquidTubs, G4Material::GetMaterial("Water"), "BufferLog");

  BufferLiquidLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto BufferLiquidPhys = new G4PVPlacement(0, G4ThreeVector(), BufferLiquidLog, "BufferLiquidPhys",
                                            BufferTankLog, false, 0, fGeomCheck);

  new G4LogicalBorderSurface("buffer_logsurf1", BufferTankPhys, BufferLiquidPhys, Stainless_opsurf);
  new G4LogicalBorderSurface("buffer_logsurf2", BufferLiquidPhys, BufferTankPhys, Stainless_opsurf);

  // Target
  G4double targetR = cm * geom_db["target_radius"];
  G4double targetH = cm * geom_db["target_height"];
  G4double targetT = cm * geom_db["target_thickness"];

  auto TargetTankTubs = new G4Tubs("TargetTankTubs", 0, targetR, targetH / 2, 0, 360 * deg);
  auto TargetTankLog =
      new G4LogicalVolume(TargetTankTubs, G4Material::GetMaterial("Acrylic"), "TargetTankLog");
  visatt = new G4VisAttributes();
  visatt->SetForceSolid(true);
  TargetTankLog->SetVisAttributes(visatt);
  auto TargetTankPhys = new G4PVPlacement(0, G4ThreeVector(), TargetTankLog, "TargetTankPhys",
                                          BufferLiquidLog, false, 0, fGeomCheck);

  auto TargetLSTubs =
      new G4Tubs("TargetLSTubs", 0, targetR - targetT, targetH / 2 - targetT, 0, 360 * deg);
  auto TargetLSLog =
      new G4LogicalVolume(TargetLSTubs, G4Material::GetMaterial("LS_LAB"), "TargetLSLog");
  TargetLSLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto TargetLSPhys = new G4PVPlacement(0, G4ThreeVector(), TargetLSLog, "TargetLSPhys",
                                        TargetTankLog, false, 0, fGeomCheck);

  if (fPMTPositionDataFile.empty()) {
    G4String msg = "Error, pmt position data file could not be opened.\n";
    G4cerr << msg << G4endl;
    G4Exception("LSCDetectorConstruction::LSCDetectorConstruction", "", FatalException, msg);
  }

  string line;
  ifstream pmtposfile;

  double coord_x, coord_y, coord_z;
  int pmtno, nring, region;

  ///////////////////////////////////////////////////////////////////////////
  // --- make the fundamental inner  PMT assembly
  ///////////////////////////////////////////////////////////////////////////
  auto _logiInnerPMT = new LSC_10inch_LogicalVolume(
      "InnerPMT", G4Material::GetMaterial("Water"), G4Material::GetMaterial("Glass"),
      Photocathode_opsurf, G4Material::GetMaterial("PMT_Vac"), G4Material::GetMaterial("Steel"),
      nullptr,
      pmtsd); // sensitive detector hook

  char PMTname[64];

  pmtposfile.open(fPMTPositionDataFile.c_str());
  while (getline(pmtposfile, line)) {
    if (line.empty() || line[0] == '#') continue;

    istringstream sline(line);
    sline >> pmtno >> coord_x >> coord_y >> coord_z >> nring >> region;

    sprintf(PMTname, "InnerPMTPhys%d", pmtno);

    G4double r = sqrt(coord_x * coord_x + coord_y * coord_y + coord_z * coord_z);
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

    // shift PMT inward by fPMTOffset from the inner wall
    G4double ox = 0, oy = 0, oz = 0;
    if (region > 0) {
      oz = -fPMTOffset;
    } else if (region < 0) {
      oz = fPMTOffset;
    } else {
      G4double rxy = std::sqrt(coord_x * coord_x + coord_y * coord_y);
      ox = -coord_x / rxy * fPMTOffset;
      oy = -coord_y / rxy * fPMTOffset;
    }
    G4ThreeVector pmtpos(coord_x + ox, coord_y + oy, coord_z + oz);

    new G4PVPlacement(PMT_rotation, pmtpos, PMTname, _logiInnerPMT, BufferLiquidPhys, false, pmtno,
                      fGeomCheck);
  }
}
