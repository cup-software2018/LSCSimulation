#include <fstream>
#include <sstream>
#include <string>

#include "G4Box.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"

#include "GLG4Sim/GLG4param.hh"
#include "LSCSim/LSCDetectorConstruction.hh"
#include "LSCSim/LSCPMTSD.hh"
#include "LSCSim/LSC_PMT_LogicalVolume.hh"

using namespace std;
using namespace CLHEP;

void LSCDetectorConstruction::ConstructDetector_Prototype(
    G4VPhysicalVolume * worldphys, LSCPMTSD * pmtsd, GLG4param & geom_db)
{
  auto WorldLog = worldphys->GetLogicalVolume();

  // Buffer
  G4double bufferR = cm * geom_db["buffer_radius"];
  G4double bufferH = cm * geom_db["buffer_height"];
  G4double bufferT = cm * geom_db["buffer_thickness"];

  auto BufferTankTubs =
      new G4Tubs("BufferTankTubs", 0, bufferR, bufferH / 2, 0, 360 * deg);
  auto BufferTankLog = new G4LogicalVolume(
      BufferTankTubs, G4Material::GetMaterial("Steel_unpolished"), "BufferTankLog");
  auto visatt = new G4VisAttributes(G4Color(1, 1, 1, 0.5));
  visatt->SetForceSolid(true);
  BufferTankLog->SetVisAttributes(visatt);
  auto BufferTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), BufferTankLog, "BufferTankPhys",
                        WorldLog, false, 0, fGeomCheck);

  auto BufferLiquidTubs = new G4Tubs("BufferLiquidTubs", 0, bufferR - bufferT,
                                     bufferH / 2 - bufferT, 0, 360 * deg);
  auto BufferLiquidLog = new G4LogicalVolume(
      BufferLiquidTubs, G4Material::GetMaterial("Water"), "BufferLog");

  BufferLiquidLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto BufferLiquidPhys =
      new G4PVPlacement(0, G4ThreeVector(), BufferLiquidLog, "BufferLiquidPhys",
                        BufferTankLog, false, 0, fGeomCheck);

  new G4LogicalBorderSurface("buffer_logsurf1", BufferTankPhys,
                             BufferLiquidPhys, Stainless_opsurf);
  new G4LogicalBorderSurface("buffer_logsurf2", BufferLiquidPhys,
                             BufferTankPhys, Stainless_opsurf);

  // Target
  G4double targetR = cm * geom_db["target_radius"];
  G4double targetH = cm * geom_db["target_height"];
  G4double targetT = cm * geom_db["target_thickness"];

  auto TargetTankTubs =
      new G4Tubs("TargetTankTubs", 0, targetR, targetH / 2, 0, 360 * deg);
  auto TargetTankLog = new G4LogicalVolume(
      TargetTankTubs, G4Material::GetMaterial("Acrylic"), "TargetTankLog");
  visatt = new G4VisAttributes();
  visatt->SetForceSolid(true);
  TargetTankLog->SetVisAttributes(visatt);
  auto TargetTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), TargetTankLog, "TargetTankPhys",
                        BufferLiquidLog, false, 0, fGeomCheck);

  auto TargetLSTubs = new G4Tubs("TargetLSTubs", 0, targetR - targetT,
                                 targetH / 2 - targetT, 0, 360 * deg);
  auto TargetLSLog = new G4LogicalVolume(
      TargetLSTubs, G4Material::GetMaterial("LS_LAB"), "TargetLSLog");
  TargetLSLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto TargetLSPhys =
      new G4PVPlacement(0, G4ThreeVector(), TargetLSLog, "TargetLSPhys",
                        TargetTankLog, false, 0, fGeomCheck);

  if (fPMTPositionDataFile.empty()) {
    G4String msg = "Error, pmt position data file could not be opened.\n";
    G4cerr << msg << G4endl;
    G4Exception("LSCDetectorConstruction::LSCDetectorConstruction", "",
                FatalException, msg);
  }

  string line;
  ifstream pmtposfile;

  double coord_x, coord_y, coord_z;
  int pmtno, nring, region;

  int doReflector = geom_db["reflector_on"];

  if (doReflector) {
    // reflector
    double reflectorR, reflectorH, PMTR;
    double reflectorT = 10.0 * mm;

    pmtposfile.open(fPMTPositionDataFile.c_str());

    getline(pmtposfile, line);
    istringstream sline(line);
    sline >> reflectorR >> reflectorH >> PMTR;

    reflectorR += reflectorT/2;
    reflectorH += reflectorT;

    // making a hollow cylinder
    G4Tubs * bodyTubs =
        new G4Tubs("bodyTubs", reflectorR, reflectorR + reflectorT,
                   reflectorH / 2.0, 0. * deg, 360. * deg);
    G4Tubs * capTubs = new G4Tubs("capTubs", 0, reflectorR + reflectorT,
                                  reflectorT / 2, 0. * deg, 360. * deg);

    G4UnionSolid * withTop =
        new G4UnionSolid("withTop", bodyTubs, capTubs, 0,
                         G4ThreeVector(0, 0, (reflectorH + reflectorT) / 2.0));

    G4UnionSolid * reflectorTubs =
        new G4UnionSolid("reflectorTubs", withTop, capTubs, 0,
                         G4ThreeVector(0, 0, -(reflectorH + reflectorT) / 2));

    auto holeTubs =
        new G4Tubs("bodyTubs", 0, PMTR, 320 * mm, 0. * deg, 360. * deg);

    // read first pmt position
    getline(pmtposfile, line);
    sline = istringstream(line);
    sline >> pmtno >> coord_x >> coord_y >> coord_z >> nring >> region;

    auto reflectorWithHoles = new G4SubtractionSolid(
        "reflectorWithHoles", reflectorTubs, holeTubs, nullptr,
        G4ThreeVector(coord_x, coord_y, coord_z));

    while (getline(pmtposfile, line)) {
      if (line.empty() || line[0] == '#') continue;
      sline = istringstream(line);
      sline >> pmtno >> coord_x >> coord_y >> coord_z >> nring >> region;

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

      reflectorWithHoles =
          new G4SubtractionSolid("reflectorWithHoles", reflectorWithHoles,
                                 holeTubs, PMT_rotation, pmtpos);
    }
    pmtposfile.close();

    auto reflectorLog = new G4LogicalVolume(
        reflectorWithHoles, G4Material::GetMaterial("Teflon"), "reflectorLog");

    visatt = new G4VisAttributes(G4Color(1, 1, 1, 0.7));
    visatt->SetForceSolid(true);
    reflectorLog->SetVisAttributes(visatt);

    new G4PVPlacement(0, G4ThreeVector(), reflectorLog, "reflectorPhys",
                      BufferLiquidLog, false, 0, fGeomCheck);
  }

  ///////////////////////////////////////////////////////////////////////////
  // --- make the fundamental inner  PMT assembly
  ///////////////////////////////////////////////////////////////////////////
  auto _logiInnerPMT = new LSC_10inch_LogicalVolume(
      "InnerPMT", G4Material::GetMaterial("Water"),
      G4Material::GetMaterial("Glass"), Photocathode_opsurf,
      G4Material::GetMaterial("PMT_Vac"), G4Material::GetMaterial("Steel"),
      nullptr,
      pmtsd); // sensitive detector hook

  char PMTname[64];

  pmtposfile.open(fPMTPositionDataFile.c_str());
  getline(pmtposfile, line);
  while (getline(pmtposfile, line)) {
    if (line.empty() || line[0] == '#') continue;

    istringstream sline(line);
    sline >> pmtno >> coord_x >> coord_y >> coord_z >> nring >> region;

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
                      BufferLiquidPhys, false, pmtno, fGeomCheck);
  }
}
