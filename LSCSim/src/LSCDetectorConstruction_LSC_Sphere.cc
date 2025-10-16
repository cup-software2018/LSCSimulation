#include <fstream>
#include <sstream>
#include <string>

#include "G4Box.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"

#include "GLG4Sim/GLG4param.hh"
#include "LSCSim/LSCDetectorConstruction.hh"
#include "LSCSim/LSCPMTSD.hh"
#include "LSCSim/LSC_PMT_LogicalVolume.hh"

using namespace std;
using namespace CLHEP;

void LSCDetectorConstruction::ConstructDetector_LSC_Sphere(
    G4VPhysicalVolume * vetoliquid, LSCPMTSD * pmtsd, GLG4param & geom_db)
{
  auto VetoLiquidLog = vetoliquid->GetLogicalVolume();

  // Buffer
  G4double bufferR = cm * geom_db["buffer_radius"];
  G4double bufferT = cm * geom_db["buffer_thickness"];
  auto BufferTankSphere = new G4Sphere("BufferTankSphere", 0, bufferR, 0,
                                       360. * deg, 0, 180. * deg);
  auto BufferTankLog =
      new G4LogicalVolume(BufferTankSphere, G4Material::GetMaterial("Steel"),
                          "BufferTankLog", 0, 0, 0);
  // BufferTankLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  BufferTankLog->SetVisAttributes(G4Colour::Black());
  auto BufferTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), BufferTankLog, "BufferTankPhys",
                        VetoLiquidLog, false, fGeomCheck);

  auto BufferLiquidSphere = new G4Sphere(
      "BufferLiquidSphere", 0, bufferR - bufferT, 0, 360. * deg, 0, 180. * deg);
  auto BufferLiquidLog =
      new G4LogicalVolume(BufferLiquidSphere, G4Material::GetMaterial("Water"),
                          "BufferLog", 0, 0, 0);
  // BufferLiquidLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  BufferLiquidLog->SetVisAttributes(G4Colour(0, 0, 1, 0.1)); // blue
  auto BufferLiquidPhys =
      new G4PVPlacement(0, G4ThreeVector(), BufferLiquidLog, "BufferLiquidPhys",
                        BufferTankLog, false, fGeomCheck);

  new G4LogicalBorderSurface("buffer_logsurf1", BufferTankPhys,
                             BufferLiquidPhys, Stainless_opsurf);
  new G4LogicalBorderSurface("buffer_logsurf2", BufferLiquidPhys,
                             BufferTankPhys, Stainless_opsurf);

  // Target
  G4double targetR = cm * geom_db["target_radius"];
  G4double targetT = cm * geom_db["target_thickness"];

  auto TargetTankSphere = new G4Sphere("TargetTankSphere", 0, targetR, 0,
                                       360. * deg, 0, 180. * deg);
  auto TargetTankLog =
      new G4LogicalVolume(TargetTankSphere, G4Material::GetMaterial("Acrylic"),
                          "TargetTankLog", 0, 0, 0);
  // TargetTankLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  TargetTankLog->SetVisAttributes(G4Colour(1, 1, 1, 0.1)); // white
  auto TargetTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), TargetTankLog, "TargetTankPhys",
                        BufferLiquidLog, false, fGeomCheck);

  auto TargetLSSphere = new G4Sphere("TargetLSSphere", 0, targetR - targetT, 0,
                                     360. * deg, 0, 180. * deg);
  auto TargetLSLog =
      new G4LogicalVolume(TargetLSSphere, G4Material::GetMaterial("Pyrene_LS"),
                          "TargetLSLog", 0, 0, 0);
  // TargetLSLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  TargetLSLog->SetVisAttributes(G4Colour(0, 1, 0, 0.1)); // green
  auto TargetLSPhys =
      new G4PVPlacement(0, G4ThreeVector(), TargetLSLog, "TargetLSPhys",
                        TargetTankLog, false, fGeomCheck);

  ///////////////////////////////////////////////////////////////////////////
  // --- make the fundamental inner  PMT assembly
  ///////////////////////////////////////////////////////////////////////////
  auto _logiInnerPMT20 = new LSC_20inch_LogicalVolume(
      "InnerPMT", G4Material::GetMaterial("Water"),
      G4Material::GetMaterial("Glass"), Photocathode_opsurf,
      G4Material::GetMaterial("PMT_Vac"), G4Material::GetMaterial("Steel"),
      nullptr,
      pmtsd); // sensitive detector hook

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

    sprintf(PMTname, "InnerPMTPhys%d", pmtno);

    G4double r =
        sqrt(coord_x * coord_x + coord_y * coord_y + coord_z * coord_z);
    G4double dx = -coord_x / r;
    G4double dy = -coord_y / r;
    G4double dz = -coord_z / r;

    double angle_y = (-1.0) * atan2(dx, dz);
    double angle_x = atan2(dy, sqrt(dx * dx + dz * dz));

    auto PMT_rotation = new G4RotationMatrix();
    PMT_rotation->rotateY(angle_y);
    PMT_rotation->rotateX(angle_x);

    G4ThreeVector pmtpos(coord_x, coord_y, coord_z);

    new G4PVPlacement(PMT_rotation, pmtpos, PMTname, _logiInnerPMT20,
                      BufferLiquidPhys, false, pmtno - 1, fGeomCheck);
  }
}
