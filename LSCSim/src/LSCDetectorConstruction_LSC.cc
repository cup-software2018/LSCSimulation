#include <fstream>
#include <sstream>
#include <string>

#include "G4Box.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "GLG4Sim/GLG4param.hh"
#include "LSCSim/LSCDetectorConstruction.hh"
#include "LSCSim/LSCPMTSD.hh"
#include "LSCSim/LSC_PMT_LogicalVolume.hh"

using namespace std;
using namespace CLHEP;

void LSCDetectorConstruction::ConstructDetector_LSC(
    G4VPhysicalVolume * worldphys, LSCPMTSD * pmtsd, GLG4param & geom_db)
{
  auto WorldLog = worldphys->GetLogicalVolume();

  /*
    // Pit
    G4double pitR = cm * geom_db["pit_radius"];
    G4double pitH = cm * geom_db["pit_height"];
    auto PitTubs = new G4Tubs("PitTubs", 0, pitR, pitH / 2, 0, 360 * deg);

    // Tunnel
    G4double tunnelR = cm * geom_db["tunnel_radius"];
    G4double tunnelH = cm * geom_db["tunnel_height"];
    auto TunnelTubs =
        new G4Tubs("TunnelTubs", 0, tunnelR, tunnelH / 2, 0 * deg, 180 * deg);

    auto TunnelRot = new G4RotationMatrix();
    TunnelRot->rotateX(-90 * deg);
    auto PitUni = new G4UnionSolid(
        "PitUni_1", PitTubs, TunnelTubs, TunnelRot,
        G4ThreeVector(0, TunnelTubs->GetZHalfLength() * 1.8, pitH / 2));
    TunnelRot->rotateX(180 * deg);
    TunnelRot->rotateZ(90 * deg);
    PitUni = new G4UnionSolid(
        "PitUni_2", PitUni, TunnelTubs,
        G4Transform3D(*TunnelRot,
                      G4ThreeVector(TunnelTubs->GetZHalfLength() * 1.8, 0,
                                    pitH / 2 - tunnelR)));
    TunnelRot->rotateY(-10 * deg);
    TunnelRot->rotateZ(-135 * deg);
    PitUni = new G4UnionSolid(
        "PitUni_3", PitUni, TunnelTubs,
        G4Transform3D(*TunnelRot, G4ThreeVector(-pitR * 1.25, -pitR * 1.25,
                                                -pitH / 2 + tunnelR / 2)));
    // Cavern (Top)
    G4double cavernR = cm * geom_db["cavern_arch_radius"];
    G4double cavernH = cm * geom_db["cavern_arch_height"];
    G4double cavernT = cm * geom_db["cavern_arch_thickness"];
    G4double cavernA = deg * geom_db["cavern_arch_angle"];
    auto CavernTubs = new G4Tubs("CavernTubs", cavernR - cavernH / 8, cavernR,
                                 cavernH / 2, 0 * deg, cavernA * 2);

    // Cavern (Bottom)
    G4double cavernX = cm * geom_db["cavern_arch_height"];
    G4double cavernY = cm * geom_db["cavern_arch_height"];
    G4double cavernZ = cm * geom_db["cavern_height"];
    auto CavernBox = new G4Box("CavernBox", cavernX / 2, cavernY / 2,
                               (cavernZ - cavernT) / 2.);
    auto CavernRot = new G4RotationMatrix();
    CavernRot->rotateZ(-cavernA);
    CavernRot->rotateY(-90 * deg);

    // Cavern (Union of pit, tunnel, top and bottom caverns)
    auto CavernUni = new G4UnionSolid(
        "CavernUni", CavernBox, CavernTubs,
        G4Transform3D(*CavernRot, G4ThreeVector(0, 0,
                                                CavernBox->GetZHalfLength() -
                                                    (cavernR - cavernH / 10))));
    CavernUni = new G4UnionSolid(
        "CavernUni", PitUni, CavernUni, 0,
        G4ThreeVector(0, 0, pitH / 2. + CavernBox->GetZHalfLength()));
    auto CavernLog = new G4LogicalVolume(
        CavernUni, G4Material::GetMaterial("Air"), "CavernLog");
    CavernLog->SetVisAttributes(G4Colour::Gray());
    auto CavernPhys =
        new G4PVPlacement(0, G4ThreeVector(), CavernLog, "CavernPhys", WorldLog,
                          false, 0, fGeomCheck);
*/
  // Veto
  G4double vetoR = cm * geom_db["veto_radius"];
  G4double vetoH = cm * geom_db["veto_height"];
  G4double vetoT = cm * geom_db["veto_thickness"];
  auto VetoTankTubs =
      new G4Tubs("VetoTankTubs", 0, vetoR, vetoH / 2, 0, 360 * deg);
  auto VetoTankLog = new G4LogicalVolume(
      VetoTankTubs, G4Material::GetMaterial("Steel"), "VetoTankLog", 0, 0, 0);
  VetoTankLog->SetVisAttributes(G4Colour::Black());
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

  // Buffer
  G4double bufferR = cm * geom_db["buffer_radius"];
  G4double bufferH = cm * geom_db["buffer_height"];
  G4double bufferT = cm * geom_db["buffer_thickness"];
  auto BufferTankTubs =
      new G4Tubs("BufferTankTubs", 0, bufferR, bufferH / 2, 0, 360 * deg);
  auto BufferTankLog =
      new G4LogicalVolume(BufferTankTubs, G4Material::GetMaterial("Steel"),
                          "BufferTankLog", 0, 0, 0);
  // BufferTankLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  BufferTankLog->SetVisAttributes(G4Colour::Black());
  auto BufferTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), BufferTankLog, "BufferTankPhys",
                        WorldLog, false, fGeomCheck);

  auto BufferLiquidTubs = new G4Tubs("BufferLiquidTubs", 0, bufferR - bufferT,
                                     bufferH / 2 - bufferT, 0, 360 * deg);
  auto BufferLiquidLog = new G4LogicalVolume(
      BufferLiquidTubs, G4Material::GetMaterial("Water"), "BufferLog", 0, 0, 0);
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
  G4double targetH = cm * geom_db["target_height"];
  G4double targetT = cm * geom_db["target_thickness"];

  auto TargetTankTubs =
      new G4Tubs("TargetTankTubs", 0, targetR, targetH / 2, 0, 360 * deg);
  auto TargetTankLog =
      new G4LogicalVolume(TargetTankTubs, G4Material::GetMaterial("Acrylic"),
                          "TargetTankLog", 0, 0, 0);
  // TargetTankLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  TargetTankLog->SetVisAttributes(G4Colour(1, 1, 1, 0.1)); // white
  auto TargetTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), TargetTankLog, "TargetTankPhys",
                        BufferLiquidLog, false, fGeomCheck);

  auto TargetLSTubs = new G4Tubs("TargetLSTubs", 0, targetR - targetT,
                                 targetH / 2 - targetT, 0, 360 * deg);
  auto TargetLSLog = new G4LogicalVolume(
      TargetLSTubs, G4Material::GetMaterial("LS_LAB"), "TargetLSLog", 0, 0, 0);
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

    new G4PVPlacement(PMT_rotation, pmtpos, PMTname, _logiInnerPMT20,
                      BufferLiquidPhys, false, pmtno - 1, fGeomCheck);
  }
}