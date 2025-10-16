#include <iostream>
#include <stdlib.h>

#include "G4RunManager.hh"
#include "G4String.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"
#include "TRandom3.h"

#include "CLHEP/Random/MTwistEngine.h"
#include "GLG4Sim/GLG4PrimaryGeneratorAction.hh"
#include "LSCSim/LSCDetectorConstruction.hh"
#include "LSCSim/LSCEventAction.hh"
#include "LSCSim/LSCPhysicsList.hh"
#include "LSCSim/LSCRootManager.hh"
#include "LSCSim/LSCRunAction.hh"
#include "LSCSim/LSCSteppingAction.hh"
#include "LSCSim/LSCTrackingAction.hh"

using namespace std;

G4int NumGamma;
G4double GammaEnergies[15];

void PrintHelp()
{
  cout << endl;
  cout << "Usage: LSCSim [-n # of event] [-o output] [-f macro]" << endl
       << "              [-g geometry] [-p pmt_position data] [-m material] "
          "[-v macro]"
       << endl;
  cout << endl;

  exit(0);
}

int main(int argc, char ** argv)
{
  if (argc < 2) {
    PrintHelp();
    return 0;
  }

  int opt;
  int nevent = 1;
  int doVis = 0;

  G4String outputFileName;
  G4String macroFileName;
  G4String materialData;
  G4String geometryData;
  G4String pmtposData;
  G4String detectorType = "LSCS";

  while ((opt = getopt(argc, argv, "o:f:m:g:p:n:vhcst")) != -1) {
    switch (opt) {
      case 'o': outputFileName = G4String(optarg); break;
      case 'f': macroFileName = G4String(optarg); break;
      case 'm': materialData = G4String(optarg); break;
      case 'g': geometryData = G4String(optarg); break;
      case 'p': pmtposData = G4String(optarg); break;
      case 'n': nevent = atoi(G4String(optarg).data()); break;
      case 'v': doVis = 1; break;
      case 'c': detectorType = "LSCC"; break;
      case 's': detectorType = "LSCS"; break;
      case 't': detectorType = "PROTO"; break;
      case 'h': PrintHelp(); break;
      default: PrintHelp();
    }
  }

  if (!doVis && nevent == 0) {
    G4cout << "LSCSim main: Number of event is not set !" << G4endl;
    PrintHelp();
  }
  if (!doVis && macroFileName.empty()) {
    G4cout << "LSCSim main: No macro file !" << G4endl;
    PrintHelp();
  }

  gRandom->SetSeed(0);
  int seed = gRandom->Integer(kMaxInt);

  CLHEP::MTwistEngine randomEngine;
  G4Random::setTheEngine(&randomEngine);
  G4Random::getTheEngine()->setSeed(seed, 0);

  G4RunManager * runManager = new G4RunManager;

  LSCDetectorConstruction * LSCDetector = new LSCDetectorConstruction();
  LSCDetector->SetDetectorType(detectorType);
  if (!geometryData.empty()) LSCDetector->SetGeometryDataFile(geometryData);
  if (!pmtposData.empty()) LSCDetector->SetPMTPositionDataFile(pmtposData);
  if (!materialData.empty()) LSCDetector->SetMaterialDataFile(materialData);

  runManager->SetUserInitialization(LSCDetector);
  runManager->SetUserInitialization(new LSCPhysicsList);

  if (outputFileName.empty()) outputFileName = "simout.root";

  LSCRootManager * rootManager = new LSCRootManager();
  rootManager->SetRootFile(outputFileName.data());

  runManager->SetUserAction(new GLG4PrimaryGeneratorAction);
  runManager->SetUserAction(new LSCRunAction(rootManager));
  runManager->SetUserAction(new LSCEventAction(rootManager));
  runManager->SetUserAction(new LSCTrackingAction(rootManager));
  runManager->SetUserAction(new LSCSteppingAction(rootManager));

  G4UImanager * UImanager = G4UImanager::GetUIpointer();

  G4UIExecutive * ui = nullptr;
  G4VisManager * visManager = nullptr;
  if (doVis) {
    ui = new G4UIExecutive(argc, argv);
    visManager = new G4VisExecutive;
    visManager->Initialize();

    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  else {
    G4String command = Form("/control/execute %s", macroFileName.c_str());
    UImanager->ApplyCommand(command);
    command = Form("/run/beamOn %d", nevent);
    UImanager->ApplyCommand(command);
  }
  
  if (rootManager) delete rootManager;
  if (visManager != nullptr) delete visManager;

  // job termination
  delete runManager;
  return 0;
}
