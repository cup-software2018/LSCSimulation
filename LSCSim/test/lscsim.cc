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
       << "              [-g geometry] [-m material] [-v]" << endl;
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

  while ((opt = getopt(argc, argv, "o:f:m:g:n:vh")) != -1) {
    switch (opt) {
      case 'o': outputFileName = G4String(optarg); break;
      case 'f': macroFileName = G4String(optarg); break;
      case 'm': materialData = G4String(optarg); break;
      case 'g': geometryData = G4String(optarg); break;
      case 'n': nevent = atoi(G4String(optarg).data()); break;
      case 'v': doVis = 1; break;
      case 'h': PrintHelp(); break;
      default: PrintHelp();
    }
  }

  if (!doVis && nevent == 0) {
    G4cout << "LSCSim main: Number of event is not set !" << G4endl;
    PrintHelp();
  }
  if (macroFileName.empty()) {
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
  if (!geometryData.empty()) LSCDetector->SetGeometryDataFile(geometryData);
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

  G4VisManager * visManager = NULL;
  if (doVis) {
    visManager = new G4VisExecutive;
    visManager->Initialize();
  }

  // get the pointer to the UI manager and set verbosities
  G4UImanager * UImanager = G4UImanager::GetUIpointer();

  G4String command = "/run/initialize";
  UImanager->ApplyCommand(command);

  command = "/control/execute ";
  command += macroFileName;
  UImanager->ApplyCommand(command);

  if (!doVis) {
    command = Form("/run/beamOn %d", nevent);
    UImanager->ApplyCommand(command);
  }
  else {
    UImanager->ApplyCommand("/control/execute vis.mac");

    G4UIsession * theSession = new G4UIterminal(new G4UItcsh);
    theSession->SessionStart();
    delete theSession;
  }

  if (rootManager) delete rootManager;
  if (visManager != NULL) delete visManager;

  // job termination
  delete runManager;
  return 0;
}
