#include "LSCSim/LSCSteppingAction.hh"
#include "LSCSim/LSCRecorderBase.hh"

#include "G4VProcess.hh"
#include "G4SteppingManager.hh"

LSCSteppingAction::LSCSteppingAction(LSCRecorderBase * r)
    : fRecorder(r)
{
}

LSCSteppingAction::~LSCSteppingAction() {}

void LSCSteppingAction::UserSteppingAction(const G4Step * theStep)
{
  G4SteppingManager * steppingManager = fpSteppingManager;
  G4VProcess * proc = steppingManager->GetfCurrentProcess();

  if (fRecorder) fRecorder->RecordStep(theStep, proc);
}
