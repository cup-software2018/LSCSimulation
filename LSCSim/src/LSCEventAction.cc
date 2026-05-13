#include "G4Event.hh"
#include "G4EventManager.hh"
#include "LSCEventAction.hh"
#include "LSCRecorderBase.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LSCEventAction::LSCEventAction(LSCRecorderBase * r)
  : fRecorder(r)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LSCEventAction::~LSCEventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LSCEventAction::BeginOfEventAction(const G4Event * anEvent)
{
  if (fRecorder) fRecorder->BeginOfEvent(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LSCEventAction::EndOfEventAction(const G4Event * anEvent)
{
  if (fRecorder) fRecorder->EndOfEvent(anEvent);
}
