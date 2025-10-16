#include "LSCSim/LSCTrackingAction.hh"

#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"

#include "LSCSim/LSCRecorderBase.hh"
#include "LSCSim/LSCTrajectory.hh"
#include "LSCSim/LSCUserTrackInformation.hh"

LSCTrackingAction::LSCTrackingAction(LSCRecorderBase * r)
    : fRecorder(r)
{
}

void LSCTrackingAction::PreUserTrackingAction(const G4Track * aTrack)
{
  // Use custom trajectory class
  fpTrackingManager->SetTrajectory(new LSCTrajectory(aTrack));

  // This user track information is only relevant to the photons
  //fpTrackingManager->SetUserTrackInformation(new LSCUserTrackInformation);
  aTrack->SetUserInformation(new LSCUserTrackInformation);
}

void LSCTrackingAction::PostUserTrackingAction(const G4Track * aTrack)
{
  LSCTrajectory * trajectory =
      (LSCTrajectory *)fpTrackingManager->GimmeTrajectory();
  LSCUserTrackInformation * trackInformation =
      (LSCUserTrackInformation *)aTrack->GetUserInformation();

  if (trackInformation->GetTrackStatus() == LSCTrackStatus::deferred) return;

  if (fRecorder) fRecorder->RecordTrack(aTrack);      

  // Let's choose to draw only the photons that hit the sphere and a pmt
  if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    const G4VProcess * creator = aTrack->GetCreatorProcess();
    if (creator && creator->GetProcessName() == "OpWLS") {
      trajectory->WLS();
      trajectory->SetDrawTrajectory(true);
    }

    //if (trackInformation->GetTrackStatus() & hitPMT)
      trajectory->SetDrawTrajectory(true);
  }
  // draw all other (not optical photon) trajectories
  else trajectory->SetDrawTrajectory(true);

  if (trackInformation->GetForceDrawTrajectory())
    trajectory->SetDrawTrajectory(true);
}
