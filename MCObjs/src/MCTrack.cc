#include <iostream>

#include "MCObjs/MCTrack.hh"
#include "MCObjs/MCStep.hh"

using namespace std;

ClassImp(MCTrack)

MCTrack::MCTrack()
    : TClonesArray("MCStep")
{
  fParticleName = "";
  fPDGCode = 0;

  fTrackId = 0;
  fParentId = 0;

  fVx = 0;
  fVy = 0;
  fVz = 0;
  fKineticEnergy = 0;
  fGlobalTime = 0;
  fLocalTime = 0;

  fProcessName = "";

  fNStep = 0;
}

MCTrack::MCTrack(const MCTrack & trk)
    : TClonesArray(trk)
{
  fParticleName = trk.GetParticleName();
  fPDGCode = trk.GetPDGCode();

  fTrackId = trk.GetTrackId();
  fParentId = trk.GetParentId();

  trk.GetVertex(fVx, fVy, fVz);
  fKineticEnergy = trk.GetKineticEnergy();
  fGlobalTime = trk.GetGlobalTime();
  fLocalTime = trk.GetLocalTime();

  fProcessName = trk.GetProcessName();
}

MCTrack::~MCTrack() {}

MCStep * MCTrack::AddStep() { return new ((*this)[fNStep++]) MCStep(); }

void MCTrack::Clear(const Option_t * opt)
{
  fNStep = 0;
  Delete();
}

void MCTrack::Print(const Option_t * opt) const
{
  cout << "+++++++++ TrackId = " << fTrackId << endl;
  cout << "         MotherId = " << fParentId << endl;
  cout << "     ParticleName = " << fParticleName << endl;
  cout << "     ProcessName  = " << fProcessName << endl;
  cout << Form("     Global Time  = %6.2f [ns]", fGlobalTime) << endl;
  cout << Form("   Kinetic Energy = %6.2f [MeV]", fKineticEnergy) << endl;
  cout << endl;
}
