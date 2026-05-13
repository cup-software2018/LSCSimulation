#include <iostream>

#include "MCObjs/MCTrack.hh"
#include "MCObjs/MCStep.hh"


ClassImp(MCTrack)

MCTrack::MCTrack()
    : TClonesArray("MCStep")
{
}

MCTrack::MCTrack(const MCTrack & trk)
    : TClonesArray(trk)
    , fPDGCode(trk.GetPDGCode())
    , fTrackId(trk.GetTrackId())
    , fParentId(trk.GetParentId())
    , fKineticEnergy(trk.GetKineticEnergy())
    , fGlobalTime(trk.GetGlobalTime())
    , fLocalTime(trk.GetLocalTime())
    , fParticleName(trk.GetParticleName())
    , fProcessName(trk.GetProcessName())
{
  trk.GetVertex(fVx, fVy, fVz);
}

MCTrack::~MCTrack() = default;

MCStep * MCTrack::AddStep()
{
  return new ((*this)[GetEntriesFast()]) MCStep();
}

void MCTrack::Clear(Option_t * opt)
{
  TClonesArray::Clear("C");
}

void MCTrack::Print(Option_t * opt) const
{
  std::cout << "+++++++++ TrackId = " << fTrackId << std::endl;
  std::cout << "         MotherId = " << fParentId << std::endl;
  std::cout << "     ParticleName = " << fParticleName << std::endl;
  std::cout << "     ProcessName  = " << fProcessName << std::endl;
  std::cout << Form("     Global Time  = %0.6f [ns]", fGlobalTime) << std::endl;
  std::cout << Form("   Kinetic Energy = %0.6f [MeV]", fKineticEnergy) << std::endl;
  std::cout << std::endl;
}
