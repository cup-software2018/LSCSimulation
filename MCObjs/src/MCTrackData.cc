#include "MCObjs/MCTrackData.hh"
#include "MCObjs/MCTrack.hh"

ClassImp(MCTrackData)

MCTrackData::MCTrackData()
    : TClonesArray("MCTrack")
{
  fN = 0;
}

MCTrackData::MCTrackData(const MCTrackData & data)
    : TClonesArray(data)
{
}

MCTrackData::~MCTrackData() {}

MCTrack * MCTrackData::Add() { return new ((*this)[fN++]) MCTrack(); }

MCTrack * MCTrackData::FindTrack(int id)
{
  MCTrack * track = nullptr;

  int n = GetN();
  for (int i = 0; i < n; i++) {
    MCTrack * rtrack = Get(i);
    if (rtrack->GetTrackId() == id) {
      track = rtrack;
      break;
    }
  }

  return track;
}

MCTrack * MCTrackData::GetParentTrack(MCTrack * track) const
{
  MCTrack * mother = nullptr;

  int pid = track->GetParentId();

  int n = GetN();
  for (int i = 0; i < n; i++) {
    MCTrack * rtrack = Get(i);
    if (rtrack->GetParentId() == pid) {
      mother = rtrack;
      break;
    }
  }

  return mother;
}

void MCTrackData::Clear(const Option_t * opt)
{
  fN = 0;
  Delete();
}

void MCTrackData::Print(const Option_t * opt) const
{
  int n = GetN();
  for (int i = 0; i < n; i++) {
    MCTrack * track = Get(i);
    track->Print();
  }
}
