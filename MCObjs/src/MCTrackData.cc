#include "MCObjs/MCTrack.hh"
#include "MCObjs/MCTrackData.hh"

ClassImp(MCTrackData)

MCTrackData::MCTrackData()
  : TClonesArray("MCTrack")
{
}

MCTrackData::MCTrackData(const MCTrackData & data)
  : TClonesArray(data)
{
}

MCTrackData::~MCTrackData() = default;

MCTrack * MCTrackData::Add() { return new ((*this)[GetEntriesFast()]) MCTrack(); }

MCTrack * MCTrackData::FindTrack(int id) const
{
  for (auto obj : *this) {
    auto track = static_cast<MCTrack *>(obj);
    if (track->GetTrackId() == id) return track;
  }
  return nullptr;
}

MCTrack * MCTrackData::GetParentTrack(MCTrack * track) const
{
  return FindTrack(track->GetParentId());
}

void MCTrackData::Clear(Option_t * opt) { TClonesArray::Clear("C"); }

void MCTrackData::Print(Option_t * opt) const
{
  for (auto obj : *this) {
    static_cast<MCTrack *>(obj)->Print();
  }
}
