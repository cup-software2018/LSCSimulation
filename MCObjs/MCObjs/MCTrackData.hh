#pragma once

#include "TClonesArray.h"

#include "MCObjs/MCTrack.hh"
class MCTrackData : public TClonesArray {
public:
  MCTrackData();
  MCTrackData(const MCTrackData & data);
  virtual ~MCTrackData();

  void Clear(Option_t * opt = "") override;

  MCTrack * Add();

  int GetN() const;
  MCTrack * Get(int i) const;
  MCTrack * FindTrack(int id) const;
  MCTrack * GetParentTrack(MCTrack * track) const;

  void Print(Option_t * opt = "") const override;

  ClassDef(MCTrackData, 1)
};

inline int MCTrackData::GetN() const { return GetEntriesFast(); }
inline MCTrack * MCTrackData::Get(int n) const { return static_cast<MCTrack *>(At(n)); }

