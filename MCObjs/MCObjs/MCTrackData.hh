#ifndef MCTrackData_HH
#define MCTrackData_HH

#include "TClonesArray.h"

class MCTrack;
class MCTrackData : public TClonesArray {
public:
  MCTrackData();
  MCTrackData(const MCTrackData & data);
  virtual ~MCTrackData();

  virtual void Clear(const Option_t * opt = "");

  MCTrack * Add();

  int GetN() const;
  MCTrack * Get(int i) const;
  MCTrack * FindTrack(int id);
  MCTrack * GetParentTrack(MCTrack * track) const;

  virtual void Print(const Option_t * opt = "") const;

private:
  int fN; //!

  ClassDef(MCTrackData, 1)
};

inline int MCTrackData::GetN() const { return GetEntriesFast(); }
inline MCTrack * MCTrackData::Get(int n) const { return (MCTrack *)At(n); }

#endif
