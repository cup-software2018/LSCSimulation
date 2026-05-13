#ifndef PMTHit_h
#define PMTHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "TClonesArray.h"

class MCPhotonHit;
class PMTHit : public G4VHit {
public:
  PMTHit();
  PMTHit(const PMTHit & pmt);
  virtual ~PMTHit();

  const PMTHit & operator=(const PMTHit & pmt);
  int operator==(const PMTHit & pmt) const;

  void * operator new(size_t);
  void operator delete(void * aHit);

  void SetPMTId(int n) { fPMTId = n; }
  int GetPMTId() const { return fPMTId; }

  MCPhotonHit * AddHit();

  int GetNHit() const;
  MCPhotonHit * GetHit(int i) const;
  TClonesArray * GetPhotonHitColl() const;

private:
  int fPMTId;
  int fNHit;

  G4ThreeVector fPMTPosition;
  TClonesArray * fPhotonHitColl;
};

inline int PMTHit::GetNHit() const { return fPhotonHitColl->GetEntriesFast(); }
inline MCPhotonHit * PMTHit::GetHit(int n) const
{
  return (MCPhotonHit *)fPhotonHitColl->At(n);
}
inline TClonesArray * PMTHit::GetPhotonHitColl() const
{
  return fPhotonHitColl;
}

typedef G4THitsCollection<PMTHit> PMTHitsCollection;
extern G4Allocator<PMTHit> PMTHitAllocator;

#endif
