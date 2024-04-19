#include "LSCSim/PMTHit.hh"

#include "MCObjs/MCPhotonHit.hh"

G4Allocator<PMTHit> PMTHitAllocator;

PMTHit::PMTHit()
    : G4VHit()
{
  fPMTId = -1;
  fNHit = 0;
  fPhotonHitColl = new TClonesArray("MCPhotonHit");
}

PMTHit::PMTHit(const PMTHit & pmt)
    : G4VHit()
{
  fPMTId = pmt.GetPMTId();
  fNHit = pmt.GetNHit();
  fPhotonHitColl = pmt.GetPhotonHitColl();
}

PMTHit::~PMTHit()
{
  fPhotonHitColl->Delete();
  delete fPhotonHitColl;
}

MCPhotonHit * PMTHit::AddHit()
{
  return new ((*fPhotonHitColl)[fNHit++]) MCPhotonHit();
}

const PMTHit & PMTHit::operator=(const PMTHit & pmt)
{
  fPMTId = pmt.GetPMTId();
  fNHit = pmt.GetNHit();

  return *this;
}

G4int PMTHit::operator==(const PMTHit & pmt) const
{
  return (fPMTId == pmt.GetPMTId());
}

void * PMTHit::operator new(size_t)
{
  void * aHit;
  aHit = (void *)PMTHitAllocator.MallocSingle();
  return aHit;
}

void PMTHit::operator delete(void * aHit)
{
  PMTHitAllocator.FreeSingle((PMTHit *)aHit);
}