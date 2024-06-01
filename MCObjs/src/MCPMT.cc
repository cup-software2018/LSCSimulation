#include <iostream>
#include "MCObjs/MCPMT.hh"
#include "MCObjs/MCPhotonHit.hh"

using namespace std;

ClassImp(MCPMT)

MCPMT::MCPMT()
  : TClonesArray("MCPhotonHit")
{
  fPMTId = 0;
  fNHit = 0;
}

MCPMT::MCPMT(const MCPMT & pmt)
  : TClonesArray(pmt)
{
  fPMTId = pmt.GetId();
}

MCPMT::~MCPMT() {}

MCPhotonHit * MCPMT::AddHit() { return new ((*this)[fNHit++]) MCPhotonHit(); }

MCPhotonHit * MCPMT::AddHit(MCPhotonHit * hit)
{
  return new ((*this)[fNHit++]) MCPhotonHit(*hit);
}

void MCPMT::Clear(const Option_t * opt)
{
  fNHit = 0;
  Delete();
}

void MCPMT::Print(const Option_t * opt) const
{
  cout << Form("PMTID=%4d, # of hit=%d", fPMTId, GetNHit()) << endl;
}
