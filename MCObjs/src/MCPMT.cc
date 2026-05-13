#include <iostream>
#include "MCObjs/MCPMT.hh"
#include "MCObjs/MCPhotonHit.hh"


ClassImp(MCPMT)

MCPMT::MCPMT()
    : TClonesArray("MCPhotonHit")
{
}

MCPMT::MCPMT(int id)
    : TClonesArray("MCPhotonHit")
    , fPMTId(id)
{
}

MCPMT::MCPMT(const MCPMT & pmt)
    : TClonesArray(pmt)
    , fPMTId(pmt.GetId())
{
}

MCPMT::~MCPMT() = default;

MCPhotonHit * MCPMT::AddHit()
{
  return new ((*this)[GetEntriesFast()]) MCPhotonHit();
}

MCPhotonHit * MCPMT::AddHit(MCPhotonHit * hit)
{
  return new ((*this)[GetEntriesFast()]) MCPhotonHit(*hit);
}

void MCPMT::Clear(Option_t * opt)
{
  TClonesArray::Clear("C");
}

int MCPMT::Compare(const TObject * object) const
{
  auto comp = static_cast<const MCPMT *>(object);
  if (this->GetId() > comp->GetId()) return 1;
  else if (this->GetId() < comp->GetId()) return -1;

  return 0;
}

void MCPMT::Print(Option_t * opt) const
{
  std::cout << Form("PMTID=%4d, # of hit=%d", fPMTId, GetNHit()) << std::endl;
}
