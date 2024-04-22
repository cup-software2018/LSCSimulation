#include "MCObjs/MCPrimaryData.hh"
#include "MCObjs/MCPrimary.hh"

ClassImp(MCPrimaryData)

MCPrimaryData::MCPrimaryData()
  : TClonesArray("MCPrimary")
{
  fN = 0;
}

MCPrimaryData::MCPrimaryData(const MCPrimaryData & data)
  : TClonesArray(data)
{
}

MCPrimaryData::~MCPrimaryData() {}

MCPrimary * MCPrimaryData::Add() { return new ((*this)[fN++]) MCPrimary(); }

void MCPrimaryData::Clear(const Option_t * opt)
{
  fN = 0;
  Delete();
}

void MCPrimaryData::Print(const Option_t * opt) const
{
  int n = GetN();
  for (int i = 0; i < n; i++) {
    MCPrimary * track = Get(i);
    track->Print();
  }
}
