#include "MCObjs/MCPMTData.hh"
#include "MCObjs/MCPMT.hh"

ClassImp(MCPMTData)

MCPMTData::MCPMTData()
  : TClonesArray("MCPMT")
{
  fN = 0;
}

MCPMTData::MCPMTData(const MCPMTData & data)
  : TClonesArray(data)
{
}

MCPMTData::~MCPMTData() {}

MCPMT * MCPMTData::Add() { return new ((*this)[fN++]) MCPMT(); }

void MCPMTData::Clear(const Option_t * opt)
{
  fN = 0;
  Delete();
}

void MCPMTData::Print(const Option_t * opt) const
{
  int n = GetN();
  for (int i = 0; i < n; i++) {
    MCPMT * track = Get(i);
    track->Print();
  }
}
