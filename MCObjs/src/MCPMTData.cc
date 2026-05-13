#include "MCObjs/MCPMTData.hh"
#include "MCObjs/MCPMT.hh"

ClassImp(MCPMTData)

MCPMTData::MCPMTData()
    : TClonesArray("MCPMT")
{
}

MCPMTData::MCPMTData(const MCPMTData & data)
    : TClonesArray(data)
{
}

MCPMTData::~MCPMTData() = default;

MCPMT * MCPMTData::Add()
{
  return new ((*this)[GetEntriesFast()]) MCPMT();
}

void MCPMTData::Clear(Option_t * opt)
{
  TClonesArray::Clear("C");
}

void MCPMTData::Print(Option_t * opt) const
{
  for (auto obj : *this) {
    static_cast<MCPMT *>(obj)->Print();
  }
}
