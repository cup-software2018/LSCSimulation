#include "MCObjs/MCScintData.hh"
#include "MCObjs/MCScint.hh"

ClassImp(MCScintData)

MCScintData::MCScintData()
    : TClonesArray("MCScint")
{
  fN = 0;
}

MCScintData::MCScintData(const MCScintData & data)
    : TClonesArray(data)
{
}

MCScintData::~MCScintData() {}

MCScint * MCScintData::Add() { return new ((*this)[fN++]) MCScint(); }

void MCScintData::Clear(const Option_t * opt)
{
  fN = 0;
  Delete();
}

void MCScintData::Print(const Option_t * opt) const
{
  int n = GetN();
  for (int i = 0; i < n; i++) {
    MCScint * scint = Get(i);
    scint->Print();
  }
}
