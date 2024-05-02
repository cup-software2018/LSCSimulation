#include "MCObjs/MCScintData.hh"

#include <iostream>

#include "MCObjs/MCScint.hh"

using namespace std;

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
MCScint * MCScintData::Add(int id) { return new ((*this)[fN++]) MCScint(id); }

void MCScintData::Clear(const Option_t * opt)
{
  fN = 0;
  Delete();
}

MCScint * MCScintData::FindScint(int id)
{
  MCScint * scint = nullptr;

  int n = GetN();
  for (int i = 0; i < n; i++) {
    MCScint * rscint = Get(i);
    if (rscint->GetVolumeId() == id) {
      scint = rscint;
      break;
    }
  }

  return scint;
}

void MCScintData::Print(const Option_t * opt) const
{
  int n = GetN();
  cout << Form("===> MCScintData: number of scintillation volume: %d", n)
       << endl;
  for (int i = 0; i < n; i++) {
    MCScint * scint = Get(i);
    scint->Print();
  }
}
