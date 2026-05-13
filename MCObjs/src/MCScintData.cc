#include <iostream>

#include "MCScint.hh"
#include "MCScintData.hh"

ClassImp(MCScintData)

MCScintData::MCScintData()
  : TClonesArray("MCScint")
{
}

MCScintData::MCScintData(const MCScintData & data)
  : TClonesArray(data)
{
}

MCScintData::~MCScintData() = default;

MCScint * MCScintData::Add() { return new ((*this)[GetEntriesFast()]) MCScint(); }

MCScint * MCScintData::Add(int id) { return new ((*this)[GetEntriesFast()]) MCScint(id); }

void MCScintData::Clear(Option_t * opt) { TClonesArray::Clear("C"); }

MCScint * MCScintData::FindScint(int id)
{
  for (auto obj : *this) {
    auto scint = static_cast<MCScint *>(obj);
    if (scint->GetVolumeId() == id) return scint;
  }
  return nullptr;
}

void MCScintData::Print(Option_t * opt) const
{
  int n = GetN();
  std::cout << Form("===> MCScintData: number of scintillation volume: %d", n) << std::endl;
  for (auto obj : *this) {
    static_cast<MCScint *>(obj)->Print();
  }
}
