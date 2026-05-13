#include "MCObjs/MCScint.hh"

#include <iostream>

#include "MCObjs/MCScintStep.hh"


ClassImp(MCScint)

MCScint::MCScint()
    : TClonesArray("MCScintStep")
{
}

MCScint::MCScint(int id)
    : TClonesArray("MCScintStep")
    , fVolumeId(id)
{
}

MCScint::MCScint(const MCScint & scint)
    : TClonesArray(scint)
    , fVolumeId(scint.GetVolumeId())
    , fNScintPhoton(scint.GetNScintPhoton())
    , fEdep(scint.GetEnergyDeposit())
    , fEdepQuenched(scint.GetEnergyVisible())
{
}

MCScint::~MCScint() = default;

MCScintStep * MCScint::AddStep()
{
  return new ((*this)[GetEntriesFast()]) MCScintStep();
}

void MCScint::Clear(Option_t * opt)
{
  fEdep = 0;
  fEdepQuenched = 0;
  fNScintPhoton = 0;
  TClonesArray::Clear("C");
}

void MCScint::Print(Option_t * opt) const
{
  std::cout << Form("Volume: %d ", fVolumeId) << std::endl;
  std::cout << Form("   Deposit Energy    = %.6f [MeV]", GetEnergyDeposit()) << std::endl;
  std::cout << Form("   Visible Energy    = %.6f [MeV]", GetEnergyVisible()) << std::endl;
  std::cout << Form("   # of Scint Photon = %d         ", GetNScintPhoton()) << std::endl;

  if (GetNStep() > 0) {
    std::cout << Form("   # of Scint Step   = %d         ", GetNStep()) << std::endl;
  }
}
