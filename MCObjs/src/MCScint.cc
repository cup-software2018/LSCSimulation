#include <iostream>

#include "MCObjs/MCScint.hh"

ClassImp(MCScint)

MCScint::MCScint()
  : TObject()
{
}

MCScint::MCScint(int id)
  : TObject(),
    fVolumeId(id)
{
}

MCScint::MCScint(const MCScint & scint)
  : TObject(scint),
    fVolumeId(scint.GetVolumeId()),
    fNScintPhoton(scint.GetNScintPhoton()),
    fEdep(scint.GetEnergyDeposit()),
    fEdepQuenched(scint.GetEnergyVisible()),
    fSteps(scint.fSteps)
{
}

MCScint::~MCScint() = default;

MCScintStep * MCScint::AddStep()
{
  fSteps.emplace_back();
  return &fSteps.back();
}

void MCScint::Clear(Option_t * opt)
{
  fEdep = 0;
  fEdepQuenched = 0;
  fNScintPhoton = 0;
  fSteps.clear();
}

void MCScint::Print(Option_t * opt) const
{
  std::cout << Form("Volume: %d ", fVolumeId) << std::endl;
  std::cout << Form("   Deposit Energy    = %.6f [MeV]", GetEnergyDeposit()) << std::endl;
  std::cout << Form("   Visible Energy    = %.6f [MeV]", GetEnergyVisible()) << std::endl;
  std::cout << Form("   # of Scint Photon = %d         ", GetNScintPhoton()) << std::endl;

  if (!fSteps.empty()) {
    std::cout << Form("   # of Scint Step   = %d         ", GetNStep()) << std::endl;
  }
}
