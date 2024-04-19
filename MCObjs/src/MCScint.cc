#include <iostream>

#include "MCObjs/MCScintStep.hh"
#include "MCObjs/MCScint.hh"

using namespace std;

ClassImp(MCScint)

    MCScint::MCScint()
    : TClonesArray("MCScintStep")
{
  fVolumeId = 0;

  // Scint
  fEdep = 0;
  fEdepQuenched = 0;
  fNScintPhoton = 0;

  fNStep = 0;
}

MCScint::MCScint(int id)
    : TClonesArray("MCScintStep")
{
  fVolumeId = id;

  // Scint
  fEdep = 0;
  fEdepQuenched = 0;
  fNScintPhoton = 0;

  fNStep = 0;
}

MCScint::MCScint(const MCScint & scint)
    : TClonesArray(scint)
{
  fVolumeId = scint.GetVolumeId();

  // Scint
  fEdep = scint.GetEnergyDeposit();
  fEdepQuenched = scint.GetEnergyVisible();
  fNScintPhoton = scint.GetNScintPhoton();
}

MCScint::~MCScint() {}

MCScintStep * MCScint::AddStep() { return new ((*this)[fNStep++]) MCScintStep(); }

void MCScint::Clear(const Option_t * opt)
{
  fNStep = 0;
  Delete();
}

void MCScint::print(const Option_t * opt) const
{
  cout << Form("   Deposit Energy    = %6.3f [MeV]", GetEnergyDeposit())
       << endl;
  cout << Form("   Visible Energy    = %6.3f [MeV]", GetEnergyVisible())
       << endl;
  cout << Form("   # of Scint Photon = %d         ", GetNScintPhoton()) << endl;

  if (GetNStep() > 0) {
    cout << Form("   # of Scint Step   = %d         ", GetNStep())
         << endl;
  }
}
