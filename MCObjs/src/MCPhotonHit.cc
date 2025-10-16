#include "MCObjs/MCPhotonHit.hh"

ClassImp(MCPhotonHit)

MCPhotonHit::MCPhotonHit()
  : TObject()
{
  fTime = 0;
  fKE = 0;
}

MCPhotonHit::MCPhotonHit(const MCPhotonHit & photon)
  : TObject()
{
  fTime = photon.GetTime();
  fKE = photon.GetKineticEnergy();
}

MCPhotonHit::~MCPhotonHit() {}


int MCPhotonHit::Compare(const TObject * object) const
{
  auto comp = (MCPhotonHit*)object;
  if (this->GetTime() < comp->GetTime()) return 1;
  else if (this->GetTime() > comp->GetTime()) return -1;

  return 0;
}
