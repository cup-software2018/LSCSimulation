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
