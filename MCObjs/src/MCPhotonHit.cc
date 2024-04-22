#include "MCObjs/MCPhotonHit.hh"

ClassImp(MCPhotonHit)

MCPhotonHit::MCPhotonHit()
  : TObject()
{
  fTime = -1;
  fKE = -1;
}

MCPhotonHit::MCPhotonHit(const MCPhotonHit & photon)
  : TObject()
{
  fTime = photon.GetTime();
  fKE = photon.GetKineticEnergy();
}

MCPhotonHit::~MCPhotonHit() {}
