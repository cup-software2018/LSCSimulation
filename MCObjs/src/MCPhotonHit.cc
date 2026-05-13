#include "MCPhotonHit.hh"

ClassImp(MCPhotonHit)

MCPhotonHit::MCPhotonHit()
  : TObject()
{
}

MCPhotonHit::MCPhotonHit(const MCPhotonHit & photon)
  : TObject(),
    fTime(photon.GetTime()),
    fKE(photon.GetKineticEnergy())
{
}

MCPhotonHit::~MCPhotonHit() = default;

int MCPhotonHit::Compare(const TObject * object) const
{
  auto comp = static_cast<const MCPhotonHit *>(object);
  if (this->GetTime() < comp->GetTime()) return -1;
  else if (this->GetTime() > comp->GetTime()) return 1;

  return 0;
}
