#ifndef MCPhotonHit_hh
#define MCPhotonHit_hh

#include "TObject.h"

class MCPhotonHit : public TObject {
public:
  MCPhotonHit();
  MCPhotonHit(const MCPhotonHit & photon);
  ~MCPhotonHit();

  void SetTime(float t);
  void SetKineticEnergy(float KE);

  float GetTime() const;
  float GetKineticEnergy() const;
  float GetWavelength() const;

  virtual bool IsSortable() const { return true; }
  virtual int Compare(const TObject * object) const;

private:
  float fTime;
  float fKE;

  ClassDef(MCPhotonHit, 1)
};

inline void MCPhotonHit::SetTime(float t) { fTime = t; }
inline void MCPhotonHit::SetKineticEnergy(float KE) { fKE = KE; }
inline float MCPhotonHit::GetTime() const { return fTime; }
inline float MCPhotonHit::GetKineticEnergy() const { return fKE; }
inline float MCPhotonHit::GetWavelength() const { return 1239.83968E-6 / fKE; }

#endif
