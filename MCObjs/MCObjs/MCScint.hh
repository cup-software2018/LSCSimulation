#ifndef MCScint_hh
#define MCScint_hh

#include "TClonesArray.h"

#include "MCObjs/MCScintStep.hh"
class MCScint : public TClonesArray {
public:
  MCScint();
  MCScint(int id);
  MCScint(const MCScint & scint);
  virtual ~MCScint();

  void Clear(Option_t * opt = "") override;

  void SetVolumeId(int val);
  void AddEnergyDeposit(float val);
  void AddEnergyVisible(float val);
  void AddScintPhotons(int val);

  int GetVolumeId() const;
  float GetEnergyDeposit() const;
  float GetEnergyVisible() const;
  int GetNScintPhoton() const;

  MCScintStep * AddStep();
  int GetNStep() const;
  MCScintStep * GetStep(int i) const;

  void Print(Option_t * opt = "") const override;

private:
  int fVolumeId = -1;
  int fNScintPhoton = 0;
  float fEdep = 0;
  float fEdepQuenched = 0;

  ClassDef(MCScint, 1)
};

inline void MCScint::SetVolumeId(int val) { fVolumeId = val; }

inline void MCScint::AddEnergyDeposit(float val) { fEdep += val; }

inline void MCScint::AddEnergyVisible(float val) { fEdepQuenched += val; }

inline void MCScint::AddScintPhotons(int val) { fNScintPhoton += val; }

inline int MCScint::GetVolumeId() const { return fVolumeId; }

inline float MCScint::GetEnergyDeposit() const { return fEdep; }

inline float MCScint::GetEnergyVisible() const { return fEdepQuenched; }

inline int MCScint::GetNScintPhoton() const { return fNScintPhoton; }

inline int MCScint::GetNStep() const { return GetEntriesFast(); }
inline MCScintStep * MCScint::GetStep(int n) const
{
  return static_cast<MCScintStep *>(At(n));
}

#endif
