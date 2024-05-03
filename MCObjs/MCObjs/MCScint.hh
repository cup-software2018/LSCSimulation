#ifndef MCScint_hh
#define MCScint_hh

#include "TClonesArray.h"

class MCScintStep;
class MCScint : public TClonesArray {
public:
  MCScint();
  MCScint(int id);
  MCScint(const MCScint & scint);
  virtual ~MCScint();

  virtual void Clear(const Option_t * opt = "");

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

  virtual void Print(const Option_t * opt = "") const;

private:
  int fVolumeId;
  int fNScintPhoton;
  float fEdep;
  float fEdepQuenched;

  // Step
  int fNStep; //!

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
  return (MCScintStep *)At(n);
}

#endif
