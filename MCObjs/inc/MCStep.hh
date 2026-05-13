#pragma once

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

class MCStep : public TObject {
public:
  MCStep();
  MCStep(const MCStep & step);
  virtual ~MCStep() = default;

  void SetStepLength(float val);
  void SetEnergyDeposit(float val);
  void SetEnergyDepositNonIonizing(float val);
  void SetGlobalTime(float val);
  void SetLocalTime(float val);
  void SetStepPoint(float x, float y, float z);
  void SetVolumeName(const char * name);

  float GetStepLength() const;
  float GetEnergyDeposit() const;
  float GetEnergyDepositNonIonizing() const;
  float GetGlobalTime() const;
  float GetLocalTime() const;
  void GetStepPoint(float & x, float & y, float & z) const;
  const char * GetVolumeName() const;

protected:
  float fStepLength = 0;
  float fEnergyDeposit = 0;
  float fEnergyDepositNonIonizing = 0;
  float fGlobalTime = 0;
  float fLocalTime = 0;
  float fX = 0, fY = 0, fZ = 0;
  TString fVolumeName;

  ClassDef(MCStep, 1)
};

inline void MCStep::SetStepLength(float val) { fStepLength = val; }

inline void MCStep::SetEnergyDeposit(float val) { fEnergyDeposit = val; }

inline void MCStep::SetEnergyDepositNonIonizing(float val) { fEnergyDepositNonIonizing = val; }

inline void MCStep::SetGlobalTime(float val) { fGlobalTime = val; }

inline void MCStep::SetLocalTime(float val) { fLocalTime = val; }

inline void MCStep::SetStepPoint(float x, float y, float z)
{
  fX = x;
  fY = y;
  fZ = z;
}

inline void MCStep::SetVolumeName(const char * val) { fVolumeName = val; }

inline float MCStep::GetStepLength() const { return fStepLength; }

inline float MCStep::GetEnergyDeposit() const { return fEnergyDeposit; }

inline float MCStep::GetEnergyDepositNonIonizing() const { return fEnergyDepositNonIonizing; }

inline float MCStep::GetGlobalTime() const { return fGlobalTime; }

inline float MCStep::GetLocalTime() const { return fLocalTime; }

inline void MCStep::GetStepPoint(float & x, float & y, float & z) const
{
  x = fX;
  y = fY;
  z = fZ;
}

inline const char * MCStep::GetVolumeName() const { return fVolumeName.Data(); }
