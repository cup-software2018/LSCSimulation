#pragma once

#include <vector>

#include "AbsSignal.hh"

class PMTSignal : public AbsSignal {
public:
  enum Model
  {
    GAUSS = 1,
    MOYAL = 2
  };

  PMTSignal();
  PMTSignal(int pid);

  void SetPID(int pid);

  void SetModel(Model mod);
  void SetMeanCharge(double pC);
  void SetChargeRMS(double pC);
  void SetBackscatterFraction(double val);
  void SetFallTime(double val);
  void SetTransitTime(double val);
  void SetTTS(double val);
  void SetDarkRate(double val);
  void IncludeDark();

  int GetPID() const;

  void AddHitTime(double val) override;
  void Initialize() override;
  void Prepare() override;
  double EvalPulseHeight(double t) override;

private:
  double PMTPulse(double t) const;
  double Moyal(double t) const;
  double Gauss(double t) const;

private:
  Model fModel = GAUSS;

  int fPID = 0;

  double fMeanCharge = 1.6;  // pC
  double fChargeRMS  = 0.0;  // pC
  double fBackscatterFraction = 0.0;  // η: fraction of partially amplified PE
  double fFallTime = 4.0;     // ns
  double fTransitTime = 50.0; // ns
  double fTTS = 10.0;         // ns
  double fDarkRate = 5000;    // Hz

  bool fIncludeDark = false;

  std::vector<double> fHitTime;
  std::vector<double> fPulseTime;
  std::vector<double> fPulseSize;

  ClassDef(PMTSignal, 2)
};

inline void PMTSignal::SetModel(Model mod) { fModel = mod; }
inline void PMTSignal::SetMeanCharge(double pC) { fMeanCharge = pC; }
inline void PMTSignal::SetChargeRMS(double pC) { fChargeRMS = pC; }
inline void PMTSignal::SetBackscatterFraction(double val) { fBackscatterFraction = val; }
inline void PMTSignal::SetFallTime(double val) { fFallTime = val; }
inline void PMTSignal::SetTransitTime(double val) { fTransitTime = val; }
inline void PMTSignal::SetTTS(double val) { fTTS = val; }
inline void PMTSignal::SetDarkRate(double val) { fDarkRate = val; }

inline void PMTSignal::IncludeDark() { fIncludeDark = true; }

inline void PMTSignal::Initialize()
{
  fHitTime.clear();
  fIsPrepared = false;
}

inline void PMTSignal::SetPID(int pid) { fPID = pid; }
inline int PMTSignal::GetPID() const { return fPID; }
