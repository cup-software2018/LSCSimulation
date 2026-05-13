#ifndef PMTSignal_hh
#define PMTSignal_hh

#include <iostream>

#include "TrgSim/AbsSignal.hh"

const int kNMAXHIT = 20000;

class PMTSignal : public AbsSignal {
public:
  enum Model { GAUSS = 1, MOYAL = 2 };

  PMTSignal();
  PMTSignal(int pid);
  virtual ~PMTSignal();

  void SetPID(int pid);

  void SetModel(Model mod);
  void SetGain(double val);
  void SetGainScale(double val);
  void SetFallTime(double val);
  void SetTransitTime(double val);
  void SetTTS(double val);
  void SetQuantumEfficiency(double val);
  void SetDarkRate(double val);
  void IncludeDark();
  
  int GetPID() const;

  virtual void AddHitTime(double val);
  virtual void Initialize();
  virtual void Prepare();
  virtual double EvalPulseHeight(double t);

private:
  double PMTPulse(double t);
  double Moyal(double t) const;
  double Gauss(double t) const;

private:
  Model fModel;

  int fPID;

  double fGain;
  double fScale;
  double fFallTime;
  double fTransitTime;
  double fTTS;
  double fDarkRate;
  double fQEff;

  bool fIncludeDark;

  int fNPE;
  double fHitTime[kNMAXHIT];
  double fPulseTime[kNMAXHIT];
  double fPulseSize[kNMAXHIT];

  ClassDef(PMTSignal, 1)
};

inline void PMTSignal::SetModel(Model mod) { fModel = mod; }

inline void PMTSignal::SetGain(double val) { fGain = val; }

inline void PMTSignal::SetGainScale(double val) { fScale = val; }

inline void PMTSignal::SetFallTime(double val) { fFallTime = val; }

inline void PMTSignal::SetTransitTime(double val) { fTransitTime = val; }

inline void PMTSignal::SetTTS(double val) { fTTS = val; }

inline void PMTSignal::SetDarkRate(double val) { fDarkRate = val; }

inline void PMTSignal::SetQuantumEfficiency(double val) { fQEff = val; }

inline void PMTSignal::IncludeDark() { fIncludeDark = true; }

inline void PMTSignal::Initialize()
{
  fNPE = 0;
  fIsPrepared = false;
}

inline void PMTSignal::SetPID(int pid) { fPID = pid; }

inline int PMTSignal::GetPID() const { return fPID; }

#endif
