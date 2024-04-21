#include "TrgSim/PMTSignal.hh"

#include <algorithm>
#include <iostream>

#include "TMath.h"

using namespace std;

ClassImp(PMTSignal)

PMTSignal::PMTSignal()
{
  fPID = 0;

  fModel = GAUSS;

  fGain = 1.0e+7;
  fScale = 1.0;
  fFallTime = 4.0;     // ns
  fTransitTime = 50.0; // ns
  fTTS = 10.0;         // ns
  fDarkRate = 5000;    // Hz
  fQEff = 0;

  fIncludeDark = false;

  fNPE = 0;
}

PMTSignal::PMTSignal(int pid)
{
  fPID = pid;

  fModel = GAUSS;

  fGain = 1.0e+7;
  fScale = 1.0;
  fFallTime = 4.0;     // ns
  fTransitTime = 50.0; // ns
  fTTS = 10.0;         // ns
  fDarkRate = 5000;    // Hz
  fQEff = 0;

  fIncludeDark = false;

  fNPE = 0;
}

PMTSignal::~PMTSignal() {}

void PMTSignal::AddHitTime(double val)
{
  if (fNPE == kNMAXHIT) {
    Warning("AddHitTime",
            "Warning! NPE should not be greater than NMAXHIT (= %d)", kNMAXHIT);
    return;
  }

  bool ishit = true;
  if (fQEff > 0) {
    if (fRandom->Rndm() > fQEff) { ishit = false; }
  }

  if (ishit) {
    fHitTime[fNPE] = val;
    fNPE += 1;
  }
}

void PMTSignal::Prepare()
{
  fMinimumTime = fMaximumTime = 0;

  if (fNPE == 0) return;

  std::sort(fHitTime, fHitTime + fNPE);

  // Dark hit
  if (fIncludeDark) {
    double lastHitTime = fHitTime[fNPE - 1] + 30 * fFallTime;
    int ndark = fRandom->Poisson(fDarkRate * lastHitTime * 1.0e-9);

    for (int i = 0; i < ndark; i++) {
      double hitTime = lastHitTime * fRandom->Rndm();
      AddHitTime(hitTime);
    }

    std::sort(fHitTime, fHitTime + fNPE);
  }

  double pC = TMath::Qe() * 1e+12;

  for (int i = 0; i < fNPE; i++) {
    fPulseTime[i] = fHitTime[i] + fRandom->Gaus(fTransitTime, fTTS);

    fPulseSize[i] = (fScale == 0)
                        ? fGain
                        : fRandom->Gaus(fGain, fScale * TMath::Sqrt(fGain));
    fPulseSize[i] = pC * fPulseSize[i]; // pC
  }

  fMinimumTime = fPulseTime[0] - 10 * fFallTime;
  if (fMinimumTime < 0) fMinimumTime = 0;

  fMaximumTime = fPulseTime[fNPE - 1] + 30 * fFallTime;

  fIsPrepared = true;
}

double PMTSignal::Moyal(double t) const
{
  double lambda = fFallTime;

  // it's enough time for normalization
  if (t > 30 * lambda || t < -5 * lambda) return 0.;

  double arg = -t / lambda;

  // TMath::Sqrt(2 * TMath::Pi()) = 2.506628274631000242
  return TMath::Exp(0.5 * (arg - TMath::Exp(arg))) /
         (lambda * 2.506628274631000242);
}

double PMTSignal::Gauss(double t) const
{
  double sigma = fFallTime;

  // it's enough time for normalization
  if (t > 10 * sigma || t < -10 * sigma) return 0.;

  return TMath::Gaus(t, 0, sigma, kTRUE);
}

double PMTSignal::PMTPulse(double t)
{
  if (!fIsPrepared) Prepare();

  double val = 0;

  switch (fModel) {
    case GAUSS: val = Gauss(t); break;
    case MOYAL: val = Moyal(t); break;
    default: break;
  }

  return val;
}

double PMTSignal::EvalPulseHeight(double t)
{
  if (fNPE == 0) return 0;

  double val = 0;
  for (int i = 0; i < fNPE; i++) {
    val += fPulseSize[i] * PMTPulse(t - fPulseTime[i]);
  }

  return val;
}
