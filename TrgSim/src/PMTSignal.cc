#include <algorithm>

#include "TMath.h"
#include "PMTSignal.hh"

ClassImp(PMTSignal)

PMTSignal::PMTSignal()
  : PMTSignal(0)
{
}

PMTSignal::PMTSignal(int pid)
  : fPID(pid)
{
}

void PMTSignal::AddHitTime(double val)
{
  fHitTime.push_back(val);
}

void PMTSignal::Prepare()
{
  fMinimumTime = fMaximumTime = 0;

  if (fHitTime.empty()) return;

  std::sort(fHitTime.begin(), fHitTime.end());

  // Dark hits bypass QE — thermal emission, not photon detection
  if (fIncludeDark) {
    double windowEnd = fHitTime.back() + 30 * fFallTime;
    int ndark = fRandom->Poisson(fDarkRate * windowEnd * 1.0e-9);
    for (int i = 0; i < ndark; i++) {
      fHitTime.push_back(fRandom->Uniform(windowEnd));
    }
    std::sort(fHitTime.begin(), fHitTime.end());
  }

  const int npe = static_cast<int>(fHitTime.size());
  fPulseTime.resize(npe);
  fPulseSize.resize(npe);

  for (int i = 0; i < npe; i++) {
    fPulseTime[i] = fHitTime[i] + fRandom->Gaus(fTransitTime, fTTS);

    if (fBackscatterFraction > 0 && fRandom->Rndm() < fBackscatterFraction) {
      // Partially amplified (PA): back-scattered photoelectron deposits only a
      // fraction of its energy at the first dynode. Charge is uniform in [0, μ].
      fPulseSize[i] = fRandom->Uniform(fMeanCharge);
    } else {
      // Fully amplified (FA): Gaussian around mean charge
      double q = (fChargeRMS == 0) ? fMeanCharge
                                   : fRandom->Gaus(fMeanCharge, fChargeRMS);
      fPulseSize[i] = std::max(0.0, q);
    }
  }

  fMinimumTime = fPulseTime.front() - 10 * fFallTime;
  if (fMinimumTime < 0) fMinimumTime = 0;
  fMaximumTime = fPulseTime.back() + 30 * fFallTime;

  fIsPrepared = true;
}

double PMTSignal::Moyal(double t) const
{
  double lambda = fFallTime;
  if (t > 30 * lambda || t < -5 * lambda) return 0.;
  double arg = -t / lambda;
  // TMath::Sqrt(2 * TMath::Pi()) = 2.506628274631000242
  return TMath::Exp(0.5 * (arg - TMath::Exp(arg))) / (lambda * 2.506628274631000242);
}

double PMTSignal::Gauss(double t) const
{
  double sigma = fFallTime;
  if (t > 10 * sigma || t < -10 * sigma) return 0.;
  return TMath::Gaus(t, 0, sigma, kTRUE);
}

double PMTSignal::PMTPulse(double t) const
{
  switch (fModel) {
    case GAUSS: return Gauss(t);
    case MOYAL: return Moyal(t);
    default: return 0;
  }
}

double PMTSignal::EvalPulseHeight(double t)
{
  if (fHitTime.empty()) return 0;
  if (!fIsPrepared) Prepare();

  double val = 0;
  for (int i = 0; i < static_cast<int>(fPulseTime.size()); i++) {
    val += fPulseSize[i] * PMTPulse(t - fPulseTime[i]);
  }
  return val;
}
