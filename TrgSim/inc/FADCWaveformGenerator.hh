#pragma once

#include <vector>

#include "TH1D.h"
#include "TMath.h"
#include "TObject.h"
#include "TRandom3.h"

#include "AbsSignal.hh"

class FADCWaveformGenerator : public TObject {
public:
  FADCWaveformGenerator();
  virtual ~FADCWaveformGenerator();

  void SetADCBits(int val);
  void SetVpp(double val);          // mV
  void SetSamplingRate(double val); // MHz
  void SetRecordLength(double val); // ns

  void SetPedestalMean(double val); // ADC counts
  void SetPedestalRMS(double val);  // ADC counts
  void SetTermination(double val);  // ohm

  void SetSignal(AbsSignal * signal);

  void Digitize();

  const unsigned short * GetWaveform() const;
  int GetNdp() const;
  TH1D * GetWaveformHist() const;

  virtual void Draw(Option_t * option = "");

private:
  void Prepare();

private:
  double fResolution = TMath::Power(2, 12); // ADC counts full scale (2^NBIT)
  double fVpp = 2500;                       // mV
  double fSamplingRate = 500;               // MHz
  double fRecordLength = 512;               // ns
  double fBinTimeWidth = 2;                 // ns, computed in Prepare()
  int fNdp = 256;                           // computed in Prepare()

  double fPedMean = 0;      // ADC counts
  double fPedRMS = 1;       // ADC counts
  double fTermination = 50; // ohm

  AbsSignal * fSignal = nullptr;

  std::vector<unsigned short> fWaveform;

  TH1D * fWaveformHist = nullptr;

  TRandom3 * fRandom = nullptr;

  ClassDef(FADCWaveformGenerator, 0)
};

inline void FADCWaveformGenerator::SetADCBits(int val) { fResolution = TMath::Power(2, val); }

inline void FADCWaveformGenerator::SetVpp(double val) { fVpp = val; }
inline void FADCWaveformGenerator::SetSamplingRate(double val) { fSamplingRate = val; }
inline void FADCWaveformGenerator::SetRecordLength(double val) { fRecordLength = val; }
inline void FADCWaveformGenerator::SetPedestalMean(double val) { fPedMean = val; }
inline void FADCWaveformGenerator::SetPedestalRMS(double val) { fPedRMS = val; }
inline void FADCWaveformGenerator::SetTermination(double val) { fTermination = val; }
inline void FADCWaveformGenerator::SetSignal(AbsSignal * signal) { fSignal = signal; }

inline const unsigned short * FADCWaveformGenerator::GetWaveform() const
{
  return fWaveform.data();
}

inline int FADCWaveformGenerator::GetNdp() const { return fNdp; }

inline TH1D * FADCWaveformGenerator::GetWaveformHist() const { return fWaveformHist; }
