#ifndef FADCWaveformGenerator_hh
#define FADCWaveformGenerator_hh

#include "TH1D.h"
#include "TMath.h"
#include "TObject.h"
#include "TRandom3.h"

#include "TrgSim/AbsSignal.hh"

class FADCWaveformGenerator : public TObject {
public:
  FADCWaveformGenerator();
  virtual ~FADCWaveformGenerator();

  void SetNBIT(int val);
  void SetVpp(double val);          // mV
  void SetSamplingRate(double val); // MHz
  void SetTimeWindow(double val);   // ns

  void SetPedOffset(double val);   // ADC count
  void SetPedRMS(double val);      // mV
  void SetTermination(double val); // ohm

  void SetSignal(AbsSignal * signal);

  void Digitize();

  const unsigned short * GetWaveform() const;
  TH1D * GetWaveformHist() const;

  virtual void Draw(Option_t * option = "");

private:
  void Prepare();

private:
  double fResolution;
  double fVpp;
  double fBinTimeWidth;
  int fNBin;

  double fPedOffset;
  double fPedRMS;
  double fTermination;

  bool fEnableEXT;
  bool fEnablePCT;
  bool fEnablePWT;
  bool fEnablePST;

  AbsSignal * fSignal;

  unsigned short * fWaveform;

  TH1D * fWaveformHist;
  TH1D * fChannelTrgHist;

  TRandom3 * fRandom;

  ClassDef(FADCWaveformGenerator, 0)
};

inline void FADCWaveformGenerator::SetNBIT(int val)
{
  fResolution = TMath::Power(2, val);
}

inline void FADCWaveformGenerator::SetVpp(double val) { fVpp = val; }

inline void FADCWaveformGenerator::SetSamplingRate(double val)
{
  fBinTimeWidth = 1 / (val / 1000.);
}

inline void FADCWaveformGenerator::SetTimeWindow(double val)
{
  fNBin = val / fBinTimeWidth;
}

inline void FADCWaveformGenerator::SetPedOffset(double val)
{
  fPedOffset = val;
}

inline void FADCWaveformGenerator::SetPedRMS(double val) { fPedRMS = val; }

inline void FADCWaveformGenerator::SetTermination(double val)
{
  fTermination = val;
}

inline void FADCWaveformGenerator::SetSignal(AbsSignal * signal)
{
  fSignal = signal;
}

inline const unsigned short * FADCWaveformGenerator::GetWaveform() const
{
  return fWaveform;
}

inline TH1D * FADCWaveformGenerator::GetWaveformHist() const
{
  return fWaveformHist;
}

#endif
