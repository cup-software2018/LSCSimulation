/*
 *
 *  Module:  FADCWaveformGenerator
 *
 *  Author:  Jaison Lee
 *
 *  Purpose: ADC digitization simulation
 *
 *  Last Update:      $Author: cupsoft $
 *  Update Date:      $Date: 2018/12/28 00:21:15 $
 *  CVS/RCS Revision: $Revision: 1.6 $
 *  Status:           $State: Exp $
 *
 */

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
  void SetVpp(double val);
  void SetBinTimeWidth(double val);
  void SetNBin(int val);

  void SetPedOffset(double val);
  void SetPedRMS(double val);
  void SetTermination(double val);

  void SetSignal(AbsSignal * signal);

  void Digitize();

  const unsigned int * GetWaveform() const;
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

  unsigned int * fWaveform;

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

inline void FADCWaveformGenerator::SetBinTimeWidth(double val)
{
  fBinTimeWidth = val;
}

inline void FADCWaveformGenerator::SetNBin(int val) { fNBin = val; }

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

inline const unsigned int * FADCWaveformGenerator::GetWaveform() const
{
  return fWaveform;
}

inline TH1D * FADCWaveformGenerator::GetWaveformHist() const
{
  return fWaveformHist;
}

#endif
