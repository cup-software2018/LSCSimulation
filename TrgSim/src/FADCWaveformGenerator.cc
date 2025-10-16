#include "TrgSim/FADCWaveformGenerator.hh"

#include <iostream>

using namespace std;

ClassImp(FADCWaveformGenerator)

    FADCWaveformGenerator::FADCWaveformGenerator()
  : TObject()
{
  fRandom = new TRandom3;
  fRandom->SetSeed(0);

  fResolution = TMath::Power(2, 12);
  fVpp = 2500;
  fBinTimeWidth = 2;
  fNBin = 256;

  fPedOffset = 0;
  fPedRMS = 1;

  fTermination = 50;

  fSignal = nullptr;
  fWaveform = nullptr;
  fWaveformHist = new TH1D();
  fWaveformHist->SetNameTitle("WaveformHist", "");
}

FADCWaveformGenerator::~FADCWaveformGenerator()
{
  delete fRandom;
  delete fWaveformHist;

  if (fWaveform) { delete[] fWaveform; }
}

void FADCWaveformGenerator::Prepare()
{
  if (!fWaveform) {
    fWaveform = new unsigned short[fNBin];
    fWaveformHist->SetBins(fNBin, 0, fNBin);
  }
  fWaveformHist->Reset();

  if (fSignal) { fSignal->Prepare(); }
}

void FADCWaveformGenerator::Digitize()
{
  Prepare();

  double covfactor = fResolution * fTermination / fVpp;

  int nbin = fWaveformHist->GetNbinsX();

  for (int j = 1; j <= nbin; j++) {
    double time = fWaveformHist->GetBinCenter(j) * fBinTimeWidth;
    double ped = fRandom->Gaus(fPedOffset, fPedRMS * fResolution / fVpp);
    double sig = 0;
    if (fSignal) { sig = covfactor * fSignal->EvalPulseHeight(time); }
    sig += ped;
    if (sig < 0) sig = 0;

    unsigned short digi = TMath::Min(TMath::Nint(sig), TMath::Nint(fResolution));

    fWaveform[j - 1] = digi;
    fWaveformHist->SetBinContent(j, digi);
  }
}

void FADCWaveformGenerator::Draw(Option_t * option)
{
  fWaveformHist->Draw(option);
}