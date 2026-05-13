#include "FADCWaveformGenerator.hh"

ClassImp(FADCWaveformGenerator)

FADCWaveformGenerator::FADCWaveformGenerator()
  : TObject()
{
  fRandom = new TRandom3;
  fRandom->SetSeed(0);

  fWaveformHist = new TH1D();
  fWaveformHist->SetNameTitle("WaveformHist", "");
}

FADCWaveformGenerator::~FADCWaveformGenerator()
{
  delete fRandom;
  delete fWaveformHist;
}

void FADCWaveformGenerator::Prepare()
{
  fBinTimeWidth = 1000.0 / fSamplingRate;
  fNdp = static_cast<int>(fRecordLength / fBinTimeWidth);

  fWaveform.resize(fNdp);
  fWaveformHist->SetBins(fNdp, 0, fNdp);
  fWaveformHist->Reset();

  if (fSignal) { fSignal->Prepare(); }
}

void FADCWaveformGenerator::Digitize()
{
  Prepare();

  const double covfactor = fResolution * fTermination / fVpp;
  const int maxADC = TMath::Nint(fResolution);

  for (int j = 1; j <= fNdp; j++) {
    double time = (j - 0.5) * fBinTimeWidth;
    double sig = fRandom->Gaus(fPedMean, fPedRMS);
    if (fSignal) { sig += covfactor * fSignal->EvalPulseHeight(time); }
    if (sig < 0) sig = 0;

    fWaveform[j - 1] =
        static_cast<unsigned short>(TMath::Min(TMath::Nint(sig), maxADC));
    fWaveformHist->SetBinContent(j, fWaveform[j - 1]);
  }
}

void FADCWaveformGenerator::Draw(Option_t * option) { fWaveformHist->Draw(option); }
