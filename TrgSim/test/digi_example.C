R__LOAD_LIBRARY(libHist)
R__LOAD_LIBRARY(libTrgSim)

void digi_example()
{
  auto pmt = new PMTSignal();
  pmt->SetModel(PMTSignal::MOYAL);
  pmt->SetTransitTime(100);
  pmt->SetTTS(1.2);
  pmt->SetMeanCharge(1.6);
  pmt->SetChargeRMS(0.4);
  pmt->SetBackscatterFraction(0.3);

  auto fadc = new FADCWaveformGenerator();
  fadc->SetADCBits(12);
  fadc->SetVpp(2500);
  fadc->SetSamplingRate(500);
  fadc->SetRecordLength(250);
  fadc->SetPedestalMean(100);
  fadc->SetPedestalRMS(0.5); // ADC counts (formerly mV)

  fadc->SetSignal(pmt);

  pmt->AddHitTime(0);
  fadc->Digitize();

  fadc->Draw();
}