R__LOAD_LIBRARY(libHist)
R__LOAD_LIBRARY(libTrgSim)

void digi_example()
{
  auto pmt = new PMTSignal();
  pmt->SetTransitTime(100);
  pmt->SetTTS(1.2);

  pmt->AddHitTime(0);
  pmt->SetModel(PMTSignal::MOYAL);
  // pmt->Draw();

  auto fadc = new FADCWaveformGenerator();
  fadc->SetADCBits(12);
  fadc->SetVpp(2500);
  fadc->SetSamplingRate(500);
  fadc->SetRecordLength(250);
  fadc->SetPedestalMean(100);
  fadc->SetPedestalRMS(0.5);  // ADC counts (formerly mV)

  fadc->SetSignal(pmt);
  fadc->Digitize();

  fadc->Draw();
}