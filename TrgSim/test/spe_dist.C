R__LOAD_LIBRARY(libHist)
R__LOAD_LIBRARY(libTrgSim)

#include <algorithm>

// ---------------------------------------------------------------------------
// SPE fit model
// ---------------------------------------------------------------------------
const double kSQRT2      = TMath::Sqrt(2);
const double kSQRTPIo2   = TMath::Sqrt(TMath::Pi() / 2);
const double kEpsilon    = std::numeric_limits<double>::epsilon();

double gFitLo, gFitHi;

double NormalizedGauss(double x, double m, double s, double lo, double up)
{
  double norm = kSQRTPIo2 * s * (TMath::Erf((m - lo) / kSQRT2 / s)
                                - TMath::Erf((m - up) / kSQRT2 / s));
  if (norm <= 0) return 0;
  return TMath::Gaus(x, m, s) / norm;
}

// Antiderivative of Erf((x-a)/(√2·s)) w.r.t. x
static double ErfAntideriv(double x, double a, double s)
{
  double u = (x - a) / (kSQRT2 * s);
  return (x - a) * TMath::Erf(u)
       + (kSQRT2 * s / TMath::Sqrt(TMath::Pi())) * TMath::Exp(-u * u);
}

// Gauss(q0, s0) ⊗ Uniform(0, q1), normalized over [lo, up]
// Models the backscatter (PA) charge distribution
double NormalizedBackscatter(double x, double q0, double s0, double q1,
                             double lo, double up)
{
  if (q1 <= 0) return 0;

  double val = 0.5 / q1 * (TMath::Erf((x - q0)        / (kSQRT2 * s0))
                          - TMath::Erf((x - q0 - q1)   / (kSQRT2 * s0)));

  double norm = 0.5 / q1 * ( ErfAntideriv(up, q0,      s0)
                            - ErfAntideriv(lo, q0,      s0)
                            - ErfAntideriv(up, q0 + q1, s0)
                            + ErfAntideriv(lo, q0 + q1, s0));
  if (norm <= 0) return 0;
  return val / norm;
}

// p[0]=q0, p[1]=s0, p[2]=q1, p[3]=s1, p[4]=mu, p[5]=eta, p[6]=nn
double fitfunc(double *v, double *p)
{
  double x   = v[0];
  double q0  = p[0];
  double s0  = p[1];
  double q1  = p[2];
  double s1  = p[3];
  double mu  = p[4];
  double eta = p[5];
  double nn  = p[6];

  double psum = 0, val = 0;

  for (int n = 0; ; n++) {
    if (n > 76) break;
    if (1 - psum <= kEpsilon) break;

    double poisson = TMath::Poisson(n, mu);
    psum += poisson;

    if (n == 0) {
      val += poisson * NormalizedGauss(x, q0, s0, gFitLo, gFitHi);
    } else {
      // FA: Gaussian peak
      val += poisson * (1 - eta)
           * NormalizedGauss(x, q0 + n * q1,
                             TMath::Sqrt(s0 * s0 + n * s1 * s1),
                             gFitLo, gFitHi);
      // PA: backscatter — Gauss(q0,s0) ⊗ Uniform(0, q1)
      val += poisson * eta
           * NormalizedBackscatter(x, q0, s0, q1, gFitLo, gFitHi);
    }
  }

  return nn * val;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void spe_dist(int nEvents = 100000, double mu = 0.2, double head = 10, double tail = 40) // ns
{
  auto pmt = new PMTSignal();
  pmt->SetModel(PMTSignal::MOYAL);
  pmt->SetTransitTime(200);
  pmt->SetTTS(1.2);
  pmt->SetMeanCharge(1.6);
  pmt->SetChargeRMS(0.4);
  pmt->SetBackscatterFraction(0.2);

  const int    kADCBits      = 12;
  const double kVpp          = 2500;  // mV
  const double kSamplingRate = 500;  // MHz
  const double kTermination  = 50;    // ohm

  auto fadc = new FADCWaveformGenerator();
  fadc->SetADCBits(kADCBits);
  fadc->SetVpp(kVpp);
  fadc->SetSamplingRate(kSamplingRate);
  fadc->SetRecordLength(500);
  fadc->SetPedestalMean(100);
  fadc->SetPedestalRMS(1.0);
  fadc->SetSignal(pmt);

  // ADC·bin → pC
  const double kBinWidth   = 1000.0 / kSamplingRate;
  const double kResolution = TMath::Power(2, kADCBits);
  const double adcToPc     = kVpp * kBinWidth / (kResolution * kTermination);

  auto hCharge = new TH1D("hCharge", "SPE Charge Distribution", 300, -0.5, 5.0);
  hCharge->GetXaxis()->SetTitle("Integrated charge (pC)");
  hCharge->GetYaxis()->SetTitle("Events");

  TRandom3 rng(0);

  for (int ev = 0; ev < nEvents; ev++) {
    pmt->Initialize();

    int npe = rng.Poisson(mu);
    for (int i = 0; i < npe; i++)
      pmt->AddHitTime(0);

    fadc->Digitize();

    const unsigned short *wf = fadc->GetWaveform();
    int ndp = fadc->GetNdp();

    // Pedestal from first 20 bins
    double ped = 0;
    for (int j = 0; j < 20; j++) ped += wf[j];
    ped /= 20.0;

    // Find peak bin
    int maxBin = 0;
    double maxVal = wf[0] - ped;
    for (int j = 1; j < ndp; j++) {
      double val = wf[j] - ped;
      if (val > maxVal) { maxVal = val; maxBin = j; }
    }

    // Integrate [maxBin - head, maxBin + tail]
    int headBins = static_cast<int>(head / kBinWidth);
    int tailBins = static_cast<int>(tail / kBinWidth);
    int start = std::max(0, maxBin - headBins);
    int stop  = std::min(ndp - 1, maxBin + tailBins);
    double charge = 0;
    for (int j = start; j <= stop; j++)
      charge += wf[j] - ped;

    hCharge->Fill(charge * adcToPc);
  }

  // Fit
  gFitLo = hCharge->GetXaxis()->GetXmin();
  gFitHi = hCharge->GetXaxis()->GetXmax();

  auto tf = new TF1("spefit", fitfunc, gFitLo, gFitHi, 7);
  tf->SetParNames("q0", "s0", "q1", "s1", "mu", "eta", "nn");
  tf->SetParameters(0,       // q0: pedestal mean
                    0.05,    // s0: pedestal sigma (pC)
                    1.6,     // q1: SPE mean charge (pC)
                    0.4,     // s1: SPE sigma (pC)
                    mu,      // mu: Poisson mean
                    0.2,     // eta: backscatter fraction
                    hCharge->Integral("width"));
  tf->SetParLimits(1, 0, 1);     // s0 > 0
  tf->SetParLimits(2, 0.5, 5);   // q1 physical range
  tf->SetParLimits(3, 0, 2);     // s1 > 0
  tf->SetParLimits(4, 0, 5);     // mu > 0
  tf->SetParLimits(5, 0, 1);     // 0 <= eta <= 1
  tf->SetParLimits(6, 0, 1e9);

  hCharge->Fit(tf, "R");
  hCharge->Draw();
  tf->Draw("same");
}
