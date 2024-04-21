#include "TrgSim/AbsSignal.hh"

ClassImp(AbsSignal)

AbsSignal::AbsSignal()
    : TObject()
{
  fIsPrepared = false;

  fRandom = new TRandom3();
  fRandom->SetSeed(0);

  fGraph = nullptr;
}

AbsSignal::~AbsSignal()
{
  delete fRandom;
  if (fGraph) delete fGraph;
}

void AbsSignal::Draw(Option_t * option)
{
  if (!fIsPrepared) Prepare();

  double timei = GetTimeMinimum();
  double timef = GetTimeMaximum();

  double dtime = (timef - timei) / 10000;

  if (fGraph) delete fGraph;

  fGraph = new TGraph();

  int np = 0;
  for (double t = timei; t < timef; t += dtime) {
    double height = EvalPulseHeight(t);
    fGraph->SetPoint(np, t, height);
    np += 1;
  }

  TString opt(option);
  if (opt.IsNull()) opt = "AL";

  fGraph->Draw(opt.Data());
}
