#ifndef AbsSignal_hh
#define AbsSignal_hh

#include "TObject.h"
#include "TRandom3.h"
#include "TGraph.h"

class AbsSignal : public TObject {
public:
  AbsSignal();
  virtual ~AbsSignal();

  virtual void AddHitTime(double val) = 0;
  virtual void Initialize() = 0;
  virtual void Prepare() = 0;
  virtual double EvalPulseHeight(double t) = 0;

  virtual double GetTimeMinimum() const;
  virtual double GetTimeMaximum() const;

  virtual void Draw(Option_t * option = "");  

protected:
  bool fIsPrepared;

  double fMinimumTime;
  double fMaximumTime;

  TGraph * fGraph;
  TRandom3 * fRandom;  

  ClassDef(AbsSignal, 0)
};

inline double AbsSignal::GetTimeMinimum() const { return fMinimumTime; }
inline double AbsSignal::GetTimeMaximum() const { return fMaximumTime; }

#endif
