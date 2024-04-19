#ifndef MCScintStep_hh
#define MCScintStep_hh

#include "MCObjs/MCStep.hh"

class MCScintStep : public MCStep {
public:
  MCScintStep();
  MCScintStep(const MCScintStep & step);
  virtual ~MCScintStep();

  inline void SetEnergyVisible(double val) { fEnergyVisible = val; }
  inline void SetNScintPhoton(int val) { fNScintPhoton = val; }

  inline float GetEnergyVisible() const { return fEnergyVisible; }
  inline int GetNScintPhoton() const { return fNScintPhoton; }

private:
  int fNScintPhoton;
  float fEnergyVisible;

  ClassDef(MCScintStep, 1)
};
#endif
