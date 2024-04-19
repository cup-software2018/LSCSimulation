#include "MCObjs/MCScintStep.hh"

ClassImp(MCScintStep)

MCScintStep::MCScintStep()
    : MCStep()
{
  fEnergyVisible = 0;
  fNScintPhoton = 0;
}

MCScintStep::MCScintStep(const MCScintStep & step)
    : MCStep(step)
{
  fEnergyVisible = step.GetEnergyVisible();
  fNScintPhoton = step.GetNScintPhoton();
}

MCScintStep::~MCScintStep() {}
