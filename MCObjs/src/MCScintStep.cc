#include "MCScintStep.hh"

ClassImp(MCScintStep)

MCScintStep::MCScintStep()
  : MCStep()
{
}

MCScintStep::MCScintStep(const MCScintStep & step)
  : MCStep(step),
    fNScintPhoton(step.GetNScintPhoton()),
    fEnergyVisible(step.GetEnergyVisible())
{
}

MCScintStep::~MCScintStep() = default;
