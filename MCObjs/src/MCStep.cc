#include "MCStep.hh"

ClassImp(MCStep)

MCStep::MCStep()
  : TObject()
{
}

MCStep::MCStep(const MCStep & step)
  : TObject(),
    fStepLength(step.GetStepLength()),
    fEnergyDeposit(step.GetEnergyDeposit()),
    fEnergyDepositNonIonizing(step.GetEnergyDepositNonIonizing()),
    fGlobalTime(step.GetGlobalTime()),
    fLocalTime(step.GetLocalTime()),
    fVolumeName(step.GetVolumeName())
{
  step.GetStepPoint(fX, fY, fZ);
}
