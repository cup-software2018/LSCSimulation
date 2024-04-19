#include "MCObjs/MCStep.hh"

ClassImp(MCStep)

MCStep::MCStep()
    : TObject()
{
  fStepLength = 0;
  fEnergyDeposit = 0;
  fEnergyDepositNonIonizing = 0;
  fGlobalTime = 0;
  fLocalTime = 0;
  fX = 0;
  fY = 0;
  fZ = 0;
}

MCStep::MCStep(const MCStep & step)
    : TObject()
{
  fStepLength = step.GetStepLength();
  fEnergyDeposit = step.GetEnergyDeposit();
  fEnergyDepositNonIonizing = step.GetEnergyDepositNonIonizing();
  fGlobalTime = step.GetGlobalTime();
  fLocalTime = step.GetLocalTime();
  step.GetStepPoint(fX, fY, fZ);
  fVolumeName = step.GetVolumeName();
}

MCStep::~MCStep() {}
