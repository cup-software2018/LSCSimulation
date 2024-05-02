#include <iostream>
#include "TString.h"
#include "MCObjs/MCEventInfo.hh"

using namespace std;

ClassImp(MCEventInfo)

MCEventInfo::MCEventInfo()
  : TObject()
{
  fEventNumber = 0;
}

MCEventInfo::MCEventInfo(const MCEventInfo & info)
  : TObject()
{
  fEventNumber = info.GetEventNumber();
}

MCEventInfo::~MCEventInfo() {}

void MCEventInfo::Print(const Option_t * opt)
{
  cout << Form("++++++++++++++++ Event No: %d ++++++++++++++++++", fEventNumber) << endl;
}
