#include "MCObjs/MCEventInfo.hh"

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
