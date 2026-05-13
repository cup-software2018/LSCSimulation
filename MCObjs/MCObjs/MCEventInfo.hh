#ifndef MCEventInfo_hh
#define MCEventInfo_hh

#include "TObject.h"

class MCEventInfo : public TObject {
public:
  MCEventInfo();
  MCEventInfo(const MCEventInfo & info);
  virtual ~MCEventInfo();

  void SetEventNumber(unsigned int n);
  unsigned int GetEventNumber() const;

  void Print(Option_t * opt = "") const override;

private:
  unsigned int fEventNumber = 0;

  ClassDef(MCEventInfo, 1)
};

inline void MCEventInfo::SetEventNumber(unsigned int n) { fEventNumber = n; }
inline unsigned int MCEventInfo::GetEventNumber() const { return fEventNumber; }

#endif
