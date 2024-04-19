#ifndef MCPrimaryData_HH
#define MCPrimaryData_HH

#include "TClonesArray.h"

class MCPrimary;
class MCPrimaryData : public TClonesArray {
public:
  MCPrimaryData();
  MCPrimaryData(const MCPrimaryData & data);
  virtual ~MCPrimaryData();

  virtual void Clear(const Option_t * opt = "");

  MCPrimary * Add();

  int GetN() const;
  MCPrimary * Get(int i) const;

  virtual void Print(const Option_t * opt = "") const;

private:
  int fN; //!

  ClassDef(MCPrimaryData, 1)
};

inline int MCPrimaryData::GetN() const { return GetEntriesFast(); }
inline MCPrimary * MCPrimaryData::Get(int n) const { return (MCPrimary *)At(n); }

#endif
