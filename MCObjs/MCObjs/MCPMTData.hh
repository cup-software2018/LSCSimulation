#ifndef MCPMTData_HH
#define MCPMTData_HH

#include "TClonesArray.h"

class MCPMT;
class MCPMTData : public TClonesArray {
public:
  MCPMTData();
  MCPMTData(const MCPMTData & data);
  virtual ~MCPMTData();

  virtual void Clear(const Option_t * opt = "");

  MCPMT * Add();

  int GetN() const;
  MCPMT * Get(int i) const;

  virtual void Print(const Option_t * opt = "") const;

private:
  int fN; //!

  ClassDef(MCPMTData, 1)
};

inline int MCPMTData::GetN() const { return GetEntriesFast(); }
inline MCPMT * MCPMTData::Get(int n) const { return (MCPMT *)At(n); }

#endif
