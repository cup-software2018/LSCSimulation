#ifndef MCPMTData_HH
#define MCPMTData_HH

#include "TClonesArray.h"

#include "MCObjs/MCPMT.hh"
class MCPMTData : public TClonesArray {
public:
  MCPMTData();
  MCPMTData(const MCPMTData & data);
  virtual ~MCPMTData();

  void Clear(Option_t * opt = "") override;

  MCPMT * Add();

  int GetN() const;
  MCPMT * Get(int i) const;

  void Print(Option_t * opt = "") const override;

  ClassDef(MCPMTData, 1)
};

inline int MCPMTData::GetN() const { return GetEntriesFast(); }
inline MCPMT * MCPMTData::Get(int n) const { return static_cast<MCPMT *>(At(n)); }

#endif
