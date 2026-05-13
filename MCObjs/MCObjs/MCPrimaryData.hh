#pragma once

#include "TClonesArray.h"

#include "MCObjs/MCPrimary.hh"
class MCPrimaryData : public TClonesArray {
public:
  MCPrimaryData();
  MCPrimaryData(const MCPrimaryData & data);
  virtual ~MCPrimaryData();

  void Clear(Option_t * opt = "") override;

  MCPrimary * Add();

  int GetN() const;
  MCPrimary * Get(int i) const;

  void Print(Option_t * opt = "") const override;

  ClassDef(MCPrimaryData, 1)
};

inline int MCPrimaryData::GetN() const { return GetEntriesFast(); }
inline MCPrimary * MCPrimaryData::Get(int n) const { return static_cast<MCPrimary *>(At(n)); }

