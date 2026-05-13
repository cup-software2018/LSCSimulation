#ifndef MCScintData_HH
#define MCScintData_HH

#include "TClonesArray.h"

#include "MCObjs/MCScint.hh"
class MCScintData : public TClonesArray {
public:
  MCScintData();
  MCScintData(const MCScintData & data);
  virtual ~MCScintData();

  void Clear(Option_t * opt = "") override;

  MCScint * Add();
  MCScint * Add(int id);

  int GetN() const;
  MCScint * Get(int i) const;
  MCScint * FindScint(int id);

  void Print(Option_t * opt = "") const override;

  ClassDef(MCScintData, 1)
};

inline int MCScintData::GetN() const { return GetEntriesFast(); }
inline MCScint * MCScintData::Get(int n) const { return static_cast<MCScint *>(At(n)); }

#endif
