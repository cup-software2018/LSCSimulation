#ifndef MCScintData_HH
#define MCScintData_HH

#include "TClonesArray.h"

class MCScint;
class MCScintData : public TClonesArray {
public:
  MCScintData();
  MCScintData(const MCScintData & data);
  virtual ~MCScintData();

  virtual void Clear(const Option_t * opt = "");

  MCScint * Add();

  int GetN() const;
  MCScint * Get(int i) const;

  virtual void Print(const Option_t * opt = "") const;

private:
  int fN; //!

  ClassDef(MCScintData, 1)
};

inline int MCScintData::GetN() const { return GetEntriesFast(); }
inline MCScint * MCScintData::Get(int n) const { return (MCScint *)At(n); }

#endif
