#ifndef MCPMT_hh
#define MCPMT_hh

#include "TClonesArray.h"

class MCPhotonHit;
class MCPMT : public TClonesArray {
public:
  MCPMT();
  MCPMT(int id);
  MCPMT(const MCPMT & pmt);
  virtual ~MCPMT();

  virtual void Clear(const Option_t * opt = "");

  void SetId(int id);
  int GetId() const;

  MCPhotonHit * AddHit();
  MCPhotonHit * AddHit(MCPhotonHit * hit);

  int GetNHit() const;
  MCPhotonHit * GetHit(int i) const;

  virtual bool IsSortable() const { return true; }
  virtual int Compare(const TObject * object) const;

  virtual void Print(const Option_t * opt = "") const;

private:
  int fPMTId;
  int fNHit; //!

  ClassDef(MCPMT, 1);
};

inline void MCPMT::SetId(int Id) { fPMTId = Id; }
inline int MCPMT::GetId() const { return fPMTId; }

inline int MCPMT::GetNHit() const { return GetEntriesFast(); }
inline MCPhotonHit * MCPMT::GetHit(int n) const { return (MCPhotonHit *)At(n); }


#endif
