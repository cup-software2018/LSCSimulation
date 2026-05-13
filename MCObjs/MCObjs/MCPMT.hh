#ifndef MCPMT_hh
#define MCPMT_hh

#include "TClonesArray.h"

#include "MCObjs/MCPhotonHit.hh"
class MCPMT : public TClonesArray {
public:
  MCPMT();
  MCPMT(int id);
  MCPMT(const MCPMT & pmt);
  virtual ~MCPMT();

  void Clear(Option_t * opt = "") override;

  void SetId(int id);
  int GetId() const;

  MCPhotonHit * AddHit();
  MCPhotonHit * AddHit(MCPhotonHit * hit);

  int GetNHit() const;
  MCPhotonHit * GetHit(int i) const;

  bool IsSortable() const override { return true; }
  int Compare(const TObject * object) const override;

  void Print(Option_t * opt = "") const override;

private:
  int fPMTId = -1;

  ClassDef(MCPMT, 1);
};

inline void MCPMT::SetId(int Id) { fPMTId = Id; }
inline int MCPMT::GetId() const { return fPMTId; }

inline int MCPMT::GetNHit() const { return GetEntriesFast(); }
inline MCPhotonHit * MCPMT::GetHit(int n) const { return static_cast<MCPhotonHit *>(At(n)); }


#endif
