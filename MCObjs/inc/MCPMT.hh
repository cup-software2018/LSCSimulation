#pragma once

#include <vector>

#include "TObject.h"

#include "MCPhotonHit.hh"

class MCPMT : public TObject {
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

  void Sort();

  bool IsSortable() const override { return true; }
  int Compare(const TObject * object) const override;

  void Print(Option_t * opt = "") const override;

private:
  int fPMTId = -1;
  std::vector<MCPhotonHit> fHits;

  ClassDef(MCPMT, 2);
};

inline void MCPMT::SetId(int id) { fPMTId = id; }
inline int MCPMT::GetId() const { return fPMTId; }

inline int MCPMT::GetNHit() const { return (int)fHits.size(); }
inline MCPhotonHit * MCPMT::GetHit(int n) const { return const_cast<MCPhotonHit *>(&fHits[n]); }
