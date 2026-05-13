#ifndef MCTrack_hh
#define MCTrack_hh

#include "TClonesArray.h"
#include "TString.h"

#include "MCObjs/MCStep.hh"
class MCTrack : public TClonesArray {
public:
  MCTrack();
  MCTrack(const MCTrack & trk);
  virtual ~MCTrack();

  void Clear(Option_t * opt = "") override;

  void SetParticleName(const char * name);
  void SetPDGCode(int code);
  void SetTrackId(int id);
  void SetParentId(int id);
  void SetVertex(double x, double y, double z);
  void SetKineticEnergy(double ke);
  void SetGlobalTime(double time);
  void SetLocalTime(double time);
  void SetProcessName(const char * name);

  const char * GetParticleName() const;
  int GetPDGCode() const;
  int GetTrackId() const;
  int GetParentId() const;
  void GetVertex(double & x, double & y, double & z) const;
  double GetKineticEnergy() const;
  double GetGlobalTime() const;
  double GetLocalTime() const;
  const char * GetProcessName() const;

  MCStep * AddStep();

  int GetNStep() const;
  MCStep * GetStep(int i) const;

  void Print(Option_t * opt = "") const override;

private:
  int fPDGCode = 0;
  int fTrackId = 0;
  int fParentId = 0;

  double fVx = 0, fVy = 0, fVz = 0;
  double fKineticEnergy = 0;
  double fGlobalTime = 0;
  double fLocalTime = 0;

  TString fParticleName;
  TString fProcessName;

  ClassDef(MCTrack, 1)
};
inline void MCTrack::SetParticleName(const char * name)
{
  fParticleName = name;
}

inline void MCTrack::SetPDGCode(int code) { fPDGCode = code; }

inline void MCTrack::SetTrackId(int id) { fTrackId = id; }

inline void MCTrack::SetParentId(int id) { fParentId = id; }

inline void MCTrack::SetVertex(double x, double y, double z)
{
  fVx = x;
  fVy = y;
  fVz = z;
}

inline void MCTrack::SetKineticEnergy(double ke) { fKineticEnergy = ke; }

inline void MCTrack::SetGlobalTime(double time) { fGlobalTime = time; }

inline void MCTrack::SetLocalTime(double time) { fLocalTime = time; }

inline void MCTrack::SetProcessName(const char * name) { fProcessName = name; }

inline const char * MCTrack::GetParticleName() const
{
  return fParticleName.Data();
}

inline int MCTrack::GetPDGCode() const { return fPDGCode; }

inline int MCTrack::GetTrackId() const { return fTrackId; }

inline int MCTrack::GetParentId() const { return fParentId; }

inline void MCTrack::GetVertex(double & x, double & y, double & z) const
{
  x = fVx;
  y = fVy;
  z = fVz;
}

inline double MCTrack::GetKineticEnergy() const { return fKineticEnergy; }

inline double MCTrack::GetGlobalTime() const { return fGlobalTime; }

inline double MCTrack::GetLocalTime() const { return fLocalTime; }

inline const char * MCTrack::GetProcessName() const
{
  return fProcessName.Data();
}

inline int MCTrack::GetNStep() const { return GetEntriesFast(); }
inline MCStep * MCTrack::GetStep(int n) const { return static_cast<MCStep *>(At(n)); }

#endif
