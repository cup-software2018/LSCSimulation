#ifndef MCPrimary_hh
#define MCPrimary_hh

#include <iostream>

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

class MCPrimary : public TObject {
public:
  MCPrimary();
  MCPrimary(const MCPrimary & prim);
  virtual ~MCPrimary();

  void SetParticleName(const char * name);
  void SetPDGCode(int val);
  void SetVertex(double x, double y, double z);
  void SetMomentum(double x, double y, double z);
  void SetKineticEnergy(double val);
  void SetT0(double val);
  void SetTrackId(int val);

  TString GetParticleName() const;
  int GetPDGCode() const;
  void GetVertex(double & x, double & y, double & z) const;
  void GetMomentum(double & x, double & y, double & z) const;
  double GetKineticEnergy() const;
  double GetT0() const;
  int GetTrackId() const;

  virtual void Print(const Option_t * opt = "") const;

private:
  TString fName;
  int fPDGCode;
  int fTrackId;
  double fVx, fVy, fVz;
  double fPx, fPy, fPz;
  double fKineticEnergy;
  double fT0;

  ClassDef(MCPrimary, 1)
};

inline void MCPrimary::SetParticleName(const char * name) { fName = name; }

inline void MCPrimary::SetPDGCode(int val) { fPDGCode = val; }

inline void MCPrimary::SetVertex(double x, double y, double z)
{
  fVx = x;
  fVy = y;
  fVz = z;  
}

inline void MCPrimary::SetMomentum(double x, double y, double z)
{
  fPx = x;
  fPy = y;
  fPz = z;  
}

inline void MCPrimary::SetKineticEnergy(double val) { fKineticEnergy = val; }

inline void MCPrimary::SetT0(double val) { fT0 = val; }

inline void MCPrimary::SetTrackId(int val) { fTrackId = val; }

inline TString MCPrimary::GetParticleName() const { return fName; }

inline int MCPrimary::GetPDGCode() const { return fPDGCode; }

inline void MCPrimary::GetVertex(double & x, double & y, double & z) const 
{ 
  x = fVx;
  y = fVy;
  z = fVz;
}

inline void MCPrimary::GetMomentum(double & x, double & y, double & z) const 
{ 
  x = fPx;
  y = fPy;
  z = fPz;
}

inline double MCPrimary::GetKineticEnergy() const { return fKineticEnergy; }

inline double MCPrimary::GetT0() const { return fT0; }

inline int MCPrimary::GetTrackId() const { return fTrackId; }
#endif
