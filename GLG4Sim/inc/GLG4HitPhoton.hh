#ifndef __GLG4HitPhoton_hh__
#define __GLG4HitPhoton_hh__
/** @file GLG4HitPhoton.hh
    Declares GLG4HitPhoton class and helper functions.
    
    This file is part of the GenericLAND software library.
    $Id: GLG4HitPhoton.hh,v 1.3 2014/04/11 01:41:02 yjko Exp $
    
    @author Glenn Horton-Smith
*/

#include <iostream>

/** GLG4HitPhoton stores information about a photon that makes a
    photoelectron in a PMT.  With count>1, it records multiple
    p.e. made at the same time at a PMT.

    The general contract for GLG4HitPhoton is as follows:
      - remember PMT ID, Time, position, KE, momentum, polarization, and count
      - provide Get/Set for all of the above, plus AddCount for count
      - initialize count=1, time & id = some invalid value

    This is almost the same "general contract" that was implemented
    for KLG4sim's KLHitPhoton by O. Tajima and G. Horton-Smith, but
    the code was rewritten for GLG4sim in December 2004.

    @author Glenn Horton-Smith
*/

class GLG4HitPhoton {
public:
  GLG4HitPhoton() { }

  void SetPMTID(int id) { fPMTID= id; }
  void SetTime(double t) { fTime= t; }
  void SetKineticEnergy(double KE);
  void SetWavelength(double wl);
  void SetPosition(double x, double y, double z);
  void SetMomentum(double x, double y, double z);
  void SetPolarization(double x, double y, double z);
  void SetCount(int count) { fCount= count; }
  void AddCount(int dcount) { fCount+= dcount; }
  void SetProcessTag(int proctag) { fProctag= proctag; } //EJ: 2007-11-06

  int GetPMTID() const { return fPMTID; }
  double GetTime() const { return fTime; }
  double GetKineticEnergy() const;
  double GetWavelength() const;
  template <class T> inline void GetPosition(T &x, T &y, T &z) const;
  template <class T> inline void GetMomentum(T &x, T &y, T &z) const;
  template <class T> inline void GetPolarization(T &x, T &y, T &z) const;
  int GetCount() const { return fCount; }
  int GetProcessTag() const { return fProctag; } //EJ: 2007-11-06
  void Print(std::ostream &) const;
  
private:
  double fTime;        /// time of hit 
  int fPMTID;          /// ID number of PMT the HitPhoton hit
  float fKE;           /// kinetic energy 
  float fPosition[3];  /// x,y,z components of position
  float fMomentum[3];  /// x,y,z components of momentum (normalized?)
  float fPolarization[3]; /// x,y,z components of polarization
  int fCount;          /// count of photons, often 1
  int fProctag;          /// process tag: 1=Cerenkov, 2=Scintillation, 3=Reemission: EJ:2007-11-06
};

template <class T> inline void 
GLG4HitPhoton::GetPosition(T &x, T &y, T &z) const {
  x= fPosition[0];
  y= fPosition[1];
  z= fPosition[2];
}

template <class T> inline void 
GLG4HitPhoton::GetMomentum(T &x, T &y, T &z) const {
  x= fMomentum[0];
  y= fMomentum[1];
  z= fMomentum[2];
}

template <class T> inline void 
GLG4HitPhoton::GetPolarization(T &x, T &y, T &z) const {
  x= fPolarization[0];
  y= fPolarization[1];
  z= fPolarization[2];
}


/** comparison function for sorting GLG4HitPhoton pointers
 */
inline bool
Compare_HitPhotonPtr_TimeAscending(const GLG4HitPhoton *a,
				   const GLG4HitPhoton *b)
{
  return a->GetTime() < b->GetTime();
}

#endif // __GLG4HitPhoton_hh__
