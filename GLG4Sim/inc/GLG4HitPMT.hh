/** @file GLG4HitPMT.hh
    Declares GLG4HitPMT class and helper functions.
    
    This file is part of the GenericLAND software library.
    $Id: GLG4HitPMT.hh,v 1.2 2013/11/11 02:39:14 jslee Exp $
    
    @author Glenn Horton-Smith, December 2004
*/
#ifndef GLG4HitPMT_hh
#define GLG4HitPMT_hh

#include <vector>
#include <cstddef>

#include "GLG4Sim/GLG4HitPhoton.hh"


/** GLG4HitPMT stores information about a PMT that detected one or more
    photoelectrons.

    The general contract for GLG4HitPMT is as follows:
      - remember ID and all photons that made photoelectrons in this PMT
      - provide Clear, DetectPhoton, SortTimeAscending, GetID, GetEntries,
        GetPhoton, and Print(ostream &) functions
      - take "ownership" of the photons registered by DetectPhoton,
        i.e., delete them when Clear or destrictor is called

    This is almost the same "general contract" that was implemented
    for KLG4sim's KLHitPMT by O. Tajima and G. Horton-Smith, but the
    code was rewritten for GLG4sim in December 2004.

    @author Glenn Horton-Smith
*/

class GLG4HitPMT {
public:
  GLG4HitPMT(int ID);
  ~GLG4HitPMT();

  void Clear();
  void DetectPhoton(GLG4HitPhoton*);
  void SortTimeAscending();

  int GetID() const { return fID; }
  int GetEntries() const { return fPhotons.size(); }
  GLG4HitPhoton* GetPhoton(int i) const { return fPhotons[i]; }

  void Print(std::ostream &, bool fullDetailsMode=false);

  static const size_t kApproxMaxIndividualHitPhotonsPerPMT;
  static const double kMergeTime;
  
private:
  int fID;
  std::vector<GLG4HitPhoton*> fPhotons;
};


/** comparison function for sorting GLG4HitPMT pointers
 */
inline bool Compare_HitPMTPtr_TimeAscending(const GLG4HitPMT * a, const GLG4HitPMT * b)
{
  // put empties at the end
  if(!a || a->GetEntries()<=0)
    return false;
  if(!b || b->GetEntries()<=0)
    return true;

  return a->GetPhoton(0)->GetTime() < b->GetPhoton(0)->GetTime();
}
#endif
