// This file is part of the GenericLAND software library.
// $Id: GLG4OpAttenuation.hh,v 1.2 2013/11/14 02:24:31 jslee Exp $
//
//  "Attenuation" (absorption or scattering) of optical photons
//
//   GenericLAND Simulation
//
//   Original: Glenn Horton-Smith, Dec 2001
//
// GLG4OpAttenuation.hh
// 

#ifndef GLG4OpAttenuation_hh
#define GLG4OpAttenuation_hh

#include "G4OpAbsorption.hh"

class GLG4OpAttenuation : public G4OpAbsorption
{
public: // Without description
  GLG4OpAttenuation(const G4String& processName = "Attenuation");
  ~GLG4OpAttenuation();

public: // With description
  // This is the method implementing attenuation of optical 
  // photons.  Fraction of photons scattered or absorbed is
  // determined by the MaterialProperyVector "OPSCATFRAC".

  G4VParticleChange * PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
};
#endif
