/** @file GLG4HitPhoton.cc
    For GLG4HitPhoton class.

    This file is part of the GenericLAND software library.
    $Id: GLG4HitPhoton.cc,v 1.4 2014/04/11 01:41:02 yjko Exp $

    @author Glenn Horton-Smith, December 2004
*/

#include "GLG4Sim/GLG4HitPhoton.hh"
#include <CLHEP/Units/PhysicalConstants.h>
//#include <CLHEP/Units/GlobalPhysicalConstants.h> //EJ

using namespace CLHEP;

/// set kinetic energy and wavelength of photon.
void GLG4HitPhoton::SetKineticEnergy(double KE) { fKE = KE; }

/// set wavelength and kinetic energy of photon
void GLG4HitPhoton::SetWavelength(double wl) { fKE = 2 * pi * hbarc / wl; }

void GLG4HitPhoton::SetPosition(double x, double y, double z)
{
  fPosition[0] = x;
  fPosition[1] = y;
  fPosition[2] = z;
}

void GLG4HitPhoton::SetMomentum(double x, double y, double z)
{
  fMomentum[0] = x;
  fMomentum[1] = y;
  fMomentum[2] = z;
}

void GLG4HitPhoton::SetPolarization(double x, double y, double z)
{
  fPolarization[0] = x;
  fPolarization[1] = y;
  fPolarization[2] = z;
}

double GLG4HitPhoton::GetKineticEnergy() const { return fKE; }

double GLG4HitPhoton::GetWavelength() const { return 2 * pi * hbarc / fKE; }
