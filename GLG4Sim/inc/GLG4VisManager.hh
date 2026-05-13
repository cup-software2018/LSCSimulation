// This file is part of the GenericLAND software library.
// $Id: GLG4VisManager.hh,v 1.1.1.1 2013/11/08 05:33:05 jslee Exp $
//
// GenericLAND visualization manager based on Geant4's "MyVisManager"
//   -- main purpose of defining our own is to reorient "up" vector
//
// See class description of G4VisManager for more details.
//
// Author:  Glenn Horton-Smith, Jan 28, 2000
//
#ifdef G4VIS_USE
#ifndef GLG4VISMANAGER_HH
#define GLG4VISMANAGER_HH

#include "G4VisManager.hh"

class GLG4VisManager: public G4VisManager {

public: // With description

  GLG4VisManager ();

private:

  void RegisterGraphicsSystems ();

};

#endif
#endif
