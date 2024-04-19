// This file is part of the GenericLAND software library.
// $Id: GLG4SteppingAction.hh,v 1.1.1.1 2013/11/08 05:33:05 jslee Exp $
//
#ifndef __GLG4SteppingAction_H__
#define __GLG4SteppingAction_H__ 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class GLG4PrimaryGeneratorAction;

class GLG4SteppingAction : public G4UserSteppingAction
{
public:
  GLG4SteppingAction();
  void UserSteppingAction(const G4Step* aStep);

private:
  GLG4PrimaryGeneratorAction* myGenerator;
};

#endif
