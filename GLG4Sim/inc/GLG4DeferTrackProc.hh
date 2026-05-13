// This file is part of the GenericLAND software library.
// $Id: GLG4DeferTrackProc.hh,v 1.2 2013/11/11 01:22:22 jslee Exp $
//
// Process to limit step length to stay within event time
// and defer long tracks (AND tracks which start after event time) using
// defered particle "generator".
//
// Written: G. Horton-Smith, 29-Oct-2001
//
#ifndef GLG4DeferTrackProc_hh
#define GLG4DeferTrackProc_hh

#include "globals.hh"
#include "G4VProcess.hh"


class GLG4PrimaryGeneratorAction;
class G4HadronCaptureProcess;
class G4UImessenger; // for G4ProcessTable.hh

class GLG4DeferTrackProc : public G4VProcess 
{
 public:  //with description     
  GLG4DeferTrackProc(const G4String & processName = "DeferTrackProc");
  ~GLG4DeferTrackProc();

  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track & track,
							G4double previousStepSize,
							G4ForceCondition * condition);

  virtual G4VParticleChange * PostStepDoIt(const G4Track&, const G4Step&);
  //  no operation in  AtRestDoIt      
  virtual G4VParticleChange * AtRestDoIt(const G4Track&, const G4Step&) { return NULL; };
  //  no operation in  AlongStepDoIt
  virtual G4VParticleChange * AlongStepDoIt(const G4Track&, const G4Step&) { return NULL; };			    
  //  no operation in  AtRestGPIL
  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition*) { return -1.0; };
  //  no operation in  AlongStepGPIL
  virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
							 G4double  ,
							 G4double  ,
							 G4double& ,
							 G4GPILSelection*) { return -1.0; };

private:
  // hide assignment operator as private 
  GLG4DeferTrackProc(GLG4DeferTrackProc&);
  GLG4DeferTrackProc& operator=(const GLG4DeferTrackProc& right);
  
private:
  GLG4PrimaryGeneratorAction * _generator;
  G4double kCarTolerance;
};
#endif

