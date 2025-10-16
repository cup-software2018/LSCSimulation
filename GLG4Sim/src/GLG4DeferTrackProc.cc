// This file is part of the GenericLAND software library.
// $Id: GLG4DeferTrackProc.cc,v 1.3 2013/11/11 01:22:25 jslee Exp $
//
// Process to limit step length to stay within event time
// and defer long tracks (AND tracks which start after event time) using
// defered particle "generator".
//
// Written: G. Horton-Smith, 29-Oct-2001
//
////////////////////////////////////////////////////////////////////////

#include "GLG4Sim/GLG4DeferTrackProc.hh"
#include "GLG4Sim/GLG4PrimaryGeneratorAction.hh"

#include "G4ios.hh"
#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"
#include "G4ProcessTable.hh"
#include "G4GeometryTolerance.hh" // for kCarTolerance

#include "TString.h"

#include "LSCSim/LSCUserTrackInformation.hh"

using namespace CLHEP;

GLG4DeferTrackProc::GLG4DeferTrackProc(const G4String & aName)
  : G4VProcess(aName)
{
  if(verboseLevel > 0){
    G4cout<<GetProcessName()<<" is created "<< G4endl;
  }
  
  _generator= GLG4PrimaryGeneratorAction::GetTheGLG4PrimaryGeneratorAction();
  
  if(_generator == 0){
    G4String msg = "GLG4DeferTrackProc:: no GLG4PrimaryGeneratorAction instance.";
    G4Exception("GLG4DeferTrackProc::GLG4DeferTrackProc", "", FatalException, msg);
  }

  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

GLG4DeferTrackProc::~GLG4DeferTrackProc()
{}

GLG4DeferTrackProc::GLG4DeferTrackProc(GLG4DeferTrackProc & right)
  : G4VProcess(right)
{}

G4double GLG4DeferTrackProc::PostStepGetPhysicalInteractionLength(const G4Track & aTrack, G4double,
								  G4ForceCondition * condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  // apply maximum time limit
  G4double dTime= (_generator->GetEventWindow() - aTrack.GetGlobalTime());
  if(dTime <= 0.0){
    return kCarTolerance;
  }

  G4double beta = aTrack.GetDynamicParticle()->GetTotalMomentum()/aTrack.GetTotalEnergy();
  return beta*c_light*dTime;
}

G4VParticleChange * GLG4DeferTrackProc::PostStepDoIt(const G4Track & aTrack, const G4Step&)
{
  _generator->DeferTrackToLaterEvent(&aTrack);

  LSCUserTrackInformation * info = (LSCUserTrackInformation*)aTrack.GetUserInformation();
  info->SetTrackStatusFlags(LSCTrackStatus::deferred);

  aParticleChange.Initialize(aTrack);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  return &aParticleChange;
}
