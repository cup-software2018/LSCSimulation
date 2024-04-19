// This file is part of the GenericLAND software library.
// $Id: GLG4SteppingAction.cc,v 1.2 2013/11/09 23:48:54 jslee Exp $
//
//
//  GenericLAND Simulation
//
//  Concrete implementation of G4UserSteppingAction
//
//  Current uses:
//    * Measure inter-step CPU time, broken down by process and particle type
//
//  Anticipated uses:
//    * Find PMT _fast_ when entering outer buffer
//
//  Author: Glenn Horton-Smith, April 7, 2000

#include "GLG4Sim/GLG4SteppingAction.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4OpticalPhoton.hh"
#include "GLG4Sim/GLG4Scint.hh"
#include "G4SteppingManager.hh"
#include "G4VisExtent.hh"
#include "GLG4Sim/GLG4PrimaryGeneratorAction.hh"
#include "G4VSolid.hh"

GLG4SteppingAction::GLG4SteppingAction()
{
   myGenerator= GLG4PrimaryGeneratorAction::GetTheGLG4PrimaryGeneratorAction();
   if (myGenerator == 0) {
     G4Exception("GLG4SteppingAction:: no GLG4PrimaryGeneratorAction instance.", "", FatalException,"");
   }
}

#ifdef G4DEBUG
#include "G4Timer.hh"
#include "map"

struct GLG4SteppingAction_time_s {
  G4double sumtime;
  G4int stepcount;
  GLG4SteppingAction_time_s() { sumtime=0.0; stepcount=0; }
  void sum(G4double time) { sumtime+= time; stepcount++; }
};

typedef G4std::map<G4String,GLG4SteppingAction_time_s> GLG4SteppingAction_time_map;
typedef G4std::map<G4String,GLG4SteppingAction_time_s>::iterator GLG4SteppingAction_time_map_iterator;

GLG4SteppingAction_time_map GLG4SteppingAction_times;
G4double GLG4SteppingAction_internal_time = 0.0;

int
GLG4SteppingAction_dump_times(void)
{
  for (GLG4SteppingAction_time_map_iterator i= GLG4SteppingAction_times.begin();
       i != GLG4SteppingAction_times.end();
       i++)
    {
      G4cout << i->first << ' '
	     << i->second.sumtime << ' '
	     << i->second.stepcount << G4endl;
      G4cout.flush();
    }
  return 1;
}

#include "G4VisAttributes.hh"
#include "G4Color.hh"

const int nrowIlluminationMap= 600;
const int ncolIlluminationMap= 600;
G4double widthIlluminationMap= -1.0;
G4double IlluminationMap[6][nrowIlluminationMap][ncolIlluminationMap][3];

static void
updateIlluminationMap(int imap, G4double srow, G4double scol, const G4Color *cp)
{
  G4double *mp= IlluminationMap[imap][(int)(srow*nrowIlluminationMap)]
    [(int)(scol*ncolIlluminationMap)];
  mp[0] += cp->GetRed();
  mp[1] += cp->GetGreen();
  mp[2] += cp->GetBlue();
}
		      
int
GLG4SteppingAction_dump_IlluminationMap(void)
{
  G4double maxval= 0.0, meanval=0.0;
  G4int nsum= 0;
  {for (int kmap=0; kmap<6; kmap++)
    for (int irow=0; irow< nrowIlluminationMap; irow++)
      for (int jcol=0; jcol< ncolIlluminationMap; jcol++)
	for (int kcol=0; kcol<3; kcol++)
	  {
	    G4double v= IlluminationMap[kmap][irow][jcol][kcol];
	    if (v > 0.0)
	      {
		if (v > maxval)
		  maxval= v;
		meanval += v;
		nsum++;
	      }
	  }
  }
  if (maxval == 0.0)
    {
      G4cout << "Empty Illumination Map" << G4endl;
      return 0;
    }
  meanval /= nsum;
  {for (int kmap=0; kmap<6; kmap++)
    {
      static char filename[]="map#.ppm";
      filename[3]= kmap+'0';
      G4std::ofstream of(filename);
      of << "P6\n# Illumination map " << kmap << "\n"
	     << nrowIlluminationMap << ' ' << nrowIlluminationMap << " 255\n";
      for (int irow=0; irow< nrowIlluminationMap; irow++)
	{
	  for (int jcol=0; jcol< ncolIlluminationMap; jcol++)
	    for (int kcol=0; kcol<3; kcol++)
	      {
		G4double v=IlluminationMap[kmap][irow][jcol][kcol];
		unsigned char byte;
		if (v < meanval)
		  byte= (unsigned char)(128*v/meanval);
		else
		  byte= (unsigned char)(128.0
					+ 127.99*(v-meanval)/(maxval-meanval));
		of.put(byte);
	      }
	}
      of.close();
    }
  }
  return 1;
}

G4double GLG4SteppingAction_totEdep = 0.0;
G4int  GLG4SteppingAction_MaxStepNumber= 100000;
#endif /* G4DEBUG */

void
GLG4SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* track= aStep->GetTrack();

#ifdef G4DEBUG
  static G4Timer timer;
  static G4int num_zero_steps_in_a_row=0;

  timer.Stop();
  G4double dut= timer.GetUserElapsed();
  timer.Start();

  // check for NULL world volume
  if (track->GetVolume() == NULL) {
    const G4VProcess* lastproc= track->GetStep()->GetPostStepPoint()
      ->GetProcessDefinedStep();
    G4cerr << "GLG4SteppingAction: Track in NULL volume, terminating!\n"
	   << " step_no=" << track->GetCurrentStepNumber()
	   << " type=" << track->GetDefinition()->GetParticleName()
	   << "\n volume=NULL"
	   << " last_process="
	   << (lastproc!=0 ? lastproc->GetProcessName() : G4String("NULL"))
	   << "\n position=" << track->GetPosition()
	   << " momentum=" << track->GetMomentum()
	   << G4endl;
    track->SetTrackStatus( fStopAndKill );    
  }

  // check for very high number of steps
  if (track->GetCurrentStepNumber() > GLG4SteppingAction_MaxStepNumber) {
    const G4VPhysicalVolume* pv= track->GetVolume();
    const G4VProcess* lastproc= track->GetStep()->GetPostStepPoint()
      ->GetProcessDefinedStep();
    G4cerr << "GLG4SteppingAction: Too many steps for this track, terminating!\n"
	   << " step_no=" << track->GetCurrentStepNumber()
	   << " type=" << track->GetDefinition()->GetParticleName()
	   << "\n volume=" << (pv!=0 ? pv->GetName() : G4String("NULL"))
	   << " last_process="
	   << (lastproc!=0 ? lastproc->GetProcessName() : G4String("NULL"))
	   << "\n position=" << track->GetPosition()
	   << " momentum=" << track->GetMomentum()
	   << G4endl;
    track->SetTrackStatus( fStopAndKill );    
  }

  // check for too many zero steps in a row
  if (aStep->GetStepLength() <= 0.0 && track->GetCurrentStepNumber() > 1) {
    ++num_zero_steps_in_a_row;
    if (num_zero_steps_in_a_row >= 4) {
      G4cerr << "GLG4SteppingAction: Too many zero steps for this track, terminating!"
	     << G4endl;
      track->SetTrackStatus( fStopAndKill );
      num_zero_steps_in_a_row= 0;
    }
  }
  else
    num_zero_steps_in_a_row= 0;

  // check total energy deposit
  GLG4SteppingAction_totEdep += aStep->GetTotalEnergyDeposit();

  // debugging (timing) info
  static G4String lastParticleName;
  if (dut > 0.0) {
    G4String key;
    G4String particleName= track->GetDefinition()->GetParticleName();
    const G4VProcess
      *preProcess=aStep->GetPreStepPoint()->GetProcessDefinedStep(),
      *postProcess= aStep->GetPostStepPoint()->GetProcessDefinedStep();
    if (preProcess)
      key += preProcess->GetProcessName();
    key += "_";
    if (postProcess)
      key += postProcess->GetProcessName();
    GLG4SteppingAction_times[key].sum(dut);
    GLG4SteppingAction_times[particleName].sum(dut);
    key += "_" + particleName;
    if (particleName != lastParticleName) {
      key += "_" + lastParticleName;
      lastParticleName= particleName;
    }
    GLG4SteppingAction_times[key].sum(dut);
  }
#endif

  // do scintillation photons, and also re-emission
    {
      // invoke scintillation process
      G4VParticleChange * pParticleChange
	= GLG4Scint::GenericPostPostStepDoIt(aStep);
      // were any secondaries defined?
      G4int iSecondary= pParticleChange->GetNumberOfSecondaries();
      if (iSecondary > 0)
	{
	  // add secondaries to the list
	  while ( (iSecondary--) > 0 )
	    {
	      G4Track * tempSecondaryTrack
		= pParticleChange->GetSecondary(iSecondary);
	      fpSteppingManager->GetfSecondary()
		->push_back( tempSecondaryTrack );
	    }
	}
      // clear ParticleChange
      pParticleChange->Clear();
    }

  // Commented out because this duplicates the function of GLG4DeferTrackProc
  // // if end step time since start of event is past event window, defer to later
  // G4TrackStatus status= track->GetTrackStatus();
  // if ( (status==fAlive || status==fSuspend) &&
  //      aStep->GetPostStepPoint()->GetGlobalTime() >
  //      myGenerator->GetEventWindow() ) {
  //   myGenerator->DeferTrackToLaterEvent(track);
  //   track->SetTrackStatus(status=fStopAndKill);
  // }
  
#ifdef G4DEBUG
  // Particle illumination map
  G4TrackStatus status= track->GetTrackStatus();
  if ( status == fStopAndKill  ||  status == fKillTrackAndSecondaries )
    {
      G4double sx, sy, sz;
      if (widthIlluminationMap <= 0.0) {
	const G4VTouchable * t= track->GetTouchable();
	int depth= t->GetHistoryDepth();
	G4VSolid* s= t->GetSolid( depth ); // this should be the world volume
	G4VisExtent vx( s->GetExtent() );
	sx= fabs(vx.GetXmin())+fabs(vx.GetXmax());
	sy= fabs(vx.GetYmin())+fabs(vx.GetYmax());
	sz= fabs(vx.GetZmin())+fabs(vx.GetZmax());
	// take the middle value
	if (sx <= sy && sy <= sz) widthIlluminationMap= sy;
	else if (sy <= sx && sx <= sz) widthIlluminationMap= sx;
	else if (sx <= sz && sz <= sy) widthIlluminationMap= sz;
	else G4Exception("GLG4SteppingAction: oh dear, transitivity doesn't work.");
      }
      G4ThreeVector postPos= aStep->GetPostStepPoint()->GetPosition();
      sx= 0.5 + postPos.x() / widthIlluminationMap;
      sy= 0.5 + postPos.y() / widthIlluminationMap;
      sz= 0.5 + postPos.z() / widthIlluminationMap;
      if (sx > 0.0 && sy > 0.0 && sz > 0.0
	  && sx < 1.0 && sy < 1.0 && sz < 1.0)
	{
	  static const G4Color defaultcolor(0.1,0.1,0.1);
	  G4VPhysicalVolume *pv;
	  const G4VisAttributes *att;
	  const G4Color *c;
	  ( ((pv=track->GetNextVolume())||(pv=track->GetVolume()))
	    && (att=pv->GetLogicalVolume()->GetVisAttributes())
	    && (c=&(att->GetColor()))
	    ) || (c=&defaultcolor);
	  updateIlluminationMap( ((sx<0.5) ? 0 : 1),  1.0-sz,  sy,  c);
	  updateIlluminationMap( ((sy<0.5) ? 2 : 3),  1.0-sz,  sx,  c);
	  updateIlluminationMap( ((sz<0.5) ? 4 : 5),  1.0-sy,  sx,  c);
	}
    }
  // reset timer; measure our own elapsed time
  timer.Stop();
  GLG4SteppingAction_internal_time += timer.GetUserElapsed();
  timer.Start();
#endif
}
