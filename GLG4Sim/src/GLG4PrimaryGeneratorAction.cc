// This file is part of the GenericLAND software library.
// $Id: GLG4PrimaryGeneratorAction.cc,v 1.2 2013/11/09 23:48:54 jslee Exp $
//
// new GLG4PrimaryGeneratorAction.cc by Glenn Horton-Smith, August 3-17, 2001

////////////////////////////////////////////////////////////////
// GLG4PrimaryGeneratorAction
////////////////////////////////////////////////////////////////

// include files (boring)
#include "GLG4Sim/GLG4PrimaryGeneratorAction.hh"

#include <stdio.h> // for sprintf

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "TString.h"
#include "globals.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "GLG4Sim/GLG4PosGen.hh" // for global position generator
#include "GLG4Sim/GLG4PrimaryGeneratorMessenger.hh"
#include "GLG4Sim/GLG4VertexGen.hh" // for vertex generator
#include "GLG4Sim/GLG4param.hh"     // for GLG4param

using namespace CLHEP;

// here are the static constants and variables (boring)
const char *
    GLG4PrimaryGeneratorAction::thePositionCodeNames[theNumPosGenCodes] = {
        "Scintillator",          // 0
        "Rock-shell",            // 1
        "Rope",                  // 2
        "Chimney",               // 3
        "LSC-vessel",            // 4
        "OB-to-vessel",          // 5
        "OD-to-vessel",          // 6
        "Rock-to-vessel",        // 7
        "Cosmic",                // 8
        "User-controlled",       // 9
        "Spare-user-controlled", // 10
        "Inner-buffer",          // 11
        "Delayed-particle"};     // 12

const char *
    GLG4PrimaryGeneratorAction::theVertexCodeNames[theNumVertexGenCodes] = {
        "reactor-anti-neutrino",
        "terrestrial-anti-neutrino",
        "solar-neutrino",
        "Uranium",
        "Thorium",
        "Radon",
        "Thoron",
        "K-40",
        "Co-60",
        "Kr-85",
        "C-10",
        "C-11",
        "C-14",
        "muons",
        "external-from-OB",
        "external-from-OD",
        "external-from-rock",
        "test-gun",
        "spare-HEPevt",
        "Delayed-particle"};

GLG4PrimaryGeneratorAction::codepair_t
    GLG4PrimaryGeneratorAction::theEventGeneratorCodes[theNumEventTypes] = {
        {0, 0},   {0, 1},   {0, 2},  {9, 17},  {2, 3},  {2, 4},  {2, 5},
        {2, 6},   {2, 7},   {3, 3},  {3, 4},   {3, 5},  {3, 6},  {3, 7},
        {3, 8},   {4, 3},   {4, 4},  {4, 5},   {4, 6},  {4, 7},  {5, 14},
        {6, 15},  {7, 16},  {8, 13}, {0, 3},   {0, 4},  {0, 5},  {0, 6},
        {0, 7},   {0, 9},   {0, 10}, {0, 11},  {0, 12}, {1, 3},  {1, 4},
        {1, 5},   {1, 6},   {1, 7},  {10, 18}, {11, 0}, {11, 1}, {11, 2},
        {11, 3},  {11, 4},  {11, 5}, {11, 6},  {11, 7}, {11, 9}, {11, 10},
        {11, 11}, {11, 12}, {12, 19}};

GLG4VVertexGen *
    GLG4PrimaryGeneratorAction::theVertexGenerators[theNumVertexGenCodes];

GLG4VPosGen *
    GLG4PrimaryGeneratorAction::thePositionGenerators[theNumPosGenCodes];

GLG4PrimaryGeneratorAction *
    GLG4PrimaryGeneratorAction::theGLG4PrimaryGeneratorAction = 0;

// static functions that are not inline (boring)
G4String GLG4PrimaryGeneratorAction::GetEventTypeName(int argEventType)
{
  G4String name(
      thePositionCodeNames[theEventGeneratorCodes[argEventType].poscode]);
  name += ":";
  name += (theVertexCodeNames[theEventGeneratorCodes[argEventType].vertexcode]);
  return name;
}

// constructor and destructor (mildly interesting because generators created)
// GLG4PrimaryGeneratorAction::GLG4PrimaryGeneratorAction(
//    GLG4DetectorConstruction * argDC)
//    : myDetector(argDC)
GLG4PrimaryGeneratorAction::GLG4PrimaryGeneratorAction()
{
  if (theGLG4PrimaryGeneratorAction == 0) {
    theGLG4PrimaryGeneratorAction = this;
  }
  else {
    G4Exception(
        "Error, more than one GLG4PrimaryGeneratorAction instantiated.\n"
        "Sorry, but this is a no-no because GLG4SteppingAction relies on\n"
        "GLG4PrimaryGeneratorAction::GetTheGLG4PrimaryGeneratorAction().\n"
        "This is yucky, I know -- please rewrite GLG4SteppingAction AND\n"
        "all main() programs so that constructor accepts a pointer to\n"
        "the GLG4PrimaryGeneratorAction you really want them to use.",
        "", FatalException, "");
  }

  // initialize messenger and time fields
  myMessenger = new GLG4PrimaryGeneratorMessenger(this);
  myUniversalTime = 0.0;
  myUniversalTimeSincePriorEvent = 0.0;

  // initialize generator state
  GLG4param & db(GLG4param::GetDB());
  myEventWindow = db.GetWithDefault("gen.eventWindow", 1000. * ns);
  myChainClip = db.GetWithDefault("gen.chainClip", 100.0 * second);

  for (int i = 0; i < theNumEventTypes; i++) {
    char key[16];
    sprintf(key, "gen.rate%d", i);
    myEventRate[i] = db.GetWithDefault(key, 0.0);
    sprintf(key, "gen.trig%d", i);
    switch (i) {
      case 32: // C-14 in scintillator
      case 33: // U-238  in rock
      case 34: // Th-232 in rock
      case 35: // Radon  in rock
      case 36: // Thoron in rock
      case 37: // K-40   in rock
        myEventTriggerCondition[i] =
            (int)(db.GetWithDefault(key, kGeneratorTriggerPileupOnly));
        break;
      case kDelayEvtIndex:
        myEventTriggerCondition[i] =
            (int)(db.GetWithDefault(key, kGeneratorTriggerDelay));
        break;
      default:
        myEventTriggerCondition[i] =
            (int)(db.GetWithDefault(key, kGeneratorTriggerNormal));
        break;
    }
    myTimeToNextEvent[i] = -1.0;
  }

  // set generator arrays -- classes and initial state are hard-coded
  // (but user can change state manually, and database does not track-FIXME!)
  if (theVertexGenerators[0] == 0) {
    // kGunEvtIndex == 3 is assumed in many scripts and such.  Don't change it.
    if (kGunEvtIndex != 3) {
      G4Exception("Error, kGunIndex != 3 in GLG4PrimaryGeneratorAction.", "",
                  FatalException, "");
    }

    // set the vertex generator array
    theVertexGenerators[17] = new GLG4VertexGen_Gun("gen.vtx17");
    theVertexGenerators[kDelayVtxIndex] =
        new GLG4VertexGen_Stack("gen.vtx19", this);
    int i;
    for (i = 0; i < theNumVertexGenCodes; i++)
      if (theVertexGenerators[i] == NULL) {
        char dbname[16];
        sprintf(dbname, "gen.vtx%d", i);
        theVertexGenerators[i] = new GLG4VertexGen_HEPEvt(dbname);
      }

    // set the position generator array
    thePositionGenerators[5] = new GLG4PosGen_null("gen.pos5"); // OB-to-balloon
    thePositionGenerators[6] = new GLG4PosGen_null("gen.pos6"); // OD-to-balloon
    thePositionGenerators[7] =
        new GLG4PosGen_null("gen.pos7"); // Rock-to-balloon
    thePositionGenerators[8] =
        new GLG4PosGen_Cosmic("gen.pos8"); // external cosmic rays
    thePositionGenerators[kDelayPosIndex] =
        new GLG4PosGen_null("gen.pos12"); // delayed track stack
    for (i = 0; i < theNumPosGenCodes; i++)
      if (thePositionGenerators[i] == NULL) {
        char dbname[16];
        sprintf(dbname, "gen.pos%d", i);
        thePositionGenerators[i] = new GLG4PosGen_PointPaintFill(dbname);
      }

    // set reasonable defaults for each position generator:
    // fill scintillator
    thePositionGenerators[0]->SetState("0 0 0 fill");
    // fill rock shell
    thePositionGenerators[1]->SetState("0 0 -10100 fill");
    // fill Kevlar within 5 mm of balloon (custom generator needed otherwise)
    thePositionGenerators[2]->SetState("0 0 0 paint ! 5 Kevlar");
    // fill the stainless within 10 cm of the balloon surface
    //  (chimney, bellmouth, flange and part of pipe on bottom)
    thePositionGenerators[3]->SetState("0 0 0 paint ! 100 StainlessSteel");
    // uniformly coat the balloon film, 100 um thick
    thePositionGenerators[4]->SetState("0 0 0 paint ! 0.1");
    //  [nothing needed for 5..7, the positions are set in the VtxGen]
    // cosmic rays in 30-meter-tall by 20-meter-wide rectangle
    thePositionGenerators[8]->SetState("20000 30000");
    //  [defaults are okay for 9,10, the User-controlled particle sources]
    // fill inner buffer
    thePositionGenerators[11]->SetState("7000 0 0 fill");
  }
}

GLG4PrimaryGeneratorAction::~GLG4PrimaryGeneratorAction() {}

// "Set" functions (non-inline because of GLG4param interaction)
void GLG4PrimaryGeneratorAction::SetEventRate(int i, double r)
{
  myEventRate[i] = r;

  GLG4param & db(GLG4param::GetDB());
  char key[16];
  sprintf(key, "gen.rate%d", i);
  db[key] = r;
}

void GLG4PrimaryGeneratorAction::SetEventTriggerCondition(int iev, int itc)
{
  myEventTriggerCondition[iev] = itc;

  GLG4param & db(GLG4param::GetDB());
  char key[16];
  sprintf(key, "gen.trig%d", iev);
  db[key] = itc;
}

void GLG4PrimaryGeneratorAction::SetEventWindow(double argEventWindow)
{
  myEventWindow = argEventWindow;

  GLG4param & db(GLG4param::GetDB());
  db["gen.eventWindow"] = argEventWindow;
}

void GLG4PrimaryGeneratorAction::SetChainClip(double argChainClip)
{
  myChainClip = argChainClip;

  GLG4param & db(GLG4param::GetDB());
  db["gen.chainClip"] = argChainClip;
}

// GeneratePrimaries (this is the interesting part!)
void GLG4PrimaryGeneratorAction::GeneratePrimaries(G4Event * argEvent)
{
  int next_event_type = -1;
  G4double min_time_to_next_event = DBL_MAX;

  // find the next event, resetting any "expired" events along the way
  for (int i = 0; i < theNumEventTypes; i++) {

    if (myTimeToNextEvent[i] < 0.0 && myEventRate[i] > 0.0 &&
        myEventTriggerCondition[i] == kGeneratorTriggerNormal) {
      myTimeToNextEvent[i] = -log(1.0 - G4UniformRand()) / myEventRate[i];
      // myTimeToNextEvent[i] = 1.0 / myEventRate[i];
      //G4cout << Form("1-- %5d %20e %20e", i, myTimeToNextEvent[i],
      //               myEventRate[i])
      //       << G4endl;
    }
    if (myTimeToNextEvent[i] >= 0.0 &&
        myTimeToNextEvent[i] < min_time_to_next_event) {
      next_event_type = i;
      min_time_to_next_event = myTimeToNextEvent[i];
      //G4cout << Form("2-- %5d %20e %20e", i, myTimeToNextEvent[i],
      //               myEventRate[i])
      //       << G4endl;
    }
  }

  // event type -1 means a commonly-made user error
  if (next_event_type == -1) {
    G4cerr << "GLG4PrimaryGeneratorAction: no non-zero event rates!\n"
           << "Setting gun rate to 1 Hz (use /generator/rates to reset)"
           << G4endl;
    next_event_type = kGunEvtIndex;
    min_time_to_next_event = 0.0;
    myEventRate[next_event_type] = 1e-9;
  }

  // update universal time and decrement all time-to-next-events
  myUniversalTimeSincePriorEvent = min_time_to_next_event;
  myUniversalTime += min_time_to_next_event;

  for (int i = 0; i < theNumEventTypes; i++) {
    myTimeToNextEvent[i] -= min_time_to_next_event;
  }

  myTimeToNextEvent[next_event_type] =
      -log(1.0 - G4UniformRand()) /
      myEventRate[next_event_type]; // new time-to-next
  // myTimeToNextEvent[next_event_type] =
  //     1.0 / myEventRate[next_event_type]; // new time-to-next
  myTypeOfCurrentEvent = next_event_type;

  //G4cout << Form("%5d %20e", next_event_type,
  //               myTimeToNextEvent[next_event_type])
  //       << G4endl;

  // generate the event!
  int next_vtx_code = theEventGeneratorCodes[next_event_type].vertexcode;
  int next_pos_code = theEventGeneratorCodes[next_event_type].poscode;
  // "vertex"
  theVertexGenerators[next_vtx_code]->GeneratePrimaryVertex(argEvent);
  G4PrimaryVertex * v = argEvent->GetPrimaryVertex(0);
  // and "position"
  if (v != 0) {
    thePositionGenerators[next_pos_code]->GenerateVertexPositions(
        v, myChainClip, myEventRate[next_event_type]);
  }
  else {
    G4cerr << "Warning, no vertex generated by vertex generator:"
           << " event_type=" << next_event_type << " vtx_code=" << next_vtx_code
           << G4endl;
    v = new G4PrimaryVertex(0., 0., 0., 0.);
    argEvent->AddPrimaryVertex(v);
  }

  // Note: if you want a fake "Informaton" in your vertex lists which
  // contains some internal generator state in the "px, py, pz" components,
  // then you can uncomment the next seven lines:
  // // add virtual particle containing time generator info
  // enum { kChronatonPseudoPDGcode = -1999999999 };
  // v->SetPrimary
  //   (new G4PrimaryParticle( kChronatonPseudoPDGcode,
  // 			    myUniversalTime,        // univeral time of event
  // 			    min_time_to_next_event, // time since last event
  // 			    (double)next_event_type ));// type of this event

  // add pileup events from normal-triggering primary events
  for (int i = 0; i < theNumEventTypes; i++) {
    while (myTimeToNextEvent[i] >= 0.0 &&
           myTimeToNextEvent[i] < myEventWindow) {
      int vtx_code = theEventGeneratorCodes[i].vertexcode;
      int pos_code = theEventGeneratorCodes[i].poscode;
      int nstart = argEvent->GetNumberOfPrimaryVertex();
      theVertexGenerators[vtx_code]->GeneratePrimaryVertex(argEvent);
      G4PrimaryVertex * vn = argEvent->GetPrimaryVertex(nstart);
      if (vn)
        thePositionGenerators[pos_code]->GenerateVertexPositions(
            vn, myChainClip, myEventRate[i], myTimeToNextEvent[i]);
      myTimeToNextEvent[i] += -log(1.0 - G4UniformRand()) / myEventRate[i];
    }
  }

  // add pileup events from pileup-triggered events (which may have rates >
  // 1/myEventWindow)
  for (int i = 0; i < theNumEventTypes; i++) {
    if (myEventTriggerCondition[i] == kGeneratorTriggerPileupOnly) {
      int vtx_code = theEventGeneratorCodes[i].vertexcode;
      int pos_code = theEventGeneratorCodes[i].poscode;
      double t = 0.0;
      while ((t += -log(1.0 - G4UniformRand()) / myEventRate[i]) <
             myEventWindow) {
        int nstart = argEvent->GetNumberOfPrimaryVertex();
        theVertexGenerators[vtx_code]->GeneratePrimaryVertex(argEvent);
        G4PrimaryVertex * vn = argEvent->GetPrimaryVertex(nstart);
        if (vn)
          thePositionGenerators[pos_code]->GenerateVertexPositions(
              vn, myChainClip, myEventRate[i], t);
      }
    }
  }
  // done!
}

void GLG4PrimaryGeneratorAction::DeferTrackToLaterEvent(const G4Track * track)
{
  NotifyTimeToNextStackedEvent(track->GetGlobalTime());
  ((GLG4VertexGen_Stack *)(theVertexGenerators[19]))->StackIt(track);
  // track->SetTrackStatus( fStopAndKill ); // must be done by caller
}

void GLG4PrimaryGeneratorAction::NotifyTimeToNextStackedEvent(double t)
{
  if (myTimeToNextEvent[kDelayEvtIndex] < 0.0 ||
      t < myTimeToNextEvent[kDelayEvtIndex])
    myTimeToNextEvent[kDelayEvtIndex] = t;
}
