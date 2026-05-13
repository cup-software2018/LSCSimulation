// This file is part of the GenericLAND software library.
// $Id: GLG4PrimaryGeneratorAction.hh,v 1.1.1.1 2013/11/08 05:33:05 jslee Exp $
//
// new GLG4PrimaryGeneratorAction.hh by Glenn Horton-Smith, August 3-17, 2001

#ifndef __GLG4PrimaryGeneratorAction_hh__
#define __GLG4PrimaryGeneratorAction_hh__ 1
////////////////////////////////////////////////////////////////
// GLG4PrimaryGeneratorAction
////////////////////////////////////////////////////////////////

#include "G4VUserPrimaryGeneratorAction.hh" // for user primary vertex gen.

class GLG4PrimaryGeneratorMessenger;
class GLG4DetectorConstruction;
class G4Event;
class G4Track;
class G4String;
class GLG4VPosGen;
class GLG4VVertexGen;

class GLG4PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  // GLG4PrimaryGeneratorAction(GLG4DetectorConstruction *argDC); // not used
  GLG4PrimaryGeneratorAction();
  ~GLG4PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event * argEvent); // generate primary particles

  void DeferTrackToLaterEvent(const G4Track * track); // postpone to later

  void NotifyTimeToNextStackedEvent(double t); // note time to stacked evt

  double GetUniversalTime() { return myUniversalTime; }

  double GetUniversalTimeSincePriorEvent()
  {
    return myUniversalTimeSincePriorEvent;
  }

  int GetTypeOfCurrentEvent() { return myTypeOfCurrentEvent; }

  double GetEventRate(int i) { return myEventRate[i]; }
  void SetEventRate(int i, double r);

  int GetEventTriggerCondition(int iev) { return myEventTriggerCondition[iev]; }
  void SetEventTriggerCondition(int iev, int itc);

  double GetEventWindow() { return myEventWindow; }
  void SetEventWindow(double argEventWindow);

  double GetChainClip() { return myChainClip; }
  void SetChainClip(double argChainClip);

  static G4String GetEventTypeName(int argEventType);

  static const char * GetPositionCodeName(int argPosCode)
  {
    return thePositionCodeNames[argPosCode];
  }
  static const char * GetVertexCodeName(int argVertexCode)
  {
    return theVertexCodeNames[argVertexCode];
  }
  static int GetPositionCodeForEventType(int argEventType)
  {
    return theEventGeneratorCodes[argEventType].poscode;
  }
  static int GetVertexCodeForEventType(int argEventType)
  {
    return theEventGeneratorCodes[argEventType].vertexcode;
  }

  static GLG4VVertexGen * GetVertexGenerator(int i)
  {
    return theVertexGenerators[i];
  }
  static GLG4VPosGen * GetPositionGenerator(int i)
  {
    return thePositionGenerators[i];
  }

  static GLG4PrimaryGeneratorAction * GetTheGLG4PrimaryGeneratorAction()
  {
    return theGLG4PrimaryGeneratorAction;
  }

  enum {
    kGunEvtIndex = 3,
    kGunPosIndex = 9,
    kGunVtxIndex = 17,
    kDelayEvtIndex = 51,
    kDelayPosIndex = 12,
    kDelayVtxIndex = 19
  };
  enum {
    theNumEventTypes = 52,
    theNumPosGenCodes = 13,
    theNumVertexGenCodes = 20
  };
  enum {
    kGeneratorTriggerNormal = 0,
    kGeneratorTriggerPileupOnly = 1,
    kGeneratorTriggerDelay = 2
  };

private:
  // GLG4DetectorConstruction *myDetector;
  GLG4PrimaryGeneratorMessenger * myMessenger;

  double myUniversalTime;
  double myUniversalTimeSincePriorEvent;
  int myTypeOfCurrentEvent;
  double myEventWindow;
  double myChainClip;
  double myEventRate[theNumEventTypes];
  int myEventTriggerCondition[theNumEventTypes];
  double myTimeToNextEvent[theNumEventTypes];

  static GLG4VVertexGen * theVertexGenerators[theNumVertexGenCodes];
  static GLG4VPosGen * thePositionGenerators[theNumPosGenCodes];
  static const char * theVertexCodeNames[theNumVertexGenCodes];
  static const char * thePositionCodeNames[theNumPosGenCodes];
  typedef struct codepair_s {
    int poscode, vertexcode;
  } codepair_t;
  static codepair_t theEventGeneratorCodes[theNumEventTypes];
  static GLG4PrimaryGeneratorAction * theGLG4PrimaryGeneratorAction;
};

#endif
