#ifndef LSCEventAction_hh
#define LSCEventAction_hh

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class LSCRecorderBase;
class LSCEventAction : public G4UserEventAction {
public:
  LSCEventAction(LSCRecorderBase *);
  virtual ~LSCEventAction();

public:
  virtual void BeginOfEventAction(const G4Event *);
  virtual void EndOfEventAction(const G4Event *);

private:
  LSCRecorderBase * fRecorder;
};

#endif
