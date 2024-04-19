#ifndef LSCPMTSD_hh
#define LSCPMTSD_hh

#include "G4DataVector.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4VSensitiveDetector.hh"

#include "LSCSim/PMTHit.hh"

class G4Step;
class G4HCofThisEvent;

class LSCPMTSD : public G4VSensitiveDetector {
public:
  LSCPMTSD(G4String name);
  virtual ~LSCPMTSD();

  virtual void Initialize(G4HCofThisEvent *);
  virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
  virtual void EndOfEvent(G4HCofThisEvent *);
  virtual void clear();
  void DrawAll();
  void PrintAll();

  void SimpleHit(G4int ipmt, G4double time, G4double kineticEnergy,
                 const G4ThreeVector & position, const G4ThreeVector & momentum,
                 const G4ThreeVector & polarization, G4int iHitPhotonCount);

private:
  PMTHitsCollection * fPMTHitCollection;
};

#endif
