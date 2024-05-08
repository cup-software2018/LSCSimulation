#include "LSCSim/LSCPMTSD.hh"

#include <iostream>

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include "MCObjs/MCPhotonHit.hh"
using namespace std;

LSCPMTSD::LSCPMTSD(G4String name)
    : G4VSensitiveDetector(name),
      fPMTHitCollection(0)
{
  collectionName.insert("pmtHitCollection");
}

LSCPMTSD::~LSCPMTSD() {}

void LSCPMTSD::Initialize(G4HCofThisEvent * hitsCE)
{
  fPMTHitCollection =
      new PMTHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Store collection with event and keep ID
  static G4int hitCID = -1;
  if (hitCID < 0) { hitCID = GetCollectionID(0); }
  hitsCE->AddHitsCollection(hitCID, fPMTHitCollection);
}

G4bool LSCPMTSD::ProcessHits(G4Step * step, G4TouchableHistory *)
{
  if (step->GetTrack()->GetDefinition() !=
      G4OpticalPhoton::OpticalPhotonDefinition())
    return false;

  step->GetTrack()->SetTrackStatus(fStopAndKill);

  int pmtNumber = 999;

  const G4VTouchable * touch = step->GetPreStepPoint()->GetTouchable();
  int nd = touch->GetHistoryDepth();

  for (int id = 0; id < nd; id++) {
    G4String pname = touch->GetVolume(id)->GetName();

    if (G4StrUtil::contains(pname, "PMTPhys"))
      pmtNumber = touch->GetReplicaNumber(id);
  }

  G4double globalTime = step->GetPostStepPoint()->GetGlobalTime();
  G4double kineticEnergy = step->GetPostStepPoint()->GetKineticEnergy();

  PMTHit * pmt = nullptr;

  G4int n = fPMTHitCollection->entries();
  for (G4int i = 0; i < n; i++) {
    if ((*fPMTHitCollection)[i]->GetPMTId() == pmtNumber) {
      pmt = (*fPMTHitCollection)[i];
      break;
    }
  }

  if (!pmt) {           // this pmt wasnt previously hit in this event
    pmt = new PMTHit(); // so create new hit
    pmt->SetPMTId(pmtNumber);
    fPMTHitCollection->insert(pmt);
  }

  MCPhotonHit * hit = pmt->AddHit();
  hit->SetTime((float)globalTime);
  hit->SetKineticEnergy((float)kineticEnergy);

  return true;
}

void LSCPMTSD::SimpleHit(G4int ipmt, G4double time, G4double kineticEnergy,
                         const G4ThreeVector & hit_position,
                         const G4ThreeVector & hit_momentum,
                         const G4ThreeVector & hit_polarization,
                         G4int iHitPhotonCount)
{
  PMTHit * pmt = nullptr;

  G4int n = fPMTHitCollection->entries();
  for (G4int i = 0; i < n; i++) {
    if ((*fPMTHitCollection)[i]->GetPMTId() == ipmt) {
      pmt = (*fPMTHitCollection)[i];
      break;
    }
  }

  if (!pmt) {           // this pmt wasnt previously hit in this event
    pmt = new PMTHit(); // so create new hit
    pmt->SetPMTId(ipmt);
    fPMTHitCollection->insert(pmt);
  }

  MCPhotonHit * hit = pmt->AddHit();
  hit->SetTime((float)time);
  hit->SetKineticEnergy((float)kineticEnergy);
}

void LSCPMTSD::EndOfEvent(G4HCofThisEvent *) {}

void LSCPMTSD::clear() {}

void LSCPMTSD::DrawAll() {}

void LSCPMTSD::PrintAll() {}
