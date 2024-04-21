#ifndef LSCScintillation_hh
#define LSCScintillation_hh

#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleMomentum.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4UImessenger.hh"
#include "G4VRestDiscreteProcess.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "templates.hh"
#include "G4EmSaturation.hh"

// #include "LSCSim/LSCEmSaturation.hh"

class G4UIdirectory;
class G4UIcmdWithABool;

class LSCScintillation : public G4VRestDiscreteProcess, public G4UImessenger {

public:
  ////////////////////////////////
  // Constructors and Destructor
  ////////////////////////////////

  LSCScintillation(const G4String & processName = "Scintillation",
                   G4ProcessType type = fElectromagnetic);
  ~LSCScintillation();

private:
  LSCScintillation(const LSCScintillation & right);

  //////////////
  // Operators
  //////////////

  LSCScintillation & operator=(const LSCScintillation & right);

public:
  ////////////
  // Methods
  ////////////

  virtual void SetNewValue(G4UIcommand *, G4String);

  // LSCScintillation Process has both PostStepDoIt (for energy
  // deposition of particles in flight) and AtRestDoIt (for energy
  // given to the medium by particles at rest)

  G4bool IsApplicable(const G4ParticleDefinition & aParticleType);
  // Returns true -> 'is applicable', for any particle type except
  // for an 'opticalphoton' and for short-lived particles

  G4double GetMeanFreePath(const G4Track & aTrack, G4double,
                           G4ForceCondition *);
  // Returns infinity; i. e. the process does not limit the step,
  // but sets the 'StronglyForced' condition for the DoIt to be
  // invoked at every step.

  G4double GetMeanLifeTime(const G4Track & aTrack, G4ForceCondition *);
  // Returns infinity; i. e. the process does not limit the time,
  // but sets the 'StronglyForced' condition for the DoIt to be
  // invoked at every step.

  G4VParticleChange * PostStepDoIt(const G4Track & aTrack,
                                   const G4Step & aStep);
  G4VParticleChange * AtRestDoIt(const G4Track & aTrack, const G4Step & aStep);

  // These are the methods implementing the scintillation process.

  void SetTrackSecondariesFirst(const G4bool state);
  // If set, the primary particle tracking is interrupted and any
  // produced scintillation photons are tracked next. When all
  // have been tracked, the tracking of the primary resumes.

  G4bool GetTrackSecondariesFirst() const;
  // Returns the boolean flag for tracking secondaries first.

  void SetScintillationYieldFactor(const G4double yieldfactor);
  // Called to set the scintillation photon yield factor, needed when
  // the yield is different for different types of particles. This
  // scales the yield obtained from the G4MaterialPropertiesTable.

  G4double GetScintillationYieldFactor() const;
  // Returns the photon yield factor.

  void SetScintillationExcitationRatio(const G4double excitationratio);
  // Called to set the scintillation exciation ratio, needed when
  // the scintillation level excitation is different for different
  // types of particles. This overwrites the YieldRatio obtained
  // from the G4MaterialPropertiesTable.

  G4double GetScintillationExcitationRatio() const;
  // Returns the scintillation level excitation ratio.

  G4PhysicsTable * GetEmitIntegralTable() const;
  // Returns the address of the fast scintillation integral table.

  void AddSaturation(G4EmSaturation * sat) { emSaturation = sat; }
  // Adds Birks Saturation to the process.

  void RemoveSaturation() { emSaturation = NULL; }
  // Removes the Birks Saturation from the process.

  G4EmSaturation * GetSaturation() const { return emSaturation; }
  // Returns the Birks Saturation.

  void SetScintillationByParticleType(const G4bool);
  // Called by the user to set the scintillation yield as a function
  // of energy deposited by particle type

  G4bool GetScintillationByParticleType() const
  {
    return scintillationByParticleType;
  }
  // Return the boolean that determines the method of scintillation
  // production

  void DumpPhysicsTable() const;
  // Prints the fast and slow scintillation integral tables.

  // step infomation
  G4double GetEnergyDeposit() const { return fTotalEdep; }
  G4double GetEnergyVisible() const { return fEvis; }
  G4double GetNScintillationPhoton() const { return fNScintillationPhoton; }

protected:
  void BuildThePhysicsTable();
  // It builds either the fast or slow scintillation integral table;
  // or both.

  ///////////////////////
  // Class Data Members
  ///////////////////////

  G4PhysicsTable * theEmitIntegralTable;

  G4bool fTrackSecondariesFirst;

  G4double YieldFactor;

  G4double ExcitationRatio;

  G4bool scintillationByParticleType;

private:
  G4UIdirectory * fScintillationDir;

  G4UIcmdWithABool * fScintOnCmd;

  G4bool fIsScintOn;

  // setep infomation
  G4double fTotalEdep;
  G4double fEvis;
  G4int fNScintillationPhoton;

  G4double single_exp(G4double t, G4double tau2);
  G4double bi_exp(G4double t, G4double tau1, G4double tau2);

  // emission time distribution when there is a finite rise time
  G4double sample_time(G4double tau1, G4double tau2);

  G4EmSaturation * emSaturation;
};

////////////////////
// Inline methods
////////////////////

inline G4bool
LSCScintillation::IsApplicable(const G4ParticleDefinition & aParticleType)
{
  if (aParticleType.GetParticleName() == "opticalphoton") return false;
  if (aParticleType.IsShortLived()) return false;

  return true;
}

inline void LSCScintillation::SetTrackSecondariesFirst(const G4bool state)
{
  fTrackSecondariesFirst = state;
}

inline G4bool LSCScintillation::GetTrackSecondariesFirst() const
{
  return fTrackSecondariesFirst;
}

inline void
LSCScintillation::SetScintillationYieldFactor(const G4double yieldfactor)
{
  YieldFactor = yieldfactor;
}

inline G4double LSCScintillation::GetScintillationYieldFactor() const
{
  return YieldFactor;
}

inline void LSCScintillation::SetScintillationExcitationRatio(
    const G4double excitationratio)
{
  ExcitationRatio = excitationratio;
}

inline G4double LSCScintillation::GetScintillationExcitationRatio() const
{
  return ExcitationRatio;
}

inline G4PhysicsTable * LSCScintillation::GetEmitIntegralTable() const
{
  return theEmitIntegralTable;
}

inline void LSCScintillation::DumpPhysicsTable() const
{
  if (theEmitIntegralTable) {
    G4int PhysicsTableSize = theEmitIntegralTable->entries();
    G4PhysicsOrderedFreeVector * v;

    for (G4int i = 0; i < PhysicsTableSize; i++) {
      v = (G4PhysicsOrderedFreeVector *)(*theEmitIntegralTable)[i];
      v->DumpValues();
    }
  }
}

inline G4double LSCScintillation::single_exp(G4double t, G4double tau2)
{
  return std::exp(-1.0 * t / tau2) / tau2;
}

inline G4double LSCScintillation::bi_exp(G4double t, G4double tau1,
                                         G4double tau2)
{
  return std::exp(-1.0 * t / tau2) * (1 - std::exp(-1.0 * t / tau1)) / tau2 /
         tau2 * (tau1 + tau2);
}

#endif /* LSCScintillation_h */
