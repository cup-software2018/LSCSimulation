/**@file
   Class for process of attenuation of optical photons

   Original by Jaison Lee, Nov. 2013.
   initiated from G4OpAbsorption, G4OpWLS and GLG4OpAttenuation

   G4OpAbsorption   : attenution
   G4OpWLS          : reemission
   GLGOpAttenuation : scattering

   @author Jaison Lee
*/

#ifndef LSCOpAttenuation_hh
#define LSCOpAttenuation_hh

#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "templates.hh"

class G4VWLSTimeGeneratorProfile;

class LSCOpAttenuation : public G4VDiscreteProcess {
public:
  LSCOpAttenuation(const G4String & processName = "Attenuation",
                   G4ProcessType type = fOptical);
  ~LSCOpAttenuation();

private:
  LSCOpAttenuation(const LSCOpAttenuation & right);
  LSCOpAttenuation & operator=(const LSCOpAttenuation & right);
  void BuildThePhysicsTable();

public:
  // Returns true -> 'is applicable' only for an optical photon.
  G4bool IsApplicable(const G4ParticleDefinition & aParticleType);

  // Returns the absorption length for bulk absorption of optical
  // photons in media with a specified attenuation length.
  G4double GetMeanFreePath(const G4Track & aTrack, G4double,
                           G4ForceCondition *);

  // This is the method implementing bulk absorption of optical
  // photons.
  G4VParticleChange * PostStepDoIt(const G4Track & aTrack,
                                   const G4Step & aStep);

  // Returns the address of the WLS integral table.
  G4PhysicsTable * GetIntegralTable() const;

  // Prints the WLS integral table.
  void DumpPhysicsTable() const;

  // Selects the time profile generator
  void UseTimeProfile(const G4String name);

protected:
  G4VWLSTimeGeneratorProfile * WLSTimeGeneratorProfile;
  G4PhysicsTable * theIntegralTable;
};

inline G4bool
LSCOpAttenuation::IsApplicable(const G4ParticleDefinition & aParticleType)
{
  return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

inline G4PhysicsTable * LSCOpAttenuation::GetIntegralTable() const
{
  return theIntegralTable;
}

inline void LSCOpAttenuation::DumpPhysicsTable() const
{
  G4int PhysicsTableSize = theIntegralTable->entries();
  G4PhysicsOrderedFreeVector * v;

  for (G4int i = 0; i < PhysicsTableSize; i++) {
    v = (G4PhysicsOrderedFreeVector *)(*theIntegralTable)[i];
    v->DumpValues();
  }
}
#endif
