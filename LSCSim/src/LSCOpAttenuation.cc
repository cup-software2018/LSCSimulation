#include "LSCSim/LSCOpAttenuation.hh"

#include "G4GeometryTolerance.hh"
#include "G4OpProcessSubType.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4WLSTimeGeneratorProfileDelta.hh"
#include "G4WLSTimeGeneratorProfileExponential.hh"
#include "G4ios.hh"

static const int N_COSTHETA_ENTRIES = 129;

static double Cos2ThetaTable[N_COSTHETA_ENTRIES];
static int TableInitialized = 0;

static void InitializeTable(void)
{
  double cos2th = 0.0;

  for (int i = 0; i < N_COSTHETA_ENTRIES - 1; i++) {
    double x = i / (double)(N_COSTHETA_ENTRIES - 1);
    double old_cos2th;
    do { // find exact root by iterating to convergence
      old_cos2th = cos2th;
      double costh = 2.0 * x / (3.0 - cos2th);
      cos2th = costh * costh;
    } while (fabs(old_cos2th - cos2th) >
             G4GeometryTolerance::GetInstance()->GetAngularTolerance());

    Cos2ThetaTable[i] = cos2th;
  }

  Cos2ThetaTable[N_COSTHETA_ENTRIES - 1] = 1.0;
  TableInitialized = 1;
}

LSCOpAttenuation::LSCOpAttenuation(const G4String & processName,
                                   G4ProcessType type)
    : G4VDiscreteProcess(processName, type)
{
  if (verboseLevel > 0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }

  SetProcessSubType(fOpAbsorption);

  if (!TableInitialized) InitializeTable();

  WLSTimeGeneratorProfile =
      new G4WLSTimeGeneratorProfileDelta("WLSTimeGeneratorProfileDelta");

  theIntegralTable = 0;
  BuildThePhysicsTable();
}

LSCOpAttenuation::~LSCOpAttenuation()
{
  if (theIntegralTable != 0) {
    theIntegralTable->clearAndDestroy();
    delete theIntegralTable;
  }

  delete WLSTimeGeneratorProfile;
}

G4VParticleChange * LSCOpAttenuation::PostStepDoIt(const G4Track & aTrack,
                                                   const G4Step & aStep)
{
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle * aParticle = aTrack.GetDynamicParticle();
  const G4Material * aMaterial = aTrack.GetMaterial();

  G4double thePhotonMomentum = aParticle->GetTotalMomentum();

  G4MaterialPropertiesTable * aMaterialPropertiesTable =
      aMaterial->GetMaterialPropertiesTable();
  if (!aMaterialPropertiesTable)
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  G4MaterialPropertyVector * OpScatFracVector =
      aMaterialPropertiesTable->GetProperty("OPSCATFRAC");

  G4double OpScatFrac = 0.0;
  if (OpScatFracVector) OpScatFrac = OpScatFracVector->Value(thePhotonMomentum);
  else {
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  // Scattering
  if (OpScatFrac > 0.0 && G4UniformRand() < OpScatFrac) {
    G4double urand = G4UniformRand() - 0.5;
    G4double Cos2Theta0 = Cos2ThetaTable[(
        int)(fabs(urand) * 2.0 * (N_COSTHETA_ENTRIES - 1) + 0.5)];
    G4double CosTheta = 4.0 * urand / (3.0 - Cos2Theta0);

#ifdef G4DEBUG
    if (fabs(CosTheta) > 1.0) {
      cerr << "GLG4OpAttenution: Warning, CosTheta=" << CosTheta
           << " urand=" << urand << endl;
      CosTheta = CosTheta > 0.0 ? 1.0 : -1.0;
    }
#endif

    G4double SinTheta = sqrt(1.0 - CosTheta * CosTheta);
    G4double Phi = (2.0 * G4UniformRand() - 1.0) * M_PI;
    G4ThreeVector e2(
        aParticle->GetMomentumDirection().cross(aParticle->GetPolarization()));

    G4ThreeVector NewMomentum =
        (CosTheta * aParticle->GetPolarization() +
         (SinTheta * cos(Phi)) * aParticle->GetMomentumDirection() +
         (SinTheta * sin(Phi)) * e2)
            .unit();

    // polarization is normal to new momentum and in same plane as
    // old new momentum and old polarization
    G4ThreeVector NewPolarization =
        (aParticle->GetPolarization() - CosTheta * NewMomentum).unit();

    aParticleChange.ProposeMomentumDirection(NewMomentum);
    aParticleChange.ProposePolarization(NewPolarization);

    if (verboseLevel > 0) G4cout << "\n** Photon scattered! **" << G4endl;
  }
  // Reemission or absoption
  else {
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    if (verboseLevel > 0) G4cout << "\n** Photon absorbed! **" << G4endl;

    const G4MaterialPropertyVector * WLS_Intensity =
        aMaterialPropertiesTable->GetProperty("WLSSPECTRUM");
    // No emission spectrum, just kill it (absorbed!)
    if (!WLS_Intensity) return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

    G4StepPoint * pPostStepPoint = aStep.GetPostStepPoint();

    G4int NumPhotons = 1;
    G4double QuntumEff = 1;
    if (aMaterialPropertiesTable->ConstPropertyExists("WLSPROBABILITY"))
      QuntumEff = aMaterialPropertiesTable->GetConstProperty("WLSPROBABILITY");

    if (G4UniformRand() > QuntumEff) {
      aParticleChange.SetNumberOfSecondaries(0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    aParticleChange.SetNumberOfSecondaries(NumPhotons);

    G4double primaryEnergy = aTrack.GetDynamicParticle()->GetKineticEnergy();

    // Retrieve the WLS Integral for this material
    // new G4PhysicsOrderedFreeVector allocated to hold CII's
    G4int materialIndex = aMaterial->GetIndex();

    G4double WLSTime = 0. * ns;
    G4PhysicsOrderedFreeVector * WLSIntegral = 0;

    if (aMaterialPropertiesTable->GetConstProperty("WLSTIMECONSTANT"))
      WLSTime = aMaterialPropertiesTable->GetConstProperty("WLSTIMECONSTANT");
    WLSIntegral =
        (G4PhysicsOrderedFreeVector *)((*theIntegralTable)(materialIndex));

    // Max WLS Integral
    G4double CIImax = WLSIntegral->GetMaxValue();
    G4double sampledEnergy;

    // Make sure the energy of the secondary is less than that of the primary
    for (G4int j = 1; j <= 100; j++) {
      // Determine photon energy
      G4double CIIvalue = G4UniformRand() * CIImax;
      sampledEnergy = WLSIntegral->GetEnergy(CIIvalue);
      if (verboseLevel > 1) {
        G4cout << "sampledEnergy = " << sampledEnergy << G4endl;
        G4cout << "CIIvalue =      " << CIIvalue << G4endl;
      }
      if (sampledEnergy <= primaryEnergy) break;
    }

    // If no such energy can be sampled, return
    if (sampledEnergy > primaryEnergy) {
      if (verboseLevel > 1)
        G4cout << " *** One less WLS photon will be returned ***" << G4endl;
      aParticleChange.SetNumberOfSecondaries(0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    // Generate random photon direction
    G4double cost = 1. - 2. * G4UniformRand();
    G4double sint = std::sqrt((1. - cost) * (1. + cost));

    G4double phi = twopi * G4UniformRand();
    G4double sinp = std::sin(phi);
    G4double cosp = std::cos(phi);

    G4double px = sint * cosp;
    G4double py = sint * sinp;
    G4double pz = cost;

    // Create photon momentum direction vector
    G4ParticleMomentum photonMomentum(px, py, pz);

    // Determine polarization of new photon
    G4double sx = cost * cosp;
    G4double sy = cost * sinp;
    G4double sz = -sint;

    G4ThreeVector photonPolarization(sx, sy, sz);

    G4ThreeVector perp = photonMomentum.cross(photonPolarization);

    phi = twopi * G4UniformRand();
    sinp = std::sin(phi);
    cosp = std::cos(phi);

    photonPolarization = cosp * photonPolarization + sinp * perp;

    photonPolarization = photonPolarization.unit();

    // Generate a new photon:
    G4DynamicParticle * aWLSPhoton =
        new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), photonMomentum);
    aWLSPhoton->SetPolarization(photonPolarization.x(), photonPolarization.y(),
                                photonPolarization.z());

    aWLSPhoton->SetKineticEnergy(sampledEnergy);

    // Generate new G4Track object:
    // Must give position of WLS optical photon
    G4double TimeDelay = WLSTimeGeneratorProfile->GenerateTime(WLSTime);
    G4double aSecondaryTime = (pPostStepPoint->GetGlobalTime()) + TimeDelay;

    G4ThreeVector aSecondaryPosition = pPostStepPoint->GetPosition();

    G4Track * aSecondaryTrack =
        new G4Track(aWLSPhoton, aSecondaryTime, aSecondaryPosition);

    aSecondaryTrack->SetTouchableHandle(aTrack.GetTouchableHandle());
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());

    aParticleChange.AddSecondary(aSecondaryTrack);

    if (verboseLevel > 0) {
      G4cout
          << "\n Exiting from LSCOpAttenuation::DoIt -- NumberOfSecondaries = "
          << aParticleChange.GetNumberOfSecondaries() << G4endl;
    }
  }

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

void LSCOpAttenuation::BuildThePhysicsTable()
{
  if (theIntegralTable) return;

  const G4MaterialTable * theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  // create new physics table
  if (!theIntegralTable) theIntegralTable = new G4PhysicsTable(numOfMaterials);

  // loop for materials
  for (G4int i = 0; i < numOfMaterials; i++) {
    G4PhysicsOrderedFreeVector * aPhysicsOrderedFreeVector =
        new G4PhysicsOrderedFreeVector();

    // Retrieve vector of WLS wavelength intensity for
    // the material from the material's optical properties table.
    G4Material * aMaterial = (*theMaterialTable)[i];

    G4MaterialPropertiesTable * aMaterialPropertiesTable =
        aMaterial->GetMaterialPropertiesTable();

    if (aMaterialPropertiesTable) {
      G4MaterialPropertyVector * theWLSVector =
          aMaterialPropertiesTable->GetProperty("WLSSPECTRUM");

      if (theWLSVector) {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*theWLSVector)[0];

        if (currentIN >= 0.0) {
          // Create first (photon energy)
          G4double currentPM = theWLSVector->Energy(0);
          G4double currentCII = 0.0;

          aPhysicsOrderedFreeVector->InsertValues(currentPM, currentCII);

          // Set previous values to current ones prior to loop
          G4double prevPM = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN = currentIN;

          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for (size_t j = 1; j < theWLSVector->GetVectorLength(); j++) {
            currentPM = theWLSVector->Energy(j);
            currentIN = (*theWLSVector)[j];

            currentCII = 0.5 * (prevIN + currentIN);
            currentCII = prevCII + (currentPM - prevPM) * currentCII;

            aPhysicsOrderedFreeVector->InsertValues(currentPM, currentCII);

            prevPM = currentPM;
            prevCII = currentCII;
            prevIN = currentIN;
          }
        }
      }
    }
    // The WLS integral for a given material
    // will be inserted in the table according to the
    // position of the material in the material table.

    theIntegralTable->insertAt(i, aPhysicsOrderedFreeVector);
  }
}

G4double LSCOpAttenuation::GetMeanFreePath(const G4Track & aTrack, G4double,
                                           G4ForceCondition *)
{
  const G4DynamicParticle * aParticle = aTrack.GetDynamicParticle();
  const G4Material * aMaterial = aTrack.GetMaterial();

  G4double thePhotonMomentum = aParticle->GetTotalMomentum();

  G4MaterialPropertiesTable * aMaterialPropertyTable;
  G4MaterialPropertyVector * AttenuationLengthVector;

  G4double AttenuationLength = DBL_MAX;

  aMaterialPropertyTable = aMaterial->GetMaterialPropertiesTable();
  if (aMaterialPropertyTable) {
    AttenuationLengthVector = aMaterialPropertyTable->GetProperty("ABSLENGTH");
    if (AttenuationLengthVector) {
      AttenuationLength = AttenuationLengthVector->Value(thePhotonMomentum);
    }
  }

  return AttenuationLength;
}

void LSCOpAttenuation::UseTimeProfile(const G4String name)
{
  if (name == "delta") {
    delete WLSTimeGeneratorProfile;
    WLSTimeGeneratorProfile = new G4WLSTimeGeneratorProfileDelta("delta");
  }
  else if (name == "exponential") {
    delete WLSTimeGeneratorProfile;
    WLSTimeGeneratorProfile =
        new G4WLSTimeGeneratorProfileExponential("exponential");
  }
  else
    G4Exception("LSCOpWLS::UseTimeProfile", "em0202", FatalException,
                "generator does not exist");
}
