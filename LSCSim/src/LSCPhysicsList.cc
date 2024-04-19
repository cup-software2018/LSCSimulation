#include "LSCSim/LSCPhysicsList.hh"

#include <iomanip>

#include "G4DeexPrecoParameters.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4NuclearLevelData.hh"
#include "G4NuclideTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UserLimits.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "GLG4Sim/GLG4DeferTrackProc.hh"
#include "GLG4Sim/GLG4param.hh"

// Constructor /////////////////////////////////////////////////////////////
LSCPhysicsList::LSCPhysicsList()
    : G4VModularPhysicsList()
{
  defaultCutValue = 1.0 * micrometer; //
  cutForGamma = defaultCutValue;
  cutForElectron = 1.0 * nanometer;
  cutForPositron = defaultCutValue;

  VerboseLevel = 1;
  OpVerbLevel = 0;

  SetVerboseLevel(VerboseLevel);

  // mandatory for G4NuclideTable
  //
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1 * picosecond);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0 * eV);

  // read new PhotonEvaporation data set
  //
  G4DeexPrecoParameters * deex =
      G4NuclearLevelData::GetInstance()->GetParameters();
  //  deex->SetCorrelatedGamma(false);
  deex->SetCorrelatedGamma(true);
  deex->SetStoreAllLevels(true);
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife() /
                       std::log(2.));

  // for messenger
  physDir = new G4UIdirectory("/LSC/phys/");
  physDir->SetGuidance("PhysicsList control");

  gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/LSC/phys/CutGamma", this);
  gammaCutCmd->SetGuidance("Set gamma cut.");
  gammaCutCmd->SetParameterName("Gcut", false);
  gammaCutCmd->SetUnitCategory("Length");
  gammaCutCmd->SetRange("Gcut>0.0");
  gammaCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  electCutCmd = new G4UIcmdWithADoubleAndUnit("/LSC/phys/CutEl", this);
  electCutCmd->SetGuidance("Set electron cut.");
  electCutCmd->SetParameterName("Ecut", false);
  electCutCmd->SetUnitCategory("Length");
  electCutCmd->SetRange("Ecut>0.0");
  electCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  posCutCmd = new G4UIcmdWithADoubleAndUnit("/LSC/phys/CutPos", this);
  posCutCmd->SetGuidance("Set positron cut.");
  posCutCmd->SetParameterName("Pcut", false);
  posCutCmd->SetUnitCategory("Length");
  posCutCmd->SetRange("Pcut>0.0");
  posCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  allCutCmd = new G4UIcmdWithADoubleAndUnit("/LSC/phys/CutsAll", this);
  allCutCmd->SetGuidance("Set cut for all.");
  allCutCmd->SetParameterName("cut", false);
  allCutCmd->SetUnitCategory("Length");
  allCutCmd->SetRange("cut>0.0");
  allCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  verboseCmd = new G4UIcmdWithAnInteger("/LSC/phys/verbose", this);
  verboseCmd->SetGuidance("set verbose for physics processes");
  verboseCmd->SetParameterName("verbose", true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("verbose>=0");
  verboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

// Destructor //////////////////////////////////////////////////////////////
LSCPhysicsList::~LSCPhysicsList()
{
  delete physDir;
  delete gammaCutCmd;
  delete electCutCmd;
  delete posCutCmd;
  delete allCutCmd;
  delete verboseCmd;
}

// Construct Particles /////////////////////////////////////////////////////
void LSCPhysicsList::ConstructParticle()
{

  // In this method, static member functions should be called
  // for all particles which yodu want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  ConstructMyBosons();
  ConstructMyLeptons();
  ConstructMyHadrons();
  ConstructMyShortLiveds();
}

// construct Bosons://///////////////////////////////////////////////////
void LSCPhysicsList::ConstructMyBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // OpticalPhotons
  G4OpticalPhoton::OpticalPhotonDefinition();
}

// construct Leptons://///////////////////////////////////////////////////
void LSCPhysicsList::ConstructMyLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4MesonConstructor.hh"

// construct Hadrons://///////////////////////////////////////////////////
void LSCPhysicsList::ConstructMyHadrons()
{
  //  mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  //  baryons
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

  //  ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}

#include "G4ShortLivedConstructor.hh"
// construct Shortliveds://///////////////////////////////////////////////////
void LSCPhysicsList::ConstructMyShortLiveds()
{
  // ShortLiveds
  G4ShortLivedConstructor slConstructor;
  slConstructor.ConstructParticle();
}

// Construct Processes //////////////////////////////////////////////////////
void LSCPhysicsList::ConstructProcess()
{

  AddTransportation();

  AddParameterisation();

  ConstructEM();
  // if (emName == "livermore") {
  //   auto a = new G4EmLivermorePhysics;
  //   a->ConstructProcess();
  // }
  // else {
  //   ConstructEM();
  G4cout << "EM Physics is ConstructEM(): default!" << G4endl;
  //}

  ConstructOp();

  // ConstructHad();

  ConstructGeneral();
}

// G4FastSimulation Processes
// //////////////////////////////////////////////////////
#include "G4FastSimulationManagerProcess.hh"

void LSCPhysicsList::AddParameterisation()
{
  auto theParticleIterator = GetParticleIterator();

  G4FastSimulationManagerProcess * theFastSimulationManagerProcess =
      new G4FastSimulationManagerProcess();
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition * particle = theParticleIterator->value();
    G4ProcessManager * pmanager = particle->GetProcessManager();
    // both postStep and alongStep action are required if the detector
    // makes use of ghost volumes. If no ghost, the postStep
    // is sufficient (and faster?).
#define Cup_USES_GHOST_VOLUMES 0
#if Cup_USES_GHOST_VOLUMES
    pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);
#else
    pmanager->AddProcess(theFastSimulationManagerProcess, -1, -1, 1);
#endif
  }
}

// Electromagnetic Processes ////////////////////////////////////////////////
// all charged particles

// gamma
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"

// e-
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"

// e+
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eplusAnnihilation.hh"

// alpha and GenericIon and deuterons, triton, He3:
#include "G4EnergyLossTables.hh"

// muon:
#include "G4MuBremsstrahlung.hh"
#include "G4MuIonisation.hh"
#include "G4MuPairProduction.hh"
// #include "G4MuonMinusCaptureAtRest.hh"
#include "G4MuonMinusCapture.hh"

// OTHERS:
#include "G4IonParametrisedLossModel.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"

// em process options to allow msc step-limitation to be switched off
// #include "G4EmProcessOptions.hh"
#include "G4LossTableManager.hh"

void LSCPhysicsList::ConstructEM()
{

  // set a finer grid of the physic tables in order to improve precision
  // former LowEnergy models have 200 bins up to 100 GeV
  // G4EmProcessOptions opt;
  // opt.SetMaxEnergy(100 * GeV);
  // opt.SetDEDXBinning(200);
  // opt.SetLambdaBinning(200);

  auto theParticleIterator = GetParticleIterator();

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition * particle = theParticleIterator->value();
    G4ProcessManager * pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    G4String particleType = particle->GetParticleType();
    G4double charge = particle->GetPDGCharge();

    if (particleName == "gamma") {
      // gamma
      G4RayleighScattering * theRayleigh = new G4RayleighScattering();
      theRayleigh->SetEmModel(
          new G4LivermoreRayleighModel()); // not strictly necessary
      pmanager->AddDiscreteProcess(theRayleigh);

      G4PhotoElectricEffect * thePhotoElectricEffect =
          new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);

      G4ComptonScattering * theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
      pmanager->AddDiscreteProcess(theComptonScattering);

      G4GammaConversion * theGammaConversion = new G4GammaConversion();
      theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
      pmanager->AddDiscreteProcess(theGammaConversion);
    }
    else if (particleName == "e-") {
      // electron
      // process ordering: AddProcess(name, at rest, along step, post step)
      // Multiple scattering
      G4eMultipleScattering * msc = new G4eMultipleScattering();
      pmanager->AddProcess(msc, -1, 1, 1);

      // Ionisation
      G4eIonisation * eIonisation = new G4eIonisation();
      eIonisation->SetEmModel(new G4LivermoreIonisationModel());
      eIonisation->SetStepFunction(0.2,
                                   100 * um); // improved precision in tracking
      pmanager->AddProcess(eIonisation, -1, 2, 2);

      // Bremsstrahlung
      G4eBremsstrahlung * eBremsstrahlung = new G4eBremsstrahlung();
      eBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel());
      pmanager->AddProcess(eBremsstrahlung, -1, -3, 3);
    }
    else if (particleName == "e+") {
      // positron
      G4eMultipleScattering * msc = new G4eMultipleScattering();
      pmanager->AddProcess(msc, -1, 1, 1);

      // Ionisation
      G4eIonisation * eIonisation = new G4eIonisation();
      eIonisation->SetStepFunction(0.2, 100 * um); //
      pmanager->AddProcess(eIonisation, -1, 2, 2);

      // Bremsstrahlung (use default, no low-energy available)
      pmanager->AddProcess(new G4eBremsstrahlung(), -1, -1, 3);

      // Annihilation
      pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1, 4);
    }
    else if (particleName == "mu+" || particleName == "mu-") {
      // muon
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation(), -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung(), -1, -1, 3);
      pmanager->AddProcess(new G4MuPairProduction(), -1, -1, 4);
      if (particleName == "mu-")
        // pmanager->AddProcess(new G4MuonMinusCaptureAtRest(), 0,-1,-1);
        pmanager->AddProcess(new G4MuonMinusCapture(), 0, -1, -1);
    }
    else if (particleName == "proton" || particleName == "pi+" ||
             particleName == "pi-") {
      // multiple scattering
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      // ionisation
      G4hIonisation * hIonisation = new G4hIonisation();
      hIonisation->SetStepFunction(0.2, 50 * um);
      pmanager->AddProcess(hIonisation, -1, 2, 2);

      // bremmstrahlung
      pmanager->AddProcess(new G4hBremsstrahlung, -1, -3, 3);
    }
    else if (particleName == "alpha" || particleName == "deuteron" ||
             particleName == "triton" || particleName == "He3") {
      // multiple scattering
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      // ionisation
      G4ionIonisation * ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 20 * um);
      pmanager->AddProcess(ionIoni, -1, 2, 2);
    }
    else if (particleName == "GenericIon") {
      // OBJECT may be dynamically created as either a GenericIon or nucleus
      // G4Nucleus exists and therefore has particle type nucleus
      // genericIon:

      // multiple scattering
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      // ionisation
      G4ionIonisation * ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 20 * um);
      pmanager->AddProcess(ionIoni, -1, 2, 2);
    }

    else if ((!particle->IsShortLived()) && (charge != 0.0) &&
             (particle->GetParticleName() != "chargedgeantino")) {
      // all others charged particles except geantino
      G4hMultipleScattering * aMultipleScattering = new G4hMultipleScattering();
      G4hIonisation * ahadronIon = new G4hIonisation();

      // multiple scattering
      pmanager->AddProcess(aMultipleScattering, -1, 1, 1);

      // ionisation
      pmanager->AddProcess(ahadronIon, -1, 2, 2);
    }
  }

  // switch on fluorescence, PIXE and Auger:
  // opt.SetFluo(true);
  // opt.SetPIXE(true);
  // opt.SetAuger(true);

  // Deexcitation
  //
  G4VAtomDeexcitation * de = new G4UAtomicDeexcitation();
  de->SetFluo(true);
  de->SetAuger(true);
  de->SetPIXE(true);
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

// Optical Processes ////////////////////////////////////////////////////////
// EJ: start
#include "G4Cerenkov.hh"
#include "G4EmSaturation.hh"
#include "G4OpBoundaryProcess.hh"

#include "LSCSim/LSCOpAttenuation.hh"
#include "LSCSim/LSCScintillation.hh"
// EJ: end

void LSCPhysicsList::ConstructOp()
{
  // EJ: start
  // scintillation process
  LSCScintillation * theScintProcessDef = new LSCScintillation("Scintillation");
  // theScintProcessDef->DumpPhysicsTable();
  theScintProcessDef->SetTrackSecondariesFirst(true);
  theScintProcessDef->SetScintillationYieldFactor(1.0);     //
  theScintProcessDef->SetScintillationExcitationRatio(0.0); //
  theScintProcessDef->SetVerboseLevel(OpVerbLevel);

  G4EmSaturation * emSaturation =
      G4LossTableManager::Instance()->EmSaturation();
  theScintProcessDef->AddSaturation(emSaturation);

  /*
    // scintillation process for alpha:
    CupScintillation * theScintProcessAlpha =
        new CupScintillation("Scintillation");
    // theScintProcessNuc->DumpPhysicsTable();
    theScintProcessAlpha->SetTrackSecondariesFirst(true);
    theScintProcessAlpha->SetScintillationYieldFactor(1.1);
    theScintProcessAlpha->SetScintillationExcitationRatio(1.0);
    theScintProcessAlpha->SetVerboseLevel(OpVerbLevel);

    theScintProcessAlpha->AddSaturation(emSaturation);

    // scintillation process for heavy nuclei
    CupScintillation * theScintProcessNuc = new
    CupScintillation("Scintillation");
    // theScintProcessNuc->DumpPhysicsTable();
    theScintProcessNuc->SetTrackSecondariesFirst(true);
    theScintProcessNuc->SetScintillationYieldFactor(0.2);
    theScintProcessNuc->SetScintillationExcitationRatio(1.0);
    theScintProcessNuc->SetVerboseLevel(OpVerbLevel);

    theScintProcessNuc->AddSaturation(emSaturation);
  */

  // optical processes
  LSCOpAttenuation * theAttenuationProcess = new LSCOpAttenuation();
  theAttenuationProcess->UseTimeProfile("exponential");
  theAttenuationProcess->SetVerboseLevel(OpVerbLevel);

  G4OpBoundaryProcess * theBoundaryProcess = new G4OpBoundaryProcess();
  theBoundaryProcess->SetVerboseLevel(OpVerbLevel);

  // Cerenkov
  G4Cerenkov * theCerenkovProcess = new G4Cerenkov();
  theCerenkovProcess->SetTrackSecondariesFirst(true);

  auto theParticleIterator = GetParticleIterator();

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition * particle = theParticleIterator->value();
    G4ProcessManager * pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (theScintProcessDef->IsApplicable(*particle)) {
      pmanager->AddProcess(theScintProcessDef);
      pmanager->SetProcessOrderingToLast(theScintProcessDef, idxAtRest);
      pmanager->SetProcessOrderingToLast(theScintProcessDef, idxPostStep);
      pmanager->AddProcess(theCerenkovProcess);
      pmanager->SetProcessOrdering(theCerenkovProcess, idxPostStep);
    }

    if (particleName == "opticalphoton") {
      pmanager->AddDiscreteProcess(theAttenuationProcess);
      pmanager->AddDiscreteProcess(theBoundaryProcess);
    }
  }
}

// Hadronic processes ////////////////////////////////////////////////////////
void LSCPhysicsList::ConstructHad() {}

// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4RadioactiveDecay.hh"

#include "GLG4Sim/GLG4DeferTrackProc.hh"

void LSCPhysicsList::ConstructGeneral()
{
  auto theParticleIterator = GetParticleIterator();

  // Add Decay Process
  G4Decay * theDecayProcess = new G4Decay();
  GLG4DeferTrackProc * theDeferProcess = new GLG4DeferTrackProc();
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition * particle = theParticleIterator->value();
    G4ProcessManager * pmanager = particle->GetProcessManager();

    if (theDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) {
      pmanager->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
    if (!particle->IsShortLived())
      pmanager->AddDiscreteProcess(theDeferProcess);
  }

  // Declare radioactive decay to the GenericIon in the IonTable.
  const G4IonTable * theIonTable =
      G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay * theRadioactiveDecay = new G4RadioactiveDecay();

  for (G4int i = 0; i < theIonTable->Entries(); i++) {
    G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
    G4String particleType = theIonTable->GetParticle(i)->GetParticleType();

    if (particleName == "GenericIon") {
      G4ProcessManager * pmanager =
          theIonTable->GetParticle(i)->GetProcessManager();
      pmanager->SetVerboseLevel(VerboseLevel);
      pmanager->AddProcess(theRadioactiveDecay);
      pmanager->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
      pmanager->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
    }

    if (particleName == "triton") {
      G4ProcessManager * pmanager =
          theIonTable->GetParticle(i)->GetProcessManager();
      pmanager->SetVerboseLevel(VerboseLevel);
      pmanager->AddProcess(theRadioactiveDecay);
      pmanager->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
      pmanager->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
    }
  }
}

// Cuts /////////////////////////////////////////////////////////////////////
void LSCPhysicsList::SetCuts()
{

  if (verboseLevel > 1) G4cout << "LSCPhysicsList::SetCuts:";
  if (verboseLevel > 0) {
    G4cout << "LSCPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << G4endl;
  }

  // special for low energy physics
  // G4double lowlimit=250*eV;
  // G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100.*GeV);

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  if (verboseLevel > 0) DumpCutValuesTable();
}

void LSCPhysicsList::SetNewValue(G4UIcommand * command, G4String newValue)
{
  if (command == gammaCutCmd)
    cutForGamma = gammaCutCmd->GetNewDoubleValue(newValue);
  else if (command == electCutCmd)
    cutForElectron = electCutCmd->GetNewDoubleValue(newValue);
  else if (command == posCutCmd)
    cutForPositron = posCutCmd->GetNewDoubleValue(newValue);
  else if (command == allCutCmd) {
    G4double cut = allCutCmd->GetNewDoubleValue(newValue);
    cutForGamma = cut;
    cutForElectron = cut;
    cutForPositron = cut;
  }
  else if (command == verboseCmd)
    VerboseLevel = verboseCmd->GetNewIntValue(newValue);
}