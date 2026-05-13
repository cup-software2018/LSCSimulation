
#ifndef LSCPhysicsList_h
#define LSCPhysicsList_h 1

#include "G4UImessenger.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcommand;

class LSCPhysicsList : public G4VModularPhysicsList, public G4UImessenger {
public:
  LSCPhysicsList();
  virtual ~LSCPhysicsList();

public:
  virtual void SetCuts();
  void SetNewValue(G4UIcommand *, G4String);

protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();

  // these methods Construct physics processes and register them
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  virtual void ConstructHad();
  virtual void ConstructOp();
  virtual void AddParameterisation();

private:
  G4bool omitHadronicProc;
  G4bool omitNeutHP;

  G4int VerboseLevel;
  G4int OpVerbLevel;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double cutForAll;

  // these methods Construct particles
  void ConstructMyBosons();
  void ConstructMyLeptons();
  void ConstructMyHadrons();
  void ConstructMyShortLiveds();

  G4UIdirectory * physDir;
  G4UIcmdWithADoubleAndUnit * gammaCutCmd;
  G4UIcmdWithADoubleAndUnit * electCutCmd;
  G4UIcmdWithADoubleAndUnit * posCutCmd;
  G4UIcmdWithADoubleAndUnit * allCutCmd;
  G4UIcmdWithAnInteger * verboseCmd;
  G4UIcmdWithAString * cherenkovOnCmd;
};

#endif