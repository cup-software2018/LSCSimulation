// This file is part of the GenericLAND software library.
// $Id: GLG4PrimaryGeneratorMessenger.cc,v 1.2 2013/11/09 23:48:54 jslee Exp $
//
// GLG4PrimaryGeneratorMessenger.cc by Glenn Horton-Smith, Feb. 1999
// updated Aug. 3-17, 2001, for new GLG4PrimaryGeneratorAction

////////////////////////////////////////////////////////////////
// GLG4PrimaryGeneratorMessenger
////////////////////////////////////////////////////////////////

#include "GLG4Sim/GLG4PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4ParticleTable.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ios.hh"
#include "fstream" // for file streams
#include "globals.hh"
#include "iomanip" // for G4std::setw(), etc..
#include "sstream" // for string streams

#include "GLG4Sim/GLG4PosGen.hh" // for GetState(), SetState()
#include "GLG4Sim/GLG4PrimaryGeneratorAction.hh"
#include "GLG4Sim/GLG4VertexGen.hh" // for GetState(), SetState()
#include "GLG4Sim/local_g4compat.hh"

using namespace std;
using namespace CLHEP;

GLG4PrimaryGeneratorMessenger::GLG4PrimaryGeneratorMessenger(
    GLG4PrimaryGeneratorAction * myaction)
    : myAction(myaction)
{
  G4UIdirectory * GenDir = new G4UIdirectory("/generator/");
  GenDir->SetGuidance("Control the primary event generator.");

  ListCmd = new G4UIcommand("/generator/list", this);

  RateCmd = new G4UIcommand("/generator/rates", this);
  RateCmd->SetGuidance("Set/show the event rates (in Hz).");
  RateCmd->SetGuidance("  Choice: index of rate to be set (omit to show all)");
  RateCmd->SetGuidance("          [can also use full name instead of index]");
  RateCmd->SetGuidance("  Value : optional new value for rate (in Hz)");
  RateCmd->SetGuidance(
      "  Flag  : optional new value of \"pileup-only\" condition");
  RateCmd->SetGuidance(
      "          If flag == 1, then events of this type will only");
  RateCmd->SetGuidance(
      "          occur overlapping with other events (within ATWD frame)");

  G4UIparameter * aParam = new G4UIparameter('s');
  RateCmd->SetParameter(aParam);
  aParam->SetParameterName("choice");
  aParam->SetOmittable(true);

  aParam = new G4UIparameter('d');
  RateCmd->SetParameter(aParam);
  aParam->SetParameterName("value");
  aParam->SetParameterRange("value >= 0.0");
  aParam->SetOmittable(true);

  aParam = new G4UIparameter('d');
  RateCmd->SetParameter(aParam);
  aParam->SetParameterName("flag");
  aParam->SetOmittable(true);

  GunCmd = new G4UIcommand("/generator/gun", this);
  GunCmd->SetGuidance("Set/show gun parameters.");
  GunCmd->SetGuidance("  particle_name:           name of particle");
  GunCmd->SetGuidance("  x_mm, y_mm, z_mm:        position of gun in mm");
  GunCmd->SetGuidance("  px_MeV, py_MeV, pz_MeV:  momentum in MeV/c");
  GunCmd->SetGuidance("  K_MeV:                   kinetic energy override");
  GunCmd->SetGuidance("  pol_x, pol_y, pol_z:     polarization [optional]");
  GunCmd->SetGuidance("  mult:                    multiplicity [optional]");
  GunCmd->SetGuidance("For isotropic, leave px,py,pz zero and set K_MeV.");
  GunCmd->SetGuidance("For random polarization, leave pol_* zero.");
  GunCmd->SetGuidance("Note: this command is an alias for the two commands:\n"
                      "   /generator/pos/set 9 \"x y z\"\n"
                      "   /generator/vtx/set 17 \"particle_name px py pz K "
                      "polx poly polz mult\"\n"
                      "Use the above two commands directly for more control.\n"
                      "Type the */set commands with only one argument for "
                      "current state and help.");

  GunCmd->SetParameter(new G4UIparameter("particle_name", 's', true));
  GunCmd->SetParameter(new G4UIparameter("x_mm", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("y_mm", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("z_mm", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("px_MeV", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("py_MeV", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("pz_MeV", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("K_MeV", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("pol_x", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("pol_y", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("pol_z", 'd', true));
  GunCmd->SetParameter(new G4UIparameter("mult", 'd', true));

  VtxSetCmd = new G4UIcommand("/generator/vtx/set", this);
  VtxSetCmd->SetGuidance("Set vertex generator state");
  VtxSetCmd->SetGuidance("Usage: /generator/vtx/set (vertex_code) [setting]");
  VtxSetCmd->SetGuidance("where (vertex_code) is given in /generator/list");
  VtxSetCmd->SetGuidance("[use either the integer code or the name]");
  VtxSetCmd->SetGuidance("and setting is the state string (in quotes).");
  VtxSetCmd->SetGuidance("Omit [setting] for generator-specific help.");
  VtxSetCmd->SetParameter(new G4UIparameter("event_type", 's', false));
  VtxSetCmd->SetParameter(new G4UIparameter("setting", 's', true));

  PosSetCmd = new G4UIcommand("/generator/pos/set", this);
  PosSetCmd->SetGuidance("Set position generator state");
  PosSetCmd->SetGuidance("Usage: /generator/pos/set (pos_code) [setting]");
  PosSetCmd->SetGuidance("where (pos_code) is given by /generator/list");
  PosSetCmd->SetGuidance("[use either the integer code or the name]");
  PosSetCmd->SetGuidance("and setting is the state string (in quotes).");
  PosSetCmd->SetGuidance("Omit [setting] for generator-specific help.");
  PosSetCmd->SetParameter(new G4UIparameter("event_type", 's', false));
  PosSetCmd->SetParameter(new G4UIparameter("setting", 's', true));

  EventWindowCmd = new G4UIcommand("/generator/event_window", this);
  EventWindowCmd->SetGuidance("Set/show event time window");
  EventWindowCmd->SetParameter(new G4UIparameter("window(ns)", 'd', true));

  ChainClipCmd = new G4UIcommand("/generator/chain_clip", this);
  ChainClipCmd->SetGuidance("Set/show chain clipping time");
  ChainClipCmd->SetParameter(new G4UIparameter("chain_clip(ns)", 'd', true));
}

GLG4PrimaryGeneratorMessenger::~GLG4PrimaryGeneratorMessenger()
{
  delete RateCmd;
}

void GLG4PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,
                                                G4String newValues)
{
  if (command == ListCmd) {
    G4cout << "Generator codes and descriptions:\n"
           << " evt  pos  vtx  name (== posname:vtxname)\n";
    for (int i = 0; i < myAction->theNumEventTypes; i++) {
      G4cout << G4std::setw(4) << i << " " << G4std::setw(4)
             << myAction->GetPositionCodeForEventType(i) << " "
             << G4std::setw(4) << myAction->GetVertexCodeForEventType(i) << "  "
             << myAction->GetEventTypeName(i) << G4endl;
    }
  }
  else if (command == RateCmd) {
    istringstream is((const char *)newValues);
    int index = -1;
    double rate = -1.0;
    int flag;
    char id[128];
    is.get(id, sizeof(id), ' ');
    if (strlen(id) == 0) {
      // no argument, just print rate table and end
      G4cout << "Rates:\n";
      for (int i = 0; i < myAction->theNumEventTypes; i++) {
        G4cout << G4std::setw(3) << i << " " << G4std::setw(7)
               << myAction->GetEventRate(i) / (1. / second)
               << (myAction->GetEventTriggerCondition(i) ? " * " : "   ")
               << myAction->GetEventTypeName(i) << G4endl;
      }
      G4cout << " * indicates \"pileup-only\" generator condition set."
             << G4endl;
    }
    else {
      // have one or more arguments, parse them
      is.clear();
      is.seekg(0);
      is >> index;
      if (is.fail()) {
        // wasn't an integer, maybe a name.  Look for it. (The stupid way.)
        for (int i = 0; i < myAction->theNumEventTypes; i++) {
          if (strcmp(myAction->GetEventTypeName(i), id) == 0) {
            index = i;
            break;
          }
        }
      }
      if (index < 0 || index >= myAction->theNumEventTypes) {
        G4cout << "Index must be non-negative and less than "
               << myAction->theNumEventTypes
               << ", or must be a name exactly matching a defined generator."
               << "  (No match for \"" << id << "\" was found.)" << G4endl;
      }
      else {
        is >> rate;
        if (is.fail()) {
          G4cout << G4std::setw(3) << index << " " << G4std::setw(7)
                 << myAction->GetEventRate(index) / (1. / second) << "  "
                 << myAction->GetEventTypeName(index) << G4endl;
        }
        else if (rate < 0.0) {
          G4cout << "Rate must be >= zero\n";
        }
        else {
          myAction->SetEventRate(index, rate * (1. / second));
        }
        flag = myAction->GetEventTriggerCondition(index);
        is >> flag;
        if (!is.fail()) myAction->SetEventTriggerCondition(index, flag);
      }
    }
  }
  else if (command == GunCmd) {
    if (newValues.length() > 0)
      do { // this "do" is just to allow "break" to exit this block
        // get UI pointer for ApplyCommand
        G4UImanager * theUI = G4UImanager::GetUIpointer();
        // parse and set new values
        istringstream is((const char *)newValues);
        ostringstream os_vtx;
        ostringstream os_pos;

        // set particle type
        G4String newParticleName;
        is >> newParticleName;
        if (is.fail() || newParticleName.length() == 0) break;
        os_vtx << "/generator/vtx/set " << myAction->kGunVtxIndex << " \""
               << newParticleName << '\t';

        // set particle origin
        G4double x, y, z;
        is >> x >> y >> z;
        if (is.fail()) goto do_vtx;
        os_pos << "/generator/pos/set " << myAction->kGunPosIndex << " \"" << x
               << ' ' << y << ' ' << z << "\"" << G4std::ends;
        theUI->ApplyCommand(os_pos.str());
        os_pos.seekp(0);
        // os_pos.freeze(0); // prevents memory leak

        // set particle momentum direction
        is >> x >> y >> z;
        if (is.fail()) goto do_vtx;
        os_vtx << x << ' ' << y << ' ' << z << '\t';

        // optional override of kinetic energy
        G4double KE;
        is >> KE;
        if (is.fail()) goto do_vtx;
        os_vtx << KE << '\t';

        // set particle polarization
        is >> x >> y >> z;
        if (is.fail()) goto do_vtx;
        os_vtx << x << ' ' << y << ' ' << z << '\t';

        G4int mult;
        is >> mult;
        if (is.fail()) goto do_vtx;
        os_vtx << mult;

      do_vtx:
        os_vtx << "\"" << G4std::ends;
        theUI->ApplyCommand(os_vtx.str());
        os_vtx.seekp(0);
        // os_vtx.freeze(0); // prevents memory leak
      } while (false);
    // print current settings
    G4cout << "Current TestGun settings:\n\t" << GetCurrentValue(command)
           << G4endl;
  }
  else if (command == VtxSetCmd) {
    istringstream is((const char *)newValues);
    int index = -1;
    is >> index;
    if (is.fail()) {
      // wasn't an integer, maybe a name.  Look for it. (The stupid way.)
      char id[64];
      is.clear();
      is.get(id, sizeof(id), ' ');
      for (int i = 0; i < myAction->theNumVertexGenCodes; i++) {
        if (strcmp(myAction->GetVertexCodeName(i), id) == 0) {
          index = i;
          break;
        }
      }
    }
    if (index < 0 || index >= myAction->theNumVertexGenCodes) {
      G4cerr << "/generator/vtx/set: invalid index or name: arguments are \""
             << newValues << "\"" << G4endl;
      return;
    }
    myAction->GetVertexGenerator(index)->SetState(newValues.substr(is.tellg()));
  }
  else if (command == PosSetCmd) {
    istringstream is((const char *)newValues);
    int index = -1;
    is >> index;
    if (is.fail()) {
      // wasn't an integer, maybe a name.  Look for it. (The stupid way.)
      char id[64];
      is.clear();
      is.get(id, sizeof(id), ' ');
      for (int i = 0; i < myAction->theNumPosGenCodes; i++) {
        if (strcmp(myAction->GetPositionCodeName(i), id) == 0) {
          index = i;
          break;
        }
      }
    }
    if (index < 0 || index >= myAction->theNumPosGenCodes) {
      G4cerr << "/generator/pos/set: invalid index or name: arguments are \""
             << newValues << "\"" << G4endl;
      return;
    }
    myAction->GetPositionGenerator(index)->SetState(
        newValues.substr(is.tellg()));
  }
  else if (command == EventWindowCmd) {
    if (newValues.length() == 0) {
      G4cout << "Current window is " << myAction->GetEventWindow() / ns << " ns"
             << G4endl;
    }
    else {
      istringstream is((const char *)newValues);
      G4double newWindow = -1.0;
      is >> newWindow;
      if (is.fail() || newWindow <= 0.0) {
        G4cerr << "Time window must be positive" << G4endl;
      }
      else {
        myAction->SetEventWindow(newWindow * ns);
      }
    }
  }
  else if (command == ChainClipCmd) {
    if (newValues.length() == 0) {
      G4cout << "Current clip time is " << myAction->GetChainClip() / ns
             << " ns" << G4endl;
    }
    else {
      istringstream is((const char *)newValues);
      G4double newWindow = -1.0;
      is >> newWindow;
      if (is.fail() || newWindow <= 0.0) {
        G4cerr << "Time window must be positive" << G4endl;
      }
      else {
        myAction->SetChainClip(newWindow * ns);
      }
    }
  }
  else {
    G4cerr << "invalid GLG4 \"set\" command";
  }
}

G4String GLG4PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
{
  if (command == RateCmd) {
    ostringstream os;
    os << "Rates:\n";
    for (int i = 0; i < myAction->theNumEventTypes; i++) {
      os << G4std::setw(3) << i << " " << G4std::setw(7)
         << myAction->GetEventRate(i) / (1. / second)
         << (myAction->GetEventTriggerCondition(i) ? " * " : "   ")
         << myAction->GetEventTypeName(i) << G4endl;
    }
    os << G4std::ends;
    G4String rv(os.str());
    // os.freeze(0);
    return rv;
  }
  if (command == GunCmd) {
    G4String vstring =
        myAction->GetVertexGenerator(myAction->kGunVtxIndex)->GetState();
    G4String pstring =
        myAction->GetPositionGenerator(myAction->kGunPosIndex)->GetState();
    vstring.insert(vstring.find_first_of(" \t") + 1, pstring + "\t");
    return vstring;
  }
  if (command == VtxSetCmd) {
    ostringstream os;
    os << "Vertex generator settings:\n";
    for (int index = 0; index < myAction->theNumVertexGenCodes; index++)
      os << index << '\t' << myAction->GetVertexGenerator(index)->GetState()
         << '\n';
    os << G4std::ends;
    G4String rv(os.str());
    // os.freeze(0);
    return rv;
  }
  else if (command == PosSetCmd) {
    ostringstream os;
    os << "Position generator settings:\n";
    for (int index = 0; index < myAction->theNumPosGenCodes; index++)
      os << index << '\t' << myAction->GetPositionGenerator(index)->GetState()
         << '\n';
    os << G4std::ends;
    G4String rv(os.str());
    // os.freeze(0);
    return rv;
  }
  else if (command == EventWindowCmd) {
    ostringstream os;
    os << myAction->GetEventWindow() << G4std::ends;
    G4String rv(os.str());
    // os.freeze(0);
    return rv;
  }
  else if (command == ChainClipCmd) {
    ostringstream os;
    os << myAction->GetChainClip() << G4std::ends;
    G4String rv(os.str());
    // os.freeze(0);
    return rv;
  }

  return G4String("invalid GLG4PrimaryGeneratorMessenger \"get\" command");
}