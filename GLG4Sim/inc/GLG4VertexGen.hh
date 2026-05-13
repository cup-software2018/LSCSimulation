#ifndef __GLG4VertexGen_h__
#define __GLG4VertexGen_h__ 1
/** @file
 Declares GenericLAND global vertex generator classes for primary events,
 (See note on GenericLAND generators for more information.)

 This file is part of the GenericLAND software library.
 $Id: GLG4VertexGen.hh,v 1.1.1.1 2013/11/08 05:33:05 jslee Exp $

 @author G.Horton-Smith, August 3, 2001
*/

#include <set>     // for multiset
#include <stdio.h> // for FILE

#include "G4PrimaryVertex.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;
class G4Track;
// class G4PrimaryVertex;
class G4ParticleDefinition;
class GLG4PrimaryGeneratorAction;

/** Virtual base class for vertex generators */
class GLG4VVertexGen {
public:
  GLG4VVertexGen(const char * arg_dbname)
      : _dbname(arg_dbname)
  {
  }
  virtual ~GLG4VVertexGen() {}
  virtual void GeneratePrimaryVertex(G4Event * argEvent) = 0;
  // adds vertices.
  virtual void SetState(G4String newValues) = 0;
  // sets filename or other information needed by vertex generator
  virtual G4String GetState() = 0;
  // returns the current state information in a form that can be understood
  // by SetState (and, hopefully, a well-informed human)
protected:
  G4String _dbname; // used for GLG4param key prefix
};

/** vertex generator that can generate a primary vertex with one more
    particles of a given type, direction, energy, and polarization.
    Allows for randomly isotropic direction and random transverse polarization
    of spin-1, mass=0 particles */
class GLG4VertexGen_Gun : public GLG4VVertexGen {
public:
  GLG4VertexGen_Gun(const char * arg_dbname);
  virtual ~GLG4VertexGen_Gun();
  virtual void GeneratePrimaryVertex(G4Event * argEvent);
  // generates a primary vertex with given particle type, direction, energy,
  // and consistent polarization.
  virtual void SetState(G4String newValues);
  // format: particle_name  dir_x dir_y dir_z  kinetic_energy  polx poly polz
  // If dir_x==dir_y==dir_z==0, the directions are isotropic.
  // If particle has mass==0 and spin==1, final polarization will be
  // projected into plane perpindicular to momentum and made a unit vector;
  // if polarization has zero magnitude, a polarization is chosen randomly.
  virtual G4String GetState();
  // returns current state formatted as above

public:
  // the following useful static const data should be universally accessable
  // (I copied it from the G4IonTable source code, where it is privatized
  // with no accessor functions.)
  enum { numberOfElements = 110 };
  static const char * theElementNames[numberOfElements];

private:
  G4ParticleDefinition * _pDef;
  G4ThreeVector _mom;
  G4double _ke;
  G4ThreeVector _pol;
  G4int _multiplicity;
};

/** vertex generator that generates a primary vertex based on
    information from an ASCII stream.
    The ASCII stream contains information in a HEPEVT-like style.
*/
class GLG4VertexGen_HEPEvt : public GLG4VVertexGen {
public:
  GLG4VertexGen_HEPEvt(const char * arg_dbname);
  virtual ~GLG4VertexGen_HEPEvt();
  virtual void GeneratePrimaryVertex(G4Event * argEvent);
  // Generates a primary vertex based on information from an ASCII stream.
  // The ASCII stream contains information in a HEPEVT style.
  virtual void SetState(G4String newValues);
  // The argument is a Perl "open" style filename argument:
  // if the filename ends with '|', the filename is interpreted as
  // a command which pipes output to us; otherwise, if the filename begins
  // with '<' or an ordinary character, the file is opened for input.
  // The "pipe" style can be used to read a gzip'ped file or to start
  // a program which feeds events to us directly.
  virtual G4String GetState();
  // returns current state formatted as above

  void Open(const char * argFilename);
  void GetDataLine(char * buffer, size_t size);
  void Close();

  enum {
    kIonCodeOffset = 9800000,        // nuclei have codes like 98zzaaa
    kPDGcodeModulus = 10000000,      // PDG codes are 7 digits long
    kISTHEP_ParticleForTracking = 1, // only ISTHEP==1 are for tracking
    kISTHEP_InformatonMin = 100,     // 100 and above are "informatons"
    kISTHEP_Max = 213                // <= MAXINT / kPDGcodeModulus
  };

private:
  G4String _filename;
  FILE * _file;
  bool _isPipe;
};

class GLG4VertexGen_Stack : public GLG4VVertexGen {
public:
  GLG4VertexGen_Stack(const char * arg_dbname,
                      GLG4PrimaryGeneratorAction * argGLG4PGA);
  virtual ~GLG4VertexGen_Stack();
  virtual void GeneratePrimaryVertex(G4Event * argEvent);
  // Generates a primary vertex from stacked tracks.
  virtual void SetState(G4String newValues);
  // Does nothing.
  virtual G4String GetState();
  // returns ""

  void StackIt(const G4Track * track);

  struct LessThan_for_VertexTime {
    bool operator()(G4PrimaryVertex * a, G4PrimaryVertex * b) const
    {
      return (a->GetT0() < b->GetT0());
    }
  };

private:
  GLG4PrimaryGeneratorAction * fGLG4PGA;
  std::multiset<G4PrimaryVertex *, LessThan_for_VertexTime> _stack;
  typedef std::multiset<G4PrimaryVertex *, LessThan_for_VertexTime>::iterator
      _stack_iterator_t;
};

#endif
