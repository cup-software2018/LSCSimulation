// This file is part of the GenericLAND software library.
// $Id: GLG4param.hh,v 1.1.1.1 2013/11/08 05:33:05 jslee Exp $
//
// GLG4param.hh
// defines a class for a "database" of numeric parameters
// (a string/value hashtable that is initialized from a file)
// Original by Glenn Horton-Smith, Jan 2001
#ifndef __GLG4param_hh__
#define __GLG4param_hh__

#include "globals.hh"
#include "map"
#include "iostream"
#include "local_g4compat.hh"

class GLG4param : public G4std::map<G4String,G4double>
{
private:
  GLG4param(); // singleton
  static GLG4param * theGLG4param;
  
public:
  enum EOverride { kKeepExistingValue, kOverrideExistingValue };
  static inline GLG4param & GetDB();
  static inline GLG4param * GetDBPtr();
  void ReadFile(const char *filename, EOverride oflag= kKeepExistingValue);
  void ReadFile(G4std::istream &is, EOverride oflag= kKeepExistingValue);
  void WriteFile(G4std::ostream &os);
  inline G4double GetWithDefault(G4String name, G4double defaultValue);
};

// inline functions
GLG4param *
GLG4param::GetDBPtr()
{
  // Get (and possibly create) the pointer to the singleton
  //----------------
  if (theGLG4param == NULL)
    theGLG4param= new GLG4param();
  return theGLG4param;
}

GLG4param &
GLG4param::GetDB()
{
  // Get (and possibly create) the pointer to the singleton
  //----------------
  if (theGLG4param == NULL)
    theGLG4param= new GLG4param();
  return *theGLG4param;
}

G4double
GLG4param::GetWithDefault(G4String name, G4double defaultValue) {
  if (count(name))
    return (*this)[name];
  else
    return (*this)[name]= defaultValue;
}

#endif
