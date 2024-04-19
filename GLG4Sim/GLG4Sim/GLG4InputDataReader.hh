// This file is part of the GenericLAND software library.
// $Id: GLG4InputDataReader.hh,v 1.2 2013/11/13 03:58:15 jslee Exp $
//
// GLG4InputDataReader.hh
// v.0 by Glenn Horton-Smith, Feb 12, 1999
// adapted for namespace support in Geant4.1, 4.2 -- 14-Jul-2000

#ifndef GLG4InputDataReader_hh
#define GLG4InputDataReader_hh

#include <iostream>

#include "globals.hh"
#include "local_g4compat.hh"

class GLG4InputDataReader {
public:

  class MyTokenizer {
  private:
    G4std::istream *isptr;
  public:
    MyTokenizer(G4std::istream &is) { isptr = &is; nval = 0.0;}
    
    enum { TT_EOF=-1, TT_STRING='a', TT_NUMBER='0' };
    
    int ttype;

    G4double nval;
    G4String sval;

    int nextToken(void);

    void dumpOn(G4std::ostream &os);
  };
  
  static int ReadMaterials(G4std::istream &is);
};
#endif
