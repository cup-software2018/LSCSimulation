// This file is part of the GenericLAND software library.
// $Id: GLG4param.cc,v 1.2 2013/11/09 23:48:54 jslee Exp $
//
// GLG4param.cc
// defines a class for a "database" of numeric parameters
// (a string/value hashtable that is initialized from a file)
// Original by Glenn Horton-Smith, Jan 2001

#include "GLG4Sim/GLG4param.hh"

#include "fstream"
#include "sstream"

using namespace std;

GLG4param * GLG4param::theGLG4param= NULL;

// Constructor
GLG4param::GLG4param() {
  // nothing to do
}


// ReadFile
void
GLG4param::ReadFile(G4std::istream &is, EOverride oflag)
{
  while (is.good())
    {
      // get a line from the file
      char linebuffer[128];
      is.getline( linebuffer, sizeof(linebuffer)-1 );
      if ( is.fail() )
	break;

      // put the line in an istrstream for convenient parsing
      istringstream lineStream(linebuffer);

      // parse out name
      G4String name;
      G4double value;
      lineStream >> name;
      
      // skip lines beginning with '#' or '\n';
      if (lineStream.fail()
	  || name.length()==0 || (name.data())[0] == '\n'
	  || (name.data())[0]=='#' )
	continue;

      // parse out value
      lineStream >> value;
      if (lineStream.fail()) {
	G4cerr << "Warning: invalid/missing value encountered"
	  "in GLG4param::ReadFile(): line was ``" << linebuffer << "''" << G4endl;
	continue;
      }

      // does name already exist in hash?
      if (count(name)) {
	if ( oflag == kKeepExistingValue ) {
	  G4cout << "Info: GLG4param::ReadFile: retaining previous setting "
		 << name << "=" << (*this)[name]
		 << ", ignoring value " << value << " specified by file.\n";
	}
	else {
	  G4cout << "Info: GLG4param::ReadFile: OVERRIDING previous setting "
		 << name << "=" << (*this)[name]
		 << ", new value is " << value << " as specified by file.\n";
	  // insert name/value into hash
	  (*this)[name]= value;
	}
      }
      else {
	// insert name/value into hash
	(*this)[name]= value;
      }
    }
}

void
GLG4param::ReadFile(const char *filename, EOverride oflag)
{
  G4std::ifstream ifs;
  ifs.open(filename);
  if (!ifs.good()) {
    G4cout << "GLG4param::ReadFile : could not open " << filename << G4endl;
    return;
  }
  ReadFile(ifs, oflag);
  ifs.close();
}

void
GLG4param::WriteFile(G4std::ostream &os)
{
  iterator i;
  for (i=begin(); i!=end(); i++)
    os << (*i).first << '\t' << (*i).second << G4endl;
  os.flush();
}
