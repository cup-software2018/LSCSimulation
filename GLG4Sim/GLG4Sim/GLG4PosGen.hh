// This file is part of the GenericLAND software library.
// $Id: GLG4PosGen.hh,v 1.2 2013/11/09 23:49:10 jslee Exp $
//
// GenericLAND global position generator for primary events,
// by G.Horton-Smith, August 3, 2001
// (See note on GenericLAND generators for more information.)
#ifndef __GLG4PosGen_h__
#define __GLG4PosGen_h__ 1

#include "globals.hh"
#include "G4GeometryTolerance.hh"

class G4PrimaryVertex;
class G4VPhysicalVolume;
class G4Material;

#include <vector>

class GLG4VPosGen {
public:
  GLG4VPosGen(const char *arg_dbname) : _dbname(arg_dbname) { 
    kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  }
  virtual ~GLG4VPosGen() { }

  virtual void GenerateVertexPositions( G4PrimaryVertex *argVertex,
				double max_chain_time,
				double event_rate,
				double dt=0.0
				);
  // offsets the positions of already-filled vertices, "splitting" chains at
  // max_chain_time ("split" means a new random position is generated and
  // a random time offset added according to event_rate), and optionally
  // adding an extra time offset to each vertex.  Default uses GeneratePosition
  // to generate the offset for the chain (or each unsplit segment thereof).
  
  virtual void GeneratePosition( G4ThreeVector *argResult ) = 0;
  // generates random position, for usual case where there is no
  // dependence on particle types, momenta, etc.
  
  virtual void SetState( G4String newValues ) = 0;
  // sets filename or other information needed by global position generator
  
  virtual G4String GetState() = 0;
  // returns the current state information in a form that can be understood
  // by SetState (and, hopefully, a well-informed human)
  
  static void Strip( G4String &s, const char *stripchars=" \t\"");
  // strips leading and trailing characters from s
  
protected:
  G4String _dbname; // used for GLG4param key prefix
  G4double kCarTolerance;
};

class GLG4PosGen_null : public GLG4VPosGen {
public:
  GLG4PosGen_null(const char *arg_dbname) : GLG4VPosGen(arg_dbname) { }
  virtual void GeneratePosition( G4ThreeVector * /*argResult*/ ) { }
  // Does nothing.
  // Useful with a vertex generator which already sets positions.
  void SetState( G4String /*newValues*/ ) { }
  G4String GetState() { return "GLG4PosGen_null has no state"; }
};

class GLG4PosGen_PointPaintFill : public GLG4VPosGen {
public:
  GLG4PosGen_PointPaintFill(const char *arg_dbname);
  virtual void GeneratePosition( G4ThreeVector *argResult );
  // Generates a position either at a fixed point in the global coordinates
  // or uniformly filling the volume which contains the given point.
  // (This approach to specifying the volume isn't my favorite, but
  // it is strongly motivated by subtleties in Geant4's geometry code.)
  // - Fixed point is fastest.
  // - A random point in a compact physical volume is also pretty fast.
  // - A volume which only sparsely fills its geometric "extent" may
  //   require many iterations to find an internal point -- this will be slow.
  void SetState( G4String newValues );
  // newValues == x y z coordinates in mm (separated by white space),
  // optionally followed by keyword "fill" for volume-filling mode,
  // optionally followed by name of physical volume expected at that position;
  // or keyword "paint" for surface-painting mode,
  // optionally followed by name of physical volume expected at that position,
  // optionally followed by thickness of coat of "paint" (external to volume).
  // optionally followed by name of material to which to restrict "paint".
  G4String GetState();
  // returns current state in format above
private:
  enum { kPoint=0, kFill, kPaint };
  G4ThreeVector  _fixedPos;
  G4int _mode;
  G4double _thickness;
  G4String _pVolumeName;
  G4VPhysicalVolume* _pVolume;
  G4String _materialName;
  G4Material* _material;
  int _ntried;
  int _nfound;
  G4double _boundingBoxVolume;
  std::vector<G4ThreeVector> _intercepts;

};

class GLG4PosGen_Cosmic : public GLG4VPosGen {
public:
  GLG4PosGen_Cosmic(const char *arg_dbname);
  virtual void GenerateVertexPositions( G4PrimaryVertex *argVertex,
					double max_chain_time,
					double event_rate,
					double dt=0.0
					);
  // external flux uniformly distributed over area normal to
  // incident direction of first track in vertex
  void GeneratePosition(G4ThreeVector *);  // (not used)
  void SetState( G4String newValues );
  // newValues == "width height"
  //  width == width of rectangular area normal to incident direction (mm)
  //  height == height of rectangular area normal to incident direction (mm)
  // Appropriate values for GenericLAND would be 20000 33000.
  // The rectangle is rotated so that the "width" direction vector lies in
  // the XY plane for non-zero polar angle of the incident track.
  // If a track is generated which completely misses the detector,
  // the PDG code of the vertex tracks are modified to make them have
  // an impossible effective ISTHEP codes of 1, so Geant4 does not attempt
  // to track them.
  // This means the generated external flux in Hz/mm**2 is always
  //     flux = (rate/(width*height)),
  // regardless of the geometry of the detector, where "rate" is the rate
  // set via /generators/rate.  The "rate" must be chosen appropriately
  // for the area of the rectangle.
  G4String GetState();
  // returns current state in format above
private:
  G4double _width, _height;
};


#endif
