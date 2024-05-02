#include <iostream>
#include "MCObjs/MCPrimary.hh"

using namespace std;

ClassImp(MCPrimary)

MCPrimary::MCPrimary()
  : TObject()
{
  fName = "";
  fPDGCode = 0;
  fTrackId = 0;
  fVx = 0; fVy = 0; fVz = 0;
  fPx = 0; fPy = 0; fPz = 0;
  fKineticEnergy = 0;
  fT0 = 0;
}

MCPrimary::MCPrimary(const MCPrimary & prim)
  : TObject()
{
  fName = prim.GetParticleName();
  fPDGCode = prim.GetPDGCode();
  prim.GetVertex(fVx, fVy, fVz);
  prim.GetMomentum(fPx, fPy, fPz);
  fKineticEnergy = prim.GetKineticEnergy();
  fT0 = prim.GetT0();
}

MCPrimary::~MCPrimary()
{
}

void MCPrimary::Print(const Option_t * opt) const
{
  cout << "          particle : " << fName << endl;
  cout << "       Global Time : " << fT0 << " [ns]" << endl;
  cout << "    Kinetic Energy : " << fKineticEnergy << " [MeV]" << endl;
  cout << "          Momentum : "
       << Form("%7.2f", fPx) << " "
       << Form("%7.2f", fPy) << " "
       << Form("%7.2f", fPz) << " [MeV]" << endl;
  cout << "            Vertex : "
       << Form("%7.2f", fVx) << " "
       << Form("%7.2f", fVy) << " "
       << Form("%7.2f", fVz) << " [mm]" << endl;
}
