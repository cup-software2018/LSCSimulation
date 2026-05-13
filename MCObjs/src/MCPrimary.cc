#include <iostream>

#include "MCObjs/MCPrimary.hh"

ClassImp(MCPrimary)

MCPrimary::MCPrimary()
  : TObject()
{
}

MCPrimary::MCPrimary(const MCPrimary & prim)
  : TObject(),
    fName(prim.GetParticleName()),
    fPDGCode(prim.GetPDGCode()),
    fTrackId(prim.GetTrackId()),
    fKineticEnergy(prim.GetKineticEnergy()),
    fT0(prim.GetT0())
{
  prim.GetVertex(fVx, fVy, fVz);
  prim.GetMomentum(fPx, fPy, fPz);
}

MCPrimary::~MCPrimary() = default;

void MCPrimary::Print(Option_t * opt) const
{
  std::cout << "          particle : " << fName << std::endl;
  std::cout << "       Global Time : " << fT0 << " [ns]" << std::endl;
  std::cout << "    Kinetic Energy : " << fKineticEnergy << " [MeV]" << std::endl;
  std::cout << "          Momentum : " << Form("%7.2f", fPx) << " " << Form("%7.2f", fPy) << " "
            << Form("%7.2f", fPz) << " [MeV]" << std::endl;
  std::cout << "            Vertex : " << Form("%7.2f", fVx) << " " << Form("%7.2f", fVy) << " "
            << Form("%7.2f", fVz) << " [mm]" << std::endl;
}
