#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "GLG4Sim/GLG4InputDataReader.hh"
#include "LSCSim/LSCDetectorConstruction.hh"

using namespace CLHEP;

void LSCDetectorConstruction::ConstructMaterials()
{
  G4NistManager * man = G4NistManager::Instance();

  // === Common Elements ================================================
  G4bool isotope = false;
  G4Element * elementH = man->FindOrBuildElement("H", isotope);
  G4Element * elementC = man->FindOrBuildElement("C", isotope);
  G4Element * elementN = man->FindOrBuildElement("N", isotope);
  G4Element * elementO = man->FindOrBuildElement("O", isotope);
  G4Element * elementAl = man->FindOrBuildElement("Al", isotope);
  G4Element * elementSi = man->FindOrBuildElement("Si", isotope);
  G4Element * elementK = man->FindOrBuildElement("K", isotope);
  G4Element * elementCr = man->FindOrBuildElement("Cr", isotope);
  G4Element * elementFe = man->FindOrBuildElement("Fe", isotope);
  G4Element * elementNi = man->FindOrBuildElement("Ni", isotope);
  G4Element * elementNa = man->FindOrBuildElement("Na", isotope);
  G4Element * elementI = man->FindOrBuildElement("I", isotope);
  G4Element * elementCs = man->FindOrBuildElement("Cs", isotope);
  G4Element * elementCa = man->FindOrBuildElement("Ca", isotope);
  G4Element * elementF = man->FindOrBuildElement("F", isotope);
  G4Element * elementCu = man->FindOrBuildElement("Cu", isotope);
  G4Element * elementPb = man->FindOrBuildElement("Pb", isotope);
  G4Element * elementNb = man->FindOrBuildElement("Nb", isotope);
  G4Element * elementAu = man->FindOrBuildElement("Au", isotope);
  G4Element * elementGe = man->FindOrBuildElement("Ge", isotope);

  // ... insert any additional material definitions here
  G4String name;
  G4double density;
  G4double mol;
  G4int nelements;
  G4int natoms; // JW: comment out unused variable (2024.02.13.)
  G4MaterialPropertiesTable * MPT;

  // --- Air  N=70% O=30% ---------
  name = "Air";
  density = 1.29e-3 * g / cm3;
  nelements = 2;

  auto air = new G4Material(name, density, nelements);
  air->AddElement(elementN, 70 * perCent);
  air->AddElement(elementO, 30 * perCent);

  // --- PMT vacuum is very dilute air -------
  name = "PMT_Vac";
  density = 1e-3 * kGasThreshold;         // from PhysicalConstants.h
  G4double temperature = STP_Temperature; // from PhysicalConstants.h
  G4double pressure = STP_Pressure * density / (1.29e-3 * g / cm3);
  auto PMT_Vac =
      new G4Material(name, density, 1, kStateGas, temperature, pressure);
  PMT_Vac->AddMaterial(air, 1.);

  // --- Rock  SiO2 ---------------
  name = "Rock";
  density = 2.7 * g / cm3;
  nelements = 2;

  auto rock = new G4Material(name, density, nelements);
  rock->AddElement(elementSi, natoms = 1);
  rock->AddElement(elementO, natoms = 2);

  // --- Glass  SiO2 ---------------
  name = "Glass";
  density = 2.2 * g / cm3; // changed 1999/12/03 (was 2.7*g/cm3) -- GAS
  nelements = 2;

  auto glass = new G4Material(name, density, nelements);
  glass->AddElement(elementSi, natoms = 1);
  glass->AddElement(elementO, natoms = 2);

  // --- Iron  Fe ----------------
  name = "Iron";
  density = 7.87 * g / cm3;
  nelements = 1;

  auto steel = new G4Material(name, density, nelements);
  steel->AddElement(elementFe, natoms = 1);

  // --- Water  H2O ---------------
  name = "Water";
  density = 1.0 * g / cm3;
  nelements = 2;

  auto water = new G4Material(name, density, nelements);
  water->AddElement(elementH, natoms = 2);
  water->AddElement(elementO, natoms = 1);

  // --- Stainless Steel  71% Fe, 19% Cr, 10% Ni ------
  name = "Steel";
  density = 7.87 * g / cm3;
  nelements = 3;

  auto stainless = new G4Material(name, density, nelements);
  stainless->AddElement(elementFe, 0.71);
  stainless->AddElement(elementCr, 0.19);
  stainless->AddElement(elementNi, 0.10);

  // --- Lead  Pb ------
  name = "Lead";
  density = 11.35 * g / cm3;
  nelements = 1;

  auto lead = new G4Material(name, density, nelements);
  lead->AddElement(elementPb, natoms = 1);

  // --- Aluminum  Al ------
  name = "Aluminum";
  density = 2.7 * g / cm3;
  nelements = 1;

  auto aluminum = new G4Material(name, density, nelements);
  aluminum->AddElement(elementAl, natoms = 1);

  // --- Copper Cu ------
  name = "Copper";
  density = 8.96 * g / cm3;
  nelements = 1;

  auto copper = new G4Material(name, density, nelements);
  copper->AddElement(elementCu, natoms = 1);

  //               H H
  // --- Acrylic  -C-C- --------------------
  //               H COOCH3
  name = "Acrylic";
  density = 1.14 * g / cm3;
  nelements = 3;

  auto acrylic = new G4Material(name, density, nelements);
  acrylic->AddElement(elementH, natoms = 6);
  acrylic->AddElement(elementC, natoms = 4);
  acrylic->AddElement(elementO, natoms = 2);

  //
  // ............................. Grease SiO2 ...........................
  //
  name = "OptGrease";
  density = 1.06 * g / cm3;
  nelements = 2;
  auto grease = new G4Material(name, density, nelements);
  grease->AddElement(elementO, 2);
  grease->AddElement(elementSi, 1);

  //
  // ............................. Teflon PTFE ...........................
  //
  name = "Teflon";
  density = 2.2 * g / cm3;
  nelements = 2;
  auto teflon = new G4Material(name, density, nelements);
  teflon->AddElement(elementC, natoms = 2);
  teflon->AddElement(elementF, natoms = 4);

  // --- Polyethylene
  name = "Polyethylene";
  density = 0.91 * g / cm3;
  nelements = 2;

  auto polyethylene = new G4Material(name, density, nelements);
  polyethylene->AddElement(elementH, natoms = 2);
  polyethylene->AddElement(elementC, natoms = 1);

  // --- Tyvek  ==  High Density Polyethylene:  (...-CH2-CH2-...)*n
  name = "Tyvek";
  density = 0.96 * g / cm3;
  nelements = 2;

  auto tyvek = new G4Material(name, density, nelements);
  tyvek->AddElement(elementH, natoms = 2);
  tyvek->AddElement(elementC, natoms = 1);

  // photocathode material, approximated as elemental cesium
  name = "photocathode";
  density = 5. * g / cm3; // true??
  auto Photocathode_mat = new G4Material(name, density, nelements = 1);
  Photocathode_mat->AddElement(elementK, 1);

  // --- Mineral Oil  (CH2)n ------
  name = "MineralOil";
  density = 0.77 * g / cm3;
  nelements = 2;

  auto mineralOil = new G4Material(name, density, nelements);
  mineralOil->AddElement(elementC, natoms = 1);
  mineralOil->AddElement(elementH, natoms = 2);

  // Use the chemical formula as a useful label
  mineralOil->SetChemicalFormula("OIL");

  //
  // ............................. LAB ................................
  //
  // LAB (CnH2n+1-C6H5, n=9~14) //added by EJJeon (2008-02-26)
  G4int num_C;
  G4int num_H;
  char Name[15];
  G4Material * LAB[6];
  density = 0.86 * g / cm3;
  nelements = 2;
  for (int i = 0; i < 6; i++) {
    num_C = i + 15;
    num_H = 2 * (i + 9) + 6;
    sprintf(Name, "LAB_n=%i", i + 9);
    LAB[i] = new G4Material(Name, density, nelements);
    LAB[i]->AddElement(elementC, num_C);
    LAB[i]->AddElement(elementH, num_H);

    // Use the chemical formula as a label
    LAB[i]->SetChemicalFormula("AROMATIC");

    // Calculate the molecular weight
    mol = elementC->GetA() * num_C + elementH->GetA() * num_H;
    // Allocate memory for a new Material Property Table
    MPT = new G4MaterialPropertiesTable();
    // Fill with the molecular weight
    MPT->AddConstProperty("MOL", mol / g, true);
    // Attach this MPT to the pseudocumene
    LAB[i]->SetMaterialPropertiesTable(MPT);
  }

  //
  // ............................. PPO ...............................
  //
  // PPO (C15 H11 N 0) -- also called DPO, 2,5-diphenyloxazole
  density = 1.06 * g / cm3; // ??? at T=?
  auto PPO = new G4Material(name = "PPO", density, nelements = 4);
  PPO->AddElement(elementC, 15);
  PPO->AddElement(elementH, 11);
  PPO->AddElement(elementN, 1);
  PPO->AddElement(elementO, 1);

  // Use the chemical formula as a label
  PPO->SetChemicalFormula("FLUOR");

  // Calculate the molecular weight
  mol = elementC->GetA() * 15 + elementH->GetA() * 11 + elementN->GetA() * 1 +
        elementO->GetA() * 1;
  // Allocate memory for a new Material Property Table
  MPT = new G4MaterialPropertiesTable();
  // Fill with the molecular weight
  MPT->AddConstProperty("MOL", mol / g, true);
  // Attach this MPT to the PC
  PPO->SetMaterialPropertiesTable(MPT);

  //
  // .............................. Bis-MSB
  // .....................................
  //
  density = 1.3 * g / cm3; // Unknown
  nelements = 2;
  auto BisMSB = new G4Material("Bis-MSB", density, nelements);

  // Use the chemical formula as a label
  BisMSB->SetChemicalFormula("WLS");

  BisMSB->AddElement(elementC, 24);
  BisMSB->AddElement(elementH, 22);

  // Calculate the molecular weight
  mol = elementC->GetA() * 24 + elementH->GetA() * 22;
  // Allocate memory for a new Material Property Table
  MPT = new G4MaterialPropertiesTable();
  // Fill with the molecular weight
  MPT->AddConstProperty("MOL", mol / g, true);
  // Attach this MPT to the Bis-MSB
  BisMSB->SetMaterialPropertiesTable(MPT);

  // --- LAB-based Liquid Scintillator
  density = 0.86 * g / cm3;
  nelements = 8;
  auto LS_LAB = new G4Material(name = "LS_LAB", density, nelements);

  G4double PPO_fraction = 3 * g / (m3 * density);      // 3 g/l
  G4double BisMSB_fraction = 30 * mg / (m3 * density); // 30 mg/l

  LS_LAB->AddMaterial(LAB[0], 0.0047 / (1.0 + PPO_fraction + BisMSB_fraction));
  LS_LAB->AddMaterial(LAB[1], 0.097 / (1.0 + PPO_fraction + BisMSB_fraction));
  LS_LAB->AddMaterial(LAB[2], 0.3385 / (1.0 + PPO_fraction + BisMSB_fraction));
  LS_LAB->AddMaterial(LAB[3], 0.3472 / (1.0 + PPO_fraction + BisMSB_fraction));
  LS_LAB->AddMaterial(LAB[4], 0.2083 / (1.0 + PPO_fraction + BisMSB_fraction));
  LS_LAB->AddMaterial(LAB[5], 0.0043 / (1.0 + PPO_fraction + BisMSB_fraction));
  LS_LAB->AddMaterial(PPO,
                      PPO_fraction / (1.0 + PPO_fraction + BisMSB_fraction));
  LS_LAB->AddMaterial(BisMSB,
                      BisMSB_fraction / (1.0 + PPO_fraction + BisMSB_fraction));

  LS_LAB->GetIonisation()->SetBirksConstant(0.117 * mm / MeV);

  // Pyrene
  density = 1.271 * g / cm3;
  nelements = 2;
  auto Pyrene = new G4Material("Pyrene", density, nelements);

  Pyrene->SetChemicalFormula("FLUOR");
  Pyrene->AddElement(elementC, 16);
  Pyrene->AddElement(elementH, 10);

  mol = elementC->GetA() * 16 + elementH->GetA() * 12;
  MPT = new G4MaterialPropertiesTable();
  MPT->AddConstProperty("MOL", mol / g, true);
  Pyrene->SetMaterialPropertiesTable(MPT);

  density = 0.86 * g / cm3;
  nelements = 7;
  auto Pyrene_LS = new G4Material(name = "Pyrene_LS", density, nelements);

  G4double Pyrene_fraction = 10 * g / (m3 * density);
  Pyrene_LS->AddMaterial(LAB[0], 0.0047 / (1.0 + Pyrene_fraction));
  Pyrene_LS->AddMaterial(LAB[1], 0.097 / (1.0 + Pyrene_fraction));
  Pyrene_LS->AddMaterial(LAB[2], 0.3385 / (1.0 + Pyrene_fraction));
  Pyrene_LS->AddMaterial(LAB[3], 0.3472 / (1.0 + Pyrene_fraction));
  Pyrene_LS->AddMaterial(LAB[4], 0.2083 / (1.0 + Pyrene_fraction));
  Pyrene_LS->AddMaterial(LAB[5], 0.0043 / (1.0 + Pyrene_fraction));
  Pyrene_LS->AddMaterial(Pyrene, Pyrene_fraction / (1.0 + Pyrene_fraction));
  Pyrene_LS->GetIonisation()->SetBirksConstant(0.117 * mm / MeV);

  //***Material properties tables
  std::ifstream ifs;
  if (fMaterialDataFile.empty()) {
    G4String msg = "Error, material properties file could not be opened.\n";
    G4cerr << msg << G4endl;
    G4Exception("LSCDetectorConstruction::LSCDetectorConstruction", "",
                FatalException, msg);
  }
  else {
    ifs.open(fMaterialDataFile.data());
  }

  // now read materials, keeping error count
  int errorCount_ReadMaterials = GLG4InputDataReader::ReadMaterials(ifs);
  // close file
  ifs.close();

  if (errorCount_ReadMaterials) {
    G4cerr << "Error count after reading material properties file is "
           << errorCount_ReadMaterials << G4endl;
    G4String msg = "Error reading material properties file.\n";
    G4Exception("LSCDetectorConstruction::LSCDetectorConstruction", "",
                FatalException, msg);
  }

  //
  // == Create optical surfaces
  Photocathode_opsurf = new G4OpticalSurface("Photocathode_opsurf");
  Photocathode_opsurf->SetType(dielectric_metal); // ignored if RINDEX defined
  Photocathode_opsurf->SetMaterialPropertiesTable(
      G4Material::GetMaterial("photocathode")->GetMaterialPropertiesTable());

  Stainless_opsurf = new G4OpticalSurface("Stainless_opsurf");
  Stainless_opsurf->SetFinish(polished);
  Stainless_opsurf->SetModel(glisur);
  Stainless_opsurf->SetType(dielectric_metal);
  Stainless_opsurf->SetPolish(0.95); // a guess -- FIXME?
  Stainless_opsurf->SetMaterialPropertiesTable(
      stainless->GetMaterialPropertiesTable());

  Polyethylene_opsurf = new G4OpticalSurface("Polyethylene_opsurf");
  Polyethylene_opsurf->SetFinish(ground);              // a guess -- FIXME?
  Polyethylene_opsurf->SetModel(glisur);               // a guess -- FIXME?
  Polyethylene_opsurf->SetType(dielectric_dielectric); // a guess -- FIXME?
  Polyethylene_opsurf->SetPolish(0.7);                 // a guess -- FIXME?
  Polyethylene_opsurf->SetMaterialPropertiesTable(
      polyethylene->GetMaterialPropertiesTable());

  Tyvek_opsurf = new G4OpticalSurface("Tyvek_opsurf");
  Tyvek_opsurf->SetFinish(ground);
  Tyvek_opsurf->SetModel(glisur);
  Tyvek_opsurf->SetType(dielectric_metal);
  Tyvek_opsurf->SetPolish(0.01); // a guess -- FIXME
  Tyvek_opsurf->SetMaterialPropertiesTable(tyvek->GetMaterialPropertiesTable());

  Teflon_opsurf = new G4OpticalSurface("Teflon_opsurf");
  Teflon_opsurf->SetFinish(ground);
  Teflon_opsurf->SetModel(glisur);
  Teflon_opsurf->SetType(dielectric_metal);
  Teflon_opsurf->SetPolish(0.01); // a guess -- FIXME
  Teflon_opsurf->SetMaterialPropertiesTable(
      teflon->GetMaterialPropertiesTable());
}