const double mm = 1;
const double cm = 10 * mm;
const double m = 100 * cm;

// 20-inch R12860
//const double pmt_radius = 50.8 / 2 * cm;
//const double pmt_height = 70 * cm;
//const double pmt_photo_r = 46.0 / 2 * cm;
//const double pmt_photo_h = 54 * cm;

// 10-inch R7801
const double pmt_radius = 25.3 / 2 * cm;
const double pmt_height = 34 * cm;
const double pmt_photo_r = 22.0 / 2 * cm;
const double pmt_photo_h = 24 * cm;

void plotgl()
{
  ifstream infile("pmtpos_prototype.dat");

  double buffer_radius = 2.2 / 2 * m;
  double buffer_height = 2.2 * m;

  TGeoManager * geom = new TGeoManager("nucleus", "Model of a nucleus");
  geom->SetNsegments(400);

  TGeoMaterial * matEmptySpace = new TGeoMaterial("EmptySpace", 0, 0, 0);
  TGeoMaterial * matPMT = new TGeoMaterial("PMT", .938, 1., 10000.);

  TGeoMedium * EmptySpace = new TGeoMedium("Empty", 1, matEmptySpace);
  TGeoMedium * PMT = new TGeoMedium("PMT", 1, matPMT);

  TGeoVolume * top = geom->MakeTube("WORLD", EmptySpace, 0, buffer_radius, buffer_height / 2);
  geom->SetTopVolume(top);

  TGeoVolume * pmt = geom->MakeSphere("pmt", PMT, 0., pmt_radius);
  pmt->SetLineColor(kRed);
  pmt->SetTransparency();

  int pmtno, ringno;
  double x, y, z;

  int n = 0;
  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    iss >> pmtno >> x >> y >> z >> ringno;
    top->AddNode(pmt, n, new TGeoTranslation(x, y, z));
    n += 1;
  }

  geom->CloseGeometry();
  geom->SetVisLevel(4);
  top->Draw("ogl");
}
