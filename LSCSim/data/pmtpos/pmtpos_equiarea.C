const double pi = TMath::Pi();

const double mm = 1;
const double cm = 10 * mm;
const double m = 100 * cm;

// 20-inch R12860
const double pmt_radius = 50.8 / 2 * cm;
const double pmt_height = 68 * cm;
const double pmt_photo_r = 46.0 / 2 * cm;
const double pmt_photo_h = 49.2 * cm;

// 10-inch R7801
// const double pmt_radius = 25.3 / 2 * cm;
// const double pmt_height = 32 * cm;
// const double pmt_photo_r = 22.0 / 2 * cm;
// const double pmt_photo_h = 22.33 * cm;

struct pmtpos {
  int id;
  double x, y, z;
  int ring;
  int where;
};

void pmtpos_equiarea()
{
  // real buffer dimension
  double buffer_radius = 17. / 2 * m;
  double buffer_height = 17. * m;
  // double buffer_radius = 2.2 / 2 * m;
  // double buffer_height = 2.2 * m;

  double buffer_offset_r = 5. * cm;
  double buffer_offset_h = 5. * cm;
  double pmt_r_offset = 1. * cm;

  // virtual photo cylinder
  double photo_cylinder_r = buffer_radius - buffer_offset_r - pmt_photo_h;
  double photo_cylinder_h = buffer_height - 2 * buffer_offset_h - 2 * pmt_photo_h;

  double pmt_area = pi * pmt_radius * pmt_radius;
  double photo_top_area = pi * photo_cylinder_r * photo_cylinder_r;

  int npmt_init = photo_top_area / pmt_area;
  double unit_area_top = photo_top_area / npmt_init;
  double unit_length_top = TMath::Sqrt(unit_area_top);

  if (unit_length_top - 2 * pmt_radius < 0) {
    unit_length_top = 2 * (pmt_radius + pmt_r_offset);
    unit_area_top = unit_length_top * unit_length_top;
  }

  int nring = 2 * photo_cylinder_r / unit_length_top;
  nring = nring % 2 == 0 ? nring - 1 : nring;

  double unit_d = 2 * photo_cylinder_r / nring;
  double unit_r = unit_d / 2;
  nring = nring / 2 + 1;

  cout << Form("number of ring=%d (ur=%f) for top and bottom", nring, unit_d) << endl;

  int ring = 0;
  int npmt = 0;
  vector<pmtpos *> pmts;

  double top_z = buffer_height / 2 - buffer_offset_h - pmt_height / 2;

  int npmt_top = 0;
  for (int i = 1; i < nring; i++) {
    ring += 1;
    double r = i * unit_d;
    double l = 2 * pi * r;
    int np = l / unit_length_top;
    for (int j = 0; j < np; j++) {
      npmt += 1;
      npmt_top += 1;
      double x = r * TMath::Cos(j * 2 * pi / np);
      double y = r * TMath::Sin(j * 2 * pi / np);
      pmtpos * pmt = new pmtpos;
      pmt->id = npmt;
      pmt->x = x;
      pmt->y = y;
      pmt->z = top_z;
      pmt->ring = ring;
      pmt->where = 1;
      pmts.push_back(pmt);
    }
  }

  //
  // barrel
  //
  double top_barrel_offset = pmt_height - pmt_photo_h;
  double photo_cylinder_h_temp = photo_cylinder_h - 2 * top_barrel_offset;

  double photo_cylinder_l = 2 * pi * photo_cylinder_r;
  double photo_barrel_area = photo_cylinder_l * photo_cylinder_h_temp;

  npmt_init = photo_barrel_area / pmt_area;
  double unit_area_barrel = photo_barrel_area / npmt_init;
  double unit_length_barrel = TMath::Sqrt(unit_area_barrel);

  if (unit_length_barrel - 2 * pmt_radius < 0) {
    unit_length_barrel = 2 * (pmt_radius + pmt_r_offset);
    unit_area_barrel = unit_length_barrel * unit_length_barrel;
  }

  int nx = photo_cylinder_l / unit_length_barrel;
  int ny = photo_cylinder_h_temp / unit_length_barrel;

  double unit_dx = photo_cylinder_l / nx;
  double unit_dy = photo_cylinder_h_temp / ny;

  cout << Form("barrel x:%d (%f) y:%d (%f)", nx, unit_dx, ny, unit_dy) << endl;

  // barrel pmts position
  double barrel_r = buffer_radius - buffer_offset_r - pmt_height / 2;

  int npmt_barrel = 0;
  for (int j = 0; j < ny; j++) {
    ring += 1;
    for (int i = 0; i < nx; i++) {
      npmt += 1;
      npmt_barrel += 1;
      double x = barrel_r * TMath::Cos(i * 2 * pi / nx);
      double y = barrel_r * TMath::Sin(i * 2 * pi / nx);
      double z = photo_cylinder_h_temp / 2 - (j + 0.5) * unit_dy;
      pmtpos * pmt = new pmtpos;
      pmt->id = npmt;
      pmt->x = x;
      pmt->y = y;
      pmt->z = z;
      pmt->ring = ring;
      pmt->where = 0;
      pmts.push_back(pmt);
    }
  }

  //
  // bottom
  //
  int npmt_bottom = 0;
  for (int i = 0; i < nring; i++) {
    ring += 1;
    double r = (nring - 1 - i) * unit_d;
    if (r == 0) {
      npmt += 1;
      npmt_bottom += 1;
      double x = 0;
      double y = 0;
      pmtpos * pmt = new pmtpos;
      pmt->id = npmt;
      pmt->x = x;
      pmt->y = y;
      pmt->z = -top_z;
      pmt->ring = ring;
      pmt->where = -1;
      pmts.push_back(pmt);
      continue;
    }
    double l = 2 * pi * r;
    int np = l / unit_length_top;
    for (int j = 0; j < np; j++) {
      npmt += 1;
      npmt_bottom += 1;
      double x = r * TMath::Cos(j * 2 * pi / np);
      double y = r * TMath::Sin(j * 2 * pi / np);
      pmtpos * pmt = new pmtpos;
      pmt->id = npmt;
      pmt->x = x;
      pmt->y = y;
      pmt->z = -top_z;
      pmt->ring = ring;
      pmt->where = -1;
      pmts.push_back(pmt);
    }
  }

  double pmt_photo_area = pi * pmt_photo_r * pmt_photo_r;
  double photo_total_area = 2 * photo_top_area + photo_barrel_area;

  cout << "- number of top pmt = " << npmt_top << endl;
  cout << "- number of bottom pmt = " << npmt_bottom << endl;
  cout << "- number of barrel pmt = " << npmt_barrel << endl;
  cout << "- total number of pmt = " << npmt << endl;
  cout << "- photo-coverage = " << npmt * pmt_photo_area / photo_total_area << endl;

  // positon file out
  ofstream out("pmtpos.dat");
  for (auto pmt : pmts) {
    out << Form("%5d %10.2f %10.2f %10.2f %5d %5d", pmt->id, pmt->x, pmt->y, pmt->z, pmt->ring, pmt->where) << endl;
  }
}
