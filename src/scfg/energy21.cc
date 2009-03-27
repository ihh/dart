#include "scfg/energy21.h"

Duplex_lookup::Duplex_lookup()
  : duplex_idx (4, 4, (int) -1)
{
  // set up duplex index lookup tables
  int sym5_tmp[6] = { 0, 1, 2, 3, 2, 3 };
  int sym3_tmp[6] = { 3, 2, 1, 0, 3, 2 };
  duplex_sym5 = vector<int> (sym5_tmp, sym5_tmp + 6);
  duplex_sym3 = vector<int> (sym3_tmp, sym3_tmp + 6);
  for (int i = 0; i < 6; ++i) duplex_idx (duplex_sym5[i], duplex_sym3[i]) = i;
}

RNA_energy::RNA_energy()
  : stacked (4, 4, array2d<Score> (4, 4, -InfinityScore)),
    hairpin (4, 4, array2d<Score> (4, 4, -InfinityScore)),
    interior (4, 4, array2d<Score> (4, 4, -InfinityScore))
{ }

void RNA_energy::affine_fit (const vector<Score>& loop, int min_len, int max_len, Score& A_ret, Score& B_ret)
{
  double ysum = 0;
  double xsum = 0;
  double x2sum = 0;
  double xysum = 0;
  double n = 0;
  for (int x = min_len; x <= max_len; ++x)
    {
      const double y = loop[x];
      ysum += y;
      xsum += x;
      x2sum += x*x;
      xysum += x*y;
      n += 1;
    }
  const double Ey = ysum / n;
  const double Ex = xsum / n;
  const double Ex2 = x2sum / n;
  const double Exy = xysum / n;
  const double covxy = Exy - Ex*Ey;
  const double varx = Ex2 - Ex*Ex;
  B_ret = (Score) (covxy / varx);
  A_ret = (Score) (Ey - B_ret * Ex);
}

Energy21::Energy21() : RNA_energy()
{
  stacked(0,3)(0,3) = 1190;
  stacked(0,3)(1,2) = 2520;
  stacked(0,3)(2,1) = 2000;
  stacked(0,3)(2,3) = 660;
  stacked(0,3)(3,0) = 1150;
  stacked(0,3)(3,2) = 1150;

  stacked(1,2)(0,3) = 2250;
  stacked(1,2)(1,2) = 3380;
  stacked(1,2)(2,1) = 2310;
  stacked(1,2)(2,3) = 1230;
  stacked(1,2)(3,0) = 2000;
  stacked(1,2)(3,2) = 2460;

  stacked(2,1)(0,3) = 2870;
  stacked(2,1)(1,2) = 3960;
  stacked(2,1)(2,1) = 3380;
  stacked(2,1)(2,3) = 1600;
  stacked(2,1)(3,0) = 2520;
  stacked(2,1)(3,2) = 2490;

  stacked(2,3)(0,3) = 1230;
  stacked(2,3)(1,2) = 2490;
  stacked(2,3)(2,1) = 2460;
  stacked(2,3)(2,3) = 570;
  stacked(2,3)(3,0) = 1150;
  stacked(2,3)(3,2) = -1280;

  stacked(3,0)(0,3) = 1460;
  stacked(3,0)(1,2) = 2870;
  stacked(3,0)(2,1) = 2250;
  stacked(3,0)(2,3) = 970;
  stacked(3,0)(3,0) = 1190;
  stacked(3,0)(3,2) = 1230;

  stacked(3,2)(0,3) = 970;
  stacked(3,2)(1,2) = 1600;
  stacked(3,2)(2,1) = 1230;
  stacked(3,2)(2,3) = 460;
  stacked(3,2)(3,0) = 660;
  stacked(3,2)(3,2) = 570;

  hairpin(0,0)(0,3) = 1150;
  hairpin(0,0)(1,2) = 1310;
  hairpin(0,0)(2,1) = 1860;
  hairpin(0,0)(2,3) = 1510;
  hairpin(0,0)(3,0) = 980;
  hairpin(0,0)(3,2) = 980;

  hairpin(0,1)(0,3) = 890;
  hairpin(0,1)(1,2) = 1410;
  hairpin(0,1)(2,1) = 1220;
  hairpin(0,1)(2,3) = 1100;
  hairpin(0,1)(3,0) = 800;
  hairpin(0,1)(3,2) = 800;

  hairpin(0,2)(0,3) = 1900;
  hairpin(0,2)(1,2) = 2550;
  hairpin(0,2)(2,1) = 2480;
  hairpin(0,2)(2,3) = 2240;
  hairpin(0,2)(3,0) = 1600;
  hairpin(0,2)(3,2) = 1600;

  hairpin(0,3)(0,0) = 980;
  hairpin(0,3)(0,1) = 1260;
  hairpin(0,3)(0,2) = 1920;
  hairpin(0,3)(0,3) = 1260;
  hairpin(0,3)(1,0) = 800;
  hairpin(0,3)(1,1) = 790;
  hairpin(0,3)(1,2) = 790;
  hairpin(0,3)(1,3) = 790;
  hairpin(0,3)(2,0) = 1600;
  hairpin(0,3)(2,1) = 1300;
  hairpin(0,3)(2,2) = 1300;
  hairpin(0,3)(2,3) = 1300;
  hairpin(0,3)(3,0) = 930;
  hairpin(0,3)(3,1) = 930;
  hairpin(0,3)(3,2) = 930;
  hairpin(0,3)(3,3) = 930;

  hairpin(1,0)(0,3) = 1080;
  hairpin(1,0)(1,2) = 1690;
  hairpin(1,0)(2,1) = 2390;
  hairpin(1,0)(2,3) = 1740;
  hairpin(1,0)(3,0) = 1260;
  hairpin(1,0)(3,2) = 1260;

  hairpin(1,1)(0,3) = 830;
  hairpin(1,1)(1,2) = 730;
  hairpin(1,1)(2,1) = 1280;
  hairpin(1,1)(2,3) = 1100;
  hairpin(1,1)(3,0) = 790;
  hairpin(1,1)(3,2) = 790;

  hairpin(1,2)(0,0) = 1860;
  hairpin(1,2)(0,1) = 2390;
  hairpin(1,2)(0,2) = 2520;
  hairpin(1,2)(0,3) = 2330;
  hairpin(1,2)(1,0) = 1220;
  hairpin(1,2)(1,1) = 1280;
  hairpin(1,2)(1,2) = 1220;
  hairpin(1,2)(1,3) = 1100;
  hairpin(1,2)(2,0) = 2480;
  hairpin(1,2)(2,1) = 2290;
  hairpin(1,2)(2,2) = 1810;
  hairpin(1,2)(2,3) = 2290;
  hairpin(1,2)(3,0) = 1750;
  hairpin(1,2)(3,1) = 1800;
  hairpin(1,2)(3,2) = 1750;
  hairpin(1,2)(3,3) = 1580;

  hairpin(1,3)(0,3) = 640;
  hairpin(1,3)(1,2) = 1020;
  hairpin(1,3)(2,1) = 1800;
  hairpin(1,3)(2,3) = 1270;
  hairpin(1,3)(3,0) = 930;
  hairpin(1,3)(3,2) = 930;

  hairpin(2,0)(0,3) = 2170;
  hairpin(2,0)(1,2) = 2190;
  hairpin(2,0)(2,1) = 2520;
  hairpin(2,0)(2,3) = 2390;
  hairpin(2,0)(3,0) = 1920;
  hairpin(2,0)(3,2) = 1920;

  hairpin(2,1)(0,0) = 1310;
  hairpin(2,1)(0,1) = 1690;
  hairpin(2,1)(0,2) = 2190;
  hairpin(2,1)(0,3) = 1690;
  hairpin(2,1)(1,0) = 1410;
  hairpin(2,1)(1,1) = 730;
  hairpin(2,1)(1,2) = 730;
  hairpin(2,1)(1,3) = 680;
  hairpin(2,1)(2,0) = 2550;
  hairpin(2,1)(2,1) = 1800;
  hairpin(2,1)(2,2) = 1650;
  hairpin(2,1)(2,3) = 1800;
  hairpin(2,1)(3,0) = 1020;
  hairpin(2,1)(3,1) = 1020;
  hairpin(2,1)(3,2) = 1020;
  hairpin(2,1)(3,3) = 960;

  hairpin(2,2)(0,3) = 1600;
  hairpin(2,2)(1,2) = 1650;
  hairpin(2,2)(2,1) = 1810;
  hairpin(2,2)(2,3) = 1710;
  hairpin(2,2)(3,0) = 1300;
  hairpin(2,2)(3,2) = 1300;

  hairpin(2,3)(0,0) = 980;
  hairpin(2,3)(0,1) = 1260;
  hairpin(2,3)(0,2) = 1920;
  hairpin(2,3)(0,3) = 1260;
  hairpin(2,3)(1,0) = 800;
  hairpin(2,3)(1,1) = 790;
  hairpin(2,3)(1,2) = 790;
  hairpin(2,3)(1,3) = 790;
  hairpin(2,3)(2,0) = 1600;
  hairpin(2,3)(2,1) = 1300;
  hairpin(2,3)(2,2) = 1300;
  hairpin(2,3)(2,3) = 1300;
  hairpin(2,3)(3,0) = 930;
  hairpin(2,3)(3,1) = 930;
  hairpin(2,3)(3,2) = 930;
  hairpin(2,3)(3,3) = 930;

  hairpin(3,0)(0,0) = 1150;
  hairpin(3,0)(0,1) = 1080;
  hairpin(3,0)(0,2) = 2170;
  hairpin(3,0)(0,3) = 1160;
  hairpin(3,0)(1,0) = 890;
  hairpin(3,0)(1,1) = 830;
  hairpin(3,0)(1,2) = 390;
  hairpin(3,0)(1,3) = 570;
  hairpin(3,0)(2,0) = 1900;
  hairpin(3,0)(2,1) = 1200;
  hairpin(3,0)(2,2) = 1600;
  hairpin(3,0)(2,3) = 1200;
  hairpin(3,0)(3,0) = 430;
  hairpin(3,0)(3,1) = 640;
  hairpin(3,0)(3,2) = 430;
  hairpin(3,0)(3,3) = 550;

  hairpin(3,1)(0,3) = 570;
  hairpin(3,1)(1,2) = 680;
  hairpin(3,1)(2,1) = 1100;
  hairpin(3,1)(2,3) = 910;
  hairpin(3,1)(3,0) = 790;
  hairpin(3,1)(3,2) = 790;

  hairpin(3,2)(0,0) = 1510;
  hairpin(3,2)(0,1) = 1740;
  hairpin(3,2)(0,2) = 2390;
  hairpin(3,2)(0,3) = 1750;
  hairpin(3,2)(1,0) = 1100;
  hairpin(3,2)(1,1) = 1100;
  hairpin(3,2)(1,2) = 850;
  hairpin(3,2)(1,3) = 910;
  hairpin(3,2)(2,0) = 2240;
  hairpin(3,2)(2,1) = 1750;
  hairpin(3,2)(2,2) = 1710;
  hairpin(3,2)(2,3) = 1750;
  hairpin(3,2)(3,0) = 1140;
  hairpin(3,2)(3,1) = 1270;
  hairpin(3,2)(3,2) = 1140;
  hairpin(3,2)(3,3) = 1040;

  hairpin(3,3)(0,3) = 400;
  hairpin(3,3)(1,2) = 960;
  hairpin(3,3)(2,1) = 1580;
  hairpin(3,3)(2,3) = 1040;
  hairpin(3,3)(3,0) = 930;
  hairpin(3,3)(3,2) = 930;

  interior(0,0)(0,3) = 1150;
  interior(0,0)(1,2) = 1690;
  interior(0,0)(2,1) = 1950;
  interior(0,0)(2,3) = 1320;
  interior(0,0)(3,0) = 1170;
  interior(0,0)(3,2) = 1640;

  interior(0,1)(0,3) = 1170;
  interior(0,1)(1,2) = 1790;
  interior(0,1)(2,1) = 1690;
  interior(0,1)(2,3) = 1200;
  interior(0,1)(3,0) = 1080;
  interior(0,1)(3,2) = 1560;

  interior(0,2)(0,3) = 2280;
  interior(0,2)(1,2) = 2930;
  interior(0,2)(2,1) = 3050;
  interior(0,2)(2,3) = 2430;
  interior(0,2)(3,0) = 2260;
  interior(0,2)(3,2) = 2740;

  interior(0,3)(0,0) = 1170;
  interior(0,3)(0,1) = 1260;
  interior(0,3)(0,2) = 2400;
  interior(0,3)(0,3) = 780;
  interior(0,3)(1,0) = 1080;
  interior(0,3)(1,1) = 1070;
  interior(0,3)(1,2) = 310;
  interior(0,3)(1,3) = 1070;
  interior(0,3)(2,0) = 2260;
  interior(0,3)(2,1) = 830;
  interior(0,3)(2,2) = 1300;
  interior(0,3)(2,3) = 830;
  interior(0,3)(3,0) = 450;
  interior(0,3)(3,1) = 1120;
  interior(0,3)(3,2) = 450;
  interior(0,3)(3,3) = 2070;

  interior(1,0)(0,3) = 1270;
  interior(1,0)(1,2) = 1880;
  interior(1,0)(2,1) = 1910;
  interior(1,0)(2,3) = 1360;
  interior(1,0)(3,0) = 1260;
  interior(1,0)(3,2) = 1730;

  interior(1,1)(0,3) = 1210;
  interior(1,1)(1,2) = 1580;
  interior(1,1)(2,1) = 1650;
  interior(1,1)(2,3) = 1200;
  interior(1,1)(3,0) = 1070;
  interior(1,1)(3,2) = 1550;

  interior(1,2)(0,0) = 1950;
  interior(1,2)(0,1) = 1910;
  interior(1,2)(0,2) = 3090;
  interior(1,2)(0,3) = 2330;
  interior(1,2)(1,0) = 1690;
  interior(1,2)(1,1) = 1650;
  interior(1,2)(1,2) = 1220;
  interior(1,2)(1,3) = 1770;
  interior(1,2)(2,0) = 3050;
  interior(1,2)(2,1) = 2290;
  interior(1,2)(2,2) = 1910;
  interior(1,2)(2,3) = 2290;
  interior(1,2)(3,0) = 1750;
  interior(1,2)(3,1) = 1800;
  interior(1,2)(3,2) = 1750;
  interior(1,2)(3,3) = 2810;

  interior(1,3)(0,3) = 1020;
  interior(1,3)(1,2) = 1680;
  interior(1,3)(2,1) = 1800;
  interior(1,3)(2,3) = 1180;
  interior(1,3)(3,0) = 1120;
  interior(1,3)(3,2) = 1590;

  interior(2,0)(0,3) = 2550;
  interior(2,0)(1,2) = 2850;
  interior(2,0)(2,1) = 3090;
  interior(2,0)(2,3) = 2580;
  interior(2,0)(3,0) = 2400;
  interior(2,0)(3,2) = 2870;

  interior(2,1)(0,0) = 1690;
  interior(2,1)(0,1) = 1880;
  interior(2,1)(0,2) = 2850;
  interior(2,1)(0,3) = 1690;
  interior(2,1)(1,0) = 1790;
  interior(2,1)(1,1) = 1580;
  interior(2,1)(1,2) = 730;
  interior(2,1)(1,3) = 1620;
  interior(2,1)(2,0) = 2930;
  interior(2,1)(2,1) = 1800;
  interior(2,1)(2,2) = 1740;
  interior(2,1)(2,3) = 1800;
  interior(2,1)(3,0) = 1020;
  interior(2,1)(3,1) = 1680;
  interior(2,1)(3,2) = 1020;
  interior(2,1)(3,3) = 2670;

  interior(2,2)(0,3) = 1410;
  interior(2,2)(1,2) = 1740;
  interior(2,2)(2,1) = 1910;
  interior(2,2)(2,3) = 1420;
  interior(2,2)(3,0) = 1300;
  interior(2,2)(3,2) = 1780;

  interior(2,3)(0,0) = 1640;
  interior(2,3)(0,1) = 1730;
  interior(2,3)(0,2) = 2870;
  interior(2,3)(0,3) = 1540;
  interior(2,3)(1,0) = 1560;
  interior(2,3)(1,1) = 1550;
  interior(2,3)(1,2) = 690;
  interior(2,3)(1,3) = 1550;
  interior(2,3)(2,0) = 2740;
  interior(2,3)(2,1) = 1780;
  interior(2,3)(2,2) = 1780;
  interior(2,3)(2,3) = 1780;
  interior(2,3)(3,0) = 930;
  interior(2,3)(3,1) = 1590;
  interior(2,3)(3,2) = 930;
  interior(2,3)(3,3) = 2540;

  interior(3,0)(0,0) = 1150;
  interior(3,0)(0,1) = 1270;
  interior(3,0)(0,2) = 2550;
  interior(3,0)(0,3) = 680;
  interior(3,0)(1,0) = 1170;
  interior(3,0)(1,1) = 1210;
  interior(3,0)(1,2) = -90;
  interior(3,0)(1,3) = 1040;
  interior(3,0)(2,0) = 2280;
  interior(3,0)(2,1) = 730;
  interior(3,0)(2,2) = 1410;
  interior(3,0)(2,3) = 730;
  interior(3,0)(3,0) = -50;
  interior(3,0)(3,1) = 1020;
  interior(3,0)(3,2) = -50;
  interior(3,0)(3,3) = 1970;

  interior(3,1)(0,3) = 1040;
  interior(3,1)(1,2) = 1620;
  interior(3,1)(2,1) = 1770;
  interior(3,1)(2,3) = 1200;
  interior(3,1)(3,0) = 1070;
  interior(3,1)(3,2) = 1550;

  interior(3,2)(0,0) = 1320;
  interior(3,2)(0,1) = 1360;
  interior(3,2)(0,2) = 2580;
  interior(3,2)(0,3) = 800;
  interior(3,2)(1,0) = 1200;
  interior(3,2)(1,1) = 1200;
  interior(3,2)(1,2) = -0;
  interior(3,2)(1,3) = 1200;
  interior(3,2)(2,0) = 2430;
  interior(3,2)(2,1) = 800;
  interior(3,2)(2,2) = 1420;
  interior(3,2)(2,3) = 800;
  interior(3,2)(3,0) = 90;
  interior(3,2)(3,1) = 1180;
  interior(3,2)(3,2) = 90;
  interior(3,2)(3,3) = 2080;

  interior(3,3)(0,3) = 1820;
  interior(3,3)(1,2) = 2670;
  interior(3,3)(2,1) = 2810;
  interior(3,3)(2,3) = 2080;
  interior(3,3)(3,0) = 2070;
  interior(3,3)(3,2) = 2540;

  Score interior_loop_tmp[31] = { -InfinityScore, -InfinityScore, -3890, -4840, -4650, -5030, -5410, -5600, -5690, -5790, -5970, -6070, -6070, -6160, -6260, -6350, -6450, -6450, -6540, -6540, -6640, -6730, -6730, -6730, -6830, -6830, -6920, -6920, -7020, -7020, -7020 };
  interior_loop = vector<Score> (interior_loop_tmp, interior_loop_tmp + 31);
  
  Score bulge_loop_tmp[31] = { -InfinityScore, -3700, -2940, -3320, -3980, -4550, -4740, -4930, -5030, -5120, -5220, -5410, -5410, -5500, -5600, -5690, -5790, -5790, -5880, -5880, -5970, -5970, -6070, -6070, -6160, -6160, -6160, -6260, -6350, -6350, -6350 };
  bulge_loop = vector<Score> (bulge_loop_tmp, bulge_loop_tmp + 31);
  
  Score hairpin_loop_tmp[31] = { -InfinityScore, -InfinityScore, -InfinityScore, -3890, -4650, -4170, -4460, -4740, -4840, -4930, -5030, -5120, -5220, -5310, -5410, -5500, -5500, -5600, -5600, -5690, -5790, -5790, -5880, -5880, -5970, -5970, -5970, -6070, -6070, -6160, -6160 };
  hairpin_loop = vector<Score> (hairpin_loop_tmp, hairpin_loop_tmp + 31);
}
