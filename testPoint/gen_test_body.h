#ifndef GEN_TEST_BODY_H
#define GEN_TEST_BODY_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iterator>
#include <cstdio>
#include <cstring>

#include <random>
#include <chrono>

#include <plplot/plstream.h>

#define PI 2.0 * asin(1.0)
#define D2R PI / 180.0
#define R2D 180.0 / PI
#define x2(x) (x * x)

unsigned seedRD = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 mySeed(seedRD);

plstream* pls; // pointer for plplotcxx

std::vector<double> temp;

double rad = 6.0;
double centerX = 0.0, centerY = 0.0, centerZ = 0.0;
double xx, yy, zz;
double epsa = 1.0e-3;
double x1max = -999999.999, x1min = -x1max, y1max = x1max, y1min = x1min, z1max = x1max, z1min = x1min;
double x2max = -999999.999, x2min = -x2max, y2max = x2max, y2min = x2min, z2max = x2max, z2min = x2min;
double delO1O2 = 2.5 * rad; // delta between two objects;
const int numPoints = 1000;
const int strayPts = 1000;

std::uniform_real_distribution<double> distR(rad - 0.25, rad + 0.25);
std::uniform_real_distribution<double> distTht(-40.0, 30.0);
std::uniform_real_distribution<double> distPhi(-180.0, 180.0);
std::uniform_real_distribution<double> distLen(-10.0, 10.0);

void obj_Sphere(double& xx, double& yy, double& zz) {
  double rr, tht, phi;
  rr = distR(mySeed);
  tht = distTht(mySeed) * D2R;
  phi = distPhi(mySeed) * D2R;
  xx = rr * cos(tht) * cos(phi);
  yy = rr * cos(tht) * sin(phi);
  zz = rr * sin(tht);
}

void obj_Cylinder(double& xx, double& yy, double& zz) {
  double rr, phi;
  rr = distR(mySeed);
  phi = distPhi(mySeed) * D2R;
  xx = rr * cos(phi);
  yy = rr * sin(phi);
  zz = distLen(mySeed);
}


#endif // GEN_TEST_BODY_H
