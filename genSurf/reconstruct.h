#ifndef RECONSTRUCT_H
#define RECONSTRUCT_H

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
#include "GL/glut.h"

#define PI 2.0 * asin(1.0)
#define D2R PI / 180.0
#define R2D 180.0 / PI
#define x2(x) (x * x)

unsigned seedRD = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 mySeed(seedRD);

plstream* pls; // pointer for plplotcxx

std::vector<std::vector<double>> pointVec; // store all points from all objects and stray points
std::vector<std::vector<double>> DELTA2;   // Dissimilarity matrix distance^2
std::vector<std::vector<double>> voxel;
std::vector<double> temp;

template <class T> class GLVector {
  public:
    T glX;
    T glY;
    T glZ;
};

double rad = 6.0;
double xx, yy, zz;
double dx = 0.0, dy = 0.0, dz = 0.0;
double epsa = 1.0e-3;
int totalPts = 0;
int xseg = 30, yseg = 30, zseg = 30;
int totalVoxel = xseg * yseg * zseg;
int icnt, jcnt, kcnt;

#endif // RECONSTRUCT_H
