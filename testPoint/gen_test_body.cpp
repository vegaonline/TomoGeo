#include "gen_test_body.h"

int main()
{

    // PLPLOT related routine jobs
    pls = new plstream();
    pls->init();

    const int NPTS = numPoints;
    PLFLT alt = 45.0;
    PLFLT az = 0.0;

    PLFLT* x1Plot = new PLFLT[NPTS]; //  X
    PLFLT* y1Plot = new PLFLT[NPTS]; //  Y
    PLFLT* z1Plot = new PLFLT[NPTS]; //  Z of Object 1

    PLFLT* x2Plot = new PLFLT[NPTS]; //  X
    PLFLT* y2Plot = new PLFLT[NPTS]; //  Y
    PLFLT* z2Plot = new PLFLT[NPTS]; //  Z of Object 2

    PLFLT* x3Plot = new PLFLT[strayPts]; //  X
    PLFLT* y3Plot = new PLFLT[strayPts]; //  Y
    PLFLT* z3Plot = new PLFLT[strayPts]; //  Z of Object 3

    // pls->dev("xwin");

    for(int ii = 0; ii < numPoints; ii++) {
        //  Object 1
        obj_Sphere(xx, yy, zz);
        xx -= 0.5 * del0102;
        x1max = std::max(x1max, xx);
        x1min = std::min(x1min, xx);
        y1max = std::max(y1max, yy);
        y1min = std::min(y1min, yy);
        z1max = std::max(z1max, zz);
        z1min = std::min(z1min, zz);
        x1Plot[ii] = xx;
        y1Plot[ii] = yy;
        z1Plot[ii] = zz;

        temp.push_back(xx);
        temp.push_back(yy);
        temp.push_back(zz);
        pointVec.push_back(temp);
        temp.clear();

        //  Object 2
        obj_Cylinder(xx, yy, zz);
        xx += 0.5 * del0102;
        x2max = std::max(x2max, xx);
        x2min = std::min(x2min, xx);
        y2max = std::max(y2max, yy);
        y2min = std::min(y2min, yy);
        z2max = std::max(z2max, zz);
        z2min = std::min(z2min, zz);
        x2Plot[ii] = xx;
        y2Plot[ii] = yy;
        z2Plot[ii] = zz;

        temp.push_back(xx);
        temp.push_back(yy);
        temp.push_back(zz);
        pointVec.push_back(temp);
    }


    double maxX = std::max(x1max + d2, x2max + d2);
    double minX = std::min(x1min - d2, x2min - d2);
    double maxY = std::max(y1max, y2max);
    double minY = std::min(y1min, y2min);
    double maxZ = std::max(z1max, z2max);
    double minZ = std::min(z1min, z2min);

    pls->env(3.0 * minX, 3.0 * maxX, 3.0 * minY, 3.0 * maxY, 0, 0);

    pls->col0(1);
    double baseX = 2.0 * std::max(std::abs(minX), maxX); // (xmax - xmin);
    double baseY = 2.0 * std::max(std::abs(minY), maxY); // (ymax - ymin);
    double Ht = (maxZ - minZ);

    // define the window
    pls->w3d(baseX, baseY, Ht, minX, maxX, minY, maxY, minZ, maxZ, alt, az);

    // Plot Object 1
    pls->col0(3);
    pls->string3(NPTS, x1Plot, y1Plot, z1Plot, ".");

    // plot Object 2
    pls->col0(4);
    pls->string3(NPTS, x2Plot, y2Plot, z2Plot, ".");

    // Now Plot some stray points around object 1 and 2
    std::uniform_real_distribution<double> strayX(minX, maxX);
    std::uniform_real_distribution<double> strayY(minY, maxY);
    std::uniform_real_distribution<double> strayZ(minZ, maxZ);i

    for(int nPts = 0; nPts < strayPts; nPts++) {
        x3Plot[nPts] = strayX(mySeed);
        y3Plot[nPts] = strayY(mySeed);
        z3Plot[nPts] = strayZ(mySeed);
        double xs = x3Plot[nPts];
        double ys = y3Plot[nPts];
        double zs = z3Plot[nPts];
        temp.push_back(xs);
        temp.push_back(ys);
        temp.push_back(zs);
        pointVec.push_back(temp);
    }
    pls->col0(5);
    pls->string3(NPTS, x3Plot, y3Plot, z3Plot, ".");

    std::ofstream filePOINT;
    filePOINT.open("points3D.dat");
    for(std::vector<std::vector<double>>::iterator itr1 = pointVec.begin(); itr1 != pointVec.end(); ++itr1) {
        for(std::vector<double>::iterator itr2 = itr1->begin(); itr2 != itr1->end(); ++itr2) {
            filePOINT << *itr2 << "  ";
        }
        filePOINT << std::endl;
    }

    filePOINT.close();
    delete[] x1Plot;
    delete[] y1Plot;
    delete[] z1Plot;
    delete[] x2Plot;
    delete[] y2Plot;
    delete[] z2Plot;
    delete[] x3Plot;
    delete[] y3Plot;
    delete[] z3Plot;
    pointVec.clear();
    delete pls;
    return 0;
}
