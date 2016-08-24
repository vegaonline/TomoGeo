#include "kdnode.cpp"

int main()
{
    int npart= 0;
     std:;string line;
    std::ifstream filePOINT("points3D.dat");
    if (!filePOINT.eof()) {
	std::getline(filePOINT, line))
	std::stringstream ss(line);
	i = 0;
	while (ss >> temp[i])
            i++;
	pointVec.push_back(temp);
	++npart;
	temp.clear();
     }


    totalPts = npart;

    double maxX = std::max(x1max + d2, x2max + d2), minX = std::min(x1min - d2, x2min - d2);
    double maxY = std::max(y1max, y2max), minY = std::min(y1min, y2min);
    double maxZ = std::max(z1max, z2max), minZ = std::min(z1min, z2min);

    dx = (maxX - minX) / (xseg - 1);
    dy = (maxY - minY) / (yseg - 1);
    dz = (maxZ - minZ) / (zseg - 1);

    // Prepare Dissimilarity Matrix
    double deltaVal = 0.0;
    double dirCosl = 0.0, dirCosm = 0.0, dirCosn = 0.0, dirCos2 = 0.0;

    /*
        for(int ii = 0; ii < totalPts; ii++) {
            deltaVal = 0.0;
            std::vector<double> temp;
            dirCosl = pointVec[ii][0];
            dirCosm = pointVec[ii][1];
            dirCosn = pointVec[ii][2];
            dirCos2 = x2(dirCosl) + x2(dirCosm) + x2(dirCosn);
            // deltaVal = dirCosl

            DELTA2.push_back(temp);
            temp.clear();
        }
        std::cout << DELTA2.size() << std::endl;
        std::cout << " DONE..... You may exit..... " << std::endl;
    */

    // Generate voxel of dvx, dvy, dvz
    temp.clear();
    icnt = 0;
    jcnt = 0;
    kcnt = 0;
    for(int ii = 0; ii < totalVoxel; ii++) {
        if(!(ii % xseg)) {
            icnt = 0;
            jcnt++;
        }
        if(!(jcnt % yseg)) {
            icnt = 0;
            jcnt = 0;
            kcnt++;
        }

        temp.push_back(ii);
        temp.push_back(minX + icnt * dx);       // xlow
        temp.push_back(minX + (icnt + 1) * dx); // xhigh

        temp.push_back(minY + jcnt * dy);       // ylow
        temp.push_back(minY + (jcnt + 1) * dy); // yhigh

        temp.push_back(minZ + kcnt * dz);       // zlow
        temp.push_back(minZ + (kcnt + 1) * dz); // zhigh

        temp.push_back(0.0); // density for this voxel
        voxel.push_back(temp);
        temp.clear();
        icnt++;
    }

    // populate voxel with density

    int ilen = 0;

    for(std::vector<std::vector<double>>::iterator itr1 = voxel.begin(); itr1 != voxel.end(); ++itr1) {
        for(std::vector<double>::iterator itr2 = itr1->begin(); itr2 != itr1->end(); ++itr2) {
            if(ilen < 10)
                std::cout << *itr2 << "  ";
        }
        ilen++;
        std::cout << std::endl;
    }

    pointVec.clear();
    DELTA2.clear();
    voxel.clear();
    return 0;
}
