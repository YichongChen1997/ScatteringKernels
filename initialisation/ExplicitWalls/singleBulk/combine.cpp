#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <time.h>
#include <cstdlib>

using namespace std;
const float PI = 3.1415;

double fRand(double fMin, double fMax);

//DATA: https://ww2.chemistry.gatech.edu/~lw26/structure/small_molecules/index.html

int main()
{
    string line;
    double x, y, z, xLo, xHi, yLo, yHi, zLo, zHi;
    double rCut, thickness, H, mG, sigmaG, mW, rho, zDynamicRange;
    int atoms, typ, i, j, k, c, m, id, mol = 0, molId = 0, nWallCount = 0, co;
    int nAtomsTot = 0;

    // reference scales (real to SI unit)
    double refLength = 1e-10;
    double refForce = 6.9477E-11;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime;
    double kB = 1.38064852e-23;
    double zRigidRange = 3; // rigid bound of the wall

    //***********************************
    int Nlines = 7;         // find out the No. of lines of parameters    "wc -l < Parameters.dat"
    ifstream Parameters("Parameters.dat", ios::in);
    for (int n = 1; n < Nlines + 1; n++)
    {
        if (n == 1)  {Parameters >> xLo >> xHi;                                  getline(Parameters, line);}
        if (n == 2)  {Parameters >> yLo >> yHi;                                  getline(Parameters, line);}
        if (n == 3)  {Parameters >> H >> rCut >> thickness;                      getline(Parameters, line);}
        if (n == 4)  {Parameters >> rho;                                         getline(Parameters, line);}
        if (n == 5)  {Parameters >> sigmaG;                                      getline(Parameters, line);}
        if (n == 6)  {Parameters >> mW;                                          getline(Parameters, line);}
        if (n == 7)  {Parameters >> mG;                                          getline(Parameters, line);}
    }

    srand(time(NULL));

    double Lx = xHi - xLo;
    double Ly = yHi - yLo;
    zLo = -thickness - zRigidRange;
    zHi = H + 2 * thickness + zRigidRange;

    // create only the bottom wall
    cout << "read in unit cell coordinates" << endl;

    ifstream data_unitCell("unitCell.dat", ios::in);
    getline(data_unitCell, line);

    int nSites;
    double unitCellSize;
    data_unitCell >> nSites;
    getline(data_unitCell, line);
    getline(data_unitCell, line);

    cout << "nSites = " << nSites << endl;

    data_unitCell >> unitCellSize;
    getline(data_unitCell, line);
    getline(data_unitCell, line);

    cout << "unitCellSize = " << unitCellSize << endl;

    vector<double> unitCellX(nSites);
    vector<double> unitCellY(nSites);
    vector<double> unitCellZ(nSites);

    for (i = 0; i < nSites; i++)
    {
        data_unitCell >> x >> y >> z;

        cout << x << " " << y << " " << z << endl;
        unitCellX[i] = x * unitCellSize;
        unitCellY[i] = y * unitCellSize;
        unitCellZ[i] = z * unitCellSize;
    }

    int nX = ceil(Lx / unitCellSize);
    int nY = ceil(Ly / unitCellSize);
    int nZ = ceil(2 * thickness / unitCellSize);
    cout << "------------------------------" << endl;
    cout << "No. of cells in each direction" << endl;
    cout << "nX = " << nX << endl;
    cout << "nY = " << nY << endl;
    cout << "nZ = " << nZ << endl;

    // read in the artificial surface
    cout << "------------------------------" << endl;
    cout << "read in the artificial surface" << endl;

    ifstream data_roughSurface("roughSurface.dat", ios::in);
    vector<vector<double> > roughSurface;

    for (j = 0; j < nY; j++)
        {
            vector<double> roughTemp(nX, 0.0);
            for (i = 0; i < nX; i++)
            {
                data_roughSurface >> roughTemp[i];
                // roughTemp[i] = roughTemp[i] / refLength;
            }
            roughSurface.push_back(roughTemp);
        }

    // define new Lx and Lz, and box geometry
    Lx = nX * unitCellSize;
    xLo = 0.0;
    xHi = Lx;

    Ly = nY * unitCellSize;
    yLo = 0.0;
    yHi = Ly;

    cout << "change of box co-ordinates" << endl;
    cout << "xLo = " << xLo << endl;
    cout << "xHi = " << xHi << endl;
    cout << "yLo = " << yLo << endl;
    cout << "yHi = " << yHi << endl;

    vector<double> xWall;
    vector<double> yWall;
    vector<double> zWall;

    double rIx, rIy, rIz, xP, yP, zP;
    double smallOffset = 0.2;

    double zMin = 0.0;

    for (k = 0; k < nZ; k++)
    {
        for (j = 0; j < nY; j++)
        {
            for (i = 0; i < nX; i++)
            {
                for (m = 0; m < nSites; m++)
                {
                    rIx = unitCellX[m];
                    rIy = unitCellY[m];
                    rIz = unitCellZ[m];

                    xP = rIx + unitCellSize * i + xLo + smallOffset;
                    yP = rIy + unitCellSize * j + yLo + smallOffset;
                    zP = thickness - (rIz + unitCellSize * k);

                    if (zP <= roughSurface[i][j])    
                    {
                        xWall.push_back(xP);
                        yWall.push_back(yP);
                        zWall.push_back(zP);

                        if (zP < zMin)
                        {
                            zMin = zP;
                        }

                    }
                }
            }
        }
    }

    int nWall = xWall.size();

    vector<double> xWall2;
    vector<double> yWall2;
    vector<double> zWall2;

    double zMax = 0.0;

    for (i = 0; i < nWall; i++)
    {
        z = -zWall[i] + H + thickness;
        xWall2.push_back(xWall[i]);
        yWall2.push_back(yWall[i]);
        zWall2.push_back(z);

        if (z > zMax)
        {
            zMax = z;
        }
        if (z < zMin)
        {
            zMin = z;
        }

        z = zWall[i];

        if (z > zMax)
        {
            zMax = z;
        }
        if (z < zMin)
        {
            zMin = z;
        }
    }

    cout << "zMin = " << zMin << endl;
    cout << "zMax = " << zMax << endl;

    nAtomsTot = nWall * 2;


    // create gas molecules
    double volume = (xHi - xLo) * (yHi - yLo) * H * refLength * refLength * refLength;

    int nGasPr = rho * volume;
    double offsetX = 3;
    double offsetY = 3;
    double offsetZ = 3;
    int nPx, nPy, nPz, nP;

    vector<double> xG(nGasPr);
    vector<double> yG(nGasPr);
    vector<double> zG(nGasPr);

    // RiGHT
    double rhoR = nGasPr / (Lx * Ly * H);
    double Nbox = rhoR * Lx * Lx * Lx;
    int nB = int(pow(Nbox, 1.0 / 3));

    cout << "nB = " << nB << endl;

    nPx = int(((Lx - offsetX) * nB / Lx) + 1);
    nPy = int(((Ly - offsetY) * nB / Lx) + 1);
    nPz = int(((H -  thickness - offsetZ) * nB / Lx) + 1);
    nP = 0;
    bool search = true;
    while (search)
    {
        nP = nPx * nPy * nPz;
        if (nP > nGasPr)
        {
            search = false;
        }
        else
        {
            nPx++;
            nPz++;
            nPy++;
        }
    }

    double uCX = (Lx - offsetX) / double(nPx);
    double uCY = (Ly - offsetY) / double(nPy);
    double uCZ = (H - thickness - 2 * offsetZ) /  double(nPz);

    cout << " nPx = " << nPx << endl
         << " nPy = " << nPy << endl
         << " nPz = " << nPz << endl;

    cout << " uCX = " << uCX << endl
         << " uCY = " << uCY << endl
         << " uCZ = " << uCZ << endl;

    co = 0;

    for (j = 0; j < nPy; j++)
    {
        for (k = 0; k < nPz; k++)
        {
            for (i = 0; i < nPx; i++)
            {
                if (co < nGasPr)
                {
                    xG[co] = xLo + i * uCX + offsetX;
                    yG[co] = yLo + j * uCY + offsetY;
                    zG[co] = thickness + k * uCZ + offsetZ;
                    co++;
                }
            }
        }
    }

    nAtomsTot += co;

    cout << "------------------------------------"  << endl
         << "nGas  mols  to insert = " << nGasPr    << endl
         << "nWall atoms to insert = " << nWall     << endl
         << "Total atoms to insert = " << nAtomsTot << endl
         << "------------------------------------"  << endl;

    ofstream writer("data.dat");

    writer << "LAMMPS data file" << endl;
    writer << endl;
    writer << nAtomsTot << '\t' << "atoms" << endl;
    writer << endl;
    writer << "4	atom types" << endl;
    writer << endl;
    writer << xLo << '\t' << xHi << "	xlo xhi" << endl;
    writer << yLo << '\t' << yHi << "	ylo yhi" << endl;
    writer << zLo << '\t' << zHi << "	zlo zhi" << endl;

    writer << endl;
    writer << "Masses" << endl;
    writer << endl;
    writer << "1 " << mG << endl;   // gas
    writer << "2 " << mW << endl;   // wall bottom
    writer << "3 " << mW << endl;   // wall bottom
    writer << "4 " << mW << endl;   // constrained wall
    writer << endl;
    writer << "Atoms" << endl;
    writer << endl;

    typ = 1;
    for (i = 0; i < nGasPr; i++)
    {
        molId++;
        mol++;
        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xG[i] << '\t' << yG[i] << '\t' << zG[i] << endl;
    }

    cout << "Number of GasPr atoms placed = " << mol << endl;
    double refNdenstiy = mol / (Lx * Ly * H * pow(refLength, 3));

    int nConstWall = 0;

    for (i = 0; i < nWall; i++)
    {
        typ = 2; // on purpose here
        molId++;
        mol++;

        if ((zWall[i] >= zMax - zRigidRange) || (zWall[i] <= zMin + zRigidRange))
        {
            typ = 4;
            nConstWall++;
        }

        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xWall[i] << '\t' << yWall[i] << '\t' << zWall[i] << endl;
    }

    for (i = 0; i < nWall; i++)
    {
        typ = 3; // on purpose here
        molId++;
        mol++;

        if ((zWall2[i] >= zMax - zRigidRange) || (zWall2[i] <= zMin + zRigidRange))
        {
            typ = 4;
            nConstWall++;
        }

        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xWall2[i] << '\t' << yWall2[i] << '\t' << zWall2[i] << endl;
    }

    cout << "Number of total atoms placed = " << mol << endl
        << "------------------------------------"  << endl;

    writer << endl;
    writer.close();


    double redNdensity = refNdenstiy * (PI * pow(sigmaG, 3) * pow(refLength, 3) / 6);
    double confRatio = H / sigmaG;
    double lambda = 1 / (sqrt(2) * PI * pow(sigmaG, 2) * pow(refLength, 2) * refNdenstiy);
    double Kn = lambda / (H * refLength);

    cout << "Reference number density = " << refNdenstiy << endl;
    cout << "Reduced number density   = " << redNdensity << endl;
    cout << "Confinemen ratio R       = " << confRatio << endl;
    cout << "Mean free path           = " << lambda << endl;
    cout << "Knudsen number (Kn)      = " << Kn << endl;
    cout << "------------------------------------"  << endl;

    return 0;
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
