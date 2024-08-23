#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <time.h>
#include <cstdlib>

using namespace std;
const float PI = 3.1415926;

double fRand(double fMin, double fMax);

// DATA: https://ww2.chemistry.gatech.edu/~lw26/structure/small_molecules/index.html

int main()
{
    string line;
    double x, y, z, vx, vy, vz, xLo, xHi, yLo, yHi, zLo, zHi, Tg, Txshift, Tyshift, Nshift;
    double rCut, thickness, H, uC, mG, sigmaG, mW, zGas;
    int typ, i, j, k, c, m, mol = 0, molId = 0, co;
    int nAtomsTot = 0;

    // reference scales (real to SI unit)
    double refLength = 1e-10;
    double refForce = 6.9477E-11;
    double refTime = 1e-15;
    double refMass = 1.66054e-27;
    double refVelocity = refLength / refTime;
    double kB = 1.38064852e-23;
    double zRigidRange = 3; // rigid bound of the wall
    double zDynamicRange;   // top layer of the wall

    int Nlines = 8;         // find out the No. of lines of parameters    "wc -l < Parameters.dat"
    ifstream Parameters("Parameters.dat", ios::in);
    for (int n = 1; n < Nlines + 1; n++)
    {
        if (n == 1)  {Parameters >> xLo >> xHi;                                 getline(Parameters, line);}
        if (n == 2)  {Parameters >> yLo >> yHi;                                 getline(Parameters, line);}
        if (n == 3)  {Parameters >> H >> rCut >> thickness >> zDynamicRange;    getline(Parameters, line);}
        if (n == 4)  {Parameters >> zGas >> uC;                                 getline(Parameters, line);}
        if (n == 5)  {Parameters >> Tg;                                         getline(Parameters, line);}
        if (n == 6)  {Parameters >> Txshift >> Tyshift >> Nshift;               getline(Parameters, line);}
        if (n == 7)  {Parameters >> mW;                                         getline(Parameters, line);}
        if (n == 8)  {Parameters >> mG;                                         getline(Parameters, line);}
    }

    srand(time(NULL));

    double Lx = xHi - xLo;
    double Ly = yHi - yLo;
    zLo = - (thickness + rCut);
    zHi = H + (thickness + rCut);
    double mGas = mG * refMass;
    double vM = sqrt(2 * kB * Tg / mGas) / refVelocity;  // most probable velocity
    double vRef = sqrt(kB * Tg / mGas) / refVelocity;    // reference velocity used in velocity initialisation

    // create top and bottom wall
    cout << "The most probable speed of gas is: " << vM * refVelocity << " m/s" << endl;
    cout << "The reference speed of gas is: " << vRef * refVelocity << " m/s" << endl;
    cout << "------------------------------------" << endl;
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
    int nZ = ceil(thickness / unitCellSize);

    // define new Lx and H, and box geometry
    Lx = nX * unitCellSize;
    xLo = 0.0;
    xHi = Lx;

    Ly = nY * unitCellSize;
    yLo = 0.0;
    yHi = Ly;

    cout << "------------------------------------" << endl;
    cout << "change of box co-ordinates" << endl;
    cout << "xLo = " << xLo << endl;
    cout << "xHi = " << xHi << endl;
    cout << "yLo = " << yLo << endl;
    cout << "yHi = " << yHi << endl;
    cout << "------------------------------------" << endl;

    vector<double> xWall;
    vector<double> yWall;
    vector<double> zWall;

    double rIx, rIy, rIz, xP, yP, zP;
    double smallOffset = 0.2;
    double roughSurface = 0.0;

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
                    zP = - (rIz + unitCellSize * k); // just be careful of this

                    xWall.push_back(xP);
                    yWall.push_back(yP);
                    zWall.push_back(zP);
                }
            }
        }
    }

    int nWall = xWall.size();

    vector<double> xWall2;
    vector<double> yWall2;
    vector<double> zWall2;

    double zMax = 0.0;
    double zMin = 0.0;

    for (i = 0; i < nWall; i++)
    {
        z = -zWall[i] + H;
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

    // no. of gas molecules in each direction
    int nPx = int(Lx / uC);
    int nPy = int(Ly / uC);
    int nPz = int(2 * zGas / uC);
    double uMx = Lx / nPx;
    double uMy = Ly / nPy;
    double uMz = 2 * zGas / nPz;

    // no. of gas molecules in the middle of the box.
    int nGasPr = nPx * nPy * nPz;

    vector<double> xM(nGasPr);
    vector<double> yM(nGasPr);
    vector<double> zM(nGasPr);

    cout << "No. of GasPr molecules = " << nGasPr << endl;

    // Specify the coordinates of GasPr
    co = 0;
    for (i = 0; i < nPx; i++)
    {
        for (j = 0; j < nPy; j++)
        {
            for (k = 0; k < nPz; k++)
            {
                xM[co] = xLo + uMx * i + uMx * 0.5;
                yM[co] = yLo + uMy * j + uMy * 0.5;
                zM[co] = 2 * rCut + uMz * k;

                co++;
            }
        }
    }

    cout << "No. of GasPr molecules added = " << co << endl;
    nAtomsTot += co;

    cout << "nGasPr mols = " << co << endl
         << "nWall atoms = " << 2 * nWall << endl
         << "Total atoms = " << nAtomsTot << endl
         << endl;

    ofstream writer("data.dat");

    writer << "LAMMPS data file" << endl;
    writer << endl;
    writer << nAtomsTot << '\t' << "atoms" << endl;
    writer << endl;
    writer << "6	atom types" << endl;
    writer << endl;
    writer << xLo << '\t' << xHi << "	xlo xhi" << endl;
    writer << yLo << '\t' << yHi << "	ylo yhi" << endl;
    writer << zLo << '\t' << zHi << "	zlo zhi" << endl;

    writer << endl;
    writer << "Masses" << endl;
    writer << endl;
    writer << "1 " << mG << endl;   // gas
    writer << "2 " << mW << endl;   // wall bottom
    writer << "3 " << mW << endl;   // wall top
    writer << "4 " << mW << endl;   // constrained wall
    writer << "5 " << mW << endl;   // bottom dynamic wall
    writer << "6 " << mW << endl;   // top    dynamic wall
    writer << endl;
    writer << "Atoms" << endl;
    writer << endl;

    vector<double> vxAtoms;
    vector<double> vyAtoms;
    vector<double> vzAtoms;

    double rnd, theta, rho;
    typ = 1;
    for (i = 0; i < nGasPr; i++)
    {
        molId++;
        mol++;
        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xM[i] << '\t' << yM[i] << '\t' << zM[i] << endl;

        // Maxwellian Impingement
        rnd = 0 + ((double)rand() / RAND_MAX) * (1 - 0);
        theta = 2 * PI * rnd;
        rnd = 0 + ((double)rand() / RAND_MAX) * (1 - 0);
        rho = sqrt(-2 * log(rnd)); 
        vx = rho * cos(theta);     
        vy = rho * sin(theta);
        rnd = 0 + ((double)rand() / RAND_MAX) * (1 - 0);
        vz = sqrt(-2 * log(rnd));

        vx = vx * vRef + Txshift * vM;
        vy = vy * vRef + Tyshift * vM;
        vz = vz * vRef + Nshift  * vM;

        vxAtoms.push_back(vx);
        vyAtoms.push_back(vy);
        vzAtoms.push_back(-vz);
    }

    cout << "Number of GasPr atoms placed = " << mol << endl;

    int nConstWall = 0;
    int nBotDynamic = 0;
    int nTopDynamic = 0;

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
        if ((zWall[i] >=  - zDynamicRange) && (zWall[i] <= 0))
        {
            typ = 5;
            nBotDynamic++;
        }

        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xWall[i] << '\t' << yWall[i] << '\t' << zWall[i] << endl;

        // randomly pick a unit vector and make it equal to vMag
        vx = 0.0;
        vy = 0.0;
        vz = 0.0;

        // cout << vx << " " << vy << " " << vz << endl;
        vxAtoms.push_back(vx);
        vyAtoms.push_back(vy);
        vzAtoms.push_back(vz);
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
        if ((zWall2[i] >=  H) && (zWall2[i] <= H + zDynamicRange))
        {
            typ = 6;
            nTopDynamic++;
        }

        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xWall2[i] << '\t' << yWall2[i] << '\t' << zWall2[i] << endl;

        // randomly pick a unit vector and make it equal to vMag
        vx = 0.0;
        vy = 0.0;
        vz = 0.0;

        // cout << vx << " " << vy << " " << vz << endl;
        vxAtoms.push_back(vx);
        vyAtoms.push_back(vy);
        vzAtoms.push_back(vz);
    }

    cout << " nConstWall = " << nConstWall << endl;

    //**** add velocities
    writer << endl;
    writer << "Velocities" << endl;
    writer << endl;

    for (i = 0; i < nAtomsTot; i++)
    {
        writer << i + 1 << '\t' << vxAtoms[i] << '\t' << vyAtoms[i] << '\t' << vzAtoms[i] << endl;
    }

    writer << endl;

    writer.close();

    return 0;
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}