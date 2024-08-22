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
    string line, typeName;
    char l;
    double x, y, z, chg, xLo, xHi, yLo, yHi, zLo, zHi, xT, xJ, yJ, zJ, vx, vy, vz, magV, rMag;
    int atoms, typ, i, j, k, c, m, id, mol = 0, molId = 0, bID, aID, nCNT = 0, nWallCount = 0, nWallK = 0, nKWallCount = 0, co;
    int nAtomsTot = 0;
    double shiftX = 0, shiftY = 0, shiftZ = 0;

    double mi = 2.66389E-26;
    double refLength = 1e-10;
    double refForce = 6.9477E-11;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime;
    double kB = 1.38064852e-23;

    double Msigma = 3.73;   // gas-gas interaction
    double Mdiameter = 3.8; // diameter of methane

    srand(time(NULL));

    //**** THINGS TO CHANGE *****

    double rCut = 15;       // gas-wall interaction distance
    double thickness = 15;  // thickness of the wall
    double zRigidRange = 5; // rigid bound of the wall

    // roughness
    double amp = 5.0;
    double period = 25;

    double X = 200;
    double Y = 200;
    double H = 40;
    double uCM = 1.5; // unit-cell spacing between methane molecules

    double atomForceSi = 8.0e-15;                       // Newtons (ensure the magnitude of velocity at decades of m/s)
    double atomForceReal = atomForceSi / (69.4786e-12); // Kcal/mole - Angstrom

    double pressure = 3.0; // MPa
    bool setPressure = true;

    double rho = 0;
    //double rho = 10; // kg/m^3 either define rho directly here or use the pressure

    //***************************
    // density to be set automatically
    if (setPressure)
    {

        double p = pressure;

        if (p < 50)
        {
            rho = 0.000009592344 * p * p * p * p + (-0.001207211400) * p * p * p + 0.022536131000 * p * p + 4.531223300000 * p + 0.021405916000;
        }
        else
        {
            rho = -0.000000399685 * p * p * p * p + 0.000213941490 * p * p * p + (-0.046716196000) * p * p + 5.698349000000 * p + (-0.360298820000);
        }
    }

    xLo = 0;
    yLo = 0;
    zLo = -thickness - rCut;
    xHi = X;
    yHi = Y;
    zHi = H + thickness + rCut;

    double Lx = xHi - xLo;
    double Ly = yHi - yLo;
    double Lz = zHi - zLo;

    // create top and bottom wall

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

    double roughSurfaceX = 0.0;
    double roughSurfaceY = 0.0;

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
                    zP = thickness - (rIz + unitCellSize * k); // just be careful of this

                    roughSurfaceX = amp * (sin(2 * PI * xP / period));
                    roughSurfaceY = amp * (sin(2 * PI * yP / period));

                    if ((zP <= roughSurfaceX) && (zP <= roughSurfaceY))
                    {
                        xWall.push_back(xP);
                        yWall.push_back(yP);
                        zWall.push_back(zP);
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

    // create methane
    double volume = (xHi - xLo) * (yHi - yLo) * H * refLength * refLength * refLength;

    int nMethane = ceil(rho * volume / mi);

    cout << "number of methane atoms to insert = " << nMethane << endl;

    double offsetX = 3;
    double offsetY = 3;
    double offsetZ = 15;
    int nPx, nPy, nPz, nP;

    Lz = H; // over-ride to ensure code below works

    vector<double> xM(nMethane);
    vector<double> yM(nMethane);
    vector<double> zM(nMethane);

    // RiGHT
    double rhoR = nMethane / (Lx * Ly * Lz);
    double Nbox = rhoR * Lx * Lx * Lx;
    int nB = int(pow(Nbox, 1.0 / 3));

    cout << "nB = " << nB << endl;

    nPx = int(((Lx - offsetX) * nB / Lx) + 1);
    nPy = int(((Ly - offsetY) * nB / Lx) + 1);
    nPz = int(((Lz - offsetZ) * nB / Lx) + 1);
    nP = 0;
    bool search = true;
    while (search)
    {
        nP = nPx * nPy * nPz;
        if (nP > nMethane)
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
    double uCZ = (Lz - offsetZ) / double(nPz);

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
                if (co < nMethane)
                {
                    xM[co] = xLo + i * uCX + offsetX;
                    yM[co] = yLo + j * uCY + offsetY;
                    zM[co] = 0 + k * uCZ + offsetZ;
                    co++;
                }
            }
        }
    }

    cout << " no. of methane molecules added = " << co << endl;
    nAtomsTot += co;

    cout << "nMethane mols = " << co << endl
         << "nWall atoms = " << nWall << endl
         << "Total atoms = " << nAtomsTot << endl
         << "------------------------------------" << endl
         << endl;

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
    writer << "1   16.043000" << endl; //CH4 - methane
    writer << "2   195.084" << endl;   //wall bottom
    writer << "3   195.084" << endl;   //wall top
    writer << "4   195.084" << endl;   //constrained wall
    writer << endl;
    writer << "Atoms" << endl;
    writer << endl;

    typ = 1;
    for (i = 0; i < nMethane; i++)
    {
        molId++;
        mol++;
        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xM[i] << '\t' << yM[i] << '\t' << zM[i] << endl;
    }

    cout << "Number of methane atoms placed = " << mol << endl;

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

    cout << "nConstWall = " << nConstWall << endl;
    cout << endl;

    writer << endl;

    writer.close();

    double refNdenstiy = co / (Lx * Ly * H * pow(refLength, 3));
    double redNdensity = refNdenstiy * (PI * pow(Msigma, 3) * pow(refLength, 3) / 6);
    double confRatio = H / Msigma;
    double lambda = 1 / (sqrt(2) * PI * pow(Mdiameter, 2) * pow(refLength, 2) * refNdenstiy);
    double Kn = lambda / (H * refLength);

    cout << "Reference number density = " << refNdenstiy << endl;
    cout << "Reduced number density = " << redNdensity << endl;
    cout << "Confinemen ratio R = " << confRatio << endl;
    cout << "Mean free path = " << lambda << endl;
    cout << "Knudsen number = " << Kn << endl;
    cout << "External force per atom in SI unit: " << atomForceSi << endl;
    cout << "External force per atom in real unit: " << atomForceReal << endl;
    cout << endl;

    return 0;
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
