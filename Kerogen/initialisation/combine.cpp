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

    double Msigma = 3.73;   // LJ \sigma, diameter of the methane molecule
    double MCsigma = 3.545; // sigma of methane - C interaction

    srand(time(NULL));

    //**** THINGS TO CHANGE ******
    double rCut = 15;       // gas-wall interaction distance
    double zRigidRange = 2; // rigid bound of the wall
    double uCM = 1.5; // unit-cell spacing between methane molecules

    double atomForceSi = 8.0e-15;                       // Newtons (ensure the magnitude of velocity at decades of m/s)
    double atomForceReal = atomForceSi / (69.4786e-12); // Kcal/mole - Angstrom

    double pressure = 2.0; // MPa
    bool setPressure = true;

    double rho = 0; // kg/m^3 either define rho directly here or use the pressure
    //double rho = 10;

    // roughness
    double amp = 2.0;
    double period = 25;
    double zCut = 50;    // NOMINAL PLANE - choose a number between 0 and 50 to take as the surface to the methane channel
    double thickness = 15;
    
    // Simulation box dimension
    double X = 50; // rough estimates
    double Y = 50;
    double H = 100;
    
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
    zLo = -thickness - rCut - zRigidRange;
    xHi = X;
    yHi = Y;
    zHi = H + thickness + rCut + zRigidRange;

    double Lx = xHi - xLo;
    double Ly = yHi - yLo;
    double Lz = zHi - zLo;

    //-- create Kerogen wall
    ifstream readerKerogen("EFK_50A_0.80.xyz"); // kerogen
    readerKerogen >> nWallK;

    getline(readerKerogen, line);
    getline(readerKerogen, line);

    vector<double> xK0;
    vector<double> yK0;
    vector<double> zK0;
    vector<int> typK0;

    // read in kerogen block and find the preferred cutting plane

    for (i = 0; i < nWallK; i++)
    {
        readerKerogen >> typeName >> x >> y >> z;

        if (typeName == "C")
        {
            typ = 3;
        }
        else if (typeName == "O")
        {
            typ = 4;
        }
        else if (typeName == "H")
        {
            typ = 5;
        }

        x = x + 25;
        y = y + 25;
        z = z + 25;

        if (z > zCut) // take part above zCut and shift it below and start channel interface at 0,0,0
        {
            z = z - 50 - zCut;
        }
        else
        {
            z = z - zCut;
        }

        xK0.push_back(x); // slab of kerogen with channel interface at z = 0
        yK0.push_back(y);
        zK0.push_back(z);
        typK0.push_back(typ);

        xK0.push_back(x); // add a slab of kerogen on top for carving out roughness
        yK0.push_back(y);
        zK0.push_back(z+50);
        typK0.push_back(typ);
    }
    
    int nPts0 = xK0.size();
    
    cout << "No. of atoms in block 0 = " << xK0.size() << endl;    
    
    vector<double> xK1;
    vector<double> yK1;
    vector<double> zK1;
    vector<int> typK1;

    for (i = 0; i < nPts0; i++)
    {    
        x = xK0[i];
        y = yK0[i];
        z = zK0[i];
        typ = typK0[i];

        if (z >= -thickness) // ensure thickness between bottom of kerogen and nominal plane
        {
            xK1.push_back(x);
            yK1.push_back(y);
            zK1.push_back(z);
            typK1.push_back(typ);
        }
    }
    
    cout << "No. of atoms in block 1 = " << xK1.size() << endl;

    // now let's replicate the kerogen in x and y directions

    vector<double> xK2;
    vector<double> yK2;
    vector<double> zK2;
    vector<int> typK2;

    int nX = ceil(Lx / 50.0); // No. of kerogen blocks in that direction
    int nY = ceil(Ly / 50.0);

    int nPts1 = xK1.size();

    for (i = 0; i < nPts1; i++)
    {
        x = xK1[i];
        y = yK1[i];
        z = zK1[i];
        typ = typK1[i];

        for (k = 0; k < nX; k++)
        {
            xJ = x + k * 50;

            for (j = 0; j < nY; j++)
            {
                yJ = y + j * 50;

                xK2.push_back(xJ);
                yK2.push_back(yJ);
                zK2.push_back(z);
                typK2.push_back(typ);
            }
        }
    }

    cout << "No. of atoms in block 2 = " << xK2.size() << endl;

    // now let's apply roughness to surface

    vector<double> xK3;
    vector<double> yK3;
    vector<double> zK3;
    vector<int> typK3;

    int nPts2 = xK2.size();

    double roughSurfaceX = 0.0;
    double roughSurfaceY = 0.0;

    for (i = 0; i < nPts2; i++)
    {
        x = xK2[i];
        y = yK2[i];
        z = zK2[i];
        typ = typK2[i];

        roughSurfaceX = amp * (sin(2 * PI * x / period) );
        roughSurfaceY = amp * (sin(2 * PI * y / period) );	
        
        if ((z <= roughSurfaceX) && (z <= roughSurfaceY))
        {
            xK3.push_back(x);
            yK3.push_back(y);
            zK3.push_back(z);
            typK3.push_back(typ);
        }
    }

    cout << "No. of atoms in block 3 = " << xK3.size() << endl;

    // now transfer what we have done to the main vectors

    vector<double> xK;
    vector<double> yK;
    vector<double> zK;
    vector<int> typK;

    int nPts3 = xK3.size();

    for (i = 0; i < nPts3; i++)
    {
        x = xK3[i];
        y = yK3[i];
        z = zK3[i];
        typ = typK3[i];

        xK.push_back(x);
        yK.push_back(y);
        zK.push_back(z);
        typK.push_back(typ);
    }

    nWallK = xK.size(); // No. of kerogen particles on the bottom

    nAtomsTot = nWallK * 2;

    cout << "No. of atoms in the kerogen: " << nAtomsTot << endl;

    // Create bottom & top barrier

    double uC = 2.0;
    int nSx = int(Lx / uC);
    int nSy = int(Ly / uC);

    int nS = nSx * nSy;
    vector<double> xS(nS);
    vector<double> yS(nS);
    vector<double> zS(nS);

    cout << "No. of atoms on one piston = " << nS << endl;
    co = 0;
    for (i = 0; i < nSx; i++)
    {
        for (j = 0; j < nSy; j++)
        {
            xS[co] = uC * i + uC * 0.5;
            yS[co] = uC * j + uC * 0.5;
            zS[co] = 0.0;
            co++;
        }
    }

    cout << "No. of atoms on the bottom barrier = " << xS.size() << endl;

    nAtomsTot += nS * 2;

    // create methane

    double volume = Lx * Ly * H * refLength * refLength * refLength;

    int nMethane = ceil(rho * volume / mi);

    cout << "No. of methane atoms to insert = " << nMethane << endl;

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

    cout << " nB = " << nB << endl;

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

    cout << "No. of methane molecules added = " << co << endl;
    nAtomsTot += co;

    cout << "Total atoms = " << nAtomsTot << endl
         << "------------------------------------" << endl
         << endl;

    ofstream writer("data.dat");

    writer << "LAMMPS data file" << endl;
    writer << endl;
    writer << nAtomsTot << '\t' << "atoms" << endl;
    writer << endl;
    writer << "5	atom types" << endl;
    writer << endl;
    writer << xLo << '\t' << xHi << "	xlo xhi" << endl;
    writer << yLo << '\t' << yHi << "	ylo yhi" << endl;
    writer << zLo << '\t' << zHi << "	zlo zhi" << endl;

    writer << endl;
    writer << "Masses" << endl;
    writer << endl;
    writer << "1   16.043000" << endl; //CH4 - methane
    writer << "2   195.084" << endl;   //Barrier - Platinum
    writer << "3   12.0107" << endl;   //CARBON - KEROGEN
    writer << "4   15.9990" << endl;   //OXYGEN - KEROGEN
    writer << "5   1.00800" << endl;   //HYDROGEN - KEROGEN
    writer << endl;
    writer << "Atoms" << endl;
    writer << endl;

    typ = 1;
    for (i = 0; i < nMethane; i++)
    {
        molId++;
        mol++;
        // YC: shifted all molecules by "amp" to make the centreline of sinusoidal curve as z = 0
        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xM[i] << '\t' << yM[i] << '\t' << zM[i] << endl; 
    }

    int nConstWall = 0;

    for (i = 0; i < nS; i++)
    {
        typ = 2; // bottom barrier atoms
        molId++;
        mol++;

        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xS[i] << '\t' << yS[i] << '\t' << zS[i] - thickness << endl;
    }

    for (i = 0; i < nS; i++)
    {
        typ = 2; // top barrier atoms
        molId++;
        mol++;

        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xS[i] << '\t' << yS[i] << '\t' << zS[i] + H + thickness  << endl;
    }

    for (i = 0; i < nWallK; i++)
    {
        mol++; // bottom kerogen particles
        nKWallCount++;
        writer << mol << '\t' << mol << '\t' << typK[i] << '\t' << 0.0 << '\t' << xK[i] << '\t' << yK[i] << '\t' << zK[i] << endl;
    }

    for (i = 0; i < nWallK; i++)
    {
        mol++; // top kerogen particles
        nKWallCount++;
        writer << mol << '\t' << mol << '\t' << typK[i] << '\t' << 0.0 << '\t' << xK[i] << '\t' << yK[i] << '\t' << -zK[i] + H << endl;
    }

    // cout << "number of kerogen atoms placed = " << nKWallCount << endl;

    writer << endl;
    writer.close();

    double refNdenstiy = co / (Lx * Ly * H * pow(refLength, 3));
    double redNdensity = refNdenstiy * (PI * pow(Msigma, 3) * pow(refLength, 3) / 6);
    double confRatio = H / Msigma;
    double lambda = 1 / (sqrt(2) * PI * pow(Msigma, 2) * pow(refLength, 2) * refNdenstiy);
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
