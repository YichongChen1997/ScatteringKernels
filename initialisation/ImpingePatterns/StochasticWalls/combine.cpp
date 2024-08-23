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

// DATA: https://ww2.chemistry.gatech.edu/~lw26/structure/small_molecules/index.html

int main()
{
    string line;
    double Mass, xLo, xHi, yLo, yHi, zLo, zHi, vx, vy, vz, Tg, Txshift, Tyshift, Nshift;
    double rCut, thickness, H, uC, mG, sigmaG, mW, zGas;
    int typ, i, j, k, mol = 0, molId = 0, co;
    int nAtomsTot = 0;

    double refLength = 1e-10;
    double refForce = 6.9477E-11;
    double refTime = 1e-15;
    double refMass = 1.66054e-27;
    double refVelocity = refLength / refTime;
    double kB = 1.38064852e-23;
    srand(time(NULL));

    int Nlines = 7;           // find out the No. of lines of parameters    "wc -l < Parameters.dat"
    ifstream Parameters("Parameters.dat", ios::in);
    for (int n = 1; n < Nlines + 1; n++)
    {
        if (n == 1)  {Parameters >> xLo >> xHi;                     getline(Parameters, line);}
        if (n == 2)  {Parameters >> yLo >> yHi;                     getline(Parameters, line);}
        if (n == 3)  {Parameters >> zLo >> zHi;                     getline(Parameters, line);}
        if (n == 4)  {Parameters >> zGas >> uC >> rCut;             getline(Parameters, line);}
        if (n == 5)  {Parameters >> Tg;                             getline(Parameters, line);}
        if (n == 6)  {Parameters >> Txshift >> Tyshift >> Nshift;   getline(Parameters, line);}
        if (n == 7)  {Parameters >> mG;                             getline(Parameters, line);}
    }
    
    //**** THINGS TO CHANGE *****
    double Lx = xHi - xLo;
    double Ly = yHi - yLo;
    double Lz = zHi - zLo;
    double mGas = mG * refMass;
    double vM = sqrt(2 * kB * Tg / mGas) / refVelocity;  // most probable velocity
    double vRef = sqrt(kB * Tg / mGas) / refVelocity;    // reference velocity used in velocity initialisation

    // no. of GasPr molecules in each direction
    int nPx = int(Lx / uC);
    int nPy = int(Ly / uC);
    int nPz = int(2 * zGas / uC);
    double uMx = Lx / nPx;
    double uMy = Ly / nPy;
    double uMz = 2 * zGas / nPz;

    // no. of GasPr molecules in the middle of the box.
    int nGasPr = nPx * nPy * nPz;

    vector<double> xG(nGasPr);
    vector<double> yG(nGasPr);
    vector<double> zG(nGasPr);

    cout << "No. of GasPr molecules to be inserted = " << nGasPr << endl;

    // Specify the coordinates of GasPr
    co = 0;
    for (i = 0; i < nPx; i++)
    {
        for (j = 0; j < nPy; j++)
        {
            for (k = 0; k < nPz; k++)
            {
                xG[co] = xLo + uMx * i + uMx * 0.5;
                yG[co] = yLo + uMy * j + uMy * 0.5;
                zG[co] = 1.5 * rCut + uMz * k;

                co++;
            }
        }
    }

    cout << "No. of GasPr molecules added = " << co << endl;
    nAtomsTot += co;

    cout << "nGasPr mols = " << co << endl
         << "Total atoms = " << nAtomsTot << endl
         << endl;

    ofstream writer("data.dat");

    writer << "LAMMPS data file" << endl;
    writer << endl;
    writer << nAtomsTot << '\t' << "atoms" << endl;
    writer << endl;
    writer << "1	atom types" << endl;
    writer << endl;
    writer << xLo << '\t' << xHi << "	xlo xhi" << endl;
    writer << yLo << '\t' << yHi << "	ylo yhi" << endl;
    writer << zLo << '\t' << zHi << "	zlo zhi" << endl;

    writer << endl;
    writer << "Masses" << endl;
    writer << endl;
    writer << "1  "<<  mG << endl;    // gas 
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
        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xG[i] << '\t' << yG[i] << '\t' << zG[i] << endl;

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