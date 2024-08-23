#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

using namespace std;

int main()
{
    string line, fName;
    int id, n, i, t, nAtoms, yStep, currentTimeStep, typ;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE;
    double xLo, xHi, yLo, yHi, zLo, zHi, Lz;
    
    double nDynamicB, nDynamicT, keDynamicB, keDynamicT, tempDynamicB, tempDynamicT;
    double nBulkB,    nBulkT,    keBulkB,    keBulkT,    tempBulkB,    tempBulkT;
    double nAllB,     nAllT,     keAllB,     keAllT,     tempAllB,     tempAllT;

    double kB = 1.38064852e-23;
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs

    //***************************
    // Specify the parameters for the post-processing
    int nTimeSteps;
    double deltaT, tSkip, skipTimeStep;
    double rCut, H;
    double Tg, Tw;
    double mG, mW, mGas, mWall;

    int Nlines = 5;        // find out the No. of lines of parameters  "wc -l < Parameters.dat"
    ifstream Specification("Specification.dat", ios::in);
    if (!Specification.is_open())
    {
        cerr << "Error opening Specification.dat file" << endl;
        return 1;
    }
    for (int n = 1; n < Nlines + 1; n++)
    {
        if (n == 1)  {Specification >> nTimeSteps;                                 getline(Specification, line);}
        if (n == 2)  {Specification >> deltaT >> tSkip >> skipTimeStep;            getline(Specification, line);}
        if (n == 3)  {Specification >> rCut >> H;                                  getline(Specification, line);}
        if (n == 4)  {Specification >> Tg >> Tw;                                   getline(Specification, line);}
        if (n == 5)  {Specification >> mG >> mW;                                   getline(Specification, line);}
    }
    mGas  = mG * refMass;
    mWall = mW * refMass;
    //***************************

    ifstream data("dump_meas_wall.lammpstrj", ios::in);

    ofstream topFile("MeasTemp_TopWall.txt", ios::out);
    ofstream bottomFile("MeasTemp_BottomWall.txt", ios::out);

    nTimeSteps /= 100;

    for (t = 0; t < nTimeSteps; t++)
    {
        for (n = 1; n < 10; n++)
        {
            if (n == 4)
            {
                data >> nAtoms;
                // cout << "nAtoms = " << nAtoms << endl;
            }

            if (n == 2)
            {
                data >> currentTimeStep;

                cout << "currentTimeStep = " << currentTimeStep
                     << "; t = " << t << " [ "
                     << 100 * float(t + 1) / float(nTimeSteps)
                     << "% ]" << endl;
            }

            getline(data, line);
        }
        
        nDynamicB = 0, nDynamicT = 0, keDynamicB = 0.0, keDynamicT = 0.0, tempDynamicB = 0.0, tempDynamicT = 0.0;
        nBulkB    = 0, nBulkT    = 0, keBulkB    = 0.0, keBulkT    = 0.0, tempBulkB    = 0.0, tempBulkT    = 0.0;
        nAllB     = 0, nAllT     = 0, keAllB     = 0.0, keAllT     = 0.0, tempAllB     = 0.0, tempAllT     = 0.0;

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz;

            if (typ == 2)
            {
                nBulkB++;
                keBulkB += mWall * (vx * vx + vy * vy + vz * vz);

                nAllB++;
                keAllB += mWall * (vx * vx + vy * vy + vz * vz);
            }
            if (typ == 3)
            {
                nBulkT++;
                keBulkT += mWall * (vx * vx + vy * vy + vz * vz);

                nAllT++;
                keAllT += mWall * (vx * vx + vy * vy + vz * vz);
            }
            if (typ == 5)
            {
                nDynamicB++;
                keDynamicB += mWall * (vx * vx + vy * vy + vz * vz);

                nAllB++;
                keAllB += mWall * (vx * vx + vy * vy + vz * vz);
            }
            if (typ == 6)
            {
                nDynamicT++;
                keDynamicT += mWall * (vx * vx + vy * vy + vz * vz);

                nAllT++;
                keAllT += mWall * (vx * vx + vy * vy + vz * vz);
            }
        }

        if ((nDynamicB > 0) && (nBulkB > 0) && (nAllB > 0))
        {
            tempDynamicB = keDynamicB * refVelocity * refVelocity / (3.0 * kB * double(nDynamicB));
            tempBulkB    = keBulkB    * refVelocity * refVelocity / (3.0 * kB * double(nBulkB));
            tempAllB     = keAllB     * refVelocity * refVelocity / (3.0 * kB * double(nAllB));
        }

        if ((nDynamicT > 0) && (nBulkT > 0) && (nAllT > 0))
        {
            tempDynamicT = keDynamicT * refVelocity * refVelocity / (3.0 * kB * double(nDynamicT));
            tempBulkT    = keBulkT    * refVelocity * refVelocity / (3.0 * kB * double(nBulkT));
            tempAllT     = keAllT     * refVelocity * refVelocity / (3.0 * kB * double(nAllT));
        }

        getline(data, line);

        if (t > 0)
        {
            bottomFile << deltaT * t * tSkip * refTime << '\t'
                       << tempDynamicB << '\t'
                       << tempBulkB  << '\t'
                       << tempAllB  << '\t'
                       << endl;        

            topFile << deltaT * t * tSkip * refTime << '\t'
                    << tempDynamicT << '\t'
                    << tempBulkT << '\t'
                    << tempAllT << '\t'
                    << endl;
        }
    }

    return 0;
}