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
    int id, n, t, tTime, nAtoms, yStep, currentTimeStep, typ;
    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r;
    double TAU1, TAU2, TAU3, pMiddle, rhoMiddle, keMiddle, tempMiddle, keKerogenTop, keKerogenBottom, tempKerogenTop, tempKerogenBottom;
    int nMethaneTop, nMethaneBottom, nMethaneMiddle, nMethaneKerogen;
    double kB = 1.38064852e-23;

    double Lx, Ly, Lz, volume;
    double xLo, xHi, yLo, yHi, zLo, zHi;

    double mi = 2.66389e-26; // mass of one molecule
    double Msigma = 3.73;    // diameter of one methane moleule
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

    int tSkip = 5000;
    double deltaT = 0.5;

    ifstream nSteps("nTimeSteps.dat", ios::in);
    int nTimeSteps; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l
    nSteps >> nTimeSteps;

    int skipTimeStep = 1;
    ifstream Height("Height.dat", ios::in);
    double H;
    Height >> H;

    ifstream data("dump_equil_gas.lammpstrj", ios::in);
    ofstream kerogenFile("equil_kerogen_vs_time.txt", ios::out);
    ofstream middleFile("equil_middle_vs_time.txt", ios::out);

    int count = 0;

    for (t = 0; t < nTimeSteps; t++)
    {
        for (n = 1; n < 10; n++)
        {
            if (n == 2)
            {
                data >> currentTimeStep;

                cout << "currentTimeStep = " << currentTimeStep
                     << "; t = " << t << " [ "
                     << 100 * float(t + 1) / float(nTimeSteps)
                     << "% ]" << endl;
            }

            if (n == 4)
            {
                data >> nAtoms;
                //                 cout << "nAtoms = " << nAtoms << endl;
            }

            if (n == 6)
            {
                data >> xLo >> xHi;
            }

            if (n == 7)
            {
                data >> yLo >> yHi;
            }

            if (n == 8)
            {
                data >> zLo >> zHi;
            }

            getline(data, line);
        }

        TAU1 = 0;
        TAU2 = 0;
        TAU3 = 0;

        nMethaneTop = 0;
        nMethaneBottom = 0;
        nMethaneMiddle = 0;
        keMiddle = 0.0;
        keKerogenTop = 0.0;
        keKerogenBottom = 0.0;

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3;

            if (typ == 1) // methane
            {
                if (z < 0)
                {
                    nMethaneBottom++;
                    keKerogenBottom += mi * (vx * vx + vy * vy + vz * vz);
                }
                else if (z > H)
                {
                    nMethaneTop++;
                    keKerogenTop += mi * (vx * vx + vy * vy + vz * vz);
                }
                //else if( (z >= 0.5*sigmaWall) && ( z <= (Hactual+0.5*sigmaWall)) )
                else
                {
                    nMethaneMiddle++;

                    TAU1 += tau1;
                    TAU2 += tau2;
                    TAU3 += tau3;
                    keMiddle += mi * (vx * vx + vy * vy + vz * vz);
                }
            }
        }

        Lx = xHi - xLo;
        Ly = yHi - yLo;
        volume = H * Lx * Ly;
        cout << "The volume of the middle: " << volume << endl;

        pMiddle = -((TAU1 + TAU2 + TAU3) / (3.0 * volume)) * refPressure / 1e6;
        rhoMiddle = mi * nMethaneMiddle / (volume * refLength * refLength * refLength);

        tempMiddle = 0.0;

        if (nMethaneMiddle > 0)
        {
            tempMiddle = keMiddle * refVelocity * refVelocity / (3.0 * kB * double(nMethaneMiddle));
        }

        tempKerogenBottom = 0.0;
        tempKerogenTop = 0.0;

        if (nMethaneBottom > 0)
        {
            tempKerogenBottom = keKerogenBottom * refVelocity * refVelocity / (3.0 * kB * double(nMethaneBottom));
        }

        if (nMethaneTop > 0)
        {
            tempKerogenTop = keKerogenTop * refVelocity * refVelocity / (3.0 * kB * double(nMethaneTop));
        }

        kerogenFile << t * deltaT * tSkip * refTime << '\t'
                    << nMethaneBottom << '\t'
                    << tempKerogenBottom << '\t'
                    << nMethaneTop << '\t'
                    << tempKerogenTop << '\t'
                    << endl;

        middleFile << t * deltaT * tSkip * refTime << '\t'
                   << pMiddle << '\t'
                   << rhoMiddle << '\t'
                   << nMethaneMiddle << '\t'
                   << tempMiddle << '\t'
                   << endl;

        getline(data, line);
    }

    return 0;
}
