#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

using namespace std;

int main()
{
    string line, fName, typeName;
    int id, i, j, k, n, m, t, tTime, nAtoms, yStep, currentTimeStep, typ;
    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, myPE, myKE, r, rMag, xJ, yJ, zJ, vxT;
    double TAU1, TAU2, TAU3, p, rho, keMiddle, tempMiddle, tempKerogen;
    double velMiddle, velKerogen, volume2, volume3, rhoN1, rhoN2, rhoN3, rhoKerogen;
    double velMiddle2, velMiddle3, massFlowRate, massFlowRate2, massFlowRate3;
    int nMethaneKerogen, nMethaneMiddle, N1, N2, N3;
    double TAU1b, TAU2b, TAU3b, p2, TAU1c, TAU2c, TAU3c, p3;
    double xLo, xHi, yLo, yHi, zLo, zHi;
    double Lx, Ly, Lz, volume;
    double kB = 1.38064852e-23;

    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;
    double refVolume = refLength * refLength * refLength;

    //********* INPUT ***********
    ifstream nSteps("nTimeSteps.dat", ios::in);
    int nTimeSteps; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l
    nSteps >> nTimeSteps;

    int skipTimeStep = 100; // when to start measuring bins -> after steady state

    int tSkip = 500; // match the output interval in lammps

    double deltaT = 0.5;
    double mi = 2.66389e-26;

    ifstream data("dump_meas_gas.lammpstrj", ios::in);
    ofstream file("massFlowRate.dat", ios::out);

    //**************************
    double vAVG = 0.0, mDot = 0.0;
    int count = 0;
    double nAvSteps = 0.0;

    int bz;
    int c = 0;
    int dN = floor(nTimeSteps / 100);

    for (t = 0; t < nTimeSteps; t++)
    {
        for (n = 1; n < 10; n++)
        {
            if (n == 2)
            {
                data >> currentTimeStep;
                c++;
                if (c >= dN)
                {
                    cout << "currentTimeStep = " << currentTimeStep
                         << "; t = " << t << " [ "
                         << 100 * float(t + 1) / float(nTimeSteps)
                         << "% ]" << endl;

                    c = 0;
                }
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

        if (t == 0)
        {
            Lx = xHi - xLo;
            Ly = yHi - yLo;
            Lz = zHi - zLo;
        }

        vAVG = 0.0;

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3;

            vAVG += vx;
        }

        mDot = vAVG * mi * refVelocity / (Lx * refLength); // mass flow rate (kg/s)

        vAVG /= double(nAtoms); // average velocity for all molecules, i.e. bulk velocity (m/s)

        file << deltaT * t * tSkip * refTime << '\t' << mDot << '\t' << vAVG * refVelocity << endl;

        getline(data, line);
    }

    return 0;
}
