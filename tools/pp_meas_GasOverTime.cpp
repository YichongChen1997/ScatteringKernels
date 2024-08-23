#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

// Output the mean velocity and total mGas flow rate of gas molecules within the channel.
using namespace std;

int main()
{
    string line;
    int id, i, j, k, n, m, t, tTime, nAtoms, currentTimeStep, typ;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, fx, fy, fz;
    double xLo, xHi, yLo, yHi, zLo, zHi;
    double Lx, Ly, Lz, volume;

    double massFlowRate, vMeanGas, keGas, tempGas, pressureGas;
    double keGasTemp, stressTemp;

    double kB = 1.38064852e-23;
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;
    double refVolume = refLength * refLength * refLength;

    //***************************
    // Specify the parameters for the post-processing
    int nTimeSteps;
    double deltaT, tSkip, skipTimeStep;
    double rCut, H;
    double Tg, Tw;
    double mG, mW, mGas, mWall;

    int Nlines = 5;        // find out the No. of lines of parameters    "wc -l < Parameters.dat"
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

    ifstream data("dump_meas_gas.lammpstrj", ios::in);
    ofstream file("gasOverTime.txt", ios::out);

    //**************************
    int count = 0;
    int nAvSteps = 0;
    int dN = floor(nTimeSteps / 100);

    for (t = 0; t < nTimeSteps; t++)
    {
        for (n = 1; n < 10; n++)
        {
            if (n == 2)
            {
                data >> currentTimeStep;
                count++;
                if (count >= dN)
                {
                    cout << "currentTimeStep = " << currentTimeStep
                         << "; t = " << t << " [ "
                         << 100 * float(t + 1) / float(nTimeSteps)
                         << "% ]" << endl;

                    count = 0;
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

        vMeanGas = 0.0;
        keGasTemp = 0.0;

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3 >> KE >> PE >> fx >> fy >> fz;

            if (t >= skipTimeStep)
            {
                vMeanGas   += vx;
                keGasTemp  += (vx * vx + vy * vy + vz * vz);
                stressTemp += (tau1 + tau2 + tau3);
            }
        }

        if (t >= skipTimeStep)
        {
            nAvSteps++;
        } 

        // mGas Flow Rate (kg/s)
        massFlowRate = vMeanGas * mGas * refVelocity / (Lx * refLength);

        // Mean Velocity (m/s)
        vMeanGas /= double(nAtoms);

        // Kinetic Energy
        keGas = mGas * (keGasTemp - vMeanGas * vMeanGas * nAtoms);

        // Temperature of gas (K)
        tempGas = keGas * refVelocity * refVelocity / (3.0 * kB * nAtoms);

        // Pressure of gas (MPa)
        pressureGas = -(stressTemp / (3.0 * Lx * Ly * H * nAvSteps)) * refPressure / 1e6;

        file << deltaT * t * tSkip * refTime << '\t'
             << vMeanGas * refVelocity       << '\t'        
             << massFlowRate                 << '\t'
             << keGas                        << '\t' 
             << tempGas                      << '\t'
             << pressureGas << endl;

        getline(data, line);
    }

    return 0;
}