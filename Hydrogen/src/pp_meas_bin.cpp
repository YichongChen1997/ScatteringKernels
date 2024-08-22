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
    double TAU1, TAU2, TAU3, p, rho, vAVG, vAVG_kerogen, keKerogen, keMiddle, tempMiddle, tempKerogen;
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
    double mi = 6.6465E-27;  // mass of one Helium molecule
    double MCsigma = 2.640;  // sigma of Helium

    // Specify the parameters for the post-processing
    double H, rCut;
    int nTimeSteps;

    ifstream Specification("Specification.dat", ios::in);
    Specification >> H >> rCut >> nTimeSteps;
    // ----------------------------------------------------

    int skipTimeStep = 101; // when to start measuring bins -> after steady state
    int tSkip = 500;
    double deltaT = 2.0;

    double binWidth = 0.5;    // coarse
    double binWidth2 = 0.1; // fine

    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    //**************************
    int nBins = 0, nBins2 = 0;

    vector<double> molFieldTemp;
    vector<double> velFieldTemp;
    vector<double> keFieldTemp;
    vector<double> molField;
    vector<double> keField;
    vector<double> velField;
    vector<double> stressField;

    vector<double> molFieldTemp2;
    vector<double> velFieldTemp2;
    vector<double> keFieldTemp2;
    vector<double> molField2;
    vector<double> keField2;
    vector<double> velField2;
    vector<double> stressField2;

    int count = 0;
    double nAvSteps = 0.0;

    int bz, bz2;
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
            // I'm trying to make profiles symmetric
            Lx = xHi - xLo;
            Ly = yHi - yLo;
            Lz = zHi - zLo;
            nBins = ceil(Lz / binWidth);

            int halfnBins = ceil(nBins / 2);
            nBins = 2 * halfnBins;
            binWidth = Lz / double(nBins);

            cout << "binwidth update = " << binWidth << endl;
            cout << "nBins = " << nBins << endl;

            molFieldTemp.resize(nBins, 0.0);
            velFieldTemp.resize(nBins, 0.0);
            keFieldTemp.resize(nBins, 0.0);

            molField.resize(nBins, 0.0);
            keField.resize(nBins, 0.0);
            velField.resize(nBins, 0.0);
            stressField.resize(nBins, 0.0);

            nBins2 = ceil(Lz / binWidth2);
            int halfnBins2 = ceil(nBins2 / 2);
            nBins2 = 2 * halfnBins2;
            binWidth2 = Lz / double(nBins2);

            cout << "binwidth update = " << binWidth2 << endl;
            cout << "nBins = " << nBins2 << endl;

            molFieldTemp2.resize(nBins2, 0.0);
            velFieldTemp2.resize(nBins2, 0.0);
            keFieldTemp2.resize(nBins, 0.0);

            molField2.resize(nBins2, 0.0);
            keField2.resize(nBins2, 0.0);
            velField2.resize(nBins2, 0.0);
            stressField2.resize(nBins2, 0.0);
        }

        // store the lagrangian properties first because we need to calculate velocity first
        //vector<double> xR(nAtoms, 0.0);
        //vector<double> yR(nAtoms, 0.0);
        vector<double> zR(nAtoms, 0.0);
        vector<double> vxR(nAtoms, 0.0);
        vector<double> vyR(nAtoms, 0.0);
        vector<double> vzR(nAtoms, 0.0);
        vector<double> tau1R(nAtoms, 0.0);
        vector<double> tau2R(nAtoms, 0.0);
        vector<double> tau3R(nAtoms, 0.0);

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3;

            //xR[n]=x;
            //yR[n]=y;
            zR[n] = z;
            vxR[n] = vx;
            vyR[n] = vy;
            vzR[n] = vz;
            tau1R[n] = tau1;
            tau2R[n] = tau2;
            tau3R[n] = tau3;
        }

        // reset fields before use
        for (i = 0; i < nBins; i++)
        {
            molFieldTemp[i] = 0.0;
            velFieldTemp[i] = 0.0;
        }

        for (i = 0; i < nBins2; i++)
        {
            molFieldTemp2[i] = 0.0;
            velFieldTemp2[i] = 0.0;
        }

        // first calculate mean velocity in bin
        for (n = 0; n < nAtoms; n++)
        {
            z = zR[n];
            vx = vxR[n];

            if (t >= skipTimeStep)
            {
                if ((z >= zLo) && (z <= zHi))
                {
                    bz = floor((z - zLo) / binWidth);

                    if (bz < 0)
                    {
                        bz = 0;
                    }
                    if (bz >= nBins)
                    {
                        bz = nBins - 1;
                    }

                    molFieldTemp[bz] += 1.0;
                    velFieldTemp[bz] += vx;

                    // fine measurements

                    bz2 = floor((z - zLo) / binWidth2);

                    if (bz2 < 0)
                    {
                        bz2 = 0;
                    }
                    if (bz2 >= nBins2)
                    {
                        bz2 = nBins2 - 1;
                    }

                    molFieldTemp2[bz2] += 1.0;
                    velFieldTemp2[bz2] += vx;
                }
            }
        }

        for (i = 0; i < nBins; i++)
        {
            if (molFieldTemp[i] > 0)
            {
                velFieldTemp[i] /= molFieldTemp[i]; // mean velocity in bin
            }
        }

        for (i = 0; i < nBins2; i++)
        {
            if (molFieldTemp2[i] > 0)
            {
                velFieldTemp2[i] /= molFieldTemp2[i];
            }
        }

        // now calculate other properties
        for (n = 0; n < nAtoms; n++)
        {
            if (t >= skipTimeStep)
            {
                z = zR[n];
                vx = vxR[n];
                vy = vyR[n];
                vz = vzR[n];
                tau1 = tau1R[n];
                tau2 = tau2R[n];
                tau3 = tau3R[n];

                if ((z >= zLo) && (z <= zHi))
                {
                    bz = floor((z - zLo) / binWidth);

                    if (bz < 0)
                    {
                        bz = 0;
                    }
                    if (bz >= nBins)
                    {
                        bz = nBins - 1;
                    }

                    molField[bz] += 1.0;
                    keField[bz] += mi * (vx * vx + vy * vy + vz * vz - velFieldTemp[bz] * velFieldTemp[bz]);

                    velField[bz] += vx;
                    stressField[bz] += (tau1 + tau2 + tau3);

                    // fine measurements
                    bz2 = floor((z - zLo) / binWidth2);

                    if (bz2 < 0)
                    {
                        bz2 = 0;
                    }
                    if (bz2 >= nBins2)
                    {
                        bz2 = nBins2 - 1;
                    }

                    molField2[bz2] += 1.0;
                    keField2[bz2] += mi * (vx * vx + vy * vy + vz * vz - velFieldTemp2[bz2] * velFieldTemp2[bz2]);

                    velField2[bz2] += vx;
                    stressField2[bz2] += (tau1 + tau2 + tau3);
                }
            }
        }

        if (t >= skipTimeStep)
        {
            nAvSteps += 1.0;
        }

        //****
        getline(data, line);
    }

    // bin measurements

    cout << "bin measurements" << endl;

    {
        ofstream binFile1("bins_meas_coarse.txt", ios::out);

        double bin, rho = 0.0, vel = 0.0, temp = 0.0, binVol, binVol2, pressure;

        for (i = 0; i < nBins; i++)
        {
            bin = zLo + binWidth * 0.5 + binWidth * i; // coordinate of bin

            binVol = binWidth * Lx * Ly * refVolume; // volume of each bin

            rho = molField[i] * mi / (binVol * nAvSteps); // average density of each bin over time

            vel = 0.0;
            temp = 0.0;

            if (molField[i] > 0.0)
            {
                vel = velField[i] * refVelocity / molField[i];                            // average velocity of each bin over time
                temp = keField[i] * refVelocity * refVelocity / (3.0 * kB * molField[i]); // average temperature
            }

            pressure = -(stressField[i] / (3.0 * binVol * nAvSteps / refVolume)) * refPressure / 1e6;

            binFile1 << bin << '\t'
                    << rho << '\t'
                    << molField[i] / (binVol * nAvSteps) << '\t' // number density
                    << vel << '\t'
                    << temp << '\t'
                    << pressure << '\t'
                    << endl;
        }
    }

    {
        ofstream binFile2("bins_meas_fine.txt", ios::out);

        double bin, rho = 0.0, vel = 0.0, temp = 0.0, binVol, binVol2, pressure;

        for (i = 0; i < nBins2; i++)
        {
            bin = zLo + binWidth2 * 0.5 + binWidth2 * i;

            binVol2 = binWidth2 * Lx * Ly * refVolume;

            rho = molField2[i] * mi / (binVol2 * nAvSteps);

            vel = 0.0;
            temp = 0.0;

            if (molField2[i] > 0.0)
            {
                vel = velField2[i] * refVelocity / molField2[i];
                temp = keField2[i] * refVelocity * refVelocity / (3.0 * kB * molField2[i]);
                // temp = (keField2[i] - velFieldTemp2[i] * velFieldTemp2[i]) * refVelocity * refVelocity / (3.0 * kB * molField2[i]);
            }

            pressure = -(stressField2[i] / (3.0 * binVol2 * nAvSteps / refVolume)) * refPressure / 1e6;

            binFile2 << bin << '\t'
                    << rho << '\t'
                    << molField2[i] / (binVol2 * nAvSteps) << '\t' // number density
                    << vel << '\t'
                    << temp << '\t'
                    << pressure << '\t'
                    << endl;
        }
    }

    return 0;
}
