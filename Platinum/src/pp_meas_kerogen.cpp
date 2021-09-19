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
    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r, rMag, xJ, yJ, zJ;
    double TAU1, TAU2, TAU3, p, rhoMiddle, vAVG_Middle, vAVG_kerogenTop, vAVG_kerogenBottom, keMiddle, keKerogenTop, keKerogenBottom, tempMiddle, tempTop, tempBottom;
    double velTop, velBottom, volume2, rhoN1, rhoN2, rhoKerogen;
    double molTempTop, molTempBottom, molTempMiddle;
    double velTempTop, velTempBottom, velTempMiddle;
    double vAVG_Middle2, massFlowRate, massFlowRate2;
    int nMethaneTop, nMethaneBottom, nMethaneMiddle, N1, nMethaneMiddle2;
    double TAU1b, TAU2b, TAU3b, p2;
    double xLo, xHi, yLo, yHi, zLo, zHi;
    double Lx, Ly, Lz, volume;
    double kB = 1.38064852e-23;

    double mi = 2.66389e-26; // mass of one molecule
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

    // ***************************
    int tSkip = 500;
    double deltaT = 0.5;

    ifstream nSteps("nTimeSteps.dat", ios::in);
    int nTimeSteps; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l
    nSteps >> nTimeSteps;

    int skipTimeStep = 1;
    ifstream Height("Height.dat", ios::in);
    double H;
    Height >> H;

    double rCut = 15;

    double Msigma = 3.73; // diameter of one methane moleule
    // ***************************

    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    ofstream kerogenFile("meas_kerogen_vs_time.txt", ios::out);
    ofstream middleFile("meas_middle_vs_tims.txt", ios::out);

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

        molTempTop = 0.0;
        molTempBottom = 0.0;
        molTempMiddle = 0.0;

        velTempTop = 0.0;
        velTempBottom = 0.0;
        velTempMiddle = 0.0;

        // first calculate the mean velocity along the channel
        for (n = 0; n < nAtoms; n++)
        {
            z = zR[n];
            vx = vxR[n];

            if (t > skipTimeStep)
            {
                if (z < 0)
                {
                    molTempBottom += 1.0;
                    velTempBottom += vx;
                }
                if (z > H)
                {
                    molTempTop += 1.0;
                    velTempTop += vx;
                }
                else
                {
                    molTempMiddle += 1.0;
                    velTempMiddle += vx;
                }
            }
        }

        TAU1 = 0;
        TAU2 = 0;
        TAU3 = 0;

        TAU1b = 0;
        TAU2b = 0;
        TAU3b = 0;

        nMethaneBottom = 0;
        nMethaneTop = 0;
        nMethaneMiddle = 0;

        vAVG_kerogenBottom = 0;
        vAVG_kerogenTop = 0;
        vAVG_Middle = 0;
        vAVG_Middle2 = 0;

        keMiddle = 0.0;
        keKerogenTop = 0.0;
        keKerogenBottom = 0.0;

        nMethaneMiddle2 = 0;

        // now calculate other properties
        for (n = 0; n < nAtoms; n++)
        {
            z = zR[n];
            vx = vxR[n];
            vy = vyR[n];
            vz = vzR[n];
            tau1 = tau1R[n];
            tau2 = tau2R[n];
            tau3 = tau3R[n];

            if (z < 0)
            {
                nMethaneBottom++;
                vAVG_kerogenBottom += vx;
                keKerogenBottom += mi * (vx * vx + vy * vy + vz * vz - (velTempBottom / molTempBottom) * (velTempBottom / molTempBottom));
            }
            else if (z > H)
            {
                nMethaneTop++;
                vAVG_kerogenTop += vx;
                keKerogenTop += mi * (vx * vx + vy * vy + vz * vz - (velTempTop / molTempTop) * (velTempTop / molTempTop));
            }
            // H actual
            else
            {
                nMethaneMiddle++;

                TAU1 += tau1;
                TAU2 += tau2;
                TAU3 += tau3;

                vAVG_Middle += vx;
                keMiddle += mi * (vx * vx + vy * vy + vz * vz - (velTempMiddle / molTempMiddle) * (velTempMiddle / molTempMiddle));
            }

            // H - Msigma
            if ((z >= (0.5 * Msigma)) && (z <= (H - 0.5 * Msigma)))
            {
                nMethaneMiddle2++;
                vAVG_Middle2 += vx;
                TAU1b += tau1;
                TAU2b += tau2;
                TAU3b += tau3;
            }
        }

        Lx = xHi - xLo;
        Ly = yHi - yLo;
        volume = H * Lx * Ly;
        volume2 = (H - Msigma) * Lx * Ly;
        cout << "The volume of the middle: " << volume << endl;

        // Pressure
        p = -((TAU1 + TAU2 + TAU3) / (3.0 * volume)) * refPressure / 1e6;
        p2 = -((TAU1b + TAU2b + TAU3b) / (3.0 * volume2)) * refPressure / 1e6;

        // Density
        rhoMiddle = mi * nMethaneMiddle / (volume * refLength * refLength * refLength);

        // Number density
        rhoN1 = nMethaneMiddle / (volume * refLength * refLength * refLength);
        rhoN2 = nMethaneMiddle2 / (volume2 * refLength * refLength * refLength);
        
        tempMiddle = 0.0;

        // Mean velocity
        if (nMethaneMiddle > 0)
        {
            tempMiddle = keMiddle * refVelocity * refVelocity / (3.0 * kB * double(nMethaneMiddle));
            vAVG_Middle *= refVelocity / double(nMethaneMiddle);
        }
        if (nMethaneMiddle2 > 0)
        {
            vAVG_Middle2 *= refVelocity / double(nMethaneMiddle2);
        }

        // Mass flow rate
        massFlowRate = vAVG_Middle * double(nMethaneMiddle) * mi / (Lx * refLength);
        massFlowRate2 = vAVG_Middle2 * double(nMethaneMiddle2) * mi / (Lx * refLength);

        // Properties of kerogen
        velBottom = 0.0;
        tempBottom = 0.0;

        if (nMethaneBottom > 0)
        {
            velBottom = vAVG_kerogenBottom * refVelocity / double(nMethaneBottom);
            tempBottom = keKerogenBottom * refVelocity * refVelocity / (3.0 * kB * double(nMethaneBottom));
        }

        velTop = 0.0;
        tempTop = 0.0;

        if (nMethaneTop > 0)
        {
            velTop = vAVG_kerogenTop * refVelocity / double(nMethaneTop);
            tempTop = keKerogenTop * refVelocity * refVelocity / (3.0 * kB * double(nMethaneTop));
        }

        middleFile << t * deltaT * tSkip * refTime << '\t'
                   << rhoN1 << '\t'
                   << rhoN2 << '\t'
                   << p << '\t'
                   << p2 << '\t'
                   << nMethaneMiddle << '\t'
                   << nMethaneMiddle2 << '\t'
                   << vAVG_Middle << '\t'
                   << vAVG_Middle2 << '\t'
                   << massFlowRate << '\t'
                   << massFlowRate2 << '\t'
                   << rhoMiddle << '\t'
                   << tempMiddle << '\t'
                   << endl;

        kerogenFile << t * deltaT * tSkip * refTime << '\t'
                    << nMethaneBottom << '\t'
                    << nMethaneTop << '\t'
                    << tempBottom << '\t'
                    << tempTop << '\t'
                    << velBottom << '\t'
                    << velTop << '\t'
                    << endl;

        getline(data, line);
    }
    return 0;
}
