#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;
const float PI = 3.1415;

vector<double> linspace(double min, double max, size_t N);

int main()
{
    string line, fName;
    int id, n, i, t, tTime, nAtoms, yStep, currentTimeStep, typ, bo;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r, zMean;
    double xLo, xHi, yLo, yHi, zLo, zHi, Lz;
    double kB = 1.38064852e-23;

    double mi = 2.1802E-25;   // mass of one Argon
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

    double mole = 6.02214076e23;
    double KcalJ = 4184;
    double refEnergy = KcalJ / mole;

    //**** CHANGE **************
    double Tw = 423;                    // temperature of wall (Kelvin)
    double vM = sqrt(2 * kB * Tw / mi); // most probable speed
    double MCsigma = 4.10;              // diameter of methane

    ifstream Height("Height.dat", ios::in);
    double H;
    Height >> H;
    double rCut = 15;

    ifstream nSteps("nTimeSteps.dat", ios::in);
    int nTimeSteps; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l
    nSteps >> nTimeSteps;

    // **************************
    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    ofstream file1("NMACs.txt", ios::out);
    ofstream file2("TMACs.txt", ios::out);
    ofstream file3("NKEACs.txt", ios::out);

    int nWallCollisions = 0;    // collect total number of collided molecules with wall
    int nStartedCollisions = 0; // collect number of collisions that have started
    vector<int> nCollisions;    // per molecule collect how many times molecule collided

    vector<bool> startFromMiddle;
    vector<bool> crossedStart;

    // keep track of molecule positions, velocities and time
    // only when they are in the middle of the domain (stop when molecules leave)
    vector<double> tM;       // time
    vector<double> xM;       // position x
    vector<double> yM;       // position y
    vector<double> zM;       // position z
    vector<double> vxM;      // velocity x
    vector<double> vyM;      // velocity y
    vector<double> vzM;      // velocity z
    vector<double> keM;      // kinetic energy
    vector<double> peM;      // potential energy
    vector<int> leftFromTop; // 1 = top; 0 = bottom

    // the following is collision data ("start" = when the molecule enters rCut and
    // "end" = when the molecule has left rCut).
    // Note that one molecule can hit a wall multiple times. so it can be in this list many times

    vector<double> tEnter;
    vector<double> xEnter;
    vector<double> yEnter;
    vector<double> zEnter;
    vector<double> vxEnter;
    vector<double> vyEnter;
    vector<double> vzEnter;

    vector<double> tStart;
    vector<double> xStart;
    vector<double> yStart;
    vector<double> zStart;
    vector<double> vxStart;
    vector<double> vyStart;
    vector<double> vzStart;

    vector<double> tEnd;
    vector<double> xEnd;
    vector<double> yEnd;
    vector<double> zEnd;
    vector<double> vxEnd;
    vector<double> vyEnd;
    vector<double> vzEnd;

    vector<int> leftFromTheTop;

    for (t = 0; t < nTimeSteps; t++)
    {
        for (n = 1; n < 10; n++)
        {
            if (n == 4)
            {
                data >> nAtoms;
            }

            if (n == 2)
            {
                data >> currentTimeStep;

                cout << "currentTimeStep = " << currentTimeStep
                     << "; t = " << t << " [ "
                     << 100 * float(t + 1) / float(nTimeSteps)
                     << "% ]" << endl;
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

        // set fields
        if (t == 0)
        {
            nCollisions.resize(nAtoms, 0);
            startFromMiddle.resize(nAtoms, false);
            crossedStart.resize(nAtoms, false);

            tM.resize(nAtoms, 0);
            xM.resize(nAtoms, 0.0);
            yM.resize(nAtoms, 0.0);
            zM.resize(nAtoms, 0.0);
            vxM.resize(nAtoms, 0.0);
            vyM.resize(nAtoms, 0.0);
            vzM.resize(nAtoms, 0.0);

            tEnter.resize(nAtoms, 0);
            xEnter.resize(nAtoms, 0.0);
            yEnter.resize(nAtoms, 0.0);
            zEnter.resize(nAtoms, 0.0);
            vxEnter.resize(nAtoms, 0.0);
            vyEnter.resize(nAtoms, 0.0);
            vzEnter.resize(nAtoms, 0.0);
        }

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3;

            if (typ == 1) // methane
            {
                if (t == 0)
                {
                    if ((z >= rCut) && (z <= H - rCut))
                    {
                        startFromMiddle[n] = true; // inside middle region
                    }
                    tM[n] = t;

                    xM[n] = x;
                    yM[n] = y;
                    zM[n] = z;

                    vxM[n] = vx * refVelocity;
                    vyM[n] = vy * refVelocity;
                    vzM[n] = vz * refVelocity;
                }
                else
                {
                    if ((z < rCut) || (z > H - rCut))
                    {
                        if (startFromMiddle[n])
                        {
                            // the collision has started
                            crossedStart[n] = true;

                            tEnter[n] = tM[n];

                            xEnter[n] = xM[n];
                            yEnter[n] = yM[n];
                            zEnter[n] = zM[n];

                            vxEnter[n] = vxM[n];
                            vyEnter[n] = vyM[n];
                            vzEnter[n] = vzM[n];

                            nStartedCollisions++;
                        }
                        // switch into false for re-use this molecule
                        startFromMiddle[n] = false;
                    }

                    else if ((z >= rCut) && (z <= H - rCut)) // YC: the middle region
                    {
                        if (crossedStart[n])
                        {
                            // push data when the molecule started the collision
                            tStart.push_back(tEnter[n]);

                            xStart.push_back(xEnter[n]);
                            yStart.push_back(yEnter[n]);
                            zStart.push_back(zEnter[n]);

                            vxStart.push_back(vxEnter[n]);
                            vyStart.push_back(vyEnter[n]);
                            vzStart.push_back(vzEnter[n]);

                            /* if the molecule is within the middle region and has started the collision,
                             flag it as the END */
                            tEnd.push_back(t);

                            xEnd.push_back(x);
                            yEnd.push_back(y);
                            zEnd.push_back(z);

                            vxEnd.push_back(vx * refVelocity);
                            vyEnd.push_back(vy * refVelocity);
                            vzEnd.push_back(vz * refVelocity);

                            // switch into false for re-use this molecule
                            crossedStart[n] = false;

                            nCollisions[n]++;
                            nWallCollisions++;
                        }
                        // switch into true for re-use this molecule
                        startFromMiddle[n] = true;
                    }

                    tM[n] = t;

                    xM[n] = x;
                    yM[n] = y;
                    zM[n] = z;

                    vxM[n] = vx * refVelocity;
                    vyM[n] = vy * refVelocity;
                    vzM[n] = vz * refVelocity;

                }
            }
        }

        getline(data, line);

        cout << "number of gas-wall collisions = " << nWallCollisions << endl;
        cout << endl;
    }

    leftFromTop.resize(nWallCollisions, 0);
    for (i = 0; i < nWallCollisions; i++)
    {
        // to indicate which direction the molecule left the middle (top or bottom)
        if (zStart[i] + zEnd[i] >= H)
        {
            leftFromTop[i] = 1;
        }
        else if (zStart[i] + zEnd[i] < H)
        {
            leftFromTop[i] = -1;
        }
    }

    cout << "Analysis Starts: " << endl;
    cout << "Number of Collisions Started: " << nStartedCollisions << endl;
    cout << "Number of Collisions Ended: " << vzEnd.size() << endl;
    cout << endl;

    vector<double> vN, vNi;
    vector<double> vT, vTi;
    vector<double> vTx, vTxi;
    vector<double> vTy, vTyi;
    double vMagT, vMagTi;
    double vMag, vMagi;

    int coTop = 0, coBottom = 0;

    for (i = 0; i < nWallCollisions; i++)
    {
        if (leftFromTop[i] == 1)
        {
            // normal
            vN.push_back(abs(vzEnd[i]));
            vNi.push_back(abs(vzStart[i]));

            // combined tangents
            vMagT = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i];
            vMagTi = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i];

            vMag = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i] + vzEnd[i] * vzEnd[i];
            vMagi = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i] + vzStart[i] + vzStart[i];

            if (vMagT > 0)
            {
                vT.push_back(sqrt(vMagT));
            }
            else
            {
                vT.push_back(0.0);
            }

            if (vMagTi > 0)
            {
                vTi.push_back(sqrt(vMagTi));
            }
            else
            {
                vTi.push_back(0.0);
            }

            vTx.push_back(vxEnd[i]);
            vTxi.push_back(vxStart[i]);

            vTy.push_back(vyEnd[i]);
            vTyi.push_back(vyStart[i]);

            coTop++;
        }

        if (leftFromTop[i] == -1)
        {
            // normal
            vN.push_back(abs(vzEnd[i]));
            vNi.push_back(abs(vzStart[i]));

            // combined tangents
            vMagT = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i];
            vMagTi = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i];

            vMag = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i] + vzEnd[i] * vzEnd[i];
            vMagi = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i] + vzStart[i] + vzStart[i];

            if (vMagT > 0)
            {
                vT.push_back(sqrt(vMagT));
            }
            else
            {
                vT.push_back(0.0);
            }

            if (vMagTi > 0)
            {
                vTi.push_back(sqrt(vMagTi));
            }
            else
            {
                vTi.push_back(0.0);
            }

            vTx.push_back(vxEnd[i]);
            vTxi.push_back(vxStart[i]);

            vTy.push_back(vyEnd[i]);
            vTyi.push_back(vyStart[i]);

            coBottom++;
        }
    }

    cout << "No. of Collisions at top: " << coTop << endl;
    cout << "No. of Collisions at bottom: " << coBottom << endl;

    int nPts = vN.size();

    vector<double> range = linspace(0.1, 3.9, 39);

    for (int n = 0; n < range.size(); n++)
    {
        double vTx_full = 0.0, vTxi_full = 0.0, vTx_meanF, vTxi_meanF;
        double vN_full = 0.0, vNi_full = 0.0, vN_meanF, vNi_meanF;

        double vTx_partial = 0.0, vTxi_partial = 0.0, vTx_meanP, vTxi_meanP;
        double vN_partial = 0.0, vNi_partial = 0.0, vN_meanP, vNi_meanP;

        double sigmaTx_f1 = 0.0, sigmaTx_f2;
        double sigmaTx_p1 = 0.0, sigmaTx_p2;
        double sigmaN_f1 = 0.0, sigmaN_f2;
        double sigmaN_p1 = 0.0, sigmaN_p2;
        double alphaN_p1 = 0.0, alphaN_p2;
        double alphaN_f1 = 0.0, alphaN_f2;

        int countN = 0, countTx = 0, countPositiveN = 0, countPositiveTx = 0;

        for (int i = 0; i < nPts; ++i)
        {

            if ((vNi[i] >= vM * 0.0) && (vNi[i] < vM * 6))
            {
                vNi_full += vNi[i];
                vN_full += vN[i];
                countPositiveN++;
            }

            if (vTxi[i] >= 0.0)
            {
                vTxi_full += vTxi[i];
                vTx_full += vTx[i];
                countPositiveTx++;
            }

            if ((vNi[i] > vM * (range[n] - 0.1)) && (vNi[i] <= vM * (range[n] + 0.1)))
            {
                vNi_partial += vNi[i];
                vN_partial += vN[i];
                countN++;
            }

            if ((vTxi[i] > vM * (range[n] - 0.1)) && (vTxi[i] <= vM * (range[n] + 0.1)))
            {
                vTxi_partial += vTxi[i];
                vTx_partial += vTx[i];
                countTx++;
            }
        }

        vNi_meanP = vNi_partial / countN; // mean value of incident normal velocities within a partial range
        vN_meanP = vN_partial / countN;
        vTxi_meanP = vTxi_partial / countTx;
        vTx_meanP = vTx_partial / countTx;

        vTxi_meanF = vTxi_full / countPositiveTx; // mean value of incident tangential velocities for full range.
        vTx_meanF = vTx_full / countPositiveTx;
        vNi_meanF = vNi_full / countPositiveN;
        vN_meanF = vN_full / countPositiveN;

        // ********Tangential Momentum Accommodation Coefficients***********

        // Partial TMAC
        /* Method One (General expression) */
        sigmaTx_p1 = (vTxi_partial - vTx_partial) / vTxi_partial;

        /* Method Two (Least-Squares fitting) */
        double beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            if ((vTxi[i] > vM * (range[n] - 0.1)) && (vTxi[i] <= vM * (range[n] + 0.1)))
            {
                beta_nu += (vTxi[i] - vTxi_meanP) * (vTx[i] - vTx_meanP);   // beta numerator
                beta_de += (vTxi[i] - vTxi_meanP) * (vTxi[i] - vTxi_meanP); // denominator
            }
        }
        sigmaTx_p2 = 1 - beta_nu / beta_de;

        // Global TMAC
        /* Method One (General expression) */
        sigmaTx_f1 = (vTxi_meanF - vTx_meanF) / vTxi_meanF;

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            beta_nu += (vTxi[i] - vTxi_meanF) * (vTx[i] - vTx_meanF);   // beta numerator
            beta_de += (vTxi[i] - vTxi_meanF) * (vTxi[i] - vTxi_meanF); // beta donominator
        }
        sigmaTx_f2 = 1 - beta_nu / beta_de;

        // ********Normal Momentum Accommodation Coefficients***************

        // Partial NMAC
        /* Method One (General expression) */
        sigmaN_p1 = (vNi_partial - vN_partial) / (vNi_partial - sqrt(PI) / 2 * vM * countN);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            if ((vNi[i] > vM * (range[n] - 0.1)) && (vNi[i] <= vM * (range[n] + 0.1)))
            {
                beta_nu += (vNi[i] - vNi_meanP) * (vN[i] - vN_meanP);   // beta numerator
                beta_de += (vNi[i] - vNi_meanP) * (vNi[i] - vNi_meanP); // denominator
            }
        }
        sigmaN_p2 = 1 - beta_nu / beta_de;

        // Global NMAC
        /* Method One (General expression) */
        sigmaN_f1 = (vNi_meanF - vN_meanF) / (vNi_meanF - sqrt(PI) * vM / 2);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            beta_nu += (vNi[i] - vNi_meanF) * (vN[i] - vN_meanF);   // beta numerator
            beta_de += (vNi[i] - vNi_meanF) * (vNi[i] - vNi_meanF); // beta donominator
        }
        sigmaN_f2 = 1 - beta_nu / beta_de;

        // ********Kinetic Energy Accommodation Coefficients (Normal component)***************

        // Partial NKEAC
        /* Method One (General expression) */
        alphaN_p1 = (vNi_partial * vNi_partial - vN_partial * vN_partial) / (vNi_partial * vNi_partial - (sqrt(PI) / 2 * vM * countN) * (sqrt(PI) / 2 * vM * countN));

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            if ((vNi[i] > vM * (range[n] - 0.1)) && (vNi[i] <= vM * (range[n] + 0.1)))
            {
                beta_nu += (vNi[i] * vNi[i] - vNi_meanP * vNi_meanP) * (vN[i] * vN[i] - vN_meanP * vN_meanP);     // beta numerator
                beta_de += (vNi[i] * vNi[i] - vNi_meanP * vNi_meanP) * (vNi[i] * vNi[i] - vNi_meanP * vNi_meanP); // denominator
            }
        }
        alphaN_p2 = 1 - beta_nu / beta_de;

        // Global NKEAC
        /* Method One (General expression) */
        alphaN_f1 = (vNi_meanF * vNi_meanF - vN_meanF * vN_meanF) / (vNi_meanF * vNi_meanF - (sqrt(PI) * vM / 2) * (sqrt(PI) * vM / 2));

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            beta_nu += (vNi[i] * vNi[i] - vNi_meanF * vNi_meanF) * (vN[i] * vN[i] - vN_meanF * vN_meanF);     // beta numerator
            beta_de += (vNi[i] * vNi[i] - vNi_meanF * vNi_meanF) * (vNi[i] * vNi[i] - vNi_meanF * vNi_meanF); // beta donominator
        }
        alphaN_f2 = 1 - beta_nu / beta_de;

        file1 << range[n] << " " << sigmaN_p1 << " " << sigmaN_p2 << " " << sigmaN_f1 << " " << sigmaN_f2 << endl;
        file2 << range[n] << " " << sigmaTx_p1 << " " << sigmaTx_p2 << " " << sigmaTx_f1 << " " << sigmaTx_f2 << endl;
        file3 << range[n] << " " << alphaN_p1 << " " << alphaN_p2 << " " << alphaN_f1 << " " << alphaN_f2 << endl;
    }

    return 0;
}

// Create a vector of evenly spaced numbers.
vector<double> linspace(double min, double max, size_t N)
{
    vector<double> linspace;
    double delta = (max - min) / double(N - 1);
    for (int i = 0; i < N; i++)
    {
        linspace.push_back(min + i * delta);
    }
    return linspace;
}
