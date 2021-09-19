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

    double mi = 2.66389e-26; // mass of one molecule
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
    double MCsigma = 3.545;             // diameter of methane

    ifstream Height("Height.dat", ios::in);
    double H;
    Height >> H;
    double rCut = 15;

    ifstream nSteps("nTimeSteps.dat", ios::in);
    int nTimeSteps; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l
    nSteps >> nTimeSteps;

    // **************************
    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    ofstream output("EACs.txt", ios::out);

    int nWallCollisions = 0;    // collect total number of collided molecules with wall
    int nStartedCollisions = 0; // collect number of collisions that have started
    vector<int> nCollisions;    // per molecule collect how many times molecule collided (collision distribution)

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

    for (t = 0; t < nTimeSteps; t++)
    {
        // Header post-process (general to all the dump files)
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

        // set initial fields
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

    vector<double> ke, kei;

    int coTop = 0, coBottom = 0;

    for (i = 0; i < nWallCollisions; i++)
    {
        if (leftFromTop[i] == 1)
        {
            ke.push_back(vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i] + vzEnd[i] * vzEnd[i]);
            kei.push_back(vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i] + vzStart[i] * vzStart[i]);

            coTop++;
        }

        if (leftFromTop[i] == -1)
        {
            ke.push_back(vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i] + vzEnd[i] * vzEnd[i]);
            kei.push_back(vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i] + vzStart[i] * vzStart[i]);

            coBottom++;
        }
    }

    cout << "Size of ke: " << ke.size() << endl;
    cout << "Size of kei: " << kei.size() << endl;
    cout << "No. of collisions at top: " << coTop << endl;
    cout << "No. of collisions at bottom: " << coBottom << endl;

    int nPts = ke.size();

    vector<double> range = linspace(0.1, 1.9, 19);

    for (int n = 0; n < range.size(); n++)
    {
        double EAC_f_v1, EAC_f_v2;
        double EAC_p_v1, EAC_p_v2;

        double KEi_f_all = 0.0, KE_f_all = 0.0;
        double KEi_f_mean, KE_f_mean;

        double KEi_p_all = 0.0, KE_p_all = 0.0;
        double KEi_p_mean, KE_p_mean;

        int countEF = 0, countEP = 0;

        double beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {

            // summation of the total kinetic energy
            KEi_f_all += kei[i];
            KE_f_all += ke[i];
            countEF++;

            if ((kei[i] > vM * vM * (range[n] - 0.1)) && (kei[i] <= vM * vM * (range[n] + 0.1)))
            {
                KEi_p_all += kei[i];
                KE_p_all += ke[i];
                countEP++;
            }
        }
        // Global EAC

        /* Method One */
        EAC_f_v1 = (KEi_f_all - KE_f_all) / (KEi_f_all - (4 * kB * Tw / mi) * countEF);
        cout << "Value of total incident energy: " << KEi_f_all << endl;
        cout << "Value of total reflected energy: " << KE_f_all << endl;
        cout << "Expected total energy: " << (4 * kB * Tw / mi) * countEF << endl;
        cout << "Energy Accommodation coefficient is: " << EAC_f_v1 << endl;

        /* Method Two */
        cout << "No. of summed data: " << countEF << endl;
        KEi_f_mean = KEi_f_all / countEF;
        KE_f_mean = KE_f_all / countEF;
        for (int i = 0; i < nPts; ++i)
        {
            beta_nu += (kei[i] - KEi_f_mean) * (ke[i] - KE_f_mean);   // beta numerator
            beta_de += (kei[i] - KEi_f_mean) * (kei[i] - KEi_f_mean); // beta denominator
        }
        EAC_f_v2 = 1 - beta_nu / beta_de;

        // Partial EAC
        /* Method One */
        cout << "countEP is: " << countEP << endl;
        KEi_p_mean = KEi_p_all / countEP;
        KE_p_mean = KE_p_all / countEP;

        EAC_p_v1 = (KEi_p_mean - KE_p_mean) / (KEi_p_mean - (4 * kB * Tw / mi));


        /* Method Two */
        beta_nu = 0.0, beta_de = 0.0;
        for (int i = 0; i < nPts; i++)
        {
            if ((kei[i] > vM * vM * (range[n] - 0.1)) && (kei[i] <= vM * vM * (range[n] + 0.1)))
            {
            beta_nu += (kei[i] - KEi_p_mean) * (ke[i] - KE_p_mean);
            beta_de += (kei[i] - KEi_p_mean) * (kei[i] - KEi_p_mean);
            }
        }
        EAC_p_v2 = 1 - beta_nu / beta_de;

        output << range[n] << " " << EAC_p_v1 << " " << EAC_p_v2 << " " << EAC_f_v1 << " " << EAC_f_v2 << endl;
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
