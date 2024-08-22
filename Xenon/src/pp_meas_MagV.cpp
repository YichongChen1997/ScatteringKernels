#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;

int main()
{
    string line, fName;
    int id, n, i, t, tTime, nAtoms, yStep, currentTimeStep, typ, bo;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r, zMean;
    double xLo, xHi, yLo, yHi, zLo, zHi, Lz;
    double kB = 1.38064852e-23;

    double mi = 2.1802E-25; // mass of one molecule
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

    //**** CHANGE **************

    double Tw = 423;                    // temperature of wall (Kelvin)
    double vM = sqrt(2 * kB * Tw / mi); // most probable speed
    
    ifstream Height("Height.dat", ios::in);
    double H;
    Height >> H;
    double rCut = 15;
    double rCut = 15;

    int nTimeSteps = 2000; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas_gas.lammpstrj | wc -l

    double binWidth = 10; // binwidth of velocity
    int maxVout = vM * 4; // max output velocity
    int nBins = ceil(maxVout / binWidth);

    // **************************

    ifstream data("dump_meas_gas.lammpstrj", ios::in);

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
    vector<int> leftFromTop; // 1 = top; 0 = bottom

    // the following is collision data ("start" = when the molecule enters rCut and
    // "end" = when the molecule has left rCut). Every row is a single atomic collision.
    // Note that one molecule can hit a wall multiple times. so it can be in this list many times
    // Yichong you have all this data you can use for analysis

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
                //                 cout << "nAtoms = " << nAtoms << endl;
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
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> KE >> PE >> tau1 >> tau2 >> tau3;

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
    vector<double> vMag, vMagi;
    double vMagT, vMagTi;

    int coTop = 0, coBottom = 0;

    for (i = 0; i < nWallCollisions; i++)
    {
        if (leftFromTop[i] == 1)
        {

            // normal
            vN.push_back(abs(vzEnd[i]));
            vNi.push_back(abs(vzStart[i]));

            // combined tangents
            vMagT = sqrt(vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i]);
            vMagTi = sqrt(vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i]);

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

            // total velocity
            vMag.push_back(sqrt(abs(vzEnd[i]) * abs(vzEnd[i]) + vMagT * vMagT));
            vMagi.push_back(sqrt(abs(vzStart[i]) * abs(vzStart[i]) + vMagTi * vMagTi));            

            coTop++;
        }

        if (leftFromTop[i] == -1)
        {

            // normal
            vN.push_back(abs(vzEnd[i]));
            vNi.push_back(abs(vzStart[i]));

            // combined tangents
            vMagT = sqrt(vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i]);
            vMagTi = sqrt(vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i]);

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

            // total velocity
            vMag.push_back(sqrt(abs(vzEnd[i]) * abs(vzEnd[i]) + vMagT * vMagT));
            vMagi.push_back(sqrt(abs(vzStart[i]) * abs(vzStart[i]) + vMagTi * vMagTi));

            coBottom++;
        }
    }

    cout << "Size of vN: " << vN.size() << endl;
    cout << "No. of Collisions at top: " << coTop << endl;
    cout << "No. of Collisions at bottom: " << coBottom << endl;

    int nPts = vN.size();

    double sum_vTx = 0.0, sum_vTxi = 0.0, ave_vTx, ave_vTxi;
    for (int i = 0; i < nPts; ++i)
    {
        sum_vTxi += vTxi[i];
        sum_vTx += vTx[i];
    }
    ave_vTxi = sum_vTxi / vTxi.size();
    ave_vTx = sum_vTx / vTx.size();

    cout << "The average incident velocity of tangential X components: " << ave_vTxi << endl;
    cout << "The average reflected velocity of tangential X components: " << ave_vTx << endl;
    cout << endl;

    // Normal 0.5
    {
        // Velocity distribution for selected range (Normal)
        vector<double> distribution(nBins, 0.0);
        vector<double> bins(nBins, 0.0);
        int count = 0;

        for (i = 0; i < nBins; i++)
        {
            bins[i] = binWidth * 0.5 + binWidth * i;
        }

        for (i = 0; i < nPts; i++)
        { 
            // Modify here for different incident velocity magnitude
            if ((vMagi[i] > 0) && (vMagi[i] <= vM * 0.5))
            {
                bo = floor(vN[i] / binWidth);

                if (bo >= nBins)
                {
                    bo = nBins - 1;
                }

                distribution[bo] += 1.0;
                count++;
            }
        }
        cout << "Number of molecules with incident velocity magnitude smaller than 0.5: " << count << endl;

        vector<double> probability(nBins, 0.0);
        double areaUnderGraph = 0.0;
        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }
        ofstream file("vMag_05_N.txt");

        for (i = 0; i < nBins; i++)
        {
            file << bins[i] << " " << distribution[i] << " " << probability[i] << endl;
        }
    }

    // Normal 1.0
    {
        // Velocity distribution for selected range (Normal)
        vector<double> distribution(nBins, 0.0);
        vector<double> bins(nBins, 0.0);
        int count = 0;

        for (i = 0; i < nBins; i++)
        {
            bins[i] = binWidth * 0.5 + binWidth * i;
        }

        for (i = 0; i < nPts; i++)
        { 
            // Modify here for different incident velocity magnitude
            if ((vMagi[i] > 0) && (vMagi[i] <= vM * 1.0))
            {
                bo = floor(vN[i] / binWidth);

                if (bo >= nBins)
                {
                    bo = nBins - 1;
                }

                distribution[bo] += 1.0;
                count++;
            }
        }
        cout << "Number of molecules with incident velocity magnitude smaller than 1.0: " << count << endl;

        vector<double> probability(nBins, 0.0);
        double areaUnderGraph = 0.0;
        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }
        ofstream file("vMag_10_N.txt");

        for (i = 0; i < nBins; i++)
        {
            file << bins[i] << " " << distribution[i] << " " << probability[i] << endl;
        }
    }    

    // Tx 0.5
    {
        // velocity distribution for selected range
        vector<double> distributionA(nBins, 0.0);
        vector<double> distributionB(nBins, 0.0);
        vector<double> binsA(nBins, 0.0);
        vector<double> binsB(nBins, 0.0);
        int count = 0;

        for (i = 0; i < nBins; i++)
        {
            binsA[i] = binWidth * 0.5 + binWidth * i;
            binsB[i] = -binsA[i];
        }

        for (i = 0; i < nPts; i++)
        {
            if ((vMagi[i] > 0) && (vMagi[i] <= vM * 0.5))
            {
                if (vTx[i] >= 0)
                {
                    bo = floor(vTx[i] / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionA[bo] += 1.0;
                }
                else
                {
                    bo = floor(abs(vTx[i]) / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionB[bo] += 1.0;
                }
                count++;
            }
        }
        cout << "Number of molecules with incident velocity magnitude smaller than 0.5: " << count << endl;
        vector<double> probabilityA(nBins, 0.0);
        vector<double> probabilityB(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (binsA[i + 1] - binsA[i]) * (distributionA[i] + distributionA[i + 1]) * 0.5;
            areaUnderGraph += abs(binsB[i + 1] - binsB[i]) * (distributionB[i] + distributionB[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probabilityA[i] = distributionA[i] / areaUnderGraph;
            probabilityB[i] = distributionB[i] / areaUnderGraph;
        }

        ofstream file("vMag_05_Tx.txt");

        for (i = 0; i < nBins; i++)
        {
            file << binsB[nBins - i - 1] << " " << distributionB[nBins - i - 1] << " " << probabilityB[nBins - i - 1] << endl;
        }

        for (i = 0; i < nBins; i++)
        {
            file << binsA[i] << " " << distributionA[i] << " " << probabilityA[i] << endl;
        }
    }

    // Tx 1.0
    {
        // velocity distribution for selected range
        vector<double> distributionA(nBins, 0.0);
        vector<double> distributionB(nBins, 0.0);
        vector<double> binsA(nBins, 0.0);
        vector<double> binsB(nBins, 0.0);
        int count = 0;

        for (i = 0; i < nBins; i++)
        {
            binsA[i] = binWidth * 0.5 + binWidth * i;
            binsB[i] = -binsA[i];
        }

        for (i = 0; i < nPts; i++)
        {
            if ((vMagi[i] > 0) && (vMagi[i] <= vM * 1.0))
            {
                if (vTx[i] >= 0)
                {
                    bo = floor(vTx[i] / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionA[bo] += 1.0;
                }
                else
                {
                    bo = floor(abs(vTx[i]) / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionB[bo] += 1.0;
                }
                count++;
            }
        }
        cout << "Number of molecules with incident velocity magnitude smaller than 1.0: " << count << endl;
        vector<double> probabilityA(nBins, 0.0);
        vector<double> probabilityB(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (binsA[i + 1] - binsA[i]) * (distributionA[i] + distributionA[i + 1]) * 0.5;
            areaUnderGraph += abs(binsB[i + 1] - binsB[i]) * (distributionB[i] + distributionB[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probabilityA[i] = distributionA[i] / areaUnderGraph;
            probabilityB[i] = distributionB[i] / areaUnderGraph;
        }

        ofstream file("vMag_10_Tx.txt");

        for (i = 0; i < nBins; i++)
        {
            file << binsB[nBins - i - 1] << " " << distributionB[nBins - i - 1] << " " << probabilityB[nBins - i - 1] << endl;
        }

        for (i = 0; i < nBins; i++)
        {
            file << binsA[i] << " " << distributionA[i] << " " << probabilityA[i] << endl;
        }
    }    


    return 0;
}
