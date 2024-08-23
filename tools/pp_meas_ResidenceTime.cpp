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
    string line;
    int id, n, i, t, nAtoms, currentTimeStep, typ, bo;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, fx, fy, fz;
    double xLo, xHi, yLo, yHi, zLo, zHi, Lz;

    double kB = 1.38064852e-23;
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

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
    double vM = sqrt(2 * kB * Tw / mGas); // most probable speed
    double binWidth = 1;   // binwidth of time
    int maxDis = 200;   
    int nBins = ceil(maxDis / binWidth);

    ifstream data("dump_meas_gas.lammpstrj", ios::in);
    // **************************

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
    // Note that one molecule can hit a wall multiple times. so it can be in this list many times.

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
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3 >> KE >> PE >> fx >> fy >> fz;

            if (typ == 1) // methane
            {
                if (t == 0)
                {
                    if (z >= rCut)
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
                    if (z < rCut)
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

                    else if (z >= rCut)
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

    cout << "Analysis Starts: " << endl;
    cout << "Number of Collisions Started: " << nStartedCollisions << endl;
    cout << "Number of Collisions Ended: " << vzEnd.size() << endl;
    cout << endl;

    vector<double> tIn, tOut;

    int coTop = 0, coBottom = 0;
    double sumTime = 0.0;

    for (i = 0; i < nWallCollisions; i++)
    {

        tIn.push_back(tStart[i]);
        tOut.push_back(tEnd[i]);

        coBottom++;
    
        sumTime += (tEnd[i] - tStart[i]);
    }

    cout << "Size of tIn:  " << tIn.size() << endl;
    cout << "Size of tOut: " << tOut.size() << endl;

    cout << "No. of Collisions at top: " << coTop << endl;
    cout << "No. of Collisions at bottom: " << coBottom << endl;

    int nPts = tOut.size();

    // ************************ Residence Time Distribution *******************************
    {
        double tMax = 0.0;
        vector<double> distribution(nBins, 0.0);
        vector<double> bins(nBins, 0.0);

        for (i = 0; i < nBins; i++)
        {
            bins[i] = binWidth * 0.5 + binWidth * i;
        }

        for (i = 0; i < nPts; i++)
        {
            if ( (tOut[i] - tIn[i]) < maxDis)
            {
                bo = floor( (tOut[i] - tIn[i]) / binWidth);

                if (bo >= nBins)
                {
                    bo = nBins - 1;
                }

                distribution[bo] += 1.0;
            }

            if ((tOut[i] - tIn[i]) > tMax)
            {
                tMax = tOut[i] - tIn[i];
            }
        }

        cout << "Maximum abosorped time = " << tMax * deltaT * tSkip * refTime << endl;

        vector<double> probability(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }

        ofstream file("timeDis.txt");

        for (i = 0; i < nBins; i++)
        {
            file << bins[i] * deltaT * tSkip * refTime << " " << distribution[i] << " " << probability[i] << " " <<  sumTime * deltaT * tSkip * refTime/nPts << endl;
        }
    }
    return 0;
}
