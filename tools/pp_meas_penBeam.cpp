#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;

vector<double> linspace(double min, double max, size_t N);

int main()
{
    string line, fName;
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
    double vM = sqrt(2 * kB * Tw / mGas); // most probable speed
    double binWidth = 10; // binwidth of velocity
    int maxVout = vM * 4; // max output velocity
    int nBins = ceil(maxVout / binWidth);
    cout << "The number of bins " << nBins << endl;

    double depWidth = 1;
    int nzBins = ceil(rCut / depWidth);
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

    vector<double> zPen; // record the penetration depth of the gas molecule
    vector<double> zMin; // record the minimum value of the single collision process

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

            zMin.resize(nAtoms, 0.0); // initiliase this variable
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

                            zMin[n] = zM[n]; //
                        }
                        // switch into false for re-use this molecule
                        startFromMiddle[n] = false;

                        if (zMin[n] >= zM[n])
                        {
                            zMin[n] = zM[n];
                        }
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

                            zPen.push_back(zMin[n]);

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

    vector<double> zDepth;
    leftFromTop.resize(nWallCollisions, 0);

    vector<double> vN, vNi;
    vector<double> vT, vTi;
    vector<double> vTx, vTxi;
    vector<double> vTy, vTyi;
    double vMagT, vMagTi;

    int coTop = 0, coBottom = 0;
    for (i = 0; i < nWallCollisions; i++)
    {
        if (zPen[i] < rCut)
        {
            zDepth.push_back(zPen[i]);

            // tangential X
            vTx.push_back(vxEnd[i]);
            vTxi.push_back(vxStart[i]);

            // normal
            vN.push_back(abs(vzEnd[i]));
            vNi.push_back(abs(vzStart[i]));
        }
    }

    cout << "Analysis Starts: " << endl;
    cout << "Number of Collisions Started: " << nStartedCollisions << endl;
    cout << "Number of Collisions Ended: " << vzEnd.size() << endl;
    cout << endl;

    int nzPts = zDepth.size();

    // ******************************************************************
    // Beam (partial) tangential velocity distributions depending on the penetration depth
    // ******************************************************************
    vector<double> nDepthID = linspace(-rCut, rCut, 10);
    vector<double> range = linspace(0.20, 2.0, 10);

    ofstream depthBeamDis_Tx("depthBeamDis_Tx.txt");

    vector<vector<vector<double> > > DepthDisPos(nDepthID.size(), vector<vector<double> >(range.size(), vector<double>(nBins)));
    vector<vector<vector<double> > > DepthDisNeg(nDepthID.size(), vector<vector<double> >(range.size(), vector<double>(nBins)));
    vector<vector<vector<double> > > DepthProPos(nDepthID.size(), vector<vector<double> >(range.size(), vector<double>(nBins)));
    vector<vector<vector<double> > > DepthProNeg(nDepthID.size(), vector<vector<double> >(range.size(), vector<double>(nBins)));

    for (int m = 0; m < nDepthID.size(); m++) // loops over the penetration depth
    {
        for (int n = 0; n < range.size(); n++) // loops over the incident velocity magnitude
        {
            vector<double> distributionPos(nBins, 0.0);
            vector<double> distributionNeg(nBins, 0.0);
            vector<double> binsPos(nBins, 0.0);
            vector<double> binsNeg(nBins, 0.0);

            for (i = 0; i < nBins; i++)
            {
                binsPos[i] = binWidth * 0.5 + binWidth * i;
                binsNeg[i] = -binsPos[i];
            }

            for (i = 0; i < nzPts; i++)
            {
                if (abs(vTx[i]) < maxVout)
                {
                    if ((zDepth[i] > nDepthID[m] - rCut / nDepthID.size()) && (zDepth[i] <= nDepthID[m] + rCut / nDepthID.size()))    // select the penetration depth range
                    {
                        if ((vTxi[i] > vM * (range[n] - 0.2)) && (vTxi[i] <= vM * (range[n] + 0.2)))                                  // select the incident velocity range
                        {
                            if (vTx[i] >= 0)
                            {
                                bo = floor(vTx[i] / binWidth);

                                if (bo >= nBins)
                                {
                                    bo = nBins - 1;
                                }

                                distributionPos[bo] += 1.0;
                            }
                            else
                            {
                                bo = floor(abs(vTx[i]) / binWidth);

                                if (bo >= nBins)
                                {
                                    bo = nBins - 1;
                                }

                                distributionNeg[bo] += 1.0;
                            }
                        }
                    }
                }
            }

            vector<double> probabilityPos(nBins, 0.0);
            vector<double> probabilityNeg(nBins, 0.0);

            // scale with area
            double areaUnderGraph = 0.0;

            for (i = 0; i < nBins - 1; i++)
            {
                areaUnderGraph += (binsPos[i + 1] - binsPos[i]) * (distributionPos[i] + distributionPos[i + 1]) * 0.5;
                areaUnderGraph += abs(binsNeg[i + 1] - binsNeg[i]) * (distributionNeg[i] + distributionNeg[i + 1]) * 0.5;
            }

            for (i = 0; i < nBins; i++)
            {
                probabilityPos[i] = distributionPos[i] / areaUnderGraph;
                probabilityNeg[i] = distributionNeg[i] / areaUnderGraph;
            }

            for (int i = 0; i < nBins; i++)
            {
                DepthDisNeg[m][n][i] = distributionNeg[i];
                DepthProNeg[m][n][i] = probabilityNeg[i];
                depthBeamDis_Tx << binsNeg[nBins - i - 1] << " " << DepthDisNeg[m][n][nBins - i - 1] << " " << DepthProNeg[m][n][nBins - i - 1] << endl;
            }

            for (int i = 0; i < nBins; i++)
            {
                DepthDisPos[m][n][i] = distributionPos[i];
                DepthProPos[m][n][i] = probabilityPos[i];
                depthBeamDis_Tx << binsPos[i] << " " << DepthDisPos[m][n][i] << " " << DepthProPos[m][n][i] << endl;
            }
        }
    }
    // ************************ Subfunction Ends ************************

    // ******************************************************************
    // Beam normal velocity distributions depending on the penetration depth
    // ******************************************************************
    ofstream depthBeamDis_N("depthBeamDis_N.txt");

    vector<vector<vector<double> > > DepthDis(nDepthID.size(), vector<vector<double> >(range.size(), vector<double>(nBins)));
    vector<vector<vector<double> > > DepthPro(nDepthID.size(), vector<vector<double> >(range.size(), vector<double>(nBins)));

    for (int m = 0; m < nDepthID.size(); m++)
    {
        for (int n = 0; n < range.size(); n++)
        {
            vector<double> distribution(nBins, 0.0);
            vector<double> bins(nBins, 0.0);

            for (i = 0; i < nBins; i++)
            {
                bins[i] = binWidth * 0.5 + binWidth * i;
            }

            for (i = 0; i < nzPts; i++)
            {
                if ((zDepth[i] > nDepthID[m] - rCut / nDepthID.size()) && (zDepth[i] <= nDepthID[m] + rCut / nDepthID.size()))
                {
                    if ((vNi[i] > vM * (range[n] - 0.2)) && (vNi[i] <= vM * (range[n] + 0.2)))
                    {
                        bo = floor(vN[i] / binWidth);

                        if (bo >= nBins)
                        {
                            bo = nBins - 1;
                        }

                        distribution[bo] += 1.0;
                    }
                }
            }

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

            for (int i = 0; i < nBins; i++)
            {
                DepthDis[m][n][i] = distribution[i];
                DepthPro[m][n][i] = probability[i];
                depthBeamDis_N << bins[i] << " " << DepthDis[m][n][i] << " " << DepthPro[m][n][i] << endl;
            }
        }
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