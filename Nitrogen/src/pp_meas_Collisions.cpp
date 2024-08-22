#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;
const float PI = 3.1415;

// Determine whether the two numbers have same sign
bool sameSign(double num1, double num2);

int main()
{
    string line, fName;
    int id, n, i, t, tTime, nAtoms, yStep, currentTimeStep, typ, bo;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r, zMean;
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
    double rCut, Tw, H, Msigma, Mass;

    ifstream Specification("Specification.dat", ios::in);
    getline(Specification, line);
    Specification >> nTimeSteps  >> rCut >> Tw  >>  H  >>  Msigma >> Mass;
    // **************************

    double vM = sqrt(2 * kB * Tw / Mass); // most probable speed

    ifstream data("dump_meas_gas.lammpstrj", ios::in);
    ofstream collisionDis("collisionDis.txt");
    ofstream collisionAlphaTx("collisionAlphaTx.txt", ios::out);
    ofstream collisionAlphaEn("collisionAlphaEn.txt", ios::out);

    int nWallCollisions = 0;    // collect total number of collided molecules with wall
    int nStartedCollisions = 0; // collect number of collisions that have started
    vector<int> nCollisions;    // number of collisions per molecule

    vector<bool> startFromMiddle;
    vector<bool> crossedStart;

    // keep track of molecule positions, velocities and time
    // over the entire simulation domain
    vector<double> tM;       // time
    vector<double> xM;       // position x
    vector<double> yM;       // position y
    vector<double> zM;       // position z
    vector<double> vxM;      // velocity x
    vector<double> vyM;      // velocity y
    vector<double> vzM;      // velocity z
    vector<int> leftFromTop; // 1 = top; 0 = bottom
    vector<int> rCutCollisions;

    // the following is collision data ("start" = when the molecule enters rCut and
    // "end" = when the molecule has left rCut). Every row is a single atomic collision.
    // Note that one molecule can hit a wall multiple times. so it can be in this list many times
    vector<double> tEnter;
    vector<double> xEnter;
    vector<double> yEnter;
    vector<double> zEnter;
    vector<double> vxEnter;
    vector<double> vyEnter;
    vector<double> vzEnter;

    vector<double> vxDir;
    vector<double> vyDir;
    vector<double> vzDir; // keep track of the direction of molecules

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
    vector<int> AtomCollisions;

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

            rCutCollisions.resize(nAtoms, 0);
            vxDir.resize(nAtoms, 0);
            vyDir.resize(nAtoms, 0);
            vzDir.resize(nAtoms, 0);

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

                            vxDir[n] = vxM[n];
                            vyDir[n] = vyM[n];
                            vzDir[n] = vzM[n];

                            nStartedCollisions++;

                            // switch into false for re-use this molecule
                            startFromMiddle[n] = false;
                        }

                        else if ((crossedStart[n]) && (!startFromMiddle[n]))
                        {
                            if ((!sameSign(vxDir[n], vxM[n])) || (!sameSign(vyDir[n], vyM[n])) || (!sameSign(vzDir[n], vzM[n])))
                            {
                                rCutCollisions[n]++;
                            }
                            if (!sameSign(vxDir[n], vxM[n]))
                            {
                                vxDir[n] = -vxDir[n];
                            }
                            if (!sameSign(vyDir[n], vyM[n]))
                            {
                                vyDir[n] = -vyDir[n];
                            }
                            if (!sameSign(vzDir[n], vzM[n]))
                            {
                                vzDir[n] = -vzDir[n];
                            }
                        }
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

                            AtomCollisions.push_back(rCutCollisions[n]);
                            rCutCollisions[n] = 0;

                            nCollisions[n]++;
                            nWallCollisions++;
                        }
                        // switch into true for using the molecules in the middle region
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

    cout << "------------------------------------------------------------" << endl;
    cout << "Analysis Starts: " << endl;
    cout << "Number of Collisions Started: " << nStartedCollisions << endl;
    cout << "Number of Collisions Ended: " << vzEnd.size() << endl;
    cout << endl;

    vector<double> vTx, vTxi;
    vector<double> vN, vNi;
    vector<double> keTx, keTxi;
    vector<double> keN, keNi;
    vector<double> ke, kei;

    int coTop = 0, coBottom = 0;

    for (i = 0; i < nWallCollisions; i++)
    {
        if (leftFromTop[i] == 1)
        {
            // tangential momentum
            vTx.push_back(vxEnd[i]);
            vTxi.push_back(vxStart[i]);

            // normal momentum
            vN.push_back(abs(vzEnd[i]));
            vNi.push_back(abs(vzStart[i]));

            // tangential kinetic energy
            keTx.push_back(vxEnd[i] * vxEnd[i]);
            keTxi.push_back(vxStart[i] * vxStart[i]);

            // normal kinetic energy
            keN.push_back(vzEnd[i] * vzEnd[i]);
            keNi.push_back(vzStart[i] * vzStart[i]);

            // total kinetic energy
            ke.push_back(vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i] + vzEnd[i] * vzEnd[i]);
            kei.push_back(vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i] + vzStart[i] * vzStart[i]);

            coTop++;            
        }

        if (leftFromTop[i] == -1)
        {
            // tangential momentum
            vTx.push_back(vxEnd[i]);
            vTxi.push_back(vxStart[i]);

            // normal momentum
            vN.push_back(abs(vzEnd[i]));
            vNi.push_back(abs(vzStart[i]));

            // tangential kinetic energy
            keTx.push_back(vxEnd[i] * vxEnd[i]);
            keTxi.push_back(vxStart[i] * vxStart[i]);

            // normal kinetic energy
            keN.push_back(vzEnd[i] * vzEnd[i]);
            keNi.push_back(vzStart[i] * vzStart[i]);

            // total kinetic energy
            ke.push_back(vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i] + vzEnd[i] * vzEnd[i]);
            kei.push_back(vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i] + vzStart[i] * vzStart[i]);

            coBottom++;
        }
    }

    cout << "No. of Collisions at top: " << coTop << endl;
    cout << "No. of Collisions at bottom: " << coBottom << endl;

    int nPts = ke.size();

    double CoBinWidth = 1;        // binwidth of collision distribution
    int CoMax = 20;               // Max times of collisions
    int CoBins = CoMax / CoBinWidth;

    // Distribution of Collision Times
    {
        vector<double> distribution(CoBins, 0.0);
        vector<double> bins(CoBins, 0.0);

        for (i = 0; i < CoBins; i++)
        {
            bins[i] = CoBinWidth * i;
        }

        for (i = 0; i < AtomCollisions.size(); i++)
        {
            bo = AtomCollisions[i] / CoBinWidth;

            if (bo >= CoBins)
            {
                bo = CoBins - 1;
            }

            distribution[bo] += 1.0;
        }

        vector<double> probability(CoBins, 0.0);
        double areaUnderGraph = 0.0;
        for (i = 0; i < CoBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < CoBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }

        for (i = 0; i < CoBins; i++)
        {
            collisionDis << bins[i] << " " << distribution[i] << " " << probability[i] << endl;
        }
    }
    cout << "Collision Distribution Finished" << endl;

    int collisionPts = 20;
    vector<double> alphaTx_v1, alphaTx_v2;
    vector<double> vTx_mean, vTxi_mean;
    alphaTx_v1.resize(collisionPts, 0.0);
    alphaTx_v2.resize(collisionPts, 0.0);
    vTx_mean.resize(collisionPts, 0.0);
    vTxi_mean.resize(collisionPts, 0.0);

    for (n = 0; n < collisionPts; n++)
    {
        double beta_nu = 0.0, beta_de = 0.0;
        int countPositiveTx = 0;

        for (int i = 0; i < nPts; i++)
        {
            if (AtomCollisions[i] == n)
            {
                if (vxStart[i] > 0)
                {
                    vTxi_mean[n] += vxStart[i];
                    vTx_mean[n] += vxEnd[i];
                    countPositiveTx++;
                }
            }
        }
        vTxi_mean[n] = vTxi_mean[n] / countPositiveTx;
        vTx_mean[n] = vTx_mean[n] / countPositiveTx;

        for (int i = 0; i < nPts; i++)
        {
            if (AtomCollisions[i] == n)
            {
                beta_nu += (vxStart[i] - vTxi_mean[n]) * (vxEnd[i] - vTx_mean[n]);
                beta_de += (vxStart[i] - vTxi_mean[n]) * (vxStart[i] - vTxi_mean[n]);
            }
        }
        alphaTx_v1[n] = 1 - beta_nu / beta_de;
        
        alphaTx_v2[n] = (vTxi_mean[n] - vTx_mean[n]) / vTxi_mean[n];

        collisionAlphaTx << n << " " << alphaTx_v1[n] << " " << alphaTx_v2[n] << endl;
    }

    vector<double> alphaEn_v1, alphaEn_v2;
    vector<double> En_mean, Eni_mean;
    alphaEn_v1.resize(collisionPts, 0.0);
    alphaEn_v2.resize(collisionPts, 0.0);
    En_mean.resize(collisionPts, 0.0);
    Eni_mean.resize(collisionPts, 0.0);

    for (n = 0; n < collisionPts; n++)
    {
        double beta_nu = 0.0, beta_de = 0.0;
        int countPositiveN = 0;

        for (int i = 0; i < nPts; i++)
        {
            if (AtomCollisions[i] == n)
            {
                if (vzStart[i] < 0)
                {
                    Eni_mean[n] += keNi[i];
                    En_mean[n]  += keN[i];
                    countPositiveN++;
                }
            }
        }

        Eni_mean[n] = Eni_mean[n] / countPositiveN;
        En_mean[n]  = En_mean[n] / countPositiveN;

        alphaEn_v1[n] = (Eni_mean[n] -  En_mean[n]) / (Eni_mean[n]  - vM * vM);

        for (int i = 0; i < nPts; i++)
        {
            if (AtomCollisions[i] == n)
            {
                beta_nu += (keNi[i] - Eni_mean[n]) * (keN[i]  - En_mean[n]);    // beta numerator
                beta_de += (keNi[i] - Eni_mean[n]) * (keNi[i] - Eni_mean[n]);   // denominator
            }
        }

        alphaEn_v2[n] = 1 - beta_nu / beta_de;

        collisionAlphaEn << n << " " << alphaEn_v1[n] << " " << alphaEn_v2[n] << endl;
    }

    cout << "Analysis Ended" << endl;
    return 0;
}

bool sameSign(double num1, double num2)
{
    return (num1 >= 0 && num2 >= 0) || (num1 < 0 && num2 < 0);
}
