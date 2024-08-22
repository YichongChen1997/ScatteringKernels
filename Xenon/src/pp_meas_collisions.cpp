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

    double mi = 2.1802E-25; // mass of one molecule
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

    //**** CHANGE ***************
    double Tw = 423;                    // temperature of wall (Kelvin)
    double vM = sqrt(2 * kB * Tw / mi); // most probable speed

    ifstream Height("Height.dat", ios::in);
    double H;
    Height >> H;
    double rCut = 15;

    ifstream nSteps("nTimeSteps.dat", ios::in);
    int nTimeSteps; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l
    nSteps >> nTimeSteps;
    // **************************

    ifstream data("dump_meas_gas.lammpstrj", ios::in);
    ofstream fileAtom("AtomCollisions.txt", ios::out);
    ofstream collisionTMAC("collision_vs_TMACs.txt", ios::out);
    ofstream collisionNKEAC("collision_vs_NKEAC.txt", ios::out);

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

    vector<double> vNS, vNiS; // velocity - NormalDiretion - SingleCollisions
    vector<double> vTS, vTiS;
    vector<double> vTxS, vTxiS;
    vector<double> vTyS, vTyiS;
    vector<double> vMagS, vMagiS;
    double vMagTS, vMagTiS;
    int S_coTop = 0, S_coBottom = 0;

    vector<double> vNM, vNiM; // velocity - NormalDiretion - MultipleCollisions
    vector<double> vTM, vTiM;
    vector<double> vTxM, vTxiM;
    vector<double> vTyM, vTyiM;
    vector<double> vMagM, vMagiM;
    double vMagTM, vMagTiM;
    int M_coTop = 0, M_coBottom = 0;

    for (i = 0; i < nWallCollisions; i++)
    {
        if (leftFromTop[i] == 1)
        {
            if (AtomCollisions[i] <= 3) // single collision cases
            // if (AtomCollisions[i] == 1)
            {
                // normal
                vNS.push_back(abs(vzEnd[i]));
                vNiS.push_back(abs(vzStart[i]));

                // combined tangents
                vMagTS = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i];
                vMagTiS = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i];

                // total velocity
                vMagS.push_back(sqrt(abs(vzEnd[i]) * abs(vzEnd[i]) + vMagTS * vMagTS));
                vMagiS.push_back(sqrt(abs(vzStart[i]) * abs(vzStart[i]) + vMagTiS * vMagTiS));

                if (vMagTS > 0)
                {
                    vTS.push_back(sqrt(vMagTS));
                }
                else
                {
                    vTS.push_back(0.0);
                }

                if (vMagTiS > 0)
                {
                    vTiS.push_back(sqrt(vMagTiS));
                }
                else
                {
                    vTiS.push_back(0.0);
                }

                vTxS.push_back(vxEnd[i]);
                vTxiS.push_back(vxStart[i]);

                vTyS.push_back(vyEnd[i]);
                vTyiS.push_back(vyStart[i]);

                S_coTop++;
            }

            else
            {
                // normal
                vNM.push_back(abs(vzEnd[i]));
                vNiM.push_back(abs(vzStart[i]));

                // combined tangents
                vMagTM = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i];
                vMagTiM = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i];

                // total velocity
                vMagM.push_back(sqrt(abs(vzEnd[i]) * abs(vzEnd[i]) + vMagTM * vMagTM));
                vMagiM.push_back(sqrt(abs(vzStart[i]) * abs(vzStart[i]) + vMagTiM * vMagTiM));

                if (vMagTM > 0)
                {
                    vTM.push_back(sqrt(vMagTM));
                }
                else
                {
                    vTM.push_back(0.0);
                }

                if (vMagTiM > 0)
                {
                    vTiM.push_back(sqrt(vMagTiM));
                }
                else
                {
                    vTiM.push_back(0.0);
                }

                vTxM.push_back(vxEnd[i]);
                vTxiM.push_back(vxStart[i]);

                vTyM.push_back(vyEnd[i]);
                vTyiM.push_back(vyStart[i]);

                M_coTop++;
            }
        }

        if (leftFromTop[i] == -1)
        {

            if (AtomCollisions[i] <= 3) // single collision cases
            // if (AtomCollisions[i] == 1)
            {
                // normal
                vNS.push_back(abs(vzEnd[i]));
                vNiS.push_back(abs(vzStart[i]));

                // combined tangents
                vMagTS = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i];
                vMagTiS = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i];

                // total velocity
                vMagS.push_back(sqrt(abs(vzEnd[i]) * abs(vzEnd[i]) + vMagTS * vMagTS));
                vMagiS.push_back(sqrt(abs(vzStart[i]) * abs(vzStart[i]) + vMagTiS * vMagTiS));

                if (vMagTS > 0)
                {
                    vTS.push_back(sqrt(vMagTS));
                }
                else
                {
                    vTS.push_back(0.0);
                }

                if (vMagTiS > 0)
                {
                    vTiS.push_back(sqrt(vMagTiS));
                }
                else
                {
                    vTiS.push_back(0.0);
                }

                vTxS.push_back(vxEnd[i]);
                vTxiS.push_back(vxStart[i]);

                vTyS.push_back(vyEnd[i]);
                vTyiS.push_back(vyStart[i]);

                S_coBottom++;
            }

            else
            {
                // normal
                vNM.push_back(abs(vzEnd[i]));
                vNiM.push_back(abs(vzStart[i]));

                // combined tangents
                vMagTM = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i];
                vMagTiM = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i];

                // total velocity
                vMagM.push_back(sqrt(abs(vzEnd[i]) * abs(vzEnd[i]) + vMagTM * vMagTM));
                vMagiM.push_back(sqrt(abs(vzStart[i]) * abs(vzStart[i]) + vMagTiM * vMagTiM));

                if (vMagTM > 0)
                {
                    vTM.push_back(sqrt(vMagTM));
                }
                else
                {
                    vTM.push_back(0.0);
                }

                if (vMagTiM > 0)
                {
                    vTiM.push_back(sqrt(vMagTiM));
                }
                else
                {
                    vTiM.push_back(0.0);
                }

                vTxM.push_back(vxEnd[i]);
                vTxiM.push_back(vxStart[i]);

                vTyM.push_back(vyEnd[i]);
                vTyiM.push_back(vyStart[i]);

                M_coBottom++;
            }
        }
    }

    cout << "No. of Single Collisions: " << vNS.size() << endl;
    cout << "No. of Single Collisions at top: " << S_coTop << endl;
    cout << "No. of Single Collisions at bottom: " << S_coBottom << endl;
    cout << endl;

    cout << "No. of Multiple Collisions: " << vNM.size() << endl;
    cout << "No. of Multiple Collisions at top: " << M_coTop << endl;
    cout << "No. of Multiple Collisions at bottom: " << M_coBottom << endl;
    cout << endl;

    cout << "No. of Total Collisions: " << AtomCollisions.size() << endl;

    int nPts = AtomCollisions.size();
    int SnPts = vNS.size();
    int MnPts = vNM.size();

    double CoBinWidth = 1; // binwidth of collision distribution
    int CoMax = 10;        // Max times of collisions
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
            fileAtom << AtomCollisions[i] << endl;
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
        ofstream file("CollisionDis.txt");

        for (i = 0; i < CoBins; i++)
        {
            file << bins[i] << " " << distribution[i] << " " << probability[i] << endl;
        }
    }
    cout << "Collision Distribution Finished" << endl;

    int collisionPts = 20;
    vector<double> TMAC, TMAC2;
    vector<double> vTx_meanF;
    vector<double> vTxi_meanF;
    TMAC.resize(collisionPts, 0.0);
    TMAC2.resize(collisionPts, 0.0);
    vTx_meanF.resize(collisionPts, 0.0);
    vTxi_meanF.resize(collisionPts, 0.0);

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
                    vTxi_meanF[n] += vxStart[i];
                    vTx_meanF[n] += vxEnd[i];
                    countPositiveTx++;
                }
            }
        }
        vTxi_meanF[n] = vTxi_meanF[n] / countPositiveTx;
        vTx_meanF[n] = vTx_meanF[n] / countPositiveTx;

        for (int i = 0; i < nPts; i++)
        {
            if (AtomCollisions[i] == n)
            {
                beta_nu += (vxStart[i] - vTxi_meanF[n]) * (vxEnd[i] - vTx_meanF[n]);
                beta_de += (vxStart[i] - vTxi_meanF[n]) * (vxStart[i] - vTxi_meanF[n]);
            }
        }
        TMAC[n] = 1 - beta_nu / beta_de;
        
        TMAC2[n] = (vTxi_meanF[n] - vTx_meanF[n]) / vTxi_meanF[n];

        collisionTMAC << n << " " << TMAC[n] << " " << TMAC2[n] << endl;
    }

    vector<double> alphaN1, alphaN2;
    vector<double> vN_meanF, vNi_meanF;

    alphaN1.resize(collisionPts, 0.0);
    alphaN2.resize(collisionPts, 0.0);
    vN_meanF.resize(collisionPts, 0.0);
    vNi_meanF.resize(collisionPts, 0.0);

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
                    vNi_meanF[n] += abs(vzStart[i]);
                    vN_meanF[n]  += abs(vzEnd[i]);
                    countPositiveN++;
                }
            }
        }

        vNi_meanF[n] = vNi_meanF[n] / countPositiveN;
        vN_meanF[n]  = vN_meanF[n] / countPositiveN;

        alphaN1[n] = (vNi_meanF[n] * vNi_meanF[n] - vN_meanF[n] * vN_meanF[n]) / (vNi_meanF[n] * vNi_meanF[n] - (sqrt(PI) * vM / 2) * (sqrt(PI) * vM / 2));

        for (int i = 0; i < nPts; i++)
        {
            if (AtomCollisions[i] == n)
            {
                beta_nu += (vzStart[i] * vzStart[i] - vNi_meanF[n] * vNi_meanF[n]) * (vzEnd[i] * vzEnd[i] - vN_meanF[n] * vN_meanF[n]);
                beta_de += (vzStart[i] * vzStart[i] - vNi_meanF[n] * vNi_meanF[n]) * (vzStart[i] * vzStart[i] - vNi_meanF[n] * vNi_meanF[n]);
            }
        }

        alphaN2[n] = 1 - beta_nu / beta_de;

        collisionNKEAC << n << " " << alphaN1[n] << " " << alphaN2[n] << endl;
    }

    cout << "Analysis Ended" << endl;
    return 0;
}

bool sameSign(double num1, double num2)
{
    return (num1 >= 0 && num2 >= 0) || (num1 < 0 && num2 < 0);
}
