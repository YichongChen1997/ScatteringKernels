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
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

    double mole = 6.02214076e23;
    double KcalJ = 4184;
    double refEnergy = KcalJ / mole;

    //***************************
    // Specify the parameters for the post-processing
    int nTimeSteps;
    double rCut, Tw, H, Msigma, Mass;

    ifstream Specification("Specification.dat", ios::in);
    getline(Specification, line);
    Specification >> nTimeSteps  >> rCut >> Tw  >>  H  >>  Msigma >> Mass;
    // **************************

    double vM = sqrt(2 * kB * Tw / Mass); // most probable speed
    double NorV = sqrt(kB * Tw / Mass);   // reference velocity

    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    ofstream file1("AlphaTx.txt", ios::out);
    ofstream file2("AlphaN.txt",  ios::out);
    ofstream file3("AlphaEx.txt", ios::out);
    ofstream file4("AlphaEn.txt", ios::out);
    ofstream file5("AlphaE.txt",  ios::out);

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

    vector<double> range = linspace(0.1, 2.9, 29);

    for (int n = 0; n < range.size(); n++)
    {
        double alphaTx_f_v1, alphaTx_f_v2, alphaTx_p_v1, alphaTx_p_v2;
        double alphaN_f_v1,  alphaN_f_v2,  alphaN_p_v1,  alphaN_p_v2;
        double alphaEx_f_v1, alphaEx_f_v2, alphaEx_p_v1, alphaEx_p_v2;
        double alphaEn_f_v1, alphaEn_f_v2, alphaEn_p_v1, alphaEn_p_v2;
        double alphaE_f_v1,  alphaE_f_v2,  alphaE_p_v1,  alphaE_p_v2;

        // AlphaTx
        double vTx_f_all = 0.0, vTx_f_mean, vTxi_f_all = 0.0, vTxi_f_mean;
        double vTx_p_all = 0.0, vTx_p_mean, vTxi_p_all = 0.0, vTxi_p_mean;

        // AlphaN
        double vN_f_all = 0.0, vN_f_mean, vNi_f_all = 0.0, vNi_f_mean;
        double vN_p_all = 0.0, vN_p_mean, vNi_p_all = 0.0, vNi_p_mean;

        // AlphaEx
        double Ex_f_all = 0.0, Ex_f_mean, Exi_f_all = 0.0, Exi_f_mean;
        double Ex_p_all = 0.0, Ex_p_mean, Exi_p_all = 0.0, Exi_p_mean;

        // AlphaEn
        double En_f_all = 0.0, En_f_mean, Eni_f_all = 0.0, Eni_f_mean;
        double En_p_all = 0.0, En_p_mean, Eni_p_all = 0.0, Eni_p_mean;

        // AlphaE
        double E_f_all = 0.0, E_f_mean, Ei_f_all = 0.0, Ei_f_mean;
        double E_p_all = 0.0, E_p_mean, Ei_p_all = 0.0, Ei_p_mean;

        // count the No. of valid data
        int cF_Tx = 0, cP_Tx = 0, cF_N = 0, cP_N = 0;
        int cF_Ex = 0, cP_Ex = 0, cF_En = 0, cP_En = 0, cF_E = 0, cP_E = 0;


        for (int i = 0; i < nPts; ++i)
        {
            // summation of momentum (global)
            if (vTxi[i] >= 0.0)
            {
                vTxi_f_all += vTxi[i];
                vTx_f_all  += vTx[i];
                cF_Tx++;
            }

            if ((vNi[i] >= vM * 0.0) && (vNi[i] < vM * 6))
            {
                vNi_f_all += vNi[i];
                vN_f_all += vN[i];
                cF_N++;
            }

            // summation of energy (global)
            Exi_f_all += keTxi[i];
            Ex_f_all  += keTx[i];

            Eni_f_all += keNi[i];
            En_f_all  += keN[i];

            Ei_f_all += kei[i];
            E_f_all  += ke[i];

            cF_Tx++;
            cF_N++;
            cF_Ex++;
            cF_En++;
            cF_E++;

            // summation of momentum (partial)
            if ((vTxi[i] > vM * (range[n] - 0.1)) && (vTxi[i] <= vM * (range[n] + 0.1)))
            {
                vTxi_p_all += vTxi[i];
                vTx_p_all  += vTx[i];
                cP_Tx++;

                Exi_p_all += keTxi[i];
                Ex_p_all  += keTx[i];
                cP_Ex++;
            }

            if ((vNi[i] > vM * (range[n] - 0.1)) && (vNi[i] <= vM * (range[n] + 0.1)))
            {
                vNi_p_all += vNi[i];
                vN_p_all  += vN[i];
                cP_N++;

                Eni_p_all += keNi[i];
                En_p_all  += keN[i];
                cP_En++;
            }

            // summation of energy (partial)
            if ((kei[i] > NorV * NorV * (range[n] - 0.1)) && (kei[i] <= NorV * NorV * (range[n] + 0.1)))
            {
                Ei_p_all += kei[i];
                E_p_all  += ke[i];
                cP_E++;
            }            
        }
        // *********************************
        vTx_f_mean = vTx_f_all / cF_Tx;
        vTx_p_mean = vTx_p_all / cP_Tx;

        vTxi_f_mean = vTxi_f_all / cF_Tx;
        vTxi_p_mean = vTxi_p_all / cP_Tx;
        // *********************************
        vN_f_mean  = vN_f_all  / cF_N;
        vN_p_mean  = vN_p_all  / cP_N;

        vNi_f_mean  = vNi_f_all  / cF_N;
        vNi_p_mean  = vNi_p_all  / cP_N;
        // *********************************
        Ex_f_mean  = Ex_f_all  / cF_Ex;
        Ex_p_mean  = Ex_p_all  / cP_Ex;

        Exi_f_mean  = Exi_f_all  / cF_Ex;
        Exi_p_mean  = Exi_p_all  / cP_Ex;
        // *********************************
        En_f_mean  = En_f_all  / cF_En;
        En_p_mean  = En_p_all  / cP_En;

        Eni_f_mean  = Eni_f_all  / cF_En;
        Eni_p_mean  = Eni_p_all  / cP_En;
        // *********************************
        E_f_mean   = E_f_all   / cF_E;
        E_p_mean   = E_p_all   / cP_E;

        Ei_f_mean  = Ei_f_all   / cF_E;
        Ei_p_mean  = Ei_p_all   / cP_E;
        // *********************************

        // ********Tangential Momentum Accommodation Coefficients***********
        // *****************************************************************
        // Partial TMAC
        /* Method One (General expression) */
        alphaTx_p_v1 = (vTxi_p_all - vTx_p_all) / vTxi_p_all;

        /* Method Two (Least-Squares fitting) */
        double beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            if ((vTxi[i] > vM * (range[n] - 0.1)) && (vTxi[i] <= vM * (range[n] + 0.1)))
            {
                beta_nu += (vTxi[i] - vTxi_p_mean) * (vTx[i]  - vTx_p_mean);     // beta numerator
                beta_de += (vTxi[i] - vTxi_p_mean) * (vTxi[i] - vTxi_p_mean);    // denominator
            }
        }
        alphaTx_p_v2 = 1 - beta_nu / beta_de;

        // Global TMAC
        /* Method One (General expression) */
        alphaTx_f_v1 = (vTxi_f_all - vTx_f_all) / vTxi_f_all;

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            beta_nu += (vTxi[i] - vTxi_f_mean) * (vTx[i]  - vTx_f_mean);   // beta numerator
            beta_de += (vTxi[i] - vTxi_f_mean) * (vTxi[i] - vTxi_f_mean);  // beta donominator
        }
        alphaTx_f_v2 = 1 - beta_nu / beta_de;
        // *******************************************************************
        // *******************************************************************

        // ******** Normal Momentum Accommodation Coefficients ***************
        // *******************************************************************
        // Partial NMAC
        /* Method One (General expression) */
        alphaN_p_v1 = (vNi_p_all - vN_p_all) / (vNi_p_all - sqrt(PI) / 2 * vM * cP_N);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            if ((vNi[i] > vM * (range[n] - 0.1)) && (vNi[i] <= vM * (range[n] + 0.1)))
            {
                beta_nu += (vNi[i] - vNi_p_mean) * (vN[i]  - vN_p_mean);   // beta numerator
                beta_de += (vNi[i] - vNi_p_mean) * (vNi[i] - vNi_p_mean);   // denominator
            }
        }
        alphaN_p_v2 = 1 - beta_nu / beta_de;

        // Global NMAC
        /* Method One (General expression) */
        alphaN_f_v1 = (vNi_f_all - vN_f_all) / (vNi_f_all - sqrt(PI) / 2 * vM * cF_N);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            beta_nu += (vNi[i] - sqrt(PI) / 2 * vM) * (vN[i]  - vN_f_mean);   // beta numerator
            beta_de += (vNi[i] - sqrt(PI) / 2 * vM) * (vNi[i] - sqrt(PI) / 2 * vM);  // beta donominator
        }
        alphaN_f_v2 = 1 - beta_nu / beta_de;
        // *******************************************************************
        // *******************************************************************

        // ******** Tangential Kinetic Energy Acccommodation Coefficients **********
        // *************************************************************************
        // Partial TKEAC
        /* Method One (General expression) */
        alphaEx_p_v1 = (Exi_p_all - Ex_p_all) / (Exi_p_all - vM * vM / 2 * cP_Ex);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            if ((vTxi[i] > vM * (range[n] - 0.1)) && (vTxi[i] <= vM * (range[n] + 0.1)))
            {
                beta_nu += (keTxi[i] - Exi_p_mean) * (keTx[i]  - Ex_p_mean);    // beta numerator
                beta_de += (keTxi[i] - Exi_p_mean) * (keTxi[i] - Exi_p_mean);   // denominator
            }
        }
        alphaEx_p_v2 = 1 - beta_nu / beta_de;

        // Global TKEAC
        /* Method One (General expression) */
        alphaEx_f_v1 = (Exi_f_all - Ex_f_all) / (Exi_f_all - vM * vM / 2 * cF_Ex);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            beta_nu += (keTxi[i] - Exi_f_mean) * (keTx[i]  - Ex_f_mean);    // beta numerator
            beta_de += (keTxi[i] - Exi_f_mean) * (keTxi[i] - Exi_f_mean);   // beta donominator
        }
        alphaEx_f_v2 = 1 - beta_nu / beta_de;
        // *************************************************************************
        // *************************************************************************

        // ******** Normal Kinetic Energy Acccommodation Coefficients***************
        // *************************************************************************
        // Partial NKEAC
        /* Method One (General expression) */
        alphaEn_p_v1 = (Eni_p_all - En_p_all) / (Eni_p_all - vM * vM * cP_En);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            if ((vNi[i] > vM * (range[n] - 0.1)) && (vNi[i] <= vM * (range[n] + 0.1)))
            {
                beta_nu += (keNi[i] - Eni_p_mean) * (keN[i]  - En_p_mean);    // beta numerator
                beta_de += (keNi[i] - Eni_p_mean) * (keNi[i] - Eni_p_mean);   // denominator
            }
        }
        alphaEn_p_v2 = 1 - beta_nu / beta_de;

        // Global NKEAC
        /* Method One (General expression) */
        alphaEn_f_v1 = (Eni_f_all - En_f_all) / (Eni_f_all - vM * vM * cF_En);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            beta_nu += (keNi[i] - Eni_f_mean) * (keN[i]  - En_f_mean);    // beta numerator
            beta_de += (keNi[i] - Eni_f_mean) * (keNi[i] - Eni_f_mean);   // beta donominator
        }
        alphaEn_f_v2 = 1 - beta_nu / beta_de;
        // *************************************************************************

        // ***************  Energy Acccommodation Coefficients *********************
        // *************************************************************************
        // Partial EAC
        /* Method One (General expression) */
        alphaE_p_v1 = (Ei_p_all - E_p_all) / (Ei_p_all - 2 * vM * vM * cP_E);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            if ((kei[i] > NorV * NorV * (range[n] - 0.1)) && (kei[i] <= NorV * NorV * (range[n] + 0.1)))
            {
                beta_nu += (kei[i] - Ei_p_mean) * (ke[i]  - E_p_mean);    // beta numerator
                beta_de += (kei[i] - Ei_p_mean) * (kei[i] - Ei_p_mean);   // denominator
            }
        }
        alphaE_p_v2 = 1 - beta_nu / beta_de;

        // Global EAC
        /* Method One (General expression) */
        alphaE_f_v1 = (Ei_f_all - E_f_all) / (Ei_f_all - 2 * vM * vM * cF_E);

        /* Method Two (Least-Squares fitting) */
        beta_nu = 0.0, beta_de = 0.0;

        for (int i = 0; i < nPts; ++i)
        {
            beta_nu += (kei[i] - Ei_f_mean) * (ke[i]  - E_f_mean);    // beta numerator
            beta_de += (kei[i] - Ei_f_mean) * (kei[i] - Ei_f_mean);   // beta donominator
        }
        alphaE_f_v2 = 1 - beta_nu / beta_de;
        // *************************************************************************
        // *************************************************************************

        file1 << range[n] << " " << alphaTx_p_v1 << " " << alphaTx_p_v2 << " " << alphaTx_f_v1 << " " << alphaTx_f_v2 << endl;
        file2 << range[n] << " " << alphaN_p_v1  << " " << alphaN_p_v2  << " " << alphaN_f_v1  << " " << alphaN_f_v2  << endl;
        file3 << range[n] << " " << alphaEx_p_v1 << " " << alphaEx_p_v2 << " " << alphaEx_f_v1 << " " << alphaEx_f_v2 << endl;
        file4 << range[n] << " " << alphaEn_p_v1 << " " << alphaEn_p_v2 << " " << alphaEn_f_v1 << " " << alphaEn_f_v2 << endl;
        file5 << range[n] << " " << alphaE_p_v1  << " " << alphaE_p_v2  << " " << alphaE_f_v1  << " " << alphaE_f_v2  << endl;
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