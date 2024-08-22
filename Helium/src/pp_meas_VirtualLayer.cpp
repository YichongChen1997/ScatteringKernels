#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

// Output the coordinate of every gas molecules within the virtual plane (near-surface layer), at any specified time.
using namespace std;

int main()
{
    string line, fName;
    int id, n, i, j, k, t, tTime, nAtoms, yStep, currentTimeStep, typ, bo;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r, xLo, xHi, yLo, yHi, zLo, zHi;

    double TAU1Left, TAU2Left, TAU3Left, TAU1Right, TAU2Right, TAU3Right, keRight;
    double pistonLeftX, pistonRightX, barrierLeftX, barrierRightX;
    int pistonLeftCount, pistonRightCount, barrierLeftCount, barrierRightCount;
    int nMethaneLeft, nMethaneRight, nMethaneMiddle;
    double heightLeft, heightRight, volumeLeft, volumeRight, rhoLeft, rhoRight, pLeft, pRight, tempRight;
    double TAU1Middle, TAU2Middle, TAU3Middle, keMiddle, pMiddle, tempMiddle;
    double xMin, xMax;

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

    int tSkip = 1;
    double deltaT = 1.0;

    // The timeIds of tracked molecules
    vector<int> timeIds;
    vector<int> nMethaneCount;

    for (int pick = 1; pick < 5; pick++)
    {
        timeIds.push_back(nTimeSteps * 2 / (pick + 2));
    }

    nMethaneCount.resize(timeIds.size(), 0);

    ifstream data("dump_meas_gas.lammpstrj", ios::in);
    ofstream dumpfile("Molecules.xyz", ios::out);

    vector<vector<double> > xP(timeIds.size());
    vector<vector<double> > yP(timeIds.size());
    vector<vector<double> > zP(timeIds.size());

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

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3;

            if (typ == 1)
            {
                for (i = 0; i < timeIds.size(); i++)
                {
                    if (t == timeIds[i])
                    {
                        if ((z < 15) || (z > 185))
                        {
                            nMethaneCount[i]++;
                            xP[i].push_back(x);
                            yP[i].push_back(y);
                            zP[i].push_back(z);
                        }
                    }
                }
            }
        }

        getline(data, line);
    }

    {
        // Cartesian coordinate of the molecules
        cout << "Create xyz dump file" << endl;

        int totalCount = 0;
        for (i = 0; i < timeIds.size(); i++)
        {
            totalCount += nMethaneCount[i];
        }
        dumpfile << " " << totalCount << endl;
        for (i = 0; i < timeIds.size(); i++)
        {
            dumpfile << " " << nMethaneCount[i];
        }
        dumpfile << endl;
        for (i = 0; i < timeIds.size(); i++)
        {
            for (j = 0; j < nMethaneCount[i]; j++)
            {
                dumpfile << i << " " << xP[i][j] << " " << yP[i][j] << " " << zP[i][j] << endl;
            }
        }
    }

    return 0;
}
