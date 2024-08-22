#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

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

    int tSkip = 1;
    double deltaT = 2.0;

    double mi = 2.66389e-26; // mass of one molecule
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

    ifstream nSteps("nTimeSteps.dat", ios::in);
    int nTimeSteps; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l
    nSteps >> nTimeSteps;

    // The Ids of tracked molecules
    vector<int> molIds;
    molIds.push_back(53);
    molIds.push_back(235);
    molIds.push_back(1000);
    molIds.push_back(5000);

    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    vector<vector<double>> xT(molIds.size());
    vector<vector<double>> yT(molIds.size());
    vector<vector<double>> zT(molIds.size());

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

            if (typ == 1) // methane
            {
                for (i = 0; i < molIds.size(); i++)
                {
                    if (id == molIds[i])
                    {
                        xT[i].push_back(x);
                        yT[i].push_back(y);
                        zT[i].push_back(z);
                    }
                }
            }
        }

        getline(data, line);
    }

    {
        // Cartesian
        cout << "Create ovito dump file" << endl;

        ofstream file("molecules.trj");

        for (j = 0; j < nTimeSteps; j++)
        {
            file << "ITEM: TIMESTEP" << endl;
            file << j << endl;
            file << "ITEM: NUMBER OF ATOMS" << endl;
            file << molIds.size() << endl;
            file << "ITEM: BOX BOUNDS pp pp pp" << endl;
            file << xLo << " " << xHi << endl;
            file << yLo << " " << yHi << endl;
            file << zLo << " " << zHi << endl;
            file << "ITEM: ATOMS id type x y z" << endl;

            for (i = 0; i < molIds.size(); i++)
            {
                file << i + 1 << " " << i << " " << xT[i][j] << " " << yT[i][j] << " " << zT[i][j] << endl;
            }
        }
    }
    
    {
        // Cartesian distance to the wall
        cout << "Distance to the wall" << endl;

        ofstream Distance("Distance.txt");

        for (j = 0; j < nTimeSteps; j++)
        {
            Distance << j;

            for (i = 0; i < molIds.size(); i++)
            {
                Distance << " " << zT[i][j];
            }

            Distance << endl;
        }

    }

    return 0;
}