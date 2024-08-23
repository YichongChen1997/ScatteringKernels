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
    string line;
    int id, n, i, j, k, t, nAtoms, currentTimeStep, typ, bo;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r, xLo, xHi, yLo, yHi, zLo, zHi, fx, fy, fz;

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

    // The timeIds of tracked molecules
    vector<int> timeIds;
    vector<int> nGasCount;

    for (int pick = 1; pick < 5; pick++)
    {
        timeIds.push_back(nTimeSteps * 2 / (pick + 2));
    }

    nGasCount.resize(timeIds.size(), 0);

    ifstream data("dump_meas_gas.lammpstrj", ios::in);
    ofstream dumpfile("nearWallMolecules.xyz", ios::out);

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
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3 >> KE >> PE >> fx >> fy >> fz;

            if (typ == 1)
            {
                for (i = 0; i < timeIds.size(); i++)
                {
                    if (t == timeIds[i])
                    {
                        if ((z < rCut) || (z > H - rCut))
                        {
                            nGasCount[i]++;
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
            totalCount += nGasCount[i];
        }
        dumpfile << " " << totalCount << endl;
        for (i = 0; i < timeIds.size(); i++)
        {
            dumpfile << " " << nGasCount[i];
        }
        dumpfile << endl;
        for (i = 0; i < timeIds.size(); i++)
        {
            for (j = 0; j < nGasCount[i]; j++)
            {
                dumpfile << i << " " << xP[i][j] << " " << yP[i][j] << " " << zP[i][j] << endl;
            }
        }
    }

    return 0;
}
