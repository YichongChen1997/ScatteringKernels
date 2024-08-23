#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

using namespace std;

int main()
{
    string line;
    int id, n, i, j, k, t, nAtoms, currentTimeStep, typ, bo;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, xLo, xHi, yLo, yHi, zLo, zHi, fx, fy, fz;

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

    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    //***************************
    // The Ids of tracked molecules
    vector<int> molIds;
    
    molIds.push_back(53);
    molIds.push_back(235);
    molIds.push_back(1000);
    molIds.push_back(5000);

    vector<vector<double> > xT(molIds.size());
    vector<vector<double> > yT(molIds.size());
    vector<vector<double> > zT(molIds.size());

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

        ofstream Distance("distanceToSurf.txt");
        
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