#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <limits> // Add this line to include numeric_limits
#include <ios>    // Add this line to include streamsize

using namespace std;

const float PI = 3.1415926;

// Function to read parameters from the specification file
void readParameters(const string& filename, int& nTimeSteps, double& deltaT, double& tSkip, double& skipTimeStep,
                    double& rCut, double& H, double& Tg, double& Tw, double& mG, double& mW, double& mGas, double& mWall) {
    ifstream specification(filename);
    if (!specification.is_open()) {
        cerr << "Error opening " << filename << " file" << endl;
        exit(1);
    }
    
    specification >> nTimeSteps;
    specification.ignore(numeric_limits<streamsize>::max(), '\n');
    specification >> deltaT >> tSkip >> skipTimeStep;
    specification.ignore(numeric_limits<streamsize>::max(), '\n');
    specification >> rCut >> H;
    specification.ignore(numeric_limits<streamsize>::max(), '\n');
    specification >> Tg >> Tw;
    specification.ignore(numeric_limits<streamsize>::max(), '\n');
    specification >> mG >> mW;
    specification.ignore(numeric_limits<streamsize>::max(), '\n');
    
    double refMass = 1.66054e-27;
    mGas = mG * refMass;
    mWall = mW * refMass;
}

// Function to initialize 3D vectors for bins
template <typename T>
vector<vector<vector<T> > > initialize3DVector(int xSize, int ySize, int zSize, T initialValue) {
    return vector<vector<vector<T> > >(zSize, vector<vector<T> >(ySize, vector<T>(xSize, initialValue)));
}

int main() {
    string line;
    int id, type, n, t, nAtoms, currentTimeStep;
    double x, y, z, vx, vy, vz, tau1, tau2, tau3, PE, KE, fx, fy, fz;
    double xLo, xHi, yLo, yHi, zLo, zHi;
    double kB = 1.38064852e-23;

    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double mole = 6.02214076e23;
    double KcalJ = 4184;
    double refEnergyK = KcalJ / mole / kB;
    double refForceN = (KcalJ / mole) / refLength;

    // Read parameters from the specification file
    int nTimeSteps;
    double deltaT, tSkip, skipTimeStep, rCut, H, Tg, Tw, mG, mW, mGas, mWall;
    readParameters("Specification.dat", nTimeSteps, deltaT, tSkip, skipTimeStep, rCut, H, Tg, Tw, mG, mW, mGas, mWall);

    int NxBins = 800, NyBins = 800, NzBins = static_cast<int>(20 * rCut); // Specify the number of bins for potential contours
    double mostProbableSpeed = sqrt(2 * kB * Tw / mGas); // Most probable speed
    double refVelocityNorm = sqrt(kB * Tw / mGas); // Reference velocity

    // Open files for writing and reading
    ofstream filePE("PES.txt");
    ofstream fileFZ("Fz.txt");
    ifstream dataFile("dump_meas_gas.lammpstrj");

    // Initialize vectors to keep track of molecule properties
    vector<double> timeVec, posXVec, posYVec, posZVec, velXVec, velYVec, velZVec, keVec, peVec, forceZVec;

    // Initialize 3D vectors for bins
    auto PES_count_Bins = initialize3DVector(NxBins, NyBins, NzBins, 0);
    auto PES_value_Bins = initialize3DVector(NxBins, NyBins, NzBins, 0.0);
    auto FZ_count_Bins = initialize3DVector(NxBins, NyBins, NzBins, 0);
    auto FZ_value_Bins = initialize3DVector(NxBins, NyBins, NzBins, 0.0);

    for (t = 0; t < nTimeSteps; ++t) {
        for (n = 1; n <= 9; ++n) {
            if (n == 4) dataFile >> nAtoms;
            if (n == 2) {
                dataFile >> currentTimeStep;
                cout << "currentTimeStep = " << currentTimeStep
                     << "; t = " << t << " [ "
                     << 100.0 * (t + 1) / nTimeSteps
                     << "% ]" << endl;
            }
            if (n == 6) dataFile >> xLo >> xHi;
            if (n == 7) dataFile >> yLo >> yHi;
            if (n == 8) dataFile >> zLo >> zHi;
            getline(dataFile, line);
        }

        // Initialize vectors on the first time step
        if (t == 0) {
            timeVec.resize(nAtoms, 0);
            posXVec.resize(nAtoms, 0.0);
            posYVec.resize(nAtoms, 0.0);
            posZVec.resize(nAtoms, 0.0);
            velXVec.resize(nAtoms, 0.0);
            velYVec.resize(nAtoms, 0.0);
            velZVec.resize(nAtoms, 0.0);
            keVec.resize(nAtoms, 0.0);
            peVec.resize(nAtoms, 0.0);
            forceZVec.resize(nAtoms, 0.0);
        }

        for (n = 0; n < nAtoms; ++n) {
            dataFile >> id >> type >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3 >> KE >> PE >> fx >> fy >> fz;

            if (type == 1) { // Gases
                timeVec[n] = t;
                posXVec[n] = x;
                posYVec[n] = y;
                posZVec[n] = z;
                velXVec[n] = vx * refVelocity;
                velYVec[n] = vy * refVelocity;
                velZVec[n] = vz * refVelocity;
                keVec[n] = KE * refEnergyK;
                peVec[n] = PE * refEnergyK;
                forceZVec[n] = fz * refForceN;

                if (z < rCut) {
                    int idx = static_cast<int>(floor(posXVec[n] / (xHi - xLo) * NxBins));
                    int idy = static_cast<int>(floor(posYVec[n] / (yHi - yLo) * NyBins));
                    int idz = static_cast<int>(floor((posZVec[n] / rCut) * NzBins + 1));

                    idx = max(0, min(NxBins - 1, idx));
                    idy = max(0, min(NyBins - 1, idy));
                    idz = max(0, min(NzBins - 1, idz));

                    // Accumulate samples
                    PES_value_Bins[idz][idy][idx] += peVec[n];
                    PES_count_Bins[idz][idy][idx]++;
                    FZ_value_Bins[idz][idy][idx] += forceZVec[n];
                    FZ_count_Bins[idz][idy][idx]++;
                }
            }
        }

        getline(dataFile, line);
    }

    // Write results to files
    for (size_t i = 0; i < PES_value_Bins.size(); ++i) {
        for (size_t j = 0; j < PES_value_Bins[i].size(); ++j) {
            for (size_t k = 0; k < PES_value_Bins[i][j].size(); ++k) {
                if (PES_count_Bins[i][j][k] == 0) {
                    PES_value_Bins[i][j][k] = 0.0;
                    FZ_value_Bins[i][j][k] = 0.0;
                } else {
                    PES_value_Bins[i][j][k] /= PES_count_Bins[i][j][k];
                    FZ_value_Bins[i][j][k] /= FZ_count_Bins[i][j][k];
                }
                filePE << PES_value_Bins[i][j][k] << " ";
                fileFZ << FZ_value_Bins[i][j][k] << " ";
            }
            filePE << endl;
            fileFZ << endl;
        }
        filePE << endl;
        fileFZ << endl;
    }
    return 0;
}