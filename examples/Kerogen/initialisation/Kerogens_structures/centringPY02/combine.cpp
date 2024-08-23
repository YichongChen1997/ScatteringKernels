#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;
const float PI=3.1415;

int main()
{
	string line, typeName;
	char l;
	double x,y,z,chg,xLo,xHi,yLo,yHi,zLo,zHi,xT;
	int atoms,typ,i,j,k,id,mol=0,molId=0,nWall,bID,aID,nCNT=0,nWallCount=0,nWallK=0,nKWallCount=0,co;
	nWall=0;
    int nAtomsTot =0;
    double shiftX = 0, shiftY = 0, shiftZ = 0;
	double mi = 2.66389E-26;    

	double refLength = 1e-10;
  	double refForce = 6.9477E-11;

    ofstream writer("PY02_50A_1.20_CENTRED.xyz");    
    ifstream readerKerogen("PY02_50A_1.20_NOTCENTRED.xyz"); // kerogen
    readerKerogen >> nWallK;
    writer<< nWallK << endl;
    writer<< endl;
    
    getline(readerKerogen,line);
    getline(readerKerogen,line);

	vector<double> xK;
	vector<double> yK;
	vector<double> zK;
	vector<int> idsK;
	
	double xMin = 10000;
	
	for(i=0;i<nWallK;i++)
    {
        readerKerogen >> typeName >> x >> y >> z;
        
        shiftX = -25;
        shiftY = -25;//<-
        shiftZ = -25;
              
		x += shiftX;
		y += shiftY;
		z += shiftZ;
		
		writer<< typeName << " " << x << " " << y << " " << z <<endl;	
    }

    writer.close();

    return 0;

}
