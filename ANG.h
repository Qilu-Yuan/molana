#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <chemfiles.hpp>
#include<iomanip>
#include<fstream>
#include<vector>

using namespace std;
using namespace chemfiles;
using std::string;
using std::stringstream;

class ANG
{
private:
    /* data */
    int NStat,NAngle,NBin;
    double DelR, SumCOSTheta,SumTheta;
    vector<double> RBin,VBin,PTheta;
public:
    void set(int *NBead,int *NPoly,fstream& Log);
    void calcluate(int *NBead,int *NPoly,const vector<Vector3D>& RX);
    void write(int *NBead,int *NPoly);
};

