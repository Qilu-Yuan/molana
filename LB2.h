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


class LB2
{
private:
        /* data */
    int NStat,NBond,NBin;
    double DelR, SumR2;
    vector<double> RBin,VBin,PR2;

public:
    void set(int *NBead,int *NPoly,fstream& Log);
    void calcluate(int* NBead,int* NPoly,const vector<Vector3D>& RX);
    void write(int *NBead);
};

