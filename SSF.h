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

class SSF
{
private:
    /* data */
    const int Qmax = 10000,QIMax = 20;
    const double Qrh = 15.0;
    int NQ,NStat;
    vector<int> Qn,Qc;
    vector<vector<double>> Qx,Qy,Qz;
    vector<double> Q, Sbb;

public:
    void set(double *Box,fstream& Log);
    void calcluate(int *NPart,const vector<Vector3D>& RX);
    void write(int *NPart);
};

