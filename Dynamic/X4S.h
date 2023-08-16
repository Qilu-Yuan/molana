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

class X4S
{
private:
    /* data */
    const double a2=0.09;
    vector<double> Qs1, Qs2;
public:
    void set(int *TMax);
    void calcluate(int *NPart,double *Box,int t,int dt,const vector<Vector3D>& RX,const vector<vector <double>> &RX0, const vector<vector <double>> &RY0, const vector<vector <double>> &RZ0);
    void write(int *NPart,double *Vol,int dt,const vector<int>& NDyna,double MDTime,fstream& file);

};


