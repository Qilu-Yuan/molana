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

class ASL
{
private:
    /* data */
    int Nmb;
    vector<int> Pmb;
    double DelR2;
    vector<double> r2;
    
    int NCluster;
    vector<int> NIT, L, Num;
    vector<double> ST;
    void Get_Pmb(vector<double>& r2,vector<int>& Pmb, int *NPart);
public:
    void set(int *NPart,int *TMax);
    void calcluate(int *NPart,double *Box,int t,int dt,const vector<Vector3D>& RX,const vector<vector <double>> &RX0, const vector<vector <double>> &RY0, const vector<vector <double>> &RZ0);
    void write(int dt,const vector<int>& NDyna,double MDTime,fstream& file);
};

