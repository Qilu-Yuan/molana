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

class EEV
{
private:
    /* data */
    vector<double> SumEEV;
public:
    void set(int *TMax);
    void calcluate(int *NPoly,int t,int dt,const vector<Vector3D>& EX,const vector<vector <double>> &EX0, const vector<vector <double>> &EY0, const vector<vector <double>> &EZ0);
    void write(int *NPoly,int dt,const vector<int>& NDyna,double MDTime,fstream& file);
};


