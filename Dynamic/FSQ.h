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

class FSQ
{
private:
    /* data */
    vector<vector<double>> Qx;
    vector<double> FCOS1, FCOS2;
    int Qn;
    const int QImax = 20;
    void SetFSQ(double *Box,int *Qn,vector<vector<double>>& Qx,const int *QImax,double *Q,double *QrSet);
public:
    void set(double *Box,int *TMax,fstream& Log);
    void calcluate(int *NPart,int t,int dt,const vector<Vector3D>& RX,const vector<vector <double>> &RX0, const vector<vector <double>> &RY0, const vector<vector <double>> &RZ0);
    void write(int *NPart,int dt,const vector<int>& NDyna,double MDTime,fstream& file);
};  

