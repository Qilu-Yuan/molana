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

void SetSSF(double *Box,vector<int> &Qn,vector<int> &Qc,vector<vector<double>>& Qx,vector<vector<double>>& Qy,vector<vector<double>>& Qz,int *NQ,int Qmax,int QIMax,double Qrh);

void SSF(int Switch,int *NPart,double *Box,double *Vol,const vector<Vector3D>& RX);
