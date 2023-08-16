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

using namespace std;
using namespace chemfiles;
using std::string;
using std::stringstream;

class BasicInfo
{   
public:
    BasicInfo(int *NPart,int *NBead,int *NPoly,double *Box,double *Vol,double *Rho,double *DeltL,double *TProdL);
    void cpu_time(string& Time);
    void Init(int *NFrame, int *NAtom,double *Delt_in,double *TProd,double *Delt,int *IT0,const int *T0Max,string &Trajname);
};