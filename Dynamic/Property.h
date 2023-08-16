#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include<fstream>
#include <chemfiles.hpp>
#include<iomanip>
#include<vector>

#include"MSD.h"
#include"FSQ.h"
#include"X4S.h"
#include"ASL.h"
#include"EEV.h"

using namespace std;
using namespace chemfiles;
using std::string;
using std::stringstream;

class Properties
{
private:
    int *TMax,*IT0,*T0Max;
    double *Delt;
    int *NPart,*NBead,*NPoly;
    double *Box,*Vol,*Rho;
    int Nt,Nt0;
    vector<int> NDyna, ttv0;
    
    double MDTime;
    vector<Vector3D> EX;
    vector<vector <double>> RX0, RY0, RZ0; 
    vector<vector<double>> EX0,EY0,EZ0;

    string name_msd,name_fsq,name_asl,name_x4s,name_eev;
    MSD MSD1;FSQ FSQ1;ASL ASL1;X4S X4S1; EEV EEV1;

public:
    Properties(int *NPart_in,int *NBead_in,int *NPoly_in,double *Box_in,double *Vol_in,double *Rho_in,int *TMax_in,int *IT0_in,int *T0Max_in,double *Delt_in,fstream& Log,char& type);
    void calculate(const vector<Vector3D>& RX,const vector<Vector3D>& PX);
    void write();
};
