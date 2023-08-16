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

#include"Init.h"

using namespace std;
using namespace chemfiles;
using std::string;
using std::stringstream;
// using chemfiles::Trajectory;

BasicInfo::BasicInfo(int *NPart,int *NBead,int *NPoly,double *Box,double *Vol,double *Rho,double *Delt_in,double *TProd_in){

    string line;
    stringstream iss;
    vector<double> ThermoI, Thermo;
/////////////////////////read basic informations of the system
    iss.clear();
    fstream Information;
    Information.open("../Conf/Basic.txt",ios::in);
    getline(Information,line);
    getline(Information,line); iss.str(line); iss >> *NPart; iss >> *NBead; iss >> *NPoly;iss.clear();
    getline(Information,line);
    getline(Information,line); iss.str(line); iss >> *Box;iss.clear();
    getline(Information,line);
    getline(Information,line); iss.str(line); iss >> Delt_in[0];iss >> TProd_in[0];iss >> Delt_in[1];iss >> TProd_in[1];iss >> Delt_in[2];iss >> TProd_in[2];iss.clear();
    Information.close();
    cout<<Delt_in[0]<<" "<<TProd_in[0]<<" "<<Delt_in[1]<<" "<<TProd_in[1]<<" "<<Delt_in[2]<<" "<<TProd_in[2]<<endl;
    *Vol = pow(*Box,3);
    *Rho = *NPart/(*Vol);
}



void BasicInfo::cpu_time(string& Time){
    time_t nowtime;
    time(&nowtime);
    stringstream ftm;
    tm *p = localtime(&nowtime);
    ftm << "Date & Time: " <<setw(2) <<setfill('0')<<p->tm_year+1900 <<"-";
    ftm<<setw(2)<<setfill('0')<<p->tm_mon+1<<"-";
    ftm<<setw(2)<<setfill('0')<<p->tm_mday<<",";
    ftm<<setw(2)<<setfill('0')<<p->tm_hour<<":";
    ftm<<setw(2)<<setfill('0')<<p->tm_min<<":";
    ftm<<setw(2)<<setfill('0')<<p->tm_sec;
    Time = ftm.str();
}

void BasicInfo::Init(int *NFrame, int *NAtom,double *Delt_in,double *TProd,double *Delt,int *IT0,const int *T0Max,string &Trajname){

    *Delt = *Delt_in;
    *IT0=int(round(0.5*(*TProd)/(*T0Max*(*Delt))));
/////////////////////////read the number of frames and atoms
    auto input = chemfiles::Trajectory(Trajname);
    *NFrame =  input.nsteps();
    auto Frame = input.read();
    *NAtom = Frame.size();
    input.close();
}

