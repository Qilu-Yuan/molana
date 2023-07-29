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
// using chemfiles::Trajectory;

void cpu_time(string& Time){
    time_t nowtime;
    time(&nowtime);
    stringstream ftm;
    tm *p = localtime(&nowtime);
    ftm << "Date & Time: "<<p->tm_year+1900 <<"-"<<p->tm_mon+1<<"-"<<p->tm_mday<<","<<p->tm_hour<<":"<<p->tm_min<<":"<<p->tm_sec;
    Time = ftm.str();
}

void Init(int *NPart,int *NBead,int *NPoly,double *Box,double *Vol,double *Rho,int *NItem,int *NStep,double *Temp,double *Pres,int *NFrame, int *NAtom){

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
    Information.close();

    *Vol = pow(*Box,3);
    *Rho = *NPart/(*Vol);

    auto input = chemfiles::Trajectory("../Conf/CoordL.dcd");
    *NFrame =  input.nsteps();
    auto Frame = input.read();
    *NAtom = Frame.size();
    input.close();

/////////////////////////read and write the basic thermodynamics properits
    fstream InfoThermo;
    InfoThermo.open("../Conf/InfoThermo",ios::in);
    getline(InfoThermo,line); iss.str(line); iss >> *NItem; iss >> *NStep;iss.clear();
    InfoThermo.close();

    int NProp = *NItem-1,Step;
    Thermo.resize(NProp);
    ThermoI.resize(NProp);
    fill(Thermo.begin(),Thermo.end(),0);

    fstream ThermoIN;
    ThermoIN.open("../Conf/Thermo.txt",ios::in);
    getline(ThermoIN,line);
    for(int i=0;i<*NStep;i++){
        getline(ThermoIN,line); iss.str(line);
        iss >>Step;
        for(int i=0;i<NProp;i++){
            iss >> ThermoI[i];
            Thermo[i] += ThermoI[i];
        } 
        iss.clear();
    }
    ThermoIN.close();

    for(int i=0;i<NProp;i++) Thermo[i] = Thermo[i]/(*NStep);

    *Temp = Thermo[0];
    *Pres = Thermo[1];


    fstream OThermo;
    OThermo.open("Thermo",ios::out);
    OThermo<<scientific<<setprecision(8)<<*Rho<<" ";
    for(int i=0;i<NProp;i++) OThermo<<Thermo[i]<<" ";
    OThermo<<endl;
    OThermo.close();

}

