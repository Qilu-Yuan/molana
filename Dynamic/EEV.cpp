#include"EEV.h"

void EEV::set(int *TMax){
    SumEEV.resize(*TMax);
    fill(SumEEV.begin(),SumEEV.end(),0.0);
}

void EEV::calcluate(int *NPoly,int t,int dt,const vector<Vector3D>&EX,const vector<vector<double>>&EX0,const vector<vector<double>>&EY0,const vector<vector<double>>&EZ0){

    for (int i=0;i<(*NPoly);i++){
        double EProd = EX[i][0]*EX0[i][t] + EX[i][1]*EY0[i][t] + EX[i][2]*EZ0[i][t];
        SumEEV[dt] += EProd;
    }
}

void EEV::write(int *NPoly,int dt,const vector<int>& NDyna,double MDTime,fstream& file){

    double eev;
    static double eev0;
    eev = SumEEV[dt]/double(NDyna[dt])/double(*NPoly);
    if (dt ==0){
        eev0=eev;
    }
    file<<scientific<<setprecision(8)<<MDTime<<" "<<eev/eev0<<endl;


}

