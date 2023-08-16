#include"MSD.h"

void MSD::set(int *TMax){
    this->FR2.resize(*TMax);
    this->FR4.resize(*TMax);
    fill(FR2.begin(),FR2.end(),0.0);
    fill(FR4.begin(),FR4.end(),0.0);
}

void MSD::calcluate(int *NPart,int t,int dt,const vector<Vector3D>& RX,const vector<vector <double>> &RX0, const vector<vector <double>> &RY0, const vector<vector <double>> &RZ0){
    double dx, dy, dz,r2, r4, Sum1, Sum2;
    Sum1 = 0.0;
    Sum2 = 0.0;

    for (int i=0;i<*NPart;i++){
        dx = RX[i][0] - RX0[i][t];
        dy = RX[i][1] - RY0[i][t];
        dz = RX[i][2] - RZ0[i][t];
        r2 = dx*dx + dy*dy + dz*dz;
        r4 = r2*r2;
        Sum1 += r2;
        Sum2 += r4;
    }
    FR2[dt] += Sum1;
    FR4[dt] += Sum2;

}

void MSD::write(int *NPart,int dt,const vector<int>& NDyna,double MDTime,fstream& file){
    double r2, r4, a2;
    r2 = FR2[dt]/(double(NDyna[dt])*double(*NPart));
    r4 = FR4[dt]/(double(NDyna[dt])*double(*NPart));
    a2 = 0.0;
    if (r2>0){
        a2 = 0.6*r4/(r2*r2)-1;
    } 
    file<<scientific<<setprecision(8)<<MDTime<<" "<<r2<<" "<<a2<<endl;
}
