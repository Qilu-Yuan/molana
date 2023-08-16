#include"X4S.h"

void X4S::set(int *TMax){
    Qs1.resize(*TMax);
    Qs2.resize(*TMax);
    fill(Qs1.begin(),Qs1.end(),0.0);
    fill(Qs2.begin(),Qs2.end(),0.0);
}

void X4S::calcluate(int *NPart,double *Box,int t,int dt,const vector<Vector3D>& RX,const vector<vector <double>> &RX0, const vector<vector <double>> &RY0, const vector<vector <double>> &RZ0){
    double dx, dy, dz,r2,SumQs;

    for (int i=0;i<*NPart;i++){
        dx = RX[i][0] - RX0[i][t];
        dy = RX[i][1] - RY0[i][t];
        dz = RX[i][2] - RZ0[i][t];
        dx = dx - *Box*round(dx/(*Box));
        dy = dy - *Box*round(dy/(*Box));
        dz = dz - *Box*round(dz/(*Box));
        r2 = dx*dx + dy*dy + dz*dz;
        if (r2<a2){
            SumQs += 1.0;
        }
    }
    Qs1[dt]=Qs1[dt]+SumQs;
    Qs2[dt]=Qs2[dt]+pow(SumQs,2);
}

void X4S::write(int *NPart,double *Vol,int dt,const vector<int>& NDyna,double MDTime,fstream& file){

    double qss = Qs1[dt]/(double(NDyna[dt])*double(*NPart));
                //X4S
    double qs1i = Qs1[dt]/(double(NDyna[dt]));
    double qs2i = Qs2[dt]/(double(NDyna[dt]));
    double x4s = *Vol*(qs2i-pow(qs1i,2))/pow(double(*NPart),2);
    file<<scientific<<setprecision(8)<<MDTime<<" "<<qss<<" "<<x4s<<endl;
}