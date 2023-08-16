#include"FSQ.h"

void FSQ::set(double *Box,int *TMax,fstream& Log){
    double Q,QrSet;

    this->FCOS1.resize(*TMax);
    this->FCOS2.resize(*TMax);
    fill(FCOS1.begin(),FCOS1.end(),0.0);
    fill(FCOS2.begin(),FCOS2.end(),0.0);

    this->Qx.resize(QImax);
    for(int i=0;i<QImax;i++){
        this->Qx[i].resize(3);
    }

    string line;
    stringstream iss;
    fstream Information;
    Information.open("./Qset.txt",ios::in);
    getline(Information,line); iss.str(line); iss >> QrSet; iss.clear();
    Information.close();

    FSQ::SetFSQ(Box,&Qn,Qx,&QImax,&Q,&QrSet);

    Log << "    Chosen Q         :  "<<Q <<endl;
    Log << "    Number of Vectors:  "<<Qn<<endl;

}

void FSQ::calcluate(int *NPart,int t,int dt,const vector<Vector3D>& RX,const vector<vector <double>> &RX0, const vector<vector <double>> &RY0, const vector<vector <double>> &RZ0){

    double dx, dy, dz, SetCOS,SumCOS;
    SumCOS = 0.0;
    for (int i=0;i<*NPart;i++){
        dx = RX[i][0] - RX0[i][t];
        dy = RX[i][1] - RY0[i][t];
        dz = RX[i][2] - RZ0[i][t];
        for (int k =0;k<Qn;k++){
            SetCOS = cos(double(Qx[k][0]*dx + Qx[k][1]*dy + Qx[k][2]*dz));
            SumCOS += SetCOS;
        }
    }
    FCOS1[dt] += SumCOS/Qn;
    FCOS2[dt] += pow(SumCOS/Qn,2);

}

void FSQ::write(int *NPart,int dt,const vector<int>& NDyna,double MDTime,fstream& file){
    double fsq,fsq1,fsq2,xsq;
    fsq = FCOS1[dt]/(double(NDyna[dt])*double(*NPart));;
    fsq1 = FCOS1[dt]/(double(NDyna[dt]));
    fsq2 = FCOS2[dt]/(double(NDyna[dt]));
    xsq = (fsq2 - fsq1*fsq1)/double(*NPart);
    file<<scientific<<setprecision(8)<<MDTime<<" "<<fsq<<" "<<xsq<<endl;
}

void FSQ::SetFSQ(double *Box,int *Qn,vector<vector<double>>& Qx,const int *QImax,double *Q,double *QrSet){

    int QiSet, K2,Q2;
    double DeltQ,DelMin;

    double unit=2*M_PI/(*Box);

    QiSet = int(*QrSet/unit)+1;

    DelMin = 1e10;

    for(int kx=0;kx <=QiSet;kx++){
        for(int ky=-QiSet;ky<=QiSet;ky++){
            for(int kz=-QiSet;kz<=QiSet;kz++){
                K2 = pow(kx,2) + pow(ky,2) + pow(kz,2);
                DeltQ = fabs(unit*sqrt(double(K2))-*QrSet);
                if (DeltQ<DelMin){
                    DelMin = DeltQ;
                    Q2 = K2;
                }
            }
        }
    }

    *Qn = 0;
    for(int kx=0;kx <=QiSet;kx++){
        for(int ky=-QiSet;ky<=QiSet;ky++){
            for(int kz=-QiSet;kz<=QiSet;kz++){
                K2 = pow(kx,2) + pow(ky,2) + pow(kz,2);
                if (K2==Q2 && *Qn<*QImax ){
                    Qx[*Qn][0] = unit*kx;
                    Qx[*Qn][1] = unit*ky;
                    Qx[*Qn][2] = unit*kz;
                    *Qn = *Qn + 1;
                }
            }
        }
    }
    *Q = unit*sqrt(double(Q2));

}

