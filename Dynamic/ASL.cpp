#include"ASL.h"

void ASL::Get_Pmb(vector<double>& r2,vector<int>& Pmb, int *NPart){

    int N,LOGNB2,NN,Pmbtmp;
    double r2tmp;

    N=*NPart;
    for(int i=0;i<N;i++){
        Pmb[i]=i;
    }
    LOGNB2 = int(log(double(N))/log(2.0)+1e-5);
    NN=N;
    for(int h=0;h<LOGNB2;h++){
        NN=NN/2;
        int k=N-NN;
        for(int j=0;j<k;j++){
            int i =j;
            do
            {
                int L = i + NN;
                if (r2[L] > r2[i])
                {
                    r2tmp = r2[i];
                    r2[i] = r2[L];
                    r2[L] = r2tmp;
                    Pmbtmp = Pmb[i];
                    Pmb[i] = Pmb[L];
                    Pmb[L] = Pmbtmp;
                    i -= NN;
                }
            } while (i >= 0);
        }
    }
}


void ASL::set(int *NPart,int *TMax){

    ST.resize((*TMax));
    fill(ST.begin(),ST.end(),0.0);
    Nmb = int(*NPart*0.065);
    DelR2 = pow(0.55,2);
    r2.resize(*NPart);Pmb.resize(*NPart);
    NIT.resize(Nmb);L.resize(Nmb);Num.resize(Nmb);
}

void ASL::calcluate(int *NPart,double *Box,int t,int dt,const vector<Vector3D>& RX,const vector<vector <double>> &RX0, const vector<vector <double>> &RY0, const vector<vector <double>> &RZ0){
    double dx, dy, dz,drjk2,drkj2;

    for (int i=0;i<(*NPart);i++){
        dx = RX[i][0] - RX0[i][t];
        dy = RX[i][1] - RY0[i][t];
        dz = RX[i][2] - RZ0[i][t];
        r2[i] = dx*dx + dy*dy + dz*dz;
    }
// 	====string-like collective motion
// 	====get mobile particles
    ASL::Get_Pmb(r2,Pmb,NPart);   
// 	====sort the clusters
    for(int i=0;i<Nmb;i++){
        L[i]=i;
    }
    for (int i=0;i<Nmb-1;i++){
        if (i == L[i]){
            int j = i;
            do{
                for (int k=i;k<Nmb;k++){
                    int LK=L[k];
                    if (LK == k){
                        dx = RX[Pmb[j]][0] - RX0[Pmb[k]][t];
                        dy = RX[Pmb[j]][1] - RY0[Pmb[k]][t];
                        dz = RX[Pmb[j]][2] - RZ0[Pmb[k]][t];
                        dx = dx - *Box*round(dx/(*Box));
                        dy = dy - *Box*round(dy/(*Box));
                        dz = dz - *Box*round(dz/(*Box));
                        drjk2 = dx*dx + dy*dy + dz*dz;
                        dx = RX[Pmb[k]][0] - RX0[Pmb[j]][t];
                        dy = RX[Pmb[k]][1] - RY0[Pmb[j]][t];
                        dz = RX[Pmb[k]][2] - RZ0[Pmb[j]][t];
                        dx = dx - *Box*round(dx/(*Box));
                        dy = dy - *Box*round(dy/(*Box));
                        dz = dz - *Box*round(dz/(*Box));
                        drkj2 = dx*dx + dy*dy + dz*dz;
                        if (min(drjk2,drkj2)<DelR2){
                            L[k] = L[j];
                            L[j] = LK;
                        }
                    }
                }
                j = L[j];
            }while (j!=i);
        }
    }
// 	====count number in a cluster containing particle IT
    for (int it=0;it<Nmb;it++){
        NIT[it] = 1;
        int LIT = L[it];
        while(LIT!=it){
            NIT[it] += 1;
            LIT = L[LIT];
        }
    }
	// ====count number of total particles for a given cluster size
    fill(Num.begin(),Num.end(),0);


    for (int it=0;it<Nmb;it++){
        for(int i=0;i<Nmb;i++){
            if (NIT[it] == i+1){
                Num[i] += 1;
            }
        }
    }

    NCluster = 0;
    for (int i=0;i<Nmb;i++){
        NCluster += Num[i]/(i+1);
    }
    ST[dt] += double(Nmb)/double(NCluster);
}

void ASL::write(int dt,const vector<int>& NDyna,double MDTime,fstream& file){
    double asl = ST[dt]/NDyna[dt];
    file<<scientific<<setprecision(8)<<MDTime<<" "<<asl<<endl;
}
