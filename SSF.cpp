#include"SSF.h"

void SetSSF(double *Box,vector<int> &Qn,vector<int> &Qc,vector<vector<double>>& Qx,vector<vector<double>>& Qy,vector<vector<double>>& Qz,int *NQ,int Qmax,int QIMax,double Qrh){
    int Qih,Nk,Kx,Ky,Kz,K2,Ks;
    vector<int> Kn,Kc;

    Qih = int(Qrh*(*Box)/(2*M_PI))+1;

    Kn.resize(Qmax);
    Kc.resize(Qmax);

    Ks = 8;
    int K = 0;

    for(int i=Ks;i<=Ks*Qih;i++){
        Kc[K] = int(round(pow((double(i)/double(Ks)),2)));
        K += 1;
    }

    Nk = 0 ;
    for(int i=1;i<K;i++){
        if (Kc[i]-Kc[i-1]>0){
            Nk = Nk + 1;
            Kc[Nk] = Kc[i];
        }
    }

    for (int i = 0;i<Nk;i++){
        Kn[i] = 0;
        for (Kx = 0;Kx<=Qih;Kx++){
            for (Ky = -Qih;Ky<=Qih;Ky++){
                for (Kz = -Qih;Kz<=Qih;Kz++){
                    K2 = pow(Kx,2) + pow(Ky,2) + pow(Kz,2);
                    if (K2==Kc[i] && Kn[i]<QIMax){
                        Kn[i] = Kn[i] + 1;
                        Qx[i][Kn[i]-1] = Kx;
                        Qy[i][Kn[i]-1] = Ky;
                        Qz[i][Kn[i]-1] = Kz;
                    }
                }
            }
        }
    }
    //	====set number of wave numbers and wave vectors for S(q)
	*NQ = 0;
    for (int i =0;i<Nk;i++){
        if (Kn[i]>0){
            *NQ = *NQ + 1;
            Qn[*NQ-1] = Kn[i];
            Qc[*NQ-1] = Kc[i];
            for (int j=0;j<Kn[i];j++){
                Qx[*NQ-1][j] = 2*M_PI/(*Box)*Qx[i][j];
                Qy[*NQ-1][j] = 2*M_PI/(*Box)*Qy[i][j];
                Qz[*NQ-1][j] = 2*M_PI/(*Box)*Qz[i][j];
            }
        }
    }

}



void SSF(int Switch,int *NPart,double *Box,double *Vol,vector<vector<double>>& RX,vector<vector<double>>& PX){
    const int Qmax = 10000, QIMax= 20; // Maximum number of Q values
    static int NQ;
    const double Qrh=15.0;
    static vector<int> Qn,Qc;
    static vector<vector<double>> Qx,Qy,Qz;

    double Kx, Ky, Kz, SumCOS, SumSIN, SSF;
    static int NStat;

    static double *Q,*Sbb;

    if (Switch == 1){
        NStat = NStat + 1;
		
		for(int i=0;i<NQ;i++){
            for(int j=0;j<Qn[i];j++){
                Kx = Qx[i][j];
                Ky = Qy[i][j];
                Kz = Qz[i][j];

                SumCOS = 0.0;
                SumSIN = 0.0;
                for(int k=0;k<*NPart;k++){
                    SumCOS = SumCOS + cos(Kx*RX[k][0] + Ky*RX[k][1] + Kz*RX[k][2]);
                    SumSIN = SumSIN + sin(Kx*RX[k][0] + Ky*RX[k][1] + Kz*RX[k][2]);
                }
                Sbb[i] = Sbb[i] + pow(SumCOS,2) + pow(SumSIN,2);
            }
        }
    }         
    else if (Switch == 0){
        NStat = 0;

        Qn.resize(Qmax);
        Qc.resize(Qmax);
        
        Qx.resize(Qmax);
        Qy.resize(Qmax);
        Qz.resize(Qmax);
        for(int i=0;i<Qmax;i++){
            Qx[i].resize(QIMax);
            Qy[i].resize(QIMax);
            Qz[i].resize(QIMax);
        }


        SetSSF(Box,Qn,Qc,Qx,Qy,Qz,&NQ,Qmax,QIMax,Qrh);

        Q = (double*) calloc(NQ, sizeof(double));
        Sbb = (double*) calloc(NQ, sizeof(double));
        
        for(int i=0;i<NQ;i++){
            Q[i] = 2*M_PI/(*Box)*sqrt(double(Qc[i]));
            Sbb[i] = 0.0;
        }

        fstream Log;
        Log.open("LogStat",ios::app);
        Log << "  To calculate SSF" <<endl;
        Log << "     Number of Q values = " << NQ <<endl;
        Log << "     Minimun & Maximum Q value = " << Q[0]<<" "<<Q[NQ-1] <<endl;
        Log.close();
    }
    else if (Switch ==2 ){
        fstream SSF;
        SSF.open("SSF",ios::out);
        for(int i=0;i<NQ;i++){
            double ssf = Sbb[i]/(double(NStat)*double(*NPart)*double(Qn[i]));

            SSF<<scientific<<setprecision(8)<<Q[i]<<" "<<ssf<<endl;
        }
        SSF.close();
    }
}
