#include "ANG.h"

void ANG(int Switch,int *NPart,int *NBead,int *NPoly,double *Box,double *Vol,vector<vector<double>>& RX,vector<vector<double>>& PX){ //to calculate the angle and angle distribution
    
    static int NStat,NAngle,NBin;
    int Bin;
    double Rmax,CosTheta,Theta;
    static double DelR, SumCOSTheta,SumTheta;
    
    double dxij[3],dxjk[3],rij2,rjk2;
    static vector<double> RBin,VBin,PTheta;
    // static double *RBin,*VBin,*PTheta;


    if (Switch == 1){
        NStat = NStat + 1;
		
		for(int p=0;p<*NPoly;p++){
            for(int b=0;b<(*NBead-2);b++){
                int i = p*(*NBead)+b;
                int j =i+1,k=i+2;
                dxij[0] = RX[i][0] - RX[j][0];
                dxij[1] = RX[i][1] - RX[j][1];
                dxij[2] = RX[i][2] - RX[j][2];
                dxjk[0] = RX[j][0] - RX[k][0];
                dxjk[1] = RX[j][1] - RX[k][1];
                dxjk[2] = RX[j][2] - RX[k][2];
                rij2 = pow(dxij[0],2) + pow(dxij[1],2) + pow(dxij[2],2);
                rjk2 = pow(dxjk[0],2) + pow(dxjk[1],2) + pow(dxjk[2],2);
                CosTheta = (dxij[0]*dxjk[0]+dxij[1]*dxjk[1]+dxij[2]*dxjk[2])/sqrt(rij2*rjk2);
                Theta = acos(CosTheta);
                SumCOSTheta = SumCOSTheta + CosTheta;
                SumTheta = SumTheta + Theta;
                Bin = int(Theta/DelR);
                if (Bin <= NBin){
                    PTheta[Bin] = PTheta[Bin] + 1.0;
                }
            }
        }         
    }
    else if (Switch == 0){
        NStat = 0;
        
        NAngle = (*NBead-2)*(*NPoly);

        Rmax = M_PI;
        DelR = Rmax/1.0e+2;

        NBin = int(Rmax/DelR);

        RBin.resize(NBin+1);
        VBin.resize(NBin+1);
        PTheta.resize(NBin+1);


        for(int i=0;i<=NBin;i++){
            RBin[i] = (i+0.5)*DelR;
            VBin[i] = 4.0/3.0*M_PI*(pow((i+1)*DelR,3)-pow(i*DelR,3));
            PTheta[i] = 0.0;
        }

        SumCOSTheta = 0.0;
        SumTheta = 0.0;

        fstream Log;
        Log.open("LogStat",ios::app);
        Log << "  To calculate ANG" <<endl;
        Log.close();
    }
    else if (Switch ==2 ){
        fstream PANG;
        PANG.open("PANG",ios::out);

        for(int i=0;i<(NBin);i++){
            double pang = 4*M_PI*pow(RBin[i],2)*PTheta[i]/double(NAngle)/VBin[i]/double(NStat);
            // cout<<NStat<<endl;
            PANG<<scientific<<setprecision(8)<<RBin[i]<<" "<<pang<<endl;
        }
        PANG.close();
        double cosang = SumCOSTheta/double(NAngle)/double(NStat);
        double ang = SumTheta/double(NAngle)/double(NStat);

        fstream ANG;
        ANG.open("ANG",ios::out);
        ANG<<scientific<<setprecision(8)<<cosang<<" "<<ang<<endl;
        ANG.close();
    }
}
