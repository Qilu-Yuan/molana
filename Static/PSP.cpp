#include "PSP.h"


void PSP(int Switch,int *NPart,int *NBead,int *NPoly,double *Box,double *Vol,vector<vector<double>>& RX,vector<vector<double>>& PX){ //to calculate polymer specific properity
    
    static int NStat;
    
    double dx,dy,dz,r2,L2;
    double BCF,RE2,RG2;
    static double SumRe2,SumRg2;
    static vector<double> Fbb;
    // static double *Fbb;
    static vector<vector<double>> LX,UX;
    // static double **LX,**UX;


    if (Switch == 1){
        NStat = NStat + 1;
		
		for(int i=0;i<*NPoly;i++){
            for(int j=0;j<(*NBead-1);j++){
                int k = i*(*NBead)+j;
                LX[j][0] = RX[k+1][0]-RX[k][0];
                LX[j][1] = RX[k+1][1]-RX[k][1];
                LX[j][2] = RX[k+1][2]-RX[k][2];
                L2 = pow(LX[j][0],2)+pow(LX[j][1],2)+pow(LX[j][2],2);
                UX[j][0] = LX[j][0]/sqrt(L2);
                UX[j][1] = LX[j][1]/sqrt(L2);
                UX[j][2] = LX[j][2]/sqrt(L2);
            }

            for (int j=0;j<(*NBead-1);j++){
                Fbb[j] = Fbb[j] + UX[j][0]*UX[0][0]+UX[j][1]*UX[0][1]+UX[j][2]*UX[0][2];
            }

            dx = RX[i*(*NBead)][0] - RX[i*(*NBead)+(*NBead-1)][0];
            dy = RX[i*(*NBead)][1] - RX[i*(*NBead)+(*NBead-1)][1];
            dz = RX[i*(*NBead)][2] - RX[i*(*NBead)+(*NBead-1)][2];
            r2 = pow(dx,2)+pow(dy,2)+pow(dz,2);
            SumRe2 = SumRe2 + r2;

            for(int j=0;j<*NBead;j++){
                dx = RX[i*(*NBead)+j][0] - PX[i][0];
                dy = RX[i*(*NBead)+j][1] - PX[i][1];
                dz = RX[i*(*NBead)+j][2] - PX[i][2];
                r2 = pow(dx,2)+pow(dy,2)+pow(dz,2);
                SumRg2 = SumRg2 + r2;
            }
        }         
    }
    else if (Switch == 0){
        NStat = 0;
        
        SumRe2 = 0.0;
        SumRg2 = 0.0;
        // cout<<"rmax is :"<<Rmax<<endl;

        Fbb.resize(*NBead);
        fill(Fbb.begin(),Fbb.end(),0.0);
        LX.resize(*NBead);
        UX.resize(*NBead);
        for(int i=0;i<*NBead;i++){
            LX[i].resize(3);
            UX[i].resize(3);
        }        
        // Fbb = new double[*NBead];

        fstream Log;
        Log.open("LogStat",ios::app);
        Log << "  To calculate PSP" <<endl;
        Log.close();
    }
    else if (Switch ==2 ){
        fstream BCF;
        BCF.open("BCF",ios::out);

        for(int i=0;i<(*NBead-1);i++){
            double bcf = Fbb[i]/double(NStat)/double(*NPoly);
            // cout<<NStat<<endl;
            BCF<<scientific<<setprecision(8)<<double(i)<<" "<<bcf<<endl;
        }
        BCF.close();
        RE2 = SumRe2/double(NStat)/double(*NPoly);
        RG2 = SumRg2/double(NStat)/double(*NPoly)/double(*NBead);
        fstream Size;
        Size.open("Size",ios::out);
        Size<<scientific<<setprecision(8)<<RE2<<" "<<RG2<<endl;
        Size.close();
    }

}
