#include "MID.h"

void MID(int Switch,int *NPart,int *NBead,int *NPoly,double *Box,double *Vol,vector<vector<double>>& RX,vector<vector<double>>& PX){ //to calculate compute mean-square internal distance
    
    static int NStat;
    double dx[3],r2,MID;
    static vector<double> Acc,Sum,Rn2;


    if (Switch == 1){
        NStat = NStat + 1;

        for(int i=0;i<(*NBead)-1;i++){
            Sum[i] = 0.0;
            Acc[i] = 0.0;
        }
		for(int p=0;p<*NPoly;p++){
            for(int b1=0;b1<(*NBead)-1;b1++){
                int i = p*(*NBead)+b1;
                for(int b2 =b1+1;b2<(*NBead);b2++ ){
                    int j = p*(*NBead)+b2;
                    dx[0] = RX[j][0] - RX[i][0];
                    dx[1] = RX[j][1] - RX[i][1];
                    dx[2] = RX[j][2] - RX[i][2];
                    r2 = pow(dx[0],2) + pow(dx[1],2) + pow(dx[2],2);
                    int n = b2-b1-1;//this is the ith element in the vector
                    Sum[n] = Sum[n] + r2;
                    Acc[n] = Acc[n] + 1.0;
                }
            }
        }
        for(int i=0;i<*NBead-1;i++){
            Rn2[i] = Rn2[i] + Sum[i]/Acc[i];
        }         
    }
    else if (Switch == 0){
        NStat = 0;
        
        Acc.resize(*NBead-1);
        Sum.resize(*NBead-1);
        Rn2.resize(*NBead-1);

        fill(Rn2.begin(),Rn2.end(),0.0);

        fstream Log;
        Log.open("LogStat",ios::app);
        Log << "  To calculate MID" <<endl;
        Log.close();
    }
    else if (Switch ==2 ){
        fstream MID;
        MID.open("MID",ios::out);

        for(int i=0;i<(*NBead)-1;i++){
            // double mid = Rn2[i]/double(NStat)/double(i+1);
            double mid = Rn2[i]/double(NStat)/double(i+1);
            MID<<scientific<<setprecision(8)<<double(i+1)<<" "<<mid<<endl;

        }
        MID.close();
    }
}
