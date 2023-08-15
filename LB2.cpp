#include "LB2.h"

void LB2::set(int *NBead,int *NPoly,fstream& Log){
    NStat = 0;
    NBond = (*NBead-1)*(*NPoly);

    double Rmax = 2.0;
    DelR = Rmax/2.0e+2;

    NBin = int(Rmax/DelR);

    RBin.resize(NBin+1);
    VBin.resize(NBin+1);
    PR2.resize(NBin+1);

    for(int i=0;i<=NBin;i++){
        RBin[i] = (i+0.5)*DelR;
        VBin[i] = 4.0/3.0*M_PI*(pow((i+1)*DelR,3)-pow(i*DelR,3));
        PR2[i] = 0.0;
    }
    SumR2 = 0.0;

    Log << "  To calculate LB2" <<endl;
}

void LB2::calcluate(int* NBead,int* NPoly,const vector<Vector3D>& RX){
    NStat = NStat + 1;
    double dx[3],r2;
    int Bin;

    for(int p=0;p<*NPoly;p++){
        for(int b=0;b<(*NBead-1);b++){
            int i = p*(*NBead)+b;
            int j =i+1;
            
            dx[0] = RX[j][0] - RX[i][0];
            dx[1] = RX[j][1] - RX[i][1];
            dx[2] = RX[j][2] - RX[i][2];
            r2 = pow(dx[0],2) + pow(dx[1],2) + pow(dx[2],2);
            SumR2 = SumR2 + r2;
            Bin = int(r2/DelR);
            if (Bin <= NBin){
                PR2[Bin] = PR2[Bin] + 1.0;
            }
        }
    }
}

void LB2::write(int *NBead){
    fstream PLB2;
    PLB2.open("PLB2",ios::out);

    for(int i=0;i<=NBin;i++){
        double plb2 = 4*M_PI*pow(RBin[i],2)*PR2[i]/double(NBond)/VBin[i]/double(NStat);
        PLB2<<scientific<<setprecision(8)<<RBin[i]<<" "<<plb2<<endl;
    }
    PLB2.close();

    double lb2 = SumR2/double(NBond)/double(NStat);

    fstream LB2;
    LB2.open("LB2",ios::out);
    LB2<<scientific<<setprecision(8)<<lb2<<endl;
    LB2.close();
}   