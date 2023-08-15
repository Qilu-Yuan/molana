#include "PSP.h"

void PSP::set(int *NBead,fstream& Log){
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

    Log<<"  To calculate PSP"<<endl;
}

void PSP::calcluate(int* NBead,int* NPoly,const vector<Vector3D>& RX,const vector<Vector3D>& PX){
    NStat = NStat + 1;

    double dx,dy,dz,r2,L;

    for(int i=0;i<*NPoly;i++){
        for(int j=0;j<(*NBead-1);j++){
            int k = i*(*NBead)+j;
            LX[j][0] = RX[k+1][0]-RX[k][0];
            LX[j][1] = RX[k+1][1]-RX[k][1];
            LX[j][2] = RX[k+1][2]-RX[k][2];
            L =sqrt(pow(LX[j][0],2)+pow(LX[j][1],2)+pow(LX[j][2],2));
            UX[j][0] = LX[j][0]/L;
            UX[j][1] = LX[j][1]/L;
            UX[j][2] = LX[j][2]/L;
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

void PSP::write(int *NBead,int *NPoly){
    fstream BCF;
    BCF.open("BCF",ios::out);

    for(int i=0;i<(*NBead-1);i++){
        double bcf = Fbb[i]/double(*NPoly)/double(NStat);
        BCF<<scientific<<setprecision(8)<<double(i)<<" "<<bcf<<endl;
    }
    BCF.close();

    double RE2 = SumRe2/double(NStat)/double(*NPoly);
    double RG2 = SumRg2/double(NStat)/double(*NPoly)/double(*NBead);
    fstream Size;
    Size.open("Size",ios::out);
    Size<<scientific<<setprecision(8)<<RE2<<" "<<RG2<<endl;
    Size.close();

}