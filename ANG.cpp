#include"ANG.h"

void ANG::set(int *NBead,int *NPoly,fstream& Log){
    this->NStat = 0;
    this->NAngle = (*NBead-2)*(*NPoly);
    double Rmax = M_PI;
    this->DelR = Rmax/1.0e+2;
    this->NBin = int(Rmax/this->DelR);
    this->RBin.resize(this->NBin+1);
    this->VBin.resize(this->NBin+1);
    this->PTheta.resize(this->NBin+1);

    for(int i=0;i<=this->NBin;i++){
        this->RBin[i] = (i+0.5)*this->DelR;
        this->VBin[i] = 4.0/3.0*M_PI*(pow((i+1)*this->DelR,3)-pow(i*this->DelR,3));
        this->PTheta[i] = 0.0;
    }

    this->SumCOSTheta = 0.0;
    this->SumTheta = 0.0;

    Log << "  To calculate ANG" <<endl;

}
void ANG::calcluate(int *NBead,int *NPoly,const vector<Vector3D>& RX){ //to calculate the angle and angle distribution
    double dxij[3],dxjk[3],rij2,rjk2;    
    double CosTheta,Theta;

    this->NStat += 1;
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
            this->SumCOSTheta = this->SumCOSTheta + CosTheta;
            this->SumTheta = this->SumTheta + Theta;

            int Bin = int(Theta/(this->DelR));
            if (Bin <= (this->NBin)){
                this->PTheta[Bin] = this->PTheta[Bin] + 1.0;
            }
        }
    }
}

void ANG::write(int *NBead,int *NPoly){
        fstream PANG;
        PANG.open("PANG",ios::out);

        for(int i=0;i<(this->NBin);i++){
            double pang = 4*M_PI*pow(this->RBin[i],2)*this->PTheta[i]/(this->NStat*(*NPoly)*(*NBead-2)*this->VBin[i]);
            PANG<<scientific<<setprecision(8) << this->RBin[i] << " " << pang << endl;
        }
        PANG.close();
        double AvgCOSTheta = this->SumCOSTheta/(this->NStat*this->NAngle);
        double AvgTheta = this->SumTheta/(this->NStat*this->NAngle);
        fstream ANG;
        ANG.open("ANG",ios::out);
        ANG<<scientific<<setprecision(8) << AvgCOSTheta <<" "<< AvgTheta << endl;
        ANG.close();

}