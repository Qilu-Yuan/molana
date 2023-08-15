#include"RDF.h"

void RDF::set(double *Box,fstream& Log){
    NStat = 0;
    double Rmax = *Box /2;

    DelR = 0.01;
    NBin = int(Rmax/DelR);
    
    RBin.resize(NBin+1);
    VBin.resize(NBin+1);
    Gbb.resize(NBin+1);

    fill(Gbb.begin(),Gbb.end(),0.0);

    for(int i=0;i<=NBin;i++){
        RBin[i]= (i+0.5)*DelR;
        VBin[i]= 4*M_PI*pow(DelR,3)*(pow(double(i+1),3)-pow(double(i),3))/3;
    }

    Log << "  To calculate RDF" <<endl;
}

void RDF::calcluate(double *Box,int *NPart,const vector<Vector3D>& RX){
    int Bin;
    double dx,dy,dz,r2;
    NStat = NStat + 1;
    
    for(int i=0;i<*NPart;i++){
        for(int j=i+1;j<*NPart;j++){
            dx = RX[i][0] - RX[j][0];
            dy = RX[i][1] - RX[j][1];
            dz = RX[i][2] - RX[j][2];
            dx -= (*Box)*round(dx/(*Box));
            dy -= (*Box)*round(dy/(*Box));
            dz -= (*Box)*round(dz/(*Box));
            r2 = pow(dx,2) + pow(dy,2) + pow(dz,2);
            Bin = int(sqrt(r2)*100);
            if (Bin<=NBin){
                Gbb[Bin] = Gbb[Bin] + 1;
            } 
        }
    }         
}

void RDF::write(double *Vol,int *NPart){
        fstream RDF;
        RDF.open("RDF",ios::out);
        for(int i=0;i<NBin;i++){
            double rdf = (2*(*Vol)*Gbb[i])/(VBin[i]*NStat*pow(*NPart,2));
            // cout<<NStat<<endl;
            RDF<<scientific<<setprecision(8)<<RBin[i]<<" "<<rdf<<endl;
        }
        RDF.close();
}
