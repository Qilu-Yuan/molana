#include "Property.h"

Properties::Properties(int *NPart_in,int *NBead_in,int *NPoly_in,double *Box_in,double *Vol_in,double *Rho_in,int *TMax_in,int *IT0_in,int *T0Max_in,double *Delt_in,fstream& Log,char& type){
    NPart = NPart_in;
    NBead = NBead_in;
    NPoly = NPoly_in;
    Box = Box_in;
    Vol = Vol_in;
    Rho = Rho_in;
    TMax = TMax_in;
    IT0 = IT0_in;
    T0Max = T0Max_in;
    Delt = Delt_in;

    Nt = 0;
    Nt0 = 0;

    NDyna.resize(*TMax);
    fill(NDyna.begin(),NDyna.end(),0.0);
    ttv0.resize(*T0Max);

    RX0.resize(*NPart);RY0.resize(*NPart);RZ0.resize(*NPart);
    for (int i=0;i<*NPart;i++){
        RX0[i].resize(*T0Max);
        RY0[i].resize(*T0Max);
        RZ0[i].resize(*T0Max);
    }

    EX0.resize(*NPoly);EY0.resize(*NPoly);EZ0.resize(*NPoly);
    EX.resize(*NPoly);
    for(int i=0;i<*NPoly;i++){
        EX0[i].resize(*T0Max);
        EY0[i].resize(*T0Max);
        EZ0[i].resize(*T0Max);
    }

    stringstream prop;
    prop<<"MSD_"<<type;
    name_msd = prop.str();prop.clear();prop.str("");
    prop<<"FSQ_"<<type;
    name_fsq = prop.str();prop.clear();prop.str("");
    prop<<"ASL_"<<type;
    name_asl = prop.str();prop.clear();prop.str("");
    prop<<"X4S_"<<type;
    name_x4s = prop.str();prop.clear();prop.str("");
    prop<<"EEV_"<<type;
    name_eev = prop.str();prop.clear();prop.str("");
//<<<<<<<<<<<------------properities------------->>>>>>>>>>>>>>>   
    Log<<"  To calculate "<<name_msd<<endl;
    MSD1.set(this->TMax);
    Log<<"  To calculate "<<name_fsq<<endl;
    FSQ1.set(Box,this->TMax,Log);
    Log<<"  To calculate "<<name_asl<<endl;
    ASL1.set(NPart,this->TMax);
    Log<<"  To calculate "<<name_x4s<<endl;
    X4S1.set(this->TMax);
    Log<<"  To calculate "<<name_eev<<endl;
    EEV1.set(this->TMax);
//<<<<<<<<<<<--------------------------------->>>>>>>>>>>>>>>
}

void Properties::calculate(const vector<Vector3D>& RX,const vector<Vector3D>& PX){
    // int ttel,dt;
    Nt+=1;
		
    for(int i=0;i<*NPoly;i++){
        EX[i][0] = RX[(i)*(*NBead)+*NBead-1][0]-RX[i*(*NBead)][0];
        EX[i][1] = RX[(i)*(*NBead)+*NBead-1][1]-RX[i*(*NBead)][1];
        EX[i][2] = RX[(i)*(*NBead)+*NBead-1][2]-RX[i*(*NBead)][2];
    }
    if (Nt%(*IT0) ==0){
        Nt0+=1;
        int ttel = (Nt0-1)%(*T0Max) ;
        ttv0[ttel] = Nt;
        for (int i=0;i<*NPart;i++){
            RX0[i][ttel] = RX[i][0];
            RY0[i][ttel] = RX[i][1];
            RZ0[i][ttel] = RX[i][2];
        }
        for (int i=0;i<*NPoly;i++){
            EX0[i][ttel] = EX[i][0];
            EY0[i][ttel] = EX[i][1];
            EZ0[i][ttel] = EX[i][2];
        }
    }
    for (int t= 0;t<min(Nt0,(*T0Max));t++){
        int dt = Nt - ttv0[t];
        if (dt>=0 && dt<(*TMax)){
            NDyna[dt] += 1;
//<<<<<<<<<<<------------properities------------->>>>>>>>>>>>>>>
            MSD1.calcluate(NPart,t,dt,RX,RX0,RY0,RZ0);
            FSQ1.calcluate(NPart,t,dt,RX,RX0,RY0,RZ0);
            ASL1.calcluate(NPart,Box,t,dt,RX,RX0,RY0,RZ0);
            X4S1.calcluate(NPart,Box,t,dt,RX,RX0,RY0,RZ0);
            EEV1.calcluate(NPoly,t,dt,EX,EX0,EY0,EZ0);
//<<<<<<<<<<--------------------------------->>>>>>>>>>>>>>>>>>

        }
    }
}

void Properties::write(){
    fstream msd,fsq,asl,x4s,eev;

    msd.open(name_msd,ios::out);
    fsq.open(name_fsq,ios::out);
    asl.open(name_asl,ios::out);
    x4s.open(name_x4s,ios::out);
    eev.open(name_eev,ios::out);
    for(int dt=0;dt<(*TMax);dt++){
        MDTime = (*Delt)*(dt);
        if (NDyna[dt] > 0){
//<<<<<<<<<<<------------properities------------->>>>>>>>>>>>>>>
            MSD1.write(NPart,dt,NDyna,MDTime,msd);
            FSQ1.write(NPart,dt,NDyna,MDTime,fsq);
            ASL1.write(dt,NDyna,MDTime,asl);
            X4S1.write(NPart,Vol,dt,NDyna,MDTime,x4s);
            EEV1.write(NPoly,dt,NDyna,MDTime,eev);
//<<<<<<<<<<--------------------------------->>>>>>>>>>>>>>>>>>
        }
    }
    msd.close();
    fsq.close();
    asl.close();
    x4s.close();
    eev.close();
}
