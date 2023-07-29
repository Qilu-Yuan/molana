#include"GetStat.h"

void Properity(int Switch,int *NPart,int *NBead,int *NPoly,double *Box,double *Vol, vector<vector<double>>& RX, vector<vector<double>>& PX){
    RDF(Switch,NPart,Box,Vol,RX,PX);
    SSF(Switch,NPart,Box,Vol,RX,PX);
    PSP(Switch,NPart,NBead,NPoly,Box,Vol,RX,PX);
    ANG(Switch,NPart,NBead,NPoly,Box,Vol,RX,PX);
    LB2(Switch,NPart,NBead,NPoly,Box,Vol,RX,PX);
    MID(Switch,NPart,NBead,NPoly,Box,Vol,RX,PX);

}


int main(int argc, char* argv[]){

    int NPart, NBead, NPoly; //Npart--number of particles, NBead--number of particles on a chain, Npoly-Number of polymers
    double Box,Vol,Rho,Temp,Pres;
	int NFrame,NAtom,NItem,NStep, FStatic, Foutput;

	time_t TimeStart, TimeFinal;
    string Time;

    TimeStart = time(NULL);

    cpu_time(Time);
    
    fstream Log;

    Log.open("LogStat",ios::out);Log<<Time<<endl<<endl;


    Init(&NPart,&NBead,&NPoly,&Box,&Vol,&Rho,&NItem,&NStep,&Temp,&Pres,&NFrame, &NAtom);

    if (NAtom != NPart){
        cout<<"NAtom from DCD file is inconsistent with NPart!"<<endl;
        return 1;
    }
    Log<<"---->>Basic information"               <<endl;
    Log<<"  Total particles               :    "<<NPart <<endl;
    Log<<"  Chain length                  :    "<<NBead <<endl;
    Log<<"  Chain number                  :    "<<NPoly <<endl;
    Log<<"  Box dimension                 :    "<<Box   <<endl;
    Log<<"  Number density                :    "<<Rho   <<endl<<endl;
    Log<<"  Number of frames in Trajectory:    "<<NFrame<<endl;
    Log<<"  Number of atoms  in Trajectory:    "<<NAtom <<endl<<endl;
	Log<<"  Data info in Thermo           :    "<<NItem<<" , "<< NStep<<endl;
	Log<<"  Temperature                   :    "<<Temp<<endl;
	Log<<"  Pressure                      :    "<<Pres<<endl<<endl;

    //Vector
    vector<vector<double>> RX(NPart, vector<double>(3, 0.0));
    vector<vector<double>> PX(NPoly, vector<double>(3, 0.0));

	FStatic = int(round(double(NFrame)/1000));
	Foutput = int(round(double(NFrame)/10));
    if (FStatic==0){
        FStatic = 1;
    } 


    Log.close();

    Properity(0,&NPart,&NBead,&NPoly,&Box,&Vol,RX,PX);

    fstream Log1;
    Log1.open("LogStat",ios::app);
    Log1<<endl<<"---->>calculation begins"<<endl;

    auto input = chemfiles::Trajectory("../Conf/CoordL.dcd",'r');

   
    for(int i = 1;i<=NFrame;i++){

        auto Frame = input.read();

        auto positions = Frame.positions();

        for (int i=0;i<NPart;i++){
            RX[i][0]=positions[i][0];
            RX[i][1]=positions[i][1];
            RX[i][2]=positions[i][2];
        }

        for (int i=0;i<NPoly;i++){
            PX[i][0] = 0;
            PX[i][1] = 0;
            PX[i][2] = 0;
            for (int j=0;j<NBead;j++){
                PX[i][0] += RX[i*NBead+j][0];
                PX[i][1] += RX[i*NBead+j][1];
                PX[i][2] += RX[i*NBead+j][2];
            }
            PX[i][0] = PX[i][0]/NBead;
            PX[i][1] = PX[i][1]/NBead;
            PX[i][2] = PX[i][2]/NBead;
        }
        
        if(i%10 ==0) cout<<i<<" "<<NFrame<<endl;

        if(i%FStatic==0) Properity(1,&NPart,&NBead,&NPoly,&Box,&Vol,RX,PX);
        
        if(i%Foutput==0){
            Properity(2,&NPart,&NBead,&NPoly,&Box,&Vol,RX,PX);
            Log1<<"#frames completed: "<< i <<endl;
        } 
    }
    input.close();
    Properity(2,&NPart,&NBead,&NPoly,&Box,&Vol,RX,PX);

    Log1<<"<<----calculation completed"<< endl <<endl;

    TimeFinal = time(NULL);

    Log1<<"Total time is "<<(TimeFinal-TimeStart)/60<<" minutes"<<endl;
    Log1.close();
}
