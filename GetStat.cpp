#include"GetStat.h"

int main(){

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


	FStatic = int(round(double(NFrame)/1000));
	Foutput = int(round(double(NFrame)/10));
    if (FStatic==0){
        FStatic = 1;
    } 

    Properties Properties1(&NPart,&NBead,&NPoly,&Box,&Vol,Log);

    Log<<endl<<"---->>calculation begins"<<endl;

    auto input = chemfiles::Trajectory("../Conf/CoordL.dcd",'r');

    vector<Vector3D> PX(NPoly, Vector3D(0.0, 0.0, 0.0));
    
    for(int i = 1;i<=NFrame;i++){

        const auto Frame = input.read();

        auto RX1 = Frame.positions();

        for (int i=0;i<NPoly;i++){
            PX[i][0] = 0;
            PX[i][1] = 0;
            PX[i][2] = 0;
            for (int j=0;j<NBead;j++){
                int k = i*NBead+j;
                PX[i][0] += RX1[k][0];
                PX[i][1] += RX1[k][1];
                PX[i][2] += RX1[k][2];
            }
            PX[i][0] = PX[i][0]/NBead;
            PX[i][1] = PX[i][1]/NBead;
            PX[i][2] = PX[i][2]/NBead;
        }

        if(i%10 ==0) cout<<i<<" "<<NFrame<<endl;
         if(i%FStatic==0) Properties1.calcluate(RX1,PX);
        
        if(i%Foutput==0){
            Properties1.write();
            // Ang1.write(&NBead,&NPoly);
            Log<<"#frames completed: "<< i <<endl;
        } 
    }
    input.close();
    // Ang1.write(&NBead,&NPoly);//
    Properties1.write();

    Log<<"<<----calculation completed"<< endl <<endl;

    TimeFinal = time(NULL);

    Log<<"Total time is "<<double(TimeFinal-TimeStart)/60.0<<" minutes"<<endl;
    Log.close();
}
