#include"GetDyna.h"

int main(){
    int NPart, NBead, NPoly; //Npart--number of particles, NBead--number of particles on a chain, Npoly-Number of polymers
    double Box,Vol,Rho;
	int NFrame,NAtom,NItem,NStep, FStatic, Foutput;
    double Delt_in[3], TProd_in[3],Delt;
    int IT0,TMax;
    int T0Max = 100;//the number of anverage times
    

    char Type[] ="LMS"; //to calculate the three types of TCF:long time(L);middle time(M);short time(S)

    BasicInfo Basic(&NPart,&NBead,&NPoly,&Box,&Vol,&Rho,Delt_in,TProd_in);
   
    for(int i=0;i<3;i++){
////////////get the time        
        time_t TimeStart, TimeFinal;
        string Time;
        TimeStart = time(NULL);
        Basic.cpu_time(Time);

////////////the name of log file and trajectory file
        stringstream fileLog;
        fileLog<<"LogDyn"<<Type[i];
        string Logname = fileLog.str();

        stringstream Traj;
        Traj << "../Conf/Coord" << Type[i] << ".dcd";
        string Trajname = Traj.str();

//////////////write the log file 

        fstream Log;
        Log.open(Logname,ios::out);Log<<Time<<endl<<endl;

////////////ITO//the number of start point of the calculation
        Basic.Init(&NFrame, &NAtom, &Delt_in[i],&TProd_in[i],&Delt,&IT0,&T0Max,Trajname);

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
        Log<<"  Time   step   of   TCF        :    " <<Delt <<endl;
        Log<<"  Maximum  time  for  TCF       :    "<<Delt*IT0*T0Max<<endl;
        Log<<"  Number of frames in Trajectory:    "<<NFrame<<endl;
        Log<<"  Number of atoms  in Trajectory:    "<<NAtom <<endl<<endl;

        TMax = NFrame;

////////////after how many frames to output
    	Foutput = int(round(double(NFrame)/10));

        Properties Properities1(&NPart,&NBead,&NPoly,&Box,&Vol,&Rho,&TMax,&IT0,&T0Max,&Delt,Log,Type[i]);
        
        Log<<endl<<"---->>calculation begins"<<endl;

        auto input = chemfiles::Trajectory(Trajname,'r');   

        vector<Vector3D> PX(NPoly,Vector3D(0.0,0.0,0.0)); 
        for(int i = 1;i<=NFrame;i++){

            const auto Frame = input.read();

            auto RX1 = Frame.positions();

/////////////calculate the center of mass of each chain
            for (int i=0;i<NPoly;i++){
                for (int j=0;j<NBead;j++){
                    int k=i*NBead+j;
                    PX[i][0] += RX1[k][0];
                    PX[i][1] += RX1[k][1];
                    PX[i][2] += RX1[k][2];
                }
                PX[i][0] = PX[i][0]/NBead;
                PX[i][1] = PX[i][1]/NBead;
                PX[i][2] = PX[i][2]/NBead;
            }

            if(i%10 ==0) {cout<<i<<" "<<NFrame<<endl;}
            
            Properities1.calculate(RX1,PX);

            if(i%Foutput==0){
                Properities1.write();
                Log<<"#frames completed: "<< i <<endl;
            } 
        }
        input.close();
        Properities1.write();

        Log<<"<<----calculation completed"<< endl <<endl;

        TimeFinal = time(NULL);

        Log<<"Total time is "<<(TimeFinal-TimeStart)/60<<" minutes"<<endl;
        Log.close();
    }
    
}
