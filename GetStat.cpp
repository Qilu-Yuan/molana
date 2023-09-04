#include"GetStat.h"

int main(){
    int NPart, NBead, NPoly; //Npart--number of particles, NBead--number of particles on a chain, Npoly-Number of polymers
    double Box,Vol,Rho,Temp,Pres;
	int NFrame,NAtom,NItem,NStep, FStatic, Foutput;
    int Nevery = 5,Nrepeat =1000,Nfreq=5000;
	time_t TimeStart, TimeFinal;
    string Time;

    // int Nevery = 2,Nrepeat =6,Nfreq=100;
    if (Nevery <= 0 || Nrepeat <= 0 || Nfreq <= 0 || Nrepeat*Nevery>Nfreq){
        cout<<"Illegal Parameter!"<<endl;
        return 1;
    }  
    TimeStart = time(NULL);

    cpu_time(Time);
    
    fstream Log;

    Log.open("LogStat",ios::out);Log<<Time<<endl<<endl;

    
    Init(&NPart,&NBead,&NPoly,&Box,&Vol,&Rho,&NItem,&NStep,&Temp,&Pres,&NFrame, &NAtom);

    int Ncycle = NFrame/Nfreq;

    int NTRE = Ncycle*(Nrepeat); //the total repeat times

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
    Log<<"  Number of atoms  in Trajectory:    "<<NAtom <<endl;
    Log<<"  Number of frames in Trajectory:    "<<NFrame<<endl;
    Log<<"  Number of frames to Be Read   :    "<<NTRE<<endl<<endl;
	Log<<"  Data info in Thermo           :    "<<NItem<<" , "<< NStep<<endl;
	Log<<"  Temperature                   :    "<<Temp<<endl;
	Log<<"  Pressure                      :    "<<Pres<<endl<<endl;

    Properties Properties1(&NPart,&NBead,&NPoly,&Box,&Vol,Log);

    Log<<endl<<"---->>calculation begins"<<endl;


    auto input = chemfiles::Trajectory("../Conf/CoordL.dcd",'r');

    vector<Vector3D> PX(NPoly, Vector3D(0.0, 0.0, 0.0));
    vector<vector<Vector3D>> RX(NTRE, vector<Vector3D>(NPart, Vector3D(0.0, 0.0, 0.0)));

///////////read the frames to be calculated////////////
    int r = 0;
    cout << "Frame ";
    for(int i = 1;i<=NFrame&&r<NTRE;i++){
        const auto Frame = input.read();
        
        int Nrepeatleft;
        if(i%Nfreq ==1) Nrepeatleft = Nrepeat;

        if((i+(Nevery*(Nrepeatleft-1)))%Nfreq==0 && Nrepeatleft>0){
            RX[r] = Frame.positions();
            r++;
            Nrepeatleft--;
            cout << i <<" ";
        }
    }
    cout<< " are read." << endl;

    input.close();
/////////////////////////////////////////////////////
    Foutput = int(round(double(NTRE)/10));
    for(int i = 0;i<NTRE;i++){
        for (int j=0;j<NPoly;j++){
            PX[j][0] = 0;
            PX[j][1] = 0;
            PX[j][2] = 0;
            for (int k=0;k<NBead;k++){
                int l = j*NBead+k;
                PX[j][0] += RX[i][l][0];
                PX[j][1] += RX[i][l][1];
                PX[j][2] += RX[i][l][2];
            }
            PX[j][0] = PX[j][0]/NBead;
            PX[j][1] = PX[j][1]/NBead;
            PX[j][2] = PX[j][2]/NBead;
        }

        Properties1.calcluate(RX[i],PX);
        cout<<i<<" "<<NTRE<<endl;    
        if((i+1)%Foutput==0){
            Properties1.write();
            Log<<"#frames completed: "<< i+1 <<endl;
        } 
    }

    Properties1.write();

    Log<<"<<----calculation completed"<< endl <<endl;

    TimeFinal = time(NULL);

    Log<<"Total time is "<<double(TimeFinal-TimeStart)/60.0<<" minutes"<<endl;
    Log.close();
}
