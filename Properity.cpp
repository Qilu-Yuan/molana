#include "Properity.h"

Properties::Properties(int *NPart_in,int *NBead_in,int *NPoly_in,double *Box_in,double *Vol_in,fstream& Log)
{

    NPart = NPart_in;
    NBead = NBead_in;
    NPoly = NPoly_in;
    Box = Box_in;
    Vol = Vol_in;

//set some initial value of the properties
    Ang1.set(NBead,NPoly,Log);
    RDF1.set(Box,Log);
    SSF1.set(Box,Log);
    MID1.set(NBead,Log);
    LB21.set(NBead,NPoly,Log);
    PSP1.set(NBead,Log);

}
void Properties::calcluate(const vector<Vector3D>& RX,const vector<Vector3D>& PX){
    Ang1.calcluate(NBead,NPoly,RX);
    RDF1.calcluate(Box,NPart,RX);
    SSF1.calcluate(NPart,RX);
    MID1.calcluate(NBead,NPoly,RX);
    LB21.calcluate(NBead,NPoly,RX);
    PSP1.calcluate(NBead,NPoly,RX,PX);
}
void Properties::write(){
    Ang1.write(NBead,NPoly);
    RDF1.write(Vol,NPart);
    SSF1.write(NPart);
    MID1.write(NBead);
    LB21.write(NBead);
    PSP1.write(NBead,NPoly);
}
