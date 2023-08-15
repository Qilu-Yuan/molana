#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include<fstream>
#include <chemfiles.hpp>
#include<iomanip>
#include<vector>

#include"RDF.h"
#include"SSF.h"
#include"PSP.h"
#include"ANG.h"
#include"LB2.h"
#include"MID.h"

using namespace std;
using namespace chemfiles;
using std::string;
using std::stringstream;

class Properties
{
private:
    /* data */
    int *NPart, *NBead, *NPoly; //Npart--number of particles, NBead--number of particles on a chain, Npoly-Number of polymers
    double *Box,*Vol; //Box--box size, Vol--volume of the box
    ANG Ang1;RDF RDF1;SSF SSF1;MID MID1;LB2 LB21;PSP PSP1;
public:
    Properties(int *NPart_in,int *NBead_in,int *NPoly_in,double *Box_in,double *Vol_in,fstream& Log);
    void calcluate(const vector<Vector3D>& RX,const vector<Vector3D>& PX);
    void write();
};
