#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <chemfiles.hpp>
#include<iomanip>
#include<fstream>
#include<vector>

using namespace std;
using namespace chemfiles;
using std::string;
using std::stringstream;

class RDF
{
private:
    /* data */
    int NStat,NBin;
    double DelR;
    vector<double> RBin,VBin,Gbb;
public:
    void set(double *Box,fstream& Log);
    void calcluate(double *Box,int *NPart,const vector<Vector3D>& RX);
    void write(double *Vol,int *NPart);
};