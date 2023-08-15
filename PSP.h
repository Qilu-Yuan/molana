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

class PSP
{
private:
    /* data */
    int NStat;
    double SumRe2, SumRg2;
    vector<double> Fbb;
    vector<vector<double>> LX,UX;
public:
    void set(int *NBead,fstream& Log);
    void calcluate(int* NBead,int* NPoly,const vector<Vector3D>& RX,const vector<Vector3D>& PX);
    void write(int *NBead,int *NPoly);
};