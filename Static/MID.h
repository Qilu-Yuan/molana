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

class MID
{
private:
        /* data */
    vector<double> Acc,Sum,Rn2;
    int NStat;
public:
    void set(int *NBead,fstream& Log);
    void calcluate(int* NBead,int* NPoly,const vector<Vector3D>& RX);
    void write(int *NBead);
};
