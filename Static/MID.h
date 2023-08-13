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

void MID(int Switch,int *NPart,int *NBead,int *NPoly,double *Box,double *Vol,const vector<Vector3D>& RX,const vector<Vector3D>& PX);
