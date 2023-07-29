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

void LB2(int Switch,int *NPart,int *NBead,int *NPoly,double *Box,double *Vol,vector<vector<double>>& RX,vector<vector<double>>& PX);
