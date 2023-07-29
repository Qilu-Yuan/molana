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

using namespace std;
using namespace chemfiles;
using std::string;
using std::stringstream;

void cpu_time(string& Time);
void Init(int* NPart, int* NBead, int* NPoly, double* Box, double* Vol, double* Rho, int* NItem, int* NStep, double* Temp, double* Pres, int* NFrame, int* NAtom);
