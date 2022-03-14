#ifndef EXP_DAT
#define EXP_DAT

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

class experimentalData{
 public:
  experimentalData(const char* scatterFile);
  bool binDataCheck(double &dMax,double &qmin,double &qmax);
  int setPhases(double &dMax,double &kmin,double &kmax);
  double gaussianFixed(double &s);
  double gaussianHydFixed(double &s);
  double hydrationScatter(double &s);
  double fitToScattering(std::vector<double> &molDists,std::vector<double> &solDists,std::vector<double> &solMolDists,int &molSz,int &solSz);
  void writeScatteringToFile(std::vector<double> &molDists,std::vector<double> &solDists,std::vector<double> &solMolDists,int &molSz,int &solSz,const char* filename);
 private:
  std::vector<double> gaussianParams;
  std::vector<double> gaussianParamsHyd;
  std::vector<std::pair<double,double> > exprDat;
  std::vector<double> qvals;
  double absKmin;
  double absKmax;
  int noDistBins;
  std::vector<std::pair<double,double> > distBins;
  std::vector<std::vector<double> > scatPhases;
  double a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,c;
  double a1H,a2H,a3H,a4H,a5H,b1H,b2H,b3H,b4H,b5H,cH;
  std::vector<double> experimentalIntensity;
  std::vector<double> experimentalSD;
  std::vector<std::pair<double,double> > scatVec;
};

#endif


