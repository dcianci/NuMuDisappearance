#ifndef MINIBOONE_DIS_H
#define MINIBOONE_DIS_H

#include "datasets.h"

class MiniBooNE_dis: public dataset{
  public:
    MiniBooNE_dis(bool _nubar){
      nubar = _nubar;
    }
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    bool nubar;
    const int nBins = 16;
  	int nFOsc;

    std::vector < float > FullData, Signal, Libdis_noosc;
    std::vector < std::vector < float > > Full_fractCovMatrix, Libdis_sinsq;

    TMatrixT <float> covMatrix, cov;
};

#endif
