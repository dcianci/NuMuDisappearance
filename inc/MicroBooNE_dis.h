#ifndef MICROBOONE_DIS_H
#define MICROBOONE_DIS_H

#include "datasets.h"

class MicroBooNE_dis{
  public:
    MicroBooNE_dis(){};
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    double Chi2(Oscillator osc, neutrinoModel nu, bool debug,double* data_external=NULL);

    TTree* Tree(){ return chi2Nt->Tree();  };
    void Write(){ chi2Nt->Write(); };


  private:
    bool shapeonly, signalInject,statsonly;
    int nBins, dm2_precalc_density, dm2i, indi;
    const int nMC = 6039;
    double dm2, sin22th, mnu_lowbound, mnu_hibound;

    std::vector < double > FullData, Background, Libdis_noosc,MCStatSquared_noosc;
    std::vector < std::vector < double > > Full_fractCovMatrix, Libdis_sinsq,MCStatSquared_lib;
    TMatrixT <double> covMatrix, cov;

  protected:
    int dof;
    OutTree * chi2Nt;
};

class MicroBooNE_dis_2d{
  public:
    MicroBooNE_dis_2d(){};
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    double Chi2(Oscillator osc, neutrinoModel nu, bool debug);

    TTree* Tree(){ return chi2Nt->Tree();  };
    void Write(){ chi2Nt->Write(); };


  private:
    bool shapeonly, signalInject,statsonly;
    int nBins, dm2_precalc_density, dm2i, indi;
    const int nMC = 6039;
    double dm2, sin22th, mnu_lowbound, mnu_hibound;

    std::vector < double > FullData, Background, Libdis_noosc,MCStatSquared_noosc;
    std::vector < std::vector < double > > Full_fractCovMatrix, Libdis_sinsq,MCStatSquared_lib;

    TMatrixT <double> covMatrix, cov;


  protected:
    int dof;
    OutTree * chi2Nt;
};

#endif
