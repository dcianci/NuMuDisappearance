#ifndef DATASETS_H
#define DATASETS_H

#include "OscTools.h"

class dataset{
  public:
    virtual int Init(std::string dataLoc, Oscillator osc, bool debug){ return 0; };
    virtual float Chi2(Oscillator osc, neutrinoModel nu, bool debug){ return 0.f; };

    TTree* Tree(){ return chi2Nt->Tree();  };
    void Write(){ chi2Nt->Write(); };

  protected:
    int dof;
    OutTree * chi2Nt;
};

#endif
