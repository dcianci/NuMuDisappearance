#ifndef OSCTOOLS_H
#define OSCTOOLS_H
#include "globalFit.h"

#define dm2VecMaxDim 601

class neutrinoModel{
private:
  double mnu, theta14, theta24, theta34;
  int Zero(){
    mnu = 0.0;  theta14 = 0.0;  theta24 = 0.0;  theta34 = 0.0;
  };

public:
  neutrinoModel(){
    Zero();
  };
  int Init(double _mnu, double _par1, double _par2, double _par3 = 0.0){
    mnu = _mnu; theta14 = _par1;  theta24 = _par2; theta34 = _par3;
    return 1;
  };

  double ProbAmp(std::string interaction){
    if(interaction == "ee"){
      return pow(sin(2*theta14),2);
    }
    else if(interaction == "mumu" || interaction == "mm"){
      return 4*pow(cos(theta14),2)*pow(sin(theta24),2)*(1-pow(cos(theta14),2)*pow(sin(theta24),2));
    }
    else if(interaction == "mue" || interaction == "me"){
      return pow(sin(2*theta14),2)*pow(sin(2*theta24),2);
    }
    else{
      std::cout << "BAD PROBAMP" << std::endl;
      return -1;
    }
  };

  std::array<double,4> OscParams(){
    std::array<double,4> ops = { mnu, theta14, theta24, theta34};
    return ops;
  };

  double Dm2(){ return pow(mnu,2); };
  double Umu4(){ return cos(theta14)*sin(theta24); };
  double Mnu(){ return mnu; };

  void Print(){
    std::cout << "MNU: " << mnu << std::endl;
    std::cout << "TH14: " << theta14 << std::endl;
    std::cout << "TH24: " << theta24 << std::endl;
    std::cout << "TH34: " << theta34 << std::endl;
  }
};


struct neutrinoModel_general{
    double mNu[3], Ue[3], Um[3], phi[3];
    double dm41Sq, dm51Sq, dm61Sq, dm54Sq, dm64Sq, dm65Sq;
    void zero(){
        for(int i = 0; i < 3; i ++){
            mNu[i] = 0; Ue[i] = 0;
            Um[i] = 0;  phi[i] = 0;
        }
    }
    void difference(){
        dm41Sq = pow(mNu[0],2);
        dm51Sq = pow(mNu[1],2);
        dm61Sq = pow(mNu[2],2);
        dm54Sq = dm51Sq - dm41Sq;
        dm64Sq = dm61Sq - dm41Sq;
        dm65Sq = dm61Sq - dm51Sq;
    }
};

struct oscContribution{
    double dm2[6], aMuE[6], aMuMu[6], aMuE_CPV[6], aEE[6];
};

class Oscillator{
  public:
    Oscillator(){};
    Oscillator(float _dm2Min, float _dm2Max, float _UMin, float _UMax, float _USqMax, float _stepSize, float _temperature, int _nSteriles, int _gridpts, bool _CPConserving, int _nmcgen, int seed);
    neutrinoModel InitializeMarkovParams();
    neutrinoModel NewModel(neutrinoModel modelOld);
    bool RejectModel(neutrinoModel model);

    int GridSize(){ return gridpts; };
    double dm2Vec[dm2VecMaxDim];

    TRandom3 RanGen;

    bool UsingAtm;

    int nMCGen, nSteriles;

    void PrintMarkovSettings(){
      std::cout << "Temperature: " << temp << " Stepsize: " << step << std::endl;
    }
		float temp,ran[13];
    int gridpts;

  private:

    float dm2Min, dm2Max, UMin, UMax, USqMax, step;
    bool CPConserving, reject1, reject2, reject3, reject4, usingUe, usingUm;
};


class OutTree{
  public:
    OutTree(std::string tag);

    void Fill(float _chi2, int _dof, neutrinoModel _nuModel);
    void Write(){ myTree->Write();  };
    TTree *Tree(){ return myTree->CloneTree(); };

  private:
    TTree *myTree;
    int dof;
    double chi2, m41, theta14, theta24, theta34;
};

// Common integral functions
double sinFunc(double x);
double sineInt(double x);
double cosFunc(double x);
double cosineInt(double x);
double IndexToValue(double _index, double _min, double _max, int _grdpts, std::string _scale="log");
int ValueToIndex(double value, double _min, double _max, int _grdpts, std::string _scale="log");

#endif
