#include "OscTools.h"

OutTree::OutTree(std::string tag){
  TString t_tag = tag;
  myTree = new TTree(t_tag,t_tag);
  myTree->Branch("chi2",&chi2,"chi2/D");
	myTree->Branch("dof",&dof,"dof/I");
  myTree->Branch("m41",&m41,"m41/D");
  myTree->Branch("theta14",&theta14,"theta14/D");
	myTree->Branch("theta24",&theta24,"theta24/D");
  myTree->Branch("theta34",&theta34,"theta34/D");
}

void OutTree::Fill(float _chi2, int _dof, neutrinoModel _nuModel){
  chi2 = _chi2; dof = _dof;
  std::array<double,4> ops = _nuModel.OscParams();
  m41 = ops[0];  theta14 = ops[1]; theta24 = ops[2]; theta34 = ops[3];
  myTree->Fill();
}

Oscillator::Oscillator(float _dm2Min, float _dm2Max, float _UMin, float _UMax, float _USqMax, float _stepSize, float _temperature, int _nSteriles, int _gridpts, bool _CPConserving, int _nmcgen, int seed){
  dm2Min = _dm2Min;   dm2Max = _dm2Max;
  UMin = _UMin;       UMax = _UMax;
  gridpts = _gridpts;
  USqMax = _USqMax;
  nSteriles = _nSteriles;
  nMCGen = _nmcgen;
  CPConserving = _CPConserving;

  usingUe = true;
  usingUm = true;

  RanGen.SetSeed(seed);
  step = RanGen.Rndm() * _stepSize;
	temp = RanGen.Rndm() * _temperature;

  for(int i = 0; i < dm2VecMaxDim; i++){
      dm2Vec[i] = pow(10,TMath::Log10(dm2Min) + double(i) / (dm2VecMaxDim-1) * TMath::Log10(dm2Max/dm2Min));
  }
}

/*
// Get your models sorted out (all of this has been tested to shit. we're good.)
neutrinoModel Oscillator::InitializeMarkovParams(){

  // Initialize new model!
  neutrinoModel modelOld;
  modelOld.zero();
  bool reject;
  int nCPFactors = nSteriles*(nSteriles-1)/2;

  do{
    // Fill up our random array
    RanGen.RndmArray(13,ran);

    // Initial Params for Mass and Mixing
    if(nSteriles > 0){
      for(int i = 0; i < nSteriles; i++){
        modelOld.mNu[i] = pow(10., (TMath::Log10(dm2Min) + ran[3*i]*TMath::Log10(dm2Max/dm2Min))/2);
        //modelOld.mNu[i] = pow(10., (TMath::Log10(.317) + ran[3*i]*TMath::Log10(3.832/.316))/2);
        modelOld.Ue[i] = pow(10., TMath::Log10(UMin) + ran[3*i + 1] * TMath::Log10(UMax/UMin));
        //modelOld.Ue[i] = pow(10., TMath::Log10(.051) + ran[3*i + 1] * TMath::Log10(.57/0.051));
        modelOld.Um[i] = pow(10., TMath::Log10(UMin) + ran[3*i + 2] * TMath::Log10(UMax/UMin));
      }
    }
    // Now, let's do the CP factors, phi
    if(nCPFactors > 0){
      for(int i = 0; i < nCPFactors; i++){
        modelOld.phi[i] = double(ran[nSteriles*3+i]*2*TMath::Pi());
        if(CPConserving){
          if(modelOld.phi[i] < TMath::Pi()) modelOld.phi[i] = 0;
          else  modelOld.phi[i] = TMath::Pi();
        }
      }
    }
    reject = RejectModel(modelOld);
  }while(reject);

  return modelOld;
}

neutrinoModel Oscillator::NewModel(neutrinoModel modelOld){

  // Initialize new model!
  neutrinoModel model;
  model.zero();
  bool reject;
  int nCPFactors = nSteriles*(nSteriles-1)/2;

  do{
    // Generate some random numbers!
    RanGen.RndmArray(13,ran);

    // Alright, let's step forward with these masses and mixing matrix elements!
    for(int i = 0; i < nSteriles; i++){
			model.mNu[i] = pow(10., (TMath::Log10(modelOld.mNu[i]) + (ran[3*i] - .5)*2*step*TMath::Log10(dm2Max/dm2Min))/2);
      //model.mNu[i] = pow(10., (TMath::Log10(modelOld.mNu[i]) + (ran[3*i] - .5)*2*step*TMath::Log10(3.832/.316))/2);
      //if(usingUe) model.Ue[i] = pow(10.,TMath::Log10(modelOld.Ue[i]) + (ran[3*i+1] - .5)*2*step*TMath::Log10(.57/0.051));
      if(usingUe) model.Ue[i] = pow(10.,TMath::Log10(modelOld.Ue[i]) + (ran[3*i+1] - .5)*2*step*TMath::Log10(UMax/UMin));
      else    model.Ue[i] = 0.;

      if(usingUm) model.Um[i] = pow(10.,TMath::Log10(modelOld.Um[i]) + (ran[3*i+2] - .5)*2*step*TMath::Log10(UMax/UMin));
      else    model.Um[i] = 0.;
    }
    if(nCPFactors > 0){
      for(int j = 0; j < nCPFactors; j++){
        model.phi[j] = modelOld.phi[j] + 2.*(ran[nSteriles*3 + j] - 0.5)*2.*TMath::Pi()*step;
        if(CPConserving == 1){
          if(model.phi[j] < TMath::Pi())    model.phi[j] = 0;
          else model.phi[j] = TMath::Pi();
        }
      }
    }
    reject = RejectModel(model);
		// for prospect, make sure ue4 is within proper bounds
		//reject = reject || model.Ue[0] < .051 || model.Ue[0] > .57;
  }while(reject);

  return model;
}

bool Oscillator::RejectModel(neutrinoModel model){

  int nCPFactors = nSteriles*(nSteriles-1)/2;
  // Now, we'll reject the model if matrix elements are too large
  reject1 = pow(model.Ue[0],2) + pow(model.Um[0],2) > USqMax || pow(model.Ue[1],2) + pow(model.Um[1],2) > USqMax || pow(model.Ue[2],2) + pow(model.Um[2],2) > USqMax || pow(model.Ue[0],2) + pow(model.Ue[1],2) + pow(model.Ue[2],2) > USqMax || pow(model.Um[0],2) + pow(model.Um[1],2) + pow(model.Um[2],2) > USqMax;

  // Another condition can be applied to avoid a negative under the square root for atmospheric neutrinos
  if(UsingAtm){
    double dmuMax = .25;

    double A = (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2)) * (pow(model.Um[0],2) + pow(model.Um[1],2) + pow(model.Um[2],2)) +
          pow((model.Um[0]*model.Um[1]),2) + pow((model.Um[0]*model.Um[2]),2) + pow((model.Um[1]*model.Um[2]),2);

    reject1 = reject1 || (1. - 4*A) < pow(1-2*dmuMax,2);
  }

  // More rejection conditions!
  reject2 = false;
  reject3 = false;
  reject4 = false;
  if(nSteriles > 1){
    // DEGENERACY FIX
    reject2 = model.mNu[1] < model.mNu[0] || abs(pow(model.mNu[1],2) - pow(model.mNu[0],2)) < dm2Min;
  }
  if(nSteriles > 2){
    // DEGENERACY FIX
    reject3 = model.mNu[2] < model.mNu[0] || model.mNu[2] < model.mNu[1] || abs(pow(model.mNu[2],2) - pow(model.mNu[0],2)) < dm2Min || abs(pow(model.mNu[2],2) - pow(model.mNu[1],2)) < dm2Min;
  }

  // For the Markov chain case, gotta check a few more things. Essentially whether or not we've stepped out of bounds.
  if(nSteriles > 0){
    for(int i = 0; i < nSteriles; i++){
      reject4 = reject4 || pow(model.mNu[i],2) < dm2Min || pow(model.mNu[i],2) > dm2Max;

      if(usingUe)  reject4 = reject4 || model.Ue[i] < UMin || model.Ue[i] > UMax;
      if(usingUm)  reject4 = reject4 || model.Um[i] < UMin || model.Um[i] > UMax;
    }
  }
  if(nCPFactors > 0){
    for(int i = 0; i < nCPFactors; i++){
      reject4 = reject4 || model.phi[i] < 0 || model.phi[i] > 2*TMath::Pi();
    }
  }

  return (reject1 || reject2 || reject3 || reject4);
}
*/

double sinFunc(double x){
    // Sine function for integration
    return sin(x)/x;
}
double sineInt(double x){
    // Sine function for integration reasons
    ROOT::Math::Functor1D wf(&sinFunc);
    ROOT::Math::Integrator ig;
    ig.SetFunction(wf);

    return ig.Integral(0.,x);
}
double cosFunc(double x){
    // cosine function for integration
    return (cos(x) - 1)/x;
}
double cosineInt(double x){
    // Sine function for integration reasons
    ROOT::Math::Functor1D wf(&cosFunc);
    ROOT::Math::Integrator ig;
    ig.SetFunction(wf);

    return TMath::EulerGamma() + log(x) + ig.Integral(0.,x);
}
double IndexToValue(double _index, double _min, double _max, int _grdpts, std::string _scale){
  if(_scale == "log"){
    return pow(10,log10(_min) + _index * log10(_max/_min)/(_grdpts));
  }
  else if(_scale == "linear"){
    return _min + _index*(_max-_min)/(_grdpts);
  }
  else{
    std::cout << "Scale not yet supported" << std::endl;
    return -999;
  }
}

int ValueToIndex(double value, double _min, double _max, int _grdpts, std::string _scale){
  int indi;
  if(_scale == "log"){
    indi = floor((_grdpts)*log10(value/_min)/log10(_max/_min));
  }
  else if(_scale == "linear"){
    indi = floor((_grdpts)*(value-_min)/(_max-_min));
  }
  else{
    std::cout << "Scale not yet supported" << std::endl;
    return -999;
  }
  if(indi < 0 || indi > _grdpts)
    return -1;
  else
    return indi;
}
