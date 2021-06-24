#include "MicroBooNE_dis.h"

int MicroBooNE_dis::Init(std::string dataLoc, Oscillator osc, bool debug){

  ///////////////////////////////////////////////////////////////////////
  shapeonly=false;
  signalInject=false;
  statsonly = false;
  double sin22th_inject(0), dm2_inject(11.0);
  nBins = 19;
  mnu_lowbound = 0.1;  mnu_hibound = 10.0;
  std::string s_variable = "Enu_1m1p";
  std::string s_datatag = "apr22";
  ///////////////////////////////////////////////////////////////////////
  dm2_precalc_density = osc.GridSize();

  Background.resize(nBins);
  FullData.resize(nBins);
  Full_fractCovMatrix.resize(nBins, std::vector<double>(nBins));
  Libdis_sinsq.resize(dm2_precalc_density, std::vector<double>(nBins));
	Libdis_noosc.resize(nBins);
  MCStatSquared_lib.resize(dm2_precalc_density,std::vector<double>(nBins));
  MCStatSquared_noosc.resize(nBins);
  double var_binedges[nBins+1];
  double nu_var[nMC];
  double nu_EnuTrue[nMC];
  double nu_LnuTrue[nMC];
  double nu_weight[nMC];
  double potweight, fullwgt;
  ///////////////////////////////////////////////////////////////////////

  // Load all our data
  ifstream file;
  std::string infile = dataLoc+"uboone/"+s_variable+"_"+s_datatag+"_MC.txt";
  file.open(infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << infile << std::endl;
	for(short i = 0; i < nMC; i++){
		file >> nu_var[i];
    file >> nu_EnuTrue[i];
    file >> nu_LnuTrue[i];
    file >> potweight;
    file >> fullwgt;
    nu_weight[i] = potweight*fullwgt;
  }
	file.close();

  infile = dataLoc+"uboone/"+s_variable+"_"+s_datatag+"_data.txt";
  file.open(infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << infile << std::endl;
  for(short i = 0; i < nBins; i++){
    file >> FullData[i];
  }
  file.close();

  infile = dataLoc+"uboone/"+s_variable+"_"+s_datatag+"_bkg.txt";
  file.open(infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << infile << std::endl;
  for(short i = 0; i < nBins; i++){
    file >> Background[i];
  }
  file.close();

  std::cout << "BINEDGES: ";
  infile = dataLoc+"uboone/"+s_variable+"_"+s_datatag+"_binedges.txt";
  file.open(infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << infile << std::endl;
  for(short i = 0; i <= nBins; i++){
      file >> var_binedges[i];
      std::cout << ", " << var_binedges[i];
  }
  file.close();
  std::cout << std::endl;

  //std::cout << "USING DIFFERENT MATRIX:" << std::endl;
  //infile = dataLoc+"uboone/"+s_variable+"_detsys_fracsysmatrix.txt";
  infile = dataLoc+"uboone/"+s_variable+"_"+s_datatag+"_fracsysmatrix.txt";
	file.open(infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << infile << std::endl;
	for(short i = 0; i < nBins; i++)
		for(short j = 0; j < nBins; j++){
			file >> Full_fractCovMatrix[i][j];
    }
	file.close();

  // To save drastically on the chi2 calculation, precompute all the sines and sine-squareds now!
	for(int mi = 0; mi < dm2_precalc_density; mi++){
    dm2 = pow(IndexToValue(mi,mnu_lowbound,mnu_hibound,dm2_precalc_density),2);

		for(int iB = 0; iB < nBins; iB++){
			Libdis_sinsq[mi][iB] = 0;
			if(mi == 0)
				Libdis_noosc[iB] = 0;
		}
		for(int imc = 0; imc < nMC; imc++){   // Loop over mc
      indi = -1;
      for(int iB = 0; iB < nBins; iB++){    // find which bin our mc fits into
	      if(nu_var[imc] > var_binedges[iB] && nu_var[imc] < var_binedges[iB+1]){
          indi = iB;
          break;
        }
      }
		  double ETru = nu_EnuTrue[imc];
			double LTru = nu_LnuTrue[imc];
			Libdis_sinsq[mi][indi] += nu_weight[imc]*pow(sin(1.267*dm2*LTru/ETru),2); // oscillated events in each bin
      MCStatSquared_lib[mi][indi] += pow(nu_weight[imc]*pow(sin(1.267*dm2*LTru/ETru),2),2);
			if(mi == 0){
				Libdis_noosc[indi] += nu_weight[imc];  // unoscillated events in each bin
        MCStatSquared_noosc[indi] += pow(nu_weight[imc],2);
      }
    }
  }
  dof = nBins;

  // If you want to do a signal injection, I won't stop you. We'll just replace our data.
  if(signalInject){
    std::cout << "INJECTING SIGNAL AT DM2: " << dm2_inject << " SIN22THETA: " << sin22th_inject << std::endl;
    for(int iB = 0; iB < nBins; iB ++){
      std::cout << "DATA: " << FullData[iB] << " ";
      FullData[iB] = Libdis_noosc[iB] + Background[iB];
      std::cout << "NULL: " << FullData[iB] << "  |  " << Libdis_noosc[iB] + Background[iB] << std::endl;
    }
    for(int iB = 0; iB < nBins; iB++){
      dm2i = ValueToIndex(dm2_inject,pow(mnu_lowbound,2),pow(mnu_hibound,2),dm2_precalc_density);
      FullData[iB] -= sin22th_inject*Libdis_sinsq[dm2i][iB];
    }
  }

  // Initialize output tree
  chi2Nt = new OutTree("MicroBooNE");

  if(debug){
    std::cout << "H:data(MC): ";
    for(int iB = 0; iB < nBins; iB++){
      std::cout << FullData[iB] << " ";
    }
    std::cout << std::endl;
  }


  if(debug){
    std::cout << "MicroBooNE initialized. Bins: " << nBins << std::endl;
  }
  return dof;
}

double MicroBooNE_dis::Chi2(Oscillator osc, neutrinoModel model,bool debug, double* data_external){

  if(data_external !=NULL){
    if(debug) std::cout << "Using External Data";
    for(int i = 0; i < nBins; i++){
      FullData[i] = data_external[i];
    }
  }

  double chi2(0.0);

  std::vector <  double > Prediction, PredictionNoOsc, MCStatSquared;
  Prediction.resize(nBins);
  PredictionNoOsc.resize(nBins);
  MCStatSquared.resize(nBins);
  covMatrix.ResizeTo(nBins, nBins);
  covMatrix.Zero();

  // Initialize contributions from osc probability
	sin22th = model.ProbAmp("mumu");

	for(int iB = 0; iB < nBins; iB++){
    dm2i = ValueToIndex(model.Mnu(),mnu_lowbound,mnu_hibound,dm2_precalc_density);
    PredictionNoOsc[iB] = Libdis_noosc[iB] + Background[iB];
    Prediction[iB] = Libdis_noosc[iB] + Background[iB] - sin22th*Libdis_sinsq[dm2i][iB];
    MCStatSquared[iB] = MCStatSquared_noosc[iB] - pow(sin22th,2) * MCStatSquared_lib[dm2i][iB];
	}

  if(shapeonly){
    double obsIntegral(0.0), mcIntegral(0.0), covIntegral(0.0), fnorm;
    for(int iB = 0; iB < nBins; iB++){
      obsIntegral += FullData[iB];
      mcIntegral += Prediction[iB];
      for(int jB = 0; jB < nBins;jB++){
        covIntegral += Full_fractCovMatrix[iB][jB] * Prediction[iB] * Prediction[jB];
      }
    }
    fnorm = covIntegral/pow(mcIntegral,2);
    if(debug) std::cout <<  "FNORM: " << fnorm  << "  ---- " << sqrt(fnorm) << std::endl;

    for(int iB = 0; iB < nBins; iB++){
      Prediction[iB] *= (obsIntegral/mcIntegral); // normalize prediction
    }
    for(int iB = 0; iB < nBins; iB++){
  		for(int jB = 0; jB < nBins; jB++){
        if(statsonly)
  			   covMatrix[iB][jB] = 0;
        else
          covMatrix[iB][jB] = (Full_fractCovMatrix[iB][jB] - fnorm) * Prediction[iB] * Prediction[jB]; // remove normalization component of cov matrix
        if(iB == jB){
          if(MCStatSquared[iB]==0)
            std::cout << "MATRIX WILL BE SINGULAR. PLEEZ FIX. NO EMPTY BINS ALLOWED" << std::endl;
          else
				    covMatrix[iB][jB] += MCStatSquared[iB] + Prediction[iB]; // Add statistical error of signal prediction
  			}
  		}
  	}
  }
  else{
    for(int iB = 0; iB < nBins; iB++){
		  for(int jB = 0; jB < nBins; jB++){
        if(statsonly)
			     covMatrix[iB][jB] = 0;
        else
           covMatrix[iB][jB] = Full_fractCovMatrix[iB][jB] * Prediction[iB] * Prediction[jB];
			  if(iB == jB){
          if(MCStatSquared[iB]==0){
            std::cout << "MATRIX WILL BE SINGULAR. PLEEZ FIX. NO EMPTY BINS ALLOWED" << std::endl;
          }
          else
				    covMatrix[iB][jB] += MCStatSquared[iB] + Prediction[iB]; // Add statistical error of signal prediction
			  }
		  }
	  }
  }

  //std::cout << std::endl;
	// Now, let's invert this bad boy
	cov.ResizeTo(nBins,nBins);

  if(debug){
    std::cout << "H:err: ";
    for(int iB = 0; iB < nBins; iB++){
		  std::cout << sqrt(covMatrix[iB][iB]) << " ";
	  }
    std::cout << std::endl;
  }

	cov = covMatrix.Invert();

	// Finally, let's put everything together and calculate the chisq
	for(int iB = 0; iB < nBins; iB++){
		for(int jB = 0; jB < nBins; jB++){
		   chi2 += (FullData[iB] - Prediction[iB]) * cov[iB][jB] * (FullData[jB] - Prediction[jB]);
		}
	}

  if(debug){
    std::cout << "H:data: ";
    for(int iB = 0; iB < nBins; iB++){
      std::cout << FullData[iB] << " ";
    }
    std::cout << std::endl;
    std::cout << "H:pred: ";
    for(int iB = 0; iB < nBins; iB++){
      std::cout << Prediction[iB] << " ";
    }
    std::cout << std::endl;

    for(int iB = 0; iB < nBins; iB++){
      std::cout << "H:Chi2: ";
  		for(int jB = 0; jB < nBins; jB++){
        std::cout << (FullData[iB] - Prediction[iB]) * cov[iB][jB] * (FullData[jB] - Prediction[jB]) << " ";
      }
      std::cout << std::endl;
    }
  }

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug){
    std::cout << "MicroBooNE Chi2: " << chi2 << std::endl;
  }

  return chi2;
}
