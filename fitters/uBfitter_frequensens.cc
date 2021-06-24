// May 13, 2021
// Here, we'll do some fake data studies. The only differences from before are that we're looping through fake data, fed in by an external source (see infile below)

// This code scans through a grid of points in 3+1 sterile neutrino oscillation space (bounds defined below) and calculates the chi2 for numu disappearance in microboone at each point
// inputs are "-n X" 	which assigns a grid density X by X (default is 500x500)
//						"-d" 		debug mode
//						"-s"		single point (rather than grid; good for testing)


#include "fitter.h"
#include "MicroBooNE_dis.h"

bool debug,singlept;
int ngrdpts;

int bruteforce(){

  MicroBooNE_dis* ub = new MicroBooNE_dis;
  const int nExp = 1001;
  const int nBins = 19;

  Oscillator osc;
  osc.gridpts = ngrdpts;

  std::cout << "Initializing"  << std::endl;
  std::string dataLoc = "/home/dcianci/Physics/MicroBooNEDisappearance/SBN_3plusN/GlobalFits_v2/data/";
  ub->Init(dataLoc,osc,debug);
  std::cout << "Dataset Initialized!" << std::endl;

  // Load in fake data (created elsewhere with sbnfit)
  std::vector<std::vector<double>> v_fakeData;
  v_fakeData.resize(nExp,std::vector<double>(nBins));

  std::ifstream file;
  //std::string infile = dataLoc+"uboone/fakedatastore_june11_highdm2.txt";
  std::string infile = dataLoc+"../build/fakedatastore_s_0.350000_d_2.000000.txt";
  file.open(infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << infile << std::endl;
	for(short i = 0; i < nExp; i++){
    for(short j = 0; j < nBins; j++){
      file >> v_fakeData[i][j];
    }
  }
	file.close();

  // Create output File
  std::string s_outfile = "fakedata_chi2.txt";
	std::cout << "Output File: " << s_outfile << std::endl;
  std::ofstream myfile;
  myfile.open(s_outfile);


  // Create a nu model
  neutrinoModel nuModel;

  int count = 0;
  double chi2;
  std::cout << "Beginning chi2 loop" << std::endl;

  double theta_lowbound(.01), theta_hibound(3.1415926/4.);
  double mnu_lowbound(.1), mnu_hibound(10.0);
  double sin22th_lowbound(0.01), sin22th_hibound(1.0);
  double* _fakedata = new double[nBins];

  for(int iexp = 0; iexp < nExp; iexp++){
    std::cout << "EXP: " << iexp << std::endl;
    for(int ib = 0; ib < nBins; ib++){
      _fakedata[ib] = v_fakeData[iexp][ib];
    }

    // first, calculate  chi2 at point.
    double si_sin22th = 0.4;
    double si_theta24 = asin(sqrt((1.0-sqrt(1.0-si_sin22th))/2.0));
    double si_mnu = sqrt(70.0);
    nuModel.Init(si_mnu,0,si_theta24,0);
    // Calculate chi2s
    double si_chi2 = ub->Chi2(osc,nuModel,debug,_fakedata);
    double chi2min = 9999;

    count = 0;
    for(int mi = 0; mi < ngrdpts; mi++){
      for(int t42i = 0; t42i < ngrdpts; t42i++){
        std::cout << "Progress: " << float(count)/(pow(ngrdpts,2)/(100.f)) << "\% \r";

        double mnu = IndexToValue(mi,mnu_lowbound,mnu_hibound,ngrdpts);
        double theta14 = 0;
        double theta34 = 0;
        double sin22th = IndexToValue(t42i,sin22th_lowbound,sin22th_hibound,ngrdpts);
        double theta24 = asin(sqrt((1.0-sqrt(1.0-sin22th))/2.0));
        nuModel.Init(mnu,theta14,theta24,theta34);
        // Calculate chi2s
        chi2 = ub->Chi2(osc,nuModel,debug,_fakedata);
        myfile << chi2 << " ";
        //chi2min = min(chi2,chi2min);

        count++;
      }
    }
    //myfile << si_chi2-chi2min << std::endl;
    std::cout << std::endl;
    myfile << std::endl;
  }

  myfile.close();


  std::string s_binsfile = "fakedata_coords.txt";
  std::ofstream binsfile;
  binsfile.open(s_binsfile);
  for(int mi = 0; mi < ngrdpts; mi++){
    for(int t42i = 0; t42i < ngrdpts; t42i++){
      double mnu = IndexToValue(mi,mnu_lowbound,mnu_hibound,ngrdpts);
      double sin22th = IndexToValue(t42i,sin22th_lowbound,sin22th_hibound,ngrdpts);
      double theta24 = asin(sqrt((1.0-sqrt(1.0-sin22th))/2.0));
      binsfile << sin22th << " " << pow(mnu,2) << std::endl;
    }
  }
  binsfile.close();

  return 0;
}

int main(int argc, char* argv[]){

  int iarg(0), index;
  opterr=1;

  ngrdpts = 500;
  debug = false;
  singlept = false;

  const struct option longopts[] = {
    {"debug", optional_argument, 0, 'd'},
    {"ngrid", optional_argument, 0, 'n'},
    {"singlept", optional_argument, 0, 's'},
	  {0,			no_argument, 		0,  0},
  };

  while(iarg != -1){
    iarg = getopt_long(argc,argv, "n:ds", longopts, &index);

    switch(iarg){
		  case 'n':
			  ngrdpts = atoi(optarg);
        std::cout << "GRIDSIZE: " << ngrdpts << std::endl;
			  break;
      case 'd':
        debug = true;
        std::cout << "DEBUG MODE" << std::endl;
        break;
      case 's':
        singlept = true;
        std::cout << "SINGLE POINT MODE" << std::endl;
        break;
      case '?':
		  case 'h':
			  std::cout<<"I'd like an input, friend."<<std::endl;
        std::cout<<"\t-n\t--ngrid\t\t Input int grid width (default 500)"<<std::endl;
        std::cout<<"\t-s\t--singlept"<<std::endl;
        std::cout<<"\t-d\t--debug"<<std::endl;
			  return 0;
	  }
  }

  bruteforce();
  return 0;
}
