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

  Oscillator osc;
  osc.gridpts = ngrdpts;

  std::cout << "Initializing"  << std::endl;
  std::string dataLoc = "/home/dcianci/Physics/MicroBooNEDisappearance/SBN_3plusN/GlobalFits_v2/data/";
  ub->Init(dataLoc,osc,debug);
  std::cout << "Dataset Initialized!" << std::endl;

  // Create output File
  std::string outfile = "ubfit.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  OutTree * chi2Nt = new OutTree("Total");

  // Create a nu model
  neutrinoModel nuModel;

  int count = 0;
  double chi2;
  std::cout << "Beginning chi2 loop" << std::endl;

  double theta_lowbound(.01), theta_hibound(3.1415926/4.);
  double mnu_lowbound(.1), mnu_hibound(10.0);
  double sin22th_lowbound(0.01), sin22th_hibound(1.0);

  double chi2min(9999.9), theta24min, mnumin, cmin;

  for(int mi = 0; mi < ngrdpts; mi++){
    for(int t42i = 0; t42i < ngrdpts; t42i++){
      if(!debug)
      	std::cout << "Progress: " << float(count)/(pow(ngrdpts,2)/(100.f)) << "\% \r";

      double mnu = IndexToValue(mi,mnu_lowbound,mnu_hibound,ngrdpts);
      double theta14 = 0;
      double theta34 = 0;
      double sin22th = IndexToValue(t42i,sin22th_lowbound,sin22th_hibound,ngrdpts);
      double theta24 = asin(sqrt((1.0-sqrt(1.0-sin22th))/2.0));
      nuModel.Init(mnu,theta14,theta24,theta34);

      // Calculate chi2s
      chi2 = ub->Chi2(osc,nuModel,debug);

      if(chi2 < chi2min){
        chi2min = chi2;
        theta24min = theta24;
        mnumin = mnu;
        cmin = count;
      }

      if(debug){
        //std::cout << "Total chi2 for model: " << mnu  << " " <<  theta14 << " " << theta24 << " " <<  theta34 << std::endl;
        std::cout << chi2 << std::endl;
      }

      chi2Nt->Fill(chi2,0,nuModel);
      count ++;

      if(singlept)
        break;
    }
    if(singlept)
      break;
  }

  std::cout << "MINS: " << theta24min << " " << mnumin << " " << chi2min << " " << cmin << std::endl;

  // Write everything to File
  std::cout << "Writing to file." << std::endl;
  chi2Nt->Write();
  ub->Write();
  f->Close();

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
