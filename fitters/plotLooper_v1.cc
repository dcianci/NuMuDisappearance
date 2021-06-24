// This code does everything you need to easily make sensitivity plots. It will:
// Marginalize your 3d points across into three different 2d osc channels
// Perform raster scans for three different osc channels
// Create a chiogram for ""
// Enjoy :)

#include "fitter.h"
#include "TH3D.h"

int ntupleProcess(std::string xml){

  bool debug = true;
  bool useIC = false;
  bool use3dplots = false;
  bool use2dplots = false;

  ProcessReader rdr;
  if(rdr.Load(xml))
    return 0;

  double m41_hibound(10.0), m41_lowbound(.1f);
  double t14_hibound(3.14/4.0), t14_lowbound(.01);
  double t24_hibound(3.14/4.0), t24_lowbound(.01);
  double t34_hibound(3.14/4.0), t34_lowbound(.01);
  double sin22th_hibound(1.0), sin22th_lowbound(0.01);

  // Now, sum up the chi2's of everything and throw them in a vector
  int dof;
  std::vector < double > v_chi2, v_m41, v_theta14, v_theta24, v_theta34;

  // Initialize input tree variables
  double _chi2, _m41, _theta14, _theta24, _theta34;
  int _dof(0);

  // Initialize output tree variables
  double chi2_min(9999.0), sin22th_mue_min, sin22th_ee_min,sin22th_mumu_min, dm2_min;
  double chi2, dm2, sin22th, sin22th_mue, sin22th_mumu, sin22th_ee;
  double m41_min, theta14_min,  theta24_min;

  // For raster scans,  i'm going to throw everything on a grid and just scan from left to right.
  std::vector < std::vector < double > > chi2grid_mue,chi2grid_mumu,chi2grid_ee;
  chi2grid_mue.assign(rdr.gridpts_sin22th, std::vector < double > (rdr.gridpts_dm2, 0.));
  chi2grid_mumu.assign(rdr.gridpts_sin22th, std::vector < double > (rdr.gridpts_dm2, 0.));
  chi2grid_ee.assign(rdr.gridpts_sin22th, std::vector < double > (rdr.gridpts_dm2, 0.));


  std::cout << "Combining the chi2s from multiple datasets into one." << std::endl;
  std::cout << "If any trees are 3d, please have those listed FIRST in your XML" << std::endl;
  std::cout << "NTrees: " << rdr.GetNTrees() << std::endl;

  // create vector containers
  v_chi2.resize(rdr.GetTree(0)->GetEntries());
  std::fill(v_chi2.begin(), v_chi2.end(), 0);

  v_m41.resize(rdr.GetTree(0)->GetEntries());
  v_theta14.resize(rdr.GetTree(0)->GetEntries());
  v_theta24.resize(rdr.GetTree(0)->GetEntries());
  v_theta34.resize(rdr.GetTree(0)->GetEntries());

  for(int i = 0; i < rdr.GetNTrees(); i++){
    if(rdr.GetName(i) == "IceCube" && rdr.GetNTrees() != 1){  // don't do this if we're looking at ONLY icecube
      useIC = true;
      continue;
    }

    rdr.GetTree(i)->SetBranchAddress("chi2",&_chi2);
    rdr.GetTree(i)->SetBranchAddress("dof",&_dof);
    rdr.GetTree(i)->SetBranchAddress("m41",&_m41);
    rdr.GetTree(i)->SetBranchAddress("theta14",&_theta14);
    rdr.GetTree(i)->SetBranchAddress("theta24",&_theta24);
    rdr.GetTree(i)->GetEntry(10); // this does nothing important but it makes root happy.
    use3dplots = use3dplots || rdr.GetOscType(i);
    use2dplots = !(use2dplots || rdr.GetOscType(i));

    dof += _dof;

    for(int j = 0; j < rdr.GetTree(i)->GetEntries(); j++){
      rdr.GetTree(i)->GetEntry(j);

      if(rdr.GetOscType(i)>0){
        v_chi2[j] += _chi2;
        if(i==0){
          v_m41[j] = _m41;
          v_theta14[j] = _theta14;
          v_theta24[j] = _theta24;
        }
        else{
          if(v_m41[j] != _m41 || v_theta14[j] != _theta14 || v_theta24[j] != _theta24){
            std::cout << "MISALIGNMENT IN 3D TREE " << i << " at  entry: " << j << std::endl;
            std::cout << v_m41[j] << " / " << _m41 << " " << v_theta14[j] << " / " << _theta14 << " " << v_theta24[j] << " / " << _theta24<< std::endl;
          }
        }
      }
      else if(use3dplots){ // if we've got a 2d case (nue dis)
        for(int x = 0; x < pow(v_chi2.size(),1.0/3.0); x++){
          v_chi2[j*100+x] += _chi2;

          if(v_m41[j*100] != _m41 || v_theta14[j*100] != _theta14){
            std::cout << "MISALIGNMENT IN 2D TREE " << i << " at  entry: " << j << std::endl;
            std::cout << v_m41[j*100] << " / " << _m41 << " " << v_theta14[j*100] << " / " << _theta14 << std::endl;
          }
        }
      }
      else{ // if we only have 2d cases (nue dis)
          v_chi2[j] += _chi2;
          if(i==0){
            v_m41[j] = _m41;
            v_theta14[j] = _theta14;
            v_theta24[j] = _theta24;
          }
          else{
              if(v_m41[j] != _m41 || v_theta14[j] != _theta14){
            std::cout << "MISALIGNMENT IN 2D TREE " << i << " at  entry: " << j << std::endl;
            std::cout << v_m41[j] << " / " << _m41 << " " << v_theta14[j] << " / " << _theta14 << std::endl;
          }
        }
      }
    }
  }

  //  Now we  worry about icecube. only a special case because it's got 4 dimensions
  if(useIC) {
    std::cout << "Load icecube, the special case for it's fourth  dimension" << std::endl;
    TFile* fProcIC = new TFile("/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/data/IC_chi2grids.root","READ");
    std::vector<TH3D*> vecLogICHisto;
    for(int i = 0; i < 10; i++){
      vecLogICHisto.push_back((TH3D*)fProcIC->Get(("loghisto_"+std::to_string(i)).c_str()));
    }

    int size1 = v_chi2.size();

    v_chi2.resize(size1*10);
    v_m41.resize(size1*10);
    v_theta14.resize(size1*10);
    v_theta24.resize(size1*10);
    v_theta34.resize(size1*10);

    for(int i =  0; i < size1;i++){
      for(int ith34 = 1; ith34 < 10; ith34++){
        v_chi2.at(ith34*size1+i)= v_chi2.at(i)+vecLogICHisto.at(ith34)->Interpolate(log10(v_m41.at(i)),log10(v_theta14.at(i)),log10(v_theta24.at(i)));
        v_m41.at(ith34*size1+i) = v_m41.at(i);
        v_theta14.at(ith34*size1+i)= v_theta14.at(i);
        v_theta24.at(ith34*size1+i)= v_theta24.at(i);
        v_theta34.at(ith34*size1+i)= ith34;
      }
      v_chi2[i] += vecLogICHisto.at(0)->Interpolate(log10(v_m41.at(i)),log10(v_theta14.at(i)),log10(v_theta24.at(i)));
    }
    std::cout << "New size = " << v_chi2.size();
  }

  int imin;
  for(int i = 0; i < v_chi2.size(); i++){
    dm2 = pow(v_m41[i],2);
    sin22th_mue = pow(sin(2*v_theta14[i]),2)*pow(sin(2*v_theta24[i]),2);
    sin22th_mumu = 4 * pow(cos(v_theta14[i]),2) * pow(sin(v_theta24[i]),2) * (1 - pow(cos(v_theta14[i]),2) * pow(sin(v_theta24[i]),2));
    sin22th_ee = pow(sin(2*v_theta14[i]),2);

    if(v_chi2[i] < chi2_min){
      chi2_min = v_chi2[i];
      dm2_min = dm2;
      sin22th_mue_min = sin22th_mue;
      sin22th_mumu_min = sin22th_mumu;
      sin22th_ee_min = sin22th_ee;
      m41_min = v_m41[i];
      theta14_min = v_theta14[i];
      theta24_min = v_theta24[i];
      imin = i;
    }
  }

  std::cout << "chi2 min: " << chi2_min << " dm2 min: " << dm2_min << std::endl;
  std::cout << "sin22th_mue min: " << sin22th_mue_min << std::endl;
  std::cout << "sin22th_mumu min: " << sin22th_mumu_min << std::endl;
  std::cout << "sin22th_ee min: " << sin22th_ee_min << std::endl;
  std::cout << "RAW MIN (m41, theta14, theta24): ( " << m41_min << ", " << theta14_min << ", " << theta24_min << " )" << std::endl;
  std::cout << "i: " << imin << std::endl;

  // Create output File
  std::string outfile = rdr.tag + "_proc_v1.root";
  std::cout << "Output File: " << outfile << std::endl;
  TFile *f = new TFile(outfile.c_str(), "RECREATE");
  if(f->IsZombie()){
    std::cout << "Error: couldn't create output file." << std::endl;
    return 0;
  }

  TTree * t_nueapp_1sig = new TTree("nueapp_1sig","nueapp_1sig");
  TTree * t_nueapp_99 = new TTree("nueapp_99","nueapp_99");
  TTree * t_nueapp_90 = new TTree("nueapp_90","nueapp_90");
  TTree * t_nueapp_90_exclusion = new TTree("nueapp_90_exclusion","nueapp_90_exclusion");
  TTree * t_nueapp_95_raster = new TTree("nueapp_95_raster","nueapp_95_raster");
  TTree * t_nuedis_1sig = new TTree("nuedis_1sig","nuedis_1sig");
  TTree * t_nuedis_99 = new TTree("nuedis_99","nuedis_99");
  TTree * t_nuedis_90 = new TTree("nuedis_90","nuedis_90");
  TTree * t_nuedis_90_exclusion = new TTree("nuedis_90_exclusion","nuedis_90_exclusion");
  TTree * t_nuedis_95_raster = new TTree("nuedis_95_raster","nuedis_95_raster");
  TTree * t_numudis_1sig = new TTree("numudis_1sig","numudis_1sig");
  TTree * t_numudis_99 = new TTree("numudis_99","numudis_99");
  TTree * t_numudis_90 = new TTree("numudis_90","numudis_90");
  TTree * t_numudis_90_exclusion = new TTree("numudis_90_exclusion","numudis_90_exclusion");
  TTree * t_numudis_95_raster = new TTree("numudis_95_raster","numudis_95_raster");
  TTree * t_numudis_90_raster = new TTree("numudis_90_raster","numudis_90_raster");
  TTree * t_numudis_chiogram = new TTree("numudis_chiogram","numudis_chiogram");


  // Nue Appearance
  t_nueapp_1sig->Branch("chi2",&chi2,"chi2/D");
  t_nueapp_99->Branch("chi2",&chi2,"chi2/D");
  t_nueapp_90->Branch("chi2",&chi2,"chi2/D");
  t_nueapp_90_exclusion->Branch("chi2",&chi2,"chi2/D");
  t_nueapp_95_raster->Branch("chi2",&chi2,"chi2/D");
  t_nueapp_1sig->Branch("sin22th",&sin22th_mue,"sin22th/D");
  t_nueapp_99->Branch("sin22th",&sin22th_mue,"sin22th/D");
  t_nueapp_90->Branch("sin22th",&sin22th_mue,"sin22th/D");
  t_nueapp_90_exclusion->Branch("sin22th",&sin22th_mue,"sin22th/D");
  t_nueapp_95_raster->Branch("sin22th",&sin22th,"sin22th/D");
  t_nueapp_1sig->Branch("dm2",&dm2,"dm2/D");
  t_nueapp_99->Branch("dm2",&dm2,"dm2/D");
  t_nueapp_90->Branch("dm2",&dm2,"dm2/D");
  t_nueapp_90_exclusion->Branch("dm2",&dm2,"dm2/D");
  t_nueapp_95_raster->Branch("dm2",&dm2,"dm2/D");

  // numu disappearance
  t_numudis_99->Branch("chi2",&chi2,"chi2/D");
  t_numudis_90->Branch("chi2",&chi2,"chi2/D");
  t_numudis_90_exclusion->Branch("chi2",&chi2,"chi2/D");
  t_numudis_95_raster->Branch("chi2",&chi2,"chi2/D");
  t_numudis_90_raster->Branch("chi2",&chi2,"chi2/D");
  t_numudis_chiogram->Branch("chi2",&chi2,"chi2/D");
  t_numudis_99->Branch("sin22th",&sin22th_mumu,"sin22th/D");
  t_numudis_90->Branch("sin22th",&sin22th_mumu,"sin22th/D");
  t_numudis_90_exclusion->Branch("sin22th",&sin22th_mumu,"sin22th/D");
  t_numudis_95_raster->Branch("sin22th",&sin22th,"sin22th/D");
  t_numudis_90_raster->Branch("sin22th",&sin22th,"sin22th/D");
  t_numudis_chiogram->Branch("sin22th",&sin22th,"sin22th/D");
  t_numudis_99->Branch("dm2",&dm2,"dm2/D");
  t_numudis_90->Branch("dm2",&dm2,"dm2/D");
  t_numudis_90_exclusion->Branch("dm2",&dm2,"dm2/D");
  t_numudis_95_raster->Branch("dm2",&dm2,"dm2/D");
  t_numudis_90_raster->Branch("dm2",&dm2,"dm2/D");
  t_numudis_chiogram->Branch("dm2",&dm2,"dm2/D");
  t_numudis_1sig->Branch("chi2",&chi2,"chi2/D");
  t_numudis_1sig->Branch("sin22th",&sin22th_mumu,"sin22th/D");
  t_numudis_1sig->Branch("dm2",&dm2,"dm2/D");



  // nue disappearance
  t_nuedis_99->Branch("chi2",&chi2,"chi2/D");
  t_nuedis_90->Branch("chi2",&chi2,"chi2/D");
  t_nuedis_90_exclusion->Branch("chi2",&chi2,"chi2/D");
  t_nuedis_95_raster->Branch("chi2",&chi2,"chi2/D");
  t_nuedis_99->Branch("sin22th",&sin22th_ee,"sin22th/D");
  t_nuedis_90->Branch("sin22th",&sin22th_ee,"sin22th/D");
  t_nuedis_90_exclusion->Branch("sin22th",&sin22th_ee,"sin22th/D");
  t_nuedis_95_raster->Branch("sin22th",&sin22th,"sin22th/D");
  t_nuedis_99->Branch("dm2",&dm2,"dm2/D");
  t_nuedis_90->Branch("dm2",&dm2,"dm2/D");
  t_nuedis_90_exclusion->Branch("dm2",&dm2,"dm2/D");
  t_nuedis_95_raster->Branch("dm2",&dm2,"dm2/D");
  t_nuedis_1sig->Branch("chi2",&chi2,"chi2/D");
  t_nuedis_1sig->Branch("sin22th",&sin22th_ee,"sin22th/D");
  t_nuedis_1sig->Branch("dm2",&dm2,"dm2/D");


  double c_90(4.71), c_99(9.21), c_1sig(2.30), c_90_1sided(3.22);  // 2dof
  int im, is;
  for(int i = 0; i < v_chi2.size(); i++){
    dm2 = pow(v_m41[i],2);
    sin22th_mue = pow(sin(2*v_theta14[i]),2)*pow(sin(2*v_theta24[i]),2);
    sin22th_mumu = 4 * pow(cos(v_theta14[i]),2) * pow(sin(v_theta24[i]),2) * (1 - pow(cos(v_theta14[i]),2) * pow(sin(v_theta24[i]),2));
    sin22th_ee = pow(sin(2*v_theta14[i]),2);
    chi2 = v_chi2[i];

    if(chi2 < chi2_min + c_1sig){
      t_nueapp_1sig->Fill();
      t_numudis_1sig->Fill();
      t_nuedis_1sig->Fill();
    }

    if(chi2 < chi2_min + c_90){
      t_nueapp_90->Fill();
      t_numudis_90->Fill();
      t_nuedis_90->Fill();
    }
    if(chi2 < chi2_min + c_99){
      t_nueapp_99->Fill();
      t_numudis_99->Fill();
      t_nuedis_99->Fill();
    }
    if(chi2 > chi2_min + c_90_1sided){
      t_nueapp_90_exclusion->Fill();
      t_numudis_90_exclusion->Fill();
      t_nuedis_90_exclusion->Fill();
    }

    // Also fill up grids for raster scans
    im = ValueToIndex(dm2,pow(m41_lowbound,2),pow(m41_hibound,2),rdr.gridpts_dm2);
    is = ValueToIndex(sin22th_mue,sin22th_lowbound,sin22th_hibound,rdr.gridpts_sin22th);
    if(is < rdr.gridpts_sin22th && im < rdr.gridpts_dm2 && is > -1){
      if(chi2grid_mue[is][im] == 0)
        chi2grid_mue[is][im] = v_chi2[i];
      else
        chi2grid_mue[is][im] = min(chi2grid_mue[is][im],v_chi2[i]);
    }
    is = ValueToIndex(sin22th_mumu,sin22th_lowbound,sin22th_hibound,rdr.gridpts_sin22th);
    if(is < rdr.gridpts_sin22th && im < rdr.gridpts_dm2 && is > -1){
      if(chi2grid_mumu[is][im] == 0)
        chi2grid_mumu[is][im] = v_chi2[i];
      else
        chi2grid_mumu[is][im] = min(chi2grid_mumu[is][im],v_chi2[i]);
    }
    is = ValueToIndex(sin22th_ee,sin22th_lowbound,sin22th_hibound,rdr.gridpts_sin22th);
    if(is < rdr.gridpts_sin22th && im < rdr.gridpts_dm2 && is > -1){
      if(chi2grid_ee[is][im] == 0)
        chi2grid_ee[is][im] = v_chi2[i];
      else
        chi2grid_ee[is][im] = min(chi2grid_ee[is][im],v_chi2[i]);
    }
  }
  t_nueapp_99->Write();
  t_nueapp_90->Write();
  t_nueapp_90_exclusion->Write();
  t_numudis_99->Write();
  t_numudis_90->Write();
  t_numudis_90_exclusion->Write();
  t_nuedis_99->Write();
  t_nuedis_90->Write();
  t_nuedis_90_exclusion->Write();
  t_nueapp_1sig->Write();
  t_numudis_1sig->Write();
  t_nuedis_1sig->Write();

  // perform raster scans
  for(int dm = 0; dm < rdr.gridpts_dm2; dm++){
    double chi2minRaster = 3000.;
    int sinsStartRaster = 0;

    // Find minimum point for this particular dm2
    for(int sins = 0; sins < rdr.gridpts_sin22th; sins++){
      if(chi2grid_mue[sins][dm] < chi2minRaster && chi2grid_mue[sins][dm] >= chi2_min){
        chi2minRaster = chi2grid_mue[sins][dm];
        sinsStartRaster = sins;
      }
    }
    for(int sins = sinsStartRaster; sins < rdr.gridpts_sin22th; sins++){
      chi2 = chi2grid_mue[sins][dm];
      if(chi2 > 2.706 + chi2minRaster){ //95%
        dm2 = IndexToValue(dm,pow(m41_lowbound,2),pow(m41_hibound,2),rdr.gridpts_dm2);
        sin22th = IndexToValue(sins,sin22th_lowbound,sin22th_hibound,rdr.gridpts_sin22th);
        t_nueapp_95_raster->Fill();
        break;
      }
    }
  }

  for(int dm = 0; dm < rdr.gridpts_dm2; dm++){  // numudis
    double chi2minRaster = 3000.;
    int sinsStartRaster = 0;

    // Find minimum point for this particular dm2
    for(int sins = 0; sins < rdr.gridpts_sin22th; sins++){
      chi2 = chi2grid_mumu[sins][dm];
      dm2 = IndexToValue(dm,pow(m41_lowbound,2),pow(m41_hibound,2),rdr.gridpts_dm2);
      sin22th = IndexToValue(sins,sin22th_lowbound,sin22th_hibound,rdr.gridpts_sin22th);
      if(chi2grid_mumu[sins][dm] < chi2minRaster && chi2grid_mumu[sins][dm] >= chi2_min){
        chi2minRaster = chi2grid_mumu[sins][dm];
        sinsStartRaster = sins;
      }
      t_numudis_chiogram->Fill();
    }
    for(int sins = sinsStartRaster; sins < rdr.gridpts_sin22th; sins++){
      chi2 = chi2grid_mumu[sins][dm];
      if(chi2 > 2.706 + chi2minRaster){ //95%
        dm2 = IndexToValue(dm,pow(m41_lowbound,2),pow(m41_hibound,2),rdr.gridpts_dm2);
        sin22th = IndexToValue(sins,sin22th_lowbound,sin22th_hibound,rdr.gridpts_sin22th);
        t_numudis_95_raster->Fill();
        break;
      }
    }
    for(int sins = sinsStartRaster; sins < rdr.gridpts_sin22th; sins++){
      chi2 = chi2grid_mumu[sins][dm];
      if(chi2 > 1.28 + chi2minRaster){ //90%
        dm2 = IndexToValue(dm,pow(m41_lowbound,2),pow(m41_hibound,2),rdr.gridpts_dm2);
        sin22th = IndexToValue(sins,sin22th_lowbound,sin22th_hibound,rdr.gridpts_sin22th);
        t_numudis_90_raster->Fill();
        break;
      }
    }
  }

  for(int dm = 0; dm < rdr.gridpts_dm2; dm++){  // nuedis
    double chi2minRaster = 3000.;
    int sinsStartRaster = 0;

    // Find minimum point for this particular dm2
    for(int sins = 0; sins < rdr.gridpts_sin22th; sins++){
      if(chi2grid_ee[sins][dm] < chi2minRaster && chi2grid_ee[sins][dm] >= chi2_min){
        chi2minRaster = chi2grid_ee[sins][dm];
        sinsStartRaster = sins;
      }
    }
    for(int sins = sinsStartRaster; sins < rdr.gridpts_sin22th; sins++){
      chi2 = chi2grid_ee[sins][dm];
      if(chi2 > 2.706 + chi2minRaster){ //95%
        dm2 = IndexToValue(dm,pow(m41_lowbound,2),pow(m41_hibound,2),rdr.gridpts_dm2);
        sin22th = IndexToValue(sins,sin22th_lowbound,sin22th_hibound,rdr.gridpts_sin22th);
        t_nuedis_95_raster->Fill();
        break;
      }
    }
  }
  std::cout << "nueapp raster has " << t_nueapp_95_raster->GetEntries() << " elements" << std::endl;
  std::cout << "nuedis raster has " << t_nuedis_95_raster->GetEntries() << " elements" << std::endl;
  std::cout << "numudis raster has " << t_numudis_95_raster->GetEntries() << " elements" << std::endl;
  t_nueapp_95_raster->Write();
  t_numudis_95_raster->Write();
  t_numudis_90_raster->Write();
  t_nuedis_95_raster->Write();
  t_numudis_chiogram->Write();

  f->Close();
  return 0;
}


int main(int argc, char* argv[]){

  std::string xml = "";
  int iarg = 0;
  opterr=1;
  int index;

  const struct option longopts[] = {
    {"xml", 		required_argument, 	0, 'x'},
	  {0,			no_argument, 		0,  0},
  };

  while(iarg != -1){
    iarg = getopt_long(argc,argv, "x:t:", longopts, &index);

    switch(iarg){
		  case 'x':
			  xml = optarg;
			  break;
      case '?':
		  case 'h':
			  std::cout<<"I need an input, friend."<<std::endl;
			  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
			  return 0;
	  }
  }
  if(xml == ""){
    std::cout << "Gimme an XML input or I won't start, I swear to god." << std::endl;
    return 0;
  }

  ntupleProcess(xml);
  return 0;
}
