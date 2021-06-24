#include "datasets.h"
#include "XMLTools.h"


int FitReader::Load(std::string xml){
	std::cout << "YOU DON'T NEED THIS FOR NUMU DISAPPEARANCE. CARRY ON" << std::endl;
  return 0;
}

int ProcessReader::Load(std::string xml){

  // Load XML file
  XMLDocument doc;
  if (doc.LoadFile(xml.c_str())){
    std::cout << "Couldn't load XML file! I quit." << std::endl;
    return 1;
  }
  XMLHandle hDoc(&doc);

  // We'll have an oscillator and several datasets
  XMLElement *pData, *pProc;
  pData = doc.FirstChildElement("dataset");
  pProc = doc.FirstChildElement("procopts");

  std::string dset,dloc;

  if (!pData){
    std::cout << "No datasets in config. Outta here." << std::endl;
    return 1;
  }
  else while(pData){
    dset = pData->Attribute("name");
    if(stoi(pData->Attribute("use"))){
      dloc = pData->Attribute("loc");
      std::cout << "Adding tree " << dset << " from " << dloc << "." << std::endl;
      data_names.push_back(dset);
      data_files.push_back(new TFile(dloc.c_str(),"READ"));
      data_trees.push_back((TTree*)data_files.back()->Get(dset.c_str()));
      data_osc.push_back(stoi(pData->Attribute("osc")));
    }

    pData = pData->NextSiblingElement("dataset");
  }
  if(data_names.size()==0){
    std::cout << "No valid datasets requested." << std::endl;
    return 1;
  }

  if (!pProc){
    std::cout << "No process options. Fuck this." << std::endl;
    return 1;
  }
  else{
    tag = pProc->Attribute("tag");
    gridpts_dm2 = atoi(pProc->Attribute("gridpts_dm2"));
    gridpts_sin22th = atoi(pProc->Attribute("gridpts_sin22th"));
  }
  return 0;
}
