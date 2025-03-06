#include "OutputHandler.h"

OutputHandler::OutputHandler() : calculator{0}, stathandler{0}
{
}
OutputHandler::OutputHandler(std::string name, std::string cont_suffix) : calculator{0}, stathandler{0}
{
  initialise(name, cont_suffix);
}
OutputHandler::~OutputHandler()
{
  delete stathandler;
  delete calculator;
}
void OutputHandler::initialise(std::string name, std::string cont_suffix)
{
  // Read input file
  TFile* infile = new TFile(name.c_str(), "READ");
  if (!infile) {
    printf("File %s not found!\n", name.c_str());
    return;
  };
  if (infile->IsZombie()) {
    printf("File %s is a zombie!\n", name.c_str());
    delete infile;
    return;
  };

  if(calculator) delete calculator;
  calculator = new Calculator();
  calculator->setContainerSuffix(cont_suffix);
  calculator->readData(infile);
  if(stathandler) delete stathandler;
  stathandler = new StatHandler();

  return;
}
TH1D* OutputHandler::getHist(Calculator::ObsFunc func, StatHandler::StatMethod method)
{
  std::vector<TH1D*> samplehists;
  for (int i = -1; i < 10; ++i){
    samplehists.push_back(calculator->calculate(func, i));
    if(!samplehists[i+1]) std::cout << "Warning: samplehist " << i << " not calculated" << std::endl;
  }
  auto reth = stathandler->getNominal(samplehists, method);
  // sysHandler->applySyst(reth);
  return reth;
}
TH1D* OutputHandler::getHist(Calculator::PtDiffFunc func, double min, double max, StatHandler::StatMethod method)
{
  std::vector<TH1D*> samplehists;
  for (int i = -1; i < 10; ++i){
    samplehists.push_back(calculator->calculate(func, i, min, max));
    if(!samplehists[i+1]) std::cout << "Warning: samplehist " << i << " not calculated" << std::endl;
/*     for(int bin = 1; bin <= samplehists[i+1]->GetNbinsX();++bin){
      std::cout << samplehists[i+1]->GetBinContent(bin) << std::endl;
    }
 */  }
  auto reth = stathandler->getNominal(samplehists, method);
  // sysHandler->applySyst(reth);
  return reth;
}

TH1D* OutputHandler::getNorm(int n, int m, int k, StatHandler::StatMethod method)
{
  std::vector<TH1D*> samplehists;
  for (int i = -1; i < 10; ++i){
    samplehists.push_back(calculator->getNormNC(n,m,k,i));
    if(!samplehists[i+1]) std::cout << "Warning: samplehist " << i << " not calculated" << std::endl;
  }
  auto reth = stathandler->getNominal(samplehists, method);
  // sysHandler->applySyst(reth);
  return reth;
}