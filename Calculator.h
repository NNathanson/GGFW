
#ifndef CALCULATOR_H_
#define CALCULATOR_H_

#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <map>
#include <string>
#include "TH1.h"
#include "TString.h"
#include "FlowContainer.h"
#include "FlowPtContainer.h"

class Calculator
{
 public:
  using ObsFunc = TH1D* (Calculator::*)(int);
  using PtDiffFunc = TH1D* (Calculator::*)(int, double, double);
  Calculator();
  Calculator(TFile* infile,bool gap);
  virtual ~Calculator();
  void readData(TFile* infile,bool gap);
  void setContainerSuffix(std::string cont_suffix) { if(!cont_suffix.empty() && cont_suffix.at(0)!='_') cont_suffix.insert(cont_suffix.begin(),'_'); current_suffix = cont_suffix; };
  void SetIDName(std::string name) { idname = name; id_changed = true;}
  TH1D* calculate(ObsFunc func, int i);
  TH1D* calculate(PtDiffFunc func, int i, double min, double max);
  void Rebin(int nbins, double* binedges);

  // Flow coefficients
  TH1D* getVNN(int nhar, int npar, std::string id, int i, bool onPt = kFALSE, double min = -1, double max = -1);
  TH1D* getV22(int i) { return getVNN(2, 2, "ChGap", i); }
  TH1D* getV32(int i) { return getVNN(3, 2, "ChGap", i); }
  TH1D* getV42(int i) { return getVNN(4, 2, "ChGap", i); }
  TH1D* getV24(int i) { return getVNN(2, 4, "ChFull", i); }
  TH1D* getV26(int i) { return getVNN(2, 6, "ChFull", i); }
  TH1D* getV28(int i) { return getVNN(2, 8, "ChFull", i); }
  TH1D* getV22pt(int i, double min, double max) { return getVNN(2, 2, "ChGap", i, kTRUE, min, max); }
  TH1D* getV32pt(int i, double min, double max) { return getVNN(3, 2, "ChGap", i, kTRUE, min, max); }
  TH1D* getV42pt(int i, double min, double max) { return getVNN(4, 2, "ChGap", i, kTRUE, min, max); }
  TH1D* getV24pt(int i, double min, double max) { return getVNN(2, 4, "ChFull", i, kTRUE, min, max); }
 
  TH1D* test0CV22pt2(int i);
  TH1D* test1CV22pt2(int i);
  TH1D* test2CV22pt2(int i);
  TH1D* test0CV22pt3(int i);
  TH1D* test1CV22pt3(int i);
  TH1D* test2CV22pt3(int i);
  TH1D* test3CV22pt3(int i);
  TH1D* test0CV22pt4(int i);
  TH1D* test1CV22pt4(int i);
  TH1D* test2CV22pt4(int i);
  TH1D* test3CV22pt4(int i);
  TH1D* test4CV22pt4(int i);
  TH1D* test0CV24pt(int i);
  TH1D* test1CV24pt(int i);
  TH1D* test0CV24pt2(int i);
  TH1D* test1CV24pt2(int i);
  TH1D* test2CV24pt2(int i);
  TH1D* v2test0CV24pt2(int i);
  TH1D* v2test1CV24pt2(int i);
  TH1D* v2test2CV24pt2(int i);
  TH1D* v2test3CV24pt2(int i);
  TH1D* v2test4CV24pt2(int i);

  // Flow correlations
  TH1D* getNN(int nhar, int npar, std::string id, int i);
  TH1D* get22(int i) { return getNN(2, 2, "ChGap", i); }
  TH1D* get24(int i) { return getNN(2, 4, "ChFull", i); }
  TH1D* get32(int i) { return getNN(3, 2, "ChGap", i); }
  TH1D* get42(int i) { return getNN(4, 2, "ChGap", i); }

  // Flow cumulants
  TH1D* getCNN(int nhar, int npar, std::string id, int i);
  TH1D* getC22(int i) { return getCNN(2, 2, "ChGap", i); }
  TH1D* getC24(int i) { return getCNN(2, 4, "ChFull", i); }
  TH1D* getC26(int i) { return getCNN(2, 6, "ChFull", i); }
  TH1D* getC28(int i) { return getCNN(2, 8, "ChFull", i); }
  TH1D* getC32(int i) { return getCNN(3, 2, "ChGap", i); }
  TH1D* getC42(int i) { return getCNN(4, 2, "ChGap", i); }

  // Flow-pt cumulants
  TH1D* getCVN2pt(int nhar, int i);
  TH1D* getCV22pt(int i) { return getCVN2pt(2, i); }
  TH1D* getCV32pt(int i) { return getCVN2pt(3, i); }
  TH1D* getCV42pt(int i) { return getCVN2pt(4, i); }
  // Could in theory be generalised to CVN2ptX, but will we ever need those? In any case trivial to change
  TH1D* getCV22pt2(int i);
  TH1D* getCV22pt3(int i);
  TH1D* getCV22pt4(int i);
  TH1D* getCV24pt(int i);
  TH1D* getCV24pt2(int i);
  // Flow-pt normalised cumulants

  // Flow-pt correlators - for checks only
  TH1D* getV22pt_Mpt0(int i);
  TH1D* getV22pt_Mpt1(int i);

  TH1D* getV32pt_Mpt0(int i);
  TH1D* getV32pt_Mpt1(int i);
  
  TH1D* getV42pt_Mpt0(int i);
  TH1D* getV42pt_Mpt1(int i);

  TH1D* getV24pt2_Mpt0(int i);
  TH1D* getV24pt2_Mpt1(int i);
  TH1D* getV24pt2_Mpt2(int i);

  // Pt-functions
  TH1D* getPtCorr(int i, int m);
  TH1D* getMpt(int i) { return getPtCorr(i, 1); }
  TH1D* getMpt2(int i) { return getPtCorr(i, 2); }
  TH1D* getMpt3(int i) { return getPtCorr(i, 3); }
  TH1D* getMpt4(int i) { return getPtCorr(i, 4); }
  TH1D* getMpt5(int i) { return getPtCorr(i, 5); }
  TH1D* getMpt6(int i) { return getPtCorr(i, 6); }
  TH1D* getMpt7(int i) { return getPtCorr(i, 7); }
  TH1D* getMpt8(int i) { return getPtCorr(i, 8); }

  TH1D* getKN(int i, int n);
  TH1D* getK2(int i) { return getKN(i, 2); }
  TH1D* getK3(int i) { return getKN(i, 3); }
  TH1D* getK4(int i) { return getKN(i, 4); }
  TH1D* getMuN(int i, int n);
  TH1D* getMu2(int i) { return getMuN(i, 2); }
  TH1D* getMu3(int i) { return getMuN(i, 3); }
  TH1D* getMu4(int i) { return getMuN(i, 4); }
 private:
  // Data members
  FlowContainer* flowcont;     //!
  FlowPtContainer* flowptcont; //!
  TList* corr_list;            //!
  std::string idname;
  bool id_changed;
  std::string current_suffix;

  ClassDef(Calculator,1);
};
namespace CalculatorFunctions {
  //Flow correlations
  extern Calculator::ObsFunc k22;
  extern Calculator::ObsFunc k32;
  extern Calculator::ObsFunc k42;
  extern Calculator::ObsFunc k24;

  //Flow cumulants
  extern Calculator::ObsFunc kC22;
  extern Calculator::ObsFunc kC24;
  extern Calculator::ObsFunc kC26;
  extern Calculator::ObsFunc kC28;
  extern Calculator::ObsFunc kC32;
  extern Calculator::ObsFunc kC42;

  //Flow coefficients
  extern Calculator::ObsFunc kV22;
  extern Calculator::ObsFunc kV32;
  extern Calculator::ObsFunc kV42;
  extern Calculator::ObsFunc kV24;
  extern Calculator::ObsFunc kV26;
  extern Calculator::ObsFunc kV28;

  extern Calculator::PtDiffFunc kV22pt;
  extern Calculator::PtDiffFunc kV32pt;
  extern Calculator::PtDiffFunc kV42pt;
  extern Calculator::PtDiffFunc kV24pt;

  //Flow-pt cumulants
  extern Calculator::ObsFunc kCV22pt;
  extern Calculator::ObsFunc kCV32pt;
  extern Calculator::ObsFunc kCV42pt;
  extern Calculator::ObsFunc kCV22pt2;
  extern Calculator::ObsFunc kCV22pt3;
  extern Calculator::ObsFunc kCV22pt4;
  extern Calculator::ObsFunc kCV24pt;
  extern Calculator::ObsFunc kCV24pt2;

  //Flow-pt correlations
  extern Calculator::ObsFunc kV22pt_Mpt0;
  extern Calculator::ObsFunc kV22pt_Mpt1;
  extern Calculator::ObsFunc kV24pt2_Mpt0;
  extern Calculator::ObsFunc kV24pt2_Mpt1;
  extern Calculator::ObsFunc kV24pt2_Mpt2;

  //pt correlations
  extern Calculator::ObsFunc kMpt;
  extern Calculator::ObsFunc kMpt2;
  extern Calculator::ObsFunc kMpt3;
  extern Calculator::ObsFunc kMpt4;

  //pt cumulants
  extern Calculator::ObsFunc kK2;
  extern Calculator::ObsFunc kK3;
  extern Calculator::ObsFunc kK4;
  extern Calculator::ObsFunc kMu2;
  extern Calculator::ObsFunc kMu3;
  extern Calculator::ObsFunc kMu4;
}

#endif // PWGCF_GENERICFRAMEWORK_CORE_CALCULATOR_H_