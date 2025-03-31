// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FlowPtContainer.h
/// \brief Class to handle angular and transverse momentum correlations
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch
#ifndef GGFW_FLOWPTCONTAINER_H
#define GGFW_FLOWPTCONTAINER_H

#include <algorithm>
#include <vector>
#include <complex>
#include <variant>
#include "BootstrapProfile.h"
#include "TNamed.h"
#include "TList.h"
#include "TCollection.h"
//#include "Framework/HistogramSpec.h"
#include "GFW.h"


namespace ggfw::eventweight
{
enum EventWeight {
  UnityWeight,
  TupleWeight
};
};

using namespace ggfw::eventweight;

class FlowPtContainer : public TNamed
{
 public:
  using FillType = std::variant<std::complex<double>, double>;
  FlowPtContainer();
  explicit FlowPtContainer(const char* name);
  ~FlowPtContainer();
  FlowPtContainer(const char* name, const char* title);
  //void initialise(Int_t nMultiBins, double* multiBins, const int& m, const GFWCorrConfigs& configs, const int& nsub = 10);
  void initialise(int nbinsx, double* xbins, const int& m, const GFWCorrConfigs& configs, const int& nsub = 10);
  void initialise(int nbinsx, double xlow, double xhigh, const int& m, const GFWCorrConfigs& configs, const int& nsub = 10);
  void initializeSubevent(int nbinsx, double* xbins, const int& m, const GFWCorrConfigs& configs, const int& nsub = 10);
  
  void fill(const double& w, const double& pt);
  void fillSub1(const double& w, const double& pt);
  void fillSub2(const double& w, const double& pt);
  void fillArray(FillType a, FillType b, double c, double d);
  
  int getVectorIndex(const int i, const int j) { return j * (mpar + 1) + i; }                                              // index for 2d array for storing pt correlations
  int getVectorIndex(const int i, const int j, const int k, const int l) { return i + j * 3 + k * 3 * 3 + l * 3 * 3 * 5; } // index for 4d array for std vnpt correlation - size 3x3x5x5
  
  void calculateCorrelations();
  void calculateSubeventCorrelations();

  void calculateCMTerms();
  void fillPtProfiles(const double& lMult, const double& rn);
  void fillSubeventPtProfiles(const double& lMult, const double& rn);


  void fillVnPtCorrProfiles(const double& lMult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask);
  void fillVnDeltaPtProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask);
  void fillVnDeltaPtStdProfiles(const double& centmult, const double& rn);
  void fillVnPtCorrStdProfiles(const double& centmult, const double& rn);
  void fillVnPtProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
  {
    if (fUseCentralMoments)
      fillVnDeltaPtProfiles(centmult, flowval, flowtuples, rn, mask);
    else
      fillVnPtCorrProfiles(centmult, flowval, flowtuples, rn, mask);
  }
  void fillVnPtStdProfiles(const double& centmult, const double& rn)
  {
    if (fUseCentralMoments)
      fillVnDeltaPtStdProfiles(centmult, rn);
    else
      fillVnPtCorrStdProfiles(centmult, rn);
  }
  


  void fillCMProfiles(const double& lMult, const double& rn);
  void fillCMSubeventProfiles(const double &lMult, const double &rn);


  TList* getCorrList() { return fCorrList; }
  TList* getCMTermList() { return fCMTermList; }
  TList* getCovList() { return fCovList; }
  TList* getSubList() { return fSubList; }
  TList* getSubCMList() { return fSubCMList; }
  
  void setEventWeight(const unsigned int& lWeight) { fEventWeight = lWeight; }
  void setUseCentralMoments(bool newval) { fUseCentralMoments = newval; }
  void setUseGapMethod(bool newval) { fUseGap = newval; }
  bool usesCentralMoments() { return fUseCentralMoments; }
  bool usesGap() { return fUseGap; }
  void rebinMulti(int nbins);
  void rebinMulti(int nbins, double* binedges);
  TH1* getCentralMomentHist(int ind, int m);
  TH1* getCumulantHist(int ind, int m);
  TH1* getCorrHist(int ind, int m);
  int getMpar() { return mpar; }
  Long64_t Merge(TCollection* collist);
  double orderedAddition(std::vector<double> vec);
  void createCentralMomentList();
  void calculateCentralMomentHists(std::vector<TH1*> inh, int ind, int m, TH1* hMpt);
  void createCumulantList();
  void calculateCumulantHists(std::vector<TH1*> inh, int ind);
  void clearVector()
  {
    sumP.clear();
    sumP.resize((mpar + 1) * (mpar + 1));
    cmVal.clear();
    cmDen.clear();
    fillCounter = 0;
    arr.clear();
    arr.resize(3 * 3 * 5 * 5, {0.0, 0.0});
    warr.clear();
    warr.resize(3 * 3 * 5 * 5, 0.0);
  };

  TList* fCMTermList;
  TList* fCorrList;
  TList* fCovList;
  TList* fSubList;
  TList* fSubCMList;
  TList* fCumulantList;
  TList* fCentralMomentList;

  int mpar;                  //!
  int fillCounter;           //!
  unsigned int fEventWeight; //!
  bool fUseCentralMoments;   //!
  bool fUseGap;              //!
  bool fUsePtSub;
  void mergeBSLists(TList* source, TList* target);
  TH1* raiseHistToPower(TH1* inh, double p);
  std::vector<double> sumP;              //!
  std::vector<double> insub1;              //!
  std::vector<double> insub2;              //!
  std::vector<double> corrNum;           //!
  std::vector<double> corrNum1;           //!
  std::vector<double> corrNum2;           //!
  std::vector<double> corrDen;           //!
  std::vector<double> corrDen1;           //!
  std::vector<double> corrDen2;           //!
  std::vector<double> cmVal;             //!
  std::vector<double> cmVal1;             //!
  std::vector<double> cmVal2;             //!
  std::vector<double> cmDen;             //!
  std::vector<double> cmDen1;             //!
  std::vector<double> cmDen2;             //!
  std::vector<std::complex<double>> arr; //!
  std::vector<double> warr;              //!
  template <typename T>
  double getStdAABBCC(T& inarr);
  template <typename T>
  double getStdAABBCD(T& inarr);
  template <typename T>
  double getStdAABBDD(T& inarr);
  template <typename T>
  double getStdAABBC(T& inarr);
  template <typename T>
  double getStdAABBD(T& inarr);
  template <typename T>
  double getStdABCC(T& inarr);
  template <typename T>
  double getStdABCD(T& inarr);
  template <typename T>
  double getStdABDD(T& inarr);
  template <typename T>
  double getStdABC(T& inarr);
  template <typename T>
  double getStdABD(T& inarr);
  template <typename T>
  double getStdABCCCC(T& inarr);
  template <typename T>
  double getStdABCCCD(T& inarr);
  template <typename T>
  double getStdABCCDD(T& inarr);
  template <typename T>
  double getStdABCDDD(T& inarr);
  template <typename T>
  double getStdABDDDD(T& inarr);
  template <typename T>
  double getStdABCCC(T& inarr);
  template <typename T>
  double getStdABCCD(T& inarr);
  template <typename T>
  double getStdABCDD(T& inarr);
  template <typename T>
  double getStdABDDD(T& inarr);

 private:
  static constexpr float FactorialArray[9] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320.};
  static constexpr int SignArray[9] = {1, -1, 1, -1, 1, -1, 1, -1, 1};
  ClassDef(FlowPtContainer, 2);
};

#endif