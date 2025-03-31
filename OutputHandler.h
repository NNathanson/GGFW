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

#ifndef OUTPUTHANDLER_H_
#define OUTPUTHANDLER_H_

#include <functional>
#include <string>
#include <vector>
#include "TH1.h"
#include "FlowContainer.h"
#include "FlowPtContainer.h"
#include "Calculator.h"
#include "StatHandler.h"

class OutputHandler : public TNamed
{
 public:
  OutputHandler();
  OutputHandler(std::string name, bool gap, std::string cont_suffix = "");

  ~OutputHandler();

  void initialise(std::string name, bool gap, std::string cont_suffix = "");

  TH1D* getHist(Calculator::ObsFunc func, StatHandler::StatMethod method = StatHandler::kVarianceOfSamples);
  TH1D* getHist(Calculator::PtDiffFunc func, double min, double max, StatHandler::StatMethod method = StatHandler::kVarianceOfSamples);
  TH1D* getNC(Calculator::ObsFunc funcNum, int n, int m, Calculator::ObsFunc funcDenomMpt, StatHandler::StatMethod method = StatHandler::kVarianceOfSamples);
  TH1D* getNCFull(Calculator::ObsFunc funcNum, int n, int m, Calculator::ObsFunc funcDenomMpt, StatHandler::StatMethod method = StatHandler::kVarianceOfSamples);

  void SetIDName(std::string name) { if(calculator) calculator->SetIDName(name); }
  void Rebin(int nbins, double* binedges);

 private:
  Calculator* calculator;   //!
  StatHandler* stathandler; //!

  ClassDef(OutputHandler, 1);
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_OUTPUTHANDLER_H_