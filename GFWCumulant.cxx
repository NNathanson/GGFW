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

/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572 by A. Bilandzic et al.)
A part of <GFW.cxx/h>
A container to store Q vectors for one subevent with an extra layer to recursively calculate particle correlations.
If used, modified, or distributed, please aknowledge the author of this code.
*/

#include "GFWCumulant.h"
#include <iostream>

using std::complex;
using std::vector;

GFWCumulant::GFWCumulant() : fQvector(0),
                             fUsed(kBlank),
                             fNEntries(-1),
                             fN(1),
                             fPow(1),
                             fPt(1),
                             fFilledPts(0),
                             fInitialized(false) {}

GFWCumulant::~GFWCumulant() {}
void GFWCumulant::FillArray(int ptin, double phi, double weight, double SecondWeight)
{
  ////std::cout<<"in gfw cumulant fill array"<<std::endl;
  ////std::cout<<"ptin "<<ptin<<std::endl;
  ////std::cout<<"fPt "<<fPt<<std::endl;
  if (!fInitialized)
    CreateComplexVectorArray(1, 1, 1);
  if (fPt == 1)
    ptin = 0; // If one bin, then just fill it straight; otherwise, if ptin is out-of-range, do not fill
  else if (ptin < 0 || ptin >= fPt)
    return;
  fFilledPts[ptin] = true;
  ////std::cout<<"fN "<<fN<<std::endl;
  for (int lN = 0; lN < fN; lN++) {
    ////std::cout<<"test"<<std::endl;
    double lSin = sin(lN * phi); // No need to recalculate for each power
    double lCos = cos(lN * phi); // No need to recalculate for each power
    ////std::cout<<"lN "<<lN<<std::endl;
    ////std::cout<<"PW(lN) "<<PW(lN)<<std::endl;
    for (int lPow = 0; lPow < PW(lN); lPow++) {
      //std::cout<<"test 2"<<std::endl;
      double lPrefactor = 0;
      // Dont calculate it twice; multiplication is cheaper that power
      // Also, if second weight is specified, then keep the first weight with power no more than 1, and us the other weight otherwise
      // this is important when POIs are a subset of REFs and have different weights than REFs
      if (SecondWeight > 0 && lPow > 1)
        lPrefactor = pow(SecondWeight, lPow - 1) * weight;
      else
        lPrefactor = pow(weight, lPow);
      double qsin = lPrefactor * lSin;
      double qcos = lPrefactor * lCos;
      //std::cout<<"!!! about to fill fQvector with ["<<ptin<<"]["<<lN<<"]["<<lPow<<"]"<<std::endl;
      fQvector[ptin][lN][lPow] += complex<double>(qcos, qsin);
    }
  }
  Inc();
};
void GFWCumulant::ResetQs()
{
  if (!fNEntries)
    return; // If 0 entries, then no need to reset. Otherwise, if -1, then just initialized and need to set to 0.
  for (int i = 0; i < fPt; i++) {
    fFilledPts[i] = false;
    for (int lN = 0; lN < fN; lN++) {
      for (int lPow = 0; lPow < PW(lN); lPow++) {
        fQvector[i][lN][lPow] = fNullQ;
      }
    }
  }
  fNEntries = 0;
};
void GFWCumulant::DestroyComplexVectorArray()
{
  if (!fInitialized)
    return;
  for (int l_n = 0; l_n < fN; l_n++) {
    for (int i = 0; i < fPt; i++) {
      delete[] fQvector[i][l_n];
    }
  }
  for (int i = 0; i < fPt; i++) {
    delete[] fQvector[i];
  }
  delete[] fQvector;
  delete[] fFilledPts;
  fInitialized = false;
  fNEntries = -1;
};

void GFWCumulant::CreateComplexVectorArray(int N, int Pow, int Pt)
{
  DestroyComplexVectorArray();
  vector<int> pwv;
  for (int i = 0; i < N; i++)
    pwv.push_back(Pow);
  CreateComplexVectorArrayVarPower(N, pwv, Pt);
};
void GFWCumulant::CreateComplexVectorArrayVarPower(int N, vector<int> PowVec, int Pt)
{
  //std::cout<<"in CreateComplexVectorArrayVarPower"<<std::endl;
  //std::cout<<N<<std::endl;
  DestroyComplexVectorArray();
  fN = N;
  fPow = 0;
  fPt = Pt;
  fFilledPts = new bool[Pt];
  fPowVec = PowVec;
  fQvector = new complex<double>**[fPt];
  for (int i = 0; i < fPt; i++) {
    fQvector[i] = new complex<double>*[fN];
  }
  for (int l_n = 0; l_n < fN; l_n++) {
    for (int i = 0; i < fPt; i++) {
      fQvector[i][l_n] = new complex<double>[PW(l_n)];
    }
  }
  ResetQs();
  fInitialized = true;
  //std::cout<<fN<<std::endl;
};
complex<double> GFWCumulant::Vec(int n, int p, int ptbin)
{
  ////std::cout<<"in vec"<<std::endl;
  
  if (!fInitialized)
    return 0;
  ////std::cout<<fInitialized<<std::endl;

  if (ptbin >= fPt || ptbin < 0)
    ptbin = 0;

  ////std::cout<<"ptbin "<<ptbin<<std::endl;
  ////std::cout<<"n "<<n<<std::endl;
  if (n >= 0){
    ////std::cout<<"in if statement"<<std::endl;
    ////std::cout<<fQvector[0]<<std::endl;
    ////std::cout<<"before return"<<std::endl;
    
    return fQvector[ptbin][n][p];
  }
  ////std::cout<<"out of if statement"<<std::endl;
  return conj(fQvector[ptbin][-n][p]);
};
bool GFWCumulant::IsPtBinFilled(int ptb)
{
  if (!fFilledPts)
    return false;
  if (ptb > 0) {
    if (fPt == 1)
      ptb = 0;
    else if (ptb >= fPt)
      return false; // This is in case we are differential and going out of range for whatever reason.
  }
  return fFilledPts[ptb];
}
