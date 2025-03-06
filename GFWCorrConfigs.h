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
Class steers the initialization and calculation of n-particle correlations. Uses recursive function, all terms are calculated only once.
Latest version includes the calculation of any number of gaps and any combination of harmonics (including eg symmetric cumulants, etc.)
If used, modified, or distributed, please acknowledge the author of this code.
*/
#ifndef GGFW_GFWCORRCONFIGS_H
#define GGFW_GFWCORRCONFIGS_H


#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <TObject.h>
#include <TMath.h>

template <typename T>
int CheckSameSize(const std::vector<T>& first)
{
  return first.size();
}

/// \class GFWCorrConfigs
/// \brief Class which implements configurable correlations for GFW
///
class GFWCorrConfigs
{
 public:
  GFWCorrConfigs(std::vector<std::string> corrs_ = {"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}",
                                                    "refFull {2 2 -2 -2}"},
                 std::vector<std::string> heads_ = {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"},
                 std::vector<int> pTDifs_ = {0, 0, 0, 0, 0}, std::vector<int> pTCorrMasks_ = {15, 1, 1, 0, 0}) : corrs{std::move(corrs_)}, heads{std::move(heads_)}, pTDifs{std::move(pTDifs_)}, pTCorrMasks{std::move(pTCorrMasks_)} {};

  // auto Print() const
  // {
  //   for (size_t i = 0; i < corrs.size(); ++i) {
  //     LOGF(info, "{%s,%s,%d,%d}", heads[i].c_str(), corrs[i].c_str(), pTDifs[i], pTCorrMasks[i]);
  //   }
  //   return;
  // }

  auto GetSize() const { return CheckSameSize(corrs, heads, pTDifs); }

  void SetCorrs(std::vector<std::string> corrs_) { corrs = std::move(corrs_); }
  const auto& GetCorrs() const { return corrs; }

  void SetHeads(std::vector<std::string> heads_) { heads = std::move(heads_); }
  const auto& GetHeads() const { return heads; }

  void SetpTDifs(std::vector<int> pTDifs_) { pTDifs = std::move(pTDifs_); }
  const auto& GetpTDifs() const { return pTDifs; }

  void SetpTCorrMasks(std::vector<int> pTCorrMasks_) { pTCorrMasks = std::move(pTCorrMasks_); }
  const auto& GetpTCorrMasks() const { return pTCorrMasks; }

 private:
  std::vector<std::string> corrs;
  std::vector<std::string> heads;
  std::vector<int> pTDifs;
  std::vector<int> pTCorrMasks;
  ClassDefNV(GFWCorrConfigs, 1);
};


#endif