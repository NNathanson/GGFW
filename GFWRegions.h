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
#ifndef GGFW_GFWREGIONS_H
#define GGFW_GFWREGIONS_H


#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <TObject.h>
#include <TMath.h>


template <typename T, typename... Args>
int CheckSameSize(const std::vector<T>& first, const std::vector<Args>&... rest)
{
  size_t size = first.size();
  bool allSameSize = ((size == rest.size()) && ...);

  return allSameSize ? size : -1;
}

/// \class GFWRegions
/// \brief Class which implements configurable regions for the GFW
///

class GFWRegions
{
 public:
  GFWRegions(std::vector<std::string> names_ = {"refN", "refP", "refFull"}, std::vector<float> etaminvals_ = {-0.8, 0.4, -0.8},
             std::vector<float> etamaxvals_ = {-0.4, 0.8, 0.8},
             std::vector<int> pTDifs_ = {0, 0, 0},
             std::vector<int> bitmasks_ = {1, 1, 1}) : names{std::move(names_)}, etaminvals{std::move(etaminvals_)}, etamaxvals{std::move(etamaxvals_)}, pTDifs{std::move(pTDifs_)}, bitmasks{std::move(bitmasks_)} {};

  // auto Print() const
  // {
  //   for (size_t i = 0; i < names.size(); ++i) {
  //     LOGF(info, "{%s, %.1f, %.1f, %d, %d}", names[i].c_str(), etaminvals[i], etamaxvals[i], pTDifs[i], bitmasks[i]);
  //   }
  //   return;
  // }

  auto GetSize() const { return CheckSameSize(names, etaminvals, etamaxvals, pTDifs, bitmasks); }

  void SetNames(std::vector<std::string> names_) { names = std::move(names_); }
  const auto& GetNames() const { return names; }

  void SetEtaMin(std::vector<float> etaminvals_) { etaminvals = std::move(etaminvals_); }
  const auto& GetEtaMin() const { return etaminvals; }

  void SetEtaMax(std::vector<float> etamaxvals_) { etamaxvals = std::move(etamaxvals_); }
  const auto& GetEtaMax() const { return etamaxvals; }

  void SetpTDifs(std::vector<int> pTDifs_) { pTDifs = std::move(pTDifs_); }
  const auto& GetpTDifs() const { return pTDifs; }

  void SetBitmasks(std::vector<int> bitmasks_) { bitmasks = std::move(bitmasks_); }
  const auto& GetBitmasks() const { return bitmasks; }

 private:
  std::vector<std::string> names;
  std::vector<float> etaminvals;
  std::vector<float> etamaxvals;
  std::vector<int> pTDifs;
  std::vector<int> bitmasks;
  ClassDefNV(GFWRegions, 1);
};

#endif