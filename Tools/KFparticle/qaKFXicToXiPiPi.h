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

/// \file qaKFXicToXiPiPi.h
/// \author Carolina Reetz <c.reetz@cern.ch>, University of Heidelberg and GSI Darmstadt

#ifndef TOOLS_KFPARTICLE_QAKFXICTOXIPIPI_H_
#define TOOLS_KFPARTICLE_QAKFXICTOXIPIPI_H_

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::track;


namespace preselectionCutsXicToXiPiPi
{
static const int nBinsPt = 10;
static const int nCutVars = 11;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double binsPt[nBinsPt + 1] = {
  0.,
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  12.,
  24.,
  36.};
auto vecBinsPt = std::vector<double>{binsPt, binsPt + nBinsPt + 1};

// default values for the cuts                 m   ptXi ptPi0 ptPi1  dL    dLXY cosp cospXY dcaXY: Xi Pi0 Pi1
constexpr double cuts[nBinsPt][nCutVars] = {{0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1},   /* 0  < pT < 1  */
                                            {0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1},   /* 1  < pT < 2  */
                                            {0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1},   /* 2  < pT < 3  */
                                            {0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1},   /* 3  < pT < 4  */
                                            {0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1},   /* 4  < pT < 5  */
                                            {0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1},   /* 5  < pT < 6  */
                                            {0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1},   /* 6  < pT < 8  */
                                            {0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1},   /* 8  < pT < 12 */
                                            {0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1},   /* 12 < pT < 24 */
                                            {0.400, 0.4, 0.4, 0.4, 0.005, 0.005, 0.8, 0.8, 0.1, 0.1, 0.1}};  /* 24 < pT < 36 */

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9"};

// column labels
static const std::vector<std::string> labelsCutVar = {"m", "pT Xi", "pT Pi0", "pT Pi1", "chi2PCA", "max decay length", "max decay length XY", "cos pointing angle", "cos pointing angle XY",  "max DCAXY Xi", "max DCAXY Pi0", "max DCAXY Pi1"};
} // namespace preselectionCutsXicToXiPiPi