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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "VarManager.h"

void o2::aod::dqhistograms::DefineHistograms(HistogramManager* hm, const char* histClass, const char* groupName, const char* subGroupName)
{
  //
  // Add a predefined group of histograms to the HistogramManager hm and histogram class histClass
  // NOTE: The subGroupName argument may contain several keywords, but the user should take care of ambiguities. TODO: fix it!
  // NOTE: All of the histograms which match any of the group or subgroup names will be added to the same histogram class !!
  //            So one has to make sure not to mix e.g. event-wise with track-wise histograms
  // NOTE: The subgroup name can be empty. In this case just a minimal set of histograms corresponding to the group name will be defined
  //
  TString groupStr = groupName;
  groupStr.ToLower();
  TString subGroupStr = subGroupName;
  subGroupStr.ToLower();
  if (!groupStr.CompareTo("event")) {
    if (!subGroupStr.Contains("generator")) {
      hm->AddHistogram(histClass, "VtxZ", "Vtx Z", false, 60, -15.0, 15.0, VarManager::kVtxZ);
      hm->AddHistogram(histClass, "VtxZ_Run", "Vtx Z", true, VarManager::GetDummyNRuns(), -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo, 60, -15.0, 15.0, VarManager::kVtxZ, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "BC", "Event per BC", false, 3564, 0.0, 3564.0, VarManager::kBCOrbit);
    }
    if (subGroupStr.Contains("trigger")) {
      hm->AddHistogram(histClass, "IsINT7", "Is INT7", false, 2, -0.5, 1.5, VarManager::kIsINT7);
      if (subGroupStr.Contains("muon") || subGroupStr.Contains("all")) {
        hm->AddHistogram(histClass, "IsINT7inMUON", "INT7inMUON", false, 2, -0.5, 1.5, VarManager::kIsINT7inMUON);
        hm->AddHistogram(histClass, "IsMuonSingleLowPt7", "Is MuonSingleLowPt7", false, 2, -0.5, 1.5, VarManager::kIsMuonSingleLowPt7);
        hm->AddHistogram(histClass, "IsMuonSingleHighPt7", "Is MuonSingleHighPt7", false, 2, -0.5, 1.5, VarManager::kIsMuonSingleHighPt7);
        hm->AddHistogram(histClass, "IsMuonUnlikeLowPt7", "Is MuonUnlikeLowPt7", false, 2, -0.5, 1.5, VarManager::kIsMuonUnlikeLowPt7);
        hm->AddHistogram(histClass, "IsMuonLikeLowPt7", "Is MuonLikeLowPt7", false, 2, -0.5, 1.5, VarManager::kIsMuonLikeLowPt7);
      }
      if (subGroupStr.Contains("up") || subGroupStr.Contains("all")) {
        hm->AddHistogram(histClass, "IsCUP8", "CUP8", false, 2, -0.5, 1.5, VarManager::kIsCUP8);
        hm->AddHistogram(histClass, "IsCUP9", "CUP9", false, 2, -0.5, 1.5, VarManager::kIsCUP9);
        hm->AddHistogram(histClass, "IsMUP10", "MUP10", false, 2, -0.5, 1.5, VarManager::kIsMUP10);
        hm->AddHistogram(histClass, "IsMUP11", "MUP11", false, 2, -0.5, 1.5, VarManager::kIsMUP11);
      }
      if (subGroupStr.Contains("emc") || subGroupStr.Contains("all")) {
        hm->AddHistogram(histClass, "IsEMC7", "EMC7", false, 2, -0.5, 1.5, VarManager::kIsEMC7);
      }
    }
    if (subGroupStr.Contains("time")) {
      hm->AddHistogram(histClass, "CollTime", "Coll. time wrt BC", false, 100, 0.0, 100.0, VarManager::kCollisionTime);
      hm->AddHistogram(histClass, "CollTimeRes", "Coll. time resolution", false, 100, 0.0, 100.0, VarManager::kCollisionTimeRes);
      hm->AddHistogram(histClass, "CollTime_VtxZ", "Coll. time wrt BC vs vtx-z", false, 50, -15.0, 15., VarManager::kVtxZ, 100, 0.0, 100.0, VarManager::kCollisionTime);
      hm->AddHistogram(histClass, "CollTimeRes_VtxZ", "Coll. time resolution ", false, 50, -15.0, 15., VarManager::kVtxZ, 100, 0.0, 100.0, VarManager::kCollisionTimeRes);
      hm->AddHistogram(histClass, "CollTimeRes_MultTPC", "Coll. time resolution ", false, 50, 0.0, 500., VarManager::kMultTPC, 100, 0.0, 100.0, VarManager::kCollisionTimeRes);
      hm->AddHistogram(histClass, "CollTimeRes_MultPV", "Coll. time resolution ", false, 50, 0.0, 500., VarManager::kVtxNcontribReal, 100, 0.0, 100.0, VarManager::kCollisionTimeRes);
    }
    if (subGroupStr.Contains("vtx")) {
      hm->AddHistogram(histClass, "VtxX", "Vtx X", false, 200, -0.1, 0.1, VarManager::kVtxX);
      hm->AddHistogram(histClass, "VtxY", "Vtx Y", false, 200, -0.1, 0.1, VarManager::kVtxY);
      hm->AddHistogram(histClass, "VtxYVtxX", "Vtx Y vs Vtx X", false, 100, -0.1, 0.1, VarManager::kVtxX, 100, -0.1, 0.1, VarManager::kVtxY);
    }
    if (subGroupStr.Contains("vtxpp")) {
      hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 100.0, VarManager::kVtxNcontrib);
    }
    if (subGroupStr.Contains("vtxPbPb")) {
      hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 20000.0, VarManager::kVtxNcontrib);
    }
    if (subGroupStr.Contains("cent")) {
      hm->AddHistogram(histClass, "CentV0M", "CentV0M", false, 100, 0., 100., VarManager::kCentVZERO);
      hm->AddHistogram(histClass, "CentV0M_vtxZ", "CentV0M vs Vtx Z", false, 60, -15.0, 15.0, VarManager::kVtxZ, 20, 0., 100., VarManager::kCentVZERO);
      hm->AddHistogram(histClass, "CentFT0C", "CentFT0C", false, 100, 0., 100., VarManager::kCentFT0C);
      hm->AddHistogram(histClass, "CentFT0C_vtxZ", "CentFT0C vs Vtx Z", false, 60, -15.0, 15.0, VarManager::kVtxZ, 20, 0., 100., VarManager::kCentFT0C);
      hm->AddHistogram(histClass, "CentFT0C_MultTPC", "CentFT0C vs MultTPC", false, 100, 0., 100., VarManager::kCentFT0C, 50, 0., 50., VarManager::kMultTPC);
      hm->AddHistogram(histClass, "CentFT0C_Run", "Cent FT0C", true, VarManager::GetDummyNRuns(), -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo, 100, 0., 100., VarManager::kCentFT0C, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
    }
    if (subGroupStr.Contains("mult")) {
      if (subGroupStr.Contains("pp")) {
        hm->AddHistogram(histClass, "MultTPC", "MultTPC", false, 300, 0.0, 500.0, VarManager::kMultTPC);
        hm->AddHistogram(histClass, "MultFV0A", "MultFV0A", false, 300, 0.0, 500.0, VarManager::kMultFV0A);
        hm->AddHistogram(histClass, "MultFT0A", "MultFT0A", false, 300, 0.0, 200.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "MultFT0C", "MultFT0C", false, 300, 0.0, 300.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultFDDA", "MultFDDA", false, 300, 0.0, 300.0, VarManager::kMultFDDA);
        hm->AddHistogram(histClass, "MultFDDC", "MultFDDC", false, 300, 0.0, 50.0, VarManager::kMultFDDC);
        hm->AddHistogram(histClass, "MultTracklets", "MultTracklets", false, 100, 0.0, 250.0, VarManager::kMultTracklets);
        hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 100.0, VarManager::kVtxNcontribReal);
        hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 100.0, VarManager::kVtxNcontrib);
        hm->AddHistogram(histClass, "MultTPC_MultFV0A", "MultTPC vs MultFV0A", false, 100, 0, 500.0, VarManager::kMultTPC, 100, 0, 500.0, VarManager::kMultFV0A);
        hm->AddHistogram(histClass, "MultTPC_MultFT0A", "MultTPC vs MultFT0A", false, 100, 0, 500.0, VarManager::kMultTPC, 100, 0, 200.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "MultTPC_MultFT0C", "MultTPC vs MultFT0C", false, 100, 0, 500.0, VarManager::kMultTPC, 100, 0, 300.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultFT0A_MultFT0C", "MultFT0A vs MultFT0C", false, 100, 0, 200.0, VarManager::kMultFT0A, 100, 0, 300.0, VarManager::kMultFT0C);
      } else {
        hm->AddHistogram(histClass, "MultTPC", "MultTPC", false, 100, 0.0, 25000.0, VarManager::kMultTPC);
        hm->AddHistogram(histClass, "MultFV0A", "MultFV0A", false, 100, 0.0, 25000.0, VarManager::kMultFV0A);
        hm->AddHistogram(histClass, "MultFT0A", "MultFT0A", false, 100, 0.0, 25000.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "MultFT0C", "MultFT0C", false, 100, 0.0, 25000.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultFDDA", "MultFDDA", false, 100, 0.0, 25000.0, VarManager::kMultFDDA);
        hm->AddHistogram(histClass, "MultFDDC", "MultFDDC", false, 100, 0.0, 25000.0, VarManager::kMultFDDC);
        hm->AddHistogram(histClass, "MultZNA", "MultZNA", false, 100, 0.0, 25000.0, VarManager::kMultZNA);
        hm->AddHistogram(histClass, "MultZNC", "MultZNC", false, 100, 0.0, 25000.0, VarManager::kMultZNC);
        hm->AddHistogram(histClass, "MultTracklets", "MultTracklets", false, 100, 0.0, 25000.0, VarManager::kMultTracklets);
        hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 10000.0, VarManager::kVtxNcontribReal);
        hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 10000.0, VarManager::kVtxNcontrib);
        hm->AddHistogram(histClass, "MultTPC_MultFV0A", "MultTPC vs MultFV0A", false, 100, 0, 25000.0, VarManager::kMultTPC, 100, 0, 25000.0, VarManager::kMultFV0A);
        hm->AddHistogram(histClass, "MultTPC_MultFT0A", "MultTPC vs MultFT0A", false, 100, 0, 25000.0, VarManager::kMultTPC, 100, 0, 25000.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "MultTPC_MultFT0C", "MultTPC vs MultFT0C", false, 100, 0, 25000.0, VarManager::kMultTPC, 100, 0, 25000.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultFT0A_MultFT0C", "MultFT0A vs MultFT0C", false, 100, 0, 25000.0, VarManager::kMultFT0A, 100, 0, 25000.0, VarManager::kMultFT0C);
      }
    }
    if (subGroupStr.Contains("ftmulpbpb")) {
      hm->AddHistogram(histClass, "MultTPC", "MultTPC", false, 100, 0.0, 50000.0, VarManager::kMultTPC);
      hm->AddHistogram(histClass, "MultFT0C", "MultFT0C", false, 100, 0.0, 60000.0, VarManager::kMultFT0C);
      hm->AddHistogram(histClass, "MultFT0A", "MultFT0A", false, 100, 0.0, 180000.0, VarManager::kMultFT0A);
      hm->AddHistogram(histClass, "VtxNContribReal", "Vtx n contributors", false, 100, 0.0, 10000.0, VarManager::kVtxNcontribReal);
      hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 10000.0, VarManager::kVtxNcontrib);
      hm->AddHistogram(histClass, "MultTPC_VtxNContrib", "MultTPC vs VtxNContrib", false, 100, 0, 50000.0, VarManager::kMultTPC, 100, 0, 10000.0, VarManager::kVtxNcontrib);
      hm->AddHistogram(histClass, "MultFT0C_VtxNContrib", "MultFT0C vs VtxNContrib", false, 100, 0, 60000.0, VarManager::kMultFT0C, 100, 0, 10000.0, VarManager::kVtxNcontrib);
      hm->AddHistogram(histClass, "MultFT0A_VtxNContrib", "MultFT0A vs VtxNContrib", false, 100, 0, 180000.0, VarManager::kMultFT0A, 100, 0, 10000.0, VarManager::kVtxNcontrib);
    }
    if (subGroupStr.Contains("mc")) {
      hm->AddHistogram(histClass, "MCVtxX_VtxX", "Vtx X (MC vs rec)", false, 100, -0.5, 0.5, VarManager::kVtxX, 100, -0.5, 0.5, VarManager::kMCVtxX);
      hm->AddHistogram(histClass, "MCVtxY_VtxY", "Vtx Y (MC vs rec)", false, 100, -0.5, 0.5, VarManager::kVtxY, 100, -0.5, 0.5, VarManager::kMCVtxY);
      hm->AddHistogram(histClass, "MCVtxZ_VtxZ", "Vtx Z (MC vs rec)", false, 75, -15.0, 15.0, VarManager::kVtxZ, 75, -15.0, 15.0, VarManager::kMCVtxZ);
      hm->AddHistogram(histClass, "MCVtxZ", "Vtx Z (MC)", false, 75, -15.0, 15.0, VarManager::kMCVtxZ);
      hm->AddHistogram(histClass, "MCImpPar_CentVZERO", "MC impact param vs CentVZERO", false, 50, 0.0, 100.0, VarManager::kCentVZERO, 20, 0.0, 20.0, VarManager::kMCEventImpParam);
    }
    if (subGroupStr.Contains("generator")) {
      hm->AddHistogram(histClass, "MCVtxX", "Vtx X", false, 1000, -0.5, 0.5, VarManager::kMCVtxX);
      hm->AddHistogram(histClass, "MCVtxY", "Vtx Y", false, 1000, -0.5, 0.5, VarManager::kMCVtxY);
      hm->AddHistogram(histClass, "MCVtxX_VtxY", "Vtx X vs Vtx Y", false, 200, -0.2, 0.2, VarManager::kMCVtxX, 200, -0.2, 0.2, VarManager::kMCVtxY);
      hm->AddHistogram(histClass, "MCVtxZ", "Vtx Z", false, 60, -15.0, 15.0, VarManager::kMCVtxZ);
      hm->AddHistogram(histClass, "MCVtxZ_VtxX", "Vtx X vs Vtx Z", false, 60, -15.0, 15.0, VarManager::kMCVtxZ, 200, -0.2, 0.2, VarManager::kMCVtxX);
      hm->AddHistogram(histClass, "MCVtxX_VtxY", "Vtx X vs Vtx Y", false, 200, 15.0, 15.0, VarManager::kMCVtxZ, 200, -0.2, 0.2, VarManager::kMCVtxY);
      hm->AddHistogram(histClass, "MCImpPar", "MC impact param", false, 20, 0.0, 20.0, VarManager::kMCEventImpParam);
    }
    if (subGroupStr.Contains("qvector")) {
      int varZNA[3] = {VarManager::kQ1ZNAX, VarManager::kQ1ZNAY, VarManager::kCentFT0C};
      int varZNC[3] = {VarManager::kQ1ZNCX, VarManager::kQ1ZNCY, VarManager::kCentFT0C};

      int bins[3] = {500, 500, 18};
      double minBins[3] = {-10, -10, 0};
      double maxBins[3] = {10, 10, 90};
      hm->AddHistogram(histClass, "Q1ZNAX_Q1ZNAY_CentFT0C", "", 3, varZNA, bins, minBins, maxBins, 0, -1, kTRUE);
      hm->AddHistogram(histClass, "Q1ZNCX_Q1ZNCY_CentFT0C", "", 3, varZNC, bins, minBins, maxBins, 0, -1, kTRUE);

      hm->AddHistogram(histClass, "IntercalibZNA_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -50.0, 50.0, VarManager::KIntercalibZNA);
      hm->AddHistogram(histClass, "IntercalibZNC_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -50.0, 50.0, VarManager::KIntercalibZNC);

      hm->AddHistogram(histClass, "Q2X0A", "", false, 500, -10.0, 10.0, VarManager::kQ2X0A);
      hm->AddHistogram(histClass, "Q2Y0A", "", false, 500, -10.0, 10.0, VarManager::kQ2Y0A);
      hm->AddHistogram(histClass, "Q2X0B", "", false, 500, -10.0, 10.0, VarManager::kQ2X0B);
      hm->AddHistogram(histClass, "Q2Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ2Y0B);
      hm->AddHistogram(histClass, "Q2X0C", "", false, 500, -10.0, 10.0, VarManager::kQ2X0C);
      hm->AddHistogram(histClass, "Q2Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ2Y0C);
      hm->AddHistogram(histClass, "Q2X0A_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2X0A);
      hm->AddHistogram(histClass, "Q2Y0A_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2Y0A);
      hm->AddHistogram(histClass, "Q2X0B_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2X0B);
      hm->AddHistogram(histClass, "Q2Y0B_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2Y0B);
      hm->AddHistogram(histClass, "Q2X0C_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2X0C);
      hm->AddHistogram(histClass, "Q2Y0C_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2Y0C);
      hm->AddHistogram(histClass, "Q2X0A_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2X0A);
      hm->AddHistogram(histClass, "Q2Y0A_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2Y0A);
      hm->AddHistogram(histClass, "Q2X0B_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2X0B);
      hm->AddHistogram(histClass, "Q2Y0B_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2Y0B);
      hm->AddHistogram(histClass, "Q2X0C_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2X0C);
      hm->AddHistogram(histClass, "Q2Y0C_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2Y0C);
      hm->AddHistogram(histClass, "Q3X0A", "", false, 500, -10.0, 10.0, VarManager::kQ3X0A);
      hm->AddHistogram(histClass, "Q3Y0A", "", false, 500, -10.0, 10.0, VarManager::kQ3Y0A);
      hm->AddHistogram(histClass, "Q3X0B", "", false, 500, -10.0, 10.0, VarManager::kQ3X0B);
      hm->AddHistogram(histClass, "Q3Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ3Y0B);
      hm->AddHistogram(histClass, "Q3X0C", "", false, 500, -10.0, 10.0, VarManager::kQ3X0C);
      hm->AddHistogram(histClass, "Q3Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ3Y0C);
      hm->AddHistogram(histClass, "Q3X0A_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3X0A);
      hm->AddHistogram(histClass, "Q3Y0A_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3Y0A);
      hm->AddHistogram(histClass, "Q3X0B_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3X0B);
      hm->AddHistogram(histClass, "Q3Y0B_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3Y0B);
      hm->AddHistogram(histClass, "Q3X0C_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3X0C);
      hm->AddHistogram(histClass, "Q3Y0C_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3Y0C);
      hm->AddHistogram(histClass, "Q3X0A_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3X0A);
      hm->AddHistogram(histClass, "Q3Y0A_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3Y0A);
      hm->AddHistogram(histClass, "Q3X0B_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3X0B);
      hm->AddHistogram(histClass, "Q3Y0B_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3Y0B);
      hm->AddHistogram(histClass, "Q3X0C_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3X0C);
      hm->AddHistogram(histClass, "Q3Y0C_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3Y0C);
      hm->AddHistogram(histClass, "Q2X0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2X0A);
      hm->AddHistogram(histClass, "Q2Y0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2Y0A);
      hm->AddHistogram(histClass, "Q2X0B_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2X0B);
      hm->AddHistogram(histClass, "Q2Y0B_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2Y0B);
      hm->AddHistogram(histClass, "Q2X0C_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2X0C);
      hm->AddHistogram(histClass, "Q2Y0C_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2Y0C);
      hm->AddHistogram(histClass, "Q3X0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ3X0A);
      hm->AddHistogram(histClass, "Q3Y0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ3Y0A);
      hm->AddHistogram(histClass, "Psi2A", "", false, 100, -2.0, 2.0, VarManager::kPsi2A);
      hm->AddHistogram(histClass, "Psi2B", "", false, 100, -2.0, 2.0, VarManager::kPsi2B);
      hm->AddHistogram(histClass, "Psi2C", "", false, 100, -2.0, 2.0, VarManager::kPsi2C);
      hm->AddHistogram(histClass, "Psi2A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 100, -2.0, 2.0, VarManager::kPsi2A);
      hm->AddHistogram(histClass, "Psi2B_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 100, -2.0, 2.0, VarManager::kPsi2B);
      hm->AddHistogram(histClass, "Psi2C_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 100, -2.0, 2.0, VarManager::kPsi2C);
      if (subGroupStr.Contains("cross")) {
        hm->AddHistogram(histClass, "Q1ZNACXX_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 4000, -2, 2, VarManager::kQ1ZNACXX);
        hm->AddHistogram(histClass, "Q1ZNACYY_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 4000, -2, 2, VarManager::kQ1ZNACYY);
        hm->AddHistogram(histClass, "Q1ZNACYX_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 4000, -2, 2, VarManager::kQ1ZNACYX);
        hm->AddHistogram(histClass, "Q1ZNACXY_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 4000, -2, 2, VarManager::kQ1ZNACXY);
        hm->AddHistogram(histClass, "Q2X0A_Q2Y0A", "", false, 500, -10.0, 10.0, VarManager::kQ2X0A, 500, -10.0, 10.0, VarManager::kQ2Y0A);
        hm->AddHistogram(histClass, "Q2X0B_Q2Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ2X0B, 500, -10.0, 10.0, VarManager::kQ2Y0B);
        hm->AddHistogram(histClass, "Q2X0C_Q2Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ2X0C, 500, -10.0, 10.0, VarManager::kQ2Y0C);
        hm->AddHistogram(histClass, "Q2X0B_Q2Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ2X0B, 500, -10.0, 10.0, VarManager::kQ2Y0C);
        hm->AddHistogram(histClass, "Q2X0C_Q2Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ2X0C, 500, -10.0, 10.0, VarManager::kQ2Y0B);
        hm->AddHistogram(histClass, "Q3X0B_Q3Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ3X0B, 500, -10.0, 10.0, VarManager::kQ3Y0C);
        hm->AddHistogram(histClass, "Q3X0C_Q3Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ3X0C, 500, -10.0, 10.0, VarManager::kQ3Y0B);
        hm->AddHistogram(histClass, "Q2YYAB_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YYAB);
        hm->AddHistogram(histClass, "Q2XXAB_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XXAB);
        hm->AddHistogram(histClass, "Q2XYAB_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XYAB);
        hm->AddHistogram(histClass, "Q2YXAB_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YXAB);
        hm->AddHistogram(histClass, "Q2YYAC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YYAC);
        hm->AddHistogram(histClass, "Q2XXAC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XXAC);
        hm->AddHistogram(histClass, "Q2XYAC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XYAC);
        hm->AddHistogram(histClass, "Q2YXAC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YXAC);
        hm->AddHistogram(histClass, "Q2YYBC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YYBC);
        hm->AddHistogram(histClass, "Q2XXBC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XXBC);
        hm->AddHistogram(histClass, "Q2XYBC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XYBC);
        hm->AddHistogram(histClass, "Q2YXBC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YXBC);
        hm->AddHistogram(histClass, "Q2YYAB_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YYAB);
        hm->AddHistogram(histClass, "Q2XXAB_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XXAB);
        hm->AddHistogram(histClass, "Q2XYAB_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XYAB);
        hm->AddHistogram(histClass, "Q2YXAB_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YXAB);
        hm->AddHistogram(histClass, "Q2YYAC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YYAC);
        hm->AddHistogram(histClass, "Q2XXAC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XXAC);
        hm->AddHistogram(histClass, "Q2XYAC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XYAC);
        hm->AddHistogram(histClass, "Q2YXAC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YXAC);
        hm->AddHistogram(histClass, "Q2YYBC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YYBC);
        hm->AddHistogram(histClass, "Q2XXBC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XXBC);
        hm->AddHistogram(histClass, "Q2XYBC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XYBC);
        hm->AddHistogram(histClass, "Q2YXBC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YXBC);
      }
    }
    if (subGroupStr.Contains("res")) {
      hm->AddHistogram(histClass, "R2SP_TPCFT0A_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2SP_AB);
      hm->AddHistogram(histClass, "R2SP_TPCFT0C_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2SP_AC);
      hm->AddHistogram(histClass, "R2SP_FT0AFT0C_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2SP_BC);
      hm->AddHistogram(histClass, "R3SP_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR3SP);
      hm->AddHistogram(histClass, "R3EP_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR3EP);
      hm->AddHistogram(histClass, "R2EP_TPCFT0A_CentV0M", "", false, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2EP_AB);
      hm->AddHistogram(histClass, "R2EP_TPCFT0C_CentV0M", "", false, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2EP_AC);
      hm->AddHistogram(histClass, "R2EP_FT0CFT0A_CentV0M", "", false, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2EP_BC);
      hm->AddHistogram(histClass, "R2SP_TPCFT0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_AB);
      hm->AddHistogram(histClass, "R2SP_TPCFT0C_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_AC);
      hm->AddHistogram(histClass, "R2SP_FT0AFT0C_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_BC);
      hm->AddHistogram(histClass, "R2SP_FT0CTPCPOS_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_FT0CTPCPOS);
      hm->AddHistogram(histClass, "R2SP_FT0CTPCNEG_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_FT0CTPCNEG);
      hm->AddHistogram(histClass, "R2SP_FT0ATPCPOS_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_FT0ATPCPOS);
      hm->AddHistogram(histClass, "R2SP_FT0ATPCNEG_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_FT0ATPCNEG);
      hm->AddHistogram(histClass, "R3SP_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR3SP);
      hm->AddHistogram(histClass, "R2EP_TPCFT0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_AB);
      hm->AddHistogram(histClass, "R2EP_TPCFT0C_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_AC);
      hm->AddHistogram(histClass, "R2EP_FT0CFT0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_BC);
      hm->AddHistogram(histClass, "R2EP_FT0CTPCPOS_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_FT0CTPCPOS);
      hm->AddHistogram(histClass, "R2EP_FT0CTPCNEG_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_FT0CTPCNEG);
      hm->AddHistogram(histClass, "R2EP_FT0ATPCPOS_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_FT0ATPCPOS);
      hm->AddHistogram(histClass, "R2EP_FT0ATPCNEG_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_FT0ATPCNEG);
      hm->AddHistogram(histClass, "R3EP_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR3EP);
    }
    if (subGroupStr.Contains("filter")) {
      hm->AddHistogram(histClass, "IsDoubleGap", "Is double gap", false, 2, -0.5, 1.5, VarManager::kIsDoubleGap);
      hm->AddHistogram(histClass, "IsSingleGapA", "Is single gap on side A", false, 2, -0.5, 1.5, VarManager::kIsSingleGapA);
      hm->AddHistogram(histClass, "IsSingleGapC", "Is single gap on side C", false, 2, -0.5, 1.5, VarManager::kIsSingleGapC);
    }
    if (subGroupStr.Contains("zdc")) {
      hm->AddHistogram(histClass, "energyCommonZNA_energyCommonZNC", "Common ZNA energy vs common ZNC energy", false, 1050, -10.0, 200.0, VarManager::kEnergyCommonZNA, 1050, -10.0, 200.0, VarManager::kEnergyCommonZNC);
      hm->AddHistogram(histClass, "energyCommonZNA_energyCommonZNC_lowRange", "Common ZNA energy vs common ZNC energy", false, 220, -2.0, 20.0, VarManager::kEnergyCommonZNA, 220, -2.0, 20.0, VarManager::kEnergyCommonZNC);
    }
  } // end "event"

  if (!groupStr.CompareTo("two-collisions")) {
    hm->AddHistogram(histClass, "DeltaZ", "z_{1} - z_{2}", false, 400, -20., 20., VarManager::kTwoEvDeltaZ);
    hm->AddHistogram(histClass, "DeltaZ_Z1", "z_{1} - z_{2} vs z_{1}", false, 24, -12., 12., VarManager::kTwoEvPosZ1, 300, -15., 15., VarManager::kTwoEvDeltaZ);
    hm->AddHistogram(histClass, "DeltaR", "r_{1} - r_{2}", false, 200, -0.1, 0.1, VarManager::kTwoEvDeltaR);
    hm->AddHistogram(histClass, "DeltaX", "x_{1} - x_{2}", false, 200, -0.1, 0.1, VarManager::kTwoEvDeltaX);
    hm->AddHistogram(histClass, "DeltaY", "y_{1} - y_{2}", false, 200, -0.1, 0.1, VarManager::kTwoEvDeltaY);
    hm->AddHistogram(histClass, "Z1vsZ2", "z_{1} vs z_{2}", false, 120, -12.0, 12.0, VarManager::kTwoEvPosZ1, 120, -12.0, 12.0, VarManager::kTwoEvPosZ2);
    hm->AddHistogram(histClass, "R1vsR2", "r_{1} vs r_{2}", false, 100, -0.1, 0.1, VarManager::kTwoEvPosR1, 100, -0.1, 0.1, VarManager::kTwoEvPosR2);
    hm->AddHistogram(histClass, "NContrib1vs2", "n.contrib 1 vs 2", false, 100, 0.0, 100.0, VarManager::kTwoEvPVcontrib1, 100, 0.0, 100.0, VarManager::kTwoEvPVcontrib2);
  }

  if (!groupStr.CompareTo("track")) {
    hm->AddHistogram(histClass, "Pt", "p_{T} distribution", false, 2000, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Eta", "#eta distribution", false, 500, -5.0, 5.0, VarManager::kEta);
    hm->AddHistogram(histClass, "Phi", "#varphi distribution", false, 500, -2. * TMath::Pi(), 2. * TMath::Pi(), VarManager::kPhi);
    hm->AddHistogram(histClass, "Phi_Pt", "#varphi distribution", false, 50, 0.0, 10.0, VarManager::kPt, 720, 0.0, TMath::TwoPi(), VarManager::kPhi);
    hm->AddHistogram(histClass, "IsPVcontrib_pt", "is PV contributor vs pt", false, 50, 0.0, 50.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kPVContributor);
    hm->AddHistogram(histClass, "IsPVcontrib_pt_prof", "is PV contributor vs pt", true, 50, 0.0, 50.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kPVContributor);
    if (subGroupStr.Contains("ambiguity") && !subGroupStr.Contains("muon")) {
      hm->AddHistogram(histClass, "AmbiguityInBunch", "in bunch collision ambiguity", false, 10, 0., 10., VarManager::kBarrelNAssocsInBunch);
      hm->AddHistogram(histClass, "AmbiguityInBunch_pt", "in bunch collision ambiguity vs p_{T}", false, 50, 0.0, 10.0, VarManager::kPt, 10, 0., 10., VarManager::kBarrelNAssocsInBunch);
      hm->AddHistogram(histClass, "AmbiguityOutOfBunch", "out of bunch collision ambiguity", false, 10, 0., 10., VarManager::kBarrelNAssocsOutOfBunch);
      hm->AddHistogram(histClass, "AmbiguityOutOfBunch_pt", "out of bunch collision ambiguity vs p_{T}", false, 50, 0.0, 10.0, VarManager::kPt, 10, 0., 10., VarManager::kBarrelNAssocsOutOfBunch);
    }
    if (subGroupStr.Contains("cent")) {
      hm->AddHistogram(histClass, "Pt_CentFT0C", "p_{T} distribution", false, 2000, 0.0, 20.0, VarManager::kPt, 20, 0.0, 100.0, VarManager::kCentFT0C);
      hm->AddHistogram(histClass, "Eta_CentFT0C", "#eta distribution", false, 500, -5.0, 5.0, VarManager::kEta, 20, 0.0, 100.0, VarManager::kCentFT0C);
      hm->AddHistogram(histClass, "Phi_CentFT0C", "#varphi distribution", false, 500, -2. * TMath::Pi(), 2. * TMath::Pi(), VarManager::kPhi, 20, 0.0, 100.0, VarManager::kCentFT0C);
    }
    if (subGroupStr.Contains("kine")) {
      hm->AddHistogram(histClass, "Phi_Eta", "#phi vs #eta distribution", false, 200, -5.0, 5.0, VarManager::kEta, 200, -2. * TMath::Pi(), 2. * TMath::Pi(), VarManager::kPhi);
      hm->AddHistogram(histClass, "Eta_Pt", "", false, 20, -1.0, 1.0, VarManager::kEta, 100, 0.0, 20.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Eta_VtxZ", "", false, 100, -1.0, 1.0, VarManager::kEta, 300, -15.0, 15.0, VarManager::kVtxZ);
      hm->AddHistogram(histClass, "Px", "p_{x} distribution", false, 200, 0.0, 20.0, VarManager::kPx);
      hm->AddHistogram(histClass, "Py", "p_{y} distribution", false, 200, 0.0, 20.0, VarManager::kPy);
      hm->AddHistogram(histClass, "Pz", "p_{z} distribution", false, 200, 0.0, 20.0, VarManager::kPz);
      hm->AddHistogram(histClass, "InvPt", "1/p_{T} distribution", false, 1500, 0.0, 15.0, VarManager::kInvPt);
    }
    if (subGroupStr.Contains("time")) {
      hm->AddHistogram(histClass, "TrackTime", "", false, 800, -10000.0, 10000.0, VarManager::kTrackTime);
      hm->AddHistogram(histClass, "TrackTimeRes", "", false, 400, 0.0, 1000.0, VarManager::kTrackTimeRes);
      hm->AddHistogram(histClass, "TrackTimeResRelative", "", false, 100, 0.0, 100.0, VarManager::kTrackTimeResRelative);
      hm->AddHistogram(histClass, "TOFExpMom", "", false, 100, 0.0, 10.0, VarManager::kTOFExpMom);
    }
    if (subGroupStr.Contains("map")) {
      hm->AddHistogram(histClass, "HasITSandTPC", "", false, 2, -0.5, 1.5, VarManager::kHasITS, 2, -0.5, 1.5, VarManager::kHasTPC);
      hm->AddHistogram(histClass, "HasITSandTOF", "", false, 2, -0.5, 1.5, VarManager::kHasITS, 2, -0.5, 1.5, VarManager::kHasTOF);
      hm->AddHistogram(histClass, "HasITSandTRD", "", false, 2, -0.5, 1.5, VarManager::kHasITS, 2, -0.5, 1.5, VarManager::kHasTRD);
      hm->AddHistogram(histClass, "HasTPCandTOF", "", false, 2, -0.5, 1.5, VarManager::kHasTPC, 2, -0.5, 1.5, VarManager::kHasTOF);
    }
    if (subGroupStr.Contains("its")) {
      hm->AddHistogram(histClass, "ITSncls", "Number of cluster in ITS", false, 8, -0.5, 7.5, VarManager::kITSncls);
      hm->AddHistogram(histClass, "ITSncls_vs_BC", "Average Number of cluster in ITS vs BC", true, 3564, 0.0, 3564, VarManager::kBCOrbit, 2, 0.0, 7.0, VarManager::kITSncls);
      hm->AddHistogram(histClass, "ITSchi2", "ITS chi2", false, 100, 0.0, 50.0, VarManager::kITSchi2);
      hm->AddHistogram(histClass, "IsITSrefit", "", false, 2, -0.5, 1.5, VarManager::kIsITSrefit);
      hm->AddHistogram(histClass, "IsSPDany", "", false, 2, -0.5, 1.5, VarManager::kIsSPDany);
      hm->AddHistogram(histClass, "IsSPDfirst", "", false, 2, -0.5, 1.5, VarManager::kIsSPDfirst);
      if (subGroupStr.Contains("cluster")) {
        hm->AddHistogram(histClass, "ITSClusterMap", "", false, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_pin", "ITSClustermap vs pin", false, 200, 0.0, 20.0, VarManager::kPin, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_SignedPin", "ITSClustermap vs SignedPin", false, 400, -20.0, 20.0, VarManager::kSignedPin, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_pt", "ITSClustermap vs pt", false, 200, 0.0, 20.0, VarManager::kPt, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_eta", "ITSClustermap vs eta", false, 100, -1.0, 1.0, VarManager::kEta, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_phi", "ITSClustermap vs phi", false, 315, 0.0, 6.3, VarManager::kPhi, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "SignedPin_P_ITSMap", "SignedPin vs P vs ITSMap", false, 400, -20.0, 20.0, VarManager::kSignedPin, 200, 0.0, 20.0, VarManager::kP, 2, -0.5, 1.5, VarManager::kHasITS);
      }
      if (subGroupStr.Contains("cent")) {
        hm->AddHistogram(histClass, "ITSncls_CentFT0C", "Number of cluster in ITS", false, 8, -0.5, 7.5, VarManager::kITSncls, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "ITSchi2_CentFT0C", "ITS chi2", false, 100, 0.0, 50.0, VarManager::kITSchi2, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "ITSClusterMap_CentFT0C", "", false, 128, -0.5, 127.5, VarManager::kITSClusterMap, 20, 0.0, 100.0, VarManager::kCentFT0C);
      }
      if (subGroupStr.Contains("clssize")) {
        hm->AddHistogram(histClass, "ITSclssize_p", "Mean ITS cluster size vs P", false, 200, 0.0, 10.0, VarManager::kP, 150, 0.0, 15., VarManager::kITSmeanClsSize);
      }
    }
    if (subGroupStr.Contains("itsvspt")) {
      hm->AddHistogram(histClass, "ITSncls_Pt", "Number of cluster in ITS vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 8, -0.5, 7.5, VarManager::kITSncls);
      hm->AddHistogram(histClass, "ITSchi2_Pt", "ITS chi2 vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 100, 0.0, 50.0, VarManager::kITSchi2);
      hm->AddHistogram(histClass, "IsITSrefit_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsITSrefit);
      hm->AddHistogram(histClass, "IsSPDany_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsSPDany);
      hm->AddHistogram(histClass, "IsSPDfirst_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsSPDfirst);
      hm->AddHistogram(histClass, "ITSncls_signedPt", "Number of cluster in ITS vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 8, -0.5, 7.5, VarManager::kITSncls);
      hm->AddHistogram(histClass, "ITSchi2_signedPt", "ITS chi2 vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 100, 0.0, 50.0, VarManager::kITSchi2);
      hm->AddHistogram(histClass, "IsITSrefit_signedPt", "", false, 400, -10.0, 10.0, VarManager::kSignedPt, 2, -0.5, 1.5, VarManager::kIsITSrefit);
    }
    if (subGroupStr.Contains("tpc")) {
      hm->AddHistogram(histClass, "TPCncls", "Number of cluster in TPC", false, 160, -0.5, 159.5, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCncls_Phi", "Number of cluster in TPC vs #varphi", true, 720, 0.0, TMath::TwoPi(), VarManager::kPhi, 10, 0.0, 200.0, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCncls_PhiPt", "Number of cluster in TPC vs p_{T} and #varphi", true, 20, 0.0, 10.0, VarManager::kPt, 720, 0.0, TMath::TwoPi(), VarManager::kPhi, 10, 0.0, 200.0, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCncls_Run", "Number of cluster in TPC", true, (VarManager::GetNRuns() > 0 ? VarManager::GetNRuns() : 1), 0.5, 0.5 + VarManager::GetNRuns(), VarManager::kRunId,
                       10, -0.5, 159.5, VarManager::kTPCncls, 10, 0., 1., VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "TPCnclsCR", "Number of crossed rows in TPC", false, 160, -0.5, 159.5, VarManager::kTPCnclsCR);
      hm->AddHistogram(histClass, "TPCncls_TPCnclsCR", "Number of TPC cluster vs Number of crossed rows in TPC", false, 160, -0.5, 159.5, VarManager::kTPCncls, 160, -0.5, 159.5, VarManager::kTPCnclsCR);
      hm->AddHistogram(histClass, "IsTPCrefit", "", false, 2, -0.5, 1.5, VarManager::kIsTPCrefit);
      hm->AddHistogram(histClass, "IsGoldenChi2", "", false, 2, -0.5, 1.5, VarManager::kIsGoldenChi2);
      hm->AddHistogram(histClass, "TPCchi2", "TPC chi2", false, 100, 0.0, 10.0, VarManager::kTPCchi2);
      hm->AddHistogram(histClass, "pin_vs_p", "", false, 400, -20.0, 20.0, VarManager::kSignedPin, 200, 0.0, 20, VarManager::kP);
      if (subGroupStr.Contains("cent")) {
        hm->AddHistogram(histClass, "TPCncls_CentFT0C", "Number of cluster in TPC", false, 160, -0.5, 159.5, VarManager::kTPCncls, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "TPCnclsCR_CentFT0C", "Number of crossed rows in TPC", false, 160, -0.5, 159.5, VarManager::kTPCnclsCR, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "TPCchi2_CentFT0C", "TPC chi2", false, 100, 0.0, 10.0, VarManager::kTPCchi2, 20, 0.0, 100.0, VarManager::kCentFT0C);
      }
    }
    if (subGroupStr.Contains("tpcvspt")) {
      hm->AddHistogram(histClass, "TPCncls_Pt", "Number of cluster in TPC vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 160, -0.5, 159.5, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCnclsCR_Pt", "Number of crossed rows in TPC vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 160, -0.5, 159.5, VarManager::kTPCnclsCR);
      hm->AddHistogram(histClass, "IsTPCrefit_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsTPCrefit);
      hm->AddHistogram(histClass, "IsGoldenChi2_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsGoldenChi2);
      hm->AddHistogram(histClass, "TPCchi2_Pt", "TPC chi2 vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 100, 0.0, 10.0, VarManager::kTPCchi2);
      hm->AddHistogram(histClass, "TPCncls_signedPt", "Number of cluster in TPC vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 160, -0.5, 159.5, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCnclsCR_signedPt", "Number of crossed rows in TPC vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 160, -0.5, 159.5, VarManager::kTPCnclsCR);
      hm->AddHistogram(histClass, "TPCchi2_signedPt", "TPC chi2 vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 100, 0.0, 10.0, VarManager::kTPCchi2);
    }
    if (subGroupStr.Contains("tpcpid")) {
      if (subGroupStr.Contains("tpcpid_fine")) {
        // fine binning for pIN: steps in 10 MeV/c from 0 to 1 GeV/c and 100 MeV/c up to 10 GeV/c
        double pIN_bins[281];
        for (int i = 0; i <= 200; i++)
          pIN_bins[i] = 0.01 * i;
        for (int i = 1; i <= 80; i++)
          pIN_bins[200 + i] = 2. + 0.1 * i;
        int nbins_pIN = sizeof(pIN_bins) / sizeof(*pIN_bins) - 1;

        double TPCdEdx_bins[201];
        for (int i = 0; i <= 200; i++)
          TPCdEdx_bins[i] = i;
        int nbins_TPCdEdx = sizeof(TPCdEdx_bins) / sizeof(*TPCdEdx_bins) - 1;

        double nSigma_bins[101];
        for (int i = 0; i <= 100; i++)
          nSigma_bins[i] = -5. + 0.1 * i;
        int nbins_nSigma = sizeof(nSigma_bins) / sizeof(*nSigma_bins) - 1;

        hm->AddHistogram(histClass, "TPCdedx_pIN", "TPC dE/dx vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_TPCdEdx, TPCdEdx_bins, VarManager::kTPCsignal);
        hm->AddHistogram(histClass, "TPCnSigEle_pIN", "TPC n-#sigma(e) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigPi_pIN", "TPC n-#sigma(#pi) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigKa_pIN", "TPC n-#sigma(K) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaKa);
        hm->AddHistogram(histClass, "TPCnSigPr_pIN", "TPC n-#sigma(p) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaPr);
        hm->AddHistogram(histClass, "TPCnSigEl_Corr_pIN", "TPC n-#sigma(e) Corr. vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaEl_Corr);
        hm->AddHistogram(histClass, "TPCnSigPi_Corr_pIN", "TPC n-#sigma(#pi) Corr. vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaPi_Corr);
        hm->AddHistogram(histClass, "TPCnSigPr_Corr_pIN", "TPC n-#sigma(p) Corr. vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaPr_Corr);
      } else {
        hm->AddHistogram(histClass, "TPCdedx_pIN", "TPC dE/dx vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 200, 0.0, 200., VarManager::kTPCsignal);
        hm->AddHistogram(histClass, "TPCnSigEle_pIN", "TPC n-#sigma(e) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigPi_pIN", "TPC n-#sigma(#pi) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigKa_pIN", "TPC n-#sigma(K) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaKa);
        hm->AddHistogram(histClass, "TPCnSigPr_pIN", "TPC n-#sigma(p) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaPr);
        hm->AddHistogram(histClass, "TPCnSigEl_Corr_pIN", "TPC n-#sigma(e) Corr. vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaEl_Corr);
        hm->AddHistogram(histClass, "TPCnSigPi_Corr_pIN", "TPC n-#sigma(#pi) Corr. vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaPi_Corr);
        hm->AddHistogram(histClass, "TPCnSigKa_Corr_pIN", "TPC n-#sigma(K) Corr. vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaKa_Corr);
        hm->AddHistogram(histClass, "TPCnSigPr_Corr_pIN", "TPC n-#sigma(p) Corr. vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaPr_Corr);
      }
      hm->AddHistogram(histClass, "TPCnSigEl_Corr_Eta", "TPC n-#sigma(e) Corr. vs Eta", false, 20, -1.0, 1.0, VarManager::kEta, 100, -5.0, 5.0, VarManager::kTPCnSigmaEl_Corr);
      hm->AddHistogram(histClass, "TPCnSigPi_Corr_Eta", "TPC n-#sigma(#pi) Corr. vs Eta", false, 20, -1.0, 1.0, VarManager::kEta, 100, -5.0, 5.0, VarManager::kTPCnSigmaPi_Corr);
      hm->AddHistogram(histClass, "TPCnSigKa_Corr_Eta", "TPC n-#sigma(K) Corr. vs Eta", false, 20, -1.0, 1.0, VarManager::kEta, 100, -5.0, 5.0, VarManager::kTPCnSigmaKa_Corr);
      hm->AddHistogram(histClass, "TPCnSigPr_Corr_Eta", "TPC n-#sigma(p) Corr. vs Eta", false, 20, -1.0, 1.0, VarManager::kEta, 100, -5.0, 5.0, VarManager::kTPCnSigmaPr_Corr);
    }
    if (subGroupStr.Contains("postcalib")) {
      const int kNvarsPID = 4;
      const int kTPCnsigmaNbins = 70;
      double tpcNsigmaBinLims[kTPCnsigmaNbins + 1];
      for (int i = 0; i <= kTPCnsigmaNbins; ++i)
        tpcNsigmaBinLims[i] = -7.0 + 0.2 * i;

      const int kPinEleNbins = 18;
      double pinEleBinLims[kPinEleNbins + 1] = {0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0};

      const int kEtaNbins = 9;
      double etaBinLimsI[kEtaNbins + 1] = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9};

      const int kTPCnClusterbins = 16;
      double tpcNclusterBinLims[kTPCnClusterbins + 1];
      for (int i = 0; i <= kTPCnClusterbins; ++i)
        tpcNclusterBinLims[i] = 10 * i;

      TArrayD nSigBinLimits[kNvarsPID];
      nSigBinLimits[0] = TArrayD(kTPCnsigmaNbins + 1, tpcNsigmaBinLims);
      nSigBinLimits[1] = TArrayD(kTPCnClusterbins + 1, tpcNclusterBinLims);
      nSigBinLimits[2] = TArrayD(kPinEleNbins + 1, pinEleBinLims);
      nSigBinLimits[3] = TArrayD(kEtaNbins + 1, etaBinLimsI);

      const int kPinKaNbins = 15;
      double pinKaBinLims[kPinKaNbins + 1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0};

      TArrayD nSigBinLimitsKa[kNvarsPID];
      nSigBinLimitsKa[0] = TArrayD(kTPCnsigmaNbins + 1, tpcNsigmaBinLims);
      nSigBinLimitsKa[1] = TArrayD(kTPCnClusterbins + 1, tpcNclusterBinLims);
      nSigBinLimitsKa[2] = TArrayD(kPinKaNbins + 1, pinKaBinLims);
      nSigBinLimitsKa[3] = TArrayD(kEtaNbins + 1, etaBinLimsI);

      if (subGroupStr.Contains("electron")) {
        int varsPIDnSigEle[kNvarsPID] = {VarManager::kTPCnSigmaEl, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        int varsPIDnSigEle_Corr[kNvarsPID] = {VarManager::kTPCnSigmaEl_Corr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        hm->AddHistogram(histClass, "nSigmaTPCelectron", "TPC n_{#sigma}(e) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigEle, nSigBinLimits);
        hm->AddHistogram(histClass, "nSigmaTPCelectron_Corr", "TPC n_{#sigma}^{Corr}(e) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigEle_Corr, nSigBinLimits);
      }
      if (subGroupStr.Contains("pion")) {
        int varsPIDnSigPion[kNvarsPID] = {VarManager::kTPCnSigmaPi, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        int varsPIDnSigPion_Corr[kNvarsPID] = {VarManager::kTPCnSigmaPi_Corr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        hm->AddHistogram(histClass, "nSigmaTPCpion", "TPC n_{#sigma}(pion) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigPion, nSigBinLimits);
        hm->AddHistogram(histClass, "nSigmaTPCpion_Corr", "TPC n_{#sigma}^{Corr}(pion) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigPion_Corr, nSigBinLimits);
      }
      if (subGroupStr.Contains("kaon")) {
        int varsPIDnSigKaon[kNvarsPID] = {VarManager::kTPCnSigmaKa, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        int varsPIDnSigKaon_Corr[kNvarsPID] = {VarManager::kTPCnSigmaKa_Corr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        hm->AddHistogram(histClass, "nSigmaTPCkaon", "TPC n_{#sigma}(kaon) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigKaon, nSigBinLimitsKa);
        hm->AddHistogram(histClass, "nSigmaTPCkaon_Corr", "TPC n_{#sigma}^{Corr}(kaon) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigKaon_Corr, nSigBinLimitsKa);
      }
      if (subGroupStr.Contains("proton")) {
        int varsPIDnSigProton[kNvarsPID] = {VarManager::kTPCnSigmaPr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        int varsPIDnSigProton_Corr[kNvarsPID] = {VarManager::kTPCnSigmaPr_Corr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        hm->AddHistogram(histClass, "nSigmaTPCproton", "TPC n_{#sigma}(proton) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigProton, nSigBinLimits);
        hm->AddHistogram(histClass, "nSigmaTPCproton_Corr", "TPC n_{#sigma}^{Corr}(proton) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigProton_Corr, nSigBinLimits);
      }
    }
    if (subGroupStr.Contains("tofpid")) {
      if (subGroupStr.Contains("tofpid_fine")) {
        // fine binning for pIN: steps in 10 MeV/c from 0 to 1 GeV/c and 100 MeV/c up to 10 GeV/c
        double pIN_bins[281];
        for (int i = 0; i <= 200; i++)
          pIN_bins[i] = 0.01 * i;
        for (int i = 1; i <= 80; i++)
          pIN_bins[200 + i] = 2. + 0.1 * i;
        int nbins_pIN = sizeof(pIN_bins) / sizeof(*pIN_bins) - 1;

        double TOFbeta_bins[241];
        for (int i = 0; i <= 240; i++)
          TOFbeta_bins[i] = 0.005 * i;
        int nbins_TOFbeta = sizeof(TOFbeta_bins) / sizeof(*TOFbeta_bins) - 1;

        double nSigma_bins[101];
        for (int i = 0; i <= 100; i++)
          nSigma_bins[i] = -5. + 0.1 * i;
        int nbins_nSigma = sizeof(nSigma_bins) / sizeof(*nSigma_bins) - 1;

        hm->AddHistogram(histClass, "TOFbeta_pIN", "TOF #beta vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_TOFbeta, TOFbeta_bins, VarManager::kTOFbeta);
        hm->AddHistogram(histClass, "TOFnSigEle_pIN", "TOF n-#sigma(e) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTOFnSigmaEl);
        hm->AddHistogram(histClass, "TOFnSigPi_pIN", "TOF n-#sigma(#pi) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTOFnSigmaPi);
        hm->AddHistogram(histClass, "TOFnSigKa_pIN", "TOF n-#sigma(K) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTOFnSigmaKa);
        hm->AddHistogram(histClass, "TOFnSigPr_pIN", "TOF n-#sigma(p) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTOFnSigmaPr);
      } else {
        hm->AddHistogram(histClass, "TOFbeta_pIN", "TOF #beta vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 240, 0.0, 1.2, VarManager::kTOFbeta);
        hm->AddHistogram(histClass, "TOFnSigEle_pIN", "TOF n-#sigma(e) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTOFnSigmaEl);
        hm->AddHistogram(histClass, "TOFnSigPi_pIN", "TOF n-#sigma(#pi) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTOFnSigmaPi);
        hm->AddHistogram(histClass, "TOFnSigKa_pIN", "TOF n-#sigma(K) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTOFnSigmaKa);
        hm->AddHistogram(histClass, "TOFnSigPr_pIN", "TOF n-#sigma(p) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTOFnSigmaPr);
      }
    }
    if (subGroupStr.Contains("pidcorre")) {
      const int kNvarsPID = 3;
      const int kNbins_pIN = 169;
      double pIN_bins[kNbins_pIN + 1];
      for (int i = 0; i <= 140; i++)
        pIN_bins[i] = 0.01 * i + 0.1;
      for (int i = 1; i <= 15; i++)
        pIN_bins[140 + i] = 1.5 + 0.1 * i;
      for (int i = 1; i <= 14; i++)
        pIN_bins[155 + i] = 3. + 0.5 * i;

      const int kNbins_pINmore = 135;
      double pIN_binsmore[kNbins_pINmore + 1];
      for (int i = 0; i <= 120; i++)
        pIN_binsmore[i] = 0.01 * i + 0.3;
      for (int i = 1; i <= 10; i++)
        pIN_binsmore[120 + i] = 1.5 + 0.2 * i;
      pIN_binsmore[131] = 4.;
      pIN_binsmore[132] = 5.;
      pIN_binsmore[133] = 6.;
      pIN_binsmore[134] = 8.;
      pIN_binsmore[135] = 10.;

      const int kNbins_nSigma = 100;
      double nSigma_bins[kNbins_nSigma + 1];
      for (int i = 0; i <= kNbins_nSigma; i++)
        nSigma_bins[i] = -5. + 0.1 * i;

      const int kNbins_nSigmamore = 50;
      double nSigma_binsmore[kNbins_nSigmamore + 1];
      for (int i = 0; i <= kNbins_nSigmamore; i++)
        nSigma_binsmore[i] = -5. + 0.2 * i;

      const int kNbins_nSigmagrob = 24;
      double nSigma_binsgrob[kNbins_nSigmagrob + 1];
      for (int i = 0; i <= kNbins_nSigmagrob; i++)
        nSigma_binsgrob[i] = -6. + 0.5 * i;

      const int kNbins_TOFbeta = 120;
      double TOFbeta_bins[kNbins_TOFbeta + 1];
      for (int i = 0; i <= kNbins_TOFbeta; i++)
        TOFbeta_bins[i] = 0.01 * i;

      const int kNbins_TPCdEdx = 140;
      double TPCdEdx_bins[kNbins_TPCdEdx + 1];
      for (int i = 0; i <= kNbins_TPCdEdx; i++)
        TPCdEdx_bins[i] = i + 20;

      TArrayD nSigmaBinLimits[kNvarsPID];
      nSigmaBinLimits[0] = TArrayD(kNbins_pIN + 1, pIN_bins);
      nSigmaBinLimits[1] = TArrayD(kNbins_nSigma + 1, nSigma_bins);
      nSigmaBinLimits[2] = TArrayD(kNbins_nSigma + 1, nSigma_bins);

      TArrayD nSignalBinLimits[kNvarsPID];
      nSignalBinLimits[0] = TArrayD(kNbins_pIN + 1, pIN_bins);
      nSignalBinLimits[1] = TArrayD(kNbins_TPCdEdx + 1, TPCdEdx_bins);
      nSignalBinLimits[2] = TArrayD(kNbins_TOFbeta + 1, TOFbeta_bins);

      int varsPIDnSignal[kNvarsPID] = {VarManager::kPin, VarManager::kTPCsignal, VarManager::kTOFbeta};
      hm->AddHistogram(histClass, "nSignalTPCTOF", "", kNvarsPID, varsPIDnSignal, nSignalBinLimits);

      if (subGroupStr.Contains("more")) {
        const int kNvarsPIDmore = 4;
        TArrayD nSigmaBinLimitsmore[kNvarsPIDmore];
        nSigmaBinLimitsmore[0] = TArrayD(kNbins_pINmore + 1, pIN_binsmore);
        nSigmaBinLimitsmore[1] = TArrayD(kNbins_nSigmamore + 1, nSigma_binsmore);
        nSigmaBinLimitsmore[2] = TArrayD(kNbins_nSigmagrob + 1, nSigma_binsgrob);
        nSigmaBinLimitsmore[3] = TArrayD(kNbins_nSigmamore + 1, nSigma_binsmore);
        int varsPIDnSigmamore[kNvarsPIDmore] = {VarManager::kPin, VarManager::kTPCnSigmaEl, VarManager::kTPCnSigmaPi, VarManager::kTOFnSigmaEl};
        hm->AddHistogram(histClass, "nSigmaTPCTOF", "", kNvarsPIDmore, varsPIDnSigmamore, nSigmaBinLimitsmore);
      } else {
        int varsPIDnSigma[kNvarsPID] = {VarManager::kPin, VarManager::kTPCnSigmaEl, VarManager::kTOFnSigmaEl};
        hm->AddHistogram(histClass, "nSigmaTPCTOF", "", kNvarsPID, varsPIDnSigma, nSigmaBinLimits);
      }
    }
    if (subGroupStr.Contains("runbyrun")) {
      hm->AddHistogram(histClass, "TPCncls_run", "Number of cluster in TPC vs RunNumber", false, (VarManager::GetNRuns() > 0 ? VarManager::GetNRuns() : 1), -0.5, -0.5 + VarManager::GetNRuns(), VarManager::kRunIndex,
                       160, -0.5, 159.5, VarManager::kTPCncls, 10, 0., 1., VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "TPCdEdx_run", "TPCdEdx vs RunNumber", false, (VarManager::GetNRuns() > 0 ? VarManager::GetNRuns() : 1), -0.5, -0.5 + VarManager::GetNRuns(), VarManager::kRunIndex,
                       300, 0., 300., VarManager::kTPCsignal, 10, 0., 1., VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "TPCchi2_run", "TPCchi2 vs RunNumber", false, (VarManager::GetNRuns() > 0 ? VarManager::GetNRuns() : 1), -0.5, -0.5 + VarManager::GetNRuns(), VarManager::kRunIndex,
                       100, 0., 10., VarManager::kTPCchi2, 10, 0., 1., VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "Pt_Run", "p_{T} distribution", false, VarManager::GetDummyNRuns(), -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo,
                       2000, 0.0, 20.0, VarManager::kPt, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "ITSncls_Run", "Number of cluster in ITS", false, 100, -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo,
                       8, -0.5, 7.5, VarManager::kITSncls, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "ITSchi2_Run", "ITS chi2", false, VarManager::GetDummyNRuns(), -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo,
                       100, 0.0, 50.0, VarManager::kITSchi2, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "TPCncls_Run", "Number of cluster in TPC", false, 100, -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo,
                       160, -0.5, 159.5, VarManager::kTPCncls, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "TPCchi2_Run", "TPC chi2", false, VarManager::GetDummyNRuns(), -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo,
                       100, 0.0, 10.0, VarManager::kTPCchi2, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "TPCdedx_Run", "TPC dE/dx", false, VarManager::GetDummyNRuns(), -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo,
                       200, 0.0, 200., VarManager::kTPCsignal, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "TPCnSigEle_Run", "TPC n-#sigma(e)", false, 100, -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo,
                       100, -5.0, 5.0, VarManager::kTPCnSigmaEl, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "DCAxy_Run", "DCA_{xy}", false, VarManager::GetDummyNRuns(), -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo,
                       400, -2.0, 2.0, VarManager::kTrackDCAxy, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
      hm->AddHistogram(histClass, "DCAz_Run", "DCA_{z}", false, VarManager::GetDummyNRuns(), -0.5 + VarManager::GetDummyFirst(), 0.5 + VarManager::GetDummyLast(), VarManager::kRunNo,
                       800, -4.0, 4.0, VarManager::kTrackDCAz, 1, 0, 1, VarManager::kNothing, VarManager::GetRunStr().Data());
    }
    if (subGroupStr.Contains("dca")) {
      hm->AddHistogram(histClass, "DCAxy", "DCA_{xy}", false, 400, -2.0, 2.0, VarManager::kTrackDCAxy);
      hm->AddHistogram(histClass, "DCAz", "DCA_{z}", false, 800, -4.0, 4.0, VarManager::kTrackDCAz);
      hm->AddHistogram(histClass, "DCAsigXY", "DCA_{XY} [#sigma]", false, 100, -10.0, 10.0, VarManager::kTrackDCAsigXY);
      hm->AddHistogram(histClass, "DCAsigZ", "DCA_{Z} [#sigma]", false, 100, -10.0, 10.0, VarManager::kTrackDCAsigZ);
      if (subGroupStr.Contains("pt")) {
        hm->AddHistogram(histClass, "Pt_DCAxy", "p_{T} vs DCA_{xy}", false, 200, 0.0, 20.0, VarManager::kPt, 400, -2.0, 2.0, VarManager::kTrackDCAxy);
        hm->AddHistogram(histClass, "Pt_DCAz", "p_{T} vs DCA_{z}", false, 200, 0.0, 20.0, VarManager::kPt, 800, -4.0, 4.0, VarManager::kTrackDCAz);
        hm->AddHistogram(histClass, "Pt_DCAsigXY", "p_{T} vs DCA_{XY} [#sigma]", false, 200, 0.0, 20.0, VarManager::kPt, 100, -10.0, 10.0, VarManager::kTrackDCAsigXY); // JJ:edit
        hm->AddHistogram(histClass, "Pt_DCAsigZ", "p_{T} vs DCA_{Z} [#sigma]", false, 200, 0.0, 20.0, VarManager::kPt, 100, -10.0, 10.0, VarManager::kTrackDCAsigZ);
        hm->AddHistogram(histClass, "Pt_DCAresXY", "p_{T} vs #DeltaDCA_{XY}", false, 200, 0.0, 10.0, VarManager::kPt, 100, -0.03, 0.03, VarManager::kTrackDCAresXY);
        hm->AddHistogram(histClass, "Pt_DCAresZ", "p_{T} vs #DeltaDCA_{Z}", false, 200, 0.0, 10.0, VarManager::kPt, 100, -0.03, 0.03, VarManager::kTrackDCAresZ);
      }
      if (subGroupStr.Contains("eta")) {
        hm->AddHistogram(histClass, "Eta_DCAxy", "#eta vs DCA_{xy}", false, 20, -1.0, 1.0, VarManager::kEta, 400, -2.0, 2.0, VarManager::kTrackDCAxy);
        hm->AddHistogram(histClass, "Eta_DCAz", "#eta vs DCA_{z}", false, 20, -1.0, 1.0, VarManager::kEta, 800, -4.0, 4.0, VarManager::kTrackDCAz);
      }
      if (subGroupStr.Contains("dca_fine")) { // Fine binning
        hm->AddHistogram(histClass, "DCAxy_fine", "DCA_{xy}", false, 1000, -0.2, 0.2, VarManager::kTrackDCAxy);
        hm->AddHistogram(histClass, "DCAz_fine", "DCA_{z}", false, 1000, -0.2, 0.2, VarManager::kTrackDCAz);
        hm->AddHistogram(histClass, "IsSPDfirst_Pt_DCAxy_fine", "IsSPDfirst vs p_{T} vs DCA_{xy}", false, 200, 0.0, 20.0, VarManager::kPt, 1000, -0.2, 0.2, VarManager::kTrackDCAxy, 2, -0.5, 1.5, VarManager::kIsSPDfirst);
        hm->AddHistogram(histClass, "IsSPDfirst_Pt_DCAz_fine", "IsSPDfirst vs p_{T} vs DCA_{z}", false, 200, 0.0, 20.0, VarManager::kPt, 1000, -0.2, 0.2, VarManager::kTrackDCAz, 2, -0.5, 1.5, VarManager::kIsSPDfirst);
        if (subGroupStr.Contains("pt")) {
          hm->AddHistogram(histClass, "PtLow_DCAxy_fine", "p_{T} vs DCA_{xy}", false, 100, 0.0, 2.0, VarManager::kPt, 1000, -0.2, 0.2, VarManager::kTrackDCAxy);
          hm->AddHistogram(histClass, "PtLow_DCAz_fine", "p_{T} vs DCA_{z}", false, 100, 0.0, 2.0, VarManager::kPt, 1000, -0.2, 0.2, VarManager::kTrackDCAz);
          hm->AddHistogram(histClass, "PtHigh_DCAxy_fine", "p_{T} vs DCA_{xy}", false, 200, 0.0, 20.0, VarManager::kPt, 1000, -0.05, 0.05, VarManager::kTrackDCAxy);
          hm->AddHistogram(histClass, "PtHigh_DCAz_fine", "p_{T} vs DCA_{z}", false, 200, 0.0, 20.0, VarManager::kPt, 1000, -0.05, 0.05, VarManager::kTrackDCAz);
        }
      }
    }
    if (subGroupStr.Contains("muon")) {
      if (!subGroupStr.Contains("ambiguity")) {
        hm->AddHistogram(histClass, "MuonNClusters", "", false, 100, 0.0, 10.0, VarManager::kMuonNClusters);
        hm->AddHistogram(histClass, "pdca", "", false, 200, 0.0, 1000., VarManager::kMuonPDca);
        hm->AddHistogram(histClass, "RAtAbsorberEnd", "", false, 100, 0.0, 200., VarManager::kMuonRAtAbsorberEnd);
        hm->AddHistogram(histClass, "Chi2", "", false, 100, 0.0, 200.0, VarManager::kMuonChi2);
        hm->AddHistogram(histClass, "Chi2MCHMID", "", false, 100, 0.0, 200.0, VarManager::kMuonChi2MatchMCHMID);
        hm->AddHistogram(histClass, "Chi2MCHMFT", "", false, 100, 0.0, 200.0, VarManager::kMuonChi2MatchMCHMFT);
        hm->AddHistogram(histClass, "Chi2MatchScoreMCHMFT", "", false, 100, 0.0, 200.0, VarManager::kMuonMatchScoreMCHMFT);
        hm->AddHistogram(histClass, "MuonCXX", "", false, 100, -1.0, 1.0, VarManager::kMuonCXX);
        hm->AddHistogram(histClass, "MuonCYY", "", false, 100, -1.0, 1.0, VarManager::kMuonCYY);
        hm->AddHistogram(histClass, "MuonCPhiPhi", "", false, 100, -1.0, 1.0, VarManager::kMuonCPhiPhi);
        hm->AddHistogram(histClass, "MuonCTglTgl", "", false, 100, -1.0, 1.0, VarManager::kMuonCTglTgl);
        hm->AddHistogram(histClass, "MuonC1Pt21Pt2", "", false, 100, -1.0, 1.0, VarManager::kMuonC1Pt21Pt2);
        hm->AddHistogram(histClass, "MCHBitMap_vs_pt", "MCH vs pt", false, 1025, 0.0, 1025.0, VarManager::kMCHBitMap, 400, 0, 100, VarManager::kPt);
        hm->AddHistogram(histClass, "MuonTime", "", false, 100, -1.0, 1.0, VarManager::kMuonTime);
        hm->AddHistogram(histClass, "MuonTimeRes", "", false, 100, -1.0, 1.0, VarManager::kMuonTimeRes);
        hm->AddHistogram(histClass, "MuonDcaX_vs_phi", "", false, 2000, -20.0, 20.0, VarManager::kMuonDCAx, 200, -2. * TMath::Pi(), 2. * TMath::Pi(), VarManager::kPhi);
        hm->AddHistogram(histClass, "MuonDcaY_vs_phi", "", false, 2000, -20.0, 20.0, VarManager::kMuonDCAy, 200, -2. * TMath::Pi(), 2. * TMath::Pi(), VarManager::kPhi);
        hm->AddHistogram(histClass, "MuonDcaX_vs_eta", "", false, 2000, -20.0, 20.0, VarManager::kMuonDCAx, 500, -5.0, 5.0, VarManager::kEta);
        hm->AddHistogram(histClass, "MuonDcaY_vs_eta", "", false, 2000, -20.0, 20.0, VarManager::kMuonDCAy, 500, -5.0, 5.0, VarManager::kEta);
      } else {
        hm->AddHistogram(histClass, "AmbiguityInBunch", "", false, 10, 0.0, 10., VarManager::kMuonNAssocsInBunch);
        hm->AddHistogram(histClass, "AmbiguityOutOfBunch", "", false, 10, 0.0, 10., VarManager::kMuonNAssocsOutOfBunch);
      }
    }

    if (subGroupStr.Contains("muon-pdca")) {
      hm->AddHistogram(histClass, "p", "p", false, 200, 0.0, 20.0, VarManager::kP);
      hm->AddHistogram(histClass, "pdca_vs_p", "pDCA vs p", false, 2000, 0.0, 20.0, VarManager::kP, 200, 0.0, 1000., VarManager::kMuonPDca);
      hm->AddHistogram(histClass, "pdca_vs_pt", "pDCA vs pt", false, 2000, 0.0, 20.0, VarManager::kPt, 200, 0.0, 1000., VarManager::kMuonPDca);
      hm->AddHistogram(histClass, "pdca_vs_Rabs", "pDCA vs R_{abs}", false, 100, 0., 200., VarManager::kMuonRAtAbsorberEnd, 200, 0.0, 1000., VarManager::kMuonPDca);
      hm->AddHistogram(histClass, "pdca_vs_Rabs_vs_p", "pDCA vs R_{abs} vs p", false, 2000, 0.0, 20.0, VarManager::kP, 100, 0., 200., VarManager::kMuonRAtAbsorberEnd, 200, 0.0, 1000., VarManager::kMuonPDca);
      hm->AddHistogram(histClass, "pdca_vs_Rabs_vs_pt", "pDCA vs R_{abs} vs pt", false, 2000, 0.0, 20.0, VarManager::kPt, 100, 0., 200., VarManager::kMuonRAtAbsorberEnd, 200, 0.0, 1000., VarManager::kMuonPDca);
    }
    if (subGroupStr.Contains("mft-pid")) {
      hm->AddHistogram(histClass, "hMftTrackEtaVsPt", "", false, 100, -5.f, -2.f, VarManager::kEta, 100, 0.f, 20.f, VarManager::kPt);
      hm->AddHistogram(histClass, "hMftNClusters", "", false, 16, -0.5f, 15.5f, VarManager::kMftNClusters);
      hm->AddHistogram(histClass, "hMftClusterSize", "", false, 101, -0.5f, 100.5f, VarManager::kMftClusterSize);
      hm->AddHistogram(histClass, "hMftMeanClusterSize", "", false, 200, 0.f, 20.f, VarManager::kMftMeanClusterSize);
    }
    if (subGroupStr.Contains("mc")) {
      hm->AddHistogram(histClass, "Pt_vs_PtMC", "pT vs MC pT", false, 200, 0.0, 20.0, VarManager::kPt, 200, 0.0, 20.0, VarManager::kMCPt);
      hm->AddHistogram(histClass, "Eta_vs_EtaMC", "#eta vs MC #eta", false, 50, -1.0, 1.0, VarManager::kEta, 50, -1.0, 1.0, VarManager::kMCEta);
      hm->AddHistogram(histClass, "Phi_vs_PhiMC", "#varphi vs MC #varphi", false, 50, 0.0, 2. * TMath::Pi(), VarManager::kPhi, 50, 0.0, 2. * TMath::Pi(), VarManager::kMCPhi);
      hm->AddHistogram(histClass, "TrackPDGcode", "PDG code of track", false, 10001, -5000, 5000, VarManager::kMCPdgCode);
    }
    if (subGroupStr.Contains("mcmother")) {
      hm->AddHistogram(histClass, "MotherPDGcode", "PDG code of mother", false, 10001, -5000, 5000, VarManager::kMCMotherPdgCode);
    }
    if (subGroupStr.Contains("dmeson")) {
      hm->AddHistogram(histClass, "Mass", "", false, 500, 0.0, 5.0, VarManager::kMass);
      hm->AddHistogram(histClass, "Rapidity", "", false, 400, -4.0, 4.0, VarManager::kRap);
    }
    if (subGroupStr.Contains("tpcpidvstofpid")) {
      hm->AddHistogram(histClass, "tpcNSigmaKa_tofNSigmaKa", "", false, 200, -10., 10., VarManager::kTPCnSigmaKa, 200, -10., 10., VarManager::kTOFnSigmaKa);
      hm->AddHistogram(histClass, "tpcNSigmaPi_tofNSigmaPi", "", false, 200, -10., 10., VarManager::kTPCnSigmaPi, 200, -10., 10., VarManager::kTOFnSigmaPi);
    }
  }

  if (!groupStr.CompareTo("mctruth_triple")) {
    hm->AddHistogram(histClass, "Eta_Pt", "", false, 100, -2.0, 2.0, VarManager::kPairEta, 200, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "Eta_Pt_lepton1", "", false, 100, -2.0, 2.0, VarManager::kEta1, 200, 0.0, 20.0, VarManager::kPt1);
    hm->AddHistogram(histClass, "Eta_Pt_lepton2", "", false, 100, -2.0, 2.0, VarManager::kEta2, 200, 0.0, 20.0, VarManager::kPt2);
    hm->AddHistogram(histClass, "Eta_Pt_Photon", "", false, 100, -2.0, 2.0, VarManager::kEta, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Phi_Eta", "#phi vs #eta distribution", false, 200, -5.0, 5.0, VarManager::kPairEta, 200, -2. * TMath::Pi(), 2. * TMath::Pi(), VarManager::kPairPhi);
    hm->AddHistogram(histClass, "Mass_Dilepton", "", false, 4500, 0.0, 4.5, VarManager::kPairMassDau);
    hm->AddHistogram(histClass, "Mass_Photon", "", false, 500, 0.0, 0.1, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Mass_Dilepton_Mass_Photon", "", false, 500, 0.0, 5.0, VarManager::kPairMassDau, 500, 0.0, 5.0, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Pt_Dilepton", "", false, 2000, 0.0, 20.0, VarManager::kPairPtDau);
    hm->AddHistogram(histClass, "Pt_Photon", "", false, 500, 0.0, 5.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Mass_DileptonPhoton", "", false, 4500, 0.0, 4.5, VarManager::kPairMass);
    hm->AddHistogram(histClass, "Pt_DileptonPhoton", "", false, 2000, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "Mass_Pt_DileptonPhoton", "", false, 500, 0.0, 5.0, VarManager::kPairMass, 200, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "DeltaMass", "", false, 1500, 0.0, 1.5, VarManager::kDeltaMass);
    hm->AddHistogram(histClass, "DeltaMass_ptdileptonphoton", "", false, 1000, 0.0, 1.0, VarManager::kDeltaMass, 3000, 0.0, 30.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "DeltaMass_Jpsi", "", false, 1500, 3, 4.5, (VarManager::kDeltaMass_jpsi));
    hm->AddHistogram(histClass, "Rapidity", "", false, 400, -4.0, 4.0, VarManager::kRap);
  }
  if (!groupStr.CompareTo("mctruth_pair")) {
    hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kMass, 40, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Pt", "", false, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Pt_Dilepton", "", false, 200, 0.0, 20.0, VarManager::kPairPtDau);
    hm->AddHistogram(histClass, "Eta_Pt_lepton1", "", false, 100, -2.0, 2.0, VarManager::kEta1, 200, 0.0, 20.0, VarManager::kPt1);
    hm->AddHistogram(histClass, "Eta_Pt_lepton2", "", false, 100, -2.0, 2.0, VarManager::kEta2, 200, 0.0, 20.0, VarManager::kPt2);
    hm->AddHistogram(histClass, "Mass", "", false, 500, 0.0, 5.0, VarManager::kMass);
    hm->AddHistogram(histClass, "Eta_Pt", "", false, 40, -2.0, 2.0, VarManager::kEta, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Phi_Eta", "#phi vs #eta distribution", false, 200, -5.0, 5.0, VarManager::kEta, 200, -2. * TMath::Pi(), 2. * TMath::Pi(), VarManager::kPhi);
  }
  if (!groupStr.CompareTo("mctruth_track")) {
    hm->AddHistogram(histClass, "PtMC", "MC pT", false, 200, 0.0, 20.0, VarManager::kMCPt);
    hm->AddHistogram(histClass, "EtaMC", "MC #eta", false, 50, -5.0, 5.0, VarManager::kMCEta);
    hm->AddHistogram(histClass, "PhiMC", "MC #phi", false, 50, -6.3, 6.3, VarManager::kMCPhi);
    hm->AddHistogram(histClass, "MCY", "MC y", false, 50, -5.0, 5.0, VarManager::kMCY);
    hm->AddHistogram(histClass, "EtaMC_PtMC", "", false, 40, -2.0, 2.0, VarManager::kMCEta, 200, 0.0, 20.0, VarManager::kMCPt);
    hm->AddHistogram(histClass, "VzMC", "MC vz", false, 100, -15.0, 15.0, VarManager::kMCVz);
    hm->AddHistogram(histClass, "VzMC_VtxZMC", "MC vz vs MC vtxZ", false, 50, -15.0, 15.0, VarManager::kMCVz, 50, -15.0, 15.0, VarManager::kMCVtxZ);
    hm->AddHistogram(histClass, "Weight", "", false, 50, 0.0, 5.0, VarManager::kMCParticleWeight);
  }

  if (!groupStr.CompareTo("pair")) {
    if (subGroupStr.Contains("barrel")) {
      hm->AddHistogram(histClass, "Mass", "", false, 500, 0.0, 5.0, VarManager::kMass);
      hm->AddHistogram(histClass, "Mass_HighRange", "", false, 375, 0.0, 15.0, VarManager::kMass);
      hm->AddHistogram(histClass, "Pt", "", false, 2000, 0.0, 20., VarManager::kPt);
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 125, 0.0, 5.0, VarManager::kMass, 40, 0.0, 20.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Eta_Pt", "", false, 40, -2.0, 2.0, VarManager::kEta, 40, 0.0, 20.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Mass_VtxZ", "", true, 30, -15.0, 15.0, VarManager::kVtxZ, 500, 0.0, 5.0, VarManager::kMass);
      if (subGroupStr.Contains("pbpb")) {
        hm->AddHistogram(histClass, "Mass_CentFT0C", "", false, 500, 0.0, 5.0, VarManager::kMass, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "Pt_CentFT0C", "", false, 500, 0.0, 1.5, VarManager::kPt, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "Mass_Pt_CentFT0C", "", false, 500, 0.0, 5.0, VarManager::kMass, 400, 0.0, 40.0, VarManager::kPt, 20, 0.0, 100.0, VarManager::kCentFT0C);
      }
      if (subGroupStr.Contains("mult")) {
        hm->AddHistogram(histClass, "Mass_Pt_MultFV0A", "", false, 200, 0.0, 5.0, VarManager::kMass, 40, 0.0, 40.0, VarManager::kPt, 100, 0.0, 25000.0, VarManager::kMultFV0A);
      }
      if (subGroupStr.Contains("polarization")) {
        hm->AddHistogram(histClass, "cosThetaHE", "", false, 100, -1., 1., VarManager::kCosThetaHE);
        hm->AddHistogram(histClass, "cosThetaCS", "", false, 100, -1., 1., VarManager::kCosThetaCS);
        hm->AddHistogram(histClass, "PhiHE", "", false, 100, -TMath::Pi(), TMath::Pi(), VarManager::kPhiHE);
        hm->AddHistogram(histClass, "PhiCS", "", false, 100, -TMath::Pi(), TMath::Pi(), VarManager::kPhiCS);
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaHE", "", false, 100, 1.0, 5.0, VarManager::kMass, 250, 0.0, 25.0, VarManager::kPt, 40, -1., 1., VarManager::kCosThetaHE);
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaCS", "", false, 100, 1.0, 5.0, VarManager::kMass, 250, 0.0, 25.0, VarManager::kPt, 40, -1., 1., VarManager::kCosThetaCS);
        hm->AddHistogram(histClass, "Mass_Pt_PhiHE", "", false, 100, 1.0, 5.0, VarManager::kMass, 250, 0.0, 25.0, VarManager::kPt, 40, -TMath::Pi(), TMath::Pi(), VarManager::kPhiHE);
        hm->AddHistogram(histClass, "Mass_Pt_PhiCS", "", false, 100, 1.0, 5.0, VarManager::kMass, 250, 0.0, 25.0, VarManager::kPt, 40, -TMath::Pi(), TMath::Pi(), VarManager::kPhiCS);
      }
      if (subGroupStr.Contains("upsilon")) {
        hm->AddHistogram(histClass, "MassUpsilon_Pt", "", false, 500, 7.0, 12.0, VarManager::kMass, 400, 0.0, 40.0, VarManager::kPt);
      }
      if (subGroupStr.Contains("dalitz")) {
        hm->AddHistogram(histClass, "MassLow", "", false, 500, 0.0, 0.05, VarManager::kMass);
        hm->AddHistogram(histClass, "PsiPair", "", false, 200, -1.5, 1.5, VarManager::kPsiPair);
        hm->AddHistogram(histClass, "PsiPair_DeltaPhi", "", false, 100, -0.5, 0.5, VarManager::kDeltaPhiPair, 100, -1.5, 1.5, VarManager::kPsiPair);
      }
      if (subGroupStr.Contains("vertexing")) {
        hm->AddHistogram(histClass, "UsedKF", "", false, 2, -0.5, 1.5, VarManager::kUsedKF);
        hm->AddHistogram(histClass, "Lz", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLz);
        hm->AddHistogram(histClass, "Lxy", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLxy);
        hm->AddHistogram(histClass, "Lxyz", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLxyz);
        hm->AddHistogram(histClass, "Tauz", "", false, 1000, -0.02, 0.02, VarManager::kVertexingTauz);
        hm->AddHistogram(histClass, "Tauxy", "", false, 1000, -0.03, 0.03, VarManager::kVertexingTauxy);
        hm->AddHistogram(histClass, "Tauxy_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -0.03, 0.03, VarManager::kVertexingTauxy);
        hm->AddHistogram(histClass, "Tauz_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -0.03, 0.03, VarManager::kVertexingTauz);
        hm->AddHistogram(histClass, "LzProj", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLzProjected);
        hm->AddHistogram(histClass, "LxyProj", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLxyProjected);
        hm->AddHistogram(histClass, "LxyzProj", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLxyzProjected);
        hm->AddHistogram(histClass, "TauzProj", "", false, 1000, -0.03, 0.03, VarManager::kVertexingTauzProjected);
        hm->AddHistogram(histClass, "TauxyProj", "", false, 1000, -0.03, 0.03, VarManager::kVertexingTauxyProjected);
        hm->AddHistogram(histClass, "TauxyzProj", "", false, 1000, -0.03, 0.03, VarManager::kVertexingTauxyzProjected);
        hm->AddHistogram(histClass, "LxyProj_Pt", "", false, 10, 0.0, 20.0, VarManager::kPt, 1000, -1.0, 1.0, VarManager::kVertexingLxyProjected);
        hm->AddHistogram(histClass, "LxyProj_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -1.0, 1.0, VarManager::kVertexingLxyProjected);
        hm->AddHistogram(histClass, "LzProj_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -1.0, 1.0, VarManager::kVertexingLzProjected);
      }

      if (subGroupStr.Contains("kalman-filter")) {
        hm->AddHistogram(histClass, "LxyErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyErr);
        hm->AddHistogram(histClass, "LxyzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyzErr);
        hm->AddHistogram(histClass, "TauxyErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingTauxyErr);
        hm->AddHistogram(histClass, "VtxingProcCode", "", false, 10, 0.0, 10.0, VarManager::kVertexingProcCode);
        hm->AddHistogram(histClass, "VtxingChi2PCA", "", false, 100, 0.0, 10.0, VarManager::kVertexingChi2PCA);
        hm->AddHistogram(histClass, "KFMass", "", false, 500, 0.0, 5.0, VarManager::kKFMass);
        hm->AddHistogram(histClass, "LxyOverDLxy", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyOverErr);
        hm->AddHistogram(histClass, "LzOverDLz", "", false, 100, 0.0, 10.0, VarManager::kVertexingLzOverErr);
        hm->AddHistogram(histClass, "LxyzOverDLxyz", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyzOverErr);
        hm->AddHistogram(histClass, "KFTauxy", "", false, 1000, -0.2, 0.2, VarManager::kVertexingTauxy);
        hm->AddHistogram(histClass, "KFTrack0DCAxyz", "", false, 400, -2.0, 2.0, VarManager::kKFTrack0DCAxyz);
        hm->AddHistogram(histClass, "KFTrack1DCAxyz", "", false, 400, -2.0, 2.0, VarManager::kKFTrack1DCAxyz);
        hm->AddHistogram(histClass, "KFTracksDCAxyzMax", "", false, 400, -2.0, 2.0, VarManager::kKFTracksDCAxyzMax);
        hm->AddHistogram(histClass, "KFDCAxyzBetweenProngs", "", false, 400, -2.0, 2.0, VarManager::kKFDCAxyzBetweenProngs);
        hm->AddHistogram(histClass, "KFTrack0DCAxy", "", false, 400, -2.0, 2.0, VarManager::kKFTrack0DCAxy);
        hm->AddHistogram(histClass, "KFTrack1DCAxy", "", false, 400, -2.0, 2.0, VarManager::kKFTrack1DCAxy);
        hm->AddHistogram(histClass, "KFTracksDCAxyMax", "", false, 400, -2.0, 2.0, VarManager::kKFTracksDCAxyMax);
        hm->AddHistogram(histClass, "KFDCAxyBetweenProngs", "", false, 400, -2.0, 2.0, VarManager::kKFDCAxyBetweenProngs);
        hm->AddHistogram(histClass, "KFChi2OverNDFGeo", "", false, 150, -5, 10, VarManager::kKFChi2OverNDFGeo);
        hm->AddHistogram(histClass, "KFTrack0DeviationFromPV", "", false, 150, 0, 15e+6, VarManager::kKFTrack0DeviationFromPV);
        hm->AddHistogram(histClass, "KFTrack1DeviationFromPV", "", false, 150, 0, 15e+6, VarManager::kKFTrack1DeviationFromPV);
        hm->AddHistogram(histClass, "KFTrack0DeviationxyFromPV", "", false, 150, 0, 15e+6, VarManager::kKFTrack0DeviationxyFromPV);
        hm->AddHistogram(histClass, "KFTrack1DeviationxyFromPV", "", false, 150, 0, 15e+6, VarManager::kKFTrack1DeviationxyFromPV);
        hm->AddHistogram(histClass, "KFPairDCAxyz", "", false, 400, -2.0, 2.0, VarManager::kKFJpsiDCAxyz);
        hm->AddHistogram(histClass, "KFPairDCAxy", "", false, 400, -2.0, 2.0, VarManager::kKFJpsiDCAxy);
        hm->AddHistogram(histClass, "KFPairDeviationFromPV", "", false, 150, 0, 15e+6, VarManager::kKFPairDeviationFromPV);
        hm->AddHistogram(histClass, "KFPairDeviationxyFromPV", "", false, 150, 0, 15e+6, VarManager::kKFPairDeviationxyFromPV);
        hm->AddHistogram(histClass, "KFCosPA", "", false, 300, -1.5, 1.5, VarManager::kKFCosPA);
        hm->AddHistogram(histClass, "KFMassGeoTop", "", false, 500, 0.0, 5.0, VarManager::kKFMassGeoTop);
        hm->AddHistogram(histClass, "KFChi2OverNDFGeoTop", "", false, 150, -5, 10, VarManager::kKFChi2OverNDFGeoTop);
        hm->AddHistogram(histClass, "KFNTrks2PV", "", false, 210, -10, 200, VarManager::kKFNContributorsPV);
        hm->AddHistogram(histClass, "Mass_DCAxyzTwoProngs", "", false, 500, 0.0, 5.0, VarManager::kMass, 400, -2.0, 2.0, VarManager::kKFDCAxyzBetweenProngs);
        hm->AddHistogram(histClass, "Mass_DCAxyTwoProngs", "", false, 500, 0.0, 5.0, VarManager::kMass, 400, -2.0, 2.0, VarManager::kKFDCAxyBetweenProngs);
        hm->AddHistogram(histClass, "Mass_KFChi2OverNDFGeo", "", false, 500, 0.0, 5.0, VarManager::kMass, 150, -5, 10, VarManager::kKFChi2OverNDFGeo);
      }
      if (subGroupStr.Contains("run2-vertexing-definitions")) {
        hm->AddHistogram(histClass, "LzProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLzProjected);
        hm->AddHistogram(histClass, "LxyProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyProjected);
        hm->AddHistogram(histClass, "LxyzProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyzProjected);
        hm->AddHistogram(histClass, "TauzProj", "", false, 4000, -0.5, 0.5, VarManager::kVertexingTauzProjected);
        hm->AddHistogram(histClass, "TauxyProjNs", "", false, 4000, -0.1, 0.1, VarManager::kVertexingTauxyProjectedNs);
        hm->AddHistogram(histClass, "TauxyProjCm", "", false, 4000, -0.5, 0.5, VarManager::kVertexingTauxyProjected);
      }
      if (subGroupStr.Contains("multidimentional-vertexing-histograms")) {
        hm->AddHistogram(histClass, "pT_TauxyProj", "", false, 1000, -0.2, 0.2, VarManager::kVertexingTauxyProjected, 20, 0.0, 20., VarManager::kPt);
        hm->AddHistogram(histClass, "InvMass_TauxyProj", "", false, 500, 0.0, 5.0, VarManager::kMass, 1000, -0.2, 0.2, VarManager::kVertexingTauxyProjected);
        hm->AddHistogram(histClass, "Eta_TauxyProj", "", false, 40, -2.0, 2.0, VarManager::kEta, 1000, -0.2, 0.2, VarManager::kVertexingTauxyProjected);
        hm->AddHistogram(histClass, "Rap_TauxyProj", "", false, 200, -1.0, 1.0, VarManager::kRap, 1000, -0.2, 0.2, VarManager::kVertexingTauxyProjected);

        const int kNvarsPair = 4;
        const int kInvMassNbins = 3;
        double InvMassBinLims[kInvMassNbins + 1] = {2.2, 2.6, 3.4, 3.6};

        const int kPtNbins = 10;
        double PtBinLims[kPtNbins + 1] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 12., 20.};

        const int kTauNBins = 500;
        double TauBinLims[kTauNBins + 1];
        for (int i = 0; i <= kTauNBins; ++i)
          TauBinLims[i] = -0.3 + (0.0015 * i);

        TArrayD nCutsBinLimits[kNvarsPair];
        nCutsBinLimits[0] = TArrayD(kInvMassNbins + 1, InvMassBinLims);
        nCutsBinLimits[1] = TArrayD(kPtNbins + 1, PtBinLims);
        nCutsBinLimits[2] = TArrayD(kTauNBins + 1, TauBinLims);
        nCutsBinLimits[3] = TArrayD(kTauNBins + 1, TauBinLims);

        int varsPair[kNvarsPair] = {VarManager::kMass, VarManager::kPt, VarManager::kVertexingTauzProjected, VarManager::kVertexingTauxyProjected};
        hm->AddHistogram(histClass, "tau_MultiD", "Invariant mass vs. pT vs. eta vs. rapidity vs. Run2 tau", kNvarsPair, varsPair, nCutsBinLimits);
      }

      if (subGroupStr.Contains("flow")) {
        hm->AddHistogram(histClass, "Mass_u2q2", "u_{2}Q_{2}^{A} vs m", true, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kU2Q2);
        hm->AddHistogram(histClass, "Mass_u3q3", "u_{3}Q_{3}^{A} vs m", true, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kU3Q3);
        hm->AddHistogram(histClass, "Mass_cos2DeltaPhi", "cos 2(#varphi-#Psi_{2}^{A}) vs m", true, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kCos2DeltaPhi);
        hm->AddHistogram(histClass, "Mass_cos3DeltaPhi", "cos 3(#varphi-#Psi_{3}^{A}) vs m", true, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kCos3DeltaPhi);
      }
    } else if (subGroupStr.Contains("dimuon")) {
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 750, 0.0, 15.0, VarManager::kMass, 120, 0.0, 30.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Mass_Rapidity", "", false, 750, 0.0, 15.0, VarManager::kMass, 150, 2.5, 4.0, VarManager::kRap);
      if (subGroupStr.Contains("dimuon-multi-diff")) {
        int varsKine[3] = {VarManager::kMass, VarManager::kPt, VarManager::kRap};
        int binsKine[3] = {250, 120, 60};
        double xminKine[3] = {0.0, 0.0, 2.5};
        double xmaxKine[3] = {5.0, 30.0, 4.0};
        hm->AddHistogram(histClass, "Mass_Pt_Rapidity", "", 3, varsKine, binsKine, xminKine, xmaxKine, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-centr")) {
        hm->AddHistogram(histClass, "Mass_CentFT0C", "", false, 750, 0.0, 15.0, VarManager::kMass, 100, 0., 100., VarManager::kCentFT0C);
      }
      if (subGroupStr.Contains("qc")) {
        hm->AddHistogram(histClass, "Mass_VtxZ", "", true, 30, -15.0, 15.0, VarManager::kVtxZ, 750, 0.0, 15.0, VarManager::kMass);
        hm->AddHistogram(histClass, "DeltaPtotTracks", "", false, 2000, -100., 100., VarManager::kDeltaPtotTracks);
        hm->AddHistogram(histClass, "Mass_DeltaPtotTracks", "", false, 150, 2.0, 5.0, VarManager::kMass, 200, -100., 100., VarManager::kDeltaPtotTracks);
      }
      if (subGroupStr.Contains("mixedevent")) {
        hm->AddHistogram(histClass, "Mass_cos2DeltaPhiMu1", "cos 2(#varphi_{#mu1}-#phi_{#mu#mu}) vs m", false, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kCos2DeltaPhiMu1);
        hm->AddHistogram(histClass, "Mass_cos2DeltaPhiMu2", "cos 2(#varphi_{#mu2}-#phi_{#mu#mu}) vs m", false, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kCos2DeltaPhiMu2);
        hm->AddHistogram(histClass, "R2SP1_CentFT0C", "mass vs centrality vs. R2SP_event1", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kTwoR2SP1);
        hm->AddHistogram(histClass, "R2SP2_CentFT0C", "mass vs centrality vs. R2SP_event2", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kTwoR2SP2);
        hm->AddHistogram(histClass, "R2EP1_CentFT0C", "mass vs centrality vs. R2EP_event1", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kTwoR2EP1);
        hm->AddHistogram(histClass, "R2EP2_CentFT0C", "mass vs centrality vs. R2EP_event2", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kTwoR2EP2);
        hm->AddHistogram(histClass, "U2Q2_CentFT0C_ev1", "mass vs. centrality vs. U2Q2_event1", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kU2Q2Ev1);
        hm->AddHistogram(histClass, "U2Q2_CentFT0C_ev2", "mass vs. centrality vs. U2Q2_event2", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kU2Q2Ev2);
      }
      if (subGroupStr.Contains("dimuon-polarization-he")) {
        int varspTHE[4] = {VarManager::kMass, VarManager::kPt, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int varsFT0CMulHE[4] = {VarManager::kMass, VarManager::kMultFT0C, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int binspT[4] = {100, 20, 20, 20};
        int binsMul[4] = {100, 20, 20, 20};
        double xminpT[4] = {1., 0., -1., -3.14};
        double xmaxpT[4] = {5., 20., 1., +3.14};
        double xminMul[4] = {1., 0., -1., -3.14};
        double xmaxMul[4] = {5., 2000., 1., +3.14};
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaHE_phiHE", "", 4, varspTHE, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_MultFT0C_cosThetaHE_phiHE", "", 4, varsFT0CMulHE, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-polarization-cs")) {
        int varspTCS[4] = {VarManager::kMass, VarManager::kPt, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int varsFT0CMulCS[4] = {VarManager::kMass, VarManager::kMultFT0C, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int binspT[4] = {100, 20, 20, 20};
        int binsMul[4] = {100, 20, 20, 20};
        double xminpT[4] = {1., 0., -1., -3.14};
        double xmaxpT[4] = {5., 20., 1., +3.14};
        double xminMul[4] = {1., 0., -1., -3.14};
        double xmaxMul[4] = {5., 2000., 1., +3.14};
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaCS_phiCS", "", 4, varspTCS, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_MultFT0C_cosThetaCS_phiCS", "", 4, varsFT0CMulCS, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("multiplicity-fvoa")) {
        int varsFV0AMulHE[4] = {VarManager::kMass, VarManager::kMultFV0A, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int varsFV0AMulCS[4] = {VarManager::kMass, VarManager::kMultFV0A, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int binsMul[4] = {100, 20, 20, 20};
        double xminMul[4] = {1., 0., -1., -3.14};
        double xmaxMul[4] = {5., 5000., 1., +3.14};
        hm->AddHistogram(histClass, "Mass_MultFV0A_cosThetaHE_phiHE", "", 4, varsFV0AMulHE, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_MultFV0A_cosThetaCS_phiCS", "", 4, varsFV0AMulCS, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("multiplicity-tpc")) {
        int varsTPCMulHE[4] = {VarManager::kMass, VarManager::kMultTPC, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int varsTPCMulCS[4] = {VarManager::kMass, VarManager::kMultTPC, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int binsMul[4] = {100, 20, 20, 20};
        double xminMul[4] = {1., 0., -1., -3.14};
        double xmaxMul[4] = {5., 2000., 1., +3.14};
        hm->AddHistogram(histClass, "Mass_MultTPC_cosThetaCS_phiCS", "", 4, varsTPCMulCS, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_MultTPC_cosThetaHE_phiHE", "", 4, varsTPCMulHE, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("vertexing-forward")) {
        hm->AddHistogram(histClass, "Lxyz", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyz);
        hm->AddHistogram(histClass, "Lz", "", false, 100, 0.0, 10.0, VarManager::kVertexingLz);
        hm->AddHistogram(histClass, "Tauz", "", false, 100, -0.01, 0.01, VarManager::kVertexingTauz);
        hm->AddHistogram(histClass, "LxyzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyzErr);
        hm->AddHistogram(histClass, "LzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLzErr);
        hm->AddHistogram(histClass, "TauzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingTauzErr);
        hm->AddHistogram(histClass, "VtxingProcCode", "", false, 10, 0.0, 10.0, VarManager::kVertexingProcCode);
        hm->AddHistogram(histClass, "VtxingChi2PCA", "", false, 100, 0.0, 10.0, VarManager::kVertexingChi2PCA);
        hm->AddHistogram(histClass, "Pz", "", false, 200, 0.0, 100.0, VarManager::kVertexingPz);
        hm->AddHistogram(histClass, "Secondary_Vertexing_z", "", false, 200, -10.0, 10.0, VarManager::kVertexingSV);
        hm->AddHistogram(histClass, "Primary_Vertexing_z", "", false, 200, -10.0, 10.0, VarManager::kVtxZ);
      }
      if (subGroupStr.Contains("pbpb")) {
        hm->AddHistogram(histClass, "Mass_CentVZERO", "", false, 750, 0.0, 15.0, VarManager::kMass, 100, 0., 100., VarManager::kCentVZERO);
        hm->AddHistogram(histClass, "Pt_CentVZERO", "", false, 120, 0.0, 30.0, VarManager::kPt, 100, 0., 100., VarManager::kCentVZERO);
        hm->AddHistogram(histClass, "Rapidity_CentVZERO", "", false, 200, 2.5, 4.0, VarManager::kRap, 100, 0., 100., VarManager::kCentVZERO);
        hm->AddHistogram(histClass, "Mass_CentFT0C", "", false, 750, 0.0, 15.0, VarManager::kMass, 100, 0., 100., VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "Pt_CentFT0C", "", false, 120, 0.0, 30.0, VarManager::kPt, 100, 0., 100., VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "Rapidity_CentFT0C", "", false, 200, 2.5, 4.0, VarManager::kRap, 100, 0., 100., VarManager::kCentFT0C);
      }
      if (subGroupStr.Contains("lowmass")) {
        hm->AddHistogram(histClass, "MassLow", "", false, 400, 0.0, 2.0, VarManager::kMass);
      }
      if (subGroupStr.Contains("lmmumu")) {
        hm->AddHistogram(histClass, "Mass_QuadDCAabsXY", "", false, 250, 0.0, 5.0, VarManager::kMass, 500, 0.0, 1, VarManager::kQuadDCAabsXY);
      }
      if (subGroupStr.Contains("flow-dimuon")) {
        int varV2[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kU2Q2, VarManager::kCos2DeltaPhi};
        int varV3[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kU3Q3, VarManager::kCos3DeltaPhi};

        int bins[5] = {250, 60, 9, 200, 200};
        double minBins[5] = {0.0, 0.0, 0.0, -10.0, -10.0};
        double maxBins[5] = {5.0, 30.0, 90.0, 10.0, 10.0};
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_V2", "", 5, varV2, bins, minBins, maxBins, 0, -1, kTRUE);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_V3", "", 5, varV3, bins, minBins, maxBins, 0, -1, kTRUE);
      }
      if (subGroupStr.Contains("res-flow-dimuon")) {
        int varV2[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kR2SP_AB, VarManager::kR2EP_AB};
        int varV3[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kR3SP, VarManager::kR3EP};

        int bins[5] = {125, 120, 9, 200, 200};
        double minBins[5] = {0.0, 0.0, 0.0, -10.0, -10.0};
        double maxBins[5] = {5.0, 30.0, 90.0, 10.0, 10.0};
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_R2", "", 5, varV2, bins, minBins, maxBins, 0, -1, kTRUE);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_R3", "", 5, varV3, bins, minBins, maxBins, 0, -1, kTRUE);
      }
      if (subGroupStr.Contains("z-boson")) {
        hm->AddHistogram(histClass, "MassZboson", "", false, 240, 20.0, 140.0, VarManager::kMass);
      }
      if (subGroupStr.Contains("dimuon-vtxncontrib")) {
        hm->AddHistogram(histClass, "MassMult", "", false, 750, 0.0, 15.0, VarManager::kMass, 301, -0.5, 300.5, VarManager::kVtxNcontrib);
        hm->AddHistogram(histClass, "MassVtxZMult", "", false, 300, 0.0, 6.0, VarManager::kMass, 60, -15.0, 15.0, VarManager::kVtxZ, 100, 0.0, 100.0, VarManager::kVtxNcontrib);
      }
    } else if (subGroupStr.Contains("electronmuon")) {
      hm->AddHistogram(histClass, "Mass", "", false, 750, 0.0, 30.0, VarManager::kMass);
      hm->AddHistogram(histClass, "Pt", "", false, 120, 0.0, 30.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Rapidity", "", false, 500, -1.0, 4.0, VarManager::kRap);
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 750, 0.0, 30.0, VarManager::kMass, 120, 0.0, 30.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Mass_Rapidity", "", false, 750, 0.0, 30.0, VarManager::kMass, 500, -1.0, 4.0, VarManager::kRap);
      hm->AddHistogram(histClass, "Mass_VtxZ", "", true, 30, -15.0, 15.0, VarManager::kVtxZ, 750, 0.0, 30.0, VarManager::kMass);
      hm->AddHistogram(histClass, "DeltaPhiPair", "", false, 130, -6.5, 6.5, VarManager::kDeltaPhiPair);
    }
    if (subGroupStr.Contains("dielectrons")) {
      if (subGroupStr.Contains("prefilter")) {
        hm->AddHistogram(histClass, "MassLow_OpeningAngle", "", false, 150, 0., 0.15, VarManager::kMass, 80, 0., 0.8, VarManager::kOpeningAngle);
      }
      if (subGroupStr.Contains("phiv")) {
        hm->AddHistogram(histClass, "Mass_Pt_PhiV", "", false, 20, 0.0, 0.2, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt, 100, 0.0, TMath::Pi(), VarManager::kPairPhiv);
      }
      if (subGroupStr.Contains("double-phi-v")) {
        hm->AddHistogram(histClass, "Mass_Pt_PhiV", "", false, 20, 0.0, 0.2, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt, 100, 0.0, TMath::Pi(), VarManager::kPairPhiv, "", "", "", -1, -1, true);
      }
      if (subGroupStr.Contains("largemass-phi-v")) {
        // binning for mee at large scales:
        // every 10 MeV from 0 to 0.2 GeV/c2
        // every 100 MeV from 0.2 to 1. GeV/c2
        // every 500 GeV from 1 to 5 GeV/c2
        double mee_bins[37];
        for (int i = 0; i <= 20; i++)
          mee_bins[i] = 0.01 * i;
        for (int i = 1; i <= 8; i++)
          mee_bins[20 + i] = 0.2 + 0.1 * i;
        for (int i = 1; i <= 8; i++)
          mee_bins[28 + i] = 1. + 0.5 * i;
        int nbins_mee = sizeof(mee_bins) / sizeof(*mee_bins) - 1;

        // binning for ptee at large scales:
        // every 0.2 GeV/c from 0 to 10 GeV/c
        double ptee_bins[51];
        for (int i = 0; i <= 50; i++)
          ptee_bins[i] = 0.2 * i;
        int nbins_ptee = sizeof(ptee_bins) / sizeof(*ptee_bins) - 1;

        // binning for phiv:
        // steps of size pi/100
        double phiv_bins[101];
        for (int i = 0; i <= 100; i++)
          phiv_bins[i] = TMath::Pi() / 100. * i;
        int nbins_phiv = sizeof(phiv_bins) / sizeof(*phiv_bins) - 1;

        // 3D histo
        hm->AddHistogram(histClass, "Mass_Pt_PhiV", "", false, nbins_mee, mee_bins, VarManager::kMass, nbins_ptee, ptee_bins, VarManager::kPt, nbins_phiv, phiv_bins, VarManager::kPairPhiv);
      }
      if (subGroupStr.Contains("meeptee")) {
        hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt);
      }
      if (subGroupStr.Contains("double-mee-ptee")) {
        hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt, 0, 0, 0, -1, "", "", "", -1, -1, true);
      }
      if (subGroupStr.Contains("lmee")) {
        hm->AddHistogram(histClass, "Mass_QuadDCAsigXY", "", false, 50, 0.0, 5.0, VarManager::kMass, 50, 0.0, 20.0, VarManager::kQuadDCAsigXY);
        hm->AddHistogram(histClass, "Mass_QuadDCAsigZ", "", false, 50, 0.0, 5.0, VarManager::kMass, 50, 0.0, 20.0, VarManager::kQuadDCAsigZ);
        hm->AddHistogram(histClass, "Mass_Pt_QuadDCAsigXYZ", "", false, 500, 0.0, 5.0, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt, 50, 0.0, 20.0, VarManager::kQuadDCAsigXYZ);
        hm->AddHistogram(histClass, "Mass_Pt_SignQuadDCAsigXY", "", false, 50, 0.0, 5.0, VarManager::kMass, 50, 0.0, 10.0, VarManager::kPt, 100, -20.0, 20.0, VarManager::kSignQuadDCAsigXY);
      }
      if (subGroupStr.Contains("largescale")) {
        // binning for mee at large scales:
        // every 10 MeV from 0 to 1.1 GeV/c2
        // every 50 MeV from 1.1 to 2.7 GeV/c2
        // every 10 MeV from 2.7 to 3.2 GeV/c2
        // every 50 MeV from 3.2 to 12 GeV/c2
        double mee_bins[369];
        for (int i = 0; i <= 110; i++)
          mee_bins[i] = 0.01 * i;
        for (int i = 1; i <= 32; i++)
          mee_bins[110 + i] = 1.1 + 0.05 * i;
        for (int i = 1; i <= 50; i++)
          mee_bins[142 + i] = 2.7 + 0.01 * i;
        for (int i = 1; i <= 176; i++)
          mee_bins[192 + i] = 3.2 + 0.05 * i;
        int nbins_mee = sizeof(mee_bins) / sizeof(*mee_bins) - 1;

        // binning for ptee at large scales:
        // every 0.1 GeV/c from 0 to 10 GeV/c
        // every 0.5 GeV/c from 10 to 30 GeV/c
        double ptee_bins[201];
        for (int i = 0; i <= 100; i++)
          ptee_bins[i] = 0.1 * i;
        for (int i = 1; i <= 100; i++)
          ptee_bins[100 + i] = 10 + 0.2 * i;
        int nbins_ptee = sizeof(ptee_bins) / sizeof(*ptee_bins) - 1;

        // binning for dca at large scales:
        // every 0.1 sigma from 0 to 5 sigma
        // every 0.5 sigma from 5 to 10 sigma
        // every 1.0 sigma from 10 to 40 sigma
        double dca_bins[91];
        for (int i = 0; i <= 50; i++)
          dca_bins[i] = 0.1 * i;
        for (int i = 1; i <= 10; i++)
          dca_bins[50 + i] = 5 + 0.5 * i;
        for (int i = 1; i <= 30; i++)
          dca_bins[60 + i] = 10 + 1 * i;
        int nbins_dca = sizeof(dca_bins) / sizeof(*dca_bins) - 1;

        hm->AddHistogram(histClass, "Mass_QuadDCAsigXY", "", false, nbins_mee, mee_bins, VarManager::kMass, nbins_dca, dca_bins, VarManager::kQuadDCAsigXY);
        hm->AddHistogram(histClass, "Mass_QuadDCAsigZ", "", false, nbins_mee, mee_bins, VarManager::kMass, nbins_dca, dca_bins, VarManager::kQuadDCAsigZ);
        hm->AddHistogram(histClass, "Mass_Pt_QuadDCAsigXYZ", "", false, nbins_mee, mee_bins, VarManager::kMass, nbins_ptee, ptee_bins, VarManager::kPt, nbins_dca, dca_bins, VarManager::kQuadDCAsigXYZ);
      }
    }
    if (subGroupStr.Contains("opencharm")) {
      if (subGroupStr.Contains("dmeson")) {
        double mD0_bins[51];
        for (int i = 0; i <= 50; i++) {
          mD0_bins[i] = 1.7 + i * 0.006;
        }
        int nbins_mD0 = sizeof(mD0_bins) / sizeof(*mD0_bins) - 1;

        double ptD0_bins[31];
        for (int i = 0; i <= 16; i++)
          ptD0_bins[i] = 0.125 * i;
        for (int i = 1; i <= 8; i++)
          ptD0_bins[16 + i] = 2 + 0.25 * i;
        for (int i = 1; i <= 4; i++)
          ptD0_bins[24 + i] = 4 + 1 * i;
        ptD0_bins[29] = 8;
        ptD0_bins[30] = 20;
        int nbins_ptD0 = sizeof(ptD0_bins) / sizeof(*ptD0_bins) - 1;
        hm->AddHistogram(histClass, "MassD0region", "", false, nbins_mD0, mD0_bins, VarManager::kMass);
        hm->AddHistogram(histClass, "MassD0region_Pt", "", false, nbins_mD0, mD0_bins, VarManager::kMass, nbins_ptD0, ptD0_bins, VarManager::kPt);
        hm->AddHistogram(histClass, "MassD0region_eta", "", false, 50, 1.7, 2.0, VarManager::kMass, 40, -2., 2., VarManager::kEta);
        hm->AddHistogram(histClass, "MassD0region_TauxyzProj", "", false, 50, 1.7, 2.0, VarManager::kMass, 1000, -0.03, 0.03, VarManager::kVertexingTauxyzProjected);
      }
      if (subGroupStr.Contains("lambdac")) {
        hm->AddHistogram(histClass, "MassLambdacRegion", "", false, 50, 2.15, 2.4, VarManager::kMass);
        hm->AddHistogram(histClass, "MassLambdacRegion_Pt", "", false, 50, 2.15, 2.4, VarManager::kMass, 40, 0.0, 20.0, VarManager::kPt);
        hm->AddHistogram(histClass, "MassLambdacRegion_eta", "", false, 50, 2.15, 2.4, VarManager::kMass, 40, -2., 2., VarManager::kEta);
        hm->AddHistogram(histClass, "MassLambdacRegion_TauxyzProj", "", false, 50, 2.15, 2.4, VarManager::kMass, 1000, -0.03, 0.03, VarManager::kVertexingTauxyzProjected);
      }
    }
  }

  if (!groupStr.CompareTo("dilepton-track")) {
    if (subGroupStr.Contains("mixedevent")) { // for mixed event
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 40, 0.0, 20.0, VarManager::kPairMass, 40, 0.0, 20.0, VarManager::kPairPt);
      hm->AddHistogram(histClass, "Mass", "", false, 750, 0.0, 30.0, VarManager::kPairMass);
      hm->AddHistogram(histClass, "Pt", "", false, 750, 0.0, 30.0, VarManager::kPairPt);
    }
    if (subGroupStr.Contains("invmass")) {
      hm->AddHistogram(histClass, "Mass_Dilepton", "", false, 125, 0.0, 5.0, VarManager::kPairMassDau);
      hm->AddHistogram(histClass, "Mass_Hadron", "", false, 125, 0.0, 5.0, VarManager::kMassDau);
      hm->AddHistogram(histClass, "Delta_Mass", "", false, 125, 0.0, 5.0, VarManager::kDeltaMass);
      hm->AddHistogram(histClass, "Mass_Dilepton_Mass_Hadron", "", false, 125, 0.0, 5.0, VarManager::kPairMassDau, 125, 0.0, 5.0, VarManager::kMassDau);
      hm->AddHistogram(histClass, "Pt_Dilepton", "", false, 120, 0.0, 30.0, VarManager::kPairPtDau);
      hm->AddHistogram(histClass, "Pt_Track", "", false, 120, 0.0, 30.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Mass", "", false, 750, 0.0, 30.0, VarManager::kPairMass);
      hm->AddHistogram(histClass, "Pt", "", false, 750, 0.0, 30.0, VarManager::kPairPt);
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 40, 0.0, 20.0, VarManager::kPairMass, 40, 0.0, 20.0, VarManager::kPairPt);
      hm->AddHistogram(histClass, "Pt_Dilepton__Pt", "", false, 40, 0.0, 20.0, VarManager::kPairPtDau, 40, 0.0, 20.0, VarManager::kPairPt);
      hm->AddHistogram(histClass, "Pt_Track__Pt", "", false, 40, 0.0, 20.0, VarManager::kPt, 40, 0.0, 20.0, VarManager::kPairPt);
    }
    if (subGroupStr.Contains("vertexing")) {
      hm->AddHistogram(histClass, "UsedKF", "", false, 2, -0.5, 1.5, VarManager::kUsedKF);
      hm->AddHistogram(histClass, "KFMass", "", false, 750, 0.0, 30.0, VarManager::kKFMass);
      hm->AddHistogram(histClass, "Lz", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLz);
      hm->AddHistogram(histClass, "Lxy", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxy);
      hm->AddHistogram(histClass, "Lxyz", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyz);
      hm->AddHistogram(histClass, "Tauz", "", false, 4000, -0.01, 0.01, VarManager::kVertexingTauz);
      hm->AddHistogram(histClass, "Tauxy", "", false, 4000, -0.01, 0.01, VarManager::kVertexingTauxy);
      hm->AddHistogram(histClass, "LxyzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyzErr);
      hm->AddHistogram(histClass, "LzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLzErr);
      hm->AddHistogram(histClass, "TauzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingTauzErr);
      hm->AddHistogram(histClass, "VtxingProcCode", "", false, 10, 0.0, 10.0, VarManager::kVertexingProcCode);
      hm->AddHistogram(histClass, "VtxingChi2PCA", "", false, 100, 0.0, 10.0, VarManager::kVertexingChi2PCA);
      hm->AddHistogram(histClass, "LzProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLzProjected);
      hm->AddHistogram(histClass, "LxyProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyProjected);
      hm->AddHistogram(histClass, "LxyzProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyzProjected);
      hm->AddHistogram(histClass, "TauzProj", "", false, 4000, -0.5, 0.5, VarManager::kVertexingTauzProjected);
      hm->AddHistogram(histClass, "TauxyProj", "", false, 4000, -0.5, 0.5, VarManager::kVertexingTauxyProjected);
      hm->AddHistogram(histClass, "CosPointingAngle", "", false, 100, 0.0, 1.0, VarManager::kCosPointingAngle);
    }
    if (subGroupStr.Contains("correlation")) {
      hm->AddHistogram(histClass, "DeltaEta_DeltaPhi", "", false, 20, -2.0, 2.0, VarManager::kDeltaEta, 50, -8.0, 8.0, VarManager::kDeltaPhi);
      hm->AddHistogram(histClass, "DeltaEta_DeltaPhiSym", "", false, 20, -2.0, 2.0, VarManager::kDeltaEta, 50, -8.0, 8.0, VarManager::kDeltaPhiSym);
    }
    if (subGroupStr.Contains("opencharm")) {
      hm->AddHistogram(histClass, "Delta_Mass_DstarD0region", "", false, 50, 0.14, 0.16, VarManager::kDeltaMass);
    }
  }

  if (!groupStr.CompareTo("dilepton-hadron-array-correlation")) {
    const int kInvMassBins = 500;
    double InvMassBinLims[kInvMassBins + 1];
    for (int i = 0; i <= kInvMassBins; i++)
      InvMassBinLims[i] = 0 + i * 0.01;

    const int kDelEtaBins = 20;
    double DelEtaBinLims[kDelEtaBins + 1];
    for (int i = 0; i <= kDelEtaBins; i++)
      DelEtaBinLims[i] = -2 + i * 0.2;

    const int kDelPhiBins = 50;
    double DelPhiBinLims[kDelPhiBins + 1];
    for (int i = 0; i <= kDelPhiBins; i++)
      DelPhiBinLims[i] = -8 + i * 0.32;

    const int kPtBins = 20;
    double PtBinLims[kPtBins + 1];
    for (int i = 0; i <= kPtBins; i++)
      PtBinLims[i] = 0 + i * 1;

    TArrayD nJPsiHadCorr[4];
    nJPsiHadCorr[0] = TArrayD(kInvMassBins + 1, InvMassBinLims);
    nJPsiHadCorr[1] = TArrayD(kDelEtaBins + 1, DelEtaBinLims);
    nJPsiHadCorr[2] = TArrayD(kDelPhiBins + 1, DelPhiBinLims);
    nJPsiHadCorr[3] = TArrayD(kPtBins + 1, PtBinLims);

    int varsJPsiHadCorr[4] = {VarManager::kPairMassDau, VarManager::kDeltaEta, VarManager::kDeltaPhi, VarManager::kPairPtDau};
    hm->AddHistogram(histClass, "InvMass_DelEta_DelPhi", "", 4, varsJPsiHadCorr, nJPsiHadCorr);
  }
  if (!groupStr.CompareTo("dilepton-charmhadron")) {
    if (subGroupStr.EqualTo("jpsitomumu")) {
      hm->AddHistogram(histClass, "hMassVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 300, 2.f, 5.f, VarManager::kMass);
      hm->AddHistogram(histClass, "hRapVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 50, -4.5f, -2.0f, VarManager::kRap);
      hm->AddHistogram(histClass, "hPhiVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 180, -constants::math::PI, constants::math::PI, VarManager::kPhi);
    } else if (subGroupStr.EqualTo("jpsitoee")) {
      hm->AddHistogram(histClass, "hMassVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 300, 2.f, 5.f, VarManager::kMass);
      hm->AddHistogram(histClass, "hRapVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 60, -1.5f, 1.5f, VarManager::kRap);
      hm->AddHistogram(histClass, "hPhiVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 180, -constants::math::PI, constants::math::PI, VarManager::kPhi);
    } else if (subGroupStr.EqualTo("dmeson")) {
      hm->AddHistogram(histClass, "hMassVsPtVsBdtDmeson", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 200, 1.7f, 2.1f, VarManager::kMassCharmHadron);
      hm->AddHistogram(histClass, "hRapVsPtVsBdtDmeson", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 60, -1.5f, 1.5f, VarManager::kRapCharmHadron);
      hm->AddHistogram(histClass, "hPhiVsPtVsBdtDmeson", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 180, 0., 2 * constants::math::PI, VarManager::kPhiCharmHadron);
    } else if (subGroupStr.EqualTo("jpsitomumudmeson")) {
      hm->AddHistogram(histClass, "hMassVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 300, 2.f, 5.f, VarManager::kMass);
      hm->AddHistogram(histClass, "hRapVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 50, -4.5f, -2.0f, VarManager::kRap);
      hm->AddHistogram(histClass, "hPhiVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 180, -constants::math::PI, constants::math::PI, VarManager::kPhi);
      hm->AddHistogram(histClass, "hMassVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 200, 1.7f, 2.1f, VarManager::kMassCharmHadron);
      hm->AddHistogram(histClass, "hRapVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 60, -1.5f, 1.5f, VarManager::kRapCharmHadron);
      hm->AddHistogram(histClass, "hPhiVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 180, 0., 2 * constants::math::PI, VarManager::kPhiCharmHadron);
    } else if (subGroupStr.EqualTo("jpsitoeedmeson")) {
      hm->AddHistogram(histClass, "hMassVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 300, 2.f, 5.f, VarManager::kMass);
      hm->AddHistogram(histClass, "hRapVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 60, -1.5f, 1.5f, VarManager::kRap);
      hm->AddHistogram(histClass, "hPhiVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 180, -constants::math::PI, constants::math::PI, VarManager::kPhi);
      hm->AddHistogram(histClass, "hMassVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 200, 1.7f, 2.1f, VarManager::kMassCharmHadron);
      hm->AddHistogram(histClass, "hRapVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 60, -1.5f, 1.5f, VarManager::kRapCharmHadron);
      hm->AddHistogram(histClass, "hPhiVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 180, 0., 2 * constants::math::PI, VarManager::kPhiCharmHadron);
    }
  }
  if (!groupStr.CompareTo("dilepton-dihadron")) {
    if (subGroupStr.EqualTo("xtojpsipipi")) {
      hm->AddHistogram(histClass, "hMass_X3872", "", false, 1000, 3.0, 5.0, VarManager::kQuadMass);
      hm->AddHistogram(histClass, "hPt_X3872", "", false, 200, 0.0, 20.0, VarManager::kQuadPt);
      hm->AddHistogram(histClass, "hMass_Pt_X3872", "", false, 1000, 3.0, 5.0, VarManager::kQuadMass, 150, 0.0, 15.0, VarManager::kQuadPt);
      hm->AddHistogram(histClass, "hCostheta_Jpsi_Dihadron", "", false, 100, -1.0, 1.0, VarManager::kCosthetaDileptonDitrack);
      hm->AddHistogram(histClass, "hPtDilepton_PtDihadron", "", false, 200, 0, 20, VarManager::kPairPt, 100, 0, 10, VarManager::kDitrackPt);
      hm->AddHistogram(histClass, "hPtDilepton_MassDihadron", "", false, 200, 0, 20, VarManager::kPairPt, 100, 3.0, 5.0, VarManager::kDitrackMass);
    }
  }
  if (!groupStr.CompareTo("dilepton-photon-mass")) {
    hm->AddHistogram(histClass, "Mass_Dilepton", "", false, 500, 0.0, 5.0, VarManager::kPairMassDau);
    hm->AddHistogram(histClass, "Mass_Photon", "", false, 500, 0.0, 0.1, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Mass_Dilepton_Mass_Photon", "", false, 250, 0.0, 5.0, VarManager::kPairMassDau, 250, 0.0, 5.0, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Pt_Dilepton", "", false, 2000, 0.0, 20.0, VarManager::kPairPtDau);
    hm->AddHistogram(histClass, "Pt_Photon", "", false, 500, 0.0, 5.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Mass_DileptonPhoton", "", false, 4500, 0.0, 4.5, VarManager::kPairMass);
    hm->AddHistogram(histClass, "Pt_DileptonPhoton", "", false, 2000, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kPairMass, 200, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "DeltaMass", "", false, 1500, 0.0, 1.5, VarManager::kDeltaMass);
    hm->AddHistogram(histClass, "DeltaMass_pt", "", false, 1000, 0.0, 1.0, VarManager::kDeltaMass, 3000, 0.0, 30.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "Eta_Pt", "", false, 100, -2.0, 2.0, VarManager::kPairEta, 200, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "DeltaMass_Jpsi", "", false, 1500, 3, 4.5, VarManager::kDeltaMass_jpsi);
    hm->AddHistogram(histClass, "Rapidity", "", false, 400, -4.0, 4.0, VarManager::kRap);
    hm->AddHistogram(histClass, "Eta_Pt_Dilepton", "", false, 100, -2.0, 2.0, VarManager::kDeltaEta, 200, 0.0, 20.0, VarManager::kPairPtDau);
    hm->AddHistogram(histClass, "Eta_Pt_Photon", "", false, 100, -2.0, 2.0, VarManager::kEta, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Eta_Pt_lepton1", "", false, 100, -2.0, 2.0, VarManager::kEta1, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Eta_Pt_lepton2", "", false, 100, -2.0, 2.0, VarManager::kEta2, 200, 0.0, 20.0, VarManager::kPt);
  }

  if (!groupStr.CompareTo("photon")) {
    hm->AddHistogram(histClass, "Pt_Photon", "p_{T} distribution", false, 4500, 0.0, 4.5, VarManager::kPt);
    hm->AddHistogram(histClass, "Eta", "#eta distribution", false, 500, -5.0, 5.0, VarManager::kEta);
    hm->AddHistogram(histClass, "Eta_Pt", "", false, 100, -2.0, 2.0, VarManager::kEta, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Phi", "#varphi distribution", false, 500, -2. * TMath::Pi(), 2. * TMath::Pi(), VarManager::kPhi);
    hm->AddHistogram(histClass, "Mass_Photon", "", false, 500, 0.0, 0.1, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kMassDau, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Rapidity", "", false, 400, -4.0, 4.0, VarManager::kRap);
  }
}
