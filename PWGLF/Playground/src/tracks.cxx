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
///
/// \brief this is a playground task to look at track PID
/// \author 
/// \since

#include <iostream>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyTracksRun3 = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe>;
using MyTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe>;

void GetChargeRatio(HistogramRegistry reg, std::shared_ptr<TH1> histo1, std::shared_ptr<TH1> histo2, std::shared_ptr<TH1> histo3, std::shared_ptr<TH1> histo4) {
  const char *label[4] = {"all tracks", "protons", "nuclei", "triton"};
  const Double_t ChargeRatio1 = histo1->GetBinContent(1)/histo1->GetBinContent(2);
  const Double_t ChargeRatio2 = histo2->GetBinContent(1)/histo2->GetBinContent(2);
  const Double_t ChargeRatio3 = histo3->GetBinContent(1)/histo3->GetBinContent(2);
  const Double_t ChargeRatio4 = histo4->GetBinContent(1)/histo4->GetBinContent(2);
  for (Int_t i=0; i<4; i++) {
    reg.get<TH1>(HIST("Tracks/hChargeRatio"))->GetXaxis()->SetBinLabel(i+1, label[i]);
  }
  reg.get<TH1>(HIST("Tracks/hChargeRatio"))->SetBinContent(1, ChargeRatio1);
  reg.get<TH1>(HIST("Tracks/hChargeRatio"))->SetBinContent(2, ChargeRatio2);
  reg.get<TH1>(HIST("Tracks/hChargeRatio"))->SetBinContent(3, ChargeRatio3);
  reg.get<TH1>(HIST("Tracks/hChargeRatio"))->SetBinContent(4, ChargeRatio4);
}

struct trackchecks {
  /// histogram configurables
  Configurable<int> nBinsLow{"nBinsLow", 200, "N bins in all histos"};
  Configurable<int> nBinsHigh{"nBinsHigh", 800, "N bins in all histos"};
  Configurable<double> minXTPCsignal{"minXTPCsignal", 0.1, "X axis min for TPCsignal histo"};
  Configurable<double> maxXTPCsignal{"maxXTPCsignal", 2.5, "X axis max for TPCsignal histo"};
  Configurable<double> maxPt{"maxPt", 2.5, "max for pT axis"};
  
  /// track configurables
  Configurable<float> eta{"eta", 0.8, "eta"}; // expression table
  Configurable<float> tpcNClsFound{"tpcNClsFound", 70, "Number of TPC clusters"}; // dynamic
  Configurable<float> tpcCrossedRowsOverFindableCls{"tpcCrossedRowsOverFindableCls", 0.8, "Ratio crossed rows over findable"}; // dynamic

  /// Filters (Cannot filter on dynamic columns -> only filter on rapidity (expression column))
  Filter rapFilter = nabs(aod::track::eta) < eta;
  
  /// histogram registry
  HistogramRegistry registry{
    "registry",
    {
      /// event
      {"Events/hVertexX", "PV position x; Vtx_{x} (cm); entries", {HistType::kTH1F, {{nBinsLow, -0.8, 0.8}}}},
      {"Events/hVertexY", "PV position y; Vtx_{y} (cm); entries", {HistType::kTH1F, {{nBinsLow, -0.8, 0.8}}}},
      {"Events/hVertexZ", "PV position z; Vtx_{z} (cm); entries", {HistType::kTH1F, {{nBinsLow, -12., 12.}}}},
      {"Events/hVertexXY", "PV position xy; Vtx_{x} (cm); Vtx_{y} (cm)", {HistType::kTH2F, {{nBinsHigh, -0.8, 0.8}, {nBinsHigh, -0.8, 0.8}}}},
      {"Events/hVertexNcontrib", "Real Number of contributors; real N contributors; entries", {HistType::kTH1F, {{80, 0., 80.}}}},
      {"Events/hMultiplicity", "Event multiplicity; multiplicity; entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"Events/hNprotons", "Number of selected protons per event; multiplicity; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
      {"Events/hNnuclei", "Number of selected deuterons and tritons per event; multiplicity; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
      {"Events/hTritons", "Number of selected tritons per event; multiplicity; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
      /// tracks
      {"Tracks/heta", "Track eta; #eta; entries", {HistType::kTH1F, {{nBinsLow, -1., 1.}}}},
      {"Tracks/hphi", "Track phi; #phi; entries", {HistType::kTH1F, {{nBinsLow, 0., 6.5}}}},
      {"Tracks/hsign", "Track charge sign; charge sign; entries", {HistType::kTH1F, {{2, -2., 2.}}}},
      {"Tracks/hpt", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{nBinsLow, 0., 10.}}}},
      {"Tracks/hDCAxy", "Track DCA xy to the PV; DCA_{xy} (cm); entries", {HistType::kTH1F, {{nBinsHigh, -2., 2.}}}},
      {"Tracks/hDCAz", "Track DCA z to the PV; DCA_{xy} (cm); entries", {HistType::kTH1F, {{nBinsHigh, -5., 5.}}}},
      {"Tracks/hNSigmaTPCPi", "nSigma pion; #it{p} (GeV/#it{c}); n#sigma(#pi)", {HistType::kTH2F, {{nBinsHigh, 0., maxPt}, {nBinsHigh, -10., 10.}}}},
      {"Tracks/hNSigmaTPCKa", "nSigma kaon; #it{p} (GeV/#it{c}); n#sigma(K)", {HistType::kTH2F, {{nBinsHigh, 0., maxPt}, {nBinsHigh, -10., 10.}}}},
      {"Tracks/hNSigmaTPCPr", "nSigma proton; #it{p} (GeV/#it{c}); n#sigma(p)", {HistType::kTH2F, {{nBinsHigh, 0., maxPt}, {nBinsHigh, -10., 10.}}}},
      {"Tracks/hNSigmaTPCDe", "nSigma deuteron; #it{p} (GeV/#it{c}); n#sigma(d)", {HistType::kTH2F, {{nBinsHigh, 0., maxPt}, {nBinsHigh, -10., 10.}}}},
      {"Tracks/hNSigmaTPCTr", "nSigma triton; #it{p} (GeV/#it{c}); n#sigma(t)", {HistType::kTH2F, {{nBinsHigh, 0., maxPt}, {nBinsHigh, -10., 10.}}}},
      {"Tracks/hNSigmaTPCHe", "nSigma helium; #it{p} (GeV/#it{c}); n#sigma(He)", {HistType::kTH2F, {{nBinsHigh, 0., 20.0}, {nBinsHigh, -30., 30.}}}},
      {"Tracks/hTPCsignal", "TPC signal; #it{p} (GeV/#it{c}); TPC signal (a. u.)", {HistType::kTH2F, {{3000, minXTPCsignal, maxXTPCsignal}, {nBinsHigh, 20., 500.}}}},
      {"Tracks/hTPCsignalCharge", "TPC signal; #it{p}/#it{z} (GeV/#it{c}); TPC signal (a. u.)", {HistType::kTH2F, {{3000, -1.0*maxXTPCsignal, maxXTPCsignal}, {nBinsHigh, 20., 500.}}}},
      /// protons
      {"Protons/heta", "Track eta; #eta; entries", {HistType::kTH1F, {{nBinsLow, -1., 1.}}}},
      {"Protons/hphi", "Track phi; #phi; entries", {HistType::kTH1F, {{nBinsLow, 0., 6.5}}}},
      {"Protons/hsign", "Track charge sign; charge sign; entries", {HistType::kTH1F, {{2, -2., 2.}}}},
      {"Protons/hDCAxy", "Track DCA xy to the PV; DCA_{xy} (cm); entries", {HistType::kTH1F, {{nBinsHigh, -2., 2.}}}},
      {"Protons/hDCAz", "Track DCA z to the PV; DCA_{xy} (cm); entries", {HistType::kTH1F, {{nBinsHigh, -5., 5.}}}},
      {"Protons/hpt", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{nBinsLow, 0., 10.}}}},
      {"Protons/hNSigmaTPCPr", "nSigma proton; #it{p} (GeV/#it{c}); n#sigma(p)", {HistType::kTH2F, {{nBinsHigh, 0., maxPt}, {nBinsHigh, -10., 10.}}}},
      {"Protons/hTPCsignal", "TPC signal; #it{p} (GeV/#it{c}); TPC signal (a. u.)", {HistType::kTH2F, {{3000, minXTPCsignal, maxXTPCsignal}, {nBinsHigh, 20., 500.}}}},
      /// nuclei (deuterons and triton)
      {"Nuclei/heta", "Track eta; #eta; entries", {HistType::kTH1F, {{nBinsLow, -1., 1.}}}},
      {"Nuclei/hphi", "Track phi; #phi; entries", {HistType::kTH1F, {{nBinsLow, 0., 6.5}}}},
      {"Nuclei/hsign", "Track charge sign; charge sign; entries", {HistType::kTH1F, {{2, -2., 2.}}}},
      {"Nuclei/hDCAxy", "Track DCA xy to the PV; DCA_{xy} (cm); entries", {HistType::kTH1F, {{nBinsHigh, -2., 2.}}}},
      {"Nuclei/hDCAz", "Track DCA z to the PV; DCA_{xy} (cm); entries", {HistType::kTH1F, {{nBinsHigh, -5., 5.}}}},
      {"Nuclei/hpt", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{nBinsLow, 0., 10.}}}},
      {"Nuclei/hNSigmaTPCDe", "nSigma deuteron; #it{p} (GeV/#it{c}); n#sigma(d)", {HistType::kTH2F, {{nBinsHigh, 0., maxPt}, {nBinsHigh, -10., 10.}}}},
      {"Nuclei/hTPCsignal", "TPC signal; #it{p} (GeV/#it{c}); TPC signal (a. u.)", {HistType::kTH2F, {{3000, minXTPCsignal, maxXTPCsignal}, {nBinsHigh, 20., 500.}}}},
      /// tritons selected with nsigma
      {"Triton/heta", "Track eta; #eta; entries", {HistType::kTH1F, {{nBinsLow, -1., 1.}}}},
      {"Triton/hphi", "Track phi; #phi; entries", {HistType::kTH1F, {{nBinsLow, 0., 6.5}}}},
      {"Triton/hsign", "Track charge sign; charge sign; entries", {HistType::kTH1F, {{2, -2., 2.}}}},
      {"Triton/hDCAxy", "Track DCA xy to the PV; DCA_{xy} (cm); entries", {HistType::kTH1F, {{nBinsHigh, -2., 2.}}}},
      {"Triton/hDCAz", "Track DCA z to the PV; DCA_{xy} (cm); entries", {HistType::kTH1F, {{nBinsHigh, -5., 5.}}}},
      {"Triton/hpt", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{nBinsLow, 0., 10.}}}},
      {"Triton/hNSigmaTPCTr", "nSigma triton; #it{p} (GeV/#it{c}); n#sigma(tr)", {HistType::kTH2F, {{nBinsHigh, 0., maxPt}, {nBinsHigh, -10., 10.}}}},
      {"Triton/hTPCsignal", "TPC signal; #it{p} (GeV/#it{c}); TPC signal (a. u.)", {HistType::kTH2F, {{3000, minXTPCsignal, maxXTPCsignal}, {nBinsHigh, 20., 500.}}}}
      
    }
  };

  void processEvents(soa::Join<aod::Collisions, aod::Mults, aod::EvSels>::iterator const& collision) {
    /// basic event selection
    if (!collision.sel8() || TMath::Abs(collision.posZ()) >= 10.) return;
    /// fill event histograms
    registry.get<TH1>(HIST("Events/hVertexNcontrib"))->Fill(collision.multNTracksPV());
    if (collision.multNTracksPV() <= 1) return;
    registry.get<TH1>(HIST("Events/hVertexX"))->Fill(collision.posX());
    registry.get<TH1>(HIST("Events/hVertexY"))->Fill(collision.posY());
    registry.get<TH1>(HIST("Events/hVertexZ"))->Fill(collision.posZ());
    registry.get<TH2>(HIST("Events/hVertexXY"))->Fill(collision.posX(), collision.posY());
  }
  PROCESS_SWITCH(trackchecks, processEvents, "Process events, plot event properties", true);


  void processTracks(soa::Join<aod::Collisions, aod::Mults, aod::EvSels>::iterator const& collision, soa::Filtered<MyTracksRun3> const& tracks) {
    /// basic event selection
    if (!collision.sel8() || TMath::Abs(collision.posZ()) >= 10. || collision.multNTracksPV() <= 1) return;
    /// loop over tracks
    int nTracks = 0; /// counter for multiplicity
    int nProtons = 0; /// counter for protons
    int nNuclei = 0; /// counter for deuterons and tritons
    for (auto& track : tracks) {
      if (!track.hasITS() || !track.hasTPC() || track.tpcNClsFound() < tpcNClsFound || track.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableCls) continue; 
      registry.get<TH1>(HIST("Tracks/heta"))->Fill(track.eta());
      registry.get<TH1>(HIST("Tracks/hphi"))->Fill(track.phi());
      registry.get<TH1>(HIST("Tracks/hsign"))->Fill(track.sign());
      registry.get<TH1>(HIST("Tracks/hpt"))->Fill(track.pt());
      registry.get<TH1>(HIST("Tracks/hDCAxy"))->Fill(track.dcaXY());
      registry.get<TH1>(HIST("Tracks/hDCAz"))->Fill(track.dcaZ());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCPi"))->Fill(track.tpcInnerParam(), track.tpcNSigmaPi());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCKa"))->Fill(track.tpcInnerParam(), track.tpcNSigmaKa());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCPr"))->Fill(track.tpcInnerParam(), track.tpcNSigmaPr());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCDe"))->Fill(track.tpcInnerParam(), track.tpcNSigmaDe());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCTr"))->Fill(track.tpcInnerParam(), track.tpcNSigmaTr());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCHe"))->Fill(track.tpcInnerParam(), track.tpcNSigmaHe());
      registry.get<TH2>(HIST("Tracks/hTPCsignal"))->Fill(track.tpcInnerParam(), track.tpcSignal());
      registry.get<TH2>(HIST("Tracks/hTPCsignalCharge"))->Fill(track.tpcInnerParam()/track.sign(), track.tpcSignal());
      nTracks++;
      /// select protons from dE/dx
      if ((track.pt()<0.45 && track.pt()>0.25 && track.tpcSignal()>350 && track.tpcSignal()<500) || (track.pt()<0.5 && track.pt()>0.29 && track.tpcSignal()>250 && track.tpcSignal()<=350) || (track.pt()<0.6 && track.pt()>0.35 && track.tpcSignal()>170 && track.tpcSignal()<=250)) {
        registry.get<TH1>(HIST("Protons/heta"))->Fill(track.eta());
        registry.get<TH1>(HIST("Protons/hphi"))->Fill(track.phi());
        registry.get<TH1>(HIST("Protons/hsign"))->Fill(track.sign());
        registry.get<TH1>(HIST("Protons/hpt"))->Fill(track.pt());
        registry.get<TH1>(HIST("Protons/hDCAxy"))->Fill(track.dcaXY());
        registry.get<TH1>(HIST("Protons/hDCAz"))->Fill(track.dcaZ());
        registry.get<TH2>(HIST("Protons/hNSigmaTPCPr"))->Fill(track.tpcInnerParam(), track.tpcNSigmaPr());
        registry.get<TH2>(HIST("Protons/hTPCsignal"))->Fill(track.tpcInnerParam(), track.tpcSignal());
        nProtons++;
      }
      /// select deuterons and tritons from dE/dx
      if (track.pt()<1.5 && ((track.pt()>0.45 && track.tpcSignal()>320 && track.tpcSignal()<500) || (track.pt()>0.6 && track.tpcSignal()>210 && track.tpcSignal()<=320) || (track.pt()>0.8 && track.tpcSignal()>150 && track.tpcSignal()<=210))) {
        registry.get<TH1>(HIST("Nuclei/heta"))->Fill(track.eta());
        registry.get<TH1>(HIST("Nuclei/hphi"))->Fill(track.phi());
        registry.get<TH1>(HIST("Nuclei/hsign"))->Fill(track.sign());
        registry.get<TH1>(HIST("Nuclei/hDCAxy"))->Fill(track.dcaXY());
        registry.get<TH1>(HIST("Nuclei/hDCAz"))->Fill(track.dcaZ());
        registry.get<TH1>(HIST("Nuclei/hpt"))->Fill(track.pt());
        registry.get<TH2>(HIST("Nuclei/hNSigmaTPCDe"))->Fill(track.tpcInnerParam(), track.tpcNSigmaDe());
        registry.get<TH2>(HIST("Nuclei/hTPCsignal"))->Fill(track.tpcInnerParam(), track.tpcSignal());
        nNuclei++;
      }
    } /// end track loop
    // LOG(info) << "Tracks: " << nTracks;
    // LOG(info) << "Protons: " << nProtons;
    // LOG(info) << "Deuterons: " << nDeuterons;
    // LOG(info) << "Tritons: " << nTritons;

    /// fill multiplicity histograms
    registry.get<TH1>(HIST("Events/hMultiplicity"))->Fill(nTracks);
    registry.get<TH1>(HIST("Events/hNprotons"))->Fill(nProtons);
    registry.get<TH1>(HIST("Events/hNnuclei"))->Fill(nNuclei);

  } /// end process
  PROCESS_SWITCH(trackchecks, processTracks, "Process tracks, PID with dE/dx signal", false);

  void processTracksPID(soa::Join<aod::Collisions, aod::Mults, aod::EvSels>::iterator const& collision, soa::Filtered<MyTracksRun3> const& tracks) {
    /// basic event selection
    if (!collision.sel8() || TMath::Abs(collision.posZ()) >= 10. || collision.multNTracksPV() <= 1) return;
    /// loop over tracks
    int nTracks = 0; /// counter for multiplicity
    int nProtons = 0; /// counter for protons
    int nDeuterons = 0; /// counter for deuterons
    int nTritons = 0; /// counter for tritons
    for (auto& track : tracks) {
      if (!track.hasITS() || !track.hasTPC() || track.tpcNClsFound() < tpcNClsFound || track.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableCls) continue; 
      registry.get<TH1>(HIST("Tracks/heta"))->Fill(track.eta());
      registry.get<TH1>(HIST("Tracks/hphi"))->Fill(track.phi());
      registry.get<TH1>(HIST("Tracks/hsign"))->Fill(track.sign());
      registry.get<TH1>(HIST("Tracks/hpt"))->Fill(track.pt());
      registry.get<TH1>(HIST("Tracks/hDCAxy"))->Fill(track.dcaXY());
      registry.get<TH1>(HIST("Tracks/hDCAz"))->Fill(track.dcaZ());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCPi"))->Fill(track.tpcInnerParam(), track.tpcNSigmaPi());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCKa"))->Fill(track.tpcInnerParam(), track.tpcNSigmaKa());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCPr"))->Fill(track.tpcInnerParam(), track.tpcNSigmaPr());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCDe"))->Fill(track.tpcInnerParam(), track.tpcNSigmaDe());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCTr"))->Fill(track.tpcInnerParam(), track.tpcNSigmaTr());
      registry.get<TH2>(HIST("Tracks/hNSigmaTPCHe"))->Fill(track.tpcInnerParam(), track.tpcNSigmaHe());
      registry.get<TH2>(HIST("Tracks/hTPCsignal"))->Fill(track.tpcInnerParam(), track.tpcSignal());
      registry.get<TH2>(HIST("Tracks/hTPCsignalCharge"))->Fill(track.tpcInnerParam()/track.sign(), track.tpcSignal());
      nTracks++;
      /// select protons with nsigma
      if (track.tpcInnerParam()<1.2 && track.tpcInnerParam()>0. && track.tpcNSigmaPr()<4. && track.tpcNSigmaPr()>-2.) {
        registry.get<TH1>(HIST("Protons/heta"))->Fill(track.eta());
        registry.get<TH1>(HIST("Protons/hphi"))->Fill(track.phi());
        registry.get<TH1>(HIST("Protons/hsign"))->Fill(track.sign());
        registry.get<TH1>(HIST("Protons/hpt"))->Fill(track.pt());
        registry.get<TH1>(HIST("Protons/hDCAxy"))->Fill(track.dcaXY());
        registry.get<TH1>(HIST("Protons/hDCAz"))->Fill(track.dcaZ());
        registry.get<TH2>(HIST("Protons/hNSigmaTPCPr"))->Fill(track.tpcInnerParam(), track.tpcNSigmaPr());
        registry.get<TH2>(HIST("Protons/hTPCsignal"))->Fill(track.tpcInnerParam(), track.tpcSignal());
        nProtons++;
      }
      /// select deuterons with nsigma
      if (track.tpcInnerParam()<1.5 && track.tpcInnerParam()>0. && track.tpcNSigmaDe()<4. && track.tpcNSigmaDe()>-2.) {
        registry.get<TH1>(HIST("Nuclei/heta"))->Fill(track.eta());
        registry.get<TH1>(HIST("Nuclei/hphi"))->Fill(track.phi());
        registry.get<TH1>(HIST("Nuclei/hsign"))->Fill(track.sign());
        registry.get<TH1>(HIST("Nuclei/hDCAxy"))->Fill(track.dcaXY());
        registry.get<TH1>(HIST("Nuclei/hDCAz"))->Fill(track.dcaZ());
        registry.get<TH1>(HIST("Nuclei/hpt"))->Fill(track.pt());
        registry.get<TH2>(HIST("Nuclei/hNSigmaTPCDe"))->Fill(track.tpcInnerParam(), track.tpcNSigmaDe());
        registry.get<TH2>(HIST("Nuclei/hTPCsignal"))->Fill(track.tpcInnerParam(), track.tpcSignal());
        nDeuterons++;
      }
      /// select tritons with nsigma
      if (track.tpcInnerParam()<2.0 && track.tpcInnerParam()>0. && track.tpcNSigmaTr()<4. && track.tpcNSigmaTr()>-2.5) {
        registry.get<TH1>(HIST("Triton/heta"))->Fill(track.eta());
        registry.get<TH1>(HIST("Triton/hphi"))->Fill(track.phi());
        registry.get<TH1>(HIST("Triton/hsign"))->Fill(track.sign());
        registry.get<TH1>(HIST("Triton/hDCAxy"))->Fill(track.dcaXY());
        registry.get<TH1>(HIST("Triton/hDCAz"))->Fill(track.dcaZ());
        registry.get<TH1>(HIST("Triton/hpt"))->Fill(track.pt());
        registry.get<TH2>(HIST("Triton/hNSigmaTPCTr"))->Fill(track.tpcInnerParam(), track.tpcNSigmaTr());
        registry.get<TH2>(HIST("Triton/hTPCsignal"))->Fill(track.tpcInnerParam(), track.tpcSignal());
        nTritons++;
      }
    } /// end track loop
    // LOG(info) << "Tracks: " << nTracks;
    // LOG(info) << "Protons: " << nProtons;
    // LOG(info) << "Deuterons: " << nDeuterons;
    // LOG(info) << "Tritons: " << nTritons;

    /// fill multiplicity histograms
    registry.get<TH1>(HIST("Events/hMultiplicity"))->Fill(nTracks);
    registry.get<TH1>(HIST("Events/hNprotons"))->Fill(nProtons);
    registry.get<TH1>(HIST("Events/hNnuclei"))->Fill(nDeuterons);
    registry.get<TH1>(HIST("Events/hTritons"))->Fill(nTritons);

  } /// end process
  PROCESS_SWITCH(trackchecks, processTracksPID, "Process tracks, PID with nsigma", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<trackchecks>(cfgc)
  };
}