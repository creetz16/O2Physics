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
/// \brief QA task to study properties of propagated tracks

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// allows for candidate-by-candidate comparison using Cascade to []CascData link table
using TracksLabeled = soa::Join<aod::StoredTracks, aod::StoredTracksCov, aod::TracksDCA, aod::TracksDCACov, aod::TracksExtra, aod::McTrackLabels>;
using Tracks = soa::Join<aod::StoredTracks, aod::StoredTracksCov, aod::TracksDCA, aod::TracksDCACov, aod::TracksExtra>;

struct lfpropStudy {

  ConfigurableAxis axisDCAxy{"axisDCAxy", {1000, -10.f, 10.f}, "DCA_{xy} (cm)"};
  ConfigurableAxis axisDCAz{"axisDCAz", {1000, -10.f, 10.f}, "DCA_{z} (cm)"};
  ConfigurableAxis axisMom{"axisMom", {1000, 0.0f, 10.f}, "p (GeV/c)"};

  Configurable<float> d_pTMin{"d_pTMin", 0.3, "minimum track mommentum"};
  Configurable<float> d_TPCrowsMin{"d_TPCrowsMin", 70, "minimum number of TPC crossed rows"};
  Configurable<float> d_TPCrowsOverFindMin{"d_TPCrowsOverFindMin", 0.8, "minimum for ratio of TPC crossed rows over findable"};
  
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Filter collisionFilter = (aod::evsel::sel8 == true);
  Filter etaFilter = (nabs(aod::track::eta) < 0.8f);
  Filter ptFilter = (aod::track::pt > d_pTMin);

  void init(InitContext const&)
  {
    histos.add("hEventCounter", "hEventCounter", kTH1F, {{1, 0.0f, 1.0f}});

    histos.add("hDCAxy", "hDCAxy", kTH1F, {axisDCAxy});
    histos.add("hDCAz", "hDCAz", kTH1F, {axisDCAz});
    histos.add("hPx", "hPx", kTH1F, {axisMom});
    histos.add("hPy", "hPy", kTH1F, {axisMom});
    histos.add("hPz", "hPz", kTH1F, {axisMom});

    histos.add("hDCAxyEl", "hDCAxyEl", kTH1F, {axisDCAxy});
    histos.add("hDCAzEl", "hDCAzEl", kTH1F, {axisDCAz});
    histos.add("hDCAxyPi", "hDCAxyPi", kTH1F, {axisDCAxy});
    histos.add("hDCAzPi", "hDCAzPi", kTH1F, {axisDCAz});
    histos.add("hDCAxyKa", "hDCAxyKa", kTH1F, {axisDCAxy});
    histos.add("hDCAzKa", "hDCAzKa", kTH1F, {axisDCAz});
    histos.add("hDCAxyPr", "hDCAxyPr", kTH1F, {axisDCAxy});
    histos.add("hDCAzPr", "hDCAzPr", kTH1F, {axisDCAz});
    histos.add("hDCAxyDe", "hDCAxyDe", kTH1F, {axisDCAxy});
    histos.add("hDCAzDe", "hDCAzDe", kTH1F, {axisDCAz});
    histos.add("hDCAxyTr", "hDCAxyTr", kTH1F, {axisDCAxy});
    histos.add("hDCAzTr", "hDCAzTr", kTH1F, {axisDCAz});
    histos.add("hDCAxyHe", "hDCAxyHe", kTH1F, {axisDCAxy});
    histos.add("hDCAzHe", "hDCAzHe", kTH1F, {axisDCAz});
    histos.add("hDCAxyAl", "hDCAxyAl", kTH1F, {axisDCAxy});
    histos.add("hDCAzAl", "hDCAzAl", kTH1F, {axisDCAz});
  }

  void processData(soa::Filtered<aod::Collision> const& collision, Tracks const& Tracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : Tracks) {

      // track selection
      if (track.hasITS()) continue; // only look at tracks which start outside the ITS
      if (track.tpcNClsCrossedRows() < d_TPCrowsMin || track.tpcCrossedRowsOverFindableCls() < d_TPCrowsOverFindMin) continue;

      histos.fill(HIST("hPx"), track.px());
      histos.fill(HIST("hPy"), track.py());
      histos.fill(HIST("hPz"), track.pz());
      histos.fill(HIST("hDCAxy"), track.dcaXY());
      histos.fill(HIST("hDCz"), track.dcaZ());
    }
  }
  PROCESS_SWITCH(lfpropStudy, processData, "process data", true);

  void processMC(soa::Filtered<aod::Collision> const& collision, TracksLabeled const& Tracks, aod::McParticles const& particlesMC)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : Tracks) {

      if (!track.has_mcParticle() || track.mcParticleId() <= -1 || track.mcParticleId() > particlesMC.size()) continue; // only look at tracks which had right PID assgined
      if (track.hasITS()) continue; // only look at tracks which start outside the ITS

      if (track.mcParticle().pdgCode() == kElectron) {
        histos.fill(HIST("hDCAxyEl"), track.dcaXY());
        histos.fill(HIST("hDCAzEl"), track.dcaZ());
      } else if (track.mcParticle().pdgCode() == kPiPlus) {
        histos.fill(HIST("hDCAxyPi"), track.dcaXY());
        histos.fill(HIST("hDCAzPi"), track.dcaZ());
      } else if (track.mcParticle().pdgCode() == kKPlus) {
        histos.fill(HIST("hDCAxyKa"), track.dcaXY());
        histos.fill(HIST("hDCAzKa"), track.dcaZ());
      } else if (track.mcParticle().pdgCode() == kProton) {
        histos.fill(HIST("hDCAxyPr"), track.dcaXY());
        histos.fill(HIST("hDCAzPr"), track.dcaZ());
      } else if (track.mcParticle().pdgCode() == 1000010020) {
        histos.fill(HIST("hDCAxyDe"), track.dcaXY());
        histos.fill(HIST("hDCAzDe"), track.dcaZ());
      } else if (track.mcParticle().pdgCode() == 1000010030) {
        histos.fill(HIST("hDCAxyTr"), track.dcaXY());
        histos.fill(HIST("hDCAzTr"), track.dcaZ());
      } else if (track.mcParticle().pdgCode() == 1000020030) {
        histos.fill(HIST("hDCAxyHe"), track.dcaXY());
        histos.fill(HIST("hDCAzHe"), track.dcaZ());
      } else if (track.mcParticle().pdgCode() == 1000020040) {
        histos.fill(HIST("hDCAxyAl"), track.dcaXY());
        histos.fill(HIST("hDCAzAl"), track.dcaZ());
      } 
    }
  }
  PROCESS_SWITCH(lfpropStudy, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lfpropStudy>(cfgc)};
}
