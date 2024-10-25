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

/// \file candidateSelectorXicToXiPiPi.cxx
/// \brief Ξc± → Ξ∓ π± π± candidate selector
///
/// \author Phil Lennart Stahlhut <phil.lennart.stahlhut@cern.ch>, Heidelberg University

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h" // findBin function

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorXicToXiPiPi {
  Produces<aod::HfSelXicToXiPiPi> hfSelXicToXiPiPiCandidate;

  // LF V0 selections
  Configurable<double> radiusV0Min{"radiusV0Min", 1.0, "Min V0 radius"};
  Configurable<double> cosPAV0Min{"cosPAV0Min", 0.95, "Min valueCosPA V0"};
  Configurable<double> dcaV0DauMax{"dcaV0DauMax", 2.0, "Max DCA V0 daughters"};
  Configurable<float> dcaV0ToPvMin{"dcaV0ToPvMin", 0.02, "DCA V0 To PV"};
  Configurable<float> dcaNegToPvMin{"dcaNegToPvMin", 0.05, "DCA Neg To PV"};
  Configurable<float> dcaPosToPvMin{"dcaPosToPvMin", 0.05, "DCA Pos To PV"};
  // Configurable<bool> applyTrkSelLf{"applyTrkSelLf", true, "Apply track selection for LF daughters"};

  // limits for Xic pT
  Configurable<float> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<float> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  
  // pT independent kinematic selections
  Configurable<double> etaTrackPionsFromXicMax{"etaTrackPionsFromXicMax", 0.8, "Max absolute value of eta for charm baryon bachelor"};
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 1.0, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptPiFromCascMin{"ptPiFromCascMin", 0.15, "Min pT pi <- casc"};

  // pT dependent topological and kindematic cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_xi_pi_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xic_to_xi_pi_pi::cuts[0], hf_cuts_xic_to_xi_pi_pi::nBinsPt, hf_cuts_xic_to_xi_pi_pi::nCutVars, hf_cuts_xic_to_xi_pi_pi::labelsPt, hf_cuts_xic_to_xi_pi_pi::labelsCutVar}, "Xicplus candidate selection per pT bin"};
  
  // PID options
  Configurable<bool> usePid{"usePid", true, "Switch for PID selection at track level"};
  Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept Status::NotApplicable [(NotApplicable for one detector) and (NotApplicable or Conditional for the other)] in PID selection"};
  
  // PID - TPC selections
  Configurable<float> ptPiPidTpcMin{"ptPiPidTpcMin", 0.15, "Lower bound of pion pT for TPC PID"};
  Configurable<float> ptPiPidTpcMax{"ptPiPidTpcMax", 20., "Upper bound of pion pT for TPC PID"};
  Configurable<float> nSigmaTpcPiMax{"nSigmaTpcPiMax", 5., "Nsigma cut on TPC only for pions"};
  Configurable<float> nSigmaTpcCombinedPiMax{"nSigmaTpcCombinedPiMax", 5., "Nsigma cut on TPC combined with TOF for pions"};

  Configurable<float> ptPrPidTpcMin{"ptPrPidTpcMin", 0.15, "Lower bound of proton pT for TPC PID"};
  Configurable<float> ptPrPidTpcMax{"ptPrPidTpcMax", 20., "Upper bound of proton pT for TPC PID"};
  Configurable<float> nSigmaTpcPrMax{"nSigmaTpcPrMax", 5., "Nsigma cut on TPC only for protons"};
  Configurable<float> nSigmaTpcCombinedPrMax{"nSigmaTpcCombinedPrMax", 5., "Nsigma cut on TPC combined with TOF for protons"};
  
  // PID - TOF selections
  Configurable<float> ptPiPidTofMin{"ptPiPidTofMin", 0.15, "Lower bound of pion pT for TOF PID"};
  Configurable<float> ptPiPidTofMax{"ptPiPidTofMax", 20., "Upper bound of pion pT for TOF PID"};
  Configurable<float> nSigmaTofPiMax{"nSigmaTofPiMax", 5., "Nsigma cut on TOF only for pions"};
  Configurable<float> nSigmaTofCombinedPiMax{"nSigmaTofCombinedPiMax", 5., "Nsigma cut on TOF combined with TPC for pions"};

  Configurable<float> ptPrPidTofMin{"ptPrPidTofMin", 0.15, "Lower bound of proton pT for TOF PID"};
  Configurable<float> ptPrPidTofMax{"ptPrPidTofMax", 20., "Upper bound of proton pT for TOF PID"};
  Configurable<float> nSigmaTofPrMax{"nSigmaTofPrMax", 5., "Nsigma cut on TOF only for protons"};
  Configurable<float> nSigmaTofCombinedPrMax{"nSigmaTofCombinedPrMax", 5., "Nsigma cut on TOF combined with TPC for protons"};

  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};

  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;

  using TracksPidWithSel = soa::Join<aod::TracksWExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TrackSelection>;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    if (usePid) {
      // pion
      selectorPion.setRangePtTpc(ptPiPidTpcMin, ptPiPidTpcMax);
      selectorPion.setRangeNSigmaTpc(-nSigmaTpcPiMax, nSigmaTpcPiMax);
      selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedPiMax, nSigmaTpcCombinedPiMax);
      selectorPion.setRangePtTof(ptPiPidTofMin, ptPiPidTofMax);
      selectorPion.setRangeNSigmaTof(-nSigmaTofPiMax, nSigmaTofPiMax);
      selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedPiMax, nSigmaTofCombinedPiMax);
      // proton
      selectorProton.setRangePtTpc(ptPrPidTpcMin, ptPrPidTpcMax);
      selectorProton.setRangeNSigmaTpc(-nSigmaTpcPrMax, nSigmaTpcPrMax);
      selectorProton.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedPrMax, nSigmaTpcCombinedPrMax);
      selectorProton.setRangePtTof(ptPrPidTofMin, ptPrPidTofMax);
      selectorProton.setRangeNSigmaTof(-nSigmaTofPrMax, nSigmaTofPrMax);
      selectorProton.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedPrMax, nSigmaTofCombinedPrMax);
    }

    if (activateQA) {
      constexpr int kNBinsSelections = 1 + SelectionStep::NSelectionSteps;
      std::string labels[kNBinsSelections];
      labels[0] = "No selection";
      labels[1 + SelectionStep::RecoSkims] = "Skims selection";
      labels[1 + SelectionStep::RecoTopol] = "Skims & Topological selections";
      labels[1 + SelectionStep::RecoPID] = "Skims & Topological & PID selections";
      labels[1 + SelectionStep::RecoMl] = "Skims & Topological & PID & ML selections";
      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
      registry.add("hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH2>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }
  }

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T1>
  bool selectionTopol(const T1& hfCandXic)
  {
    auto candpT = hfCandXic.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false;
    }

    // check candidate mass is within a defined mass window
    if (std::abs(hfCandXic.invMassXic() - o2::constants::physics::MassXiCPlus) > cuts->get(pTBin, "m")) {
      return false;
    }

    // cosine of pointing angle
    if (hfCandXic.cpa() <= cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }

    // cosine of pointing angle XY
    if (hfCandXic.cpaXY() <= cuts->get(pTBin, "cos pointing angle XY")) {
      return false;
    }

    // candidate maximum decay length
    if (hfCandXic.decayLength() > cuts->get(pTBin, "max decay length")) {
      return false;
    }

    // candidate maximum decay length XY
    if (hfCandXic.decayLengthXY() > cuts->get(pTBin, "max decay length XY")) {
      return false;
    }

    // candidate chi2PC
    if (hfCandXic.chi2PCA() > cuts->get(pTBin, "chi2PCA")) {
      return false;
    }

    // maximum DCA of daughters
    if ((std::abs(hfCandXic.impactParameter0()) > cuts->get(pTBin, "max impParXY Xi")) ||
        (std::abs(hfCandXic.impactParameter1()) > cuts->get(pTBin, "max impParXY Pi0")) ||
        (std::abs(hfCandXic.impactParameter2()) > cuts->get(pTBin, "max impParXY Pi1"))) {
      return false;
    }

    // cut on daughter pT
    if (hfCandXic.ptProng0() < cuts->get(pTBin, "pT Xi") ||
        hfCandXic.ptProng1() < cuts->get(pTBin, "pT Pi0") ||
        hfCandXic.ptProng2() < cuts->get(pTBin, "pT Pi1")) {
      return false;
    }

    // cut on sum of pion pT
    if (hfCandXic.ptProng1() + hfCandXic.ptProng2() < cuts->get(pTBin, "pT Pi0+Pi1")) {
      return false;
    }

    return true;
  }

  /// Apply PID selection
  /// \param pidTrackPi0   PID status of trackPi0 (prong1 of Xic candidate)
  /// \param pidTrackPi1   PID status of trackPi1 (prong2 of Xic candidate)
  /// \param pidTrackPr    PID status of trackPr (positive daughter of V0 candidate)
  /// \param pidTrackPiLam PID status of trackPiLam (negative daughter of V0 candidate)
  /// \param pidTrackPiXi  PID status of trackPiXi (Bachelor of cascade candidate)
  /// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  /// \return true if prongs of Xic candidate pass all selections
  bool selectionPid(TrackSelectorPID::Status const pidTrackPi0,
                    TrackSelectorPID::Status const pidTrackPi1,
                    TrackSelectorPID::Status const pidTrackPr,
                    TrackSelectorPID::Status const pidTrackPiLam,
                    TrackSelectorPID::Status const pidTrackPiXi,
                    bool const acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && (pidTrackPi0 != TrackSelectorPID::Accepted || pidTrackPi1 != TrackSelectorPID::Accepted || pidTrackPr != TrackSelectorPID::Accepted || pidTrackPiLam != TrackSelectorPID::Accepted || pidTrackPiXi != TrackSelectorPID::Accepted)) {
      return false;
    }
    if (acceptPIDNotApplicable && (pidTrackPi0 == TrackSelectorPID::Rejected || pidTrackPi1 == TrackSelectorPID::Rejected || pidTrackPr == TrackSelectorPID::Rejected || pidTrackPiLam == TrackSelectorPID::Rejected || pidTrackPiXi == TrackSelectorPID::Rejected)) {
      return false;
    }
    return true;
  }

  void process(aod::HfCandXic const& hfCandsXic,
               TracksPidWithSel const&)
  {
    for (const auto& hfCandXic : hfCandsXic) {
      int statusXicToXiPiPi = 0;
      auto ptCandXic = hfCandXic.pt();

      if (activateQA) {
        registry.fill(HIST("hSelections"), 1, ptCandXic);
      }

      // No hfFlag -> by default skim selected
      SETBIT(statusXicToXiPiPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusXicToXiPiPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandXic);
      }

      // topological cuts
      if (!selectionTopol(hfCandXic)) {
        hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
        continue;
      }
      SETBIT(statusXicToXiPiPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusXicToXiPiPi = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandXic);
      }

      // track-level PID selection
      if (usePid) {
        auto trackPi0 = hfCandXic.pi0_as<TracksPidWithSel>();
        auto trackPi1 = hfCandXic.pi1_as<TracksPidWithSel>();
        auto trackV0PosDau = hfCandXic.posTrack_as<TracksPidWithSel>();
        auto trackV0NegDau = hfCandXic.negTrack_as<TracksPidWithSel>();
        auto trackPiFromXi = hfCandXic.bachelor_as<TracksPidWithSel>();
        // assign proton and pion hypothesis to V0 daughters
        auto trackPr = trackV0PosDau;
        auto trackPiFromLam = trackV0NegDau;
        if (hfCandXic.sign() < 0) {
          trackPr = trackV0NegDau;
          trackPiFromLam = trackV0PosDau;
        }
        // PID info
        TrackSelectorPID::Status pidTrackPi0 = selectorPion.statusTpcAndTof(trackPi0);
        TrackSelectorPID::Status pidTrackPi1 = selectorPion.statusTpcAndTof(trackPi1);
        TrackSelectorPID::Status pidTrackPr = selectorProton.statusTpcAndTof(trackPr);
        TrackSelectorPID::Status pidTrackPiLam = selectorPion.statusTpcAndTof(trackPiFromLam);
        TrackSelectorPID::Status pidTrackPiXi = selectorPion.statusTpcAndTof(trackPiFromXi);

        if (!selectionPid(pidTrackPi0, pidTrackPi1, pidTrackPr, pidTrackPiLam, pidTrackPiXi, acceptPIDNotApplicable.value)) {
          hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
          continue;
        }
        SETBIT(statusXicToXiPiPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusXicToXiPiPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandXic);
        }
      }

      hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorXicToXiPiPi>(cfgc)};
}
