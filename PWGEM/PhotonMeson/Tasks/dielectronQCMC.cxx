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
// ========================
//
// This code runs loop over dalitz ee table for dalitz QC.
//    Please write to: daiki.sekihata@cern.ch

#include "TString.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "CCDB/BasicCCDBManager.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMEventCut.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::mcutil;
using namespace o2::aod::pwgem::photon;
using namespace o2::aod::pwgem::dilepton::mcutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyMCTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBit, aod::EMPrimaryElectronMCLabels>;
using MyMCTrack = MyMCTracks::iterator;

struct dielectronQCMC {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};

  ConfigurableAxis ConfMeeBins{"ConfMeeBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00}, "mee bins for output histograms"};
  ConfigurableAxis ConfPteeBins{"ConfPteeBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00}, "pTee bins for output histograms"};
  ConfigurableAxis ConfDCAeeBins{"ConfDCAeeBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCAee bins for output histograms"};

  EMEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", -1, "min. occupancy"};
    Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 1000000000, "max. occupancy"};
  } eventcuts;

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 1e+10, "max pair dca3d in sigma"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", true, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.9, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 70, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_min_TPCNsigmaMu{"cfg_min_TPCNsigmaMu", -0.0, "min. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaMu{"cfg_max_TPCNsigmaMu", +0.0, "max. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -1e+10, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3.0, "min. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3.0, "max. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3.0, "min. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3.0, "max. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};

    // CCDB configuration for PID ML
    Configurable<std::string> BDTLocalPathGamma{"BDTLocalPathGamma", "pid_ml_xgboost.onnx", "Path to the local .onnx file"};
    Configurable<std::string> BDTPathCCDB{"BDTPathCCDB", "Users/d/dsekihat/pwgem/pidml/", "Path on CCDB"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
  } dileptoncuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  struct : ConfigurableGroup {
    std::string prefix = "mctrackcut_group";
    Configurable<float> min_mcPt{"min_mcPt", 0.05, "min. MC pT"};
    Configurable<float> max_mcPt{"max_mcPt", 1e+10, "max. MC pT"};
    Configurable<float> max_mcEta{"max_mcEta", 0.9, "max. MC eta"};
  } mctrackcuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view ele_source_types[9] = {"lf/", "Photon/", "PromptJPsi/", "NonPromptJPsi/", "PromptPsi2S/", "NonPromptPsi2S/", "c2e/", "b2e/", "b2c2e/"};

  ~dielectronQCMC() {}

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry, cfgDoFlow);

    const AxisSpec axis_mass{ConfMeeBins, "m_{ee} (GeV/c^{2})"};
    const AxisSpec axis_pt{ConfPteeBins, "p_{T,ee} (GeV/c)"};
    const AxisSpec axis_dca{ConfDCAeeBins, "DCA_{ee}^{3D} (#sigma)"};

    // generated info
    fRegistry.add("Generated/sm/Pi0/hMvsPt", "m_{ee} vs. p_{T,ee} ULS", kTH2F, {axis_mass, axis_pt}, true);
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Eta/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/EtaPrime/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Rho/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Omega/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Phi/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/PromptJPsi/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/NonPromptJPsi/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/PromptPsi2S/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/NonPromptPsi2S/");

    fRegistry.add("Generated/ccbar/c2e_c2e/hadron_hadron/hMvsPt", "m_{ee} vs. p_{T,ee}", kTH2F, {axis_mass, axis_pt}, true);
    fRegistry.addClone("Generated/ccbar/c2e_c2e/hadron_hadron/", "Generated/ccbar/c2e_c2e/meson_meson/");
    fRegistry.addClone("Generated/ccbar/c2e_c2e/hadron_hadron/", "Generated/ccbar/c2e_c2e/baryon_baryon/");
    fRegistry.addClone("Generated/ccbar/c2e_c2e/hadron_hadron/", "Generated/ccbar/c2e_c2e/meson_baryon/");
    fRegistry.addClone("Generated/ccbar/c2e_c2e/", "Generated/bbbar/b2e_b2e/");
    fRegistry.addClone("Generated/ccbar/c2e_c2e/", "Generated/bbbar/b2c2e_b2c2e/");
    fRegistry.addClone("Generated/ccbar/c2e_c2e/", "Generated/bbbar/b2c2e_b2e_sameb/");
    fRegistry.addClone("Generated/ccbar/c2e_c2e/", "Generated/bbbar/b2c2e_b2e_diffb/"); // LS

    // reconstructed pair info
    fRegistry.add("Pair/sm/Photon/hs", "hs pair", kTHnSparseF, {axis_mass, axis_pt, axis_dca}, true);
    fRegistry.add("Pair/sm/Photon/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, false);
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Pi0/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Eta/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/EtaPrime/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Rho/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Omega/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Phi/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/PromptJPsi/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/NonPromptJPsi/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/PromptPsi2S/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/NonPromptPsi2S/");

    fRegistry.add("Pair/ccbar/c2e_c2e/hadron_hadron/hs", "hs pair", kTHnSparseF, {axis_mass, axis_pt, axis_dca}, true);
    fRegistry.addClone("Pair/ccbar/c2e_c2e/hadron_hadron/", "Pair/ccbar/c2e_c2e/meson_meson/");
    fRegistry.addClone("Pair/ccbar/c2e_c2e/hadron_hadron/", "Pair/ccbar/c2e_c2e/baryon_baryon/");
    fRegistry.addClone("Pair/ccbar/c2e_c2e/hadron_hadron/", "Pair/ccbar/c2e_c2e/meson_baryon/");
    fRegistry.addClone("Pair/ccbar/c2e_c2e/", "Pair/bbbar/b2e_b2e/");
    fRegistry.addClone("Pair/ccbar/c2e_c2e/", "Pair/bbbar/b2c2e_b2c2e/");
    fRegistry.addClone("Pair/ccbar/c2e_c2e/", "Pair/bbbar/b2c2e_b2e_sameb/");
    fRegistry.addClone("Pair/ccbar/c2e_c2e/", "Pair/bbbar/b2c2e_b2e_diffb/"); // LS

    // track info
    fRegistry.add("Track/lf/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/lf/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/lf/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, false);
    fRegistry.add("Track/lf/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/lf/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/lf/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {100, 0., 1000}}, false);
    fRegistry.add("Track/lf/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {100, 0., 1000}}, false);
    fRegistry.add("Track/lf/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/lf/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/lf/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/lf/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFbeta", "TOF beta;p_{in} (GeV/c);TOF #beta", kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaMu", "TOF n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/lf/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/lf/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/lf/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/lf/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Track/lf/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH1F, {{32, 0, 16}}, false);
    fRegistry.addClone("Track/lf/", "Track/Photon/");
    fRegistry.addClone("Track/lf/", "Track/PromptJPsi/");
    fRegistry.addClone("Track/lf/", "Track/NonPromptJPsi/");
    fRegistry.addClone("Track/lf/", "Track/PromptPsi2S/");
    fRegistry.addClone("Track/lf/", "Track/NonPromptPsi2S/");
    fRegistry.addClone("Track/lf/", "Track/c2e/");
    fRegistry.addClone("Track/lf/", "Track/b2e/");
    fRegistry.addClone("Track/lf/", "Track/b2c2e/");
  }

  bool cfgDoFlow = false;
  void init(InitContext&)
  {
    DefineEMEventCut();
    DefineDileptonCut();
    addhistograms();
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetOccupancyRange(eventcuts.cfgOccupancyMin, eventcuts.cfgOccupancyMax);
  }

  void DefineDileptonCut()
  {
    fDileptonCut = DalitzEECut("fDileptonCut", "fDileptonCut");

    // for pair
    fDileptonCut.SetMeeRange(dileptoncuts.cfg_min_mass, dileptoncuts.cfg_max_mass);
    fDileptonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dileptoncuts.cfg_phiv_intercept) / dileptoncuts.cfg_phiv_slope; });
    fDileptonCut.SetPairDCARange(dileptoncuts.cfg_min_pair_dca3d, dileptoncuts.cfg_max_pair_dca3d); // in sigma
    fDileptonCut.ApplyPhiV(dileptoncuts.cfg_apply_phiv);
    fDileptonCut.ApplyPrefilter(dileptoncuts.cfg_apply_pf);
    fDileptonCut.RequireITSibAny(dileptoncuts.cfg_require_itsib_any);
    fDileptonCut.RequireITSib1st(dileptoncuts.cfg_require_itsib_1st);

    // for track
    fDileptonCut.SetTrackPtRange(dileptoncuts.cfg_min_pt_track, 1e+10f);
    fDileptonCut.SetTrackEtaRange(-dileptoncuts.cfg_max_eta_track, +dileptoncuts.cfg_max_eta_track);
    fDileptonCut.SetMinNClustersTPC(dileptoncuts.cfg_min_ncluster_tpc);
    fDileptonCut.SetMinNCrossedRowsTPC(dileptoncuts.cfg_min_ncrossedrows);
    fDileptonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDileptonCut.SetChi2PerClusterTPC(0.0, dileptoncuts.cfg_max_chi2tpc);
    fDileptonCut.SetChi2PerClusterITS(0.0, dileptoncuts.cfg_max_chi2its);
    fDileptonCut.SetNClustersITS(dileptoncuts.cfg_min_ncluster_its, 7);
    fDileptonCut.SetMeanClusterSizeITSob(0, 16);
    fDileptonCut.SetMaxDcaXY(dileptoncuts.cfg_max_dcaxy);
    fDileptonCut.SetMaxDcaZ(dileptoncuts.cfg_max_dcaz);

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaMuRange(dileptoncuts.cfg_min_TPCNsigmaMu, dileptoncuts.cfg_max_TPCNsigmaMu);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTPCNsigmaKaRange(dileptoncuts.cfg_min_TPCNsigmaKa, dileptoncuts.cfg_max_TPCNsigmaKa);
    fDileptonCut.SetTPCNsigmaPrRange(dileptoncuts.cfg_min_TPCNsigmaPr, dileptoncuts.cfg_max_TPCNsigmaPr);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);

    if (dileptoncuts.cfg_pid_scheme == static_cast<int>(DalitzEECut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      o2::ml::OnnxModel* eid_bdt = new o2::ml::OnnxModel();
      if (dileptoncuts.loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        std::map<std::string, std::string> metadata;
        bool retrieveSuccessGamma = ccdbApi.retrieveBlob(dileptoncuts.BDTPathCCDB.value, ".", metadata, dileptoncuts.timestampCCDB.value, false, dileptoncuts.BDTLocalPathGamma.value);
        if (retrieveSuccessGamma) {
          eid_bdt->initModel(dileptoncuts.BDTLocalPathGamma.value, dileptoncuts.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the Gamma model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      } else {
        eid_bdt->initModel(dileptoncuts.BDTLocalPathGamma.value, dileptoncuts.enableOptimizations.value);
      }

      fDileptonCut.SetPIDModel(eid_bdt);
    } // end of PID ML
  }

  template <typename TTrack, typename TMCParticles>
  int FindLF(TTrack const& posmc, TTrack const& elemc, TMCParticles const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 111, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 221, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 331, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 113, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 223, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 333, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 443, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 100443, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <typename T>
  bool isInAcceptance(T const& t1)
  {
    if ((mctrackcuts.min_mcPt < t1.pt() && t1.pt() < mctrackcuts.max_mcPt) && abs(t1.eta()) < mctrackcuts.max_mcEta) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TCollision, typename TTrack1, typename TTrack2, typename TMCParticles>
  bool fillTruePairInfo(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, TMCParticles const& mcparticles)
  {
    if (dileptoncuts.cfg_pid_scheme == static_cast<int>(DalitzEECut::PIDSchemes::kPIDML)) {
      if (!fDileptonCut.IsSelectedTrack<true>(t1, collision) || !fDileptonCut.IsSelectedTrack<true>(t2, collision)) {
        return false;
      }
    } else { // cut-based
      if (!fDileptonCut.IsSelectedTrack(t1) || !fDileptonCut.IsSelectedTrack(t2)) {
        return false;
      }
    }

    if (!fDileptonCut.IsSelectedPair(t1, t2, collision.bz())) {
      return false;
    }

    auto t1mc = t1.template emmcparticle_as<TMCParticles>();
    auto t2mc = t2.template emmcparticle_as<TMCParticles>();

    int mother_id = FindLF(t1mc, t2mc, mcparticles);
    int hfee_type = IsHF(t1mc, t2mc, mcparticles);
    if (mother_id < 0 && hfee_type < 0) {
      return false;
    }
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    if (abs(v12.Rapidity()) > maxY) {
      return false;
    }

    float dca_t1_3d = t1.dca3DinSigma();
    float dca_t2_3d = t2.dca3DinSigma();
    float dca_ee_3d = std::sqrt((dca_t1_3d * dca_t1_3d + dca_t2_3d * dca_t2_3d) / 2.);
    float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), collision.bz());

    if (mother_id > -1 && t1mc.pdgCode() * t2mc.pdgCode() < 0) {
      auto mcmother = mcparticles.iteratorAt(mother_id);
      if (mcmother.isPhysicalPrimary() || mcmother.producedByGenerator()) {
        if ((t1mc.isPhysicalPrimary() || t1mc.producedByGenerator()) && (t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
          switch (abs(mcmother.pdgCode())) {
            case 111:
              fRegistry.fill(HIST("Pair/sm/Pi0/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              fRegistry.fill(HIST("Pair/sm/Pi0/hMvsPhiV"), phiv, v12.M());
              fillTrackInfo<0>(t1);
              fillTrackInfo<0>(t2);
              break;
            case 221:
              fRegistry.fill(HIST("Pair/sm/Eta/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              fRegistry.fill(HIST("Pair/sm/Eta/hMvsPhiV"), phiv, v12.M());
              fillTrackInfo<0>(t1);
              fillTrackInfo<0>(t2);
              break;
            case 331:
              fRegistry.fill(HIST("Pair/sm/EtaPrime/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              fRegistry.fill(HIST("Pair/sm/EtaPrime/hMvsPhiV"), phiv, v12.M());
              fillTrackInfo<0>(t1);
              fillTrackInfo<0>(t2);
              break;
            case 113:
              fRegistry.fill(HIST("Pair/sm/Rho/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              fRegistry.fill(HIST("Pair/sm/Rho/hMvsPhiV"), phiv, v12.M());
              fillTrackInfo<0>(t1);
              fillTrackInfo<0>(t2);
              break;
            case 223:
              fRegistry.fill(HIST("Pair/sm/Omega/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              fRegistry.fill(HIST("Pair/sm/Omega/hMvsPhiV"), phiv, v12.M());
              fillTrackInfo<0>(t1);
              fillTrackInfo<0>(t2);
              break;
            case 333:
              fRegistry.fill(HIST("Pair/sm/Phi/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              fRegistry.fill(HIST("Pair/sm/Phi/hMvsPhiV"), phiv, v12.M());
              fillTrackInfo<0>(t1);
              fillTrackInfo<0>(t2);
              break;
            case 443: {
              if (IsFromBeauty(mcmother, mcparticles) > 0) {
                fRegistry.fill(HIST("Pair/sm/NonPromptJPsi/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fRegistry.fill(HIST("Pair/sm/NonPromptJPsi/hMvsPhiV"), phiv, v12.M());
                fillTrackInfo<3>(t1);
                fillTrackInfo<3>(t2);
              } else {
                fRegistry.fill(HIST("Pair/sm/PromptJPsi/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fRegistry.fill(HIST("Pair/sm/PromptJPsi/hMvsPhiV"), phiv, v12.M());
                fillTrackInfo<2>(t1);
                fillTrackInfo<2>(t2);
              }
              break;
            }
            case 100443: {
              if (IsFromBeauty(mcmother, mcparticles) > 0) {
                fRegistry.fill(HIST("Pair/sm/NonPromptPsi2S/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fRegistry.fill(HIST("Pair/sm/NonPromptPsi2S/hMvsPhiV"), phiv, v12.M());
                fillTrackInfo<5>(t1);
                fillTrackInfo<5>(t2);
              } else {
                fRegistry.fill(HIST("Pair/sm/PromptPsi2S/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fRegistry.fill(HIST("Pair/sm/PromptPsi2S/hMvsPhiV"), phiv, v12.M());
                fillTrackInfo<4>(t1);
                fillTrackInfo<4>(t2);
              }
              break;
            }
            default:
              break;
          }

        } else if (!(t1mc.isPhysicalPrimary() || t1mc.producedByGenerator()) && !(t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
          switch (abs(mcmother.pdgCode())) {
            case 22:
              fRegistry.fill(HIST("Pair/sm/Photon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              fRegistry.fill(HIST("Pair/sm/Photon/hMvsPhiV"), phiv, v12.M());
              fillTrackInfo<1>(t1);
              fillTrackInfo<1>(t2);
              break;
            default:
              break;
          }
        } // end of primary/secondary selection
      }   // end of primary selection for same mother
    } else if (hfee_type > -1) {
      if ((t1mc.isPhysicalPrimary() || t1mc.producedByGenerator()) && (t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
        auto mp1 = mcparticles.iteratorAt(t1mc.mothersIds()[0]);
        auto mp2 = mcparticles.iteratorAt(t2mc.mothersIds()[0]);
        if (t1mc.pdgCode() * t2mc.pdgCode() < 0) { // ULS
          switch (hfee_type) {
            case static_cast<int>(EM_HFeeType::kCe_Ce): {
              fRegistry.fill(HIST("Pair/ccbar/c2e_c2e/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
                fRegistry.fill(HIST("Pair/ccbar/c2e_c2e/meson_meson/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fillTrackInfo<6>(t1);
                fillTrackInfo<6>(t2);
              } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
                fRegistry.fill(HIST("Pair/ccbar/c2e_c2e/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fillTrackInfo<6>(t1);
                fillTrackInfo<6>(t2);
              } else {
                fRegistry.fill(HIST("Pair/ccbar/c2e_c2e/meson_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fillTrackInfo<6>(t1);
                fillTrackInfo<6>(t2);
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBe_Be): {
              fRegistry.fill(HIST("Pair/bbbar/b2e_b2e/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              if (isBeautyMeson(mp1) && isBeautyMeson(mp2)) {
                fRegistry.fill(HIST("Pair/bbbar/b2e_b2e/meson_meson/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fillTrackInfo<7>(t1);
                fillTrackInfo<7>(t2);
              } else if (isBeautyBaryon(mp1) && isBeautyBaryon(mp2)) {
                fRegistry.fill(HIST("Pair/bbbar/b2e_b2e/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fillTrackInfo<7>(t1);
                fillTrackInfo<7>(t2);
              } else {
                fRegistry.fill(HIST("Pair/bbbar/b2e_b2e/meson_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fillTrackInfo<7>(t1);
                fillTrackInfo<7>(t2);
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_BCe): {
              fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2c2e/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2c2e/meson_meson/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fillTrackInfo<8>(t1);
                fillTrackInfo<8>(t2);
              } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2c2e/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fillTrackInfo<8>(t1);
                fillTrackInfo<8>(t2);
              } else {
                fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2c2e/meson_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
                fillTrackInfo<8>(t1);
                fillTrackInfo<8>(t2);
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): { // ULS
              fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2e_sameb/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2e_sameb/meson_meson/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2e_sameb/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              } else {
                fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2e_sameb/meson_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              }
              if ((isCharmMeson(mp1) || isCharmBaryon(mp1)) && (isBeautyMeson(mp2) || isBeautyBaryon(mp2))) {
                fillTrackInfo<7>(t1);
                fillTrackInfo<8>(t2);
              } else {
                fillTrackInfo<8>(t1);
                fillTrackInfo<7>(t2);
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): // LS
              LOGF(info, "You should not see kBCe_Be_DiffB in ULS. Good luck.");
              break;
            default:
              break;
          }
        } else { // LS
          switch (hfee_type) {
            case static_cast<int>(EM_HFeeType::kCe_Ce):
              LOGF(info, "You should not see kCe_Ce in LS. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBe_Be):
              LOGF(info, "You should not see kBe_Be in LS. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBCe_BCe):
              LOGF(info, "You should not see kBCe_BCe in LS. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): // ULS
              LOGF(info, "You should not see kBCe_Be_SameB in LS. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): { // LS
              fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2e_diffb/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2e_diffb/meson_meson/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2e_diffb/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              } else {
                fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2e_diffb/meson_baryon/hs"), v12.M(), v12.Pt(), dca_ee_3d);
              }
              if ((isCharmMeson(mp1) || isCharmBaryon(mp1)) && (isBeautyMeson(mp2) || isBeautyBaryon(mp2))) {
                fillTrackInfo<7>(t1);
                fillTrackInfo<8>(t2);
              } else {
                fillTrackInfo<8>(t1);
                fillTrackInfo<7>(t2);
              }
              break;
            }
            default:
              break;
          }
        }
      }
    } // end of HF evaluation
    return true;
  }

  template <int e_source_id, typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    // fill track info that belong to true pairs.
    if (std::find(used_trackIds.begin(), used_trackIds.end(), track.globalIndex()) == used_trackIds.end()) {
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hPt"), track.pt());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hEtaPhi"), track.phi(), track.eta());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hMeanClusterSizeITS"), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFbeta"), track.tpcInnerParam(), track.beta());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
      used_trackIds.emplace_back(track.globalIndex());
    }
  }

  std::vector<int> used_trackIds;
  SliceCache cache;
  Preslice<MyMCTracks> perCollision_track = aod::emprimaryelectron::emeventId;
  Filter trackFilter = static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && o2::aod::track::tpcChi2NCl < static_cast<float>(dileptoncuts.cfg_max_chi2tpc) && o2::aod::track::itsChi2NCl < static_cast<float>(dileptoncuts.cfg_max_chi2its) && nabs(o2::aod::track::dcaXY) < static_cast<float>(dileptoncuts.cfg_max_dcaxy) && nabs(o2::aod::track::dcaZ) < static_cast<float>(dileptoncuts.cfg_max_dcaz);
  Filter pidFilter = (static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl)) && (o2::aod::pidtpc::tpcNSigmaPi < static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaPi) || static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaPi) < o2::aod::pidtpc::tpcNSigmaPi) && ((0.96f < o2::aod::pidtofbeta::beta && o2::aod::pidtofbeta::beta < 1.04f) || o2::aod::pidtofbeta::beta < 0.f);
  using FilteredMyMCTracks = soa::Filtered<MyMCTracks>;
  Partition<FilteredMyMCTracks> posTracks = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyMCTracks> negTracks = o2::aod::emprimaryelectron::sign < int8_t(0);

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processQCMC(FilteredMyCollisions const& collisions, FilteredMyMCTracks const& tracks, aod::EMMCParticles const& mcparticles, aod::EMMCEvents const&)
  {
    used_trackIds.reserve(tracks.size());

    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, cfgDoFlow);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, cfgDoFlow);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      // LOGF(info, "centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        fillTruePairInfo(collision, pos, ele, mcparticles);
      } // end of ULS pair loop

      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        fillTruePairInfo(collision, pos1, pos2, mcparticles);
      } // end of ULS pair loop

      for (auto& [ele1, ele2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS__
        fillTruePairInfo(collision, ele1, ele2, mcparticles);
      } // end of ULS pair loop

    } // end of collision loop

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();

  } // end of process
  PROCESS_SWITCH(dielectronQCMC, processQCMC, "run Dalitz QC", true);

  Partition<aod::EMMCParticles> posTracksMC = o2::aod::mcparticle::pdgCode == -11; // e+
  Partition<aod::EMMCParticles> negTracksMC = o2::aod::mcparticle::pdgCode == +11; // e-
  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  void processGen(MyCollisions const& collisions, aod::EMMCEvents const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      auto mccollision = collision.emmcevent_as<aod::EMMCEvents>();

      auto posTracks_per_coll = posTracksMC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);
      auto negTracks_per_coll = negTracksMC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);

      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        // LOGF(info, "pdg1 = %d, pdg2 = %d", t1.pdgCode(), t2.pdgCode());

        if (!isInAcceptance(t1) || !isInAcceptance(t2)) {
          continue;
        }

        if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
          continue;
        }
        if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
          continue;
        }

        int mother_id = FindLF(t1, t2, mcparticles);
        int hfee_type = IsHF(t1, t2, mcparticles);
        if (mother_id < 0 && hfee_type < 0) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        if (abs(v12.Rapidity()) > maxY) {
          continue;
        }

        if (mother_id > -1) {
          auto mcmother = mcparticles.iteratorAt(mother_id);
          if (mcmother.isPhysicalPrimary() || mcmother.producedByGenerator()) {

            switch (abs(mcmother.pdgCode())) {
              case 111:
                fRegistry.fill(HIST("Generated/sm/Pi0/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 221:
                fRegistry.fill(HIST("Generated/sm/Eta/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 331:
                fRegistry.fill(HIST("Generated/sm/EtaPrime/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 113:
                fRegistry.fill(HIST("Generated/sm/Rho/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 223:
                fRegistry.fill(HIST("Generated/sm/Omega/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 333:
                fRegistry.fill(HIST("Generated/sm/Phi/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 443: {
                if (IsFromBeauty(mcmother, mcparticles) > 0) {
                  fRegistry.fill(HIST("Generated/sm/NonPromptJPsi/hMvsPt"), v12.M(), v12.Pt());
                } else {
                  fRegistry.fill(HIST("Generated/sm/PromptJPsi/hMvsPt"), v12.M(), v12.Pt());
                }
                break;
              }
              case 100443: {
                if (IsFromBeauty(mcmother, mcparticles) > 0) {
                  fRegistry.fill(HIST("Generated/sm/NonPromptPsi2S/hMvsPt"), v12.M(), v12.Pt());
                } else {
                  fRegistry.fill(HIST("Generated/sm/PromptPsi2S/hMvsPt"), v12.M(), v12.Pt());
                }
                break;
              }
              default:
                break;
            }
          }
        } else if (hfee_type > -1) {
          auto mp1 = mcparticles.iteratorAt(t1.mothersIds()[0]);
          auto mp2 = mcparticles.iteratorAt(t2.mothersIds()[0]);
          switch (hfee_type) {
            case static_cast<int>(EM_HFeeType::kCe_Ce): {
              fRegistry.fill(HIST("Generated/ccbar/c2e_c2e/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
                fRegistry.fill(HIST("Generated/ccbar/c2e_c2e/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
                fRegistry.fill(HIST("Generated/ccbar/c2e_c2e/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/ccbar/c2e_c2e/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBe_Be): {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2e/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if (isBeautyMeson(mp1) && isBeautyMeson(mp2)) {
                fRegistry.fill(HIST("Generated/bbbar/b2e_b2e/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if (isBeautyBaryon(mp1) && isBeautyBaryon(mp2)) {
                fRegistry.fill(HIST("Generated/bbbar/b2e_b2e/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2e_b2e/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_BCe): {
              fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2c2e/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
                fRegistry.fill(HIST("Generated/bbbar/b2e_b2e/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
                fRegistry.fill(HIST("Generated/bbbar/b2e_b2e/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2e_b2e/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): { // ULS
              fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_sameb/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_sameb/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_sameb/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_sameb/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): // LS
              LOGF(info, "You should not see kBCe_Be_DiffB in ULS. Good luck.");
              break;
            default:
              break;
          }
        } // end of HF evaluation
      }   // end of true ULS pair loop

      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) {
        // LOGF(info, "pdg1 = %d, pdg2 = %d", t1.pdgCode(), t2.pdgCode());

        if (!isInAcceptance(t1) || !isInAcceptance(t2)) {
          continue;
        }

        if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
          continue;
        }
        if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
          continue;
        }

        int hfee_type = IsHF(t1, t2, mcparticles);
        if (hfee_type < 0) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (abs(v12.Rapidity()) > maxY) {
          continue;
        }
        if (hfee_type > -1) {
          auto mp1 = mcparticles.iteratorAt(t1.mothersIds()[0]);
          auto mp2 = mcparticles.iteratorAt(t2.mothersIds()[0]);
          switch (hfee_type) {
            case static_cast<int>(EM_HFeeType::kCe_Ce):
              LOGF(info, "You should not see kCe_Ce in LS++. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBe_Be):
              LOGF(info, "You should not see kBe_Be in LS++. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBCe_BCe):
              LOGF(info, "You should not see kBCe_BCe in LS++. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): // ULS
              LOGF(info, "You should not see kBCe_Be_SameB in LS++. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): { // LS
              fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_diffb/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_diffb/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_diffb/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_diffb/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            default:
              break;
          }
        }
      } // end of true LS++ pair loop

      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) {
        // LOGF(info, "pdg1 = %d, pdg2 = %d", t1.pdgCode(), t2.pdgCode());

        if (!isInAcceptance(t1) || !isInAcceptance(t2)) {
          continue;
        }

        if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
          continue;
        }
        if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
          continue;
        }

        int hfee_type = IsHF(t1, t2, mcparticles);
        if (hfee_type < 0) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (abs(v12.Rapidity()) > maxY) {
          continue;
        }
        if (hfee_type > -1) {
          auto mp1 = mcparticles.iteratorAt(t1.mothersIds()[0]);
          auto mp2 = mcparticles.iteratorAt(t2.mothersIds()[0]);
          switch (hfee_type) {
            case static_cast<int>(EM_HFeeType::kCe_Ce):
              LOGF(info, "You should not see kCe_Ce in LS--. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBe_Be):
              LOGF(info, "You should not see kBe_Be in LS--. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBCe_BCe):
              LOGF(info, "You should not see kBCe_BCe in LS--. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): // ULS
              LOGF(info, "You should not see kBCe_Be_SameB in LS--. Good luck.");
              break;
            case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): { // LS
              fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_diffb/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_diffb/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_diffb/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2e_diffb/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            default:
              break;
          }
        }
      } // end of true LS++ pair loop

    } // end of collision loop
  }
  PROCESS_SWITCH(dielectronQCMC, processGen, "run genrated info", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(dielectronQCMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dielectronQCMC>(cfgc, TaskName{"dielectron-qc-mc"})};
}
