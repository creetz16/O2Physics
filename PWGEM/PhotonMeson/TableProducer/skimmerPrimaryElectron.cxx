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

/// \brief write relevant information for dalitz ee analysis to an AO2D.root file. This file is then the only necessary input to perform pcm analysis.
/// \author daiki.sekihata@cern.ch

#include <map>
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::pwgem::photonmeson;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTrack = MyTracks::iterator;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;
using MyTrackMC = MyTracksMC::iterator;

struct skimmerPrimaryElectron {
  enum class EM_EEPairType : int {
    kULS = 0,
    kLSpp = +1,
    kLSnn = -1,
  };

  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;
  Produces<aod::EMPrimaryElectrons> emprimaryelectrons;
  Produces<aod::EMPrimaryElectronsCov> emprimaryelectronscov;
  // Produces<o2::aod::EMEventsBz> events_bz;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // Operation and minimisation criteria
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 10, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<float> max_mean_itsob_cluster_size{"max_mean_itsob_cluster_size", 16.f, "max. <ITSob cluster size> x cos(lambda)"}; // this is to suppress random combination. default 4 + 1 for skimming.
  Configurable<int> minitsncls{"minitsncls", 4, "min. number of ITS clusters"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max. chi2/NclsITS"};
  Configurable<float> minpt{"minpt", 0.05, "min pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1.0f, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 1.0f, "max DCAz in cm"};
  Configurable<float> dca_3d_sigma_max{"dca_3d_sigma_max", 1e+10, "max DCA 3D in sigma"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -3.0, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 4.0, "max. TPC n sigma for electron inclusion"};
  Configurable<float> maxTOFNsigmaEl{"maxTOFNsigmaEl", 4.0, "max. TOF n sigma for electron inclusion"};
  Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", 0.0, "min. TPC n sigma for pion exclusion"}; // set to -2 for lowB, -999 for nominalB
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 0.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> maxTPCNsigmaKa{"maxTPCNsigmaKa", 2.0, "max. TPC n sigma for kaon exclusion"};
  Configurable<float> maxTPCNsigmaPr{"maxTPCNsigmaPr", 2.0, "max. TPC n sigma for proton exclusion"};
  Configurable<float> maxMee{"maxMee", 1e+10, "max. mee to store ee pairs"};
  Configurable<bool> storeLS{"storeLS", false, "flag to store LS pairs"};
  Configurable<bool> applyKaRej_TPC{"applyKaRej_TPC", false, "flag to apply Kaon rejection in TPC at the skimming level"}; // used for CEFP only
  Configurable<bool> applyPrRej_TPC{"applyPrRej_TPC", false, "flag to apply Kaon rejection in TPC at the skimming level"}; // used for CEFP only

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hNpairs", "hNpairs;pair type;Number of Pairs", {HistType::kTH1F, {{3, -1.5f, +1.5f}}}},
      {"Track/hTPCdEdx_Pin_before", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0.f, 200.f}}}},
      {"Track/hTPCdEdx_Pin_after", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0.f, 200.f}}}},
      {"Track/hTOFbeta_Pin_before", "TOF beta vs. p_{in};p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {240, 0.f, 1.2f}}}},
      {"Track/hTOFbeta_Pin_after", "TOF beta vs. p_{in};p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {240, 0.f, 1.2f}}}},
      {"Track/hTPCNsigmaEl_before", "TPC n sigma e vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaEl_before", "TOF n sigma e vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaMu_before", "TPC n sigma #mu vs. p_{in};p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaMu_before", "TOF n sigma #mu vs. p_{in};p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaPi_before", "TPC n sigma #pi vs. p_{in};p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaPi_before", "TOF n sigma #pi vs. p_{in};p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaEl_after", "TPC n sigma e vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaEl_after", "TOF n sigma e vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaMu_after", "TPC n sigma #mu vs. p_{in};p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaMu_after", "TOF n sigma #mu vs. p_{in};p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaPi_after", "TPC n sigma #pi vs. p_{in};p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaPi_after", "TOF n sigma #pi vs. p_{in};p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"Pair/hMeePtee_ULS", "ULS m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", {HistType::kTH2F, {{400, 0.f, 4.f}, {100, 0.f, 10.f}}}},
      {"Pair/hMeePtee_LSpp", "LS++ m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", {HistType::kTH2F, {{400, 0.f, 4.f}, {100, 0.f, 10.f}}}},
      {"Pair/hMeePtee_LSnn", "LS-- m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", {HistType::kTH2F, {{400, 0.f, 4.f}, {100, 0.f, 10.f}}}},
    },
  };

  std::pair<int8_t, std::set<uint8_t>> itsRequirement = {1, {0, 1, 2}}; // any hits on 3 ITS ib layers.

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  template <bool isMC, typename TCollision, typename TTrack>
  bool checkTrack(TCollision const& collision, TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (track.tpcChi2NCl() > maxchi2tpc) {
      return false;
    }

    if (track.itsChi2NCl() > maxchi2its) {
      return false;
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }
    if (track.itsNCls() < minitsncls) {
      return false;
    }

    auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    if (hits < itsRequirement.first) {
      return false;
    }

    uint32_t itsClusterSizes = track.itsClusterSizes();
    int total_cluster_size = 0, nl = 0;
    for (unsigned int layer = 3; layer < 7; layer++) {
      int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
      if (cluster_size_per_layer > 0) {
        nl++;
      }
      total_cluster_size += cluster_size_per_layer;
    }
    if (static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(track.tgl())) > max_mean_itsob_cluster_size) {
      return false;
    }

    if (track.tpcNClsFound() < min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < min_tpc_cr_findable_ratio) {
      return false;
    }

    if ((0.0 < track.beta() && track.beta() < 0.95) || 1.05 < track.beta()) {
      return false;
    }

    gpu::gpustd::array<float, 2> dcaInfo;
    auto track_par_cov_recalc = getTrackParCov(track);
    std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
    getPxPyPz(track_par_cov_recalc, pVec_recalc);
    float dcaXY = dcaInfo[0];
    float dcaZ = dcaInfo[1];

    if (track_par_cov_recalc.getPt() < minpt || abs(track_par_cov_recalc.getEta()) > maxeta) {
      return false;
    }

    float dca_3d = 999.f;
    float det = track_par_cov_recalc.getSigmaY2() * track_par_cov_recalc.getSigmaZ2() - track_par_cov_recalc.getSigmaZY() * track_par_cov_recalc.getSigmaZY();
    if (det < 0) {
      dca_3d = 999.f;
    } else {
      float chi2 = (dcaXY * dcaXY * track_par_cov_recalc.getSigmaZ2() + dcaZ * dcaZ * track_par_cov_recalc.getSigmaY2() - 2. * dcaXY * dcaZ * track_par_cov_recalc.getSigmaZY()) / det;
      dca_3d = std::sqrt(std::abs(chi2) / 2.);
    }
    if (dca_3d > dca_3d_sigma_max) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    return isElectron_TPChadrej(track) || isElectron_TOFrecovery(track);
  }

  template <typename TTrack>
  bool isElectron_TPChadrej(TTrack const& track)
  {
    if (track.tpcNSigmaEl() < minTPCNsigmaEl || maxTPCNsigmaEl < track.tpcNSigmaEl()) {
      return false;
    }
    if (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi) {
      return false;
    }
    if (applyKaRej_TPC && abs(track.tpcNSigmaKa()) < maxTPCNsigmaKa) {
      return false;
    }
    if (applyPrRej_TPC && abs(track.tpcNSigmaPr()) < maxTPCNsigmaPr) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron_TOFrecovery(TTrack const& track)
  {
    if (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi) {
      return false;
    }
    return minTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < maxTPCNsigmaEl && abs(track.tofNSigmaEl()) < maxTOFNsigmaEl;
  }

  template <bool isMC, typename TCollision, typename TTracks>
  void fillTrackHistogram(TCollision const& collision, TTracks const& tracks)
  {
    for (auto& track : tracks) {
      if (!checkTrack<isMC>(collision, track) || !isElectron(track)) {
        continue;
      }
      fRegistry.fill(HIST("Track/hTPCdEdx_Pin_before"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/hTOFbeta_Pin_before"), track.tpcInnerParam(), track.beta());
      fRegistry.fill(HIST("Track/hTPCNsigmaEl_before"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/hTOFNsigmaEl_before"), track.tpcInnerParam(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/hTPCNsigmaMu_before"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/hTOFNsigmaMu_before"), track.tpcInnerParam(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/hTPCNsigmaPi_before"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/hTOFNsigmaPi_before"), track.tpcInnerParam(), track.tofNSigmaPi());
    }
  }

  template <typename TCollision, typename TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& track)
  {
    if (std::find(stored_trackIds.begin(), stored_trackIds.end(), std::make_pair(collision.globalIndex(), track.globalIndex())) == stored_trackIds.end()) {
      gpu::gpustd::array<float, 2> dcaInfo;
      auto track_par_cov_recalc = getTrackParCov(track);
      std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
      getPxPyPz(track_par_cov_recalc, pVec_recalc);
      float dcaXY = dcaInfo[0];
      float dcaZ = dcaInfo[1];

      float pt_recalc = track_par_cov_recalc.getPt();
      float eta_recalc = track_par_cov_recalc.getEta();
      float phi_recalc = track_par_cov_recalc.getPhi();
      // float tgl_recalc = track_par_cov_recalc.getTgl();
      // float cYY_recalc = track_par_cov_recalc.getSigmaY2();
      // float cZZ_recalc = track_par_cov_recalc.getSigmaZ2();
      // float cZY_recalc = track_par_cov_recalc.getSigmaZY();

      // if (collision.globalIndex() != track.collisionId()) {
      //   LOGF(info, "dca_xy: before = %f , after = %f | pT: before = %f , after = %f | cYY : before = %f , after = %f", track.dcaXY(), dcaInfo[0], track.pt(), pt_recalc, track.cYY(), cYY_recalc);
      // }

      emprimaryelectrons(collision.globalIndex(), track.globalIndex(), track.sign(),
                         pt_recalc, eta_recalc, phi_recalc, dcaXY, dcaZ,
                         track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                         track.tpcChi2NCl(), track.tpcInnerParam(),
                         track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                         track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                         track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap(),
                         track_par_cov_recalc.getX(), track_par_cov_recalc.getAlpha(), track_par_cov_recalc.getY(), track_par_cov_recalc.getZ(), track_par_cov_recalc.getSnp(), track_par_cov_recalc.getTgl());

      emprimaryelectronscov(
        track_par_cov_recalc.getSigmaY2(),
        track_par_cov_recalc.getSigmaZY(),
        track_par_cov_recalc.getSigmaZ2(),
        track_par_cov_recalc.getSigmaSnpY(),
        track_par_cov_recalc.getSigmaSnpZ(),
        track_par_cov_recalc.getSigmaSnp2(),
        track_par_cov_recalc.getSigmaTglY(),
        track_par_cov_recalc.getSigmaTglZ(),
        track_par_cov_recalc.getSigmaTglSnp(),
        track_par_cov_recalc.getSigmaTgl2(),
        track_par_cov_recalc.getSigma1PtY(),
        track_par_cov_recalc.getSigma1PtZ(),
        track_par_cov_recalc.getSigma1PtSnp(),
        track_par_cov_recalc.getSigma1PtTgl(),
        track_par_cov_recalc.getSigma1Pt2());

      fRegistry.fill(HIST("Track/hTPCdEdx_Pin_after"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/hTOFbeta_Pin_after"), track.tpcInnerParam(), track.beta());
      fRegistry.fill(HIST("Track/hTPCNsigmaEl_after"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/hTOFNsigmaEl_after"), track.tpcInnerParam(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/hTPCNsigmaMu_after"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/hTOFNsigmaMu_after"), track.tpcInnerParam(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/hTPCNsigmaPi_after"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/hTOFNsigmaPi_after"), track.tpcInnerParam(), track.tofNSigmaPi());
      stored_trackIds.emplace_back(std::make_pair(collision.globalIndex(), track.globalIndex()));
    }
  }

  template <bool isMC, EM_EEPairType pairtype, typename TCollision, typename TTracks1, typename TTracks2>
  void fillPairInfo(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    if constexpr (pairtype == EM_EEPairType::kULS) { // ULS
      for (auto& t1 : tracks1) {
        for (auto& t2 : tracks2) {
          if (!checkTrack<isMC>(collision, t1) || !checkTrack<isMC>(collision, t2)) {
            continue;
          }
          if (!isElectron(t1) || !isElectron(t2)) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          fRegistry.fill(HIST("Pair/hMeePtee_ULS"), v12.M(), v12.Pt());

          if (v12.M() > maxMee) { // don't store
            continue;
          }
          fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
          fillTrackTable(collision, t1);
          fillTrackTable(collision, t2);
        }    // end of t2
      }      // end of t1
    } else { // LS
      std::vector<std::pair<int64_t, int64_t>> used_pairs;
      used_pairs.reserve(tracks1.size() * tracks2.size());

      for (auto& t1 : tracks1) {
        for (auto& t2 : tracks2) {
          if (t1.globalIndex() == t2.globalIndex()) {
            continue;
          }

          if (std::find(used_pairs.begin(), used_pairs.end(), std::make_pair(t1.globalIndex(), t2.globalIndex())) == used_pairs.end()) {
            if (!checkTrack<isMC>(collision, t1) || !checkTrack<isMC>(collision, t2)) {
              continue;
            }
            if (!isElectron(t1) || !isElectron(t2)) {
              continue;
            }

            ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

            if (pairtype == EM_EEPairType::kLSpp) {
              fRegistry.fill(HIST("Pair/hMeePtee_LSpp"), v12.M(), v12.Pt());
            } else if (pairtype == EM_EEPairType::kLSnn) {
              fRegistry.fill(HIST("Pair/hMeePtee_LSnn"), v12.M(), v12.Pt());
            }

            if (v12.M() > maxMee) { // don't store
              continue;
            }
            fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
            fillTrackTable(collision, t1);
            fillTrackTable(collision, t2);

            used_pairs.emplace_back(std::make_pair(t1.globalIndex(), t2.globalIndex()));
            used_pairs.emplace_back(std::make_pair(t2.globalIndex(), t1.globalIndex()));
          }
        } // end of t2
      }   // end of t1
      used_pairs.clear();
      used_pairs.shrink_to_fit();
    }
  }

  std::vector<std::pair<int64_t, int64_t>> stored_trackIds;
  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& nabs(o2::aod::track::dcaXY) < dca_xy_max&& nabs(o2::aod::track::dcaZ) < dca_z_max&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& o2::aod::track::itsChi2NCl < maxchi2its&& ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  Filter pidFilter = minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl && ((0.96f < o2::aod::pidtofbeta::beta && o2::aod::pidtofbeta::beta < 1.04f) || o2::aod::pidtofbeta::beta < 0.f);
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  // ---------- for data ----------

  void processRec_SA(aod::Collisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    stored_trackIds.reserve(tracks.size());

    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      // events_bz(d_bz);

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      fillTrackHistogram<false>(collision, posTracks_per_coll);
      fillTrackHistogram<false>(collision, negTracks_per_coll);

      fillPairInfo<false, EM_EEPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        fillPairInfo<false, EM_EEPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        fillPairInfo<false, EM_EEPairType::kLSnn>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processRec_SA, "process reconstructed info only", true); // standalone

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  void processRec_TTCA(aod::Collisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks, aod::TrackAssoc const& trackIndices)
  {
    stored_trackIds.reserve(tracks.size() * 2);

    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      // events_bz(d_bz);

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      std::vector<MyTrack> posTracks_per_coll;
      std::vector<MyTrack> negTracks_per_coll;
      posTracks_per_coll.reserve(trackIdsThisCollision.size());
      negTracks_per_coll.reserve(trackIdsThisCollision.size());

      for (auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyFilteredTracks>();
        if (!checkTrack<false>(collision, track) || !isElectron(track)) {
          continue;
        }

        if (track.sign() > 0) {
          posTracks_per_coll.emplace_back(track);
        } else {
          negTracks_per_coll.emplace_back(track);
        }
      }

      fillTrackHistogram<false>(collision, posTracks_per_coll);
      fillTrackHistogram<false>(collision, negTracks_per_coll);

      fillPairInfo<false, EM_EEPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        fillPairInfo<false, EM_EEPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        fillPairInfo<false, EM_EEPairType::kLSnn>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }

      posTracks_per_coll.clear();
      negTracks_per_coll.clear();
      posTracks_per_coll.shrink_to_fit();
      negTracks_per_coll.shrink_to_fit();
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processRec_TTCA, "process reconstructed info only", false); // with TTCA

  // ---------- for MC ----------

  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Partition<MyFilteredTracksMC> posTracksMC = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracksMC = o2::aod::track::signed1Pt < 0.f;
  void processMC_SA(soa::Join<aod::McCollisionLabels, aod::Collisions> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks)
  {
    stored_trackIds.reserve(tracks.size());

    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      // events_bz(d_bz);

      auto posTracks_per_coll = posTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      fillTrackHistogram<true>(collision, posTracks_per_coll);
      fillTrackHistogram<true>(collision, negTracks_per_coll);

      fillPairInfo<true, EM_EEPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        fillPairInfo<true, EM_EEPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        fillPairInfo<true, EM_EEPairType::kLSnn>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processMC_SA, "process reconstructed and MC info ", false);

  void processMC_TTCA(soa::Join<aod::McCollisionLabels, aod::Collisions> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks, aod::TrackAssoc const& trackIndices)
  {
    stored_trackIds.reserve(tracks.size() * 2);

    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      // events_bz(d_bz);

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      std::vector<MyTrackMC> posTracks_per_coll;
      std::vector<MyTrackMC> negTracks_per_coll;
      posTracks_per_coll.reserve(trackIdsThisCollision.size());
      negTracks_per_coll.reserve(trackIdsThisCollision.size());

      for (auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyFilteredTracksMC>();
        if (!checkTrack<true>(collision, track) || !isElectron(track)) {
          continue;
        }

        if (track.sign() > 0) {
          posTracks_per_coll.emplace_back(track);
        } else {
          negTracks_per_coll.emplace_back(track);
        }
      }

      fillTrackHistogram<true>(collision, posTracks_per_coll);
      fillTrackHistogram<true>(collision, negTracks_per_coll);

      fillPairInfo<true, EM_EEPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        fillPairInfo<true, EM_EEPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        fillPairInfo<true, EM_EEPairType::kLSnn>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }

      posTracks_per_coll.clear();
      negTracks_per_coll.clear();
      posTracks_per_coll.shrink_to_fit();
      negTracks_per_coll.shrink_to_fit();
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processMC_TTCA, "process reconstructed info only", false); // with TTCA
};

struct prefilterPrimaryElectron {
  enum class EM_Electron_PF : int {
    kElFromPC = 0, // electron from photon conversion
  };

  Produces<aod::EMPrimaryElectronsPrefilterBit> ele_pfb;

  SliceCache cache;
  Preslice<aod::Tracks> perCol_track = o2::aod::track::collisionId;
  PresliceUnsorted<aod::EMPrimaryElectrons> perCol_ele = o2::aod::emprimaryelectron::collisionId;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  Configurable<float> max_pt_itsonly{"max_pt_itsonly", 0.15, "max pT for ITSonly tracks at PV"};
  Configurable<float> min_dcatopv{"min_dcatopv", -1.f, "DCAxy To PV"};
  Configurable<float> minpt{"minpt", 0.05, "min pt for track for loose track sample"};
  Configurable<float> maxeta{"maxeta", 1.2, "eta acceptance for loose track sample"};
  Configurable<int> mincrossedrows{"mincrossedrows", 40, "min crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max chi2/NclsITS"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -4.0, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 4.0, "max. TPC n sigma for electron inclusion"};
  Configurable<float> slope{"slope", 0.0185, "slope for m vs. phiv"};
  Configurable<float> intercept{"intercept", -0.0280, "intercept for m vs. phiv"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hMvsPhiV_PV", "mass vs. phiv;#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", {HistType::kTH2F, {{90, .0f, TMath::Pi()}, {100, 0, 0.1}}}},
    },
  };

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  template <typename TCollision, typename TTrack>
  bool checkTrack(TCollision const& collision, TTrack const& track)
  {
    if (!track.hasITS()) {
      return false;
    }
    if (track.itsChi2NCl() > maxchi2its) {
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcNSigmaEl() < minTPCNsigmaEl || maxTPCNsigmaEl < track.tpcNSigmaEl()) {
        return false;
      }
      if (track.tpcNClsCrossedRows() < mincrossedrows) {
        return false;
      }

      if (track.tpcChi2NCl() > maxchi2its) {
        return false;
      }
    }

    gpu::gpustd::array<float, 2> dcaInfo;
    auto track_par_recalc = getTrackPar(track);
    std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_recalc, 2.f, matCorr, &dcaInfo);
    getPxPyPz(track_par_recalc, pVec_recalc);

    if (abs(dcaInfo[0]) < min_dcatopv) {
      return false;
    }

    if (isITSonlyTrack(track) && track_par_recalc.getPt() > max_pt_itsonly) {
      return false;
    }
    if (abs(track_par_recalc.getEta()) > maxeta) {
      return false;
    }

    return true;
  }

  template <int loose_track_sign, typename TCollision, typename TTrack1, typename TTrack2>
  bool reconstructPC(TCollision const& collision, TTrack1 const& ele, TTrack2 const& pos)
  {
    if (pos.sign() * ele.sign() > 0) { // reject same sign pair
      return false;
    }

    float mee, phiv = 0;
    gpu::gpustd::array<float, 2> dcaInfo;
    std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz

    if constexpr (loose_track_sign > 0) { // positive track is loose track
      auto track_par_recalc = getTrackPar(pos);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_recalc, 2.f, matCorr, &dcaInfo);
      getPxPyPz(track_par_recalc, pVec_recalc);

      ROOT::Math::PtEtaPhiMVector v1(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
      ROOT::Math::PtEtaPhiMVector v2(track_par_recalc.getPt(), track_par_recalc.getEta(), track_par_recalc.getPhi(), o2::constants::physics::MassElectron);
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      mee = v12.M();
      phiv = getPhivPair(pVec_recalc[0], pVec_recalc[1], pVec_recalc[2], ele.px(), ele.py(), ele.pz(), pos.sign(), ele.sign(), d_bz);
    } else {
      auto track_par_recalc = getTrackPar(ele);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_recalc, 2.f, matCorr, &dcaInfo);
      getPxPyPz(track_par_recalc, pVec_recalc);

      ROOT::Math::PtEtaPhiMVector v1(track_par_recalc.getPt(), track_par_recalc.getEta(), track_par_recalc.getPhi(), o2::constants::physics::MassElectron);
      ROOT::Math::PtEtaPhiMVector v2(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      mee = v12.M();
      phiv = getPhivPair(pos.px(), pos.py(), pos.pz(), pVec_recalc[0], pVec_recalc[1], pVec_recalc[2], pos.sign(), ele.sign(), d_bz);
    }

    fRegistry.fill(HIST("hMvsPhiV_PV"), phiv, mee);

    if (mee < slope * phiv + intercept) {
      return true;
    } else {
      return false;
    }
  }

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& min_dcatopv < nabs(o2::aod::track::dcaXY) && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true;
  using MyFilteredTracks = soa::Filtered<MyTracks>;
  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  Partition<aod::EMPrimaryElectrons> positrons = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<aod::EMPrimaryElectrons> electrons = o2::aod::emprimaryelectron::sign < int8_t(0);

  void processPrefilter_TTCA(aod::Collisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const&, aod::EMPrimaryElectrons const& primaryelectrons, aod::TrackAssoc const& trackIndices)
  {
    std::map<int, uint8_t> pfb_map; // map track.globalIndex -> prefilter bit

    for (auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      auto positrons_per_coll = positrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample
      auto electrons_per_coll = electrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      std::vector<MyTrack> posTracks_per_coll;
      std::vector<MyTrack> negTracks_per_coll;
      posTracks_per_coll.reserve(trackIdsThisCollision.size());
      negTracks_per_coll.reserve(trackIdsThisCollision.size());

      for (auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyTracks>();
        if (!checkTrack(collision, track)) {
          continue;
        }

        if (track.sign() > 0) {
          posTracks_per_coll.emplace_back(track);
        } else {
          negTracks_per_coll.emplace_back(track);
        }
      }

      for (auto& empos : positrons_per_coll) {
        // auto pos = tracks.rawIteratorAt(empos.trackId()); // use rawIterator, if the table is filtered.
        for (auto& ele : negTracks_per_coll) {
          if (!checkTrack(collision, ele)) {
            continue;
          }
          if (empos.trackId() == ele.globalIndex()) {
            continue;
          }
          bool isPC = reconstructPC<-1>(collision, ele, empos);
          if (isPC) {
            pfb_map[empos.globalIndex()] |= (uint8_t(1) << static_cast<int>(EM_Electron_PF::kElFromPC));
          }
        } // end of loose electron loop
      }   // end of signal positon loop

      for (auto& emele : electrons_per_coll) {
        // auto ele = tracks.rawIteratorAt(emele.trackId()); // use rawIterator, if the table is filtered.
        for (auto& pos : posTracks_per_coll) {
          if (!checkTrack(collision, pos)) { // track cut is applied to loose sample
            continue;
          }
          if (emele.trackId() == pos.globalIndex()) {
            continue;
          }
          bool isPC = reconstructPC<+1>(collision, emele, pos);
          if (isPC) {
            pfb_map[emele.globalIndex()] |= (uint8_t(1) << static_cast<int>(EM_Electron_PF::kElFromPC));
          }
        } // end of loose positon loop
      }   // end of signal electron loop

      posTracks_per_coll.clear();
      negTracks_per_coll.clear();
      posTracks_per_coll.shrink_to_fit();
      negTracks_per_coll.shrink_to_fit();
    } // end of collision loop

    for (auto& ele : primaryelectrons) {
      ele_pfb(pfb_map[ele.globalIndex()]);
    }

    pfb_map.clear();
  }
  PROCESS_SWITCH(prefilterPrimaryElectron, processPrefilter_TTCA, "process prefilter with TTCA", false);

  void processPrefilter_SA(aod::Collisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const&, aod::EMPrimaryElectrons const& primaryelectrons)
  {
    std::map<int, uint8_t> pfb_map; // map track.globalIndex -> prefilter bit

    for (auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache); // loose track sample
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache); // loose track sample

      auto positrons_per_coll = positrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample
      auto electrons_per_coll = electrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample

      for (auto& [ele, empos] : combinations(CombinationsFullIndexPolicy(negTracks_per_coll, positrons_per_coll))) {
        // auto pos = tracks.rawIteratorAt(empos.trackId()); // use rawIterator, if the table is filtered.
        if (!checkTrack(collision, ele)) { // track cut is applied to loose sample
          continue;
        }
        if (empos.trackId() == ele.globalIndex()) {
          continue;
        }

        bool isPC = reconstructPC<-1>(collision, ele, empos);
        if (isPC) {
          pfb_map[empos.globalIndex()] |= (uint8_t(1) << static_cast<int>(EM_Electron_PF::kElFromPC));
        }
      }

      for (auto& [pos, emele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, electrons_per_coll))) {
        // auto ele = tracks.rawIteratorAt(emele.trackId()); // use rawIterator, if the table is filtered.
        if (!checkTrack(collision, pos)) { // track cut is applied to loose sample
          continue;
        }
        if (emele.trackId() == pos.globalIndex()) {
          continue;
        }
        bool isPC = reconstructPC<+1>(collision, emele, pos);
        if (isPC) {
          pfb_map[emele.globalIndex()] |= (uint8_t(1) << static_cast<int>(EM_Electron_PF::kElFromPC));
        }
      }
    } // end of collision loop

    for (auto& ele : primaryelectrons) {
      ele_pfb(pfb_map[ele.globalIndex()]);
    }

    pfb_map.clear();
  }
  PROCESS_SWITCH(prefilterPrimaryElectron, processPrefilter_SA, "process prefilter standalone", false);

  void processDummy(aod::EMPrimaryElectrons const& primaryelectrons)
  {
    for (int i = 0; i < primaryelectrons.size(); i++) {
      ele_pfb(0);
    }
  }
  PROCESS_SWITCH(prefilterPrimaryElectron, processDummy, "process dummy", true);
};

struct associateAmbiguousElectron {
  Produces<aod::EMAmbiguousElectronSelfIds> em_amb_ele_ids;

  SliceCache cache;
  PresliceUnsorted<aod::EMPrimaryElectrons> perTrack = o2::aod::emprimaryelectron::trackId;
  std::vector<int> ambele_self_Ids;

  void process(aod::EMPrimaryElectrons const& electrons)
  {
    for (auto& electron : electrons) {
      auto electrons_with_same_trackId = electrons.sliceBy(perTrack, electron.trackId());
      ambele_self_Ids.reserve(electrons_with_same_trackId.size());
      for (auto& amp_ele : electrons_with_same_trackId) {
        if (amp_ele.globalIndex() == electron.globalIndex()) { // don't store myself.
          continue;
        }
        ambele_self_Ids.emplace_back(amp_ele.globalIndex());
      }
      em_amb_ele_ids(ambele_self_Ids);
      ambele_self_Ids.clear();
      ambele_self_Ids.shrink_to_fit();
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<skimmerPrimaryElectron>(cfgc, TaskName{"skimmer-primary-electron"}),
    adaptAnalysisTask<prefilterPrimaryElectron>(cfgc, TaskName{"prefilter-primary-electron"}),
    adaptAnalysisTask<associateAmbiguousElectron>(cfgc, TaskName{"associate-ambiguous-electron"})};
}
