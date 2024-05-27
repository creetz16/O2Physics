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

/// \file qaKFXicToXiPiPi.cxx
/// \brief Task to reconstruct XicPlus --> (Xi --> (Lam--> p pi) pi) pi  with the KFParticle software package
///
/// \author Carolina Reetz <c.reetz@cern.ch>, University of Heidelberg

#ifndef HomogeneousField
#define HomogeneousField
#endif

// KFParticle
#include <KFParticleBase.h>
#include <KFParticle.h>
#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFVertex.h>

// Root
#include <TPDGCode.h>

// O2
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"

// O2Physics
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/KFparticle/KFUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "Tools/KFparticle/qaKFXicToXiPiPi.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;
using namespace o2::constants::physics;

using KFCascadesLinked = soa::Join<aod::Cascades, aod::KFCascDataLink>;
using KFCascFull = soa::Join<aod::KFCascDatas, aod::KFCascCovs>;
using TaggedCascades = soa::Join<aod::Cascades, aod::CascTags>;

struct qaKFXicToXiPiPibuilder {

  // ============ Tables =============
  Produces<aod::> row;
  Produces<aod::> row;

  // =========== Filtering and preslicing ==========
  Preslice<KFCascFull> cascadesPerCollision = aod::cascdata::collisionId;
  Preslice<TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  // =========== Data members ==========
  Service<ccdb::BasicCCDBManager> ccdb;o2::vertexing::DCAFitterN<3> fitter3body;
  base::MatLayerCylSet* lut = nullptr;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int mRunNumber;
  float mBz;
  // PDG masses
  double massP{0.};
  double massPi{0.};
  double massXi{0.};
  double massLambda{0.};

  // ============ Configurables ============
  // general
  Configurable<bool> isRun3{"isRun3", true, "Run3 dataset"};
  Configurable<bool> fillQAHistograms{"fillQAHistograms", true, "fill histograms"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  Configurable<bool> doDCAFitterPreMinimum{"doDCAFitterPreMinimum", true, "Do DCAFitter pre-optimization before KF fit to include material corrections"};
  Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices in DCA fitter"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs in DCA fitter"};

  // CCDB options
  struct : ConfigurableGroup {
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } ccdbConfigurations;

  // track selections
  struct : ConfigurableGroup {
    Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
    Configurable<float> d_pTMin{"d_pTMin", 0., "Minimum momentum for tracks"};
    Configurable<int> d_crossedrows{"d_crossedrows", 70, "Minimum number of crossed rows in the TPC"};
    Configurable<float> d_eta{"d_eta", 0.8, "Maximum eta"};
    Configurable<float> d_dcaXYTrackPV{"d_dcaXYTrackPV", 2., "DCA XY of the daughter tracks to the PV"};
    Configurable<float> d_dcaZTrackPV{"d_dcaZTrackPV", 10., "DCA Z of the daughter tracks to the PV"};
  } trackConfigurations;

  // TPC PID selections
  struct : ConfigurableGroup {
    Configurable<double> d_nSigmaTPCpi{"nSigmaTPCpi", 5., "TPC Nsigma pion cut"};
    Configurable<double> d_nSigmaTPCpr{"nSigmaTPCpr", 5., "TPC Nsigma proton cut"};
  } pidConfigurations;

  // cascadebuilder selections
  struct : ConfigurableGroup {
    Configurable<float> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "Minimum cascade cosPA to the PV"};
    Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.9, "Minimum cascade decay radius"};
    Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "Maximum DCA between cascade daughters"};
    Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.05, "Minimum DCA bach to PV"};
    Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "Maximum distance from V0 mass"};
  }

  // candidate preselections
  struct : ConfigurableGroup {
    Configurable<float> v0TransvRadius{"v0TransvRadius", 1.0, "Minimum V0 radius in xy plane"};           // 1.2 (xi) and 1.1 (omega) in run2
    Configurable<float> cascTransvRadius{"cascTransvRadius", 0.4, "Minimum cascade radius in xy plane"};  // 0.5 cm (xi) and 0.6 (omega) in run2
    Configurable<float> dcaBachToPv{"dcaBachToPv", 0.03, "Minimum DCA bach to PV"};                       // 0.04 in run2
    Configurable<float> dcaV0ToPv{"dcaV0ToPv", 0.02, "Minimum DCA V0 to PV"};                             // 0.03 in run2
    Configurable<double> v0CosPA{"v0CosPA", 0.95, "Minimum V0 cosPA to the PV"};                          // 0.97 in run2 - KEEP LOSE to re-cut after PVRefit! - double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<double> cascCosPA{"cascCosPA", 0.95, "Minimum cascade cosPA to the PV"};                 // 0.97 in run2 - KEEP LOSE to re-cut after PVRefit! - double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> dcaV0Dau{"dcaV0Dau", 1.0, "Maximum DCA betweem V0 daughters"};                    // conservative, a cut ar 1.0 should also be fine
    Configurable<float> dcaCascDau{"dcaCascDau", 1.0, "Maximum DCA between cascade daughters"};           // conservative, a cut ar 1.0 should also be fine
    Configurable<float> dcaNegToPv{"dcaNegToPv", 0.05, "Minimum DCA Neg To PV"};                          // 0.06 in run2
    Configurable<float> dcaPosToPv{"dcaPosToPv", 0.05, "Minimum DCA Pos To PV"};                          // 0.06 in run2
    Configurable<float> v0MassWindow{"v0MassWindow", 0.01, "Maximum distance from V0 mass"};              // 0.008 in run2
    Configurable<float> cascadeMassWindow{"cascadeMassWindow", 0.01, "Maximum distance from cascade mass"};
  }

  // Helper struct to pass candidate information
  /// TODO: update for Xic candidates
  struct {
    int v0Id;
    int positiveId;
    int negativeId;
    int bachelorId;
    int pi0Id;
    int pi1Id;
    float positiveX;
    float negativeX;
    float bachelorX;
    float pi0Id;
    float pi1Id;
    int charge;
    std::array<float, 3> pos;
    std::array<float, 3> bachP;
    float dcacascdau;
    float bachDCAxy;
    float cosPA;
    float cascradius;
    float cascDCAxy; // cascade DCA xy (with bending)
    float cascDCAz;  // cascade DCA z
    std::array<float, 3> v0pos;
    std::array<float, 3> v0mompos;
    std::array<float, 3> v0momneg;
    std::array<float, 3> v0pospos;
    std::array<float, 3> v0posneg;
    std::array<float, 3> cascademom;
    std::array<float, 3> kfv0mom;
    float v0dcadau;
    float v0dcapostopv;
    float v0dcanegtopv;
    float mXi;
    float mOmega;
    float yXi;
    float yOmega;
    float bachBaryonCosPA;
    float bachBaryonDCAxyToPV;
    float kfMLambda;
    float kfV0Chi2;
    float kfCascadeChi2;
    std::array<float, 21> kfCascadeCov;
    std::array<float, 21> kfV0Cov;
    std::array<float, 21> kfV0DauPosCov;
    std::array<float, 21> kfV0DauNegCov;
  } xiccandidate;

  void init(InitContext const&)
  {
    mRunNumber = 0;
    mBz = 0.;
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // Particle masses
    massP = o2::constants::physics::MassProton;
    massPi = o2::constants::physics::MassPiPlus;
    massXi = o2::constants::physics::MassXiMinus;
    massLambda = o2::constants::physics::MassLambda0;

    if (fillQAHistograms)
    {
      /// QA histogram axis definitions
      const AxisSpec axVtxXY{1000, -0.1f, 0.1f};
      const AxisSpec axVtxZ{500, -15.f, 15.f};
      const AxisSpec axVtxCov{500, -0.00007f, 0.00007f};
      const AxisSpec axMomentum{2000, 0.f, 5.0f};
      const AxisSpec axNsigma{1200, -20.0f, 6.0f};
      const AxisSpec axDCA{2500, -0.6f, 0.6f};

      /// QA histogram definitions: events
      auto hEventCounter = hist.add<TH1>("Events/hEventCounter", "Event selection flow", kTH1F, {{3, 0.5f, 3.5f}});
      hEventCounter->GetXaxis()->SetBinLabel(hEventCounter->FindBin(1), "Sel8");
      hEventCounter->GetXaxis()->SetBinLabel(hEventCounter->FindBin(2), "|Vtx_{z}|<10cm");
      hist.add("Events/hVertexX", "PV position x; Vtx_{x} (cm); entries", kTH1F, {axVtxXY});
      hist.add("Events/hVertexY", "PV position y; Vtx_{y} (cm); entries", kTH1F, {axVtxXY});
      hist.add("Events/hVertexXY", "PV position xy; Vtx_{x} (cm); Vtx_{y} (cm)", kTH2F, {axVtxXY, axVtxXY});
      hist.add("Events/hVertexZ", "PV position z; Vtx_{z} (cm); entries", kTH1F, {axVtxZ});
      hist.add("Events/hVertexCovXX", "PV Cov_{xx}; Cov_{xx} (cm^{2}); entries", kTH1F, {axVtxCov});
      hist.add("Events/hVertexCovYY", "PV Cov_{yy}; Cov_{yy} (cm^{2}); entries", kTH1F, {axVtxCov});
      hist.add("Events/hVertexCovZZ", "PV Cov_{zz}; Cov_{zz} (cm^{2}); entries", kTH1F, {axVtxCov});
      hist.add("Events/hVertexCovXY", "PV Cov_{xy}; Cov_{xy} (cm^{2}); entries", kTH1F, {axVtxCov});
      hist.add("Events/hVertexCovXZ", "PV Cov_{xz}; Cov_{xz} (cm^{2}); entries", kTH1F, {axVtxCov});
      hist.add("Events/hVertexCovYZ", "PV Cov_{yz}; Cov_{yz} (cm^{2}); entries", kTH1F, {axVtxCov});

      /// QA histogram definitions: tracks
      auto hCandidateCounter = hist.add<TH1>("hCandidateCounter", "Candidate selecton flow", kTH1F, {{10, 0.5f, 10.5f}});
      hCandidateCounter->GetXaxis()->SetBinLabel(hCandidateCounter->FindBin(1), "In acceptance");
      hCandidateCounter->GetXaxis()->SetBinLabel(hCandidateCounter->FindBin(2), "Collision ID");
      hCandidateCounter->GetXaxis()->SetBinLabel(hCandidateCounter->FindBin(3), "Charge");
      hCandidateCounter->GetXaxis()->SetBinLabel(hCandidateCounter->FindBin(4), "Selected daughters");
      hCandidateCounter->GetXaxis()->SetBinLabel(hCandidateCounter->FindBin(5), "Has SV");
      hCandidateCounter->GetXaxis()->SetBinLabel(hCandidateCounter->FindBin(6), "Sel Xic geo");
      hCandidateCounter->GetXaxis()->SetBinLabel(hCandidateCounter->FindBin(7), "Sel Xic topo");
      hist.add("Tracks/hTPCNSigmaPr", "nSigma TPC proton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(p)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hTPCNSigmaPi", "nSigma TPC pion from #Lambda; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(pi)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hTPCNSigmaPiFromXi", "nSigma TPC pion from #Xi; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(pi)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hTPCNSigmaPi0", "nSigma TPC pion0 from #Xi_{c}^{+}; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(pi)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hTPCNSigmaPi1", "nSigma TPC pion1 from Xic; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(pi)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hTOFNSigmaPr", "nSigma TOF proton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TOF}(p)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hTOFNSigmaPi", "nSigma TOF pion from #Lambda; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TOF}(pi)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hTOFNSigmaPiFromXi", "nSigma TOF pion from #Xi; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TOF}(pi)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hTOFNSigmaPi0", "nSigma TOF pion0 from #Xi_{c}^{+}; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TOF}(pi)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hTOFNSigmaPi1", "nSigma TOF pion1 from #Xi_{c}^{+}; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TOF}(pi)", kTH2D, {axMomentum, axNsigma});
      hist.add("Tracks/hDCAxy", "DCA_{xy} to PV; DCA_{xy} (cm); entries", kTH1D, {axDCA});
      hist.add("Tracks/hDCAz", "DCA_{z} to PV; DCA_{z} (cm); entries", kTH1D, {axDCA});

      // Cascade mass spectra
      hist.add("Cascade/hMassXiMinus", "hMassXiMinus", {HistType::kTH1F, {{400, 1.122f, 1.522f, "Inv. Mass (GeV/c^{2})"}}});
      hist.add("Cascade/hMassXiPlus", "hMassXiPlus", {HistType::kTH1F, {{400, 1.122f, 1.522f, "Inv. Mass (GeV/c^{2})"}}});
      hist.add("Cascade/hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH1F, {{400, 1.472f, 1.872f, "Inv. Mass (GeV/c^{2})"}}});
      hist.add("Cascade/hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH1F, {{400, 1.472f, 1.872f, "Inv. Mass (GeV/c^{2})"}}});

      // Cascade topology
      hist.add("Cascade/hV0Radius", "hV0Radius", {HistType::kTH1F, {{500, 0.0, 100.0, "radius (cm)"}}});
      hist.add("Cascade/hCascRadius", "hCascRadius", {HistType::kTH1F, {{500, 0.0, 100.0, "radius (cm)"}}});
      hist.add("Cascade/hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f, "cos(PA)"}}});
      hist.add("Cascade/hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f, "cos(PA)"}}});
      hist.add("Cascade/hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "DCA to PV (cm)"}}});
      hist.add("Cascade/hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "DCA to PV (cm)"}}});
      hist.add("Cascade/hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "DCA to PV (cm)"}}});
      hist.add("Cascade/hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "DCA to PV (cm)"}}});
      hist.add("Cascade/hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{500, 0.0f, 5.0f, "DCA(V0 daughters) (cm)"}}});
      hist.add("Cascade/hDCACascDau", "hDCACascDau", {HistType::kTH1F, {{500, 0.0f, 5.0f, "DCA (cascade daughters) (cm)"}}});
      hist.add("Cascade/hLambdaMass", "hLambdaMass", {HistType::kTH1F, {{400, 0.916f, 1.316f, "Inv. Mass V0 (GeV/c^{2})"}}});
    } // if fillQAHistograms

  } // end init

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    auto timestamp = bc.timestamp();
    o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, timestamp);

    if (!grpmag) {
      LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);

    // Fetch magnetic field from ccdb for current collision
    mBz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << mBz << " kG";
    mRunNumber = bc.runNumber();
    LOG(info) << "Magnetic field has been initialised";

    /// Set magnetic field for KF vertexing
    #ifdef HomogeneousField
        KFParticle::SetField(mBz);
    #endif
    LOG(info) << "Magnetic field has been set for KF";

    if (useMatCorrType == 2) {
      base::Propagator::Instance()->setMatLUT(lut);
      LOG(info) << "LUT has been set for propagator";
    }
  } // end initCCDB

  template <typename TCascade>
  bool isPreselectedCascade(const TCascade& cascade, const float& pvx, const float& pvy, const float& pvz)
  {
    if (cascade.v0cosPA(pvx, pvy, pvz) > v0CosPA &&
        cascade.casccosPA(pvx, pvy, pvz) > cascCosPA &&
        cascade.dcacascdaughters() < dcaCascDau &&
        cascade.dcaV0daughters() < dcaV0Dau &&
        std::abs(cascade.dcanegtopv()) > dcaNegToPv &&
        std::abs(cascade.dcapostopv()) > dcaPosToPv &&
        std::abs(cascade.dcabachtopv()) > dcaBachToPv &&
        std::abs(cascade.dcav0topv(pvx, pvy, pvz)) > dcaV0ToPv &&
        cascade.v0radius() > v0TransvRadius &&
        cascade.cascradius() > cascTransvRadius &&
        std::abs(cascade.mLambda() - massLambda) < v0MassWindow) {

      if (fillQAHistograms) {
        hist.fill(HIST("Cascade/hV0Radius"), cascade.v0radius());
        hist.fill(HIST("Cascade/hCascRadius"), cascade.cascradius());
        hist.fill(HIST("Cascade/hV0CosPA"), cascade.v0cosPA(pvx, pvy, pvz));
        hist.fill(HIST("Cascade/hCascCosPA"), cascade.casccosPA(pvx, pvy, pvz));
        hist.fill(HIST("Cascade/hDCAPosToPV"), cascade.dcapostopv());
        hist.fill(HIST("Cascade/hDCANegToPV"), cascade.dcanegtopv());
        hist.fill(HIST("Cascade/hDCABachToPV"), cascade.dcabachtopv());
        hist.fill(HIST("Cascade/hDCAV0ToPV"), cascade.dcav0topv(pvx, pvy, pvz));
        hist.fill(HIST("Cascade/hDCAV0Dau"), cascade.dcaV0daughters());
        hist.fill(HIST("Cascade/hDCACascDau"), cascade.dcacascdaughters());
        hist.fill(HIST("Cascade/hLambdaMass"), cascade.mLambda());
        if (cascade.sign() < 0) {
          hist.fill(HIST("Cascade/hMassXiMinus"), cascade.mXi());
          hist.fill(HIST("Cascade/hMassOmegaMinus"), cascade.mOmega());
        } else {
          hist.fill(HIST("Cascade/hMassXiPlus"), cascade.mXi());
          hist.fill(HIST("Cascade/hMassOmegaPlus"), cascade.mOmega());
        }
      }
      return true;
    }
    
    return false;
  }

  template <class TTrackTo, typename TCascade, typename TCollision>
  bool buildXicCandidate(TCascade const& cascade, TCollision const& collision)
  {
    // Track casting of cascade daughter tracks
    auto bachTrack = cascade.template bachelor_as<TTrackTo>();
    auto posTrack = cascade.template posTrack_as<TTrackTo>();
    auto negTrack = cascade.template negTrack_as<TTrackTo>();

    // Overall Xic candidate charge
    xiccandidate.charge = bachTrack.signed1Pt() > 0 ? -1 : +1;

    hist.fill(HIST("hCandidateCounter"), 0.5); // all cascade candidates

    // check collision IDs of cascade daughters
    if (posTrack.collisionId() != negTrack.collisionId() || 
        posTrack.collisionId() != bachTrack.collisionId() || 
        negTrack.collisionId() != bachTrack.collisionId()) {
      return false;
    }
    hist.fill(HIST("hCandidateCounter"), 1.5);

    // check global IDs of cascade daughters (exclude duplicates)
    if (posTrack.globalIndex() == negTrack.globalIndex() || 
        posTrack.globalIndex() == bachTrack.globalIndex() || 
        negTrack.globalIndex() == bachTrack.globalIndex()) {
      return false;
    }
    hist.fill(HIST("hCandidateCounter"), 2.5);

    // preselect cascade
    if (!isPreselectedCascade(cascade, collision.posX(), collision.posY(), collision.posZ()))
      return false;
    hist.fill(HIST("hCandidateCounter"), 3.5);

    // create KFParticle cascade object from cascdata table info
    std::array<float, 3> xyzCasc = {cascade.x(), cascade.y(), cascade.z()};
    std::array<float, 3> pxpypzCasc = {cascade.px(), cascade.py(), cascade.pz()};
    float xyzpxpypzCasc[6];
    for (int i{0}; i < 3; ++i) {
      xyzpxpypzCasc[i] = xyz[i];
      xyzpxpypzCasc[i + 3] = pxpypz[i];
    }
    std::array<float, 21> covCasc;
    for (int i{0}; i < 21; ++i) {
      covCasc[i] = cascade.kfTrackCovMat()[i];
    }
    KFParticle kfCascXi;
    float Mini, SigmaMini, M, SigmaM;
    kfCascXi.GetMass(Mini, SigmaMini);
    LOG(debug) << "Cascade KFParticle mass before creation: " << Mini << " +- " << SigmaMini;
    try {
      kfCascXi.Create(xyzpxpypzCasc, covCasc.data(), xiccandidate.charge, massXi);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create cascade KFParticle from cascdata table" << e.what();
    }
    kfCascXi.GetMass(M, SigmaM);
    LOG(debug) << "Cascade KFParticle mass after creation: " << M << " +- " << SigmaM;

    // combine cascade with pion tracks
    /// TODO: continue here!
    auto groupedPionTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
    for (auto trackIdPion0 = groupedPionTrackIndices.begin(); )





    return true;
  }

  template <typename TCollision>
  void buildCandidateTable(TCollision const& collision)
  {
    /// TODO: fill candidate table here

  }
  }

  void process(aod::Collision const& collision, 
              soa::Filtered<KFCascFull> const& cascades, 
              FullTracksExtIU const&, 
              aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    hist.fill(HIST("hEventCounter"), 0.5);

    // grouped cascades for this collision
    const uint64_t collIdx = collision.globalIndex();
    auto CascadeTable_thisCollision = cascades.sliceBy(cascadesPerCollision, collIdx);
    // loop over grouped cascdatas
    for (auto& cascade : cascades) {

      // create Xic+ candidates starting from cascades
      bool validXicCandidate = buildXicCandidate<TTrackTo>(cascade, collision);
      if (!validXicCandidate)
        continue;
      
      // fill candidate table
      buildCandidateTable(collision);

    }

  }






}
