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
/// \file   nucleiRecoQA.cxx
/// \author 
/// \brief  task to study nuclei momentum reconstruction

#include <PWGLF/DataModel/NucleiRecoQA.h>
#include <CCDB/BasicCCDBManager.h>
#include <string>
#include "TableHelper.h"

/// O2
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

/// O2Physics
#include "CommonDataFormat/InteractionRecord.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
//#include "Tools/KFparticle/KFUtilities.h"

/// KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

#ifndef HomogeneousField
#define HomogeneousField
#endif

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"He3"};
} // namespace

std::string trackdirs[] = {"Tracks", "Deuterons", "Tritons", "Helium3"};
//const char* trackdirs[4] = {"Tracks", "Deuterons", "Tritons", "Helium3"};

struct NucleiRecoQA {
    /// create output tables
    Produces<o2::aod::TableCollisions> tableCollisions;
    Produces<o2::aod::TableTracks> tableTracks;

    /// general
    int runNumber = 0.;
    double magneticField = 0.;
    Configurable<bool> isRun3{"isRun3", true, "Flag to run with Run 3 data"};
    Configurable<bool> isMC{"isMC", false, "Flag to run over MC"};
    Configurable<bool> writeCollTable{"writeCollTable", false, "flag to write collision properties to a table"};
    Configurable<bool> writeTrackTable{"writeTrackTable", false, "flag to write track properties to a table"};

    // bethe bloch parameters
    Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
    Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

    // PDG codes
    Configurable<int> hyperPdg{"hyperPDG", 1010010030, "PDG code of the hyper-mother (could be 3LamH or 4LamH)"};
    Configurable<int> heDauPdg{"heDauPDG", 1000020030, "PDG code of the helium (could be 3He or 4He)"};

    /// settings for CCDB
    Service<o2::ccdb::BasicCCDBManager> ccdb;
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"}; /// what is this lut??? -> material budget distribution?
    Configurable<std::string> ccdbPathGeo{"ccdbPathGeo", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"}; /// grp? -> needed for magnet field in Run 3?
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
    o2::base::MatLayerCylSet* lut; 

    /// event selection
    Configurable<bool> eventSel{"eventSel", true, "basic event selection"};

    /// track selection
    Configurable<int> trackSelType{"trackSelType", 1, "option for track selection: 0 = no cut, 1 = kGlobalTrack, 2 = kGlobalTrackWoPtEta, 3 = kGlobalTrackWoDCA, 4 = kQualityTracks, 5 = kInAcceptanceTracks"};
    Configurable<bool> trackSel{"trackSel", true, "basic track selection"};
    Configurable<float> eta{"eta", 0.8, "eta"}; // expression table
    Configurable<float> tpcNClsFound{"tpcNClsFound", 100, "Number of TPC clusters"}; // dynamic
    Configurable<float> dEdxMin{"dEdxMin", 240, "Minimum dE/dx"}; 
    Configurable<float> ptMin{"ptMin", 1.5, "Minimum track pT"}; 

    /// histogram configurables
    Configurable<double> maxPt{"maxPt", 5.0, "max for pT axis"};

    Filter trackFilter = (trackSelType.value == 0) || 
                         ((trackSelType.value == 1) && requireGlobalTrackInFilter()) ||
                         ((trackSelType.value == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                         ((trackSelType.value == 3) && requireGlobalTrackWoDCAInFilter()) ||
                         ((trackSelType.value == 4) && requireQualityTracksInFilter()) ||
                         ((trackSelType.value == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

    /// histogram registry
    HistogramRegistry hist;

    /// initialise magnetic field 
    void initMagneticFieldCCDB(o2::aod::BCsWithTimestamps::iterator const& bc, int& mRunNumber,
                             o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, std::string ccdbPathGrp, o2::base::MatLayerCylSet* lut,
                             bool isRun3) {      
        LOGF(info, "====== initCCDB function called (isRun3==%d)", isRun3);
        if (!isRun3) { /// Run 2 GRP object
            o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp()); /// get GRP object from specified CCDB path for current timestamp
            if (grpo == nullptr) {
                LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
            }
            o2::base::Propagator::initFieldFromGRP(grpo); /// intialises field using the fetched object
            o2::base::Propagator::Instance()->setMatLUT(lut); /// sets the material budget instance to the above specified defined MatLayerCylSet (do I need this???)
            LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());
        } else { /// Run 3 GRP object
            o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrp, bc.timestamp()); /// gets magnetic field for Run 3 for current timestamp
            if (grpo == nullptr) {
                LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
            }
            o2::base::Propagator::initFieldFromGRP(grpo);
            o2::base::Propagator::Instance()->setMatLUT(lut);
            LOGF(info, "Setting magnetic field to current %f A for run %d from its GRP CCDB object (type o2::parameters::GRPMagField)", grpo->getL3Current(), bc.runNumber());
        }
        mRunNumber = bc.runNumber(); /// set mRunNumber variable to current run number
    } /// end initMagneticFieldCCDB

    void init(InitContext const&) { /// what is actually passed to the function ?
        //if ((processData == true) && (processMC == true)) LOGF(fatal, "Cannot enable processData and processMC at the same time. Please choose one.");

        ccdb->setURL(ccdbUrl.value); /// set CCDB url
        ccdb->setCaching(true); /// enabling object chaching
        ccdb->setLocalObjectValidityChecking(); /// what happens here?
        lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut)); /// "fill" lut pointer with actual material budget info from specified CCDB path
        if (!o2::base::GeometryManager::isGeometryLoaded()) {
            ccdb->get<TGeoManager>(ccdbPathGeo); /// load geometry from specified CCDB path (what does this mean? what happens?)
        }

        /// set runNumber variable to 0 (WHY???)
        runNumber = 0; 

        /// axis definitions
        const AxisSpec axVtxXY{1000, -0.1, 0.1};
        const AxisSpec axVtxZ{500, -15., 15.};
        const AxisSpec axVtxCov{500, -0.00007, 0.00007};
        const AxisSpec axNContrib{100, -0.5, 99.5};
        const AxisSpec axMultiplicity{100, -0.5, 99.5};
        const AxisSpec axSignedMomentum{2000, -maxPt, maxPt};
        const AxisSpec axMomentum{2000, 0., maxPt};
        const AxisSpec axNsigma{1200, -20.0, 6.0};
        const AxisSpec axDCA{2500, -0.6, 0.6};

        /// histos definitions: events
        if (isMC) {
            auto hEventSels = hist.add<TH1>("Events/hEventSels", "Event Selections", kTH1D, {{5, 0.5, 5.5}});
            hEventSels->GetXaxis()->SetBinLabel(hEventSels->FindBin(1), "All events");
            hEventSels->GetXaxis()->SetBinLabel(hEventSels->FindBin(2), "with MC collision");
            hEventSels->GetXaxis()->SetBinLabel(hEventSels->FindBin(3), "Sel8");
            hEventSels->GetXaxis()->SetBinLabel(hEventSels->FindBin(4), "|Vtx_{z}|<10cm");
            hEventSels->GetXaxis()->SetBinLabel(hEventSels->FindBin(5), "numContrib>1");
        } else {
            auto hEventSels = hist.add<TH1>("Events/hEventSels", "Event Selections", kTH1D, {{4, 0.5, 4.5}});
            hEventSels->GetXaxis()->SetBinLabel(hEventSels->FindBin(1), "All events");
            hEventSels->GetXaxis()->SetBinLabel(hEventSels->FindBin(2), "Sel8");
            hEventSels->GetXaxis()->SetBinLabel(hEventSels->FindBin(3), "|Vtx_{z}|<10cm");
            hEventSels->GetXaxis()->SetBinLabel(hEventSels->FindBin(4), "numContrib>1");
        }
        hist.add("Events/hVertexX", "PV position x; Vtx_{x} (cm); entries", kTH1D, {axVtxXY});
        hist.add("Events/hVertexY", "PV position y; Vtx_{y} (cm); entries", kTH1D, {axVtxXY});
        hist.add("Events/hVertexXY", "PV position xy; Vtx_{x} (cm); Vtx_{y} (cm)", kTH2D, {axVtxXY, axVtxXY});
        hist.add("Events/hVertexZ", "PV position z; Vtx_{z} (cm); entries", kTH1D, {axVtxZ});
        hist.add("Events/hVertexCovXX", "PV Cov_{xx}; Cov_{xx} (cm^{2}); entries", kTH1D, {axVtxCov});
        hist.add("Events/hVertexCovYY", "PV Cov_{yy}; Cov_{yy} (cm^{2}); entries", kTH1D, {axVtxCov});
        hist.add("Events/hVertexCovZZ", "PV Cov_{zz}; Cov_{zz} (cm^{2}); entries", kTH1D, {axVtxCov});
        hist.add("Events/hVertexCovXY", "PV Cov_{xy}; Cov_{xy} (cm^{2}); entries", kTH1D, {axVtxCov});
        hist.add("Events/hVertexCovXZ", "PV Cov_{xz}; Cov_{xz} (cm^{2}); entries", kTH1D, {axVtxCov});
        hist.add("Events/hVertexCovYZ", "PV Cov_{yz}; Cov_{yz} (cm^{2}); entries", kTH1D, {axVtxCov});
        hist.add("Events/hVertexNcontrib", "Real Number of contributors; real N contributors; entries", kTH1D, {axNContrib});
        hist.add("Events/hMultiplicity", "Event multiplicity; multiplicity; entries", kTH1D, {axMultiplicity});
        hist.add("Events/hMultiplicityDeuterons", "Event multiplicity; # of deuterons; entries", kTH1D, {axMultiplicity});
        hist.add("Events/hMultiplicityTritons", "Event multiplicity; # of tritons; entries", kTH1D, {axMultiplicity});
        hist.add("Events/hMultiplicityHelium3", "Event multiplicity; # of helium3; entries", kTH1D, {axMultiplicity});

        /// histos definitions: tracks
        // for (unsigned int i = 0; i < sizeof(trackdirs); i++) {
        //     auto& dir = trackdirs[i];
        hist.add("Tracks/hpx", "Track momentum x; #it{p}_{x} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tracks/hpy", "Track momentum y; #it{p}_{y} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tracks/hpz", "Track momentum z; #it{p}_{z} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tracks/hpt", "Transverse momentum; #it{p}_{T} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tracks/hp", "Track momentum; #it{p} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tracks/hSign", "Track sign; ; entries", kTH1D, {{2, -2., 2.}});
        hist.add("Tracks/hEta", "Track eta; #eta; entries", kTH1D, {{500, -1., 1.}});
        hist.add("Tracks/hPhi", "Track phi; #phi; entries", kTH1D, {{500, 0, 6.5}});
        hist.add("Tracks/hEtaPhi", "Track eta vs phi; #eta; #phi", kTH2D, {{500, -1., 1.}, {500, 0, 6.5}});
        hist.add("Tracks/hTPCinnerParam", "TPC mommentum; #it{p}_{TPC} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tracks/hTPCSignal", "TPC signal vs TPC momentum; #it{p}_{TPC}/Z (GeV/#it{c}); d#it{E}/d#it{x} (A.U.)", kTH2D, {axSignedMomentum, {1000, 20., 500.}});
        hist.add("Tracks/hTPCNSigmaPr", "nSigma TPC proton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(p)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tracks/hTPCNSigmaDe", "nSigma TPC deuteron; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(de)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tracks/hTPCNSigmaTr", "nSigma TPC triton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(tr)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tracks/hTPCNSigmaHe", "nSigma TPC helium; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(he)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tracks/hTOFNSigmaPr", "nSigma TOF proton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TOF}(p)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tracks/hTOFNSigmaDe", "nSigma TOF deuteron; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TOF}(de)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tracks/hTOFNSigmaTr", "nSigma TOF triton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TOF}(tr)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tracks/hTOFNSigmaHe", "nSigma TOF helium; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TOF}(he)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tracks/hDCAxy", "DCA_{xy} to PV; DCA_{xy} (cm); entries", kTH1D, {axDCA});
        hist.add("Tracks/hDCAz", "DCA_{z} to PV; DCA_{z} (cm); entries", kTH1D, {axDCA});
        hist.add("Tracks/hDCAxyVsDCAz", "DCA_{xy} vs DCA_{z} to PV; DCA_{xy} (cm); DCA_{z} (cm)", kTH2D, {axDCA, axDCA});
        hist.add("Tracks/hDCAxyVsPt", "DCA_{xy} vs pT; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2D, {axMomentum, axDCA});
        hist.add("Tracks/hDCAzVsPt", "DCA_{z} vs pT; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2D, {axMomentum, axDCA});
        hist.add("Tracks/hBeta", "TOF beta vs signed momentum; #it{p}/Z (GeV/#it{c}); #beta", kTH2D, {axSignedMomentum, {1200, 0.0, 1.2}});
        hist.add("Tracks/hTPCNcrossedRows", "TPC number of crossed rows; number of crossed rows TPC; entries", kTH1D, {{110, 60, 170}});
        hist.add("Tracks/hTPCNclusters", "TPC number of clusters; number of clusters TPC; entries", kTH1D, {{200, 0., 200.}});
        hist.add("Tracks/hTPCchi2", "Chi2 / clusters for the TPC track segment; #chi^{2}_{TPC}/n_{TPC}", kTH1D, {{11, -0.5, 10.5}});
        hist.add("Tracks/hITSchi2", "Chi2 / clusters for the ITS track segment; #chi^{2}_{ITS}/n_{ITS}; entries", kTH1D, {{51, -0.5, 50.5}});
        hist.add("Tracks/hTPCchi2vsPt", "Chi2 / clusters for the TPC track segment vs pT; #it{p}_{T} (GeV/#it{c}); #chi^{2}_{TPC}/n_{TPC}; entries", kTH2D, {axMomentum, {11, -0.5, 10.5}});
        hist.add("Tracks/hITSchi2vsPt", "Chi2 / clusters for the ITS track segment vs pT; #it{p}_{T} (GeV/#it{c}); #chi^{2}_{ITS}/n_{ITS}", kTH2D, {axMomentum, {51, -0.5, 50.5}});
        hist.add("Tracks/hPtResolution", "pT resolution; #it{p}^{gen}_{T} (GeV/#it{c}); (#it{p}^{gen}_{T}-#it{p}^{reco}_{T})/#it{p}^{gen}_{T}", kTH2D, {axMomentum, {1000, -10, 10}});
        hist.add("Tracks/hPhiResolution", "phi resolution; #phi^{gen}; (#phi^{gen}-#phi^{reco})/#phi^{gen}", kTH2D, {{500, 0, 6.5}, {100, -0.005, 0.05}});
        hist.add("Tracks/hGenPt", "Generated pT; #it{p}^{gen}_{T} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tracks/hDeltaPt", "#phi^{gen}-#it{p}^{reco}_{T}; (#phi^{gen}-#it{p}^{reco}_{T}) (GeV/#it{c}); entries", kTH1D, {{1000, -10, 10}});

        hist.add("Deuterons/hpx", "Track momentum x; #it{p}_{x} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Deuterons/hpy", "Track momentum y; #it{p}_{y} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Deuterons/hpz", "Track momentum z; #it{p}_{z} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Deuterons/hpt", "Transverse momentum; #it{p}_{T} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Deuterons/hp", "Track momentum; #it{p} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Deuterons/hSign", "Track sign; ; entries", kTH1D, {{2, -2., 2.}});
        hist.add("Deuterons/hEta", "Track eta; #eta; entries", kTH1D, {{500, -1., 1.}});
        hist.add("Deuterons/hPhi", "Track phi; #phi; entries", kTH1D, {{500, 0, 6.5}});
        hist.add("Deuterons/hEtaPhi", "Track eta vs phi; #eta; #phi", kTH2D, {{500, -1., 1.}, {500, 0, 6.5}});
        hist.add("Deuterons/hTPCinnerParam", "TPC mommentum; #it{p}_{TPC} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Deuterons/hTPCSignal", "TPC signal vs TPC momentum; #it{p}_{TPC}/Z (GeV/#it{c}); d#it{E}/d#it{x} (A.U.)", kTH2D, {axSignedMomentum, {1000, 20., 500.}});
        hist.add("Deuterons/hTPCNSigmaPr", "nSigma TPC proton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(p)", kTH2D, {axMomentum, axNsigma});
        hist.add("Deuterons/hTPCNSigmaDe", "nSigma TPC deuteron; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(de)", kTH2D, {axMomentum, axNsigma});
        hist.add("Deuterons/hTPCNSigmaTr", "nSigma TPC triton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(tr)", kTH2D, {axMomentum, axNsigma});
        hist.add("Deuterons/hTPCNSigmaHe", "nSigma TPC helium; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(he)", kTH2D, {axMomentum, axNsigma});
        hist.add("Deuterons/hDCAxy", "DCA_{xy} to PV; DCA_{xy} (cm); entries", kTH1D, {axDCA});
        hist.add("Deuterons/hDCAz", "DCA_{z} to PV; DCA_{z} (cm); entries", kTH1D, {axDCA});
        hist.add("Deuterons/hDCAxyVsDCAz", "DCA_{xy} vs DCA_{z} to PV; DCA_{xy} (cm); DCA_{z} (cm)", kTH2D, {axDCA, axDCA});
        hist.add("Deuterons/hDCAxyVsPt", "DCA_{xy} vs pT; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2D, {axMomentum, axDCA});
        hist.add("Deuterons/hDCAzVsPt", "DCA_{z} vs pT; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2D, {axMomentum, axDCA});
        hist.add("Deuterons/hPtResolution", "pT resolution; #it{p}^{gen}_{T} (GeV/#it{c}); (#it{p}^{gen}_{T}-#it{p}^{reco}_{T})/#it{p}^{gen}_{T}", kTH2D, {axMomentum, {1000, -10, 10}});
        hist.add("Deuterons/hPhiResolution", "phi resolution; #phi^{gen}; (#phi^{gen}-#phi^{reco})/#phi^{gen}", kTH2D, {{500, 0, 6.5}, {100, -0.005, 0.05}});
        hist.add("Deuterons/hGenPt", "Generated pT; #it{p}^{gen}_{T} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Deuterons/hDeltaPt", "#phi^{gen}-#it{p}^{reco}_{T}; (#phi^{gen}-#it{p}^{reco}_{T}) (GeV/#it{c}); entries", kTH1D, {{1000, -10, 10}});

        hist.add("Tritons/hpx", "Track momentum x; #it{p}_{x} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tritons/hpy", "Track momentum y; #it{p}_{y} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tritons/hpz", "Track momentum z; #it{p}_{z} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tritons/hpt", "Transverse momentum; #it{p}_{T} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tritons/hp", "Track momentum; #it{p} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tritons/hSign", "Track sign; ; entries", kTH1D, {{2, -2., 2.}});
        hist.add("Tritons/hEta", "Track eta; #eta; entries", kTH1D, {{500, -1., 1.}});
        hist.add("Tritons/hPhi", "Track phi; #phi; entries", kTH1D, {{500, 0, 6.5}});
        hist.add("Tritons/hEtaPhi", "Track eta vs phi; #eta; #phi", kTH2D, {{500, -1., 1.}, {500, 0, 6.5}});
        hist.add("Tritons/hTPCinnerParam", "TPC mommentum; #it{p}_{TPC} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tritons/hTPCSignal", "TPC signal vs TPC momentum; #it{p}_{TPC}/Z (GeV/#it{c}); d#it{E}/d#it{x} (A.U.)", kTH2D, {axSignedMomentum, {1000, 20., 500.}});
        hist.add("Tritons/hTPCNSigmaPr", "nSigma TPC proton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(p)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tritons/hTPCNSigmaDe", "nSigma TPC deuteron; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(de)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tritons/hTPCNSigmaTr", "nSigma TPC triton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(tr)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tritons/hTPCNSigmaHe", "nSigma TPC helium; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(he)", kTH2D, {axMomentum, axNsigma});
        hist.add("Tritons/hDCAxy", "DCA_{xy} to PV; DCA_{xy} (cm); entries", kTH1D, {axDCA});
        hist.add("Tritons/hDCAz", "DCA_{z} to PV; DCA_{z} (cm); entries", kTH1D, {axDCA});
        hist.add("Tritons/hDCAxyVsDCAz", "DCA_{xy} vs DCA_{z} to PV; DCA_{xy} (cm); DCA_{z} (cm)", kTH2D, {axDCA, axDCA});
        hist.add("Tritons/hDCAxyVsPt", "DCA_{xy} vs pT; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2D, {axMomentum, axDCA});
        hist.add("Tritons/hDCAzVsPt", "DCA_{z} vs pT; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2D, {axMomentum, axDCA});
        hist.add("Tritons/hPtResolution", "pT resolution; #it{p}^{gen}_{T} (GeV/#it{c}); (#it{p}^{gen}_{T}-#it{p}^{reco}_{T})/#it{p}^{gen}_{T}", kTH2D, {axMomentum, {1000, -10, 10}});
        hist.add("Tritons/hPhiResolution", "phi resolution; #phi^{gen}; (#phi^{gen}-#phi^{reco})/#phi^{gen}", kTH2D, {{500, 0, 6.5}, {100, -0.005, 0.05}});
        hist.add("Tritons/hGenPt", "Generated pT; #it{p}^{gen}_{T} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Tritons/hDeltaPt", "#phi^{gen}-#it{p}^{reco}_{T}; (#phi^{gen}-#it{p}^{reco}_{T}) (GeV/#it{c}); entries", kTH1D, {{1000, -10, 10}});

        hist.add("Helium3/hpx", "Track momentum x; #it{p}_{x} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Helium3/hpy", "Track momentum y; #it{p}_{y} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Helium3/hpz", "Track momentum z; #it{p}_{z} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Helium3/hpt", "Transverse momentum; #it{p}_{T} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Helium3/hp", "Track momentum; #it{p} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Helium3/hSign", "Track sign; ; entries", kTH1D, {{2, -2., 2.}});
        hist.add("Helium3/hEta", "Track eta; #eta; entries", kTH1D, {{500, -1., 1.}});
        hist.add("Helium3/hPhi", "Track phi; #phi; entries", kTH1D, {{500, 0, 6.5}});
        hist.add("Helium3/hEtaPhi", "Track eta vs phi; #eta; #phi", kTH2D, {{500, -1., 1.}, {500, 0, 6.5}});
        hist.add("Helium3/hTPCinnerParam", "TPC mommentum; #it{p}_{TPC} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Helium3/hTPCSignal", "TPC signal vs TPC momentum; #it{p}_{TPC}/Z (GeV/#it{c}); d#it{E}/d#it{x} (A.U.)", kTH2D, {axSignedMomentum, {1000, 20., 500.}});
        hist.add("Helium3/hTPCNSigmaPr", "nSigma TPC proton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(p)", kTH2D, {axMomentum, axNsigma});
        hist.add("Helium3/hTPCNSigmaDe", "nSigma TPC deuteron; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(de)", kTH2D, {axMomentum, axNsigma});
        hist.add("Helium3/hTPCNSigmaTr", "nSigma TPC triton; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(tr)", kTH2D, {axMomentum, axNsigma});
        hist.add("Helium3/hTPCNSigmaHe", "nSigma TPC helium; #it{p}_{T} (GeV/#it{c}); n#it{#sigma}_{TPC}(he)", kTH2D, {axMomentum, axNsigma});
        hist.add("Helium3/hDCAxy", "DCA_{xy} to PV; DCA_{xy} (cm); entries", kTH1D, {axDCA});
        hist.add("Helium3/hDCAz", "DCA_{z} to PV; DCA_{z} (cm); entries", kTH1D, {axDCA});
        hist.add("Helium3/hDCAxyVsDCAz", "DCA_{xy} vs DCA_{z} to PV; DCA_{xy} (cm); DCA_{z} (cm)", kTH2D, {axDCA, axDCA});
        hist.add("Helium3/hDCAxyVsPt", "DCA_{xy} vs pT; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2D, {axMomentum, axDCA});
        hist.add("Helium3/hDCAzVsPt", "DCA_{z} vs pT; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2D, {axMomentum, axDCA});
        hist.add("Helium3/hPtResolution", "pT resolution; #it{p}^{gen}_{T} (GeV/#it{c}); (#it{p}^{gen}_{T}-#it{p}^{reco}_{T})/#it{p}^{gen}_{T}", kTH2D, {axMomentum, {1000, -10, 10}});
        hist.add("Helium3/hPhiResolution", "phi resolution; #phi^{gen}; (#phi^{gen}-#phi^{reco})/#phi^{gen}", kTH2D, {{500, 0, 6.5}, {100, -0.005, 0.05}});
        hist.add("Helium3/hGenPt", "Generated pT; #it{p}^{gen}_{T} (GeV/#it{c}); entries", kTH1D, {axMomentum});
        hist.add("Helium3/hDeltaPt", "#phi^{gen}-#it{p}^{reco}_{T}; (#phi^{gen}-#it{p}^{reco}_{T}) (GeV/#it{c}); entries", kTH1D, {{1000, -10, 10}});
        
    } /// end init function

    /// TODO: event selections: sel8, nContr>1, Vtx_z<10
    /// event selection function (to be used later on)
    template <typename CollisionType>
    bool isSelectedColl(const CollisionType& collision) { /// here I have to pass a single collision
        if (isRun3) {
            if (!collision.sel8()) return false;
        }
        else {
            if (!collision.sel7()) return false;
        }
        return true;
    }

    template <typename TrackType>
    bool isTrackInAcceptance(const TrackType& track) { /// here I have to pass a single track
        if (abs(track.eta()) > eta) return false;
        return true;
    }

    /// function to fill event histograms
    template <typename CollisionType>
    void fillCollisionHistos(const CollisionType& collision) {
        hist.fill(HIST("Events/hVertexX"), collision.posX());
        hist.fill(HIST("Events/hVertexY"), collision.posY());
        hist.fill(HIST("Events/hVertexXY"), collision.posX(), collision.posY());
        hist.fill(HIST("Events/hVertexZ"), collision.posZ());
        hist.fill(HIST("Events/hVertexCovXX"), collision.covXX());
        hist.fill(HIST("Events/hVertexCovYY"), collision.covYY());
        hist.fill(HIST("Events/hVertexCovZZ"), collision.covZZ());
        hist.fill(HIST("Events/hVertexCovXY"), collision.covXY());
        hist.fill(HIST("Events/hVertexCovXZ"), collision.covXZ());
        hist.fill(HIST("Events/hVertexCovYZ"), collision.covYZ());
        hist.fill(HIST("Events/hVertexNcontrib"), collision.multNTracksPV());
    }

    /// function to fill event histograms
    template <typename TrackType>
    void fillTrackHistos(const TrackType& track) {
        hist.fill(HIST("Tracks/hpx"), track.px());
        hist.fill(HIST("Tracks/hpy"), track.py());
        hist.fill(HIST("Tracks/hpz"), track.pz());
        hist.fill(HIST("Tracks/hpt"), track.pt());
        hist.fill(HIST("Tracks/hp"), track.p());
        hist.fill(HIST("Tracks/hSign"), track.sign());
        hist.fill(HIST("Tracks/hEta"), track.eta());
        hist.fill(HIST("Tracks/hPhi"), track.phi());
        hist.fill(HIST("Tracks/hEtaPhi"), track.eta(), track.phi());
        hist.fill(HIST("Tracks/hTPCinnerParam"), track.tpcInnerParam());
        hist.fill(HIST("Tracks/hTPCSignal"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
        hist.fill(HIST("Tracks/hTPCNSigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        hist.fill(HIST("Tracks/hTPCNSigmaDe"), track.tpcInnerParam(), track.tpcNSigmaDe());
        hist.fill(HIST("Tracks/hTPCNSigmaTr"), track.tpcInnerParam(), track.tpcNSigmaTr());
        hist.fill(HIST("Tracks/hTPCNSigmaHe"), track.tpcInnerParam(), track.tpcNSigmaHe());
        hist.fill(HIST("Tracks/hTOFNSigmaPr"), track.p(), track.tofNSigmaPr());
        hist.fill(HIST("Tracks/hTOFNSigmaDe"), track.p(), track.tofNSigmaDe());
        hist.fill(HIST("Tracks/hTOFNSigmaTr"), track.p(), track.tofNSigmaTr());
        hist.fill(HIST("Tracks/hTOFNSigmaHe"), track.p(), track.tofNSigmaHe());
        hist.fill(HIST("Tracks/hDCAxy"), track.dcaXY());
        hist.fill(HIST("Tracks/hDCAz"), track.dcaZ());
        hist.fill(HIST("Tracks/hDCAxyVsDCAz"), track.dcaXY(), track.dcaZ());
        hist.fill(HIST("Tracks/hDCAxyVsPt"), track.pt(), track.dcaXY());
        hist.fill(HIST("Tracks/hDCAzVsPt"), track.pt(), track.dcaZ());
        hist.fill(HIST("Tracks/hBeta"), track.p()/track.sign(), track.beta());
        hist.fill(HIST("Tracks/hTPCNcrossedRows"), track.tpcNClsCrossedRows());
        hist.fill(HIST("Tracks/hTPCNclusters"), track.tpcNClsFound());
        hist.fill(HIST("Tracks/hTPCchi2"), track.tpcChi2NCl());
        hist.fill(HIST("Tracks/hITSchi2"), track.itsChi2NCl());
        hist.fill(HIST("Tracks/hTPCchi2vsPt"), track.pt(), track.tpcChi2NCl());
        hist.fill(HIST("Tracks/hITSchi2vsPt"), track.pt(), track.itsChi2NCl());
    }

    template <typename TrackType>
    void fillDeuteronHistos(const TrackType& track) {
        hist.fill(HIST("Deuterons/hpx"), track.px());
        hist.fill(HIST("Deuterons/hpy"), track.py());
        hist.fill(HIST("Deuterons/hpz"), track.pz());
        hist.fill(HIST("Deuterons/hpt"), track.pt());
        hist.fill(HIST("Deuterons/hp"), track.p());
        hist.fill(HIST("Deuterons/hSign"), track.sign());
        hist.fill(HIST("Deuterons/hEta"), track.eta());
        hist.fill(HIST("Deuterons/hPhi"), track.phi());
        hist.fill(HIST("Deuterons/hEtaPhi"), track.eta(), track.phi());
        hist.fill(HIST("Deuterons/hTPCinnerParam"), track.tpcInnerParam());
        hist.fill(HIST("Deuterons/hTPCSignal"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
        hist.fill(HIST("Deuterons/hTPCNSigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        hist.fill(HIST("Deuterons/hTPCNSigmaDe"), track.tpcInnerParam(), track.tpcNSigmaDe());
        hist.fill(HIST("Deuterons/hTPCNSigmaTr"), track.tpcInnerParam(), track.tpcNSigmaTr());
        hist.fill(HIST("Deuterons/hTPCNSigmaHe"), track.tpcInnerParam(), track.tpcNSigmaHe());
        hist.fill(HIST("Deuterons/hDCAxy"), track.dcaXY());
        hist.fill(HIST("Deuterons/hDCAz"), track.dcaZ());
        hist.fill(HIST("Deuterons/hDCAxyVsDCAz"), track.dcaXY(), track.dcaZ());
        hist.fill(HIST("Deuterons/hDCAxyVsPt"), track.pt(), track.dcaXY());
        hist.fill(HIST("Deuterons/hDCAzVsPt"), track.pt(), track.dcaZ());
    }

    template <typename TrackType>
    void fillTritonHistos(const TrackType& track) {
        hist.fill(HIST("Tritons/hpx"), track.px());
        hist.fill(HIST("Tritons/hpy"), track.py());
        hist.fill(HIST("Tritons/hpz"), track.pz());
        hist.fill(HIST("Tritons/hpt"), track.pt());
        hist.fill(HIST("Tritons/hp"), track.p());
        hist.fill(HIST("Tritons/hSign"), track.sign());
        hist.fill(HIST("Tritons/hEta"), track.eta());
        hist.fill(HIST("Tritons/hPhi"), track.phi());
        hist.fill(HIST("Tritons/hEtaPhi"), track.eta(), track.phi());
        hist.fill(HIST("Tritons/hTPCinnerParam"), track.tpcInnerParam());
        hist.fill(HIST("Tritons/hTPCSignal"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
        hist.fill(HIST("Tritons/hTPCNSigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        hist.fill(HIST("Tritons/hTPCNSigmaDe"), track.tpcInnerParam(), track.tpcNSigmaDe());
        hist.fill(HIST("Tritons/hTPCNSigmaTr"), track.tpcInnerParam(), track.tpcNSigmaTr());
        hist.fill(HIST("Tritons/hTPCNSigmaHe"), track.tpcInnerParam(), track.tpcNSigmaHe());
        hist.fill(HIST("Tritons/hDCAxy"), track.dcaXY());
        hist.fill(HIST("Tritons/hDCAz"), track.dcaZ());
        hist.fill(HIST("Tritons/hDCAxyVsDCAz"), track.dcaXY(), track.dcaZ());
        hist.fill(HIST("Tritons/hDCAxyVsPt"), track.pt(), track.dcaXY());
        hist.fill(HIST("Tritons/hDCAzVsPt"), track.pt(), track.dcaZ());
    }

    template <typename TrackType>
    void fillHeliumHistos(const TrackType& track) {
        hist.fill(HIST("Helium3/hpx"), track.px());
        hist.fill(HIST("Helium3/hpy"), track.py());
        hist.fill(HIST("Helium3/hpz"), track.pz());
        hist.fill(HIST("Helium3/hpt"), track.pt());
        hist.fill(HIST("Helium3/hp"), track.p());
        hist.fill(HIST("Helium3/hSign"), track.sign());
        hist.fill(HIST("Helium3/hEta"), track.eta());
        hist.fill(HIST("Helium3/hPhi"), track.phi());
        hist.fill(HIST("Helium3/hEtaPhi"), track.eta(), track.phi());
        hist.fill(HIST("Helium3/hTPCinnerParam"), track.tpcInnerParam());
        hist.fill(HIST("Helium3/hTPCSignal"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
        hist.fill(HIST("Helium3/hTPCNSigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        hist.fill(HIST("Helium3/hTPCNSigmaDe"), track.tpcInnerParam(), track.tpcNSigmaDe());
        hist.fill(HIST("Helium3/hTPCNSigmaTr"), track.tpcInnerParam(), track.tpcNSigmaTr());
        hist.fill(HIST("Helium3/hTPCNSigmaHe"), track.tpcInnerParam(), track.tpcNSigmaHe());
        hist.fill(HIST("Helium3/hDCAxy"), track.dcaXY());
        hist.fill(HIST("Helium3/hDCAz"), track.dcaZ());
        hist.fill(HIST("Helium3/hDCAxyVsDCAz"), track.dcaXY(), track.dcaZ());
        hist.fill(HIST("Helium3/hDCAxyVsPt"), track.pt(), track.dcaXY());
        hist.fill(HIST("Helium3/hDCAzVsPt"), track.pt(), track.dcaZ());
    }

    /// function to fill the event tree
    template <typename CollisionType>
    void fillCollisionTable(const CollisionType& collision) {
        /// filling event properties
        tableCollisions(collision.bc().runNumber(),
                        collision.posX(),
                        collision.posY(),
                        collision.posZ(),
                        collision.covXX(),
                        collision.covYY(),
                        collision.covZZ(),
                        collision.multNTracksPV(),
                        collision.chi2(),
                        collision.sel8());
    }

    /// collision and track tables
    using CollisionTable = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
    using TrackTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTOFbeta, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe>;

    /// process function data --> only MC process function further updated, this is obsolete for now!!
    void processData(CollisionTable::iterator const& collision, 
                   soa::Filtered<TrackTable> const& tracks, /// Filter: Filterbit Global track, ...
                   aod::BCsWithTimestamps const&) {
        auto bc = collision.bc_as<aod::BCsWithTimestamps>();
        if (runNumber!=bc.runNumber()) {
            initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag:ccdbPathGrp, lut, isRun3);
            magneticField = o2::base::Propagator::Instance()->getNominalBz(); /// what is happening here??
        }
        /// fill first bin of event selection hist
        hist.fill(HIST("Events/hEventSels"), 1.f);
        /// event selection
        if (eventSel) {
            /// select collisions according to sel8 trigger
            if (!isSelectedColl(collision)) return;
            hist.fill(HIST("Events/hEventSels"), 2.f);
            /// select collisions in acceptance
            if (!(abs(collision.posZ()) < 10)) return;
            hist.fill(HIST("Events/hEventSels"), 3.f);
            /// select PVs with more than one contributor
            if (!(collision.multNTracksPV() > 1)) return;
            hist.fill(HIST("Events/hEventSels"), 4.f);
        }

        /// fill event histograms and tree
        fillCollisionHistos(collision);
        if (writeCollTable) fillCollisionTable(collision);

        /// track loop to get event multiplicity
        int nTracks = 0.;
        for (auto& track : tracks) {
            /// apply track selection
            if (!isTrackInAcceptance(track)) return; /// eta filter
            if (track.tpcNClsFound() < tpcNClsFound) return;
            fillTrackHistos(track);
            nTracks++;
        }
        hist.fill(HIST("Events/hMultiplicity"), nTracks);
    }
    PROCESS_SWITCH(NucleiRecoQA, processData, "process Data", true);

    /// process function MC
    void processMC(soa::Join<CollisionTable, aod::McCollisionLabels>::iterator const& collision, 
                   soa::Filtered<soa::Join<TrackTable, aod::McTrackLabels>> const& tracks, /// Filter: Filterbit Global track, ...
                   aod::BCsWithTimestamps const&, aod::McCollisions const& mcCollisions, 
                   aod::McParticles const& particlesMC) {
        auto bc = collision.bc_as<aod::BCsWithTimestamps>();
        if (runNumber!=bc.runNumber()) {
            initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag:ccdbPathGrp, lut, isRun3);
            magneticField = o2::base::Propagator::Instance()->getNominalBz(); /// what is happening here??
        }
        /// fill first bin of event selection hist
        hist.fill(HIST("Events/hEventSels"), 1.f);
        /// remove collisions without a MC collision
        if (!collision.has_mcCollision()) return;
        hist.fill(HIST("Events/hEventSels"), 2.f);
        /// event selection
        if (eventSel) {
            /// select collisions according to sel8 trigger
            if (!isSelectedColl(collision)) return;
            hist.fill(HIST("Events/hEventSels"), 3.f);
            /// select collisions in acceptance
            if (!(abs(collision.posZ()) < 10)) return;
            hist.fill(HIST("Events/hEventSels"), 4.f);
            /// select PVs with more than one contributor
            if (!(collision.multNTracksPV() > 1)) return;
            hist.fill(HIST("Events/hEventSels"), 5.f);
        }

        /// fill event histograms and tree
        fillCollisionHistos(collision);
        if (writeCollTable) fillCollisionTable(collision);

        /// track loop to get event multiplicity
        int nTracks = 0.;
        int nDeuterons = 0.;
        int nTritons = 0.;
        int nHelium3 = 0.;
        for (auto& track : tracks) {
            /// apply track selection
            if (!isTrackInAcceptance(track)) continue; /// eta filter
            /// fill histograms
            fillTrackHistos(track);
            nTracks++;

            /// select deuterons
            if (track.tpcNSigmaDe() > 0 && track.tpcNSigmaDe() < 4 && track.pt() > 0.3) {
                fillDeuteronHistos(track);
                nDeuterons++;
                // access MC truth information with mcCollision() and mcParticle() methods
                if (track.mcParticleId() >= -1 && track.mcParticleId() <= particlesMC.size()) {
                    auto GenPt = track.mcParticle().pt();
                    auto RecoPt = track.pt();
                    auto deltaPt = GenPt - RecoPt;
                    auto deltaPhi = track.mcParticle().phi() - track.phi();
                    hist.fill(HIST("Deuterons/hPtResolution"), GenPt, deltaPt);
                    hist.fill(HIST("Deuterons/hGenPt"), GenPt);
                    hist.fill(HIST("Deuterons/hDeltaPt"), deltaPt);
                    hist.fill(HIST("Deuterons/hPhiResolution"), track.mcParticle().phi(), deltaPhi);
                };
            };

            /// select tritons
            if (track.tpcNSigmaTr() > -0.5 && track.tpcNSigmaTr() < 4 && track.pt() > 0.9) {
                fillTritonHistos(track);
                nTritons++;
                // access MC truth information with mcCollision() and mcParticle() methods
                if (track.mcParticleId() >= -1 && track.mcParticleId() <= particlesMC.size()) {
                    auto GenPt = track.mcParticle().pt();
                    auto RecoPt = track.pt();
                    auto deltaPt = GenPt - RecoPt;
                    auto deltaPhi = track.mcParticle().phi() - track.phi();
                    hist.fill(HIST("Tritons/hPtResolution"), GenPt, deltaPt);
                    hist.fill(HIST("Tritons/hGenPt"), GenPt);
                    hist.fill(HIST("Tritons/hDeltaPt"), deltaPt);
                    hist.fill(HIST("Tritons/hPhiResolution"), track.mcParticle().phi(), deltaPhi);
                };
            };

            /// select helium3 with high dE/dx
            if (track.pt() > ptMin && track.tpcSignal() > dEdxMin) {
                fillHeliumHistos(track);
                nHelium3++;
                // access MC truth information with mcCollision() and mcParticle() methods
                if (track.mcParticleId() >= -1 && track.mcParticleId() <= particlesMC.size()) {
                    auto GenPt = track.mcParticle().pt();
                    auto RecoPt = track.pt();
                    auto deltaPt = GenPt - RecoPt;
                    auto deltaPhi = track.mcParticle().phi() - track.phi();
                    hist.fill(HIST("Helium3/hPtResolution"), GenPt, deltaPt);
                    hist.fill(HIST("Helium3/hGenPt"), GenPt);
                    hist.fill(HIST("Helium3/hDeltaPt"), deltaPt);
                    hist.fill(HIST("Helium3/hPhiResolution"), track.mcParticle().phi(), deltaPhi);
                };
            };
        }
        hist.fill(HIST("Events/hMultiplicity"), nTracks);
        hist.fill(HIST("Events/hMultiplicityDeuterons"), nDeuterons);
        hist.fill(HIST("Events/hMultiplicityTritons"), nTritons);
        hist.fill(HIST("Events/hMultiplicityHelium3"), nHelium3);
    }
    PROCESS_SWITCH(NucleiRecoQA, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleiRecoQA>(cfgc)
  };
}