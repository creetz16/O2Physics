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
// This code runs loop over v0 photons for PCM QC.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include "TString.h"
#include "THashList.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzRotation.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

using MyV0Photons = soa::Join<aod::V0Photons, aod::V0RecalculationAndKF>;
using MyV0Photon = MyV0Photons::iterator;

struct PCMQCMC {
  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "analysis,qc,nocut", "Comma separated list of v0 photon cuts"};

  std::vector<V0PhotonCut> fPCMCuts;

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputTrack{"Track"};
  OutputObj<THashList> fOutputV0{"V0"};
  OutputObj<THashList> fOutputGen{"Generated"};
  THashList* fMainList = new THashList();

  // static constexpr std::string_view ambtracktypes[2] = {"NonAmb", "Amb"};
  void addhistograms()
  {

    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));
    o2::aod::emphotonhistograms::DefineHistograms(list_ev, "Event");

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Track");
    THashList* list_tr = reinterpret_cast<THashList*>(fMainList->FindObject("Track"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "V0");
    THashList* list_v0 = reinterpret_cast<THashList*>(fMainList->FindObject("V0"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::emphotonhistograms::DefineHistograms(list_gen, "Generated", "ConversionStudy");

    for (const auto& cut : fPCMCuts) {
      const char* cutname = cut.GetName();
      o2::aod::emphotonhistograms::AddHistClass(list_tr, cutname);
      o2::aod::emphotonhistograms::AddHistClass(list_v0, cutname);
    }

    // for single tracks
    for (auto& cut : fPCMCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("Track")->FindObject(cutname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list, "Track");
    }

    // for V0s
    for (auto& cut : fPCMCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("V0")->FindObject(cutname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list, "V0");
    }
  }

  void DefineCuts()
  {
    TString cutNamesStr = fConfigPCMCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPCMCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of PCM cuts = %d", fPCMCuts.size());
  }

  void init(InitContext& context)
  {
    DefineCuts();
    addhistograms(); // please call this after DefinCuts();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputTrack.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Track")));
    fOutputV0.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0")));
    fOutputGen.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Generated")));
  }

  template <typename T>
  void fillHistosLeg(const T& leg, const char* cutname)
  {
    reinterpret_cast<TH1F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hPt"))->Fill(leg.pt());
    reinterpret_cast<TH2F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hEtaPhi"))->Fill(leg.phi(), leg.eta());
    reinterpret_cast<TH2F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hDCAxyz"))->Fill(leg.dcaXY(), leg.dcaZ());
    reinterpret_cast<TH1F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hNclsTPC"))->Fill(leg.tpcNClsFound());
    reinterpret_cast<TH1F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hNclsITS"))->Fill(leg.itsNCls());
    reinterpret_cast<TH1F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hNcrTPC"))->Fill(leg.tpcNClsCrossedRows());
    reinterpret_cast<TH1F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hTPCNcr2Nf"))->Fill(leg.tpcCrossedRowsOverFindableCls());
    reinterpret_cast<TH1F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hChi2TPC"))->Fill(leg.tpcChi2NCl());
    reinterpret_cast<TH1F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hChi2ITS"))->Fill(leg.itsChi2NCl());
    reinterpret_cast<TH2F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hTPCdEdx"))->Fill(leg.tpcInnerParam(), leg.tpcSignal());
    reinterpret_cast<TH2F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hTPCNsigmaEl"))->Fill(leg.tpcInnerParam(), leg.tpcNSigmaEl());
    reinterpret_cast<TH2F*>(fMainList->FindObject("Track")->FindObject(cutname)->FindObject("hTPCNsigmaPi"))->Fill(leg.tpcInnerParam(), leg.tpcNSigmaPi());
  }

  template <typename T>
  void fillHistosV0(const T& v0, const char* cutname)
  {
    reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hPt"))->Fill(v0.pt());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hEtaPhi"))->Fill(v0.phi(), v0.eta());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hRadius"))->Fill(v0.vz(), v0.v0radius());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hRadius_recalc"))->Fill(v0.recalculatedVtxZ(), v0.recalculatedVtxR());
    reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hCosPA"))->Fill(abs(v0.cospa()));
    reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hPCA"))->Fill(v0.pca());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hAPplot"))->Fill(v0.alpha(), v0.qtarm());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hMassGamma"))->Fill(v0.v0radius(), v0.mGamma());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hMassGamma_recalc"))->Fill(v0.recalculatedVtxR(), v0.mGamma());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hGammaPsiPair"))->Fill(v0.psipair(), v0.mGamma());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hGammaRxy"))->Fill(v0.vx(), v0.vy());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hGammaRxy_recalc"))->Fill(v0.recalculatedVtxX(), v0.recalculatedVtxY());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hKFChi2vsR_recalc"))->Fill(v0.recalculatedVtxR(), v0.chiSquareNDF());
    reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cutname)->FindObject("hKFChi2vsZ_recalc"))->Fill(v0.recalculatedVtxZ(), v0.chiSquareNDF());
  }

  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int FindCommonMotherFrom2Prongs(TMCParticle1 const& p1, TMCParticle2 const& p2, const int expected_pdg1, const int expected_pdg2, const int expected_mother_pdg, TMCParticles const& mcparticles)
  {
    if (p1.globalIndex() == p2.globalIndex())
      return -1; // mc particle p1 and p2 is identical. reject.

    if (p1.pdgCode() != expected_pdg1)
      return -1;
    if (p2.pdgCode() != expected_pdg2)
      return -1;

    if (!p1.has_mothers())
      return -1;
    if (!p2.has_mothers())
      return -1;

    // LOGF(info,"original motherid1 = %d , motherid2 = %d", p1.mothersIds()[0], p2.mothersIds()[0]);

    int motherid1 = p1.mothersIds()[0];
    auto mother1 = mcparticles.iteratorAt(motherid1);
    int mother1_pdg = mother1.pdgCode();

    int motherid2 = p2.mothersIds()[0];
    auto mother2 = mcparticles.iteratorAt(motherid2);
    int mother2_pdg = mother2.pdgCode();

    // LOGF(info,"motherid1 = %d , motherid2 = %d", motherid1, motherid2);

    if (motherid1 != motherid2)
      return -1;
    if (mother1_pdg != mother2_pdg)
      return -1;
    if (mother1_pdg != expected_mother_pdg)
      return -1;
    return motherid1;
  }

  Preslice<MyV0Photons> perCollision = aod::v0photon::collisionId;

  using MyMCV0Legs = soa::Join<aod::V0Legs, aod::EMMCParticleLabels>;
  void processQCMC(soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels> const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const& v0legs, aod::EMMCParticles const& mcparticles)
  {
    for (auto& collision : collisions) {
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hZvtx_before"))->Fill(collision.posZ());
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(1.0);
      if (!collision.sel8()) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(2.0);

      if (collision.numContrib() < 0.5) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(3.0);

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(4.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hZvtx_after"))->Fill(collision.posZ());
      auto V0Photons_coll = v0photons.sliceBy(perCollision, collision.collisionId());

      for (const auto& cut : fPCMCuts) {
        int ng = 0;
        for (auto& g : V0Photons_coll) {
          auto pos = g.posTrack_as<MyMCV0Legs>();
          auto ele = g.negTrack_as<MyMCV0Legs>();
          auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
          auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

          int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
          if (photonid < 0) { // check swap, true electron is reconstructed as positron and vice versa.
            photonid = FindCommonMotherFrom2Prongs(posmc, elemc, 11, -11, 22, mcparticles);
          }

          if (photonid < 0) {
            continue;
          }

          if (cut.IsSelected<MyMCV0Legs>(g)) {
            fillHistosV0(g, cut.GetName());
            ng++;
            for (auto& leg : {pos, ele}) {
              fillHistosLeg(leg, cut.GetName());
            }
          }
        } // end of v0 loop
        reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hNgamma"))->Fill(ng);
      } // end of cut loop
    }   // end of collision loop
  }     // end of process

  Preslice<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emreducedmceventId;
  void processGen(soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels> const& collisions, aod::EMReducedMCEvents const& mccollisions, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : collisions) {
      auto mccollision = collision.emreducedmcevent();
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(1.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_before"))->Fill(mccollision.mcPosZ());
      if (!collision.sel8()) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(2.0);

      if (collision.numContrib() < 0.5) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(3.0);

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(4.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_after"))->Fill(mccollision.mcPosZ());

      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (auto& mctrack : mctracks_coll) {

        if (IsEleFromPC(mctrack, mcparticles) > 0) {
          float rxy = sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2));
          reinterpret_cast<TH2F*>(fMainList->FindObject("Generated")->FindObject("hGammaRZ"))->Fill(mctrack.vz(), rxy);

          if (abs(mctrack.eta()) > 0.9) {
            continue;
          }
          reinterpret_cast<TH2F*>(fMainList->FindObject("Generated")->FindObject("hGammaRxy"))->Fill(mctrack.vx(), mctrack.vy());
          if (rxy < 6) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy0_6cm"))->Fill(mctrack.phi());
          } else if (rxy < 10) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy6_10cm"))->Fill(mctrack.phi());
          } else if (rxy < 20) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy10_20cm"))->Fill(mctrack.phi());
          } else if (rxy < 30) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy20_30cm"))->Fill(mctrack.phi());
          } else if (rxy < 40) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy30_40cm"))->Fill(mctrack.phi());
          } else if (rxy < 50) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy40_50cm"))->Fill(mctrack.phi());
          } else if (rxy < 60) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy50_60cm"))->Fill(mctrack.phi());
          } else if (rxy < 70) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy60_70cm"))->Fill(mctrack.phi());
          } else if (rxy < 80) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy70_80cm"))->Fill(mctrack.phi());
          } else if (rxy < 90) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hConvPhi_Rxy80_90cm"))->Fill(mctrack.phi());
          }
        }
      }
    }
  }

  void processDummy(aod::EMReducedEvents::iterator const& collision)
  {
    // do nothing
  }

  PROCESS_SWITCH(PCMQCMC, processQCMC, "run PCM QC in MC", true);
  PROCESS_SWITCH(PCMQCMC, processGen, "run generated information", false);
  PROCESS_SWITCH(PCMQCMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PCMQCMC>(cfgc, TaskName{"pcm-qc-mc"})};
}
