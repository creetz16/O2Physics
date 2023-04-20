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

/// \file NucleiRecoQA.h

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_NUCLEIRECOQATABLES_H_
#define PWGLF_DATAMODEL_NUCLEIRECOQATABLES_H_

using namespace o2;

namespace o2::aod
{
namespace eventtrack
{
/// (Name, Getter, Type)
/// events
DECLARE_SOA_COLUMN(IsEventRejected, isEventRejected, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(XVtx, xVtx, float);
DECLARE_SOA_COLUMN(YVtx, yVtx, float);
DECLARE_SOA_COLUMN(ZVtx, zVtx, float);
DECLARE_SOA_COLUMN(CovXX, covXX, float);
DECLARE_SOA_COLUMN(CovYY, covYY, float);
DECLARE_SOA_COLUMN(CovZZ, covZZ, float);
DECLARE_SOA_COLUMN(NContrib, nContrib, float);
DECLARE_SOA_COLUMN(Chi2, chi2, float);
/// tracks
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float);
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float);
DECLARE_SOA_COLUMN(TPCNSigmaTr, tpcNSigmaTr, float);
DECLARE_SOA_COLUMN(TPCNSigmaHe, tpcNSigmaHe, float);
DECLARE_SOA_COLUMN(TPCNSigmaAl, tpcNSigmaAl, float);
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNSigmaPr, float);
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float);
DECLARE_SOA_COLUMN(TOFNSigmaTr, tofNSigmaTr, float);
DECLARE_SOA_COLUMN(TOFNSigmaHe, tofNSigmaHe, float);
DECLARE_SOA_COLUMN(TOFNSigmaAl, tofNSigmaAl, float);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);
DECLARE_SOA_COLUMN(Beta, beta, float);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int16_t);
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);
DECLARE_SOA_COLUMN(TPCChi2Ncl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);
DECLARE_SOA_COLUMN(GenPt, genPt, float);           // Momentum of the candidate (x direction)
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);         // Momentum of the candidate (y direction)
DECLARE_SOA_COLUMN(GenEta, genEta, float);         // Momentum of the candidate (z direction)
} // namespace eventtrack

DECLARE_SOA_TABLE(TableCollisions, "AOD", "TableCollisions",
                  eventtrack::RunNumber,
                  eventtrack::XVtx,
                  eventtrack::YVtx,
                  eventtrack::ZVtx,
                  eventtrack::CovXX,
                  eventtrack::CovYY,
                  eventtrack::CovZZ,
                  eventtrack::NContrib,
                  eventtrack::Chi2,
                  eventtrack::IsEventRejected);

DECLARE_SOA_TABLE(TableTracks, "AOD", "TableTracks",
                  eventtrack::Px,
                  eventtrack::Py,
                  eventtrack::Pz,
                  eventtrack::Pt,
                  eventtrack::P,
                  eventtrack::Sign,
                  eventtrack::Eta,
                  eventtrack::Phi,
                  eventtrack::TPCInnerParam,
                  eventtrack::TPCSignal,
                  eventtrack::TPCNSigmaPr,
                  eventtrack::TPCNSigmaDe,
                  eventtrack::TPCNSigmaTr,
                  eventtrack::TPCNSigmaHe,
                  eventtrack::TPCNSigmaAl,
                  eventtrack::TOFNSigmaPr,
                  eventtrack::TOFNSigmaDe,
                  eventtrack::TOFNSigmaTr,
                  eventtrack::TOFNSigmaHe,
                  eventtrack::TOFNSigmaAl,
                  eventtrack::HasTOF,
                  eventtrack::HasTRD,
                  eventtrack::DcaXY,
                  eventtrack::DcaZ,
                  eventtrack::Beta,
                  eventtrack::TPCNClsCrossedRows,
                  eventtrack::TPCCrossedRowsOverFindableCls,
                  eventtrack::TPCNClsFound,
                  eventtrack::GenPt,
                  eventtrack::GenPhi,
                  eventtrack::GenEta);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_NUCLEIRECOQATABLES_H_
