#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "FCTTrack.h"
#include "Field/MagneticField.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "ITSMFTSimulation/Hit.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTree.h"
#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TProfile.h>
#include <TStyle.h>

#endif

// FCTTools

#include "fcttools/HistosHelpers.C"
#include "fcttools/MagField.C"

using o2::MCTrackT;
using o2::fct::FCTTrack;
using o2::itsmft::Hit;
using eventFoundTracks = std::vector<bool>;
using std::vector;

bool DEBUG_VERBOSE = false;
bool EXPORT_HISTOS_IMAGES = false;
Int_t minHitsPerTrack = 10;

bool InnerBorder(float_t tanl)
{ // -3.6
  auto abstanl = std::abs(tanl);
  return (abstanl > 18.1 && abstanl < 18.2855);
}

bool InnerRegion(float_t tanl)
{ // 3.5 < |eta| < 3.6
  auto abstanl = std::abs(tanl);
  return (abstanl > 16.542646 && abstanl < 18.2855);
}

bool OuterBorder(float_t tanl)
{ // -2.8
  auto abstanl = std::abs(tanl);
  return (abstanl > 8.1919179 && abstanl < 8.2748525);
}

bool OuterRegion(float_t tanl)
{ // 2.8 < |eta| < 2.9
  auto abstanl = std::abs(tanl);
  return (abstanl > 8.1919179 && abstanl < 9.0595683);
}

bool pt_1(float_t pt)
{
  return ((pt > 0.9) && (pt < 1.1));
}

bool pt_4(float_t pt)
{
  return ((pt > 3.9) && (pt < 4.1));
}

//_________________________________________________________________________________________________
int FCTTrackerChecker(const Char_t* trkFile = "fcttracks.root",
                      const Char_t* o2sim_KineFile = "o2sim_Kine.root",
                      const Char_t* HitsFCTFile = "o2sim_HitsFCT.root")
{

  // Histos parameters
  Double_t pMin = 0.0;
  Double_t pMax = 100.0;
  Double_t deltaetaMin = -.1;
  Double_t deltaetaMax = +.1;
  Double_t etaMin = -5.0;
  Double_t etaMax = +5.0;
  Double_t deltaphiMin = -.2;
  Double_t deltaphiMax = .2;
  Double_t deltatanlMin = -2.0;
  Double_t deltatanlMax = 2.0;

  // Seed configuration
  std::string seed_cfg{trkFile};
  std::string trk_start{"fcttracks_"};
  std::string trk_ext{".root"};
  std::string trk_trk{"fcttracks"};
  if (seed_cfg.find(trk_start) < seed_cfg.length())
    seed_cfg.replace(seed_cfg.find(trk_start), trk_start.length(), "");
  if (seed_cfg.find(trk_ext) < seed_cfg.length())
    seed_cfg.replace(seed_cfg.find(trk_ext), trk_ext.length(), "");
  if (seed_cfg.find(trk_trk) < seed_cfg.length())
    seed_cfg.replace(seed_cfg.find(trk_trk), trk_trk.length(), "");
  std::cout << seed_cfg << std::endl;

  // histos
  gStyle->SetOptStat("emr");
  gStyle->SetStatW(.28);
  gStyle->SetStatH(.26);
  gStyle->SetPalette(1, 0);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(4);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(3);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(3);
  gStyle->SetLabelSize(0.06, "xyz");
  gStyle->SetLabelOffset(0.01, "y");
  gStyle->SetTitleSize(0.06, "xyz");
  gStyle->SetTitleSize(0.08, "o");
  gStyle->SetTitleOffset(.95, "Y");
  gStyle->SetTitleFillColor(10);
  gStyle->SetStatColor(10);

  enum TH3HistosCodes {
    kFCTTrackDeltaXVertexPtEta,
    kFCTTrackDeltaYVertexPtEta,
    kFCTTrackPtResolutionPtEta,
    kFCTTrackInvPtResolutionPtEta,
    kFCTTrackInvQPtResolutionPtEta,
    kFCTTrackInvQPtPullPtEta
  };

  std::map<int, const char*> TH3Names{
    {kFCTTrackDeltaXVertexPtEta, "FCTTrackDeltaXVertexPtEta"},
    {kFCTTrackDeltaYVertexPtEta, "FCTTrackDeltaYVertexPtEta"},
    {kFCTTrackPtResolutionPtEta, "FCTTrackPtResolutionPtEta"},
    {kFCTTrackInvPtResolutionPtEta, "FCTTrackInvPtResolutionPtEta"},
    {kFCTTrackInvQPtResolutionPtEta, "FCTTrackInvQPtResolutionPtEta"},
    {kFCTTrackInvQPtPullPtEta, "FCTTrackInvQPtPullPtEta"}};
  //
  std::map<int, const char*> TH3Titles{
    {kFCTTrackDeltaXVertexPtEta, "FCTTrackDeltaXVertexPtEta"},
    {kFCTTrackDeltaYVertexPtEta, "FCTTrackDeltaYVertexPtEta"},
    {kFCTTrackPtResolutionPtEta, "FCTTrackPtResolutionPtEta"},
    {kFCTTrackInvQPtPullPtEta, "FCTTrackInvQPtPullPtEta"},
    {kFCTTrackInvQPtResolutionPtEta, "FCTTrackInvQPtResolutionPtEta"},
    {kFCTTrackInvPtResolutionPtEta, "FCTTrackInvPtResolutionPtEta"}};

  std::map<int, std::array<double, 9>> TH3Binning{
    {kFCTTrackDeltaYVertexPtEta, {105, 0, 21, 62, 1.7, 4.5, 2e4, -1e4, 1e4}},
    {kFCTTrackDeltaXVertexPtEta, {105, 0, 21, 62, 1.7, 4.5, 2e4, -1e4, 14}},
    {kFCTTrackPtResolutionPtEta, {105, 0, 21, 62, 1.7, 4.5, 1000, -2, 50}},
    {kFCTTrackInvQPtPullPtEta, {105, 0, 21, 62, 1.7, 4.5, 200, -5, 5}},
    {kFCTTrackInvQPtResolutionPtEta, {105, 0, 21, 62, 1.7, 4.5, 2000, -2, 2}},
    {kFCTTrackInvPtResolutionPtEta, {105, 0, 21, 62, 1.7, 4.5, 2500, -5, 150}}};

  std::map<int, const char*> TH3XaxisTitles{
    {kFCTTrackDeltaXVertexPtEta, "p_t"},
    {kFCTTrackDeltaYVertexPtEta, "p_t"},
    {kFCTTrackPtResolutionPtEta, "p_t"},
    {kFCTTrackInvQPtPullPtEta, "p_t"},
    {kFCTTrackInvQPtResolutionPtEta, "p_t"},
    {kFCTTrackInvPtResolutionPtEta, "p_t"}};

  //
  std::map<int, const char*> TH3YaxisTitles{
    {kFCTTrackDeltaXVertexPtEta, "\\eta"},
    {kFCTTrackDeltaYVertexPtEta, "\\eta"},
    {kFCTTrackPtResolutionPtEta, "\\eta"},
    {kFCTTrackInvQPtPullPtEta, "\\eta"},
    {kFCTTrackInvQPtResolutionPtEta, "\\eta"},
    {kFCTTrackInvPtResolutionPtEta, "\\eta"}};

  std::map<int, const char*> TH3ZaxisTitles{
    {kFCTTrackDeltaXVertexPtEta, "X residual at vertex (um)"},
    {kFCTTrackDeltaYVertexPtEta, "Y residual at vertex (um)"},
    {kFCTTrackPtResolutionPtEta, "(p_t residual)/pt"},
    {kFCTTrackInvQPtPullPtEta, "(\\Delta q/p_t)/\\sigma_{q/pt}"},
    {kFCTTrackInvQPtResolutionPtEta, "(q/p_t residual)/(q/pt)"},
    {kFCTTrackInvPtResolutionPtEta, "(1/p_t residual)/(1/pt)"}};

  enum TH2HistosCodes {
    kFCTTrackDeltaXYVertex,
    kFCTTrackDeltaXYVertex0_1,
    kFCTTrackDeltaXYVertex1_4,
    kFCTTrackDeltaXYVertex4plus,
    kFCTTrackQPRec_MC,
    kFCTTrackPtResolution,
    kFCTTrackPtResolutionInner,
    kFCTTrackPtResolutionOuter,
    kFCTTrackInvPtResolution,
    kFCTTrackInvPtResolutionInner,
    kFCTTrackInvPtResolutionOuter,
    kMCTracksEtaZ
  };

  std::map<int, const char*> TH2Names{
    {kFCTTrackDeltaXYVertex, "FCT Tracks at vertex"},
    {kFCTTrackDeltaXYVertex0_1, "FCT Tracks Vertex at Z = 0 Pt0_1"},
    {kFCTTrackDeltaXYVertex1_4, "FCT Tracks Vertex at Z = 0 Pt1_4"},
    {kFCTTrackDeltaXYVertex4plus, "FCT Tracks Vertex at Z = 0 Pt4plus"},
    {kFCTTrackQPRec_MC, "FCT Track QP FITxMC"},
    {kFCTTrackPtResolution, "FCT Track Pt Resolution"},
    {kFCTTrackPtResolutionInner, "FCT Track Pt Resolution Inner"},
    {kFCTTrackPtResolutionOuter, "FCT Track Pt Resolution Outer"},
    {kFCTTrackInvPtResolution, "FCT Track InvPt Resolution"},
    {kFCTTrackInvPtResolutionInner, "FCT Track InvPt ResolutionInner"},
    {kFCTTrackInvPtResolutionOuter, "FCT Track InvPt ResolutionOuter"},
    {kMCTracksEtaZ, "MCTracks_eta_z"}};

  std::map<int, const char*> TH2Titles{
    {kFCTTrackDeltaXYVertex, "FCT Tracks at Z_vertex"},
    {kFCTTrackDeltaXYVertex0_1, "FCT Tracks at Z_vertex (pt < 1)"},
    {kFCTTrackDeltaXYVertex1_4, "FCT Tracks at Z_vertex (1 < pt < 4)"},
    {kFCTTrackDeltaXYVertex4plus, "FCT Tracks at Z_vertex (pt > 4)"},
    {kFCTTrackQPRec_MC, "Charged Momentum: Reconstructed vs MC"},
    {kFCTTrackPtResolution, "Pt Resolution"},
    {kFCTTrackPtResolutionInner, "\\text{Pt Resolution }(3.5 < \\eta < 3.6)"},
    {kFCTTrackPtResolutionOuter, "\\text{Pt Resolution }(2.8 < \\eta < 2.9 )"},
    {kFCTTrackInvPtResolution, "InvPt Resolution"},
    {kFCTTrackInvPtResolutionInner, "\\text{InvPt Resolution }(3.5 < \\eta < 3.6)"},
    {kFCTTrackInvPtResolutionOuter, "\\text{InvPt Resolution }(2.8 < \\eta < 2.9 )"},
    {kMCTracksEtaZ, "MC Tracks: Pseudorapidity vs zVertex"}};

  std::map<int, std::array<double, 6>> TH2Binning{
    {kFCTTrackDeltaXYVertex, {100, -.5, .5, 100, -.5, .5}},
    {kFCTTrackDeltaXYVertex0_1, {100, -.5, .5, 100, -.5, .5}},
    {kFCTTrackDeltaXYVertex1_4, {100, -.5, .5, 100, -.5, .5}},
    {kFCTTrackDeltaXYVertex4plus, {100, -.5, .5, 100, -.5, .5}},
    {kFCTTrackQPRec_MC, {100, -10, 10, 100, -10, 10}},
    {kFCTTrackPtResolution, {10, 0, 10, 1000, -2, 50}},
    {kFCTTrackPtResolutionInner, {10, 0, 10, 1000, -2, 50}},
    {kFCTTrackPtResolutionOuter, {10, 0, 10, 1000, -2, 50}},
    {kFCTTrackInvPtResolution, {10, 0, 10, 2500, -5, 150}},
    {kFCTTrackInvPtResolutionInner, {10, 0, 10, 2500, -5, 150}},
    {kFCTTrackInvPtResolutionOuter, {10, 0, 10, 2500, -5, 150}},
    {kMCTracksEtaZ, {31, -15, 16, 25, etaMin, etaMax}}};

  std::map<int, const char*> TH2XaxisTitles{
    {kFCTTrackDeltaXYVertex, "\\Delta x \\text{ [mm]}"},
    {kFCTTrackDeltaXYVertex0_1, "\\Delta x \\text{ [mm]}"},
    {kFCTTrackDeltaXYVertex1_4, "\\Delta x \\text{ [mm]}"},
    {kFCTTrackDeltaXYVertex4plus, "\\Delta x \\text{ [mm]}"},
    {kFCTTrackQPRec_MC, "(q.p)_{MC} [GeV]"},
    {kFCTTrackPtResolution, "pt_{MC} [GeV]"},
    {kFCTTrackPtResolutionInner, "pt_{MC} [GeV]"},
    {kFCTTrackPtResolutionOuter, "pt_{MC} [GeV]"},
    {kFCTTrackInvPtResolution, "pt_{MC} [GeV]"},
    {kFCTTrackInvPtResolutionInner, "pt_{MC} [GeV]"},
    {kFCTTrackInvPtResolutionOuter, "pt_{MC} [GeV]"},
    {kMCTracksEtaZ, "Vertex PosZ [cm]"}};

  std::map<int, const char*> TH2YaxisTitles{
    {kFCTTrackDeltaXYVertex, "\\Delta y \\text{ [mm]}"},
    {kFCTTrackDeltaXYVertex0_1, "\\Delta y \\text{ [mm]}"},
    {kFCTTrackDeltaXYVertex1_4, "\\Delta y \\text{ [mm]}"},
    {kFCTTrackDeltaXYVertex4plus, "\\Delta y \\text{ [mm]}"},
    {kFCTTrackQPRec_MC, "(q.p)_{fit} [GeV]"},
    {kFCTTrackPtResolution, "(pt_{fit} - pt_{MC}) / pt_{MC}"},
    {kFCTTrackPtResolutionInner, "(pt_{fit} - pt_{MC}) / pt_{MC}"},
    {kFCTTrackPtResolutionOuter, "(pt_{fit} - pt_{MC}) / pt_{MC}"},
    {kFCTTrackInvPtResolution, "(1/(p_t)_{fit} - 1/(p_t)_{MC})/(1/(p_t)_{MC})"},
    {kFCTTrackInvPtResolutionInner, "(1/(p_t)_{fit} - 1/(p_t)_{MC})/(1/(p_t)_{MC})"},
    {kFCTTrackInvPtResolutionOuter, "(1/(p_t)_{fit} - 1/(p_t)_{MC})/(1/(p_t)_{MC})"},
    {kMCTracksEtaZ, "\\eta"}};

  enum TH1HistosCodes {
    kFCTTrackDeltaXErr,
    kFCTTrackDeltaYErr,
    kFCTTrackDeltaPhiErr,
    kFCTTrackDeltaTanLErr,
    kFCTTrackDeltainvQPtErr,
    kFCTTracksP,
    kFCTTrackDeltaTanl,
    kFCTTrackDeltaTanl0_1,
    kFCTTrackDeltaTanl1_4,
    kFCTTrackDeltaTanl4plus,
    kFCTTrackDeltaPhi,
    kFCTTrackDeltaPhi0_1,
    kFCTTrackDeltaPhi1_4,
    kFCTTrackDeltaPhi4plus,
    kFCTTrackDeltaPhiDeg,
    kFCTTrackDeltaPhiDeg0_1,
    kFCTTrackDeltaPhiDeg1_4,
    kFCTTrackDeltaPhiDeg4plus,
    kFCTTrackDeltaInvQPt,
    kFCTTrackDeltaInvQPtSeed,
    kFCTTrackDeltaX,
    kFCTTrackDeltaX0_1,
    kFCTTrackDeltaX1_4,
    kFCTTrackDeltaX4plus,
    kFCTTrackDeltaY,
    kFCTTrackR,
    kFCTTrackQ,
    kFCTTrackQ0_1,
    kFCTTrackQ1_4,
    kFCTTrackQ4plus,
    kFCTTrackXPull1_innerBorder,
    kFCTTrackXPull1_OuterBorder,
    kFCTTrackXPull4_innerBorder,
    kFCTTrackXPull4_OuterBorder,
    kFCTTrackYPull1_innerBorder,
    kFCTTrackYPull1_OuterBorder,
    kFCTTrackYPull4_innerBorder,
    kFCTTrackYPull4_OuterBorder,
    kFCTTrackPhiPull1_innerBorder,
    kFCTTrackPhiPull1_OuterBorder,
    kFCTTrackPhiPull4_innerBorder,
    kFCTTrackPhiPull4_OuterBorder,
    kFCTTrackTanlPull1_innerBorder,
    kFCTTrackTanlPull1_OuterBorder,
    kFCTTrackTanlPull4_innerBorder,
    kFCTTrackTanlPull4_OuterBorder,
    kFCTTrackInvQPtPull1_innerBorder,
    kFCTTrackInvQPtPull1_OuterBorder,
    kFCTTrackInvQPtPull4_innerBorder,
    kFCTTrackInvQPtPull4_OuterBorder,
    kFCTTrackChi2,
    kMCTrackspT,
    kMCTracksp,
    kMCTrackEta
  };

  std::map<int, const char*> TH1Names{
    {kFCTTracksP, "FCT Tracks Fitted p"},
    {kFCTTrackDeltaXErr, "Delta X / SigmaX"},
    {kFCTTrackDeltaYErr, "Delta Y / SigmaY"},
    {kFCTTrackDeltaPhiErr, "Delta Phi at Vertex / SigmaPhi"},
    {kFCTTrackDeltaTanLErr, "Delta_Tanl / SigmaTanl"},
    {kFCTTrackDeltainvQPtErr, "Delta_InvQPt / Sigma_{q/pt}"},
    {kFCTTrackDeltaTanl, "FCT Tracks Fitted Delta_tanl"},
    {kFCTTrackDeltaTanl0_1, "FCT Tracks tanl (pt < 1)"},
    {kFCTTrackDeltaTanl1_4, "FCT Tracks tanl (1 < pt < 4)"},
    {kFCTTrackDeltaTanl4plus, "FCT Tracks tanl (pt > 4)"},
    {kFCTTrackDeltaPhi, "FCT Tracks Fitted Phi at Vertex"},
    {kFCTTrackDeltaPhi0_1, "FCT Tracks Fitted Phi at Vertex [rad] (pt < 1)"},
    {kFCTTrackDeltaPhi1_4,
     "FCT Tracks Fitted Phi at Vertex [rad] (1 < pt < 4)"},
    {kFCTTrackDeltaPhi4plus,
     "FCT Tracks Fitted Phi at Vertex [rad] (pt > 4)"},
    {kFCTTrackDeltaPhiDeg, "FCT Tracks Fitted Phi at Vertex [deg]"},
    {kFCTTrackDeltaPhiDeg0_1,
     "FCT Tracks Fitted Phi at Vertex [deg] (pt < 1)"},
    {kFCTTrackDeltaPhiDeg1_4,
     "FCT Tracks Fitted Phi at Vertex [deg] (1 < pt < 4)"},
    {kFCTTrackDeltaPhiDeg4plus,
     "FCT Tracks Fitted Phi at Vertex [deg] (pt > 4)"},
    {kFCTTrackDeltaInvQPt, "FCT Tracks invQPt"},
    {kFCTTrackDeltaInvQPtSeed, "FCT Tracks invQPt Seed"},
    {kFCTTrackDeltaX, "FCT Tracks Delta X"},
    {kFCTTrackDeltaX0_1, "FCT Tracks Delta X (pt < 1)"},
    {kFCTTrackDeltaX1_4, "FCT Tracks Delta X (1 < pt < 4)"},
    {kFCTTrackDeltaX4plus, "FCT Tracks Delta X (pt > 4)"},
    {kFCTTrackDeltaY, "FCT Tracks Delta Y"},
    {kFCTTrackR, "FCT Tracks Delta R"},
    {kFCTTrackQ, "Charge Match"},
    {kFCTTrackQ0_1, "Charge Match (pt < 1)"},
    {kFCTTrackQ1_4, "Charge Match (1 < pt < 4)"},
    {kFCTTrackQ4plus, "Charge Match (pt > 4)"},
    {kFCTTrackXPull1_innerBorder, "TrackXPull1_innerBorder"},
    {kFCTTrackXPull1_OuterBorder, "TrackXPull1_OuterBorder"},
    {kFCTTrackXPull4_innerBorder, "TrackXPull4_innerBorder"},
    {kFCTTrackXPull4_OuterBorder, "TrackXPull4_OuterBorder"},
    {kFCTTrackYPull1_innerBorder, "TrackYPull1_innerBorder"},
    {kFCTTrackYPull1_OuterBorder, "TrackYPull1_OuterBorder"},
    {kFCTTrackYPull4_innerBorder, "TrackYPull4_innerBorder"},
    {kFCTTrackYPull4_OuterBorder, "TrackYPull4_OuterBorder"},
    {kFCTTrackPhiPull1_innerBorder, "TrackPhiPull1_innerBorder"},
    {kFCTTrackPhiPull1_OuterBorder, "TrackPhiPull1_OuterBorder"},
    {kFCTTrackPhiPull4_innerBorder, "TrackPhiPull4_innerBorder"},
    {kFCTTrackPhiPull4_OuterBorder, "TrackPhiPull4_OuterBorder"},
    {kFCTTrackTanlPull1_innerBorder, "TrackTanlPull1_innerBorder"},
    {kFCTTrackTanlPull1_OuterBorder, "TrackTanlPull1_OuterBorder"},
    {kFCTTrackTanlPull4_innerBorder, "TrackTanlPull4_innerBorder"},
    {kFCTTrackTanlPull4_OuterBorder, "TrackTanlPull4_OuterBorder"},
    {kFCTTrackInvQPtPull1_innerBorder, "TrackInvQPtPull1_innerBorder"},
    {kFCTTrackInvQPtPull1_OuterBorder, "TrackInvQPtPull1_OuterBorder"},
    {kFCTTrackInvQPtPull4_innerBorder, "TrackInvQPtPull4_innerBorder"},
    {kFCTTrackInvQPtPull4_OuterBorder, "TrackInvQPtPull4_OuterBorder"},
    {kFCTTrackChi2, "Tracks Chi2"},
    {kMCTrackspT, "MC Tracks p_T"},
    {kMCTracksp, "MC Tracks p"},
    {kMCTrackEta, "MC Tracks eta"}};

  std::map<int, const char*> TH1Titles{
    {kFCTTracksP, "Standalone FCT Tracks P"},
    {kFCTTrackDeltaXErr, "\\Delta X / \\sigma_X"},
    {kFCTTrackDeltaYErr, "\\Delta Y / \\sigma_Y"},
    {kFCTTrackDeltaPhiErr, "\\Delta \\phi / \\sigma_\\phi"},
    {kFCTTrackDeltaTanLErr, "\\Delta TanL / \\sigma_{TanL} "},
    {kFCTTrackDeltainvQPtErr, "\\Delta(q/Pt) / \\sigma_{q/pt}"},
    {kFCTTrackDeltaTanl, "tanl_{Fit} - tanl_{MC} "},
    {kFCTTrackDeltaTanl0_1, "tanl_{Fit} - tanl_{MC} (pt < 1)"},
    {kFCTTrackDeltaTanl1_4, "tanl_{Fit} - tanl_{MC} (1 < p_t < 4)"},
    {kFCTTrackDeltaTanl4plus, "tanl_{Fit} - tanl_{MC} (p_t > 4)"},
    {kFCTTrackDeltaPhi, "\\phi _{Fit} - \\phi_{MC}"},
    {kFCTTrackDeltaPhi0_1, "\\phi _{Fit} - \\phi_{MC}"},
    {kFCTTrackDeltaPhi1_4, "\\phi _{Fit} - \\phi_{MC}"},
    {kFCTTrackDeltaPhi4plus, "\\phi _{Fit} - \\phi_{MC}"},
    {kFCTTrackDeltaPhiDeg, "\\phi _{Fit} - \\phi_{MC}"},
    {kFCTTrackDeltaPhiDeg0_1, "\\phi _{Fit} - \\phi_{MC}"},
    {kFCTTrackDeltaPhiDeg1_4, "\\phi _{Fit} - \\phi_{MC}"},
    {kFCTTrackDeltaPhiDeg4plus, "\\phi _{Fit} - \\phi_{MC}"},
    {kFCTTrackDeltaInvQPt, "FCT Tracks \\Delta invQPt"},
    {kFCTTrackDeltaInvQPtSeed, "FCT Tracks \\Delta invQPt Seed"},
    {kFCTTrackDeltaX, "FCT Tracks Delta X at Z_vertex"},
    {kFCTTrackDeltaX0_1, "FCT Tracks Delta X at Z_vertex"},
    {kFCTTrackDeltaX1_4, "FCT Tracks Delta X at Z_vertex"},
    {kFCTTrackDeltaX4plus, "FCT Tracks Delta X at Z_vertex"},
    {kFCTTrackDeltaY, "FCT Tracks Delta Y at Z_vertex"},
    {kFCTTrackR, "FCT Tracks Delta R at Z_vertex"},
    {kFCTTrackQ, "FCT Tracks Charge Match"},
    {kFCTTrackQ0_1, "FCT Tracks Charge Match (pt < 1)"},
    {kFCTTrackQ1_4, "FCT Tracks Charge Match (1 < pt < 4)"},
    {kFCTTrackQ4plus, "FCT Tracks Charge Match (pt > 4)"},
    {kFCTTrackChi2, "FCT Tracks ~ \\chi^2"},
    {kFCTTrackXPull1_innerBorder, "\\text{TrackXPull1GeV } \\eta = -3.6"},
    {kFCTTrackXPull1_OuterBorder, "\\text{TrackXPull1GeV } \\eta = -2.8"},
    {kFCTTrackXPull4_innerBorder, "\\text{TrackXPull4GeV } \\eta = -3.6"},
    {kFCTTrackXPull4_OuterBorder, "\\text{TrackXPull4GeV } \\eta = -2.8"},
    {kFCTTrackYPull1_innerBorder, "\\text{TrackYPull1GeV } \\eta = -3.6"},
    {kFCTTrackYPull1_OuterBorder, "\\text{TrackYPull1GeV } \\eta = -2.8"},
    {kFCTTrackYPull4_innerBorder, "\\text{TrackYPull4GeV } \\eta = -3.6"},
    {kFCTTrackYPull4_OuterBorder, "\\text{TrackYPull4GeV } \\eta = -2.8"},
    {kFCTTrackPhiPull1_innerBorder, "\\text{TrackPhiPull1GeV } \\eta = -3.6"},
    {kFCTTrackPhiPull1_OuterBorder, "\\text{TrackPhiPull1GeV } \\eta = -2.8"},
    {kFCTTrackPhiPull4_innerBorder, "\\text{TrackPhiPull4GeV } \\eta = -3.6"},
    {kFCTTrackPhiPull4_OuterBorder, "\\text{TrackPhiPull4GeV } \\eta = -2.8"},
    {kFCTTrackTanlPull1_innerBorder, "\\text{TrackTanlPull1GeV } \\eta = -3.6"},
    {kFCTTrackTanlPull1_OuterBorder, "\\text{TrackTanlPull1GeV } \\eta = -2.8"},
    {kFCTTrackTanlPull4_innerBorder, "\\text{TrackTanlPull4GeV } \\eta = -3.6"},
    {kFCTTrackTanlPull4_OuterBorder, "\\text{TrackTanlPull4GeV } \\eta = -2.8"},
    {kFCTTrackInvQPtPull1_innerBorder, "\\text{TrackInvQPtPull1GeV } \\eta = -3.6"},
    {kFCTTrackInvQPtPull1_OuterBorder, "\\text{TrackInvQPtPull1GeV } \\eta = -2.8"},
    {kFCTTrackInvQPtPull4_innerBorder, "\\text{TrackInvQPtPull4GeV } \\eta = -3.6"},
    {kFCTTrackInvQPtPull4_OuterBorder, "\\text{TrackInvQPtPull4GeV } \\eta = -2.8"},
    {kMCTrackspT, "MC Tracks p_T"},
    {kMCTracksp, "MC Tracks p"},
    {kMCTrackEta, "MC Tracks Pseudorapidity"}};

  std::map<int, std::array<double, 3>> TH1Binning{
    {kFCTTracksP, {500, pMin, pMax}},
    {kFCTTrackDeltaXErr, {500, -5, 5}},
    {kFCTTrackDeltaYErr, {500, -5, 5}},
    {kFCTTrackDeltaPhiErr, {200, -5, 5}},
    {kFCTTrackDeltaTanLErr, {200, -5, 5}},
    {kFCTTrackDeltainvQPtErr, {200, -5, 5}},
    {kFCTTrackDeltaTanl, {200, deltatanlMin, deltatanlMax}},
    {kFCTTrackDeltaTanl0_1, {200, deltatanlMin, deltatanlMax}},
    {kFCTTrackDeltaTanl1_4, {200, deltatanlMin, deltatanlMax}},
    {kFCTTrackDeltaTanl4plus, {200, deltatanlMin, deltatanlMax}},
    {kFCTTrackDeltaPhi, {200, deltaphiMin, deltaphiMax}},
    {kFCTTrackDeltaPhi0_1, {200, deltaphiMin, deltaphiMax}},
    {kFCTTrackDeltaPhi1_4, {200, deltaphiMin, deltaphiMax}},
    {kFCTTrackDeltaPhi4plus, {200, deltaphiMin, deltaphiMax}},
    {kFCTTrackDeltaPhiDeg,
     {1000, TMath::RadToDeg() * deltaphiMin,
      TMath::RadToDeg() * deltaphiMax}},
    {kFCTTrackDeltaPhiDeg0_1,
     {1000, TMath::RadToDeg() * deltaphiMin,
      TMath::RadToDeg() * deltaphiMax}},
    {kFCTTrackDeltaPhiDeg1_4,
     {1000, TMath::RadToDeg() * deltaphiMin,
      TMath::RadToDeg() * deltaphiMax}},
    {kFCTTrackDeltaPhiDeg4plus,
     {1000, TMath::RadToDeg() * deltaphiMin,
      TMath::RadToDeg() * deltaphiMax}},
    {kFCTTrackDeltaInvQPt, {1000, -10., 10.}},
    {kFCTTrackDeltaInvQPtSeed, {1000, -0.5, 0.5}},
    {kFCTTrackDeltaX, {1000, -5, 5}},
    {kFCTTrackDeltaX0_1, {1000, -5, 5}},
    {kFCTTrackDeltaX1_4, {1000, -5, 5}},
    {kFCTTrackDeltaX4plus, {1000, -5, 5}},
    {kFCTTrackDeltaY, {1000, -5, 5}},
    {kFCTTrackR, {250, 0, 5}},
    {kFCTTrackQ, {5, -2.1, 2.1}},
    {kFCTTrackQ0_1, {5, -2.1, 2.1}},
    {kFCTTrackQ1_4, {5, -2.1, 2.1}},
    {kFCTTrackQ4plus, {5, -2.1, 2.1}},
    {kFCTTrackChi2, {10000, 0, 1000}},
    {kFCTTrackXPull1_innerBorder, {200, -5, 5}},
    {kFCTTrackXPull1_OuterBorder, {200, -5, 5}},
    {kFCTTrackXPull4_innerBorder, {200, -5, 5}},
    {kFCTTrackXPull4_OuterBorder, {200, -5, 5}},
    {kFCTTrackYPull1_innerBorder, {200, -5, 5}},
    {kFCTTrackYPull1_OuterBorder, {200, -5, 5}},
    {kFCTTrackYPull4_innerBorder, {200, -5, 5}},
    {kFCTTrackYPull4_OuterBorder, {200, -5, 5}},
    {kFCTTrackPhiPull1_innerBorder, {200, -5, 5}},
    {kFCTTrackPhiPull1_OuterBorder, {200, -5, 5}},
    {kFCTTrackPhiPull4_innerBorder, {200, -5, 5}},
    {kFCTTrackPhiPull4_OuterBorder, {200, -5, 5}},
    {kFCTTrackTanlPull1_innerBorder, {200, -5, 5}},
    {kFCTTrackTanlPull1_OuterBorder, {200, -5, 5}},
    {kFCTTrackTanlPull4_innerBorder, {200, -5, 5}},
    {kFCTTrackTanlPull4_OuterBorder, {200, -5, 5}},
    {kFCTTrackInvQPtPull1_innerBorder, {200, -5, 5}},
    {kFCTTrackInvQPtPull1_OuterBorder, {200, -5, 5}},
    {kFCTTrackInvQPtPull4_innerBorder, {200, -5, 5}},
    {kFCTTrackInvQPtPull4_OuterBorder, {200, -5, 5}},
    {kMCTrackspT, {5000, 0, 50}},
    {kMCTracksp, {1000, pMin, pMax}},
    {kMCTrackEta, {1000, etaMin, etaMax}}};

  std::map<int, const char*> TH1XaxisTitles{
    {kFCTTracksP, "p [GeV]"},
    {kFCTTrackDeltaXErr, "\\Delta x  /\\sigma_{x}"},
    {kFCTTrackDeltaYErr, "\\Delta y  /\\sigma_{y}"},
    {kFCTTrackDeltaPhiErr, "\\Delta \\phi  /\\sigma_{\\phi}"},
    {kFCTTrackDeltaTanLErr, "\\Delta tanl /\\sigma_{tanl}"},
    {kFCTTrackDeltainvQPtErr, "\\Delta (q/p_t)/\\sigma_{q/Pt}"},
    {kFCTTrackDeltaTanl, "\\Delta tanl"},
    {kFCTTrackDeltaTanl0_1, "\\Delta tanl"},
    {kFCTTrackDeltaTanl1_4, "\\Delta tanl"},
    {kFCTTrackDeltaTanl4plus, "\\Delta tanl"},
    {kFCTTrackDeltaPhi, "\\Delta \\phi ~[rad]"},
    {kFCTTrackDeltaPhi0_1, "\\Delta \\phi ~[rad]"},
    {kFCTTrackDeltaPhi1_4, "\\Delta \\phi ~[rad]"},
    {kFCTTrackDeltaPhi4plus, "\\Delta \\phi ~[rad]"},
    {kFCTTrackDeltaPhiDeg, "\\Delta \\phi ~[deg]"},
    {kFCTTrackDeltaPhiDeg0_1, "\\Delta \\phi ~[deg]"},
    {kFCTTrackDeltaPhiDeg1_4, "\\Delta \\phi ~[deg]"},
    {kFCTTrackDeltaPhiDeg4plus, "\\Delta \\phi ~[deg]"},
    {kFCTTrackDeltaInvQPt, "\\Delta invQPt"},
    {kFCTTrackDeltaInvQPtSeed, "\\Delta invQPt Seed"},
    {kFCTTrackDeltaX, "\\Delta x \\text{ [mm]}"},
    {kFCTTrackDeltaX0_1, "\\Delta x \\text{ [mm]}"},
    {kFCTTrackDeltaX1_4, "\\Delta x \\text{ [mm]}"},
    {kFCTTrackDeltaX4plus, "\\Delta x \\text{ [mm]}"},
    {kFCTTrackDeltaY, "\\Delta y \\text{ [mm]}"},
    {kFCTTrackR, "\\Delta r \\text{ [mm]}"},
    {kFCTTrackQ, "q_{fit}-q_{MC}"},
    {kFCTTrackQ0_1, "q_{fit}-q_{MC}"},
    {kFCTTrackQ1_4, "q_{fit}-q_{MC}"},
    {kFCTTrackQ4plus, "q_{fit}-q_{MC}"},
    {kFCTTrackChi2, "\\chi^2"},
    {kFCTTrackXPull1_innerBorder, "x_pull"},
    {kFCTTrackXPull1_OuterBorder, "x_pull"},
    {kFCTTrackXPull4_innerBorder, "x_pull"},
    {kFCTTrackXPull4_OuterBorder, "x_pull"},
    {kFCTTrackYPull1_innerBorder, "y_pull"},
    {kFCTTrackYPull1_OuterBorder, "y_pull"},
    {kFCTTrackYPull4_innerBorder, "y_pull"},
    {kFCTTrackYPull4_OuterBorder, "y_pull"},
    {kFCTTrackPhiPull1_innerBorder, "phi_pull"},
    {kFCTTrackPhiPull1_OuterBorder, "phi_pull"},
    {kFCTTrackPhiPull4_innerBorder, "phi_pull"},
    {kFCTTrackPhiPull4_OuterBorder, "phi_pull"},
    {kFCTTrackTanlPull1_innerBorder, "tanl_pull"},
    {kFCTTrackTanlPull1_OuterBorder, "tanl_pull"},
    {kFCTTrackTanlPull4_innerBorder, "tanl_pull"},
    {kFCTTrackTanlPull4_OuterBorder, "tanl_pull"},
    {kFCTTrackInvQPtPull1_innerBorder, "q/pt_pull"},
    {kFCTTrackInvQPtPull1_OuterBorder, "q/pt_pull"},
    {kFCTTrackInvQPtPull4_innerBorder, "q/pt_pull"},
    {kFCTTrackInvQPtPull4_OuterBorder, "q/pt_pull"},
    {kMCTrackspT, "p_t [GeV]"},
    {kMCTracksp, "p [GeV]"},
    {kMCTrackEta, " \\eta"}};

  // Create histograms
  const int nTH1Histos = TH1Names.size();
  std::vector<std::unique_ptr<TH1F>> TH1Histos(nTH1Histos);
  auto nHisto = 0;
  for (auto& h : TH1Histos) {
    h = std::make_unique<TH1F>(TH1Names[nHisto], TH1Titles[nHisto],
                               (int)TH1Binning[nHisto][0],
                               TH1Binning[nHisto][1], TH1Binning[nHisto][2]);
    h->GetXaxis()->SetTitle(TH1XaxisTitles[nHisto]);
    ++nHisto;
  }

  const int nTH2Histos = TH2Names.size();
  std::vector<std::unique_ptr<TH2F>> TH2Histos(nTH2Histos);
  auto n2Histo = 0;
  for (auto& h : TH2Histos) {
    h = std::make_unique<TH2F>(TH2Names[n2Histo], TH2Titles[n2Histo],
                               (int)TH2Binning[n2Histo][0],
                               TH2Binning[n2Histo][1], TH2Binning[n2Histo][2],
                               (int)TH2Binning[n2Histo][3],
                               TH2Binning[n2Histo][4], TH2Binning[n2Histo][5]);
    h->GetXaxis()->SetTitle(TH2XaxisTitles[n2Histo]);
    h->GetYaxis()->SetTitle(TH2YaxisTitles[n2Histo]);

    h->SetOption("COLZ");
    ++n2Histo;
  }

  const int nTH3Histos = TH3Names.size();
  std::vector<std::unique_ptr<TH3F>> TH3Histos(nTH3Histos);
  auto n3Histo = 0;
  for (auto& h : TH3Histos) {
    h = std::make_unique<TH3F>(TH3Names[n3Histo], TH3Titles[n3Histo],
                               (int)TH3Binning[n3Histo][0],
                               TH3Binning[n3Histo][1],
                               TH3Binning[n3Histo][2],
                               (int)TH3Binning[n3Histo][3],
                               TH3Binning[n3Histo][4],
                               TH3Binning[n3Histo][5],
                               (int)TH3Binning[n3Histo][6],
                               TH3Binning[n3Histo][7],
                               TH3Binning[n3Histo][8]);
    h->GetXaxis()->SetTitle(TH3XaxisTitles[n3Histo]);
    h->GetYaxis()->SetTitle(TH3YaxisTitles[n3Histo]);
    h->GetZaxis()->SetTitle(TH3ZaxisTitles[n3Histo]);

    //h->SetOption("COLZ");
    ++n3Histo;
  }

  // Profiles histograms
  auto PtRes_Profile = new TProfile("Pt_res_prof", "Profile of pt{fit}/pt{MC}",
                                    10, 0, 10, 0, 20, "s");
  PtRes_Profile->GetXaxis()->SetTitle("pt_{MC}");
  PtRes_Profile->GetYaxis()->SetTitle("mean(Pt_{Fit}/Pt_{MC})");

  auto DeltaX_Profile = new TProfile("DeltaX_prof", "Vertexing resolution", 14,
                                     0, 7, -10000., 10000., "s");
  DeltaX_Profile->GetXaxis()->SetTitle("pt_{MC} [GeV]");
  DeltaX_Profile->GetYaxis()->SetTitle("\\sigma_x ~[\\mu m]");

  // TEfficiency histogram
  TEfficiency* qMatchEff = new TEfficiency(
    "QMatchEff", "Charge Match;p_t [GeV];#epsilon", 10, 0, 10);

  // Counters
  Int_t nChargeMatch = 0;
  Int_t nChargeMiss = 0;
  Int_t nChargeMatch0_1 = 0;
  Int_t nChargeMiss0_1 = 0;
  Int_t nChargeMatch1_4 = 0;
  Int_t nChargeMiss1_4 = 0;
  Int_t nChargeMatch4plus = 0;
  Int_t nChargeMiss4plus = 0;

  // Files & Trees
  // MC
  TFile* o2sim_KineFileIn = new TFile(o2sim_KineFile);
  TTree* o2SimKineTree = (TTree*)o2sim_KineFileIn->Get("o2sim");

  vector<MCTrackT<float>>* mcTr = nullptr;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  // FCT Hits
  TFile* HitsFCTFileIn = new TFile(HitsFCTFile);
  TTree* o2FCTHitsTree = (TTree*)HitsFCTFileIn->Get("o2sim");
  vector<Hit>* fcthit = nullptr;
  o2FCTHitsTree->SetBranchAddress("FCTHit", &fcthit);

  Int_t numberOfEvents = o2SimKineTree->GetEntries();
  Int_t numberOfFCTEvents = o2FCTHitsTree->GetEntries();
  if (numberOfEvents == numberOfFCTEvents)
    std::cout << "numberOfEvents = " << numberOfEvents << std::endl;
  else {
    std::cout << "ERROR: Inconsistent number of entries on " << o2sim_KineFile
              << " and " << HitsFCTFile << std::endl;
    return -1;
  }

  // FCT Tracks
  TFile* trkFileIn = new TFile(trkFile);
  TTree* fctTrackTree = (TTree*)trkFileIn->Get("o2sim");
  std::vector<o2::fct::FCTTrackExt> trackFCTVec, *trackFCTVecP = &trackFCTVec;
  fctTrackTree->SetBranchAddress("FCTTrack", &trackFCTVecP);

  vector<Int_t>* recoTrackIDs = nullptr;
  fctTrackTree->SetBranchAddress("FCTTrackID", &recoTrackIDs);

  fctTrackTree->GetEntry(0);
  o2SimKineTree->GetEntry(0);

  auto field_z = getZField(0, 0, 0); // Get field at Center of ALICE

  std::string outfilename = "Fittercheck_" + std::string(trkFile);

  TFile outFile(outfilename.c_str(), "RECREATE");

  // Reconstructed FCT Tracks
  std::cout << "Loop over events and reconstructed FCT Tracks!" << std::endl;
  // TracksFCT - Identify reconstructed tracks
  auto totalTracks = 0;
  for (int iEvent = 0; iEvent < numberOfEvents; iEvent++) {
    fctTrackTree->GetEntry(iEvent);
    o2SimKineTree->GetEntry(iEvent);
    auto iTrack = 0;
    //if (DEBUG_VERBOSE)
    std::cout << "Processing Event # " << iEvent << std::endl;
    o2SimKineTree->GetEntry(iEvent);
    for (auto& trackFCT : trackFCTVec) {
      auto trackID = recoTrackIDs->at(iTrack);
      if (trackFCT.getNumberOfPoints() < minHitsPerTrack) {
        iTrack++;
        continue;
      }
      if (1) {
        if (DEBUG_VERBOSE) {

          std::cout << "  Track #" << iTrack << ": TrackID = " << trackID
                    << std::endl;
        }

        MCTrackT<float>* thisTrack = &(*mcTr).at(trackID);
        auto vx_MC = thisTrack->GetStartVertexCoordinatesX();
        auto vy_MC = thisTrack->GetStartVertexCoordinatesY();
        auto vz_MC = thisTrack->GetStartVertexCoordinatesZ();
        auto Pt_MC = thisTrack->GetPt();
        auto P_MC = thisTrack->GetP();
        auto phi_MC = TMath::ATan2(thisTrack->Py(), thisTrack->Px());
        auto eta_MC = atanh(thisTrack->GetStartVertexMomentumZ() / P_MC);
        auto tanl_MC = thisTrack->Pz() / thisTrack->GetPt();
        auto pdgcode_MC = thisTrack->GetPdgCode();

        int Q_MC;
        if (TDatabasePDG::Instance()->GetParticle(pdgcode_MC)) {
          Q_MC =
            TDatabasePDG::Instance()->GetParticle(pdgcode_MC)->Charge() / 3;
        }

        else {
          iTrack++;
          continue;
          Q_MC = 0;
          // std::cout << " => pdgcode ERROR " << Q_MC <<  "\n";
        }
        auto invQPt_MC = 1.0 * Q_MC / thisTrack->GetPt();

        trackFCT.propagateToZhelix(vz_MC, field_z);

        auto Q_fit = trackFCT.getCharge();
        auto dx = trackFCT.getX() - vx_MC;
        auto dy = trackFCT.getY() - vy_MC;
        auto d_eta = trackFCT.getEta() - eta_MC;
        auto d_tanl = trackFCT.getTanl() - tanl_MC;
        auto Pt_fit = trackFCT.getPt();
        auto invQPt_Fit = trackFCT.getInvQPt();
        auto invQPt_seed = trackFCT.getInvQPtSeed();
        auto d_invQPt = Q_fit / Pt_fit - Q_MC / Pt_MC;
        auto d_invQPtSeed = invQPt_seed - Q_fit / Pt_fit;
        auto P_fit = trackFCT.getP();
        auto P_res = P_fit / P_MC;
        auto Pt_res = Pt_fit / Pt_MC;
        auto d_Phi = trackFCT.getPhi() - phi_MC;
        auto d_Charge = Q_fit - Q_MC;
        auto trackChi2 = trackFCT.getTrackChi2();

        TH3Histos[kFCTTrackDeltaXVertexPtEta]->Fill(Pt_MC, std::abs(eta_MC), 1e4 * dx);
        TH3Histos[kFCTTrackDeltaYVertexPtEta]->Fill(Pt_MC, std::abs(eta_MC), 1e4 * dy);
        TH3Histos[kFCTTrackPtResolutionPtEta]->Fill(Pt_MC, std::abs(eta_MC), (Pt_fit - Pt_MC) / Pt_MC);
        TH3Histos[kFCTTrackInvQPtPullPtEta]->Fill(Pt_MC, std::abs(eta_MC), d_invQPt / sqrt(trackFCT.getCovariances()(4, 4)));
        TH3Histos[kFCTTrackInvPtResolutionPtEta]->Fill(Pt_MC, std::abs(eta_MC), (1.0 / Pt_fit - 1.0 / Pt_MC) * Pt_MC);
        TH3Histos[kFCTTrackInvQPtResolutionPtEta]->Fill(Pt_MC, std::abs(eta_MC), (invQPt_Fit - invQPt_MC) / invQPt_MC);

        TH1Histos[kFCTTracksP]->Fill(trackFCT.getP());
        TH1Histos[kFCTTrackDeltaTanl]->Fill(d_tanl);
        TH1Histos[kFCTTrackDeltaPhi]->Fill(d_Phi);
        TH1Histos[kFCTTrackDeltaInvQPt]->Fill(d_invQPt);
        TH1Histos[kFCTTrackDeltaInvQPtSeed]->Fill(d_invQPtSeed);
        TH1Histos[kFCTTrackDeltaPhiDeg]->Fill(TMath::RadToDeg() * d_Phi);
        TH1Histos[kFCTTrackDeltaX]->Fill(10. * dx);

        TH1Histos[kFCTTrackDeltaXErr]->Fill(
          dx / sqrt(trackFCT.getCovariances()(0, 0)));
        TH1Histos[kFCTTrackDeltaYErr]->Fill(
          dy / sqrt(trackFCT.getCovariances()(1, 1)));
        TH1Histos[kFCTTrackDeltaPhiErr]->Fill(
          d_Phi / sqrt(trackFCT.getCovariances()(2, 2)));
        TH1Histos[kFCTTrackDeltaTanLErr]->Fill(
          d_tanl / sqrt(trackFCT.getCovariances()(3, 3)));
        TH1Histos[kFCTTrackDeltainvQPtErr]->Fill(
          d_invQPt / sqrt(trackFCT.getCovariances()(4, 4)));

        DeltaX_Profile->Fill(Pt_MC, dx * 1e4);
        TH1Histos[kFCTTrackDeltaY]->Fill(10 * dy);
        TH1Histos[kFCTTrackR]->Fill(10.0 * sqrt(dx * dx + dy * dy));
        TH1Histos[kFCTTrackQ]->Fill(d_Charge);
        TH1Histos[kFCTTrackChi2]->Fill(trackChi2);
        TH2Histos[kFCTTrackDeltaXYVertex]->Fill(10.0 * dx, 10.0 * dy);
        TH2Histos[kFCTTrackQPRec_MC]->Fill(P_MC * Q_MC, P_fit * Q_fit);
        TH2Histos[kFCTTrackPtResolution]->Fill(Pt_MC, (Pt_fit - Pt_MC) / Pt_MC);
        PtRes_Profile->Fill(Pt_MC, Pt_fit / Pt_MC);
        TH2Histos[kFCTTrackInvPtResolution]->Fill(
          Pt_MC, (1.0 / Pt_fit - 1.0 / Pt_MC) * Pt_MC);

        // MC histos
        TH1Histos[kMCTrackspT]->Fill(Pt_MC);
        TH1Histos[kMCTracksp]->Fill(P_MC);
        TH1Histos[kMCTrackEta]->Fill(eta_MC);
        TH2Histos[kMCTracksEtaZ]->Fill(vz_MC, eta_MC);

        // Differential histos

        if (InnerRegion(tanl_MC)) {
          TH2Histos[kFCTTrackPtResolutionInner]->Fill(Pt_MC, (Pt_fit - Pt_MC) / Pt_MC);
          TH2Histos[kFCTTrackInvPtResolutionInner]->Fill(
            Pt_MC, (1.0 / Pt_fit - 1.0 / Pt_MC) * Pt_MC);
        }

        if (OuterRegion(tanl_MC)) {
          TH2Histos[kFCTTrackPtResolutionOuter]->Fill(Pt_MC, (Pt_fit - Pt_MC) / Pt_MC);
          TH2Histos[kFCTTrackInvPtResolutionOuter]->Fill(
            Pt_MC, (1.0 / Pt_fit - 1.0 / Pt_MC) * Pt_MC);
        }

        if (Pt_MC <= 1.0) {
          TH2Histos[kFCTTrackDeltaXYVertex0_1]->Fill(10.0 * dx, 10.0 * dy);
          TH1Histos[kFCTTrackDeltaTanl0_1]->Fill(d_tanl);
          TH1Histos[kFCTTrackDeltaPhi0_1]->Fill(d_Phi);
          TH1Histos[kFCTTrackDeltaPhiDeg0_1]->Fill(TMath::RadToDeg() * d_Phi);
          TH1Histos[kFCTTrackDeltaX0_1]->Fill(10.0 * dx);
          TH1Histos[kFCTTrackQ0_1]->Fill(d_Charge);
          d_Charge ? nChargeMiss0_1++ : nChargeMatch0_1++;
        }

        if (InnerBorder(tanl_MC) and pt_1(Pt_MC)) {
          TH1Histos[kFCTTrackXPull1_innerBorder]->Fill(dx / sqrt(trackFCT.getCovariances()(0, 0)));
          TH1Histos[kFCTTrackYPull1_innerBorder]->Fill(dy / sqrt(trackFCT.getCovariances()(1, 1)));
          TH1Histos[kFCTTrackPhiPull1_innerBorder]->Fill(d_Phi / sqrt(trackFCT.getCovariances()(2, 2)));
          TH1Histos[kFCTTrackTanlPull1_innerBorder]->Fill(d_tanl / sqrt(trackFCT.getCovariances()(3, 3)));
          TH1Histos[kFCTTrackInvQPtPull1_innerBorder]->Fill(d_invQPt / sqrt(trackFCT.getCovariances()(4, 4)));
        }

        if (InnerBorder(tanl_MC) and pt_4(Pt_MC)) {
          TH1Histos[kFCTTrackXPull4_innerBorder]->Fill(dx / sqrt(trackFCT.getCovariances()(0, 0)));
          TH1Histos[kFCTTrackYPull4_innerBorder]->Fill(dy / sqrt(trackFCT.getCovariances()(1, 1)));
          TH1Histos[kFCTTrackPhiPull4_innerBorder]->Fill(d_Phi / sqrt(trackFCT.getCovariances()(2, 2)));
          TH1Histos[kFCTTrackTanlPull4_innerBorder]->Fill(d_tanl / sqrt(trackFCT.getCovariances()(3, 3)));
          TH1Histos[kFCTTrackInvQPtPull4_innerBorder]->Fill(d_invQPt / sqrt(trackFCT.getCovariances()(4, 4)));
        }

        if (OuterBorder(tanl_MC) and pt_1(Pt_MC)) {
          TH1Histos[kFCTTrackXPull1_OuterBorder]->Fill(dx / sqrt(trackFCT.getCovariances()(0, 0)));
          TH1Histos[kFCTTrackYPull1_OuterBorder]->Fill(dy / sqrt(trackFCT.getCovariances()(1, 1)));
          TH1Histos[kFCTTrackPhiPull1_OuterBorder]->Fill(d_Phi / sqrt(trackFCT.getCovariances()(2, 2)));
          TH1Histos[kFCTTrackTanlPull1_OuterBorder]->Fill(d_tanl / sqrt(trackFCT.getCovariances()(3, 3)));
          TH1Histos[kFCTTrackInvQPtPull1_OuterBorder]->Fill(d_invQPt / sqrt(trackFCT.getCovariances()(4, 4)));
        }

        if (OuterBorder(tanl_MC) and pt_4(Pt_MC)) {
          TH1Histos[kFCTTrackXPull4_OuterBorder]->Fill(dx / sqrt(trackFCT.getCovariances()(0, 0)));
          TH1Histos[kFCTTrackYPull4_OuterBorder]->Fill(dy / sqrt(trackFCT.getCovariances()(1, 1)));
          TH1Histos[kFCTTrackPhiPull4_OuterBorder]->Fill(d_Phi / sqrt(trackFCT.getCovariances()(2, 2)));
          TH1Histos[kFCTTrackTanlPull4_OuterBorder]->Fill(d_tanl / sqrt(trackFCT.getCovariances()(3, 3)));
          TH1Histos[kFCTTrackInvQPtPull4_OuterBorder]->Fill(d_invQPt / sqrt(trackFCT.getCovariances()(4, 4)));
        }

        if (Pt_MC > 1.0 and Pt_MC <= 4) {
          TH2Histos[kFCTTrackDeltaXYVertex1_4]->Fill(10.0 * dx, 10.0 * dy);
          TH1Histos[kFCTTrackDeltaTanl1_4]->Fill(d_tanl);
          TH1Histos[kFCTTrackDeltaPhi1_4]->Fill(d_Phi);
          TH1Histos[kFCTTrackDeltaPhiDeg1_4]->Fill(TMath::RadToDeg() * d_Phi);
          TH1Histos[kFCTTrackDeltaX1_4]->Fill(10.0 * dx);
          TH1Histos[kFCTTrackQ1_4]->Fill(d_Charge);
          d_Charge ? nChargeMiss1_4++ : nChargeMatch1_4++;
        }
        if (Pt_MC > 4.0) {
          TH2Histos[kFCTTrackDeltaXYVertex4plus]->Fill(10.0 * dx, 10.0 * dy);
          TH1Histos[kFCTTrackDeltaTanl4plus]->Fill(d_tanl);
          TH1Histos[kFCTTrackDeltaPhi4plus]->Fill(d_Phi);
          TH1Histos[kFCTTrackDeltaPhiDeg4plus]->Fill(TMath::RadToDeg() * d_Phi);
          TH1Histos[kFCTTrackDeltaX4plus]->Fill(10.0 * dx);
          TH1Histos[kFCTTrackQ4plus]->Fill(d_Charge);
          d_Charge ? nChargeMiss4plus++ : nChargeMatch4plus++;
        }

        d_Charge ? nChargeMiss++ : nChargeMatch++;
        qMatchEff->Fill(!d_Charge, Pt_MC);
      }
      iTrack++;
    } // Loop on TracksFCT
    totalTracks += iTrack;
  } // Loop over events

  // Customize histograms
  TH1Histos[kFCTTrackQ]->SetTitle(
    Form("nChargeMatch = %d (%.2f%%)", nChargeMatch,
         100. * nChargeMatch / (nChargeMiss + nChargeMatch)));
  TH1Histos[kFCTTrackQ0_1]->SetTitle(
    Form("nChargeMatch = %d (%.2f%%)", nChargeMatch0_1,
         100. * nChargeMatch0_1 / (nChargeMiss0_1 + nChargeMatch0_1)));
  TH1Histos[kFCTTrackQ1_4]->SetTitle(
    Form("nChargeMatch = %d (%.2f%%)", nChargeMatch1_4,
         100. * nChargeMatch1_4 / (nChargeMiss1_4 + nChargeMatch1_4)));
  TH1Histos[kFCTTrackQ4plus]->SetTitle(
    Form("nChargeMatch = %d (%.2f%%)", nChargeMatch4plus,
         100. * nChargeMatch4plus / (nChargeMiss4plus + nChargeMatch4plus)));

  qMatchEff->SetTitle(Form("Charge match = %.2f%%",
                           100. * nChargeMatch / (nChargeMiss + nChargeMatch)));

  // Remove stat boxes
  TH2Histos[kFCTTrackQPRec_MC]->SetStats(0);
  TH2Histos[kFCTTrackPtResolution]->SetStats(0);
  TH2Histos[kFCTTrackPtResolutionInner]->SetStats(0);
  TH2Histos[kFCTTrackPtResolutionOuter]->SetStats(0);
  TH2Histos[kFCTTrackInvPtResolution]->SetStats(0);
  TH2Histos[kFCTTrackInvPtResolutionInner]->SetStats(0);
  TH2Histos[kFCTTrackInvPtResolutionOuter]->SetStats(0);
  TH2Histos[kMCTracksEtaZ]->SetStats(0);
  PtRes_Profile->SetStats(0);
  DeltaX_Profile->SetStats(0);
  TH1Histos[kFCTTrackQ]->SetStats(0);

  // Fit Slices: Pt resolution
  FitSlicesy(*TH2Histos[kFCTTrackInvPtResolution], *TH2Histos[kFCTTrackQPRec_MC], "E((1/pt_{fit} - 1.pt_{MC}) / (1/pt_{MC}))", "1/Pt Resolution");
  FitSlicesy(*TH2Histos[kFCTTrackInvPtResolutionInner], *TH2Histos[kFCTTrackQPRec_MC], "E((1/pt_{fit} - 1.pt_{MC}) / (1/pt_{MC}))", "\\text{1/Pt Resolution }(3.5 < \\eta < 3.6)");
  FitSlicesy(*TH2Histos[kFCTTrackInvPtResolutionOuter], *TH2Histos[kFCTTrackQPRec_MC], "E((1/pt_{fit} - 1.pt_{MC}) / (1/pt_{MC}))", "\\text{1/Pt Resolution }(2.8 < \\eta < 2.9 )");
  FitSlicesy(*TH2Histos[kFCTTrackPtResolution], *TH2Histos[kFCTTrackQPRec_MC], "E((pt_{fit} - pt_{MC}) / pt_{MC})", "Pt Resolution");
  FitSlicesy(*TH2Histos[kFCTTrackPtResolutionInner], *TH2Histos[kFCTTrackQPRec_MC], "E((pt_{fit} - pt_{MC}) / pt_{MC})", "\\text{Pt Resolution }(3.5 < \\eta < 3.6)");
  FitSlicesy(*TH2Histos[kFCTTrackPtResolutionOuter], *TH2Histos[kFCTTrackQPRec_MC], "E((pt_{fit} - pt_{MC}) / pt_{MC})", "\\text{Pt Resolution }(2.8 < \\eta < 2.9 )");

  // sigmaX resolution Profile
  TH1D* DeltaX_Error = new TH1D();
  DeltaX_Error = DeltaX_Profile->ProjectionX("DeltaX_Error", "C=E");

  // pt resolution Profile
  TH1D* pt_resolution_from_profile = new TH1D();
  pt_resolution_from_profile = PtRes_Profile->ProjectionX("Pt Resolution from Profile", "C=E");

  // Summary Canvases
  auto param_resolution = summary_report_3x2(
    *TH2Histos[kFCTTrackDeltaXYVertex], *TH2Histos[kFCTTrackPtResolution],
    *PtRes_Profile, *DeltaX_Error, *TH2Histos[kFCTTrackQPRec_MC], *qMatchEff,
    "Param Summary", seed_cfg, 0, 0, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kFCTTrackDeltaXYVertex]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackPtResolution]->Integral() /
                     TH2Histos[kFCTTrackPtResolution]->GetEntries()),
    "-", "-",
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackQPRec_MC]->Integral() /
                     TH2Histos[kFCTTrackQPRec_MC]->GetEntries()),
    "-");

  auto covariances_summary = summary_report_3x2(
    *TH1Histos[kFCTTrackDeltaXErr], *TH1Histos[kFCTTrackDeltaPhiErr],
    *TH1Histos[kFCTTrackDeltainvQPtErr], *TH1Histos[kFCTTrackDeltaYErr],
    *TH1Histos[kFCTTrackDeltaTanLErr], *TH2Histos[kFCTTrackQPRec_MC],
    "Covariances Summary", seed_cfg, 1, 1, 1, 1, 1, 0,
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaXErr]->Integral() /
                     TH1Histos[kFCTTrackDeltaXErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaPhiErr]->Integral() /
                     TH1Histos[kFCTTrackDeltaPhiErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltainvQPtErr]->Integral() /
                     TH1Histos[kFCTTrackDeltainvQPtErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaYErr]->Integral() /
                     TH1Histos[kFCTTrackDeltaYErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaTanLErr]->Integral() /
                     TH1Histos[kFCTTrackDeltaTanLErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackQPRec_MC]->Integral() /
                     TH2Histos[kFCTTrackQPRec_MC]->GetEntries()));

  auto long_summary = summary_report_3x3(
    *TH2Histos[kFCTTrackDeltaXYVertex], *TH1Histos[kFCTTrackDeltaXErr],
    *TH1Histos[kFCTTrackDeltaYErr], *DeltaX_Error,
    *TH2Histos[kFCTTrackQPRec_MC], *TH1Histos[kFCTTrackDeltaPhiErr],
    *qMatchEff, *TH1Histos[kFCTTrackDeltainvQPtErr],
    *TH1Histos[kFCTTrackDeltaTanLErr], "Summary3x3", seed_cfg, 0, 1, 1, 0, 0,
    1, 0, 1, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kFCTTrackDeltaXYVertex]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaXErr]->Integral() /
                     TH1Histos[kFCTTrackDeltaXErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaYErr]->Integral() /
                     TH1Histos[kFCTTrackDeltaYErr]->GetEntries()),
    "-",
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackQPRec_MC]->Integral() /
                     TH2Histos[kFCTTrackQPRec_MC]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaPhiErr]->Integral() /
                     TH1Histos[kFCTTrackDeltaPhiErr]->GetEntries()),
    "-",
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltainvQPtErr]->Integral() /
                     TH1Histos[kFCTTrackDeltainvQPtErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaTanLErr]->Integral() /
                     TH1Histos[kFCTTrackDeltaTanLErr]->GetEntries()));

  auto param_summary_diff_pt = summary_report_3x3(
    *TH1Histos[kFCTTrackDeltaX0_1], *TH1Histos[kFCTTrackDeltaTanl0_1],
    *TH1Histos[kFCTTrackDeltaPhiDeg0_1], *TH1Histos[kFCTTrackDeltaX1_4],
    *TH1Histos[kFCTTrackDeltaTanl1_4], *TH1Histos[kFCTTrackDeltaPhiDeg1_4],
    *TH1Histos[kFCTTrackDeltaX4plus], *TH1Histos[kFCTTrackDeltaTanl4plus],
    *TH1Histos[kFCTTrackDeltaPhiDeg4plus], "ParamSummaryVsPt", seed_cfg, 1, 1,
    1, 1, 1, 1, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaX0_1]->Integral() /
                     TH1Histos[kFCTTrackDeltaX0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaTanl0_1]->Integral() /
                     TH1Histos[kFCTTrackDeltaTanl0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaPhiDeg0_1]->Integral() /
                     TH1Histos[kFCTTrackDeltaPhiDeg0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaX1_4]->Integral() /
                     TH1Histos[kFCTTrackDeltaX1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaTanl1_4]->Integral() /
                     TH1Histos[kFCTTrackDeltaTanl1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaPhiDeg1_4]->Integral() /
                     TH1Histos[kFCTTrackDeltaPhiDeg1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaX4plus]->Integral() /
                     TH1Histos[kFCTTrackDeltaX4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaTanl4plus]->Integral() /
                     TH1Histos[kFCTTrackDeltaTanl4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaPhiDeg4plus]->Integral() /
                     TH1Histos[kFCTTrackDeltaPhiDeg4plus]->GetEntries()));

  auto pt_resolution = summary_report(
    *TH2Histos[kFCTTrackPtResolution], *TH2Histos[kFCTTrackQPRec_MC],
    *PtRes_Profile, *qMatchEff, "Pt Summary", seed_cfg, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackPtResolution]->Integral() /
                     TH2Histos[kFCTTrackPtResolution]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackQPRec_MC]->Integral() /
                     TH2Histos[kFCTTrackQPRec_MC]->GetEntries()));

  auto pt_resolution_2 = summary_report(
    *TH2Histos[kFCTTrackPtResolution],
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kFCTTrackPtResolutionInner]->GetName()) +
       std::string("_2"))
        .c_str()),
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kFCTTrackPtResolution]->GetName()) +
       std::string("_1"))
        .c_str()),
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kFCTTrackPtResolutionOuter]->GetName()) +
       std::string("_2"))
        .c_str()),
    "Pt Resolution Summary", seed_cfg, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackPtResolution]->Integral() /
                     TH2Histos[kFCTTrackPtResolution]->GetEntries()));

  auto invpt_resolution = summary_report(
    *TH2Histos[kFCTTrackInvPtResolution], *TH2Histos[kFCTTrackQPRec_MC],
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kFCTTrackInvPtResolution]->GetName()) +
       std::string("_1"))
        .c_str()),
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kFCTTrackInvPtResolution]->GetName()) +
       std::string("_2"))
        .c_str()),
    "InvPt Summary", seed_cfg, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackInvPtResolution]->Integral() /
                     TH2Histos[kFCTTrackInvPtResolution]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackQPRec_MC]->Integral() /
                     TH2Histos[kFCTTrackQPRec_MC]->GetEntries()));

  auto invpt_resolution_2 = summary_report(
    *TH2Histos[kFCTTrackInvPtResolution],
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kFCTTrackInvPtResolutionInner]->GetName()) +
       std::string("_2"))
        .c_str()),
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kFCTTrackInvPtResolution]->GetName()) +
       std::string("_1"))
        .c_str()),
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kFCTTrackInvPtResolutionOuter]->GetName()) +
       std::string("_2"))
        .c_str()),
    "InvPt Resolution Summary", seed_cfg, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackInvPtResolution]->Integral() /
                     TH2Histos[kFCTTrackInvPtResolution]->GetEntries()));

  auto vertexing_resolution = summary_report(
    *TH2Histos[kFCTTrackDeltaXYVertex], *TH1Histos[kFCTTrackDeltaX],
    *DeltaX_Error, *TH1Histos[kFCTTrackDeltaPhiDeg], "Vertexing Summary",
    seed_cfg, 0, 1, 0, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kFCTTrackDeltaXYVertex]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaX]->Integral() /
                     TH1Histos[kFCTTrackDeltaX]->GetEntries()),
    Form("-"),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaPhiDeg]->Integral() /
                     TH1Histos[kFCTTrackDeltaPhiDeg]->GetEntries()));

  auto vertexing_resolution0_1 = summary_report(
    *TH2Histos[kFCTTrackDeltaXYVertex0_1], *TH1Histos[kFCTTrackDeltaX0_1],
    *TH1Histos[kFCTTrackDeltaTanl0_1], *TH1Histos[kFCTTrackDeltaPhiDeg0_1],
    "Vertexing Summary pt < 1", seed_cfg, 0, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackDeltaXYVertex0_1]->Integral() /
                     TH2Histos[kFCTTrackDeltaXYVertex0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaX0_1]->Integral() /
                     TH1Histos[kFCTTrackDeltaX0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaTanl0_1]->Integral() /
                     TH1Histos[kFCTTrackDeltaTanl0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaPhiDeg0_1]->Integral() /
                     TH1Histos[kFCTTrackDeltaPhiDeg0_1]->GetEntries()));

  auto vertexing_resolution1_4 = summary_report(
    *TH2Histos[kFCTTrackDeltaXYVertex1_4], *TH1Histos[kFCTTrackDeltaX1_4],
    *TH1Histos[kFCTTrackDeltaTanl1_4], *TH1Histos[kFCTTrackDeltaPhiDeg1_4],
    "Vertexing Summary 1 < p_t < 4", seed_cfg, 0, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kFCTTrackDeltaXYVertex1_4]->Integral() /
                     TH2Histos[kFCTTrackDeltaXYVertex1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaX1_4]->Integral() /
                     TH1Histos[kFCTTrackDeltaX1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaTanl1_4]->Integral() /
                     TH1Histos[kFCTTrackDeltaTanl1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaPhiDeg1_4]->Integral() /
                     TH1Histos[kFCTTrackDeltaPhiDeg1_4]->GetEntries()));

  auto vertexing_resolution4plus = summary_report(
    *TH2Histos[kFCTTrackDeltaXYVertex4plus], *TH1Histos[kFCTTrackDeltaX4plus],
    *TH1Histos[kFCTTrackDeltaTanl4plus],
    *TH1Histos[kFCTTrackDeltaPhiDeg4plus], "Vertexing Summary p_t > 4",
    seed_cfg, 0, 1, 1, 1,
    Form("%.2f%%", 100.0 *
                     TH2Histos[kFCTTrackDeltaXYVertex4plus]->Integral() /
                     TH2Histos[kFCTTrackDeltaXYVertex4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaX4plus]->Integral() /
                     TH1Histos[kFCTTrackDeltaX4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaTanl4plus]->Integral() /
                     TH1Histos[kFCTTrackDeltaTanl4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackDeltaPhiDeg4plus]->Integral() /
                     TH1Histos[kFCTTrackDeltaPhiDeg4plus]->GetEntries()));
  // Pulls summaries

  auto XpullSummary = summary_report(
    *TH1Histos[kFCTTrackXPull1_innerBorder], *TH1Histos[kFCTTrackXPull4_innerBorder],
    *TH1Histos[kFCTTrackXPull1_OuterBorder], *TH1Histos[kFCTTrackXPull4_OuterBorder],
    "XpullSummary", seed_cfg, 1, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackXPull1_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackXPull1_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackXPull4_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackXPull4_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackXPull1_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackXPull1_OuterBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackXPull4_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackXPull4_OuterBorder]->GetEntries()));
  //
  auto YpullSummary = summary_report(
    *TH1Histos[kFCTTrackYPull1_innerBorder], *TH1Histos[kFCTTrackYPull4_innerBorder],
    *TH1Histos[kFCTTrackYPull1_OuterBorder], *TH1Histos[kFCTTrackYPull4_OuterBorder],
    "YpullSummary", seed_cfg, 1, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackYPull1_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackYPull1_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackYPull4_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackYPull4_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackYPull1_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackYPull1_OuterBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackYPull4_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackYPull4_OuterBorder]->GetEntries()));
  //
  auto PhipullSummary = summary_report(
    *TH1Histos[kFCTTrackPhiPull1_innerBorder], *TH1Histos[kFCTTrackPhiPull4_innerBorder],
    *TH1Histos[kFCTTrackPhiPull1_OuterBorder], *TH1Histos[kFCTTrackPhiPull4_OuterBorder],
    "PhipullSummary", seed_cfg, 1, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackPhiPull1_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackPhiPull1_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackPhiPull4_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackPhiPull4_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackPhiPull1_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackPhiPull1_OuterBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackPhiPull4_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackPhiPull4_OuterBorder]->GetEntries()));
  //
  auto TanlpullSummary = summary_report(
    *TH1Histos[kFCTTrackTanlPull1_innerBorder], *TH1Histos[kFCTTrackTanlPull4_innerBorder],
    *TH1Histos[kFCTTrackTanlPull1_OuterBorder], *TH1Histos[kFCTTrackTanlPull4_OuterBorder],
    "TanlpullSummary", seed_cfg, 1, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackTanlPull1_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackTanlPull1_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackTanlPull4_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackTanlPull4_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackTanlPull1_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackTanlPull1_OuterBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackTanlPull4_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackTanlPull4_OuterBorder]->GetEntries()));
  //
  auto InvQPtpullSummary = summary_report(
    *TH1Histos[kFCTTrackInvQPtPull1_innerBorder], *TH1Histos[kFCTTrackInvQPtPull4_innerBorder],
    *TH1Histos[kFCTTrackInvQPtPull1_OuterBorder], *TH1Histos[kFCTTrackInvQPtPull4_OuterBorder],
    "InvQPtpullSummary", seed_cfg, 1, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackInvQPtPull1_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackInvQPtPull1_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackInvQPtPull4_innerBorder]->Integral() /
                     TH1Histos[kFCTTrackInvQPtPull4_innerBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackInvQPtPull1_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackInvQPtPull1_OuterBorder]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kFCTTrackInvQPtPull4_OuterBorder]->Integral() /
                     TH1Histos[kFCTTrackInvQPtPull4_OuterBorder]->GetEntries()));

  // Write histograms to file and export images

  outFile.mkdir("MoreHistos");
  outFile.cd("MoreHistos");

  for (auto& h : TH3Histos) {
    h->Write();
    if (EXPORT_HISTOS_IMAGES)
      exportHisto(*h);
  }

  for (auto& h : TH2Histos) {
    h->Write();
    if (EXPORT_HISTOS_IMAGES)
      exportHisto(*h);
  }

  for (auto& h : TH1Histos) {
    h->Write();
    if (EXPORT_HISTOS_IMAGES)
      exportHisto(*h);
  }

  PtRes_Profile->Write();
  DeltaX_Profile->Write();
  DeltaX_Error->Write();
  pt_resolution_from_profile->Write();
  qMatchEff->Write();
  outFile.Close();

  std::cout << "---------------------------------------------------"
            << std::endl;
  std::cout << "-------------   Fitting Summary   -----------------"
            << std::endl;
  std::cout << "---------------------------------------------------"
            << std::endl;
  std::cout << " P_mean = " << TH1Histos[kFCTTracksP]->GetMean() << std::endl;
  std::cout << " P_StdDev = " << TH1Histos[kFCTTracksP]->GetStdDev()
            << std::endl;
  std::cout << " Tanl_mean = " << TH1Histos[kFCTTrackDeltaTanl]->GetMean()
            << std::endl;
  std::cout << " Tanl_StdDev = " << TH1Histos[kFCTTrackDeltaTanl]->GetStdDev()
            << std::endl;
  std::cout << " Tanl_StdDev(pt<1) = "
            << TH1Histos[kFCTTrackDeltaTanl0_1]->GetStdDev() << std::endl;
  std::cout << " Tanl_StdDev(1<pt<4) = "
            << TH1Histos[kFCTTrackDeltaTanl1_4]->GetStdDev() << std::endl;
  std::cout << " Tanl_StdDev(pt>4) = "
            << TH1Histos[kFCTTrackDeltaTanl4plus]->GetStdDev() << std::endl;
  std::cout << " Phi_mean = " << TH1Histos[kFCTTrackDeltaPhi]->GetMean()
            << std::endl;
  std::cout << " Phi_StdDev = " << TH1Histos[kFCTTrackDeltaPhi]->GetStdDev()
            << std::endl;
  std::cout << " Phi_StdDev(pt<1) = "
            << TH1Histos[kFCTTrackDeltaPhi0_1]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDev(1<pt<4) = "
            << TH1Histos[kFCTTrackDeltaPhi1_4]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDev(pt>4) = "
            << TH1Histos[kFCTTrackDeltaPhi4plus]->GetStdDev() << std::endl;
  std::cout << " Phi_meanDeg = " << TH1Histos[kFCTTrackDeltaPhiDeg]->GetMean()
            << std::endl;
  std::cout << " Phi_StdDevDeg = "
            << TH1Histos[kFCTTrackDeltaPhiDeg]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDevDeg(pt<1) = "
            << TH1Histos[kFCTTrackDeltaPhiDeg0_1]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDevDeg(1<pt<4) = "
            << TH1Histos[kFCTTrackDeltaPhiDeg1_4]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDevDeg(pt>4) = "
            << TH1Histos[kFCTTrackDeltaPhiDeg4plus]->GetStdDev() << std::endl;
  std::cout << " DeltaX_mean = " << TH1Histos[kFCTTrackDeltaX]->GetMean()
            << std::endl;
  std::cout << " DeltaX_StdDev = " << TH1Histos[kFCTTrackDeltaX]->GetStdDev()
            << std::endl;
  std::cout << " DeltaX_StdDev(pt<1) = "
            << TH1Histos[kFCTTrackDeltaX0_1]->GetStdDev() << std::endl;
  std::cout << " DeltaX_StdDev(1<pt<4) = "
            << TH1Histos[kFCTTrackDeltaX1_4]->GetStdDev() << std::endl;
  std::cout << " DeltaX_StdDev(pt>4) = "
            << TH1Histos[kFCTTrackDeltaX4plus]->GetStdDev() << std::endl;
  std::cout << " DeltaY_mean = " << TH1Histos[kFCTTrackDeltaY]->GetMean()
            << std::endl;
  std::cout << " DeltaY_StdDev = " << TH1Histos[kFCTTrackDeltaY]->GetStdDev()
            << std::endl;
  std::cout << " R_mean = " << TH1Histos[kFCTTrackR]->GetMean() << std::endl;
  std::cout << " R_StdDev = " << TH1Histos[kFCTTrackR]->GetStdDev()
            << std::endl;
  std::cout << " Charge_mean = " << TH1Histos[kFCTTrackDeltaY]->GetMean()
            << std::endl;
  std::cout << " nChargeMatch = " << nChargeMatch << " ("
            << 100. * nChargeMatch / (nChargeMiss + nChargeMatch) << "%)"
            << std::endl;
  std::cout << "---------------------------------------------------"
            << std::endl;

  return 0;
}
