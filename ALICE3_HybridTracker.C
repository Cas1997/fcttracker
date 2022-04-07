#include "FCTTrack.h"
#include "TFile.h"
#include "TTree.h"
#include "TrackFitter.h"
#include "fcttools/MagField.C"

#include <map>

#pragma link C++ class o2::fct::FCTTrack + ;
#pragma link C++ class std::vector < o2::fct::FCTTrack> + ;
#pragma link C++ class o2::fct::FCTTrackExt + ;
#pragma link C++ class std::vector < o2::fct::FCTTrackExt> + ;

using o2::fct::FCTTrack;
using o2::itsmft::Hit;
using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;

bool DEBUG_VERBOSE = true;
auto PrintTracksTotal = 10;
//std::vector<Double_t> layersX2X0 = {};

std::vector<Double_t> loadx2X0fromFile(std::string configFileName);
void setSeedCovariances(FCTTrack& track);

void ALICE3_HybridTracker(Double_t clResolution = 10.0e-4, bool verbose = true) // clResolution = 8.44e-4
{
  DEBUG_VERBOSE = verbose;

  o2::fct::TrackFitter fitter;
  fitter.mVerbose = DEBUG_VERBOSE;
  auto field_z = getZField(0, 0, 0); // Get field at Center of ALICE
  fitter.setBz(field_z);
  fitter.setLayersx2X0(loadx2X0fromFile("FCT_layout.cfg"));

  std::string FCThitfile = "o2sim_HitsFCT.root";
  std::string TRKhitfile = "o2sim_HitsTRK.root";

  std::string ALICE3HibridTracksfile = "ALICE3tracks.root";

  vector<map<int, FCTTrack>> allFwdTracks;
  vector<vector<o2::fct::FCTTrackExt>> recoFwdTracks;
  vector<vector<Int_t>> recoFwdTrackIDs;

  TFile* FCTHitFileIn = new TFile(FCThitfile.c_str());
  TTree* FCThitTree = (TTree*)FCTHitFileIn->Get("o2sim");

  TFile* TRKHitFileIn = new TFile(TRKhitfile.c_str());
  TTree* TRKhitTree = (TTree*)TRKHitFileIn->Get("o2sim");

  vector<Hit>* FCThit = nullptr;
  FCThitTree->SetBranchAddress("FCTHit", &FCThit);

  vector<Hit>* TRKhit = nullptr;
  TRKhitTree->SetBranchAddress("TRKHit", &TRKhit);

  Int_t nEvents = FCThitTree->GetEntries();

  allFwdTracks.resize(nEvents);
  recoFwdTracks.resize(nEvents);
  recoFwdTrackIDs.resize(nEvents);

  auto nPrintTracks = 0;

  for (Int_t event = 0; event < nEvents; event++) {

    FCThitTree->GetEntry(event);
    TRKhitTree->GetEntry(event);

    Int_t nFCTHits = FCThit->size();
    Int_t nTRKHits = TRKhit->size();

    if (1) {
      std::cout << " Building tracks at event " << event << std::endl;
    }

    auto& FwdTracksMap = allFwdTracks.at(event);
    for (Int_t iHit = 0; iHit < nFCTHits; iHit++) { // Attach FCT hits to tracks
      Hit* thisHit = &(*FCThit)[iHit];
      Int_t myID = thisHit->GetTrackID();
      FwdTracksMap[myID] = FwdTracksMap[myID];
      FwdTracksMap[myID].addHit(*thisHit, iHit, clResolution, true);
      FwdTracksMap[myID].setEventID(event);
      FwdTracksMap[myID].setTrackID(myID);
    }

    for (Int_t iHit = 0; iHit < nTRKHits; iHit++) { // Attach TRK hits to tracks
      Hit* thisHit = &(*TRKhit)[iHit];
      Int_t myID = thisHit->GetTrackID();
      FwdTracksMap[myID] = FwdTracksMap[myID];
      FwdTracksMap[myID].addTRKHit(*thisHit, iHit, clResolution, true);
      FwdTracksMap[myID].setEventID(event);
      FwdTracksMap[myID].setTrackID(myID);
    }
  }

  TFile* FCTTracksFileOut = new TFile(ALICE3HibridTracksfile.c_str(), "RECREATE");
  TTree* FCTTree = new TTree("o2sim", "Tree with FCT Tracks");
  vector<o2::fct::FCTTrackExt>* recoTracks;
  vector<Int_t>* recoTrackIDs;
  FCTTree->Branch("FCTTrack", &recoTracks);
  FCTTree->Branch("FCTTrackID", &recoTrackIDs);

  for (Int_t event = 0; event < nEvents; event++) { // Fit and save tracks
    if (1) {
      std::cout << " Fitting tracks at event " << event << std::endl;
    }
    recoTracks = &recoFwdTracks[event];
    recoTrackIDs = &recoFwdTrackIDs[event];
    auto& FwdTracksMap = allFwdTracks.at(event);

    auto iter = FwdTracksMap.begin();
    while (iter != FwdTracksMap.end()) {
      auto& trackID = iter->first;
      auto& track = iter->second;
      auto nHits = track.getNumberOfPoints();
      if (nHits > 6) {
        track.sort();
        fitter.initTrack(track);
        setSeedCovariances(track);
        fitter.fit(track);
        recoTracks->emplace_back(track);
        recoTrackIDs->emplace_back(trackID);
        if (fitter.mVerbose) {
          if (nPrintTracks > PrintTracksTotal)
            fitter.mVerbose = false;
        }
      }
      ++iter;
      nPrintTracks++;
    }

    FCTTree->Fill();
  }
  FCTTracksFileOut->Write();
}

std::vector<Double_t> loadx2X0fromFile(std::string configFileName = "FCT_layout.cfg")
{
  std::vector<Double_t> Layersx2X0;
  std::ifstream ifs(configFileName.c_str());
  if (!ifs.good()) {
    LOG(FATAL) << " Invalid FCTBase.configFile!";
  }
  std::string tempstr;
  Double_t z_layer, r_in, r_out, Layerx2X0;
  char delimiter;
  int layerNumber = 0;
  while (std::getline(ifs, tempstr)) {
    if (tempstr[0] == '#') {
      LOG(INFO) << " Comment: " << tempstr;
      continue;
    }
    LOG(INFO) << " Line: " << tempstr;
    std::istringstream iss(tempstr);
    iss >> z_layer;
    iss >> r_in;
    iss >> r_out;
    iss >> Layerx2X0;

    Layersx2X0.push_back(Layerx2X0);
    LOG(INFO) << " loadx2X0fromFile z =  " << z_layer << " ; x/X0 = " << Layerx2X0 << std::endl;
  }
  return Layersx2X0;
}

void setSeedCovariances(FCTTrack& track)
{

  auto tan_k = 10.0;
  auto q2pt_c = 10.;
  auto phi_c = 1 / 4.;

  SMatrix55Sym Covariances{}; ///< \brief Covariance matrix of track parameters
  Double_t q2ptcov = TMath::Max(std::abs(track.getInvQPt()), .5);
  Double_t tanlerr = TMath::Max(std::abs(track.getTanl()), .5);

  Covariances(0, 0) = 1;
  Covariances(1, 1) = 1;
  Covariances(2, 2) = phi_c * TMath::Pi() * TMath::Pi();
  Covariances(3, 3) = tan_k * tanlerr * tanlerr;
  Covariances(4, 4) = q2pt_c * q2ptcov * q2ptcov;
  track.setCovariances(Covariances);
}
