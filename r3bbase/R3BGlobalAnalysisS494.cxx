/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ------------------------------------------------------------
// -----                  R3BGlobalAnalysisS494                -----
// -----          Created April 13th 2016 by M.Heil       -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaHitData.h"
#include "R3BCalifaMappedData.h"

#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BLosMappedData.h"

#include "R3BBeamMonitorMappedData.h"

#include "R3BGlobalAnalysisS494.h"

#include "R3BSci8CalData.h"
#include "R3BSci8MappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BRoluCalData.h"
#include "R3BRoluMappedData.h"

#include "R3BPaddleCalData.h"

#include "R3BPspxCalData.h"
#include "R3BPspxMappedData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberMappedData.h"

#include "R3BMCTrack.h"
#include "R3BTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

#include "tracker_routines.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

R3BGlobalAnalysisS494::R3BGlobalAnalysisS494()
    : R3BGlobalAnalysisS494("GlobalAnalysis", 1)
{
}

R3BGlobalAnalysisS494::R3BGlobalAnalysisS494(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fCuts(0)
    , fGhost(0)
    , fPairs(1)
    , fSimu(0)
    , fGraphCuts(0)
    , fB(-1710)
    , fcut_chiX(5000)
    , fcut_chiY(5000)
    , fNEvents(0)
    , fMappedItemsCalifa(NULL)
    , fCalItemsCalifa(NULL)
    , fHitItemsCalifa(NULL)
{
}

R3BGlobalAnalysisS494::~R3BGlobalAnalysisS494() {}

InitStatus R3BGlobalAnalysisS494::Init()
{
    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    cout << "R3BGlobalAnalysisS494::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
    FairRunOnline* run = FairRunOnline::Instance();

    // Get objects for detectors on all levels
    fMCTrack = (TClonesArray*)mgr->GetObject("MCTrack");
    fTrack = (TClonesArray*)mgr->GetObject("Track");
    maxevent = mgr->CheckMaxEventNo();

    // get access to data of Califa
    fMappedItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaMappedData");
    if (!fMappedItemsCalifa)
        LOG(WARNING) << "GlobalAnalysis: CalifaMappedData not found";
    fCalItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaCrystalCalData");
    if (!fCalItemsCalifa)
        LOG(WARNING) << "GlobalAnalysis: CalifaCrystalCalData not found";
    fHitItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaHitData");
    if (!fHitItemsCalifa)
        LOG(WARNING) << "GlobalAnalysis: CalifaHitData not found";

    cout << "R3BGlobalAnalysisS494::Max num events: " << maxevent << endl;

    //------------------------------------------------------------------------
    // graphical cuts
    //------------------------------------------------------------------------
    if (fGraphCuts)
    {
        cout << "READING GRAPHICAL CUTS " << endl;
        cut_EHe_EC = NULL;

        TFile* fevse = TFile::Open("Cut_global_evse.root", "read");
        cut_EHe_EC = dynamic_cast<TCutG*>(fevse->Get("Cut_EHeEC"));

        cout << "GRAPHICAL CUTS ARE READ" << endl;
    }
    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    //    char strNameC[255];
    //    sprintf(strNameC, "Cave_C_position");
    //    fh_Cave_position = new TH2F(strNameC, "", 2100, -100., 2000., 1000, -500., 500.);

    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    //

    fh_px_He = new TH1F("px_He", " px He", 500, -500., 500);
    fh_px_He->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_He->GetYaxis()->SetTitle("counts");

    fh_py_He = new TH1F("py_He", " py He", 500, -500., 500);
    fh_py_He->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_He->GetYaxis()->SetTitle("counts");

    fh_pz_He = new TH1F("pz_He", " pz He", 1000, 2000., 5000.);
    fh_pz_He->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_He->GetYaxis()->SetTitle("counts");

    fh_p_He = new TH1F("p_He", " p He", 1000, 2000., 5000.);
    fh_p_He->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_He->GetYaxis()->SetTitle("counts");

    fh_px_C = new TH1F("px_C", " px C", 500, -500., 500);
    fh_px_C->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_C->GetYaxis()->SetTitle("counts");

    fh_py_C = new TH1F("py_C", " py C", 500, -500., 500);
    fh_py_C->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_C->GetYaxis()->SetTitle("counts");

    fh_pz_C = new TH1F("pz_C", " pz C", 2000, 10000., 16000.);
    fh_pz_C->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_C->GetYaxis()->SetTitle("counts");

    fh_p_C = new TH1F("p_C", " p C", 200, 10000., 16000.);
    fh_p_C->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_C->GetYaxis()->SetTitle("counts");

    fh_px_O = new TH1F("px_O", " px O", 500, -500., 500);
    fh_px_O->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_O->GetYaxis()->SetTitle("counts");

    fh_py_O = new TH1F("py_O", " py O", 500, -500., 500);
    fh_py_O->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_O->GetYaxis()->SetTitle("counts");

    fh_pz_O = new TH1F("pz_O", " pz O", 1000, 16000., 19000.);
    fh_pz_O->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_O->GetYaxis()->SetTitle("counts");

    fh_p_O = new TH1F("p_O", " p O", 1000, 16000., 19000.);
    fh_p_O->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_O->GetYaxis()->SetTitle("counts");

    fh_target_xy = new TH2F("target_xy", "target xy ", 500, -1, 1, 500, -1, 1);
    fh_target_xy->GetXaxis()->SetTitle("x / cm");
    fh_target_xy->GetYaxis()->SetTitle("y / cm");

    fh_px_He_nc = new TH1F("px_He_nc", " px He all chi2", 500, -500., 500);
    fh_px_He_nc->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_He_nc->GetYaxis()->SetTitle("counts");

    fh_py_He_nc = new TH1F("py_He_nc", " py He all chi2", 500, -500., 500);
    fh_py_He_nc->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_He_nc->GetYaxis()->SetTitle("counts");

    fh_pz_He_nc = new TH1F("pz_He_nc", " pz He all chi2", 1000, 2000., 5000.);
    fh_pz_He_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_He_nc->GetYaxis()->SetTitle("counts");

    fh_p_He_nc = new TH1F("p_He_nc", " p He all chi2", 1000, 2000., 5000.);
    fh_p_He_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_He_nc->GetYaxis()->SetTitle("counts");

    fh_px_C_nc = new TH1F("px_C_nc", " px C all chi2", 500, -500., 500);
    fh_px_C_nc->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_C_nc->GetYaxis()->SetTitle("counts");

    fh_py_C_nc = new TH1F("py_C_nc", " py C all chi2", 500, -500., 500);
    fh_py_C_nc->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_C_nc->GetYaxis()->SetTitle("counts");

    fh_pz_C_nc = new TH1F("pz_C_nc", " pz C all chi2", 2000, 10000., 16000.);
    fh_pz_C_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_C_nc->GetYaxis()->SetTitle("counts");

    fh_p_C_nc = new TH1F("p_C_nc", " p C all chi2", 2000, 10000., 16000.);
    fh_p_C_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_C_nc->GetYaxis()->SetTitle("counts");

    fh_px_O_nc = new TH1F("px_O_nc", " px O all chi2", 500, -500., 500);
    fh_px_O_nc->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_O_nc->GetYaxis()->SetTitle("counts");

    fh_py_O_nc = new TH1F("py_O_nc", " py O all chi2", 500, -500., 500);
    fh_py_O_nc->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_O_nc->GetYaxis()->SetTitle("counts");

    fh_pz_O_nc = new TH1F("pz_O_nc", " pz O all chi2", 1000, 16000., 19000.);
    fh_pz_O_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_O_nc->GetYaxis()->SetTitle("counts");

    fh_p_O_nc = new TH1F("p_O_nc", " p O all chi2", 1000, 16000., 19000.);
    fh_p_O_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_O_nc->GetYaxis()->SetTitle("counts");

    fh_psum_vs_r_nc = new TH2F("psum_vs_targetR_nc", "psum vs Rtarget all chi2", 250, 0, 2.5, 1000, 16000., 19000.);
    fh_psum_vs_r_nc->GetXaxis()->SetTitle("Radius / cm");
    fh_psum_vs_r_nc->GetYaxis()->SetTitle("psum / MeV/c");

    fh_px_vs_r = new TH2F("px_vs_targetR", "px vs Rtarget ", 250, 0, 2.5, 500, -500., 500.);
    fh_px_vs_r->GetXaxis()->SetTitle("Radius / cm");
    fh_px_vs_r->GetYaxis()->SetTitle("px / MeV/c");

    fh_py_vs_r = new TH2F("py_vs_targetR", "py vs Rtarget ", 250, 0, 2.5, 500, -500., 500.);
    fh_py_vs_r->GetXaxis()->SetTitle("Radius / cm");
    fh_py_vs_r->GetYaxis()->SetTitle("py / MeV/c");

    fh_pz_vs_r = new TH2F("pz_vs_targetR", "pz vs Rtarget ", 250, 0, 2.5, 1500, 16000., 19000.);
    fh_pz_vs_r->GetXaxis()->SetTitle("Radius / cm");
    fh_pz_vs_r->GetYaxis()->SetTitle("pz / MeV/c");

    fh_target_xy_nc = new TH2F("target_xy_nc", "target xy all chi2", 300, -1.5, 1.5, 300, -1.5, 1.5);
    fh_target_xy_nc->GetXaxis()->SetTitle("x / cm");
    fh_target_xy_nc->GetYaxis()->SetTitle("y / cm");

    fh_target_xy_MC = new TH2F("target_xy_MC", "target xy MC", 300, -1.5, 1.5, 300, -1.5, 1.5);
    fh_target_xy_MC->GetXaxis()->SetTitle("x / cm");
    fh_target_xy_MC->GetYaxis()->SetTitle("y / cm");

    fh_target_xx = new TH1F("target_xHe_xC", "diff target xC vs xHe ", 300, -1.5, 1.5);
    fh_target_xx->GetXaxis()->SetTitle("xHe-xC / cm");
    fh_target_xx->GetYaxis()->SetTitle("counts");

    fh_target_yy = new TH1F("target_yHe_yC", "target yC vs yHe ", 300, -1.5, 1.5);
    fh_target_yy->GetXaxis()->SetTitle("yHe-yC / cm");
    fh_target_yy->GetYaxis()->SetTitle("counts");

    fh_chi2 = new TH1F("chi2", "chi2 ", 5000, 0., 500);
    fh_chi2->GetXaxis()->SetTitle("Chi2");
    fh_chi2->GetYaxis()->SetTitle("counts");

    fh_dx = new TH1F("tracker_dx", "tracker dx ", 2000, -5., 5);
    fh_dx->GetXaxis()->SetTitle("dx / cm");
    fh_dx->GetYaxis()->SetTitle("counts");

    fh_dy = new TH1F("tracker_dy", "tracker dy ", 2000, -5., 5.);
    fh_dy->GetXaxis()->SetTitle("dy / cm");
    fh_dy->GetYaxis()->SetTitle("counts");

    fh_dz = new TH1F("tracker_dz", "tracker dz ", 2000, -5., 5.);
    fh_dz->GetXaxis()->SetTitle("dz / cm");
    fh_dz->GetYaxis()->SetTitle("counts");

    fh_dpxO = new TH1F("tracker_dpxO", "tracker dpx O", 1000, -50, 50);
    fh_dpxO->GetXaxis()->SetTitle("dpx / percent");
    fh_dpxO->GetYaxis()->SetTitle("counts");

    fh_dpyO = new TH1F("tracker_dpyO", "tracker dpy O", 1000, -50, 50);
    fh_dpyO->GetXaxis()->SetTitle("dpy / percent");
    fh_dpyO->GetYaxis()->SetTitle("counts");

    fh_dpzO = new TH1F("tracker_dpzO", "tracker dpz O", 1000, -50, 50);
    fh_dpzO->GetXaxis()->SetTitle("dpz / percent");
    fh_dpzO->GetYaxis()->SetTitle("counts");

    fh_dpO = new TH1F("tracker_dpO", "tracker dp O", 1000, -50, 50);
    fh_dpO->GetXaxis()->SetTitle("dp / percent");
    fh_dpO->GetYaxis()->SetTitle("counts");

    fh_dpxC = new TH1F("tracker_dpxC", "tracker dpx C", 2000, -100, 100);
    fh_dpxC->GetXaxis()->SetTitle("dpx / percent");
    fh_dpxC->GetYaxis()->SetTitle("counts");

    fh_dpyC = new TH1F("tracker_dpyC", "tracker dpy C", 2000, -100, 100);
    fh_dpyC->GetXaxis()->SetTitle("dpy / percent");
    fh_dpyC->GetYaxis()->SetTitle("counts");

    fh_dpzC = new TH1F("tracker_dpzC", "tracker dpzC", 2000, -100, 100);
    fh_dpzC->GetXaxis()->SetTitle("dpz / percent");
    fh_dpzC->GetYaxis()->SetTitle("counts");

    fh_dpC = new TH1F("tracker_dpC", "tracker dpC", 400, -20, 20);
    fh_dpC->GetXaxis()->SetTitle("dp / percent");
    fh_dpC->GetYaxis()->SetTitle("counts");

    fh_dp = new TH1F("tracker_dp", "tracker dp", 400, -20, 20);
    fh_dp->GetXaxis()->SetTitle("dp / percent");
    fh_dp->GetYaxis()->SetTitle("counts");

    fh_dpxHe = new TH1F("tracker_dpxHe", "tracker Hedpx ", 2000, -100, 100);
    fh_dpxHe->GetXaxis()->SetTitle("dpx / percent");
    fh_dpxHe->GetYaxis()->SetTitle("counts");

    fh_dpyHe = new TH1F("tracker_dpyHe", "tracker Hedpy ", 2000, -100, 100);
    fh_dpyHe->GetXaxis()->SetTitle("dpy / percent");
    fh_dpyHe->GetYaxis()->SetTitle("counts");

    fh_dpzHe = new TH1F("tracker_dpzHe", "tracker dpz He", 400, -20, 20);
    fh_dpzHe->GetXaxis()->SetTitle("dpz / percent");
    fh_dpzHe->GetYaxis()->SetTitle("counts");

    fh_dpHe = new TH1F("tracker_dpHe", "tracker dp He", 400, -20, 20);
    fh_dpHe->GetXaxis()->SetTitle("dp / percent");
    fh_dpHe->GetYaxis()->SetTitle("counts");

    fh_thetax_px = new TH2F("thethax_px", "theta_x vs px ", 500, -500, 500, 200, -100., 100.);
    fh_thetax_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_thetax_px->GetYaxis()->SetTitle("theta_x ");

    fh_thetay_py = new TH2F("thethay_py", "theta_y vs py ", 500, -500, 500, 200, -100., 100.);
    fh_thetay_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_thetay_py->GetYaxis()->SetTitle("theta_y ");

    fh_thetax_py = new TH2F("thethax_py", "theta_x vs py ", 200, -100, 100, 200, -100., 100.);
    fh_thetax_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_thetax_py->GetYaxis()->SetTitle("theta_x ");

    fh_thetay_px = new TH2F("thethay_px", "theta_y vs px ", 200, -100, 100, 200, -100., 100.);
    fh_thetay_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_thetay_px->GetYaxis()->SetTitle("theta_y ");

    fh_dpy_dpx = new TH2F("tracker_dpy_dpx", "tracker dpy vs dpx ", 200, -100, 100, 200, -100., 100.);
    fh_dpy_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_dpy_dpx->GetYaxis()->SetTitle("dpy / percent");

    fh_thetax_dpx_C = new TH2F("tracker_thethax_dpx_C", "tracker theta_x vs dpx C", 2000, -200, 200, 200, -50., 50.);
    fh_thetax_dpx_C->GetXaxis()->SetTitle("dpx / MeV/c");
    fh_thetax_dpx_C->GetYaxis()->SetTitle("theta_x / mrad ");

    fh_thetay_dpy_C = new TH2F("tracker_thethay_dpy_C", "tracker theta_y vs dpy_C ", 2000, -200, 200, 200, -50., 50.);
    fh_thetay_dpy_C->GetXaxis()->SetTitle("dpy / MeV/c");
    fh_thetay_dpy_C->GetYaxis()->SetTitle("theta_y / mrad ");

    fh_thetax_dpx_He = new TH2F("tracker_thethax_dpx_He", "tracker theta_x vs dpx He", 2000, -200, 200, 200, -50., 50.);
    fh_thetax_dpx_He->GetXaxis()->SetTitle("dpx / MeV/c");
    fh_thetax_dpx_He->GetYaxis()->SetTitle("theta_x / mrad ");

    fh_thetay_dpy_He =
        new TH2F("tracker_thethay_dpy_He", "tracker theta_y vs dpy_He ", 2000, -200, 200, 200, -50., 50.);
    fh_thetay_dpy_He->GetXaxis()->SetTitle("dpy / MeV/c");
    fh_thetay_dpy_He->GetYaxis()->SetTitle("theta_y / mrad ");

    fh_thetax_dpy = new TH2F("tracker_thethax_dpy", "tracker theta_x vs dpy ", 2000, -200, 200, 200, -50., 50.);
    fh_thetax_dpy->GetXaxis()->SetTitle("dpy / MeV/c");
    fh_thetax_dpy->GetYaxis()->SetTitle("theta_x / mrad ");

    fh_thetay_dpx = new TH2F("tracker_thethay_dpx", "tracker theta_y vs dpx ", 2000, -200, 200, 200, -50., 50.);
    fh_thetay_dpx->GetXaxis()->SetTitle("dpx / MeV/c");
    fh_thetay_dpx->GetYaxis()->SetTitle("theta_y / mrad ");

    fh_x_dpx = new TH2F("tracker_px_dpx", "tracker px vs dpx ", 1000, -50, 50, 400, -500., 500.);
    fh_x_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_x_dpx->GetYaxis()->SetTitle("px / MeV/c");

    fh_y_dpy = new TH2F("tracker_py_dpy", "tracker py vs dpy ", 1000, -50, 50, 400, -500., 500.);
    fh_y_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_y_dpy->GetYaxis()->SetTitle("py / MeV/c");

    fh_px_px = new TH2F("px_px", "px He vs. px C", 400, -500, 500, 400, -500., 500.);
    fh_px_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_px->GetYaxis()->SetTitle("px / MeV/c");

    fh_py_py = new TH2F("py_py", "py He vs. py C", 400, -500, 500, 400, -500., 500.);
    fh_py_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_py->GetYaxis()->SetTitle("py / MeV/c");

    fh_pz_pz = new TH2F("pz_pz", "pz He vs. pz C", 400, 0, 10000, 400, 10000., 20000.);
    fh_pz_pz->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_pz->GetYaxis()->SetTitle("pz / MeV/c");

    fh_p_p = new TH2F("p_p", "p He vs. p C", 400, 0, 10000, 400, 10000., 20000.);
    fh_p_p->GetXaxis()->SetTitle("p / MeV/c");
    fh_p_p->GetYaxis()->SetTitle("p / MeV/c");

    fh_px_px_nc = new TH2F("px_px_nc", "px He vs. px C all chi2", 400, -500, 500, 400, -500., 500.);
    fh_px_px_nc->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_px_nc->GetYaxis()->SetTitle("px / MeV/c");

    fh_py_py_nc = new TH2F("py_py_nc", "py He vs. py C all chi2", 400, -500, 500, 400, -500., 500.);
    fh_py_py_nc->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_py_nc->GetYaxis()->SetTitle("py / MeV/c");

    fh_pz_pz_nc = new TH2F("pz_pz_nc", "pz He vs. pz C all chi2", 400, 2000, 6000, 500, 10000., 20000.);
    fh_pz_pz_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_pz_nc->GetYaxis()->SetTitle("pz / MeV/c");

    fh_px_px_mc = new TH2F("px_px_MC", "px He vs. px C MC", 400, -500, 500, 400, -500., 500.);
    fh_px_px_mc->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_px_mc->GetYaxis()->SetTitle("px / MeV/c");

    fh_py_py_mc = new TH2F("py_py_MC", "py He vs. py C MC", 400, -500, 500, 400, -500., 500.);
    fh_py_py_mc->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_py_mc->GetYaxis()->SetTitle("py / MeV/c");

    fh_pz_pz_mc = new TH2F("pz_pz_MC", "pz He vs. pz C MC", 400, 2000, 6000, 1000, 10000., 20000.);
    fh_pz_pz_mc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_pz_mc->GetYaxis()->SetTitle("pz / MeV/c");

    // Analysis
    fh_chiy_vs_chix_He_nc = new TH2F("chiC_vs_chiHe", "chi2 C vs chi2 He", 1000, 0., 50000., 1000, 0., 50000.);
    fh_chiy_vs_chix_He_nc->GetXaxis()->SetTitle("chi2 He");
    fh_chiy_vs_chix_He_nc->GetYaxis()->SetTitle("chi2 C");

    fh_chiy_vs_chix_He =
        new TH2F("chiC_vs_chiHe_after_cut", "chi2 C vs chi2 He after cut", 1000, 0., 50000., 1000, 0., 50000.);
    fh_chiy_vs_chix_He->GetXaxis()->SetTitle("chi2 He");
    fh_chiy_vs_chix_He->GetYaxis()->SetTitle("chi2 C");

    fh_chiy_vs_chix_C = new TH1F("chi2summ", "chi2sum", 5000, 0., 1000.);
    fh_chiy_vs_chix_C->GetXaxis()->SetTitle("chi2");
    fh_chiy_vs_chix_C->GetYaxis()->SetTitle("counts");

    fh_theta26 = new TH1F("theta26", "theta between alpha and carbon ", 125, 0., 5);
    fh_theta26->GetXaxis()->SetTitle("angle / degree");
    fh_theta26->GetYaxis()->SetTitle("counts");

    fh_theta26_nc = new TH1F("theta26_nc", "theta between alpha and carbon all chi2 ", 125, 0., 5);
    fh_theta26_nc->GetXaxis()->SetTitle("angle / degree");
    fh_theta26_nc->GetYaxis()->SetTitle("counts");

    fh_theta26_vs_chi = new TH2F("theta_vs_chi", "theta vs chi", 1000, 0., 500., 125, 0., 5);
    fh_theta26_vs_chi->GetXaxis()->SetTitle("chi ");
    fh_theta26_vs_chi->GetYaxis()->SetTitle("angle / degree");

    fh_theta_16O = new TH1F("theta_16O", "theta of excited 16O ", 125, 0., 5);
    fh_theta_16O->GetXaxis()->SetTitle("angle / degree");
    fh_theta_16O->GetYaxis()->SetTitle("counts");

    fh_phi_16O = new TH1F("phi_16O", "phi of excited 16O ", 360, 0., 360);
    fh_phi_16O->GetXaxis()->SetTitle("angle / degree");
    fh_phi_16O->GetYaxis()->SetTitle("counts");

    fh_theta26_cm = new TH1F("theta26_cm", "theta between alpha and carbon in cm-system ", 360, 0, 360.);
    fh_theta26_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta26_cm->GetYaxis()->SetTitle("counts");

    fh_phi26_cm = new TH1F("phi26_cm", "phi between alpha and carbon in cm-system", 360, 0., 360.);
    fh_phi26_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi26_cm->GetYaxis()->SetTitle("counts");

    fh_theta_4He_cm = new TH1F("theta_4He_cm", "Theta 4He in cm-system", 180, 0., 180.);
    fh_theta_4He_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_4He_cm->GetYaxis()->SetTitle("counts");

    fh_phi_4He_cm = new TH1F("phi_4He_cm", "phi 4He in cm-system", 360, 0., 360.);
    fh_phi_4He_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_4He_cm->GetYaxis()->SetTitle("counts");

    fh_theta_bc_cm = new TH1F("theta_bc_cm", "Theta bc in cm-system", 180, 0., 180.);
    fh_theta_bc_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_bc_cm->GetYaxis()->SetTitle("counts");

    fh_theta_bc_cm_nc = new TH1F("theta_bc_cm_nc", "Theta bc in cm-system all chi2", 180, 0., 180.);
    fh_theta_bc_cm_nc->GetXaxis()->SetTitle("angle / degree");
    fh_theta_bc_cm_nc->GetYaxis()->SetTitle("counts");

    fh_phi_bc_cm_nc = new TH1F("phi_bc_cm_nc", "phi bc in cm-system all chi2", 360, 0., 360.);
    fh_phi_bc_cm_nc->GetXaxis()->SetTitle("angle / degree");
    fh_phi_bc_cm_nc->GetYaxis()->SetTitle("counts");

    fh_phi_bc_cm = new TH1F("phi_bc_cm", "phi bc in cm-system", 360, 0., 360.);
    fh_phi_bc_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_bc_cm->GetYaxis()->SetTitle("counts");

    fh_phiMC_bc_cm = new TH1F("phiMC_bc_cm", "phiMC bc in cm-system", 360, 0., 360.);
    fh_phiMC_bc_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phiMC_bc_cm->GetYaxis()->SetTitle("counts");

    fh_phi_bc_cm_polar = new TH2F("phi_bc_cm_polar", "phi_bc_cm_polar", 720, 0., 360., 100, 0., 10.);
    fh_phi_bc_cm_polar->GetXaxis()->SetTitle("angle / degree");
    fh_phi_bc_cm_polar->GetYaxis()->SetTitle("counts");

    fh_phiMC_bc_cm_polar = new TH2F("phiMC_bc_cm_polar", "phiMC_bc_cm_polar", 720, 0., 360., 100, 0., 10.);
    fh_phiMC_bc_cm_polar->GetXaxis()->SetTitle("angle / degree");
    fh_phiMC_bc_cm_polar->GetYaxis()->SetTitle("counts");

    fh_theta_12C_cm = new TH1F("theta_12C_cm", "Theta 12C in cm-system", 180, 0., 180.);
    fh_theta_12C_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_12C_cm->GetYaxis()->SetTitle("counts");

    fh_phi_12C_cm = new TH1F("phi_12C_cm", "phi 12C in cm-system", 360, 0., 360.);
    fh_phi_12C_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_12C_cm->GetYaxis()->SetTitle("counts");

    fh_theta26_simu = new TH1F("theta26_simu", "theta 26 simulation", 500, 0., 5);
    fh_theta26_simu->GetXaxis()->SetTitle("angle / degree");
    fh_theta26_simu->GetYaxis()->SetTitle("counts");

    fh_minv_simu = new TH1F("minv_MC", "minv_MC ", 500, 14900, 14920);
    fh_minv_simu->GetXaxis()->SetTitle("minv / MeV");
    fh_minv_simu->GetYaxis()->SetTitle("counts");

    fh_minv = new TH1F("minv", "minv ", 2000, 10, 30);
    fh_minv->GetXaxis()->SetTitle("minv / MeV");
    fh_minv->GetYaxis()->SetTitle("counts");

    fh_Erel = new TH1F("Erel", "Erel ", 500, 0., 50);
    fh_Erel->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel->GetYaxis()->SetTitle("counts");

    fh_Erel_vs_chi = new TH2F("Erel_vs_chi", "Erel vs chi", 100, 0, 50, 200, 0., 20);
    fh_Erel_vs_chi->GetYaxis()->SetTitle("Erel / MeV");
    fh_Erel_vs_chi->GetXaxis()->SetTitle("chi2");

    fh_phi26_vs_chi = new TH2F("phi26_vs_chi", "phi26 vs chi", 1000, 0, 500, 360, 0., 360);
    fh_phi26_vs_chi->GetYaxis()->SetTitle("phi / deg");
    fh_phi26_vs_chi->GetXaxis()->SetTitle("chi2");

    fh_psum_vs_theta26 = new TH2F("psum_vs_theta26", "psum_vs_theta26", 125, 0., 5, 600, 16000., 19000);
    fh_psum_vs_theta26->GetYaxis()->SetTitle("theta26 / deg");
    fh_psum_vs_theta26->GetXaxis()->SetTitle("psum / MeV/c");

    fh_Erel_nc = new TH1F("Erel_nc", "Erel all chi2", 1000, 0., 100);
    fh_Erel_nc->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel_nc->GetYaxis()->SetTitle("counts");

    fh_ErelL = new TH1F("ErelL", "Erel links", 500, 0., 50);
    fh_ErelL->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelL->GetYaxis()->SetTitle("counts");

    fh_ErelR = new TH1F("ErelR", "Erel rechts ", 500, 0., 50);
    fh_ErelR->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelR->GetYaxis()->SetTitle("counts");

    fh_ErelB = new TH1F("ErelB", "Erel both sides ", 500, 0., 50);
    fh_ErelB->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelB->GetYaxis()->SetTitle("counts");

    fh_ErelB_nc = new TH1F("ErelB_nc", "Erel both sides all chi2", 1000, 0., 100);
    fh_ErelB_nc->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelB_nc->GetYaxis()->SetTitle("counts");

    fh_Erel_simu = new TH1F("Erel_simu", "Erel simulation", 500, 0., 50);
    fh_Erel_simu->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel_simu->GetYaxis()->SetTitle("counts");

    fh_Erel_dErel = new TH2F("dpsum_vs_dErel", "dpsum vs dErel", 500, -10, 10, 500, -500., 500);
    fh_Erel_dErel->GetXaxis()->SetTitle("dErel / MeV");
    fh_Erel_dErel->GetYaxis()->SetTitle("dpsum / MeV/c");

    fh_phibcm_vs_px = new TH2F("phibcm_vs_px", "phibcm vs. px", 400, -500., 500, 360, 0, 360.);
    fh_phibcm_vs_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_phibcm_vs_px->GetYaxis()->SetTitle("phibcm / deg");

    fh_phibcm_vs_py = new TH2F("phibcm_vs_py", "phibcm vs. py", 400, -500., 500, 360, 0, 360.);
    fh_phibcm_vs_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_phibcm_vs_py->GetYaxis()->SetTitle(" phibcm / deg");

    fh_phibcm_vs_pz = new TH2F("phibcm_vs_pz", "phibcm vs. pz", 600, 2000., 18000, 360, 0, 360);
    fh_phibcm_vs_pz->GetXaxis()->SetTitle("pz / MeV/c");
    fh_phibcm_vs_pz->GetYaxis()->SetTitle(" phibcm / deg");

    fh_phibcm_vs_psum = new TH2F("phibcm_vs_psum", "phibcm vs. pz", 600, 16000., 19000, 360, 0, 360);
    fh_phibcm_vs_psum->GetXaxis()->SetTitle("psum / MeV/c");
    fh_phibcm_vs_psum->GetYaxis()->SetTitle(" phibcm / deg");

    fh_phibcm_vs_px_MC = new TH2F("phibcm_vs_px_MC", "phibcm vs. px MC", 400, -500., 500, 360, 0, 360.);
    fh_phibcm_vs_px_MC->GetXaxis()->SetTitle("px / MeV/c");
    fh_phibcm_vs_px_MC->GetYaxis()->SetTitle("phibcm / deg");

    fh_phibcm_vs_py_MC = new TH2F("phibcm_vs_py_MC", "phibcm vs. py MC", 400, -500., 500, 360, 0, 360.);
    fh_phibcm_vs_py_MC->GetXaxis()->SetTitle("py / MeV/c");
    fh_phibcm_vs_py_MC->GetYaxis()->SetTitle(" phibcm / deg");

    fh_phibcm_vs_pz_MC = new TH2F("phibcm_vs_pz_MC", "phibcm vs. pz MC", 600, 2000., 18000, 360, 0, 360);
    fh_phibcm_vs_pz_MC->GetXaxis()->SetTitle("pz / MeV/c");
    fh_phibcm_vs_pz_MC->GetYaxis()->SetTitle(" phibcm / deg");

    fh_phibcm_vs_psum_MC = new TH2F("phibcm_vs_psum_MC", "phibcm vs. pz MC", 600, 16000., 19000, 360, 0, 360);
    fh_phibcm_vs_psum_MC->GetXaxis()->SetTitle("psum / MeV/c");
    fh_phibcm_vs_psum_MC->GetYaxis()->SetTitle(" phibcm / deg");

    fh_Erel_vs_psum = new TH2F("Erel_vs_psum", "Erel vs. psum", 500, 16000., 18000., 100, 0, 10.);
    fh_Erel_vs_psum->GetXaxis()->SetTitle("psum / MeV/c");
    fh_Erel_vs_psum->GetYaxis()->SetTitle(" Erel / MeV");

    fh_mass_nc = new TH2F("mass_nc", "mass He vs mass C all chi2", 2000, 0, 20, 500, 0, 5);
    fh_mass_nc->GetXaxis()->SetTitle("mass C");
    fh_mass_nc->GetYaxis()->SetTitle(" mass H");

    fh_mass = new TH2F("mass", "mass He vs mass C selected chi2", 2000, 0, 20, 500, 0, 5);
    fh_mass->GetXaxis()->SetTitle("mass C");
    fh_mass->GetYaxis()->SetTitle(" mass H");

    fh_energy_nc = new TH2F("energy_nc", "energy He vs energy C all chi2", 1000, 10, 20, 1000, 0, 10);
    fh_energy_nc->GetXaxis()->SetTitle("energy C / GeV");
    fh_energy_nc->GetYaxis()->SetTitle(" energy H / GeV");

    fh_energy = new TH2F("energy", "energy He vs energy C selected chi2", 1000, 10, 20, 1000, 0, 10);
    fh_energy->GetXaxis()->SetTitle("energy C / GeV");
    fh_energy->GetYaxis()->SetTitle(" energy H / GeV");

    fh_psum = new TH1F("sum_p", "ptot", 500, 16000., 18000.);
    fh_psum->GetXaxis()->SetTitle("p / MeV/c");
    fh_psum->GetYaxis()->SetTitle("count");

    fh_pzsum = new TH1F("sum_pz", "pztot", 500, 16000., 18000.);
    fh_pzsum->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pzsum->GetYaxis()->SetTitle("count");

    fh_psum_nc = new TH1F("sum_p_nc", "ptot all chi2", 500, 16000., 18000.);
    fh_psum_nc->GetXaxis()->SetTitle("p / MeV/c");
    fh_psum_nc->GetYaxis()->SetTitle("count");

    fh_psum_vs_chi = new TH2F("sum_p_vs_chi", "ptot vs chi2", 1000, 0, 500, 400, 16800., 17800.);
    fh_psum_vs_chi->GetYaxis()->SetTitle("p / MeV/c");
    fh_psum_vs_chi->GetXaxis()->SetTitle("chi2");

    fh_psum_vs_event = new TH2F("psum_vs_event", "psum vs. event number", 10000, 0, 1e6, 80, 17000., 17800.);
    fh_psum_vs_event->GetXaxis()->SetTitle("Event number");
    fh_psum_vs_event->GetYaxis()->SetTitle("p / MeV/c");

    fh_psum_vs_event_nc =
        new TH2F("psum_vs_event_nc", "psum vs. event number all chi2", 10000, 0, 1e6, 80, 17000., 17800.);
    fh_psum_vs_event_nc->GetXaxis()->SetTitle("Event number");
    fh_psum_vs_event_nc->GetYaxis()->SetTitle("p / MeV/c");

    fh_psum_MC = new TH1F("sum_p_MC", "sum_p_MC", 2000, 16000., 18000.);
    fh_psum_MC->GetXaxis()->SetTitle("p / MeV/c");
    fh_psum_MC->GetYaxis()->SetTitle("count");

    fh_pzsum_MC = new TH1F("sum_pz_MC", "sum_pz_MC", 2000, 16000., 18000.);
    fh_pzsum_MC->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pzsum_MC->GetYaxis()->SetTitle("count");

    fh_dErel = new TH1F("dErel", "delta Erel ", 400, -100, 100.);
    fh_dErel->GetYaxis()->SetTitle("counts");
    fh_dErel->GetXaxis()->SetTitle("Delta Erel / %");

    fh_dtheta = new TH1F("dtheta26_cm", "delta theta ", 1000, -5., 5.);
    fh_dtheta->GetYaxis()->SetTitle("counts");
    fh_dtheta->GetXaxis()->SetTitle("Delta theta / deg");

    fh_Erel_vs_theta26 = new TH2F("Erel_vs_theta", "Erel both sides vs. theta", 125, 0., 5, 200, 0, 20.);
    fh_Erel_vs_theta26->GetXaxis()->SetTitle("angle / deg");
    fh_Erel_vs_theta26->GetYaxis()->SetTitle("Erel / MeV");

    fh_Erel_vs_theta26_nc = new TH2F("Erel_vs_theta_nc", "Erel both sides vs. theta all ci2", 250, 0., 5, 200, 0, 20.);
    fh_Erel_vs_theta26_nc->GetXaxis()->SetTitle("angle / deg");
    fh_Erel_vs_theta26_nc->GetYaxis()->SetTitle("Erel / MeV");

    fh_Erel_vs_thetaMC = new TH2F("Erel_vs_thetaMC", "Erel both sides vs. theta MC", 250, 0., 5, 200, 0, 20.);
    fh_Erel_vs_thetaMC->GetXaxis()->SetTitle("angle / deg");
    fh_Erel_vs_thetaMC->GetYaxis()->SetTitle("Erel / MeV");

    fh_Erel_vs_phibc = new TH2F("Erel_vs_phibccm", "Erel both sides vs. phi_bc_cm", 200, 0, 20., 400, 0., 400.);
    fh_Erel_vs_phibc->GetYaxis()->SetTitle("phi angle bc / deg");
    fh_Erel_vs_phibc->GetXaxis()->SetTitle("Erel / MeV");

    fh_Erel_vs_phibcMC =
        new TH2F("Erel_vs_phibccmMC", "Erel both sides vs. phi_bc_cm MCtrack", 200, 0, 20., 400, 0., 400.);
    fh_Erel_vs_phibcMC->GetYaxis()->SetTitle("phi angle bc / deg");
    fh_Erel_vs_phibcMC->GetXaxis()->SetTitle("Erel / MeV");

    fh_pHe_vs_theta26 = new TH2F("pHe_vs_theta", "theta vs pHe", 400, 2000, 6000, 200, 0, 100.);
    fh_pHe_vs_theta26->GetYaxis()->SetTitle("angle / mrad");
    fh_pHe_vs_theta26->GetXaxis()->SetTitle("p / MeV/c");

    fh_psum_vs_theta26_nc = new TH2F("psum_vs_theta_nc", "theta vs psum nc", 125, 0., 5, 600, 16000., 19000);
    fh_psum_vs_theta26_nc->GetYaxis()->SetTitle("angle / deg");
    fh_psum_vs_theta26_nc->GetXaxis()->SetTitle("p / MeV/c");

    if (fMappedItemsCalifa)
    {
        fh_califa_hitenergy =
            new TH2F("fh_califa_hit_energy", "Califa hitE vs crystal id", 100, 0, 100, 2000, 0., 10000.);
        fh_califa_hitenergy->GetYaxis()->SetTitle("Energy / keV");
        fh_califa_hitenergy->GetXaxis()->SetTitle("Cluster #");

        fh_califa_calenergy =
            new TH2F("fh_califa_cal_energy", "Califa calE vs crystal id", 5000, 0, 5000, 2000, 0., 10000.);
        fh_califa_calenergy->GetYaxis()->SetTitle("Energy / keV");
        fh_califa_calenergy->GetXaxis()->SetTitle("Crystal #");

        fh_califa_energy = new TH2F("fh_califa_energy", "Califa E vs crystal id", 5000, 0, 5000, 2000, 0., 10000.);
        fh_califa_energy->GetYaxis()->SetTitle("Energy / keV");
        fh_califa_energy->GetXaxis()->SetTitle("Crystal #");

        fh_califa_energy_nc = new TH2F(
            "fh_califa_energy_nc", "Califa E vs crystal id no pair selection", 5000, 0, 5000, 1000, 0., 10000.);
        fh_califa_energy_nc->GetYaxis()->SetTitle("Energy / keV");
        fh_califa_energy_nc->GetXaxis()->SetTitle("Crystal #");
    }

    TCanvas* check = new TCanvas("CheckingGlobal", "CheckingGlobal", 10, 10, 900, 900);
    if (fPairs)
    {
        check->Divide(4, 4);
        check->cd(1);
        fh_chi2->Draw();
        check->cd(2);
        gPad->SetLogz();
        fh_psum_vs_chi->Draw("colz");
        check->cd(3);
        gPad->SetLogz();
        fh_theta26_vs_chi->Draw("colz");
        check->cd(4);
        gPad->SetLogz();
        fh_Erel_vs_chi->Draw("colz");
        check->cd(5);
        gPad->SetLogz();
        fh_target_xy_nc->Draw("colz");
        check->cd(6);
        gPad->SetLogz();
        fh_px_px_nc->Draw("colz");
        check->cd(7);
        gPad->SetLogz();
        fh_py_py_nc->Draw("colz");
        check->cd(8);
        gPad->SetLogz();
        fh_pz_pz_nc->Draw("colz");
        check->cd(9);
        fh_p_He_nc->Draw();
        check->cd(10);
        fh_p_C_nc->Draw();
        check->cd(11);
        fh_psum_nc->Draw();
        check->cd(12);
        fh_Erel_nc->Draw();
        check->cd(13);
        gPad->SetLogz();
        fh_Erel_vs_theta26_nc->Draw("colz");
        check->cd(14);
        fh_theta26_nc->Draw();
        check->cd(15);
        fh_theta_bc_cm_nc->Draw();
        check->cd(16);
        fh_phi_bc_cm_nc->Draw();

        TCanvas* check1 = new TCanvas("SelectedChi2", "SelectedChi2", 10, 10, 900, 900);
        check1->Divide(4, 3);
        check1->cd(1);
        gPad->SetLogz();
        fh_target_xy->Draw("colz");
        check1->cd(2);
        gPad->SetLogz();
        fh_px_px->Draw("colz");
        check1->cd(3);
        gPad->SetLogz();
        fh_py_py->Draw("colz");
        check1->cd(4);
        gPad->SetLogz();
        fh_pz_pz->Draw("colz");
        check1->cd(5);
        fh_p_He->Draw();
        check1->cd(6);
        fh_p_C->Draw();
        check1->cd(7);
        fh_psum->Draw();
        check1->cd(8);
        fh_Erel->Draw();
        check1->cd(9);
        gPad->SetLogz();
        fh_Erel_vs_theta26->Draw("colz");
        check1->cd(10);
        fh_theta26->Draw();
        check1->cd(11);
        fh_theta_bc_cm->Draw();
        check1->cd(12);
        fh_phi_bc_cm->Draw();
    }
    else
    {

        check->Divide(4, 4);
        check->cd(1);
        fh_chi2->Draw();
        check->cd(2);
        gPad->SetLogz();
        fh_psum_vs_chi->Draw("colz");
        check->cd(3);
        gPad->SetLogz();
        fh_target_xy_nc->Draw("colz");
        check->cd(4);
        gPad->SetLogz();
        fh_psum_vs_event_nc->Draw("colz");
        check->cd(5);
        fh_p_O_nc->Draw();
        check->cd(6);
        fh_px_O_nc->Draw();
        check->cd(7);
        fh_py_O_nc->Draw();
        check->cd(8);
        fh_pz_O_nc->Draw();
        check->cd(11);
        gPad->SetLogz();
        fh_target_xy->Draw("colz");
        check->cd(12);
        gPad->SetLogz();
        fh_psum_vs_event->Draw("colz");
        check->cd(13);
        fh_p_O->Draw();
        check->cd(14);
        fh_px_O->Draw();
        check->cd(15);
        fh_py_O->Draw();
        check->cd(16);
        fh_pz_O->Draw();
    }

    if (fSimu)
    {
        TCanvas* checkMC = new TCanvas("CheckingGlobalMC", "CheckingGlobalMC", 10, 10, 1200, 900);

        if (fPairs)
        {
            checkMC->Divide(4, 3);
            checkMC->cd(1);
            fh_dx->Draw();
            checkMC->cd(2);
            fh_dy->Draw();
            checkMC->cd(3);
            fh_dp->Draw();
            checkMC->cd(4);
            gPad->SetLogz();
            fh_dpy_dpx->Draw("colz");
            checkMC->cd(5);
            fh_x_dpx->Draw("colz");
            checkMC->cd(6);
            fh_y_dpy->Draw("colz");
            checkMC->cd(7);
            fh_dErel->Draw();
            checkMC->cd(8);
            fh_target_xy->Draw("colz");
            checkMC->cd(9);
            fh_target_xy_MC->Draw("colz");
            checkMC->cd(10);
            gPad->SetLogz();
            fh_thetay_dpy_C->Draw("colz");
            checkMC->cd(11);
            gPad->SetLogz();
            fh_thetax_dpx_He->Draw("colz");
            checkMC->cd(12);
            gPad->SetLogz();
            fh_thetay_dpy_He->Draw("colz");
        }
        else
        {
            checkMC->Divide(3, 3);
            checkMC->cd(1);
            fh_dx->Draw();
            checkMC->cd(2);
            fh_dy->Draw();
            checkMC->cd(3);
            fh_dpxO->Draw();
            checkMC->cd(4);
            fh_dpyO->Draw();
            checkMC->cd(5);
            fh_dpzO->Draw();
            checkMC->cd(6);
            fh_dpO->Draw();
            checkMC->cd(7);
            gPad->SetLogz();
            fh_dpy_dpx->Draw("colz");
        }
    }

    return kSUCCESS;
}

void R3BGlobalAnalysisS494::Exec(Option_t* option)
{
    if (fNEvents / 100000. == (int)fNEvents / 100000)
        std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
                  << " %) " << std::flush;

    // cout << "New event " << fNEvents << endl;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    pHex = -100.;
    pHey = -100.;
    pHez = -100.;
    pCx = -100.;
    pCy = -100.;
    pCz = -100.;
    XHe = -100;
    YHe = -100;
    ZHe = -100;
    XC = -100;
    YC = -100;
    ZC = -100;
    mtrackHe = -100;
    mtrackC = -100;
    pHex_mc = -100.;
    pHey_mc = -100.;
    pHez_mc = -100.;
    pCx_mc = -100.;
    pCy_mc = -100.;
    pCz_mc = -100.;
    XHe_mc = -100;
    YHe_mc = -100;
    ZHe_mc = -100;
    XC_mc = -100;
    YC_mc = -100;
    ZC_mc = -100;
    Double_t px, py, pz;
    Double_t theta_16, theta_26;
    Double_t costh26;
    Double_t chiHex = 1.e+36, chiHey = 1.e+36, chiCx = 1.e+36, chiCy = 1.e+36;
    Double_t chix = 1.e+36, chiy = 1.e+36;

    Bool_t is_alpha = false;
    Bool_t is_carbon = false;
    Bool_t is_oxygen = false;
    Bool_t is_tracked = false;

    Int_t nHitsTrack = 0;
    Int_t nHitsMCTrack = 0;

    if (fTrack)
    {

        // read in Track parameter
        nHitsTrack = fTrack->GetEntriesFast();

        //  if(nHitsTrack < 1) return;

        // cout << "Track hits: " << nHitsTrack << endl;
        for (Int_t l = 0; l < nHitsTrack; l++)
        {
            R3BTrack* aTrack = (R3BTrack*)fTrack->At(l);

            //  LOG(DEBUG) << "Charge " << aTrack->GetQ() << endl;

            if (aTrack->GetQ() == 2)
            {
                // 4He
                XHe = aTrack->GetX(); // cm
                YHe = aTrack->GetY(); // cm
                ZHe = aTrack->GetZ(); // cm

                pHex = aTrack->GetPx(); //  + ((std::rand() / (float)RAND_MAX) - 0.5)/100.* aTrack->GetPx(); // MeV/c
                pHey = aTrack->GetPy(); // + ((std::rand() / (float)RAND_MAX) - 0.5)/100.* aTrack->GetPx(); // MeV/c
                pHez = aTrack->GetPz(); // MeV/c
                mtrackHe = aTrack->GetAoZ() * aTrack->GetQ();

                chiHex = aTrack->GetChix();
                chiHey = aTrack->GetChiy();

                is_alpha = true;
                alpha.SetPxPyPzE(pHex, pHey, pHez, sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2) + pow(mHe, 2)));

                LOG(DEBUG) << "******************************************" << endl;
                LOG(DEBUG) << "Track In 4He"
                           << "x " << XHe << " y " << YHe << " z " << ZHe << endl;
                LOG(DEBUG) << "px " << pHex << " py " << pHey << " z " << pHez << endl;
                LOG(DEBUG) << "chiHex " << chiHex << ", " << chiHey << endl;
            }
            if (aTrack->GetQ() == 6)
            {
                // 12C
                XC = aTrack->GetX();
                YC = aTrack->GetY();
                ZC = aTrack->GetZ();
                mtrackC = aTrack->GetAoZ() * aTrack->GetQ();

                pCx = aTrack->GetPx(); // + ((std::rand() / (float)RAND_MAX) - 0.5)/100.* aTrack->GetPx();
                pCy = aTrack->GetPy(); // + ((std::rand() / (float)RAND_MAX) - 0.5)/100.* aTrack->GetPx();
                pCz = aTrack->GetPz();

                chiCx = aTrack->GetChix();
                chiCy = aTrack->GetChiy();

                is_carbon = true;
                carbon.SetPxPyPzE(pCx, pCy, pCz, sqrt(pow(pCx, 2) + pow(pCy, 2) + pow(pCz, 2) + pow(mC, 2)));

                LOG(DEBUG) << "******************************************" << endl;
                LOG(DEBUG) << "Track In 12C"
                           << "x " << XC << " y " << YC << " z " << ZC << endl;
                LOG(DEBUG) << "px " << pCx << " py " << pCy << " z " << pCz << endl;
                LOG(DEBUG) << "chiCx " << chiCx << ", " << chiCy << endl;
            }
            if (aTrack->GetQ() == 8)
            {
                // 160
                Xf = aTrack->GetX();
                Yf = aTrack->GetY();
                Zf = aTrack->GetZ();

                Pxf = aTrack->GetPx();
                Pyf = aTrack->GetPy();
                Pzf = aTrack->GetPz();
                is_oxygen = true;
                Pf_tot = sqrt((Pxf * Pxf) + (Pyf * Pyf) + (Pzf * Pzf));

                chix = aTrack->GetChix();
                chiy = aTrack->GetChiy();

                LOG(DEBUG) << "******************************************" << endl;
                LOG(DEBUG) << "Track In 16O"
                           << "x " << Xf << " y " << Yf << " z " << Zf << endl;
                LOG(DEBUG) << "px " << Pxf << " py " << Pyf << " z " << Pzf << endl;
            }
        }
    }
    if (fMCTrack && fSimu)
    {
        // read in Monte Carlo Track parameter

        nHitsMCTrack = fMCTrack->GetEntriesFast();
        //  cout << "MCTrack hits: " << nHitsMCTrack << endl;

        for (Int_t l = 0; l < nHitsMCTrack; l++)
        {
            //  cout << "Original MC Data is analyzed" << endl;
            R3BMCTrack* aTrack = (R3BMCTrack*)fMCTrack->At(l);

            Int_t PID = aTrack->GetPdgCode();
            Int_t mother = aTrack->GetMotherId();
            //  LOG(DEBUG) << "PID " << PID << endl;
            if (mother < 0)
            {
                if (PID == 1000020040)
                {
                    // 4He
                    XHe_mc = aTrack->GetStartX(); // cm
                    YHe_mc = aTrack->GetStartY(); // cm
                    ZHe_mc = aTrack->GetStartZ(); // cm
                    THe_mc = aTrack->GetStartT(); // cm

                    pHex_mc = aTrack->GetPx() * 1000.; // MeV/c
                    pHey_mc = aTrack->GetPy() * 1000.;
                    pHez_mc = aTrack->GetPz() * 1000.;

                    is_alpha = true;
                    alphaMC.SetPxPyPzE(pHex_mc,
                                       pHey_mc,
                                       pHez_mc,
                                       sqrt(pow(pHex_mc, 2) + pow(pHey_mc, 2) + pow(pHez_mc, 2) + pow(mHe, 2)));

                    //  LOG(DEBUG) << "******************************************" << endl;
                    //  LOG(DEBUG) << "MC Track In 4He"
                    //             << "x " << XHe_mc << " y " << YHe_mc << " z " << ZHe_mc << endl;
                    // LOG(DEBUG) << "px " << pHex_mc << " py " << pHey_mc << " pz " << pHez_mc << endl;
                }
                if (PID == 1000060120)
                {
                    // 12C
                    XC_mc = aTrack->GetStartX();
                    YC_mc = aTrack->GetStartY();
                    ZC_mc = aTrack->GetStartZ();
                    TC_mc = aTrack->GetStartT();

                    pCx_mc = aTrack->GetPx() * 1000.;
                    pCy_mc = aTrack->GetPy() * 1000.;
                    pCz_mc = aTrack->GetPz() * 1000.;

                    is_carbon = true;
                    carbonMC.SetPxPyPzE(
                        pCx_mc, pCy_mc, pCz_mc, sqrt(pow(pCx_mc, 2) + pow(pCy_mc, 2) + pow(pCz_mc, 2) + pow(mC, 2)));

                    //  LOG(DEBUG) << "******************************************" << endl;
                    //  LOG(DEBUG) << "MC Track In 12C"
                    //             << "x " << XC_mc << " y " << YC_mc << " z " << ZC_mc << endl;
                    //  LOG(DEBUG) << "px " << pCx_mc << " py " << pCy_mc << " pz " << pCz_mc << endl;
                }
                if (PID == 1000080160)
                {
                    // 16O
                    Xf_mc = aTrack->GetStartX();
                    Yf_mc = aTrack->GetStartY();
                    Zf_mc = aTrack->GetStartZ();

                    Pxf_mc = aTrack->GetPx() * 1000.;
                    Pyf_mc = aTrack->GetPy() * 1000.;
                    Pzf_mc = aTrack->GetPz() * 1000.;

                    is_oxygen = true;
                    Pf_tot_mc = sqrt((Pxf_mc * Pxf_mc) + (Pyf_mc * Pyf_mc) + (Pzf_mc * Pzf_mc));

                    LOG(DEBUG) << "******************************************" << endl;
                    LOG(DEBUG) << "MC Track In 16O"
                               << "x " << Xf_mc << " y " << Yf_mc << " z " << Zf_mc << endl;
                    LOG(DEBUG) << "px " << Pxf_mc << " py " << Pyf_mc << " z " << Pzf_mc << endl;
                }
            }
        }
    }

    if (fMappedItemsCalifa)
    {
        // CALIFA
        auto detCalifa = fMappedItemsCalifa;
        Int_t nHitsCalifa = detCalifa->GetEntriesFast();
        // cout<<"Califa hits: "<<nHitsCalifa<<endl;

        for (Int_t ihit = 0; ihit < nHitsCalifa; ihit++)
        {
            R3BCalifaMappedData* hitCalifa = (R3BCalifaMappedData*)detCalifa->At(ihit);
            if (!hitCalifa)
                continue;

            Int_t Crystal = hitCalifa->GetCrystalId();
            Int_t Energy = hitCalifa->GetEnergy();
            // cout << "Califa: " << Crystal << " Energy: " << Energy << endl;
            if (Energy > 0)
            {
                fh_califa_energy_nc->Fill(Crystal, Energy);
            }
        }
    }

    // new version of analysis
    if (fTrack && is_carbon && is_alpha && fPairs && chiHex < 1.e+36) // 17364.240
    {
        fNeventstrack += 1;

        fh_chiy_vs_chix_C->Fill(chiCx + chiHex);
        fh_chi2->Fill(chiHex);
        fh_chi2->Fill(chiCx);
        Double_t sq_chi = sqrt(chiHex * chiHex + chiCx * chiCx);
        fh_chiy_vs_chix_He_nc->Fill(chiHex, chiCx);

        LOG(DEBUG) << "Entering Pair analysis***" << endl;
        // if (chiHex < fcut_chiX && chiHey < fcut_chiY && chiCx < fcut_chiX && chiCy < fcut_chiY)
        //		&&((alpha.Px() > 0 && carbon.Px() < 0) || (alpha.Px() < 0 && carbon.Px() > 0)))
        if (alpha.Pz() > 0 && carbon.Pz() > 0)
        {
            //   if(abs(XHe-XC)*100. > 0.03) return;
            //   if(abs(YHe-YC)*100. > 0.03) return;

            //  if(abs(XHe-XC)*100. > 0.2) return;
            //  if(abs(YHe-YC)*100. > 0.2) return;

            pa = alpha.Vect();
            pc = carbon.Vect();

            // Calculate angle between alpha and C
            if (alpha.Pz() == 0 || carbon.Pz() == 0)
                return;
            /*
                        if (alpha.Pz() < 3900 || alpha.Pz() > 4800)
                           return;
                       if (carbon.Pz() < 12600 || carbon.Pz() > 13500)
                           return;
           */
            theta_26 = alpha.Angle(carbon.Vect()) * TMath::RadToDeg(); // angle alpha carbon (theta)
            Double_t phi_26 = (alpha.Phi() - carbon.Phi()) * TMath::RadToDeg();

            // if(theta_26 < 0.5) return;

            /* * Calculate invariant mass and relative energy * */
            m_inva = (alpha + carbon).M(); // invariant mass

            Erel = m_inva - mHe - mC; // relative Energy

            /*
                        Double_t costh26_test = (alpha.Px() * carbon.Px() + alpha.Py() * carbon.Py() + alpha.Pz() *
               carbon.Pz()) / (sqrt(alpha.Px() * alpha.Px() + alpha.Py() * alpha.Py() + alpha.Pz() * alpha.Pz()) *
                               sqrt(carbon.Px() * carbon.Px() + carbon.Py() * carbon.Py() + carbon.Pz() * carbon.Pz()));

                        Double_t test_theta = acos(costh26_test) * TMath::RadToDeg(); // opening angle
                        Double_t m_inv_test = sqrt(mHe * mHe + mC * mC + 2. * sqrt(pc.Mag() * pc.Mag() + mC * mC) *
                                              sqrt(pa.Mag() * pa.Mag() + mHe * mHe) -
                                              2. * pa.Mag() * pc.Mag() * cos(test_theta * TMath::DegToRad()));

                        cout<<"Checking minv: "<<m_inva<<", "<<m_inv_test<<endl;
            */

            /* * transfer to cm system and rotate * */

            oxygen = alpha + carbon; // excited oxygen
            // cout<<"gamma "<<oxygen.Gamma()<<endl;
            Double_t theta_16O = oxygen.Theta() * TMath::RadToDeg();
            Double_t phi_16O = oxygen.Phi() * TMath::RadToDeg();
            if (phi_16O < 0)
                phi_16O += 360.;

            // photon.SetPxPyPzE(oxygen.Px, oxygen.Py, 0., sqrt(pow(oxygen.Px, 2) + pow(oxygen.Py, 2)));

            // transfer to cm system and make some rotations
            // boost them to centre of mass
            alpha_cm = alpha;
            carbon_cm = carbon;
            oxygen_cm = oxygen;
            alpha_cm.Boost(-oxygen.BoostVector());
            carbon_cm.Boost(-oxygen.BoostVector());
            oxygen_cm.Boost(-oxygen.BoostVector());
            // cout << " px: " << alpha_cm.Px() + carbon_cm.Px() << "  py:  " << alpha_cm.Py() + carbon_cm.Py()
            //     << "  pz:  " << alpha_cm.Pz() + carbon_cm.Pz() << endl;
            /*
                        cout << "lab " << alpha.Px() << "  " << alpha.Py() << "  " << alpha.Pz() << endl;
                        cout << "lab " << carbon.Px() << "  " << carbon.Py() << "  " << carbon.Pz() << endl;
                        cout << "lab " << oxygen.Px() << "  " << oxygen.Py() << "  " << oxygen.Pz() << endl;

                        cout << "cm " << alpha_cm.Px() << "  " << alpha_cm.Py() << "  " << alpha_cm.Pz() << endl;
                        cout << "cm " << carbon_cm.Px() << "  " << carbon_cm.Py() << "  " << carbon_cm.Pz() << endl;
                        cout << "cm " << oxygen_cm.Px() << "  " << oxygen_cm.Py() << "  " << oxygen_cm.Pz() << endl;
            */

            Double_t phi = (alpha_cm.Phi() - carbon_cm.Phi()) * TMath::RadToDeg();
            Double_t theta = alpha_cm.Angle(carbon_cm.Vect()) * TMath::RadToDeg();

            phi_bc_cm = (oxygen.Phi() - alpha_cm.Phi()) * TMath::RadToDeg();

            if (phi_bc_cm < 0)
                phi_bc_cm += 360.;

            Double_t tthetaxC = pCx / pCz;    // * 1000.;
            Double_t tthetayC = pCy / pCz;    // * 1000.;
            Double_t tthetaxHe = pHex / pHez; // * 1000.;
            Double_t tthetayHe = pHey / pHez; // * 1000.;

            Double_t thetaxC = atan(tthetaxC) * 1000.;
            Double_t thetayC = atan(tthetayC) * 1000.;
            Double_t thetaxHe = atan(tthetaxHe) * 1000.;
            Double_t thetayHe = atan(tthetayHe) * 1000.;

            Double_t Ea = mHe * sqrt(1. + (pa.Mag() / mHe) * (pa.Mag() / mHe)) / 1000.;
            Double_t Ec = mC * sqrt(1. + (pc.Mag() / mC) * (pc.Mag() / mC)) / 1000.;
            /*
                   cout<<setprecision(10)<<"12C einzelne: "<<pCx<<" - "<<carbon.Px()<<"; "<<pCy<<" - "<<carbon.Py()<<";
               "<<pCz<<" - "<<carbon.Pz()<<endl; cout<<setprecision(10)<<"12C total: "<<pc.Mag()<<" -
               "<<sqrt(pCx*pCx+pCy*pCy+pCz*pCz)<<endl; cout<<setprecision(10)<<"4He einzelne: "<<pHex<<" -
               "<<alpha.Px()<<"; "<<pHey<<" - "<<alpha.Py()<<"; "<<pHez<<" - "<<alpha.Pz()<<endl;
                   cout<<setprecision(10)<<"4He total: "<<pa.Mag()<<" - "<<sqrt(pHex*pHex+pHey*pHey+pHez*pHez)<<endl;
                   cout<<setprecision(10)<<"Sum of two: "<<(pa+pc).Mag()<<",
               "<<sqrt((pHex+pCx)*(pHex+pCx)+(pHey+pCy)*(pHey+pCy)+(pHez+pCz)*(pHez+pCz))<<", "<<
                                       sqrt(pa.Mag()*pa.Mag()+pc.Mag()*pc.Mag()+2. * pa.Mag() * pc.Mag() * cos(theta_26
               * TMath::DegToRad()))<<endl;
             */

            fh_target_xy_nc->Fill(XHe, YHe);
            fh_target_xy_nc->Fill(XC, YC);

            fh_mass_nc->Fill(mtrackC, mtrackHe);
            fh_energy_nc->Fill(Ec, Ea);
            fh_px_He_nc->Fill(alpha.Px());
            fh_py_He_nc->Fill(alpha.Py());
            fh_pz_He_nc->Fill(alpha.Pz());
            fh_p_He_nc->Fill(pa.Mag());

            fh_px_C_nc->Fill(carbon.Px());
            fh_py_C_nc->Fill(carbon.Py());
            fh_pz_C_nc->Fill(carbon.Pz());
            fh_p_C_nc->Fill(pc.Mag());

            fh_px_px_nc->Fill(alpha.Px(), carbon.Px());
            fh_py_py_nc->Fill(alpha.Py(), carbon.Py());
            fh_pz_pz_nc->Fill(alpha.Pz(), carbon.Pz());

            fh_psum_nc->Fill((pa + pc).Mag());
            fh_theta26_nc->Fill(theta_26);
            fh_Erel_nc->Fill(Erel); // relative Energy plot

            fh_psum_vs_event_nc->Fill(fNeventstrack, (pa + pc).Mag());
            fh_psum_vs_chi->Fill(sq_chi, (pa + pc).Mag());

            fh_theta26_vs_chi->Fill(sq_chi, theta_26);
            fh_Erel_vs_theta26_nc->Fill(theta_26, Erel);
            fh_Erel_vs_chi->Fill(sq_chi, Erel);
            fh_phi26_vs_chi->Fill(sq_chi, phi_26);
            fh_psum_vs_theta26_nc->Fill(theta_26, (pa + pc).Mag());

            if ((alpha.Px() > 0. && carbon.Px() < 0.) || (alpha.Px() < 0. && carbon.Px() > 0.))
            {

                fh_ErelB_nc->Fill(Erel);
            }

            fh_theta_bc_cm_nc->Fill(alpha_cm.Theta() * TMath::RadToDeg());
            fh_phi_bc_cm_nc->Fill(phi_bc_cm);

            fh_chiy_vs_chix_He->Fill(chiHex, chiCx);

            is_tracked = false;

            //  Double_t Eacheck = 17.1514-0.665768*Ec;
            Double_t Eacheck = -75.6469 + 10.2155 * Ec - 0.318913 * Ec * Ec;
            Bool_t becut = false;
            if (fGraphCuts)
            {
                if (cut_EHe_EC->IsInside(Ec, Ea))
                    becut = true;
            }
            else
            {
                if (abs(Ea - Eacheck) < 0.12)
                    becut = true;
            }
            // Dima's tracker only one chi2;
            if (sq_chi < (fcut_chiX + fcut_chiY) && becut && int(mtrackC + 0.5) == 12 &&
                int(mtrackHe + 0.5) == 4) //  abs(carbon.Pz()+alpha.Pz()-17390.) < 400. )
            {
                //  if ((alpha.Px() > 0. && carbon.Px() > 0.) || (alpha.Px() < 0. && carbon.Px() < 0.)) return;
                //  if ((alpha.Py() > 0. && carbon.Py() > 0.) || (alpha.Py() < 0. && carbon.Py() < 0.)) return;
                //   if((pa + pc).Mag() < 17352. || (pa + pc).Mag()> 17432.) return;

                is_tracked = true;

                fh_energy->Fill(Ec, Ea);

                fh_minv->Fill(Ec / 1000. + Ea / 1000.);

                //  fh_chiy_vs_chix_He->Fill(chiHex,chiCx);
                fh_target_xy->Fill(XHe, YHe);
                fh_target_xy->Fill(XC, YC);

                fh_target_xx->Fill(XHe - XC);
                fh_target_yy->Fill(YHe - YC);

                fh_mass->Fill(mtrackC, mtrackHe);

                fh_psum_vs_event->Fill(fNeventstrack, (pa + pc).Mag());

                fh_px_He->Fill(alpha.Px());
                fh_py_He->Fill(alpha.Py());
                fh_pz_He->Fill(alpha.Pz());
                fh_p_He->Fill(pa.Mag());

                fh_px_C->Fill(carbon.Px());
                fh_py_C->Fill(carbon.Py());
                fh_pz_C->Fill(carbon.Pz());
                fh_p_C->Fill(pc.Mag());

                fh_px_px->Fill(alpha.Px(), carbon.Px());
                fh_py_py->Fill(alpha.Py(), carbon.Py());
                fh_pz_pz->Fill(alpha.Pz(), carbon.Pz());
                fh_p_p->Fill(pa.Mag(), pc.Mag());

                fh_psum->Fill((pa + pc).Mag());
                fh_pzsum->Fill(alpha.Pz() + carbon.Pz());

                fh_theta26->Fill(theta_26);
                fh_Erel->Fill(Erel); // relative Energy plot

                fh_phi_bc_cm->Fill(phi_bc_cm);

                if (alpha.Px() < 0. && carbon.Px() < 0.)
                    fh_ErelR->Fill(Erel);
                if (alpha.Px() > 0. && carbon.Px() > 0.)
                    fh_ErelL->Fill(Erel);

                fh_Erel_vs_theta26->Fill(theta_26, Erel);
                fh_pHe_vs_theta26->Fill(pa.Mag(), theta_26 * TMath::DegToRad() * 1000.);
                fh_psum_vs_theta26->Fill(theta_26, (pa + pc).Mag());

                if ((alpha.Px() > 0. && carbon.Px() < 0.) || (alpha.Px() < 0. && carbon.Px() > 0.))
                {
                    // if (theta_26 > 0.7)
                    fh_ErelB->Fill(Erel);
                }
                fh_phibcm_vs_px->Fill(carbon.Px(), phi_bc_cm);
                fh_phibcm_vs_py->Fill(carbon.Py(), phi_bc_cm);
                fh_phibcm_vs_pz->Fill(carbon.Pz(), phi_bc_cm);
                fh_phibcm_vs_px->Fill(alpha.Px(), phi_bc_cm);
                fh_phibcm_vs_py->Fill(alpha.Py(), phi_bc_cm);
                fh_phibcm_vs_pz->Fill(alpha.Pz(), phi_bc_cm);
                fh_phibcm_vs_psum->Fill((pa + pc).Mag(), phi_bc_cm);
                fh_Erel_vs_psum->Fill((pa + pc).Mag(), Erel);
                fh_theta_16O->Fill(theta_16O); // theta oxygen
                fh_phi_16O->Fill(phi_16O);

                fh_phi26_cm->Fill(phi);     // phi between alpha + carbon in cm
                fh_theta26_cm->Fill(theta); // angle between alpha + carbon in cm

                fh_theta_4He_cm->Fill(alpha_cm.Theta() * TMath::RadToDeg());
                fh_phi_4He_cm->Fill(alpha_cm.Phi() * TMath::RadToDeg());
                fh_theta_12C_cm->Fill(carbon_cm.Theta() * TMath::RadToDeg());
                fh_phi_12C_cm->Fill(carbon_cm.Phi() * TMath::RadToDeg());

                fh_theta_bc_cm->Fill(alpha_cm.Theta() * TMath::RadToDeg());

                fh_phi_bc_cm_polar->Fill(phi_bc_cm, 1);
                fh_Erel_vs_phibc->Fill(Erel, phi_bc_cm);

                fh_thetax_px->Fill(pHex, thetaxHe);
                fh_thetax_px->Fill(pCx, thetaxC);

                fh_thetay_py->Fill(pHey, thetayHe);
                fh_thetay_py->Fill(pCy, thetayC);

                fh_thetax_py->Fill(pHey, thetaxHe);
                fh_thetax_py->Fill(pCy, thetaxC);

                fh_thetay_px->Fill(pHex, thetayHe);
                fh_thetay_px->Fill(pCx, thetayC);
            }

            /*
                        // rotation in oxygen direction
                        //TLorentzRotation transform;
                        //transform.RotateX(oxygen.Px()).RotateY(oxygen.Py()).RotateZ(0.);
                        TVector3 transform = oxygen.Vect(); /// TODO: check if it needs to be normed
                        alpha *= transform; // apply transformation vector
                        carbon *= transform;

                        if(abs(alpha.Phi()*TMath::RadToDeg()) > 88 && abs(alpha.Phi()*TMath::RadToDeg()) < 92)
                        {
                        fh_target_xy->Fill(XHe * 100., YHe * 100.);

                        fh_px_He->Fill(pHex);
                        fh_py_He->Fill(pHey);
                        fh_pz_He->Fill(pHez);

                        fh_px_C->Fill(pCx);
                        fh_py_C->Fill(pCy);
                        fh_pz_C->Fill(pCz);

                            fh_theta26->Fill(theta26); ///angle alpha carbon (theta)
                            fh_theta_16O->Fill(theta16O); // theta oxygen


                            fh_Erel->Fill(Erel); // relative Energy plot
                        fh_cs_cd->Fill(Erel,Erel);


            */

        } // end if chi2
    }     // end if trackHits>1

    if (fTrack && is_oxygen && !fPairs)
    {

        fh_chi2->Fill(chix);
        fh_psum_vs_chi->Fill(chix, Pf_tot);
        fh_target_xy_nc->Fill(Xf, Yf);
        fh_px_O_nc->Fill(Pxf);
        fh_py_O_nc->Fill(Pyf);
        fh_pz_O_nc->Fill(Pzf);
        fh_p_O_nc->Fill(Pf_tot);
        fh_psum_vs_event_nc->Fill(fNEvents, Pf_tot);
        Double_t Radius = sqrt(Xf * Xf + Yf * Yf);
        fh_psum_vs_r_nc->Fill(Radius, Pf_tot);

        if (chix < fcut_chiX && chiy < fcut_chiY)
        {
            is_tracked = true;

            fh_target_xy->Fill(Xf, Yf);

            fh_px_O->Fill(Pxf);
            fh_py_O->Fill(Pyf);
            fh_pz_O->Fill(Pzf);
            fh_p_O->Fill(Pf_tot);
            fh_psum_vs_event->Fill(fNEvents, Pf_tot);
            fh_pz_vs_r->Fill(Radius, Pzf);
            fh_px_vs_r->Fill(Radius, Pxf);
            fh_py_vs_r->Fill(Radius, Pyf);

        } // end if chi2
    }     // end if trackHits>1

    if (fSimu && nHitsMCTrack > 0)
    {
        if (fPairs) //&& is_tracked)
        {

            paMC = alphaMC.Vect();
            pcMC = carbonMC.Vect();
            Double_t theta_26MC = alphaMC.Angle(carbonMC.Vect()) * TMath::RadToDeg(); /// angle alpha carbon (theta)
            fh_theta26_simu->Fill(theta_26MC);
            // Double_t randtheta = (std::rand() / (float)RAND_MAX) - 0.5;
            //** Calculate invariant mass and relative energy **/
            m_invaMC = (alphaMC + carbonMC).M(); // invariant mass
            ErelMC = m_invaMC - mHe - mC;        // relative Energy
            fh_Erel_vs_thetaMC->Fill(theta_26MC, ErelMC);

            /* * transfer to cm system and rotate * */

            oxygenMC = alphaMC + carbonMC; // excited oxygen
            // cout<<"gamma "<<oxygen.Gamma()<<endl;
            Double_t thetaMC_16O = oxygenMC.Theta() * TMath::RadToDeg();
            Double_t phiMC_16O = oxygenMC.Phi() * TMath::RadToDeg();
            if (phiMC_16O < 0)
                phiMC_16O += 360.;

            // transfer to cm system and make some rotations
            // boost them to centre of mass
            alphaMC_cm = alphaMC;
            carbonMC_cm = carbonMC;
            oxygenMC_cm = oxygenMC;
            alphaMC_cm.Boost(-oxygenMC.BoostVector());
            carbonMC_cm.Boost(-oxygenMC.BoostVector());
            oxygenMC_cm.Boost(-oxygenMC.BoostVector());

            Double_t phiMC = (alphaMC_cm.Phi() - carbonMC_cm.Phi()) * TMath::RadToDeg();
            Double_t thetaMC = alphaMC_cm.Angle(carbonMC_cm.Vect()) * TMath::RadToDeg();

            fh_psum_MC->Fill((paMC + pcMC).Mag());
            fh_pzsum_MC->Fill(alphaMC.Pz() + carbonMC.Pz());

            phiMC_bc_cm = (oxygenMC.Phi() - alphaMC_cm.Phi()) * TMath::RadToDeg();
            if (phiMC_bc_cm < 0)
                phiMC_bc_cm += 360.;
            fh_phiMC_bc_cm->Fill(phiMC_bc_cm);

            fh_phiMC_bc_cm_polar->Fill(phiMC_bc_cm, 1);

            fh_minv_simu->Fill(m_invaMC);
            fh_Erel_vs_phibcMC->Fill(ErelMC, phiMC_bc_cm);
            fh_Erel_simu->Fill(ErelMC);

            fh_px_px_mc->Fill(pHex_mc, pCx_mc);
            fh_py_py_mc->Fill(pHey_mc, pCy_mc);
            fh_pz_pz_mc->Fill(pHez_mc, pCz_mc);

            fh_phibcm_vs_px_MC->Fill(pCx_mc, phiMC_bc_cm);
            fh_phibcm_vs_py_MC->Fill(pCy_mc, phiMC_bc_cm);
            fh_phibcm_vs_pz_MC->Fill(pCz_mc, phiMC_bc_cm);
            fh_phibcm_vs_px_MC->Fill(pHex_mc, phiMC_bc_cm);
            fh_phibcm_vs_py_MC->Fill(pHey_mc, phiMC_bc_cm);
            fh_phibcm_vs_pz_MC->Fill(pHez_mc, phiMC_bc_cm);
            fh_phibcm_vs_psum_MC->Fill((paMC + pcMC).Mag(), phiMC_bc_cm);

            if (is_tracked)
            {
                fh_dErel->Fill((ErelMC - Erel) / ErelMC * 100.);
                fh_dtheta->Fill(theta_26MC - theta_26);

                fh_dx->Fill(XHe_mc - XHe);
                fh_dx->Fill(XC_mc - XC);
                fh_dy->Fill(YHe_mc - YHe);
                fh_dy->Fill(YC_mc - YC);
                fh_target_xy_MC->Fill(XHe_mc, YHe_mc);
                fh_target_xy_MC->Fill(XC_mc, YC_mc);

                //   cout << "X-momenta: MC_He: " << pHex_mc << ", tracked He: " << pHex << endl;
                //   cout << "X-momenta: MC_C: " << pCx_mc << ", tracked C: " << pCx << endl;
                //   cout << "Y-momenta: MC_He: " << pHey_mc << ", tracked He: " << pHey << endl;
                //   cout << "Y-momenta: MC_C: " << pCy_mc << ", tracked C: " << pCy << endl;
                //   cout << "Momenta: MC_He: " << paMC.Mag() << ", tracked He: " << pa.Mag() << endl;
                //   cout << "Momenta: MC_C: " << pcMC.Mag() << ", tracked C: " << pc.Mag() << endl;

                fh_dpxC->Fill((pCx_mc - pCx) / pCx_mc * 100.);
                fh_dpxHe->Fill((pHex_mc - pHex) / pHex_mc * 100.);
                fh_dpyC->Fill((pCy_mc - pCy) / pCy_mc * 100.);
                fh_dpyHe->Fill((pHey_mc - pHey) / pHey_mc * 100.);
                fh_dpzC->Fill((pCz_mc - pCz) / pCz_mc * 100.);
                fh_dpzHe->Fill((pHez_mc - pHez) / pHez_mc * 100.);
                fh_dpC->Fill((pcMC.Mag() - pc.Mag()) / pcMC.Mag() * 100.);
                fh_dpHe->Fill((paMC.Mag() - pa.Mag()) / paMC.Mag() * 100.);
                fh_dp->Fill(((paMC + pcMC).Mag() - (pa + pc).Mag()) / ((paMC + pcMC).Mag()) * 100.);
                fh_Erel_dErel->Fill(ErelMC - Erel, (paMC + pcMC).Mag() - (pa + pc).Mag());

                Double_t diff =
                    sqrt(pCx_mc * pCx_mc + pCy_mc * pCy_mc + pCz_mc * pCz_mc) - sqrt(pCx * pCx + pCy * pCy + pCz * pCz);
                //	cout<<"Difference: "<<pcMC.Mag() - pc.Mag()<<", "<<diff<<endl;

                fh_x_dpx->Fill((pCx_mc - pCx) / pCx_mc * 100., pCx);
                fh_x_dpx->Fill((pHex_mc - pHex) / pHex_mc * 100., pHex);
                fh_y_dpy->Fill((pCy_mc - pCy) / pCy_mc * 100., pCy);
                fh_y_dpy->Fill((pHey_mc - pHey) / pHey_mc * 100., pHey);
                fh_dpy_dpx->Fill((pCx_mc - pCx) / pCx_mc * 100., (pCy_mc - pCy) / pCy_mc * 100.);
                fh_dpy_dpx->Fill((pHex_mc - pHex) / pHex_mc * 100., (pHey_mc - pHey) / pHey_mc * 100.);

                // Comparison with simulation
                Double_t tthetaxC_simu = pCx_mc / pCz_mc;    // * 1000.;
                Double_t tthetayC_simu = pCy_mc / pCz_mc;    // * 1000.;
                Double_t tthetaxHe_simu = pHex_mc / pHez_mc; // * 1000.;
                Double_t tthetayHe_simu = pHey_mc / pHez_mc; // *  1000.;

                Double_t thetaxC_simu = atan(tthetaxC_simu) * 1000.;
                Double_t thetayC_simu = atan(tthetayC_simu) * 1000.;
                Double_t thetaxHe_simu = atan(tthetaxHe_simu) * 1000.;
                Double_t thetayHe_simu = atan(tthetayHe_simu) * 1000.;

                fh_thetax_dpx_He->Fill((pHex_mc - pHex), thetaxHe_simu);
                fh_thetax_dpx_C->Fill((pCx_mc - pCx), thetaxC_simu);

                fh_thetay_dpy_He->Fill((pHey_mc - pHey), thetayHe_simu);
                fh_thetay_dpy_C->Fill((pCy_mc - pCy), thetayC_simu);

                fh_thetax_dpy->Fill((pHey_mc - pHey), thetaxHe_simu);
                fh_thetax_dpy->Fill((pCy_mc - pCy), thetaxC_simu);

                fh_thetay_dpx->Fill((pHex_mc - pHex), thetayHe_simu);
                fh_thetay_dpx->Fill((pCx_mc - pCx), thetayC_simu);
            }
        }
        else
        {

            if (is_tracked)
            {

                fh_dx->Fill((Xf_mc - Xf));
                fh_dy->Fill((Yf_mc - Yf));

                fh_dpxO->Fill((Pxf_mc - Pxf) / Pxf_mc * 100.);
                fh_dpyO->Fill((Pyf_mc - Pyf) / Pyf_mc * 100.);
                fh_dpzO->Fill((Pzf_mc - Pzf) / Pzf_mc * 100.);
                // fh_dpO->Fill((Pf_tot_mc - Pf_tot) / Pf_tot_mc * 100.);

                fh_x_dpx->Fill((Pxf_mc - Pxf) / Pxf_mc * 100., Xf);
                fh_y_dpy->Fill((Pyf_mc - Pyf) / Pyf_mc * 100., Yf);
                // fh_dpy_dpx->Fill((Pxf_mc - Pxf) / Pxf_mc * 100., (Pyf_mc - Pyf) / Pyf_mc * 100.);
            }
        }
    }
    if (fMappedItemsCalifa && is_tracked)
    {
        // CALIFA
        auto detCalifa = fMappedItemsCalifa;
        Int_t nHitsCalifa = detCalifa->GetEntriesFast();
        // cout<<"Califa hits: "<<nHitsCalifa<<endl;

        for (Int_t ihit = 0; ihit < nHitsCalifa; ihit++)
        {
            R3BCalifaMappedData* hitCalifa = (R3BCalifaMappedData*)detCalifa->At(ihit);
            if (!hitCalifa)
                continue;

            Int_t Crystal = hitCalifa->GetCrystalId();
            Int_t Energy = hitCalifa->GetEnergy();
            // cout << "Califa: " << Crystal << " Energy: " << Energy << endl;
            if (Energy > 0)
            {
                fh_califa_energy->Fill(Crystal, Energy);
            }
        }
    }
    if (fCalItemsCalifa && is_tracked)
    {
        // CALIFA
        auto detCalifa = fCalItemsCalifa;
        Int_t nHitsCalifa = detCalifa->GetEntriesFast();
        // cout<<"Califa hits: "<<nHitsCalifa<<endl;

        for (Int_t ihit = 0; ihit < nHitsCalifa; ihit++)
        {
            R3BCalifaCrystalCalData* hitCalifa = (R3BCalifaCrystalCalData*)detCalifa->At(ihit);
            if (!hitCalifa)
                continue;

            Int_t Crystal = hitCalifa->GetCrystalId();
            Int_t Energy = hitCalifa->GetEnergy();
            // cout << "Califa: " << Crystal << " Energy: " << Energy << endl;
            if (Energy > 0)
            {
                fh_califa_calenergy->Fill(Crystal, Energy);
            }
        }
    }
    if (fHitItemsCalifa && is_tracked)
    {
        // CALIFA
        auto detCalifa = fHitItemsCalifa;
        Int_t nHitsCalifa = detCalifa->GetEntriesFast();
        // cout<<"Califa hits: "<<nHitsCalifa<<endl;

        for (Int_t ihit = 0; ihit < nHitsCalifa; ihit++)
        {
            R3BCalifaHitData* hitCalifa = (R3BCalifaHitData*)detCalifa->At(ihit);
            if (!hitCalifa)
                continue;

            Int_t Crystal = hitCalifa->GetClusterId();
            Int_t Energy = hitCalifa->GetEnergy();
            // cout << "Califa: " << Crystal << " Energy: " << Energy << endl;
            if (Energy > 0)
            {
                fh_califa_hitenergy->Fill(Crystal, Energy);
            }
        }
    }

    // second, old version of analysis
    if (is_carbon && is_alpha && 1 == 0)
    {
        fh_chiy_vs_chix_He->Fill(chiHex, chiHey);
        fh_chiy_vs_chix_C->Fill(chiCx + chiCy);

        if (chiHex < 1e10 && chiHey < 1e10 && chiCx < 1e10 && chiCy < 1e10)
        {
            fh_target_xy->Fill(XHe * 100., YHe * 100.);

            fh_px_He->Fill(pHex);
            fh_py_He->Fill(pHey);
            fh_pz_He->Fill(pHez);

            fh_px_C->Fill(pCx);
            fh_py_C->Fill(pCy);
            fh_pz_C->Fill(pCz);

            fh_px_px->Fill(pHex, pCx);
            fh_py_py->Fill(pHey, pCy);
            fh_pz_pz->Fill(pHez, pCz);

            // Calculate angle between alphs and C
            if (pCz == 0 || pHez == 0)
                return;

            costh26 = (pHex * pCx + pHey * pCy + pHez * pCz) /
                      (sqrt(pHex * pHex + pHey * pHey + pHez * pHez) * sqrt(pCx * pCx + pCy * pCy + pCz * pCz));

            theta_26 = acos(costh26) * 180. / 3.14159; // opening angle

            fh_theta26->Fill(theta_26);

            Double_t pHe = sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2));
            Double_t pC = sqrt(pow(pCx, 2) + pow(pCy, 2) + pow(pCz, 2));
            Double_t eHe = sqrt(pow(pHe, 2) + pow(mHe, 2)) - mHe;
            Double_t eC = sqrt(pow(pC, 2) + pow(mC, 2)) - mC;

            fh_p_He->Fill(pHe);
            fh_p_C->Fill(pC);

            Double_t Erela = sqrt(pow((mHe + mC + eHe + eC), 2) - pow(pHe, 2) - pow(pC, 2) -
                                  2 * pHe * pC * cos(theta_26 * 3.1415 / 180.)) -
                             (mHe + mC); // Erel

            Double_t ptota = sqrt(pow(pHe, 2) + pow(pC, 2) + 2 * pHe * pC * cos(theta_26 * 3.1415 / 180.));

            Double_t m_inv = sqrt(mHe * mHe + mC * mC + 2. * sqrt(pC * pC + mC * mC) * sqrt(pHe * pHe + mHe * mHe) -
                                  2. * pHe * pC * cos(theta_26 * 3.1415 / 180.));
            Double_t Erelb = m_inv - mHe - mC;

            fh_Erel->Fill(Erela);
            if (pHex < 0. && pCx < 0.)
                fh_ErelR->Fill(Erela);
            if (pHex > 0. && pCx > 0.)
                fh_ErelL->Fill(Erela);
            if ((pHex > 0. && pCx < 0.) || (pHex < 0. && pCx > 0.))
                fh_ErelB->Fill(Erela);

            //	fh_dErel_vs_x->Fill(xTest * 100., Erela - 1.);
            //	fh_dErel_vs_y->Fill(yTest * 100., Erela - 1.);
            // cout << "Theta 26: " << theta_26 << " Erel: " << Erela << " " << Erelb << endl;

            // transfer to cm system and make some rotations

            px = pHex + pCx;
            py = pHey + pCy;
            pz = pHez + pCz;

            if (pz != 0.)
            {
                theta_16 = atan2(sqrt(pow(px, 2) + pow(py, 2)), pz) * 180. / 3.14159265;
            }

            fh_theta_16O->Fill(theta_16);

            Double_t n2[10], eaa, vaa, g2, v2[10], lt2[10][10];
            Double_t tmp1[10], tmp2[10], tmp3[10], tmp4[10];
            Double_t pxO_1, pyO_1, pzO_1, esum_3;
            Double_t psum_1, thcm, phcm, tanphicm;

            esum_3 = sqrt(pow((sqrt(pow(pHe, 2) + pow(pC, 2) + 2 * pHe * pC * cos(theta_26 * 3.14159 / 180.))), 2) +
                          pow(mO, 2)) -
                     mO;

            psum_1 = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
            thcm = acos(pz / psum_1);
            tanphicm = py / px;

            if (px > 0)
            {
                phcm = atan(tanphicm) * 180. / 3.14159265;
            }

            if (px < 0. && py > 0.)
            {
                phcm = (atan(tanphicm) * 180. / 3.14159265) + 180.;
            }

            if (px < 0. && py < 0.)
            {
                phcm = (atan(tanphicm) * 180. / 3.14159265) - 180.;
            }

            if (phcm < 0)
            {
                phcm = phcm + 360.;
            }
            phcm = phcm * 3.14159265 / 180.;

            //            fh_theta_16O_cm->Fill(thcm);
            //            fh_phi_16O_cm->Fill(phcm);

            n2[1] = px / ptota;
            n2[2] = py / ptota;
            n2[3] = pz / ptota;

            //       eaa   = esum_3;
            //     velocity of 6Li
            //       vaa=sqrt(eaa*eaa+2*mO*eaa)/(eaa+mO);

            //       cout<< " va1 "<< vaa<<endl;

            eaa = sqrt(mO * mO + psum_1 * psum_1);
            vaa = psum_1 / eaa;

            //       cout<< " va2 "<< vaa<<endl;

            g2 = 1.0 / sqrt(1.0 - vaa * vaa);

            v2[1] = vaa * n2[1];
            v2[2] = vaa * n2[2];
            v2[3] = vaa * n2[3];

            // Lorentz transformation
            lt2[1][1] = g2;
            lt2[1][2] = -1. * g2 * v2[1];
            lt2[1][3] = -1. * g2 * v2[2];
            lt2[1][4] = -1. * g2 * v2[3];
            lt2[2][1] = lt2[1][2];
            lt2[2][2] = 1.0 + (g2 - 1.0) * n2[1] * n2[1];
            lt2[2][3] = (g2 - 1.0) * n2[1] * n2[2];
            lt2[2][4] = (g2 - 1.0) * n2[1] * n2[3];
            lt2[3][1] = lt2[1][3];
            lt2[3][2] = lt2[2][3];
            lt2[3][3] = 1.0 + (g2 - 1.0) * n2[2] * n2[2];
            lt2[3][4] = (g2 - 1.0) * n2[2] * n2[3];
            lt2[4][1] = lt2[1][4];
            lt2[4][2] = lt2[2][4];
            lt2[4][3] = lt2[3][4];
            lt2[4][4] = 1.0 + (g2 - 1.0) * n2[3] * n2[3];

            //      four momentum of nucleus b
            tmp1[1] = sqrt(mHe * mHe + pHe * pHe);
            tmp1[2] = pHex;
            tmp1[3] = pHey;
            tmp1[4] = pHez;

            //     four momentum of nucleus c
            tmp2[1] = sqrt(mC * mC + pC * pC);
            tmp2[2] = pCx;
            tmp2[3] = pCy;
            tmp2[4] = pCz;

            //     four momentum of nucleus a for a check of the transformation
            tmp4[1] = sqrt(mO * mO + psum_1 * psum_1);
            tmp4[2] = px;
            tmp4[3] = py;
            tmp4[4] = pz;
            //       cout <<"compare " << tmp4[1]<<"  "<< esum_3+mO<<endl;

            /*
                            cout <<"in " << tmp4[1]<<endl;
                            cout <<"in " << tmp4[2]<<endl;
                            cout <<"in " << tmp4[3]<<endl;
                            cout <<"in " << tmp4[4]<<endl;
            */
            for (Int_t i1 = 1; i1 < 5; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 5; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp1[i2];
                }
            }

            Double_t pHex_1, pHey_1, pHez_1, pHe_1, pCx_1, pCy_1, pCz_1, pC_1;

            pHex_1 = tmp3[2];
            pHey_1 = tmp3[3];
            pHez_1 = tmp3[4];

            for (Int_t i1 = 1; i1 < 5; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 5; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp2[i2];
                }
            }

            pCx_1 = tmp3[2];
            pCy_1 = tmp3[3];
            pCz_1 = tmp3[4];

            for (Int_t i1 = 1; i1 < 5; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 5; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp4[i2];
                }
            }

            pxO_1 = tmp3[2];
            pyO_1 = tmp3[3];
            pzO_1 = tmp3[4];

            /*
                            cout <<"out " << tmp3[1]<<endl;
                            cout <<"out " << tmp3[2]<<endl;
                            cout <<"out " << tmp3[3]<<endl;
                            cout <<"out " << tmp3[4]<<endl;
            */

            Double_t pHex_cm, pHey_cm, pHez_cm, pHe_cm;
            Double_t pCx_cm, pCy_cm, pCz_cm, pC_cm;
            Double_t thet_cm_He, thet_cm_C, thet;
            Double_t theta_cm_He, theta_cm_C, theta;
            Double_t pxyd_cm, costhetacm_He, tanphi_d, phicm_He, pxyhe_cm;
            Double_t costhetacm_C, tanphi_he, phicm_C, theta_cm_26;
            Double_t phi_cm_26, pxO_cm, pyO_cm, pzO_cm;

            //     Vector for rotation
            n2[1] = py / sqrt(px * px + py * py);
            n2[2] = -px / sqrt(px * px + py * py);
            n2[3] = 0;

            //       cout <<"winkel " << thcm<<"  "<<theta6/180.*3.14159265<<endl;

            lt2[1][1] = cos(thcm) + n2[1] * n2[1] * (1. - cos(thcm));
            lt2[1][2] = n2[1] * n2[2] * (1. - cos(thcm)) - n2[3] * sin(thcm);
            lt2[1][3] = n2[1] * n2[3] * (1. - cos(thcm)) + n2[2] * sin(thcm);
            lt2[2][1] = n2[2] * n2[1] * (1. - cos(thcm)) + n2[3] * sin(thcm);
            lt2[2][2] = cos(thcm) + n2[2] * n2[2] * (1. - cos(thcm));
            lt2[2][3] = n2[2] * n2[3] * (1. - cos(thcm)) - n2[1] * sin(thcm);
            lt2[3][1] = n2[3] * n2[1] * (1. - cos(thcm)) - n2[2] * sin(thcm);
            lt2[3][2] = n2[3] * n2[2] * (1. - cos(thcm)) + n2[1] * sin(thcm);
            lt2[3][3] = cos(thcm) + n2[3] * n2[3] * (1. - cos(thcm));

            //      four momentum of nucleus b
            tmp1[1] = pHex_1;
            tmp1[2] = pHey_1;
            tmp1[3] = pHez_1;

            //     four momentum of nucleus c
            tmp2[1] = pCx_1;
            tmp2[2] = pCy_1;
            tmp2[3] = pCz_1;

            //     four momentum of nucleus a for a check of the transformation
            tmp4[1] = px;
            tmp4[2] = py;
            tmp4[3] = pz;
            //       cout <<"compare " << tmp4[1]<<"  "<< esum_3+mO<<endl;

            /*
                            cout <<"in1 " << tmp4[1]<<endl;
                            cout <<"in1 " << tmp4[2]<<endl;
                            cout <<"in1 " << tmp4[3]<<endl;
            */
            Double_t pHex_2, pHey_2, pHez_2, pHe_2, pCx_2, pCy_2, pCz_2, pC_2;
            Double_t pxO_2, pyO_2, pzO_2;

            for (Int_t i1 = 1; i1 < 4; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 4; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp1[i2];
                }
            }

            pHex_cm = tmp3[1];
            pHey_cm = tmp3[2];
            pHez_cm = tmp3[3];

            for (Int_t i1 = 1; i1 < 4; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 4; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp2[i2];
                }
            }

            pCx_cm = tmp3[1];
            pCy_cm = tmp3[2];
            pCz_cm = tmp3[3];

            for (Int_t i1 = 1; i1 < 4; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 4; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp4[i2];
                }
            }

            pxO_cm = tmp3[1];
            pyO_cm = tmp3[2];
            pzO_cm = tmp3[3];

            /*
                            cout <<"out1 " << tmp3[1]<<endl;
                            cout <<"out1 " << tmp3[2]<<endl;
                            cout <<"out1 " << tmp3[3]<<endl;
            */

            // 4He fragment

            pHe_cm = sqrt(pHex_cm * pHex_cm + pHey_cm * pHey_cm + pHez_cm * pHez_cm);

            costhetacm_He = pHez_cm / pHe_cm;
            theta_cm_He = acos(costhetacm_He) * 180. / 3.14159265;
            tanphi_d = pHey_cm / pHex_cm;

            if (pHex_cm > 0)
            {
                phicm_He = atan(tanphi_d) * 180. / 3.14159265;
            }

            if (pHex_cm < 0. && pHey_cm > 0.)
            {
                phicm_He = (atan(tanphi_d) * 180. / 3.14159265) + 180.;
            }

            if (pHex_cm < 0. && pHey_cm < 0.)
            {
                phicm_He = (atan(tanphi_d) * 180. / 3.14159265) - 180.;
            }

            if (phicm_He < 0)
            {
                phicm_He = phicm_He + 360.;
            }
            fh_theta_4He_cm->Fill(theta_cm_He);
            fh_phi_4He_cm->Fill(phicm_He);

            // 12C fragment

            pC_cm = sqrt(pCx_cm * pCx_cm + pCy_cm * pCy_cm + pCz_cm * pCz_cm);

            costhetacm_C = pCz_cm / pC_cm;
            theta_cm_C = acos(costhetacm_C) * 180. / 3.14159265;
            tanphi_he = pCy_cm / pCx_cm;

            if (pCx_cm > 0.)
            {
                phicm_C = atan(tanphi_he) * 180. / 3.14159265;
            }

            if (pCx_cm < 0. && pCy_cm > 0.)
            {
                phicm_C = (atan(tanphi_he) * 180. / 3.14159265) + 180.;
            }

            if (pCx_cm < 0. && pCy_cm < 0.)
            {
                phicm_C = (atan(tanphi_he) * 180. / 3.14159265) - 180.;
            }

            if (phicm_C < 0.)
            {
                phicm_C = phicm_C + 360.;
            }
            fh_theta_12C_cm->Fill(theta_cm_C);
            fh_phi_12C_cm->Fill(phicm_C);

            theta_cm_26 = theta_cm_He + theta_cm_C;
            phi_cm_26 = phicm_He - phicm_C;

            fh_theta26_cm->Fill(theta_cm_26);
            fh_phi26_cm->Fill(phi_cm_26);

            fh_px_px_mc->Fill(pHex_cm, pCx_cm);
            fh_py_py_mc->Fill(pHey_cm, pCy_cm);
            fh_pz_pz_mc->Fill(pHez_cm, pCz_cm);

        } // end if chi2
    }     // end if trackHits>1

    fNEvents += 1;

    return;
}

void R3BGlobalAnalysisS494::FinishEvent()
{
    fTrack->Clear();
    if (fMCTrack)
        fMCTrack->Clear();
}

void R3BGlobalAnalysisS494::FinishTask()
{
    fh_chi2->Write();
    fh_psum_vs_chi->Write();
    fh_Erel_vs_chi->Write();
    fh_phi26_vs_chi->Write();
    fh_psum_vs_theta26->Write();
    fh_chiy_vs_chix_He->Write();
    fh_chiy_vs_chix_He_nc->Write();
    fh_chiy_vs_chix_C->Write();
    fh_theta26_vs_chi->Write();

    fh_target_xy_nc->Write();
    fh_px_He_nc->Write();
    fh_py_He_nc->Write();
    fh_pz_He_nc->Write();
    fh_p_He_nc->Write();
    fh_px_C_nc->Write();
    fh_py_C_nc->Write();
    fh_pz_C_nc->Write();
    fh_p_C_nc->Write();
    fh_px_O_nc->Write();
    fh_py_O_nc->Write();
    fh_pz_O_nc->Write();
    fh_p_O_nc->Write();
    fh_px_px_nc->Write();
    fh_py_py_nc->Write();
    fh_pz_pz_nc->Write();
    fh_psum_vs_event_nc->Write();
    fh_psum_nc->Write();

    fh_theta26_nc->Write();
    fh_ErelB_nc->Write();
    fh_Erel_nc->Write();
    fh_theta_bc_cm_nc->Write();
    fh_phi_bc_cm_nc->Write();
    fh_Erel_vs_theta26_nc->Write();
    fh_psum_vs_r_nc->Write();

    fh_target_xy->Write();
    fh_target_xx->Write();
    fh_target_yy->Write();

    fh_px_He->Write();
    fh_py_He->Write();
    fh_pz_He->Write();
    fh_p_He->Write();
    fh_px_C->Write();
    fh_py_C->Write();
    fh_pz_C->Write();
    fh_p_C->Write();
    fh_px_O->Write();
    fh_py_O->Write();
    fh_pz_O->Write();
    fh_p_O->Write();
    fh_px_px->Write();
    fh_py_py->Write();
    fh_pz_pz->Write();
    fh_p_p->Write();
    fh_psum->Write();
    fh_pzsum->Write();
    fh_psum_vs_event->Write();
    fh_pHe_vs_theta26->Write();
    fh_psum_vs_theta26_nc->Write();
    fh_pz_vs_r->Write();
    fh_px_vs_r->Write();
    fh_py_vs_r->Write();

    fh_theta26->Write();
    fh_minv->Write();
    fh_Erel->Write();
    fh_ErelB->Write();
    fh_ErelR->Write();
    fh_ErelL->Write();
    fh_Erel_vs_psum->Write();
    fh_Erel_vs_theta26->Write();
    fh_Erel_vs_phibc->Write();
    fh_phi_bc_cm->Write();
    fh_phibcm_vs_px->Write();
    fh_phibcm_vs_py->Write();
    fh_phibcm_vs_pz->Write();
    fh_phibcm_vs_psum->Write();

    fh_thetax_px->Write();
    fh_thetay_py->Write();
    fh_thetax_py->Write();
    fh_thetay_px->Write();

    fh_theta26_cm->Write();
    fh_phi26_cm->Write();
    fh_theta_16O->Write();
    fh_phi_16O->Write();
    fh_theta_4He_cm->Write();
    fh_phi_4He_cm->Write();
    fh_theta_bc_cm->Write();
    fh_phi_bc_cm_polar->Write();
    fh_theta_12C_cm->Write();
    fh_phi_12C_cm->Write();

    fh_mass_nc->Write();
    fh_mass->Write();
    fh_energy_nc->Write();
    fh_energy->Write();

    if (fMCTrack)
    {
        fh_target_xy_MC->Write();
        fh_dx->Write();
        fh_dy->Write();
        fh_dz->Write();
        fh_dpxO->Write();
        fh_dpyO->Write();
        fh_dpzO->Write();
        fh_dpO->Write();
        fh_dpxC->Write();
        fh_dpyC->Write();
        fh_dpzC->Write();
        fh_dpC->Write();
        fh_dpxHe->Write();
        fh_dpyHe->Write();
        fh_dpzHe->Write();
        fh_dpHe->Write();
        fh_dp->Write();
        fh_thetax_dpx_C->Write();
        fh_thetay_dpy_C->Write();
        fh_thetax_dpx_He->Write();
        fh_thetay_dpy_He->Write();
        fh_thetax_dpy->Write();
        fh_thetay_dpx->Write();
        fh_x_dpx->Write();
        fh_y_dpy->Write();
        fh_dpy_dpx->Write();
        fh_Erel_vs_thetaMC->Write();
        fh_dErel->Write();
        fh_Erel_dErel->Write();
        fh_dtheta->Write();
        fh_psum_MC->Write();
        fh_pzsum_MC->Write();
        fh_theta26_simu->Write();
        fh_Erel_simu->Write();
        fh_phiMC_bc_cm->Write();
        fh_phiMC_bc_cm_polar->Write();
        fh_Erel_vs_phibcMC->Write();
        fh_phibcm_vs_px_MC->Write();
        fh_phibcm_vs_py_MC->Write();
        fh_phibcm_vs_pz_MC->Write();
        fh_phibcm_vs_psum_MC->Write();
        fh_px_px_mc->Write();
        fh_py_py_mc->Write();
        fh_pz_pz_mc->Write();
        fh_minv_simu->Write();
    }

    if (fMappedItemsCalifa)
    {
        fh_califa_hitenergy->Write();
        fh_califa_calenergy->Write();
        fh_califa_energy->Write();
        fh_califa_energy_nc->Write();
    }
}
ClassImp(R3BGlobalAnalysisS494)
