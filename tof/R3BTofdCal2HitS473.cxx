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
// -----                 R3BTofdCal2Hit                   -----
// -----            Created May 2016 by M.Heil            -----
// -----           Modified Dec 2019 by L.Bott            -----
// ------------------------------------------------------------

#include "R3BTofdCal2HitS473.h"
#include "R3BEventHeader.h"
#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BTCalEngine.h"
#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdHitModulePar.h"
#include "R3BTofdHitPar.h"
#include "R3BTofdGlobalParS473.h"

#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "FairRootManager.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"

#include "TClonesArray.h"
#include "TMath.h"

#include <iostream>

using namespace std;
#define IS_NAN(x) TMath::IsNaN(x)

#define N_TOFD_HIT_PLANE_MAX 4
#define N_TOFD_HIT_PADDLE_MAX 44

namespace
{
	//initPars
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 20; // nanoseconds.
	Long_t maxN = 0;
	Long_t counter = 0;
	Int_t chargeCounter[4] = {0,0,0,0};
	Int_t fact0Counter[4] = {0,0,0,0};
	Int_t negSqrtCounter = 0;
	Int_t losZCounter[6] = {0,0,0,0,0,0,};
} // namespace

R3BTofdCal2HitS473::R3BTofdCal2HitS473()
    : FairTask("TofdCal2HitS473", 1)
    , fCalItems(NULL)
    , fCalItemsLos(NULL)
    , fHitItemsLos(NULL)
    , fHitItems(new TClonesArray("R3BTofdHitData"))
    , fNofHitItems(0)
    , fNofHitPars(0)
    , fHitPar(NULL)
	, fGlobPar(NULL)
    , fTrigger(-1)
    , fTpat(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fTofdQ(1)
    , fTofdHisto(true)
    , fTofdTotPos(true)
    , fnEvents(0)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , maxevent(0)
    , countloshit(0)
    , wrongtrigger(0)
    , wrongtpat(0)
    , headertpat(0)
    , events_in_cal_level(0)
    , inbarcoincidence(0)
    , countreset(0)
    , hitsbeforereset(0)
    , eventstore(0)
    , incoincidence(0)
    , singlehit(0)
    , multihit(0)
    , bars_with_multihit(0)
    , events_wo_tofd_hits(0)
	, event(NULL)
	, fGlobal(true)
	, fCutLosMult(false)
{
    fhLosXYP = NULL;
    fhChargeLosTofD = NULL;
    fh_los_pos = NULL;
    if (fTofdHisto)
    {
        fhxy12 = NULL;
        fhxy12tot = NULL;
        fhxy34 = NULL;
        fhxy34tot = NULL;
        fhCharge = NULL;
        fhAverageCharge = NULL;
		fhChargevsTof = NULL;
		fhChargevsPos = NULL;
		fhChargeLos = NULL;

        for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
        {
            fhQ[i] = NULL;
            fhxy[i] = NULL;
            fhQvsEvent[i] = NULL;
            fhQM[i] = NULL;
            fhMvsQ[i] = NULL;
            fhTdiff[i] = NULL;
			fhTsync[i] = NULL;
			fhEloss[i] = NULL;
			fhChargePlane[i] = NULL;
			fhElossVsBar[i] = NULL;
			fhPosToT[i] = NULL;
			fhPosDiff[i] = NULL;

            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {
                fhQvsPos[i][j] = NULL;
                fhTdiffvsQ[i][2 * j] = NULL;
                fhTdiffvsQ[i][2 * j + 1] = NULL;
                fhQvsQ[i][2 * j] = NULL;
                fhQvsQ[i][2 * j + 1] = NULL;
                fhQvsTof[i][j] = NULL;
                fhTvsTof[i][j] = NULL;
                fhToTvsTofw[i][j] = NULL;
				fhSqrtQvsPosToT[i][j] = NULL;	
				fhLuminosity[i][j] = NULL;
				fhSqrtToT[i][j] = NULL;
            }
        }
    }
}

R3BTofdCal2HitS473::R3BTofdCal2HitS473(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems(NULL)
    , fCalItemsLos(NULL)
    , fHitItemsLos(NULL)
    , fHitItems(new TClonesArray("R3BTofdHitData"))
    , fNofHitItems(0)
    , fNofHitPars(0)
    , fHitPar(NULL)
	, fGlobPar(NULL)
    , fTrigger(-1)
    , fTpat(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fTofdQ(1)
    , fTofdHisto(true)
    , fTofdTotPos(true)
    , fnEvents(0)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , maxevent(0)
    , countloshit(0)
    , wrongtrigger(0)
    , wrongtpat(0)
    , headertpat(0)
    , events_in_cal_level(0)
    , inbarcoincidence(0)
    , countreset(0)
    , hitsbeforereset(0)
    , eventstore(0)
    , incoincidence(0)
    , singlehit(0)
    , multihit(0)
    , bars_with_multihit(0)
    , events_wo_tofd_hits(0)
	, event(NULL)
	, fGlobal(true)
	, fCutLosMult(false)
{
    fhLosXYP = NULL;
    fhChargeLosTofD = NULL;
    fh_los_pos = NULL;
    if (fTofdHisto)
    {
        fhxy12 = NULL;
        fhxy12tot = NULL;
        fhxy34 = NULL;
        fhxy34tot = NULL;
        fhCharge = NULL;
        fhAverageCharge = NULL;
		fhChargevsTof = NULL;
		fhChargevsPos = NULL;
		fhChargeLos = NULL;

        for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
        {
            fhQ[i] = NULL;
            fhxy[i] = NULL;
            fhQvsEvent[i] = NULL;
            fhQM[i] = NULL;
            fhMvsQ[i] = NULL;
            fhTdiff[i] = NULL;
			fhTsync[i] = NULL;
			fhEloss[i] = NULL;
			fhChargePlane[i] = NULL;
			fhElossVsBar[i] = NULL;
			fhPosToT[i] = NULL;
			fhPosDiff[i] = NULL;

            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {
                fhQvsPos[i][j] = NULL;
                fhTdiffvsQ[i][2 * j] = NULL;
                fhTdiffvsQ[i][2 * j + 1] = NULL;
                fhQvsQ[i][2 * j] = NULL;
                fhQvsQ[i][2 * j + 1] = NULL;
                fhQvsTof[i][j] = NULL;
                fhTvsTof[i][j] = NULL;
                fhToTvsTofw[i][j] = NULL;
				fhSqrtQvsPosToT[i][j] = NULL;
				fhLuminosity[i][j] = NULL;
				fhSqrtToT[i][j] = NULL;
            }
        }
    }
}

R3BTofdCal2HitS473::~R3BTofdCal2HitS473()
{
    if (fhLosXYP)
        delete fhLosXYP;
    if (fhChargeLosTofD)
        delete fhChargeLosTofD;
    if (fh_los_pos)
        delete fh_los_pos;
    if (fTofdHisto)
    {
        if (fhxy12)
            delete fhxy12;
        if (fhxy12tot)
            delete fhxy12tot;
        if (fhxy34)
            delete fhxy34;
        if (fhxy34tot)
            delete fhxy34tot;
        if (fhAverageCharge)
            delete fhAverageCharge;
        if (fhChargevsTof)
			delete fhChargevsTof;
        if (fhChargevsPos)
			delete fhChargevsPos;
        if (fhCharge)
            delete fhCharge;
		if (fhChargeLos)
			delete fhChargeLos;
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            if (fhQ[i])
                delete fhQ[i];
            if (fhxy[i])
                delete fhxy[i];
            if (fhQvsEvent[i])
                delete fhQvsEvent[i];
            if (fhQM[i])
                delete fhQM[i];
            if (fhMvsQ[i])
                delete fhMvsQ[i];
            if (fhTdiff[i])
                delete fhTdiff[i];
			if (fhTsync[i])
				delete fhTsync[i];
			if (fhEloss[i])
				delete fhEloss[i];
			if (fhChargePlane[i])
				delete fhChargePlane[i];
			if (fhElossVsBar[i])
				delete fhElossVsBar[i];
			if(fhPosToT[i])
				delete fhPosToT[i];
			if(fhPosDiff[i])
				delete fhPosDiff[i];
            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {
                if (fhQvsPos[i][j])
                    delete fhQvsPos[i][j];
                if (fhTdiffvsQ[i][2 * j])
                    delete fhTdiffvsQ[i][2 * j];
                if (fhTdiffvsQ[i][2 * j + 1])
                    delete fhTdiffvsQ[i][2 * j + 1];
                if (fhQvsQ[i][2 * j])
                    delete fhQvsQ[i][2 * j];
                if (fhQvsQ[i][2 * j + 1])
                    delete fhQvsQ[i][2 * j + 1];
                if (fhQvsTof[i][j])
                    delete fhQvsTof[i][j];
                if (fhTvsTof[i][j])
                    delete fhTvsTof[i][j];
                if (fhToTvsTofw[i][j])
                    delete fhQvsTof[i][j];
				if (fhSqrtQvsPosToT[i][j])
					delete fhSqrtQvsPosToT[i][j];
				if(fhLuminosity[i][j])
					delete fhLuminosity[i][j];
				if(fhSqrtToT[i][j])
					delete fhSqrtToT[i][j];
            }
        }
    }
    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
}

InitStatus R3BTofdCal2HitS473::Init()
{
    fHitPar = (R3BTofdHitPar*)FairRuntimeDb::instance()->getContainer("TofdHitPar");
    if (!fHitPar)
    {
        LOG(error) << "Could not get access to TofdHitPar-Container.";
        fNofHitPars = 0;
        return kFATAL;
    }
    fNofHitPars = fHitPar->GetNumModulePar();
    if (fNofHitPars == 0)
    {
        LOG(error) << "There are no Hit parameters in container TofdHitPar";
        return kFATAL;
    }
	fGlobPar = (R3BTofdGlobalParS473*)FairRuntimeDb::instance()->getContainer("TofdGlobalPar");
	if(!fGlobPar)
	{
		LOG(error) << "Could not get access to TofdGlobalHitPar-Container.";
		return kFATAL;
	}
    // get access to Cal data
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
    fCalItems = (TClonesArray*)mgr->GetObject("TofdCal");
    if (NULL == fCalItems)
        LOG(fatal) << "Branch TofdCal not found";
    maxevent = mgr->CheckMaxEventNo();
    fCalItemsLos = (TClonesArray*)mgr->GetObject("LosCal");
    if (NULL == fCalItemsLos)
        LOG(warning) << "Branch LosCal not found";
    fHitItemsLos = (TClonesArray*)mgr->GetObject("LosHit");
    if (NULL == fHitItemsLos)
        LOG(warning) << "Branch LosHit not found";
    // request storage of Hit data in output tree
    mgr->Register("TofdHit", "Land", fHitItems, kTRUE);

	maxN = maxevent;
    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BTofdCal2HitS473::SetParContainers()
{
    fHitPar = (R3BTofdHitPar*)FairRuntimeDb::instance()->getContainer("TofdHitPar");
    if (!fHitPar)
    {
        LOG(error) << "Could not get access to TofdHitPar-Container.";
        fNofHitPars = 0;
        return;
    }
	fGlobPar = (R3BTofdGlobalParS473*)FairRuntimeDb::instance()->getContainer("TofdGlobalPar");
	if(!fGlobPar)
	{
		LOG(error) << "Could not get access to TofdGlobalPar-Container.";
		return;
	}
}

InitStatus R3BTofdCal2HitS473::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BTofdCal2HitS473::Exec(Option_t* option)
{
	counter++;
	event = new std::vector<hit>();
    if (fnEvents / 100000. == (int)fnEvents / 100000)
        std::cout << "\rEvents: " << fnEvents << " / " << maxevent << " (" << (int)(fnEvents * 100. / maxevent)
                  << " %) " << std::flush;

    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
    {
        wrongtrigger++;
        return;
    }
    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit = fTpat - 1;
    if (fTpat_bit >= 0)
    {
        Int_t itpat = header->GetTpat();
        Int_t tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        if ((header) && (tpatvalue == 0))
        {
            wrongtpat++;
            return;
        }
    }
    headertpat++;
    Double_t timeRef = 0.;
    Double_t timeLos = 0;
    Double_t timeP0 = 0.;
    Double_t LosTresM = 0;
    Double_t LosQ = 0;
    Double_t xLosP = 1000;
    Double_t yLosP = 1000;
    Double_t randx;
    std::vector<std::vector<std::vector<Double_t>>> q;
    std::vector<std::vector<std::vector<Double_t>>> tof;
    std::vector<std::vector<std::vector<Double_t>>> x;
    std::vector<std::vector<std::vector<Double_t>>> y;
    std::vector<std::vector<std::vector<Double_t>>> yToT;
    UInt_t vmultihits[N_PLANE_MAX + 1][N_TOFD_HIT_PADDLE_MAX * 2 + 1];
    for (Int_t i = 0; i <= fNofPlanes; i++)
    {
        q.push_back(std::vector<std::vector<Double_t>>());
        tof.push_back(std::vector<std::vector<Double_t>>());
        x.push_back(std::vector<std::vector<Double_t>>());
        y.push_back(std::vector<std::vector<Double_t>>());
        yToT.push_back(std::vector<std::vector<Double_t>>());
        for (Int_t j = 0; j <= 2 * N_TOFD_HIT_PADDLE_MAX; j++)
        {
            vmultihits[i][j] = 0;
            q[i].push_back(std::vector<Double_t>());
            tof[i].push_back(std::vector<Double_t>());
            x[i].push_back(std::vector<Double_t>());
            y[i].push_back(std::vector<Double_t>());
            yToT[i].push_back(std::vector<Double_t>());
        }
    }
    if (fHitItemsLos)
    {
        Int_t nHits = fHitItemsLos->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            countloshit++;
            LOG(warning) << "LOS Ihit  " << ihit << " " << nHits;
            R3BLosHitData* hitData = (R3BLosHitData*)fHitItemsLos->At(ihit);
            if (ihit == 0)
                timeLos = hitData->fTime_ns;
            LOG(warning) << "LOS Time " << timeLos;
            if (std::isnan(timeLos))
                return; /// in s444 run2930 Event 22333208
            if (NULL == fh_los_pos)
            {
                char strName[255];
                sprintf(strName, "LOS_X_vs_Y_MCFD");
                //     fh_los_pos = new TH2F(strName, "", 2000, -10., 10., 2000, -10., 10.);
                fh_los_pos = new TH2F(strName, "", 2000, -10., 10., 2000, -10., 10.);
            }
            fh_los_pos->Fill(hitData->fX_cm, hitData->fY_cm);

            if (ihit == 0)
                LosQ = hitData->fZ;
			if(LosQ < 0)
			{
				if(LosQ < -10)
				{ 
					losZCounter[0]++;
				} else 
				{
					losZCounter[1]++;
				}
			} else
			{
				if(LosQ == 0.)
				{
					losZCounter[2]++;
				} else
				{
					if(LosQ < 1)
					{
						losZCounter[3]++;
					} else
					{
						if(LosQ < 10)
						{
							losZCounter[4]++;
						} else
						{
							losZCounter[5]++;
						}
					}
				}
			}

            if (NULL == fhChargeLosTofD)
            {
                char strName[255];
                sprintf(strName, "LosQ_vs_TofdQ");
                fhChargeLosTofD = new TH2F(strName, "", 400, 0., 100., 800, 0., 200.);
                fhChargeLosTofD->GetYaxis()->SetTitle("Charge LOS");
                fhChargeLosTofD->GetXaxis()->SetTitle("Charge ToFD");
            }
			//Cut on Los Multiplicity
			if (fCutLosMult)
			{
				if (fhChargeLos)
				{
					//TODO: Check for pileup
				}
				if (nHits > 1) return;
			}
			if (NULL == fhChargeLos)	
			{
				fhChargeLos = new TH1F("LosQ","LosQ", 440,-10.,100.);
				fhChargeLos->GetXaxis()->SetTitle("Charge LOS");
				fhChargeLos->GetYaxis()->SetTitle("Entries");
			}
			fhChargeLos->Fill(LosQ);
        }
    }

    // std::cout<<"new event!*************************************\n";

    Int_t nHits = fCalItems->GetEntries();
    Int_t nHitsEvent = 0;
    // Organize cals into bars.
    struct Entry
    {
        std::vector<R3BTofdCalData*> top;
        std::vector<R3BTofdCalData*> bot;
    };
    std::map<size_t, Entry> bar_map;
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        auto* hit0 = (R3BTofdCalData*)fCalItems->At(ihit);
        size_t idx = hit0->GetDetectorId() * fPaddlesPerPlane * hit0->GetBarId();
        // std::cout << "Hits: " << hit->GetDetectorId() << ' ' << hit->GetBarId() << ' ' << hit->GetSideId() << ' '
        //          << hit->GetTimeLeading_ns() << ' ' << hit->GetTimeTrailing_ns() << '\n';
        auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
        auto& vec = 1 == hit0->GetSideId() ? ret.first->second.top : ret.first->second.bot;
        vec.push_back(hit0);
        events_in_cal_level++;
    }

    // Find coincident PMT hits.
    // std::cout << "Print:\n";
    for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
    {
    reset:
        // for (auto it2 = it->second.top.begin(); it->second.top.end() != it2; ++it2) {
        // std::cout << "Top: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        // }
        // for (auto it2 = it->second.bot.begin(); it->second.bot.end() != it2; ++it2) {
        // std::cout << "Bot: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        // }
        auto const& top_vec = it->second.top;
        auto const& bot_vec = it->second.bot;
        size_t top_i = 0;
        size_t bot_i = 0;
        for (; top_i < top_vec.size() && bot_i < bot_vec.size();)
        {
            auto top = top_vec.at(top_i);
            auto bot = bot_vec.at(bot_i);
            auto top_ns = top->GetTimeLeading_ns();
            auto bot_ns = bot->GetTimeLeading_ns();
            auto dt = top_ns - bot_ns;
            // Handle wrap-around.
            auto dt_mod = fmod(dt + c_range_ns, c_range_ns);
            if (dt < 0)
            {
                // We're only interested in the short time-differences, so we
                // want to move the upper part of the coarse counter range close
                // to the lower range, i.e. we cut the middle of the range and
                // glue zero and the largest values together.
                dt_mod -= c_range_ns;
            }
            // std::cout << top_i << ' ' << bot_i << ": " << top_ns << ' ' << bot_ns << " = " << dt << ' ' <<
            // std::abs(dt_mod) << '\n';
            if (std::abs(dt_mod) < c_bar_coincidence_ns)
            {
                inbarcoincidence++;
                // Hit!
                // std::cout << "Hit!\n";
                Int_t iPlane = top->GetDetectorId(); // 1..n
                Int_t iBar = top->GetBarId();        // 1..n
                if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                {
                    LOG(error) << "R3BTofdCal2HitPar::Exec() : more detectors than expected! Det: " << iPlane
                               << " allowed are 1.." << fNofPlanes;
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    LOG(error) << "R3BTofdCal2HitPar::Exec() : more bars then expected! Det: " << iBar
                               << " allowed are 1.." << fPaddlesPerPlane;
                    continue;
                }

                auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);

				if(top_tot == 0.){cout << "top_tot = 0" << endl;}
				if(bot_tot == 0.){cout << "bot_tot = 0" << endl;}

                // register multi hits
                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
                if (iPlane == 1 || iPlane == 3)
                    vmultihits[iPlane][iBar * 2 - 2] += 1;
                if (iPlane == 2 || iPlane == 4)
                    vmultihits[iPlane][iBar * 2] += 1;
                vmultihits[iPlane][iBar * 2 - 1] += 1;

                nHitsEvent += 1;
                R3BTofdHitModulePar* par = fHitPar->GetModuleParAt(iPlane, iBar);
				R3BTofdGlobalParS473* gPar = fGlobPar;
                if (!par)
                {
                    LOG(info) << "R3BTofdCal2Hit::Exec : Hit par not found, Plane: " << top->GetDetectorId()
                              << ", Bar: " << top->GetBarId();
                    continue;
                }
                // walk corrections
                if (par->GetPar1Walk() == 0. || par->GetPar2Walk() == 0. || par->GetPar3Walk() == 0. ||
                    par->GetPar4Walk() == 0. || par->GetPar5Walk() == 0.)
                    LOG(info) << "Walk correction not found!";
                bot_ns = bot_ns - walk(bot_tot,
                                                 par->GetPar1Walk(),
                                                 par->GetPar2Walk(),
                                                 par->GetPar3Walk(),
                                                 par->GetPar4Walk(),
                                                 par->GetPar5Walk());
                top_ns = top_ns - walk(top_tot,
                                                 par->GetPar1Walk(),
                                                 par->GetPar2Walk(),
                                                 par->GetPar3Walk(),
                                                 par->GetPar4Walk(),
                                                 par->GetPar5Walk());
                
		// calculate tdiff
                auto tdiff = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2()));

                // calculate time-of-flight
				// calcToF
                if (timeLos == 0)
                    LOG(warning) << "Los Time is zero! ";
                Double_t ToF0 = (bot_ns + top_ns) / 2. - timeLos - par->GetSync();
				Double_t ToF = ToF0 - gPar->GetToFOffset(iPlane);
                if (std::isnan(ToF))
                {
		    cout << "bot " << bot_ns << ", top " << top_ns << ", timeLos " << timeLos << ", sync" << par->GetSync() << endl;
                    break;
		    LOG(fatal) << "ToFD ToF not found";

                }
                if (timeP0 == 0.)
                    timeP0 = ToF;

                if (timeLos == 0)
                { // no LOS in s454
                    /// What to do here:
                    /// check if all ToF in one event are in a range of 3000ns (readout window) and shift times
                    /// according to that
                    ///
                    ///       this is the first hit
                    ///       I
                    /// e.g. 171; 9439; 179; 1117; 175 -->> 171+c_range_ns; 9439; 179+c_range_ns; 1117+c_range_ns;
                    /// 175+c_range_ns
                    ///             I
                    ///             this should be the first hit -> counter resets -> other hits follow
                    if (ToF - timeP0 < -3000.)
                    {
                        ToF += c_range_ns;
                    }
                    if (ToF - timeP0 > 3000.)
                    {
                        timeP0 = ToF;
                        it = bar_map.begin();
                        countreset++;
                        hitsbeforereset += nHitsEvent;
                        for (Int_t i = 0; i <= fNofPlanes; i++)
                        {
                            for (Int_t j = 0; j <= 2 * N_TOFD_HIT_PADDLE_MAX; j++)
                            {
                                tof[i][j].clear();
                                x[i][j].clear();
                                y[i][j].clear();
                                yToT[i][j].clear();
                                q[i][j].clear();
                                vmultihits[i][j] = 0;
                                nHitsEvent = 0;
                            }
                        }
                        LOG(warning) << "Found new first hit -> will reset";
                        goto reset; /// TODO: how to do without goto?
                    }
                }

                if (timeLos != 0)
                {
                    LOG(debug) << "Found LOS detector";
                    LOG(debug) << "check for coincidence range: TOF: " << ToF << " c_range: " << c_range_ns << "\n";
                    while (ToF < -c_range_ns / 2)
                    {
                        ToF += c_range_ns;
                        LOG(debug) << "Shift up\n";
                    }
                    while (ToF > c_range_ns / 2)
                    {
                        ToF -= c_range_ns;
                        LOG(debug) << "Shift down\n";
                    }
                }

				hit h;

				h.p = iPlane;


                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
				/*
                if (iPlane == 1 || iPlane == 3)
                    tof[iPlane][iBar * 2 - 2].push_back(ToF);
                if (iPlane == 2 || iPlane == 4)
                    tof[iPlane][iBar * 2].push_back(ToF);
                tof[iPlane][iBar * 2 - 1].push_back(ToF);
				*/
				h.t = ToF;


                // calculate y-position
                auto pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();

                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
				/*
                if (iPlane == 1 || iPlane == 3)
                    y[iPlane][iBar * 2 - 2].push_back(pos);
                if (iPlane == 2 || iPlane == 4)
                    y[iPlane][iBar * 2].push_back(pos);
                y[iPlane][iBar * 2 - 1].push_back(pos);
				*/
				h.y = pos;

                // calculate y-position from ToT
                auto posToT =
                    par->GetLambda() * log((top_tot * par->GetToTOffset2()) / (bot_tot * par->GetToTOffset1()));
                /*
				if (iPlane == 1 || iPlane == 3)
                    yToT[iPlane][iBar * 2 - 2].push_back(posToT);
                if (iPlane == 2 || iPlane == 4)
                    yToT[iPlane][iBar * 2].push_back(posToT);
                yToT[iPlane][iBar * 2 - 1].push_back(posToT);
				*/

                if (fTofdTotPos){
                    pos = posToT;
					h.y = posToT;
				}

                // calculate x-position
				Double_t paddle_width = 2.7;
				Double_t paddle_thickness = 0.5;
				Double_t air_gap_paddles = 0.04;
				Double_t air_gap_layer = 4.5;
				Int_t number_paddles = 44;
				Int_t number_layers = 4;

				Double_t detector_width = number_paddles * paddle_width + (number_paddles - 1) * air_gap_paddles + paddle_width / 2;
				Double_t detector_thickness = (number_layers - 1) * air_gap_layer + number_layers * paddle_thickness;
                randx = (std::rand() / (float)RAND_MAX);
                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
				/*
                if (iPlane == 1 || iPlane == 3)
                    x[iPlane][iBar * 2 - 2].push_back(iBar * 2.8 - 21. * 2.8 - 1.4 - 1.4 * randx);
                if (iPlane == 2 || iPlane == 4)
                    x[iPlane][iBar * 2].push_back(iBar * 2.8 - 21. * 2.8 + 1.4 - 1.4 * randx);
                x[iPlane][iBar * 2 - 1].push_back(iBar * 2.8 - 21. * 2.8 - 1.4 * randx);
				*/

				Double_t xp = -1000.;
				Double_t xp2 = -1000.;
				if(iPlane == 1 || iPlane == 3)
				{
					xp = -detector_width / 2 + (iBar - 1) * (paddle_width + air_gap_paddles) + paddle_width * randx;
					h.vb1 = iBar * 2 - 2;
				}
				if(iPlane == 2 || iPlane == 4)
				{
					xp = -detector_width / 2 + paddle_width / 2  + (iBar - 1) * (paddle_width + air_gap_paddles) + paddle_width * randx;
				    h.vb1 = iBar * 2 - 1;
				}
				h.x1 = xp;
				h.vb2 = h.vb1 + 1;
				h.x2 = xp;			



                // correct for position dependence and calculate nuclear charge Z
                Double_t para[4];
                para[0] = par->GetPar1a();
                para[1] = par->GetPar1b();
                para[2] = par->GetPar1c();
                para[3] = par->GetPar1d();

				//calcQ
                Double_t qb = 0.;
				Double_t q0 = 0.;
				Double_t tot = 0.;
                if (fTofdTotPos)
                {

                	para[0] = par->GetPola();
                	para[1] = par->GetPolb();
                	para[2] = par->GetPolc();
                	para[3] = par->GetPold();
                    // via pol3
					Double_t factor =  (para[0] + para[1] * posToT + para[2] * pow(posToT, 2) + para[3] * pow(posToT, 3));
					if(factor == 0) fact0Counter[iPlane - 1]++;
					tot = TMath::Sqrt(top_tot * bot_tot);
                    q0 = tot / factor;
                    qb = q0 * gPar->GetEloss(iPlane);
                }
                else
                {
                    // via double exponential:
                    auto q1 =
                        bot_tot / (para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3]);
                    para[0] = par->GetPar2a();
                    para[1] = par->GetPar2b();
                    para[2] = par->GetPar2c();
                    para[3] = par->GetPar2d();
                    auto q2 =
                        top_tot / (para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3]);
                    q1 = q1 * gPar->GetEloss(iPlane);
                    q2 = q2 * gPar->GetEloss(iPlane);
                    q0 = (q1 + q2) / 2.;
					qb = q0 * gPar->GetEloss(iPlane);
                }
				
				h.eloss = qb;

				event->push_back(h);

				nHitsEvent++;

		//Apply velocity correction

				
                Double_t parz[3];
                parz[0] = par->GetPar1za();
                parz[1] = par->GetPar1zb();
                parz[2] = par->GetPar1zc();
				


                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
				/*
                if (parz[0] > 0 && parz[2] > 0)
                {
                    if (iPlane == 1 || iPlane == 3)
                        q[iPlane][iBar * 2 - 2].push_back(parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                    if (iPlane == 2 || iPlane == 4)
                        q[iPlane][iBar * 2].push_back(parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                    q[iPlane][iBar * 2 - 1].push_back(parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                }
                else
                {
                    if (iPlane == 1 || iPlane == 3)
                        q[iPlane][iBar * 2 - 2].push_back(qb);
                    if (iPlane == 2 || iPlane == 4)
                        q[iPlane][iBar * 2].push_back(qb);
                    q[iPlane][iBar * 2 - 1].push_back(qb);
                    parz[0] = 1.;
                    parz[1] = 0.;
                    parz[2] = 1.;
                }
				*/


                if (fTofdHisto)
                {
                    // fill control histograms
                    CreateHistograms(iPlane, iBar);
                    // fhTof[iPlane-1]->Fill(iBar,ToF);
                    fhTsync[iPlane-1]->Fill(iBar,ToF);
                    fhTdiff[iPlane - 1]->Fill(iBar, tdiff);
                    fhQvsPos[iPlane - 1][iBar - 1]->Fill(pos, q0);
                    fhQvsTof[iPlane - 1][iBar - 1]->Fill(qb, ToF);
                    fhTvsTof[iPlane - 1][iBar - 1]->Fill(dt_mod, ToF);
                    fhToTvsTofw[iPlane - 1][iBar - 1]->Fill((bot_tot + top_tot) / 2.,
                	                                            ToF); // needed to get TOF w/o walk correction
					fhEloss[iPlane - 1] -> Fill(qb);
					fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->Fill(posToT, qb);
					fhElossVsBar[iPlane - 1]->Fill(iBar,qb);
					fhLuminosity[iPlane - 1][iBar - 1]->Fill(q0);
					fhSqrtToT[iPlane - 1][iBar - 1]->Fill(tot);
					fhPosToT[iPlane - 1]->Fill(iBar,posToT);
					fhPosDiff[iPlane - 1]->Fill(iBar,pos);
                }


                // Time reference in case on has the master signal in one of the TDC channels.
                // Not used at the moment.
                timeRef = 0;

                ++top_i;
                ++bot_i;
            }
            else if (dt < 0 && dt > -c_range_ns / 2)
            {
                ++top_i;
                LOG(warning) << "Not in bar coincidence increase top counter";
            }
            else
            {
                ++bot_i;
                LOG(warning) << "Not in bar coincidence increase bot counter";
            }
        }
    }

    // Now all hits in this event are analyzed

    Double_t hit_coinc = 500.; // coincidence window for hits in one event in ns. physics says max 250 ps
    Double_t maxChargeDiff = 1.;   // maximum charge difference between two planes for averaged hits

    LOG(debug) << "Hits in this event: " << nHitsEvent;
    if (nHitsEvent == 0) events_wo_tofd_hits++;

	//time order of hits in this event
	//std::sort(event->begin(), event->end(), [](hit const& a, hit const& b) {return a.t < b.t; });

    // init arrays to store hits
	/*
    Double_t tArrQ[2 * nHitsEvent + 1];
    Double_t tArrT[2 * nHitsEvent + 1];
    Double_t tArrX[2 * nHitsEvent + 1];
    Double_t tArrY[2 * nHitsEvent + 1];
    Double_t tArrYT[2 * nHitsEvent + 1];
    Double_t tArrP[2 * nHitsEvent + 1];
    Double_t tArrB[2 * nHitsEvent + 1];
    Bool_t tArrU[2 * nHitsEvent + 1];
    for (int i = 0; i < (2 * nHitsEvent + 1); i++)
    {
        tArrQ[i] = -1.;
        tArrT[i] = -1.;
        tArrX[i] = -1.;
        tArrY[i] = -1.;
        tArrYT[i] = -1.;
        tArrP[i] = -1.;
        tArrB[i] = -1.;
        tArrU[i] = kFALSE;
    }

    for (Int_t i = 1; i <= fNofPlanes; i++)
    {
        for (Int_t j = 0; j < fPaddlesPerPlane * 2 + 1; j++)
        {
            if (vmultihits[i][j] > 1)
            {
                bars_with_multihit++;
                multihit += vmultihits[i][j] - 1;
            }
        }
    }

    // order events for time
    for (Int_t i = 1; i <= fNofPlanes; i++)
    { // loop over planes i
        for (Int_t j = 0; j < fPaddlesPerPlane * 2 + 1; j++)
        { // loop over virtual paddles j
            if (tof[i][j].empty() == false)
            { // check paddle for entries
                for (Int_t m = 0; m < tof[i][j].size(); m++)
                { // loop over multihits m
                    Int_t p = 0;
                    if (tArrT[0] == -1.)
                    { // first entry
                        LOG(debug) << "First entry plane/bar " << i << "/" << j;
                        tArrQ[0] = q[i][j].at(m);
                        tArrT[0] = tof[i][j].at(m);
                        tArrX[0] = x[i][j].at(m);
                        tArrY[0] = y[i][j].at(m);
                        tArrYT[0] = yToT[i][j].at(m);
                        tArrP[0] = i;
                        tArrB[0] = j;
                    }
                    else
                    {
                        if (tof[i][j].at(m) < tArrT[0])
                        { // new first entry with smaller time
                            LOG(debug) << "Insert new first " << i << " " << j;
                            insertX(2 * nHitsEvent, tArrQ, q[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrT, tof[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrX, x[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrY, y[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrYT, yToT[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrP, i, 1);
                            insertX(2 * nHitsEvent, tArrB, j, 1);
                        }
                        else
                        {
                            while (tof[i][j].at(m) > tArrT[p] && tArrT[p] != -1.)
                            {
                                p++; // find insert position
                                if (p > 2 * nHitsEvent + 1)
                                    LOG(fatal) << "Insert position oor"; // should not happen
                            }

                            LOG(debug) << "Will insert at " << p;
                            if (p > 0 && tof[i][j].at(m) > tArrT[p - 1] && tof[i][j].at(m) != tArrT[p])
                            { // insert at right position
                                LOG(debug) << "Insert at " << p << " " << i << " " << j;
                                insertX(2 * nHitsEvent, tArrQ, q[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrT, tof[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrX, x[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrY, y[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrYT, yToT[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrP, i, p + 1);
                                insertX(2 * nHitsEvent, tArrB, j, p + 1);
                            }
                            else
                            {
                                if (tof[i][j].at(m) == tArrT[p])
                                { // handle virtual bars
                                    LOG(debug) << "Insert virtual bar " << i << " " << j;
                                    insertX(2 * nHitsEvent, tArrQ, q[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrT, tof[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrX, x[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrY, y[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrYT, yToT[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrP, i, p + 2);
                                    insertX(2 * nHitsEvent, tArrB, j, p + 2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
	*/

    // Now we can analyze the hits in this event

	/*
    if (fTofdHisto)
    {
        for (Int_t a = 0; a < 2 * nHitsEvent; a++)
        { // loop over all hits
            eventstore++;
            fhQ[((Int_t)tArrP[a]) - 1]->Fill(tArrB[a], tArrQ[a]);        // charge per plane
            fhQvsEvent[((Int_t)tArrP[a]) - 1]->Fill(fnEvents, tArrQ[a]); // charge vs event #
            if (fTofdTotPos)
            {
                fhxy[((Int_t)tArrP[a]) - 1]->Fill(tArrB[a], tArrYT[a]); // xy of plane
            }
            else
            {
                fhxy[((Int_t)tArrP[a]) - 1]->Fill(tArrB[a], tArrY[a]); // xy of plane
            }
        }
    }
	*/
		//std::vector<hit>* vec = new std::vector<hit>();

		/*
        while (tArrT[ihit] < time0 + hit_coinc)
        { // check if in coincidence window
            incoincidence++;
		
            hitscoinc++; // number of hits in coincidence

            if (fTofdHisto)
            {
                if (tArrP[ihit] == plane0 && charge0 != tArrQ[ihit])
                {
                    fhQM[(Int_t)tArrP[ihit] - 1]->Fill(charge0, tArrQ[ihit]);
                }
            }

            LOG(debug) << "Hit in coincidence window: " << tArrP[ihit] << " " << tArrB[ihit] << " " << tArrT[ihit]
                       << " " << tArrQ[ihit];

			//vec->push_back(hit());
			hit h;
			h.t = tArrT[ihit];
			h.p = tArrP[ihit];
			h.vb = tArrB[ihit];
			h.eloss = tArrQ[ihit];
			h.x = tArrX[ihit];
			h.y = tArrY[ihit];
			if(fTofdTotPos) h.y = tArrYT[ihit];

			fhEloss[h.p-1]->Fill(h.eloss);

			event->push_back(h);

            ihit++;
            if (ihit >= 2 * nHitsEvent)
                break;
        }
		*/

		//calcBeta


		Double_t c0 = 299.792458; // c0 in mm/ns

		for(Int_t cHit = 0; cHit < event->size(); cHit++)
		{
			hit* h = &event->at(cHit);
			R3BTofdGlobalParS473* gPar = fGlobPar;
			//TODO: calculate beta from ToFs
			Double_t dist = gPar->GetToFSlope(h->p)*c0;
			Double_t beta = (dist/h->t)/c0;
			//TODO: calculate Q from Eloss and beta
			Double_t factor = gPar->GetBetheSlope() * h->eloss + gPar->GetBetheOffset();
			if((gPar->GetBetheSlope() * h->eloss + gPar->GetBetheOffset()) < 0)
			{ 
				//cout << "sqrt < 0" << endl;
				negSqrtCounter++;
				continue;
			}
			Double_t charge = beta * TMath::Sqrt(gPar->GetBetheSlope() * h->eloss + gPar->GetBetheOffset());
			if((charge < 0 || charge > 80) && h->eloss < 1e7) 
			{	
			//	cout << beta << " " << h->eloss << " " << gPar->GetBetheSlope() << " " << 
			//		" " << gPar->GetBetheOffset() << " " << factor << " " << charge << endl;
				chargeCounter[h->p - 1]++;
			}
			fhChargePlane[h->p - 1]->Fill(charge);
			fhCharge->Fill(charge);
			//TODO: save Hits in file
			new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData(h->t,
															h->x1,
															h->y,
															charge,
															h->t,
															h->eloss,
															h->p,
															h->vb1,
															h->t,
															h->t
												);
			//hit++;
			if(fTofdHisto)
			{
				fhxy[h->p - 1]->Fill(h->x1, h->y);
			}
		}

		



	/*
    for (Int_t hit = 0; hit < 2 * nHitsEvent; hit++)
    { // loop over not averaged hits
        if (tArrU[hit] == false)
        {
            LOG(debug) << "Single Hit for Plane " << tArrP[hit] << " " << tArrB[hit];
            tArrU[hit] = tArrU[hit + 1] = true;
            // store single hits only seen in planes
            singlehit++;
            if (fTofdTotPos)
            {
                new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData(tArrT[hit],
                                                                  (tArrX[hit] + tArrX[hit + 1]) / 2.,
                                                                  tArrYT[hit],
                                                                  tArrQ[hit],
                                                                  -1.,
                                                                  tArrQ[hit],
                                                                  tArrP[hit]);
            }
            else
            {
                new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData(tArrT[hit],
                                                                  (tArrX[hit] + tArrX[hit + 1]) / 2.,
                                                                  tArrY[hit],
                                                                  tArrQ[hit],
                                                                  -1.,
                                                                  tArrQ[hit],
                                                                  tArrP[hit]);
            }
            hit++;
        }
    }
	*/
    // std::cout<<"Used up hits in this event:\n";
    /*
	for (Int_t a = 0; a < 2 * nHitsEvent; a++)
    {
        // std::cout << tArrU[a] << " ";
        if (tArrU[a] != true)
            LOG(fatal);

    }
	*/
    // std::cout << "\n";
    //LOG(debug) << "------------------------------------------------------\n";

    fnEvents++;
}

void R3BTofdCal2HitS473::CreateHistograms(Int_t iPlane, Int_t iBar)
{
    Double_t max_charge = 80.;
    // create histograms if not already existing

    if (NULL == fhTdiff[iPlane - 1])
    {
        char strName1[255];
        char strName2[255];
        sprintf(strName1, "Time_Diff_Plane_%d", iPlane);
        sprintf(strName2, "Time Diff Plane %d", iPlane);
        fhTdiff[iPlane - 1] = new TH2F(strName1, strName2, 50, 0, 50, 400, -8., 8.);
        fhTdiff[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTdiff[iPlane - 1]->GetYaxis()->SetTitle("Time difference (PM1 - PM2) in ns");
    }

	if (NULL == fhTsync[iPlane - 1])
	{
		char strName1[255];
		char strName2[255];
		sprintf(strName1, "Time_Sync_Plane_%d", iPlane);
		sprintf(strName2, "Time Sync Plane %d", iPlane);
		fhTsync[iPlane - 1] = new TH2F(strName1, strName2, 50,0,50,800,30,110);
		fhTsync[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
		fhTsync[iPlane - 1]->GetYaxis()->SetTitle("ToF in ns");
	}

	if(NULL == fhEloss[iPlane - 1])
	{
		char strName1[255];
		char strName2[255];
		sprintf(strName1, "Eloss_Plane_%d", iPlane);
		sprintf(strName2, "Eloss Plane %d; Eloss / MeV;", iPlane);
		fhEloss[iPlane - 1] = new TH1F(strName1, strName2, 8000., 2000.,10000.);
	}

    if (NULL == fhQvsTof[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Q_vs_ToF_Plane_%d_Bar_%d", iPlane, iBar);
        fhQvsTof[iPlane - 1][iBar - 1] = new TH2F(strName, "", 1000, 0., 10000, 1000, -10, 40);
        fhQvsTof[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToF in ns");
        fhQvsTof[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Charge");
    }

    if (NULL == fhTvsTof[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "T_vs_ToF_Plane_%d_Bar_%d", iPlane, iBar);
        fhTvsTof[iPlane - 1][iBar - 1] = new TH2F(strName, "", 625, -25, 25, 1000, -10, 40);
        fhTvsTof[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToF in ns");
        fhTvsTof[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("T1-T2 in ns");
    }

    if (NULL == fhToTvsTofw[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "ToT_vs_ToF_Plane_%d_Bar_%d_w", iPlane, iBar);
        fhToTvsTofw[iPlane - 1][iBar - 1] = new TH2F(strName, "", 1000, 0., 200, 1000, -10, 40);
        fhToTvsTofw[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT in ns");
        fhToTvsTofw[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToF in ns");
    }

    if (NULL == fhQvsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Q_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        fhQvsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 200, -100, 100, 600, 0., 1.2);
        fhQvsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Charge");
        fhQvsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Position in cm");
    }

    if (NULL == fhQ[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "Q_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Q Plane %d ", iPlane);
        fhQ[iPlane - 1] = new TH2F(strName1, strName2, 90, 0, 90, max_charge * 10, 0., max_charge);
        fhQ[iPlane - 1]->GetYaxis()->SetTitle("Charge");
        fhQ[iPlane - 1]->GetXaxis()->SetTitle("Paddle number");
    }
    /*
    if (NULL == fhxy[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "xy_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "xy of Plane %d ", iPlane);
        fhxy[iPlane - 1] = new TH2F(strName1, strName2, 160, -80, 80, 400, -100., 100.);
        fhxy[iPlane - 1]->GetYaxis()->SetTitle("y-position in cm");
        fhxy[iPlane - 1]->GetXaxis()->SetTitle("x-position in cm");
    }
    */
    if (NULL == fhxy[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "xy_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "xy of Plane %d ", iPlane);
        fhxy[iPlane - 1] = new TH2F(strName1, strName2, 125, -62.5, 62.5, 400, -100., 100.);
        fhxy[iPlane - 1]->GetYaxis()->SetTitle("y-position in cm");
        fhxy[iPlane - 1]->GetXaxis()->SetTitle("y-position in cm");
    }
    if (NULL == fhQvsEvent[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "QvsEvent_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Charge vs Event # Plane %d ", iPlane);
        fhQvsEvent[iPlane - 1] = new TH2F(strName1, strName2, 2e2, 0, 2e9, max_charge * 10, 0., max_charge); // 2e4 2e9
        fhQvsEvent[iPlane - 1]->GetYaxis()->SetTitle("Charge");
        fhQvsEvent[iPlane - 1]->GetXaxis()->SetTitle("Event #");
    }
    if (NULL == fhChargePlane[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "Charge_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Charge Plane %d ", iPlane);
        fhChargePlane[iPlane - 1] = new TH1F(strName1, strName2, max_charge * 10, 0., max_charge);
        fhChargePlane[iPlane - 1]->GetXaxis()->SetTitle("Charge");
        fhChargePlane[iPlane - 1]->GetYaxis()->SetTitle("Count");
    }

    // Multiplicity
    if (NULL == fhQM[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "QvsQt0_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Q vs Q_time0 Plane %d ", iPlane);
        fhQM[iPlane - 1] =
            new TH2F(strName1, strName2, max_charge * 10, 0., max_charge, max_charge * 10, 0., max_charge);
        fhQM[iPlane - 1]->GetYaxis()->SetTitle("Charge particle i");
        fhQM[iPlane - 1]->GetXaxis()->SetTitle("Charge first particle");
    }

    if (iPlane == 1 || iPlane == 3)
    {
        if (NULL == fhMvsQ[iPlane - 1])
        {
            char strName1[255];
            sprintf(strName1, "QvsHits_Plane_%d", iPlane);
            char strName2[255];
            sprintf(strName2, "Q vs Hit # Plane %d ", iPlane);
            char strName3[255];
            sprintf(strName3, "Hits in planes %d %d in coincidence window", iPlane, iPlane + 1);
            fhMvsQ[iPlane - 1] = new TH2F(strName1, strName2, max_charge * 10, 0., max_charge, 20, 0., 20);
            fhMvsQ[iPlane - 1]->GetYaxis()->SetTitle(strName3);
            fhMvsQ[iPlane - 1]->GetXaxis()->SetTitle("#sum Charge");
        }
    }

    if (iPlane == 1 || iPlane == 3)
    {
        if (NULL == fhTdiffvsQ[iPlane - 1][2 * iBar - 2])
        {
            char strName[255];
            sprintf(strName, "Tdiff_Plane_%dand%d_Bar_%dvsQ", iPlane, iPlane + 1, iBar * 2 - 1);
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 2] = new TH2F(strName, "", 1000, -10, 10, 1200, 0., 60.);
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 2]->GetYaxis()->SetTitle("charge");
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 2]->GetXaxis()->SetTitle("dt in ns");
        }
        if (NULL == fhTdiffvsQ[iPlane - 1][2 * iBar - 3])
        {
            char strName[255];
            sprintf(strName, "Tdiff_Plane_%dand%d_Bar_%dvsQ", iPlane, iPlane + 1, iBar * 2 - 2);
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 3] = new TH2F(strName, "", 1000, -10, 10, 1200, 0., 60.);
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 3]->GetYaxis()->SetTitle("charge");
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 3]->GetXaxis()->SetTitle("dt in ns");
        }
    }

    if (NULL == fhQvsQ[iPlane - 1][iBar * 2 - 2])
    {
        char strName[255];
        sprintf(strName, "Q_vs_Q_Plane_%d_Bar_%d", iPlane, iBar * 2 - 1);
        fhQvsQ[iPlane - 1][iBar * 2 - 2] = new TH2F(strName, "", 1000, 0, max_charge, 1000, 0., max_charge);
        char strNamex[255];
        char strNamey[255];
        if (iPlane % 2 == 0)
        {
            sprintf(strNamex, "Charge plane %d", iPlane);
            sprintf(strNamey, "Charge plane %d", iPlane + 1);
        }
        else
        {
            sprintf(strNamex, "Charge plane %d", iPlane + 1);
            sprintf(strNamey, "Charge plane %d", iPlane);
        }
        fhQvsQ[iPlane - 1][iBar * 2 - 2]->GetYaxis()->SetTitle(strNamey);
        fhQvsQ[iPlane - 1][iBar * 2 - 2]->GetXaxis()->SetTitle(strNamex);
    }

	if (NULL == fhSqrtQvsPosToT[iPlane - 1][iBar - 1])
	{
        char strName1[255];
		char strName2[255];
        sprintf(strName1, "SqrtQ_vs_PosToT_Plane_%d_Bar_%d", iPlane, iBar);
		sprintf(strName2, "SqrtQ vs PosToT Plane %d Bar %d", iPlane, iBar);
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1] = new TH2F(strName1, strName2, 400, -100, 100, 4000, 2000., 6000);
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("y-position from ToT in cm");
		fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Energy loss in a.u.");
	}

    if (iPlane == 1 || iPlane == 3)
    {
        if (NULL == fhQvsQ[iPlane - 1][iBar * 2 - 3])
        {
            char strName[255];
            sprintf(strName, "Q_vs_Q_Plane_%d_Bar_%d", iPlane, iBar * 2 - 2);
            fhQvsQ[iPlane - 1][iBar * 2 - 3] = new TH2F(strName, "", 1000, 0, max_charge, 1000, 0., max_charge);
            char strNamex[255];
            char strNamey[255];
            if (iPlane % 2 == 0)
            {
                sprintf(strNamex, "Charge plane %d", iPlane);
                sprintf(strNamey, "Charge plane %d", iPlane + 1);
            }
            else
            {
                sprintf(strNamex, "Charge plane %d", iPlane + 1);
                sprintf(strNamey, "Charge plane %d", iPlane);
            }
            fhQvsQ[iPlane - 1][iBar * 2 - 3]->GetYaxis()->SetTitle(strNamey);
            fhQvsQ[iPlane - 1][iBar * 2 - 3]->GetXaxis()->SetTitle(strNamex);
        }
    }

    if (iPlane == 2 || iPlane == 4)
    {
        if (NULL == fhQvsQ[iPlane - 1][iBar * 2 - 1])
        {
            char strName[255];
            sprintf(strName, "Q_vs_Q_Plane_%d_Bar_%d", iPlane, iBar * 2);
            fhQvsQ[iPlane - 1][iBar * 2 - 1] = new TH2F(strName, "", 1000, 0, max_charge, 1000, 0., max_charge);
            char strNamex[255];
            char strNamey[255];
            if (iPlane % 2 == 0)
            {
                sprintf(strNamex, "Charge plane %d", iPlane);
                sprintf(strNamey, "Charge plane %d", iPlane + 1);
            }
            else
            {
                sprintf(strNamex, "Charge plane %d", iPlane + 1);
                sprintf(strNamey, "Charge plane %d", iPlane);
            }
            fhQvsQ[iPlane - 1][iBar * 2 - 1]->GetYaxis()->SetTitle(strNamey);
            fhQvsQ[iPlane - 1][iBar * 2 - 1]->GetXaxis()->SetTitle(strNamex);
        }
    }


    if (NULL == fhCharge)
    {
        char strName[255];
        sprintf(strName, "Charge_of_TofD");
        fhCharge = new TH1F(strName, "Charge of ToFD", 1000, 0., max_charge);
        fhCharge->GetYaxis()->SetTitle("Counts");
        fhCharge->GetXaxis()->SetTitle("Charge");
    }

    if (NULL == fhAverageCharge)
    {
        char strName[255];
        sprintf(strName, "Average_Charge");
        fhAverageCharge = new TH2F(strName, "", 90, 0, 90, 1200, 0., 60.);
        fhAverageCharge->GetYaxis()->SetTitle("Q");
        fhAverageCharge->GetXaxis()->SetTitle("Bar #");
    }

    if (NULL == fhChargevsTof)
    {
        char strName[255];
        sprintf(strName, "Charge_vs_Tof_of_TofD");
        fhChargevsTof = new TH2F(strName, "",  2000, -10., 40.,1000, 0, 100);
        fhChargevsTof->GetXaxis()->SetTitle("ToF in ns");
        fhChargevsTof->GetYaxis()->SetTitle("Charge");
    }

	if (NULL == fhElossVsBar[iPlane - 1])
	{
		fhElossVsBar[iPlane - 1] = new TH2F(Form("ElossVsBar_Plane_%i",iPlane),
						Form("Eloss Vs Bar No Plane %i;Bar No; Eloss / MeV",iPlane),
						44,1,44,10000,0,10000);
	}
	if (NULL ==  fhLuminosity[iPlane - 1][iBar - 1])
	{
		fhLuminosity[iPlane - 1][iBar - 1] = new TH1F(Form("Luminosity_Plane_%i_Bar_%i", iPlane, iBar),
						Form("Luminosity Plane %i Bar %i; Luminosity in a.u.; Entries", iPlane, iBar),
						150,0,1.5);
	}
	if (NULL == fhSqrtToT[iPlane - 1][iBar - 1])
	{
		fhSqrtToT[iPlane - 1][iBar - 1] = new TH1F(Form("SqrtToT_Plane_%i_Bar_%i", iPlane, iBar),
						Form("Sqrt ToT Plane %i, Bar %i", iPlane, iBar),
						300,0.,300.);
	}
	if(NULL == fhPosToT[iPlane - 1])
	{
		fhPosToT[iPlane - 1] = new TH2F(Form("PosToTVsBar_Plane_%i",iPlane),
						Form("PosToT Vs Bar Plane %i; Bar #; PosToT in cm",iPlane),
						44,0,44,1600,-80,80);
	}
	if(NULL == fhPosDiff[iPlane - 1])
	{
		fhPosDiff[iPlane - 1] = new TH2F(Form("PosDiffVsBar_Plane_%i",iPlane),
						Form("PosDiff Vs Bar Plane %i; Bar #; PosDiff in cm",iPlane),
						44,0,44,1600,-80,80);
	}
    
}
void R3BTofdCal2HitS473::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
        fNofHitItems = 0;
    }
    if (fCalItemsLos)
    {
        fCalItemsLos->Clear();
    }
    if (fHitItemsLos)
    {
        fHitItemsLos->Clear();
    }
}

void R3BTofdCal2HitS473::FinishTask()
{
	cout << "Ended Task R3BTofdCal2HitS473." << endl;
	cout << counter << "/" << maxN << " Events processed." << endl;
    if (fhLosXYP)
        fhLosXYP->Write();
    // if (fhChargeLosTofD) fhChargeLosTofD->Write();
    if (fh_los_pos)
        fh_los_pos->Write();
    if (fTofdHisto)
    {
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            if (fhQ[i])
                fhQ[i]->Write();
            if (fhxy[i])
                fhxy[i]->Write();
            if (fhQvsEvent[i])
                fhQvsEvent[i]->Write();
            if (fhQM[i])
                fhQM[i]->Write();
            if (fhMvsQ[i])
                fhMvsQ[i]->Write();
            // if (fhTof[i]) fhTof[i]->Write();
            if (fhTdiff[i])
                fhTdiff[i]->Write();
            if (fhTsync[i]) 
				fhTsync[i]->Write();
			if(fhEloss[i])
				fhEloss[i]->Write();
			if(fhChargePlane[i])
				fhChargePlane[i]->Write();
			if(fhElossVsBar[i])
				fhElossVsBar[i]->Write();
			if(fhPosToT[i])
				fhPosToT[i]->Write();
			if(fhPosDiff[i])
				fhPosDiff[i]->Write();
            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {

                // control histogram time particles
                if (fhQvsPos[i][j])
                    fhQvsPos[i][j]->Write();
                /*
				if (fhTdiffvsQ[i][2 * j])
                    fhTdiffvsQ[i][2 * j]->Write();
                if (fhTdiffvsQ[i][2 * j + 1])
                    fhTdiffvsQ[i][2 * j + 1]->Write();
                
				if (fhQvsQ[i][2 * j])
                    fhQvsQ[i][2 * j]->Write();
                if (fhQvsQ[i][2 * j + 1])
                    fhQvsQ[i][2 * j + 1]->Write();
				*/
                if (fhQvsTof[i][j])
                    fhQvsTof[i][j]->Write();
                if (fhTvsTof[i][j])
                    fhTvsTof[i][j]->Write();
                if (fhToTvsTofw[i][j])
                    fhToTvsTofw[i][j]->Write();
				if(fhSqrtQvsPosToT[i][j])
					fhSqrtQvsPosToT[i][j]->Write();
				if(fhLuminosity[i][j])
					fhLuminosity[i][j]->Write();
				if(fhSqrtToT[i][j])
					fhSqrtToT[i][j]->Write();
            }
        }
        if (fhxy12)
            fhxy12->Write();
        if (fhxy12tot)
            fhxy12tot->Write();
        if (fhxy34)
            fhxy34->Write();
        if (fhxy34tot)
            fhxy34tot->Write();
        if (fhCharge)
            fhCharge->Write();
        if (fhAverageCharge)
            fhAverageCharge->Write();
        if (fhChargevsTof)
			fhChargevsTof->Write();
        if (fhChargevsPos)
			fhChargevsPos->Write();
		if (fhChargeLos)
			fhChargeLos->Write();
        // if (fhQp12) fhQp12->Write();
        // if (fhQp34) fhQp34->Write();
    }
    std::cout << "\n\nSome statistics:\n"
              << "Total number of events in tree  " << maxevent << "\n"
              << "Max Event analyzed              " << fnEvents + wrongtrigger + wrongtpat << "\n"
              << "Events in LOS                   " << countloshit << "\n"
              << "Events skipped due to trigger   " << wrongtrigger << "\n"
              << "Events skipped due to tpat      " << wrongtpat << "\n"
              << "Events with correct header&tpat " << headertpat << "\n"
              << "Events without ToFd hits        " << events_wo_tofd_hits << "\n"
              << "Events in cal level             " << events_in_cal_level << "\n"
              << "Hits in bar coincidence         " << inbarcoincidence << "\n"
              << "Number of resets                " << countreset << "\n"
              << "Hits before reset               " << hitsbeforereset << "\n"
              << "Bars with multihits             " << bars_with_multihit << "\n"
              << "Multihits                       " << multihit << "\n"
              << "Events stored                   " << eventstore / 2 << " <-> " << inbarcoincidence - hitsbeforereset
              << " = " << inbarcoincidence << " - " << hitsbeforereset
              << " (Events in bar coincidence - Hits before reset)\n"
              << "Events in coincidence window    " << incoincidence / 2 << "\n"
              << "Events in single planes         " << singlehit << "\n"
			  << "Events with a factor of 0       " << fact0Counter[0] << ";" << fact0Counter[1] <<
			 		";" << fact0Counter[2] << ";" << fact0Counter[3] << "\n"
			  << "Events with charge ooB          " << chargeCounter[0] << ";" << chargeCounter[1]  <<
			 		";" << chargeCounter[2] << ";" << chargeCounter[3] <<  "\n"
			  << "Hits with radicand < 0           " << negSqrtCounter << "\n"
			  << "LosQ [-inf,-100)/[-100,-0)/{0}/(+0,100]/[100,inf]" << "\n"
			  << losZCounter[0] << "/" << losZCounter[1] << "/" << losZCounter[2] << "/"
			  << losZCounter[3] << "/" << losZCounter[4] << "/" << losZCounter[5] << "\n"
			  << "\n";
}

Double_t R3BTofdCal2HitS473::betaCorr(Double_t delta)
{
    //    Double_t corr=-3.*delta;  //correction for Ag

    Double_t corr = -1. * delta; // correction for 12C
    corr = 0.;
    return corr;
}

Double_t R3BTofdCal2HitS473::walk(Double_t Q,
                              Double_t par1,
                              Double_t par2,
                              Double_t par3,
                              Double_t par4,
                              Double_t par5) // new method
{
    Double_t y = 0;
    y = -30.2 + par1 * TMath::Power(Q, par2) + par3 / Q + par4 * Q + par5 * Q * Q;
    return y;
}

Double_t R3BTofdCal2HitS473::saturation(Double_t x) // not used
{
    Double_t kor;
    Int_t voltage = 600;
    if (voltage == 600)
    {
        if (x < 173)
        {
            kor = 0.;
        }
        else if (x > 208)
        {
            kor = -1.73665e+03 + 2.82009e+01 * 208. - 1.53846e-01 * (208. * 208.) + 2.82425e-04 * (208. * 208. * 208.);
        }
        else
        {
            kor = -1.73665e+03 + 2.82009e+01 * x - 1.53846e-01 * (x * x) + 2.82425e-04 * (x * x * x);
        }
    }
    if (voltage == 500)
    {
        if (x < 95.5)
        {
            kor = 0.;
        }
        else if (x > 124)
        {
            kor = 1.08 * x - 112.44;
        }
        else
        {
            kor = 643.257 - 16.7823 * x + 0.139822 * (x * x) - 0.000362154 * (x * x * x);
        }
    }
    if (voltage == 700)
    {
        if (x < 198)
        {
            kor = 0.;
        }
        else if (x > 298)
        {
            kor = 0.21 * x - 45.54;
        }
        else
        {
            kor = -19067 + 383.93 * x - 3.05794 * (x * x) + 0.0120429 * (x * x * x) - 2.34619e-05 * (x * x * x * x) +
                  1.81076e-08 * (x * x * x * x * x);
        }
    }
    return kor;
}

Double_t* R3BTofdCal2HitS473::insertX(Int_t n, Double_t arr[], Double_t x, Int_t pos)
{
    Int_t i;

    // increase the size by 1
    n++;

    // shift elements forward
    for (i = n; i >= pos; i--)
        arr[i] = arr[i - 1];

    // insert x at pos
    arr[pos - 1] = x;

    return arr;
}

//beta

Double_t R3BTofdCal2HitS473::preBeta(Double_t beta, Double_t Eloss)
{
	return beta;
	/*	Only viable for 120Sn?		
	Double_t E3 = 1/TMath::Sqrt(1-pow(beta,2.))*931.494*119.87477;
	Double_t E4 = E3 + Eloss;
	Double_t beta0 = TMath::Sqrt(1-pow(931.494*119.87477/E4,2.));

	return beta0;
	*/
}

Int_t R3BTofdCal2HitS473::SameVB(hit* h, hit* g)
{
	if(abs(h->p - g->p) == 1)
	{
		if(h->x1 == g->x1) return 1;
		if(h->x1 == g->x2) return 2;
		if(h->x2 == g->x1) return 3;
		if(h->x2 == g->x2) return 4;
	}

	return -1;
}

void R3BTofdCal2HitS473::FindZ(Double_t Z, Double_t* pars)
{
	Double_t offsets[4] = {0,0,0,0};
	Double_t times[4] = {0,0,0,0};
	R3BTofdGlobalParS473* gPar = fGlobPar;
	for(Int_t iPlane = 1; iPlane < 5; iPlane++)
	{
		//Eloss per plane
		TH1F* hist = (TH1F*) fhEloss[iPlane - 1]->Clone();
		hist->Rebin(4);
		Int_t binmax = hist->GetMaximumBin();
		Double_t Max = hist->GetXaxis()->GetBinCenter(binmax);
		TF1* fgaus = new TF1("fgaus", "gaus(0)", Max - 0.3, Max + 0.3);
		hist->Fit("fgaus", "QR0");
		offsets[iPlane - 1] = fgaus->GetParameter(1);		
	
		//Times per plane
		TH2F* histo = (TH2F*) fhTsync[iPlane - 1]->Clone();
		TH1F* histo_py = dynamic_cast<TH1F*>(histo->ProjectionY(Form("ToF_Plane_%i",iPlane),0,-1,""));
		histo_py->Rebin(4);
		binmax = histo_py->GetMaximumBin();
		Max = histo_py->GetXaxis()->GetBinCenter(binmax);
		TF1* ggaus = new TF1("ggaus", "gaus(0)", Max - 0.3, Max + 0.3);
		histo_py->Fit("ggaus", "QR0");
		times[iPlane - 1] = ggaus->GetParameter(1);

		//Charge Factor per plane
		Double_t beta = gPar->GetToFSlope(iPlane) / times[iPlane - 1];
		Double_t charge = beta * TMath::Sqrt(offsets[iPlane - 1] * gPar->GetBetheSlope() + gPar->GetBetheOffset());

		pars[iPlane - 1] = Z / charge;
	}


}

ClassImp(R3BTofdCal2HitS473)
