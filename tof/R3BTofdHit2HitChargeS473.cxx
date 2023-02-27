
#include "R3BTofdHit2HitChargeS473.h"
#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"
#include "R3BTofdHitData.h"
#include "R3BTofdHitChargeData.h"

#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "FairRootManager.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"

#include "TClonesArray.h"
#include "TMath.h"

#include <iostream>

using namespace std;
#ifndef IS_NAN
#define IS_NAN(x) TMath::IsNaN(x)

#define N_TOFD_HIT_PLANE_MAX 4
#define N_TOFD_HIT_PADDLE_MAX 44
#endif

namespace
{
	TFile* hfilename;
}

R3BTofdHit2HitChargeS473::R3BTofdHit2HitChargeS473()
	: FairTask("TofdHit2HitChargeS473",1)
	, fHitItems(NULL)
	, fHitChargeItems(new TClonesArray("R3BTofdHitChargeData"))
	, fNofHitChargeItems(0)
	, fTpat(-1)
	, fTrigger(-1)
    , fNofPlanes(4)
    , fPaddlesPerPlane(44)
    , fTofdQ(1)
    , fTofdHisto(true)
	, fnEvents(0)
	, maxevent(0)
{
	fhCharge = NULL;
	for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
	{
		fPeakPos[i] = 0;
		if(fTofdHisto)
		{
			fhCharges[i] = NULL;
		}
	}
}

R3BTofdHit2HitChargeS473::R3BTofdHit2HitChargeS473(const char* name, Int_t iVerbose)
	: FairTask(name,iVerbose)
	, fHitItems(NULL)
	, fHitChargeItems(new TClonesArray("R3BTofdHitChargeData"))
	, fNofHitChargeItems(0)
	, fTpat(-1)
	, fTrigger(-1)
    , fNofPlanes(4)
    , fPaddlesPerPlane(44)
    , fTofdQ(1)
    , fTofdHisto(true)
	, fnEvents(0)
	, maxevent(0)
{
	fhCharge = NULL;
	for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
	{
		fPeakPos[i] = 0;
		if(fTofdHisto)
		{
			fhCharges[i] = NULL;
		}
	}
}

R3BTofdHit2HitChargeS473::~R3BTofdHit2HitChargeS473()
{
	if(fhCharge)
		delete fhCharge;

	if(fTofdHisto)
	{	
		for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
		{
			if(fhCharges[i])
				delete fhCharges[i];
		}
	}
}

InitStatus R3BTofdHit2HitChargeS473::Init()
{
    // get access to Hit data
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
    fHitItems = (TClonesArray*)mgr->GetObject("TofdHit");
    if (NULL == fHitItems)
        LOG(fatal) << "Branch TofdHit not found";
	hfilename = TFile::Open(fHitFile);
	if(hfilename == 0)
	{
		LOG(fatal) << "Cannot open TofdHit file";
		return kFATAL;
	}
    maxevent = mgr->CheckMaxEventNo();
    mgr->Register("TofdHitCharge", "Land", fHitChargeItems, kTRUE);

	return kSUCCESS;
}

void R3BTofdHit2HitChargeS473::SetParContainers()
{
	FindZ(fTofdQ, fPeakPos);	
}

InitStatus R3BTofdHit2HitChargeS473::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BTofdHit2HitChargeS473::Exec(Option_t* option)
{	
    if (fnEvents / 100000. == (int)fnEvents / 100000)
        std::cout << "\rEvents: " << fnEvents << " / " << maxevent << " (" << (int)(fnEvents * 100. / maxevent)
                  << " %) " << std::flush;

    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
    {
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
            return;
        }
	}

    Int_t nHits = fHitItems->GetEntries();
    Int_t nHitsEvent = 0;

	struct hit
	{
		Double_t x,y,t, Eloss;
		Int_t p;
	};

	for(Int_t iHit = 0; iHit < nHits; iHit++)
	{	
        auto* hit0 = (R3BTofdHitData*)fHitItems->At(iHit);

		hit h;

		h.x = hit0->GetX();
		h.y = hit0->GetY();
		h.t = hit0->GetTime();
		h.p = hit0->GetDetId();
		h.Eloss = hit0->GetEloss();

		Double_t charge = h.Eloss * fPeakPos[h.p - 1];

		fhCharge->Fill(charge);

		if(fTofdHisto)
		{
			fhCharges[h.p - 1]->Fill(charge);
		}

		//Save in new file
		new ((*fHitChargeItems)[fNofHitChargeItems++]) R3BTofdHitChargeData(h.t,
															h.x,
															h.y,
															charge,
															h.Eloss,
															h.p
												);
	}
}

void R3BTofdHit2HitChargeS473::CreateHistograms(Int_t iPlane, Int_t iBar)
{
	if(fhCharge == NULL)
	{
		fhCharge = new TH1F("Charge_Tofd",
						"Charge Tofd;Charge in e;Count",
						320,0.,80.);
	}

	if(fhCharges[iPlane - 1] == NULL)
	{
		fhCharges[iPlane - 1] = new TH1F(Form("Charge_Tofd_Plane_%i",iPlane),
						Form("Charge Tofd Plane %i;Charge in e;Count",iPlane),
						320,0.,80.);
	}
}

void R3BTofdHit2HitChargeS473::FinishEvent()
{
    if (fHitChargeItems)
    {
        fHitChargeItems->Clear();
    }
}

void R3BTofdHit2HitChargeS473::FinishTask()
{
	//Console Output
	if(fhCharge)
		fhCharge->Write();

	if(fTofdHisto)
	{
		for(Int_t i = 0; i < 4; i++)
		{
			if(fhCharges[i])
				fhCharges[i]->Write();
		}
	}
	//Statistics
}

void R3BTofdHit2HitChargeS473::FindZ(Double_t Z, Double_t* pars)
{
	Double_t offsets[4] = {0,0,0,0};
	for(Int_t iPlane = 1; iPlane < 5; iPlane++)
	{
		if(hfilename->Get(Form("Eloss_Plane_%i",iPlane)))
		{
			//Eloss per plane
			TH1F* hist = (TH1F*) hfilename->Get(Form("Eloss_Plane_%i",iPlane))->Clone();
			hist->Rebin(4);
			Int_t binmax = hist->GetMaximumBin();
			Double_t Max = hist->GetXaxis()->GetBinCenter(binmax);
			TF1* fgaus = new TF1("fgaus", "gaus(0)", Max - 0.3, Max + 0.3);
			hist->Fit("fgaus", "QR0");
			offsets[iPlane - 1] = fgaus->GetParameter(1);

			pars[iPlane - 1] = Z/offsets[iPlane - 1];
		}
	}
}

ClassImp(R3BTofdHit2HitChargeS473)
