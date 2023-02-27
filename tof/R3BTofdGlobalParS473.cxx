#include "R3BTofdGlobalParS473.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TF1.h"
#include "TH1F.h"
#include "TPad.h"

using namespace std;

ClassImp(R3BTofdGlobalParS473);

R3BTofdGlobalParS473::R3BTofdGlobalParS473(const char* name, const char* title, const char* context, Bool_t own)
	: FairParGenericSet(name, title, context, own)
{
	//Reset all parameters
	clear();
}

R3BTofdGlobalParS473::~R3BTofdGlobalParS473() {}

void R3BTofdGlobalParS473::putParams(FairParamList* list)
{
	LOG(info) << "R3BTofdGlobalParS473::putParams() called";
	if (!list)
	{
		return;
	}
	list->add("BetheSlope", fBetheSlope);
	list->add("BetheOffset",fBetheOffset);
	//list->add("fToT", fTTofdOffset);
	for(Int_t i = 0; i < 4; i++)
	{
		TString name = Form("Eloss%01d",i+1);
		list->add(name, fEloss[i]);

		name = Form("ZDist%01d",i+1);
		list->add(name, fZDist[i]);

		name = Form("TOffset%01d",i+1);
		list->add(name, fTOffset[i]);

	}
}

Bool_t R3BTofdGlobalParS473::getParams(FairParamList* list)
{
	if (!list)
	{
		return kFALSE;
	}
	if (!list->fill("BetheSlope", &fBetheSlope)) return kFALSE;
	if (!list->fill("BetheOffset", &fBetheOffset)) return kFALSE;
	//if (!list->fill("fToT", &fTTofdOffset)) return kFALSE;
	for (Int_t i = 0; i < 4; i++)
	{
		TString name = Form("Eloss%01d",i+1);
		if(list->fill(name, &fEloss[i])) return kFALSE;

		name = Form("ZDist%01d",i+1);
		if(list->fill(name, &fZDist[i])) return kFALSE;

		name = Form("TOffset%01d",i+1);
		if(list->fill(name, &fTOffset[i])) return kFALSE;

	}

	return kTRUE;
}

void R3BTofdGlobalParS473::clear()
{
	fBetheSlope = fBetheOffset = 0;//fToFOffset = fToFSlope = 0;//fTTofdOffset = 0.;
	for(Int_t i = 0; i < 4; i++)
	{
		fEloss[i] = fZDist[i] = fTOffset[i] = fToFOffset[i] = fToFSlope[i] =  0.;
	}
}


void R3BTofdGlobalParS473::printParams()
{
	LOG(info) << "   R3BTofdGlobalParS473:  ToFD Velocity Parameters: ";
	//LOG(info) << "   fToF: " << fTTofdOffset;
	for(Int_t i = 0; i<4; i++){
		LOG(info) << "   fTOffset" << i+1 << ": "  << fTOffset[i];
	}
	for(Int_t i = 0; i<4; i++){
		LOG(info) << "   fZDist" << i+1 << ": "  << fZDist[i];
	}
	for(Int_t i = 0; i<4; i++){
		LOG(info) << "   fEloss" << i+1 << ": "  << fEloss[i];
	}
	for(Int_t i = 0; i<4; i++){
		LOG(info) << "   fToFOffset" << i+1 << ": "  << fToFOffset[i];
	}
	for(Int_t i = 0; i<4; i++){
		LOG(info) << "   fToFSlope" << i+1 << ": "  << fToFSlope[i];
	}

	LOG(info) << "   fBetheSlope: " << fBetheSlope;
	LOG(info) << "   fBetheOffset: " << fBetheOffset;
	//LOG(info) << "	 fToFSlope: " << fToFSlope;
	//LOG(info) << "	 fToFOffset: " << fToFOffset;
}

void R3BTofdGlobalParS473::DrawParams()
{
	// Is this even necessary?
}
