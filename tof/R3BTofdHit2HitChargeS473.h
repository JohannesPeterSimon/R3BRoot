
#ifndef R3BTOFDHIT2HITCHARGES473
#define R3BTOFDHIT2HITCHARGES473

#define N_TOFD_HIT_PLANE_MAX 4
#define N_TOFD_HIT_PADDLE_MAX 44

#include <map>

#include "FairTask.h"
#include "THnSparse.h"

class TClonesArray;
class R3BTofdHitModulePar;
class R3BTofdHitPar;
class R3BTofdGlobalParS473;
class R3BEventHeader;
class TH1F;
class TH2F;


/**
 * An analysis task to set the highest peak to the correct charge.
 */
class R3BTofdHit2HitChargeS473 : public FairTask
{
  public:

	/**
	 * Default constructor.
	 * Creates an instance of the task with default parameters.
	 */
	R3BTofdHit2HitChargeS473();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTofdHit2HitChargeS473(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTofdHit2HitChargeS473();
  
    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for re-initialization of parameter containers
     * in case the Run ID has changed.
     */
    virtual InitStatus ReInit();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    virtual void SetParContainers();

	/**
	 * Creates all necessary histograms for the given plane and bar.
	 */
    virtual void CreateHistograms(Int_t iPlane, Int_t iBar);
  
    /**
     * Method for setting the nuclear charge of main beam
     */
    inline void SetTofdQ(Double_t Q) { fTofdQ = Q; }

    /**
     * Method for setting histograms
     */
    inline void SetTofdHisto(Bool_t Histo) { fTofdHisto = Histo; }
  
	/**
	 * Method fits highest peak to main charge
	 */
	virtual void FindZ(Double_t Charge, Double_t*);

    /**
     * Method for selecting events with certain trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }
    
	/**
     * Methods for setting number of planes and paddles
     */
    inline void SetNofModules(Int_t planes, Int_t ppp)
    {
        fNofPlanes = planes;
        fPaddlesPerPlane = ppp;
    }
    inline void ReadHitFile(TString file) { fHitFile = file; }

  private:

    TClonesArray* fHitItems;    			/**< Array with Hit items - input data. */
	TClonesArray* fHitChargeItems;			/**< Array with HitCharge items - output data */
    UInt_t fNofHitItems;		        	/**< Number of hit items for cur event. */
	UInt_t fNofHitChargeItems = 0;
    R3BEventHeader* header;  		    	/**< Event header - input data. */
    Double_t fClockFreq;       				/**< Clock cycle in [ns]. */
    Int_t fTrigger;            				/**< Trigger value. */
    Int_t fTpat;
    Double_t fTofdQ;
    Bool_t fTofdHisto;
    UInt_t fnEvents;
    UInt_t fNofPlanes;
    UInt_t fPaddlesPerPlane; 				/**< Number of paddles per plane. */
	Double_t fPeakPos[N_TOFD_HIT_PLANE_MAX];
	TString fHitFile;
	UInt_t maxevent;

	//arrays of control histograms
	TH1F* fhCharges[N_TOFD_HIT_PLANE_MAX];	/**< Charges in single plane */
	TH1F* fhCharge;							/**< Charges all planes */

  public:
    ClassDef(R3BTofdHit2HitChargeS473, 1)
};

#endif
