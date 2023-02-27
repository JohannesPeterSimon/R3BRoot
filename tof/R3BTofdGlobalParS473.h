#ifndef R3BTOFDGLOBALPARS473_H
#define R3BTOFDGLOBALPARS473_H

#include "FairParGenericSet.h"

/**
 * Container for velocity calibration for ToFD. This class is used to store TofD velocity parameters.
 */

class R3BTofdGlobalParS473 : public FairParGenericSet
{
	public: 
	/**
	 * Standard Constructor
	 * @param name a name of container
	 * @param title a title of container
	 * @param context context/purpose for parameters and conditions
	 * @param own class ownership, if flag is kTRUE FairDB has the par
	 */
	R3BTofdGlobalParS473(const char* name = "TofdGlobalPar",
					const char* title = "Tofd Velocity Calibration Parameters",
					const char* context = "TestDefaultContext",
					Bool_t own = kTRUE);

	/**
	 * Destructor
	 * Frees the memory allocated by this object
	 */
	virtual ~R3BTofdGlobalParS473(void);

	/**
	 * A method to reset all parameter values.
	 */
	void clear(void);

	/**
	 * A method to write parameters using RuntimeDB
	 * @param list a list of parameters
	 */
	void putParams(FairParamList* list);

	/** 
	 * A method to read parameters using RuntimeDB
	 * @param list a list of parameters
	 * @return kTRUE if successfull, else kFALSE
	 */
	Bool_t getParams(FairParamList* list);

	/**
	 * A method to print values of parameters to the standard output using FairLogger
	 */
	void printParams();

	/**
	 * A method to Draw values of parameters to the current canvas
	 */
	void DrawParams();

	/** Accesssor Functions **/
	Double_t GetZDist(Int_t plane) const { return fZDist[plane - 1]; }
	Double_t GetTOffset(Int_t plane) const { return fTOffset[plane -1]; }
	Double_t GetEloss(Int_t plane) const { return fEloss[plane - 1]; }
	Double_t GetBetheSlope() const { return fBetheSlope; }
	Double_t GetBetheOffset() const { return fBetheOffset; }
	Double_t GetToFOffset(Int_t plane) const {return fToFOffset[plane - 1]; }
	Double_t GetToFSlope(Int_t plane) const {return fToFSlope[plane - 1]; }

	void SetBetheOffset(Double_t offset) { fBetheOffset = offset; }
	void SetBetheSlope(Double_t slope) { fBetheSlope = slope; }
	void SetToFSlope(Int_t plane, Double_t tofSlope) { fToFSlope[plane - 1] = tofSlope; }
	void SetToFOffset(Int_t plane, Double_t tofOffset) { fToFOffset[plane - 1] = tofOffset; }
	void SetEloss(Int_t plane, Double_t eloss) { fEloss[plane -1] = eloss; }
	void SetTOffset(Int_t plane, Double_t t) { fTOffset[plane -1] = t; }
	void SetZDist(Int_t plane, Double_t z) {fZDist[plane -1] = z; }


   private:

	Double_t fZDist[4]; 	/**< Distance of a plane to plane 1. */
	Double_t fTOffset[4]; 	/**< ToF between a plane and plane 1. */
	Double_t fEloss[4];		/**< Atima Eloss in a plane. */
	Double_t fBetheSlope;	/**< Slope of the Bethe fit. */
	Double_t fBetheOffset;	/**< Offset of the Bethe fit. */
	Double_t fToFSlope[4];		/**< Slope of the velocity dependent ToF offset. */
	Double_t fToFOffset[4];		/**< Offset of the velocity dependent ToF offset. */

	ClassDef(R3BTofdGlobalParS473,1);
};

#endif
