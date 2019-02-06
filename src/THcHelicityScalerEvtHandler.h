#ifndef THcHelicityScalerEvtHandler_
#define THcHelicityScalerEvtHandler_

/////////////////////////////////////////////////////////////////////
//
//   THcHelicityScalerEvtHandler
//   Class to handle Hall C scaler events (type 0)
//   author  Edward Brash (brash@jlab.org)
//   based on THaScalerEvtHandler by Robert Michaels (rom@jlab.org)
//
/////////////////////////////////////////////////////////////////////

#include "THaEvtTypeHandler.h"
#include "Decoder.h"
#include <string>
#include <vector>
#include <set>
#include "TTree.h"
#include "TString.h"
#include <cstring>


class THcHelicityScalerEvtHandler : public THaEvtTypeHandler {

public:

   THcHelicityScalerEvtHandler(const char*, const char*);
   virtual ~THcHelicityScalerEvtHandler();

   Int_t Analyze(THaEvData *evdata);
   Int_t AnalyzeBuffer(UInt_t *rdata, Bool_t onlysync);
   virtual EStatus Init( const TDatime& run_time);
   virtual Int_t   ReadDatabase(const TDatime& date );
   virtual Int_t End( THaRunBase* r=0 );

   virtual void SetUseFirstEvent(Bool_t b = kFALSE) {fUseFirstEvent = b;}
   virtual void SetDelayedType(int evtype);
   virtual void SetROC(Int_t roc) {fROC=roc;}
   virtual void SetBankID(Int_t bankid) {fBankID=bankid;}

private:

   static size_t FindNoCase(const std::string& sdata, const std::string& skey);

   Int_t fNumBCMs;
   Double_t *fBCM_Gain;
   Double_t *fBCM_Offset;
   Double_t *fBCM_delta_charge;
   
   Int_t fROC;
   UInt_t fBankID;

   Double_t fTotalTime;
   Double_t fDeltaTime;
   Double_t fPrevTotalTime;
   Double_t fbcm_Current_Threshold;
   Double_t fClockFreq;
   Int_t fbcm_Current_Threshold_Index;
   std::vector <std::string> fBCM_Name;
   UInt_t evcount;
   Double_t evcountR;
   UInt_t evNumber;
   Int_t Nvars, ifound, fNormIdx, fNormSlot, nscalers;
   Bool_t fUseFirstEvent;
   Bool_t fOnlySyncEvents;
   Bool_t fOnlyBanks;
   Int_t fDelayedType;
   Int_t fClockChan;
   UInt_t fLastClock;
   Int_t fClockOverflows;
   std::vector<UInt_t*> fDelayedEvents;
   std::set<UInt_t> fRocSet;

   THcHelicityScalerEvtHandler(const THcHelicityScalerEvtHandler& fh);
   THcHelicityScalerEvtHandler& operator=(const THcHelicityScalerEvtHandler& fh);

   ClassDef(THcHelicityScalerEvtHandler,0)  // Scaler Event handler

};

#endif
