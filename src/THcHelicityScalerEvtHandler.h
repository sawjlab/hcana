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
#include <algorithm>
#include <set>
#include "TTree.h"
#include "TString.h"
#include <cstring>

class THcHelicity;

class THcHelicityScalerEvtHandler : public THaEvtTypeHandler {

public:

   THcHelicityScalerEvtHandler(const char*, const char*);
   virtual ~THcHelicityScalerEvtHandler();

   Int_t Analyze(THaEvData *evdata);
   Int_t AnalyzeBuffer(UInt_t *rdata, Bool_t onlysync);
   Int_t AnalyzeHelicityScaler(UInt_t *p);
	
   virtual EStatus Init( const TDatime& run_time);
   virtual Int_t   ReadDatabase(const TDatime& date );
   virtual Int_t End( THaRunBase* r=0 );

   virtual void SetUseFirstEvent(Bool_t b = kFALSE) {fUseFirstEvent = b;}
   virtual void SetDelayedType(int evtype);
   virtual void SetROC(Int_t roc) {fROC=roc;}
   virtual void SetBankID(Int_t bankid) {fBankID=bankid;}
   virtual void SetHelicityDetector(THcHelicity *f) {fglHelicityDetector = f;}

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
   std::vector <int> eventnumbers; // running storage of event numbers
   Int_t nevents; // # of helicity windows in each helicity bank
   Int_t DAQ_rep_hel_windows(Int_t nevents); // value of helicity in each window of DAQ reported helicity
   std::vector <int> DAQ_rep_hel_bank; // running storage of DAQ reported helicity values
   std::vector <int> DAQ_pred_hel_bank; // running storage of DAQ predicted helicity values
   std::vector <int> DAQ_act_hel_bank; // running storage of DAQ actual helicity values
   std::vector <Int_t> random_seed; // random seed for 30-bit (120 windows) shift-register
   Int_t pos_accu_scaler;
   Int_t neg_accu_scaler;
   std::vector <int> actual_helicity; // actual helicity aligned with DAQ (beginning at the appropriate event)
   std::vector <int> positive_hel_scalers; // running storage of accumulated scalers (+ helicity) for each event
   std::vector <int> negative_hel_scalers; // running storage of accumulated scalers (- helicity) for each event
   Int_t bit1;
   Int_t bit7;
   Int_t bit28;
   Int_t bit29;
   Int_t bit30;
   Int_t newbit;
   Int_t next_quartet[4];
   std::vector<UInt_t*> fDelayedEvents;
   std::set<UInt_t> fRocSet;

   THcHelicityScalerEvtHandler(const THcHelicityScalerEvtHandler& fh);
   THcHelicityScalerEvtHandler& operator=(const THcHelicityScalerEvtHandler& fh);

   ClassDef(THcHelicityScalerEvtHandler,0)  // Scaler Event handler
   Int_t quartet_1[4];  // for finding and checking quartet pattern
   Int_t quartet_2[4];

   THcHelicity *fglHelicityDetector;
   

};

#endif
