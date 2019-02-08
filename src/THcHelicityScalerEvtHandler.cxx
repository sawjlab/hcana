/** \class THcHelicityScalerEvtHandler
   \ingroup Base

\brief Event handler for Hall C helicity scalers


~~~
     gHaEvtHandlers->Add (new THcHelicityScalerEvtHandler("H","HC helicity scalers"));
~~~
To enable debugging you may try this in the setup script
~~~
     THcScalerEvtHandler *hhelscaler = new THcScalerEvtHandler("H","HC helicity scalers");
     hscaler->SetDebugFile("HHelScaler.txt");
     gHaEvtHandlers->Add (hhelscaler);
~~~
\author  
*/

#include "THaEvtTypeHandler.h"
#include "THcHelicityScalerEvtHandler.h"
#include "THaCodaData.h"
#include "THaEvData.h"
#include "THcParmList.h"
#include "THcGlobals.h"
#include "THaGlobals.h"
#include "TNamed.h"
#include "TMath.h"
#include "TString.h"
#include "TROOT.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <iterator>
#include "THaVarList.h"
#include "VarDef.h"
#include "Helper.h"

using namespace std;
using namespace Decoder;

static const UInt_t ICOUNT    = 1;
static const UInt_t IRATE     = 2;
static const UInt_t ICURRENT = 3;
static const UInt_t ICHARGE   = 4;
static const UInt_t ITIME   = 5;
static const UInt_t ICUT = 6;
static const UInt_t MAXCHAN   = 32;
static const UInt_t defaultDT = 4;

THcHelicityScalerEvtHandler::THcHelicityScalerEvtHandler(const char *name, const char* description)
  : THaEvtTypeHandler(name,description),
    fBCM_Gain(0), fBCM_Offset(0), fBCM_delta_charge(0),
    fBankID(9801),
    evcount(0), evcountR(0.0), ifound(0), fNormIdx(-1),
    fNormSlot(-1),
    fUseFirstEvent(kTRUE),
    fOnlySyncEvents(kFALSE), fOnlyBanks(kFALSE), fDelayedType(-1),
    fClockChan(-1), fLastClock(0), fClockOverflows(0)
{
  fRocSet.clear();
}

THcHelicityScalerEvtHandler::~THcHelicityScalerEvtHandler()
{
  delete [] fBCM_Gain;
  delete [] fBCM_Offset;
  delete [] fBCM_delta_charge;

  for( vector<UInt_t*>::iterator it = fDelayedEvents.begin();
       it != fDelayedEvents.end(); ++it )
    delete [] *it;
  fDelayedEvents.clear();
}

Int_t THcHelicityScalerEvtHandler::End( THaRunBase* )
{
  // Process any delayed events in order received

  cout << "THcHelicityScalerEvtHandler::End Analyzing " << fDelayedEvents.size() << " delayed scaler events" << endl;
  for(std::vector<UInt_t*>::iterator it = fDelayedEvents.begin();
      it != fDelayedEvents.end(); ++it) {
    UInt_t* rdata = *it;
    AnalyzeBuffer(rdata,kFALSE);
  }

  for( vector<UInt_t*>::iterator it = fDelayedEvents.begin();
       it != fDelayedEvents.end(); ++it )
    delete [] *it;
  fDelayedEvents.clear();

  return 0;
}


Int_t THcHelicityScalerEvtHandler::ReadDatabase(const TDatime& date )
{
  char prefix[2];
  prefix[0]='g';
  prefix[1]='\0';
  fNumBCMs = 0;
  DBRequest list[]={
    {"NumBCMs",&fNumBCMs, kInt, 0, 1},
    {0}
  };
  gHcParms->LoadParmValues((DBRequest*)&list, prefix);
  //cout << " NUmber of BCMs = " << fNumBCMs << endl;
  //
  if(fNumBCMs > 0) {
    fBCM_Gain = new Double_t[fNumBCMs];
    fBCM_Offset = new Double_t[fNumBCMs];
    fBCM_delta_charge= new Double_t[fNumBCMs];
    string bcm_namelist;
    DBRequest list2[]={
      {"BCM_Gain",      fBCM_Gain,         kDouble, (UInt_t) fNumBCMs},
      {"BCM_Offset",     fBCM_Offset,       kDouble,(UInt_t) fNumBCMs},
      {"BCM_Names",     &bcm_namelist,       kString},
      {"BCM_Current_threshold",     &fbcm_Current_Threshold,       kDouble,0, 1},
      {"BCM_Current_threshold_index",     &fbcm_Current_Threshold_Index,       kInt,0,1},
      {0}
    };
    fbcm_Current_Threshold = 0.0;
    fbcm_Current_Threshold_Index = 0;
    gHcParms->LoadParmValues((DBRequest*)&list2, prefix);
    vector<string> bcm_names = vsplit(bcm_namelist);
    for(Int_t i=0;i<fNumBCMs;i++) {
      fBCM_Name.push_back(bcm_names[i]+".scal");
      fBCM_delta_charge[i]=0.;
    }
  }
  fTotalTime=0.;
  fPrevTotalTime=0.;
  fDeltaTime=-1.;
  //
  //
  return kOK;
}
void THcHelicityScalerEvtHandler::SetDelayedType(int evtype) {
  /**
   * \brief Delay analysis of this event type to end.
   *
   * Final scaler events generated in readout list end routines may not
   * come in order in the data stream.  If the event type of a end routine
   * scaler event is set, then the event contents will be saved and analyzed
   * at the end of the analysis so that time ordering of scaler events is preserved.
   */
  fDelayedType = evtype;
}
  
Int_t THcHelicityScalerEvtHandler::Analyze(THaEvData *evdata)
{

  if ( !IsMyEvent(evdata->GetEvType()) ) return -1;

  if (fDebugFile) {
    *fDebugFile << endl << "---------------------------------- "<<endl<<endl;
    *fDebugFile << "\nEnter THcHelicityScalerEvtHandler  for fName = "<<fName<<endl;
    EvDump(evdata);
  }

  UInt_t *rdata = (UInt_t*) evdata->GetRawDataBuffer();

  if( evdata->GetEvType() == fDelayedType) { // Save this event for processing later
    Int_t evlen = evdata->GetEvLength();
    
    UInt_t *datacopy = new UInt_t[evlen];
    fDelayedEvents.push_back(datacopy);
    memcpy(datacopy,rdata,evlen*sizeof(UInt_t));
    return 1;
  } else { 			// A normal event
    if (fDebugFile) *fDebugFile<<"\n\nTHcHelicityScalerEvtHandler :: Debugging event type "<<dec<<evdata->GetEvType()<< " event num = " << evdata->GetEvNum() << endl<<endl;
    evNumber=evdata->GetEvNum();
    Int_t ret;
    if((ret=AnalyzeBuffer(rdata,fOnlySyncEvents))) {
      //
    }
    return ret;

  }

}
Int_t THcHelicityScalerEvtHandler::AnalyzeBuffer(UInt_t* rdata, Bool_t onlysync)
{

  // Parse the data, load local data arrays.
  UInt_t *p = (UInt_t*) rdata;

  UInt_t *plast = p+*p;		// Index to last word in the bank
  Int_t roc = -1;
  Int_t evlen = *p+1;

  ifound=0;
  while(p<plast) {
    Int_t banklen = *p;
    p++;			  // point to header
    if (fDebugFile) {
      *fDebugFile << "Bank: " << hex << *p << dec << " len: " << *(p-1) << endl;
    }
    if((*p & 0xff00) == 0x1000) {	// Bank Containing banks
      if(evlen-*(p-1) > 1) { // Don't use overall event header
        roc = (*p>>16) & 0xf;
	if(fDebugFile) *fDebugFile << "ROC: " << roc << " " << evlen << " " << *(p-1) << hex << " " << *p << dec << endl;
	//	cout << "ROC: " << roc << " " << evlen << " " << *(p-1) << hex << " " << *p << dec << endl;
	if(fRocSet.find(roc)==fRocSet.end()) { // Not a ROC with helicity scaler
	  p+=*(p-1)-1;		// Skip to next ROC
	}
      }
      p++;				// Now pointing to a bank in the bank
    } else if (((*p & 0xff00) == 0x100) && (*p != 0xC0000100)) {
      // Bank containing integers.  Look for scalers
      // This is either ROC bank containing integers or
      // a bank within a ROC containing data from modules of a single type
      // Look for banks with the helicity scaler bank ID (9801)
      // Assume that very first word is a scaler header
      // At any point in the bank where the word is not a matching
      // header, we stop.
      UInt_t tag = (*p>>16) & 0xffff; // Bank ID (ROC #)
      //      UInt_t num = (*p) & 0xff;
      UInt_t *pnext = p+banklen;	// Next bank
      p++;			// First data word

      // If the bank is not a helicity scaler bank
      // or it is not one of the ROC containing helcity scaler data
      // skip to the next bank
      //cout << "BankID=" << tag << endl;
      if (tag != fBankID) {
	p = pnext;		// Fall through to end of the above else if
	//	cout << "  Skipping to next bank" << endl;

      } else {
	// This is a helicity scaler bank
	cout << "ROC " << roc << "  Scaler bank " << tag << ": " << banklen << endl;
	cout << hex;
	for(Int_t i=0;i<banklen-1;i++) {
	  if(i%5==0) {
	    cout << dec << i << ":";
	  }
	  cout << " " << hex << p[i];
	  if(i%5==4) {
	    cout << endl;
	  }
	}
	cout << dec << endl;
      }
      while(p < pnext) {
	Int_t nskip = 0;
	if(fDebugFile) {
	  *fDebugFile << "Scaler Header: " << hex << *p << dec;
	}
	if(nskip == 0) {
	  if(fDebugFile) {
	    *fDebugFile << endl;
	  }
	  break;	// Didn't find a matching header
	}
	p = p + nskip;
      }
      p = pnext;
    } else {
      p = p+*(p-1);		// Skip to next bank
    }
  }

  if (fDebugFile) {
    *fDebugFile << "Finished with decoding.  "<<endl;
    *fDebugFile << "   Found flag   =  "<<ifound<<endl;
  }

  // HMS has headers which are different from SOS, but both are
  // event type 0 and come here.  If you found no headers, return.

  if (!ifound) return 0;

  evcount = evcount + 1;
  evcountR = evcount;
  //
  
  return 1;
}


THaAnalysisObject::EStatus THcHelicityScalerEvtHandler::Init(const TDatime& date)
{
  //
  ReadDatabase(date);

  fStatus = kOK;
  fNormIdx = -1;

  for( vector<UInt_t*>::iterator it = fDelayedEvents.begin();
       it != fDelayedEvents.end(); ++it )
    delete [] *it;
  fDelayedEvents.clear();

  cout << "Howdy !  We are initializing THcHelicityScalerEvtHandler !!   name =   "
        << fName << endl;

  if(eventtypes.size()==0) {
    eventtypes.push_back(0);  // Default Event Type
  }

  fRocSet.insert(5);		// List ROCs that have helicity scalers
  fRocSet.insert(8);		// Should make configurable

  //
  return kOK;
}

size_t THcHelicityScalerEvtHandler::FindNoCase(const string& sdata, const string& skey)
{
  // Find iterator of word "sdata" where "skey" starts.  Case insensitive.
  string sdatalc, skeylc;
  sdatalc = "";  skeylc = "";
  for (string::const_iterator p =
	 sdata.begin(); p != sdata.end(); ++p) {
    sdatalc += tolower(*p);
  }
  for (string::const_iterator p =
	 skey.begin(); p != skey.end(); ++p) {
    skeylc += tolower(*p);
  }
  if (sdatalc.find(skeylc,0) == string::npos) return -1;
  return sdatalc.find(skeylc,0);
};

ClassImp(THcHelicityScalerEvtHandler)
