/** \class THcHelicityScaler
   \ingroup Base

\brief Event handler for Hall C helicity scalers

~~~
~~~
     THcHelcityScaler *hhelscaler = new THcHelicityScaler("H","HC helicity scalers");
     // hscaler->SetDebugFile("HHelScaler.txt");
     hhelscaler->SetROC(8);   // 5 for HMS defaults to 8 for SHMS
     hhelscaler->SetBankID(0x9801); // Will default to this
     gHaEvtHandlers->Add (hhelscaler);
~~~
\author  
*/

//#include "THaEvtTypeHandler.h"
#include "THcHelicityScaler.h"
#include "THaCodaData.h"
#include "THaEvData.h"
#include "THcGlobals.h"
#include "THaGlobals.h"
#include "THcParmList.h"
#include "THcHelicity.h"
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
#include <bitset>
#include <iterator>
#include "THaVarList.h"
#include "VarDef.h"
#include "Helper.h"

using namespace std;
using namespace Decoder;

THcHelicityScaler::THcHelicityScaler(const char *name, const char* description)
  : THaEvtTypeHandler(name,description),
    fBankID(9801),
    fUseFirstEvent(kTRUE),
    fDelayedType(-1),
    fBCM_Gain(0), fBCM_Offset(0)
{
  fROC=-1;
  fNScalerChannels = 32;

  AddEvtType(1);
  AddEvtType(2);
  AddEvtType(4);
  AddEvtType(5);
  AddEvtType(6);
  AddEvtType(7);
  SetDelayedType(129);
  
}

THcHelicityScaler::~THcHelicityScaler()
{
  delete [] fBCM_Gain;
  delete [] fBCM_Offset;

  for( vector<UInt_t*>::iterator it = fDelayedEvents.begin();
       it != fDelayedEvents.end(); ++it )
    delete [] *it;
  fDelayedEvents.clear();
}

Int_t THcHelicityScaler::End( THaRunBase* )
{
  // Process any delayed events in order received

  cout << "THcHelicityScaler::End Analyzing " << fDelayedEvents.size() << " delayed helicity scaler events" << endl;
    for(std::vector<UInt_t*>::iterator it = fDelayedEvents.begin();
      it != fDelayedEvents.end(); ++it) {
    UInt_t* rdata = *it;
    AnalyzeBuffer(rdata);
  }

  for( vector<UInt_t*>::iterator it = fDelayedEvents.begin();
       it != fDelayedEvents.end(); ++it )
    delete [] *it;
  fDelayedEvents.clear();

  //  cout << " -- Helicity Scalers -- " << endl;
  for(Int_t i=0;i<fNScalerChannels;i++) {
    if(fScalerSums[i]>0.5) {
      fAsymmetry[i] = (fHScalers[0][i]-fHScalers[1][i])/fScalerSums[i];
      fAsymmetryError[i] = 2*TMath::Sqrt(fHScalers[0][i]*fHScalers[1][i]
					*fScalerSums[i])
	/(fScalerSums[i]*fScalerSums[i]);
    } else {
      fAsymmetry[i] = 0.0;
      fAsymmetryError[i] = 0.0;
    }
    //    printf("%2d %12.0f %12.0f %12.0f %12.8f\n",i,fScalerSums[i],
    //	   fHScalers[0][i],fHScalers[1][i],
    //	   fAsymmetry[i]);
  }
  //  cout << " ---------------------- " << endl;

  // Compute Charge Asymmetries

  Double_t pclock = fHScalers[0][fClockIndex];
  Double_t mclock = fHScalers[1][fClockIndex];
  cout << endl << "---------------------- Beam Charge Asymmetries ---------------------- " << endl;
  cout << "  BCM        Total     Charge        Beam ON     Beam ON      Asymmetry" << endl;
  cout << " Name       Charge    Asymmetry       Charge    Asymmetry        Error"    << endl;
  for(Int_t i=0;i<fNumBCMs;i++) {
    Int_t index = fBCMtoScalerIndex[i];
    if(index>=0) {
      Double_t pcounts = fHScalers[0][index];
      Double_t mcounts = fHScalers[1][index];
      //      cout << index << " " << fBCM_Name[i] << " " << pcounts << " " << mcounts
      //	   << " " << fBCM_Gain[i]
      //      	   << " " << fBCM_Offset[i] << endl;
      Double_t pcharge = (pcounts - (pclock/fClockFreq)*fBCM_Offset[i])
	/fBCM_Gain[i];
      Double_t mcharge = (mcounts - (mclock/fClockFreq)*fBCM_Offset[i])
	/fBCM_Gain[i];
      fCharge[i] = pcharge+mcharge;
      if(fCharge[i]>0.0) {
	fChargeAsymmetry[i] = (pcharge-mcharge)/fCharge[i];
      } else {
	fChargeAsymmetry[i] = 0.0;
      }
      Double_t asy, asyerr;
      if(fAsymmetryCount[i] <= 1) {
	asy = -100;
	asyerr = 0.0;
      } else {
	asy = fAsymmetrySum[i]/fAsymmetryCount[i];
	if(fAsymmetrySum2[i] >= fAsymmetryCount[i]*asy*asy) {
	  asyerr = TMath::Sqrt((fAsymmetrySum2[i] -
					 fAsymmetryCount[i]*asy*asy) /
					(fAsymmetryCount[i]*(fAsymmetryCount[i]-1)));
	} else {
	  asyerr = 0.0;
	}
      }
      printf("%6s %12.2f %12.8f %12.2f %12.8f %12.8f\n",fBCM_Name[i].c_str(),fCharge[i],
	     fChargeAsymmetry[i],fChargeSum[i],asy,asyerr);
    }
  }
  fTime = (pclock+mclock)/fClockFreq;
  if(pclock+mclock>0) {
    fTimeAsymmetry = (pclock-mclock)/(pclock+mclock);
  } else {
    fTimeAsymmetry = 0.0;
  }
  printf("TIME(s)%12.2f %12.8f %12.2f\n",fTime, fTimeAsymmetry, fTimeSum);
  if(fNTriggersPlus+fNTriggersMinus > 0) {
    fTriggerAsymmetry = ((Double_t) (fNTriggersPlus-fNTriggersMinus))/(fNTriggersPlus+fNTriggersMinus);
  } else {
    fTriggerAsymmetry = 0.0;
  }
  cout << endl << "--------------------------------------------------------------------- " << endl;
  return 0;
}


Int_t THcHelicityScaler::ReadDatabase(const TDatime& date )
{
  char prefix[2];
  prefix[0]='g'; prefix[1]='\0';

  fNumBCMs = 0;
  string bcm_namelist;
  fImin = 2.5;			// Minimum current to calculate a charge asymmetry
  fIminBCM_index = 0;		// Which BCM to use
  DBRequest list[]={
		    {"NumBCMs",&fNumBCMs, kInt, 0, 1},
		    {"Imin_charge_asymmetry", &fImin, kDouble, 0, 1},
		    {"Imin_BCM", &fIminBCM_index, kInt, 0, 1},
		    {"BCM_Names",     &bcm_namelist,       kString},
		    {0}
  };
  gHcParms->LoadParmValues((DBRequest*)&list, prefix);
  if(fNumBCMs > 0) {
    fBCM_Gain = new Double_t[fNumBCMs];
    fBCM_Offset = new Double_t[fNumBCMs];
    DBRequest list2[]={
      {"BCM_Gain",      fBCM_Gain,         kDouble, (UInt_t) fNumBCMs},
      {"BCM_Offset",     fBCM_Offset,       kDouble,(UInt_t) fNumBCMs},
      {0}
    };
    gHcParms->LoadParmValues((DBRequest*)&list2, prefix);
    fBCM_Name = vsplit(bcm_namelist);
  }
  std::map<std::string, Int_t> bcmnametoscalerindex;
  bcmnametoscalerindex["BCM1"] = 0;
  bcmnametoscalerindex["BCM2"] = 2;
  bcmnametoscalerindex["Unser"] = 6;
  bcmnametoscalerindex["BCM4A"] = 10;
  bcmnametoscalerindex["BCM4B"] = 4;
  bcmnametoscalerindex["BCM4C"] = 12;
  //  bcmindex["1MHz"] = 8;
  fClockIndex=8;
  fClockFreq=1000000.0;
  fBCMtoScalerIndex = new Int_t[fNumBCMs];
  for(Int_t i=0;i<fNumBCMs;i++) {
    if(bcmnametoscalerindex.find(fBCM_Name[i]) != bcmnametoscalerindex.end()) {
      Int_t index=bcmnametoscalerindex[fBCM_Name[i]];
      fBCMtoScalerIndex[i] = index;
    } else {
      fBCMtoScalerIndex[i] = -1;
    }
  }
  
  return kOK;
}
void THcHelicityScaler::SetDelayedType(int evtype) {
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
  
Int_t THcHelicityScaler::Analyze(THaEvData *evdata)
{

  if ( !IsMyEvent(evdata->GetEvType()) ) return -1;

  if (fDebugFile) {
    *fDebugFile << endl << "---------------------------------- "<<endl<<endl;
    *fDebugFile << "\nEnter THcHelicityScaler  for fName = "<<fName<<endl;
    EvDump(evdata);
  }

  UInt_t *rdata = (UInt_t*) evdata->GetRawDataBuffer();
  
  if(evdata->GetEvType() == fDelayedType) { // Save this event for processing later
    Int_t evlen = evdata->GetEvLength();
    UInt_t *datacopy = new UInt_t[evlen];
    fDelayedEvents.push_back(datacopy);
    memcpy(datacopy,rdata,evlen*sizeof(UInt_t));
    return 1;
  } else { 			// A normal event
    if (fDebugFile) *fDebugFile<<"\n\nTHcHelicityScaler :: Debugging event type "<<dec<<evdata->GetEvType()<< " event num = " << evdata->GetEvNum() << endl<<endl;
    evNumber=evdata->GetEvNum();
    Int_t ret;
    if((ret=AnalyzeBuffer(rdata))) {
      //
    }
    return ret;
  }

}
Int_t THcHelicityScaler::AnalyzeBuffer(UInt_t* rdata)
{
  fNTrigsInBuf = 0;

  // Parse the data, load local data arrays.
  UInt_t *p = (UInt_t*) rdata;

  UInt_t *plast = p+*p;		// Index to last word in the bank
  Int_t roc = -1;
  Int_t evlen = *p+1;

  Int_t ifound=0;
  
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
//		cout << "ROC: " << roc << " " << evlen << " " << *(p-1) << hex << " " << *p << dec << endl;
	if(roc != fROC) { // Not a ROC with helicity scaler
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
  //          UInt_t num = (*p) & 0xff;
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
	if (roc == fROC) {
	  Int_t nevents = (banklen-2)/fNScalerChannels;
	  //cout << "# of helicity events in bank:" << " " << nevents << endl;
	  if (nevents > 100) {
	    cout << "Error! Beam off for too long" << endl;	
	  }
	  
	  fNTrigsInBuf = 0;
	  // Save helcitiy and quad info for THcHelicity
	  for (Int_t iev = 0; iev < nevents; iev++) {  // find number of helicity events in each bank
	    Int_t index = fNScalerChannels*iev+1;
	    AnalyzeHelicityScaler(p+index);
	    //	    cout << "H: " << evNumber << endl;
	  }
	}
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

  if (!ifound) return 0;

  return 1;

 	
}

Int_t THcHelicityScaler::AnalyzeHelicityScaler(UInt_t *p)
{
  Int_t hbits = (p[0]>>30) & 0x3; // quartet and helcity bits in scaler word
  Bool_t isquartet = (hbits&2) != 0;
  Int_t ispos = hbits&1;
  Int_t actualhelicity = 0;
  fHelicityHistory[fNTrigsInBuf] = hbits;
  fNTrigsInBuf++;
  fNTriggers++;

  Int_t quartetphase = (fNTriggers-fFirstCycle)%4;
  if(fFirstCycle >= -10) {
    if(quartetphase == 0) {
      Int_t predicted = RanBit30(fRingSeed_reported);
      fRingSeed_reported = ((fRingSeed_reported<<1) | ispos) & 0x3FFFFFFF;
      // Check if ringseed_predicted agrees with reported if(fNBits>=30)
      if(fNBits >= 30 && predicted != fRingSeed_reported) {
	cout << "THcHelicityScaler: Helicity Prediction Failed" << endl;
	cout << "Reported  " << bitset<32>(fRingSeed_reported) << endl;
	cout << "Predicted " << bitset<32>(predicted) << endl;
      }
      fNBits++;
      if(fNBits==30) {
	cout << "THcHelicityScaler: A " << bitset<32>(fRingSeed_reported) <<
	  " found at cycle " << fNTriggers << endl;
      }
    } else if (quartetphase == 3) {
      if(!isquartet) {
	cout << "THcHelicityScaler: Quartet bit expected but not set (" <<
	  fNTriggers << ")" << endl;
	fNBits = 0;
	fRingSeed_reported = 0;
	fRingSeed_actual = 0;
	fFirstCycle = -100;
      }
    }
  } else { 			// First cycle not yet identified
    if(isquartet) { // Helicity and quartet signal for next set of scalers
      fFirstCycle = fNTriggers - 3;
      quartetphase = (fNTriggers-fFirstCycle)%4;
      //// Helicity at start of quartet is same as last of quartet, so we can start filling the seed
      fRingSeed_reported = ((fRingSeed_reported<<1) | ispos) & 0x3FFFFFFF;
      fNBits++;
      if(fNBits==30) {
	cout << "THcHelicityScaler: B " << bitset<32>(fRingSeed_reported) <<
	  " found at cycle " << fNTriggers << endl;
      }
    }
  }

  if(fNBits>=30) {
    fRingSeed_actual = RanBit30(fRingSeed_reported);
    fRingSeed_actual = RanBit30(fRingSeed_actual);

#define DELAY9
#ifdef DELAY9
    if(quartetphase == 3) {
      fRingSeed_actual = RanBit30(fRingSeed_actual);
      actualhelicity = (fRingSeed_actual&1)?+1:-1;
    } else {
      actualhelicity = (fRingSeed_actual&1)?+1:-1;
      if(quartetphase == 0 || quartetphase == 1) {
	actualhelicity = -actualhelicity;
      }
    }
    quartetphase = (quartetphase+1)%4;
#else
    actualhelicity = (fRingSeed_actual&1)?+1:-1;
    if(quartetphase == 1 || quartetphase == 2) {
      actualhelicity = -actualhelicity;
    }
#endif
  } else {
    fRingSeed_actual = 0;
  }

  if(actualhelicity!=0) {
    Int_t hindex = (actualhelicity>0)?0:1;
    (actualhelicity>0)?(fNTriggersPlus++):(fNTriggersMinus++);
    Int_t countarray[fNScalerChannels];
    for(Int_t i=0;i<fNScalerChannels;i++) {
      Int_t count = p[i]&0xFFFFFF; // Bottom 24 bits
      countarray[i] = count;
      fHScalers[hindex][i] += count;
      fScalerSums[i] += count;
    }
    // Compute current
    //
    Double_t time = countarray[fClockIndex]/fClockFreq;
    Double_t current = (countarray[fBCMtoScalerIndex[fIminBCM_index]]
			/time-fBCM_Offset[fIminBCM_index])
      /fBCM_Gain[fIminBCM_index];
    //    cout << "Time: " << time << "  Current: " << current << endl;
    if(quartetphase==0) {
      fHaveCycle[0] = fHaveCycle[1] = fHaveCycle[2] = fHaveCycle[3] = kFALSE;
    }
    if(current >= fImin && (quartetphase==0 || fHaveCycle[max(quartetphase-1,0)])) {
      fHaveCycle[quartetphase] = kTRUE;
      for(Int_t i=0;i<fNumBCMs;i++) {
	Int_t index = fBCMtoScalerIndex[i];
	Double_t charge = (countarray[index]
			   -time*fBCM_Offset[i])/fBCM_Gain[i];
	fTimeCycle[quartetphase] = time;
	fChargeCycle[quartetphase][i] = charge;
      }
    }
    //    cout << quartetphase << " " << fHaveCycle[0] << " " << fHaveCycle[1] << " " << fHaveCycle[2] << " " << fHaveCycle[3] << endl;
    if(quartetphase == 3 && fHaveCycle[3]) {	// Compute charge asymmetries for this quartet
      for(Int_t i=0;i<fNumBCMs;i++) {
	Double_t asy = actualhelicity*(fChargeCycle[0][i]+fChargeCycle[3][i]
			-fChargeCycle[1][i]-fChargeCycle[2][i]) /
	  (fChargeCycle[0][i]+fChargeCycle[3][i]+fChargeCycle[1][i]+fChargeCycle[2][i]);
	fChargeSum[i] += fChargeCycle[0][i]+fChargeCycle[1][i]
	  +fChargeCycle[2][i]+fChargeCycle[3][i];
	fAsymmetrySum[i] += asy;
	fAsymmetrySum2[i] += asy*asy;
	fAsymmetryCount[i]++;
      }
      fTimeSum += fTimeCycle[0]+fTimeCycle[1]
	  +fTimeCycle[2]+fTimeCycle[3];
    }
  }
  return(0);
}
//_____________________________________________________________________________
Int_t  THcHelicityScaler::RanBit30(Int_t ranseed)
{
  
  UInt_t bit7    = (ranseed & 0x00000040) != 0;
  UInt_t bit28   = (ranseed & 0x08000000) != 0;
  UInt_t bit29   = (ranseed & 0x10000000) != 0;
  UInt_t bit30   = (ranseed & 0x20000000) != 0;

  UInt_t newbit = (bit30 ^ bit29 ^ bit28 ^ bit7) & 0x1;

  ranseed =  ( (ranseed<<1) | newbit ) & 0x3FFFFFFF;

  return ranseed;

}
//_____________________________________________________________________________
THaAnalysisObject::EStatus THcHelicityScaler::Init(const TDatime& date)
{
  
  ReadDatabase(date);

  fStatus = kOK;

  for( vector<UInt_t*>::iterator it = fDelayedEvents.begin();
       it != fDelayedEvents.end(); ++it )
    delete [] *it;
  fDelayedEvents.clear();

  cout << "Howdy !  We are initializing THcHelicityScaler !!   name =   "
        << fName << endl;

  if(eventtypes.size()==0) {
    eventtypes.push_back(0);  // Default Event Type
  }

  if(fROC < 0) {
    fROC = 8;			// Default to SHMS crate
  }

  fNTriggers = 0;
  fNTrigsInBuf = 0;
  fFirstCycle = -100;
  fRingSeed_reported = 0;
  fRingSeed_actual = 0;
  fNBits = 0;
  fNTriggersPlus = fNTriggersMinus = 0;
  fHScalers[0] = new Double_t[fNScalerChannels];
  fHScalers[1] = new Double_t[fNScalerChannels];
  fScalerSums = new Double_t[fNScalerChannels];
  fAsymmetry = new Double_t[fNScalerChannels];
  fAsymmetryError = new Double_t[fNScalerChannels];
  for(Int_t i=0;i<4;i++) {
    fChargeCycle[i] = new Double_t[fNScalerChannels];
    fHaveCycle[i] = kFALSE;
    for(Int_t j=0;j<fNScalerChannels;j++) {
      fChargeCycle[i][j] = 0.0;
    }
    fTimeCycle[i] = 0.0;
  }
  fFirstHelicity = 0;		// Helicity of first cycle in quartet
  for(Int_t i=0;i<fNScalerChannels;i++) {
    fHScalers[0][i] = 0.0;
    fHScalers[1][i] = 0.0;
    fScalerSums[i] = 0.0;
    fAsymmetry[i] = 0.0;
    fAsymmetryError[i] = 0.0;
  }
  fChargeSum = new Double_t[fNScalerChannels];
  fAsymmetrySum = new Double_t[fNScalerChannels];
  fAsymmetrySum2 = new Double_t[fNScalerChannels];
  fAsymmetryCount = new Int_t[fNScalerChannels];
  for(Int_t i=0;i<fNumBCMs;i++) {
    fChargeSum[i] = 0.0;
    fAsymmetrySum[i] = 0.0;
    fAsymmetrySum2[i] = 0.0;
    fAsymmetryCount[i] = 0;
  }
  fTimeSum = 0.0;

  fCharge = new Double_t[fNumBCMs];
  fChargeAsymmetry = new Double_t[fNumBCMs];

  fTime = fTimeAsymmetry = 0;
  fTriggerAsymmetry = 0.0;

  MakeParms();

  return kOK;
}

void THcHelicityScaler::MakeParms()
{
  /**
     Put Various helicity scaler results in gHcParms so they can be included in results.
  */
  gHcParms->Define(Form("g%s_hscaler_plus[%d]",fName.Data(),fNScalerChannels),
		   "Plus Helcity Scalers",*(fHScalers[0]));
  gHcParms->Define(Form("g%s_hscaler_minus[%d]",fName.Data(),fNScalerChannels),
		   "Minus Helcity Scalers",*(fHScalers[1]));
  gHcParms->Define(Form("g%s_hscaler_sum[%d]",fName.Data(),fNScalerChannels),
		   "Helcity Scalers Sum",*fScalerSums);
  gHcParms->Define(Form("g%s_hscaler_asy[%d]",fName.Data(),fNScalerChannels),
		   "Helicity Scaler Asymmetry[%d]",*fAsymmetry);
  gHcParms->Define(Form("g%s_hscaler_asyerr[%d]",fName.Data(),fNScalerChannels),
		   "Helicity Scaler Asymmetry Error[%d]",*fAsymmetryError);
  gHcParms->Define(Form("g%s_hscaler_triggers",fName.Data()),
		   "Total Helicity Scaler Triggers",fNTriggers);
  gHcParms->Define(Form("g%s_hscaler_triggers_plus",fName.Data()),
		   "Positive Helicity Scaler Triggers",fNTriggersPlus);
  gHcParms->Define(Form("g%s_hscaler_triggers_minus",fName.Data()),
		   "Negative Helicity Scaler Triggers",fNTriggersMinus);
  gHcParms->Define(Form("g%s_hscaler_charge[%d]",fName.Data(),fNumBCMs),
		   "Helicity Gated Charge",*fCharge);
  gHcParms->Define(Form("g%s_hscaler_charge_asy[%d]",fName.Data(),fNumBCMs),
		   "Helicity Gated Charge Asymmetry",*fChargeAsymmetry);
  gHcParms->Define(Form("g%s_hscaler_time",fName.Data()),
		   "Helicity Gated Time (sec)",fTime);
  gHcParms->Define(Form("g%s_hscaler_time_asy",fName.Data()),
		   "Helicity Gated Time Asymmetry",fTimeAsymmetry);
  gHcParms->Define(Form("g%s_hscaler_trigger_asy",fName.Data()),
		   "Helicity Trigger Asymmetry",fTriggerAsymmetry);
}

ClassImp(THcHelicityScaler)

/*

Better charge asymmetry calc.

Compute charge asymmetry by pairs.  Require average current of
pair to be > i_min

error = sqrt( sum_i  (xi-xave)^2) / sqrt(N*N-1)
s(xi^2) - 2 * s(xi) xave + N xave^2
s(xi^2) - N xave*xave

xave = s(xi)/N
s(xi^2) - 2 s(xi) * s(xi)/N + N (s(xi)/N)^2
s(xi^2) - 2 s(xi)^2 / N + s(xi)^2/N
s(xi^2) - s(xi)^2 / N

Accumulate sum of xi, xi**2 
sqrt( sxi2 - sxi*sxi/N ) /sqrt(N*(N-1))

*/
