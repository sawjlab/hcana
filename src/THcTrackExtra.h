#ifndef ROOT_THcTrackExtra
#define ROOT_THcTrackExtra

//////////////////////////////////////////////////////////////////////////
//
// THcTrackExtra
//
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TClonesArray;
class THaTrackingDetector;
class THaTrack;

class THcTrackExtra : public TObject {

public:

  // Default constructor
  THcTrackExtra() 
    : TObject(),
      fTrack(0), fDEDX(kBig), fEnergy(kBig)
  { }

  // Constructor with fp coordinates
  // FIXME: this really should be setting detector coordinates
  THcTrackExtra( THaTrack *track, Double_t dedx=0.0, Double_t energy=0.0)
    : TObject(),
      fTrack(track), fDEDX(dedx), fEnergy(energy)
  { }

  virtual ~THcTrackExtra();

  void              Clear( Option_t* opt="" );

  THaTrack*         GetTrack()   const { return fTrack;}
  Double_t          GetDEDX()    const { return fDEDX; }
  Double_t          GetEnergy()  const { return fEnergy; }

  void              Print( Option_t* opt="" ) const;

  void              SetDEDX(Double_t dedx)     { fDEDX = dedx; }
  void              SetEnergy(Double_t energy) { fEnergy = energy; }

  virtual Bool_t    IsSortable() const { return kFALSE; }
#if 0
  virtual Int_t	    Compare(const TObject* obj) const;
#endif

protected:

  THaTrack*         fTrack;          // Track this information belongs to
  Double_t          fDEDX;           // dEdX from hodoscopes
  Double_t          fEnergy;         // Energy from calorimeter

  static const Double_t kBig;
  
  ClassDef(THcTrackExtra,0)       // A generic particle track
};

#endif
