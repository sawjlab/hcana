//*-- Author :    Ole Hansen   29 March 2001

//////////////////////////////////////////////////////////////////////////
//
// THcTrackExtra
//
// Goes into a TClonesArray parallel with the TClonesArray of THaTracks
// Holds Hall C specific information needed for selecting the golden track
// that doesn't fit in the THaTrack class.
// This can go away either when THaTrack is extended with additional
// data members or a way is made to generate spectrometer specific track
// classes
//
//////////////////////////////////////////////////////////////////////////

#include "THcTrackExtra.h"
#include "THaTrack.h"
#include <iostream>

using namespace std;

const Double_t THcTrackExtra::kBig = 1e38;

//_____________________________________________________________________________
THcTrackExtra::~THcTrackExtra()
{
  // Destructor. Delete the track ID associated with this track.

}

//_____________________________________________________________________________
void THcTrackExtra::Clear( const Option_t* opt )
{
  // If *opt == 'F' then reset all track quantities, else just
  // delete memory managed by this track.
  // (We need this behavior so we can Clear("C") the track TClonesArray
  // without the overhead of clearing everything.)
 
  //FIXME: too complicated. Do we really need to reallocate the trackID?

  if( opt && (*opt == 'F') ) {
    // Initialize data members
  }
}

//_____________________________________________________________________________
void THcTrackExtra::Print( Option_t* opt ) const
{
  // Print track parameters
  TObject::Print( opt );
}

#if 0
// Will need to figure out if we need sort and how it fits in with
// sorting of THaTrack
//_____________________________________________________________________________
Int_t THcTrackExtra::Compare(const TObject * obj) const
{
  // compare two tracks by chi2/ndof
  // for track array sorting

  const THcTrackExtra* tr = dynamic_cast<const THcTrackExtra*>(obj);
  if (!tr) return 0;

  Double_t v1 = GetChi2() / SafeNDoF( GetNDoF() );
  Double_t v2 = tr->GetChi2()/ SafeNDoF( tr->GetNDoF() );

  if( v1<v2 ) return -1;
  else if( v1==v2 ) return 0;
  else return 1;
}
#endif

//_____________________________________________________________________________

ClassImp(THcTrackExtra)

