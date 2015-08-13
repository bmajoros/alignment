/****************************************************************
 Banding.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Banding_H
#define INCL_Banding_H

enum BandingType { 
  NO_BANDING, 
  FIXED_WIDTH_BANDING,    // fixed-width margin around diagonal

  // These are not implemented yet:
  STEPPING_STONE_BANDING, // similar to DoubleScan
  TWAIN_BANDING,          // similar to TWAIN
  STEINER_BANDING         // using a Steiner graph
};


#endif
