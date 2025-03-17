//==========================================================================================
//
// Copyright (C) 2021 Dr. Ninad Joshi <email@ninadjoshi.de>  
// P. Auffahrt
// All Rights Reserved
//
//===========================================================================
//==========================================================================================
//
// Copyright (C) 2021 Dr. Ninad Joshi <email@ninadjoshi.de>  
// All Rights Reserved
//
//==========================================================================================

#ifndef FIELDCYL_H_INCLUDED
#define FIELDCYL_H_INCLUDED

#include "vectorcyl.h"
#include <vector>
#include <complex>
#include <string>

//FieldCyl class objects represent scalar fields in cylindrical coordinates with double field values
template <typename T>
class FieldCyl
{
private:
	
	//number of grid points in radial (r), azimuthal (theta) and longitudinal (z) direction
    unsigned int nr{};
    unsigned int nth{};
    unsigned int nz{};
    //vector of values stored in scalar field
    std::vector<T> fieldVal{};

public:

    //constructor with number of grid points in respective directions as parameters
    FieldCyl( unsigned int n1 = 1, unsigned int n2 = 1, unsigned int n3 = 1 );

    //Access functions to number of grid points in r-direction, theta-direction and z-direction
    unsigned int getnr() { return nr; }
    unsigned int getnth(){ return nth; }
    unsigned int getnz() { return nz; }

    //gives value of scalar field at grid point denoted by parameters
    T get( unsigned int i, unsigned int j, unsigned int k );

    //sets the scalar field value at grid point denoted by parameters to parameter val
    void set( unsigned int i, unsigned int j, unsigned int k , T val );
    
    //gives a reference to the value of fieldval for indices i,j,k
    T& ref ( unsigned int i, int j, unsigned int k );

    //creates a copy of the field
    FieldCyl copy();
    
    //sets all the values in fieldVal to zero
    void clear();

    //gives a reference to the value of the field at the grid point denoted by parameters
    T& operator() ( unsigned int i, int j, unsigned int k );
    
    //returns a FieldCyl with all elements in fieldVal sign flipped
    FieldCyl operator- ();
    
    //Adds the values of two fields and returns a FieldCyl
    FieldCyl operator+ ( FieldCyl& f );
    
    //set all elements of fieldVal to one value
    void setVal( T val );
    
    //gives the maximal absolute value stored on any grid point
	double getmax();
	
	//scatters value density of val at position pos onto field
	void scatter ( VectorCyl<double> pos, T val, double dr, double dth, double dz, int nth, bool rshift = false );
	
	//scatters value density of val at position pos onto field for 2D-field (nth = 1)
	void scatter2D ( VectorCyl<double> pos, T val, double dr, double dth, double dz, bool rshift = false );
	
	//scatter value val at position pos onto field defined on staggered grid
	void scatterVal ( VectorCyl<double> pos, T val, double dr, double dth, double dz, int nth, bool rshift = false );
	
	//scatter value val at position pos onto field defined on staggered grid for 2D-field (nth = 1)
	void scatterVal2D ( VectorCyl<double> pos, T val, double dr, double dth, double dz, bool rshift = false );
	
	//gathers value of field at pos and returns that value
	T gather ( VectorCyl<double> pos, double dr, double dth, double dz, int nth );
	
	//gathers value of field at pos and returns that value for 2D-field
	T gather2D ( VectorCyl<double> pos, double dr, double dz );
	
	//writes field data to file with given filepath
	void fieldToFile( std::string filepath, int nrD, int nthD, int nzMin, int nzMax, double dr, double dth, double dz );
	
	//writes 2d-field data in r-z-plane given by j to file
	void fieldToFileRZ( std::string filepath, int nrD, int nthD, int nzMin, int nzMax, double dr, double dth, double dz, int j, double rshift = 0 );
	
	//writes field data of process 0 to file with given filepath from point nz0 to nzD in z-direction
	void fieldToFileSingleProcess( std::string filepath, int w_rank, int nrD, int nthD, int nz0, int nzD, double dr, double dth, double dz );
	
	//writes surfaces of field domain to file
	void fieldToFileSurface( std::string filepath, int w_rank, int nrD, int nthD, int nz0, int nzD, double dr, double dth, double dz );
	
	//writes 2d-field data in r-z-plane given by j to file
    void fieldToFileZR( std::string filepath, int nrD, int nthD, int nzMin, int nzMax, double dr, double dth, double dz, int j );
    
    //writes 2d-field data in r-theta-plane given by k to file
    void fieldToFileThetaR( std::string filepath, int nrD, int nthD, int nzMin, int nzMax, double dr, double dth, double dz, int k );
};

#include "fieldcyl.inl"

#endif // FIELDCYL_H_INCLUDED
