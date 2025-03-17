//==========================================================================================
//
// Copyright (C) 2021 Dr. Ninad Joshi <email@ninadjoshi.de>  
// All Rights Reserved
//
//==========================================================================================

#ifndef VECTORFIELDCYL_H_INCLUDED
#define VECTORFIELDCYL_H_INCLUDED

#include "fieldcyl.h"
#include "vectorcyl.h"
#include "constants.h"
#include <string>


//VectorFieldCyl class objects are a 3D vector fields making use of the FieldCyl class
template <typename T>
class VectorFieldCyl
{
private:

	//number of grid points in radial (r), azimuthal (theta) and longitudinal (z) direction
	unsigned int nr{};
	unsigned int nth{};
	unsigned int nz{};
	//simulated part of domain
	double part{};
	//components of vector field in radial, azimuthal and longitudinal direction for every grid point
	FieldCyl<T> vr;
	FieldCyl<T> vth;
	FieldCyl<T> vz;
	//theta-coordinate for every grid point
	FieldCyl<double> angle;
	
public:	
	
	//Constructor
	VectorFieldCyl( unsigned int n1 = 1, unsigned int n2 = 1, unsigned int n3 = 1, double partial = 1 );
	
    //Access functions to number of grid points in r-direction, theta-direction and z-direction
    unsigned int getnr() { return nr; }
    unsigned int getnth() { return nth; }
    unsigned int getnz() { return nz; }
	
	//gives a reference to the r-component of the field at the grid point denoted by parameters
	T& r ( unsigned int i, int j, unsigned int k ) { return vr(i,j,k);	}

	//gives a reference to the theta-component of the field at the grid point denoted by parameters	
	T& th ( unsigned int i, int j, unsigned int k ) { return vth(i,j,k); }

	//gives a reference to the z-component of the field at the grid point denoted by parameters	
	T& z ( unsigned int i, int j, unsigned int k )	{ return vz(i,j,k);	}
	
	//gives a reference to the angle at the grid point denoted by parameters	
	double& a ( unsigned int i, int j, unsigned int k )	{ return angle(i,j,k); }
	
	//returns a vector field normalized to given parameter
	VectorFieldCyl<double> normalize( double normR, double normTh, double normZ, int w_size , int w_rank );
	
	//sets all components to zero for all grid points
	void clear();
	
	//Returns a VectorFieldCyl object with all elements sign flipped
	VectorFieldCyl operator- ();
	
	//Returns a vector consisting of the field components of the field at i,j,k tranformed to angle = 0
	VectorCyl<T> angleZero ( int i, int j, int k );
	
	//Scatter components of vector vec at position pos to field (density)
	void scatter ( VectorCyl<double>& pos, VectorCyl<T>& vec, double dr, double dth, double dz, int nth );
	
	//Scatter value of components of vector vec at position pos to field
	void scatterVal ( VectorCyl<double>& pos, VectorCyl<T>& vec, double dr, double dth, double dz, int nth );
	
	//Gathers the values for components of vector field of circumvening grid points of position pos and returns a vector with correctly transformed components
	VectorCyl<T> gather ( VectorCyl<double>& pos, double dr, double dth, double dz, int nth );
	
	//Gather-function for 2D-fields
	VectorCyl<T> gather2D ( VectorCyl<double>& pos, double dr, double dz );
	
	//Gathers the values for components of staggered grid vector field at b-field definition points and returns a vector with correctly transformed components
	VectorCyl<T> gatherStaggeredB ( VectorCyl<double>& pos, double dr, double dth, double dz, int nth );
	
	//GatherStaggerdB-function for 2D-fields
	VectorCyl<T> gatherStaggeredB2D ( VectorCyl<double>& pos, double dr, double dz );
	
	//Gathers the values for components of staggered grid vector field at dynamic e-field definition points and returns a vector with correctly transformed components
	VectorCyl<T> gatherStaggeredE ( VectorCyl<double>& pos, double dr, double dth, double dz, int nth );
	
	//GatherStaggeredE-function for 2D-fields
	VectorCyl<T> gatherStaggeredE2D ( VectorCyl<double>& pos, double dr, double dz );
	
	//Writes vector field to file in r-theta-plane
	void fieldToFileThetaR( std::string filepath, int nrD, int nthD, int nzMin, int nzMax, double dr, double dth, double dz, int k );
	
	//Writes vector field to file in r-z-plane
	void fieldToFileRZ( std::string filepath, int nrD, int nthD, int nzMin, int nzMax, double dr, double dth, double dz, int j, bool rePart = true );
	
	//Writes Theta-component of fiedl to file in r-z-plane
	void fieldToFileThetaRZ( std::string filepath, int nrD, int nzMin, int nzMax, double dr, double dz, int j = 0, bool rePart = true );
};

#include "vectorfieldcyl.inl"


#endif // VECTORFIELDCYL_H_INCLUDED
