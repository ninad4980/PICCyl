//==========================================================================================
//
// Copyright (C) 2021 Dr. Ninad Joshi <email@ninadjoshi.de>  
// All Rights Reserved
//
//==========================================================================================

#include "mathoperations.h"
#include "vectorcyl.h"
#include <cmath>


VectorCyl<double> cylCordTrans( VectorCyl<double> v, double r, double theta, double z )
{
	using std::sin;
	using std::cos;
	using std::sqrt;
	using std::atan2;
	
	//compute Cartesian coordinates with origin at origin of X'
	double xPr{ v.r*cos( v.a ) };
	double yPr{ v.r*sin( v.a ) };
	
	//transform to Cartesian coordinates with origin at origin X
	double x{ xPr + r*cos( theta ) };
	double y{ yPr + r*sin( theta ) };
	double zN{ v.z + z };
	
	//compute cylindrical coordinate vector in coordinate system X
	VectorCyl<double> vec{ sqrt( x*x + y*y ), 0, zN, atan2( y, x ) };
	
	return vec;
}

VectorCyl<double> cylCordTrans2( VectorCyl<double> v, double r, double theta, double z )
{
	using std::sin;
	using std::cos;
	using std::sqrt;
	using std::atan2;
	
	//compute Cartesian coordinates with origin at origin of X
	double x{ v.r*cos( v.a ) };
	double y{ v.r*sin( v.a ) };
	
	//transform to Cartesian coordinates with origin at origin of X'
	double xPr{ x - r*cos( theta ) };
	double yPr{ y - r*sin( theta ) };
	double zPr{ v.z - z };
	
	//compute cylindrical coordinate vector in coordinate sytem X'
	VectorCyl<double> vec{ sqrt( xPr*xPr + yPr*yPr ), 0, zPr, atan2( yPr, xPr ) };
	
	return vec;
}

double numIntAreaCircle( double z1, double z2, double R, int slices )
{
	using std::sqrt;
	double zDiff{ z2 - z1 };
	double deltaZ{ zDiff / slices };
	double reVal{ 0 };
	
	for ( int i{ 0 }; i < slices; ++i )
	{
		double r1{ sqrt( R*R - ( i*deltaZ + z1 )*( i*deltaZ + z1 ) ) };
		double r2{ sqrt( R*R - ( (i+1)*deltaZ + z1 )*( (i+1)*deltaZ + z1 ) ) };
		reVal += deltaZ*r2 + deltaZ*( r1 - r2 )/2; 
	}
	
	return reVal;
}

