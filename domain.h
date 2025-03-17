//==========================================================================================
//
// Copyright (C) 2021 Dr. Ninad Joshi <email@ninadjoshi.de>  
// All Rights Reserved
//
//==========================================================================================

#ifndef DOMAIN_H_INCLUDED
#define DOMAIN_H_INCLUDED

#include "vectorcyl.h"
#include "vectorfieldcyl.h"
#include "fieldcyl.h"
#include "species.h"
#include "timer.h"
#include <string>
#include <complex>

using Complex = std::complex<double>;

//forward declaration of Species class for member function declaration of Domain class
class Species;

//Domain class objects are storing field values for all defined grid points and provide functions to evaluate and access those fields
class Domain
{
private:
	
	//set to true for a 2D-simulation
	bool dim2{ false };
	//number of grid points in radial (r), azimuthal (theta) and longitudinal (z) direction
	int nr{};
	int nth{};
	int nz{};
	//number of first and last grid point of process for z-decomposition of domain
	int nZmin{};
	int nZmax{};
	//radius and height of cylindrical domain
	double R{};
	double H{};
	//radius and length of plasma chamber
	double chamberRadius{};
	double chamberLength{};
	//simulated part of domain
	double part{};
	//distance between two grid points in radial (in m), azimuthal (in rad) and longitudinal (in m) direction
	double dr{};
	double dth{};
	double dz{};
	//time step of simulation
	double dt{};
	//time step count
	int it{ 0 };
	//time steps per rf-cycle
	int rfsteps{};
	//minimal and maximal z-value covered by process
	double zMin{};
	double zMax{};
	//rank of process and total number of processes in mpi-decomposition
	int w_rank{};
	int w_size{};
	//coil current frequency
	double omega{};
	//smallest radial distance considered to be unequal zero
	double zeroDistanceLimit{ 1e-10 };
	//maximal materialIndex for computation of potential
	int evIndex{ 2 };
	//constant magnetic and electric field in cylindrical components for test purposes
	VectorCyl<double> B{};
	VectorCyl<double> E{};
	//material index for every grid point
	FieldCyl<int> materialIndex{};
	//el, potential, charge density and curent density for every grid point
	FieldCyl<double> phi{};
	FieldCyl<double> rho{};
	VectorFieldCyl<Complex> jCurr{};
	VectorFieldCyl<Complex> coilCurr{};
	//charge density of electrons and ions for every grid point
	FieldCyl<double> elChargeDensity{};
	FieldCyl<double> ioChargeDensity{};
	//density and drift velocity of neutral gas for every grid point
	FieldCyl<double> neutralGasDensity{};
	VectorFieldCyl<double> neutralGasDriftVel{};
	//electric and magnetic field for every grid point
	VectorFieldCyl<double> Efield{};
	VectorFieldCyl<Complex> Bfield{};
	VectorFieldCyl<Complex> Edynfield{};
	//timer for domain
	Timer dTimer{};
	//test fields
	FieldCyl<double> phi1{};
	FieldCyl<double> phi2{};
	FieldCyl<double> phi3{};
	FieldCyl<double> phi4{};
	FieldCyl<double> phi5{};
	FieldCyl<double> phi6{};
	
public:

	//Constructor of cylindrical domain, taking node numbers in r-, theta-, z-direction, radius and height of domain as input
	Domain ( int n1, int n2, int n3, double radius, double height, double partial, double dT, double frequency );
	
	//Access functions
	bool get2D() { return dim2; }
	int getnr() { return nr; }
	int getnth() { return nth; }
	int getnz() { return nz; }
	int getnZmin() { return nZmin; }
	int getnZmax() { return nZmax; }
	double getR() { return R; }
	double getH() { return H; }
	double getChRadius() { return chamberRadius; }
	double getChLength() { return chamberLength; }
	double getpart() { return part; }
	double getdr() { return dr; }
	double getdth() { return dth; }
	double getdz() { return dz; }
	double getdt() { return dt; }
	int getit() { return it; }
	int getrfsteps() { return rfsteps; }
	double getzMin() { return zMin; }
	double getzMax() { return zMax; }
	int getrank() { return w_rank; }
	int getsize() { return w_size; }
	double getomega() { return omega; }
	double getzDL() { return zeroDistanceLimit; }
	VectorCyl<double> getB() { return B; }
	VectorCyl<double> getE() { return E; }
	VectorFieldCyl<double> getEfield() { return Efield; }
	
	//References to field values
	int& matIndexVal( int i, int j, int k ) { return materialIndex(i,j,k); }
	double& phiVal( int i, int j, int k ) { return phi(i,j,k); }
	double& rhoVal( int i, int j, int k ) { return rho(i,j,k); }
	double& Er( int i, int j, int k ) { return Efield.r(i,j,k); }
	double& Eth( int i, int j, int k ) { return Efield.th(i,j,k); }
	double& Ez( int i, int j, int k ) { return Efield.z(i,j,k); }
	
	//References to fields
	FieldCyl<double>& rhoField() { return rho; }
	FieldCyl<double>& elChDen() { return elChargeDensity; }
	FieldCyl<double>& ioChDen() { return ioChargeDensity; }
	FieldCyl<double>& nGasDen() { return neutralGasDensity; }
	VectorFieldCyl<Complex>& currDen() { return jCurr; }
	VectorFieldCyl<Complex>& cCurrent() { return coilCurr; }
	VectorFieldCyl<double>& nGasVel() { return neutralGasDriftVel; }
	VectorFieldCyl<double>& eField() { return Efield; }
	VectorFieldCyl<Complex>& bField() { return Bfield; }
	VectorFieldCyl<Complex>& dynEField() { return Edynfield; }
	Timer& timer() { return dTimer; }
	
	//Value Setting
	void setChamberDimension( double radius, double length );
	void setzDL( double zDL ) { zeroDistanceLimit = zDL; }
	void setB( VectorCyl<double> bF ) { B = bF; }
	void setE( VectorCyl<double> eF ) { E = eF; }
	
	//Split domain over different processes along z direction
	void splitDomainMPI ( int world_rank, int world_size );
	
	//Output information about simulation
	void endSimulation();
	
	//Checks if given position is within domain
	bool inBounds( VectorCyl<double>& pos );
	
	//Advance time step count
	void advanceDomain() { ++it; }
	
	//Compute exp() of current phase of coil current
	Complex computeExpPhase();
	
	//Compute exp() of negative coil current phase
	Complex computeExpNegPhase();
	
	//Solves Poisson equation for electrostatic potential with Gauss-Seidel method
	void poissonGSSolver( double tol );
	
	//Solves Poisson equation for electrostatic potential with Gauss-Seidel method in 2 dimensions
	void poissonGSSolver2D( double tol );
		
	//Computes the potential phi via Jacobi method with charge density rho and clears rho after that
	void jacobiPotentialSolver ( double tol , int potSet = 0 );
	//for potSet = 1 a double disc potential is preset (test potential)
	void discPot();				//not in use anymore
	//for potSet = 2 a double cylinder potential is preset (test potential)
	void cylPot();				//not in use anymore
	
	//set up pot for Multigrid
	void setUpMultiGrid();
	
	//compute tol for multigrid
	double computeTol( double refTol );
	
	//test function for multigrid solver
	void testMultigrid( double tol );
	
	//Computes complex amplitude of dynamic e-field 
	void dynFieldSolver( double tol );
	
	//Computes the gradient of the given scalar field and returns it as a VectorFieldCyl class object
	VectorFieldCyl<double> grad( FieldCyl<double>& sField );
	
	//Computes gradient of a scalar field and gives out a 2D-fiel in r-z-plane
	VectorFieldCyl<double> grad2D( FieldCyl<double>& sField );
	
	//Computes curl of E-field defined on staggered grid (Yee grid)
	template <typename T>
	VectorFieldCyl<T> curlE( VectorFieldCyl<T>& vField );
	
	//Computes curl of E-field in 2D-r-z-plane
	template <typename T>
	VectorFieldCyl<T> curlE2D( VectorFieldCyl<T>& vField );
	
	//Computes curl of E-field in 3D
	template <typename T>
	VectorFieldCyl<T> curlE3D( VectorFieldCyl<T>& vField );
	
	//Computes curl of B-field defined on staggered grid
	template <typename T>
	VectorFieldCyl<T> curlB( VectorFieldCyl<T>& vField );
	
	//Computes curl of B-field in 2D-r-z-plane
	template <typename T>
	VectorFieldCyl<T> curlB2D( VectorFieldCyl<T>& vField );
	
	//Computes curl of B-field in 3D
	template <typename T>
	VectorFieldCyl<T> curlB3D( VectorFieldCyl<T>& vField );
	
	//Computes complex amplitude of dynamic e-field for boundary points
	void computeEdynBiotSavart();
	
	//Computes e-field with vector potential solution at given grid poin
	VectorCyl<Complex> biotSavartIntegral( int iPr, int jPr, int kPr );
	
	//Shares the current density data between different processes
	void shareCurrentData();
	
	//Compute electromagnetic fields
	void computeEMFields( double tolStatic, double tolDyn );
	
	//Computes complex amplitude of b-field
	void computeBComplex();
	
	//Applies boundary conditions for E-field
	void applyBoundaries();
	
	//Computes B-field with Maxwell-Faraday equation
	void computeBFaradaysLaw();
	
	//Computes E-field with Ampère-Maxwell equation
	void computeEAmperesLaw();
	
	//Computes initial B-field 
	void computeBInit();
	
	//Computes initial EM-fields
	void initialEMSolve();
	
	//Computes E-field and B-field for new time step
	void advanceFields();

	//Create electrons and ions to have quasi neutrality in plasma and maintain minDensity
	void retainDensityQNRIT4( Species& electrons, Species& ions, double Te, double Ti, double minDensity );
	
	//Create electrons and ions to have quasi neutrality in plasma and maintain minDensity in RIT2.5
	void retainDensityQNRIT25( Species& electrons, Species& ions, double Te, double Ti, double minDensity );
		
	//Computes the total charge density rho
	void computeRho( Species& electrons, Species& ions );
	
	//Compute total charge density with new evaluation of ion charge density only
	void computeRhoIonsOnly( Species& ions );
	
	//Compute total charge density with new evaluation of electron charge density only
	void computeRhoElectronsOnly( Species& electrons );
	
	//Computes quasi-neutral charge density for simulation without electron species
	void computeRhoQN( Species& ions, double nEl, double phiPlasma );		//currently not implemented correctly
	
	//Computes E-field as negative divergence of phi
	void computeEField();
	
	//Computes normalization of electric field
	VectorFieldCyl<double> normEfield();
	
	//Transformation between logical coordniates (i,j,k) and position vector
	VectorCyl<double> lcToPos( double i, double j, double k );
	VectorCyl<double> posToLc( VectorCyl<double>& pos );
	
	//Evaluates cell index of given position for given process
	int posToCell( VectorCyl<double>& pos );
	//Evaluates volume of cell with given cellIndex
	double cellVolume( int cellIndex );
	
	//Writes fields to files
	void writeToFileRho( int t );
	void writeToFileElDen( int t );
	void writeToFileIonDen( int t );
	void writeToFilePhi( int t );
	void writeToFileNormCurrDen( int t );
	void writeToFileCurrDen( int t );
	void writeToFileNGas( int t );
	void writeToFileEfield( int t );
	void writeToFileBfield( int t );
	void writeDynEfieldThR( int t, int k );
	void writeBfieldThR( int t, int k );
	void writeDynEfieldZR( int t, int j = 0 );
	void writeBfieldZR( int t, int j = 0, bool rePart = true );
	void writeDynEfieldThetaZR( int t, int j = 0, bool rePart = true );
	void writeToFileMatIndex( int nz0, int nzMax, double partial );
	void writeToFileMatIndexSurface( int nz0, int nzMax, double partial );
	void writeToFileMatIndexZR();
	void writeToFileMatIndexThetaR( double z );
	
	//Input of neutral Gas properties from file
	void inputNeutralGasDen( std::string filepath );
	void inputNeutralGasDen2D( std::string filepath );
	void inputNeutralGasVel( std::string filepath );
	void inputNeutralGasVel2D( std::string filepath );
	
};


//template member functions definitions
template <typename T>
VectorFieldCyl<T> Domain::curlE( VectorFieldCyl<T>& vField )
{
	//get number of grid points for fields
	int nrf{ static_cast<int>( vField.getnr() ) };
	int nthf{ static_cast<int>( vField.getnth() ) };
	int nzf{ static_cast<int>( vField.getnz() ) };
	
	//initialize return vector field
	VectorFieldCyl<T> reField( nrf, nthf, nzf, getpart() );
	
	if ( nthf == 1 )
	{
		reField = curlE2D( vField );
	}
	else
	{
		reField = curlE3D( vField );
	}
	
	return reField;
	
}

template <typename T>
VectorFieldCyl<T> Domain::curlE2D( VectorFieldCyl<T>& vField )
{
	//get number of grid points for fields
	int nrf{ static_cast<int>( vField.getnr() ) };
	int nzf{ static_cast<int>( vField.getnz() ) };
	
	int numZmin{ nZmin };
    int numZmax{ nZmax };
    if ( w_rank == w_size - 1 )
    {
        numZmax = nzf - 1;
    }
    
	//initialize return vector field
	VectorFieldCyl<T> reField( nrf, 1, nzf, getpart() );
	
	//bulk case
	for ( int i{ 0 }; i < nrf - 1; ++i )
	{
		for ( int k{ numZmin }; k < numZmax; ++k )
		{
			if ( materialIndex(i,0,k) == 1 && materialIndex(i,0,k+1) == 1 )
			{
				//compute field components on staggered grid with finite difference method or stokes theorem
				reField.r(i,0,k) = ( vField.th(i,0,k+1) - vField.th(i,0,k) )/dz;
				reField.th(i,0,k) = ( vField.r(i,0,k+1) - vField.r(i,0,k) )/dz - ( vField.z(i+1,0,k) - vField.z(i,0,k) )/dr;
				reField.z(i,0,k) = ( vField.th(i+1,0,k) - vField.th(i,0,k) )/dr;
			}
		}
	}
	//set open boundaries
	//case i == nrf - 1
	for ( int k{ numZmin }; k < numZmax; ++k )
	{
		reField.r(nrf-1,0,k) = ( vField.th(nrf-1,0,k+1) - vField.th(nrf-1,0,k) )/dz;
		reField.th(nrf-1,0,k) = reField.th(nrf-2,0,k);
		reField.z(nrf-1,0,k) = reField.z(nrf-2,0,k);
	}
	//case k == nzf - 1
	if ( w_rank == w_size - 1 )
	{
		for ( int i{ 1 }; i < nrf - 1; ++i )
		{
			reField.r(i,0,nzf-1) = reField.r(i,0,nzf-2);
			reField.th(i,0,nzf-1) = reField.th(i,0,nzf-2);
			reField.z(i,0,nzf-1) = ( vField.th(i+1,0,nzf-1) - vField.th(i,0,nzf-1) )/dr;
		}
		//case i == 0 && k == nzf - 1
		reField.r(0,0,nzf-1) = ( reField.r(0,0,nzf-2) + reField.r(1,0,nzf-1) )/2.0;
		reField.th(0,0,nzf-1) = ( reField.th(0,0,nzf-2) + reField.th(1,0,nzf-1) )/2.0;
		reField.z(0,0,nzf-1) = ( reField.z(0,0,nzf-2) + reField.z(1,0,nzf-1) )/2.0;
		
		//case i == nrf - 1 && k == nzf - 1
		reField.r(nrf-1,0,nzf-1) = ( reField.r(nrf-1,0,nzf-2) + reField.r(nrf-2,0,nzf-1) )/2.0;
		reField.th(nrf-1,0,nzf-1) = ( reField.th(nrf-1,0,nzf-2) + reField.th(nrf-2,0,nzf-1) )/2.0;
		reField.z(nrf-1,0,nzf-1) = ( reField.z(nrf-1,0,nzf-2) + reField.z(nrf-2,0,nzf-1) )/2.0;
	}
	
	if ( w_rank == 0 )
	{
		//case i == nrf - 1 && k == 0
		reField.r(nrf-1,0,0) = ( reField.r(nrf-1,0,1) + reField.r(nrf-2,0,0) )/2.0;
		reField.th(nrf-1,0,0) = ( reField.th(nrf-1,0,1) + reField.th(nrf-2,0,0) )/2.0;
		reField.z(nrf-1,0,0) = ( reField.z(nrf-1,0,1) + reField.z(nrf-2,0,0) )/2.0;
	}

	return reField;
	
}

template <typename T>
VectorFieldCyl<T> Domain::curlE3D( VectorFieldCyl<T>& vField )
{
	//get number of grid points for fields
	int nrf{ static_cast<int>( vField.getnr() ) };
	int nthf{ static_cast<int>( vField.getnth() ) };
	int nzf{ static_cast<int>( vField.getnz() ) };
	
	int numZmin{ nZmin };
    int numZmax{ nZmax };
    if ( w_rank == w_size - 1 )
    {
        numZmax = nzf - 1;
    }
    
	//initialize return vector field
	VectorFieldCyl<T> reField( nrf, nthf, nzf, getpart() );
	
	//bulk case, i != 0
	for ( int i{ 1 }; i < nrf - 1; ++i )
	{
		for ( int j{ 0 }; j < nthf; ++j )
		{
			for ( int k{ numZmin }; k < numZmax; ++k )
			{
				if ( materialIndex(i,j,k) == 1 )
                {
				//compute field components on staggered grid with finite difference method or stokes theorem
				reField.r(i,j,k) = ( 1.0 /(i*dr) )*( vField.z(i,j+1,k) - vField.z(i,j,k) )/dth - ( vField.th(i,j,k+1) - vField.th(i,j,k) )/dz;
				reField.th(i,j,k) = ( vField.r(i,j,k+1) - vField.r(i,j,k) )/dz - ( vField.z(i+1,j,k) - vField.z(i,j,k) )/dr;
				reField.z(i,j,k) = ( vField.th(i+1,j,k) - vField.th(i,j,k) )/dr - ( 1.0 /( ( i + 0.5 )*dr ) )*( vField.r(i,j+1,k) - vField.r(i,j,k) )/dth;
				}
			}
		}
	}
	//bulk case, i == 0
	for ( int j{ 0 }; j < nthf; ++j )
	{
		for ( int k{ numZmin }; k < numZmax; ++k )
		{
			if ( materialIndex(0,j,k) == 1 )
            {
			//compute field components on staggered grid with stokes theorem
			int j1{ static_cast<int>( ( ( reField.a(0,j,k) + constants::pi/2 )/( 2*constants::pi ) )*( nthf/part ) + 0.5 ) }; 
			int j2{ static_cast<int>( ( ( reField.a(0,j,k) + 3*constants::pi/2 )/( 2*constants::pi ) )*( nthf/part ) + 0.5 ) };
			
			reField.r(0,j,k) = ( vField.z(1,j1,k) - vField.z(1,j2,k) )/( 2*dr ) + ( vField.r(0,j2,k+1) - vField.r(0,j1,k+1) + vField.r(0,j1,k) - vField.r(0,j2,k) )/( 2*dz );
			reField.th(0,j,k) = ( vField.r(0,j,k+1) - vField.r(0,j,k) )/dz - ( vField.z(1,j,k) - vField.z(0,j,k) )/dr;
			reField.z(0,j,k) = ( vField.th(1,j,k) - vField.th(0,j,k) )/dr - ( 1.0 /( ( 0.5 )*dr ) )*( vField.r(0,j+1,k) - vField.r(0,j,k) )/dth;
			}
		}
	}
	//set open boundaries
	//case k == nzf - 1
	if ( w_rank == w_size - 1 )
	{
		for ( int i{ 0 }; i < nrf - 1; ++i )
		{
			for ( int j{ 0 }; j < nthf; ++j )
			{
				reField.r(i,j,nzf-1) = reField.r(i,j,nzf-2);
				reField.th(i,j,nzf-1) = reField.th(i,j,nzf-2);
				reField.z(i,j,nzf-1) = ( vField.th(i+1,j,nzf-1) - vField.th(i,j,nzf-1) )/dr - ( 1.0 /( ( i + 0.5 )*dr ) )*( vField.r(i,j+1,nzf-1) - vField.r(i,j,nzf-1) )/dth;
			}
		}
	}
	//case i == nrf - 1
	for ( int j{ 0 }; j < nthf; ++j )
	{
		for ( int k{ numZmin }; k < numZmax; ++k )
		{
			reField.r(nrf-1,j,k) = ( 1.0 /( (nrf-1)*dr) )*( vField.z(nrf-1,j+1,k) - vField.z(nrf-1,j,k) )/dth - ( vField.th(nrf-1,j,k+1) - vField.th(nrf-1,j,k) )/dz;
			reField.th(nrf-1,j,k) = reField.th(nrf-2,j,k);
			reField.z(nrf-1,j,k) = reField.z(nrf-2,j,k);
		}
	}
	//case i == nrf - 1 && k == nzf - 1
	if ( w_rank == w_size - 1 )
	{
		for ( int j{ 0 }; j < nthf; ++j )
		{
			reField.r(nrf-1,j,nzf-1) = ( reField.r(nrf-1,j,nzf-2) + reField.r(nrf-2,j,nzf-1) )/2.0;
			reField.th(nrf-1,j,nzf-1) = ( reField.th(nrf-1,j,nzf-2) + reField.th(nrf-2,j,nzf-1) )/2.0;
			reField.z(nrf-1,j,nzf-1) = ( reField.z(nrf-1,j,nzf-2) + reField.z(nrf-2,j,nzf-1) )/2.0;
		}
	}
	
	return reField;
}

template <typename T>
VectorFieldCyl<T> Domain::curlB( VectorFieldCyl<T>& vField )
{
	//get number of grid points for fields
	int nrf{ static_cast<int>( vField.getnr() ) };
	int nthf{ static_cast<int>( vField.getnth() ) };
	int nzf{ static_cast<int>( vField.getnz() ) };
	
	//initialize return vector field
	VectorFieldCyl<T> reField( nrf, nthf, nzf, getpart() );
	
	if ( nthf == 1 )
	{
		reField = curlB2D( vField );
	}
	else
	{
		reField = curlB3D( vField );
	}
	
	return reField;
}

template <typename T>
VectorFieldCyl<T> Domain::curlB2D( VectorFieldCyl<T>& vField )
{
	//get number of grid points for fields
	int nrf{ static_cast<int>( vField.getnr() ) };
	int nzf{ static_cast<int>( vField.getnz() ) };
	
	int numZmin{ nZmin };
    int numZmax{ nZmax };
	if ( w_rank == 0 )
    {
    	numZmin = 1;
    }
    
	//initialize return vector field
	VectorFieldCyl<T> reField( nrf, 1, nzf, getpart() );
	
	//bulk case, i != 0
	for ( int i{ 1 }; i < nrf; ++i )
	{
		for ( int k{ numZmin }; k < numZmax; ++k )
		{
			//compute field components on staggered grid with finite difference method or stokes theorem
			reField.r(i,0,k) = ( vField.th(i,0,k) - vField.th(i,0,k-1) )/dz;
			reField.th(i,0,k) = ( vField.r(i,0,k) - vField.r(i,0,k-1) )/dz - ( vField.z(i,0,k) - vField.z(i-1,0,k) )/dr;
			reField.z(i,0,k) = ( vField.th(i,0,k) - vField.th(i-1,0,k) )/dr;
		}
	}
	//set open boundaries
	//case i == 0
	for ( int k{ numZmin }; k < numZmax; ++k )
	{		
		reField.r(0,0,k) = ( vField.th(0,0,k) - vField.th(0,0,k-1) )/dz;
		reField.th(0,0,k) = reField.th(1,0,k);
		reField.z(0,0,k) = reField.z(1,0,k);
	}
	//case k == 0
	if ( w_rank == 0 )
	{
		for ( int i{ 1 }; i < nrf; ++i )
		{
			reField.r(i,0,0) = reField.r(i,0,1);
			reField.th(i,0,0) = reField.th(i,0,1);
			reField.z(i,0,0) = ( vField.th(i,0,0) - vField.th(i-1,0,0) )/dr;
		}
		//case i == 0 && k == 0
		reField.r(0,0,0) = ( reField.r(0,0,1) + reField.r(1,0,0) )/2.0;
		reField.th(0,0,0) = ( reField.th(0,0,1) + reField.th(1,0,0) )/2.0;
		reField.z(0,0,0) = ( reField.z(0,0,1) + reField.z(1,0,0) )/2.0;
		
		//case i == nrf - 1 && k == 0
		reField.r(nrf-1,0,0) = ( reField.r(nrf-1,0,1) + reField.r(nrf-2,0,0) )/2.0;
		reField.th(nrf-1,0,0) = ( reField.th(nrf-1,0,1) + reField.th(nrf-2,0,0) )/2.0;
		reField.z(nrf-1,0,0) = ( reField.z(nrf-1,0,1) + reField.z(nrf-2,0,0) )/2.0;
	}
	
	if ( w_rank == w_size - 1 )
	{
		//case i == 0 && k == nzf - 1
		reField.r(0,0,nzf-1) = ( reField.r(0,0,nzf-2) + reField.r(1,0,nzf-1) )/2.0;
		reField.th(0,0,nzf-1) = ( reField.th(0,0,nzf-2) + reField.th(1,0,nzf-1) )/2.0;
		reField.z(0,0,nzf-1) = ( reField.z(0,0,nzf-2) + reField.z(1,0,nzf-1) )/2.0;
	}
	
	return reField;
	
}

template <typename T>
VectorFieldCyl<T> Domain::curlB3D( VectorFieldCyl<T>& vField )
{
	//get number of grid points for fields
	int nrf{ static_cast<int>( vField.getnr() ) };
	int nthf{ static_cast<int>( vField.getnth() ) };
	int nzf{ static_cast<int>( vField.getnz() ) };
	
	int numZmin{ nZmin };
    int numZmax{ nZmax };
	if ( w_rank == 0 )
    {
    	numZmin = 1;
    }
    
	//initialize return vector field
	VectorFieldCyl<T> reField( nrf, nthf, nzf, getpart() );
	
	//bulk case, i != 0
	for ( int i{ 1 }; i < nrf; ++i )
	{
		for ( int j{ 0 }; j < nthf; ++j )
		{
			//for ( int k{ 1 }; k < nz; ++k )
			for ( int k{ numZmin }; k < numZmax; ++k )
			{
				//compute field components on staggered grid with finite difference method or stokes theorem
				reField.r(i,j,k) = ( 1.0 /( ( i + 0.5 )*dr) )*( vField.z(i,j,k) - vField.z(i,j-1,k) )/dth - ( vField.th(i,j,k) - vField.th(i,j,k-1) )/dz;
				reField.th(i,j,k) = ( vField.r(i,j,k) - vField.r(i,j,k-1) )/dz - ( vField.z(i,j,k) - vField.z(i-1,j,k) )/dr;
				reField.z(i,j,k) = ( vField.th(i,j,k) - vField.th(i-1,j,k) )/dr - ( 1.0 /( ( i )*dr ) )*( vField.r(i,j,k) - vField.r(i,j-1,k) )/dth;
			}
		}
	}
	//bulk case, i == 0
	for ( int j{ 0 }; j < nthf; ++j )
	{
		//for ( int k{ 1 }; k < nz; ++k )
		for ( int k{ numZmin }; k < numZmax; ++k )
		{
			//compute field components on staggered grid with stokes theorem
			int j2{ static_cast<int>( ( ( reField.a(0,j,k) + constants::pi )/( 2*constants::pi ) )*( nthf/part ) + 0.5 ) };
			
			reField.r(0,j,k) = ( 1.0 /( ( 0.5 )*dr) )*( vField.z(0,j,k) - vField.z(0,j-1,k) )/dth - ( vField.th(0,j,k) - vField.th(0,j,k-1) )/dz;
			reField.th(0,j,k) = ( vField.r(0,j,k) - vField.r(0,j,k-1) )/dz  - ( vField.z(0,j,k) - vField.z(0,j2,k) )/dr;
			
			for ( int jP{ 0 }; jP < nthf; ++jP )
			{
				reField.z(0,j,k) += 4.0*vField.th(0,jP,k)/( dr*nthf*part );		//area of regular n-polygon A=nthf*(dr*dth/2)*(dr/2)/2 (for partial simulation *part)
			}
		}
	}
	//set open boundaries
	//case k == 0
	if ( w_rank == 0 )
	{
		for ( int i{ 1 }; i < nrf; ++i )
		{
			for ( int j{ 0 }; j < nthf; ++j )
			{
				reField.r(i,j,0) = reField.r(i,j,1);
				reField.th(i,j,0) = reField.th(i,j,1);
				reField.z(i,j,0) = ( vField.th(i,j,0) - vField.th(i-1,j,0) )/dr - ( 1.0 /( ( i )*dr ) )*( vField.r(i,j,0) - vField.r(i,j-1,0) )/dth;
			}
		}
		for ( int j{ 0 }; j < nthf; ++j )
		{
			reField.r(0,j,0) = reField.r(0,j,1);
			reField.th(0,j,0) = reField.th(0,j,1);
			
			for ( int jP{ 0 }; jP < nthf; ++jP )
			{
				reField.z(0,j,0) += 4.0*vField.th(0,jP,0)/( dr*nthf*part );		//area of regular n-polygon A=nthf*(dr*dth/2)*(dr/2)/2
			}
		}
	}
	
	return reField;
}


#endif // DOMAIN_H_INCLUDED
