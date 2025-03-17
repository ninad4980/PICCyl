//==========================================================================================
//
// Copyright (C) 2021 Dr. Ninad Joshi <email@ninadjoshi.de>  
//
// All Rights Reserved
//
//==========================================================================================

#include "collisions.h"
#include "domain.h"
#include "species.h"
#include "vectorcyl.h"
#include "constants.h"
#include "parameters.h"
#include <cmath>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <complex>
#include <mpi.h>

using Complex = std::complex<double>;

//Constructor for MCC interaction class, sets correct mass for neutral gas type if available
MCC::MCC ( Species& sp, Domain& domaind, std::string gasT ): species{ sp }, domain{ domaind }, gas{ gasT }
{
	//checks if chosen gas type is available and assigns atom mass accordingly
	if ( gasT == "krypton" )
	{
		M = constants::m_kr;
	}
	else if ( gasT == "neon" )
	{
		M = constants::m_ne;
	}
	else if ( gasT == "xenon" )
	{
		M = constants::m_xe;
	}
	else
	{
		std::cerr << "Invalid input for gas type! Please use neon, krypton or xenon." << std::endl;
	}
}

//Constructor for ElectronNeutralCollision interaction class, sets correct values for available neutral gas types
ElectronNeutralCollision::ElectronNeutralCollision ( Species& electron, Species& ion, Domain& domaind, std::string gasT ): MCC{ electron, domaind, gasT }, ions{ ion }
{
	//checks selected gas type and assigns corresponding values for ionization energy, electron-energy distribution factor and maximal collision cross section
	if ( gasT == "krypton" )
	{
		ionizationEn = electronNeutralIonization::Kr::I;
		B = electronNeutralIonization::Kr::S;
		sigmaMax = 3.5e-20;
	}
	else if ( gasT == "neon" )
	{
		ionizationEn = electronNeutralIonization::Ne::I;
		B = electronNeutralIonization::Ne::S;
		sigmaMax = 6.74e-21;
	}
	else if ( gasT == "xenon" )
	{
		ionizationEn = electronNeutralIonization::Xe::I;
		excitationEn = electronNeutralExcitation::Xe::E;
		B = electronNeutralIonization::Xe::S;
		sigmaMax = 2.9e-19;
	}
}

void MCC::computeVMax ()
{
	vMax = 0;
	//loops through all particles
	for ( long long p{ 0 }; p < static_cast<long long>( species.size() ); ++p )
	{
		//computes particle velocity norm and assigns it to vMax, if it is larger than current vMax
		Particle& part{ species.at(p) };
		part.vel.realNorm();
		if ( part.vel.norm > vMax ) { vMax = part.vel.norm; }
	}
	
	//check if there are multiple processes
	if ( domain.getsize() > 1 )
	{
		//get maximum velocity across all processes
		double vMaxSingle{ vMax };
		MPI_Allreduce( &vMaxSingle, &vMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
		vMax = vMaxSingle;
	}
}

int MCC::computeNMax ( double dt )
{
	using std::exp;
	long long N{ static_cast<long long>( species.size() ) };
	long long NTot{ N };
	
	//check if there are multiple processes
	if ( domain.getsize() > 1 )
	{
		//compute total number of particles
		MPI_Allreduce( &N, &NTot, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD );
	} 
	//compute maximum collision frequency
	nuMax = nNeutralMax*sigmaMax*vMax;
	//compute maximal number of collisions
	long long NMax{ static_cast<long long>( NTot*( 1 - exp( -nuMax*dt ) ) + 0.5 ) };
	//check if there are multiple processes
	if ( domain.getsize() > 1 )
	{
		//compute maximal number of colliding particles for single process
		long long NMaxSingle{ static_cast<long long>( NMax*N*1.0/NTot + 0.5 ) };
		NMax = NMaxSingle;
	}
	//set minimal value of NMax to 1
	if ( NMax == 0 ) { NMax = 1; }
	//check if NMax is larger than number of particles and adjust NMax accordingly
	if ( NMax > N ) { NMax = N; }
	return NMax;
}

void ElectronNeutralCollision::apply ( double dt )
{
	domain.timer().startTimer( "collision module" );
	//seed prng for starting values
	std::mt19937 mersenne{ static_cast<std::mt19937::result_type>(std::time(nullptr)) };

	//evaluate maximal number of collisions taking place Nmax
	computeVMax();
	int Nmax{ computeNMax( dt ) };
	//std::cout << "process " << domain.getrank() << " Nmax is " << Nmax << '\n';
	//std::cout << "process " << domain.getrank() << " particleNum is " << species.size() << '\n';
	
	std::uniform_int_distribution<int> rN1( 0, species.size() - 1 );
	std::uniform_real_distribution<double> rN2( 0.0, 1.0 );
	
	//loop through Nmax random particles stored in species particles vector
	for ( int p{ 0 }; p < Nmax; ++p )
	{
		//select a radnom particle form species particles vector
		Particle& part{ species.at( rN1(mersenne) ) };
		
		//get neutral gas density at particle position
		//double nNeutral{ domain.nGasDen().gather( part.pos, domain.getdr(), domain.getdth(), domain.getdz(), domain.getnth() ) };
		
		//evaluate cross sections for particle
		evalSigmaSc( part );
		evalSigmaEx( part );
		evalSigmaIo( part );
		//evaluate collision frequencies for particle
		double nuSc{ nNeutral*sigmaSc*part.vel.norm };
		double nuEx{ nNeutral*sigmaEx*part.vel.norm };
		double nuIo{ nNeutral*sigmaIo*part.vel.norm };
		//get random number between 0 and 1 for "null collision" procedure
		double ctProb{ rN2(mersenne) };
		
		//test if scattering, excitation, ionization or null collision takes place and apply respective function
		if ( ctProb < (nuSc / nuMax) )
		{
			elasticScatter( part, nuSc );
			++eNElasticCollisions;
		}
		else if ( ctProb < ( (nuSc + nuEx) / nuMax ) )
		{
			excitation( part, nuEx );
			++eNExcitations;
		}
		else if ( ctProb < ( (nuSc + nuEx + nuIo) / nuMax ) )
		{
			ionization( part, nuIo );
			++eNIonizations;
		}
		else
		{
			++eNNullCollisions;
		}
		
	}
	
	domain.timer().stopTimer( "collision module" );
}

void ElectronNeutralCollision::evalSigmaSc ( Particle& part )
{
	using std::pow;
	using std::log;
	using std::exp;
	
	//check if selected gas type is xenon
	if ( gas == "xenon" )
	{
		//compute particle's kinetic energy (in eV)
		double kinEn{ ( constants::m_e*part.vel.norm*part.vel.norm / 2 ) / constants::q_e };
		double x{ kinEn };
		
		//check which fitting formula has to be applied depending on energy
		if ( x > 16.5 )
		{
			//initialize parameter values (stored in parametrs.h)
			double F{ electronNeutralElasticScatter::Xe::F };
			double G{ electronNeutralElasticScatter::Xe::G };
			
			//compute cross section with appropriate fitting formula
			sigmaSc = F*pow( x, G );
		}
		else if ( x > 0.63 )
		{
			//initialize parameter values (stored in parametrs.h)
			double A{ electronNeutralElasticScatter::Xe::A };
			double B{ electronNeutralElasticScatter::Xe::B };
			double C{ electronNeutralElasticScatter::Xe::C };
			double D{ electronNeutralElasticScatter::Xe::D };
			double E{ electronNeutralElasticScatter::Xe::E };
			
			//compute cross section with appropriate fitting formula
			sigmaSc = A*exp( B*pow( log( x ) - C, 2 ) - E ) + D;
		}
		else if ( x > 0 )
		{
			//initialize parameter values (stored in parametrs.h)
			double a{ electronNeutralElasticScatter::Xe::a };
			double b{ electronNeutralElasticScatter::Xe::b };
			double c{ electronNeutralElasticScatter::Xe::c };
			double d{ electronNeutralElasticScatter::Xe::d };
			
			//compute cross section with appropriate fitting formula
			sigmaSc = a*exp( b*pow( x - d, c ) );
		}
	}
	else
	{
		std::cerr << "Evaluation of elastic scatter cross section implemented for Xenon only!" << std::endl;
	}
	
	//for negative computed cross section set sigmaSc to zero
	if ( sigmaSc < 0 ) { sigmaSc = 0; }
}

void ElectronNeutralCollision::evalSigmaEx ( Particle& part )
{
	using std::log;
	
	double A{};
	double B{};
	double C{};
	double D{};
	double E{};
	
	//compute kinetic energy of particle (in eV)
	double kinEn{ ( constants::m_e*part.vel.norm*part.vel.norm / 2 ) / constants::q_e };
	
	//check if gas type is xenon
	if ( gas == "xenon" )
	{
		//assign parameter values for xenon (stored in parametrs.h)
		A = electronNeutralExcitation::Xe::A;
		B = electronNeutralExcitation::Xe::B;
		C = electronNeutralExcitation::Xe::C;
		D = electronNeutralExcitation::Xe::D;
		E = electronNeutralExcitation::Xe::E;
		
		double x{ kinEn / E };
		
		//check if kinetic energy is larger than minimum excitation energy
		if ( x > 1 )
		{
			//compute cross section with fitting formula
			sigmaEx = ( 1 / x )*( A*log( x ) + ( B*log( x ) + C*( x - 1 ) ) / ( x + D ) );
		}
		else
		{
			//in case kinEn < E set cross section to zero
			sigmaEx = 0.0;
		}	
	}
	else
	{
		std::cerr << "Evaluation of excitation cross section implemented for Xenon only!" << std::endl;
	}
	
	//for negative computed cross section set sigmaEx to zero
	if ( sigmaEx < 0 ) { sigmaEx = 0; }
}

void ElectronNeutralCollision::evalSigmaIo ( Particle& part )
{
	using std::log;
	double A{}; 
	double B{};
	double C{};
	double D{};
	double I{};
	
	//compute kinetic energy of particle (in eV)
	double kinEn{ ( constants::m_e*part.vel.norm*part.vel.norm / 2 ) / constants::q_e };
	
	//compute ionization cross section for given gas type and energy kinEn
	if ( gas == "neon" )
	{
		//assign parameter values for nenon (stored in parametrs.h)
		A = electronNeutralIonization::Ne::A;
		B = electronNeutralIonization::Ne::B;
		C = electronNeutralIonization::Ne::C;
		D = electronNeutralIonization::Ne::D;
		I = electronNeutralIonization::Ne::I;	
		
		double x{ kinEn / I };
		//check if kinEn is larger than ionization energy
		if ( x > 1 )
		{
			//compute cross section with fitting formula
			sigmaIo = ( 1 / x )*( A*log( x ) + ( B*log( x ) + C*( x - 1 ) ) )/( x + D )*(1e-4);
		}
		else
		{
			//if kinEn < I set cross section to zero
			sigmaIo = 0.0;
		}
	}
	else if ( gas == "xenon" )
	{
		//assign parameter values for xenon (stored in parametrs.h)
		A = electronNeutralIonization::Xe::A;
		B = electronNeutralIonization::Xe::B;
		C = electronNeutralIonization::Xe::C;
		D = electronNeutralIonization::Xe::D;
		I = electronNeutralIonization::Xe::I;
		
		double x{ kinEn / I };
		//check if kinEn is larger than ionization energy
		if ( x > 1 )
		{
			//compute cross section with fitting formula
			sigmaIo = ( 1 / x )*( A*log( x ) + ( B*log( x ) + C*( x - 1 ) ) )/( x + D );
		}
		else
		{
			//if kinEn < I set cross section to zero
			sigmaIo = 0.0;
		}
	}
	else if ( gas == "krypton" )
	{
		//assign parameter values for krypton (stored in parametrs.h)
		double a{ electronNeutralIonization::Kr::a };
		double b{ electronNeutralIonization::Kr::b };
		double c{ electronNeutralIonization::Kr::c };
		I = electronNeutralIonization::Kr::I;
	
		double x{ kinEn / I };
		//check if kinEn is larger than ionization energy
		if ( x > 1 )
		{
			//compute cross section with fitting formula
			sigmaIo = ( 1 / x )*( a*log( x ) + b + c*( 1 / x ) )*(1e-4);
		}
		else
		{
			//if kinEn < I set cross section to zero
			sigmaIo = 0.0;
		}
	}
	
	//for negative computed cross section set sigmaIo to zero (cross section can not be smaller than zero)
	if ( sigmaIo < 0 ) { sigmaIo = 0.0; }
}

void ElectronNeutralCollision::elasticScatter ( Particle& part, double nuSc, double energyLoss )
{
	using std::cos;
	using std::sin;
	using std::acos;
	using std::pow;
	using std::sqrt;
	using std::arg;
	using std::tan;
	
	//seed prng for starting values
	std::mt19937 mersenne{ static_cast<std::mt19937::result_type>(std::time(nullptr)) };
	std::uniform_real_distribution<double> rN( 0.0, 1.0 );
	
	//compute normalized vector in particle velocity direction
	VectorCyl<double> vInc{ part.vel.real().normalize() };
	VectorCyl<double> eR{ 1, 0, 0, part.pos.a };
	
	//kinetic energy of particle in eV
	double kinEn{ ( constants::m_e*part.vel.norm*part.vel.norm / 2 ) / constants::q_e - energyLoss };
	
	//compute scattering angles
	double chi{ acos( ( 2 + kinEn - 2*pow( 1 + kinEn, rN(mersenne) ) ) / kinEn ) };
	double phi{ 2*constants::pi*rN(mersenne) };
	double alpha{ acos( vInc*eR ) };
	
	//velocity direction of scattered electron
	VectorCyl<double> vScatter{ vInc*cos( chi ) + ( vInc % eR )*( sin( chi )*sin( phi )/sin( alpha ) ) + (vInc % ( eR % vInc ))*( sin( chi )*cos( phi )/sin( alpha ) ) };
	
	//kinetic energy loss in collision
	double deltaEn{ ( 2*constants::m_e / M )*( 1 - cos( chi ) )*kinEn };
	
	//new velocity norm
	double vNorm{ sqrt( ( kinEn - deltaEn )*2 / constants::m_e ) };
	
	//compute imaginary part of new complex velocity
	VectorCyl<double> phase{ arg( part.vel.r ), arg( part.vel.theta ), arg( part.vel.z ), part.vel.a };
	VectorCyl<Complex> imag{ { 0, tan( phase.r )*vScatter.r }, { 0, tan( phase.theta )*vScatter.theta }, { 0, tan( phase.z )*vScatter.z }, phase.a };
	 
	//velocity after scattering
	part.vel = (vScatter + imag)*vNorm;	
	
}

void ElectronNeutralCollision::excitation ( Particle& part, double nuEx )
{
	//scattering after excitation energy substraction
	elasticScatter( part, nuEx, excitationEn );
}

void ElectronNeutralCollision::ionization ( Particle& part, double nuIo )
{
	using std::tan;
	using std::atan2;
	using std::sqrt;
	
	//seed prng for starting values
	std::mt19937 mersenne{ static_cast<std::mt19937::result_type>(std::time(nullptr)) };
	
	//random number distribution for kinetic energy distribution between electrons
	std::uniform_real_distribution<double> rN( 0.0, 1.0 );
	
	//kinetic energy of incoming electron
	double kinEnInc{ ( constants::m_e*part.vel.norm*part.vel.norm / 2 ) / constants::q_e };
	
	//compute energy of new electron
	double kinEnNew{ B*tan( rN(mersenne)*atan2( kinEnInc - ionizationEn, 2*B ) ) };
	//velocity norm of new electron
	double vNewNorm{ sqrt( kinEnNew*2*constants::q_e / constants::m_e ) };
	//velocity norm of incident electron after energy loss 
	double vNorm{ sqrt( ( kinEnInc - ionizationEn - kinEnNew )*2*constants::q_e / constants::m_e ) };
	
	//update particle of incident electron after energy loss
	part.vel = part.vel.normalize()*vNorm;
	
	//create new electron particle
	Particle partNew{ part.pos, part.vel.normalize()*vNewNorm, -1, true };
	
	//scattering of both electrons
	elasticScatter( part, nuIo );
	elasticScatter( partNew, nuIo );
	
	//add new electron and ion to respective species class objects
	species.addParticle( partNew );
	ions.addMaxwellianParticle( part.pos, T );
	
}

void ElectronNeutralCollision::writeLogfile( int t )
{
	//write logfile for collisions of each process
	//set up mpi-file
	std::string filepath{ "logfiles/ElMCCLogfilet" + std::to_string( t ) + ".dat" };
	MPI_Status status;
	MPI_File fh;
	MPI_File_open( MPI_COMM_WORLD, filepath.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );
	
	//set up output buffers
	char buf[410];
	char buf4[35];
	char buf5[44];
	char buf6[44];
	char buf7[44];
	char buf8[44];
	
	//write number of collision events to buffer
	snprintf( buf, 25, "Collisions process %i \n", domain.getrank() );
	snprintf( buf4, 34, "MCC Electron-Neutral-Collisions \n" );
	strcat( buf, buf4 );
	snprintf( buf5, 43, "Null Collisions: %f \n", eNNullCollisions );
	strcat( buf, buf5 );
	snprintf( buf6, 43, "Elastic Scatter Collisions: %f \n", eNElasticCollisions );
	strcat( buf, buf6 );
	snprintf( buf7, 43, "Excitation Collisions: %f \n", eNExcitations );
	strcat( buf, buf7 );
	snprintf( buf8, 43, "Ionization Collisions: %f \n", eNIonizations );
	strcat( buf, buf8 );
	
	//write output buffer of every process to file
	MPI_File_write_ordered( fh, buf, strlen(buf), MPI_CHAR, &status );
	
	MPI_File_close( &fh );
	
	//write logfile for total number of collisions
	
	double eNNullCollisionsTot{ eNNullCollisions };
	double eNElasticCollisionsTot{ eNElasticCollisions };
	double eNExcitationsTot{ eNExcitations };
	double eNIonizationsTot{ eNIonizations };
	
	//for multiple processes evaluate total number of events for all collisions types
	if ( domain.getsize() != 1 )
	{
		MPI_Reduce( &eNNullCollisions, &eNNullCollisionsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &eNElasticCollisions, &eNElasticCollisionsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &eNExcitations, &eNExcitationsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &eNIonizations, &eNIonizationsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );	
	}
	
	//write total number of collisions to file
	if ( domain.getrank() == 0 )
	{
		std::ofstream lFile{ "logfiles/ElMCCLogfileTotalt" + std::to_string( t ) + ".dat" };
		if ( !lFile )
		{
			std::cerr << "Logfile electrons could not be opened!" << std::endl;
		}
		
		lFile << "Total MCC Electron-Neutral-Collisions\n";
		lFile << "Null Collisions: " << eNNullCollisionsTot << '\n';
		lFile << "Elastic Scatter Collisions: " << eNElasticCollisionsTot << '\n';
		lFile << "Excitation Collisions: " << eNExcitationsTot << '\n';
		lFile << "Ionization Collisions: " << eNIonizationsTot << '\n';
		
		lFile.close();
	}
}

void ElectronNeutralCollision::writeLog( std::string filepath, int access_mode )
{
	//write logfile for collisions of each process
	//set up mpi-file
	MPI_Status status;
	MPI_File fh;
	MPI_File_open( MPI_COMM_WORLD, filepath.c_str(), access_mode, MPI_INFO_NULL, &fh );
	
	//set up output buffers
	char buf[235];
	char buf4[35];
	char buf5[50];
	char buf6[50];
	char buf7[50];
	char buf8[50];
	
	//write number of collision events to buffer
	snprintf( buf, 25, "Collisions process %i\n", domain.getrank() );
	snprintf( buf4, 34, "MCC Electron-Neutral-Collisions\n" );
	strcat( buf, buf4 );
	snprintf( buf5, 49, "Null Collisions: %.6g\n", eNNullCollisions );
	strcat( buf, buf5 );
	snprintf( buf6, 49, "Elastic Scatter Collisions: %.6g\n", eNElasticCollisions );
	strcat( buf, buf6 );
	snprintf( buf7, 49, "Excitation Collisions: %.6g\n", eNExcitations );
	strcat( buf, buf7 );
	snprintf( buf8, 49, "Ionization Collisions: %.6g\n\n", eNIonizations );
	strcat( buf, buf8 );
	
	//write output buffer of every process to file
	MPI_File_write_ordered( fh, buf, strlen(buf), MPI_CHAR, &status );
	
	MPI_File_close( &fh );

}

void ElectronNeutralCollision::writeLogTot( std::string filepath, int access_mode )
{
	//write logfile for collisions of each process
	//set up mpi-file
	MPI_Status status;
	MPI_File fh;
	MPI_File_open( MPI_COMM_WORLD, filepath.c_str(), access_mode, MPI_INFO_NULL, &fh );
	
	
	//write log for total number of collisions
	double eNNullCollisionsTot{ eNNullCollisions };
	double eNElasticCollisionsTot{ eNElasticCollisions };
	double eNExcitationsTot{ eNExcitations };
	double eNIonizationsTot{ eNIonizations };
	
	//for multiple processes evaluate total number of events for all collisions types
	if ( domain.getsize() != 1 )
	{
		MPI_Reduce( &eNNullCollisions, &eNNullCollisionsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &eNElasticCollisions, &eNElasticCollisionsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &eNExcitations, &eNExcitationsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &eNIonizations, &eNIonizationsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );	
	}
	
	//write total number of collisions to file
	if ( domain.getrank() == 0 )
	{
		//set up output buffers
		char buf[230];
		
		//write number of collision events to buffer
		snprintf( buf, 39 , "Total MCC Electron-Neutral-Collisions\n" );
		char buf1[40];
		snprintf( buf1, 39, "Null Collisions: %.6g\n", eNNullCollisionsTot );
		strcat( buf, buf1 );
		char buf2[50];
		snprintf( buf2, 49, "Elastic Scatter Collisions: %.6g\n", eNElasticCollisionsTot );
		strcat( buf, buf2 );
		char buf3[50];
		snprintf( buf3, 49, "Excitation Collisions: %.6g\n", eNExcitationsTot );
		strcat( buf, buf3 );
		char buf4[50];
		snprintf( buf4, 49, "Ionization Collisions: %.6g\n\n", eNIonizationsTot );
		strcat( buf, buf4 );
		
		MPI_File_write( fh, buf, strlen(buf), MPI_CHAR, &status );

	}
	
	MPI_File_close( &fh );

}

void IonNeutralCollision::apply ( double dt )
{
	domain.timer().startTimer( "collision module" );
	using std::sqrt;
		
	//seed prng for starting values
	std::mt19937 mersenne{ static_cast<std::mt19937::result_type>(std::time(nullptr)) };
	
	//standard deviation of velocity distributions of neutral atom
	double stDev{ sqrt( constants::k_B*T / M ) };
	
	//set up normal Distributions for every velocity component (Maxwellian velocity distribution)
	std::normal_distribution<double> v1( 0.0, stDev );
	std::normal_distribution<double> v2( 0.0, stDev );
	std::normal_distribution<double> v3( 0.0, stDev );
	
	std::uniform_real_distribution<double> rN2( 0.0, 1.0 );
	
	//loop through all particles in species particles vector
	for ( long long p{ 0 }; p < static_cast<long long>( species.size() ); ++p )
	{
		Particle& part{ species.at( p ) };
		
		//get neutral gas density at particle position
		//double nNeutral{ domain.nGasDen().gather( part.pos, domain.getdr(), domain.getdth(), domain.getdz(), domain.getnth() ) };
		
		//create random velocity for neutral gas atom from Maxwellian distribution
		VectorCyl<Complex> vNeutral{ v1(mersenne), v2(mersenne), v3(mersenne), 0.0 };
		vNeutral.angleShift( part.pos.a );
		
		//include drift velocity of neutral particles
		//vNeutral += nGasVel().gather( part.pos, domain.getdr(), domain.getdth(), domain.getnz(), domain.getnth() );
		
		//transform particle velocity to neutral atom reference frame
		part.vel -= vNeutral.complexConversion();
		
		//compute norm of particle velocity
		part.vel.realNorm();
	
		//evaluate cross sections for particle
		evalSigmaSc( part );
		evalSigmaChEx( part );
		//evaluate collision frequencies for particle
		double nuSc{ nNeutral*sigmaSc*part.vel.norm };
		double nuChEx{ nNeutral*sigmaChEx*part.vel.norm };
		
		//get random number between 0 and 1
		double ctProb1{ rN2(mersenne) };
		double ctProb2{ rN2(mersenne) };
		
		//check if collision occurs and apply respective function
		if ( ctProb1 < ( 1 - exp( -nuSc*dt ) ) )
		{
			elasticScatter( part );
			++iNElasticCollisions;
		}
		
		if ( ctProb2 < ( 1 - exp( -nuChEx*dt ) ) )
		{
			chargeExchange( part );
			++iNChargeExchanges;
		}
		
		//transform particle velocity back to laboratory frame
		part.vel += vNeutral.complexConversion();
		
	}	
	
	domain.timer().stopTimer( "collision module" );
}

void IonNeutralCollision::evalSigmaSc( Particle& part )
{
	using std::sqrt;
	using std::pow;
	
	double a{};
	double b{};
	
	//compute kinetic energy of particle
	double kinEn{ ( M*part.vel.norm*part.vel.norm / 2 ) / constants::q_e };

	//check if selected gas type is xenon
	if ( gas == "xenon" )
	{
		//assign parameter values for xenon (stored in parameters.h)
		a = ionNeutralElasticScatter::Xe::a;
		b = ionNeutralElasticScatter::Xe::b;
		
		//evaluate elastic scattering cross section with fitting formula
		sigmaSc = a*pow( kinEn, b );
	}
	else
	{
		std::cerr << "Evaluation of charge exchange cross section implemented for xenon only!" << std::endl;
	}
}

void IonNeutralCollision::evalSigmaChEx( Particle& part )
{
	using std::sqrt;
	
	//compute kinetic energy of particle
	double kinEn{ ( M*part.vel.norm*part.vel.norm / 2 ) / constants::q_e };
	
	//cross section charge exchange formula for neon from Plasmapic code
	if ( gas == "neon" )
	{
		sigmaChEx = 5.98e-20 / ( sqrt( kinEn ) + 1e-30 );
	}
	//cross section charge exchange formula for xenon from Plasmapic code
	else if ( gas == "xenon" )
	{
		sigmaChEx = 4.541e-20 / ( sqrt( kinEn ) + 1e-30 );
	}
	else
	{
		std::cerr << "Evaluation of charge exchange cross section implemented for neon and xenon only!" << std::endl;
	}
}

void IonNeutralCollision::elasticScatter( Particle& part )
{	
	using std::sqrt;
	using std::cos;
	using std::sin;
	using std::acos;
	using std::pow;
	using std::arg;
	using std::tan;
	
	//seed prng for starting values
	std::mt19937 mersenne{ static_cast<std::mt19937::result_type>(std::time(nullptr)) };	
	std::uniform_real_distribution<double> rN( 0.0, 1.0 );
	
	//compute normalized vector in particle velocity direction
	VectorCyl<double> vInc{ part.vel.real().normalize() };
	VectorCyl<double> eR{ 1, 0, 0, part.pos.a };
	
	//kinetic energy of particle in eV
	double kinEn{ ( constants::m_e*part.vel.norm*part.vel.norm / 2 ) / constants::q_e };
	
	//compute scattering angles
	double Theta{ acos( 1 - 2*rN(mersenne) ) };
	double phi{ 2*constants::pi*rN(mersenne) };
	double alpha{ acos( vInc*eR ) };
	
	//velocity direction of scattered electron
	VectorCyl<double> vScatter{ vInc*cos( Theta ) + ( vInc % eR )*( sin( Theta )*sin( phi )/sin( alpha ) ) + (vInc % ( eR % vInc ))*( sin( Theta )*cos( phi )/sin( alpha ) ) };
	
	//kinetic energy after scattering
	double kinEnSc{ kinEn*cos( Theta / 2 )*cos( Theta / 2) };
	
	//new velocity norm
	double vNorm{ sqrt( kinEnSc*2 / constants::m_e ) };
	
	//compute imaginary part of new complex velocity
	VectorCyl<double> phase{ arg( part.vel.r ), arg( part.vel.theta ), arg( part.vel.z ), part.vel.a };
	VectorCyl<Complex> imag{ { 0, tan( phase.r )*vScatter.r }, { 0, tan( phase.theta )*vScatter.theta }, { 0, tan( phase.z )*vScatter.z }, phase.a };
	 
	//velocity after scattering
	part.vel = (vScatter + imag)*vNorm;	
		
}

void IonNeutralCollision::chargeExchange( Particle& part )
{
	//in neutral atom reference frame the particle velocity of the new ion is zero
	part.vel.r = 0.0;
	part.vel.theta = 0.0;
	part.vel.z = 0.0;
}

void IonNeutralCollision::writeLogfile( int t )
{
	//write logfile for collisions of each process
	//set up mpi-file
	std::string filepath{ "logfiles/IoMCCLogfilet" + std::to_string( t ) + ".dat" };
	MPI_Status status;
	MPI_File fh;
	MPI_File_open( MPI_COMM_WORLD, filepath.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );
	
	//set up output buffers
	char buf[144];
	char buf1[30];
	char buf2[44];
	char buf3[44];

	//write number of collision events to buffer
	snprintf( buf, 25, "Collisions process %i \n", domain.getrank() );
	snprintf( buf1, 29, "MCC Ion-Neutral-Collisions \n" );
	strcat( buf, buf1 );
	snprintf( buf2, 43, "Elastic Scatter Collisions: %f \n", iNElasticCollisions );
	strcat( buf, buf2 );
	snprintf( buf3, 43, "Charge Exchange Collisions: %f \n", iNChargeExchanges );
	strcat( buf, buf3 );
	
	//write output buffer of every process to file
	MPI_File_write_ordered( fh, buf, strlen(buf), MPI_CHAR, &status );
	
	MPI_File_close( &fh );
	
	//write logfile for total number of collisions
	double iNElasticCollisionsTot{ iNElasticCollisions };
	double iNChargeExchangesTot{ iNChargeExchanges };

	//for multiple processes evaluate total number of events for all collisions types
	if ( domain.getsize() != 1 )
	{
		MPI_Reduce( &iNElasticCollisions, &iNElasticCollisionsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &iNChargeExchanges, &iNChargeExchangesTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	}

	//write total number of collisions to file
	if ( domain.getrank() == 0 )
	{
		std::ofstream lFile{ "logfiles/IoMCCLogfileTotalt" + std::to_string( t ) + ".dat" };
		if ( !lFile )
		{
			std::cerr << "Logfile ions could not be opened!" << std::endl;
		}
		
		lFile << "Total MCC Ion-Neutral-Collisions\n";
		lFile << "Elastic Scatter Collisions: " << iNElasticCollisionsTot << '\n';
		lFile << "Charge Exchange Collisions: " << iNChargeExchangesTot << '\n';
		
		lFile.close();
	}
}

void IonNeutralCollision::writeLog( std::string filepath, int access_mode )
{
	//write logfile for collisions of each process
	//set up mpi-file
	MPI_Status status;
	MPI_File fh;
	MPI_File_open( MPI_COMM_WORLD, filepath.c_str(), access_mode, MPI_INFO_NULL, &fh );
	
	//set up output buffers
	char buf[164];
	char buf1[30];
	char buf2[50];
	char buf3[50];

	//write number of collision events to buffer
	snprintf( buf, 25, "Collisions process %i\n", domain.getrank() );
	snprintf( buf1, 29, "MCC Ion-Neutral-Collisions\n" );
	strcat( buf, buf1 );
	snprintf( buf2, 49, "Elastic Scatter Collisions: %.6g\n", iNElasticCollisions );
	strcat( buf, buf2 );
	snprintf( buf3, 49, "Charge Exchange Collisions: %.6g\n\n", iNChargeExchanges );
	strcat( buf, buf3 );
	
	//write output buffer of every process to file
	MPI_File_write_ordered( fh, buf, strlen(buf), MPI_CHAR, &status );
	
	MPI_File_close( &fh );
	
}

void IonNeutralCollision::writeLogTot( std::string filepath, int access_mode )
{
	//write logfile for collisions of each process
	//set up mpi-file
	MPI_Status status;
	MPI_File fh;
	MPI_File_open( MPI_COMM_WORLD, filepath.c_str(), access_mode, MPI_INFO_NULL, &fh );
	
	//write log for total number of collisions
	double iNElasticCollisionsTot{ iNElasticCollisions };
	double iNChargeExchangesTot{ iNChargeExchanges };

	//for multiple processes evaluate total number of events for all collisions types
	if ( domain.getsize() != 1 )
	{
		MPI_Reduce( &iNElasticCollisions, &iNElasticCollisionsTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &iNChargeExchanges, &iNChargeExchangesTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	}

	//write total number of collisions to file
	if ( domain.getrank() == 0 )
	{
		char buf[135];
		
		snprintf( buf, 34 , "Total MCC Ion-Neutral-Collisions\n" );
		char buf1[50];
		snprintf( buf1, 49, "Elastic Scatter Collisions: %.6g\n", iNElasticCollisionsTot );
		strcat( buf, buf1 );
		char buf2[50];
		snprintf( buf2, 49, "Charge Exchange Collisions: %.6g\n\n", iNChargeExchangesTot );
		strcat( buf, buf2 );
		
		MPI_File_write( fh, buf, strlen(buf), MPI_CHAR, &status );
	
	}
	
	MPI_File_close( &fh );
	
}

DSMC::DSMC ( Species& neutral, Domain& domaind, std::string gasT ): neutrals{ neutral }, domain{ domaind }, gas{ gasT }, Fn{ neutral.getspwt() }
{
	//checks if chosen gas type is available and assigns atom mass and cross section parameters accordingly
	if ( gasT == "krypton" )
	{
		M = constants::m_kr;
		Tref = dsmcNeutralNeutralVHSScatter::Kr::T_ref;
		dref = dsmcNeutralNeutralVHSScatter::Kr::d_ref;
		omega = dsmcNeutralNeutralVHSScatter::Kr::omega;		
	}
	else if ( gasT == "neon" )
	{
		M = constants::m_ne;
		Tref = dsmcNeutralNeutralVHSScatter::Ne::T_ref;
		dref = dsmcNeutralNeutralVHSScatter::Ne::d_ref;
		omega = dsmcNeutralNeutralVHSScatter::Ne::omega;
	}
	else if ( gasT == "xenon" )
	{
		M = constants::m_xe;
		Tref = dsmcNeutralNeutralVHSScatter::Xe::T_ref;
		dref = dsmcNeutralNeutralVHSScatter::Xe::d_ref;
		omega = dsmcNeutralNeutralVHSScatter::Xe::omega;
	}
	else
	{
		std::cerr << "Invalid input for gas type! Please use neon, krypton or xenon." << std::endl;
	}
	
	//compute reduced mass
	mr = M / 2;
}

void DSMC::apply ( double dt )
{
	//seed prng for starting values
	std::mt19937 mersenne{ std::random_device{}() };
	
	//set up uniform random number distribution
	std::uniform_real_distribution<double> rN2( 0.0, 1.0 );
	
	//compute number of cells for given process
	int nr{ domain.getnr() };
	int nth{ domain.getnth() };
	int nZmin{ domain.getnZmin() };
	int nZmax{ domain.getnZmax() };
	int numCells{};
	if ( domain.getrank() == domain.getsize() - 1 )
	{
		numCells = ( nr - 1 )*nth*( nZmax - nZmin - 1 );
	}
	else
	{
		numCells = ( nr - 1 )*nth*( nZmax - nZmin );
	}
	
	//set up vector of particle adresses for every cell associated to given process
	std::vector<std::vector<Particle*>> partCell( numCells );
	
	//loop through all particles in neutrals species particle vector
	for ( auto& part : neutrals.part() )
	{
		//evaluate cell index of particle position and add particle adress to cell particle vector
		int cInd{ domain.posToCell( part.pos ) };
		
		/*
		if ( cInd > numCells - 1 )
		{
			std::cerr << "Error: Cell index of particle out of range!" << std::endl;
			std::cerr << "maximum index: " << numCells - 1 << " particle cell index: " << cInd << std::endl;
			std::cerr << "process range " << domain.getzMin() << " - " << domain.getzMax() << std::endl;
			std::cerr << " particle position: " << part.pos.r << ", " << part.pos.a << ", " << part.pos.z << std::endl;
			std::cerr << "grid number r: " << static_cast<int>( part.pos.r / domain.getdr() ) << std::endl;
			std::cerr << "grid number theta: " << static_cast<int>( part.pos.a / domain.getdth() ) << std::endl;
			std::cerr << "grid number z: " << static_cast<int>( part.pos.z / domain.getdz() ) << std::endl;
		}
		*/
		
		partCell.at(cInd).emplace_back( &part );
	}
	
	//loop through all cells
	for ( int c{ 0 }; c < numCells; ++c )
	{
		//initialize a reference to the particle adress vector of the cell
		std::vector<Particle*>& partAdress{ partCell.at(c) };
		int numPart{ static_cast<int>( partAdress.size() ) };
		double V{ domain.cellVolume( c ) };
		
		//check if there is more than one particle in cell
		if ( numPart > 1 )
		{
			//set up integer uniform random number distribution
			std::uniform_int_distribution<int> rN1( 0, numPart - 1 );
			
			//compute number of particle pairs to check for collision
			int numPairs{ static_cast<int>( 0.5*numPart*numPart*Fn*sigmaCrMax*dt/V + 0.5 ) };
			
			//loop through particle pairs
			for ( int p{ 0 }; p < numPairs; ++p )
			{
				//find random number for first particle
				int p1{ rN1( mersenne ) };
				int p2{};
				
				//find random number for second particle (must be unequal to number of first particle)
				do
				{
					p2 = rN1( mersenne );
				}
				while ( p2 == p1 );
				
				//initialize reference to particles
				Particle& part1{ *partAdress.at(p1) };
				Particle& part2{ *partAdress.at(p2) };
				
				//compute relative velocity of particles
				VectorCyl<double> vel1{ part1.vel.real() };
				VectorCyl<double> vel2{ part2.vel.real() };
				vel1.angleShift( -part1.vel.a );
				vel2.angleShift( -part2.vel.a );
				VectorCyl<double> relVel{ vel1 - vel2 };
				relVel.norma();
				double cr{ relVel.norm };
				
				//evaluate cross section
				double sigma{ evalSigma( cr ) };
				
				//update maximum cr times sigma
				if ( sigma*cr > sigmaCrMax ) { sigmaCrMax = sigma*cr; }
				
				//probability of collision
				double P{ sigma*cr / sigmaCrMax };
				
				//check if collision occurs
				if ( P > rN2( mersenne ) )
				{
					collide( part1, part2, vel1, vel2, cr, mersenne );
				}
				
			}
		}
	}
	
}

double DSMC::evalSigma( double cr )
{
	using std::pow;
	using std::tgamma;
	using constants::pi;
	using constants::k_B;
	
	//return cross section for variable hard sphere model (VHS) according to Bird
	return pi*dref*dref*pow( ( 2*k_B*Tref )/( mr*cr*cr ), omega - 0.5 )/tgamma( 2.5 - omega );
}

void DSMC::collide( Particle& p1, Particle& p2, VectorCyl<double>& vel1, VectorCyl<double> vel2, double cr, std::mt19937& mersenne )
{
	using std::sqrt;
	using std::cos;
	using std::sin;
	using constants::pi;
	
	//set up uniform random number distribution
	std::uniform_real_distribution<double> rN( 0.0, 1.0 );
	
	VectorCyl<double> cm{ vel1 + vel2 / 2 };				//for two particles of different mass cm = vel1*m1 + vel2*m2 / ( m1 + m2 )
	
	//compute angles for new random relative velocity
	double cosTh{ 2*rN(mersenne) - 1 };
	double sinTh{ sqrt( 1 - cosTh*cosTh ) };
	double phi{ 2*pi*rN(mersenne) };
	
	//compute new relative velocity
	VectorCyl<double> relVelNew{ cosTh*cr, sinTh*cos( phi )*cr, sinTh*sin( phi )*cr, 0 };
	//for 2D-simulation set theta-component of velocity to zero
	if ( domain.getnth() == 1 ) { relVelNew.theta = 0; }
	
	//compute new particle velocities
	VectorCyl<double> vel1New{ cm + relVelNew / 2 };
	VectorCyl<double> vel2New{ cm - relVelNew / 2 };
	vel1New.angleShift( p1.pos.a );
	vel2New.angleShift( p2.pos.a );
	
	//assign new velocities to particles
	p1.vel = vel1New.complexConversion();
	p2.vel = vel2New.complexConversion();
	
}

void DSMC::writeLogfile( int t )
{
	//write logfile for collisions of each process
	//set up mpi-file
	std::string filepath{ "logfiles/DSMCLogfilet" + std::to_string( t ) + ".dat" };
	MPI_Status status;
	MPI_File fh;
	MPI_File_open( MPI_COMM_WORLD, filepath.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );
	
	//set up output buffers
	char buf[71];
	char buf1[44];

	//write number of collision events to buffer
	snprintf( buf, 25, "Collisions process %i \n", domain.getrank() );
	snprintf( buf1, 43, "Elastic Scatter Collisions: %f \n", nColl );
	strcat( buf, buf1 );
	
	//write output buffer of every process to file
	MPI_File_write_ordered( fh, buf, strlen(buf), MPI_CHAR, &status );
	
	MPI_File_close( &fh );
	
	//write logfile for total number of collisions
	double nCollTot{ nColl };

	//for multiple processes evaluate total number of events for all collisions types
	if ( domain.getsize() != 1 )
	{
		MPI_Reduce( &nColl, &nCollTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	}

	//write total number of collisions to file
	if ( domain.getrank() == 0 )
	{
		std::ofstream lFile{ "logfiles/DSMCLogfileTotalt" + std::to_string( t ) + ".dat" };
		if ( !lFile )
		{
			std::cerr << "Logfile dsmc could not be opened!" << std::endl;
		}
		
		lFile << "Total Neutral-Neutral-Collisions\n";
		lFile << "Elastic Scatter Collisions: " << nCollTot << '\n';
		
		lFile.close();
	}
}

void ChemistryIonize::apply ( double dt ) {

	//seed prng for starting values
	std::mt19937 mersenne{ static_cast<std::mt19937::result_type>(std::time(nullptr)) };

	std::uniform_real_distribution<double> ranNum( 0.0, 1.0 );
	int nr{ domain.getnr() };
	int nth{ domain.getnth() };
	int nZmin{ domain.getnZmin() };
	int nZmax{ domain.getnZmax() };
	double dr{ domain.getdr() };
	double dth{ domain.getdth() };
	double dz{ domain.getdz() };
	
	for ( int i{ 0 }; i < nr - 1; ++i )
	{
		for ( int j{ 0 }; j < nth; ++j )
		{
			for ( int k{ nZmin }; k < nZmax - 1; ++k )
			{
				//volume of cell
				double dV{ constants::pi*((i+1)*(i+1)-i*i)*dr*dr*dz/nth };
				//neutral gas density
				double nN{ 1e14 };
				//electron density at the center of the cell
				double ne{ -domain.rhoField().gather( { (i+0.5)*dr, 0, (k+0.5)*dz, (j+0.5)*dth }, dr, dth, dz, nth )/constants::q_e };
				
				double dni = rate*nN*ne*dt;
				
				//number of macroparticles
				int numIons{ static_cast<int>( dni*dV/ions.getspwt() + ranNum(mersenne) ) };
				//add ions to species vector
				for ( int p{ 0 }; p < numIons; ++p )
				{
					VectorCyl<double> pos{ (i+ranNum(mersenne))*dr, 0, (k+ranNum(mersenne))*dz, (j+ranNum(mersenne))*dth };
					VectorCyl<double> vel{ 0, 0, 0, pos.a };
					ions.addParticleToAll( pos, vel );
				}
			}
		}
	}
	
	
}
 

