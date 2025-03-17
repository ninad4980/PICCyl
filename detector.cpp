//==========================================================================================
//
// Copyright (C) 2021 Dr. Ninad Joshi <email@ninadjoshi.de>  
// All Rights Reserved
//
//==========================================================================================

#include "detector.h"
#include "domain.h"
#include "species.h"
#include "vectorcyl.h"
#include <mpi.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>

void Detector::saveParticles()
{
	domain.timer().startTimer( "particle save" );
	
	//initialize variable to count number of detected particles
	int numP{ 0 };
	
	//check if there are particles of given species in detector
	if ( species.getnDetector() > 0 )
	{
		//iterate through all particles of given species
		for ( Particle& p : species.part() )
		{
			//check if particle p has id corresponding to detector grid point 
			if ( p.id == 3 )
			{
				//save particles p in particles vector of detector
				Particle part{ p };
				particles.emplace_back( part );
				//set particle p.id to zero and unset p.alive in species particles vector
				p.id = 0;
				p.alive = false;
				++numP;
			}
		}
	}
	
	//initialize variable for total number of detected particles
	int numPAll{ 0 };
	
	//evaluate total number of detected particles
	MPI_Allreduce( &numP, &numPAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	
	//update maximum of detected particles
	if ( numPAll > maxPartFlow ) { maxPartFlow = numPAll; }
	
	//add number of total detected particles to numParticlesPerTime vector
	numParticlesPerTime.emplace_back( numPAll );
	
	domain.timer().stopTimer( "particle save" );
}

bool Detector::constParticleFlow( double relDiffMax, int avgNum )
{
	using std::abs;
	
	//evaluate iteration step
	int itStep{ static_cast<int>( numParticlesPerTime.size() ) };
	
	//check if iteration step is larger than 2 times avgNum
	if ( itStep > 3*avgNum )
	{		
		//initialize variables for avarage detected particles
		double avgPartFlowLast{ 0 };
		double avgPartFlowPrevious{ 0 };
		double avgPartFlowPrevious2{ 0 };
		
		//evaluate avarage detected particles in last avgNum iteration steps and the avgNum steps previous to that
		for ( int st{ itStep - avgNum }; st < itStep; ++st )
		{
			avgPartFlowLast += static_cast<double>( numParticlesPerTime.at( st ) )/ avgNum;
			avgPartFlowPrevious += static_cast<double>( numParticlesPerTime.at( st - avgNum ) )/ avgNum;
			avgPartFlowPrevious2 += static_cast<double>( numParticlesPerTime.at( st - 2*avgNum ) )/ avgNum;
		}
		
		//std::cout << "avgPartFlowLast: " << avgPartFlowLast << '\n';
		//std::cout << "avgPartFlowPrevious: " << avgPartFlowPrevious << '\n';
		
		
		//check if avarage detected particle numbers are larger than zero
		if ( ( avgPartFlowLast > 0 ) && ( avgPartFlowPrevious > 0 ) )
		{
			//compute relative difference between avarage detected particle numbers
			double relDiff{ abs( avgPartFlowLast - avgPartFlowPrevious ) / avgPartFlowLast };
			double relDiff2{ abs( avgPartFlowLast - avgPartFlowPrevious2 ) / avgPartFlowLast };
			
			//std::cout << "relDiff: " << relDiff << '\n';
			
			//return true if relative difference is smaller than specified value
			if ( ( relDiff < relDiffMax ) && ( relDiff2 < relDiffMax ) ) { return true; }
		}
	}
	
	//return false in all other cases
	return false;
}

bool Detector::maxParticleFlowReached( int numSteps )
{
	int timeSteps{ static_cast<int>( numParticlesPerTime.size() ) };
	
	if ( ( timeSteps > numSteps ) && ( maxPartFlow > 0 ) )
	{
		for ( int it{ timeSteps - numSteps }; it < timeSteps; ++it )
		{
			if ( numParticlesPerTime.at( it ) > 0.99*maxPartFlow ) { return false; }
		}
		
		return true;
	}
	
	return false;
}

void Detector::calcEmittance()
{
	using std::atan2;
	using std::sqrt;
	using std::fabs;
	
	long unsigned int Np{ particles.size() };
	double rAvg{ 0 };
	double r2Avg{ 0 };
	double rpAvg{ 0 };
	double rp2Avg{ 0 };
	double rrpSum{ 0 };
	double rSum{ 0 };
	double rpSum{ 0 };
	
	if ( Np > 0 )
	{
		for ( auto& p : particles )
		{
			rAvg += p.pos.r / Np;
			r2Avg += p.pos.r*p.pos.r / Np;
			rpAvg += atan2( p.vel.r.real(), p.vel.z.real() ) / Np;
			rp2Avg += atan2( p.vel.r.real(), p.vel.z.real() )*atan2( p.vel.r.real(), p.vel.z.real() ) / Np;
			rrpSum += p.pos.r*atan2( p.vel.r.real(), p.vel.z.real() );
			rSum += p.pos.r;
			rpSum += atan2( p.vel.r.real(), p.vel.z.real() );
		}
		
		double r2{ r2Avg - rAvg*rAvg };
		double rp2{ rp2Avg - rpAvg*rpAvg };
		
		double rskew{ ( rrpSum - rSum*rpSum / Np ) / Np };
		
		species.setEmittance( sqrt( fabs( r2*rp2 - rskew*rskew ) ) );		
	}
}

void Detector::calcyEmittance()
{
	using std::atan2;
	using std::sqrt;
	using std::fabs;
	
	//get number of particles stored in detector
	long unsigned int Np{ particles.size() };
	double yAvg{ 0 };
	double y2Avg{ 0 };
	double ypAvg{ 0 };
	double yp2Avg{ 0 };
	double yypSum{ 0 };
	double ySum{ 0 };
	double ypSum{ 0 };
	
	//check if particles number in detector is greater zero
	if ( Np > 0 )
	{
		//loop through all particles stored in detector
		for ( auto& p : particles )
		{
			//compute Cartesian components of position and velocity vectors
			VectorCyl<double> posCar{ p.pos };
			posCar.angleShift( -posCar.a );
			VectorCyl<double> velCar{ p.vel.real() };
			velCar.angleShift( -velCar.a );
			
			//compute average y position
			yAvg += posCar.theta / Np;
			//compute average of the square of y position
			y2Avg += posCar.theta*posCar.theta / Np;
			//compute average of the angle between y- and z-component of velocity
			ypAvg += atan2( velCar.theta, velCar.z ) / Np;
			//compute average of the square of the angle between y- and z-component of velocity
			yp2Avg += atan2( velCar.theta, velCar.z )*atan2( velCar.theta, velCar.z ) / Np;
			//compute sum of product between y position and angle between y- and z-componet of velocity
			yypSum += posCar.theta*atan2( velCar.theta, velCar.z );
			//compute sum of y position
			ySum += posCar.theta;
			//compute sum of angle between y- and z-component of velocity
			ypSum += atan2( velCar.theta, velCar.z );
		}
		
		//compute variance of y position
		double y2{ y2Avg - yAvg*yAvg };
		//compute variance of angle between y- and z-component of velocity
		double yp2{ yp2Avg - ypAvg*ypAvg };
		
		//compute angle-position correlation of particles
		double yskew{ ( yypSum - ySum*ypSum / Np ) / Np };
		
		//compute and set rms-y-emittance
		species.setyEmittance( sqrt( fabs( y2*yp2 - yskew*yskew ) ) );		
	}
}

void Detector::writeParticleFlowToFile( int stepNum, double dt )
{
	//check if process rank is zero
	if ( domain.getrank() == 0 )
	{
		//compute number of iteration steps to avarage over
		int avgNum{ static_cast<int>( numParticlesPerTime.size() ) / stepNum };
		
		//open file stream fout for particle flow in dependence of time
		std::ofstream fout{ "results/particleFlow.dat" };
		if ( !fout )
		{
			std::cerr << "particleFlow.dat could not ne opened!" << std::endl;
		}
		
		//iterate over number of time steps to write down in file
		for ( int it{ 0 }; it < stepNum; ++ it )
		{
			//initialize variable for avarage number of particles
			double partAvg{ 0 };
			
			//iterate over iteration steps to avarage over
			for ( int itS{ it*avgNum }; itS < ( it + 1 )*avgNum; ++itS )
			{
				//compute avarage number of particles
				partAvg += static_cast<double>( numParticlesPerTime.at( itS ) )/ avgNum;
			}
			
			//write down avarage number of detected particles per time to file
			fout << ( it + 1 )*avgNum*dt << '\t' << partAvg << '\n';
		}
		
		fout.close();
		
	}
}

void Detector::writePhasespaceToFile( int t )
{
	using std::atan2;
	
	//get number of particles in detector
	int partNum{ static_cast<int>( particles.size() ) };
	int maxNum{};
	
	//evaluate maximum number of particles in detector across all processes
	MPI_Allreduce( &partNum, &maxNum, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	
	//set up mpi-file
	std::string filepath{ "results/phasespaceDetectort" + std::to_string( t ) + ".dat" };
	MPI_Status status;
	MPI_File fh;
	MPI_File_open( MPI_COMM_WORLD, filepath.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );
	
	//set up buffer for particles phasespace
	char buf[ (maxNum+1)*74 ];
	
	//write number of process to buffer
	snprintf( buf, 47, "phasespace process %i \n", domain.getrank() );
	//loop through all particles in detector particles vector
	for ( Particle& p : particles )
	{
		//write particle position and velocity to buffer
		char buf1[74];
		snprintf( buf1, 73, "%.4e %.4e %.4e %.4e %.4e %.4e \n", p.pos.r, p.pos.a, p.pos.z, p.vel.r.real(), p.vel.theta.real(), p.vel.z.real() );
		strcat( buf, buf1 );
	}
		
	//write buffer to file in order of process
	MPI_File_write_ordered( fh, buf, strlen(buf), MPI_CHAR, &status );
	
	MPI_File_close( &fh );
}

void Detector::writePositionsRTheta( int t )
{
	using std::atan2;
	
	//get number of particles in detector
	int partNum{ static_cast<int>( particles.size() ) };
	int maxNum{};
	
	//evaluate maximum number of particles in detector across all processes
	MPI_Allreduce( &partNum, &maxNum, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	
	//set up mpi-file
	std::string filepath{ "results/particlePosRTheta" + std::to_string( t ) + ".dat" };
	MPI_Status status;
	MPI_File fh;
	MPI_File_open( MPI_COMM_WORLD, filepath.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );
	
	//set up buffer for particles phasespace
	char buf[ (maxNum+1)*27 ];
	
	//loop through all particles in detector particles vector
	for ( Particle& p : particles )
	{
		//write particle position and velocity to buffer
		char buf1[27];
		snprintf( buf1, 26, "%.4e %.4e\n", p.pos.r, p.pos.a );
		strcat( buf, buf1 );
	}
		
	//write buffer to file in order of process
	MPI_File_write_ordered( fh, buf, strlen(buf), MPI_CHAR, &status );
	
	MPI_File_close( &fh );
}
