#include "kernels.hpp"
#include <cmath>

// The comparison would not be fair if the same kernel weren't used for both data layouts. So we write in a way that suits the old way more.
void computeBondSumAny(const int DIMENSION, const int NUM_NEIGHBORHOODS, const int* __restrict__ NEIGHBORHOOD_LENGTHS, const int* __restrict__ OWNED_INDICES, const int* __restrict__ NEIGHBOR_INDICES, const double* __restrict__ CURR_COORDS, double* __restrict__ bond_sum, double* __restrict__ neighbor_reactions){
	std::vector<double> components(DIMENSION);
	double *componentsPtr, *ownedSum, *ownedCurrCoords;
	double mag(0);
	int *nid(NEIGHBOR_INDICES), ownedID; 

	for(int neighborhood(0); neighborhood < NUM_NEIGHBORS; ++ neighborhood){
		ownedID = OWNED_INDICES[neighborhood];
		ownedDofID = ownedID*DIMENSION;
		ownedSum = &bond_sum[ownedDofID];
		ownedCurrCoords = &CURR_COORDS[ownedDofID];

		for(int neighbor(0); neighbor < NEIGHBORHOOD_LENGTHS[neighborhood]; ++ neighbor){
			
			const int NEIGHBOR_SCALED_INDEX(*(++ nid)*DIMENSION);
			mag = 0.0;
			componentsPtr = &components[0];

			// Compute current length components and magnitude
			for(int dof(0); dof < DIMENSION; ++dof, ++ componentsPtr){
				*componentsPtr = CURR_COORDS[NEIGHBOR_SCALED_INDEX + dof] - *(ownedCurrCoords + dof); 
				mag += *componentsPtr * (*componentsPtr);
			}
			mag = sqrt(mag);

			// Compute sum of bond magnitudes as well as 'reactions'
			componentsPtr = &components[0];
			for(int dof(0); dof < DIMENSION; ++dof, ++ componentsPtr){
				bond_sum[ownedDofID + dof] += mag/(*componentsPtr);
				neighbor_reactions[NEIGHBOR_SCALED_INDEX + dof] -= mag/(*componentsPtr);
			}
		}

		componentsPtr = &components[0];
		for(int dof(0); dof < DIMENSION; ++ dof){
			*componentsPtr = CURR_COORDS[ dof] - CURR_COORDS[+ dof];

		}	
	}	
};
//! \file elastic.cxx 

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include <cmath>

void computeInternalForceLinearElasticSimplified
(
		const double* xOverlap,
		const double* yOverlap,
		double* fInternalOverlap,
		double* fReactionsOverlap,
		const int*  localIndexList,
		const int*  neighborhoodPreLengths,
		const int numOwnedPoints
)
{

	/*
	 * Compute processor local contribution to internal force
	 */
	const double *xOwned;
	const double *yOwned;
	double *fOwned;

	const int *neighPtr = localIndexList;
	const int *neighLengths = neighborhoodPreLengths;
	double X_dx, X_dy, X_dz, zeta;
	double Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz;
	for(int p=0;p<numOwnedPoints;p++, xOwned = (, yOwned +=3, fOwned+=3,  ++ neighLengths){
		xOwned = &xOverlap[*neighLengths * 3];	
		yOwned = &yOverlap[*neighLengths * 3];


		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const double *Y = yOwned;
		alpha = 15.0*MU/(*m);
		double selfCellVolume = v[p];
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const double *YP = &yOverlap[3*localId];
			X_dx = XP[0]-X[0];
			X_dy = XP[1]-X[1];
			X_dz = XP[2]-X[2];
			zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
			Y_dx = YP[0]-Y[0];
			Y_dy = YP[1]-Y[1];
			Y_dz = YP[2]-Y[2];
			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
      e = dY - zeta;

			t = zeta * e;
			fx = t * Y_dx / dY;
			fy = t * Y_dy / dY;
			fz = t * Y_dz / dY;

			*(fOwned+0) += fx;
			*(fOwned+1) += fy;
			*(fOwned+2) += fz;
			fReactionsOverlap[3*localId+0] -= fx;
			fReactionsOverlap[3*localId+1] -= fy;
			fReactionsOverlap[3*localId+2] -= fz;

		}

	}
}


