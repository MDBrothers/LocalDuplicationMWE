#include "kernels.hpp"
#include <cmath>


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


void computeInternalForceLinearElasticSimplifiedNew
(
		const double* __restrict__ xOverlap,
		const double* __restrict__ yOverlap,
		double* __restrict__ fInternalOverlap,
		double* __restrict__ fReactionsOverlap,
		const int* __restrict__ localIndexList,
		const int* __restrict__ neighborhoodLengths,
		const int numOwnedPoints
)
{

	/*
	 * Compute processor local contribution to internal force
	 */
	const double *xOwned, *yOwned;
	double *fOwned;

	const int *neighPtr = localIndexList;
	const int *neighLengths = neighborhoodLengths;
	double X_dx, X_dy, X_dz, zeta;
	double Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz;

	for(int p=0;p<numOwnedPoints;++ p, ++ neighLengths){
		// The first 'neighbor' is always the pivot
		xOwned = &xOverlap[*neighPtr * 3];	
		yOwned = &yOverlap[*neighPtr * 3];
		fOwned = &fInternalOverlap[*neighPtr * 3];

		int numNeigh = *neighLengths; neighPtr++;
		const double *X = xOwned;
		const double *Y = yOwned;

		for(int n=0;n<numNeigh;n++,neighPtr++){
			int localId = *neighPtr;
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

void computeInternalForceLinearElasticSimplifiedOld
(
		const double* __restrict__ xOverlap,
		const double* __restrict__ yOverlap,
		double* __restrict__ fInternalOverlap,
		double* __restrict__ fReactionsOverlap,
		const int* __restrict__ localIndexList,
		const int* __restrict__ neighborhoodLengths,
		const int numOwnedPoints
)
{

	/*
	 * Compute processor local contribution to internal force
	 */
	const double *xOwned, *yOwned;
	double *fOwned;

	const int *neighPtr = localIndexList;
	const int *neighLengths = neighborhoodLengths;
	double X_dx, X_dy, X_dz, zeta;
	double Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz;

	for(int p=0;p<numOwnedPoints;++ p, ++ neighLengths){
		// The first 'neighbor' is always the pivot
		xOwned = &xOverlap[*neighPtr * 3];	
		yOwned = &yOverlap[*neighPtr * 3];
		fOwned = &fInternalOverlap[*neighPtr * 3];

		int numNeigh = *neighLengths; neighPtr++;
		const double *X = xOwned;
		const double *Y = yOwned;

		for(int n=0;n<numNeigh;n++,neighPtr++){
			int localId = *neighPtr;
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
			fInternalOverlap[3*localId+0] -= fx;
			fInternalOverlap[3*localId+1] -= fy;
			fInternalOverlap[3*localId+2] -= fz;
		}

	}
}


