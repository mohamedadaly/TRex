/*
-----------------------------------------------------------------------
Copyright: 2010-2015, iMinds-Vision Lab, University of Antwerp
           2014-2015, CWI, Amsterdam

Contact: astra@uantwerpen.be
Website: http://sf.net/projects/astra-toolbox

This file is part of the ASTRA Toolbox.


The ASTRA Toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The ASTRA Toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the ASTRA Toolbox. If not, see <http://www.gnu.org/licenses/>.

-----------------------------------------------------------------------
$Id$
*/

#ifndef _INC_ASTRA_TREXALGORITHMTERMS
#define _INC_ASTRA_TREXALGORITHMTERMS

#include "Globals.h"
#include "Config.h"

#include "AstraObjectManager.h"
#include "Float32ProjectionData2D.h"
#include "Float32VolumeData2D.h"
#include "Float32VolumeData3DMemory.h"
#include "Logging.h"
#include "ReconstructionAlgorithm2D.h"


namespace astra {

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
//							Priors and Data terms for TRex
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
	
// ----------------------------------------------------------------------------
// Priors
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Base class for priors

class _AstraExport CTRexPrior {
public:
	float32 m_fSigma;
	CTRexPrior() : m_fSigma(1.f) {}	
	CTRexPrior(float32 _fSigma) : m_fSigma(_fSigma) {}	
	virtual ~CTRexPrior() {}

	// Multiply the volume by the matrix K
	virtual void K(const CFloat32VolumeData2D* pX,  
		CFloat32VolumeData3DMemory* pKx) = 0;

	// Multiply the volume by the transpose of K
	virtual void Kt(const CFloat32VolumeData3DMemory* pKx, 
		CFloat32VolumeData2D* pX) = 0;

	// Apply the proximal operator on the input u (of size Kx).
	virtual void prox(const CFloat32VolumeData3DMemory* pU, float32 rho,
		CFloat32VolumeData3DMemory* pV) = 0;

	// Return the squared norm of K ||K||^2_2
	virtual float32 norm() = 0;

	// Return the depth of the object returned by applying K
	virtual int depth() = 0;
};

// ----------------------------------------------------------------------------
// Anisotropic TV
class _AstraExport CTRexPriorATV : public CTRexPrior {
public:
	CTRexPriorATV() : CTRexPrior() {}
	CTRexPriorATV(float _fSigma) : CTRexPrior(_fSigma) {}

	// Multiply the volume by the matrix K * sigma
	virtual void K(const CFloat32VolumeData2D* pX, 
		CFloat32VolumeData3DMemory* pKx);

	// Multiply the volume by the transpose of K * sigma
	virtual void Kt(const CFloat32VolumeData3DMemory* pKx, 
		CFloat32VolumeData2D* pX);

	// Apply the proximal operator on the input u (of size Kx).
	virtual void prox(const CFloat32VolumeData3DMemory* pU, float32 rho,
		CFloat32VolumeData3DMemory* pV);

	// Return the squared norm of K ||K||^2_2
	virtual float32 norm();

	virtual int depth() { return 2; }
};

// ----------------------------------------------------------------------------
// Isotropic TV
class _AstraExport CTRexPriorITV : public CTRexPriorATV {
public:
	CTRexPriorITV() : CTRexPriorATV() {}
	CTRexPriorITV(float _fSigma) : CTRexPriorATV(_fSigma) {}

	// Apply the proximal operator on the input u (of size Kx).
	virtual void prox(const CFloat32VolumeData3DMemory* pU, float32 rho,
		CFloat32VolumeData3DMemory* pV);

};

// ----------------------------------------------------------------------------
// Sum of Absolute Differences
class _AstraExport CTRexPriorSAD : public CTRexPriorATV {
public:
	CTRexPriorSAD() : CTRexPriorATV() {}
	CTRexPriorSAD(float _fSigma) : CTRexPriorATV(_fSigma) {}

	// Multiply the volume by the matrix K * sigma
	virtual void K(const CFloat32VolumeData2D* pX, 
		CFloat32VolumeData3DMemory* pKx);

	// Multiply the volume by the transpose of K * sigma
	virtual void Kt(const CFloat32VolumeData3DMemory* pKx, 
		CFloat32VolumeData2D* pX);

	// Return the squared norm of K ||K||^2_2
	virtual float32 norm();

	virtual int depth() { return 8; }

};
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Data Terms
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Base class for data terms
class _AstraExport CTRexData {
public:
	CReconstructionAlgorithm2D* m_pAlg;
	Config* m_pCfg;
	// ProxInput for the prox operator
	CFloat32VolumeData2D* m_pProxInput;

	CTRexData();
	virtual ~CTRexData();

	// Initialize the algorithm
	virtual bool init(const Config& cfg, CVolumeGeometry2D* geom,
					   float32 fLambda);
	// Runs the algorithm for iter iterations
	virtual void prox(int32 iter);
protected:
	// Parses the input config and sets it up to be initialized
	virtual bool parse(const Config& cfg, CVolumeGeometry2D* geom,
					   float32 fLambda) = 0;
};

// Least Squares data term
class _AstraExport CTRexDataLS : public CTRexData {
protected:
	virtual bool parse(const Config& cfg, CVolumeGeometry2D* geom,
					   float32 fLambda);
};

// Weighted least squares
class _AstraExport CTRexDataWLS : public CTRexDataLS {
public:
	CTRexDataWLS() 
	{
		m_pW = NULL;
		m_iWlsRoot = 1;
		ASTRA_INFO("WlsRoot const = %d", m_iWlsRoot);
	}
protected:
	// Weights for WLS that multiply the sinogram and the projection matrix.
	// They are processed here during init to have the correct form e.g. sqrt
	// if WlsRoot = 1.
	CFloat32ProjectionData2D* m_pW;

	// WLS nth root to apply before feeding to the prox operator for WLS data 
	// term
	int m_iWlsRoot;

	virtual bool parse(const Config& cfg, CVolumeGeometry2D* geom,
					   float32 fLambda);
};

} // namespace

#endif