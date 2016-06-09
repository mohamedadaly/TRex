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

#include "astra/TRexAlgorithmTerms.h"

#include "astra/ReconstructionAlgorithm2D.h"
#include "astra/AstraObjectFactory.h"

using namespace std;

namespace astra {

// ---------------------------------------------------------------------------
// Priors
// ---------------------------------------------------------------------------

void CTRexPriorATV::K(const CFloat32VolumeData2D* pX, 
					  CFloat32VolumeData3DMemory* pKx) 
{
	// Forward discrete difference in 2D
	//
	// assert the size of Kx (width, height, 2)
	ASTRA_ASSERT(pKx->getHeight() == pX->getHeight() &&
		pKx->getWidth() == pX->getWidth() &&
		pKx->getDepth() == this->depth());

	int rows = pX->getHeight();
	int cols = pX->getWidth();

	// Init
	float32* pDataOut = pKx->getData();
	const float32* pDataIn = pX->getDataConst();

	// Compute horizontal difference in first slice
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols - 1; ++j) {
			// The data is written in row-major order
			int iindx = i * cols + j;			
			// first slice
			int oindx = iindx;
			pDataOut[oindx] = pDataIn[iindx + 1] - pDataIn[iindx];
		}
	}

	// Compute vertical difference in second slice
	for (int i = 0; i < rows - 1; ++i) {
		for (int j = 0; j < cols; ++j) {
			int iindx = i * cols + j;		
			// second slice
			int oindx = iindx + rows * cols; 
			pDataOut[oindx] = pDataIn[iindx + cols] - pDataIn[iindx];
		}
	}

	//pKx->printInfo("Kx before");

	// Multiply by sigma
	if (m_fSigma != 1.f) {
		*pKx *= m_fSigma;
	}

	//ASTRA_INFO("Sigma: %f", sigma);
	//pX->printInfo("X");
	//pKx->printInfo("Kx after");
}

void CTRexPriorATV::Kt(const CFloat32VolumeData3DMemory* pKx, 
					   CFloat32VolumeData2D* pX)
{
	// Negative divergence in 2D (transpose of forward difference)
	//
	ASTRA_ASSERT(pKx);
	ASTRA_ASSERT(pX);
	ASTRA_ASSERT(pKx->getHeight() == pX->getHeight() &&
		pKx->getWidth() == pX->getWidth() &&
		pKx->getDepth() == this->depth());


	// Init
	pX->setData(0.f);
	float32* pDataOut = pX->getData();
	const float32* pDataIn = pKx->getDataConst();
	
	int rows = pX->getHeight();
	int cols = pX->getWidth();

	// First slice: horizontal difference. -ve in place and +ve to the left
	//
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			int oindx = i * cols + j;
			// first slice
			int iindx = oindx;
			//ASTRA_ASSERT(indx < pX->getSize());
			//ASTRA_ASSERT(indx-1 >= 0 && indx-1 < pSlice->getSize());
			if (j > 0)
				// add the value on the left and subtract in place
				pDataOut[oindx] += pDataIn[iindx - 1] - pDataIn[iindx];
			else
				pDataOut[oindx] -= pDataIn[iindx];
		}
	}

	// Second slice: vertical difference. -ve in place and +ve to the top
	//
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			int oindx = i * cols + j;
			// second slice
			int iindx = oindx + cols*rows;
			//ASTRA_ASSERT(indx < pX->getSize());
			//ASTRA_ASSERT(indx-cols >= 0 && indx-cols < pSlice->getSize());
			if (i > 0)
				// add the value on the top and subtract in place
				pDataOut[oindx] += pDataIn[iindx - cols] - pDataIn[iindx];
			else
				pDataOut[oindx] -= pDataIn[iindx];
		}
	}
	//pX->printInfo("X before");

	// Multiply by sigma
	if (m_fSigma != 1.f) {
		*pX *= m_fSigma;
	}

	//pKx->printInfo("Kx");
	//pX->printInfo("X after");
}

void CTRexPriorATV::prox(const CFloat32VolumeData3DMemory* pU, float32 rho,
						 CFloat32VolumeData3DMemory* pV)
{
	ASTRA_ASSERT(pU->getDepth() == pV->getDepth() &&
		pU->getSize() == pV->getSize());

	const float32* pDataIn = pU->getDataConst();
	float32* pDataOut = pV->getData();
	// loop on values
	int sz = pU->getSize();
	for (int i = 0; i < sz; ++i) {
		pDataOut[i] = max(0.f, pDataIn[i] - rho) - 
			max(0.f, -pDataIn[i] - rho);
	}

	//pU->printInfo("U");
	//pV->printInfo("V");
}

float32 CTRexPriorATV::norm()
{
	return 8.f * m_fSigma * m_fSigma;
}

void CTRexPriorITV::prox(const CFloat32VolumeData3DMemory* pU, float32 rho,
						 CFloat32VolumeData3DMemory* pV)
{
	//ASTRA_INFO("Depth: %d & %d", pU->getSize(), pV->getSize());
	ASTRA_ASSERT(pU->getDepth() == pV->getDepth() &&
		pU->getDepth() == this->depth() &&
		pU->getSize() == pV->getSize());

	int rows = pU->getHeight();
	int cols = pU->getWidth();
	int depth = pU->getDepth();

	const float32* pDataIn = pU->getDataConst();
	float32* pDataOut = pV->getData();

	// Compute norm across third dimension
	CVolumeGeometry2D geo(cols, rows);
	CFloat32VolumeData2D* pNorm = new CFloat32VolumeData2D(&geo);
	float32* pNormData = pNorm->getData();

	int offset = rows * cols;
	#pragma omp for
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			int indx = i * cols + j;
			float32 val = pDataIn[indx] * pDataIn[indx] +
				pDataIn[indx + offset] * pDataIn[indx + offset];
			pNormData[indx] = sqrtf(val);
		}
	}

	// Compute output: V = U - U*rho / max(rho, norm)
	for (int k = 0; k < this->depth(); ++k) {
		offset = k * rows * cols;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				int indx2 = i * cols + j;
				int indx3 = indx2 + offset;
				pDataOut[indx3] = pDataIn[indx3] -
					(rho * pDataIn[indx3] / max(rho, pNormData[indx2]));
			}
		}
	}

	// clear
	ASTRA_DELETE(pNorm);

	//pU->printInfo("U");
	//pV->printInfo("V");
}

void CTRexPriorSAD::K(const CFloat32VolumeData2D* pX, 
					  CFloat32VolumeData3DMemory* pKx) 
{
	// Forward discrete difference in 2D in 8 directions
	//
	// assert the size of Kx (width, height, 8)
	ASTRA_ASSERT(pKx->getHeight() == pX->getHeight() &&
		pKx->getWidth() == pX->getWidth() &&
		pKx->getDepth() == this->depth());

	int r = pX->getHeight();
	int c = pX->getWidth();

	// Init
	pKx->setData(0.f);
	float32* pDataOut = pKx->getData();
	const float32* pDataIn = pX->getDataConst();

	// Differences in this order
	//                     6     7     8
	//                     5     x     1
	//                     4     3     2
	//                   1    2    3   4     5     6   7  8
	int32 iStart[] =	{0,	  0,   0,   0,   0,    1,  1, 1   };
	int32 iEnd[] =		{r,   r-1, r-1, r-1, r,    r,  r, r   };
	int32 jStart[] =	{0,   0,   0,   1,   1,    1,  0, 0   };
	int32 jEnd[] =		{c-1, c-1, c,   c,   c,    c,  c, c-1 };
	int32 off[] =       {1,   c+1, c, c-1,  -1, -c-1, -c, -c+1};

	// Compute horizontal difference in first slice
	for (int k = 0; k < pKx->getDepth(); ++k) {
		for (int i = iStart[k]; i < iEnd[k]; ++i) {
			for (int j = jStart[k]; j < jEnd[k]; ++j) {
				// The data is written in row-major order
				int iindx = i * c + j;			
				// k th slice
				int oindx = iindx + k * r*c;
				// Put the value in the correct slice
				pDataOut[oindx] = pDataIn[iindx + off[k]] - pDataIn[iindx];
			}
		}
	}

	//pKx->printInfo("Kx before");

	// Multiply by sigma
	if (m_fSigma != 1.f) {
		*pKx *= m_fSigma;
	}

	//ASTRA_INFO("Sigma: %f", sigma);
	//pX->printInfo("X");
	//pKx->printInfo("Kx after");
}

void CTRexPriorSAD::Kt(const CFloat32VolumeData3DMemory* pKx, 
					   CFloat32VolumeData2D* pX)
{
	// Negative divergence in 2D (transpose of forward difference)
	//
	ASTRA_ASSERT(pKx);
	ASTRA_ASSERT(pX);
	ASTRA_ASSERT(pKx->getHeight() == pX->getHeight() &&
		pKx->getWidth() == pX->getWidth() &&
		pKx->getDepth() == this->depth());


	// Init
	pX->setData(0.f);
	float32* pDataOut = pX->getData();
	const float32* pDataIn = pKx->getDataConst();
	
	int rows = pX->getHeight();
	int cols = pX->getWidth();

	// Compute the sum across the slices
	//
	for (int k = 0; k < pKx->getDepth(); ++k) {
		int offset = k * rows * cols;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				int oindx = i * cols + j;
				int iindx = oindx + offset;
				pDataOut[oindx] += pDataIn[iindx];
			}
		}
	}
	// Multiply by -2
	*pX *= -2.f;

	// Multiply by sigma
	if (m_fSigma != 1.f) {
		*pX *= m_fSigma;
	}

	//pKx->printInfo("Kx");
	//pX->printInfo("X after");
}

float32 CTRexPriorSAD::norm()
{
	return 24.f * m_fSigma * m_fSigma;
}


CTRexData::CTRexData() {
	m_pAlg = NULL;
	m_pCfg = NULL;
	m_pProxInput = NULL;
}

CTRexData::~CTRexData() 
{
	ASTRA_DELETE(m_pAlg);
	ASTRA_DELETE(m_pCfg);
	CData2DManager::getSingleton().remove(
		CData2DManager::getSingleton().getIndex(m_pProxInput));
	m_pProxInput = NULL;
}

bool CTRexData::init(const Config& cfg, CVolumeGeometry2D* geom,
					float32 fLambda) 
{
	// parse
	if (!parse(cfg, geom, fLambda))
		return false;
	// Init
	ASTRA_ASSERT(m_pAlg);
	ASTRA_ASSERT(m_pCfg);
	// Initialize the prox algorithm
	return m_pAlg->initialize(*m_pCfg);

}
void CTRexData::prox(int32 iter) 
{
	ASTRA_ASSERT(m_pAlg);
	// run for the required iterations
	m_pAlg->run(iter);
}

// ----------------------------------------------------------------------------
bool CTRexDataLS::parse(const Config& cfg, CVolumeGeometry2D* geom,
					   float32 fLambda)
{
	ASTRA_ASSERT(geom);

	// copy config
	m_pCfg = new Config(cfg);
	ASTRA_ASSERT(m_pCfg);

	// Get the data proximal operator algorithm
	string sProx = m_pCfg->self.getOption("DataProx","");
	m_pAlg = dynamic_cast<CReconstructionAlgorithm2D*>(
		CAlgorithmFactory::getSingleton().create(sProx));
	ASTRA_CONFIG_CHECK(m_pAlg, "CTRexAlgorithm", 
		"Error initializing: unknown DataProx.");

	// Turning it off for now.
	// TODO: check.
	m_pCfg->self.removeOption("IterationMetricsId");
	m_pCfg->self.removeOption("ComputeIterationMetrics");
	//if (m_bComputeIterationMetrics) {
	//	// Create a matrix for that
	//	m_pProxMetrics = new CFloat32VolumeData2D();
	//	// Store in manager and save Id to pass to prox solver
	//	id = CData2DManager::getSingleton().store(m_pProxMetrics);
	//	// Remove the old value and add a new one
	//	m_pData->m_pCfg->self.removeOption("IterationMetricsId");
	//	m_pData->m_pCfg->self.addOption("IterationMetricsId", 
	//		static_cast<float32>(id));
	//}

	// Add ProxInput to the prox config
	// Create object and add to manager
	m_pProxInput = new CFloat32VolumeData2D(geom);
	int id = CData2DManager::getSingleton().store(m_pProxInput);
	// Add the id to config
	m_pCfg->self.removeChildNode("ProxInputDataId");
	m_pCfg->self.addChildNode("ProxInputDataId", 
		static_cast<float32>(id));

	// Set Lambda
	m_pCfg->self.removeOption("Lambda");
	m_pCfg->self.addOption("Lambda", fLambda);

	return true;
}

bool CTRexDataWLS::parse(const Config& cfg, CVolumeGeometry2D* geom,
					   float32 fLambda)
{
	// Call the parent function to parse the common stuff.
	if (!CTRexDataLS::parse(cfg, geom, fLambda))
		return false;

	// WLS root
	m_iWlsRoot = static_cast<int>(
		m_pCfg->self.getOptionNumerical("WlsRoot", m_iWlsRoot));
	//CC.markOptionParsed("WlsRoot");
	ASTRA_CONFIG_CHECK(m_iWlsRoot > 0, "CTRexAlgorithm", 
		"Error initializing: WlsRoot <= 0");

	// WLS weights
	int32 id = static_cast<int>(
		m_pCfg->self.getOptionNumerical("WlsWeightDataId", -1));
	m_pW = dynamic_cast<CFloat32ProjectionData2D*>(
		CData2DManager::getSingleton().get(id));
	ASTRA_CONFIG_CHECK(m_pW, "CTRexAlgorithm", 
		"Error initializing: Invalid WlsWeightDataId");	
	//CC.markOptionParsed("WlsWeightDataId");
	//ASTRA_CONFIG_CHECK(m_pW, "CTRexAlgorithm", 
	//	"Error initializing: unknown DataProx.");

	// Put in the correct form by taking sqrt of W
	for (int i = m_iWlsRoot; i > 0; --i) {
		m_pW->sqrt();
	}

	return true;
}

} // namespace
