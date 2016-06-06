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

#include "astra/TRexAlgorithm.h"

#include <boost/lexical_cast.hpp>

#include "astra/AstraObjectManager.h"
#include "astra/DataProjectorPolicies.h"
#include "astra/Logging.h"
#include "astra/Float32VolumeData3DMemory.h"


using namespace std;

namespace astra {

// ---------------------------------------------------------------------------
// Priors
// ---------------------------------------------------------------------------

void CTRexPriorATV::K(const CFloat32VolumeData2D* pX, 
					  CFloat32VolumeData3D* pKx, float32 sigma) 
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
	CFloat32VolumeData2D* pSlice;
	float32* pDataOut;
	const float32* pDataIn = pX->getDataConst();

	// Compute horizontal difference in first slice
	pSlice = pKx->fetchSliceZ(0);
	pSlice->setData(0.f);
	pDataOut = pSlice->getData();	
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols - 1; ++j) {
			// The data is written in row-major order
			int indx = i * cols + j;			
			pDataOut[indx] = pDataIn[indx + 1] - pDataIn[indx];
		}
	}
	pKx->returnSliceZ(0, pSlice);

	// Compute vertical difference in second slice
	pSlice = pKx->fetchSliceZ(1);
	pSlice->setData(0.f);
	pDataOut = pSlice->getData();	
	for (int i = 0; i < rows - 1; ++i) {
		for (int j = 0; j < cols; ++j) {
			int indx = i * cols + j;			
			pDataOut[indx] = pDataIn[indx + cols] - pDataIn[indx];
		}
	}
	pKx->returnSliceZ(1, pSlice);

	// Multiply by sigma
	if (sigma != 1.f) {
		*pKx *= sigma;
	}
}

void CTRexPriorATV::Kt(const CFloat32VolumeData3D* pKx, 
					   CFloat32VolumeData2D* pX, float32 sigma)
{
	// Negative divergence in 2D (transpose of forward difference)
	//
	ASTRA_ASSERT(pKx->getHeight() == pX->getHeight() &&
		pKx->getWidth() == pX->getWidth() &&
		pKx->getDepth() == this->depth());


	// Init
	pX->setData(0.f);
	float32* pDataOut = pX->getData;
	CFloat32VolumeData2D* pSlice;
	const float32* pDataIn;
	
	int rows = pX->getHeight();
	int cols = pX->getWidth();

	// First slice: hirozontal difference. -ve in place and +ve to the left
	//
	pSlice = pKx->fetchSliceZ(0);
	pDataIn = pSlice->getDataConst();
	// subtract the slice first
	*pX -= *pSlice;
	for (int i = 0; i < rows; ++i) {
		for (int j = 1; j < cols; ++j) {
			int indx = i * cols + j;
			// add the value on the left
			pDataOut[indx] += pDataIn[indx - 1];
		}
	}

	// Second slice: vertical difference. -ve in place and +ve to the top
	//
	pSlice = pKx->fetchSliceZ(1);
	pDataIn = pSlice->getDataConst();
	// subtract the slice first
	*pX -= *pSlice;
	for (int i = 1; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			int indx = i * cols + j;
			// add the value on the top
			pDataOut[indx] += pDataIn[indx - cols];
		}
	}

	// Multiply by sigma
	if (sigma != 1.f) {
		*pX *= sigma;
	}
}

void CTRexPriorATV::prox(const CFloat32VolumeData3D* pU, float32 rho,
						 CFloat32VolumeData3D* pV)
{
	ASTRA_ASSERT(pU->getDepth() == pV->getDepth() &&
		pU->getSize() == pV->getSize());

	// loop on slices
	for (int k = 0; k < pU->getDepth(); ++k) {
		// Get  slice
		CFloat32VolumeData2D* pSliceIn = pU->fetchSliceZ(k);
		CFloat32VolumeData2D* pSliceOut = pV->fetchSliceZ(k);
		const float32* pDataIn = pSliceIn->getDataConst();
		float32* pDataOut = pSliceOut->getData();
		int sz = pSliceIn->getSize();

		// loop on slice
		for (int i = 0; i < sz; ++i) {
			pDataOut[i] = max(0.f, pDataIn[i] - rho) - 
				max(0.f, -pDataIn[i] - rho);
		}
		// Put back slice
		pV->returnSliceZ(k, pSliceOut);
	}

}

float32 CTRexPriorATV::norm(float32 sigma)
{
	return 8.f * sigma * sigma;
}

#include "astra/Projector2DImpl.inl"

// type of the algorithm, needed to register with CAlgorithmFactory
std::string CTRexAlgorithm::type = "TREX";


//---------------------------------------------------------------------------------------
// Clear - Constructors
void CTRexAlgorithm::_clear()
{
	//CReconstructionAlgorithm2D::_clear();
	CSartAlgorithm::_clear();
	m_fLambda = 1.0f;
	m_pC = NULL;
	m_pY = NULL;
	m_pW = NULL;
}

//---------------------------------------------------------------------------------------
// Clear - Public
void CTRexAlgorithm::clear()
{
	//CReconstructionAlgorithm2D::clear();
	CSartAlgorithm::clear();
	//if (m_piProjectionOrder) {
	//	delete[] m_piProjectionOrder;
	//	m_piProjectionOrder = NULL;
	//}
	//m_iProjectionCount = 0;
	//m_iCurrentProjection = 0;
	//m_bIsInitialized = false;
	//m_iIterationCount = 0;
	ASTRA_DELETE(m_pY);
	ASTRA_DELETE(m_pC);
}

//----------------------------------------------------------------------------------------
// Constructor
CTRexAlgorithm::CTRexAlgorithm() 
{
	_clear();
}

//----------------------------------------------------------------------------------------
// Constructor
CTRexAlgorithm::CTRexAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction) 
{
	_clear();
	CSartAlgorithm::CSartAlgorithm(_pProjector, _pSinogram, _pReconstruction);
	//initialize(_pProjector, _pSinogram, _pReconstruction);
}

//----------------------------------------------------------------------------------------
// Constructor
CTRexAlgorithm::CTRexAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction,
							   int* _piProjectionOrder, 
							   int _iProjectionCount)
{
	_clear();
	CSartAlgorithm::CSartAlgorithm(_pProjector, _pSinogram, _pReconstruction, 
		_piProjectionOrder, _iProjectionCount);
	//initialize(_pProjector, _pSinogram, _pReconstruction, _piProjectionOrder, _iProjectionCount);
}

//----------------------------------------------------------------------------------------
// Destructor
CTRexAlgorithm::~CTRexAlgorithm() 
{
	clear();
}

//---------------------------------------------------------------------------------------
// Initialize - Config
bool CTRexAlgorithm::initialize(const Config& _cfg)
{
	assert(_cfg.self);
	ConfigStackCheck<CAlgorithm> CC("SartProxOperator", this, _cfg);
	
	// if already initialized, clear first
	if (m_bIsInitialized) {
		clear();
	}

	// initialization of parent class
	if (!CSartAlgorithm::initialize(_cfg)) {
		return false;
	}

	//// projection order
	//m_iCurrentProjection = 0;
	//m_iProjectionCount = m_pProjector->getProjectionGeometry()->getProjectionAngleCount();
	//string projOrder = _cfg.self.getOption("ProjectionOrder", "sequential");
	//CC.markOptionParsed("ProjectionOrder");
	//if (projOrder == "sequential") {
	//	m_piProjectionOrder = new int[m_iProjectionCount];
	//	for (int i = 0; i < m_iProjectionCount; i++) {
	//		m_piProjectionOrder[i] = i;
	//	}
	//} else if (projOrder == "random") {
	//	m_piProjectionOrder = new int[m_iProjectionCount];
	//	for (int i = 0; i < m_iProjectionCount; i++) {
	//		m_piProjectionOrder[i] = i;
	//	}
	//	for (int i = 0; i < m_iProjectionCount-1; i++) {
	//		int k = (rand() % (m_iProjectionCount - i));
	//		int t = m_piProjectionOrder[i];
	//		m_piProjectionOrder[i] = m_piProjectionOrder[i + k];
	//		m_piProjectionOrder[i + k] = t;
	//	}
	//} else if (projOrder == "custom") {
	//	vector<float32> projOrderList = _cfg.self.getOptionNumericalArray("ProjectionOrderList");
	//	m_piProjectionOrder = new int[projOrderList.size()];
	//	for (int i = 0; i < m_iProjectionCount; i++) {
	//		m_piProjectionOrder[i] = static_cast<int>(projOrderList[i]);
	//	}
	//	CC.markOptionParsed("ProjectionOrderList");
	//}

	// Lambda
	m_fLambda = _cfg.self.getOptionNumerical("Lambda", m_fLambda);
	CC.markOptionParsed("Lambda");

	//// Alpha
	//m_fAlpha = _cfg.self.getOptionNumerical("Alpha", m_fAlpha);
	//CC.markOptionParsed("Alpha");

	//// Clear RaySum after each sweep. Defaults to true.
	//m_bClearRayLength = _cfg.self.getOptionBool("ClearRayLength", m_bClearRayLength);
	//CC.markOptionParsed("ClearRayLength");

	// Input volume
	XMLNode node = _cfg.self.getSingleNode("ProxInputDataId");
	ASTRA_CONFIG_CHECK(node, "SartProxOperator", "No Proximal Input tag specified.");
	int id = boost::lexical_cast<int>(node.getContent());
	m_pProxInput = dynamic_cast<CFloat32VolumeData2D*>(CData2DManager::getSingleton().get(id));
	CC.markNodeParsed("ProxInputDataId");

	// WLS weights
	id = static_cast<int>(_cfg.self.getOptionNumerical("WlsWeightDataId", -1));
	m_pW = dynamic_cast<CFloat32ProjectionData2D*>(CData2DManager::getSingleton().get(id));
	CC.markOptionParsed("WlsWeightDataId");

	// create data objects
	//m_pTotalRayLength = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	//m_pTotalPixelWeight = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	//m_pDiffSinogram = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pY = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pC = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}


//----------------------------------------------------------------------------------------
bool CTRexAlgorithm::_check()
{
	// check base class
	ASTRA_CONFIG_CHECK(CSartAlgorithm::_check(), "SartProxOperator",
		"Error in ReconstructionAlgorithm2D initialization");

	return true;
}

//---------------------------------------------------------------------------------------
// Information - All
map<string,boost::any> CTRexAlgorithm::getInformation() 
{
	map<string, boost::any> res;
	res["ProjectionOrder"] = getInformation("ProjectionOrder");
	return mergeMap<string,boost::any>(CReconstructionAlgorithm2D::getInformation(), res);
};

//---------------------------------------------------------------------------------------
// Information - Specific
boost::any CTRexAlgorithm::getInformation(std::string _sIdentifier) 
{
	if (_sIdentifier == "ProjectionOrder") {
		vector<float32> res;
		for (int i = 0; i < m_iProjectionCount; i++) {
			res.push_back(m_piProjectionOrder[i]);
		}
		return res;
	}
	return CAlgorithm::getInformation(_sIdentifier);
};

//----------------------------------------------------------------------------------------
// Iterate
void CTRexAlgorithm::run(int _iNrIterations)
{
	// check initialized
	ASTRA_ASSERT(m_bIsInitialized);

	m_bShouldAbort = false;

	// Clear reconstruction (not using the input initial value).
	if (m_bClearReconstruction) {
		m_pReconstruction->setData(0.f);
	}
	// alias
	CFloat32VolumeData2D* pX = m_pReconstruction;

	// Initialiazation
	CVolumeGeometry3D geom(m_pReconstruction->getWidth(), 
		m_pReconstruction->getHeight(), m_pPrior->depth());
	CFloat32VolumeData3D* pZ = new CFloat32VolumeData3DMemory(&geom);
	CFloat32VolumeData3D* pU = new CFloat32VolumeData3DMemory(&geom);
	CFloat32VolumeData3D* pKx = new CFloat32VolumeData3DMemory(&geom);
	CFloat32VolumeData2D* pT = new CFloat32VolumeData2D(pX->getGeometry());

	// Outer loop
	for (int iIteration = 0; iIteration < _iNrIterations; ++iIteration) {
		// start timer
		m_ulTimer = CPlatformDepSystemCode::getMSCount();

		// x-step: data term
		//
		// Kx - z + u
		m_pPrior->K(pX, pKx, m_fSigma);
		*pKx -= *pZ;
		*pKx += *pU;
		// t = Kt * Kx
		m_pPrior->Kt(pKx, pT, m_fSigma);
		// t = -rho*mu*t
		*pT *= -m_fRho * m_fMu;
		// t = x - rho*mu * t
		*pT += *pX;


		// z-step: prior
		//
		// Kx + u
		m_pPrior->K(pX, pKx, m_fSigma);
		*pKx += *pU;
		// z = prox(Kx + u)
		m_pPrior->prox(pKx, m_fSigma / m_fRho, pZ);

		// u-step: dual variable
		//
		// u = u + Kx - z
		m_pPrior->K(pX, pKx, m_fSigma);
		*pU += *pKx;
		*pU -= *pZ;

		// end timer
		m_ulTimer = CPlatformDepSystemCode::getMSCount() - m_ulTimer;

		// Compute metrics.
		computeIterationMetrics(iIteration, _iNrIterations, m_pDiffSinogram);
	}

	// Clear
	ASTRA_DELETE(pZ);
	ASTRA_DELETE(pU);
	ASTRA_DELETE(pKx);
	ASTRA_DELETE(pT);

	//for (int i=0; i < m_pReconstruction->getSize(); ++i) {
	//	ASTRA_INFO("voxel=%d val=%f", i, m_pReconstruction->getData()[i]);
	//}
}
//----------------------------------------------------------------------------------------

} // namespace astra
