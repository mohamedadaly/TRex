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

#include "astra/AstraObjectFactory.h"
#include "astra/AstraObjectManager.h"
#include "astra/DataProjectorPolicies.h"
#include "astra/Logging.h"
#include "astra/Float32VolumeData3DMemory.h"
#include "astra/ArtProxOperatorAlgorithm.h"
#include "astra/BicavProxOperatorAlgorithm.h"
#include "astra/SartProxOperatorAlgorithm.h"
#include "astra/OrderedSubsetSQSProxOperatorAlgorithm.h"

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
	float32* pDataOut = pX->getData();
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


// ----------------------------------------------------------------------------
void CTRexDataLS::init() 
{

}
void CTRexDataLS::prox(const CFloat32VolumeData2D* pU, CFloat32VolumeData2D* pX)
{

}

// ----------------------------------------------------------------------------

#include "astra/Projector2DImpl.inl"

// type of the algorithm, needed to register with CAlgorithmFactory
std::string CTRexAlgorithm::type = "TREX";


//---------------------------------------------------------------------------------------
// Clear - Constructors
void CTRexAlgorithm::_clear()
{
	CReconstructionAlgorithm2D::_clear();
	//CSartAlgorithm::_clear();
	m_fSigma = 1.f;
	m_fRho = 100.f;
	m_fMu = -1.f;
	m_fWlsRoot = 1;
	m_iInnterIter = 2;
	m_pPrior = NULL;
	m_pData = NULL;
	m_pProxMetrics = NULL;
	m_pProxInput = NULL;
	//m_pDataProxOperator = NULL;
	//m_pTomoProxOperatorConfig = NULL;
}

//---------------------------------------------------------------------------------------
// Clear - Public
void CTRexAlgorithm::clear()
{
	CReconstructionAlgorithm2D::clear();
	//CSartAlgorithm::clear();
	//if (m_piProjectionOrder) {
	//	delete[] m_piProjectionOrder;
	//	m_piProjectionOrder = NULL;
	//}
	//m_iProjectionCount = 0;
	//m_iCurrentProjection = 0;
	//m_bIsInitialized = false;
	//m_iIterationCount = 0;
	ASTRA_DELETE(m_pPrior);
	ASTRA_DELETE(m_pData);
	ASTRA_DELETE(m_pProxInput);
	ASTRA_DELETE(m_pProxMetrics);
	//CData2DManager::getSingleton().remove(m_iProxMetricsId);
	//CData2DManager::getSingleton().remove(m_iProxInputId);
	//ASTRA_DELETE(m_pDataProxOperator);
	//ASTRA_DELETE(m_pTomoProxOperatorConfig);
	
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
	//initialize(_pProjector, _pSinogram, _pReconstruction);
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
	ConfigStackCheck<CAlgorithm> CC("TRexAlgorithm", this, _cfg);
	
	// if already initialized, clear first
	if (m_bIsInitialized) {
		clear();
	}

	// initialization of parent class which reads SART's options
	if (!CReconstructionAlgorithm2D::initialize(_cfg)) {
		return false;
	}

	// Prior
	string sPrior = _cfg.self.getOption("Prior","");
	if (sPrior == "ATV")  {
		m_pPrior = new CTRexPriorATV();
	//} else if (sPrior == "ITV") {
	//	m_pPrior = new CTRexPriorITV();
	//} else if (sPrior == "SAD") {
	//	m_pPrior = new CTRexPriorSAD();
	}
	CC.markOptionParsed("Prior");
	ASTRA_CONFIG_CHECK(m_pPrior, "CTRexAlgorithm", 
		"Error initializing: unimplemented prior");

	// rho, sigma, mu
	m_fSigma = _cfg.self.getOptionNumerical("Sigma", m_fSigma);
	CC.markOptionParsed("Sigma");
	ASTRA_CONFIG_CHECK(m_fSigma > 0, "CTRexAlgorithm", 
		"Error initializing: Sigma <= 0");

	m_fRho = _cfg.self.getOptionNumerical("Rho", m_fRho);
	CC.markOptionParsed("Rho");
	ASTRA_CONFIG_CHECK(m_fRho > 0, "CTRexAlgorithm", 
		"Error initializing: Rho <= 0");

	m_fMu = _cfg.self.getOptionNumerical("Mu", m_fMu);
	CC.markOptionParsed("Mu");
	if (m_fMu <= 0) {
		m_fMu = 1.f / (m_fRho * m_pPrior->norm(m_fSigma));
	}

	// WLS root
	m_fWlsRoot = _cfg.self.getOptionNumerical("WlsRoot", m_fWlsRoot);
	CC.markOptionParsed("WlsRoot");
	ASTRA_CONFIG_CHECK(m_fWlsRoot > 0, "CTRexAlgorithm", 
		"Error initializing: WlsRoot <= 0");

	// Inner iterations
	m_iInnterIter = static_cast<int>(
		_cfg.self.getOptionNumerical("InnerIter", m_iInnterIter));
	CC.markOptionParsed("InnerIter");
	ASTRA_CONFIG_CHECK(m_iInnterIter > 0, "CTRexAlgorithm", 
		"Error initializing: InnerIter <= 0");
	
	// Data term
	string sData = _cfg.self.getOption("Data","");
	if (sData == "LS")  {
		m_pData = new CTRexDataLS();
	//} else if (sData == "WLS") {
	//	m_pData = new CTRexDataWLS();
	}
	CC.markOptionParsed("Data");
	ASTRA_CONFIG_CHECK(m_pData, "CTRexAlgorithm", 
		"Error initializing: unimplemented Data term.");

	// Get the data proximal operator algorithm
	string sProx = _cfg.self.getOption("DataProx","");
	m_pData->m_pAlg = dynamic_cast<CReconstructionAlgorithm2D*>(
		CAlgorithmFactory::getSingleton().create(sProx));
	ASTRA_CONFIG_CHECK(m_pData->m_pAlg, "CTRexAlgorithm", 
		"Error initializing: unknown DataProx.");
	CC.markOptionParsed("DataProx");

	// Get the config for the Data Prox as a copy of this config
	m_pData->m_pCfg = new Config(_cfg);

	// Compute metrics?
	if (m_bComputeIterationMetrics) {
		// Create a matrix for that
		m_pProxMetrics = new CFloat32VolumeData2D();
		// Store in manager and save Id to pass to prox solver
		int id = CData2DManager::getSingleton().store(m_pProxMetrics);
		m_pData->m_pCfg->self.getSingleNode("IterationMetricsId").
			setContent(static_cast<float32>(id));
	}

	// Add ProxInput to the prox config
	{
		// Create object and add to manager
		CFloat32VolumeData2D* m_pProxInput = new 
			CFloat32VolumeData2D(m_pReconstruction->getGeometry());
		int id = CData2DManager::getSingleton().store(m_pProxInput);
		// Add the id to config
		XMLNode node = _cfg.self.getSingleNode("ProxInputDataId");
		if (node) {
			node.setContent(static_cast<float32>(id));
		} else {
			m_pData->m_pCfg->self.addChildNode("ProxInputDataId", 
				static_cast<float32>(id));
		}
	}

	// Initialize the prox algorithm
	m_pData->m_pAlg->initialize(*m_pData->m_pCfg);

	//XMLNode node = _cfg.self.getSingleNode("ProxInputDataId");
	//ASTRA_CONFIG_CHECK(node, "SartProxOperator", "No Proximal Input tag specified.");
	//int id = boost::lexical_cast<int>(node.getContent());
	//m_pProxInput = dynamic_cast<CFloat32VolumeData2D*>(CData2DManager::getSingleton().get(id));
	//CC.markNodeParsed("ProxInputDataId");

	//// WLS weights
	//int id = static_cast<int>(_cfg.self.getOptionNumerical("WlsWeightDataId", -1));
	//m_pW = dynamic_cast<CFloat32ProjectionData2D*>(CData2DManager::getSingleton().get(id));
	//CC.markOptionParsed("WlsWeightDataId");


	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}


//----------------------------------------------------------------------------------------
bool CTRexAlgorithm::_check()
{
	// check base class
	ASTRA_CONFIG_CHECK(CReconstructionAlgorithm2D::_check(), "TRexAlgorithm",
		"Error in TRexAlgorithm initialization");

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
	//if (_sIdentifier == "ProjectionOrder") {
	//	vector<float32> res;
	//	for (int i = 0; i < m_iProjectionCount; i++) {
	//		res.push_back(m_piProjectionOrder[i]);
	//	}
	//	return res;
	//}
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
	// alias to ProxInput
	CFloat32VolumeData2D* pT = m_pProxInput;
	//CFloat32VolumeData2D* pT = new CFloat32VolumeData2D(pX->getGeometry());

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
		// run the data proximal operator: input in pT (ProxInput) and
		// output in pX (Reconstruction)
		m_pData->m_pAlg->run(m_iInnterIter);

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
		computeIterationMetrics(iIteration, _iNrIterations, NULL);

		// Get the inner iteration metrics
		// clear the inner iteration metrics		

	}

	// Clear
	ASTRA_DELETE(pZ);
	ASTRA_DELETE(pU);
	ASTRA_DELETE(pKx);
	//ASTRA_DELETE(pT);

	//for (int i=0; i < m_pReconstruction->getSize(); ++i) {
	//	ASTRA_INFO("voxel=%d val=%f", i, m_pReconstruction->getData()[i]);
	//}
}
//----------------------------------------------------------------------------------------

} // namespace astra
