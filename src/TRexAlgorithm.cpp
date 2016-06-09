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

#include <iostream>
#include <boost/lexical_cast.hpp>

#include "astra/AstraObjectFactory.h"
#include "astra/AstraObjectManager.h"
#include "astra/DataProjectorPolicies.h"
#include "astra/Logging.h"
#include "astra/Float32VolumeData3DMemory.h"

using namespace std;

namespace astra {

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
	//m_fWlsRoot = 1;
	m_iInnterIter = 2;
	m_pPrior = NULL;
	m_pData = NULL;
	//m_pProxMetrics = NULL;
	//m_pProxInput = NULL;
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
	//ASTRA_DELETE(m_pProxInput);
	//ASTRA_DELETE(m_pProxMetrics);

	// Delete these from the data manager since they are stored there
	//CData2DManager::getSingleton().remove(
	//	CData2DManager::getSingleton().getIndex(m_pProxInput));
	//m_pProxInput = NULL;
	//CData2DManager::getSingleton().remove(
	//	CData2DManager::getSingleton().getIndex(m_pProxMetrics));
	//m_pProxMetrics = NULL;
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

	int id;

	// Sigma
	m_fSigma = _cfg.self.getOptionNumerical("Sigma", m_fSigma);
	CC.markOptionParsed("Sigma");
	ASTRA_CONFIG_CHECK(m_fSigma > 0, "CTRexAlgorithm", 
		"Error initializing: Sigma <= 0");

	// Prior
	string sPrior = _cfg.self.getOption("Prior","");
	if (sPrior == "ATV")  {
		m_pPrior = new CTRexPriorATV(m_fSigma);
	} else if (sPrior == "ITV") {
		m_pPrior = new CTRexPriorITV(m_fSigma);
	} else if (sPrior == "SAD") {
		m_pPrior = new CTRexPriorSAD(m_fSigma);
	}
	CC.markOptionParsed("Prior");
	ASTRA_CONFIG_CHECK(m_pPrior, "CTRexAlgorithm", 
		"Error initializing: unimplemented prior");

	// Rho and Mu
	m_fRho = _cfg.self.getOptionNumerical("Rho", m_fRho);
	CC.markOptionParsed("Rho");
	ASTRA_CONFIG_CHECK(m_fRho > 0, "CTRexAlgorithm", 
		"Error initializing: Rho <= 0");

	m_fMu = _cfg.self.getOptionNumerical("Mu", m_fMu);
	CC.markOptionParsed("Mu");
	if (m_fMu <= 0) {
		m_fMu = 1.f / (m_fRho * m_pPrior->norm());
	}	

	// Inner iterations
	m_iInnterIter = static_cast<int>(
		_cfg.self.getOptionNumerical("InnerIter", 
									 static_cast<float32>(m_iInnterIter)));
	CC.markOptionParsed("InnerIter");
	ASTRA_CONFIG_CHECK(m_iInnterIter > 0, "CTRexAlgorithm", 
		"Error initializing: InnerIter <= 0");

	//// Debugging
	//Config c = Config(_cfg);
	//ASTRA_INFO("%s\n", c.self.toString().c_str());
	//c.self.removeOption("Alpha");
	//c.self.addOption("Alpha", 1.234);
	//c.self.removeChildNode("ProxInputDataId");
	//c.self.addChildNode("ProxInputDataId", 1234);
	//ASTRA_INFO("%s\n", c.self.toString().c_str());
	//return false;

	// Data term
	string sData = _cfg.self.getOption("Data","");
	if (sData == "LS")  {
		m_pData = new CTRexDataLS();
	} else if (sData == "WLS") {
		m_pData = new CTRexDataWLS();
	}
	CC.markOptionParsed("Data");
	ASTRA_CONFIG_CHECK(m_pData, "CTRexAlgorithm", 
		"Error initializing: unimplemented Data term.");

	// Parse the input and init the data term
	ASTRA_CONFIG_CHECK(
		m_pData->init(_cfg, m_pReconstruction->getGeometry(), m_fMu),
		"CTRexAlgorithm", "Error initializing: Data term");

	// -> Moved to Data class
	// Get the data proximal operator algorithm
	//string sProx = _cfg.self.getOption("DataProx","");
	//m_pData->m_pAlg = dynamic_cast<CReconstructionAlgorithm2D*>(
	//	CAlgorithmFactory::getSingleton().create(sProx));
	//ASTRA_CONFIG_CHECK(m_pData->m_pAlg, "CTRexAlgorithm", 
	//	"Error initializing: unknown DataProx.");
	//CC.markOptionParsed("DataProx");

	//// Get the config for the Data Prox as a copy of this config
	//m_pData->m_pCfg = new Config(_cfg);
	////ASTRA_INFO("%s\n", m_pData->m_pCfg->self.toString().c_str());

	//// Compute metrics?
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

	//// Add ProxInput to the prox config
	//// Create object and add to manager
	//m_pProxInput = new 
	//	CFloat32VolumeData2D(m_pReconstruction->getGeometry());
	//id = CData2DManager::getSingleton().store(m_pProxInput);
	//// Add the id to config
	//m_pData->m_pCfg->self.removeChildNode("ProxInputDataId");
	//m_pData->m_pCfg->self.addChildNode("ProxInputDataId", 
	//	static_cast<float32>(id));

	//// Set Lambda to Mu
	//m_pData->m_pCfg->self.removeOption("Lambda");
	//m_pData->m_pCfg->self.addOption("Lambda", m_fMu);

	//// Initialize the prox algorithm
	//m_pData->m_pAlg->initialize(*m_pData->m_pCfg);

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
	ASTRA_ASSERT(pX);

	// Initialiazation
	CVolumeGeometry3D geom(m_pReconstruction->getWidth(), 
		m_pReconstruction->getHeight(), m_pPrior->depth());
	CFloat32VolumeData3DMemory* pZ = new CFloat32VolumeData3DMemory(&geom);
	CFloat32VolumeData3DMemory* pU = new CFloat32VolumeData3DMemory(&geom);
	CFloat32VolumeData3DMemory* pKx = new CFloat32VolumeData3DMemory(&geom);
	// alias to ProxInput
	CFloat32VolumeData2D* pT = m_pData->m_pProxInput;
	ASTRA_ASSERT(pT);
	//CFloat32VolumeData2D* pT = new CFloat32VolumeData2D(pX->getGeometry());

	// Initialize
	m_pPrior->K(pX, pZ);
	m_pPrior->K(pX, pU);

	// Outer loop
	for (int iIteration = 0; iIteration < _iNrIterations; ++iIteration) {
		// start timer
		m_ulTimer = CPlatformDepSystemCode::getMSCount();
		//ASTRA_INFO("Iteration: %d", iIteration);

		// x-step: data term
		//
		//ASTRA_INFO("X-step");
		// Kx - z + u
		m_pPrior->K(pX, pKx);
		*pKx -= *pZ;
		*pKx += *pU;
		// t = Kt * Kx
		m_pPrior->Kt(pKx, pT);
		// t = -rho*mu*t
		*pT *= -m_fRho * m_fMu;
		// t = x - rho*mu * t
		*pT += *pX;
		// run the data proximal operator: input in pT (ProxInput) and
		// output in pX (Reconstruction)
		//m_pData->m_pAlg->run(m_iInnterIter);
		m_pData->prox(m_iInnterIter);
		// Clamp
		if (m_bUseMinConstraint)
			pX->clampMin(m_fMinValue);
		//pX->printInfo("X");

		// z-step: prior
		//
		//ASTRA_INFO("Z-step");
		// Kx + u
		//pZ->printInfo("Z before");
		m_pPrior->K(pX, pKx);
		*pKx += *pU;
		// z = prox(Kx + u)
		m_pPrior->prox(pKx, m_fSigma / m_fRho, pZ);
		//pZ->printInfo("Z");

		// u-step: dual variable
		//
		//ASTRA_INFO("U-step");
		// u = u + Kx - z
		//pU->printInfo("U before");
		//pU->opera *pKx;
		//pU->operator=(*pKx);
		pU->copyData(pKx->getDataConst());
		*pU -= *pZ;
		//m_pPrior->K(pX, pKx, m_fSigma);
		//*pU += *pKx;
		//*pU -= *pZ;
		//pU->printInfo("U");

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
