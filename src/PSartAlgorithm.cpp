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

#include "astra/PSartAlgorithm.h"

#include <boost/lexical_cast.hpp>

#include "astra/AstraObjectManager.h"
#include "astra/DataProjectorPolicies.h"
#include "astra/Logging.h"


using namespace std;

namespace astra {

#include "astra/Projector2DImpl.inl"

// type of the algorithm, needed to register with CAlgorithmFactory
std::string CPSartAlgorithm::type = "PSART";


//---------------------------------------------------------------------------------------
// Clear - Constructors
void CPSartAlgorithm::_clear()
{
	//CReconstructionAlgorithm2D::_clear();
	CSartAlgorithm::_clear();
	m_fLambda = 1.0f;
}

//---------------------------------------------------------------------------------------
// Clear - Public
void CPSartAlgorithm::clear()
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
}

//----------------------------------------------------------------------------------------
// Constructor
CPSartAlgorithm::CPSartAlgorithm() 
{
	_clear();
}

//----------------------------------------------------------------------------------------
// Constructor
CPSartAlgorithm::CPSartAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction) 
{
	CSartAlgorithm::CSartAlgorithm(_pProjector, _pSinogram, _pReconstruction);
	//_clear();
	//initialize(_pProjector, _pSinogram, _pReconstruction);
}

//----------------------------------------------------------------------------------------
// Constructor
CPSartAlgorithm::CPSartAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction,
							   int* _piProjectionOrder, 
							   int _iProjectionCount)
{
	CSartAlgorithm::CSartAlgorithm(_pProjector, _pSinogram, _pReconstruction, 
		_piProjectionOrder, _iProjectionCount);
	//_clear();
	//initialize(_pProjector, _pSinogram, _pReconstruction, _piProjectionOrder, _iProjectionCount);
}

//----------------------------------------------------------------------------------------
// Destructor
CPSartAlgorithm::~CPSartAlgorithm() 
{
	clear();
}

//---------------------------------------------------------------------------------------
// Initialize - Config
bool CPSartAlgorithm::initialize(const Config& _cfg)
{
	assert(_cfg.self);
	ConfigStackCheck<CAlgorithm> CC("PSartAlgorithm", this, _cfg);
	
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
	ASTRA_CONFIG_CHECK(node, "PSART", "No Proximal Input tag specified.");
	int id = boost::lexical_cast<int>(node.getContent());
	m_pProxInput = dynamic_cast<CFloat32VolumeData2D*>(CData2DManager::getSingleton().get(id));
	CC.markNodeParsed("ProxInputDataId");

	// create data objects
	//m_pTotalRayLength = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	//m_pTotalPixelWeight = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	//m_pDiffSinogram = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pY = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}


//----------------------------------------------------------------------------------------
bool CPSartAlgorithm::_check()
{
	// check base class
	ASTRA_CONFIG_CHECK(CSartAlgorithm::_check(), "PSART", "Error in ReconstructionAlgorithm2D initialization");

	return true;
}

//---------------------------------------------------------------------------------------
// Information - All
map<string,boost::any> CPSartAlgorithm::getInformation() 
{
	map<string, boost::any> res;
	res["ProjectionOrder"] = getInformation("ProjectionOrder");
	return mergeMap<string,boost::any>(CReconstructionAlgorithm2D::getInformation(), res);
};

//---------------------------------------------------------------------------------------
// Information - Specific
boost::any CPSartAlgorithm::getInformation(std::string _sIdentifier) 
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
void CPSartAlgorithm::run(int _iNrIterations)
{
	// check initialized
	ASTRA_ASSERT(m_bIsInitialized);

	m_bShouldAbort = false;

	// data projectors
	CDataProjectorInterface* pForwardProjector;
	CDataProjectorInterface* pBackProjector;
	CDataProjectorInterface* pFirstForwardProjector;


	// Init Y = 0
	m_pY->setData(0.f);
	// Init Reconstruction = ProxInput
	m_pReconstruction->copyData(m_pProxInput->getData());
	m_pReconstruction->updateStatistics();
	//ASTRA_INFO("Initialized from ProxInput: max=%f min=%f id=%d", 
	//	m_pReconstruction->getGlobalMax(), m_pReconstruction->getGlobalMin(),
	//	CData2DManager::getSingleton().getIndex(m_pReconstruction));
	//for (int i=0; i < m_pReconstruction->getSize(); ++i) {
	//	ASTRA_INFO("voxel=%d val=%f", i, m_pReconstruction->getData()[i]);
	//}	

	// constant
	float32 fSqrt2Lambda = sqrtf(m_fLambda * 2.f);
	//ASTRA_INFO("Sqrt2Lambda = %f", fSqrt2Lambda);
	//ASTRA_INFO("Alpha = %f", m_fAlpha);

	//ASTRA_INFO("UseMinConst=%d MinVal=%f", m_bUseMinConstraint, m_fMinValue);
	//ASTRA_INFO("UseMaxConst=%d MaxVal=%f", m_bUseMaxConstraint, m_fMaxValue);
	// Scale projections by sqrt(2 * lambda)
	// m_pSinogram->operator*=(fSqrt_2_lambda);

	// backprojection data projector
	pBackProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			PSARTBPPolicy(m_pReconstruction, m_pDiffSinogram, 
						 m_pTotalPixelWeight, m_pTotalRayLength, 
						 m_pY, m_fAlpha, fSqrt2Lambda),			// PSART backprojection
			m_bUseSinogramMask, m_bUseReconstructionMask, true // options on/off
		); 

	// first time forward projection data projector,
	// also computes total pixel weight and total ray length
	pForwardProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			CombinePolicy<DiffFPPolicy, TotalPixelWeightPolicy>(					// 3 basic operations
				DiffFPPolicy(m_pReconstruction, m_pDiffSinogram, m_pSinogram),								// forward projection with difference calculation
				TotalPixelWeightPolicy(m_pTotalPixelWeight)),												// calculate the total pixel weights
			m_bUseSinogramMask, m_bUseReconstructionMask, true											 // options on/off
		);

	// first time forward projection data projector,
	// computes total ray length
	pFirstForwardProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			TotalRayLengthPolicy(m_pTotalRayLength),													// calculate the total ray lengths
			m_bUseSinogramMask, m_bUseReconstructionMask, true											 // options on/off
		);

	// Perform the first forward projection to compute ray lengths
	m_pTotalRayLength->setData(0.0f);
	m_pTotalPixelWeight->setData(0.0f);
	pFirstForwardProjector->project();

	// iteration loop, each iteration loops over all available projections
	for (int iIteration = 0; iIteration < _iNrIterations && !m_bShouldAbort; ++iIteration) {
		//ASTRA_INFO("Iteration %d", iIteration);
		// Clear RayLength before another loop over projections. This is needed so that
		// RayLength is correct, because updating RayLength with the forward projection
		// again will multiply the RayLength when processing the same ray in the next
		// iteration.
		//if (m_bClearRayLength) {
		//	m_pTotalRayLength->setData(0.f);
		//}

		// start timer
		m_ulTimer = CPlatformDepSystemCode::getMSCount();

		// loop over projections
		for (int iP = 0; iP < m_iProjectionCount; ++iP) {
			// projection id
			// int iProjection = m_piProjectionOrder[m_iIterationCount % m_iProjectionCount];
			int iProjection = m_piProjectionOrder[iP % m_iProjectionCount];
			//ASTRA_INFO(" Projection %d", iProjection);

			// forward projection and difference calculation
			m_pTotalPixelWeight->setData(0.0f);
			pForwardProjector->projectSingleProjection(iProjection);
			// backprojection
			pBackProjector->projectSingleProjection(iProjection);
			// update iteration count
			m_iIterationCount++;

			if (m_bUseMinConstraint)
				m_pReconstruction->clampMin(m_fMinValue);
			if (m_bUseMaxConstraint)
				m_pReconstruction->clampMax(m_fMaxValue);
		}

		// end timer
		m_ulTimer = CPlatformDepSystemCode::getMSCount() - m_ulTimer;

		// Compute metrics.
		computeIterationMetrics(iIteration, _iNrIterations);
	}

	ASTRA_DELETE(pForwardProjector);
	ASTRA_DELETE(pBackProjector);
	ASTRA_DELETE(pFirstForwardProjector);

	//for (int i=0; i < m_pReconstruction->getSize(); ++i) {
	//	ASTRA_INFO("voxel=%d val=%f", i, m_pReconstruction->getData()[i]);
	//}
}
//----------------------------------------------------------------------------------------

} // namespace astra
