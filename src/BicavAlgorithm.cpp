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

#include "astra/BicavAlgorithm.h"

#include <boost/lexical_cast.hpp>

#include "astra/AstraObjectManager.h"
#include "astra/DataProjectorPolicies.h"

using namespace std;

namespace astra {

#include "astra/Projector2DImpl.inl"

// type of the algorithm, needed to register with CAlgorithmFactory
std::string CBicavAlgorithm::type = "BICAV";


//---------------------------------------------------------------------------------------
// Clear - Constructors
void CBicavAlgorithm::_clear()
{
	CReconstructionAlgorithm2D::_clear();
	m_piProjectionOrder = NULL;
	m_iProjectionCount = 0;
	m_iCurrentProjection = 0;
	m_bIsInitialized = false;
	m_iIterationCount = 0;
	m_fAlpha = 1.0f;
	m_bClearRayLength = true;
}

//---------------------------------------------------------------------------------------
// Clear - Public
void CBicavAlgorithm::clear()
{
	CReconstructionAlgorithm2D::clear();
	if (m_piProjectionOrder) {
		delete[] m_piProjectionOrder;
		m_piProjectionOrder = NULL;
	}
	m_iProjectionCount = 0;
	m_iCurrentProjection = 0;
	m_bIsInitialized = false;
	m_iIterationCount = 0;
}

//----------------------------------------------------------------------------------------
// Constructor
CBicavAlgorithm::CBicavAlgorithm() 
{
	_clear();
}

//----------------------------------------------------------------------------------------
// Constructor
CBicavAlgorithm::CBicavAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction) 
{
	_clear();
	initialize(_pProjector, _pSinogram, _pReconstruction);
}

//----------------------------------------------------------------------------------------
// Constructor
CBicavAlgorithm::CBicavAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction,
							   int* _piProjectionOrder, 
							   int _iProjectionCount)
{
	_clear();
	initialize(_pProjector, _pSinogram, _pReconstruction, _piProjectionOrder, _iProjectionCount);
}

//----------------------------------------------------------------------------------------
// Destructor
CBicavAlgorithm::~CBicavAlgorithm() 
{
	clear();
}

//---------------------------------------------------------------------------------------
// Initialize - Config
bool CBicavAlgorithm::initialize(const Config& _cfg)
{
	assert(_cfg.self);
	ConfigStackCheck<CAlgorithm> CC("BicavAlgorithm", this, _cfg);
	
	// if already initialized, clear first
	if (m_bIsInitialized) {
		clear();
	}

	// initialization of parent class
	if (!CReconstructionAlgorithm2D::initialize(_cfg)) {
		return false;
	}

	// projection order
	m_iCurrentProjection = 0;
	m_iProjectionCount = m_pProjector->getProjectionGeometry()->getProjectionAngleCount();
	string projOrder = _cfg.self.getOption("ProjectionOrder", "sequential");
	CC.markOptionParsed("ProjectionOrder");
	if (projOrder == "sequential") {
		m_piProjectionOrder = new int[m_iProjectionCount];
		for (int i = 0; i < m_iProjectionCount; i++) {
			m_piProjectionOrder[i] = i;
		}
	} else if (projOrder == "random") {
		m_piProjectionOrder = new int[m_iProjectionCount];
		for (int i = 0; i < m_iProjectionCount; i++) {
			m_piProjectionOrder[i] = i;
		}
		for (int i = 0; i < m_iProjectionCount-1; i++) {
			int k = (rand() % (m_iProjectionCount - i));
			int t = m_piProjectionOrder[i];
			m_piProjectionOrder[i] = m_piProjectionOrder[i + k];
			m_piProjectionOrder[i + k] = t;
		}
	} else if (projOrder == "custom") {
		vector<float32> projOrderList = _cfg.self.getOptionNumericalArray("ProjectionOrderList");
		m_piProjectionOrder = new int[projOrderList.size()];
		for (int i = 0; i < m_iProjectionCount; i++) {
			m_piProjectionOrder[i] = static_cast<int>(projOrderList[i]);
		}
		CC.markOptionParsed("ProjectionOrderList");
	}

	// Alpha
	m_fAlpha = _cfg.self.getOptionNumerical("Alpha", m_fAlpha);
	CC.markOptionParsed("Alpha");

	// Clear RaySum after each sweep. Defaults to true.
	m_bClearRayLength = _cfg.self.getOptionBool("ClearRayLength", m_bClearRayLength);
	CC.markOptionParsed("ClearRayLength");

	// create data objects
	m_pTotalRayLength = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pTotalPixelWeight = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	m_pDiffSinogram = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}

//---------------------------------------------------------------------------------------
// Initialize - C++
bool CBicavAlgorithm::initialize(CProjector2D* _pProjector, 
								CFloat32ProjectionData2D* _pSinogram, 
								CFloat32VolumeData2D* _pReconstruction)
{
	// if already initialized, clear first
	if (m_bIsInitialized) {
		clear();
	}

	// required classes
	m_pProjector = _pProjector;
	m_pSinogram = _pSinogram;
	m_pReconstruction = _pReconstruction;

	// ray order
	m_iCurrentProjection = 0;
	m_iProjectionCount = _pProjector->getProjectionGeometry()->getProjectionAngleCount();
	m_piProjectionOrder = new int[m_iProjectionCount];
	for (int i = 0; i < m_iProjectionCount; i++) {
		m_piProjectionOrder[i] = i;
	}

	// create data objects
	m_pTotalRayLength = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pTotalPixelWeight = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	m_pDiffSinogram = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}

//---------------------------------------------------------------------------------------
// Initialize - C++
bool CBicavAlgorithm::initialize(CProjector2D* _pProjector, 
								CFloat32ProjectionData2D* _pSinogram, 
								CFloat32VolumeData2D* _pReconstruction,
								int* _piProjectionOrder, 
								int _iProjectionCount)
{
	// required classes
	m_pProjector = _pProjector;
	m_pSinogram = _pSinogram;
	m_pReconstruction = _pReconstruction;

	// ray order
	m_iCurrentProjection = 0;
	m_iProjectionCount = _iProjectionCount;
	m_piProjectionOrder = new int[m_iProjectionCount];
	for (int i = 0; i < m_iProjectionCount; i++) {
		m_piProjectionOrder[i] = _piProjectionOrder[i];
	}

	// create data objects
	m_pTotalRayLength = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pTotalPixelWeight = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	m_pDiffSinogram = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}

//----------------------------------------------------------------------------------------
bool CBicavAlgorithm::_check()
{
	// check base class
	ASTRA_CONFIG_CHECK(CReconstructionAlgorithm2D::_check(), "BICAV", "Error in ReconstructionAlgorithm2D initialization");

	// check projection order all within range
	for (int i = 0; i < m_iProjectionCount; ++i) {
		ASTRA_CONFIG_CHECK(0 <= m_piProjectionOrder[i] && m_piProjectionOrder[i] < m_pProjector->getProjectionGeometry()->getProjectionAngleCount(), "BICAV", "Projection Order out of range.");
	}

	return true;
}

//---------------------------------------------------------------------------------------
// Information - All
map<string,boost::any> CBicavAlgorithm::getInformation() 
{
	map<string, boost::any> res;
	res["ProjectionOrder"] = getInformation("ProjectionOrder");
	return mergeMap<string,boost::any>(CReconstructionAlgorithm2D::getInformation(), res);
};

//---------------------------------------------------------------------------------------
// Information - Specific
boost::any CBicavAlgorithm::getInformation(std::string _sIdentifier) 
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
void CBicavAlgorithm::run(int _iNrIterations)
{
	// check initialized
	ASTRA_ASSERT(m_bIsInitialized);

	m_bShouldAbort = false;

	// data projectors
	CDataProjectorInterface* pForwardProjector;
	CDataProjectorInterface* pBackProjector;
	CDataProjectorInterface* pFirstForwardProjector;

	// Initialize m_pReconstruction to zero.
	m_pReconstruction->setData(0.f);

	// backprojection data projector
	pBackProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			SIRTBPPolicy(m_pReconstruction, m_pDiffSinogram, 
			m_pTotalPixelWeight, m_pTotalRayLength, m_fAlpha),	// SIRT backprojection
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
				TotalPixelWeightPolicy(m_pTotalPixelWeight, true)),												// calculate the total pixel weights
			m_bUseSinogramMask, m_bUseReconstructionMask, true											 // options on/off
		);

	// first time forward projection data projector,
	// computes total ray length
	pFirstForwardProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			TotalRayLengthPolicy(m_pTotalRayLength, true),													// calculate the total ray lengths
			m_bUseSinogramMask, m_bUseReconstructionMask, true											 // options on/off
		);

	// Perform the first forward projection to compute ray lengths
	m_pTotalRayLength->setData(0.0f);
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
	}


	ASTRA_DELETE(pForwardProjector);
	ASTRA_DELETE(pBackProjector);
	ASTRA_DELETE(pFirstForwardProjector);
}
//----------------------------------------------------------------------------------------

} // namespace astra
