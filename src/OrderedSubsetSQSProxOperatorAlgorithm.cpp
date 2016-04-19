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

#include "astra/OrderedSubsetSQSProxOperatorAlgorithm.h"

#include <boost/lexical_cast.hpp>

#include "astra/AstraObjectManager.h"
#include "astra/DataProjectorPolicies.h"
#include "astra/Logging.h"


using namespace std;

namespace astra {

#include "astra/Projector2DImpl.inl"

// type of the algorithm, needed to register with CAlgorithmFactory
std::string COrderedSubsetSQSProxOperatorAlgorithm::type = "OS-SQS-PROX";


//---------------------------------------------------------------------------------------
// Clear - Constructors
void COrderedSubsetSQSProxOperatorAlgorithm::_clear()
{
	//CReconstructionAlgorithm2D::_clear();
	CSartProxOperatorAlgorithm::_clear();
	//m_fLambda = 1.0f;
	m_pTempVol = NULL;
}

//---------------------------------------------------------------------------------------
// Clear - Public
void COrderedSubsetSQSProxOperatorAlgorithm::clear()
{
	//CReconstructionAlgorithm2D::clear();
	CSartProxOperatorAlgorithm::clear();

	ASTRA_DELETE(m_pTempVol);
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
COrderedSubsetSQSProxOperatorAlgorithm::COrderedSubsetSQSProxOperatorAlgorithm() 
{
	_clear();
}

//----------------------------------------------------------------------------------------
// Constructor
COrderedSubsetSQSProxOperatorAlgorithm::COrderedSubsetSQSProxOperatorAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction) 
{
	_clear();
	CSartProxOperatorAlgorithm::CSartProxOperatorAlgorithm(_pProjector, _pSinogram, _pReconstruction);
	//initialize(_pProjector, _pSinogram, _pReconstruction);
}

//----------------------------------------------------------------------------------------
// Constructor
COrderedSubsetSQSProxOperatorAlgorithm::COrderedSubsetSQSProxOperatorAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction,
							   int* _piProjectionOrder, 
							   int _iProjectionCount)
{
	_clear();
	CSartProxOperatorAlgorithm::CSartProxOperatorAlgorithm(_pProjector, _pSinogram, _pReconstruction, 
		_piProjectionOrder, _iProjectionCount);
	//initialize(_pProjector, _pSinogram, _pReconstruction, _piProjectionOrder, _iProjectionCount);
}

//----------------------------------------------------------------------------------------
// Destructor
COrderedSubsetSQSProxOperatorAlgorithm::~COrderedSubsetSQSProxOperatorAlgorithm() 
{
	clear();
}

//---------------------------------------------------------------------------------------
// Initialize - Config
bool COrderedSubsetSQSProxOperatorAlgorithm::initialize(const Config& _cfg)
{
	assert(_cfg.self);
	ConfigStackCheck<CAlgorithm> CC("OS-SQS-PROX", this, _cfg);
	
	// if already initialized, clear first
	if (m_bIsInitialized) {
		clear();
	}

	// initialization of parent class
	if (!CSartProxOperatorAlgorithm::initialize(_cfg)) {
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

	//// Lambda
	//m_fLambda = _cfg.self.getOptionNumerical("Lambda", m_fLambda);
	//CC.markOptionParsed("Lambda");

	//// Alpha
	//m_fAlpha = _cfg.self.getOptionNumerical("Alpha", m_fAlpha);
	//CC.markOptionParsed("Alpha");

	//// Clear RaySum after each sweep. Defaults to true.
	//m_bClearRayLength = _cfg.self.getOptionBool("ClearRayLength", m_bClearRayLength);
	//CC.markOptionParsed("ClearRayLength");

	//// Input volume
	//XMLNode node = _cfg.self.getSingleNode("ProxInputDataId");
	//ASTRA_CONFIG_CHECK(node, "PSART", "No Proximal Input tag specified.");
	//int id = boost::lexical_cast<int>(node.getContent());
	//m_pProxInput = dynamic_cast<CFloat32VolumeData2D*>(CData2DManager::getSingleton().get(id));
	//CC.markNodeParsed("ProxInputDataId");

	// create data objects
	m_pTempVol = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	//m_pTotalRayLength = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	//m_pTotalPixelWeight = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	//m_pDiffSinogram = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}


//----------------------------------------------------------------------------------------
bool COrderedSubsetSQSProxOperatorAlgorithm::_check()
{
	// check base class
	ASTRA_CONFIG_CHECK(CSartProxOperatorAlgorithm::_check(), "OS-SQS-PROX", "Error in ReconstructionAlgorithm2D initialization");

	return true;
}

//---------------------------------------------------------------------------------------
// Information - All
map<string,boost::any> COrderedSubsetSQSProxOperatorAlgorithm::getInformation() 
{
	map<string, boost::any> res;
	res["ProjectionOrder"] = getInformation("ProjectionOrder");
	return mergeMap<string,boost::any>(CReconstructionAlgorithm2D::getInformation(), res);
};

//---------------------------------------------------------------------------------------
// Information - Specific
boost::any COrderedSubsetSQSProxOperatorAlgorithm::getInformation(std::string _sIdentifier) 
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
void COrderedSubsetSQSProxOperatorAlgorithm::run(int _iNrIterations)
{
	// check initialized
    ASTRA_ASSERT(m_bIsInitialized);

    m_bShouldAbort = false;

    // data projectors
    CDataProjectorInterface* pForwardProjector;
    CDataProjectorInterface* pBackProjector;
    CDataProjectorInterface* pFirstForwardProjector;
    CDataProjectorInterface* pFirstBackProjector;

    // Initialize m_pReconstruction to zero.
    if (m_bClearReconstruction) {
        m_pReconstruction->setData(0.f);
    }

    // backprojection data projector
    pBackProjector = dispatchDataProjector(
            m_pProjector, 
            SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
            ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			SIRTBPPolicy(m_pReconstruction, m_pDiffSinogram, 
				m_pTotalPixelWeight, m_pTotalRayLength, m_pPreconditioner, 
				m_fAlpha * m_fLambda * 2, m_pW),  // OS-SQS Prox backprojection
            m_bUseSinogramMask, m_bUseReconstructionMask, true // options on/off
        ); 

    // also computes total pixel weight and total ray length
    pForwardProjector = dispatchDataProjector(
            m_pProjector, 
            SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
            ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
            DiffFPPolicy(m_pReconstruction, m_pDiffSinogram, 
			  m_pSinogram, m_pW),								// forward projection with difference calculation
            m_bUseSinogramMask, m_bUseReconstructionMask, true											 // options on/off
        );

    // first time forward projection data projector,
    // computes total ray length (sum of rows) and total pixel weights (sum of columns)
    pFirstForwardProjector = dispatchDataProjector(
            m_pProjector, 
            SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
            ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
            TotalRayLengthPolicy(m_pTotalRayLength, false, m_pW),
            m_bUseSinogramMask, m_bUseReconstructionMask, true 											 // options on/off
        );

	// Backproject the TotalRayLength to compute the required column sums for SQS
	pFirstBackProjector = dispatchDataProjector(
            m_pProjector, 
            SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
            ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			DefaultBPPolicy(m_pTotalPixelWeight, m_pTotalRayLength, m_pW),
            m_bUseSinogramMask, m_bUseReconstructionMask, true 											 // options on/off
        );

	// We assume the input weights are already square-rooted.
	////WLS?
	//if (m_pW != NULL) {
	//	// Take square root to prepare
	//	m_pW->sqrt();
	//	// Scale sinogram
	//	*m_pSinogram *= *m_pW;
	//}

	// Perform the first forward projection to compute ray lengths.
	m_pTotalRayLength->setData(0.0f);    
    pFirstForwardProjector->project();

	// Backproject the row sums to compute the pixel weights (column sums) i.e. A^T * A
	m_pTotalPixelWeight->setData(0.0f);
	pFirstBackProjector->project();
	// Scale by 2 lambda and add 1
	*m_pTotalPixelWeight *= (2. * m_fLambda);
	*m_pTotalPixelWeight += 1.f;
	// Scale by the number of projections
	*m_pTotalPixelWeight *= 1.0 / m_iProjectionCount;

	// The row sums are just set to 1, and untouched again.
	m_pTotalRayLength->setData(1.0f);

	//// end of init.
    //ttime = CPlatformDepSystemCode::getMSCount() - timer;

	// Number of voxels to update per projection for the second term (u - x)
	int iVoxelsPerProjection = m_pReconstruction->getSize() / m_iProjectionCount + 1;

    // iteration loop, each iteration loops over all available projections
    for (int iIteration = 0; iIteration < _iNrIterations && !m_bShouldAbort; ++iIteration) {
        // start timer
        m_ulTimer = CPlatformDepSystemCode::getMSCount();

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

			// Update tempVol = ProxIn - Reconstruction
			m_pTempVol->copyData(m_pProxInput->getData());
			*m_pTempVol -= *m_pReconstruction;

			// forward projection and difference calculation
            pForwardProjector->projectSingleProjection(iProjection);
            // backprojection
            pBackProjector->projectSingleProjection(iProjection);


			//Update all the voxels, not sure if wrong, but gives better results ...
			// Add alpha * tempVol / TotalPixel
			*m_pTempVol /= *m_pTotalPixelWeight;
			// CHECK if division is correct!!
			//*m_pTempVol *= m_fAlpha / m_iProjectionCount;
			*m_pTempVol *= m_fAlpha;
			*m_pReconstruction += *m_pTempVol;

			//// Update only a subset of the Reconstruction corresponding to the current subset
			//for (int i = iVoxelsPerProjection * iP; 
			//	i < min(iVoxelsPerProjection * (iP + 1), m_pReconstruction->getSize()); 
			//	++i) {
			//	// Get tempVol[i]
			//	float32 val = m_pTempVol->getData()[i];
			//	// Update alpha * tempVol / TotalPixel
			//	val /= m_pTotalPixelWeight->getData()[i];
			//	// CHECK if division here is correct!!
			//	val *= m_fAlpha / m_iProjectionCount;  
			//	// Update voxel
			//	m_pReconstruction->getData()[i] += val;
			//}			

            // update iteration count
            m_iIterationCount++;

            // We need to check here, as checking inside the BP (as in ART)
            // is not correct.
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
	ASTRA_DELETE(pFirstBackProjector);

	//for (int i=0; i < m_pReconstruction->getSize(); ++i) {
	//	ASTRA_INFO("voxel=%d val=%f", i, m_pReconstruction->getData()[i]);
	//}
}
//----------------------------------------------------------------------------------------

} // namespace astra
