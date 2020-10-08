#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "elasticTwistingForce.h"
#include "externalGravityForce.h"
#include "inertialForce.h"

// include external force
#include "RegularizedStokeslet.h" // RSS-based viscous force
#include "dampingForce.h" // RFT-based viscous force
#include "externalHeadForce.h" // Force on head
#include "externalContactForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"

//
// clamp function is copied from BASim code (Columbia Univ.)
/** Clamps scalar to the range [min,max]. */
template <typename T> inline T clamp( const T& scalar, const T& min, const T& max) 
{
  if (scalar < min) return min;
  if (scalar > max) return max;
  return scalar;
}

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	void setRodStepper();
	void updateTimeStep();
	int simulationRunning();
	int numPoints();
	double getScaledCoordinate(int j);
	double getCurrentTime();
	double getTotalTime();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);
		
private:

	// Physical parameters
	double RodLength;
	double helixradius, helixpitch;
	double rodRadius;
	int numVertices;
	double youngM;
	double Poisson;
	double shearM;
	double deltaTime;
	double totalTime;
	double density;
	Vector3d gVector;
	double viscosity;
	double epsilon;
	double distf;
	double nTurn;
	double flagellaLength;
	
	// Damping force parameters
	double eta_per;
	double eta_par;
	double headSize;
	double C_translation, C_rotation; // numerical prefactor to multiply the force against head rotation and translation
	
	double tol, stol;
	int maxIter; // maximum number of iterations
	int maxIterContact;
	double characteristicForce;
	double forceTol;
	
	// Geometry
	MatrixXd vertices;
	double currentTime;
	
	// Rod
	elasticRod *rod;

	// std::vector<elasticRod*> rodsVector;
	
	// set up the time stepper
	timeStepper *stepper;
	double *totalForce;

//	std::vector<timeStepper*> stepperVector;
//	std::vector<double*> totalForceVector;
	
	// declare the forces
	elasticStretchingForce *m_stretchForce;
	elasticBendingForce *m_bendingForce;
	elasticTwistingForce *m_twistingForce;

	bool useRSS; // if true, RSS will be used. Otherwise, RFT will be used.
	dampingForce *m_dampingForce; // Resistive force theory based model
	RegularizedStokeslet *m_RegularizedStokeslet; // Regularized Stokeslet Segment (RSS) based model
	externalHeadForce *m_externalHeadForce; // Force on head
	
	bool includeContact; // if true, contact will be handled. Otherwise, it will be ignored.
	externalContactForce *m_externalContactForce; // contact force

	inertialForce *m_inertialForce;	// inertial force
	externalGravityForce *m_gravityForce; // gravity
	
//	std::vector<elasticStretchingForce*> v_stretchForce;
//	std::vector<elasticBendingForce*> v_bendingForce;
//	std::vector<elasticTwistingForce*> v_twistingForce;
//	std::vector<inertialForce*> v_inertialForce;
//	std::vector<externalGravityForce*> v_gravityForce;
//	std::vector<externalContactForce*> v_externalContactForce;

	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry();
	void rodBoundaryCondition();

	// Variables about angular velocity
	double deltaTwist; // actuation
	// double deltaTheta;
	string inputName; // name of file containing omega
	std::vector<double> timeSeries, omegaSeries;
    int numOmegaPoints, currentOmegaIndex;
    void ReadOmegaData();
    
	bool render; // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?

	void updateEachRod();

	double axisLength;

	double ClosestPtSegmentSegment( const Vector3d& p1, const Vector3d& q1, const Vector3d& p2, const Vector3d& q2, double& s, double& t, Vector3d& c1, Vector3d& c2 );
	double smallDist2;

	void prepareForContact();
	int ne, nv;

	int contactNum;

	double minDistace;

	double axisLengthInput;
	double deltaLengthInput;

	VectorXd reactionForce;

	void computeReactionForce();

	Vector3d f_0_rod_0;
	Vector3d f_1_rod_0;
	Vector3d f_2_rod_0;
	Vector3d f_rod_0;

	Vector3d f_0_rod_1;
	Vector3d f_1_rod_1;
	Vector3d f_2_rod_1;
	Vector3d f_rod_1;
};

#endif
