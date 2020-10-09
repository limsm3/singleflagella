#include "externalHeadForce.h"
#include <iostream>

externalHeadForce::externalHeadForce(elasticRod &m_rod, timeStepper &m_stepper, 
	double m_viscosity, double m_headSize, double m_C_translation, double m_C_rotation)
{
	rod = &m_rod;
	stepper = &m_stepper;
	viscosity = m_viscosity;			   
	dt = rod->dt;
	headSize = m_headSize;
	C_translation = m_C_translation;
	C_rotation = m_C_rotation;

	Id3<<1,0,0,
		0,1,0,
        0,0,1;
	
	f.setZero(3);
	ForceVec = VectorXd::Zero(rod->ndof); // added
	f_head = 0;
	headNode = rod->headNode; // index of head node
	delta_head = headSize;
}

externalHeadForce::~externalHeadForce()
{
	;
}

void externalHeadForce::computeFh()
{
	ForceVec = VectorXd::Zero(rod->ndof);

	//
	// Force on head against translation (Stokes Law)
	//
	u = rod->getVelocity(headNode);
	f = - C_translation * 6.0 * M_PI * viscosity * headSize * u;
	for (int k=0; k < 3; k++)
	{
		ind = 4 * headNode + k;
		stepper->addForce(ind, - f[k]); // subtracting external force
		ForceVec(ind) = ForceVec(ind) + f[k]; // added
	}
 
	//
	// Force on head against rotation (can be improved)
	//
	omega_head = rod->getOmega(headNode);	
	f_head = - C_rotation * 8.0 * M_PI * viscosity * pow(headSize, 3.0) * omega_head ; // force on headNode+1

	ind = 4 * (headNode) + 3;
	stepper->addForce(ind, - f_head); // subtracting external force
	ForceVec(ind) = ForceVec(ind) + f_head; // added	


}

void externalHeadForce::computeJh()
{
	// Remember that dF/dx = 1/dt * dF/dv

	//
	// Force on head against translation (Stokes Law)
	//
	jac = - C_translation * 6.0 * M_PI * viscosity * headSize / dt * Id3;
	for (int kx = 0; kx < 3; kx++)
	{
		indx = 4 * headNode + kx;
		for (int ky = 0; ky < 3; ky++)
		{
			indy = 4 * headNode + ky;
			stepper->addJacobian(indx, indy, - jac(kx,ky)); // subtracting external force
		}
	}
	
	//
	// Force on head against rotation (can be improved)
	//
	double jac_head = - C_rotation * 8.0 * M_PI * viscosity * pow(headSize, 3.0)  / dt; // force on headNode+1
	

		indx = 4 * (headNode) + 3;
		int ky = 3;
		indy = 4 * (headNode) + ky;
	    stepper->addJacobian(indx, indy, - jac_head); // subtracting external force
		
	

}
