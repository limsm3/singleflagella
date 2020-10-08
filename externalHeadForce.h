#ifndef EXTERNALHEADFORCE_H
#define EXTERNALHEADFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class externalHeadForce
{
public:
	externalHeadForce(elasticRod &m_rod, timeStepper &m_stepper, double m_viscosity, 
		double m_headSize, double m_C_translation, double m_C_rotation);
	~externalHeadForce();
	void computeFh();
	void computeJh();

	VectorXd ForceVec; // added

private:
	elasticRod *rod;
	timeStepper *stepper;
	double viscosity;
    Vector3d t, u, u_i, u_f, du, f;
    int ind, indx, indy;
    Matrix3d Id3, jac;
    double dt;
    double C_translation, C_rotation; // numerical prefactors for translation and rotation
    int headNode;
    double headSize;
	double delta_head;
	double omega_head;
};

#endif
