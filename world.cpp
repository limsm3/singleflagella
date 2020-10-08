#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean
    inputName = m_inputData.GetStringOpt("input-file");		// string
	
	// Physical parameters
//	RodLength = m_inputData.GetScalarOpt("RodLength");      // meter -> replaced with helix parameter and axial length
    helixpitch = m_inputData.GetScalarOpt("helixpitch");    // meter
    helixradius = m_inputData.GetScalarOpt("helixradius");  // meter
    gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
    maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
	maxIterContact = m_inputData.GetIntOpt("maxIterContact");				// maximum number of iterations
	rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter
//	numVertices = m_inputData.GetIntOpt("numVertices");     // int_num -> replaced with delta_l
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s
	epsilon = m_inputData.GetScalarOpt("epsilon");
	
	headSize = m_inputData.GetScalarOpt("headSize");
	C_translation = m_inputData.GetScalarOpt("C-translation");
	C_rotation = m_inputData.GetScalarOpt("C-rotation");

	useRSS = m_inputData.GetBoolOpt("use-RSS"); // If true, RSS will be used. Otherwise, RFT.
	includeContact = m_inputData.GetBoolOpt("include-contact");
	
	eta_per = 4.0 * M_PI * viscosity / ( log( 2.0 * helixpitch / rodRadius) + 0.5 );
	eta_par = 2.0 * M_PI * viscosity / ( log( 2.0 * helixpitch / rodRadius) - 0.5 );
	
	axisLengthInput = m_inputData.GetScalarOpt("axisLengthInput");
	deltaLengthInput = m_inputData.GetScalarOpt("deltaLengthInput");

	//////////////// needs to be fixed - world global variable new variable flagellaLength, headLength

	nTurn = axisLengthInput / helixpitch;
	flagellaLength = nTurn * sqrt(((2 * M_PI * helixradius) * (2 * M_PI * helixradius)) + (helixpitch * helixpitch));

	// geometry of helix
	distf = 0.028; // right now it is fixed but should be an input
	double newRodLength = (flagellaLength * 1) + helixradius + (deltaLengthInput * 2); // RODLENGTH CORRECTED
	RodLength = newRodLength;
	//int newNe = RodLength / deltaLengthInput;
	/*
	double nTurn = axisLengthInput / helixpitch;
    double newRodLength = nTurn * sqrt( (2 * M_PI * helixradius) * (2 * M_PI * helixradius) + helixpitch * helixpitch );

    RodLength = 2*newRodLength;

    int newNe = RodLength / deltaLengthInput;*/

    //numVertices = newNe + 1; 

	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus
	
    // Read input file to get angular velocity
    ReadOmegaData();

    smallDist2 = pow( rodRadius/1000.0, 2);

    minDistace = 1000.0; 

    contactNum = 0; // Initialize number of contacts to 0
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::ReadOmegaData()
{
	ifstream infile;

    infile.open (inputName.c_str() );
    if (!infile.is_open())
    {
		cout << "Unable to open file to read omega";
		timeStep = Nstep; // we are exiting
	}

	numOmegaPoints = 0;
    double a, b;
    while (infile >> a >> b)
    {
		numOmegaPoints++;
        timeSeries.push_back(a);
        omegaSeries.push_back(b);
    }
    infile.close();
    
    currentOmegaIndex = 0; // keeps track of current angular velocity
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	ReadOmegaData();

	ostringstream name;
	name.precision(4);
	name << fixed;
    name << "datafiles/simDER";
	name << "_numvertex_" << rod->nv;
	name << "_numhead_"<< rod->headNode;
	name << "_crotation" << C_rotation;
    name << "_axisLength_" << axisLengthInput;
    name << "_helixPitch_" << helixpitch;
    name << "_helixRadius_" << helixradius;
	name << "_omega_" << omegaSeries[currentOmegaIndex];
	name << "_totalTime_" << totalTime;
	name << ".txt";

    outfile.open(name.str().c_str());
    outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) return;
	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	//  data output every 0.1 seconds.
	
	if (fmod(timeStep,1000)==0) 
	{
		for (int i = 0; i < rod->nv; i++)
		{
			Vector3d xCurrent = 100 * rod->getVertex(i);

			outfile << currentTime << ", 0, " << xCurrent(0) << ", " << xCurrent(1) << ", " << xCurrent(2) << ", " << endl;
		}
		/*
		Vector3d xCurrent = rod->getVertex(rod->headNode);
		Vector3d xCurrent2 = rod->getVertex(rod->headNode-1);
		Vector3d headForce = f_0_rod_0;
		Vector3d leftflagella = f_1_rod_0;
		Vector3d rightflagella = f_2_rod_0;
		//double Forcehead = computeReactionForce();

		outfile << currentTime << " " << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << " " << xCurrent2(0) << " " << xCurrent2(1) << " " << xCurrent2(2) << " " << headForce(0) << " " << headForce(1) << " " << headForce(2) < < " " << xCurrent2(2) << " " << headForce(0) << " " << headForce(1) << " " << headForce(2) << " " << leftflagella(0) << " " << leftflagella(1) << " " << leftflagella(2) < <  " " << rightflagella(0) << " " << rightflagella(1) << " " << rightflagella(2) << endl;
	*/
	}
		/*
	for (int i = 0; i < rod->nv; i++)
	{
		Vector3d xCurrent = rod->getVertex(i);

		outfile << currentTime << " " << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << endl;
	}
*/
	}

void world::setRodStepper()
{
	rodGeometry();

	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
		youngM, shearM, RodLength, distf);

	rodBoundaryCondition();

	rod->setup();

	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce( *rod, *stepper);
	m_bendingForce = new elasticBendingForce( *rod, *stepper);
	m_twistingForce = new elasticTwistingForce( *rod, *stepper);
	m_inertialForce = new inertialForce( *rod, *stepper);
	m_gravityForce = new externalGravityForce( *rod, *stepper, gVector);
	
	if (includeContact==true) // If contact should be included, declare that force
		m_externalContactForce =new externalContactForce( *rod, *stepper);

	// dampingForce added
	if (useRSS == true)
		m_RegularizedStokeslet = new RegularizedStokeslet(*rod, *stepper, viscosity, epsilon);
	else
		m_dampingForce = new dampingForce(*rod, *stepper, viscosity, eta_per, eta_par);
	
	m_externalHeadForce = new externalHeadForce(*rod, *stepper, viscosity, headSize, C_translation, C_rotation); // force on head
	
	rod->updateTimeStep();
	
	timeStep = 0;
	currentTime = 0.0;

	Nstep = totalTime/deltaTime;

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = M_PI * pow(rodRadius ,4)/4.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;
	int totVertices = 2 * (numVertices) + 7;
	ne = totVertices - 1;
	nv = totVertices;

	reactionForce = VectorXd::Zero(rod->ndof);
}

// Setup geometry
void world::rodGeometry()
{
	double helixA = helixradius;
	double helixB = helixpitch / (2.0 * M_PI);

	double nTurn = axisLengthInput / helixpitch;
	double flagellaLength = nTurn * sqrt((2 * M_PI * helixA) * (2 * M_PI * helixA) + helixpitch * helixpitch);

	int newNe = flagellaLength / deltaLengthInput;

	numVertices = newNe + 1;
	int totVertices = 1* (numVertices) + 3; // NOTE THIS CHANGE
	vertices = MatrixXd(totVertices, 3);

	double helixT = flagellaLength / sqrt(helixA * helixA + helixB * helixB);
	double delta_t = helixT / (numVertices - 1); // step for t->[0, T]

	// geometry of helix
	double delta_l = flagellaLength / (numVertices - 1);
	double RodLength = flagellaLength*1 + helixradius + delta_l*2; //RODLENGTH CORRECTED

	vertices(0, 0) = 0.0;
	vertices(0, 1) = 0.0;
	vertices(0, 2) = 0.0;

	vertices(1, 0) = - delta_l;
	vertices(1, 1) = 0.0;
	vertices(1, 2) = 0.0;

	vertices(2, 0) = - 2.0 * delta_l;
	vertices(2, 1) = 0.0;
	vertices(2, 2) = 0.0;
	
	int j = 3;
	for (double tt = 0.0; j < numVertices; tt += delta_t)
	{
		vertices(j, 0) = helixB * tt;
		vertices(j, 1) = helixA * cos(tt);
		vertices(j, 2) = -helixA * sin(tt);
		j++;
	}

}

void world::rodBoundaryCondition()
{
/*
	rod->setVertexBoundaryCondition(rod->getVertex(numVertices+2), numVertices+2);
	rod->setVertexBoundaryCondition(rod->getVertex(numVertices+1), numVertices+1);
	rod->setVertexBoundaryCondition(rod->getVertex(numVertices+3), numVertices+3);
*/
}
	

void world::updateTimeStep()
{
	if (currentOmegaIndex < numOmegaPoints-1 && timeSeries[currentOmegaIndex + 1] <= currentTime)
	{
		currentOmegaIndex++;
	}
	
	deltaTwist = omegaSeries[currentOmegaIndex] * (2.0*M_PI/60.0) * deltaTime;

    rod->setTwistBoundaryCondition(rod->undeformedTwist(1) + deltaTwist, 1); // NOTE THIS CHANGE
    /*rod->setTwistBoundaryCondition(rod->undeformedTwist(numVertices+5) + deltaTwist, (numVertices+5));
    /* rod->setTwistBoundaryCondition(rod->undeformedTwist(numVertices-1) - deltaTwist, (numVertices-1));
	rod->setTwistBoundaryCondition(rod->undeformedTwist(numVertices+3) - deltaTwist, (numVertices+3)); */
	
	// Start with a trial solution for our solution x
	rod->updateGuess(); // x = x0 + u * dt

	if (includeContact==true) m_externalContactForce->setZeroForce();
		

	// compute hydrodynamic force - viscous force appllied. make change on here. Write your own. Just add viscous force.
	if (useRSS == true) m_RegularizedStokeslet->prepareForViscousForce();

	// solve the EOM
	updateEachRod();

	// compute contact
	if (includeContact == true)	
	{
		prepareForContact();
		int contiter = 0;
		// resolve contact
		while (contactNum != 0 && contiter <= maxIterContact)
		{
			if (render == 1)
			{
				cout << "Contact detected, contact number = " << contactNum << endl;
			}

			rod->updateGuess();
			updateEachRod();
			
			contiter = contiter+1;
		}
	}

	computeReactionForce(); // Use this function if we need to calculate reaction force

	// update time step
	rod->updateTimeStep();

	currentTime += deltaTime;
		
	timeStep++;
}

void world::updateEachRod()
{
	double normf = forceTol * 10.0;	
	double normf0 = 0;
	
	bool solved = false;
	
	iter = 0;
		
	while (solved == false)
	{
		rod->prepareForIteration();
		
		stepper->setZero();

		// Compute the forces and the jacobians
		m_inertialForce->computeFi();
		m_inertialForce->computeJi();
			
		m_stretchForce->computeFs();
		m_stretchForce->computeJs();
			
		m_bendingForce->computeFb();
		m_bendingForce->computeJb();
		
		m_twistingForce->computeFt();
		m_twistingForce->computeJt();

		m_gravityForce->computeFg();
		m_gravityForce->computeJg();

		if (includeContact==true) m_externalContactForce->computeFc();

		if (useRSS==true)
			m_RegularizedStokeslet->computeFrs();
		else
		{
			m_dampingForce->computeFd(); // DampingForce added
			m_dampingForce->computeJd(); // DampingForce added
		}
		
		m_externalHeadForce->computeFh();
		m_externalHeadForce->computeJh();
		
		// Compute norm of the force equations.
		normf = 0.0;
		for (int i=0; i < rod->uncons; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}

		normf = sqrt(normf);

		if (iter == 0) normf0 = normf;
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}
		
		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			rod->updateNewtonX(totalForce);

			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}

	if (render) 
	{
		cout << "Time: " << currentTime << " iter=" << iter << endl;
	}
	
	if (solved == false)
	{
		timeStep = Nstep; // we are exiting
	}
}

int world::simulationRunning()
{
	if (timeStep<Nstep) 
		return 1;
	else 
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int j)
{
	return rod->x[j];
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}

void world::prepareForContact()
{
	contactNum = 0;

	for (int j = 0; j < ne; j++)
	{
		const Vector3d x_1 = rod->getVertex(j);
		const Vector3d x_2 = rod->getVertex(j+1);

		for (int l = 0; l < ne; l++)
		{
			// Two edges side by side are always going to contact. We will ignore it.
			if ( abs(l-j) <=1 ) continue;
			
			const Vector3d x_3 = rod->getVertex(l);
			const Vector3d x_4 = rod->getVertex(l+1);

			// compute min length of two segements 
			Vector3d c1;
   			Vector3d c2;
   			double s, t; 

   			double sqrdist = ClosestPtSegmentSegment(x_1,x_2,x_3,x_4, s, t, c1, c2);

    		//	if (minDistace > sqrt(sqrdist))
    		//	{
    		//		minDistace = sqrt(sqrdist);
    		//	}
		

			if (sqrdist < (2.0*rodRadius) * (2.0*rodRadius) )
			{
   				double pen = rodRadius + rodRadius - sqrt(sqrdist);

				Vector3d n = c2-c1; // contact normal
				n.normalize();

				double wi = s;
				double wj = t;
				double d =  2.0 * rodRadius;
				double mdij = d - pen;
	
				double del_r_i   = 0.5 * (mdij - d) * wi;
				double del_r_ip1 = 0.5 * (mdij - d) * (1.0 - wi);	

				double del_r_j   = 0.5 * (d - mdij) * wj;
				double del_r_jp1 = 0.5 * (d - mdij) * (1.0 - wj); 	

				double mi   = rod->massArray(4*j);
				double mip1 = rod->massArray(4*(j+1));

				double mj   = rod->massArray(4*l);
				double mjp1 = rod->massArray(4*(l+1));

				Vector3d f1 = n * del_r_i   * mi   / (deltaTime*deltaTime);
				Vector3d f2 = n * del_r_ip1 * mip1 / (deltaTime*deltaTime);

				Vector3d f3 = n * del_r_j   * mj   / (deltaTime*deltaTime);
				Vector3d f4 = n * del_r_jp1 * mjp1 / (deltaTime*deltaTime);

				m_externalContactForce->getContactForce(j + 0, f1);
				m_externalContactForce->getContactForce(j + 1, f2);
				m_externalContactForce->getContactForce(l + 0, f3);
				m_externalContactForce->getContactForce(l + 1, f4);

				contactNum = contactNum + 1;
   			}
		}
	}
}

double world::ClosestPtSegmentSegment( const Vector3d& p1, const Vector3d& q1, const Vector3d& p2, const Vector3d& q2, double& s, double& t, Vector3d& c1, Vector3d& c2 )
{ 
  Vector3d d1 = q1 - p1; // Direction vector of segment S1
  Vector3d d2 = q2 - p2; // Direction vector of segment S2
  Vector3d r = p1 - p2;
  double a = d1.dot(d1); // Squared length of segment S1, always nonnegative
  double e = d2.dot(d2); // Squared length of segment S2, always nonnegative
  double f = d2.dot(r);
  
  // Check if either or both segments degenerate into points
  if (a <= smallDist2 && e <= smallDist2) 
  {
    // Both segments degenerate into points
    s = t = 0.0;
    c1 = p1;
    c2 = p2;
    return (c1-c2).dot(c1-c2);
  }

  if (a <= smallDist2) 
  {
    // First segment degenerates into a point
    s = 0.0;
    t = f / e; // s = 0 => t = (b*s + f) / e = f / e
    t = clamp(t, 0.0, 1.0);
  }
  else 
  {
    double c = d1.dot(r);
    if (e <= smallDist2) 
    {
      // Second segment degenerates into a point
      t = 0.0;
      s = clamp(-c/a, 0.0, 1.0); // t = 0 => s = (b*t - c) / a = -c / a
    } 
    else 
    {
      // The general nondegenerate case starts here
      double b = d1.dot(d2);
      double denom = a*e-b*b; // Always nonnegative
      
      // If segments not parallel, compute closest point on L1 to L2, and
      // clamp to segment S1. Else pick arbitrary s (here 0)
      if (denom != 0.0) 
      {
        s = clamp((b*f - c*e) / denom, 0.0, 1.0);
      } 
      else s = 0.0;
      
      // Compute point on L2 closest to S1(s) using
      // t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e
      t = (b*s + f) / e;
      
      // If t in [0,1] done. Else clamp t, recompute s for the new value
      // of t using s = Dot((P2+D2*t)-P1,D1) / Dot(D1,D1)= (t*b - c) / a
      // and clamp s to [0, 1]
      if (t < 0.0) 
      {
        t = 0.0;
        s = clamp(-c / a, 0.0, 1.0);
      } 
      else if (t > 1.0) 
      {
        t = 1.0;
        s = clamp((b - c) / a, 0.0, 1.0);
      }
    }
  }
  
  c1 = p1 + d1 * s;
  c2 = p2 + d2 * t;
  return (c1 - c2).dot(c1 - c2);
}

void world::computeReactionForce()
{
	// THIS FUNCTION IS OUTDATED AND SHOULD NOT BE USED
	
	if (useRSS == true) // RSS
		reactionForce = m_inertialForce->ForceVec - m_bendingForce->ForceVec - m_twistingForce->ForceVec - m_stretchForce->ForceVec - m_RegularizedStokeslet->ForceVec;
	else // RFT
		reactionForce = m_inertialForce->ForceVec - m_bendingForce->ForceVec - m_twistingForce->ForceVec - m_stretchForce->ForceVec - m_dampingForce->ForceVec;
	
	int k = rod->headNode;

	f_1_rod_0(0) = reactionForce(4 * (k - 1));
	f_1_rod_0(1) = reactionForce(4 * (k - 1) + 1);
	f_1_rod_0(2) = reactionForce(4 * (k - 1) + 2);

	f_0_rod_0(0) = reactionForce(4*k); 
	f_0_rod_0(1) = reactionForce(4*k+1); 
	f_0_rod_0(2) = reactionForce(4*k+2); 

	f_2_rod_0(0) = reactionForce(4*(k+1));
	f_2_rod_0(1) = reactionForce(4 * (k + 1)+1);
	f_2_rod_0(2) = reactionForce(4 * (k + 1) + 2);

	//f_rod_0 = f_0_rod_0 + f_1_rod_0;
}
