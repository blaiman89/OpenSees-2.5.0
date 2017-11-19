/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 6223 $
// $Date: 2016-03-01 16:32:50 -0800 (Tue, 01 Mar 2016) $
// $URL: svn://peera.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/analysis/integrator/AlphaOS.h $

#ifndef AlphaOS_h
#define AlphaOS_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/05
// Revision: A
//
// Description: This file contains the class definition for AlphaOS.
// AlphaOS is an algorithmic class for performing a transient analysis
// using the Alpha-Operator-Splitting integration scheme.
// The parameter alpha corresponds to 1+alpha_{HHT}.

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class AlphaOS : public TransientIntegrator
{
public:
    // constructors
    AlphaOS();
    AlphaOS(double alpha,
        bool updElemDisp = false);
    AlphaOS(double alpha, double beta, double gamma,
        bool updElemDisp = false);
    
    // destructor
    ~AlphaOS();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    
    int domainChanged(void);
    int newStep(double deltaT);
    int revertToLastStep(void);
    int update(const Vector &deltaU);
    int commit(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:
    virtual int formElementResidual(void);
    
private:
    double alpha;
    double beta;
    double gamma;
    bool updElemDisp;  // a flag indicating if element displacements are updated during commit
    double deltaT;
    
    int updateCount;                // method should only have one update per step
    double c1, c2, c3;              // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot;  // response quantities at time t
    Vector *U, *Udot, *Udotdot;     // response quantities at time t+deltaT
    Vector *Ualpha, *Ualphadot;     // response quantities at time t+alpha*deltaT
    Vector *Upt;                    // predictor displacements at time t
};

#endif
