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

// $Revision: 1.15 $
// $Date: 2009/08/19 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/LinearCrdTransf3d1414.cpp,v $


// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Modified: 04/2005 Andreas Schellenberg (getBasicTrialVel, getBasicTrialAccel)
//
// Purpose: This file contains the implementation for the 
// LinearCrdTransf3d14 class. LinearCrdTransf3d14 is a linear
// transformation for a planar frame between the global 
// and basic coordinate systems
// Sensitivity: QGu UCSD 2009

#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>
#include <elementAPI.h>
#include <string>
#include <LinearCrdTransf3d14.h>

// initialize static variables
Matrix LinearCrdTransf3d14::Tlg(12,12);
Matrix LinearCrdTransf3d14::kg(12,12);

void* OPS_LinearCrdTransf3d14()
{
	if(OPS_GetNumRemainingInputArgs() < 4) {
		opserr<<"insufficient arguments for LinearCrdTransf3d\n";
		return 0;
	}

	// get tag
	int tag;
	int numData = 1;
	if(OPS_GetIntInput(&numData,&tag) < 0) return 0;

	// get vector
	Vector vec(3);
	double* vptr = &vec(0);
	numData = 3;
	if(OPS_GetDoubleInput(&numData,vptr) < 0) return 0;

	// get option
	Vector jntOffsetI(3), jntOffsetJ(3);
	double *iptr=&jntOffsetI(0), *jptr=&jntOffsetJ(0);
	while(OPS_GetNumRemainingInputArgs() > 6) {
		std::string type = OPS_GetString();
		if(type == "-jntOffset") {
			if(OPS_GetDoubleInput(&numData,iptr) < 0) return 0;
			if(OPS_GetDoubleInput(&numData,jptr) < 0) return 0;
		}
	}

	return new LinearCrdTransf3d14(tag,vec,jntOffsetI,jntOffsetJ);
}

// constructor:
LinearCrdTransf3d14::LinearCrdTransf3d14(int tag, const Vector &vecInLocXZPlane):
CrdTransf(tag, CRDTR_TAG_LinearCrdTransf3d14),
nodeIPtr(0), nodeJPtr(0), 
L(0), nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] = 0.0;
        
        R[2][0] = vecInLocXZPlane(0);
        R[2][1] = vecInLocXZPlane(1);
        R[2][2] = vecInLocXZPlane(2);
        
        // Does nothing
}


// constructor:
LinearCrdTransf3d14::LinearCrdTransf3d14(int tag, const Vector &vecInLocXZPlane,
										 const Vector &rigJntOffset1,
										 const Vector &rigJntOffset2):
CrdTransf(tag, CRDTR_TAG_LinearCrdTransf3d14),
	nodeIPtr(0), nodeJPtr(0),
	nodeIOffset(0), nodeJOffset(0), L(0),
	nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 3; j++)
			R[i][j] = 0.0;

	R[2][0] = vecInLocXZPlane(0);
	R[2][1] = vecInLocXZPlane(1);
	R[2][2] = vecInLocXZPlane(2);

	// check rigid joint offset for node I
	if (&rigJntOffset1 == 0 || rigJntOffset1.Size() != 3 ) 
	{
		opserr << "LinearCrdTransf3d::LinearCrdTransf3d:  Invalid rigid joint offset vector for node I\n";
		opserr << "Size must be 3\n";      
	}
	else if (rigJntOffset1.Norm() > 0.0) 
	{
		nodeIOffset = new double[3];
		nodeIOffset[0] = rigJntOffset1(0);
		nodeIOffset[1] = rigJntOffset1(1);
		nodeIOffset[2] = rigJntOffset1(2);
	}

	// check rigid joint offset for node J
	if (&rigJntOffset2 == 0 || rigJntOffset2.Size() != 3 )
	{
		opserr << "LinearCrdTransf3d::LinearCrdTransf3d:  Invalid rigid joint offset vector for node J\n";
		opserr << "Size must be 3\n";      
	}
	else if (rigJntOffset2.Norm() > 0.0) 
	{
		nodeJOffset = new double[3];
		nodeJOffset[0] = rigJntOffset2(0);
		nodeJOffset[1] = rigJntOffset2(1);
		nodeJOffset[2] = rigJntOffset2(2);
	}
}


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
LinearCrdTransf3d14::LinearCrdTransf3d14():
CrdTransf(0, CRDTR_TAG_LinearCrdTransf3d14),
nodeIPtr(0), nodeJPtr(0),
L(0), nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] = 0.0;
}


// destructor:
LinearCrdTransf3d14::~LinearCrdTransf3d14() 
{
    if (nodeIInitialDisp != 0)
        delete [] nodeIInitialDisp;
    if (nodeJInitialDisp != 0)
        delete [] nodeJInitialDisp;
}


int
LinearCrdTransf3d14::commitState(void)
{
    return 0;
}


int
LinearCrdTransf3d14::revertToLastCommit(void)
{
    return 0;
}


int
LinearCrdTransf3d14::revertToStart(void)
{
    return 0;
}


int 
LinearCrdTransf3d14::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
    int error;
    
    nodeIPtr = nodeIPointer;
    nodeJPtr = nodeJPointer;
    
    if ((!nodeIPtr) || (!nodeJPtr))
    {
        opserr << "\nLinearCrdTransf3d14::initialize";
        opserr << "\ninvalid pointers to the element nodes\n";
        return -1;
    }
    
    // see if there is some initial displacements at nodes
    if (initialDispChecked == false) {
        const Vector &nodeIDisp = nodeIPtr->getDisp();
        const Vector &nodeJDisp = nodeJPtr->getDisp();
        for (int i=0; i<7; i++)
            if (nodeIDisp(i) != 0.0) {
                nodeIInitialDisp = new double [7];
                for (int j=0; j<7; j++)
                    nodeIInitialDisp[j] = nodeIDisp(j);
                i = 7;
            }
            
            for (int j=0; j<7; j++)
                if (nodeJDisp(j) != 0.0) {
                    nodeJInitialDisp = new double [7];
                    for (int i=0; i<7; i++)
                        nodeJInitialDisp[i] = nodeJDisp(i);
                    j = 7;
                }
                
                initialDispChecked = true;
    }
    
    // get element length and orientation
    if ((error = this->computeElemtLengthAndOrient()))
        return error;
    
    static Vector XAxis(3);
    static Vector YAxis(3);
    static Vector ZAxis(3);
    
    // get 3by3 rotation matrix
    if ((error = this->getLocalAxes(XAxis, YAxis, ZAxis)))
        return error;
    
    return 0;
}


int
LinearCrdTransf3d14::update(void)
{
    return 0;
}


int 
LinearCrdTransf3d14::computeElemtLengthAndOrient()
{
    // element projection
    static Vector dx(3);
    
    const Vector &ndICoords = nodeIPtr->getCrds();
    const Vector &ndJCoords = nodeJPtr->getCrds();
    
    dx(0) = ndJCoords(0) - ndICoords(0);
    dx(1) = ndJCoords(1) - ndICoords(1);
    dx(2) = ndJCoords(2) - ndICoords(2);
    
    if (nodeIInitialDisp != 0) {
        dx(0) -= nodeIInitialDisp[0];
        dx(1) -= nodeIInitialDisp[1];
        dx(2) -= nodeIInitialDisp[2];
    }
    
    if (nodeJInitialDisp != 0) {
        dx(0) += nodeJInitialDisp[0];
        dx(1) += nodeJInitialDisp[1];
        dx(2) += nodeJInitialDisp[2];
    }
    
    // calculate the element length
    L = dx.Norm();
    
    if (L == 0.0) {
        opserr << "\nLinearCrdTransf3d14::computeElemtLengthAndOrien: 0 length\n";
        return -2;  
    }
    
    // calculate the element local x axis components (direction cossines)
    // wrt to the global coordinates
    R[0][0] = dx(0)/L;
    R[0][1] = dx(1)/L;
    R[0][2] = dx(2)/L;
    
    return 0;
}


void
LinearCrdTransf3d14::compTransfMatrixLocalGlobal(Matrix &Tlg)
{
	// setup transformation matrix from local to global
	Tlg.Zero();

	Tlg(0, 0) = Tlg(3, 3) = Tlg(6, 6) = Tlg(9, 9) = R[0][0];
	Tlg(0, 1) = Tlg(3, 4) = Tlg(6, 7) = Tlg(9, 10) = R[0][1];
	Tlg(0, 2) = Tlg(3, 5) = Tlg(6, 8) = Tlg(9, 11) = R[0][2];
	Tlg(1, 0) = Tlg(4, 3) = Tlg(7, 6) = Tlg(10, 9) = R[1][0];
	Tlg(1, 1) = Tlg(4, 4) = Tlg(7, 7) = Tlg(10, 10) = R[1][1];
	Tlg(1, 2) = Tlg(4, 5) = Tlg(7, 8) = Tlg(10, 11) = R[1][2];
	Tlg(2, 0) = Tlg(5, 3) = Tlg(8, 6) = Tlg(11, 9) = R[2][0];
	Tlg(2, 1) = Tlg(5, 4) = Tlg(8, 7) = Tlg(11, 10) = R[2][1];
	Tlg(2, 2) = Tlg(5, 5) = Tlg(8, 8) = Tlg(11, 11) = R[2][2];
}


int
LinearCrdTransf3d14::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
    // Compute y = v cross x
    // Note: v(i) is stored in R[2][i]
    static Vector vAxis(3);
    vAxis(0) = R[2][0];	vAxis(1) = R[2][1];	vAxis(2) = R[2][2];
    
    static Vector xAxis(3);
    xAxis(0) = R[0][0];	xAxis(1) = R[0][1];	xAxis(2) = R[0][2];
    XAxis(0) = xAxis(0);    XAxis(1) = xAxis(1);    XAxis(2) = xAxis(2);
    
    static Vector yAxis(3);
    yAxis(0) = vAxis(1)*xAxis(2) - vAxis(2)*xAxis(1);
    yAxis(1) = vAxis(2)*xAxis(0) - vAxis(0)*xAxis(2);
    yAxis(2) = vAxis(0)*xAxis(1) - vAxis(1)*xAxis(0);
    
    double ynorm = yAxis.Norm();
    
    if (ynorm == 0) {
        opserr << "\nLinearCrdTransf3d14::getLocalAxes";
        opserr << "\nvector v that defines plane xz is parallel to x axis\n";
        return -3;
    }
    
    yAxis /= ynorm;
    
    YAxis(0) = yAxis(0);    YAxis(1) = yAxis(1);    YAxis(2) = yAxis(2);
    
    // Compute z = x cross y
    static Vector zAxis(3);
    
    zAxis(0) = xAxis(1)*yAxis(2) - xAxis(2)*yAxis(1);
    zAxis(1) = xAxis(2)*yAxis(0) - xAxis(0)*yAxis(2);
    zAxis(2) = xAxis(0)*yAxis(1) - xAxis(1)*yAxis(0);
    ZAxis(0) = zAxis(0);    ZAxis(1) = zAxis(1);    ZAxis(2) = zAxis(2);
    
    // Fill in transformation matrix
    R[1][0] = yAxis(0);
    R[1][1] = yAxis(1);
    R[1][2] = yAxis(2);
    
    R[2][0] = zAxis(0);
    R[2][1] = zAxis(1);
    R[2][2] = zAxis(2);
    
    return 0;
}


double 
LinearCrdTransf3d14::getInitialLength(void)
{
    return L;
}


double 
LinearCrdTransf3d14::getDeformedLength(void)
{
    return L;
}


const Vector &
LinearCrdTransf3d14::getBasicTrialDisp (void)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getTrialDisp();
    const Vector &disp2 = nodeJPtr->getTrialDisp();
    
    static double ug[14];
    for (int i = 0; i < 7; i++) {
        ug[i]   = disp1(i);
        ug[i+7] = disp2(i);
    }
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<7; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<7; j++)
            ug[j+7] -= nodeJInitialDisp[j];
    }
    
    double oneOverL = 1.0/L;
    
    static Vector ub(9);
    
    static double ul[14];
// these need to be changed for 14 dofs, keep 7 and 14 the same
    ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
    ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];
    
    ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
    ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
    ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];
    
    ul[6]=ug[6]; // do not transform warping

    ul[7]  = R[0][0]*ug[7] + R[0][1]*ug[8] + R[0][2]*ug[9];
    ul[8]  = R[1][0]*ug[7] + R[1][1]*ug[8] + R[1][2]*ug[9];
    ul[9]  = R[2][0]*ug[7] + R[2][1]*ug[8] + R[2][2]*ug[9];
    
    ul[10] = R[0][0]*ug[10] + R[0][1]*ug[11] + R[0][2]*ug[12];
    ul[11] = R[1][0]*ug[10] + R[1][1]*ug[11] + R[1][2]*ug[12];
    ul[12] = R[2][0]*ug[10] + R[2][1]*ug[11] + R[2][2]*ug[12];
	
    ul[13] = ug[13]; // do not transform warping
    
    ub(8) = ul[7] - ul[0];
    double tmp;
    tmp = oneOverL*(ul[1]-ul[8]);
    ub(1) = ul[5] + tmp;
    ub(5) = ul[12] + tmp;
    tmp = oneOverL*(ul[9]-ul[2]);
    ub(2) = ul[4] + tmp;
    ub(6) = ul[11] + tmp;
    ub(0) = (-ul[10] + ul[3])/2;
	ub(4) = -ub(0);
	ub(3) = ul[6];
	ub(7) = ul[13];
    
    return ub;
}


const Vector &
LinearCrdTransf3d14::getBasicIncrDisp (void)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getIncrDisp();
    const Vector &disp2 = nodeJPtr->getIncrDisp();
    
    static double ug[14];
    for (int i = 0; i < 7; i++) {
        ug[i]   = disp1(i);
        ug[i+7] = disp2(i);
    }
    
    double oneOverL = 1.0/L;
    
    static Vector ub(9);
    
    static double ul[14];
// these need to be changed for 14 dofs, keep 7 and 14 the same
    ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
    ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];
    
    ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
    ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
    ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];
    
    ul[6]=ug[6]; // do not transform warping

    ul[7]  = R[0][0]*ug[7] + R[0][1]*ug[8] + R[0][2]*ug[9];
    ul[8]  = R[1][0]*ug[7] + R[1][1]*ug[8] + R[1][2]*ug[9];
    ul[9]  = R[2][0]*ug[7] + R[2][1]*ug[8] + R[2][2]*ug[9];
    
    ul[10] = R[0][0]*ug[10] + R[0][1]*ug[11] + R[0][2]*ug[12];
    ul[11] = R[1][0]*ug[10] + R[1][1]*ug[11] + R[1][2]*ug[12];
    ul[12] = R[2][0]*ug[10] + R[2][1]*ug[11] + R[2][2]*ug[12];
	
    ul[13] = ug[13]; // do not transform warping
    
    ub(8) = ul[7] - ul[0];
    double tmp;
    tmp = oneOverL*(ul[1]-ul[8]);
    ub(1) = ul[5] + tmp;
    ub(5) = ul[12] + tmp;
    tmp = oneOverL*(ul[9]-ul[2]);
    ub(2) = ul[4] + tmp;
    ub(6) = ul[11] + tmp;
    ub(0) = (-ul[10] + ul[3])/2;
    ub(4) = -ub(0);
    ub(3) = ul[6];
    ub(7) = ul[13];
    
    return ub;
}


const Vector &
LinearCrdTransf3d14::getBasicIncrDeltaDisp(void)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getIncrDeltaDisp();
    const Vector &disp2 = nodeJPtr->getIncrDeltaDisp();
    
    static double ug[14];
    for (int i = 0; i < 7; i++) {
        ug[i]   = disp1(i);
        ug[i+7] = disp2(i);
    }
    
    double oneOverL = 1.0/L;
    
    static Vector ub(9);
    
    static double ul[14];
    
// these need to be changed for 14 dofs, keep 7 and 14 the same
    ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
    ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];
    
    ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
    ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
    ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];
    
    ul[6]=ug[6]; // do not transform warping

    ul[7]  = R[0][0]*ug[7] + R[0][1]*ug[8] + R[0][2]*ug[9];
    ul[8]  = R[1][0]*ug[7] + R[1][1]*ug[8] + R[1][2]*ug[9];
    ul[9]  = R[2][0]*ug[7] + R[2][1]*ug[8] + R[2][2]*ug[9];
    
    ul[10] = R[0][0]*ug[10] + R[0][1]*ug[11] + R[0][2]*ug[12];
    ul[11] = R[1][0]*ug[10] + R[1][1]*ug[11] + R[1][2]*ug[12];
    ul[12] = R[2][0]*ug[10] + R[2][1]*ug[11] + R[2][2]*ug[12];
	
    ul[13] = ug[13]; // do not transform warping
    
    ub(8) = ul[7] - ul[0];
    double tmp;
    tmp = oneOverL*(ul[1]-ul[8]);
    ub(1) = ul[5] + tmp;
    ub(5) = ul[12] + tmp;
    tmp = oneOverL*(ul[9]-ul[2]);
    ub(2) = ul[4] + tmp;
    ub(6) = ul[11] + tmp;
    ub(0) = (-ul[10] + ul[3])/2;
    ub(4) = -ub(0);
    ub(3) = ul[6];
    ub(7) = ul[13];
    return ub;
}


const Vector &
LinearCrdTransf3d14::getBasicTrialVel(void)
{
	// determine global velocities
	const Vector &vel1 = nodeIPtr->getTrialVel();
	const Vector &vel2 = nodeJPtr->getTrialVel();
	
	static double vg[14];
	for (int i = 0; i < 7; i++) {
		vg[i]   = vel1(i);
		vg[i+7] = vel2(i);
	}
	
	double oneOverL = 1.0/L;
	
	static Vector vb(9);
	
	static double vl[14];
	
// these need to be changed for 14 dofs, keep 7 and 14 the same
    vl[0]  = R[0][0]*vg[0] + R[0][1]*vg[1] + R[0][2]*vg[2];
    vl[1]  = R[1][0]*vg[0] + R[1][1]*vg[1] + R[1][2]*vg[2];
    vl[2]  = R[2][0]*vg[0] + R[2][1]*vg[1] + R[2][2]*vg[2];
    
    vl[3]  = R[0][0]*vg[3] + R[0][1]*vg[4] + R[0][2]*vg[5];
    vl[4]  = R[1][0]*vg[3] + R[1][1]*vg[4] + R[1][2]*vg[5];
    vl[5]  = R[2][0]*vg[3] + R[2][1]*vg[4] + R[2][2]*vg[5];
    
    vl[6]=vg[6]; // do not transform warping

    vl[7]  = R[0][0]*vg[7] + R[0][1]*vg[8] + R[0][2]*vg[9];
    vl[8]  = R[1][0]*vg[7] + R[1][1]*vg[8] + R[1][2]*vg[9];
    vl[9]  = R[2][0]*vg[7] + R[2][1]*vg[8] + R[2][2]*vg[9];
    
    vl[10] = R[0][0]*vg[10] + R[0][1]*vg[11] + R[0][2]*vg[12];
    vl[11] = R[1][0]*vg[10] + R[1][1]*vg[11] + R[1][2]*vg[12];
    vl[12] = R[2][0]*vg[10] + R[2][1]*vg[11] + R[2][2]*vg[12];
	
    vl[13] = vg[13]; // do not transform warping
	
    vb(8) = vl[7] - vl[0];
    double tmp;
    tmp = oneOverL*(vl[1]-vl[8]);
    vb(1) = vl[5] + tmp;
    vb(5) = vl[12] + tmp;
    tmp = oneOverL*(vl[9]-vl[2]);
    vb(2) = vl[4] + tmp;
    vb(6) = vl[11] + tmp;
    vb(0) = (-vl[10] + vl[3])/2;
    vb(4) = -vb(0);
    vb(3) = vl[6];
    vb(7) = vl[13];
	
    return vb;
}


const Vector &
LinearCrdTransf3d14::getBasicTrialAccel(void)
{
	// determine global accelerations
	const Vector &accel1 = nodeIPtr->getTrialAccel();
	const Vector &accel2 = nodeJPtr->getTrialAccel();
	
	static double ag[14];
	for (int i = 0; i < 7; i++) {
		ag[i]   = accel1(i);
		ag[i+7] = accel2(i);
	}
	
	double oneOverL = 1.0/L;
	
	static Vector ab(9);
	
	static double al[14];
	
// these need to be changed for 14 dofs, keep 7 and 14 the same
    al[0]  = R[0][0]*ag[0] + R[0][1]*ag[1] + R[0][2]*ag[2];
    al[1]  = R[1][0]*ag[0] + R[1][1]*ag[1] + R[1][2]*ag[2];
    al[2]  = R[2][0]*ag[0] + R[2][1]*ag[1] + R[2][2]*ag[2];
    
    al[3]  = R[0][0]*ag[3] + R[0][1]*ag[4] + R[0][2]*ag[5];
    al[4]  = R[1][0]*ag[3] + R[1][1]*ag[4] + R[1][2]*ag[5];
    al[5]  = R[2][0]*ag[3] + R[2][1]*ag[4] + R[2][2]*ag[5];
    
    al[6]=ag[6]; // do not transform warping

    al[7]  = R[0][0]*ag[7] + R[0][1]*ag[8] + R[0][2]*ag[9];
    al[8]  = R[1][0]*ag[7] + R[1][1]*ag[8] + R[1][2]*ag[9];
    al[9]  = R[2][0]*ag[7] + R[2][1]*ag[8] + R[2][2]*ag[9];
    
    al[10] = R[0][0]*ag[10] + R[0][1]*ag[11] + R[0][2]*ag[12];
    al[11] = R[1][0]*ag[10] + R[1][1]*ag[11] + R[1][2]*ag[12];
    al[12] = R[2][0]*ag[10] + R[2][1]*ag[11] + R[2][2]*ag[12];
	
    al[13] = ag[13]; // do not transform warping
	
	
    ab(8) = al[7] - al[0];
    double tmp;
    tmp = oneOverL*(al[1]-al[8]);
    ab(1) = al[5] + tmp;
    ab(5) = al[12] + tmp;
    tmp = oneOverL*(al[9]-al[2]);
    ab(2) = al[4] + tmp;
    ab(6) = al[11] + tmp;
    ab(0) = (-al[10] + al[3])/2;
    ab(4) = -ab(0);
    ab(3) = al[6];
    ab(7) = al[13];
	
    return ab;
}


const Vector &
LinearCrdTransf3d14::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
    // transform resisting forces from the basic system to local coordinates
    static double pl[14];
    
    double q0 = pb(0); //T1
    double q1 = pb(1); //Mz1
    double q2 = pb(2); //My1
    double q3 = pb(3); //Bx1
    double q4 = pb(4); //T2
    double q5 = pb(5); //Mz2
    double q6 = pb(6); //My2
    double q7 = pb(7); //Bx2
    double q8 = pb(8); //P
    
    double oneOverL = 1.0/L;
    
    pl[0]  = -q8; //P1
    pl[1]  =  oneOverL*(q1+q5); //Vy1
    pl[2]  = -oneOverL*(q2+q6); //Vz1
    pl[3]  =  0.5*(q0-q4); //T1
    pl[4]  =  q2; //My1
    pl[5]  =  q1; //Mz1
    pl[6]  =  q3; //Bx1
    pl[7]  =  q8; //P2
    pl[8]  = -pl[1]; //Vy2
    pl[9]  = -pl[2]; //Vz2
    pl[10] = -pl[3]; //T2
    pl[11] =  q6; //My2
    pl[12] =  q5;  //Mz2
    pl[13] =  q7; //Bx2
    
    // transform resisting forces  from local to global coordinates
    static Vector pg(14);
    
    pg(0) = R[0][0]*pl[0] + R[1][0]*pl[1] + R[2][0]*pl[2];
    pg(1) = R[0][1]*pl[0] + R[1][1]*pl[1] + R[2][1]*pl[2];
    pg(2) = R[0][2]*pl[0] + R[1][2]*pl[1] + R[2][2]*pl[2];
    
    pg(3) = R[0][0]*pl[3] + R[1][0]*pl[4] + R[2][0]*pl[5];
    pg(4) = R[0][1]*pl[3] + R[1][1]*pl[4] + R[2][1]*pl[5];
    pg(5) = R[0][2]*pl[3] + R[1][2]*pl[4] + R[2][2]*pl[5];
    
    pg(6) = pl[6];

    pg(7) = R[0][0]*pl[7] + R[1][0]*pl[8] + R[2][0]*pl[9];
    pg(8) = R[0][1]*pl[7] + R[1][1]*pl[8] + R[2][1]*pl[9];
    pg(9) = R[0][2]*pl[7] + R[1][2]*pl[8] + R[2][2]*pl[9];
    
    pg(10) = R[0][0]*pl[10] + R[1][0]*pl[11] + R[2][0]*pl[12];
    pg(11) = R[0][1]*pl[10] + R[1][1]*pl[11] + R[2][1]*pl[12];
    pg(12) = R[0][2]*pl[10] + R[1][2]*pl[11] + R[2][2]*pl[12];

    pg(13) = pl[13];
	
    return pg;
}

const Matrix &
LinearCrdTransf3d14::getGlobalStiffMatrix (const Matrix &KB, const Vector &pb)
{
    static Matrix kg(14,14);	// Global stiffness for return
    static double kb[9][9];	// Basic stiffness
    static double kl[14][14];	// Local stiffness
    static double tmp[14][14];	// Temporary storage
    
    double oneOverL = 1.0/L;
    
    int i,j;
    for (i = 0; i < 9; i++)
        for (j = 0; j < 9; j++)
            kb[i][j] = KB(i,j);
        
        // Transform basic stiffness to local system
        // First compute kb*T_{bl}
        for (i = 0; i < 9; i++) {
            tmp[i][0]  = -kb[i][8];
            tmp[i][1]  =  oneOverL*(kb[i][1]+kb[i][5]);
            tmp[i][2]  = -oneOverL*(kb[i][2]+kb[i][6]);
            tmp[i][3]  =  0.5*(kb[i][0]-kb[i][4]);
            tmp[i][4]  =  kb[i][2];
            tmp[i][5]  =  kb[i][1];
            tmp[i][6]  =  kb[i][3];
            tmp[i][7]  =  kb[i][8];
            tmp[i][8]  = -tmp[i][1];
            tmp[i][9]  = -tmp[i][2];
            tmp[i][10] = -tmp[i][3];
            tmp[i][11] =  kb[i][6];
	    tmp[i][12] =  kb[i][5];
	    tmp[i][13] =  kb[i][7];
        }
        
        // Now compute T'_{bl}*(kb*T_{bl})
        for (i = 0; i < 14; i++) {
            kl[0][i]  = -tmp[8][i];
            kl[1][i]  =  oneOverL*(tmp[1][i]+tmp[5][i]);
            kl[2][i]  = -oneOverL*(tmp[2][i]+tmp[6][i]);
            kl[3][i]  =  0.5*(tmp[0][i]-tmp[4][i]);
            kl[4][i]  =  tmp[2][i];
            kl[5][i]  =  tmp[1][i];
            kl[6][i]  =  tmp[3][i];
            kl[7][i]  =  tmp[8][i];
            kl[8][i]  = -kl[1][i];
            kl[9][i]  = -kl[2][i];
            kl[10][i] = -kl[3][i];
            kl[11][i] =  tmp[6][i];
	    kl[12][i] =  tmp[5][i];
	    kl[13][i] =  tmp[7][i];
        }
        
        // Transform local stiffness to global system
        // First compute kl*T_{lg}
        int m;
        for (m = 0; m < 14; m++) {
            tmp[m][0] = kl[m][0]*R[0][0] + kl[m][1]*R[1][0]  + kl[m][2]*R[2][0];
            tmp[m][1] = kl[m][0]*R[0][1] + kl[m][1]*R[1][1]  + kl[m][2]*R[2][1];
            tmp[m][2] = kl[m][0]*R[0][2] + kl[m][1]*R[1][2]  + kl[m][2]*R[2][2];
            
            tmp[m][3] = kl[m][3]*R[0][0] + kl[m][4]*R[1][0]  + kl[m][5]*R[2][0];
            tmp[m][4] = kl[m][3]*R[0][1] + kl[m][4]*R[1][1]  + kl[m][5]*R[2][1];
            tmp[m][5] = kl[m][3]*R[0][2] + kl[m][4]*R[1][2]  + kl[m][5]*R[2][2];
            
	    tmp[m][6] = kl[m][6];

            tmp[m][7] = kl[m][7]*R[0][0] + kl[m][8]*R[1][0]  + kl[m][9]*R[2][0];
            tmp[m][8] = kl[m][7]*R[0][1] + kl[m][8]*R[1][1]  + kl[m][9]*R[2][1];
            tmp[m][9] = kl[m][7]*R[0][2] + kl[m][8]*R[1][2]  + kl[m][9]*R[2][2];
            
            tmp[m][10] = kl[m][10]*R[0][0] + kl[m][11]*R[1][0] + kl[m][12]*R[2][0];
            tmp[m][11] = kl[m][10]*R[0][1] + kl[m][11]*R[1][1] + kl[m][12]*R[2][1];
            tmp[m][12] = kl[m][10]*R[0][2] + kl[m][11]*R[1][2] + kl[m][12]*R[2][2];
            
            tmp[m][13] = kl[m][13];
        }
        
        // Now compute T'_{lg}*(kl*T_{lg})
        for (m = 0; m < 14; m++) {
            kg(0,m) = R[0][0]*tmp[0][m] + R[1][0]*tmp[1][m]  + R[2][0]*tmp[2][m];
            kg(1,m) = R[0][1]*tmp[0][m] + R[1][1]*tmp[1][m]  + R[2][1]*tmp[2][m];
            kg(2,m) = R[0][2]*tmp[0][m] + R[1][2]*tmp[1][m]  + R[2][2]*tmp[2][m];
            
            kg(3,m) = R[0][0]*tmp[3][m] + R[1][0]*tmp[4][m]  + R[2][0]*tmp[5][m];
            kg(4,m) = R[0][1]*tmp[3][m] + R[1][1]*tmp[4][m]  + R[2][1]*tmp[5][m];
            kg(5,m) = R[0][2]*tmp[3][m] + R[1][2]*tmp[4][m]  + R[2][2]*tmp[5][m];
            
	    kg(6,m) = tmp[6][m];

            kg(7,m) = R[0][0]*tmp[7][m] + R[1][0]*tmp[8][m]  + R[2][0]*tmp[9][m];
            kg(8,m) = R[0][1]*tmp[7][m] + R[1][1]*tmp[8][m]  + R[2][1]*tmp[9][m];
            kg(9,m) = R[0][2]*tmp[7][m] + R[1][2]*tmp[8][m]  + R[2][2]*tmp[9][m];
            
            kg(10,m)  = R[0][0]*tmp[10][m] + R[1][0]*tmp[11][m] + R[2][0]*tmp[12][m];
            kg(11,m) = R[0][1]*tmp[10][m] + R[1][1]*tmp[11][m] + R[2][1]*tmp[12][m];
            kg(12,m) = R[0][2]*tmp[10][m] + R[1][2]*tmp[11][m] + R[2][2]*tmp[12][m];
            
	    kg(13,m) = tmp[13][m];
        }
        
        return kg;
}


const Matrix &
LinearCrdTransf3d14::getInitialGlobalStiffMatrix (const Matrix &KB)
{
    static Matrix kg(14,14);	// Global stiffness for return
    static double kb[9][9];		// Basic stiffness
    static double kl[14][14];	// Local stiffness
    static double tmp[14][14];	// Temporary storage
    
    double oneOverL = 1.0/L;
    
    int i,j;
    for (i = 0; i < 9; i++)
        for (j = 0; j < 9; j++)
            kb[i][j] = KB(i,j);
        
        // Transform basic stiffness to local system
        // First compute kb*T_{bl}
        for (i = 0; i < 9; i++) {
            tmp[i][0]  = -kb[i][8];
            tmp[i][1]  =  oneOverL*(kb[i][1]+kb[i][5]);
            tmp[i][2]  = -oneOverL*(kb[i][2]+kb[i][6]);
            tmp[i][3]  =  0.5*(kb[i][0]-kb[i][4]);
            tmp[i][4]  =  kb[i][2];
            tmp[i][5]  =  kb[i][1];
            tmp[i][6]  =  kb[i][3];
            tmp[i][7]  =  kb[i][8];
            tmp[i][8]  = -tmp[i][1];
            tmp[i][9]  = -tmp[i][2];
            tmp[i][10] = -tmp[i][3];
            tmp[i][11] =  kb[i][6];
	    tmp[i][12] =  kb[i][5];
	    tmp[i][13] =  kb[i][7];
        }
        
        // Now compute T'_{bl}*(kb*T_{bl})
        for (i = 0; i < 14; i++) {
            kl[0][i]  = -tmp[8][i];
            kl[1][i]  =  oneOverL*(tmp[1][i]+tmp[5][i]);
            kl[2][i]  = -oneOverL*(tmp[2][i]+tmp[6][i]);
            kl[3][i]  =  0.5*(tmp[0][i]-tmp[4][i]);
            kl[4][i]  =  tmp[2][i];
            kl[5][i]  =  tmp[1][i];
            kl[6][i]  =  tmp[3][i];
            kl[7][i]  =  tmp[8][i];
            kl[8][i]  = -kl[1][i];
            kl[9][i]  = -kl[2][i];
            kl[10][i] = -kl[3][i];
            kl[11][i] =  tmp[6][i];
	    kl[12][i] =  tmp[5][i];
	    kl[13][i] =  tmp[7][i];
        }
        
        // Transform local stiffness to global system
        // First compute kl*T_{lg}
        int m;
        for (m = 0; m < 14; m++) {
            tmp[m][0] = kl[m][0]*R[0][0] + kl[m][1]*R[1][0]  + kl[m][2]*R[2][0];
            tmp[m][1] = kl[m][0]*R[0][1] + kl[m][1]*R[1][1]  + kl[m][2]*R[2][1];
            tmp[m][2] = kl[m][0]*R[0][2] + kl[m][1]*R[1][2]  + kl[m][2]*R[2][2];
            
            tmp[m][3] = kl[m][3]*R[0][0] + kl[m][4]*R[1][0]  + kl[m][5]*R[2][0];
            tmp[m][4] = kl[m][3]*R[0][1] + kl[m][4]*R[1][1]  + kl[m][5]*R[2][1];
            tmp[m][5] = kl[m][3]*R[0][2] + kl[m][4]*R[1][2]  + kl[m][5]*R[2][2];
            
	    tmp[m][6] = kl[m][6];

            tmp[m][7] = kl[m][7]*R[0][0] + kl[m][8]*R[1][0]  + kl[m][9]*R[2][0];
            tmp[m][8] = kl[m][7]*R[0][1] + kl[m][8]*R[1][1]  + kl[m][9]*R[2][1];
            tmp[m][9] = kl[m][7]*R[0][2] + kl[m][8]*R[1][2]  + kl[m][9]*R[2][2];
            
            tmp[m][10] = kl[m][10]*R[0][0] + kl[m][11]*R[1][0] + kl[m][12]*R[2][0];
            tmp[m][11] = kl[m][10]*R[0][1] + kl[m][11]*R[1][1] + kl[m][12]*R[2][1];
            tmp[m][12] = kl[m][10]*R[0][2] + kl[m][11]*R[1][2] + kl[m][12]*R[2][2];
            
            tmp[m][13] = kl[m][13];
            
        }
        
        // Now compute T'_{lg}*(kl*T_{lg})
        for (m = 0; m < 14; m++) {
            kg(0,m) = R[0][0]*tmp[0][m] + R[1][0]*tmp[1][m]  + R[2][0]*tmp[2][m];
            kg(1,m) = R[0][1]*tmp[0][m] + R[1][1]*tmp[1][m]  + R[2][1]*tmp[2][m];
            kg(2,m) = R[0][2]*tmp[0][m] + R[1][2]*tmp[1][m]  + R[2][2]*tmp[2][m];
            
            kg(3,m) = R[0][0]*tmp[3][m] + R[1][0]*tmp[4][m]  + R[2][0]*tmp[5][m];
            kg(4,m) = R[0][1]*tmp[3][m] + R[1][1]*tmp[4][m]  + R[2][1]*tmp[5][m];
            kg(5,m) = R[0][2]*tmp[3][m] + R[1][2]*tmp[4][m]  + R[2][2]*tmp[5][m];
            
	    kg(6,m) = tmp[6][m];

            kg(7,m) = R[0][0]*tmp[7][m] + R[1][0]*tmp[8][m]  + R[2][0]*tmp[9][m];
            kg(8,m) = R[0][1]*tmp[7][m] + R[1][1]*tmp[8][m]  + R[2][1]*tmp[9][m];
            kg(9,m) = R[0][2]*tmp[7][m] + R[1][2]*tmp[8][m]  + R[2][2]*tmp[9][m];
            
            kg(10,m)  = R[0][0]*tmp[10][m] + R[1][0]*tmp[11][m] + R[2][0]*tmp[12][m];
            kg(11,m) = R[0][1]*tmp[10][m] + R[1][1]*tmp[11][m] + R[2][1]*tmp[12][m];
            kg(12,m) = R[0][2]*tmp[10][m] + R[1][2]*tmp[11][m] + R[2][2]*tmp[12][m];
            
	    kg(13,m) = tmp[13][m];
        }
        
        return kg;
}


CrdTransf *
LinearCrdTransf3d14::getCopy3d(void)
{
    // create a new instance of LinearCrdTransf3d14 
    
    LinearCrdTransf3d14 *theCopy;
    
    static Vector xz(3);
    xz(0) = R[2][0];
    xz(1) = R[2][1];
    xz(2) = R[2][2];
    
    theCopy = new LinearCrdTransf3d14(this->getTag(), xz);
    
    theCopy->nodeIPtr = nodeIPtr;
    theCopy->nodeJPtr = nodeJPtr;
    theCopy->L = L;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            theCopy->R[i][j] = R[i][j];
        
        return theCopy;
}


int 
LinearCrdTransf3d14::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    
    static Vector data(19);
    data(0) = this->getTag();
    data(1) = L;
    
    if (nodeIInitialDisp != 0) {
        data(2) = nodeIInitialDisp[0];
        data(3) = nodeIInitialDisp[1];
        data(4) = nodeIInitialDisp[2];
        data(5) = nodeIInitialDisp[3];
        data(6) = nodeIInitialDisp[4];
        data(7) = nodeIInitialDisp[5];
	data(8) = nodeIInitialDisp[6];
    } else {
        data(2)  = 0.0;
        data(3)  = 0.0;
        data(4) = 0.0;
        data(5) = 0.0;
        data(6) = 0.0;
        data(7) = 0.0;
	data(8) = 0.0;
    }
    
    if (nodeJInitialDisp != 0) {
        data(9)  = nodeJInitialDisp[0];
        data(10) = nodeJInitialDisp[1];
        data(11) = nodeJInitialDisp[2];
        data(12) = nodeJInitialDisp[3];
        data(13) = nodeJInitialDisp[4];
        data(14) = nodeJInitialDisp[5];
	data(15) = nodeJInitialDisp[6];
    } else {
        data(9)  = 0.0;
        data(10) = 0.0;
        data(11) = 0.0;
        data(12) = 0.0;
        data(13) = 0.0;
	data(14) = 0.0;
	data(15) = 0.0;
    }
    
    data(16) = R[2][0];
    data(17) = R[2][1];
    data(18) = R[2][2];
    
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
        opserr << "LinearCrdTransf3d14::sendSelf - failed to send Vector\n";
        
        return res;
    }
    
    return res;
}


int 
LinearCrdTransf3d14::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    
    static Vector data(19);
    
    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
        opserr << "LinearCrdTransf3d14::recvSelf - failed to receive Vector\n";
        
        return res;
    }
    
    this->setTag((int)data(0));
    L = data(1);
    data(0) = this->getTag();
    data(1) = L;
    
    
    int flag;
    int i,j;
    
            
    flag = 0;
    for (i=2; i<=8; i++)
	if (data(i) != 0.0)
	    flag = 1;

    if (flag == 1) {
	if (nodeIInitialDisp == 0)
	    nodeIInitialDisp = new double[7];
	for (i=2, j=0; i<=8; i++, j++)
	    nodeIInitialDisp[j] = data(i);
    }
	
    flag = 0;
    for (i=9; i<=15; i++)
	if (data(i) != 0.0)
	    flag = 1;

    if (flag == 1) {
	if (nodeJInitialDisp == 0)
	    nodeJInitialDisp = new double [7];
	for (i=9, j=0; i<=15; i++, j++)
	    nodeJInitialDisp[j] = data(i);
    }
	
    R[2][0] = data(16);
    R[2][1] = data(17);
    R[2][2] = data(18);
    
    
    initialDispChecked = true;
    return res;
}


const Matrix &
LinearCrdTransf3d14::getGlobalMatrixFromLocal(const Matrix &ml)
{
	this->compTransfMatrixLocalGlobal(Tlg);  // OPTIMIZE LATER
	kg.addMatrixTripleProduct(0.0, Tlg, ml, 1.0);  // OPTIMIZE LATER

	return kg;
}


const Vector &
LinearCrdTransf3d14::getPointGlobalCoordFromLocal(const Vector &xl)
{
    static Vector xg(3);
    
    //xg = nodeIPtr->getCrds() + nodeIOffset;
    xg = nodeIPtr->getCrds();
    
    if (nodeIInitialDisp != 0) {
        xg(0) -= nodeIInitialDisp[0];
        xg(1) -= nodeIInitialDisp[1];
        xg(2) -= nodeIInitialDisp[2];
    }
    
    
    // xg = xg + Rlj'*xl
    //xg.addMatrixTransposeVector(1.0, Rlj, xl, 1.0);
    xg(0) += R[0][0]*xl(0) + R[1][0]*xl(1) + R[2][0]*xl(2);
    xg(1) += R[0][1]*xl(0) + R[1][1]*xl(1) + R[2][1]*xl(2);
    xg(2) += R[0][2]*xl(0) + R[1][2]*xl(1) + R[2][2]*xl(2);
    
    return xg;  
}


const Vector &
LinearCrdTransf3d14::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getTrialDisp();
    const Vector &disp2 = nodeJPtr->getTrialDisp();
    
    static double ug[14];
    for (int i = 0; i < 7; i++)
    {
        ug[i]   = disp1(i);
        ug[i+7] = disp2(i);
    }
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<7; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<7; j++)
            ug[j+7] -= nodeJInitialDisp[j];
    }
    
    
    
    // transform global end displacements to local coordinates
    //ul.addMatrixVector(0.0, Tlg,  ug, 1.0);       //  ul = Tlg *  ug;
    static double ul[14];
    
    ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
    ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];
    
    ul[8]  = R[1][0]*ug[7] + R[1][1]*ug[8] + R[1][2]*ug[9];
    ul[9]  = R[2][0]*ug[7] + R[2][1]*ug[8] + R[2][2]*ug[9];
    
    // compute displacements at point xi, in local coordinates
    static double uxl[3];
    static Vector uxg(3);
    
    uxl[0] = uxb(0) +        ul[0];
    uxl[1] = uxb(1) + (1-xi)*ul[1] + xi*ul[8];
    uxl[2] = uxb(2) + (1-xi)*ul[2] + xi*ul[9];
    
    // rotate displacements to global coordinates
    // uxg = Rlj'*uxl
    //uxg.addMatrixTransposeVector(0.0, Rlj, uxl, 1.0);
    uxg(0) = R[0][0]*uxl[0] + R[1][0]*uxl[1] + R[2][0]*uxl[2];
    uxg(1) = R[0][1]*uxl[0] + R[1][1]*uxl[1] + R[2][1]*uxl[2];
    uxg(2) = R[0][2]*uxl[0] + R[1][2]*uxl[1] + R[2][2]*uxl[2];
    
    return uxg;  
}


void
LinearCrdTransf3d14::Print(OPS_Stream &s, int flag)
{
    s << "\nCrdTransf: " << this->getTag() << " Type: LinearCrdTransf3d14";
}
