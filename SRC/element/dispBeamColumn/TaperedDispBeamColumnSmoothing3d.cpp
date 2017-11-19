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

// $Revision: 1.25 $
// $Date: 2008/11/04 21:32:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/TaperedDispBeamColumnSmoothing3d.cpp,v $

// Modified by Xi Zhang from University of Sydney, Australia (include warping degrees of freedom). Refer to 
// Formulation and Implementation of Three-dimensional Doubly Symmetric Beam-Column Analyses with Warping Effects in OpenSees
// Research Report R917, School of Civil Engineering, University of Sydney.
// Description: This file contains the class definition for TaperedDispBeamColumnSmoothing3d which include warping degrees of freedom.
// Modified by Xi Zhang from University of Sydney, Australia (include warping degrees of freedom). Refer to 
// Formulation and Implementation of Three-dimensional Doubly Symmetric Beam-Column Analyses with Warping Effects in OpenSees
// Research Report R917, School of Civil Engineering, University of Sydney.

//#define DEBUG

#include <TaperedDispBeamColumnSmoothing3d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>
#include <Parameter.h>
using std::string;
using namespace std;


Matrix TaperedDispBeamColumnSmoothing3d::K(14, 14);
Vector TaperedDispBeamColumnSmoothing3d::P(14);
double TaperedDispBeamColumnSmoothing3d::workArea[200];

TaperedDispBeamColumnSmoothing3d::TaperedDispBeamColumnSmoothing3d(int tag, int nd1,
                                                 int nd2, int numSec,
                                                 SectionForceDeformation
                                                 ** s,
                                                 BeamIntegration & bi,
                                                 CrdTransf &
                                                 coordTransf,
                                                 double r)
:Element(tag, ELE_TAG_TaperedDispBeamColumnSmoothing3d), numSections(numSec),
theSections(0), crdTransf(0), beamInt(0), connectedExternalNodes(2),
Q(14), q(9), rho(r), parameterID(0)
{
    // Allocate arrays of pointers to SectionForceDeformations
    theSections = new SectionForceDeformation *[numSections];

    if (theSections == 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::TaperedDispBeamColumnSmoothing3d - failed to allocate section model pointer\n";
        exit(-1);
    }

    for (int i = 0; i < numSections; i++) {

        // Get copies of the material model for each integration point
        theSections[i] = s[i]->getCopy();

        // Check allocation
        if (theSections[i] == 0) {
            opserr <<
                "TaperedDispBeamColumnSmoothing3d::DispBeamColumn3d -- failed to get a copy of section model\n";
            exit(-1);
        }
    }

    beamInt = bi.getCopy();

    if (beamInt == 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::DispBeamColumn3d - failed to copy beam integration\n";
        exit(-1);
    }

    crdTransf = coordTransf.getCopy3d();

    if (crdTransf == 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::DispBeamColumn3d - failed to copy coordinate transformation\n";
        exit(-1);
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;


    theNodes[0] = 0;
    theNodes[1] = 0;

    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;
    q0[3] = 0.0;
    q0[4] = 0.0;
    q0[5] = 0.0;
    q0[6] = 0.0;
    q0[7] = 0.0;
    q0[8] = 0.0;
    q0[9] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    p0[3] = 0.0;
    p0[4] = 0.0;
}

TaperedDispBeamColumnSmoothing3d::TaperedDispBeamColumnSmoothing3d()
:  
Element(0, ELE_TAG_TaperedDispBeamColumnSmoothing3d),
numSections(0), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2), Q(14), q(9), rho(0.0), parameterID(0)
{

    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;
    q0[3] = 0.0;
    q0[4] = 0.0;
    q0[5] = 0.0;
    q0[6] = 0.0;
    q0[7] = 0.0;
    q0[8] = 0.0;
    q0[9] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    p0[3] = 0.0;
    p0[4] = 0.0;

    theNodes[0] = 0;
    theNodes[1] = 0;
}

TaperedDispBeamColumnSmoothing3d::~TaperedDispBeamColumnSmoothing3d()
{
    for (int i = 0; i < numSections; i++) {
        if (theSections[i])
            delete theSections[i];
    }

    // Delete the array of pointers to SectionForceDeformation pointer arrays
    if (theSections)
        delete[]theSections;

    if (crdTransf)
        delete crdTransf;

    if (beamInt != 0)
        delete beamInt;
}

int TaperedDispBeamColumnSmoothing3d::getNumExternalNodes() const 
{

    return 2;
}

const ID & TaperedDispBeamColumnSmoothing3d::getExternalNodes()
{
    
    return connectedExternalNodes;
}

Node **TaperedDispBeamColumnSmoothing3d::getNodePtrs()
{

    return theNodes;
}

int TaperedDispBeamColumnSmoothing3d::getNumDOF()
{

    return 14;
}

void TaperedDispBeamColumnSmoothing3d::setDomain(Domain * theDomain)
{
    // Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    if (theNodes[0] == 0 || theNodes[1] == 0) {
        opserr <<
            "FATAL ERROR TaperedDispBeamColumnSmoothing3d (tag: %d), node not found in domain",
            this->getTag();

        return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    // release this restriction
    /*if (dofNd1 != 6 || dofNd2 != 6) {
       //opserr << "FATAL ERROR TaperedDispBeamColumnSmoothing3d (tag: %d), has differing number of DOFs at its nodes",
       //   this->getTag());

       return;
       } */

    if (crdTransf->initialize(theNodes[0], theNodes[1])) {
        // Add some error check
    }

    double L = crdTransf->getInitialLength();

    if (L == 0.0) {
        // Add some error check
    }

    this->DomainComponent::setDomain(theDomain);


    this->update();
}

int TaperedDispBeamColumnSmoothing3d::commitState()
{

    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::commitState () - failed in base class";
    }

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
        retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}

int TaperedDispBeamColumnSmoothing3d::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
        retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int TaperedDispBeamColumnSmoothing3d::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
        retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int TaperedDispBeamColumnSmoothing3d::update(void)
{
    int err = 0;

    #ifdef DEBUG
    /*
    static int ntimes = 1;
	fprintf(stdout, "Entering TaperedDispBeamColumnSmoothing3d for the %d-th times\n", ntimes);
	++ntimes;
    */
    #endif

    // Update the transformation
    crdTransf->update();

    // Get basic deformations
    // this is collecting qn (Chang) | n (Syd), which is the natural
    // displacements (no rigid body)
    const Vector & v = crdTransf->getBasicTrialDisp();

    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;
    double oneOverLsquare = 1.0 / L / L;

    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);

    // Loop over the integration points
 
    for (int i = 0; i < numSections; i++) {

        int order = theSections[i]->getOrder();
        const ID & code = theSections[i]->getType();

        // this is about to take the qn vector and create v
        // this is the same e as used in fibersection3d
        Vector e(workArea, 12);  // Without smoothing e(..., 8)

        double xi6 = 6.0 * xi[i];
        double xi12 = 12.0 * xi[i];
        double xi1 = xi[i];
        double x3 = xi[i] * xi[i] * xi[i];
        double x2 = xi[i] * xi[i];   

        int j;

        // here total strain (not incremental) is used
        // e is v in Sydney's code
        // NOTE: This probably means that FiberSection3d needs reordering


  // SHAPE FUNCTIONS
        // N means shape function, dN means 1st derivative, ddN 2nd 
        // N* - * is the order of the shape function 3-cubic 1-linear
        // N*x -  x is the dof of the shape function
        // xi is bringing in x/L already
        double N31   = 1 - 3 * x2 + 2 * x3;
        double N32   = xi[i] * L * (1 - xi[i]) * (1 - xi[i]);
        double N33   = 3 * x2 - 2 * x3;
        double N34   = -xi[i] * L * (xi[i] - x2);
        double dN31  = 6 * oneOverL * (x2 - xi[i]); //modified
        double dN32  = 1 + 3 * x2 - 4 * xi[i];
        double dN33  = 6 * oneOverL * (xi[i] - x2);
        double dN34  = 3 * x2 - 2 * xi[i];
        double ddN31 = oneOverLsquare * (xi12 - 6);
        double ddN32 = oneOverL * (xi6 - 4);
        double ddN33 = oneOverLsquare * (6 - xi12);
        double ddN34 = oneOverL * (xi6 - 2);
        double dN12  = oneOverL;

	// The v referred to is the nodal displacements in the correct
	// frame. It has this order [rx1 rz1 ry1 rx1’ rx2 rz2 ry2 rx2’ e]
	// This is the new e:

	e(0)  = dN12 * v(8); // u prime
	e(1)  = v(0); // rx1
	e(2) = v(1); // rz1
	e(3) = v(2);  // ry1
	e(4) = v(4); // rx2
	e(5) = v(5); // rz2
	e(6) = v(6); // ry2
	e(7) = ddN32 * v(1) + ddN34 * v(5); // v dprime
	e(8) = -ddN32 * v(2) - ddN34 * v(6); // w dprime
	e(9) = N31 * v(0) + N32 * v(3) + N33 * v(4) + N34 * v(7);     // Theta
	e(10) = dN31 * v(0) + dN32 * v(3) + dN33 * v(4) + dN34 * v(7); // Theta prime
	e(11) = ddN31*v(0) + ddN32*v(3) + ddN33*v(4) + ddN34*v(7);     // Theta dprime
	
        // Set the section deformations
	err += theSections[i]->setTrialSectionDeformation(e);
    }

    if (err != 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::update() - failed setTrialSectionDeformations()\n";
        return err;
    }

    return 0;
}

const Matrix & TaperedDispBeamColumnSmoothing3d::getTangentStiff()
{
    static Matrix kb(9, 9);
    static Matrix N1(6,12);
    //static Matrix N1(6, 8);
    static Matrix N2(12,9);
    //static Matrix N2(8, 9);
    static Matrix N3(12,12);
    //static Matrix N3(8, 8);
    static Matrix kbPart1(9, 9);
    static Matrix Gmax(12,12); 
    //static Matrix Gmax(8, 8);
    static Matrix kbPart2(9, 9);


    const Vector & v = crdTransf->getBasicTrialDisp();

    // Zero for integral
    kb.Zero();

    q.Zero();

    double L = crdTransf->getInitialLength();
    double oneOverL       = 1.0 / L;
    double oneOverLsquare = oneOverL / L;
    double oneOverLcube   = oneOverLsquare / L;

    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        int order = theSections[i]->getOrder();
        const ID & code = theSections[i]->getType();

	// N1 is Q in Chang and Ndeld1 in Sydney
	// N2 is B in Chang and Ndeld2 in Sydney
	// N3 is made up here for calc purposes and is N1_transpose * ks * N1
	// kbPart1 is the elastic stiffness matrix part of kb
	// kb is the element stiffness matrix, natural
	// kbPart2 is the geometric stiffness matrix
	
        N1.Zero();
        N2.Zero();
        N3.Zero();
        kbPart1.Zero();
        Gmax.Zero();
        kbPart2.Zero();

	// SHAPE FUNCTIONS
        double xi6 = 6.0 * xi[i];
        double xi12 = 12.0 * xi[i];
        double xi1 = xi[i];
        double x3 = xi[i] * xi[i] * xi[i];
        double x2 = xi[i] * xi[i];    

  // SHAPE FUNCTIONS
        // N means shape function, dN means 1st derivative, ddN 2nd 
        // N* - * is the order of the shape function 3-cubic 1-linear
        // N*x -  x is the dof of the shape function
        // xi is bringing in x/L already
        double N31   = 1 - 3 * x2 + 2 * x3;
        double N32   = xi[i] * L * (1 - xi[i]) * (1 - xi[i]);
        double N33   = 3 * x2 - 2 * x3;
        double N34   = -xi[i] * L * (xi[i] - x2);
        double dN31  = 6 * oneOverL * (x2 - xi[i]); //modified
        double dN32  = 1 + 3 * x2 - 4 * xi[i];
        double dN33  = 6 * oneOverL * (xi[i] - x2);
        double dN34  = 3 * x2 - 2 * xi[i];
        double ddN31 = oneOverLsquare * (xi12 - 6);
        double ddN32 = oneOverL * (xi6 - 4);
        double ddN33 = oneOverLsquare * (6 - xi12);
        double ddN34 = oneOverL * (xi6 - 2);
        double dN12  = oneOverL;

	N1(0,0)  = 1.0;
        N1(0,2)  = (4*v(1)-v(5))/30;
        N1(0,3)  = (4*v(2)-v(6))/30;
        N1(0,5)  = (4*v(5)-v(1))/30;
        N1(0,6)  = (4*v(6)-v(2))/30;
        N1(1,7)  = 1.0;
        N1(1,8)  = N31 * v(0) + N32 * v(3) + N33 * v(4) + N34 * v(7); // Theta
        N1(1,9)  = -ddN32 * v(2) - ddN34 * v(6); // w dprime
        N1(2,7)  = N1(1,8);
        N1(2,8)  = -1.0;
        N1(2,9)  = ddN32 * v(1) + ddN34 * v(5); // v dprime
        N1(3,10) = dN31 * v(0) + dN32 * v(3) + dN33 * v(4) + dN34 * v(7); // Theta prime
        N1(4,11) = 1.0;
        N1(5,10) = 1.0; 

       	N2(0,8)  = dN12;
        N2(1,0)  = 1;
        N2(2,1)  = 1;
        N2(3,2)  = 1;
        N2(4,4)  = 1;
        N2(5,5)  = 1;
        N2(6,6)  = 1;
        N2(7,1)  = ddN32;  // rz1->v dprime
        N2(7,5)  = ddN34;  // rz2->v dprime
        N2(8,2)  = -ddN32; // ry1->w dprime
        N2(8,6)  = -ddN34; // ry2->w dprime
        N2(9,0)  = N31;    // rx1->theta
        N2(9,3)  = N32;    // rx1 prime->theta
        N2(9,4)  = N33;    // rx2->theta
        N2(9,7)  = N34;    // rx2 prime->theta
        N2(10,0) = dN31;   // rx1->theta prime
        N2(10,3) = dN32;   // rx1 prime->theta prime
        N2(10,4) = dN33;   // rx2->theta prime
        N2(10,7) = dN34;   // rx2 prime->theta prime
        N2(11,0) = ddN31;  // rx1->theta dprime
        N2(11,3) = ddN32;  // rx1 prime->theta dprime
        N2(11,4) = ddN33;  // rx2->theta dprime
        N2(11,7) = ddN34;  // rx2 prime->theta dprime 
	
	#ifdef DEBUG
	/*
	    static int counter = 1;

	    if (counter == 1) {
		for (int a = 0; i < 8; ++i) {
		    for (int b = 0; j < 9; ++j) { 
			fprintf(stdout, "%.4f ", N2(a,b));
		    }
		    fprintf(stdout, "\n");
		}
		fprintf(stdout, "\n");
	    }
	    ++counter;
	    */
	#endif

/*
    fprintf(stdout, "About to obtain the stiffness tangent dimensions: \n")
    int rows = ks.noRows();
    int cols = ks.noCols();
    printf("The number of rows: %i\n", rows);
    printf("The number of columns: %i\n\n", cols);
*/

        // Get the section tangent stiffness and stress resultant
        const Matrix & ks = theSections[i]->getSectionTangent();

	#ifdef DEBUG
	/*
           int rows = ks.noRows();
           int cols = ks.noCols();
	    // Print Statements for N1
	    for (int i = 0; i < rows; ++i) {
	    	for (int j = 0; j < cols; ++j) { 
    		    //fprintf(stdout, "N2(%d,%d) = %.4f ", i, j, N2(i,j));
		    fprintf(stdout, " %.4f ", ks(i,j));
		}
		fprintf(stdout, "\n");
	    }
	    fprintf(stdout, "\n"); */
	#endif

        //calculate kb, refer to Alemdar

        N3.addMatrixTripleProduct(0.0, N1, ks, 1.0);
        kbPart1.addMatrixTripleProduct(0.0, N2, N3, 1.0);
        const Vector &s = theSections[i]->getStressResultant();

	#ifdef DEBUG
	/*
	    fprintf(stdout, "getTangentStiff\n");
	    for(int c = 0; c < 6; ++c) 
		fprintf(stdout, "The value of s[%d] = %.4f\n", c, s(c));

	    fprintf(stdout, "\n");
	*/
	#endif

	// s needs to be in the following order
	// [P Mz My W B Tsv]
	 Gmax(2,2)   = (4/30)*s(0);
	 Gmax(3,3)   = Gmax(2,2);
	 Gmax(5,5)   = Gmax(2,2);
	 Gmax(6,6)   = Gmax(2,2);
	 Gmax(2,5)   = -s(0)/30;
	 Gmax(5,2)   = Gmax(2,5);
	 Gmax(3,6)   = Gmax(2,5);
	 Gmax(6,3)   = Gmax(2,5);
	 Gmax(7,9)   = s(2);
	 Gmax(9,7)   = Gmax(7,9);
	 Gmax(8,9)   = s(1);
	 Gmax(9,8)   = Gmax(8,9);
	 Gmax(10,10) = s(3); 
	 kbPart2.addMatrixTripleProduct(0.0, N2, Gmax, 1.0); 

        // Perform numerical integration
        //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
        double wti = wt[i];

        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 9; k++) {
                kb(j, k) += (kbPart1(j, k) + kbPart2(j, k)) * L * wti;
            }
        }

        static Vector qProduct1(12);
        static Vector qProduct2(9);
        qProduct1.Zero();
        qProduct2.Zero();
        qProduct1.addMatrixTransposeVector(0.0, N1, s, 1.0);
        qProduct2.addMatrixTransposeVector(0.0, N2, qProduct1, 1.0);

        for (int j = 0; j < 9; j++) {
            q(j) += qProduct2(j) * L * wti;
        }
    }

    // Transform to global stiffness
    K = crdTransf->getGlobalStiffMatrix(kb, q);

    return K;
}

const Matrix &TaperedDispBeamColumnSmoothing3d::getInitialBasicStiff()
{
    static Matrix kb(9, 9);
    static Matrix N1(6,12);
    //static Matrix N1(6, 8);
    static Matrix N2(12,9);
    //static Matrix N2(8, 9);
    static Matrix N3(12,12);
    //static Matrix N3(8, 8);
    static Matrix kbPart1(9, 9);
    static Matrix Gmax(12,12); 
    //static Matrix Gmax(8, 8);
    static Matrix kbPart2(9, 9);

    // Zero for integral
    kb.Zero();
    const Vector & v = crdTransf->getBasicTrialDisp();

    double L = crdTransf->getInitialLength();
    double oneOverL       = 1.0 / L;
    double oneOverLsquare = oneOverL / L;
    double oneOverLcube   = oneOverLsquare / L;

    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        int order = theSections[i]->getOrder();
        const ID & code = theSections[i]->getType();

        N1.Zero();
        N2.Zero();
        N3.Zero();
        kbPart1.Zero();
        Gmax.Zero();
        kbPart2.Zero();

	// SHAPE FUNCTIONS
        double xi6 = 6.0 * xi[i];
        double xi12 = 12.0 * xi[i];
        double xi1 = xi[i];
        double x3 = xi[i] * xi[i] * xi[i];
        double x2 = xi[i] * xi[i];    

        // SHAPE FUNCTIONS
        // N means shape function, dN means 1st derivative, ddN 2nd 
        // N* - * is the order of the shape function 3-cubic 1-linear
        // N*x -  x is the dof of the shape function
        // xi is bringing in x/L already
        double N31   = 1 - 3 * x2 + 2 * x3;
        double N32   = xi[i] * L * (1 - xi[i]) * (1 - xi[i]);
        double N33   = 3 * x2 - 2 * x3;
        double N34   = -xi[i] * L * (xi[i] - x2);
        double dN31  = 6 * oneOverL * (x2 - xi[i]); //modified
        double dN32  = 1 + 3 * x2 - 4 * xi[i];
        double dN33  = 6 * oneOverL * (xi[i] - x2);
        double dN34  = 3 * x2 - 2 * xi[i];
        double ddN31 = oneOverLsquare * (xi12 - 6);
        double ddN32 = oneOverL * (xi6 - 4);
        double ddN33 = oneOverLsquare * (6 - xi12);
        double ddN34 = oneOverL * (xi6 - 2);
        double dN12  = oneOverL;

	N1(0,0)  = 1.0;
	N1(0,2)  = (4*v(1)-v(5))/30;
        N1(0,3)  = (4*v(2)-v(6))/30;
        N1(0,5)  = (4*v(5)-v(1))/30;
        N1(0,6)  = (4*v(6)-v(2))/30;
        N1(1,7)  = 1.0;
        N1(1,8)  = N31 * v(0) + N32 * v(3) + N33 * v(4) + N34 * v(7); // Theta
        N1(1,9)  = -ddN32 * v(2) - ddN34 * v(6); // w dprime
        N1(2,7)  = N1(1,8);
        N1(2,8)  = -1.0;
        N1(2,9)  = ddN32 * v(1) + ddN34 * v(5); // v dprime
        N1(3,10) = dN31 * v(0) + dN32 * v(3) + dN33 * v(4) + dN34 * v(7); // Theta prime
        N1(4,11) = 1.0;
        N1(5,10) = 1.0; 

	N2(0,8) = dN12;
        N2(1,0) = 1;
        N2(2,1) = 1;
        N2(3,2) = 1;
        N2(4,4) = 1;
        N2(5,5) = 1;
        N2(6,6) = 1;
        N2(7,1) = ddN32;  // rz1->v dprime
        N2(7,5) = ddN34;  // rz2->v dprime
        N2(8,2) = -ddN32; // ry1->w dprime
        N2(8,6) = -ddN34; // ry2->w dprime
        N2(9,0) = N31;    // rx1->theta
        N2(9,3) = N32;    // rx1 prime->theta
        N2(9,4) = N33;    // rx2->theta
        N2(9,7) = N34;    // rx2 prime->theta
        N2(10,0) = dN31;  // rx1->theta prime
        N2(10,3) = dN32;  // rx1 prime->theta prime
        N2(10,4) = dN33;  // rx2->theta prime
        N2(10,7) = dN34;  // rx2 prime->theta prime
        N2(11,0) = ddN31; // rx1->theta dprime
        N2(11,3) = ddN32; // rx1 prime->theta dprime
        N2(11,4) = ddN33; // rx2->theta dprime
        N2(11,7) = ddN34; // rx2 prime->theta dprime 

	#ifdef DEBUG
	    /*
	    if (counter == 2) {
		for (int a = 0; a < 8; ++a) {
		    for (int b = 0; b < 9; ++b) { 
			//fprintf(stdout, "N1(%d,%d) = %.3f ", i, j, N1(i,j));
			fprintf(stdout, "%.4f ", N1(a,b));
		    }
		    fprintf(stdout, "\n");
		}
	    }
	    */
	#endif

        // Get the section tangent stiffness and stress resultant
        const Matrix & ks = theSections[i]->getInitialTangent();

        N3.addMatrixTripleProduct(0.0, N1, ks, 1.0);
        kbPart1.addMatrixTripleProduct(0.0, N2, N3, 1.0);
        const Vector & s = theSections[i]->getStressResultant();

	#ifdef DEBUG
	fprintf(stdout, "getInitialBasicStiff\n");
	for(int c = 0; c < 6; ++c) 
	    fprintf(stdout, "The value of s[%d] = %.4f\n", c, s(c));

	fprintf(stdout, "\n");
	#endif

	// s needs to be in the following order:
	// [P Mz My W B Tsv]
	Gmax(2,2)   = (4/30)*s(0);
        Gmax(3,3)   = Gmax(2,2);
        Gmax(5,5)   = Gmax(2,2);
        Gmax(6,6)   = Gmax(2,2);
        Gmax(2,5)   = -s(0)/30;
        Gmax(5,2)   = Gmax(2,5);
        Gmax(3,6)   = Gmax(2,5);
        Gmax(6,3)   = Gmax(2,5);
        Gmax(7,9)   = s(2);
        Gmax(9,7)   = Gmax(7,9);
        Gmax(8,9)   = s(1);
        Gmax(9,8)   = Gmax(8,9);
        Gmax(10,10) = s(3); 
        kbPart2.addMatrixTripleProduct(0.0, N2, Gmax, 1.0); 

	#ifdef DEBUG
	    static int counter2 = 1;
	    /*
	    if (counter == 2) {
		for (int a = 0; a < 8; ++a) {
		    for (int b = 0; b < 9; ++b) { 
			//fprintf(stdout, "N1(%d,%d) = %.3f ", i, j, N1(i,j));
			fprintf(stdout, "%.4f ", N1(a,b));
		    }
		    fprintf(stdout, "\n");
		}
	    }
	    */

	    if (counter2 == 2) {
		for (int a = 0; a < 8; ++a) {
		    for (int b = 0; b < 8; ++b) { 
			//fprintf(stdout, "N1(%d,%d) = %.3f ", i, j, N1(i,j));
			fprintf(stdout, "%.4f ", Gmax(a,b));
		    }
		    fprintf(stdout, "\n");
		}
		fprintf(stdout, "\n");
	    }

	    for(int c = 0; c < 6; ++c) 
		fprintf(stdout, "s[%d] = %.4f\n", c, s(c));

	    ++counter2;
	#endif

        // Perform numerical integration
        double wti = wt[i];

        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 9; k++) {
                kb(j, k) += kbPart1(j, k) * L * wti + kbPart2(j, k) * L * wti;
            }
        }
    }

    return kb;
}

const Matrix &TaperedDispBeamColumnSmoothing3d::getInitialStiff()
{
    const Matrix &kb = this->getInitialBasicStiff();

    // Transform to global stiffness
    K = crdTransf->getInitialGlobalStiffMatrix(kb);

    return K;
}

const Matrix &TaperedDispBeamColumnSmoothing3d::getMass()
{

    K.Zero();

    if (rho == 0.0)
        return K;

    double L = crdTransf->getInitialLength();
    double m = 0.5 * rho * L; // this just applies half half the mass 
                              // each node

    K(0, 0) = K(1, 1) = K(2, 2) = K(7, 7) = K(8, 8) = K(9, 9) = m;

    return K;
}

void TaperedDispBeamColumnSmoothing3d::zeroLoad(void)
{
    Q.Zero();

    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;
    q0[3] = 0.0;
    q0[4] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    p0[3] = 0.0;
    p0[4] = 0.0;

    return;
}

int TaperedDispBeamColumnSmoothing3d::addLoad(ElementalLoad * theLoad,
                                     double loadFactor)
{

    int type;
    const Vector &data = theLoad->getData(type, loadFactor);
    double L = crdTransf->getInitialLength();

    if (type == LOAD_TAG_Beam3dUniformLoad) {
        double wy = data(0) * loadFactor;       // Transverse
        double wz = data(1) * loadFactor;       // Transverse
        double wx = data(2) * loadFactor;       // Axial (+ve from node I to J)

        double Vy = 0.5 * wy * L;
        double Mz = Vy * L / 6.0;       // wy*L*L/12
        double Vz = 0.5 * wz * L;
        double My = Vz * L / 6.0;       // wz*L*L/12
        double P = wx * L;

        // Reactions in basic system
        p0[0] -= P;
        p0[1] -= Vy;
        p0[2] -= Vy;
        p0[3] -= Vz;
        p0[4] -= Vz;

        // Fixed end forces in basic system
        q0[0] -= 0.5 * P;
        q0[1] -= Mz;
        q0[2] += Mz;
        q0[3] += My;
        q0[4] -= My;
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
        double Py = data(0) * loadFactor;
        double Pz = data(1) * loadFactor;
        double N = data(2) * loadFactor;
        double aOverL = data(3);

        if (aOverL < 0.0 || aOverL > 1.0)
            return 0;

        double a = aOverL * L;
        double b = L - a;

        // Reactions in basic system
        p0[0] -= N;
        double V1, V2;
        V1 = Py * (1.0 - aOverL);
        V2 = Py * aOverL;
        p0[1] -= V1;
        p0[2] -= V2;
        V1 = Pz * (1.0 - aOverL);
        V2 = Pz * aOverL;
        p0[3] -= V1;
        p0[4] -= V2;

        double L2 = 1.0 / (L * L);
        double a2 = a * a;
        double b2 = b * b;

        // Fixed end forces in basic system
        q0[0] -= N * aOverL;
        double M1, M2;
        M1 = -a * b2 * Py * L2;
        M2 = a2 * b * Py * L2;
        q0[1] += M1;
        q0[2] += M2;
        M1 = -a * b2 * Pz * L2;
        M2 = a2 * b * Pz * L2;
        q0[3] -= M1;
        q0[4] -= M2;
    }
    else {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::addLoad() -- load type unknown for element with tag: "
            << this->getTag() << endln;
        return -1;
    }

    return 0;
}

int TaperedDispBeamColumnSmoothing3d::
addInertiaLoadToUnbalance(const Vector &accel)
{
    // Check for a quick return
    if (rho == 0.0)
        return 0;

    // Get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
        return -1;
    }

    double L = crdTransf->getInitialLength();
    double m = 0.5 * rho * L;

    // Want to add ( - fact * M R * accel ) to unbalance
    // Take advantage of lumped mass matrix
    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);
    Q(2) -= m * Raccel1(2);
    Q(7) -= m * Raccel2(0);
    Q(8) -= m * Raccel2(1);
    Q(9) -= m * Raccel2(2);

    return 0;
}

const Vector &TaperedDispBeamColumnSmoothing3d::getResistingForce()
{
    #ifdef DEBUG
    static int count = 1;
	fprintf(stdout, "Inside function TaperedDispBeamColumnSmoothing3d::getResistingForce %d-th time\n", count);
	++count; 
    #endif

    double L = crdTransf->getInitialLength();
    double oneOverL       = 1.0 / L;
    double oneOverLsquare = oneOverL / L;
    double oneOverLcube   = oneOverLsquare / L;
    static Matrix N1(6,12);
    //static Matrix N1(6, 8);
    static Matrix N2(12,9);
    //static Matrix N2(8, 9);
    static Matrix N3(12,12);
    static Matrix Gmax(12,12); 

    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);
    const Vector & v = crdTransf->getBasicTrialDisp();
    // Zero for integration
    q.Zero();

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        int order = theSections[i]->getOrder();
        const ID & code = theSections[i]->getType();
        N1.Zero();
        N2.Zero();

	// SHAPE FUNCTIONS
        double xi6 = 6.0 * xi[i];
        double xi12 = 12.0 * xi[i];
        double xi1 = xi[i];
        double x3 = xi[i] * xi[i] * xi[i];
        double x2 = xi[i] * xi[i]; 

      	// SHAPE FUNCTIONS
        // N means shape function, dN means 1st derivative, ddN 2nd 
        // N* - * is the order of the shape function 3-cubic 1-linear
        // N*x -  x is the dof of the shape function
        // xi is bringing in x/L already
        double N31   = 1 - 3 * x2 + 2 * x3;
        double N32   = xi[i] * L * (1 - xi[i]) * (1 - xi[i]);
        double N33   = 3 * x2 - 2 * x3;
        double N34   = -xi[i] * L * (xi[i] - x2);
        double dN31  = 6 * oneOverL * (x2 - xi[i]); //modified
        double dN32  = 1 + 3 * x2 - 4 * xi[i];
        double dN33  = 6 * oneOverL * (xi[i] - x2);
        double dN34  = 3 * x2 - 2 * xi[i];
        double ddN31 = oneOverLsquare * (xi12 - 6);
        double ddN32 = oneOverL * (xi6 - 4);
        double ddN33 = oneOverLsquare * (6 - xi12);
        double ddN34 = oneOverL * (xi6 - 2);
        double dN12  = oneOverL;

        N1(0,0) = 1.0;
        N1(0,2) = (4*v(1)-v(5))/30;
        N1(0,3) = (4*v(2)-v(6))/30;
        N1(0,5) = (4*v(5)-v(1))/30;
        N1(0,6) = (4*v(6)-v(2))/30;
	N1(1,7) = 1.0;
	N1(1,8) = N31*v(0) + N32*v(3) + N33*v(4) + N34*v(7); // Theta	
	N1(1,9) = -(ddN32*v(2) + ddN34 * v(6)); // w dprime
	N1(2,7) = N1(1,8);
	N1(2,8) = -1.0;
	N1(2,9) = ddN32*v(1) + ddN34*v(5); // v dprime
	N1(3,10) = dN31*v(0) + dN32*v(3) + dN33*v(4) + dN34*v(7); // Theta prime
	N1(4,11) = 1.0;
	N1(5,10) = 1.0; 

        N2(0,8) = dN12;
        N2(1,0)  = 1;
        N2(2,1)  = 1;
        N2(3,2)  = 1;
        N2(4,4)  = 1;
        N2(5,5)  = 1;
        N2(6,6)  = 1;
        N2(7,1)  = ddN32;  // rz1->v dprime
	N2(7,5)  = ddN34;  // rz2->v dprime
	N2(8,2)  = -ddN32; // ry1->w dprime
	N2(8,6)  = -ddN34; // ry2->w dprime
	N2(9,0)  = N31;    // rx1->theta
	N2(9,3)  = N32;    // rx1 prime->theta
	N2(9,4)  = N33;    // rx2->theta
	N2(9,7)  = N34;    // rx2 prime->theta
	N2(10,0) = dN31;   // rx1->theta prime
	N2(10,3) = dN32;   // rx1 prime->theta prime
	N2(10,4) = dN33;   // rx2->theta prime
	N2(10,7) = dN34;   // rx2 prime->theta prime
	N2(11,0) = ddN31;  // rx1->theta dprime
	N2(11,3) = ddN32;  // rx1 prime->theta dprime
	N2(11,4) = ddN33;  // rx2->theta dprime
	N2(11,7) = ddN34;  // rx2 prime->theta dprime 

	#ifdef DEBUG
	    /*
	    static int counter = 1;
	    // Print Statements for N1
	    fprintf(stdout, "Section %d: N1 Matrix::getResistingForce\n", i);
	    for (int i = 0; i < 6; ++i) {
	    	for (int j = 0; j < 8; ++j) { 
    		    //fprintf(stdout, "N2(%d,%d) = %.4f ", i, j, N2(i,j));
		    fprintf(stdout, " %5.2e ", N1(i,j));
		}
		fprintf(stdout, "\n");
	    }
	    fprintf(stdout, "\n");

	    // Print Statements for N2
	    fprintf(stdout, "N2 matrix in Section %d\n", i);
	    for (int i = 0; i < 8; ++i) {
    		for (int j = 0; j < 9; ++j) { 
		    //fprintf(stdout, "N2(%d,%d) = %.4f ", i, j, N2(i,j));
		    fprintf(stdout, " %.3f ", N2(i,j));
		}
		fprintf(stdout, "\n");
	    }
	    fprintf(stdout, "\n");
	    ++counter;
	    */
	#endif
	
        // Get section stress resultant
        const Vector &s = theSections[i]->getStressResultant();
        //opserr<<"s"<<s;

	#ifdef DEBUG
	fprintf(stdout, "getResistingForce\n");
	for(int c = 0; c < 6; ++c) 
	    fprintf(stdout, "The value of s[%d] = %.4f\n", c, s(c));

	fprintf(stdout, "\n");
	#endif

        double wti = wt[i];

        static Vector qProduct1(12);
        static Vector qProduct2(9);
        qProduct1.Zero();
        qProduct2.Zero();
        qProduct1.addMatrixTransposeVector(0.0, N1, s, 1.0);
        qProduct2.addMatrixTransposeVector(0.0, N2, qProduct1, 1.0);

        for (int j = 0; j < 9; j++) {
            q(j) += qProduct2(j) * L * wti;
        }
    }

    // Transform forces
    Vector p0Vec(p0, 5);
    P = crdTransf->getGlobalResistingForce(q, p0Vec);

    // Subtract other external nodal loads ... P_res = P_int - P_ext
    P.addVector(1.0, Q, -1.0);

    return P;
}

const Vector & TaperedDispBeamColumnSmoothing3d::getResistingForceIncInertia()
{

    this->getResistingForce();

    if (rho != 0.0) {
        const Vector & accel1 = theNodes[0]->getTrialAccel();
        const Vector & accel2 = theNodes[1]->getTrialAccel();

        // Compute the current resisting force
        this->getResistingForce();

        double L = crdTransf->getInitialLength();
        double m = 0.5 * rho * L;

        P(0) += m * accel1(0);
        P(1) += m * accel1(1);
        P(2) += m * accel1(2);
        P(7) += m * accel2(0);
        P(8) += m * accel2(1);
        P(9) += m * accel2(2);

        // add the damping forces if rayleigh damping
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0
            || betaKc != 0.0)
            P += this->getRayleighDampingForces();

    }
    else {

        // add the damping forces if rayleigh damping
        if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            P += this->getRayleighDampingForces();
    }

    
    return P;
}

int TaperedDispBeamColumnSmoothing3d::sendSelf(int commitTag,
                                      Channel & theChannel)
{

    // place the integer data into an ID

    int dbTag = this->getDbTag();
    int i, j;
    int loc = 0;

    static ID idData(7);        // one bigger than needed so no clash later
    idData(0) = this->getTag();
    idData(1) = connectedExternalNodes(0);
    idData(2) = connectedExternalNodes(1);
    idData(3) = numSections;
    idData(4) = crdTransf->getClassTag();
    int crdTransfDbTag = crdTransf->getDbTag();
    if (crdTransfDbTag == 0) {
        crdTransfDbTag = theChannel.getDbTag();
        if (crdTransfDbTag != 0)
            crdTransf->setDbTag(crdTransfDbTag);
    }
    idData(5) = crdTransfDbTag;

    if (alphaM != 0 || betaK != 0 || betaK0 != 0 || betaKc != 0)
        idData(6) = 1;
    else
        idData(6) = 0;


    if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to send ID data\n";
        return -1;
    }

    if (idData(6) == 1) {
        // send damping coefficients
        static Vector dData(4);
        dData(0) = alphaM;
        dData(1) = betaK;
        dData(2) = betaK0;
        dData(3) = betaKc;
        if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
            opserr <<
                "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to send double data\n";
            return -1;
        }
    }

    // send the coordinate transformation
    if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to send crdTranf\n";
        return -1;
    }


    //
    // send an ID for the sections containing each sections dbTag and classTag
    // if section ha no dbTag get one and assign it
    //

    ID idSections(2 * numSections);
    loc = 0;
    for (i = 0; i < numSections; i++) {
        int sectClassTag = theSections[i]->getClassTag();
        int sectDbTag = theSections[i]->getDbTag();
        if (sectDbTag == 0) {
            sectDbTag = theChannel.getDbTag();
            theSections[i]->setDbTag(sectDbTag);
        }

        idSections(loc) = sectClassTag;
        idSections(loc + 1) = sectDbTag;
        loc += 2;
    }

    if (theChannel.sendID(dbTag, commitTag, idSections) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to send ID data\n";
        return -1;
    }

    //
    // send the sections
    //

    for (j = 0; j < numSections; j++) {
        if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
            opserr << "TaperedDispBeamColumnSmoothing3d::sendSelf() - section "
                << j << "failed to send itself\n";
            return -1;
        }
    }

    return 0;
}

int TaperedDispBeamColumnSmoothing3d::recvSelf(int commitTag,
                                      Channel & theChannel,
                                      FEM_ObjectBroker & theBroker)
{

    //
    // receive the integer data containing tag, numSections and coord transformation info
    //
    int dbTag = this->getDbTag();
    int i;

    static ID idData(7);        // one bigger than needed so no clash with section ID

    if (theChannel.recvID(dbTag, commitTag, idData) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::recvSelf() - failed to recv ID data\n";
        return -1;
    }

    this->setTag(idData(0));
    connectedExternalNodes(0) = idData(1);
    connectedExternalNodes(1) = idData(2);

    int crdTransfClassTag = idData(4);
    int crdTransfDbTag = idData(5);

    if (idData(6) == 1) {
        // recv damping coefficients
        static Vector dData(4);
        if (theChannel.recvVector(dbTag, commitTag, dData) < 0) {
            opserr <<
                "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to recv double data\n";
            return -1;
        }
        alphaM = dData(0);
        betaK = dData(1);
        betaK0 = dData(2);
        betaKc = dData(3);
    }

    // create a new crdTransf object if one needed
    if (crdTransf == 0
        || crdTransf->getClassTag() != crdTransfClassTag) {
        if (crdTransf != 0)
            delete crdTransf;

        crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

        if (crdTransf == 0) {
            opserr << "TaperedDispBeamColumnSmoothing3d::recvSelf() - " <<
                "failed to obtain a CrdTrans object with classTag" <<
                crdTransfClassTag << endln;
            return -2;
        }
    }

    crdTransf->setDbTag(crdTransfDbTag);

    // invoke recvSelf on the crdTransf object
    if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to recv crdTranf\n";
        return -3;
    }

    //
    // recv an ID for the sections containing each sections dbTag and classTag
    //

    ID idSections(2 * idData(3));
    int loc = 0;

    if (theChannel.recvID(dbTag, commitTag, idSections) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::recvSelf() - failed to recv ID data\n";
        return -1;
    }

    //
    // now receive the sections
    //

    if (numSections != idData(3)) {

        //
        // we do not have correct number of sections, must delete the old and create
        // new ones before can recvSelf on the sections
        //

        // delete the old
        if (numSections != 0) {
            for (int i = 0; i < numSections; i++)
                delete theSections[i];
            delete[]theSections;
        }

        // create a new array to hold pointers
        theSections = new SectionForceDeformation *[idData(3)];
        if (theSections == 0) {
            opserr <<
                "TaperedDispBeamColumnSmoothing3d::recvSelf() - out of memory creating sections array of size"
                << idData(3) << endln;
            exit(-1);
        }

        // create a section and recvSelf on it
        numSections = idData(3);
        loc = 0;

        for (i = 0; i < numSections; i++) {
            int sectClassTag = idSections(loc);
            int sectDbTag = idSections(loc + 1);
            loc += 2;
            theSections[i] = theBroker.getNewSection(sectClassTag);
            if (theSections[i] == 0) {
                opserr <<
                    "TaperedDispBeamColumnSmoothing3d::recvSelf() - Broker could not create Section of class type"
                    << sectClassTag << endln;
                exit(-1);
            }
            theSections[i]->setDbTag(sectDbTag);
            if (theSections[i]->
                recvSelf(commitTag, theChannel, theBroker) < 0) {
                opserr <<
                    "TaperedDispBeamColumnSmoothing3d::recvSelf() - section "
                    << i << "failed to recv itself\n";
                return -1;
            }
        }

    }
    else {

        // 
        // for each existing section, check it is of correct type
        // (if not delete old & create a new one) then recvSelf on it
        //

        loc = 0;
        for (i = 0; i < numSections; i++) {
            int sectClassTag = idSections(loc);
            int sectDbTag = idSections(loc + 1);
            loc += 2;

            // check of correct type
            if (theSections[i]->getClassTag() != sectClassTag) {
                // delete the old section[i] and create a new one
                delete theSections[i];
                theSections[i] =
                    theBroker.getNewSection(sectClassTag);
                if (theSections[i] == 0) {
                    opserr <<
                        "TaperedDispBeamColumnSmoothing3d::recvSelf() - Broker could not create Section of class type"
                        << sectClassTag << endln;
                    exit(-1);
                }
            }

            // recvSelf on it
            theSections[i]->setDbTag(sectDbTag);
            if (theSections[i]->
                recvSelf(commitTag, theChannel, theBroker) < 0) {
                opserr <<
                    "TaperedDispBeamColumnSmoothing3d::recvSelf() - section "
                    << i << "failed to recv itself\n";
                return -1;
            }
        }
    }

    return 0;
}

// This may need to be updated, we need to be able to extract 
// forces to get moment axial and shear force diagrams.
void TaperedDispBeamColumnSmoothing3d::Print(OPS_Stream &s, int flag)
{
    s << "\nTaperedDispBeamColumnSmoothing3d, element id:  " << this->
        getTag() << endln;
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tmass density:  " << rho << endln;

    double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    N   = q(8);
    Mz1 = q(1);
    Mz2 = q(5);
    Vy  = (Mz1 + Mz2) * oneOverL;
    My1 = q(2);
    My2 = q(6);
    Vz  = -(My1 + My2) * oneOverL;
    T   = q(4);

    s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
        << -N + p0[0] << ' ' << Mz1 << ' ' << Vy +
        p0[1] << ' ' << My1 << ' ' << Vz +
        p0[3] << ' ' << -T << endln;
    s << "\tEnd 2 Forces (P Mz Vy My Vz T): " << N << ' ' << Mz2 <<
        ' ' << -Vy + p0[2] << ' ' << My2 << ' ' << -Vz +
        p0[4] << ' ' << T << endln;

    //for (int i = 0; i < numSections; i++)
    //theSections[i]->Print(s,flag);
}


int TaperedDispBeamColumnSmoothing3d::displaySelf(Renderer & theViewer,
                                         int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector & end1Crd = theNodes[0]->getCrds();
    const Vector & end2Crd = theNodes[1]->getCrds();

    static Vector v1(3);
    static Vector v2(3);

    if (displayMode >= 0) {
        const Vector & end1Disp = theNodes[0]->getDisp();
        const Vector & end2Disp = theNodes[1]->getDisp();

        for (int i = 0; i < 3; i++) {
            v1(i) = end1Crd(i) + end1Disp(i) * fact;
            v2(i) = end2Crd(i) + end2Disp(i) * fact;
        }
    }
    else {
        int mode = displayMode * -1;
        const Matrix & eigen1 = theNodes[0]->getEigenvectors();
        const Matrix & eigen2 = theNodes[1]->getEigenvectors();
        if (eigen1.noCols() >= mode) {
            for (int i = 0; i < 3; i++) {
                v1(i) = end1Crd(i) + eigen1(i, mode - 1) * fact;
                v2(i) = end2Crd(i) + eigen2(i, mode - 1) * fact;
            }

        }
        else {
            for (int i = 0; i < 3; i++) {
                v1(i) = end1Crd(i);
                v2(i) = end2Crd(i);
            }
        }
    }
    return theViewer.drawLine(v1, v2, 1.0, 1.0);
}

// Might need to update this upcoming function.
Response *TaperedDispBeamColumnSmoothing3d::setResponse(const char **argv,
                                               int argc,
                                               OPS_Stream & output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "TaperedDispBeamColumnSmoothing3d");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes[0]);
    output.attr("node2", connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0], "forces") == 0
        || strcmp(argv[0], "force") == 0
        || strcmp(argv[0], "globalForce") == 0
        || strcmp(argv[0], "globalForces") == 0) {

        output.tag("ResponseType", "Px_1");
        output.tag("ResponseType", "Py_1");
        output.tag("ResponseType", "Pz_1");
        output.tag("ResponseType", "Bx_1");
        output.tag("ResponseType", "Mx_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "Px_2");
        output.tag("ResponseType", "Py_2");
        output.tag("ResponseType", "Pz_2");
        output.tag("ResponseType", "Mx_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");
	output.tag("ResponseType", "Bx_2");


        theResponse = new ElementResponse(this, 1, P);

        // local force -
    }
    else if (strcmp(argv[0], "localForce") == 0
             || strcmp(argv[0], "localForces") == 0) {

// I changed this but I don't know if it worked worth anything
// because I'm not sure what ElementResponse(this, 2, ...) does...2 probably puts P into local coordinates
// and therefore I'd need to see that code.  Otherwise this may crash things BEWARE

        output.tag("ResponseType", "N_1");
        output.tag("ResponseType", "Vy_1");
        output.tag("ResponseType", "Vz_1");
        output.tag("ResponseType", "T_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "Bx_1");
        output.tag("ResponseType", "N_2");
        output.tag("ResponseType", "Vy_2");
        output.tag("ResponseType", "Vz_2");
        output.tag("ResponseType", "T_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");
        output.tag("ResponseType", "Bx_2");

        theResponse = new ElementResponse(this, 2, P);

        // chord rotation -
    }
    else if (strcmp(argv[0], "chordRotation") == 0
             || strcmp(argv[0], "chordDeformation") == 0
             || strcmp(argv[0], "basicDeformation") == 0) {

        output.tag("ResponseType", "eps");
        output.tag("ResponseType", "thetaZ_1");
        output.tag("ResponseType", "thetaZ_2");
        output.tag("ResponseType", "thetaY_1");
        output.tag("ResponseType", "thetaY_2");
        output.tag("ResponseType", "thetaX");

        theResponse = new ElementResponse(this, 3, Vector(6));

        // plastic rotation -
    }
    else if (strcmp(argv[0], "plasticRotation") == 0
             || strcmp(argv[0], "plasticDeformation") == 0) {

        output.tag("ResponseType", "epsP");
        output.tag("ResponseType", "thetaZP_1");
        output.tag("ResponseType", "thetaZP_2");
        output.tag("ResponseType", "thetaYP_1");
        output.tag("ResponseType", "thetaYP_2");
        output.tag("ResponseType", "thetaXP");

        theResponse = new ElementResponse(this, 4, Vector(6));

        // section response -
    }
    else if (strstr(argv[0], "sectionX") != 0) {
        if (argc > 2) {
            float sectionLoc = atof(argv[1]);

            double xi[maxNumSections];
            double L = crdTransf->getInitialLength();
            beamInt->getSectionLocations(numSections, L, xi);

            sectionLoc /= L;

            float minDistance = fabs(xi[0] - sectionLoc);
            int sectionNum = 0;
            for (int i = 1; i < numSections; i++) {
                if (fabs(xi[i] - sectionLoc) < minDistance) {
                    minDistance = fabs(xi[i] - sectionLoc);
                    sectionNum = i;
                }
            }

            output.tag("GaussPointOutput");
            output.attr("number", sectionNum + 1);
            output.attr("eta", xi[sectionNum] * L);

            theResponse =
                theSections[sectionNum]->setResponse(&argv[2],
                                                     argc - 2,
                                                     output);
        }
    }
    else if (strcmp(argv[0], "section") == 0) {
        if (argc > 2) {

            int sectionNum = atoi(argv[1]);
            if (sectionNum > 0 && sectionNum <= numSections) {

                double xi[maxNumSections];
                double L = crdTransf->getInitialLength();
                beamInt->getSectionLocations(numSections, L, xi);

                output.tag("GaussPointOutput");
                output.attr("number", sectionNum);
                output.attr("eta", xi[sectionNum - 1] * L);

                theResponse =
                    theSections[sectionNum - 1]->setResponse(&argv[2],
                                                             argc - 2,
                                                             output);

                output.endTag();
            }
        }
    }

    output.endTag();
    return theResponse;
}

// Update this???
int TaperedDispBeamColumnSmoothing3d::getResponse(int responseID,
                                         Information &eleInfo)
{

    double N, V, M1, M2, T1, T2, Bi1, Bi2;
    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    if (responseID == 1)
        return eleInfo.setVector(this->getResistingForce());

    else if (responseID == 2) {
	// NOTE: INSTEAD OF "P" WE MUST USE "Plocal" -> 14 total local forces
	// Axial
	N = q(8);
        P(0) = -N;
        P(7) = N;

        // Torsion
	T1 = q(0);
	T2 = q(4);
        P(3) = T1;
        P(10) = T2;

        // Moments about z and shears along y
        M1 = q(1);
        M2 = q(5);
        P(5) = M1;
        P(12) = M2;

        V = (M1 + M2) * oneOverL;
        P(1) = V;
        P(8) = -V;

        // Moments about y and shears along z
        M1 = q(2);
        M2 = q(6);
        P(4) = M1;
        P(11) = M2;

        V = -(M1 + M2) * oneOverL;
        P(2) = V;
        P(9) = -V;

	// Bi-Moments about x
	Bi1 = q(3);
	Bi2 = q(7);
	P(6) = Bi1;
	P(13) = Bi2;

        return eleInfo.setVector(P);
    }

    // Chord rotation
    else if (responseID == 3)
        return eleInfo.setVector(crdTransf->getBasicTrialDisp());

    // Plastic rotation
    else if (responseID == 4) {
        static Vector vp(6);
        static Vector ve(6);
        const Matrix & kb = this->getInitialBasicStiff();
        kb.Solve(q, ve);
        vp = crdTransf->getBasicTrialDisp();
        vp -= ve;
        return eleInfo.setVector(vp);
    }

    else
        return -1;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int TaperedDispBeamColumnSmoothing3d::setParameter(const char **argv, int argc,
                                          Parameter &param)
{
    if (argc < 1)
        return -1;

    // If the parameter belongs to the element itself
    if (strcmp(argv[0], "rho") == 0)
        return param.addObject(1, this);

    if (strstr(argv[0], "sectionX") != 0) {
        if (argc < 3)
            return -1;

        float sectionLoc = atof(argv[1]);

        double xi[maxNumSections];
        double L = crdTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        sectionLoc /= L;

        float minDistance = fabs(xi[0] - sectionLoc);
        int sectionNum = 0;
        for (int i = 1; i < numSections; i++) {
            if (fabs(xi[i] - sectionLoc) < minDistance) {
                minDistance = fabs(xi[i] - sectionLoc);
                sectionNum = i;
            }
        }
        return theSections[sectionNum]->setParameter(&argv[2],
                                                     argc - 2, param);
    }
    // If the parameter belongs to a section or lower
    if (strstr(argv[0], "section") != 0) {

        if (argc < 3)
            return -1;

        // Get section and material tag numbers from user input
        int paramSectionTag = atoi(argv[1]);

        // Find the right section and call its setParameter method
        int ok = 0;
        for (int i = 0; i < numSections; i++)
            if (paramSectionTag == theSections[i]->getTag())
                ok +=
                    theSections[i]->setParameter(&argv[2], argc - 2,
                                                 param);

        return ok;
    }

    else if (strstr(argv[0], "integration") != 0) {

        if (argc < 2)
            return -1;

        return beamInt->setParameter(&argv[1], argc - 1, param);
    }

    // Default, send to every object
    int ok = 0;
    for (int i = 0; i < numSections; i++)
        ok += theSections[i]->setParameter(argv, argc, param);
    ok += beamInt->setParameter(argv, argc, param);
    return ok;
}

int TaperedDispBeamColumnSmoothing3d::updateParameter(int parameterID,
                                             Information & info)
{
    if (parameterID == 1) {
        rho = info.theDouble;
        return 0;
    }
    else
        return -1;
}


int TaperedDispBeamColumnSmoothing3d::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;

    return 0;
}

const Matrix &
    TaperedDispBeamColumnSmoothing3d::getKiSensitivity(int gradNumber)
{
    K.Zero();
    return K;
}

const Matrix &
    TaperedDispBeamColumnSmoothing3d::getMassSensitivity(int gradNumber)
{
    K.Zero();
    return K;
}



const Vector &
    TaperedDispBeamColumnSmoothing3d::
getResistingForceSensitivity(int gradNumber)
{
    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);

    // Zero for integration
    static Vector dqdh(6);
    dqdh.Zero();

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        int order = theSections[i]->getOrder();
        const ID & code = theSections[i]->getType();

        //double xi6 = 6.0*pts(i,0);
        double xi6 = 6.0 * xi[i];
        //double wti = wts(i);
        double wti = wt[i];

        // Get section stress resultant gradient
        const Vector & dsdh =
            theSections[i]->getStressResultantSensitivity(gradNumber,
                                                          true);

        // Perform numerical integration on internal force gradient
        double sensi;
        for (int j = 0; j < order; j++) {
            sensi = dsdh(j) * wti;
            switch (code(j)) {
            case SECTION_RESPONSE_P:
                dqdh(0) += sensi;
                break;
            case SECTION_RESPONSE_MZ:
                dqdh(1) += (xi6 - 4.0) * sensi;
                dqdh(2) += (xi6 - 2.0) * sensi;
                break;
            case SECTION_RESPONSE_MY:
                dqdh(3) += (xi6 - 4.0) * sensi;
                dqdh(4) += (xi6 - 2.0) * sensi;
                break;
            case SECTION_RESPONSE_T:
                dqdh(5) += sensi;
                break;
            default:
                break;
            }
        }
    }

    // Transform forces
    static Vector dp0dh(6);     // No distributed loads

    P.Zero();

    //////////////////////////////////////////////////////////////

    if (crdTransf->isShapeSensitivity()) {

        // Perform numerical integration to obtain basic stiffness matrix
        // Some extra declarations
        static Matrix kbmine(6, 6);
        kbmine.Zero();
        q.Zero();

        double tmp;

        int j, k;

        for (int i = 0; i < numSections; i++) {

            int order = theSections[i]->getOrder();
            const ID & code = theSections[i]->getType();

            //double xi6 = 6.0*pts(i,0);
            double xi6 = 6.0 * xi[i];
            //double wti = wts(i);
            double wti = wt[i];

            const Vector & s = theSections[i]->getStressResultant();
            const Matrix & ks = theSections[i]->getSectionTangent();

            Matrix ka(workArea, order, 6);
            ka.Zero();

            double si;
            for (j = 0; j < order; j++) {
                si = s(j) * wti;
                switch (code(j)) {
                case SECTION_RESPONSE_P:
                    q(0) += si;
                    for (k = 0; k < order; k++) {
                        ka(k, 0) += ks(k, j) * wti;
                    }
                    break;
                case SECTION_RESPONSE_MZ:
                    q(1) += (xi6 - 4.0) * si;
                    q(2) += (xi6 - 2.0) * si;
                    for (k = 0; k < order; k++) {
                        tmp = ks(k, j) * wti;
                        ka(k, 1) += (xi6 - 4.0) * tmp;
                        ka(k, 2) += (xi6 - 2.0) * tmp;
                    }
                    break;
                case SECTION_RESPONSE_MY:
                    q(3) += (xi6 - 4.0) * si;
                    q(4) += (xi6 - 2.0) * si;
                    for (k = 0; k < order; k++) {
                        tmp = ks(k, j) * wti;
                        ka(k, 3) += (xi6 - 4.0) * tmp;
                        ka(k, 4) += (xi6 - 2.0) * tmp;
                    }
                    break;
                case SECTION_RESPONSE_T:
                    q(5) += si;
                    for (k = 0; k < order; k++) {
                        ka(k, 5) += ks(k, j) * wti;
                    }
                    break;
                default:
                    break;
                }
            }
            for (j = 0; j < order; j++) {
                switch (code(j)) {
                case SECTION_RESPONSE_P:
                    for (k = 0; k < 6; k++) {
                        kbmine(0, k) += ka(j, k);
                    }
                    break;
                case SECTION_RESPONSE_MZ:
                    for (k = 0; k < 6; k++) {
                        tmp = ka(j, k);
                        kbmine(1, k) += (xi6 - 4.0) * tmp;
                        kbmine(2, k) += (xi6 - 2.0) * tmp;
                    }
                    break;
                case SECTION_RESPONSE_MY:
                    for (k = 0; k < 6; k++) {
                        tmp = ka(j, k);
                        kbmine(3, k) += (xi6 - 4.0) * tmp;
                        kbmine(4, k) += (xi6 - 2.0) * tmp;
                    }
                    break;
                case SECTION_RESPONSE_T:
                    for (k = 0; k < 6; k++) {
                        kbmine(5, k) += ka(j, k);
                    }
                    break;
                default:
                    break;
                }
            }
        }

        const Vector & A_u = crdTransf->getBasicTrialDisp();
        double dLdh = crdTransf->getdLdh();
        double d1overLdh = -dLdh / (L * L);
        // a^T k_s dadh v
        dqdh.addMatrixVector(1.0, kbmine, A_u, d1overLdh);

        // k dAdh u
        const Vector & dAdh_u =
            crdTransf->getBasicTrialDispShapeSensitivity();
        dqdh.addMatrixVector(1.0, kbmine, dAdh_u, oneOverL);

        // dAdh^T q
        P += crdTransf->getGlobalResistingForceShapeSensitivity(q,
                                                                dp0dh,
                                                                gradNumber);
    }

    // A^T (dqdh + k dAdh u)
    P += crdTransf->getGlobalResistingForce(dqdh, dp0dh);

    return P;
}



// NEW METHOD
int TaperedDispBeamColumnSmoothing3d::commitSensitivity(int gradNumber,
                                               int numGrads)
{
    // Get basic deformation and sensitivities
    const Vector & v = crdTransf->getBasicTrialDisp();

    static Vector dvdh(6);
    dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;
    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);

    // Some extra declarations
    double d1oLdh = crdTransf->getd1overLdh();

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        int order = theSections[i]->getOrder();
        const ID & code = theSections[i]->getType();

        Vector e(workArea, order);

        //double xi6 = 6.0*pts(i,0);
        double xi6 = 6.0 * xi[i];

        for (int j = 0; j < order; j++) {
            switch (code(j)) {
            case SECTION_RESPONSE_P:
                e(j) = oneOverL * dvdh(0)
                    + d1oLdh * v(0);
                break;
            case SECTION_RESPONSE_MZ:
                e(j) =
                    oneOverL * ((xi6 - 4.0) * dvdh(1) +
                                (xi6 - 2.0) * dvdh(2))
                    + d1oLdh * ((xi6 - 4.0) * v(1) +
                                (xi6 - 2.0) * v(2));
                break;
            case SECTION_RESPONSE_MY:
                e(j) =
                    oneOverL * ((xi6 - 4.0) * dvdh(3) +
                                (xi6 - 2.0) * dvdh(4))
                    + d1oLdh * ((xi6 - 4.0) * v(3) +
                                (xi6 - 2.0) * v(4));
                break;
            case SECTION_RESPONSE_T:
                e(j) = oneOverL * dvdh(5)
                    + d1oLdh * v(5);
                break;
            default:
                e(j) = 0.0;
                break;
            }
        }

        // Set the section deformations
        theSections[i]->commitSensitivity(e, gradNumber, numGrads);
    }

    return 0;
}
// AddingSensitivity:END /////////////////////////////////////////////
