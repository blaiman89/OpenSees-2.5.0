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

// $Revision: 1.31 $
// $Date: 2009/09/28 22:48:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/TaperedFiberSectionSmoothing3d.cpp,v $

// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of TaperedFiberSectionSmoothing3d.
// Modified by Xi Zhang from University of Sydney, Australia (include warping degrees of freedom). Refer to 
// Formulation and Implementation of Three-dimensional Doubly Symmetric Beam-Column Analyses with Warping Effects in OpenSees
// Research Report R917, School of Civil Engineering, University of Sydney.

//#define DEBUG

#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <TaperedFiberSectionSmoothing3d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <iostream>
#include <fstream>
using std::string;
using namespace std;


ID TaperedFiberSectionSmoothing3d::code(SEC_TAG_TaperedFiberSectionSmoothing3d);

// constructors:
TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d(int tag, int num,
                               Fiber ** fibers, double hratio, double g):SectionForceDeformation(tag,
                                                               SEC_TAG_TaperedFiberSectionSmoothing3d),
numFibers(num), theMaterials(0), matData(0), yBar(0.0), zBar(0.0), G(g), hRatio(hratio),
e(12), eCommit(12), s(0), ks(0)
{
    #ifdef DEBUG
	fprintf(stdout, "Inside function TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d(...)\n");
    #endif

    if (numFibers != 0) {
        theMaterials = new UniaxialMaterial *[numFibers];

        if (theMaterials == 0) {
            opserr <<
                "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to allocate Material pointers\n";
            exit(-1);
        }

        matData = new double[numFibers * 5];

        if (matData == 0) {
            opserr <<
                "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to allocate double array for material data\n";
            exit(-1);
        }

        double Qz = 0.0;
        double Qy = 0.0;
        double A  = 0.0;

        for (int i = 0; i < numFibers; i++) {
            Fiber *theFiber = fibers[i];
            double yLoc, zLoc, Area, tP;
            theFiber->getFiberLocation(yLoc, zLoc);
	    tP = theFiber->gettP();                    // MODIFIED FOR TAPERED STUDY
            Area = theFiber->getArea();

            Qz += yLoc * Area;
            Qy += zLoc * Area;
            A  += Area;

            matData[i * 5]     = yLoc;
            matData[i * 5 + 1] = zLoc;
            matData[i * 5 + 2] = Area;
	    matData[i * 5 + 3] = hRatio;
	    matData[i * 5 + 4] = tP; // plate thickness for the fiber...this assumes fibers are located at mid-thickness only!
            UniaxialMaterial *theMat = theFiber->getMaterial();
            theMaterials[i] = theMat->getCopy();

            if (theMaterials[i] == 0) {
                opserr <<
                    "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to get copy of a Material\n";
                exit(-1);
            }
        }

        yBar = Qz / A;
        zBar = Qy / A;
    }

    s  = new Vector(sData, 6);
    ks = new Matrix(kData, 6, 6);

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;
    sData[3] = 0.0;
    sData[4] = 0.0;
    sData[5] = 0.0;

    for (int i = 0; i < 36; i++)
        kData[i] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;

    // AddingSensitivity:BEGIN ////////////////////////////////////
    parameterID = 0;
    SHVs = 0;
    // AddingSensitivity:END //////////////////////////////////////
}

// constructor for blank object that recvSelf needs to be invoked upon
TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d():
SectionForceDeformation(0, SEC_TAG_TaperedFiberSectionSmoothing3d),
numFibers(0), theMaterials(0), matData(0),
yBar(0.0), zBar(0.0), G(0.0), hRatio(0.0), e(12), eCommit(4), s(0), ks(0)
{
    #ifdef DEBUG
	fprintf(stdout, "Inside function TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d()\n");
    #endif
    s = new Vector(sData, 6);
    ks = new Matrix(kData, 6, 6);

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;
    sData[3] = 0.0;
    sData[4] = 0.0;
    sData[5] = 0.0;

    for (int i = 0; i < 36; i++)
        kData[i] = 0.0;

	// what does this next piece do?
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;

    // AddingSensitivity:BEGIN ////////////////////////////////////
    parameterID = 0;
    SHVs = 0;
    // AddingSensitivity:END //////////////////////////////////////
}

// Looks like we do not really care much about this since it doesn't get touched
int TaperedFiberSectionSmoothing3d::addFiber(Fiber & newFiber)
{
    #ifdef DEBUG
	fprintf(stdout, "Inside function TaperedFiberSectionSmoothing3d::addFiber\n");
    #endif
    // need to create a larger array
    int newSize = numFibers + 1;

    UniaxialMaterial **newArray = new UniaxialMaterial *[newSize];
    double *newMatData = new double[5 * newSize];

    if (newArray == 0 || newMatData == 0) {
        opserr <<
            "TaperedFiberSectionSmoothing3d::addFiber -- failed to allocate Fiber pointers\n";
        exit(-1);
    }

    // copy the old pointers
    int i;
    for (i = 0; i < numFibers; i++) {
        newArray[i] = theMaterials[i];
        newMatData[5 * i] = matData[5 * i];
        newMatData[5 * i + 1] = matData[5 * i + 1];
        newMatData[5 * i + 2] = matData[5 * i + 2];
        newMatData[5 * i + 3] = matData[5 * i + 3];
	newMatData[5 * i + 4] = matData[5 * i + 4];
    }
    // set the new pointers
    double yLoc, zLoc, Area, hRatio, tP;
    newFiber.getFiberLocation(yLoc, zLoc);
    tP = newFiber.gettP();
    Area = newFiber.getArea();


    newMatData[numFibers * 5]     = yLoc;
    newMatData[numFibers * 5 + 1] = zLoc;
    newMatData[numFibers * 5 + 2] = Area;
    newMatData[numFibers * 5 + 3] = hRatio;
    newMatData[numFibers * 5 + 4] = tP;

    UniaxialMaterial *theMat = newFiber.getMaterial();
    newArray[numFibers] = theMat->getCopy();

    if (newArray[numFibers] == 0) {
        opserr <<
            "TaperedFiberSectionSmoothing3d::addFiber -- failed to get copy of a Material\n";
        exit(-1);

        delete[]newArray;
        delete[]newMatData;
        return -1;
    }

    numFibers++;

    if (theMaterials != 0) {
        delete[]theMaterials;
        delete[]matData;
    }

    theMaterials = newArray;
    matData = newMatData;

    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
        yLoc   = matData[5 * i];
        zLoc   = matData[5 * i + 1];
        Area   = matData[5 * i + 2];

        A  += Area;
        Qz += yLoc * Area;
        Qy += zLoc * Area;
    }

    yBar = Qz / A;
    zBar = Qy / A;

    return 0;
}



// destructor:
TaperedFiberSectionSmoothing3d::~TaperedFiberSectionSmoothing3d()
{
    if (theMaterials != 0) {
        for (int i = 0; i < numFibers; i++)
            if (theMaterials[i] != 0)
                delete theMaterials[i];

        delete[]theMaterials;
    }

    if (matData != 0)
        delete[]matData;

    if (s != 0)
        delete s;

    if (ks != 0)
        delete ks;
}

int TaperedFiberSectionSmoothing3d::setTrialSectionDeformation(const Vector &deforms)
{
    #ifdef DEBUG
	fprintf(stdout, "Inside function TaperedFiberSectionSmoothing3d::setTrialSectionDeformation\n");
	fprintf(stdout, "The number of fibers = %d\n", numFibers);
	fprintf(stdout, "The tag is = %d\n", code);
	fprintf(stdout, "The number of G = %lf\n", G);
	fprintf(stdout, "The value of hRatio = %lf\n", hRatio);
	fprintf(stdout, "The value of yBar = %lf\n", yBar);
	fprintf(stdout, "The value of zBar = %lf\n", zBar);
    #endif
    int res = 0;
    e = deforms; // this is the section deformations

    kData[0] = 0.0;
    kData[1] = 0.0;
    kData[2] = 0.0;
    kData[3] = 0.0;
    kData[4] = 0.0;
    kData[5] = 0.0;
    kData[6] = 0.0;
    kData[7] = 0.0;
    kData[8] = 0.0;
    kData[9] = 0.0;
    kData[10] = 0.0;
    kData[11] = 0.0;
    kData[12] = 0.0;
    kData[13] = 0.0;
    kData[14] = 0.0;
    kData[15] = 0.0;
    kData[16] = 0.0;
    kData[17] = 0.0;
    kData[18] = 0.0;
    kData[19] = 0.0;
    kData[20] = 0.0;
    kData[21] = 0.0;
    kData[22] = 0.0;
    kData[23] = 0.0;
    kData[24] = 0.0;
    kData[25] = 0.0;
    kData[26] = 0.0;
    kData[27] = 0.0;
    kData[28] = 0.0;
    kData[29] = 0.0;
    kData[30] = 0.0;
    kData[31] = 0.0;
    kData[32] = 0.0;
    kData[33] = 0.0;
    kData[34] = 0.0;
    kData[35] = 0.0; //GJ changing to dealing with tP and G; // Originally 0.0

    
    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;
    sData[3] = 0.0;
    sData[4] = 0.0;
    sData[5] = 0.0;

    int loc = 0;

    // deforms is a vector of generalized strains
    // this is e coming from DispBeamColumn3d
    // order is [u', v', w', v'', w'', theta, theta', theta''] 
    double d0 = deforms(0);
    double d1 = deforms(1);
    double d2 = deforms(2);
    double d3 = deforms(3);
    double d4 = deforms(4);
    double d5 = deforms(5);
    double d6 = deforms(6);
    double d7 = deforms(7);
    // For smoothing add these in:
    // for smoothing the order is [u' rx1 rz1 ry1 rx2 rz2 ry2 v'' w'' theta theta' theta'']
    double d8 = deforms(8);
    double d9 = deforms(9);
    double d10 = deforms(10);
    double d11 = deforms(11);


    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat = theMaterials[i];
        double y = matData[loc++];
        double z = matData[loc++];
        double A = matData[loc++];
        double yP = y * matData[loc++]; //where matData[loc++] is hRatio
	double tP = matData[loc++];     // I want this to be plate thickness

	// calculate psi function
	double psi = 0.0;
	psi = 2 * yP * z;
	
        // calculate sectorial area
        double omig = 0.0;
	omig = y * z; // NOTE: This has to have the right sign from input

	#ifdef DEBUG
	    //fprintf(stdout, "The value of psi in setTrialSectionDeformation = %.4f\n", psi);
	    //fprintf(stdout, "The value of y in setTrialSectionDeformation = %.4f\n", y);
	    //fprintf(stdout, "The value of z in setTrialSectionDeformation = %.4f\n", z);
	    //fprintf(stdout, "The value of Area in setTrialSectionDeformation = %.4f\n", A);
	    //fprintf(stdout, "The value of yP in setTrialSectionDeformation = %.4f\n", yP);
	    //fprintf(stdout, "The value of omig in setTrialSectionDeformation = %.4f\n", omig);
	    //fprintf(stdout, "The value of psi in setTrialSectionDeformation = %.4f\n", psi);
	#endif

        // determine material strain and set it, include second order terms
	// this is only the axial strain, shear strain to be handled seperately
	// NOTE: These are taking into consideration smoothing.
	double strain = d0 - y * d7 - z * d8 -omig*d11 - psi * d10 + z * d7 * d9
		-y * d8 * d9 + 0.5 * (y * y + z * z) * d10 * d10 + d3 * (d3/15 - d6/60)
		+ d6 * (d6/15 - d3/60) + d2 * (d2/15 - d5/60) + d5 * (d5/15 - d2/60);

	double gamma =
		-tP * e[10];

	/*
	 * NOTE: For normal TaperedFiberSectionSmoothing3d.
        double strain =
            d0 - y * d3 - z * d4 - omig * d7 + 0.5 * (y * y + z * z) *
	    d6 * d6 - y * d5 * d4 + z * d5 * d3 - psi * d6 + 0.5*(d1 * d1 + d2 * d2);
		
	// here is the shear strain, max for the plate edge
	double gamma = -tP * e[6];
	*/

	// If smoothing is used, the nodal rotations are going to have to 
	// be retrieved from  some other class to redefine strain as:
	// deforms has 12 entries
	// [u' rx1 rz1 ry1 rx2 rz2 ry2 v'' w'' theta theta' theta'']
	/*
	double strain = d0 - y * d7 - z * d8 -omig*d11 - psi * d10 + z * d7 * d9
		-y * d8 * d9 + 0.5 * (y * y + z * z) * d10 * d10 + d3 * (d3/15 - d6/60)
		+ d6 * (d6/15 - d3/60) + d2 * (d2/15 - d5/60) + d5 * (d5/15 - d2/60);

	double gamma =
		-tP * d10;
	*/

	// retrieve axial tangent modulus and stress
        double tangent, stress;
        res += theMat->setTrial(strain, stress, tangent);

	// shear assumed to behave elastically always, J2 plasticity can be added if needed
	double tau = G * gamma;   // NEED TO ADD G, same for whole section and element

        double value = tangent * A;
        double vas1 = y * value;
        double vas2 = z * value;
        double vas1as2 = vas1 * z;

        // section stiffness matrix k, refer to Alemdar

	kData[0]  += value;
	kData[1]  += -1 * vas1;
	kData[2]  += vas2;
	kData[3]  += (y * y + z * z) * value;
	kData[4]  += -omig * value;
	kData[5]  += -psi * value;
	kData[6]  = kData[1];
	kData[7]  += y * vas1;
	kData[8]  += -1 * vas1as2;
	kData[9]  += -1 * vas1 * (y * y + z * z);
	kData[10] += vas1 * omig;
	kData[11] += vas1 * psi;
	kData[12] = kData[2];
	kData[13] = kData[8];
	kData[14] += z * vas2; 
	kData[15] += vas2 * (y * y + z * z);
	kData[16] += -vas2 * omig;
	kData[17] += -vas2 * psi;
	kData[18] = kData[3];
	kData[19] = kData[9];
	kData[20] = kData[15];
	kData[21] += (y * y + z * z) * (y * y + z * z) * value;
	kData[22] += -omig * value * (y * y + z * z);
	kData[23] += -psi * value * (y * y + z * z);
	kData[24] = kData[4];
	kData[25] = kData[10];
	kData[26] = kData[16];
	kData[27] = kData[22];
	kData[28] += omig * omig * value;
	kData[29] += omig * psi * value;
	kData[30] = kData[5];
	kData[31] = kData[11];
	kData[32] = kData[17];
	kData[33] = kData[23];
	kData[34] = kData[29];
	kData[35] += psi * psi * value + (G * A * tP * tP)/3;  // the second term replaces GJ

        // section force vector D, refer to Alemdar
        double fs0 = stress * A;
        sData[0] += fs0;
        sData[1] += -1.0 * fs0 * y;
        sData[2] += 1.0 * fs0 * z;
        sData[3] += fs0 * (y * y + z * z);
        sData[4] += -fs0 * omig;
	sData[5] += -fs0 * psi - (tau * A * tP)/3;  // second term adds pure st. venant torsion, it was missing before!!!

	}
	#ifdef DEBUG
	   /* 
    	    fprintf(stdout, "The value of kData[0] = %.2f\n", kData[0]);
    	    fprintf(stdout, "The value of kData[1] = %.2f\n", kData[1]);
    	    fprintf(stdout, "The value of kData[2] = %.2f\n", kData[2]);
    	    fprintf(stdout, "The value of kData[3] = %.2f\n", kData[3]);
    	    fprintf(stdout, "The value of kData[4] = %.2f\n", kData[4]);
    	    fprintf(stdout, "The value of kData[5] = %.2f\n", kData[5]);
    	    fprintf(stdout, "The value of kData[6] = %.2f\n", kData[6]);
    	    fprintf(stdout, "The value of kData[7] = %.2f\n", kData[7]);
    	    fprintf(stdout, "The value of kData[8] = %.2f\n", kData[8]);
    	    fprintf(stdout, "The value of kData[9] = %.2f\n", kData[9]);
    	    fprintf(stdout, "The value of kData[10] = %.2f\n", kData[10]);
    	    fprintf(stdout, "The value of kData[11] = %.2f\n", kData[11]);
    	    fprintf(stdout, "The value of kData[12] = %.2f\n", kData[12]);
    	    fprintf(stdout, "The value of kData[13] = %.2f\n", kData[13]);
    	    fprintf(stdout, "The value of kData[14] = %.2f\n", kData[14]);
    	    fprintf(stdout, "The value of kData[15] = %.2f\n", kData[15]);
    	    fprintf(stdout, "The value of kData[16] = %.2f\n", kData[16]);
    	    fprintf(stdout, "The value of kData[17] = %.2f\n", kData[17]);
    	    fprintf(stdout, "The value of kData[18] = %.2f\n", kData[18]);
    	    fprintf(stdout, "The value of kData[19] = %.2f\n", kData[19]);
    	    fprintf(stdout, "The value of kData[20] = %.2f\n", kData[20]);
    	    fprintf(stdout, "The value of kData[21] = %.2f\n", kData[21]);
    	    fprintf(stdout, "The value of kData[22] = %.2f\n", kData[22]);
    	    fprintf(stdout, "The value of kData[23] = %.2f\n", kData[23]);
    	    fprintf(stdout, "The value of kData[24] = %.2f\n", kData[24]);
    	    fprintf(stdout, "The value of kData[25] = %.2f\n", kData[25]);
    	    fprintf(stdout, "The value of kData[26] = %.2f\n", kData[26]);
    	    fprintf(stdout, "The value of kData[27] = %.2f\n", kData[27]);
    	    fprintf(stdout, "The value of kData[28] = %.2f\n", kData[28]);
    	    fprintf(stdout, "The value of kData[29] = %.2f\n", kData[29]);
    	    fprintf(stdout, "The value of kData[30] = %.2f\n", kData[30]);
    	    fprintf(stdout, "The value of kData[31] = %.2f\n", kData[31]);
    	    fprintf(stdout, "The value of kData[32] = %.2f\n", kData[32]);
    	    fprintf(stdout, "The value of kData[33] = %.2f\n", kData[33]);
    	    fprintf(stdout, "The value of kData[34] = %.2f\n", kData[34]);
    	    fprintf(stdout, "The value of kData[35] = %.2f\n", kData[35]);
	    */
	#endif


    return res;
}

const Matrix &TaperedFiberSectionSmoothing3d::getInitialTangent(void)
{
    #ifdef DEBUG
	fprintf(stdout, "Inside function TaperedFiberSectionSmoothing3d::getInitialTangent\n");
    #endif
    
    static double kInitialData[36];
    static Matrix kInitial(kInitialData, 6, 6);
    for (int i = 0; i < 36; i++)
        kInitialData[i] = 0.0;

    int loc = 0;

    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat = theMaterials[i];
        double y  = matData[loc++];
        double z  = matData[loc++];
        double A  = matData[loc++];
        double yP = y * matData[loc++];
	double tP = matData[loc++];

        // calculate sectorial area
        double omig = 0.0;
	omig = y * z;

	// calculate psi function
	double psi = 0.0;
	psi = 2 * yP * z;
	#ifdef DEBUG
	    fprintf(stdout, "The value of psi in getInitialTangent = %.4f\n", psi);
	#endif
	

        double tangent = theMat->getInitialTangent();

        double value = tangent * A;
        double vas1  = y * value;
        double vas2  = z * value;
        double vas1as2 = vas1 * z;

        // section stiffness matrix k, refer to Alemdar
	kInitialData[0]  += value;
	kInitialData[1]  += -1*vas1;
	kInitialData[2]  += vas2;
	kInitialData[3]  += (y*y+z*z)*value;
	kInitialData[4]  += -omig*value;
	kInitialData[5]  += -psi*value;
	kInitialData[6]  = kInitialData[1];
	kInitialData[7]  += y*vas1;
	kInitialData[8]  += -1*vas1as2;
	kInitialData[9]  += -1*vas1*(y*y+z*z);
	kInitialData[10] += vas1*omig;
	kInitialData[11] += vas1*psi;
	kInitialData[12] = kInitialData[2];
	kInitialData[13] = kInitialData[8];
	kInitialData[14] += z*vas2; 
	kInitialData[15] += vas2*(y*y+z*z);
	kInitialData[16] += -vas2*omig;
	kInitialData[17] += -vas2*psi;
	kInitialData[18] = kInitialData[3];
	kInitialData[19] = kInitialData[9];
	kInitialData[20] = kInitialData[15];
	kInitialData[21] += (y*y+z*z)*(y*y+z*z)*value;
	kInitialData[22] += -omig*value*(y*y+z*z);
	kInitialData[23] += -psi*value*(y*y+z*z);
	kInitialData[24] = kInitialData[4];
	kInitialData[25] = kInitialData[10];
	kInitialData[26] = kInitialData[16];
	kInitialData[27] = kInitialData[22];
	kInitialData[28] += omig*omig*value;
	kInitialData[29] += omig*psi*value;
	kInitialData[30] = kInitialData[5];
	kInitialData[31] = kInitialData[11];
	kInitialData[32] = kInitialData[17];
	kInitialData[33] = kInitialData[23];
	kInitialData[34] = kInitialData[29];
	kInitialData[35] += psi * psi * value + (G * A * tP * tP)/3;
    }


    return kInitial; // Initially was only "kInitial"
}

const Vector & TaperedFiberSectionSmoothing3d::getSectionDeformation(void)
{
    return e;
}

const Matrix & TaperedFiberSectionSmoothing3d::getSectionTangent(void)
{
    return *ks;
}

const Vector & TaperedFiberSectionSmoothing3d::getStressResultant(void)
{
    return *s;
}

SectionForceDeformation *TaperedFiberSectionSmoothing3d::getCopy(void)
{
    #ifdef DEBUG
	fprintf(stdout, "Inside function TaperedFiberSectionSmoothing3d::getCopy\n");
    #endif
    TaperedFiberSectionSmoothing3d *theCopy = new TaperedFiberSectionSmoothing3d();
    theCopy->setTag(this->getTag());

    theCopy->numFibers = numFibers;

    if (numFibers != 0) {
        theCopy->theMaterials = new UniaxialMaterial *[numFibers];

        if (theCopy->theMaterials == 0) {
            opserr <<
                "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to allocate Material pointers\n";
            exit(-1);
        }

        theCopy->matData = new double[numFibers * 5];

        if (theCopy->matData == 0) {
            opserr <<
                "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to allocate double array for material data\n";
            exit(-1);
        }


        for (int i = 0; i < numFibers; i++) {
            theCopy->matData[i * 5]     = matData[i * 5];
            theCopy->matData[i * 5 + 1] = matData[i * 5 + 1];
            theCopy->matData[i * 5 + 2] = matData[i * 5 + 2];
            theCopy->matData[i * 5 + 3] = matData[i * 5 + 3];
            theCopy->matData[i * 5 + 4] = matData[i * 5 + 4];
            theCopy->theMaterials[i] = theMaterials[i]->getCopy();

            if (theCopy->theMaterials[i] == 0) {
                opserr <<
                    "TaperedFiberSectionSmoothing3d::getCopy -- failed to get copy of a Material\n";
                exit(-1);
            }
        }
    }

    theCopy->eCommit = eCommit;
    theCopy->e = e;
    theCopy->G = G;
    theCopy->hRatio = hRatio;
    theCopy->yBar = yBar;
    theCopy->zBar = zBar;

    for (int i = 0; i < 36; i++)
        theCopy->kData[i] = kData[i];

    /*
    #ifdef DEBUG
	for (int i = 0; i < 36; i++) {
	fprintf(stdout, "The value of kData[%d] = %.4f\n", i, kData[i]);
	}
    #endif
    */

    theCopy->sData[0] = sData[0];
    theCopy->sData[1] = sData[1];
    theCopy->sData[2] = sData[2];
    theCopy->sData[3] = sData[3];
    theCopy->sData[4] = sData[4];
    theCopy->sData[5] = sData[5];

    return theCopy;
}

const ID & TaperedFiberSectionSmoothing3d::getType()
{
    return code;
}

int TaperedFiberSectionSmoothing3d::getOrder() const const
{
    return 6; // Initially was a return value of "5"
}

int TaperedFiberSectionSmoothing3d::commitState(void)
{
    #ifdef DEBUG
	fprintf(stdout, "In function TaperedFiberSectionSmoothing3d::commitState\n");
    #endif
    
    int err = 0;

    for (int i = 0; i < numFibers; i++)
        err += theMaterials[i]->commitState();

    eCommit = e;

    return err;
}

int TaperedFiberSectionSmoothing3d::revertToLastCommit(void)
{
    #ifdef DEBUG
	fprintf(stdout, "In function TaperedFiberSectionSmoothing3d::revertToLastCommit\n");
    #endif
    
    int err = 0;

    // Last committed section deformations
    e = eCommit;


    kData[0] = 0.0;
    kData[1] = 0.0;
    kData[2] = 0.0;
    kData[3] = 0.0;
    kData[4] = 0.0;
    kData[5] = 0.0;
    kData[6] = 0.0;
    kData[7] = 0.0;
    kData[8] = 0.0;
    kData[9] = 0.0;
    kData[10] = 0.0;
    kData[11] = 0.0;
    kData[12] = 0.0;
    kData[13] = 0.0;
    kData[14] = 0.0;
    kData[15] = 0.0;
    kData[16] = 0.0;
    kData[17] = 0.0;
    kData[18] = 0.0;
    kData[19] = 0.0;
    kData[20] = 0.0;
    kData[21] = 0.0;
    kData[22] = 0.0;
    kData[23] = 0.0;
    kData[24] = 0.0;
    kData[25] = 0.0;
    kData[26] = 0.0;
    kData[27] = 0.0;
    kData[28] = 0.0;
    kData[29] = 0.0;
    kData[30] = 0.0;
    kData[31] = 0.0;
    kData[32] = 0.0;
    kData[33] = 0.0;
    kData[34] = 0.0;
    kData[35] = 0; //GJ; //modified from 0.0

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;
    sData[3] = 0.0;
    sData[4] = 0.0;
    sData[5] = 0.0;


    int loc = 0;

    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat = theMaterials[i];
        double y  = matData[loc++];
        double z  = matData[loc++];
        double A  = matData[loc++];
        double yP = y * matData[loc++];
        double tP = matData[loc++];

        // invoke revertToLast on the material
        err += theMat->revertToLastCommit();

        double tangent = theMat->getTangent();
        double stress = theMat->getStress();

	double gamma = -tP * e[10];
	double tau = G * gamma;   // NEED TO ADD G, same for whole section and element

        double value = tangent * A;
        double vas1 = y * value;
        double vas2 = z * value;
        double vas1as2 = vas1 * z;
        double omig = 0.0;
	omig = y * z;

	// calculate psi function
	double psi = 0.0;
	psi = 2 * yP * z;
	#ifdef DEBUG
	    fprintf(stdout, "The value of psi in revertToLastCommit = %.4f\n", psi);
	#endif

	kData[0]  += value;
	kData[1]  += -1*vas1;
	kData[2]  += vas2;
	kData[3]  += (y*y+z*z)*value;
	kData[4]  += -omig*value;
	kData[5]  += -psi*value;
	kData[6]  = kData[1];
	kData[7]  += y*vas1;
	kData[8]  += -1*vas1as2;
	kData[9]  += -1*vas1*(y*y+z*z);
	kData[10] += vas1*omig;
	kData[11] += vas1*psi;
	kData[12] = kData[2];
	kData[13] = kData[8];
	kData[14] += z*vas2; 
	kData[15] += vas2*(y*y+z*z);
	kData[16] += -vas2*omig;
	kData[17] += -vas2*psi;
	kData[18] = kData[3];
	kData[19] = kData[9];
	kData[20] = kData[15];
	kData[21] += (y*y+z*z)*(y*y+z*z)*value;
	kData[22] += -omig*value*(y*y+z*z);
	kData[23] += -psi*value*(y*y+z*z);
	kData[24] = kData[4];
	kData[25] = kData[10];
	kData[26] = kData[16];
	kData[27] = kData[22];
	kData[28] += omig*omig*value;
	kData[29] += omig*psi*value;
	kData[30] = kData[5];
	kData[31] = kData[11];
	kData[32] = kData[17];
	kData[33] = kData[23];
	kData[34] = kData[29];
	kData[35] += psi*psi*value+ (G * A * tP * tP)/3;

        double fs0 = stress * A;  // I'm not sure where stress came from but you need tau like above as well for sData[5]

        sData[0] += fs0;
        sData[1] += -1.0 * fs0 * y;
        sData[2] += 1.0 * fs0 * z;
        sData[3] += fs0 * (y * y + z * z);
        sData[4] += -fs0 * omig;
	sData[5] += -fs0 * psi - (tau * A * tP)/3;
    }

    return err;
}

int TaperedFiberSectionSmoothing3d::revertToStart(void)
{
    // revert the fibers to start    
    int err = 0;


    kData[0] = 0.0;
    kData[1] = 0.0;
    kData[2] = 0.0;
    kData[3] = 0.0;
    kData[4] = 0.0;
    kData[5] = 0.0;
    kData[6] = 0.0;
    kData[7] = 0.0;
    kData[8] = 0.0;
    kData[9] = 0.0;
    kData[10] = 0.0;
    kData[11] = 0.0;
    kData[12] = 0.0;
    kData[13] = 0.0;
    kData[14] = 0.0;
    kData[15] = 0.0;
    kData[16] = 0.0;
    kData[17] = 0.0;
    kData[18] = 0.0;
    kData[19] = 0.0;
    kData[20] = 0.0;
    kData[21] = 0.0;
    kData[22] = 0.0;
    kData[23] = 0.0;
    kData[24] = 0.0;
    kData[25] = 0.0;
    kData[26] = 0.0;
    kData[27] = 0.0;
    kData[28] = 0.0;
    kData[29] = 0.0;
    kData[30] = 0.0;
    kData[31] = 0.0;
    kData[32] = 0.0;
    kData[33] = 0.0;
    kData[34] = 0.0;
    kData[35] = 0.0; // modified from 0.0

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;
    sData[3] = 0.0;
    sData[4] = 0.0;
    sData[5] = 0.0;

    int loc = 0;

    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat = theMaterials[i];
        double y = matData[loc++];
        double z = matData[loc++];
        double A = matData[loc++];
        double yP = y * matData[loc++];
        double tP = matData[loc++];

        double omig = 0.0;
	omig = y * z;

	// calculate psi function 
	double psi = 0.0;
	psi = 2 * yP * z;
	#ifdef DEBUG
	    fprintf(stdout, "The value of psi in revertToStart = %.4f\n", psi);
	#endif

        // invoke revertToStart on the material
        err += theMat->revertToStart();

        double tangent = theMat->getTangent();
        double stress = theMat->getStress();

	double gamma = -tP * e[10];

	double tau = G * gamma;   // NEED TO ADD G, same for whole section and element

        double value = tangent * A;
        double vas1 = y * value;
        double vas2 = z * value;
        double vas1as2 = vas1 * z;

	kData[0]  += value;
	kData[1]  += -1*vas1;
	kData[2]  += vas2;
	kData[3]  += (y*y+z*z)*value;
	kData[4]  += -omig*value;
	kData[5]  += -psi*value;
	kData[6]  = kData[1];
	kData[7]  += y*vas1;
	kData[8]  += -1*vas1as2;
	kData[9]  += -1*vas1*(y*y+z*z);
	kData[10] += vas1*omig;
	kData[11] += vas1*psi;
	kData[12] = kData[2];
	kData[13] = kData[8];
	kData[14] += z*vas2; 
	kData[15] += vas2*(y*y+z*z);
	kData[16] += -vas2*omig;
	kData[17] += -vas2*psi;
	kData[18] = kData[3];
	kData[19] = kData[9];
	kData[20] = kData[15];
	kData[21] += (y*y+z*z)*(y*y+z*z)*value;
	kData[22] += -omig*value*(y*y+z*z);
	kData[23] += -psi*value*(y*y+z*z);
	kData[24] = kData[4];
	kData[25] = kData[10];
	kData[26] = kData[16];
	kData[27] = kData[22];
	kData[28] += omig*omig*value;
	kData[29] += omig*psi*value;
	kData[30] = kData[5];
	kData[31] = kData[11];
	kData[32] = kData[17];
	kData[33] = kData[23];
	kData[34] = kData[29];
	kData[35] += psi*psi*value+ (G * A * tP * tP)/3;

        double fs0 = stress * A;  // need tau for sData[5]

        sData[0] += fs0;
        sData[1] += -1.0 * fs0 * y;
        sData[2] += 1.0 * fs0 * z;
        sData[3] += fs0 * (y * y + z * z);
        sData[4] += -fs0 * omig;
	sData[5] += -fs0 * psi - (tau * A * tP)/3;
    }


    return err;
}

int TaperedFiberSectionSmoothing3d::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    // create an id to send objects tag and numFibers, 
    //     size 3 so no conflict with matData below if just 1 fiber
    static ID data(3);
    data(0) = this->getTag();
    data(1) = numFibers;
    int dbTag = this->getDbTag();
    res += theChannel.sendID(dbTag, commitTag, data);
    if (res < 0) {
        opserr <<
            "FiberSection2d::sendSelf - failed to send ID data\n";
        return res;
    }

    if (numFibers != 0) {

        // create an id containingg classTag and dbTag for each material & send it
        ID materialData(2 * numFibers);
        for (int i = 0; i < numFibers; i++) {
            UniaxialMaterial *theMat = theMaterials[i];
            materialData(2 * i) = theMat->getClassTag();
            int matDbTag = theMat->getDbTag();
            if (matDbTag == 0) {
                matDbTag = theChannel.getDbTag();
                if (matDbTag != 0)
                    theMat->setDbTag(matDbTag);
            }
            materialData(2 * i + 1) = matDbTag;
        }

        res += theChannel.sendID(dbTag, commitTag, materialData);
        if (res < 0) {
            opserr <<
                "FiberSection2d::sendSelf - failed to send material data\n";
            return res;
        }

        // send the fiber data, i.e. area and loc
        Vector fiberData(matData, 5 * numFibers);
        res += theChannel.sendVector(dbTag, commitTag, fiberData);
        if (res < 0) {
            opserr <<
                "FiberSection2d::sendSelf - failed to send material data\n";
            return res;
        }

        // now invoke send(0 on all the materials
        for (int j = 0; j < numFibers; j++)
            theMaterials[j]->sendSelf(commitTag, theChannel);
    }

    return res;
}

int TaperedFiberSectionSmoothing3d::recvSelf(int commitTag, Channel & theChannel,
                             FEM_ObjectBroker & theBroker)
{
    int res = 0;

    static ID data(3);

    int dbTag = this->getDbTag();
    res += theChannel.recvID(dbTag, commitTag, data);

    if (res < 0) {
        opserr <<
            "FiberSection2d::sendSelf - failed to recv ID data\n";
        return res;
    }

    this->setTag(data(0));

    // recv data about materials objects, classTag and dbTag
    if (data(1) != 0) {
        ID materialData(2 * data(1));
        res += theChannel.recvID(dbTag, commitTag, materialData);
        if (res < 0) {
            opserr <<
                "FiberSection2d::sendSelf - failed to send material data\n";
            return res;
        }

        // if current arrays not of correct size, release old and resize
        if (theMaterials == 0 || numFibers != data(1)) {
            // delete old stuff if outa date
            if (theMaterials != 0) {
                for (int i = 0; i < numFibers; i++)
                    delete theMaterials[i];
                delete[]theMaterials;
                if (matData != 0)
                    delete[]matData;
                matData = 0;
                theMaterials = 0;
            }

            // create memory to hold material pointers and fiber data
            numFibers = data(1);
            if (numFibers != 0) {

                theMaterials = new UniaxialMaterial *[numFibers];

                if (theMaterials == 0) {
                    opserr <<
                        "FiberSection2d::recvSelf -- failed to allocate Material pointers\n";
                    exit(-1);
                }

                for (int j = 0; j < numFibers; j++)
                    theMaterials[j] = 0;

                matData = new double[numFibers * 5];

                if (matData == 0) {
                    opserr <<
                        "FiberSection2d::recvSelf  -- failed to allocate double array for material data\n";
                    exit(-1);
                }
            }
        }

        Vector fiberData(matData, 5 * numFibers);
        res += theChannel.recvVector(dbTag, commitTag, fiberData);
        if (res < 0) {
            opserr <<
                "FiberSection2d::sendSelf - failed to send material data\n";
            return res;
        }

        int i;
        for (i = 0; i < numFibers; i++) {
            int classTag = materialData(2 * i);
            int dbTag = materialData(2 * i + 1);

            // if material pointed to is blank or not of corrcet type, 
            // release old and create a new one
            if (theMaterials[i] == 0)
                theMaterials[i] =
                    theBroker.getNewUniaxialMaterial(classTag);
            else if (theMaterials[i]->getClassTag() != classTag) {
                delete theMaterials[i];
                theMaterials[i] =
                    theBroker.getNewUniaxialMaterial(classTag);
            }

            if (theMaterials[i] == 0) {
                opserr <<
                    "FiberSection2d::recvSelf -- failed to allocate double array for material data\n";
                exit(-1);
            }

            theMaterials[i]->setDbTag(dbTag);
            res +=
                theMaterials[i]->recvSelf(commitTag, theChannel,
                                          theBroker);
        }

        double Qz = 0.0;
        double Qy = 0.0;
        double A = 0.0;
        double yLoc, zLoc, Area;

        // Recompute centroid
        for (i = 0; i < numFibers; i++) {
            yLoc   = matData[5 * i];
            zLoc   = matData[5 * i + 1];
            Area   = matData[5 * i + 2];

            A  += Area;
            Qz += yLoc * Area;
            Qy += zLoc * Area;
        }

        yBar = Qz / A;
        zBar = Qy / A;
    }

    return res;
}

void TaperedFiberSectionSmoothing3d::Print(OPS_Stream & s, int flag)
{
    if (flag == 2) {
        for (int i = 0; i < numFibers; i++) {
            s << -matData[5 * i] << " " << matData[5 * i +
                                                   1] << " " <<
                matData[5 * i + 2] << " ";
            s << theMaterials[i]->
                getStress() << " " << theMaterials[i]->
                getStrain() << endln;
        }
    }
    else {
        s << "\nTaperedFiberSectionSmoothing3d, tag: " << this->getTag() << endln;
        s << "\tSection code: " << code;
        s << "\tNumber of Fibers: " << numFibers << endln;
        s << "\tCentroid: (" << yBar << ", " << zBar << ')' << endln;

        if (flag == 1) {
            for (int i = 0; i < numFibers; i++) {
                s << "\nLocation (y, z) = (" << matData[5 *
                                                         i] << ", " <<
                    matData[5 * i + 1] << ")";
                s << "\nArea = " << matData[5 * i + 2] << endln;
                theMaterials[i]->Print(s, flag);
            }
        }
    }
}

Response *TaperedFiberSectionSmoothing3d::setResponse(const char **argv, int argc,
                                      OPS_Stream & output)
{

    const ID & type = this->getType();
    int typeSize = this->getOrder();

    Response *theResponse = 0;

    output.tag("SectionOutput");
    output.attr("secType", this->getClassType());
    output.attr("secTag", this->getTag());

    // deformations
    if (strcmp(argv[0], "deformations") == 0
        || strcmp(argv[0], "deformation") == 0) {
        for (int i = 0; i < typeSize; i++) {
            int code = type(i);
            switch (code) {
            case SECTION_RESPONSE_MZ:
                output.tag("ResponseType", "kappaZ");
                break;
            case SECTION_RESPONSE_P:
                output.tag("ResponseType", "eps");
                break;
            case SECTION_RESPONSE_VY:
                output.tag("ResponseType", "gammaY");
                break;
            case SECTION_RESPONSE_MY:
                output.tag("ResponseType", "kappaY");
                break;
            case SECTION_RESPONSE_VZ:
                output.tag("ResponseType", "gammaZ");
                break;
            case SECTION_RESPONSE_T:
                output.tag("ResponseType", "theta");
                break;
            default:
                output.tag("ResponseType", "Unknown");
            }
        }
        theResponse =
            new MaterialResponse(this, 1,
                                 this->getSectionDeformation());

        // forces
    }
    else if (strcmp(argv[0], "forces") == 0
             || strcmp(argv[0], "force") == 0) {
        for (int i = 0; i < typeSize; i++) {
            int code = type(i);
            switch (code) {
            case SECTION_RESPONSE_MZ:
                output.tag("ResponseType", "Mz");
                break;
            case SECTION_RESPONSE_P:
                output.tag("ResponseType", "P");
                break;
            case SECTION_RESPONSE_VY:
                output.tag("ResponseType", "Vy");
                break;
            case SECTION_RESPONSE_MY:
                output.tag("ResponseType", "My");
                break;
            case SECTION_RESPONSE_VZ:
                output.tag("ResponseType", "Vz");
                break;
            case SECTION_RESPONSE_T:
                output.tag("ResponseType", "T");
                break;
            default:
                output.tag("ResponseType", "Unknown");
            }
        }
        theResponse =
            new MaterialResponse(this, 2, this->getStressResultant());

        // force and deformation
    }
    else if (strcmp(argv[0], "forceAndDeformation") == 0) {
        for (int i = 0; i < typeSize; i++) {
            int code = type(i);
            switch (code) {
            case SECTION_RESPONSE_MZ:
                output.tag("ResponseType", "kappaZ");
                break;
            case SECTION_RESPONSE_P:
                output.tag("ResponseType", "eps");
                break;
            case SECTION_RESPONSE_VY:
                output.tag("ResponseType", "gammaY");
                break;
            case SECTION_RESPONSE_MY:
                output.tag("ResponseType", "kappaY");
                break;
            case SECTION_RESPONSE_VZ:
                output.tag("ResponseType", "gammaZ");
                break;
            case SECTION_RESPONSE_T:
                output.tag("ResponseType", "theta");
                break;
            default:
                output.tag("ResponseType", "Unknown");
            }
        }
        for (int j = 0; j < typeSize; j++) {
            int code = type(j);
            switch (code) {
            case SECTION_RESPONSE_MZ:
                output.tag("ResponseType", "Mz");
                break;
            case SECTION_RESPONSE_P:
                output.tag("ResponseType", "P");
                break;
            case SECTION_RESPONSE_VY:
                output.tag("ResponseType", "Vy");
                break;
            case SECTION_RESPONSE_MY:
                output.tag("ResponseType", "My");
                break;
            case SECTION_RESPONSE_VZ:
                output.tag("ResponseType", "Vz");
                break;
            case SECTION_RESPONSE_T:
                output.tag("ResponseType", "T");
                break;
            default:
                output.tag("ResponseType", "Unknown");
            }
        }

        theResponse =
            new MaterialResponse(this, 4,
                                 Vector(2 * this->getOrder()));

    }

    else {
        if (argc > 2 || strcmp(argv[0], "TaperedFiber") == 0) {

            int key = numFibers;
            int passarg = 2;


            if (argc <= 3) {    // fiber number was input directly

                key = atoi(argv[1]);

            }
            else if (argc > 4) {        // find fiber closest to coord. with mat tag
                int matTag = atoi(argv[3]);
                double yCoord = atof(argv[1]);
                double zCoord = atof(argv[2]);
                double closestDist;
                double ySearch, zSearch, dy, dz;
                double distance;
                int j;

                // Find first fiber with specified material tag
                for (j = 0; j < numFibers; j++) {
                    if (matTag == theMaterials[j]->getTag()) {
                        ySearch = matData[5 * j];
                        zSearch = matData[5 * j + 1];
                        dy = ySearch - yCoord;
                        dz = zSearch - zCoord;
                        closestDist = sqrt(dy * dy + dz * dz);
                        key = j;
                        break;
                    }
                }

                // Search the remaining fibers
                for (; j < numFibers; j++) {
                    if (matTag == theMaterials[j]->getTag()) {
                        ySearch = matData[5 * j];
                        zSearch = matData[5 * j + 1];
                        dy = ySearch - yCoord;
                        dz = zSearch - zCoord;
                        distance = sqrt(dy * dy + dz * dz);
                        if (distance < closestDist) {
                            closestDist = distance;
                            key = j;
                        }
                    }
                }
                passarg = 4;
            }

            else {              // fiber near-to coordinate specified
                double yCoord = atof(argv[1]);
                double zCoord = atof(argv[2]);
                double closestDist;
                double ySearch, zSearch, dy, dz;
                double distance;
                ySearch = matData[0];
                zSearch = matData[1];
                dy = ySearch - yCoord;
                dz = zSearch - zCoord;
                closestDist = sqrt(dy * dy + dz * dz);
                key = 0;
                for (int j = 1; j < numFibers; j++) {
                    ySearch = matData[5 * j];
                    zSearch = matData[5 * j + 1];
                    dy = ySearch - yCoord;
                    dz = zSearch - zCoord;
                    distance = sqrt(dy * dy + dz * dz);
                    if (distance < closestDist) {
                        closestDist = distance;
                        key = j;
                    }
                }
                passarg = 3;
            }

            if (key < numFibers && key >= 0) {
                output.tag("FiberOutput");
                output.attr("yLoc", -matData[5 * key]);
                output.attr("zLoc", matData[5 * key + 1]);
                output.attr("area", matData[5 * key + 2]);

                theResponse =
                    theMaterials[key]->setResponse(&argv[passarg],
                                                   argc - passarg,
                                                   output);

                output.endTag();
            }
        }
    }

    output.endTag();
    return theResponse;
}


int TaperedFiberSectionSmoothing3d::getResponse(int responseID,
                                Information & sectInfo)
{
    // Just call the base class method ... don't need to define
    // this function, but keeping it here just for clarity
    return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int TaperedFiberSectionSmoothing3d::setParameter(const char **argv, int argc,
                                 Parameter & param)
{
    if (argc < 3)
        return -1;


    int result = 0;

    // A material parameter
    if (strstr(argv[0], "material") != 0) {

        // Get the tag of the material
        int paramMatTag = atoi(argv[1]);

        // Loop over fibers to find the right material(s)
        int ok = 0;
        for (int i = 0; i < numFibers; i++)
            if (paramMatTag == theMaterials[i]->getTag()) {
                ok = theMaterials[i]->setParameter(&argv[2], argc - 2,
                                                   param);
                if (ok != -1)
                    result = ok;
            }

        return result;
    }

    int ok = 0;

    // loop over every material
    for (int i = 0; i < numFibers; i++) {
        ok = theMaterials[i]->setParameter(argv, argc, param);
        if (ok != -1)
            result = ok;
    }

    return result;
}

const Vector &
    TaperedFiberSectionSmoothing3d::getSectionDeformationSensitivity(int gradIndex)
{
    static Vector dummy(3);
    dummy.Zero();
    if (SHVs != 0) {
        dummy(0) = (*SHVs) (0, gradIndex);
        dummy(1) = (*SHVs) (1, gradIndex);
        dummy(2) = (*SHVs) (2, gradIndex);
    }
    return dummy;
}


const Vector &
    TaperedFiberSectionSmoothing3d::getStressResultantSensitivity(int gradIndex,
                                                  bool conditional)
{

    static Vector ds(3);

    ds.Zero();

    double stressGradient;
    int loc = 0;


    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat = theMaterials[i];
        double y = matData[loc++] - yBar;
        double z = matData[loc++] - zBar;
        double A = matData[loc++];
        stressGradient =
            theMaterials[i]->getStressSensitivity(gradIndex,
                                                  conditional);
        stressGradient *= A;
        ds(0) += stressGradient;
        ds(1) += stressGradient * y;
        ds(2) += stressGradient * z;

    }

    return ds;
}

const Matrix &
    TaperedFiberSectionSmoothing3d::getSectionTangentSensitivity(int gradIndex)
{
    static Matrix something(2, 2);

    something.Zero();

    return something;
}

int TaperedFiberSectionSmoothing3d::commitSensitivity(const Vector & defSens,
                                      int gradIndex, int numGrads)
{

    // here add SHVs to store the strain sensitivity.

    if (SHVs == 0) {
        SHVs = new Matrix(3, numGrads);
    }

    (*SHVs) (0, gradIndex) = defSens(0);
    (*SHVs) (1, gradIndex) = defSens(1);
    (*SHVs) (2, gradIndex) = defSens(2);

    int loc = 0;

    double d0 = defSens(0);
    double d1 = defSens(1);
    double d2 = defSens(2);

    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat = theMaterials[i];
        double y = matData[loc++] - yBar;
        double z = matData[loc++] - zBar;
        loc++;                  // skip A data.

        double strainSens = d0 + y * d1 + z * d2;



        theMat->commitSensitivity(strainSens, gradIndex, numGrads);
    }

    return 0;
}

// AddingSensitivity:END ///////////////////////////////////
