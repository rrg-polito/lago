/*
 *	Created on: Mar 27, 2012
 *  Authors: Luca Carlone, Solomon Jingchun Yin, Stefano Rosa
 *			Robotics Research Group, Politecnico di Torino, Torino, Italy
 *
 * This software is licenced under the Common Creative License,
 * Attribution-NonCommercial-ShareAlike 3.0
 *
 * You are free:
 *   - to Share - to copy, distribute and transmit the work
 *   - to Remix - to adapt the work
 *
 * Under the following conditions:
 *
 *   - Attribution. You must attribute the work in the manner specified
 *     by the author or licensor (but not in any way that suggests that
 *     they endorse you or your use of the work).
 *  
 *   - Noncommercial. You may not use this work for commercial purposes.
 *  
 *   - Share Alike. If you alter, transform, or build upon this work,
 *     you may distribute the resulting work only under the same or
 *     similar license to this one.
 *
 * Any of the above conditions can be waived if you get permission
 * from the copyright holder.  Nothing in this license impairs or
 * restricts the author's moral rights.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  
 **********************************************************************/

extern "C" {
#include "cs.h"
}
#include "utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
#include <limits>
#include <sys/stat.h>
#include <cstdlib>
using namespace std;

// Main code
int main(int argc, char ** argv) {
	if (argc == 1) {
		cerr
				<< "Please use the following syntax:  \"./graphSLAM nameOfInputFile.graph (add d for debug is desired\"."
				<< endl;
		return 0;
	}

	bool debug = 0;
	if ((argc == 3) && (*argv[2] == 'd'))
		debug = 1;
	/*****************************************************************************************
	 * retrieve the odometry data;
	 ***************************************************************************************/
	cout << "Retrieving the edges data ..." << endl;
	ifstream ifs;
	string filename = argv[1];
	ifs.open(filename.c_str());
	if (!ifs.is_open()) {
		cout << "Can't open file " << filename << endl;
		return -1;
	}
	ifs.precision(numeric_limits<double>::digits10 + 1);
	vector < pose > vPoses;
	int id;
	double a, b, c;
	pose tmpPose;
	int id1, id2;
	double meas1, meas2, meas3;
	double t_ff, t_fs, t_ss, t_rr, t_fr, t_sr;
	edge tmpEdge;
	vector < edge > vEdges;
	string tmpStr;

	string strTmp;
	while (getline(ifs, tmpStr)) {
		stringstream ss;
		ss << tmpStr;
		ss >> strTmp;
		if (strTmp == "VERTEX2") {
			ss >> id >> a >> b >> c;
			tmpPose.id = id;
			tmpPose.x = a;
			tmpPose.y = b;
			tmpPose.theta = c;
			vPoses.push_back(tmpPose);
		} else if (strTmp == "EDGE2") {
			ss >> id1 >> id2 >> meas1 >> meas2 >> meas3 >> t_ff >> t_fs >> t_ss
					>> t_rr >> t_fr >> t_sr;
			tmpEdge.id1 = id1;
			tmpEdge.id2 = id2;
			tmpEdge.meas1 = meas1;
			tmpEdge.meas2 = meas2;
			tmpEdge.meas3 = meas3;
			tmpEdge.t_ff = t_ff;
			tmpEdge.t_fs = t_fs;
			tmpEdge.t_ss = t_ss;
			tmpEdge.t_rr = t_rr;
			tmpEdge.t_fr = t_fr;
			tmpEdge.t_sr = t_sr;
			vEdges.push_back(tmpEdge);
		}
	}
	cout << "Data retrieved!" << endl;

	/****************************************************************************************/

	int m = vEdges.size(); // number of edges (measurements)
	int np1 = vPoses.size(); // n+1 : number of nodes
	int n = np1 - 1; // number of nodes excluding the anchor node

	double **edges = new double *[m]; // each row is "id1 id2 dx dy dth"
	double *Delta_l = new double[2 * m]; //relative position measurements
	double *delta = new double[m]; // relative orientation measurements (non regularized)
	for (int i = 0; i < m; i++) {
		edges[i] = new double[5];
		edges[i][0] = vEdges.at(i).id1;
		edges[i][1] = vEdges.at(i).id2;
		edges[i][2] = vEdges.at(i).meas1;
		edges[i][3] = vEdges.at(i).meas2;
		edges[i][4] = vEdges.at(i).meas3;
		Delta_l[2 * i] = vEdges.at(i).meas1;
		Delta_l[2 * i + 1] = vEdges.at(i).meas2;
		delta[i] = vEdges.at(i).meas3;
	}
	//  construct information matrix of the relative measurements: Omega = diag(Omega_Delta, Omega_delta)
	cs * Omega_ccf = constructOmega(vEdges);

	// check the correctness of the first odometric constraints section of the EDGES data;
	/****************************************************************************************/
	int cntCk = 0;
	for (int idx = 0; idx < n; idx++) {
		if ((edges[idx][1] - edges[idx][0]) == 1 && edges[idx][0] == idx)
			cntCk++;
	}
	if (cntCk == n)
		cout << "Check on the format of input data successfully completed!"
				<< endl;
	else {
		cout << "Invalid format for the input file:" << endl;
		cout << "- the first node needs to have id 0" << endl;
		cout
				<< "- the first n constraints in EDGES needs to be the odometric contraints"
				<< endl;
		return 0;
	}
	/****************************************************************************************/

	/*
	 * to verify the correctness of the data retrieved (only in debug mode);
	 */
	/****************************************************************************************/
	ofstream ofs;
	if (debug == 1) {
		//	to print the pose vector data to file;
		cout << "printing initial poses" << endl;
		string fn1 = "initial_poses.txt";
		ofs.open(fn1.c_str());
		ofs.precision(numeric_limits<double>::digits10 + 1);
		for (unsigned int i = 0; i < vPoses.size(); i++) {

			ofs << vPoses.at(i).id << "\t" << vPoses.at(i).x << "\t"
					<< vPoses.at(i).y << "\t" << vPoses.at(i).theta << "\r\n";
		}
		ofs.close();

		//		        to print the edge vector data to file;
		cout << "printing initial edges" << endl;
		string fn2 = "initial_edges.txt";
		ofs.open(fn2.c_str());
		ofs.precision(numeric_limits<double>::digits10 + 1);
		for (unsigned int i = 0; i < vEdges.size(); i++) {

			ofs << vEdges.at(i).id1 << "\t" << vEdges.at(i).id2 << "\t"
					<< vEdges.at(i). meas1 << "\t" << vEdges.at(i).meas2
					<< "\t" << vEdges.at(i).meas3 << "\t" << vEdges.at(i).t_ff
					<< "\t" << vEdges.at(i).t_fs << "\t" << vEdges.at(i).t_ss
					<< "\t" << vEdges.at(i).t_rr << "\t" << vEdges.at(i).t_fr
					<< "\t" << vEdges.at(i).t_sr << "\r\n";
		}
		ofs.close();
		cout << "printing complete" << endl;
	}
	/****************************************************************************************/

	/***************************************************************************************
	 *    Phase 1: estimate nodes' orientations
	 ***************************************************************************************/

	/******************************************************************************************
	 * start the timer;
	 ******************************************************************************************/
	struct timeval start;
	gettimeofday(&start, NULL);
	double startseconds;
	startseconds = calculateSeconds(start);
	//========================================================================================

	/******************************************************************************************
	 * Regularize the orientation column vector from the original sensor data;
	 *******************************************************************************************/
	//    construct the measurement vector MeasVector = [(Delta^l)'   delta']'
	cs * MeasVector_ccf = constructMeasVector(edges, m, n);

	if (debug == 1) {
		cs * Unity = cs_spalloc(3 * vEdges.size(), 1, 3 * vEdges.size(), 1, 1);
		for (unsigned int i = 0; i < 3 * vEdges.size(); i++) {
			cs_entry(Unity, i, 0, 1);
		}
		cs * Unity_ccf = cs_compress(Unity);
		cs * R = cs_multiply(Omega_ccf, Unity_ccf);
		double * RArr = fromCompressedColumnVectorToDoubleArrayPointer(R);
		writeDoubleArrayToDisk(RArr, R->m, "summationForEachRowOfOmega.txt");
	}

	//  construct matrix: Z = diag(I_2m, A') \in 3m X (2m+n)
	cs * Z_ccf = constructZ(edges, m, n);
	cs * ZTran_ccf = cs_transpose(Z_ccf, 1);

	//    construct matrix: OmegaZ = Z' *  Omega *  Z;
	cs* ZTranOmega_ccf = cs_multiply(ZTran_ccf, Omega_ccf);

	cs* OmegaZ_ccf = cs_multiply(ZTranOmega_ccf, Z_ccf);

	//cs_free(ZTran_ccf);   //LC: insert if needed
	//cs_free(Z_ccf);	//LC: insert if needed

	//    construct vector: bZ = Z' *  Omega *  MeasVector;
	cs* bZ_ccf = cs_multiplyMatVect(ZTranOmega_ccf, MeasVector_ccf);

	double * bZ = fromCompressedColumnVectorToDoubleArrayPointer(bZ_ccf);

	// solve linear system (OmegaZ) z = bZ 
	int order = 1;
	//double tolerance = 0.001;

	//cs_lusol(order, OmegaZ_ccf, bZ, tolerance);
	cs_cholsol(order, OmegaZ_ccf, bZ);
	double * z_hat = bZ; // solution z_hat = [Delta^l   theta_hat]
	cs * z_hat_ccf = fromDoubleArrayPointerToCompressedColumnVector(z_hat,
			2 * m + n);

	if (debug == 1) {
		cout << "phase 1 completed" << endl;
		double * z_hat_deb = fromCompressedColumnVectorToDoubleArrayPointer(
				z_hat_ccf);
		writeDoubleArrayToDisk(z_hat_deb, z_hat_ccf->m, "z_hat.txt");

		cs* ones3 = constructOneMatrix(2 * m + n, 1);
		cs* OmegaZ_vect = cs_multiply(OmegaZ_ccf, ones3);
		double * OmegaZ_vect_deb =
				fromCompressedColumnVectorToDoubleArrayPointer(OmegaZ_vect);
		writeDoubleArrayToDisk(OmegaZ_vect_deb, OmegaZ_vect->m,
				"OmegaZ_vect.txt");
	}

	/***************************************************************************************
	 *    Phase 2: express measurements in global frame
	 ***************************************************************************************/

	///////////////////// CONSTRUCT SUPPORT MATRICES  /////////////////////
	//    construct rotation matrix;  R

	//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
	cs* R_ccf = constructRotationMatrix(z_hat, edges, m);
	cs * RTran_ccf = cs_transpose(R_ccf, 1);

	//    construct transformation matrix: T =  [R    0;  0   I_n]
	cs * T_ccf = constructT(R_ccf, m, n);
	cs* y_ccf = cs_multiply(T_ccf, z_hat_ccf);
	/////////////////////  /////////////////////  /////////////////////  

	// cs * JTinv_ccf = constructTinv(z_hat, edges, m, n);

	//    construct (negative) Jacobian Matrix;
	cs * nJ_ccf = constructNegJacobianMatrix(z_hat, edges, m, n);

	//    construct Identity matrix of size n X 2m
	cs * Zero_ccf = constructZeroMatrix(n, 2 * m); // (row, col)

	//    construct Identity matrix of size n
	cs * I_n_ccf = constructIdentity(n);

	//  construct inverse of the Jacobian of the transformation matrix: JTinv =  [R'    -R' *J;  0   I_n]
	cs* nRTranJ = cs_multiply(RTran_ccf, nJ_ccf);
	cs * U_JTinv = combineTwoSparseMatrixLeftAndRight(RTran_ccf, nRTranJ);
	cs * D_JTinv = combineTwoSparseMatrixLeftAndRight(Zero_ccf, I_n_ccf);
	cs * JTinv_ccf = combineTwoSparseMatrixUpAndLow(U_JTinv, D_JTinv);

	cs * JTinvTran_ccf = cs_transpose(JTinv_ccf, 1);
	//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

	//  obtain matrix: OmegaY =  Tinv' * OmegaZ * Tinv
	cs* JTinvTranOmega = cs_multiply(JTinvTran_ccf, OmegaZ_ccf);
	cs* OmegaY_ccf = cs_multiply(JTinvTranOmega, JTinv_ccf);

	if (debug == 1) {
		cout << "phase 2 completed" << endl;
		double * y_ccf_deb = fromCompressedColumnVectorToDoubleArrayPointer(
				y_ccf);
		writeDoubleArrayToDisk(y_ccf_deb, y_ccf->m, "y_ccf.txt");
	}

	/***************************************************************************************
	 *    Phase 3: solve the second linear system to obtain the configuration estimate
	 ***************************************************************************************/

	//    contruct Matrix B = [A_2t  zeros;  zeros  eye(n)]
	cs * B_ccf = constructB(edges, m, n);
	cs * BTran_ccf = cs_transpose(B_ccf, 1);

	// construct vector bX = B' * OmegaY * y
	cs * BTranOmegaY = cs_multiply(BTran_ccf, OmegaY_ccf);
	cs * bX_cff = cs_multiply(BTranOmegaY, y_ccf);

	if (debug == 1) {
		cout << "printing data for phase 3" << endl;
		double * bX_cff_deb = fromCompressedColumnVectorToDoubleArrayPointer(
				bX_cff);
		writeDoubleArrayToDisk(bX_cff_deb, bX_cff->m, "bX_cff.txt");
	}

	double * bX = fromCompressedColumnVectorToDoubleArrayPointer(bX_cff);

	// construct matrix omegaX = B' * OmegaY * B
	cs * OmegaX_cff = cs_multiply(BTranOmegaY, B_ccf);

	//   solve the linear equation: OmegaX x = bX
	//	cs_lusol(order, OmegaX_cff, bX, tolerance);
	cs_cholsol(order, OmegaX_cff, bX);
	double * poseEstimated = bX;

	/*********************************************************************************
	 * Terminate the timer;
	 *********************************************************************************/
	struct timeval end;
	gettimeofday(&end, NULL);
	double endSeconds;
	endSeconds = calculateSeconds(end);

	double elapsedSeconds;
	elapsedSeconds = endSeconds - startseconds;
	cout << "Execution time:" << elapsedSeconds << " seconds." << endl;
	//========================================================================================

	if (debug == 1) {
		cout << "phase 3 completed" << endl;
		int row = 3 * (vPoses.size() - 1);
		writeDoubleArrayToDisk(poseEstimated, row, "pose_estimated.txt");
	}
	//***************************************************

    //***************************************************
    // write results on file
    ofs.open("output_poses.txt");
    ofs.precision(numeric_limits<double>::digits10 + 1);
    ofs << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\r\n";
    for (unsigned int idx = 0; idx < vPoses.size() - 1; idx++) {
        ofs << poseEstimated[2 * idx] << "\t" << poseEstimated[2 * idx + 1]
		<< "\t" << standardizeOrientation(
										  poseEstimated[idx + 2 * (vPoses.size() - 1)]) << "\r\n";
    }
    ofs.close();
	
    // write results in TORO format for evaluation of the accuracy and other benchmarking metrics
	
    //    *******************************************
    ofs.open("output_graph.txt");
    ofs.precision(numeric_limits<double>::digits10 + 1);
    ofs << "VERTEX2" << "\t" << 0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0
	<< "\r\n";
    for (unsigned int idx = 0; idx < vPoses.size() - 1; idx++) {
        ofs << "VERTEX2" << "\t" << idx + 1 << "\t" << poseEstimated[2 * idx]
		<< "\t" << poseEstimated[2 * idx + 1] << "\t"
		<< standardizeOrientation(
								  poseEstimated[idx + 2 * (vPoses.size() - 1)]) << "\r\n";
    }
    for (unsigned int i = 0; i < vEdges.size(); i++) {
		
        ofs << "EDGE2" << "\t" << vEdges.at(i).id1 << "\t" << vEdges.at(i).id2
		<< "\t" << vEdges.at(i). meas1 << "\t" << vEdges.at(i).meas2
		<< "\t" << vEdges.at(i).meas3 << "\t" << vEdges.at(i).t_ff
		<< "\t" << vEdges.at(i).t_fs << "\t" << vEdges.at(i).t_ss
		<< "\t" << vEdges.at(i).t_rr << "\t" << vEdges.at(i).t_fr
		<< "\t" << vEdges.at(i).t_sr << "\r\n";
    }
    ofs.close();


    return 0; 
}
