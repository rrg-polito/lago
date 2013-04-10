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
#include <cstdlib>
#include <sys/stat.h>
using namespace std;

// Main code
int main(int argc, char ** argv) {
	if (argc == 1) {
		cerr
				<< "Please use the following syntax:  \"./graphSLAM_sc nameOfInputFile.graph (add d for debug if desired) \"."
				<< endl;
		return 0;
	}

	bool debug = 0;
	if ((argc == 3) && (*argv[2] == 'd'))
		debug = 1;
	/*****************************************************************************************
	 * retrieve the odometry data;
	 ***************************************************************************************/

	cout << "Retrieving the original data ..." << endl;
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
	cout << "Data retrieved." << endl;

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
	//  construct information matrix of the relative measurements
	//    construct covariance matrix of the relative orientation measurements
	cs * PDeltaInvMatrix_ccf = constructPdeltaInvMatrix(vEdges);
	//    construct the covariance matrix of the relative position measurements
	cs * PDeltaLInverseMatrix_ccf = constructPDeltaLInverseMatrix(vEdges);

	// check the correctness of the first odometric constraints section of the EDGES data;
	/****************************************************************************************/
	unsigned int cntCk = 0;
	for (int idx = 0; idx < (int) (vPoses.size() - 1); idx++) {
		if ((vEdges.at(idx).id1 == idx) && (vEdges.at(idx).id2 == idx + 1)) {
			cntCk++;
		}
	}
	if (cntCk == vPoses.size() - 1) {
		cout << "Check on the format of the input data successfully completed!"
				<< endl;
	} else {
		cout << "number of odometric constraints in the EDGES data: " << cntCk
				<< endl;
		cout << "number of poses: " << vPoses.size() << endl;
		cout << "The first odometric section of the Edges data is NOT correct."
				<< endl;
		return 0;
	}

	/****************************************************************************************/

	/*
	 * to verify the correctness of the data retrieved;
	 */

	ofstream ofs;
	if (debug == 1) {
		cout << "number of poses: " << vPoses.size() << endl;
		cout << "number of odometric constraints in the EDGES data: " << cntCk
				<< endl;
		cout << "number of loop-closing constraints: " << vEdges.size() - cntCk
				<< endl;
		cout << "number of edges: " << vEdges.size() << endl;

		//	to print the pose vector data to file;

		string fn1 = "initial_poses.txt";
		ofs.open(fn1.c_str());
		ofs.precision(numeric_limits<double>::digits10 + 1);
		for (unsigned int i = 0; i < vPoses.size(); i++) {

			ofs << vPoses.at(i).id << "\t" << vPoses.at(i).x << "\t"
					<< vPoses.at(i).y << "\t" << vPoses.at(i).theta << "\r\n";
		}
		ofs.close();

		//		        to print the edge vector data to file;

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

		string fn3 = "data.graph";
		ofs.open(fn3.c_str());
		ofs.precision(numeric_limits<double>::digits10 + 1);
		for (unsigned int i = 0; i < vPoses.size(); i++) {

			ofs << "VERTEX2" << "\t" << vPoses.at(i).id + 1 << "\t"
					<< vPoses.at(i).x << "\t" << vPoses.at(i).y << "\t"
					<< vPoses.at(i).theta << "\r\n";
		}

		for (unsigned int i = 0; i < vEdges.size(); i++) {
			ofs << "EDGE2" << "\t" << vEdges.at(i).id1 << "\t"
					<< vEdges.at(i).id2 << "\t" << vEdges.at(i). meas1 << "\t"
					<< vEdges.at(i).meas2 << "\t" << vEdges.at(i).meas3 << "\t"
					<< vEdges.at(i).t_ff << "\t" << vEdges.at(i).t_fs << "\t"
					<< vEdges.at(i).t_ss << "\t" << vEdges.at(i).t_rr << "\t"
					<< vEdges.at(i).t_fr << "\t" << vEdges.at(i).t_sr << "\r\n";
		}

		ofs.close();

		// for original positions plot;
		string fn4 = "Positions_original.txt";
		ofs.open(fn4.c_str());
		ofs.precision(numeric_limits<double>::digits10 + 1);
		for (unsigned int i = 0; i < vPoses.size(); i++) {

			ofs << vPoses.at(i).x << "\t" << vPoses.at(i).y << "\r\n";
		}
		ofs.close();

		//	standardize the POSE and the EDGES data;
		if (vPoses.at(0).id == 1) {
			for (unsigned int idx = 0; idx < vPoses.size(); idx++) {
				vPoses.at(idx).id -= 1;
			}
			for (unsigned int i = 0; i < vEdges.size(); i++) {
				vEdges.at(i).id1 -= 1;
				vEdges.at(i).id2 -= 1;
			}
			cout << "original data standardized." << endl;
		}
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

	//  construct the vector of relative orientation measurements  delta
	cs * orientationVector = constructOrientationVector(vEdges);
	cs * orientationVector_ccf = cs_compress(orientationVector);
	cs * v = regularization(vEdges, vPoses);
	cs * v_ccf = cs_compress(v);
	cs * v_r = roundTo2Pi(v_ccf);
	cs * v_r_ccf = cs_compress(v_r);
	cs * regularizedOrientationVector = cs_add(orientationVector_ccf, v_r_ccf,
			1, -1);

	//cs * regularizedOrientationVector = constructRegularizedOrientationVector(edges, m, n);

	if (debug == 1) {
		double * delta_r = fromCompressedColumnVectorToDoubleArrayPointer(
				regularizedOrientationVector);
		writeDoubleArrayToDisk(delta_r, regularizedOrientationVector->m,
				"delta_regularized.txt");
		cout << "orientations of the original pose data has been regularized."
				<< endl;

		cs * Unity1 = cs_spalloc(vEdges.size(), 1, vEdges.size(), 1, 1);
		for (unsigned int i = 0; i < vEdges.size(); i++) {
			cs_entry(Unity1, i, 0, 1);
		}
		cs * Unity1_ccf = cs_compress(Unity1);
		cs * R1 = cs_multiply(PDeltaInvMatrix_ccf, Unity1_ccf);
		double * RArr1 = fromCompressedColumnVectorToDoubleArrayPointer(R1);
		writeDoubleArrayToDisk(RArr1, R1->m, "summationForPDeltaInv.txt");
	}
	/***************************************************************************************
	 * estimate nodes' orientations
	 ***************************************************************************************/

	//    construct reduced Incidence matrix
	cs * reducedIncidenceMatrixTran_ccf = constructReducedIncidenceMatrixTran(
			edges, m, n);
	cs* reducedIncidenceMatrix =
			cs_transpose(reducedIncidenceMatrixTran_ccf, 1);

	// A * Omega_delta
	cs* IncidenceAndPdeltaInvMatrix = cs_multiply(reducedIncidenceMatrix,
			PDeltaInvMatrix_ccf);

	// A * Omega_delta * delta
	cs* IncidencePdeltaInvAndOrientationVector = cs_multiply(
			IncidenceAndPdeltaInvMatrix, regularizedOrientationVector);
	double * doubleIncidencePdeltaInvAndOrientationVector =
			fromCompressedColumnVectorToDoubleArrayPointer(
					IncidencePdeltaInvAndOrientationVector);

	// A * Omega_delta * A'
	cs * PThetaEstimatedInv = cs_multiply(IncidenceAndPdeltaInvMatrix,
			reducedIncidenceMatrixTran_ccf);

	// solve linear system (A * Omega_delta * A') theta = A * Omega_delta * delta
	int order = 1;
	//double tolerance = 0.001;
	//cs_lusol(order, PThetaEstimatedInv,
	//		doubleIncidencePdeltaInvAndOrientationVector, tolerance);
			cs_cholsol(order, PThetaEstimatedInv,
					doubleIncidencePdeltaInvAndOrientationVector);
	double * theta_estimated = doubleIncidencePdeltaInvAndOrientationVector;

	if (debug == 1) {
		writeDoubleArrayToDisk(theta_estimated, vPoses.size() - 1,
				"theta_estimated.txt");
		cout << ".................the estimation of theta done." << endl;
	}

	/***************************************************************************************
	 *    Phase 2: express measurements in global frame
	 ***************************************************************************************/

	//  construct the vector of relative position measurements  Delta
	cs * pV_ccf = constructPositionColumnVector(edges, m);
	double * position_meas = fromCompressedColumnVectorToDoubleArrayPointer(
			pV_ccf);

	//    construct rotation matrix;  R
	cs* rotationMatrix_ccf = constructRotationMatrixWithEstimatedOrientations(
			theta_estimated, edges, m);

	if (debug == 1) {
		writeDoubleArrayToDisk(position_meas, pV_ccf->m, "Delta.txt");
		cout << "vector of relative position measurement done." << endl;
		cout << ".................the estimation of [rho;theta]^T starts."
				<< endl;
		cout << "rotation matrix done." << endl;
	}

	//    construct z vector, with z = [Delta^g  hat{theta}];
	cs * z_position = cs_multiply(rotationMatrix_ccf, pV_ccf); //Delta^g
	cs * z_ccf = constructZVector(z_position, theta_estimated, n);

	if (debug == 1) {
		double* z_d = fromCompressedColumnVectorToDoubleArrayPointer(z_ccf);
		writeDoubleArrayToDisk(z_d, z_ccf->m, "z.txt");
		cout << "[R_theta*Delta^l;theta_hat] has been constructed." << endl;
		cout << "P_Delta^l constructed." << endl;

		cs * Unity2 = cs_spalloc(2 * vEdges.size(), 1, 2 * vEdges.size(), 1, 1);
		for (unsigned int i = 0; i < 2 * vEdges.size(); i++) {
			cs_entry(Unity2, i, 0, 1);
		}
		cs * Unity2_ccf = cs_compress(Unity2);
		cs * R2 = cs_multiply(PDeltaLInverseMatrix_ccf, Unity2_ccf);
		double * RArr2 = fromCompressedColumnVectorToDoubleArrayPointer(R2);
		writeDoubleArrayToDisk(RArr2, R2->m, "summationForPDeltaLInverse.txt");
	}

	//    transform the covariance in the global frame: PDeltaGInverseMatrix = R * PDeltaL * R'
	cs * rotationAndPDeltaLInverseMatrix = cs_multiply(rotationMatrix_ccf,
			PDeltaLInverseMatrix_ccf);
	cs * RotationMatrixTran = cs_transpose(rotationMatrix_ccf, 1);
	cs* PDeltaGInverseMatrix = cs_multiply(rotationAndPDeltaLInverseMatrix,
			RotationMatrixTran);

	//    construct Jacobian Matrix;
	cs * Jacobian_ccf = constructNegJacobianMatrix(edges, theta_estimated,
			position_meas, m, n);
	cs * JacobianTran = cs_transpose(Jacobian_ccf, 1);

	//construct PZInverse Matrix;
	// note that the minus is included in the definition of the Jacobian
	cs * UR = cs_multiply(PDeltaGInverseMatrix, Jacobian_ccf); // - Omega_Delta^g * J

	cs * U_t = combineTwoSparseMatrixLeftAndRight(PDeltaGInverseMatrix, UR);// [Omega_Delta^g  - Omega_Delta^g * J]
	cs * U = cs_compress(U_t);

	// note that the minus is included in the definition of the Jacobian
	cs * LL = cs_multiply(JacobianTran, PDeltaGInverseMatrix);

	cs * JacobianTranPDeltaGInverseJacobian = cs_multiply(LL, Jacobian_ccf);

	cs * LR = cs_add(PThetaEstimatedInv, JacobianTranPDeltaGInverseJacobian, 1,
			1);

	cs * L_t = combineTwoSparseMatrixLeftAndRight(LL, LR);
	cs * L = cs_compress(L_t);

	cs * PZInv_t = combineTwoSparseMatrixUpAndLow(U, L);
	cs * PZInv = cs_compress(PZInv_t);

	//    contruct Matrix B' = [A_2t  zeros;  zeros  eye(n)]
	cs * BTran_ccf = constructBT(edges, m, n);
	cs * B_ccf = cs_transpose(BTran_ccf, 1);

	//   solve the linear equation: EL p = ER
	cs * BPZInv = cs_multiply(B_ccf, PZInv);
	cs * EL = cs_multiply(BPZInv, BTran_ccf);
	cs * ER = cs_multiply(BPZInv, z_ccf);
	double * dER = fromCompressedColumnVectorToDoubleArrayPointer(ER);

	//cs_lusol(order, EL, dER, tolerance);
cs_cholsol(order, EL,dER);
	double * poseEstimated = dER;

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
	//*********************************************************************************/

	if (debug == 1) {
		int row = 3 * (vPoses.size() - 1);
		writeDoubleArrayToDisk(poseEstimated, row, "pose_estimated.txt");
		cout
				<< ".................the estimation of [rho;theta]^T has finished."
				<< endl;
	}

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

