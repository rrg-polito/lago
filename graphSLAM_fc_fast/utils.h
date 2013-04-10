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

#ifndef THREEPHASE_H_
#define THREEPHASE_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <limits>
#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
using namespace std;

using namespace std;

const double PI = 3.14159265358979;

struct edge {
	int id1;
	int id2;
	double meas1;
	double meas2;
	double meas3;
	double t_ff;
	double t_fs;
	double t_ss;
	double t_rr;
	double t_fr;
	double t_sr;
};

struct pose {
	int id;
	double x;
	double y;
	double theta;
};


// this function is useful for computing statistics on the computational effort
double calculateSeconds(struct timeval time) {
	double seconds;
	seconds = time.tv_sec + 1e-6 * time.tv_usec;
	return seconds;
}


double standardizeOrientation(double delta) {
	double tmp = fmod(delta, 2 * PI);
	if (tmp < -PI) {
		tmp += 2 * PI;
	} else if (tmp > PI) {
		tmp -= 2 * PI;
	}
	return tmp;
}


double * fromCompressedColumnVectorToDoubleArrayPointer(cs* T) {
	//    declare a double array with all the elements set to be zeros;
	int size = T->m;

	double * arrPtr = new double[size];

	for (int i = 0; i < size; i++) {
		arrPtr[i] = 0;
	}
	int ptr = 0;
	int indexLow = *(T->p + ptr);
	int indexHigh = *(T->p + ptr + 1) - 1;
	for (int idx = indexLow; idx <= indexHigh; idx++) {
		int row = *(T->i + idx);
		double val = *(T->x + idx);
		arrPtr[row] = val;
	}
	return arrPtr;
}


cs * constructOrientationVector(vector<edge> vEdges) {
	int rowNum = vEdges.size();
	int colNum = 1;
	int nzmax = rowNum * colNum;
	int values = 1;
	int triplet = 1;
	cs *orientationVector = cs_spalloc(rowNum, colNum, nzmax, values, triplet);
	for (unsigned int i = 0; i < vEdges.size(); i++) {
		cs_entry(orientationVector, i, 0, vEdges.at(i).meas3);
	}

	return orientationVector;
}


cs * constructMeasVector(double ** edges, int m, int n){

	int rowNum = 3 * m;
	int colNum = 1;
	int nzmax = rowNum * colNum;
	int values = 1;
	int triplet = 1;
	cs *MeasVector = cs_spalloc(rowNum, colNum, nzmax, values, triplet);
	double tot_orientation;

	for (int i = 0; i < n; i++) {
		// Relative position measurements Delta^l
		cs_entry(MeasVector, 2*i, 0, edges[i][2]); 
		cs_entry(MeasVector, 2*i+1, 0, edges[i][3]); 
		//  these do not need regularization
		cs_entry(MeasVector, 2*m+i, 0, edges[i][4]); 
	}

	for (int i = n; i < m; i++) { // these have to be regularized
		{
		// Relative position measurements Delta^l
		cs_entry(MeasVector, 2*i, 0, edges[i][2]); 
		cs_entry(MeasVector, 2*i+1, 0, edges[i][3]); 
		// Regularization for relative orientation measurements
		tot_orientation = edges[i][4]; // we first add the loop closing constraints
		int id1 = edges[i][0];
		int id2 = edges[i][1];
		//	assign -1s for forward and 1s on the contrary;
		if (id1 < id2) { // we compute the total orientation over cycles
			for (int j = id1; j <= id2 - 1; j++) {
				tot_orientation = tot_orientation - edges[j][4];
			}
		} else {  // we compute the total orientation over cycles
			for (int j = id2; j <= id1 - 1; j++) {
				tot_orientation = tot_orientation + edges[j][4];
			}
		}
		// here tot_orientation contains the sum of the orientations along a cycle
		tot_orientation = 2 * PI * round(tot_orientation / (2 * PI)); // this is the regularization term
		cs_entry(MeasVector, 2*m+i, 0, edges[i][4] - tot_orientation); 
		}
	}
	cs * MeasVector_ccf = cs_compress(MeasVector);
	return MeasVector_ccf;
}


cs * constructCycleBasisMatrix(vector<edge> vEdges, vector<pose> vPoses) {
	int m = vEdges.size();
	int n = vPoses.size() - 1;
	int numLoop = m - n;
	int rowLen = numLoop;
	int colLen = m;

	cs * CycleBasisMatrix = cs_spalloc(rowLen, colLen, rowLen * colLen, 1, 1);
	for (int idx = 0; idx < rowLen; idx++) { // for each loop
		int id1 = vEdges.at(n + idx).id1;
		int id2 = vEdges.at(n + idx).id2;
		cs_entry(CycleBasisMatrix, idx, n + idx, 1);
		//	assign -1s for forward and 1s on the contrary;
		if (id1 < id2) {
			for (int i = id1; i <= id2 - 1; i++) {
				cs_entry(CycleBasisMatrix, idx, i, -1);
			}
		} else {
			for (int i = id2; i <= id1 - 1; i++) {
				cs_entry(CycleBasisMatrix, idx, i, 1);
			}
		}

	}
	return CycleBasisMatrix;
}


cs * regularization(vector<edge> vEdges, vector<pose> vPoses) {
	int m = vEdges.size();
	int n = vPoses.size() - 1;
	int numLoop = m - n;
	int rowLen = numLoop;
	int colLen = 1;
	cs * v = cs_spalloc(rowLen, colLen, rowLen * colLen, 1, 1);

	double * tmp_orientations = new double[m];
	for (int i = 0; i < m; i++) {
		tmp_orientations[i] = vEdges.at(i).meas3;
	}

	for (int idx = 0; idx < rowLen; idx++) { // for each loop
		double tot_orientation = tmp_orientations[n + idx]; // we first add the loop closing constraints
		int id1 = vEdges.at(n + idx).id1;
		int id2 = vEdges.at(n + idx).id2;
		//	assign -1s for forward and 1s on the contrary;
		if (id1 < id2) { // we compute the total orientation over cycles
			for (int i = id1; i <= id2 - 1; i++) {
				tot_orientation = tot_orientation - tmp_orientations[i];
			}
		} else {  // we compute the total orientation over cycles
			for (int i = id2; i <= id1 - 1; i++) {
				tot_orientation = tot_orientation + tmp_orientations[i];
			}
		}
		double val=tot_orientation;
		cs_entry(v, n + idx, 0, val);
	}
	return v;
}


cs * constructV(vector<edge> vEdges, vector<pose> vPoses, cs* CDelta) {
	int m = vEdges.size();
	int n = vPoses.size() - 1;
	cs * v = cs_spalloc(m, 1, m - n, 1, 1);

	int ptr = 0;
	int indexLow = *(CDelta->p + ptr);
	int indexHigh = *(CDelta->p + ptr + 1) - 1;
	for (int idx = indexLow; idx <= indexHigh; idx++) {
		int r = *(CDelta->i + idx);
		double val = *(CDelta->x + idx);
		cs_entry(v, n + r, 0, val);
	}

	return v;
}


cs * roundTo2Pi(cs* v) {
	int m = v->m;
	cs * v_r = cs_spalloc(m, 1, m, 1, 1);

	for (int idx = 0; idx < v->nzmax; idx++) {
		double tmp = *(v->x + idx);
		tmp = 2 * PI * round(tmp / (2 * PI));
		int i = *(v->i + idx);
		cs_entry(v_r, i, 0, tmp);
	}
	return v_r;
}


cs * constructZ(double ** edges, int m, int n) {

	int rowLen = 3*m;
	int colLen = 2* m + n;
	cs * Z = cs_spalloc(rowLen, colLen, 4*m, 1, 1);

	for (int i = 0; i < m; i++) { // for each edge
		cs_entry(Z, 2*i, 2*i, 1);
		cs_entry(Z, 2*i+1, 2*i+1, 1);
		
		int id1 = edges[i][0]; // we pick the tail of the i-th edge
		int id2 = edges[i][1]; // we pick the head of the i-th edge
		
		if (id1 > 0) {
			cs_entry(Z, 2*m + i, 2*m + id1 - 1, -1);
		}
		if (id2 > 0) {
			cs_entry(Z, 2*m + i, 2*m +  id2 - 1, 1);
		}
	}
	cs * Z_ccf = cs_compress(Z);
	return Z_ccf;
}


void writeDoubleArrayToDisk(double * doubleArr, int length, string filename) {
	ofstream ofs;
	ofs.open(filename.c_str());
	ofs.precision(numeric_limits<double>::digits10 + 1);
	for (int i = 0; i < length; i++) {
		ofs << doubleArr[i] << "\r\n";
	}
	ofs.close();
}


cs * fromDoubleArrayPointerToCompressedColumnVector(double * z_hat, int length_z_hat){
	cs * z_hat_triplet = cs_spalloc(length_z_hat, 1, length_z_hat * 1, 1, 1);

	for (int k = 0; k < length_z_hat; k++) 
		cs_entry(z_hat_triplet, k, 0, z_hat[k]);

	cs * z_hat_ccf = cs_compress(z_hat_triplet);

	return z_hat_ccf;
}


cs * constructRotationMatrix(double * z_hat, double **edges, int m) {
	int row = 2 * m;
	int col = 2 * m;
	int id;
	double tmpOrientation;
	double c_tmp, s_tmp;
	cs * rotationMatrix = cs_spalloc(row, col, row * 2, 1, 1); // at most 2 nonzero elements for each row
	for (int k = 0; k < m; k++) {
		id = edges[k][0];
		if (id == 0) { // it is the anchor
			tmpOrientation = 0;// the anchor has orientation 0 by convention
		} else { // it is not the anchor
			tmpOrientation = z_hat[row + id - 1];
		}
		//tmpOrientation_normalized = standardizeOrientation(tmpOrientation);
		c_tmp =  cos(tmpOrientation);
		s_tmp =  sin(tmpOrientation);
		cs_entry(rotationMatrix, 2 * k, 2 * k, 		c_tmp);
		cs_entry(rotationMatrix, 2 * k, 2 * k + 1, -1 * s_tmp);
		cs_entry(rotationMatrix, 2 * k + 1, 2 * k,	s_tmp);
		cs_entry(rotationMatrix, 2 * k + 1, 2 * k + 1,	c_tmp);
	}
	cs * rotationMatrix_ccf = cs_compress(rotationMatrix);
	return rotationMatrix_ccf;
}


cs* constructIdentity(int u) {
	cs* eye = cs_spalloc(u, u, u, 1, 1);
	for (int i = 0; i < u; i++) {
		cs_entry(eye, i, i, 1);
	}
	cs * eye_ccf = cs_compress(eye);
	return eye_ccf;
}


cs* constructZeroMatrix(int row, int col) {
	cs* zero = cs_spalloc(row, col, 1, 1, 1);
	cs_entry(zero, 0, 0, 0);
	cs * zero_ccf = cs_compress(zero);
	return zero_ccf;
}

cs* constructOneMatrix(int row, int col) {
	cs* ones = cs_spalloc(row, col, row * col, 1, 1);
	for(int i=0; i < row; i++){
		for(int j=0; j < col; j++){
			cs_entry(ones, i, j, 1);
		}
	}
	cs * ones_ccf = cs_compress(ones);
	return ones_ccf;
}


cs * constructOmega(vector<edge> vEdges) {
	unsigned int m = vEdges.size();
	int m3 = 3*m;
	cs* Omega = cs_spalloc(m3, m3, m3 * 3, 1, 1); // at most 3 nonzero elements for each row
	
	for (unsigned int i = 0; i < m; i++) {
		//Relative position information
		cs_entry(Omega, 2 * i, 2 * i, vEdges.at(i).t_ff);
		cs_entry(Omega, 2 * i + 1, 2 * i, vEdges.at(i).t_fs);
		cs_entry(Omega, 2 * i, 2 * i + 1, vEdges.at(i).t_fs);
		cs_entry(Omega, 2 * i + 1, 2 * i + 1, vEdges.at(i).t_ss);
		//Relative orientation information
		cs_entry(Omega, 2*m + i, 2*m + i, vEdges.at(i).t_rr);
		// Correlation terms
		cs_entry(Omega, 2 * i, 2*m + i, vEdges.at(i).t_fr);
		cs_entry(Omega, 2 * i + 1, 2*m + i, vEdges.at(i).t_sr);
		
		cs_entry(Omega,  2*m + i, 2 * i,  vEdges.at(i).t_fr);
		cs_entry(Omega, 2*m + i, 2 * i + 1, vEdges.at(i).t_sr);

	}
	cs * Omega_ccf = cs_compress(Omega);
	return Omega_ccf;
}


cs * constructB(double ** edges, int m, int n) {

	int row = 2 * m + n;
	int col = 3 * n;

	cs* B = cs_spalloc(row, col, 4 * m + n, 1, 1);
	// max nz = 4m (for the augmented incidence matrix) + n (for the identity matrix)

	for (int k = 0; k < m; k++) { // for each edge
		int id1 = edges[k][0];
		int id2 = edges[k][1];
		if (id1>0){ // it is not the anchor node
			cs_entry(B, 2 * k, 2 * (id1-1), -1);
			cs_entry(B, 2 * k + 1, 2 * (id1-1) + 1, -1);
		}
		if (id2>0){ // it is not the anchor node
		cs_entry(B, 2 * k, 2 * (id2 -1) , 1);
		cs_entry(B, 2 * k + 1, 2 * (id2-1) + 1, 1);
		}
	}
	for (int k = 0; k < n; k++) { // we fill the lower right part with I_n
		cs_entry(B, 2 * m + k, 2 * n + k, 1);
	}
	cs * B_ccf = cs_compress(B);
	return B_ccf;
}


cs * constructNegJacobianMatrix(double * z_hat,  double ** edges, int m, int n) {
	int row = 2 * m;
	int col = n;
	// Jacobian has dimension 2m times n, where n is the number of poses-1
	int id;
	double tmp1, tmp2;
	double eleTheta, c_th, s_th; // this will contain the angle of the node labeled with id
	int idxTheta;
	cs* partJacobianMatrix = cs_spalloc(row, col, row, 1, 1); // at most 1 nonzero element for each row

	for (int i = 0; i < m; i++) { // for each edge
		id = edges[i][0]; // id of the tail of the edge

		if (id != 0) { // if the considered node is not the anchor node

			idxTheta = id - 1;
			eleTheta = z_hat[row + idxTheta];
			c_th = cos(eleTheta);
			s_th = sin(eleTheta);

			tmp1 = -s_th * z_hat[2 * i] - c_th * z_hat[2 * i + 1];
			cs_entry(partJacobianMatrix, 2 * i, idxTheta, -tmp1); //

			tmp2 = c_th * z_hat[2 * i] - s_th * z_hat[2 * i + 1];
			cs_entry(partJacobianMatrix, 2 * i + 1, idxTheta, -tmp2); //tmp2
		} 
		//else: nothing - the jacobian does not include the column corresponding to the anchor node,
		// therefore here we have to do nothing	
	}
	cs * partJacobianMatrix_ccf = cs_compress(partJacobianMatrix);
	return partJacobianMatrix_ccf;
}


cs * constructTinv(double * z_hat, double **edges, int m, int n) {
	int row = 2 * m + n;
	int col = 2 * m + n;
	int id;
	double tmp1, tmp2, c1, s1;
	double tmpOrientation;

	cs * Tinv = cs_spalloc(row, col, 6*m+n, 1, 1); 
	// at most 3 nonzero elements for the first 2m rows, then 1 for the last n rows
	for (int k = 0; k < m; k++) {
		id = edges[k][0];
		if (id == 0) { // it is the anchor
			tmpOrientation = 0;// the anchor has orientation 0 by convention
			// upper left R'
			c1 = cos(tmpOrientation);
			s1 = sin(tmpOrientation); 
			cs_entry(Tinv, 2 * k, 2 * k, 		c1);
			cs_entry(Tinv, 2 * k, 2 * k + 1,  	s1);
			cs_entry(Tinv, 2 * k + 1, 2 * k, 	-s1);
			cs_entry(Tinv, 2 * k + 1, 2 * k + 1,	c1);

		} else { // it is not the anchor
			tmpOrientation = z_hat[row + id - 1];
			// upper left R'
			c1 = cos(tmpOrientation);
			s1 = sin(tmpOrientation); 
			double s_tmp = sin(tmpOrientation);
			double c_tmp = cos(tmpOrientation);
			cs_entry(Tinv, 2 * k, 2 * k, 		c1);
			cs_entry(Tinv, 2 * k, 2 * k + 1,  	s1);
			cs_entry(Tinv, 2 * k + 1, 2 * k, 	-s1);
			cs_entry(Tinv, 2 * k + 1, 2 * k + 1,	c1);
			// Now retrive the Negative Jacobian
			tmp1 = s_tmp * z_hat[2 * k] + c_tmp * z_hat[2 * k + 1];
			tmp2 = -c_tmp * z_hat[2 * k] + s_tmp * z_hat[2 * k + 1];
			cs_entry(Tinv, 2 * k, 		2*m + id, c1 * tmp1 + s1 * tmp2); 
			cs_entry(Tinv, 2 * k + 1, 	2*m + id, -s1 * tmp1 + c1 * tmp2); 
		}
	}

	for (int k = 0; k < n; k++) 
		cs_entry(Tinv, 2*m + k, 2*m + k, 1); 

	cs * Tinv_ccf = cs_compress(Tinv);
	return Tinv_ccf;
}



cs * combineTwoSparseMatrixLeftAndRight(cs*L, cs*R) {
	int rowL = L->m;
	int colL = L->n;
	int rowR = R->m;
	int colR = R->n;
	int row = rowL;
	int col = colL + colR;

	if (rowL != rowR) {
		cerr << "can't be combined...row number should match." << endl;
	}

	int nzmax = L->nzmax + R->nzmax;  // changed S
	cs * M = cs_spalloc(row, col, nzmax, 1, 1);

	for (int ptr = 0; ptr < colL; ptr++) {
		int indexLow = *(L->p + ptr);
		int indexHigh = *(L->p + ptr + 1) - 1;
		for (int idx = indexLow; idx <= indexHigh; idx++) {
			int row = *(L->i + idx);
			double val = *(L->x + idx);
			cs_entry(M, row, ptr, val);
		}
	}
	for (int ptr = 0; ptr < colR; ptr++) {
		int indexLow = *(R->p + ptr);
		int indexHigh = *(R->p + ptr + 1) - 1;
		for (int idx = indexLow; idx <= indexHigh; idx++) {
			int row = *(R->i + idx);
			double val = *(R->x + idx);
			cs_entry(M, row, colL + ptr, val);
		}
	}

	cs * M_ccf = cs_compress(M);
	return M_ccf;
}


cs * combineTwoSparseMatrixUpAndLow(cs* U, cs* D) {
	int rowU = U->m;
	int colU = U->n;
	int rowD = D->m;
	int colD = D->n;
	int row = rowU + rowD;
	int col = colU;

	if (colU != colD) {
		cerr << "can't be combined...column number should match." << endl;
	}

	int nzmax = U->nzmax + D->nzmax; // changed S
	cs * M = cs_spalloc(row, col, nzmax, 1, 1);

	for (int ptr = 0; ptr < colU; ptr++) {
		int indexLow = *(U->p + ptr);
		int indexHigh = *(U->p + ptr + 1) - 1;
		for (int idx = indexLow; idx <= indexHigh; idx++) {
			int row = *(U->i + idx);
			if (row >= rowU)
				cout << "error: row= " << row << endl;
			double val = *(U->x + idx);
			cs_entry(M, row, ptr, val);
		}
	}
	for (int ptr = 0; ptr < colD; ptr++) {
		int indexLow = *(D->p + ptr);
		int indexHigh = *(D->p + ptr + 1) - 1;
		for (int idx = indexLow; idx <= indexHigh; idx++) {
			int row = *(D->i + idx);
			double val = *(D->x + idx);
			cs_entry(M, row + rowU, ptr, val);
		}
	}

	cs * M_ccf = cs_compress(M);
	return M_ccf;
}


cs * combineTwoSparseMatrixDiag(cs* U, cs* D) {
	int rowU = U->m;
	int colU = U->n;
	int rowD = D->m;
	int colD = D->n;
	int row = rowU + rowD;
	int col = colU + colD;

	int nzmax = U->nzmax + D->nzmax; // changed S
	cs * M = cs_spalloc(row, col, nzmax, 1, 1);

	for (int ptr = 0; ptr < colU; ptr++) {
		int indexLow = *(U->p + ptr);
		int indexHigh = *(U->p + ptr + 1) - 1;
		for (int idx = indexLow; idx <= indexHigh; idx++) {
			int row = *(U->i + idx);
			double val = *(U->x + idx);
			cs_entry(M, row, ptr, val);
		}
	}
	for (int ptr = 0; ptr < colD; ptr++) {
		int indexLow = *(D->p + ptr);
		int indexHigh = *(D->p + ptr + 1) - 1;
		for (int idx = indexLow; idx <= indexHigh; idx++) {
			int row = *(D->i + idx);
			double val = *(D->x + idx);
			cs_entry(M, row + rowU, ptr + colU, val);
		}
	}
	cs * M_ccf = cs_compress(M);
	return M_ccf;
}


cs * constructT(cs* R, int m, int n) {

	int rowR = 2*m;
	int colR = 2*m;
	int row = rowR + n;
	int col = colR + n;

	cs * M = cs_spalloc(row, col, 2*rowR + n, 1, 1);

	for (int ptr = 0; ptr < colR; ptr++) {
		int indexLow = *(R->p + ptr);
		int indexHigh = *(R->p + ptr + 1) - 1;
		for (int idx = indexLow; idx <= indexHigh; idx++) {
			int row = *(R->i + idx);
			double val = *(R->x + idx);
			cs_entry(M, row, ptr, val);
		}
	}
	for (int ptr = 0; ptr < n; ptr++)  // this builds the lower right identity matrix
		cs_entry(M, rowR + ptr, rowR + ptr, 1);

	cs * M_ccf = cs_compress(M);
	return M_ccf;
}



/* C = A*vect */
cs *cs_multiplyMatVect(const cs *A, const cs *B)
{
    csi p, j, nz = 0, *Cp, *Ci, *Bp, m, n, *w, values, *Bi ;
    double *x, *Bx, *Cx ;
    cs *C ;
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;      /* check inputs */
    if (A->n != B->m) return (NULL) ;
    m = A->m ;
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ;
    w = (ptrdiff_t*) cs_calloc (m, sizeof (csi)) ;                    /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? (double*) cs_malloc (m, sizeof (double)) : NULL ; /* get workspace */

    C = cs_spalloc (m, n, m, values, 0) ;        /* allocate result */

    if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !cs_sprealloc (C, 2*(C->nzmax)+m))
        {
            return (cs_done (C, w, x, 0)) ;             /* out of memory */
        }
        Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
        Cp [j] = nz ;                   /* column j of C starts here */
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            nz = cs_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
        }
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    cs_sprealloc (C, 0) ;               /* remove extra space from C */
    return (cs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}


cs *cs_multiplyFlexible(const cs *A, const cs *B, int nzmaxC)
// is the same as the standard cs sparse matrix multiplication, but allows to fix the
// number of nonzero elements of C = A*B, for optimized memory allocation
{
    csi p, j, nz = 0, *Cp, *Ci, *Bp, m, n, *w, values, *Bi ;
    double *x, *Bx, *Cx ;
    cs *C ;
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;      /* check inputs */
    if (A->n != B->m) return (NULL) ;
    m = A->m ;
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ;
    w = (ptrdiff_t*) cs_calloc (m, sizeof (csi)) ;                    /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? (double*) cs_malloc (m, sizeof (double)) : NULL ; /* get workspace */

    C = cs_spalloc (m, n, nzmaxC, values, 0) ;        /* allocate result */

    if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !cs_sprealloc (C, 2*(C->nzmax)+m))
        {
            return (cs_done (C, w, x, 0)) ;             /* out of memory */
        }
        Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
        Cp [j] = nz ;                   /* column j of C starts here */
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            nz = cs_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
        }
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    cs_sprealloc (C, 0) ;               /* remove extra space from C */
    return (cs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}


#endif /* THREEPHASE_H_ */
