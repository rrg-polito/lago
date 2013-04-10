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
#include <cstddef>
#include <limits>

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

double standardizeOrientation(double delta) {
	double tmp = fmod(delta, 2 * PI);
	if (tmp < -PI) {
		tmp += 2 * PI;
	} else if (tmp > PI) {
		tmp -= 2 * PI;
	}
	return tmp;
}

// this function is useful for computing statistics on the computational effort
double calculateSeconds(struct timeval time) {
	double seconds;
	seconds = time.tv_sec + 1e-6 * time.tv_usec;
	return seconds;
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

cs * constructRegularizedOrientationVector(double ** edges, int m, int n){

	int rowNum = m;
	int colNum = 1;
	int nzmax = rowNum * colNum;
	int values = 1;
	int triplet = 1;
	cs *OrientationVector = cs_spalloc(rowNum, colNum, nzmax, values, triplet);
	double tot_orientation;

	for (int i = 0; i < n; i++) {
		//  these do not need regularization
		cs_entry(OrientationVector, i, 0, edges[i][4]); 
	}

	for (int i = n; i < m; i++) { // these have to be regularized
		{
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
		cs_entry(OrientationVector, i, 0, edges[i][4] - tot_orientation); 
		}
	}
	cs * OrientationVector_ccf = cs_compress(OrientationVector);
	return OrientationVector_ccf;
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
        //    assign -1s for forward and 1s on the contrary;
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

cs * constructPdeltaInvMatrix(vector<edge> vEdges) {
	int u = vEdges.size();
	cs* PDeltaInv = cs_spalloc(u, u, u, 1, 1);
	for (int i = 0; i < u; i++) {
		double val = vEdges.at(i).t_rr;
		cs_entry(PDeltaInv, i, i, val);
	}
	cs * PDeltaInv_ccf = cs_compress(PDeltaInv);
	return PDeltaInv_ccf;
}

cs * constructReducedIncidenceMatrixTran(double ** edges, int m, int n) {
	int rowLen = m;
	int colLen = n; //it is the reduced incidence matrix
	cs * reducedIncidenceMatrixTran = cs_spalloc(rowLen, colLen, rowLen * 2, 1, 1);
	for (int i = 0; i < m; i++) { // for each edge
		int id1 = edges[i][0]; // we pick the tail of the i-th edge
		int id2 = edges[i][1]; // we pick the head of the i-th edge
		if (id1 > 0) {
			cs_entry(reducedIncidenceMatrixTran, i, id1 - 1, -1);
		}
		if (id2 > 0) {
			cs_entry(reducedIncidenceMatrixTran, i, id2 - 1, 1);
		}
	}
	cs * reducedIncidenceMatrixTran_ccf = cs_compress(reducedIncidenceMatrixTran);
	return reducedIncidenceMatrixTran_ccf;
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

cs * constructPositionColumnVector(double ** edges, int m) {
	int row = 2 * m;
	cs * positionVector = cs_spalloc(row, 1, row * 1, 1, 1);

	for (int k = 0; k < m; k++) {
		cs_entry(positionVector, 2 * k, 0, edges[k][2]);
		cs_entry(positionVector, 2 * k + 1, 0, edges[k][3]);
	}
	cs * positionVector_ccf = cs_compress(positionVector);
	return positionVector_ccf;
}

cs * constructRotationMatrixWithEstimatedOrientations(double * orientationVector, double ** edges, int m) {
	int row = 2*m;
	int col = 2*m;
	int id;
	double tmpOrientation, c_tmp, s_tmp;

	cs * rotationMatrix = cs_spalloc(row, col, row * 2, 1, 1);
	for (int k = 0; k < m; k++) {
		id = edges[k][0];
		if (id == 0) { // it is the anchor
			tmpOrientation = 0;// the anchor has orientation 0 by convention
		} else { // it is not the anchor
			tmpOrientation = orientationVector[id - 1];
		}
		s_tmp = sin(tmpOrientation);
		c_tmp = cos(tmpOrientation);

		cs_entry(rotationMatrix, 2 * k, 2 * k, 		c_tmp);
		cs_entry(rotationMatrix, 2 * k, 2 * k + 1, 	-s_tmp);
		cs_entry(rotationMatrix, 2 * k + 1, 2 * k,	s_tmp);
		cs_entry(rotationMatrix, 2 * k + 1, 2 * k + 1,   c_tmp);
	}
	cs * rotationMatrix_ccf = cs_compress(rotationMatrix);
	return rotationMatrix_ccf;
}

cs * constructZVector(cs * z_position, double* theta_estimated, int length) {
	int rowZ = z_position->m;
	int row = rowZ + length;
	cs * z = cs_spalloc(row, 1, row, 1, 1);

	for (int idx = 0; idx < rowZ; idx++) {
		cs_entry(z, idx, 0, *(z_position->x + idx));
	}
	for (int i = 0; i < length; i++) {
		double theta = theta_estimated[i];
		cs_entry(z, rowZ + i, 0, theta);
	}
	cs* z_ccf = cs_compress(z);
	return z_ccf;
}

cs* constructPDeltaLInverseMatrix(vector<edge> vEdges) {
	int u = 2 * vEdges.size();
	cs * PDeltaLInv = cs_spalloc(u, u, 4 * u, 1, 1);

	for (unsigned int i = 0; i < vEdges.size(); i++) {
		cs_entry(PDeltaLInv, 2 * i, 2 * i, vEdges.at(i).t_ff);
		cs_entry(PDeltaLInv, 2 * i + 1, 2 * i, vEdges.at(i).t_fs);
		cs_entry(PDeltaLInv, 2 * i, 2 * i + 1, vEdges.at(i).t_fs);
		cs_entry(PDeltaLInv, 2 * i + 1, 2 * i + 1, vEdges.at(i).t_ss);
	}
	if (!PDeltaLInv)
		cout << "null matrix" << endl;

	cs * PDeltaLInv_ccf = cs_compress(PDeltaLInv);
	return PDeltaLInv_ccf;
}


cs * constructBT(double ** edges, int m, int n) {
	int row = 2 * m + n;
	int col = 2 * n + n;

	cs * BT = cs_spalloc(row, col, 4*m + n, 1, 1);

	for (int k = 0; k < m; k++) {
		int id1 = edges[k][0];
		int id2 = edges[k][1];
		if (id1>0){
		cs_entry(BT, 2 * k, 2 * (id1-1), -1);
		cs_entry(BT, 2 * k + 1, 2 * (id1-1) + 1, -1);
		}
		if (id2>0){
		cs_entry(BT, 2 * k, 2 * (id2-1), 1);
		cs_entry(BT, 2 * k + 1, 2 * (id2-1) + 1, 1);
		}
	}

	for (int k = 0; k < n; k++) {
		cs_entry(BT, 2 * m + k, 2 * n + k, 1);
	}
	cs * BT_ccf = cs_compress(BT);
	return BT_ccf;
}

cs * constructJacobianMatrix(vector<edge> vEdges, double * theta,
		double *position_meas, vector<pose> vPoses) {
	int row = 2 * vEdges.size();
	int col = vPoses.size() - 1;
	// Jacobian has dimension 2m times n, where n is the number of poses-1
	int id;
	double tmp1, tmp2;
	double eleTheta; // this will contain the angle of the node labeled with id
	int idxTheta;
	cs* partJacobianMatrix = cs_spalloc(row, col, row, 1, 1);

	for (unsigned int i = 0; i < vEdges.size(); i++) { // for each edge
		id = vEdges.at(i).id1; // id of the tail of the edge

		if (id != vPoses.at(0).id) { // if the considered node is not the anchor node

			idxTheta = id - 1;
			eleTheta = theta[idxTheta];
			tmp1 = -sin(eleTheta) * position_meas[2 * i] - cos(eleTheta)
					* position_meas[2 * i + 1];
			cs_entry(partJacobianMatrix, 2 * i, idxTheta, -tmp1); //

			tmp2 = cos(eleTheta) * position_meas[2 * i] - sin(eleTheta)
					* position_meas[2 * i + 1];
			cs_entry(partJacobianMatrix, 2 * i + 1, idxTheta, -tmp2); //tmp2
		} else {
			// the jacobian does not include the column corresponding to the anchor node,
			// therefore here we have to do nothing
		}
	}
	cs * partJacobianMatrix_ccf = cs_compress(partJacobianMatrix);
	return partJacobianMatrix_ccf;
}

cs * constructNegJacobianMatrix(double ** edges, double * theta, double *position_meas, int m, int n) {
	int row = 2 * m;
	int col = n;
	// Jacobian has dimension 2m times n, where n is the number of poses-1
	int id;
	double tmp1, tmp2, c1, s1;
	double eleTheta; // this will contain the angle of the node labeled with id
	int idxTheta;
	cs* partJacobianMatrix = cs_spalloc(row, col, row, 1, 1);

	for (int i = 0; i < m; i++) { // for each edge
		id = edges[i][0]; // id of the tail of the edge

		if (id != 0) { // if the considered node is not the anchor node

			idxTheta = id - 1;
			eleTheta = theta[idxTheta];
			s1 = sin(eleTheta);
			c1 = cos(eleTheta);
			tmp1 = s1 * position_meas[2 * i] + c1 * position_meas[2 * i + 1];
			cs_entry(partJacobianMatrix, 2 * i, idxTheta, tmp1); 

			tmp2 = - c1 * position_meas[2 * i] + s1 * position_meas[2 * i + 1];
			cs_entry(partJacobianMatrix, 2 * i + 1, idxTheta, tmp2); 
		} 
		//else: do nothing - the jacobian does not include the column corresponding to the anchor node,
		// therefore here we have to do nothing

	}
	cs * partJacobianMatrix_ccf = cs_compress(partJacobianMatrix);
	return partJacobianMatrix_ccf;
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

	int nzmax = L->nzmax + R->nzmax;
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

	return M;
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

	int nzmax = U->nzmax + D->nzmax;
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

	return M;
}


/* C = A*vect */
cs *cs_multiplyMatVect (const cs *A, const cs *B)
{
    csi p, j, nz = 0, *Cp, *Ci, *Bp, m, n, *w, values, *Bi ;
    double *x, *Bx, *Cx ;
    cs *C ;
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;      /* check inputs */
    if (A->n != B->m) return (NULL) ;
    m = A->m ;
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x;
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

#endif /* THREEPHASE_H_ */
