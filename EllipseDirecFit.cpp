#include<iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
void main()
{
	/** Ellipse Direct Fit */

	int XY[][2] = { { 330, 275 }, { 331, 275 }, { 328, 276 }, { 329, 276 }, { 332, 276 }, { 333, 276 }, { 334, 276 }, { 326, 277 }, { 327, 277 }, { 335, 277 }, { 325, 278 }, { 336, 278 }, { 325, 279 }, { 337, 279 }, { 325, 280 }, { 337, 280 }, { 325, 281 }, { 337, 281 }, { 325, 282 }, { 337, 282 }, { 325, 283 }, { 337, 283 }, { 325, 284 }, { 337, 284 }, { 326, 285 }, { 336, 285 }, { 326, 286 }, { 335, 286 }, { 327, 287 }, { 328, 287 }, { 333, 287 }, { 334, 287 }, { 329, 288 }, { 330, 288 }, { 331, 288 }, { 332, 288 } };

	const int RowSize = sizeof(XY) / sizeof(XY[0]);
	//int col = sizeof(XY[0])/sizeof(int);

	double centroid[2];

	double xSum = 0.0, ySum = 0.0;

	for (int i = 0; i < RowSize; i++)
	{
		xSum += XY[i][0];
		ySum += XY[i][1];
	}

	double xMean = xSum / RowSize;
	double yMean = ySum / RowSize;

	centroid[0] = xMean;
	centroid[1] = yMean;

	double D1[RowSize][3], D2[RowSize][3]; /***************************************/

	for (int i = 0; i < RowSize; i++)
	{
		D1[i][0] = pow((XY[i][0] - centroid[0]), 2);
		D1[i][1] = (XY[i][0] - centroid[0])*(XY[i][1] - centroid[1]);
		D1[i][2] = pow((XY[i][1] - centroid[1]), 2);

		D2[i][0] = (XY[i][0] - centroid[0]);
		D2[i][1] = (XY[i][1] - centroid[1]);
		D2[i][2] = 1;
	}

	double transposeD1[3][RowSize];/******************************************************/
	double transposeD2[3][RowSize]; /************************************/

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < RowSize; j++) /****************************************/
		{
			transposeD1[i][j] = D1[j][i];
			transposeD2[i][j] = D2[j][i];
		}
	}

	double S1[3][3] = { 0 };
	double S2[3][3] = { 0 };
	double S3[3][3] = { 0 };

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < RowSize; k++) /************************************************/
			{
				S1[i][j] += transposeD1[i][k] * D1[k][j];

				S2[i][j] += transposeD1[i][k] * D2[k][j];

				S3[i][j] += transposeD2[i][k] * D2[k][j];
			}

		}
	}

	double inverseS3[3][3];
	double determinant = 0;

	for (int i = 0; i < 3; i++)
	{
		determinant = determinant + (S3[0][i] * (S3[1][(i + 1) % 3] * S3[2][(i + 2) % 3] - S3[1][(i + 2) % 3] * S3[2][(i + 1) % 3]));
	}

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			inverseS3[i][j] = ((S3[(j + 1) % 3][(i + 1) % 3] * S3[(j + 2) % 3][(i + 2) % 3]) - (S3[(j + 1) % 3][(i + 2) % 3] * S3[(j + 2) % 3][(i + 1) % 3])) / determinant;
		}
	}

	double transposeS2[3][3];

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			transposeS2[i][j] = S2[j][i];
		}
	}


	double T[3][3] = { 0 };

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++) /*************************************/
			{
				T[i][j] += -inverseS3[i][k] * transposeS2[k][j];
			}

		}
	}

	double M[3][3], tmp[3][3] = { 0 };

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				tmp[i][j] += S2[i][k] * T[k][j];
			}
			M[i][j] = S1[i][j] + tmp[i][j];
		}
	}

	double newM[3][3];

	for (int i = 0; i < 3; i++)
	{

		newM[0][i] = M[2][i] / 2;
		newM[1][i] = -M[1][i];
		newM[2][i] = M[0][i] / 2;
	}

	MatrixXd m(3, 3);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			m(i, j) = newM[i][j];
		}
	}

	std::cout << m << std::endl;

	EigenSolver<MatrixXd> es(m);

	cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;

	MatrixXd D = es.pseudoEigenvalueMatrix();
	MatrixXd V = es.pseudoEigenvectors();

	cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
	cout << "The pseudo-eigenvector matrix V is:" << endl << V << endl;
	cout << "***********************" << endl << V(0, 1) << endl;
	cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;

	double cound[3] = { 0 };

	int test1 = sizeof(cound) / sizeof(cound[0]);

	for (int i = 0; i < sizeof(cound) / sizeof(cound[0]); i++)
	{
		cound[i] += 4*V(0, i)*V(2,i)-pow(V(1,i),2);
	
	}

	int index;
	for (int i = 0; i < sizeof(cound) / sizeof(cound[0]); i++)
	{
		if (cound[i]>0)
		{
			index = i;
		}
	}

	double A1[3];

	for (int i = 0; i < sizeof(A1) / sizeof(A1[0]); i++)
	{
		A1[i] = V(i, index);
	}

	double A[6];

	A[0] = A1[0];
	A[1] = A1[1];
	A[2] = A1[2];

	double tmp2[3] = {0};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			
			tmp2[i] += T[i][j]*A1[j] ;

		}
	}

	A[3] = tmp2[0];
	A[4] = tmp2[1];
	A[5] = tmp2[2];

	double A3 = A[3] - 2 * A[0] * centroid[0] - A[1] * centroid[1];
	double A4 = A[4] - 2 * A[2] * centroid[1] - A[1] * centroid[0];
	double A5 = A[5] + A[0] * centroid[0] * centroid[0] + A[2] * centroid[1] * centroid[1] + A[1] * centroid[0] * centroid[1] - A[3] * centroid[0] - A[4] * centroid[1];

	A[3] = A3;
	A[4] = A4;
	A[5] = A5;

	double normA = 0;

	for (int i = 0; i < 6; i++)
	{
		normA = pow(A[i], 2);
	}

	normA = sqrt(normA);

	for (int i = 0; i < 6; i++)
	{
		A[i] = A[i] / normA;
	}

	double a = A[0];
	double b = A[1];
	double c = A[2];
	double d = A[3];
	double e = A[4];
	double f = A[5];

	/** https://math.stackexchange.com/a/820896 */

	double q = 64 * ((f*(4 * a*c - b*b) - a*e*e + b*d*e - c*d*d) / pow((4 * a*c - b *b), 2)); //coefficient normalizing factor

	double s = 1.0 / 4.0* (sqrt(abs(q)*sqrt(b *b + (a - c) *(a - c)))); //distance between center and focal point

	double r_max = 1.0 / 8.0 *(sqrt(2 * abs(q)*sqrt(b *b + (a - c) *(a - c)) - 2 * q*(a + c))); //semi - major axis length

	double r_min = sqrt(r_max*r_max - s *s); //semi - minor axis length

	int x0 = round((b*e - 2 * c*d) / (4 * a*c - b*b));
	int y0 = round((b*d - 2 * a*e) / (4 * a*c - b*b));

	cout << "Rmax: " << r_max << " Rmin: " << r_min << " Center:(" << x0 << "," << y0 << ")";

}
