#include "pch.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <omp.h>

using namespace std;

#define PI 3.14159265
double al = 0.0,
vy = 0.5,
u = 0.93,
gam = pow((1 - u * u), 2),
hz = 0.06,
ht = 0.005,
bet = 0.03;

const int vx = 1,
kk = 15,
mm = 13,
z0 = -10,
kz = 800,
j = 5;

double aa1 = pow((ht / hz), 2);
double aa2 = pow(ht, 2);
int nn = 1, lll = 200, ll = 0;

//double* f = new double[kk];
double* H = new double[kk];
double A[mm][kk];
double HH[kk];
double* b = new double[mm];


int main()
{
	omp_set_num_threads(4);
	double f[16] = { 0, 1, 0.001672, 0.189, 0.0002537, 0.083, 0.00004931, 0.044, 0.00001042, 0.027, 0.000002295, 0.017, 0.000005185, 0.011, 0.0000001191, 0.007831 };
	//double f[16] = { 0, 1, 0.1672, 0.189, 0.02537, 0.083, 0.004931, 0.044, 0.001042, 0.027, 0.0002295, 0.017, 0.0005185, 0.011, 0.00001191, 0.007831 };

	std::cout << " program" << std::endl;
	float** qqx1 = new float*[lll];
	float** qqx2 = new float*[lll];
	float** qqy1 = new float*[lll];
	float** qqy2 = new float*[lll];
	float** qqx = new float*[lll];
	float** qqy = new float*[lll];
	float** eeex = new float *[lll];
	float** per = new float*[lll];
	float sum1 = 0.0,
		sum2 = 0.0,
		sum3 = 0.0,
		d = 0.0;

	for (int i = 0; i < kz; i++)
	{
		qqx1[i] = new float[kz];
		qqx2[i] = new float[kz];
		qqy1[i] = new float[kz];
		qqy2[i] = new float[kz];
		qqx[i] = new float[kz];
		qqy[i] = new float[kz];
		eeex[i] = new float[kz];
		per[i] = new float[kz];
	}

	float* ks = new float[kz];
	for (int i = 0; i <= kz; i++)
	{
		ks[i] = z0 + i * hz;
		//std::cout << "ks[i] = "<<ks[i] << std::endl
	}
#pragma omp parallel
	for (int k = 1; k <= lll; k++)
		for (int i = 0; i <= kz; i++)
		{
			qqx1[k][i] = 4 * exp(-pow((ks[i] / gam), 2))*exp(-bet * pow((k - ll), 2));
			qqx2[k][i] = 4 * exp(-pow(((ks[i] - u * ht) / gam), 2))*exp(-bet * pow((k - ll), 2));
		}
	for (int k = 1; k <= lll; k++)
		for (int i = 0; i <= kz; i++)
		{
			qqy1[k][i] = 0.0;
			qqy2[k][i] = 0.0;
		}

	ofstream file;
	//file.open("C://Users//82G//Desktop//учеба//магистратура//2 семестр///Математическое моделирование для решения междисциплинарных задач//n.csv", ios::out);

	while (nn < 15)
	{
	#pragma omp parallel
		for (int k = 1; k <= lll; k++)
		{
			qqx[k][0] = (4 * qqx2[k][1] - qqx2[k][2]) / 3;
			qqx[k][kz] = (4 * qqx2[k][kz - 1] - qqx2[k][kz - 2]) / 3;
			qqy[k][0] = (4 * qqy2[k][1] - qqy2[k][2]) / 3;
			qqy[k][kz] = (4 * qqy2[k][kz - 1] - qqy2[k][kz - 2]) / 3;
		}
	#pragma omp parallel
		for (int i = 1; i <= kz - 1; i++)
		{
			for (int p = 2; p <= kk; p++)
			{
				sum1 += f[p] * sin(p*(qqx2[1][i] * cos(al) + qqy[1][i] * sin(al)));
				sum2 += f[p] * sin(p*(qqx2[lll][i] * cos(al) + qqy[lll][i] * sin(al)));
			}
			qqx[1][i] = 2 * qqx2[1][i] - qqx1[1][i] + aa1 * vx*(qqx2[1][i + 1] - 2 * qqx2[1][i] + qqx2[1][i - 1]) +
				1.5*aa1*(qqx2[2][i] - qqx2[1][i]) +
				aa2 * cos(al)*(sin(qqx2[1][i] * cos(al) + qqy[1][i] * sin(al))) + sum1;

			qqx[lll][i] = 2 * qqx2[lll][i] - qqx1[lll][i] + aa1 * vx*(qqx2[lll][i + 1] - 2 * qqx2[lll][i] + qqx2[lll][i - 1]) +
				aa1 * (qqx2[lll][i] - qqx2[lll - 1][i])*((-lll + 1 / 2) / lll) +
				aa2 * cos(al)*(sin(qqx2[lll][i] * cos(al) + qqy[lll][i] * sin(al))) + sum2;

			qqy[1][i] = 2 * qqy2[1][i] - qqy1[1][i] + aa1 * vy*(qqy2[1][i + 1] - 2 * qqy2[1][i] + qqy2[1][i - 1]) +
				1.5*aa1*(qqy2[2][i] - qqy2[1][i]) +
				aa2 * sin(al)*(sin(qqx2[1][i] * cos(al) + qqy[1][i] * sin(al))) + sum1;

			qqy[lll][i] = 2 * qqy2[lll][i] - qqy1[lll][i] + aa1 * vy*(qqy2[lll][i + 1] - 2 * qqy2[lll][i] + qqy2[lll][i - 1]) +
				aa1 * (qqy2[lll][i] - qqy2[lll - 1][i])*((-lll + 1 / 2) / lll) +
				aa2 * sin(al)*(sin(qqx2[lll][i] * cos(al) + qqy[lll][i] * sin(al))) + sum2;
		}
	#pragma omp parallel
		for (int k = 2; k <= lll - 1; k++)
			for (int i = 1; i <= kz - 1; i++)
			{
				for (int p = 2; p <= kk; p++)
				{
					sum3 += f[p] * sin(p*(qqx2[k][i] * cos(al) + qqy2[k][i] * sin(al)));
				}
				d = ((k + 1 / 2)*(qqx2[k + 1][i] - qqx2[k][i]) - (k - 1 / 2)*(qqx2[k][i] - qqx2[k - 1][i])) / k;

				qqx[k][i] = 2 * qqx2[k][i] - qqx1[k][i] + aa1 * vx*(qqx2[k][i + 1] - 2 * qqx2[k][i] + qqx2[k][i - 1]) +
					aa1 * d + aa2 * cos(al)*(sin(qqx2[k][i] * cos(al) + qqy2[k][i] * sin(al))) + sum3;
				//cout << "qqx = " << qqx[k][i] << endl;

				qqy[k][i] = 2 * qqy2[k][i] - qqy1[k][i] + aa1 * vy*(qqy2[k][i + 1] - 2 * qqy2[k][i] + qqy2[k][i - 1]) +
					aa1 * d + aa2 * sin(al)*(sin(qqx2[k][i] * cos(al) + qqy2[k][i] * sin(al))) + sum3;
			}
	#pragma omp parallel
		for (int k = 1; k <= lll; k++)
			for (int i = 1; i <= kz; i++)
			{
				qqx1[k][i] = qqx2[k][i];
				qqx2[k][i] = qqx[k][i];
				 
				qqy1[k][i] = qqy2[k][i];
				qqy2[k][i] = qqy[k][i];
				//cout << qqy2[k][i] << endl;
			}
		nn = nn + 1;
		cout << nn << std::endl;
		file << nn << std::endl;
	}

	ofstream file2;
	file2.open("C://Users//kulbina//Downloads//ImpulsNew//15_00.dat", ios::out);

	//qqx

	double max = per[0][0];
	double norm = 0.0;
	//double per = 0.0;


#pragma omp parallel
	for (int k = 1; k <= lll; k++)
		for (int i = 2; i <= kz - 2; i++)
		{
			eeex[k][i] = (-qqx[k][i + 1] + qqx[k][i]) / ht;
			per[k][i] = pow(eeex[k][i], 2);

			//file2 << setw(6)<< per<< std::endl;
			//cout << "eeex["<<k<<"]["<<i<<"] = " << eeex[k][i] << endl;
		}
#pragma omp parallel
	for (int k = 1; k <= lll; k++)
		for (int i = 2; i <= kz - 2; i++)
		{
			if (per[k][i] > max) {
				max = per[k][i]; 
			}
			//cout << "max = " << max << 	
		}

#pragma omp parallel
	for (int k = 0; k <= lll; k++)
	{
		for (int i = 0 ; i <= kz - 2; i++)
		{
			norm = per[k][i] / max;
			if (norm >= 0)
			{
				file2 << setw(6) << norm<< ";";
			}
		}
		file2 << std::endl;
	}

	std::cout << "End of program" << std::endl;

	return 0;
}
