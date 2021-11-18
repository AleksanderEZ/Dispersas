#include <mkl.h>
#include "mkl_dss.h"
#include <cstdio>

#pragma warning(disable : 4996)

void printArray(int* data, int length) {
	printf("[ ");
	for (int i = 0; i < length-1; i++)
	{
		printf("%d, ", data[i]);
	}
	printf("%d ]\n", data[length-1]);
}

void printArray(double* data, int length) {
	printf("[ ");
	for (int i = 0; i < length - 1; i++)
	{
		printf("%lf, ", data[i]);
	}
	printf("%lf ]\n", data[length - 1]);
}

void actividadPractica4ejercicio1() {
	//Resolver el siguiente sistema Ax=b aplicando la metodología DSS

	//int A[25] = {
	//	1, 0, 0, 0, 0,
	//	2, 3, 0, 0, 0,
	//	4, 0, 5, 0, 0,
	//	0, 6, 0, 7, 0,
	//	0, 0, 8, 0, 9
	//};

	int A[9] = {
		1, 2, 3, 4, 5, 6, 7, 8, 9
	};

	int B[5] = { 1, 2, 3, 4, 5 };

	int rowIndex[5 + 1] = { 0, 1, 3, 5, 7, 9 };
	int colIndex[9] = { 0, 0, 1, 0, 2, 1, 3, 2, 4 };

	const MKL_INT size = 5;
	const MKL_INT nonZeroAmount = 9;

	_MKL_DSS_HANDLE_t handle;
	MKL_INT opt = MKL_DSS_ZERO_BASED_INDEXING + MKL_DSS_SINGLE_PRECISION + MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR;
	int info = dss_create(handle, opt);
	printf("Create: %d\n", info);

	opt = MKL_DSS_NON_SYMMETRIC;
	info = dss_define_structure(handle, opt, rowIndex, size, size, colIndex, nonZeroAmount);
	printf("Structure: %d\n", info);

	int* perm = (int*) malloc(5 * sizeof(int));
	opt = MKL_DSS_AUTO_ORDER;
	info = dss_reorder(handle, opt, 0);
	printf("Reorder: %d\n", info);

	opt = MKL_DSS_INDEFINITE;
	info = dss_factor_real(handle, opt, A);
	printf("Factor: %d\n", info);

	printArray(B, 5);
	opt = MKL_DSS_FORWARD_SOLVE;
	info = dss_solve_real(handle, opt, B, size, B);
	printArray(B, 5);
	printf("Solve real: %d\n", info);

	opt = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR;
	info = dss_delete(handle, opt);
	printf("Delete: %d\n", info);
}

void actividadPractica2ejercicio1() {
	double A[49] = {
		 0,         0,         0,         0,   64.6376,         0,         0,
	75.2755,   80.7525,   98.9732,  108.0870,  101.2359,         0,         0,
		 0,   99.8230,         0,         0,         0,         0,         0,
		 0,         0,         0,         0,         0,         0,         0,
		 0,  114.9331,         0,         0,  104.7711,         0,         0,
	71.8266,         0,         0,         0,   83.9254,         0,   99.7240,
		 0,   72.0374,         0,   61.9563,   69.9837,   92.4787,         0
	};
	
	// CSR
	int job[6] = { 0,0,0,2,40,3 };
	int ia[8];
	int ja[16];
	double acsr[16];
	MKL_INT info;
	const int m = 7;
	int nnz = 16;

	mkl_ddnscsr(job, &m, &m, A, &m, acsr, ja, ia, &info);
	printArray(acsr, 16);
	printArray(ja, 16);
	printArray(ia, 8);
	printf("info: %d\n", info);

	// COO
	int job2[6] = { 0,0,0,2,40,3 };
	int rowind[16];
	int colind[16];
	double acoo[16];

	mkl_dcsrcoo(job2, &m, acsr, ja, ia, &nnz, acoo, rowind, colind, &info);
	printArray(acoo, 16);
	printArray(rowind, 16);
	printArray(colind, 16);
	printf("info: %d\n", info);

	// BSR XD no se puede porque matrix 7x7

	// CSC
	int job3[6] = { 0,0,0,0,0,3 };
	double absr[16];
	int jab[16];
	int iab[8];
	mkl_dcsrcsc(job3, &m, acsr, ja, ia, absr, jab, iab, &info);
	printArray(absr, 16);
	printArray(jab, 16);
	printArray(iab, 8);
	printf("info: %d\n", info);
}

int main() {
	printf("Ejercicio 2-1:\n");
	actividadPractica2ejercicio1();
	printf("Ejercicio 4-1:\n");
	actividadPractica4ejercicio1();
	return 0;
}