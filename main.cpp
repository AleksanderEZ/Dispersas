#include <mkl.h>
#include "mkl_dss.h"
#include <cstdio>

void printArray(int* data, int length) {
	printf("[ ");
	for (int i = 0; i < length-1; i++)
	{
		printf("%d, ", data[i]);
	}
	printf("%d ]\n", data[length-1]);
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

	opt = 0;
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

int main() {
	actividadPractica4ejercicio1();
	return 0;
}