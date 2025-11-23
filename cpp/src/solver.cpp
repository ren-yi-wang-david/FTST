#include <mkl.h>
#include <vector>
#include <iostream>
#include "mkl_solver.hpp"
void thermal_solver::mkl_solver_csr(
	std::vector<int> &ArowNEntry,
	std::vector<std::vector<int>> &ArowCols,
	std::vector<std::vector<double>> &ArowVals,
	std::vector<double> &x,
	std::vector<double> &b)
{
	int n = ArowNEntry.size();
	int nnz = 0;

	// Build row_ptr
	std::vector<int> ia(n + 1);
	ia[0] = 1; // Fortran-style 1-based
	for (int i = 0; i < n; i++)
	{
		nnz += ArowNEntry[i];
		ia[i + 1] = nnz + 1;
	}

	// Build ja + a
	std::vector<int> ja(nnz);
	std::vector<double> a(nnz);

	int idx = 0;
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < ArowNEntry[i]; k++)
		{
			ja[idx] = ArowCols[i][k] + 1; // 1-based
			a[idx] = ArowVals[i][k];
			idx++;
		}
	}

	void *pt[64] = {0};
	int iparm[64] = {0};
	for (int i = 0; i < 64; i++)
		iparm[i] = 0;

	int maxfct = 1, mnum = 1, phase, error = 0, msglvl = 0;
	int mtype = 11;
	int nrhs = 1;

	iparm[0] = 1;

	x.resize(n);

	// Phase 11
	phase = 11;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
			&n, a.data(), ia.data(), ja.data(),
			nullptr, &nrhs, iparm, &msglvl,
			nullptr, nullptr, &error);

	// Phase 22
	phase = 22;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
			&n, a.data(), ia.data(), ja.data(),
			nullptr, &nrhs, iparm, &msglvl,
			nullptr, nullptr, &error);

	// Phase 33
	phase = 33;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
			&n, a.data(), ia.data(), ja.data(),
			nullptr, &nrhs, iparm, &msglvl,
			b.data(), x.data(), &error);

	// Phase -1
	phase = -1;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
			&n, a.data(), ia.data(), ja.data(),
			nullptr, &nrhs, iparm, &msglvl,
			nullptr, nullptr, &error);
}
