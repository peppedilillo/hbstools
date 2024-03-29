#ifndef POISSONFOCUSSES_H
#define POISSONFOCUSSES_H

#if defined(_WIN32)
#  define DLL00_EXPORT_API __declspec(dllexport)
#else
#  define DLL00_EXPORT_API
#endif

#include "poissonfocus.h"

enum pfs_errors
{
	PFS_NO_ERRORS = 0,
	PFS_ERROR_INVALID_ALLOCATION,
	PFS_ERROR_INVALID_INPUT,
};

typedef struct pfs PoissonFocusSES;

/**
 * This is the API for online mode applications.
 * For an example of usage, see implementation of pf_interface, in `focus.c`.
 */

PoissonFocusSES* pfs_init(enum pfs_errors* e, double threshold_std,
		double mu_min, double alpha, int m, int sleep);

void pfs_terminate(PoissonFocusSES* f);

enum pfs_errors pfs_step(PoissonFocusSES* f, bool* trigflag, count_t x_t);

struct pf_change pfs_get_change(PoissonFocusSES* f);

/**
 * This is the API for offline mode applications
 */

DLL00_EXPORT_API enum pfs_errors
pfs_interface(struct pf_changepoint* cp, count_t* xs, size_t len,
		double threshold_std, double mu_min, double alpha, int m, int sleep);

/**
 * Utilities.
 */

DLL00_EXPORT_API enum pfs_errors pfs_check_init_parameters(double threshold_std, double mu_min,
		double alpha, int m, int sleep);

#endif //POISSONFOCUSSES_H
