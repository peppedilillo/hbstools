#ifndef POISSONFOCUS_H
#define POISSONFOCUS_H

#if defined(_WIN32)
#  define DLL00_EXPORT_API __declspec(dllexport)
#else
#  define DLL00_EXPORT_API
#endif

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define PF_MAXCURVES 64

typedef int64_t count_t;

enum pf_errors
{
	PF_NO_ERRORS = 0,
	PF_ERROR_INVALID_ALLOCATION,
	PF_ERROR_INVALID_INPUT,
};

typedef struct pf PoissonFocus;

/**
 * This is the API for online mode applications.
 * For an example of usage, see implementation of pf_interface, in `focus.c`.
 */

struct pf_change
{
	double significance_std;
	int offset;
};

PoissonFocus* pf_init(enum pf_errors* err, double threshold_std, double mu_min);

void pf_terminate(PoissonFocus* f);

enum pf_errors pf_step(PoissonFocus* f, bool* trigflag, count_t x, double b);

struct pf_change pf_get_change(PoissonFocus* f);

/**
 * This is the API for offline mode applications
 */

struct pf_changepoint
{
	double significance_std;
	size_t changepoint;
	size_t triggertime;
};

DLL00_EXPORT_API enum pf_errors
pf_interface(struct pf_changepoint* cp, count_t* xs, double* bs, size_t len,
		double threshold, double mu_min);

/**
 * Utilities.
 */

struct pf_changepoint
pf_change2changepoint(struct pf_change c, size_t t);

DLL00_EXPORT_API enum pf_errors pf_check_init_parameters(double threshold_std, double mu_min);

#endif //POISSONFOCUS_H