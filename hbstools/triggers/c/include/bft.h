#ifndef BFT_H
#define BFT_H

#if defined(_WIN32)
#  define DLL00_EXPORT_API __declspec(dllexport)
#else
#  define DLL00_EXPORT_API
#endif

#include "poissonfocusses.h"

#define DETECTORS_NUMBER 4

enum bft_errors
{
	BFT_NO_ERRORS = 0,
	BFT_ERROR_INVALID_ALLOCATION,
	BFT_ERROR_INVALID_INPUT
};

typedef struct bft Bft;

/**
 * This is the API for online mode applications.
 * For an example of usage, see implementation of pf_interface, in `focus.c`.
 */

struct bft_changes
{
	struct pf_change changes[DETECTORS_NUMBER];
};

Bft* bft_init(enum bft_errors* err, double threshold_std, double mu_min,
		double alpha, int m, int sleep, int majority);

void bft_terminate(Bft* bft);

enum bft_errors bft_step(Bft* bft, bool* t, count_t* xs);

struct bft_changes bft_get_changes(Bft* bft);

/**
 * This is the API for offline mode applications
 */

struct bft_changepoints
{
	struct pf_changepoint changepoints[DETECTORS_NUMBER];
};

DLL00_EXPORT_API enum bft_errors bft_interface(struct bft_changepoints* cps,
		const count_t* xss, size_t len,
		double threshold_std, double mu_min,
		double alpha, int m, int sleep,
		int majority);

/**
 * Utilities
 */

DLL00_EXPORT_API enum bft_errors bft_check_init_parameters(double threshold_std,
		double mu_min, double alpha, int m, int sleep, int majority);

#endif //BFT_H