#include "bft.h"
#include <string.h>

#define BITARRAY_BYTE_SIZE ((DETECTORS_NUMBER + 7) / 8)


/**
 * BFT stands for Big Focus Trigger :^).
 * It is a manager of multiple FOCuS algorithms, with autonoumous
 * background estimate, working independently.
 */
struct bft
{
	double threshold_std;
	double mu_min;
	double alpha;
	int m;
	int sleep;
	int majority;
	unsigned char bitarray[BITARRAY_BYTE_SIZE];
	PoissonFocusSES* fs[DETECTORS_NUMBER];
};

/**
 *  Sets on the n-th bit of the BFT's dead quadrant bit array.
 */
void bitarray_set(unsigned char* bitarray, int n)
{
	int index = n / 8;  // byte index
	int offset = n % 8;  // bit offset
	unsigned char mask = (unsigned char)(1 << offset); // little-endian

	bitarray[index] |= mask;
}

/**
 *  Counts the number of _zeros_ in BFT's bitarray, i.e. the number of
 *  working detectors/algorithms.
 */
int bitarray_count(const unsigned char* bitarray)
{
	int counts = DETECTORS_NUMBER;
	for (int i = 0; i < BITARRAY_BYTE_SIZE; i++)
	{ // iterate over bytes
		unsigned char byte = bitarray[i];
		int maxoffset = i < BITARRAY_BYTE_SIZE - 1 ? 8 : DETECTORS_NUMBER % 8;
		for (int offset = 0; offset < maxoffset; offset++)
		{
			counts -= byte & 1;
			byte >>= 1;
		}
	}
	return counts;
}

/**
 * Defines and checks the domain of the init function arguments.
 * Returns an error if the arguments are invalid.
 */
enum bft_errors bft_check_init_parameters(double threshold_std, double mu_min,
		double alpha, int m, int sleep, int majority)
{
	// @formatter:off
    if (
            pfs_check_init_parameters(
                    threshold_std, mu_min,
                    alpha, m, sleep)
            == PFS_ERROR_INVALID_INPUT ||
			              majority < 1 ||
                          majority > DETECTORS_NUMBER
            )
        return BFT_ERROR_INVALID_INPUT;
    // @formatter:on
	return BFT_NO_ERRORS;
}

/**
 * We have this helper function separate from the main interface.
 * The reason for this is to easen conversion of this program to
 * an implementation without dynamic allocation.
 */
static struct bft* init_helper(struct bft* bft, PoissonFocusSES** fs,
		double threshold_std, double mu_min,
		double alpha, int m, int sleep, int majority)
{
	bft->threshold_std = threshold_std;
	bft->mu_min = mu_min;
	bft->alpha = alpha;
	bft->m = m;
	bft->sleep = sleep;
	bft->majority = majority;
	memset(bft->bitarray, 0, sizeof(bft->bitarray));
	for (int i = 0; i < DETECTORS_NUMBER; i++)
		bft->fs[i] = fs[i];
	return bft;
}

/**
 * Dynamically each algorithms with their own data structures.
 * If any error is encountered at this point, stuff allocated prior to the
 * error are freed, then we return NULL.
 * Wraps init_helper.
 */
struct bft* bft_init(enum bft_errors* e, double threshold_std, double mu_min,
		double alpha, int m, int sleep, int majority)
{
	if (bft_check_init_parameters(threshold_std, mu_min, alpha, m, sleep, majority))
	{
		*e = BFT_ERROR_INVALID_INPUT;
		return NULL;
	}

	struct bft* bft = malloc(sizeof(struct bft));
	if (bft == NULL)
	{
		*e = BFT_ERROR_INVALID_ALLOCATION;
		return NULL;
	}

	PoissonFocusSES* fs[DETECTORS_NUMBER];
	enum pfs_errors focusexp_err;
	for (int i = 0; i < DETECTORS_NUMBER; i++)
	{
		fs[i] = pfs_init(&focusexp_err, threshold_std, mu_min, alpha, m, sleep);
		if (fs[i] == NULL)
		{
			free(bft);
			// terminate instances initialized before failure
			for (int j = i - 1; j >= 0; j--)
				pfs_terminate(fs[j]);
			*e = BFT_ERROR_INVALID_ALLOCATION;
			return NULL;
		}
	}

	*e = BFT_NO_ERRORS;
	return init_helper(bft, fs, threshold_std, mu_min, alpha, m, sleep, majority);
}

/**
 * Terminates freeing memory.
 */
void bft_terminate(struct bft* bft)
{
	if (bft != NULL)
		for (int i = 0; i < DETECTORS_NUMBER; i++)
			free(bft->fs[i]);
	free(bft);
}

/**
 * The base step. Will return an error if background estimate is not positive.
 * After using this you are supposed to check the bool flag to see
 * if you got a trigger.
 *
 * @param bft : a pointer to an initialized structure
 * @param got_trigger : this is a bool flag, we will write to this if
 * we found any trigger during the update
 * @param xs : this is a pointer to an array of `DETECTORS_NUMBER` ordered counts.
 * @return : error code
 */
enum bft_errors bft_step(Bft* bft, bool* got_trigger, count_t* xs)
{
	bool focusdes_triggered = false;
	int triggered_detectors = 0;
	enum pfs_errors err = PFS_NO_ERRORS;
	for (int i = 0; i < DETECTORS_NUMBER; i++)
	{
		err = pfs_step(bft->fs[i], &focusdes_triggered, xs[i]);
		if (err == PFS_ERROR_INVALID_INPUT) bitarray_set(bft->bitarray, i);
		else if (focusdes_triggered) triggered_detectors++;
	}
	*got_trigger = triggered_detectors >= bft->majority;

	return bitarray_count(bft->bitarray) < bft->majority ? BFT_ERROR_INVALID_INPUT
														 : BFT_NO_ERRORS;
}

/**
 * The step function returns YES/NO information on wether a trigger happened.
 * If you need more information on how and when a trigger actually happened you
 * call `get_change`, which will return a `changes` structure.
 * This structure bears changes for all detectors. See `pfs_get_change` for
 * more information on changes.
 *
 * Within the present implementation a non-trivial change is returned only if
 * the algorithm triggered, otherwise changes will be (0.0, 0).
 */
struct bft_changes bft_get_changes(Bft* bft)
{
	struct bft_changes r;
	for (int i = 0; i < DETECTORS_NUMBER; i++)
		r.changes[i] = pfs_get_change(bft->fs[i]);
	return r;
}

/**
 * Changes are what you get out of the algorithms when they run online.
 * Changepoints are what you get out of the algorithms when they run offline.
 * The difference between the two is in how the modes deal with time: changes will
 * report time as an offset index from present iteration; changepoints will
 * return the actual step index of the anomaly.
 * This implies that in online mode the user is responsible for keeping track
 * of the time passed, while in offline model that is on us.
 *
 * This is an utility function for converting between changes and changepoints.
 */
struct bft_changepoints bft_changes2changepoints(struct bft_changes c, size_t t)
{
	struct bft_changepoints cps;
	for (int i = 0; i < DETECTORS_NUMBER; i++)
	{
		cps.changepoints[i] = pf_change2changepoint(c.changes[i], t);
	}
	return cps;
}

/**
 * An interface example. This is intended to either be used for offline
 * applications, or to show how to use the functions of this library.
 * It supposes we are dealing with four detectors (DETECTORS_NUMBER == 4).
 *
 * @param cps : a pointer, here we will store the results
 * @param xss : an array of counts with shape (DETECTOR_NUMBER, len).
 * (t-th count, from detector i-th, is accessed with xss[i * len + t]).
 * @param len : length of input array.
 * @param threshold_std : threshold value in standard deviations
 * @param mu_min : will keep memory usage low and constant at cost of losing
 * older chagepoint. See Ward 2023, Dilillo 2024 for more info.
 * @param alpha : simple exponential smoothing factor.
 * @param m : only counts gathered up to m time-steps in the past are used
 * for assessing background
 * @param sleep : testing for anomalies starts only after `sleep + m` time-steps.
 * @param majority : the number of detectors simultaneously over threshold should
 * be equal or greater than this number for a trigger to pass through.
 * @return : an error code.
 */
enum bft_errors bft_interface(struct bft_changepoints* cps,
		const count_t* xss, size_t len,
		double threshold_std, double mu_min,
		double alpha, int m, int sleep,
		int majority)
{
	// inititalization can fail either because of wrong inputs,
	// or because failed allocation for curve stack.
	// we return error code, setting the changepoint to (0.0, 0, 0)
	enum bft_errors err = BFT_NO_ERRORS;
	Bft* bft = bft_init(&err, threshold_std, mu_min, alpha, m, sleep, majority);

	if (err == BFT_ERROR_INVALID_ALLOCATION ||
		err == BFT_ERROR_INVALID_INPUT)
	{
		*cps = (struct bft_changepoints){ 0 };
		bft_terminate(bft);
		return err;
	}

	bool got_trigger = false;
	size_t t;
	for (t = 0; t < len; t++)
	{
		count_t xs[DETECTORS_NUMBER];
		size_t i;
		for (i = 0; i < DETECTORS_NUMBER; i++) xs[i] = xss[i * len + t];

		err = bft_step(bft, &got_trigger, xs);

		if (err == BFT_ERROR_INVALID_INPUT)
		{
			// the focus_step fails if it is provided with a non-positive background.
			// if this happens at iteration `t` for any of the detectors,
			// we return immediately.
			break;
		}
		else if (got_trigger)
		{
			// do your trigger things, then return with changepoint:
			// (significance_std, changepoint, triggertime)
			// where changepoint < triggertime and significance_std > 0.
			break;
		}
		else
			// continue operations, log if needed.
			// if no trigger is found by the end of the time series,
			// return changepoint (0.0, len + 1, len).
			;
	}

	*cps = bft_changes2changepoints(bft_get_changes(bft), t == len ? t - 1 : t);
	bft_terminate(bft);
	return err;
}