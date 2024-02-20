#include <assert.h>
#include "poissonfocusses.h"
#include "poissonfocus.h"

/**
 * A queue of int implemented over a circular buffer.
 * Most recent counts are parked here for while, until they are used
 * to update the background estimate. Then they are discarded.
 * This prevents pollution of background estimate with source photons.
 * Since we'll (..should..) never enqueue onto full queue,
 * or dequeue from an empty one, there are no guard-rails against these.
 */
struct queue
{
	int head;
	int tail;
	int len;
	count_t* arr;
};

/**
 * Wrapped helper to easen conversion to a static memory implementation.
 */
static struct queue*
queue_init_helper(struct queue* q, count_t* buffer, int len)
{
	q->head = 0;
	q->tail = 0;
	q->len = len;
	q->arr = buffer;
	return q;
}

static struct queue* queue_init(int len)
{
	struct queue* q = malloc(sizeof(struct queue));
	if (q == NULL)
		return NULL;

	count_t* buffer = malloc(sizeof(count_t) * (len + 1));
	if (buffer == NULL)
	{
		free(q);
		return NULL;
	}

	return queue_init_helper(q, buffer, len);
}

static void queue_terminate(struct queue* q)
{
	if (q != NULL)
		free(q->arr);
	free(q);
}

static bool queue_empty(struct queue* q)
{
	return q->head == q->tail;
}

static bool queue_full(struct queue* q)
{
	return (q->tail + 1) % (q->len + 1) == q->head;
}

static void queue_enqueue(struct queue* q, count_t n)
{
	assert(!queue_full(q));

	*(q->arr + q->tail++) = n;
	if (q->tail > q->len)
		q->tail = 0;
}

static count_t queue_dequeue(struct queue* q)
{
	assert(!queue_empty(q));

	count_t t = *(q->arr + q->head++);
	if (q->head > q->len)
		q->head = 0;
	return t;
}

/**
 * Now the actual implementation of focus with automatic background estimate
 * through single exponential smoothin.
 */
enum status_codes
{
	COLLECT = 0,
	UPDATE,
	TEST,
	STOP,
};

struct status
{
	enum status_codes code;
	enum pfs_errors latest_error;
};

struct pfs
{
	PoissonFocus* focus;
	struct queue* queue;
	double alpha;
	int m;
	int sleep;
	int t;
	double lambda_t;
	struct status status;
};

/**
 * Defines and checks the domain of the init function arguments.
 * Will return an error if the arguments are invalid.
 */
enum pfs_errors pfs_check_init_parameters(double threshold_std, double mu_min,
	double alpha, int m, int sleep)
{
	// @formatter:off
  if (
	  pf_check_init_parameters(threshold_std, mu_min)
	     == PF_ERROR_INVALID_INPUT ||
	              alpha     <  0.0 ||
				  alpha     >  1.0 ||
				  m         <  1   ||
				  sleep     <  0
	  )
    // @formatter:on
	  return PFS_ERROR_INVALID_INPUT;
	return PFS_NO_ERRORS;
}

/**
 * Dynamically allocate an algorithm with its own data structures.
 * If any error is encountered at this point, stuff allocated prior to the
 * error are freed, then we return NULL.
 * Wraps init_helper.
 */
static struct pfs* init_helper(struct pfs* f,
	struct pf* focus,
	struct queue* q, double alpha,
	int m, int sleep)
{
	f->focus = focus;
	f->queue = q;
	f->alpha = alpha;
	f->m = m;
	f->sleep = sleep;
	f->t = sleep + m;
	f->status = (struct status){
		COLLECT,
		PFS_NO_ERRORS
	};
	return f;
}

/**
 * Initializes FOCuS with automatic background estimate by exponential smoothing.
 *
 * @param e : a pointer to an error.
 * @param threshold_std : threshold value in standard deviations
 * @param mu_min : will keep memory usage low and constant at cost of losing
 * older chagepoint. See FOCuS Ward 2023 or Dilillo 2024 for more info.
 * @param alpha : simple exponential smoothing factor.
 * @param m : only counts gathered up to m time-steps in the past are used
 * for assessing background
 * @param sleep : testing for anomalies starts only after `sleep + m` time-steps.
 * @return a pointer
 */
struct pfs* pfs_init(enum pfs_errors* e, double threshold_std,
	double mu_min, double alpha, int m, int sleep)
{
	if (pfs_check_init_parameters(threshold_std, mu_min, alpha, m, sleep))
	{
		*e = PFS_ERROR_INVALID_INPUT;
		return NULL;
	}

	struct pfs* f = malloc(sizeof(struct pfs));
	if (f == NULL)
	{
		*e = PFS_ERROR_INVALID_ALLOCATION;
		return NULL;
	}

	struct queue* q = queue_init(m);
	if (q == NULL)
	{
		free(f);
		*e = PFS_ERROR_INVALID_ALLOCATION;
		return NULL;
	}

	enum pf_errors status;
	PoissonFocus* focus = pf_init(&status, threshold_std, mu_min);
	if (status == PF_ERROR_INVALID_ALLOCATION)
	{
		free(f);
		queue_terminate(q);
		*e = PFS_ERROR_INVALID_ALLOCATION;
		return NULL;
	}

	*e = PFS_NO_ERRORS;
	return init_helper(f, focus, q, alpha, m, sleep);
}

/**
 * Terminates freeing memory.
 */
void pfs_terminate(struct pfs* f)
{
	if (f != NULL)
	{
		queue_terminate(f->queue);
		pf_terminate(f->focus);
	}
	free(f);
}

/**
 * Automatically sets initial value for exponential smoothing.
 * Smoothed observation is set to the mean of the counts in the queue buffer.
 * Intended to be used over a full queue.
 * Trend initial value is always set to 0.0.
 * See: https://www.itl.nist.gov/div898/handbook/pmc/section4/pmc433.htm
 */
static void set_initial_bkg(struct pfs* f)
{
	struct queue* q = f->queue;
	assert(queue_full(q));

	count_t total = 0;
	int i = q->head;
	while (i != q->tail)
	{
		total += *(q->arr + i++);
		if (i > q->len)
			i = 0;
	}

	f->lambda_t = (double)total / f->m;
}

/**
 * Updates background estimate through exponential smoothing.
 * Returns a background estimate based on past observations.
 */
static double update_bkg(struct pfs* f, count_t x_t)
{
	double old_lambda = f->lambda_t;
	double lambda_t = f->alpha * x_t + (1 - f->alpha) * old_lambda;
	return lambda_t;
}

/**
 * A quality control function.
 * To have a trigger we require a changepoint earlier than m.
 * When this happens, we return true.
 */
static inline bool triggered(struct pfs* f, const bool* focus_triggered)
{
	if (*focus_triggered)
	{
		struct pf_change c = pf_get_change(f->focus);
		if (c.offset < f->m)
			return true;
	}
	return false;
}

/**
 * Gets the oldest count out of queue, and use it to update the background
 * estimate. Then, feed the background estimate and the most recent count to
 * focus. Will return an error if focus is fed a negative background or count.
 */
static inline enum pf_errors
step_test(struct pfs* f, bool* trigflag, count_t x_t)
{
	// updates background and queue
	count_t x_t_m = queue_dequeue(f->queue);
	f->lambda_t = update_bkg(f, x_t_m);
	queue_enqueue(f->queue, x_t);

	// feeds focus and checks if we got triggers
	bool focus_triggered;
	enum pf_errors err;
	err = pf_step(f->focus, &focus_triggered, x_t, f->lambda_t);
	*trigflag = triggered(f, &focus_triggered);
	return err;
}

/**
 * Gets oldest count out of queue and use it to update the background estimate.
 */
static inline void step_update(struct pfs* f, count_t x_t)
{
	count_t x_t_m = queue_dequeue(f->queue);
	f->lambda_t = update_bkg(f, x_t_m);
	queue_enqueue(f->queue, x_t);
}

/**
 * The base step. Will return an error if background estimate is not positive.
 * After using this you are supposed to check the bool flag to see
 * if you got a trigger.
 *
 * @param f : a pointer to an initialized structure.
 * @param trigflag : this is a bool flag, we will write here if we found a trigger.
 * @param x_t : latest count.
 * @return : error code
 */
enum pfs_errors pfs_step(PoissonFocusSES* f, bool* trigflag, count_t x_t)
{
	switch (f->status.code)
	{
	case TEST:
	{
		// update bkg and run focus test.
		if (step_test(f, trigflag, x_t) == PF_ERROR_INVALID_INPUT)
		{
			// focus stopped because of an invalid input. we stop too.
			f->status.code = STOP;
			f->status.latest_error = PFS_ERROR_INVALID_INPUT;
			return PFS_ERROR_INVALID_INPUT;
		}
		break;
	}
	case UPDATE:
	{
		// updates background estimate without testing.
		// not accessed if sleep == 0.
		step_update(f, x_t);
		// `t` acts as a countdown, when it gets to 0 we start operations.
		if (--f->t == 0)
			f->status.code = TEST;
		break;
	}
	case COLLECT:
	{
		// adds to the queue until full, then initialize background.
		queue_enqueue(f->queue, x_t);
		if (--f->t == f->sleep)
		{
			assert(queue_full(f->queue));
			set_initial_bkg(f);
			// if sleep is 0 we start testing at next iteration.
			f->status.code = f->sleep ? UPDATE : TEST;
		}
		break;
	}
	case STOP:
	{
		// an error happened. do nothing and return the error.
		return f->status.latest_error;
	}
	}
	return PFS_NO_ERRORS;
}

/**
 * The step function returns YES/NO information on wether a trigger happened.
 * If you need more information on how and when a trigger actually happened you
 * call `get_change`, which will return a `change` structure.
 * This structure bears the significance of the trigger in units of standard
 * deviations and its time offset as a step index.
 * If the algorithm did not resolve a trigger, will return (0.0, 0)
 */
struct pf_change pfs_get_change(struct pfs* f)
{
	struct pf_change c = pf_get_change(f->focus);
	if (c.offset < f->m)
		return c;
	return (struct pf_change){ 0 };
}

/**
 * Logs time-step, most recent count, background estimate and curve lists.
 */
void pfs_print(struct pfs* f, size_t t, count_t x_t)
{
	pf_print(f->focus, t, x_t, f->lambda_t);
}

/**
 * An interface example. This is intended to either be used for offline
 * applications, or to show how to use the functions of this library.
 *
 * @param cp : a pointer to a changepoint, where we store the results.
 * @param xs : an array of counts with length len.
 * @param len : length of input array.
 * @param threshold_std : threshold value in standard deviations
 * @param mu_min : will keep memory usage low and constant at cost of losing
 * older chagepoint. See Ward 2023, Dilillo 2024 for more info.
 * @param alpha : simple exponential smoothing factor.
 * @param m : only counts gathered up to m time-steps in the past are used
 * for assessing background
 * @param sleep : testing for anomalies starts only after `sleep + m` time-steps.
 * @return : an error code.
 */
enum pfs_errors
pfs_interface(struct pf_changepoint* cp, count_t* xs, size_t len,
	double threshold_std, double mu_min,
	double alpha, int m, int sleep)
{
	enum pfs_errors err = PFS_NO_ERRORS;
	PoissonFocusSES* focusexp = pfs_init(&err, threshold_std, mu_min, alpha, m, sleep);

	// inititalization can fail either because of wrong inputs,
	// or because failed allocation for curve stack.
	// we return error code, setting the changepoint to (0.0, 1, 0)
	if (err == PFS_ERROR_INVALID_ALLOCATION ||
		err == PFS_ERROR_INVALID_INPUT)
	{
		*cp = (struct pf_changepoint){ 0 };
		pfs_terminate(focusexp);
		return err;
	}

	bool got_trigger;
	size_t t;
	for (t = 0; t < len; t++)
	{
		err = pfs_step(focusexp, &got_trigger, xs[t]);

		if (err == PFS_ERROR_INVALID_INPUT)
		{
			// the focus_step fails if it is provided with a non-positive background.
			// if this happens at iteration `t`, we return immediately with
			// changepoint (0.0, t + 1, t).
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
		{
			// continue operations, log if needed.
			// if no trigger is found by the end of the time series,
			// return changepoint (0.0, len + 1, len).
			;
		}
	}

	*cp = pf_change2changepoint(pfs_get_change(focusexp),
		t == len ? t - 1 : t);
	pfs_terminate(focusexp);
	return err;
}
