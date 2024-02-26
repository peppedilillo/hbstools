#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "poissonfocus.h"

#define STACK_LEN (PF_MAXCURVES + 1)

#define COUNT_MAX INT64_MAX
typedef int64_t acc_count_t;
typedef int64_t acc_timestep_t;

struct curve
{
	acc_count_t x;
	double b;
	acc_timestep_t t;
	double m;
};

static struct curve NULL_CURVE;
static struct curve TAIL_CURVE;

/**
 * A stack implementation over a circular buffer, with curve elements.
 */
struct stack
{
	int head;
	int tail;
	int capacity;
	struct curve* arr;
};

static struct stack*
stack_init_helper(struct stack* s, struct curve* curve_buffer, int capacity)
{
	s->head = 0;
	s->tail = 0;
	s->capacity = capacity;
	s->arr = curve_buffer;
	return s;
}

static struct stack* stack_init()
{
	struct stack* s = malloc(sizeof(struct stack));
	if (s == NULL) return NULL;
	struct curve* curve_buffer = malloc(sizeof(struct curve) * STACK_LEN);
	if (curve_buffer == NULL)
	{
		free(s);
		return NULL;
	}
	return stack_init_helper(s, curve_buffer, PF_MAXCURVES);
}

static void stack_terminate(struct stack* s)
{
	if (s != NULL) free(s->arr);
	free(s);
}

static int stack_empty(struct stack* s)
{
	return s->head == s->tail;
}

static int stack_full(struct stack* s)
{
	return s->head == s->capacity ? s->tail == 0 : s->head + 1 == s->tail;
}

static void stack_push(struct stack* s, struct curve* q)
{
	if (stack_full(s))
	{
		// we just throw away the oldest curve.
		// TODO: Since we are using dynamic allocation anyway, this should be
		//       changed to call buffer resizing..
		s->tail == s->capacity ? s->tail = 0 : s->tail++;
		*(s->arr + s->tail) = TAIL_CURVE;
	}
	*(s->arr + s->head) = *q;
	s->head == s->capacity ? s->head = 0 : s->head++;
}

static struct curve* stack_pop(struct stack* s)
{
	assert(!stack_empty(s));
	s->head == 0 ? s->head = s->capacity : s->head--;
	return s->arr + s->head;
}

static struct curve* stack_peek(struct stack* s)
{
	assert(!stack_empty(s));
	return s->arr + (s->head == 0 ? s->capacity - 1 : s->head - 1);
}

static void stack_reset(struct stack* s)
{
	s->head = 0;
	s->tail = 0;
}

static double curve_max(struct curve* c, struct curve* acc)
{
	count_t x = (count_t)(acc->x - c->x);
	double b = (acc->b - c->b);
	assert (x > b);
	return x * log(x / b) - (x - b);
}

static int curve_dominate(struct curve* p, struct curve* q, struct curve* acc)
{
	double p_x = (double)(acc->x - p->x);
	double p_b = acc->b - p->b;
	double q_x = (double)(acc->x - q->x);
	double q_b = acc->b - q->b;
	if (p_x * q_b - q_x * p_b > 0)
		return +1;
	return -1;
}

/**
 * This is private. For the public equivalent see `focus_change`.
 * The difference between `change` and `focus_change` is that the latter
 * express significance in units of standard deviations, while the former
 * express it as loglikelihood ratio. Offsets are non-negative.
 */
struct change
{
	double significance_llr;
	int offset;
};

enum status_codes
{
	TEST = 0,
	STOP
};

struct status
{
	enum status_codes code;
	enum pf_errors latest_error;
};

struct pf
{
	struct status status;
	struct stack* curves;
	struct change change;
	double mu_crit;
	double threshold_llr;
};

/**
 * Defines and checks the domain of the init function arguments.
 * Returns an error if the arguments are invalid.
 */
enum pf_errors pf_check_init_parameters(double threshold_std, double mu_min)
{
	if (
		threshold_std <= 0.0 ||
		mu_min < 1.0
		)
		return PF_ERROR_INVALID_INPUT;
	return PF_NO_ERRORS;
}

/**
 * We have this helper function separate from the main interface.
 * The reason for this is to easen conversion of this program to
 * an implementation without dynamic allocation.
 */
static struct pf*
init_helper(struct pf* f, struct stack* curves, double threshold_std, double mu_min)
{
	f->status = (struct status){
		TEST,
		PF_NO_ERRORS
	};
	f->curves = curves;
	f->change = (struct change){ 0.0, 0 };
	f->threshold_llr = threshold_std * threshold_std / 2;
	f->mu_crit = (mu_min == 1. ? 1.0 : (mu_min - 1) / log(mu_min));

	NULL_CURVE = (struct curve){ 0 };
	TAIL_CURVE = (struct curve){ COUNT_MAX, 0., 0, 0. };
	stack_push(curves, &TAIL_CURVE);
	stack_push(curves, &NULL_CURVE);
	return f;
}

/**
 * Dynamically allocate a focus structure with its curve stack.
 * If any error is encountered at this point, stuff allocated prior to the
 * error are freed, then we return NULL.
 * Wraps init_helper.
 */
struct pf*
pf_init(enum pf_errors* err, double threshold_std, double mu_min)
{
	if (threshold_std <= 0.0 || mu_min < 1.0)
	{
		*err = PF_ERROR_INVALID_INPUT;
		return NULL;
	}
	struct pf* f = malloc(sizeof(struct pf));
	if (f == NULL)
	{
		*err = PF_ERROR_INVALID_ALLOCATION;
		return NULL;
	}
	struct stack* curves = stack_init();
	if (curves == NULL)
	{
		free(f);
		*err = PF_ERROR_INVALID_ALLOCATION;
		return NULL;
	}

	*err = PF_NO_ERRORS;
	init_helper(f, curves, threshold_std, mu_min);
	return f;
}

/**
 * Terminates focus and frees memory.
 */
void pf_terminate(struct pf* f)
{
	if (f != NULL) stack_terminate(f->curves);
	free(f);
}

/**
 * Fast FOCuS maximizer, implementing the strategy of the paper
 * https://doi.org/10.48550/arXiv.2302.04743, by Kes Ward.
 */
static void maximize(struct pf* f, struct curve* p, struct curve* acc)
{
	struct stack* curves = f->curves;
	double m = acc->m - p->m;
	int i = curves->head;
	while (m + p->m >= f->threshold_llr)
	{
		if (m >= f->threshold_llr)
		{
			f->change.significance_llr = m;
			f->change.offset = (int)(acc->t - p->t);
			break;
		}
		i == 0 ? i = curves->capacity : i--;
		p = (curves->arr + i);
		m = curve_max(p, acc);
	}
}

/**
 * Fast FOCuS updater, see Dilillo 2024.
 */
static struct change step_helper(struct pf* f, count_t x_t, double b_t)
{
	// b_t here is supposed to be greater than zero, and you are supposed
	// to check for that in wrapper.
	f->change = (struct change){ 0.0, 0 };

	struct stack* curves = f->curves;
	struct curve* p = stack_pop(curves);
	struct curve acc = { p->x + x_t, p->b + b_t, p->t + 1, p->m };
	while (curve_dominate(p, stack_peek(curves), &acc) < 0)
		p = stack_pop(curves);

	if ((double)(acc.x - p->x) > f->mu_crit * (acc.b - p->b))
	{
		double m = curve_max(p, &acc);
		acc.m = p->m + m;
		maximize(f, p, &acc);
		stack_push(curves, p);
		stack_push(curves, &acc);
	}
	else
	{
		stack_reset(curves);
		stack_push(curves, &TAIL_CURVE);
		stack_push(curves, &NULL_CURVE);
	}

	return f->change;
}

static inline bool triggered(struct pf* f)
{
	if (f->change.significance_llr > 0)
		return true;
	return false;
}

/**
 * A wrapper to the focus's update step.
 * The algorithm can be in one of two states:
 * 1. Running. Being in this state guarantees no error occurred in the past.
 *             If a valid background value is give, the algorithm updates,
 *             notifies if it got a trigger through the `got_trigger` flag
 *             and returns a no-error code. If an invalid background is given,
 *             the algorithm stops, the `got_trigger` is set to false, and
 *             an appropriate error is returned.
 * 2. Stopped. This happens if we successfully initialized the algorithm but
 *             ended up with an error during runtime: we do nothing and return
 *             the latest error code. The algorithm is _not_ stopped if it
 *             encounters a trigger.
 *
 * @param f : a pointer to an initialized structure.
 * @param trigflag : this is a bool flag, we will write to this if
 * we found any trigger during the update.
 * @param x : latest count.
 * @param b : latest background rate value.
 * @return : error code
 */
enum pf_errors
pf_step(struct pf* f, bool* trigflag, count_t x, double b)
{
    *trigflag = false;

    switch (f->status.code)
	{
	case TEST :
	{
		if (b <= 0 || x < 0)
		{
			f->status.code = STOP;
			f->status.latest_error = PF_ERROR_INVALID_INPUT;
			f->change = (struct change){ 0 };
			return PF_ERROR_INVALID_INPUT;
		}
		f->change = step_helper(f, x, b);
		*trigflag = triggered(f);
		break;
	}
	case STOP :
	{
		return f->status.latest_error;
	}
	}
	return PF_NO_ERRORS;
}

/**
 * The step function returns YES/NO information on wether a trigger happened.
 * If you need more information on how and when a trigger actually happened you
 * call `get_change`, which will return a `change` structure.
 * This structure bears the significance of the trigger in units of standard
 * deviations and its time offset as a non-negative step index.
 *
 * Within the present implementation a non-trivial change is returned only if
 * the algorithm triggered, otherwise the change will be (0.0, 0).
 */
struct pf_change pf_get_change(struct pf* f)
{
	return (struct pf_change){ sqrt(2 * f->change.significance_llr),
							   f->change.offset };
}

/**
 * Changes are what you get out of the algorithm when it runs online.
 * Changepoints are what you get out of the algorithms when it runs offline.
 * The difference between the two is in how the modes deal with time: changes will
 * report time as an offset index from present iteration; changepoints will
 * return the actual step index of the anomaly.
 * This implies that in online mode the user is responsible for keeping track
 * of the time passed, while in offline model that is on us.
 *
 * This is an utility function for converting between changes and changepoints.
 */
struct pf_changepoint
pf_change2changepoint(struct pf_change c, size_t t)
{
	return (struct pf_changepoint){
		c.significance_std,
		t - c.offset + 1,
		t
	};
}

/**
 * An interface example. This is intended to either be used for offline
 * applications, or to show how to use the functions of this library.
 *
 * @param cp : a changepoint point, here we will store the result
 * @param xs : an array of counts with length len.
 * @param bs : an array of background rates, as expected per bin-step,
 * with length len.
 * @param len : number of time-steps in the time series.
 * @param threshold : threshold, in units of standard deviation.
 * @param mu_min : will keep memory usage low and constant at cost of losing
 * older chagepoint. See Ward 2023, Dilillo 2024 for more info.
 * @return
 */
enum pf_errors
pf_interface(struct pf_changepoint* cp, count_t* xs, double* bs, size_t len, double threshold, double mu_min)
{
	enum pf_errors err = PF_NO_ERRORS;
	PoissonFocus* focus = pf_init(&err, threshold, mu_min);

	// inititalization can fail either because of wrong inputs,
	// or because failed allocation for curve stack.
	// we return error code, setting the changepoint to (0.0, 0, 0)
	if (err == PF_ERROR_INVALID_ALLOCATION
		|| err == PF_ERROR_INVALID_INPUT)
	{
		*cp = (struct pf_changepoint){ 0 };
		pf_terminate(focus);
		return err;
	}

	// loop over data and keep track of steps in `t`.
	bool got_trigger;
	size_t t;
	for (t = 0; t < len; t++)
	{
		err = pf_step(focus, &got_trigger, xs[t], bs[t]);

		if (err == PF_ERROR_INVALID_INPUT)
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

	*cp = pf_change2changepoint(pf_get_change(focus),
		t == len ? t - 1 : t);
	pf_terminate(focus);
	return err;
}
