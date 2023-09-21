#include <math.h>
#include <stdio.h>

#define PRINT_LEVEL 0
#define	MAX_ORDER 16
#define MAXPOW 32		
#define	SMALL_ENOUGH 1.0e-18

/*
 * structure type for representing a polynomial
 */
typedef  	struct	p {
  int	ord;
  double	coef[MAX_ORDER+1];
} poly;

double RELERROR;
int MAXIT, MAX_ITER_SECANT;

/* Forward declaration */
void sbisect(int np, poly* sseq, double min, double max, int atmin, int atmax, double* roots, int roots_offset);
int modrf(int ord, double *coef, double	a, double b, double *val);
int numroots(int np, poly* sseq, int* atneg, int* atpos);
int numchanges(int np, poly* sseq, double a);
int buildsturm(int ord, poly* sseq);
double evalpoly (int ord, double* coef, double x);
double hyper_tan(double a, double x);

/* set termination criteria for polynomial solver */
void initialize_sturm_(double *tol_secant, int *max_iter_sturm, int *max_iter_secant)
{
  RELERROR = *tol_secant;
  MAXIT = *max_iter_sturm;
  MAX_ITER_SECANT = *max_iter_secant;
}

void solve_sturm_(int *p_order, int *n_root, double *poly_coeffs, double *roots)
{
  poly sseq[MAX_ORDER*2];
  double min, max, dt1;
  int order, i, j, nroots, nchanges, np, atmin, atmax;

  order = *p_order;

  for (i = order; i >= 0; i--) 
    {
      sseq[0].coef[i] = poly_coeffs[i];
    }

  if (PRINT_LEVEL > 0) 
    {
      for (i = order; i >= 0; i--) 
	{
	  printf("coefficients in Sturm solver\n");
	  printf("%d %lf\n", i, sseq[0].coef[i]);
	}
    }

  /*
   * build the Sturm sequence
   */
  np = buildsturm(order, sseq);
  
  if (PRINT_LEVEL > 0) 
    {
      printf("Sturm sequence for:\n");
      for (i = order; i >= 0; i--)
	printf("%lf ", sseq[0].coef[i]);
      printf("\n\n");
      for (i = 0; i <= np; i++) 
	{
	  for (j = sseq[i].ord; j >= 0; j--)
	    printf("%lf ", sseq[i].coef[j]);
	  printf("\n");
	}
      printf("\n");
    }

  /* 
   * get the number of real roots
   */

  nroots = numroots(np, sseq, &atmin, &atmax);


  if (nroots == 0) 
    {
      // printf("solve: no real roots\n");
      *n_root = nroots;
      return;
    }

  if (PRINT_LEVEL > 0) 
    printf("Number of real roots: %d\n", nroots);
 
  /*
   * calculate the bracket that the roots live in
   */
  min = -1.0;
  nchanges = numchanges(np, sseq, min);

  for (i = 0; nchanges != atmin && i != MAXPOW; i++) { 
    min *= 10.0;
    nchanges = numchanges(np, sseq, min);
  }

  if (nchanges != atmin) {
    printf("solve: unable to bracket all negative roots\n");
    atmin = nchanges;
  }

  max = 1.0;
  nchanges = numchanges(np, sseq, max);
  for (i = 0; nchanges != atmax && i != MAXPOW; i++) { 
    max *= 10.0;
    nchanges = numchanges(np, sseq, max);
  }
  
  if (nchanges != atmax) {
    printf("solve: unable to bracket all positive roots\n");
    atmax = nchanges;
  }

  nroots = atmin - atmax;

  /*
   * perform the bisection.
   */

  sbisect(np, sseq, min, max, atmin, atmax, roots, 0);

  *n_root = nroots;

  /*
   * write out the roots...
   */
  if (PRINT_LEVEL > 0) 
    {
      if (nroots == 1) {
	printf("\n1 distinct real root at x = %f\n", roots[0]);
      } else {
	printf("\n%d distinct real roots for x: \n", nroots);
	
	for (i = 0; i != nroots; i++)
	  {
	    printf("%f\n", roots[i]);
	  }
      }
    }

  return;
}

/*
 * modp
 *
 *	calculates the modulus of u(x) / v(x) leaving it in r, it
 *  returns 0 if r(x) is a constant.
 *  note: this function assumes the leading coefficient of v 
 *	is 1 or -1
 */

double hyper_tan(double a, double x)
{
  double exp_x1, exp_x2, ax;

  ax = a*x;
  if (ax > 100.0)
    {
      return(1.0);
    }
  else if (ax < -100.0)
    {
      return(-1.0);
    }
  else    
    {
      exp_x1 = exp(ax);
      exp_x2 = exp(-ax);
      return (exp_x1 - exp_x2)/(exp_x1 + exp_x2);
    }
}

static int modp(poly *u, poly *v, poly *r)
{
	int		k, j;
	double	*nr, *end, *uc;

	nr = r->coef;
	end = &u->coef[u->ord];

	uc = u->coef;
	while (uc <= end)
			*nr++ = *uc++;

	if (v->coef[v->ord] < 0.0) {


			for (k = u->ord - v->ord - 1; k >= 0; k -= 2)
				r->coef[k] = -r->coef[k];

			for (k = u->ord - v->ord; k >= 0; k--)
				for (j = v->ord + k - 1; j >= k; j--)
					r->coef[j] = -r->coef[j] - r->coef[v->ord + k]
					* v->coef[j - k];
	} else {
			for (k = u->ord - v->ord; k >= 0; k--)
				for (j = v->ord + k - 1; j >= k; j--)
				r->coef[j] -= r->coef[v->ord + k] * v->coef[j - k];
	}

	k = v->ord - 1;
	while (k >= 0 && fabs(r->coef[k]) < SMALL_ENOUGH) {
		r->coef[k] = 0.0;
		k--;
	}

	r->ord = (k < 0) ? 0 : k;

	return(r->ord);
}

/*
 * buildsturm
 *
 *	build up a sturm sequence for a polynomial in smat, returning
 * the number of polynomials in the sequence
 */
int buildsturm(int ord, poly* sseq)
{
	int		i;
	double	f, *fp, *fc;
	poly	*sp;

	sseq[0].ord = ord;
	sseq[1].ord = ord - 1;

	/*
	 * calculate the derivative and normalise the leading
	 * coefficient.
	 */
	f = fabs(sseq[0].coef[ord]*ord);

	fp = sseq[1].coef;
	fc = sseq[0].coef + 1;

	for (i = 1; i <= ord; i++)
	  {
	    *fp++ = *fc++ * i / f;
	  }


	/*
	 * construct the rest of the Sturm sequence
	 */
	//	for (sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++) {
	for (sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++) {

	  /*
	   * reverse the sign and normalise
	   */

	  f = -fabs(sp->coef[sp->ord]);
	  for (fp = &sp->coef[sp->ord]; fp >= sp->coef; fp--)
		  *fp /= f;
	}

	sp->coef[0] = -sp->coef[0];	/* reverse the sign */

	return(sp - sseq);
}

/*
 * numroots
 *
 *	return the number of distinct real roots of the polynomial
 * described in sseq.
 */
int numroots(int np, poly* sseq, int* atneg, int* atpos)
{
		int		atposinf, atneginf;
		poly	*s;
		double	f, lf;

		atposinf = atneginf = 0;


	/*
	 * changes at positive infinity
	 */
	lf = sseq[0].coef[sseq[0].ord];

	for (s = sseq + 1; s <= sseq + np; s++) {
			f = s->coef[s->ord];
			if (lf == 0.0 || lf * f < 0)
				atposinf++;
		lf = f;
	}

	/*
	 * changes at negative infinity
	 */
	if (sseq[0].ord & 1)
			lf = -sseq[0].coef[sseq[0].ord];
	else
			lf = sseq[0].coef[sseq[0].ord];

	for (s = sseq + 1; s <= sseq + np; s++) {
			if (s->ord & 1)
				f = -s->coef[s->ord];
			else
				f = s->coef[s->ord];
			if (lf == 0.0 || lf * f < 0)
				atneginf++;
			lf = f;
	}

	*atneg = atneginf;
	*atpos = atposinf;
	//	printf("atneginf, atposinf = %d %d\n", atneginf, atposinf);
	return(atneginf - atposinf);
}

/*
 * numchanges
 *
 *	return the number of sign changes in the Sturm sequence in
 * sseq at the value a.
 */
int numchanges(int np, poly* sseq, double a)
{
	int		changes;
	double	f, lf;
	poly	*s;

	changes = 0;

	lf = evalpoly(sseq[0].ord, sseq[0].coef, a);

	for (s = sseq + 1; s <= sseq + np; s++) {
			f = evalpoly(s->ord, s->coef, a);
			if (lf == 0.0 || lf * f < 0)
				changes++;
			lf = f;
//			printf("lf %lf %d \n", f, changes);
	}

	//	printf("%d \n", changes);
	return(changes);
}

/*
 * sbisect
 *
 *	uses a bisection based on the sturm sequence for the polynomial
 * described in sseq to isolate intervals in which roots occur,
 * the roots are returned in the roots array in order of magnitude.
 */
void sbisect(int np, poly* sseq, double min, double max, int atmin, int atmax, double* roots, int roots_offset)
{
  double	mid;
  int		n1 = 0, n2 = 0, its, atmid, nroot;

  if ((roots_offset >= 16) || (roots_offset < 0)) // Yep, this happens...
    return;

  if ((nroot = atmin - atmax) == 1) 
    {

      /*
       * first try a less expensive technique.
       */
      //      printf("min max %lf, %lf \n", min, max);
      if (modrf(sseq->ord, sseq->coef, min, max, &roots[roots_offset]))
	return;

      //      printf("try hard way\n");
      /*
       * if we get here we have to evaluate the root the hard
       * way by using the Sturm sequence.
       */
      for (its = 0; its < MAXIT; its++) 
	{
	  mid = (min + max) / 2;
	  
	  atmid = numchanges(np, sseq, mid);
	  
	  if (fabs(mid) > RELERROR) 
	    {
	      if (fabs((max - min) / mid) < RELERROR) 
		{
		  roots[roots_offset] = mid;
		  return;
		}
	    } 
	  else if (fabs(max - min) < RELERROR) 
	    {
	      roots[roots_offset] = mid;
	      return;
	    }
	  
	  if ((atmin - atmid) == 0)
	    min = mid;
	  else
	    max = mid;
	}
      
      if (its == MAXIT) 
	{
	  fprintf(stderr, "sbisect: overflow min %f max %f\
					diff %e nroot %d n1 %d n2 %d\n",
		  min, max, max - min, nroot, n1, n2);
	  roots[roots_offset] = mid;
	}
      return;
    }

  /*
   * more than one root in the interval, we have to bisect...
   */

  for (its = 0; its < MAXIT; its++) 
    {

      mid = (min + max) / 2;
      
      atmid = numchanges(np, sseq, mid);

      n1 = atmin - atmid;
      n2 = atmid - atmax;

      if (n1 != 0 && n2 != 0) 
	{
	  sbisect(np, sseq, min, mid, atmin, atmid, roots, roots_offset);
	  sbisect(np, sseq, mid, max, atmid, atmax, roots, roots_offset + n1);
	  break;
	}
      
      if (n1 == 0)
	min = mid;
      else
	max = mid;
    }

  if (its == MAXIT) 
    {
      /*
      fprintf(stderr, "sbisect: roots too close together\n");
      fprintf(stderr, "sbisect: overflow min %f max %f diff %e\
				nroot %d n1 %d n2 %d\n",
	      min, max, max - min, nroot, n1, n2);
      */
      for (n1 = atmax; n1 < atmin; n1++) {
		const int scary_idx = roots_offset + n1 - atmax;
		if ((scary_idx >= 0) && (scary_idx < 16))
			roots[scary_idx] = mid;
	  }
    }
}

/*
 * evalpoly
 *
 *	evaluate polynomial defined in coef returning its value.
 */
double evalpoly (int ord, double* coef, double x)
{
	double	*fp, f;

	fp = &coef[ord];
	f = *fp;

	for (fp--; fp >= coef; fp--)
	f = x * f + *fp;

	return(f);
}


/*
 * modrf
 *
 *	uses the modified regula-falsi method to evaluate the root
 * in interval [a,b] of the polynomial described in coef. The
 * root is returned is returned in *val. The routine returns zero
 * if it can't converge.
 */
int modrf(int ord, double *coef, double	a, double b, double *val)
{
  int its;
  double fa, fb, x, fx, lfx;
  double *fp, *scoef, *ecoef;

  scoef = coef;
  ecoef = &coef[ord];

  fb = fa = *ecoef;
  for (fp = ecoef - 1; fp >= scoef; fp--) {
    fa = a * fa + *fp;
    fb = b * fb + *fp;
  }

  /*
   * if there is no sign difference the method won't work
   */
  if (fa * fb > 0.0)
    return(0);

  /*  commented out to avoid duplicate solutions when the bounds are close to the roots

  if (fabs(fa) < RELERROR) 
    {
      *val = a;
      return(1);
    }
  
  if (fabs(fb) < RELERROR) 
    {
      *val = b;
      return(1);
    }
  */

  lfx = fa;

  for (its = 0; its < MAX_ITER_SECANT; its++) 
    {
      x = (fb * a - fa * b) / (fb - fa);

      // constrain that x stays in the bounds
      if (x < a || x > b)
	x = 0.5 * (a+b);

      fx = *ecoef;
      for (fp = ecoef - 1; fp >= scoef; fp--)
	fx = x * fx + *fp;

      if (fabs(x) > RELERROR) 
	{
	  if (fabs(fx / x) < RELERROR) 
	    {
	      *val = x;
	      //	      printf(" x, fx %lf %lf\n", x, fx);
	      return(1);
	    }
	} 
      else if (fabs(fx) < RELERROR) 
	{
	  *val = x;
	  //	  printf(" x, fx %lf %lf\n", x, fx);
	  return(1);
	}

      if ((fa * fx) < 0) 
	{
	  b = x;
	  fb = fx;
	  if ((lfx * fx) > 0)
	    fa /= 2;
	} 
      else 
	{
	  a = x;
	  fa = fx;
	  if ((lfx * fx) > 0)
	    fb /= 2;
	}

      lfx = fx;
    }

  //fprintf(stderr, "modrf overflow %f %f %f\n", a, b, fx);

  return(0);
}

