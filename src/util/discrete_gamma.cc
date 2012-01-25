#include "util/discrete_gamma.h"
#include "util/logfile.h"
#include "util/math_fn.h"
#include "util/vector_output.h"

// All code in this file was taken from MrBayes (with superficial C->C++ conversion)
// by Ian Holmes on January 25, 2012

#define POINTGAMMA(prob,alpha,beta) 		PointChi2(prob,2.0*(alpha))/(2.0*(beta))

double  PointChi2 (double prob, double v);
double  LnGamma (double alp);
double  PointNormal (double prob);
double  IncompleteGamma (double x, double alpha, double LnGamma_alpha);

/* For comparing floating points: two values are the same if the absolute difference is less then 
   this value.
*/
#ifndef ETA
#define ETA (1E-30)
#endif

/*---------------------------------------------------------------------------------
|
|   DiscreteGamma
|
|   Discretization of gamma distribution with equal proportions in each
|   category.
|
---------------------------------------------------------------------------------*/
Discrete_gamma::Discrete_gamma (double alfa, double beta, int K, bool median)
  : vector<double> (K, (double) 0),
    alpha (alpha),
    beta (beta),
    K (K),
    median (median)
{
  vector<double>& rK (*this);

	int 			i;
	double 			gap05 = 1.0/(2.0*K), t, factor = alfa/beta*K, lnga1;

	if (median) 
		{
		for (i=0; i<K; i++) 
			rK[i] = POINTGAMMA((i*2.0+1.0)*gap05, alfa, beta);
		for (i=0,t=0; i<K; i++) 
			t += rK[i];
		for (i=0; i<K; i++)     
			rK[i] *= factor / t;
		}
	else 
		{
		lnga1 = LnGamma(alfa+1);
		/* calculate the points in the gamma distribution */
		for (i=0; i<K-1; i++) 
			rK[i] = POINTGAMMA((i+1.0)/K, alfa, beta);
		/* calculate the cumulative values */
		for (i=0; i<K-1; i++) 
			rK[i] = IncompleteGamma(rK[i] * beta, alfa + 1.0, lnga1);
		rK[K-1] = 1.0;
		/* calculate the relative values and rescale */
		for (i=K-1; i>0; i--)
			{
			rK[i] -= rK[i-1];
			rK[i] *= factor;
			}
		rK[0] *= factor;
		}
}

/*---------------------------------------------------------------------------------
|
|   PointChi2
|
|   Returns z so that Prob(x < z) = prob where x is Chi2 distributed with df=v. 
|   Returns -1 if in error.   0.000002 < prob < 0.999998.
|
---------------------------------------------------------------------------------*/
double PointChi2 (double prob, double v)

{

	double 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
					x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;

	if (p < 0.000002 || p > 0.999998 || v <= 0.0) 
		return (-1.0);
	g = LnGamma (v/2.0);
	xx = v/2.0;   
	c = xx - 1.0;
	if (v >= -1.24*log(p)) 
		goto l1;
	ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
	if (ch-e<0) 
		return (ch);
	goto l4;
	l1:
		if (v > 0.32) 
			goto l3;
		ch = 0.4;   
		a = log(1.0-p);
	l2:
		q = ch;  
		p1 = 1.0+ch*(4.67+ch);  
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = PointNormal (p);
		p1 = 0.222222/v;   
		ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v+6.0)  
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
		q = ch;   
		p1 = 0.5*ch;
		if ((t = IncompleteGamma (p1, xx, g)) < 0.0) 
			{
			CLOGERR << "   Error: Problem in PointChi2\n";
			return (-1.0);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));   
		b = t/ch;  
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a)) / 2520.0;
		s6 = (120.0+c*(346.0+127.0*c)) / 5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e) 
			goto l4;
		return (ch);

}


/*---------------------------------------------------------------------------------
|
|   LnGamma
|
|   Calculates the log of the gamma function. The Gamma function is equal
|   to:
|
|      Gamma(alp) = {integral from 0 to infinity} t^{alp-1} e^-t dt
|
|   The result is accurate to 10 decimal places. Stirling's formula is used
|   for the central polynomial part of the procedure.
|
|   Pike, M. C. and I. D. Hill.  1966.  Algorithm 291: Logarithm of the gamma
|      function.  Communications of the Association for Computing
|      Machinery, 9:684.
|      
---------------------------------------------------------------------------------*/
double LnGamma (double alp)

{

	double 		x = alp, f = 0.0, z;
	
	if (x < 7) 
		{
		f = 1.0;
		z = x-1.0;
		while (++z < 7.0)  
			f *= z;
		x = z;   
		f = -log(f);
		}
	z = 1.0 / (x*x);
	return  (f + (x-0.5)*log(x) - x + 0.918938533204673 + 
			(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
			0.083333333333333)/x);  

}


/*---------------------------------------------------------------------------------
|
|   PointNormal
|
|   Returns z so That Prob{x<z} = prob where x ~ N(0,1) and
|   (1e-12) < prob < 1-(1e-12).  Returns (-9999) if in error.
|
|   Odeh, R. E. and J. O. Evans.  1974.  The percentage points of the normal
|     distribution.  Applied Statistics, 22:96-97 (AS70).
|
|   Newer methods:
|
|   Wichura, M. J.  1988.  Algorithm AS 241: The percentage points of the
|      normal distribution.  37:477-484.
|   Beasley, JD & S. G. Springer.  1977.  Algorithm AS 111: The percentage
|      points of the normal distribution.  26:118-121.
|
---------------------------------------------------------------------------------*/
double PointNormal (double prob)

{

	double 		a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
 					a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
 					b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
 					y, z = 0, p = prob, p1;

	p1 = (p<0.5 ? p : 1-p);
	if (p1<1e-20) 
	   return (-9999);
	y = sqrt (log(1/(p1*p1)));   
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	
	return (p<0.5 ? -z : z);

}


/*---------------------------------------------------------------------------------
|
|   IncompleteGamma
|
|   Returns the incomplete gamma ratio I(x,alpha) where x is the upper
|   limit of the integration and alpha is the shape parameter.  Returns (-1)
|   if in error.   
|
|   Bhattacharjee, G. P.  1970.  The incomplete gamma integral.  Applied
|      Statistics, 19:285-287 (AS32)
|
---------------------------------------------------------------------------------*/
double IncompleteGamma (double x, double alpha, double LnGamma_alpha)

{

	int 			i;
	double 		p = alpha, g = LnGamma_alpha,
					accurate = 1e-8, overflow = 1e30,
					factor, gin = 0.0, rn = 0.0, a = 0.0, b = 0.0, an = 0.0, 
					dif = 0.0, term = 0.0, pn[6];

	if (fabs(x) < ETA) 
		return (0.0);
	if (x < 0 || p <= 0) 
		return (-1.0);

	factor = exp(p*log(x)-x-g);   
	if (x>1 && x>=p) 
		goto l30;
	gin = 1.0;  
	term = 1.0;  
	rn = p;
	l20:
		rn++;
		term *= x/rn;   
		gin += term;
		if (term > accurate) 
			goto l20;
		gin *= factor/p;
		goto l50;
	l30:
		a = 1.0-p;   
		b = a+x+1.0;  
		term = 0.0;
		pn[0] = 1.0;  
		pn[1] = x;  
		pn[2] = x+1;  
		pn[3] = x*b;
		gin = pn[2]/pn[3];
	l32:
		a++;  
		b += 2.0;  
		term++;   
		an = a*term;
		for (i=0; i<2; i++) 
			pn[i+4] = b*pn[i+2]-an*pn[i];
		if (fabs(pn[5]) < ETA) 
			goto l35;
		rn = pn[4]/pn[5];   
		dif = fabs(gin-rn);
		if (dif>accurate) 
			goto l34;
		if (dif<=accurate*rn) 
			goto l42;
	l34:
		gin = rn;
	l35:
		for (i=0; i<4; i++) 
			pn[i] = pn[i+2];
		if (fabs(pn[4]) < overflow) 
			goto l32;
		for (i=0; i<4; i++) 
			pn[i] /= overflow;
		goto l32;
	l42:
		gin = 1.0-factor*gin;
	l50:
		return (gin);

}
