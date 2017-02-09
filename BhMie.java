
import Complex;

public class BhMie
{

private Complex[] cxs1;//perpendicular electric field
private Complex[] cxs2;//parallel electric field
private double gsca;//assymetry factor
private double qsca;//scattering cross section
private double qext;//extinction cross section
private double qback;//backscattering cross section
private static Complex CXONE=new Complex(1.0, 0.0);


/**
 * Adaptation in java of the FORTRAN code (see https://github.com/hyperion-rt/bhmie)
 * itself derived from the Bohren-Huffman Mie scattering subroutine to calculate 
 * scattering and absorption by a homogenous isotropic sphere.
 * @param cxs1 perpendicular electric field
 * @param cxs2 parallel electric field
 * @param nang number of angles calculated
 * @param x reduced particle size: x=2*Math.PI*diameter/2/lambda*fluidIndex;
 * @param cxref complex index of refraction of the particle material
 */
public void  bhmie(Complex[] cxs1,Complex[] cxs2,int nang,double x,Complex cxref)
{
/* .. Array Arguments .. */
/* .. Local Scalars ..*/
Complex cxan=new Complex();
Complex cxbn=new Complex();
Complex cxan1=new Complex();
Complex cxbn1=new Complex();
Complex cxxi, cxxi0, cxy, cxxi1;
Complex cxtemp;
double apsi, apsi0, apsi1, chi, chi0, chi1, dang, fn, p, pii,rn, t, theta, xstop, ymod;
double  dn, dx, psi, psi0, psi1;
int  j, jj, n, nmx, nn, nstop;
/* .. Local Arrays ..*/
double[] amu=new double[nang+1];
double[] pi=new double[nang+1];
double[] pi0=new double[nang+1];
double[] pi1=new double[nang+1];
double[] tau=new double[nang+1];

pii = Math.PI;
dx = x;
cxy = Complex.mul(new Complex(x,0.0),cxref);

/* Series expansion terminated after NSTOP terms */
xstop = x + 4.E0*Math.pow(x,0.3333) + 2.0;
nstop = (int)xstop;
ymod = Complex.abs(cxy);
nmx = (int)(Math.max(xstop,ymod) + 15);

Complex[] cxd=new Complex[nmx+1];

dang = .5E0*pii/ (double)(nang-1);
for (j = 1; j<=nang; j++) 
	{
	theta = (double)(j-1)*dang;
	amu[j] = Math.cos(theta);
	}



/* Logarithmic derivative D(J) calculated by downward recurrence
    beginning with initial value (0.,0.) at J=NMX */

cxd[nmx] =new  Complex(0.E0,0.E0);
nn = nmx - 1;

for (n = 1; n<= nn; n++) 
	{
	rn = nmx - n + 1;
	cxtemp=Complex.add(cxd[nmx-n+1],Complex.div(new Complex(rn,0.0),cxy));
	cxtemp=Complex.div(CXONE,cxtemp);
	cxd[nmx-n]=Complex.sub(Complex.div(new Complex(rn,0.0),cxy),cxtemp);
	}

for ( j = 1; j <= nang; j++) 
	{
	pi0[j] = 0.E0;
	pi1[j] = 1.E0;
	}
nn = 2*nang - 1;
for(j = 1; j<= nn; j++) 
	{
	cxs1[j] = new Complex(0.E0,0.E0);
	cxs2[j] = new Complex(0.E0,0.E0);
	}



/* Riccati-Bessel functions with real argument X
    calculated by upward recurrence */

psi0 = Math.cos(dx);
psi1 = Math.sin(dx);
chi0 = -Math.sin(x);
chi1 = Math.cos(x);
apsi0 = psi0;
apsi1 = psi1;
cxxi0 = new Complex(apsi0,-chi0);
cxxi1 = new Complex(apsi1,-chi1);
qsca = 0.E0;
gsca = 0.E0;

for ( n = 1; n <= nstop; n++) 
	{  
	dn = n;
	rn = n;
	fn = (2.E0*rn+1.E0)/(rn*(rn+1.E0));
	psi = (2.E0*dn-1.E0)*psi1/dx - psi0;
	apsi = psi;
	chi = (2.E0*rn-1.E0)*chi1/x - chi0;
	cxxi = new Complex(apsi,-chi);
/* Store previous values of AN and BN for use
    in computation of g=<cos(theta)> */
	if (n>1) 
		{
		cxan1 = cxan;
		cxbn1 = cxbn;
		}

/* Compute AN and BN:*/

cxan=Complex.div(cxd[n],cxref);
cxan=Complex.add(cxan,new Complex(rn/x,0.0));
cxan=Complex.mul(cxan,new Complex(apsi,0.0));
cxan=Complex.sub(cxan,new Complex(apsi1,0.0));

cxtemp=Complex.div(cxd[n],cxref);
cxtemp=Complex.add(cxtemp,new Complex(rn/x,0.0));
cxtemp=Complex.mul(cxtemp,cxxi);
cxtemp=Complex.sub(cxtemp,cxxi1);
cxan=Complex.div(cxan,cxtemp);

cxbn=Complex.mul(cxref,cxd[n]);
cxbn=Complex.add(cxbn,new Complex(rn/x,0.0));
cxbn=Complex.mul(cxbn,new Complex(apsi,0.0));
cxbn=Complex.sub(cxbn,new Complex(apsi1,0.0));

cxtemp=Complex.mul(cxref,cxd[n]);
cxtemp=Complex.add(cxtemp,new Complex(rn/x,0.0));
cxtemp=Complex.mul(cxtemp,cxxi);
cxtemp=Complex.sub(cxtemp,cxxi1);
cxbn=Complex.div(cxbn,cxtemp);

/* Augment sums for *qsca and g=<cos(theta)> */
qsca = qsca + (2.*rn+1.)*(Complex.abs(cxan)*Complex.abs(cxan)+Complex.abs(cxbn)*Complex.abs(cxbn));
gsca = gsca + ((2.*rn+1.)/(rn*(rn+1.)))*(cxan.re()*cxbn.re()+cxan.im()*cxbn.im());

if (n>1) 
	{
	gsca = gsca + ((rn-1.)*(rn+1.)/rn)*(cxan1.re()*cxan.re()+
	cxan1.im()*cxan.im()+cxbn1.re()*cxbn.re()+cxbn1.im()*cxbn.im());
 	}

for ( j = 1; j<= nang; j++) 
	{
	jj = 2*nang - j;
	pi[j] = pi1[j];
	tau[j] = rn*amu[j]*pi[j] - (rn+1.E0)*pi0[j];
	p = Math.pow(-1.0,n-1);

	cxtemp=Complex.mul(cxan,new Complex(pi[j],0.0));
	cxtemp=Complex.add(cxtemp,Complex.mul(cxbn,new Complex(tau[j],0.0)));
	cxtemp=Complex.mul(new Complex(fn,0.0),cxtemp);
	cxs1[j]=Complex.add(cxs1[j],cxtemp);
	t = Math.pow(-1.0,n);

	cxtemp=Complex.mul(cxan,new Complex(tau[j],0.0));
	cxtemp=Complex.add(cxtemp,Complex.mul(cxbn,new Complex(pi[j],0.0)));
	cxtemp=Complex.mul(new Complex(fn,0.0),cxtemp);
	cxs2[j]=Complex.add(cxs2[j],cxtemp);

	if (j!=jj) 
		{
		cxtemp=Complex.mul(cxan,new Complex(pi[j]*p,0.0));
		cxtemp=Complex.add(cxtemp,Complex.mul(cxbn,new Complex(tau[j]*t,0.0)));
		cxtemp=Complex.mul(new Complex(fn,0.0),cxtemp);
		cxs1[jj]=Complex.add(cxs1[jj],cxtemp);

		cxtemp=Complex.mul(cxan,new Complex(tau[j]*t,0.0));
		cxtemp=Complex.add(cxtemp,Complex.mul(cxbn,new Complex(pi[j]*p,0.0)));
		cxtemp=Complex.mul(new Complex(fn,0.0),cxtemp);
		cxs2[jj]=Complex.add(cxs2[jj],cxtemp);
		}
	}

psi0 = psi1;
psi1 = psi;
apsi1 = psi1;
chi0 = chi1;
chi1 = chi;
cxxi1 = new Complex(apsi1,-chi1);

/*  For each angle J, compute pi_n+1
    from PI = pi_n , PI0 = pi_n-1 */

for ( j = 1; j<= nang; j++) 
	{
	pi1[j] = ((2.*rn+1.)*amu[j]*pi[j]-(rn+1.)*pi0[j])/rn;
	pi0[j] = pi[j];
	}
} /*end of big for */

/*  Have summed sufficient terms.
     Now compute *qsca,*qext,*qback,and *gsca */
gsca = 2.* gsca/ qsca;
qsca = (2.E0/(x*x))* qsca;
qext = (4.E0/(x*x))*cxs1[1].re();
qback = (4.E0/(x*x))*Complex.abs(cxs1[2*nang-1])*Complex.abs(cxs1[2*nang-1]);

return;
}


public double getQsca()
{
return qsca;
}


public double getQext()
{
return qext;
}


public double getGsca()
{
return gsca;
}


public double getQback()
{
return qback;
}



}
