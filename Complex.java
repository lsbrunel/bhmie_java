
// Author Laurent Brunel
package com.cionin.math;


/**
complex number
*/
public class Complex
{
private double re,im;

/**
@param re real part
@param im imaginary part
*/
public Complex(double re,double im)
{
this.re=re;
this.im=im;
}

public Complex()
{
this.re=0;
this.im=0;
}

public double re()
{
return re;
}

public double im()
{
return im;
}

public static double re(Complex z)
{
return z.re;
}

public double im(Complex z)
{
return z.im;
}

/**
return a new complex with this phase and norm
@param norm norm
@param phase phase
*/
public static Complex ComplexPolar(double norm,double phase)
{
return new Complex(norm*Math.cos(phase),norm*Math.sin(phase));
}


public double phase()
{
//if (re()==0) 
//	if (im()>=0) return Math.PI/2;
//	if (im()<0) return -Math.PI/2;
//return Math.atan(im()/re());
return Math.atan2(this.im, this.re);
}

public void setre(double r)
{
re=r;
}
public void setim(double r)
{
im=r;
}

/**return a new Complex addition of c1 and c2*/
public static Complex add(Complex c1,Complex c2)
{
Complex c=new Complex();
c.re=c1.re+c2.re;
c.im=c1.im+c2.im;
return c;
}

/**return a new Complex: c added to r*/
public static Complex add(Complex c,double r)
{
Complex c1=new Complex();
c1.re=c.re+r;
c1.im=c.im;
return c1;
}

/**return a new Complex added to r*/
public static Complex add(double r,Complex c)
{
return add(c,r);
}

/**return a new Complex added to r*/
public  Complex add(double r)
{
return add(this,r);
}


/**return a new Complex added to c*/
public  Complex add(Complex c)
{
return add(this,c);
}


public static Complex sub(Complex c1,Complex c2)
{
Complex c=new Complex();
c.re=c1.re-c2.re;
c.im=c1.im-c2.im;
return c;
}

public static Complex sub(double d1,Complex c2)
{
Complex c=new Complex();
c.re=d1-c2.re;
c.im=-c2.im;
return c;
}

public static Complex sub(Complex c1,double d2)
{
Complex c=new Complex();
c.re=c1.re-d2;
c.im=c1.im;
return c;
}

public  Complex sub(Complex c1)
{
return sub(this,c1);
}
public  Complex sub(double r)
{
return add(this,-r);
}


public static Complex mul(Complex c,double r)
{
Complex c1=new Complex();
c1.re=c.re*r;
c1.im=c.im*r;
return c1;
}

public static Complex mul(double r,Complex c)
{
return mul(c,r);
}

/**return a new Complex product of c1 and c2*/
public static Complex mul(Complex c1,Complex c2)
{
Complex c=new Complex();
c.re=c1.re*c2.re-c1.im*c2.im;
c.im=c1.re*c2.im+c2.re*c1.im;
return c;
}

/**return a new Complex product of c1 and c2*/
public static Complex mul(Complex c1,Complex c2,Complex c3)
{
return mul(c1,c2).mul(c3);
}
/**return a new Complex after multiplication by c*/
public Complex mul(Complex c)
{
return mul(this,c);
}
public Complex mul(double r)
{
return mul(this,r);
}


public static Complex div(Complex c,double r)
{
Complex c1=new Complex();
c1.re=c.re/r;
c1.im=c.im/r;
return c1;
}

public  Complex div(double r)
{
return div(this,r);
}

public  Complex div(Complex c)
{
return div(this,c);
}

public static Complex div(double r,Complex c)
{
return mul(inv(c),r);
}

public static Complex div(Complex c1,Complex c2)
{
return mul(c1,inv(c2));
}


public static Complex inv(Complex c)
{
double r=c.re*c.re+c.im*c.im;
Complex c1=new Complex();
c1.re=c.re/r;
c1.im=-c.im/r;
return c1;
}

public static Complex conj(Complex c)
{
Complex c1=new Complex();
c1.re=c.re;
c1.im=-c.im;
return c1;
}
public Complex conj()
{
return conj(this);
}



public static double norme(Complex c)
{
double r=Math.sqrt(c.re*c.re+c.im*c.im);
return r;
}

public static double abs(Complex c)
{
return norme(c);
}

public  double norme()
{
return norme(this);
}

public String toString()
{
String s="";
s+=re+" "+im;
return s;
}

/**return a new Complex copy of this one*/
public Complex copy() {
	
	return new Complex(this.re,this.im);
}


/**
 * complex sin function
 * sin(x+iy) = sin(x)cosh(y) + i cos(x)sinh(y)
 * @param z
 * @return
 */
public static Complex SIN(Complex z)
{
double x=z.re;
double y=z.im;
return new Complex(Math.sin(x)*Math.cosh(y),Math.cos(x)*Math.sinh(y));
}

/**
 * complex cos function
 * cos(x+iy) = cos(x)cosh(y) - i sin(x)sinh(y)
 * @param z
 * @return
 */
public static Complex COS(Complex z)
{
double x=z.re;
double y=z.im;
return new Complex(Math.cos(x)*Math.cosh(y),-Math.sin(x)*Math.sinh(y));
}


}
