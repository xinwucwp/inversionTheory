/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package hw6;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Nonlinear inversion with Gauss-Newton method.
 * <p>
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.04.26
 */
public class NonlinearInverter {

  /**
   * Sets dh for perturbing.
   * @param dh perturbation for computing Jacobian matrix.
   */
  public void setPerturb(double dh) {
    _dh = (float)dh;
  }
  /**
   * Sets forward.
   * @param xo the x-coordinates of observation locations.
   * @param zo the z-coordinates of observation locations.
   * @param xc central x-coordinates of dykes.
   * @param zt z-coordinates of the top of dykes.
   * @param wd width of each dyke, all dykes have the same width.
   * @param dc density contrast of each dyke, all dykes have the same dc.
   */
  public void setForward(float[] xo, float[] zo, 
    float[] xc, float[] zt, float wd, float dc) 
  {
    _xo = xo;
    _zo = zo;
    _xc = xc;
    _zt = zt;
    _wd = wd;
    _dc = dc;
  }

  /**
   * Sets scale factors for the smallest model and the derivative term 
   * in the model objective function.
   * @param alphaX scale factor for the smallest model.
   * @param alphaS scale factor for the derivative term.
   */
  public void setWmWd(double alphaX, double alphaS, double dSigma) {
    _alphaX = (float)alphaX;
    _alphaS = (float)alphaS;
    _dSigma = (float)dSigma;
  }
  /**
   * Sets a half-width for smoothing.
   * The smoothing serves as a preconditioner; they accelerate convergence
   * of the iterative solver used to compute model perturbation.
   * @param sigma half-width for smoothing, in samples.
   */
  public void setSmoothing(double sigma) {
    _sigma = (float)sigma;
  }
  
  /**
   * Sets parameters that control the number of CG solver iterations.
   * @param small stop iterations when error norm is reduced by this fraction.
   * @param niter maximum number of CG solver iterations.
   */
  public void setForCG(double small, int niter) {
    _small = (float)small;
    _niter = niter;
  }

  /**
   * Gets mappings computed from specified slopes and linearities.
   * @param sd sampling of data.
   * @param sm sampling of model.
   * @param d array of observed data.
   * @param m array of model to be updated.
   */
  public void inverter( 
    Sampling sm, float beta, float[] d, float[] m) 
  {
    // Sampling parameters.
    int nd = d.length;
    int nm = sm.getCount();
    float dm = (float)sm.getDelta();
    float[] dk = zerofloat(nd);
    float[] dp = zerofloat(nd);
    // Gauss-Newtom iterates
    boolean bool = true;
    forward(m,dk); // modeling data 
    float[] dok = sub(d,dk);// difference between observed and modeled data
    float dmc = dataMisfit(dok);// current data misfit
    float mof = modelNorm2(dm,m); 
    float tofc = dmc+beta*mof;
    while (bool){ 
      float[] b = new float[nm]; // right-hand side
      float[] p = new float[nm]; // updating direction
      VecArrayFloat1 vb = new VecArrayFloat1(b);
      VecArrayFloat1 vp = new VecArrayFloat1(p);
      Smoother smoother = new Smoother(nm,_sigma);
      A a = new A(smoother,dm,beta,dk,m);
      CgSolver cs = new CgSolver(_small,_niter);
      makeRhs(dm,beta,dk,dok,m,b);
      smoother.applyTranspose(b);
      cs.solve(a,vb,vp);
      smoother.apply(p);
      float alpha = lineSearch(dm,beta,d,p,m);
      add(m,mul(alpha,p),m);
      forward(m,dk);
      sub(d,dk,dok);
      float dmu = dataMisfit(dok);// updated data misfit
      float mou = modelNorm2(dm,m); 
      float tofu = dmu+beta*mou;
      System.out.println("tofc="+tofc); 
      System.out.println("tofu="+tofu); 
      if ((tofc-tofu)/tofc <=0.02f) {bool = false;}
      tofc = tofu;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _dh = 10.0f; 
  private float _alphaX = 1.0f; // coefficient of the derivative term for model
  private float _alphaS = 0.0f; //coefficient of smallest model 
  private float _dSigma = 0.0f; //standard deviation of observed data
  private float _sigma = 8.0f; // precon smoothing extent for 1st dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  // private for forward
  private static float[] _xo,_zo,_xc,_zt;
  private static float _wd,_dc;
  // Conjugate-gradient operators.
  private class A implements CgSolver.A {
    A(Smoother sf, float dm, float ba, float[] dk, float[] mk) 
    {
      _sf = sf;
      _dm = dm;
      _ba = ba;
      _dk = dk;
      _mk = mk;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat1 v1x = (VecArrayFloat1)vx;
      VecArrayFloat1 v1y = (VecArrayFloat1)vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      float[] z = copy(x);
      _sf.apply(z);
      applyLhs(_dm,_ba,_dk,_mk,z,y);
      _sf.applyTranspose(y);
    }
    private Smoother _sf;
    private float _dm,_ba;
    private float[] _dk,_mk;
  }

  // Smoother used as a preconditioner.
  private static class Smoother {
    public Smoother(
      int nm, float sigma) 
    {
      _sigma = sigma;
    }
    public void apply(float[] x) {
      smooth1(_sigma,x);
    }
    public void applyTranspose(float[] x) {
      smooth1(_sigma,x);
    }
    private float _sigma;
  }

  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply(x,x);
  }

  // Line search for each Gauss-Newtom iterate
  private float lineSearch(final float dm, final float beta, 
    final float[] d, final float[] p, final float[] m)
  {
    int nd = d.length;
    final float[] dk  = new float[nd];
    LineSearch.Function func = new LineSearch.Function() {
      public double[] evaluate(double s) {
        float sm = (float)s;
        float sp = (float)s;
        float[] mk = add(mul(sm,p),m);
        double f = (double) objFunc(beta,dm,d,mk); 
        sm -= sc;
        sp += sc;
        float[] mm = add(mul(sm,p),m);
        float[] mp = add(mul(sp,p),m);
        double fp = (double) objFunc(beta,dm,d,mp); 
        double fm = (double) objFunc(beta,dm,d,mm); 
        double g = (fp-fm)/(sc*2.0);
        return new double[] {f,g};
      }
      private float sc = 0.001f;
    };
    double s = 1.0e-3;
    LineSearch ls = new LineSearch(func,1.0e-10,1.0e-3,1.0e-3);
    double[] fg = func.evaluate(0.0);
    double f = fg[0];
    double g = fg[1];
    if(g>=0.0) {g=-1.0;}
    System.out.println("f="+f);
    System.out.println("g="+g);
    double smin = 0.0;
    double smax = 4.0*max(1.0,s);
    LineSearch.Result lsr = ls.search(s,f,g,smin,smax);
    System.out.println("s="+lsr.s);
    return (float) lsr.s;
  }
  // compute the total object function
  private float objFunc(float beta, float dm, float[] d, float[] mk) {
    int nd = d.length;
    float[] dk = new float[nd];
    forward(mk,dk);
    float dmk = dataMisfit(sub(d,dk));
    float mnk = modelNorm2(dm,mk);
    return (dmk+beta*mnk);
  } 

  private float dataMisfit(float[] dok) {
    int nd = dok.length;
    float wd = 1.0f/_dSigma;
    float ws = wd*wd;
    float sum = 0.0f;
    for (int id=0; id<nd; ++id) {
      sum += ws*dok[id]*dok[id];
    } 
    return sum;
  }

  private float modelNorm2(float dm, float[] m) {
    float sum = 0.0f;
    int nm = m.length;
    float[] ym = new float[nm]; 
    applyWmWm(dm,m,ym);
    for (int im=0; im<nm; ++im) {
      sum += ym[im]*ym[im];
    } 
    return sum;
  }

  private void makeRhs(
    float dm, float beta, float[] dk, float[] dok, float[] mk, float[] y) 
  {
    int nm = mk.length;
    int nd = dk.length;
    float[] yd = new float[nd]; 
    float[] y1 = new float[nm]; 
    float[] y2 = new float[nm]; 
    applyWdWd(dok,yd);
    applyJacbT(dk,mk,yd,y1);
    applyWmWm(dm,mk,y2);
    sub(y1,mul(beta,y2),y);
  }
  
  private void applyLhs(
    float dm, float beta, float[] dk, float[] mk, float[] x, float[] y) 
  {
    int nm = mk.length;
    int nd = dk.length;
    float[] y1 = new float[nm];
    float[] y2 = new float[nm];
    float[] yd = new float[nd];
    applyJacb(dk,mk,x,yd);
    applyWdWd(yd,yd);
    applyJacbT(dk,mk,yd,y1);
    applyWmWm(dm,x,y2);
    add(y1,mul(beta,y2),y);
  }

  // apply W'dWd operator
  private void applyWdWd(float[] x, float[] y) {
    float wd = 1.0f/_dSigma; 
    float ws = wd*wd;
    mul(ws,x,y);
  }

  // apply W'mWm operator
  private void applyWmWm(float dm, float[] x, float[] y) {
    int n = x.length;
    float ws = _alphaS*dm;
    float dx = _alphaX/dm;
    float d11 = dx;
    float d12 = dx;
    float d22 = dx;
    for (int i=1; i<n;++i) {
      float xa = 0.0f;
      float xb = 0.0f;
      xa += x[i  ];
      xb -= x[i-1];
      float ya = d11*xa+d12*xb;
      float yb = d12*xa+d22*xb;
      y[i  ] += ya;
      y[i-1] -= yb;
    }
    add(mul(ws,x),y,y);
  }

  // apply J operator 
  private void applyJacb(
    float[] dk, float[] mk, float[] x, float[] y) 
  {
    zero(y);
    int nd = dk.length;
    int nm = mk.length;
    float[] dp = new float[nd];
    for (int im=0; im<nm; ++im) {
      mk[im] += _dh; 
      forward(mk,dp);
      mk[im] -= _dh; 
      for (int id=0; id<nd; ++id) 
        y[id] += x[im]*(dp[id]-dk[id])/_dh; 
    }
  }

  // apply J' operator 
  private void applyJacbT(
    float[] dk, float[] mk, float[] x, float[] y) 
  {
    zero(y);
    int nd = dk.length;
    int nm = mk.length;
    float[] dp = new float[nd];
    for (int im=0; im<nm; ++im) { 
      mk[im] += _dh;
      forward(mk,dp);
      mk[im] -= _dh;
      for (int id=0; id<nd; ++id) 
        y[im] +=x[id]*(dp[id]-dk[id])/_dh; 
    }
  }

  // compute predicted data
  public static void forward(float[] mk, float[] dk) {
    zero(dk);
    int n = mk.length;
    for (int i=0; i<n; ++i) {
      float zbi =  mk[i];
      float zti = _zt[i];
      float xci = _xc[i];
      float[] dki = forward(xci,zti,zbi);
      add(dki,dk,dk); 
    } 
  }
  // forward for each prism
  public static float[] forward(float xc, float zt, float zb) {
    // construct the "polygon" representing the vertical dyke
    int np = 4;
    float[] xp = new float[np];  
    float[] zp = new float[np];  
    float swd1 = 0.5f*_wd;
    float swd2 = 0.0001f*_wd;
    xp[0] = xc-swd1; zp[0] = zt-swd2;
    xp[1] = xc+swd1; zp[1] = zt+swd2;
    xp[2] = xp[1];   zp[2] = zb+swd2;
    xp[3] = xp[0];   zp[3] = zb-swd2;
    float gcons = 0.006672f; 
    int nd = _xo.length;
    float[] sum = zerofloat(nd);
    for (int ip=0; ip<np; ++ip) {
      int ipp = ip+1;
      if (ip==np-1) ipp = 0;
      float xpi = xp[ip];
      float zpi = zp[ip];
      float xpe = xp[ipp];
      float zpe = zp[ipp];
      float dxi = abs(xpe-xpi);
      float dzi = abs(zpe-zpi);
      if (dzi<=0.000001f) zpi +=0.00001f*dxi; 
      float[] x1 = sub(xpi,_xo);
      float[] x2 = sub(xpe,_xo);
      float[] z1 = sub(zpi,_zo);
      float[] z2 = sub(zpe,_zo);
      float[] r1 = add(mul(x1,x1),mul(z1,z1));
      float[] r2 = add(mul(x2,x2),mul(z2,z2));
      float[] bt = sub(z2,z1);
      float[] alpha  = div(sub(x2,x1),bt);
      float[] beta   = div(sub(mul(x1,z2),mul(x2,z1)),bt);
      float[] factor = div(beta, add(1.0f,mul(alpha,alpha)));
      float[] term1  = mul(0.5f,log(div(r2,r1)));
      float[] term2  = sub(atan(z2,x2), atan(z1,x1)); 
      float[] update = sub(term1, mul(alpha,term2));
      sum = add(sum, mul(factor,update));
    }
    float sca = 2.0f*_dc*gcons;
    return mul(sca,sum); 
  }

  private static float[] atan(float[] y, float[] x) {
    int n = x.length;
    float[] r = new float[n];
    for (int i=0; i<n; ++i)
      r[i] = atan2(y[i],x[i]); 
    return r;
  }
  private static void printStats(String s, int i1, float[][] a) {
    int n2 = a.length;
    float amin = a[0][i1];
    float amax = a[0][i1];
    for (int i2=1; i2<n2; ++i2) {
      if (a[i2][i1]<amin) amin = a[i2][i1];
      if (a[i2][i1]>amax) amax = a[i2][i1];
    }
    System.out.println(s+": i1="+i1+" min="+amin+" max="+amax);
  }
}
