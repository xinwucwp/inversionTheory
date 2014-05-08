/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package hv;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * Flattens and unflattens locally planar features in a 3D image.
 * <p>
 * We assume that a 3D image f(x1,x2,x3) has coherent features that are
 * locally planar with finite (non-vertical) slope. This class computes and
 * applies coordinate mappings for flattening and unflattening such features.
 * In this version, the mappings change only the vertical coordinate.
 * <p> 
 * Let x1 denote the vertical coordinate for the image f(x1,x2,x3). This class
 * computes a mapping function x1(u1,u2,u3) such that features in a
 * transformed image g(u1,u2,u3) = f(x1(u1,u2),u2,u3) are approximately
 * horizontal. This process is often called "flattening", and the transformed
 * image g(u1,u2,u3) is "the flattened image." Likewise, for any constant u1,
 * the curve defined by the function x1(u1,x2,x3) is called a "horizon".
 * <p>
 * If the coordinates u1, u2 and u3 are sampled finely enough, then the
 * mapping function x1(u1,u2,u3) is invertible. Specifically, there exists an
 * inverse mapping u1(x1,x2,x3) such that x1 = x1(u1(x1,x2,x3),x2,u3) for all
 * x2 and x3. Currently, samplings of u1, u2 and u3 are set to be identical to
 * those for x1, x2 and x3. If used directly, this sampling of the mapping
 * x1(u1,u2,u3) may cause aliasing of the flattened image g(u1,u2,u3), so that
 * f(x1,x2,x3) cannot be recovered by unflattening.
 * <p>
 * Without additional constraints, the description above of the mapping
 * function x1(u1,u2) is ambiguous. For example, one could add any constant to
 * this mapping function, and features in the transformed image g(u1,u2) would
 * still be horizontal. The mapping function x1(u1,u2) is here constrained so
 * that x1(u1,u2) = u1 for the first and last sampled values of u1. In other
 * words, features at the top or bottom of an image are not shifted by
 * flattening or unflattening.
 * @author Dave Hale and Xinming Wu, Colorado School of Mines
 * @version 2013.01.29
 */
public class RelativeGeologicTime3C {

  /** Coordinate mappings u1(x1,x2,x3) and x1(u1,u2,u3). */
  public static class Mappings {
    
    /** Sampling for the 1st dimension (the vertical coordinate). */
    public Sampling s1;
    
    /** Sampling for the 2nd dimension. */
    public Sampling s2;
    
    /** Sampling for the 3rd dimension. */
    public Sampling s3;

    /** Array of sampled u1(x1,x2,x3). */
    public float[][][] u1;
    
    /** Array of sampled x1(u1,u2,u3). */
    public float[][][] x1;

    /**
     * Uses these mappings to flatten the specified image.
     * @param f the image to flatten.
     * @return the flattened image.
     */
    public float[][][] flatten(float[][][] f) {
      return apply(x1,f);
    }

    /**
     * Uses these mappings to unflatten the specified image.
     * @param f the image to unflatten.
     * @return the unflattened image.
     */
    public float[][][] unflatten(float[][][] f) {
      return apply(u1,f);
    }

    /**
     * Gets the flattening shifts s(u1,u2,u3) = u1 - x1(u1,u2,u3).
     * @return the flattening shifts.
     */
    public float[][][] getShiftsS() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n3 = s3.getCount();
      float[][][] s = new float[n3][n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float u1i = f1+i1*d1;
            s[i3][i2][i1] = u1i-x1[i3][i2][i1];
          }
        }
      }
      return s;
    }

    /**
     * Gets the unflattening shifts r(x1,x2,x3) = u1(x1,x2,x3) - x1.
     * @return the unflattening shifts.
     */
    public float[][][] getShiftsR() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n3 = s3.getCount();
      float[][][] r = new float[n3][n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float x1i = f1+i1*d1;
            r[i3][i2][i1] = u1[i3][i2][i1]-x1i;
          }
        }
      }
      return r;
    }

    private Mappings(
      Sampling s1, Sampling s2, Sampling s3, 
      float[][][] u1, float[][][] x1) 
    {
      this.s1 = s1;
      this.s2 = s2;
      this.s3 = s3;
      this.u1 = u1;
      this.x1 = x1;
    }

    private float[][][] apply(final float[][][] ux, final float[][][] f) {
      final int n1 = s1.getCount();
      final int n2 = s2.getCount();
      final int n3 = s3.getCount();
      final double d1 = s1.getDelta();
      final double f1 = s1.getFirst();
      final SincInterp si = new SincInterp();
      final float[][][] g = new float[n3][n2][n1];
      Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2)
          si.interpolate(n1,d1,f1,f[i3][i2],n1,ux[i3][i2],g[i3][i2]);
      }});
      return g;
    }
  }

  /**
   * Sets the relative weight of the PDE dr(x1,x2)/dx1 = 0.
   * Increasing this weight will cause shifts r(x1,x2) and s(u1,u2) to vary
   * more slowly with vertical coordinates x1 and u1, respectively. A weight
   * of 1.0 will cause this equation to get as much weight as other PDEs that
   * cause contours of constant u1 = u1(x1,x2) to be aligned with coherent
   * planar image features.
   * @param w1 the weight.
   */
  public void setWeight1(double w1) {
    _weight1 = (float)w1;
  }
  
  public void setScale(double sc) {
    _scale = (float)sc;
  }

  /**
   * Sets half-widths for smoothings in 1st and 2nd dimensions.
   * These smoothings serve as a preconditioner; they accelerate convergence
   * of the iterative solver used to compute mappings.
   * @param sigma1 half-width for smoothing in 1st dimension, in samples.
   * @param sigma2 half-width for smoothing in 2nd dimension, in samples.
   */
  public void setSmoothings(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  /**
   * Sets parameters that control the number of solver iterations.
   * @param small stop iterations when error norm is reduced by this fraction.
   * @param niter maximum number of solver iterations.
   */
  public void setIterations(double small, int niter) {
    _small = (float)small;
    _niter = niter;
  }

  /**
   * Gets mappings computed from specified slopes and planarities.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param p2 array of slopes of image features.
   * @param ep array of planarities of image features.
   */
  public Mappings getMappingsFromSlopes(
    Sampling s1, Sampling s2, Sampling s3,
    float[][][] uz, float[][][] uy, float[][][] ux, 
    float[][][] wp, float[][][] ws, float[][] rd, 
    float[][] k1, float[][] k2, float[][] k3) 
  {
    // Sampling parameters.
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();
    final int n3 = s3.getCount();
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    float d3 = (float)s3.getDelta();
    float f1 = (float)s1.getFirst();
    int nc = k1.length;
    int[][] k1a = new int[nc][];
    int[][] k2a = new int[nc][];
    int[][] k3a = new int[nc][];
    float[][][] wc = zerofloat(n1,n2,n3);
    add(wc,0.08f,wc);
    for (int ic=0; ic<nc; ++ic) {
      int np = k1[ic].length;
      k1a[ic] = new int[np];
      k2a[ic] = new int[np];
      k3a[ic] = new int[np];
      for (int ip=0; ip<np; ++ip) {
        int i1 = (int)k1[ic][ip];
        int i2 = (int)k2[ic][ip];
        int i3 = (int)k3[ic][ip];
        k1a[ic][ip] = i1;
        k2a[ic][ip] = i2;
        k3a[ic][ip] = i3;
        wc[i3][i2][i1] = 0.08f;
        if(i1+1<n1) wc[i3][i2][i1+1] = 0.08f;
        if(i1-1>=0) wc[i3][i2][i1-1] = 0.08f;
      }
    }

    // Compute shifts r(x1,x2,x3), in samples.
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] u = new float[n3][n2][n1]; // shifts, in samples
    initializeRGTs(rd,k1a,k2a,k3a,u); // initial shifts to satisfy constraints
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vu = new VecArrayFloat3(u);
    A3 a3 = new A3(_scale,_weight1,wp,uz,uy,ux);
    M3 m3 = new M3(_sigma1,_sigma2,_sigma3,ws,k1a,k2a,k3a);
    //testSpd("a2",n1,n2,a2);
    //testSpd("m2",n1,n2,m2);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(_weight1,b);
    cs.solve(a3,m3,vb,vu);
    final float[][][] r = u;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1i = f1+i1*d1;
          r[i3][i2][i1] = (u[i3][i2][i1]-x1i)/d1;
        }
      }
    }

    //checkShifts(k1a,k2a,k3a,r);
    cleanShifts(r);
    // Compute u1(x1,x2,x3).
    final float[][][] u1 = r;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1i = f1+i1*d1;
          u1[i3][i2][i1] = x1i+r[i3][i2][i1]*d1;
        }
      }
    }

    // Compute x1(u1,u2).
    final float[][][] x1 = b;
    final InverseInterpolator ii = new InverseInterpolator(s1,s1);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) 
        ii.invert(u1[i3][i2],x1[i3][i2]);
    }});

    return new Mappings(s1,s2,s3,u1,x1);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _scale = 0.0f;    // scale the curvature term
  private float _weight1 = 0.01f; // weight for dr/d1 = 0 equation
  private float _sigma1 = 4.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 4.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 4.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 1000; // maximum number of CG iterations

  // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(float sc, float w1, float[][][] wp,
       float[][][] uz, float[][][] uy, float[][][] ux) 
    {
      _wp = wp;
      _uz = uz;
      _uy = uy;
      _ux = ux;
      _sc = sc;
      _w1 = w1;
      //testSpd();
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      float[][][] y1 = copy(x);
      float[][][] y2 = copy(x);
      VecArrayFloat3 v3yy = new VecArrayFloat3(y2);
      v3yy.zero();
      zero(y);
      applyLhs(_w1,_wp,_uz,_uy,_ux,z,y);
      if (_sc>0.0f) {
        applyLaplace(z ,y1);
        applyLaplace(y1,y2);
        v3y.add(1.f,v3yy,_sc);
      }
    }
    private float _sc;
    private float _w1;
    private float[][][] _wp;
    private float[][][] _uz;
    private float[][][] _uy;
    private float[][][] _ux;
    public void testSpd() {
      // symmetric: y'Ax = x'(A'y) = x'Ay
      // positive-semidefinite: x'Ax >= 0
      int n1 = _wp[0][0].length;
      int n2 = _wp[0].length;
      int n3 = _wp.length;
      float[][][] x = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] y = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] ax = zerofloat(n1,n2,n3);
      float[][][] ay = zerofloat(n1,n2,n3);
      VecArrayFloat3 vx = new VecArrayFloat3(x);
      VecArrayFloat3 vy = new VecArrayFloat3(y);
      VecArrayFloat3 vax = new VecArrayFloat3(ax);
      VecArrayFloat3 vay = new VecArrayFloat3(ay);
      apply(vx,vax);
      apply(vy,vay);
      double yax = vy.dot(vax);
      double xay = vx.dot(vay);
      double xax = vx.dot(vax);
      double yay = vy.dot(vay);
      System.out.println("A3: yax="+yax+" xay="+xay);
      System.out.println("A3: xax="+xax+" yay="+yay);
    }
  }

  // Preconditioner; includes smoothers and (optional) constraints.
  private static class M3 implements CgSolver.A {
    M3(float sigma1, float sigma2, float sigma3, float[][][] wp, 
       int[][] k1, int[][] k2, int[][] k3) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
      _wp = wp;
      if (k1!=null && k2!=null && k3!=null) {
        _k1 = copy(k1);
        _k2 = copy(k2);
        _k3 = copy(k3);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      copy(x,y);
      constrain(_k1,_k2,_k3,y);
      removeAverage(y);
      smooth3(_sigma3,_wp,y);
      smooth2(_sigma2,_wp,y);
      smooth1(2.0f*_sigma1,_wp,y);
      smooth2(_sigma2,_wp,y);
      smooth3(_sigma3,_wp,y);
      removeAverage(y);
      constrain(_k1,_k2,_k3,y);
    }
    private float _sigma1,_sigma2,_sigma3;
    private float[][][] _wp;
    private int[][] _k1,_k2,_k3;
  }

  public static void initializeRGTs
    (float[][] rd, int[][] k1, int[][] k2, 
     int[][] k3, float[][][] u) 
  {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          u[i3][i2][i1] = (float)i1;
    if (k1!=null && k2!=null &&k3!=null) {
      int nc = k1.length;
      for (int ic=0; ic<nc; ++ic) {
        int nk = k1[ic].length;
        int i1 = k1[ic][0];
        int i2 = k2[ic][0];
        int i3 = k3[ic][0];
        u[i3][i2][i1] = (float)i1;
        float ur = (float)i1 + rd[ic][0]; 
        for (int ik=1; ik<nk; ++ik) {
          i1  = k1[ic][ik];
          i2  = k2[ic][ik];
          i3  = k3[ic][ik];
          u[i3][i2][i1] = ur-rd[ic][ik];
        }
      }
    }
  }

  public static void checkShifts(int[][] k1, int[][] k2, int[][] k3, float[][][] u) {
    if (k1!=null && k2!=null &&k3!=null) {
      int nc = k1.length;
      for (int ic=0; ic<nc; ++ic) {
        trace("ic="+ic);
        int nk = k1[ic].length;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = k1[ic][ik];
          int i2 = k2[ic][ik];
          int i3 = k3[ic][ik];
          trace("  i1="+i1+" i2="+i2+" i3="+i3+" u="+u[i3][i2][i1]);
          //assert r[i2][i1]==rp+ip-i1:"shifts r satisfy constraints";
        }
      }
    }
  }

  public static void constrain(int[][] k1, int[][] k2, int[][] k3, float[][][] x) {
    if (k1!=null && k2!=null &&k3!=null) {
      int nc = k1.length;
      for (int ic=0; ic<nc; ++ic) {
        int nk = k1[ic].length;
        float sum = 0.0f;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = k1[ic][ik];
          int i2 = k2[ic][ik];
          int i3 = k3[ic][ik];
          sum += x[i3][i2][i1];
        }
        float avg = sum/(float)nk;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = k1[ic][ik];
          int i2 = k2[ic][ik];
          int i3 = k3[ic][ik];
          x[i3][i2][i1] = avg;
        }
      }
    }
  }
  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }
  
  private static void smooth1(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    int n2 = x.length;
    int n1 = x[0].length;
    float c = 0.5f*sigma*sigma;
    float[] st = fillfloat(1.0f,n1);
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      if (s!=null) {
        for (int i1=0; i1<n1; ++i1)
          st[i1] = s[i2][i1];
      }
      for (int i1=0; i1<n1; ++i1)
        xt[i1] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }
  private static void smooth1(final float sigma, final float[][][] s, final float[][][] x) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] x3 = x[i3];
      float[][] s3 = (s!=null)?s[i3]:null;
      smooth1(sigma,s3,x3);
    }});
  }
  


  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] st = fillfloat(1.0f,n2);
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      if (s!=null) {
        for (int i2=0; i2<n2; ++i2)
          st[i2] = s[i2][i1];
      }
      for (int i2=0; i2<n2; ++i2)
        xt[i2] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }
  private static void smooth2(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] s3 = (s!=null)?s[i3]:null;
      float[][] x3 = x[i3];
      smooth2(sigma,s3,x3);
    }});
  }

  // Smoothing for dimension 3.
  private static void smooth3(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] s2 = (s!=null)?new float[n3][]:null;
      float[][] x2 = new float[n3][];
      for (int i3=0; i3<n3; ++i3) {
        if (s!=null)
          s2[i3] = s[i3][i2];
        x2[i3] = x[i3][i2];
      }
      smooth2(sigma,s2,x2);
    }});
  }

  private static void removeAverage(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float nh = (float)(n2*n3);
    for (int i1=0; i1<n1; ++i1) {
      float sumx = 0.0f;
      for (int i3=0; i3<n3; ++i3)  
        for (int i2=0; i2<n2; ++i2)  
          sumx += x[i3][i2][i1];
      float avgx = sumx/nh;
      for (int i3=0; i3<n3; ++i3) 
        for (int i2=0; i2<n2; ++i2) 
          x[i3][i2][i1] -= avgx; 
    }
  }

  private static void makeRhs(
    float w1, float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    float w1s = w1*w1;
    float ya = 0.25f*w1s;
    float yb = 0.25f*w1s;
    float yc = 0.25f*w1s;
    float yd = 0.25f*w1s;
    for (int i3=1,i3m=0; i3<n3; ++i3,++i3m) {
      for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
        for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
          y[i3 ][i2 ][i1 ] += ya;
          y[i3 ][i2 ][i1m] -= yd;
          y[i3 ][i2m][i1 ] += yb;
          y[i3 ][i2m][i1m] -= yc;
          y[i3m][i2 ][i1 ] += yc;
          y[i3m][i2 ][i1m] -= yb;
          y[i3m][i2m][i1 ] += yd;
          y[i3m][i2m][i1m] -= ya;
        }
      }
    }
  }

  private static void applyLaplace( final float[][][]x, final float[][][] y)
  {   
    final int n3 = x.length; 
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLaplaceSlice3(i3,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLaplaceSlice3(i3,x,y);
    }});
  }
  private static void applyLaplaceSlice3(
    int i3,float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x01 = x[i3  ][i2-1];
      float[] x10 = x[i3-1][i2  ];
      float[] x11 = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y01 = y[i3  ][i2-1];
      float[] y10 = y[i3-1][i2  ];
      float[] y11 = y[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float d11 = 1.0f;
        float d22 = 1.0f;
        float d33 = 1.0f;
        float xa = 0.0f;
        float xb = 0.0f;
        float xc = 0.0f;
        float xd = 0.0f;
        xa += x00[i1 ];
        xd -= x00[i1m];
        xb += x01[i1 ];
        xc -= x01[i1m];
        xc += x10[i1 ];
        xb -= x10[i1m];
        xd += x11[i1 ];
        xa -= x11[i1m];
        float x1 = 0.25f*(xa+xb+xc+xd);
        float x2 = 0.25f*(xa-xb+xc-xd);
        float x3 = 0.25f*(xa+xb-xc-xd);
        float y1 = d11*x1;
        float y2 = d22*x2;
        float y3 = d33*x3;
        float ya = 0.25f*(y1+y2+y3);
        float yb = 0.25f*(y1-y2+y3);
        float yc = 0.25f*(y1+y2-y3);
        float yd = 0.25f*(y1-y2-y3);
        y00[i1 ] += ya;
        y00[i1m] -= yd;
        y01[i1 ] += yb;
        y01[i1m] -= yc;
        y10[i1 ] += yc;
        y10[i1m] -= yb;
        y11[i1 ] += yd;
        y11[i1m] -= ya;
      }
    }
  }

  private static void applyLhs(
    final float w1, final float[][][] wp, 
    final float[][][] uz, final float[][][] uy, final float[][][] ux,
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,w1,wp,uz,uy,ux,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,w1,wp,uz,uy,ux,x,y);
    }});
  }
  private static void applyLhsSlice3(
    int i3, float w1, float[][][] wp, 
    float[][][] uz,float[][][] uy, float[][][] ux,
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    float w1s = w1*w1;
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x01 = x[i3  ][i2-1];
      float[] x10 = x[i3-1][i2  ];
      float[] x11 = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y01 = y[i3  ][i2-1];
      float[] y10 = y[i3-1][i2  ];
      float[] y11 = y[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
        float uzi = uz[i3][i2][i1];
        float uyi = uy[i3][i2][i1];
        float uxi = ux[i3][i2][i1];
        float wps = wpi*wpi;
        float uzs = uzi*uzi;
        float uys = uyi*uyi;
        float uxs = uxi*uxi;
        float d11 = w1s+wps*(uys+uxs);
        float d12 = -wps*uyi*uzi;
        float d13 = -wps*uxi*uzi;
        float d22 = wps*uzs;
        float d33 = wps*uzs;
        float xa = 0.0f;
        float xb = 0.0f;
        float xc = 0.0f;
        float xd = 0.0f;
        xa += x00[i1 ];
        xd -= x00[i1m];
        xb += x01[i1 ];
        xc -= x01[i1m];
        xc += x10[i1 ];
        xb -= x10[i1m];
        xd += x11[i1 ];
        xa -= x11[i1m];
        float x1 = 0.25f*(xa+xb+xc+xd);
        float x2 = 0.25f*(xa-xb+xc-xd);
        float x3 = 0.25f*(xa+xb-xc-xd);
        float y1 = d11*x1+d12*x2+d13*x3;
        float y2 = d12*x1+d22*x2       ;
        float y3 = d13*x1       +d33*x3;
        float ya = 0.25f*(y1+y2+y3);
        float yb = 0.25f*(y1-y2+y3);
        float yc = 0.25f*(y1+y2-y3);
        float yd = 0.25f*(y1-y2-y3);
        y00[i1 ] += ya;
        y00[i1m] -= yd;
        y01[i1 ] += yb;
        y01[i1m] -= yc;
        y10[i1 ] += yc;
        y10[i1m] -= yb;
        y11[i1 ] += yd;
        y11[i1m] -= ya;
      }
    }
  }
  private static void trace(String s) {
    System.out.println(s);
  }
  // Post-processing of computed shifts to ensure monotonic u1.
  private static void cleanShifts(float[][][] r) {
    int n1 = r[0][0].length;
    int n2 = r[0].length;
    int n3 = r.length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          if (r[i3][i2][i1]<=r[i3][i2][i1-1]-0.99f)
            r[i3][i2][i1] = r[i3][i2][i1-1]-0.99f;
        }
      }
    }
  }
}
