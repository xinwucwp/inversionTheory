/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package hv;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * Generate an RGT volume using seismic normal vectors estimated from a 3D seismic image.

 * Contours of constant RGT values are seismic horizons. Therefore, RGT countours should 
 * have the same structures as seismic reflectors, and hence the gradient vectors 
 * of an RGT volume, that are perpendicular to RGT contours, should be parallel to 
 * seismic normal vectors, that are perpendicular to seismic reflectors.
 *
 * To generate a more accurate RGT volume for a complicated seismic image, some control 
 * points are incorporated to the method, and are implemented as a simple preconditioner 
 * for the conjugate gradient method that is used to generate an RGT volume.
 *
 * From an RGT volume t(x,z) with RGT t monotonically increasing with its vertical 
 * coordinate z, a horizon volume z(x,t) can be obtained by inverse linear interpolation, 
 * and then be used to flatten a seismic image.
 * @author Dave Hale and Xinming Wu, Colorado School of Mines
 * @version 2014.04.05
 */
public class RelativeGeologicTime2C {

  /** Coordinate mappings u1(x1,x2) and x1(u1,u2). */
  public static class Mappings {
    
    /** Sampling for the 1st dimension (the vertical coordinate). */
    public Sampling s1;
    
    /** Sampling for the 2nd dimension. */
    public Sampling s2;

    /** Array of sampled u1(x1,x2). */
    public float[][] u1;
    
    /** Array of sampled x1(u1,u2). */
    public float[][] x1;

    /**
     * Uses these mappings to flatten the specified image.
     * @param f the image to flatten.
     * @return the flattened image.
     */
    public float[][] flatten(float[][] f) {
      return apply(x1,f);
    }

    /**
     * Uses these mappings to unflatten the specified image.
     * @param f the image to unflatten.
     * @return the unflattened image.
     */
    public float[][] unflatten(float[][] f) {
      return apply(u1,f);
    }

    /**
     * Gets the flattening shifts s(u1,u2) = u1 - x1(u1,u2).
     * @return the flattening shifts.
     */
    public float[][] getShiftsS() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      float[][] s = new float[n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = f1+i1*d1;
          s[i2][i1] = u1i-x1[i2][i1];
        }
      }
      return s;
    }

    /**
     * Gets the unflattening shifts r(x1,x2) = u1(x1,x2) - x1.
     * @return the unflattening shifts.
     */
    public float[][] getShiftsR() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      float[][] r = new float[n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1i = f1+i1*d1;
          r[i2][i1] = u1[i2][i1]-x1i;
        }
      }
      return r;
    }

    private Mappings(Sampling s1, Sampling s2, float[][] u1, float[][] x1) {
      this.s1 = s1;
      this.s2 = s2;
      this.u1 = u1;
      this.x1 = x1;
    }

    private float[][] apply(float[][] ux, float[][] f) {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      double d1 = s1.getDelta();
      double f1 = s1.getFirst();
      SincInterp si = new SincInterp();
      float[][] g = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2)
        si.interpolate(n1,d1,f1,f[i2],n1,ux[i2],g[i2]);
      return g;
    }
  }

  /**
   * Sets the relative weight of the PDE dr(x1,x2)/dx1 = 0.
   * Increasing this weight will cause shifts r(x1,x2) and s(u1,u2) to vary
   * more slowly with vertical coordinates x1 and u1, respectively. A weight
   * of 1.0 will cause this equation to get as much weight as other PDEs that
   * cause contours of constant u1 = u1(x1,x2) to be aligned with coherent
   * linear image features.
   * @param w1 the weight.
   */
  public void setWeight1(double w1) {
    _weight1 = (float)w1;
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
   * Gets mappings computed from specified slopes and linearities.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param p2 array of slopes of image features.
   * @param el array of linearities of image features.
   */
  public Mappings getMappingsFromSlopes(
    Sampling s1, Sampling s2, float[][] nx, float[][] nz, float[][] el) 
  {
    return getMappingsFromSlopes(s1,s2,nx,nz,el,null);
  }

  /**
   * Gets mappings computed from specified slopes and constraints.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param p2 array of slopes of image features.
   * @param el array of linearities of image features.
   * @param cs constraints, arrays of (x1,x2) for which u1 is constant.
   *  Contents of the array cs are as follows:
   *  <pre>
   *  {{{x1_00,x1_01,x1_02,...},{x1_10,x1_11,...},...},
   *   {{x2_00,x2_01,x2_02,...},{x2_10,x2_11,...},...}}
   *  </pre>
   *  Each constraint point has coordinates (x1,x2). In the array cs, the
   *  arrays of arrays of x1 coordinates precede the arrays of arrays of x2
   *  coordinates. Each array of x1 coordinates has a corresponding array of
   *  x2 coordinates, and together this pair of arrays of (x1,x2) represents
   *  one constraint. All points in a constraint will be constrained to have
   *  the same u1; in other words, they will all lie on the same horizon.
   */
  public Mappings getMappingsFromSlopes(
    Sampling s1, Sampling s2, float[][] uz, float[][] ux, float[][] el, float[][][] c)
  {
    // Sampling parameters.
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    float f1 = (float)s1.getFirst();

    // Convert any constraints in cs to sample indices k1 and k2.
    int[][] k1 = null;
    int[][] k2 = null;
    if (c!=null) {
      float[][] c1 = c[0];
      float[][] c2 = c[1];
      int nc = c1.length;
      k1 = new int[nc][];
      k2 = new int[nc][];
      for (int ic=0; ic<nc; ++ic) {
        int nk = c1[ic].length;
        k1[ic] = new int[nk];
        k2[ic] = new int[nk];
        for (int ik=0; ik<nk; ++ik) {
          k1[ic][ik] = s1.indexOfNearest(c1[ic][ik]);
          k2[ic][ik] = s2.indexOfNearest(c2[ic][ik]);
        }
      }
    }


    // Compute shifts r(x1,x2), in samples.
    float[][] b = new float[n2][n1]; // for right-hand side b
    float[][] u = new float[n2][n1]; // RGTs u, in samples
    initializeRGTs(k1,k2,u); // initial shifts to satisfy constraints
    checkRGTs(k1,k2,u);
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vu = new VecArrayFloat2(u);
    A2 a2 = new A2(_weight1,el,uz,ux);
    M2 m2 = new M2(_sigma1,_sigma2,el,k1,k2);
    //testSpd("a2",n1,n2,a2);
    //testSpd("m2",n1,n2,m2);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(_weight1,b);
    cs.solve(a2,m2,vb,vu);
    checkRGTs(k1,k2,u);
    float[][] r = u;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        r[i2][i1] = (u[i2][i1]-f1)/d1-i1;
      }
    }
    cleanShifts(r);

    // Compute u1(x1,x2) from shifts r.
    float[][] u1 = r;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        u1[i2][i1] = f1+(i1+r[i2][i1])*d1;
      }
    }

    // Compute x1(u1,u2) using inverse linear interpolation.
    float[][] x1 = b;
    InverseInterpolator ii = new InverseInterpolator(s1,s1);
    for (int i2=0; i2<n2; ++i2) 
      ii.invert(u1[i2],x1[i2]);
    printStats("u1",u1);
    printStats("x1",x1);
    return new Mappings(s1,s2,u1,x1);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _weight1 = 0.010f; // weight for dr/d1 = 0 equation
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 1000; // maximum number of CG iterations

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(float w1, float[][] wp, float[][] uz, float[][] ux) {
      _w1 = w1;
      _wp = wp;
      _ux = ux;
      _uz = uz;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      applyLhs(_w1,_wp,_uz,_ux,x,y);
    }
    private float _w1;
    private float[][] _wp;
    private float[][] _uz;
    private float[][] _ux;
  }

  // Preconditioner; includes smoothers and (optional) constraints.
  private static class M2 implements CgSolver.A {
    M2(float sigma1, float sigma2, float[][] wp, int[][] k1, int[][] k2) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _wp = wp;
      if (k1!=null && k2!=null) {
        _k1 = copy(k1);
        _k2 = copy(k2);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      copy(x,y);
      constrain(_k1,_k2,y);
      //removeAverage(y);
      smooth2(_sigma2,_wp,y);
      smooth1(2.0f*_sigma1,y);
      smooth2(_sigma2,_wp,y);
      //removeAverage(y);
      constrain(_k1,_k2,y);
    }
    private float _sigma1,_sigma2;
    private float[][] _wp;
    private int[][] _k1,_k2;
  }

  // Initializes RGT values u to satisfy constraints that u1 = u2=... .
  public static void initializeRGTs(int[][] k1, int[][] k2, float[][] u) {
    int n2 = u.length; 
    int n1 = u[0].length; 
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        u[i2][i1] = (float)i1;
      }
    }
    if (k1!=null && k2!=null) {
      int nc = k1.length;
      for (int ic=0; ic<nc; ++ic) {
        int nk = k1[ic].length;
        float ur = k1[ic][0]; 
        for (int ik=0; ik<nk; ++ik) {
          int i1 = k1[ic][ik];
          int i2 = k2[ic][ik];
          u[i2][i1] = ur;
        }
      }
    }
  }

  // Asserts that RGT values u satisfy constraints that u is constant.
  public static void checkRGTs(int[][] k1, int[][] k2, float[][] u) {
    if (k1!=null && k2!=null) {
      int nc = k1.length;
      for (int ic=0; ic<nc; ++ic) {
        trace("ic="+ic);
        int nk = k1[ic].length;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = k1[ic][ik];
          int i2 = k2[ic][ik];
          trace("  i1="+i1+" i2="+i2+" u="+u[i2][i1]);
          //assert r[i2][i1]==rp+ip-i1:"shifts r satisfy constraints";
        }
      }
    }
  }

  // Projects x into the null space of constraints. For each constraint (if
  // any), replaces values at constrained samples with the averages of those
  // values.
  public static void constrain(int[][] k1, int[][] k2, float[][] x) {
    if (k1!=null && k2!=null) {
      int nc = k1.length;
      for (int ic=0; ic<nc; ++ic) {
        int nk = k1[ic].length;
        float sum = 0.0f;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = k1[ic][ik];
          int i2 = k2[ic][ik];
          sum += x[i2][i1];
        }
        float avg = sum/nk;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = k1[ic][ik];
          int i2 = k2[ic][ik];
          x[i2][i1] = avg;
        }
      }
    }
  }

  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
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

  private static void removeAverage(float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    float nh = (float)(n2);
    for (int i1=0; i1<n1; ++i1) {
      float sumx = 0.0f;
      for (int i2=0; i2<n2; ++i2)  
        sumx += x[i2][i1];
      float avgx = sumx/nh;
      for (int i2=0; i2<n2; ++i2) 
        x[i2][i1] -= avgx; 
    }
  }
  /* 
  private static void makeRhs(float w1, float[][] y) {
    zero(y);
    //add(y,w1*w1,y);
  }*/
  private static void makeRhs(float w1, float[][] y) {
    zero(y);
    int n1 = y[0].length;
    int n2 = y.length;
    float w1s = w1*w1;
    float ya = 0.5f*w1s;
    float yb = 0.5f*w1s;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }
  
  private static void applyLhs(
    float w1, float[][] wp, float[][] uz, float[][] ux, float[][] x, float[][] y) 
  {
    zero(y);
    int n1 = x[0].length;
    int n2 = x.length;
    float w1s = w1*w1;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        float uxi = ux[i2][i1];
        float uzi = uz[i2][i1];
        float wps = wpi*wpi;
        float uxs = uxi*uxi;
        float uzs = uzi*uzi;
        float d11 = wps*uxs+w1s;
        float d12 = -wps*uxi*uzi;
        float d22 = wps*uzs;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = 0.5f*(xa+xb);
        float x2 = 0.5f*(xa-xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  // Post-processing of computed shifts to ensure monotonic u1.
  private static void cleanShifts(float[][] r) {
    int n1 = r[0].length;
    int n2 = r.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if (r[i2][i1]<=r[i2][i1-1]-0.99f)
          r[i2][i1] = r[i2][i1-1]-0.99f;
      }
    }
  }

  private static void printStats(String s, float[][] a) {
    trace(s+": min="+min(a)+" max="+max(a));
  }

  private static void trace(String s) {
    System.out.println(s);
  }
  private void testSpd(String s, int n1, int n2, CgSolver.A a) {
    // symmetric: y'Ax = x'(A'y) = x'Ay
    // positive-semidefinite: x'Ax >= 0
    float[][] x = sub(randfloat(n1,n2),0.5f);
    float[][] y = sub(randfloat(n1,n2),0.5f);
    float[][] ax = zerofloat(n1,n2);
    float[][] ay = zerofloat(n1,n2);
    VecArrayFloat2 vx = new VecArrayFloat2(x);
    VecArrayFloat2 vy = new VecArrayFloat2(y);
    VecArrayFloat2 vax = new VecArrayFloat2(ax);
    VecArrayFloat2 vay = new VecArrayFloat2(ay);
    a.apply(vx,vax);
    a.apply(vy,vay);
    double yax = vy.dot(vax);
    double xay = vx.dot(vay);
    double xax = vx.dot(vax);
    double yay = vy.dot(vay);
    System.out.println(s+": yax="+yax+" xay="+xay+" should be equal");
    System.out.println(s+": xax="+xax+" yay="+yay+" should be positive");
  }
}
