package hw2;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.la.*;
import static edu.mines.jtk.util.ArrayMath.*;
 // Xinming Wu, Colorado School of Mines
 // 2014.02.09

public class MinNormModel {

  // constructed and true model
  public void trueAndConstructed(
    Sampling s,  double[] alpha, 
    double[] mt, double[] mc) 
  {
    int nm = s.getCount();
    int na = alpha.length; 
    double fx = s.getFirst();
    double dx = s.getDelta(); 
    double pi = Math.PI;
    for (int im=0; im<nm; ++im) {
      double xi = fx + im*dx;
      mt[im] = 1.0 - 0.5*cos(2.0*pi*xi);
      for (int ia=1; ia<=na; ++ia) {
        double ci = 1.0 - (double)ia;
        mc[im] += alpha[ia-1]*exp(ci*xi);
      } 
    }
  }  

  // produce accurate data
  public void applyForData(double[] d) {
    int n = d.length;
    double pi = Math.PI;
    d[0] = 1.0f;
    for (int i=2; i<=n; ++i) {
      double ci = 1.0-(double)i;
      double ei = Math.exp(ci) - 1.0;
      d[i-1] = ei/ci - ci*ei/(2.0*ci*ci+8.0*pi*pi);  
    }
  }

  // compute coefficients alpha 
  public void coefficients(
    double[] d, double[] e, 
    double[][] v, double[] alpha) {
    int n2 = v.length; 
    int n1 = v[0].length; 
    double[] mc = new double[n1];
    for(int i2=0; i2<n2; ++i2) {
      double[] vi2 = v[i2];
      double rd = dot(vi2,d);
      double ei = e[i2];     
      for(int i1=0; i1<n1; ++i1) {
        double vi = vi2[i1];
        alpha[i1] += vi*rd/ei; 
      }
    }
  } 
  // eigenvalue decomposition of the inner product matrix
  public void eigenDecompo(double[][] g, double[] e, double[][] v) {
    int n2 = g.length;
    int n1 = g[0].length;
    DMatrix dmg = new DMatrix(g);
    DMatrixEvd evd = new DMatrixEvd(dmg);
    DMatrix dmv = evd.getV();
    double[] de = evd.getRealEigenvalues();
    double[][] dv = dmv.get();
    for (int i2=0; i2<n2; ++i2) {
      int i2r = n2-1-i2;
      e[i2r] = de[i2];
      for (int i1=0; i1<n1; ++i1) {
        v[i2r][i1] = dv[i1][i2];
      }
    }
  } 

  // compute the inner product of the kernel functions
  public void kernelMatrix(double[][] g) {
    int n2 = g.length;
    int n1 = g.length;
    for (int i2=1; i2<=n2; ++i2) {
      for (int i1=1; i1<=n1; ++i1) {
        double ci1 = 1.0 - (double)i1;
        double ci2 = 1.0 - (double)i2;
        double c12 = ci1+ci2;
        if(c12==0.0) {
          g[i2-1][i1-1] = 1.0;
        } else {
          double e12 = Math.exp(c12) - 1.0;
          g[i2-1][i1-1] = e12/c12;
        }
      }
    }
  }
  
  private double dot(double[] a, double[] b) {
    int n = a.length;
    double v = 0.0;
    for (int i=0; i<n; ++i) {
      v += a[i]*b[i]; 
    }
    return v;
  }

}
