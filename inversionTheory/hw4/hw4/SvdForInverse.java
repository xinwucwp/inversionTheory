package hw4;

import java.util.Random;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;
 // Xinming Wu, Colorado School of Mines
 // 2014.03.15

public class SvdForInverse {

  public void setForNoise(double mean, double vari) {
    _mean = mean;
    _vari = vari;
  } 
  public void constructG(Sampling sx, double[][] G) {
    int nk = G.length;
    int nx = G[0].length; 
    double fx = sx.getFirst();
    double dx = sx.getDelta();
    double pi = Math.PI;
    for (int ik=1; ik<=nk; ++ik) { 
      for (int ix=0; ix<nx; ++ix) { 
        double xi = fx+dx*ix;
        double ki = (double)ik-1.0;
        double x1 = 0.25*ki*xi;
        double x2 = 2.00*pi*x1;
        G[ik-1][ix] = cos(x2)*exp(-x1)*dx; 
      }
    }
  }
  public void normalizeG(double sigma, double[][] G, double[][] NG) { 
    int n2 = G.length;
    int n1 = G[0].length;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        NG[i2][i1] = G[i2][i1]/sigma; 
  }

  public void computeData(double[][] G, double[] m, double[] dt, double[] dn) {
    int n2 = G.length;
    int n1 = G[0].length;
    for (int i2=0; i2<n2; ++i2) { 
      double gm = 0.0;
      for (int i1=0; i1<n1; ++i1) 
        gm += G[i2][i1]*m[i1];
      dt[i2] = gm;
      dn[i2] = gm;
    }
    addNoise(_mean,_vari,dn);
  } 

  public void normalizeData(double sigma, double[] d, double[] nd) {
    int n = d.length;
    for (int i=0; i<n; ++i)
      nd[i] = d[i]/sigma;
  }

  public void svdForG(double[][] G, double[] s, double[][] u, double[][] v) {
    DMatrix dG = new DMatrix(G);
    DMatrixSvd svd = dG.svd();
    DMatrix dU = svd.getU();
    DMatrix dV = svd.getV();
    double[] S = svd.getSingularValues();
    copy(S,s);
    dU.get(u);
    dV.get(v);
  } 
   
  public void modelCoefficients(
    double[] s,  double[][] u, double[] d,
    double[] rd, double[] ra) 
  {
    int nk = s.length;
    for (int i1=0; i1<nk; ++i1) { 
      double vi = 0.0;
      for (int i2=0; i2<nk; ++i2) { 
        vi += u[i2][i1]*d[i2];
      }
      rd[i1] = vi;
      ra[i1] = vi/s[i1];
    }
  }
  
  public void computePhiD(double[] rd, double[] phiD) {
    int n = rd.length;
    zero(phiD);
    double sum = 0.0; 
    for (int i=0; i<n; ++i)
      sum += rd[i]*rd[i]; 
    for (int q=0; q<n; ++q) {
      double sumq = 0.0;
      for (int i=0; i<q; ++i)
        sumq += rd[i]*rd[i]; 
      phiD[q] = sum-sumq;
    }
  }

  public void computePhiM(double[] ra, double[] phiM) {
    int n = ra.length;
    zero(phiM);
    for (int q=0; q<n; ++q) {
      double sum = 0.0;
      for (int i=0; i<q; ++i)
        sum += ra[i]*ra[i];
      phiM[q] = sum;
    }
  }

  public int findQ(double[] phiD, int N) {
    double Nd = (double) N; 
    int n = phiD.length;
    int q = 0;
    double dMin = Double.POSITIVE_INFINITY;
    for (int i=0; i<n; ++i) {
      double di = abs(phiD[i]-Nd);
      if (di<dMin) {
        q = i;
        dMin = di;
      }
    }
    return q;
  }

  public void modelConstruct(
    int[] t, double[] s,  double[][] u,  double[][] v, 
    double[] d, double[] mc) 
  {
    int nk = s.length;
    int nx = mc.length;
    zero(mc);
    int i1b = t[0];
    int i1e = t[1];
    for (int i1=i1b; i1<i1e; ++i1) { 
      double vi = 0.0;
      for (int i2=0; i2<nk; ++i2) { 
        vi += u[i2][i1]*d[i2];
      }
      double rai = vi/s[i1];
      for (int ix=0; ix<nx; ++ix)
        mc[ix] += v[ix][i1]*rai;
    }
  }

  //////////////////////////////////////////////////////////////
  // private
  private void addNoise(double mean, double vari, double[] d) {
    int n = d.length;
    Random rm = new Random();
    for (int i=0; i<n; ++i) {
      double ni = mean+rm.nextGaussian()*vari;
      d[i] += ni;
    } 
  }

  private double _mean = 0.0;
  private double _vari = 0.05;
}
