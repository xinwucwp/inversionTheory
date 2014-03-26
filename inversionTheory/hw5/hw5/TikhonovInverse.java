package hw5;

import java.util.Random;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;
 // Xinming Wu, Colorado School of Mines
 // 2014.03.20

public class TikhonovInverse {

  public enum Method {
    SMALLEST,FLATTEST
  }

  public void setForNoise(double mean, double vari) {
    _mean = mean;
    _vari = vari;
  } 

  public void setForFlattest(double dx, double alphaS, double alphaX) {
    _dx = dx;
    _alphaS = alphaS;
    _alphaX = alphaX;
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

  public void computeData(double[][] G, double[] m, double[] d) {
    int n2 = G.length;
    int n1 = G[0].length;
    for (int i2=0; i2<n2; ++i2) { 
      double gm = 0.0;
      for (int i1=0; i1<n1; ++i1) 
        gm += G[i2][i1]*m[i1];
      d[i2] = gm;
    }
  } 

  public void addNoise(double[] d) {
    addNoise(_mean,_vari,d);
  }
  
  public void solveForPhi(
    Method method, double wd, Sampling sb, double[][] G, double[] d, 
    double[] phiD, double[] phiM) 
  {
    int n2 = G.length;
    int n1 = G[0].length;
    int nb = sb.getCount();
    double bb = sb.getFirst();
    double db = sb.getDelta();
    double[] b   = zerodouble(n1);
    double[] mc  = zerodouble(n1);
    double[] dp  = zerodouble(n2);
    double[][] A  = zerodouble(n1,n1);
    double[][] Gs = zerodouble(n1,n1);
    double[][] Wm = zerodouble(n1,n1);
    computeWmWm(_dx,_alphaS,_alphaX,Wm);
    computeGWWG(G,wd,Gs);
    applyForRhs(G,wd,d,b); 
    double beta = bb;
    for (int i=0; i<nb; i++, beta+=db) {
      A = copy(Gs);
      if(method==Method.SMALLEST)
        addDiagonal(beta,A);
      if(method==Method.FLATTEST)
        applyWmWm(beta,Wm,A);
      DMatrix dA  = new DMatrix(A);
      DMatrix dAi = dA.inverse();
      applyInverse(dAi.get(),b,mc);
      computeData(G,mc,dp);
      phiM[i] = computePhiM(mc);
      phiD[i] = computePhiD(wd,d,dp);
    }
  }
 
  
  public void computeGcv(
    Method method, double wd, Sampling sb, double[][] G, double[] d, 
    double[] gcv) 
  {
    int n2 = G.length;
    int n1 = G[0].length;
    int nb = sb.getCount();
    double bb = sb.getFirst();
    double db = sb.getDelta();
    double[] b   = zerodouble(n1);
    double[] mc  = zerodouble(n1);
    double[] dp  = zerodouble(n2);
    double[][] A  = zerodouble(n1,n1);
    double[][] Gs = zerodouble(n1,n1);
    double[][] Wm = zerodouble(n1,n1);
    double[][] Gw = zerodouble(n1,n1);
    computeWmWm(_dx,_alphaS,_alphaX,Wm);
    computeGWWG(G,wd,Gs);
    applyForRhs(G,wd,d,b); 
    double beta = bb;
    mul(G,1.0/(wd*wd),Gw);
    DMatrix dG  = new DMatrix(G);
    DMatrix dGT = dG.transpose();
    dGT = dGT.timesEquals(1.0/(wd*wd));
    for (int i=0; i<nb; i++, beta+=db) {
      A = copy(Gs);
      if(method==Method.SMALLEST)
        addDiagonal(beta,A);
      if(method==Method.FLATTEST)
        applyWmWm(beta,Wm,A);
      DMatrix dA  = new DMatrix(A);
      DMatrix dAi = dA.inverse();
      DMatrix dAG = dAi.times(dGT);
      DMatrix dGAG = dG.times(dAG);
      applyInverse(dAi.get(),b,mc);
      computeData(G,mc,dp);
      double numi = computePhiD(wd,d,dp);
      double deni = (double)n2-dGAG.trace();
      gcv[i] = (double)(n2*n2)*numi/(deni*deni);
    }
  }

  public void optimalSolver(
    Method method, double wd, double beta, double[][] G, double[] d, 
    double[] mc) 
  {
    int n2 = G.length;
    int n1 = G[0].length;
    double[] b   = zerodouble(n1);
    double[][] A  = zerodouble(n1,n1);
    double[][] Wm = zerodouble(n1,n1);
    computeWmWm(_dx,_alphaS,_alphaX,Wm);
    computeGWWG(G,wd,A);
    applyForRhs(G,wd,d,b); 
    if(method==Method.SMALLEST)
      addDiagonal(beta,A);
    if(method==Method.FLATTEST)
      applyWmWm(beta,Wm,A);
    DMatrix dA = new DMatrix(A);
    DMatrix dAi = dA.inverse();
    applyInverse(dAi.get(),b,mc);
  }

  public void computeCurvature(double[] phiD, double[] phiM, double[] cv) {
    int n = phiD.length;
    double[] dd1 = zerodouble(n);
    double[] dd2 = zerodouble(n);
    double[] dm1 = zerodouble(n);
    double[] dm2 = zerodouble(n);
    firstDerivative(log10(phiD),dd1);
    firstDerivative(log10(phiM),dm1);
    secondDerivative(log10(phiD),dd2);
    secondDerivative(log10(phiM),dm2);
    for (int i=0; i<n; ++i) {
      double dd1i = dd1[i];
      double dd2i = dd2[i];
      double dm1i = dm1[i];
      double dm2i = dm2[i];
      double numi = dm1i*dd2i-dd1i*dm2i;
      double deni = pow((dd1i*dd1i+dm1i*dm1i),1.5);
      cv[i] = -numi/deni;
    }
  } 

  public double betaFromPhiD(Sampling sb, double r, double[] phiD) {
    double beta = 0.0;
    int n = sb.getCount();
    double bb = sb.getFirst();
    double db = sb.getDelta();
    double dm = Double.POSITIVE_INFINITY; 
    for (int i=0; i<n; ++i) {
      double di = abs(phiD[i]-r); 
      if (di<dm) {
        dm = di;
        beta = bb+db*(double)i; 
      }
    }
    return beta; 
  }
  
  public double betaFromLcurve(Sampling sb, double[] cv) {
    int n = cv.length;
    double vm = Double.NEGATIVE_INFINITY; 
    int ind = 0;
    for (int i=0; i<n; ++i) {
      if (cv[i]>vm){
        ind = i; vm = cv[i];
      }
    }
    double fb = sb.getFirst();
    double db = sb.getDelta();
    double bi = fb+db*(double)ind;
    double ci = cv[ind];
    double cm = cv[ind-1];
    double cp = cv[ind+1];
    return parabolicPeak(db,bi,cm,ci,cp);
  } 

  public double betaFromGcv(Sampling sb, double[] gcv) {
    int n = gcv.length;
    double vm = Double.POSITIVE_INFINITY; 
    int ind = 0;
    for (int i=0; i<n; ++i) {
      if (gcv[i]<vm){
        ind = i; vm = gcv[i];
      }
    }
    double fb = sb.getFirst();
    double db = sb.getDelta();
    double bi = fb+db*(double)ind;
    double ci = gcv[ind];
    double cm = gcv[ind-1];
    double cp = gcv[ind+1];
    return parabolicPeak(db,bi,cm,ci,cp);
  } 

  public void computeWmWm(double dx, double alphaS, double alphaX,double[][] Wm) {
    int n = Wm.length;
    double wsi = alphaS*dx;
    double dx1 = alphaX/dx;
    double dx2 = 2.0*dx1;
    for (int i=1; i<n; ++i){
      Wm[i][i-1] = -dx1;
      Wm[i-1][i] = -dx1;
      Wm[i][i  ] =  dx2+wsi;
    }
    Wm[0  ][0  ] = dx1+wsi;
    Wm[n-1][n-1] = dx1+wsi;
  }   
  private void applyWmWm(double beta, double[][] Wm, double[][] x) {
    int n = Wm.length;
    x[0][0] += beta*Wm[0][0];
    for (int i=1; i<n; ++i){
      x[i][i  ] += beta*Wm[i  ][i];
      x[i][i-1] += beta*Wm[i][i-1];
      x[i-1][i] += beta*Wm[i-1][i];
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

  private void randomVector(double[] v) {
    int n = v.length;
    Random rm = new Random();
    double mean = 0.0;
    double vari = 50.0;
    for (int i=0; i<n; ++i) {
      double ni = mean+rm.nextGaussian()*vari;
      if (ni-mean<0.0) v[i] = -1.0;
      else v[i] = 1.0;
    } 
  }

  public void computeGWWG(double[][] G, double wd, double[][] A) {
    DMatrix dG = new DMatrix(G);
    DMatrix dA = new DMatrix(A);
    double dwd = 1.0/wd;
    dG = dG.timesEquals(dwd);
    dA = dG.transposeTimes(dG); 
    dA.get(A); 
  }
  
  private void applyForRhs(double[][] G, double wd, double[] d, double[] b) {
    int n2 = G.length;
    int n1 = G[0].length;
    double wds = 1.0/(wd*wd);
    for (int i1=0; i1<n1; ++i1) { 
      double sum = 0.0;
      for (int i2=0; i2<n2; ++i2)  
        sum += G[i2][i1]*d[i2];
      b[i1] = sum*wds;
    }
  } 

  private void applyInverse(double[][] A, double[] b, double[] x) {
    int n2 = A.length;
    int n1 = A[0].length;
    for (int i2=0; i2<n2; ++i2) {
      double sum = 0.0;
      for (int i1=0; i1<n1; ++i1)
        sum += A[i2][i1]*b[i1];
      x[i2] = sum;
    }
  }
  private double computePhiD(double wd, double[] d, double[] dp) {
    int n = d.length;
    double sum = 0.0;
    double wds = 1.0/(wd*wd);
    for (int i=0; i<n; ++i) {
      double dd = d[i]-dp[i];
      sum += dd*dd; 
    }
    return sum*wds;
  }

  private double computePhiM(double[] mc) {
    int n = mc.length;
    double sum = 0.0;
    for (int i=0; i<n; ++i) {
      double mi = mc[i];
      double dm = mi*mi;
      sum += dm;
    }
    return sum;
  }

  private void addDiagonal(double a, double[][] A) {
    int n2 = A.length;
    int n1 = A[0].length;
    int n = min(n1,n2);
    for (int i=0; i<n; ++i)
      A[i][i] += a;
  }
  
  private void firstDerivative(double[] x, double[] y) {
    int n = x.length;
    for (int i=1; i<n-1; ++i) 
      y[i] = 0.5*(x[i+1]-x[i-1]);
    y[0  ]=y[1  ];
    y[n-1]=y[n-2];
  }
 
  private void secondDerivative(double[] x, double[] y) {
    int n = x.length;
    for (int i=1; i<n-1; ++i) 
      y[i] = x[i+1]+x[i-1]-2.0*x[i];
    y[0  ]=y[1  ];
    y[n-1]=y[n-2];
  }

  private double innerProduct(double[] u, double[] v) {
    int n = u.length;
    double p = 0.0;
    for (int i=0; i<n; ++i)
      p += u[i]*v[i];
    return p;
  }
  private double parabolicPeak(double dz, double zi, double um, double ui, double up) {
    double a = um-up;
    double b = 2.0f*(um+up)-4.0f*ui;
    return (zi+dz*a/b);
  }
  
  private double _mean = 0.0;
  private double _vari = 0.05;
  private double _dx = 1.0;
  private double _alphaX = 1.0;
  private double _alphaS = 1.0;
}
