package hv;

import java.util.ArrayList;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import he.*;

/**
 * Set up control points for generating a horizon volume
 * Method rearrange:
 *   Compute the average depth of all the control points of each set, and set 
 *   the point with the average depth as the reference point for that set.
 * Method extend:
 *   Optional but usually useful: extend scattered control points to a control 
 *   surfaces (each set produces one surface) by using the "HorizonExtractorC" method.
 *
 * @author Xinming Wu
 * @version 2014.03.13
 */

public class SetupConstraints {
  public SetupConstraints(float[] k1, float[] k2, float[] k3) {
    _k1 = k1;
    _k2 = k2;
    _k3 = k3;
  }
 
  public void setForExtension (float sigma1, float sigma2, float scale) {
    _scale  = scale ;
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  } 

  public float[][] rearrange() {
    int ai = 0;
    int np = _k1.length;
    float[][] k = zerofloat(np,4);
    float k1a = sum(_k1)/(float)np;
    float min = Float.POSITIVE_INFINITY;
    for (int ip=0; ip<np; ip++) {
      int k1i = round(_k1[ip]);
      k[0][ip] = (float)k1i; 
      k[1][ip] = _k2[ip];
      k[2][ip] = _k3[ip];
      k[3][ip] = _k1[ip]-k1i;
      float dk = abs(_k1[ip]-k1a);
      if (dk<min){min = dk; ai=ip;}
    }
    float t0 = k[0][0];
    float t1 = k[1][0];
    float t2 = k[2][0];
    float t3 = k[3][0];
    k[0][0] = k[0][ai];  k[0][ai] = t0;
    k[1][0] = k[1][ai];  k[1][ai] = t1;
    k[2][0] = k[2][ai];  k[2][ai] = t2;
    k[3][0] = k[3][ai];  k[3][ai] = t3;
    return k;
  }

  public float[][] extend(
    int w2, int w3, float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] u, float[][][] wp) {
    int np = _k1.length;
    int n3  = p2.length;
    int n2  = p2[0].length;
    int n1  = p2[0][0].length;
    int ib2 = (int)min(_k2)-w2; if(ib2<0   ) {ib2=0;   }
    int ib3 = (int)min(_k3)-w3; if(ib3<0   ) {ib3=0;   }
    int ie2 = (int)max(_k2)+w2; if(ie2>n2-1) {ie2=n2-1;} 
    int ie3 = (int)max(_k3)+w3; if(ie3>n3-1) {ie3=n3-1;} 
    int n2m = ie2-ib2+1;
    int n3m = ie3-ib3+1;
    float[][][] um  = copy(n1,n2m,n3m,0,ib2,ib3,u );
    float[][][] p2m = copy(n1,n2m,n3m,0,ib2,ib3,p2);
    float[][][] p3m = copy(n1,n2m,n3m,0,ib2,ib3,p3);
    float[][][] epm = copy(n1,n2m,n3m,0,ib2,ib3,ep);
    sub(_k2,ib2,_k2);
    sub(_k3,ib3,_k3);
    SurfaceExtractorC se = new SurfaceExtractorC();
    se.setSmoothings(_sigma1,_sigma2);
    se.setWeights(_scale);
    float lmt = n1-1.f;
    float[][] surf = se.surfaceInitialization(n2m,n3m,lmt,_k1,_k2,_k3);
    se.surfaceUpdateFromSlopes(epm,p2m,p3m,_k1,_k2,_k3,surf);
    int ci = 0;
    int ai = 0;
    float[][] k = zerofloat(n2m*n3m,4);
    float min = Float.POSITIVE_INFINITY;
    float avg = heightAvg(lmt,surf);
    for (int i3=0; i3<n3m; i3++) {
      for (int i2=0; i2<n2m; i2++) {
        int i1m = round(surf[i3][i2]);
        int i2m = i2+ib2;     
        int i3m = i3+ib3;     
        if (i1m<n1-1){
          k[0][ci] = i1m; 
          k[1][ci] = i2m;
          k[2][ci] = i3m;
          k[3][ci] = surf[i3][i2]-i1m;
        float df = abs(surf[i3][i2]-avg);
        if (df<min){min = df; ai=ci;}
          ci++;
        }
      }
    }
    float t0 = k[0][0];
    float t1 = k[1][0];
    float t2 = k[2][0];
    float t3 = k[3][0];
    k[0][0] = k[0][ai]; k[0][ai] = t0;
    k[1][0] = k[1][ai]; k[1][ai] = t1;
    k[2][0] = k[2][ai]; k[2][ai] = t2;
    k[3][0] = k[3][ai]; k[3][ai] = t3;
    return copy(ci,4,0,0,k);
  } 
  private void removeFaults(float[][][] f) {
    int n = _k1.length;
    ArrayList<float[]> kk = new ArrayList<float[]>();
    for (int i=0; i<n; ++i) {
      int i1 = (int)_k1[i];
      int i2 = (int)_k2[i];
      int i3 = (int)_k3[i];
      float fi = f[i3][i2][i1];
      if (fi>0.0f) {
        float[] kki = {_k1[i],_k2[i],_k3[i]};
        kk.add(kki);
      }
    }
    int np = kk.size();
    _k1 = new float[np];
    _k2 = new float[np];
    _k3 = new float[np];
    for (int ip=0; ip<np; ++ip) {
      _k1[ip] = kk.get(ip)[0];
      _k2[ip] = kk.get(ip)[1];
      _k3[ip] = kk.get(ip)[2];
    }
  } 
  private float heightAvg(float lmt, float[][] x) {
    int n2=x.length;
    int n1=x[0].length;
    float sum = 0.0f;
    float ci = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = x[i2][i1];
        if(xi<lmt) {sum +=xi; ci +=1.0f;}
      }
    }
    return sum/ci;
  }
  private float[] _k1;
  private float[] _k2;
  private float[] _k3;
  private float _scale  = 0.0f;
  private float _sigma1 = 12.0f;
  private float _sigma2 = 12.0f;
}
