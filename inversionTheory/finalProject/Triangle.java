package hv;

import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
public class Triangle {
  public float[] trianglesForSurface(float[][] surf, float top, float down) {
    int n1 = surf[0].length;
    int n2 = surf.length;
    float[] xyz = new float[n1*n2*18];
    int i = 0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float surfi = surf[i2][i1];
        if (surfi<down && surfi>top && i1<n1-1 && i2<n2-1) {
          xyz[i   ] = (float) i2;
          xyz[i+1 ] = (float) i1;
          xyz[i+2 ] = surf[i2][i1];
          xyz[i+3 ] = (float) i2+1;
          xyz[i+4 ] = (float) i1;
          xyz[i+5 ] = surf[i2+1][i1];
          xyz[i+6 ] = (float) i2+1;
          xyz[i+7 ] = (float) i1+1;
          xyz[i+8 ] = surf[i2+1][i1+1];

          xyz[i+9 ] = (float) i2;
          xyz[i+10] = (float) i1;
          xyz[i+11] = surf[i2][i1];
          xyz[i+12] = (float) i2+1;
          xyz[i+13] = (float) i1+1;
          xyz[i+14] = surf[i2+1][i1+1];
          xyz[i+15] = (float) i2;
          xyz[i+16] = (float) i1+1;
          xyz[i+17] = surf[i2][i1+1];
          i = i+18;
        } 
      } 
    }
    return copy(i,0,xyz);
  }

  public float[] trianglesForS3(int n1, int n2, int i3, float[] top) {
    float[] xyz = new float[n1*n2*18];
    int i = 0;
    for (int i2=0; i2<n2-1; i2++) {
      for (int i1=round(top[i2]); i1<n1-1; i1++) {
        xyz[i   ] = (float) i3;
	xyz[i+1 ] = (float) i2;
	xyz[i+2 ] = (float) i1;
	xyz[i+3 ] = (float) i3;
	xyz[i+4 ] = (float) i2+1;
	xyz[i+5 ] = (float) i1;
	xyz[i+6 ] = (float) i3;
	xyz[i+7 ] = (float) i2+1;
	xyz[i+8 ] = (float) i1+1;

	xyz[i+9 ] = (float) i3;
	xyz[i+10] = (float) i2;
	xyz[i+11] = (float) i1;
	xyz[i+12] = (float) i3;
	xyz[i+13] = (float) i2+1;
	xyz[i+14] = (float) i1+1;
	xyz[i+15] = (float) i3;
	xyz[i+16] = (float) i2;
	xyz[i+17] = (float) i1+1;
	i = i+18;
      }
    }
    return copy(i,0,xyz);
  }
}
