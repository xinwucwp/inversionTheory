import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from hv import *
from fault import *
from util import *

seismicDir = "../../data/"
ffile = "fs"
p2file = "p2s"
wpfile = "wps"
s1 = Sampling(20,1.0,0.0)
s2 = Sampling(50,1.0,0.0)
sp1 = Sampling(20,1.0,0.000)
sp2 = Sampling(50,1.0,0.000)
dx,dt = 1.0,1.0
fx,ft = 0.0,0.0
n1,n2 = s1.count,s2.count
d1,d2 = s1.delta,s2.delta
f1,f2 = sp1.first,sp2.first
def main(args):
  flatten()
  #horizonExtract()
def flatten():
  f = readImage(n1,n2,ffile)
  p2 = readImage(n1,n2,p2file)
  wp = readImage(n1,n2,wpfile)
  k1s = [[12,10]]
  k2s = [[41,10]]
  #f = FakeData.seismic2d2011A(n1,n2,30)
  plot(sp1,sp2,f,clab="Amplitude",vlabel="Time (s)",hlabel="Inline (km)", 
       cmin=-2.5, cmax=2.5,png="f")
  plot(sp1,sp2,f,c=(k1s,k2s),clab="Amplitude",vlabel="Time (s)",hlabel="Inline (km)",
       cmin=-2.5,cmax=2.5,png="c")
  fl = Flattener2C()
  fl.setWeight1(0.05)
  fl.setIterations(0.01,1000)
  fl.setSmoothings(4.0,8.0)
  for cs in [[k1s,k2s]]:
    if cs:
      nc = len(cs[0])
      psuffix = str(nc)
      tsuffix = " ("+str(nc)+" constraints)"
    else:
      psuffix = "0"
      tsuffix = " (no constraints)"
    fm = fl.getMappingsFromSlopes(sp1,sp2,p2,wp,cs)
    g = fm.flatten(f)
    h = fm.unflatten(g)
    s = fm.getShiftsS()
    u = fm.u1
    mul(u,sp1.delta,u)
    add(u,sp1.first,u)
    rc = readImage(n1,n2,"r1c")
    mrc = readImage(n1,n2,"Mr1c")
    sc0= readImage(n1,n2,"sc0")
    scf= readImage(n1,n2,"scf")
    scfm= readImage(n1,n2,"scfm")
    print "shifts:"
    print sc0[10][10]
    print sc0[41][12]
    print u[10][10]
    print u[41][12]
    #plot(sp1,sp2,rc,c=cs,cmap=jet,cmin=min(mrc),cmax=max(mrc),title="r",png="rc"+psuffix)
    #plot(sp1,sp2,mrc,c=cs,cmap=jet,cmin=min(mrc),cmax=max(mrc),title="r",png="mrc"+psuffix)
    plot(sp1,sp2,sc0,c=cs,cmap=jet,cmin=-3.2,cmax=3,title="r",png="initialShifts"+psuffix)
    plot(sp1,sp2,scf,c=cs,cmap=jet,cmin=-3.2,cmax=3.8,title="r",png="updateShifts"+psuffix)
    plot(sp1,sp2,u,c=cs,cmap=jet,clab="RGT",vlabel="Time (s)",hlabel="Inline (km)",png="rgt"+psuffix)
    #plot(sp1,sp2,rc,c=cs,cmap=jet,clab="RGT",vlabel="Time (s)",hlabel="Inline (km)",title="r",png="rgt"+psuffix)
    #plot(sp1,sp2,mrc,c=cs,cmap=jet,clab="RGT",vlabel="Time (s)",hlabel="Inline (km)",title="Pr",png="rgt"+psuffix)
    if cs:
      nc = len(cs[0])
      for ic in range(nc):
        np = len(cs[0][ic])
        for ip in range(np):
          i1 = cs[0][ic][ip]
          i2 = cs[1][ic][ip]
          cs[0][ic][ip]=i1+s[i2][i1]
    #plot(sp1,sp2,f,u=u,clab="RGT",vlabel="Time (s)",hlabel="Inline (km)",png="hs"+psuffix)
    plot(sp1,sp2,g,c=cs,clab="Amplitude",vlabel="RGT",
         hlabel="Inline (km)",cmin=-2.5,cmax=2.5,png="g"+psuffix)
    writeImage("rgt"+psuffix,u)
    writeImage("shift"+psuffix,s)
    writeImage("hv"+psuffix,fm.x1)
    writeImage("g"+psuffix,g)
    #plot(sp1,sp2,h,title="Unflattened"+tsuffix,png="h"+psuffix)
    #plot(s1,s2,s,cmap=jet,title="Shifts"+tsuffix,png="s"+psuffix)
    print "average shift =",sum(s)/(n1*n2),"samples"

def flip2(f):
  n2 = len(f)
  for i2 in range(n2/2):
    fi2 = f[i2]
    f[i2] = f[n2-1-i2]
    f[n2-1-i2] = fi2

def faults(f):
  fs = FaultSemblance()
  p2 = fs.slopes(f)
  snd = fs.semblanceNumDen(p2,f)
  fs = FaultScanner2(20,snd,FaultScanner2.Smoother.SHEAR)
  (fl,ft) = fs.scan(-20,20)
  #(fl,ft) = fs.thin((fl,ft))
  return fl
def gain(n1,n2,x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
pngDir = None
pngDir = "../../png/hv/2d/slides/"

def plot(s1,s2,x,u=None,c=None,cmap=ColorMap.GRAY,clab=None,vlabel=None,hlabel=None,
  cmin=0,cmax=0,title=None,png=None):
  #sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.NONE)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if c: 
    k1s,k2s=c
    tp1 = PointsView.Mark.HOLLOW_CIRCLE
    pv = sp.addPoints(k1s,k2s)
    pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(tp1)
    pv.setMarkColor(Color.BLUE)
    pv.setMarkSize(45)
    pv.setLineWidth(6)
    #sp.add(pv)
  # sp.setFontSizeForSlide(0.5,0.9,16.0/9.0)
  frame = PlotFrame(sp) 
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setSize(1000+120,400)
  frame.setVisible(True)
  if pngDir and png:
    frame.paintToPng(360,3.33,pngDir+png+".png")

#############################################################################
# utilities

def readImage(n1,n2,name):
  fileName = seismicDir+name+".dat"
  #n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
