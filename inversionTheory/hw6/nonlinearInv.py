import sys
from math import *
from java.io import *
from java.nio import *
from java.awt import *
from java.lang import *
from java.util import *
from javax.swing import *
from java.util.Random import *

from edu.mines.jtk.io import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from hw6 import*
#############################################################################
# parameters
nd = 48
nm = 120
dd,dm=250.0,100.0
fd,fm=-6000.0,-6000.0
sd = Sampling(nd,dd,fd)
sm = Sampling(nm,dm,fm)
alphaS,alphaX=0.00002,1.0
dSigma = 0.2
zo = fillfloat(-2.0,nd)
zt = zerofloat(nm)
xc = zerofloat(nm)
dk = zerofloat(nd)
mk = fillfloat(200.0,nm)
rd = Random()
for i in range(120):
  xc[i] = sm.getValue(i)+50.0
wd,dc=100.0,-1.0
#pngDir = "../../../HW6/images/"
#dataDir = "../../../HW6/data/"
pngDir = "./images/"
dataDir = "./data/"
#############################################################################

def main(args):
  #beta =[1.0e-1,1.0e-2,1.0e-3,1.0e-4,1.0e-5,1.0e-6]
  beta = [1.0e-4]
  da = readData(48,4,'hw06')
  xo = da[0]
  do = da[2]
  ni = NonlinearInverter()
  ni.setSmoothing(12.0)
  ni.setForCG(0.01,200)
  ni.setForward(xo,zo,xc,zt,wd,dc)  
  ni.setPerturb(35.0)
  ni.setWmWd(alphaX,alphaS,dSigma)
  ms = zerofloat(nm,len(beta))
  for i in range(len(beta)):
    print "beta="
    print beta[i]
    print "======================================"
    mk = fillfloat(200.0,nm)
    ni.inverter(sm,beta[i],do,mk)
    copy(mk,ms[i]) 
    plot1D(sm,mk,"Model","x(m)",title="model")
  #writeData("models",ms)
  ms = readData(nm,6,"models")
  cs = [Color.blue,Color.cyan,Color.green,Color.yellow,Color.orange,Color.red]
  plotMultiple1D(sm,ms,cs,"Recovered models","x(m)",png="models")
  plot1D(sd,do,"Observed data","x(m)",title="data")
  '''
  phi=[321130.03,4815.872,150.72491,30.443737,19.774208]
  sh = Sampling(5,1.0,0.0)
  plot1D(sh,phi,"Object function","Gauss-Newton iterates",title="gs")
  '''
##################################################################
# plots
jet = ColorMap.JET
gray = ColorMap.GRAY
def plot1D(s,x1,vlabel,hlabel,x2=None,title=None):
  sp = SimplePlot()
  pv = sp.addPoints(s,x1)
  if x2:
   pv1 = sp.addPoints(s,x1)
   pv2 = sp.addPoints(s,x2)
   pv1.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
   pv2.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
   pv1.setLineColor(Color.blue)
   pv1.setMarkColor(Color.blue)
   pv2.setLineColor(Color.red)
   pv2.setMarkColor(Color.red)
   pv1.setMarkSize(8.0)
   pv2.setMarkSize(8.0)
   pv1.setLineWidth(2.0)
   pv2.setLineWidth(2.0)
  else:
    #pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setLineColor(Color.black)
    pv.setLineWidth(2.0)
    #pv.setMarkSize(8.0)
  sp.setSize(1000+80,400)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  sp.setFontSize(18)
  if title:
    sp.paintToPng(720,3.3,pngDir+title+".png")

def plotMultiple1D(s,v,cs,vlabel,hlabel,png=None):
  sp = SimplePlot()
  n = len(cs)
  for i in range(n):
    pv = sp.addPoints(s,v[i])
    pv.setLineStyle(PointsView.Line.SOLID)
    pv.setLineColor(cs[i])
    pv.setLineWidth(2.0)
    pv.setMarkSize(8.0)
  sp.setVLimits(0,1500)
  sp.setHLimits(-6000,6000)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  sp.setSize(1000+80,400)
  if png:
    sp.paintToPng(720,3.3,pngDir+png+".png")

def plot2D(s1,s2,x,wd,vlabel,hlabel,cbar,cmap=jet,png=None,orient=True):
  n1 = s1.getCount()
  n2 = s2.getCount()
  if orient:
    sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
    pv = sp.addPixels(s1,s2,x)
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  else:
    sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
    pv = sp.addPixels(s1,s2,x)
    pv.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cmap)
  dims = [n1,n2]
  sp.setSize(wd+80,400)
  sp.setVInterval(2)
  sp.setHInterval(10)
  sp.setVInterval(10)
  cb=sp.addColorBar()
  cb.setLabel(cbar)
  cb.setWidthMinimum(80)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  sp.setFontSize(18)
  if png:
    sp.paintToPng(720,3.3,pngDir+png+".png")
#############################################################################
def readData(n1,n2,name):
  fileName = dataDir+name+".dat"
  data = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.BIG_ENDIAN)
  ais.readFloats(data)
  ais.close()
  return data

def readData1D(n1,name):
  fileName = dataDir+name+".dat"
  data = zerodouble(n1)
  ais = ArrayInputStream(fileName)
  ais.readDoubles(data)
  ais.close()
  return data

def writeData(name,data):
  fileName = dataDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(data)
  aos.close()
  return data


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
