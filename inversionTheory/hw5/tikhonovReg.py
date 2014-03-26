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

from hw5 import *

#############################################################################
# parameters
nk = 20
nx = 100
nb = 1000
bb,db = 0.01,0.01
dx = 0.01
fx = dx*0.5
s1 = Sampling(nx,1,0)
s2 = Sampling(nk,1,0)
sb = Sampling(nb,db,bb)
sx = Sampling(nx,dx,fx)
xb = zerodouble(nb)
for i in range(nb):
  xb[i] = log10(bb+i*db)
mean,sigma = 0.0,0.05
alphaS,alphaX=1.0,0.25
#pngDir = "../../../HW5/images/"
#dataDir = "../../../HW5/data/"
pngDir = "./images/"
dataDir = "./data/"
#############################################################################

def main(args):
  m = readData(nx,2,"model")
  mt = m[1]
  ti = TikhonovInverse()
  g = zerodouble(nx,nk) 
  ti.constructG(sx,g)
  showWmWm()
  computeData(g,mt)
  applySmallest(g,mt)
  applyFlattest(g,mt)
def showWmWm():
  ti = TikhonovInverse()
  w = zerodouble(nx,nx) 
  ti.computeWmWm(dx,alphaS,alphaX,w)
  plot2D(s1,s1,w,600,"nx","nx","Wm_ij",cmap=jet,png="Wm",orient=False)

  #########################data
def computeData(g,mt):
  dt = zerodouble(nk)
  ti = TikhonovInverse()
  ti.setForNoise(mean,sigma)
  ti.computeData(g,mt,dt)
  #copy(dt,dn)
  #ti.addNoise(dn)
  #writeData("dn3",dn)
  dn=readData1D(nk,"dn2")
  plot1D(s2,dt,"Data","nk",x2=dn,title="data")
  #########################Smallest
def applySmallest(g,mt):
  dn=readData1D(nk,"dn2")
  ti = TikhonovInverse()
  cv = zerodouble(nb)
  mc1 = zerodouble(nx)
  mc2 = zerodouble(nx)
  mc3 = zerodouble(nx)
  phiD = zerodouble(nb)
  phiM = zerodouble(nb)
  ti.solveForPhi(TikhonovInverse.Method.SMALLEST,sigma,sb,g,dn,phiD,phiM)
  ti.computeCurvature(phiD,phiM,cv)
  beta1 = ti.betaFromPhiD(sb,20.0,phiD)
  beta2 = ti.betaFromLcurve(sb,cv)
  print "beta1=";print beta1
  print "beta2=";print beta2
  ti.optimalSolver(TikhonovInverse.Method.SMALLEST,sigma,beta1,g,dn,mc1)
  #plot1D(sx,mt,"SmallestModel","x",x2=mc,title="smallestModel1")
  ti.optimalSolver(TikhonovInverse.Method.SMALLEST,sigma,beta2,g,dn,mc2)
  #plot1D(sx,mt,"SmallestModel","x",x2=mc,title="smallestModel2")
  plot1D(xb,log10(phiD),"log10(phiD)","Beta",title="smallestPhiD")
  plot1D(xb,log10(phiM),"log10(phiM)","Beta",title="smallestPhiM")
  plotLogScale(xb,cv,"curvature","Beta",title="smallestCurvature")
  plotLogScale(log10(phiM),log10(phiD),"Log10(phiD)","Log10(phiM)",title="smallestT")
  #########################Gcv for smallest
  gcv = zerodouble(nb)
  ti.setForFlattest(dx,alphaS,alphaX)
  ti.computeGcv(TikhonovInverse.Method.SMALLEST,sigma,sb,g,dn,gcv)
  plotLogScale(xb,log10(gcv),"GCV","Beta",title="smallestGcv")
  beta3 = ti.betaFromGcv(sb,gcv)
  print "beta3=";print beta3
  ti.optimalSolver(TikhonovInverse.Method.SMALLEST,sigma,beta3,g,dn,mc3)
  ms = zerodouble(4,nx) 
  ms[0],ms[1] = mt, mc1;
  ms[2],ms[3] = mc2,mc3;
  cs = [Color.blue,Color.red,Color.green,Color.orange]
  plotMultiple1D(sx,ms,cs,"Models","x",png="smallestModels")

  #########################Flattest
def applyFlattest(g,mt):
  dn=readData1D(nk,"dn2")
  ti = TikhonovInverse()
  cv = zerodouble(nb)
  mc1 = zerodouble(nx)
  mc2 = zerodouble(nx)
  mc3 = zerodouble(nx)
  phiD = zerodouble(nb)
  phiM = zerodouble(nb)
  ti.setForFlattest(dx,alphaS,alphaX)
  ti.solveForPhi(TikhonovInverse.Method.FLATTEST,sigma,sb,g,dn,phiD,phiM)
  ti.computeCurvature(phiD,phiM,cv)
  beta1 = ti.betaFromPhiD(sb,20.0,phiD)
  beta2 = ti.betaFromLcurve(sb,cv)
  print "beta1=";print beta1
  print "beta2=";print beta2
  ti.setForFlattest(dx,alphaS,alphaX)
  ti.optimalSolver(TikhonovInverse.Method.FLATTEST,sigma,beta1,g,dn,mc1)
  #plot1D(sx,mt,"FlattestModel","x",x2=mc1,title="flattestModel1")
  ti.optimalSolver(TikhonovInverse.Method.FLATTEST,sigma,beta2,g,dn,mc2)
  #plot1D(sx,mt,"FlattestModel","x",x2=mc2,title="flattestModel2")
  plot1D(xb,log10(phiD),"log10(phiD)","Beta",title="flattestPhiD")
  plot1D(xb,log10(phiM),"log10(phiM)","Beta",title="flattestPhiM")
  plotLogScale(xb,cv,"Curvature","Beta",title="flattestCurvature")
  plotLogScale(log10(phiM),log10(phiD),"Log10(phiD)","Log10(phiM)",title="flattestT")
  #########################Gcv for flattest
  gcv = zerodouble(nb)
  ti.setForFlattest(dx,alphaS,alphaX)
  ti.computeGcv(TikhonovInverse.Method.FLATTEST,sigma,sb,g,dn,gcv)
  plotLogScale(xb,log10(gcv),"GCV","Beta",title="flattestGcv")
  beta3 = ti.betaFromGcv(sb,gcv)
  print "beta3=";print beta3
  ti.optimalSolver(TikhonovInverse.Method.FLATTEST,sigma,beta3,g,dn,mc3)
  #plot1D(sx,mt,"FlattestModel","x",x2=mc3,title="flattestModel3")
  ms = zerodouble(4,nx) 
  ms[0],ms[1] = mt, mc1;
  ms[2],ms[3] = mc2,mc3;
  cs = [Color.blue,Color.red,Color.green,Color.orange]
  plotMultiple1D(sx,ms,cs,"Models","x",png="flattestModels")
  
##################################################################
# plots
jet = ColorMap.JET
gray = ColorMap.GRAY
def plotLogScale(x,y,vlabel,hlabel,title=None):
  sp = SimplePlot()
  pv = sp.addPoints(x,y)
  #pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  #pv.setLineColor(Color.black)
  pv.setLineWidth(2.0)
  #pv.setMarkSize(8.0)
  sp.setSize(1000+80,400)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  sp.setFontSize(18)
  if title:
    sp.paintToPng(720,3.3,pngDir+title+".png")

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
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(cs[i])
    pv.setLineColor(cs[i])
    pv.setLineWidth(2.0)
    pv.setMarkSize(8.0)
  #sp.setVLimits(-14,10)
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
  data = zerodouble(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.BIG_ENDIAN)
  ais.readDoubles(data)
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
  aos.writeDoubles(data)
  aos.close()
  return data


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
