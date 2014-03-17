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

from hw4 import *

#############################################################################
# parameters
nk = 20
nx = 100
dx = 0.01
fx = dx*0.5
s1 = Sampling(nx,1,0)
s2 = Sampling(nk,1,0)
sf = Sampling(nk,1,0)
sx = Sampling(nx,dx,fx)
sigma = 0.05
pngDir = "../../../HW4/images/"
dataDir = "../../../HW4/data/"
#############################################################################

def main(args):
  applyForAll()
def applyForAll():
  g  = zerodouble(nx,nk) 
  a = zerodouble(nx,nk) 
  sfi = SvdForInverse()
  sfi.constructG(sx,g)
  div(g,sigma,a)
  plot2D(s1,s2,a,1000,"nk","nx","A_ij",cmap=jet)#,png="A",orient=False)
  #########################data
  dt = zerodouble(nk) 
  dn = zerodouble(nk) 
  b  = zerodouble(nk) 
  m = readData(nx,2,"model")
  mt = m[1]
  sfi.setForNoise(0.0,sigma)
  sfi.computeData(g,mt,dt,dn)
  #sfi.normalizeData(sigma,dn,nd) 
  div(dn,sigma,b)
  div(dt,sigma,dt)
  plot1D(s2,dt,"NormalizedData","nk",x2=b)#,title="NormalizedData")
  #########################svd
  s = zerodouble(nk) 
  u = zerodouble(nk,nk) 
  v = zerodouble(nk,nx) 
  sfi.svdForG(a,s,u,v)
  plot1D(s2,s,"Singular values","Index")#,title="ss")
  #plot2D(s2,s1,v,600,"v_i","Index (j)","v_ij",cmap=jet,png="vs")
  #plot2D(s2,s2,u,400,"u_i","Index (j)","u_ij",cmap=jet,png="us")
  ########################construct coefficients
  rd = zerodouble(nk) 
  ra = zerodouble(nk) 
  vsl = zerodouble(3,nk)
  sfi.modelCoefficients(s,u,b,rd,ra)
  vsl[0]  = log10(abs(rd)) 
  vsl[1]  = log10(abs(s)) 
  vsl[2]  = log10(abs(ra))
  plotMultiple1D(s2,vsl)#,png="coefficients")
  ######################## dataMisfit, modelObjectFunction and Tikhonov curves
  phiD = zerodouble(nk)
  phiM = zerodouble(nk)
  sfi.computePhiD(rd,phiD)
  sfi.computePhiM(ra,phiM)
  q = sfi.findQ(phiD,nk)
  print q
  plot1D(sf,log10(phiD),"Log10(phiD)","Index")#,title="dm")
  phiD = log10(phiD)
  for i in range(nk-1):
    phiM[i+1] = log10(phiM[i+1])
  plot1D(sf,phiM,"Log10(phiM)","Index")#,title="mof")
  plotLogScale(phiM,phiD,"Log10(phiD)","Log10(phiM)")#,title="tc")
  plot1D(sf,phiD,"Log10(phiM)&Log10(phiD)","x",x2=phiM)#,title="dm&mof")
  ########################model construct
  mc = zerodouble(nx)
  mc1 = zerodouble(nx)
  mc2 = zerodouble(nx)
  mc3 = zerodouble(nx)
  sfi.modelConstruct([0,nk ],s,u,v,b,mc)
  sfi.modelConstruct([0,q+1],s,u,v,b,mc1)
  sfi.modelConstruct([q+1,14],s,u,v,b,mc2)
  sfi.modelConstruct([14,nk],s,u,v,b,mc3)
  plot1D(sx,mt,"Model1","x",x2=mc,title="model")
  plot1D(sx,mt,"Model","x",x2=mc1,title="model1")
  plot1D(sx,mt,"Model","x",x2=mc2,title="model2")
  plot1D(sx,mt,"Model","x",x2=mc3,title="model3")
  #plot2D(v,cmap=jet)
##################################################################
# plots
jet = ColorMap.JET
gray = ColorMap.GRAY
def plotLogScale(x,y,vlabel,hlabel,title=None):
  sp = SimplePlot()
  pv = sp.addPoints(x,y)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  #pv.setLineColor(Color.magenta)
  pv.setLineColor(Color.black)
  pv.setLineWidth(2.0)
  pv.setMarkSize(8.0)
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
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    #pv.setLineColor(Color.magenta)
    pv.setLineColor(Color.black)
    pv.setLineWidth(2.0)
    pv.setMarkSize(8.0)
  sp.setSize(1000+80,400)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  sp.setFontSize(18)
  if title:
    sp.paintToPng(720,3.3,pngDir+title+".png")

def plotMultiple1D(s,v,png=None):
  sp = SimplePlot()
  colors = [Color.red,Color.green,Color.blue]
  '''
  colors = [Color.red,Color.magenta,Color.blue,
            Color.cyan,Color.green,Color.yellow]
  '''
  for i in range(3):
    pv = sp.addPoints(s,v[i])
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(colors[i])
    pv.setLineColor(colors[i])
    #pv.setLineWidth(2.0)
    pv.setMarkSize(8.0)
  sp.setVLimits(-14,10)
  sp.setVInterval(4)
  sp.setVLabel("Log10(v)")
  sp.setHLabel("Index")
  sp.setSize(1000+80,400)
  if png:
    sp.paintToPng(720,3.3,pngDir+png+".png")

def plot2D(s1,s2,x,wd,vlabel,hlabel,cbar,cmap=jet,png=None,orient=True):
  n1 = s1.getCount()
  n2 = s2.getCount()
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  if orient:
    sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
    pv = sp.addPixels(s1,s2,x)
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  else:
    sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
    pv = sp.addPixels(s1,s2,x)
    pv.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
  pv.setColorModel(cmap);
  dims = [n1,n2]
  sp.setSize(wd+80,400);
  sp.setVInterval(2)
  if (wd>400):
    sp.setHInterval(10)
  else:
    sp.setHInterval(2)
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

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
