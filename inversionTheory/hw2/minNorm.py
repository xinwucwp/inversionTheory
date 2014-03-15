import sys
from math import *
from java.awt import *
from java.lang import *
from java.util import *
from javax.swing import *
from java.util.Random import *

from edu.mines.jtk.la import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from hw2 import *

#############################################################################
# parameters
n = 20
dm = 0.05
fm = 0.00
nm = int(1.0/dm)+1
s = Sampling(n,1.0,1.0)
s1 = Sampling(n,1.0,1.0)
s2 = Sampling(n,1.0,1.0)
sm = Sampling(nm,dm,fm)
pngDir = "../../HW2/images/"
#############################################################################

def main(args):
  applyForAll()
def applyForAll():
  mnm = MinNormModel()
  e = zerodouble(n)
  d = zerodouble(n)
  v = zerodouble(n,n)
  g = zerodouble(n,n)
  alpha = zerodouble(n)
  mnm.applyForData(d) 
  mnm.kernelMatrix(g)
  mnm.eigenDecompo(g,e,v)
  mnm.coefficients(d,e,v,alpha)
  mt = zerodouble(nm)
  mc = zerodouble(nm)
  mnm.trueAndConstructed(sm,alpha,mt,mc)
  print e[0]/e[n-1]
  print e[n-1]
  #plot1D(sm,mt,"models","x",x2=mc,title="model_6")
  #plot1D(sm,sub(mc,mt),"difference","x",title="diff_06")
  #plot1D(s,e,"eigenvalues","index",title="eigenvalues")
  #plotEigenvectors(s,v)
  #plot1D(s,d,"data","index",title="data")
def plotMatrix(n):
  e = zerodouble(n)
  v = zerodouble(n,n)
  g = zerodouble(n,n)
  mnm = MinNormModel()
  mnm.kernelMatrix(g)
  #evd = DMatrixEvd(DMatrix(g))
  #dme = evd.getD()
  #dmv = evd.getV()
  mnm.eigenDecompo(g,e,v)
  plot2D(s1,s2,g)
  plot2D(s1,s2,e)
  plot2D(s1,s2,v)
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
   pv1.setLineColor(Color.blue)
   pv2.setLineColor(Color.red)
   pv1.setLineWidth(2.0)
   pv2.setLineWidth(2.0)
  else:
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    #pv.setLineColor(Color.magenta)
    pv.setLineColor(Color.black)
    pv.setLineWidth(2.0)
    pv.setMarkSize(8.0)
  sp.setSize(1000,400)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  if title:
    sp.paintToPng(720,3.3,pngDir+title+".png")
def plotEigenvectors(s,v):
  sp = SimplePlot()
  colors = [Color.red,Color.magenta,Color.blue,
            Color.cyan,Color.green,Color.yellow]
  for i in range(6):
    pv = sp.addPoints(s,v[i])
    pv.setLineColor(colors[i])
    pv.setLineWidth(2.0)
    pv.setMarkSize(8.0)
  sp.setVLabel("eigenvectors")
  sp.setHLabel("index")
  sp.setSize(1000,500)
  sp.paintToPng(720,3.3,pngDir+"eigenvectors.png")
#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
