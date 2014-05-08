import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from hv import *
from util import *

#pngDir = None
pngDir = "./png/"
seismicDir = "../../data/f3d/"
s1 = Sampling(90 ,1.0,0)
s2 = Sampling(420,1.0,0)
s3 = Sampling(335,1.0,0)
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
#k1,k2,k3 = n1/2,n2/2,n3/2
#k1,k2,k3 = 173,n2/2,n3/2
#k1,k2,k3 = 173,0,120
h1 = 45
h2 = 65
hs = [28,36,48,55,65,70]
hs = [70,65,55,48,36,28]
#k1,k2,k3 = 69,419,0; azimuth=130; elevation=40 # for 3D view of horizon 1
#k1,k2,k3 = 82,419,0; azimuth=100; elevation=40 # for 3D view of horizon 2
#k1,k2,k3 = 88,60,160; azimuth=285; elevation=11 # for 3D view of all horizons
k1,k2,k3 = 69,415,8; azimuth=130; elevation=30 # for 3D view of strips
fmin,fmax = -5.5,5.5
k1f,k2f,k3f = 65,406,114
k1f,k2f,k3f = 48,406,114
k1f,k2f,k3f = 48,406,0
gmin,gmax,gint,glab = -2.0,2.0,0.5,"Amplitude"
background = Color.WHITE

def main(args):
  #slopes()
  flatten()
  #display("gmep")


def slopes():
  f = readImage("gm")
  p2 = copy(f)
  p3 = copy(f)
  ep = copy(f)
  lsf = LocalSlopeFinder(4.0,4.0)
  lsf.findSlopes(f,p2,p3,ep);
  writeImage("gmp2",p2)
  writeImage("gmp3",p3)
  writeImage("gmep",ep)
  for g in [p2,p3,ep]:
    world = World()
    addImage2ToWorld(world,f,g)
    makeFrame(world)

def flatten():
  f = readImage("gm")
  p2 = readImage("gmp2")
  p3 = readImage("gmp3")
  fl = readImage("flm")
  ep = readImage("gmep")
  
  wp = readImage("wp")
  p2 = mul(d1/d2,p2)
  p3 = mul(d1/d3,p3)
  ep = pow(ep,14.0)
  k11 = [ 50, 67, 50, 52, 50, 41, 31]
  k12 = [ 40,370,391,373,365,201,283]
  k13 = [150,270,165,150,115,195,160]

  k21 = [ 77, 51, 41, 42, 42, 77, 58, 48, 40, 50, 35]
  k22 = [ 65,169, 93,160,243,110,359,261,236,269,221]
  k23 = [130, 90, 20, 50, 30,210, 65, 50,  5, 10, 10]

  k31 = [ 31, 51, 52, 33, 40]
  k32 = [394,300, 10,400,410]
  k33 = [106,298,211, 92, 39]

  k41 = [ 48, 48, 35, 27, 26,10,18]
  k42 = [ 56,408, 21,363,405,50,10]
  k43 = [300,290,190,140, 58,21,48]
  kk = [k11,k12,k13]
  sc1 = SetupConstraints(k11,k12,k13)
  sc2 = SetupConstraints(k21,k22,k23)
  sc3 = SetupConstraints(k31,k32,k33)
  sc4 = SetupConstraints(k41,k42,k43)
  kk1 = sc1.extend(500,500,p2,p3,ep,f,ep)
  kk2 = sc2.extend(500,500,p2,p3,ep,f,ep)
  kk3 = sc3.extend(500,500,p2,p3,ep,f,ep)
  kk4 = sc4.extend(500,500,p2,p3,ep,f,ep)
  #kk1 = sc1.rearrange()
  #kk2 = sc2.rearrange()
  #kk3 = sc3.rearrange()
  #kk4 = sc4.rearrange()
  k1 = [kk1[0],kk2[0],kk3[0],kk4[0]]
  k2 = [kk1[1],kk2[1],kk3[1],kk4[1]]
  k3 = [kk1[2],kk2[2],kk3[2],kk4[2]]
  k4 = [kk1[3],kk2[3],kk3[3],kk4[3]]
  #k1,k2,k3,k4=None,None,None,None
  fl = Flattener3C()
  fl.setIterations(0.01,100)
  fl.setSmoothings(8.0,2.0)
  fl.setWeight1(0.06)
  fl.setScale(0.005)
  fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep,k4,k1,k2,k3)
  gc = fm.flatten(f)
  writeImage("gmfcSeg",gc)
  s = fm.getShiftsS()
  writeImage("gmscSeg",s)
  u = fm.u1
  print "s min =",min(s),"max =",max(s),"avg =",sum(s)/n1/n2/n3
  print "wp min =",min(wp),"max =",max(wp)
  world = World()
  addImageToWorld(world,f)
  addImageToWorld(world,gc)
  makeFrame(world)
def rgbFromAmplitude(f,h,r,g,b):
  amp = zerofloat(n2*n3)
  si = SincInterpolator()
  si.setUniform(n1,1,0,n2,1,0,n3,1,0,f)
  i = 0
  for i3 in range(n3):
    for i2 in range(n2):
      amp[i] = si.interpolate(h[i3][i2],i2,i3)
      i = i+1
  aMin = -1.0
  aMax = 1.0
  mp = ColorMap(aMin,aMax,ColorMap.RED_WHITE_BLUE)
  ampRGB = mp.getRgbFloats(amp)
  i = 0
  for i3 in range(n3):
    for i2 in range(n2):
      r[i3][i2] = ampRGB[i  ] 
      g[i3][i2] = ampRGB[i+1] 
      b[i3][i2] = ampRGB[i+2] 
      i = i+3

def display(filename):
  f = readImage(filename)
  world = World()
  ipg = addImageToWorld(world,f,cmap=gray)
  ipg.setSlices(k1,k2,k3)
  makeFrame(world)

def displayHorizon(h):
  f = readImage("gm")
  r = zerofloat(n2,n3)
  g = zerofloat(n2,n3)
  b = zerofloat(n2,n3)
  rgbFromAmplitude(f,h,r,g,b)
  world = World()
  ipg = addImageToWorld(world,f,cmap=jet,cmin=-0.6,cmax=0.6)
  ipg.setSlices(k1,k2,k3)
  tg  = TriangleGroup(True,s3,s2,add(h,0.0),r,g,b)
  world.addChild(tg)
  makeFrame(world)

#############################################################################
# read/write files
def readImage2(name,n2,n3):
  fileName = seismicDir+name+".dat"
  image = zerofloat(n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image
  
def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2,n3 = s1.count,s2.count,s3.count
  #n1,n2,n3 = 120,480,430
  image = zerofloat(n1,n2,n3)
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

def readSlice3(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s


#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  ipg.setColorModel2(ColorMap.getJet())
  #ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg

def makeFrame(world):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  lightPosition=[-0.18,-0.4,0.8,0.0] # good for horizons 1 and 2
  lightPosition=[0.,0.,1.0,0.0] #default position
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  #view.setLightPosition(lightPosition)
  zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.1)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  #frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1000,900)
  frame.setVisible(True)
  frame.paintToFile("height.png")
  return frame
"""
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
"""
 
def display2(s,png=None):
  pp = PlotPanel(1,1,
    PlotPanel.Orientation.X1RIGHT_X2UP,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  pp.setHInterval(2.0)
  pp.setVInterval(2.0)
  pp.setHLabel("Crossline (km)")
  pp.setVLabel("Inline (km)")
  pv = pp.addPixels(s2,s3,slice23(k1,s))
  pv.setClips(fmin,fmax)
  pf = PlotFrame(pp)
  pf.setFontSizeForSlide(1.0,1.0)
  pf.setSize(926,510)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")
 
def display3(s,c=None,clabel="",cmin=0,cmax=0,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,s)
  pp.setSlices(k1,k2,k3)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Crossline (km)")
  pp.setLabel3("Inline (km)")
  pp.setClips(fmin,fmax)
  if c:
    cb = pp.addColorBar(clabel)
    #cb.setInterval(1.0)
    pp.setColorBarWidthMinimum(140)
    pp.setLineColor(Color.BLACK)
  else:
    pp.setLineColor(Color.YELLOW)
    #cb = pp.addColorBar("Amplitude")
    #cb.setInterval(5.0)
  pp.setInterval1(0.5)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1,200)
  if c:
    pv12 = PixelsView(s1,s2,slice12(k3,c))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(s1,s3,slice13(k2,c))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(s2,s3,slice23(k1,c))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setFontSizeForSlide(1.0,1.0)
  if c:
    pf.setSize(1036,814)
  else:
    pf.setSize(859,814)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")

def plot3f(g,a=None,amin=None,amax=None,amap=None,alab=None,aint=None,
           png=1):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X3RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  pp.setSlices(k1f,k2f,k3f)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Inline (km)")
  pp.setLabel3("Crossline (km)")
  pp.mosaic.setHeightElastic(0,180)
  pp.mosaic.setHeightElastic(1, 70)
  pp.setClips(gmin,gmax)
  pp.setColorModel(gray)
  if a:
    pp.setLineColor(Color.WHITE)
    cb = pp.addColorBar(alab)
    if aint:
      cb.setInterval(aint)
  else:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(0.5)
  pp.setInterval1(0.1)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  if a:
    pv12 = PixelsView(s1,s2,slice12(k3,a))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv12.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv13 = PixelsView(s1,s3,slice13(k2,a))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv23 = PixelsView(s2,s3,slice23(k1f,a))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    pv23.setInterpolation(PixelsView.Interpolation.NEAREST)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if amin!=amax:
        pv.setClips(amin,amax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  pp.setColorBarWidthMinimum(170)
  #pf.setFontSizeForSlide(1.0,0.8)
  pf.setSize(1200,700)
  pf.setVisible(True)
  if png and pngDir:
    png = pngDir+"f"+str(k1f)
    pf.paintToPng(360,7.0,png+".png")

def plotFrame(s1,s2,f,h,i3t):
  orient = PlotPanel.Orientation.X1DOWN_X2RIGHT
  panel  = PlotPanel(1,1,orient)
  pxv    = panel.addPixels(0,0,s1,s2,f)
  pxv.setColorModel(ColorMap.GRAY)
  ptv = panel.addPoints(0,0,h[0],h[1])
  ptv.setStyle("r-")
  panel.setTitle("section "+i3t)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
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
