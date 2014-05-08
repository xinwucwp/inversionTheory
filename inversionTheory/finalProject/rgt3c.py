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
pngDir = "../../png/hv/rgt/"
seismicDir = "../../data/f3d/"
n1,n2,n3=90,420,335
f1,f2,f3=1.480,0.000,5.250
d1,d2,d3=0.004,0.025,0.025
s1 = Sampling(90 ,d1,f1)
s2 = Sampling(420,d2,f2)
s3 = Sampling(335,d3,f3)
h1,h2 =48,68#for slide 
h1,h2 =46,68#for print
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
k11=[ 50, 67, 50, 52, 50, 41, 35]
k12=[ 40,370,391,373,365,201,280]
k13=[150,270,165,150,115,195,160]
k21=[ 77, 51,41, 42, 42, 77, 58, 48, 40, 50, 35, 57]
k22=[ 65,169,93,160,243,110,359,261,236,269,221,270]
k23=[130, 90,20, 50, 30,210, 65, 50,  5, 10, 10,105]
k31 = [ 31, 51, 52, 33, 40]
k32 = [394,300, 10,400,410]
k33 = [106,298,211, 92, 39]
k41 = [ 48, 48, 35, 27, 26,10,18]
k42 = [ 56,408, 21,363,405,50,10]
k43 = [300,290,190,140, 58,21,48]
cp1  = [k11,k12,k13] 
cp2  = [k21,k22,k23] 

def main(args):
  #slopes()
  #flatten()
  slices()
  pg1 = setPointGroup(cp1,8.0)  
  pg2 = setPointGroup(cp2,8.0)  
  k1,k2,k3 = 69,419,5
  azimuth,elevation=130,40 
  horizonExtraction(k1,k2,k3,azimuth,elevation,[48,46],1.0,pg1,"surf1")
  k1,k2,k3 = 85,419,5
  azimuth,elevation=100,40
  horizonExtraction(k1,k2,k3,azimuth,elevation,[h2,h2],1.5,pg2,"surf2")

def slopes():
  f = readImage("gm")
  u1 = copy(f)
  u2 = copy(f)
  u3 = copy(f)
  p2 = copy(f)
  p3 = copy(f)
  ep = copy(f)
  lof = LocalOrientFilter(2.0,1.0)
  lof.applyForNormalPlanar(f,u1,u2,u3,ep)
  lsf = LocalSlopeFinder(2.0,1.0)
  lsf.findSlopes(f,p2,p3,ep);
  writeImage("gmu1",u1)
  writeImage("gmu2",u2)
  writeImage("gmu3",u3)
  writeImage("gmp2",p2)
  writeImage("gmp3",p3)
  #writeImage("gmep",ep)
  for g in [u1,u2,u3,ep]:
    world = World()
    addImage2ToWorld(world,f,g)
    makeFrame(s1,s2,s3,azimuth,elevation,world)
def flatten():
  d1,d2,d3=1.0,1.0,1.0
  f1,f2,f3=0.0,0.0,0.0
  s1 = Sampling(90 ,d1,f1)
  s2 = Sampling(420,d2,f2)
  s3 = Sampling(335,d3,f3)
  f = readImage("gm")
  u1 = readImage("gmu1")
  u2 = readImage("gmu2")
  u3 = readImage("gmu3")
  p2 = readImage("gmp2")
  p3 = readImage("gmp3")
  fl = readImage("flm")
  ep = readImage("gmep")
  p2 = mul(d1/d2,p2)
  p3 = mul(d1/d3,p3)
  wh = pow(ep,10.0)
  ws = pow(ep,6.0)
  wp = pow(ep,6.0)
  thin = Thin();
  ut1 = zerofloat(n1,n2,n3)
  ut2 = zerofloat(n1,n2,n3)
  thin.applyForHorizontal(fl,ut1,ut2)
  sub(1.0,ut2,ut2)
  ut2=pow(ut2,2.0)
  mul(wp,ut2,wp)
  fl3 = Flattener3()
  fl3.setIterations(0.01,100)
  fl3.setSmoothings(2.0,4.0)
  fl3.setWeight1(0.08)
  flm = fl3.getMappingsFromSlopes(s1,s2,s3,p2,p3,wp)
  g = flm.flatten(f)
  writeImage("gmfSeg",g)
  s = flm.getShiftsS()
  writeImage("gmsSeg",s)
  u = flm.u1
  writeImage("rgtSeg",u)
  print "s min =",min(s),"max =",max(s),"avg =",sum(s)/n1/n2/n3
  print "wp min =",min(wp),"max =",max(wp)
  world = World()
  addImageToWorld(s1,s2,s3,world,f)
  addImageToWorld(s1,s2,s3,world,g)
  makeFrame(s1,s2,s3,azimuth,elevation,world)

def flattenC():
  d1,d2,d3=1.0,1.0,1.0
  f1,f2,f3=0.0,0.0,0.0
  s1 = Sampling(90 ,d1,f1)
  s2 = Sampling(420,d2,f2)
  s3 = Sampling(335,d3,f3)
  f = readImage("gm")
  u1 = readImage("gmu1")
  u2 = readImage("gmu2")
  u3 = readImage("gmu3")
  p2 = readImage("gmp2")
  p3 = readImage("gmp3")
  fl = readImage("flm")
  ep = readImage("gmep")
  p2 = mul(d1/d2,p2)
  p3 = mul(d1/d3,p3)
  wh = pow(ep,10.0)
  ws = pow(ep,6.0)
  wp = pow(ep,6.0)
  thin = Thin();
  ut1 = zerofloat(n1,n2,n3)
  ut2 = zerofloat(n1,n2,n3)
  thin.applyForHorizontal(fl,ut1,ut2)
  sub(1.0,ut2,ut2)
  ut2=pow(ut2,2.0)
  mul(wp,ut2,wp)
  sigma1,sigma2=6.0,6.0
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
  sc1.setForExtension(sigma1,sigma2,0.0)
  sc2.setForExtension(sigma1,sigma2,0.0)
  sc3.setForExtension(sigma1,sigma2,0.0)
  sc4.setForExtension(sigma1,sigma2,0.0)
  kk1 = sc1.extend(500,500,p2,p3,wh,f,wh)
  kk2 = sc2.extend(500,500,p2,p3,wh,f,wh)
  kk3 = sc3.extend(500,500,p2,p3,wh,f,wh)
  kk4 = sc4.extend(500,500,p2,p3,wh,f,wh)
  '''
  np1 = 140485
  np2 = 99547
  np3 = 140700
  np4 = 140700
  kk1 = readImage2("c1",np1,4)
  kk2 = readImage2("c2",np2,4)
  kk3 = readImage2("c3",np3,4)
  kk4 = readImage2("c4",np4,4)
  '''
  k51 = [ 37, 41, 46, 49, 59, 67, 67, 71]
  k52 = [185,205,223,239,256,276,298,312]
  k53 = [  0,  0,  0,  0,  0,  0,  0,  0]
  sc5 = SetupConstraints(k51,k52,k53)
  sc5.setForExtension(sigma1,sigma2,0.0)
  kk5 = sc5.extend(15,25,p2,p3,wh,f,wh)
  k1 = [kk1[0],kk2[0],kk3[0],kk4[0],kk5[0]]
  k2 = [kk1[1],kk2[1],kk3[1],kk4[1],kk5[1]]
  k3 = [kk1[2],kk2[2],kk3[2],kk4[2],kk5[2]]
  k4 = [kk1[3],kk2[3],kk3[3],kk4[3],kk5[3]]
  w1 = zerofloat(n1,n2,n3)
  rgt = RelativeGeologicTime3C()
  rgt.setIterations(0.01,100)
  rgt.setSmoothings(2.0,4.0)
  rgt.setWeight1(0.08)
  rgt.setScale(0.00001)
  rgtm = rgt.getMappingsFromSlopes(s1,s2,s3,u1,u2,u3,wp,ws,k4,k1,k2,k3)
  gc = rgtm.flatten(f)
  writeImage("gmfcSeg",gc)
  s = rgtm.getShiftsS()
  writeImage("gmscSeg",s)
  u = rgtm.u1
  writeImage("gmsuSeg",u)
  print "s min =",min(s),"max =",max(s),"avg =",sum(s)/n1/n2/n3
  print "wp min =",min(wp),"max =",max(wp)
  world = World()
  addImageToWorld(s1,s2,s3,world,f)
  addImageToWorld(s1,s2,s3,world,gc)
  makeFrame(s1,s2,s3,azimuth,elevation,world)

def slices():
  fname,gname,rgtname="gm","gmfcSeg","gmsuSeg"
  #fname,gname,rgtname="gm","gmfSeg","rgtSeg"
  fc = readImage(gname)
  fu = readImage(rgtname)
  gm = readImage(fname)
  mul(fu,d1,fu)
  add(fu,1.48,fu)
  fc = gain(fc)
  gm = gain(gm)
  fmin,fmax= -3.5,3.0
  rmin,rmax=1.42,2.0
  k2,k3=406,5
  k1= h2
  label1,cbar = "Time (s)","Amplitude"
  plot3f(k1,k2,k3,gm,cbar,label1,gmap=gray,cint=1.0,gmin=fmin,gmax=fmax,png=fname)
  label1,cbar = "RGT","Amplitude"
  k1 = h1
  plot3f(k1,k2,k3,fc,cbar,label1,gmap=gray,cint=1.0,gmin=fmin,gmax=fmax,png=gname+"fl1")
  k1 = h2
  plot3f(k1,k2,k3,fc,cbar,label1,gmap=gray,cint=1.0,gmin=fmin,gmax=fmax,png=gname+"fl2")
  label1,cbar = "Time (s)","RGT"
  k1 = h2
  plot3f(k1,k2,k3,fu,cbar,label1,gmap=jet,cint=0.2,gmin=rmin,gmax=rmax,png=rgtname)
def setPointGroup(kk1,size):
  np  = len(kk1[0])
  xyz = zerofloat(np*3)
  rgb = zerofloat(np*3)
  ki = 0
  for i in range(np):
    xyz[ki  ] = kk1[2][i]
    xyz[ki+1] = kk1[1][i]
    xyz[ki+2] = kk1[0][i]+1.0
    rgb[ki  ]  = 0.2
    rgb[ki+1]  = 0.8
    rgb[ki+2]  = 0.2 
    ki = ki+3
  pg = PointGroup(size,xyz,rgb);
  states = StateSet();
  cs = ColorState();
  cs.setColor(Color.GREEN);
  states.add(cs);
  lms = LightModelState();
  lms.setTwoSide(True);
  states.add(lms);
  ms = MaterialState();
  ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
  ms.setShininess(100.0);
  states.add(ms);
  pg.setStates(states);
  return pg;
def rgbFromHeight(h,r,g,b):
  n1 = len(h[0])
  n2 = len(h)
  ht = zerofloat(n1*n2)
  mp = ColorMap(-max(h),-min(h),ColorMap.JET)
  i = 0
  for i1 in range(n1):
    for i2 in range(n2):
      ht[i] = -h[i2][i1]
      i=i+1
  htRGB = mp.getRgbFloats(ht)
  i = 0
  for i1 in range(n1):
    for i2 in range(n2):
      r[i2][i1] = htRGB[i  ] 
      g[i2][i1] = htRGB[i+1] 
      b[i2][i1] = htRGB[i+2] 
      i = i+3
def horizonExtraction(k1,k2,k3,azimuth,elevation,tau,shift,pg,png):
  f = readImage("gm")
  su = readImage("gmsSeg")
  sc = readImage("gmscSeg")
  f = gain(f)
  i = 0
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  s3 = Sampling(n3,1.0,0.0)
  for s in [su,sc]:
    r = rampfloat(0.0,1.0,0.0,0.0,n1,n2,n3)
    t = add(s1.first,mul(s1.delta,sub(r,s)))
    h  = slice23(tau[i],t)
    h  = add(h,shift) 
    for i2 in range(n2):
      for i3 in range(n3):
        if h[i3][i2]>=n1:
          h[i3][i2] = (n1-1)
    hr = zerofloat(n2,n3)
    hg = zerofloat(n2,n3)
    hb = zerofloat(n2,n3)
    rgbFromHeight(h,hr,hg,hb)
    world = World()
    ipg = addImageToWorld(s1,s2,s3,world,f,cmap=gray,cmin=-3.5,cmax=3)
    ipg.setSlices(k1,k2,k3)
    tg = TriangleGroup(True,s3,s2,h,hr,hg,hb)
    if(i==1):
      world.addChild(pg)
      png = png+"c"
    world.addChild(tg)
    makeFrame(s1,s2,s3,azimuth,elevation,world,png=png)
    i = i+1
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

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET

def addImageToWorld(s1,s2,s3,world,image,cmap=gray,cmin=0,cmax=0):
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

def makeFrame(s1,s2,s3,azimuth,elevation,world,png=None):
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
  if png and pngDir:
    png = pngDir+png
    frame.paintToFile(png+".png")
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
def plot3f(k1,k2,k3,g,cbar,label1,gmap=gray,cint=0.5,gmin=None,gmax=None,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  pp.setSlices(k1,k2,k3)
  pp.setInterpolation(PixelsView.Interpolation.LINEAR)
  pp.setLabel1(label1)
  pp.setLabel2("Inline (km)")
  pp.setLabel3("Crossline (km)")
  pp.mosaic.setHeightElastic(0,180)
  pp.mosaic.setHeightElastic(0, 220)
  if gmin !=gmax:
    pp.setClips(gmin,gmax)
  pp.setColorModel(gmap)
  pp.setLineColor(Color.YELLOW)
  cb = pp.addColorBar(cbar)
  cb.setInterval(cint)
  #pp.setInterval1(0.1)# for print
  pp.setInterval1(0.2)# for slide
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  #pp.setColorBarWidthMinimum(45)
  #pf.setFontSizeForPrint(7.0,480)
  pp.setColorBarWidthMinimum(110)
  pf.setFontSizeForSlide(1.0,1.0,16.0/9.0)
  pf.setSize(1200,750)
  pf.setVisible(True)
  if png and pngDir:
    png = pngDir+png
    pf.paintToPng(360,6.66,png+".png")
 
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
