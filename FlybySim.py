#!/usr/bin/env python3
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# (c) Maisie Rashman, Department of Physics, University of Oxford

import argparse
from astropy.io import ascii
from astropy.table import QTable, MaskedColumn
from astropy import units as u
import imageio
import math 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.patheffects as pe
from matplotlib.offsetbox import AnchoredText, OffsetImage, AnnotationBbox
import numpy as np
import os
import pandas as pd
import pathlib
from PyQt5.QtWidgets import QWidget, QHBoxLayout
import seaborn as sns
from shapely import LineString
from shapely.geometry.polygon import Polygon
import sys
from urllib.request import Request, urlopen
import visvis as vv
from visvis import Pointset
import time

__version__ = '4.0'
__date__ = '2024-03-11'
__author__ = 'Maisie Rashman'
__all__ = ['arc','downloadFile', 'parseArguments', 'progress', 'rotateArc', \
'rotateFilters', 'MainWindow', 'view', 'cometMesh', 'coma', 'window', 'FPA', 'sim']

"""
#########################################################################################
										GLOBAL DEFINITIONS
#########################################################################################
"""
verts = [
(50., 125.), 
(100., 87.5),  
(100., 125.),  
(150., 87.5), 
(100., 50.),
(100., 87.5),  
(50., 50.),
(50., 125.)]
codes = [
Path.MOVETO,
Path.LINETO,
Path.LINETO,
Path.LINETO,
Path.LINETO,
Path.LINETO,
Path.LINETO,
Path.CLOSEPOLY]

global path
path = Path(verts, codes)

global wdth
wdth = os.get_terminal_size()[0]
"""
#########################################################################################
										FUNCTIONS
#########################################################################################
"""
def arc(ox, oy):
	"""
	Returns the pixel coordinates of a ~9 deg circular arc centred on the TIRI FOV.		
	
	####  INPUT	 ####
	ox		: float
		x-coordinate about which the rotation should occur. Almost always set to the 
		centre of the pixel array.
				  
	oy		: float
		y-coordinate about which the rotation should occur. Almost always set to the 
		centre of the pixel array.
			
	####  OUTPUT  ####	   
	x		: float
		x-coordinates of arc. 
		   
	y		: float
		y-coordinates of arc.
		
	"""
	r = ox/np.sin((ox*0.26e-3)/2)
	theta = np.linspace(-1*np.radians(9), np.radians(9), 1280)
	y = r*np.cos(theta) - (r-oy)
	x = r*np.sin(theta) + ox 
	return x,y
	
def downloadFile(dl_path, dl_url):
	"""
	Downloads files from url and saves to directory.
	
	####  PARAMETERS ####
	DL_PATH : str
			  Download path on the local machine, relative to this function.
	
	DL_URL	: str
			  Download url of the requested file.		
	"""
	file_name = dl_url.split('/')[-1]
	pathlib.Path(dl_path).mkdir(parents=True, exist_ok=True)
	if not os.path.isfile(dl_path + file_name):
		opener=urllib.request.build_opener()
		opener.addheaders=[('User-Agent','Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/36.0.1941.0 Safari/537.36')]
		urllib.request.install_opener(opener)
		urllib.request.urlretrieve(dl_url, dl_path + file_name)
	return
	
	
def parseArguments():
	"""
	Defines parser for command line interface
	
	####  RETURNS  ####	
	parser	: argparse.ArgumentParser
		instance of argparse.ArgumentParser
	"""
	# Create argument parser
	parser = argparse.ArgumentParser(
					prog='FlybySim',
					description="""FLYBYSIM runs a simulation of the comet interceptor \
					encounter for the TIRI sub-module of the instrument MIRMIS.""",
					epilog="Comet download and display written using Thomas Albin's \
					Space Science Tutorials - Part 12: https://github.com/ThomasAlbin/\
					SpaceScienceTutorial\n\nComet meshes have to be found and downloaded \
					prior to using.\nMost available comet shape models can be found \
					here: https://sbn.psi.edu/pds/shape-models\n\nExample download of \
					c67P mesh:\ndirectory_path = 'path/to/local'\ndownload_url = \
					'https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos\
					/rosetta/raw/misc/models/ROS_CG_M001_OSPCLPS_N_V1.OBJ'\n\
					download_file(directory_path, download_url)", formatter_class = \
					argparse.ArgumentDefaultsHelpFormatter)

	# Positional mandatory arguments
	parser.add_argument("mesh", type=str, metavar="path to mesh", help = \
	"path to comet shape model .OBJ file stored on local disk.")


	# Optional arguments
	parser.add_argument('-rest', type=str, default='tiri', dest ='rFrame', \
	choices=['tiri', 'sc'], help="Rest frame of simulation view")

	parser.add_argument("-ff", type=int, default = 20, dest="upSpeed", \
	help='Fast forward non-datacube simulation frames by this amount.')
	parser.add_argument("-V", type=float, default=70, help='Encounter velocity \
	vector in km/s.')
	parser.add_argument("-Vc", type=float, default=15, help='Comet velocity in \
	km/s.')
	parser.add_argument("-ca", type=float, default=1000, help="Closest approach \
	distance in km.")
	parser.add_argument("-o", type=str, dest='outLink', default=os.getcwd(), \
	help="Path to local directory to save output. If not provided, defaults to current \
	working directory.")
	parser.add_argument("-readfreq", type=float, dest='tRes', default=1/20, \
	help="Readout frequency of the TIRI detector in Hz.")
	parser.add_argument('-fps', type=float, dest = 'fRate', default=20, \
	help="Number of frames per second in output gif.")
	parser.add_argument('--nosave', action='store_false', dest='saveOutput', help= \
	"Don't save model gif and table output.")	
	parser.add_argument("--noann", action='store_false', dest='annotate', \
	help="Don't annotate simulation frames with TIRI filters, etc.")
	parser.add_argument("--saveA", action='store_true', dest='saveAnn', \
	help="Don't save annotated simulation frames as .npy arrays.")		
	parser.add_argument('--saveraw', action='store_true', dest='saveRaw', help \
	= 'Save raw simulation frames to local disk.')
	return parser

def progress(txt=None,p=False):
	if txt != None:
		print(txt,end="", flush=True)
		time.sleep(0.2) 
	else:
		if p == 1:
			print('.\n', end="\n")
			time.sleep(0.2)
		else:
			print('.',end="", flush=True)
			time.sleep(0.2)
	return

def rotateArc(angle, ox, oy):
	"""
	Returns the pixel coordinates of a ~9 deg circular arc rotated about the centre of 
	TIRI FOV by a given angle.
	
	####  INPUT	 ####
	angle	: float
		Angle of rotation in deg
				  
	ox		: float
		x-coordinate about which the rotation should occur. Almost always set to the 
		centre of the pixel array.
				  
	oy		: float
		y-coordinate about which the rotation should occur. Almost always set to the 
		centre of the pixel array.	  
				
	####  OUTPUT  ####	   
	xr		: float
		x-coordinates of rotated arc.
		 
	yr		: float
		y-coordinates of rotated arc.
			
	"""
	r = ox/np.sin((ox*0.26e-3)/2)
	theta = np.linspace(-1*np.radians(12), np.radians(12), 1000) #set to 12 to account for rotation
	y = r*np.cos(theta) - (r-oy)
	x = r*np.sin(theta) + ox
	cosang, sinang = np.cos(np.radians(angle)), np.sin(np.radians(angle))	
	xr = [ (x[i]-ox)*cosang-(y[i]-oy)*sinang+ox for i in range(0,len(x))]
	yr = [ (x[i]-ox)*sinang-(y[i]-oy)*cosang+oy for i in range(0,len(x))]	
	return xr,yr

def rotateFilters(border, angle, ox, oy):
	"""
	Takes the boundary edge of a filter and rotates it in the FOV by a given angle.

	####  INPUT	 ####
	border		: array-like, float
		x- and y- coordinates of the straight edge of the filter in the format 
		[[x0, x1, ..., xn],[y0, y1, ..., yn]]
		
	angle		: float
		Angle of rotation in degrees.

	ox			: float
		x-coordinate about which the rotation should occur. Almost always set 
		to the centre of the pixel array.
		  
	oy			: float
		y-coordinate about which the rotation should occur. Almost always set 
		to the centre of the pixel array. 
			
	####  OUTPUT  ####	
	line		: array-like, float
		x- and y- coordinates at the extremes of the newly rotated filter edge
		in the format [[x0_r, xn_r], [y0_r, yn_r]]. 
					
	"""	 
	
	cosang, sinang = np.cos(np.radians(angle)), np.sin(np.radians(angle)) 
	p1x, p1y = border[0][0], border[1][0]
	p2x, p2y = border[0][-1], border[1][-1]
	q1x, q1y = ox + cosang * (p1x - ox) - sinang * (p1y - oy), oy + sinang * (p1x - ox) + cosang * (p1y - oy)
	q2x, q2y = ox + cosang * (p2x - ox) - sinang * (p2y - oy), oy + sinang * (p2x - ox)+ cosang * (p2y - oy)
	line = [[q1x, q2x], [q1y, q2y]]
	return line

"""
#########################################################################################
										CLASSES
#########################################################################################
"""
class coma: 
	"""
	Creates a very inauthentic representation of a coma. Just for visuals, not 
	scientifically relevant.
	"""
	
	__slots__ = ("cloud")
	
	def __init__(self): 
		self.cloud = self.generateComa()
		
	def generateComa(self):
		"""
		Generates a 29000 point cloud to describe the rough shape of an idealised coma. 
		The cloud is composed of four different parts (for no purpose but
		~aesthetics~):
			1)	a 3d cone of 6000 points with coordinates in the range 
				(-38 < x < 35)km, (-38 < y < 880)km, (-41 < z < 37)km.
			2)	a 3D normal distribution of 10000 points with coordinates in 
				the range (-1900 < x < 2150)km, (-3100 < y < 5200)km, 
				(-1900 < z < 2000)km.
			3)	a 3D cone of 3000 points with coordinates in the range 
				(-8 < x < 8)km, (-5 < y < 260)km, (8 < z < 8)km.
			4)	a 3D normal distribution of 10000 points with coordinates in 
				the range (-200 < x < 200)km, (-3000 < y < 5300)km, 
				(-200 < z < 200)km.

		#### OUTPUT ####
		pp		=  Point set with 29000 points		
		"""
		
		a = 1000
		b = 15
		h = a * (np.random.random(3000)) ** (1/3)
		r = ((b/a) * h * np.sqrt(np.random.random(3000)))

		t =	 2*np.pi * np.random.random(3000)
		x = r * np.cos(t)*np.random.random(3000)*2
		y = (a-h)
		z = r * np.sin(t)*np.random.random(3000)*2

		a = np.concatenate([(np.asarray([x,y,z])).T, np.random.normal(size=(3000,3))*10])
		a[:,1] = a[:,1]
		a2 = np.squeeze(np.dstack([np.random.normal(size=(10000))*500, np.random.normal(size=(10000))*1000, np.random.normal(size=(10000))*500]))
		a2[:,1] = a2[:,1]+980

		a4 = np.squeeze(np.dstack([np.random.normal(size=(10000))*50, np.random.normal(size=(10000))*100, np.random.normal(size=(10000))*50]))
		a4[:,1] = a2[:,1]+100
	
		a3 = (np.asarray([x,y,z])*0.3).T
		a3[:,1] = a3[:,1]-5

		pp = Pointset(np.concatenate([a, a2, a3, a4]))
		return pp  
		
class cometMesh:
	
	"""
	Loads and displays 3D shape model as a visvis mesh object which can be translated 
	and rotated interactively in the window. 
	"""
	
	def __init__(self, faces, vertices):
		"""		
		#####  INPUT ####
		faces		=  Array containing the location of faces in the shape model.
	
		vertices	=  Array containing the location of vertices in the shape model.
		""" 
		self.faces = faces
		self.vertices = vertices
		self.obj = vv.mesh(vertices=self.vertices, faces=self.faces, verticesPerFace=3)
		self.obj.specular = 0.2
		self.obj.diffuse = 0.9
		self.Rot = vv.Transform_Rotate(angle=100, ax=-1, ay=1, az=-1)
		self.obj.transformations.insert(0,self.Rot)		  

class FPA:
	""" 
	Routines to simulate the projection of the TIRI filter array on the detector. 
	"""
	
	def __init__(self):
		self.loadInfo()
		self.factor = 1
		
	def loadInfo(self):
		"""
		Loads and stores info on the layout of the TIRI filter array 
		"""
		self.nbN = 8
		self.bbN = 1
		self.nb_xy = [32,480]
		self.bb_xy = [200,140]
		self.gap_nb = 8
		self.gap_bb = 57
		self.m_h = 22
		self.m_v = 20
		self.order = ['broadband', 'A', 'B', 'C', 'D', 'broadband', 'E', 'F', 'G', 'H']
		self.integrate = {'A': 0.0552, 'B': 0.0448, 'C': 0.0400, 'D': 0.0435, 'E': 0.0520, 'F': 0.0743, 'G': 0.1452, 'H': 0.0710, 'broadband': 0.0002}
		return
		
	def setScale(self, n):
		"""
		Can be used to adjust the scale of the filter layout to account for different
		dpi values in the simulation window. Base simulation has twice as many pixels 
		and therefore needs the filter array scaling by a factor of 2.
		"""
		self.factor = n
		return

	def getFPAdims(self):
		"""
		Returns the x,y pixel dimensions of the simulated TIRI focal plane array. 
		"""
		xy_dim = [640,480]
		xy_dim = [xy_dim[0]*self.factor, xy_dim[1]*self.factor]
		return xy_dim
	
	def verticalCentres(self):
		"""
		Returns the pixel location of the centre of each filter in the x-axis
		"""
		n = self.factor
		centres = []
		ind = np.arange(0,4,1)
		for i in ind:
			centres.append(((self.m_h*0.5)+(0.5*self.nb_xy[0])+(i*self.nb_xy[0])+(i*self.gap_nb))*n)
		end_blc = centres[-1]+(0.5*self.nb_xy[0]*n)
		start_blc = end_blc+(self.gap_bb*2*n)+(self.bb_xy[0]*n)
		for i in ind:
			centres.append((start_blc+(0.5*self.nb_xy[0]*n)+(i*self.nb_xy[0]*n)+(i*self.gap_nb*n)))
		return centres
		
	def horizontalCentres(self):
		"""
		Returns the pixel location of the centre of each filter in the y-axis		
		"""
		n = self.factor
		centres = []
		ind = np.arange(0,4,1)
		for i in ind:
			centres.append(((self.m_v*0.5)+(0.5*self.nb_xy[0])+(i*self.nb_xy[0])+(i*self.gap_nb))*n)
		end_blc = centres[-1]+(0.5*self.nb_xy[0]*n)
		start_blc = end_blc+(self.gap_nb*2*n)+(self.bb_xy[1]*n)
		for i in ind:
			centres.append((start_blc+(0.5*self.nb_xy[0]*n)+(i*self.nb_xy[0]*n)+(i*self.gap_nb*n)))
		return centres
	
	def filterOverlay(self):
		"""
		Maps out the location of all the filters as they appear on the TIRI fpa 
		and returns (x,y) line arrays of the filter edges.
		"""
		n = self.factor
		if not hasattr(self, "xy_dim"):
			self.xy_dim = self.getFPAdims() 
		fEdges = {}
		vc = self.verticalCentres()
		hc = self.horizontalCentres()
		for i, name in enumerate(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'broadband']):
			xhL = np.arange(self.xy_dim[0]/2 - self.bb_xy[0]*n / 2, self.xy_dim[0]/2+ self.bb_xy[0]*n / 2+1, 1) 
			if name == 'broadband':
				yhL = np.asarray([self.xy_dim[1]/2-self.bb_xy[1]*n / 2]*len(xhL))
				yhR = np.asarray([self.xy_dim[1]/2+ self.bb_xy[1]*n / 2]*len(xhL))
				yvL = np.arange(self.xy_dim[1]/2- self.bb_xy[1]*n / 2, self.xy_dim[1]/2+ self.bb_xy[1]*n / 2, 1)
				xvL = np.asarray([self.xy_dim[0]/2- self.bb_xy[0]*n / 2]*len(yvL))
				xvR = np.asarray([self.xy_dim[0]/2+ self.bb_xy[0]*n / 2]*len(yvL))
			else:
				yhL = np.asarray([hc[i]- self.nb_xy[0]*n / 2]*len(xhL))
				yhR = np.asarray([hc[i]+ self.nb_xy[0]*n / 2]*len(xhL)) 
				yvL = np.arange(0,self.xy_dim[1], 1)
				xvL = np.asarray([vc[i]- self.nb_xy[0]*n / 2]*len(yvL))
				xvR = np.asarray([vc[i]+ self.nb_xy[0]*n / 2]*len(yvL))						
			fEdges.update({name: {'vL' : [xvL,yvL], 'vR': [xvR,yvL], 'hL': [xhL,yhL], 'hR': [xhL,yhR]}})
		return fEdges

class MainWindow(QWidget):
	"""
	Creates application object for the pyqt5 backend and embeds a visvis figure in the
	underlying GUI app.
	"""
	def __init__(self, *args):
		QWidget.__init__(self, *args)
		self.fig = vv.backends.backend_pyqt5.Figure(self)
		self.sizer = QHBoxLayout(self)	   
		self.sizer.setContentsMargins(0,0,0,0)
		self.sizer.addWidget(self.fig._widget, 2)
		self.setWindowTitle('Comet Interceptor')
		self.show()

class window:
	"""
	Parent class wrapper for all things visual.
	"""
	
	def __init__(self, V=70, Vc=15, Ca = 1000):
		"""		
		##### INPUT ##### 
		Vc		 =	Float or integer of the desired velocity the simulated 
								comet will travel at. Informs the view of the comet for 
								each timestep.
								Default is 15km/s
							
		V = Float or integer of the desired flyby relative velocity, 
								defined in ESA-COMET-SYS-RS-001 as the modulus of the 
								spacecraft heliocentric velocity - comet heliocentric 
								velocity. For simplicity, the spacecraft velocity will be 
								taken as interceptorVelocity minus cometVelocity for all 
								calculations. Informs the view of the comet for each 
								timestep.
								Default is 70km/s
								Can be adjusted with X.interceptorVelocity = float
		"""
		
		self.cometVelocity = Vc
		self.interceptorVelocity = V
		self.closeApproach = Ca
		
	def whichComet(self, cometPath):
		"""
		Assigns path to comet shape model on local disk
		"""
		
		progress(f'Preparing flyby window')
		self.cometLink = cometPath
		progress()
		progress(p=True)
		return
		
	def setup(self):
		"""
		Sets up the 3D visual window for the simulation. Calls several internal methods 
		that define and control the appearance of the window, comet, coma, camera FOV and 
		lighting. Assigns class attributes in preparation for running the simulation.	
		"""
		
		self.setWindow()
		self.loadComet()
		self.setLight()
		self.setAxes()
		self.setCamera()
		progress('Prepping Simulation')
		self.draw()
		progress()
		self.startZoom = self.axes.camera.zoom
		progress()
		self.startView = self.axes.GetView().get('loc')
		progress(p=True)
		self.pixSCX = self.getPS(321,1, 320, 1)[0]
		print('Ready to go!')
		return

	def setWindow(self):
		"""
		Instantiates the window the simulation will run in, acts as a wrapper for visvis 
		functions and classes. Sets sub-classes as attributes so they can be called for 
		use easily later.	
		"""

		progress('Opening simulation window')
		self.app = vv.use()
		progress()
		self.app.Create()
		progress()
		self.mainW = MainWindow()
		progress()
		self.mainW.resize(661, 501)
		progress()
		self.figure = vv.gcf()
		progress()
		self.axes = vv.gca()
		
		progress(p=True)
		return								

	def loadComet(self):
		"""
		Loads comet shape model from user provided path, converts model faces and 
		vertices into a visvis mesh object and stores as a class attribute. Instantiates 
		a coma() class object and stores the resulting point set object as a class 
		attribute. Displays comet mesh and point cloud coma in the current axes instance.	
		"""
		progress(f'Loading comet and coma')
		df = pd.read_csv(os.path.join(self.link, self.cometLink), delim_whitespace=True, names=['TYPE','X1', 'X2', 'X3'], engine='python', header=None, on_bad_lines='warn')									
		progress()		
		vertices = df.loc[df['TYPE'] == 'v'][['X1', 'X2', 'X3']].values.tolist()
		progress()
		faces = df.loc[df['TYPE'] == 'f'][['X1', 'X2', 'X3']].values
		progress()		
		faces = faces - 1
		progress()
		faces = faces.astype(int)
		progress()
		faces = faces.tolist()	
		progress()
		comet = cometMesh(faces, vertices)
		progress()
		coma_cl = coma()
		progress()
		self.comet = comet.obj
		progress()
		self.coma = coma_cl.cloud
		progress()
		del comet
		progress()
		l = vv.plot(self.coma, ms='.', mc = 'w', mw='1', mew = 0, ls='', axesAdjust=False , alpha=0.2)
		progress()
		self.axes.Draw()
		progress()
		self.figure.DrawNow()
		progress(p=True)
		return				
 
	def setLight(self):
		"""
		Interfaces with the visvis .axes.light class to turn off default window view 
		lighting, creates a new instance of the light class and stores as a class 
		attribute. Sets the light to be directional, coming from behind the "spacecraft".
		"""
		
		progress('Turning on the lights')
		self.axes.light0.Off()
		progress()
		self.lightObj = self.axes.lights[1]
		progress()
		self.lightObj.On()
		progress()
		self.lightObj.isDirectional = True
		progress()
		self.lightObj.position = (0.3,-0.8,-0.5, 0.0)
		progress(p=True)
		return					
		
	def setAxes(self):
		"""
		Adjusts .axes parameters to make the simulation look like is is happening in 
		space, i.e. setting background colour to black, hiding the axes, turning off 
		automatic aspect adjustments and setting the axes camera to 3D to allow x,y and 
		z translations of the simulation.		
		"""
		
		progress('Preparing the Axes')
		self.axes.bgcolor = (0, 0, 0)
		progress()
		self.axes.axis.showGrid = False
		progress()
		self.axes.axis.visible = False
		progress()
		self.axes.daspectAuto = False
		progress()
		self.axes.camera = '3d' 
		progress(p=True)
		return						

	def setCamera(self, f=0, az = 90, el=0):
		"""
		Adjusts the simulation "camera", i.e. where the spacecraft appears to be in 
		location w.r.t. the comet.	
		"""
		
		self.axes.camera.fov = f
		self.axes.camera.azimuth = az
		self.axes.camera.elevation = el
		return								

	def getPS(self, x, y, x2, y2):
		"""
		Determines the pixel scale of the current axes instance in the WCS from pixel 
		coordinates in the local coordinate system.
		"""
		
		return np.asarray(self.axes.camera.ScreenToWorld((x,y)))-np.asarray(self.axes.camera.ScreenToWorld((x2, y2)))
		
	def draw(self):
		"""
		Updates the current window.
		"""
		
		self.axes.Draw()
		self.figure.DrawNow()
		return
		
	def setView(self, time):
		"""
		Calls view() class for a given time in the encounter.
		"""
		if hasattr(self, 'view'):
			self.view.updateTime(time)
		else:
			self.view = view(time, self.interceptorVelocity, self.closeApproach)
		return	

	def movingView(self, time, step, params, rFrame='tiri', first=None, scanning=False):
		"""
		Advances the simulation in time for a given time step. An internal call of 
		.setView() returns the ifov and delta angle. For the first frame in the 
		simulation, internal calls to start\_zoom() and .setCamera() return the 
		baseline view parameters which are stored for use throughout. At all other 
		times, the comet and coma translations from the previous iteration are removed 
		and the axes reset to the default view. The encounter velocity and time step are 
		used to return the camera (i.e. spacecraft/tiri) location and comet translation, 
		respectively, for simulations run in the `tiri' rest frame, an internal call to 
		activeView translates the comet in the FOV to simulates scanning of the mirror 
		when operating mode recorded as 'moving'. The .view.angle and .view.ifov are used 
		to scale the zoom parameter to simulate the expected view.
		"""

		time_el = step
		ifov, ang, notes = params
		scx = ifov/(self.pixSCX*2)
		if first is True:
			z2 = self.startZoom/abs(scx)
			self.axes.camera.zoom = z2
			self.setCamera(az=ang)
		else:
			old_ang = self.angle
			self._removeTransformations()
			CI_dist = time_el*(self.interceptorVelocity-self.cometVelocity)
			com_dist = time_el*self.cometVelocity
			z2 = self.startZoom/abs(scx)
			self.axes.camera.Reset()
			if rFrame == 'tiri':
				self.axes.camera.roll = 90-ang
			self.axes.SetView({'zoom':z2})
			self.axes._wobjects[1]._transformations.insert(1, vv.Transform_Translate(dy = com_dist))
			self.axes._wobjects[2]._transformations.insert(1, vv.Transform_Translate(dy = com_dist))
			self.axes.Draw()
			self.setCamera(az = ang, el=0)
			loc = np.asarray(self.axes.GetView().get('loc'))
			loc[1] = CI_dist	
			self.axes.SetView({'loc': tuple(loc)})
			self.axes.Draw()
			if rFrame == 'tiri' and notes:
				if notes.get('comet loc'):
					self.activeView(notes, ang)
		self.figure.DrawNow()
		self.angle = self.axes.camera.azimuth
		return
	
	def activeView(self, notes, ang):
		"""
		When rest frame set to `tiri' and operating mode classed as 'moving', translates 
		the comet in the FOV to simulate movement of the pointing mirror to scan the 
		filter array across the comet nucleus. When operating mode classed as 'waiting', 
		translates the comet to simulate the comet entering the FOV from the 
		right hand side and moving into the centre of the FOV.
		"""
		self.axes.Draw()	
		offsety, offsetz = self.axes.camera.ScreenToWorld((320, 240))
		fmov = np.asarray(notes.get('comet loc'))/2
		y,z = self.axes.camera.ScreenToWorld(fmov)	
		self.axes._wobjects[1]._transformations.insert(0, vv.Transform_Translate(dy=y+offsety+3))
		self.axes._wobjects[2]._transformations.insert(0, vv.Transform_Translate(dy=y+offsety+3))
		if notes.get('mode') == 'waiting':
			self.axes.Draw()
			return
		self.axes._wobjects[2]._transformations.insert(0, vv.Transform_Translate(dx=z-offsetz))
		self.axes._wobjects[1]._transformations.insert(0, vv.Transform_Translate(dx=z-offsetz))
		self.axes.Draw()
		return

	def _removeTransformations(self):
		"""
		Removes transformations previously applied to the nucleus and coma models.
		"""
		for i in range(0,len(self.axes._wobjects[1]._transformations)):
			try:
				self.axes._wobjects[1]._transformations.remove(self.axes._wobjects[1]._transformations[0])
			except:
				pass
			try:
				self.axes._wobjects[2]._transformations.remove(self.axes._wobjects[2]._transformations[0])
			except:
				pass	
		self.axes.Draw()
		self.axes._wobjects[1]._transformations.insert(0,vv.Transform_Rotate(angle=100, ax=-1, ay=1, az=-1))
		self.axes.Draw()
		return	
	
	def closeWindow(self):
		"""
		Breaks the connection to the GUI app and closes the simulation window.
		"""
		self.mainW.close()
		return
		
class sim(window):
	
	"""
	Simulation interface, derivative of .window() class.	
	"""
	
	def __init__(self, cometObj, Opath, encounterV=70, cometV=15,temporalR=1/20, \
	rest = 'tiri', closestA=1000):
		"""		   
		#####  INPUT  ####
		cometObj	: str
			Path to comet shape model.
			   
		encounterV	: float or int
			Flyby relative velocity, defined in ESA-COMET-SYS-RS-001 as the modulus of 
			the spacecraft heliocentric velocity - comet heliocentric velocity. Stored 
			in class property eV.
			
		cometV	: float or int
			The desired comet velocity in km/s.
						
		closestA	: float or int
			LoS distance in km between the spacecraft and comet at closest approach. 
			Default is 1000km.	 

		temporalR	:	float
			Frame rate in Hz of the TIRI dectector during the encounter. 
			   
		Opath		: str 
			Path to directory where comet shape models are stored on disk. Simulation 
			output will be saved to this directory. 
			Stored in class property .link
		
		rest		: str
			Rest frame for simulation. The two available options are:
			1)	View defined in spacecraft rest frame, i.e. no rotation or scanning of 
				comet in FOV.
			2)	View defined in TIRI detector rest frame, i.e. with rotation and scanning
				of comet in FOV.
				
		"""
			

		self.cA = closestA
		self.eV = encounterV
		self.cV = cometV
		self.delta ='\u03B4'
		self.cpath = cometObj
		self.link = Opath
		self.timeRes = temporalR
		self.rFrame = rest
						
	def windowSetup(self):
		"""
		Instantiates the parent class and sets up the 3D visual window for the simulation.		
		"""
		super().__init__(self.eV, self.cV)
		self.whichComet(self.cpath)
		self.setup()
		return

	def setupDirs(self, rFrame, flag):
		"""
		Checks for and creates directories for simulation outputs and/or intermediate
		products as desired.
		"""
		folders = {
		'general': {'tiri': "TIRI_restframe", 'sc': "SC_restframe"}, 
		'saveRaw': {'tiri': "TIRI_restframe/RAW", 'sc': "SC_restframe/RAW"},
		'saveAnnotate' : {'tiri': "TIRI_restframe/PROCESSED", 'sc': "SC_restframe/PROCESSED"}
		}
	
		topLevel = os.path.join(self.link, f'simulationData_{np.round(self.timeRes, 4)}s')
		if not os.path.isdir(topLevel):
			os.mkdir(topLevel)
		if flag == 'general':
			midLevel = os.path.join(topLevel, folders.get(flag).get(rFrame))
			if not os.path.isdir(midLevel):
				os.mkdir(midLevel)		
			return midLevel
		if flag == 'saveRaw' or flag == 'saveAnnotate':
			subLevel = os.path.join(topLevel, folders.get(flag).get(rFrame))
			if not os.path.isdir(subLevel):
				os.mkdir(subLevel)		
			return subLevel

	def filterCentre(self, whichFilt):
		"""
		Constructs dictionary with line arrays at centre pixel location of a given filter.
		"""
		Edges = [self.fEdges.get(whichFilt).get(i) for i in ['vL', 'vR', 'hL', 'hR']]
		MidV = np.median([Edges[0], Edges[1]], axis=0)
		MidH = np.median([Edges[2], Edges[3]], axis=0)
		return {whichFilt: {'V': MidV, 'H': MidH}}

	def mapEncounter(self):
		"""
		Calulates geometry of encounter for each timestep 
		"""
		self.angleList = []
		self.distList = []
		self.ifovList = []
		self.indexList = []
		self.timeList = []
		
		encounterD = (self.cA/np.tan(np.radians(5)))*2
		tTotal = encounterD/self.eV
		timeSteps = np.round(np.arange(0, tTotal,self.timeRes),4)

		self.setView(0)
		for i,t in enumerate(timeSteps):
			self.setView(t)
			if self.view.angle >= 80:
				self.indexList.append(i)
				self.timeList.append(t)
				self.ifovList.append(self.view.ifov)
				self.distList.append(self.view.b)
				self.angleList.append(self.view.angle)
		return
		
	def mapCubes(self):
		"""
		Identifies which timesteps fall within the angular limits of the datacube 
		requirements.
		"""
		# DC1 @ 80-100 (90 +/- 10)
		# DC2 @ 95 - 125 (110 +/- 15)
		# DC3 @ 120 - 130 (125 +/- 5)
		# DC4 @ 160 - 170 (165 +/- 5)
		# DC 5 @ > 173
		# DC 6 @ > 173
		self.notes = [None]*len(self.angleList)
		self.defineFilts()
		self.filterOrder = self.fpa.order
		it = self.fpa.integrate

		self.ft = {}
		for key, i in it.items():
			self.ft[key]=math.ceil(i/self.timeRes)
		if self.rFrame == 'tiri':
			self.dcflag = None
		else:
			self.dcflag = 'DC1'
		self.lastAdjust = [0,0]
		self.vhflag = 'v'
		self.cometLoc = None
		j = 0 
		for i, a in enumerate(self.angleList):
			if i <= j and i != 0:
				continue	
			if a >= 80 and self.dcflag == None:
				j,a = self.firstEntry(j, a)	
				self.dcflag ='DC1'		
				continue
			
			if a >= 80 and self.dcflag == 'DC1':
				j, a = self.dataCube(i, a, j)
				self.dcflag = 'DC2'
				continue

			if a >= 95 and self.dcflag == 'DC2':	
				j, a = self.dataCube(i, a, j)
				self.dcflag = 'DC3'
				continue
		
			if a >= 120 and self.dcflag == 'DC3':
				j, a = self.dataCube(i, a, j)
				self.dcflag = 'DC4'
				continue
	
			if a >= 160 and self.dcflag == 'DC4':
				j, a = self.dataCube(i, a, j)
				self.dcflag = 'DC5'
				continue
	
			if a >= 173 and self.dcflag == 'DC5':
				j, a = self.dataCube(i, a, j)
				self.dcflag = 'DC6'
				continue

			if a >= 173 and self.dcflag == 'DC6':
				j, a = self.dataCube(i, a, j)
				self.dcflag = 'end'
				continue
		return		
		
	def dataCube(self,i, a, j):
		"""
		Maps out the operating mode of TIRI at each timestep within a datacube. Possible
		active modes are `waiting', `exposing' and `moving'. Returns None if TIRI is idle.
		"""
		for n, f in enumerate(self.filterOrder):
			if n == 0:
				self.cometLoc = [640, 480]
				j = i
				j, a = self.exposing(f, j)
			else:
				if self.rFrame == 'tiri':
					whichmove = [self.filterOrder[n-1], f]
					parkmove = [f, 'broadband']
				else:
					whichmove = [f,self.filterOrder[n-1]]
					parkmove = ['broadband', f]
				j, a = self.moving(whichmove, j, a)
				j, a = self.exposing(f, j)
				if n == len(self.filterOrder)-1:
					j, a = self.moving(parkmove, j, a, park=True)
		return j, a

	def firstEntry(self, j, a):
		"""
		Maps out at which timesteps TIRI is in `waiting' operating mode, i.e. the
		time between the comet first entering the FOV and reaching the centre of
		the FOV. The encounter velocity and detector read out rate are used to determine
		how many simulation frames it takes for the comet to reach the central pixel in 
		the FPA and the apparent pixel location of the comet in each frame as
		it travels. This mode is only available if the rest frame is set to `tiri'
		"""
		movingList = []
		frames = math.ceil(self.eV*self.timeRes)
		inc = (1280-640)/frames
		_ = [movingList.append(foo) for foo in np.arange(frames+1)+j]
		for nm, dm in enumerate(movingList):
			loc = 1280-(inc*(nm))
			self.notes[dm] = {'cube':None, 'filter':None, 'mode':'waiting', 'move f': None, 'comet loc': [loc, 960], 'orientation': None}
		return dm+1, self.angleList[dm+1]
		
	def moving(self, mov, j, a, park=False):
		"""
		Maps out at which timesteps TIRI is in `moving' operating mode, i.e. when the 
		pointing mirror is in motion. The encounter velocity and detector read out rate 
		are used to determine how many simulation frames it takes for the mirror to scan
		between filters, i.e. move the comet nucleus from the centre (in x-axis) of one
		filter to another, in the order defined in the FPA object. The apparent pixel 
		location of the comet in each frame as the mirror moves is also stored for later.
		"""
		movingList = []
		mft = (360/16)*self.timeRes #no. of degrees of movement in 1 frame.
		if self.vhflag == 'v':
			orient = 'V'
		else:
			orient = 'H'
		centres = [self.filterCentre(f)for f in mov]
		if self.rFrame == 'tiri':
			Arc = np.asarray(rotateArc(a+90, self.ox, self.oy))
			startMid = np.asarray(centres[0][mov[0]].get(orient))
			endMid = np.asarray(centres[1][mov[1]].get(orient))
		else:
			Arc = np.asarray(arc(self.ox, self.oy))
			startMid = np.asarray(self.rotateTime(centres[0][mov[0]].get(orient), a))
			endMid = np.asarray(self.rotateTime(centres[1][mov[1]].get(orient), a))
		
		arcLine =  LineString(Arc.T)
		lineStart = LineString(startMid.T)
		lineEnd = LineString(endMid.T)
		if self.checkIntersects(lineStart, arcLine) and self.checkIntersects(lineEnd, arcLine):
			startLoc = lineStart.intersection(arcLine)
			endLoc = lineEnd.intersection(arcLine)
		else:
			orient = 'H'
			self.vhflag = 'h'
			if self.rFrame == 'tiri':
				Arc = np.asarray(rotateArc(a+90, self.ox, self.oy))
				startMid = np.asarray(centres[0][mov[0]].get(orient))
				endMid = np.asarray(centres[1][mov[1]].get(orient))
			else:
				Arc = np.asarray(arc(self.ox, self.oy))
				startMid = np.asarray(self.rotateTime(centres[0][mov[0]].get(orient), a))
				endMid = np.asarray(self.rotateTime(centres[1][mov[1]].get(orient), a))
				arcLine =  LineString(Arc.T)
			lineStart = LineString(startMid.T)
			lineEnd = LineString(endMid.T)
			if self.checkIntersects(lineStart, arcLine) and self.checkIntersects(lineEnd, arcLine):
				startLoc = lineStart.intersection(arcLine)
				endLoc = lineEnd.intersection(arcLine)
			else:
				print(f'found no matches, j={j}, a={a}, i={i}, vhflag={self.vhflag}')	
		r = self.ox/np.sin((self.ox*0.26e-3)/2)
		hyp_pix = r*(np.arcsin(startLoc.x/r)-np.arcsin(endLoc.x/r))
		if self.rFrame == 'tiri':
			if orient == 'H':			
				hyp_pix = r*(np.arcsin(startLoc.y/r)-np.arcsin(endLoc.y/r))		
		hyp_deg = abs(hyp_pix)*(np.rad2deg(0.26e-3/2))
		hyp_frames = math.ceil(hyp_deg/mft)
		inc = hyp_pix/hyp_frames
		_ = [movingList.append(foo) for foo in np.arange(hyp_frames+1)+j]
		for nm, dm in enumerate(movingList):
			newa = self.angleList[dm]
			if self.rFrame == 'tiri':
				if orient == 'H':
					loc = startLoc.y-(inc*(nm))
					incL = np.squeeze(np.mgrid[centres[0][mov[0]].get('H')[0][0]:centres[0][mov[0]].get('H')[0][-1], loc:loc+1])						
				else:			
					loc = startLoc.x-(inc*(nm))
					incL = np.squeeze(np.mgrid[loc:loc+1, 0:(self.oy*2)])		
				arcN = np.asarray(rotateArc(newa+90, self.ox, self.oy))	
			else:
				arcN = Arc
				newMid = np.asarray(self.rotateTime(centres[0][mov[0]].get(orient), newa))
				incY = (endLoc.y-startLoc.y)/hyp_frames
				locX = newMid[0]-(inc*(nm))
				locY = newMid[1]+(incY*nm)
				incL = np.squeeze(np.stack([locX, locY]))
			arcLineN =	LineString(arcN.T)
			lineL = LineString(incL.T)			
			if self.checkIntersects(lineL, arcLineN):
				newLoc = lineL.intersection(arcLineN)
			else:
				print(f'no match for j= {j}, a= {a}, filter move = {mov}')
				break
			if self.rFrame == 'tiri':
				comL = [newLoc.x, newLoc.y]
				fChange = mov
				fAdjust = [0,0]
			else:
				comL = [640, 480]
				fChange = mov[::-1]
				miniAdjust = 0
				prevAdjust = self.lastAdjust
				if 'broadband' in mov:
					if 'A' in mov:
						prevAdjust = [0,0]
						miniAdjust = (9/self.ifovList[dm])	
				fAdjust = [newLoc.x-startLoc.x-miniAdjust+prevAdjust[0], newLoc.y-startLoc.y+prevAdjust[1]]
			if park:
				self.notes[dm] = {'cube':self.dcflag, 'filter':None, 'mode':'moving', 'filter change': 'parking', 'comet loc': comL, 'orientation': self.vhflag, 'filter adjust': fAdjust}
			else:
				self.notes[dm] = {'cube':self.dcflag, 'filter':None, 'mode':'moving', 'filter change': fChange, 'comet loc': comL, 'orientation': self.vhflag, 'filter adjust': fAdjust}
		self.lastAdjust = fAdjust
		self.cometLoc = comL
		return dm+1, self.angleList[dm+1]

	def exposing(self, f, j):
		"""
		Maps out at which timesteps TIRI is in `exposing' operating mode, i.e. when the 
		detector is integrating. Individual dwell times for each filter are 
		defined in the FPA object as the total integration time required to reach 
		SNR=100 and this is translated into  using the detector read out 
		rate and rounding up to an integer value for all partial frames to determine how 
		many simulation frames each filter requires in exposing mode and which timesteps
		each filter is exposing over.
		"""
		exposingList = []
		dwell = self.ft.get(f)
		_ = [exposingList.append(foo) for foo in np.arange(dwell)+j]
		for di in exposingList:
			self.notes[di] = {'cube':self.dcflag, 'filter':f, 'mode':'exposing', 'filter change': None, 'comet loc': self.cometLoc, 'orientation': self.vhflag, 'filter adjust': self.lastAdjust}
		return di+1, self.angleList[di+1]
		
	def addAnnotate(self, a, t, ifov, notes, ind, ff):
		"""
		Annotates simulation frames with TIRI filter array and visual info if desired.		
		"""
		rArc = np.asarray(rotateArc(a+90, self.ox, self.oy))
		ax.set_xticks([])
		ax.set_yticks([])
		fig.tight_layout(pad=0)
		ax.set_xlim(0,1280)
		ax.set_ylim(0,960)
		mxy = [0,0]
		if self.rFrame == 'tiri':
			rAngle = 90
			labelL = [0.90,0.95]
		else:
			rAngle = a
			labelL = [0.90,0.15]
			if notes:
				 mxy = notes.get('filter adjust')
		for i in self.filterOrder:
			if i == 'broadband':
				p = self.createPoly(*mxy, i, rAngle, 'h')
				ax.plot(*p.exterior.xy, lw=0.9)
			else:
				for j in ['v', 'h']:
					p = self.createPoly(*mxy, i, rAngle, j)
					ax.plot(*p.exterior.xy, lw=0.9)
		ax.text(*labelL, f't = {np.round(t)}s\n{self.delta} = {np.round(a,1)}\n{np.round(ifov, 2)}km/pix', color='k', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'), va='top', ha='center', transform = ax.transAxes)
		if notes:
			mode = notes.get('mode')	
			if mode == 'moving':
				at = AnchoredText(f'MIRROR MOVING', loc = 'upper center', frameon=False, pad=0,prop=dict(color='b', size=20))
				ax.add_artist(at)
			elif mode == 'exposing':
				sp = self.createPoly(*mxy, notes.get('filter'), rAngle, 'h' if notes.get('filter')=='broadband' else notes.get('orientation'))
				ax.fill(*sp.exterior.xy, color='pink', alpha=0.3)		
			if mode != 'waiting' and self.rFrame == 'tiri':
				sns.lineplot(x = rArc[0], y=rArc[1], ax=ax, alpha=0.2, color='yellow', linestyle='--')
		else:
			if ff!=0:
				patch = patches.PathPatch(path, facecolor='yellow', edgecolor='white', lw=1.5, zorder=10)
				ax.add_patch(patch)
				ax.text(170,87.5, f'x{ff}', color='yellow', fontsize='x-large', family='fantasy', va='center', ha='left', variant='small-caps', path_effects=[pe.withStroke(linewidth=3, foreground="k")])
		fig.show()
		fig.tight_layout(pad=0)
		return

	def flyby(self, fast=20, annotate=True, saveRaw=False, saveAnnotate=True):
		"""
		Runs the simulation.
		"""
		global fig
		global ax
		fig, ax = plt.subplots()
	
		print(f'+{"-"*(wdth-2)}+')
		print('Simulating encounter...'.center(wdth))

		if saveRaw or saveAnnotate:
			self.outputPath = self.setupDirs(self.rFrame, 'general')
		if saveRaw:
			rawPath = self.setupDirs(self.rFrame, 'saveRaw')
		if saveAnnotate:
			annotatePath = self.setupDirs(self.rFrame, 'saveAnnotate')
	
		self.mapEncounter()
		self.mapCubes()
		self.dropflag = 0
		self.outputImages = []
		tempView = view(0, 70, 1000)
		self.movingView(0, 0, [tempView.ifov, tempView.angle, None], rFrame=self.rFrame, first=True)
		self.draw()
		dTracker = 'DC1'	
		for i, t in enumerate(self.timeList):
			if self.notes[i] and self.notes[i].get('cube') == dTracker:
				print(f'Data Cube at {self.delta} = {np.round(self.angleList[i])}')
				dTracker = f'DC{int(dTracker[-1])+1}'
			params = [self.ifovList[i], self.angleList[i], self.notes[i]]
			if self.notes[i]:
				self.movingView(t, self.timeRes, params, rFrame=self.rFrame, scanning=True)
			elif fast!= 0 and self.dropflag%fast == 0:
				self.movingView(t, self.timeRes, params, rFrame=self.rFrame)		
				self.dropflag+=1
			else:
				self.dropflag+=1
				continue
			self.draw() 
			temp_image = vv.getframe(vv.gca())

			if saveRaw:
				np.save(os.path.join(rawPath, f'{np.round(self.angleList[i], 3)}raw.npy'), np.asarray(temp_image))
				if saveOutput and not annotate:
					self.outputImages.append(temp_image)
		
			if annotate:
				ax.clear()
				ax.imshow(temp_image)
				ax.set_ylim(0,960)
				ax.set_xlim(0,1280)
				self.addAnnotate(self.angleList[i], t, self.ifovList[i], self.notes[i], i, fast)
				fig.canvas.draw()		
				data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
				data = data.reshape(tuple(np.asarray(fig.canvas.get_width_height()[::-1])*2) + (3,))
				if saveAnnotate:
					np.save(os.path.join(annotatePath, f'{np.round(self.angleList[i], 3)}A.npy'), np.asarray(data))
				self.outputImages.append(data)
		
		self.closeWindow()
		plt.close(fig)
		plt.close('all')
		print('...Encounter Complete'.center(wdth))
		return
	
	def outputTable(self):
		"""
		Creates table with the following entries: time since start of encounter, distance 
		between comet interceptor and target comet along the line of sight, delta angle, 
		pixel scale and operating mode.	
		"""		
		datatab = QTable()					
		datatab['Time'] = self.timeList
		datatab['Time'].unit = u.s
		datatab['Delta Angle'] = np.round(self.angleList,3)
		datatab['Delta Angle'].unit = u.deg
		datatab['Distance'] = np.round(self.distList, 3)
		datatab['Distance'].unit = u.km
		datatab['Pixel Scale'] = np.round(self.ifovList,4)
		datatab['Pixel Scale'].unit = u.km/u.pix
		mode = MaskedColumn([m.get('mode').upper() if m else np.nan for m in self.notes], mask=[False if m else True for m in self.notes])
		datatab['Mode'] = mode
		return datatab
						
	def output(self, fRate):
		"""
		Outputs .csv and .gif products to local disk if desired. 		
		"""
		self.setupDirs(self.rFrame, 'general')
		print(f'+{"-"*(wdth-2)}+')
		print('| Assembling Outputs |', end='')
		saveName = f'flyby-{(self.rFrame).upper()}-{np.round(self.timeRes, 4)}'
		arr = self.outputTable()
		aIms = self.outputImages
		print(' Saving Outputs|') 
		ascii.write(arr, f'{os.path.join(self.outputPath,saveName)}.csv', format='csv', overwrite=True)
		imageio.mimsave(os.path.join(self.outputPath, f'{saveName}.gif'), aIms, duration=fRate)
		return
	
	def defineFilts(self):
		"""
		Assigns class attributes containing the locations of filter boundaries in pixels.		
		"""
		self.fpa = FPA()
		self.fpa.setScale(2)
		self.fEdges = self.fpa.filterOverlay()
		self.ox, self.oy = np.asarray(self.fpa.xy_dim)/2
		return
											
	def rotateTime(self, filterEdge, angle):
		"""
		Takes in an array of coordinates defining the pixel locations of the TIRI filter 
		boundaries edges and returns coordinates rotated to compensate for pointing
		mirror rotation at the given time.	
		"""
		
		angle = 90 - angle
		rFilter = rotateFilters(filterEdge, angle, self.ox, self.oy)
		return np.asarray(rFilter)
		
	@staticmethod
	def checkIntersects(l1, a1):
		"""
		Checks for an intersection between a straight line and an arc.
		"""		
		return l1.intersects(a1)

	def createPoly(self, mx, my, whichFilt, angle, flag):
		"""
		Returns a polygon object with the spatial parameters (edge, rotation, location, 
		etc) of a specified filter at a given time. Takes into account rotation and (x,y) 
		translation in spacecraft rest frame as a result of the movement of the TIRI 
		pointing mirror.		
		"""
			
		if flag == 'h':
			LEdge = self.fEdges.get(whichFilt).get('hL')
			REdge = self.fEdges.get(whichFilt).get('hR')
		else:
			LEdge = self.fEdges.get(whichFilt).get('vL')
			REdge = self.fEdges.get(whichFilt).get('vR')
		rfiltL = self.rotateTime(LEdge, angle)
		rfiltR = self.rotateTime(REdge, angle)
		if self.rFrame == 'tiri':
			sign = 1
		else:
			sign = -1
		lx,ly = np.asarray(rfiltL[0])-sign*mx, np.asarray(rfiltL[1])-sign*my
		rx, ry = np.asarray(rfiltR[0])-sign*mx, np.asarray(rfiltR[1])-sign*my
		polygon = Polygon([(lx[0], ly[0]), (lx[1], ly[1]),(rx[::-1][0], ry[::-1][0]),\
		(rx[::-1][1], ry[::-1][1])])
		return polygon

class view:
	def __init__(self, time, encounterVelocity, closestApproach):
		self.time = time
		self.eV = encounterVelocity
		self.ca = closestApproach
		self.initialDist = self.ca/np.tan(np.radians(5))
		self.getView()

	def updateTime(self, time):
		"""
		Updates the view to reflect a different time in the encounter.
		"""
		self.time = time
		self.getView()
		return

	def getView(self):
		"""
		Returns geometry of encounter for the stored time.
		"""
		self.travelDistance()
		self.timeView()
		self.viewDistance() 
		self.viewPixscale() 
		return 

		
	def travelDistance(self):
		"""
		Calculates the distance between the comet and spacecraft along the direction
		of travel.
		"""
		ci_dist = self.time*(self.eV)
		self.a = self.initialDist-ci_dist	
		return

	def timeView(self):
		"""
		Calculates the current delta angle between the spacecraft and the comet.	
		"""
		if self.a == 0:
			ang = 90
		else:
			ang = np.rad2deg(np.arctan(self.ca/self.a))
			if self.a < 0:
				ang = (90+(90+ang))
		self.angle = ang
		return

	def viewDistance(self):
		"""
		Calculates the current distance between the spacecraft and the comet. 		
		"""
		self.b = self.ca/np.sin(np.radians(self.angle))
		return	

	def viewPixscale(self):
		"""
		Calculates the ifov, i.e. the km/pixel scale observed by TIRI at the comet 
		nucleus surface for a given delta angle and distance.		
		"""
		self.ifov = np.tan((0.26e-3)/2)*self.b*2
		return 


"""
#########################################################################################
										MAIN
#########################################################################################
"""

if __name__ == "__main__":
	_p = parseArguments()
	strhelp = _p.format_help()
	if len(sys.argv) < 1:
		raise SystemExit(print('Incorrect usage, see below:\n' +strhelp))
	else:
		args = _p.parse_args()
		print('_' * wdth)
		print("{:^{w}}".format(
	r""" 
 _ \    \    _ \    \     \  |  __| __ __|  __|  _ \   __| _) 
 __/   _ \     /   _ \   |\/ |  _|     |    _|     / \__ \    
_|   _/  _\ _|_\ _/  _\ _|  _| ___|   _|   ___| _|_\ ____/ _)"""   , w=wdth))
														   
		print(f'+{"-"*(wdth-2)}+')
		print("{:<{w1}} {:<2} {:<{w2}}".format("REST FRAME", ' ', f'{"spacecraft" if args.rFrame == "sc" else "tiri"}',w1=int(0.3*wdth), w2=int(0.4*wdth)))
		print("{:<{w1}} {:<2} {:<{w2}}".format("ENCOUNTER VELOCITY", ' ', f"{args.V} km/s", w1=int(0.3*wdth), w2=int(0.4*wdth)))
		print("{:<{w1}} {:<2} {:<{w2}}".format("COMET VELOCITY", ' ', f"{args.Vc} km/s", w1=int(0.3*wdth), w2=int(0.4*wdth)))
		print("{:<{w1}} {:<2} {:<{w2}}".format("CLOSEST APPROACH DISTANCE", ' ', f"{args.ca} km", w1=int(0.3*wdth), w2=int(0.4*wdth)))
		print("{:<{w1}} {:<2} {:<{w2}}".format("TIRI DETECTOR READOUT FREQUENCY", ' ', f"{args.tRes} Hz", w1=int(0.3*wdth), w2=int(0.4*wdth)))
		print("{:<{w1}} {:<2} {:<{w2}}".format("FAST FORWARD", ' ', f"{args.upSpeed}x", w1=int(0.3*wdth), w2=int(0.4*wdth)))
		print("{:<{w1}} {:<2} {:<{w2}}".format("OUTPUT DIRECTORY", ' ',args.outLink, w1=int(0.3*wdth), w2=int(0.4*wdth)))
	
		print(f'+{"-"*(wdth-2)}+')
		print("{:^{w}}".format(
	r""" 
  _ \   _ \ __ __| _ _|   _ \    \ |   __| _) 
 (   |  __/    |     |   (   |  .  | \__ \    
\___/  _|     _|   ___| \___/  _|\_| ____/ _)"""   , w=wdth))
		print("{:<{w1}} {:<2} {:<{w2}}".format("FRAME ANNOTATIONS", ' ', f"{'ON' if args.annotate else 'OFF'}{', saved to disk' if args.saveAnn else ''}", w1=int(0.3*wdth), w2=int(0.4*wdth)))			
		if args.saveOutput:
			print("{:<{w1}} {:<2} {:<{w2}}".format("GIF FPS", ' ', f'{args.fRate} fps', w1=int(0.3*wdth), w2=int(0.4*wdth)))
			print(f'+{"-"*(wdth-2)}+')
			print("{:^{w}}".format(
r""" 		
  _ \   |  | __ __|  _ \  |  | __ __|   __| _) 
 (   |  |  |    |    __/  |  |    |   \__ \    
\___/  \__/    _|   _|   \__/    _|   ____/ _)""", w=wdth)) 			
			print("{:<{w1}} {:<2} {:<{w2}}".format("GIF FILENAME", ' ', f"flyby-{(args.rFrame).upper()}-{np.round(args.tRes, 4)}.gif", w1=int(0.3*wdth), w2=int(0.4*wdth)))
			print("{:<{w1}} {:<2} {:<{w2}}".format("TABLE FILENAME", ' ', f"flyby-{(args.rFrame).upper()}-{np.round(args.tRes, 4)}.csv", w1=int(0.3*wdth), w2=int(0.4*wdth)))
		if args.saveRaw:
			print("{:<{w1}} {:<2} {:<{w2}}".format("RAW FRAMES", ' ', f"{'TIRI_restframe' if args.rFrame == 'tiri' else 'SC_restframe'}/RAW/%%%.npy", w1=int(0.3*wdth), w2=int(0.4*wdth)))
		if args.saveAnn:
			print("{:<{w1}} {:<2} {:<{w2}}".format("ANNOTATED FRAMES", ' ', f"{'TIRI_restframe' if args.rFrame == 'tiri' else 'SC_restframe'}/PROCESSED/%%%.npy", w1=int(0.3*wdth), w2=int(0.4*wdth)))						
		print(f'+{"-"*(wdth-2)}+')
		
	model = sim(args.mesh, args.V, args.Vc, os.path.join(args.outLink,''), args.tRes, \
	rest=args.rFrame, closestA=args.ca)
	model.windowSetup()
	model.flyby(args.upSpeed, args.annotate, args.saveRaw, args.saveAnn)
	if args.saveOutput:
		model.output(args.fRate)

