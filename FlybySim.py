#!/usr/bin/env python3
# +=====+========+==============================================+======================================+
# |Issue| Date   | Description Of Change                        | Comment                              |
# +=====+========+==============================================+======================================+
# | 1.0 |12/06/23| Initial release for i-CDR                    |                                      |
# +-----+--------+----------------------------------------------+--------------------------------------+
# | 2.0 |23/10/23| Included expected projection of filters on   | Docstrings for classes and functions |
# |     |        | FPA and scanning motion with respect to data-| have been updated where necessary.   |
# |     |        | cube requirements and estimates of exposure  |                                      |
# |     |        | time per filter.                             |                                      |
# |     |        | Adjusted output product options to include:  |                                      |
# |     |        | 1) gif including datacube scanning motion    |                                      |
# |     |        |      w.r.t spacecraft rest frame.            |                                      |
# |     |        | 2) csv file with tabulated time, delta       |                                      |
# |     |        |      angle, distance and pixel scale values  |                                      |
# |     |        |      associated with each simulation frame.  |                                      | 
# +-----+--------+----------------------------------------------+--------------------------------------+
# |     |        |                                              |                                      |
# +-----+--------+----------------------------------------------+--------------------------------------+
# |     |        |                                              |                                      |
# +-----+--------+----------------------------------------------+--------------------------------------+
# |     |        |                                              |                                      |
# +-----+--------+----------------------------------------------+--------------------------------------+
# |     |        |                                              |                                      |
# +-----+--------+----------------------------------------------+--------------------------------------+

from astropy.io import ascii
from astropy.table import Table
from astropy import units as u
import imageio
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.patheffects as pe
import numpy as	 np
import os
import pandas as pd
import pathlib
from PyQt5.QtWidgets import QWidget, QHBoxLayout
import seaborn as sns
from shapely import LineString
from shapely.geometry.polygon import Polygon
import urllib.request
import visvis as vv
from visvis import Pointset
import time

__version__ = '2.0'
__date__ = '2023-10-23'
__author__ = 'Maisie Rashman'

"""
#########################################################################################
									USAGE INFO
#########################################################################################
Comet download and display written using Thomas Albin's Space Science Tutorials - Part 12:	
https://github.com/ThomasAlbin/SpaceScienceTutorial

Code explicitly written to be run interactively, ipython in particular but should work in 
a jupyter notebook too. Uses a deprecated visual engine (visvis) that crashes all the 
time and takes up a lot of RAM

Comet 67P/Churyumov–Gerasimenko is the default
Comet velocity default = 15km/s
Flyby velocity vector default = 70km/s

Comet meshes have to be found and downloaded prior to using. Most available comet shape
models can be found here: https://sbn.psi.edu/pds/shape-models/
Example download of c67P mesh:

directory_path = 'directory to wherever you want to keep the file'
download_url = 'https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/rosetta
				/raw/misc/models/ROS_CG_M001_OSPCLPS_N_V1.OBJ'
download_file(directory_path, download_url)


#########################################################################################
										FUNCTIONS
#########################################################################################
"""

def downloadFile(dl_path, dl_url):
	"""
	Downloads files from url and saves to directory.
	
	####  INPUT	 ####
	DL_PATH : str
			  Download path on the local machine, relative to this function.
	
	DL_URL	: str
			  Download url of the requested file.		
	"""
	file_name = dl_url.split('/')[-1]
	pathlib.Path(dl_path).mkdir(parents=True, exist_ok=True)
	if not os.path.isfile(dl_path + file_name):
		urllib.request.urlretrieve(dl_url, dl_path + file_name)
	return
	
def makecDict(cometNames, cometApehelion, cometPerihelion, cometFiles):
	"""
	Creates a dictionary of comets with available shape models, alongside their velocity
	at perihelion or 1AU, whichever is greatest, and the name of the shape model file.
	
	####  INPUT	 ####
	cometNames		: list, str
			Name/identifiers of comets.
	
	cometAphelion	: list, float
			Aphelion of comets in AU.
				
	cometPerihelion : list, float
			Perihelion of comets in AU.
	
	cometFiles		: list, str
			File names of comet shape models.
				
	####  OUTPUT  ####	   
	cometDict		: dict
			Nested dictionary in form {name: {velocity, file_name}}.
			
	"""
	
	if len(cometNames) == len(cometApehelion) == len(cometPerihelion) == len(cometFiles):
		cometDict = {}		
		for i in range(0, len(cometNames)):
			cp = cometPerihelion[i]
			ap = cometApehelion[i]
			if cp < 1:
				r = 1
			else:
				r = cp
			sub_dict = {'velocity': comVel(ap, cp, r), 'file': cometFiles[i]}
			cometDict[cometNames[i]] = sub_dict
	else:
		print('Missing info - parameter lists are not the same length')
	return cometDict

def comVel(apehelion, perihelion, sun_comet_distance):
	"""
	Calculates the velocity of a comet at a given distance from the sun.
	
	####  INPUT	 ####
	aphelion			: float
			Aphelion of comet in AU
	
	perihelion			: float
			Perihelion of comet in AU.
				
	sun_comet_distance	: float
			Distance of comet from the sun in AU.
				
	####  OUTPUT  ####	   
	v			: float
			Velocity of comet in km/s at sun_comet_distance. 
					
	"""
	AU = 1.495978707e11
	ap = apehelion*AU 
	per = perihelion*AU 
	a = (per+ap)/2
	M = 2e30
	G = 6.67430e-11
	v = np.sqrt(G*M*((2/(sun_comet_distance*AU))-(1/a)))
	return v*1e-3

def progress(txt=None,p=None, q=None):
	if txt != None:
		print(txt,end="", flush=True)
		time.sleep(0.2) 
	elif q != None:
		ind = q[0]
		total = q[1]-1
		percent = np.round((q[0]/q[1])*100, 2)
		if ind != total:
			print(f'\r {percent}%', end='') 
		else:
			print(f'\r {percent}%') 
	else:
		if p == 1:
			print('.\n', end="\n")
			time.sleep(0.2)		
		else:
			print('.',end="", flush=True)
			time.sleep(0.2)
	return

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
	theta = np.linspace(-1*np.radians(12), np.radians(12), 1000)
	y = r*np.cos(theta) - (r-oy)
	x = r*np.sin(theta) + ox
	cosang, sinang = np.cos(np.radians(angle)), np.sin(np.radians(angle))	
	xr = [ (x[i]-ox)*cosang-(y[i]-oy)*sinang+ox for i in range(0,len(x))]
	yr = [ (x[i]-ox)*sinang-(y[i]-oy)*cosang+oy for i in range(0,len(x))]	
	return xr,yr
	
def arc(ox, oy):
	"""
	Returns the pixel coordinates of a ~9 deg circular arc centred on the TIRI FOV.		
	
	####  INPUT  ####
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
	theta = np.linspace(-1*np.radians(12), np.radians(12), 1280)
	y = r*np.cos(theta) - (r-oy)
	x = r*np.sin(theta) + ox 
	return x,y
	
"""
#########################################################################################
										CLASSES
#########################################################################################
"""
class MainWindow(QWidget):
	def __init__(self, *args):
		QWidget.__init__(self, *args)
		self.fig = vv.backends.backend_pyqt5.Figure(self)
		self.sizer = QHBoxLayout(self)	   
		self.sizer.setContentsMargins(0,0,0,0)
		self.sizer.addWidget(self.fig._widget, 2)
		self.setWindowTitle('Comet Interceptor')
		self.show()

class view:
	"""
	Calculates the location of the spacecraft w.r.t. the comet at a certain point in time, 
	for use to set the view of the comet in the simulation
	"""
	
	def __init__(self, time, cometVelocity, interceptorVelocity):
	
		"""
		#####  INPUT  ####
		time			=  Float value in seconds of time elapsed since the start of fly-by.  
				   
		cometVelocity	=  Float or integer in km/s of the desired velocity the simulated 
						   comet is travelling at. 
	
		interceptorV	=  Float or integer in km/s of the desired flyby relative velocity, 
						   defined in ESA-COMET-SYS-RS-001 as the modulus of the spacecraft 
						   heliocentric velocity - comet heliocentric velocity. For 
						   simplicity, the spacecraft velocity will be taken as 
						   interceptorV-cometVelocity.
	
		#####  PROPERTIES ###
		.time			= set to time input
		.cometV			= set to cometVelocity input
		.iV				= set to interceptorVelocity input	 
		.angle			= Delta angle in degrees returned by .time_view()	
		.distance		= Seperation in km returned by .view_distance()
		.ifov			= ifov in km/pixel returned by .view_ifov().	
		"""
	
		self.time = time
		self.cometV = cometVelocity
		self.iV = interceptorVelocity
		self.angle = self.timeView()
		self.distance = self.viewDistance()
		self.ifov = self.viewPixscale()

	def timeView(self):
		"""
		Calculates the current delta angle between the spacecraft and the comet. Does 
		this over four calculations:
		   1) the maximum distance between the two objects, i.e. at the start 
			  of flyby defined as when the delta angle = 5 degrees. Distance 
			  defined as the adjacent side to the delta angle in a right- 
			  angled triangle with shortest (opposite) side equal to closest 
			  approach - 1000km. 
		   2) the total distance moved by each object since the start of flyby 
			  using time elapsed and object velocity. 
		   3) the distance between the two objects by reducing the maximum 
			  distance by the total distance moved by each object. 
		   4) the delta angle using a right angled triangle defined with 
			  opposite side equal to closest approach - 1000km and adjacent 
			  side equal to the value calculated in 3).
		Angles after closest approach are corrected to be continuous.
		   
		##### EQUATION ####
		tan(theta) = opposite/adjacent
		delta = arctan(1000/distance)

		#### OUTPUT ####
		ang = Float in degrees
		"""
		
		ini_dist = 1000/np.tan(np.radians(5))
		ci_dist = self.time*(self.iV-self.cometV)
		com_dist = self.time*self.cometV
		curr_dist = ini_dist-ci_dist-com_dist
		ang = np.rad2deg(np.arctan(1000/curr_dist))
		if curr_dist < 0:
			ang = (90+(90+ang))
		return ang
		
	def viewDistance(self):
		"""
		Calculates the current distance between the spacecraft and the comet. Where 
		distance is defined as the hypotenuse of the right-angled triangle with shortest 
		side (opposite) equal to closes approach - 1000km and angle is equal to the 
		delta angle returned by timeview()

		##### EQUATION ####
		sin(theta) = opposite/hypotenuse
		distance = 1000/sin(delta)
									
		#### OUTPUT ####
		Float in km along the line of sight between the s/c and comet
		"""
		
		a = self.angle
		return 1000/np.sin(np.radians(a))	
		
	def viewPixscale(self):
		"""
		Calculates the ifov, i.e. the km/pixel scale observed by TIRI at the comet 
		nucleus surface for a given delta angle and .viewDistance().

		##### EQUATION ####
		tan(theta) = opposite/adjacent
		ifov = 2*(tan(delta/2)*seperation
		
		#### OUTPUTS ####
		float in km/pixel 
		"""
		
		return np.tan((0.26e-3)/2)*self.viewDistance()*2	
 
class cometMesh:
	
	"""
	Creates a 3D mesh from a shape model and applies transformations to get the best 
	viewing angle (for no purpose but ~aesthetics~). Mostly just a redundant wrapper for
	visvis classes and methods. 
	"""
	
	def __init__(self, faces, vertices):
		"""		
		#####  INPUTS ####
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

class coma: 
	"""
	Creates a very inauthentic representation of a coma. Just for visuals, not 
	scientifically relevant.
	"""
	
	__slots__ = ("cloud")
	
	def __init__(self): 
		"""
		##### PROPERTIES ####
		.cloud		=  Point cloud representing the coma returned by .generateComa()
		"""
		
		self.cloud = self.generateComa()
		
	def generateComa(self):
		"""
		Creates a 29000 point cloud composed of four different parts (for no purpose but
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
		
class window:
	"""
	Parent class wrapper for all things simulation.
	"""
	
	def __init__(self):
		"""		
		##### DEFAULTS #####
		.name				 =	String containing the name of the comet to use in the 
								simulation, this must match the name in the comet 
								dictionary.
								Default is 'c67p'.
								Can be adjusted with X.name = string	   

		.cometVelocity		 =	Float or integer of the desired velocity the simulated 
								comet will travel at. Informs the view of the comet for 
								each timestep.
								Default is 15km/s
								Can be adjusted with X.cometVelocity = float
							
		.interceptorVelocity =	Float or integer of the desired flyby relative velocity, 
								defined in ESA-COMET-SYS-RS-001 as the modulus of the 
								spacecraft heliocentric velocity - comet heliocentric 
								velocity. For simplicity, the spacecraft velocity will be 
								taken as interceptorVelocity minus cometVelocity for all 
								calculations. Informs the view of the comet for each 
								timestep.
								Default is 70km/s
								Can be adjusted with X.interceptorVelocity = float
		"""
		
		self.name = 'c67p'
		self.cometVelocity = 15
		self.interceptorVelocity = 70
		
	def whichComet(self, cometDict):
		"""
		Assigns class attributes based on a search of the given comet dictionary using 
		.name attribute as the search identifier.

		#####  INPUT ####
		cometDict		=  Dictionary of available comet shape models. Needed to identify 
						   the parameters and locate the shape models to create the comet 
						   mesh for the simulation.
			
		#### OUTPUT ####
		.cometLink		=  String containing the file name of the comet shape model for 
						   use in the simulation.
		.cometVelocity	=  Float or integer of the desired velocity the simulated comet 
						   will travel at. Informs the view of the comet for each 
						   timestep.		
		"""
		
		progress(f'Preparing flyby window for {self.name}')
		self.cometLink = cometDict.get(self.name).get('file')
		progress()
		self.cometVelocity = cometDict.get(self.name).get('velocity')
		progress(p=1)
		print(f'Location of shape model: {self.link+self.cometLink}')
		print(f'Setting comet velocity as {self.cometVelocity}')
		return
		
	def setup(self):
		"""
		Sets up the 3D visual window for the simulation. Calls several internal methods 
		that define and control the appearance of the window, comet, coma, camera FOV and 
		lighting. Assigns class attributes in preparation for running the simulation.
		
		#### OUTPUT ####
		.startZoom	=	Stores a baseline scaling for simulated spaceraft to comet 
						distance by querying the default camera view set in local 
						coordinates by the visvis .axes.camera class. This is used in the 
						simulation to determine the needed wcs adjustment of the viewing 
						parameters (i.e. camera, location, zoom, angle) for each timestep 
						to simulate the flyby.
				
		.startView	=	Stores a baseline location of the simulated spaceraft w.r.t. the 
						comet by querying the default camera location set by the visvis 
						.axes._view class. This is used in the simulation to determine 
						the needed wcs adjustment of the viewing parameters (i.e. camera, 
						location, zoom, angle) for each timestep to simulate the flyby.
				
		.pixSCX		=	Stores the default pixel scale of the window in local 
						coordinates. This is used in the simulation to determine the 
						needed wcs adjustment of the viewing parameters (i.e. camera, 
						location, zoom, angle) for each timestep to simulate the flyby.		
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
		progress(p=1)
		self.pixSCX = self.getPS(321,1, 320, 1)[0]
		print('Ready to go!')
		return

	def setWindow(self):
		"""
		Instantiates the window the simulation will run in, acts as a wrapper for visvis 
		functions and classes. Sets sub-classes as attributes so they can be called for 
		use easily later.
		
		#### OUTPUT ####
		.axes	=	Stores the visvis .axes class.
		.figure =	Stores the visvis window wrapper.		
		"""
		%gui qt
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
		progress(p=1)
		return								

	def loadComet(self):
		"""
		Loads comet shape model from .link and .cometLink, converts model faces and 
		vertices into a mesh object, stored as a class attribute. Instantiates a coma() 
		class object and stores the resulting point set object as a class attribute. 
		Displays comet mesh and point cloud coma in the current axes instance.
		
		#### OUTPUT ####  
		.coma	=	random point cloud representation of a coma. Not accurate
					just for visuals. 
		.comet	=	3D comet mesh object.			
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
		progress(p=1)
		return				
 
	def setLight(self):
		"""
		Interfaces with the visvis .axes.light class to turn off default window view 
		lighting, creates a new instance of the light class and stores as a class 
		attribute. Sets the light to be directional, coming from behind the "spacecraft".
		
		#### OUTPUT ####
		.lightObj	=  Stores an instance of axes.lights to control the lighting in the 
					   simulation.	
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
		progress(p=1)
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
		progress(p=1)
		return						

	def setCamera(self, f=0, az = 90, el=0):
		"""
		Adjusts the simulation "camera", i.e. where the spacecraft appears to be in 
		location w.r.t. the comet.	
		
		#####  INPUT  ####
		f	=	Integer value of the projection of the axes, can be set to adjust the 
				apparent focal length of the simulation, i.e. the depth of field.
				Default is 0

		az	=	The azimuth angle of the "spacecraft". Can be set to adjust the face-on 
				viewing angle of the comet.
				Default is 90

		el	=	The elevation of the "spacecraft". Can be set to adjust the face-on 
				viewing angle of the comet. Default is 0	
		"""
		
		self.axes.camera.fov = f
		self.axes.camera.azimuth = az
		self.axes.camera.elevation = el
		return								

	def getPS(self, x, y, x2, y2):
		"""
		Determines the pixel scale of the current axes instance in the WCS from pixel 
		coordinates in the local coordinate system.
		
		####  INPUT	 ####
		x	=	x-coordinate of a pixel.

		y	=	y-coordinate of a pixel.

		x2	=	x-coordinate of neighbouring pixel.

		y2	=	y-coordinate of neighbouring pixel.
		
		####  OUTPUT  ####
		The pixel scale of the current axes instance in WCS 
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
					
		#####  INPUT  ####
		time	=  Float value of time elapsed since beginning of fly-by in seconds.
			
		####  OUTPUT  ####
		.view	= instance of view() class
		"""		
		self.view = view(time, self.cometVelocity, self.interceptorVelocity)
		return		

	def movingView(self, time, step, first=None):
		"""
		Advances the simulation in time for a given time step. An internal call of 
		.setView() returns the ifov and delta angle. For the first frame in the 
		simulation, .internal calls to start_zoom() and .setCamera() return the baseline 
		view parameters which are stored for use throughout. At all other times, the 
		comet and coma translations from the previous iteration are removed and the axes
		reset to the default view. The encounter velocity and time step are used to 
		return the camera (i.e. spacecraft) location and comet translation, respectively.
		The .view.angle and .view.ifov are used to scale the zoom parameter to simulate 
		the expected view.	

		#####  INPUT  ####
		time	=  float value of time elapsed since the start of fly-by. 

		step	=  float value in seconds of time elapsed since previous 
				   frame. Time steps can be stochastic.

		first	=  Boolean flag that when set to true indicates the first
				   frame in the simulation for special handling.
				   Default = False
			
		####  OUTPUT  ####
		.angle	=  Float value of the delta angle represented in the current 
				   axes instance. Overwrites previous angle.

		"""
		
		time_el = step
		self.setView(time)
		ifov = self.view.ifov
		ang = self.view.angle
		scx = ifov/(self.pixSCX*2)
		if first is True:
			z2 = self.startZoom/abs(scx)
			self.axes.camera.zoom = z2
			self.setCamera(az=ang)
		else:
			old_ang = self.angle
			if len(self.axes._wobjects[1]._transformations) > 1:
				self.axes._wobjects[1]._transformations.remove(self.axes._wobjects[1]._transformations[0])
			if len(self.axes._wobjects[2]._transformations) > 0:
				self.axes._wobjects[2]._transformations.remove(self.axes._wobjects[2]._transformations[0])
			CI_dist = time_el*(self.interceptorVelocity-self.cometVelocity)
			com_dist = time_el*self.cometVelocity
			z2 = self.startZoom/abs(scx)
			self.axes.camera.Reset()
			self.axes.SetView({'zoom':z2})
			self.axes._wobjects[1]._transformations.insert(0, vv.Transform_Translate(dy = com_dist))
			self.axes._wobjects[2]._transformations.insert(0, vv.Transform_Translate(dy = com_dist))
			self.setCamera(az = ang, el=0)
			loc = np.asarray(self.axes.GetView().get('loc'))
			loc[1] = CI_dist	
			self.axes.SetView({'loc': tuple(loc)})
		self.axes.Draw()
		self.figure.DrawNow()
		self.angle = self.axes.camera.azimuth
		return
	
class sim(window):
	"""
	Simulation interface, derivative of .window() class.	
	"""
	def __init__(self, cometDict, encounterV, temporalR, Opath, closestA=1000, saveMemory=True):
		"""		   
		#####  INPUT  ####
		cometDict	: dict
			Available comet shape models. Needed to identify the parameters and locate 
			the shape models to create the comet mesh for the simulation. 
			Stored in class property .cdict
			   
		encounterV	: float or int
			Flyby relative velocity, defined in ESA-COMET-SYS-RS-001 as the modulus of 
			the spacecraft heliocentric velocity - comet heliocentric velocity.
						
		closestA	: float or int
			LoS distance in km between the spacecraft and comet at closest approach. 
			Default is 1000km.	 

		temporalR	:	float
			Frame rate in Hz of the TIRI dectector during the encounter. 
			Stored in class property .timeRes
			   
		Opath		: str 
			Path to directory where comet shape models are stored on disk. Simulation 
			output will be saved to this directory. 
			Stored in class property .link
			   
		saveMemory	: bool
			When set to True, runs some RAM saving evasive manoeuvres. Default is True. 
			Can be adjusted with .flag = `bool`
			Stored in class property .flag		
		"""
				
		self.timeline(encounterV, closestA, temporalR)
		self.delta ='\u03B4'
		self.cdict = cometDict
		self.link = Opath
		self.timeRes = temporalR
		self.flag = saveMemory
		if self.flag:
			self.setupDirs()
			
	def setupDirs(self):
		"""
		If .flag is set, output directories for raw and processed simulation frames 
		are created on disk in the location provided by .link.
		
		####  OUTPUT  ####
		.outputPath : @property, str
			Path to parent folder (naming structure `simulationData_Xs', where X is 
			.timeRes rounded to 4 significant figures) containing new output directories.
			
		"""
		
		self.outputPath = os.path.join(self.link, f'simulationData_{np.round(self.timeRes, 4)}s')
		rawPath = os.path.join(self.outputPath, "RAW")
		annPath = os.path.join(self.outputPath, "PROCESSED")
		if not os.path.isdir(self.outputPath):
			os.mkdir(self.outputPath)
			os.mkdir(rawPath)
			os.mkdir(annPath)
		elif os.path.isdir(self.outputPath):
			if not os.path.isdir(rawPath):
				os.mkdir(rawPath)
			if not os.path.isdir(annPath):
				os.mkdir(annPath)
		return

	def timeline(self, eVel, closeA, timeR):
		"""
		Calculates total distance covered by the spacecraft in the direction of travel 
		over the encounter and returns an array of timesteps for the encounter at a given 
		frame rate.
		
		####  INPUT	 ####
		eVel		: float or int 
			The encounter relative velocity.
		
		closeA		: float or int
			Distance between the spacecraft and comet at closest approach.
		
		timeR		: float 
			The frame rate during the encounter.
			
		####  OUTPUT  ####
		.timeSteps	: @property, arr-like
			Array of timesteps for duration of encounter at given resolution.
		
		"""
		
		encounterD = (closeA/np.tan(np.radians(5)))*2
		tTotal = encounterD/eVel
		self.timeSteps = np.round(np.arange(0, tTotal,timeR),4)
		return 
						
	def windowSetup(self, name='c67p'):
		"""
		Instantiates the parent class and sets up the 3D visual window for the simulation.		
		
		####  INPUT	 ####
		name	: str
			String containing the name of the comet to use in the simulation, this must 
			match the name in the comet dictionary. Default is `c67p`
			
		####  OUTPUT  ####
		Interactive window pop-up, set for simulation with models displayed. See .window() 
		class for more information on customising the visuals.	
		
		"""
		super().__init__()
		self.name = name
		self.whichComet(self.cdict)
		self.setup()
		return
				
	def flyBy(self):
		"""
		Runs the simulation, which can be viewed live in the interactive pop-up window.
		
		####  OUTPUT  ####
		.images		: @property, arr-like
			If .flag is set to False, list containing the still frames at each time step.
			Can be adjusted or re-assigned in the event of a crash with .images(`list`).
			
		.angleList	: @property, arr-like
			List of the delta angles at each time step. Can be adjusted or re-assigned 
			in the event of a crash with .angleList(`list`).
			
		.indexList	: @property, arr-like
			List of the indexes associated with timesteps that result in output 
			simulation frames. This is needed in this current iteration of the sim as 
			time periods within the encounter during which TIRI either cannot observe or
			does not observe a datacube are "sped up" to a rate of 1fps, by 
			dropping ever other frame.
			
		.timeList	: @property, arr-like
			List of the timesteps for each output simulation frame.
			
		.ifovList	: @property, arr-like
			List of the pixel scale in km/pix at for each output simulation frame.	
		"""
		
		plt.ion()
		if self.flag:
			outputPath = os.path.join(self.outputPath, "RAW")
		else:
			self.images = []
		self.angleList = []
		self.indexList = []
		self.timeList = []
		self.ifovList = []
		self.distList = []				
		for i,t in enumerate(self.timeSteps):
			validFlag = False
			if t == 0:	
				self.movingView(t, 0, first=True)
				self.draw()
				self.indexList.append(i)
				self.timeList.append(t)
				self.ifovList.append(self.view.ifov)
				self.distList.append(self.view.distance)
				validFlag = True
			else:
				self.view.time = t
				a = self.view.timeView()
				if t%1 == 0 or (a > 82 and a < 100) or (a > 94 and a <125) or \
					(a > 122.99 and a < 128) or (a > 159.99 and a < 161.5) or \
					(a > 172.99 and a < 173.2):
					validFlag = True
					self.indexList.append(i)
					self.timeList.append(t)
					self.movingView(t, t-self.timeSteps[i-1])
					self.ifovList.append(self.view.ifov)
					self.distList.append(self.view.distance)
					self.draw()
			if validFlag:
				temp_image = vv.getframe(vv.gca())
				self.angleList.append(self.view.angle)
				if self.flag == True:
					print(f'Saving t={t}s simulation frame to {outputPath}/{i}.npy')
					np.save(os.path.join(outputPath, f'{i}.npy'), np.asarray(temp_image))
				else:
					self.images.append(temp_image)	
		return
		
	def defineFilts(self):
		"""
		Assigns class attributes containing the locations of filter boundaries in pixels 
		(all pixel locations multiplied by two to account for screen dpi, i.e. the number 
		of pixels in the output simulation images is double the number TIRI actually has) 
		for plotting.

		####  OUTPUT  ####
		.fpa		: @property, child
			Instance of FPA() class.
			
		.fEdges		: @property, dict
			Dictionary containing the pixel location of the boundary edges of the TIRI 
			filters in the TIRI fov.
			
		.ox		: @property, int
			Central pixel x-axis of TIRI fov.
			
		.oy		: @property, int
			Central pixel y-axis of TIRI fov.
			
		"""

		self.fpa = FPA()
		self.fpa.factor=2
		self.fEdges = self.fpa.filterOverlay()
		self.ox, self.oy = np.asarray(self.fpa.xy_dim)/2
		return
									
	def ind(self, lims):
		"""
		Returns index of simulation frames that occur within the angular limits
		of a datacube requirement.
		   
		#####  INPUT  ####
		lims	: Tuple 
			 Delta angles, in the form (start, stop), between which TIRI is expected to
			 actively observe.
						   
		####  OUTPUT  ####
		.lims	: @property, Tuple 
			Indexes, in the form (start, stop), which locate the subset of simulation 
			frames that occur between the given delta angles.
								
		"""
				
		ab = self.angleList
		rnge = [i for i,j in enumerate(ab) if lims[0]<j<lims[1]]
		self.lims = [min(rnge), max(rnge)]
		return 
		
	def rerun(self, num):
		"""
		A roundabout way of getting around having requirements for a variable number of
		datacubes in each observing period (as recorded in observing plan/instrument 
		requirement documents). Calls datacube annotation function directly
			   
		#####  INPUT  ####
		num		: int
			Number of datacubes required.
		 
		"""
		
		for i in np.arange(num):
			if i == 1:
				self.lims[0] = self.Done[-1]+1
			self.dataCube(self.lims)
		return
		
	def rotateTime(self, filterEdge, angle):
		"""
		Takes in an array of coordinates defining the pixel locations of the TIRI filter 
		boundaries edges and returns coordinates rotated to compensate for pointing
		mirror rotation at the given time.

		#####  INPUT  ####
		filterEdge		:  Arr-like 
			(x,y) coordinates of the filter boundaries.
		
		angle			: float
			Delta angle.

		####  OUTPUT  ####
		rFilter		: Arr-like
			(x,y) coordinates of the filter boundaries, rotated about the centre of the 
			FOV by an angle determined by the delta angle at the timestep.		
		"""
		
		angle = 90 - angle
		rFilter = rotateFilters(filterEdge, angle, self.ox, self.oy)
		return np.asarray(rFilter)
		
	@staticmethod
	def checkIntersects(l1, a1):
		"""
		Checks for an intersection between a straight line and an arc.
		
		####  OUTPUT  ####
		bool set to True for intersection, False for no-intersection
		"""
		
		return l1.intersects(a1)
				
	def scanMove(self, angle, whichFilt):
		"""
		Calculates the x- and y-axis translation required to centre a given filter on, 
		or as close as possible, the comet nucleus.
		The apparent translation of the filters is not linear as a result of the 
		rotation of the TIRI pointing mirror, which introduces an apparent trajectory
		along an arc of ~ 9 deg.

		#####  INPUT  ####
		angle		: float
			Delta angle.
		
		whichFilt	: str
			Letter identifier for filter, e.g. 'A'.

		####  OUTPUT  ####
		Mxy			: arr-like 
			Magnitude of (x, y) translation in pixels.

		.vhflag		: @property, str  
			Updates the vhflag attribute with a 'v' if the trajectory arc intersects with
			a vertical filter, and 'h' if the trajectory arc doesn't intersect with
			a vertical filter and does intersect with a horizontal filter. 
				
		"""

		if whichFilt == 'broadband':			
			return [0,0], [0,0]
		validE = ['vL', 'vR', 'hL', 'hR']
		Edges = [self.fEdges.get(whichFilt).get(i) for i in validE] 
		options = {}	
		_ = [options.update({i: self.rotateTime(Edges[j], angle)}) for j,i in enumerate(validE)]
		rArc = np.asarray(arc(self.ox, self.oy))
			
		a1 = LineString(rArc.T)
		isct = {}
		fa = None
		for edge, arr in options.items():
			l1 = LineString(np.asarray(arr).T)
			if self.checkIntersects(l1, a1):
				inters = l1.intersection(a1)
				isct.update({edge: inters})
				if fa == edge[0]:
					self.vhflag = fa
					break
				else:
					fa = edge[0]					
		if len(isct) == 0:
			return False
		halfx = (options.get(f'{fa}R')[0][0]-options.get(f'{fa}L')[0][0])/2
		if fa:
			if f'{fa}L' in isct:
				slitcorr = self.ox-halfx
				Mxy = [isct.get(f'{fa}L').x-slitcorr, isct.get(f'{fa}L').y-self.oy]
			else:
				slitcorr = self.ox-(halfx*2)
				Mxy = [isct.get(f'{fa}R').x-slitcorr, isct.get(f'{fa}R').y-self.oy]
		return Mxy

	def createPoly(self, mx, my, whichFilt, angle, flag):
		"""
		Returns a polygon object with the spatial parameters (edge, rotation, location, 
		etc) of a specified filter at a given time. Takes into account rotation and (x,y) 
		translation in spacecraft rest frame as a result of the movement of the TIRI 
		pointing mirror.

		####  INPUT	 ####
		mx			: float 
			x-axis translation of filter in spacecraft rest frame.	
			   
		my			: float 
			y-axis translation of filter in spacecraft rest frame.
			
		whichFilt	: str
			Letter identifier for filter, e.g. 'A'.
		
		angle		:
			Delta angle.
		
		flag		: str
			Letter identifier for filter subset, either vertical or horizontal, e.g. 'h'.

		####  OUTPUT  ####
		polygon		: obj 
			An interior set of (x,y) coordinates consisting of the infinitely many points 
			within a boundary set of (x,y) coordinates consisting of one or more Curves, 
			and an exterior set of of (x,y) coordinates with boundaries defined by the 
			left edges of filter i and filter i+1, with any rotation and (x,y) scan 
			translation applied.		
		"""
			
		if flag == 'h':
			LEdge = self.fEdges.get(whichFilt).get('hL')
			REdge = self.fEdges.get(whichFilt).get('hR')
		else:
			LEdge = self.fEdges.get(whichFilt).get('vL')
			REdge = self.fEdges.get(whichFilt).get('vR')
		rfiltL = self.rotateTime(LEdge, angle)
		rfiltR = self.rotateTime(REdge, angle)
		lx,ly = np.asarray(rfiltL[0])-mx, np.asarray(rfiltL[1])-my
		rx, ry = np.asarray(rfiltR[0])-mx, np.asarray(rfiltR[1])-my
		polygon = Polygon([(lx[0], ly[0]), (lx[1], ly[1]),(rx[::-1][0], ry[::-1][0]),\
		(rx[::-1][1], ry[::-1][1])])
		return polygon
					
	def dataCube(self, lims):
		"""
		Takes index values of a subset of simulation frames within a datacube 
		requirement and annotates appropriately. Considers integration times per filter 
		and trajectory of comet/filters.
				   
		#####  INPUT  ####
		lims	: tuple 
			(start, stop) index of the simulation frames that occur between the delta 
			angle limits of a datacube requirement.
			
		####  OUTPUT  ####
		.SaveImages		: @property, arr-like
			Appends the annotated simulation frame to a list of stored images
			
		.Done			:@property, arr-like
			Appends the index of the annotated simulation frame to the end of the list.
			
		"""
		
		filterOrder = self.fpa.order
		it = {'A': 0.0552, 'B': 0.0448, 'C': 0.0400, 'D': 0.0435, 'E': 0.0520, 'F': 0.0743, 'G': 0.1452, 'H': 0.0710, 'broadband': 0.0002}
		print(f'data cube at {self.delta} = {np.round(self.angleList[self.lims[0]],1)}')
		a = 0
		b = 0
		for e in range(lims[0],lims[1]):
			if a > 9:
				self.vhflag = 'v'
				ax.clear()
				break
			if self.flag:					
				inputPath = os.path.join(self.outputPath, "RAW")
				im = np.load(os.path.join(inputPath, f'{self.indexList[e]}.npy'), mmap_mode=None)
			else:
				im = self.images[e]
			whichFilt = filterOrder[a]
			angle = self.angleList[e]
			mxy = self.scanMove(angle, whichFilt)
			if not mxy:
				a+=1
				b=0
				continue		
			ax.imshow(im)		
			ax.set_xlim(0,1280)
			ax.set_ylim(0,960)
			ax.text(1000,40, f't = {np.round(self.timeList[e])}s\n{self.delta} = {np.round(self.angleList[e],1)}\n{np.round(self.ifovList[e], 2)}km/pix', color='k', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
			p = self.createPoly(mxy[0], mxy[1], whichFilt, angle, self.vhflag)
			ax.fill(*p.exterior.xy, color='pink', alpha=0.3)	
			ax.set_xticks([])
			ax.set_yticks([])
			for i in filterOrder:
				if i == 'broadband':
					p = self.createPoly(mxy[0], mxy[1], i, angle,  'h')
					ax.plot(*p.exterior.xy)
				else:
					for j in ['v', 'h']:
						p = self.createPoly(mxy[0], mxy[1], i, angle, j)
						ax.plot(*p.exterior.xy)		
			plt.show()
			fig.tight_layout(pad=0)
			fig.canvas.draw()				
			data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
			data = data.reshape(tuple(np.asarray(fig.canvas.get_width_height()[::-1])*2) \
			+ (3,))
			if self.flag:
				outputPath = os.path.join(self.outputPath, "PROCESSED")
				np.save(os.path.join(outputPath, f'{self.indexList[e]}A.npy'), np.asarray(data))
			else:
				self.saveImages[e] = data
			ax.clear()
			b+= self.timeRes/2
			if b >= it.get(whichFilt):
				a+=1
				b = 0
			self.Done.append(e)
		return
	
	def notDatacube(self, e):
		"""
		Annotates simulation frames that are not part of datacubes, and are "sped up" to 
		a rate of 1fps.
		
		####  INPUT	 ####
		e		: int
			Index of simulaton frame.
			
		####  OUTPUT  ####
		.SaveImages		: @property, arr-like
			Appends the annotated simulation frame to a list of stored images.
	
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

		path = Path(verts, codes)
		filterOrder = self.fpa.order
		if self.flag:
			inputPath = os.path.join(self.outputPath, "RAW")
			im = np.load(os.path.join(inputPath, f'{self.indexList[e]}.npy'), mmap_mode=None)
		else:
			im = self.images[e]
		ax.imshow(im)
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_xlim(0,1280)
		ax.set_ylim(0,960)
		ax.text(1000,40, f't = {np.round(self.timeList[e])}s\n{self.delta} = {np.round(self.angleList[e],1)}\n{np.round(self.ifovList[e], 2)}km/pix', color='k', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))		
		if e == 0 or (self.timeList[e] - self.timeList[e-1]) == 1:
			patch = patches.PathPatch(path, facecolor='yellow', lw=1.5)
			ax.add_patch(patch)
			ax.text(170,87.5, f'x{int(1/self.timeRes)}', color='yellow', fontsize='x-large', family='fantasy', va='center', ha='left', variant='small-caps', path_effects=[pe.withStroke(linewidth=3, foreground="k")])
		if self.angleList[e] < 80:
			plt.show()
			fig.tight_layout(pad=0)
			fig.canvas.draw()
			data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
			data = data.reshape(tuple(np.asarray(fig.canvas.get_width_height()[::-1])*2) \
			+ (3,))
		else:
			for i in filterOrder:
				if i == 'broadband':
					p = self.createPoly(0,0, i, self.angleList[e],	'h')
					ax.plot(*p.exterior.xy)
				else:
					for j in ['v', 'h']:
						p = self.createPoly(0,0, i, self.angleList[e], j)
						ax.plot(*p.exterior.xy)
			plt.show()
			fig.tight_layout(pad=0)
			fig.canvas.draw()
			data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
			data = data.reshape(tuple(np.asarray(fig.canvas.get_width_height()[::-1])*2) \
			+ (3,))
		if self.flag:
			outputPath = os.path.join(self.outputPath, "PROCESSED")
			np.save(os.path.join(outputPath, f'{self.indexList[e]}A.npy'), np.asarray(data))
		else:
			self.saveImages[e] = data
		ax.clear()
		return	
		
	def outputTable(self):
		"""
		Creates table with entries of time since start of encounter, distance between 
		comet interceptor and target comet along the line of sight, delta angle and 
		pixel scale, alongside their respective units, for each simulation frame.
		
		####  OUTPUT  ####
		datatab		: tab, ascii format
			Table of dimensions 4xN where N is the total number of frames in the 
			simulation. Columns are [Time, Delta Angle, Distance, Pixel Scale] with units
			[s, deg, km, km/pixel] respectively.
			
		"""
		
		datatab = Table()					
		datatab['Time'] = self.timeList
		datatab['Time'].unit = u.s
		datatab['Delta Angle'] = self.angleList
		datatab['Delta Angle'].unit = u.deg
		datatab['Distance'] = self.distList
		datatab['Distance'].unit = u.km
		datatab['Pixel Scale'] = self.ifovList
		datatab['Pixel Scale'].unit = u.km/u.pix
		return datatab
	
	def Output(self, save=True, fRate=40):
		"""
		Processes simulation frames, adds annotations and representative movements.
		
		#####  INPUT  ####
		save	: bool
			When set to true, saves the annotated images as a gif in the output directory.
			Default is True.
			
		fRate	: int
		   Number of frames per second in output gif. Default is 40.
			
		####  OUTPUT  ####
		If .flag is True, writes annotated simulation frames to file.
		
		If save is True, writes gif video of annotated simulation frames and csv
		containing key parameters for each simulation frame (see .outputTable() for more 
		details) to file.
			
		"""
		if not self.flag:
			self.saveImages = [0]* len(self.indexList)
		global ax
		global fig
		print('Defining filter positions')
		self.defineFilts()
		self.vhflag = 'v'
		self.Done = [0] 
			
		print('Adding annotations at angles where there is a requirement to observe')
		fig = plt.figure(num=1, clear=True)
		ax = fig.add_subplot()
		# requirement: 1 datacube and 3 thermal images at 90 +/- 10 deg
		anglz = np.asarray(self.angleList)
		self.ind((anglz[np.squeeze(np.where(anglz >= 85.))[0]],anglz[np.squeeze \
		(np.where(anglz >= 100.))[0]]))
		self.rerun(1)
		# requirement, 1 datacube and 3 thermal images at 110 +/- 15 deg
		self.ind((anglz[np.squeeze(np.where(anglz >= 95.))[0]],anglz[np.squeeze \
		(np.where(anglz >= 125.))[0]]))
		self.rerun(1)
		# requirement, 1 datacube at 128 +/- 5 deg
		self.ind((anglz[np.squeeze(np.where(anglz >= 123.))[0]],anglz[np.squeeze \
		(np.where(anglz >= 128))[0]]))
		self.rerun(1)
		# requirement, 1 datacube at 160 +/- 5 deg
		self.ind((anglz[np.squeeze(np.where(anglz >=160.))[0]],anglz[np.squeeze \
		(np.where(anglz >= 162.))[0]]))
		self.rerun(1)
		# requirement, 2 datacubes at > 173		
		self.ind((anglz[np.squeeze(np.where(anglz >= 173.))[0]],anglz[np.squeeze \
		(np.where(anglz >= 173.5))[0]]))
		self.rerun(2)			
		print('Adding annotations at angles where there is no requirement to observe:')
		leftover = [i for i in np.arange(0,len(self.indexList)) if i not in self.Done]
		leftover = [0] + leftover
		for i in leftover:
			self.notDatacube(i) 
		if save:
			saveName = f'TIRIencounter_{np.round(self.timeRes, 4)}'
			print('Assembling Table')
			arr = self.outputTable()
			print('Assembling Gif')
			if self.flag:
				inputPath = os.path.join(self.outputPath, "PROCESSED")
				aIms = []
				files = sorted([i for i in os.listdir(inputPath) if i.endswith('.npy')], key=lambda \
						x: float(x.strip('A.npy'))) 
				for im in files:
					aIms.append(np.load(os.path.join(inputPath, im), mmap_mode=None))
			else:
				aIms = self.saveImages
			print(f'Saving outputs as {saveName}.gif and {saveName}.csv')
			ascii.write(arr, f'{os.path.join(self.outputPath,saveName)}.csv', format='csv')
			imageio.mimsave(os.path.join(self.outputPath, f'{saveName}.gif'), aIms, duration=fRate)
		plt.close(fig)
		plt.close('all')
		return

class FPA:
	""" 
	Routines to simulate the projection of the TIRI filter array on the detector. 
	"""
	
	def __init__(self):
		self.loadLayoutI()
		self.factor = 1
		
	def loadLayoutI(self):
		self.nbN = 8
		self.bbN = 1
		self.nb_xy = [32,480]
		self.bb_xy = [200,140]
		self.gap_nb = 8
		self.gap_bb = 57
		self.m_h = 22
		self.m_v = 20
		self.order = ['broadband', 'A', 'B', 'C', 'D', 'broadband', 'E', 'F', 'G', 'H']
		return
		
	def setScale(self, n):
		self.factor = n
		return

	def getFPAdims(self):
		xy_dim = [640,480]
		xy_dim = [xy_dim[0]*self.factor, xy_dim[1]*self.factor]
		return xy_dim
	
	def verticalCentres(self):
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


"""
#########################################################################################
									USAGE EXAMPLE
#########################################################################################
"""
plt.ioff()

# Make nested dictionary of available comet models by calling function
# makecDict([names], [apehelion], [perihelion], [file names]):
cometDict = makecDict(['c67p', 'halley', 'tempel', 'c81p'], [5.68, 35.14, 4.748, \
5.308],[1.243, 0.59278, 1.542, 1.592],['ROS_CG_M001_OSPCLPS_N_V1.OBJ', \
'Halley_Giotto_Vega_Stooke_Model_1.obj', 'TEMPEL1_9P_K032_THO_V01.OBJ', \
'wild2_cart_full.tab.obj'])

# define path to the directory with comet shape models and where any output will be saved
outLink = '/path/to/dir' 


# Instantiate sim class object:
threeD = sim(cometDict, 70, 1/60, outLink)

# Simulation runs by default in save memory mode (i.e. simulation frames are saved to disk 
# instead of stored in memory, and loaded into memory one-by-one when needed). To turn 
# this off run:
# threeD.flag = False

# Instantiate window and display the 3D model:
threeD.windowSetup('c67p')

# Run closest approach fly-by "simulation" which produces an array of images and an 
# array of delta angles. 
threeD.flyBy()


# Save flyby data just in case of crash, so you don't have to rerun the simulation:
spareAngles = threeD.angleList
spareIndexes = threeD.indexList
spareTimes = threeD.timeList
spareIfov = threeD.ifovList
spareDist = threeD.distList

# In the event of having to instantiate a new gif and window object after simulation
# these lists can be used to skip the need to start from scratch by setting:
threeD.angleList = spareAngles
threeD.indexList = spareIndexes
threeD.timeList = spareTimes
threeD.ifovList = spareIfov
threeD.distList = spareDist

# Annotate the simulation frames and create gif and table output:
threeD.Output(fRate=50)



