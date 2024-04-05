#!/usr/bin/env python3

"""
#########################################################################################
									USAGE INFO
#########################################################################################
Comet download and display written using Thomas Albin's Space Science Tutorials - Part 12:  
https://github.com/ThomasAlbin/SpaceScienceTutorial

Code explicitly written to be run interactively, ipython in particular but should work in 
a jupyter notebook too. Uses a deprecated visual engine (visvis) that crashes all the 
time and takes up a lot of RAM, too lazy to make it more efficient or learn how to 
implement maintained versions of comparable packages. Fight me.

Comet 67P/Churyumov–Gerasimenko is the default
Comet velocity default = 34km/s
Flyby velocity vector default = 70km/s

Comet meshes have to be found and downloaded prior to using. Most available comet shape
models can be found here: https://sbn.psi.edu/pds/shape-models/
Example download of c67P mesh:

directory_path = 'directory to wherever you want to keep the file'
download_url = 'https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/rosetta
				/raw/misc/models/ROS_CG_M001_OSPCLPS_N_V1.OBJ'
download_file(directory_path, download_url)

#########################################################################################
										CLASSES
#########################################################################################
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	view() - 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Calculates the location of the spacecraft w.r.t. the comet at a certain point in time, 
	for use to set the view of the comet in the simulation 

	##### CALL ####
	X = view(time, cometV, interceptorV)
	
------------------------------------------------------------------------------------------				   		
	Class attributes:
------------------------------------------------------------------------------------------				   
	#####  INPUTS ####
	time  			=  Float value in seconds of time elapsed since the start of fly-by.  
				   
	cometV			=  Float or integer of the desired velocity the simulated comet is 
					   travelling at. 
	
	interceptorV	=  Float or integer of the desired flyby relative velocity, defined 
					   in ESA-COMET-SYS-RS-001 as the modulus of the spacecraft 
					   heliocentric velocity - comet heliocentric velocity. For 
					   simplicity, the spacecraft velocity will be taken as 
					   interceptorV-cometV

------------------------------------------------------------------------------------------				   		
	Class properties:
------------------------------------------------------------------------------------------	
	.angle			= Delta angle in degrees returned by time_view()
	
	.distance		= Seperation in km returned by view_distance().

	.ifov			= ifov in km/pixel returned by view_ifov().
		   
------------------------------------------------------------------------------------------				           
	Class methods:	
------------------------------------------------------------------------------------------				   	  		  
	timeView		=  Calculates the current delta angle between the spacecraft and the 
					   comet. Does this over four calculations:
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

					   ##### CALL ####
					   X.timeView()
				   
					   ##### INTERNAL CALLS #####
					   None
				   
					   #####  INPUTS ####
					   None
								
					   #### OUTPUTS ####
					   ang = Delta angle in degrees 

------------------------------------------------------------------------------------------				   	  		  
	viewDistance	=  Calculates the current distance between the spacecraft and the 
				   	   comet. Where distance is defined as the hypotenuse of the right-
				   	   angled triangle with shortest side (opposite) equal to closest 
				   	   approach (1000km) and angle is equal to the delta angle returned by
				   	   timeview()

					   ##### EQUATION ####
					   sin(theta) = opposite/hypotenuse
					   distance = 1000/sin(delta)
				   				
				   				   
					   ##### CALL ####
				   	   X.viewDistance()
				   
					   ##### INTERNAL CALLS #####
				   	   None
				   
					   #####  INPUTS ####
					   None
								
					   #### OUTPUTS ####
					   The distance in km between the s/c and comet
			   				
------------------------------------------------------------------------------------------
	viewPixscale	=  Calculates the ifov, i.e. the number of km per pixel of the TIRI 
					   detector given a viewing angle returned by timeView()/2 and a 
					   seperation returned by viewDistance().
				   
					   ##### EQUATION ####
					   tan(theta) = opposite/adjacent
					   ifov = 2*(tan(delta/2)*seperation
				   
					   ##### CALL ####
					   X.viewPixscale()
				   
					   ##### INTERNAL CALLS #####
					   .viewDistance()
				   
					   #####  INPUTS ####
					   None
								
					   #### OUTPUTS ####
					   The ifov of TIRI in km/pixel 
		   				
------------------------------------------------------------------------------------------				   
------------------------------------------------------------------------------------------	


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	cometMesh() - 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Creates a 3D mesh from a shape model and applies transformations to get the best 
	viewing angle (for no purpose but ~aesthetics~). Mostly just a redundant wrapper for
	visvis classes and methods.

	##### CALL ####
	X = cometMesh(faces, vertices)
	
------------------------------------------------------------------------------------------				   		
	Class attributes:
------------------------------------------------------------------------------------------				   
	#####  INPUTS ####
	faces		=  Array containing the location of faces in the shape model
				   Can be adjusted with X.faces = array
	
	vertices	=  Array containing the location of vertices in the shape model
				   Can be adjusted with X.vertices = array
							   		   	   				   
------------------------------------------------------------------------------------------				   
------------------------------------------------------------------------------------------


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	coma() - 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Creates a very inauthentic representation of a coma. Just for visuals, not 
	scientifically relevant.

	##### CALL ####
	X = coma()

------------------------------------------------------------------------------------------				   		
	Class properties:
------------------------------------------------------------------------------------------	
	.cloud		=  Point cloud representing the coma returned by generateComa()				  
			   
------------------------------------------------------------------------------------------				           
	Class methods:	
------------------------------------------------------------------------------------------				   	  		  
	generateComa	=  Creates a 29000 point cloud composed of four different parts (for
						no purpose but ~aesthetics~):
					   1) a 3d cone of 6000 points with coordinates in the range 
					      (-38 < x < 35)km, (-38 < y < 880)km, (-41 < z < 37)km.
					   2) a 3D normal distribution of 10000 points with coordinates in 
					      the range (-1900 < x < 2150)km, (-3100 < y < 5200)km, 
					      (-1900 < z < 2000)km.
					   3) a 3D cone of 3000 points with coordinates in the range 
					      (-8 < x < 8)km, (-5 < y < 260)km, (8 < z < 8)km.
					   4) a 3D normal distribution of 10000 points with coordinates in 
					      the range (-200 < x < 200)km, (-3000 < y < 5300)km, 
					      (-200 < z < 200)km.
				   
					   ##### CALL ####
				   	   X.generateComa()
				   
					   ##### INTERNAL CALLS #####
				   	   None
				   
					   #####  INPUTS ####
					   None
								
					   #### OUTPUTS ####
					   pp		=  Point set with 29000 points
     				   		   	   				   
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------	       


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	window() -  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	class wrapper for all things simulation.

	##### CALL ####
	X = window()
	
------------------------------------------------------------------------------------------				   		
	Class attributes:
------------------------------------------------------------------------------------------				   
	#####  INPUTS #####
	None
				  
	##### DEFAULTS #####
	.name				 =  String containing the name of the comet to use in the 
						    simulation, this must match the name in the comet dictionary.
				            Default is 'c67p'.
				            Can be adjusted with X.name = string	   

	.cometVelocity		 =  Float or integer of the desired velocity the simulated comet 
						    will travel at. Informs the view of the comet for each 
						    timestep.
					        Default is 15km/s
					        Can be adjusted with X.cometVelocity = float
					        
	.interceptorVelocity =  Float or integer of the desired flyby relative velocity, defined 
					   		in ESA-COMET-SYS-RS-001 as the modulus of the spacecraft 
					   		heliocentric velocity - comet heliocentric velocity. For 
					   		simplicity, the spacecraft velocity will be taken as 
					   		interceptorVelocity minus cometVelocity for all calculations. 
					   		Informs the view of the comet for each timestep.
						    Default is 70km/s
				   			Can be adjusted with X.interceptorVelocity = float
------------------------------------------------------------------------------------------				           
	Class methods:	
------------------------------------------------------------------------------------------				   	  		  			   	  		  
	whichComet	=  Assigns class attributes based on a search of the given comet 
				   dictionary using .name attribute as the search identifier.
				   
				   ##### CALL ####
				   X.whichComet(cometDict)
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   cometDict		=  Dictionary of available comet shape models. Needed 
				   					   to identify the parameters and locate the shape 
				   					   models to create the comet mesh for the simulation.
				   				
				   #### OUTPUTS ####
				   .cometLink		=  String containing the file name of the comet shape 
				   				   	   model for use in the simulation.
				   .cometVelocity	=  Float or integer of the desired velocity the 
				   					   simulated comet will travel at. Informs the view 
				   					   of the comet for each timestep. 

------------------------------------------------------------------------------------------				   	  		  
	setup		=  Sets up the 3D visual window for the simulation. Calls several 
				   internal methods that define and control the appearance of the window, 
				   comet, coma, camera FOV and lighting. Assigns class attributes in 
				   preparation for running the simulation.
				      
				   ##### CALL ####
				   X.setup()
				   
				   ##### INTERNAL CALLS #####
				   .setWindow(), .loadComet(), .setLight(), .setAxes(), .setCamera(), 
				   .draw(), .getPS()
				   
				   #####  INPUTS ####
				   None 	
							
				   #### OUTPUTS ####
				   .startZoom	=   Stores a baseline scaling for simulated spaceraft to
				   					comet distance by querying the default camera view 
				   					set in local coordinates by the visvis .axes.camera 
				   					class. This is used in the simulation to determine the
				   					needed wcs adjustment of the viewing parameters (i.e. 
				   					camera, location, zoom, angle) for each timestep to 
				   					simulate the flyby.
				   					
				   .startView	=	Stores a baseline location of the simulated spaceraft
				   					w.r.t. the comet by querying the default camera 
				   					location set by the visvis .axes._view class. This is
				   					used in the simulation to determine the needed wcs 
				   					adjustment of the viewing parameters (i.e. camera, 
				   					location, zoom, angle) for each timestep to simulate 
				   					the flyby.
				   					
				   .pixSCX		=	Stores the default pixel scale of the window in local 
				   					coordinates. This is used in the simulation to 
				   					determine the needed wcs adjustment of the viewing 
				   					parameters (i.e. camera, location, zoom, angle) for 
				   					each timestep to simulate the flyby.
		   		   				
------------------------------------------------------------------------------------------	
	setWindow	=  Instantiates the window the simulation will run in, acts as a wrapper
				   for visvis functions and classes. Sets sub-classes as attributes so
				   they can be called for use easily later.
				   
				   ##### CALL ####
				   X.setWindow()
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   None
				   				
				   #### OUTPUTS ####
				   .axes	=   Stores the visvis .axes class.
				   .figure	=	Stores the visvis window wrapper.			
------------------------------------------------------------------------------------------	
	loadComet	=  Loads comet shape model from .link and .cometLink, converts model
				   faces and vertices into a mesh object, stored as a class attribute. 
				   Instantiates a coma() class object and stores the resulting point 
				   set object as a class attribute. Displays comet mesh and point cloud 
				   coma in the current axes instance.
				   
				   ##### CALL ####
				   X.loadComet()
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   None
				   				
				   #### OUTPUTS ####  
				   .coma	=   random point cloud representation of a coma. Not accurate
				   				just for visuals. 
				   .comet	=   3D comet mesh object. 
				   				
------------------------------------------------------------------------------------------	
	setLight	=  Interfaces with the visvis .axes.light class to turn off default 
				   window view lighting, creates a new instance of the light class and 
				   stores as a class attribute. Sets the light to be directional and as 
				   close as i could get it to coming from behind the "spacecraft" to 
				   simulate the direction of light coming from the sun.
				   
				   ##### CALL ####
				   X.setLight()
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   None
				   				
				   #### OUTPUTS ####
				   .lightObj	=  Stores an instance of axes.lights to control the
				   				   lighting in the simulation.			   				
------------------------------------------------------------------------------------------	
	setAxes		=  Adjusts .axes parameters to make the simulation look like is is 
				   happening in space, i.e. setting background colour to black, hiding 
				   the axes, turning off automatic aspect adjustments and setting the 
				   axes camera to 3D to allow x,y and z translations of the simulation.
				   
				   ##### CALL ####
				   X.setAxes()
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   None
				   				
				   #### OUTPUTS ####
				   None
				   				
------------------------------------------------------------------------------------------	
	setCamera	=  Adjusts the simulation "camera", i.e. where the spacecraft appears to 
				   be in location w.r.t. the comet.
				   
				   ##### CALL ####
				   X.setCamera(f=0, az=90, el=0)
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   f	=   Integer value of the projection of the axes, can be set to 
				   			adjust the apparent focal length of the simulation, i.e. the
				   			depth of field.
				   			Default is 0
				   	
				   az	=   The azimuth angle of the "spacecraft". Can be set to adjust the
				   			face-on viewing angle of the comet.
				   			Default is 90
				   			
				   el	=	The elevation of the "spacecraft". Can be set to adjust the
				   			face-on viewing angle of the comet.
				   			Default is 0
				   					
				   #### OUTPUTS ####
				   None
				   				
------------------------------------------------------------------------------------------	
	getPS		=  Determines the pixel scale of the current axes instance in the WCS from
				   pixel coordinates given in the local coordinate system. Can be used to
				   either return the pixel scale in x or y (should be the same) by varying
				   either the x-coordinate inputs by +/- 1 and keeping the y-coordinate 
				   inputs the same or keeping the x-coordinates the same and varying the 
				   y-coordinate inputs by +/- 1.
				   
				   ##### CALL ####
				   X.getPS(x,y,x2,y2)
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   x	=   x-coordinate of a pixel.
				   
				   y	=   y-coordinate of a pixel.
				   
				   x2	=   x-coordinate of neighbouring pixel.
				   
				   y2	=   y-coordinate of neighbouring pixel.
				   				
				   #### OUTPUTS ####
				   The pixel scale of the current axes instance in WCS
				   				
------------------------------------------------------------------------------------------	    
	draw		=  Re-draws the axes and updates everything in it.
				   
				   ##### CALL ####
				   X.draw()
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   None
				   				
				   #### OUTPUTS ####
				   None
				   				
------------------------------------------------------------------------------------------	
	getSV		=  Gets the initial location of the simulated spaceraft w.r.t. the comet 
				   by querying the default camera location set by the visvis .axes._view 
				   class. 
				   
				   ##### CALL ####
				   X.getSV()
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   None
				   				
				   #### OUTPUTS ####
				   startView	=  Initial view (i.e. spacecraft location) of the 
				   				   simulated comet.
				   				
------------------------------------------------------------------------------------------	
	setView		=  Creates a new instance of the view() class and stores as a class 
				   attribute for easy use.
				   
				   ##### CALL ####
				   X.setView(time)
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   time		=  Float value of time elapsed since beginning of fly-by
				   			   in seconds.
				   				
				   #### OUTPUTS ####
				   .view  = instance of view() class
				   				
------------------------------------------------------------------------------------------	
	movingView	=  The method that advances the simulation in time. It takes in the time-
				   step, calls setView() to calculate the ifov and delta angle at the new 
				   timestep. If this is the first frame in the simulation, then 
				   start_zoom() and setCamera() are called to store the baseline view 
				   parameters for the rest of the simulation. 
				   For all other steps in the simulation, the comet and coma translations 
				   from the previous timestep are removed and the axes are reset to the 
				   default view. 
				   To account for non-regular time steps, the time elapsed since the 
				   last simulation frames and the velocity attributes are used to 
				   calculate the distance moved by the spacecraft the comet are set since
				   the last frame. These are set as the location of the axes camera and 
				   the comet translation respectively. The angle and ifov are used to set 
				   and correctly scale the zoom parameter to the appropriate distance 
				   between the spacecraft and the comet
				   
				   ##### CALL ####
				   X.movingView(time, step, first=None)
				   
				   ##### INTERNAL CALLS #####
				   .setView(), .startZoom(), .setCamera()
				   
				   #####  INPUTS ####
				   time		=  float value of time elapsed since the start of fly-by.
				   
				   step		=  float value in seconds of time elapsed since previous 
				   			   frame.
				   
				   first	=  Boolean flag that when set to true indicates the first
				   			   frame in the simulation for special handling.
				   			   Default = False
				   				
				   #### OUTPUTS ####
				   .angle	=  Float value of the delta angle represented in the current 
				   			   axes instance. Overwrites previous angle.
				   				
------------------------------------------------------------------------------------------	
------------------------------------------------------------------------------------------	


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	sim() - 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	the main class to interface with for the simulation, instantiated with a dictionary 
	of available comet shape models.
	
	##### CALL #####
	X = sim(dictionary, saveMemory)
	
------------------------------------------------------------------------------------------				   		
	Class attributes:
------------------------------------------------------------------------------------------				   
	#####  INPUTS ####
	cdict       =  Dictionary of available comet shape models. Needed to identify the 
				   parameters and locate the shape models to create the comet mesh for
				   the simulation.
				   Can be adjusted with X.cdict = dictionary
				   
	saveMemory	=  Boolean flag that when set to True runs some RAM saving evasive 
				   manoeuvres.
				   Default is True.
				   Can be adjusted with X.flag = Bool
				  
	##### DEFAULTS #####
	.timeSteps  =  Numpy array of time steps for the simulation in seconds, rounded to 
				   4 sig fig.
				   Default value is 0 to 350 seconds in steps of 0.0625s.
				   Can be adjusted with X.timeSteps = array
	
				  
	.saveImages =  Array to hold all the annotated simulation images prior to 
				   saving them as a gif.
				   Default is list of 0s the same length as timeSteps. 
				   Can be adjusted with X.saveImages = list
				  
	.saveName   =  String containing the desired name of the gif output.
				   Default is "TIRI_closest_approach".
				   Can be adjusted with X.saveName = string
				   
	.link		=  A string containing the path to directory containing the comet shape 
				   models and where the gif will be saved.
				   Default is the current working directory.
				   Can be adjusted with X.link = string
				   
------------------------------------------------------------------------------------------				           
	Class methods:	
------------------------------------------------------------------------------------------				   	  		  
	windowSetup =  Creates an instantiation of the window class and sets up the 3D visual 
				   window for the simulation.
				   
				   ##### METHOD #####
				   1)

				   ##### EQUATIONS #####
				   None		
				   		   
				   ##### CALL ####
				   X.windowSetup(name)
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   #####  INPUTS ####
				   name		=	String containing the name of the comet to use in the 
				   				simulation, this must match the name in the comet 
				   				dictionary.
				   				Default is 'c67p'.
				   				
				   #### OUTPUTS ####
				   Interactive window pop-up, ready for simulation with 3D comet model and
				   coma representative. See window() class for more information on 
				   customising the visuals.	
				   		   	   				   
------------------------------------------------------------------------------------------				   
	flyBy       =  Runs the simulation, which can be viewed live in the interactive pop-up
	               window and assigns two new class attributes to store the outputs for 
	               later use.
	               
				   ##### METHOD #####
				   1)
				   
				   ##### EQUATIONS #####
				   None
				   				   
				   ##### CALL #####
				   X.flyby()
				   
				   ##### INTERNAL CALLS #####	
				   None
				   	   
				   ##### INPUTS ####
				   None
				   
				   #### OUTPUTS ####
				   If .flag is set to True, saves all simulation frames into a directory
				   called sim/ in .link.
				   .images	  = If .flag is set to False, a list containing the still 
				   				frames at each time step.
				   			    Can be adjusted or re-assigned in the event of a crash
				   			    with X.images(list).
				   .angleList = A list of the delta angles at each time step.
				   				Can be adjusted or re-assigned in the event of a crash 
				   				with X.angleList(list).
				   .distList  = A list of the LoS distance to the comet at each time step.
				   				Can be adjusted or re-assigned in the event of a crash 
				   				with X.distList(list).
				   				
------------------------------------------------------------------------------------------				   				   
	makeOutput	=  Takes the simulation results, annotates them with timestamps and delta
				   angles and stores them by populating X.saveImages. Called 
				   with a flag to turn off/on saving images to the local disk in gif 
				   format.  
				   
				   ##### METHOD #####
				   None
				   
				   ##### EQUATIONS #####
				   None	
				   			   
				   ##### CALL #####
				   X.makeOutput(s=True, f=30)
				   
				   ##### INTERNAL CALLS #####
				   .outputData(), .outputTable()
				   
				   ##### INPUTS ####
				   s	    =   When set to true, saves the annotated images as a gif at 
				   				the directory specified in X.link, with the file name 
				   				specified in X.saveName.
				   				Default is True.
				   f		=   Integer of number of frames per second in output gif.
				   			    Default is 30
				   			    
				   #### OUTPUTS ####
				   If s is set to True, saves a gif with f frames per second 
				   and an array
				   If s is set to True, saves a csv file containing an ascii table 
				   of the physical parameters from the simulation.		   			   			   		        

------------------------------------------------------------------------------------------				   
	outputData	=  Takes in simulation frames and annotates them with time and angle. 
				   
				   ##### METHOD ####
				   None

				   ##### EQUATIONS ####	
				   None			   			   
				   	  
				   ##### CALL ####
				   X.outputData(lims)
				   
				   ##### INTERNAL CALLS #####
				   None
				   
				   ##### INPUTS ####
				   lims			=  Index to locate the simulation frame parameters.
				   
				   #### OUTPUTS ####
				   .SaveImages	=  Updates with the annotated simulation frame.

------------------------------------------------------------------------------------------				   	
------------------------------------------------------------------------------------------				   
	outputTable	=  Takes in simulation timeSteps and outputs an ascii table with columns
				   of the following arrays:
				   
				   1) Time elapsed(s) 
				   2) Delta angle (deg)
				   3) Distance along LoS to comet (km)
				   4) Pixel scale (km)
				   
				   ##### METHOD ####
				   None

				   ##### EQUATIONS ####	
				   None			   			   
				   	  
				   ##### CALL ####
				   X.outputTable()
				   
				   ##### INTERNAL CALLS #####
				   view()
				   
				   ##### INPUTS ####
				   None
				   				   
				   #### OUTPUTS ####
				   dataTab	=  (4 x length(time array)) ascii table

#########################################################################################
										FUNCTIONS
#########################################################################################

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	download_file(DL_PATH, DL_URL) - 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Downloads files from url and saves to directory.
	
	#####  INPUTS ####
	DL_PATH	: str
        	  Download path on the local machine, relative to this function.
    
    DL_URL	: str
              Download url of the requested file.	
             
	##### OUTPUTS #####	   
	None
------------------------------------------------------------------------------------------				   	
------------------------------------------------------------------------------------------				   

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	makecDict(cometNames, cometAphelion, cometPerihelion, cometFiles) - 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Creates a dictionary of comets with available shape models, alongside their velocity
	at perihelion or 1AU, whichever is greatest, and the name of the shape model file.
	
	#####  INPUTS ####
	cometNames 		: list, str
					  Name/identifiers of comets.
    
    cometAphelion	: list, float
    				  Aphelion of comets in AU.
                
    cometPerihelion	: list, float
    				  Perihelion of comets in AU.
    
    cometFiles		: list, str
    				  File names of comet shape models.
    	        
	##### OUTPUTS #####	   
	cometDict		: dict
					  Nested dictionary of form {name: {velocity, file name}}
------------------------------------------------------------------------------------------				   	
------------------------------------------------------------------------------------------				   

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	comVel(aphelion, perihelion, sun_comet_distance) - 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Calculates the velocity of a comet with at a given distance from the sun.
	
	#####  INPUTS ####
	aphelion 			: float
					  	  Aphelion of comet in AU
    
    perihelion			: float
    				  	  Perihelion of comet in AU.
                
    sun_comet_distance	: float
    				  	  Distance of comet from the sun in AU
    	        
	##### OUTPUTS #####	   
	v			: float
	              velocity of comet in km/s at sun_comet_distance        
					                                   			  
------------------------------------------------------------------------------------------				   	
------------------------------------------------------------------------------------------

"""
from astropy.io import ascii
from astropy.table import Table
from astropy import units as u
import imageio
import matplotlib.pyplot as plt
import numpy as	 np
import os
import pandas as pd
import pathlib
from PyQt5.QtWidgets import QWidget, QHBoxLayout
from shapely import Point, LineString, MultiPoint
from shapely.geometry.polygon import Polygon
import time
import urllib.request
import visvis as vv
from visvis import Point, Pointset

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	DEFINE FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

def downloadFile(dl_path, dl_url):
	file_name = dl_url.split('/')[-1]
	pathlib.Path(dl_path).mkdir(parents=True, exist_ok=True)
	if not os.path.isfile(dl_path + file_name):
		urllib.request.urlretrieve(dl_url, dl_path + file_name)
	return
	
def makecDict(cometNames, cometApehelion, cometPerihelion, cometFiles):
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

		
"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	DEFINE CLASSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
	def __init__(self, time, cometVelocity, interceptorVelocity):
		self.time = time
		self.cometV = cometVelocity
		self.iV = interceptorVelocity
		self.angle = self.timeView()
		self.distance = self.viewDistance()
		self.ifov = self.viewPixscale()

	def timeView(self):
		ini_dist = 1000/np.tan(np.radians(5))
		ci_dist = self.time*(self.iV-self.cometV)
		com_dist = self.time*self.cometV
		curr_dist = ini_dist-ci_dist-com_dist
		ang = np.rad2deg(np.arctan(1000/curr_dist))
		if curr_dist < 0:
			ang = (90+(90+ang)) 
		return ang
		
	def viewDistance(self):
		a = self.angle
		return 1000/np.sin(np.radians(a))	
		
	def viewPixscale(self):
		return np.tan((0.26e-3)/2)*self.viewDistance()*2
		
 
class cometMesh:
	def __init__(self, faces, vertices):
		self.faces = faces
		self.vertices = vertices
		self.obj = vv.mesh(vertices=self.vertices, faces=self.faces, verticesPerFace=3)
		self.obj.specular = 0.2
		self.obj.diffuse = 0.9
		self.Rot = vv.Transform_Rotate(angle=100, ax=-1, ay=1, az=-1)
		self.obj.transformations.insert(0,self.Rot)		  

class coma:
	__slots__ = ("cloud")
	def __init__(self):
		self.cloud = self.generateComa()
		
	def generateComa(self):
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
	%gui qt 
	def __init__(self):
		self.name = 'c67p'
		self.cometVelocity = 15
		self.interceptorVelocity = 70
		
	def whichComet(self, cometDict):
		progress(f'Preparing flyby window for {self.name}')
		self.cometLink = cometDict.get(self.name).get('file')
		progress()
		self.cometVelocity = cometDict.get(self.name).get('velocity')
		progress(p=1)
		print(f'Location of shape model: {self.link+self.cometLink}')
		print(f'Setting comet velocity as {self.cometVelocity}')
		return
		
	def setup(self):
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
		progress('Loading comet and coma')
		df = pd.read_csv(self.link+self.cometLink, delim_whitespace=True, names=['TYPE','X1', 'X2', 'X3'], engine='python', header=None, on_bad_lines='warn')									
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
		self.axes.camera.fov = f
		self.axes.camera.azimuth = az
		self.axes.camera.elevation = el
		return
	
	def getPS(self, x, y, x2, y2):
		return np.asarray(self.axes.camera.ScreenToWorld((x,y)))-np.asarray(self.axes.camera.ScreenToWorld((x2, y2)))
		
	def draw(self):
		self.axes.Draw()
		self.figure.DrawNow()
		return
		
	def setView(self, time):
		self.view = view(time, self.cometVelocity, self.interceptorVelocity)
		return		

	def movingView(self, time, step, first=None):
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
	def __init__(self, cometDict, saveMemory=True):
		self.timeSteps = np.round(np.arange(0, 350,0.125),4)
		self.saveName = "comet_encounter"
		self.delta ='\u03B4'
		self.cdict = cometDict
		self.link = os.getcwd()
		self.flag = saveMemory
		self.scanShift = [(0,0)]*len(self.timeSteps)
				
	def windowSetup(self, name='c67p'):
		self.window = window()
		self.window.name = name
		self.window.link = self.link
		self.window.whichComet(self.cdict)
		self.window.setup()
		return	
				
	def flyBy(self):
		plt.ion()
		if self.flag == True:
			print('Memory saving mode is on')
			try:
				os.mkdir(f'{self.link}sim')
			except:
				pass			  
		else:
			print('Memory saving mode is off')
			comet_images = []
		aa = []
		print("Running Simulation")
		for i,t in enumerate(self.timeSteps):
			progress(q=[i,len(self.timeSteps)])
			if t == 0:
				self.window.movingView(t, 0, first=True)
			else:
				self.window.movingView(t, t-self.timeSteps[i-1])
			self.window.draw()
			self.window.axes.Draw()
			self.window.figure.DrawNow()			
			temp_image = vv.getframe(vv.gca())
			if self.flag == True:
				np.save(self.link+'sim/'+str(i), np.asarray(temp_image))
			elif self.flag!=True:
				comet_images.append(temp_image)
			aa.append(self.window.view.angle)
		if self.flag != True:
			self.images = comet_images
		self.angleList = aa
		return		
		
	def makeOutput(self, s=True, f=40):	
		self.saveImages = [0]* len(self.timeSteps)
		global ax
		global fig
		print('Adding annotations:')
		fig = plt.figure(num=1, clear=True)
		ax = fig.add_subplot()
		for i in range(0, len(self.timeSteps)):
			progress(q=[i,len(self.timeSteps)])
			self.outputData(i)
		print('Assembling Table:')
		arr = self.outputTable()	
		if s is True:	
			print(f'Saving outputs as {self.saveName}.gif and {self.saveName}.csv')
			progress()
			imageio.mimsave(f'{self.link+self.saveName}.gif', self.saveImages, fps=f)
			progress()
			ascii.write(arr, f'{threeD.link+threeD.saveName}.csv', format='csv')
			progress(p=1)
		plt.close(fig)
		plt.close('all')
		return
		
	def outputData(self, e):			
		if self.flag == True:
			im = np.load(self.link+'sim/'+str(e)+'.npy', mmap_mode=None)
		else:
			im = self.images[e]
		ax.imshow(im)
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_xlim(0,1280)
		ax.set_ylim(0,960)
		ax.text(1000,40, f't = {np.round(self.timeSteps[e])}s\n{self.delta} = {np.round(self.angleList[e],1)}', color='k', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))				
		plt.show()
		fig.tight_layout(pad=0)
		fig.canvas.draw()		
		data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
		data = data.reshape(tuple(np.asarray(fig.canvas.get_width_height()[::-1])*2) \
		+ (3,))		
		self.saveImages[e] = data
		ax.clear()
		return

	def outputTable(self):
		datatab =Table()
		threeD.distance = [0]*len(threeD.timeSteps)
		threeD.pixscale = [0]*len(threeD.timeSteps)
		threeD.deltaangle = [0]*len(threeD.timeSteps)
		for e,i in enumerate(threeD.timeSteps):
			progress(q=[e,len(self.timeSteps)])
			v = view(i, threeD.window.cometVelocity, threeD.window.interceptorVelocity)
			threeD.distance[e] = v.distance
			threeD.deltaangle[e] = v.angle
			threeD.pixscale[e] = v.ifov										
		datatab['Time'] = threeD.timeSteps
		datatab['Time'].unit = u.s
		datatab['Delta Angle'] = threeD.angleList
		datatab['Delta Angle'].unit = u.deg
		datatab['Distance'] = threeD.distance
		datatab['Distance'].unit = u.km
		datatab['Pixel Scale'] = threeD.pixscale
		datatab['Pixel Scale'].unit = u.km/u.pix	
		return datatab
"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	RUN SIMULATION AND SAVE GIF AND TABLE OUTPUT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
# Make nested dictionary of available comet models by calling function
# makecDict([names], [apehelion], [perihelion], [file names]):

	cometDict = makecDict(['c67p', 'halley', 'tempel', 'c81p'], [5.68, 35.14, 4.748, \
	5.308],[1.243, 0.59278, 1.542, 1.592],['ROS_CG_M001_OSPCLPS_N_V1.OBJ', \
	'Halley_Giotto_Vega_Stooke_Model_1.obj', 'TEMPEL1_9P_K032_THO_V01.OBJ', \
	'wild2_cart_full.tab.obj'])
	
# Instantiate sim class object: 

	threeD = sim(cometDict)	

# set some necessary parameters:
# set the name of the output objects if different from default ['comet_encounter']:

	threeD.saveName = 'TIRI_newrot_sgl' 

# set path to the directory with comet shape models and where gif will be saved,
# must end with "/":

	threeD.link = '/path/to/dir' 

# Set time resolution of simulation if different from default [0 - 350s, step = 0.125s]:
	
	threeD.timeSteps = np.round(np.arange(0, 350,1),4)

# Simulation runs by default in save memory mode (i.e. simulation frames are saved to disk 
# instead of stored in memory, and loaded into memory one-by-one when needed). To turn 
# this off:

	threeD.flag = False
	
# Instantiate window and display the 3D model:

	threeD.windowSetup('c67p')

# Run closest approach fly-by "simulation" which produces an array of images and an 
# array of delta angles:

	threeD.flyBy()	

# Save flyby data just in case of crash, so you don't have to rerun the simulation.
# In the event of having to instantiate a new gif and window object after simulation
# these lists can be used to skip the need to start from scratch by setting
# threeD.images = spareImages and threeD.angleList = spareAngles:

	spareAngles = threeD.angleList

# # Annotate the simulation frames and save as a gif:
	threeD.makeOutput(f=20)












	




	

    























		
			


		
		

	

	
		





    

