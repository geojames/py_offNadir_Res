#!/usr/bin/env python
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
__author__ = 'James T. Dietrich'
__contact__ = 'james.t.dietrich@dartmouth.edu'
__copyright__ = '(c) James Dietrich 2016'
__license__ = 'MIT'
__date__ = 'Fri Jul 29 08:46:22 2016'
__version__ = '1.0'
__status__ = "initial release"
__url__ = "https://github.com/geojames/py_offNadir_Res"

"""
Name:           offNadir_resolution_v1-0.py
Compatibility:  Python 3.5
Description:    This program calculates the pixel resolution and instanteous
                    field of view for a camera from an input file provided by
                    the user. The file nees to be a comma-delimited text file 
                    (*.CSV) with a header row:
                    
                    Name,focal,sensor_x,sensor_y,pixel_x,pixel_y,flyH,angle
                    P3,3.61,6.24,4.71,6000,4000,100,10
                    
                    where, Name is camera name, focal is the focal length of 
                    camera (mm), sensor_x is the sensor's long dimension (mm),
                    sensor_y is the sensor's short dimension (mm), pixel_x is
                    the number of pixels along the long dimension, pixel_y is
                    the number of pixels along the short dimension, flyH is the
                    flying height/altitude of the camera, angle is the 
                    off-nadir angle of the camera (0 = nadir, 90 = horizontal)
                    *the software will only work for low-oblique images (i.e. 
                    no horizon visable) and will spit out a error if the
                    off-nadir angle puts the top part of the verical field of
                    view over the horizon.
                    
                    Run:
                        - Run the software from a Python shell or editor
                        - A file chooser will pop-up (it may try to hide)
                            - Choose your CSV file
                        - The program will run and create a formated text file
                            with the output in the same location as the input
                            file. The name will be the same with the suffix
                            '_resolution'


URL:            https://github.com/geojames/py_offNadir_Res

Requires:       tkinter, numpy, pandas, sympy, matplotlib

Dev ToDo:       N/A

AUTHOR:         James T. Dietrich
ORGANIZATION:   Dartmouth College
Contact:        james.t.dietrich@dartmouth.edu
Copyright:      (c) James Dietrich 2016

"""
#------------------------------------------------------------------------------
# Imports
import os
import sys
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
import pandas as pd
import sympy.geometry as spg
import matplotlib.path as mplPath
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt

#-----------

# footprint

def footprint(sensor):
    
    '''
    Caculates the foot print of the off nadir camera by projecting rays from
    the sensor corners through the "lens" (focal length) out onto the ground.
    It's a lot of fun linear algebra that the SYMPY library handles.
    '''
    
    # Setup DF to house camera footprint polygons
    footprint = pd.DataFrame(np.zeros((1,5)), 
                             columns=['fov_h','fov_v','path','pp_x','pp_y'])
    
    # convert sensor dimensions to meters, divide x/y for corner coord calc
    f = sensor.focal[0] * 0.001
    sx = sensor.sensor_x[0] / 2 * 0.001
    sy = sensor.sensor_y[0] / 2 * 0.001

    # calculate the critical pitch (in degrees) where the horizon will be 
    #   visible with the horizon viable, the ray projections go backward 
    #   and produce erroneous IFOV polygons (90 - 0.5*vert_fov)
    #   exit with error message if critical pitch is exceeded
    
    crit_pitch = 90 - np.rad2deg(np.arctan(sy / f))
    
    if sensor.angle[0] >= crit_pitch:
        print('!!! The provided parameters indicate that the vertical field')
        print('\t of view extends above the horizon. Please start over and')
        print('\t try a shallower camera angle. The maximum angle for this')
        print('\t camera is %0.2f' %(crit_pitch))
        sys.exit()
    
    # calculate horz and vert field of view angles
    footprint.fov_h = 2 * np.rad2deg(np.arctan(sx / f))
    footprint.fov_v = 2 * np.rad2deg(np.arctan(sy / f))
    
    # sensor corners (UR,LR,LL,UL), north-oriented and zero pitch
    corners = np.array([[0+sx,0-f,sensor.flyH[0]+sy],
                       [0+sx,0-f,sensor.flyH[0]-sy],
                       [0-sx,0-f,sensor.flyH[0]-sy],
                       [0-sx,0-f,sensor.flyH[0]+sy]])
    
    # offset corner points by cam x,y,z for rotation
    cam_pt = np.atleast_2d(np.array([0, 0, sensor.flyH[0]]))
    corner_p = corners - cam_pt
    
    # convert off nadir angle to radians
    pitch = np.deg2rad(90.0-sensor.angle[0])
    
    # setup pitch rotation matrix (r_x)
    r_x = np.matrix([[1.0,0.0,0.0],
                     [0.0,np.cos(pitch),-1*np.sin(pitch)],
                     [0.0,np.sin(pitch),np.cos(pitch)]])
    
    # rotate corner_p by r_x, add back cam x,y,z offsets
    p_out = np.matmul(corner_p, r_x) + cam_pt

    # GEOMETRY
    # Set Sympy 3D point for the camera and a 3D plane for intersection
    cam_sp = spg.Point3D(0, 0, sensor.flyH[0])
    plane = spg.Plane(spg.Point3D(0, 0, 0),
                              normal_vector=(0,0,1))
    
    # blank array for footprint intersection coords
    inter_points = np.zeros((corners.shape[0],2))
    
    # for each sensor corner point
    idx_b = 0
    for pt in np.asarray(p_out):
        
        # create a Sympy 3D point and create a Sympy 3D ray from 
        #   corner point through camera point
        pt_sp = spg.Point3D(pt[0],pt[1],pt[2])
        ray = spg.Ray3D(pt_sp,cam_sp)
        
        # calculate the intersection of the ray with the plane                
        inter_pt = plane.intersection(ray)
        
        # Extract out the X,Y coords fot eh intersection point
        #   ground intersect points will be in this order (LL,UL,UR,LR)
        inter_points[idx_b,0] = inter_pt[0].x.evalf()
        inter_points[idx_b,1] = inter_pt[0].y.evalf()
        
        idx_b += 1
    
        # append inter_points to footprints as a matplotlib path object
        footprint.path[0] = mplPath.Path(inter_points)
    
    # calculate the principle point by intersecting the corners of the ifov path
    ll_pt = spg.Point(inter_points[0,0],inter_points[0,1])
    ul_pt = spg.Point(inter_points[1,0],inter_points[1,1])
    ur_pt = spg.Point(inter_points[2,0],inter_points[2,1])
    lr_pt = spg.Point(inter_points[3,0],inter_points[3,1])
    line_ll_ur = spg.Line(ll_pt,ur_pt)
    line_lr_ul = spg.Line(lr_pt,ul_pt)
    pp_inter = line_ll_ur.intersection(line_lr_ul)
    footprint.pp_x = pp_inter[0].x.evalf()
    footprint.pp_y = pp_inter[0].y.evalf()  
        
    return footprint
# END - def footprint

def graphics(sensor,ifov):
    '''
    def grapics - displays 2 plots containing a verticle profile view of the 
    camera angle and a planimetric view of the footprint
    '''
    
    # copy the path object containing the footprint coordinates, 
    path = ifov.path[0]

    # create a filled triagle object for the vert cross-section plot
    fov = np.array([[0,sensor.flyH],[path.vertices[0][1],0],[path.vertices[1][1],0]])
    tri = mpatches.Polygon(fov,facecolor='b', alpha=0.25)
    
    # plot setup
    fig = plt.figure(figsize = plt.figaspect(0.5))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    
    # build the verticle profile subplot
    ax1.add_patch(tri)
    ax1.scatter(0,sensor.flyH,label='Camera')
    ax1.scatter(ifov.pp_y,0, color='r',label='Pric. Point')
    ax1.hlines(0,path.vertices[0][1]-50,path.vertices[1][1]+50)
    ax1.vlines(0,0,sensor.flyH,'r','--')
    ax1.axis('equal')
    ax1.grid()
    ax1.set_title("Vertical Cross-section")
    ax1.legend()
    
    # build the planimetric subplot
    ifov_patch = mpatches.PathPatch(path, facecolor='b', alpha=0.25)
    ax2.add_patch(ifov_patch)
    ax2.scatter(0,0,label='Camera')
    ax2.scatter(ifov.pp_x,ifov.pp_y, color='r',label='Pric. Point')
    ax2.axis('equal')
    ax2.grid()
    ax2.set_title("Planimetric View (Top)")
# END - graphics 

def resolution(sensor,ifov):

    # create blank data frame for resolutions and break out ifov points
    res = pd.DataFrame(np.zeros((2,4)),columns=['near','mid','far','area'])
    ifov_path = ifov.path[0]
    
    # far field resolution, calc ground distance between far IFOV points
    #   divided by x pixel count
    res.far[0] = euclid_dist(ifov_path.vertices[1][0],ifov_path.vertices[1][1],
                           ifov_path.vertices[2][0],ifov_path.vertices[2][1])
    
    res.far[1] = res.far[0] / sensor.pixel_x[0]

    # near field resolution
    res.near[0] = euclid_dist(ifov_path.vertices[0][0],ifov_path.vertices[0][1],
                           ifov_path.vertices[3][0],ifov_path.vertices[3][1])
    
    res.near[1] = res.near[0] / sensor.pixel_x[0]
    
    # mid-field (principle point) resolution
    #   trig to calculate width of the ifov at the principle point
    pp_slantdist = euclid_dist(0.0,float(sensor.flyH[0]),float(ifov.pp_y[0]),0.0)
    
    res.mid[0] = 2 * (pp_slantdist * np.tan(0.5*np.radians(ifov.fov_h[0])))
    
    res.mid[1] = res.mid[0] / sensor.pixel_x[0]

    # calculate area covered by the fov
    h = ifov_path.vertices[2][1] - ifov_path.vertices[0][1]
    res.area[0] = ((res.near[0]+res.far[0])/2) * h    

    return res
# END - resolution
    
def euclid_dist(x1,y1,x2,y2):
    
    # Simple euclidean distance calculator
    return np.sqrt((x1-x2)**2+(y1-y2)**2)
    
# END - euclid_dist

def write_output(file,sensor,ifov,res):
    
    f = open(file, mode = 'w')
    
    f.write('Resolution Calculations Report\n\n')
    f.write('Camera Type:\t%s\nFocal Length (mm):\t%0.2f\n\n'
            %(sensor.Name[0],sensor.focal[0]))
    f.write('Resolutions calculated for...\n')
    f.write('Flying Height (m):\t%0.1f\nOff-nadir Angle(deg):\t%0.2f\n\n'
            %(sensor.flyH[0], sensor.angle[0]))
    f.write('The IFOV covers:\n\t%0.2f sq. meters / %0.2f ha / %0.2f sq km\n\n' 
            %(res.area[0],res.area[0] * 0.0001, res.area[0] * 1e-6))
    f.write('The IFOV pixel resolutions (m | cm)...\n')
    f.write('Near Field =\t%0.3f | %0.2f\n' %(res.near[1],res.near[1]/0.01))
    f.write('Mid Field  =\t%0.3f | %0.2f\n' %(res.mid[1],res.mid[1]/0.01))
    f.write('Far Field  =\t%0.3f | %0.2f\n' %(res.far[1],res.far[1]/0.01))
    
    f.close()
# END - write_output
    
# MAIN
def main():
    
    app = tk.Tk() # setup tkinter window
    
    # target points - as CSV point cloud (x,y,z,w_surf,r,g,b) from CloudCompare
    #   will be read in 7500 point chunks for memory managment purposes   
    sensor_file = fd.askopenfilename(title='Open Point Cloud',
                                     filetypes=[('Comma-Delimited Files (*.csv)',
                                     '*.csv')],initialdir=os.getcwd())
    sensor = pd.read_csv(sensor_file)
    
    # clean up the tkinter app
    app.destroy()
    
    # call footprint func to calc the corner coordinates of the IFOV
    ifov = footprint(sensor)
    
    # display graphs of the ifov 
    graphics(sensor,ifov)
    
    # call resolution func to calculate the resolution
    res = resolution(sensor, ifov)

    # set output file name
    out_file = sensor_file[:-4] + "_resolution" + ".txt"
    
    # func to write output into a nicely formatted text file
    write_output(out_file, sensor, ifov, res)
    
    print("...PROCESSING COMPLETE...")
    

if __name__ == "__main__":
    main()
