# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 23:54:12 2019

@author: user
"""

from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import numpy as np
import math
import time
from collections import Counter
import sys    
def fEDC(part, model, singlePlane, multiPlane, *args, **kwargs):
    unit = kwargs.get('unit')
    mat = kwargs.get('mat')
    cusMat = kwargs.get('cusMat')
    N = kwargs.get('N')
    g0 = kwargs.get('g0')
    g90 = kwargs.get('g90')
    cusSigma = kwargs.get('cusSigma')
    selectAxis = kwargs.get('selectAxis')
    axisGlobal = kwargs.get('axisGlobal')
    if part == '':
        raise Exception('Part not defined please make a part before using software')
    if N == None:
        raise Exception('Enter number of fracture planes')
    if g0 == None or g90 == None:
        raise Exception('Enter properties of material')
    p = mdb.models[model].parts[part]
    
    # timing script
    t1 = time.time()
    # creating plane
    area2 = {}
    U = {}
    F ={}
    pos = {}
    fileNumber = 1
    axis = []
    plotU = []
    plotF = []
    x = []
    y = []
    z = []
    nodeSet = p.nodes
    elementSet = p.elements
    if len(elementSet) == 0 or len(elementSet) == None:
        raise Exception('No mesh defined please define mesh before using software')
    #setting material properties
    if cusMat == 'Material from database':
        if mat == 'T800s/M21':
            G90 = 0.255
            G0 = 209
            sigma1 = 3066.96
            
        elif mat == 'T300/913':
            G90 = 0.211
            G0 = 133
            sigma1 = 3530
        elif mat == 'T300/920':
            G90 = 0.456
            G0 = 132
            sigma1 = 3530
        elif cusMat == True:
            G90 = g90
            G0 = g0
    if unit == 'millimeter':
        G90 = G90/1000000
        G0 = G0/1000000
        sigma1 = sigma1/1000000
    
    if cusMat == 'Custom Material':
        G90 = g90
        G0 = g0
        sigma1 = cusSigma

    
    if multiPlane == True and selectAxis == 'Use global axis':
        for n in range(0,len(nodeSet)):
            x.append( nodeSet[n].coordinates[0])
            y.append(nodeSet[n].coordinates[1])
            z.append(nodeSet[n].coordinates[2])
    
        #finding maximum and minimum coordinates of the model
        minX = min(x)
        minY = min(y)
        minZ = min(z)
        maxX = max(x)
        maxY = max(y)
        maxZ = max(z)
        
        if axisGlobal == 'X':
            start = minX
            end = maxX
        elif axisGlobal == 'Y':
            start = minY
            end = maxY
        elif axisGlobal == 'Z': 
            start = minZ
            end = maxZ
        
        for z in np.arange(start,end,abs(start-end)/N):
            if axisGlobal == 'X':
                s = np.array([z,0,0])
                q = np.array([z,0,1])
                r = np.array([z,1,0])
                ax = 'x'
            elif axisGlobal == 'Y':
                s = np.array([0,z,0])
                q = np.array([0,z,1])
                r = np.array([1,z,0])
                ax = 'y'
            elif axisGlobal == 'Z':
                s = np.array([0,0,z])
                q = np.array([0,1,z])
                r = np.array([1,0,z])
                ax = 'z'
            lambRecord = []
        
            p.DatumPlaneByThreePoints(point1=s, 
                point2=q, 
                point3=r) 
            
            pq = q - s
            pr = r - s
            
            perp = np.cross(pq,pr)
            d = -np.dot(perp,-s)
            #select elements
            vec = np.zeros((12,3))
            vertexply = np.zeros((8,3))
            vertexplyi = np.zeros((12,3))
            vertexplyj = np.zeros((12,3))
            vertexsubply = np.zeros((8,3))
            #finding area of intercepts in the general case of any plane interception on the element
            for e in range(0,len(elementSet)):
        #        if e == 3:
                connected = []
                    
                connected = elementSet[e].connectivity
                label = elementSet[e].label
                area2[label] = []
                    # sort node index to keep track of lines and vertex
                sort = np.sort(connected)
                side1 = np.array(nodeSet[sort[0]].coordinates)-np.array(nodeSet[sort[1]].coordinates)
                side2 = np.array(nodeSet[sort[0]].coordinates)-np.array(nodeSet[sort[2]].coordinates)
                side3 = np.array(nodeSet[sort[0]].coordinates)-np.array(nodeSet[sort[4]].coordinates)
                size = [np.linalg.norm(side1)]
                size.append(np.linalg.norm(side2))
                size.append(np.linalg.norm(side3))
                maxSize = max(size)
                diag = math.sqrt(maxSize**2+maxSize**2)
                # calculate the elements near the plane taking account element mesh size
                if float(abs(np.dot(perp,np.array(nodeSet[sort[0]].coordinates)) - d)/(np.linalg.norm(perp))) <= float(math.sqrt(diag**2+maxSize**2)):
                    region = []
                    nLayups = len(p.compositeLayups.keys())
                    if nLayups == 0:
                        raise Exception('Composite layup not defined please define layup before using software')
                    # number of layups
                    for n in range(0,nLayups):
                        key = p.compositeLayups.keys()[n]
                        nPlies = len(p.compositeLayups[key].plies)
                        #ply stack direction
                        stackAxis = p.compositeLayups[key].orientation.stackDirection
                        plyThickness = []
                        # calculate of the total relative thickness of plies add to one otherwise 
                        #terminate the program
                        for j in range(0,nPlies):
                            plyThickness.append(p.compositeLayups[key].plies[j].thickness)
                            region.append(p.compositeLayups[key].plies[j].region[0])
                        if sum(plyThickness) != 1 and region != [region[0]] * nPlies:
                            raise Exception('Defined relative thiknesses do not add up to one and/or regions defined for compositelayup are different')
                        plyRelThickness = 0  
                        # Three possible directions for the ply stack
                        # Vertex change depending on stack direction
                        if str(stackAxis) == 'STACK_1':
                            if np.array(nodeSet[sort[4]].coordinates)[0] < np.array(nodeSet[sort[0]].coordinates)[0]:
                                vertexsubply[0,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[5]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[6]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[7]].coordinates)
                            
                                vertexsubply[4,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[3]].coordinates)
                            else:
                                vertexsubply[4,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[5]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[6]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[7]].coordinates)
                            
                                vertexsubply[0,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[3]].coordinates)
                        elif str(stackAxis) == 'STACK_2':
                            if np.array(nodeSet[sort[1]].coordinates)[1] < np.array(nodeSet[sort[0]].coordinates)[1]:
                                vertexsubply[0,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[3]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[5]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[7]].coordinates)
                                
                                vertexsubply[4,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[6]].coordinates)
                            else:
                                vertexsubply[4,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[3]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[5]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[7]].coordinates)
                                
                                vertexsubply[0,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[6]].coordinates)
                        elif str(stackAxis) == 'STACK_3':
                            if np.array(nodeSet[sort[2]].coordinates)[2] < np.array(nodeSet[sort[0]].coordinates)[2]:
                                vertexsubply[0,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[3]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[6]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[7]].coordinates)
                                
                                vertexsubply[4,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[5]].coordinates)
                            
                            else:
                                vertexsubply[4,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[3]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[6]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[7]].coordinates)
                                
                                vertexsubply[0,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[5]].coordinates)
                       # assign vertex to different variable for further calculation
                        for i in range(0,4):     
                                vertexply[i,:] = vertexsubply[i,:]
                        for j in range(0,nPlies):
                            plyRelThickness += p.compositeLayups[key].plies[j].thickness
                            lamb = []
                            
                            #assign ply coordinates according to ply thickness using equation of a line 
                            for i in range(0,4):
                                vertexply[i+4,:] = vertexsubply[i,:] + plyRelThickness * (vertexsubply[i+4,:] - vertexsubply[i,:])
                                vertexplyi[i,:] = vertexply[i,:]
                                vertexplyj[i,:] = vertexply[i+4,:]
                                
                            vertexplyi[4,:] = vertexply[0,:] 
                            vertexplyj[4,:] = vertexply[1,:]
                            
                            vertexplyi[5,:] = vertexply[3,:]
                            vertexplyj[5,:] = vertexply[1,:]
                            
                            vertexplyi[6,:] = vertexply[2,:]
                            vertexplyj[6,:] = vertexply[0,:]
                            
                            vertexplyi[7,:] = vertexply[2,:] 
                            vertexplyj[7,:] = vertexply[3,:] 
                            
                            vertexplyi[8,:] = vertexply[4,:]
                            vertexplyj[8,:] = vertexply[5,:]
                            
                            vertexplyi[9,:] = vertexply[7,:]
                            vertexplyj[9,:] = vertexply[5,:]
                            
                            vertexplyi[10,:] = vertexply[6,:]
                            vertexplyj[10,:] = vertexply[4,:]
                            
                            vertexplyi[11,:] = vertexply[6,:]
                            vertexplyj[11,:] = vertexply[7,:]
                            
            #                # Algorithm for lines and finding lambda (gradient of intercept)
            #                # refer to paper
                            
                            for i in range(0,len(vec)):
                                vec[i,:] = vertexplyi[i,:] - vertexplyj[i,:]
                                if np.dot(perp, vec[i,:]) != 0:
                                    lamb.append(round((d - np.dot(perp, vertexplyj[i,:]))/(np.dot(perp, vec[i,:])),3)) 
                                else:
                                    lamb.append(-100)
                            k = 0  
                            intercept = np.zeros((12,3))
                            count = Counter(lamb)
            #                print(vec, vertexplyj)
            #                print(vertexplyi)
            #                print(lamb)
                            # checking for every value of lambda if there is an intercept
                            for i in range(0,len(lamb)):
                                if 0 < lamb[i] < 1 and count[0] >= 4 or count[1] >= 4:
                                    continue
                                elif lamb[i] < 1 and lamb[i] >= 0:
                                    #if there is an intercept find the coordinates using parametric method
                                    temp = (np.multiply(vec[i,:],lamb[i]) + vertexplyj[i,:])
                                    intercept[k,:] = temp
                                    #if lambda = 1 check if there are over lapping intersections
                                    if lamb[i] == 1:
                                        for j in range(0,len(intercept)):
                                            if j < k and np.array_equal(intercept[j,:], intercept[k,:]):
                                                k -= 1
                                    k += 1
                                # check if intersection is at the nodes and if the itersections lie
                                # diagonally directly at the nodes so we have lambda = 1 and 0 
                                # the same number of times
                                elif lamb[i] == 1 and count[0] == count[1]:
                                    temp = (np.multiply(vec[i,:],lamb[i]) + vertexplyj[i,:])
                                    intercept[k,:] = temp
                                    if lamb[i] == 0:
                                        for j in range(0,len(intercept)):
                                            if j < k and np.array_equal(intercept[j,:], intercept[k,:]):
                                                k -= 1
                                    k += 1
                                    
                                if lamb[i] <=1 and lamb[i] >= 0:
                                    lambRecord.append(lamb[i])
            #                print(intercept)
                            #find the areas by sorting the coordinates of intersections
                            norm = 0
                            angles = {}
            
            #                        print(lamb)
                            #for just three points of intersection find area of triangle
                            if k == 3:
                                pq = intercept[0,:] - intercept[1,:]
                                pr = intercept[0,:] - intercept[2,:]
                                cros = np.cross(pq,pr)
                                norm += np.linalg.norm(cros)/2
                                area2[label].append(norm)
                            #things become harder with more than 3 points
                            elif 4 <= k <= 6:
                                #define reference vector
                                ref = intercept[0,:] - intercept[1,:]
                                #set the reference angle of the reference vector to zero
                                angles[1] = 0
                                for i in range(2,k):
                                    #find all other vectors with respect to the first point
                                    #find their angles with respect to the reference angle
                                    #then arrange from smallest to largest angle
                                    ref2 = intercept[0,:] - intercept[i,:]
                                    cros = np.cross(ref,ref2)
                                    refNorm = np.linalg.norm(cros)
                                    product = np.linalg.norm(ref)*np.linalg.norm(ref2)
                                    dot = np.dot(ref,ref2)
                                    #if statement added as to test which side the point lies
                                    #on reference vector. If the z componenet of the cross
                                    #product of the reference vector and the new vector
                                    #is negative then the point lies on one side, if positive
                                    #point lies on other side
                                    if cros[2] < 0:
                                        sin = round(-refNorm/product,5)
                                    else:
                                        sin = round(refNorm/product,5)
                                    cos = round(dot/product,5)
                                    #between 0 and 90
                                    if sin >= 0 and cos >= 0:
                                        theta = (180*math.acos(cos))/math.pi
                                        angles[i] = theta
                                    #between -90 and 0
                                    elif sin <= 0 and cos >= 0:
                                        theta = (180*math.asin(sin))/math.pi
                                        angles[i] = theta
                                    #between 90 and 180
                                    elif sin >= 0 and cos <= 0:
                                        theta = 180-((180*math.asin(sin))/math.pi)
                                        angles[i] = theta
                                    #between -180 and -90 
                                    elif sin <= 0 and cos <= 0:
                                        theta = -180 + ((180*math.asin(sin))/math.pi)
                                        angles[i] = theta
                                sorted_angles = sorted(angles.items(), key=lambda x: x[1])
                                sorted_index = np.array(sorted_angles)[:,0]
                                for i in range(0,len(sorted_index)-1):
                                    pq = intercept[0,:] - intercept[sorted_index[i],:]
                                    pr = intercept[0,:] - intercept[sorted_index[i+1],:]
                                    cros = np.cross(pq,pr)
                                    norm += np.linalg.norm(cros)/2
                                area2[label].append(norm)
                            else:
                                area2[label].append(0)
                            
                            #assign coordinates of last four vertex to first four vertex
                            #in order to calcualte area of the new ply
                            for i in range(0,4):
                                vertexply[i,:] = vertexply[i+4,:]
                           #5.7 seconds when added to only consider elements near the plane and 28.4 seconds without modifycation 
            if lambRecord[0:] == [1] * len(lambRecord) or lambRecord[0:] == [0] * len(lambRecord):
                        area2 = {}
            for label in area2.keys():
                if area2[label] == [0] * nPlies or area2[label] == []:
                    del area2[label]
            #    print(vertexsubply)        
            #    print('Area calculated from intersecting planes: {0}'.format(area2))
                #print('Area calculated for special case where plane intersected nodes: {0}'.format(area))
        
    
            if area2 != {}:
            # extracting material orientation
                Name = 'csys-plane'
                if Name in p.features.keys():
                    del p.features[Name]
                #creating coordinate system for plane
                p.DatumCsysByThreePoints( origin = s, name = Name,coordSysType = CARTESIAN,point1 = q, point2 = s + perp)
               #extracting coordinate axis
                planeAxis1 = p.datums[p.features[Name].id].axis1
                planeAxis2 = p.datums[p.features[Name].id].axis2
                planeAxis3 = p.datums[p.features[Name].id].axis3
                
                
            
        #        Rx = np.zeros((3,3))
                Ry = np.zeros((3,3))
                Rz = np.zeros((3,3))
                materialAngle = {}
                #for each defined ply layup
                for i in range(0,nLayups):
                    key = p.compositeLayups.keys()[i]
                    csys = p.compositeLayups[key].orientation.localCsys
                    refAxis = str(p.compositeLayups[key].orientation.axis)
                    addAngle = p.compositeLayups[key].orientation.angle
                    #extract layup axis
                    plyaxis1 = p.datums[csys].axis1
        #            plyaxis2 = p.datums[sys].axis2
        #            plyaxis3 = p.datums[sys].axis3
                    #for every ply
                    
                    for j in range(0,nPlies):
                        plyOrientation = p.compositeLayups[key].plies[j].orientation
                        plyAngle = p.compositeLayups[key].plies[j].orientationValue
                        plyOrient = str(p.compositeLayups[key].plies[j].axis)
                        plyOrientType = str(p.compositeLayups[key].plies[j].orientationType)
                        
                        #if no material orientation is defined for each individual py then 
                        #orientatio will be that of the layup
                        if plyOrientation == None:
                            direction1 = plyaxis1.direction
        #                    direction2 = plyaxis2.direction
        #                    direction3 = plyaxis3.direction
                            totalAngle = addAngle + plyAngle
                            
                            if plyOrientType == 'ANGLE_45':
                                totalAngle = addAngle + 45
                            elif plyOrientType == 'ANGLE_90':
                                totalAngle = addAngle + 90
                            elif plyOrientType == 'ANGLLE_NEG45':
                                totalAngle = addAngle - 45
                                
                                
                            if totalAngle != 0: 
        #                        if refAxis == 'AXIS_1':
        #                    #rotation wrt x axis
        #                            Rx = np.array([[1,0,0],[0,math.cos(totalAngle*math.pi/180),math.sin(totalAngle*math.pi/180)],[0,-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180)]])
        #                            direction2 = np.dot(Rx,plyaxis2.direction)
        #                            direction3 = np.dot(Rx,plyaxis3.direction)
                                   
                                if refAxis == 'AXIS_2':
                                    Ry = np.array([[math.cos(totalAngle*math.pi/180),0,-math.sin(totalAngle*math.pi/180)],[0,1,0],[math.sin(totalAngle*math.pi/180),0,math.cos(45*math.pi/180)]])
                                    direction1 = np.dot(Ry,plyaxis1.direction)
        #                            direction3 = np.dot(Ry,plyaxis3.direction)
                               
                                elif refAxis == 'AXIS_3':
                                    Rz = np.array([[math.cos(totalAngle*math.pi/180),-math.sin(totalAngle*math.pi/180),0],[-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180),0],[0,0,1]])
                                    direction1 = np.dot(Rz,plyaxis1.direction)
        #                            direction2 = np.dot(Rz,plyaxis2.direction)
                        else:
                            #if there is a defines Csys for each ply then change the axis defined
                            plyaxis1 = p.datums[plyOrientation].axis1
        #                    plyaxis2 = p.datums[plyOrientation].axis2
        #                    plyaxis3 = p.datums[plyOrientation].axis3
                        
                            direction1 = plyaxis1.direction
        #                    direction2 = plyaxis2.direction
        #                    direction3 = plyaxis3.direction
            #               
                            totalAngle = plyAngle
                            if plyOrientType == 'ANGLE_45':
                                totalAngle = 45
                            elif plyOrientType == 'ANGLE_90':
                                totalAngle = 90
                            elif plyOrientType == 'ANGLLE_NEG45':
                                totalAngle = -45
                                
        #                    if plyOrient == 'AXIS_1':
        #                    #rotation wrt x axis
        #                        Rx = np.array([[1,0,0],[0,math.cos(totalAngle*math.pi/180),math.sin(totalAngle*math.pi/180)],[0,-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180)]])
        #                        direction2 = np.dot(Rx,plyaxis2.direction)
        #                        direction3 = np.dot(Rx,plyaxis3.direction)
                               
                            elif plyOrient == 'AXIS_2':
                                Ry = np.array([[math.cos(totalAngle*math.pi/180),0,-math.sin(totalAngle*math.pi/180)],[0,1,0],[math.sin(totalAngle*math.pi/180),0,math.cos(45*math.pi/180)]])
                                direction1 = np.dot(Ry,plyaxis1.direction)
        #                        direction3 = np.dot(Ry,plyaxis3.direction)
                                
                            elif plyOrient == 'AXIS_3':
                                Rz = np.array([[math.cos(totalAngle*math.pi/180),-math.sin(totalAngle*math.pi/180),0],[-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180),0],[0,0,1]])
                                direction1 = np.dot(Rz,plyaxis1.direction)
        #                        direction2 = np.dot(Rz,plyaxis2.direction)
            
            #            
                        #finding angle between direction of plane axis and fibre axis    
                        
        #                dotx = np.dot(direction1,planeAxis1.direction)
                        doty = np.dot(direction1,planeAxis2.direction)
        #                dotz = np.dot(direction1,planeAxis3.direction)
        #                productx = np.linalg.norm(direction1)*np.linalg.norm(planeAxis1.direction)
                        producty = np.linalg.norm(direction1)*np.linalg.norm(planeAxis2.direction)
        
        
                        cosy = doty/producty
        
          
                        theta = 180*math.acos(cosy)/math.pi
                        
        
                        if theta > 90:
                            theta = 180 - theta
                                
                        #assigning angles to elements
                        plyRegion = p.compositeLayups[key].plies[j].region[0]
                        plyElements = p.sets[plyRegion].elements
                        for k in range(0,len(plyElements)):
                            label = plyElements[k].label
                            if label in area2.keys() and label in materialAngle.keys(): #check if another ply is assigned to the same element
                                materialAngle[label].append(theta) #append to the orientation angle of the element
                                    #assume plys have same relative thickness
                            elif label in area2.keys():
                                materialAngle[label] = [theta]
                                    #if only one ply is assigned to the element
                
            
        #        print('Angle output for each element with respect to x, y and z axis of cut plane: {0}'.format(materialAngle))
                
                #calculating toughness of composite in relation to fibre orientation
                #Then storing the toughness of composite in a database if existing 
                #orientation already exist
                materialG = {}
                totalU = 0
                materialSigma1 = {}
                totalF = 0
                axis.append(z)
                for element in area2.keys():
                    theta = materialAngle[element]
                    #writing  and reading to database
                    materialG[element] = G0 * np.cos(np.multiply(theta,np.pi/180)) + G90 * np.sin(np.multiply(theta,np.pi/180))
                    materialSigma1[element] = sigma1 * np.cos(np.multiply(theta,np.pi/180))
                    #calculate energy dissipation
                    
                    totalU += np.dot(materialG[element],area2[element])
                    totalF += np.dot(materialSigma1[element],area2[element])
                if  totalU not in U.keys():  
                    U[totalU] = [p.features.keys()[-2]]
                    F[totalU] = [totalF]
                    pos[totalU] = [z]
                else:
                    U[totalU].append(p.features.keys()[-2])
                    F[totalU].append(totalF)
                    pos[totalU].append(z)
                print('Total energy dissipated from failure: {0} kJ at {1}, = {2}, Critical force: {3} MN'.format(float(totalU),ax,z,float(totalF)))
                plotU.append(totalU)
                plotF.append(totalF)
        critU = min(U.keys())
        planes = U[critU]
        force = F[critU]
        critPos = pos[critU]      
        print('Critical energy is: {0} kJ. Critical Force: {1} MN at the following planes: {2}, at {3} = {4}'.format(critU,force,planes,ax,critPos))
        t2 = time.time()
        print('Run time: {0}'.format(t2-t1))
        
        #plotting output energies with increasing co-ordinates
        #getting existing xyplots
        outputFiles = session.xyPlots.keys()
        fileName = 'Energy Output{0}'.format(fileNumber)
        while fileName in outputFiles:
            fileNumber += 1
            fileName = 'Energy Output{0}'.format(fileNumber)
        #creating new data plot
        xyp = session.XYPlot('Energy Output{0}'.format(fileNumber))
        #creating new chart
        chartName = xyp.charts.keys()[0]
        chart = xyp.charts[chartName]
        #x and y label
        yQuantity = visualization.QuantityType(type=ENERGY)
        c = np.empty((len(axis),2))
        c[:,0] = axis
        c[:,1] = plotU
        #creating data plot
        xy2 = xyPlot.XYData(data=(c), sourceDescription='Entered from keyboard', 
             axis2QuantityType=yQuantity, )
        #creating curve
        c2 = session.Curve(xyData=xy2)
        #plotting curve
        chart.setValues(curvesToPlot=(c2, ), )
        session.charts[chartName].axes1[0].axisData.setValues(useSystemTitle=False, 
        title='{0} coordinate'.format(ax))
        #changing view to graph
        session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
 
    
    # force plot
        fileNameF = 'Force Output{0}'.format(fileNumber)
        while fileName in outputFiles:
            fileNumber += 1
            fileNameF = 'Force Output{0}'.format(fileNumber)
        xyp = session.XYPlot(fileNameF)
    
        chartName = xyp.charts.keys()[0]
        chart = xyp.charts[chartName]
        c = np.empty((len(axis),2))
        c[:,0] = axis
        c[:,1] = plotF
        xy2 = xyPlot.XYData(data=(c), sourceDescription='Entered from keyboard', 
            )
        c2 = session.Curve(xyData=xy2)
        chart.setValues(curvesToPlot=(c2, ), )
        session.charts[chartName].axes2[0].axisData.setValues(useSystemTitle=False, 
        title='Force')
        session.charts[chartName].axes1[0].axisData.setValues(useSystemTitle=False, 
        title='{0} coordinate'.format(ax))
        session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
        
        
    #if ownaxis selected 
    if multiPlane == True and selectAxis == 'Use two points':
        axisPoints = kwargs.get('axisPoints')
        # timing script
        t1 = time.time()
        # creating plane
        area2 = {}
        U = {}    
        pos = {}
        axis = []
        plotU = []
        if axisPoints == 'Select 2 Points (Nodes In Mesh) To Create Axis':
            axisPointStart = kwargs.get('axisPointStart') 
            axisPointEnd = kwargs.get('axisPointEnd')
            if axisPointStart == None or axisPointEnd == None:
                raise Exception('Points not selected')
            s = np.array(axisPointStart.coordinates)
            q = np.array(axisPointEnd.coordinates)
        elif axisPoints == 'Define Points Manually':
            coordAxisStart = kwargs.get('coordAxisStart')
            coordAxisEnd = kwargs.get('coordAxisEnd')
            print(coordAxisStart)
            if coordAxisStart == () or coordAxisEnd == ():
                raise Exception('Points not defined')
            elif len(coordAxisStart[0]) != 3 or len(coordAxisEnd[0]) != 3:
                raise Exception('Axis points not defined properly')
            s = np.array([coordAxisStart[0][0],coordAxisStart[0][1],coordAxisStart[0][2]])
            q = np.array([coordAxisEnd[0][0],coordAxisEnd[0][1],coordAxisEnd[0][2]])
            
            
        perp = np.subtract(q,s)
        axisLength = np.linalg.norm(perp)
        p.DatumAxisByTwoPoint(point1=s, 
            point2=q)
        d1 = p.datums
        planeName = p.features.keys()[-1]
        inc = 1/float(N)
        for z in np.arange(0,1, inc):
            lambRecord = []
            pt2 = s + np.multiply(perp,z)
            d = np.dot(perp,pt2)
    
            p.DatumPlaneByPointNormal(normal=d1[p.features[planeName].id], point=s + np.multiply(perp,z))
            
            #select elements
            vec = np.zeros((12,3))
            vertexply = np.zeros((8,3))
            vertexplyi = np.zeros((12,3))
            vertexplyj = np.zeros((12,3))
            vertexsubply = np.zeros((8,3))
            #finding area of intercepts in the general case of any plane interception on the element
            for e in range(0,len(elementSet)):
    
                connected = []
                    
                connected = elementSet[e].connectivity
                label = elementSet[e].label
                area2[label] = []
                    # sort node index to keep track of lines and vertex
                sort = np.sort(connected)
                side1 = np.array(nodeSet[sort[0]].coordinates)-np.array(nodeSet[sort[1]].coordinates)
                side2 = np.array(nodeSet[sort[0]].coordinates)-np.array(nodeSet[sort[2]].coordinates)
                side3 = np.array(nodeSet[sort[0]].coordinates)-np.array(nodeSet[sort[4]].coordinates)
                size = [np.linalg.norm(side1)]
                size.append(np.linalg.norm(side2))
                size.append(np.linalg.norm(side3))
                maxSize = max(size)
                diag = math.sqrt(maxSize**2+maxSize**2)
                # calculate the elements near the plane taking account element mesh size
                if float(abs(np.dot(perp,np.array(nodeSet[sort[0]].coordinates)) - d)/(np.linalg.norm(perp))) <= float(math.sqrt(diag**2+maxSize**2)):
                    region = []
                    nLayups = len(p.compositeLayups.keys())
                    if nLayups == 0:
                        sys.exit('Composite layup not defined please define layup before using software')

                    # number of layups
                    for n in range(0,nLayups):
                        key = p.compositeLayups.keys()[n]
                        nPlies = len(p.compositeLayups[key].plies)
                        #ply stack direction
                        stackAxis = p.compositeLayups[key].orientation.stackDirection
                        plyThickness = []
                        # calculate of the total relative thickness of plies add to one otherwise 
                        #terminate the program
                        for j in range(0,nPlies):
                            plyThickness.append(p.compositeLayups[key].plies[j].thickness)
                            region.append(p.compositeLayups[key].plies[j].region[0])
                        if sum(plyThickness) != 1 and region != [region[0]] * nPlies:
                            sys.exit('Defined relative thiknesses do not add up to one and/or regions defined for compositelayup are different')
                        plyRelThickness = 0  
                        # Three possible directions for the ply stack
                        # Vertex change depending on stack direction
                        if str(stackAxis) == 'STACK_1':
                            if np.array(nodeSet[sort[4]].coordinates)[0] < np.array(nodeSet[sort[0]].coordinates)[0]:
                                vertexsubply[0,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[5]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[6]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[7]].coordinates)
                            
                                vertexsubply[4,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[3]].coordinates)
                            else:
                                vertexsubply[4,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[5]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[6]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[7]].coordinates)
                            
                                vertexsubply[0,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[3]].coordinates)
                        elif str(stackAxis) == 'STACK_2':
                            if np.array(nodeSet[sort[1]].coordinates)[1] < np.array(nodeSet[sort[0]].coordinates)[1]:
                                vertexsubply[0,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[3]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[5]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[7]].coordinates)
                                
                                vertexsubply[4,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[6]].coordinates)
                            else:
                                vertexsubply[4,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[3]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[5]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[7]].coordinates)
                                
                                vertexsubply[0,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[6]].coordinates)
                        elif str(stackAxis) == 'STACK_3':
                            if np.array(nodeSet[sort[2]].coordinates)[2] < np.array(nodeSet[sort[0]].coordinates)[2]:
                                vertexsubply[0,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[3]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[6]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[7]].coordinates)
                                
                                vertexsubply[4,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[5]].coordinates)
                            
                            else:
                                vertexsubply[4,:] = np.array(nodeSet[sort[2]].coordinates)
                                vertexsubply[5,:] = np.array(nodeSet[sort[3]].coordinates)
                                vertexsubply[6,:] = np.array(nodeSet[sort[6]].coordinates)
                                vertexsubply[7,:] = np.array(nodeSet[sort[7]].coordinates)
                                
                                vertexsubply[0,:] = np.array(nodeSet[sort[0]].coordinates)
                                vertexsubply[1,:] = np.array(nodeSet[sort[1]].coordinates)
                                vertexsubply[2,:] = np.array(nodeSet[sort[4]].coordinates)
                                vertexsubply[3,:] = np.array(nodeSet[sort[5]].coordinates)
                       # assign vertex to different variable for further calculation
                        for i in range(0,4):     
                                vertexply[i,:] = vertexsubply[i,:]
                        for j in range(0,nPlies):
                            plyRelThickness += p.compositeLayups[key].plies[j].thickness
                            lamb = []
                            
                            #assign ply coordinates according to ply thickness using equation of a line 
                            for i in range(0,4):
                                vertexply[i+4,:] = vertexsubply[i,:] + plyRelThickness * (vertexsubply[i+4,:] - vertexsubply[i,:])
                                vertexplyi[i,:] = vertexply[i,:]
                                vertexplyj[i,:] = vertexply[i+4,:]
                                
                            vertexplyi[4,:] = vertexply[0,:] 
                            vertexplyj[4,:] = vertexply[1,:]
                            
                            vertexplyi[5,:] = vertexply[3,:]
                            vertexplyj[5,:] = vertexply[1,:]
                            
                            vertexplyi[6,:] = vertexply[2,:]
                            vertexplyj[6,:] = vertexply[0,:]
                            
                            vertexplyi[7,:] = vertexply[2,:] 
                            vertexplyj[7,:] = vertexply[3,:] 
                            
                            vertexplyi[8,:] = vertexply[4,:]
                            vertexplyj[8,:] = vertexply[5,:]
                            
                            vertexplyi[9,:] = vertexply[7,:]
                            vertexplyj[9,:] = vertexply[5,:]
                            
                            vertexplyi[10,:] = vertexply[6,:]
                            vertexplyj[10,:] = vertexply[4,:]
                            
                            vertexplyi[11,:] = vertexply[6,:]
                            vertexplyj[11,:] = vertexply[7,:]
                            
            #                # Algorithm for lines and finding lambda (gradient of intercept)
            #                # refer to paper
                            
                            for i in range(0,len(vec)):
                                vec[i,:] = vertexplyi[i,:] - vertexplyj[i,:]
                                if np.dot(perp, vec[i,:]) != 0:
                                    lamb.append(round((d - np.dot(perp, vertexplyj[i,:]))/(np.dot(perp, vec[i,:])),3)) 
                                else:
                                    lamb.append(-100)
                            k = 0  
                            intercept = np.zeros((12,3))
                            count = Counter(lamb)
            #               
                            # checking for every value of lambda if there is an intercept
                            for i in range(0,len(lamb)):
                                if 0 < lamb[i] < 1 and count[0] >= 4 or count[1] >= 4:
                                    continue
                                elif lamb[i] < 1 and lamb[i] >= 0:
                                    #if there is an intercept find the coordinates using parametric method
                                    temp = (np.multiply(vec[i,:],lamb[i]) + vertexplyj[i,:])
                                    intercept[k,:] = temp
                                    #if lambda = 1 check if there are over lapping intersections
                                    if lamb[i] == 1:
                                        for j in range(0,len(intercept)):
                                            if j < k and np.array_equal(intercept[j,:], intercept[k,:]):
                                                k -= 1
                                    k += 1
                                # check if intersection is at the nodes and if the itersections lie
                                # diagonally directly at the nodes so we have lambda = 1 and 0 
                                # the same number of times
                                elif lamb[i] == 1 and count[0] == count[1]:
                                    temp = (np.multiply(vec[i,:],lamb[i]) + vertexplyj[i,:])
                                    intercept[k,:] = temp
                                    if lamb[i] == 0:
                                        for j in range(0,len(intercept)):
                                            if j < k and np.array_equal(intercept[j,:], intercept[k,:]):
                                                k -= 1
                                    k += 1
                                    
                                if abs(lamb[i]) <=1 and lamb[i] >= 0:
                                    lambRecord.append(lamb[i])
                            #find the areas by sorting the coordinates of intersections
                            norm = 0
                            angles = {}
    
                            #for just three points of intersection find area of triangle
                            if k == 3:
                                pq = intercept[0,:] - intercept[1,:]
                                pr = intercept[0,:] - intercept[2,:]
                                cros = np.cross(pq,pr)
                                norm += np.linalg.norm(cros)/2
                                area2[label].append(norm)
                            #things become harder with more than 3 points
                            elif 4<= k <= 6:
                                #define reference vector
                                ref = intercept[0,:] - intercept[1,:]
                                #set the reference angle of the reference vector to zero
                                angles[1] = 0
                                for i in range(2,k):
                                    #find all other vectors with respect to the first point
                                    #find their angles with respect to the reference angle
                                    #then arrange from smallest to largest angle
                                    ref2 = intercept[0,:] - intercept[i,:]
                                    cros = np.cross(ref,ref2)
                                    refNorm = np.linalg.norm(cros)
                                    product = np.linalg.norm(ref)*np.linalg.norm(ref2)
                                    dot = np.dot(ref,ref2)
                                    #if statement added as to test which side the point lies
                                    #on reference vector. If the z componenet of the cross
                                    #product of the reference vector and the new vector
                                    #is negative then the point lies on one side, if positive
                                    #point lies on other side
                                    if cros[2] < 0:
                                        sin = round(-refNorm/product,5)
                                    else:
                                        sin = round(refNorm/product,5)
                                    cos = round(dot/product,5)
                                    #between 0 and 90
                                    if sin >= 0 and cos >= 0:
                                        theta = (180*math.acos(cos))/math.pi
                                        angles[i] = theta
                                    #between -90 and 0
                                    elif sin <= 0 and cos >= 0:
                                        theta = (180*math.asin(sin))/math.pi
                                        angles[i] = theta
                                    #between 90 and 180
                                    elif sin >= 0 and cos <= 0:
                                        theta = 180-((180*math.asin(sin))/math.pi)
                                        angles[i] = theta
                                    #between -180 and -90 
                                    elif sin <= 0 and cos <= 0:
                                        theta = -180 + ((180*math.asin(sin))/math.pi)
                                        angles[i] = theta
                                sorted_angles = sorted(angles.items(), key=lambda x: x[1])
                                sorted_index = np.array(sorted_angles)[:,0]
                                for i in range(0,len(sorted_index)-1):
                                    pq = intercept[0,:] - intercept[sorted_index[i],:]
                                    pr = intercept[0,:] - intercept[sorted_index[i+1],:]
                                    cros = np.cross(pq,pr)
                                    norm += np.linalg.norm(cros)/2
                                area2[label].append(norm)
                            else:
                                area2[label].append(0)
                            
                            #assign coordinates of last four vertex to first four vertex
                            #in order to calcualte area of the new ply
                            for i in range(0,4):
                                vertexply[i,:] = vertexply[i+4,:]
                           #5.7 seconds when added to only consider elements near the plane and 28.4 seconds without modifycation 
            if lambRecord[0:] == [1] * len(lambRecord) or lambRecord[0:] == [0] * len(lambRecord):
                            area2 = {}
            for label in area2.keys():
                if area2[label] == [0] * nPlies or area2[label] == []:
                    del area2[label]
            fileNumber = 0
            if area2 != {}:
            # extracting material orientation
                Name = 'csys-plane'
                if Name in p.features.keys():
                    del p.features[Name]
                #creating coordinate system for plane   
        #        Rx = np.zeros((3,3))
                Ry = np.zeros((3,3))
                Rz = np.zeros((3,3))
                materialAngle = {}
                #for each defined ply layup
                for i in range(0,nLayups):
                    key = p.compositeLayups.keys()[i]
                    csys = p.compositeLayups[key].orientation.localCsys
                    refAxis = str(p.compositeLayups[key].orientation.axis)
                    addAngle = p.compositeLayups[key].orientation.angle
                    #extract layup axis
                    plyaxis1 = p.datums[csys].axis1

                    #for every ply
                    
                    for j in range(0,nPlies):
                        plyOrientation = p.compositeLayups[key].plies[j].orientation
                        plyAngle = p.compositeLayups[key].plies[j].orientationValue
                        plyOrient = str(p.compositeLayups[key].plies[j].axis)
                        plyOrientType = str(p.compositeLayups[key].plies[j].orientationType)
                        
                        #if no material orientation is defined for each individual py then 
                        #orientatio will be that of the layup
                        if plyOrientation == None:
                            direction1 = plyaxis1.direction

                            totalAngle = addAngle + plyAngle
                            
                            if plyOrientType == 'ANGLE_45':
                                totalAngle = addAngle + 45
                            elif plyOrientType == 'ANGLE_90':
                                totalAngle = addAngle + 90
                            elif plyOrientType == 'ANGLLE_NEG45':
                                totalAngle = addAngle - 45
                                
                                
                            if totalAngle != 0: 
        #                        if refAxis == 'AXIS_1':
        #                    #rotation wrt x axis
        #                            Rx = np.array([[1,0,0],[0,math.cos(totalAngle*math.pi/180),math.sin(totalAngle*math.pi/180)],[0,-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180)]])
        #                            direction2 = np.dot(Rx,plyaxis2.direction)
        #                            direction3 = np.dot(Rx,plyaxis3.direction)
                                   
                                if refAxis == 'AXIS_2':
                                    Ry = np.array([[math.cos(totalAngle*math.pi/180),0,-math.sin(totalAngle*math.pi/180)],[0,1,0],[math.sin(totalAngle*math.pi/180),0,math.cos(45*math.pi/180)]])
                                    direction1 = np.dot(Ry,plyaxis1.direction)
        #                            direction3 = np.dot(Ry,plyaxis3.direction)
                               
                                elif refAxis == 'AXIS_3':
                                    Rz = np.array([[math.cos(totalAngle*math.pi/180),-math.sin(totalAngle*math.pi/180),0],[-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180),0],[0,0,1]])
                                    direction1 = np.dot(Rz,plyaxis1.direction)
        #                            direction2 = np.dot(Rz,plyaxis2.direction)
                        else:
                            #if there is a defines Csys for each ply then change the axis defined
                            plyaxis1 = p.datums[plyOrientation].axis1

                            direction1 = plyaxis1.direction

            #               
                            totalAngle = plyAngle
                            if plyOrientType == 'ANGLE_45':
                                totalAngle = 45
                            elif plyOrientType == 'ANGLE_90':
                                totalAngle = 90
                            elif plyOrientType == 'ANGLLE_NEG45':
                                totalAngle = -45
                                
   
                            if plyOrient == 'AXIS_2':
                                Ry = np.array([[math.cos(totalAngle*math.pi/180),0,-math.sin(totalAngle*math.pi/180)],[0,1,0],[math.sin(totalAngle*math.pi/180),0,math.cos(45*math.pi/180)]])
                                direction1 = np.dot(Ry,plyaxis1.direction)
        #                        direction3 = np.dot(Ry,plyaxis3.direction)
                                
                            elif plyOrient == 'AXIS_3':
                                Rz = np.array([[math.cos(totalAngle*math.pi/180),-math.sin(totalAngle*math.pi/180),0],[-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180),0],[0,0,1]])
                                direction1 = np.dot(Rz,plyaxis1.direction)
        #                        direction2 = np.dot(Rz,plyaxis2.direction)
            
            #            
                        #finding angle between direction of plane axis and fibre axis    
                        

                        doty = np.dot(direction1,perp)
                        producty = np.linalg.norm(direction1)*np.linalg.norm(perp)
        
                        cosy = doty/producty
        
                        theta = 180*math.acos(cosy)/math.pi
                        
                        if theta > 90:
                            theta = 180 - theta
                                
                        #assigning angles to elements
                        plyRegion = p.compositeLayups[key].plies[j].region[0]
                        plyElements = p.sets[plyRegion].elements
                        for k in range(0,len(plyElements)):
                            label = plyElements[k].label
                            if label in area2.keys() and label in materialAngle.keys(): #check if another ply is assigned to the same element
                                materialAngle[label].append(theta) #append to the orientation angle of the element
                                    #assume plys have same relative thickness
                            elif label in area2.keys():
                                materialAngle[label] = [theta]
                                    #if only one ply is assigned to the element
                
                #calculating toughness of composite in relation to fibre orientation
                #Then storing the toughness of composite in a database if existing 
                #orientation already exist
                materialG = {}
                totalU = 0
                materialSigma1 = {}
                totalF = 0
                axis.append(z*axisLength)
                for element in area2.keys():
                    theta = materialAngle[element]
                    #writing  and reading to database
                    materialG[element] = G0 * np.cos(np.multiply(theta,np.pi/180)) + G90 * np.sin(np.multiply(theta,np.pi/180))
                    materialSigma1[element] = sigma1 * np.cos(np.multiply(theta, np.pi/180))
                    #calculate energy dissipation
                    
                    totalU += np.dot(materialG[element],area2[element])
                    totalF += np.dot(materialSigma1[element],area2[element])
                if  totalU not in U.keys():  
                    U[totalU] = [p.features.keys()[-2]]
                    F[totalU] = [totalF]
                    pos[totalU] = [z]
                else:
                    U[totalU].append(p.features.keys()[-2])
                    F[totalU].append(totalF)
                    pos[totalU].append(z)
                print('Total energy dissipated from failure: {0} kJ, Critical force: {1} MN at Point: {2}'.format(float(totalU),float(totalF),pt2))
                plotU.append(totalU)
                plotF.append(totalF)
        critU = min(U.keys())
        planes = U[critU]
        force = F[totalU]
        print('Critical energy is: {0} kJ. Critical Force {1} MN at the following planes: {2}, at Point {3}'.format(critU,force,planes,pt2))
        t2 = time.time()
        print('Run time: {0}'.format(t2-t1))
        
        outputFiles = session.xyPlots.keys()
        fileName = 'Energy Output{0}'.format(fileNumber)
        while fileName in outputFiles:
            fileNumber += 1
            fileName = 'Energy Output{0}'.format(fileNumber)
        xyp = session.XYPlot(fileName)
    
        chartName = xyp.charts.keys()[0]
        chart = xyp.charts[chartName]
        yQuantity = visualization.QuantityType(type=ENERGY)
        c = np.empty((len(axis),2))
        c[:,0] = axis
        c[:,1] = plotU
        xy2 = xyPlot.XYData(data=(c), sourceDescription='Entered from keyboard', 
             axis2QuantityType=yQuantity, )
        c2 = session.Curve(xyData=xy2)
        chart.setValues(curvesToPlot=(c2, ), )
        session.charts[chartName].axes1[0].axisData.setValues(useSystemTitle=False, 
        title='Distance along axis path')
        session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
    
        fileNameF = 'Force Output{0}'.format(fileNumber)
        while fileName in outputFiles:
            fileNumber += 1
            fileNameF = 'Force Output{0}'.format(fileNumber)
        xyp = session.XYPlot(fileNameF)
    
        chartName = xyp.charts.keys()[0]
        chart = xyp.charts[chartName]
        c = np.empty((len(axis),2))
        c[:,0] = axis
        c[:,1] = plotF
        xy2 = xyPlot.XYData(data=(c), sourceDescription='Entered from keyboard', 
            )
        c2 = session.Curve(xyData=xy2)
        chart.setValues(curvesToPlot=(c2, ), )
        session.charts[chartName].axes2[0].axisData.setValues(useSystemTitle=False, 
        title='Force')
        session.charts[chartName].axes1[0].axisData.setValues(useSystemTitle=False, 
        title='Distance along axis path')
        session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
        
        
    if singlePlane == True:
        area2 = {}
        selectPsingle = kwargs.get('selectPsingle')
        if selectPsingle == 'Select points in viewport:':
            coord1 = kwargs.get('coord1')
            coord2 = kwargs.get('coord2')
            coord3 = kwargs.get('coord3')
            if coord1 == None or coord2 == None or coord3 == None:
                raise Exception('Point not picked please pick points')
            s = np.array(coord1.coordinates)
            q = np.array(coord2.coordinates)
            r = np.array(coord3.coordinates)
        elif selectPsingle == 'Define Points Manually:':
            enterCoord1 = kwargs.get('enterCoord1')
            enterCoord2 = kwargs.get('enterCoord2')
            enterCoord3 = kwargs.get('enterCoord3')
            if enterCoord1 == () or enterCoord2 == () or enterCoord3 == ():
                raise Exception('Points not defined please define point')
            elif len(enterCoord1[0]) != 3 or len(enterCoord2[0]) != 3 or len(enterCoord3[0]) != 3:
                raise Exception('Axis points not defined properly')
            s = np.array([enterCoord1[0][0],enterCoord1[0][1],enterCoord1[0][2]])
            q = np.array([enterCoord2[0][0],enterCoord2[0][1],enterCoord2[0][2]])
            r = np.array([enterCoord3[0][0],enterCoord3[0][1],enterCoord3[0][2]])
            
        lambRecord = []
            
        p.DatumPlaneByThreePoints(point1=s, 
            point2=q, 
            point3=r) 
        
        pq = q - s
        pr = r - s
        
        perp = np.cross(pq,pr)
        d = -np.dot(perp,-s)
        #select elements
        vec = np.zeros((12,3))
        vertexply = np.zeros((8,3))
        vertexplyi = np.zeros((12,3))
        vertexplyj = np.zeros((12,3))
        vertexsubply = np.zeros((8,3))
        #finding area of intercepts in the general case of any plane interception on the element
        for e in range(0,len(elementSet)):
            connected = []
                
            connected = elementSet[e].connectivity
            label = elementSet[e].label
            area2[label] = []
                # sort node index to keep track of lines and vertex
            sort = np.sort(connected)
            side1 = np.array(nodeSet[sort[0]].coordinates)-np.array(nodeSet[sort[1]].coordinates)
            side2 = np.array(nodeSet[sort[0]].coordinates)-np.array(nodeSet[sort[2]].coordinates)
            side3 = np.array(nodeSet[sort[0]].coordinates)-np.array(nodeSet[sort[4]].coordinates)
            size = [np.linalg.norm(side1)]
            size.append(np.linalg.norm(side2))
            size.append(np.linalg.norm(side3))
            maxSize = max(size)
            diag = math.sqrt(maxSize**2+maxSize**2)
            # calculate the elements near the plane taking account element mesh size
            if float(abs(np.dot(perp,np.array(nodeSet[sort[0]].coordinates)) - d)/(np.linalg.norm(perp))) <= float(math.sqrt(diag**2+maxSize**2)):
                region = []
                nLayups = len(p.compositeLayups.keys())
                if nLayups == 0:
                        sys.exit('Composite layup not defined please define layup before using software')

                # number of layups
                for n in range(0,nLayups):
                    key = p.compositeLayups.keys()[n]
                    nPlies = len(p.compositeLayups[key].plies)
                    #ply stack direction
                    stackAxis = p.compositeLayups[key].orientation.stackDirection
                    plyThickness = []
                    # calculate of the total relative thickness of plies add to one otherwise 
                    #terminate the program
                    for j in range(0,nPlies):
                        plyThickness.append(p.compositeLayups[key].plies[j].thickness)
                        region.append(p.compositeLayups[key].plies[j].region[0])
                    if sum(plyThickness) != 1 and region != [region[0]] * nPlies:
                        sys.exit('Defined relative thiknesses do not add up to one and/or regions defined for compositelayup are different')
                    plyRelThickness = 0  
                    # Three possible directions for the ply stack
                    # Vertex change depending on stack direction
                    if str(stackAxis) == 'STACK_1':
                        if np.array(nodeSet[sort[4]].coordinates)[0] < np.array(nodeSet[sort[0]].coordinates)[0]:
                            vertexsubply[0,:] = np.array(nodeSet[sort[4]].coordinates)
                            vertexsubply[1,:] = np.array(nodeSet[sort[5]].coordinates)
                            vertexsubply[2,:] = np.array(nodeSet[sort[6]].coordinates)
                            vertexsubply[3,:] = np.array(nodeSet[sort[7]].coordinates)
                        
                            vertexsubply[4,:] = np.array(nodeSet[sort[0]].coordinates)
                            vertexsubply[5,:] = np.array(nodeSet[sort[1]].coordinates)
                            vertexsubply[6,:] = np.array(nodeSet[sort[2]].coordinates)
                            vertexsubply[7,:] = np.array(nodeSet[sort[3]].coordinates)
                        else:
                            vertexsubply[4,:] = np.array(nodeSet[sort[4]].coordinates)
                            vertexsubply[5,:] = np.array(nodeSet[sort[5]].coordinates)
                            vertexsubply[6,:] = np.array(nodeSet[sort[6]].coordinates)
                            vertexsubply[7,:] = np.array(nodeSet[sort[7]].coordinates)
                        
                            vertexsubply[0,:] = np.array(nodeSet[sort[0]].coordinates)
                            vertexsubply[1,:] = np.array(nodeSet[sort[1]].coordinates)
                            vertexsubply[2,:] = np.array(nodeSet[sort[2]].coordinates)
                            vertexsubply[3,:] = np.array(nodeSet[sort[3]].coordinates)
                    elif str(stackAxis) == 'STACK_2':
                        if np.array(nodeSet[sort[1]].coordinates)[1] < np.array(nodeSet[sort[0]].coordinates)[1]:
                            vertexsubply[0,:] = np.array(nodeSet[sort[1]].coordinates)
                            vertexsubply[1,:] = np.array(nodeSet[sort[3]].coordinates)
                            vertexsubply[2,:] = np.array(nodeSet[sort[5]].coordinates)
                            vertexsubply[3,:] = np.array(nodeSet[sort[7]].coordinates)
                            
                            vertexsubply[4,:] = np.array(nodeSet[sort[0]].coordinates)
                            vertexsubply[5,:] = np.array(nodeSet[sort[2]].coordinates)
                            vertexsubply[6,:] = np.array(nodeSet[sort[4]].coordinates)
                            vertexsubply[7,:] = np.array(nodeSet[sort[6]].coordinates)
                        else:
                            vertexsubply[4,:] = np.array(nodeSet[sort[1]].coordinates)
                            vertexsubply[5,:] = np.array(nodeSet[sort[3]].coordinates)
                            vertexsubply[6,:] = np.array(nodeSet[sort[5]].coordinates)
                            vertexsubply[7,:] = np.array(nodeSet[sort[7]].coordinates)
                            
                            vertexsubply[0,:] = np.array(nodeSet[sort[0]].coordinates)
                            vertexsubply[1,:] = np.array(nodeSet[sort[2]].coordinates)
                            vertexsubply[2,:] = np.array(nodeSet[sort[4]].coordinates)
                            vertexsubply[3,:] = np.array(nodeSet[sort[6]].coordinates)
                    elif str(stackAxis) == 'STACK_3':
                        if np.array(nodeSet[sort[2]].coordinates)[2] < np.array(nodeSet[sort[0]].coordinates)[2]:
                            vertexsubply[0,:] = np.array(nodeSet[sort[2]].coordinates)
                            vertexsubply[1,:] = np.array(nodeSet[sort[3]].coordinates)
                            vertexsubply[2,:] = np.array(nodeSet[sort[6]].coordinates)
                            vertexsubply[3,:] = np.array(nodeSet[sort[7]].coordinates)
                            
                            vertexsubply[4,:] = np.array(nodeSet[sort[0]].coordinates)
                            vertexsubply[5,:] = np.array(nodeSet[sort[1]].coordinates)
                            vertexsubply[6,:] = np.array(nodeSet[sort[4]].coordinates)
                            vertexsubply[7,:] = np.array(nodeSet[sort[5]].coordinates)
                        
                        else:
                            vertexsubply[4,:] = np.array(nodeSet[sort[2]].coordinates)
                            vertexsubply[5,:] = np.array(nodeSet[sort[3]].coordinates)
                            vertexsubply[6,:] = np.array(nodeSet[sort[6]].coordinates)
                            vertexsubply[7,:] = np.array(nodeSet[sort[7]].coordinates)
                            
                            vertexsubply[0,:] = np.array(nodeSet[sort[0]].coordinates)
                            vertexsubply[1,:] = np.array(nodeSet[sort[1]].coordinates)
                            vertexsubply[2,:] = np.array(nodeSet[sort[4]].coordinates)
                            vertexsubply[3,:] = np.array(nodeSet[sort[5]].coordinates)
                   # assign vertex to different variable for further calculation
                    for i in range(0,4):     
                            vertexply[i,:] = vertexsubply[i,:]
                    for j in range(0,nPlies):
                        plyRelThickness += p.compositeLayups[key].plies[j].thickness
                        lamb = []
                        
                        #assign ply coordinates according to ply thickness using equation of a line 
                        for i in range(0,4):
                            vertexply[i+4,:] = vertexsubply[i,:] + plyRelThickness * (vertexsubply[i+4,:] - vertexsubply[i,:])
                            vertexplyi[i,:] = vertexply[i,:]
                            vertexplyj[i,:] = vertexply[i+4,:]
                            
                        vertexplyi[4,:] = vertexply[0,:] 
                        vertexplyj[4,:] = vertexply[1,:]
                        
                        vertexplyi[5,:] = vertexply[3,:]
                        vertexplyj[5,:] = vertexply[1,:]
                        
                        vertexplyi[6,:] = vertexply[2,:]
                        vertexplyj[6,:] = vertexply[0,:]
                        
                        vertexplyi[7,:] = vertexply[2,:] 
                        vertexplyj[7,:] = vertexply[3,:] 
                        
                        vertexplyi[8,:] = vertexply[4,:]
                        vertexplyj[8,:] = vertexply[5,:]
                        
                        vertexplyi[9,:] = vertexply[7,:]
                        vertexplyj[9,:] = vertexply[5,:]
                        
                        vertexplyi[10,:] = vertexply[6,:]
                        vertexplyj[10,:] = vertexply[4,:]
                        
                        vertexplyi[11,:] = vertexply[6,:]
                        vertexplyj[11,:] = vertexply[7,:]
                        
        #                # Algorithm for lines and finding lambda (gradient of intercept)
        #                # refer to paper
                        
                        for i in range(0,len(vec)):
                            vec[i,:] = vertexplyi[i,:] - vertexplyj[i,:]
                            if np.dot(perp, vec[i,:]) != 0:
                                lamb.append(round((d - np.dot(perp, vertexplyj[i,:]))/(np.dot(perp, vec[i,:])),3)) 
                            else:
                                lamb.append(-100)
                        k = 0  
                        intercept = np.zeros((12,3))
                        count = Counter(lamb)
            # checking for every value of lambda if there is an intercept
                        for i in range(0,len(lamb)):
                            if 0 < lamb[i] < 1 and count[0] >= 4 or count[1] >= 4:
                                continue
                            elif lamb[i] < 1 and lamb[i] >= 0:
                                #if there is an intercept find the coordinates using parametric method
                                temp = (np.multiply(vec[i,:],lamb[i]) + vertexplyj[i,:])
                                intercept[k,:] = temp
                                #if lambda = 1 check if there are over lapping intersections
                                if lamb[i] == 1:
                                    for j in range(0,len(intercept)):
                                        if j < k and np.array_equal(intercept[j,:], intercept[k,:]):
                                            k -= 1
                                k += 1
                            # check if intersection is at the nodes and if the itersections lie
                            # diagonally directly at the nodes so we have lambda = 1 and 0 
                            # the same number of times
                            elif lamb[i] == 1 and count[0] == count[1]:
                                temp = (np.multiply(vec[i,:],lamb[i]) + vertexplyj[i,:])
                                intercept[k,:] = temp
                                if lamb[i] == 0:
                                    for j in range(0,len(intercept)):
                                        if j < k and np.array_equal(intercept[j,:], intercept[k,:]):
                                            k -= 1
                                k += 1
                                
                            if abs(lamb[i]) <=1 and lamb[i] >= 0:
                                lambRecord.append(lamb[i])
                        #find the areas by sorting the coordinates of intersections
                        norm = 0
                        angles = {}

                        #for just three points of intersection find area of triangle
                        if k == 3:
                            pq = intercept[0,:] - intercept[1,:]
                            pr = intercept[0,:] - intercept[2,:]
                            cros = np.cross(pq,pr)
                            norm += np.linalg.norm(cros)/2
                            area2[label].append(norm)
                        #things become harder with more than 3 points
                        elif 4 <= k <= 6:
                            #define reference vector
                            ref = intercept[0,:] - intercept[1,:]
                            #set the reference angle of the reference vector to zero
                            angles[1] = 0
                            for i in range(2,k):
                                #find all other vectors with respect to the first point
                                #find their angles with respect to the reference angle
                                #then arrange from smallest to largest angle
                                ref2 = intercept[0,:] - intercept[i,:]
                                cros = np.cross(ref,ref2)
                                refNorm = np.linalg.norm(cros)
                                product = np.linalg.norm(ref)*np.linalg.norm(ref2)
                                dot = np.dot(ref,ref2)
                                #if statement added as to test which side the point lies
                                #on reference vector. If the z componenet of the cross
                                #product of the reference vector and the new vector
                                #is negative then the point lies on one side, if positive
                                #point lies on other side
                                if cros[2] < 0:
                                    sin = round(-refNorm/product,5)
                                else:
                                    sin = round(refNorm/product,5)
                                cos = round(dot/product,5)
                                #between 0 and 90
                                if sin >= 0 and cos >= 0:
                                    theta = (180*math.acos(cos))/math.pi
                                    angles[i] = theta
                                #between -90 and 0
                                elif sin <= 0 and cos >= 0:
                                    theta = (180*math.asin(sin))/math.pi
                                    angles[i] = theta
                                #between 90 and 180
                                elif sin >= 0 and cos <= 0:
                                    theta = 180-((180*math.asin(sin))/math.pi)
                                    angles[i] = theta
                                #between -180 and -90 
                                elif sin <= 0 and cos <= 0:
                                    theta = -180 + ((180*math.asin(sin))/math.pi)
                                    angles[i] = theta
                            sorted_angles = sorted(angles.items(), key=lambda x: x[1])
                            sorted_index = np.array(sorted_angles)[:,0]
                            for i in range(0,len(sorted_index)-1):
                                pq = intercept[0,:] - intercept[sorted_index[i],:]
                                pr = intercept[0,:] - intercept[sorted_index[i+1],:]
                                cros = np.cross(pq,pr)
                                norm += np.linalg.norm(cros)/2
                            area2[label].append(norm)
                        else:
                            area2[label].append(0)
                        
                        #assign coordinates of last four vertex to first four vertex
                        #in order to calcualte area of the new ply
                        for i in range(0,4):
                            vertexply[i,:] = vertexply[i+4,:]
                       #5.7 seconds when added to only consider elements near the plane and 28.4 seconds without modifycation 
        if lambRecord[0:] == [1] * len(lambRecord) or lambRecord[0:] == [0] * len(lambRecord):
                    area2 = {}
        for label in area2.keys():
            if area2[label] == [0] * nPlies or area2[label] == []:
                del area2[label]

        if area2 != {}:
        # extracting material orientation
            Name = 'csys-plane'
            if Name in p.features.keys():
                del p.features[Name]
            #creating coordinate system for plane
            p.DatumCsysByThreePoints( origin = s, name = Name,coordSysType = CARTESIAN,point1 = q, point2 = s + perp)
           #extracting coordinate axis
            planeAxis1 = p.datums[p.features[Name].id].axis1
            planeAxis2 = p.datums[p.features[Name].id].axis2
            planeAxis3 = p.datums[p.features[Name].id].axis3
            
            
        
    #        Rx = np.zeros((3,3))
            Ry = np.zeros((3,3))
            Rz = np.zeros((3,3))
            materialAngle = {}
            #for each defined ply layup
            for i in range(0,nLayups):
                key = p.compositeLayups.keys()[i]
                csys = p.compositeLayups[key].orientation.localCsys
                refAxis = str(p.compositeLayups[key].orientation.axis)
                addAngle = p.compositeLayups[key].orientation.angle
                #extract layup axis
                plyaxis1 = p.datums[csys].axis1

                #for every ply
                
                for j in range(0,nPlies):
                    plyOrientation = p.compositeLayups[key].plies[j].orientation
                    plyAngle = p.compositeLayups[key].plies[j].orientationValue
                    plyOrient = str(p.compositeLayups[key].plies[j].axis)
                    plyOrientType = str(p.compositeLayups[key].plies[j].orientationType)
                    
                    #if no material orientation is defined for each individual py then 
                    #orientatio will be that of the layup
                    if plyOrientation == None:
                        direction1 = plyaxis1.direction

                        totalAngle = addAngle + plyAngle
                        
                        if plyOrientType == 'ANGLE_45':
                            totalAngle = addAngle + 45
                        elif plyOrientType == 'ANGLE_90':
                            totalAngle = addAngle + 90
                        elif plyOrientType == 'ANGLLE_NEG45':
                            totalAngle = addAngle - 45
                            
                            
                        if totalAngle != 0: 

                            if refAxis == 'AXIS_2':
                                Ry = np.array([[math.cos(totalAngle*math.pi/180),0,-math.sin(totalAngle*math.pi/180)],[0,1,0],[math.sin(totalAngle*math.pi/180),0,math.cos(45*math.pi/180)]])
                                direction1 = np.dot(Ry,plyaxis1.direction)
    #                            direction3 = np.dot(Ry,plyaxis3.direction)
                           
                            elif refAxis == 'AXIS_3':
                                Rz = np.array([[math.cos(totalAngle*math.pi/180),-math.sin(totalAngle*math.pi/180),0],[-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180),0],[0,0,1]])
                                direction1 = np.dot(Rz,plyaxis1.direction)
    #                            direction2 = np.dot(Rz,plyaxis2.direction)
                    else:
                        #if there is a defines Csys for each ply then change the axis defined
                        plyaxis1 = p.datums[plyOrientation].axis1
  
                        direction1 = plyaxis1.direction

        #               
                        totalAngle = plyAngle
                        if plyOrientType == 'ANGLE_45':
                            totalAngle = 45
                        elif plyOrientType == 'ANGLE_90':
                            totalAngle = 90
                        elif plyOrientType == 'ANGLLE_NEG45':
                            totalAngle = -45
                            
                        if plyOrient == 'AXIS_2':
                            Ry = np.array([[math.cos(totalAngle*math.pi/180),0,-math.sin(totalAngle*math.pi/180)],[0,1,0],[math.sin(totalAngle*math.pi/180),0,math.cos(45*math.pi/180)]])
                            direction1 = np.dot(Ry,plyaxis1.direction)
    #                        direction3 = np.dot(Ry,plyaxis3.direction)
                            
                        elif plyOrient == 'AXIS_3':
                            Rz = np.array([[math.cos(totalAngle*math.pi/180),-math.sin(totalAngle*math.pi/180),0],[-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180),0],[0,0,1]])
                            direction1 = np.dot(Rz,plyaxis1.direction)
        
        #            
                    #finding angle between direction of plane axis and fibre axis    
                    
                    doty = np.dot(direction1,planeAxis2.direction)
                    producty = np.linalg.norm(direction1)*np.linalg.norm(planeAxis2.direction)
    
    
                    cosy = doty/producty
    
                    
                    theta = 180*math.acos(cosy)/math.pi
                    
    
                    if theta > 90:
                        theta = 180 - theta
                            
                    #assigning angles to elements
                    plyRegion = p.compositeLayups[key].plies[j].region[0]
                    plyElements = p.sets[plyRegion].elements
                    for k in range(0,len(plyElements)):
                        label = plyElements[k].label
                        if label in area2.keys() and label in materialAngle.keys(): #check if another ply is assigned to the same element
                            materialAngle[label].append(theta) #append to the orientation angle of the element
                                #assume plys have same relative thickness
                        elif label in area2.keys():
                            materialAngle[label] = [theta]
                                #if only one ply is assigned to the element
            
 
            #calculating toughness of composite in relation to fibre orientation
            #Then storing the toughness of composite in a database if existing 
            #orientation already exist
            materialG = {}
            totalU = 0
            totalF = 0
            materialSigma1 = {}
            for element in area2.keys():
                theta = materialAngle[element]
                #writing  and reading to database
                materialG[element] = G0 * np.cos(np.multiply(theta,np.pi/180)) + G90 * np.sin(np.multiply(theta,np.pi/180))
                materialSigma1[element] = sigma1 * np.cos(np.multiply(theta,np.pi/180))
      
                #calculate energy dissipation
                
                totalU += np.dot(materialG[element],area2[element])
                totalF += np.dot(materialSigma1[element],area2[element])    
                
            print('Total energy dissipated from failure: {0} kJ. Critical Force: {1} MN for plane selected at: {2}, {3}, {4}'.format(float(totalU),float(totalF), s,q,r))
            t2 = time.time()
            print('Run time: {0}'.format(t2-t1))
    if multiPlane == False and singlePlane == False:
        raise Exception('No checkboxes were ticked please tick at least one checkbox')