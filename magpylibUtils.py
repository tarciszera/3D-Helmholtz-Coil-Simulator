# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 11:01:56 2019

@author: Tarcis feat. Lucas

Helmholtz project auxiliar library

air permeability = 1.25663753e10−6

"""
import numpy as np
import magpylib as magpy
import matplotlib.pyplot as plt
from numpy import array,linspace,arange,sin,cos,tan,sqrt,r_,c_,pi,exp,sign,amax,amin,arctan,meshgrid,mean,polyfit,poly1d,std,dot,log,log10,logspace,arccos
#import gxPackage_v03 as gx
#import plotConfig_v18 as pconfig 
from numpy.linalg import norm
import datetime
import time

#%% coil function
def coil(color = 'k', I = 1, coilIntDiameter = 1, referencePosition=[0,0,0], wireDiameter=1, 
         loopsInEachEvenLayer=1, evenLayers = 1, oddLayers = 1, collectionToAdd = magpy.Collection()):
    """ 
    Return the model of a rolled coil around the 'z' axis.
    
    Parameters
    ----------
    I : float
        The Coil Current in Amperes.
    coilIntDiameter : {non-zero int}
        The Diameter of the coil's inner loops, in mm.
    referencePosition : List, default [0,0,0]
        A List containing 3 coordinates.
    wireDiameter : float, default 1
        Diameter of the coil wire in mm.
    loopsInEachEvenLayer : non-zero integer, default 1 
        How many loops in a single even layer
    evenLayers : integer, default 1
        Integer number of even layers.
    oddLayers : integer, default 1
        Integer number of odd layers.
    collectionToAdd : Magpy Collection, default empty
        Append coil source to a magpy collection.
    
    Example
    ---------
    test = coil(I = 4, coilIntDiameter = 77.26 , referencePosition=[0,0,0], wireDiameter=1.59, 
                loopsInEachEvenLayer=3, evenLayers = 3, oddLayers = 3, collectionToAdd = magpy.Collection())
    test.displaySystem()
    
    """
    
     
        
    # Initialization for the loop
    even = True                                                     # Flag if the coil rolling is in even format
    singleCoilWidth = wireDiameter*loopsInEachEvenLayer   
    centerPosition = referencePosition[2]                                    # The center position of the coil in Z
    initPosEven = centerPosition - singleCoilWidth/2 + wireDiameter/2     # Initial position in the 'even' form of rolling the coil
    initPosOdd = centerPosition - singleCoilWidth/2 + wireDiameter        # Initial position in the 'odd' form of rolling the coil
    stepPos = wireDiameter                                          # The position increment
    finalPosEven = centerPosition + singleCoilWidth/2 - wireDiameter/2 + stepPos    # Final position in the 'even' form of rolling the coil
    finalPosOdd = centerPosition + singleCoilWidth/2 - wireDiameter + stepPos     # Final position in the 'odd' form of rolling the coil
    stepDiameter = 2*(wireDiameter**2 - (wireDiameter/2)**2)**0.5   # The diameter increment
    subLoops = np.int64(singleCoilWidth/wireDiameter)                     # The max value of the loops in the 'even' rolling, the odd max value is (subLoops - 1)
    coilIntDiameter += wireDiameter
    
    if loopsInEachEvenLayer == 1:
        nLoops = evenLayers + oddLayers
    else:
        nLoops = loopsInEachEvenLayer*evenLayers + (loopsInEachEvenLayer-1)*oddLayers        # Total number of loops in each coil   
    # ----------------------------
    
    
    # Starting to roll the coil
    while nLoops > 0:                                   # loop that counts the number of wireLoops already done
        if (even == 1) | (singleCoilWidth == wireDiameter):                                                    # The statement that says if the rolling is in the 'even' or 'odd' form
            zs = np.arange(initPosEven, finalPosEven, stepPos)      # The 'zs' is the variable that will be used for rolling loop position
            for n in range(0,subLoops):                             
                coilLoop = magpy.source.current.Circular(curr=I, dim=coilIntDiameter,pos=[0,0,zs[n]])
                collectionToAdd.addSources(coilLoop)
                nLoops = nLoops - 1
                if nLoops <= 0:
                    break
            even = False
            coilIntDiameter = coilIntDiameter + stepDiameter
            
        else:                                              # The statement that says if the rolling is in the 'even' or 'odd' form
            zs = np.arange(initPosOdd, finalPosOdd, stepPos)      # The 'zs' is the variable that will be used for rolling loop position
            for n in range(0,subLoops-1):
                coilLoop = magpy.source.current.Circular(curr=I, dim=coilIntDiameter,pos=[0,0,zs[n]])            
                collectionToAdd.addSources(coilLoop)
                nLoops = nLoops - 1
                if nLoops <= 0:
                    break
            even = True
            coilIntDiameter = coilIntDiameter + stepDiameter
            
    return collectionToAdd
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
#%% ConstrucAcess function 
def constructAccess( I = 4, centralPointBetweenHelmCoils = [0,0,0], coilExtDiameter = 1, relativeDistance = 1, 
                    wireDiameter = 1, singleCoilWidth = 1, lengthBetweenCoils = 1, helmCoil = magpy.Collection()):
    
    """
    Returns a collection with the wires access in 'x' axis orientation
    to the coils of a helmholtz collection with 'z' field orientation
    
    Parameters
    ----------
    I : float
        The Coil Current in Amperes.
    centralPointBetweenHelmCoils : integer, default 1
        Integer number of odd layers.
    coilExtDiameter : float
        The external diameter of one of the helmholtz coils, in mm.
    relativeDistance : float 0~1, default 1
        The porcentage of the 1.5*coilExtDiameter distance from the 
        centralPointBetweenHelmCoils that the wire will split to the coils
    wireDiameter : float, default 1
        Diameter of the coil wire in mm.
    singleCoilWidth : float, default 1 
        Width of a single coil of helmholtz coil, in mm
    lengthBetweenCoils : float, default 1
        The internal length between the two coils of the helmholtz coil
    helmCoil : Magpy Collection, default empty
        Append the wires to a magpy collection.
    
    Example
    ---------
    test = constructAccess( I = 4, centralPointBetweenHelmCoils = [0,0,0], coilExtDiameter = 1, relativeDistance = 0.5, 
                    wireDiameter = 1, singleCoilWidth = 1, lengthBetweenCoils = 1, helmCoil = magpy.Collection()):
    
    test.displaySystem()
    """
    
    # distance from the centralPointBetweenHelmCoils to the center of one coil of Helmholtz coil
    halfCoillengthBetweenCoils = (lengthBetweenCoils + singleCoilWidth)/2
                                 
    # Initializing the statingPoint of the wires with the central point of the coil + some distance relative to the coilIntDiameter
    wiresStartingPoint = centralPointBetweenHelmCoils                                        
    wiresStartingPoint[0] += 1.5*coilExtDiameter   
    
    #  |
    #
    #
    vector = [[wiresStartingPoint[0],wiresStartingPoint[1], wireDiameter +  wiresStartingPoint[2]], [wiresStartingPoint[0]*relativeDistance,0, wireDiameter] ]      
    line = magpy.source.current.Line(curr=I,vertices=vector)
    helmCoil.addSources(line)
    
    #  _|
    # 
    # 
    vector = [vector[1], [vector[1][0], vector[1][1], vector[1][2] + halfCoillengthBetweenCoils]]
    line = magpy.source.current.Line(curr=I,vertices=vector)
    helmCoil.addSources(line)
    
    #  _|
    # |
    # |
    vector = [vector[1], [coilExtDiameter/2 - wireDiameter, vector[1][1], vector[1][2]]]
    line = magpy.source.current.Line(curr=I,vertices=vector)
    helmCoil.addSources(line)
    
    #  _|
    # | 
    # ||
    vector = [[vector[1][0], vector[1][1], vector[1][2] - wireDiameter], [vector[0][0] -wireDiameter, vector[0][1], vector[0][2] - wireDiameter]]
    line = magpy.source.current.Line(curr=I,vertices=vector)
    helmCoil.addSources(line)
    
    # going to the other side
    #  _|
    # | __
    # ||
    vector = [vector[1], [vector[1][0], vector[1][1], -vector[1][2]]]
    line = magpy.source.current.Line(curr=I,vertices=vector)
    helmCoil.addSources(line)
    
    #  _|
    # | __
    # || |
    vector = [vector[1], [coilExtDiameter/2 - wireDiameter, vector[1][1], vector[1][2]]]
    line = magpy.source.current.Line(curr=I,vertices=vector)
    helmCoil.addSources(line)
    
    #  _|
    # | __ |
    # || ||
    vector = [[vector[1][0], vector[1][1], vector[1][2] - wireDiameter], [vector[0][0] + wireDiameter, vector[0][1], vector[0][2] - wireDiameter]]
    line = magpy.source.current.Line(curr=I,vertices=vector)
    helmCoil.addSources(line)
    
    #  _|  _
    # | __ |
    # || ||
    vector = [vector[1], [wiresStartingPoint[0]*relativeDistance,0, -wireDiameter]]
    line = magpy.source.current.Line(curr=I,vertices=vector)
    helmCoil.addSources(line)
    
    #  _| |_
    # | __ |
    # || ||
    vector = [[wiresStartingPoint[0]*relativeDistance,0, -wireDiameter], [wiresStartingPoint[0],wiresStartingPoint[1], -wireDiameter +  wiresStartingPoint[2]]]
    line = magpy.source.current.Line(curr=I,vertices=vector)
    helmCoil.addSources(line)
    
    return helmCoil

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
#%% helmholtzCoil function 
def helmholtzCoil(color = 'k', collectionToAdd = 'NULL', I = 1, coilIntDiameter = 1, wireDiameter=1, loopsInEachEvenLayer=1, evenLayers=1, oddLayers=1, 
                  coilsDistance = 1, orientation = 'z', relativeDistance = 0.5,):
    ### (_TO DO: make the helmholtz coil do not touch each other_)
    """ 
    Generate and add 2 coils to a collection , with r distance (in mm) between both coils.
    Return a colection with the helmholtz coil inside
    
    Parameters
    ----------
    I : float
        The Helmhotz Coil Current in Amperes.
    coilIntDiameter : float
        The internal diameter of the coil in mm
    wireDiameter : float, optimizationional
        Diameter of the coil wire in mm.
    loopsInEachEvenLayer : non-zero integer, default 1 
        How many loops in a single even layer
    evenLayers : integer, default 1
        Integer number of even layers.
    oddLayers : integer, default 1
        Integer number of odd layers.
    collectionToAdd : Magpy Collection, default empty
        Append coil source to a magpy collection.
    coilsDistance : float , default 1
        Distance between the coil of the helmholtz coil
    orientation : {'x', 'y', 'z'}, optimizationional
        Change the orientation of the resulting coil having its center parallel with Plane X, Y, or Z.
    
    *OBS: the odd layers always must be evenLayers -1 or = evenLayers 
    
    Example
    ---------
    test = helmholtzCoil(I = 4, coilIntDiameter = 77.26, wireDiameter=1.59, 
                     loopsInEachEvenLayer=3, evenLayers=3, oddLayers=3, 
                     collectionToAdd=magpy.Collection(), coilsDistance = 77.26/2, orientation = 'x')
    test.displaySystem()
    """
    # seting some values to do the calculations after that
    if collectionToAdd == 'NULL':
        collectionToAdd=magpy.Collection()
        
    helmoltzRotationCenterXYZ = [0,0,0]
    singleCoilWidth = wireDiameter*loopsInEachEvenLayer 
    
    # doing the inverse rotations to change the helmhotz coil orientation (because this function only create the coil in 'z' orientation)
    if orientation == 'y':
        collectionToAdd.rotate(90,[1,0,0],anchor=helmoltzRotationCenterXYZ)
        collectionToAdd.rotate(90,[0,0,1],anchor=helmoltzRotationCenterXYZ)
    elif orientation == 'x':
        collectionToAdd.rotate(-90,[0,0,1],anchor=helmoltzRotationCenterXYZ)
        collectionToAdd.rotate(-90,[1,0,0],anchor=helmoltzRotationCenterXYZ)
    # ---------------------------------------------------------------------
    
    
    # construct the two coils of a helmholtz coil ----------------------------------------------------------------------------
    collectionToAdd = coil(I = I, color = color, coilIntDiameter = coilIntDiameter , referencePosition = [0,0,-coilsDistance/2-singleCoilWidth/2], 
                           wireDiameter=wireDiameter, loopsInEachEvenLayer=loopsInEachEvenLayer, 
                           evenLayers = evenLayers, oddLayers = oddLayers, collectionToAdd = collectionToAdd)
    collectionToAdd = coil(I = I, color = color, coilIntDiameter = coilIntDiameter , referencePosition = [0,0, coilsDistance/2+singleCoilWidth/2], 
                           wireDiameter=wireDiameter, loopsInEachEvenLayer=loopsInEachEvenLayer, 
                           evenLayers = evenLayers, oddLayers = oddLayers, collectionToAdd = collectionToAdd)
    # ------------------------------------------------------------------------------------------------------------------------
    '''
    # construct the access of supply wires to the two coils of a helmholtz coil ----------------------------------------------
    stepRadius = (wireDiameter**2 - (wireDiameter/2)**2)**0.5
    coilExtDiameter = coilIntDiameter + (wireDiameter + (evenLayers + oddLayers -1)*stepRadius)*2    # calculate the external diameter of the coil
    collectionToAdd = constructAccess(I = I, coilExtDiameter = coilExtDiameter, wireDiameter = wireDiameter, 
                                      centralPointBetweenHelmCoils = [0,0,0], lengthBetweenCoils = coilsDistance, 
                                      singleCoilWidth = singleCoilWidth, helmCoil = collectionToAdd, 
                                      relativeDistance = relativeDistance)
    # ------------------------------------------------------------------------------------------------------------------------
    '''
    
    # doing the rotations to change the helmhotz coil orientation
    if orientation == 'y':
        collectionToAdd.rotate(-90,[0,0,1],anchor=helmoltzRotationCenterXYZ)
        collectionToAdd.rotate(-90,[1,0,0],anchor=helmoltzRotationCenterXYZ)
    elif orientation == 'x':
        collectionToAdd.rotate(90,[1,0,0],anchor=helmoltzRotationCenterXYZ)
        collectionToAdd.rotate(90,[0,0,1],anchor=helmoltzRotationCenterXYZ)
    # ---------------------------------------------------------------------
        
    return collectionToAdd

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
#%%
def plotBxyz(collectionToPlot, collectionToCompare = magpy.Collection(), plotBounds = [0,10,0,10], orientation = 'z', orderMagnitude = 'mT',
          fieldDif = False, figureSize = [10,8], nPlotPoints = 40, xyz0 = 0, figureTittle = 'Figure', compareToCenter = False, centerToCompare = [0,0,0]):
    
    """ 
    Plot in a 2d (x,z) in a y0 = constant graphic with 2 screens the [module %, angles º(related to the z direction)] 
    diference of the B in [0,0,0] and the other points of the graphic
    
    Parameters
    ----------
    collectionToPlot : Magpy Collection
        The collection used to caculate the plots
    plotBounds : float list [4]
        its the initial and final X and Z boundaries
        plotBound[0] = initialX, plotBound[1] = finalX
        plotBound[2] = initialZ, plotBound[3] = finalZ
    percent : bool
        If True, plot the percent diference relative to the B = [0,0,0] field
        If False, plot the true field in every point of the graphic
    optimization : bool
        If True, plot based in an eliptical line in the screen
        If False, plot every point related to the field in each point
    figureSize : list(integers) [2]
        Define the size of the entire figure
    nPlotPoints : integer
        Number of point that will be used to plot the figure in each axis
    y0 : float
        The plot made is an 2d 'z' by 'x', this select where in 'y' the plot will be made 
    figureTittle : string
        Define the entire figure tittle
        
    Example
    ---------
    test = helmholtzCoil(I = 4, coilIntDiameter = 77.26, wireDiameter=1.59, 
                         loopsInEachEvenLayer=3, evenLayers=3, oddLayers=3, 
                         collectionToAdd=magpy.Collection(), coilsDistance = 77.26/2, orientation = 'z')
    plotB(collectionToPlot = test, optimization= True, figureTittle = 'Test figure', nPlotPoints = 10, percent = True)
        
    test.displaySystem()
    """
    if orderMagnitude == 'uT':
        multiplier = 1000
    else:
        multiplier = 1
    # Here creates the screens and name the axis
    fig = plt.figure(figsize=(figureSize[0], figureSize[1]), dpi=80)
    fig.suptitle(figureTittle, fontsize=16)
    if compareToCenter:
        AXS = [fig.add_subplot(2,2,i, axisbelow=True) for i in range(1,5)]
        subPlot1,subPlot2,subPlot3,subPlot4 = AXS
    else:
        subPlot1 = fig.add_subplot(1,1,1, axisbelow=True)
        
    
    if orientation == 'y':
        subPlot1.set_xlabel('x axis (mm)')
        subPlot1.set_ylabel('z axis (mm)')
        
    elif orientation == 'x':
        subPlot1.set_xlabel('z axis (mm)')
        subPlot1.set_ylabel('y axis (mm)')
        collectionToPlot.rotate(90,[0,1,0],anchor=[0,0,0])
        collectionToCompare.rotate(90,[0,1,0],anchor=[0,0,0])
        collectionToPlot.rotate(90,[1,0,0],anchor=[0,0,0])
        collectionToCompare.rotate(90,[1,0,0],anchor=[0,0,0])
        
    elif orientation == 'z':
        subPlot1.set_xlabel('y axis (mm)')
        subPlot1.set_ylabel('x axis (mm)')
        collectionToPlot.rotate(-90,[0,1,0],anchor=[0,0,0])
        collectionToCompare.rotate(-90,[0,1,0],anchor=[0,0,0])
        collectionToPlot.rotate(-90,[0,0,1],anchor=[0,0,0])
        collectionToCompare.rotate(-90,[0,0,1],anchor=[0,0,0])

    # Here the process to create the plot arrays
    xs = linspace(plotBounds[0],plotBounds[1],nPlotPoints)
    zs = linspace(plotBounds[2],plotBounds[3],nPlotPoints)
        
    X,Z = meshgrid(xs,zs)
    
    # Here the process to take the points of the amplitude error and degree error
    Bs = array([[collectionToPlot.getB([xs[x],xyz0,zs[z]]) for x in range(0,len(xs))] for z in range(0,len(zs))])
    Bs = Bs*multiplier
    amp = Bs[:,:,2]
    
    if fieldDif:
        Bs2 = array([[collectionToCompare.getB([xs[x],xyz0,zs[z]]) for x in range(0,len(xs))] for z in range(0,len(zs))])
        Bs2 = Bs2*multiplier
        
        if orientation == 'x':
            amp2 = Bs2[:,:,2]
            #subPlot1.set_title('shows the % difference\n B1tc and B2tp, relative to B2tp \n(comparition in "y" orientation)')
            # shows the 'y' comparation
        elif orientation == 'z':
            amp2 = Bs2[:,:,2]
            #subPlot1.set_title('shows the % difference\n B1tc and B2tp, relative to B2tp \n(comparition in "x" orientation)')
            # shows the 'x' comparation
        elif orientation == 'y':
            amp2 = Bs2[:,:,2]
            #subPlot1.set_title('shows the % difference\n B1tc and B2tp, relative to B2tp \n(comparition in "z" orientation)')
            # shows the 'z' comparation
        ampDif = np.abs(amp2-amp)/amp
        subPlot1.contourf( X, Z, ampDif,10,cmap=plt.cm.jet)
        cp=subPlot1.contour( X, Z, ampDif,30,colors=('k','k'),linestyles = ('-','-'),linewidths=(1,1))
        subPlot1.clabel(cp, fontsize=10, inline=True,fmt= '%4.3f') 
    
    elif compareToCenter:
        # calculating the diffence of the the comparing point to the whole map
        Bcenter = collectionToPlot.getB(centerToCompare)
        Bcenter = Bcenter*multiplier
        Bs = array([[collectionToPlot.getB([xs[x],xyz0,zs[z]]) for x in range(0,len(xs))] for z in range(0,len(zs))])
        Bs = Bs*multiplier
        BsNorm = ((Bs[:,:,0])**2 + (Bs[:,:,1])**2 + (Bs[:,:,2])**2 )**0.5
        BcenterNorm = ((Bcenter[0])**2 + (Bcenter[1])**2 + (Bcenter[2])**2 )**0.5
        
        normDif = (BsNorm-BcenterNorm)/BcenterNorm*100
        # ---------------------------------------------------------------------------
        # Here the process to actualy plot int the screen
        subPlot1.set_title('B vector amplitude Error realative to the "center" (%)')
        subPlot1.contourf( X, Z, normDif,10,cmap=plt.cm.jet)
        cp=subPlot1.contour( X, Z, normDif,30,colors=('k','k'),linestyles = ('-','-'),linewidths=(1,1))
        subPlot1.clabel(cp, fontsize=10, inline=True,fmt= '%4.3f')
        
        # Here calculate the angle to plot
        
        if orientation == 'y':
            subPlot1.set_xlabel('.')
            subPlot2.set_xlabel('.')
            subPlot2.set_ylabel('z axis (mm)')
            subPlot3.set_xlabel('x axis (mm)')
            subPlot3.set_ylabel('z axis (mm)')
            subPlot4.set_xlabel('x axis (mm)')
            subPlot4.set_ylabel('z axis (mm)')
            Bsx = Bs[:,:,0]
            Bsy = Bs[:,:,1]
            Bsz = Bs[:,:,2]
        
        elif orientation == 'x':
            subPlot1.set_xlabel('.')
            subPlot2.set_xlabel('.')
            subPlot2.set_ylabel('y axis (mm)')
            subPlot3.set_xlabel('z axis (mm)')
            subPlot3.set_ylabel('y axis (mm)')
            subPlot4.set_xlabel('z axis (mm)')
            subPlot4.set_ylabel('y axis (mm)')
            Bsx = Bs[:,:,1]
            Bsy = Bs[:,:,2]
            Bsz = Bs[:,:,0]
        
        elif orientation == 'z':
            subPlot1.set_xlabel('.')
            subPlot2.set_xlabel('.')
            subPlot2.set_ylabel('x axis (mm)')
            subPlot3.set_xlabel('y axis (mm)')
            subPlot3.set_ylabel('x axis (mm)')
            subPlot4.set_xlabel('y axis (mm)')
            subPlot4.set_ylabel('x axis (mm)')
            Bsx = Bs[:,:,2]
            Bsy = Bs[:,:,0]
            Bsz = Bs[:,:,1]
        
        if orderMagnitude == 'uT':
            subPlot2.set_title('x component (uT)')
            subPlot3.set_title('y component (uT)')
            subPlot4.set_title('z component (uT)')
            
        else:
            subPlot2.set_title('x component (mT)')
            subPlot3.set_title('y component (mT)')
            subPlot4.set_title('z component (mT)')
        #print(Bs)
        subPlot2.contourf( X, Z, Bsx,10,cmap=plt.cm.jet)
        cp=subPlot2.contour( X, Z, Bsx,30,colors=('k','k'),linestyles = ('-','-'),linewidths=(1,1))
        subPlot2.clabel(cp, fontsize=10, inline=True,fmt= '%4.3f')    
                
        
        subPlot3.contourf( X, Z, Bsy,10,cmap=plt.cm.jet)
        cp=subPlot3.contour( X, Z, Bsy,30,colors=('k','k'),linestyles = ('-','-'),linewidths=(1,1))
        subPlot3.clabel(cp, fontsize=10, inline=True,fmt= '%4.3f')       
        
        
        subPlot4.contourf( X, Z, Bsz,10,cmap=plt.cm.jet)
        cp=subPlot4.contour( X, Z, Bsz,30,colors=('k','k'),linestyles = ('-','-'),linewidths=(1,1))
        subPlot4.clabel(cp, fontsize=10, inline=True,fmt= '%4.3f')         
    else:
        if orderMagnitude == 'uT':
            subPlot1.set_title('B vector amplitude (uT)')
        else:
            subPlot1.set_title('B vector amplitude (mT)')
        subPlot1.contourf( X, Z, amp,10,cmap=plt.cm.jet)
        cp=subPlot1.contour( X, Z, amp,30,colors=('k','k'),linestyles = ('-','-'),linewidths=(1,1))
        subPlot1.clabel(cp, fontsize=10, inline=True,fmt= '%4.3f')
    
    if orientation == 'x':
        collectionToPlot.rotate(-90,[1,0,0],anchor=[0,0,0])
        collectionToCompare.rotate(-90,[1,0,0],anchor=[0,0,0])
        collectionToPlot.rotate(-90,[0,1,0],anchor=[0,0,0])
        collectionToCompare.rotate(-90,[0,1,0],anchor=[0,0,0])
        
    elif orientation == 'z':
        collectionToPlot.rotate(90,[0,0,1],anchor=[0,0,0])
        collectionToCompare.rotate(90,[0,0,1],anchor=[0,0,0])
        collectionToPlot.rotate(90,[0,1,0],anchor=[0,0,0])
        collectionToCompare.rotate(90,[0,1,0],anchor=[0,0,0])
        
        
#%% nextCoilMinimalDiameter function 
def nextCoilMinimalDiameter(coilExtDiameter = 1, coilExtlengthBetweenCoils = 1):
    """ 
    In more than 1d axis Helmholtz Coil context, it returns the minimal diameter 
    that the next axis Coil must have to fit the system, mechanicaly speaking
    
    Parameters
    ----------
    coilExtDiameter : float
        It's the previous helmholtz Coil External Diameter
    coilExtlengthBetweenCoils : float
        It's the length between the two coils of the previous helmholtz coil
        
    Example
    ---------
    test = helmholtzCoil(I = 4, coilIntDiameter = 77.26, wireDiameter=1.59, 
                         loopsInEachEvenLayer=3, evenLayers=3, oddLayers=3, 
                         collectionToAdd=magpy.Collection(), coilsDistance = 77.26/2, orientation = 'z')
    plotB(collectionToPlot = test, optimization= True, figureTittle = 'Test figure', nPlotPoints = 10, percent = True)
        
    test.displaySystem()
    """    
    
    teta = pi/6
    RecBase = np.cos(teta)*coilExtDiameter
    RecHeigth = coilExtlengthBetweenCoils
    newCoilIntDiameter = (RecHeigth**2 + RecBase**2)**0.5
        
    return newCoilIntDiameter

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
#%% spiral function 
def spiral(diameter = 2, spiralExtremities = [0,1], nPoints = 40, nTurns = 4, I = 1, startingSpiralAng = 0, 
           rotationSense = 1, referencePosition = [0,0,0]):
    
    """ 
    Returns an object that contains a spiral line field
    
    Parameters
    ----------
    diameter : float
        Diameter of the spiral
    spiralExtremities : list (float) [2]
        The beginning [0], and the end [1] in 'z' orientation of the spiral
    nPoints : float
        Number of points in the spiral, greater the number, more "realistic" model
        of a spiral you get
    nTurns : positive integer 
        Number of the turns that will have in your spiral
    I : float
        The current amplitude that will create the magnetic field
    startingSpiralAng : float (0~180)
        It's the starting angle of the spiral
    rotationSense : integer (-1,1)
        This multiply the angle, so the sense will be positive increasing or negative increasing
    referencePosition : list (float) [3]
        is the starting point where the spiral will be created, in 'z' orientation
          
    Example
    ---------
    test = spiral(diameter = 2, spiralExtremities = [0,4], nPoints = 40, nTurns = 4, I = 1, startingSpiralAng = 0, 
           rotationSense = 1, referencePosition = [1,1,1])
    test = magpy.Collection(test)
    test.displaySystem()
    """    
    
    # Creating the linear space wheres that will creates the spiral
    spiralExtremities = linspace(spiralExtremities[0],spiralExtremities[1],nPoints)
    phi = linspace(0,rotationSense*2*pi*nTurns,nPoints)
    
    # Creating the spiral around the 'z' axis
    vertices = [[diameter/2*cos(phi[i]),diameter/2*sin(phi[i]),spiralExtremities[i]] for i in range(0,nPoints)]
    spi = magpy.source.current.Line(curr=I,vertices=vertices, pos=referencePosition)
    spi.rotate(startingSpiralAng,[0,0,1],anchor=referencePosition)
    
    return spi

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
#%% twistedPair function 
def twistedPair(diameter = 2, nPoints = 40, nTurns = 10, I = 4, collectionToAdd = magpy.Collection(), referencePosition= [0,0,0], end = [0,0,1]):
    
    """ 
    Returns a Collection that is a previous Collection with a twisted pair plus that
    *OBS: still have to fix the rotations
    
    Parameters
    ----------
    diameter : float
        Diameter of the spiral
    spiralExtremities : list (float) [2]
        The beginning [0], and the end [1] in 'z' orientation of the spiral
    nPoints : float
        Number of points in the spiral, greater the number, more "realistic" model
        of a spiral you get
    nTurns : positive integer 
        Number of the turns that will have in your spiral
    I : float
        The current amplitude that will create the magnetic field
    startingSpiralAng : float (0~180)
        It's the starting angle of the spiral
    rotationSense : integer (-1,1)
        This multiply the angle, so the sense will be positive increasing or negative increasing
    referencePosition : list (float) [3]
        is the starting point where the spiral will be created, in 'z' orientation
          
    Example
    ---------
    test = twistedPair(diameter = 2, nPoints = 40, nTurns = 4, I = 4, collectionToAdd = magpy.Collection(), referencePosition= [0,0,0], end = [0,0,3])
    test.displaySystem()
    
    """    
    
    # Find the rotations angles for the right positioning of the twisted pair ------------------------------
    spiralExtremities = np.zeros(2)
    spiralExtremities[0] = 0
    spiralExtremities[1] = ((end[0] - referencePosition[0])**2 + (end[1] - referencePosition[1])**2 + (end[2] - referencePosition[2])**2 )**0.5

    print(spiralExtremities[1])         # debug
    # *OBS: still have to fix the rotations
    """
    phi = -np.arccos((end[2] - referencePosition[2])/spiralExtremities[1])
    print(phi)                          # debug
    teta = np.arccos((end[0] - referencePosition[0])/spiralExtremities[1]/np.sin(phi))*180/pi
    print(teta)                         # debug
    phi = phi*180/pi
    
    # -------------------------------------------------------------------------------------
    
    # Doing the inverse rotations -------------------------------------
    collectionToAdd.rotate(-teta,[0,0,1],anchor=referencePosition)
    collectionToAdd.rotate(-phi,[1,0,0],anchor=referencePosition)
    # -------------------------------------------------------------------------------------
    """
    # Doing one spiral in the 'z' positive orientation
    spi = spiral(diameter = diameter, spiralExtremities = [spiralExtremities[0],spiralExtremities[1]], nPoints = nPoints, 
                 nTurns = nTurns, I = I, startingSpiralAng = 0, referencePosition = referencePosition)
    collectionToAdd.addSources(spi)

    # Doing the second spiral in the 'z' negative orientation
    spi = spiral(diameter = diameter, spiralExtremities = [spiralExtremities[1],spiralExtremities[0]], nPoints = nPoints, 
                 nTurns = nTurns, I = I, startingSpiralAng = 180, rotationSense = -1, referencePosition = referencePosition)
    collectionToAdd.addSources(spi)
    """
    # Doing the rotations -------------------------------------
    collectionToAdd.rotate(phi,[1,0,0],anchor=referencePosition)
    collectionToAdd.rotate(teta,[0,0,1],anchor=referencePosition)
    # -------------------------------------------------------------------------------------
    """
    return collectionToAdd

# -----------------------------------------------------------------------------------------------------------------------------------------------------------

def WhatTimeIsIt():
    justTime = str(datetime.datetime.now())
    #return justTime[0:22]
    return time.time()
#%% desiredCoilDiameter function 
def desiredCoilDiameter(nIteration = 1000, breakCondition = 0.005, desiredField = 1, initialValue = 0, 
                        I = 1, coilIntDiameter = 1, wireDiameter=1, loopsInEachEvenLayer=3, evenLayers=3, oddLayers=3, 
                        orientation='z', relativeDistance = 0.5):
    """ 
    Returns a Collection that is a previous Collection with a twisted pair plus that
    *OBS: still have to fix the rotations
    
    Parameters
    ----------
    nIteration : integer
        Determine how many iterations the loop will do
    breakCondition : float
        It's the B - Ba break condition
    coilIntDiameterStep : float
        The step that the coil diameter tested will increase or decrease
    desiredField : float
        The field amplitude value that the function will search for
    I : float
        The Helmhotz Coil Current in Amperes.
    coilIntDiameter : float
        The internal diameter of the coil in mm (initial value)
    wireDiameter : float, optimizationional
        Diameter of the coil wire in mm.
    loopsInEachEvenLayer : non-zero integer, default 1 
        How many loops in a single even layer
    evenLayers : integer, default 1
        Integer number of even layers.
    oddLayers : integer, default 1
        Integer number of odd layers.
    collectionToAdd : Magpy Collection, default empty
        Append coil source to a magpy collection.
    coilsDistance : float , default 1
        Distance between the coil of the helmholtz coil
    orientation : {'x', 'y', 'z'}, optimizationional
        Change the orientation of the resulting coil having its center parallel with Plane X, Y, or Z.
    
    *OBS: the odd layers always must be evenLayers -1 or = evenLayers 
    
    Example
    ---------
    test = twistedPair(diameter = 2, nPoints = 40, nTurns = 4, I = 4, collectionToAdd = magpy.Collection(), referencePosition= [0,0,0], end = [0,0,3])
    test.displaySystem()
    
    """    
    
    # Setting the initial values of B and Ba
    samplingTime = 100
    previousError = 0
    previousCoilIntDiameter = initialValue
    i = 0
    # The iteration loop that will find the correct HelmoltzCoil to get the I desired with B desired together
    while (1):
        
        
        coilIntDiameter = -previousError*samplingTime + previousCoilIntDiameter
        if coilIntDiameter < 0.1:
            coilIntDiameter = 0.1
        #b = WhatTimeIsIt()
        # Creates a new helmholtzCoil, to search for the right coilintDiameter value
        EmptyCol = magpy.Collection()
        Coil = helmholtzCoil(I = I, coilIntDiameter = coilIntDiameter, wireDiameter=wireDiameter, 
                             loopsInEachEvenLayer=loopsInEachEvenLayer, evenLayers=evenLayers, oddLayers=oddLayers, 
                             collectionToAdd=EmptyCol, coilsDistance = coilIntDiameter/2, orientation = orientation,
                             relativeDistance = relativeDistance)
        # --------------------------------------------------------------------------------------------
        #a = WhatTimeIsIt()
        #print(a-b)
        # Store the B in Ba and calculate the new B
        field = Coil.getB([0,0,0])
        if(orientation == 'z'):
            field = field[2]
        elif(orientation == 'y'):
            field = field[1]
        elif(orientation == 'x'):
            field = field[0]
        else:
            field = norm(field)
        
        error = desiredField - field
        #print(error)
        #print(coilIntDiameter)
        previousError  = error
        previousCoilIntDiameter = coilIntDiameter
        
        i += 1
        if i > nIteration:
            break
        if np.linalg.norm(error) < breakCondition:
            break
        #print(B)
    # -----------------------------------------------------------------------------------------------------
    
    return coilIntDiameter

 # Setting input values -------------------------------------------------------------------------------------------------------------------------
 
def coilDesign(nIteration = 40,
               initialValue = 0.1,
               desiredField = 0.1, 
               I = 0.1, 
               wireDiameter=0.6, 
               loopsInEachEvenLayer=5, 
               evenLayers=5, oddLayers=5, 
               orientation='z',
               relativeDistance = 1/3,
               maxError = 0.1,
               color = 'k',
               wireResistence = 0.5):
    
    
    # ----------------------------------------------------------------------------------------------------------------------------------------------
    # Finding the right coil internal diameter to have the desired field with the zCoilCurrent in A
    newCoilIntDiameter = desiredCoilDiameter(nIteration = nIteration, 
                                              breakCondition = desiredField*maxError/100, 
                                              initialValue = initialValue,
                                              desiredField = desiredField, 
                                              I = I, 
                                              coilIntDiameter = initialValue, 
                                              wireDiameter=wireDiameter, 
                                              loopsInEachEvenLayer=loopsInEachEvenLayer, 
                                              evenLayers=evenLayers, oddLayers=oddLayers, 
                                              orientation=orientation,
                                              relativeDistance = relativeDistance)
    # ---------------------------------------------------------------------------------------------
    
    
    # Constructing the desired coil ---------------------------------------------------------------
    Coil = helmholtzCoil(color = color,I = I, 
                           coilIntDiameter = newCoilIntDiameter, 
                           wireDiameter=wireDiameter, 
                           loopsInEachEvenLayer=loopsInEachEvenLayer, 
                           evenLayers=evenLayers, oddLayers=oddLayers, 
                           coilsDistance=newCoilIntDiameter/2, 
                           orientation=orientation, 
                           relativeDistance = relativeDistance)
    # --------------------------------------------------------------------------------------------- 
    
    # calculating the energy dissipation by joule effect --------------------------------------------------------------------------------------------
    coilExtDiameter = newCoilIntDiameter
    wireLenght = 0
    
    for i in range(0,evenLayers + oddLayers):
        if (i % 2) == 0:
            wireLenght += coilExtDiameter*pi*loopsInEachEvenLayer
        elif (i % 2) == 1:
            wireLenght += coilExtDiameter*pi*(loopsInEachEvenLayer-1)
        coilExtDiameter += ((3)**(1/2))*wireDiameter/2
    
    wireLenghtMeters = 2*wireLenght/1000
    coilResistence = wireResistence*wireLenghtMeters
    coilMaxPower = coilResistence*(I**2)

    
    return Coil, newCoilIntDiameter, coilExtDiameter, coilResistence, coilMaxPower, wireLenghtMeters
