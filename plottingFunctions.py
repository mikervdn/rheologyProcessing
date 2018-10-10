#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 09:17:41 2018

@author: mike
"""

import numpy as np 
import os
from matplotlib import pyplot
import fittingFunctions


    
def shearStressVsViscosityPlotter(topDir,titleOfPlot="normal",rescaledViscosity=False,noErrors=False):
    
    
    #Get a list of all the subdirectories
    listOfSubDirs = next(os.walk(topDir))[1]
    listOfSubDirs = [float(i) for i in listOfSubDirs]
    listOfSubDirs.sort()
    listOfSubDirs = [str(i) for i in listOfSubDirs]
    listOfSubDirs = listOfSubDirs[::-1]
    
    f, ax = pyplot.subplots()
    
    
    for packingFraction in listOfSubDirs:
        
        if rescaledViscosity == False:
            loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedData.npy" ))
        else:
            loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedDataVR.npy" ))
            
            
        if noErrors == False:
            if rescaledViscosity ==False:
                loadInErrordata = np.load(os.path.join(topDir, packingFraction , r"errorsData.npy" ))
            else:
                loadInErrordata = np.load(os.path.join(topDir, packingFraction , r"errorsDataVR.npy" ))
        
                
        #First the plot as a function of shear stress.
        if noErrors == False:
            ax.errorbar(loadInAveragedData[:,4],loadInAveragedData[:,2],xerr= loadInErrordata[:,4],yerr=loadInErrordata[:,2],label=packingFraction,marker='o',markersize=5)
        else:
            ax.plot(loadInAveragedData[:,4],loadInAveragedData[:,2],label=packingFraction,marker='o')
        #Next plot the viscosity as a function of shear rate,


    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.xlabel("Shear Stress [Pa]")
    
    if rescaledViscosity==False:
        pyplot.ylabel("Viscosity [Pa s]")
    else:
        pyplot.ylabel(r"Rescaled Viscosity $\frac{\eta}{\eta_0}$")
        
    pyplot.title(titleOfPlot)
        
    pyplot.grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    return f, ax

def tauMinTauMaxPlotter(topDir,tauMin,tauMax,titleOfPlot="normal",rescaledViscosity=False,noErrors=False):
    
    
    listOfSubDirs = next(os.walk(topDir))[1] 
    
    minViscosityList = np.zeros(len(listOfSubDirs))
    maxViscosityList = np.zeros(len(listOfSubDirs))
    
    (phiVsMinViscosity,phiVsMaxViscosity,phiVsMinViscosityErrors,phiVsMaxViscosityErrors) = fittingFunctions.newtonianPlateauFinder(topDir,minViscosityList,maxViscosityList,rescaledViscosity,tauMin,tauMax)
    
    f, ax = shearStressVsViscosityPlotter(topDir,titleOfPlot,rescaledViscosity,noErrors)
    
    (numSystems,_) = np.shape(phiVsMaxViscosity)
    
    for i in range(0,numSystems):
        
        ax.plot(phiVsMinViscosity[i,:][2],phiVsMinViscosity[i,:][1],'o',markerSize=15,color="red")
        ax.plot(phiVsMaxViscosity[i,:][2],phiVsMaxViscosity[i,:][1],'o',markerSize=15,color="blue")
        
    return f,ax
    
    

    
def shearRateVsViscosityPlotter(topDir,solvent="normal",rescaleViscosity=False):
    
    
    #Get a list of all the subdirectories
    listOfSubDirs = next(os.walk(topDir))[1]
    listOfSubDirs = [float(i) for i in listOfSubDirs]
    listOfSubDirs.sort()
    listOfSubDirs = [str(i) for i in listOfSubDirs]
    listOfSubDirs = listOfSubDirs[::-1]
    
        
    ax = pyplot.subplot(111)
    
    for packingFraction in listOfSubDirs:
        
        if rescaleViscosity ==False:
            loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedData.npy" ))
            loadInErrordata = np.load(os.path.join(topDir, packingFraction , r"errorsData.npy" ))
        else:
            loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedDataVR.npy" ))
            loadInErrordata = np.load(os.path.join(topDir, packingFraction , r"errorsDataVR.npy" ))
            
        
                
        #First the plot as a function of shear stress.
        ax.errorbar(loadInAveragedData[:,3],loadInAveragedData[:,2],xerr= loadInErrordata[:,3],yerr=loadInErrordata[:,2],label=packingFraction,marker='o')
        #Next plot the viscosity as a function of shear rate,
    
    
    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.xlabel(r"Shear Rate [$\frac{1}{s}$]")
    
    if rescaleViscosity==False:
        pyplot.ylabel("Viscosity [Pa s]")
    else:
        pyplot.ylabel(r"Rescaled Viscosity $\frac{\eta}{\eta_0}$")
        
    if solvent == "GdCl":
        pyplot.title("Flow Curves for GdCl Cornstarch Suspensions")
    if solvent =="normal":
        pyplot.title("Flow Curves for Normal Cornstarch Suspensions")
        
    pyplot.grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    

        

def viscosityDivergencePlot(viscosityValues,packingFractionValues,viscosityErrors,labelsForCurves,fixAlpha=False,titleOfPlot ="Newtonian Plateau's Divergences"):
    
    
    shapeOfVis = np.shape(viscosityValues)
    
    if np.shape(shapeOfVis)[0]==1:
        numSystems=1
    else:
        numSystems = shapeOfVis[0]
    
    ax = pyplot.subplot(111)
    
    if fixAlpha== False:
        def divergingViscosityFunc(x,phiJ, alpha):
            return ( 1 - (x/phiJ) )**( -alpha )
    else:
        def divergingViscosityFunc(x,phiJ):
            return ( 1 - (x/phiJ) )**(-fixAlpha)
    
    
        
        
    
    phiJHolder = np.zeros(numSystems)
    phiJErrorHolder = np.zeros(numSystems)
    
    if fixAlpha == False:
        alphaHolder = np.zeros(numSystems)
        alphaErrorHolder = np.zeros(numSystems)
    
    
    for i in range(0,numSystems):
        print(numSystems)
        if numSystems == 1:
            currentViscosityErrors = viscosityErrors
            currentPhiValues = packingFractionValues
            currentViscosityValues = viscosityValues
        else:
            currentViscosityErrors = viscosityErrors[i,:]
            currentPhiValues = packingFractionValues[i,:]
            currentViscosityValues = viscosityValues[i,:]
            
            
        currentViscosityValues = currentViscosityValues[~np.isnan(currentViscosityValues)]
        currentViscosityErrors = currentViscosityErrors[~np.isnan(currentViscosityValues)] 
        currentPhiValues = currentPhiValues[~np.isnan(currentPhiValues)]
        
        if fixAlpha == False:
            (phiJHolder[i],alphaHolder[i],phiJErrorHolder[i],alphaErrorHolder[i]) = fittingFunctions.viscosityDivergenceFit(currentPhiValues,currentViscosityValues,currentViscosityErrors,False,False)
            ax.plot(np.linspace(0,phiJHolder[i]-.001,1000),divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i],alphaHolder[i]),label=labelsForCurves[i] +" Fit" )

            minimumError = divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]-phiJErrorHolder[i],alphaHolder[i]+alphaErrorHolder[i])
            maximumError = divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]+phiJErrorHolder[i],alphaHolder[i]-alphaErrorHolder[i])
            #pyplot.plot(np.linspace(0,phiJHolder[i]-.001,1000),minimumError,label="minimumError")
            #pyplot.plot(np.linspace(0,phiJHolder[i]-.001,1000), maximumError,label="maximumError")
            
            ax.fill_between(np.linspace(0,phiJHolder[i]-.001,1000), minimumError, maximumError, where=maximumError <= minimumError, facecolor='lightblue')
        else:
            (phiJHolder[i],phiJErrorHolder[i]) = fittingFunctions.viscosityDivergenceFit(currentPhiValues,currentViscosityValues,currentViscosityErrors,False,fixAlpha)
            minimumError = divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]-phiJErrorHolder[i])
            maximumError = divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]+phiJErrorHolder[i])
            ax.fill_between(np.linspace(0,phiJHolder[i]-.001,1000), minimumError, maximumError, where=maximumError <= minimumError, facecolor='lightblue')
            
            ax.plot(np.linspace(0,phiJHolder[i]-.001,1000),divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]),label=labelsForCurves[i] +" Fit" )
        print(currentViscosityErrors)
        ax.errorbar(currentPhiValues, currentViscosityValues,fmt='o',yerr=currentViscosityErrors,label=labelsForCurves[i])

        
        
    ax.legend()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    pyplot.title(titleOfPlot)
    pyplot.xlabel("Packing Fraction $\phi$")
    pyplot.ylabel("Rescaled Viscosity $\eta_R$")
    pyplot.xlim(np.nanmin(packingFractionValues)-.05,np.nanmax(packingFractionValues)+.04)
    pyplot.ylim(0,np.nanmax(viscosityValues)+20)
    if fixAlpha == False:
        return (phiJHolder,alphaHolder,phiJErrorHolder,alphaErrorHolder)
    else:
        return (phiJHolder,phiJErrorHolder)
    

def logarithmicPlotter(viscosityValues0,viscosityValues1,packingFractionValues0,packingFractionValues1,alpha0,phiJ0,alpha1,phiJ1,highOrLow):
    
    logViscosityValues0 = np.log10(viscosityValues0)
    logPackingFraction0 = np.log10(1-packingFractionValues0/phiJ0)
    
    logViscosityValues1 = np.log10(viscosityValues1)
    logPackingFraction1 = np.log10(1-packingFractionValues1/phiJ1)
    
    
    pyplot.plot(logPackingFraction0,logViscosityValues0,'o')
    
    pyplot.plot(logPackingFraction1,logViscosityValues1,'o')
    #pyplot.title()
    if highOrLow == "high":
        pyplot.xlabel(r"$\log(1-\frac{\phi}{\phi_M})$")
        pyplot.title("Frictional Newtonian Regimes")
        
    if highOrLow == "low":
        pyplot.xlabel(r"$\log(1-\frac{\phi}{\phi_0})$")
        pyplot.title("Frictionless Newtonian Regimes")
    
        
    pyplot.legend(("Pure Cornstarch Suspension","GdCl Cornstarch Suspension"))
    pyplot.ylabel(r"$\log(\eta)$")
    
    font = {'family' : 'normal',

        'size'   : 10}

    pyplot.rc('font', **font)
    #pyplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    
def logarithmicPlotterAllLines(phiVsMinViscosityG,phiVsMinViscosityC,phiVsMaxViscosityG,phiVsMaxViscosityC,phiMG,alphaMG,phi0G,alpha0G,phiMC,alphaMC,phi0C,alpha0C,suspendViscosityC,suspendViscosityG):
     
    
    logViscosityValuesC_0= np.log10(phiVsMinViscosityC[:,1])
    logPackingFractionC_0 = np.log10(1-phiVsMinViscosityC[:,0]/phi0C)
    logViscosityValuesG_0= np.log10(phiVsMinViscosityG[:,1])
    logPackingFractionG_0 = np.log10(1-phiVsMinViscosityG[:,0]/phi0G)
    
    logViscosityValuesC_M= np.log10(phiVsMaxViscosityC[:,1])
    logPackingFractionC_M = np.log10(1-phiVsMaxViscosityC[:,0]/phiMC)
    logViscosityValuesG_M= np.log10(phiVsMaxViscosityG[:,1])
    logPackingFractionG_M = np.log10(1-phiVsMaxViscosityG[:,0]/phiMG)
    
    
    
    pyplot.plot(logPackingFractionC_0,logViscosityValuesC_0,'o',label="Frictionless Cornstarch")
    pyplot.plot(logPackingFractionG_0,logViscosityValuesG_0,'o',label="Frictionless Cornstarch/GdCl")
    pyplot.plot(logPackingFractionC_M,logViscosityValuesC_M,'o',label="Frictional Cornstarch")
    pyplot.plot(logPackingFractionG_M,logViscosityValuesG_M,'o',label="Frictional Cornstarch/GdCl")
    
    
    pyplot.legend()
    pyplot.xlabel(r"$\log(1-\frac{\phi}{\phi_J})$")
    pyplot.ylabel(r"$\log(\eta)$")



def sideBySideCurvePlot(topDir1,topDir2,topDir3,outputDir,rescaledViscosity=False,allPlots=False,evenCurves="none"):
    
    #get list of all of the subdirectories in topDir1 and topDir2
    dirsIn1 = next(os.walk(topDir1))[1]
    dirsIn2 = next(os.walk(topDir2))[1]
    dirsIn3 = next(os.walk(topDir3))[1]
    #They need to be sorted so we have to convert them from strings to floats and then back
    dirsIn1 = np.sort([float(i) for i in dirsIn1])
    dirsIn2 = np.sort([float(i) for i in dirsIn2])
    dirsIn3 = np.sort([float(i) for i in dirsIn3])
    #Back to strings!
    dirsIn1 = [str(i) for i in dirsIn1]
    dirsIn2 = [str(i) for i in dirsIn2]
    dirsIn3 = [str(i) for i in dirsIn3]
    
    pyplot.figure(figsize=(14.0, 5.0))

    
    markersForPlots = ["o","x","D","h","+","P"]
    
    for i in range(0,len(dirsIn1)):
        
        if allPlots == True:
            if evenCurves == "even":
                if i%2!=0:
                    continue
            if evenCurves == "odd":
                if i%2==0:
                    continue
                

            
        

        if rescaledViscosity == True:
            dir1Data = np.load(os.path.join(topDir1,dirsIn1[i], r"averagedDataVR.npy"))
            dir1Errors = np.load(os.path.join(topDir1,dirsIn1[i], r"errorsDataVR.npy"))
            dir2Data = np.load(os.path.join(topDir2,dirsIn2[i], r"averagedDataVR.npy"))
            dir2Errors = np.load(os.path.join(topDir2,dirsIn2[i], r"errorsDataVR.npy"))
            dir3Data = np.load(os.path.join(topDir3,dirsIn3[i], r"averagedDataVR.npy"))
            dir3Errors = np.load(os.path.join(topDir3,dirsIn3[i], r"errorsDataVR.npy"))
        else:
            dir1Data = np.load(os.path.join(topDir1,dirsIn1[i], r"averagedData.npy"))
            dir1Errors = np.load(os.path.join(topDir1,dirsIn1[i], r"errorsData.npy"))
            dir2Data = np.load(os.path.join(topDir2,dirsIn2[i], r"averagedData.npy"))
            dir2Errors = np.load(os.path.join(topDir2,dirsIn2[i], r"errorsData.npy"))
            dir3Data = np.load(os.path.join(topDir3,dirsIn3[i], r"averagedData.npy"))
            dir3Errors = np.load(os.path.join(topDir3,dirsIn3[i], r"errorsData.npy"))
        
        
        pyplot.errorbar(dir1Data[:,4],dir1Data[:,2],yerr=dir1Errors[:,2],label=dirsIn1[i]+"C",marker=markersForPlots[i],color="black", markersize=10)
        pyplot.errorbar(dir2Data[:,4],dir2Data[:,2],yerr=dir2Errors[:,2],label=dirsIn2[i]+"GdCl",marker=markersForPlots[i],color="red", markersize=10)
        pyplot.errorbar(dir3Data[:,4],dir3Data[:,2],yerr=dir3Errors[:,2],label=dirsIn3[i]+"NH4",marker=markersForPlots[i],color="blue", markersize=10)
        
        figName = str(round(float(dirsIn1[i])))
        pyplot.yscale('log')
        pyplot.xscale('log')
        pyplot.xlabel(r"Shear Stress [Pa]")
        
        
        if rescaledViscosity == True:
            pyplot.ylabel(r"Rescaled Viscosity $\frac{\eta}{\eta_0}$")
        else:
            pyplot.ylabel(r"Viscosity $\eta$ [Pa s]")
            
        pyplot.legend()
        
        
        if allPlots == False:
            if rescaledViscosity == True:
                pyplot.savefig(os.path.join(outputDir,figName+"rescaled"))
            else:
                pyplot.savefig(os.path.join(outputDir,figName))
            pyplot.clf()
        else:
            pyplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                    
        
    if allPlots == True:
        pyplot.savefig(os.path.join(outputDir,"all"+evenCurves+"plots"))
        
        

def shearJammingPhaseDiagram(phi0,phiM,tauStar,alpha):
    '''
    This generates a plot of the shear jamming phase diagram using the Wyart-Cates 
    model.  This model is outlined in the paper found here:
    
    https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.112.098302
    
    This model gives the shear rate as a function of the density of your suspension
    and the as a function of the shear stress:
        gammaDot = (1/eta_0)*tau*(1+phi/(phiM+(phi0-phiM)*e^(-tau/tauStar)))^alpha
        where gammaDot is the shear rate, eta_0 is the viscosity of the suspending solvent
        tau is the shear stress, phi is the densty, phi0 is the frictionless jamming
        density, phiM is the frictional jamming density, and alpha is a free parameter in the
        model.
        
    
    '''
    
    #This is the y-values, or the stresses used for each curve
    stressValues = np.linspace(10**-5,10**5,100000)
    
    
    #Find the DST line as a function of the stressValues defined above.  This function
    #returns packing fractions as a function of shear stress
    def DSTLine(tau):
        return (np.exp(-tau/tauStar)*((phi0-phiM)+phiM*np.exp(tau/tauStar))**2)/(phiM*(np.exp(tau/tauStar)-1-alpha*tau/tauStar)+phi0*(1+alpha*tau/tauStar))
    
    #Find the Shear Jamming line as a function of the stressValues defined above.  
    #This function returns packing fractions as a function of shear stress
    def shearJammingLine(tau):
        return phiM+(phi0-phiM)*np.exp(-tau/tauStar)
    
    #This is the x-values, or the densities that define the DST line and SJ line
    phiValuesDST = DSTLine(stressValues)
    phiValuesSJ = shearJammingLine(stressValues)
    
    #This generates the vertical line that defines frictionless jamming
    phiValuesJamming = phi0*np.ones(100000)
    
    #This is another vertical line that is not seen in the plot but acts as boundary for the fill_between function 
    phiValuesBeyondJamming = phiValuesJamming+0.06
    
    #Plot all four curves, the DST curve, the SJ curve, and the two vertical lines.
    pyplot.semilogy(phiValuesDST,stressValues,phiValuesSJ,stressValues,phiValuesJamming,stressValues,phiValuesBeyondJamming,stressValues,label = ['DST', 'SJ','J'],color="black")
    
    #This shades the region between DST and SJ green (This is the DST region)
    pyplot.fill_betweenx(stressValues,phiValuesDST,phiValuesSJ,where=phiValuesDST <= phiValuesSJ, facecolor='lightblue', interpolate=True)
    #This shades the region between SJ and frictionless jamming blue (This is the SJ region)
    pyplot.fill_betweenx(stressValues,phiValuesSJ,phiValuesJamming,where=phiValuesSJ <= phiValuesJamming, facecolor='lightgreen', interpolate=True)
    #This shades the region between frictionless jamming and vertical line out of sight grey (This is the frictionlesss jamming region)
    pyplot.fill_betweenx(stressValues,phiValuesJamming,phiValuesBeyondJamming,where=phiValuesJamming <= phiValuesBeyondJamming, facecolor='olive', interpolate=True)
        
    
    pyplot.ylim(10**-1, 1000)
    pyplot.xlim(phiM-.1,phi0+.05)
    
    pyplot.xlabel(r"Packing Fraction $\phi$")
    pyplot.ylabel(r"Shear Stress $ \tau $ [Pa]")
    pyplot.title("Shear Jamming Phase Diagram")
    
    
    
    
    
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]



def shearJammingPhaseDiagramTwoSystems(phi01,phiM1,tauStar1,alpha1,phi02,phiM2,tauStar2,alpha2,labels,DST=False,SJ=False):
    '''
    This generates a plot of the shear jamming phase diagram using the Wyart-Cates 
    model.  This model is outlined in the paper found here:
    
    https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.112.098302
    
    This model gives the shear rate as a function of the density of your suspension
    and the as a function of the shear stress:
        gammaDot = (1/eta_0)*tau*(1+phi/(phiM+(phi0-phiM)*e^(-tau/tauStar)))^alpha
        where gammaDot is the shear rate, eta_0 is the viscosity of the suspending solvent
        tau is the shear stress, phi is the densty, phi0 is the frictionless jamming
        density, phiM is the frictional jamming density, and alpha is a free parameter in the
        model.
        
    
    '''
    
    #This is the y-values, or the stresses used for each curve
    stressValues = np.linspace(10**-5,10**5,100000)
    
    
    #Find the DST line as a function of the stressValues defined above.  This function
    #returns packing fractions as a function of shear stress
    def DSTLine(tau,tauStar,phi0,phiM,alpha):
        return (np.exp(-tau/tauStar)*((phi0-phiM)+phiM*np.exp(tau/tauStar))**2)/(phiM*(np.exp(tau/tauStar)-1-alpha*tau/tauStar)+phi0*(1+alpha*tau/tauStar))
    
    #Find the Shear Jamming line as a function of the stressValues defined above.  
    #This function returns packing fractions as a function of shear stress
    def shearJammingLine(tau,tauStar,phi0,phiM,alpha):
        return phiM+(phi0-phiM)*np.exp(-tau/tauStar)
    
    
    
    
    
    #This is the x-values, or the densities that define the DST line and SJ line for system 1.
    phiValuesDST1 = DSTLine(stressValues,tauStar1,phi01,phiM1,alpha1)
    phiValuesSJ1 = shearJammingLine(stressValues,tauStar1,phi01,phiM1,alpha1)
    #This is the x-values, or the densities that define the DST line and SJ line for system 2.
    phiValuesDST2 = DSTLine(stressValues,tauStar2,phi02,phiM2,alpha2)
    phiValuesSJ2 = shearJammingLine(stressValues,tauStar2,phi02,phiM2,alpha2)
    
    #This generates the vertical line that defines frictionless jamming
    phiValuesJamming1 = phi01*np.ones(100000)
    phiValuesJamming2 = phi02*np.ones(100000)

 

    if DST ==True:
        pyplot.semilogy(phiValuesDST1,stressValues,label = 'DST-'+ labels[0],color="black")
        #Plot all four curves, the DST curve, the SJ curve, and the two vertical lines for the first system.
        pyplot.semilogy(phiValuesDST2,stressValues,label = 'DST-'+ labels[1],color="red")
        pyplot.fill_betweenx(stressValues,phiValuesDST1,phiValuesSJ1,where=phiValuesDST1 <= phiValuesSJ1, facecolor='lightblue', interpolate=True)
        pyplot.fill_betweenx(stressValues,phiValuesDST2,phiValuesSJ2,where=phiValuesDST2 <= phiValuesSJ2, facecolor='lightgreen', interpolate=True)
        pyplot.semilogy(phiValuesSJ1,stressValues,label = 'SJ-'+ labels[0],color="black")
        pyplot.semilogy(phiValuesSJ2,stressValues,label = 'SJ-'+ labels[1],color="red")
    if SJ == True:
        pyplot.semilogy(phiValuesSJ1,stressValues,label = 'SJ-'+ labels[0],color="black")
        pyplot.semilogy(phiValuesSJ2,stressValues,label = 'SJ-'+ labels[1],color="red")
        pyplot.semilogy(phiValuesJamming1,stressValues,label = 'J-'+ labels[0],color="black")
        pyplot.semilogy(phiValuesJamming2,stressValues,label = 'J-'+ labels[1],color="red")
        pyplot.fill_betweenx(stressValues,phiValuesSJ1,phiValuesJamming1,where=phiValuesSJ1 <= phiValuesJamming1, facecolor='lightgreen', interpolate=True)
        pyplot.fill_betweenx(stressValues,phiValuesSJ2,phiValuesJamming2,where=phiValuesSJ2 <= phiValuesJamming2, facecolor='lightgreen', interpolate=True)
        
        
                                                                                                                                            
    pyplot.legend()


    phiM = min(phiM1,phiM2)
    phi0 = max(phi01,phi02)
    
    pyplot.ylim(10**-1, 1000)
    pyplot.xlim(phiM-.05,phi0+.05)
    
    pyplot.xlabel(r"Packing Fraction $\phi$")
    pyplot.ylabel(r"Shear Stress $ \tau $ [Pa]")
    5
    if DST==True and SJ ==True:
        pyplot.title("Shear Jamming Phase Diagram")
    
    if DST==True and SJ==False:
        pyplot.title("Discontinuous Shear Thickening Regions")
        
    if DST==False and SJ==True:
        pyplot.title("Shear Jamming Regions")
        
    
    
    
def WCModelComparison(topDir,outputDir,phi0,phiM,tauStar,alpha,beta=False):
    
    
    
    listOfSubDirs = next(os.walk(topDir))[1]

    
    stressValues = np.linspace(.0001,1e3,10000)
    
    if beta==False:
        def rescaledNewtonianViscosity(stress,phi):
            return (1-phi/(phiM+(phi0-phiM)*np.exp(-stress/tauStar)))**(-alpha)
    else:
        def rescaledNewtonianViscosity(stress,phi):
            return (1-phi/(phi0+(phiM-phi0)*np.exp((-tauStar/stress))**beta))**(-alpha)
    
    for packingFraction in listOfSubDirs:
        loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedDataVR.npy" ))
        loadInErrors = np.load(os.path.join(topDir, packingFraction , r"errorsDataVR.npy" ))
        
        phiValues = (float(packingFraction)/100)*np.ones(10000)
        
        
        
        #Plot Wyart-Cates Prediction
        pyplot.plot(stressValues,rescaledNewtonianViscosity(stressValues,phiValues))
        #Plot Actual Data
        pyplot.errorbar(loadInAveragedData[:,4],loadInAveragedData[:,2],yerr=loadInErrors[:,2],label=packingFraction,marker='o',markersize=5)
        
        pyplot.title(packingFraction)   
        pyplot.xlabel("Shear Stress [Pa]")
        pyplot.ylabel("Rescaled Viscosity")
        pyplot.yscale('log')
        pyplot.xscale('log')
        pyplot.ylim(1,1e3)
        pyplot.xlim(1e-1,1e3)
        pyplot.savefig(os.path.join(outputDir,str(float(packingFraction))+".png"))
        pyplot.clf()
        
        
def WCModelComparisonTwoSystems(topDir1,topDir2,outputDir,phi01,phiM1,tauStar1,alpha1,phi02,phiM2,tauStar2,alpha2,phi01E,phi02E,phiM1E,phiM2E,tauStar1E,tauStar2E,alpha1E,alpha2E,stretchedExponential=False,beta1=False,beta2=False):
    
    
    
    dirsIn1 = next(os.walk(topDir1))[1]
    dirsIn2 = next(os.walk(topDir2))[1]
    
    #They need to be sorted so we have to convert them from strings to floats and then back
    dirsIn1 = np.sort([float(i) for i in dirsIn1])
    dirsIn2 = np.sort([float(i) for i in dirsIn2])
    #Back to strings!
    dirsIn1 = [str(i) for i in dirsIn1]
    dirsIn2 = [str(i) for i in dirsIn2]

    
    stressValues = np.linspace(.0001,1e3,10000)
    
    if stretchedExponential == False:
        def rescaledNewtonianViscosity1(stress,phi):
            return (1-phi/(phiM1+(phi01-phiM1)*np.exp(-stress/tauStar1)))**(-alpha1)
        
        def rescaledNewtonianViscosity2(stress,phi):
            return (1-phi/(phiM2+(phi02-phiM2)*np.exp(-stress/tauStar2)))**(-alpha2)
    else:
        def rescaledNewtonianViscosity1(stress,phi):
            return (1-phi/(phi01+(phiM1-phi01)*np.exp((-stress/tauStar1)**beta1)))**(-alpha1)
        
        def rescaledNewtonianViscosity2(stress,phi):
            return (1-phi/(phi02+(phiM2-phi02)*np.exp((-stress/tauStar2)**beta2)))**(-alpha2)
        
    
    
    
    
    for i in range(0,len(dirsIn1)):
        
        loadInAveragedData1 = np.load(os.path.join(topDir1,dirsIn1[i], r"averagedDataVR.npy"))
        loadInErrors1 = np.load(os.path.join(topDir1,dirsIn1[i], r"errorsDataVR.npy"))
        
        loadInAveragedData2 = np.load(os.path.join(topDir2,dirsIn2[i], r"averagedDataVR.npy"))
        loadInErrors2 = np.load(os.path.join(topDir2,dirsIn2[i], r"errorsDataVR.npy"))
        
        phiValues1 = (float(dirsIn1[i])/100)*np.ones(10000)
        phiValues2 = (float(dirsIn2[i])/100)*np.ones(10000)
        
        
        #Plot Wyart-Cates Prediction for 1
        pyplot.plot(stressValues,rescaledNewtonianViscosity1(stressValues,phiValues1),label="Cornstarch Fit")
        
        if stretchedExponential ==False:
            yMinE1 = (1-(float(dirsIn1[i])/100)/((phiM1-phiM1E)+((phi01-phi01E)-(phiM1-phiM1E))*np.exp(-stressValues/(tauStar1-tauStar1E))))**(-(alpha1-alpha1E))
            yMaxE1 = (1-(float(dirsIn1[i])/100)/((phiM1+phiM1E)+((phi01+phi01E)-(phiM1+phiM1E))*np.exp(-stressValues/(tauStar1+tauStar1E))))**(-(alpha1+alpha1E))
 


        yMinE1 = np.squeeze(yMinE1)
        
        
        yMaxE1 = np.squeeze(yMaxE1)
        
        
        pyplot.fill_between(stressValues,yMinE1,yMaxE1,yMinE1>=yMaxE1)
        #Plot Wyart-Cates Prediction for 2
        pyplot.plot(stressValues,rescaledNewtonianViscosity2(stressValues,phiValues2),label="GdCl Fit")
        if stretchedExponential==False:
            yMinE2 = (1-(float(dirsIn2[i])/100)/((phiM2-phiM2E)+((phi02-phi02E)-(phiM2-phiM2E))*np.exp(-stressValues/(tauStar2-tauStar2E))))**(-(alpha2-alpha2E))
            yMaxE2 = (1-(float(dirsIn2[i])/100)/((phiM2+phiM2E)+((phi02+phi02E)-(phiM2+phiM2E))*np.exp(-stressValues/(tauStar2+tauStar2E))))**(-(alpha2+alpha2E))

            
        
        yMinE2 = np.squeeze(yMinE2)
        yMaxE2 = np.squeeze(yMaxE2)
        pyplot.fill_between(stressValues,yMinE2,yMaxE2,yMinE2>=yMaxE2)
        
        pyplot.errorbar(loadInAveragedData1[:,4],loadInAveragedData1[:,2],yerr=loadInErrors1[:,2],label=dirsIn1[i]+"C",marker='o')
        pyplot.errorbar(loadInAveragedData2[:,4],loadInAveragedData2[:,2],yerr=loadInErrors2[:,2],label=dirsIn2[i]+"GdCl",marker="D")
        
        figName = str(round(float(dirsIn1[i])))
        
        pyplot.title(figName)   
        pyplot.xlabel("Shear Stress [Pa]")
        pyplot.ylabel("Rescaled Viscosity")
        pyplot.yscale('log')
        pyplot.xscale('log')
        pyplot.ylim(1,5e3)
        pyplot.xlim(1e-1,1e3)
        pyplot.legend()
        pyplot.savefig(os.path.join(outputDir,str(float(figName))+".png"))
        pyplot.clf()
        
        

        
        