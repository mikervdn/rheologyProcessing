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

def shearStressVsViscosityPlotter(topDir, maxMinPtsOn=False,minViscosityList=False,maxViscosityList=False,solvent="normal",rescaleViscosity=True,tauMin=False,noErrors=False,suspendingViscosity=False):
    
    
    #Get a list of all the subdirectories
    listOfSubDirs = next(os.walk(topDir))[1]
    listOfSubDirs = [float(i) for i in listOfSubDirs]
    listOfSubDirs.sort()
    listOfSubDirs = [str(i) for i in listOfSubDirs]
    listOfSubDirs = listOfSubDirs[::-1]
    
    if maxViscosityList != False:
        if (len(minViscosityList)!=len(listOfSubDirs)) or (len(maxViscosityList)!=len(listOfSubDirs)):
            return print("The length of minViscosityList and the length of maxViscosityList must be the same as the number of packing fractions in topDir.")
        


    if maxMinPtsOn==True:
        phiVsMaxViscosity = np.array([])
        phiVsMinViscosity = np.array([])
        
        
    ax = pyplot.subplot(111)
    
    for packingFraction in listOfSubDirs:
        
        if rescaleViscosity == False:
            loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedData.npy" ))
        else:
            loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedDataVR.npy" ))
            
            
        if noErrors == False:
            if rescaleViscosity ==False:
                loadInErrordata = np.load(os.path.join(topDir, packingFraction , r"errorsData.npy" ))
            else:
                loadInErrordata = np.load(os.path.join(topDir, packingFraction , r"errorsDataVR.npy" ))
        
                
        #First the plot as a function of shear stress.
        if noErrors == False:
            ax.errorbar(loadInAveragedData[:,4],loadInAveragedData[:,2],xerr= loadInErrordata[:,4],yerr=loadInErrordata[:,2],label=packingFraction,marker='o',markersize=5)
        else:
            ax.plot(loadInAveragedData[:,4],loadInAveragedData[:,2],label=packingFraction,marker='o')
        #Next plot the viscosity as a function of shear rate,




    #THIS IS THE ONLY THING THAT SHOULD REMAIN OF maxMinPtsOn
    if maxMinPtsOn==True:
        if tauMin==False:
            if (suspendingViscosity==False) and (rescaleViscosity!=False):
                (phiVsMinViscosity,phiVsMaxViscosity,phiM,alphaM,phi0,alpha0,phiJError0,alphaError0,phiJErrorM,alphaErrorM)=fittingFunctions.newtonianPlateauFinder(topDir,minViscosityList,maxViscosityList,False,True,False)
            else:
                (phiVsMinViscosity,phiVsMaxViscosity,phiM,alphaM,phi0,alpha0,phiJError0,alphaError0,phiJErrorM,alphaErrorM)=fittingFunctions.newtonianPlateauFinder(topDir,minViscosityList,maxViscosityList,suspendingViscosity,False,False)
        else:
            if (suspendingViscosity==False) and (rescaleViscosity!=False):
                (phiVsMinViscosity,phiVsMaxViscosity,phiM,alphaM,phi0,alpha0,phiJError0,alphaError0,phiJErrorM,alphaErrorM)=fittingFunctions.newtonianPlateauFinder(topDir,minViscosityList,maxViscosityList,False,True,tauMin)
            else:
                (phiVsMinViscosity,phiVsMaxViscosity,phiM,alphaM,phi0,alpha0,phiJError0,alphaError0,phiJErrorM,alphaErrorM)=fittingFunctions.newtonianPlateauFinder(topDir,minViscosityList,maxViscosityList,suspendingViscosity,False,tauMin)
                

        for packingFraction in listOfSubDirs:
            
            if rescaleViscosity == False:
                loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedData.npy" ))
            else:
                loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedDataVR.npy" ))
                

            if float(packingFraction)/100 in phiVsMinViscosity:
                (x,y) = np.where(phiVsMinViscosity==float(packingFraction)/100)
                minViscosity = phiVsMinViscosity[x,1]
     
                if tauMin==False:
                    shearStressForMinViscosity = loadInAveragedData[np.where(loadInAveragedData[:,2]==minViscosity)[0][0],4]
                    ax.plot(shearStressForMinViscosity,minViscosity,marker='o',markersize=10,color='blue',mfc='none' , zorder=10 )
                else:
                    tauMinOnset = find_nearest(loadInAveragedData[:,4],tauMin)
                    ax.plot(tauMinOnset,minViscosity[0],marker='o',markersize=10,color='blue',mfc='none' , zorder=10 )

            if float(packingFraction)/100 in phiVsMaxViscosity:
                (x,y) = np.where(phiVsMaxViscosity==float(packingFraction)/100)
                maxViscosity = phiVsMaxViscosity[x,1]
                
                shearStressForMaxViscosity = loadInAveragedData[np.where(loadInAveragedData[:,2]==maxViscosity)[0][0],4]
                ax.plot(shearStressForMaxViscosity,maxViscosity,marker='o',markersize=10,color='red',mfc='none', zorder=10  )
            
            
        ax.text(1.3, 0.9, r"$\phi_0=$ "+ str(round(phi0,4)) + " $\pm$ " + str(round(phiJError0,4)), transform=ax.transAxes, fontsize=14,verticalalignment='top')
        ax.text(1.3, 0.83, r"$\phi_M=$ "+ str(round(phiM,4)) + " $\pm$ " + str(round(phiJErrorM,4)), transform=ax.transAxes, fontsize=14,verticalalignment='top')
        ax.text(1.3, 0.76, r"$\alpha_0=$ "+ str(round(alpha0,4)) + " $\pm$ " + str(round(alphaError0,4)), transform=ax.transAxes, fontsize=14,verticalalignment='top')
        ax.text(1.3, 0.69, r"$\alpha_M=$ "+ str(round(alphaM,4)) + " $\pm$ " + str(round(alphaErrorM,4)), transform=ax.transAxes, fontsize=14,verticalalignment='top')
        
            
        


            

    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.xlabel("Shear Stress [Pa]")
    
    if rescaleViscosity==False:
        pyplot.ylabel("Viscosity [Pa s]")
    else:
        pyplot.ylabel(r"Rescaled Viscosity $\frac{\eta}{\eta_0}$")
        
    if solvent == "GdCl":
        pyplot.title("Flow Curves for GdCl Cornstarch Suspensions")
    if solvent =="normal":
        pyplot.title("Flow Curves for Normal Cornstarch Suspensions")
    if solvent =="NH4":
        pyplot.title("Flow Curves for $(NH_4)_2 SO_4$ Cornstarch Suspensions")
        
    pyplot.grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    #if maxMinPtsOn == True:
        #pyplot.text(1700, 10.2, "$\phi_0$ = " + str(round(phi0,4)) + "$\pm$" + str(round(phiJError0,4)),fontsize=14)
        #pyplot.text(1700, 7.8, "$\phi_M$ = " + str(round(phiM,4))+ "$\pm$" + str(round(phiJErrorM,4)),fontsize=14)
        #pyplot.text(1700, .0035, r"$\alpha_M$ = " + str(round(alphaM,4))+ "$\pm$" + str(round(alphaErrorM,4)),fontsize=14)
        #pyplot.text(1700, .0065, r"$\alpha_0$ = " + str(round(alpha0,4))+ "$\pm$" + str(round(alphaError0,4)),fontsize=14)
#        ax.text(2, 1, "$\phi_0$ = " + str(round(phi0,2)))
#        ax.text(1600, 4, "$alpha_0$ = " + str(round(alpha0,2)))
        #ax.text(1600, 0, "$\phi_M$ = " + str(round(phiM,2)))
        #ax.text(1600, 0, "$alpha_M$ = " + str(round(alphaM,2)))
    # Put a legend to the right of the current axis
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    if maxMinPtsOn == True:
        if rescaleViscosity != False:
            return (phiVsMinViscosity,phiVsMaxViscosity,phiM,alphaM,phi0,alpha0,phiJError0,alphaError0,phiJErrorM,alphaErrorM)
        else:
            return (phiVsMinViscosity,phiVsMaxViscosity)
    

        
    
    
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
    

        

def viscosityDivergencePlot(viscosityValues,packingFractionValues,labelsForCurves,fixAlpha=False,titleOfPlot ="Newtonian Plateau's Divergences"):
    
    (numSystems,_) = np.shape(viscosityValues)
    
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
        currentViscosityValues = viscosityValues[i,:]
        currentViscosityValues = currentViscosityValues[~np.isnan(currentViscosityValues)]
        
        currentPhiValues = packingFractionValues[i,:]
        currentPhiValues = currentPhiValues[~np.isnan(currentPhiValues)]
        if fixAlpha == False:
            (phiJHolder[i],alphaHolder[i],phiJErrorHolder[i],alphaErrorHolder[i]) = fittingFunctions.fitAlphaPhiJRescaled(currentViscosityValues, currentPhiValues)
            ax.plot(np.linspace(0,phiJHolder[i]-.001,1000),divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i],alphaHolder[i]),label=labelsForCurves[i] +" Fit" )

            
            minimumError = divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]-phiJErrorHolder[i],alphaHolder[i]+alphaErrorHolder[i])
            maximumError = divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]+phiJErrorHolder[i],alphaHolder[i]-alphaErrorHolder[i])
            #pyplot.plot(np.linspace(0,phiJHolder[i]-.001,1000),minimumError,label="minimumError")
            #pyplot.plot(np.linspace(0,phiJHolder[i]-.001,1000), maximumError,label="maximumError")
            
            ax.fill_between(np.linspace(0,phiJHolder[i]-.001,1000), minimumError, maximumError, where=maximumError <= minimumError, facecolor='lightblue')
        else:
            (phiJHolder[i],phiJErrorHolder[i]) = fittingFunctions.fitPhiJRescaledFixedAlpha(currentViscosityValues, currentPhiValues,fixAlpha)
            minimumError = divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]-phiJErrorHolder[i])
            maximumError = divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]+phiJErrorHolder[i])
            ax.fill_between(np.linspace(0,phiJHolder[i]-.001,1000), minimumError, maximumError, where=maximumError <= minimumError, facecolor='lightblue')
            
            ax.plot(np.linspace(0,phiJHolder[i]-.001,1000),divergingViscosityFunc(np.linspace(0,phiJHolder[i]-.001,1000),phiJHolder[i]),label=labelsForCurves[i] +" Fit" )
            
        ax.plot(currentPhiValues, currentViscosityValues,'o',label=labelsForCurves[i])

        
        
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



def sideBySideCurvePlot(topDir1,topDir2,outputDir,rescaledViscosity=False,allPlots=False,evenCurves="none"):
    
    #get list of all of the subdirectories in topDir1 and topDir2
    dirsIn1 = next(os.walk(topDir1))[1]
    dirsIn2 = next(os.walk(topDir2))[1]
    
    #They need to be sorted so we have to convert them from strings to floats and then back
    dirsIn1 = np.sort([float(i) for i in dirsIn1])
    dirsIn2 = np.sort([float(i) for i in dirsIn2])
    #Back to strings!
    dirsIn1 = [str(i) for i in dirsIn1]
    dirsIn2 = [str(i) for i in dirsIn2]
    
    pyplot.figure(figsize=(14.0, 5.0))


    
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
        else:
            dir1Data = np.load(os.path.join(topDir1,dirsIn1[i], r"averagedData.npy"))
            dir1Errors = np.load(os.path.join(topDir1,dirsIn1[i], r"errorsData.npy"))
            dir2Data = np.load(os.path.join(topDir2,dirsIn2[i], r"averagedData.npy"))
            dir2Errors = np.load(os.path.join(topDir2,dirsIn2[i], r"errorsData.npy"))
        
        
        pyplot.errorbar(dir1Data[:,4],dir1Data[:,2],yerr=dir1Errors[:,2],label=dirsIn1[i]+"C",marker='o')
        pyplot.errorbar(dir2Data[:,4],dir2Data[:,2],yerr=dir2Errors[:,2],label=dirsIn2[i]+"GdCl",marker="D")
        
        
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



def shearJammingPhaseDiagramTwoSystems(phi01,phiM1,tauStar1,alpha1,phi02,phiM2,tauStar2,alpha2,DST=False,SJ=False):
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
        pyplot.semilogy(phiValuesDST1,stressValues,label = 'DST-1',color="black")
        #Plot all four curves, the DST curve, the SJ curve, and the two vertical lines for the first system.
        pyplot.semilogy(phiValuesDST2,stressValues,label = 'DST-2',color="black")
        pyplot.fill_betweenx(stressValues,phiValuesDST1,phiValuesSJ1,where=phiValuesDST1 <= phiValuesSJ1, facecolor='lightblue', interpolate=True)
        pyplot.fill_betweenx(stressValues,phiValuesDST2,phiValuesSJ2,where=phiValuesDST2 <= phiValuesSJ2, facecolor='lightgreen', interpolate=True)
        pyplot.semilogy(phiValuesSJ1,stressValues,label = 'SJ-1',color="black")
        pyplot.semilogy(phiValuesSJ2,stressValues,label = 'SJ-2',color="black")
    if SJ == True:
        pyplot.semilogy(phiValuesSJ1,stressValues,label = 'SJ-1',color="black")
        pyplot.semilogy(phiValuesSJ2,stressValues,label = 'SJ-2',color="black")
        pyplot.semilogy(phiValuesJamming1,stressValues,label = 'J-1',color="black")
        pyplot.semilogy(phiValuesJamming2,stressValues,label = 'J-2',color="black")
        pyplot.fill_betweenx(stressValues,phiValuesSJ1,phiValuesJamming1,where=phiValuesSJ1 <= phiValuesJamming1, facecolor='lightgreen', interpolate=True)
        pyplot.fill_betweenx(stressValues,phiValuesSJ2,phiValuesJamming2,where=phiValuesSJ2 <= phiValuesJamming2, facecolor='lightgreen', interpolate=True)
        
        
                                                                                                                                            
    


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
        
        
def WCModelComparisonTwoSystems(topDir1,topDir2,outputDir,phi01,phiM1,tauStar1,alpha1,phi02,phiM2,tauStar2,alpha2):
    
    
    
    dirsIn1 = next(os.walk(topDir1))[1]
    dirsIn2 = next(os.walk(topDir2))[1]
    
    #They need to be sorted so we have to convert them from strings to floats and then back
    dirsIn1 = np.sort([float(i) for i in dirsIn1])
    dirsIn2 = np.sort([float(i) for i in dirsIn2])
    #Back to strings!
    dirsIn1 = [str(i) for i in dirsIn1]
    dirsIn2 = [str(i) for i in dirsIn2]

    
    stressValues = np.linspace(.0001,1e3,10000)
    

    def rescaledNewtonianViscosity1(stress,phi):
        return (1-phi/(phiM1+(phi01-phiM1)*np.exp(-stress/tauStar1)))**(-alpha1)
    
    def rescaledNewtonianViscosity2(stress,phi):
        return (1-phi/(phiM2+(phi02-phiM2)*np.exp(-stress/tauStar2)))**(-alpha2)

    
    for i in range(0,len(dirsIn1)):
        
        loadInAveragedData1 = np.load(os.path.join(topDir1,dirsIn1[i], r"averagedDataVR.npy"))
        loadInErrors1 = np.load(os.path.join(topDir1,dirsIn1[i], r"errorsDataVR.npy"))
        
        loadInAveragedData2 = np.load(os.path.join(topDir2,dirsIn2[i], r"averagedDataVR.npy"))
        loadInErrors2 = np.load(os.path.join(topDir2,dirsIn2[i], r"errorsDataVR.npy"))
        
        phiValues1 = (float(dirsIn1[i])/100)*np.ones(10000)
        phiValues2 = (float(dirsIn2[i])/100)*np.ones(10000)
        
        
        #Plot Wyart-Cates Prediction for 1
        pyplot.plot(stressValues,rescaledNewtonianViscosity1(stressValues,phiValues1),label="Cornstarch Fit")
        #Plot Wyart-Cates Prediction for 2
        pyplot.plot(stressValues,rescaledNewtonianViscosity2(stressValues,phiValues2),label="GdCl Fit")
        
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
        
        

        
        