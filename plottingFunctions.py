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

def shearStressVsViscosityPlotter(topDir, maxMinPtsOn=False,minViscosityList=False,maxViscosityList=False,suspendingViscosity=False,solvent="normal",rescaleViscosity=False):
    
    
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
    
    counter = 0
    for packingFraction in listOfSubDirs:
        
        
        loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedData.npy" ))
        loadInErrordata = np.load(os.path.join(topDir, packingFraction , r"errorsData.npy" ))
        
        #Rescale viscosities to be dimensionless
        if rescaleViscosity == True:
            loadInAveragedData[:,2] = (1/suspendingViscosity)*loadInAveragedData[:,2]
            
                
        #First the plot as a function of shear stress.
        ax.errorbar(loadInAveragedData[:,4],loadInAveragedData[:,2],xerr= loadInErrordata[:,4],yerr=loadInErrordata[:,2],label=packingFraction,marker='o')
        #Next plot the viscosity as a function of shear rate,

        if maxMinPtsOn==True:
            sortedViscosities = np.sort(loadInAveragedData[:,2])
            
            if minViscosityList[counter] == -1:
                minViscosity = np.nan
            else:
                minViscosity = sortedViscosities[minViscosityList[counter]]
                
            if maxViscosityList[counter] == -1:
                maxViscosity = np.nan
            else:
                maxViscosity = sortedViscosities[-(maxViscosityList[counter]+1)]
                    

            if minViscosityList[counter] != -1:
                shearStressForMinViscosity = loadInAveragedData[np.where(loadInAveragedData[:,2]==minViscosity)[0][0],4]
                ax.plot(shearStressForMinViscosity,minViscosity,marker='o',markersize=10,color='blue',mfc='none' , zorder=10 )
            
            if maxViscosityList[counter] != -1:
                shearStressForMaxViscosity = loadInAveragedData[np.where(loadInAveragedData[:,2]==maxViscosity)[0][0],4]
                ax.plot(shearStressForMaxViscosity,maxViscosity,marker='o',markersize=10,color='red',mfc='none', zorder=10  )
            counter = counter +1
            
            if phiVsMinViscosity.size == 0:
                
                phiVsMinViscosity = np.array([float(packingFraction)/100 ,minViscosity])
                phiVsMaxViscosity = np.array([float(packingFraction)/100 ,maxViscosity])   
            else:
                
                phiVsMinViscosity = np.vstack( (phiVsMinViscosity,np.array([float(packingFraction)/100 ,minViscosity])))
                phiVsMaxViscosity = np.vstack( (phiVsMaxViscosity,np.array([float(packingFraction)/100 ,maxViscosity])) )


    if maxMinPtsOn == True:
        #Sort the viscosities based no the entry in the first column.
        phiVsMinViscosity = phiVsMinViscosity[phiVsMinViscosity[:,0].argsort()]
        phiVsMaxViscosity = phiVsMaxViscosity[phiVsMaxViscosity[:,0].argsort()]
        
        phiVsMinViscosity = phiVsMinViscosity[~np.isnan(phiVsMinViscosity).any(axis=1)]
        phiVsMaxViscosity = phiVsMaxViscosity[~np.isnan(phiVsMaxViscosity).any(axis=1)]
    if suspendingViscosity != False:
        if rescaleViscosity ==False:
            
            (phiM,alphaM,phiJErrorM,alphaErrorM) = fittingFunctions.fitAlphaPhiJNoRescale(phiVsMaxViscosity[:,1], phiVsMaxViscosity[:,0], suspendingViscosity)
            (phi0,alpha0,phiJError0,alphaError0) = fittingFunctions.fitAlphaPhiJNoRescale(phiVsMinViscosity[:,1], phiVsMinViscosity[:,0], suspendingViscosity)
            #textstr = '$\phi_0=%.2f$\n$\alpha_0=%.2f$\n$\phi_m=%.2f$\n$\alpha_m=%.2f$' % (phi0, alpha0, phiM,alphaM)
        else:
            (phiM,alphaM,phiJErrorM,alphaErrorM) = fittingFunctions.fitAlphaPhiJRescaled(phiVsMaxViscosity[:,1], phiVsMaxViscosity[:,0])
            (phi0,alpha0,phiJError0,alphaError0) = fittingFunctions.fitAlphaPhiJRescaled(phiVsMinViscosity[:,1], phiVsMinViscosity[:,0])
            
             
    
    
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
    
        
    if suspendingViscosity != False:
        return (phiVsMinViscosity,phiVsMaxViscosity,phiM,alphaM,phi0,alpha0,phiJError0,alphaError0,phiJErrorM,alphaErrorM)
    
    if maxMinPtsOn == True:
        return (phiVsMinViscosity,phiVsMaxViscosity)
        
    
    
def shearRateVsViscosityPlotter(topDir,suspendingViscosity=False,solvent="normal",rescaleViscosity=False):
    
    
    #Get a list of all the subdirectories
    listOfSubDirs = next(os.walk(topDir))[1]
    listOfSubDirs = [float(i) for i in listOfSubDirs]
    listOfSubDirs.sort()
    listOfSubDirs = [str(i) for i in listOfSubDirs]
    listOfSubDirs = listOfSubDirs[::-1]
    
        
    ax = pyplot.subplot(111)
    
    for packingFraction in listOfSubDirs:
        
        
        loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedData.npy" ))
        loadInErrordata = np.load(os.path.join(topDir, packingFraction , r"errorsData.npy" ))
        
        #Rescale viscosities to be dimensionless
        if rescaleViscosity == True:
            loadInAveragedData[:,2] = (1/suspendingViscosity)*loadInAveragedData[:,2]
            
                
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
    

        

def viscosityDivergencePlot(lowViscosityValues,lowViscosityValuesGdCl, packingFractionValuesGdCl,packingFractionValuesCornstarch, suspendingViscosityCornstarch,suspendingViscosityGdCl):
    
    
    def divergingViscosityFuncGdCl(x,phiJ, alpha):
        return suspendingViscosityGdCl * ( 1 - (x/phiJ) )**( -alpha )
    
    def divergingViscosityFuncNormal(x,phiJ, alpha):
        return suspendingViscosityCornstarch * ( 1 - (x/phiJ) )**( -alpha )
    
    (phiJGdCl,alphaGdCl,phiJErrorG,alphaErrorG) = fittingFunctions.fitAlphaPhiJRescaled(lowViscosityValuesGdCl, packingFractionValuesGdCl)
    (phiJ0,alpha0,phiJError0,alphaError0) = fittingFunctions.fitAlphaPhiJRescaled(lowViscosityValues, packingFractionValuesCornstarch)
    
#    phiJLineYValues = np.linspace(0,max(highViscosityValues),50)
#    phiJMLineXValues = phiJGdCl*np.ones((50))
#    phiJ0LineXValues = phiJ0*np.ones((50))
    
    
    packingFractionValuesGdClXValues = np.linspace(0,.37,100)
    packingFractionValuesCornstarchXValues = np.linspace(0,.37,100)
    
    
    pyplot.plot(packingFractionValuesGdCl,lowViscosityValuesGdCl,'o',color='C0')
    pyplot.plot(packingFractionValuesCornstarch,lowViscosityValues,'o',color='C1')
    pyplot.plot(packingFractionValuesGdClXValues,divergingViscosityFuncGdCl(packingFractionValuesGdClXValues,phiJGdCl,alphaGdCl),color='C0')
    pyplot.plot(packingFractionValuesCornstarchXValues,divergingViscosityFuncNormal(packingFractionValuesCornstarchXValues,phiJ0,alpha0),color='C1')
#    pyplot.plot(phiJMLineXValues,phiJLineYValues,'.')
#    pyplot.plot(phiJ0LineXValues,phiJLineYValues,'.')
    pyplot.legend(("GdCl Cornstarch Suspension","Normal Cornstarch Suspension","GdCl Fit","Normal Cornstarch Fit"))
    pyplot.title("Frictional Newtonian Viscosities Vs Packing Fraction")
    pyplot.xlabel("Packing Fraction $\phi$")
    pyplot.ylabel("Viscosity [Pa s]")
    

def logarithmicPlotter(viscosityValues0,viscosityValues1,packingFractionValues0,packingFractionValues1,suspendingViscosity0,alpha0,phiJ0,suspendingViscosity1,alpha1,phiJ1,highOrLow):
    
    #logEta0 = np.log10(suspendingViscosity0)*np.ones(len(viscosityValues0))
    logViscosityValues0 = np.log10(viscosityValues0)
    logPackingFraction0 = np.log10(1-packingFractionValues0/phiJ0)
    
    #logEta1 = np.log10(suspendingViscosity1)*np.ones(len(viscosityValues1))
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
    
    
#    logEta0C_0 = np.log10(suspendViscosityC)*np.ones(len(phiVsMinViscosityC[:,1]))
#    logEta0G_0 = np.log10(suspendViscosityG)*np.ones(len(phiVsMinViscosityG[:,1]))
#    
#    
#    logEta0C_M = np.log10(suspendViscosityC)*np.ones(len(phiVsMaxViscosityC[:,1]))
#    logEta0G_M = np.log10(suspendViscosityG)*np.ones(len(phiVsMaxViscosityG[:,1]))    
#    
    
    
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



def sideBySideCurvePlot(topDir1,topDir2,outputDir,baseViscosity1,baseViscosity2,rescaledViscosity=False,allPlots=False,evenCurves="none"):
    
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
                

            
        
        #if i%2==0:
        #   continue

        dir1Data = np.load(os.path.join(topDir1,dirsIn1[i], r"averagedData.npy"))
        dir1Errors = np.load(os.path.join(topDir1,dirsIn1[i], r"errorsData.npy"))
        
        dir2Data = np.load(os.path.join(topDir2,dirsIn2[i], r"averagedData.npy"))
        dir2Errors = np.load(os.path.join(topDir2,dirsIn2[i], r"errorsData.npy"))
        
        if rescaledViscosity == True:
            dir1Data[:,2] = dir1Data[:,2]/baseViscosity1
            dir2Data[:,2] = dir2Data[:,2]/baseViscosity2
            
        
        pyplot.errorbar(dir1Data[:,4],dir1Data[:,2],xerr= dir1Errors[:,4],yerr=dir1Errors[:,2],label=dirsIn1[i]+"C",marker='o')
        pyplot.errorbar(dir2Data[:,4],dir2Data[:,2],xerr= dir2Errors[:,4],yerr=dir2Errors[:,2],label=dirsIn2[i]+"GdCl",marker='1')
        
        
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
    pyplot.fill_betweenx(stressValues,phiValuesJamming,phiValuesBeyondJamming,where=phiValuesJamming <= phiValuesBeyondJamming, facecolor='beige', interpolate=True)
        
    
    pyplot.ylim(10**-1, 10000)
    pyplot.xlim(phiM-.1,phi0+.05)
    
    pyplot.xlabel(r"Packing Fraction $\phi$")
    pyplot.ylabel(r"Shear Stress $ \tau $ [Pa]")
    pyplot.title("Shear Jamming Phase Diagram")
    
    
    
    
    
    
    
    
    
    

    
    
        
        
