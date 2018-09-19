#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 16:17:35 2018

@author: mike
"""
from scipy.optimize import curve_fit
import numpy as np
import os
import plottingFunctions

def tauStarFit(packingFractionValues,stressValues,shearRateValues,phi0,phiM,alpha,suspendingViscosity):
    
    #X is a tuple that is X = (packingFractionValues,stressValues)
    def strainRateFunc(X,tauStar):
        packingFraction,shearStress = X
        
        return (1/suspendingViscosity)*shearStress*(1- packingFraction/(phiM+(phi0-phiM)*np.exp(-shearStress/tauStar)))**alpha
    
    (popt, pcov) = curve_fit(strainRateFunc, (packingFractionValues,stressValues), shearRateValues,bounds=(0, [10000]))
    
    return popt,pcov
    

def allParameterFit(packingFractionValues,stressValues,shearRateValues,suspendingViscosity):
    
    
    #X is a tuple that is X = (packingFractionValues,stressValues)
    def strainRateFunc(X,tauStar,phiM,phi0,alpha):
        packingFraction,shearStress = X
        
        return (1/suspendingViscosity)*shearStress*(1- packingFraction/(phiM+(phi0-phiM)*np.exp(-shearStress/tauStar)))**alpha
    
    (popt, pcov) = curve_fit(strainRateFunc, (packingFractionValues,stressValues), shearRateValues,bounds=(0, [100,.8,.8,5]))
    
    return popt,np.sqrt(np.diag(pcov))
    

def tauStarFitStretchExponential(packingFractionValues,stressValues,shearRateValues,phi0,phiM,alpha,suspendingViscosity):
    
    #X is a tuple that is X = (packingFractionValues,stressValues)
    def strainRateFunc(X,tauStar,beta):
        packingFraction,shearStress = X
        
        return (1/suspendingViscosity)*shearStress*(1- packingFraction/(phiM+(phi0-phiM)*np.exp((-shearStress/tauStar)**beta)))**alpha
    
    (popt, pcov) = curve_fit(strainRateFunc, (packingFractionValues,stressValues), shearRateValues,bounds=([0,.999], [10000,1]))
    
    return popt,pcov

def errorInPhiDataPoint(cornstarchMass,solventMass,cornstarchRho,solventRho,errorCornstarchMass,errorSolventMass,errorSolventRho):
    
    
    densityDenom = (solventMass/solventRho) + (cornstarchMass/cornstarchRho)
    
    dPhi_dCornstarchMass = 1/(cornstarchRho*densityDenom) - cornstarchMass/(cornstarchRho*densityDenom)**2
    
    dPhi_dSolventMass = -cornstarchMass/(cornstarchRho*solventRho*densityDenom**2)
    
    dPhi_dSolventDensity = cornstarchMass*solventMass/(cornstarchRho*(solventRho*densityDenom)**2)
    
    return np.sqrt((dPhi_dCornstarchMass*errorCornstarchMass)**2 + (dPhi_dSolventMass*errorSolventMass)**2 + (dPhi_dSolventDensity*errorSolventRho)**2)
    
    
def newtonianPlateauFinder(topDir,minViscosityList,maxViscosityList,suspendingViscosity=False,rescaleViscosity=False,tauMin=False,fixAlpha=False):
    
    listOfSubDirs = next(os.walk(topDir))[1]
    listOfSubDirs = [float(i) for i in listOfSubDirs]
    listOfSubDirs.sort()
    listOfSubDirs = [str(i) for i in listOfSubDirs]
    listOfSubDirs = listOfSubDirs[::-1]
    
    if (len(minViscosityList)!=len(listOfSubDirs)) or (len(maxViscosityList)!=len(listOfSubDirs)):
        return print("The length of minViscosityList and the length of maxViscosityList must be the same as the number of packing fractions in topDir.")
    
    phiVsMaxViscosity = np.array([])
    phiVsMinViscosity = np.array([])
    
    counter = 0
    for packingFraction in listOfSubDirs:
        
        if rescaleViscosity == False:
            loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedData.npy" ))
            loadInErrorData = np.load(os.path.join(topDir, packingFraction , r"errorsData.npy" ))
        else:
            loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedDataVR.npy" ))
            loadInErrorData = np.load(os.path.join(topDir, packingFraction , r"errorsDataVR.npy" ))
                    
        sortedViscosities = np.sort(loadInAveragedData[:,2])
                 
        if minViscosityList[counter] == -1:
            minViscosity = np.nan
            errorOnMinViscosity = np.nan
        else:
            if tauMin == False:
                minViscosity = sortedViscosities[minViscosityList[counter]]
                indexOfMinViscosity = np.where(loadInAveragedData[:,2]==minViscosity)
                errorOnMinViscosity = loadInErrorData[indexOfMinViscosity[0][0],2]
            else:
                tauMinOnset = plottingFunctions.find_nearest(loadInAveragedData[:,4],tauMin)
                minViscosity = loadInAveragedData[np.where(loadInAveragedData[:,4]==tauMinOnset),2][0]
                indexOfMinViscosity = np.where(loadInAveragedData[:,2]==minViscosity)
                errorOnMinViscosity = loadInErrorData[indexOfMinViscosity[0][0],2]
        if maxViscosityList[counter] == -1:
            maxViscosity = np.nan
            errorOnMaxViscosity = np.nan
        else:
            maxViscosity = sortedViscosities[-(maxViscosityList[counter]+1)]
            indexOfMaxViscosity = np.where(loadInAveragedData[:,2]==maxViscosity)
            errorOnMaxViscosity = loadInErrorData[indexOfMaxViscosity[0][0],2]
            
        counter = counter +1
        if phiVsMinViscosity.size == 0:
            
            phiVsMinViscosity = np.array([float(packingFraction)/100 ,minViscosity])
            phiVsMaxViscosity = np.array([float(packingFraction)/100 ,maxViscosity])
            phiVsMinViscosityErrors = np.array([float(packingFraction)/100 ,errorOnMinViscosity])
            phiVsMaxViscosityErrors = np.array([float(packingFraction)/100 ,errorOnMaxViscosity])
            
        else:
            if tauMin==False:
                phiVsMinViscosity = np.vstack( (phiVsMinViscosity,np.array([float(packingFraction)/100 ,minViscosity])))
            else:
                phiVsMinViscosity = np.vstack( (phiVsMinViscosity,np.array([float(packingFraction)/100 ,minViscosity[0]])))
            
            phiVsMaxViscosity = np.vstack( (phiVsMaxViscosity,np.array([float(packingFraction)/100 ,maxViscosity])) )
            phiVsMinViscosityErrors = np.vstack( (phiVsMinViscosityErrors,np.array([float(packingFraction)/100 ,errorOnMinViscosity])))
            phiVsMaxViscosityErrors = np.vstack( (phiVsMaxViscosityErrors,np.array([float(packingFraction)/100 ,errorOnMaxViscosity])))
    
    
    phiVsMinViscosity = phiVsMinViscosity[phiVsMinViscosity[:,0].argsort()]
    phiVsMaxViscosity = phiVsMaxViscosity[phiVsMaxViscosity[:,0].argsort()]
    phiVsMaxViscosityErrors = phiVsMaxViscosityErrors[phiVsMaxViscosityErrors[:,0].argsort()]
    phiVsMinViscosityErrors = phiVsMinViscosityErrors[phiVsMinViscosityErrors[:,0].argsort()]
        
    phiVsMinViscosity = phiVsMinViscosity[~np.isnan(phiVsMinViscosity).any(axis=1)]
    phiVsMaxViscosity = phiVsMaxViscosity[~np.isnan(phiVsMaxViscosity).any(axis=1)]
    phiVsMinViscosityErrors = phiVsMinViscosityErrors[~np.isnan(phiVsMinViscosityErrors).any(axis=1)]
    phiVsMaxViscosityErrors = phiVsMaxViscosityErrors[~np.isnan(phiVsMaxViscosityErrors).any(axis=1)]
    
    
    
    if (rescaleViscosity == False) and (fixAlpha==False):
        (phiM,alphaM,phiJErrorM,alphaErrorM) = viscosityDivergenceFit(phiVsMaxViscosity[:,0], phiVsMaxViscosity[:,1],phiVsMaxViscosityErrors[:,1], suspendingViscosity,False)
        (phi0,alpha0,phiJError0,alphaError0) = viscosityDivergenceFit(phiVsMinViscosity[:,0], phiVsMinViscosity[:,1],phiVsMinViscosityErrors[:,1], suspendingViscosity,False)
        return (phiVsMinViscosity,phiVsMaxViscosity,phiVsMinViscosityErrors,phiVsMaxViscosityErrors,phi0,alpha0,phiJError0,alphaError0,phiM,alphaM,phiJErrorM,alphaErrorM)
    if (rescaleViscosity!=False) and (fixAlpha==False):
        (phiM,alphaM,phiJErrorM,alphaErrorM) = viscosityDivergenceFit(phiVsMaxViscosity[:,0], phiVsMaxViscosity[:,1],phiVsMaxViscosityErrors[:,1], False,False)
        (phi0,alpha0,phiJError0,alphaError0) = viscosityDivergenceFit(phiVsMinViscosity[:,0], phiVsMinViscosity[:,1],phiVsMinViscosityErrors[:,1], False,False)
        return (phiVsMinViscosity,phiVsMaxViscosity,phiVsMinViscosityErrors,phiVsMaxViscosityErrors,phi0,alpha0,phiJError0,alphaError0,phiM,alphaM,phiJErrorM,alphaErrorM)
    if (rescaleViscosity != False) and (fixAlpha!=False):
        (phiM,phiJErrorM) = viscosityDivergenceFit(phiVsMaxViscosity[:,0], phiVsMaxViscosity[:,1],phiVsMaxViscosityErrors[:,1], False,fixAlpha)
        (phi0,phiJError0) = viscosityDivergenceFit(phiVsMinViscosity[:,0], phiVsMinViscosity[:,1],phiVsMinViscosityErrors[:,1], False,fixAlpha)
        return (phiVsMinViscosity,phiVsMaxViscosity,phiVsMinViscosityErrors,phiVsMaxViscosityErrors,phi0,phiJError0,phiM,phiJErrorM)
    if (rescaleViscosity == False) and (fixAlpha!=False):
        (phiM,phiJErrorM) = viscosityDivergenceFit(phiVsMaxViscosity[:,0], phiVsMaxViscosity[:,1],phiVsMaxViscosityErrors[:,1], suspendingViscosity,fixAlpha)
        (phi0,phiJError0) = viscosityDivergenceFit(phiVsMinViscosity[:,0], phiVsMinViscosity[:,1],phiVsMinViscosityErrors[:,1], suspendingViscosity,fixAlpha)
        return (phiVsMinViscosity,phiVsMaxViscosity,phiVsMinViscosityErrors,phiVsMaxViscosityErrors,phi0,phiJError0,phiM,phiJErrorM)
    

    
def viscosityDivergenceFit(packingFractions,viscosities,viscosityErrors,suspendingViscosity=False,fixAlpha=False):
    
    if (suspendingViscosity!=False) and (fixAlpha!=False):
        def divergingViscosityFunc(x,phiJ):
            return  suspendingViscosity*( 1 - (x/phiJ) )**( -fixAlpha )
        
        (phiJ, phiJError) =curve_fit(divergingViscosityFunc,packingFractions,viscosities,sigma=viscosityErrors)  
        return (phiJ, phiJError)
        
    if (suspendingViscosity!=False) and (fixAlpha==False):
        def divergingViscosityFunc(x,phiJ,alpha):
            return  suspendingViscosity*( 1 - (x/phiJ) )**( -alpha )
        (popt, pcov) =curve_fit(divergingViscosityFunc,packingFractions,viscosities,sigma=viscosityErrors)  
        errors = np.sqrt(np.diag(pcov))
        phiJ = popt[0]
        alpha = popt[1]
        phiJError = errors[0]
        alphaError = errors[1]
        return (phiJ, alpha,phiJError,alphaError)
    
    
    if (suspendingViscosity==False) and (fixAlpha==False):
        def divergingViscosityFunc(x,phiJ,alpha):
            return  ( 1 - (x/phiJ) )**( -alpha )
        (popt, pcov) =curve_fit(divergingViscosityFunc,packingFractions,viscosities,sigma=viscosityErrors)  
        errors = np.sqrt(np.diag(pcov))
        phiJ = popt[0]
        alpha = popt[1]
        phiJError = errors[0]
        alphaError = errors[1]
        return (phiJ, alpha,phiJError,alphaError)
    

    if (suspendingViscosity==False) and (fixAlpha!=False):
        def divergingViscosityFunc(x,phiJ):
            return  ( 1 - (x/phiJ) )**( -fixAlpha )
        
        (phiJ, phiJError) =curve_fit(divergingViscosityFunc,packingFractions,viscosities,sigma=viscosityErrors)  
        return (phiJ, phiJError)
        

