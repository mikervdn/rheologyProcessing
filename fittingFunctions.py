#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 16:17:35 2018

@author: mike
"""
from scipy.optimize import curve_fit
import numpy as np


def fitAlphaPhiJNoRescale(viscosityValues, packingFractionValues, suspendingViscosity):
    
    def divergingViscosityFunc(x,phiJ, alpha):
        
        return suspendingViscosity * ( 1 - (x/phiJ) )**( -alpha )
    
    
    (popt, pcov) = curve_fit(divergingViscosityFunc, packingFractionValues, viscosityValues,bounds=(0, [1, 100]))
    
    alpha = popt[1]
    phiJ = popt[0]
    errors = np.sqrt(np.diag(pcov))
    alphaError = errors[1]
    phiJError = errors[0]
    
    
    return (phiJ,alpha,phiJError,alphaError)

def fitAlphaPhiJRescaled(viscosityValues, packingFractionValues):
    
    def divergingViscosityFunc(x,phiJ, alpha):
        
        return  ( 1 - (x/phiJ) )**( -alpha )
    
    
    (popt, pcov) = curve_fit(divergingViscosityFunc, packingFractionValues, viscosityValues,bounds=(0, [1, 100]))
    
    alpha = popt[1]
    phiJ = popt[0]
    errors = np.sqrt(np.diag(pcov))
    alphaError = errors[1]
    phiJError = errors[0]
    
    
    return (phiJ,alpha,phiJError,alphaError)

def tauStarFit(packingFractionValues,stressValues,shearRateValues,phi0,phiM,alpha,suspendingViscosity):
    
    #X is a tuple that is X = (packingFractionValues,stressValues)
    def strainRateFunc(X,tauStar):
        packingFraction,shearStress = X
        
        return (1/suspendingViscosity)*shearStress*(1- packingFraction/(phiM+(phi0-phiM)*np.exp(shearStress/tauStar)))**alpha
    
    (popt, pcov) = curve_fit(strainRateFunc, (packingFractionValues,stressValues), stressValues,bounds=(0, [10000]))
    