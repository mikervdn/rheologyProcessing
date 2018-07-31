#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 10:28:17 2018

@author: Mike van der Naald
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 16:28:27 2018

This package is a bunch of functions that first parses the data from the Anton-Parr
rheometer into a series of numpy arrays and saves them (rheologyDataParserOld).

Once the data is saved into a series of numpy arrays all of the arrays can be averaged together
and saved into an averaged numpy array using dataSetPackager.

All of the data 
"""

import numpy as np 
import os
from matplotlib import pyplot
from scipy.optimize import curve_fit



def rheologyDataParserOld( fileName, outputDirectory ):    
    """
    This function takes the data from the old format of the Anton-Parr rheometer and parses it into a series of text files that has
    less information but also makes the data easier to deal with.  To get a feeling
    for what this function is doing take a look at the old format and see how it is arranged first.
    
    The fileName is the fileName of the textfile that was outputted from the Anton-Parr program that is in the format 
    that is obtained when one does file-->export-->contents of active window-->tab separated text file.  The file exported
    should all be of the same packing fraction.  The output directory is where each data set is outputted as a numpy 
    saved file.  Each numpy file has a [N,6] array in it where N is the number of data points taken and the 6 correpsonds to 
    the 6 quantities measured during a constant shear rate experiment [point # , time, viscosity, shear rate, shear stress, torque].
    
    """
    
    
    #Loading in the file to parse.
    #We read the whole file as a list of lines to make the parsing easier.    
    with open(fileName, "rb") as f1:
         lines = f1.readlines()
         
    
    #We need to go through and find all of the lines that have "Data Series Information\n"
    #in them because these lines are the ones that segment our data into different experiments.
    linesWithDataSeriesInfo = np.zeros([0])
    counter = 0
    for line in lines:
        if "Data Series Information\r\n" in line.decode(encoding="utf-8", errors='ignore'):
            linesWithDataSeriesInfo= np.append(linesWithDataSeriesInfo , counter)
        counter = counter + 1
        
        
    #Now we know the total number of datasets we have is:
    totalNumOfExperiments = len(linesWithDataSeriesInfo)

    for i in range(0,totalNumOfExperiments):
        #This picks out one of the blocks of data and now we need to parse it into 
        #the separate experiments and the name of the experiement.  These will go into 
        #a text file that is in a directory that is named after the experiment 

        if i == totalNumOfExperiments-1:
            currentDataSet = lines[ int(linesWithDataSeriesInfo[i]) : int(len(lines)) ]
        else:
            currentDataSet = lines[ int(linesWithDataSeriesInfo[i]) : int(linesWithDataSeriesInfo[i+1]) ]
        
        #Now we need to pick through this data set and get all of the information we want out of it.
        #First we"ll get the name of the experiment by looping through all of the rows of the data set
        #and finding the name.
        for j in range(0, len(currentDataSet)):
            if "Name:" in currentDataSet[j].decode(encoding="utf-8", errors='ignore'):
                nameOfCurrentDataSet = currentDataSet[j][8:]
        
        nameOfCurrentDataSet = nameOfCurrentDataSet.decode(encoding="utf-8", errors='ignore').replace("\r\n","")
        #Now that we have the name of the dataset we can make a new directory 
        #inside the output directory with the name of the experiment.
        filePath = (outputDirectory+os.sep+nameOfCurrentDataSet)
        os.mkdir(filePath)
        
        
        
        
        
        #Now we want to loop through and find all the lines in our data set that have
        # "\t[s]\t[Pa·s]\t[1/s]\t[Pa]\t[mNm]\t[]\n"
        #because that is the beginning of a dataset.
        linesWhereDataStarts = np.zeros([0])
        for k in range(0, len(currentDataSet)):
            if "\t[s]\t[Pas]\t[1/s]\t[Pa]\t[mNm]\t[]\r\n" in currentDataSet[k].decode(encoding="utf-8", errors='ignore'):
                linesWhereDataStarts = np.append(linesWhereDataStarts, k)
                
                
        linesWhereDataStarts=[int(i) for i in linesWhereDataStarts]     
        #Now that we know where the data starts we can start to pick out individual datasets

        linesWhereDataEnds = np.zeros([0])
        for j in linesWhereDataStarts:
            for k in range(j,len(currentDataSet)):
                if currentDataSet[k].decode(encoding="utf-8", errors='ignore') == "\r\n":
                    linesWhereDataEnds = np.append(linesWhereDataEnds, k)
                    break
        
        linesWhereDataEnds=[int(i) for i in linesWhereDataEnds]  
        
        if len(linesWhereDataEnds)!=len(linesWhereDataStarts):
            linesWhereDataEnds = np.append(linesWhereDataEnds, len(currentDataSet) )
            
        
        for j in range(0,len(linesWhereDataEnds)):

            filePathFinal = (outputDirectory+os.sep+nameOfCurrentDataSet)
            dataSetToSave = currentDataSet[linesWhereDataStarts[j]+1:linesWhereDataEnds[j]]
            finalDataSet = np.zeros((len(dataSetToSave),6))

            for k in range(0,len(dataSetToSave)):
                finalDataSet[k,:]=np.fromstring(dataSetToSave[k].decode(encoding="utf-8", errors='ignore').replace('\t',' ').replace(',', '').replace('\n',''),sep=' ')
            finalNameForData = r"dataSet-" + str(j) + r".txt"
            filePathFinal = (outputDirectory+os.sep+nameOfCurrentDataSet+os.sep+finalNameForData)
            np.save(filePathFinal,finalDataSet)
            

def dataSetPackager(topDir):
    """
    This function takes the path topDir and finds all of the data files that end in .npy, these files should have been generated by
    the above function "oldRheologyDataParser".  That function takes the raw rheometer data and packages it into an array that has rows 
    that are arranged as [pt#, time, viscosity, shear rate, shear stress, torque].  This function loads in each individual flow curve generated
    by "oldRheologyDataParser" and averages together the xvalues (shear stress) and averages together the yvalues (viscosity) and generates
    errors for averaged data point.  
    
    YOU MUST DELETE YOUR PRESHEAR .npy FILES FROM YOUR topDir DIRECTORY OR ELSE THIS FUNCTION WILL TRY TO AVERAGE THOSE IN
    
    inputs:
        topDir is the directory in which all of the flow curves that you would like to average together.  Delete the preshears from this
        directories MAKE SURE TO PUT r BEFORE THE STRING, FOR EXAMPLE topDir = r"blah"
        
    outputs:
        It

        
    """
    
    #Intialize the containers that will hold all of the data sets.
    dataSetHolder = np.array([])
    
    #This loop will find all files in topDir that end with .npy
    for dirpath, dirnames, filenames in os.walk(topDir):
        for filename in [f for f in filenames if f.endswith(r".npy")]:
            
            if (filename == r"averagedData.npy") or (filename == r"errorsData.npy"):
                break
            #This is the current dataset that is organized as [pt#, time, viscosity, shear rate, shear stress, torque]
            currentDataSet = np.expand_dims(np.load(os.path.join(dirpath, filename)),axis=2)
            
            #The current data set should be ordered from lowest controlled shear rate to highest shear rate so we check and then flip if need be.
            if currentDataSet[0,4,0] > currentDataSet[-1,4,0]:
                currentDataSet = np.flip(currentDataSet,axis=0)
                        
            if dataSetHolder.any():
                dataSetHolder = np.append(dataSetHolder , currentDataSet , axis=2 )
            else:
                dataSetHolder = currentDataSet
                
    
    averagedData = np.mean(dataSetHolder,axis=2)
    standardDeviationData = np.std(dataSetHolder,axis=2)
    
    np.save(os.path.join(topDir, r"averagedData"), averagedData )
    np.save(os.path.join(topDir, r"errorsData"), standardDeviationData )

    return (averagedData,standardDeviationData)
            
            
    
def shearStressVsViscosityPlotter(topDir, maxMinPtsOn=False):
    
    
    #Get a list of all the subdirectories
    listOfSubDirs = next(os.walk(topDir))[1]
    
    if maxMinPtsOn==True:
        phiVsMaxViscosity = np.array([])
        phiVsMinViscosity = np.array([])
    
    for packingFraction in listOfSubDirs:
        
        
        loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedData.npy" ))
        loadInErrordata = np.load(os.path.join(topDir, packingFraction , r"errorsData.npy" ))
        
        #We want to plot the viscosity both as a function of shear stress and shear rate.
        
        #First the plot as a function of shear stress.
        pyplot.errorbar(loadInAveragedData[:,4],loadInAveragedData[:,2],xerr= loadInErrordata[:,4],yerr=loadInErrordata[:,2],label=packingFraction,marker='o')
        #Next plot the viscosity as a function of shear rate,

        if maxMinPtsOn==True:

            pyplot.plot(max(loadInAveragedData[:,4]),max(loadInAveragedData[:,2]),marker='o',markersize=12,color='red' )
            pyplot.plot(min(loadInAveragedData[:,4]),min(loadInAveragedData[:,2]),marker='o',markersize=12,color='red' )
            
            
            if phiVsMinViscosity.size == 0:
                
                phiVsMinViscosity = np.array([float(packingFraction) ,min(loadInAveragedData[:,2])])
                phiVsMaxViscosity = np.array([float(packingFraction) ,max(loadInAveragedData[:,2])])   
            else:
                
                phiVsMinViscosity = np.vstack( (phiVsMinViscosity,np.array([float(packingFraction) ,min(loadInAveragedData[:,1])])))
                phiVsMaxViscosity = np.vstack( (phiVsMaxViscosity,np.array([float(packingFraction) ,max(loadInAveragedData[:,1])])) )
                
    
    
    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.xlabel("Shear Rate [1/s]")
    pyplot.ylabel("Viscosity [Pa X s]")
    pyplot.grid(True)
    pyplot.legend()
    
    
    if maxMinPtsOn == True:
        #Sort the viscosities based no the entry in the first column.
        phiVsMinViscosity = phiVsMinViscosity[phiVsMinViscosity[:,0].argsort()]
        phiVsMaxViscosity = phiVsMaxViscosity[phiVsMaxViscosity[:,0].argsort()]
        return (phiVsMinViscosity,phiVsMaxViscosity)
        

def fitAlphaPhiJ(viscosityValues, packingFractionValues, suspendingViscosity):
    
    def divergingViscosityFunc(x,phiJ, alpha):
        
        return suspendingViscosity * ( 1 - x / phiJ )**( -alpha )
    
    
    (popt, pcov) = curve_fit(divergingViscosityFunc, packingFractionValues, viscosityValues)
    
    #alpha = popt[0]
    #phiJ = popt[1]
    
    return popt

    
    

