# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 16:28:27 2018

This package is a bunch of functions that generates a shear-jamming phase diagram.
"""

import numpy as np 
import os
from matplotlib import pyplot
from scipy.optimize import curve_fit



def oldRheologyDataParser(fileName,outputDirectory):    
    """
    This function takes the data from the old format of the Anton-Parr rheometer and parses it into a series of text files that has
    less information but also makes the data easier to deal with.  To get a feeling
    for what this function is doing take a look at the old format and see how it is arranged first.
    
    This package is a bunch of functions that generates a shear-jamming phase diagram.
    """
    
    
    #Loading in the file to parse.
    #We read the whole file as a list of lines to make the parsing easier.    
    with open(fileName, "r") as f1:
         lines = f1.readlines()
         
    
    #We need to go through and find all of the lines that have "Data Series Information\n"
    #in them because these lines are the ones that segment our data into different experiments.
    linesWithDataSeriesInfo = np.zeros([0])
    counter = 0
    for line in lines:
        if "Data Series Information\n" in line:
            linesWithDataSeriesInfo= np.append(linesWithDataSeriesInfo , counter)
        counter = counter + 1
        
        
    #Now we know the total number of datasets we have is:
    totalNumOfExperiments = len(linesWithDataSeriesInfo)

    for i in range(0,totalNumOfExperiments):
        print("i equals " + str(i))
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
            if "Name:" in currentDataSet[j]:
                nameOfCurrentDataSet = currentDataSet[j][8:].strip("\n").replace(" ", "")
        
        #Now that we have the name of the dataset we can make a new directory 
        #inside the output directory with the name of the experiment.
        filePath = (outputDirectory+os.sep+nameOfCurrentDataSet)
        os.mkdir(filePath)
        
        
        
        
        
        #Now we want to loop through and find all the lines in our data set that have
        # "\t[s]\t[Pa·s]\t[1/s]\t[Pa]\t[mNm]\t[]\n"
        #because that is the beginning of a dataset.
        linesWhereDataStarts = np.zeros([0])
        for k in range(0, len(currentDataSet)):
            if "\t[s]\t[Pa·s]\t[1/s]\t[Pa]\t[mNm]\t[]\n" in currentDataSet[k]:
                linesWhereDataStarts = np.append(linesWhereDataStarts, k)
                
                
        linesWhereDataStarts=[int(i) for i in linesWhereDataStarts]     
        #Now that we know where the data starts we can start to pick out individual datasets

        linesWhereDataEnds = np.zeros([0])
        for j in linesWhereDataStarts:
            for k in range(j,len(currentDataSet)):
                if currentDataSet[k]=="\n":
                    linesWhereDataEnds = np.append(linesWhereDataEnds, k)
                    break
                
        if len(linesWhereDataEnds)!=len(linesWhereDataStarts):
            linesWhereDataEnds = np.append(linesWhereDataEnds, len(currentDataSet) )
                                
        linesWhereDataEnds=[int(i) for i in linesWhereDataEnds]  
        
        for j in range(0,len(linesWhereDataEnds)):

            filePathFinal = (outputDirectory+os.sep+nameOfCurrentDataSet)
            dataSetToSave = currentDataSet[linesWhereDataStarts[j]+1:linesWhereDataEnds[j]]
            finalDataSet = np.zeros((len(dataSetToSave),6))

            for k in range(0,len(dataSetToSave)):
                finalDataSet[k,:]=np.fromstring(dataSetToSave[k].replace('\t',' ').replace(',', '').replace('\n',''),sep=' ')
            finalNameForData = r"dataSet-" + str(j) + r".txt"
            filePathFinal = (outputDirectory+os.sep+nameOfCurrentDataSet+os.sep+finalNameForData)
            np.save(filePathFinal,finalDataSet)
            
            
    
def dataSetAverager(x,y):
    """This function takes in multiple data sets and averages them and returns their error as well.  
    x should be of the shape (N,numPoints) where N is the number of experiments done and numPoints is the data points collected
    
    """

    xAverage = np.mean(x,axis=0)    
    yAverage = np.mean(y,axis=0) 
    
    xError = np.std(x,axis=0)
    yError = np.std(y,axis=0)
    
    return (xAverage,yAverage,xError,yError)


def stressVsViscosityDataSetPackager(topDir):
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
        averagedShearSress the averaged shear stress from all the flow curves found in topDir 
        
        averagedViscosity the averaged viscosity from all the flow curves found in topDir 
        
        errorsInShearSress is the standard deviation of the shear stresses averaged from the flow curves found in topDir 
        
        errorsInViscosity is the standard deviation of the averaged viscosities from the flow curves found in topDir
        
    """
    
    #Intialize the containers that will hold all of the data sets.
    shearStressHolder = []
    viscosityHolder = []
    
    #This loop will find all files in topDir that end with .npy
    for dirpath, dirnames, filenames in os.walk(topDir):
        for filename in [f for f in filenames if f.endswith(r".npy")]:
            
            if filename == r"averagedShearStressVsViscosity.npy":
                break
            #This is the current dataset that is organized as [pt#, time, viscosity, shear rate, shear stress, torque]
            currentDataSet = np.load(os.path.join(dirpath, filename))
            
            #The current data set should be ordered from lowest controlled shear rate to highest shear rate so we check and then flip if need be.
            if currentDataSet[0,4] > currentDataSet[-1,4]:
                currentDataSet = np.flip(currentDataSet,axis=0)
                        
            if viscosityHolder == []:
                viscosityHolder = np.expand_dims(currentDataSet[:,2],axis=1)
            else:
                viscosityHolder = np.append(viscosityHolder , np.expand_dims(currentDataSet[:,2],axis=1) , axis=1 )
                 
                
            if shearStressHolder == []:
                shearStressHolder = np.expand_dims(currentDataSet[:,4],axis=1)
            else:
                shearStressHolder = np.append(shearStressHolder , np.expand_dims(currentDataSet[:,4],axis=1) , axis=1 )
                

    
    
    
    (averagedShearSress, averagedViscosity, errorsInShearSress,errorsInViscosity) = dataSetAverager( np.transpose(shearStressHolder) , np.transpose(viscosityHolder) )
    
    np.save(os.path.join(topDir, r"averagedShearStressVsViscosity"), np.concatenate( (np.expand_dims(averagedShearSress,axis=1) ,np.expand_dims(averagedViscosity,axis=1),np.expand_dims(errorsInShearSress,axis=1),np.expand_dims(errorsInViscosity,axis=1)), axis=1 ) )

    return (averagedShearSress, averagedViscosity, errorsInShearSress,errorsInViscosity)
            
            
    
def shearStressVsViscosityPlotter(topDir, maxMinPtsOn=False):
    
    
    #Get a list of all the subdirectories
    listOfSubDirs = next(os.walk(topDir))[1]
    
    if maxMinPtsOn==True:
        phiVsMaxViscosity = np.array([])
        phiVsMinViscosity = np.array([])
    
    for packingFraction in listOfSubDirs:
        
        loadInAveragedData = np.load(os.path.join(topDir, packingFraction , r"averagedShearStressVsViscosity.npy" ))
        
        pyplot.errorbar(loadInAveragedData[:,0],loadInAveragedData[:,1],loadInAveragedData[:,2],loadInAveragedData[:,3],label=packingFraction,marker='o')
        
        if maxMinPtsOn==True:

            pyplot.plot(max(loadInAveragedData[:,0]),max(loadInAveragedData[:,1]),marker='o',markersize=12,color='red' )
            pyplot.plot(min(loadInAveragedData[:,0]),min(loadInAveragedData[:,1]),marker='o',markersize=12,color='red' )
            
            
            if phiVsMinViscosity.size == 0:
                phiVsMinViscosity = np.array([float(packingFraction) ,min(loadInAveragedData[:,1])])
                phiVsMaxViscosity = np.array([float(packingFraction) ,max(loadInAveragedData[:,1])])   
            else:
                print(phiVsMinViscosity)
                phiVsMinViscosity = np.vstack( (phiVsMinViscosity,np.array([float(packingFraction) ,min(loadInAveragedData[:,1])])))
                phiVsMaxViscosity = np.vstack( (phiVsMaxViscosity,np.array([float(packingFraction) ,max(loadInAveragedData[:,1])])) )
                
    
    
    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.xlabel("Shear Stress [Pa]")
    pyplot.ylabel("Viscosity [Pa X s]")
    pyplot.grid(True)
    pyplot.legend()
    if maxMinPtsOn == True:
        return (phiVsMinViscosity,phiVsMaxViscosity)
        
    
    
        


def fitAlphaPhiJ(viscosityValues, packingFractionValues, suspendingViscosity):
    
    def divergingViscosityFunc(x,phiJ, alpha):
        
        return suspendingViscosity * ( 1 - x / phiJ )**(-alpha)
    
    
    (popt, pcov) = curve_fit(divergingViscosityFunc, packingFractionValues, viscosityValues)
    
    #alpha = popt[0]
    #phiJ = popt[1]
    
    return popt

    
    

