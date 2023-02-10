#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:16:34 2023

@author: Aaron Mooney
"""

# References:
# Sauro, J., &; Lewis, J. R. (2016). 
# In Quantifying the user experience: Practical statistics for user research 
# (pp. 19–60). essay, Morgan Kaufmann. 

import pandas as pd
import numpy as np
import statistics as stats
import scipy.stats
from scipy.stats import sem
import scipy.stats
import math

# *** Descriptive Statistics Block ***

def getNumScores(scores):
    """Reads in a list of scores. Calculates the 'n' value.
    Returns the 'n' value as numScores."""
    
    numScores = len(scores)
    
    return numScores

def getSumScores(scores):
    """Reads in a list of scores. Calculates the sum of the values.
    Returns the sum as sumScores."""

    sumScores = sumScores = round((sum(scores)), 4)
    
    return sumScores

def getRangeVals(scores):
    """Reads in a list of scores. Calculates the min, max, and range
    values. 
    Returns the min, max, and range values as a list object."""
    
    minVal =  minVal = round((min(scores)), 4)
    maxVal = round((max(scores)), 4)
    rangeVal = round((maxVal - minVal), 4)
     
    rangeValList = [minVal, maxVal, rangeVal]
    
    return rangeValList

def getMeanVal(scores):
    """Reads in a list of scores. Calculates the mean value.
    Returns the mean value as meanScores."""
    
    meanScores = round(stats.mean(scores), 4)
    
    return meanScores

def getStdDev(scores):
    """Reads in a list of scores. Calculates the stdDev value.
    Returns the stdDev value."""
    
    stdDev = round(stats.stdev(scores), 4)
    
    return stdDev

def getStdErr(scores):
    """Reads in a list of scores. Calculates the stdErr value.
    Returns the stdErr value."""
    
    stdErr = round(sem(scores), 4)
    
    return stdErr

def getMedian(scores):
    """Reads in a list of scores. Calculates the median value.
    Returns the median value."""
    
    medianVal = stats.median(scores)
    
    return medianVal

def getDescriptStats(scores):
    """Reads in a list of scores. Calls the functions necessary to calculate
    values for descriptive statistics.
    Returns list object containing the descriptive statistics values.""" 
    
    numScores = getNumScores(scores)
    sumScores = getSumScores(scores)
    rangeValList = getRangeVals(scores)
    meanVal = getMeanVal(scores)
    medianVal = getMedian(scores)
    stdDev = getStdDev(scores)
    stdErr = getStdErr(scores)
    
    descriptStats = [numScores, sumScores, meanVal, medianVal, stdDev, stdErr]
    
    return descriptStats

def reportDescriptStats(descriptStats):
    """Reads in the list object generated from getDescripStats function.
    Breaks the list object value into variables, prints the variable label
    along with the variable value.
    Returns None."""
    
    numScores = descriptStats[0]
    sumScores = descriptStats[1]
    
    #minVal = descriptStats[2][0]
    #maxVal = descriptStats[2][1]
    #rangeVals = descriptStats[2][2]
        
    meanVal = round((descriptStats[2]), 2)
    medianVal = round((descriptStats[3]), 2)
        
    rawStdDev = descriptStats[4]
    stdDev = round((rawStdDev), 2)
    
    rawStdErr = descriptStats[5]
    stdErr = round((rawStdErr), 2)
        
    print('\n*** Descriptive Statistics ***\n')
    print('n:', numScores)
    print('')
        
    # print('Range:', rangeVals)
    # print('Min Value:', minVal)
    # print('Max Value:', maxVal)
    # print('')
    
    print('Mean:', meanVal)
    print('Median:', medianVal)
    print('')
    
    print('Standard deviation:', stdDev)
    print('Standard Error:', stdErr)
    
    return None

# *** End Descriptive Statistic Block ***

# *** Sample Sizes < 25 Block ***

def getDegFree(scores):
    """Reads in a list of scores. Calculates the degrees of freedom for a 
    single sample test.
    Returns the degrees of freedom."""
    
    numScores = len(scores)
    doF = numScores - 1
    print('Degrees of Freedom:', doF)
    
    return doF

def getLogScore(scores):
    """Takes a list of scored values and calculates the log value for
    those scores. Returns a list containing log values of scores."""
    
    sortedScores = selSort(scores)
    logScores = np.log(sortedScores)
    
    logScoreList = []
                
    for elem in logScores:
        rndLogScores = round(elem, 2)
        logScoreList.append(rndLogScores)
        
    return logScoreList

def expLogScores(logConfIntVals):
    """Reads in the confidence interval list object calculated earlier.
    Calculates the exp value to return the log value to the standard value.
    Returns the standard confidence interval values."""
    
    logLCL = logConfIntVals[0]
    logUCL = logConfIntVals[1]
    logMean = logConfIntVals[2]
    
    LCL = round((np.exp(logLCL)), 3)
    UCL = round((np.exp(logUCL)), 3)
    meanVal = round((np.exp(logMean)), 3)
    
    confIntVals = [LCL, UCL]
    
    return confIntVals

def getConfIntvals(scores, alpha):
    """Calculates the confidence intervals for a list of data. 
    Takes in a list of data and the alpha level, calculates tabled-T values
    and confidence level values.
    Returns a list object with the mean and confidence interval values."""
    
    alphaTwoTail = (1 - (alpha / 2))    
    
    descriptStats = []
    descriptStats = getDescriptStats(scores)
    doF = descriptStats[0] - 1
       
    tCrit = round(scipy.stats.t.ppf(alphaTwoTail, doF), 4)
    stdErr = descriptStats[5]
    critVal = (stdErr * tCrit)
    
    meanVal = descriptStats[2]
    LCL = np.round((meanVal - critVal), 4)
    UCL = np.round((meanVal + critVal), 4)
    confIntVals = [LCL, UCL, meanVal]
    
    return confIntVals

def reportConfInt(confIntVals, alpha):
    """Reads the upper and lower confidence interval along with the mean.
    Prints the interval. 
    Return a list object with the LCL and UCL values."""
    
    LCL = int(confIntVals[0])
    UCL = int(confIntVals[1])
    #meanVal = int(confIntVals[2])
    
    print('\n** Confidence Interval **\n')
    
    #print('Confidence Interval:', LCL, '<', meanVal, '<', UCL, '\n')
    
    rptPct = int((1 - alpha) * 100)
    valType = 'median'
        
    print('The', rptPct, '% confidence interval around the', valType,'is between',
          LCL, 'seconds and', UCL, 'seconds.')
    
    rptConfIntVals = [LCL, UCL]
    
    return rptConfIntVals

# *** End Samples < 25 Block ***

# *** Sample Sizes > 25 Block ***

def selSort(scores):
    """Assumes that L is a list of elements that can be compared using >.
    Sorts L in ascending order.
    Returns a list object containing the sorted scores."""
    
    # Selection sort - order n^2
    sortedScores = []
    
    for i in range(len(scores) - 1):
        #Invariant: the list scores[:i] is sorted
        minIndx = i
        minVal = scores[i]
        j = i + 1
        while j < len(scores):
            if minVal > scores[j]:
                minIndx = j
                minVal = scores[j]
            j += 1
        temp = scores[i]
        scores[i] = scores[minIndx]
        scores[minIndx] = temp
        sortedScores = scores
        
    return sortedScores

def setSampleSizeFlag(scores):
    """Reads in a list of scores. Counts scores and sets a flag for the
    test type to use.
    Returns the large sample flag as lrgSampFlag."""
    
    sampSizeFlag = 0
    numScores = getNumScores(scores)
    
    if numScores >= 25:
        sampSizeFlag = 1
    
    return sampSizeFlag

def calcLrgSampIndx(scores, alpha, pctile):
    """Reads in a list of scores, then calculates the z-score formula
    for a specific percentile. This determines the index values for 
    identifying the values that comprise the interval. 
    Returns the index values."""
    
    numScores = getNumScores(scores)
    
    npVal = (numScores * pctile)
    qVal = (1 - pctile)
   
    zVal = scipy.stats.norm.ppf(1 - (alpha / 2))
    
    confIntAdj = zVal * ((npVal * qVal) ** 0.5)
    
    rawLCLIndex = (npVal - confIntAdj)
    rawUCLIndex = (npVal + confIntAdj)
            
    LCLIndex = math.ceil(rawLCLIndex)
    UCLIndex = math.ceil(rawUCLIndex)
            
    lrgSampIndx = [LCLIndex, UCLIndex]
        
    return lrgSampIndx

def getIndxCIVals(sortedScores, lrgSampIndx):
    """Reads in a list of sorted scores and the indexes for large samples.
    Uses the index to select the original interval values from the sorted 
    list. 
    Returns the values at the index."""
    
    LCLIndex = lrgSampIndx[0]
    UCLIndex = lrgSampIndx[1]
    
    LCLVal = sortedScores[LCLIndex - 1]
    UCLVal = sortedScores[UCLIndex - 1]
        
    indxConfIntVals = [LCLVal, UCLVal]
    
    return indxConfIntVals
            
def scoreTaskTimes(scores, alpha, pctile):
    """Reads in a list of scores and the alpha level. Runs the procedures
    necessary to calculate the geometric mean of the task time.
    Returns the confidence interval values for use in test harnesses."""
    
    sortedScores = selSort(scores)
    #print(sortedScores)
    
    descriptStats = getDescriptStats(sortedScores)
    reportDescriptStats(descriptStats)
    sampSizeFlag = setSampleSizeFlag(scores)
    
    if sampSizeFlag == 1:
        lrgSampIndx = calcLrgSampIndx(scores, alpha, pctile)
        indxConfIntVals = getIndxCIVals(sortedScores, lrgSampIndx)
        rptConfIntVals = reportConfInt(indxConfIntVals, alpha)
    
    else:
        logScores = getLogScore(sortedScores)
        print(logScores)
        
        logConfIntVals = getConfIntvals(logScores, alpha)
        print(logConfIntVals)
    
        confIntVals = expLogScores(logConfIntVals)
        rptConfIntVals = reportConfInt(confIntVals, alpha)
    
    
    
    return rptConfIntVals


### *******************************************************
### Use this block if you have data as a list of values. 
### *******************************************************
# If you have a list of values, enter them in these brackets. 

#scores = [40, 36, 53, 56, 110, 48, 34, 44, 30, 40, 80]
#alpha = 0.05
#pctile = 0.5 # Percentile. Median = 0.5
#scoreTaskTimes(scores, alpha, pctile)

### *******************************************************


### *******************************************************
### Use this block if you have data in a .csv file.
### *******************************************************
# If you have a data table to score, use this block.
# Place your .csv file in the same directory as this file. 
# File Name Example: 'FakeUXData.csv'

inFile = 'Fake_Task_Time_Data.csv'
df = pd.read_csv(inFile, index_col=False)

# Include the column name you wish to calculate the scores on.
# Example: 'TaskTimes_5'

scores = df['TT_1']
alpha = 0.05
pctile = .5
scoreTaskTimes(scores, alpha, pctile)
### *******************************************************














