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
from scipy.stats import sem
import scipy.stats

def getDescriptStats(scores):
    """Reads in a list of scores. Calculates the descriptive statistics
    for the list of values.
    Returns a list object containing the descriptive statistics-related values."""
    
    # Calculate (n) based on length of data list.
    print('')
    numScores = len(scores)
    print('n:', numScores)
    
    # Sum of scores in data list.
    sumScores = round((sum(scores)), 2)
    print('Sum:', sumScores)
    print('')
    
    # Range of values.
    minVal = round((min(scores)), 2)
    maxVal = round((max(scores)), 2)
    rangeVals = round((maxVal - minVal), 2)
    print('Range:', rangeVals)
    print('Min Value:', minVal)
    print('Max Value:', maxVal)
    print('')
    
    # Mean of scores in data list.
    meanScores = round(stats.mean(scores), 2)
    print('Mean:', meanScores)
    
    # Standard deviation of scores in data list.
    stdDev = round(stats.stdev(scores), 2)
    print('Standard deviation:', stdDev)
    
    # Standard error of scores in data list.
    stdErr = round(sem(scores), 2)
    print('Standard Error:', stdErr)
    
    return numScores, sumScores, meanScores, stdDev, stdErr

def getDegFree(scores):
    """Reads in a list of scores. Calculates the degrees of freedom for a 
    single sample test.
    Returns the degrees of freedom."""
    
    numScores = len(scores)
    doF = numScores - 1
    print('Degrees of Freedom:', doF)
    
    return doF

def getConfIntvals(scores, alpha):
    """Calculates the confidence intervals for a list of data. 
    Takes in a list of data and the alpha level, calculates tabled-T values
    and confidence level values.
    Returns a list object with the mean and confidence interval values."""
    
    # Calculates the value for use in the tCrit calculation.
    alphaTwoTail = (1 - (alpha / 2))    
    
    # Creates an empty list object for getDescriptStats return values
    descriptStats = []
    
    # This call returns a list of values calculated in getDescriptStats
    descriptStats = getDescriptStats(scores)
    
    # Calculates degrees of freedom based on index value of (n), [0] from getDescriptStats
    doF = descriptStats[0] - 1
       
    # Calculates the tabled t value, or the t-critical value for a two-tailed t-test.
    tCrit = round(scipy.stats.t.ppf(alphaTwoTail, doF), 2)
    
    # Calculates the critical values based on StdError (index 4) from getDescriptStats. 
    stdErr = descriptStats[4]
    critVal = (stdErr * tCrit)
    
    # Calculates the values for the confidence interval.
    meanVal = descriptStats[2]
    LCL = round((meanVal - critVal), 2)
    UCL = round((meanVal + critVal), 2)
    
    confIntVals = [LCL, UCL, meanVal]
    
    return confIntVals

def reportConfInt(confIntVals):
    """Reads the upper and lower confidence interval along with the mean.
    Prints the interval. 
    Return None."""
    
    LCL = confIntVals[0]
    UCL = confIntVals[1]
    meanVal = confIntVals[2]
            
    print('Confidence Interval:', LCL, '<', meanVal, '<', UCL, '\n')
    print('The population value of the mean at the alpha', alpha, 'level is between', LCL,
          'seconds and', UCL, 'seconds.')
    
    return None

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

def getLogScore(sortedScores):
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
    
    LCL = round((np.exp(logLCL)), 2)
    UCL = round((np.exp(logUCL)), 2)
    meanVal = round((np.exp(logMean)), 2)
    
    confIntVals = [LCL, UCL, meanVal]
    
    return confIntVals
    
def scoreTaskTimes(scores, alpha):
    """Reads in a list of scores and the alpha level. Runs the procedures
    necessary to calculate the geometric mean of the task time.
    Returns None."""
    
    sortedScores = selSort(scores)
    logScores = getLogScore(sortedScores)
    logConfIntVals = getConfIntvals(logScores, alpha)
    confIntVals = expLogScores(logConfIntVals)
    reportConfInt(confIntVals)
    
    return None


### *******************************************************
### Use this block if you have data as a list of values. 
### *******************************************************
# If you have a list of values, enter them in these brackets. 
#scores = [62, 56, 71, 83, 52, 114, 127, 191, 78, 94, 67, 99, 102, 153]

#alpha = 0.05
#scoreTaskTimes(scores, alpha)
### *******************************************************


### *******************************************************
### Use this block if you have data in a .csv file.
### *******************************************************
# If you have a data table to score, use this block.
# Place your .csv file in the same directory as this file. 
# File Name Example: 'FakeUXData.csv'

inFile = 'FakeTaskTimeData.csv'
df = pd.read_csv(inFile, index_col=False)

# Include the column name you wish to calculate the scores on.
# Example: 'TaskTimes_5'

scores = df['TT_1']
alpha = 0.05
scoreTaskTimes(scores, alpha)
### *******************************************************














