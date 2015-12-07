from __future__ import division
import csv
import sys
import math

def mpcTree(guessFileName, data, burnin, mpcTreeName, farms):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)
    readData.next()
    readData.next()
    readData.next()
    farms.insert(0, 'state')
    stateData = list()
    currentLine = readData.next()
    totalStates = 0
    while currentLine!=None:
        if int(currentLine[0])>burnin:
            referenceDict = dict()
            for i in range(len(farms)):
                referenceDict[farms[i]] = currentLine[i]
            stateData.append(referenceDict)
            totalStates = totalStates+1
        try:
            currentLine = readData.next()
        except StopIteration:
            currentLine=None
    for network in stateData:
        totalParentCredibility = 0
        totalLogParentCredibility = 0
        for i in range(1,len(farms)):
            currentFarm = farms[i]
            parentCredibility = 0
            logParentCredibility = 0
            for guessDict in data:
                if guessDict["name"]==currentFarm:
                    parentCredibility = (guessDict[network[currentFarm]])/totalStates
                    if parentCredibility == 0:
                        print "zero credibility"
                    logParentCredibility = math.log((guessDict[network[currentFarm]])/totalStates)
            totalParentCredibility = totalParentCredibility + parentCredibility
            totalLogParentCredibility = totalLogParentCredibility + logParentCredibility
        network['TPC']=totalParentCredibility
        network['TLPC']=totalLogParentCredibility
    currentHighestTPC = 0
    currentHighestTLPC = float('-inf')
    currentBestNetworkMPC = None
    currentBestNetworkMSPC = None
    for network in stateData:
        if network['TPC']>currentHighestTPC:
            print 'State '+network['state']+": new highest additive credibility, " + str(network['TPC'])
            currentHighestTPC = network['TPC']
            currentBestNetworkMSPC = network
        if network['TLPC']>currentHighestTLPC:
            print 'State '+network['state']+": new highest multiplicative credibility, exp(" + str(network['TLPC']) + ")"
            currentHighestTLPC = network['TLPC']
            currentBestNetworkMPC = network
    print
    writeData = csv.writer(open(mpcTreeName, 'w'))
    mpcBestGuessLine = list()
    mpcBestGuessLine.append('MPC tree')
    mspcBestGuessLine = list()
    mspcBestGuessLine.append('MSPC tree')
    for i in range(1,len(farms)):
        mpcBestGuessLine.append(currentBestNetworkMPC[farms[i]])
        mspcBestGuessLine.append(currentBestNetworkMSPC[farms[i]]) 
    individualMultCredLine = list()
    individualMultCredLine.append('MPC credibility')
    individualAddCredLine = list()
    individualAddCredLine.append('MSPC credibility')
    for i in range(1,len(farms)):
        for guessDict in data:
            if guessDict["name"]==farms[i]:
                individualMultCredLine.append(str(guessDict[currentBestNetworkMPC[farms[i]]]/totalStates))
                individualAddCredLine.append(str(guessDict[currentBestNetworkMSPC[farms[i]]]/totalStates))
    for case in range(len(farms)):
        line = [farms[case], mpcBestGuessLine[case], individualMultCredLine[case]]    
        writeData.writerow(line)
    print
    
def getFarmList(guessFileName, includeStart):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)
    readData.next()
    readData.next()
    headers = readData.next()
    headers.remove("state")
    for index, item in enumerate(headers):
        if item.endswith('_infector'):
            headers[index] = item[:-9]
    if includeStart:
        headers.append("Start")
    return headers

def getGuesses(guessFileName, headers, burnin):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)
    data = list()
    readData.next()
    readData.next()
    readData.next()
    for farmName in headers:
        currentDict = dict()
        currentDict['name']=farmName
        data.append(currentDict)
    for currentDict in data:
        for farmName in headers:
            currentDict[farmName]=0
        currentDict['Start']=0
    currentLine = readData.next()
    while currentLine!=None:
        if int(currentLine[0])>burnin:           
            for i in range(1, len(currentLine)):
                for currentDict in data:
                    if currentDict['name']==headers[i-1]:
                        currentDict[currentLine[i]]=currentDict[currentLine[i]]+1
        try:
            currentLine=readData.next()
        except StopIteration:
            currentLine=None
    return data

def main():
    print "Reading list of child farms"
    childFarmList = getFarmList(sys.argv[1], False)
    print "Reading BEAST output"
    guesses = getGuesses(sys.argv[1], childFarmList, 1000000)
    print "Finding maximum parent credibility network (output to "+sys.argv[2]+")"
    mpcTree(sys.argv[1], guesses, 1000000, sys.argv[2], childFarmList)

if __name__ == '__main__':
    main()