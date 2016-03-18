from __future__ import division
import csv
import math
import argparse

def mpcTree(guessFileName, data, burnin, mpcTreeName, farms, burninInLines):

    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)

    currentLine = readData.next()

    while currentLine[0].startswith("#"):
        currentLine = readData.next()

    # no header line

    currentLine = readData.next()

    farms.insert(0, 'state')
    stateData = list()
    totalStates = 0
    totalCountedStates = 0
    while currentLine is not None:
        totalStates += 1
        if burninInLines:
            burningIn = totalStates < burnin
        else:
            burningIn = int(currentLine[0]) < burnin

        if not burningIn:
            totalCountedStates += 1
            referenceDict = dict()
            for i in range(len(farms)):
                referenceDict[farms[i]] = currentLine[i]
            stateData.append(referenceDict)

        try:
            currentLine = readData.next()
        except StopIteration:
            currentLine = None
    for network in stateData:
        totalLogParentCredibility = 0
        for i in range(1,len(farms)):
            currentFarm = farms[i]
            logParentCredibility = 0
            for guessDict in data:
                if guessDict["name"] == currentFarm:
                    logParentCredibility = math.log((guessDict[network[currentFarm]])/totalStates)
            totalLogParentCredibility = totalLogParentCredibility + logParentCredibility
        network['TLPC']=totalLogParentCredibility
    currentHighestTLPC = float('-inf')
    currentBestNetworkMPC = None
    for network in stateData:
        if network['TLPC']>currentHighestTLPC:
            print 'State '+network['state']+": new highest credibility, exp(" + str(network['TLPC']) + ")"
            currentHighestTLPC = network['TLPC']
            currentBestNetworkMPC = network
    print
    writeData = csv.writer(open(mpcTreeName, 'w'))
    mpcBestGuessLine = list()
    mpcBestGuessLine.append('MPC tree')
    for i in range(1,len(farms)):
        mpcBestGuessLine.append(currentBestNetworkMPC[farms[i]])
    individualMultCredLine = list()
    individualMultCredLine.append('MPC credibility')
    for i in range(1,len(farms)):
        for guessDict in data:
            if guessDict["name"]==farms[i]:
                individualMultCredLine.append(str(guessDict[currentBestNetworkMPC[farms[i]]]/totalStates))
    for case in range(len(farms)):
        line = [farms[case], mpcBestGuessLine[case], individualMultCredLine[case]]    
        writeData.writerow(line)
    print


def getHostList(guessFileName):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)

    currentLine = readData.next()

    while currentLine[0].startswith("#"):
        currentLine = readData.next()

    headers = currentLine
    headers.remove("state")
    for index, item in enumerate(headers):
        if item.endswith('_infector'):
            headers[index] = item[:-9]
    return headers


def getGuesses(guessFileName, headers, burnin, burninInLines):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)

    currentLine = readData.next()

    while currentLine[0].startswith("#"):
        currentLine = readData.next()

    data = list()

    for farmName in headers:
        currentDict = dict()
        currentDict['name']=farmName
        data.append(currentDict)
    for currentDict in data:
        for farmName in headers:
            currentDict[farmName]=0
        currentDict['Start']=0

    # no header line

    currentLine = readData.next()

    totalStates = 0

    while currentLine!=None:
        totalStates += 1

        if burninInLines:
            burningIn = totalStates < burnin
        else:
            burningIn = int(currentLine[0]) < burnin

        if not burningIn:
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
    parser = argparse.ArgumentParser("Processes a BEAST transmission tree output file to obtain the maximum parent"
                                     " credibility transmisstion tree")
    parser.add_argument("beastOutput", help="Name of a BEAST transmission tree output file")
    parser.add_argument("outputName", help="CSV file to which the MPC tree is to be written")
    parser.add_argument("-b", "--burnin", default=1000, help="Integer: burnin (default=1000 lines)")
    parser.add_argument("-t", "--burninInLines", default=True, help="Boolean: whether the burnin refers to the number "
                                                                    "of lines in the input file. If false, it refers "
                                                                    "to the number of MCMC states instead")

    args = parser.parse_args()

    print "Reading list of child hosts"
    childFarmList = getHostList(args.beastOutput)
    print "Reading BEAST output"
    guesses = getGuesses(args.beastOutput, childFarmList, args.burnin, args.burninInLines)
    print "Finding maximum parent credibility network (output to "+args.outputName+")"
    mpcTree(args.beastOutput, guesses, args.burnin, args.outputName, childFarmList, args.burninInLines)

if __name__ == '__main__':
    main()