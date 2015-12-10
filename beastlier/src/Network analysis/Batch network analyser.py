from __future__ import division
import csv
import sys
import math

def getCorrectParents(correctFileName):
    readData = csv.reader(open(correctFileName, 'rU'))
    correctAnswers=dict()
    for row in readData:
        if row[1]=="start":
            row[1]="Start"
        if row[0]!="Child":
            correctAnswers[row[0]]=row[1]
    return correctAnswers

def getCorrectDS(farmList, correctParents):
    correctDescendantSet=dict()
    for farm in farmList:
        correctDescendantSet[farm]=findDescendantSet(farm, farmList, correctParents)
    return correctDescendantSet

def getCorrectSM(farmList, correctParents):
    correctSubtrees=set()
    for farm in farmList:
        correctSubtrees.add(findSubtreeMembers(farm, farmList, correctParents))
    return correctSubtrees
    
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

def getDSGuesses(guessFileName, headers, burnin):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)
    data = list()
    readData.next()
    readData.next()
    readData.next()
    for farmName in headers:
        currentDict = dict()
        currentDict['name']=farmName
        data.append(currentDict)
    currentLine = readData.next()
    while currentLine!=None:
        if int(currentLine[0])>burnin:
            tempDict = dict()
            for i in range(0, len(headers)):
                tempDict[headers[i]]=currentLine[i+1]
            for i in range(0, len(headers)):
                currentDS = findDescendantSet(headers[i],headers,tempDict)
                for currentDict in data:
                    if currentDict['name']==headers[i]:
                        if currentDS in currentDict:
                            currentDict[currentDS]=currentDict[currentDS]+1
                        else:
                            currentDict[currentDS]=1
        try:
            currentLine=readData.next()
        except StopIteration:
            currentLine=None
    return data

def getSMGuesses(guessFileName, headers, burnin):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)
    data = dict()
    readData.next()
    readData.next()
    readData.next()
    currentLine = readData.next()
    while currentLine!=None:
        if int(currentLine[0])>burnin:
            tempDict = dict()
            for i in range(0, len(headers)):
                tempDict[headers[i]]=currentLine[i+1]
            for i in range(0, len(headers)):
                currentSM = findSubtreeMembers(headers[i],headers,tempDict)
                if currentSM in data:
                    data[currentSM] = data[currentSM]+1
                else:
                    data[currentSM] = 1
        try:
            currentLine=readData.next()
        except StopIteration:
            currentLine=None
    return data 

def findDescendantSet(farmString, allFarmsList, parentDict):
    descendants = set()
    for farm in allFarmsList:
        if farm!=farmString and farm!='Start' and farm!='state':
            currentFarm = farm
            startReached = False
            farmReached = False
            while (not startReached and not farmReached):
                currentParent = parentDict[currentFarm]
                if currentParent == farmString:
                    farmReached = True
                if currentParent == 'Start':
                    startReached = True
                currentFarm = currentParent
            if farmReached:
                descendants.add(farm)
    return frozenset(descendants)
    
def findSubtreeMembers(farmString, allFarmsList, parentDict):
    members = set()
    members.add(farmString)
    for farm in allFarmsList:
        if farm!=farmString and farm!='Start' and farm!='state':
            currentFarm = farm
            startReached = False
            farmReached = False
            while (not startReached and not farmReached):
                currentParent = parentDict[currentFarm]
                if currentParent == farmString:
                    farmReached = True
                if currentParent == 'Start':
                    startReached = True
                currentFarm = currentParent
            if farmReached:
                members.add(farm)
    return frozenset(members)

def analyse(headers, correctAnswers, guesses, outputFileName, briefOutputFileName, rowstring):
    outputFile = open(outputFileName, 'w')
    briefOutputFile = open(briefOutputFileName, 'a')
    briefWriter = csv.writer(briefOutputFile)
    totalCorrectGuesses=0
    totalFarms=0
    for farm in headers:
        farmDict = dict()
        total=0
        correctGuess=""
        if farm!="Start":
            totalFarms=totalFarms+1
            outputFile.write("Farm "+farm+" (correct: "+correctAnswers[farm]+")\n")
            for currentDict in guesses:
                for farm2 in headers:
                    if currentDict['name']==farm:
                        total=total+currentDict[farm2]
                for farm2 in headers:
                    if currentDict['name']==farm:
                        farmDict[farm2] = currentDict[farm2]
                        if correctAnswers[farm]==farm2:
                            outputFile.write("**"+farm2 +"** : "+str(currentDict[farm2]) + " (" + str(round((currentDict[farm2]/total)*100,1)) + "%)\n")
                            correctGuess=farm2
                        else:
                            if currentDict[farm2]!=0:
                                outputFile.write(farm2 +": "+str(currentDict[farm2]) + " (" + str(round((currentDict[farm2]/total)*100,1)) + "%)\n")
            outputFile.write("Total: " + str(total) + '\n')
            guessedRight=True
            currentMaxNumber = farmDict[correctGuess]
            currentMaxName = correctGuess
            for farm2 in headers:
                if farmDict[farm2]>currentMaxNumber:
                    currentMaxNumber = farmDict[farm2]
                    currentMaxName = farm2
                    guessedRight=False                
            if guessedRight:
                totalCorrectGuesses = totalCorrectGuesses+1
            if guessedRight:
                outputFile.write("Guessed correctly: "+currentMaxName +" ("+str(round((currentMaxNumber/total)*100,1))+"%)\n")
                briefWriter.writerow([rowstring, str(currentMaxNumber/total), 'TRUE'])
            else:
                outputFile.write("Guessed incorrectly: "+currentMaxName +" ("+str(round((currentMaxNumber/total)*100,1))+"%)\n")
                briefWriter.writerow([rowstring, str(currentMaxNumber/total), 'FALSE'])
                outputFile.write("Correct guesses: "+correctGuess + " ("+str(round((farmDict[correctGuess]/total)*100,1))+"%)\n")                
            outputFile.write('\n')
    outputFile.write("Correctly guessed farms: "+str(totalCorrectGuesses)+" ("+str(round((totalCorrectGuesses/totalFarms)*100,1))+"%)\n")
    outputFile.close()
    
def analyseDS(childFarms, correctAnswers, guesses, outputFileName):
    outputFile = open(outputFileName, 'w')
    totalCorrectGuesses=0
    totalFarms=0
    for farm in childFarms:
        farmDict = dict()
        total=0
        correctGuess=None
        totalFarms=totalFarms+1
        outputFile.write("Farm "+farm+" (correct: "+str(correctAnswers[farm])+")\n")
        for currentDict in guesses:
            if currentDict['name']==farm:
                tempDescendantSets = currentDict.keys()
                for key in tempDescendantSets:
                    if key!='name':
                        total=total+currentDict[key]
                for key in tempDescendantSets:
                    if key!='name':
                        farmDict[key] = currentDict[key]
                        if(correctAnswers[farm]==key):
                            outputFile.write("**"+str(key)+"** : "+str(currentDict[key]) + " (" + str(round((currentDict[key]/total)*100,1)) + "%)\n")
                            correctGuess=key
                        else:
                            outputFile.write(str(key) +": "+str(currentDict[key]) + " (" + str(round((currentDict[key]/total)*100,1)) + "%)\n")
        if not farmDict.has_key(correctAnswers[farm]):
            correctGuess=correctAnswers[farm]
            farmDict[correctAnswers[farm]]=0
        outputFile.write("Total: " + str(total) + '\n')
        guessedRight=True
        currentMaxNumber = farmDict[correctGuess]
        currentMaxName = correctGuess
        for key in farmDict.keys():
            if key!='name':
                if farmDict[key]>currentMaxNumber:
                    currentMaxNumber = farmDict[key]
                    currentMaxName = key
                    guessedRight=False                
        if guessedRight:
            totalCorrectGuesses = totalCorrectGuesses+1
        if guessedRight:
            outputFile.write("Guessed correctly: "+str(currentMaxName)+" ("+str(round((currentMaxNumber/total)*100,1))+"%)\n")
        else:
            outputFile.write("Guessed incorrectly: "+str(currentMaxName)+" ("+str(round((currentMaxNumber/total)*100,1))+"%)\n")
            outputFile.write("Correct guesses: "+str(correctGuess)+" ("+str(round((farmDict[correctGuess]/total)*100,1))+"%)\n")                
        outputFile.write('\n')
    outputFile.write("Correctly guessed descendant sets: "+str(totalCorrectGuesses)+" ("+str(round((totalCorrectGuesses/totalFarms)*100,1))+"%)\n")
    outputFile.close()
    
def mdscTree(guessFileName, data, burnin, mdscTreeName, correctAnswers):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)
    readData.next()
    readData.next()
    farms = readData.next()
    parentStateData = list()
    dSStateData = list()
    currentLine = readData.next()
    totalStates = 0
    while currentLine!=None:
        if int(currentLine[0])>burnin:
            parentReferenceDict = dict()
            for i in range(len(farms)):
                parentReferenceDict[farms[i]] = currentLine[i]
            parentStateData.append(parentReferenceDict)
            dSReferenceDict = dict()
            dSReferenceDict['state'] = currentLine[0]
            for i in range(1,len(farms)):
                dSReferenceDict[farms[i]] = findDescendantSet(farms[i], farms, parentReferenceDict)
            dSStateData.append(dSReferenceDict)
            totalStates = totalStates+1          
        try:
            currentLine = readData.next()
        except StopIteration:
            currentLine=None
    for i in range(len(dSStateData)):
        currentDSData = dSStateData[i]
        currentParentData = parentStateData[i]
        totalDSCredibility = 0
        totalLogDSCredibility = 0
        for j in range(1,len(farms)):
            currentFarm = farms[j]
            dSCredibility = 0
            logDSCredibility = 0
            for guessDict in data:
                if guessDict["name"]==currentFarm:
                    dSCredibility = (guessDict[currentDSData[currentFarm]])/totalStates
                    logDSCredibility = math.log((guessDict[currentDSData[currentFarm]])/totalStates)
            totalDSCredibility = totalDSCredibility + dSCredibility
            totalLogDSCredibility = totalLogDSCredibility + logDSCredibility
            currentParentData['TDSC']=totalDSCredibility
            currentParentData['TLDSC']=totalLogDSCredibility
    currentHighestTDSC = 0
    currentHighestTLDSC = float('-inf')
    currentBestNetworkMDSC = None
    currentBestNetworkMSDSC = None
    for network in parentStateData:
        if network['TDSC']>currentHighestTDSC:
            print 'State '+network['state']+": new highest additive credibility, " + str(network['TDSC'])
            currentHighestTDSC = network['TDSC']
            currentBestNetworkMSDSC = network
        if network['TLDSC']>currentHighestTLDSC:
            print 'State '+network['state']+": new highest multiplicative credibility, exp(" + str(network['TLDSC']) + ")"
            currentHighestTLDSC = network['TLDSC']
            currentBestNetworkMDSC = network
    print
    writeData = csv.writer(open(mdscTreeName, 'w'))
    writeData.writerow(farms)
    correctLine = list()
    correctLine.append('Actual parent')
    mdscBestGuessLine = list()
    mdscBestGuessLine.append('MDSC network')
    msdscBestGuessLine = list()
    msdscBestGuessLine.append('MSDSC network')
    mdscCorrectParents = 0
    msdscCorrectParents = 0
    for i in range(1,len(farms)):
        correctLine.append(correctAnswers[farms[i]])
        mdscBestGuessLine.append(currentBestNetworkMDSC[farms[i]])
        msdscBestGuessLine.append(currentBestNetworkMSDSC[farms[i]])
        if correctAnswers[farms[i]]==currentBestNetworkMDSC[farms[i]]:
            mdscCorrectParents = mdscCorrectParents+1
        if correctAnswers[farms[i]]==currentBestNetworkMSDSC[farms[i]]:
            msdscCorrectParents = msdscCorrectParents+1
    writeData.writerow(correctLine)       
    writeData.writerow(mdscBestGuessLine)
    writeData.writerow(msdscBestGuessLine)
    print "MDSC: Parent is correct in "+str(mdscCorrectParents)+" out of "+str(len(farms)-1)+" cases (" + str(round((mdscCorrectParents/(len(farms)-1))*100,1)) + "%)"
    print "MSDSC: Parent is correct in "+str(msdscCorrectParents)+" out of "+str(len(farms)-1)+" cases (" + str(round((msdscCorrectParents/(len(farms)-1))*100,1)) + "%)"
    print
    
def mpcTree(guessFileName, data, burnin, mpcTreeName, correctAnswers, farms):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)
    readData.next()
    readData.next()
    readData.next()
    stateData = list()
    farms.insert(0, 'state')
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
        if network['TLPC']>currentHighestTLPC:
            print 'State '+network['state']+": new highest multiplicative credibility, exp(" + str(network['TLPC']) + ")"
            currentHighestTLPC = network['TLPC']
            currentBestNetworkMPC = network
    print
    writeData = csv.writer(open(mpcTreeName, 'w'))
    correctLine = list()
    correctLine.append('Actual parent')
    mpcBestGuessLine = list()
    mpcBestGuessLine.append('MPC tree')
    mpcCorrectParents = 0
    mspcCorrectParents = 0
    for i in range(1,len(farms)):
        correctLine.append(correctAnswers[farms[i]])
        mpcBestGuessLine.append(currentBestNetworkMPC[farms[i]])
        if correctAnswers[farms[i]]==currentBestNetworkMPC[farms[i]]:
            mpcCorrectParents = mpcCorrectParents+1

    individualMultCredLine = list()
    individualMultCredLine.append('MPC credibility')
    for i in range(1,len(farms)):
        for guessDict in data:
            if guessDict["name"]==farms[i]:
                individualMultCredLine.append(str(guessDict[currentBestNetworkMPC[farms[i]]]/totalStates))
    for i in range(0, len(farms)):
        row = list()
        row.append(farms[i])
        row.append(correctLine[i])
        row.append(mpcBestGuessLine[i])
        row.append(individualMultCredLine[i])
        writeData.writerow(row)
    print "MPC: Parent is correct in "+str(mpcCorrectParents)+" out of "+str(len(farms)-1)+" cases (" + str(round((mpcCorrectParents/(len(farms)-1))*100,1)) + "%)"
    print
    return mpcCorrectParents/(len(farms)-1)
    
def msmcTree(guessFileName, data, burnin, msmcTreeName, correctAnswers):
    readData = csv.reader(open(guessFileName), dialect=csv.excel_tab)
    readData.next()
    readData.next()
    farms = readData.next()
    parentStateData = list()
    subtreeStateData = list()
    currentLine = readData.next()
    totalStates = 0
    while currentLine!=None:
        if int(currentLine[0])>burnin:
            referenceDict = dict()
            for i in range(len(farms)):
                referenceDict[farms[i]] = currentLine[i]
            parentStateData.append(referenceDict)
            referenceSet = set()
            for i in range(1, len(farms)):
                referenceSet.add(findSubtreeMembers(farms[i], farms, referenceDict))   
            subtreeStateData.append(referenceSet)
            totalStates = totalStates+1
        try:
            currentLine = readData.next()
        except StopIteration:
            currentLine=None
    for i in range(len(subtreeStateData)):
        totalSubtreeCredibility = 0
        totalLogSubtreeCredibility = 0
        for subtree in subtreeStateData[i]:
            totalSubtreeCredibility = totalSubtreeCredibility + data[subtree]/totalStates
            totalLogSubtreeCredibility = totalLogSubtreeCredibility + math.log(data[subtree]/totalStates)
        parentStateData[i]['TSMC']=totalSubtreeCredibility
        parentStateData[i]['TLSMC']=totalLogSubtreeCredibility
    currentHighestTSMC = 0
    currentHighestTLSMC = float('-inf')
    currentBestNetworkMSMC = None
    currentBestNetworkMSSMC = None
    for network in parentStateData:
        if network['TSMC']>currentHighestTSMC:
            print 'State '+network['state']+": new highest additive credibility, " + str(network['TSMC'])
            currentHighestTSMC = network['TSMC']
            currentBestNetworkMSSMC = network
        if network['TLSMC']>currentHighestTLSMC:
            print 'State '+network['state']+": new highest multiplicative credibility, exp(" + str(network['TLSMC']) + ")"
            currentHighestTLSMC = network['TLSMC']
            currentBestNetworkMSMC = network
    print
    writeData = csv.writer(open(msmcTreeName, 'w'))
    writeData.writerow(farms)
    correctLine = list()
    correctLine.append('Actual parent')
    msmcBestGuessLine = list()
    msmcBestGuessLine.append('MSMC tree')
    mssmcBestGuessLine = list()
    mssmcBestGuessLine.append('MSSMC tree')
    msmcCorrectParents = 0
    mssmcCorrectParents = 0
    for i in range(1,len(farms)):
        correctLine.append(correctAnswers[farms[i]])
        msmcBestGuessLine.append(currentBestNetworkMSMC[farms[i]])
        mssmcBestGuessLine.append(currentBestNetworkMSSMC[farms[i]])
        if correctAnswers[farms[i]]==currentBestNetworkMSMC[farms[i]]:
            msmcCorrectParents = msmcCorrectParents+1
        if correctAnswers[farms[i]]==currentBestNetworkMSSMC[farms[i]]:
            mssmcCorrectParents = mssmcCorrectParents+1
    writeData.writerow(correctLine)       
    writeData.writerow(msmcBestGuessLine)
    writeData.writerow(mssmcBestGuessLine)
    print "MSMC: Parent is correct in "+str(msmcCorrectParents)+" out of "+str(len(farms)-1)+" cases (" + str(round((msmcCorrectParents/(len(farms)-1))*100,1)) + "%)"
    print "MSSMC: Parent is correct in "+str(mssmcCorrectParents)+" out of "+str(len(farms)-1)+" cases (" + str(round((mssmcCorrectParents/(len(farms)-1))*100,1)) + "%)"
    print    
    

    
    
    

def main():
    
    with open('mpcOutputGamma.csv', 'wb') as outputFile:
        writer = csv.writer(outputFile)
        writer.writerow(['Simulation','Proportion_correct'])
      
        for i in range(1,11):
            for j in range(1,6):
                for k in range(1,4):
                    print 'Doing simulation '+str(i)+":"+str(j)+":"+str(k)
                    
                    mcmcOutput = '3LT/pTree_G_'+str(i)+'_'+str(j)+'_'+str(k)+'_fixedTT_3LT_actual.net.txt'
                    
                    trueNetwork = 'network_G_'+str(i)+'_'+str(j)+'_'+str(k)+'.csv'
                    detail = '3LT/detail_'+str(i)+'_'+str(j)+'_'+str(k)+'.txt'
                    briefdetail = '3LT/detail_brief.txt'
                    mpc = '3LT/mpc_'+str(i)+'_'+str(j)+'_'+str(k)+'.csv'
        
                    burnin = 1000000
                              
                    print "Reading list of child farms"
                    childFarmList = getFarmList(mcmcOutput, False)
                    print "Reading list of parent farms"
                    parentFarmList = getFarmList(mcmcOutput, True)
                    print "Reading correct parents from tree simulation"
                    correctParents = getCorrectParents(trueNetwork)
                    print "Reading BEAST output"
                    guesses = getGuesses(mcmcOutput, childFarmList, burnin)
                    print "Analysing parents (output to "+detail+")"
                    analyse(parentFarmList, correctParents, guesses, detail, briefdetail, str(i)+"_"+str(j)+"_"+str(k))
                    print "Finding maximum parent credibility network (output to mpcTree.csv)"
                    percent = mpcTree(mcmcOutput, guesses, burnin, mpc, correctParents, childFarmList)
        
                    
                    writer.writerow([str(i), str(j), str(k), str(percent)])
                    
                    print


if __name__ == '__main__':
    main()