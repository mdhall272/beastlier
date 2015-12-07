# Turn a CSV file of farm information into XML blocks for BEAST. Version 2.0,

import csv
import xml.dom.minidom as mini
from datetime import datetime
import time
import math
import argparse

# Main method

def generateXML(csvreader, dayOne, tnMean, tnSD, standardDevInDays,
                riemannSteps, outputFileName, chainLength, logEvery, tree,
                networkLogFileName, parameterLogFileName,
                startingNetworkFileName, estimateIncubationPeriods, 
                incubationDistribution):
    xmlBlock = mini.Document()
    beastElement = xmlBlock.createElement('beast')   
    taxaComment = xmlBlock.createComment('Taxa block.')
    taxaElement = xmlBlock.createElement('taxa')
    taxaElement.setAttribute('id','taxa')
    likelihoodComment = xmlBlock.createComment('CaseToCaseLikelihood block.')
    rootElement = xmlBlock.createElement('caseToCaseTransmissionLikelihood')
    rootElement.appendChild(createNestedReferenceBlock(xmlBlock,'virusTree','newick','startingTree'))      
    farmsetname = 'lesionDatedFarmCaseSet'
    farmsetid = 'lesionDatedFarmCases'
    farmname = 'lesionDatedFarmCase'
    if incubationDistribution=='T':
        distString = 'truncatedNormalDistributionModel'
    elif incubationDistribution=='L':
        distString = 'logNormalDistributionModel'
    elif incubationDistribution=='G':
        distString = 'gammaDistributionModel'
    farmSetElement = xmlBlock.createElement(farmsetname)
    rootElement.setAttribute('id', farmsetid)
    # Information common to all farms
    if estimateIncubationPeriods:
        if incubationDistribution=="T":
            farmSetElement.appendChild(createLowerTruncNormalBlock(xmlBlock, 'incubationPeriodDistribution', 'overallIncubationPeriodDistribution', '1', '1', '0', 'inc_TN_mean', 'inc_TN_stdev'))
        elif incubationDistribution=="G":
            farmSetElement.appendChild(createGammaBlock(xmlBlock, 'incubationPeriodDistribution', 'overallIncubationPeriodDistribution', '1', '1', 'inc_gamma_scale', 'inc_gamma_shape'))
        elif incubationDistribution=="L":
            farmSetElement.appendChild(createLogNormalBlock(xmlBlock, 'incubationPeriodDistribution', 'overallIncubationPeriodDistribution', '5', '2', 'inc_LN_mean', 'inc_LN_stdev'))
    else:
        farmSetElement.appendChild(createLowerTruncNormalBlock(xmlBlock, 'incubationPeriodDistribution', 'overallIncubationPeriodDistribution', tnMean, tnSD, "0", "inc_TN_mean", "inc_TN_stdev"))
    farmSetElement.appendChild(createParameterBlockNoID(xmlBlock, 'riemannSampleSize', int(riemannSteps), True))
    if startingNetworkFileName!=None:
        rootElement.appendChild(createStringBlock(xmlBlock, 'startingNetwork', startingNetworkFileName))      
    # Individual farm data processing starts here; read from a CSV file. The file must have the column headings "Farm_ID", "Exam_date", "Cull_date", "Oldest_lesion" and "Taxon".
    headerRow=csvreader.next()
    farmIDColumn=None
    examDateColumn=None
    cullDateColumn=None
    oldestLesionColumn=None
    taxonColumn=None
    # Find the correct columns for the required information
    for i in range(0,len(headerRow)):
        if headerRow[i]=='Farm_ID':
            farmIDColumn=i
        elif headerRow[i]=='Exam_date':
            examDateColumn=i
        elif headerRow[i]=='Cull_date':
            cullDateColumn=i
        elif headerRow[i]=='Oldest_lesion':
            oldestLesionColumn=i
        elif headerRow[i]=='Taxon':
            taxonColumn=i
    if(farmIDColumn==None or examDateColumn==None or cullDateColumn==None or oldestLesionColumn==None or taxonColumn==None):
        raise Exception('Not all required columns are present')
    # Go through the rows of the CSV file and make a farm XML and a taxon block for each
    currentRow=csvreader.next()
    while currentRow!=None:
        processedExamDate = datetime.strptime(currentRow[examDateColumn], '%d/%m/%Y')
        processedCullDate = datetime.strptime(currentRow[cullDateColumn], '%d/%m/%Y')
        intExamDate = (processedExamDate - dayOne).days
        intCullDate = (processedCullDate - dayOne).days
        taxa=currentRow[taxonColumn].split('$')
        taxonElement = xmlBlock.createElement('taxon')
        taxonElement.setAttribute('id', currentRow[taxonColumn])
        taxaElement.appendChild(taxonElement)
        samplingDateElement = xmlBlock.createElement('date')
        # Events (sample and cull) now happen at the _end_ of the day. Easier that way. Hopefully. Still need to add 1, since time 0 is the beginning
        # of day 0, and thus if you subtract day 0's time value from another you get the number of days elapsed between the ends.
        samplingDateElement.setAttribute('value', str(intExamDate+1))
        samplingDateElement.setAttribute('direction','forwards')
        samplingDateElement.setAttribute('units','days')
        samplingDateElement.setAttribute('origin', dayOne.strftime('%d/%m/%Y'))
        taxonElement.appendChild(samplingDateElement)
        # Want to be adding +1 to dates here to get the end of the required days.
        farmSetElement.appendChild(createLesionDatedFarmCaseElement(xmlBlock, currentRow[farmIDColumn], intExamDate+1, intCullDate+1, float(currentRow[oldestLesionColumn]), dayOne, standardDevInDays, taxa, distString))
        try:
            currentRow=csvreader.next()
        except StopIteration:
            currentRow=None
    beastElement.appendChild(taxaComment)
    beastElement.appendChild(taxaElement)
    
    newickComment = xmlBlock.createComment('The fixed tree in Newick format')
    beastElement.appendChild(newickComment)
    newickElement = xmlBlock.createElement('newick')
    newickElement.setAttribute('id', 'startingTree')
    newickElement.setAttribute('usingDates', 'true')
    treePlaceholder = xmlBlock.createTextNode(tree)
    newickElement.appendChild(treePlaceholder)
    beastElement.appendChild(newickElement)
    
    beastElement.appendChild(likelihoodComment)
    rootElement.appendChild(farmSetElement)
    beastElement.appendChild(rootElement)
    
    
    operatorComment = xmlBlock.createComment('Operators; at the moment there is only one')
    operatorsBlock = xmlBlock.createElement('operators')
    operatorsBlock.setAttribute('id','operators')
    branchLSOElement = xmlBlock.createElement('nodePaintingSwitchOperator')
    branchLSOElement.setAttribute('weight', str(1))
    operatorsBlock.appendChild(branchLSOElement)
    ftlElement = xmlBlock.createElement('caseToCaseTransmissionLikelihood')
    ftlElement.setAttribute('idref',farmsetid)
    branchLSOElement.appendChild(ftlElement)
    beastElement.appendChild(operatorComment)
    beastElement.appendChild(operatorsBlock)
    if estimateIncubationPeriods:
        if incubationDistribution=="T" or incubationDistribution=="L":
            meanOperator = xmlBlock.createElement('scaleOperator')
            stdevOperator = xmlBlock.createElement('scaleOperator')
            meanOperator.setAttribute('scaleFactor', "0.75")
            meanOperator.setAttribute('weight', "1")
            stdevOperator.setAttribute('scaleFactor', "0.75")
            stdevOperator.setAttribute('weight', "1")
            stdevOperator.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_'+incubationDistribution+'N_mean'))
            meanOperator.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_'+incubationDistribution+'N_stdev'))
            operatorsBlock.appendChild(meanOperator)
            operatorsBlock.appendChild(stdevOperator)
        elif incubationDistribution=="G":
            scaleOperator = xmlBlock.createElement('scaleOperator')
            shapeOperator = xmlBlock.createElement('scaleOperator')
            scaleOperator.setAttribute('scaleFactor', "0.75")
            scaleOperator.setAttribute('weight', "1")
            shapeOperator.setAttribute('scaleFactor', "0.75")
            shapeOperator.setAttribute('weight', "1")
            scaleOperator.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_scale'))
            shapeOperator.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_shape'))
            operatorsBlock.appendChild(scaleOperator)
            operatorsBlock.appendChild(shapeOperator)
            
    
    mcmcComment = xmlBlock.createComment('MCMC element')
    beastElement.appendChild(mcmcComment)
    mcmcBlock = xmlBlock.createElement('mcmc')
    mcmcBlock.setAttribute('id', 'mcmc')
    mcmcBlock.setAttribute('chainLength', str(chainLength))
    mcmcBlock.setAttribute('autoOptimize', 'false')
    beastElement.appendChild(mcmcBlock)
    posteriorBlock = xmlBlock.createElement('posterior')
    posteriorBlock.setAttribute('id', 'posterior')
    mcmcBlock.appendChild(posteriorBlock)
    priorBlock = xmlBlock.createElement('prior')
    if estimateIncubationPeriods:
        if incubationDistribution=="T" or incubationDistribution=="L":
            meanPrior = xmlBlock.createElement('uniformPrior')
            meanPrior.setAttribute('lower', '1')
            meanPrior.setAttribute('upper', '10')
            meanPrior.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_'+incubationDistribution+'N_mean'))
            priorBlock.appendChild(meanPrior)
            stdevPrior = xmlBlock.createElement('uniformPrior')
            stdevPrior.setAttribute('lower', '0')
            stdevPrior.setAttribute('upper', '5')
            stdevPrior.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_'+incubationDistribution+'N_stdev'))
            priorBlock.appendChild(stdevPrior) 
        elif incubationDistribution=="G":
            scalePrior = xmlBlock.createElement('uniformPrior')
            scalePrior.setAttribute('lower', '0')
            scalePrior.setAttribute('upper', '25')
            scalePrior.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_scale'))
            priorBlock.appendChild(scalePrior)
            shapePrior = xmlBlock.createElement('uniformPrior')
            shapePrior.setAttribute('lower', '0.04')
            shapePrior.setAttribute('upper', '1E100')
            shapePrior.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_shape'))
            priorBlock.appendChild(shapePrior) 
                                  
    priorBlock.setAttribute('id', 'prior')
    posteriorBlock.appendChild(priorBlock)
    likelihoodBlock = xmlBlock.createElement('likelihood')
    likelihoodBlock.setAttribute('id', 'likelihood')
    posteriorBlock.appendChild(likelihoodBlock)
    likelihoodBlock.appendChild(createReferenceBlock(xmlBlock, 'caseToCaseTransmissionLikelihood', farmsetid))
    mcmcBlock.appendChild(createReferenceBlock(xmlBlock, 'operators', 'operators'))
    loggerComment = xmlBlock.createComment('Log elements; the network log is a separate file')
    mcmcBlock.appendChild(loggerComment)
    networkLogBlock = xmlBlock.createElement('log')
    networkLogBlock.setAttribute('id', 'networkLog')
    networkLogBlock.setAttribute('logEvery', str(logEvery))
    networkLogBlock.setAttribute('fileName', networkLogFileName)
    networkLogBlock.setAttribute('overwrite', 'true')
    networkLogBlock.appendChild(createReferenceBlock(xmlBlock, 'caseToCaseTransmissionLikelihood', farmsetid))    
    mcmcBlock.appendChild(networkLogBlock)
    parameterLogBlock = xmlBlock.createElement('log')
    parameterLogBlock.setAttribute('id', 'fileLog')
    parameterLogBlock.setAttribute('logEvery', str(logEvery))
    parameterLogBlock.setAttribute('fileName', parameterLogFileName)
    parameterLogBlock.setAttribute('overwrite', 'true')
    parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'posterior', 'posterior'))
    parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'prior', 'prior'))
    parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'likelihood', 'likelihood'))
    if estimateIncubationPeriods:
        if incubationDistribution=="T" or incubationDistribution=="L":
            parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_'+incubationDistribution+'N_mean'))
            parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_'+incubationDistribution+'N_stdev'))
        elif incubationDistribution=="G":
            parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_scale'))
            parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_shape'))
            
    mcmcBlock.appendChild(parameterLogBlock)
    screenLogBlock = xmlBlock.createElement('log')
    screenLogBlock.setAttribute('id', 'screenLog')
    screenLogBlock.setAttribute('logEvery', str(logEvery))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'posterior', 'Posterior', 'posterior', 4, 12))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'prior', 'Prior', 'prior', 4, 12))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'likelihood', 'Likelihood', 'likelihood', 4, 12))
    if estimateIncubationPeriods:
        if incubationDistribution=="T" or incubationDistribution=="L":
            screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'parameter', 'Inc. mean', 'inc_'+incubationDistribution+'N_mean', 4, 12))
            screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'parameter', 'Inc. stdev', 'inc_'+incubationDistribution+'N_stdev', 4, 12))
        elif incubationDistribution=="G":
            screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'parameter', 'Inc. scale', 'inc_gamma_scale', 4, 12))
            screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'parameter', 'Inc. shape', 'inc_gamma_shape', 4, 12))
    mcmcBlock.appendChild(screenLogBlock)

    reportBlock = xmlBlock.createElement('report')
    beastElement.appendChild(reportBlock)
    reportPropertyBlock = xmlBlock.createElement('property')
    reportPropertyBlock.setAttribute('name', 'timer')
    reportPropertyBlock.appendChild(createReferenceBlock(xmlBlock, 'mcmc', 'mcmc'))
    reportBlock.appendChild(reportPropertyBlock)
    
    xmlBlock.appendChild(beastElement)
    
    # Write to file
    xmlFile = open(outputFileName, 'w')
    xmlFile.write(xmlBlock.toprettyxml('  ', newl = '\n', encoding = 'utf-8'))
    
# Create a block with an enclosing element with name "name", and a parameter element within it with id "id". If greaterThanZero then it has a lower bound of zero.
    
def createParameterBlock(xmlBlock, name, idstring, value, greaterThanZero):
    enclosingElement = xmlBlock.createElement(name)
    parameterElement = xmlBlock.createElement('parameter')
    parameterElement.setAttribute('value', str(value))
    parameterElement.setAttribute('id', idstring)
    if greaterThanZero:
        parameterElement.setAttribute('lower', '0.0')
    enclosingElement.appendChild(parameterElement)
    return enclosingElement

# Create a block with an enclosing element with name "name" and a gamma distribution element within it.

def createScreenLogColumnBlock(xmlBlock, name, label, reference, dp, width):
    columnElement = xmlBlock.createElement('column')
    columnElement.setAttribute('label', label)
    columnElement.setAttribute('dp', str(dp))
    columnElement.setAttribute('width', str(width))
    innerElement = xmlBlock.createElement(name)
    innerElement.setAttribute('idref', reference)
    columnElement.appendChild(innerElement)
    return columnElement

def createGammaBlock(xmlBlock, name, idstring, scale, shape, scaleID, shapeID):
    enclosingElement = xmlBlock.createElement(name)
    gammaElement = xmlBlock.createElement('gammaDistributionModel')
    gammaElement.appendChild(createParameterBlock(xmlBlock, 'shape', shapeID, shape, True))
    gammaElement.appendChild(createParameterBlock(xmlBlock, 'scale', scaleID, scale, True))
    gammaElement.setAttribute('id', idstring)
    enclosingElement.appendChild(gammaElement)
    return enclosingElement

def createUpperTruncNormalBlock(xmlBlock, name, idstring, mean, stdev, maximum, meanID, sdID):
    enclosingElement = xmlBlock.createElement(name)
    truncNormalElement = xmlBlock.createElement('truncatedNormalDistributionModel')
    truncNormalElement.appendChild(createParameterBlock(xmlBlock, 'mean', meanID, mean, False))
    truncNormalElement.appendChild(createParameterBlock(xmlBlock, 'stdev', sdID, stdev, True))
    truncNormalElement.appendChild(createParameterBlockNoID(xmlBlock, 'maximum', maximum, False))
    truncNormalElement.setAttribute('id', idstring)
    enclosingElement.appendChild(truncNormalElement)
    return enclosingElement

def createLowerTruncNormalBlock(xmlBlock, name, idstring, mean, stdev, minimum, meanID, sdID):
    enclosingElement = xmlBlock.createElement(name)
    truncNormalElement = xmlBlock.createElement('truncatedNormalDistributionModel')
    truncNormalElement.appendChild(createParameterBlock(xmlBlock, 'mean', meanID, mean, False))
    truncNormalElement.appendChild(createParameterBlock(xmlBlock, 'stdev', sdID, stdev, True))
    truncNormalElement.appendChild(createParameterBlockNoID(xmlBlock, 'minimum', minimum, False))
    truncNormalElement.setAttribute('id', idstring)
    enclosingElement.appendChild(truncNormalElement)
    return enclosingElement

def createLogNormalBlock(xmlBlock, name, idstring, mean, stdev, meanID, sdID):
    enclosingElement = xmlBlock.createElement(name)
    logNormalElement = xmlBlock.createElement('logNormalDistributionModel')
    logNormalElement.setAttribute('meanInRealSpace','true')
    logNormalElement.setAttribute('stdevInRealSpace','true')   
    logNormalElement.appendChild(createParameterBlock(xmlBlock, 'mean', meanID, mean, False))
    logNormalElement.appendChild(createParameterBlock(xmlBlock, 'stdev', sdID, stdev, True))
    logNormalElement.setAttribute('id', idstring)
    enclosingElement.appendChild(logNormalElement)
    return enclosingElement
    

# Create a block with an enclosing element with name "name", and a parameter element within it. If greaterThanZero then it has a lower bound of zero.

def createParameterBlockNoID(xmlBlock, name, value, greaterThanZero):
    enclosingElement = xmlBlock.createElement(name)
    parameterElement = xmlBlock.createElement('parameter')
    parameterElement.setAttribute('value', str(value))
    if greaterThanZero:
        parameterElement.setAttribute('lower', '0.0')
    enclosingElement.appendChild(parameterElement)
    return enclosingElement

# Create a block with an enclosing element with name "name", and a date element within it with id "id". 

def createEnclosedDateBlock(xmlBlock, name, idstring, value, direction, units, origin):
    enclosingElement = xmlBlock.createElement(name)
    dateElement = xmlBlock.createElement('date')
    dateElement.setAttribute('value', str(value))
    dateElement.setAttribute('id', idstring)
    dateElement.setAttribute('direction', direction)
    dateElement.setAttribute('units', units)
    dateElement.setAttribute('origin', str(origin))
    enclosingElement.appendChild(dateElement)
    return enclosingElement

# Create a block with an enclosing element with name "name", and a date element within it. 

def createEnclosedDateBlockNoID(xmlBlock, name, value, direction, units, origin):
    enclosingElement = xmlBlock.createElement(name)
    dateElement = xmlBlock.createElement('date')
    dateElement.setAttribute('value', str(value))
    dateElement.setAttribute('direction', direction)
    dateElement.setAttribute('units', units)
    dateElement.setAttribute('origin', str(origin))
    enclosingElement.appendChild(dateElement)
    return enclosingElement

# Create a block with an enclosing element with name "name" and an integer within it.

def createIntegerBlock(xmlBlock, name, value):
    enclosingElement = xmlBlock.createElement(name)
    integerElement = xmlBlock.createElement('integer')
    integerValue = xmlBlock.createTextNode(str(value))
    integerElement.appendChild(integerValue)
    enclosingElement.appendChild(integerElement)
    return enclosingElement

def createStringBlock(xmlBlock, name, value):
    enclosingElement = xmlBlock.createElement(name)
    stringElement = xmlBlock.createElement('string')
    stringValue = xmlBlock.createTextNode(str(value))
    stringElement.appendChild(stringValue)
    enclosingElement.appendChild(stringElement)
    return enclosingElement

def createIntegerBlockID(xmlBlock, name, myid, value):
    enclosingElement = xmlBlock.createElement(name)
    integerElement = xmlBlock.createElement('integer')
    integerElement.setAttribute('id',myid)
    integerValue = xmlBlock.createTextNode(str(value))
    integerElement.appendChild(integerValue)
    enclosingElement.appendChild(integerElement)
    return enclosingElement

# Create a block of name "name" and an idref to "reference".

def createReferenceBlock(xmlBlock, name, reference):
    refBlock = xmlBlock.createElement(name)
    refBlock.setAttribute('idref', reference)
    return refBlock

# Create an enclosed reference block; "parent" is the enclosing block.

def createNestedReferenceBlock(xmlBlock, parentName, childName, reference):
    parentBlock = xmlBlock.createElement(parentName)
    parentBlock.appendChild(createReferenceBlock(xmlBlock, childName, reference))
    return parentBlock

# Create a farm element

def createLesionDatedFarmCaseElement(xmlBlock, caseID, examinationDay, cullDay, oldestLesionAge, dayOne, standardDevInDays, taxa, incubationType):
    farmElement = xmlBlock.createElement('lesionDatedFarmCase')
    farmElement.setAttribute('caseID', caseID)
    farmElement.setAttribute('id', 'Farm'+caseID)
    farmElement.appendChild(createEnclosedDateBlockNoID(xmlBlock, 'examinationDay', examinationDay, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    farmElement.appendChild(createEnclosedDateBlockNoID(xmlBlock, 'cullDay', cullDay, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    farmElement.appendChild(createParameterBlockNoID(xmlBlock, 'oldestLesionAge', oldestLesionAge, True))
    # examinationDay-oldestLesionAge is the best guess for the time of infectiousness
    farmElement.appendChild(createUpperTruncNormalBlock(xmlBlock, 'infectiousDateDistribution', 'farm'+caseID+'InfectiousDate', examinationDay-oldestLesionAge, standardDevInDays, examinationDay, caseID+"_mean", caseID+"_stdev"))
    farmElement.appendChild(createNestedReferenceBlock(xmlBlock, 'incubationPeriodDistribution', incubationType, 'overallIncubationPeriodDistribution'))
    for taxon in taxa:
        farmElement.appendChild(createReferenceBlock(xmlBlock, 'taxon', taxon))
    return farmElement
    

def pickGammaParameters(mode, variance):
    scale = (-mode+math.sqrt(math.pow(mode,2)+(4*variance)))/2
    shape = variance/math.pow(scale,2)    
    return [shape,scale]
    
# Turn a date to a fractional year. Probably deprecated at this point.
    
def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def argOrNone(argument):
    try:
        variable = argument
    except:
        variable = None
    return variable



# Main method; command line entry for data that isn't in the CSV file. If the command has a second argument, it is the XML filename to be written to. Otherwise, it writes to "eieio.xml" ;-) 

def main():
    startingNetworkFileName = None
    dayOne = None
    estimateIncubationPeriods = None
    incubationDistribution = None
    incubationMean = None
    incubationStdev = None
    standardDevInDays = None
    riemannSteps = None
    chainLength = None
    logEvery = None
    
    parser = argparse.ArgumentParser(description='Write XML for the CaseToCaseTransmissionLikelihood class in BEAST', prog="oldMacDonald")
    parser.add_argument('epiFile', help='The path of a CSV file of epidemiological information; the file must have columns "Farm_ID", "Taxon", "Exam_date", "Cull_date" and "Oldest_lesion"')
    parser.add_argument('treeFile', help='The path of a Newick file containing the single fixed tree')
    parser.add_argument('outputXMLFile', help='The name of the output XML file')
    parser.add_argument('networkOut', help='The file that BEAST will write the network output to')
    parser.add_argument('paramOut', help='The file that BEAST will write the parameter output to')
    parser.add_argument('-d', '--dayOne', help='A date to be used as day 1 of the model (before any plausible TMRCA, use dd/mm/yyyy format)')
    parser.add_argument('-i', '--estimateInc', help='Estimate the parameters of the incubation period distribution simultaneously', type=int)
    parser.add_argument('-D', '--incubationPeriodType', help='The type of probability distribution for the incubation period (T=truncated normal (>0), G=gamma, L=lognormal')
    parser.add_argument('-M', '--incubationMean', help='Mean of incubation period truncated normal distribution (the incubation period is in days)', type=float)
    parser.add_argument('-S', '--incubationSD', help='Standard deviation of incubation period truncated normal distribution (the incubation period is in days)', type=float)
    parser.add_argument('-L', '--lesionDatingSD', help='The estimated standard deviation of error in lesion dating (in days)', type=float)
    parser.add_argument('-r', '--riemannSteps', help='The number of Riemann numerical integrator steps to use', type=int)
    parser.add_argument('-l', '--chainLength', help='The length of the MCMC chain', type=int)
    parser.add_argument('-s', '--samplingFreq', help='The sampling frequency of MCMC states', type=int)
    parser.add_argument('-N', '--startingNetwork', help='A CSV file describing a starting network; should have a header line followed by a line for each taxon giving its ID in column 1 and the ID of its parent in column 2. If a node has no parent, column 2 should be "Start".')
    
    arguments = parser.parse_args()


    
    rawName = arguments.epiFile
    #    At some point you're going to have to make this check the CSV file for the correct column headings. But not now.       
    csvreader = csv.reader(open(rawName, 'rU'))
    treeFile = arguments.treeFile
    treeFileReader = open(treeFile, 'rU')
    tree = treeFileReader.read()
    outputFileName = arguments.outputXMLFile
    networkFileName = arguments.networkOut
    parameterFileName = arguments.paramOut
    try:
        startingNetworkFileName = arguments.startingNetwork
    except:
        startingNetworkFileName = None
    try:
        dayOne = arguments.dayOne  
    except:
        while dayOne==None:
            try:
                dayOne = raw_input('Enter a date to be used as day 1 of the model (before any plausible TMRCA, use dd/mm/yyyy format): ')
            except:
                print "Please enter a valid date (dd/mm/yy format)"
    try:
        estimateIncubationPeriods = bool(arguments.estimateInc)
    except:
        estimateQ = raw_input('Simultaneously estimate parameters of incubation period gamma distribution? (Y/N): ')
        if estimateQ == "Y":
            estimateIncubationPeriods = True
    if estimateIncubationPeriods:
        incubationDistribution = arguments.incubationPeriodType
        if incubationDistribution == None:
            incubationDistribution = "T"
    else:
        try:
            incubationMean = arguments.incubationMean
        except:    
            while incubationMean == None:
                try:
                    incubationMean = float(raw_input('Enter mean of incubation period (truncated) normal distribution (incubation period in days): '))
                except:
                    print "Please enter a valid decimal number"  
        try:
            incubationStdev = arguments.incubationSD
        except:
            while incubationStdev == None:
                try:
                    incubationStdev = float(raw_input('Enter standard deviation of incubation period (truncated) normal distribution (incubation period in days): '))
                except:
                    print "Please enter a valid decimal number"        
    try:
        standardDevInDays = arguments.lesionDatingSD
    except:
        while standardDevInDays == None:
            try:
                standardDevInDays = float(raw_input('Enter estimated standard deviation of error in lesion dating (in days): '))
            except:
                print "Please enter a valid decimal number"
    try:
        riemannSteps = arguments.riemannSteps
    except:
        while riemannSteps == None:
            try:
                riemannSteps = float(raw_input('Enter number of Riemann numerical integrator steps to use: '))
            except:
                print "Please enter a valid decimal number"
    try:
        chainLength = arguments.chainLength
    except:
        while chainLength == None:
            try:
                chainLength = int(raw_input('Enter length of MCMC chain: '))
            except:
                print "Please enter an integer > 0"
    try:
        logEvery = arguments.samplingFreq
    except:
        while logEvery == None:
            try:
                logEvery = int(raw_input('Enter the frequency to log states: '))
            except:
                print "Please enter an integer > 0"
    generateXML(csvreader, datetime.strptime(dayOne, '%d/%m/%Y'), incubationStdev, incubationMean, standardDevInDays,
                riemannSteps, outputFileName, chainLength, logEvery, tree,
                networkFileName, parameterFileName, startingNetworkFileName,
                estimateIncubationPeriods, incubationDistribution) 

if __name__ == '__main__':
    main()

