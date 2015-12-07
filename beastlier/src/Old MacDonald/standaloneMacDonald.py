# Turn a CSV file of farm information into XML blocks for BEAST. Version 3.0,

import csv
import xml.dom.minidom as mini
from datetime import datetime
import time
import math
import argparse

# Main method

def generateXML(csvreader, dayOne, riemannSteps, outputFileName, chainLength, logEvery, tree, networkLogFileName, parameterLogFileName,
                startingNetworkFileName, felsenstein):
    xmlBlock = mini.Document()
    beastElement = xmlBlock.createElement('beast')   
    taxaComment = xmlBlock.createComment('Taxa block.')
    taxaElement = xmlBlock.createElement('taxa')
    taxaElement.setAttribute('id','taxa')
    likelihoodComment = xmlBlock.createComment('CaseToCaseLikelihood block.')
    rootElement = xmlBlock.createElement('caseToCaseTransmissionLikelihood')
    rootElement.setAttribute('extended', 'false')
    rootElement.setAttribute('felsenstein', str(felsenstein).lower())
    rootElement.appendChild(createNestedReferenceBlock(xmlBlock,'virusTree','newick','startingTree'))      
    farmsetname = 'morelli12Outbreak'
    farmsetid = 'morelli12Cases'
    farmname = 'morelli12Case'
    farmSetElement = xmlBlock.createElement(farmsetname)
    rootElement.setAttribute('id', farmsetid)
    # Information common to all farms
    
    farmSetElement.appendChild(createGammaBlock(xmlBlock, 'incubationPeriodDistribution', 'overallIncubationPeriodDistribution', '0.1', '100', 'inc_gamma_scale', 'inc_gamma_shape'))
    
    meanProductElement = xmlBlock.createElement('productStatistic')
    meanProductElement.setAttribute('id','incubation_mean')
    meanProductElement.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_shape'))
    meanProductElement.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_scale'))
    farmSetElement.appendChild(meanProductElement)
    
    stdevProductElement = xmlBlock.createElement('productStatistic')
    stdevProductElement.setAttribute('id','incubation_var')
    stdevProductElement.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_shape'))
    stdevProductElement.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_scale'))
    stdevProductElement.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_scale'))
    farmSetElement.appendChild(stdevProductElement)

    farmSetElement.appendChild(createParameterBlockNoID(xmlBlock, 'riemannSampleSize', int(riemannSteps), True))
    farmSetElement.appendChild(createParameterBlock(xmlBlock, 'sqrtInfectiousScale', 'sqrt_inf_scale', 1, True))
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
        elif headerRow[i]=='Clinical_guess':
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
        farmSetElement.appendChild(createFarmCaseElement(xmlBlock, farmname, currentRow[farmIDColumn], intExamDate+1, intCullDate+1, 
                                                         float(currentRow[oldestLesionColumn]), dayOne, 'sqrt_inf_scale', taxa))
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
    
    
    operatorComment = xmlBlock.createComment('Operators')
    operatorsBlock = xmlBlock.createElement('operators')
    operatorsBlock.setAttribute('id','operators')
    if not felsenstein:
        branchLSOElement = xmlBlock.createElement('nodePaintingSwitchOperator')
        branchLSOElement.setAttribute('weight', str(1))
        operatorsBlock.appendChild(branchLSOElement)
        ftlElement = xmlBlock.createElement('caseToCaseTransmissionLikelihood')
        ftlElement.setAttribute('idref',farmsetid)
        branchLSOElement.appendChild(ftlElement)
    beastElement.appendChild(operatorComment)
    beastElement.appendChild(operatorsBlock)
    incShapeOperator = xmlBlock.createElement('scaleOperator')
    incScaleOperator = xmlBlock.createElement('scaleOperator')
    dOperator = xmlBlock.createElement('scaleOperator')
    incShapeOperator.setAttribute('scaleFactor', '0.75')
    incShapeOperator.setAttribute('weight', '1')
    incShapeOperator.setAttribute('autoOptimize', 'true')
    incScaleOperator.setAttribute('scaleFactor', '0.75')
    incScaleOperator.setAttribute('weight', '1')
    incScaleOperator.setAttribute('autoOptimize', 'true')
    dOperator.setAttribute('scaleFactor', '0.75')
    dOperator.setAttribute('weight', '1')
    incShapeOperator.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_shape'))
    incScaleOperator.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_scale'))
    dOperator.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'sqrt_inf_scale'))
    operatorsBlock.appendChild(incShapeOperator)
    operatorsBlock.appendChild(incScaleOperator)
    operatorsBlock.appendChild(dOperator)
            
    
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
    meanPrior = xmlBlock.createElement('uniformPrior')
    meanPrior.setAttribute('lower', '5')
    meanPrior.setAttribute('upper', '10')
    meanPrior.appendChild(createReferenceBlock(xmlBlock, 'statistic', 'incubation_mean'))
    priorBlock.appendChild(meanPrior)
    varPrior = xmlBlock.createElement('uniformPrior')
    varPrior.setAttribute('lower', '0')
    varPrior.setAttribute('upper', '4')
    varPrior.appendChild(createReferenceBlock(xmlBlock, 'statistic', 'incubation_var'))
    priorBlock.appendChild(varPrior) 
    dPrior = xmlBlock.createElement('uniformPrior')
    dPrior.setAttribute('lower', '0')
    dPrior.setAttribute('upper', '2')
    dPrior.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'sqrt_inf_scale'))
    priorBlock.appendChild(dPrior) 
                                  
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
    parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_scale'))
    parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'inc_gamma_shape'))
    parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'parameter', 'sqrt_inf_scale'))
    parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'productStatistic', 'incubation_mean'))
    parameterLogBlock.appendChild(createReferenceBlock(xmlBlock, 'productStatistic', 'incubation_var'))
            
    mcmcBlock.appendChild(parameterLogBlock)
    screenLogBlock = xmlBlock.createElement('log')
    screenLogBlock.setAttribute('id', 'screenLog')
    screenLogBlock.setAttribute('logEvery', str(logEvery))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'posterior', 'Posterior', 'posterior', 4, 12))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'prior', 'Prior', 'prior', 4, 12))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'likelihood', 'Likelihood', 'likelihood', 4, 12))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'parameter', 'Inc. shape', 'inc_gamma_shape', 4, 12))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'parameter', 'Inc. scale', 'inc_gamma_scale', 4, 12))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'productStatistic', 'Inc. mean', 'incubation_mean', 4, 12))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'productStatistic', 'Inc. var', 'incubation_var', 4, 12))
    screenLogBlock.appendChild(createScreenLogColumnBlock(xmlBlock, 'parameter', 'd', 'sqrt_inf_scale', 4, 12))
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

def createFarmCaseElement(xmlBlock, name, caseID, examinationDay, cullDay, oldestLesionAge, dayOne, sqrtInfScaleName, taxa):
    farmElement = xmlBlock.createElement(name)
    farmElement.setAttribute('caseID', caseID)
    farmElement.setAttribute('id', 'Farm'+caseID)
    farmElement.appendChild(createEnclosedDateBlockNoID(xmlBlock, 'examinationDay', examinationDay, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    farmElement.appendChild(createEnclosedDateBlockNoID(xmlBlock, 'cullDay', cullDay, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    farmElement.appendChild(createParameterBlockNoID(xmlBlock, 'oldestLesionAge', oldestLesionAge, True))
    # examinationDay-oldestLesionAge is the best guess for the time of infectiousness
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
    parser.add_argument('-r', '--riemannSteps', help='The number of Riemann numerical integrator steps to use', type=int)
    parser.add_argument('-l', '--chainLength', help='The length of the MCMC chain', type=int)
    parser.add_argument('-s', '--samplingFreq', help='The sampling frequency of MCMC states', type=int)
    parser.add_argument('-f', '--felsenstein', help='Use the Felsenstein pruning algorithm to calculate overall likelihood of the phylogenetic tree under all transmission trees (warning - may be slow)')
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
    try:
        felsenstein = bool(arguments.felsenstein)
    except:
        estimateF = raw_input('Use the Felsenstein pruning algorithm to calculate overall likelihood of the phylogenetic tree under all transmission trees (warning - may be slow)? ')
        if estimateF == "Y":
            felsenstein = True
    generateXML(csvreader, datetime.strptime(dayOne, '%d/%m/%Y'), riemannSteps, outputFileName, chainLength, logEvery, tree,
                networkFileName, parameterFileName, startingNetworkFileName, felsenstein) 

if __name__ == '__main__':
    main()

