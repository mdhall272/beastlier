# Take a file generated by BEAUTi and a CSV file of farm information and make a caseToCase XML file for BEAST. Version 2.0

import csv
from lxml import etree
from datetime import datetime
import time
import math
import argparse

# Main method

def generateXML(csvreader, dayOne, riemannSteps, outputFileName, beautiFileName, fileNameRoot, sampleTTs, model, extended):
    parser = etree.XMLParser(remove_blank_text=True)
    
    fullTree = etree.parse(open(beautiFileName), parser)  
    
    beastElement = fullTree.getroot()
    
    # Remove the old taxa elements
    taxaElement = beastElement.find('taxa')
    for taxon in taxaElement:
        taxaElement.remove(taxon)
    
    operatorsElement = beastElement.find('operators')
    
    # Add the C2CL block
    c2cElement = etree.Element('caseToCaseTransmissionLikelihood')
    beastElement.insert(beastElement.index(operatorsElement)-1,c2cElement)  
    c2cElement.set('extended', str(extended).lower())
    c2cElement.set('sampleTTs', str(sampleTTs).lower())
    
    c2cElement.append(createNestedReferenceBlock(etree,'virusTree','treeModel','treeModel')) 
    if(model=='m'):     
        farmsetname = 'morelli12Outbreak'
        farmsetid = 'morelli12Cases'
        farmname = 'morelli12Case'
    elif(model=='s'):
        farmsetname = 'simpleOutbreak'
        farmsetid = 'simpleCases'
        farmname = 'simpleCase'
    else:
        raise BaseException('Unknown model')
    outbreakElement = etree.SubElement(c2cElement,farmsetname)
    c2cElement.set('id', farmsetid)
    
    # Information common to all farms
    
    if model=='m':
        outbreakElement.append(createGammaBlock(etree, 'latentPeriodDistribution', 'overallLatentPeriodDistribution', '0.5', '20', 'inc_gamma_scale', 'inc_gamma_shape'))
        
        meanProductElement = etree.SubElement(outbreakElement,'productStatistic')
        meanProductElement.set('id','latent_mean')
        meanProductElement.append(createReferenceBlock(etree, 'parameter', 'inc_gamma_shape'))
        meanProductElement.append(createReferenceBlock(etree, 'parameter', 'inc_gamma_scale'))
        
        stdevProductElement = etree.SubElement(outbreakElement,'productStatistic')
        stdevProductElement.set('id','latent_var')
        stdevProductElement.append(createReferenceBlock(etree, 'parameter', 'inc_gamma_shape'))
        stdevProductElement.append(createReferenceBlock(etree, 'parameter', 'inc_gamma_scale'))
        stdevProductElement.append(createReferenceBlock(etree, 'parameter', 'inc_gamma_scale'))

    outbreakElement.append(createParameterBlockNoID(etree, 'riemannSampleSize', int(riemannSteps), True, 1))
    outbreakElement.append(createParameterBlock(etree, 'sqrtInfectiousScale', 'sqrt_inf_scale', "0.5", True, 1))
  
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
        taxonElement = etree.SubElement(taxaElement,'taxon')
        taxonElement.set('id', currentRow[taxonColumn])
        samplingDateElement = etree.SubElement(taxonElement,'date')
        # Events (sample and cull) now happen at the _end_ of the day. Easier that way. Hopefully. Still need to add 1, since time 0 is the beginning
        # of day 0, and thus if you subtract day 0's time value from another you get the number of days elapsed between the ends.
        samplingDateElement.set('value', str(intExamDate+1))
        samplingDateElement.set('direction','forwards')
        samplingDateElement.set('units','days')
        samplingDateElement.set('origin', dayOne.strftime('%d/%m/%Y'))
        # Want to be adding +1 to dates here to get the end of the required days.
        outbreakElement.append(createCaseElement(etree, farmname, currentRow[farmIDColumn], intExamDate+1, intCullDate+1, 
                                                         float(currentRow[oldestLesionColumn]), dayOne, 'sqrt_inf_scale', taxa))
        try:
            currentRow=csvreader.next()
        except StopIteration:
            currentRow=None
    
    
    operatorsBlock = beastElement.find('operators')
    if sampleTTs:
        treeOperators = {'narrowExchange', 'wideExchange', 'wilsonBalding', 'subtreeSlide'}
        for operator in operatorsBlock.iterchildren():
            if operator.tag in treeOperators:
                convertOperatorElement(etree, operatorsBlock, operator, farmsetid)
        branchLSOElement = etree.SubElement(operatorsBlock,'nodePaintingSwitchOperator')
        branchLSOElement.set('weight', str(5))
        ftlElement = etree.SubElement(branchLSOElement,'caseToCaseTransmissionLikelihood')
        ftlElement.set('idref',farmsetid)
    if model=='m':
        incShapeOperator = etree.SubElement(operatorsBlock,'scaleOperator')
        incScaleOperator = etree.SubElement(operatorsBlock,'scaleOperator')
        incShapeOperator.set('scaleFactor', '0.75')
        incShapeOperator.set('weight', '5')
        incShapeOperator.set('autoOptimize', 'true')
        incScaleOperator.set('scaleFactor', '0.75')
        incScaleOperator.set('weight', '5')
        incScaleOperator.set('autoOptimize', 'true')
    dOperator = etree.SubElement(operatorsBlock,'scaleOperator')
    dOperator.set('scaleFactor', '0.75')
    dOperator.set('weight', '5')
    if model=='m':
        incShapeOperator.append(createReferenceBlock(etree, 'parameter', 'inc_gamma_shape'))
        incScaleOperator.append(createReferenceBlock(etree, 'parameter', 'inc_gamma_scale'))
    dOperator.append(createReferenceBlock(etree, 'parameter', 'sqrt_inf_scale'))
#     if doGeography:
#         alphasOperator = etree.SubElement(operatorsBlock, 'scaleOperator')
#         alphasOperator.set('scaleFactor', '0.75')
#         alphasOperator.set('weight', '1')
#         alphasOperator.set('autoOptimize', 'true')
#         alphasOperator.set('scaleAllIndependently', 'true')
#         alphasOperator.append(createReferenceBlock(etree, 'parameter', 'kernelAlphas'))
    
    
            
    mcmcBlock = beastElement.find('mcmc')
    posteriorBlock = mcmcBlock.find('posterior')
    priorBlock = posteriorBlock.find('prior')
    if model=='m':
        meanPrior = etree.SubElement(priorBlock,'gammaPrior')
        meanPrior.set('shape', '1')
        meanPrior.set('scale', '1')    
        meanPrior.append(createReferenceBlock(etree, 'statistic', 'latent_mean'))
#     varPrior = etree.SubElement(priorBlock,'uniformPrior')
#     varPrior.set('lower', '0')
#     varPrior.set('upper', '6.25')
#     varPrior.append(createReferenceBlock(etree, 'statistic', 'latent_var'))
#     if doGeography:
#         alphasPrior = etree.SubElement(priorBlock,'exponentialPrior')
#         alphasPrior.set('mean', '100')
#         alphasPrior.set('offset', '0')
#         alphasPrior.append(createReferenceBlock(etree, 'parameter', 'kernelAlphas'))
#     
#     dPrior = etree.SubElement(priorBlock,'uniformPrior')
#     dPrior.set('lower', '0')
#     dPrior.set('upper', '2')
#     dPrior.append(createReferenceBlock(etree, 'parameter', 'sqrt_inf_scale'))
                                  
    likelihoodBlock = posteriorBlock.find('likelihood')
    likelihoodBlock.append(createReferenceBlock(etree, 'caseToCaseTransmissionLikelihood', farmsetid))
    
    for logBlock in mcmcBlock.iterchildren('log'):
        if logBlock.get('id')=='screenLog':
            if model=='m':
                logBlock.append(createScreenLogColumnBlock(etree, 'productStatistic', 'Inc. mean', 'latent_mean', 4, 12))
                logBlock.append(createScreenLogColumnBlock(etree, 'productStatistic', 'Inc. var', 'latent_var', 4, 12))           
            logBlock.append(createScreenLogColumnBlock(etree, 'parameter', 'd', 'sqrt_inf_scale', 4, 12))
        elif logBlock.get('id')=='fileLog':
            logBlock.set('fileName',fileNameRoot+".log.txt")
            logEvery = logBlock.get('logEvery')
            if model=='m':
                logBlock.append(createReferenceBlock(etree, 'parameter', 'inc_gamma_scale'))
                logBlock.append(createReferenceBlock(etree, 'parameter', 'inc_gamma_shape'))
                logBlock.append(createReferenceBlock(etree, 'productStatistic', 'latent_mean'))
                logBlock.append(createReferenceBlock(etree, 'productStatistic', 'latent_var'))
            logBlock.append(createReferenceBlock(etree, 'parameter', 'sqrt_inf_scale'))

#             if doGeography:
#                 logBlock.append(createReferenceBlock(etree, 'parameter', 'kernelAlphas'))
            networkLogBlock = etree.Element('log')
            networkLogBlock.set('id', 'networkLog')
            networkLogBlock.set('logEvery', logEvery)
            networkLogBlock.set('fileName', fileNameRoot+".net.txt")
            networkLogBlock.set('overwrite', 'true')
            networkLogBlock.append(createReferenceBlock(etree, 'caseToCaseTransmissionLikelihood', farmsetid))
            mcmcBlock.insert(mcmcBlock.index(logBlock)+1, networkLogBlock)
    
    treeLogBlock = mcmcBlock.find('logTree')
    treeLogBlock.set('fileName', fileNameRoot+".trees.txt")
    firstTrait = treeLogBlock.find('trait')
    paintingTrait = etree.Element('trait')
    paintingTrait.set('name', 'paintings')
    paintingTrait.set('tag', 'paintings')
    paintingTrait.append(createReferenceBlock(etree, 'caseToCaseTransmissionLikelihood', farmsetid))
    treeLogBlock.insert(treeLogBlock.index(firstTrait), paintingTrait)
                  
    # Write to file

    et = etree.ElementTree(beastElement)
    
    
    et.write(outputFileName, pretty_print=True)
    
# Create a block with an enclosing element with name "name", and a parameter element within it with id "id". If greaterThanZero then it has a lower bound of zero.
    
def createParameterBlock(tree, name, idstring, value, greaterThanZero, dim):
    enclosingElement = tree.Element(name)
    parameterElement = tree.SubElement(enclosingElement,'parameter')
    parameterElement.set('value', str(value))
    parameterElement.set('id', idstring)
    if greaterThanZero:
        parameterElement.set('lower', '0.0')
    if dim>1:
        parameterElement.set('dimension', str(dim))
    return enclosingElement

# Create a block with an enclosing element with name "name" and a gamma distribution element within it.

def createScreenLogColumnBlock(tree, name, label, reference, dp, width):
    columnElement = tree.Element('column')
    columnElement.set('label', label)
    columnElement.set('dp', str(dp))
    columnElement.set('width', str(width))
    innerElement = tree.SubElement(columnElement,name)
    innerElement.set('idref', reference)
    return columnElement

def createGammaBlock(tree, name, idstring, scale, shape, scaleID, shapeID):
    enclosingElement = tree.Element(name)
    gammaElement = tree.SubElement(enclosingElement,'gammaDistributionModel')
    gammaElement.append(createParameterBlock(tree, 'shape', shapeID, shape, True, 1))
    gammaElement.append(createParameterBlock(tree, 'scale', scaleID, scale, True, 1))
    gammaElement.set('id', idstring)
    return enclosingElement

def createUpperTruncNormalBlock(xmlBlock, name, idstring, mean, stdev, maximum, meanID, sdID):
    enclosingElement = xmlBlock.createElement(name)
    truncNormalElement = xmlBlock.createElement('truncatedNormalDistributionModel')
    truncNormalElement.appendChild(createParameterBlock(xmlBlock, 'mean', meanID, mean, False, 1))
    truncNormalElement.appendChild(createParameterBlock(xmlBlock, 'stdev', sdID, stdev, True, 1))
    truncNormalElement.appendChild(createParameterBlockNoID(xmlBlock, 'maximum', maximum, False, 1))
    truncNormalElement.setAttribute('id', idstring)
    enclosingElement.appendChild(truncNormalElement)
    return enclosingElement

def createLowerTruncNormalBlock(xmlBlock, name, idstring, mean, stdev, minimum, meanID, sdID):
    enclosingElement = xmlBlock.createElement(name)
    truncNormalElement = xmlBlock.createElement('truncatedNormalDistributionModel')
    truncNormalElement.appendChild(createParameterBlock(xmlBlock, 'mean', meanID, mean, False, 1))
    truncNormalElement.appendChild(createParameterBlock(xmlBlock, 'stdev', sdID, stdev, True, 1))
    truncNormalElement.appendChild(createParameterBlockNoID(xmlBlock, 'minimum', minimum, False, 1))
    truncNormalElement.setAttribute('id', idstring)
    enclosingElement.appendChild(truncNormalElement)
    return enclosingElement

def createLogNormalBlock(xmlBlock, name, idstring, mean, stdev, meanID, sdID):
    enclosingElement = xmlBlock.createElement(name)
    logNormalElement = xmlBlock.createElement('logNormalDistributionModel')
    logNormalElement.setAttribute('meanInRealSpace','true')
    logNormalElement.setAttribute('stdevInRealSpace','true')   
    logNormalElement.appendChild(createParameterBlock(xmlBlock, 'mean', meanID, mean, False, 1))
    logNormalElement.appendChild(createParameterBlock(xmlBlock, 'stdev', sdID, stdev, True, 1))
    logNormalElement.setAttribute('id', idstring)
    enclosingElement.appendChild(logNormalElement)
    return enclosingElement
    

# Create a block with an enclosing element with name "name", and a parameter element within it. If greaterThanZero then it has a lower bound of zero.

def createParameterBlockNoID(tree, name, value, greaterThanZero, dim):
    enclosingElement = etree.Element(name)
    parameterElement = etree.SubElement(enclosingElement,'parameter')
    parameterElement.set('value', str(value))
    if greaterThanZero:
        parameterElement.set('lower', '0.0')
    if dim>1:
        parameterElement.set('dimension', str(dim))
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

def createEnclosedDateBlockNoID(tree, name, value, direction, units, origin):
    enclosingElement = tree.Element(name)
    dateElement = tree.SubElement(enclosingElement,'date')
    dateElement.set('value', str(value))
    dateElement.set('direction', direction)
    dateElement.set('units', units)
    dateElement.set('origin', str(origin))
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

def createReferenceBlock(tree, name, reference):
    refBlock = tree.Element(name)
    refBlock.set('idref', reference)
    return refBlock

# Create an enclosed reference block; "parent" is the enclosing block.

def createNestedReferenceBlock(tree, parentName, childName, reference):
    parentBlock = tree.Element(parentName)
    parentBlock.append(createReferenceBlock(tree, childName, reference)) 
    return parentBlock

# Create a farm element

def createCaseElement(tree, name, caseID, examinationDay, cullDay, oldestLesionAge, dayOne, sqrtInfScaleName, taxa):
    caseElement = tree.Element(name)
    caseElement.set('caseID', caseID)
    caseElement.set('id', 'Farm'+caseID)
    caseElement.append(createEnclosedDateBlockNoID(tree, 'examinationDay', examinationDay, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    caseElement.append(createEnclosedDateBlockNoID(tree, 'cullDay', cullDay, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    caseElement.append(createParameterBlockNoID(tree, 'oldestLesionAge', oldestLesionAge, True, 1))
                                                   
    # examinationDay-oldestLesionAge is the best guess for the time of infectiousness
    for taxon in taxa:
        caseElement.append(createReferenceBlock(tree, 'taxon', taxon))
    return caseElement

def convertOperatorElement(tree, operatorsElement, operatorElement, c2cReference):
    operatorsElement.remove(operatorElement)
    ttoElement = tree.SubElement(operatorsElement, 'transmissionTreeOperator')
    ttoElement.append(operatorElement)
    ttoElement.append( createReferenceBlock(tree, 'caseToCaseTransmissionLikelihood', c2cReference))
    

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
    dayOne = None
    riemannSteps = None
    sampleTTs = False
    model = 'm'
    extended = False

    parser = argparse.ArgumentParser(description='Write XML for the CaseToCaseTransmissionLikelihood class in BEAST', prog="oldMacDonald")
    parser.add_argument('beautiFile', help='The path of an existing XML file generated by BEAUTi or similar, ready to have the elements for transmission tree reconstruction added')
    parser.add_argument('epiFile', help='The path of a CSV file of epidemiological information; the file must have columns "Farm_ID", "Taxon", "Exam_date", "Cull_date" and "Oldest_lesion"')
    parser.add_argument('outputXMLFile', help='The name of the output XML file')
    parser.add_argument('fileNameRoot', help='The root of the file names that BEAST will write the output to')
    parser.add_argument('-d', '--dayOne', help='A date to be used as day 1 of the model (before any plausible TMRCA, use dd/mm/yyyy format)')
    parser.add_argument('-r', '--riemannSteps', help='The number of Riemann numerical integrator steps to use', type=int)
    parser.add_argument('-s', '--sampleTTs', help='Sample transmission trees individually')
    parser.add_argument('-m', '--model', help='Case model (m=Morelli12Outbreak, s=SimpleOutbreak)')
    parser.add_argument('-e', '--extended', help='Extended version (TT>=TMRCA)')
    
    arguments = parser.parse_args()


    
    rawName = arguments.epiFile
    #    At some point you're going to have to make this check the CSV file for the correct column headings. But not now.       
    csvreader = csv.reader(open(rawName, 'rU'))
    outputFileName = arguments.outputXMLFile
    fileNameRoot = arguments.fileNameRoot
    beautiFileName = arguments.beautiFile
    model = arguments.model
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
        sampleTTs = bool(arguments.sampleTTs=="True") 
    except:
        estimateF = raw_input('Use the Felsenstein pruning algorithm to calculate overall likelihood of the phylogenetic tree under all transmission trees? ')
        if estimateF == "Y":
            sampleTTs = True
    try:
        extended = bool(arguments.extended=="True") 
    except:
        estimateE = raw_input('Use the extended version of tree-paintings (TT>=TMRCA)?')
        if estimateE == "Y":
            extended = True
    generateXML(csvreader, datetime.strptime(dayOne, '%d/%m/%Y'), riemannSteps, outputFileName, beautiFileName, fileNameRoot, sampleTTs, model, extended) 

if __name__ == '__main__':
    main()

