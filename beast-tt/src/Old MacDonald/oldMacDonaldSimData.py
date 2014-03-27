# Turn a CSV file of farm information into XML blocks for BEAST

import csv
import xml.dom.minidom as mini
from datetime import datetime
import time
import sys
import math

# Main method

def generateXML(csvreader, averageIncubationPeriod, minInfectionToExam,
                distribution, standardDevInDays, earliestPossibleInfectionDate,
                dayOne, gammaScale, gammaShape, outputFileName, logEvery,
                networkLogFileName):
    xmlBlock = mini.Document()
    beastElement = xmlBlock.createElement('beast')
    
    # Using cull dates as sampling dates in this case (because the trees do in the
    # current incarnation)
    
    taxaComment = xmlBlock.createComment('Taxa block.')
    taxaElement = xmlBlock.createElement('taxa')
    taxaElement.setAttribute('id','taxa')
    likelihoodComment = xmlBlock.createComment('FarmTransmissionLikelihood block.')
    rootElement = xmlBlock.createElement('farmTransmissionLikelihood')
    rootElement.appendChild(createNestedReferenceBlock(xmlBlock,'virusTree','newick','startingTree'))
    farmsetname=''
    farmsetid=''
    if distribution=='b':
        farmsetname = 'cottam2008FarmSet'
        farmsetid = 'cottam2008Farms'
    elif distribution=='g':
        farmsetname = 'cottam2008GammaModFarmSet'
        farmsetid = 'cottam2008GammaModFarms'
    farmSetElement = xmlBlock.createElement(farmsetname)
    rootElement.setAttribute('id', farmsetid)
    # Information common to all farms
    farmSetElement.appendChild(createParameterBlock(xmlBlock, 'averageIncubationPeriod', 'overallAverageIncubationPeriod', averageIncubationPeriod, True))
    farmSetElement.appendChild(createParameterBlock(xmlBlock, 'minInfectionToExam', 'overallMinInfectionToExam', minInfectionToExam, True))
    farmSetElement.appendChild(createParameterBlock(xmlBlock, 'standardDevInDays', 'overallStandardDevInDays', standardDevInDays, True))
    
    # Date explanations, part 1: The decimal time value in BEAST is the start of the day. dayOne begins at time zero and ends at time 1. You just want the date in the XML in a format BEAST can read. 
    if farmSetElement=='b':
        intEarliestPossibleInfectionDate = (earliestPossibleInfectionDate - dayOne).days
    farmSetElement.appendChild(createEnclosedDateBlock(xmlBlock, 'dayOne', 'dayOne', 0, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    if distribution=='b':   
        farmSetElement.appendChild(createEnclosedDateBlock(xmlBlock, 'earliestInfectionDate', 'overallEarliestInfectionDate', intEarliestPossibleInfectionDate, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    farmSetElement.appendChild(createGammaBlock(xmlBlock, 'incubationPeriodDistribution', 'overallIncubationPeriodDistribution', gammaScale, gammaShape))
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
        taxa=currentRow[taxonColumn].split('$')
        taxonElement = xmlBlock.createElement('taxon')
        taxonElement.setAttribute('id', currentRow[taxonColumn])
        taxaElement.appendChild(taxonElement)

        samplingDateElement = xmlBlock.createElement('date')
        samplingDateElement.setAttribute('value', str(currentRow[cullDateColumn]))
        samplingDateElement.setAttribute('direction','forwards')
        samplingDateElement.setAttribute('units','days')
        taxonElement.appendChild(samplingDateElement)
        processedExamDate = datetime.strptime(currentRow[examDateColumn], '%d/%m/%Y')
        processedCullDate = datetime.strptime(currentRow[cullDateColumn], '%d/%m/%Y')
        # Time zero is the start of dayOne. Events are assumed to occur at the END of a given day (at least until it occurs to me why this is flawed). Hence the +1 here.
        intExamDate = (processedExamDate - dayOne).days
        intCullDate = (processedCullDate - dayOne).days
        farmSetElement.appendChild(createFarmElement(xmlBlock, currentRow[farmIDColumn], intExamDate, intCullDate, currentRow[oldestLesionColumn], dayOne, taxa, distribution))
        try:
            currentRow=csvreader.next()
        except StopIteration:
            currentRow=None
    beastElement.appendChild(taxaComment)
    beastElement.appendChild(taxaElement)
    beastElement.appendChild(likelihoodComment)
    rootElement.appendChild(farmSetElement)
    beastElement.appendChild(rootElement)
    
    
    operatorComment = xmlBlock.createComment('BranchLocationSwitchOperator elmement; place with other operators')
    operatorsBlock = xmlBlock.createElement('operators')
    operatorsBlock.setAttribute('id','operators')
    branchLSOElement = xmlBlock.createElement('branchLocationSwitchOperator')
    branchLSOElement.setAttribute('weight', str(1))
    operatorsBlock.appendChild(branchLSOElement)
    ftlElement = xmlBlock.createElement('farmTransmissionLikelihood')
    ftlElement.setAttribute('idref','cottam2008Farms')
    branchLSOElement.appendChild(ftlElement)
    beastElement.appendChild(operatorComment)
    beastElement.appendChild(operatorsBlock)
   
    loggerComment = xmlBlock.createComment('Log element; this should not occur with the other logs, because that file would be mental')
    logBlock = xmlBlock.createElement('log')
    logBlock.setAttribute('id','networkLog')
    logBlock.setAttribute('logEvery', str(logEvery))
    logBlock.setAttribute('fileName', networkLogFileName)
    logBlock.setAttribute('overwrite', 'true')
    ftlElement2 = xmlBlock.createElement('farmTransmissionLikelihood')
    ftlElement2.setAttribute('idref','cottam2008Farms')
    logBlock.appendChild(ftlElement2)
    beastElement.appendChild(loggerComment)
    beastElement.appendChild(logBlock)
    
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

def createGammaBlock(xmlBlock, name, idstring, scale, shape):
    enclosingElement = xmlBlock.createElement(name)
    gammaElement = xmlBlock.createElement('gammaDistributionModel')
    gammaElement.appendChild(createParameterBlockNoID(xmlBlock, 'shape', shape, False))
    gammaElement.appendChild(createParameterBlockNoID(xmlBlock, 'scale', scale, False))
    gammaElement.setAttribute('id', idstring)
    enclosingElement.appendChild(gammaElement)
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
 
def createFarmElement(xmlBlock, id, examinationDay, cullDay, oldestLesionAge, dayOne, taxa, distribution):
    farmElement = xmlBlock.createElement('cottam2008Farm')
    farmElement.setAttribute('name', id)
    farmElement.setAttribute('id', 'Farm'+id)
    farmElement.appendChild(createNestedReferenceBlock(xmlBlock, 'averageIncubationPeriod', 'parameter', 'overallAverageIncubationPeriod'))
    farmElement.appendChild(createNestedReferenceBlock(xmlBlock, 'minInfectionToExam', 'parameter', 'overallMinInfectionToExam'))
    farmElement.appendChild(createNestedReferenceBlock(xmlBlock, 'standardDevInDays', 'parameter', 'overallStandardDevInDays'))
    if distribution=='b':
        farmElement.appendChild(createNestedReferenceBlock(xmlBlock, 'earliestInfectionDate', 'date', 'overallEarliestInfectionDate'))
    farmElement.appendChild(createEnclosedDateBlockNoID(xmlBlock, 'examinationDay', examinationDay, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    farmElement.appendChild(createEnclosedDateBlockNoID(xmlBlock, 'cullDay', cullDay, 'forwards', 'days', dayOne.strftime('%d/%m/%Y')))
    farmElement.appendChild(createParameterBlockNoID(xmlBlock, 'oldestLesionAge', oldestLesionAge, True))
    farmElement.appendChild(createNestedReferenceBlock(xmlBlock, 'incubationPeriodDistribution', 'gammaDistributionModel', 'overallIncubationPeriodDistribution'))
    for taxon in taxa:
        farmElement.appendChild(createReferenceBlock(xmlBlock, 'taxon', taxon))
    return farmElement

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

# If you want a beta distribution with a given mode and variance, this function gives you its parameters

def pickBetaParameters(desiredMode, desiredVariance):
    Z = desiredMode/(1-desiredMode)
    a1 = math.pow(Z,3)
    a2 = 3*math.pow(Z,2)
    a3 = 3*Z
    a4 = 1
    a = desiredVariance*(a1+a2+a3+a4)
    Y = ((2*desiredMode)-1)
    X = 1-desiredMode
    b1 = (-3*Y*math.pow(desiredMode,2))/math.pow(X,3)
    b2 = (-6*Y)*desiredMode/math.pow(X,2)
    b3 = (-3*Y)/X
    b4 = math.pow(Z,2)
    b5 = 2*Z
    b6 = 1
    b7 = -Z
    b = desiredVariance*(b1+b2+b3+b4+b5+b6) + b7
    c1 = (3*desiredMode*math.pow(Y,2))/math.pow(X,3)
    c2 = (3*math.pow(Y,2))/math.pow(X,2)
    c3 = ((-2*desiredMode*Y)/math.pow(X,2))
    c4 = ((-2*Y)/X)
    c5 = Y/X
    c = desiredVariance*(c1+c2+c3+c4) + c5
    d1 = -math.pow(Y/X,3)
    d2 = math.pow(Y/X,2)
    d = desiredVariance*(d1 + d2)
    options = numpy.roots([a,b,c,d])
    bestOption = float('inf')
    for i in range(0,3):
        if options[i]>1 and (options[i]*desiredMode - 2*desiredMode +1)/(1-desiredMode)>1 and options[i]<bestOption and scipy.imag(options[i])==0:
            bestOption = options[i]
    if(bestOption==float('inf')):
        raise Exception('No suitable beta parameters')
    else:
        output=([0,0])
        output[0]=(bestOption*desiredMode - 2*desiredMode + 1)/(1-desiredMode)
        output[1]=bestOption
    return output

# Main method; command line entry for data that isn't in the CSV file. If the command has a second argument, it is the XML filename to be written to. Otherwise, it writes to "eieio.xml" ;-) 

def main():
    if len(sys.argv)>1:
        rawName = sys.argv[1]
        try:
            csvreader = csv.reader(open(rawName, 'rU'))
        except IOError:
            print "File not found"
    else:
        validFile = False
        while validFile==False:
            rawName = raw_input('Enter path of CSV file of epidemiological information:')
            try:
                csvreader = csv.reader(open(rawName, 'rU'))
                validFile=True
            except IOError:
                print "File not found"
    if len(sys.argv)>2:
        outputFileName = sys.argv[2]
    else:
        outputFileName = 'eieio.xml'
#    At some point you're going to have to make this check the CSV file for the correct column headings. But not now.            
    averageIncubationPeriod = ''
    while averageIncubationPeriod =='':
        try:
            averageIncubationPeriod = float(raw_input('Enter estimated mean incubation period in days: '))
        except:
            print "Please enter a valid decimal number"
    minInfectionToExam= ''
    while minInfectionToExam =='':
        try:
            minInfectionToExam = float(raw_input('Enter minimum time period in days that could elapse between infection and a positive examination: '))
        except:
            print "Please enter a valid decimal number"
    distribution=''
    while distribution !='g' and distribution!='b':
        try:
            distribution=raw_input('Distribution to use for date of infection (b=beta, g=gamma): ')
        except:
            print "Please enter \'b\' or \'g\'"
    standardDevInDays=''
    while standardDevInDays == '':
        try:
            standardDevInDays = float(raw_input('Enter standard deviation of estimated date of infection period distribution in days: '))
        except:
            print "Please enter a valid decimal number"
    earliestInfDate = ''
    if distribution=='b':
        while earliestInfDate=='':
            try:
                rawEarliestInfDate = raw_input('Enter the earliest possible date for the first infection in the epidemic (use dd/mm/yyyy format): ')
                earliestInfDate = datetime.strptime(rawEarliestInfDate, '%d/%m/%Y')
            except:
                print "Please enter a valid date (dd/mm/yy format)"
    dayOne = ''
    while dayOne=='' or (distribution=='b' and dayOne>=earliestInfDate):
        try:
            rawDayOne = raw_input('Enter a date to be used as day 1 of the model (before the earliest possible infection date if using a beta distribution for date of infection, otherwise arbitrary; best to pick one before any possible TMRCA, use dd/mm/yyyy format): ')
            dayOne = datetime.strptime(rawDayOne, '%d/%m/%Y')
        except:
            print "Please enter a valid date (dd/mm/yy format)"
    gammaScale=''
    while gammaScale == '':
        try:
            gammaScale = float(raw_input('Enter scale parameter of incubation period gamma distribution (incubation period in days): '))
        except:
            print "Please enter a valid decimal number"
    gammaShape=''
    while gammaShape== '':
        try:
            gammaShape = float(raw_input('Enter shape parameter of incubation period gamma distribution (incubation period in days): '))
        except:
            print "Please enter a valid decimal number"  
    generateXML(csvreader, averageIncubationPeriod, minInfectionToExam, distribution, standardDevInDays, earliestInfDate, dayOne, gammaScale, gammaShape, outputFileName, 1000, "network.log.txt") 

if __name__ == '__main__':
    main()
