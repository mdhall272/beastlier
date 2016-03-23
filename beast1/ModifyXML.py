from __future__ import division
import imp
import csv

try:
    imp.find_module('lxml')
    lxml = True
    import lxml.etree as ET
except ImportError:
    lxml = False
    import xml.etree.ElementTree as ET

from xml.dom import minidom
from datetime import datetime
import math
import argparse


dayOne = datetime.strptime("01/01/1900", '%d/%m/%Y')

coalescentLikelihoodNames = {'gmrfSkyrideLikelihood', 'coalescentLikelihood', 'generalizedSkyLineLikelihood',
                   'gmrfSkyGridLikelihood'}

coalescentModelNames = {'constantSize', 'exponentialGrowth', 'logisticGrowth', 'expansion'}

startingTrees = {'coalescentSimulator', 'coalescentTree', 'upgmaTree'}
kernelLookup = {'e': 'exponential', 'p': 'powerLaw', 'g': 'gaussian', 'l': "logistic", 'x': 'none'}

def modifyXML(epiFileName, taxaFileName, outputFileName, beautiFileName, fileNameRoot, startingNewick, startingTT,
              kernel, fixedPT, fixedTT, infType, infArgs, latType, latArgs, indexPrior, chainLength, sampleEvery,
              dateFormat, locationsInTable):

    epiReader = csv.reader(open(epiFileName, "rb"))

    parseAsDates = dateFormat != None

    latentPeriods = latType != 'x'

    infHyperpriorLog = infType == "l"

    outbreakString = 'categoryOutbreak'
    caseString = 'categoryCase'

    if lxml:
        parser = ET.XMLParser(remove_blank_text=True)

        fullTree = ET.parse(open(beautiFileName), parser)
    else :

        fullTree = ET.parse(open(beautiFileName))

    beastElement = fullTree.getroot()

    oldCoalescentParameters = list()

    # Remove the old tree prior
    for element in list(beastElement):
        if element.tag in coalescentLikelihoodNames:
            for parameterElement in element.iter("parameter"):
                oldCoalescentParameters.append(parameterElement.get("id"))
            beastElement.remove(element)
        # Non sky- models also generate the initial tree and should be kept for this purpose only
        if element.tag in coalescentModelNames:
            for parameterElement in element.iter("parameter"):
                oldCoalescentParameters.append(parameterElement.get("id"))

    # Remove the old taxa elements
    taxaElement = beastElement.find('taxa')

    for taxon in list(taxaElement.getchildren()):
        taxaElement.remove(taxon)

    alignmentElement = beastElement.find('alignment')

    outbreakElement = ET.Element(outbreakString)
    outbreakElement.set('id', 'outbreak')
    outbreakElement.set('hasLatentPeriods', str(latentPeriods).lower())
    outbreakElement.set('hasGeography', str(kernel!='x').lower())

    beastElement.insert(beastElement.getchildren().index(alignmentElement)+1, outbreakElement)

    # Infectious periods

    infCategoryName = 'inf1'

    infWrapperElement = ET.SubElement(outbreakElement, 'infectiousPeriodPrior')

    infHyperpriorType = 'normalPeriodPriorDistribution' if infType == "ng" else 'individualPrior'

    if infHyperpriorType=='normalPeriodPriorDistribution':
        infDistributionPriorElement = ET.SubElement(infWrapperElement, infHyperpriorType)
        infDistributionPriorElement.set("mu", infArgs[0])
        infDistributionPriorElement.set("lambda", infArgs[1])
        infDistributionPriorElement.set("alpha", infArgs[2])
        infDistributionPriorElement.set("beta", infArgs[3])
    else:
        infDistributionPriorElement = ET.SubElement(infWrapperElement, infHyperpriorType)

        if infType=="l" or infType=="n":
            distElement = createNormalBlock(ET, 'distribution', infArgs[0], math.sqrt(float(infArgs[1])))
        if infType=="e":
            distElement = createExponentialBlock(ET, 'distribution', float(infArgs[0]))
        if infType=="g":
            distElement = createGammaBlock(ET, 'distribution', float(infArgs[0]), float(infArgs[1]))

        infDistributionPriorElement.insert(0, distElement)



    infDistributionPriorElement.set("id", infCategoryName+'dist')
    infDistributionPriorElement.set("log", str(infHyperpriorLog).lower())

    # Latent periods

    if latentPeriods:
        latCategoryName = 'lat1' if latentPeriods else None

        latentPeriodStart = 0.0001

        latWrapperElement = ET.SubElement(outbreakElement, 'latentPeriods')
        latParameterElement = ET.SubElement(latWrapperElement, "parameter")
        latParameterElement.set("value", str(latentPeriodStart))
        latParameterElement.set("lower", "0")
        latParameterElement.set("id", latCategoryName+".latPeriod")

    # Taxa reference

    outbreakElement.append(createReferenceBlock(ET, 'taxa', 'taxa'))

    # Case data...

    # Gather all the branch position elements together (to be added as a child of the tree likelihood element)
    branchPositionsElement = ET.Element('compoundParameter')
    branchPositionsElement.set('id','infectionTimeBranchPositions')

    epiHeaderRow=epiReader.next()
    hostIDColumnEpi=None
    endTimeColumn=None
    longColumn=None
    latColumn=None

    for i in range(0, len(epiHeaderRow)):
        if epiHeaderRow[i] == 'Host_ID':
            hostIDColumnEpi = i
        elif epiHeaderRow[i] == 'End_date':
            endTimeColumn = i
        elif locationsInTable and epiHeaderRow[i] == 'Longitude' and kernel!='x':
            longColumn = i
        elif locationsInTable and epiHeaderRow[i] == 'Latitude' and kernel!='x':
            latColumn = i

    taxaReader = csv.reader(open(taxaFileName, "rb"))

    taxaHeaderRow=taxaReader.next()
    hostIDColumnTaxa=None
    taxaNameColumn=None

    for i in range(0, len(taxaHeaderRow)):
        if taxaHeaderRow[i] == 'Host_ID':
            hostIDColumnTaxa = i
        if taxaHeaderRow[i] == 'Taxon_ID':
            taxaNameColumn = i
        if taxaHeaderRow[i] == 'Exam_date':
            examTimeColumn = i

    currentRow = epiReader.next()
    branchPositionNames = list()

    caseCount = 0

    while currentRow is not None:

        taxaReader = csv.reader(open(taxaFileName, "rb"))
        taxaReader.next()

        hostName = currentRow[hostIDColumnEpi]


        latitude = None
        longitude = None
        if kernel != 'x':
            latitude = currentRow[latColumn]
            longitude = currentRow[longColumn]

        if currentRow[endTimeColumn] != "NA":

            if(parseAsDates):
                processedEndDate = datetime.strptime(currentRow[endTimeColumn], dateFormat)
                endTime = (processedEndDate - dayOne).days
            else:
                endTime = currentRow[endTimeColumn]

            taxa = list()

            currentTaxaRow = taxaReader.next()
            while currentTaxaRow is not None:

                if(currentTaxaRow[hostIDColumnTaxa] == hostName):

                    if(parseAsDates):
                        processedExamDate = datetime.strptime(currentTaxaRow[examTimeColumn], dateFormat)
                        examTime = (processedExamDate - dayOne).days
                    else:
                        examTime = currentTaxaRow[examTimeColumn]

                    if examTime > endTime:
                        raise Exception ("For host ID "+hostName+" at least one examination time is post end-of "
                                                                 "infection.")

                    taxon=currentTaxaRow[taxaNameColumn]

                    taxa.append(taxon)

                    taxonElement = ET.SubElement(taxaElement,'taxon')
                    taxonElement.set('id', currentTaxaRow[taxaNameColumn])
                    samplingDateElement = ET.SubElement(taxonElement,'date')
                    samplingDateElement.set('value', str(examTime))
                    samplingDateElement.set('direction','forwards')
                    samplingDateElement.set('units','days')
                    attributeElement = ET.SubElement(taxonElement,'attr')
                    attributeElement.set('name', 'hostID')
                    attributeElement.text = hostName

                try:
                    currentTaxaRow = taxaReader.next()
                except StopIteration:
                    currentTaxaRow = None


            latString = None

            if latentPeriods:
                latString = latCategoryName+".latPeriod"

            outbreakElement.append(createCategoryCaseElement(ET, caseString, hostName, endTime, longitude, latitude,
                                                             0.5, taxa, infCategoryName +'dist', latString))
            branchPositionsElement.append(createReferenceBlock(ET,'parameter', hostName+'_bp'))
            branchPositionNames.append(hostName+'_bp')

            caseCount = caseCount + 1

        else:

            outbreakElement.append(createNeverInfectedCategoryCaseElement(ET, caseString, hostName, longitude,
                                                                          latitude))
        try:
            currentRow = epiReader.next()
        except StopIteration:
            currentRow = None

    # Tree model

    treeModelElement = beastElement.find('treeModel')
    treeModelElement.tag = 'partitionedTreeModel'

    treeModelElement.append(createNestedReferenceBlock(ET, "outbreak", "parameter", "outbreak"))

    # replace the old starting tree if told to do so
    if startingNewick!=None:
        for element in list(beastElement):
            if element.tag in startingTrees:
                beastElement.remove(element)
        for element in list(treeModelElement):
            if element.tag in startingTrees:
                treeModelElement.insert(treeModelElement.getchildren().index(element),
                                        createReferenceBlock(ET, 'newick', 'startingTree'))
                treeModelElement.remove(element)
        newickElement = ET.Element('newick')
        newickElement.set('id', 'startingTree')
        newickElement.set('usingDates', 'true')
        newickFile = open(startingNewick)
        newickString = newickFile.read()
        newickElement.text = newickString

        beastElement.insert(beastElement.getchildren().index(treeModelElement), newickElement)

    operatorsElement = beastElement.find('operators')

    c2cTransElement = ET.Element('caseToCaseTransmissionLikelihood')
    c2cTransElement.set('id', 'c2cTransLikelihood')
    beastElement.insert(beastElement.getchildren().index(operatorsElement), c2cTransElement)
    beastElement.insert(beastElement.getchildren().index(c2cTransElement), ET.Comment(" Probability of the transmission "
                                                                                   "tree and epidemiological "
                                                                                   "parameters given the phylogenetic "
                                                                                   "tree"))

    # Transmission rate and geography

    if kernel != 'x':
        kernelElement = ET.SubElement(c2cTransElement, 'spatialKernelFunction')
        kernelElement.set('type', kernelLookup[kernel])
        parametersBlock = createParameterBlock(ET, 'parameters', 'kernel.alpha', 1, True, 1)

        kernelElement.append(parametersBlock)
        if kernel == 'g' or 'l':
            kernelElement.set('integratorSteps', str(20))
            if kernel == "l":
                parameterElement = ET.SubElement(parametersBlock,'parameter')
                parameterElement.set('value', str(1))
                parameterElement.set('id', "kernel.r0")
                parameterElement.set('lower', '0.0')



    trElement = createParameterBlock(ET, 'transmissionRate', 'transmission_rate', 0.05, True, 1)

    c2cTransElement.append(trElement)

    # Prior on the time of first infection

    if indexPrior is not None:
        indexPriorElement = ET.SubElement(c2cTransElement, 'initialInfectionTimePrior')
        normalDistElement = ET.SubElement(indexPriorElement, 'normalDistributionModel')
        meanElement = ET.SubElement(normalDistElement, 'mean')
        sdElement = ET.SubElement(normalDistElement, 'stdev')

        if(parseAsDates):
            processedMeanIndexTime = datetime.strptime(indexPrior[0], dateFormat)
            priorMeanIndexInfectionTime = (processedMeanIndexTime - dayOne).days
        else:
            priorMeanIndexInfectionTime = str(indexPrior[0])
        meanElement.text = priorMeanIndexInfectionTime
        sdElement.text = str(indexPrior[1])

    # Tree likelihood

    c2cTreeElement = ET.SubElement(c2cTransElement, 'withinCaseCoalescent')

    c2cTreeElement.set("truncate", "true")

    if startingTT is not None:
        startingTTElement = ET.SubElement(c2cTreeElement, 'startingNetwork')
        startingTTElement.text = startingTT

    c2cTreeElement.append(createReferenceBlock(ET, outbreakString, 'outbreak'))
    c2cTreeElement.append(createReferenceBlock(ET,'treeModel','treeModel'))
    c2cTreeElement.set('id', 'withinCaseCoalescent')
    c2cTreeElement.append(createParameterBlock(ET, 'maxFirstInfToRoot', 'max_first_inf_to_root', 30, True, 1))

    infTimesWrapper = ET.Element('infectionTimeBranchPositions')
    c2cTreeElement.append(infTimesWrapper)
    infTimesWrapper.append(branchPositionsElement)

    # Within-host demographic model (fixed to logistic for this script)

    demoModelWrapper = ET.Element('demographicModel')

    logisticElement = ET.SubElement(demoModelWrapper, 'logisticGrowthN0')
    logisticElement.set('units', 'days')
    logisticElement.append(createParameterBlock(ET, 'populationSize', 'logistic.startingNe', 0.1, True, 1))
    logisticElement.append(createParameterBlock(ET, 'growthRate', 'logistic.growthRate', 1, False, 1))
    logisticElement.append(createParameterBlock(ET, 't50', 'logistic.t50', -4, False, 1))
    c2cTreeElement.append(demoModelWrapper)
    minusStatistic = ET.SubElement(c2cTreeElement, 'negativeStatistic')
    minusStatistic.set('id', 'minusT50')
    minusStatistic.append(createReferenceBlock(ET, 'parameter', 'logistic.t50'))

    ratioStatistic = ET.SubElement(c2cTreeElement, 'sumStatistic')
    ratioStatistic.set('id', 'logistic.ratio')
    oneElement = ET.SubElement(ratioStatistic, 'parameter')
    oneElement.set('value', '1')
    expElement = ET.SubElement(ratioStatistic, 'exponentialStatistic')
    productElement = ET.SubElement(expElement, 'productStatistic')
    productElement.append(createReferenceBlock(ET, 'parameter', 'logistic.growthRate'))
    productElement.append(createReferenceBlock(ET, 'negativeStatistic', 'minusT50'))

    # Operators

    operatorsBlock = beastElement.find('operators')

    # First get rid of the old coalescent model

    for operator in list(operatorsBlock.getchildren()):
        if operator.tag == "gmrfBlockUpdateOperator":
            operatorsBlock.remove(operator)
        if operator[0].get("idref") is not None and operator[0].get("idref") in oldCoalescentParameters:
            operatorsBlock.remove(operator)

    if not fixedPT:
        treeOperatorsToRemove = {'narrowExchange', 'wideExchange', 'wilsonBalding', 'subtreeSlide'}
        if fixedTT:
            treeOperatorsToAdd = {'transmissionExchangeOperatorA', 'transmissionWilsonBaldingA',
                                  'transmissionSubtreeSlideA'}
        else:
            treeOperatorsToAdd = {'transmissionExchangeOperatorA', 'transmissionExchangeOperatorB',
                              'transmissionWilsonBaldingA', 'transmissionWilsonBaldingB', 'transmissionSubtreeSlideA',
                              'transmissionSubtreeSlideB'}

        for operator in list(operatorsBlock.getchildren()):
            if operator.tag in treeOperatorsToRemove:
                operatorsBlock.remove(operator)

        for operatorName in treeOperatorsToAdd:
            newOperatorElement = ET.SubElement(operatorsBlock, operatorName)
            newOperatorElement.set('weight', str(8))
            newOperatorElement.set('resampleInfectionTimes', 'true')
            newOperatorElement.append(createReferenceBlock(ET, "withinCaseCoalescent", "withinCaseCoalescent"))
            if operatorName.startswith('transmissionSubtreeSlide'):
                newOperatorElement.set('gaussian', 'true')
    else:
        # anything that's left at this point has to do with the tree
        for operator in list(operatorsBlock):
            operatorsBlock.remove(operator)

    if not fixedTT:
        ibmElement = ET.SubElement(operatorsBlock,'infectionBranchMovementOperator')
        ibmElement.set('weight', str(15))
        ibmElement.set('resampleInfectionTimes', 'true')
        ibmElement.append(createReferenceBlock(ET, "withinCaseCoalescent", "withinCaseCoalescent"))

    bpUniformOpElement = ET.SubElement(operatorsBlock, 'uniformOperator')
    bpUniformOpElement.set('weight', str(5))
    bpUniformOpElement.set('lower', str(0))
    bpUniformOpElement.set('upper', str(1))
    bpUniformOpElement.append(createReferenceBlock(ET, 'parameter', 'infectionTimeBranchPositions'))

    transRateScaleElement = ET.SubElement(operatorsBlock, 'scaleOperator')
    transRateScaleElement.set('weight', str(5))
    transRateScaleElement.set('scaleFactor', str(0.75))
    transRateScaleElement.append(createReferenceBlock(ET, 'parameter', 'transmission_rate'))

    if kernel != 'x':
        kernelaScaleElement = ET.SubElement(operatorsBlock, 'scaleOperator')
        kernelaScaleElement.set('weight', str(5))
        kernelaScaleElement.set('scaleFactor', str(0.75))
        kernelaScaleElement.append(createReferenceBlock(ET, 'parameter', 'kernel.alpha'))

        if kernel == 'l':
            kernelaScaleElement = ET.SubElement(operatorsBlock, 'scaleOperator')
            kernelaScaleElement.set('weight', str(5))
            kernelaScaleElement.set('scaleFactor', str(0.75))
            kernelaScaleElement.append(createReferenceBlock(ET, 'parameter', 'kernel.r0'))

    if latentPeriods:
        latentPeriodScaleElement = ET.SubElement(operatorsBlock, 'scaleOperator')
        latentPeriodScaleElement.set('weight', str(5))
        latentPeriodScaleElement.set('scaleFactor', str(0.75))
        latentPeriodScaleElement.append(createReferenceBlock(ET, 'parameter', latCategoryName+".latPeriod"))

    growthRateScaleElement = ET.SubElement(operatorsBlock, 'scaleOperator')
    growthRateScaleElement.set('weight', str(5))
    growthRateScaleElement.set('scaleFactor', str(0.75))
    growthRateScaleElement.append(createReferenceBlock(ET, 'parameter', 'logistic.growthRate'))

    startingNeScaleElement = ET.SubElement(operatorsBlock, 'scaleOperator')
    startingNeScaleElement.set('weight', str(5))
    startingNeScaleElement.set('scaleFactor', str(0.75))
    startingNeScaleElement.append(createReferenceBlock(ET, 'parameter', 'logistic.startingNe'))

    ttScaleElement = ET.SubElement(operatorsBlock, 'scaleOperator')
    ttScaleElement.set('weight', str(5))
    ttScaleElement.set('scaleFactor', str(0.75))
    ttScaleElement.append(createReferenceBlock(ET, 'parameter', 'logistic.t50'))

    mcmcBlock = beastElement.find('mcmc')
    mcmcBlock.set('operatorAnalysis', fileNameRoot+".ops.txt")
    if(chainLength!=None):
        mcmcBlock.set('chainLength', str(chainLength))
    posteriorBlock = mcmcBlock.find('posterior')

    # Priors

    priorBlock = posteriorBlock.find('prior')
    for element in list(priorBlock):
        if element.tag in coalescentLikelihoodNames:
            priorBlock.remove(element)
        if len(element)>0 \
                and element[0].get('idref') is not None \
                and element[0].get('idref') in oldCoalescentParameters:
            priorBlock.remove(element)

    minusTTPriorBlock = ET.SubElement(priorBlock, 'gammaPrior')
    minusTTPriorBlock.set('shape', '10')
    minusTTPriorBlock.set('scale', '0.5')
    minusTTPriorBlock.append(createReferenceBlock(ET, 'negativeStatistic', 'minusT50'))

    ratioPriorBlock = ET.SubElement(priorBlock, 'logNormalPrior')
    ratioPriorBlock.set('mean', '4')
    ratioPriorBlock.set('stdev', '0.5')
    ratioPriorBlock.set('offset', '0')
    ratioPriorBlock.set('meanInRealSpace', 'false')
    ratioPriorBlock.append(createReferenceBlock(ET, 'parameter', 'logistic.ratio'))

    if kernel != "x":
        kaPriorBlock = ET.SubElement(priorBlock, 'exponentialPrior')
        kaPriorBlock.set('mean', '1')
        kaPriorBlock.append(createReferenceBlock(ET, 'parameter', 'kernel.alpha'))

        if kernel == "l":
            kr0PriorBlock = ET.SubElement(priorBlock, 'exponentialPrior')
            kr0PriorBlock.set('mean', '1')
            kr0PriorBlock.append(createReferenceBlock(ET, 'parameter', 'kernel.r0'))

    if latentPeriods:
        if latType=="n":
            latDistributionPriorElement = createNestedReferenceBlock(ET, "normalPrior", "parameter",
                                                                     latCategoryName+".latPeriod")
            latDistributionPriorElement.set("mean", latArgs[0])
            latDistributionPriorElement.set("stdev", str(math.sqrt(float(latArgs[1]))))
        if latType=="l":
            latDistributionPriorElement = createNestedReferenceBlock(ET, "logNormalPrior", "parameter",
                                                                     latCategoryName+".latPeriod")
            latDistributionPriorElement.set("mean", latArgs[0])
            latDistributionPriorElement.set("stdev", str(math.sqrt(float(latArgs[1]))))
        if latType=="e":
            latDistributionPriorElement = createNestedReferenceBlock(ET, "exponentialPrior", "parameter",
                                                                     latCategoryName+".latPeriod")
            latDistributionPriorElement.set("rate", latArgs[0])
        if latType=="g":
            latDistributionPriorElement = createNestedReferenceBlock(ET, "gammaPrior", "parameter",
                                                                     latCategoryName+".latPeriod")
            latDistributionPriorElement.set("shape", latArgs[0])
            latDistributionPriorElement.set("scale", latArgs[1])

    priorBlock.append(latDistributionPriorElement)

    priorBlock.append(createReferenceBlock(ET, 'caseToCaseTransmissionLikelihood', 'c2cTransLikelihood'))

    if fixedPT:
        likelihoodBlock = posteriorBlock.find('prior')
        for irrelevant in list(likelihoodBlock):
            likelihoodBlock.remove(irrelevant)


    # Logs

    for block in mcmcBlock:
        if block.tag == "log":
            if block.get('id')=='screenLog':
                if sampleEvery is not None:
                    block.set('logEvery', str(sampleEvery))

                block.append(createScreenLogColumnBlock(ET, 'parameter', 'trans_rate', 'transmission_rate', 6, 12))


                if kernel != 'x':
                    block.append(createScreenLogColumnBlock(ET, 'parameter', 'kernel.alpha', 'kernel.alpha', 4, 12))
                    if kernel == "l":
                        block.append(createScreenLogColumnBlock(ET, 'parameter', 'kernel.r0', 'kernel.r0', 4, 12))

                if latentPeriods:
                    block.append(createScreenLogColumnBlock(ET, 'parameter', 'lat.period', 'lat1.latPeriod', 4, 12))

                block.append(createScreenLogColumnBlock(ET, 'parameter', 'lg.ratio', 'logistic.ratio', 4, 12))
                block.append(createScreenLogColumnBlock(ET, 'productStatistic', 'lg.t50', 'logistic.t50', 4, 12))
                block.append(createScreenLogColumnBlock(ET, 'parameter', 'lg.gr', 'logistic.growthRate', 4, 12))

            elif block.get('id')=='fileLog':

                for child in list(block):
                    if child.get('idref') in oldCoalescentParameters \
                            or child.tag in coalescentLikelihoodNames:
                        block.remove(child)

                block.set('fileName',fileNameRoot+".log.txt")

                if sampleEvery is not None:
                    block.set('logEvery', str(sampleEvery))
                else:
                    oldSampleEvery = block.get("logEvery")


                if kernel != 'x':
                    block.append(createReferenceBlock(ET, 'parameter', 'kernel.alpha'))
                    if kernel == "l":
                        block.append(createReferenceBlock(ET, 'parameter', 'kernel.r0'))

                block.append(createReferenceBlock(ET, 'parameter', 'transmission_rate'))
                block.append(createReferenceBlock(ET, 'caseToCaseTransmissionLikelihood', 'c2cTransLikelihood'))

                if latentPeriods:
                    block.append(createReferenceBlock(ET, 'parameter', 'lat1.latPeriod'))

                block.append(createReferenceBlock(ET, 'parameter', 'logistic.ratio'))
                block.append(createReferenceBlock(ET, 'productStatistic', 'logistic.t50'))
                block.append(createReferenceBlock(ET, 'parameter', 'logistic.growthRate'))




    networkLogBlock = ET.Element('log')
    networkLogBlock.set('id', 'networkLog')
    if sampleEvery is not None:
        networkLogBlock.set('logEvery', str(sampleEvery))
    else:
        networkLogBlock.set('logEvery', str(oldSampleEvery))

    networkLogBlock.set('fileName', fileNameRoot+".net.txt")
    networkLogBlock.set('overwrite', 'true')

    networkLogBlock.append(createReferenceBlock(ET, "withinCaseCoalescent", "withinCaseCoalescent"))
    mcmcBlock.insert(mcmcBlock.getchildren().index(posteriorBlock)+4, networkLogBlock)

    treeLogBlock = mcmcBlock.find('logTree')

    treeLogBlock.set('fileName', fileNameRoot+".trees.txt")
    if sampleEvery is not None:
        treeLogBlock.set('logEvery', str(sampleEvery))
    else:
        treeLogBlock.set('logEvery', str(oldSampleEvery))
    firstTrait = treeLogBlock.find('trait')
    partitionTrait = ET.Element('trait')
    partitionTrait.set('name', 'partition')
    partitionTrait.set('tag', 'partition')
    partitionTrait.append(createReferenceBlock(ET,  "withinCaseCoalescent",  "withinCaseCoalescent"))
    treeLogBlock.insert(treeLogBlock.getchildren().index(firstTrait), partitionTrait)

    if lxml:
        outString = ET.tostring(beastElement, pretty_print = True)

        outFile = open(outputFileName, "w+")
        outFile.write(outString)

    else:
        outString = ET.tostring(beastElement)

        outFile = open(outputFileName, "w+")
        outFile.write(minidom.parseString(outString).toprettyxml(newl='\n'))


def createExponentialBlock(tree, name, mean, idstring,  meanID):
    enclosingElement = tree.Element(name)
    expElement = tree.SubElement(enclosingElement, 'exponentialDistributionModel')
    expElement.append(createParameterBlock(tree, 'mean', meanID, mean, True, 1))
    if idstring is not None:
        expElement.set('id', idstring)
    return enclosingElement

def createGammaBlock(tree, name, scale, shape, idstring, scaleID, shapeID):
    enclosingElement = tree.Element(name)
    gammaElement = tree.SubElement(enclosingElement,'gammaDistributionModel')
    gammaElement.append(createParameterBlock(tree, 'shape', shapeID, shape, True, 1))
    gammaElement.append(createParameterBlock(tree, 'scale', scaleID, scale, True, 1))
    if idstring is not None:
        gammaElement.set('id', idstring)
    return enclosingElement

def createNormalBlock(tree, name, mean, stdev, idstring = None, scaleID = None, shapeID = None):
    enclosingElement = tree.Element(name)
    gammaElement = tree.SubElement(enclosingElement,'normalDistributionModel')
    gammaElement.append(createParameterBlock(tree, 'mean', shapeID, mean, True, 1))
    gammaElement.append(createParameterBlock(tree, 'stdev', scaleID, stdev, True, 1))
    if idstring is not None:
        gammaElement.set('id', idstring)
    return enclosingElement

def createParameterBlock(tree, name, idstring, value, greaterThanZero, dim):
    enclosingElement = tree.Element(name)
    parameterElement = tree.SubElement(enclosingElement,'parameter')
    parameterElement.set('value', str(value))
    if idstring is not None:
        parameterElement.set('id', idstring)
    if greaterThanZero:
        parameterElement.set('lower', '0.0')
    if dim>1:
        parameterElement.set('dimension', str(dim))
    return enclosingElement

def createParameterBlockWithMax(tree, name, idstring, value, maxVal, greaterThanZero, dim):
    enclosingElement = tree.Element(name)
    parameterElement = tree.SubElement(enclosingElement,'parameter')
    parameterElement.set('value', str(value))
    parameterElement.set('id', idstring)
    parameterElement.set('upper', str(maxVal))
    if greaterThanZero:
        parameterElement.set('lower', '0.0')
    if dim>1:
        parameterElement.set('dimension', str(dim))
    return enclosingElement

def createReferenceBlock(tree, name, reference):
    refBlock = tree.Element(name)
    refBlock.set('idref', reference)
    return refBlock

def createNestedReferenceBlock(tree, parentName, childName, reference):
    parentBlock = tree.Element(parentName)
    parentBlock.append(createReferenceBlock(tree, childName, reference))
    return parentBlock

def createCategoryCaseElement(tree, name, caseID, endTime, longitude, latitude, bp, taxa, infCategoryName,
                              latCategoryName):
    caseElement = tree.Element(name)
    caseElement.set('wasEverInfected', 'true')
    caseElement.set('caseID', caseID)
    caseElement.set('id', caseID)
    caseElement.set('infectiousCategory', infCategoryName)
    if latCategoryName!=None:
        caseElement.set('latentCategory', latCategoryName)
    caseElement.set('endTime', str(endTime))
    caseElement.append(createParameterBlockWithMax(tree, 'infectionTimeBranchPosition', caseID+"_bp", bp, 1, True, 1))

    if longitude!=None and latitude!=None:
        coordsString = str(longitude)+" "+str(latitude)
        caseElement.append(createParameterBlock(tree, 'spatialCoordinates', None, coordsString, False, 2))

    for taxon in taxa:
        caseElement.append(createReferenceBlock(tree, 'taxon', taxon))
    return caseElement

def createNeverInfectedCategoryCaseElement(tree, name, caseID, longitude, latitude):
    caseElement = tree.Element(name)

    caseElement.set('wasEverInfected', 'false')

    caseElement.set('caseID', caseID)
    caseElement.set('id', caseID)

    if longitude!=None and latitude!=None:
        coordsString = str(longitude)+" "+str(latitude)
        caseElement.append(createParameterBlock(tree, 'spatialCoordinates', None, coordsString, False, 2))

    return caseElement

def createScreenLogColumnBlock(tree, name, label, reference, sf, width):
    columnElement = tree.Element('column')
    columnElement.set('label', label)
    columnElement.set('sf', str(sf))
    columnElement.set('width', str(width))
    innerElement = tree.SubElement(columnElement,name)
    innerElement.set('idref', reference)
    return columnElement


# Main method; command line entry for data that isn't in the CSV file or provided as arguments

def main():

    validKernels = {'e', 'p', 'g', 'l', 'x'}
    validInfHyperpriors = {'ng', 'n', 'g', 'l', 'e'}
    infArgumentCounts = {'ng' : 4, 'n' : 2, 'g' : 2, 'l' : 2, 'e' : 1}
    infArgumentStrings = {'ng' : ["mu", "lambda", "alpha", "beta"],
                    'n' : ["mean", "variance"],
                    'l' : ["mean", "variance"],
                    'g' : ["shape", "scale"],
                    'e' : ["rate"]
                    }
    validLatPriors = {'n', 'g', 'l', 'e', 'x'}

    latArgumentCounts = {'n' : 2, 'g' : 2, 'l' : 2, 'e' : 1, 'x' : 0}
    latArgumentStrings = {'n' : ["mean", "variance"],
                    'l' : ["mean", "variance"],
                    'g' : ["shape", "scale"],
                    'e' : ["rate"]
                    }



    parser = argparse.ArgumentParser(description='Modifies a BEAST 1 input file (e.g. generated by BEAUTi) for '
                                                 'transmission tree reconstruction (as per Hall et al., PLOS Comp'
                                                 ' Biol, 2016).',
                                     prog="ModifyXML")
    parser.add_argument('beautiFile', help='The path of an existing XML file (e.g. generated by BEAUTi), ready to '
                                           'have the elements for transmission tree reconstruction added.')
    parser.add_argument('epiFile', help='The path of a CSV file of epidemiological information; the file must have '
                                        'columns "Host_ID" and "End_date" and may have "longitude" and'
                                        ' "latitude". This file can contain details of known never-infected hosts'
                                        ' which should be given an end date of NA.')
    parser.add_argument('taxaFile', help='The path of a CSV file of matching each taxon to its host. Must have columns'
                                         ' "Taxon_ID", "Exam_date" and "Host_ID".')
    parser.add_argument('outputXMLFile', help='The name of the output XML file.')
    parser.add_argument('fileNameRoot', help='The root of the file names that BEAST will write the output to')
    parser.add_argument('-d', '--dateFormat', help='The date format (see documentation for the datetime Python library '
                                                   'for formats) for the date entries in the data table. If absent, '
                                                   'times of sampling and noninfectiousness are parsed as floating '
                                                   'point numbers.')
    parser.add_argument('-k', '--kernel', default="n", help='The type of spatial kernel to use (e=exponential, p=power'
                                                            ' law, g=Gaussian, l=logistic, x=none (no geography))')
    parser.add_argument('-i', '--infectiousPeriods', help='The prior distribution on the length of infectious periods.'
                                                          ' Accepts variable arguments; the first is the type of'
                                                          ' distribution and the following are its parameters.'
                                                          ' Choices: Normal-Gamma ("ng", four arguments: mu, lambda,'
                                                          ' alpha, beta), Normal ("n", two arguments: mean, variance),'
                                                          ' Lognormal ("l", two arguments, mean, variance (on log '
                                                          'scale)), Gamma ("g", two arguments, shape, scale),'
                                                          ' exponential ("e", one argument, rate).',
                        nargs='+')
    parser.add_argument('-l', '--latentPeriods', help='The prior distribution on the length of latent periods. See -i '
                                                      'for explanation. Choices: Normal ("n", two arguments: mean, '
                                                      'variance), Lognormal ("l", two arguments, mean, variance (on '
                                                      'log scale)), Gamma ("g", two arguments, shape, scale), '
                                                      'Exponential ("e", one argument, rate), no latent periods ("x", '
                                                      'no further arguments).',
                        nargs='+')
    parser.add_argument('-s', '--startingPTree', help='If a Newick filename is given here, use this as the starting '
                                                     'phylogenetic tree.')
    parser.add_argument('-t', '--startingTTree', help='If a CSV filename is given here, use this as the starting'
                                                      ' transmission tree (WARNING: this may not be compatible'
                                                      ' with the starting phylogeny).')
    parser.add_argument('-f', '--fixedPT', default=False, help='Run on a fixed phylogenetic tree (requires specified '
                                                               '-startingPTree).')
    parser.add_argument('-g', '--fixedTT', default=False, help='Run on a fixed transmission tree (requires specified '
                                                               '-startingTTree)')
    parser.add_argument('-c', '--chainLength', help="Length of the MCMC chain (if unspecified, will keep what's in the"
                                                    " original XML file).")
    parser.add_argument('-e', '--sampleEvery', help="Sampling frequency (if unspecified, will keep what's in the"
                                                    " original XML file).")
    parser.add_argument('-x', '--indexDatePrior', help="Parameters of the Normal prior distribution for the date of the"
                                                       "index infection. Two arguments; mean and stdev. The mean is"
                                                       " in the format determined by -d and the stdev is in days. If"
                                                       " absent, an uninformative prior will be used.", nargs=2)

    arguments = parser.parse_args()

    epiFileName = arguments.epiFile

    epiReader = csv.reader(open(epiFileName, 'rU'))

    headers = epiReader.next()

    for necessaryHeader in {"Host_ID", "End_date"}:
        if necessaryHeader not in headers:
            raise Exception("File "+epiFileName+" does not have required column named "+necessaryHeader)

    hasLongitude = "Longitude" in headers
    hasLatitude = "Latitude" in headers

    if hasLatitude != hasLongitude:
        raise Exception("File "+epiFileName+" has a column for only one spatial dimension")
    else:
        locationsInTable = hasLongitude

    taxaFileName = arguments.taxaFile

    taxaReader = csv.reader(open(taxaFileName, 'rU'))

    headers = taxaReader.next()

    for necessaryHeader in {"Taxon_ID", "Host_ID", "Exam_date"}:
        if necessaryHeader not in headers:
            raise Exception("File "+taxaFileName+" does not have required column named "+necessaryHeader)

    # todo write this file to a dict to send to the XML writer

    outputFileName = arguments.outputXMLFile
    fileNameRoot = arguments.fileNameRoot
    beautiFileName = arguments.beautiFile
    startingNewick = arguments.startingPTree

    if startingNewick is None:
        print "Using a random starting phylogeny"

    startingTT = arguments.startingTTree

    if startingTT is None:
        print "Using a random starting transmission tree"

    chainLength = arguments.chainLength

    fixPT = arguments.fixedPT
    if fixPT:
        if startingNewick is None:
            raise Exception("You asked for a fixed phylogeny but did not specify it.")
        print "Fixing the phylogeny"

    fixTT = arguments.fixedTT
    if fixTT:
        if startingTT is None:
            raise Exception("You asked for a fixed transmission tree but did not specify it.")
        print "Fixing the transmission tree"

    if chainLength is None:
        print "Keeping the chain length"

    sampleEvery = arguments.sampleEvery

    if sampleEvery is None:
        print "Keeping the sampling frequency"

    dateFormat = arguments.dateFormat

    if dateFormat is None:
        print "Parsing dates as numbers"

    indexPrior = arguments.indexDatePrior

    if indexPrior is None:
        print "Using a noninformative prior on the date of the index infection"

    if arguments.kernel != 'x':
        if arguments.kernel not in validKernels:
            raise Exception("Invalid spatial transmission kernel")
        if not hasLongitude:
            raise Exception("Spatial transmission kernel specified but host data lists no coordinates")
        else:
            kernelType = arguments.kernel
    else:
        kernelType = "x"
        print "Assuming no geography"

    if arguments.infectiousPeriods is not None:
        if arguments.infectiousPeriods[0] not in validInfHyperpriors:
            raise Exception("Invalid infectious period prior distribution")
        elif len(arguments.infectiousPeriods) != (infArgumentCounts.get(arguments.infectiousPeriods[0]) + 1):
            raise Exception("Wrong number of arguments for infectious period prior distribution")
        else:
            infType = arguments.infectiousPeriods[0]
            infArguments = arguments.infectiousPeriods[1:len(arguments.infectiousPeriods)]
    else:
        infType = ""
        while infType not in validInfHyperpriors:
            infType = raw_input("Please specify type of prior distribution for infectious periods (ng=Normal-Gamma,"
                                " n=normal, l=lognormal, g=gamma, e=exponential): ")
        infArguments = list()
        if infType != "x":
            for i in range(0,infArgumentCounts.get(infType)):
                input = ""
                while not input.replace('.','',1).isdigit():
                    input = raw_input("Enter "+infArgumentStrings.get(infType)[i]+" parameter of prior distribution: ")
                infArguments.append(float(input))

    if arguments.latentPeriods is not None:
        if arguments.latentPeriods[0] not in validLatPriors:
            raise Exception("Invalid latent period prior distribution")
        elif len(arguments.latentPeriods) != (latArgumentCounts.get(arguments.latentPeriods[0]) + 1):
            raise Exception("Wrong number of arguments for latent period prior distribution")
        else:
            latType = arguments.latentPeriods[0]
            latArguments = arguments.latentPeriods[1:len(arguments.latentPeriods)]
    else:
        latType = ""
        while latType not in validLatPriors:
            latType = raw_input("Please specify type of prior distribution for latent periods ("
                                "n=normal, l=lognormal, g=gamma, e=exponential, x=no latent periods): ")
        latArguments = list()
        if latType != "x":
            for i in range(0,latArgumentCounts.get(latType)):
                input = ""
                while not input.replace('.','',1).isdigit():
                    input = raw_input("Enter "+latArgumentStrings.get(latType)[i]+" parameter of prior distribution: ")
                latArguments.append(float(input))

    modifyXML(epiFileName, taxaFileName, outputFileName, beautiFileName, fileNameRoot, startingNewick, startingTT,
              kernelType, fixPT, fixTT, infType, infArguments, latType, latArguments, indexPrior, chainLength,
              sampleEvery, dateFormat, locationsInTable)


if __name__ == '__main__':
    main()
