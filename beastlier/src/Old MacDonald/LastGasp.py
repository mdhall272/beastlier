from lxml import etree

def run(inFileName, outFileName, fileNameRoot):

    parser = etree.XMLParser(remove_blank_text=True)

    fullTree = etree.parse(open(inFileName), parser)

    beastElement = fullTree.getroot()

    mcmcElement = beastElement.find("mcmc")
    mcmcElement.set("operatorAnalysis", fileNameRoot+".ops.txt")
    mcmcElement.set("chainLength", "250000000")


    priorElement = mcmcElement.find("posterior").find("prior")

    kernelParameterPriorElement = etree.SubElement(priorElement, "exponentialPrior")
    kernelParameterPriorElement.set("mean", "10")

    referenceElement = etree.SubElement(kernelParameterPriorElement, "parameter")
    referenceElement.set("idref", "kernel.alpha")

    wccBlock = beastElement.find("caseToCaseTransmissionLikelihood").find("withinCaseCoalescent")
    wccBlock.set("truncate", "true")



    for logBlock in beastElement.iter('log'):
        logBlock.set('logEvery', "25000")

        if logBlock.get('id')=='fileLog':
            logBlock.set('fileName',fileNameRoot+".log.txt")
        elif logBlock.get('id')=='networkLog':
            logBlock.set('fileName',fileNameRoot+".net.txt")

    for treeLogBlock in beastElement.iter('logTree'):
        treeLogBlock.set('logEvery', '25000')

        treeLogBlock.set('fileName', fileNameRoot+".trees.txt")




    et = etree.ElementTree(beastElement)

    et.write(outFileName, pretty_print=True)






def main():
    for i in range(1,26):
        inFileName = "prior/wcc_X_" + str(i) + "_prior.xml"
        outFileName = "prior_3/wcc_X_" + str(i) + "_prior_3.xml"
        fileNameRoot = "wcc_X_" + str(i) + "_prior_3"
        run(inFileName, outFileName, fileNameRoot)



if __name__ == '__main__':
    main()


