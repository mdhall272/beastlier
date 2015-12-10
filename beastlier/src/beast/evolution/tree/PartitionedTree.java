
package beast.evolution.tree;

import beast.app.beauti.BeautiDoc;
import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.parameter.IntegerParameter;
import beast.util.Randomizer;
import beast.util.TreeParser;

import java.util.*;

import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Matthew Hall
 * @author Nicola De Maio
 */
@Description("A tree with node partitions as defined by Hall et al")
public class PartitionedTree extends Tree {

    /*
     * Plugin inputs:
     */
//    public Input<Integer> nTypesInput = new Input<Integer>(
//            "nTypes",
//            "Number of distinct types to consider.", Input.Validate.REQUIRED);

    public Input<String> elementLabelInput = new Input<String>(
            "elementLabel",
            "Label for type traits (default 'element')", "element");
    public Input<Boolean> allowCreepInput = new Input<Boolean>(
            "allowCreep", "Allow nodes that aren't ancestral under partition?");
    public Input<IntegerParameter> leafElementsInput = new Input<IntegerParameter>(
            "leafElements",
            "Elements containing leaf nodes.");


    /*
     * Non-input fields:
     */
    protected String elementLabel;

    protected TraitSet elementTraitSet;

    //allowCreep = true is within-host diversity. = false means lineages split at transmission if there is only one
    //sample per host. It's incoherent if not, and sampled ancestors would be needed to resolve this

    protected Boolean allowCreep;

    protected List<String> elementList;

    private List<Integer> leafElements;
    private List<String> leafNames;

    public PartitionedTree() { };

    public PartitionedTree(PartitionedTreeNode rootNode) {
        //public MultiTypeTreeVolz(Node rootNode) {

        if (!(rootNode instanceof PartitionedTreeNode)) {
            throw new IllegalArgumentException("Attempted to instantiate partitioned tree with regular root node.");
        }

        setRoot(rootNode);
        initArrays();
    }

    @Override
    public void initAndValidate() throws Exception {

        super.initAndValidate();

        leafElements = new ArrayList<>();
        leafNames = new ArrayList<>();
        if (leafElementsInput.get() != null) {
            for (int i=0; i<leafElementsInput.get().getDimension(); i++) {
                leafElements.add(leafElementsInput.get().getValue(i));
                leafNames.add(String.valueOf(i));
            }
        } else {
            if (!hasElementTrait())
                throw new IllegalArgumentException("Either leafColours or element set (with name '" + elementLabel
                        + "') must be provided.");

            // Fill leaf colour array:
            for (String taxonName : elementTraitSet.taxaInput.get().asStringList()) {
                leafElements.add(getElementList().indexOf(elementTraitSet.getStringValue(taxonName)));
                leafNames.add(taxonName);
            }
        }

        for (Node node : getExternalNodes()) {
            ((PartitionedTreeNode) node).setPartitionElementNumber(
                    getElementList().indexOf(elementTraitSet.getStringValue(getTaxonId(node))));
        }


        // all partition elements containing more than one tip must be assigned first

        for(int elementNo = 0; elementNo < elementList.size(); elementNo++){
            if(getTipsInElement(elementNo).size()>1){
                PartitionedTreeNode mrca = getElementMRCA(elementNo);


                for(PartitionedTreeNode node : getTipsInElement(elementNo)){
                    PartitionedTreeNode currentNode = node;
                    while(currentNode!=mrca){
                        currentNode = (PartitionedTreeNode)node.getParent();
                        if(currentNode.getPartitionElementNumber()==-1){
                            break;
                        }
                        currentNode.setPartitionElementNumber(elementNo);
                    }

                }

                //if this procedure has failed to reach the MRCA, then the tree can't be partitioned

                if(mrca.getPartitionElementNumber()!=elementNo){
                    throw new Exception("The starting tree is not compatible with this leaf partition");
                }

            }

        }

        partitionDown((PartitionedTreeNode)getRoot());

        if(!checkPartitionIntegrity()){
            throw new Exception("Partition integrity check failed on initial tree");
        }


    }



    protected void partitionDown(PartitionedTreeNode node){
        for(Node child : node.getChildren()){
            partitionDown((PartitionedTreeNode)child);
        }

        if(node.getPartitionElementNumber()==-1){
            ArrayList<Integer> childPartitions = new ArrayList<>();

            for(Node child : node.getChildren()){

                childPartitions.add(((PartitionedTreeNode)child).getPartitionElementNumber());

            }


            Integer choice = childPartitions.get(Randomizer.nextInt(2));
            node.setPartitionElementNumber(choice);

        }

    }




    @Override
    protected void processTraits(List<TraitSet> traitList) {
        super.processTraits(traitList);

        // Record trait set associated with leaf types.
        for (TraitSet traitSet : traitList) {
            if (traitSet.getTraitName().equals(elementLabel)) {
                elementTraitSet = traitSet;
                break;
            }
        }

        // Construct type list.
        if (elementTraitSet == null) {
            if (getTaxonset() != null) {
                TraitSet dummyTraitSet = new TraitSet();

                StringBuilder sb = new StringBuilder();
                for (int i=0; i<getTaxonset().getTaxonCount(); i++) {
                    if (i>0)
                        sb.append(",\n");
                    sb.append(getTaxonset().getTaxonId(i)).append("=NOT_SET");
                }
                try {
                    dummyTraitSet.initByName(
                            "traitname", "type",
                            "taxa", getTaxonset(),
                            "value", sb.toString());
                    dummyTraitSet.setID("elementTraitSet.t:"
                            + BeautiDoc.parsePartition(getID()));
                    setElementTrait(dummyTraitSet);
                } catch (Exception ex) {
                    System.out.println("Error setting default type trait.");
                }
            }
        } else {

            Set<String> typeSet = new HashSet<>();

            int nTaxa = elementTraitSet.taxaInput.get().asStringList().size();
            for (int i = 0; i < nTaxa; i++)
                typeSet.add(elementTraitSet.getStringValue(i));

            elementList = new ArrayList<>(typeSet);
            Collections.sort(elementList);

            System.out.println("Type trait with the following types detected:");
            for (int i = 0; i < elementList.size(); i++)
                System.out.println(elementList.get(i) + " (" + i + ")");

        }
    }






    /**
     * @return TraitSet with same name as elementLabel.
     */
    public TraitSet getElementTrait() {
        if (!traitsProcessed)
            processTraits(m_traitList.get());

        return elementTraitSet;
    }

    /**
     * @return true if TraitSet with same name as elementLabel exists.
     */
    public boolean hasElementTrait() {
        if (getElementTrait() != null)
            return true;
        else
            return false;
    }

    /**
     * Specifically set the type trait set for this tree. A null value simply
     * removes the existing trait set.
     *
     * @param traitSet
     */
    public void setElementTrait(TraitSet traitSet) {
        if (hasElementTrait()) {
            m_traitList.get().remove(elementTraitSet);
        }

        //deal with this when BEAUTi is involved
//        if (traitSet != null) {
//
//            elementTraitInput.setValue(traitSet, this);
//        }

        elementTraitSet = traitSet;
    }

    /**
     * Retrieve the list of unique types identified by the type trait.
     * @return List of unique type trait value strings.
     */
    public List<String> getElementList() {
        if (!traitsProcessed)
            processTraits(m_traitList.get());

        return elementList;
    }


    /**
     * @param type
     * @return string name of given type
     */
    public String getTypeString(int type) {
        if (!traitsProcessed)
            processTraits(m_traitList.get());

        return elementList.get(type);
    }

    /**
     * @param typeString
     * @return integer type corresponding to given type string
     */
    public int getTypeFromString(String typeString) {
        if (!traitsProcessed)
            processTraits(m_traitList.get());

        return elementList.indexOf(typeString);
    }

    /**
     * @return type label to be used in logging.
     */
    public String getElementLabel() {
        return elementLabel;
    }



    @Override
    protected void initArrays() {
        // initialise tree-as-array representation + its stored variant
        m_nodes = new PartitionedTreeNode[nodeCount];
        listNodes(root, m_nodes);
        m_storedNodes = new PartitionedTreeNode[nodeCount];
        Node copy = root.copy();
        listNodes(copy, m_storedNodes);
    }

    /**
     * Convert Volz-type tree to array representation.
     *
     * @param node Root of sub-tree to convert.
     * @param nodes Array to populate with tree nodes.
     */
    private void listNodes(PartitionedTreeNode node, PartitionedTreeNode[] nodes) {
        nodes[node.getNr()] = node;
        node.m_tree = this;
        if (!node.isLeaf()) {
            listNodes(node.getLeft(), nodes);
            if (node.getRight()!=null)
                listNodes(node.getRight(), nodes);
        }
    }

//    /**
//     * Deep copy, returns a completely new Volz-type tree.
//     *
//     * @return a deep copy of this Volz-type tree
//     */
//    @Override
//    public MultiTypeTreeVolz copy() {
//        MultiTypeTreeVolz tree = new MultiTypeTreeVolz();
//        tree.ID = ID;
//        tree.index = index;
//        tree.root = root.copy();
//        tree.nodeCount = nodeCount;
//        tree.internalNodeCount = internalNodeCount;
//        tree.leafNodeCount = leafNodeCount;
//        tree.nTypes = nTypes;
//        tree.elementLabel = elementLabel;
//        return tree;
//    }

    /**
     * Deep copy, returns a completely new multi-type tree.
     *
     * @return a deep copy of this multi-type tree
     */
    @Override
    public PartitionedTree copy() {
        PartitionedTree tree = new PartitionedTree();
        tree.ID = ID;
        tree.index = index;
        tree.root = root.copy();
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
        tree.elementLabel = elementLabel;
        return tree;
    }

    /**
     * Copy all values from an existing Volz-type tree.
     *
     * @param other
     */
    @Override
    public void assignFrom(StateNode other) {
        PartitionedTree pTree = (PartitionedTree) other;

        PartitionedTreeNode[] mtNodes = new PartitionedTreeNode[pTree.getNodeCount()];
        for (int i=0; i<pTree.getNodeCount(); i++)
            mtNodes[i] = new PartitionedTreeNode();

        ID = pTree.ID;
        root = mtNodes[pTree.root.getNr()];
        root.assignFrom(mtNodes, pTree.root);
        root.parent = null;

        nodeCount = pTree.nodeCount;
        internalNodeCount = pTree.internalNodeCount;
        leafNodeCount = pTree.leafNodeCount;
        initArrays();
    }

    /**
     * Copy all values aside from IDs from an existing Volz-type tree.
     *
     * @param other
     */
    @Override
    public void assignFromFragile(StateNode other) {
        PartitionedTree pTree = (PartitionedTree) other;
        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[pTree.root.getNr()];
        Node[] otherNodes = pTree.m_nodes;
        int iRoot = root.getNr();
        assignFromFragileHelper(0, iRoot, otherNodes);
        root.height = otherNodes[iRoot].height;
        root.parent = null;

        PartitionedTreeNode mtRoot = (PartitionedTreeNode)root;
        mtRoot.partitionElementNumber = ((PartitionedTreeNode)(otherNodes[iRoot])).partitionElementNumber;
        mtRoot.cacheIndex = ((PartitionedTreeNode)(otherNodes[iRoot])).cacheIndex;

        if (otherNodes[iRoot].getLeft() != null) {
            root.setLeft(m_nodes[otherNodes[iRoot].getLeft().getNr()]);
        } else {
            root.setLeft(null);
        }
        if (otherNodes[iRoot].getRight() != null) {
            root.setRight(m_nodes[otherNodes[iRoot].getRight().getNr()]);
        } else {
            root.setRight(null);
        }
        assignFromFragileHelper(iRoot + 1, nodeCount, otherNodes);
    }

    /**
     * helper to assignFromFragile *
     */
    private void assignFromFragileHelper(int iStart, int iEnd, Node[] otherNodes) {
        for (int i = iStart; i < iEnd; i++) {
            PartitionedTreeNode sink = (PartitionedTreeNode)m_nodes[i];
            PartitionedTreeNode src = (PartitionedTreeNode)otherNodes[i];
            sink.height = src.height;
            sink.parent = m_nodes[src.parent.getNr()];

            sink.cacheIndex = src.cacheIndex;
            sink.partitionElementNumber = src.partitionElementNumber;

            if (src.getLeft() != null) {
                sink.setLeft(m_nodes[src.getLeft().getNr()]);
                if (src.getRight() != null) {
                    sink.setRight(m_nodes[src.getRight().getNr()]);
                } else {
                    sink.setRight(null);
                }
            }
        }
    }

//    /**
//     * Retrieve total number of allowed types on tree.
//     *
//     * @return total type/deme count.
//     */
//    public int getNTypes() {
//        return nTypes;
//    }

    /**
     * Obtain the number of types defined for this MultiTypeTreeVolz.
     *
     * @return number of types defined for MultiTypeTreeVolz
     */
    public int getNTypes() {
        if (!traitsProcessed)
            processTraits(m_traitList.get());

        return elementList.size();
    }

    /**
     * Check whether typing and timing of tree are sensible.
     *
     * @return true if types and times are "valid"
     */
    public boolean isValid() {
        return timesAreValid( root);
    }

    private boolean timesAreValid(Node node) {
        for (Node child : node.getChildren()) {
            double lastHeight = node.getHeight();
            if (child.getHeight()>lastHeight)
                return false;

            if (!timesAreValid(child))
                return false;
        }

        return true;
    }


    /**
     * Return string representation of Volz-type tree.  We use reflection
     * here to determine whether this is being called as part of writing
     * the state file.
     *
     * @return Volz-type tree string in Newick format.
     */
    @Override
    public String toString() {

        // Behaves differently if writing a state file
        StackTraceElement[] ste = Thread.currentThread().getStackTrace();
        if (ste[2].getMethodName().equals("toXML")) {
            // Use toShortNewick to generate Newick string without taxon labels
            String string = getNormalTree().getRoot().toShortNewick(true);

            // Sanitize ampersands if this is destined for a state file.
            return string.replaceAll("&", "&amp;");
        } else{
            return getNormalTree().getRoot().toSortedNewick(new int[1], true);
        }
    }


    /////////////////////////////////////////////////
    //           StateNode implementation          //
    /////////////////////////////////////////////////
    @Override
    protected void store() {
        storedRoot = m_storedNodes[root.getNr()];
        int iRoot = root.getNr();

        storeNodes(0, iRoot);

        storedRoot.height = m_nodes[iRoot].height;
        storedRoot.parent = null;

        if (root.getLeft()!=null)
            storedRoot.setLeft(m_storedNodes[root.getLeft().getNr()]);
        else
            storedRoot.setLeft(null);
        if (root.getRight()!=null)
            storedRoot.setRight(m_storedNodes[root.getRight().getNr()]);
        else
            storedRoot.setRight(null);

        PartitionedTreeNode mtStoredRoot = (PartitionedTreeNode)storedRoot;
        mtStoredRoot.cacheIndex = ((PartitionedTreeNode)m_nodes[iRoot]).cacheIndex;
        mtStoredRoot.partitionElementNumber = ((PartitionedTreeNode)m_nodes[iRoot]).partitionElementNumber;

        storeNodes(iRoot+1, nodeCount);
    }

    /**
     * helper to store *
     */
    private void storeNodes(int iStart, int iEnd) {
        for (int i = iStart; i<iEnd; i++) {
            PartitionedTreeNode sink = (PartitionedTreeNode)m_storedNodes[i];
            PartitionedTreeNode src = (PartitionedTreeNode)m_nodes[i];
            sink.height = src.height;
            sink.parent = m_storedNodes[src.parent.getNr()];
            if (src.getLeft()!=null) {
                sink.setLeft(m_storedNodes[src.getLeft().getNr()]);
                if (src.getRight()!=null)
                    sink.setRight(m_storedNodes[src.getRight().getNr()]);
                else
                    sink.setRight(null);
            }

            sink.cacheIndex = src.cacheIndex;
            sink.partitionElementNumber = src.partitionElementNumber;
        }
    }

    /////////////////////////////////////////////////
    // Methods implementing the Loggable interface //
    /////////////////////////////////////////////////
    @Override
    public void init(PrintStream printStream) throws Exception {

        printStream.println("#NEXUS\n");
        printStream.println("Begin taxa;");
        printStream.println("\tDimensions ntax="+getLeafNodeCount()+";");
        printStream.println("\t\tTaxlabels");
        for (int i = 0; i<getLeafNodeCount(); i++)
            printStream.println("\t\t\t"+getNodesAsArray()[i].getID());
        printStream.println("\t\t\t;");
        printStream.println("End;");

        printStream.println("Begin trees;");
        printStream.println("\tTranslate");
        for (int i = 0; i<getLeafNodeCount(); i++) {
            printStream.print("\t\t\t"+(getNodesAsArray()[i].getNr()+1)
                    +" "+getNodesAsArray()[i].getID());
            if (i<getLeafNodeCount()-1)
                printStream.print(",");
            printStream.print("\n");
        }
        printStream.print("\t\t\t;");
    }

    @Override
    public void log(int i, PrintStream printStream) {
        printStream.print("tree STATE_"+i+" = ");
        printStream.print(toString());
        printStream.print(";");


    }

    @Override
    public void close(PrintStream printStream) {
        printStream.println("End;");
    }


    /////////////////////////////////////////////////
    // Methods for flat trees rewritten            //
    /////////////////////////////////////////////////

    /**
     * Generates a new tree in which the colours are indicated as meta-data.
     *
     * Caveat: assumes more than one node exists on tree (i.e. leaf != root)
     *
     * @return Normal tree.
     */
    public Tree getNormalTree() {


        Tree flatTree = copy();
        flatTree.initArrays();



        for (Node node : getNodesAsArray()) {

            PartitionedTreeNode pNode = (PartitionedTreeNode)node;

            int nodeNum = node.getNr();

            if (node.isRoot()) {
                Node startNode = flatTree.getNode(nodeNum);
                startNode.setMetaData(elementLabel, ((PartitionedTreeNode)node).getPartitionElementNumber());
                startNode.metaDataString = String.format("%s=%d", elementLabel, pNode.getPartitionElementNumber());
                continue;
            }

            Node startNode = flatTree.getNode(nodeNum);
            startNode.setMetaData(elementLabel, ((PartitionedTreeNode)node).getPartitionElementNumber());
            startNode.metaDataString = String.format("%s=%d", elementLabel, pNode.getPartitionElementNumber());

        }

        return flatTree;
    }


    /**
     * Initialise colours and tree topology from Tree object in which colours
     * are marked by meta-data tags.
     *
     * @param tree
     * @param takeNrsFromTree
     * @throws java.lang.Exception
     */
    public void initFromNormalTree(Tree tree, boolean takeNrsFromTree) throws Exception {

        // Build new coloured tree:

        List<Node> activeTreeNodes = new ArrayList<Node>();
        List<Node> nextActiveTreeNodes = new ArrayList<Node>();
        List<PartitionedTreeNode> activePartTreeNodes = new ArrayList<PartitionedTreeNode>();
        List<PartitionedTreeNode> nextActivePartTreeNodes = new ArrayList<PartitionedTreeNode>();

        // Populate active node lists with root:
        activeTreeNodes.add(tree.getRoot());
        PartitionedTreeNode newRoot = new PartitionedTreeNode();
        activePartTreeNodes.add(newRoot);

        // Initialise counter used to number leaves when takeNrsFromFlatTree
        // is false:
        int nextNr = 0;

        while (!activeTreeNodes.isEmpty()) {

            nextActiveTreeNodes.clear();
            nextActivePartTreeNodes.clear();

            for (int idx = 0; idx<activeTreeNodes.size(); idx++) {
                Node treeNode = activeTreeNodes.get(idx);
                PartitionedTreeNode partTreeNode = activePartTreeNodes.get(idx);

                List<Integer> colours = new ArrayList<>();
                List<Double> times = new ArrayList<>();

                while (treeNode.getChildCount()==1) {
                    int col = (int) Math.round(
                            (Double) treeNode.getMetaData(elementLabel));
                    colours.add(col);
                    times.add(treeNode.getHeight());

                    treeNode = treeNode.getLeft();
                }

                // Order changes from youngest to oldest:
                Collections.reverse(colours);
                Collections.reverse(times);

                switch (treeNode.getChildCount()) {
                    case 0:
                        // Leaf at base of branch
                        if (takeNrsFromTree) {
                            partTreeNode.setNr(treeNode.getNr());
                            partTreeNode.setID(String.valueOf(treeNode.getID()));
                        } else {
                            partTreeNode.setNr(nextNr);
                            partTreeNode.setID(String.valueOf(nextNr));
                            nextNr += 1;
                        }
                        break;

                    case 2:
                        // Non-leaf at base of branch
                        nextActiveTreeNodes.add(treeNode.getLeft());
                        nextActiveTreeNodes.add(treeNode.getRight());

                        PartitionedTreeNode daughter = new PartitionedTreeNode();
                        PartitionedTreeNode son = new PartitionedTreeNode();
                        partTreeNode.addChild(daughter);
                        partTreeNode.addChild(son);
                        nextActivePartTreeNodes.add(daughter);
                        nextActivePartTreeNodes.add(son);

                        break;
                }

                // Set node type at base of multi-type tree branch:
                int nodeType = (int) Math.round(
                        (Double) treeNode.getMetaData(elementLabel));
                partTreeNode.setPartitionElementNumber(nodeType);

                // Set node height:
                partTreeNode.setHeight(treeNode.getHeight());
            }

            // Replace old active node lists with new:
            activeTreeNodes.clear();
            activeTreeNodes.addAll(nextActiveTreeNodes);

            activePartTreeNodes.clear();
            activePartTreeNodes.addAll(nextActivePartTreeNodes);

        }


        // Number internal nodes:
        numberInternalNodes(newRoot, newRoot.getAllLeafNodes().size());

        // Assign tree topology:
        assignFromWithoutID(new PartitionedTree(newRoot));
        initArrays();

    }

    /**
     * Helper method to assign sensible node numbers
     * to each internal node.  This is a post-order traversal, meaning the
     * root is given the largest number.
     *
     * @param node
     * @param nextNr
     * @return
     */
    protected int numberInternalNodes(Node node, int nextNr) {
        if (node.isLeaf())
            return nextNr;

        for (Node child : node.getChildren())
            nextNr = numberInternalNodes(child, nextNr);

        node.setNr(nextNr);
        node.setID(String.valueOf(nextNr));

        return nextNr+1;
    }

    // Check that partition rules are obeyed

    public boolean checkPartitionIntegrity(){
        boolean ok = true;

        for(int i=0; i<nodeCount; i++){
            PartitionedTreeNode node = (PartitionedTreeNode)getNode(i);
            int elementNumber = node.getPartitionElementNumber();

            boolean linked = false;
            for(Node child : node.getChildren()){
                if (((PartitionedTreeNode)child).getPartitionElementNumber()==elementNumber){
                    linked = true;
                }
            }

            if(allowCreep){
                Node parent = node.getParent();
                if (((PartitionedTreeNode)parent).getPartitionElementNumber()==elementNumber){
                    linked = true;
                }

            }

            if(!linked){
                return false;
            }

        }
        return true;
    }

    protected PartitionedTreeNode getElementMRCA(int elementNo){
        HashSet<String> elementTips = new HashSet<>();

        for(Node tip : getExternalNodes()){
            if(((PartitionedTreeNode)tip).getPartitionElementNumber()==elementNo){
                elementTips.add(getTaxonId(tip));
            }

        }

        return (PartitionedTreeNode)TreeUtils.getCommonAncestorNode(this, elementTips);

    }

    protected boolean isRootBlockedBy(int elementNo, int maybeBlockedBy){
        PartitionedTreeNode elementMRCA = getElementMRCA(elementNo);
        PartitionedTreeNode potentialBlockingMRCA = getElementMRCA(maybeBlockedBy);

        PartitionedTreeNode currentNode = elementMRCA;

        while(currentNode!=null){
            currentNode = (PartitionedTreeNode)currentNode.getParent();
            if(currentNode == potentialBlockingMRCA){
                return true;
            }

        }

        return false;

    }

    protected boolean isRootBlocked(int elementNo){
        for(int i=0; i<elementList.size(); i++){
            if(isRootBlockedBy(elementNo, i)){
                return true;
            }
        }
        return false;

    }


    protected boolean isAncestral(PartitionedTreeNode node){
        if(!allowCreep){
            return true;
        }

        int elementNo = node.getPartitionElementNumber();

        Set<String> descendentLeaves = TreeUtils.getDescendantLeaves(this, node);

        for(String id : descendentLeaves){
            for(Node tip : getExternalNodes()){
                if(getTaxonId(tip).equals(id)){
                    if(((PartitionedTreeNode)tip).getPartitionElementNumber()==elementNo){
                        return true;
                    }
                }
            }
        }

        return false;

    }

    public HashSet<PartitionedTreeNode> getTipsInElement(int elementNo){
        HashSet<PartitionedTreeNode> out = new HashSet<>();

        for(Node node : getExternalNodes()){
            if(((PartitionedTreeNode)node).getPartitionElementNumber()==elementNo){
                out.add((PartitionedTreeNode)node);
            }
        }

        return out;
    }

    public HashSet<PartitionedTreeNode> getNodesInElement(int elementNo){
        HashSet<PartitionedTreeNode> out = new HashSet<>();

        for(Node node : m_nodes){
            if(((PartitionedTreeNode)node).getPartitionElementNumber()==elementNo){
                out.add((PartitionedTreeNode)node);
            }
        }

        return out;

    }

    public HashSet<PartitionedTreeNode> getNodesInSameElement(PartitionedTreeNode node){
        int elementNo = node.getPartitionElementNumber();
        return getNodesInElement(elementNo);
    }




    /////////////////////////////////////////////////
    // Serialization and deserialization for state //
    /////////////////////////////////////////////////

    /**
     * reconstruct tree from XML fragment in the form of a DOM node *
     * @param node
     */
    @Override
    public void fromXML(org.w3c.dom.Node node) {
        try {
            String sNewick = node.getTextContent().replace("&", "");

            TreeParser parser = new TreeParser();
            parser.initByName(
                    "IsLabelledNewick", false,
                    "offset", 0,
                    "adjustTipHeights", false,
                    "singlechild", true,
                    "newick", sNewick);
            //parser.m_nThreshold.setValue(1e-10, parser);
            //parser.m_nOffset.setValue(0, parser);

            initFromNormalTree(parser, true);

            initArrays();
        } catch (Exception ex) {
            Logger.getLogger(PartitionedTree.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
