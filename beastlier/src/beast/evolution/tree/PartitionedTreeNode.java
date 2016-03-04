
package beast.evolution.tree;


import beast.core.Description;


/**
 *
 * @author Nicola De Maio
 */
@Description("A node in a tree with node partitions")
public class PartitionedTreeNode extends Node {


    int partitionElementNumber = -1;
    int cacheIndex = -1;


    /**
     * @return position of the lineage in the probability calculation cache
     */
    public int getCacheIndex() {
        return cacheIndex;
    }

    /**
     * Sets the index of the node in the cache
     *
     * @param cacheIndex New cache index.
     */
    public void setCacheIndex(int cacheIndex) {
        this.cacheIndex = cacheIndex;
    }

    /**
     * Obtain partition element number for.
     *
     * @return type value at this node
     */
    public int getPartitionElementNumber() {
        return partitionElementNumber;
    }

    /**
     * Set partition element number.
     *
     * @param partitionElementNumber New node type.
     */
    public void setPartitionElementNumber(int partitionElementNumber) {
        this.partitionElementNumber = partitionElementNumber;
    }

    /**
     * Remove all type changes from branch above node.
     */
    public void clearChanges() {
        startEditing();
    }

    /**
     * @return shallow copy of node
     */
    public PartitionedTreeNode shallowCopy() {
        PartitionedTreeNode node = new PartitionedTreeNode();
        node.height = height;
        node.parent = parent;
        node.children.addAll(children);

        node.partitionElementNumber = partitionElementNumber;
//        if (likes==null) {node.likes = null;}else{
//        	node.likes = Arrays.copyOf(likes,likes.length);
//        }
        node.cacheIndex = cacheIndex;

        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.ID = ID;

        return node;
    }

    /**
     * **************************
     * Methods ported from Node *
     ***************************
     */


    /**
     * methods from Node.java modified to include type structure
     */


    /**
     * @return (deep) copy of node
     */
    @Override
    public PartitionedTreeNode copy() {
        PartitionedTreeNode node = new PartitionedTreeNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.parent = null;
        node.ID = ID;
        node.cacheIndex = cacheIndex;
        node.partitionElementNumber = partitionElementNumber;
//        if (likes==null) {node.likes = null;}else{
//        	node.likes = Arrays.copyOf(likes,likes.length);
//        }

        if (getLeft()!=null) {
            node.setLeft(getLeft().copy());
            node.getLeft().setParent(node);
            if (getRight()!=null) {
                node.setRight(getRight().copy());
                node.getRight().setParent(node);
            }
        }

        return node;
    }

    /**
     * assign values from a tree in array representation *
     * @param nodes
     * @param node
     */
    @Override
    public void assignFrom(Node[] nodes, Node node) {
        height = node.height;
        labelNr = node.labelNr;
        metaDataString = node.metaDataString;
        parent = null;
        ID = node.getID();

        PartitionedTreeNode mtNode = (PartitionedTreeNode)node;
        cacheIndex = mtNode.cacheIndex;
        partitionElementNumber = mtNode.partitionElementNumber;
//        if (mtNode.likes==null) {likes = null;}else{
//        	likes = Arrays.copyOf(mtNode.likes,mtNode.likes.length);
//        }

        if (node.getLeft()!=null) {
            setLeft(nodes[node.getLeft().getNr()]);
            getLeft().assignFrom(nodes, node.getLeft());
            getLeft().parent = this;
            if (node.getRight()!=null) {
                setRight(nodes[node.getRight().getNr()]);
                getRight().assignFrom(nodes, node.getRight());
                getRight().parent = this;
            }
        }
    }
}
