package tree;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;

import java.util.List;

/**
 * Created by twoseventwo on 07/12/2015.
 */

@Description("A tree whose nodes are partitioned as in Hall et al (2015). No assumption is made in this class" +
        "about what the partition elements mean")
public class PartitionedTree extends Tree {

    protected String elementLabel;

    // the information that determines the partitions. In an outbreak these will be hosts/cases

    protected TraitSet pElementTraitSet;

    protected List<String> pElementList;

    public PartitionedTree() { };

    public PartitionedTree(PElementMemberNode rootNode) {
        //public MultiTypeTreeVolz(Node rootNode) {

        if (!(rootNode instanceof PElementMemberNode)) {
            throw new IllegalArgumentException("Attempted to instantiate "
                    + "multi-type tree with regular root node.");
        }

        setRoot(rootNode);
        initArrays();
    }

    @Override
    public void initAndValidate() throws Exception {

        super.initAndValidate();

        for(int element = 0; element<pElementList.size(); element++ ){
            int tipCount = 0;

            for(int tipNo = 0; tipNo < getLeafNodeCount(); tipNo++){
                if(pElementTraitSet.)


            }

        }



    }




}
