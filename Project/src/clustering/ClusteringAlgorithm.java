package clustering;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;
import exceptions.BioInfException;
import pairwisealignment.NeedlemanWunsch;
import pairwisealignment.PairwiseAlignment;
import pairwisealignment.ScoringFunction;
import utilities.ProteinSequence;
import utilities.ProteinSequenceList;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;

public abstract class ClusteringAlgorithm extends BioinfAlgorithm {

  public ProteinSequenceList psList;
  public NeedlemanWunsch nw;
  protected DistanceFunction df;
  /**
   * A String of length 1 that will be used to name the next cluster. The
   * function genNewClusterName will iterate through the alphabet and then start
   * using integers.
   */
  public String nextClusterName = "A";
  /**
   * The String here contains a concatenation of both cluster names with a space
   * in between. The double is the distance.
   */
  protected ConcurrentHashMap<Pair, Double> dm;

  /**
   * A simple doubly linked binary tree.
   * 
   * @author Max
   * 
   */
  public class Cluster {

    public String id;
    public Cluster lChild;
    public Cluster rChild;
    public Cluster parent;
    public double height;

    public Cluster(String id, Cluster lChild, Cluster rChild, Double height) {
      this.id = id;
      this.lChild = lChild;
      this.rChild = rChild;
      this.height = height;
    }

    public Cluster(String id) {
      this.id = id;
      this.lChild = null;
      this.rChild = null;
      this.height = 0.0;
    }

    public void setParent(Cluster parent) {
      this.parent = parent;
    }

    @Override
    public boolean equals(Object obj) {
      if (!(obj instanceof Cluster)) {
        return false;
      }
      Cluster c = (Cluster) obj;
      return c.id.equals(this.id);
    }

    @Override
    public int hashCode() {
      return super.hashCode();
    }

    public boolean isLeaf() {
      return (this.lChild == null && this.rChild == null);
    }

    /**
     * Recursive descent to find subTreeSize
     * 
     * @return the number of nodes in the rooted three beginning at this
     */
    public int subTreeSize() {
      int lDaughters = 0, rDaughters = 0;
      if (this.lChild != null) {
        lDaughters = this.lChild.subTreeSize();
      }
      if (this.rChild != null) {
        rDaughters = this.rChild.subTreeSize();
      }

      return 1 + lDaughters + rDaughters;
    }

    private Cluster() {
      this.id = null;
    }

    public Cluster findCluster(Cluster c) {
      // Check the current node
      if (this.equals(c)) {
        return this;
      } else if (this.isLeaf()) {
        return new Cluster();
      } else {
        Cluster result = new Cluster();
        // Search the left branch
        if (this.lChild != null) {
          result = this.lChild.findCluster(c);
        }
        // Search the right branch
        if (result.id == null && this.rChild != null) {
          result = this.rChild.findCluster(c);
        }
        return result;
      }

    }
  }

  /**
   * Simple Pair class to index the DistanceMatrix
   * 
   * @author Max
   * 
   * @param <X>
   * @param <Y>
   */
  protected final class Pair<X, Y> {
    public X left;
    public Y right;

    public Pair(X left, Y right) {
      this.left = left;
      this.right = right;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj)
        return true;
      if (obj == null)
        return false;
      if (getClass() != obj.getClass())
        return false;
      Pair other = (Pair) obj;
      if (!getOuterType().equals(other.getOuterType()))
        return false;
      // If l = l and r = r OR l = r and r = l
      if (other.left.equals(this.left) && other.right.equals(this.right)
          || other.right.equals(this.left) && other.left.equals(this.right)) {
        return true;
      } else {
        return false;
      }
    }

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + getOuterType().hashCode();
      String lString = (String) this.left;
      String rString = (String) this.right;

      if (lString.compareTo(rString) <= 0) {
        result = prime * result + ((left == null) ? 0 : left.hashCode());
        result = prime * result + ((right == null) ? 0 : right.hashCode());
      } else {
        result = prime * result + ((right == null) ? 0 : right.hashCode());
        result = prime * result + ((left == null) ? 0 : left.hashCode());
      }

      return result;
    }

    private ClusteringAlgorithm getOuterType() {
      return ClusteringAlgorithm.this;
    }
  }

  public ArrayList<Cluster> clusterTree;

  @Override
  public abstract String getName();

  @Override
  public abstract String getDescription();

  @Override
  public Vector<AlgorithmParameter> getInputParameters() {
    return super.parameters;
  }

  /*
   * (non-Javadoc)
   * 
   * @see gui.BioinfAlgorithm#run(java.util.Vector)
   */
  @Override
  public String run(Vector<AlgorithmParameter> params) {

    if (params.size() != 3) {
      return "This class expects 3 parameters, <StringList of Sequences>, "
          + " <ScoringFunctionName> <Integer - Gap Penalty>  ";
    }

    // Initialize state variables
    try {
      // Parse the first argument -- the FASTA sequences
      this.parseSequences(params.get(0));
      // Remove it, so that we can feed the remaining arguments to NW
      params.remove(0);
      nw = new NeedlemanWunsch();
      nw.internal_setup(params);
      // Use the same ScoringFunction and GapPenalty
      this.df = new DistanceFunction(nw.sf, nw.gapPenalty);
      // Build a symmetric matrix with distance between all possible alignments
      this.initializeDistanceMatrix();
      // For every sequence create a cluster
      this.initializeClusters();
      // Until every cluster is part of one big tree
      while (clusterTree.size() > 1) {
        // Find the two closest clusters
        Cluster[] clusterPair = findClosestClusterPair();
        // Make a new cluster from the union of these two.
        // and set the height as the distance between them / 2.
        Cluster union =
            new Cluster(getClusterName(), clusterPair[0], clusterPair[1],
                getDistance(clusterPair[0].id, clusterPair[1].id) / 2);
        // Cluster is doubly linked, to allow for easier backwards traversal.
        clusterPair[0].setParent(union);
        clusterPair[1].setParent(union);
        clusterTree.add(union);
        // Set the values in the distance matrix to be equal to whatever
        // formula is defined in the child-class
        // and then remove the clusters from the distance matrix
        updateDistanceMatrix(clusterPair[0].id, clusterPair[1].id, union.id);
        // remove the clusters from the cluster tree
        clusterTree.remove(clusterPair[0]);
        clusterTree.remove(clusterPair[1]);
      }

    } catch (BioInfException e) {

      return e.getMessage().toString();
    }
    
    String output = this.buildParameterOutputString();
    output += treeToNewick(this.clusterTree.get(0)) + ";" +  "\n";
    // Create a set of clusters, one for every node in the distance matrix.
    return output;
  }

  /**
   * @return A string in Newick format of clusterTree
   */
  public String treeToNewick(Cluster tree) {

    // Equivalent to depth first search
    // For every node
    String output = "";
    output = tree.id + output;
    // If the cluster is not a leaf (has children)
    if (!tree.isLeaf()) {
      output = ")" + output;
      if (tree.lChild != null) {
        output = treeToNewick(tree.lChild) + output;
      }
      if (tree.lChild != null && tree.rChild != null) {
        output = "," + output;
      }
      if (tree.rChild != null) {
        output = treeToNewick(tree.rChild) + output;
      }
      output = "(" + output;
    }

    return output;
  }

  /**
   * For each proteinSequence in psList, make a cluster and add it to the tree.
   */
  public void initializeClusters() {
    this.clusterTree = new ArrayList<Cluster>();
    for (ProteinSequence ps : this.psList.toArray()) {
      Cluster c = new Cluster(ps.id);
      this.clusterTree.add(c);
    }
  }

  /**
   * For every unique combination of DIFFERENT protein sequences determine the
   * distance and put it into the distance matrix.
   * 
   * @throws BioInfException
   */
  public void initializeDistanceMatrix() throws BioInfException {
    ProteinSequence[] arr = psList.toArray();
    this.dm = new ConcurrentHashMap<Pair, Double>();
    // for each sequence in psList
    for (int s1 = 0; s1 < arr.length - 1; s1++) {
      // for each other sequence
      for (int s2 = s1 + 1; s2 < arr.length; s2++) {
        // Find optimal pairwise alignment using Needleman-Wunsch
        PairwiseAlignment opt_alignment = nw.internal_run(arr[s1], arr[s2]);
        String s1id = arr[s1].id;
        String s2id = arr[s2].id;
        // Find distance using Feng-Doolittle
        Double distance =
            this.df.computeDistance(opt_alignment.a1, opt_alignment.a2);
        // Enter this distance in the distance matrix and make it symmetric
        dm.put(new Pair(s1id, s2id), distance);
      }
    }
  }

  /**
   * Parses any number of FASTA sequences.
   * 
   * @param param
   * @throws BioInfException
   */
  public void parseSequences(AlgorithmParameter param) throws BioInfException {
    if (param.type != gui.StringList.class
        || param.data.getClass() != gui.StringList.class) {
      throw new BioInfException("Expected StringList", "but found "
          + param.type.toString());
    }
    String[] lines = param.data.toString().trim().split("\n");
    this.psList = ProteinSequenceList.parse(lines);
  }

  public int parseGapPenalty(AlgorithmParameter gp) throws BioInfException {
    if (!(gp.type == Integer.class && gp.data.getClass() == Integer.class)) {
      throw new BioInfException("Expected Integer", "but found "
          + gp.data.getClass().toString());
    }
    if ((Integer) gp.data < 0) {
      throw new BioInfException("Expected Integer >= 0", "but found "
          + (Integer) gp.data);
    } else {
      return (Integer) gp.data;
    }
  }

  public ScoringFunction parseScoringFunction(AlgorithmParameter p)
    throws BioInfException {
    if (!(p.type == Boolean.class)) {
      throw new BioInfException("Expected Boolean", "but found "
          + p.type.toString());
    }

    return new ScoringFunction((Boolean) p.data);
  }

  /**
   * Look up the string c1Namec2Name in the distance matrix and return the
   * result.
   * 
   * @param s1Name
   * @param s2Name
   * @return
   */
  public Double getDistance(String c1Name, String c2Name) {
    return this.dm.get(new Pair(c1Name, c2Name));
  }

  /**
   * Searches the distance matrix for any reference to clusterName and removes
   * those entries.
   * 
   * @param clusterName
   */
  protected void removeFromDM(String clusterName) {
    ArrayList<Pair> keys = Collections.list(this.dm.keys());
    for (Pair key : keys) {
      String clusterId1 = (String) key.left;
      String clusterId2 = (String) key.right;
      if (clusterId1.equals(clusterName) || clusterId2.equals(clusterName)) {
        this.dm.remove(key);
      }
    }
  }

  /**
   * When implemented in WPGMA and UPGMA, this will use different means.
   * Responsible for updating the distance matrix with new values and removing
   * old entries.
   */
  public abstract void updateDistanceMatrix(String oldCluster1,
    String oldCluster2, String unionId);

  /**
   * Iterates through all values in the distance matrix and finds the smallest.
   * 
   * @return String array of length 2. Elts are names of two closest clusters.
   */
  protected Cluster[] findClosestClusterPair() {
    double minDistance = Double.MAX_VALUE;
    Pair minKey = new Pair("", "");

    // For entry in the distance matrix,
    // Which means, every pair of clusters
    for (Entry<Pair, Double> clusterPair : this.dm.entrySet()) {
      if (clusterPair.getValue() < minDistance) {
        minDistance = clusterPair.getValue();
        minKey = clusterPair.getKey();
      }
    }
    // Parse the key, which is "id1 id2"
    String clusterId1 = (String) minKey.left;
    String clusterId2 = (String) minKey.right;

    Cluster[] clusterPair = new Cluster[2];
    // Iterate through the clusterTree in order to find the cluster objects.
    for (Cluster c : clusterTree) {
      if (c.id.equals(clusterId1)) {
        clusterPair[0] = c;
      } else if (c.id.equals(clusterId2)) {
        clusterPair[1] = c;
      }
    }

    return clusterPair;
  }

  /**
   * Iterates through the alphabet before using the natural numbers to generate
   * unique names for clusters.
   * 
   * @return a unique name for the next cluster
   */
  public String getClusterName() {
    String oldName = this.nextClusterName;

    if (Character.isAlphabetic(this.nextClusterName.charAt(0))) {
      char newChar = (char) (this.nextClusterName.charAt(0) + 1);
      this.nextClusterName = newChar + "";
      if (this.nextClusterName.equals("[")) {
        this.nextClusterName = "1";
      }
    } else {
      this.nextClusterName = Integer.valueOf(this.nextClusterName) + 1 + "";
    }

    return oldName;
  }
  
  private String buildParameterOutputString() {
    StringBuilder s = new StringBuilder();
    for (AlgorithmParameter ap : this.parameters) {
      s.append(ap.name + ": " + ap.data + "\n");
    }
    s.append("\n");
    return s.toString();
  }

}
