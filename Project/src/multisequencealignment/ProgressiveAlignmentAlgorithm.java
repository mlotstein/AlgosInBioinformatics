package multisequencealignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

import pairwisealignment.PairwiseAlignment;
import clustering.ClusteringAlgorithm.Cluster;
import clustering.ClusteringAlgorithm;
import clustering.UPGMA;
import clustering.WPGMA;
import exceptions.BioInfException;
import utilities.Constants;
import utilities.ProteinSequence;
import utilities.ProteinSequenceList;
import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;

public class ProgressiveAlignmentAlgorithm extends BioinfAlgorithm {
  /**
   * The list of sequences from input
   */
  ProteinSequenceList psList;

  /**
   * The aligned sequences, with gap characters replaced.
   */
  ProteinSequenceList g;
  ClusteringAlgorithm ca;

  public ProgressiveAlignmentAlgorithm() {

    super.parameters.add(new AlgorithmParameter("Use UPGMA",
        "By unchecking this option, the program will use WPGMA.",
        Boolean.class, new Boolean(false)));

    super.parameters
        .add(new AlgorithmParameter(
            "Sequence Data",
            "Any number of correctly formatted FASTA sequences."
                + " The correct format for each is the sequence is:"
                + " '><sequence name><newline><sequence data><newline>"
                + " and can be repeated consecutively any number of times."
                + " Beware though: this algorithm requires O(m^n) time and"
                + " space to complete, where m = the average sequence "
                + " length and n = the number of sequences.",
            StringList.class,
            new StringList(
                "> S1\nMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA\n>"
                    + " S2\nMTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLPQDVEEFFEGPSEA\n> "
                    + " S3\nEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEEFFEGPSEA")));

    super.parameters
        .add(new AlgorithmParameter(
            "UseBlosum62Matrix",
            "Unchecking this box makes the program use the pam250 scoring matrix.",
            Boolean.class, new Boolean(true)));

    super.parameters.add(new AlgorithmParameter("Gap Penalty",
        "This value will be SUBTRACTED during the scoring process"
            + " every time a gap is inserted into an alignment.",
        Integer.class, 2));
  }

  @Override
  public String getName() {
    return "Progressive Alignment using Feng-Doolittle.";
  }

  @Override
  public String getDescription() {
    return new String("This class excepts any number of FASTA sequences,"
        + " and then uses a clustering algorithm to build a guide tree. Then"
        + " it iterates through the guide tree, performing pairwise"
        + " alignments with the most similar sequences.");
  }

  @Override
  public Vector<AlgorithmParameter> getInputParameters() {
    return super.parameters;
  }

  @Override
  public String run(Vector<AlgorithmParameter> params) {
    if (params.size() != 4) {
      return "This class expects 4 parameters, <ClusteringAlgorithmName> "
          + "<StringList of Sequences>, "
          + " <ScoringFunctionName> <Integer - Gap Penalty>  ";
    }

    // Parse the clustering algorithm
    try {
      this.parseClusteringAlgorithm(params.get(0));
    } catch (BioInfException e) {
      e.printStackTrace();
    }

    // Parse the second argument -- the FASTA sequences
    try {
      this.parseSequences(params.get(1));
    } catch (BioInfException e) {
      e.printStackTrace();
    }

    // Remove the name of the clustering algorithm
    params.remove(0);

    // Generate the guide tree.
    this.generateGuideTree(params);

    // The guide tree was constructed such that the cluster-nodes were
    // generated alphabetically. Thus, the two closest clusters have the
    // name "A", the next "B" and so forth.
    this.sortProteinSequenceList(this.parseGuideTree());

    g = new ProteinSequenceList();

    // Use Needleman-Wunsch to align the top two sequences
    PairwiseAlignment pa = this.ca.nw.internal_run(psList.pop(), psList.pop());
    // Take the two alignments and replace gap characters with Neutral Char
    // and then add them to G.
    g.add(OAGAAG(pa.a1, pa.id1));
    g.add(OAGAAG(pa.a2, pa.id2));

    while (!psList.isEmpty()) {
      // Take the next sequence to be aligned from the guide tree.
      ProteinSequence ps = psList.pop();
      // Find the best pairwise alignment between remaining sequences
      PairwiseAlignment bestAlignment = new PairwiseAlignment();
      Float bestScore = Float.NEGATIVE_INFINITY;

      for (ProteinSequence as : g.toArray()) {
        // NOTE: ps -- new sequence -- is in FIRST position
        // and as -- sequence taken g -- is in SECOND position
        PairwiseAlignment newAlignment = this.ca.nw.internal_run(ps, as);
        float newScore =
            this.ca.nw.sf.scoreAlignment(newAlignment.a1, newAlignment.a2,
                this.ca.nw.gapPenalty);
        if (newScore > bestScore) {
          bestScore = newScore;
          bestAlignment = newAlignment;
        }
      }
      // The new alignment will likely have gaps. Update the alignments in
      // g to reflect this.
      this.insertGapsIntoG(bestAlignment.a2, bestAlignment.id2);
      // add both new sequences to g, after replacing gaps
      g.add(OAGAAG(bestAlignment.a1, bestAlignment.id1));
      g.add(OAGAAG(bestAlignment.a2, bestAlignment.id2));
    }

    replaceStars();

    return buildOutputString();
  }

  public Double sumOfSquares() {
    Double sumOfSquares = 0.0;
    // For every column
    // For every sequence
    // For every other sequence
    for (int seq1 = 0; seq1 < this.g.size() - 1; seq1++) {
      for (int seq2 = seq1 + 1; seq2 < this.g.size(); seq2++) {
        sumOfSquares +=
            this.ca.nw.sf.scoreAlignment(
                this.g.get(seq1).seq.toArray(new Character[0]),
                this.g.get(seq2).seq.toArray(new Character[0]),
                this.ca.nw.gapPenalty);
      }
    }

    return sumOfSquares;
  }
  
  private String buildParameterOutputString() {
    StringBuilder s = new StringBuilder();
    for (AlgorithmParameter ap : this.parameters) {
      s.append(ap.name + ": " + ap.data + "\n");
    }
    s.append("\n");
    return s.toString();
  }

  public String buildOutputString() {

    StringBuilder output = new StringBuilder();
    output.append(buildParameterOutputString());
    
    Double sumOfSquares = this.sumOfSquares();
    output.append("Sum of Squares: " + sumOfSquares + "\n");
    // For each alignment of the form String1, String2
    // ... The first line is <String 1 id> <alignment characters>
    // ... The second line is a space buffer
    int maxIdLength = 0;
    for (int seqNum = 0; seqNum < this.g.size(); seqNum++) {
      if (this.g.get(seqNum).id.length() + 1 > maxIdLength) {
        maxIdLength = this.g.get(seqNum).id.length() + 1;
      }
    }

    // Each chunk will be of size 80 - maxIdLength
    // Total number of chunks = length of a sequence / chunkSize
    int chunkSize = 80 - maxIdLength;
    int seqLength = g.get(0).seq.size();
    int chunkStart = 0;
    int chunkEnd = Math.min(chunkSize - 1, seqLength - 1);
    boolean reachedEnd = false;
    while (!reachedEnd) {
      // For every sequence, find the next chunk and add it to output
      output.append("\n");
      for (int sNum = 0; sNum < this.g.size(); sNum++) {
        String id = this.g.get(sNum).id;
        char[] idPad = new char[maxIdLength - this.g.get(sNum).id.length()];
        Arrays.fill(idPad, ' ');
        StringBuilder curChunk = new StringBuilder();
        for (int cNum = chunkStart; cNum <= chunkEnd; cNum++) {
          curChunk.append(this.g.get(sNum).seq.get(cNum));
        }
        output.append(id + String.valueOf(idPad) + curChunk.toString() + "\n");
      }
      chunkStart = chunkEnd + 1;
      chunkEnd = Math.min(chunkStart + chunkSize, seqLength - 1);
      if (chunkStart >= seqLength) {
        reachedEnd = true;
      }
    }
    return output.toString();
  }

  /**
   * Run the clustering algorithm in order to generate a guide tree.
   * @param params
   */
  protected void generateGuideTree(Vector<AlgorithmParameter> params) {
    // Run the clustering algorithm to generate the guide tree
    this.ca.run(params);
  }

  /**
   * Because the naming convention for clusters during guide tree creation is to
   * progress through the alphabet (and then, natural numbers), to get the
   * order, simply recursively search for ids until the cluster tree is empty.
   * NOTE: this is DESTRUCTIVE
   * 
   * @return a string of ProteinSequence ids, corresponding to those in psList,
   *         ordered by proximity in the guide tree
   */
  public String[] parseGuideTree() {
    String[] ids = new String[this.psList.size()];
    int numIds = 0;
    // The first cluster name we will search for.
    this.ca.nextClusterName = "A";
    Cluster root = this.ca.clusterTree.get(0);
    // While the cluster tree's root node still has children
    while (!root.isLeaf()) {
      // Search for the next node and increment name counter
      Cluster c =
          root.findCluster(this.ca.new Cluster(this.ca.getClusterName()));
      // This node's children are closest clusters
      // Add their ids next
      if (c.lChild != null) {
        ids[numIds] = c.lChild.id;
        numIds++;
      }
      if (c.rChild != null) {
        ids[numIds] = c.rChild.id;
        numIds++;
      }
      // Then, remove this node from the tree.
      // There are 2 possibilities.
      // Either, we are at the root node, in which case make sure both
      // children are null
      if (c.parent == null) {
        c.lChild = null;
        c.rChild = null;
        // other wise, prune the tree above the current node
      } else {
        Cluster parent = c.parent;
        if (parent.lChild.equals(c)) {
          parent.lChild = null;
        } else {
          parent.rChild = null;
        }
      }
    }

    return ids;
  }

  /**
   * A simple selection sort, using the ids array as the guide
   * 
   * @param ids
   *          proper order of psList ids
   */
  public void sortProteinSequenceList(String[] ids) {
    for (int idNum = 0; idNum < ids.length - 1; idNum++) {
      ProteinSequence[] psArray = psList.toArray();
      boolean notFound = true;
      for (int psNum = idNum; notFound && psNum < psArray.length; psNum++) {
        if (psArray[psNum].id.equals(ids[idNum])) {
          notFound = false;
          ProteinSequence temp = psList.get(idNum);
          psList.set(idNum, psList.get(psNum));
          psList.set(psNum, temp);
        }
      }
    }

  }

  /**
   * Parses input in order to determine whether to instantiate UPGMA or WPGMA.
   * 
   * @param clusterType
   * @throws BioInfException
   */
  public void parseClusteringAlgorithm(AlgorithmParameter clusterType)
    throws BioInfException {
    if (clusterType.type != Boolean.class
        || clusterType.data.getClass() != Boolean.class) {
      throw new BioInfException("Expected Boolean", "but found "
          + clusterType.type.toString());
    }

    Boolean isUPGMA = (Boolean) clusterType.data;
    if (isUPGMA) {
      this.ca = new UPGMA();
    } else {
      this.ca = new WPGMA();
    }
  }

  /**
   * The group of proteinSequences g needs to have a gap inserted into it for
   * every gap in the alignment s.
   * 
   * @param s
   *          an array of Characters
   * @param id
   *          a string id of a sequence in g that will be removed.
   */
  private void insertGapsIntoG(Character[] s, String id) {
    ArrayList<Integer> gapPositions = new ArrayList<Integer>();
    for (int cNum = 0; cNum < s.length; cNum++) {
      if (s[cNum] == Constants.GAP_SYMBOL) {
        gapPositions.add(cNum);
      }
    }

    // Remove the id that was matched
    g.removeById(id);
    // For gap in the alignment s
    for (Integer gapPos : gapPositions) {
      // For every sequence in the psList g
      for (int psNum = 0; psNum < g.size(); psNum++) {
        // Insert a neutral character at that position
        g.get(psNum).seq.add(gapPos, Constants.NEUTRAL_CHAR);

      }
    }
  }

  /**
   * "Once a Gap, Always a Gap"
   * 
   * @param arr
   * @return the same alignment as a ProteinSequence, but with gap characters
   *         replaced
   */
  private static ProteinSequence OAGAAG(Character[] arr, String id) {
    ProteinSequence ps = new ProteinSequence();
    ps.id = id;
    try {
      for (Character c : arr) {
        if (c.equals(Constants.GAP_SYMBOL) || c.equals(Constants.NEUTRAL_CHAR)) {
          ps.addStar();
        } else {
          ps.add(c);
        }
      }
    } catch (BioInfException e) {
      e.printStackTrace();
    }
    return ps;
  }

  /**
   * Replace the '*' character with the Gap Symbol
   * 
   */
  private void replaceStars() {
    ProteinSequence[] psArr = this.g.toArray();
    for (int psNum = 0; psNum < psArr.length; psNum++) {
      ProteinSequence ps = psArr[psNum];
      Character[] cArr = ps.seq.toArray(new Character[0]);
      for (int cNum = 0; cNum < cArr.length; cNum++) {
        if (cArr[cNum].equals(Constants.NEUTRAL_CHAR)) {
          ps.seq.set(cNum, Constants.GAP_SYMBOL);
        }
      }
      this.g.set(psNum, ps);
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

  /**
   * Creates an instance of this class and calls the run method using the
   * default parameters.
   * 
   * @param args
   *          program parameters (completely ignored)
   */
  public static void main(String[] args) {
    // create an instance of this class
    ProgressiveAlignmentAlgorithm myInstance =
        new ProgressiveAlignmentAlgorithm();

    // run the example the instance with the default parameters
    BioinfAlgorithm.runAlgorithmDefaults(myInstance);
  }

}
