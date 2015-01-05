package clustering;

import java.util.ArrayList;
import java.util.Collections;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;

public class UPGMA extends ClusteringAlgorithm {

  public UPGMA() {
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
            ">s1\nILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN"
                + "\n>s2\nRDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT\n"
                + ">s3\nISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA\n"
                + ">s4\nRRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD")));

    super.parameters
    .add(new AlgorithmParameter(
        "UseBlosum62Matrix",
        "Unchecking this box makes the program use the pam250 scoring matrix.",
        Boolean.class, new Boolean(true)));

    super.parameters.add(new AlgorithmParameter("Gap Penalty",
    "This value will be SUBTRACTED during the scoring process"
        + " every time a gap is inserted into an alignment.",
    Integer.class, 6));
  }

  @Override
  public String getName() {
    return "UPGMA";
  }

  @Override
  public String getDescription() {
    return new String("This class excepts any number of FASTA sequences,"
        + " and then uses the Feng-Doolittle distance algorithm"
        + " to progressively find the closest clusters until"
        + " a binary tree is created.");
  }

  /*
   * (non-Javadoc) UPGMA uses unweighted distance measure by multiplying by
   * cardinality of clusters and then dividing by the sum of cardinality.
   * 
   * @see
   * multisequencealignment.ClusteringAlgorithm#updateDistanceMatrix(java.lang
   * .String, java.lang.String, java.lang.String)
   */
  @Override
  public void updateDistanceMatrix(String cId1, String cId2, String unionId) {
    // UPGMA uses magnitude of two clusters
    int magC1 = 0, magC2 = 0;
    for (Cluster c : this.clusterTree) {
      if (c.id.equals(cId1)) {
        magC1 = c.subTreeSize();
      } else if (c.id.equals(cId2)) {
        magC2 = c.subTreeSize();
      }
    }

    ArrayList<Pair> keys = Collections.list(this.dm.keys());
    for (Pair key : keys) {
      String clusterId1 = (String) key.left;
      String clusterId2 = (String) key.right;
      // If clusterId1 is the matching cluster id
      if ((clusterId1.equals(cId1) && !clusterId2.equals(cId2))
          || (clusterId1.equals(cId2) && !clusterId2.equals(cId1))) {
        // Then clusterId2 must be part of the new distance matrix entry
        // If the entry is null, make one, using the average
        if (this.dm.get(new Pair(unionId, clusterId2)) == null) {
          this.dm.put(new Pair(unionId, clusterId2), this.dm.get(key) * magC1
              / (magC1 + magC2));
          // Otherwise, add to the existing value
        } else {
          Double newValue =
              this.dm.get(new Pair(unionId, clusterId2)) + this.dm.get(key)
                  * magC1 / (magC1 + magC2);
          this.dm.put(new Pair(unionId, clusterId2), newValue);
        }
        // Else if clusterId2 is the matching cluster
      } else if ((clusterId2.equals(cId1) && !clusterId1.equals(cId2))
          || (clusterId2.equals(cId2) && !clusterId1.equals(cId1))) {
        // Then clusterId2 must be part of the new dist. matrix entry
        if (this.dm.get(new Pair(unionId, clusterId1)) == null) {
          this.dm.put(new Pair(unionId, clusterId1), this.dm.get(key) * magC2
              / (magC1 + magC2));
        } else {
          Double newValue =
              this.dm.get(new Pair(unionId, clusterId1)) + this.dm.get(key)
                  * magC2 / (magC1 + magC2);
          this.dm.put(new Pair(unionId, clusterId1), newValue);
        }
      }
    }

    this.removeFromDM(cId1);
    this.removeFromDM(cId2);

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
    UPGMA myInstance = new UPGMA();

    // run the example the instance with the default parameters
    BioinfAlgorithm.runAlgorithmDefaults(myInstance);
  }

}
