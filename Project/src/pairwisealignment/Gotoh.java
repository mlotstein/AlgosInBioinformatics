package pairwisealignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Stack;
import java.util.Vector;

import exceptions.BioInfException;
import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;
import utilities.Constants;

public class Gotoh extends PairwiseAlignmentAlgorithm {
  /**
   * A matrix for choosing bet. aligning or putting gaps in either sequence.
   */
  int[][] m;
  /**
   * A matrix showing value of prefixes ending with a gap in seq2.
   */
  int[][] p;
  
  /**
   *  A matrix showing value of prefixes ending with a gap in seq1.
   */
  int[][] q;
  TRACE[][][] t_m;
  TRACE[][][] t_p;
  TRACE[][][] t_q;
  int beta;

  public Gotoh() {
    super.parameters.add(new AlgorithmParameter("UseBlosumScoringMatrix",
        "Unchecking this option will cause the program to use the pam250"
            + "scoring matrx.", Boolean.class, new Boolean(false)));

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
                + "\n>s2\nRDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT\n")));

    super.parameters.add(new AlgorithmParameter("Gap Opening Penalty",
        "This value will be SUBTRACTED during the scoring process"
            + " every time a gap is inserted into an alignment.",
        Integer.class, 11));

    super.parameters.add(new AlgorithmParameter("UseCompleteTraceBack",
        "This determines the number and sort of alignments returned."
            + " Uncheck to return a random optimal alignment.", Boolean.class,
        new Boolean(true)));

    super.parameters.add(new AlgorithmParameter("Gape Extension Penalty",
        "This value will be SUBTRACTED during the scoring process"
            + " every time a gap is extended in an alignment.", Integer.class,
        1));
  }

  @Override
  public String getName() {
    return new String("Needleman-Wunsch Algorithm");
  }

  @Override
  public String getDescription() {
    return new String(
        "This class accepts any number of sequences in the"
            + " FASTA format and uses the Gotoh"
            + " Algorithm to find the global optimal pairwise"
            + " alignment. To do so, it uses a user-selected"
            + " character alignment cost function that is either PAM250 or BLOSUM62 and"
            + " and affine gap insertion/extension function, which requires defining"
            + " beta and k.");
  }

  @Override
  public Vector<AlgorithmParameter> getInputParameters() {
    return super.parameters;
  }

  @Override
  public String run(Vector<AlgorithmParameter> params) {

    if (params.size() != 5) {
      return "This class expects 5 parameters, <ScoringFunctionName>, "
          + " <StringList of Sequences> <Integer - Gap Penalty> "
          + "<Traceback Type> <Integer - Beta>";
    }

    try {
      // Parse scoring function from first argument
      this.sf = this.parseScoringFunction(params.elementAt(0));
      // Parse sequences from second argument
      this.parseSequences(params.elementAt(1));
      // Parse the gap penalty
      this.gapPenalty = this.parseGapPenalty(params.elementAt(2));
      // Parse the traceback type -- complete, or randomOptimal
      this.tbType = this.parseTracebackType(params.elementAt(3));
      // Parse the value for beta
      this.parseBeta(params.elementAt(4));
    } catch (BioInfException e) {
      return e.getMessage().toString();
    }

    // Perform Alignment and fill out traceback matrix
    initializeMatrices();
    // Recursively fill in remaining cells in all matrices
    try {
      filloutMatrices();
    } catch (BioInfException e) {
      e.printStackTrace();
    }

    // Perform Complete Traceback and use that to build output
    String output = buildOutputString(findAlignments());
    output += "Maximum Score: " + this.m[m.length - 1][m[0].length - 1];
    return output;
  }

  /**
   * Captures the integer value from @param beta and sets interval state var
   * @param beta An Bio.AlgorithmParameter of type Integer
   * @throws BioInfException
   */
  private void parseBeta(AlgorithmParameter beta) throws BioInfException {
    if (beta.type == Integer.class && beta.data.getClass() == Integer.class) {
      this.beta = (Integer) beta.data;
    } else {
      throw new BioInfException("Expected Integer", " but found "
          + beta.data.getClass().toString());
    }
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
    Gotoh myInstance = new Gotoh();

    // run the example the instance with the default parameters
    BioinfAlgorithm.runAlgorithmDefaults(myInstance);
  }

  @Override
  void initializeMatrices() {
    this.m = new int[this.seq1.seq.size() + 1][this.seq2.seq.size() + 1];
    this.p = new int[this.seq1.seq.size() + 1][this.seq2.seq.size() + 1];
    this.q = new int[this.seq1.seq.size() + 1][this.seq2.seq.size() + 1];
    this.t_m =
        new TRACE[this.seq1.seq.size() + 1][this.seq2.seq.size() + 1][3];
    this.t_p =
        new TRACE[this.seq1.seq.size() + 1][this.seq2.seq.size() + 1][3];
    this.t_q =
        new TRACE[this.seq1.seq.size() + 1][this.seq2.seq.size() + 1][3];

    this.m[0][0] = 0;
    this.p[0][0] = 0;
    this.t_p[0][0] = new TRACE[] { TRACE.NULL };
    this.t_m[0][0] = new TRACE[] { TRACE.NULL };
    this.t_q[0][0] = new TRACE[] { TRACE.NULL };

    for (int i = 1; i < this.seq1.seq.size() + 1; ++i) {
      this.p[i][0] = Integer.MIN_VALUE + 1;
      this.m[i][0] = -gapPenalty - (i - 1) * beta;
      this.t_p[i][0] = new TRACE[] { TRACE.M };
      this.t_m[i][0] = new TRACE[] { TRACE.X };
      this.t_q[i][0] = new TRACE[] { TRACE.X }; // we should never reach here
    }
    for (int j = 1; j < seq2.seq.size() + 1; ++j) {
      this.m[0][j] = -gapPenalty - (j - 1) * beta;
      this.q[0][j] = Integer.MIN_VALUE + 1;
      this.t_m[0][j] = new TRACE[] { TRACE.Y };
      this.t_p[0][j] = new TRACE[] { TRACE.Y }; // we should never reach
                                                // here
      this.t_q[0][j] = new TRACE[] { TRACE.M };
    }

  }

  @Override
  void filloutMatrices() throws BioInfException {
    Random r = new Random();
    for (int i = 1; i <= this.seq1.seq.size(); ++i) {
      for (int j = 1; j <= this.seq2.seq.size(); ++j) {
        // Fill in p and t_p
        int[] values =
            new int[] { m[i][j - 1] - gapPenalty, p[i][j - 1] - beta };
        int[] indices = new int[] { TRACE.M.ordinal(), TRACE.Y.ordinal() };
        // Selection sort on BOTH arrays
        for (int outer = 0; outer < values.length - 1; outer++) {
          for (int inner = outer + 1; inner < values.length; inner++) {
            if (values[outer] <= values[inner]) {
              // Exchange elements in BOTH arrays
              int temp = values[outer];
              values[outer] = values[inner];
              values[inner] = temp;
              temp = indices[outer];
              indices[outer] = indices[inner];
              indices[inner] = temp;
            }
          }
        }
        p[i][j] = values[0];
        // If the two options are the same,
        // either put them all in the traceback
        // or, choose one randomly
        if (values[0] == values[1]) {
          if (tbType == TracebackType.RANDOM) {
            t_p[i][j] = new TRACE[] { TRACE.values()[r.nextInt(2)] };
          } else {
            t_p[i][j] = Arrays.copyOfRange(TRACE.values(), 0, 1);
          }
          // Otherwise, just take the max
        } else {
          t_p[i][j] = new TRACE[] { TRACE.values()[indices[0]] };
        }

        // Fill in q and t_q
        values = new int[] { m[i - 1][j] - gapPenalty, q[i - 1][j] - beta };
        indices = new int[] { TRACE.M.ordinal(), TRACE.X.ordinal() };
        // Selection sort on BOTH arrays
        for (int outer = 0; outer < values.length - 1; outer++) {
          for (int inner = outer + 1; inner < values.length; inner++) {
            if (values[outer] <= values[inner]) {
              // Exchange elements in BOTH arrays
              int temp = values[outer];
              values[outer] = values[inner];
              values[inner] = temp;
              temp = indices[outer];
              indices[outer] = indices[inner];
              indices[inner] = temp;
            }
          }
        }
        q[i][j] = values[0];
        // If the two options are the same,
        // either put them all in the traceback
        // or, choose one randomly
        if (values[0] == values[1]) {
          if (tbType == TracebackType.RANDOM) {
            t_q[i][j] = new TRACE[] { TRACE.values()[r.nextInt(2)] };
          } else {
            t_q[i][j] = Arrays.copyOfRange(TRACE.values(), 0, 1);
          }
          // Otherwise, just take the max
        } else {
          t_q[i][j] = new TRACE[] { TRACE.values()[indices[0]] };
        }

        // Using two arrays, track both the values and indices of all
        // all possible directions
        values =
            new int[] {
                p[i][j],
                q[i][j],
                m[i - 1][j - 1]
                    + this.sf.lookupScore(this.seq1.seq.get(i - 1),
                        seq2.seq.get(j - 1)) };
        indices =
            new int[] { TRACE.P.ordinal(), TRACE.Q.ordinal(),
                TRACE.XY.ordinal() };
        // Selection sort on BOTH arrays
        for (int outer = 0; outer < values.length - 1; outer++) {
          for (int inner = outer + 1; inner < values.length; inner++) {
            if (values[outer] <= values[inner]) {
              // Exchange elements in BOTH arrays
              int temp = values[outer];
              values[outer] = values[inner];
              values[inner] = temp;
              temp = indices[outer];
              indices[outer] = indices[inner];
              indices[inner] = temp;
            }
          }
        }
        m[i][j] = values[0];
        // If all three options are the same,
        // either put them all in the traceback
        // or, choose one randomly
        if (values[0] == values[1] && values[1] == values[2]) {
          if (tbType == TracebackType.RANDOM) {
            t_m[i][j] = new TRACE[] { TRACE.values()[r.nextInt(3)] };
          } else {
            t_m[i][j] = Arrays.copyOfRange(TRACE.values(), 0, 2);
          }
          // If two are equal
          // either take both
          // or pick one randomly
        } else if (values[0] == values[1]) {
          if (tbType == TracebackType.RANDOM) {
            t_m[i][j] = new TRACE[] { TRACE.values()[r.nextInt(2)] };
          } else {
            t_m[i][j] =
                new TRACE[] { TRACE.values()[indices[0]],
                    TRACE.values()[indices[1]] };
          }
          // Otherwise, just take the max
        } else {
          t_m[i][j] = new TRACE[] { TRACE.values()[indices[0]] };
        }
      }
    }

  }

  @Override
  ArrayList<Character[][]> findAlignments() {
    int curX = t_m.length - 1;
    int curY = t_m[0].length - 1;
    ArrayList<Character[][]> paths = new ArrayList<Character[][]>();
    Stack<SearchState> stack = new Stack<SearchState>();
    ArrayList<Character[]> curPath = new ArrayList<Character[]>();
    TRACE[][][] curTracebackMatrix = t_m;
    char curTbMName = 'M';

    // This is basically a depth-first graph traversal
    // The traceback matrix can be thought of as an acyclic direct graph
    boolean notDone = true;
    // while we are not done
    while (notDone) {
      // for every direction in the current traceback matrix
      for (int i = 0; i < curTracebackMatrix[curX][curY].length; i++) {
        // if we are in the [0,0] space
        // ... add the completed path to the list of paths
        if (curTracebackMatrix[curX][curY][i] == TRACE.NULL) {
          Character[][] newPath =
              curPath.toArray(new Character[2][curPath.size()]);
          paths.add(newPath);
          // if the current cell has a NORTH arrow
          // we must be in the t_q
          // and the option points to t_q
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are a gap in seq1 and a character in seq2
        } else if (curTracebackMatrix[curX][curY][i] == TRACE.Y) {
          int newY = curY - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { Constants.GAP_SYMBOL, seq2.seq.get(curY - 1) };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(curX, newY, 'P', newPath);
          stack.push(newJunction);
          // if the current cell has a WEST arrow
          // we must be in the t_p
          // and the option points to t_p
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are a gap in seq2 and a character in seq1
        } else if (curTracebackMatrix[curX][curY][i] == TRACE.X) {
          int newX = curX - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { this.seq1.seq.get(curX - 1),
                  Constants.GAP_SYMBOL };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(newX, curY, 'Q', newPath);
          stack.push(newJunction);
          // if the current cell has a NORTHWEST arrow
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are characters from both sequences
        } else if (curTracebackMatrix[curX][curY][i] == TRACE.XY) {
          int newX = curX - 1;
          int newY = curY - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { this.seq1.seq.get(curX - 1),
                  seq2.seq.get(curY - 1) };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(newX, newY, 'M', newPath);

          stack.push(newJunction);
          // if the current cell has an M
          // we could be in P or Q
          // and the option points to t_m
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are characters from both sequences
        } else if (curTracebackMatrix[curX][curY][i] == TRACE.M) {
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment = new Character[0];
          int newX = curX;
          int newY = curY;
          if (curTbMName == 'P') { // we are in t_p
            newY = curY - 1;
            newAlignment =
                new Character[] { Constants.GAP_SYMBOL, seq2.seq.get(curY - 1) };
          } else if (curTbMName == 'Q') { // we are in t_q
            newX = curX - 1;
            newAlignment =
                new Character[] { this.seq1.seq.get(curX - 1),
                    Constants.GAP_SYMBOL };
          }
          newPath.add(0, newAlignment);

          SearchState newJunction = new SearchState(newX, newY, 'M', newPath);
          stack.push(newJunction);
          // if the current cell has an P
          // we must be in t_m
          // and the option points to t_p
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... BUT THERE ARE NO NEWLY ALIGNED CHARACTERS
        } else if (curTracebackMatrix[curX][curY][i] == TRACE.P) {
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          SearchState newJunction = new SearchState(curX, curY, 'P', newPath);
          stack.push(newJunction);
        } else if (curTracebackMatrix[curX][curY][i] == TRACE.Q) {
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          SearchState newJunction = new SearchState(curX, curY, 'Q', newPath);
          stack.push(newJunction);
        }
      }
      // When there are no more SearchStates to explore, it means we have
      // explored the entire graph
      if (stack.isEmpty()) {
        notDone = false;
        // Otherwise, pop off the next SearchState and start exploring
      } else {
        SearchState backtrack = stack.pop();
        curX = backtrack.curX;
        curY = backtrack.curY;
        curPath = backtrack.curPath;
        curTbMName = backtrack.mName;
        if (curTbMName == 'M') {
          curTracebackMatrix = t_m;
        } else if (curTbMName == 'P') {
          curTracebackMatrix = t_p;
        } else if (curTbMName == 'Q') {
          curTracebackMatrix = t_q;
        }
      }
    }
    return paths;
  }

}
