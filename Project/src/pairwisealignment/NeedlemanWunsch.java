package pairwisealignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Stack;
import java.util.Vector;

import utilities.Constants;
import utilities.ProteinSequence;
import exceptions.BioInfException;
import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;

public class NeedlemanWunsch extends PairwiseAlignmentAlgorithm {

  int[][] m;
  TRACE[][][] t;

  public NeedlemanWunsch() {
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
                    + "\n>s2\nRRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD\n")));

    super.parameters.add(new AlgorithmParameter("Gap Penalty",
        "This value will be SUBTRACTED during the scoring process"
            + " every time a gap is inserted into an alignment.",
        Integer.class, 8));

    super.parameters.add(new AlgorithmParameter("UseCompleteTraceBack",
        "This determines the number and sort of alignments returned."
            + " Check to return the full set of optimal alignments.", Boolean.class,
        new Boolean(true)));
  }

  @Override
  public String getName() {
    return new String("Needleman-Wunsch Algorithm");
  }

  @Override
  public String getDescription() {
    return new String("This class accepts 2 sequences in the"
        + " FASTA format and uses the Needleman-Wunsch"
        + " Algorithm to find the global optimal pairwise"
        + " alignment. To do so, it uses a user-selected"
        + " utility function that is either PAM250 or BLOSUM62.");
  }

  @Override
  public Vector<AlgorithmParameter> getInputParameters() {
    return super.parameters;
  }

  @Override
  public String run(Vector<AlgorithmParameter> params) {

    if (params.size() != 4) {
      return "This class expects 4 parameters, <ScoringFunctionName>, "
          + " <StringList of Sequences> <Integer - Gap Penalty> "
          + "<Traceback Type>";
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
    } catch (BioInfException e) {
      return e.getMessage().toString();
    }

    // Initialize Matrices
    initializeMatrices();

    // Perform Alignment and fill out traceback matrix
    try {
      filloutMatrices();
    } catch (BioInfException e) {
      e.printStackTrace();
    }

    // Perform Complete Traceback
    // and use that to build output;
    String output = buildOutputString(findAlignments());
    output += "Maximum Score: " + this.m[m.length - 1][m[0].length - 1];
    return output;
  }

  /**
   * Unlike the normal setup routine, this function expects only 2 arguments --
   * a type for the scoring function, and a gap penalty. Sequences are defined
   * manually! TracebackType is always set to Random.
   * 
   * @param params
   * @throws WrongParameterTypeException
   * @throws WrongParameterValueException
   */
  public void internal_setup(Vector<AlgorithmParameter> params)
    throws BioInfException {
    // Parse scoring function from first argument
    this.sf = this.parseScoringFunction(params.elementAt(0));
    // Parse the gap penalty
    this.gapPenalty = this.parseGapPenalty(params.elementAt(1));
    // We only want to look at ONE top alignment
    this.tbType = TracebackType.RANDOM;
  }

  /**
   * Unlike run(), this function returns the optimal alignment It also assumes
   * that setup_no_output() has been called previously
   * 
   * @param params
   * @return PairwiseAlignment an optimal pairwise alignment
   */
  public PairwiseAlignment internal_run(ProteinSequence seq1,
    ProteinSequence seq2) {

    this.seq1 = seq1;
    this.seq2 = seq2;

    PairwiseAlignment pa = new PairwiseAlignment();

    // Initialize Matrices
    initializeMatrices();

    // Perform Alignment and fill out traceback matrix
    try {
      filloutMatrices();
    } catch (BioInfException e) {
      e.printStackTrace();
    }

    Character[][] arr = findAlignments().get(0);

    // Array really should be transposed. We want arr[][0]
    ArrayList<Character> a1 = new ArrayList<Character>();
    ArrayList<Character> a2 = new ArrayList<Character>();
    for (int charNum = 0; charNum < arr.length; charNum++) {
      a1.add(arr[charNum][0]);
      a2.add(arr[charNum][1]);
    }
    pa.a1 = a1.toArray(new Character[0]);
    pa.a2 = a2.toArray(new Character[0]);
    pa.id1 = this.seq1.id;
    pa.id2 = this.seq2.id;

    return pa;
  }

  @Override
  void initializeMatrices() {
    this.m = new int[this.seq1.seq.size() + 1][this.seq2.seq.size() + 1];
    this.t = new TRACE[this.seq1.seq.size() + 1][this.seq2.seq.size() + 1][3];
    // Initialize m and t
    for (int i = 0; i < this.seq1.seq.size() + 1; ++i) {
      this.m[i][0] = i * this.gapPenalty * -1;
      this.t[i][0] = new TRACE[] { TRACE.X };
    }
    for (int j = 0; j < this.seq2.seq.size() + 1; ++j) {
      this.m[0][j] = j * this.gapPenalty * -1;
      if (j == 0) {
        this.t[0][j] = new TRACE[] { TRACE.NULL };
      } else {
        this.t[0][j] = new TRACE[] { TRACE.Y };
      }
    }

  }

  @Override
  void filloutMatrices() throws BioInfException {
    Random r = new Random();
    for (int i = 1; i <= this.seq1.seq.size(); ++i) {
      for (int j = 1; j <= this.seq2.seq.size(); ++j) {
        // Using two arrays, track both the values and indices of all
        // all possible directions
        int[] values =
            new int[] {
                m[i - 1][j] - gapPenalty,
                m[i][j - 1] - gapPenalty,
                m[i - 1][j - 1]
                    + this.sf.lookupScore(this.seq1.seq.get(i - 1),
                        this.seq2.seq.get(j - 1)) };
        int[] indices =
            new int[] { TRACE.X.ordinal(), TRACE.Y.ordinal(),
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
        // TODO: CHECK THIS
        if (values[0] == values[1] && values[1] == values[2]) {
          if (tbType == TracebackType.RANDOM) {
            t[i][j] = new TRACE[] { TRACE.values()[r.nextInt(3)] };
          } else {
            t[i][j] = Arrays.copyOfRange(TRACE.values(), 0, 2);
          }
          // If two are equal
          // either take both
          // or pick one randomly
        } else if (values[0] == values[1]) {
          if (tbType == TracebackType.RANDOM) {
            // TODO: shouldn't this refer to indices? check nw3d
            t[i][j] = new TRACE[] { TRACE.values()[indices[r.nextInt(2)]] };
          } else {
            t[i][j] =
                new TRACE[] { TRACE.values()[indices[0]],
                    TRACE.values()[indices[1]] };
          }
          // Otherwise, just take the max
        } else {
          t[i][j] = new TRACE[] { TRACE.values()[indices[0]] };
        }
      }
    }
  }

  @Override
  ArrayList<Character[][]> findAlignments() {
    // Since the traceback is built with only one value in each cell when the
    // tracebackType is Random, we can simply find all the paths.
    int curX = t.length - 1;
    int curY = t[0].length - 1;
    ArrayList<Character[][]> paths = new ArrayList<Character[][]>();
    Stack<SearchState> stack = new Stack<SearchState>();
    ArrayList<Character[]> curPath = new ArrayList<Character[]>();

    // This is basically a depth-first graph traversal
    // The traceback matrix can be thought of as an acyclic direct graph
    boolean notDone = true;
    // while we are not done
    while (notDone) {
      for (int i = 0; i < t[curX][curY].length; i++) {
        // if we are in the [0,0] space
        // ... add the completed path to the list of paths
        if (t[curX][curY][i] == TRACE.NULL) {
          Character[][] newPath =
              curPath.toArray(new Character[curPath.size()][2]);
          paths.add(newPath);
          // if the current cell has a NORTH arrow
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are a gap in seq1 and a character in seq2
        } else if (t[curX][curY][i] == TRACE.Y) {
          int newY = curY - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { Constants.GAP_SYMBOL,
                  this.seq2.seq.get(curY - 1) };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(curX, newY, newPath);
          stack.push(newJunction);
          // if the current cell has a WEST arrow
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are a gap in seq2 and a character in seq1
        } else if (t[curX][curY][i] == TRACE.X) {
          int newX = curX - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { this.seq1.seq.get(curX - 1),
                  Constants.GAP_SYMBOL };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(newX, curY, newPath);
          stack.push(newJunction);
          // if the current cell has a NORTHWEST arrow
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are characters from both sequences
        } else if (t[curX][curY][i] == TRACE.XY) {
          int newX = curX - 1;
          int newY = curY - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { this.seq1.seq.get(curX - 1),
                  this.seq2.seq.get(curY - 1) };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(newX, newY, newPath);
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
      }
    }
    return paths;
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
    NeedlemanWunsch myInstance = new NeedlemanWunsch();

    // run the example the instance with the default parameters
    BioinfAlgorithm.runAlgorithmDefaults(myInstance);
  }

}
