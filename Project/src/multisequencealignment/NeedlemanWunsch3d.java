package multisequencealignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Stack;
import java.util.Vector;

import exceptions.BioInfException;

import pairwisealignment.ScoringFunction;
import utilities.Constants;
import utilities.ProteinSequence;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;

public class NeedlemanWunsch3d extends BioinfAlgorithm {

  public ScoringFunction sf;
  public ProteinSequence seqx;
  public ProteinSequence seqy;
  public ProteinSequence seqz;
  public int[][][] m;
  public TRACE[][][][] t;

  TracebackType tbType;
  public int gapPenalty;

  static enum TracebackType {
    RANDOM, COMPLETE;
  };

  static enum TRACE {
    X, Y, Z, XY, XZ, YZ, XYZ, NULL;

    String[] toString = new String[] { "X", "Y", "Z", "XY", "XZ", "YZ", "XYZ",
        "∅" };

    @Override
    public String toString() {
      return this.toString[this.ordinal()];
    }

  };

  protected class SearchState {

    public int curX;
    public int curY;
    public int curZ;
    public ArrayList<Character[]> curPath;
    public char mName;

    public SearchState(int curX, int curY, int curZ,
        ArrayList<Character[]> curPath) {
      this.curX = curX;
      this.curY = curY;
      this.curZ = curZ;
      this.curPath = curPath;
    }

    public SearchState(int curX, int curY, int curZ, char mName,
        ArrayList<Character[]> curPath) {
      this.curX = curX;
      this.curY = curY;
      this.curZ = curZ;
      this.mName = mName;
      this.curPath = curPath;
    }
  }

  public NeedlemanWunsch3d() {
    super.parameters.add(new AlgorithmParameter("UseBlosumScoringMatrix",
        "Unchecking this option will cause the program to use the pam250"
            + "scoring matrx.", Boolean.class, new Boolean(true)));

    super.parameters.add(new AlgorithmParameter("Sequence Data",
        "Any number of correctly formatted FASTA sequences."
            + " The correct format for each is the sequence is:"
            + " '><sequence name><newline><sequence data><newline>"
            + " and can be repeated consecutively any number of times."
            + " Beware though: this algorithm requires O(m^n) time and"
            + " space to complete, where m = the average sequence "
            + " length and n = the number of sequences.", StringList.class,
        new StringList(">s1\nTTTA\n>s2\nAAAA\n>s3\nAAAT")));

    super.parameters.add(new AlgorithmParameter("Gap Penalty",
        "This value will be SUBTRACTED during the scoring process"
            + " every time a gap is inserted into an alignment.",
        Integer.class, 5));

    super.parameters
        .add(new AlgorithmParameter(
            "Traceback Type",
            "This determines the number and sort of alignments returned."
                + " If 'Random' is used, then the algorithm returns one randomly selected "
                + "optimal alignment. If 'Complete' is used, all optimal "
                + "alignments are returned.", String.class, "Complete"));
  }

  @Override
  public String getName() {
    return "Needleman-Wunsch 3D";
  }

  @Override
  public String getDescription() {
    return new String("This class accepts 3 sequences in the"
        + " FASTA format and uses the Needleman-Wunsch"
        + " Algorithm to find the global optimal "
        + " alignment. To do so, it uses a user-selected"
        + " utility function that is either PAM250 or BLOSUM62.");
  }

  @Override
  public Vector<AlgorithmParameter> getInputParameters() {
    return super.parameters;
  }

  public TracebackType parseTracebackType(AlgorithmParameter tbType)
    throws BioInfException {
    if (!(tbType.type == String.class && tbType.data.getClass() == String.class)) {
      throw new BioInfException("Expected String", "but found "
          + tbType.data.getClass().toString());
    }
    if (!(((String) tbType.data).equals("Random") || ((String) tbType.data)
        .equals("Complete"))) {
      throw new BioInfException("Expedted String value Random or Complete",
          "but found " + (String) tbType.data);
    } else if (((String) tbType.data).equals("Random")) {
      return TracebackType.RANDOM;
    } else {
      return TracebackType.COMPLETE;
    }
  }

  public int parseGapPenalty(AlgorithmParameter gp) throws BioInfException {
    if (!(gp.type == Integer.class && gp.data.getClass() == Integer.class)) {
      throw new BioInfException("Expected Integer", "but found "
          + gp.data.getClass().toString());
    }
    if ((Integer) gp.data < 0) {
      throw new BioInfException("Expected Integer >= 0", "but found"
          + (Integer) gp.data);
    } else {
      return (Integer) gp.data;
    }
  }

  public void parseSequences(AlgorithmParameter param) throws BioInfException {
    if (param.type != gui.StringList.class
        || param.data.getClass() != gui.StringList.class) {
      throw new BioInfException("Expected StringList", "but found "
          + param.type.toString());
    }
    String[] lines = param.data.toString().trim().split("\n");
    if (lines.length != 6) {
      throw new BioInfException("Expected StringList"
          + "precisely 3 sequences in 6 lines", "but found " + lines.length);
    }
    ProteinSequence[] ps = ProteinSequence.parse(lines, 3);
    this.seqx = ps[0];
    this.seqy = ps[1];
    this.seqz = ps[2];
  }

  public ScoringFunction parseScoringFunction(AlgorithmParameter p)
    throws BioInfException {
    if (!(p.type == Boolean.class || p.data.getClass() != Boolean.class)) {
      throw new BioInfException("Expected Boolean", "but found "
          + p.type.toString());
    }

    return new ScoringFunction((Boolean) p.data);
  }

  /**
   * Accepts a bunch of paths and prepares a single string for output
   * 
   * @param paths
   * @return outputString
   */
  public String buildOutputString(ArrayList<Character[][]> paths) {
    // Paths are are basically N x 2 arrays of characters, such as:
    // Path 1: { {A, T, G} {G, _, A} {_, A, _}...}
    // Path 2: { {A, _, _},{G, T, _},{_, T, A}...}
    // We are interested in a 'rows', so we need to iterate to find them.
    // So, this converts from above to:
    // { {A G _...}, {T _ A...}, {G, A, _ ...}
    // { {A G _ ..}, {_ T A...}}
    // Which are 2 DIFFERENT alignments
    StringBuilder output = new StringBuilder();
    String[][] alignments = new String[paths.size()][3];
    int pathNum = 0;
    for (Character[][] path : paths) {
      alignments[pathNum][0] = "";
      alignments[pathNum][1] = "";
      alignments[pathNum][2] = "";
      for (int c = 0; c < path.length; c++) {
        for (int aNum = 0; aNum < alignments[0].length; aNum++) {
          alignments[pathNum][aNum] =
              alignments[pathNum][aNum] + path[c][aNum];
        }
      }
      pathNum++;
    }

    // For each alignment of the form String1, String2
    // ... The first line is <String 1 id> <alignment characters>
    // ... The second line is a space buffer
    for (int alignment = 0; alignment < alignments.length; alignment++) {

      if (alignment >= 1) {
        output.append("\n****************\n");
      }

      // Each set of alignments has two strings -- the top and the bottom
      // alignment
      // For each pair, we first find the characters that go in between them
      // Then we split all three of these strings into pieces
      // And layer them together -- top, match, bottom
      // And repeat, if there is another alignment.

      // Because every line has a pad in the front eql to the seqId, we split
      // into chunks
      // of size 80 - seqId.length
      char[] pad =
          new char[Math.max(this.seqx.id.length() + 1,
              Math.max(this.seqy.id.length() + 1, this.seqz.id.length() + 1))];
      char[] seq1Pad = new char[pad.length - this.seqx.id.length()];
      Arrays.fill(seq1Pad, ' ');
      char[] seq2Pad = new char[pad.length - this.seqy.id.length()];
      Arrays.fill(seq2Pad, ' ');
      char[] seq3Pad = new char[pad.length - this.seqz.id.length()];
      Arrays.fill(seq3Pad, ' ');
      String[] topPieces =
          splitIntoChunks(alignments[alignment][0], (80 - pad.length));
      // alignments[alignment][0].split(".{1,80}");
      String[] middlePieces =
          splitIntoChunks(alignments[alignment][1], (80 - pad.length));
      String[] bottomPieces =
          splitIntoChunks(alignments[alignment][2], (80 - pad.length));

      for (int pieceNum = 0; pieceNum < topPieces.length; pieceNum++) {
        output.append(this.seqx.id + String.valueOf(seq1Pad)
            + topPieces[pieceNum] + "\n");
        output.append(this.seqy.id + String.valueOf(seq2Pad)
            + middlePieces[pieceNum] + "\n");
        output.append(this.seqz.id + String.valueOf(seq3Pad)
            + bottomPieces[pieceNum] + "\n");
      }
    }
    return output.toString();
  }

  private String[] splitIntoChunks(String s, int chunkSize) {
    ArrayList<String> as = new ArrayList<String>();
    int begin = 0;
    while (begin < s.length()) {
      as.add(s.substring(begin, Math.min(s.length(), begin + chunkSize)));
      begin += chunkSize + 1;
    }
    return as.toArray(new String[0]);
  }

  protected static char getAlignmentMatchChar(char c1, char c2) {
    if (c1 == Constants.GAP_SYMBOL || c2 == Constants.GAP_SYMBOL) {
      return ' ';
    }
    return (c1 == c2) ? '*' : ':';
  }

  protected static String getAlignmentString(String s1, String s2) {
    StringBuilder matchString = new StringBuilder();
    for (int charNum = 0; charNum < s1.length() && charNum < s2.length(); charNum++) {
      matchString.append(getAlignmentMatchChar(s1.charAt(charNum),
          s2.charAt(charNum)));
    }
    return matchString.toString();
  }

  protected static ArrayList<Character[]> deepCloneList(
    ArrayList<Character[]> al) {
    ArrayList<Character[]> clone = new ArrayList<Character[]>(al.size());
    for (Character[] a : al) {
      clone.add(a.clone());
    }
    return clone;
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
    return buildOutputString(findAlignments());
  }

  void initializeMatrices() {
    this.m =
        new int[this.seqx.seq.size() + 1][this.seqy.seq.size() + 1][this.seqz.seq
            .size() + 1];
    this.t =
        new TRACE[this.seqx.seq.size() + 1][this.seqy.seq.size() + 1][this.seqz.seq
            .size() + 1][7];
    // Initialize m and t
    for (int x = 0; x < this.seqx.seq.size() + 1; ++x) {
      this.m[x][0][0] = x;
      this.t[x][0][0] = new TRACE[] { TRACE.X };
    }
    for (int y = 0; y < this.seqy.seq.size() + 1; ++y) {
      this.m[0][y][0] = y;
      if (y == 0) {
        this.t[0][y][0] = new TRACE[] { TRACE.NULL };
      } else {
        this.t[0][y][0] = new TRACE[] { TRACE.Y };
      }
    }

    for (int z = 0; z < this.seqz.seq.size() + 1; ++z) {
      this.m[0][0][z] = z;
      if (z == 0) {
        this.t[0][0][z] = new TRACE[] { TRACE.NULL };
      } else {
        this.t[0][0][z] = new TRACE[] { TRACE.Z };
      }
    }

  }

  void filloutMatrices() throws BioInfException {
    Random r = new Random();
    for (int x = 1; x <= this.seqx.seq.size(); ++x) {
      for (int y = 1; y <= this.seqy.seq.size(); ++y) {
        for (int z = 1; z <= this.seqz.seq.size(); ++z) {
          // Using two arrays, track both the values and indices of all
          // all possible directions
          int[] values =
              new int[] {
                  m[x - 1][y][z] - gapPenalty * 2,
                  m[x][y - 1][z] - gapPenalty * 2,
                  m[x][y][z - 1] - gapPenalty * 2,
                  m[x - 1][y - 1][z]
                      - gapPenalty
                      + this.sf.lookupScore(this.seqx.seq.get(x - 1),
                          this.seqy.seq.get(y - 1)),
                  m[x - 1][y][z - 1]
                      - gapPenalty
                      + this.sf.lookupScore(this.seqx.seq.get(x - 1),
                          this.seqz.seq.get(z - 1)),
                  m[x][y - 1][z - 1]
                      - gapPenalty
                      + this.sf.lookupScore(this.seqy.seq.get(y - 1),
                          this.seqz.seq.get(z - 1)),
                  m[x - 1][y - 1][z - 1]
                      + this.sf.lookupScore(this.seqx.seq.get(x - 1),
                          this.seqy.seq.get(y - 1))
                      + this.sf.lookupScore(this.seqx.seq.get(x - 1),
                          this.seqz.seq.get(z - 1))
                      + this.sf.lookupScore(this.seqy.seq.get(y - 1),
                          this.seqz.seq.get(z - 1)) };
          int[] indices =
              new int[] { TRACE.X.ordinal(), TRACE.Y.ordinal(),
                  TRACE.Z.ordinal(), TRACE.XY.ordinal(), TRACE.XZ.ordinal(),
                  TRACE.YZ.ordinal(), TRACE.XYZ.ordinal() };
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
          m[x][y][z] = values[0];
          // Count how many are equivalent to the top alignment
          int firstSubOpt = values.length;
          for (int vNum = 1; firstSubOpt == values.length
              && vNum < values.length; vNum++) {
            if (values[vNum] != values[0]) {
              firstSubOpt = vNum;
            }
          }
          if (tbType == TracebackType.RANDOM) {
            t[x][y][z] =
                new TRACE[] { TRACE.values()[indices[r.nextInt(firstSubOpt)]] };
          } else {
            TRACE[] temp = new TRACE[firstSubOpt];
            for (int optNum = 0; optNum < firstSubOpt; optNum++) {
              temp[optNum] = TRACE.values()[indices[optNum]];
            }
            t[x][y][z] = temp;
          }
        }
      }
    }
  }

  ArrayList<Character[][]> findAlignments() {
    // Since the traceback is built with only one value in each cell when the
    // tracebackType is Random, we can simply find all the paths.
    int curX = t.length - 1;
    int curY = t[0].length - 1;
    int curZ = t[0][0].length - 1;
    ArrayList<Character[][]> paths = new ArrayList<Character[][]>();
    Stack<SearchState> stack = new Stack<SearchState>();
    ArrayList<Character[]> curPath = new ArrayList<Character[]>();

    // This is basically a depth-first graph traversal
    // The traceback matrix can be thought of as an acyclic direct graph
    boolean notDone = true;
    // while we are not done
    while (notDone) {
      for (int i = 0; i < t[curX][curY][curZ].length; i++) {
        // if we are in the [0,0] space
        // ... add the completed path to the list of paths
        if (t[curX][curY][curZ][i] == TRACE.NULL) {
          Character[][] newPath =
              curPath.toArray(new Character[2][curPath.size()]);
          paths.add(newPath);
          // if the current cell has a NORTH arrow
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are a gap in seq1 and a character in seq2
        } else if (t[curX][curY][curZ][i] == TRACE.Y) {
          int newY = curY - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { Constants.GAP_SYMBOL,
                  this.seqy.seq.get(curY - 1), Constants.GAP_SYMBOL };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(curX, newY, curZ, newPath);
          stack.push(newJunction);
          // if the current cell has a WEST arrow
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are a gap in seq2 and a character in seq1
        } else if (t[curX][curY][curZ][i] == TRACE.X) {
          int newX = curX - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { this.seqx.seq.get(curX - 1),
                  Constants.GAP_SYMBOL, Constants.GAP_SYMBOL };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(newX, curY, curZ, newPath);
          stack.push(newJunction);
          // if the current cell has a NORTHWEST arrow
          // ... build a new SearchState and put it on the stack
          // ... ... this SearchState will consist of the current path
          // ... ... plus the newly aligned characters -- which in this case
          // ... ... are characters from both sequences
        } else if (t[curX][curY][curZ][i] == TRACE.Z) {
          int newZ = curZ - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { Constants.GAP_SYMBOL, Constants.GAP_SYMBOL,
                  this.seqz.seq.get(curZ - 1) };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(curX, curY, newZ, newPath);
          stack.push(newJunction);
        } else if (t[curX][curY][curZ][i] == TRACE.XY) {
          int newX = curX - 1;
          int newY = curY - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { this.seqx.seq.get(curX - 1),
                  this.seqy.seq.get(curY - 1), Constants.GAP_SYMBOL };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(newX, newY, curZ, newPath);
          stack.push(newJunction);
        } else if (t[curX][curY][curZ][i] == TRACE.XZ) {
          int newX = curX - 1;
          int newZ = curZ - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { this.seqx.seq.get(curX - 1),
                  Constants.GAP_SYMBOL, this.seqz.seq.get(curZ - 1) };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(newX, curY, newZ, newPath);
          stack.push(newJunction);
        } else if (t[curX][curY][curZ][i] == TRACE.YZ) {
          int newY = curY - 1;
          int newZ = curZ - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { Constants.GAP_SYMBOL,
                  this.seqy.seq.get(curY - 1), this.seqz.seq.get(curZ - 1) };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(curX, newY, newZ, newPath);
          stack.push(newJunction);
        } else if (t[curX][curY][curZ][i] == TRACE.XYZ) {
          int newX = curX - 1;
          int newY = curY - 1;
          int newZ = curZ - 1;
          ArrayList<Character[]> newPath = deepCloneList(curPath);
          Character[] newAlignment =
              new Character[] { this.seqx.seq.get(curX - 1),
                  this.seqy.seq.get(curY - 1), this.seqz.seq.get(curZ - 1) };
          newPath.add(0, newAlignment);
          SearchState newJunction = new SearchState(newX, newY, newZ, newPath);
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
        curZ = backtrack.curZ;
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
    NeedlemanWunsch3d myInstance = new NeedlemanWunsch3d();

    // run the example the instance with the default parameters
    BioinfAlgorithm.runAlgorithmDefaults(myInstance);
  }

}
