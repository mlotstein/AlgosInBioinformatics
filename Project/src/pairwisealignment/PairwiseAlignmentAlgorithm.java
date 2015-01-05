package pairwisealignment;

import exceptions.BioInfException;
import gui.AlgorithmParameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;
import utilities.Constants;
import utilities.ProteinSequence;

public abstract class PairwiseAlignmentAlgorithm extends gui.BioinfAlgorithm {

  /**
   * A function for evaluating pairs of characters.
   */
  public ScoringFunction sf;
  public ProteinSequence seq1;
  public ProteinSequence seq2;

  TracebackType tbType;
  public int gapPenalty;

  static enum TracebackType {
    RANDOM, COMPLETE;
  };

  static enum TRACE {
    X, Y, XY, P, Q, M, NULL;

    char[] toString = new char[] { '←', '↑', '↖', 'P', 'Q', 'M', '∅' };

    @Override
    public String toString() {
      return "" + this.toString[this.ordinal()];
    }

  };

  /**
   * A class used during traceback that logs the state of the search.
   * 
   * @author Max
   * 
   */
  protected class SearchState {

    public int curX;
    public int curY;
    public ArrayList<Character[]> curPath;
    public char mName;

    public SearchState(int curX, int curY, ArrayList<Character[]> curPath) {
      this.curX = curX;
      this.curY = curY;
      this.curPath = curPath;
    }

    public SearchState(int curX, int curY, char mName,
        ArrayList<Character[]> curPath) {
      this.curX = curX;
      this.curY = curY;
      this.mName = mName;
      this.curPath = curPath;
    }
  }

  @Override
  public abstract String getName();

  @Override
  public abstract String getDescription();

  @Override
  public Vector<AlgorithmParameter> getInputParameters() {
    return super.parameters;
  }

  @Override
  public abstract String run(Vector<AlgorithmParameter> params);

  /**
   * Parses the information necessary to initialize a ScoringFunction from the
   * input parameter
   * 
   * @param tbType
   *          , a string with value 'Random' or 'Complete'
   * @return TracebackType
   * @throws BioInfException
   */
  public TracebackType parseTracebackType(AlgorithmParameter tbType)
    throws BioInfException {
    if (!(tbType.type == Boolean.class && tbType.data.getClass() == Boolean.class)) {
      throw new BioInfException("Expected a Boolean", "but found "
          + tbType.data.getClass().toString());
    }
    if ((Boolean) tbType.data) {
      return TracebackType.COMPLETE;
    } else {
      return TracebackType.RANDOM;
    }
  }

  /**
   * Parses the gapPenalty from the AlgorithmParameter input
   * 
   * @param gp
   * @return Integer
   * @throws BioInfException
   */
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

  /**
   * Parses exactly two sequences in FASTA format from StringList
   * AlgorithmParameter
   * 
   * @param param
   * @throws BioInfException
   *           when input contains inappropriate characters or is
   *           inappropriately formatted.
   */
  public void parseSequences(AlgorithmParameter param) throws BioInfException {
    if (param.type != gui.StringList.class
        || param.data.getClass() != gui.StringList.class) {
      throw new BioInfException("Expected StringList", "but found "
          + param.type.toString());
    }
    String[] lines = param.data.toString().trim().split("\n");
    ProteinSequence[] ps = ProteinSequence.parse(lines, 2);
    this.seq1 = ps[0];
    this.seq2 = ps[1];
  }

  /**
   * Uses a String to initialize a ScoringFunction and return it
   * @param p A String with value 'blosum62' or 'pam250'
   * @return an instantiation of ScoringFunction with type set as appropriate
   * @throws BioInfException
   */
  public ScoringFunction parseScoringFunction(AlgorithmParameter p)
    throws BioInfException {
    if (!(p.type == Boolean.class || p.data.getClass() != Boolean.class)) {
      throw new BioInfException("Expected type Boolean", "but found "
          + p.type.toString());
    }

    return new ScoringFunction((Boolean) p.data);
  }

  abstract void initializeMatrices();

  abstract void filloutMatrices() throws BioInfException;

  abstract ArrayList<Character[][]> findAlignments();

  private String buildParameterOutputString() {
    StringBuilder s = new StringBuilder();
    for (AlgorithmParameter ap : this.parameters) {
      s.append(ap.name + ": " + ap.data + "\n");
    }
    s.append("\n");
    return s.toString();
  }

  /**
   * Accepts a bunch of paths and prepares a single string for output
   * @param paths See comments in function
   * @return outputString formatted as specified in pdf
   */
  public String buildOutputString(ArrayList<Character[][]> paths) {
    // Paths are are basically N x 2 arrays of characters, such as:
    // Path 1: { {A, T} {G, _} {_, A}...}
    // Path 2: { {A, _},{G, T},{_, A}...}
    // We are interested in a 'rows', so we need to iterate to find them.
    // So, this converts from above to:
    // { {A G _...}, {T _ A...}}
    // { {A G _ ..}, {_ T A...}}
    // Which are 2 DIFFERENT alignments
    StringBuilder output = new StringBuilder();
    output.append(buildParameterOutputString());
    
    String[][] alignments = new String[paths.size()][2];
    int pathNum = 0;
    for (Character[][] path : paths) {
      alignments[pathNum][0] = "";
      alignments[pathNum][1] = "";
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
      output.append("Alignment Number " + (alignment + 1) + "\n");

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
          new char[Math.max(this.seq1.id.length(), this.seq2.id.length()) + 1];
      char[] seq1Pad = new char[pad.length - this.seq1.id.length()];
      Arrays.fill(seq1Pad, ' ');
      char[] matchPad = new char[pad.length];
      Arrays.fill(matchPad, ' ');
      char[] seq2Pad = new char[pad.length - this.seq2.id.length()];
      Arrays.fill(seq2Pad, ' ');
      String[] topPieces =
          splitIntoChunks(alignments[alignment][0], (80 - pad.length));
      // alignments[alignment][0].split(".{1,80}");
      String[] bottomPieces =
          splitIntoChunks(alignments[alignment][1], (80 - pad.length));
      String[] matchPieces =
          splitIntoChunks(
              getAlignmentString((String) alignments[alignment][0],
                  (String) alignments[alignment][1]), (80 - pad.length));

      for (int pieceNum = 0; pieceNum < topPieces.length; pieceNum++) {
        output.append(this.seq1.id + String.valueOf(seq1Pad)
            + topPieces[pieceNum] + "\n");
        output.append(String.valueOf(matchPad) + matchPieces[pieceNum] + "\n");
        output.append(this.seq2.id + String.valueOf(seq2Pad)
            + bottomPieces[pieceNum] + "\n");
      }
      output.append("\n");
    }
    return output.toString();
  }

  /**
   * Given a String s, splits into an Array of String, where each elt
   * except the last is equal to chunkSize and the last elt is the
   * remaining characters.
   * @param s
   * @param chunkSize
   * @return
   */
  private String[] splitIntoChunks(String s, int chunkSize) {
    ArrayList<String> as = new ArrayList<String>();
    int begin = 0;
    while (begin < s.length()) {
      as.add(s.substring(begin, Math.min(s.length(), begin + chunkSize)));
      begin += chunkSize + 1;
    }
    return as.toArray(new String[0]);
  }

  /**
   * As per spec pdf, when one of the characters is a Constants.GAP_SYMBOL,
   * return a ' ', when characters match in the alignment, show a '*',
   * otherwise return a ':'
   * @param c1 character
   * @param c2 character
   * @return either ' ','*' or ':'
   */
  protected static char getAlignmentMatchChar(char c1, char c2) {
    if (c1 == Constants.GAP_SYMBOL || c2 == Constants.GAP_SYMBOL) {
      return ' ';
    }
    return (c1 == c2) ? '*' : ':';
  }

  /**
   * Given an alignment, append the result of calling getAlignmentMatchChar
   * on every pair of chars with the same index.
   * @param s1
   * @param s2
   * @return A string comprised of returned values from getAlignmentMatchChar
   */
  protected static String getAlignmentString(String s1, String s2) {
    StringBuilder matchString = new StringBuilder();
    for (int charNum = 0; charNum < s1.length() && charNum < s2.length(); charNum++) {
      matchString.append(getAlignmentMatchChar(s1.charAt(charNum),
          s2.charAt(charNum)));
    }
    return matchString.toString();
  }

  /**
   * Clones the values in an ArrayList
   * @param al
   * @return cloned ArrayList.
   */
  protected static ArrayList<Character[]> deepCloneList(
    ArrayList<Character[]> al) {
    ArrayList<Character[]> clone = new ArrayList<Character[]>(al.size());
    for (Character[] a : al) {
      clone.add(a.clone());
    }
    return clone;
  }

}
