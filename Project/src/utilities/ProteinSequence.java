package utilities;

import java.util.ArrayList;
import java.util.Arrays;
import exceptions.BioInfException;

public class ProteinSequence {

  public String id;
  public ArrayList<Character> seq;

  public ProteinSequence(String inID, String proteinString)
    throws BioInfException {
    id = inID;
    seq = new ArrayList<Character>();
    for (int i = 0; i < proteinString.length(); i++) {
      seq.add(proteinString.charAt(i));
    }
  }

  public ProteinSequence(String proteinString)
    throws BioInfException {
    seq = new ArrayList<Character>();
    for (int i = 0; i < proteinString.length(); i++) {
      seq.add(proteinString.charAt(i));
    }
  }

  public ProteinSequence() {
    seq = new ArrayList<Character>();
  }

  /**
   * Used to add any valid character other than '*' to a ProteinSequence.
   * @param c
   * @throws BioInfException
   */
  public void add(char c) throws BioInfException {
    // The -1 here makes it so the last elt of Constants.SIGMA, the '*'
    // character, is not considered.
    for (int cNum = 0; cNum < Constants.SIGMA.length - 1; cNum++) {
      if (Constants.SIGMA[cNum] == c) {
        seq.add(c);
        return;
      }
    }
    throw new BioInfException("Found inappropriate character: " + c);
  }
  
  public void addStar() {
    seq.add(Constants.SIGMA[20]);
  }

  @Override
  public String toString() {
    StringBuilder result = new StringBuilder("[");
    if (this.seq.size() > 0) {
      result.append(this.seq.get(0));
    }
    for (int i = 1; i < this.seq.size(); ++i) {
      result.append(", " + this.seq.get(i));
    }
    result.append("]");
    return result.toString();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof ProteinSequence)) {
      return false;
    }
    ProteinSequence ps = (ProteinSequence) obj;
    return Arrays.equals(this.seq.toArray(), ps.seq.toArray());
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  public static ProteinSequence[] parse(String[] lines, int numPS)
    throws BioInfException {
    ProteinSequence[] ps = new ProteinSequence[numPS];
    for (int lineNum = 0; lineNum / 2 < numPS && lineNum < lines.length; lineNum +=
        2) {
      if (lines[lineNum].charAt(0) != '>') {
        throw new BioInfException("Improper FASTA format");
      }
      String identifier = lines[lineNum].substring(1);
      ps[lineNum / 2] = new ProteinSequence(identifier, lines[lineNum + 1]);
    }

    return ps;
  }

}
