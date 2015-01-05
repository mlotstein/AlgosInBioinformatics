package clustering;

import java.util.Arrays;

import exceptions.BioInfException;
import pairwisealignment.ScoringFunction;

import utilities.Constants;

public class DistanceFunction {

  ScoringFunction sf;
  int gapPenalty;

  public DistanceFunction(ScoringFunction sf, int gapPenalty) {
    this.sf = sf;
    this.gapPenalty = gapPenalty;
  }

  /**
   *  Note: this function assumes that input arrays constitute a valid alignment
   * and thus have the same number of characters and consist only of valid
   * characters.
   * @param s1a
   * @param s2a
   * @return double the distance between them
   * @throws BioInfException
   */
  public double computeDistance(Character[] s1a, Character[] s2a) throws BioInfException  {
    double Sobs = 0.0, Srand = 0.0, Smax = 0.0;

    if (s1a.length != s2a.length) {
      throw new BioInfException("Alignments of difference length!");
    }

    // Compute Sobs, the sum of the observed score of the alignment
    // using the scoring function
    // And Smax, the average of each seq. aligned with itself
    for (int cNum = 0; cNum < s1a.length; cNum++) {
      // If the alignment contains a gap
      if (s1a[cNum] == Constants.GAP_SYMBOL || s2a[cNum] == Constants.GAP_SYMBOL) {
        Sobs += (double) -1 *  this.gapPenalty;
      } else {
        Sobs += (double) this.sf.lookupScore(s1a[cNum], s2a[cNum]);
      }
      if (s1a[cNum] != Constants.GAP_SYMBOL) {
        Smax += (double) this.sf.lookupScore(s1a[cNum], s1a[cNum]);
      }   
      if (s2a[cNum] != Constants.GAP_SYMBOL) {
        Smax += (double) this.sf.lookupScore(s2a[cNum], s2a[cNum]);
      }
    }
    Smax /= 2.0;

    // Compute Srand
    double L = (double) s1a.length;
    double summation = 0.0;
    Arrays.sort(s1a);
    Arrays.sort(s2a);
    // For every possible combination of unique characters in a and b
    for (int aChar = 0; aChar < s1a.length; aChar++) {
      // If we have seen this character before, skip it
      if (!s1a[aChar].equals(Constants.GAP_SYMBOL) && (aChar == 0 || !s1a[aChar].equals(s1a[aChar - 1]))) {
        double freqA = (double) countCharFreq(s1a, s1a[aChar]);
        for (int bChar = 0; bChar < s2a.length; bChar++) {
          // If we have seen this character before, skip it
          if (!s2a[bChar].equals(Constants.GAP_SYMBOL) && (bChar == 0 || !s2a[bChar].equals(s2a[bChar - 1]))) {
            double freqB = (double) countCharFreq(s2a, s2a[bChar]);
            // If we encounter a gap, use gapPenalty, otherwise look up score
            double substitutionScore = (double) sf.lookupScore(s1a[aChar], s2a[bChar]);

            summation +=
                substitutionScore * freqA * freqB;
          }
        }
      }

    }
    double numGaps =
        (double) countCharFreq(s1a, Constants.GAP_SYMBOL)
            + countCharFreq(s2a, Constants.GAP_SYMBOL);
    Srand = summation / L - numGaps * this.gapPenalty;

    // TODO: Ask Robert how to handle case where Smax = Srand or Sobs = Srand?
    return -1 * Math.log((Math.max(Sobs - Srand,.000001)) / (Math.max(Smax - Srand,.000001)));
  }

  /**
   * Returns the number of occurrences of c in arr.
   * @param arr
   * @param c
   * @return
   */
  protected int countCharFreq(Character[] arr, Character c) {
    int freq = 0;
    for (Character ch : arr) {
      if (ch.equals(c)) {
        freq++;
      }
    }
    return freq;
  }
}
