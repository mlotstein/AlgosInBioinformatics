package pairwisealignment;
import utilities.Constants;
import exceptions.BioInfException;


public class ScoringFunction {

  public Boolean useBlosum;

  public ScoringFunction(Boolean type) throws BioInfException {
    this.useBlosum = type;
  }

  public Integer lookupScore(char c1, char c2)
    throws BioInfException {
    int tableIndex1 = Constants.AMINO_ACIDS.indexOf(c1);
    int tableIndex2 = Constants.AMINO_ACIDS.indexOf(c2);
    // If either of these characters does not occur in string of acceptable
    // amino acids, throw a BadProteinCharacterException
    if (tableIndex1 < 0 || tableIndex2 < 0) {
      throw new BioInfException(c1 + " or " + c2 + " invalid Amino Acids.");
    }

    if (!this.useBlosum) {
      return Constants.PAM_MATRIX[tableIndex1][tableIndex2];
    } else {
      return Constants.BLOSUM_MATRIX[tableIndex1][tableIndex2];
    } /*else if (this.useBlosum.equals(Constants.LEVENSHTEIN)) {
      return (c1 == c2) ? 1 : 0;
    }*/
    // Exception here?
  }

  public float scoreAlignment(Character[] arr1, Character[] arr2,
    int gapPenalty) {
    float score = 0;

    for (int cNum = 0; cNum < arr1.length & cNum < arr2.length; cNum++) {
      if (arr1[cNum].equals(Constants.GAP_SYMBOL)
          || arr2[cNum].equals(Constants.GAP_SYMBOL)) {
        score -= gapPenalty;
      } else {
        try {
          score += this.lookupScore(arr1[cNum], arr2[cNum]);
        } catch (BioInfException e) {
          e.printStackTrace();
        }
      }
    }

    return score;
  }
}
