package exceptions;

public class BioInfException extends Exception {
  
  public String expected, found;
  
  public BioInfException(String expected, String found) {
    super();
    this.expected = expected;
    this.found = found;
  }
  
  public BioInfException(String msg) {
    super();
    this.expected = msg;
    this.found = "";
  }
  
  public String getError() {
    return expected + " " + found;
  }

}
