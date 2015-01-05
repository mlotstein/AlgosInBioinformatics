//Copyright 2012, University of Freiburg,
//Chair of Algorithms and Data Structures.
//Author: Max Lotstein <losteim@informatik.uni-freiburg.de>
//Extended from code from Florian  Burle

package utilities;

import java.util.Arrays;

import exceptions.BioInfException;

/**
 * A resizable-array implementation for ProteinSequences.
 */
public class ProteinSequenceList {

  /**
   * The default capacity for an integer list is 10.
   */
  private static final int DEFAULT_CAPACITY = 10;

  private ProteinSequence[] list;
  private int size;

  /**
   * Constructs an empty list with an initial capacity of ten.
   */
  public ProteinSequenceList() {
    this.initialize(DEFAULT_CAPACITY);
  }

  /**
   * Constructs a list containing the elements of the specified array.
   */
  public ProteinSequenceList(ProteinSequence[] array) {
    if (array == null) {
      throw new NullPointerException("The array must not be null!");
    }
    this.list = Arrays.copyOf(array, Math.max(array.length, DEFAULT_CAPACITY));
    this.size = array.length;
  }

  /**
   * Constructs an empty list with the specified initial capacity.
   */
  public ProteinSequenceList(int initialCapacity) {
    if (initialCapacity < 0) {
      throw new IllegalArgumentException("The initial capacity must be >= 0!");
    }
    this.initialize(initialCapacity);
  }

  /**
   * Initialized this list with the specified initial capacity.
   */
  private void initialize(int initialCapacity) {
    this.list = new ProteinSequence[initialCapacity];
    this.size = 0;
  }

  /**
   * Doubles the internal array storage capacity. If the current capacity is 0
   * it is increased to the {@link #DEFAULT_CAPACITY}.
   */
  private void increaseStorage() {
    int newCapacity = this.list.length * 2;
    if (newCapacity == 0) {
      newCapacity = DEFAULT_CAPACITY;
    }
    this.list = Arrays.copyOf(this.list, newCapacity);
  }

  /**
   * Appends the given integer to the end of this list.
   */
  public void add(ProteinSequence value) {
    if (this.size == this.list.length) {
      this.increaseStorage();
    }

    this.list[this.size] = value;
    ++this.size;
  }

  /**
   * Removes all of the elements from this list.
   */
  public void clear() {
    this.size = 0;
  }

  /**
   * Returns the element at the specified position in this list.
   */
  public ProteinSequence get(int index) {
    if (index >= this.size) {
      throw new IndexOutOfBoundsException(
          "The given index was outside of the size of the list.");
    }
    return this.list[index];
  }

  /**
   * Returns the element at the specified position in this list.<br/>
   * <b>Caution: This method makes no range check for the index</b> for faster
   * access! Only use this method if you are sure that there is a valid element
   * at this index!
   */
  public ProteinSequence getUnsafe(int index) {
    return this.list[index];
  }

  /**
   * Replaces the element at the specified position in this list with the
   * specified value.
   */
  public void set(int index, ProteinSequence value) {
    if (index >= this.size) {
      throw new IndexOutOfBoundsException(
          "The given index was outside of the size of the list.");
    }
    this.list[index] = value;
  }

  /**
   * Removes the last element from the list and returns it.
   */
  public ProteinSequence pop() {
    if (this.size == 0) {
      throw new IllegalStateException("The ProteinSequence list is empty!");
    }

    --this.size;
    return this.list[this.size];
  }

  /**
   * Returns the last element from the list.
   */
  public ProteinSequence last() {
    if (this.size < 1) {
      throw new IllegalStateException("The integer list is empty!");
    }
    return this.list[this.size - 1];
  }

  /**
   * Returns the current capacity of this list.
   */
  protected int capacity() {
    return this.list.length;
  }

  /**
   * Returns the number of elements in this list.
   */
  public int size() {
    return this.size;
  }

  /**
   * Finds proteinSequence with id and removes it from array.
   * 
   * @param id
   */
  public void removeById(String id) {
    int pos = -1;

    for (int idPos = 0; pos < 0 && idPos < this.size; idPos++) {
      if (this.list[idPos].id.equals(id)) {
        pos = idPos;
      }
    }

    if (pos != -1) {
      for (int psNum = pos; psNum < this.size; psNum++) {
        this.list[psNum] = this.list[psNum + 1];
      }
      this.size = this.size - 1;
    }
  }

  /**
   * Returns <code>true</code> if this list contains no elements.
   */
  public boolean isEmpty() {
    return this.size == 0;
  }

  /**
   * Trims the capacity of this {@link StringList} instance to be the list's
   * current size.
   */
  public void trimToSize() {
    if (this.size < this.list.length) {
      this.list = Arrays.copyOfRange(this.list, 0, this.size);
    }
  }

  /**
   * Sorts this list into ascending numerical order.
   */
  public void sort() {
    if (this.size == this.capacity()) {
      Arrays.sort(this.list);
    } else {
      // Extract an array that only contains the real content of the list and
      // sort it.
      ProteinSequence[] sortedList = Arrays.copyOf(this.list, this.size);
      Arrays.sort(sortedList);
      // Copy the extracted sorted array back into the underlying list with it's
      // original size.
      this.list = Arrays.copyOf(sortedList, this.list.length);
    }
  }

  /**
   * Returns an array containing all of the elements in this list in proper
   * sequence (from first to last element).
   */
  public ProteinSequence[] toArray() {
    return Arrays.copyOf(this.list, this.size);
  }

  @Override
  public String toString() {
    StringBuilder result = new StringBuilder("[");
    if (this.size > 0) {
      result.append(this.list[0]);
    }
    for (int i = 1; i < this.size; ++i) {
      result.append(", " + this.list[i]);
    }
    result.append("]");
    return result.toString();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof ProteinSequenceList)) {
      return false;
    }
    ProteinSequenceList l = (ProteinSequenceList) obj;
    return Arrays.equals(this.list, l.list);
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  public static ProteinSequenceList parse(String[] lines)
    throws BioInfException {
    ProteinSequenceList psList = new ProteinSequenceList();
    for (int lineNum = 0; lineNum < lines.length; lineNum += 2) {
      if (lines[lineNum].charAt(0) != '>') {
        throw new BioInfException("Improper FASTA format");
      }

      String identifier = lines[lineNum].substring(1).trim();
      // Give the sequence a name if it doesn't have one.
      if (identifier.length() == 0) {
        identifier = "sequence" + (psList.size() + 1);
      }

      ProteinSequence ps = new ProteinSequence(identifier, lines[lineNum + 1]);
      psList.add(ps);

    }

    return psList;
  }
}
