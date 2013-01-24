package org.broadinstitute;

/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;


public class ParserIterator implements Iterator<TestRow> {

  private XReadLines m_reader;
  private TestRow m_nextEntry;
  private boolean m_debug;
  
  public ParserIterator(String filePath, boolean debug) throws FileNotFoundException {
    m_reader = new XReadLines(new File(filePath));
    m_debug = debug;
  }
  
  @Override
  public boolean hasNext() {
    if(m_reader.hasNext() == false){
      return false;
    }
    
    String line = m_reader.next();
    String[] p = line.split(" ");
    boolean reachedRead;
    if (p[7].trim().equals("true")) {
      reachedRead = true;
    } else {
      reachedRead = false;
    }

    TestRow testRow = new TestRow(p[0].getBytes(), p[1].getBytes(),
      convertToByteArray(p[2]), convertToByteArray(p[3]),
      convertToByteArray(p[4]), convertToByteArray(p[5]),
      Integer.parseInt(p[6]), reachedRead,
      Double.parseDouble(p[8]));

    m_nextEntry = testRow;
    
    if (m_debug) {
      System.out.println("Haplotype Bases:" + p[0]);
      System.out.println("Read Bases:" + p[1]);
      System.out.println("Read Quals: " + p[2]);
      System.out.println("Read Ins Quals: " + p[3]);
      System.out.println("Read Del Quals: " + p[4]);
      System.out.println("Overall GCP:" + p[5]);
      System.out.println("Haplotype start:" + p[6]);
      System.out.println("Boolean (reached read values):" + p[7]);
      System.out.println("result, likelihood:" + p[8]);
    }
    return true;
  }
  
  private byte[] convertToByteArray(String quals) {
    byte[] output = quals.getBytes();

    for (int i = 0; i < output.length; i++) {
      output[i] -= 33;
    }

    return output;
  }

  @Override
  public TestRow next() {
    return m_nextEntry;
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException("Not supported yet.");
  }
  
}
