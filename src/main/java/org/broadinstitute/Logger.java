package org.broadinstitute;

public class Logger {

  private Level m_level;

  public void setLevel(Level level){
    m_level = level;
  }

  public void info(String str){
    System.out.println(str);
  }

  public void debug(String str){
    if(m_level == Level.DEBUG){
      System.out.println(str);
    }
  }
}
