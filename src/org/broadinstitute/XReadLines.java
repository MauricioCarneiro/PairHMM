
import java.io.*;
import java.util.Iterator;

/*************************************************************************
 * Classes and Functions used to parse and apply the tests to the
 * PairHMM class
 ************************************************************************/

class XReadLines implements Iterator<String>, Iterable<String> {
    private final BufferedReader in;      // The stream we're reading from
    private String nextLine = null;       // Return value of next call to next()
    private final boolean trimWhitespace;
    private final String commentPrefix;

    public XReadLines(final File filename) throws FileNotFoundException {
        this(new FileReader(filename), true, null);
    }

    /**
     * Creates a new xReadLines object to read lines from an bufferedReader
     *
     * @param reader         file name
     * @param trimWhitespace trim whitespace
     * @param commentPrefix  prefix for comments or null if no prefix is set
     */
    public XReadLines(final Reader reader, final boolean trimWhitespace, final String commentPrefix) {
        this.in = (reader instanceof BufferedReader) ? (BufferedReader) reader : new BufferedReader(reader);
        this.trimWhitespace = trimWhitespace;
        this.commentPrefix = commentPrefix;
        try {
            this.nextLine = readNextLine();
        } catch (IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    /**
     * I'm an iterator too...
     *
     * @return an iterator
     */
    public Iterator<String> iterator() {
        return this;
    }

    public boolean hasNext() {
        return this.nextLine != null;
    }

    /**
     * Actually reads the next line from the stream, not accessible publicly
     *
     * @return the next line or null
     * @throws java.io.IOException if an error occurs
     */
    private String readNextLine() throws IOException {
        String nextLine;
        while ((nextLine = this.in.readLine()) != null) {
            if (this.trimWhitespace) {
                nextLine = nextLine.trim();
                if (nextLine.length() == 0)
                    continue;
            }
            if (this.commentPrefix != null)
                if (nextLine.startsWith(this.commentPrefix))
                    continue;
            break;
        }
        return nextLine;
    }

    /**
     * Returns the next line (optionally minus whitespace)
     *
     * @return the next line
     */
    public String next() {
        try {
            String result = this.nextLine;
            this.nextLine = readNextLine();

            // If we haven't reached EOF yet
            if (this.nextLine == null) {
                in.close();             // And close on EOF
            }

            // Return the line we read last time through.
            return result;
        } catch (IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    // The file is read-only; we don't allow lines to be removed.
    public void remove() {
        throw new UnsupportedOperationException();
    }

}
