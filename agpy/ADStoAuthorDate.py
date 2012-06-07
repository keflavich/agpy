import pybtex
from pybtex.database.input import bibtex
import pybtex.database.output.bibtex
import itertools

class Writer(pybtex.database.output.bibtex.Writer):

    def quote(self, s):
        """
        >>> w = Writer()
        >>> print w.quote('The World')
        "The World"
        >>> print w.quote(r'The \emph{World}')
        "The \emph{World}"
        >>> print w.quote(r'The "World"')
        {The "World"}
        >>> try:
        ...     print w.quote(r'The {World')
        ... except BibTeXError, error:
        ...     print error
        String has unmatched braces: The {World
        """

        self.check_braces(s)
        # REMOVED if '"' not in s:
        # REMOVED     return '"%s"' % s
        # REMOVED else:
        return '{%s}' % s


def ADStoAuthorDate(infilename,outfilename):
    """
    Changes the citation key from whatever it is (e.g., 2012MNRAS.416..465L) to
    AuthorDATE (e.g. Longmore2011)
    """

    parser = bibtex.Parser()

    bib_data = parser.parse_file(infilename)
    new_bib_data = pybtex.database.BibliographyData()

    for key in bib_data.entries.keys():
        entry = bib_data.entries[key]
        entry.fields['key'] = key
        new_key = entry.persons['author'][0].last()[0].strip("{}") + entry.fields['year']
        entry.key = new_key

        # add lower-case letter suffix for repeated name/date combinations
        alphabet = itertools.cycle("abcdefghijklmnopqrstuvqxyz")
        while new_key in new_bib_data.entries:
            new_key += alphabet.next()

        new_bib_data.add_entry(new_key, entry)

    writer = Writer()

    writer.write_file(new_bib_data,outfilename)

if __name__ == "__main__":

    import optparse

    parser = optparse.OptionParser()

    options,args = parser.parse_args()

    ADStoAuthorDate(*args)
