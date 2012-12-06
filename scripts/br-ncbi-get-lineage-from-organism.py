#!/usr/bin/env python
from HTMLParser import HTMLParser
from urllib import urlopen
import argparse
import warnings
import sys

NCBI_TAXONOMY_BROWSER_ID = "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={tax_id}"
NCBI_TAXONOMY_BROWSER_NAME = "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name={name}"
TAX_ID_LINE_START = '<em>Taxonomy ID: </em>'
TAX_ID_END = '<br'
LINEAGE_LINE_START = '<dd><a ALT="no rank"'
TOPNAME_LINE_START = '<table width="100%"><tr><td valign="top"><h2>'
TOPNAME_END = '</h2>'
# Tax levels vary between organisms on the NCBI Taxonomy browser, to get the same number of columns for
# all given species this global constant is used.
TAX_LEVELS = ['no rank', 'superkingdom', 'superphylum', 'phylum', 'class', 'order', 'family', 'genus', 'species']


class NCBILinageParser(HTMLParser):
    def __init__(self):
        self.taxon_a = False
        self.taxon_lvl = None
        self.phylogeny = list()
        HTMLParser.__init__(self)

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for attr in attrs:
                if attr[0] == "title":
                    self.taxon_a = True
                    self.taxon_lvl = attr[1]

    def handle_data(self, data):
        if self.taxon_a:
            self.phylogeny.append((self.taxon_lvl, " ".join(data.split())))
            self.taxon_a = False


def main(tax_id, is_name=False, print_header=False, print_table=False, print_taxid=False):
    if not is_name:
        filehandle = urlopen(NCBI_TAXONOMY_BROWSER_ID.format(tax_id=tax_id))
    else:
        filehandle = urlopen(NCBI_TAXONOMY_BROWSER_NAME.format(name=tax_id).replace(' ', '+'))

    for line in filehandle:
        if print_taxid and line.startswith(TAX_ID_LINE_START):
            taxid_line = line[len(TAX_ID_LINE_START):]
            sys.stdout.write(taxid_line[:taxid_line.find(TAX_ID_END)].strip() + "\t")
        elif line.startswith(TOPNAME_LINE_START):
            topname_line = line[len(TOPNAME_LINE_START):]
            topname = topname_line[:topname_line.find(TOPNAME_END)].strip()
        elif line.startswith(LINEAGE_LINE_START):
            nparser = NCBILinageParser()
            nparser.feed(line)

            if print_table:
                assert(len(set(TAX_LEVELS)) == len(TAX_LEVELS))

                # Print warning if each taxa does not appear only once, take
                # the one highest in the hierarchy
                # NOTE: no rank tends to appear twice, e.g. for id 290317
                taxa = [t for (t, v) in nparser.phylogeny]
                if not len(taxa) == len(set(taxa)):
                    warnings.warn("Non-unique taxa level name."
                                  "Taking highest in hierarchy. See: {taxa}"
                                  "".format(taxa=taxa))
                phyldict = dict(reversed(nparser.phylogeny))

                sys.stdout.write("\t".join(phyldict.setdefault(tl, "-") for tl in TAX_LEVELS))
                sys.stdout.write("\t" + topname + "\n")
            else:
                if print_header:
                    if print_taxid:
                        sys.stdout.write("taxonomy_id\t")
                    print("\t".join([t for (t, v) in nparser.phylogeny]) + "\ttopname")
                print("\t".join([v for (t, v) in nparser.phylogeny]) + "\t" + topname)
            return

    raise Exception('Line with lineage information not found.')


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Output full linage from an NCBI Taxonomy ID as tab "
                                        "separated values. Last value is the name at the top of the page,"
                                        "usually a strain.")
    argparser.add_argument('tax_ids', metavar='TaxonomyID', type=str, nargs='+',
                           help='One or more NCBI Taxonomy IDs. Can also be names if --name is supplied.')
    argparser.add_argument("--header", action="store_true", help="Print header for every record")
    argparser.add_argument("--table", action="store_true",
                           help="Print with equal number of taxonomic levels."
                           "Insert '-' if level is not available. NOTE this "
                           "does not include all possible taxonomic levels. "
                           "Prints {num} columns, i.e. {ranks}"
                           "".format(num=len(TAX_LEVELS), ranks=", ".join(TAX_LEVELS)))
    argparser.add_argument("--name", action="store_true", help="Given ids are names")
    argparser.add_argument("--printid", action="store_true", help="Prints Taxonomy ID as the first value")
    args = argparser.parse_args()

    if args.header and args.table:
        if args.printid:
            sys.stdout.write("taxonomy_id\t")
        print("\t".join(TAX_LEVELS) + "\ttopname")
    for tid in args.tax_ids:
        main(tid, is_name=args.name, print_header=args.header,
             print_table=args.table, print_taxid=args.printid)
