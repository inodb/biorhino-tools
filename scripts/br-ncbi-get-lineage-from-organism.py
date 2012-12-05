#!/usr/bin/env python
from HTMLParser import HTMLParser
from urllib import urlopen
import argparse

NCBI_TAXONOMY_BROWSER_ID = "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={tax_id}"
NCBI_TAXONOMY_BROWSER_NAME = "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name={name}"


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


def main(tax_id, is_name=False, print_header=False):
    if not is_name:
        filehandle = urlopen(NCBI_TAXONOMY_BROWSER_ID.format(tax_id=tax_id))
    else:
        filehandle = urlopen(NCBI_TAXONOMY_BROWSER_NAME.format(name=tax_id).replace(' ', '+'))

    for line in filehandle:
        if line.startswith('<dd><a ALT="no rank"'):
            nparser = NCBILinageParser()
            nparser.feed(line)
            if print_header:
                print("\t".join([tax_lvl for (tax_lvl, v) in nparser.phylogeny]))
            print("\t".join([v for (tax_lvl, v) in nparser.phylogeny]))
            return

    raise Exception('Line with lineage information not found.')


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description='Output full linage from an NCBI Taxonomy ID as tab separated values.')
    argparser.add_argument('tax_ids', metavar='TaxonomyID', type=str, nargs='+',
                           help='One or more NCBI Taxonomy IDs. Can also be names if --name is supplied.')
    argparser.add_argument("--header", action="store_true", help="Print header first")
    argparser.add_argument("--name", action="store_true", help="Given ids are names")
    args = argparser.parse_args()

    if args.header:
        print('no rank\tsuperkingdom\tsuperphylum\tphylum\tclass\torder\tfamily\tgenus\tspecies')
    for tid in args.tax_ids:
        main(tid, is_name=args.name)
