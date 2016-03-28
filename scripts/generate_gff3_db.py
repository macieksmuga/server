"""
Generate a SQLite DB file from a sequence annotation file in GFF3 format
(http://sequenceontology.org/resources/gff3.html)
The resulting database file will be mapped to a FeatureSet in the current
GA4GH server implementation, and each row in the 'feature' table will be
mapped to a Feature.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import sys
import json
import sqlite3

import utils
import ga4gh.gff3Parser as gff3

# TODO: Shift this to use the Gff3DbBackend class.

# The columns of the FEATURE table correspond to the columns of a GFF3,
# with three additional columns prepended representing the ID of this feature,
# the ID of its parent (if any), and a whitespace separated array
# of its child IDs.
_dbTableSQL = (
                "CREATE TABLE FEATURE( "
                "id TEXT PRIMARY KEY NOT NULL, "
                "parent_id TEXT, "
                "child_ids TEXT, "
                "reference_name TEXT, "
                "source TEXT, "
                "feature_type TEXT, "
                "start INT, "
                "end INT, "
                "score REAL, "
                "strand TEXT, "
                "attributes TEXT);"
                )


def _db_serialize(pyData):
    return json.dumps(pyData, separators=(',', ':'))


class Gff32Db(object):
    """
    Represents a unit of work for this script: Parse a GFF3 file
    using an external GFF3 parser, and create a corresponding SQLite DB file
    by iterating through the resulting parsed dictionary object
    (gff3Data.byFeatureId).
    """
    def __init__(self, inputFile, outputFile):
        """
        :param inputFile: source GFF3 filename (can be a full path)
        :param outputFile: destination sqlite filename (ditto)
        """
        self.gff3File = inputFile
        self.dbFile = outputFile
        if os.path.exists(outputFile):
            print("DB output file already exists, please remove or rename.",
                  file=sys.stderr)
            exit()

            # # FIXME: No current way to get childIDs in correct order.
            # childIds = [child.featureId for child in feature.children]
            # rowInsertSQL += "'{}', ".format(_db_serialize(childIds))

    def run(self):
        print("Parsing file {}...".format(self.gff3File), file=sys.stderr)
        dbconn = sqlite3.connect(self.dbFile)
        dbcur = dbconn.cursor()
        dbcur.execute(_dbTableSQL)  # create table
        gff3Parser = gff3.Gff3Parser(self.gff3File)
        for feature in gff3Parser.parse():
            parentIds = [parent.featureId for parent in feature.parents]
            parentId = parentIds[0] if len(parentIds) > 0 else ''
            dbcur.execute(
                "INSERT INTO FEATURE VALUES(?,?,?,?,?,?,?,?,?,?,?)", (
                    feature.featureId,
                    parentId,
                    _db_serialize([]),  # FIXME childIds
                    feature.seqname,
                    feature.source,
                    feature.type,
                    feature.start,
                    feature.end,
                    feature.score,
                    feature.strand,
                    _db_serialize(feature.attributes)))
            dbconn.commit()
        dbcur.close()
        dbconn.close()
        print("Done.", file=sys.stderr)


@utils.Timed()
def main():
    parser = argparse.ArgumentParser(
        description="Script to generate SQLite database corresponding to "
        "an input GFF3 sequence annotations file.")
    parser.add_argument(
        "--outputFile", "-o", default="annotations.db",
        help="The file to output the server-ready database to.")
    parser.add_argument(
        "--inputFile", "-i",
        help="Path to input GFF3 file.",
        default='.')
    parser.add_argument('--verbose', '-v', action='count', default=0)
    args = parser.parse_args()
    g2d = Gff32Db(args.inputFile, args.outputFile)
    g2d.run()


if __name__ == "__main__":
    main()
