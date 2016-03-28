"""
Tests of the gff3 Extract-Transform-Load (ETL) utility script.

Since the output of each run of the script is a DB file, the tests
themselves involve opening the resulting files in sqlite3 and
making specific value queries.

*** This test must be run from the project root directory, ex:
<$PROJECT_ROOT_DIR>: nosetests tests.integration.test_generate_gff3_db
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import tempfile
import shutil
import subprocess
import sqlite3
import json

class TestGff3EtlScript(unittest.TestCase):
    gff3Dir = 'tests/data/datasets/dataset1/sequenceAnnotations'
    etlScript ='python scripts/generate_gff3_db.py'
    testFiles = (
        'gencodeV21Set1',
        'sacCerTest',
        'specialCasesTest',
        'discontinuous'
    )

    def callEtlScript(self, fileName):
        """
        Calls the script to test on the input file, produces output DB
        file <fileName>.db in the temporary output directory.
        :param fileName: filename without extension.
        """
        subprocess.call(
            "{} -i {}/{}.gff3 -o {}/{}.db".format(
                self.etlScript,
                self.gff3Dir, fileName,
                self.tmpDbDir, fileName),
            shell=True)

    def setUp(self):
        self.tmpDbDir = tempfile.mkdtemp()
        print("temporary directory: {}".format(self.tmpDbDir))
        self.dbs = {}
        for testFile in self.testFiles:
            self.callEtlScript(testFile)
            dbConn = sqlite3.connect(
                "{}/{}.db".format(self.tmpDbDir, testFile))
            self.dbs[testFile] = dbConn

    def tearDown(self):
        for dbConn in self.dbs.values():
            dbConn.close()
        shutil.rmtree(self.tmpDbDir)

    def testGencodeAllRowsRead(self):
        sql = "SELECT count(*) FROM Feature"
        query = self.dbs['gencodeV21Set1'].execute(sql)
        ret = query.fetchone()
        self.assertEqual(ret[0], 575)

    def testGencodeGetTranscriptById(self):
        sql = (
            "SELECT * FROM Feature "
            "WHERE id='ENST00000456328.2'")
        query = self.dbs['gencodeV21Set1'].execute(sql)
        ret = query.fetchone()
        self.assertIsNotNone(ret)
        self.assertEqual(ret[1], "ENSG00000223972.5")  # parentID
        self.assertItemsEqual(json.loads(ret[2]), [
            "exon:ENST00000456328.2:1",
            "exon:ENST00000456328.2:2",
            "exon:ENST00000456328.2:3"])  # children
        self.assertEqual(ret[5], "transcript")  # feature type
        self.assertEqual(ret[6], 11869)  # start
        self.assertEqual(ret[7], 14409)  # end

    # def testDiscontinuousResolvedCorrectly(self):
