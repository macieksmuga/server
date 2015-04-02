"""
Unit tests for sidegraph.py, my toy side graph library.
Uses local data for now. 
TODO: Generalize to take any old SQL/FASTA encoded sidegraph.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import numpy
import pyfasta
import sqlite3
import os

import ga4gh.sidegraph as sidegraph


class TestCase(unittest.TestCase):
    def setUp(self):
        self._files = "tests/data/graph_references/graph_1"
        self._db = os.path.join(self._files,"refgraph1.db")
        self._conn = sqlite3.connect(self._db)

    def tearDown(self):
        self._conn.close()

    def testGetRefNames(self):
        c = self._conn.cursor()
        q = c.execute("SELECT name FROM Reference")
        with sidegraph.SideGraph(self._db, self._files) as sg:
            self.assertEquals(
                q.fetchall().sort(),
                sg.getRefNames().sort())

    def testGetJoinsDatasetLength(self):
        c = self._conn.cursor()
        q = c.execute("SELECT COUNT(*) FROM GraphJoin")
        
        with sidegraph.SideGraph(self._db, self._files) as sg:
            self.assertEquals(q.fetchone()[0], len(sg.getJoins()))

    def testGetSubgraphLength(self):
        with sidegraph.SideGraph(self._db, self._files) as sg:
            (segments, joins) = sg.getSubgraph('chr1',75,10)
            self.assertEquals(2, len(segments))
            self.assertEquals(2, len(joins))


if __name__ == '__main__':
    unittest.main()