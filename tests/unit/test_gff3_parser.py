"""
Data-driven tests for rna quantification.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import ga4gh.gff3 as gff3
import unittest


class TestGff3ParserOnTypicalFile(unittest.TestCase):
    """
    Data driven unit tests for the GFF3 parser
    """
    def setUp(self):
        self.testDataFile = "tests/data/datasets/dataset1/sequenceAnnotation/gencodeV21Set1.gff3"
        self.gff3Parser = gff3.Gff3Parser(self.testDataFile)
        self.gff3Data = self.gff3Parser.parse()

    def testFileParsedHasSomeRootFeatures(self):
        self.assertIsNotNone(self.gff3Data.roots, "No root features")
        self.assertNotEqual(len(self.gff3Data.roots),0,"No root features")

    def testSomeFeatureIsWellFormed(self):
        featId = self.gff3Data.byFeatureId.keys()[0]
        feat = self.gff3Data.byFeatureId[featId][0]
        self.assertEqual(feat.gff3Set, self.gff3Data, "feature gff3 set mismatch")
        self.assertEqual(featId, feat.featureId, "featureId mismatch")
        self.assertIsNotNone(feat.seqname, "sequence name is not populated")
        self.assertGreaterEqual(feat.end, feat.start, "end less than start")
        self.assertIn(feat.strand, u"+-", "strand is neither + nor -")
        self.assertIsNotNone(feat.source, "source is unspecified")
        self.assertIsNotNone(feat.type, "feature type is unspecified")
        self.assertGreater(feat.lineNumber, 0, "line number invalid")
        self.assertIsInstance(feat.parents, set, "parents not a set")
        self.assertIsInstance(feat.children, set, "children not a set")


    def testRootFeaturesHaveNoParents(self):
        for root in self.gff3Data.roots:
            self.assertEqual(len(root.parents), 0, "root feature has a parent")


    def testAllFeaturesContainAllRootFeatures(self):
        for root in self.gff3Data.roots:
            feat = self.gff3Data.byFeatureId[root.featureId]
            self.assertGreaterEqual(len(feat), 1,
                "root feature not in list of all features")

    def testInvalidFeatureIdKeyQueryFails(self):
        badFeatureId = "987654"
        badFeat = self.gff3Data.byFeatureId[badFeatureId]
        self.assertEqual(len(badFeat), 0, "invalid feature ID returned valid object")

    def testAllChildrenFeaturesArePresentInSet(self):
        for featList in self.gff3Data.byFeatureId.values():
            for feat in featList:
                for child in feat.children:
                    childLookup = self.gff3Data.byFeatureId[child.featureId]
                    self.assertGreaterEqual(len(childLookup), 1,
                                            "child feature not in set")


