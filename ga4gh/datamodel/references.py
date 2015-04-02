"""
Module responsible for translating reference sequence data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob

import numpy
import pyfasta

import ga4gh.protocol as protocol
import ga4gh.datamodel as datamodel
import ga4gh.sidegraph as sidegraph

class AbstractReferenceSet(datamodel.DatamodelObject):
    """
    Class representing either linear or graph ReferenceSets.
    """
    def __init__(self, id_):
        self._id = id_

    def toProtocolElement(self):
        """
        Returns a default GA4GH protocol representation of this 
        ReferenceSet.
        """
        ret = protocol.ReferenceSet()
        ret.id = self._id
        ret.description = "TODO"
        ret.sourceURI = None
        ret.assemblyId = None
        ret.md5checksum = "TODO"
        ret.ncbiTaxonId = None
        ret.sourceAccessions = []
        return ret

    def _createReference(self, mdict={}):
        """
        Returns the GA4GH protocol representation of this Reference.
        mdict is an optional metadata dictionary (as returned from 
        SideGraph.getReference(), for example)
        """
        reference = protocol.Reference()
        reference.id = self._id
        reference.sequenceId = mdict.get("sequenceId",self._id)
        reference.start = mdict.get("start", 0)
        reference.length = mdict.get("length",0)
        reference.md5checksum = mdict.get("md5checksum", "") 
        reference.name = mdict.get("name", "")
        reference.sourceAccessions = mdict.get("sourceAccessions", [])
        reference.isDerived = mdict.get("isDerived", False)
        reference.sourceDivergence = mdict.get("sourceDivergence", None)
        reference.ncbiTaxonId = mdict.get("ncbiTaxonId", None)
        reference.isPrimary = mdict.get("isPrimary",True)
        return reference


class LinearReferenceSet(AbstractReferenceSet):
    """
    Class implementing linear ReferenceSets. Such a ReferenceSet is
    a set of LinearReferences which typically comprise a reference
    assembly, such as GRCh38.
    """
    def __init__(self, id_, dataDir):
        self._id = id_
        self._referenceIdMap = {}
        self._dataDir = dataDir
        # TODO get metadata from a file within dataDir? How else will we
        # fill in the fields like ncbiTaxonId etc?
        for relativePath in glob.glob(os.path.join(self._dataDir, "*.fa.gz")):
            filename = os.path.split(relativePath)[1]
            localId = filename.split(".")[0]
            referenceId = "{}:{}".format(self._id, localId)
            # TODO What was here is wrong. It's still wrong.
            self._referenceIdMap[referenceId] = filename
        self._referenceIds = sorted(self._referenceIdMap.keys())

    def getReferences(self):
        """
        Returns the References in this ReferenceSet.
        """
        return self._referenceIdMap.values()


class GraphReferenceSet(AbstractReferenceSet):
    """
    Implements graph reference sets: Reference sequences together with
    joins, as well as some helper methods to deal with complex queries.
    Makes assumption that dataDir used for initialization contains one
    graph topology Sqlite database (.db file) together with any named 
    FASTA files containing the actual sequence string data.
    """
    def __init__(self, id_, dataDir):
        self._id = id_
        self._dataDir = dataDir
        self._referenceIdMap = {}
        # TODO: Throw exception if directory contains != 1 database file
        self._dbFile = glob.glob(os.path.join(self._dataDir, "*.db"))[0]

    def getReferences(self):
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            for refName in sg.getRefNames():
                sgRefDict = sg.getReference(refName)
                yield self._createReference(sgRefDict)

        




