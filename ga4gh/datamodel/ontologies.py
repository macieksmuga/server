"""
Support for Ontologies.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import fnmatch
import os
import ga4gh.protocol as protocol
import ga4gh.datamodel as datamodel


class FileSystemOntology(object):
    """
    The base class of storing an ontology
    At this moment an "Ontology" is just a map between names and IDs (e.g.
    in Sequence Ontology we would have "SO:0001583 <-> missense_variant")
    This is a tempotrary solution adn must be replaced by more comprehensive
    ontology object.
    """

    def __init__(self):
        self._nameIdMap = dict()
        self._idNameMap = dict()

    def add(self, id_, name):
        self._nameIdMap[id_] = name
        self._idNameMap[name] = id_

    def getId(self, name):
        return self._idNameMap[name]

    def getName(self, id_):
        return self._nameIdMap[id_]

    def readOntology(self, filename):
        with open(filename) as f:
            for line in f:
                # File format: id \t name
                fields = line.rstrip().split("\t")
                self.add(fields[0], fields[1])


class FileSystemOntologies(object):
    """
    An ontology read from the file system.
    This implementation uses a tab separated TXT file: "id\tname"
    """

    def __init__(self, localId, dataDir, backend):
        self._ontologyNameMap = dict()
        self.readOntologies(dataDir)

    def add(self, ontologyName, ontology):
        self._ontologyNameMap[ontologyName] = ontology

    def get(self, ontologyName):
        return self._ontologyNameMap[ontologyName]

    def keys(self):
        return self._ontologyNameMap.keys()

    def len(self):
        return len(self._ontologyNameMap)

    def readOntologies(self, dataDir):
        self._dataDir = dataDir
        # Find TXT files
        for filename in os.listdir(dataDir):
            if fnmatch.fnmatch(filename, '*.txt'):
                ontologyName, _ = os.path.splitext(filename)
                path = os.path.join(dataDir, filename)
                ontology = FileSystemOntology()
                ontology.readOntology(path)
                self.add(ontologyName, ontology)


class OntologyTermSet(datamodel.DatamodelObject):
    """
    A set of related ontology terms.
    """
    def __init__(self, ontologySource):
        self._ontologySource = ontologySource
        self._ontologyTermIdMap = dict()
        self._ontologyTermNameMap = dict()

    def add(self, ontologyTerm):
        """
        add an ontology term to the object
        """
        if ontologyTerm.id in self._ontologyTermIdMap:
            raise ValueError(
                "OntologyTerm id already exists: {}".format(str(ontologyTerm)))
        self._ontologyTermIdMap[ontologyTerm.id] = ontologyTerm
        if ontologyTerm.name is not None:
            if ontologyTerm.name in self._ontologyTermNameMap:
                raise ValueError(
                    "OntologyTerm already exists: {}".format(str(
                                                             ontologyTerm)))

    def create(self, id, name):
        """
        create a new ontology term
        """
        self.add(OntologyTerm(self._ontologySource, id, name))


class OntologyTerm(datamodel.DatamodelObject):
    """
    A specific ontology term.
    """
    def __init__(self, ontologySource, id, name):
        self._ontologySource = ontologySource
        self._id = id
        self._name = name

    def __str__(self):
        return "{}/{}/{}".format(self.ontologySource, self.id, str(self.name))

    def toProtocolElement(self):
        """
        Returns the representation of this OntologyTerm as the corresponding
        ProtocolElement.
        """
        gaOntologyTerm = protocol.OntologyTerm()
        gaOntologyTerm.ontologySource = self._ontologySource
        gaOntologyTerm.id = self._id
        gaOntologyTerm.name = self._name
        return gaOntologyTerm


class SequenceOntologyTermSet(OntologyTermSet):
    # TODO: tmp hack to encode the sequence ontology terms used by the BRCA
    # gene sets
    def __init__(self):
        super(SequenceOntologyTermSet, self).__init__(
            "http://www.sequenceontology.org/")
        self.create("SO:0000316", "CDS")
        self.create("SO:0000147", "exon")
        self.create("SO:0000704", "gene")
        self.create("SO:0000318", "start_codon")
        self.create("SO:0000319", "stop_codon")
        self.create("SO:0000710", "stop_codon_redefined_as_selenocysteine")
        self.create("SO:0000673", "transcript")
        self.create("SO:0000203", "UTR")

    _singleton = None

    @staticmethod
    def singleton():
        """
        obtain singleton instances of this class
        """
        if SequenceOntologyTermSet._singleton is None:
            SequenceOntologyTermSet._singleton = SequenceOntologyTermSet()
        return SequenceOntologyTermSet._singleton
