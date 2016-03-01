"""
Module responsible for translating sequence annotation data
into GA4GH native objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import ga4gh.exceptions

import ga4gh.protocol as protocol
import ga4gh.datamodel as datamodel
import ga4gh.sqliteBackend as sqliteBackend

# Note to self: There's the Feature ID as understood in a GFF3 file,
# the Feature ID that is its server-assigned compoundId, and the
# ID of the feature's row in the DB FEATURE table.
# I need to be careful about which one is which.

"""
    For this implementation, `featureSetId` is required, while `parentId`
    is optional, and filters the features within the requested `featureSetId`
    by their parent.

    Only return features on the reference with this name. Genomic positions
    are non-negative integers less than reference length.
    Requests spanning the join of circular genomes are represented as two
    requests one on each side of the join (position 0) end is also required
    If specified, this query matches only annotations which match one of the
    provided feature types.
    For now do not use the features array in search

    GFF3 data is represented by rows in a single table, named FEATURE.
    The columns of the FEATURE table correspond to the columns of a GFF3,
    with three additional columns prepended representing the ID
    of this feature, the ID of its parent (if any), and a whitespace
    separated array of its child IDs.

    _featureColumns pairs represent the ordered (column_name, column_type).
"""
_featureColumns = [('id', 'TEXT'),
                   ('parent_id', 'TEXT'),
                   ('child_ids', 'TEXT'),
                   ('reference_name', 'TEXT'),
                   ('source', 'TEXT'),
                   ('ontology_term', 'TEXT'),
                   ('start', 'INT'),
                   ('end', 'INT'),
                   ('score', 'REAL'),
                   ('strand', 'TEXT'),  # limited to one of '+'/'-' or none.
                   ('attributes', 'TEXT')  # JSON encoding of attributes dict
                   ]


class Gff3DbBackend(sqliteBackend.SqliteBackedDataSource):
    """
    Notes about the current implementation:
    For this implementation, `featureSetId` is required, while `parentId`
    is optional, and filters the features within the requested `featureSetId`
    by their parent.

    Genomic positions are non-negative integers less than reference length.
    Requests spanning the join of circular genomes are represented as two
    requests one on each side of the join (position 0)
    """

    def __init__(self, dbFile):
        super(Gff3DbBackend, self).__init__(dbFile)
        self.featureColumnNames = [f[0] for f in _featureColumns]
        self.featureColumnTypes = [f[1] for f in _featureColumns]

    def searchFeaturesInDb(self, pageToken=0, pageSize=None, **query):
        """
        :param pageToken: int representing first record to return
        :param pageSize: int representing number of records to return
        :param query: dictionary of query terms and values (as strings)
        to translate into 'WHERE' clauses.
        Keys must exactly match entries in self.featureColumnNames,
        and must be converted according to their featureColumnType.
        All text values must first be checked for illegal characters,
        and for overflow, to prevent SQL injection attacks.
        :return: an array of dictionaries, each representing a feature's
        worth of data, keyed by the column names.
        """

        # TODO: Recast this using sqlite's execute diction.
        # Ex: cur.execute("insert into people values (?, ?)", (who, age))

        sql = "SELECT * FROM FEATURE "
        whereClauses = []
        for col in query:
            if col in self.featureColumnNames:
                colIdx = self.featureColumnNames.index(col)
                colType = self.featureColumnTypes[colIdx]
                colVal = query[col]
                # simple input sanity check
                if colType is "INT":
                    colVal = int(colVal)
                    # start and end need to be handled properly for ranges
                    # featureStart < queryEnd and featureEnd > queryStart
                    if col is "start":
                        whereClauses.append("end > {}".format(colVal))
                    elif col is "end":
                        whereClauses.append("start < {}".format(colVal))
                    else:
                        whereClauses.append("{} = {}".format(col, colVal))
                else:  # TEXT of some sort
                    if "'" in colVal:
                        raise(ga4gh.exceptions.BadRequestException)
                    whereClauses.append("{} = '{}'".format(col, colVal))
        if len(whereClauses) > 0:
            sql += "WHERE {} ".format(" AND ".join(whereClauses))
        sql += "ORDER BY start, id "
        sql += sqliteBackend.limitsSql(pageToken, pageSize)
        query = self._dbconn.execute(sql)
        return sqliteBackend.sqliteRows2dicts(query.fetchall())


class AbstractFeatureSet(datamodel.DatamodelObject):
    """
    A set of sequence features annotations
    """
    compoundIdClass = datamodel.FeatureSetCompoundId

    def __init__(
            self, parentContainer, localId, referenceSetId=None):
        super(AbstractFeatureSet, self).__init__(parentContainer, localId)
        self._referenceSetId = referenceSetId
        self._name = localId
        self._sourceUri = ""
        self._attributes = protocol.Attributes()

    def toProtocolElement(self):
        """
        Returns the representation of this FeatureSet as the corresponding
        ProtocolElement.
        """
        gaFeatureSet = protocol.FeatureSet()
        gaFeatureSet.id = self.getId()
        gaFeatureSet.datasetId = self.getParentContainer().getId()
        gaFeatureSet.referenceSetId = self._referenceSetId
        gaFeatureSet.name = self._name
        gaFeatureSet.sourceUri = self._sourceUri
        gaFeatureSet.attributes = self._attributes
        return gaFeatureSet

    def getFeatureId(self, gaFeature):
        """
        :param gaFeature: protocol.Feature object
        :return: string representing ID for the specified protocol
            Feature object in this FeatureSet.
        """
        compoundId = datamodel.FeatueCompoundId(
            self.getCompoundId(), gaFeature.id)
        return str(compoundId)


class SimulatedFeatureSet(AbstractFeatureSet):
    def __init__(self, parentContainer, localId, randomSeed=1):
        super(SimulatedFeatureSet, self).__init__(
            parentContainer, localId, None)
        self._randomSeed = randomSeed


class Gff3DbFeatureSet(AbstractFeatureSet):
    """
    Stub class to directly read sequence annotation features from GFF3 files.
    Tests basic access, not to be used in production.
    """
    def __init__(self, parentContainer, localId, filePath, dataRepository):
        super(Gff3DbFeatureSet, self).__init__(
            parentContainer, localId, None)

        self._dbFilePath = filePath  # the full file path of the
        self._dataRepository = dataRepository
        self._db = Gff3DbBackend(self._dbFilePath)

    def _gaFeatureForFeatureDbRecord(self, feature):
        """
        :param feature: The DB Row representing a feature
        :return: the corresponding GA4GH protocol.Feature object
        """
        gaFeature = protocol.Feature()
        gaFeature.id = feature['id']
        gaFeature.parentId = feature['parentId']
        gaFeature.featureSetId = self._compoundId
        gaFeature.referenceName = feature['reference_name']
        gaFeature.start = feature['start']
        gaFeature.end = feature['end']
        gaFeature.featureType = feature['ontology_term']
        gaFeature.attributes = json.loads(
            feature['attributes']).toProtocolElement()
        return gaFeature


    def featureObjectGenerator(self, request):
        """
        method passed to runSearchRequest to fulfill the request
        :param request: protocol.FeatureSearchRequest
        :return: yields a protocol.Feature at a time
        """
        # parse out the various query parameters from the request.
        parentId = request.parentId
        referenceName = request.referenceName
        start = 0
        if request.get(start):
            start = int(request.start)
        end = None
        if request.get(end):
            end = int(request.end)
        ontologyTerms = request.ontologyTerms

        with self._db as dataSource:
            # featuresCount is needed to ensure that once the
            # request is fulfilled, no nextPageTokens past the
            # end of the actual dataset range are returned.
            featuresCount = dataSource.countFeaturesSearchInDb(
                parent_id=parentId, reference_name=referenceName,
                start=start, end=end, ontology_terms=ontologyTerms)
            featuresReturned = dataSource.searchFeaturesInDb(
                request.pageToken, request.pageSize,
                parent_id=parentId, reference_name=referenceName,
                start=start, end=end, ontology_terms=ontologyTerms)

        nextPageToken = request.pageToken
        for featureRecord in featuresReturned:
            gaFeature = self._gaFeatureForFeatureDbRecord(featureRecord)
            # pagination logic: None if last feature was returned,
            # else row number.
            if nextPageToken < featuresCount:
                nextPageToken += 1
            else:
                nextPageToken = None

            yield gaFeature, nextPageToken
