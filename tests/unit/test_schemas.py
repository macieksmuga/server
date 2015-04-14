"""
Tests for auto generated schemas and conversion to and from JSON.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import string
import random
import unittest

import avro.schema

import ga4gh.protocol as protocol
import tests.utils as utils


def randomString():
    """
    Returns a randomly generated short string.
    """
    randInt = random.randint(0, 10)
    randStr = ''.join(random.choice(
        string.ascii_letters) for _ in range(randInt))
    return randStr


class SchemaTest(unittest.TestCase):
    """
    Superclass of schema tests.
    """
    typicalValueMap = {
        "string": "string value",
        "int": 1000,
        "long": 10000,
        "boolean": True,
        "double": 0.125,
        "float": 0.25
    }

    instanceGenerator = utils.InstanceGenerator()

    def getAvroSchema(self, cls, fieldName):
        """
        Returns the avro schema for the specified field.
        """
        field = None
        for fld in cls.schema.fields:
            if fld.name == fieldName:
                field = fld
        return field

    def getInvalidValue(self, cls, fieldName):
        """
        Returns a value that should trigger a schema validation failure.
        """
        fieldType = self.getAvroSchema(cls, fieldName).type
        if isinstance(fieldType, avro.schema.UnionSchema):
            types = list(t.type for t in fieldType.schemas)
            val = self.instanceGenerator.generateInvalidateTypeValue(*types)
        else:
            val = self.instanceGenerator.generateInvalidateTypeValue(fieldType)
        return val

    def getTypicalValue(self, cls, fieldName):
        """
        Returns a typical value for the specified field on the specified
        Protocol class.
        """
        # We make some simplifying assumptions about how the schema is
        # structured which fits the way the GA4GH protocol is currently
        # designed but may break in the future. We try to at least flag
        # this fact here.
        err = "Schema structure assumptions violated"
        field = self.getAvroSchema(cls, fieldName)
        typ = field.type
        if isinstance(typ, avro.schema.UnionSchema):
            t0 = typ.schemas[0]
            if (isinstance(t0, avro.schema.PrimitiveSchema) and
                    t0.type == "null"):
                typ = typ.schemas[1]
            else:
                raise Exception(err)
        ret = None
        if isinstance(typ, avro.schema.MapSchema):
            ret = {"key": ["value1", "value2"]}
            if not isinstance(typ.values, avro.schema.ArraySchema):
                raise Exception(err)
        elif isinstance(typ, avro.schema.ArraySchema):
            if cls.isEmbeddedType(field.name):
                embeddedClass = cls.getEmbeddedType(field.name)
                ret = [self.getTypicalInstance(embeddedClass)]
            else:
                ret = [self.typicalValueMap[typ.items.type]]
        elif isinstance(typ, avro.schema.EnumSchema):
            ret = typ.symbols[0]
        elif isinstance(typ, avro.schema.RecordSchema):
            self.assertTrue(cls.isEmbeddedType(fieldName))
            embeddedClass = cls.getEmbeddedType(fieldName)
            ret = self.getTypicalInstance(embeddedClass)
        elif typ.type in self.typicalValueMap:
            ret = self.typicalValueMap[typ.type]
        else:
            raise Exception("schema assumptions violated")
        return ret

    def getTypicalInstance(self, cls):
        """
        Returns a typical instance of the specified protocol class.
        """
        instance = cls()
        for field in cls.schema.fields:
            setattr(instance, field.name,
                    self.getTypicalValue(cls, field.name))
        return instance

    def getRandomInstance(self, cls):
        """
        Returns an instance of the specified class with randomly generated
        values conforming to the schema.
        """
        return self.instanceGenerator.generateInstance(cls)

    def setRequiredValues(self, instance):
        """
        Sets the required values in the specified instance to typical values.
        """
        for key in instance.__slots__:
            if key in instance.requiredFields:
                value = self.getTypicalValue(type(instance), key)
                setattr(instance, key, value)

    def getDefaultInstance(self, cls):
        """
        Returns a new instance with the required values set.
        """
        instance = cls()
        self.setRequiredValues(instance)
        return instance


class EqualityTest(SchemaTest):
    """
    Tests equality is correctly calculated for different protocol elements.
    """
    def verifyEqualityOperations(self, i1, i2):
        self.assertEqual(i1, i1)
        self.assertTrue(i1 == i1)
        self.assertFalse(i1 != i1)
        self.assertEqual(i1, i2)
        self.assertTrue(i1 == i2)
        self.assertFalse(i1 != i2)
        for val in [None, {}, [], object, ""]:
            self.assertFalse(i1 == val)
            self.assertTrue(i1 != val)
            self.assertFalse(val == i1)
            self.assertTrue(val != i1)
        # Now change an attribute on one and check if equality fails.
        for field in i1.schema.fields:
            setattr(i1, field.name, "value unique to i1")
            self.assertFalse(i1 == i2)
            self.assertTrue(i1 != i2)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testSameClasses(self):
        factories = [self.getDefaultInstance, self.getTypicalInstance,
                     self.getRandomInstance]
        for cls in protocol.getProtocolClasses():
            for factory in factories:
                i1 = factory(cls)
                i2 = cls.fromJsonDict(i1.toJsonDict())
                self.verifyEqualityOperations(i1, i2)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testDifferentValues(self):
        def factory(cls):
            return cls()
        factories = [factory, self.getTypicalInstance, self.getDefaultInstance,
                     self.getRandomInstance]
        classes = list(protocol.getProtocolClasses())
        c1 = classes[0]
        for c2 in classes[1:]:
            for factory in factories:
                i1 = factory(c1)
                i2 = factory(c2)
                self.assertFalse(i1 == i2)
                self.assertTrue(i1 != i2)

    def testDifferentLengthArrays(self):
        i1 = self.getTypicalInstance(protocol.CallSet)
        i2 = protocol.CallSet.fromJsonDict(i1.toJsonDict())
        i2.variantSetIds.append("extra")
        self.assertFalse(i1 == i2)


class SerialisationTest(SchemaTest):
    """
    Tests the serialisation and deserialisation code for the schema classes
    """
    def validateClasses(self, factory):
        for cls in protocol.getProtocolClasses():
            instance = factory(cls)
            jsonStr = instance.toJsonString()
            otherInstance = cls.fromJsonString(jsonStr)
            self.assertEqual(instance, otherInstance)

            jsonDict = instance.toJsonDict()
            otherInstance = cls.fromJsonDict(jsonDict)
            self.assertEqual(instance, otherInstance)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testSerialiseDefaultValues(self):
        self.validateClasses(self.getDefaultInstance)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testSerialiseTypicalValues(self):
        self.validateClasses(self.getTypicalInstance)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testSerialiseRandomValues(self):
        self.validateClasses(self.getRandomInstance)


class ValidatorTest(SchemaTest):
    """
    Tests the validator to see if it will correctly identify instances
    that do not match the schema and also that it correctly identifies
    instances that do match the schema
    """
    def validateClasses(self, factory):
        for cls in protocol.getProtocolClasses():
            instance = factory(cls)
            jsonDict = instance.toJsonDict()
            self.assertTrue(cls.validate(jsonDict))

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testValidateDefaultValues(self):
        self.validateClasses(self.getDefaultInstance)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testValidateTypicalValues(self):
        self.validateClasses(self.getTypicalInstance)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testValidateRandomValues(self):
        self.validateClasses(self.getRandomInstance)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testValidateBadValues(self):
        for cls in protocol.getProtocolClasses():
            instance = self.getTypicalInstance(cls)
            jsonDict = instance.toJsonDict()
            self.assertFalse(cls.validate(None))
            self.assertFalse(cls.validate([]))
            self.assertFalse(cls.validate(1))
            # setting values to bad values should be invalid
            for key in jsonDict.keys():
                dct = dict(jsonDict)
                dct[key] = self.getInvalidValue(cls, key)
                self.assertFalse(cls.validate(dct))
            for c in utils.powerset(jsonDict.keys(), 10):
                if len(c) > 0:
                    dct = dict(jsonDict)
                    for f in c:
                        dct[f] = self.getInvalidValue(cls, f)
                    self.assertFalse(cls.validate(dct))


class GetProtocolClassesTest(SchemaTest):
    """
    Tests the protocol.getProtocolClasses() function to ensure it
    works correctly.
    """
    def testAllClasses(self):
        classes = protocol.getProtocolClasses()
        assert len(classes) > 0
        for class_ in classes:
            self.assertTrue(issubclass(class_, protocol.ProtocolElement))

    def testRequestAndResponseClasses(self):
        requestClasses = protocol.getProtocolClasses(protocol.SearchRequest)
        responseClasses = protocol.getProtocolClasses(protocol.SearchResponse)
        self.assertEqual(len(requestClasses), len(responseClasses))
        self.assertGreater(len(requestClasses), 0)
        for class_ in requestClasses:
            self.assertTrue(issubclass(class_, protocol.SearchRequest))
        for class_ in responseClasses:
            self.assertTrue(issubclass(class_, protocol.SearchResponse))
            valueListName = class_.getValueListName()
            self.assertGreater(len(valueListName), 0)


class SearchResponseBuilderTest(SchemaTest):
    """
    Tests the SearchResponseBuilder class to ensure that it behaves
    correctly.
    """
    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testIntegrity(self):
        # Verifies that the values we put in are exactly what we get
        # back across all subclasses of SearchResponse
        for class_ in protocol.getProtocolClasses(protocol.SearchResponse):
            instances = [
                self.getTypicalInstance(class_),
                self.getRandomInstance(class_)]
            for instance in instances:
                valueList = getattr(instance, class_.getValueListName())
                builder = protocol.SearchResponseBuilder(
                    class_, len(valueList), 2**32)
                for value in valueList:
                    builder.addValue(value)
                builder.setNextPageToken(instance.nextPageToken)
                otherInstance = class_.fromJsonString(builder.getJsonString())
                self.assertEqual(instance,  otherInstance)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testPageSizeOverflow(self):
        # Verifies that the page size behaviour is correct when we keep
        # filling after full is True.
        responseClass = protocol.SearchVariantsResponse
        valueClass = protocol.Variant
        for pageSize in range(1, 10):
            builder = protocol.SearchResponseBuilder(
                responseClass, pageSize, 2**32)
            self.assertEqual(builder.getPageSize(), pageSize)
            self.assertFalse(builder.isFull())
            for listLength in range(1, 2 * pageSize):
                builder.addValue(self.getTypicalInstance(valueClass))
                instance = responseClass.fromJsonString(
                    builder.getJsonString())
                valueList = getattr(
                    instance, responseClass.getValueListName())
                self.assertEqual(len(valueList), listLength)
                if listLength < pageSize:
                    self.assertFalse(builder.isFull())
                else:
                    self.assertTrue(builder.isFull())

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testPageSizeExactFill(self):
        responseClass = protocol.SearchVariantsResponse
        valueClass = protocol.Variant
        for pageSize in range(1, 10):
            builder = protocol.SearchResponseBuilder(
                responseClass, pageSize, 2**32)
            self.assertEqual(builder.getPageSize(), pageSize)
            while not builder.isFull():
                builder.addValue(self.getTypicalInstance(valueClass))
            instance = responseClass.fromJsonString(builder.getJsonString())
            valueList = getattr(instance, responseClass.getValueListName())
            self.assertEqual(len(valueList), pageSize)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testMaxResponseLengthOverridesPageSize(self):
        responseClass = protocol.SearchVariantsResponse
        valueClass = protocol.Variant
        typicalValue = self.getTypicalInstance(valueClass)
        typicalValueLength = len(typicalValue.toJsonString())
        for numValues in range(1, 10):
            maxResponseLength = numValues * typicalValueLength
            builder = protocol.SearchResponseBuilder(
                responseClass, 1000, maxResponseLength)
            self.assertEqual(
                maxResponseLength, builder.getMaxResponseLength())
            while not builder.isFull():
                builder.addValue(typicalValue)
            instance = responseClass.fromJsonString(builder.getJsonString())
            valueList = getattr(instance, responseClass.getValueListName())
            self.assertEqual(len(valueList), numValues)

    def testNextPageToken(self):
        responseClass = protocol.SearchVariantsResponse
        builder = protocol.SearchResponseBuilder(
            responseClass, 100, 2**32)
        # If not set, pageToken should be None
        self.assertIsNone(builder.getNextPageToken())
        instance = responseClass.fromJsonString(builder.getJsonString())
        self.assertIsNone(instance.nextPageToken)
        # page tokens can be None or any string.
        for nextPageToken in [None, "", "string"]:
            builder.setNextPageToken(nextPageToken)
            self.assertEqual(nextPageToken, builder.getNextPageToken())
            instance = responseClass.fromJsonString(builder.getJsonString())
            self.assertEqual(nextPageToken, instance.nextPageToken)
