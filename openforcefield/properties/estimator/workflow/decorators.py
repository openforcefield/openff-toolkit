# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Protocol helper decorators.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from enum import Enum


# =============================================================================================
# Custom decorators
# =============================================================================================

class MergeBehaviour(Enum):
    """A enum which describes how attributes should be handled when
    attempting to merge similar protocols.

    Notes
    -----
    Any attributes marked with a merge behavior of `ExactlyEqual`
    must be exactly for two protocols to merge.
    """

    ExactlyEqual = 0,
    SmallestValue = 1,
    GreatestValue = 2,


class BaseProtocolInputObject:
    """A custom decorator used to mark class attributes as either
    a required input, or output, of a protocol.

    Notes
    -----
    This decorator expects the protocol to have a matching private field
    in addition to the public attribute. For example if a protocol has
    an attribute `substance`, by default the protocol must also have a
    `_substance` field.
   """

    def __init__(self, class_attribute):

        documentation = class_attribute.__doc__ or None

        self.__doc__ = documentation

        self.attribute = '_' + class_attribute.__name__
        self.value_type = None

    def __get__(self, instance, owner=None):

        if instance is None:
            # Added in to fix an issue where RTD tries to call this
            # on a class, rather than an instance.
            return self

        if not hasattr(instance, self.attribute):
            raise ValueError('Missing {} attribute.'.format(self.attribute))

        return getattr(instance, self.attribute)

    def __set__(self, instance, value):

        if instance is None:
            raise ValueError('Unexpected ProtocolArgumentDecorator set use case.')

        if not hasattr(instance, self.attribute):
            raise ValueError('Missing {} attribute.'.format(self.attribute))

        from openforcefield.properties.estimator.workflow.protocols import ProtocolPath

        if not isinstance(value, self.value_type) and not isinstance(value, ProtocolPath) and value is not None:

            # Handle the special case where the decimal has been lost on float types...
            if not (self.value_type is float and isinstance(value, int)):

                raise ValueError('The {} attribute can only accept values '
                                 'of type {}'.format(self.attribute, self.value_type))

        setattr(instance, self.attribute, value)


def protocol_input(value_type, merge_behavior=MergeBehaviour.ExactlyEqual):
    """A custom decorator used to mark a protocol attribute as a possible input.

    Examples
    ----------
    To mark an attribute as an input:

    >>> @protocol_input
    >>> def substance(self, Substance):
    >>>     pass
    """

    class ProtocolInputObject(BaseProtocolInputObject):

        def __init__(self, class_attribute):
            super().__init__(class_attribute)

            self.value_type = value_type
            self.merge_behavior = merge_behavior

    return ProtocolInputObject


def protocol_output(value_type):
    """A custom decorator used to mark a protocol attribute as
    an output of the protocol.

    Examples
    ----------
    To mark a property as an output:

    >>> @protocol_output(str)
    >>> def positions(self):
    >>>     pass
    """

    class ProtocolOutputObject(BaseProtocolInputObject):

        def __init__(self, class_attribute):
            super().__init__(class_attribute)

            self.value_type = value_type

    return ProtocolOutputObject