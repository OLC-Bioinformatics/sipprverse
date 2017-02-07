#!/usr/bin/env python
from customtargets import Custom
__author__ = 'adamkoziol'

class SixteenS(Custom):

    def reporter(self):
        pass

    def __init__(self, inputobject, analysistype):
        Custom.__init__(self, inputobject, analysistype)
        self.targetpath = inputobject.customtargetpath if inputobject.customtargetpath else inputobject.targetpath \
            + 'sixteenS/'
