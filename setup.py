#!/usr/bin/env python
from setuptools import setup, find_packages
import os
__author__ = 'adamkoziol'

setup(
    name="sipprverse",
    version="0.2.29",
    packages=find_packages(),
    include_package_data=True,
    scripts=[
        os.path.join('sipprverse', 'sippr.py'),
        os.path.join('sipprverse', 'method.py')
    ],
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    description='Object oriented raw read typing software',
    url='https://github.com/OLC-Bioinformatics/sipprverse',
    long_description=open('README.md').read(),
)
