#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
''' blah '''


setup(
    name='fcctx_replicate_agreement',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['fcctx_replicate_agreement'],
    package_dir={'fcctx_replicate_agreement': 'fcctx_replicate_agreement'},
    entry_points={
        'console_scripts': ['repagree = fcctx_replicate_agreement.repagree:main']
    },
    url='https://github.com/bjpop/fcctx_replicate_agreement',
    license='LICENSE',
    description=('blah'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["networkx", "intervaltree"],
)
