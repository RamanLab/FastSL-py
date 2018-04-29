# -*- coding: utf-8 -*-

from __future__ import absolute_import

from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name="fast_sl",
    version='0.1.0',
    install_requires=[
        "cobra",
        "joblib",
        "tqdm",
        "lxml"
    ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    test_suite="fast_sl.test.suite",
    description=("FastSL-py is an efficient algorithm to identify "
                 "synthetic lethal gene/reaction sets in genome-scale "
                 "metabolic models."),
    long_description=readme(),
    url='https://github.com/RamanLab/FastSL-py',
    entry_points={
        'console_scripts': ['fast-sl=fast_sl.fast_sl:main']
    },
    author='',
    author_email='',
    license="LGPL/GPL v2+",
    keywords=("synthetic lethals flux balance analysis linear programming"
              "computational systems biology"),
    packages=['fast_sl'],
    platforms="GNU/Linux, macOS >= 10.7, Microsoft Windows >= 7",
    zip_safe=False)
