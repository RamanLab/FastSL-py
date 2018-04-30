# -*- coding: utf-8 -*-

from __future__ import absolute_import

from setuptools import setup

import fast_sl


def readme():
    with open('README.rst') as file:
        return file.read()


setup(
    name="fast_sl",
    version=fast_sl.__version__,
    install_requires=[
        "cobra",
        "joblib",
        "tqdm",
        "lxml"
    ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    description=("FastSL-py is an efficient algorithm to identify "
                 "synthetic lethal gene/reaction sets in genome-scale "
                 "metabolic models."),
    long_description=readme(),
    url='https://github.com/RamanLab/FastSL-py',
    entry_points={
        'console_scripts': ['fast-sl=fast_sl_cli:main']
    },
    author='',
    author_email='',
    license="GPL v3",
    keywords=("synthetic lethals flux balance analysis linear programming"
              "computational systems biology"),
    packages=['fast_sl'],
    platforms="GNU/Linux, macOS >= 10.7, Microsoft Windows >= 7",
)
