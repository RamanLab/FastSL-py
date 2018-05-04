# -*- coding: utf-8 -*-

from __future__ import absolute_import

from setuptools import setup, find_packages

import fastsl


def readme():
    with open('README.rst') as file:
        return file.read()


setup(
    name='fastsl',
    version=fastsl.__version__,
    project_urls={
        'Source': 'https://github.com/RamanLab/FastSL-py'
    },
    install_requires=[
        'cobra',
        'joblib',
        'tqdm',
        'lxml',
        'Click'
    ],
    python_requires='~=3.5',
    setup_requires=['pytest-runner'],
    tests_require=['pytest', 'tox'],
    description=('FastSL-py is an efficient algorithm to identify '
                 'synthetic lethal gene/reaction sets in genome-scale '
                 'metabolic models.'),
    long_description=readme(),
    url='https://github.com/RamanLab/FastSL-py',
    entry_points='''
        [console_scripts]
        fast-sl=fastsl.fastsl_cli:main
    ''',
    author='Synchon Mandal, Karthik Raman',
    author_email='fast-sl@ramanlab.groups.io',
    maintainer='Synchon Mandal',
    maintainer_email='synchonmandal@gmail.com',
    license='GPL v3',
    include_package_data=True,
    zip_safe=False,
    keywords=('synthetic-lethals flux-balance-analysis linear-programming '
              'computational-systems-biology'),
    packages=find_packages(),
    platforms='GNU/Linux, macOS >= 10.7, Microsoft Windows >= 7',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'
    ]
)
