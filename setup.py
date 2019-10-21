#!/usr/bin/env python

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='meclass2',
        version='0.2.0',
        description="Software to find 5mC and 5hmC associated differentially expressed genes",
        long_description=readme(),
        url='http://github.com/jedwards-wustl/me-classv2',
        author='Chris Schlosberg and John Edwards',
        author_email="jredwards@wustl.edu",
        license='GPLv3',
        packages=['meclass2'],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Operating System :: POSIX',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.6',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.2',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
        ],
        install_requires=[
            'numpy',
            'scipy',
            'scikit-learn>=0.18',
            'matplotlib',
            'seaborn',
            'pandas',
        ],
        #python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
        entry_points = {
            "console_scripts": [
                'meclass2_classifier = meclass2.classifier_main:main',
                'meclass2_interpolation = meclass2.interpolation_main:main',
                'meclass2_clustering = meclass2.clustering_main:main',
                'meclass2_reporting = meclass2.reporting_main:main',
                ],
            },
        include_package_data=True,
        zip_safe=False)


#packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
