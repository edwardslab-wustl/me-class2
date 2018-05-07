ME-Class2
=============

Installation
============
Using pip::

    download meclass2-1.0.tar.gz
    pip install meclass2-1.0.tar.gz

For a user-specific install use --user flag::

    pip install --user meclass2-1.0.tar.gz


Required packages
=================
* sklearn >= 0.18
* numpy
* scipy
* seaborn
* matplotlib
* pandas

We  recommend installing anaconda or a similar all in one scientific package which will include all of these.


Quick start using example data
==============================
1) Setup example\_data
-----------------------------------------------
::

    download example_data.tgz
    tar -xzf example_data.tgz
    cd example_data

2) Run interpolation
--------------------
::

    meclass2_interpolation --autosome-only -g refGene_sample.txt -z sample_data_5hmC.bg sample_data_5mC.bg sample_data.expr HRPS_test HRPS

Example run time: ~20-30 minutes

3) Run Classifier
-----------------
::

    meclass2_classifier --num_trees 1001 -t 5hmC . 
    meclass2_classifier --num_trees 1001 -t 5mC . 
    meclass2_classifier --num_trees 1001 -t 5mC_5hmC . 

Example run time: ~3-5 minutes each or faster if using multiple threads (-j option)

4) Report Results
-----------------
::

    meclass2_reporting --plot_results *.pred

Example run time: <1 minute

5) Run clustering
-----------------
::

    meclass2_clustering . cluster --tag 5mC_5hmC --numClusters=4 --lowerPredBound=0.75

Example run time: ~1 minute


Known Issues
============
1) There is a warning thrown in the clustering script about a depreciated function (axisbg). This can be safely ignored.


