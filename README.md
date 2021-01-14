HotspotFinder
================

This repository contains HotspotFinder, a tool for the identification and annotation of hotspots across the genome.


License
-------

HotspotFinder is available to the general public subject to certain conditions described in its [LICENSE](other_file.md)


Installation
------------

HotspotFinder is implemented in Python 3.6.

You can obtain the latest code from the repository and install it with pip::

    $ git clone git@bitbucket.org:bbglab/hotspotfinder.git
    $ cd hotspotfinder
    $ pip install .


The first time that you run HotspotFinder it will download data from our servers. By default the
downloaded datasets go to ``~/.bgdata``. If you want to move these datasets to another folder you have to define the
system environment variable BGDATA_LOCAL with an export command.

The following command will show you the help::

    $ hotspotfinder --help

