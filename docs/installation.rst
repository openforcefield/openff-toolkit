.. _installation:

Installation
************

Installing via `conda`
======================

The simplest way to install the Open Forcefield Toolkit is via the `conda <http://www.continuum.io/blog/conda>`_  package manager.
Packages are provided on the `omnia Anaconda Cloud channel <http://anaconda.org/omnia>`_ for Linux, OS X, and Win platforms.
The `openforcefield Anaconda Cloud page <https://anaconda.org/omnia/openforcefield>`_ has useful instructions and `download statistics <https://anaconda.org/omnia/openforcefield/files>`_.

If you are using the `anaconda <https://www.continuum.io/downloads/>`_ scientific Python distribution, you already have the ``conda`` package manager installed.
If not, the quickest way to get started is to install the `miniconda <http://conda.pydata.org/miniconda.html>`_ distribution, a lightweight minimal installation of Anaconda Python.

On ``linux``, you can install the Python 3 version into ``$HOME/miniconda3`` with (on ``bash`` systems):

.. code-block:: bash

   $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   $ bash ./Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
   $ source ~/miniconda3/etc/profile.d/conda.sh
   $ conda activate base


On ``osx``, you want to use the ``osx`` binary

.. code-block:: bash

   $ curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O
   $ bash ./Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda3
   $ source ~/miniconda3/etc/profile.d/conda.sh
   $ conda activate base


You may want to add the new ``source ~/miniconda3/etc/profile.d/conda.sh`` line to your ``~/.bashrc`` file to ensure Anaconda Python can enabled in subsequent terminal sessions.
``conda activate base`` will need to be run in each subsequent terminal session to return to the environment where the toolkit will be installed.


Note that ``openforcefield`` will be installed into this local Python installation, so that you will not need to worry about disrupting existing Python installations.

.. note:: Installation via the conda package manager is the preferred method since all dependencies are automatically fetched and installed for you.

|

Required dependencies
=======================

The ``openforcefield`` toolkit makes use of the `Omnia <http://www.omnia.md>`_ and `Conda Forge <https://conda-forge.org/>`_ free and open source community package repositories:

.. code-block:: bash

   $ conda config --add channels omnia --add channels conda-forge
   $ conda update --all

This only needs to be done once.

.. note ::

   If automation is required, provide the ``--yes`` argument to ``conda update`` and ``conda install`` comamnds.
   More information on the ``conda`` command-line API can be found in the `conda online documentation <https://conda.io/docs/commands.html>`_.

|

Release build
-------------

You can install the latest stable release build of ``openforcefield`` via the ``conda`` package with

.. code-block:: none

   $ conda config --add channels omnia --add channels conda-forge
   $ conda install openforcefield

This version is recommended for all users not actively developing new forcefield parameterization algorithms.

.. note:: The conda package manager will install dependencies from binary packages automatically, including difficult-to-install packages such as OpenMM, numpy, and scipy. This is really the easiest way to get started.

|

Upgrading your installation
---------------------------

To update an earlier ``conda`` installation of ``openforcefield`` to the latest release version, you can use ``conda update``:

.. code-block:: bash

   $ conda update openforcefield

|

Optional dependencies
---------------------

This toolkit can optionally make use of the `OpenEye toolkit <https://www.eyesopen.com/toolkit-development>`_ if the user has a license key installed.
Academic laboratories intending to release results into the public domain can `obtain a free license key <https://www.eyesopen.com/licensing-philosophy>`_, while other users (including academics intending to use the software for purposes of generating protected intellectual property) must `pay to obtain a license <https://www.eyesopen.com/pricing>`_.

To install the OpenEye toolkits (provided you have a valid license file):

.. code-block:: none

   $ conda install --yes -c openeye openeye-toolkits

No essential ``openforcefield`` release capabilities *require* the OpenEye toolkit, but the Open Force Field developers make use of it in parameterizing new open source force fields.
It is known that there are certain differences in toolkit behavior between RDKit and OpenEye when reading a small fraction of molecules, and we encourage you to report any unexpected behavior that may be caused by toolkit differences to our `issue tracker <https://github.com/openforcefield/openforcefield/issues>`_.


