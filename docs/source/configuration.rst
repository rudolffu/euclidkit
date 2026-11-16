Credentials and Configuration
=============================

Commands that query the Euclid archive log in with your ESA Cosmos account.
This page shows how to give ``euclidkit`` your credentials safely.

Create a credentials file
-------------------------

1. Create a private file in your home directory:

   .. code-block:: bash

      mkdir -p ~/.euclidkit
      touch ~/.euclidkit/.cred.txt
      chmod 600 ~/.euclidkit/.cred.txt

2. Open ``~/.euclidkit/.cred.txt`` in a text editor and add two lines:

   .. code-block:: text

      your_cosmos_username
      your_cosmos_password

Do not write the password with a shell command such as ``echo``, because the
shell records it in its history.

``~/.euclidkit/.cred.txt`` is the default location. If your file is there, you
do not need any further configuration.

Use a different credentials file
--------------------------------

You can point to another file in three ways. ``euclidkit`` uses the first one
that it finds:

1. The ``--credentials`` (``-c``) option of a command, or the
   ``credentials_file`` argument of ``EuclidArchive.login()``:

   .. code-block:: bash

      euclidkit crossmatch -i sources.fits -o matches.fits -c /path/to/cred.txt

2. The ``EUCLIDKIT_CREDENTIALS_FILE`` or ``EUCLIDKIT_CREDENTIALS`` environment
   variable:

   .. code-block:: bash

      export EUCLIDKIT_CREDENTIALS=/path/to/cred.txt

3. The ``data.credentials_file`` setting in your user configuration file,
   described below.

On ESA Datalabs, a common location is ``/media/user/cred.txt``.

Log in without a file (Python only)
-----------------------------------

If no credentials file exists, ``EuclidArchive.login()`` tries these sources
in order:

1. The ``EUCLID_USER`` and ``EUCLID_PASSWORD`` environment variables.
2. The operating-system keyring, under the service name ``euclidkit``, if the
   ``keyring`` package is installed and you pass ``user``.
3. A password prompt, if you pass ``user`` and ``prompt=True``.

.. code-block:: python

   from euclidkit.core.data_access import EuclidArchive

   archive = EuclidArchive()
   archive.login(user="your_cosmos_username", prompt=True)

User configuration file
-----------------------

The configuration file is optional. To create one with the default settings:

.. code-block:: bash

   euclidkit init-config

This writes ``~/.euclidkit/euclidkit_config.yaml``, which ``euclidkit`` reads
at startup. Use ``--output`` to write the template somewhere else.

To change the credentials path, edit the ``data`` section:

.. code-block:: yaml

   data:
     credentials_file: /home/<user>/.euclidkit/.cred.txt

.. note::

   The template also contains ``analysis`` and ``output`` sections, and
   ``init-config --template`` offers ``advanced`` and ``pipeline`` variants.
   The commands described in this documentation do not read those settings.

Next step
---------

Run a complete example: :doc:`quickstart`.
