Configuration
=============

Credentials file (recommended)
------------------------------

Create a private credentials file and lock file permissions:

.. code-block:: bash

   mkdir -p ~/.euclidkit
   touch ~/.euclidkit/.cred.txt
   chmod 600 ~/.euclidkit/.cred.txt

Edit ``~/.euclidkit/.cred.txt`` manually and store:

1. COSMOS username
2. COSMOS password

User config file
----------------

Generate a config template:

.. code-block:: bash

   euclidkit init-config --output ~/.euclidkit/euclidkit_config.yaml --template basic

Set the credentials path in ``~/.euclidkit/euclidkit_config.yaml``:

.. code-block:: yaml

   data:
     credentials_file: /home/<user>/.euclidkit/.cred.txt

Environment override
--------------------

You can override the configured credential path with:

.. code-block:: bash

   export EUCLIDKIT_CREDENTIALS=/home/<user>/.euclidkit/.cred.txt

