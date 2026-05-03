Installation and docs build
===========================

1. Create and activate a virtual environment (recommended):

.. code-block:: bash

   python3 -m venv .venv
   source .venv/bin/activate

   If you already use a named environment (for example ``Env_TG``), activate that instead.

2. Install the project in editable mode with documentation extras:

.. code-block:: bash

   pip install -e ".[docs]"

3. Build HTML documentation with Sphinx:

.. code-block:: bash

   sphinx-build -b html docs docs/_build/html

   For a clean rebuild from scratch:

.. code-block:: bash

   sphinx-build -E -a -b html docs docs/_build/html

4. Open ``docs/_build/html/index.html`` in a browser.

The ``docs/_build/`` directory is excluded from version control; CI rebuilds the site on pushes to ``main`` (see ``.github/workflows/docs.yml``).
