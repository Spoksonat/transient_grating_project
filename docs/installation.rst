Installation and docs build
===========================

1. Activate the project environment:

.. code-block:: bash

   source Env_TG/bin/activate

2. Install the project and documentation dependencies:

.. code-block:: bash

   pip install -e ".[docs]"

3. Build HTML documentation with Sphinx:

.. code-block:: bash

   sphinx-build -b html docs docs/_build/html

4. Open the generated documentation at:

.. code-block:: text

   docs/_build/html/index.html
