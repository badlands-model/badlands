Documentation is built with Sphinx. To install:    

    pip install -r requirements.txt

You should install Sphinx in the same virtualenv or Python environment as Badlands.

To build documentation:

    READTHEDOCS=True make clean html

Documentation is hosted on badlands.readthedocs.io. RTD's builder can't build
the FORTRAN components of Badlands, so we stub them out with MagicMock, per
the instructions at http://docs.readthedocs.io/en/latest/faq.html#i-get-import-errors-on-libraries-that-depend-on-c-modules
