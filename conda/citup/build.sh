#!/bin/bash

cd src
CPLUS_INCLUDE_PATH=$PREFIX/include LIBRARY_PATH=$PREFIX/lib BOOST_DIRECTORY=$PREFIX/src/boost/ make install

cd ..
python setup.py install
