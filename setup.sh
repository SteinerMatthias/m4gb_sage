#!/bin/bash

git clone https://github.com/cr-marcstevens/m4gb.git
cd ./m4gb
autoreconf --install
./configure
make check

