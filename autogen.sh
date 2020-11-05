#!/bin/sh
set -e
autoreconf -if --warnings=all
autoheader
libtoolize
aclocal
automake --add-missing
