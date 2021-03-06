#!/bin/sh

# This little script re-generates MANIFEST.in to include files that
# wouldn't otherwise get picked up by distutils. This is relevant when
# building source distributions using "python setup.py sdist", which
# would otherwise leave these files out.

MANIFEST_IN=MANIFEST.in

cat <<EOF > $MANIFEST_IN
# This automatically generated by new_manifest.sh
#
# Note: for files required to build extensions, you probably want to
# place them in the Extension depend list, rather than in this file.
#
include setup.py
include *.txt
include pyamg/*.py
EOF

