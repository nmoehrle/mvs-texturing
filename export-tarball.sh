#!/bin/bash
DEST=/tmp/texrecon.tar.gz
echo exporting the currently checked out repo revision to $DEST ...
git archive HEAD -o $DEST --prefix=texrecon/

