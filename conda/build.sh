#!/bin/bash
set -e
set -u
mkdir -p $PREFIX/bin
install -m 755 bin/pcne $PREFIX/bin/
install -m 644 bin/PCNE.R $PREFIX/bin/
echo "Installation script finished."
echo "Files installed in $PREFIX/bin:"
ls -l $PREFIX/bin