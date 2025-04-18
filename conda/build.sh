#!/bin/bash
# This script is executed by conda-build when building the package.

# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u

# Create the directory structure within the conda environment where binaries are stored.
# $PREFIX is a variable provided by conda-build, pointing to the root of the installation environment.
mkdir -p $PREFIX/bin

# Copy the executable scripts from the source directory (your project's bin/)
# into the environment's bin directory ($PREFIX/bin).
# 'install' is preferred over 'cp' as it can set permissions directly.
# -m 755 sets read/write/execute permissions for owner, read/execute for group/others (standard for executables)
install -m 755 bin/pcne $PREFIX/bin/

# -m 644 sets read/write for owner, read-only for group/others (standard for non-executable scripts/data)
install -m 644 bin/PCNE.R $PREFIX/bin/

echo "Installation script finished."
echo "Files installed in $PREFIX/bin:"
ls -l $PREFIX/bin