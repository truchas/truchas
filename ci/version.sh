#!/bin/bash

# This script extracts the project's current version from git.

set -e
set -x

git describe --tags --dirty > version

cat version
