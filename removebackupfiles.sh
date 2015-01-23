#!/bin/sh
echo "recursively removing backup files from"
pwd
find ./ -name 'CMakeLists.txt' -exec rm '{}' \; -print -or -name ".*~" -exec rm {} \; -print
