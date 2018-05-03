#!/usr/bin/env sh

cd tmp/
if [ -e "igraph_1.1.2.tar.gz" ]
then
  tar -xzf igraph*.tar.gz
  cp ../igraph_patched_file/foreign-graphml.c igraph/src/
  R CMD INSTALL ./igraph
else
  echo "An igraph version newer than 1.1.2 has been released. Please fall back to standard R install."
  exit 1
fi

