#!/bin/bash

mkdocs build
base=site
rsync -avz ${base}/* nitweb:httpdocs/contents/nap_docs/
