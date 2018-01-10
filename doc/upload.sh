#!/bin/bash

base=_build/html
rsync -avz ${base}/* nitweb:httpdocs/contents/pmd_usage/
