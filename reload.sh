#!/bin/bash

rm -f profdpm_1.1.tar.gz
R CMD REMOVE profdpm
R CMD build pkg
R CMD INSTALL profdpm_1.1.tar.gz

