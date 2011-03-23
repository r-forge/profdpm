#!/bin/bash

R CMD REMOVE profdpm
R CMD build pkg
R CMD INSTALL profdpm_2.0.tar.gz

