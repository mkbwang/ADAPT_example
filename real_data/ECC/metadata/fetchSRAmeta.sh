#!/bin/bash

source ~/.bashrc # need to activate the path to the e-utilities commands
esearch -db sra -query PRJNA752888 | efetch -format runinfo > runinfo_all.csv

