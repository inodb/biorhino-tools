#!/bin/bash
awk '{if ($0 ~ "^@HWI") { print $1"/"substr($2,0,1) } else {print $0} }' $1
