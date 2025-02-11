#!/usr/bin/env bash

output=$1
ref=$2
grep 'TEST' $output > $ref.out
diff -s $ref $ref.out
