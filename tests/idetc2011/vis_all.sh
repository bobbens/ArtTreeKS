#!/usr/bin/env bash

DIRNAME="res_index_middle_thumb_3_333"

for FILE in $DIRNAME/*.lua; do
   lua vis.lua test_run_base.lua "$FILE"
done


