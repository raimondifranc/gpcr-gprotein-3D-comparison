#!/bin/bash
 
for mmaker_nogui in `ls -1 *_mmaker_nogui`
  do
    chimera --nogui cmd:$mmaker_nogui > $mmaker_nogui\_output.txt
  done
