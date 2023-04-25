#!/bin/bash
set -eo pipefail

while true
do
  ps -o pid,user,%mem,command ax | sort -b -k3 -r  | grep python >> /app/fireatlas_nrt/memory.log
  sleep 5
done
