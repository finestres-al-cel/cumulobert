#!/usr/bin/env bash

for file in bin/*py planet_orbits/*py planet_orbits/*/*py
do
  echo "yapf --style google $file -i"
  yapf --style google $file -i
done
