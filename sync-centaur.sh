#!/bin/bash
rsync -av --exclude=build** --exclude=ntl** --exclude=cube_** --exclude=callgrind** ./ centaur:hap/
