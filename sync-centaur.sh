#!/bin/bash
rsync -av --exclude=build** --exclude=dist** --exclude=ntl** --exclude=cube_** --exclude=callgrind** ./ centaur:hap/
