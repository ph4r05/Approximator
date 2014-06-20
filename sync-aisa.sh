#!/bin/bash
rsync -av --exclude=build** --exclude=dist** --exclude=ntl** --exclude=cube_** --exclude=callgrind** ./ aisa.fi.muni.cz:hap/
