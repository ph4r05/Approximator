#!/bin/bash
rsync -av --exclude="*.bin" --exclude=ntl** --exclude=cube_** --exclude=callgrind** ./ aisa.fi.muni.cz:hap/
