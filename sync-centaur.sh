#!/bin/bash
rsync -av --exclude=build** --exclude=dist** --exclude=ntl** ./ centaur:hap/
