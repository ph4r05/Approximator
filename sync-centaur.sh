#!/bin/bash
rsync -av --exclude=build** --exclude=dist** ./ centaur:hap/
