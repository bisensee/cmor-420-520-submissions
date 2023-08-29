#!/bin/bash

read -p "Project name: " project
if [ ! -d "./.build" ]; then
	mkdir ./.build
fi
pdflatex $project.tex
mv $project*.aux ./.build
