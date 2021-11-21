#!/usr/bin/zsh

for file in CFIS*.fits
do
	name=$(echo $file | cut -d'.' -f1,2,3,4)
	mkdir -p $name
	cp actual.{params,sex} $name
	cp default.{nnw,conv} $name
	#echo $name.fits
	cd $name
	sex -c actual.sex ../$name.fits
	cd ..
done
