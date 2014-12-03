AUTOMAKE_OPTIONS = foreign
SUBDIRS = src scripts 
copy:
	if [ -d "bin" ]; then echo bin exists; else mkdir bin; fi
	cp scripts/strip.pl bin 
	cp scripts/b2c.pl bin 
	cp scripts/build_complement bin 
	cp scripts/filter.pl bin
	cp scripts/build_indexX bin 
	cp scripts/ann2loc.pl bin 
	cp scripts/build_location.pl bin 
	cp scripts/build_revcmp bin 
	cp scripts/complement.pl bin
	cp src/bwtformatdb bin 
	cp src/reverse bin 
	cp src/penguin bin 
	cp src/batalign bin
	cp src/batalignindex bin
	cp src/filter bin
	cp src/penguin-s bin 
	cp src/penguin-l bin 
	cp src/batalign bin 
