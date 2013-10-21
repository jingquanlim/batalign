AUTOMAKE_OPTIONS = foreign
#SUBDIRS = src scripts 
SUBDIRS = src  
copy:
	if [ -d "bin" ]; then echo bin exists; else mkdir bin; fi
	cp scripts/strip.pl bin 
	cp scripts/b2c.pl bin 
	cp scripts/build_complement bin 
	cp scripts/filter.pl bin
	cp scripts/build_index bin 
	cp scripts/ann2loc.pl bin 
	cp scripts/build_location.pl bin 
	cp scripts/build_revcmp bin 
	cp scripts/complement.pl bin
	cp src/bwtformatdb bin 
	cp src/reverse bin 
	cp src/batman bin 
	cp src/batdecode bin 
	cp src/filter bin
