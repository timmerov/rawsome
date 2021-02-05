# rawsome

## info

reads a raw .CR2 file from a canon rebel ti 7 using libraw.
does some processing.
writes a png file.
gory details are in main.cc.

## build

one time only:

$ ./gen-unix64.sh

to build:

$ ./make-unix64.sh

build for windows is not supported.

## run

$ ./bin/rawsome_d [in.CR2 out.png]
