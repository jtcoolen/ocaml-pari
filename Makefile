INCLUDEDIR=`llvm-config --includedir`
LIBDIR=`llvm-config --ldflags`

generatebindings: stubs_generation/generatebindings.cpp
	g++ -std=c++14 -Wall -Wextra -Wpointer-arith -I$(INCLUDEDIR) $(LIBDIR) -o $@ $< -lclang

gen-stubs: generatebindings
	./generatebindings -f src/libpari/src/headers/pari.h && mv *.ml src/ && dune build @fmt --auto-promote

clean:
	rm generatebindings
