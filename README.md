OCaml bindings to the [PARI](https://pari.math.u-bordeaux.fr/) library -- a cross platform and open-source computer algebra system designed for fast computations in number theory.

The API is documented in the [Users' Guide to the PARI library](https://pari.math.u-bordeaux.fr/pub/pari/manuals/2.15.1/libpari.pdf).

The best way to get the dependencies is through [Nix](https://nixos.org/).
Without Nix, you will need opam 2.0 and GMP installed on your system.
To generate the stubs you will need LLVM and libclang.

1. With OPAM: `opam update && opam switch create . -y --deps-only` (or with Nix: `nix develop -c $SHELL`)
2. `dune build`

To regenerate the stubs run `make gen-stubs` (after `nix-shell -p llvm libclang`).

Execute a code sample: `dune exec examples/pohlig_hellman.exe`, or view some of the available Jupyter notebooks.
