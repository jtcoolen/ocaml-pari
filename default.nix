{ system ? builtins.currentSystem, pkgs ? import <nixpkgs> { inherit system; }
}:
let
  pkgs = import (fetchTarball {
    url =
      "https://github.com/NixOS/nixpkgs/archive/0034715d703ee02226c5708c263583a82539f2ea.tar.gz";
    sha256 = "02jg1qwr3f8xagfwv5j8asyajwwgl1zmj6xaj5yfxd30l33d62y5";
  }) { inherit system; };
in let
  ocamlPackages = pkgs.ocaml-ng.ocamlPackages_5_2; #.overrideScope' (self: super: {
#    ocaml = super.ocaml.override { flambdaSupport = true; }; });

  ctypes-foreign = pkgs.lib.overrideDerivation (ocamlPackages.ctypes-foreign)
    (old: {
      NIX_CFLAGS_COMPILE = "-Wno-error=incompatible-function-pointer-types";
    });
in pkgs.clangStdenv.mkDerivation {
  name = "ocaml_pari";
  nativeBuildInputs = (with ocamlPackages; [
    ocaml
    findlib
    dune_3
    merlin
    utop
#    ppx_cstubs
    ctypes
    containers
    mdx
    odoc
    ocaml-lsp
    hacl-star
    hex
    iter
    qcheck
    memtrace
    ocamlformat
  ]) ++ (with pkgs; [
    bison
    gnumake
    pkg-config
    #gcc
    gmp
    gmpxx
    #libcxx
    # for linux
    #glibc
    #glibc.static
    clang # for i-shiny objects
    #nixfmt-rfc-style
    ocamlformat
    llvmPackages.clang-unwrapped
    llvmPackages.llvm
  ]);
  dontDetectOcamlConflicts = true;
  buildInputs =
    (with ocamlPackages; [ core ppx_expect ctypes odoc ctypes-foreign ]);
  #LD_LIBRARY_PATH = "${pkgs.glibc}/lib:${pkgs.glibc.static}/lib";
  NIX_LDFLAGS = [ "-lc -lm" ]
    ++ pkgs.lib.optionals pkgs.stdenv.isDarwin [ "-no_compact_unwind" ];
  NIX_CFLAGS_COMPILE =
    # Silence errors (-Werror) for unsupported flags on MacOS.
    pkgs.lib.optionals pkgs.stdenv.isDarwin [
      "-Wno-unused-command-line-argument"
      "-Wmacro-redefined"
      "-Wincompatible-pointer-types-discards-qualifiers"
    ];
}
