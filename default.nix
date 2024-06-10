{ system ? builtins.currentSystem, pkgs ? import <nixpkgs> { inherit system; }
}:
let
  pkgs = import (fetchTarball {
    url =
      "https://github.com/NixOS/nixpkgs/archive/90e4a08f48e2b3cd9f827dca60300c93b1fea7a3.tar.gz";
    sha256 = "1hr9xic73ciycmgxcl0d446xk1lfbybb6p2q1gkw1zd4rxnbzkih";
  }) { inherit system; };
in let
  ocamlPackages = pkgs.ocaml-ng.ocamlPackages_5_1;
  #.overrideScope' (self: super: {
  #  ocaml = super.ocaml.override { flambdaSupport = true; };

  ctypes-foreign = pkgs.lib.overrideDerivation (ocamlPackages.ctypes-foreign)
    (old: {
      NIX_CFLAGS_COMPILE = "-Wno-error=incompatible-function-pointer-types";
    });
in pkgs.fastStdenv.mkDerivation {
  name = "ocaml_pari";
  nativeBuildInputs = (with ocamlPackages; [
    ocaml
    findlib
    dune_3
    merlin
    utop
    ppx_cstubs
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
    gcc
    gmp
    gmpxx
    libcxx
    # for linux
    #glibc
    #glibc.static
    clang # for i-shiny objects
    nixfmt-rfc-style
    ocamlformat
  ]);
  dontDetectOcamlConflicts = true;
  buildInputs =
    (with ocamlPackages; [ core ppx_expect ctypes odoc ctypes-foreign ]);
  #LD_LIBRARY_PATH = "${pkgs.glibc}/lib:${pkgs.glibc.static}/lib";
  NIX_LDFLAGS = [ "-lc -lm" ]
    ++ pkgs.lib.optionals pkgs.stdenv.isDarwin [ "-Wl,-no_compact_unwind" ];
  NIX_CFLAGS_COMPILE =
    # Silence errors (-Werror) for unsupported flags on MacOS.
    pkgs.lib.optionals pkgs.stdenv.isDarwin [
      "-Wno-unused-command-line-argument"
      "-Wmacro-redefined"
      "-Wincompatible-pointer-types-discards-qualifiers"
    ];
}
