{ system ? builtins.currentSystem, pkgs ? import <nixpkgs> { inherit system; } }:
let pkgs = import (fetchTarball {
  url = "https://github.com/NixOS/nixpkgs/archive/1f0e8ac1f9a783c4cfa0515483094eeff4315fe2.tar.gz";
  sha256 = "1mdnn0fj81pgvhzmzxh0g54g6yqxfqd2fim4h4c7cf7yskcp8g48";
}) {inherit system; }; in
let ocamlPackages = pkgs.ocaml-ng.ocamlPackages_5_0;
in pkgs.fastStdenv.mkDerivation {
  name = "ocaml_pari";
  nativeBuildInputs =
    (with ocamlPackages; [ ocaml findlib dune_3 merlin utop ppx_cstubs mdx odoc ])
    ++ (with pkgs; [
      bison
      gnumake
      pkg-config
      gcc
      gmp
      gmpxx
      libcxx
      glibc
      glibc.static
      nixfmt
      ocamlformat
    ]);
  buildInputs = (with ocamlPackages; [ core ppx_expect ctypes odoc ]);
  LD_LIBRARY_PATH = "${pkgs.glibc}/lib:${pkgs.glibc.static}/lib";
  NIX_LDFLAGS = "-lc -lm";
}
