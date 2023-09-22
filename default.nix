{ pkgs ? import <nixpkgs> { } }:
let ocamlPackages = pkgs.ocaml-ng.ocamlPackages;
in pkgs.fastStdenv.mkDerivation {
  name = "ocaml_pari";
  nativeBuildInputs =
    (with ocamlPackages; [ ocaml findlib dune_3 merlin utop ppx_cstubs mdx odoc ])
    ++ (with pkgs; [
      bison
      gnumake
      pkg-config
      pkgconfig
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
#  LD_LIBRARY_PATH = "${pkgs.glibc}/lib:${pkgs.glibc.static}/lib";
  NIX_LDFLAGS = "-lc -lm";
}
