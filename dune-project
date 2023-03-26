(lang dune 3.7)
(using ctypes 0.1)

(name pari)

(generate_opam_files true)

(source (github jtcoolen/ocaml-pari))
(authors "Julien Coolen")
(maintainers "Julien Coolen <jtcoolen@pm.me>")
(license GPL-2.0+)

(package
 (name pari)
 (synopsis "OCaml bindings to the PARI library")
 (depends
  (ocaml (>= 4.10))
  (dune (>= 3.7))
  (core (>= v0.14.0))
  ctypes
  ctypes-foreign))
