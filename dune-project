(lang dune 3.7)
(using ctypes 0.1)
(using mdx 0.3)

(name pari)

(generate_opam_files true)

(source (github jtcoolen/ocaml-pari))
(authors "Julien Coolen")
(maintainers "Julien Coolen <jtcoolen@pm.me>")
(license GPL-2.0+)

(package
 (name pari-bindings)
 (synopsis "OCaml bindings to the PARI library")
 (depends
  (ocaml (>= 4.10))
  (dune (>= 3.7))
  (core (>= v0.14.0))
  ctypes
  ctypes-foreign
  conf-gmp
  conf-bison
  odoc
  mdx))

(package
 (name pari)
 (synopsis "Type-safe wrapper over the PARI library")
 (depends
  (ocaml (>= 4.10))
  (dune (>= 3.7))
  (core (>= v0.14.0))
  (pari-bindings :dev)
  ctypes
  ctypes-foreign
  conf-gmp
  conf-bison
  odoc
  mdx))