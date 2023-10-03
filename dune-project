(lang dune 3.7)
(using ctypes 0.3)
(using mdx 0.3)

(name pari)

(generate_opam_files true)

(source (github jtcoolen/ocaml-pari))
(authors "Julien Coolen")
(maintainers "Julien Coolen <jtcoolen@pm.me>")
(license GPL-2.0-or-later)

(package
 (name pari-bindings)
 (synopsis "OCaml bindings to the PARI library")
 (depends
  (ocaml (>= 4.10))
  (dune (>= 3.7))
  (core (>= v0.14.0))
  (ctypes (>= 0.13.0))
  (ctypes-foreign (>= 0.13.0))
  conf-gmp
  conf-bison
  conf-perl
  odoc
  (mdx (>= 2.0.0))))

(package
 (name pari)
 (synopsis "Type-safe wrapper over the PARI library")
 (depends
  (ocaml (>= 4.10))
  (dune (>= 3.7))
  (core (>= v0.14.0))
  (pari-bindings (= :version))
  (ctypes (>= 0.13.0))
  (ctypes-foreign (>= 0.13.0))
  conf-gmp
  conf-bison
  conf-perl
  odoc
  (mdx (>= 2.0.0))))
