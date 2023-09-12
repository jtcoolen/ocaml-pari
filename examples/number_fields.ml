open Pari

(* This line is mandatory to initialize the PARI stack;heap;table of primes... *)
let () = pari_init 500_000_000 (Unsigned.ULong.of_int 500_000)
let p = Polynomial.create [ (Integer.of_int 1, 1) ]
let () = Printf.eprintf "%s\n" (gentostr (Obj.magic p))
let p' = Polynomial.cyclotomic (Signed.Long.of_int 8)
let () = Printf.eprintf "%s\n" (gentostr (Obj.magic p'))

let q =
  Polynomial.create
    [
      (Integer.of_int 1, 3);
      (Integer.of_int (-111), 2);
      (Integer.of_int 6064, 1);
      (Integer.of_int (-189804), 0);
    ]

let () = Printf.eprintf "%s\n" (gentostr (Obj.magic q))
let () = Printf.eprintf "%b\n" (Polynomial.is_irreducible q)
let qmin : Integer.t Polynomial.t = Polynomial.minimal q
let () = Printf.eprintf "%s\n" (gentostr (Obj.magic qmin))

let () =
  Printf.eprintf "%b\n" NumberField.(are_isomorphic (create q) (create qmin))