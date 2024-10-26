open Pari

(* This line is mandatory to initialize the PARI stack;heap;table of primes... *)
let () = pari_init 500_000_000 (Unsigned.ULong.of_int 500_000)
let p = Polynomial.create [| Integer.of_int 1; Integer.of_int 0 |]
let () = Printf.eprintf "%s\n" (Polynomial.to_string p)
let p' = Polynomial.cyclotomic (Signed.Long.of_int 8)
let () = Printf.eprintf "%s\n" (Polynomial.to_string p')

let q =
  Polynomial.create
    [|
      Integer.of_int 1;
      Integer.of_int (-111);
      Integer.of_int 6064;
      Integer.of_int (-189804);
    |]

let zero = Polynomial.create [| Integer.of_int 0 |]
let qq = Polynomial.create [| q; q; zero; zero |]
let () = Printf.eprintf "%s\n" (Polynomial.to_string qq)
let () = Printf.eprintf "%b\n" (Polynomial.is_irreducible q)
let qmin : Integer.t Polynomial.t = Polynomial.minimal q
let () = Printf.eprintf "%s\n" (Polynomial.to_string qmin)
let inj_rat = Polynomial.inj_base_ring ~inj:(fun x -> x |> Integer.inj_rat)

let () =
  Printf.eprintf "%b\n"
    Number_field.(are_isomorphic (create (inj_rat q)) (create (inj_rat qmin)))

(* Gaussian integers: the ring Z[i] (here we work in the field Q(i)) *)
let gaussian_integers =
  (* Q(i) = Q[X]/(X^2+1) *)
  Number_field.create
    (Polynomial.create
       [| Integer.of_int 1; Integer.of_int 0; Integer.of_int 1 |]
    |> inj_rat)

(* Euclidean division of 6 + 8i by 1 + 5i. *)
let a =
  Number_field.elt
    [| Integer.(inj_rat (of_int 6)); Integer.(inj_rat (of_int 8)) |]

let b =
  Number_field.elt
    [| Integer.(inj_rat (of_int 1)); Integer.(inj_rat (of_int 5)) |]

let q, r = Number_field.divrem gaussian_integers a b

let () =
  Printf.eprintf "%b\n"
    Number_field.(equal a (add gaussian_integers (mul gaussian_integers b q) r))

let _nf2 =
  Number_field.create (Option.get @@ Polynomial.of_string "x^4-2" |> inj_rat)

let a = Vector.((Number_field.z_basis _nf2).%[2])
let a2 = Number_field.mul _nf2 a a
let _ = Printf.eprintf "%s\n" (gentostr a2)
let a4 = Number_field.mul _nf2 a2 a2
let _ = Printf.eprintf "a4=%s\n" (gentostr a4)

let mp =
  nf_to_scalar_or_alg (Obj.magic _nf2)
    (Number_field.add _nf2 a2
       (Number_field.elt
          [|
            Integer.of_int (-3) |> Integer.inj_rat;
            Integer.of_int 0 |> Integer.inj_rat;
            Integer.of_int 0 |> Integer.inj_rat;
            Integer.of_int 0 |> Integer.inj_rat;
          |]))

let _ = Printf.eprintf "mp=%s\n" (gentostr mp)
let _ = Printf.eprintf "mp=%s\n" (gentostr (Number_field.splitting (`F _nf2)))
let _ = Printf.eprintf "mp=%s\n" (gentostr (polred (Obj.magic mp)))
