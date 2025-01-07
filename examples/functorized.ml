open Pari
module P = F.Polynomial (Integer)

let a = P.create Integer.([| 1; 2; 3 |] |> Array.map of_int)
let () = Printf.eprintf "%s\n" (gentostr (P.to_gen a))
let b = P.create Integer.([| 1; 2; 3 |] |> Array.map of_int)
let a_plus_b = P.add a b
let () = Printf.eprintf "%s" (gentostr (P.to_gen a_plus_b))

module P2 = F.Polynomial (Integer)

let c = P2.create Integer.([| 1; 2; 3 |] |> Array.map of_int)
let a_plus_c = P.add a c
let () = Printf.eprintf "%s\n" (gentostr (P.to_gen a_plus_c))

(********)

module P3 = F.Polynomial (Integer_mod)

let _d =
  P3.create
    Integer_mod.(
      [| 1; 2; 3 |]
      |> Array.map @@ fun x ->
         create (Integer.of_int x) ~modulo:(Integer.of_int 11))

(* echoue *)
(* let a_plus_d = P3.add a d *)
(* let () = Printf.eprintf "%s\n" (gentostr (P.to_gen a_plus_d)) *)

module F11 = F.IntegerMod (struct
  let modulus = 11
end)

module F7 = F.IntegerMod (struct
  let modulus = 7
end)

module P4 = F.Polynomial (F11)
module P5 = F.Polynomial (F7)
module P6 = F.Polynomial (F.Polynomial (F7))

let d = P4.create F11.([| 1; 2; 3 |] |> Array.map create)
let d' = P5.create F7.([| 1; 2; 3 |] |> Array.map create)
let pol = P6.create [| d'; d'; d' |]
let () = Printf.eprintf "pol = %s\n" (gentostr (P6.to_gen pol))

(* echoue *)
(* let _ = P5.add d d' *)

module P4Ring : F.Ring = P4

let _ = P4Ring.add
let _ = P4Ring.mul
let sqr (type t) (module R : F.Ring with type t = t) (x : t) = R.mul x x
let p4 = (module P4 : F.Ring with type t = P4.t)
let s1 = sqr p4 d
let s2 = sqr (module P4) d

let () =
  Printf.eprintf "%s == %s, equal ? %b\n"
    (gentostr (P4.to_gen s1))
    (gentostr (P4.to_gen s2))
    (gequal (P4.to_gen s1) (P4.to_gen s2))
