open Pari

(* This line is mandatory to initialize the PARI stack, heap, table of primes... *)
let () = pari_init 2_000_000_000 (Unsigned.ULong.of_int 500_000)
let p = Integer.(random_prime ~bits_amount:100)
let a = Finite_field.generator ~order:p
let () = Printf.printf "\n%s\n" (gentostr a)

let e1 =
  Elliptic_curve.create
    ~a4:(Finite_field.pow a (Integer.of_int 4))
    ~a6:(Finite_field.pow a (Integer.of_int 6))
    ()
  |> Option.get

let () = Printf.eprintf "%s\n" (Integer.to_string (Elliptic_curve.order e1))
let n = Integer.(pow (of_int 5) (of_int 4))
let g = Finite_field.generator ~order:n

let e2 =
  Elliptic_curve.create
    ~a4:(Finite_field.pow g (Integer.of_int 4))
    ~a6:(Finite_field.pow g (Integer.of_int 6))
    ()
  |> Option.get

let () = Printf.printf "%s\n" (Integer.to_string (Elliptic_curve.order e2))

