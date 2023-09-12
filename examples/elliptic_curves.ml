open Pari_bindings

let int i = stoi (Signed.Long.of_int i)
let randomprime n = randomprime (int n)

let mkintmodu n modulus =
  mkintmodu (Unsigned.ULong.of_int n) (Unsigned.ULong.of_int modulus)

(* This line is mandatory to initialize the PARI stack, heap, table of primes... *)
let () = pari_init 50_000_000 (Unsigned.ULong.of_int 500_000)
let i = mkintmodu 2 101
let () = Printf.printf "\n%s\n" (gentostr i)
let ellinit2 ~a4 ~a6 ~p = ellinit (mkvec2 a4 a6) p (Signed.Long.of_int 38)

let ellinit ~a1 ~a2 ~a3 ~a4 ~a6 ~p =
  ellinit (mkvec5 a1 a2 a3 a4 a6) p (Signed.Long.of_int 38)

let p = randomprime (1 lsl 10)
let a = mkintmod (int 2) p

let e1 =
  ellinit2
    ~a4:(gpowgs a (Signed.Long.of_int 4))
    ~a6:(gpowgs a (Signed.Long.of_int 6))
    ~p

let () = Printf.eprintf "%s\n%s\n" (gentostr e1) (gentostr (ellcard e1 p))
let n = int 625 (* 5^4 *)

let g =
  ff_primroot
    (ffgen n (Signed.Long.of_int 0))
    Ctypes.(coerce (ptr void) (ptr gen) null)

let e2 =
  ellinit ~a1:(ff_zero g) ~a2:(ff_zero g) ~a3:(ff_zero g)
    ~a4:(ff_pow g (int 4))
    ~a6:(ff_pow g (int 6))
    ~p:(ff_p g)

let () = Printf.printf "%s\n%s" (gentostr e2) (gentostr (ellcard e2 (ff_p g)))