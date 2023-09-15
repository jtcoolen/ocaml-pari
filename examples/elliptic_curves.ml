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
    ~a6:(ff_pow g (Integer.of_int 6))
    ()
  |> Option.get

let () = Printf.printf "%s\n" (Integer.to_string (Elliptic_curve.order e2))

module type Polynomial_commitment = sig
  type common_input
  type polynomial
  type scalar
  type evaluation
  type commitment
  type proof

  val commit : common_input -> polynomial -> commitment
  val prove : common_input -> polynomial -> scalar -> proof

  val verify :
    common_input -> commitment -> scalar -> evaluation -> proof -> bool
end

let images_from_abscissa ell x =
  let open Elliptic_curve in
  let a1 = get_a1 ell in
  let a2 = get_a2 ell in
  let a3 = get_a3 ell in
  let a4 = get_a4 ell in
  let a6 = get_a6 ell in
  let open Finite_field in
  let open Finite_field.Infix in
  (*Y^2 + a_1 XY + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6*)
  Polynomial.roots_ff
    (Polynomial.create
       [
         (ff_1 x, 2);
         (inj_ring ((a1 * x) + a3), 1);
         ( inj_ring
             Infix.(
               ~-((x ^ Integer.of_int 3)
                 + (a2 * (x ^ Integer.of_int 2))
                 + (a4 * x) + a6)),
           0 );
       ])

(* l must be coprime with the characteristic of the field over which the curve is defined.
   Note that if l is prime then the Weil pairing is trivial (it always evaluates to 1).

   Inefficient way to find the elements of an l-torsion subgroup: find the abscissa of such
   a elements by looking at the roots of the l-division polynomial. *)
let _l_torsion_subgroup ell ~l =
  let l_div = Elliptic_curve.l_division_polynomial ell ~l in
  let l_div_roots = Polynomial.roots_ff l_div in
  let len = Vector.length l_div_roots in
  let c = ref 0 in
  let sg = Array.make (2 * len) (Elliptic_curve.zero ell) in
  for i = 1 to len do
    let x = Vector.(l_div_roots.%[i]) in
    let ys = images_from_abscissa ell x in
    for j = 1 to Vector.length ys do
      let y = Vector.(ys.%[j]) in
      let p = Elliptic_curve.of_coordinates ~x ~y in
      sg.(!c) <- p;
      c := !c + 1
    done
  done;
  Array.init !c (fun i -> sg.(i))

type kzg_common_input = {
  srs : (Finite_field.prime_field Elliptic_curve.t, [ `ROW ]) Vector.t;
  finite_field_generator : Finite_field.prime_field;
  curve : Finite_field.prime_field Elliptic_curve.structure;
  curve_subgroup_order : Integer.t;
  g1 : Finite_field.prime_field Elliptic_curve.t;
  g2 : Finite_field.prime_field Elliptic_curve.t;
}

module BLS = struct
  let p =
    Option.get
    @@ Integer.of_string
         "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab"

  let r =
    Option.get
    @@ Integer.of_string
         "0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"

  let g1_x =
    Option.get
    @@ Integer.of_string
         "0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb"

  let g1_y =
    Option.get
    @@ Integer.of_string
         "0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1"

  let fp = Finite_field.generator ~order:p

  let quad_ext_p =
    Finite_field.extend fp
      (`Quotient
        (Polynomial.create
           [
             (Finite_field.(inj_ring @@ ff_1 fp), 2);
             (Finite_field.(inj_ring @@ ff_1 fp), 0);
           ]))

  let () = Printf.eprintf "%s" (gentostr quad_ext_p)

  let g2_x =
    (*Finite_field.prime_field_element
      (Option.get
      @@ Integer.of_string
           "0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8"
      )
      ~p*)
    Finite_field.finite_field_element
      [|
        Option.get
        @@ Integer.of_string
             "0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8";
        Option.get
        @@ Integer.of_string
             "0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e";
      |]
      quad_ext_p

  let g2_y =
    Finite_field.finite_field_element
      [|
        Option.get
        @@ Integer.of_string
             "0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801";
        Option.get
        @@ Integer.of_string
             "0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be";
      |]
      quad_ext_p
  (*Finite_field.prime_field_element
    (Option.get
    @@ Integer.of_string
         "0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8"
    )
    ~p*)

  let g1 =
    Elliptic_curve.of_coordinates
      ~x:(Finite_field.prime_field_element g1_x ~p)
      ~y:(Finite_field.prime_field_element g1_y ~p)

  let g2 = Elliptic_curve.of_coordinates ~x:g2_x ~y:g2_y
end

module KZG :
  Polynomial_commitment
    with type common_input = kzg_common_input
     and type polynomial = Finite_field.prime Finite_field.ring Polynomial.t
     and type scalar = Finite_field.prime_field
     and type evaluation = Finite_field.prime_field
     and type commitment = Finite_field.prime_field Elliptic_curve.t = struct
  type common_input = kzg_common_input
  type polynomial = Finite_field.prime Finite_field.ring Polynomial.t
  type scalar = Finite_field.prime_field
  type evaluation = Finite_field.prime_field
  type commitment = Finite_field.prime_field Elliptic_curve.t
  type proof = Finite_field.prime_field Elliptic_curve.t

  let commit c p =
    assert (Polynomial.degree p < Vector.length c.srs);
    let cm = ref (Elliptic_curve.zero c.curve) in
    Printf.eprintf "\n\n\ncommit p(x)=%s, g1=%s\n" (gentostr p) (gentostr c.g1);
    for i = 0 to Polynomial.degree p do
      let n =
        Finite_field.(residue_class_prime (inj_field Polynomial.(p.%[i])))
      in
      let p = Vector.(c.srs.%[i + 1]) in
      Printf.eprintf "\nn=%s ;; p=%s\n" (gentostr n) (gentostr p);
      cm := Elliptic_curve.(add c.curve !cm (mul c.curve ~n ~p))
    done;
    Printf.eprintf "\nend commit\n\n\n";
    !cm

  let prove c (p : _ Finite_field.ring Polynomial.t) (x : _ Finite_field.t) :
      _ Finite_field.t Elliptic_curve.t =
    let y = Polynomial.(create [ (eval p (Finite_field.inj_ring x), 0) ]) in
    let n = Polynomial.(sub p y) in
    let d =
      Polynomial.create
        [
          (Finite_field.(inj_ring (ff_1 c.finite_field_generator)), 1);
          (Finite_field.(inj_ring (Infix.( ~- ) x)), 0);
        ]
    in
    let q = Polynomial.div n d in
    commit c q

  let verify c cm x y pi =
    let d =
      Elliptic_curve.(
        sub c.curve
          Vector.(c.srs.%[2])
          (mul c.curve ~n:(Finite_field.residue_class_prime x) ~p:c.g2))
    in
    let n =
      Elliptic_curve.(
        sub c.curve cm
          (mul c.curve ~n:(Finite_field.residue_class_prime y) ~p:c.g1))
    in
    let lhs =
      Elliptic_curve.weil_pairing_ff c.curve ~l:c.curve_subgroup_order ~p:pi
        ~q:d
    in
    let rhs =
      Elliptic_curve.weil_pairing_ff c.curve ~l:c.curve_subgroup_order ~p:n
        ~q:c.g2
    in
    Printf.eprintf
      "\npi=%s ;; d=%s ;; n=%s ;; g=%s ;; lhs=%s ;; rhs=%s ;; p=%s\n"
      (gentostr
      @@ ff_ellorder (Obj.magic c.curve) pi (Elliptic_curve.order c.curve))
      (gentostr d) (gentostr n) (gentostr c.g1) (gentostr lhs) (gentostr rhs)
      (gentostr
      @@ ff_order c.finite_field_generator Ctypes.(coerce (ptr void) t null));
    Finite_field.equal lhs rhs
end

let c =
  let g = Finite_field.generator ~order:BLS.p in
  let secret = Finite_field.random g in
  let curve =
    Elliptic_curve.create
      ~a6:(Finite_field.finite_field_element [|(Integer.of_int 3)|] BLS.quad_ext_p)
      ()
    |> Option.get
  in
  {
    srs =
      Vector.of_array
        (Array.init 100 (fun i ->
             Elliptic_curve.mul curve
               ~n:
                 Finite_field.(
                   residue_class_prime (pow secret (Integer.of_int i)))
               ~p:BLS.g1));
    finite_field_generator = g;
    curve = Obj.magic curve;
    curve_subgroup_order = BLS.r;
    g1 = BLS.g1;
    g2 = BLS.g2;
  }

let p =
  Polynomial.create
    [
      (Finite_field.(inj_ring (random c.finite_field_generator)), 2);
      (Finite_field.(inj_ring (random c.finite_field_generator)), 0);
    ]

let () = Printf.eprintf "\nP=%s\n" (gentostr p)
let () = Printf.eprintf "\nP[1]=%s\n" (gentostr Polynomial.(p.%[1]))
let () = Printf.eprintf "\nP[0]=%s\n" (gentostr Polynomial.(p.%[0]))
let () = Printf.eprintf "\nP[deg]=%s\n" (gentostr Polynomial.(p.%[degree p]))

let () =
  Printf.eprintf "\nP[deg-1]=%s\n" (gentostr Polynomial.(p.%[degree p - 1]))

let () =
  Printf.eprintf "\nP[deg+1]=%s\n" (gentostr Polynomial.(p.%[degree p + 1]))

let cm = KZG.commit c p
let x = Finite_field.random c.finite_field_generator
let y = Finite_field.inj_field (Polynomial.eval p (Finite_field.inj_ring x))
let pi = KZG.prove c p x

let () =
  Printf.eprintf "\nverify=%b\n"
    (KZG.verify c cm x y (*Finite_field.random c.finite_field_generator*) pi)

let p' =
  Polynomial.create
    [
      (Finite_field.(inj_ring (random c.finite_field_generator)), 1);
      (Finite_field.(inj_ring (random c.finite_field_generator)), 0);
    ]

let cm' = KZG.commit c p'
let lhs = Elliptic_curve.(add c.curve cm cm')
let rhs = KZG.commit c (Polynomial.add p p')
let () = Printf.eprintf "\nlhs=%s ;; rhs=%s\n" (gentostr lhs) (gentostr rhs)
let () = assert (Elliptic_curve.(equal lhs rhs))
let x = Finite_field.random c.finite_field_generator
let y = Finite_field.inj_field (Polynomial.eval p (Finite_field.inj_ring x))
let pi = KZG.prove c p' x
let () = Printf.eprintf "\nverify=%b\n" (KZG.verify c cm' x y pi)
