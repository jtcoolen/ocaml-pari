open Pari

(* This line is mandatory to initialize the PARI stack;heap;table of primes... *)
let () = pari_init 500_000_000 (Unsigned.ULong.of_int 500_000)

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

type kzg_common_input = {
  srs_g1 : (Finite_field.t Elliptic_curve.elt, [ `ROW ]) Vector.t;
  srs_g2 : (Finite_field.t Elliptic_curve.elt, [ `ROW ]) Vector.t;
  finite_field_generator : Finite_field.t;
  curve : Finite_field.t Elliptic_curve.t;
  curve_subgroup_order : Integer.t;
  g1 : Finite_field.t Elliptic_curve.elt;
  g2 : Finite_field.t Elliptic_curve.elt;
}

module ToyCurve = struct
  let p = Integer.of_int 7691
  let r = Integer.of_int 641
  let fp = Finite_field.generator ~order:p

  let quad_ext_p =
    Finite_field.extend fp
      (`Quotient
        (Polynomial.create
           [|
             Finite_field.(inj_ring (one fp));
             Finite_field.(inj_ring (zero fp));
             Finite_field.(inj_ring (one fp));
           |]))

  let g1_x =
    Finite_field.finite_field_element [| Integer.of_int 2693 |] quad_ext_p

  let g1_y =
    Finite_field.finite_field_element [| Integer.of_int 4312 |] quad_ext_p

  let g2_x =
    Finite_field.finite_field_element
      [| Integer.of_int 633; Integer.of_int 6145 |]
      quad_ext_p

  let g2_y =
    Finite_field.finite_field_element
      [| Integer.of_int 7372; Integer.of_int 109 |]
      quad_ext_p

  let g1 = Elliptic_curve.of_coordinates ~x:g1_x ~y:g1_y
  let g2 = Elliptic_curve.of_coordinates ~x:g2_x ~y:g2_y

  let curve =
    Option.get
      (Elliptic_curve.create
         ~a6:
           (Finite_field.finite_field_element [| Integer.of_int 1 |] quad_ext_p)
         ~dom:quad_ext_p ())
end

module KZG :
  Polynomial_commitment
    with type common_input = kzg_common_input
     and type polynomial = (finite_field, ring) typ Polynomial.t
     and type scalar = Finite_field.t
     and type evaluation = Finite_field.t
     and type commitment = Finite_field.t Elliptic_curve.elt = struct
  type common_input = kzg_common_input
  type polynomial = (finite_field, ring) typ Polynomial.t
  type scalar = Finite_field.t
  type evaluation = Finite_field.t
  type commitment = Finite_field.t Elliptic_curve.elt
  type proof = Finite_field.t Elliptic_curve.elt

  let commit c p =
    assert (Polynomial.degree p < Vector.length c.srs_g1);
    let f x p cm =
      let n = Polynomial.(Finite_field.(residue_class (inj_field x)).%[0]) in
      Elliptic_curve.(add c.curve cm (mul c.curve ~n ~p))
    in
    let acc = Elliptic_curve.zero c.curve in
    Polynomial.fold_left2_vec ~f ~acc p c.srs_g1

  let prove c p x =
    let y = Polynomial.(create [| eval p (Finite_field.inj_ring x) |]) in
    let n = Polynomial.(sub p y) in
    let d =
      Polynomial.create
        [|
          Finite_field.(inj_ring (one c.finite_field_generator));
          Finite_field.(inj_ring (neg x));
        |]
    in
    let q = Polynomial.div n d in
    commit c q

  let verify c cm x y pi =
    let denominator =
      Elliptic_curve.(
        sub c.curve
          Vector.(c.srs_g2.%[2])
          (mul c.curve
             ~n:Polynomial.((Finite_field.residue_class x).%[0])
             ~p:c.g2))
    in
    let numerator =
      Elliptic_curve.(
        sub c.curve cm
          (mul c.curve
             ~n:Polynomial.((Finite_field.residue_class y).%[0])
             ~p:c.g1))
    in
    let lhs =
      Elliptic_curve.weil_pairing c.curve ~l:c.curve_subgroup_order pi
        denominator
    in
    let rhs =
      Elliptic_curve.weil_pairing c.curve ~l:c.curve_subgroup_order numerator
        c.g2
    in
    Finite_field.equal lhs rhs
end

let c =
  let g = Finite_field.generator ~order:ToyCurve.r in
  let secret = Finite_field.random g in
  let coeff p i =
    Elliptic_curve.mul ToyCurve.curve
      ~n:
        Polynomial.(
          Finite_field.(residue_class (pow secret (Integer.of_int i))).%[0])
      ~p
  in
  let srs_g1 = Vector.init 100 ~f:(coeff ToyCurve.g1) in
  let srs_g2 = Vector.init 100 ~f:(coeff ToyCurve.g2) in
  {
    srs_g1;
    srs_g2;
    finite_field_generator = g;
    curve = ToyCurve.curve;
    curve_subgroup_order = ToyCurve.r;
    g1 = ToyCurve.g1;
    g2 = ToyCurve.g2;
  }

let p =
  Polynomial.create
    [|
      Finite_field.(inj_ring (random c.finite_field_generator));
      Finite_field.(inj_ring (zero c.finite_field_generator));
      Finite_field.(inj_ring (random c.finite_field_generator));
    |]

let cm = KZG.commit c p
let x = Finite_field.random c.finite_field_generator
let y = Finite_field.inj_field (Polynomial.eval p (Finite_field.inj_ring x))
let pi = KZG.prove c p x

let () =
  assert (
    not (KZG.verify c cm x (Finite_field.add y c.finite_field_generator) pi))

let p' =
  Polynomial.create
    [|
      Finite_field.(inj_ring (random c.finite_field_generator));
      Finite_field.(inj_ring (random c.finite_field_generator));
    |]

let cm' = KZG.commit c p'
let lhs = Elliptic_curve.(add c.curve cm cm')
let rhs = KZG.commit c (Polynomial.add p p')

(* KZG commitments are homomorphic *)
let () = assert (Elliptic_curve.(equal lhs rhs))
let x' = Finite_field.random c.finite_field_generator
let y' = Finite_field.inj_field (Polynomial.eval p' (Finite_field.inj_ring x'))
let pi' = KZG.prove c p' x'
let () = assert (KZG.verify c cm' x' y' pi')

(* test bilinearity of the Weil pairing *)
let r =
  Elliptic_curve.mul c.curve ~n:(Integer.of_int (Random.int 100000)) ~p:c.g2

let lhs =
  Finite_field.mul
    (Elliptic_curve.weil_pairing c.curve ~l:ToyCurve.r c.g1 r)
    (Elliptic_curve.weil_pairing c.curve ~l:ToyCurve.r c.g2 r)

let rhs =
  Elliptic_curve.weil_pairing c.curve ~l:ToyCurve.r
    (Elliptic_curve.add c.curve c.g1 c.g2)
    r

let () = assert (Finite_field.(equal lhs rhs))
