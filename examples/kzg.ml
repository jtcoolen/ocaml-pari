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

let images_from_abscissa (ell : 'a Elliptic_curve.structure) (x : 'a) =
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
       [|
         inj_ring @@ one x;
         inj_ring ((a1 * x) + a3);
         inj_ring
           Infix.(
             ~-((x ^ Integer.of_int 3)
               + (a2 * (x ^ Integer.of_int 2))
               + (a4 * x) + a6));
       |])

(* l must be prime different from the characteristic of the field over which the curve is defined.

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
  srs_g1 : (Finite_field.t Elliptic_curve.t, [ `ROW ]) Vector.t;
  srs_g2 : (Finite_field.t Elliptic_curve.t, [ `ROW ]) Vector.t;
  finite_field_generator : Finite_field.t;
  curve : Finite_field.t Elliptic_curve.structure;
  curve_subgroup_order : Integer.t;
  g1 : Finite_field.t Elliptic_curve.t;
  g2 : Finite_field.t Elliptic_curve.t;
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
             Finite_field.(inj_ring @@ one fp);
             Finite_field.(inj_ring @@ zero fp);
             Finite_field.(inj_ring @@ one fp);
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
    Elliptic_curve.create
      ~a6:(Finite_field.finite_field_element [| Integer.of_int 1 |] quad_ext_p)
      ~dom:quad_ext_p ()
    |> Option.get
end

module KZG :
  Polynomial_commitment
    with type common_input = kzg_common_input
     and type polynomial = (finite_field, ring) t Polynomial.t
     and type scalar = Finite_field.t
     and type evaluation = Finite_field.t
     and type commitment = Finite_field.t Elliptic_curve.t = struct
  type common_input = kzg_common_input
  type polynomial = (finite_field, ring) t Polynomial.t
  type scalar = Finite_field.t
  type evaluation = Finite_field.t
  type commitment = Finite_field.t Elliptic_curve.t
  type proof = Finite_field.t Elliptic_curve.t

  let commit c p =
    assert (Polynomial.degree p < Vector.length c.srs_g1);
    Polynomial.fold_left2_vec
      ~f:(fun x p cm ->
        let n = Polynomial.(Finite_field.(residue_class (inj_field x)).%[0]) in
        Elliptic_curve.(add c.curve cm (mul c.curve ~n ~p)))
      ~acc:(Elliptic_curve.zero c.curve)
      p c.srs_g1

  let prove c p x =
    let y = Polynomial.(create [| eval p (Finite_field.inj_ring x) |]) in
    let n = Polynomial.(sub p y) in
    let d =
      Polynomial.create
        [|
          Finite_field.(inj_ring (one c.finite_field_generator));
          Finite_field.(inj_ring (Infix.( ~- ) x));
        |]
    in
    let q = Polynomial.div n d in
    commit c q

  let verify c cm x y pi =
    let d =
      Elliptic_curve.(
        sub c.curve
          Vector.(c.srs_g2.%[2])
          (mul c.curve
             ~n:Polynomial.((Finite_field.residue_class x).%[0])
             ~p:c.g2))
    in
    let n =
      Elliptic_curve.(
        sub c.curve cm
          (mul c.curve
             ~n:Polynomial.((Finite_field.residue_class y).%[0])
             ~p:c.g1))
    in
    let lhs =
      Elliptic_curve.weil_pairing_ff c.curve ~l:c.curve_subgroup_order ~p:pi
        ~q:d
    in
    let rhs =
      Elliptic_curve.weil_pairing_ff c.curve ~l:c.curve_subgroup_order ~p:n
        ~q:c.g2
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
  (*let sg = l_torsion_subgroup ToyCurve.curve ~l:(Signed.Long.of_int 641) in
    Printf.eprintf "\nsubgroup\n";
    Array.iter (fun e -> Printf.eprintf "\n%s\n" @@ gentostr e) sg;
    Printf.eprintf "\nend subgroup\n";*)
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
    (Elliptic_curve.weil_pairing_ff c.curve ~l:ToyCurve.r ~p:c.g1 ~q:r)
    (Elliptic_curve.weil_pairing_ff c.curve ~l:ToyCurve.r ~p:c.g2 ~q:r)

let rhs =
  Elliptic_curve.weil_pairing_ff c.curve ~l:ToyCurve.r
    ~p:(Elliptic_curve.add c.curve c.g1 c.g2)
    ~q:r

let () = assert (Finite_field.(equal lhs rhs))
