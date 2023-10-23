include Pari_bindings

type ('kind, 'structure) t = gen

let t = gen

type ring = gen
type field = gen

let register_gc v =
  Gc.finalise_last (fun () -> pari_free Ctypes.(coerce gen (ptr void) v)) v

module Complex = struct
  type complex = private Complex
  type t = gen

  let create ~re ~im =
    let c = mkcomplex re im in
    register_gc c;
    c

  let inv = ginv
  let add = gadd
  let to_string = gentostr
end

module Real = struct
  type real = private Real
  type t = gen

  external inj_complex : t -> Complex.t = "%identity"

  let of_signed = stor
  let shift n s = mpshift n (Signed.Long.of_int s)
  let sqrt = sqrtr
  let ceil = gceil
  let add = gadd
  let inv = ginv
end

module Rational = struct
  type rational = private Rational
  type t = gen
  type ring = gen

  external inj_ring : t -> ring = "%identity"
  external inj_real : t -> Real.t = "%identity"
  external inj_complex : t -> Complex.t = "%identity"

  let shift x s = mpshift x (Signed.Long.of_int s)
end

module Integer = struct
  type integer = private Integer
  type t = gen

  external inj_rat : t -> Rational.t = "%identity"
  external inj_real : t -> Rational.t = "%identity"
  external inj_complex : t -> Rational.t = "%identity"

  let equal x y = equalii x y = 1
  let of_int i = stoi (Signed.Long.of_int i)
  let of_signed = stoi
  let shift x s = mpshift x (Signed.Long.of_int s)
  let sqrt i = inj_complex (sqrti i)
  let zero = stoi Signed.Long.zero
  let mul = mulii
  let add = addii
  let sub = subii
  let neg = negi
  let pow = powii
  let modulo = modii
  let of_string s = (* todo error handling *) Some (gp_read_str s)
  let to_string = gentostr

  let random_prime ~bits_amount =
    randomprime (mpshift (of_int 1) (Signed.Long.of_int bits_amount))

  module Infix = struct
    let ( * ) = mul
    let ( + ) = add
    let ( - ) = sub
    let ( ~- ) = neg
    let ( mod ) = modulo
    let ( = ) = equal
  end
end

module Matrix = struct
  type nonrec 'a t = gen constraint 'a = gen

  let id n =
    let m = matid (Signed.Long.of_int n) in
    register_gc m;
    m

  let inv = ginv
  let mul = gmul
  let lll = lll

  (* Pointer dereference operator *)
  let ( !@ ) = Ctypes.( !@ )

  (* Multi-index operators *)
  let getcoeff m = function
    | [| i; j |] -> safegcoeff m (Signed.Long.of_int i) (Signed.Long.of_int j)
    | [| i |] -> safegel m (Signed.Long.of_int i)
    | _ -> assert false

  let ( .%[] ) m i = Ctypes.( !@ ) (safegel m (Signed.Long.of_int i))
  let ( .%[]<- ) m i v = Ctypes.(safegel m (Signed.Long.of_int i) <-@ v)
  let ( .%[;..] ) m i = !@(getcoeff m i)
  let ( .%[;..]<- ) m i v = Ctypes.(getcoeff m i <-@ v)
  let inj x ~inj:_ = x
end

module Set = struct
  type nonrec 'a t = gen constraint 'a = gen

  let length (s : 'a t) = glength s
  let search (s : 'a t) (e : 'a) = setsearch s e
end

module Vector = struct
  type ('a, 'b) p = ('a, 'b) t
  type ('a, 'b) t = gen constraint 'a = gen constraint 'b = [< `COL | `ROW ]

  let length x = glength x |> Signed.Long.to_int
  let ( .%[] ) m i = Ctypes.( !@ ) (safegel m (Signed.Long.of_int i))
  let ( .%[]<- ) m i v = Ctypes.(safegel m (Signed.Long.of_int i) <-@ v)

  let of_array a =
    let l = Array.length a in
    let v = zerovec_block (Signed.Long.of_int l) in
    for i = 1 to l do
      v.%[i] <- a.(i - 1)
    done;
    register_gc v;
    v

  let init l ~f =
    let v = zerovec_block (Signed.Long.of_int l) in
    for i = 1 to l do
      v.%[i] <- f (i - 1)
    done;
    register_gc v;
    v

  let equal x y = gequal x y == 1

  let slice m ~start ~stop =
    vecslice m (Signed.Long.of_int start) (Signed.Long.of_int stop)

  let mul x y = gmul x y
  let add = gadd
  let sub = gsub
  let neg = gneg
  let transpose_row = gtrans
  let transpose_column = gtrans
  let to_set = gtoset
  let singleton x = mkvec x
  let concat = gconcat
  let inj x ~inj:_ = x
  let to_string = gentostr

  module Infix = struct
    let ( = ) = equal
    let ( * ) = mul
    let ( + ) = add
    let ( - ) = sub
    let ( ~- ) = neg
  end
end

let () = pari_init 50_000_000 (Unsigned.ULong.of_int 500_000)

module Polynomial = struct
  type polynomial
  type 'a t = gen constraint 'a = gen

  let equal x y = gequal x y = 1
  let to_string = gentostr
  let mul = gmul
  let div = gdiv
  let equal = equal
  let add = gadd
  let sub = gsub
  let neg = gneg
  let eval p x = poleval p x
  let degree t = degree t |> Signed.Long.to_int

  let get_coeff t i =
    gcopy @@ Ctypes.(!@(coerce gen (ptr gen) t +@ Int.add i 2))

  let create (p : 'a array) : 'a t =
    let len = Array.length p in
    let size = Signed.Long.of_int (Int.add len 2) in
    let p' = cgetg size (Signed.Long.of_int 10 (* t_POL *)) in
    let p'' = Ctypes.(coerce gen (ptr gen) p') in
    (* TODO: support 32-bit arch *)
    let typ = gentostr (type0 p.(0)) in
    let v =
      if typ = "\"t_POL\"" then Signed.Long.(succ (gvar p.(0)))
      else Signed.Long.zero
    in
    Ctypes.(p' +@ 1 <-@ Signed.Long.(shift_left v Stdlib.(64 - 2 - 16)));
    for i = 0 to len - 1 do
      Ctypes.(p'' +@ (i + 2) <-@ p.(len - 1 - i))
    done;
    let p = normalizepol_lg p' size in
    register_gc p;
    p

  let var t = function Some var -> Signed.Long.(of_int var) | None -> gvar t
  let deriv ?indeterminate t = deriv t (var t indeterminate)

  let derivn ?indeterminate t n =
    derivn t (Signed.Long.of_int n) (var t indeterminate)

  let cyclotomic n = polcyclo n (Signed.Long.of_int 0)
  let is_irreducible p = polisirreducible p = Signed.Long.one
  let minimal p = polredbest p (Signed.Long.of_int 0)
  let ( .%[] ) m i = gcopy @@ truecoef m (Signed.Long.of_int i)
  let roots_ff p = polrootsmod p Ctypes.(coerce (ptr void) t null)

  let fold_left2 ~f ~acc p p' =
    assert (degree p = degree p');
    let acc = ref acc in
    for i = 0 to degree p do
      acc := f p.%[i] p'.%[i] !acc
    done;
    !acc

  let fold_left2_vec ~f ~acc p v =
    assert (degree p < Vector.length v);
    let acc = ref acc in
    for i = 0 to degree p do
      acc := f p.%[i] Vector.(v.%[i + 1]) !acc
    done;
    !acc
end

module Number_field = struct
  type number_field = private Number_field
  type nonrec t = (number_field, field) t
  type structure = gen

  let create p =
    let nf = nfinit p Signed.Long.(of_int 4) in
    register_gc nf;
    nf

  let are_isomorphic a b = gequal0 (nfisincl a b) <> 1

  let sign nf =
    let a = Ctypes.(allocate long Signed.Long.zero) in
    let b = Ctypes.(allocate long Signed.Long.zero) in
    nf_get_sign nf a b;
    Ctypes.(!@a, !@b)

  let discriminant nf = nf_get_disc nf
  let z_basis nf = nf_get_zk nf
  let elt a = Vector.(transpose_row (of_array a))
  let add nf a b = nfadd nf a b
  let mul nf a b = nfmul nf a b

  let divrem nf a b =
    let qr = nfdivrem nf a b in
    Vector.(qr.%[1], qr.%[2])

  let equal a b = gequal a b = 1
  let ideal_norm nf a = idealnorm nf a

  let splitting = function
    | `F nf -> nfsplitting nf Ctypes.(coerce (ptr void) gen null)
    | `P p -> nfsplitting p Ctypes.(coerce (ptr void) gen null)

  module Infix = struct
    let ( = ) = equal
  end
end

type group

type 'a group_structure = {
  mul : ('a, group) t -> ('a, group) t -> ('a, group) t;
  pow : ('a, group) t -> Integer.t -> ('a, group) t;
  rand : unit -> ('a, group) t;
  hash : ('a, group) t -> Unsigned.ULong.t;
  equal : ('a, group) t -> ('a, group) t -> bool;
  equal_identity : ('a, group) t -> bool;
  bb_group : bb_group Ctypes.structure option;
}

module Fp = struct
  type t = Integer.t

  let add a b ~modulo = fp_add a b modulo
  let pow x ~exponent ~modulo = fp_pow x exponent modulo
end

type finite_field = private Finite_field

module Finite_field = struct
  type ('a, 'b) p = ('a, 'b) t
  type t = (finite_field, field) p

  external inj_ring : t -> (finite_field, ring) p = "%identity"
  external inj_field : (finite_field, ring) p -> t = "%identity"

  let generator ~order =
    ff_primroot
      (ffgen order Signed.Long.zero)
      Ctypes.(coerce (ptr void) (ptr gen) null)

  let prime_field_element x ~p = ff_z_mul (ff_1 (generator ~order:p)) x

  let finite_field_element coeffs a =
    let len = Array.length coeffs in
    fst
      (Array.fold_left
         (fun (acc, i) e ->
           (gadd (gmul e (ff_pow a (Integer.of_int (len - i - 1)))) acc, i + 1))
         (ff_zero a, 0)
         coeffs)

  let generator_from_irreducible_polynomial p = ffgen p Signed.Long.zero
  let residue_class x = ff_to_fpxq_i x

  let create ~p ~degree =
    ffinit (Integer.of_int p) (Signed.Long.of_int degree) Signed.Long.zero

  let equal a b = ff_equal a b = 1
  let add = ff_add
  let mul = ff_mul
  let pow x n = ff_pow x n
  let random = genrand
  let zero g = ff_zero g
  let one e = ff_1 e

  let extend base_field_elt = function
    | `Degree degree ->
        Vector.(
          (ffextend base_field_elt
             (ffinit (ff_p base_field_elt)
                (Signed.Long.of_int degree)
                Signed.Long.zero)
             Signed.Long.zero).%[1])
    | `Quotient modulo ->
        Vector.((ffextend base_field_elt modulo Signed.Long.zero).%[1])

  let fpxq_star ~(p : pari_ulong) ~(quotient : Fp.t Polynomial.t) :
      finite_field group_structure =
    let open Ctypes in
    let q = powuu p (Unsigned.ULong.of_int (Polynomial.degree quotient)) in
    let ret = allocate (ptr void) (from_voidp void null) in
    let bb_group = !@(get_flxq_star ret quotient (itou q)) in
    let get_fct field typ =
      let fn = getf bb_group field in
      coerce (static_funptr typ) (Foreign.funptr typ) fn
    in
    {
      bb_group = Some bb_group;
      mul = (fun a b -> flxq_mul a b quotient p);
      pow = (fun x n -> flxq_pow x n quotient p);
      rand =
        (fun () -> (get_fct bb_group_rand (ptr void @-> returning gen)) null);
      hash = hash_gen;
      equal = (fun a b -> flx_equal a b = 1);
      equal_identity = (fun a -> flx_equal1 a = 1);
    }

  let to_string = gentostr

  module Infix = struct
    let ( ~- ) = ff_neg
    let ( + ) = ff_add
    let ( - ) = ff_sub
    let ( * ) = ff_mul
    let ( ^ ) = ff_pow
  end
end

module Elliptic_curve = struct
  type elliptic_curve
  type nonrec 'a t = (elliptic_curve, group) t constraint 'a = gen
  type 'a structure = gen constraint 'a = gen

  let create ?a1 ?a2 ?a3 ?a4 ?a6 ?(dom = Ctypes.(coerce (ptr void) gen null)) ()
      =
    let to_coeff = function Some c -> c | None -> stoi Signed.Long.zero in
    let ell =
      ellinit
        (mkvec5 (to_coeff a1) (to_coeff a2) (to_coeff a3) (to_coeff a4)
           (to_coeff a6))
        dom (Signed.Long.of_int 38)
    in
    if Vector.length ell = 0 then None
    else (
      register_gc ell;
      Some ell)

  let get_a1 = ell_get_a1
  let get_a2 = ell_get_a2
  let get_a3 = ell_get_a3
  let get_a4 = ell_get_a4
  let get_a6 = ell_get_a6
  let of_coordinates ~x ~y = mkvec2 x y
  let order ell = ellcard ell Ctypes.(coerce (ptr void) gen null)
  let discriminant ell = ell_get_disc ell
  let j_invariant ell = ell_get_j ell
  let random = ellrandom
  let l_division_polynomial ell ~l = elldivpol ell l Signed.Long.zero
  let to_string = gentostr
  let weil_pairing_ff ell ~l ~p ~q = ellweilpairing ell p q l
  let add ell p q = elladd ell p q
  let sub ell p q = ellsub ell p q
  let mul ell ~n ~p = ellmul ell p n
  let equal a b = gequal a b = 1
  let generators_ff ell = ff_ellgens ell
  let zero _ = ellinf ()
end
