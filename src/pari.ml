include Pari_bindings

type 'kind ty = gen

let t = gen

type group = Group
type ring = Ring
type field = Field
type unique_factorization_domain = Unique_factorization_domain
type complex = Complex
type real = Real
type rational = Rational
type integer = Integer
type 'a polynomial = Polynomial of 'a
type integer_mod = Integer_mod
type finite_field = Finite_field
type number_field = Number_field
type 'a elliptic_curve = Elliptic_curve of 'a

let register_gc v =
  Gc.finalise_last (fun () -> pari_free Ctypes.(coerce gen (ptr void) v)) v

let gentostr = gentostr_raw
let gentobytes x = gentostr x |> String.to_bytes

module Complex = struct
  type t = gen

  let create ~re ~im = mkcomplex re im
  let inv = ginv
  let add = gadd
  let to_string = gentostr
end

module Real = struct
  type t = gen

  let[@inline] inj_complex x = Fun.id x
  let of_signed = stor
  let shift n s = mpshift n (Signed.Long.of_int s)
  let sqrt = sqrtr
  let ceil = gceil
  let add = gadd
  let inv = ginv
end

module Rational = struct
  type t = gen

  let of_int i = stoi (Signed.Long.of_int i)
  let[@inline] inj_real x = Fun.id x
  let[@inline] inj_complex x = Fun.id x
  let shift x s = mpshift x (Signed.Long.of_int s)
end

module Integer = struct
  type t = gen

  let[@inline] inj_rat x = Fun.id x
  let[@inline] inj_real x = Fun.id x
  let[@inline] inj_complex x = Fun.id x
  let equal x y = equalii x y = 1
  let of_int i = stoi (Signed.Long.of_int i)
  let to_int i = Signed.Long.to_int (itos i)
  let of_signed = stoi

  let of_hex (`Hex h) =
    fromdigits (gtovec (gtovecsmall (strtogenstr h))) (of_int 16)

  let shift x s = mpshift x (Signed.Long.of_int s)
  let sqrt i = inj_complex (sqrti i)
  let zero () = of_int 0
  let mul = mulii
  let add = addii
  let sub = subii
  let neg = negi
  let pow = powii
  let modulo = modii
  let of_string s = (* todo error handling *) Some (gp_read_str s)
  let to_string = gentostr

  let gcdext a b =
    let u = Ctypes.allocate t a in
    let v = Ctypes.allocate t a in
    let gcd = bezout a b u v in
    (gcd, Ctypes.( !@ ) u, Ctypes.( !@ ) v)

  let gcd = gcdii
  let divexact = diviiexact

  let random_prime ~bits_amount =
    randomprime (mpshift (of_int 1) (Signed.Long.of_int bits_amount))

  let random = randomi

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
  let[@inline] inj x ~inj:_ = x

  let dimensions m =
    let sz = matsize m in
    (Signed.Long.to_int (itos sz.%[1]), Signed.Long.to_int (itos sz.%[2]))
end

module Set = struct
  type nonrec 'a t = gen constraint 'a = gen

  let length (s : 'a t) = glength s
  let search (s : 'a t) (e : 'a) = setsearch s e
end

module Vector = struct
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

  let array_map ~f a =
    let l = Array.length a in
    let v = zerovec_block (Signed.Long.of_int l) in
    for i = 1 to l do
      v.%[i] <- f a.(i - 1)
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
  let[@inline] inj x ~inj:_ = x
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

  let of_string s =
    let g = gp_read_str s in
    let typ = gentostr (type0 g) in
    if typ = "\"t_POL\"" then Some g else None

  let var t = function Some var -> Signed.Long.(of_int var) | None -> gvar t
  let deriv ?indeterminate t = deriv t (var t indeterminate)

  let derivn ?indeterminate t n =
    derivn t (Signed.Long.of_int n) (var t indeterminate)

  let cyclotomic n = polcyclo n (Signed.Long.of_int 0)
  let is_irreducible p = polisirreducible p = Signed.Long.one
  let minimal p = polredbest p (Signed.Long.of_int 0)
  let ( .%[] ) m i = gcopy @@ truecoef m (Signed.Long.of_int i)
  let roots_ff p = polrootsmod p Ctypes.(coerce (ptr void) t null)

  let fold_left ~f ~acc p =
    let acc = ref acc in
    for i = 0 to degree p do
      acc := f p.%[i] !acc
    done;
    !acc

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

  let inj_base_ring ~inj:_ p = p
end

module Integer_mod = struct
  type t = gen

  let[@inline] inj_group x = Fun.id x

  let create : Integer.t -> modulo:Integer.t -> t =
   fun x ~modulo -> mkintmod x modulo

  let create_assume_prime_modulus x ~modulo = mkintmod x modulo
  let lift = lift

  let inverse x =
    let x = Ctypes.(coerce gen (ptr gen) x) in
    let modulo = Ctypes.(!@(x +@ 1)) in
    let res = Ctypes.allocate t (Integer.zero ()) in
    if invmod Ctypes.(!@(x +@ 2)) modulo res = 1 then
      Some (create Ctypes.(!@res) ~modulo)
    else None

  let mul = gmul
  let pow = powgi
  let chinese : (t, [ `ROW ]) Vector.t -> t = chinese1
  let to_string = gentostr

  let get_modulo x =
    let x = Ctypes.(coerce gen (ptr gen) x) in
    gcopy Ctypes.(!@(x +@ 1))

  let order x =
    let x' = Ctypes.(coerce gen (ptr gen) x) in
    let modulo = Ctypes.(!@(x' +@ 1)) in
    znorder x (eulerphi modulo)

  let log ~base x =
    let log = znlog x base (order base) in
    let typ = gentostr (type0 log) in
    if typ = "\"t_VEC\"" then None else Some log
end

module Number_field = struct
  type t = gen
  type elt = gen

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

type 'a group_structure = {
  mul : 'a ty -> 'a ty -> 'a ty;
  pow : 'a ty -> Integer.t -> 'a ty;
  rand : unit -> 'a ty;
  hash : 'a ty -> Unsigned.ULong.t;
  equal : 'a ty -> 'a ty -> bool;
  equal_identity : 'a ty -> bool;
  bb_group : bb_group Ctypes.structure option;
}

module Fp = struct
  type t = Integer.t

  let add a b ~modulo = fp_add a b modulo
  let pow x ~exponent ~modulo = fp_pow x exponent modulo
end

module Finite_field = struct
  type t = gen

  let generator ~order =
    ff_primroot
      (ffgen order Signed.Long.zero)
      Ctypes.(coerce (ptr void) (ptr gen) null)

  let prime_field_element x ~p = ff_z_mul (ff_1 (generator ~order:p)) x

  let inj_prime_field x =
    let p = ff_to_fpxq_i x in

    if
      glength p = Signed.Long.one
      || Polynomial.degree p = 0
      || Polynomial.degree p = 1
    then Some Polynomial.(p.%[0])
    else None

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
  let sub = ff_sub
  let neg = ff_neg
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

  let fpxq_star ~p ~quotient =
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
  type 'a t = gen constraint 'a = gen
  type 'a elt = gen constraint 'a = gen

  let create ?a1 ?a2 ?a3 ?a4 ?a6 ?(dom = Ctypes.(coerce (ptr void) gen null)) ()
      =
    let to_coeff = function Some c -> c | None -> stoi Signed.Long.zero in
    let ell =
      ellinit
        (mkvec5 (to_coeff a1) (to_coeff a2) (to_coeff a3) (to_coeff a4)
           (to_coeff a6))
        dom (Signed.Long.of_int 38)
    in
    if Vector.length ell = 0 then None else Some ell

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
  let to_string ell = gentostr ell
  let weil_pairing ell ~l p q = ellweilpairing ell p q l
  let add = elladd
  let sub = ellsub
  let mul ell ~n ~p = ellmul ell p n
  let equal a b = gequal a b = 1
  let generators ell = ellgenerators ell
  let zero _ = ellinf ()

  let get_coordinates p =
    if Vector.length p = 1 then `inf
    else `point (gcopy Vector.(p.%[1]), gcopy Vector.(p.%[2]))

  let order_elt ell elt = ellorder ell elt Ctypes.(coerce (ptr void) gen null)
  let to_string_elt = gentostr

  let log ell ~base x =
    let log = elllog ell x base (order_elt ell base) in
    let typ = gentostr (type0 log) in
    if typ = "\"t_VEC\"" then None else Some log
end

let factor e =
  let m = factor e in
  let n_rows, n_cols = Matrix.dimensions m in
  assert (n_cols = 2);
  Array.init n_rows (fun i ->
      (m.Matrix.%[i + 1; 1], Integer.to_int m.Matrix.%[i + 1; 2]))

let gequal x y = gequal x y = 1

let[@inline] with_stack_clean f =
  let ltop = get_avma () in
  let res = f () in
  gerepilecopy ltop res

let[@inline] with_stack_clean_opt f =
  let ltop = get_avma () in
  match f () with
  | Some res -> Some (gerepilecopy ltop res)
  | None ->
      set_avma ltop;
      None

let[@inline] with_stack_clean6 ?av f =
  (* TODO add gc_needed check *)
  let ltop = match av with Some av -> av | None -> get_avma () in
  let a1, a2, a3, a4, a5, a6 = f () in
  let a1' = copy_bin a1 in
  let a2' = copy_bin a2 in
  let a3' = copy_bin a3 in
  let a4' = copy_bin a4 in
  let a5' = copy_bin a5 in
  let a6' = copy_bin a6 in
  set_avma ltop;
  ( bin_copy a1',
    bin_copy a2',
    bin_copy a3',
    bin_copy a4',
    bin_copy a5',
    bin_copy a6' )
