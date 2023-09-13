include Pari_bindings

type _ t = gen

let t = gen

type 'a ring = gen
type 'a field = gen

let register_gc =
  Gc.finalise (fun v -> pari_free Ctypes.(coerce gen (ptr void) v))

module Complex = struct
  type complex
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
  type real
  type t = gen

  external inj_complex : t -> Complex.t = "%identity"

  let of_signed = stor
  let shift = mpshift
  let sqrt = sqrtr
  let ceil = gceil
  let add = gadd
  let inv = ginv
end

module Rational = struct
  type rational
  type t = gen
  type ring = gen

  external inj_ring : t -> ring = "%identity"
  external inj_real : t -> Real.t = "%identity"
  external inj_complex : t -> Complex.t = "%identity"

  let shift = mpshift
end

module Integer = struct
  type integer
  type t = gen

  external inj_rat : t -> Rational.t = "%identity"
  external inj_real : t -> Rational.t = "%identity"
  external inj_complex : t -> Rational.t = "%identity"

  let equal x y = equalii x y == 1
  let of_int i = stoi (Signed.Long.of_int i)
  let of_signed = stoi
  let shift = mpshift
  let sqrt i = inj_complex (sqrti i)
  let mul = mulii
  let add = addii
  let sub = subii
  let neg = negi
  let modulo = modii
  let of_string s = (* todo error handling *) Some (gp_read_str s)

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
  type matrix
  type nonrec 'a t = gen constraint 'a = gen

  let id n =
    let m = matid n in
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
  let ( .%[;..] ) m i = Obj.magic @@ !@(getcoeff m i)
  let ( .%[;..]<- ) m i v = Ctypes.(getcoeff m i <-@ Obj.magic @@ v)
  let inj x ~inj:_ = x
end

module Set = struct
  type nonrec 'a t = gen constraint 'a = gen

  let length (s : 'a t) = glength s
  let search (s : 'a t) (e : 'a) = setsearch s (Obj.magic e)
end

module Vector = struct
  type ('a, 'b) t = gen constraint 'a = gen constraint 'b = [< `COL | `ROW ]

  let length = glength

  let ( .%[] ) m i =
    Obj.magic @@ Ctypes.( !@ ) (safegel m (Signed.Long.of_int i))

  let ( .%[]<- ) m i v =
    Ctypes.(safegel m (Signed.Long.of_int i) <-@ Obj.magic v)

  let of_array a =
    let l = Array.length a in
    let v = zerovec_block (Signed.Long.of_int l) in
    for i = 1 to l do
      v.%[i] <- a.(i - 1)
    done;
    register_gc v;
    v

  let equal x y = gequal x y == 1
  let slice m ~start ~stop = vecslice m start stop
  let mul x y = Obj.magic @@ gmul x y
  let add = gadd
  let sub = gsub
  let neg = gneg
  let transpose_row = gtrans
  let transpose_column = gtrans
  let to_set = gtoset

  let singleton x =
    let s = mkvec (Obj.magic x) in
    register_gc s;
    s

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
  let one = stoi (Signed.Long.of_int 1)
  let mul = gmul
  let equal = equal
  let zero = stoi (Signed.Long.of_int 0)
  let add = gadd
  let sub = gsub
  let neg = gneg
  let eval p x = Obj.magic (poleval p (Obj.magic x))
  let degree t = degree t

  let get_coeff t i =
    Obj.magic Ctypes.(!@(coerce gen (ptr gen) t +@ Int.add i 2))

  let create (p : ('a * int) list) : 'a t =
    let len =
      succ
        (List.fold_left
           (fun acc (_, i) -> if i > acc then i else acc)
           (Int.neg 1) p)
    in
    let size = Signed.Long.of_int (Int.add len 2) in
    let p' = cgetg size (Signed.Long.of_int 10 (* t_POL *)) in
    let p'' = Ctypes.(coerce gen (ptr gen) p') in
    (* TODO: support 32-bit arch *)
    let typ = gentostr (type0 (fst (List.nth p 0))) in
    let v =
      if typ = "\"t_POL\"" then Signed.Long.(succ (gvar (fst (List.nth p 0))))
      else Signed.Long.zero
    in
    Ctypes.(p' +@ 1 <-@ Signed.Long.(shift_left v Stdlib.(64 - 2 - 16)));
    for i = 2 to Int.succ len do
      Ctypes.(p'' +@ i <-@ zero)
    done;
    List.iter (fun (c, idx) -> Ctypes.(p'' +@ Int.add idx 2 <-@ c)) p;
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
end

module Number_field = struct
  type number_field
  type 'a p = 'a t
  type nonrec t = number_field t
  type elt

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
  let z_basis nf : (elt p, [ `ROW ]) Vector.t = nf_get_zk nf
  let elt a : elt p = Vector.(transpose_row (of_array a))
  let add nf a b = nfadd nf a b
  let mul nf a b = nfmul nf a b

  let divrem nf a b =
    let qr = nfdivrem nf a b in
    Vector.(qr.%[1], qr.%[2])

  let ideal_norm nf a = idealnorm nf a

  let splitting = function
    | `F nf -> nfsplitting nf Ctypes.(coerce (ptr void) gen null)
    | `P p -> nfsplitting p Ctypes.(coerce (ptr void) gen null)

  module Infix = struct
    let ( = ) a b = gequal a b = 1
  end
end

type 'a group

type 'a group_structure = {
  mul : 'a group t -> 'a group t -> 'a group t;
  pow : 'a group t -> Integer.t -> 'a group t;
  rand : unit -> 'a group t;
  hash : 'a group t -> Unsigned.ULong.t;
  equal : 'a group t -> 'a group t -> bool;
  equal_identity : 'a group t -> bool;
  bb_group : bb_group Ctypes.structure option;
}

module Fp = struct
  type t = Integer.t

  let add a b ~modulo = fp_add a b modulo
  let pow x ~exponent ~modulo = fp_pow x exponent modulo
end

module Finite_field = struct
  type finite_field
  type ring = gen
  type nonrec t = finite_field t

  external inj_ring : t -> ring = "%identity"

  let generator ~order =
    ff_primroot
      (ffgen (stoi (Signed.Long.of_int order)) (Signed.Long.of_int 0))
      Ctypes.(coerce (ptr void) (ptr gen) null)

  let create ~p ~degree =
    ffinit (Integer.of_int p) (Signed.Long.of_int degree) Signed.Long.zero

  let fpxq_star ~(p : pari_ulong) ~(quotient : Fp.t Polynomial.t) :
      finite_field group_structure =
    let open Ctypes in
    let q =
      powuu p
        (Unsigned.ULong.of_int
           (Signed.Long.to_int (Polynomial.degree quotient)))
    in
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
end
