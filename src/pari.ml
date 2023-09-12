include Pari_bindings

type _ t = gen

let t = gen

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
  let qflll0 = qflll0

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
  let singleton x = mkvec (Obj.magic x)
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
      if typ = "\"t_POL\"" then Signed.Long.succ (gvar (fst (List.nth p 0)))
      else Signed.Long.zero
    in
    Ctypes.(p' +@ 1 <-@ Signed.Long.(shift_left v Int.(sub 64 @@ sub 2 16)));
    for i = 2 to Int.succ len do
      Ctypes.(p'' +@ i <-@ zero)
    done;
    List.iter (fun (c, idx) -> Ctypes.(p'' +@ Int.add idx 2 <-@ c)) p;
    normalizepol_lg p' size

  let var t = function Some var -> Signed.Long.(of_int var) | None -> gvar t
  let deriv ?indeterminate t = deriv t (var t indeterminate)

  let derivn ?indeterminate t n =
    derivn t (Signed.Long.of_int n) (var t indeterminate)

  let cyclotomic n = polcyclo n (Signed.Long.of_int 0)
  let is_irreducible p = polisirreducible p = Signed.Long.one
  let minimal p = polredbest p (Signed.Long.of_int 0)
end

module NumberField = struct
  type number_field
  type nonrec t = number_field t

  let create p = nfinit p Signed.Long.(of_int 4)
  let are_isomorphic a b = gequal0 (nfisincl a b) <> 1

  let sign nf =
    let a = Ctypes.(allocate long Signed.Long.zero) in
    let b = Ctypes.(allocate long Signed.Long.zero) in
    nf_get_sign nf a b;
    Ctypes.(!@a, !@b)

  let discriminant nf = nf_get_disc nf

  let splitting = function
    | `F nf -> nfsplitting nf Ctypes.(coerce (ptr void) gen null)
    | `P p -> nfsplitting p Ctypes.(coerce (ptr void) gen null)
end
