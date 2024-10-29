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
type 'a matrix = private Matrix of 'a

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

  (* *)

  let pol1_f2xx = pol1_f2xx
  let polx_f2xx = polx_f2xx
  let pol1_flxx = pol1_flxx
  let polx_flxx = polx_flxx
  let polisclass = polisclass
  let polrootsff = polrootsff
  let polteichmuller = polteichmuller
  let polhensellift = polhensellift
  let polcyclofactors = polcyclofactors
  let poliscyclo = poliscyclo
  let poliscycloprod = poliscycloprod
  let polredord = polredord
  let polred = polred
  let polred0 = polred0
  let polred2 = polred2
  let polredabs = polredabs
  let polredabs0 = polredabs0
  let polredabs2 = polredabs2
  let polredabsall = polredabsall
  let poltomonic = poltomonic
  let polcompositum0 = polcompositum0
  let poldiscfactors = poldiscfactors
  let polmod_nffix = polmod_nffix
  let polmod_nffix2 = polmod_nffix2
  let polcyclo_eval = polcyclo_eval
  let polhermite = polhermite
  let polhermite_eval0 = polhermite_eval0
  let polhermite_eval = polhermite_eval
  let pollaguerre = pollaguerre
  let pollaguerre_eval = pollaguerre_eval
  let pollaguerre_eval0 = pollaguerre_eval0
  let pollegendre = pollegendre
  let pollegendre_reduced = pollegendre_reduced
  let pollegendre_eval = pollegendre_eval
  let pollegendre_eval0 = pollegendre_eval0
  let polint = polint
  let polint_i = polint_i
  let polintspec = polintspec
  let polchebyshev = polchebyshev
  let polchebyshev_eval = polchebyshev_eval
  let polchebyshev1 = polchebyshev1
  let polchebyshev2 = polchebyshev2
  let polrecip = polrecip
  let polgalois = polgalois
  let polcoef = polcoef
  let polcoef_i = polcoef_i
  let poldegree = poldegree
  let pollead = pollead
  let polfnf = polfnf
  let poldivrem = poldivrem
  let polrootspadic = polrootspadic
  let poldisc0 = poldisc0
  let polresultant0 = polresultant0
  let polsym = polsym
  let polresultantext0 = polresultantext0
  let polresultantext = polresultantext
  let pol0_f2x = pol0_f2x
  let pol1_f2x = pol1_f2x
  let polx_f2x = polx_f2x
  let polx_zx = polx_zx
  let pol_x_powers = pol_x_powers
  let polclass = polclass
  let polmodular = polmodular
  let polmodular_zm = polmodular_zm
  let polmodular_zxx = polmodular_zxx
  let polgraeffe = polgraeffe
  let polmod_to_embed = polmod_to_embed
  let polrootsbound = polrootsbound
  let polsubcyclo = polsubcyclo
  let polsubcyclofast = polsubcyclofast
  let polzag = polzag
  let polylog0 = polylog0
  let polylogmult = polylogmult
  let polylogmult_interpolate = polylogmult_interpolate
  let pol_x = pol_x
  let pol_xn = pol_xn
  let pol_xnall = pol_xnall
  let polxn_flx = polxn_flx
  let pol_1 = pol_1
  let pol_0 = pol_0
  let pol0_flx = pol0_flx
  let pol1_flx = pol1_flx
  let polx_flx = polx_flx
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

  let nf_get_allroots = nf_get_allroots
  let nf_get_prec = nf_get_prec
  let nfmaxord_to_nf = nfmaxord_to_nf
  let nfcertify = nfcertify
  let nfgaloismatrix = nfgaloismatrix
  let nfgaloismatrixapply = nfgaloismatrixapply
  let nfgaloispermtobasis = nfgaloispermtobasis
  let nfinit_basic = nfinit_basic
  let nfinit_complete = nfinit_complete
  let nfinit0 = nfinit0
  let nfinitred = nfinitred
  let nfinitred2 = nfinitred2
  let nfisincl0 = nfisincl0
  let nfisisom = nfisisom
  let nfnewprec = nfnewprec
  let nfnewprec_shallow = nfnewprec_shallow
  let nfpoleval = nfpoleval
  let nfsplitting0 = nfsplitting0
  let nftyp = nftyp
  let nfgrunwaldwang = nfgrunwaldwang
  let nfgwkummer = nfgwkummer
  let mf_get_m = mf_get_m
  let mf_get_mindex = mf_get_mindex
  let mf_get_minv = mf_get_minv
  let mf_get_basis = mf_get_basis
  let mf_get_dim = mf_get_dim
  let mf_get_e = mf_get_e
  let mf_get_fields = mf_get_fields
  let mf_get_newforms = mf_get_newforms
  let mf_get_space = mf_get_space
  let mf_get_s = mf_get_s
  let mfcusp_get_vmjd = mfcusp_get_vmjd
  let mfnew_get_vj = mfnew_get_vj
  let nf_to_fq_init = nf_to_fq_init
  let nf_to_fq = nf_to_fq
  let nfm_to_fqm = nfm_to_fqm
  let nfv_to_fqv = nfv_to_fqv
  let nfx_to_fqx = nfx_to_fqx
  let nfx_to_monic = nfx_to_monic
  let nfbasis = nfbasis
  let nfcompositum = nfcompositum
  let nfdiscfactors = nfdiscfactors
  let nfmaxord = nfmaxord
  let nfmodpr = nfmodpr
  let nfmodprinit = nfmodprinit
  let nfmodprinit0 = nfmodprinit0
  let nfmodprlift = nfmodprlift
  let nfreducemodpr = nfreducemodpr
  let nfdisc = nfdisc
  let nf_to_scalar_or_alg = nf_to_scalar_or_alg
  let nf_to_scalar_or_basis = nf_to_scalar_or_basis
  let nf_cxlog = nf_cxlog
  let nfv_cxlog = nfv_cxlog
  let nfchecksigns = nfchecksigns
  let nfdiv = nfdiv
  let nfdiveuc = nfdiveuc
  let nfembed = nfembed
  let nfeltembed = nfeltembed
  let nfeltembed_i = nfeltembed_i
  let nfeltsign = nfeltsign
  let nffactorback = nffactorback
  let nfinv = nfinv
  let nfinvmodideal = nfinvmodideal
  let nfissquare = nfissquare
  let nfispower = nfispower
  let nflogembed = nflogembed
  let nfm_det = nfm_det
  let nfm_inv = nfm_inv
  let nfm_ker = nfm_ker
  let nfm_mul = nfm_mul
  let nfm_nfc_mul = nfm_nfc_mul
  let nfmod = nfmod
  let nfmuli = nfmuli
  let nfnorm = nfnorm
  let nfpolsturm = nfpolsturm
  let nfpow = nfpow
  let nfpow_u = nfpow_u
  let nfpowmodideal = nfpowmodideal
  let nfsign = nfsign
  let nfsign_arch = nfsign_arch
  let nfsign_from_logarch = nfsign_from_logarch
  let nfsqr = nfsqr
  let nfsqri = nfsqri
  let nfsub = nfsub
  let nftrace = nftrace
  let nfval = nfval
  let nfvalrem = nfvalrem
  let nf_get_varn = nf_get_varn
  let nf_get_pol = nf_get_pol
  let nf_get_degree = nf_get_degree
  let nf_get_r1 = nf_get_r1
  let nf_get_r2 = nf_get_r2
  let nf_get_index = nf_get_index
  let nf_get_m = nf_get_m
  let nf_get_g = nf_get_g
  let nf_get_roundg = nf_get_roundg
  let nf_get_tr = nf_get_tr
  let nf_get_diff = nf_get_diff
  let nf_get_ramified_primes = nf_get_ramified_primes
  let nf_get_roots = nf_get_roots
  let nf_get_zkprimpart = nf_get_zkprimpart
  let nf_get_zkden = nf_get_zkden
  let nf_get_invzk = nf_get_invzk
  let nfv_to_scalar_or_alg = nfv_to_scalar_or_alg
  let nfc_multable_mul = nfc_multable_mul
  let nfc_nf_mul = nfc_nf_mul
  let nf_get_gtwist = nf_get_gtwist
  let nf_get_gtwist1 = nf_get_gtwist1
  let nf_to_fp_coprime = nf_to_fp_coprime
  let nfdetint = nfdetint
  let nfdivmodpr = nfdivmodpr
  let nfhnf = nfhnf
  let nfhnf0 = nfhnf0
  let nfhnfmod = nfhnfmod
  let nfkermodpr = nfkermodpr
  let nfmulmodpr = nfmulmodpr
  let nfpowmodpr = nfpowmodpr
  let nfreduce = nfreduce

  let smith_normal_form nf sipm =
    (*square integral invertible pseudo-matrix sipm*)
    let d, _ = Matrix.dimensions sipm in
    let sipm =
      Vector.of_array
        [|
          sipm;
          Vector.init d ~f:(fun _ -> Integer.of_int 1);
          Vector.init d ~f:(fun _ -> Integer.of_int 1);
        |]
    in
    let res = nfsnf0 nf sipm Signed.Long.one in
    ( Vector.(res.%[1]),
      Vector.(res.%[2]) |> matbasistoalg nf,
      Vector.(res.%[3]) |> matbasistoalg nf )

  let nfsnf0 = nfsnf0
  let nfsolvemodpr = nfsolvemodpr
  let nfsubfields = nfsubfields
  let nfsubfields0 = nfsubfields0
  let nfsubfieldscm = nfsubfieldscm
  let nfsubfieldsmax = nfsubfieldsmax
  let nflist = nflist
  let nfresolvent = nfresolvent
  let nf_pv_to_prv = nf_pv_to_prv
  let nf_rnfeq = nf_rnfeq
  let nf_rnfeqsimple = nf_rnfeqsimple
  let nf_nfzk = nf_nfzk
  let nfeltup = nfeltup
  let nfislocalpower = nfislocalpower
  let nf_cxlog_normalize = nf_cxlog_normalize
  let nfcyclotomicunits = nfcyclotomicunits
  let nfsign_units = nfsign_units
  let nfsign_tu = nfsign_tu
  let nfsign_fu = nfsign_fu
  let nf_deg1_prime = nf_deg1_prime
  let nfarchstar = nfarchstar
  let nf_hyperell_locally_soluble = nf_hyperell_locally_soluble
  let nfhilbert = nfhilbert
  let nfhilbert0 = nfhilbert0
  let nfhyperellpadicfrobenius = nfhyperellpadicfrobenius
  let nffactor = nffactor
  let nffactormod = nffactormod
  let nfgcd = nfgcd
  let nfgcd_all = nfgcd_all
  let nfissquarefree = nfissquarefree
  let nfroots = nfroots
  let nfroots_if_split = nfroots_if_split
  let nfrootsof1 = nfrootsof1
  let nfrootsq = nfrootsq
  let nfx_disc = nfx_disc
  let nfx_resultant = nfx_resultant

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

  (*  *)

  let ellanalyticrank ell =
    let prec = Signed.Long.of_int 10 in
    let eps = Ctypes.(coerce (ptr void) gen null) in
    ellanalyticrank ell eps prec

  let ellanalyticrank_bitprec = ellanalyticrank_bitprec
  let ellanal_globalred_all = ellanal_globalred_all
  let ellsupersingularj_fpxq = ellsupersingularj_fpxq
  let elltrace_extension = elltrace_extension
  let ellheegner = ellheegner
  let elll1 = elll1
  let elll1_bitprec = elll1_bitprec
  let ellconvertname = ellconvertname
  let elldatagenerators = elldatagenerators
  let ellidentify = ellidentify
  let ellsearch = ellsearch
  let ellsearchcurve = ellsearchcurve
  let ell_is_integral = ell_is_integral
  let ellq_get_cm = ellq_get_cm
  let ellq_get_n = ellq_get_n
  let ellq_get_nfa = ellq_get_nfa
  let ellqp_tate_uniformization = ellqp_tate_uniformization
  let ellqp_agm = ellqp_agm
  let ellqp_u = ellqp_u
  let ellqp_u2 = ellqp_u2
  let ellqp_q = ellqp_q
  let ellqp_ab = ellqp_ab
  let ellqp_l = ellqp_l
  let ellqp_root = ellqp_root
  let ellqtwist_bsdperiod = ellqtwist_bsdperiod
  let ellr_area = ellr_area
  let ellr_ab = ellr_ab
  let ellr_eta = ellr_eta
  let ellr_omega = ellr_omega
  let ellr_roots = ellr_roots
  let ellan = ellan
  let ellanq_zv = ellanq_zv
  let ellanal_globalred = ellanal_globalred
  let ellap = ellap
  let ellap_cm_fast = ellap_cm_fast
  let ellbasechar = ellbasechar
  let ellbsd = ellbsd
  let ellchangecurve = ellchangecurve
  let ellchangeinvert = ellchangeinvert
  let ellchangepoint = ellchangepoint
  let ellchangepointinv = ellchangepointinv
  let elleisnum = elleisnum
  let elleta = elleta
  let elleulerf = elleulerf
  let ellff_get_card = ellff_get_card
  let ellff_get_gens = ellff_get_gens
  let ellff_get_group = ellff_get_group
  let ellff_get_o = ellff_get_o
  let ellff_get_p = ellff_get_p
  let ellff_get_m = ellff_get_m
  let ellff_get_d = ellff_get_d
  let ellfromj = ellfromj
  let ellgenerators = ellgenerators
  let ellglobalred = ellglobalred
  let ellgroup = ellgroup
  let ellgroup0 = ellgroup0
  let ellheight0 = ellheight0
  let ellheight = ellheight
  let ellheightmatrix = ellheightmatrix
  let ellheightoo = ellheightoo
  let ellintegralmodel = ellintegralmodel
  let ellintegralmodel_i = ellintegralmodel_i
  let elliscm = elliscm
  let ellisoncurve = ellisoncurve
  let ellisotree = ellisotree
  let ellissupersingular = ellissupersingular
  let elljissupersingular = elljissupersingular
  let elllseries = elllseries
  let elllocalred = elllocalred
  let elllog = elllog
  let ellminimaldisc = ellminimaldisc
  let ellminimalmodel = ellminimalmodel
  let ellminimaltwist = ellminimaltwist
  let ellminimaltwist0 = ellminimaltwist0
  let ellminimaltwistcond = ellminimaltwistcond
  let ellnf_vecarea = ellnf_vecarea
  let ellnf_veceta = ellnf_veceta
  let ellnf_vecomega = ellnf_vecomega
  let ellneg = ellneg
  let ellorder = ellorder
  let ellorder_q = ellorder_q
  let ellordinate = ellordinate
  let ellpadicheight0 = ellpadicheight0
  let ellpadicheightmatrix = ellpadicheightmatrix
  let ellperiods = ellperiods
  let ellrootno = ellrootno
  let ellrootno_global = ellrootno_global
  let ellsaturation = ellsaturation
  let ellsea = ellsea
  let ellsigma = ellsigma
  let ellsupersingularj = ellsupersingularj
  let elltamagawa = elltamagawa
  let elltaniyama = elltaniyama
  let elltatepairing = elltatepairing
  let elltors = elltors
  let elltors0 = elltors0
  let elltors_psylow = elltors_psylow
  let elltrace = elltrace
  let elltwist = elltwist
  let ellwp = ellwp
  let ellwp0 = ellwp0
  let ellwpseries = ellwpseries
  let ellxn = ellxn
  let ellzeta = ellzeta
  let ellfromeqn = ellfromeqn
  let ellformaldifferential = ellformaldifferential
  let ellformalexp = ellformalexp
  let ellformallog = ellformallog
  let ellformalpoint = ellformalpoint
  let ellformalw = ellformalw
  let ellnonsingularmultiple = ellnonsingularmultiple
  let ellpadicl = ellpadicl
  let ellpadicbsd = ellpadicbsd
  let ellpadicfrobenius = ellpadicfrobenius
  let ellpadicheight = ellpadicheight
  let ellpadiclog = ellpadiclog
  let ellpadicregulator = ellpadicregulator
  let ellpadics2 = ellpadics2
  let ell2cover = ell2cover
  let ellrank = ellrank
  let ellrankinit = ellrankinit
  let ellisdivisible = ellisdivisible
  let ellisogenyapply = ellisogenyapply
  let ellisogeny = ellisogeny
  let ellisomat = ellisomat
  let ellweilcurve = ellweilcurve
  let ell_get_a1 = ell_get_a1
  let ell_get_a2 = ell_get_a2
  let ell_get_a3 = ell_get_a3
  let ell_get_a4 = ell_get_a4
  let ell_get_a6 = ell_get_a6
  let ell_get_b2 = ell_get_b2
  let ell_get_b4 = ell_get_b4
  let ell_get_b6 = ell_get_b6
  let ell_get_b8 = ell_get_b8
  let ell_get_c4 = ell_get_c4
  let ell_get_c6 = ell_get_c6
  let ell_get_type = ell_get_type
  let ellff_get_field = ellff_get_field
  let ellff_get_a4a6 = ellff_get_a4a6
  let ellqp_get_zero = ellqp_get_zero
  let ellqp_get_prec = ellqp_get_prec
  let ellqp_get_p = ellqp_get_p
  let ellr_get_prec = ellr_get_prec
  let ellr_get_sign = ellr_get_sign
  let ellnf_get_nf = ellnf_get_nf
  let ellnf_get_bnf = ellnf_get_bnf
  let ellmodulareqn = ellmodulareqn
  let ellmoddegree = ellmoddegree
  let ellratpoints = ellratpoints
  let elle = elle
  let ellk = ellk
  let ellpadiclambdamu = ellpadiclambdamu
  let ell_is_inf = ell_is_inf
end

let isdiagonal e = isdiagonal e = 1

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
