type pari_ulong = Unsigned.ULong.t

val pari_ulong : pari_ulong Ctypes.typ

type 'a ty

val t : 'a ty Ctypes.typ

type byteptr = Unsigned.uchar Ctypes.ptr

val byteptr : byteptr Ctypes.typ

type pari_sp = pari_ulong

val pari_sp : pari_sp Ctypes.typ

type pari_logstyles =
  | LOGSTYLE_NONE
  | LOGSTYLE_PLAIN
  | LOGSTYLE_COLOR
  | LOGSTYLE_TEX

type err_list =
  | E_SYNTAX
  | E_BUG
  | E_ALARM
  | E_FILE
  | E_MISC
  | E_FLAG
  | E_IMPL
  | E_ARCH
  | E_PACKAGE
  | E_NOTFUNC
  | E_PREC
  | E_TYPE
  | E_DIM
  | E_VAR
  | E_PRIORITY
  | E_USER
  | E_STACK
  | E_STACKTHREAD
  | E_OVERFLOW
  | E_DOMAIN
  | E_COMPONENT
  | E_MAXPRIME
  | E_CONSTPOL
  | E_IRREDPOL
  | E_COPRIME
  | E_PRIME
  | E_MODULUS
  | E_ROOTS0
  | E_OP
  | E_TYPE2
  | E_INV
  | E_MEM
  | E_SQRTN
  | E_FILEDESC
  | E_NONE

type pari_timer
type pari_str
type pari_sieve
type forprime_t
type forcomposite_t
type forvec_t
type forpart_t
type forperm_t
type forsubset_t
type pari_plot
type genbin
type pari_mainstack
type entree
type pari_parsestate
type pari_compilestate
type pari_mtstate
type pari_evalstate
type pari_varstate
type pari_global_state
type pari_thread
type mt_state
type pari_mt
type parfor_iter
type parfor_t
type parforeach_t
type parforprime_t
type parforvec_t
type hashentry
type hashtable
type gp_path
type pariout_t
type nfmaxord_t
type qfr_data
type fp_chk_fun
type zlog_s
type bb_group
type bb_field
type bb_algebra
type bb_ring
type ring = private Ring
type field = private Field
type unique_factorization_domain = private Unique_factorization_domain
type complex = private Complex
type real = private Real
type rational = private Rational
type integer = private Integer
type 'a polynomial = private Polynomial of 'a
type integer_mod = private Integer_mod
type finite_field = private Finite_field
type number_field = private Number_field
type 'a elliptic_curve = private Elliptic_curve of 'a

val factor : 'a ty -> ('a ty * int) array

module rec Complex : sig
  type t = complex ty

  val inv : t -> t
  val add : t -> t -> t
  val create : re:Real.t -> im:Real.t -> t
  val to_string : t -> string
end

and Real : sig
  type t = real ty

  val inj_complex : t -> Complex.t
  (** {@ocaml[
    # let x1 = Integer.(inj_real (of_int 1));;
    val x1 : Real.t = <abstr>
    # let x2 = Integer.(inj_real (of_int 1));;
    val x2 : Real.t = <abstr>
    # let a = Pari.Complex.create  ~re:x1 ~im:x2;;
    val a : Pari.Complex.t = <abstr>
    # Complex.add (Real.inj_complex x1) a |> Complex.to_string;;
    - : string = "2 + I"
    ]} *)

  val of_signed : Signed.long -> Signed.long -> t
  val shift : t -> int -> t
  val sqrt : t -> t
  val ceil : t -> Integer.t
  val add : t -> t -> t
  val inv : t -> t
end

and Rational : sig
  type t = rational ty

  val of_int : int -> t
  val inj_real : t -> Real.t
  val inj_complex : t -> Complex.t
  val shift : t -> int -> t
end

and Integer : sig
  type t = integer ty

  val inj_rat : t -> Rational.t
  val inj_real : t -> Real.t
  val inj_complex : t -> Complex.t
  val of_int : int -> t
  val to_int : t -> int
  val of_hex : Hex.t -> t
  val of_signed : Signed.long -> t
  val equal : t -> t -> bool
  val shift : t -> int -> t
  val sqrt : t -> Real.t
  val zero : unit -> t
  val mul : t -> t -> t
  val add : t -> t -> t
  val sub : t -> t -> t
  val neg : t -> t
  val pow : t -> t -> t
  val modulo : t -> t -> t
  val of_string : string -> t option
  val to_string : t -> string
  val random_prime : bits_amount:int -> t
  val gcdext : t -> t -> t * t * t
  val gcd : t -> t -> t
  val divexact : t -> t -> t
  val random : t -> t

  module Infix : sig
    val ( * ) : t -> t -> t
    val ( + ) : t -> t -> t
    val ( - ) : t -> t -> t
    val ( ~- ) : t -> t
    val ( mod ) : t -> t -> t
    val ( = ) : t -> t -> bool
  end
end

type group = private Group

type 'a group_structure = {
  mul : 'a ty -> 'a ty -> 'a ty;
  pow : 'a ty -> Integer.t -> 'a ty;
  rand : unit -> 'a ty;
  hash : 'a ty -> Unsigned.ULong.t;
  equal : 'a ty -> 'a ty -> bool;
  equal_identity : 'a ty -> bool;
  bb_group : bb_group Ctypes.structure option;
}

module Set : sig
  type 'a t constraint 'a = 'b ty

  val length : 'a t -> Signed.Long.t
  val search : 'a t -> 'a -> Signed.Long.t -> Signed.Long.t
end

module Vector : sig
  type ('a, 'b) t constraint 'a = 'c ty constraint 'b = [< `COL | `ROW ]

  val length : ('a, 'b) t -> int
  val of_array : 'a array -> ('a, [ `ROW ]) t
  val array_map : f:('a -> 'b ty) -> 'a array -> ('b ty, [ `ROW ]) t
  val init : int -> f:(int -> 'a ty) -> ('a ty, [ `ROW ]) t
  val equal : ('a, 'b) t -> ('a, 'b) t -> bool
  val slice : ('a, 'b) t -> start:int -> stop:int -> ('a, 'b) t
  val ( .%[] ) : ('a, 'b) t -> int -> 'a
  val ( .%[]<- ) : ('a, 'b) t -> int -> 'a -> unit
  val mul : ('a, [ `ROW ]) t -> ('a, [ `COL ]) t -> 'a
  val add : ('a, 'b) t -> ('a, 'b) t -> ('a, 'b) t
  val sub : ('a, 'b) t -> ('a, 'b) t -> ('a, 'b) t
  val neg : ('a, 'b) t -> ('a, 'b) t
  val transpose_row : ('a, [ `ROW ]) t -> ('a, [ `COL ]) t
  val transpose_column : ('a, [ `COL ]) t -> ('a, [ `ROW ]) t
  val to_set : ('a, 'b) t -> 'a Set.t
  val singleton : 'a -> ('a, 'b) t
  val concat : ('a, 'b) t -> ('a, 'b) t -> ('a, 'b) t
  val inj : ('a, 'b) t -> inj:('a -> 'c) -> ('c, 'b) t
  val to_string : ('a, 'b) t -> string

  module Infix : sig
    val ( = ) : ('a, 'b) t -> ('a, 'b) t -> bool
    val ( * ) : ('a, [ `COL ]) t -> ('a, [ `ROW ]) t -> 'a
    val ( + ) : ('a, 'b) t -> ('a, 'b) t -> ('a, 'b) t
    val ( - ) : ('a, 'b) t -> ('a, 'b) t -> ('a, 'b) t
    val ( ~- ) : ('a, 'b) t -> ('a, 'b) t
  end
end

module Matrix : sig
  type 'a t constraint 'a = 'b ty

  val dimensions : 'a t -> int * int
  val id : int -> Integer.t t
  val inv : 'a t -> 'a ty
  val mul : 'a t -> 'a t -> 'a t
  val lll : 'a t -> 'a t
  val ( .%[] ) : 'a t -> int -> ('a, [ `COL ]) Vector.t
  val ( .%[]<- ) : 'a t -> int -> ('a, [ `COL ]) Vector.t -> unit
  val ( .%[;..] ) : 'a t -> int array -> 'a t
  val ( .%[;..]<- ) : 'a t -> int array -> 'a -> unit
  val inj : 'a t -> inj:('a -> 'b) -> 'b t
end

module Polynomial : sig
  type 'a t = 'a polynomial ty constraint 'a = 'b ty

  val to_string : 'a t -> string
  val mul : 'a t -> 'a t -> 'a t
  val div : 'a t -> 'a t -> 'a t
  val equal : 'a t -> 'a t -> bool
  val add : 'a t -> 'a t -> 'a t
  val sub : 'a t -> 'a t -> 'a t
  val neg : 'a t -> 'a t
  val eval : 'a t -> 'a -> 'a
  val degree : 'a t -> int

  val create : 'a array -> 'a t
  (** [create a] returns {m a_{n-1} X^{n-1} + ... + a_0} for array [a] of length [n].

      {@ocaml[
      # let q = Polynomial.create
        [|
          Integer.of_int 1;
          Integer.of_int (-111);
          Integer.of_int 6064;
          Integer.of_int (-189804);
        |];;
      val q : Integer.t Polynomial.t = <abstr>
      # Polynomial.to_string q;;
      - : string = "x^3 - 111*x^2 + 6064*x - 189804"
      # let zero = Polynomial.create [| Integer.of_int 0 |];;
      val zero : Integer.t Polynomial.t = <abstr>
      # let qq = Polynomial.create [| q; q; zero; zero |];;
      val qq : Integer.t Polynomial.t Polynomial.t = <abstr>
      # Polynomial.to_string qq;;
      - : string =
      "(x^3 - 111*x^2 + 6064*x - 189804)*y^3 + (x^3 - 111*x^2 + 6064*x - 189804)*y^2"
      ]}
  *)

  val of_string : string -> 'a t option
  val deriv : ?indeterminate:int -> 'a t -> 'a t
  val derivn : ?indeterminate:int -> 'a t -> int -> 'a t
  val cyclotomic : Signed.long -> Integer.t t
  val is_irreducible : 'a t -> bool

  val minimal : 'a t -> 'a t
  (** [minimal p] reduces [p] to be the minimal polynomial of the
      roots of [p] over the field of the rational numbers.
   
      {@ocaml[
      # let q = Polynomial.create
        [|
          (Rational.of_int 1);
          (Rational.of_int (-111));
          (Rational.of_int 6064);
          (Rational.of_int (-189804));
        |];;
      val q : Rational.t Polynomial.t = <abstr>
      # Polynomial.to_string q;;
      - : string = "x^3 - 111*x^2 + 6064*x - 189804"
      # let qmin = Polynomial.minimal q;;
      val qmin : Rational.t Polynomial.t = <abstr>
      # Polynomial.to_string qmin;;
      - : string = "x^3 - x^2 - 60*x - 364"
      # Number_field.(are_isomorphic (create q) (create qmin));
      - : bool = true
      ]} *)

  val ( .%[] ) : 'a t -> int -> 'a
  val roots_ff : finite_field ty t -> (finite_field ty, [ `ROW ]) Vector.t
  val fold_left : f:('b -> 'a -> 'a) -> acc:'a -> 'b t -> 'a
  val fold_left2 : f:('a -> 'b -> 'c -> 'c) -> acc:'c -> 'a t -> 'b t -> 'c

  val fold_left2_vec :
    f:('a -> 'b -> 'c -> 'c) -> acc:'c -> 'a t -> ('b, _) Vector.t -> 'c

  val inj_base_ring : inj:('a -> 'b) -> 'a t -> 'b t

  (* *)

  val pol1_f2xx : Signed.long -> Signed.long -> 'a ty
  val polx_f2xx : Signed.long -> Signed.long -> 'a ty
  val pol1_flxx : Signed.long -> Signed.long -> 'a ty
  val polx_flxx : Signed.long -> Signed.long -> 'a ty
  val polisclass : 'a ty -> Signed.long
  val polrootsff : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val polteichmuller : 'a ty -> pari_ulong -> Signed.long -> 'a ty
  val polhensellift : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
  val polcyclofactors : 'a ty -> 'a ty
  val poliscyclo : 'a ty -> Signed.long
  val poliscycloprod : 'a ty -> Signed.long
  val polredord : 'a ty -> 'a ty
  val polred : 'a ty -> 'a ty
  val polred0 : 'a ty -> Signed.long -> 'a ty -> 'a ty
  val polred2 : 'a ty -> 'a ty
  val polredabs : 'a ty -> 'a ty
  val polredabs0 : 'a ty -> Signed.long -> 'a ty
  val polredabs2 : 'a ty -> 'a ty
  val polredabsall : 'a ty -> Signed.long -> 'a ty
  val poltomonic : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
  val polcompositum0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val poldiscfactors : 'a ty -> Signed.long -> 'a ty
  val polmod_nffix : string -> 'a ty -> 'a ty -> int -> 'a ty
  val polmod_nffix2 : string -> 'a ty -> 'a ty -> 'a ty -> int -> 'a ty
  val polcyclo_eval : Signed.long -> 'a ty -> 'a ty
  val polhermite : Signed.long -> Signed.long -> 'a ty
  val polhermite_eval0 : Signed.long -> 'a ty -> Signed.long -> 'a ty
  val polhermite_eval : Signed.long -> 'a ty -> 'a ty
  val pollaguerre : Signed.long -> 'a ty -> Signed.long -> 'a ty
  val pollaguerre_eval : Signed.long -> 'a ty -> 'a ty -> 'a ty
  val pollaguerre_eval0 : Signed.long -> 'a ty -> 'a ty -> Signed.long -> 'a ty
  val pollegendre : Signed.long -> Signed.long -> 'a ty
  val pollegendre_reduced : Signed.long -> Signed.long -> 'a ty
  val pollegendre_eval : Signed.long -> 'a ty -> 'a ty
  val pollegendre_eval0 : Signed.long -> 'a ty -> Signed.long -> 'a ty
  val polint : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

  val polint_i :
    'a ty -> 'a ty -> 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty

  val polintspec :
    'a ty ->
    'a ty ->
    'a ty ->
    Signed.long ->
    Signed.long Ctypes_static.ptr ->
    'a ty

  val polchebyshev : Signed.long -> Signed.long -> Signed.long -> 'a ty
  val polchebyshev_eval : Signed.long -> Signed.long -> 'a ty -> 'a ty
  val polchebyshev1 : Signed.long -> Signed.long -> 'a ty
  val polchebyshev2 : Signed.long -> Signed.long -> 'a ty
  val polrecip : 'a ty -> 'a ty
  val polgalois : 'a ty -> Signed.long -> 'a ty
  val polcoef : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val polcoef_i : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val poldegree : 'a ty -> Signed.long -> Signed.long
  val pollead : 'a ty -> Signed.long -> 'a ty
  val polfnf : 'a ty -> 'a ty -> 'a ty
  val poldivrem : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
  val polrootspadic : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val poldisc0 : 'a ty -> Signed.long -> 'a ty
  val polresultant0 : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
  val polsym : 'a ty -> Signed.long -> 'a ty
  val polresultantext0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val polresultantext : 'a ty -> 'a ty -> 'a ty
  val pol0_f2x : Signed.long -> 'a ty
  val pol1_f2x : Signed.long -> 'a ty
  val polx_f2x : Signed.long -> 'a ty
  val polx_zx : Signed.long -> 'a ty
  val pol_x_powers : Signed.long -> Signed.long -> 'a ty
  val polclass : 'a ty -> Signed.long -> Signed.long -> 'a ty

  val polmodular :
    Signed.long -> Signed.long -> 'a ty -> Signed.long -> Signed.long -> 'a ty

  val polmodular_zm : Signed.long -> Signed.long -> 'a ty

  val polmodular_zxx :
    Signed.long -> Signed.long -> Signed.long -> Signed.long -> 'a ty

  val polgraeffe : 'a ty -> 'a ty
  val polmod_to_embed : 'a ty -> Signed.long -> 'a ty
  val polrootsbound : 'a ty -> 'a ty -> 'a ty
  val polsubcyclo : Signed.long -> Signed.long -> Signed.long -> 'a ty

  val polsubcyclofast :
    'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty

  val polzag : Signed.long -> Signed.long -> 'a ty
  val polylog0 : Signed.long -> 'a ty -> Signed.long -> Signed.long -> 'a ty
  val polylogmult : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val polylogmult_interpolate : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
  val pol_x : Signed.long -> 'a ty
  val pol_xn : Signed.long -> Signed.long -> 'a ty
  val pol_xnall : Signed.long -> Signed.long -> 'a ty
  val polxn_flx : Signed.long -> Signed.long -> 'a ty
  val pol_1 : Signed.long -> 'a ty
  val pol_0 : Signed.long -> 'a ty
  val pol0_flx : Signed.long -> 'a ty
  val pol1_flx : Signed.long -> 'a ty
  val polx_flx : Signed.long -> 'a ty
end

module Fp : sig
  type t = Integer.t

  val add : t -> t -> modulo:t -> t
  val pow : t -> exponent:t -> modulo:t -> t
end

module Finite_field : sig
  type t = finite_field ty

  val generator : order:Integer.t -> t
  val prime_field_element : Integer.t -> p:Integer.t -> t
  val inj_prime_field : t -> Fp.t option
  val finite_field_element : Integer.t array -> t -> t

  val create : p:int -> degree:int -> finite_field ty Polynomial.t
  (** [create p degree] returns a monic irreducible polynomial of the
      given [degree] over F_p[X]. *)

  val generator_from_irreducible_polynomial : finite_field ty Polynomial.t -> t
  val residue_class : t -> finite_field ty Polynomial.t
  val equal : t -> t -> bool
  val add : t -> t -> t
  val sub : t -> t -> t
  val neg : t -> t
  val mul : t -> t -> t
  val pow : t -> Integer.t -> t
  val random : t -> t
  val zero : t -> t
  val one : t -> t

  val extend :
    finite_field ty ->
    [< `Degree of int | `Quotient of finite_field ty Polynomial.t ] ->
    t

  val fpxq_star :
    p:pari_ulong -> quotient:Fp.t Polynomial.t -> finite_field group_structure

  val to_string : t -> string

  module Infix : sig
    val ( ~- ) : t -> t
    val ( + ) : t -> t -> t
    val ( - ) : t -> t -> t
    val ( * ) : t -> t -> t
    val ( ^ ) : t -> Integer.t -> t
  end
end

module Integer_mod : sig
  type t = integer_mod ty

  val inj_group : t -> integer_mod ty
  val create : Integer.t -> modulo:Integer.t -> t

  val create_assume_prime_modulus :
    Integer.t -> modulo:Integer.t -> integer_mod ty

  val lift : integer_mod ty -> Integer.t
  val inverse : integer_mod ty -> integer_mod ty option
  val mul : integer_mod ty -> integer_mod ty -> integer_mod ty
  val pow : integer_mod ty -> Integer.t -> integer_mod ty
  val chinese : (t, [ `ROW ]) Vector.t -> t
  val to_string : integer_mod ty -> string
  val get_modulo : integer_mod ty -> Integer.t
  val order : integer_mod ty -> Integer.t
  val log : base:integer_mod ty -> integer_mod ty -> Integer.t option
end

module Number_field : sig
  type t
  type elt = number_field ty

  val create : Rational.t Polynomial.t -> t
  (** [create p] returns the number field Q(X)/(p) for a monic
      irreducible polynomial [p] over the field Q of the rationals. *)

  val are_isomorphic : t -> t -> bool
  (** [are_isomorphic a b] returns true if and only if number
      fields [a] and [b] are isomorphic.

      {@ocaml[
      # let q =
        Polynomial.create
          [|
            Rational.of_int 1;
            Rational.of_int (-111);
            Rational.of_int 6064;
            Rational.of_int (-189804);
          |];;
      val q : Rational.t Polynomial.t = <abstr>
      # let zero = Polynomial.create [| Rational.of_int 0 |];;
      val zero : Rational.t Polynomial.t = <abstr>
      # let qq = Polynomial.create [| q; q; zero; zero |];;
      val qq : Rational.t Polynomial.t Polynomial.t = <abstr>
      # Polynomial.to_string qq;;
      - : string =
      "(x^3 - 111*x^2 + 6064*x - 189804)*y^3 + (x^3 - 111*x^2 + 6064*x - 189804)*y^2"
      # Polynomial.is_irreducible q;;
      - : bool = true
      # let qmin = Polynomial.minimal q;;
      val qmin : Rational.t Polynomial.t = <abstr>
      # Polynomial.to_string qmin;;
      - : string = "x^3 - x^2 - 60*x - 364"
      # Number_field.(are_isomorphic (create q) (create qmin));;
      - : bool = true
      ]} *)

  val sign : t -> Signed.Long.t * Signed.Long.t
  val discriminant : t -> Integer.t
  val z_basis : t -> (elt, [ `ROW ]) Vector.t
  val elt : Rational.t array -> elt
  val add : t -> elt -> elt -> elt
  val mul : t -> elt -> elt -> elt
  val equal : elt -> elt -> bool

  val divrem : t -> elt -> elt -> elt * elt
  (** [divrem nf a b] returns the pair (divisor, remainder) from
      the euclidean division of a by b.

      {@ocaml[
      # let gaussian_integers =
        (* Gaussian integers: the ring Z[i] (here we work in the field Q(i)) *)
        Number_field.create
        (Polynomial.create [| Rational.of_int 1; Rational.of_int 0; Rational.of_int 1 |]);;
      val gaussian_integers : Number_field.t = <abstr>
      # let a = Number_field.elt [| Rational.of_int 6; Rational.of_int 8 |];;
      val a : Number_field.elt = <abstr>
      # let b = Number_field.elt [| Rational.of_int 1; Rational.of_int 5 |];;
      val b : Number_field.elt = <abstr>
      # let q, r = (* Euclidean division of 6 + 8i by 1 + 5i. *)
        Number_field.divrem gaussian_integers a b;;
      val q : Number_field.elt = <abstr>
      val r : Number_field.elt = <abstr>
      # Number_field.(equal a (add gaussian_integers (mul gaussian_integers b q) r));;
      - : bool = true
      ]} *)

  val ideal_norm : t -> elt -> Integer.t

  val splitting :
    [< `F of t | `P of Integer.t Polynomial.t ] -> Integer.t Polynomial.t
  (** [splitting (nf|p)] given the number field [nf = Q(x)/(p)]
      or polynomial [p], returns the polynomial over Q for the
      splitting field of [p], that is the smallest field over
      which [p] is totally split. *)

  val nf_get_allroots : 'a ty -> 'a ty
  val nf_get_prec : 'a ty -> Signed.long

  val nfmaxord_to_nf :
    nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
    'a ty ->
    Signed.long ->
    'a ty

  val nfcertify : 'a ty -> 'a ty
  val nfgaloismatrix : 'a ty -> 'a ty -> 'a ty
  val nfgaloismatrixapply : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfgaloispermtobasis : 'a ty -> 'a ty -> 'a ty

  val nfinit_basic :
    nfmaxord_t Ctypes.structure Ctypes_static.ptr -> 'a ty -> unit

  val nfinit_complete :
    nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
    Signed.long ->
    Signed.long ->
    'a ty

  val nfinit0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val nfinitred : 'a ty -> Signed.long -> 'a ty
  val nfinitred2 : 'a ty -> Signed.long -> 'a ty
  val nfisincl0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfisisom : 'a ty -> 'a ty -> 'a ty
  val nfnewprec : 'a ty -> Signed.long -> 'a ty
  val nfnewprec_shallow : 'a ty -> Signed.long -> 'a ty
  val nfpoleval : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfsplitting0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nftyp : 'a ty -> Signed.long
  val nfgrunwaldwang : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfgwkummer : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
  val mf_get_m : 'a ty -> 'a ty
  val mf_get_mindex : 'a ty -> 'a ty
  val mf_get_minv : 'a ty -> 'a ty
  val mf_get_basis : 'a ty -> 'a ty
  val mf_get_dim : 'a ty -> Signed.long
  val mf_get_e : 'a ty -> 'a ty
  val mf_get_fields : 'a ty -> 'a ty
  val mf_get_newforms : 'a ty -> 'a ty
  val mf_get_space : 'a ty -> Signed.long
  val mf_get_s : 'a ty -> 'a ty
  val mfcusp_get_vmjd : 'a ty -> 'a ty
  val mfnew_get_vj : 'a ty -> 'a ty

  val nf_to_fq_init :
    'a ty ->
    'a ty Ctypes_static.ptr ->
    'a ty Ctypes_static.ptr ->
    'a ty Ctypes_static.ptr ->
    'a ty

  val nf_to_fq : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfm_to_fqm : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfv_to_fqv : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfx_to_fqx : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfx_to_monic : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
  val nfbasis : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
  val nfcompositum : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfdiscfactors : 'a ty -> 'a ty

  val nfmaxord :
    nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
    'a ty ->
    Signed.long ->
    unit

  val nfmodpr : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfmodprinit : 'a ty -> 'a ty -> 'a ty
  val nfmodprinit0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfmodprlift : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfreducemodpr : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfdisc : 'a ty -> 'a ty
  val nf_to_scalar_or_alg : 'a ty -> 'a ty -> 'a ty
  val nf_to_scalar_or_basis : 'a ty -> 'a ty -> 'a ty
  val nf_cxlog : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfv_cxlog : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfchecksigns : 'a ty -> 'a ty -> 'a ty -> int
  val nfdiv : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfdiveuc : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfembed : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfeltembed : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

  val nfeltembed_i :
    'a ty Ctypes_static.ptr -> 'a ty -> 'a ty -> Signed.long -> 'a ty

  val nfeltsign : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nffactorback : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfinv : 'a ty -> 'a ty -> 'a ty
  val nfinvmodideal : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfissquare : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

  val nfispower :
    'a ty -> 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr -> Signed.long

  val nflogembed :
    'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long -> 'a ty

  val nfm_det : 'a ty -> 'a ty -> 'a ty
  val nfm_inv : 'a ty -> 'a ty -> 'a ty
  val nfm_ker : 'a ty -> 'a ty -> 'a ty
  val nfm_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfm_nfc_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfmod : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfmuli : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfnorm : 'a ty -> 'a ty -> 'a ty
  val nfpolsturm : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfpow : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfpow_u : 'a ty -> 'a ty -> pari_ulong -> 'a ty
  val nfpowmodideal : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfsign : 'a ty -> 'a ty -> 'a ty
  val nfsign_arch : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfsign_from_logarch : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfsqr : 'a ty -> 'a ty -> 'a ty
  val nfsqri : 'a ty -> 'a ty -> 'a ty
  val nfsub : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nftrace : 'a ty -> 'a ty -> 'a ty
  val nfval : 'a ty -> 'a ty -> 'a ty -> Signed.long

  val nfvalrem :
    'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

  val nf_get_varn : 'a ty -> Signed.long
  val nf_get_pol : 'a ty -> 'a ty
  val nf_get_degree : 'a ty -> Signed.long
  val nf_get_r1 : 'a ty -> Signed.long
  val nf_get_r2 : 'a ty -> Signed.long
  val nf_get_index : 'a ty -> 'a ty
  val nf_get_m : 'a ty -> 'a ty
  val nf_get_g : 'a ty -> 'a ty
  val nf_get_roundg : 'a ty -> 'a ty
  val nf_get_tr : 'a ty -> 'a ty
  val nf_get_diff : 'a ty -> 'a ty
  val nf_get_ramified_primes : 'a ty -> 'a ty
  val nf_get_roots : 'a ty -> 'a ty
  val nf_get_zkprimpart : 'a ty -> 'a ty
  val nf_get_zkden : 'a ty -> 'a ty
  val nf_get_invzk : 'a ty -> 'a ty
  val nfv_to_scalar_or_alg : 'a ty -> 'a ty -> 'a ty
  val nfc_multable_mul : 'a ty -> 'a ty -> 'a ty
  val nfc_nf_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nf_get_gtwist : 'a ty -> 'a ty -> 'a ty
  val nf_get_gtwist1 : 'a ty -> Signed.long -> 'a ty
  val nf_to_fp_coprime : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfdetint : 'a ty -> 'a ty -> 'a ty
  val nfdivmodpr : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfhnf : 'a ty -> 'a ty -> 'a ty
  val nfhnf0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfhnfmod : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfkermodpr : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfmulmodpr : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfpowmodpr : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfreduce : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfsnf : 'a ty -> 'a ty -> 'a ty
  val nfsnf0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfsolvemodpr : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfsubfields : 'a ty -> Signed.long -> 'a ty
  val nfsubfields0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val nfsubfieldscm : 'a ty -> Signed.long -> 'a ty
  val nfsubfieldsmax : 'a ty -> Signed.long -> 'a ty
  val nflist : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
  val nfresolvent : 'a ty -> Signed.long -> 'a ty
  val nf_pv_to_prv : 'a ty -> 'a ty -> 'a ty
  val nf_rnfeq : 'a ty -> 'a ty -> 'a ty
  val nf_rnfeqsimple : 'a ty -> 'a ty -> 'a ty
  val nf_nfzk : 'a ty -> 'a ty -> 'a ty
  val nfeltup : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfislocalpower : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long
  val nf_cxlog_normalize : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val nfcyclotomicunits : 'a ty -> 'a ty -> 'a ty
  val nfsign_units : 'a ty -> 'a ty -> int -> 'a ty
  val nfsign_tu : 'a ty -> 'a ty -> 'a ty
  val nfsign_fu : 'a ty -> 'a ty -> 'a ty
  val nf_deg1_prime : 'a ty -> 'a ty
  val nfarchstar : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nf_hyperell_locally_soluble : 'a ty -> 'a ty -> 'a ty -> Signed.long
  val nfhilbert : 'a ty -> 'a ty -> 'a ty -> Signed.long
  val nfhilbert0 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long

  val nfhyperellpadicfrobenius :
    'a ty -> 'a ty -> pari_ulong -> Signed.long -> 'a ty

  val nffactor : 'a ty -> 'a ty -> 'a ty
  val nffactormod : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val nfgcd : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

  val nfgcd_all :
    'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

  val nfissquarefree : 'a ty -> 'a ty -> int
  val nfroots : 'a ty -> 'a ty -> 'a ty
  val nfroots_if_split : 'a ty Ctypes_static.ptr -> 'a ty -> 'a ty
  val nfrootsof1 : 'a ty -> 'a ty
  val nfrootsq : 'a ty -> 'a ty
  val nfx_disc : 'a ty -> 'a ty -> 'a ty
  val nfx_resultant : 'a ty -> 'a ty -> 'a ty -> 'a ty

  module Infix : sig
    val ( = ) : elt -> elt -> bool
  end
end

module Elliptic_curve : sig
  type 'a t constraint 'a = 'b ty
  type 'a elt = 'a elliptic_curve ty constraint 'a = 'b ty

  val create :
    ?a1:'a ->
    ?a2:'a ->
    ?a3:'a ->
    ?a4:'a ->
    ?a6:'a ->
    ?dom:'a ->
    unit ->
    'a t option
  (** [create ?a1 ?a2 ?a3 ?a4 ?a6] defines the curve
        {math Y^2 + a_1 XY + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6}
        Returns [None] if the input coefficients do not define an elliptic curve
        over the field from the coefficients. *)

  val get_a1 : 'a t -> 'a
  val get_a2 : 'a t -> 'a
  val get_a3 : 'a t -> 'a
  val get_a4 : 'a t -> 'a
  val get_a6 : 'a t -> 'a
  val of_coordinates : x:'a -> y:'a -> 'a elt

  val order : 'a t -> Integer.t
  (** {@ocaml[
        # let g = Finite_field.generator ~order:(Integer.of_int 625) (* 5^4 *);;
        val g : Finite_field.t = <abstr>
        # let ell = Option.get (Elliptic_curve.create ~a4:(Finite_field.pow g (Integer.of_int 4)) ~a6:(Finite_field.pow g (Integer.of_int 6)) ());;
        val ell : Finite_field.t Elliptic_curve.t = <abstr>
        # Integer.(to_string (Elliptic_curve.order ell));;
        - : string = "675"
        ]} *)

  val discriminant : 'a t -> 'a
  val j_invariant : 'a t -> 'a
  val random : 'a t -> 'a elt

  val weil_pairing : 'a t -> l:Integer.t -> 'a elt -> 'a elt -> 'a
  (** [weil_pairing ell ~l p q] returns the Weil pairing of the two points
        of [l]-torsion [p] and [q] on the elliptic curve [ell].

        {@ocaml[
        # let l = Integer.of_int 3;;
        val l : Integer.t = <abstr>
        # let ord = Integer.of_int 103;;
        val ord : Integer.t = <abstr>
        # let ell = Option.get (Elliptic_curve.create ~a3:(Finite_field.prime_field_element (Integer.of_int 1) ~p:ord) ());;
        val ell : Finite_field.t Elliptic_curve.t = <abstr>
        # let (p, q) = Elliptic_curve.(of_coordinates ~x:(Finite_field.prime_field_element (Integer.of_int 0) ~p:ord) ~y:(Finite_field.prime_field_element (Integer.of_int 0) ~p:ord), of_coordinates ~x:(Finite_field.prime_field_element (Integer.of_int 57) ~p:ord) ~y:(Finite_field.prime_field_element (Integer.of_int 46) ~p:ord));;
        val p : Finite_field.t Elliptic_curve.elt = <abstr>
        val q : Finite_field.t Elliptic_curve.elt = <abstr>
        # let scalar = Elliptic_curve.weil_pairing ell ~l p q;;
        val scalar : Finite_field.t = <abstr>
        # Finite_field.to_string scalar (* 56 mod 103 *);;
        - : string = "56"
        # Finite_field.(to_string (pow scalar l));;
        - : string = "1"
        ]}
        *)

  val l_division_polynomial : 'a t -> l:Signed.Long.t -> 'a ty Polynomial.t
  (** {@ocaml[
        # let g = Finite_field.generator ~order:(Integer.of_int 625) (* 5^4 *);;
        val g : Finite_field.t = <abstr>
        # let ell = Option.get (Elliptic_curve.create ~a6:(Finite_field.pow g (Integer.of_int 6)) ());;
        val ell : Finite_field.t Elliptic_curve.t = <abstr>
        # let pdiv7 = (Elliptic_curve.l_division_polynomial ell ~l:(Signed.Long.of_int 7));;
        val pdiv7 : Finite_field.t ty Polynomial.t = <abstr>
        # Polynomial.to_string pdiv7;;
        - : string =
        "2*x^24 + (2*x^3 + x^2 + x + 2)*x^21 + (x^3 + x^2 + 3*x + 2)*x^18 + (2*x^3 + x^2 + 3*x)*x^15 + (2*x^3 + 4*x^2 + 3)*x^12 + (4*x^3 + x^2 + 2*x)*x^9 + (3*x^3 + 3*x^2)*x^6 + (4*x^3 + 2*x^2 + 3)*x^3 + (3*x^3 + 3*x^2 + 3*x + 1)"
        # Polynomial.degree pdiv7 = ((7 * 7 - 1) / 2);;
        - : bool = true
        ]} *)

  val to_string : 'a t -> string
  val add : 'a t -> 'a elt -> 'a elt -> 'a elt
  val sub : 'a t -> 'a elt -> 'a elt -> 'a elt
  val mul : 'a t -> n:Integer.t -> p:'a elt -> 'a elt
  val equal : 'a elt -> 'a elt -> bool
  val generators : 'a t -> ('a elt, [ `ROW ]) Vector.t
  val zero : 'a t -> 'a elt
  val get_coordinates : 'a t -> [> `inf | `point of 'a * 'a ]
  val order_elt : 'a t -> 'a elt -> Integer.t
  val to_string_elt : 'a elt -> string
  val log : 'a t -> base:'a elt -> 'a elt -> Integer.t option

  val ellanalyticrank : Rational.t t -> (Rational.t, [ `ROW ]) Vector.t
  (** {@ocaml[
        # let ell = Option.get (Elliptic_curve.create ~a6:(Rational.of_int 6) ());;
        val ell : Rational.t Elliptic_curve.t = <abstr>
        # Elliptic_curve.ellanalyticrank ell |> Vector.to_string
        - : string = "[0, 3.1205690727970642238215206887060527526]"
        ]} *)

  val ellanalyticrank_bitprec : 'a ty -> 'a ty -> Signed.long -> 'a ty

  val ellanal_globalred_all :
    'a ty ->
    'a ty Ctypes_static.ptr ->
    'a ty Ctypes_static.ptr ->
    'a ty Ctypes_static.ptr ->
    'a ty

  val ellsupersingularj_fpxq : 'a ty -> 'a ty -> 'a ty
  val elltrace_extension : 'a ty -> Signed.long -> 'a ty -> 'a ty
  val ellheegner : 'a ty -> 'a ty
  val elll1 : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val elll1_bitprec : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellconvertname : 'a ty -> 'a ty
  val elldatagenerators : 'a ty -> 'a ty
  val ellidentify : 'a ty -> 'a ty
  val ellsearch : 'a ty -> 'a ty
  val ellsearchcurve : 'a ty -> 'a ty
  val ell_is_integral : 'a ty -> int
  val ellq_get_cm : 'a ty -> Signed.long
  val ellq_get_n : 'a ty -> 'a ty

  val ellq_get_nfa :
    'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

  val ellqp_tate_uniformization : 'a ty -> Signed.long -> 'a ty
  val ellqp_agm : 'a ty -> Signed.long -> 'a ty
  val ellqp_u : 'a ty -> Signed.long -> 'a ty
  val ellqp_u2 : 'a ty -> Signed.long -> 'a ty
  val ellqp_q : 'a ty -> Signed.long -> 'a ty
  val ellqp_ab : 'a ty -> Signed.long -> 'a ty
  val ellqp_l : 'a ty -> Signed.long -> 'a ty
  val ellqp_root : 'a ty -> Signed.long -> 'a ty
  val ellqtwist_bsdperiod : 'a ty -> Signed.long -> 'a ty
  val ellr_area : 'a ty -> Signed.long -> 'a ty
  val ellr_ab : 'a ty -> Signed.long -> 'a ty
  val ellr_eta : 'a ty -> Signed.long -> 'a ty
  val ellr_omega : 'a ty -> Signed.long -> 'a ty
  val ellr_roots : 'a ty -> Signed.long -> 'a ty
  val ellan : 'a ty -> Signed.long -> 'a ty
  val ellanq_zv : 'a ty -> Signed.long -> 'a ty
  val ellanal_globalred : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
  val ellap : 'a ty -> 'a ty -> 'a ty
  val ellap_cm_fast : 'a ty -> pari_ulong -> Signed.long -> Signed.long
  val ellbasechar : 'a ty -> 'a ty
  val ellbsd : 'a ty -> Signed.long -> 'a ty
  val ellchangecurve : 'a ty -> 'a ty -> 'a ty
  val ellchangeinvert : 'a ty -> 'a ty
  val ellchangepoint : 'a ty -> 'a ty -> 'a ty
  val ellchangepointinv : 'a ty -> 'a ty -> 'a ty
  val elleisnum : 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty
  val elleta : 'a ty -> Signed.long -> 'a ty
  val elleulerf : 'a ty -> 'a ty -> 'a ty
  val ellff_get_card : 'a ty -> 'a ty
  val ellff_get_gens : 'a ty -> 'a ty
  val ellff_get_group : 'a ty -> 'a ty
  val ellff_get_o : 'a ty -> 'a ty
  val ellff_get_p : 'a ty -> 'a ty
  val ellff_get_m : 'a ty -> 'a ty
  val ellff_get_d : 'a ty -> 'a ty
  val ellfromj : 'a ty -> 'a ty
  val ellgenerators : 'a ty -> 'a ty
  val ellglobalred : 'a ty -> 'a ty
  val ellgroup : 'a ty -> 'a ty -> 'a ty
  val ellgroup0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val ellheight0 : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
  val ellheight : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val ellheightmatrix : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val ellheightoo : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val ellintegralmodel : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
  val ellintegralmodel_i : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
  val elliscm : 'a ty -> Signed.long
  val ellisoncurve : 'a ty -> 'a ty -> 'a ty
  val ellisotree : 'a ty -> 'a ty
  val ellissupersingular : 'a ty -> 'a ty -> int
  val elljissupersingular : 'a ty -> int
  val elllseries : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
  val elllocalred : 'a ty -> 'a ty -> 'a ty
  val elllog : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
  val ellminimaldisc : 'a ty -> 'a ty
  val ellminimalmodel : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
  val ellminimaltwist : 'a ty -> 'a ty
  val ellminimaltwist0 : 'a ty -> Signed.long -> 'a ty
  val ellminimaltwistcond : 'a ty -> 'a ty
  val ellnf_vecarea : 'a ty -> Signed.long -> 'a ty
  val ellnf_veceta : 'a ty -> Signed.long -> 'a ty
  val ellnf_vecomega : 'a ty -> Signed.long -> 'a ty
  val ellneg : 'a ty -> 'a ty -> 'a ty
  val ellorder : 'a ty -> 'a ty -> 'a ty -> 'a ty
  val ellorder_q : 'a ty -> 'a ty -> Signed.long
  val ellordinate : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val ellpadicheight0 : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
  val ellpadicheightmatrix : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
  val ellperiods : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellrootno : 'a ty -> 'a ty -> Signed.long
  val ellrootno_global : 'a ty -> Signed.long
  val ellsaturation : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellsea : 'a ty -> Signed.long -> 'a ty
  val ellsigma : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellsupersingularj : 'a ty -> 'a ty
  val elltamagawa : 'a ty -> 'a ty
  val elltaniyama : 'a ty -> Signed.long -> 'a ty
  val elltatepairing : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
  val elltors : 'a ty -> 'a ty
  val elltors0 : 'a ty -> Signed.long -> 'a ty
  val elltors_psylow : 'a ty -> pari_ulong -> 'a ty
  val elltrace : 'a ty -> 'a ty -> 'a ty
  val elltwist : 'a ty -> 'a ty -> 'a ty
  val ellwp : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val ellwp0 : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellwpseries : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellxn : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellzeta : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val ellfromeqn : 'a ty -> 'a ty
  val ellformaldifferential : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellformalexp : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellformallog : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellformalpoint : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellformalw : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellnonsingularmultiple : 'a ty -> 'a ty -> 'a ty

  val ellpadicl :
    'a ty -> 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty -> 'a ty

  val ellpadicbsd : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
  val ellpadicfrobenius : 'a ty -> pari_ulong -> Signed.long -> 'a ty
  val ellpadicheight : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
  val ellpadiclog : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
  val ellpadicregulator : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
  val ellpadics2 : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val ell2cover : 'a ty -> Signed.long -> 'a ty
  val ellrank : 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty
  val ellrankinit : 'a ty -> Signed.long -> 'a ty

  val ellisdivisible :
    'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

  val ellisogenyapply : 'a ty -> 'a ty -> 'a ty

  val ellisogeny :
    'a ty -> 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty

  val ellisomat : 'a ty -> Signed.long -> Signed.long -> 'a ty
  val ellweilcurve : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
  val ell_get_a1 : 'a ty -> 'a ty
  val ell_get_a2 : 'a ty -> 'a ty
  val ell_get_a3 : 'a ty -> 'a ty
  val ell_get_a4 : 'a ty -> 'a ty
  val ell_get_a6 : 'a ty -> 'a ty
  val ell_get_b2 : 'a ty -> 'a ty
  val ell_get_b4 : 'a ty -> 'a ty
  val ell_get_b6 : 'a ty -> 'a ty
  val ell_get_b8 : 'a ty -> 'a ty
  val ell_get_c4 : 'a ty -> 'a ty
  val ell_get_c6 : 'a ty -> 'a ty
  val ell_get_type : 'a ty -> Signed.long
  val ellff_get_field : 'a ty -> 'a ty
  val ellff_get_a4a6 : 'a ty -> 'a ty
  val ellqp_get_zero : 'a ty -> 'a ty
  val ellqp_get_prec : 'a ty -> Signed.long
  val ellqp_get_p : 'a ty -> 'a ty
  val ellr_get_prec : 'a ty -> Signed.long
  val ellr_get_sign : 'a ty -> Signed.long
  val ellnf_get_nf : 'a ty -> 'a ty
  val ellnf_get_bnf : 'a ty -> 'a ty
  val ellmodulareqn : Signed.long -> Signed.long -> Signed.long -> 'a ty
  val ellmoddegree : 'a ty -> 'a ty
  val ellratpoints : 'a ty -> 'a ty -> Signed.long -> 'a ty
  val elle : 'a ty -> Signed.long -> 'a ty
  val ellk : 'a ty -> Signed.long -> 'a ty

  val ellpadiclambdamu :
    'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty

  val ell_is_inf : 'a ty -> int
end

val logstyle_none : int64
val logstyle_plain : int64
val logstyle_color : int64
val logstyle_tex : int64
val pari_logstyles : pari_logstyles Ctypes.typ
val e_syntax : int64
val e_bug : int64
val e_alarm : int64
val e_file : int64
val e_misc : int64
val e_flag : int64
val e_impl : int64
val e_arch : int64
val e_package : int64
val e_notfunc : int64
val e_prec : int64
val e_type : int64
val e_dim : int64
val e_var : int64
val e_priority : int64
val e_user : int64
val e_stack : int64
val e_stackthread : int64
val e_overflow : int64
val e_domain : int64
val e_component : int64
val e_maxprime : int64
val e_constpol : int64
val e_irredpol : int64
val e_coprime : int64
val e_prime : int64
val e_modulus : int64
val e_roots0 : int64
val e_op : int64
val e_type2 : int64
val e_inv : int64
val e_mem : int64
val e_sqrtn : int64
val e_filedesc : int64
val e_none : int64
val err_list : err_list Ctypes.typ
val pari_timer : pari_timer Ctypes.structure Ctypes.typ
val pari_timer_s : (Signed.long, pari_timer Ctypes.structure) Ctypes.field
val pari_timer_us : (Signed.long, pari_timer Ctypes.structure) Ctypes.field
val pari_str : pari_str Ctypes.structure Ctypes.typ
val pari_str_string : (string, pari_str Ctypes.structure) Ctypes.field
val pari_str_end : (string, pari_str Ctypes.structure) Ctypes.field
val pari_str_cur : (string, pari_str Ctypes.structure) Ctypes.field
val pari_str_size : (int, pari_str Ctypes.structure) Ctypes.field
val pari_str_use_stack : (int, pari_str Ctypes.structure) Ctypes.field
val pari_sieve : pari_sieve Ctypes.structure Ctypes.typ
val pari_sieve_start : (pari_ulong, pari_sieve Ctypes.structure) Ctypes.field
val pari_sieve_end : (pari_ulong, pari_sieve Ctypes.structure) Ctypes.field
val pari_sieve_maxpos : (pari_ulong, pari_sieve Ctypes.structure) Ctypes.field
val pari_sieve_c : (pari_ulong, pari_sieve Ctypes.structure) Ctypes.field
val pari_sieve_q : (pari_ulong, pari_sieve Ctypes.structure) Ctypes.field

val pari_sieve_sieve :
  (Unsigned.uchar Ctypes_static.ptr, pari_sieve Ctypes.structure) Ctypes.field

val forprime_t : forprime_t Ctypes.structure Ctypes.typ
val forprime_t_strategy : (int, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_bb : ('a ty, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_c : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_q : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_d : (byteptr, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_p : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_b : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field

val forprime_t_psieve :
  ( pari_sieve Ctypes.structure Ctypes_static.ptr,
    forprime_t Ctypes.structure )
  Ctypes.field

val forprime_t_sieve :
  (Unsigned.uchar Ctypes_static.ptr, forprime_t Ctypes.structure) Ctypes.field

val forprime_t_isieve :
  (Unsigned.uchar Ctypes_static.ptr, forprime_t Ctypes.structure) Ctypes.field

val forprime_t_cache :
  (pari_ulong Ctypes_static.carray, forprime_t Ctypes.structure) Ctypes.field

val forprime_t_chunk : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_a : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_end : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_sieveb : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_pos : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_maxpos : (pari_ulong, forprime_t Ctypes.structure) Ctypes.field
val forprime_t_pp : ('a ty, forprime_t Ctypes.structure) Ctypes.field
val forcomposite_t : forcomposite_t Ctypes.structure Ctypes.typ
val forcomposite_t_first : (int, forcomposite_t Ctypes.structure) Ctypes.field
val forcomposite_t_b : ('a ty, forcomposite_t Ctypes.structure) Ctypes.field
val forcomposite_t_n : ('a ty, forcomposite_t Ctypes.structure) Ctypes.field
val forcomposite_t_p : ('a ty, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_T :
  (forprime_t Ctypes.structure, forcomposite_t Ctypes.structure) Ctypes.field

val forvec_t : forvec_t Ctypes.structure Ctypes.typ
val forvec_t_first : (Signed.long, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_a :
  ('a ty Ctypes_static.ptr, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_m :
  ('a ty Ctypes_static.ptr, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_M :
  ('a ty Ctypes_static.ptr, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_n : (Signed.long, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_next :
  ( (forvec_t Ctypes.structure Ctypes_static.ptr -> 'a ty)
    Ctypes_static.static_funptr,
    forvec_t Ctypes.structure )
  Ctypes.field

val forpart_t : forpart_t Ctypes.structure Ctypes.typ
val forpart_t_k : (Signed.long, forpart_t Ctypes.structure) Ctypes.field
val forpart_t_amax : (Signed.long, forpart_t Ctypes.structure) Ctypes.field
val forpart_t_amin : (Signed.long, forpart_t Ctypes.structure) Ctypes.field
val forpart_t_nmin : (Signed.long, forpart_t Ctypes.structure) Ctypes.field
val forpart_t_nmax : (Signed.long, forpart_t Ctypes.structure) Ctypes.field
val forpart_t_strip : (Signed.long, forpart_t Ctypes.structure) Ctypes.field
val forpart_t_v : ('a ty, forpart_t Ctypes.structure) Ctypes.field
val forperm_t : forperm_t Ctypes.structure Ctypes.typ
val forperm_t_k : (Signed.long, forperm_t Ctypes.structure) Ctypes.field
val forperm_t_first : (Signed.long, forperm_t Ctypes.structure) Ctypes.field
val forperm_t_v : ('a ty, forperm_t Ctypes.structure) Ctypes.field
val forsubset_t : forsubset_t Ctypes.structure Ctypes.typ
val forsubset_t_n : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_k : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_all : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_first : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_v : ('a ty, forsubset_t Ctypes.structure) Ctypes.field
val pari_plot : pari_plot Ctypes.structure Ctypes.typ

val pari_plot_draw :
  ( (pari_plot Ctypes.structure Ctypes_static.ptr ->
    'a ty ->
    'a ty ->
    'a ty ->
    unit)
    Ctypes_static.static_funptr,
    pari_plot Ctypes.structure )
  Ctypes.field

val pari_plot_width : (Signed.long, pari_plot Ctypes.structure) Ctypes.field
val pari_plot_height : (Signed.long, pari_plot Ctypes.structure) Ctypes.field
val pari_plot_hunit : (Signed.long, pari_plot Ctypes.structure) Ctypes.field
val pari_plot_vunit : (Signed.long, pari_plot Ctypes.structure) Ctypes.field
val pari_plot_fwidth : (Signed.long, pari_plot Ctypes.structure) Ctypes.field
val pari_plot_fheight : (Signed.long, pari_plot Ctypes.structure) Ctypes.field
val pari_plot_dwidth : (Signed.long, pari_plot Ctypes.structure) Ctypes.field
val pari_plot_dheight : (Signed.long, pari_plot Ctypes.structure) Ctypes.field
val genbin : genbin Ctypes.structure Ctypes.typ
val genbin_len : (int, genbin Ctypes.structure) Ctypes.field
val genbin_x : ('a ty, genbin Ctypes.structure) Ctypes.field
val genbin_base : ('a ty, genbin Ctypes.structure) Ctypes.field

val genbin_rebase :
  ( ('a ty -> Signed.long -> unit) Ctypes_static.static_funptr,
    genbin Ctypes.structure )
  Ctypes.field

val pari_mainstack : pari_mainstack Ctypes.structure Ctypes.typ
val pari_mainstack_top : (pari_sp, pari_mainstack Ctypes.structure) Ctypes.field
val pari_mainstack_bot : (pari_sp, pari_mainstack Ctypes.structure) Ctypes.field

val pari_mainstack_vbot :
  (pari_sp, pari_mainstack Ctypes.structure) Ctypes.field

val pari_mainstack_size : (int, pari_mainstack Ctypes.structure) Ctypes.field
val pari_mainstack_rsize : (int, pari_mainstack Ctypes.structure) Ctypes.field
val pari_mainstack_vsize : (int, pari_mainstack Ctypes.structure) Ctypes.field
val pari_mainstack_memused : (int, pari_mainstack Ctypes.structure) Ctypes.field
val entree : entree Ctypes.structure Ctypes.typ
val entree_name : (string, entree Ctypes.structure) Ctypes.field
val entree_valence : (pari_ulong, entree Ctypes.structure) Ctypes.field

val entree_value :
  (unit Ctypes_static.ptr, entree Ctypes.structure) Ctypes.field

val entree_menu : (Signed.long, entree Ctypes.structure) Ctypes.field
val entree_code : (string, entree Ctypes.structure) Ctypes.field
val entree_help : (string, entree Ctypes.structure) Ctypes.field

val entree_pvalue :
  (unit Ctypes_static.ptr, entree Ctypes.structure) Ctypes.field

val entree_arity : (Signed.long, entree Ctypes.structure) Ctypes.field
val entree_hash : (pari_ulong, entree Ctypes.structure) Ctypes.field

val entree_next :
  ( entree Ctypes.structure Ctypes_static.ptr,
    entree Ctypes.structure )
  Ctypes.field

val pari_parsestate : pari_parsestate Ctypes.structure Ctypes.typ

val pari_parsestate_node :
  (Signed.long, pari_parsestate Ctypes.structure) Ctypes.field

val pari_parsestate_once : (int, pari_parsestate Ctypes.structure) Ctypes.field

val pari_parsestate_discarded :
  (Signed.long, pari_parsestate Ctypes.structure) Ctypes.field

val pari_parsestate_lex_start :
  (string, pari_parsestate Ctypes.structure) Ctypes.field

val pari_parsestate_lasterror :
  ('a ty, pari_parsestate Ctypes.structure) Ctypes.field

val pari_compilestate : pari_compilestate Ctypes.structure Ctypes.typ

val pari_compilestate_opcode :
  (Signed.long, pari_compilestate Ctypes.structure) Ctypes.field

val pari_compilestate_operand :
  (Signed.long, pari_compilestate Ctypes.structure) Ctypes.field

val pari_compilestate_accesslex :
  (Signed.long, pari_compilestate Ctypes.structure) Ctypes.field

val pari_compilestate_data :
  (Signed.long, pari_compilestate Ctypes.structure) Ctypes.field

val pari_compilestate_localvars :
  (Signed.long, pari_compilestate Ctypes.structure) Ctypes.field

val pari_compilestate_frames :
  (Signed.long, pari_compilestate Ctypes.structure) Ctypes.field

val pari_compilestate_dbginfo :
  (Signed.long, pari_compilestate Ctypes.structure) Ctypes.field

val pari_compilestate_offset :
  (Signed.long, pari_compilestate Ctypes.structure) Ctypes.field

val pari_compilestate_nblex :
  (Signed.long, pari_compilestate Ctypes.structure) Ctypes.field

val pari_compilestate_dbgstart :
  (string, pari_compilestate Ctypes.structure) Ctypes.field

val pari_mtstate : pari_mtstate Ctypes.structure Ctypes.typ

val pari_mtstate_pending_threads :
  (Signed.long, pari_mtstate Ctypes.structure) Ctypes.field

val pari_mtstate_is_thread :
  (Signed.long, pari_mtstate Ctypes.structure) Ctypes.field

val pari_mtstate_trace_level :
  (Signed.long, pari_mtstate Ctypes.structure) Ctypes.field

val pari_evalstate : pari_evalstate Ctypes.structure Ctypes.typ

val pari_evalstate_avma :
  (pari_sp, pari_evalstate Ctypes.structure) Ctypes.field

val pari_evalstate_sp :
  (Signed.long, pari_evalstate Ctypes.structure) Ctypes.field

val pari_evalstate_rp :
  (Signed.long, pari_evalstate Ctypes.structure) Ctypes.field

val pari_evalstate_var :
  (Signed.long, pari_evalstate Ctypes.structure) Ctypes.field

val pari_evalstate_lvars :
  (Signed.long, pari_evalstate Ctypes.structure) Ctypes.field

val pari_evalstate_locks :
  (Signed.long, pari_evalstate Ctypes.structure) Ctypes.field

val pari_evalstate_prec :
  (Signed.long, pari_evalstate Ctypes.structure) Ctypes.field

val pari_evalstate_trace :
  (Signed.long, pari_evalstate Ctypes.structure) Ctypes.field

val pari_evalstate_mt :
  (pari_mtstate Ctypes.structure, pari_evalstate Ctypes.structure) Ctypes.field

val pari_evalstate_comp :
  ( pari_compilestate Ctypes.structure,
    pari_evalstate Ctypes.structure )
  Ctypes.field

val pari_varstate : pari_varstate Ctypes.structure Ctypes.typ

val pari_varstate_nvar :
  (Signed.long, pari_varstate Ctypes.structure) Ctypes.field

val pari_varstate_max_avail :
  (Signed.long, pari_varstate Ctypes.structure) Ctypes.field

val pari_varstate_min_priority :
  (Signed.long, pari_varstate Ctypes.structure) Ctypes.field

val pari_varstate_max_priority :
  (Signed.long, pari_varstate Ctypes.structure) Ctypes.field

val pari_global_state : pari_global_state Ctypes.structure Ctypes.typ

val pari_global_state_bitprec :
  (Signed.long, pari_global_state Ctypes.structure) Ctypes.field

val pari_global_state_primetab :
  ('a ty, pari_global_state Ctypes.structure) Ctypes.field

val pari_global_state_seadata :
  ('a ty, pari_global_state Ctypes.structure) Ctypes.field

val pari_global_state_varpriority :
  ( Signed.long Ctypes_static.ptr,
    pari_global_state Ctypes.structure )
  Ctypes.field

val pari_global_state_varstate :
  ( pari_varstate Ctypes.structure,
    pari_global_state Ctypes.structure )
  Ctypes.field

val pari_thread : pari_thread Ctypes.structure Ctypes.typ

val pari_thread_st :
  (pari_mainstack Ctypes.structure, pari_thread Ctypes.structure) Ctypes.field

val pari_thread_gs :
  ( pari_global_state Ctypes.structure,
    pari_thread Ctypes.structure )
  Ctypes.field

val pari_thread_data : ('a ty, pari_thread Ctypes.structure) Ctypes.field
val mt_state : mt_state Ctypes.structure Ctypes.typ
val mt_state_worker : ('a ty, mt_state Ctypes.structure) Ctypes.field
val mt_state_pending : ('a ty, mt_state Ctypes.structure) Ctypes.field
val mt_state_workid : (Signed.long, mt_state Ctypes.structure) Ctypes.field
val pari_mt : pari_mt Ctypes.structure Ctypes.typ

val pari_mt_mt :
  (mt_state Ctypes.structure, pari_mt Ctypes.structure) Ctypes.field

val pari_mt_get :
  ( (mt_state Ctypes.structure Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    'a ty)
    Ctypes_static.static_funptr,
    pari_mt Ctypes.structure )
  Ctypes.field

val pari_mt_submit :
  ( (mt_state Ctypes.structure Ctypes_static.ptr -> Signed.long -> 'a ty -> unit)
    Ctypes_static.static_funptr,
    pari_mt Ctypes.structure )
  Ctypes.field

val pari_mt_end :
  (unit Ctypes_static.static_funptr, pari_mt Ctypes.structure) Ctypes.field

val parfor_iter : parfor_iter Ctypes.structure Ctypes.typ

val parfor_iter_pending :
  (Signed.long, parfor_iter Ctypes.structure) Ctypes.field

val parfor_iter_worker : ('a ty, parfor_iter Ctypes.structure) Ctypes.field

val parfor_iter_pt :
  (pari_mt Ctypes.structure, parfor_iter Ctypes.structure) Ctypes.field

val parfor_t : parfor_t Ctypes.structure Ctypes.typ
val parfor_t_a : ('a ty, parfor_t Ctypes.structure) Ctypes.field
val parfor_t_b : ('a ty, parfor_t Ctypes.structure) Ctypes.field

val parfor_t_iter :
  (parfor_iter Ctypes.structure, parfor_t Ctypes.structure) Ctypes.field

val parforeach_t : parforeach_t Ctypes.structure Ctypes.typ
val parforeach_t_x : ('a ty, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_W : ('a ty, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_i : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_l : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field

val parforeach_t_iter :
  (parfor_iter Ctypes.structure, parforeach_t Ctypes.structure) Ctypes.field

val parforprime_t : parforprime_t Ctypes.structure Ctypes.typ
val parforprime_t_v : ('a ty, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_forprime :
  (forprime_t Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_iter :
  (parfor_iter Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforvec_t : parforvec_t Ctypes.structure Ctypes.typ
val parforvec_t_v : ('a ty, parforvec_t Ctypes.structure) Ctypes.field

val parforvec_t_forvec :
  (forvec_t Ctypes.structure, parforvec_t Ctypes.structure) Ctypes.field

val parforvec_t_iter :
  (parfor_iter Ctypes.structure, parforvec_t Ctypes.structure) Ctypes.field

val hashentry : hashentry Ctypes.structure Ctypes.typ

val hashentry_key :
  (unit Ctypes_static.ptr, hashentry Ctypes.structure) Ctypes.field

val hashentry_val :
  (unit Ctypes_static.ptr, hashentry Ctypes.structure) Ctypes.field

val hashentry_hash : (pari_ulong, hashentry Ctypes.structure) Ctypes.field

val hashentry_next :
  ( hashentry Ctypes.structure Ctypes_static.ptr,
    hashentry Ctypes.structure )
  Ctypes.field

val hashtable : hashtable Ctypes.structure Ctypes.typ
val hashtable_len : (pari_ulong, hashtable Ctypes.structure) Ctypes.field

val hashtable_table :
  ( hashentry Ctypes.structure Ctypes_static.ptr Ctypes_static.ptr,
    hashtable Ctypes.structure )
  Ctypes.field

val hashtable_nb : (pari_ulong, hashtable Ctypes.structure) Ctypes.field
val hashtable_maxnb : (pari_ulong, hashtable Ctypes.structure) Ctypes.field
val hashtable_pindex : (pari_ulong, hashtable Ctypes.structure) Ctypes.field

val hashtable_hash :
  ( (unit Ctypes_static.ptr -> pari_ulong) Ctypes_static.static_funptr,
    hashtable Ctypes.structure )
  Ctypes.field

val hashtable_eq :
  ( (unit Ctypes_static.ptr -> unit Ctypes_static.ptr -> int)
    Ctypes_static.static_funptr,
    hashtable Ctypes.structure )
  Ctypes.field

val hashtable_use_stack : (int, hashtable Ctypes.structure) Ctypes.field
val gp_path : gp_path Ctypes.structure Ctypes.typ
val gp_path_PATH : (string, gp_path Ctypes.structure) Ctypes.field

val gp_path_dirs :
  (string Ctypes_static.ptr, gp_path Ctypes.structure) Ctypes.field

val pariout_t : pariout_t Ctypes.structure Ctypes.typ
val pariout_t_format : (char, pariout_t Ctypes.structure) Ctypes.field
val pariout_t_sigd : (Signed.long, pariout_t Ctypes.structure) Ctypes.field
val pariout_t_sp : (int, pariout_t Ctypes.structure) Ctypes.field
val pariout_t_prettyp : (int, pariout_t Ctypes.structure) Ctypes.field
val pariout_t_TeXstyle : (int, pariout_t Ctypes.structure) Ctypes.field
val nfmaxord_t : nfmaxord_t Ctypes.structure Ctypes.typ
val nfmaxord_t_T : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dT : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_T0 : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_unscale : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dK : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_index : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_basis : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_r1 : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_basden : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dTP : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dTE : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dKP : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dKE : ('a ty, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_certify : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field
val qfr_data : qfr_data Ctypes.structure Ctypes.typ
val qfr_data_D : ('a ty, qfr_data Ctypes.structure) Ctypes.field
val qfr_data_sqrtD : ('a ty, qfr_data Ctypes.structure) Ctypes.field
val qfr_data_isqrtD : ('a ty, qfr_data Ctypes.structure) Ctypes.field
val fp_chk_fun : fp_chk_fun Ctypes.structure Ctypes.typ

val fp_chk_fun_f :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_init :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_post :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_data :
  (unit Ctypes_static.ptr, fp_chk_fun Ctypes.structure) Ctypes.field

val fp_chk_fun_skipfirst :
  (Signed.long, fp_chk_fun Ctypes.structure) Ctypes.field

val zlog_s : zlog_s Ctypes.structure Ctypes.typ
val zlog_s_bid : ('a ty, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_P : ('a ty, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_k : ('a ty, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_sprk : ('a ty, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_archp : ('a ty, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_mod : ('a ty, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_U : ('a ty, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_hU : (Signed.long, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_no2 : (int, zlog_s Ctypes.structure) Ctypes.field
val bb_group : bb_group Ctypes.structure Ctypes.typ

val bb_group_mul :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_pow :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_rand :
  ( (unit Ctypes_static.ptr -> 'a ty) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_hash :
  ( ('a ty -> pari_ulong) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal :
  ( ('a ty -> 'a ty -> int) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal1 :
  ( ('a ty -> int) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_easylog :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_field : bb_field Ctypes.structure Ctypes.typ

val bb_field_red :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_add :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_mul :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_neg :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_inv :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_equal0 :
  ( ('a ty -> int) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_s :
  ( (unit Ctypes_static.ptr -> Signed.long -> 'a ty) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_algebra : bb_algebra Ctypes.structure Ctypes.typ

val bb_algebra_red :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_add :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sub :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_mul :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sqr :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_one :
  ( (unit Ctypes_static.ptr -> 'a ty) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_zero :
  ( (unit Ctypes_static.ptr -> 'a ty) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_ring : bb_ring Ctypes.structure Ctypes.typ

val bb_ring_add :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_mul :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
    Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_sqr :
  ( (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val buchimag : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val buchreal : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zidealstar : 'a ty -> 'a ty -> 'a ty
val zidealstarinit : 'a ty -> 'a ty -> 'a ty
val zidealstarinitgen : 'a ty -> 'a ty -> 'a ty
val factmod : 'a ty -> 'a ty -> 'a ty
val mpbern : Signed.long -> Signed.long -> unit
val simplefactmod : 'a ty -> 'a ty -> 'a ty
val listkill : 'a ty -> unit
val isprincipalforce : 'a ty -> 'a ty -> 'a ty
val isprincipalgen : 'a ty -> 'a ty -> 'a ty
val isprincipalgenforce : 'a ty -> 'a ty -> 'a ty
val f2ms_ker : 'a ty -> Signed.long -> 'a ty
val f2ms_to_f2m : 'a ty -> Signed.long -> 'a ty
val f2c_to_zc : 'a ty -> 'a ty
val f2c_to_mod : 'a ty -> 'a ty
val f2m_f2c_gauss : 'a ty -> 'a ty -> 'a ty
val f2m_f2c_invimage : 'a ty -> 'a ty -> 'a ty
val f2m_f2c_mul : 'a ty -> 'a ty -> 'a ty
val f2m_deplin : 'a ty -> 'a ty
val f2m_det : 'a ty -> pari_ulong
val f2m_det_sp : 'a ty -> pari_ulong
val f2m_gauss : 'a ty -> 'a ty -> 'a ty
val f2m_inv : 'a ty -> 'a ty
val f2m_invimage : 'a ty -> 'a ty -> 'a ty
val f2m_ker : 'a ty -> 'a ty
val f2m_ker_sp : 'a ty -> Signed.long -> 'a ty
val f2m_mul : 'a ty -> 'a ty -> 'a ty
val f2m_powu : 'a ty -> pari_ulong -> 'a ty
val f2m_rank : 'a ty -> Signed.long
val f2m_row : 'a ty -> Signed.long -> 'a ty
val f2m_rowslice : 'a ty -> Signed.long -> Signed.long -> 'a ty
val f2m_to_f2ms : 'a ty -> 'a ty
val f2m_to_flm : 'a ty -> 'a ty
val f2m_to_zm : 'a ty -> 'a ty
val f2m_to_mod : 'a ty -> 'a ty
val f2m_transpose : 'a ty -> 'a ty
val f2v_add_inplace : 'a ty -> 'a ty -> unit
val f2v_and_inplace : 'a ty -> 'a ty -> unit
val f2v_dotproduct : 'a ty -> 'a ty -> pari_ulong
val f2v_equal0 : 'a ty -> int
val f2v_hamming : 'a ty -> pari_ulong
val f2v_negimply_inplace : 'a ty -> 'a ty -> unit
val f2v_or_inplace : 'a ty -> 'a ty -> unit
val f2v_slice : 'a ty -> Signed.long -> Signed.long -> 'a ty
val f2v_subset : 'a ty -> 'a ty -> int
val f2v_to_flv : 'a ty -> 'a ty
val matid_f2m : Signed.long -> 'a ty
val f2x_f2xq_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2x_f2xqv_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2x_frobenius : 'a ty -> 'a ty
val f2x_1_add : 'a ty -> 'a ty
val f2x_add : 'a ty -> 'a ty -> 'a ty
val f2x_deflate : 'a ty -> Signed.long -> 'a ty
val f2x_degfact : 'a ty -> 'a ty
val f2x_degree : 'a ty -> Signed.long
val f2x_deriv : 'a ty -> 'a ty
val f2x_divrem : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val f2x_eval : 'a ty -> pari_ulong -> pari_ulong

val f2x_even_odd :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val f2x_extgcd :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val f2x_gcd : 'a ty -> 'a ty -> 'a ty
val f2x_get_red : 'a ty -> 'a ty
val f2x_halfgcd : 'a ty -> 'a ty -> 'a ty
val f2x_issquare : 'a ty -> int
val f2x_matfrobenius : 'a ty -> 'a ty
val f2x_mul : 'a ty -> 'a ty -> 'a ty
val f2x_recip : 'a ty -> 'a ty
val f2x_rem : 'a ty -> 'a ty -> 'a ty
val f2x_shift : 'a ty -> Signed.long -> 'a ty
val f2x_sqr : 'a ty -> 'a ty
val f2x_sqrt : 'a ty -> 'a ty
val f2x_to_f2v : 'a ty -> Signed.long -> 'a ty
val f2x_to_f2xx : 'a ty -> Signed.long -> 'a ty
val f2x_to_flx : 'a ty -> 'a ty
val f2x_to_zx : 'a ty -> 'a ty
val f2x_valrem : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val f2xc_to_flxc : 'a ty -> 'a ty
val f2xc_to_zxc : 'a ty -> 'a ty
val f2xv_to_f2m : 'a ty -> Signed.long -> 'a ty
val f2xv_to_flxv_inplace : 'a ty -> unit
val f2xv_to_zxv_inplace : 'a ty -> unit
val f2xx_f2x_add : 'a ty -> 'a ty -> 'a ty
val f2xx_f2x_mul : 'a ty -> 'a ty -> 'a ty
val f2xx_add : 'a ty -> 'a ty -> 'a ty
val f2xx_deriv : 'a ty -> 'a ty
val f2xx_renormalize : 'a ty -> Signed.long -> 'a ty
val f2xx_to_kronecker : 'a ty -> Signed.long -> 'a ty
val f2xx_to_flxx : 'a ty -> 'a ty
val f2xx_to_zxx : 'a ty -> 'a ty
val f2xx_to_f2xc : 'a ty -> Signed.long -> Signed.long -> 'a ty
val f2xxv_to_f2xm : 'a ty -> Signed.long -> Signed.long -> 'a ty
val f2xxc_to_zxxc : 'a ty -> 'a ty
val f2xy_f2xq_evalx : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xy_f2xqv_evalx : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xy_degreex : 'a ty -> Signed.long
val f2xn_div : 'a ty -> 'a ty -> Signed.long -> 'a ty
val f2xn_inv : 'a ty -> Signed.long -> 'a ty
val f2xn_red : 'a ty -> Signed.long -> 'a ty
val f2xq_artin_schreier : 'a ty -> 'a ty -> 'a ty
val f2xq_autpow : 'a ty -> Signed.long -> 'a ty -> 'a ty
val f2xq_conjvec : 'a ty -> 'a ty -> 'a ty
val f2xq_div : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xq_inv : 'a ty -> 'a ty -> 'a ty
val f2xq_invsafe : 'a ty -> 'a ty -> 'a ty
val f2xq_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xq_matrix_pow : 'a ty -> Signed.long -> Signed.long -> 'a ty -> 'a ty
val f2xq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xq_order : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xq_pow_init : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val f2xq_pow_table : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xq_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val f2xq_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty
val f2xq_sqr : 'a ty -> 'a ty -> 'a ty
val f2xq_sqrt : 'a ty -> 'a ty -> 'a ty
val f2xq_sqrt_fast : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xq_sqrtn : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val f2xq_trace : 'a ty -> 'a ty -> pari_ulong
val f2xqx_f2xq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqx_f2xq_mul_to_monic : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqx_f2xqxq_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqx_f2xqxqv_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqx_disc : 'a ty -> 'a ty -> 'a ty
val f2xqx_divrem : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val f2xqx_extgcd :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val f2xqx_gcd : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqx_get_red : 'a ty -> 'a ty -> 'a ty
val f2xqx_halfgcd : 'a ty -> 'a ty -> 'a ty -> 'a ty

val f2xqx_halfgcd_all :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val f2xqx_invbarrett : 'a ty -> 'a ty -> 'a ty

val f2xqx_ispower :
  'a ty -> Signed.long -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

val f2xqx_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqx_normalize : 'a ty -> 'a ty -> 'a ty
val f2xqx_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val f2xqx_red : 'a ty -> 'a ty -> 'a ty
val f2xqx_rem : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqx_resultant : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqx_sqr : 'a ty -> 'a ty -> 'a ty
val f2xqxq_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqxq_invsafe : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqxq_sqr : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqxq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqxq_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val f2xqxq_autpow : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val f2xqxq_auttrace : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val f2xqxqv_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val flm_to_f2m : 'a ty -> 'a ty
val flv_to_f2v : 'a ty -> 'a ty
val flx_to_f2x : 'a ty -> 'a ty
val flxc_to_f2xc : 'a ty -> 'a ty
val flxx_to_f2xx : 'a ty -> 'a ty
val flxxc_to_f2xxc : 'a ty -> 'a ty
val kronecker_to_f2xqx : 'a ty -> 'a ty -> 'a ty
val rg_to_f2xq : 'a ty -> 'a ty -> 'a ty
val rgm_to_f2m : 'a ty -> 'a ty
val rgv_to_f2v : 'a ty -> 'a ty
val rgx_to_f2x : 'a ty -> 'a ty
val z_to_f2x : 'a ty -> Signed.long -> 'a ty
val zm_to_f2m : 'a ty -> 'a ty
val zv_to_f2v : 'a ty -> 'a ty
val zx_to_f2x : 'a ty -> 'a ty
val zxx_to_f2xx : 'a ty -> Signed.long -> 'a ty
val const_f2v : Signed.long -> 'a ty
val gener_f2xq : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val get_f2xq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  bb_field Ctypes.structure Ctypes_static.ptr

val monomial_f2x : Signed.long -> Signed.long -> 'a ty
val random_f2xqx : Signed.long -> Signed.long -> 'a ty -> 'a ty
val f2x_teichmuller : 'a ty -> Signed.long -> 'a ty
val f2xq_ellcard : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xq_ellgens : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val f2xq_ellgroup :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val f2xq_elltwist :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  unit

val f2xqe_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_changepoint : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_changepointinv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_dbl : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_neg : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_order : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_tatepairing : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqe_weilpairing : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val get_f2xqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_f2xqe : 'a ty -> 'a ty -> 'a ty
val random_f2xqe : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f3c_to_mod : 'a ty -> 'a ty
val f3c_to_zc : 'a ty -> 'a ty
val f3m_ker : 'a ty -> 'a ty
val f3m_ker_sp : 'a ty -> Signed.long -> 'a ty
val f3m_mul : 'a ty -> 'a ty -> 'a ty
val f3m_row : 'a ty -> Signed.long -> 'a ty
val f3m_to_flm : 'a ty -> 'a ty
val f3m_to_zm : 'a ty -> 'a ty
val f3m_to_mod : 'a ty -> 'a ty
val f3m_transpose : 'a ty -> 'a ty
val f3v_to_flv : 'a ty -> 'a ty
val f3v_coeff : 'a ty -> Signed.long -> pari_ulong
val f3v_clear : 'a ty -> Signed.long -> unit
val f3v_set : 'a ty -> Signed.long -> pari_ulong -> unit
val flm_to_f3m : 'a ty -> 'a ty
val flv_to_f3v : 'a ty -> 'a ty
val rgm_to_f3m : 'a ty -> 'a ty
val rgv_to_f3v : 'a ty -> 'a ty
val zm_to_f3m : 'a ty -> 'a ty
val zv_to_f3v : 'a ty -> 'a ty
val zero_f3m_copy : Signed.long -> Signed.long -> 'a ty
val zero_f3v : Signed.long -> 'a ty
val fl_elldisc : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_elldisc_pre :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_ellj : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_ellj_to_a4a6 :
  pari_ulong ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  unit

val fl_ellptors :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val fl_elltwist :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  unit

val fl_elltwist_disc :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  unit

val fle_add : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val fle_dbl : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val fle_changepoint : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val fle_changepointinv : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val fle_log : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val fle_mul : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val fle_mulu : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty
val fle_order : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val fle_sub : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val fle_tatepairing :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fle_to_flj : 'a ty -> 'a ty

val fle_weilpairing :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val flj_add_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val flj_changepointinv_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flj_dbl_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val flj_mulu_pre :
  'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val flj_neg : 'a ty -> pari_ulong -> 'a ty
val flj_to_fle : 'a ty -> pari_ulong -> 'a ty
val flj_to_fle_pre : 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val fljv_factorback_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val random_fle : pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val random_fle_pre :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val random_flj_pre :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val flc_to_zc : 'a ty -> 'a ty
val flc_to_zc_inplace : 'a ty -> 'a ty
val flm_flc_gauss : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_flc_invimage : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_adjoint : 'a ty -> pari_ulong -> 'a ty
val flm_deplin : 'a ty -> pari_ulong -> 'a ty
val flm_det : 'a ty -> pari_ulong -> pari_ulong
val flm_det_sp : 'a ty -> pari_ulong -> pari_ulong
val flm_gauss : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_intersect : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_intersect_i : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_inv : 'a ty -> pari_ulong -> 'a ty
val flm_invimage : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_ker : 'a ty -> pari_ulong -> 'a ty
val flm_ker_sp : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val flm_rank : 'a ty -> pari_ulong -> Signed.long
val flm_to_zm : 'a ty -> 'a ty
val flm_to_zm_inplace : 'a ty -> 'a ty
val flv_to_zv : 'a ty -> 'a ty
val fl_to_flx : pari_ulong -> Signed.long -> 'a ty
val fl2_equal1 : 'a ty -> int
val fl2_inv_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val fl2_mul_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val fl2_norm_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl2_pow_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val fl2_sqr_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty
val fl2_sqrt_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val fl2_sqrtn_pre :
  'a ty ->
  'a ty ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flm_to_flxv : 'a ty -> Signed.long -> 'a ty
val flm_to_flxx : 'a ty -> Signed.long -> Signed.long -> 'a ty
val flv_flm_polint : 'a ty -> 'a ty -> pari_ulong -> Signed.long -> 'a ty
val flv_inv : 'a ty -> pari_ulong -> 'a ty
val flv_inv_inplace : 'a ty -> pari_ulong -> unit
val flv_inv_pre_inplace : 'a ty -> pari_ulong -> pari_ulong -> unit
val flv_inv_pre : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flv_invvandermonde : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flv_polint : 'a ty -> 'a ty -> pari_ulong -> Signed.long -> 'a ty
val flv_prod : 'a ty -> pari_ulong -> pari_ulong
val flv_prod_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong
val flv_roots_to_pol : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val flv_to_flx : 'a ty -> Signed.long -> 'a ty
val flx_fl_add : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_fl_mul : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_fl_mul_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty
val flx_fl_mul_to_monic : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_fl_sub : 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flx_fl2_eval_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty

val flx_flv_multieval : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_flxq_eval : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flx_flxq_eval_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flx_flxqv_eval : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flx_flxqv_eval_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flx_frobenius : 'a ty -> pari_ulong -> 'a ty
val flx_frobenius_pre : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_laplace : 'a ty -> pari_ulong -> 'a ty
val flx_newton : 'a ty -> Signed.long -> pari_ulong -> 'a ty
val flx_add : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_blocks : 'a ty -> Signed.long -> Signed.long -> 'a ty
val flx_composedprod : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_composedsum : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_convol : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_deflate : 'a ty -> Signed.long -> 'a ty
val flx_deriv : 'a ty -> pari_ulong -> 'a ty
val flx_diff1 : 'a ty -> pari_ulong -> 'a ty
val flx_digits : 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flx_div_by_x_x :
  'a ty -> pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> 'a ty

val flx_divrem :
  'a ty -> 'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> 'a ty

val flx_divrem_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty Ctypes_static.ptr -> 'a ty

val flx_double : 'a ty -> pari_ulong -> 'a ty
val flx_eval : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong

val flx_eval_powers_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong

val flx_eval_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val flx_extgcd :
  'a ty ->
  'a ty ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flx_extgcd_pre :
  'a ty ->
  'a ty ->
  pari_ulong ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flx_extresultant :
  'a ty ->
  'a ty ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  pari_ulong

val flx_extresultant_pre :
  'a ty ->
  'a ty ->
  pari_ulong ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  pari_ulong

val flx_fromnewton : 'a ty -> pari_ulong -> 'a ty
val flx_gcd : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_gcd_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_get_red : 'a ty -> pari_ulong -> 'a ty
val flx_get_red_pre : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_halfgcd : 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flx_halfgcd_all :
  'a ty ->
  'a ty ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flx_halfgcd_all_pre :
  'a ty ->
  'a ty ->
  pari_ulong ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flx_halfgcd_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_halve : 'a ty -> pari_ulong -> 'a ty
val flx_inflate : 'a ty -> Signed.long -> 'a ty
val flx_integ : 'a ty -> pari_ulong -> 'a ty
val flx_invbarrett : 'a ty -> pari_ulong -> 'a ty
val flx_invlaplace : 'a ty -> pari_ulong -> 'a ty
val flx_is_squarefree : 'a ty -> pari_ulong -> int
val flx_is_smooth : 'a ty -> Signed.long -> pari_ulong -> int
val flx_is_smooth_pre : 'a ty -> Signed.long -> pari_ulong -> pari_ulong -> int
val flx_matfrobenius : 'a ty -> pari_ulong -> 'a ty
val flx_matfrobenius_pre : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_mod_xn1 : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_mod_xnm1 : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_mul : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_mul_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_neg : 'a ty -> pari_ulong -> 'a ty
val flx_neg_inplace : 'a ty -> pari_ulong -> 'a ty
val flx_normalize : 'a ty -> pari_ulong -> 'a ty
val flx_powu : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_powu_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty
val flx_recip : 'a ty -> 'a ty
val flx_red : 'a ty -> pari_ulong -> 'a ty
val flx_rem : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_rem_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_renormalize : 'a ty -> Signed.long -> 'a ty
val flx_rescale : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_resultant : 'a ty -> 'a ty -> pari_ulong -> pari_ulong
val flx_resultant_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong
val flx_shift : 'a ty -> Signed.long -> 'a ty
val flx_splitting : 'a ty -> Signed.long -> 'a ty
val flx_sqr : 'a ty -> pari_ulong -> 'a ty
val flx_sqr_pre : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_sub : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_translate1 : 'a ty -> pari_ulong -> 'a ty
val flx_translate1_basecase : 'a ty -> pari_ulong -> 'a ty
val flx_to_flv : 'a ty -> Signed.long -> 'a ty
val flx_to_flxx : 'a ty -> Signed.long -> 'a ty
val flx_to_zx : 'a ty -> 'a ty
val flx_to_zx_inplace : 'a ty -> 'a ty
val flx_triple : 'a ty -> pari_ulong -> 'a ty
val flx_val : 'a ty -> Signed.long
val flx_valrem : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val flxc_flxqv_eval : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxc_flxqv_eval_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxc_flxq_eval : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxc_flxq_eval_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxc_eval_powers_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxc_neg : 'a ty -> pari_ulong -> 'a ty
val flxc_sub : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxc_to_zxc : 'a ty -> 'a ty
val flxm_flx_add_shallow : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxm_eval_powers_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxm_neg : 'a ty -> pari_ulong -> 'a ty
val flxm_sub : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxm_to_flxxv : 'a ty -> Signed.long -> 'a ty
val flxm_to_zxm : 'a ty -> 'a ty
val flxt_red : 'a ty -> pari_ulong -> 'a ty
val flxv_flc_mul : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxv_flv_multieval : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxv_flx_fromdigits : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxv_composedsum : 'a ty -> pari_ulong -> 'a ty
val flxv_prod : 'a ty -> pari_ulong -> 'a ty
val flxv_red : 'a ty -> pari_ulong -> 'a ty
val flxv_to_flm : 'a ty -> Signed.long -> 'a ty
val flxv_to_flxx : 'a ty -> Signed.long -> 'a ty
val flxv_to_zxv : 'a ty -> 'a ty
val flxv_to_zxv_inplace : 'a ty -> unit
val flxn_div : 'a ty -> 'a ty -> Signed.long -> pari_ulong -> 'a ty

val flxn_div_pre :
  'a ty -> 'a ty -> Signed.long -> pari_ulong -> pari_ulong -> 'a ty

val flxn_exp : 'a ty -> Signed.long -> pari_ulong -> 'a ty
val flxn_expint : 'a ty -> Signed.long -> pari_ulong -> 'a ty
val flxn_inv : 'a ty -> Signed.long -> pari_ulong -> 'a ty
val flxn_mul : 'a ty -> 'a ty -> Signed.long -> pari_ulong -> 'a ty

val flxn_mul_pre :
  'a ty -> 'a ty -> Signed.long -> pari_ulong -> pari_ulong -> 'a ty

val flxn_sqr : 'a ty -> Signed.long -> pari_ulong -> 'a ty
val flxn_sqr_pre : 'a ty -> Signed.long -> pari_ulong -> pari_ulong -> 'a ty
val flxn_red : 'a ty -> Signed.long -> 'a ty
val flxq_autpow : 'a ty -> pari_ulong -> 'a ty -> pari_ulong -> 'a ty

val flxq_autpow_pre :
  'a ty -> pari_ulong -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxq_autpowers : 'a ty -> pari_ulong -> 'a ty -> pari_ulong -> 'a ty
val flxq_autsum : 'a ty -> pari_ulong -> 'a ty -> pari_ulong -> 'a ty
val flxq_auttrace : 'a ty -> pari_ulong -> 'a ty -> pari_ulong -> 'a ty

val flxq_auttrace_pre :
  'a ty -> pari_ulong -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxq_charpoly : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_conjvec : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_div : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_div_pre : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxq_inv : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_inv_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxq_invsafe : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_invsafe_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxq_issquare : 'a ty -> 'a ty -> pari_ulong -> int
val flxq_is2npower : 'a ty -> Signed.long -> 'a ty -> pari_ulong -> int
val flxq_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_lroot : 'a ty -> 'a ty -> Signed.long -> 'a ty
val flxq_lroot_pre : 'a ty -> 'a ty -> Signed.long -> pari_ulong -> 'a ty
val flxq_lroot_fast : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val flxq_lroot_fast_pre :
  'a ty -> 'a ty -> 'a ty -> Signed.long -> pari_ulong -> 'a ty

val flxq_matrix_pow :
  'a ty -> Signed.long -> Signed.long -> 'a ty -> pari_ulong -> 'a ty

val flxq_matrix_pow_pre :
  'a ty ->
  Signed.long ->
  Signed.long ->
  'a ty ->
  pari_ulong ->
  pari_ulong ->
  'a ty

val flxq_minpoly : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_minpoly_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxq_mul_pre : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxq_norm : 'a ty -> 'a ty -> pari_ulong -> pari_ulong
val flxq_order : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_pow_pre : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxq_pow_init :
  'a ty -> 'a ty -> Signed.long -> 'a ty -> pari_ulong -> 'a ty

val flxq_pow_init_pre :
  'a ty -> 'a ty -> Signed.long -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxq_pow_table_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxq_pow_table : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_powu : 'a ty -> pari_ulong -> 'a ty -> pari_ulong -> 'a ty

val flxq_powu_pre :
  'a ty -> pari_ulong -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxq_powers : 'a ty -> Signed.long -> 'a ty -> pari_ulong -> 'a ty

val flxq_powers_pre :
  'a ty -> Signed.long -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxq_sqr : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_sqr_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxq_sqrt : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_sqrt_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxq_sqrtn :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> 'a ty

val flxq_trace : 'a ty -> 'a ty -> pari_ulong -> pari_ulong
val flxqc_flxq_mul : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_flxq_mul : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqv_dotproduct : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqv_dotproduct_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val rg_to_f2 : 'a ty -> pari_ulong
val rg_to_fl : 'a ty -> pari_ulong -> pari_ulong
val rg_to_flxq : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val rgx_to_flx : 'a ty -> pari_ulong -> 'a ty
val rgxv_to_flxv : 'a ty -> pari_ulong -> 'a ty
val z_to_flx : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val zxv_to_flxv : 'a ty -> pari_ulong -> 'a ty
val zxt_to_flxt : 'a ty -> pari_ulong -> 'a ty
val gener_flxq : 'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> 'a ty

val get_flxq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  pari_ulong ->
  bb_field Ctypes.structure Ctypes_static.ptr

val monomial_flx : pari_ulong -> Signed.long -> Signed.long -> 'a ty
val random_flx : Signed.long -> Signed.long -> pari_ulong -> 'a ty
val zero_flxc : Signed.long -> Signed.long -> 'a ty
val zero_flxm : Signed.long -> Signed.long -> Signed.long -> 'a ty
val zlx_translate1 : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val zx_to_flx : 'a ty -> pari_ulong -> 'a ty
val flxx_fl_mul : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxx_flx_add : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxx_flx_mul : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxx_flx_sub : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxx_laplace : 'a ty -> pari_ulong -> 'a ty
val flxx_add : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxx_blocks : 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val flxx_deriv : 'a ty -> pari_ulong -> 'a ty
val flxx_double : 'a ty -> pari_ulong -> 'a ty
val flxx_invlaplace : 'a ty -> pari_ulong -> 'a ty
val flxx_neg : 'a ty -> pari_ulong -> 'a ty
val flxx_renormalize : 'a ty -> Signed.long -> 'a ty
val flxx_shift : 'a ty -> Signed.long -> Signed.long -> 'a ty
val flxx_sub : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxx_swap : 'a ty -> Signed.long -> Signed.long -> 'a ty
val flxx_to_flm : 'a ty -> Signed.long -> 'a ty
val flxx_to_flx : 'a ty -> 'a ty
val flxx_to_flxc : 'a ty -> Signed.long -> Signed.long -> 'a ty
val flxx_to_zxx : 'a ty -> 'a ty
val flxx_translate1 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val flxx_triple : 'a ty -> pari_ulong -> 'a ty
val flxxc_sub : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxxc_to_zxxc : 'a ty -> 'a ty
val flxxm_to_zxxm : 'a ty -> 'a ty
val flxxv_to_flxm : 'a ty -> Signed.long -> Signed.long -> 'a ty
val flxxn_red : 'a ty -> Signed.long -> 'a ty
val flxy_flx_div : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxy_flx_translate : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxy_flxqv_evalx : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxy_flxqv_evalx_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxy_flxq_evalx : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxy_flxq_evalx_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxy_evalx : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxy_evalx_powers_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxy_evalx_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty
val flxyqq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqv_roots_to_pol : 'a ty -> 'a ty -> pari_ulong -> Signed.long -> 'a ty

val flxqxc_flxqxqv_eval :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxc_flxqxqv_eval_pre :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxc_flxqxq_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_autpow :
  'a ty -> Signed.long -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_autpow_pre :
  'a ty -> Signed.long -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_autsum :
  'a ty -> Signed.long -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_autsum_pre :
  'a ty -> Signed.long -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_auttrace :
  'a ty -> pari_ulong -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_auttrace_pre :
  'a ty -> pari_ulong -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_div : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_div_pre :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_inv : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_inv_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_invsafe : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_invsafe_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_matrix_pow :
  'a ty -> Signed.long -> Signed.long -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_minpoly : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_minpoly_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_mul_pre :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_pow_pre :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_powers :
  'a ty -> Signed.long -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_powers_pre :
  'a ty -> Signed.long -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_powu_pre :
  'a ty -> pari_ulong -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxq_sqr : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqxq_sqr_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxv_prod : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_flxqxqv_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqx_flxqxqv_eval_pre :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqx_flxqxq_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqx_flxqxq_eval_pre :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqx_flxq_mul : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqx_flxq_mul_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqx_flxq_mul_to_monic : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqx_flxq_mul_to_monic_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqx_newton : 'a ty -> Signed.long -> 'a ty -> pari_ulong -> 'a ty

val flxqx_newton_pre :
  'a ty -> Signed.long -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqx_composedsum : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_disc : 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqx_div_by_x_x :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> 'a ty

val flxqx_div_by_x_x_pre :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flxqx_divrem :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> 'a ty

val flxqx_divrem_pre :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  Signed.long ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flxqx_dotproduct : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqx_extgcd :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flxqx_extgcd_pre :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flxqx_fromnewton : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_fromnewton_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqx_gcd : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_gcd_pre : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqx_get_red : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_get_red_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqx_halfgcd : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqx_halfgcd_all :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flxqx_halfgcd_all_pre :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flxqx_halfgcd_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqx_invbarrett : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_invbarrett_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqx_mul : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_mul_pre : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqx_normalize : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_normalize_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqx_powu : 'a ty -> pari_ulong -> 'a ty -> pari_ulong -> 'a ty

val flxqx_powu_pre :
  'a ty -> pari_ulong -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqx_red : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_red_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqx_rem : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_rem_pre : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqx_resultant : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqx_resultant_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqx_safegcd : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_saferesultant : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_sqr : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_sqr_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqxn_expint : 'a ty -> Signed.long -> 'a ty -> pari_ulong -> 'a ty

val flxqxn_expint_pre :
  'a ty -> Signed.long -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxn_inv : 'a ty -> Signed.long -> 'a ty -> pari_ulong -> 'a ty

val flxqxn_inv_pre :
  'a ty -> Signed.long -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxn_mul : 'a ty -> 'a ty -> Signed.long -> 'a ty -> pari_ulong -> 'a ty

val flxqxn_mul_pre :
  'a ty -> 'a ty -> Signed.long -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqxn_sqr : 'a ty -> Signed.long -> 'a ty -> pari_ulong -> 'a ty

val flxqxn_sqr_pre :
  'a ty -> Signed.long -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxy_degreex : 'a ty -> Signed.long

val flxy_eval_powers_pre :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong

val fly_to_flxy : 'a ty -> Signed.long -> 'a ty
val kronecker_to_flxqx : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val kronecker_to_flxqx_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val rgx_to_flxqx : 'a ty -> 'a ty -> pari_ulong -> 'a ty

val get_flxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_flxqx : Signed.long -> Signed.long -> 'a ty -> pari_ulong -> 'a ty

val zlxx_translate1 :
  'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty

val zxx_to_kronecker : 'a ty -> 'a ty -> 'a ty
val flxq_ellcard : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxq_ellgens :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxq_ellgroup :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty

val flxq_elltwist :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  unit

val flxq_ellj : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxq_ellj_to_a4a6 :
  'a ty ->
  'a ty ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  unit

val flxqe_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqe_changepoint : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqe_changepointinv : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqe_dbl : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqe_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqe_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqe_neg : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqe_order : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqe_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqe_tatepairing :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqe_weilpairing :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqe_weilpairing_pre :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val zxx_to_flxx : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val zxxt_to_flxxt : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val zxxv_to_flxxv : 'a ty -> pari_ulong -> Signed.long -> 'a ty

val get_flxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_flxqe : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val random_flxqe : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val fl_elltrace : pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fl_elltrace_cm :
  Signed.long -> pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fp_ellcard : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_elldivpol : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val fp_ellgens : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fp_ellgroup :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val fp_ellj : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fp_ellj_to_a4a6 :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val fp_elljissupersingular : 'a ty -> 'a ty -> int

val fp_elltwist :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  unit

val fp_ffellcard : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpe_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpe_changepoint : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpe_changepointinv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpe_dbl : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpe_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpe_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpe_neg : 'a ty -> 'a ty -> 'a ty
val fpe_order : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpe_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpe_to_fpj : 'a ty -> 'a ty
val fpe_to_mod : 'a ty -> 'a ty -> 'a ty
val fpe_tatepairing : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpe_weilpairing : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpj_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpj_dbl : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpj_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpj_neg : 'a ty -> 'a ty -> 'a ty
val fpj_to_fpe : 'a ty -> 'a ty -> 'a ty
val fpxq_ellcard : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_ellcard_supersingular : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_elldivpol : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty

val fpxq_ellgens :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxq_ellgroup :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val fpxq_ellj : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_elljissupersingular : 'a ty -> 'a ty -> 'a ty -> int

val fpxq_elltwist :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  unit

val fpxqe_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqe_changepoint : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqe_changepointinv : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqe_dbl : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqe_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqe_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqe_neg : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqe_order : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqe_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxqe_tatepairing :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxqe_weilpairing :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fq_elljissupersingular : 'a ty -> 'a ty -> 'a ty -> int
val fq_ellcard_supersingular : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val rge_to_fpe : 'a ty -> 'a ty -> 'a ty
val rge_to_fpxqe : 'a ty -> 'a ty -> 'a ty -> 'a ty

val get_fpe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  bb_group Ctypes.structure Ctypes_static.ptr

val random_fpe : 'a ty -> 'a ty -> 'a ty -> 'a ty
val random_fpxqe : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_issquare : 'a ty -> 'a ty -> int
val fp_fpx_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_fpxq_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpv_fpm_polint : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val fpv_inv : 'a ty -> 'a ty -> 'a ty
val fpv_invvandermonde : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpv_polint : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val fpv_roots_to_pol : 'a ty -> 'a ty -> Signed.long -> 'a ty
val fpx_fp_add : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fp_add_shallow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fp_div : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fp_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fp_mul_to_monic : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fp_mulspec : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val fpx_fp_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fp_sub_shallow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fpv_multieval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fpxq_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fpxqv_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_fpxv_multirem : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_frobenius : 'a ty -> 'a ty -> 'a ty
val fpx_laplace : 'a ty -> 'a ty -> 'a ty
val fpx_newton : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpx_add : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_center : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_center_i : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpx_chinese_coprime :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpx_composedprod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_composedsum : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_convol : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_deriv : 'a ty -> 'a ty -> 'a ty
val fpx_digits : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_disc : 'a ty -> 'a ty -> 'a ty
val fpx_div_by_x_x : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val fpx_divrem : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val fpx_divu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val fpx_dotproduct : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpx_extgcd :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val fpx_extresultant :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val fpx_fromnewton : 'a ty -> 'a ty -> 'a ty
val fpx_gcd : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_gcd_check : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_get_red : 'a ty -> 'a ty -> 'a ty
val fpx_halve : 'a ty -> 'a ty -> 'a ty
val fpx_halfgcd : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpx_halfgcd_all :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val fpx_integ : 'a ty -> 'a ty -> 'a ty
val fpx_invbarrett : 'a ty -> 'a ty -> 'a ty
val fpx_invlaplace : 'a ty -> 'a ty -> 'a ty
val fpx_is_squarefree : 'a ty -> 'a ty -> int
val fpx_matfrobenius : 'a ty -> 'a ty -> 'a ty
val fpx_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_mulspec : 'a ty -> 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val fpx_mulu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val fpx_neg : 'a ty -> 'a ty -> 'a ty
val fpx_normalize : 'a ty -> 'a ty -> 'a ty
val fpx_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val fpx_red : 'a ty -> 'a ty -> 'a ty
val fpx_rem : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_rescale : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_resultant : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_sqr : 'a ty -> 'a ty -> 'a ty
val fpx_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpx_valrem :
  'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

val fpxc_fpxq_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxc_fpxqv_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxm_fpxqv_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_autpow : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fpxq_autpowers : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fpxq_autsum : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fpxq_auttrace : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fpxq_charpoly : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_conjvec : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_div : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_invsafe : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_issquare : 'a ty -> 'a ty -> 'a ty -> int
val fpxq_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxq_matrix_pow :
  'a ty -> Signed.long -> Signed.long -> 'a ty -> 'a ty -> 'a ty

val fpxq_minpoly : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_norm : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_order : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fpxq_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fpxq_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_sqr : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_sqrt : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxq_sqrtn :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val fpxq_trace : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqc_to_mod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqm_autsum : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fpxt_red : 'a ty -> 'a ty -> 'a ty
val fpxv_fpx_fromdigits : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxv_chinese : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val fpxv_composedsum : 'a ty -> 'a ty -> 'a ty
val fpxv_factorback : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val fpxv_prod : 'a ty -> 'a ty -> 'a ty
val fpxv_red : 'a ty -> 'a ty -> 'a ty
val fpxn_div : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpxn_exp : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpxn_expint : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpxn_inv : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpxn_mul : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpxn_sqr : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fq_issquare : 'a ty -> 'a ty -> 'a ty -> int
val fq_ispower : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long
val fq_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqc_to_mod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_to_mod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqv_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val z_to_fpx : 'a ty -> 'a ty -> Signed.long -> 'a ty
val gener_fpxq : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val gener_fpxq_local : 'a ty -> 'a ty -> 'a ty -> 'a ty

val get_fpxq_star :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpx : Signed.long -> Signed.long -> 'a ty -> 'a ty
val f2x_ddf : 'a ty -> 'a ty
val f2x_factor : 'a ty -> 'a ty
val f2x_factor_squarefree : 'a ty -> 'a ty
val f2x_is_irred : 'a ty -> int
val flx_ddf : 'a ty -> pari_ulong -> 'a ty
val flx_ddf_pre : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_is_irred : 'a ty -> pari_ulong -> int
val flx_is_totally_split : 'a ty -> pari_ulong -> int

val flx_ispower :
  'a ty -> pari_ulong -> pari_ulong -> 'a ty Ctypes_static.ptr -> Signed.long

val flx_degfact : 'a ty -> pari_ulong -> 'a ty
val flx_factor : 'a ty -> pari_ulong -> 'a ty
val flx_factor_squarefree : 'a ty -> pari_ulong -> 'a ty
val flx_factor_squarefree_pre : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_nbfact : 'a ty -> pari_ulong -> Signed.long
val flx_nbfact_pre : 'a ty -> pari_ulong -> pari_ulong -> Signed.long
val flx_nbfact_frobenius : 'a ty -> 'a ty -> pari_ulong -> Signed.long

val flx_nbfact_frobenius_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> Signed.long

val flx_nbfact_by_degree :
  'a ty -> Signed.long Ctypes_static.ptr -> pari_ulong -> 'a ty

val flx_nbroots : 'a ty -> pari_ulong -> Signed.long
val flx_oneroot : 'a ty -> pari_ulong -> pari_ulong
val flx_oneroot_split : 'a ty -> pari_ulong -> pari_ulong
val flx_oneroot_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong
val flx_oneroot_split_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong
val flx_roots : 'a ty -> pari_ulong -> 'a ty
val flx_roots_pre : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flx_rootsff : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val fpx_ddf : 'a ty -> 'a ty -> 'a ty
val fpx_ddf_degree : 'a ty -> 'a ty -> 'a ty -> Signed.long
val fpx_degfact : 'a ty -> 'a ty -> 'a ty
val fpx_factor : 'a ty -> 'a ty -> 'a ty
val fpx_factor_squarefree : 'a ty -> 'a ty -> 'a ty
val fpx_is_irred : 'a ty -> 'a ty -> int
val fpx_is_totally_split : 'a ty -> 'a ty -> int

val fpx_ispower :
  'a ty -> pari_ulong -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

val fpx_nbfact : 'a ty -> 'a ty -> Signed.long
val fpx_nbfact_frobenius : 'a ty -> 'a ty -> 'a ty -> Signed.long
val fpx_nbroots : 'a ty -> 'a ty -> Signed.long
val fpx_oneroot : 'a ty -> 'a ty -> 'a ty
val fpx_oneroot_split : 'a ty -> 'a ty -> 'a ty
val fpx_roots : 'a ty -> 'a ty -> 'a ty
val fpx_roots_mult : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpx_rootsff : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_split_part : 'a ty -> 'a ty -> 'a ty
val f2xqx_ddf : 'a ty -> 'a ty -> 'a ty
val f2xqx_degfact : 'a ty -> 'a ty -> 'a ty
val f2xqx_factor : 'a ty -> 'a ty -> 'a ty
val f2xqx_factor_squarefree : 'a ty -> 'a ty -> 'a ty
val f2xqx_roots : 'a ty -> 'a ty -> 'a ty
val flx_factorff_irred : 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flx_ffintersect :
  'a ty ->
  'a ty ->
  Signed.long ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  unit

val flx_ffisom : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_ffisom_inv : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_frobenius : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_frobenius_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flxqx_ddf : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_ddf_degree : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> Signed.long
val flxqx_degfact : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_factor : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_factor_squarefree : 'a ty -> 'a ty -> pari_ulong -> 'a ty

val flxqx_factor_squarefree_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flxqx_ispower :
  'a ty ->
  pari_ulong ->
  'a ty ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  Signed.long

val flxqx_is_squarefree : 'a ty -> 'a ty -> pari_ulong -> Signed.long
val flxqx_nbfact : 'a ty -> 'a ty -> pari_ulong -> Signed.long

val flxqx_nbfact_frobenius :
  'a ty -> 'a ty -> 'a ty -> pari_ulong -> Signed.long

val flxqx_nbfact_by_degree :
  'a ty -> Signed.long Ctypes_static.ptr -> 'a ty -> pari_ulong -> 'a ty

val flxqx_nbroots : 'a ty -> 'a ty -> pari_ulong -> Signed.long
val flxqx_roots : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqxq_halffrobenius : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val fpx_factorff : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_factorff_irred : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpx_ffintersect :
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  unit

val fpx_ffisom : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_ffisom_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_frobenius : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_ddf : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_ddf_degree : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long
val fpxqx_degfact : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_factor : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_factor_squarefree : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxqx_ispower :
  'a ty ->
  pari_ulong ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  Signed.long

val fpxqx_nbfact : 'a ty -> 'a ty -> 'a ty -> Signed.long
val fpxqx_nbfact_frobenius : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long
val fpxqx_nbroots : 'a ty -> 'a ty -> 'a ty -> Signed.long
val fpxqx_roots : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_split_part : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_halffrobenius : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_is_squarefree : 'a ty -> 'a ty -> 'a ty -> Signed.long

val fqx_ispower :
  'a ty ->
  pari_ulong ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  Signed.long

val fqx_nbfact : 'a ty -> 'a ty -> 'a ty -> Signed.long
val fqx_nbroots : 'a ty -> 'a ty -> 'a ty -> Signed.long
val factorff : 'a ty -> 'a ty -> 'a ty -> 'a ty
val factormod0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val factormodddf : 'a ty -> 'a ty -> 'a ty
val factormodsqf : 'a ty -> 'a ty -> 'a ty

val ff_parse_tp :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long ->
  int

val rootmod0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val fpxqx_fpxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_fpxqxqv_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_fpxqxq_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_digits : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_disc : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxqx_div_by_x_x :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val fpxqx_divrem :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val fpxqx_dotproduct : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxqx_extgcd :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val fpxqx_gcd : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_get_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_halfgcd : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxqx_halfgcd_all :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val fpxqx_invbarrett : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fpxqx_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_rem : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_resultant : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_sqr : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_to_mod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_div : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_invsafe : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxqxq_matrix_pow :
  'a ty -> Signed.long -> Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fpxqxq_minpoly : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_sqr : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_autpow : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_autsum : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxq_auttrace : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxt_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxv_fpxqx_fromdigits : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxv_prod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxv_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqxn_div : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fpxqxn_exp : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fpxqxn_expint : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fpxqxn_inv : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fpxqxn_mul : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fpxqxn_sqr : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fpxx_fp_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxx_fpx_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxx_add : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxx_deriv : 'a ty -> 'a ty -> 'a ty
val fpxx_halve : 'a ty -> 'a ty -> 'a ty
val fpxx_integ : 'a ty -> 'a ty -> 'a ty
val fpxx_mulu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val fpxx_neg : 'a ty -> 'a ty -> 'a ty
val fpxx_red : 'a ty -> 'a ty -> 'a ty
val fpxx_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxy_fpxq_evalx : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxy_fpxqv_evalx : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxy_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxy_evalx : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxy_evaly : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val fpxyqq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxc_to_mod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxm_to_mod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val kronecker_to_fpxqx : 'a ty -> 'a ty -> 'a ty -> 'a ty

val get_fpxqx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpxqx : Signed.long -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val flc_flv_mul : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flc_to_mod : 'a ty -> pari_ulong -> 'a ty
val flm_fl_add : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flm_fl_mul : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flm_fl_mul_inplace : 'a ty -> pari_ulong -> pari_ulong -> unit
val flm_fl_mul_pre : 'a ty -> pari_ulong -> pari_ulong -> pari_ulong -> 'a ty
val flm_fl_sub : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flm_flc_mul : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_flc_mul_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty

val flm_flc_mul_pre_flx :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> Signed.long -> 'a ty

val flm_add : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_center : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flm_mul : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_mul_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flm_neg : 'a ty -> pari_ulong -> 'a ty
val flm_powers : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flm_powu : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flm_sub : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flm_to_mod : 'a ty -> pari_ulong -> 'a ty
val flm_transpose : 'a ty -> 'a ty
val flv_fl_div : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flv_fl_div_inplace : 'a ty -> pari_ulong -> pari_ulong -> unit
val flv_fl_mul : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flv_fl_mul_inplace : 'a ty -> pari_ulong -> pari_ulong -> unit

val flv_fl_mul_part_inplace :
  'a ty -> pari_ulong -> pari_ulong -> Signed.long -> unit

val flv_add : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flv_add_inplace : 'a ty -> 'a ty -> pari_ulong -> unit
val flv_center : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val flv_dotproduct : 'a ty -> 'a ty -> pari_ulong -> pari_ulong

val flv_dotproduct_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong

val flv_neg : 'a ty -> pari_ulong -> 'a ty
val flv_neg_inplace : 'a ty -> pari_ulong -> unit
val flv_sub : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flv_sub_inplace : 'a ty -> 'a ty -> pari_ulong -> unit
val flv_sum : 'a ty -> pari_ulong -> pari_ulong
val flx_dotproduct : 'a ty -> 'a ty -> pari_ulong -> pari_ulong

val flx_dotproduct_pre :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong -> pari_ulong

val fp_to_mod : 'a ty -> 'a ty -> 'a ty
val fpc_fpv_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpc_fp_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpc_center : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpc_center_inplace : 'a ty -> 'a ty -> 'a ty -> unit
val fpc_red : 'a ty -> 'a ty -> 'a ty
val fpc_to_mod : 'a ty -> 'a ty -> 'a ty
val fpm_add : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_fp_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_fpc_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_fpc_mul_fpx : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val fpm_center : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_center_inplace : 'a ty -> 'a ty -> 'a ty -> unit
val fpm_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val fpm_red : 'a ty -> 'a ty -> 'a ty
val fpm_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_to_mod : 'a ty -> 'a ty -> 'a ty
val fpms_fpc_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpms_fpcs_solve : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpms_fpcs_solve_safe : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpms_leftkernel_elt : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fpc_add : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpc_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpv_fpms_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpv_add : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpv_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpv_dotproduct : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpv_dotsquare : 'a ty -> 'a ty -> 'a ty
val fpv_red : 'a ty -> 'a ty -> 'a ty
val fpv_to_mod : 'a ty -> 'a ty -> 'a ty
val fpvv_to_mod : 'a ty -> 'a ty -> 'a ty
val fpx_to_mod : 'a ty -> 'a ty -> 'a ty
val fpxc_to_mod : 'a ty -> 'a ty -> 'a ty
val fpxm_to_mod : 'a ty -> 'a ty -> 'a ty
val zabm_ker : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zabm_indexrank : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zabm_inv : 'a ty -> 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr -> 'a ty

val zabm_inv_ratlift :
  'a ty -> 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr -> 'a ty

val zabm_pseudoinv :
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val zv_zms_mul : 'a ty -> 'a ty -> 'a ty

val zpms_zpcs_solve :
  'a ty -> 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty

val gen_fpm_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  'a ty

val gen_zpm_dixon_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val gen_matid :
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  'a ty

val matid_flm : Signed.long -> 'a ty
val matid_f2xqm : Signed.long -> 'a ty -> 'a ty
val matid_flxqm : Signed.long -> 'a ty -> pari_ulong -> 'a ty
val random_flv : Signed.long -> pari_ulong -> 'a ty
val random_fpc : Signed.long -> 'a ty -> 'a ty
val random_fpv : Signed.long -> 'a ty -> 'a ty
val scalar_flm : Signed.long -> Signed.long -> 'a ty
val zcs_to_zc : 'a ty -> Signed.long -> 'a ty
val zms_to_zm : 'a ty -> Signed.long -> 'a ty
val zms_zc_mul : 'a ty -> 'a ty -> 'a ty
val zmv_to_flmv : 'a ty -> pari_ulong -> 'a ty
val flx_teichmuller : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val z2_sqrt : 'a ty -> Signed.long -> 'a ty
val zp_div : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zp_exp : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val zp_inv : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zp_invlift : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zp_log : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val zp_sqrt : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zp_sqrtlift : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zp_sqrtnlift : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpm_invlift : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpx_frobenius : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpx_zpxq_liftroot : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zpx_zpxq_liftroot_ea :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val zpx_liftfact : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpx_liftroot : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpx_liftroots : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpx_roots : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpxq_div : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpxq_inv : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpxq_invlift : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpxq_log : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpxq_sqrt : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zpxq_sqrtnlift :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zpxqm_prodfrobenius : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zpxqx_digits :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zpxqx_divrem :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty Ctypes_static.ptr ->
  'a ty

val zpxqx_liftfact :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zpxqx_liftroot : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zpxqx_liftroot_vald :
  'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zpxqx_liftroots : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpxqx_roots : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zpxqx_zpxqxq_liftroot :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zq_sqrtnlift :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zqx_zqxq_liftroot :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zqx_liftfact :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val zqx_liftroot : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zqx_roots : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val gen_zpm_dixon :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty

val gen_zpm_newton :
  'a ty ->
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_zpx_dixon :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty

val gen_zpx_newton :
  'a ty ->
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val quadratic_prec_mask : Signed.long -> pari_ulong
val qx_factor : 'a ty -> 'a ty
val zx_factor : 'a ty -> 'a ty
val zx_is_irred : 'a ty -> Signed.long
val zx_squff : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val rg_rgc_sub : 'a ty -> 'a ty -> 'a ty
val rgc_rg_add : 'a ty -> 'a ty -> 'a ty
val rgc_rg_div : 'a ty -> 'a ty -> 'a ty
val rgc_rg_mul : 'a ty -> 'a ty -> 'a ty
val rgc_rg_sub : 'a ty -> 'a ty -> 'a ty
val rgc_rgm_mul : 'a ty -> 'a ty -> 'a ty
val rgc_rgv_mul : 'a ty -> 'a ty -> 'a ty
val rgc_add : 'a ty -> 'a ty -> 'a ty
val rgc_is_ei : 'a ty -> Signed.long
val rgc_neg : 'a ty -> 'a ty
val rgc_sub : 'a ty -> 'a ty -> 'a ty
val rgm_rg_add : 'a ty -> 'a ty -> 'a ty
val rgm_rg_add_shallow : 'a ty -> 'a ty -> 'a ty
val rgm_rg_div : 'a ty -> 'a ty -> 'a ty
val rgm_rg_mul : 'a ty -> 'a ty -> 'a ty
val rgm_rg_sub : 'a ty -> 'a ty -> 'a ty
val rgm_rg_sub_shallow : 'a ty -> 'a ty -> 'a ty
val rgm_rgc_mul : 'a ty -> 'a ty -> 'a ty
val rgm_rgv_mul : 'a ty -> 'a ty -> 'a ty
val rgm_add : 'a ty -> 'a ty -> 'a ty
val rgm_det_triangular : 'a ty -> 'a ty
val rgm_is_qm : 'a ty -> int
val rgm_is_zm : 'a ty -> int
val rgm_isdiagonal : 'a ty -> int
val rgm_isidentity : 'a ty -> int
val rgm_isscalar : 'a ty -> 'a ty -> int
val rgm_mul : 'a ty -> 'a ty -> 'a ty
val rgm_multosym : 'a ty -> 'a ty -> 'a ty
val rgm_neg : 'a ty -> 'a ty
val rgm_powers : 'a ty -> Signed.long -> 'a ty
val rgm_sqr : 'a ty -> 'a ty
val rgm_sub : 'a ty -> 'a ty -> 'a ty
val rgm_sumcol : 'a ty -> 'a ty
val rgm_transmul : 'a ty -> 'a ty -> 'a ty
val rgm_transmultosym : 'a ty -> 'a ty -> 'a ty
val rgmrow_zc_mul : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgm_zc_mul : 'a ty -> 'a ty -> 'a ty
val rgm_zm_mul : 'a ty -> 'a ty -> 'a ty
val rgmrow_rgc_mul : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgv_rgm_mul : 'a ty -> 'a ty -> 'a ty
val rgv_rgc_mul : 'a ty -> 'a ty -> 'a ty
val rgv_rg_mul : 'a ty -> 'a ty -> 'a ty
val rgv_add : 'a ty -> 'a ty -> 'a ty
val rgv_dotproduct : 'a ty -> 'a ty -> 'a ty
val rgv_dotsquare : 'a ty -> 'a ty
val rgv_is_zmv : 'a ty -> int
val rgv_kill0 : 'a ty -> 'a ty
val rgv_neg : 'a ty -> 'a ty
val rgv_prod : 'a ty -> 'a ty
val rgv_sub : 'a ty -> 'a ty -> 'a ty
val rgv_sum : 'a ty -> 'a ty
val rgv_sumpart : 'a ty -> Signed.long -> 'a ty
val rgv_sumpart2 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgv_zc_mul : 'a ty -> 'a ty -> 'a ty
val rgv_zm_mul : 'a ty -> 'a ty -> 'a ty
val rgx_rgm_eval : 'a ty -> 'a ty -> 'a ty
val rgx_rgmv_eval : 'a ty -> 'a ty -> 'a ty
val isdiagonal : 'a ty -> int
val scalarcol : 'a ty -> Signed.long -> 'a ty
val scalarcol_shallow : 'a ty -> Signed.long -> 'a ty
val scalarmat : 'a ty -> Signed.long -> 'a ty
val scalarmat_shallow : 'a ty -> Signed.long -> 'a ty
val scalarmat_s : Signed.long -> Signed.long -> 'a ty
val kronecker_to_mod : 'a ty -> 'a ty -> 'a ty
val qx_zxqv_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxq_charpoly : 'a ty -> 'a ty -> Signed.long -> 'a ty
val qxq_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty
val qxq_to_mod_shallow : 'a ty -> 'a ty -> 'a ty
val qxqc_to_mod_shallow : 'a ty -> 'a ty -> 'a ty
val qxqm_to_mod_shallow : 'a ty -> 'a ty -> 'a ty
val qxqv_to_mod : 'a ty -> 'a ty -> 'a ty
val qxqx_homogenous_evalpow : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxqx_to_mod_shallow : 'a ty -> 'a ty -> 'a ty
val qxqxv_to_mod : 'a ty -> 'a ty -> 'a ty
val qxv_qxq_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxy_qxq_evalx : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rg_rgx_sub : 'a ty -> 'a ty -> 'a ty
val rg_get_0 : 'a ty -> 'a ty
val rg_get_1 : 'a ty -> 'a ty
val rg_to_rgc : 'a ty -> Signed.long -> 'a ty
val rgm_to_rgxv : 'a ty -> Signed.long -> 'a ty
val rgm_to_rgxv_reverse : 'a ty -> Signed.long -> 'a ty
val rgm_to_rgxx : 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgv_to_rgx : 'a ty -> Signed.long -> 'a ty
val rgv_to_rgm : 'a ty -> Signed.long -> 'a ty
val rgv_to_rgx_reverse : 'a ty -> Signed.long -> 'a ty
val rgx_rgxq_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgx_rgxqv_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgx_rgxn_eval : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgx_rgxnv_eval : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgx_rg_add : 'a ty -> 'a ty -> 'a ty
val rgx_rg_add_shallow : 'a ty -> 'a ty -> 'a ty
val rgx_rg_div : 'a ty -> 'a ty -> 'a ty
val rgx_rg_divexact : 'a ty -> 'a ty -> 'a ty
val rgx_rg_eval_bk : 'a ty -> 'a ty -> 'a ty
val rgx_rg_mul : 'a ty -> 'a ty -> 'a ty
val rgx_rg_sub : 'a ty -> 'a ty -> 'a ty
val rgx_rgv_eval : 'a ty -> 'a ty -> 'a ty
val rgx_add : 'a ty -> 'a ty -> 'a ty
val rgx_addmulxn_shallow : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgx_addmulxn : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgx_addspec : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgx_addspec_shallow : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgx_affine : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgx_blocks : 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgx_deflate : 'a ty -> Signed.long -> 'a ty
val rgx_deriv : 'a ty -> 'a ty
val rgx_digits : 'a ty -> 'a ty -> 'a ty
val rgx_div_by_x_x : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val rgx_divrem : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val rgx_divs : 'a ty -> Signed.long -> 'a ty
val rgx_equal : 'a ty -> 'a ty -> Signed.long

val rgx_even_odd :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val rgx_homogenize : 'a ty -> Signed.long -> 'a ty
val rgx_homogenous_evalpow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgx_inflate : 'a ty -> Signed.long -> 'a ty
val rgx_mul : 'a ty -> 'a ty -> 'a ty
val rgx_mul_i : 'a ty -> 'a ty -> 'a ty
val rgx_mul_normalized : 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty
val rgx_mul2n : 'a ty -> Signed.long -> 'a ty
val rgx_mulxn : 'a ty -> Signed.long -> 'a ty
val rgx_mulhigh_i : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgx_muls : 'a ty -> Signed.long -> 'a ty
val rgx_mulspec : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgx_neg : 'a ty -> 'a ty
val rgx_normalize : 'a ty -> 'a ty
val rgx_pseudodivrem : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val rgx_pseudorem : 'a ty -> 'a ty -> 'a ty
val rgx_recip : 'a ty -> 'a ty
val rgx_recip_i : 'a ty -> 'a ty
val rgx_recip_shallow : 'a ty -> 'a ty
val rgx_rem : 'a ty -> 'a ty -> 'a ty
val rgx_renormalize_lg : 'a ty -> Signed.long -> 'a ty
val rgx_rescale : 'a ty -> 'a ty -> 'a ty
val rgx_rotate_shallow : 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgx_shift : 'a ty -> Signed.long -> 'a ty
val rgx_shift_shallow : 'a ty -> Signed.long -> 'a ty
val rgx_splitting : 'a ty -> Signed.long -> 'a ty
val rgx_sqr : 'a ty -> 'a ty
val rgx_sqr_i : 'a ty -> 'a ty
val rgx_sqrhigh_i : 'a ty -> Signed.long -> 'a ty
val rgx_sqrspec : 'a ty -> Signed.long -> 'a ty
val rgx_sub : 'a ty -> 'a ty -> 'a ty
val rgx_to_rgc : 'a ty -> Signed.long -> 'a ty
val rgx_translate : 'a ty -> 'a ty -> 'a ty
val rgx_unscale : 'a ty -> 'a ty -> 'a ty
val rgxq_matrix_pow : 'a ty -> Signed.long -> Signed.long -> 'a ty -> 'a ty
val rgxq_norm : 'a ty -> 'a ty -> 'a ty
val rgxq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgxq_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty
val rgxq_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val rgxq_trace : 'a ty -> 'a ty -> 'a ty
val rgxqc_red : 'a ty -> 'a ty -> 'a ty
val rgxqm_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgxqm_red : 'a ty -> 'a ty -> 'a ty
val rgxqv_rgxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgxqv_factorback : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgxqv_red : 'a ty -> 'a ty -> 'a ty
val rgxqx_rgxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgxqx_divrem : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val rgxqx_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgxqx_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty

val rgxqx_pseudodivrem :
  'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val rgxqx_pseudorem : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgxqx_red : 'a ty -> 'a ty -> 'a ty
val rgxqx_sqr : 'a ty -> 'a ty -> 'a ty
val rgxqx_translate : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgxv_rgv_eval : 'a ty -> 'a ty -> 'a ty
val rgxv_prod : 'a ty -> 'a ty
val rgxv_rescale : 'a ty -> 'a ty -> 'a ty
val rgxv_to_rgm : 'a ty -> Signed.long -> 'a ty
val rgxv_unscale : 'a ty -> 'a ty -> 'a ty
val rgxx_to_rgm : 'a ty -> Signed.long -> 'a ty
val rgxy_degreex : 'a ty -> Signed.long
val rgxy_derivx : 'a ty -> 'a ty
val rgxy_swap : 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgxy_swapspec : 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val rgxn_div : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgxn_div_i : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgxn_eval : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgxn_exp : 'a ty -> Signed.long -> 'a ty
val rgxn_expint : 'a ty -> Signed.long -> 'a ty
val rgxn_inv : 'a ty -> Signed.long -> 'a ty
val rgxn_inv_i : 'a ty -> Signed.long -> 'a ty
val rgxn_mul : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgxn_powers : 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgxn_recip_shallow : 'a ty -> Signed.long -> 'a ty
val rgxn_red_shallow : 'a ty -> Signed.long -> 'a ty
val rgxn_reverse : 'a ty -> Signed.long -> 'a ty
val rgxn_sqr : 'a ty -> Signed.long -> 'a ty
val rgxn_sqrt : 'a ty -> Signed.long -> 'a ty
val rgxnv_red_shallow : 'a ty -> Signed.long -> 'a ty
val rgxn_powu : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val rgxn_powu_i : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val zx_translate : 'a ty -> 'a ty -> 'a ty
val zx_unscale2n : 'a ty -> Signed.long -> 'a ty
val zx_unscale : 'a ty -> 'a ty -> 'a ty
val zx_unscale_div : 'a ty -> 'a ty -> 'a ty
val zx_unscale_divpow : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zx_z_unscale : 'a ty -> Signed.long -> 'a ty
val zxq_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty
val zxq_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val zxqx_dvd : 'a ty -> 'a ty -> 'a ty -> int
val brent_kung_optpow : Signed.long -> Signed.long -> Signed.long -> Signed.long

val gen_bkeval :
  'a ty ->
  Signed.long ->
  'a ty ->
  int ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_bkeval_powers :
  'a ty ->
  Signed.long ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val get_rg_algebra : unit -> bb_algebra Ctypes.structure Ctypes_static.ptr
val rfrac_deflate_order : 'a ty -> Signed.long
val rfrac_deflate_max : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val rfrac_deflate : 'a ty -> Signed.long -> 'a ty
val zgc_g_mul_inplace : 'a ty -> 'a ty -> unit
val zgcs_add : 'a ty -> 'a ty -> 'a ty
val g_zgc_mul : 'a ty -> 'a ty -> 'a ty
val g_zg_mul : 'a ty -> 'a ty -> 'a ty
val zgc_g_mul : 'a ty -> 'a ty -> 'a ty
val zgc_z_mul : 'a ty -> 'a ty -> 'a ty
val zg_g_mul : 'a ty -> 'a ty -> 'a ty
val zg_z_mul : 'a ty -> 'a ty -> 'a ty
val zg_add : 'a ty -> 'a ty -> 'a ty
val zg_mul : 'a ty -> 'a ty -> 'a ty
val zg_neg : 'a ty -> 'a ty
val zg_normalize : 'a ty -> 'a ty
val zg_sub : 'a ty -> 'a ty -> 'a ty
val flc_lincomb1_inplace : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> unit
val vecsmall_prod : 'a ty -> 'a ty
val qm_qc_mul : 'a ty -> 'a ty -> 'a ty
val qm_det : 'a ty -> 'a ty
val qm_ker : 'a ty -> 'a ty
val qm_mul : 'a ty -> 'a ty -> 'a ty
val qm_sqr : 'a ty -> 'a ty
val rgm_check_zm : 'a ty -> string -> unit
val rgv_check_zv : 'a ty -> string -> unit
val z_zc_sub : 'a ty -> 'a ty -> 'a ty
val zv_zc_mul : 'a ty -> 'a ty -> 'a ty
val zc_q_mul : 'a ty -> 'a ty -> 'a ty
val zc_z_add : 'a ty -> 'a ty -> 'a ty
val zc_z_div : 'a ty -> 'a ty -> 'a ty
val zc_z_divexact : 'a ty -> 'a ty -> 'a ty
val zc_z_sub : 'a ty -> 'a ty -> 'a ty
val zc_zv_mul : 'a ty -> 'a ty -> 'a ty
val zc_divexactu : 'a ty -> pari_ulong -> 'a ty
val zc_add : 'a ty -> 'a ty -> 'a ty
val zc_copy : 'a ty -> 'a ty
val zc_hnfremdiv : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zc_is_ei : 'a ty -> Signed.long
val zc_lincomb : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val zc_lincomb1_inplace : 'a ty -> 'a ty -> 'a ty -> unit
val zc_lincomb1_inplace_i : 'a ty -> 'a ty -> 'a ty -> Signed.long -> unit
val zc_neg : 'a ty -> 'a ty
val zc_reducemodlll : 'a ty -> 'a ty -> 'a ty
val zc_reducemodmatrix : 'a ty -> 'a ty -> 'a ty
val zc_sub : 'a ty -> 'a ty -> 'a ty
val zc_z_mul : 'a ty -> Signed.long -> 'a ty
val zm_q_mul : 'a ty -> 'a ty -> 'a ty
val zm_z_div : 'a ty -> 'a ty -> 'a ty
val zm_z_divexact : 'a ty -> 'a ty -> 'a ty
val zm_z_mul : 'a ty -> 'a ty -> 'a ty
val zm_add : 'a ty -> 'a ty -> 'a ty
val zm_det_triangular : 'a ty -> 'a ty
val zm_diag_mul : 'a ty -> 'a ty -> 'a ty
val zm_divexactu : 'a ty -> pari_ulong -> 'a ty
val zm_equal : 'a ty -> 'a ty -> int
val zm_equal0 : 'a ty -> int
val zm_hnfdivrem : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zm_ishnf : 'a ty -> int
val zm_isdiagonal : 'a ty -> int
val zm_isidentity : 'a ty -> int
val zm_isscalar : 'a ty -> 'a ty -> int
val zm_max_lg : 'a ty -> Signed.long
val zm_mul_diag : 'a ty -> 'a ty -> 'a ty
val zm_multosym : 'a ty -> 'a ty -> 'a ty
val zm_neg : 'a ty -> 'a ty
val zm_nm_mul : 'a ty -> 'a ty -> 'a ty
val zm_pow : 'a ty -> 'a ty -> 'a ty
val zm_powu : 'a ty -> pari_ulong -> 'a ty
val zm_reducemodlll : 'a ty -> 'a ty -> 'a ty
val zm_reducemodmatrix : 'a ty -> 'a ty -> 'a ty
val zm_sqr : 'a ty -> 'a ty
val zm_sub : 'a ty -> 'a ty -> 'a ty
val zm_supnorm : 'a ty -> 'a ty
val zm_transmul : 'a ty -> 'a ty -> 'a ty
val zm_transmultosym : 'a ty -> 'a ty -> 'a ty
val zm_togglesign : 'a ty -> unit
val zm_zm_mul : 'a ty -> 'a ty -> 'a ty
val zmrow_zc_mul : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zmrow_equal0 : 'a ty -> Signed.long -> int
val zv_abscmp : 'a ty -> 'a ty -> int
val zv_cmp : 'a ty -> 'a ty -> int
val zv_dotsquare : 'a ty -> 'a ty
val zv_max_lg : 'a ty -> Signed.long
val zv_to_nv : 'a ty -> 'a ty
val zv_togglesign : 'a ty -> unit
val gram_matrix : 'a ty -> 'a ty
val nm_z_mul : 'a ty -> 'a ty -> 'a ty
val zm_mul : 'a ty -> 'a ty -> 'a ty
val zm_to_flm : 'a ty -> pari_ulong -> 'a ty
val zm_to_zm : 'a ty -> 'a ty
val zm_zc_mul : 'a ty -> 'a ty -> 'a ty
val zmv_to_zmv : 'a ty -> 'a ty
val zv_abs : 'a ty -> 'a ty
val zv_content : 'a ty -> Signed.long
val zv_dotproduct : 'a ty -> 'a ty -> Signed.long
val zv_equal : 'a ty -> 'a ty -> int
val zv_equal0 : 'a ty -> int
val zv_neg : 'a ty -> 'a ty
val zv_neg_inplace : 'a ty -> 'a ty
val zv_prod : 'a ty -> Signed.long
val zv_prod_z : 'a ty -> 'a ty
val zv_sum : 'a ty -> Signed.long
val zv_sumpart : 'a ty -> Signed.long -> Signed.long
val zv_to_flv : 'a ty -> pari_ulong -> 'a ty
val zv_z_mul : 'a ty -> Signed.long -> 'a ty
val zv_zm_mul : 'a ty -> 'a ty -> 'a ty
val zvv_equal : 'a ty -> 'a ty -> int
val kronecker_to_zxqx : 'a ty -> 'a ty -> 'a ty
val kronecker_to_zxx : 'a ty -> Signed.long -> Signed.long -> 'a ty
val qx_zx_rem : 'a ty -> 'a ty -> 'a ty
val qx_mul : 'a ty -> 'a ty -> 'a ty
val qx_sqr : 'a ty -> 'a ty
val qxqm_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxqm_sqr : 'a ty -> 'a ty -> 'a ty
val qxqx_qxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxqx_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxqx_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty
val qxqx_sqr : 'a ty -> 'a ty -> 'a ty
val rgx_check_qx : 'a ty -> string -> unit
val rgx_check_zx : 'a ty -> string -> unit
val rgx_check_zxx : 'a ty -> string -> unit
val z_zx_sub : 'a ty -> 'a ty -> 'a ty
val zx_z_add : 'a ty -> 'a ty -> 'a ty
val zx_z_add_shallow : 'a ty -> 'a ty -> 'a ty
val zx_z_eval : 'a ty -> 'a ty -> 'a ty
val zx_z_mul : 'a ty -> 'a ty -> 'a ty
val zx_z_sub : 'a ty -> 'a ty -> 'a ty
val zx_add : 'a ty -> 'a ty -> 'a ty
val zx_affine : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zx_copy : 'a ty -> 'a ty
val zx_deriv : 'a ty -> 'a ty
val zx_digits : 'a ty -> 'a ty -> 'a ty
val zxv_zx_fromdigits : 'a ty -> 'a ty -> 'a ty
val zx_div_by_x_1 : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zx_divuexact : 'a ty -> pari_ulong -> 'a ty
val zx_equal : 'a ty -> 'a ty -> int
val zx_eval1 : 'a ty -> 'a ty
val zx_max_lg : 'a ty -> Signed.long
val zx_mod_xnm1 : 'a ty -> pari_ulong -> 'a ty
val zx_mul : 'a ty -> 'a ty -> 'a ty
val zx_mulspec : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val zx_mulu : 'a ty -> pari_ulong -> 'a ty
val zx_neg : 'a ty -> 'a ty
val zx_rem : 'a ty -> 'a ty -> 'a ty
val zx_remi2n : 'a ty -> Signed.long -> 'a ty
val zx_rescale2n : 'a ty -> Signed.long -> 'a ty
val zx_rescale : 'a ty -> 'a ty -> 'a ty
val zx_rescale_lt : 'a ty -> 'a ty
val zx_shifti : 'a ty -> Signed.long -> 'a ty
val zx_sqr : 'a ty -> 'a ty
val zx_sqrspec : 'a ty -> Signed.long -> 'a ty
val zx_sub : 'a ty -> 'a ty -> 'a ty
val zx_val : 'a ty -> Signed.long
val zx_valrem : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val zxc_to_flxc : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val zxm_to_flxm : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val zxqm_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zxqm_sqr : 'a ty -> 'a ty -> 'a ty
val zxqx_zxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zxqx_sqr : 'a ty -> 'a ty -> 'a ty
val zxqx_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zxt_remi2n : 'a ty -> Signed.long -> 'a ty
val zxv_z_mul : 'a ty -> 'a ty -> 'a ty
val zxv_dotproduct : 'a ty -> 'a ty -> 'a ty
val zxv_equal : 'a ty -> 'a ty -> int
val zxv_remi2n : 'a ty -> Signed.long -> 'a ty
val zxx_z_divexact : 'a ty -> 'a ty -> 'a ty
val zxx_z_mul : 'a ty -> 'a ty -> 'a ty
val zxx_z_add_shallow : 'a ty -> 'a ty -> 'a ty
val zxx_evalx0 : 'a ty -> 'a ty
val zxx_max_lg : 'a ty -> Signed.long
val zxx_mul_kronecker : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zxx_renormalize : 'a ty -> Signed.long -> 'a ty
val zxx_sqr_kronecker : 'a ty -> Signed.long -> 'a ty
val rgxx_to_kronecker : 'a ty -> Signed.long -> 'a ty
val rgxx_to_kronecker_spec : 'a ty -> Signed.long -> Signed.long -> 'a ty
val zxn_mul : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zxn_sqr : 'a ty -> Signed.long -> 'a ty
val scalar_zx : 'a ty -> Signed.long -> 'a ty
val scalar_zx_shallow : 'a ty -> Signed.long -> 'a ty
val zx_to_zx : 'a ty -> 'a ty
val zx_z_divexact : 'a ty -> Signed.long -> 'a ty
val alg_centralproj : 'a ty -> 'a ty -> Signed.long -> 'a ty
val alg_complete : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val alg_csa_table : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val alg_cyclic : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val alg_get_absdim : 'a ty -> Signed.long
val alg_get_abssplitting : 'a ty -> 'a ty
val alg_get_aut : 'a ty -> 'a ty
val algaut : 'a ty -> 'a ty
val alg_get_auts : 'a ty -> 'a ty
val alg_get_b : 'a ty -> 'a ty
val algb : 'a ty -> 'a ty
val algcenter : 'a ty -> 'a ty
val alg_get_center : 'a ty -> 'a ty
val alg_get_char : 'a ty -> 'a ty
val algchar : 'a ty -> 'a ty
val alg_get_degree : 'a ty -> Signed.long
val algdegree : 'a ty -> Signed.long
val alg_get_dim : 'a ty -> Signed.long
val algdim : 'a ty -> Signed.long -> Signed.long
val alg_get_hasse_f : 'a ty -> 'a ty
val alghassef : 'a ty -> 'a ty
val alg_get_hasse_i : 'a ty -> 'a ty
val alghassei : 'a ty -> 'a ty
val alg_get_invbasis : 'a ty -> 'a ty
val alginvbasis : 'a ty -> 'a ty
val alg_get_multable : 'a ty -> 'a ty
val alg_get_basis : 'a ty -> 'a ty
val algbasis : 'a ty -> 'a ty
val alg_get_relmultable : 'a ty -> 'a ty
val algrelmultable : 'a ty -> 'a ty
val alg_get_splitpol : 'a ty -> 'a ty
val alg_get_splittingfield : 'a ty -> 'a ty
val algsplittingfield : 'a ty -> 'a ty
val alg_get_splittingbasis : 'a ty -> 'a ty
val alg_get_splittingbasisinv : 'a ty -> 'a ty
val alg_get_splittingdata : 'a ty -> 'a ty
val algsplittingdata : 'a ty -> 'a ty
val alg_get_tracebasis : 'a ty -> 'a ty

val alg_hasse :
  'a ty -> Signed.long -> 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty

val alg_hilbert : 'a ty -> 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val alg_matrix : 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val alg_model : 'a ty -> 'a ty -> Signed.long
val alg_quotient : 'a ty -> 'a ty -> Signed.long -> 'a ty
val algradical : 'a ty -> 'a ty
val algsimpledec : 'a ty -> Signed.long -> 'a ty
val algsimpledec_ss : 'a ty -> Signed.long -> 'a ty
val algsubalg : 'a ty -> 'a ty -> 'a ty
val alg_type : 'a ty -> Signed.long
val algadd : 'a ty -> 'a ty -> 'a ty -> 'a ty
val algalgtobasis : 'a ty -> 'a ty -> 'a ty
val algbasistoalg : 'a ty -> 'a ty -> 'a ty
val algcharpoly : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val algdisc : 'a ty -> 'a ty
val algdivl : 'a ty -> 'a ty -> 'a ty -> 'a ty
val algdivr : 'a ty -> 'a ty -> 'a ty -> 'a ty
val alggroup : 'a ty -> 'a ty -> 'a ty
val alggroupcenter : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val alghasse : 'a ty -> 'a ty -> 'a ty
val alginit : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val algindex : 'a ty -> 'a ty -> Signed.long
val alginv : 'a ty -> 'a ty -> 'a ty
val algisassociative : 'a ty -> 'a ty -> int
val algiscommutative : 'a ty -> int
val algisdivision : 'a ty -> 'a ty -> int
val algisramified : 'a ty -> 'a ty -> int
val algissemisimple : 'a ty -> int
val algissimple : 'a ty -> Signed.long -> int
val algissplit : 'a ty -> 'a ty -> int
val algisdivl : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> int
val algisinv : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> int
val algmakeintegral : 'a ty -> Signed.long -> 'a ty
val algmul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val algmultable : 'a ty -> 'a ty
val alglat_get_primbasis : 'a ty -> 'a ty
val alglat_get_scalar : 'a ty -> 'a ty
val alglatadd : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val alglatcontains : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> int
val alglatelement : 'a ty -> 'a ty -> 'a ty -> 'a ty
val alglathnf : 'a ty -> 'a ty -> 'a ty -> 'a ty
val alglatindex : 'a ty -> 'a ty -> 'a ty -> 'a ty
val alglatinter : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val alglatmul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val alglatlefttransporter : 'a ty -> 'a ty -> 'a ty -> 'a ty
val alglatrighttransporter : 'a ty -> 'a ty -> 'a ty -> 'a ty
val alglatsubset : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> int
val algneg : 'a ty -> 'a ty -> 'a ty
val algnorm : 'a ty -> 'a ty -> Signed.long -> 'a ty
val algpoleval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val algpow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val algprimesubalg : 'a ty -> 'a ty
val algramifiedplaces : 'a ty -> 'a ty
val algrandom : 'a ty -> 'a ty -> 'a ty
val algsplit : 'a ty -> Signed.long -> 'a ty
val algtomatrix : 'a ty -> 'a ty -> Signed.long -> 'a ty
val algsqr : 'a ty -> 'a ty -> 'a ty
val algsub : 'a ty -> 'a ty -> 'a ty -> 'a ty
val algtableinit : 'a ty -> 'a ty -> 'a ty
val algtensor : 'a ty -> 'a ty -> Signed.long -> 'a ty
val algtrace : 'a ty -> 'a ty -> Signed.long -> 'a ty
val algtype : 'a ty -> Signed.long
val bnfgwgeneric : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val checkalg : 'a ty -> unit
val checkhasse : 'a ty -> 'a ty -> 'a ty -> Signed.long -> unit
val checklat : 'a ty -> 'a ty -> unit
val conjclasses_algcenter : 'a ty -> 'a ty -> 'a ty
val galoischardet : 'a ty -> 'a ty -> Signed.long -> 'a ty
val galoischarpoly : 'a ty -> 'a ty -> Signed.long -> 'a ty
val galoischartable : 'a ty -> 'a ty
val f2ms_colelim : 'a ty -> Signed.long -> 'a ty
val f2m_image : 'a ty -> 'a ty
val f2m_indexrank : 'a ty -> 'a ty
val f2m_suppl : 'a ty -> 'a ty
val f2xqm_f2xqc_gauss : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqm_f2xqc_invimage : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqm_f2xqc_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqm_deplin : 'a ty -> 'a ty -> 'a ty
val f2xqm_det : 'a ty -> 'a ty -> 'a ty
val f2xqm_gauss : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqm_ker : 'a ty -> 'a ty -> 'a ty
val f2xqm_image : 'a ty -> 'a ty -> 'a ty
val f2xqm_indexrank : 'a ty -> 'a ty -> 'a ty
val f2xqm_inv : 'a ty -> 'a ty -> 'a ty
val f2xqm_invimage : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqm_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val f2xqm_rank : 'a ty -> 'a ty -> Signed.long
val f2xqm_suppl : 'a ty -> 'a ty -> 'a ty
val flm_image : 'a ty -> pari_ulong -> 'a ty
val flm_indexrank : 'a ty -> pari_ulong -> 'a ty
val flm_suppl : 'a ty -> pari_ulong -> 'a ty
val flxqm_flxqc_gauss : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_flxqc_invimage : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_flxqc_mul : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_deplin : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_det : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_gauss : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_ker : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_image : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_indexrank : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_inv : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_invimage : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_mul : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqm_rank : 'a ty -> 'a ty -> pari_ulong -> Signed.long
val flxqm_suppl : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val fpm_fpc_gauss : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_fpc_invimage : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_deplin : 'a ty -> 'a ty -> 'a ty
val fpm_det : 'a ty -> 'a ty -> 'a ty
val fpm_gauss : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_image : 'a ty -> 'a ty -> 'a ty
val fpm_indexrank : 'a ty -> 'a ty -> 'a ty
val fpm_intersect : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_intersect_i : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_inv : 'a ty -> 'a ty -> 'a ty
val fpm_invimage : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_ker : 'a ty -> 'a ty -> 'a ty
val fpm_rank : 'a ty -> 'a ty -> Signed.long
val fpm_suppl : 'a ty -> 'a ty -> 'a ty
val fqm_fqc_gauss : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_fqc_invimage : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_fqc_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_deplin : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_det : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_gauss : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_ker : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_image : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_indexrank : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_invimage : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqm_rank : 'a ty -> 'a ty -> 'a ty -> Signed.long
val fqm_suppl : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qm_image_shallow : 'a ty -> 'a ty
val qm_image : 'a ty -> 'a ty
val qm_gauss : 'a ty -> 'a ty -> 'a ty
val qm_gauss_i : 'a ty -> 'a ty -> Signed.long -> 'a ty
val qm_indexrank : 'a ty -> 'a ty
val qm_inv : 'a ty -> 'a ty
val qm_rank : 'a ty -> Signed.long
val rgm_fp_init : 'a ty -> 'a ty -> pari_ulong Ctypes_static.ptr -> 'a ty
val rgm_hadamard : 'a ty -> 'a ty
val rgm_rgc_invimage : 'a ty -> 'a ty -> 'a ty
val rgm_diagonal : 'a ty -> 'a ty
val rgm_diagonal_shallow : 'a ty -> 'a ty
val rgm_inv : 'a ty -> 'a ty
val rgm_inv_upper : 'a ty -> 'a ty
val rgm_invimage : 'a ty -> 'a ty -> 'a ty
val rgm_solve : 'a ty -> 'a ty -> 'a ty
val rgm_solve_realimag : 'a ty -> 'a ty -> 'a ty

val rgms_structelim :
  'a ty ->
  Signed.long ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  unit

val zm_det : 'a ty -> 'a ty
val zm_detmult : 'a ty -> 'a ty
val zm_gauss : 'a ty -> 'a ty -> 'a ty
val zm_ker : 'a ty -> 'a ty
val zm_imagecompl : 'a ty -> 'a ty
val zm_indeximage : 'a ty -> 'a ty
val zm_indexrank : 'a ty -> 'a ty
val zm_inv : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zm_inv_ratlift : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val zm_pseudoinv :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val zm_rank : 'a ty -> Signed.long
val zlm_gauss : 'a ty -> 'a ty -> pari_ulong -> Signed.long -> 'a ty -> 'a ty
val closemodinvertible : 'a ty -> 'a ty -> 'a ty
val deplin : 'a ty -> 'a ty
val det : 'a ty -> 'a ty
val det0 : 'a ty -> Signed.long -> 'a ty
val det2 : 'a ty -> 'a ty
val detint : 'a ty -> 'a ty
val eigen : 'a ty -> Signed.long -> 'a ty
val gauss : 'a ty -> 'a ty -> 'a ty
val gaussmodulo : 'a ty -> 'a ty -> 'a ty -> 'a ty
val gaussmodulo2 : 'a ty -> 'a ty -> 'a ty -> 'a ty

val gen_gauss :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_gauss_pivot :
  'a ty ->
  Signed.long Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_det :
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_ker :
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_matcolinvimage :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_matcolmul :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_matinvimage :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_matmul :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  'a ty

val image : 'a ty -> 'a ty
val image2 : 'a ty -> 'a ty
val imagecompl : 'a ty -> 'a ty
val indexrank : 'a ty -> 'a ty
val inverseimage : 'a ty -> 'a ty -> 'a ty
val ker : 'a ty -> 'a ty
val mateigen : 'a ty -> Signed.long -> Signed.long -> 'a ty
val matimage0 : 'a ty -> Signed.long -> 'a ty
val matker0 : 'a ty -> Signed.long -> 'a ty
val rank : 'a ty -> Signed.long
val reducemodinvertible : 'a ty -> 'a ty -> 'a ty
val reducemodlll : 'a ty -> 'a ty -> 'a ty
val split_realimag : 'a ty -> Signed.long -> Signed.long -> 'a ty
val suppl : 'a ty -> 'a ty
val flm_charpoly : 'a ty -> pari_ulong -> 'a ty
val flm_hess : 'a ty -> pari_ulong -> 'a ty
val fpm_charpoly : 'a ty -> 'a ty -> 'a ty
val fpm_hess : 'a ty -> 'a ty -> 'a ty
val frobeniusform : 'a ty -> Signed.long -> 'a ty
val qm_minors_coprime : 'a ty -> 'a ty -> 'a ty
val qm_imz : 'a ty -> 'a ty

val qm_imz_all :
  'a ty -> 'a ty Ctypes_static.ptr -> Signed.long -> Signed.long -> 'a ty

val qm_imz_hnf : 'a ty -> 'a ty
val qm_imz_hnfall : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long -> 'a ty
val qm_imq : 'a ty -> 'a ty

val qm_imq_all :
  'a ty -> 'a ty Ctypes_static.ptr -> Signed.long -> Signed.long -> 'a ty

val qm_imq_hnf : 'a ty -> 'a ty
val qm_imq_hnfall : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long -> 'a ty
val qm_charpoly_zx : 'a ty -> 'a ty
val qm_charpoly_zx_bound : 'a ty -> Signed.long -> 'a ty
val zm_charpoly : 'a ty -> 'a ty
val adj : 'a ty -> 'a ty
val adjsafe : 'a ty -> 'a ty
val caract : 'a ty -> Signed.long -> 'a ty
val caradj : 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr -> 'a ty
val carberkowitz : 'a ty -> Signed.long -> 'a ty
val carhess : 'a ty -> Signed.long -> 'a ty
val charpoly : 'a ty -> Signed.long -> 'a ty
val charpoly0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val gnorm : 'a ty -> 'a ty
val gnorml1 : 'a ty -> Signed.long -> 'a ty
val gnorml1_fake : 'a ty -> 'a ty
val gnormlp : 'a ty -> 'a ty -> Signed.long -> 'a ty
val gnorml2 : 'a ty -> 'a ty
val gsupnorm : 'a ty -> Signed.long -> 'a ty

val gsupnorm_aux :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long ->
  unit

val gtrace : 'a ty -> 'a ty
val hess : 'a ty -> 'a ty
val intersect : 'a ty -> 'a ty -> 'a ty
val jacobi : 'a ty -> Signed.long -> 'a ty
val matadjoint0 : 'a ty -> Signed.long -> 'a ty
val matcompanion : 'a ty -> 'a ty
val matrixqz0 : 'a ty -> 'a ty -> 'a ty
val minpoly : 'a ty -> Signed.long -> 'a ty
val qfgaussred : 'a ty -> 'a ty
val qfgaussred_positive : 'a ty -> 'a ty
val qfsign : 'a ty -> 'a ty
val apply0 : 'a ty -> 'a ty -> 'a ty
val diagonal : 'a ty -> 'a ty
val diagonal_shallow : 'a ty -> 'a ty
val extract0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fold0 : 'a ty -> 'a ty -> 'a ty

val genapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val genfold :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val genindexselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val genselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val gtomat : 'a ty -> 'a ty
val gtrans : 'a ty -> 'a ty
val matmuldiagonal : 'a ty -> 'a ty -> 'a ty
val matmultodiagonal : 'a ty -> 'a ty -> 'a ty

val matslice0 :
  'a ty -> Signed.long -> Signed.long -> Signed.long -> Signed.long -> 'a ty

val parapply : 'a ty -> 'a ty -> 'a ty

val parfor :
  'a ty ->
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parfor_init :
  parfor_t Ctypes.structure Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty -> unit

val parfor_next : parfor_t Ctypes.structure Ctypes_static.ptr -> 'a ty
val parfor_stop : parfor_t Ctypes.structure Ctypes_static.ptr -> unit

val parforeach :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforeach_init :
  parforeach_t Ctypes.structure Ctypes_static.ptr -> 'a ty -> 'a ty -> unit

val parforeach_next : parforeach_t Ctypes.structure Ctypes_static.ptr -> 'a ty
val parforeach_stop : parforeach_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprime :
  'a ty ->
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprime_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  unit

val parforprime_next : parforprime_t Ctypes.structure Ctypes_static.ptr -> 'a ty
val parforprime_stop : parforprime_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprimestep :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprimestep_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  unit

val parforvec :
  'a ty ->
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforvec_init :
  parforvec_t Ctypes.structure Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  unit

val parforvec_next : parforvec_t Ctypes.structure Ctypes_static.ptr -> 'a ty
val parforvec_stop : parforvec_t Ctypes.structure Ctypes_static.ptr -> unit
val parselect : 'a ty -> 'a ty -> Signed.long -> 'a ty
val select0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val shallowextract : 'a ty -> 'a ty -> 'a ty
val shallowmatextract : 'a ty -> 'a ty -> 'a ty -> 'a ty
val shallowtrans : 'a ty -> 'a ty

val vecapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val veccatapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val veccatselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val vecrange : 'a ty -> 'a ty -> 'a ty
val vecrangess : Signed.long -> Signed.long -> 'a ty

val vecselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val vecselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val vecslice0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val vecsum : 'a ty -> 'a ty
val zv_diagonal : 'a ty -> 'a ty
val addhelp : string -> string -> unit
val arity0 : 'a ty -> 'a ty
val alias0 : string -> string -> unit
val compile_str : string -> 'a ty
val delete_var : unit -> Signed.long
val fetch_user_var : string -> Signed.long
val fetch_var : unit -> Signed.long
val fetch_var_higher : unit -> Signed.long
val fetch_var_value : Signed.long -> 'a ty -> 'a ty
val gp_embedded : string -> string
val gp_embedded_init : Signed.long -> Signed.long -> unit
val gp_read_str : string -> 'a ty
val gp_read_str_bitprec : string -> Signed.long -> 'a ty
val gp_read_str_prec : string -> Signed.long -> 'a ty

val install :
  unit Ctypes_static.ptr ->
  string ->
  string ->
  entree Ctypes.structure Ctypes_static.ptr

val is_entry : string -> entree Ctypes.structure Ctypes_static.ptr
val kill0 : string -> unit
val pari_var_close : unit -> unit
val pari_var_init : unit -> unit
val pari_var_next : unit -> Signed.long
val pari_var_next_temp : unit -> Signed.long
val pari_var_create : entree Ctypes.structure Ctypes_static.ptr -> Signed.long
val name_var : Signed.long -> string -> unit
val readseq : string -> 'a ty
val safeel : 'a ty -> Signed.long -> Signed.long Ctypes_static.ptr
val safelistel : 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr
val strtoi : string -> 'a ty
val strtor : string -> Signed.long -> 'a ty
val varhigher : string -> Signed.long -> 'a ty
val varlower : string -> Signed.long -> 'a ty
val divisorslenstra : 'a ty -> 'a ty -> 'a ty -> 'a ty
val isprimeaprcl : 'a ty -> Signed.long
val qfb0 : 'a ty -> 'a ty -> 'a ty -> 'a ty

val check_quaddisc :
  'a ty ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  string ->
  unit

val check_quaddisc_imag :
  'a ty -> Signed.long Ctypes_static.ptr -> string -> unit

val check_quaddisc_real :
  'a ty -> Signed.long Ctypes_static.ptr -> string -> unit

val cornacchia :
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long

val cornacchia2 :
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long

val cornacchia2_sqrt :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long

val nucomp : 'a ty -> 'a ty -> 'a ty -> 'a ty
val nudupl : 'a ty -> 'a ty -> 'a ty
val nupow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val primeform : 'a ty -> 'a ty -> 'a ty
val primeform_u : 'a ty -> pari_ulong -> 'a ty
val qfb_1 : 'a ty -> 'a ty
val qfbcomp : 'a ty -> 'a ty -> 'a ty
val qfbcomp_i : 'a ty -> 'a ty -> 'a ty
val qfbcompraw : 'a ty -> 'a ty -> 'a ty
val qfbcompraw_i : 'a ty -> 'a ty -> 'a ty
val qfbcornacchia : 'a ty -> 'a ty -> 'a ty
val qfbpow : 'a ty -> 'a ty -> 'a ty
val qfbpow_i : 'a ty -> 'a ty -> 'a ty
val qfbpowraw : 'a ty -> Signed.long -> 'a ty
val qfbpows : 'a ty -> Signed.long -> 'a ty
val qfbred : 'a ty -> 'a ty
val qfbred_i : 'a ty -> 'a ty
val qfbred0 : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val qfbredsl2 : 'a ty -> 'a ty -> 'a ty
val qfbsolve : 'a ty -> 'a ty -> Signed.long -> 'a ty
val qfbsqr : 'a ty -> 'a ty
val qfbsqr_i : 'a ty -> 'a ty
val qfisolvep : 'a ty -> 'a ty -> 'a ty

val qfr3_comp :
  'a ty -> 'a ty -> qfr_data Ctypes.structure Ctypes_static.ptr -> 'a ty

val qfr3_compraw : 'a ty -> 'a ty -> 'a ty

val qfr3_pow :
  'a ty -> 'a ty -> qfr_data Ctypes.structure Ctypes_static.ptr -> 'a ty

val qfr3_red : 'a ty -> qfr_data Ctypes.structure Ctypes_static.ptr -> 'a ty
val qfr3_rho : 'a ty -> qfr_data Ctypes.structure Ctypes_static.ptr -> 'a ty
val qfr3_to_qfr : 'a ty -> 'a ty -> 'a ty

val qfr5_comp :
  'a ty -> 'a ty -> qfr_data Ctypes.structure Ctypes_static.ptr -> 'a ty

val qfr5_compraw : 'a ty -> 'a ty -> 'a ty
val qfr5_dist : 'a ty -> 'a ty -> Signed.long -> 'a ty

val qfr5_pow :
  'a ty -> 'a ty -> qfr_data Ctypes.structure Ctypes_static.ptr -> 'a ty

val qfr5_red : 'a ty -> qfr_data Ctypes.structure Ctypes_static.ptr -> 'a ty
val qfr5_rho : 'a ty -> qfr_data Ctypes.structure Ctypes_static.ptr -> 'a ty
val qfr5_to_qfr : 'a ty -> 'a ty -> 'a ty -> 'a ty

val qfr_data_init :
  'a ty -> Signed.long -> qfr_data Ctypes.structure Ctypes_static.ptr -> unit

val qfr_to_qfr5 : 'a ty -> Signed.long -> 'a ty
val qfrsolvep : 'a ty -> 'a ty -> 'a ty
val quadgen : 'a ty -> 'a ty
val quadgen0 : 'a ty -> Signed.long -> 'a ty
val quadpoly : 'a ty -> 'a ty
val quadpoly_i : 'a ty -> 'a ty
val quadpoly0 : 'a ty -> Signed.long -> 'a ty
val fl_2gener_pre : pari_ulong -> pari_ulong -> pari_ulong
val fl_2gener_pre_i : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val fl_log : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_log_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fl_order : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val fl_powers : pari_ulong -> Signed.long -> pari_ulong -> 'a ty

val fl_powers_pre :
  pari_ulong -> Signed.long -> pari_ulong -> pari_ulong -> 'a ty

val fl_powu : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_powu_pre :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_sqrt : pari_ulong -> pari_ulong -> pari_ulong
val fl_sqrt_pre : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_sqrt_pre_i :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_sqrtl : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_sqrtl_pre :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_sqrtn :
  pari_ulong ->
  Signed.long ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong

val fl_sqrtn_pre :
  pari_ulong ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong

val fp_2gener : 'a ty -> 'a ty
val fp_2gener_i : 'a ty -> 'a ty -> 'a ty
val fp_factored_order : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_ispower : 'a ty -> 'a ty -> 'a ty -> int
val fp_log : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_order : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_pow_init : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val fp_pow_table : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fp_pows : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fp_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val fp_sqrt : 'a ty -> 'a ty -> 'a ty
val fp_sqrt_i : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_sqrtn : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val fpv_prod : 'a ty -> 'a ty -> 'a ty
val z_zv_mod : 'a ty -> 'a ty -> 'a ty
val z_zv_mod_tree : 'a ty -> 'a ty -> 'a ty -> 'a ty
val z_chinese : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val z_chinese_all :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val z_chinese_coprime : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val z_chinese_post : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val z_chinese_pre :
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  unit

val z_factor_listp : 'a ty -> 'a ty -> 'a ty
val z_nv_mod : 'a ty -> 'a ty -> 'a ty
val zm_nv_mod_tree : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zv_allpnqn : 'a ty -> 'a ty
val zv_chinese : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zv_chinese_tree : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val zv_chinesetree : 'a ty -> 'a ty -> 'a ty
val zv_nv_mod_tree : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zv_producttree : 'a ty -> 'a ty
val zx_nv_mod_tree : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zxc_nv_mod_tree : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zxm_nv_mod_tree : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zxx_nv_mod_tree : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zideallog : 'a ty -> 'a ty -> 'a ty
val bestappr : 'a ty -> 'a ty -> 'a ty
val bestapprpade : 'a ty -> Signed.long -> 'a ty
val chinese : 'a ty -> 'a ty -> 'a ty
val chinese1 : 'a ty -> 'a ty
val chinese1_coprime_z : 'a ty -> 'a ty
val contfrac0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val contfracpnqn : 'a ty -> Signed.long -> 'a ty
val fibo : Signed.long -> 'a ty
val gboundcf : 'a ty -> Signed.long -> 'a ty
val gcf : 'a ty -> 'a ty
val gcf2 : 'a ty -> 'a ty -> 'a ty

val get_fp_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  bb_field Ctypes.structure Ctypes_static.ptr

val hilbert : 'a ty -> 'a ty -> 'a ty -> Signed.long
val hilbertii : 'a ty -> 'a ty -> 'a ty -> Signed.long
val istotient : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val krois : 'a ty -> Signed.long -> Signed.long
val kroiu : 'a ty -> pari_ulong -> Signed.long
val kronecker : 'a ty -> 'a ty -> Signed.long
val krosi : Signed.long -> 'a ty -> Signed.long
val kross : Signed.long -> Signed.long -> Signed.long
val kroui : pari_ulong -> 'a ty -> Signed.long
val krouu : pari_ulong -> pari_ulong -> Signed.long
val lcmii : 'a ty -> 'a ty -> 'a ty
val fp_invgen : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val logint0 : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val logintall : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val mpfact : Signed.long -> 'a ty
val factorial_fl : Signed.long -> pari_ulong -> pari_ulong
val factorial_fp : Signed.long -> 'a ty -> 'a ty
val muls_interval : Signed.long -> Signed.long -> 'a ty
val mulu_interval : pari_ulong -> pari_ulong -> 'a ty
val mulu_interval_step : pari_ulong -> pari_ulong -> pari_ulong -> 'a ty
val ncv_chinese_center : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val ncv_chinese_center_tree : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val nmv_chinese_center : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val nmv_chinese_center_tree : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val nonsquare_fl : pari_ulong -> pari_ulong
val nxcv_chinese_center : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val nxcv_chinese_center_tree : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val nxmv_chinese_center : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val nxv_chinese_center : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val nxv_chinese_center_tree : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val zv_chinese_center : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val odd_prime_divisors : 'a ty -> 'a ty
val pgener_fl : pari_ulong -> pari_ulong
val pgener_fl_local : pari_ulong -> 'a ty -> pari_ulong
val pgener_fp : 'a ty -> 'a ty
val pgener_fp_local : 'a ty -> 'a ty -> 'a ty
val pgener_zl : pari_ulong -> pari_ulong
val pgener_zp : 'a ty -> 'a ty
val pnqn : 'a ty -> 'a ty
val ramanujantau : 'a ty -> Signed.long -> 'a ty
val rootsof1_fl : pari_ulong -> pari_ulong -> pari_ulong
val rootsof1_fp : 'a ty -> 'a ty -> 'a ty
val rootsof1u_fp : pari_ulong -> 'a ty -> 'a ty

val u_chinese_coprime :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val znlog : 'a ty -> 'a ty -> 'a ty -> 'a ty
val znorder : 'a ty -> 'a ty -> 'a ty
val znprimroot : 'a ty -> 'a ty
val znstar : 'a ty -> 'a ty
val znstar0 : 'a ty -> Signed.long -> 'a ty
val rgv_is_zvpos : 'a ty -> int
val rgv_is_zvnon0 : 'a ty -> int
val rgv_is_prv : 'a ty -> int
val z_issquarefree_fact : 'a ty -> Signed.long

val z_lsmoothen :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val z_smoothen :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val bigomega : 'a ty -> Signed.long
val bigomegau : pari_ulong -> Signed.long
val boundfact : 'a ty -> pari_ulong -> 'a ty
val check_arith_pos : 'a ty -> string -> 'a ty
val check_arith_non0 : 'a ty -> string -> 'a ty
val check_arith_all : 'a ty -> string -> 'a ty
val clean_z_factor : 'a ty -> 'a ty
val core : 'a ty -> 'a ty

val coredisc2_fact :
  'a ty ->
  Signed.long ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val coredisc2u_fact :
  'a ty ->
  Signed.long ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  pari_ulong

val corepartial : 'a ty -> Signed.long -> 'a ty
val core0 : 'a ty -> Signed.long -> 'a ty
val core2 : 'a ty -> 'a ty
val core2partial : 'a ty -> Signed.long -> 'a ty
val coredisc : 'a ty -> 'a ty
val coredisc0 : 'a ty -> Signed.long -> 'a ty
val coredisc2 : 'a ty -> 'a ty
val corediscs : Signed.long -> pari_ulong Ctypes_static.ptr -> Signed.long
val divisors : 'a ty -> 'a ty
val divisors_factored : 'a ty -> 'a ty
val divisors0 : 'a ty -> Signed.long -> 'a ty
val divisorsu : pari_ulong -> 'a ty
val divisorsu_moebius : 'a ty -> 'a ty
val divisorsu_fact : 'a ty -> 'a ty
val divisorsu_fact_factored : 'a ty -> 'a ty
val eulerphi : 'a ty -> 'a ty
val eulerphiu : pari_ulong -> pari_ulong
val eulerphiu_fact : 'a ty -> pari_ulong
val factor_pn_1 : 'a ty -> pari_ulong -> 'a ty
val factor_pn_1_limit : 'a ty -> Signed.long -> pari_ulong -> 'a ty
val factoru_pow : pari_ulong -> 'a ty
val fuse_z_factor : 'a ty -> 'a ty -> 'a ty
val is_z_factor : 'a ty -> int
val is_z_factornon0 : 'a ty -> int
val is_z_factorpos : 'a ty -> int
val is_nf_factor : 'a ty -> int
val is_nf_extfactor : 'a ty -> int
val issquarefree : 'a ty -> Signed.long
val numdiv : 'a ty -> 'a ty
val numdivu : Signed.long -> Signed.long
val numdivu_fact : 'a ty -> Signed.long
val omega : 'a ty -> Signed.long
val omegau : pari_ulong -> Signed.long
val sumdiv : 'a ty -> 'a ty
val sumdivk : 'a ty -> Signed.long -> 'a ty
val uissquarefree : pari_ulong -> Signed.long
val uissquarefree_fact : 'a ty -> Signed.long
val usumdiv_fact : 'a ty -> 'a ty
val usumdivk_fact : 'a ty -> pari_ulong -> 'a ty
val fpx_fpc_nfpoleval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val embed_t2 : 'a ty -> Signed.long -> 'a ty
val embednorm_t2 : 'a ty -> Signed.long -> 'a ty
val embed_norm : 'a ty -> Signed.long -> 'a ty
val check_zkmodule_i : 'a ty -> int
val check_zkmodule : 'a ty -> string -> unit
val checkbid : 'a ty -> unit
val checkbid_i : 'a ty -> 'a ty
val checkbnf : 'a ty -> 'a ty
val checkbnf_i : 'a ty -> 'a ty
val checkbnr : 'a ty -> unit
val checkbnr_i : 'a ty -> 'a ty
val checkabgrp : 'a ty -> unit
val checksqmat : 'a ty -> Signed.long -> unit
val checknf : 'a ty -> 'a ty
val checknf_i : 'a ty -> 'a ty
val checknfelt_mod : 'a ty -> 'a ty -> string -> 'a ty
val checkprid : 'a ty -> unit
val checkprid_i : 'a ty -> int
val checkrnf : 'a ty -> unit
val checkrnf_i : 'a ty -> int
val factoredpolred : 'a ty -> 'a ty -> 'a ty
val factoredpolred2 : 'a ty -> 'a ty -> 'a ty
val galoisapply : 'a ty -> 'a ty -> 'a ty -> 'a ty
val get_bnf : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty

val get_bnfpol :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val get_nf : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val get_nfpol : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val get_prid : 'a ty -> 'a ty
val idealfrobenius : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealfrobenius_aut : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealramfrobenius : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealramfrobenius_aut : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealramgroups : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealramgroups_aut : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val rnfpolredabs : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfpolredbest : 'a ty -> 'a ty -> Signed.long -> 'a ty
val smallpolred : 'a ty -> 'a ty
val smallpolred2 : 'a ty -> 'a ty
val tschirnhaus : 'a ty -> 'a ty
val zx_q_mul : 'a ty -> 'a ty -> 'a ty
val zx_q_normalize : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zx_z_normalize : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zx_to_monic : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zx_primitive_to_monic : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zxx_q_mul : 'a ty -> 'a ty -> 'a ty
val fq_to_nf : 'a ty -> 'a ty -> 'a ty
val fqm_to_nfm : 'a ty -> 'a ty -> 'a ty
val fqv_to_nfv : 'a ty -> 'a ty -> 'a ty
val fqx_to_nfx : 'a ty -> 'a ty -> 'a ty
val rg_nffix : string -> 'a ty -> 'a ty -> int -> 'a ty
val rgv_nffix : string -> 'a ty -> 'a ty -> int -> 'a ty
val rgx_nffix : string -> 'a ty -> 'a ty -> int -> 'a ty
val zx_composedsum : 'a ty -> 'a ty -> 'a ty
val zx_compositum : 'a ty -> 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val zpx_disc_val : 'a ty -> 'a ty -> Signed.long
val zpx_gcd : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val zpx_monic_factor : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpx_primedec : 'a ty -> 'a ty -> 'a ty
val zpx_reduced_resultant : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val zpx_reduced_resultant_fast : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zpx_resultant_val : 'a ty -> 'a ty -> 'a ty -> Signed.long -> Signed.long
val checkmodpr : 'a ty -> unit
val compositum : 'a ty -> 'a ty -> 'a ty
val compositum2 : 'a ty -> 'a ty -> 'a ty
val get_modpr : 'a ty -> 'a ty
val indexpartial : 'a ty -> 'a ty -> 'a ty
val modpr_genfq : 'a ty -> 'a ty
val idealprimedec : 'a ty -> 'a ty -> 'a ty
val idealprimedec_galois : 'a ty -> 'a ty -> 'a ty
val idealprimedec_degrees : 'a ty -> 'a ty -> 'a ty
val idealprimedec_kummer : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val idealprimedec_limit_f : 'a ty -> 'a ty -> Signed.long -> 'a ty
val idealprimedec_limit_norm : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rnfbasis : 'a ty -> 'a ty -> 'a ty
val rnfdedekind : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfdet : 'a ty -> 'a ty -> 'a ty
val rnfdisc_factored : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val rnfdiscf : 'a ty -> 'a ty -> 'a ty
val rnfequation : 'a ty -> 'a ty -> 'a ty
val rnfequation0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfequation2 : 'a ty -> 'a ty -> 'a ty

val rnfequationall :
  'a ty ->
  'a ty ->
  Signed.long Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val rnfhnfbasis : 'a ty -> 'a ty -> 'a ty
val rnfisfree : 'a ty -> 'a ty -> Signed.long
val rnflllgram : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfpolred : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfpseudobasis : 'a ty -> 'a ty -> 'a ty
val rnfsimplifybasis : 'a ty -> 'a ty -> 'a ty
val rnfsteinitz : 'a ty -> 'a ty -> 'a ty
val factorial_lval : pari_ulong -> pari_ulong -> Signed.long

val zk_to_fq_init :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val zk_to_fq : 'a ty -> 'a ty -> 'a ty
val qxqv_to_fpm : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zkmodprinit : 'a ty -> 'a ty -> 'a ty
val idealstar : 'a ty -> 'a ty -> Signed.long -> 'a ty
val idealstarprk : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgc_to_nfc : 'a ty -> 'a ty -> 'a ty
val rgm_rgx_mul : 'a ty -> 'a ty -> 'a ty
val rgm_to_nfm : 'a ty -> 'a ty -> 'a ty
val rgx_to_nfx : 'a ty -> 'a ty -> 'a ty
val zc_nfval : 'a ty -> 'a ty -> Signed.long
val zc_nfvalrem : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val zc_prdvd : 'a ty -> 'a ty -> int
val zm_zx_mul : 'a ty -> 'a ty -> 'a ty
val zv_snf_gcd : 'a ty -> 'a ty -> 'a ty
val algtobasis : 'a ty -> 'a ty -> 'a ty
val basistoalg : 'a ty -> 'a ty -> 'a ty
val ei_multable : 'a ty -> Signed.long -> 'a ty

val get_nf_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  bb_field Ctypes.structure Ctypes_static.ptr

val famat_nfvalrem : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val gpnfvalrem : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val idealfactorback : 'a ty -> 'a ty -> 'a ty -> int -> 'a ty
val ideallist : 'a ty -> Signed.long -> 'a ty
val ideallist0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val gideallist : 'a ty -> 'a ty -> Signed.long -> 'a ty
val ideallistarch : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ideallog : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ideallogmod : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val ideallog_units : 'a ty -> 'a ty -> 'a ty
val ideallog_units0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealprincipalunits : 'a ty -> 'a ty -> Signed.long -> 'a ty
val idealstar0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val idealstarmod : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val indices_to_vec01 : 'a ty -> Signed.long -> 'a ty
val matalgtobasis : 'a ty -> 'a ty -> 'a ty
val matbasistoalg : 'a ty -> 'a ty -> 'a ty
val multable : 'a ty -> 'a ty -> 'a ty
val pr_basis_perm : 'a ty -> 'a ty -> 'a ty
val pr_equal : 'a ty -> 'a ty -> int
val rnfalgtobasis : 'a ty -> 'a ty -> 'a ty
val rnfbasistoalg : 'a ty -> 'a ty -> 'a ty
val rnfeltnorm : 'a ty -> 'a ty -> 'a ty
val rnfelttrace : 'a ty -> 'a ty -> 'a ty
val set_sign_mod_divisor : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val tablemul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val tablemul_ei : 'a ty -> 'a ty -> Signed.long -> 'a ty
val tablemul_ei_ej : 'a ty -> Signed.long -> Signed.long -> 'a ty
val tablemulvec : 'a ty -> 'a ty -> 'a ty -> 'a ty
val tablesqr : 'a ty -> 'a ty -> 'a ty
val vec01_to_indices : 'a ty -> 'a ty
val vecsmall01_to_indices : 'a ty -> 'a ty
val zk_inv : 'a ty -> 'a ty -> 'a ty
val zk_multable : 'a ty -> 'a ty -> 'a ty
val zk_scalar_or_multable : 'a ty -> 'a ty -> 'a ty
val zkchinese : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zkchinese1 : 'a ty -> 'a ty -> 'a ty
val zkchineseinit : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val zkmultable_capz : 'a ty -> 'a ty
val zkmultable_inv : 'a ty -> 'a ty

val fl_invgen :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val rm_round_maxrank : 'a ty -> 'a ty
val zm_famat_limit : 'a ty -> 'a ty -> 'a ty
val zv_cba : 'a ty -> 'a ty
val zv_cba_extend : 'a ty -> 'a ty -> 'a ty
val z_cba : 'a ty -> 'a ty -> 'a ty
val z_ppgle : 'a ty -> 'a ty -> 'a ty
val z_ppio : 'a ty -> 'a ty -> 'a ty
val z_ppo : 'a ty -> 'a ty -> 'a ty
val famatv_factorback : 'a ty -> 'a ty -> 'a ty
val famatv_zv_factorback : 'a ty -> 'a ty -> 'a ty
val famat_z_gcd : 'a ty -> 'a ty -> 'a ty
val famat_div : 'a ty -> 'a ty -> 'a ty
val famat_div_shallow : 'a ty -> 'a ty -> 'a ty
val famat_idealfactor : 'a ty -> 'a ty -> 'a ty
val famat_inv : 'a ty -> 'a ty
val famat_inv_shallow : 'a ty -> 'a ty

val famat_makecoprime :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val famat_mul : 'a ty -> 'a ty -> 'a ty
val famat_mul_shallow : 'a ty -> 'a ty -> 'a ty
val famat_mulpow_shallow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val famat_mulpows_shallow : 'a ty -> 'a ty -> Signed.long -> 'a ty
val famat_pow : 'a ty -> 'a ty -> 'a ty
val famat_pow_shallow : 'a ty -> 'a ty -> 'a ty
val famat_pows_shallow : 'a ty -> Signed.long -> 'a ty
val famat_reduce : 'a ty -> 'a ty
val famat_remove_trivial : 'a ty -> 'a ty
val famat_sqr : 'a ty -> 'a ty
val famat_to_nf : 'a ty -> 'a ty -> 'a ty
val famat_to_nf_moddivisor : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val famat_to_nf_modideal_coprime :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val famatsmall_reduce : 'a ty -> 'a ty
val gpidealfactor : 'a ty -> 'a ty -> 'a ty -> 'a ty
val gpidealval : 'a ty -> 'a ty -> 'a ty -> 'a ty

val idealhnf_z_factor :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val idealhnf_z_factor_i :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val idealhnf_inv : 'a ty -> 'a ty -> 'a ty
val idealhnf_inv_z : 'a ty -> 'a ty -> 'a ty
val idealhnf_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealadd : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealaddmultoone : 'a ty -> 'a ty -> 'a ty
val idealaddtoone : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealaddtoone0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealaddtoone_i : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealaddtoone_raw : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealappr : 'a ty -> 'a ty -> 'a ty
val idealappr0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val idealapprfact : 'a ty -> 'a ty -> 'a ty
val idealchinese : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealcoprime : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealcoprimefact : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealdiv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealdiv0 : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val idealdivexact : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealdivpowprime : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealdown : 'a ty -> 'a ty -> 'a ty
val idealfactor : 'a ty -> 'a ty -> 'a ty
val idealfactor_limit : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val idealfactor_partial : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealhnf : 'a ty -> 'a ty -> 'a ty
val idealhnf0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealhnf_principal : 'a ty -> 'a ty -> 'a ty
val idealhnf_shallow : 'a ty -> 'a ty -> 'a ty
val idealhnf_two : 'a ty -> 'a ty -> 'a ty
val idealintersect : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealinv : 'a ty -> 'a ty -> 'a ty
val idealismaximal : 'a ty -> 'a ty -> 'a ty

val idealispower :
  'a ty -> 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr -> Signed.long

val idealmin : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealmul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealmul0 : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val idealmulpowprime : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealmulred : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealnumden : 'a ty -> 'a ty -> 'a ty
val idealpow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealpow0 : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val idealpowred : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealpows : 'a ty -> 'a ty -> Signed.long -> 'a ty
val idealprod : 'a ty -> 'a ty -> 'a ty
val idealprodprime : 'a ty -> 'a ty -> 'a ty
val idealprodval : 'a ty -> 'a ty -> 'a ty -> Signed.long
val idealpseudomin : 'a ty -> 'a ty -> 'a ty
val idealpseudomin_nonscalar : 'a ty -> 'a ty -> 'a ty
val idealpseudominvec : 'a ty -> 'a ty -> 'a ty
val idealpseudored : 'a ty -> 'a ty -> 'a ty
val idealred0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealred_elt : 'a ty -> 'a ty -> 'a ty
val idealredmodpower : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val idealsqr : 'a ty -> 'a ty -> 'a ty
val idealtwoelt : 'a ty -> 'a ty -> 'a ty
val idealtwoelt0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealtwoelt2 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val idealtyp : 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> Signed.long
val idealval : 'a ty -> 'a ty -> 'a ty -> Signed.long
val isideal : 'a ty -> 'a ty -> Signed.long
val matreduce : 'a ty -> 'a ty
val prv_lcm_capz : 'a ty -> 'a ty
val prv_primes : 'a ty -> 'a ty
val pr_hnf : 'a ty -> 'a ty -> 'a ty
val pr_inv : 'a ty -> 'a ty
val pr_inv_p : 'a ty -> 'a ty
val pr_uniformizer : 'a ty -> 'a ty -> 'a ty
val sunits_makecoprime : 'a ty -> 'a ty -> 'a ty -> 'a ty
val to_famat : 'a ty -> 'a ty -> 'a ty
val to_famat_shallow : 'a ty -> 'a ty -> 'a ty
val u_ppo : pari_ulong -> pari_ulong -> pari_ulong
val vecdiv : 'a ty -> 'a ty -> 'a ty
val vecinv : 'a ty -> 'a ty
val vecmul : 'a ty -> 'a ty -> 'a ty
val vecpow : 'a ty -> 'a ty -> 'a ty
val vecsqr : 'a ty -> 'a ty
val zkc_multable_mul : 'a ty -> 'a ty -> 'a ty
val eltreltoabs : 'a ty -> 'a ty -> 'a ty
val eltabstorel : 'a ty -> 'a ty -> 'a ty
val eltabstorel_lift : 'a ty -> 'a ty -> 'a ty
val rnf_build_nfabs : 'a ty -> Signed.long -> 'a ty
val rnf_zkabs : 'a ty -> 'a ty
val rnfcomplete : 'a ty -> unit
val rnfeltabstorel : 'a ty -> 'a ty -> 'a ty
val rnfeltdown : 'a ty -> 'a ty -> 'a ty
val rnfeltdown0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfeltreltoabs : 'a ty -> 'a ty -> 'a ty
val rnfeltup : 'a ty -> 'a ty -> 'a ty
val rnfeltup0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfidealabstorel : 'a ty -> 'a ty -> 'a ty
val rnfidealdown : 'a ty -> 'a ty -> 'a ty
val rnfidealfactor : 'a ty -> 'a ty -> 'a ty
val rnfidealhnf : 'a ty -> 'a ty -> 'a ty
val rnfidealmul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rnfidealnormabs : 'a ty -> 'a ty -> 'a ty
val rnfidealnormrel : 'a ty -> 'a ty -> 'a ty
val rnfidealprimedec : 'a ty -> 'a ty -> 'a ty
val rnfidealreltoabs : 'a ty -> 'a ty -> 'a ty
val rnfidealreltoabs0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfidealtwoelement : 'a ty -> 'a ty -> 'a ty
val rnfidealup : 'a ty -> 'a ty -> 'a ty
val rnfidealup0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfinit : 'a ty -> 'a ty -> 'a ty
val rnfinit0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val get_arith_zzm : 'a ty -> 'a ty
val get_arith_z : 'a ty -> 'a ty

val gen_ph_log :
  'a ty ->
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_shanks_init :
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_shanks :
  'a ty ->
  'a ty ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_shanks_sqrtn :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_gener :
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_ellgens :
  'a ty ->
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_ellgroup :
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_factored_order :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_order :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_select_order :
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_plog :
  'a ty ->
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  'a ty

val gen_pow :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_pow_fold :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty

val gen_pow_fold_i :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty

val gen_pow_i :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_pow_init :
  'a ty ->
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_pow_table :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_powers :
  'a ty ->
  Signed.long ->
  int ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty) Ctypes_static.static_funptr ->
  'a ty

val gen_powu :
  'a ty ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_powu_fold :
  'a ty ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty

val gen_powu_fold_i :
  'a ty ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty

val gen_powu_i :
  'a ty ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_product :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val matdetmod : 'a ty -> 'a ty -> 'a ty
val matimagemod : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val matinvmod : 'a ty -> 'a ty -> 'a ty
val matkermod : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val matsolvemod : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val bernfrac : Signed.long -> 'a ty
val bernpol : Signed.long -> Signed.long -> 'a ty
val bernreal : Signed.long -> Signed.long -> 'a ty
val bernvec : Signed.long -> 'a ty
val constbern : Signed.long -> unit
val eulerfrac : Signed.long -> 'a ty
val eulerpol : Signed.long -> Signed.long -> 'a ty
val eulerreal : Signed.long -> Signed.long -> 'a ty
val eulervec : Signed.long -> 'a ty
val harmonic : pari_ulong -> 'a ty
val harmonic0 : pari_ulong -> 'a ty -> 'a ty

val qr_init :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long ->
  int

val r_from_qr : 'a ty -> Signed.long -> 'a ty
val rgm_babai : 'a ty -> 'a ty -> 'a ty

val rgm_qr_init :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long ->
  int

val rgm_gram_schmidt : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val algdep : 'a ty -> Signed.long -> 'a ty
val algdep0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val bestapprnf : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val forqfvec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> float -> Signed.long)
  Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  unit

val forqfvec0 : 'a ty -> 'a ty -> 'a ty -> unit

val forqfvec1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  unit

val gaussred_from_qr : 'a ty -> Signed.long -> 'a ty
val lindep : 'a ty -> 'a ty
val lindep_xadic : 'a ty -> 'a ty
val lindep_bit : 'a ty -> Signed.long -> 'a ty
val lindep_padic : 'a ty -> 'a ty
val lindep0 : 'a ty -> Signed.long -> 'a ty
val lindep2 : 'a ty -> Signed.long -> 'a ty
val lindepfull_bit : 'a ty -> Signed.long -> 'a ty
val mathouseholder : 'a ty -> 'a ty -> 'a ty
val matqr : 'a ty -> Signed.long -> Signed.long -> 'a ty
val minim : 'a ty -> 'a ty -> 'a ty -> 'a ty
val minim_raw : 'a ty -> 'a ty -> 'a ty -> 'a ty
val minim_zm : 'a ty -> 'a ty -> 'a ty -> 'a ty
val minim2 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qfminim0 : 'a ty -> 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val qfperfection : 'a ty -> 'a ty
val qfrep0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val seralgdep : 'a ty -> Signed.long -> Signed.long -> 'a ty
val serdiffdep : 'a ty -> Signed.long -> Signed.long -> 'a ty
val vandermondeinverse : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val vandermondeinverseinit : 'a ty -> 'a ty
val zncoppersmith : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxq_reverse : 'a ty -> 'a ty -> 'a ty
val vec_equiv : 'a ty -> 'a ty
val rgv_polint : 'a ty -> 'a ty -> Signed.long -> 'a ty
val vec_reduce : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val rgxq_reverse : 'a ty -> 'a ty -> 'a ty
val zc_union_shallow : 'a ty -> 'a ty -> 'a ty
val zv_indexsort : 'a ty -> 'a ty
val zv_sort : 'a ty -> 'a ty
val zv_sort_inplace : 'a ty -> unit
val zv_sort_shallow : 'a ty -> 'a ty
val zv_sort_uniq : 'a ty -> 'a ty
val zv_sort_uniq_shallow : 'a ty -> 'a ty
val zv_union_shallow : 'a ty -> 'a ty -> 'a ty
val binomial : 'a ty -> Signed.long -> 'a ty
val binomial0 : 'a ty -> 'a ty -> 'a ty
val binomialuu : pari_ulong -> pari_ulong -> 'a ty
val cmp_flx : 'a ty -> 'a ty -> int
val cmp_rgx : 'a ty -> 'a ty -> int
val cmp_nodata : unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int
val cmp_prime_ideal : 'a ty -> 'a ty -> int
val cmp_prime_over_p : 'a ty -> 'a ty -> int
val cmp_universal : 'a ty -> 'a ty -> int
val convol : 'a ty -> 'a ty -> 'a ty
val gen_cmp_rgx : unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int
val dirdiv : 'a ty -> 'a ty -> 'a ty
val dirmul : 'a ty -> 'a ty -> 'a ty
val eulerianpol : Signed.long -> Signed.long -> 'a ty
val gprec_wensure : 'a ty -> Signed.long -> 'a ty

val gen_indexsort :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  'a ty

val gen_indexsort_uniq :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  'a ty

val gen_search :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  Signed.long

val gen_setminus :
  'a ty -> 'a ty -> ('a ty -> 'a ty -> int) Ctypes_static.static_funptr -> 'a ty

val gen_sort :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  'a ty

val gen_sort_inplace :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  'a ty Ctypes_static.ptr ->
  unit

val gen_sort_shallow :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  'a ty

val gen_sort_uniq :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  'a ty

val getstack : unit -> Signed.long
val gettime : unit -> Signed.long
val getabstime : unit -> Signed.long
val getwalltime : unit -> 'a ty
val gprec : 'a ty -> Signed.long -> 'a ty
val gprec_wtrunc : 'a ty -> Signed.long -> 'a ty
val gprec_w : 'a ty -> Signed.long -> 'a ty
val indexlexsort : 'a ty -> 'a ty
val indexsort : 'a ty -> 'a ty
val indexvecsort : 'a ty -> 'a ty -> 'a ty
val laplace : 'a ty -> 'a ty
val lexsort : 'a ty -> 'a ty
val mathilbert : Signed.long -> 'a ty
val matqpascal : Signed.long -> 'a ty -> 'a ty

val merge_factor :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  'a ty

val merge_sort_uniq :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  'a ty

val modreverse : 'a ty -> 'a ty
val setbinop : 'a ty -> 'a ty -> 'a ty -> 'a ty
val setdelta : 'a ty -> 'a ty -> 'a ty
val setintersect : 'a ty -> 'a ty -> 'a ty
val setisset : 'a ty -> Signed.long
val setminus : 'a ty -> 'a ty -> 'a ty
val setunion : 'a ty -> 'a ty -> 'a ty
val setunion_i : 'a ty -> 'a ty -> 'a ty
val sort : 'a ty -> 'a ty

val sort_factor :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  'a ty

val stirling : Signed.long -> Signed.long -> Signed.long -> 'a ty
val stirling1 : pari_ulong -> pari_ulong -> 'a ty
val stirling2 : pari_ulong -> pari_ulong -> 'a ty

val tablesearch :
  'a ty ->
  'a ty ->
  ('a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  Signed.long

val vecbinomial : Signed.long -> 'a ty
val vecsearch : 'a ty -> 'a ty -> 'a ty -> Signed.long
val vecsort : 'a ty -> 'a ty -> 'a ty
val vecsort0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zv_search : 'a ty -> Signed.long -> Signed.long
val bits_to_int : 'a ty -> Signed.long -> 'a ty
val bits_to_u : 'a ty -> Signed.long -> pari_ulong
val binaire : 'a ty -> 'a ty
val binary_2k : 'a ty -> Signed.long -> 'a ty
val binary_2k_nv : 'a ty -> Signed.long -> 'a ty
val binary_zv : 'a ty -> 'a ty
val bittest : 'a ty -> Signed.long -> Signed.long
val fromdigits_2k : 'a ty -> Signed.long -> 'a ty
val gbitand : 'a ty -> 'a ty -> 'a ty
val gbitneg : 'a ty -> Signed.long -> 'a ty
val gbitnegimply : 'a ty -> 'a ty -> 'a ty
val gbitor : 'a ty -> 'a ty -> 'a ty
val gbittest : 'a ty -> Signed.long -> 'a ty
val gbitxor : 'a ty -> 'a ty -> 'a ty
val hammingl : pari_ulong -> Signed.long
val hammingweight : 'a ty -> Signed.long
val ibitand : 'a ty -> 'a ty -> 'a ty
val ibitnegimply : 'a ty -> 'a ty -> 'a ty
val ibitor : 'a ty -> 'a ty -> 'a ty
val ibitxor : 'a ty -> 'a ty -> 'a ty
val nv_fromdigits_2k : 'a ty -> Signed.long -> 'a ty
val bnflogef : 'a ty -> 'a ty -> 'a ty
val bnflog : 'a ty -> 'a ty -> 'a ty
val bnflogdegree : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rnfislocalcyclo : 'a ty -> Signed.long
val bnfisunit : 'a ty -> 'a ty -> 'a ty
val bnfissunit : 'a ty -> 'a ty -> 'a ty -> 'a ty
val bnfsunit : 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnfunits : 'a ty -> 'a ty -> 'a ty
val bnfisunit0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val sunits_mod_units : 'a ty -> 'a ty -> 'a ty
val buchquad : 'a ty -> float -> float -> Signed.long -> 'a ty
val quadclassno : 'a ty -> 'a ty
val quadclassnos : Signed.long -> Signed.long
val quadclassunit0 : 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty
val buchall : 'a ty -> Signed.long -> Signed.long -> 'a ty

val buchall_param :
  'a ty -> float -> float -> Signed.long -> Signed.long -> Signed.long -> 'a ty

val bnf_build_cheapfu : 'a ty -> 'a ty
val bnf_build_cycgen : 'a ty -> 'a ty
val bnf_build_matalpha : 'a ty -> 'a ty
val bnf_build_units : 'a ty -> 'a ty
val bnf_compactfu : 'a ty -> 'a ty
val bnf_compactfu_mat : 'a ty -> 'a ty
val bnf_has_fu : 'a ty -> 'a ty
val bnfinit0 : 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty
val bnfisprincipal0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnfnewprec : 'a ty -> Signed.long -> 'a ty
val bnfnewprec_shallow : 'a ty -> Signed.long -> 'a ty
val bnftestprimes : 'a ty -> 'a ty -> unit
val bnrnewprec : 'a ty -> Signed.long -> 'a ty
val bnrnewprec_shallow : 'a ty -> Signed.long -> 'a ty
val isprincipalfact : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val isprincipalfact_or_fail : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val isprincipal : 'a ty -> 'a ty -> 'a ty
val signunits : 'a ty -> 'a ty
val hermite_bound : Signed.long -> Signed.long -> 'a ty

val bnr_subgroup_sanitize :
  'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val bnr_char_sanitize :
  'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val abc_to_bnr :
  'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> int -> 'a ty

val buchray : 'a ty -> 'a ty -> Signed.long -> 'a ty
val buchraymod : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val bnrautmatrix : 'a ty -> 'a ty -> 'a ty
val bnr_subgroup_check : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val bnrchar : 'a ty -> 'a ty -> 'a ty -> 'a ty
val bnrchar_primitive : 'a ty -> 'a ty -> 'a ty -> 'a ty
val bnrclassno : 'a ty -> 'a ty -> 'a ty
val bnrclassno0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val bnrclassnolist : 'a ty -> 'a ty -> 'a ty
val bnrchar_primitive_raw : 'a ty -> 'a ty -> 'a ty -> 'a ty
val bnrconductor_factored : 'a ty -> 'a ty -> 'a ty
val bnrconductor_raw : 'a ty -> 'a ty -> 'a ty
val bnrconductormod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val bnrconductor0 : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnrconductor : 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnrconductor_i : 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnrconductorofchar : 'a ty -> 'a ty -> 'a ty
val bnrdisc0 : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnrdisc : 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnrdisclist0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val bnrgaloismatrix : 'a ty -> 'a ty -> 'a ty
val bnrgaloisapply : 'a ty -> 'a ty -> 'a ty -> 'a ty
val bnrinit0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnrinitmod : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val bnrisconductor0 : 'a ty -> 'a ty -> 'a ty -> Signed.long
val bnrisconductor : 'a ty -> 'a ty -> Signed.long
val bnrisgalois : 'a ty -> 'a ty -> 'a ty -> Signed.long
val bnrisprincipalmod : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnrisprincipal : 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnrmap : 'a ty -> 'a ty -> 'a ty
val bnrsurjection : 'a ty -> 'a ty -> 'a ty
val bnfnarrow : 'a ty -> 'a ty
val bnfcertify : 'a ty -> Signed.long
val bnfcertify0 : 'a ty -> Signed.long -> Signed.long
val bnrcompositum : 'a ty -> 'a ty -> 'a ty
val decodemodule : 'a ty -> 'a ty -> 'a ty
val discrayabslist : 'a ty -> 'a ty -> 'a ty
val discrayabslistarch : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val idealmoddivisor : 'a ty -> 'a ty -> 'a ty
val isprincipalray : 'a ty -> 'a ty -> 'a ty
val isprincipalraygen : 'a ty -> 'a ty -> 'a ty
val rnfconductor : 'a ty -> 'a ty -> 'a ty
val rnfconductor0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfnormgroup : 'a ty -> 'a ty -> 'a ty
val subgrouplist0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val bnfisnorm : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfisnorm : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rnfisnorminit : 'a ty -> 'a ty -> int -> 'a ty
val coprimes_zv : pari_ulong -> 'a ty
val char_check : 'a ty -> 'a ty -> int
val charconj : 'a ty -> 'a ty -> 'a ty
val charconj0 : 'a ty -> 'a ty -> 'a ty
val chardiv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val chardiv0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val chareval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val chargalois : 'a ty -> 'a ty -> 'a ty
val charker : 'a ty -> 'a ty -> 'a ty
val charker0 : 'a ty -> 'a ty -> 'a ty
val charmul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val charmul0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val charorder : 'a ty -> 'a ty -> 'a ty
val charorder0 : 'a ty -> 'a ty -> 'a ty
val charpow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val charpow0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val char_denormalize : 'a ty -> 'a ty -> 'a ty -> 'a ty
val char_normalize : 'a ty -> 'a ty -> 'a ty
val char_simplify : 'a ty -> 'a ty -> 'a ty
val checkznstar_i : 'a ty -> int
val cyc_normalize : 'a ty -> 'a ty
val ncharvecexpo : 'a ty -> 'a ty -> 'a ty
val znchar : 'a ty -> 'a ty
val znchar_quad : 'a ty -> 'a ty -> 'a ty
val zncharcheck : 'a ty -> 'a ty -> int
val zncharconductor : 'a ty -> 'a ty -> 'a ty
val zncharconj : 'a ty -> 'a ty -> 'a ty
val znchardecompose : 'a ty -> 'a ty -> 'a ty -> 'a ty
val znchardiv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val znchareval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val zncharinduce : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zncharisodd : 'a ty -> 'a ty -> Signed.long
val zncharker : 'a ty -> 'a ty -> 'a ty
val zncharmul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zncharorder : 'a ty -> 'a ty -> 'a ty
val zncharpow : 'a ty -> 'a ty -> 'a ty -> 'a ty
val znchartokronecker : 'a ty -> 'a ty -> Signed.long -> 'a ty
val znchartoprimitive : 'a ty -> 'a ty -> 'a ty
val znconrey_check : 'a ty -> 'a ty -> int
val znconrey_normalized : 'a ty -> 'a ty -> 'a ty
val znconreychar : 'a ty -> 'a ty -> 'a ty
val znconreyfromchar_normalized : 'a ty -> 'a ty -> 'a ty
val znconreyconductor : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val znconreyexp : 'a ty -> 'a ty -> 'a ty
val znconreyfromchar : 'a ty -> 'a ty -> 'a ty
val znconreylog : 'a ty -> 'a ty -> 'a ty
val znconreylog_normalize : 'a ty -> 'a ty -> 'a ty
val znlog0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zv_cyc_minimal : 'a ty -> 'a ty -> 'a ty -> Signed.long
val zv_cyc_minimize : 'a ty -> 'a ty -> 'a ty -> Signed.long
val closure_deriv : 'a ty -> 'a ty
val closure_derivn : 'a ty -> Signed.long -> 'a ty

val localvars_find :
  'a ty -> entree Ctypes.structure Ctypes_static.ptr -> Signed.long

val localvars_read_str : string -> 'a ty -> 'a ty
val snm_closure : entree Ctypes.structure Ctypes_static.ptr -> 'a ty -> 'a ty
val strtoclosure : string -> Signed.long -> 'a ty
val strtofunction : string -> 'a ty
val gconcat : 'a ty -> 'a ty -> 'a ty
val gconcat1 : 'a ty -> 'a ty
val matconcat : 'a ty -> 'a ty
val shallowconcat : 'a ty -> 'a ty -> 'a ty
val shallowconcat1 : 'a ty -> 'a ty
val shallowmatconcat : 'a ty -> 'a ty
val vconcat : 'a ty -> 'a ty -> 'a ty
val default0 : string -> string -> 'a ty
val getrealprecision : unit -> Signed.long
val pari_is_default : string -> entree Ctypes.structure Ctypes_static.ptr
val sd_texstyle : string -> Signed.long -> 'a ty
val sd_colors : string -> Signed.long -> 'a ty
val sd_compatible : string -> Signed.long -> 'a ty
val sd_datadir : string -> Signed.long -> 'a ty
val sd_debug : string -> Signed.long -> 'a ty
val sd_debugfiles : string -> Signed.long -> 'a ty
val sd_debugmem : string -> Signed.long -> 'a ty
val sd_factor_add_primes : string -> Signed.long -> 'a ty
val sd_factor_proven : string -> Signed.long -> 'a ty
val sd_format : string -> Signed.long -> 'a ty
val sd_histsize : string -> Signed.long -> 'a ty
val sd_log : string -> Signed.long -> 'a ty
val sd_logfile : string -> Signed.long -> 'a ty
val sd_nbthreads : string -> Signed.long -> 'a ty
val sd_new_galois_format : string -> Signed.long -> 'a ty
val sd_output : string -> Signed.long -> 'a ty
val sd_parisize : string -> Signed.long -> 'a ty
val sd_parisizemax : string -> Signed.long -> 'a ty
val sd_path : string -> Signed.long -> 'a ty
val sd_plothsizes : string -> Signed.long -> 'a ty
val sd_prettyprinter : string -> Signed.long -> 'a ty
val sd_primelimit : string -> Signed.long -> 'a ty
val sd_realbitprecision : string -> Signed.long -> 'a ty
val sd_realprecision : string -> Signed.long -> 'a ty
val sd_secure : string -> Signed.long -> 'a ty
val sd_seriesprecision : string -> Signed.long -> 'a ty
val sd_simplify : string -> Signed.long -> 'a ty
val sd_sopath : string -> int -> 'a ty
val sd_strictargs : string -> Signed.long -> 'a ty
val sd_strictmatch : string -> Signed.long -> 'a ty

val sd_string :
  string -> Signed.long -> string -> string Ctypes_static.ptr -> 'a ty

val sd_threadsize : string -> Signed.long -> 'a ty
val sd_threadsizemax : string -> Signed.long -> 'a ty

val sd_intarray :
  string -> Signed.long -> 'a ty Ctypes_static.ptr -> string -> 'a ty

val sd_toggle :
  string -> Signed.long -> string -> int Ctypes_static.ptr -> 'a ty

val sd_ulong :
  string ->
  Signed.long ->
  string ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  pari_ulong ->
  string Ctypes_static.ptr ->
  'a ty

val setdefault : string -> string -> Signed.long -> 'a ty

val setrealprecision :
  Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val digits : 'a ty -> 'a ty -> 'a ty
val fromdigits : 'a ty -> 'a ty -> 'a ty
val fromdigitsu : 'a ty -> 'a ty -> 'a ty

val gen_digits :
  'a ty ->
  'a ty ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty

val gen_fromdigits :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  'a ty

val sumdigits : 'a ty -> 'a ty
val sumdigits0 : 'a ty -> 'a ty -> 'a ty
val sumdigitsu : pari_ulong -> pari_ulong
val ecpp : 'a ty -> 'a ty
val ecpp0 : 'a ty -> Signed.long -> 'a ty
val ecppexport : 'a ty -> Signed.long -> 'a ty
val ecppisvalid : 'a ty -> Signed.long
val isprimeecpp : 'a ty -> Signed.long
val sd_breakloop : string -> Signed.long -> 'a ty
val sd_echo : string -> Signed.long -> 'a ty
val sd_graphcolormap : string -> Signed.long -> 'a ty
val sd_graphcolors : string -> Signed.long -> 'a ty
val sd_help : string -> Signed.long -> 'a ty
val sd_histfile : string -> Signed.long -> 'a ty
val sd_lines : string -> Signed.long -> 'a ty
val sd_linewrap : string -> Signed.long -> 'a ty
val sd_prompt : string -> Signed.long -> 'a ty
val sd_prompt_cont : string -> Signed.long -> 'a ty
val sd_psfile : string -> Signed.long -> 'a ty
val sd_readline : string -> Signed.long -> 'a ty
val sd_recover : string -> Signed.long -> 'a ty
val sd_timer : string -> Signed.long -> 'a ty
val pari_hit_return : unit -> unit
val gp_load_gprc : unit -> unit
val gp_meta : string -> int -> int
val gphelp_keyword_list : unit -> string Ctypes_static.ptr
val pari_center : string -> unit
val pari_community : unit -> Signed.long
val pari_print_version : unit -> unit
val gp_format_time : Signed.long -> string
val gp_format_prompt : string -> string
val pari_alarm : Signed.long -> unit
val gp_alarm : Signed.long -> 'a ty -> 'a ty
val gp_input : unit -> 'a ty
val gp_allocatemem : 'a ty -> unit
val gp_handle_exception : Signed.long -> int
val gp_alarm_handler : int -> unit
val gp_sigint_fun : unit -> unit
val gp_help : string -> Signed.long -> unit
val gp_echo_and_log : string -> string -> unit
val print_fun_list : string Ctypes_static.ptr -> Signed.long -> unit
val strtime : Signed.long -> 'a ty

val direuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty

val dirpowers : Signed.long -> 'a ty -> Signed.long -> 'a ty
val dirpowerssum : pari_ulong -> 'a ty -> Signed.long -> Signed.long -> 'a ty

val dirpowerssumfun :
  pari_ulong ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> pari_ulong -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  'a ty

val vecpowuu : Signed.long -> pari_ulong -> 'a ty
val vecpowug : Signed.long -> 'a ty -> Signed.long -> 'a ty

val forell :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  unit

val akell : 'a ty -> 'a ty -> 'a ty
val bilhell : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val checkell : 'a ty -> unit
val checkell_fq : 'a ty -> unit
val checkell_q : 'a ty -> unit
val checkell_qp : 'a ty -> unit
val checkellisog : 'a ty -> unit
val checkellpt : 'a ty -> unit
val checkell5 : 'a ty -> unit
val cxredsl2 : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val cxredsl2_i :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val ec_2divpol_evalx : 'a ty -> 'a ty -> 'a ty
val ec_3divpol_evalx : 'a ty -> 'a ty -> 'a ty
val ec_bmodel : 'a ty -> Signed.long -> 'a ty
val ec_f_evalx : 'a ty -> 'a ty -> 'a ty
val ec_h_evalx : 'a ty -> 'a ty -> 'a ty
val ec_dfdx_evalq : 'a ty -> 'a ty -> 'a ty
val ec_dfdy_evalq : 'a ty -> 'a ty -> 'a ty
val ec_dmfdy_evalq : 'a ty -> 'a ty -> 'a ty
val ec_half_deriv_2divpol : 'a ty -> Signed.long -> 'a ty
val ec_half_deriv_2divpol_evalx : 'a ty -> 'a ty -> 'a ty
val ec_phi2 : 'a ty -> Signed.long -> 'a ty
val oncurve : 'a ty -> 'a ty -> int
val orderell : 'a ty -> 'a ty -> 'a ty
val pointell : 'a ty -> 'a ty -> Signed.long -> 'a ty
val point_to_a4a6 : 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val point_to_a4a6_fl :
  'a ty -> 'a ty -> pari_ulong -> pari_ulong Ctypes_static.ptr -> 'a ty

val zell : 'a ty -> 'a ty -> Signed.long -> 'a ty
val qp_agm2_sequence : 'a ty -> 'a ty -> 'a ty

val qp_ascending_landen :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val qp_descending_landen :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val hyperell_locally_soluble : 'a ty -> 'a ty -> Signed.long

val flxq_elldivpolmod :
  'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty -> pari_ulong -> 'a ty

val fp_ellcard_sea : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val fq_ellcard_sea :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty

val fq_elldivpolmod :
  'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val externstr : string -> 'a ty
val gp_filter : string -> string
val gpextern : string -> 'a ty
val gpsystem : string -> Signed.long
val readstr : string -> 'a ty
val gentogenstr_nospace : 'a ty -> 'a ty
val gentogenstr : 'a ty -> 'a ty
val gentotexstr : 'a ty -> string
val gentostr : 'a ty -> string
val gentostr_raw : 'a ty -> string
val gentostr_unquoted : 'a ty -> string
val str : 'a ty -> 'a ty
val strexpand : 'a ty -> 'a ty
val strtex : 'a ty -> 'a ty
val brute : 'a ty -> char -> Signed.long -> unit
val dbggen : 'a ty -> Signed.long -> unit
val error0 : 'a ty -> unit
val dbg_pari_heap : unit -> unit
val err_flush : unit -> unit
val err_printf : string -> unit
val gp_getenv : string -> 'a ty
val gp_fileclose : Signed.long -> unit
val gp_fileextern : string -> Signed.long
val gp_fileflush : Signed.long -> unit
val gp_fileflush0 : 'a ty -> unit
val gp_fileopen : string -> string -> Signed.long
val gp_fileread : Signed.long -> 'a ty
val gp_filereadstr : Signed.long -> 'a ty
val gp_filewrite : Signed.long -> string -> unit
val gp_filewrite1 : Signed.long -> string -> unit
val gp_read_file : string -> 'a ty
val gp_read_str_multiline : string -> string -> 'a ty
val gp_readvec_file : string -> 'a ty
val gpinstall : string -> string -> string -> string -> unit
val gsprintf : string -> 'a ty
val itostr : 'a ty -> string
val matbrute : 'a ty -> char -> Signed.long -> unit
val os_getenv : string -> string
val uordinal : pari_ulong -> string
val outmat : 'a ty -> unit
val output : 'a ty -> unit
val rgv_to_str : 'a ty -> Signed.long -> string
val pari_add_hist : 'a ty -> Signed.long -> Signed.long -> unit
val pari_ask_confirm : string -> unit
val pari_flush : unit -> unit
val pari_get_hist : Signed.long -> 'a ty
val pari_get_histrtime : Signed.long -> Signed.long
val pari_get_histtime : Signed.long -> Signed.long
val pari_get_homedir : string -> string
val pari_histtime : Signed.long -> 'a ty
val pari_is_dir : string -> int
val pari_is_file : string -> int
val pari_last_was_newline : unit -> int
val pari_set_last_newline : int -> unit
val pari_nb_hist : unit -> pari_ulong
val pari_printf : string -> unit
val pari_putc : char -> unit
val pari_puts : string -> unit
val pari_sprintf : string -> string
val pari_stdin_isatty : unit -> int
val pari_unique_dir : string -> string
val pari_unique_filename : string -> string
val pari_unique_filename_suffix : string -> string -> string
val pari_unlink : string -> unit
val path_expand : string -> string
val pari_sprint0 : string -> 'a ty -> Signed.long -> string
val print : 'a ty -> unit
val printp : 'a ty -> unit
val print1 : 'a ty -> unit
val printf0 : string -> 'a ty -> unit
val printsep : string -> 'a ty -> unit
val printsep1 : string -> 'a ty -> unit
val printtex : 'a ty -> unit
val stack_sprintf : string -> string
val str_init : pari_str Ctypes.structure Ctypes_static.ptr -> int -> unit
val str_printf : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val str_putc : pari_str Ctypes.structure Ctypes_static.ptr -> char -> unit
val str_puts : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val strftime_expand : string -> string -> Signed.long -> unit
val strprintf : string -> 'a ty -> 'a ty
val term_color : Signed.long -> unit
val term_get_color : string -> Signed.long -> string
val texe : 'a ty -> char -> Signed.long -> unit
val warning0 : 'a ty -> unit
val write0 : string -> 'a ty -> unit
val write1 : string -> 'a ty -> unit
val writebin : string -> 'a ty -> unit
val writetex : string -> 'a ty -> unit
val bincopy_relink : 'a ty -> 'a ty -> unit
val bitprecision0 : 'a ty -> Signed.long -> 'a ty
val bitprecision00 : 'a ty -> 'a ty -> 'a ty
val break0 : Signed.long -> 'a ty
val call0 : 'a ty -> 'a ty -> 'a ty
val closure_callgen0prec : 'a ty -> Signed.long -> 'a ty
val closure_callgen1 : 'a ty -> 'a ty -> 'a ty
val closure_callgen1prec : 'a ty -> 'a ty -> Signed.long -> 'a ty
val closure_callgen2 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val closure_callgenall : 'a ty -> Signed.long -> 'a ty
val closure_callgenvec : 'a ty -> 'a ty -> 'a ty
val closure_callgenvecdef : 'a ty -> 'a ty -> 'a ty -> 'a ty
val closure_callgenvecdefprec : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val closure_callgenvecprec : 'a ty -> 'a ty -> Signed.long -> 'a ty
val closure_callvoid1 : 'a ty -> 'a ty -> unit
val closure_context : Signed.long -> Signed.long -> Signed.long
val closure_disassemble : 'a ty -> unit
val closure_err : Signed.long -> unit
val closure_evalbrk : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val closure_evalgen : 'a ty -> 'a ty
val closure_evalnobrk : 'a ty -> 'a ty
val closure_evalres : 'a ty -> 'a ty
val closure_evalvoid : 'a ty -> unit
val closure_func_err : unit -> string
val closure_trapgen : 'a ty -> Signed.long -> 'a ty
val copybin_unlink : 'a ty -> 'a ty
val getlocalprec : Signed.long -> Signed.long
val getlocalbitprec : Signed.long -> Signed.long
val get_lex : Signed.long -> 'a ty
val get_localprec : unit -> Signed.long
val get_localbitprec : unit -> Signed.long
val gp_call : unit Ctypes_static.ptr -> 'a ty -> 'a ty
val gp_callprec : unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty
val gp_call2 : unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty
val gp_callbool : unit Ctypes_static.ptr -> 'a ty -> Signed.long
val gp_callvoid : unit Ctypes_static.ptr -> 'a ty -> Signed.long
val gp_eval : unit Ctypes_static.ptr -> 'a ty -> 'a ty
val gp_evalbool : unit Ctypes_static.ptr -> 'a ty -> Signed.long
val gp_evalprec : unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty
val gp_evalupto : unit Ctypes_static.ptr -> 'a ty -> 'a ty
val gp_evalvoid : unit Ctypes_static.ptr -> 'a ty -> Signed.long
val localprec : 'a ty -> unit
val localbitprec : 'a ty -> unit
val loop_break : unit -> Signed.long
val next0 : Signed.long -> 'a ty
val pareval : 'a ty -> 'a ty
val pari_self : unit -> 'a ty
val parsum : 'a ty -> 'a ty -> 'a ty -> 'a ty
val parvector : Signed.long -> 'a ty -> 'a ty
val pop_lex : Signed.long -> unit
val pop_localprec : unit -> unit
val precision0 : 'a ty -> Signed.long -> 'a ty
val precision00 : 'a ty -> 'a ty -> 'a ty
val push_lex : 'a ty -> 'a ty -> unit
val push_localbitprec : Signed.long -> unit
val push_localprec : Signed.long -> unit
val return0 : 'a ty -> 'a ty
val set_lex : Signed.long -> 'a ty -> unit

val forcomposite_init :
  forcomposite_t Ctypes.structure Ctypes_static.ptr -> 'a ty -> 'a ty -> int

val forcomposite_next :
  forcomposite_t Ctypes.structure Ctypes_static.ptr -> 'a ty

val forprime_next : forprime_t Ctypes.structure Ctypes_static.ptr -> 'a ty

val forprime_init :
  forprime_t Ctypes.structure Ctypes_static.ptr -> 'a ty -> 'a ty -> int

val forprimestep_init :
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  int

val initprimes :
  pari_ulong ->
  Signed.long Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  byteptr

val initprimetable : pari_ulong -> unit

val init_primepointer_geq :
  pari_ulong -> byteptr Ctypes_static.ptr -> pari_ulong

val init_primepointer_gt : pari_ulong -> byteptr Ctypes_static.ptr -> pari_ulong

val init_primepointer_leq :
  pari_ulong -> byteptr Ctypes_static.ptr -> pari_ulong

val init_primepointer_lt : pari_ulong -> byteptr Ctypes_static.ptr -> pari_ulong
val maxprime : unit -> pari_ulong
val maxprimen : unit -> pari_ulong
val maxprime_check : pari_ulong -> unit
val maxprimelim : unit -> pari_ulong
val pari_init_primes : pari_ulong -> unit
val prodprimes : unit -> 'a ty

val u_forprime_next :
  forprime_t Ctypes.structure Ctypes_static.ptr -> pari_ulong

val u_forprime_init :
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  pari_ulong ->
  pari_ulong ->
  int

val u_forprime_restrict :
  forprime_t Ctypes.structure Ctypes_static.ptr -> pari_ulong -> unit

val u_forprime_arith_init :
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  int

val ff_frobenius : 'a ty -> Signed.long -> 'a ty
val ff_z_z_muldiv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ff_q_add : 'a ty -> 'a ty -> 'a ty
val ff_z_add : 'a ty -> 'a ty -> 'a ty
val ff_add : 'a ty -> 'a ty -> 'a ty
val ff_charpoly : 'a ty -> 'a ty
val ff_conjvec : 'a ty -> 'a ty
val ff_div : 'a ty -> 'a ty -> 'a ty
val ff_ellcard : 'a ty -> 'a ty
val ff_ellcard_sea : 'a ty -> Signed.long -> 'a ty
val ff_ellgroup : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val ff_elllog : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val ff_ellmul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ff_ellorder : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ff_elltwist : 'a ty -> 'a ty
val ff_ellrandom : 'a ty -> 'a ty
val ff_elltatepairing : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val ff_ellweilpairing : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val ff_equal0 : 'a ty -> int
val ff_equal1 : 'a ty -> int
val ff_equalm1 : 'a ty -> int
val ff_f : 'a ty -> Signed.long
val ff_gen : 'a ty -> 'a ty
val ff_inv : 'a ty -> 'a ty
val ff_issquare : 'a ty -> Signed.long
val ff_issquareall : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val ff_ispower : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val ff_log : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ff_map : 'a ty -> 'a ty -> 'a ty
val ff_minpoly : 'a ty -> 'a ty
val ff_mod : 'a ty -> 'a ty
val ff_mul2n : 'a ty -> Signed.long -> 'a ty
val ff_neg : 'a ty -> 'a ty
val ff_neg_i : 'a ty -> 'a ty
val ff_norm : 'a ty -> 'a ty
val ff_order : 'a ty -> 'a ty -> 'a ty
val ff_p_i : 'a ty -> 'a ty
val ff_pow : 'a ty -> 'a ty -> 'a ty
val ff_q : 'a ty -> 'a ty
val ff_samefield : 'a ty -> 'a ty -> int
val ff_sqr : 'a ty -> 'a ty
val ff_sqrt : 'a ty -> 'a ty
val ff_sqrtn : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val ff_sub : 'a ty -> 'a ty -> 'a ty
val ff_to_f2xq : 'a ty -> 'a ty
val ff_to_f2xq_i : 'a ty -> 'a ty
val ff_to_flxq : 'a ty -> 'a ty
val ff_to_flxq_i : 'a ty -> 'a ty
val ff_to_fpxq : 'a ty -> 'a ty
val ff_trace : 'a ty -> 'a ty
val ff_var : 'a ty -> Signed.long
val ffm_ffc_invimage : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffm_ffc_gauss : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffm_ffc_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffm_deplin : 'a ty -> 'a ty -> 'a ty
val ffm_det : 'a ty -> 'a ty -> 'a ty
val ffm_gauss : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffm_image : 'a ty -> 'a ty -> 'a ty
val ffm_indexrank : 'a ty -> 'a ty -> 'a ty
val ffm_inv : 'a ty -> 'a ty -> 'a ty
val ffm_invimage : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffm_ker : 'a ty -> 'a ty -> 'a ty
val ffm_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffm_rank : 'a ty -> 'a ty -> Signed.long
val ffm_suppl : 'a ty -> 'a ty -> 'a ty
val ffx_add : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffx_ddf : 'a ty -> 'a ty -> 'a ty
val ffx_degfact : 'a ty -> 'a ty -> 'a ty
val ffx_disc : 'a ty -> 'a ty -> 'a ty

val ffx_extgcd :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val ffx_factor : 'a ty -> 'a ty -> 'a ty
val ffx_factor_squarefree : 'a ty -> 'a ty -> 'a ty
val ffx_gcd : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffx_halfgcd : 'a ty -> 'a ty -> 'a ty -> 'a ty

val ffx_halfgcd_all :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val ffx_ispower :
  'a ty -> Signed.long -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

val ffx_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffx_preimage : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffx_preimagerel : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffx_rem : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffx_resultant : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffx_roots : 'a ty -> 'a ty -> 'a ty
val ffx_sqr : 'a ty -> 'a ty -> 'a ty
val ffxq_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffxq_minpoly : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffxq_sqr : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_to_ffx : 'a ty -> 'a ty -> 'a ty
val fq_to_ff : 'a ty -> 'a ty -> 'a ty
val z_ff_div : 'a ty -> 'a ty -> 'a ty
val ffembed : 'a ty -> 'a ty -> 'a ty
val fffrobenius : 'a ty -> Signed.long -> 'a ty
val ffinvmap : 'a ty -> 'a ty
val fflog : 'a ty -> 'a ty -> 'a ty -> 'a ty
val ffmap : 'a ty -> 'a ty -> 'a ty
val ffmaprel : 'a ty -> 'a ty -> 'a ty
val ffcompomap : 'a ty -> 'a ty -> 'a ty
val fforder : 'a ty -> 'a ty -> 'a ty
val ffprimroot : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val ffrandom : 'a ty -> 'a ty
val rg_is_ff : 'a ty -> 'a ty Ctypes_static.ptr -> int
val rgc_is_ffc : 'a ty -> 'a ty Ctypes_static.ptr -> int
val rgm_is_ffm : 'a ty -> 'a ty Ctypes_static.ptr -> int
val p_to_ff : 'a ty -> Signed.long -> 'a ty
val tp_to_ff : 'a ty -> 'a ty -> 'a ty
val flx_factcyclo : pari_ulong -> pari_ulong -> pari_ulong -> 'a ty
val fpx_factcyclo : pari_ulong -> 'a ty -> pari_ulong -> 'a ty
val factormodcyclo : Signed.long -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val checkgal : 'a ty -> 'a ty
val checkgroup : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val checkgroupelts : 'a ty -> 'a ty
val embed_disc : 'a ty -> Signed.long -> Signed.long -> 'a ty
val embed_roots : 'a ty -> Signed.long -> 'a ty
val galois_group : 'a ty -> 'a ty
val galoisconj : 'a ty -> 'a ty -> 'a ty
val galoisconj0 : 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty
val galoisconjclasses : 'a ty -> 'a ty
val galoisexport : 'a ty -> Signed.long -> 'a ty
val galoisfixedfield : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val galoisidentify : 'a ty -> 'a ty
val galoisinit : 'a ty -> 'a ty -> 'a ty
val galoisisabelian : 'a ty -> Signed.long -> 'a ty
val galoisisnormal : 'a ty -> 'a ty -> Signed.long
val galoispermtopol : 'a ty -> 'a ty -> 'a ty
val galoissplittinginit : 'a ty -> 'a ty -> 'a ty
val galoissubgroups : 'a ty -> 'a ty
val galoissubfields : 'a ty -> Signed.long -> Signed.long -> 'a ty
val numberofconjugates : 'a ty -> Signed.long -> Signed.long
val galoisnbpol : Signed.long -> 'a ty
val galoisgetgroup : Signed.long -> Signed.long -> 'a ty
val galoisgetname : Signed.long -> Signed.long -> 'a ty
val galoisgetpol : Signed.long -> Signed.long -> Signed.long -> 'a ty
val conj_i : 'a ty -> 'a ty
val conjvec : 'a ty -> Signed.long -> 'a ty
val divrunextu : 'a ty -> pari_ulong -> 'a ty
val gadd : 'a ty -> 'a ty -> 'a ty
val gaddsg : Signed.long -> 'a ty -> 'a ty
val gconj : 'a ty -> 'a ty
val gdiv : 'a ty -> 'a ty -> 'a ty
val gdivgs : 'a ty -> Signed.long -> 'a ty
val gdivgu : 'a ty -> pari_ulong -> 'a ty
val gdivgunextu : 'a ty -> pari_ulong -> 'a ty
val ginv : 'a ty -> 'a ty
val gmul : 'a ty -> 'a ty -> 'a ty
val gmul2n : 'a ty -> Signed.long -> 'a ty
val gmulsg : Signed.long -> 'a ty -> 'a ty
val gmulug : pari_ulong -> 'a ty -> 'a ty
val gsqr : 'a ty -> 'a ty
val gsub : 'a ty -> 'a ty -> 'a ty
val gsubsg : Signed.long -> 'a ty -> 'a ty
val mulcxi : 'a ty -> 'a ty
val mulcxmi : 'a ty -> 'a ty
val mulcxpowis : 'a ty -> Signed.long -> 'a ty
val qdivii : 'a ty -> 'a ty -> 'a ty
val qdiviu : 'a ty -> pari_ulong -> 'a ty
val qdivis : 'a ty -> Signed.long -> 'a ty
val ser_normalize : 'a ty -> 'a ty

val gassoc_proto :
  ('a ty -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  'a ty

val map_proto_g : ('a ty -> 'a ty) Ctypes_static.static_funptr -> 'a ty -> 'a ty

val map_proto_lg :
  ('a ty -> Signed.long) Ctypes_static.static_funptr -> 'a ty -> 'a ty

val map_proto_lgl :
  ('a ty -> Signed.long -> Signed.long) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val q_lval : 'a ty -> pari_ulong -> Signed.long
val q_lvalrem : 'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> Signed.long
val q_pval : 'a ty -> 'a ty -> Signed.long
val q_pvalrem : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val rgx_val : 'a ty -> Signed.long
val rgx_valrem : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val rgx_valrem_inexact : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val rgxv_maxdegree : 'a ty -> Signed.long
val zv_z_dvd : 'a ty -> 'a ty -> int
val zv_pval : 'a ty -> 'a ty -> Signed.long
val zv_pvalrem : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val zv_lval : 'a ty -> pari_ulong -> Signed.long
val zv_lvalrem : 'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> Signed.long
val zx_lvalrem : 'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> Signed.long
val zx_pval : 'a ty -> 'a ty -> Signed.long
val zx_pvalrem : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

val z_lvalrem_stop :
  'a ty Ctypes_static.ptr -> pari_ulong -> int Ctypes_static.ptr -> Signed.long

val cgetp : 'a ty -> 'a ty
val cvstop2 : Signed.long -> 'a ty -> 'a ty
val cvtop : 'a ty -> 'a ty -> Signed.long -> 'a ty
val cvtop2 : 'a ty -> 'a ty -> 'a ty
val cx_approx_equal : 'a ty -> 'a ty -> int
val cx_approx0 : 'a ty -> 'a ty -> int
val gabs : 'a ty -> Signed.long -> 'a ty
val gaffect : 'a ty -> 'a ty -> unit
val gaffsg : Signed.long -> 'a ty -> unit
val gcmp : 'a ty -> 'a ty -> int
val gequal0 : 'a ty -> int
val gequal1 : 'a ty -> int
val gequalx : 'a ty -> int
val gequalm1 : 'a ty -> int
val gcmpsg : Signed.long -> 'a ty -> int
val gcvtop : 'a ty -> 'a ty -> Signed.long -> 'a ty
val gequal : 'a ty -> 'a ty -> bool
val gequalsg : Signed.long -> 'a ty -> int
val gexpo : 'a ty -> Signed.long
val gexpo_safe : 'a ty -> Signed.long
val gpexponent : 'a ty -> 'a ty
val gpvaluation : 'a ty -> 'a ty -> 'a ty
val gvaluation : 'a ty -> 'a ty -> Signed.long
val gidentical : 'a ty -> 'a ty -> int
val gmax : 'a ty -> 'a ty -> 'a ty
val gmaxgs : 'a ty -> Signed.long -> 'a ty
val gmin : 'a ty -> 'a ty -> 'a ty
val gmings : 'a ty -> Signed.long -> 'a ty
val gneg : 'a ty -> 'a ty
val gneg_i : 'a ty -> 'a ty
val gsigne : 'a ty -> int
val gtolist : 'a ty -> 'a ty
val gtolong : 'a ty -> Signed.long
val lexcmp : 'a ty -> 'a ty -> int
val listinsert : 'a ty -> 'a ty -> Signed.long -> 'a ty
val listpop : 'a ty -> Signed.long -> unit
val listpop0 : 'a ty -> Signed.long -> unit
val listput : 'a ty -> 'a ty -> Signed.long -> 'a ty
val listput0 : 'a ty -> 'a ty -> Signed.long -> unit
val listsort : 'a ty -> Signed.long -> unit
val matsize : 'a ty -> 'a ty
val mklist : unit -> 'a ty
val mklist_typ : Signed.long -> 'a ty
val mklistcopy : 'a ty -> 'a ty
val mkmap : unit -> 'a ty
val normalizeser : 'a ty -> 'a ty
val normalizepol : 'a ty -> 'a ty
val normalizepol_approx : 'a ty -> Signed.long -> 'a ty
val padic_to_fl : 'a ty -> pari_ulong -> pari_ulong
val padic_to_fp : 'a ty -> 'a ty -> 'a ty
val quadtofp : 'a ty -> Signed.long -> 'a ty
val sizedigit : 'a ty -> Signed.long
val u_lval : pari_ulong -> pari_ulong -> Signed.long

val u_lvalrem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val u_lvalrem_stop :
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  int Ctypes_static.ptr ->
  Signed.long

val u_pval : pari_ulong -> 'a ty -> Signed.long

val u_pvalrem :
  pari_ulong -> 'a ty -> pari_ulong Ctypes_static.ptr -> Signed.long

val vecindexmax : 'a ty -> Signed.long
val vecindexmin : 'a ty -> Signed.long
val vecmax0 : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val vecmax : 'a ty -> 'a ty
val vecmin0 : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val vecmin : 'a ty -> 'a ty
val z_lval : Signed.long -> pari_ulong -> Signed.long

val z_lvalrem :
  Signed.long -> pari_ulong -> Signed.long Ctypes_static.ptr -> Signed.long

val z_pval : Signed.long -> 'a ty -> Signed.long

val z_pvalrem :
  Signed.long -> 'a ty -> Signed.long Ctypes_static.ptr -> Signed.long

val zx_lval : 'a ty -> Signed.long -> Signed.long
val hgmcyclo : 'a ty -> 'a ty
val hgmalpha : 'a ty -> 'a ty
val hgmgamma : 'a ty -> 'a ty
val hgminit : 'a ty -> 'a ty -> 'a ty
val hgmparams : 'a ty -> 'a ty

val hgmeulerfactor :
  'a ty -> 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr -> 'a ty

val hgmcoef : 'a ty -> 'a ty -> 'a ty -> 'a ty
val hgmcoefs : 'a ty -> 'a ty -> Signed.long -> 'a ty
val hgmtwist : 'a ty -> 'a ty
val hgmissymmetrical : 'a ty -> Signed.long
val hgmbydegree : Signed.long -> 'a ty
val lfunhgm : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val qp_zeta : 'a ty -> 'a ty
val lerchphi : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val lerchzeta : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zetahurwitz : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val rgx_to_ser : 'a ty -> Signed.long -> 'a ty
val rgx_to_ser_inexact : 'a ty -> Signed.long -> 'a ty
val gtoser : 'a ty -> Signed.long -> Signed.long -> 'a ty
val gtoser_prec : 'a ty -> Signed.long -> Signed.long -> 'a ty
val rfrac_to_ser : 'a ty -> Signed.long -> 'a ty
val rfrac_to_ser_i : 'a ty -> Signed.long -> 'a ty
val rfracrecip_to_ser_absolute : 'a ty -> Signed.long -> 'a ty

val rfracrecip :
  'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> Signed.long

val scalarser : 'a ty -> Signed.long -> Signed.long -> 'a ty
val sertoser : 'a ty -> Signed.long -> 'a ty
val toser_i : 'a ty -> 'a ty
val rgv_to_ser : 'a ty -> Signed.long -> Signed.long -> 'a ty
val ser0 : 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty
val padic_to_q : 'a ty -> 'a ty
val padic_to_q_shallow : 'a ty -> 'a ty
val qpv_to_qv : 'a ty -> 'a ty
val rgc_rgv_mulrealsym : 'a ty -> 'a ty -> 'a ty
val rgm_mulreal : 'a ty -> 'a ty -> 'a ty
val rgx_cxeval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgx_deflate_max : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val rgx_deflate_order : 'a ty -> Signed.long
val rgx_degree : 'a ty -> Signed.long -> Signed.long
val rgx_integ : 'a ty -> 'a ty
val rgxy_cxevalx : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zx_deflate_order : 'a ty -> Signed.long
val zx_deflate_max : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val ceil_safe : 'a ty -> 'a ty
val ceilr : 'a ty -> 'a ty
val centerlift : 'a ty -> 'a ty
val centerlift0 : 'a ty -> Signed.long -> 'a ty
val compo : 'a ty -> Signed.long -> 'a ty
val deg1pol : 'a ty -> 'a ty -> Signed.long -> 'a ty
val deg1pol_shallow : 'a ty -> 'a ty -> Signed.long -> 'a ty
val deg2pol_shallow : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val denom : 'a ty -> 'a ty
val denom_i : 'a ty -> 'a ty
val denominator : 'a ty -> 'a ty -> 'a ty
val derivser : 'a ty -> 'a ty
val diffop : 'a ty -> 'a ty -> 'a ty -> 'a ty
val diffop0 : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val diviiround : 'a ty -> 'a ty -> 'a ty
val divrem : 'a ty -> 'a ty -> Signed.long -> 'a ty
val floor_safe : 'a ty -> 'a ty
val gceil : 'a ty -> 'a ty
val gcvtoi : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val gdeflate : 'a ty -> Signed.long -> Signed.long -> 'a ty
val gdivent : 'a ty -> 'a ty -> 'a ty
val gdiventgs : 'a ty -> Signed.long -> 'a ty
val gdiventsg : Signed.long -> 'a ty -> 'a ty
val gdiventres : 'a ty -> 'a ty -> 'a ty
val gdivmod : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val gdivround : 'a ty -> 'a ty -> 'a ty
val gdvd : 'a ty -> 'a ty -> int
val geq : 'a ty -> 'a ty -> 'a ty
val geval : 'a ty -> 'a ty
val gfloor : 'a ty -> 'a ty
val gtrunc2n : 'a ty -> Signed.long -> 'a ty
val gfrac : 'a ty -> 'a ty
val gge : 'a ty -> 'a ty -> 'a ty
val ggrando : 'a ty -> Signed.long -> 'a ty
val ggt : 'a ty -> 'a ty -> 'a ty
val gimag : 'a ty -> 'a ty
val gisexactzero : 'a ty -> 'a ty
val gle : 'a ty -> 'a ty -> 'a ty
val glt : 'a ty -> 'a ty -> 'a ty
val gmod : 'a ty -> 'a ty -> 'a ty
val gmodgs : 'a ty -> Signed.long -> 'a ty
val gmodsg : Signed.long -> 'a ty -> 'a ty
val gmodulo : 'a ty -> 'a ty -> 'a ty
val gmodulsg : Signed.long -> 'a ty -> 'a ty
val gmodulss : Signed.long -> Signed.long -> 'a ty
val gne : 'a ty -> 'a ty -> 'a ty
val gnot : 'a ty -> 'a ty
val gpolvar : 'a ty -> 'a ty
val gppadicprec : 'a ty -> 'a ty -> 'a ty
val gppoldegree : 'a ty -> Signed.long -> 'a ty
val gprecision : 'a ty -> Signed.long
val gpserprec : 'a ty -> Signed.long -> 'a ty
val greal : 'a ty -> 'a ty
val grndtoi : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val ground : 'a ty -> 'a ty
val gshift : 'a ty -> Signed.long -> 'a ty
val gsubst : 'a ty -> Signed.long -> 'a ty -> 'a ty
val gsubstpol : 'a ty -> 'a ty -> 'a ty -> 'a ty
val gsubstvec : 'a ty -> 'a ty -> 'a ty -> 'a ty
val gtocol : 'a ty -> 'a ty
val gtocol0 : 'a ty -> Signed.long -> 'a ty
val gtocolrev : 'a ty -> 'a ty
val gtocolrev0 : 'a ty -> Signed.long -> 'a ty
val gtopoly : 'a ty -> Signed.long -> 'a ty
val gtopolyrev : 'a ty -> Signed.long -> 'a ty
val gtovec : 'a ty -> 'a ty
val gtovec0 : 'a ty -> Signed.long -> 'a ty
val gtovecrev : 'a ty -> 'a ty
val gtovecrev0 : 'a ty -> Signed.long -> 'a ty
val gtovecsmall : 'a ty -> 'a ty
val gtovecsmall0 : 'a ty -> Signed.long -> 'a ty
val gtrunc : 'a ty -> 'a ty
val gvar : 'a ty -> Signed.long
val gvar2 : 'a ty -> Signed.long
val hqfeval : 'a ty -> 'a ty -> 'a ty
val imag_i : 'a ty -> 'a ty
val integ : 'a ty -> Signed.long -> 'a ty
val integser : 'a ty -> 'a ty
val ser_inv : 'a ty -> 'a ty
val iscomplex : 'a ty -> int
val isexactzero : 'a ty -> int
val isrationalzeroscalar : 'a ty -> int
val isinexact : 'a ty -> int
val isinexactreal : 'a ty -> int
val isint : 'a ty -> 'a ty Ctypes_static.ptr -> int
val isrationalzero : 'a ty -> int
val issmall : 'a ty -> Signed.long Ctypes_static.ptr -> int
val lift : 'a ty -> 'a ty
val lift_shallow : 'a ty -> 'a ty
val lift0 : 'a ty -> Signed.long -> 'a ty
val liftall : 'a ty -> 'a ty
val liftall_shallow : 'a ty -> 'a ty
val liftint : 'a ty -> 'a ty
val liftint_shallow : 'a ty -> 'a ty
val liftpol : 'a ty -> 'a ty
val liftpol_shallow : 'a ty -> 'a ty
val mkcoln : Signed.long -> 'a ty
val mkintn : Signed.long -> 'a ty
val mkpoln : Signed.long -> 'a ty
val mkvecn : Signed.long -> 'a ty
val mkvecsmalln : Signed.long -> 'a ty
val modrr_safe : 'a ty -> 'a ty -> 'a ty
val modrr_i : 'a ty -> 'a ty -> 'a ty -> 'a ty
val mulreal : 'a ty -> 'a ty -> 'a ty
val numer : 'a ty -> 'a ty
val numer_i : 'a ty -> 'a ty
val numerator : 'a ty -> 'a ty -> 'a ty
val padicprec : 'a ty -> 'a ty -> Signed.long
val padicprec_relative : 'a ty -> Signed.long
val precision : 'a ty -> Signed.long
val qf_apply_rgm : 'a ty -> 'a ty -> 'a ty
val qf_apply_zm : 'a ty -> 'a ty -> 'a ty
val qfb_apply_zm : 'a ty -> 'a ty -> 'a ty
val qfbil : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qfeval : 'a ty -> 'a ty -> 'a ty
val qfeval0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qfevalb : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qfnorm : 'a ty -> 'a ty -> 'a ty
val real_i : 'a ty -> 'a ty
val round0 : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val roundr : 'a ty -> 'a ty
val roundr_safe : 'a ty -> 'a ty
val scalarpol : 'a ty -> Signed.long -> 'a ty
val scalarpol_shallow : 'a ty -> Signed.long -> 'a ty
val ser_unscale : 'a ty -> 'a ty -> 'a ty
val serprec : 'a ty -> Signed.long -> Signed.long
val serreverse : 'a ty -> 'a ty
val simplify : 'a ty -> 'a ty
val simplify_shallow : 'a ty -> 'a ty
val tayl : 'a ty -> Signed.long -> Signed.long -> 'a ty
val trunc0 : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val uu32toi : pari_ulong -> pari_ulong -> 'a ty
val uu32toineg : pari_ulong -> pari_ulong -> 'a ty
val vars_sort_inplace : 'a ty -> 'a ty
val vars_to_rgxv : 'a ty -> 'a ty
val variables_vecsmall : 'a ty -> 'a ty
val variables_vec : 'a ty -> 'a ty
val genus2red : 'a ty -> 'a ty -> 'a ty
val genus2igusa : 'a ty -> Signed.long -> 'a ty
val gchar_conductor : 'a ty -> 'a ty -> 'a ty
val gchar_identify : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val gcharalgebraic : 'a ty -> 'a ty -> 'a ty
val gcharduallog : 'a ty -> 'a ty -> 'a ty
val gchareval : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val gchari_lfun : 'a ty -> 'a ty -> 'a ty -> 'a ty
val gcharinit : 'a ty -> 'a ty -> Signed.long -> 'a ty
val gcharisalgebraic : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> int

val gcharlocal :
  'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr -> 'a ty

val gcharlog : 'a ty -> 'a ty -> Signed.long -> 'a ty
val gcharnewprec : 'a ty -> Signed.long -> 'a ty
val is_gchar_group : 'a ty -> int
val lfungchar : 'a ty -> 'a ty -> 'a ty
val vecan_gchar : 'a ty -> Signed.long -> Signed.long -> 'a ty
val eulerf_gchar : 'a ty -> 'a ty -> Signed.long -> 'a ty
val group_ident : 'a ty -> 'a ty -> Signed.long
val group_ident_trans : 'a ty -> 'a ty -> Signed.long

val hash_create_ulong :
  pari_ulong -> Signed.long -> hashtable Ctypes.structure Ctypes_static.ptr

val hash_create_str :
  pari_ulong -> Signed.long -> hashtable Ctypes.structure Ctypes_static.ptr

val hash_create :
  pari_ulong ->
  (unit Ctypes_static.ptr -> pari_ulong) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> unit Ctypes_static.ptr -> int)
  Ctypes_static.static_funptr ->
  int ->
  hashtable Ctypes.structure Ctypes_static.ptr

val hash_dbg : hashtable Ctypes.structure Ctypes_static.ptr -> unit

val hash_haskey_gen :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  'a ty

val hash_haskey_long :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  int

val hash_init :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  pari_ulong ->
  (unit Ctypes_static.ptr -> pari_ulong) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> unit Ctypes_static.ptr -> int)
  Ctypes_static.static_funptr ->
  int ->
  unit

val hash_init_gen :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  pari_ulong ->
  ('a ty -> 'a ty -> int) Ctypes_static.static_funptr ->
  int ->
  unit

val hash_init_ulong :
  hashtable Ctypes.structure Ctypes_static.ptr -> pari_ulong -> int -> unit

val hash_insert :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  unit

val hash_insert_long :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  Signed.long ->
  unit

val hash_insert2 :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  pari_ulong ->
  unit

val hash_keys : hashtable Ctypes.structure Ctypes_static.ptr -> 'a ty
val hash_values : hashtable Ctypes.structure Ctypes_static.ptr -> 'a ty

val hash_search :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  hashentry Ctypes.structure Ctypes_static.ptr

val hash_search2 :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  pari_ulong ->
  hashentry Ctypes.structure Ctypes_static.ptr

val hash_select :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  hashentry Ctypes.structure Ctypes_static.ptr ->
  int)
  Ctypes_static.static_funptr ->
  hashentry Ctypes.structure Ctypes_static.ptr

val hash_remove :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  hashentry Ctypes.structure Ctypes_static.ptr

val hash_remove_select :
  hashtable Ctypes.structure Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  hashentry Ctypes.structure Ctypes_static.ptr ->
  int)
  Ctypes_static.static_funptr ->
  hashentry Ctypes.structure Ctypes_static.ptr

val hash_destroy : hashtable Ctypes.structure Ctypes_static.ptr -> unit
val hash_gen : 'a ty -> pari_ulong
val hash_zv : 'a ty -> pari_ulong
val zx_hyperellred : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val hyperellcharpoly : 'a ty -> 'a ty
val hyperellchangecurve : 'a ty -> 'a ty -> 'a ty
val hyperelldisc : 'a ty -> 'a ty
val hyperellisoncurve : 'a ty -> 'a ty -> int
val hyperellminimaldisc : 'a ty -> 'a ty -> 'a ty
val hyperellminimalmodel : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty -> 'a ty
val hyperellpadicfrobenius0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val hyperellpadicfrobenius : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val hyperellred : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val hypergeom : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val airy : 'a ty -> Signed.long -> 'a ty
val rgm_hnfall : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long -> 'a ty
val zm_hnf : 'a ty -> 'a ty
val zm_hnf_knapsack : 'a ty -> 'a ty
val zm_hnfall : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long -> 'a ty
val zm_hnfall_i : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long -> 'a ty
val zm_hnfcenter : 'a ty -> 'a ty
val zm_hnflll : 'a ty -> 'a ty Ctypes_static.ptr -> int -> 'a ty
val zv_extgcd : 'a ty -> 'a ty

val zv_snfall :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val zv_snf_group :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val zv_snf_rank_u : 'a ty -> pari_ulong -> Signed.long
val zv_snf_trunc : 'a ty -> unit
val zm_hnfmod : 'a ty -> 'a ty -> 'a ty
val zm_hnfmodall : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zm_hnfmodall_i : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zm_hnfmodid : 'a ty -> 'a ty -> 'a ty
val zm_hnfmodprime : 'a ty -> 'a ty -> 'a ty

val zm_hnfperm :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val zm_snfclean : 'a ty -> 'a ty -> 'a ty -> unit
val zm_snf : 'a ty -> 'a ty

val zm_snf_group :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val zm_snfall :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val zm_snfall_i :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long ->
  'a ty

val zv_snfclean : 'a ty -> 'a ty
val zpm_echelon : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val gsmith : 'a ty -> 'a ty
val gsmithall : 'a ty -> 'a ty
val hnf : 'a ty -> 'a ty
val hnf_divscale : 'a ty -> 'a ty -> 'a ty -> 'a ty
val hnf_invscale : 'a ty -> 'a ty -> 'a ty
val hnf_solve : 'a ty -> 'a ty -> 'a ty
val hnf_invimage : 'a ty -> 'a ty -> 'a ty
val hnfall : 'a ty -> 'a ty
val hnfdivide : 'a ty -> 'a ty -> int
val hnflll : 'a ty -> 'a ty
val hnfmerge_get_1 : 'a ty -> 'a ty -> 'a ty
val hnfmod : 'a ty -> 'a ty -> 'a ty
val hnfmodid : 'a ty -> 'a ty -> 'a ty
val hnfperm : 'a ty -> 'a ty
val matfrobenius : 'a ty -> Signed.long -> Signed.long -> 'a ty
val mathnf0 : 'a ty -> Signed.long -> 'a ty
val matsnf0 : 'a ty -> Signed.long -> 'a ty
val smith : 'a ty -> 'a ty
val smithall : 'a ty -> 'a ty
val smithclean : 'a ty -> 'a ty
val snfrank : 'a ty -> 'a ty -> Signed.long
val zlm_echelon : 'a ty -> Signed.long -> pari_ulong -> pari_ulong -> 'a ty
val zv_snf_rank : 'a ty -> pari_ulong -> Signed.long
val z_ecm : 'a ty -> Signed.long -> Signed.long -> pari_ulong -> 'a ty
val z_factor : 'a ty -> 'a ty
val z_factor_limit : 'a ty -> pari_ulong -> 'a ty
val z_factor_until : 'a ty -> 'a ty -> 'a ty
val z_issmooth : 'a ty -> pari_ulong -> Signed.long
val z_issmooth_fact : 'a ty -> pari_ulong -> 'a ty
val z_issquarefree : 'a ty -> Signed.long
val z_pollardbrent : 'a ty -> Signed.long -> Signed.long -> 'a ty
val absz_factor : 'a ty -> 'a ty
val absz_factor_limit : 'a ty -> pari_ulong -> 'a ty

val absz_factor_limit_strict :
  'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> 'a ty

val coreu : pari_ulong -> pari_ulong
val coreu_fact : 'a ty -> pari_ulong
val factorint : 'a ty -> Signed.long -> 'a ty
val factoru : pari_ulong -> 'a ty
val tridiv_boundu : pari_ulong -> pari_ulong
val ifac_isprime : 'a ty -> int

val ifac_next :
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  int

val ifac_read :
  'a ty -> 'a ty Ctypes_static.ptr -> Signed.long Ctypes_static.ptr -> int

val ifac_skip : 'a ty -> unit
val ifac_start : 'a ty -> int -> 'a ty

val is_357_power :
  'a ty -> 'a ty Ctypes_static.ptr -> pari_ulong Ctypes_static.ptr -> int

val is_pth_power :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  pari_ulong ->
  int

val ispowerful : 'a ty -> Signed.long
val maxomegau : pari_ulong -> Signed.long
val maxomegaoddu : pari_ulong -> Signed.long
val moebius : 'a ty -> Signed.long
val moebiusu : pari_ulong -> Signed.long
val moebiusu_fact : 'a ty -> Signed.long
val nextprime : 'a ty -> 'a ty
val precprime : 'a ty -> 'a ty
val radicalu : pari_ulong -> pari_ulong
val tridiv_bound : 'a ty -> pari_ulong

val uis_357_power :
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  int

val uis_357_powermod : pari_ulong -> pari_ulong Ctypes_static.ptr -> int
val unextprime : pari_ulong -> pari_ulong
val uprecprime : pari_ulong -> pari_ulong
val vecfactorsquarefreeu : pari_ulong -> pari_ulong -> 'a ty
val vecfactorsquarefreeu_coprime : pari_ulong -> pari_ulong -> 'a ty -> 'a ty
val vecfactoru_i : pari_ulong -> pari_ulong -> 'a ty
val vecfactoru : pari_ulong -> pari_ulong -> 'a ty
val vecfactoroddu_i : pari_ulong -> pari_ulong -> 'a ty
val vecfactoroddu : pari_ulong -> pari_ulong -> 'a ty
val vecsquarefreeu : pari_ulong -> pari_ulong -> 'a ty
val chk_gerepileupto : 'a ty -> int
val copy_bin : 'a ty -> genbin Ctypes.structure Ctypes_static.ptr
val copy_bin_canon : 'a ty -> genbin Ctypes.structure Ctypes_static.ptr
val dbg_gerepile : pari_ulong -> unit
val dbg_gerepileupto : 'a ty -> unit
val errname : 'a ty -> 'a ty
val gclone : 'a ty -> 'a ty
val gcloneref : 'a ty -> 'a ty
val gclone_refc : 'a ty -> unit
val gcopy : 'a ty -> 'a ty
val gcopy_avma : 'a ty -> pari_ulong Ctypes_static.ptr -> 'a ty
val gcopy_lg : 'a ty -> Signed.long -> 'a ty
val gerepile : pari_ulong -> pari_ulong -> 'a ty -> 'a ty
val gerepileallsp : pari_ulong -> pari_ulong -> int -> unit

val gerepilecoeffssp :
  pari_ulong -> pari_ulong -> Signed.long Ctypes_static.ptr -> int -> unit

val gerepilemanysp :
  pari_ulong ->
  pari_ulong ->
  'a ty Ctypes_static.ptr Ctypes_static.ptr ->
  int ->
  unit

val getheap : unit -> 'a ty
val gsizeword : 'a ty -> Signed.long
val gsizebyte : 'a ty -> Signed.long
val gunclone : 'a ty -> unit
val gunclone_deep : 'a ty -> unit
val listcopy : 'a ty -> 'a ty
val listinit : 'a ty -> 'a ty
val msgtimer : string -> unit
val name_numerr : string -> Signed.long
val new_chunk_resize : int -> unit
val newblock : int -> 'a ty
val numerr_name : Signed.long -> string
val obj_check : 'a ty -> Signed.long -> 'a ty

val obj_checkbuild :
  'a ty -> Signed.long -> ('a ty -> 'a ty) Ctypes_static.static_funptr -> 'a ty

val obj_checkbuild_padicprec :
  'a ty ->
  Signed.long ->
  ('a ty -> Signed.long -> 'a ty) Ctypes_static.static_funptr ->
  Signed.long ->
  'a ty

val obj_checkbuild_realprec :
  'a ty ->
  Signed.long ->
  ('a ty -> Signed.long -> 'a ty) Ctypes_static.static_funptr ->
  Signed.long ->
  'a ty

val obj_checkbuild_prec :
  'a ty ->
  Signed.long ->
  ('a ty -> Signed.long -> 'a ty) Ctypes_static.static_funptr ->
  ('a ty -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  'a ty

val obj_free : 'a ty -> unit
val obj_init : Signed.long -> Signed.long -> 'a ty
val obj_insert : 'a ty -> Signed.long -> 'a ty -> 'a ty
val obj_insert_shallow : 'a ty -> Signed.long -> 'a ty -> 'a ty
val obj_reinit : 'a ty -> 'a ty
val pari_add_function : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_defaults_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_close : unit -> unit
val pari_close_opts : pari_ulong -> unit
val pari_compile_str : string -> 'a ty
val pari_daemon : unit -> int
val pari_err : int -> unit
val pari_err_last : unit -> 'a ty
val pari_err2str : 'a ty -> string
val pari_init_opts : int -> pari_ulong -> pari_ulong -> unit
val pari_init : int -> pari_ulong -> unit
val pari_stackcheck_init : unit Ctypes_static.ptr -> unit
val pari_sighandler : int -> unit
val pari_sig_init : (int -> unit) Ctypes_static.static_funptr -> unit

val pari_thread_alloc :
  pari_thread Ctypes.structure Ctypes_static.ptr -> int -> 'a ty -> unit

val pari_thread_close : unit -> unit
val pari_thread_free : pari_thread Ctypes.structure Ctypes_static.ptr -> unit
val pari_thread_init : unit -> unit
val pari_thread_start : pari_thread Ctypes.structure Ctypes_static.ptr -> 'a ty

val pari_thread_valloc :
  pari_thread Ctypes.structure Ctypes_static.ptr -> int -> int -> 'a ty -> unit

val pari_version : unit -> 'a ty
val pari_warn : int -> unit
val paristack_newrsize : pari_ulong -> unit
val paristack_resize : pari_ulong -> unit
val paristack_setsize : int -> int -> unit
val parivstack_resize : pari_ulong -> unit
val parivstack_reset : unit -> unit
val setalldebug : Signed.long -> unit
val setdebug : string -> Signed.long -> 'a ty
val shiftaddress : 'a ty -> Signed.long -> unit
val shiftaddress_canon : 'a ty -> Signed.long -> unit
val timer : unit -> Signed.long
val timer_delay : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long
val timer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val timer_printf :
  pari_timer Ctypes.structure Ctypes_static.ptr -> string -> unit

val timer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit
val timer2 : unit -> Signed.long
val trap0 : string -> 'a ty -> 'a ty -> 'a ty

val traverseheap :
  ('a ty -> unit Ctypes_static.ptr -> unit) Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  unit

val walltimer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit

val walltimer_delay :
  pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val walltimer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long
val contfraceval : 'a ty -> 'a ty -> Signed.long -> 'a ty
val contfracinit : 'a ty -> Signed.long -> 'a ty

val intcirc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val intfuncinit :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  'a ty

val intnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val intnumgauss :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val intnumgaussinit : Signed.long -> Signed.long -> 'a ty
val intnuminit : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty

val intnumosc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty ->
  Signed.long ->
  'a ty

val intnumromb :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  'a ty

val intnumromb_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  'a ty

val prodeulerrat : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val prodnumrat : 'a ty -> Signed.long -> Signed.long -> 'a ty
val quodif : 'a ty -> Signed.long -> 'a ty
val sumeulerrat : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty

val sumnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val sumnumap :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val sumnumapinit : 'a ty -> Signed.long -> 'a ty
val sumnuminit : 'a ty -> Signed.long -> 'a ty
val sumnumlagrangeinit : 'a ty -> 'a ty -> Signed.long -> 'a ty

val sumnumlagrange :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val sumnummonien :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val sumnummonieninit : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val sumnumrat : 'a ty -> 'a ty -> Signed.long -> 'a ty

val sumnumsidi :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  float ->
  Signed.long ->
  'a ty

val z_isanypower : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val z_ispow2 : 'a ty -> Signed.long
val z_ispowerall : 'a ty -> pari_ulong -> 'a ty Ctypes_static.ptr -> Signed.long
val z_issquareall : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

val zn_ispower :
  'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long

val zn_issquare : 'a ty -> 'a ty -> Signed.long
val zp_issquare : 'a ty -> 'a ty -> Signed.long
val gisanypower : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val gissquare : 'a ty -> 'a ty
val gissquareall : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val ispolygonal : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val ispower : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val isprimepower : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val ispseudoprimepower : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val issquare : 'a ty -> Signed.long
val issquareall : 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long
val sqrtint : 'a ty -> 'a ty
val sqrtint0 : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val uisprimepower : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long
val uissquare : pari_ulong -> Signed.long
val uissquareall : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val ulogintall :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val padicfields0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val padicfields : 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val bnrclassfield : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val rnfkummer : 'a ty -> 'a ty -> Signed.long -> 'a ty
val is_linit : 'a ty -> Signed.long
val ldata_get_an : 'a ty -> 'a ty
val ldata_get_dual : 'a ty -> 'a ty
val ldata_get_gammavec : 'a ty -> 'a ty
val ldata_get_degree : 'a ty -> Signed.long
val ldata_get_k : 'a ty -> 'a ty
val ldata_get_k1 : 'a ty -> 'a ty
val ldata_get_conductor : 'a ty -> 'a ty
val ldata_get_rootno : 'a ty -> 'a ty
val ldata_get_residue : 'a ty -> 'a ty
val ldata_get_type : 'a ty -> Signed.long
val ldata_isreal : 'a ty -> Signed.long
val linit_get_type : 'a ty -> Signed.long
val linit_get_ldata : 'a ty -> 'a ty
val linit_get_tech : 'a ty -> 'a ty
val lfun_get_domain : 'a ty -> 'a ty
val lfun_get_dom : 'a ty -> 'a ty
val lfun_get_factgammavec : 'a ty -> 'a ty
val lfun_get_step : 'a ty -> 'a ty
val lfun_get_pol : 'a ty -> 'a ty
val lfun_get_residue : 'a ty -> 'a ty
val lfun_get_k2 : 'a ty -> 'a ty
val lfun_get_w2 : 'a ty -> 'a ty
val lfun_get_expot : 'a ty -> 'a ty
val lfun_get_bitprec : 'a ty -> Signed.long
val lfun : 'a ty -> 'a ty -> Signed.long -> 'a ty
val lfun0 : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfuncheckfeq : 'a ty -> 'a ty -> Signed.long -> Signed.long
val lfunconductor : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfuncost : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfuncost0 : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfuncreate : 'a ty -> 'a ty
val lfundual : 'a ty -> Signed.long -> 'a ty
val lfuneuler : 'a ty -> 'a ty -> Signed.long -> 'a ty
val lfunparams : 'a ty -> Signed.long -> 'a ty
val lfunan : 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfunhardy : 'a ty -> 'a ty -> Signed.long -> 'a ty
val lfuninit : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfuninit0 : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfuninit_make : Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val lfunlambda : 'a ty -> 'a ty -> Signed.long -> 'a ty
val lfunlambda0 : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfunmisc_to_ldata : 'a ty -> 'a ty
val lfunmisc_to_ldata_shallow : 'a ty -> 'a ty
val lfunmisc_to_ldata_shallow_i : 'a ty -> 'a ty
val lfunorderzero : 'a ty -> Signed.long -> Signed.long -> Signed.long
val lfunprod_get_fact : 'a ty -> 'a ty
val lfunrootno : 'a ty -> Signed.long -> 'a ty
val lfunrootres : 'a ty -> Signed.long -> 'a ty
val lfunrtopoles : 'a ty -> 'a ty
val lfunshift : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfuntwist : 'a ty -> 'a ty -> Signed.long -> 'a ty
val lfuntheta : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfunthetacost0 : 'a ty -> 'a ty -> Signed.long -> Signed.long -> Signed.long
val lfunthetacost : 'a ty -> 'a ty -> Signed.long -> Signed.long -> Signed.long
val lfunthetainit : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfunthetacheckinit : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfunzeros : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val sdomain_isincl : float -> 'a ty -> 'a ty -> int
val theta_get_an : 'a ty -> 'a ty
val theta_get_k : 'a ty -> 'a ty
val theta_get_r : 'a ty -> 'a ty
val theta_get_bitprec : 'a ty -> Signed.long
val theta_get_m : 'a ty -> Signed.long
val theta_get_tdom : 'a ty -> 'a ty
val theta_get_isqrtn : 'a ty -> 'a ty
val vgaeasytheta : 'a ty -> int
val znchargauss : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val dirzetak : 'a ty -> 'a ty -> 'a ty
val eta_zxn : Signed.long -> Signed.long -> 'a ty
val eta_product_zxn : 'a ty -> Signed.long -> 'a ty

val etaquotype :
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val galois_get_conj : 'a ty -> 'a ty
val ldata_vecan : 'a ty -> Signed.long -> Signed.long -> 'a ty
val ldata_newprec : 'a ty -> Signed.long -> 'a ty

val lfunabelianrelinit :
  'a ty -> 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty

val lfunartin : 'a ty -> 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val lfundiv : 'a ty -> 'a ty -> Signed.long -> 'a ty
val lfunellmfpeters : 'a ty -> Signed.long -> 'a ty
val lfunetaquo : 'a ty -> 'a ty
val lfungenus2 : 'a ty -> 'a ty
val lfunmfspec : 'a ty -> Signed.long -> 'a ty
val lfunmul : 'a ty -> 'a ty -> Signed.long -> 'a ty
val lfunqf : 'a ty -> Signed.long -> 'a ty
val lfunsympow : 'a ty -> pari_ulong -> 'a ty
val lfunzetakinit : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val qfiseven : 'a ty -> Signed.long
val lfunquadneg : Signed.long -> Signed.long -> 'a ty

val zm_lll_norms :
  'a ty -> float -> Signed.long -> 'a ty Ctypes_static.ptr -> 'a ty

val kerint : 'a ty -> 'a ty
val lllfp : 'a ty -> float -> Signed.long -> 'a ty
val lllgen : 'a ty -> 'a ty
val lllgram : 'a ty -> 'a ty
val lllgramgen : 'a ty -> 'a ty
val lllgramint : 'a ty -> 'a ty
val lllgramkerim : 'a ty -> 'a ty
val lllgramkerimgen : 'a ty -> 'a ty
val lllint : 'a ty -> 'a ty
val lllintpartial : 'a ty -> 'a ty
val lllintpartial_inplace : 'a ty -> 'a ty
val lllkerim : 'a ty -> 'a ty
val lllkerimgen : 'a ty -> 'a ty
val matkerint0 : 'a ty -> Signed.long -> 'a ty
val qflll0 : 'a ty -> Signed.long -> 'a ty
val qflllgram0 : 'a ty -> Signed.long -> 'a ty
val gtomap : 'a ty -> 'a ty
val mapdelete : 'a ty -> 'a ty -> unit
val mapdomain : 'a ty -> 'a ty
val mapdomain_shallow : 'a ty -> 'a ty
val mapget : 'a ty -> 'a ty -> 'a ty
val mapisdefined : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> int
val mapput : 'a ty -> 'a ty -> 'a ty -> unit
val maptomat : 'a ty -> 'a ty
val maptomat_shallow : 'a ty -> 'a ty
val matpermanent : 'a ty -> 'a ty
val zm_permanent : 'a ty -> 'a ty
val dbllemma526 : float -> float -> float -> float -> float
val dblcoro526 : float -> float -> float -> float
val gammamellininv : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val gammamellininvasymp : 'a ty -> Signed.long -> Signed.long -> 'a ty
val gammamellininvinit : 'a ty -> Signed.long -> Signed.long -> 'a ty
val gammamellininvrt : 'a ty -> 'a ty -> Signed.long -> 'a ty
val member_a1 : 'a ty -> 'a ty
val member_a2 : 'a ty -> 'a ty
val member_a3 : 'a ty -> 'a ty
val member_a4 : 'a ty -> 'a ty
val member_a6 : 'a ty -> 'a ty
val member_area : 'a ty -> 'a ty
val member_b2 : 'a ty -> 'a ty
val member_b4 : 'a ty -> 'a ty
val member_b6 : 'a ty -> 'a ty
val member_b8 : 'a ty -> 'a ty
val member_bid : 'a ty -> 'a ty
val member_bnf : 'a ty -> 'a ty
val member_c4 : 'a ty -> 'a ty
val member_c6 : 'a ty -> 'a ty
val member_clgp : 'a ty -> 'a ty
val member_codiff : 'a ty -> 'a ty
val member_cyc : 'a ty -> 'a ty
val member_diff : 'a ty -> 'a ty
val member_disc : 'a ty -> 'a ty
val member_e : 'a ty -> 'a ty
val member_eta : 'a ty -> 'a ty
val member_f : 'a ty -> 'a ty
val member_fu : 'a ty -> 'a ty
val member_gen : 'a ty -> 'a ty
val member_group : 'a ty -> 'a ty
val member_index : 'a ty -> 'a ty
val member_j : 'a ty -> 'a ty
val member_mod : 'a ty -> 'a ty
val member_nf : 'a ty -> 'a ty
val member_no : 'a ty -> 'a ty
val member_omega : 'a ty -> 'a ty
val member_orders : 'a ty -> 'a ty
val member_p : 'a ty -> 'a ty
val member_pol : 'a ty -> 'a ty
val member_polabs : 'a ty -> 'a ty
val member_reg : 'a ty -> 'a ty
val member_r1 : 'a ty -> 'a ty
val member_r2 : 'a ty -> 'a ty
val member_roots : 'a ty -> 'a ty
val member_sign : 'a ty -> 'a ty
val member_t2 : 'a ty -> 'a ty
val member_tate : 'a ty -> 'a ty
val member_tu : 'a ty -> 'a ty
val member_zk : 'a ty -> 'a ty
val member_zkst : 'a ty -> 'a ty
val qab_tracerel : 'a ty -> Signed.long -> 'a ty -> 'a ty
val qabm_tracerel : 'a ty -> Signed.long -> 'a ty -> 'a ty
val qabv_tracerel : 'a ty -> Signed.long -> 'a ty -> 'a ty
val qab_trace_init : Signed.long -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val checkmf : 'a ty -> 'a ty
val checkmf_i : 'a ty -> int
val getcache : unit -> 'a ty
val hclassno6u : pari_ulong -> pari_ulong
val hclassno6u_no_cache : pari_ulong -> pari_ulong
val lfunmf : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mfdelta : unit -> 'a ty
val mfeh : 'a ty -> 'a ty
val mfek : Signed.long -> 'a ty
val mftheta : 'a ty -> 'a ty
val mf_get_chi : 'a ty -> 'a ty
val mf_get_n : 'a ty -> Signed.long
val mf_get_nk : 'a ty -> 'a ty
val mf_get_field : 'a ty -> 'a ty
val mf_get_gn : 'a ty -> 'a ty
val mf_get_gk : 'a ty -> 'a ty
val mf_get_k : 'a ty -> Signed.long
val mf_get_r : 'a ty -> Signed.long
val mf_get_type : 'a ty -> Signed.long
val mfatkin : 'a ty -> 'a ty -> 'a ty
val mfatkineigenvalues : 'a ty -> Signed.long -> Signed.long -> 'a ty
val mfatkininit : 'a ty -> Signed.long -> Signed.long -> 'a ty
val mfbasis : 'a ty -> Signed.long -> 'a ty
val mfbd : 'a ty -> Signed.long -> 'a ty
val mfbracket : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mfcharorder : 'a ty -> Signed.long
val mfcharmodulus : 'a ty -> Signed.long
val mfcharpol : 'a ty -> 'a ty
val mfcoef : 'a ty -> Signed.long -> 'a ty
val mfcoefs : 'a ty -> Signed.long -> Signed.long -> 'a ty
val mfconductor : 'a ty -> 'a ty -> Signed.long
val mfcosets : 'a ty -> 'a ty
val mfcuspdim : Signed.long -> Signed.long -> 'a ty -> Signed.long
val mfcuspisregular : 'a ty -> 'a ty -> Signed.long
val mfcusps : 'a ty -> 'a ty
val mfcuspval : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val mfcuspwidth : 'a ty -> 'a ty -> Signed.long
val mfderiv : 'a ty -> Signed.long -> 'a ty
val mfderive2 : 'a ty -> Signed.long -> 'a ty
val mfdescribe : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val mfdim : 'a ty -> Signed.long -> 'a ty
val mfdiv : 'a ty -> 'a ty -> 'a ty
val mfdiv_val : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mfeigenbasis : 'a ty -> 'a ty
val mfeigensearch : 'a ty -> 'a ty -> 'a ty
val mfeisenstein : Signed.long -> 'a ty -> 'a ty -> 'a ty
val mfeisensteindim : Signed.long -> Signed.long -> 'a ty -> Signed.long
val mfembed : 'a ty -> 'a ty -> 'a ty
val mfembed0 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mfeval : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val mffields : 'a ty -> 'a ty
val mffromell : 'a ty -> 'a ty
val mffrometaquo : 'a ty -> Signed.long -> 'a ty
val mffromlfun : 'a ty -> Signed.long -> 'a ty
val mffromqf : 'a ty -> 'a ty -> 'a ty
val mffulldim : Signed.long -> Signed.long -> 'a ty -> Signed.long
val mfgaloisprojrep : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mfgaloistype : 'a ty -> 'a ty -> 'a ty
val mfhecke : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mfheckemat : 'a ty -> 'a ty -> 'a ty
val mfinit : 'a ty -> Signed.long -> 'a ty
val mfiscm : 'a ty -> 'a ty
val mfiscuspidal : 'a ty -> 'a ty -> Signed.long
val mfisequal : 'a ty -> 'a ty -> Signed.long -> Signed.long
val mfisetaquo : 'a ty -> Signed.long -> 'a ty
val mfkohnenbasis : 'a ty -> 'a ty
val mfkohnenbijection : 'a ty -> 'a ty
val mfkohneneigenbasis : 'a ty -> 'a ty -> 'a ty
val mflinear : 'a ty -> 'a ty -> 'a ty
val mfmanin : 'a ty -> Signed.long -> 'a ty
val mfmatembed : 'a ty -> 'a ty -> 'a ty
val mfmul : 'a ty -> 'a ty -> 'a ty
val mfnewdim : Signed.long -> Signed.long -> 'a ty -> Signed.long
val mfolddim : Signed.long -> Signed.long -> 'a ty -> Signed.long
val mfparams : 'a ty -> 'a ty
val mfperiodpol : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val mfperiodpolbasis : Signed.long -> Signed.long -> 'a ty
val mfpetersson : 'a ty -> 'a ty -> 'a ty
val mfpow : 'a ty -> Signed.long -> 'a ty
val mfsearch : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mfshift : 'a ty -> Signed.long -> 'a ty
val mfshimura : 'a ty -> 'a ty -> Signed.long -> 'a ty

val mfslashexpansion :
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  'a ty Ctypes_static.ptr ->
  Signed.long ->
  'a ty

val mfspace : 'a ty -> 'a ty -> Signed.long
val mfsplit : 'a ty -> Signed.long -> Signed.long -> 'a ty
val mfsturm : 'a ty -> Signed.long
val mfsturmngk : Signed.long -> 'a ty -> Signed.long
val mfsturmnk : Signed.long -> Signed.long -> Signed.long
val mfsturm_mf : 'a ty -> Signed.long
val mfsymboleval : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val mfsymbol : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mftaylor : 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val mftobasis : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mftobasises : 'a ty -> 'a ty -> 'a ty
val mftocol : 'a ty -> Signed.long -> Signed.long -> 'a ty
val mftocoset : pari_ulong -> 'a ty -> 'a ty -> 'a ty
val mftonew : 'a ty -> 'a ty -> 'a ty
val mftraceform : 'a ty -> Signed.long -> 'a ty
val mftwist : 'a ty -> 'a ty -> 'a ty
val mfvecembed : 'a ty -> 'a ty -> 'a ty
val mfvectomat : 'a ty -> Signed.long -> Signed.long -> 'a ty
val fl_inv : pari_ulong -> pari_ulong -> pari_ulong
val fl_invsafe : pari_ulong -> pari_ulong -> pari_ulong

val fp_ratlift :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  int

val zm2_mul : 'a ty -> 'a ty -> 'a ty
val abscmpii : 'a ty -> 'a ty -> int
val abscmprr : 'a ty -> 'a ty -> int
val absequalii : 'a ty -> 'a ty -> int
val addii_sign : 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty
val addir_sign : 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty
val addmulii : 'a ty -> 'a ty -> 'a ty -> 'a ty
val addmulii_inplace : 'a ty -> 'a ty -> 'a ty -> 'a ty
val addrr_sign : 'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty
val addsi_sign : Signed.long -> 'a ty -> Signed.long -> 'a ty
val addsr : Signed.long -> 'a ty -> 'a ty
val addui_sign : pari_ulong -> 'a ty -> Signed.long -> 'a ty
val addumului : pari_ulong -> pari_ulong -> 'a ty -> 'a ty
val affir : 'a ty -> 'a ty -> unit
val affrr : 'a ty -> 'a ty -> unit

val cbezout :
  Signed.long ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val cgcd : Signed.long -> Signed.long -> Signed.long
val clcm : Signed.long -> Signed.long -> Signed.long
val cmpii : 'a ty -> 'a ty -> int
val cmprr : 'a ty -> 'a ty -> int
val dblexpo : float -> Signed.long
val dblmantissa : float -> pari_ulong
val dbltor : float -> 'a ty
val divir : 'a ty -> 'a ty -> 'a ty
val divis : 'a ty -> Signed.long -> 'a ty
val divis_rem : 'a ty -> Signed.long -> Signed.long Ctypes_static.ptr -> 'a ty
val absdiviu_rem : 'a ty -> pari_ulong -> pari_ulong Ctypes_static.ptr -> 'a ty
val diviuuexact : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val diviuexact : 'a ty -> pari_ulong -> 'a ty
val divri : 'a ty -> 'a ty -> 'a ty
val divrr : 'a ty -> 'a ty -> 'a ty
val divrs : 'a ty -> Signed.long -> 'a ty
val divru : 'a ty -> pari_ulong -> 'a ty
val divsi : Signed.long -> 'a ty -> 'a ty
val divsr : Signed.long -> 'a ty -> 'a ty
val divur : pari_ulong -> 'a ty -> 'a ty
val dvmdii : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val equalrr : 'a ty -> 'a ty -> int
val floorr : 'a ty -> 'a ty
val halfgcdii : 'a ty -> 'a ty -> 'a ty
val int2n : Signed.long -> 'a ty
val int2u : pari_ulong -> 'a ty
val int2um1 : pari_ulong -> 'a ty
val int_normalize : 'a ty -> Signed.long -> 'a ty
val invmod : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> int
val invmod2bil : pari_ulong -> pari_ulong
val invr : 'a ty -> 'a ty
val mantissa_real : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val modiiz : 'a ty -> 'a ty -> 'a ty -> unit
val mulir : 'a ty -> 'a ty -> 'a ty
val mulrr : 'a ty -> 'a ty -> 'a ty
val mulsi : Signed.long -> 'a ty -> 'a ty
val mulsr : Signed.long -> 'a ty -> 'a ty
val mulss : Signed.long -> Signed.long -> 'a ty
val mului : pari_ulong -> 'a ty -> 'a ty
val mulur : pari_ulong -> 'a ty -> 'a ty
val muluu : pari_ulong -> pari_ulong -> 'a ty
val muluui : pari_ulong -> pari_ulong -> 'a ty -> 'a ty
val pari_kernel_close : unit -> unit
val pari_kernel_init : unit -> unit
val pari_kernel_version : unit -> string
val remi2n : 'a ty -> Signed.long -> 'a ty
val rtodbl : 'a ty -> float
val shifti : 'a ty -> Signed.long -> 'a ty
val sqri : 'a ty -> 'a ty
val sqrr : 'a ty -> 'a ty
val sqrs : Signed.long -> 'a ty
val sqrtr_abs : 'a ty -> 'a ty
val sqrtremi : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val sqru : pari_ulong -> 'a ty
val subsr : Signed.long -> 'a ty -> 'a ty
val truedvmdii : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val truedvmdis : 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr -> 'a ty
val truedvmdsi : Signed.long -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val trunc2nr : 'a ty -> Signed.long -> 'a ty
val mantissa2nr : 'a ty -> Signed.long -> 'a ty
val truncr : 'a ty -> 'a ty
val ugcd : pari_ulong -> pari_ulong -> pari_ulong
val ulcm : pari_ulong -> pari_ulong -> pari_ulong
val umodiu : 'a ty -> pari_ulong -> pari_ulong
val vals : pari_ulong -> Signed.long
val fpc_ratlift : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpm_ratlift : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_ratlift : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxqx_gcd : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zxqx_gcd : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rnfabelianconjgen : 'a ty -> 'a ty -> 'a ty
val rnfisabelian : 'a ty -> 'a ty -> Signed.long

val forpart :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  'a ty ->
  'a ty ->
  unit

val forpart_init :
  forpart_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  'a ty ->
  'a ty ->
  unit

val forpart_next : forpart_t Ctypes.structure Ctypes_static.ptr -> 'a ty
val forpart_prev : forpart_t Ctypes.structure Ctypes_static.ptr -> 'a ty
val numbpart : 'a ty -> 'a ty
val partitions : Signed.long -> 'a ty -> 'a ty -> 'a ty

val forperm :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  'a ty ->
  unit

val forperm_init : forperm_t Ctypes.structure Ctypes_static.ptr -> 'a ty -> unit
val forperm_next : forperm_t Ctypes.structure Ctypes_static.ptr -> 'a ty

val forallsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> Signed.long -> unit

val forksubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  unit

val forsubset_next : forsubset_t Ctypes.structure Ctypes_static.ptr -> 'a ty

val forsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> 'a ty -> unit

val glambertw : 'a ty -> Signed.long -> Signed.long -> 'a ty
val mplambertw : 'a ty -> Signed.long -> 'a ty
val mplambertx : 'a ty -> Signed.long -> 'a ty
val mplambertx_logx : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mplambertxlogx_x : 'a ty -> 'a ty -> Signed.long -> 'a ty
val z_to_perm : Signed.long -> 'a ty -> 'a ty
val abelian_group : 'a ty -> 'a ty
val conjclasses_repr : 'a ty -> Signed.long -> 'a ty
val cyc_pow : 'a ty -> Signed.long -> 'a ty
val cyc_pow_perm : 'a ty -> Signed.long -> 'a ty
val cyclicgroup : 'a ty -> Signed.long -> 'a ty
val dicyclicgroup : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val group_abelianhnf : 'a ty -> 'a ty -> 'a ty
val group_abeliansnf : 'a ty -> 'a ty -> 'a ty
val group_domain : 'a ty -> Signed.long
val group_elts : 'a ty -> Signed.long -> 'a ty
val group_export : 'a ty -> Signed.long -> 'a ty
val group_export_gap : 'a ty -> 'a ty
val group_export_magma : 'a ty -> 'a ty
val group_isa4s4 : 'a ty -> Signed.long
val group_isabelian : 'a ty -> Signed.long
val group_leftcoset : 'a ty -> 'a ty -> 'a ty
val group_order : 'a ty -> Signed.long
val group_perm_normalize : 'a ty -> 'a ty -> Signed.long
val group_quotient : 'a ty -> 'a ty -> 'a ty
val group_rightcoset : 'a ty -> 'a ty -> 'a ty
val group_set : 'a ty -> Signed.long -> 'a ty
val group_subgroup_is_faithful : 'a ty -> 'a ty -> int
val group_subgroup_isnormal : 'a ty -> 'a ty -> Signed.long
val group_subgroups : 'a ty -> 'a ty
val groupelts_solvablesubgroups : 'a ty -> 'a ty
val group_to_cc : 'a ty -> 'a ty
val groupelts_abelian_group : 'a ty -> 'a ty
val groupelts_center : 'a ty -> 'a ty
val groupelts_conj_set : 'a ty -> 'a ty -> 'a ty
val groupelts_conjclasses : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val groupelts_exponent : 'a ty -> Signed.long
val groupelts_quotient : 'a ty -> 'a ty -> 'a ty
val groupelts_set : 'a ty -> Signed.long -> 'a ty
val groupelts_to_group : 'a ty -> 'a ty
val numtoperm : Signed.long -> 'a ty -> 'a ty
val perm_commute : 'a ty -> 'a ty -> int
val perm_cycles : 'a ty -> 'a ty
val perm_order : 'a ty -> 'a ty
val perm_orderu : 'a ty -> pari_ulong
val perm_pow : 'a ty -> 'a ty -> 'a ty
val perm_powu : 'a ty -> pari_ulong -> 'a ty
val perm_sign : 'a ty -> Signed.long
val perm_to_gap : 'a ty -> 'a ty
val perm_to_z : 'a ty -> 'a ty
val permcycles : 'a ty -> 'a ty
val permorder : 'a ty -> 'a ty
val permsign : 'a ty -> Signed.long
val permtonum : 'a ty -> 'a ty
val quotient_group : 'a ty -> 'a ty -> 'a ty
val quotient_groupelts : 'a ty -> 'a ty
val quotient_perm : 'a ty -> 'a ty -> 'a ty
val quotient_subgroup_lift : 'a ty -> 'a ty -> 'a ty -> 'a ty
val subgroups_tableset : 'a ty -> Signed.long -> 'a ty
val tableset_find_index : 'a ty -> 'a ty -> Signed.long
val trivialgroup : unit -> 'a ty
val vec_insert : 'a ty -> Signed.long -> 'a ty -> 'a ty
val vec_is1to1 : 'a ty -> int
val vec_isconst : 'a ty -> int
val vecperm_orbits : 'a ty -> Signed.long -> 'a ty
val vecsmall_duplicate : 'a ty -> Signed.long
val vecsmall_duplicate_sorted : 'a ty -> Signed.long
val vecsmall_indexsort : 'a ty -> 'a ty
val vecsmall_is1to1 : 'a ty -> int
val vecsmall_isconst : 'a ty -> int
val vecsmall_sort : 'a ty -> unit
val vecsmall_uniq : 'a ty -> 'a ty
val vecsmall_uniq_sorted : 'a ty -> 'a ty
val vecsmall_counting_indexsort : 'a ty -> Signed.long -> 'a ty
val vecsmall_counting_sort : 'a ty -> Signed.long -> unit
val vecsmall_counting_uniq : 'a ty -> Signed.long -> 'a ty
val vecvecsmall_indexsort : 'a ty -> 'a ty
val vecvecsmall_max : 'a ty -> Signed.long
val vecvecsmall_search : 'a ty -> 'a ty -> Signed.long
val vecvecsmall_sort : 'a ty -> 'a ty
val vecvecsmall_sort_inplace : 'a ty -> 'a ty Ctypes_static.ptr -> unit
val vecvecsmall_sort_shallow : 'a ty -> 'a ty
val vecvecsmall_sort_uniq : 'a ty -> 'a ty
val mt_broadcast : 'a ty -> unit
val mt_nbthreads : unit -> Signed.long
val mt_queue_end : pari_mt Ctypes.structure Ctypes_static.ptr -> unit

val mt_queue_get :
  pari_mt Ctypes.structure Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  'a ty

val mt_queue_start : pari_mt Ctypes.structure Ctypes_static.ptr -> 'a ty -> unit

val mt_queue_start_lim :
  pari_mt Ctypes.structure Ctypes_static.ptr -> 'a ty -> Signed.long -> unit

val mt_queue_submit :
  pari_mt Ctypes.structure Ctypes_static.ptr -> Signed.long -> 'a ty -> unit

val mt_sigint_block : unit -> unit
val mt_sigint_unblock : unit -> unit
val pari_mt_init : unit -> unit
val pari_mt_close : unit -> unit
val subcyclopclgp : 'a ty -> 'a ty -> Signed.long -> 'a ty
val subcycloiwasawa : 'a ty -> 'a ty -> Signed.long -> 'a ty
val subcyclohminus : 'a ty -> 'a ty -> 'a ty

val color_to_rgb :
  'a ty ->
  int Ctypes_static.ptr ->
  int Ctypes_static.ptr ->
  int Ctypes_static.ptr ->
  unit

val colorname_to_rgb :
  string ->
  int Ctypes_static.ptr ->
  int Ctypes_static.ptr ->
  int Ctypes_static.ptr ->
  unit

val long_to_rgb :
  Signed.long ->
  int Ctypes_static.ptr ->
  int Ctypes_static.ptr ->
  int Ctypes_static.ptr ->
  unit

val pari_plot_by_file : string -> string -> string -> unit

val pari_set_plot_engine :
  (pari_plot Ctypes.structure Ctypes_static.ptr -> unit)
  Ctypes_static.static_funptr ->
  unit

val pari_kill_plot_engine : unit -> unit

val parploth :
  'a ty -> 'a ty -> 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty

val parplothexport :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  'a ty

val plotbox : Signed.long -> 'a ty -> 'a ty -> Signed.long -> unit
val plotclip : Signed.long -> unit
val plotcolor : Signed.long -> 'a ty -> 'a ty

val plotcopy :
  Signed.long -> Signed.long -> 'a ty -> 'a ty -> Signed.long -> unit

val plotcursor : Signed.long -> 'a ty
val plotdraw : 'a ty -> Signed.long -> unit
val plotexport : 'a ty -> 'a ty -> Signed.long -> 'a ty

val ploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  'a ty

val plothexport :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  'a ty

val plothraw : 'a ty -> 'a ty -> Signed.long -> 'a ty
val plothrawexport : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val plothsizes : Signed.long -> 'a ty
val plotinit : Signed.long -> 'a ty -> 'a ty -> Signed.long -> unit
val plotkill : Signed.long -> unit
val plotline : Signed.long -> 'a ty -> 'a ty -> unit
val plotlines : Signed.long -> 'a ty -> 'a ty -> Signed.long -> unit
val plotlinetype : Signed.long -> Signed.long -> unit
val plotmove : Signed.long -> 'a ty -> 'a ty -> unit
val plotpoints : Signed.long -> 'a ty -> 'a ty -> unit
val plotpointsize : Signed.long -> 'a ty -> unit
val plotpointtype : Signed.long -> Signed.long -> unit
val plotrbox : Signed.long -> 'a ty -> 'a ty -> Signed.long -> unit

val plotrecth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  Signed.long ->
  'a ty ->
  'a ty ->
  pari_ulong ->
  Signed.long ->
  Signed.long ->
  'a ty

val plotrecthraw : Signed.long -> 'a ty -> Signed.long -> 'a ty
val plotrline : Signed.long -> 'a ty -> 'a ty -> unit
val plotrmove : Signed.long -> 'a ty -> 'a ty -> unit
val plotrpoint : Signed.long -> 'a ty -> 'a ty -> unit
val plotscale : Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> unit
val plotstring : Signed.long -> string -> Signed.long -> unit
val psdraw : 'a ty -> Signed.long -> unit

val psploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  'a ty

val psplothraw : 'a ty -> 'a ty -> Signed.long -> 'a ty

val rect2ps :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  string

val rect2ps_i :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  int ->
  string

val rect2svg :
  'a ty ->
  'a ty ->
  'a ty ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  string

val pariplot :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  unit

val zx_zp_root : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val zp_appr : 'a ty -> 'a ty -> 'a ty
val cmp_padic : 'a ty -> 'a ty -> int
val factorpadic : 'a ty -> 'a ty -> Signed.long -> 'a ty
val gdeuc : 'a ty -> 'a ty -> 'a ty
val grem : 'a ty -> 'a ty -> 'a ty
val padicappr : 'a ty -> 'a ty -> 'a ty
val flv_factorback : 'a ty -> 'a ty -> pari_ulong -> pari_ulong
val flxqv_factorback : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val fpv_factorback : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqv_factorback : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val q_content : 'a ty -> 'a ty
val q_content_safe : 'a ty -> 'a ty
val q_denom : 'a ty -> 'a ty
val q_denom_safe : 'a ty -> 'a ty
val q_div_to_int : 'a ty -> 'a ty -> 'a ty
val q_gcd : 'a ty -> 'a ty -> 'a ty
val q_mul_to_int : 'a ty -> 'a ty -> 'a ty
val q_muli_to_int : 'a ty -> 'a ty -> 'a ty
val q_primitive_part : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val q_primpart : 'a ty -> 'a ty
val q_remove_denom : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val q_factor : 'a ty -> 'a ty
val q_factor_limit : 'a ty -> pari_ulong -> 'a ty

val rg_type :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rgc_type :
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rescale_to_int : 'a ty -> 'a ty

val rgm_type :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_type2 :
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type2 :
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_rg_type :
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_chinese_coprime : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgx_disc : 'a ty -> 'a ty

val rgx_extgcd :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val rgx_extgcd_simple :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val rgx_gcd : 'a ty -> 'a ty -> 'a ty
val rgx_gcd_simple : 'a ty -> 'a ty -> 'a ty
val rgx_halfgcd : 'a ty -> 'a ty -> 'a ty

val rgx_halfgcd_all :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val rgx_rescale_to_int : 'a ty -> 'a ty
val rgx_resultant_all : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val rgx_sturmpart : 'a ty -> 'a ty -> Signed.long
val rgx_sylvestermatrix : 'a ty -> 'a ty -> 'a ty

val rgx_type :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type2 :
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type3 :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type_decode :
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val rgx_type_is_composite : Signed.long -> int
val rgxq_charpoly : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgxq_inv : 'a ty -> 'a ty -> 'a ty
val rgxq_minpoly : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rgxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty

val rgxq_ratlift :
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  int

val rgxq_sqr : 'a ty -> 'a ty -> 'a ty
val z_content : 'a ty -> 'a ty
val zx_content : 'a ty -> 'a ty
val centermod : 'a ty -> 'a ty -> 'a ty
val centermod_i : 'a ty -> 'a ty -> 'a ty -> 'a ty
val centermodii : 'a ty -> 'a ty -> 'a ty -> 'a ty
val content : 'a ty -> 'a ty
val content0 : 'a ty -> 'a ty -> 'a ty
val deg1_from_roots : 'a ty -> Signed.long -> 'a ty
val factor0 : 'a ty -> 'a ty -> 'a ty
val factorback : 'a ty -> 'a ty
val factorback2 : 'a ty -> 'a ty -> 'a ty

val gbezout :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val gdivexact : 'a ty -> 'a ty -> 'a ty

val gen_factorback :
  'a ty ->
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty -> 'a ty)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> 'a ty) Ctypes_static.static_funptr ->
  'a ty

val ggcd : 'a ty -> 'a ty -> 'a ty
val ggcd0 : 'a ty -> 'a ty -> 'a ty
val ghalfgcd : 'a ty -> 'a ty -> 'a ty
val ginvmod : 'a ty -> 'a ty -> 'a ty
val glcm : 'a ty -> 'a ty -> 'a ty
val glcm0 : 'a ty -> 'a ty -> 'a ty
val newtonpoly : 'a ty -> 'a ty -> 'a ty
val primitive_part : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val primpart : 'a ty -> 'a ty
val reduceddiscsmith : 'a ty -> 'a ty
val resultant2 : 'a ty -> 'a ty -> 'a ty
val resultant : 'a ty -> 'a ty -> 'a ty
val rnfcharpoly : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val roots_from_deg1 : 'a ty -> 'a ty
val roots_to_pol : 'a ty -> Signed.long -> 'a ty
val roots_to_pol_r1 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val sturmpart : 'a ty -> 'a ty -> 'a ty -> Signed.long

val subresext :
  'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> 'a ty

val sylvestermatrix : 'a ty -> 'a ty -> 'a ty
val trivial_fact : unit -> 'a ty
val gcdext0 : 'a ty -> 'a ty -> 'a ty
val prime_fact : 'a ty -> 'a ty
val row_q_primpart : 'a ty -> 'a ty
val vec_q_primpart : 'a ty -> 'a ty
val vecprod : 'a ty -> 'a ty
val zv_lcm : 'a ty -> 'a ty
val flx_flxy_resultant : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxx_resultant : 'a ty -> 'a ty -> pari_ulong -> Signed.long -> 'a ty
val fpx_fpxy_resultant : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_translate : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxqx_normalize : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxv_fpc_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxy_fpxq_evaly : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val fpxc_center : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxm_center : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_fp_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_div : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_halve : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_invsafe : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_mulu : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fq_neg : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_neg_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fq_sqr : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_sqrt : 'a ty -> 'a ty -> 'a ty -> 'a ty

val fq_sqrtn :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val fq_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqc_fq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqc_fqv_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqc_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqc_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqv_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqv_roots_to_pol : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val fqx_fq_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_fq_mul_to_monic : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_fq_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_translate : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fqxq_matrix_pow :
  'a ty -> Signed.long -> Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fqxq_powers : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxy_eval : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxy_evalx : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val qx_disc : 'a ty -> 'a ty
val qx_gcd : 'a ty -> 'a ty -> 'a ty
val qx_resultant : 'a ty -> 'a ty -> 'a ty
val qxq_div : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxq_intnorm : 'a ty -> 'a ty -> 'a ty
val qxq_inv : 'a ty -> 'a ty -> 'a ty
val qxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qxq_norm : 'a ty -> 'a ty -> 'a ty
val qxq_sqr : 'a ty -> 'a ty -> 'a ty
val rg_is_fp : 'a ty -> 'a ty Ctypes_static.ptr -> int

val rg_is_fpxq :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> int

val rg_to_fp : 'a ty -> 'a ty -> 'a ty
val rg_to_fpxq : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgc_to_fpc : 'a ty -> 'a ty -> 'a ty
val rgc_to_fqc : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgm_is_fpm : 'a ty -> 'a ty Ctypes_static.ptr -> int
val rgm_to_flm : 'a ty -> pari_ulong -> 'a ty
val rgm_to_fpm : 'a ty -> 'a ty -> 'a ty
val rgm_to_fqm : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgv_is_fpv : 'a ty -> 'a ty Ctypes_static.ptr -> int
val rgv_to_flv : 'a ty -> pari_ulong -> 'a ty
val rgv_to_fpv : 'a ty -> 'a ty -> 'a ty
val rgx_is_fpx : 'a ty -> 'a ty Ctypes_static.ptr -> int
val rgx_to_fpx : 'a ty -> 'a ty -> 'a ty

val rgx_is_fpxqx :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> int

val rgx_to_fpxqx : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgx_to_fqx : 'a ty -> 'a ty -> 'a ty -> 'a ty

val z_incremental_crt :
  'a ty Ctypes_static.ptr ->
  pari_ulong ->
  'a ty Ctypes_static.ptr ->
  pari_ulong ->
  int

val z_init_crt : pari_ulong -> pari_ulong -> 'a ty

val zm_incremental_crt :
  'a ty Ctypes_static.ptr ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  pari_ulong ->
  int

val zm_init_crt : 'a ty -> pari_ulong -> 'a ty
val zx_zxy_resultant : 'a ty -> 'a ty -> 'a ty

val zx_zxy_rnfequation :
  'a ty -> 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty

val zx_disc : 'a ty -> 'a ty
val zx_gcd : 'a ty -> 'a ty -> 'a ty
val zx_gcd_all : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val zx_incremental_crt :
  'a ty Ctypes_static.ptr ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  pari_ulong ->
  int

val zx_init_crt : 'a ty -> pari_ulong -> Signed.long -> 'a ty
val zx_is_squarefree : 'a ty -> int
val zx_radical : 'a ty -> 'a ty
val zx_resultant : 'a ty -> 'a ty -> 'a ty

val zxm_incremental_crt :
  'a ty Ctypes_static.ptr ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  pari_ulong ->
  int

val zxm_init_crt : 'a ty -> Signed.long -> pari_ulong -> 'a ty
val zxq_minpoly : 'a ty -> 'a ty -> Signed.long -> pari_ulong -> 'a ty
val zxq_charpoly : 'a ty -> 'a ty -> Signed.long -> 'a ty
val characteristic : 'a ty -> 'a ty
val ffnbirred : 'a ty -> Signed.long -> 'a ty
val ffnbirred0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val ffsumnbirred : 'a ty -> Signed.long -> 'a ty

val get_fq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  'a ty ->
  'a ty ->
  bb_field Ctypes.structure Ctypes_static.ptr

val init_flxq : pari_ulong -> Signed.long -> Signed.long -> 'a ty
val init_fq : 'a ty -> Signed.long -> Signed.long -> 'a ty
val residual_characteristic : 'a ty -> 'a ty
val fp_modinv_to_j : 'a ty -> Signed.long -> 'a ty -> 'a ty

val fp_polmodular_evalx :
  Signed.long -> Signed.long -> 'a ty -> 'a ty -> Signed.long -> int -> 'a ty

val check_modinv : Signed.long -> unit
val disc_best_modinv : Signed.long -> Signed.long
val modinv_height_factor : Signed.long -> Signed.long
val modinv_good_disc : Signed.long -> Signed.long -> int
val modinv_good_prime : Signed.long -> Signed.long -> int
val modinv_is_weber : Signed.long -> int
val modinv_is_double_eta : Signed.long -> int
val bpsw_isprime : 'a ty -> Signed.long
val bpsw_psp : 'a ty -> Signed.long
val addprimes : 'a ty -> 'a ty
val check_ecppcert : 'a ty -> Signed.long
val gisprime : 'a ty -> Signed.long -> 'a ty
val gispseudoprime : 'a ty -> Signed.long -> 'a ty
val gprimepi_upper_bound : 'a ty -> 'a ty
val gprimepi_lower_bound : 'a ty -> 'a ty
val isprime : 'a ty -> Signed.long
val ispseudoprime : 'a ty -> Signed.long -> Signed.long
val millerrabin : 'a ty -> Signed.long -> Signed.long
val prime : Signed.long -> 'a ty
val primecert : 'a ty -> Signed.long -> 'a ty
val primecert0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val primecertexport : 'a ty -> Signed.long -> 'a ty
val primecertisvalid : 'a ty -> Signed.long
val primepi : 'a ty -> 'a ty
val primepi_upper_bound : float -> float
val primepi_lower_bound : float -> float
val primes : Signed.long -> 'a ty
val primes_interval : 'a ty -> 'a ty -> 'a ty
val primes_interval_zv : pari_ulong -> pari_ulong -> 'a ty
val primes_upto_zv : pari_ulong -> 'a ty
val primes0 : 'a ty -> 'a ty
val primes_zv : Signed.long -> 'a ty
val randomprime0 : 'a ty -> 'a ty -> 'a ty
val removeprimes : 'a ty -> 'a ty
val uis2psp : pari_ulong -> int
val uispsp : pari_ulong -> pari_ulong -> int
val uislucaspsp : pari_ulong -> int
val uisprime : pari_ulong -> int
val uisprime_101 : pari_ulong -> int
val uisprime_661 : pari_ulong -> int
val uprime : Signed.long -> pari_ulong
val uprimepi : pari_ulong -> pari_ulong
val qfauto : 'a ty -> 'a ty -> 'a ty
val qfauto0 : 'a ty -> 'a ty -> 'a ty
val qfautoexport : 'a ty -> Signed.long -> 'a ty
val qfisom : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val qfisom0 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val qfisominit : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qfisominit0 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qforbits : 'a ty -> 'a ty -> 'a ty
val qfminimize : 'a ty -> 'a ty
val qfparam : 'a ty -> 'a ty -> Signed.long -> 'a ty
val qfsolve : 'a ty -> 'a ty
val z_isfundamental : 'a ty -> Signed.long
val classno : 'a ty -> 'a ty
val classno2 : 'a ty -> 'a ty
val hclassnof_fact : 'a ty -> 'a ty -> 'a ty -> 'a ty
val hclassno : 'a ty -> 'a ty
val hclassno6 : 'a ty -> 'a ty
val isfundamental : 'a ty -> Signed.long
val qfb_equal1 : 'a ty -> int
val qfbclassno0 : 'a ty -> Signed.long -> 'a ty
val qfi_shanks : 'a ty -> 'a ty -> Signed.long -> 'a ty
val qfi_log : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qfi_order : 'a ty -> 'a ty -> 'a ty
val quadclassnof : 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val quadclassnof_fact : 'a ty -> 'a ty -> 'a ty -> 'a ty
val quaddisc : 'a ty -> 'a ty
val quadregulator : 'a ty -> Signed.long -> 'a ty
val quadunit : 'a ty -> 'a ty
val quadunit0 : 'a ty -> Signed.long -> 'a ty
val quadunitindex : 'a ty -> 'a ty -> 'a ty
val quadunitnorm : 'a ty -> Signed.long
val sisfundamental : Signed.long -> Signed.long
val uhclassnof_fact : 'a ty -> Signed.long -> Signed.long
val unegisfundamental : pari_ulong -> Signed.long
val unegquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val uposisfundamental : pari_ulong -> Signed.long
val uposquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val uquadclassnof_fact :
  pari_ulong -> Signed.long -> 'a ty -> 'a ty -> pari_ulong

val zn_quad_roots : 'a ty -> 'a ty -> 'a ty -> 'a ty
val getrand : unit -> 'a ty
val pari_rand : unit -> pari_ulong
val randomr : Signed.long -> 'a ty
val random_f2x : Signed.long -> Signed.long -> 'a ty
val random_fl : pari_ulong -> pari_ulong
val random_bits : Signed.long -> Signed.long
val random_zv : Signed.long -> 'a ty
val setrand : 'a ty -> unit
val hyperellratpoints : 'a ty -> 'a ty -> Signed.long -> 'a ty
val qx_complex_roots : 'a ty -> Signed.long -> 'a ty
val fft : 'a ty -> 'a ty -> 'a ty
val fftinv : 'a ty -> 'a ty -> 'a ty
val cleanroots : 'a ty -> Signed.long -> 'a ty
val fujiwara_bound : 'a ty -> float
val fujiwara_bound_real : 'a ty -> Signed.long -> float
val isrealappr : 'a ty -> Signed.long -> int
val roots : 'a ty -> Signed.long -> 'a ty
val realroots : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zx_graeffe : 'a ty -> 'a ty
val zx_realroots_irred : 'a ty -> Signed.long -> 'a ty
val zx_sturm : 'a ty -> Signed.long
val zx_sturm_irred : 'a ty -> Signed.long
val zx_sturmpart : 'a ty -> 'a ty -> Signed.long
val zx_uspensky : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val factor_aurifeuille : 'a ty -> Signed.long -> 'a ty
val factor_aurifeuille_prime : 'a ty -> Signed.long -> 'a ty
val galoissubcyclo : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val znsubgroupgenerators : 'a ty -> Signed.long -> 'a ty
val subgrouplist : 'a ty -> 'a ty -> 'a ty

val forsubgroup :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  unit

val abmap_kernel : 'a ty -> 'a ty
val abmap_subgroup_image : 'a ty -> 'a ty -> 'a ty
val bnrl1 : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val bnrrootnumber : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val bnrstark : 'a ty -> 'a ty -> Signed.long -> 'a ty
val cyc2elts : 'a ty -> 'a ty
val qfbforms : 'a ty -> 'a ty
val quadhilbert : 'a ty -> Signed.long -> 'a ty
val quadray : 'a ty -> 'a ty -> Signed.long -> 'a ty
val chartogenstr : char -> 'a ty
val pari_strdup : string -> string
val pari_strndup : string -> Signed.long -> string
val stack_strcat : string -> string -> string
val stack_strdup : string -> string
val pari_strchr : 'a ty -> 'a ty
val strjoin : 'a ty -> 'a ty -> 'a ty
val strntogenstr : string -> Signed.long -> 'a ty
val strsplit : 'a ty -> 'a ty -> 'a ty
val strtogenstr : string -> 'a ty
val type_name : Signed.long -> string

val asympnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val asympnumraw :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  Signed.long ->
  'a ty ->
  Signed.long ->
  'a ty

val derivnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val derivnumk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val derivfun :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val derivfunk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val forvec_init :
  forvec_t Ctypes.structure Ctypes_static.ptr -> 'a ty -> Signed.long -> int

val forvec_next : forvec_t Ctypes.structure Ctypes_static.ptr -> 'a ty

val laurentseries :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  'a ty

val limitnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val prodeuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val prodinf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val prodinf1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val solvestep :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  'a ty

val sumalt :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val sumalt2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val sumpos :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val sumpos2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val suminf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val suminf_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val sumdivmultexpr :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty

val zbrent :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  'a ty ->
  Signed.long ->
  'a ty

val bnfisintnorm : 'a ty -> 'a ty -> 'a ty
val bnfisintnormabs : 'a ty -> 'a ty -> 'a ty
val ideals_by_norm : 'a ty -> 'a ty -> 'a ty
val thue : 'a ty -> 'a ty -> 'a ty -> 'a ty
val thueinit : 'a ty -> Signed.long -> Signed.long -> 'a ty
val pi2n : Signed.long -> Signed.long -> 'a ty
val pii2 : Signed.long -> 'a ty
val pii2n : Signed.long -> Signed.long -> 'a ty
val qp_exp : 'a ty -> 'a ty
val qp_exp_prec : 'a ty -> Signed.long
val qp_log : 'a ty -> 'a ty
val qp_sqrt : 'a ty -> 'a ty
val qp_sqrtn : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val zn_sqrt : 'a ty -> 'a ty -> 'a ty
val zp_teichmuller : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty
val agm : 'a ty -> 'a ty -> Signed.long -> 'a ty
val constcatalan : Signed.long -> 'a ty
val consteuler : Signed.long -> 'a ty
val constlog2 : Signed.long -> 'a ty
val constpi : Signed.long -> 'a ty
val cxexpm1 : 'a ty -> Signed.long -> 'a ty
val expir : 'a ty -> 'a ty
val exp1r_abs : 'a ty -> 'a ty
val gcos : 'a ty -> Signed.long -> 'a ty
val gcotan : 'a ty -> Signed.long -> 'a ty
val gcotanh : 'a ty -> Signed.long -> 'a ty
val gexp : 'a ty -> Signed.long -> 'a ty
val gexpm1 : 'a ty -> Signed.long -> 'a ty
val glog : 'a ty -> Signed.long -> 'a ty
val glog1p : 'a ty -> Signed.long -> 'a ty
val gpow : 'a ty -> 'a ty -> Signed.long -> 'a ty
val gpowers : 'a ty -> Signed.long -> 'a ty
val gpowers0 : 'a ty -> Signed.long -> 'a ty -> 'a ty
val gpowgs : 'a ty -> Signed.long -> 'a ty
val grootsof1 : Signed.long -> Signed.long -> 'a ty
val gsin : 'a ty -> Signed.long -> 'a ty
val gsinc : 'a ty -> Signed.long -> 'a ty

val gsincos :
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  Signed.long ->
  unit

val gsqrpowers : 'a ty -> Signed.long -> 'a ty
val gsqrt : 'a ty -> Signed.long -> 'a ty
val gsqrtn : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> Signed.long -> 'a ty
val gtan : 'a ty -> Signed.long -> 'a ty
val logr_abs : 'a ty -> 'a ty
val mpcos : 'a ty -> 'a ty
val mpeuler : Signed.long -> 'a ty
val mpcatalan : Signed.long -> 'a ty

val mpsincosm1 :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val mpexp : 'a ty -> 'a ty
val mpexpm1 : 'a ty -> 'a ty
val mplog : 'a ty -> 'a ty
val mplog2 : Signed.long -> 'a ty
val mppi : Signed.long -> 'a ty
val mpsin : 'a ty -> 'a ty

val mpsincos :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val pow2pis : 'a ty -> Signed.long -> 'a ty
val powpis : 'a ty -> Signed.long -> 'a ty
val powcx : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val powcx_prec : Signed.long -> 'a ty -> Signed.long -> Signed.long
val powersr : 'a ty -> Signed.long -> 'a ty
val powiu : 'a ty -> pari_ulong -> 'a ty
val powrfrac : 'a ty -> Signed.long -> Signed.long -> 'a ty
val powrs : 'a ty -> Signed.long -> 'a ty
val powrshalf : 'a ty -> Signed.long -> 'a ty
val powru : 'a ty -> pari_ulong -> 'a ty
val powruhalf : 'a ty -> pari_ulong -> 'a ty
val powgi : 'a ty -> 'a ty -> 'a ty
val rootsof1_cx : 'a ty -> Signed.long -> 'a ty
val rootsof1u_cx : pari_ulong -> Signed.long -> 'a ty
val rootsof1q_cx : Signed.long -> Signed.long -> Signed.long -> 'a ty
val rootsof1powinit : Signed.long -> Signed.long -> Signed.long -> 'a ty
val rootsof1pow : 'a ty -> Signed.long -> 'a ty
val serchop : 'a ty -> Signed.long -> 'a ty
val serchop_i : 'a ty -> Signed.long -> 'a ty
val serchop0 : 'a ty -> 'a ty
val sqrtnint : 'a ty -> Signed.long -> 'a ty
val sqrtnr_abs : 'a ty -> Signed.long -> 'a ty
val teich : 'a ty -> 'a ty
val teichmullerinit : Signed.long -> Signed.long -> 'a ty
val teichmuller : 'a ty -> 'a ty -> 'a ty

val trans_eval :
  string ->
  ('a ty -> Signed.long -> 'a ty) Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val trans_evalgen :
  string ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long -> 'a ty)
  Ctypes_static.static_funptr ->
  'a ty ->
  Signed.long ->
  'a ty

val upowuu : pari_ulong -> pari_ulong -> pari_ulong
val upowers : pari_ulong -> Signed.long -> 'a ty
val usqrtn : pari_ulong -> pari_ulong -> pari_ulong
val usqrt : pari_ulong -> pari_ulong
val qp_gamma : 'a ty -> 'a ty
val atanhuu : pari_ulong -> pari_ulong -> Signed.long -> 'a ty
val atanhui : pari_ulong -> 'a ty -> Signed.long -> 'a ty
val gacosh : 'a ty -> Signed.long -> 'a ty
val gacos : 'a ty -> Signed.long -> 'a ty
val garg : 'a ty -> Signed.long -> 'a ty
val gasinh : 'a ty -> Signed.long -> 'a ty
val gasin : 'a ty -> Signed.long -> 'a ty
val gatan : 'a ty -> Signed.long -> 'a ty
val gatanh : 'a ty -> Signed.long -> 'a ty
val gcosh : 'a ty -> Signed.long -> 'a ty
val ggammah : 'a ty -> Signed.long -> 'a ty
val ggamma : 'a ty -> Signed.long -> 'a ty
val ggamma1m1 : 'a ty -> Signed.long -> 'a ty
val glngamma : 'a ty -> Signed.long -> 'a ty
val gpsi : 'a ty -> Signed.long -> 'a ty
val gsinh : 'a ty -> Signed.long -> 'a ty
val gtanh : 'a ty -> Signed.long -> 'a ty
val mpfactr : Signed.long -> Signed.long -> 'a ty

val mpsinhcosh :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val psi1series : Signed.long -> Signed.long -> Signed.long -> 'a ty
val sumformal : 'a ty -> Signed.long -> 'a ty

val rgv_is_arithprog :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> int

val besseljzero : 'a ty -> Signed.long -> Signed.long -> 'a ty
val besselyzero : 'a ty -> Signed.long -> Signed.long -> 'a ty
val constzeta : Signed.long -> Signed.long -> 'a ty
val cxek : 'a ty -> Signed.long -> Signed.long -> 'a ty
val dblmodulus : 'a ty -> float
val dilog : 'a ty -> Signed.long -> 'a ty
val eint1 : 'a ty -> Signed.long -> 'a ty
val expipir : 'a ty -> Signed.long -> 'a ty
val expipic : 'a ty -> Signed.long -> 'a ty
val expixy : 'a ty -> 'a ty -> Signed.long -> 'a ty
val eta : 'a ty -> Signed.long -> 'a ty
val eta0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val gerfc : 'a ty -> Signed.long -> 'a ty
val gpolylog : Signed.long -> 'a ty -> Signed.long -> 'a ty
val gzeta : 'a ty -> Signed.long -> 'a ty
val hbessel1 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val hbessel2 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val hyperu : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val ibessel : 'a ty -> 'a ty -> Signed.long -> 'a ty
val incgam : 'a ty -> 'a ty -> Signed.long -> 'a ty
val incgam0 : 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val incgamc : 'a ty -> 'a ty -> Signed.long -> 'a ty
val jbessel : 'a ty -> 'a ty -> Signed.long -> 'a ty
val jbesselh : 'a ty -> 'a ty -> Signed.long -> 'a ty
val jell : 'a ty -> Signed.long -> 'a ty
val kbessel : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mpeint1 : 'a ty -> 'a ty -> 'a ty
val mpveceint1 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val sumdedekind : 'a ty -> 'a ty -> 'a ty
val sumdedekind_coprime : 'a ty -> 'a ty -> 'a ty
val szeta : Signed.long -> Signed.long -> 'a ty
val theta : 'a ty -> 'a ty -> Signed.long -> 'a ty
val thetanullk : 'a ty -> Signed.long -> Signed.long -> 'a ty
val trueeta : 'a ty -> Signed.long -> 'a ty
val u_sumdedekind_coprime : Signed.long -> Signed.long -> 'a ty
val upper_to_cx : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val veceint1 : 'a ty -> 'a ty -> Signed.long -> 'a ty
val vecthetanullk : 'a ty -> Signed.long -> Signed.long -> 'a ty
val vecthetanullk_tau : 'a ty -> Signed.long -> Signed.long -> 'a ty
val veczeta : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val weber0 : 'a ty -> Signed.long -> Signed.long -> 'a ty
val weberf : 'a ty -> Signed.long -> 'a ty
val weberf1 : 'a ty -> Signed.long -> 'a ty
val weberf2 : 'a ty -> Signed.long -> 'a ty
val ybessel : 'a ty -> 'a ty -> Signed.long -> 'a ty
val sl2_inv_shallow : 'a ty -> 'a ty
val qevproj_apply : 'a ty -> 'a ty -> 'a ty
val qevproj_apply_vecei : 'a ty -> 'a ty -> Signed.long -> 'a ty
val qevproj_down : 'a ty -> 'a ty -> 'a ty
val qevproj_init : 'a ty -> 'a ty
val rgx_act_gl2q : 'a ty -> Signed.long -> 'a ty
val rgx_act_zgl2q : 'a ty -> Signed.long -> 'a ty
val checkms : 'a ty -> unit
val checkmspadic : 'a ty -> unit
val mfnumcusps : 'a ty -> 'a ty
val mfnumcusps_fact : 'a ty -> 'a ty
val mfnumcuspsu_fact : 'a ty -> pari_ulong
val mfnumcuspsu : pari_ulong -> pari_ulong
val msfromcusp : 'a ty -> 'a ty -> 'a ty
val msfromell : 'a ty -> Signed.long -> 'a ty
val msfromhecke : 'a ty -> 'a ty -> 'a ty -> 'a ty
val msdim : 'a ty -> Signed.long
val mseval2_ooq : 'a ty -> 'a ty -> 'a ty -> 'a ty
val msgetlevel : 'a ty -> Signed.long
val msgetsign : 'a ty -> Signed.long
val msgetweight : 'a ty -> Signed.long
val msatkinlehner : 'a ty -> Signed.long -> 'a ty -> 'a ty
val mscuspidal : 'a ty -> Signed.long -> 'a ty
val mseisenstein : 'a ty -> 'a ty
val mseval : 'a ty -> 'a ty -> 'a ty -> 'a ty
val mshecke : 'a ty -> Signed.long -> 'a ty -> 'a ty
val msinit : 'a ty -> 'a ty -> Signed.long -> 'a ty
val msissymbol : 'a ty -> 'a ty -> 'a ty
val mslattice : 'a ty -> 'a ty -> 'a ty
val msomseval : 'a ty -> 'a ty -> 'a ty -> 'a ty

val mspadic_parse_chi :
  'a ty -> 'a ty Ctypes_static.ptr -> 'a ty Ctypes_static.ptr -> unit

val mspadic_unit_eigenvalue :
  'a ty -> Signed.long -> 'a ty -> Signed.long -> 'a ty

val mspadicinit : 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val mspadicl : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mspadicmoments : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mspadicseries : 'a ty -> Signed.long -> 'a ty
val mspathgens : 'a ty -> 'a ty
val mspathlog : 'a ty -> 'a ty -> 'a ty
val msnew : 'a ty -> 'a ty
val mspetersson : 'a ty -> 'a ty -> 'a ty -> 'a ty
val mspolygon : 'a ty -> Signed.long -> 'a ty
val msstar : 'a ty -> 'a ty -> 'a ty
val msqexpansion : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val mssplit : 'a ty -> 'a ty -> Signed.long -> 'a ty
val mstooms : 'a ty -> 'a ty -> 'a ty
val mscosets0 : 'a ty -> 'a ty -> 'a ty

val mscosets :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  'a ty

val msfarey :
  'a ty ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> 'a ty -> Signed.long) Ctypes_static.static_funptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val msfarey0 : 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val checkfarey_i : 'a ty -> int
val zetamult : 'a ty -> Signed.long -> 'a ty
val zetamultdual : 'a ty -> 'a ty
val zetamult_interpolate : 'a ty -> 'a ty -> Signed.long -> 'a ty
val zetamultall : Signed.long -> Signed.long -> Signed.long -> 'a ty
val zetamultconvert : 'a ty -> Signed.long -> 'a ty
val fl_add : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_addmul_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fl_addmulmul_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fl_center : pari_ulong -> pari_ulong -> pari_ulong -> Signed.long
val fl_div : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val fl_double : pari_ulong -> pari_ulong -> pari_ulong

val fl_ellj_pre :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_halve : pari_ulong -> pari_ulong -> pari_ulong
val fl_mul : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_mul_pre :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl_neg : pari_ulong -> pari_ulong -> pari_ulong
val fl_sqr : pari_ulong -> pari_ulong -> pari_ulong
val fl_sqr_pre : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val fl_sub : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val fl_triple : pari_ulong -> pari_ulong -> pari_ulong
val abscmpiu : 'a ty -> pari_ulong -> int
val abscmpui : pari_ulong -> 'a ty -> int
val absequaliu : 'a ty -> pari_ulong -> int
val absi : 'a ty -> 'a ty
val absi_shallow : 'a ty -> 'a ty
val absr : 'a ty -> 'a ty
val absrnz_equal1 : 'a ty -> int
val absrnz_equal2n : 'a ty -> int
val addiiz : 'a ty -> 'a ty -> 'a ty -> unit
val addir : 'a ty -> 'a ty -> 'a ty
val addirz : 'a ty -> 'a ty -> 'a ty -> unit
val addis : 'a ty -> Signed.long -> 'a ty
val addri : 'a ty -> 'a ty -> 'a ty
val addriz : 'a ty -> 'a ty -> 'a ty -> unit
val addrr : 'a ty -> 'a ty -> 'a ty
val addrrz : 'a ty -> 'a ty -> 'a ty -> unit
val addrs : 'a ty -> Signed.long -> 'a ty
val addsi : Signed.long -> 'a ty -> 'a ty
val addsiz : Signed.long -> 'a ty -> 'a ty -> unit
val addsrz : Signed.long -> 'a ty -> 'a ty -> unit
val addss : Signed.long -> Signed.long -> 'a ty
val addssz : Signed.long -> Signed.long -> 'a ty -> unit
val adduu : pari_ulong -> pari_ulong -> 'a ty
val affii : 'a ty -> 'a ty -> unit
val affiz : 'a ty -> 'a ty -> unit
val affrr_fixlg : 'a ty -> 'a ty -> unit
val affsi : Signed.long -> 'a ty -> unit
val affsr : Signed.long -> 'a ty -> unit
val affsz : Signed.long -> 'a ty -> unit
val affui : pari_ulong -> 'a ty -> unit
val affur : pari_ulong -> 'a ty -> unit
val cgetg_block : Signed.long -> Signed.long -> 'a ty
val cgetg_copy : 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val cgeti : Signed.long -> 'a ty
val cgetineg : Signed.long -> 'a ty
val cgetipos : Signed.long -> 'a ty
val cgetr : Signed.long -> 'a ty
val cgetr_block : Signed.long -> 'a ty
val cmpir : 'a ty -> 'a ty -> int
val cmpis : 'a ty -> Signed.long -> int
val cmpiu : 'a ty -> pari_ulong -> int
val cmpri : 'a ty -> 'a ty -> int
val cmprs : 'a ty -> Signed.long -> int
val cmpsi : Signed.long -> 'a ty -> int
val cmpsr : Signed.long -> 'a ty -> int
val cmpss : Signed.long -> Signed.long -> int
val cmpui : pari_ulong -> 'a ty -> int
val cmpuu : pari_ulong -> pari_ulong -> int
val divii : 'a ty -> 'a ty -> 'a ty
val diviiz : 'a ty -> 'a ty -> 'a ty -> unit
val divirz : 'a ty -> 'a ty -> 'a ty -> unit
val divisz : 'a ty -> Signed.long -> 'a ty -> unit
val divriz : 'a ty -> 'a ty -> 'a ty -> unit
val divrrz : 'a ty -> 'a ty -> 'a ty -> unit
val divrsz : 'a ty -> Signed.long -> 'a ty -> unit
val divsi_rem : Signed.long -> 'a ty -> Signed.long Ctypes_static.ptr -> 'a ty
val divsiz : Signed.long -> 'a ty -> 'a ty -> unit
val divsrz : Signed.long -> 'a ty -> 'a ty -> unit
val divss : Signed.long -> Signed.long -> 'a ty

val divss_rem :
  Signed.long -> Signed.long -> Signed.long Ctypes_static.ptr -> 'a ty

val divssz : Signed.long -> Signed.long -> 'a ty -> unit
val dvdii : 'a ty -> 'a ty -> int
val dvdiiz : 'a ty -> 'a ty -> 'a ty -> int
val dvdis : 'a ty -> Signed.long -> int
val dvdisz : 'a ty -> Signed.long -> 'a ty -> int
val dvdiu : 'a ty -> pari_ulong -> int
val dvdiuz : 'a ty -> pari_ulong -> 'a ty -> int
val dvdsi : Signed.long -> 'a ty -> int
val dvdui : pari_ulong -> 'a ty -> int
val dvmdiiz : 'a ty -> 'a ty -> 'a ty -> 'a ty -> unit
val dvmdis : 'a ty -> Signed.long -> 'a ty Ctypes_static.ptr -> 'a ty
val dvmdisz : 'a ty -> Signed.long -> 'a ty -> 'a ty -> unit
val dvmdsbil : Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long
val dvmdsi : Signed.long -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty
val dvmdsiz : Signed.long -> 'a ty -> 'a ty -> 'a ty -> unit
val dvmdss : Signed.long -> Signed.long -> 'a ty Ctypes_static.ptr -> 'a ty
val dvmdssz : Signed.long -> Signed.long -> 'a ty -> 'a ty -> unit
val dvmdubil : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val equalis : 'a ty -> Signed.long -> int
val equalsi : Signed.long -> 'a ty -> int
val equalui : pari_ulong -> 'a ty -> int
val equaliu : 'a ty -> pari_ulong -> int
val absequalui : pari_ulong -> 'a ty -> int
val ceildivuu : pari_ulong -> pari_ulong -> pari_ulong
val evalexpo : Signed.long -> Signed.long
val evallg : Signed.long -> Signed.long
val evalprecp : Signed.long -> Signed.long
val evalvalp : Signed.long -> Signed.long
val evalvalser : Signed.long -> Signed.long
val expi : 'a ty -> Signed.long
val expu : pari_ulong -> Signed.long
val fixlg : 'a ty -> Signed.long -> unit
val fractor : 'a ty -> Signed.long -> 'a ty
val gc_bool : pari_ulong -> int -> int
val gc_const : pari_ulong -> 'a ty -> 'a ty
val gc_double : pari_ulong -> float -> float
val gc_int : pari_ulong -> int -> int
val gc_long : pari_ulong -> Signed.long -> Signed.long
val gc_stoi : pari_ulong -> Signed.long -> 'a ty
val gc_ulong : pari_ulong -> pari_ulong -> pari_ulong
val gc_utoi : pari_ulong -> pari_ulong -> 'a ty
val gc_utoipos : pari_ulong -> pari_ulong -> 'a ty
val gc_null : pari_ulong -> 'a ty
val icopy : 'a ty -> 'a ty
val icopyspec : 'a ty -> Signed.long -> 'a ty
val int_bit : 'a ty -> Signed.long -> pari_ulong
val itor : 'a ty -> Signed.long -> 'a ty
val itos : 'a ty -> Signed.long
val itos_or_0 : 'a ty -> Signed.long
val itou : 'a ty -> pari_ulong
val itou_or_0 : 'a ty -> pari_ulong
val leafcopy : 'a ty -> 'a ty
val maxdd : float -> float -> float
val maxss : Signed.long -> Signed.long -> Signed.long
val maxuu : pari_ulong -> pari_ulong -> Signed.long
val mindd : float -> float -> float
val minss : Signed.long -> Signed.long -> Signed.long
val minuu : pari_ulong -> pari_ulong -> Signed.long
val mod16 : 'a ty -> Signed.long
val mod2 : 'a ty -> Signed.long
val mod2bil : 'a ty -> pari_ulong
val mod32 : 'a ty -> Signed.long
val mod4 : 'a ty -> Signed.long
val mod64 : 'a ty -> Signed.long
val mod8 : 'a ty -> Signed.long
val modis : 'a ty -> Signed.long -> 'a ty
val modisz : 'a ty -> Signed.long -> 'a ty -> unit
val modsi : Signed.long -> 'a ty -> 'a ty
val modsiz : Signed.long -> 'a ty -> 'a ty -> unit
val modss : Signed.long -> Signed.long -> 'a ty
val modssz : Signed.long -> Signed.long -> 'a ty -> unit
val mpabs : 'a ty -> 'a ty
val mpabs_shallow : 'a ty -> 'a ty
val mpadd : 'a ty -> 'a ty -> 'a ty
val mpaddz : 'a ty -> 'a ty -> 'a ty -> unit
val mpaff : 'a ty -> 'a ty -> unit
val mpceil : 'a ty -> 'a ty
val mpcmp : 'a ty -> 'a ty -> int
val mpcopy : 'a ty -> 'a ty
val mpdiv : 'a ty -> 'a ty -> 'a ty
val mpexpo : 'a ty -> Signed.long
val mpfloor : 'a ty -> 'a ty
val mpmul : 'a ty -> 'a ty -> 'a ty
val mpmulz : 'a ty -> 'a ty -> 'a ty -> unit
val mpneg : 'a ty -> 'a ty
val mpodd : 'a ty -> int
val mpround : 'a ty -> 'a ty
val mpsqr : 'a ty -> 'a ty
val mpsub : 'a ty -> 'a ty -> 'a ty
val mpsubz : 'a ty -> 'a ty -> 'a ty -> unit
val mptrunc : 'a ty -> 'a ty
val muliiz : 'a ty -> 'a ty -> 'a ty -> unit
val mulirz : 'a ty -> 'a ty -> 'a ty -> unit
val mulis : 'a ty -> Signed.long -> 'a ty
val muliu : 'a ty -> pari_ulong -> 'a ty
val mulri : 'a ty -> 'a ty -> 'a ty
val mulriz : 'a ty -> 'a ty -> 'a ty -> unit
val mulrrz : 'a ty -> 'a ty -> 'a ty -> unit
val mulrs : 'a ty -> Signed.long -> 'a ty
val mulru : 'a ty -> pari_ulong -> 'a ty
val mulsiz : Signed.long -> 'a ty -> 'a ty -> unit
val mulsrz : Signed.long -> 'a ty -> 'a ty -> unit
val mulssz : Signed.long -> Signed.long -> 'a ty -> unit
val negr : 'a ty -> 'a ty
val new_chunk : int -> 'a ty
val rcopy : 'a ty -> 'a ty
val rdivii : 'a ty -> 'a ty -> Signed.long -> 'a ty
val rdiviiz : 'a ty -> 'a ty -> 'a ty -> unit
val rdivis : 'a ty -> Signed.long -> Signed.long -> 'a ty
val rdivsi : Signed.long -> 'a ty -> Signed.long -> 'a ty
val rdivss : Signed.long -> Signed.long -> Signed.long -> 'a ty
val real2n : Signed.long -> Signed.long -> 'a ty
val real_m2n : Signed.long -> Signed.long -> 'a ty
val real_0 : Signed.long -> 'a ty
val real_0_bit : Signed.long -> 'a ty
val real_1 : Signed.long -> 'a ty
val real_1_bit : Signed.long -> 'a ty
val real_m1 : Signed.long -> 'a ty
val remii : 'a ty -> 'a ty -> 'a ty
val remiiz : 'a ty -> 'a ty -> 'a ty -> unit
val remis : 'a ty -> Signed.long -> 'a ty
val remisz : 'a ty -> Signed.long -> 'a ty -> unit

val remlll_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val remsi : Signed.long -> 'a ty -> 'a ty
val remsiz : Signed.long -> 'a ty -> 'a ty -> unit
val remss : Signed.long -> Signed.long -> 'a ty
val remssz : Signed.long -> Signed.long -> 'a ty -> unit
val rtor : 'a ty -> Signed.long -> 'a ty
val sdivsi : Signed.long -> 'a ty -> Signed.long

val sdivsi_rem :
  Signed.long -> 'a ty -> Signed.long Ctypes_static.ptr -> Signed.long

val sdivss_rem :
  Signed.long -> Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val get_avma : unit -> pari_ulong
val set_avma : pari_ulong -> unit

val uabsdiviu_rem :
  'a ty -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val udivuu_rem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val umodi2n : 'a ty -> Signed.long -> pari_ulong
val setabssign : 'a ty -> unit

val shift_left :
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  unit

val shift_right :
  'a ty ->
  'a ty ->
  Signed.long ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  unit

val shiftl : pari_ulong -> pari_ulong -> pari_ulong
val shiftlr : pari_ulong -> pari_ulong -> pari_ulong
val shiftr : 'a ty -> Signed.long -> 'a ty
val shiftr_inplace : 'a ty -> Signed.long -> unit
val smodis : 'a ty -> Signed.long -> Signed.long
val smodss : Signed.long -> Signed.long -> Signed.long
val stackdummy : pari_ulong -> pari_ulong -> unit
val stack_malloc : int -> string
val stack_malloc_align : int -> Signed.long -> string
val stack_calloc : int -> string
val stack_calloc_align : int -> Signed.long -> string
val subiiz : 'a ty -> 'a ty -> 'a ty -> unit
val subir : 'a ty -> 'a ty -> 'a ty
val subirz : 'a ty -> 'a ty -> 'a ty -> unit
val subis : 'a ty -> Signed.long -> 'a ty
val subisz : 'a ty -> Signed.long -> 'a ty -> unit
val subri : 'a ty -> 'a ty -> 'a ty
val subriz : 'a ty -> 'a ty -> 'a ty -> unit
val subrr : 'a ty -> 'a ty -> 'a ty
val subrrz : 'a ty -> 'a ty -> 'a ty -> unit
val subrs : 'a ty -> Signed.long -> 'a ty
val subrsz : 'a ty -> Signed.long -> 'a ty -> unit
val subsi : Signed.long -> 'a ty -> 'a ty
val subsiz : Signed.long -> 'a ty -> 'a ty -> unit
val subsrz : Signed.long -> 'a ty -> 'a ty -> unit
val subss : Signed.long -> Signed.long -> 'a ty
val subssz : Signed.long -> Signed.long -> 'a ty -> unit
val subuu : pari_ulong -> pari_ulong -> 'a ty
val togglesign : 'a ty -> unit
val togglesign_safe : 'a ty Ctypes_static.ptr -> unit
val affectsign : 'a ty -> 'a ty -> unit
val affectsign_safe : 'a ty -> 'a ty Ctypes_static.ptr -> unit
val truedivii : 'a ty -> 'a ty -> 'a ty
val truedivis : 'a ty -> Signed.long -> 'a ty
val truedivsi : Signed.long -> 'a ty -> 'a ty

val uabsdivui_rem :
  pari_ulong -> 'a ty -> pari_ulong Ctypes_static.ptr -> pari_ulong

val umodsu : Signed.long -> pari_ulong -> pari_ulong
val umodui : pari_ulong -> 'a ty -> pari_ulong
val ugcdiu : 'a ty -> pari_ulong -> pari_ulong
val ugcdui : pari_ulong -> 'a ty -> pari_ulong
val umuluu_le : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val umuluu_or_0 : pari_ulong -> pari_ulong -> pari_ulong
val utoi : pari_ulong -> 'a ty
val utoineg : pari_ulong -> 'a ty
val utoipos : pari_ulong -> 'a ty
val utor : pari_ulong -> Signed.long -> 'a ty
val uutoi : pari_ulong -> pari_ulong -> 'a ty
val uutoineg : pari_ulong -> pari_ulong -> 'a ty
val vali : 'a ty -> Signed.long
val varncmp : Signed.long -> Signed.long -> int
val varnmax : Signed.long -> Signed.long -> Signed.long
val varnmin : Signed.long -> Signed.long -> Signed.long
val pari_err_component : string -> string -> 'a ty -> 'a ty -> unit
val pari_err_dim : string -> unit
val pari_err_domain : string -> string -> string -> 'a ty -> 'a ty -> unit
val pari_err_file : string -> string -> unit
val pari_err_filedesc : string -> Signed.long -> unit
val pari_err_flag : string -> unit
val pari_err_impl : string -> unit
val pari_err_inv : string -> 'a ty -> unit
val pari_err_irredpol : string -> 'a ty -> unit
val pari_err_maxprime : pari_ulong -> unit
val pari_err_modulus : string -> 'a ty -> 'a ty -> unit
val pari_err_op : string -> 'a ty -> 'a ty -> unit
val pari_err_overflow : string -> unit
val pari_err_package : string -> unit
val pari_err_prec : string -> unit
val pari_err_prime : string -> 'a ty -> unit
val pari_err_priority : string -> 'a ty -> string -> Signed.long -> unit
val pari_err_sqrtn : string -> 'a ty -> unit
val pari_err_type : string -> 'a ty -> unit
val pari_err_type2 : string -> 'a ty -> 'a ty -> unit
val pari_err_var : string -> 'a ty -> 'a ty -> unit
val pari_err_roots0 : string -> unit
val mkintmod : 'a ty -> 'a ty -> 'a ty
val mkintmodu : pari_ulong -> pari_ulong -> 'a ty
val mkpolmod : 'a ty -> 'a ty -> 'a ty
val mkfrac : 'a ty -> 'a ty -> 'a ty
val mkfracss : Signed.long -> Signed.long -> 'a ty

val qtoss :
  'a ty ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val sstoq : Signed.long -> Signed.long -> 'a ty
val uutoq : pari_ulong -> pari_ulong -> 'a ty
val mkfraccopy : 'a ty -> 'a ty -> 'a ty
val mkrfrac : 'a ty -> 'a ty -> 'a ty
val mkrfraccopy : 'a ty -> 'a ty -> 'a ty
val gen_i : unit -> 'a ty
val cgetc : Signed.long -> 'a ty
val mkquad : 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkvecsmall : Signed.long -> 'a ty
val mkvecsmall2 : Signed.long -> Signed.long -> 'a ty
val mkvecsmall3 : Signed.long -> Signed.long -> Signed.long -> 'a ty

val mkvecsmall4 :
  Signed.long -> Signed.long -> Signed.long -> Signed.long -> 'a ty

val mkvecsmall5 :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  'a ty

val mkqfb : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkvec3 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkvec4 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkvec5 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkvecs : Signed.long -> 'a ty
val mkvec2s : Signed.long -> Signed.long -> 'a ty
val mkvec3s : Signed.long -> Signed.long -> Signed.long -> 'a ty
val mkvec4s : Signed.long -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val mkveccopy : 'a ty -> 'a ty
val mkvec2copy : 'a ty -> 'a ty -> 'a ty
val mkcol : 'a ty -> 'a ty
val mkcol2 : 'a ty -> 'a ty -> 'a ty
val mkcol3 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkcol4 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkcol5 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkcol6 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkcols : Signed.long -> 'a ty
val mkcol2s : Signed.long -> Signed.long -> 'a ty
val mkcol3s : Signed.long -> Signed.long -> Signed.long -> 'a ty
val mkcol4s : Signed.long -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val mkcolcopy : 'a ty -> 'a ty
val mkmat : 'a ty -> 'a ty
val mkmat2 : 'a ty -> 'a ty -> 'a ty
val mkmat3 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkmat4 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkmat5 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val mkmatcopy : 'a ty -> 'a ty
val mkerr : Signed.long -> 'a ty
val mkoo : unit -> 'a ty
val mkmoo : unit -> 'a ty
val inf_get_sign : 'a ty -> Signed.long
val mkmat22s : Signed.long -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val mkmat22 : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val const_vec : Signed.long -> 'a ty -> 'a ty
val const_col : Signed.long -> 'a ty -> 'a ty
val const_vecsmall : Signed.long -> Signed.long -> 'a ty
val zeropadic : 'a ty -> Signed.long -> 'a ty
val zeropadic_shallow : 'a ty -> Signed.long -> 'a ty
val zeroser : Signed.long -> Signed.long -> 'a ty
val ser_isexactzero : 'a ty -> int
val zeropol : Signed.long -> 'a ty
val zerocol : Signed.long -> 'a ty
val zerovec : Signed.long -> 'a ty
val zeromat : Signed.long -> Signed.long -> 'a ty
val zero_flx : Signed.long -> 'a ty
val zero_flv : Signed.long -> 'a ty
val zero_flm : Signed.long -> Signed.long -> 'a ty
val zero_flm_copy : Signed.long -> Signed.long -> 'a ty
val zero_f2v : Signed.long -> 'a ty
val zero_f2m : Signed.long -> Signed.long -> 'a ty
val zero_f2m_copy : Signed.long -> Signed.long -> 'a ty
val zeromatcopy : Signed.long -> Signed.long -> 'a ty
val col_ei : Signed.long -> Signed.long -> 'a ty
val vec_ei : Signed.long -> Signed.long -> 'a ty
val f2v_ei : Signed.long -> Signed.long -> 'a ty
val vecsmall_ei : Signed.long -> Signed.long -> 'a ty
val rg_col_ei : 'a ty -> Signed.long -> Signed.long -> 'a ty
val shallowcopy : 'a ty -> 'a ty
val vectrunc_init : Signed.long -> 'a ty
val coltrunc_init : Signed.long -> 'a ty
val lg_increase : 'a ty -> unit
val vectrunc_append : 'a ty -> 'a ty -> unit
val vectrunc_append_batch : 'a ty -> 'a ty -> unit
val vecsmalltrunc_init : Signed.long -> 'a ty
val vecsmalltrunc_append : 'a ty -> Signed.long -> unit
val hash_str : string -> pari_ulong
val hash_str_len : string -> Signed.long -> pari_ulong
val vec_shorten : 'a ty -> Signed.long -> 'a ty
val vec_lengthen : 'a ty -> Signed.long -> 'a ty
val vec_append : 'a ty -> 'a ty -> 'a ty
val vec_prepend : 'a ty -> 'a ty -> 'a ty
val vec_setconst : 'a ty -> 'a ty -> 'a ty
val vecsmall_shorten : 'a ty -> Signed.long -> 'a ty
val vecsmall_lengthen : 'a ty -> Signed.long -> 'a ty
val vec_to_vecsmall : 'a ty -> 'a ty
val vecsmall_to_vec : 'a ty -> 'a ty
val vecsmall_to_vec_inplace : 'a ty -> 'a ty
val vecsmall_to_col : 'a ty -> 'a ty
val vecsmall_lexcmp : 'a ty -> 'a ty -> int
val vecsmall_prefixcmp : 'a ty -> 'a ty -> int
val vecsmall_prepend : 'a ty -> Signed.long -> 'a ty
val vecsmall_append : 'a ty -> Signed.long -> 'a ty
val vecsmall_concat : 'a ty -> 'a ty -> 'a ty
val vecsmall_coincidence : 'a ty -> 'a ty -> Signed.long
val vecsmall_isin : 'a ty -> Signed.long -> Signed.long
val vecsmall_pack : 'a ty -> Signed.long -> Signed.long -> Signed.long
val vecsmall_indexmax : 'a ty -> Signed.long
val vecsmall_max : 'a ty -> Signed.long
val vecsmall_indexmin : 'a ty -> Signed.long
val vecsmall_min : 'a ty -> Signed.long
val zv_isscalar : 'a ty -> int
val qv_isscalar : 'a ty -> int
val rgv_isscalar : 'a ty -> int
val rgx_isscalar : 'a ty -> int
val rgx_equal_var : 'a ty -> 'a ty -> Signed.long
val rgx_to_rgv : 'a ty -> Signed.long -> 'a ty
val rgx_is_rational : 'a ty -> int
val rgx_is_zx : 'a ty -> int
val rgx_is_qx : 'a ty -> int
val rgx_is_monomial : 'a ty -> int
val rgv_is_zv : 'a ty -> int
val rgv_is_qv : 'a ty -> int
val rgv_isin_i : 'a ty -> 'a ty -> Signed.long -> Signed.long
val rgv_isin : 'a ty -> 'a ty -> Signed.long
val vecslicepermute : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val rowslicepermute : 'a ty -> 'a ty -> Signed.long -> Signed.long -> 'a ty
val rowslice : 'a ty -> Signed.long -> Signed.long -> 'a ty

val matslice :
  'a ty -> Signed.long -> Signed.long -> Signed.long -> Signed.long -> 'a ty

val rowsplice : 'a ty -> Signed.long -> 'a ty
val vecsplice : 'a ty -> Signed.long -> 'a ty
val rgm_minor : 'a ty -> Signed.long -> Signed.long -> 'a ty
val row : 'a ty -> Signed.long -> 'a ty
val flm_row : 'a ty -> Signed.long -> 'a ty
val rowcopy : 'a ty -> Signed.long -> 'a ty
val row_i : 'a ty -> Signed.long -> Signed.long -> Signed.long -> 'a ty
val vecreverse : 'a ty -> 'a ty
val vecsmall_reverse : 'a ty -> 'a ty
val vecreverse_inplace : 'a ty -> unit
val vecsmallpermute : 'a ty -> 'a ty -> 'a ty
val vecpermute : 'a ty -> 'a ty -> 'a ty
val rowpermute : 'a ty -> 'a ty -> 'a ty
val identity_zv : Signed.long -> 'a ty
val identity_perm : Signed.long -> 'a ty
val cyclic_perm : Signed.long -> Signed.long -> 'a ty
val perm_mul : 'a ty -> 'a ty -> 'a ty
val perm_sqr : 'a ty -> 'a ty
val perm_inv : 'a ty -> 'a ty
val perm_conj : 'a ty -> 'a ty -> 'a ty
val pari_free : unit Ctypes_static.ptr -> unit
val pari_malloc : int -> unit Ctypes_static.ptr
val pari_realloc : unit Ctypes_static.ptr -> int -> unit Ctypes_static.ptr
val pari_realloc_ip : unit Ctypes_static.ptr Ctypes_static.ptr -> int -> unit
val pari_calloc : int -> unit Ctypes_static.ptr
val cgetalloc : int -> Signed.long -> 'a ty
val icopy_avma : 'a ty -> pari_ulong -> 'a ty
val leafcopy_avma : 'a ty -> pari_ulong -> 'a ty
val gerepileuptoleaf : pari_ulong -> 'a ty -> 'a ty
val gerepileuptoint : pari_ulong -> 'a ty -> 'a ty
val gerepileupto : pari_ulong -> 'a ty -> 'a ty
val gerepilecopy : pari_ulong -> 'a ty -> 'a ty
val gunclonenull : 'a ty -> unit
val gunclonenull_deep : 'a ty -> unit

val gerepilemany :
  pari_ulong -> 'a ty Ctypes_static.ptr Ctypes_static.ptr -> int -> unit

val gerepileall : pari_ulong -> int -> unit
val gc_all : pari_ulong -> int -> 'a ty
val gerepilecoeffs : pari_ulong -> 'a ty -> int -> unit
val bin_copy : genbin Ctypes.structure Ctypes_static.ptr -> 'a ty
val genbinbase : genbin Ctypes.structure Ctypes_static.ptr -> 'a ty
val cgiv : 'a ty -> unit
val killblock : 'a ty -> unit
val is_universal_constant : 'a ty -> int
val cxcompotor : 'a ty -> Signed.long -> 'a ty
val cxtofp : 'a ty -> Signed.long -> 'a ty
val cxtoreal : 'a ty -> 'a ty
val gtodouble : 'a ty -> float
val gisdouble : 'a ty -> float Ctypes_static.ptr -> int
val gtos : 'a ty -> Signed.long
val gtou : 'a ty -> pari_ulong
val absfrac : 'a ty -> 'a ty
val absfrac_shallow : 'a ty -> 'a ty
val q_abs : 'a ty -> 'a ty
val q_abs_shallow : 'a ty -> 'a ty
val r_abs_shallow : 'a ty -> 'a ty
val r_abs : 'a ty -> 'a ty
val gtofp : 'a ty -> Signed.long -> 'a ty
val gtomp : 'a ty -> Signed.long -> 'a ty
val rgx_gtofp : 'a ty -> Signed.long -> 'a ty
val rgc_gtofp : 'a ty -> Signed.long -> 'a ty
val rgv_gtofp : 'a ty -> Signed.long -> 'a ty
val rgm_gtofp : 'a ty -> Signed.long -> 'a ty
val rgc_gtomp : 'a ty -> Signed.long -> 'a ty
val rgm_gtomp : 'a ty -> Signed.long -> 'a ty
val rgx_fpnorml2 : 'a ty -> Signed.long -> 'a ty
val rgc_fpnorml2 : 'a ty -> Signed.long -> 'a ty
val rgm_fpnorml2 : 'a ty -> Signed.long -> 'a ty
val affgr : 'a ty -> 'a ty -> unit
val affc_fixlg : 'a ty -> 'a ty -> 'a ty
val trunc_safe : 'a ty -> 'a ty
val ndec2nlong : Signed.long -> Signed.long
val ndec2prec : Signed.long -> Signed.long
val ndec2nbits : Signed.long -> Signed.long
val nbits2nlong : Signed.long -> Signed.long
val nbits2extraprec : Signed.long -> Signed.long
val nbits2prec : Signed.long -> Signed.long
val nbits2lg : Signed.long -> Signed.long
val nchar2nlong : Signed.long -> Signed.long
val prec2nbits : Signed.long -> Signed.long
val bit_accuracy_mul : Signed.long -> float -> float
val prec2nbits_mul : Signed.long -> float -> float
val bit_prec : 'a ty -> Signed.long
val bit_accuracy : Signed.long -> Signed.long
val prec2ndec : Signed.long -> Signed.long
val nbits2ndec : Signed.long -> Signed.long
val precdbl : Signed.long -> Signed.long
val divsbil : Signed.long -> Signed.long
val remsbil : Signed.long -> Signed.long
val fp_red : 'a ty -> 'a ty -> 'a ty
val fp_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_neg : 'a ty -> 'a ty -> 'a ty
val fp_halve : 'a ty -> 'a ty -> 'a ty
val fp_center : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_center_i : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_addmul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_sqr : 'a ty -> 'a ty -> 'a ty
val fp_mulu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val fp_muls : 'a ty -> Signed.long -> 'a ty -> 'a ty
val fp_inv : 'a ty -> 'a ty -> 'a ty
val fp_invsafe : 'a ty -> 'a ty -> 'a ty
val fp_div : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fp_divu : 'a ty -> pari_ulong -> 'a ty -> 'a ty
val flx_mulu : 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val get_f2x_mod : 'a ty -> 'a ty
val get_f2x_var : 'a ty -> Signed.long
val get_f2x_degree : 'a ty -> Signed.long
val get_f2xqx_mod : 'a ty -> 'a ty
val get_f2xqx_var : 'a ty -> Signed.long
val get_f2xqx_degree : 'a ty -> Signed.long
val get_flx_mod : 'a ty -> 'a ty
val get_flx_var : 'a ty -> Signed.long
val get_flx_degree : 'a ty -> Signed.long
val get_flxqx_mod : 'a ty -> 'a ty
val get_flxqx_var : 'a ty -> Signed.long
val get_flxqx_degree : 'a ty -> Signed.long
val get_fpx_mod : 'a ty -> 'a ty
val get_fpx_var : 'a ty -> Signed.long
val get_fpx_degree : 'a ty -> Signed.long
val get_fpxqx_mod : 'a ty -> 'a ty
val get_fpxqx_var : 'a ty -> Signed.long
val get_fpxqx_degree : 'a ty -> Signed.long
val submulii : 'a ty -> 'a ty -> 'a ty -> 'a ty
val mulsubii : 'a ty -> 'a ty -> 'a ty -> 'a ty
val submuliu : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val addmuliu : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val submuliu_inplace : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val addmuliu_inplace : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val lincombii : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val is_const_t : Signed.long -> int
val is_extscalar_t : Signed.long -> int
val is_intreal_t : Signed.long -> int
val is_matvec_t : Signed.long -> int
val is_noncalc_t : Signed.long -> int
val is_qfb_t : Signed.long -> int
val is_rational_t : Signed.long -> int
val is_real_t : Signed.long -> int
val is_recursive_t : Signed.long -> int
val is_scalar_t : Signed.long -> int
val is_vec_t : Signed.long -> int
val qfb_is_qfi : 'a ty -> int
val cbrtr_abs : 'a ty -> 'a ty
val cbrtr : 'a ty -> 'a ty
val sqrtnr : 'a ty -> Signed.long -> 'a ty
val logint : 'a ty -> 'a ty -> Signed.long
val ulogint : pari_ulong -> pari_ulong -> pari_ulong
val ismpzero : 'a ty -> int
val isintzero : 'a ty -> int
val isint1 : 'a ty -> int
val isintm1 : 'a ty -> int
val equali1 : 'a ty -> int
val equalim1 : 'a ty -> int
val is_pm1 : 'a ty -> int
val is_bigint : 'a ty -> int
val odd : Signed.long -> int
val both_odd : Signed.long -> Signed.long -> int
val isonstack : 'a ty -> int
val dbllog2r : 'a ty -> float
val mul_content : 'a ty -> 'a ty -> 'a ty
val inv_content : 'a ty -> 'a ty
val div_content : 'a ty -> 'a ty -> 'a ty
val mul_denom : 'a ty -> 'a ty -> 'a ty
val constant_coeff : 'a ty -> 'a ty
val leading_coeff : 'a ty -> 'a ty
val flx_lead : 'a ty -> pari_ulong
val flx_constant : 'a ty -> pari_ulong
val degpol : 'a ty -> Signed.long
val lgpol : 'a ty -> Signed.long
val lgcols : 'a ty -> Signed.long
val nbrows : 'a ty -> Signed.long
val truecoef : 'a ty -> Signed.long -> 'a ty
val zxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zxq_sqr : 'a ty -> 'a ty -> 'a ty
val rgx_copy : 'a ty -> 'a ty
val rgx_coeff : 'a ty -> Signed.long -> 'a ty
val rgx_renormalize : 'a ty -> 'a ty
val rgx_div : 'a ty -> 'a ty -> 'a ty
val rgxqx_div : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rgxqx_rem : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpx_div : 'a ty -> 'a ty -> 'a ty -> 'a ty
val flx_div : 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flx_div_pre : 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val f2x_div : 'a ty -> 'a ty -> 'a ty
val fpv_fpc_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty
val zero_zx : Signed.long -> 'a ty
val zx_shift : 'a ty -> Signed.long -> 'a ty
val zero_f2x : Signed.long -> 'a ty
val f2x_equal1 : 'a ty -> int
val f2x_equal : 'a ty -> 'a ty -> int
val f2x_copy : 'a ty -> 'a ty
val f2v_copy : 'a ty -> 'a ty
val flv_copy : 'a ty -> 'a ty
val flx_copy : 'a ty -> 'a ty
val vecsmall_copy : 'a ty -> 'a ty
val zx_equal1 : 'a ty -> int
val zx_is_monic : 'a ty -> int
val zx_renormalize : 'a ty -> Signed.long -> 'a ty
val fpx_renormalize : 'a ty -> Signed.long -> 'a ty
val fpxx_renormalize : 'a ty -> Signed.long -> 'a ty
val fpxqx_renormalize : 'a ty -> Signed.long -> 'a ty
val f2x_renormalize : 'a ty -> Signed.long -> 'a ty
val f2xx_shift : 'a ty -> Signed.long -> Signed.long -> 'a ty
val f2v_to_f2x : 'a ty -> Signed.long -> 'a ty
val sturm : 'a ty -> Signed.long
val gval : 'a ty -> Signed.long -> Signed.long
val rgx_shift_inplace_init : Signed.long -> unit
val rgx_shift_inplace : 'a ty -> Signed.long -> 'a ty
val zc_to_zc : 'a ty -> 'a ty
val zv_to_zv : 'a ty -> 'a ty
val zx_to_zv : 'a ty -> Signed.long -> 'a ty
val zv_to_zx : 'a ty -> Signed.long -> 'a ty
val zm_to_zxv : 'a ty -> Signed.long -> 'a ty
val zero_zm : Signed.long -> Signed.long -> 'a ty
val zero_zv : Signed.long -> 'a ty
val zm_transpose : 'a ty -> 'a ty
val zm_copy : 'a ty -> 'a ty
val zv_copy : 'a ty -> 'a ty
val zm_row : 'a ty -> Signed.long -> 'a ty
val zc_hnfrem : 'a ty -> 'a ty -> 'a ty
val zm_hnfrem : 'a ty -> 'a ty -> 'a ty
val zm_lll : 'a ty -> float -> Signed.long -> 'a ty

val rgm_dimensions :
  'a ty ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val rgm_shallowcopy : 'a ty -> 'a ty
val f2m_copy : 'a ty -> 'a ty
val f3m_copy : 'a ty -> 'a ty
val flm_copy : 'a ty -> 'a ty
val zv_dvd : 'a ty -> 'a ty -> int
val zm_zv_mod : 'a ty -> 'a ty -> 'a ty
val zv_zv_mod : 'a ty -> 'a ty -> 'a ty
val vecmodii : 'a ty -> 'a ty -> 'a ty
val vecmoduu : 'a ty -> 'a ty -> 'a ty
val fq_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fq_to_fpxq : 'a ty -> 'a ty -> 'a ty -> 'a ty
val rg_to_fq : 'a ty -> 'a ty -> 'a ty -> 'a ty
val gener_fq_local : 'a ty -> 'a ty -> 'a ty -> 'a ty
val random_fq : 'a ty -> 'a ty -> 'a ty
val fpxqx_div : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val flxqx_div : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxqx_div_pre : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> pari_ulong -> 'a ty
val f2xqx_div : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxy_fq_evaly : 'a ty -> 'a ty -> 'a ty -> 'a ty -> Signed.long -> 'a ty
val fqx_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_neg : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_fp_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_fq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_mulu : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fqx_sqr : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_powu : 'a ty -> pari_ulong -> 'a ty -> 'a ty -> 'a ty
val fqx_halve : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_div : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_get_red : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_rem : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fqx_divrem :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val fqx_div_by_x_x :
  'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty Ctypes_static.ptr -> 'a ty

val fqx_halfgcd : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_gcd : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty

val fqx_extgcd :
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty ->
  'a ty Ctypes_static.ptr ->
  'a ty Ctypes_static.ptr ->
  'a ty

val fqx_normalize : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_deriv : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_integ : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_factor : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_factor_squarefree : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_ddf : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_degfact : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_roots : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqx_to_mod : 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxq_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxq_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxq_div : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxq_inv : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxq_invsafe : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxq_mul : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxq_sqr : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxq_pow : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fqxn_expint : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fqxn_exp : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fqxn_inv : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fqxn_mul : 'a ty -> 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fqxn_sqr : 'a ty -> Signed.long -> 'a ty -> 'a ty -> 'a ty
val fpxq_add : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val fpxq_sub : 'a ty -> 'a ty -> 'a ty -> 'a ty -> 'a ty
val flxq_add : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val flxq_sub : 'a ty -> 'a ty -> 'a ty -> pari_ulong -> 'a ty
val f2x_coeff : 'a ty -> Signed.long -> pari_ulong
val f2x_clear : 'a ty -> Signed.long -> unit
val f2x_set : 'a ty -> Signed.long -> unit
val f2x_flip : 'a ty -> Signed.long -> unit
val f2v_coeff : 'a ty -> Signed.long -> pari_ulong
val f2v_clear : 'a ty -> Signed.long -> unit
val f2v_set : 'a ty -> Signed.long -> unit
val f2v_flip : 'a ty -> Signed.long -> unit
val f2m_coeff : 'a ty -> Signed.long -> Signed.long -> pari_ulong
val f2m_clear : 'a ty -> Signed.long -> Signed.long -> unit
val f2m_set : 'a ty -> Signed.long -> Signed.long -> unit
val f2m_flip : 'a ty -> Signed.long -> Signed.long -> unit
val f3m_coeff : 'a ty -> Signed.long -> Signed.long -> pari_ulong
val f3m_set : 'a ty -> Signed.long -> Signed.long -> pari_ulong -> unit
val matpascal : Signed.long -> 'a ty
val z_issquare : 'a ty -> Signed.long
val z_ispower : 'a ty -> pari_ulong -> Signed.long
val sqrti : 'a ty -> 'a ty
val gaddgs : 'a ty -> Signed.long -> 'a ty
val gcmpgs : 'a ty -> Signed.long -> int
val gequalgs : 'a ty -> Signed.long -> int
val gmaxsg : Signed.long -> 'a ty -> 'a ty
val gminsg : Signed.long -> 'a ty -> 'a ty
val gmulgs : 'a ty -> Signed.long -> 'a ty
val gmulgu : 'a ty -> pari_ulong -> 'a ty
val gsubgs : 'a ty -> Signed.long -> 'a ty
val gdivsg : Signed.long -> 'a ty -> 'a ty
val gmax_shallow : 'a ty -> 'a ty -> 'a ty
val gmin_shallow : 'a ty -> 'a ty -> 'a ty
val cxnorm : 'a ty -> 'a ty
val quadnorm : 'a ty -> 'a ty
val quad_disc : 'a ty -> 'a ty
val qfb_disc3 : 'a ty -> 'a ty -> 'a ty -> 'a ty
val qfb_disc : 'a ty -> 'a ty
val sqrfrac : 'a ty -> 'a ty
val normalize_frac : 'a ty -> unit
val powis : Signed.long -> 'a ty
val mpexpz : 'a ty -> 'a ty -> unit
val mplogz : 'a ty -> 'a ty -> unit
val mpcosz : 'a ty -> 'a ty -> unit
val mpsinz : 'a ty -> 'a ty -> unit
val gnegz : 'a ty -> 'a ty -> unit
val gabsz : 'a ty -> Signed.long -> 'a ty -> unit
val gaddz : 'a ty -> 'a ty -> 'a ty -> unit
val gsubz : 'a ty -> 'a ty -> 'a ty -> unit
val gmulz : 'a ty -> 'a ty -> 'a ty -> unit
val gdivz : 'a ty -> 'a ty -> 'a ty -> unit
val gdiventz : 'a ty -> 'a ty -> 'a ty -> unit
val gmodz : 'a ty -> 'a ty -> 'a ty -> unit
val gmul2nz : 'a ty -> Signed.long -> 'a ty -> unit
val gshiftz : 'a ty -> Signed.long -> 'a ty -> unit
val checkell_i : 'a ty -> int
val modpr_get_pr : 'a ty -> 'a ty
val modpr_get_p : 'a ty -> 'a ty
val modpr_get_t : 'a ty -> 'a ty
val pr_get_p : 'a ty -> 'a ty
val pr_get_gen : 'a ty -> 'a ty
val pr_get_e : 'a ty -> Signed.long
val pr_get_f : 'a ty -> Signed.long
val pr_get_tau : 'a ty -> 'a ty
val pr_is_inert : 'a ty -> int
val pr_norm : 'a ty -> 'a ty
val upr_norm : 'a ty -> pari_ulong
val cyc_get_expo : 'a ty -> 'a ty
val abgrp_get_no : 'a ty -> 'a ty
val abgrp_get_cyc : 'a ty -> 'a ty
val abgrp_get_gen : 'a ty -> 'a ty
val bnf_get_nf : 'a ty -> 'a ty
val bnf_get_clgp : 'a ty -> 'a ty
val bnf_get_no : 'a ty -> 'a ty
val bnf_get_cyc : 'a ty -> 'a ty
val bnf_get_gen : 'a ty -> 'a ty
val bnf_get_reg : 'a ty -> 'a ty
val bnf_get_logfu : 'a ty -> 'a ty
val bnf_get_sunits : 'a ty -> 'a ty
val bnf_get_tuu : 'a ty -> 'a ty
val bnf_get_tun : 'a ty -> Signed.long
val bnf_get_fu_nocheck : 'a ty -> 'a ty
val bnf_get_fu : 'a ty -> 'a ty
val bnr_get_bnf : 'a ty -> 'a ty
val bnr_get_bid : 'a ty -> 'a ty
val bnr_get_mod : 'a ty -> 'a ty
val bnr_get_nf : 'a ty -> 'a ty
val bnr_get_clgp : 'a ty -> 'a ty
val bnr_get_no : 'a ty -> 'a ty
val bnr_get_cyc : 'a ty -> 'a ty
val bnr_get_gen_nocheck : 'a ty -> 'a ty
val bnr_get_gen : 'a ty -> 'a ty
val locs_get_cyc : 'a ty -> 'a ty
val locs_get_lsprk : 'a ty -> 'a ty
val locs_get_lgenfil : 'a ty -> 'a ty
val locs_get_mod : 'a ty -> 'a ty
val locs_get_famod : 'a ty -> 'a ty
val locs_get_m_infty : 'a ty -> 'a ty
val gchar_get_basis : 'a ty -> 'a ty
val gchar_get_bnf : 'a ty -> 'a ty
val gchar_get_nf : 'a ty -> 'a ty
val gchar_get_zm : 'a ty -> 'a ty
val gchar_get_mod : 'a ty -> 'a ty
val gchar_get_modp : 'a ty -> 'a ty
val gchar_get_s : 'a ty -> 'a ty
val gchar_get_dldata : 'a ty -> 'a ty
val gchar_get_sfu : 'a ty -> 'a ty
val gchar_get_cyc : 'a ty -> 'a ty
val gchar_get_hnf : 'a ty -> 'a ty
val gchar_get_u : 'a ty -> 'a ty
val gchar_get_ui : 'a ty -> 'a ty
val gchar_get_m0 : 'a ty -> 'a ty
val gchar_get_u0 : 'a ty -> 'a ty
val gchar_get_r1 : 'a ty -> Signed.long
val gchar_get_r2 : 'a ty -> Signed.long
val gchar_get_loccyc : 'a ty -> 'a ty
val gchar_get_nc : 'a ty -> Signed.long
val gchar_get_ns : 'a ty -> Signed.long
val gchar_get_nm : 'a ty -> Signed.long
val gchar_get_evalprec : 'a ty -> Signed.long
val gchar_get_prec : 'a ty -> Signed.long
val gchar_get_nfprec : 'a ty -> Signed.long
val gchar_set_evalprec : 'a ty -> Signed.long -> unit
val gchar_set_prec : 'a ty -> Signed.long -> unit
val gchar_copy_precs : 'a ty -> 'a ty -> unit
val gchar_set_nfprec : 'a ty -> Signed.long -> unit
val gchar_get_ntors : 'a ty -> Signed.long
val gchar_get_nfree : 'a ty -> Signed.long
val gchar_get_nalg : 'a ty -> Signed.long
val gchar_set_basis : 'a ty -> 'a ty -> unit
val gchar_set_nf : 'a ty -> 'a ty -> unit
val gchar_set_ntors : 'a ty -> Signed.long -> unit
val gchar_set_nfree : 'a ty -> Signed.long -> unit
val gchar_set_nalg : 'a ty -> Signed.long -> unit
val gchar_set_cyc : 'a ty -> 'a ty -> unit
val gchar_set_huui : 'a ty -> 'a ty -> 'a ty -> 'a ty -> unit
val gchar_set_m0 : 'a ty -> 'a ty -> unit
val gchar_set_u0 : 'a ty -> 'a ty -> unit
val bid_get_mod : 'a ty -> 'a ty
val bid_get_ideal : 'a ty -> 'a ty
val bid_get_arch : 'a ty -> 'a ty
val bid_get_grp : 'a ty -> 'a ty
val bid_get_fact : 'a ty -> 'a ty
val bid_get_fact2 : 'a ty -> 'a ty
val bid_get_sprk : 'a ty -> 'a ty
val bid_get_sarch : 'a ty -> 'a ty
val bid_get_archp : 'a ty -> 'a ty
val bid_get_u : 'a ty -> 'a ty
val bid_get_no : 'a ty -> 'a ty
val bid_get_cyc : 'a ty -> 'a ty
val bid_get_gen_nocheck : 'a ty -> 'a ty
val bid_get_gen : 'a ty -> 'a ty
val znstar_get_n : 'a ty -> 'a ty
val znstar_get_fan : 'a ty -> 'a ty
val znstar_get_no : 'a ty -> 'a ty
val znstar_get_cyc : 'a ty -> 'a ty
val znstar_get_gen : 'a ty -> 'a ty
val znstar_get_conreycyc : 'a ty -> 'a ty
val znstar_get_conreygen : 'a ty -> 'a ty
val znstar_get_ui : 'a ty -> 'a ty
val znstar_get_u : 'a ty -> 'a ty
val znstar_get_pe : 'a ty -> 'a ty
val gal_get_pol : 'a ty -> 'a ty
val gal_get_p : 'a ty -> 'a ty
val gal_get_e : 'a ty -> 'a ty
val gal_get_mod : 'a ty -> 'a ty
val gal_get_roots : 'a ty -> 'a ty
val gal_get_invvdm : 'a ty -> 'a ty
val gal_get_den : 'a ty -> 'a ty
val gal_get_group : 'a ty -> 'a ty
val gal_get_gen : 'a ty -> 'a ty
val gal_get_orders : 'a ty -> 'a ty
val rnf_get_degree : 'a ty -> Signed.long
val rnf_get_nfdegree : 'a ty -> Signed.long
val rnf_get_absdegree : 'a ty -> Signed.long
val rnf_get_idealdisc : 'a ty -> 'a ty
val rnf_get_k : 'a ty -> 'a ty
val rnf_get_alpha : 'a ty -> 'a ty
val rnf_get_nf : 'a ty -> 'a ty
val rnf_get_nfzk : 'a ty -> 'a ty
val rnf_get_polabs : 'a ty -> 'a ty
val rnf_get_pol : 'a ty -> 'a ty
val rnf_get_disc : 'a ty -> 'a ty
val rnf_get_index : 'a ty -> 'a ty
val rnf_get_ramified_primes : 'a ty -> 'a ty
val rnf_get_varn : 'a ty -> Signed.long
val rnf_get_nfpol : 'a ty -> 'a ty
val rnf_get_nfvarn : 'a ty -> Signed.long
val rnf_get_zk : 'a ty -> 'a ty
val rnf_get_map : 'a ty -> 'a ty
val rnf_get_invzk : 'a ty -> 'a ty
val idealred : 'a ty -> 'a ty -> 'a ty
val idealchineseinit : 'a ty -> 'a ty -> 'a ty
val closure_arity : 'a ty -> Signed.long
val closure_is_variadic : 'a ty -> Signed.long
val closure_codestr : 'a ty -> string
val closure_get_code : 'a ty -> 'a ty
val closure_get_oper : 'a ty -> 'a ty
val closure_get_data : 'a ty -> 'a ty
val closure_get_dbg : 'a ty -> 'a ty
val closure_get_text : 'a ty -> 'a ty
val closure_get_frame : 'a ty -> 'a ty
val err_get_num : 'a ty -> Signed.long
val err_get_compo : 'a ty -> Signed.long -> 'a ty
val pari_err_bug : string -> unit
val pari_err_constpol : string -> unit
val pari_err_coprime : string -> 'a ty -> 'a ty -> unit
val with_stack_clean : (unit -> 'a ty) -> 'a ty
val with_stack_clean_opt : (unit -> 'a ty option) -> 'a ty option

val with_stack_clean6 :
  ?av:pari_ulong ->
  (unit -> 'k1 ty * 'k2 ty * 'k3 ty * 'k4 ty * 'k5 ty * 'k6 ty) ->
  'k1 ty * 'k2 ty * 'k3 ty * 'k4 ty * 'k5 ty * 'k6 ty

val gentobytes : 'a ty -> bytes
