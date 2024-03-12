type pari_ulong = Unsigned.ULong.t

val pari_ulong : pari_ulong Ctypes.typ

type ('kind, 'structure) typ

val t : ('kind, 'structure) typ Ctypes.typ

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
type elliptic_curve = private Elliptic_curve

val factor :
  ('kind, unique_factorization_domain) typ ->
  (('kind, 'structure) typ * int) array

module rec Complex : sig
  type t = (complex, field) typ

  val inv : t -> t
  val add : t -> t -> t
  val create : re:Real.t -> im:Real.t -> t
  val to_string : t -> string
end

and Real : sig
  type t = (real, field) typ

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
  type t = (rational, field) typ

  val inj_ring : t -> (rational, ring) typ
  val inj_real : t -> Real.t
  val inj_complex : t -> Complex.t
  val shift : t -> int -> t
end

and Integer : sig
  type t = (integer, ring) typ

  val inj_rat : t -> Rational.t
  val inj_real : t -> Real.t
  val inj_complex : t -> Complex.t

  val inj_unique_factorization_domain :
    t -> (integer, unique_factorization_domain) typ

  val to_integer : (integer, _) typ -> t
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
  mul : ('a, group) typ -> ('a, group) typ -> ('a, group) typ;
  pow : ('a, group) typ -> Integer.t -> ('a, group) typ;
  rand : unit -> ('a, group) typ;
  hash : ('a, group) typ -> Unsigned.ULong.t;
  equal : ('a, group) typ -> ('a, group) typ -> bool;
  equal_identity : ('a, group) typ -> bool;
  bb_group : bb_group Ctypes.structure option;
}

module Set : sig
  type 'a t constraint 'a = ('b, 'c) typ

  val length : 'a t -> Signed.Long.t
  val search : 'a t -> 'a -> Signed.Long.t -> Signed.Long.t
end

module Vector : sig
  type ('a, 'b) t constraint 'a = ('c, 'd) typ constraint 'b = [< `COL | `ROW ]

  val length : ('a, 'b) t -> int
  val of_array : 'a array -> ('a, [ `ROW ]) t

  val array_map :
    f:('a -> ('b, 'c) typ) -> 'a array -> (('b, 'c) typ, [ `ROW ]) t

  val init : int -> f:(int -> ('a, 'b) typ) -> (('a, 'b) typ, [ `ROW ]) t
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
  type 'a t constraint 'a = ('b, 'c) typ

  val dimensions : 'a t -> int * int
  val id : int -> Integer.t t
  val inv : 'a t -> 'a t
  val mul : 'a t -> 'a t -> 'a t
  val lll : 'a t -> 'a t
  val ( .%[] ) : 'a t -> int -> ('a, [ `COL ]) Vector.t
  val ( .%[]<- ) : 'a t -> int -> ('a, [ `COL ]) Vector.t -> unit
  val ( .%[;..] ) : 'a t -> int array -> 'a t
  val ( .%[;..]<- ) : 'a t -> int array -> 'a -> unit
  val inj : 'a t -> inj:('a -> 'b) -> 'b t
end

module rec Polynomial : sig
  type 'a t = ('a polynomial, ring) typ constraint 'a = ('b, ring) typ

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
          (Integer.of_int 1);
          (Integer.of_int (-111));
          (Integer.of_int 6064);
          (Integer.of_int (-189804));
        |];;
      val q : Integer.t Polynomial.t = <abstr>
      # Polynomial.to_string q;;
      - : string = "x^3 - 111*x^2 + 6064*x - 189804"
      # let qmin = Polynomial.minimal q;;
      val qmin : ('a, ring) typ Polynomial.t = <abstr>
      # Polynomial.to_string qmin;;
      - : string = "x^3 - x^2 - 60*x - 364"
      # Number_field.(are_isomorphic (create q) (create qmin));
      - : bool = true
      ]} *)

  val ( .%[] ) : 'a t -> int -> 'a

  val roots_ff :
    (finite_field, ring) typ t -> ((finite_field, field) typ, [ `ROW ]) Vector.t

  val fold_left : f:('b -> 'a -> 'a) -> acc:'a -> 'b t -> 'a

  val fold_left2 :
    f:('a -> 'b -> ('c, 'd) typ -> ('c, 'd) typ) ->
    acc:('c, 'd) typ ->
    'a t ->
    'b t ->
    ('c, 'd) typ

  val fold_left2_vec :
    f:('a -> 'b -> ('c, 'd) typ -> ('c, 'd) typ) ->
    acc:('c, 'd) typ ->
    'a t ->
    ('b, _) Vector.t ->
    ('c, 'd) typ

  val inj_base_ring : inj:('a -> 'b) -> 'a t -> 'b t
end

and Fp : sig
  type t = Integer.t

  val add : t -> t -> modulo:t -> t
  val pow : t -> exponent:t -> modulo:t -> t
end

and Finite_field : sig
  type t = (finite_field, field) typ

  val inj_ring : t -> (finite_field, ring) typ
  val inj_field : (finite_field, ring) typ -> t
  val generator : order:Integer.t -> t
  val prime_field_element : Integer.t -> p:Integer.t -> t
  val inj_prime_field : t -> Fp.t option
  val finite_field_element : Integer.t array -> t -> t

  val create : p:int -> degree:int -> (finite_field, ring) typ Polynomial.t
  (** [create p degree] returns a monic irreducible polynomial of the
      given [degree] over F_p[X]. *)

  val generator_from_irreducible_polynomial :
    (finite_field, ring) typ Polynomial.t -> t

  val residue_class : t -> (finite_field, ring) typ Polynomial.t
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
    (finite_field, field) typ ->
    [< `Degree of int | `Quotient of (finite_field, ring) typ Polynomial.t ] ->
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
  type t = (integer_mod, ring) typ

  val inj_group : t -> (integer_mod, group) typ
  val create : Integer.t -> modulo:Integer.t -> t

  val create_assume_prime_modulus :
    Integer.t -> modulo:Integer.t -> (integer_mod, field) typ

  val lift : (integer_mod, _) typ -> Integer.t
  val inverse : (integer_mod, 'a) typ -> (integer_mod, 'a) typ option

  val mul :
    (integer_mod, 'a) typ -> (integer_mod, 'a) typ -> (integer_mod, 'a) typ

  val pow : (integer_mod, 'a) typ -> Integer.t -> (integer_mod, 'a) typ
  val chinese : (t, [ `ROW ]) Vector.t -> t
  val to_string : (integer_mod, _) typ -> string
  val get_modulo : (integer_mod, _) typ -> Integer.t
  val order : (integer_mod, _) typ -> Integer.t

  val log :
    base:(integer_mod, _) typ -> (integer_mod, _) typ -> Integer.t option
end

module Number_field : sig
  type t
  type elt = (number_field, field) typ

  val create : (rational, ring) typ Polynomial.t -> t
  (** [create p] returns the number field Q(X)/(p) for a monic
      irreducible polynomial [p] over the field Q of the rationals. *)

  val are_isomorphic : t -> t -> bool
  (** [are_isomorphic a b] returns true if and only if number
      fields [a] and [b] are isomorphic.

      {@ocaml[
      # let q =
        Polynomial.create
          [|
            Integer.of_int 1;
            Integer.of_int (-111);
            Integer.of_int 6064;
            Integer.of_int (-189804);
          |];;
      val q : Integer.t Polynomial.t = <abstr>
      # let zero = Polynomial.create [| Integer.of_int 0 |];;
      val zero : Integer.t Polynomial.t = <abstr>
      # let qq = Polynomial.create [| q; q; zero; zero |];;
      val qq : Integer.t Polynomial.t Polynomial.t = <abstr>
      # Polynomial.to_string qq;;
      - : string =
      "(x^3 - 111*x^2 + 6064*x - 189804)*y^3 + (x^3 - 111*x^2 + 6064*x - 189804)*y^2"
      # Polynomial.is_irreducible q;;
      - : bool = true
      # let qmin : Integer.t Polynomial.t = Polynomial.minimal q;;
      val qmin : Integer.t Polynomial.t = <abstr>
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
        (Polynomial.create [| Integer.of_int 1; Integer.of_int 0; Integer.of_int 1 |]);;
      val gaussian_integers : Number_field.t = <abstr>
      # let a = Number_field.elt [| Integer.(inj_rat (of_int 6)); Integer.(inj_rat (of_int 8)) |];;
      val a : Number_field.elt = <abstr>
      # let b = Number_field.elt [| Integer.(inj_rat (of_int 1)); Integer.(inj_rat (of_int 5)) |];;
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

  module Infix : sig
    val ( = ) : elt -> elt -> bool
  end
end

module Elliptic_curve : sig
  type 'a t constraint 'a = ('b, field) typ
  type 'a elt = (elliptic_curve, group) typ constraint 'a = ('b, field) typ

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

  val l_division_polynomial :
    ('a, field) typ t -> l:Signed.Long.t -> ('a, ring) typ Polynomial.t
  (** {@ocaml[
      # let g = Finite_field.generator ~order:(Integer.of_int 625) (* 5^4 *);;
      val g : Finite_field.t = <abstr>
      # let ell = Option.get (Elliptic_curve.create ~a6:(Finite_field.pow g (Integer.of_int 6)) ());;
      val ell : Finite_field.t Elliptic_curve.t = <abstr>
      # let pdiv7 = (Elliptic_curve.l_division_polynomial ell ~l:(Signed.Long.of_int 7));;
      val pdiv7 : (finite_field, ring) typ Polynomial.t = <abstr>
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

val forprime_t_bb :
  (('kind, 'structure) typ, forprime_t Ctypes.structure) Ctypes.field

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

val forprime_t_pp :
  (('kind, 'structure) typ, forprime_t Ctypes.structure) Ctypes.field

val forcomposite_t : forcomposite_t Ctypes.structure Ctypes.typ
val forcomposite_t_first : (int, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_b :
  (('kind, 'structure) typ, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_n :
  (('kind, 'structure) typ, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_p :
  (('kind, 'structure) typ, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_T :
  (forprime_t Ctypes.structure, forcomposite_t Ctypes.structure) Ctypes.field

val forvec_t : forvec_t Ctypes.structure Ctypes.typ
val forvec_t_first : (Signed.long, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_a :
  ( ('kind, 'structure) typ Ctypes_static.ptr,
    forvec_t Ctypes.structure )
  Ctypes.field

val forvec_t_m :
  ( ('kind, 'structure) typ Ctypes_static.ptr,
    forvec_t Ctypes.structure )
  Ctypes.field

val forvec_t_M :
  ( ('kind, 'structure) typ Ctypes_static.ptr,
    forvec_t Ctypes.structure )
  Ctypes.field

val forvec_t_n : (Signed.long, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_next :
  ( (forvec_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ)
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

val forpart_t_v :
  (('kind, 'structure) typ, forpart_t Ctypes.structure) Ctypes.field

val forperm_t : forperm_t Ctypes.structure Ctypes.typ
val forperm_t_k : (Signed.long, forperm_t Ctypes.structure) Ctypes.field
val forperm_t_first : (Signed.long, forperm_t Ctypes.structure) Ctypes.field

val forperm_t_v :
  (('kind, 'structure) typ, forperm_t Ctypes.structure) Ctypes.field

val forsubset_t : forsubset_t Ctypes.structure Ctypes.typ
val forsubset_t_n : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_k : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_all : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_first : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field

val forsubset_t_v :
  (('kind, 'structure) typ, forsubset_t Ctypes.structure) Ctypes.field

val pari_plot : pari_plot Ctypes.structure Ctypes.typ

val pari_plot_draw :
  ( (pari_plot Ctypes.structure Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
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
val genbin_x : (('kind, 'structure) typ, genbin Ctypes.structure) Ctypes.field

val genbin_base :
  (('kind, 'structure) typ, genbin Ctypes.structure) Ctypes.field

val genbin_rebase :
  ( (('kind, 'structure) typ -> Signed.long -> unit) Ctypes_static.static_funptr,
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
  (('kind, 'structure) typ, pari_parsestate Ctypes.structure) Ctypes.field

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
  (('kind, 'structure) typ, pari_global_state Ctypes.structure) Ctypes.field

val pari_global_state_seadata :
  (('kind, 'structure) typ, pari_global_state Ctypes.structure) Ctypes.field

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

val pari_thread_data :
  (('kind, 'structure) typ, pari_thread Ctypes.structure) Ctypes.field

val mt_state : mt_state Ctypes.structure Ctypes.typ

val mt_state_worker :
  (('kind, 'structure) typ, mt_state Ctypes.structure) Ctypes.field

val mt_state_pending :
  (('kind, 'structure) typ, mt_state Ctypes.structure) Ctypes.field

val mt_state_workid : (Signed.long, mt_state Ctypes.structure) Ctypes.field
val pari_mt : pari_mt Ctypes.structure Ctypes.typ

val pari_mt_mt :
  (mt_state Ctypes.structure, pari_mt Ctypes.structure) Ctypes.field

val pari_mt_get :
  ( (mt_state Ctypes.structure Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    pari_mt Ctypes.structure )
  Ctypes.field

val pari_mt_submit :
  ( (mt_state Ctypes.structure Ctypes_static.ptr ->
    Signed.long ->
    ('kind, 'structure) typ ->
    unit)
    Ctypes_static.static_funptr,
    pari_mt Ctypes.structure )
  Ctypes.field

val pari_mt_end :
  (unit Ctypes_static.static_funptr, pari_mt Ctypes.structure) Ctypes.field

val parfor_iter : parfor_iter Ctypes.structure Ctypes.typ

val parfor_iter_pending :
  (Signed.long, parfor_iter Ctypes.structure) Ctypes.field

val parfor_iter_worker :
  (('kind, 'structure) typ, parfor_iter Ctypes.structure) Ctypes.field

val parfor_iter_pt :
  (pari_mt Ctypes.structure, parfor_iter Ctypes.structure) Ctypes.field

val parfor_t : parfor_t Ctypes.structure Ctypes.typ

val parfor_t_a :
  (('kind, 'structure) typ, parfor_t Ctypes.structure) Ctypes.field

val parfor_t_b :
  (('kind, 'structure) typ, parfor_t Ctypes.structure) Ctypes.field

val parfor_t_iter :
  (parfor_iter Ctypes.structure, parfor_t Ctypes.structure) Ctypes.field

val parforeach_t : parforeach_t Ctypes.structure Ctypes.typ

val parforeach_t_x :
  (('kind, 'structure) typ, parforeach_t Ctypes.structure) Ctypes.field

val parforeach_t_W :
  (('kind, 'structure) typ, parforeach_t Ctypes.structure) Ctypes.field

val parforeach_t_i : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_l : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field

val parforeach_t_iter :
  (parfor_iter Ctypes.structure, parforeach_t Ctypes.structure) Ctypes.field

val parforprime_t : parforprime_t Ctypes.structure Ctypes.typ

val parforprime_t_v :
  (('kind, 'structure) typ, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_forprime :
  (forprime_t Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_iter :
  (parfor_iter Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforvec_t : parforvec_t Ctypes.structure Ctypes.typ

val parforvec_t_v :
  (('kind, 'structure) typ, parforvec_t Ctypes.structure) Ctypes.field

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

val nfmaxord_t_T :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dT :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_T0 :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_unscale :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dK :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_index :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_basis :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_r1 : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_basden :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dTP :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dTE :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dKP :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dKE :
  (('kind, 'structure) typ, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_certify : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field
val qfr_data : qfr_data Ctypes.structure Ctypes.typ

val qfr_data_D :
  (('kind, 'structure) typ, qfr_data Ctypes.structure) Ctypes.field

val qfr_data_sqrtD :
  (('kind, 'structure) typ, qfr_data Ctypes.structure) Ctypes.field

val qfr_data_isqrtD :
  (('kind, 'structure) typ, qfr_data Ctypes.structure) Ctypes.field

val fp_chk_fun : fp_chk_fun Ctypes.structure Ctypes.typ

val fp_chk_fun_f :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_init :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_post :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_data :
  (unit Ctypes_static.ptr, fp_chk_fun Ctypes.structure) Ctypes.field

val fp_chk_fun_skipfirst :
  (Signed.long, fp_chk_fun Ctypes.structure) Ctypes.field

val zlog_s : zlog_s Ctypes.structure Ctypes.typ
val zlog_s_bid : (('kind, 'structure) typ, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_P : (('kind, 'structure) typ, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_k : (('kind, 'structure) typ, zlog_s Ctypes.structure) Ctypes.field

val zlog_s_sprk :
  (('kind, 'structure) typ, zlog_s Ctypes.structure) Ctypes.field

val zlog_s_archp :
  (('kind, 'structure) typ, zlog_s Ctypes.structure) Ctypes.field

val zlog_s_mod : (('kind, 'structure) typ, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_U : (('kind, 'structure) typ, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_hU : (Signed.long, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_no2 : (int, zlog_s Ctypes.structure) Ctypes.field
val bb_group : bb_group Ctypes.structure Ctypes.typ

val bb_group_mul :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_pow :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_rand :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_hash :
  ( (('kind, 'structure) typ -> pari_ulong) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal :
  ( (('kind, 'structure) typ -> ('kind, 'structure) typ -> int)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal1 :
  ( (('kind, 'structure) typ -> int) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_easylog :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_field : bb_field Ctypes.structure Ctypes.typ

val bb_field_red :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_add :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_mul :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_neg :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_inv :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_equal0 :
  ( (('kind, 'structure) typ -> int) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_s :
  ( (unit Ctypes_static.ptr -> Signed.long -> ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_algebra : bb_algebra Ctypes.structure Ctypes.typ

val bb_algebra_red :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_add :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sub :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_mul :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sqr :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_one :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_zero :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_ring : bb_ring Ctypes.structure Ctypes.typ

val bb_ring_add :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_mul :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_sqr :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) typ ->
    ('kind, 'structure) typ)
    Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val buchimag :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val buchreal :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zidealstar :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zidealstarinit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zidealstarinitgen :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val factmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpbern : Signed.long -> Signed.long -> unit

val simplefactmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val listkill : ('kind, 'structure) typ -> unit

val isprincipalforce :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val isprincipalgen :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val isprincipalgenforce :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2ms_ker : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2ms_to_f2m :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2c_to_zc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2c_to_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2m_f2c_gauss :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2m_f2c_invimage :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2m_f2c_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2m_deplin : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2m_det : ('kind, 'structure) typ -> pari_ulong
val f2m_det_sp : ('kind, 'structure) typ -> pari_ulong

val f2m_gauss :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2m_inv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2m_invimage :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2m_ker : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2m_ker_sp :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2m_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2m_powu : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val f2m_rank : ('kind, 'structure) typ -> Signed.long
val f2m_row : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2m_rowslice :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val f2m_to_f2ms : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2m_to_flm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2m_to_zm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2m_to_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2m_transpose : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2v_add_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val f2v_and_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val f2v_dotproduct :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong

val f2v_equal0 : ('kind, 'structure) typ -> int
val f2v_hamming : ('kind, 'structure) typ -> pari_ulong

val f2v_negimply_inplace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val f2v_or_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val f2v_slice :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val f2v_subset : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val f2v_to_flv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val matid_f2m : Signed.long -> ('kind, 'structure) typ

val f2x_f2xq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2x_f2xqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2x_frobenius : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2x_1_add : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_deflate :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2x_degfact : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2x_degree : ('kind, 'structure) typ -> Signed.long
val f2x_deriv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val f2x_eval : ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val f2x_even_odd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val f2x_extgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val f2x_gcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_get_red : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_halfgcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_issquare : ('kind, 'structure) typ -> int
val f2x_matfrobenius : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_recip : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_rem :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_shift :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2x_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2x_sqrt : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_to_f2v :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2x_to_f2xx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2x_to_flx : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2x_to_zx : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2x_valrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val f2xc_to_flxc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2xc_to_zxc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xv_to_f2m :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2xv_to_flxv_inplace : ('kind, 'structure) typ -> unit
val f2xv_to_zxv_inplace : ('kind, 'structure) typ -> unit

val f2xx_f2x_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xx_f2x_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xx_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xx_deriv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xx_renormalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2xx_to_kronecker :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2xx_to_flxx : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2xx_to_zxx : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xx_to_f2xc :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val f2xxv_to_f2xm :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val f2xxc_to_zxxc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xy_f2xq_evalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xy_f2xqv_evalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xy_degreex : ('kind, 'structure) typ -> Signed.long

val f2xn_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val f2xn_inv : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val f2xn_red : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2xq_artin_schreier :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xq_autpow :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_conjvec :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xq_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xq_invsafe :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xq_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_matrix_pow :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_pow_init :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_pow_table :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xq_sqrt :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xq_sqrt_fast :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_sqrtn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val f2xq_trace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong

val f2xqx_f2xq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_f2xq_mul_to_monic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_f2xqxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_f2xqxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_disc :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val f2xqx_extgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val f2xqx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_get_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_halfgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_halfgcd_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val f2xqx_invbarrett :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_ispower :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val f2xqx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_normalize :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_rem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqx_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqxq_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqxq_invsafe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqxq_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqxq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqxq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqxq_autpow :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqxq_auttrace :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqxqv_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val flm_to_f2m : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flv_to_f2v : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flx_to_f2x : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flxc_to_f2xc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flxx_to_f2xx : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flxxc_to_f2xxc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val kronecker_to_f2xqx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rg_to_f2xq :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_to_f2m : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgv_to_f2v : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgx_to_f2x : ('kind, 'structure) typ -> ('kind, 'structure) typ
val z_to_f2x : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val zm_to_f2m : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_to_f2v : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zx_to_f2x : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxx_to_f2xx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val const_f2v : Signed.long -> ('kind, 'structure) typ

val gener_f2xq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val get_f2xq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  bb_field Ctypes.structure Ctypes_static.ptr

val monomial_f2x : Signed.long -> Signed.long -> ('kind, 'structure) typ
val pol1_f2xx : Signed.long -> Signed.long -> ('kind, 'structure) typ
val polx_f2xx : Signed.long -> Signed.long -> ('kind, 'structure) typ

val random_f2xqx :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2x_teichmuller :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2xq_ellcard :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_ellgens :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xq_ellgroup :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val f2xq_elltwist :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val f2xqe_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_changepoint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_changepointinv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_dbl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_neg :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_tatepairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqe_weilpairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val get_f2xqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_f2xqe :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val random_f2xqe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f3c_to_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f3c_to_zc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f3m_ker : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f3m_ker_sp :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f3m_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f3m_row : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val f3m_to_flm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f3m_to_zm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f3m_to_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f3m_transpose : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f3v_to_flv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f3v_coeff : ('kind, 'structure) typ -> Signed.long -> pari_ulong
val f3v_clear : ('kind, 'structure) typ -> Signed.long -> unit
val f3v_set : ('kind, 'structure) typ -> Signed.long -> pari_ulong -> unit
val flm_to_f3m : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flv_to_f3v : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgm_to_f3m : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgv_to_f3v : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zm_to_f3m : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_to_f3v : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zero_f3m_copy : Signed.long -> Signed.long -> ('kind, 'structure) typ
val zero_f3v : Signed.long -> ('kind, 'structure) typ
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
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

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

val fle_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fle_dbl :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val fle_changepoint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val fle_changepointinv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val fle_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fle_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fle_mulu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fle_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fle_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fle_tatepairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fle_to_flj : ('kind, 'structure) typ -> ('kind, 'structure) typ

val fle_weilpairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val flj_add_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flj_changepointinv_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flj_dbl_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flj_mulu_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flj_neg : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flj_to_fle :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flj_to_fle_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val fljv_factorback_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val random_fle :
  pari_ulong -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val random_fle_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val random_flj_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flc_to_zc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flc_to_zc_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ

val flm_flc_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_flc_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_adjoint :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_deplin :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_det : ('kind, 'structure) typ -> pari_ulong -> pari_ulong
val flm_det_sp : ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flm_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_intersect :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_intersect_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_inv : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_ker : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_ker_sp :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val flm_rank : ('kind, 'structure) typ -> pari_ulong -> Signed.long
val flm_to_zm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flm_to_zm_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flv_to_zv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val fl_to_flx : pari_ulong -> Signed.long -> ('kind, 'structure) typ
val fl2_equal1 : ('kind, 'structure) typ -> int

val fl2_inv_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fl2_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fl2_norm_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fl2_pow_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fl2_sqr_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fl2_sqrt_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val fl2_sqrtn_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flm_to_flxv :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flm_to_flxx :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val flv_flm_polint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val flv_inv : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val flv_inv_inplace : ('kind, 'structure) typ -> pari_ulong -> unit

val flv_inv_pre_inplace :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> unit

val flv_inv_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flv_invvandermonde :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flv_polint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val flv_prod : ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flv_prod_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> pari_ulong

val flv_roots_to_pol :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val flv_to_flx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flx_fl_add :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_fl_mul :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_fl_mul_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_fl_mul_to_monic :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_fl_sub :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_fl2_eval_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_flv_multieval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_flxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_flxq_eval_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_flxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_flxqv_eval_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_frobenius :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_frobenius_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_laplace :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_newton :
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_blocks :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val flx_composedprod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_composedsum :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_convol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_deflate :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flx_deriv : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val flx_diff1 : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_digits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_div_by_x_x :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flx_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flx_divrem_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flx_double :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_eval : ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> pari_ulong

val flx_eval_powers_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val flx_eval_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val flx_extgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flx_extgcd_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flx_extresultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong

val flx_extresultant_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong

val flx_fromnewton :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_gcd_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_get_red :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_get_red_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_halfgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_halfgcd_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flx_halfgcd_all_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flx_halfgcd_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_halve : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_inflate :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flx_integ : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_invbarrett :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_invlaplace :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_is_squarefree : ('kind, 'structure) typ -> pari_ulong -> int
val flx_is_smooth : ('kind, 'structure) typ -> Signed.long -> pari_ulong -> int

val flx_is_smooth_pre :
  ('kind, 'structure) typ -> Signed.long -> pari_ulong -> pari_ulong -> int

val flx_matfrobenius :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_matfrobenius_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_mod_xn1 :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_mod_xnm1 :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_neg : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_neg_inplace :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_normalize :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_powu :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_powu_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_recip : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flx_red : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_rem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_rem_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_renormalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flx_rescale :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_resultant :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flx_resultant_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val flx_shift :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flx_splitting :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flx_sqr : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_sqr_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_translate1 :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_translate1_basecase :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_to_flv :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flx_to_flxx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flx_to_zx : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flx_to_zx_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ

val flx_triple :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_val : ('kind, 'structure) typ -> Signed.long

val flx_valrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val flxc_flxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxc_flxqv_eval_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxc_flxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxc_flxq_eval_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxc_eval_powers_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxc_neg : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxc_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxc_to_zxc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val flxm_flx_add_shallow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxm_eval_powers_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxm_neg : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxm_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxm_to_flxxv :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flxm_to_zxm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flxt_red : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxv_flc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxv_flv_multieval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxv_flx_fromdigits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxv_composedsum :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxv_prod : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val flxv_red : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxv_to_flm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flxv_to_flxx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flxv_to_zxv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flxv_to_zxv_inplace : ('kind, 'structure) typ -> unit

val flxn_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxn_div_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxn_exp :
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxn_expint :
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxn_inv :
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxn_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxn_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxn_sqr :
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxn_sqr_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxn_red : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flxq_autpow :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_autpow_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_autpowers :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_autsum :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_auttrace :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_auttrace_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_charpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_conjvec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_div_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_inv_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_invsafe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_invsafe_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_issquare :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong -> int

val flxq_is2npower :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  int

val flxq_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_lroot :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val flxq_lroot_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_lroot_fast :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val flxq_lroot_fast_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_matrix_pow :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_matrix_pow_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_minpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_minpoly_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_norm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flxq_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_pow_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_pow_init :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_pow_init_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_pow_table_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_pow_table :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_powu_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_powers_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_sqr_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_sqrt :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_sqrt_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_sqrtn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxq_trace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flxqc_flxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_flxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqv_dotproduct :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqv_dotproduct_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val rg_to_f2 : ('kind, 'structure) typ -> pari_ulong
val rg_to_fl : ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val rg_to_flxq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val rgx_to_flx :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val rgxv_to_flxv :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val z_to_flx :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val zxv_to_flxv :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val zxt_to_flxt :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val gener_flxq :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val get_flxq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  bb_field Ctypes.structure Ctypes_static.ptr

val monomial_flx :
  pari_ulong -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val random_flx :
  Signed.long -> Signed.long -> pari_ulong -> ('kind, 'structure) typ

val zero_flxc : Signed.long -> Signed.long -> ('kind, 'structure) typ

val zero_flxm :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val zlx_translate1 :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val zx_to_flx : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxx_fl_mul :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flxx_flx_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxx_flx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxx_flx_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxx_laplace :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxx_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxx_blocks :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val flxx_deriv :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxx_double :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxx_invlaplace :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxx_neg : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxx_renormalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flxx_shift :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val flxx_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxx_swap :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val flxx_to_flm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flxx_to_flx : ('kind, 'structure) typ -> ('kind, 'structure) typ

val flxx_to_flxc :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val flxx_to_zxx : ('kind, 'structure) typ -> ('kind, 'structure) typ

val flxx_translate1 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val flxx_triple :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxxc_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxxc_to_zxxc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flxxm_to_zxxm : ('kind, 'structure) typ -> ('kind, 'structure) typ

val flxxv_to_flxm :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val flxxn_red :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val flxy_flx_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxy_flx_translate :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxy_flxqv_evalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxy_flxqv_evalx_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxy_flxq_evalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxy_flxq_evalx_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxy_evalx :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flxy_evalx_powers_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxy_evalx_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxyqq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqv_roots_to_pol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val flxqxc_flxqxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxc_flxqxqv_eval_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxc_flxqxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_autpow :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_autpow_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_autsum :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_autsum_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_auttrace :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_auttrace_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_div_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_inv_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_invsafe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_invsafe_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_matrix_pow :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_minpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_minpoly_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_pow_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_powers_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_powu_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_sqr_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxv_prod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_flxqxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_flxqxqv_eval_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_flxqxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_flxqxq_eval_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_flxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_flxq_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_flxq_mul_to_monic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_flxq_mul_to_monic_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_newton :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_newton_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_composedsum :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_disc :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_div_by_x_x :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxqx_div_by_x_x_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxqx_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxqx_divrem_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxqx_dotproduct :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_extgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxqx_extgcd_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxqx_fromnewton :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_fromnewton_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_gcd_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_get_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_get_red_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_halfgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_halfgcd_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxqx_halfgcd_all_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxqx_halfgcd_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_invbarrett :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_invbarrett_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_normalize :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_normalize_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_powu_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_red_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_rem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_rem_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_resultant_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_safegcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_saferesultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_sqr_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxn_expint :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxn_expint_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxn_inv :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxn_inv_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxn_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxn_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxn_sqr :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxn_sqr_pre :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxy_degreex : ('kind, 'structure) typ -> Signed.long

val flxy_eval_powers_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fly_to_flxy :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val kronecker_to_flxqx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val kronecker_to_flxqx_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val rgx_to_flxqx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val get_flxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val pol1_flxx : Signed.long -> Signed.long -> ('kind, 'structure) typ
val polx_flxx : Signed.long -> Signed.long -> ('kind, 'structure) typ

val random_flxqx :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val zlxx_translate1 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val zxx_to_kronecker :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val flxq_ellcard :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_ellgens :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_ellgroup :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxq_elltwist :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val flxq_ellj :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_ellj_to_a4a6 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val flxqe_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_changepoint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_changepointinv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_dbl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_neg :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_tatepairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_weilpairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqe_weilpairing_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val zxx_to_flxx :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val zxxt_to_flxxt :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val zxxv_to_flxxv :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val get_flxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_flxqe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val random_flxqe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val polisclass : ('kind, 'structure) typ -> Signed.long
val fl_elltrace : pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fl_elltrace_cm :
  Signed.long -> pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fp_ellcard :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_elldivpol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_ellgens :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_ellgroup :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fp_ellj :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_ellj_to_a4a6 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val fp_elljissupersingular :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val fp_elltwist :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val fp_ffellcard :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_changepoint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_changepointinv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_dbl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_neg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpe_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_to_fpj : ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpe_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpe_tatepairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpe_weilpairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpj_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpj_dbl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpj_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpj_neg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpj_to_fpe :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxq_ellcard :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_ellcard_supersingular :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_elldivpol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_ellgens :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_ellgroup :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpxq_ellj :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_elljissupersingular :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val fpxq_elltwist :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val fpxqe_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_changepoint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_changepointinv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_dbl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_neg :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_tatepairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqe_weilpairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_elljissupersingular :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val fq_ellcard_supersingular :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rge_to_fpe :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rge_to_fpxqe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val get_fpe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  bb_group Ctypes.structure Ctypes_static.ptr

val ellsupersingularj_fpxq :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val elltrace_extension :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val random_fpe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val random_fpxqe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_issquare : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val fp_fpx_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_fpxq_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpv_fpm_polint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fpv_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpv_invvandermonde :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpv_polint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fpv_roots_to_pol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fpx_fp_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fp_add_shallow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fp_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fp_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fp_mul_to_monic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fp_mulspec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fpx_fp_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fp_sub_shallow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fpv_multieval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fpxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fpxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_fpxv_multirem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_frobenius :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_laplace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_newton :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_center_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_chinese_coprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_composedprod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_composedsum :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_convol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_deriv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_digits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_disc :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_div_by_x_x :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpx_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpx_divu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_dotproduct :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_extgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpx_extresultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpx_fromnewton :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_gcd_check :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_get_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_halve :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_halfgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_halfgcd_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpx_integ :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_invbarrett :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_invlaplace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_is_squarefree :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val fpx_matfrobenius :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_mulspec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val fpx_mulu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_neg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_normalize :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_rem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_rescale :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_valrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val fpxc_fpxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxc_fpxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxm_fpxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_autpow :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_autpowers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_autsum :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_auttrace :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_charpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_conjvec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_invsafe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_issquare :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val fpxq_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_matrix_pow :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_minpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_norm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_sqrt :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_sqrtn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpxq_trace :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqc_to_mod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqm_autsum :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxt_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxv_fpx_fromdigits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxv_chinese :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpxv_composedsum :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxv_factorback :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fpxv_prod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxv_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxn_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxn_exp :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxn_expint :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxn_inv :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxn_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxn_sqr :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_issquare :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val fq_ispower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fq_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqc_to_mod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_to_mod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqv_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val z_to_fpx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gener_fpxq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gener_fpxq_local :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val get_fpxq_star :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpx :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2x_ddf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2x_factor : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2x_factor_squarefree : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2x_is_irred : ('kind, 'structure) typ -> int
val flx_ddf : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_ddf_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_is_irred : ('kind, 'structure) typ -> pari_ulong -> int
val flx_is_totally_split : ('kind, 'structure) typ -> pari_ulong -> int

val flx_ispower :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val flx_degfact :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_factor :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_factor_squarefree :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_factor_squarefree_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_nbfact : ('kind, 'structure) typ -> pari_ulong -> Signed.long

val flx_nbfact_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> Signed.long

val flx_nbfact_frobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long

val flx_nbfact_frobenius_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  Signed.long

val flx_nbfact_by_degree :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_nbroots : ('kind, 'structure) typ -> pari_ulong -> Signed.long
val flx_oneroot : ('kind, 'structure) typ -> pari_ulong -> pari_ulong
val flx_oneroot_split : ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flx_oneroot_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> pari_ulong

val flx_oneroot_split_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> pari_ulong

val flx_roots : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_roots_pre :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flx_rootsff :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val fpx_ddf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_ddf_degree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fpx_degfact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_factor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_factor_squarefree :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_is_irred : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val fpx_is_totally_split :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val fpx_ispower :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val fpx_nbfact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val fpx_nbfact_frobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fpx_nbroots :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val fpx_oneroot :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_oneroot_split :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_roots :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_roots_mult :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_rootsff :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_split_part :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_ddf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_degfact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_factor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_factor_squarefree :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqx_roots :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val flx_factorff_irred :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_ffintersect :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val flx_ffisom :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_ffisom_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_frobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_frobenius_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_ddf :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_ddf_degree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long

val flxqx_degfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_factor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_factor_squarefree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_factor_squarefree_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_ispower :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val flxqx_is_squarefree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long

val flxqx_nbfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long

val flxqx_nbfact_frobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long

val flxqx_nbfact_by_degree :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_nbroots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long

val flxqx_roots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqxq_halffrobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val fpx_factorff :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_factorff_irred :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_ffintersect :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val fpx_ffisom :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_ffisom_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_frobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_ddf :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_ddf_degree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fpxqx_degfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_factor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_factor_squarefree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_ispower :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val fpxqx_nbfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fpxqx_nbfact_frobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fpxqx_nbroots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fpxqx_roots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_split_part :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_halffrobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_is_squarefree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fqx_ispower :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val fqx_nbfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fqx_nbroots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val factorff :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val factormod0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val factormodddf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val factormodsqf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_parse_tp :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  int

val polrootsff :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rootmod0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fpxqx_fpxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_fpxqxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_fpxqxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_digits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_disc :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_div_by_x_x :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpxqx_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpxqx_dotproduct :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_extgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpxqx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_get_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_halfgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_halfgcd_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpxqx_invbarrett :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_rem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_to_mod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_invsafe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_matrix_pow :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_minpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_autpow :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_autsum :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxq_auttrace :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxt_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxv_fpxqx_fromdigits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxv_prod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxv_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxn_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxn_exp :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxn_expint :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxn_inv :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxn_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqxn_sqr :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxx_fp_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxx_fpx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxx_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxx_deriv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxx_halve :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxx_integ :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxx_mulu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxx_neg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxx_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxx_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxy_fpxq_evalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxy_fpxqv_evalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxy_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxy_evalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxy_evaly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fpxyqq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxc_to_mod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxm_to_mod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val kronecker_to_fpxqx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val get_fpxqx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpxqx :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val flc_flv_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flc_to_mod :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_fl_add :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flm_fl_mul :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flm_fl_mul_inplace :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> unit

val flm_fl_mul_pre :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_fl_sub :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flm_flc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_flc_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_flc_mul_pre_flx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val flm_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_center :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_mul_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_neg : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_powers :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flm_powu :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flm_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flm_to_mod :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_transpose : ('kind, 'structure) typ -> ('kind, 'structure) typ

val flv_fl_div :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flv_fl_div_inplace :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> unit

val flv_fl_mul :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flv_fl_mul_inplace :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> unit

val flv_fl_mul_part_inplace :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> Signed.long -> unit

val flv_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flv_add_inplace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong -> unit

val flv_center :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val flv_dotproduct :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flv_dotproduct_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val flv_neg : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val flv_neg_inplace : ('kind, 'structure) typ -> pari_ulong -> unit

val flv_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flv_sub_inplace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong -> unit

val flv_sum : ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flx_dotproduct :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flx_dotproduct_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fp_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpc_fpv_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpc_fp_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpc_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpc_center_inplace :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val fpc_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpc_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpm_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_fp_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_fpc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_fpc_mul_fpx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fpm_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_center_inplace :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val fpm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpm_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpms_fpc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpms_fpcs_solve :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpms_fpcs_solve_safe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpms_leftkernel_elt :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpc_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpc_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpv_fpms_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpv_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpv_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpv_dotproduct :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpv_dotsquare :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpv_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpv_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpvv_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpx_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxc_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxm_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zabm_ker :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zabm_indexrank :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zabm_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zabm_inv_ratlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zabm_pseudoinv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zv_zms_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zpms_zpcs_solve :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gen_fpm_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gen_zpm_dixon_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gen_matid :
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val matid_flm : Signed.long -> ('kind, 'structure) typ

val matid_f2xqm :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val matid_flxqm :
  Signed.long ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val random_flv : Signed.long -> pari_ulong -> ('kind, 'structure) typ

val random_fpc :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val random_fpv :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val scalar_flm : Signed.long -> Signed.long -> ('kind, 'structure) typ

val zcs_to_zc :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zms_to_zm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zms_zc_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zmv_to_flmv :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flx_teichmuller :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val z2_sqrt : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zp_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zp_exp :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val zp_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zp_invlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zp_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val zp_sqrt :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zp_sqrtlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zp_sqrtnlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpm_invlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpx_frobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpx_zpxq_liftroot :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpx_zpxq_liftroot_ea :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val zpx_liftfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpx_liftroot :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpx_liftroots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpx_roots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxq_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxq_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxq_invlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxq_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxq_sqrt :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxq_sqrtnlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxqm_prodfrobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxqx_digits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxqx_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zpxqx_liftfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxqx_liftroot :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxqx_liftroot_vald :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxqx_liftroots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxqx_roots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpxqx_zpxqxq_liftroot :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zq_sqrtnlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zqx_zqxq_liftroot :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zqx_liftfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zqx_liftroot :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zqx_roots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gen_zpm_dixon :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_zpm_newton :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_zpx_dixon :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_zpx_newton :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val polteichmuller :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val polhensellift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val quadratic_prec_mask : Signed.long -> pari_ulong
val qx_factor : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zx_factor : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zx_is_irred : ('kind, 'structure) typ -> Signed.long

val zx_squff :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val polcyclofactors : ('kind, 'structure) typ -> ('kind, 'structure) typ
val poliscyclo : ('kind, 'structure) typ -> Signed.long
val poliscycloprod : ('kind, 'structure) typ -> Signed.long

val rg_rgc_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_rg_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_rg_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_rg_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_rg_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_rgm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_rgv_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_is_ei : ('kind, 'structure) typ -> Signed.long
val rgc_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rg_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rg_add_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rg_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rg_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rg_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rg_sub_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rgc_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rgv_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_det_triangular : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgm_is_qm : ('kind, 'structure) typ -> int
val rgm_is_zm : ('kind, 'structure) typ -> int
val rgm_isdiagonal : ('kind, 'structure) typ -> int
val rgm_isidentity : ('kind, 'structure) typ -> int
val rgm_isscalar : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val rgm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_multosym :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_powers :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_sumcol : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_transmul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_transmultosym :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgmrow_zc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgm_zc_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_zm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgmrow_rgc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgv_rgm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_rgc_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_rg_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_dotproduct :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_dotsquare : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgv_is_zmv : ('kind, 'structure) typ -> int
val rgv_kill0 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgv_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgv_prod : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_sum : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_sumpart :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgv_sumpart2 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgv_zc_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_zm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rgm_eval :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rgmv_eval :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val isdiagonal : ('kind, 'structure) typ -> int

val scalarcol :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val scalarcol_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val scalarmat :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val scalarmat_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val scalarmat_s : Signed.long -> Signed.long -> ('kind, 'structure) typ

val kronecker_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qx_zxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxq_charpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qxq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxq_to_mod_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxqc_to_mod_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxqm_to_mod_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxqv_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxqx_homogenous_evalpow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxqx_to_mod_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxqxv_to_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxv_qxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxy_qxq_evalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rg_rgx_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rg_get_0 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rg_get_1 : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rg_to_rgc :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_to_rgxv :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_to_rgxv_reverse :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_to_rgxx :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgv_to_rgx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgv_to_rgm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgv_to_rgx_reverse :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_rgxq_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgx_rgxqv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgx_rgxn_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_rgxnv_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_rg_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rg_add_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rg_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rg_divexact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rg_eval_bk :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rg_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rg_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rgv_eval :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_addmulxn_shallow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_addmulxn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_addspec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_addspec_shallow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_affine :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgx_blocks :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_deflate :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_deriv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_digits :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_div_by_x_x :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgx_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgx_divs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_equal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val rgx_even_odd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val rgx_homogenize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_homogenous_evalpow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgx_inflate :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_mul_i :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_mul_normalized :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_mul2n :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_mulxn :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_mulhigh_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_muls : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_mulspec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgx_normalize : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_pseudodivrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgx_pseudorem :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_recip : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgx_recip_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgx_recip_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rem :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_renormalize_lg :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_rescale :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_rotate_shallow :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_shift :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_shift_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_splitting :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgx_sqr_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_sqrhigh_i :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_sqrspec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_to_rgc :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_translate :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_unscale :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxq_matrix_pow :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxq_norm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxq_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxq_trace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxqc_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxqm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxqm_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxqv_rgxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxqv_factorback :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxqv_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxqx_rgxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxqx_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgxqx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxqx_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxqx_pseudodivrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgxqx_pseudorem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxqx_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxqx_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxqx_translate :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxv_rgv_eval :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxv_prod : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxv_rescale :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxv_to_rgm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxv_unscale :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxx_to_rgm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxy_degreex : ('kind, 'structure) typ -> Signed.long
val rgxy_derivx : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxy_swap :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxy_swapspec :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxn_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxn_div_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxn_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxn_exp : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxn_expint :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxn_inv : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxn_inv_i :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxn_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxn_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxn_recip_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxn_red_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxn_reverse :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxn_sqr : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxn_sqrt :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxnv_red_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxn_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxn_powu_i :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val zx_translate :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_unscale2n :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zx_unscale :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_unscale_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_unscale_divpow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zx_z_unscale :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zxq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zxq_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zxqx_dvd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val brent_kung_optpow : Signed.long -> Signed.long -> Signed.long -> Signed.long

val gen_bkeval :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  int ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_bkeval_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val get_rg_algebra : unit -> bb_algebra Ctypes.structure Ctypes_static.ptr
val rfrac_deflate_order : ('kind, 'structure) typ -> Signed.long

val rfrac_deflate_max :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rfrac_deflate :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zgc_g_mul_inplace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val zgcs_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val g_zgc_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val g_zg_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zgc_g_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zgc_z_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zg_g_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zg_z_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zg_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zg_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zg_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zg_normalize : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zg_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val flc_lincomb1_inplace :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  unit

val vecsmall_prod : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_qc_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_det : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qm_ker : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgm_check_zm : ('kind, 'structure) typ -> string -> unit
val rgv_check_zv : ('kind, 'structure) typ -> string -> unit

val z_zc_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_zc_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_q_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_z_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_z_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_z_divexact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_z_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_zv_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_divexactu :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val zc_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_hnfremdiv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zc_is_ei : ('kind, 'structure) typ -> Signed.long

val zc_lincomb :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zc_lincomb1_inplace :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val zc_lincomb1_inplace_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val zc_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_reducemodlll :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_reducemodmatrix :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_z_mul : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zm_q_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_z_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_z_divexact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_z_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_det_triangular : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_diag_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_divexactu :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val zm_equal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val zm_equal0 : ('kind, 'structure) typ -> int

val zm_hnfdivrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zm_ishnf : ('kind, 'structure) typ -> int
val zm_isdiagonal : ('kind, 'structure) typ -> int
val zm_isidentity : ('kind, 'structure) typ -> int
val zm_isscalar : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val zm_max_lg : ('kind, 'structure) typ -> Signed.long

val zm_mul_diag :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_multosym :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_nm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_pow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_powu : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val zm_reducemodlll :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_reducemodmatrix :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_supnorm : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_transmul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_transmultosym :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_togglesign : ('kind, 'structure) typ -> unit

val zm_zm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zmrow_zc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zmrow_equal0 : ('kind, 'structure) typ -> Signed.long -> int
val zv_abscmp : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val zv_cmp : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val zv_dotsquare : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_max_lg : ('kind, 'structure) typ -> Signed.long
val zv_to_nv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_togglesign : ('kind, 'structure) typ -> unit
val gram_matrix : ('kind, 'structure) typ -> ('kind, 'structure) typ

val nm_z_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_to_flm : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val zm_to_zm : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_zc_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zmv_to_zmv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_abs : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_content : ('kind, 'structure) typ -> Signed.long

val zv_dotproduct :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val zv_equal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val zv_equal0 : ('kind, 'structure) typ -> int
val zv_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_neg_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_prod : ('kind, 'structure) typ -> Signed.long
val zv_prod_z : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_sum : ('kind, 'structure) typ -> Signed.long
val zv_sumpart : ('kind, 'structure) typ -> Signed.long -> Signed.long
val zv_to_flv : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val zv_z_mul : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zv_zm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zvv_equal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val kronecker_to_zxqx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val kronecker_to_zxx :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val qx_zx_rem :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qx_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qx_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxqm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxqm_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxqx_qxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxqx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxqx_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxqx_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_check_qx : ('kind, 'structure) typ -> string -> unit
val rgx_check_zx : ('kind, 'structure) typ -> string -> unit
val rgx_check_zxx : ('kind, 'structure) typ -> string -> unit

val z_zx_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_z_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_z_add_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_z_eval :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_z_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_z_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_affine :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zx_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zx_deriv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_digits :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxv_zx_fromdigits :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_div_by_x_1 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zx_divuexact :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val zx_equal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val zx_eval1 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zx_max_lg : ('kind, 'structure) typ -> Signed.long

val zx_mod_xnm1 :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val zx_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_mulspec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val zx_mulu : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val zx_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_rem :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_remi2n :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zx_rescale2n :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zx_rescale :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_rescale_lt : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_shifti :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zx_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_sqrspec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zx_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_val : ('kind, 'structure) typ -> Signed.long

val zx_valrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val zxc_to_flxc :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val zxm_to_flxm :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val zxqm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zxqm_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxqx_zxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zxqx_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxqx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zxt_remi2n :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zxv_z_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxv_dotproduct :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxv_equal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val zxv_remi2n :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zxx_z_divexact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxx_z_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxx_z_add_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxx_evalx0 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zxx_max_lg : ('kind, 'structure) typ -> Signed.long

val zxx_mul_kronecker :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zxx_renormalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zxx_sqr_kronecker :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxx_to_kronecker :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgxx_to_kronecker_spec :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val zxn_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zxn_sqr : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val scalar_zx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val scalar_zx_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zx_to_zx : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_z_divexact :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val alg_centralproj :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val alg_complete :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val alg_csa_table :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val alg_cyclic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val alg_get_absdim : ('kind, 'structure) typ -> Signed.long
val alg_get_abssplitting : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_aut : ('kind, 'structure) typ -> ('kind, 'structure) typ
val algaut : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_auts : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_b : ('kind, 'structure) typ -> ('kind, 'structure) typ
val algb : ('kind, 'structure) typ -> ('kind, 'structure) typ
val algcenter : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_center : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_char : ('kind, 'structure) typ -> ('kind, 'structure) typ
val algchar : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_degree : ('kind, 'structure) typ -> Signed.long
val algdegree : ('kind, 'structure) typ -> Signed.long
val alg_get_dim : ('kind, 'structure) typ -> Signed.long
val algdim : ('kind, 'structure) typ -> Signed.long -> Signed.long
val alg_get_hasse_f : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alghassef : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_hasse_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alghassei : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_invbasis : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alginvbasis : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_multable : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_basis : ('kind, 'structure) typ -> ('kind, 'structure) typ
val algbasis : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_relmultable : ('kind, 'structure) typ -> ('kind, 'structure) typ
val algrelmultable : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_splitpol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_splittingfield : ('kind, 'structure) typ -> ('kind, 'structure) typ
val algsplittingfield : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_splittingbasis : ('kind, 'structure) typ -> ('kind, 'structure) typ

val alg_get_splittingbasisinv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ

val alg_get_splittingdata : ('kind, 'structure) typ -> ('kind, 'structure) typ
val algsplittingdata : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alg_get_tracebasis : ('kind, 'structure) typ -> ('kind, 'structure) typ

val alg_hasse :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val alg_hilbert :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val alg_matrix :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val alg_model :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val alg_quotient :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val algradical : ('kind, 'structure) typ -> ('kind, 'structure) typ

val algsimpledec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val algsimpledec_ss :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val algsubalg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val alg_type : ('kind, 'structure) typ -> Signed.long

val algadd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val algalgtobasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val algbasistoalg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val algcharpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val algdisc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val algdivl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val algdivr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val alggroup :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val alggroupcenter :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val alghasse :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val alginit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val algindex : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val alginv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val algisassociative : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val algiscommutative : ('kind, 'structure) typ -> int
val algisdivision : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val algisramified : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val algissemisimple : ('kind, 'structure) typ -> int
val algissimple : ('kind, 'structure) typ -> Signed.long -> int
val algissplit : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val algisdivl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val algisinv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val algmakeintegral :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val algmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val algmultable : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alglat_get_primbasis : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alglat_get_scalar : ('kind, 'structure) typ -> ('kind, 'structure) typ

val alglatadd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val alglatcontains :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val alglatelement :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val alglathnf :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val alglatindex :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val alglatinter :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val alglatmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val alglatlefttransporter :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val alglatrighttransporter :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val alglatsubset :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val algneg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val algnorm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val algpoleval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val algpow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val algprimesubalg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val algramifiedplaces : ('kind, 'structure) typ -> ('kind, 'structure) typ

val algrandom :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val algsplit : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val algtomatrix :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val algsqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val algsub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val algtableinit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val algtensor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val algtrace :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val algtype : ('kind, 'structure) typ -> Signed.long

val bnfgwgeneric :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val checkalg : ('kind, 'structure) typ -> unit

val checkhasse :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val checklat : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val conjclasses_algcenter :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoischardet :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val galoischarpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val galoischartable : ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfgrunwaldwang :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfgwkummer :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val f2ms_colelim :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2m_image : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2m_indexrank : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2m_suppl : ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqm_f2xqc_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqm_f2xqc_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqm_f2xqc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqm_deplin :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqm_det :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqm_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqm_ker :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqm_image :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqm_indexrank :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqm_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val f2xqm_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val f2xqm_rank :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val f2xqm_suppl :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val flm_image : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_indexrank :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_suppl : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flxqm_flxqc_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_flxqc_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_flxqc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_deplin :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_det :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_ker :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_image :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_indexrank :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqm_rank :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long

val flxqm_suppl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val fpm_fpc_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_fpc_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_deplin :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpm_det :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpm_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_image :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpm_indexrank :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpm_intersect :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_intersect_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpm_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_ker :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpm_rank : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val fpm_suppl :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fqm_fqc_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_fqc_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_fqc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_deplin :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_det :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_ker :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_image :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_indexrank :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqm_rank :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val fqm_suppl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qm_image_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qm_image : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_gauss :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_gauss_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qm_indexrank : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qm_inv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qm_rank : ('kind, 'structure) typ -> Signed.long

val rgm_fp_init :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgm_hadamard : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rgc_invimage :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_diagonal : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgm_diagonal_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgm_inv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rgm_inv_upper : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_invimage :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_solve :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_solve_realimag :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgms_structelim :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val zm_det : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zm_detmult : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_gauss :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_ker : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zm_imagecompl : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zm_indeximage : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zm_indexrank : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zm_inv_ratlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zm_pseudoinv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zm_rank : ('kind, 'structure) typ -> Signed.long

val zlm_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val closemodinvertible :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val deplin : ('kind, 'structure) typ -> ('kind, 'structure) typ
val det : ('kind, 'structure) typ -> ('kind, 'structure) typ
val det0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val det2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val detint : ('kind, 'structure) typ -> ('kind, 'structure) typ
val eigen : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gauss :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gaussmodulo :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gaussmodulo2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gen_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_gauss_pivot :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_det :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_ker :
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_matcolinvimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_matcolmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_matinvimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_matmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val image : ('kind, 'structure) typ -> ('kind, 'structure) typ
val image2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val imagecompl : ('kind, 'structure) typ -> ('kind, 'structure) typ
val indexrank : ('kind, 'structure) typ -> ('kind, 'structure) typ

val inverseimage :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ker : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mateigen :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val matimage0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val matker0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val rank : ('kind, 'structure) typ -> Signed.long

val reducemodinvertible :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val reducemodlll :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val split_realimag :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val suppl : ('kind, 'structure) typ -> ('kind, 'structure) typ

val flm_charpoly :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val flm_hess : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val fpm_charpoly :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpm_hess :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val frobeniusform :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val qm_minors_coprime :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_imz : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_imz_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val qm_imz_hnf : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_imz_hnfall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ

val qm_imq : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_imq_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val qm_imq_hnf : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_imq_hnfall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ

val qm_charpoly_zx : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qm_charpoly_zx_bound :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zm_charpoly : ('kind, 'structure) typ -> ('kind, 'structure) typ
val adj : ('kind, 'structure) typ -> ('kind, 'structure) typ
val adjsafe : ('kind, 'structure) typ -> ('kind, 'structure) typ
val caract : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val caradj :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val carberkowitz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val carhess : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val charpoly : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val charpoly0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val gnorm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gnorml1 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gnorml1_fake : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gnormlp :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gnorml2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gsupnorm : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gsupnorm_aux :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  unit

val gtrace : ('kind, 'structure) typ -> ('kind, 'structure) typ
val hess : ('kind, 'structure) typ -> ('kind, 'structure) typ

val intersect :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val jacobi : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val matadjoint0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val matcompanion : ('kind, 'structure) typ -> ('kind, 'structure) typ

val matrixqz0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val minpoly : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val qfgaussred : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qfgaussred_positive : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qfsign : ('kind, 'structure) typ -> ('kind, 'structure) typ

val apply0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val diagonal : ('kind, 'structure) typ -> ('kind, 'structure) typ
val diagonal_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ

val extract0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fold0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val genapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val genfold :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val genindexselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val genselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gtomat : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gtrans : ('kind, 'structure) typ -> ('kind, 'structure) typ

val matmuldiagonal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val matmultodiagonal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val matslice0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val parapply :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val parfor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parfor_init :
  parfor_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val parfor_next :
  parfor_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val parfor_stop : parfor_t Ctypes.structure Ctypes_static.ptr -> unit

val parforeach :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforeach_init :
  parforeach_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val parforeach_next :
  parforeach_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val parforeach_stop : parforeach_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprime_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val parforprime_next :
  parforprime_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val parforprime_stop : parforprime_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprimestep :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprimestep_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val parforvec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforvec_init :
  parforvec_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val parforvec_next :
  parforvec_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val parforvec_stop : parforvec_t Ctypes.structure Ctypes_static.ptr -> unit

val parselect :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val select0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val shallowextract :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val shallowmatextract :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val shallowtrans : ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val veccatapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val veccatselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val vecrange :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecrangess : Signed.long -> Signed.long -> ('kind, 'structure) typ

val vecselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val vecselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val vecslice0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val vecsum : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_diagonal : ('kind, 'structure) typ -> ('kind, 'structure) typ
val addhelp : string -> string -> unit
val arity0 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val alias0 : string -> string -> unit
val compile_str : string -> ('kind, 'structure) typ
val delete_var : unit -> Signed.long
val fetch_user_var : string -> Signed.long
val fetch_var : unit -> Signed.long
val fetch_var_higher : unit -> Signed.long

val fetch_var_value :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gp_embedded : string -> string
val gp_embedded_init : Signed.long -> Signed.long -> unit
val gp_read_str : string -> ('kind, 'structure) typ
val gp_read_str_bitprec : string -> Signed.long -> ('kind, 'structure) typ
val gp_read_str_prec : string -> Signed.long -> ('kind, 'structure) typ

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
val readseq : string -> ('kind, 'structure) typ

val safeel :
  ('kind, 'structure) typ -> Signed.long -> Signed.long Ctypes_static.ptr

val safelistel :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr

val strtoi : string -> ('kind, 'structure) typ
val strtor : string -> Signed.long -> ('kind, 'structure) typ
val varhigher : string -> Signed.long -> ('kind, 'structure) typ
val varlower : string -> Signed.long -> ('kind, 'structure) typ

val divisorslenstra :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val isprimeaprcl : ('kind, 'structure) typ -> Signed.long

val qfb0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val check_quaddisc :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  string ->
  unit

val check_quaddisc_imag :
  ('kind, 'structure) typ -> Signed.long Ctypes_static.ptr -> string -> unit

val check_quaddisc_real :
  ('kind, 'structure) typ -> Signed.long Ctypes_static.ptr -> string -> unit

val cornacchia :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val cornacchia2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val cornacchia2_sqrt :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val nucomp :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nudupl :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nupow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val primeform :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val primeform_u :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val qfb_1 : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbcomp :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbcomp_i :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbcompraw :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbcompraw_i :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbcornacchia :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbpow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbpow_i :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbpowraw :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val qfbpows : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val qfbred : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qfbred_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbred0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfbredsl2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbsolve :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qfbsqr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qfbsqr_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfisolvep :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfr3_comp :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val qfr3_compraw :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfr3_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val qfr3_red :
  ('kind, 'structure) typ ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val qfr3_rho :
  ('kind, 'structure) typ ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val qfr3_to_qfr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfr5_comp :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val qfr5_compraw :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfr5_dist :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qfr5_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val qfr5_red :
  ('kind, 'structure) typ ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val qfr5_rho :
  ('kind, 'structure) typ ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val qfr5_to_qfr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfr_data_init :
  ('kind, 'structure) typ ->
  Signed.long ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  unit

val qfr_to_qfr5 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val qfrsolvep :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val quadgen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val quadgen0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val quadpoly : ('kind, 'structure) typ -> ('kind, 'structure) typ
val quadpoly_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val quadpoly0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

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

val fl_powers :
  pari_ulong -> Signed.long -> pari_ulong -> ('kind, 'structure) typ

val fl_powers_pre :
  pari_ulong ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

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

val fp_2gener : ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_2gener_i :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_factored_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_ispower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val fp_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_pow_init :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_pow_table :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_pows :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_sqrt :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_sqrt_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_sqrtn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fpv_prod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_zv_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_zv_mod_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val z_chinese :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val z_chinese_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val z_chinese_coprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val z_chinese_post :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val z_chinese_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val z_factor_listp :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_nv_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_nv_mod_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zv_allpnqn : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_chinese :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zv_chinese_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zv_chinesetree :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_nv_mod_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zv_producttree : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_nv_mod_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zxc_nv_mod_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zxm_nv_mod_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zxx_nv_mod_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zideallog :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bestappr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bestapprpade :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val chinese :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val chinese1 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val chinese1_coprime_z : ('kind, 'structure) typ -> ('kind, 'structure) typ

val contfrac0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val contfracpnqn :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val fibo : Signed.long -> ('kind, 'structure) typ
val gboundcf : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gcf : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gcf2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val get_fp_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  bb_field Ctypes.structure Ctypes_static.ptr

val hilbert :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val hilbertii :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val istotient :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val krois : ('kind, 'structure) typ -> Signed.long -> Signed.long
val kroiu : ('kind, 'structure) typ -> pari_ulong -> Signed.long

val kronecker :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val krosi : Signed.long -> ('kind, 'structure) typ -> Signed.long
val kross : Signed.long -> Signed.long -> Signed.long
val kroui : pari_ulong -> ('kind, 'structure) typ -> Signed.long
val krouu : pari_ulong -> pari_ulong -> Signed.long

val lcmii :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_invgen :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val logint0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val logintall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val mpfact : Signed.long -> ('kind, 'structure) typ
val factorial_fl : Signed.long -> pari_ulong -> pari_ulong

val factorial_fp :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val muls_interval : Signed.long -> Signed.long -> ('kind, 'structure) typ
val mulu_interval : pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val mulu_interval_step :
  pari_ulong -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val ncv_chinese_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ncv_chinese_center_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nmv_chinese_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val nmv_chinese_center_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nonsquare_fl : pari_ulong -> pari_ulong

val nxcv_chinese_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val nxcv_chinese_center_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nxmv_chinese_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val nxv_chinese_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val nxv_chinese_center_tree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zv_chinese_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val odd_prime_divisors : ('kind, 'structure) typ -> ('kind, 'structure) typ
val pgener_fl : pari_ulong -> pari_ulong
val pgener_fl_local : pari_ulong -> ('kind, 'structure) typ -> pari_ulong
val pgener_fp : ('kind, 'structure) typ -> ('kind, 'structure) typ

val pgener_fp_local :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val pgener_zl : pari_ulong -> pari_ulong
val pgener_zp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val pnqn : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ramanujantau :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rootsof1_fl : pari_ulong -> pari_ulong -> pari_ulong

val rootsof1_fp :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rootsof1u_fp :
  pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val u_chinese_coprime :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val znlog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val znorder :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znprimroot : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val rgv_is_zvpos : ('kind, 'structure) typ -> int
val rgv_is_zvnon0 : ('kind, 'structure) typ -> int
val rgv_is_prv : ('kind, 'structure) typ -> int
val z_issquarefree_fact : ('kind, 'structure) typ -> Signed.long

val z_lsmoothen :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val z_smoothen :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val bigomega : ('kind, 'structure) typ -> Signed.long
val bigomegau : pari_ulong -> Signed.long
val boundfact : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val check_arith_pos :
  ('kind, 'structure) typ -> string -> ('kind, 'structure) typ

val check_arith_non0 :
  ('kind, 'structure) typ -> string -> ('kind, 'structure) typ

val check_arith_all :
  ('kind, 'structure) typ -> string -> ('kind, 'structure) typ

val clean_z_factor : ('kind, 'structure) typ -> ('kind, 'structure) typ
val core : ('kind, 'structure) typ -> ('kind, 'structure) typ

val coredisc2_fact :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val coredisc2u_fact :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong

val corepartial :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val core0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val core2 : ('kind, 'structure) typ -> ('kind, 'structure) typ

val core2partial :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val coredisc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val coredisc0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val coredisc2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val corediscs : Signed.long -> pari_ulong Ctypes_static.ptr -> Signed.long
val divisors : ('kind, 'structure) typ -> ('kind, 'structure) typ
val divisors_factored : ('kind, 'structure) typ -> ('kind, 'structure) typ

val divisors0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val divisorsu : pari_ulong -> ('kind, 'structure) typ
val divisorsu_moebius : ('kind, 'structure) typ -> ('kind, 'structure) typ
val divisorsu_fact : ('kind, 'structure) typ -> ('kind, 'structure) typ
val divisorsu_fact_factored : ('kind, 'structure) typ -> ('kind, 'structure) typ
val eulerphi : ('kind, 'structure) typ -> ('kind, 'structure) typ
val eulerphiu : pari_ulong -> pari_ulong
val eulerphiu_fact : ('kind, 'structure) typ -> pari_ulong

val factor_pn_1 :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val factor_pn_1_limit :
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val factoru_pow : pari_ulong -> ('kind, 'structure) typ

val fuse_z_factor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val is_z_factor : ('kind, 'structure) typ -> int
val is_z_factornon0 : ('kind, 'structure) typ -> int
val is_z_factorpos : ('kind, 'structure) typ -> int
val is_nf_factor : ('kind, 'structure) typ -> int
val is_nf_extfactor : ('kind, 'structure) typ -> int
val issquarefree : ('kind, 'structure) typ -> Signed.long
val numdiv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val numdivu : Signed.long -> Signed.long
val numdivu_fact : ('kind, 'structure) typ -> Signed.long
val omega : ('kind, 'structure) typ -> Signed.long
val omegau : pari_ulong -> Signed.long
val sumdiv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val sumdivk : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val uissquarefree : pari_ulong -> Signed.long
val uissquarefree_fact : ('kind, 'structure) typ -> Signed.long
val usumdiv_fact : ('kind, 'structure) typ -> ('kind, 'structure) typ

val usumdivk_fact :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val fpx_fpc_nfpoleval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val embed_t2 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val embednorm_t2 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val embed_norm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val check_zkmodule_i : ('kind, 'structure) typ -> int
val check_zkmodule : ('kind, 'structure) typ -> string -> unit
val checkbid : ('kind, 'structure) typ -> unit
val checkbid_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val checkbnf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val checkbnf_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val checkbnr : ('kind, 'structure) typ -> unit
val checkbnr_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val checkabgrp : ('kind, 'structure) typ -> unit
val checksqmat : ('kind, 'structure) typ -> Signed.long -> unit
val checknf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val checknf_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val checknfelt_mod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  string ->
  ('kind, 'structure) typ

val checkprid : ('kind, 'structure) typ -> unit
val checkprid_i : ('kind, 'structure) typ -> int
val checkrnf : ('kind, 'structure) typ -> unit
val checkrnf_i : ('kind, 'structure) typ -> int

val factoredpolred :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val factoredpolred2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoisapply :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val get_bnf :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val get_bnfpol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val get_nf :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val get_nfpol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val get_prid : ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealfrobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealfrobenius_aut :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealramfrobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealramfrobenius_aut :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealramgroups :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealramgroups_aut :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nf_get_allroots : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_prec : ('kind, 'structure) typ -> Signed.long

val nfmaxord_to_nf :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfcertify : ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfgaloismatrix :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfgaloismatrixapply :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfgaloispermtobasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfinit_basic :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  unit

val nfinit_complete :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val nfinit0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val nfinitred :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val nfinitred2 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val nfisincl0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfisisom :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfnewprec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val nfnewprec_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val nfpoleval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfsplitting0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nftyp : ('kind, 'structure) typ -> Signed.long
val polredord : ('kind, 'structure) typ -> ('kind, 'structure) typ
val polred : ('kind, 'structure) typ -> ('kind, 'structure) typ

val polred0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val polred2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val polredabs : ('kind, 'structure) typ -> ('kind, 'structure) typ

val polredabs0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val polredabs2 : ('kind, 'structure) typ -> ('kind, 'structure) typ

val polredabsall :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val poltomonic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rnfpolredabs :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfpolredbest :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val smallpolred : ('kind, 'structure) typ -> ('kind, 'structure) typ
val smallpolred2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val tschirnhaus : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_q_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_q_normalize :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zx_z_normalize :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zx_to_monic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zx_primitive_to_monic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zxx_q_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fq_to_nf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fqm_to_nfm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fqv_to_nfv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fqx_to_nfx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rg_nffix :
  string ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int ->
  ('kind, 'structure) typ

val rgv_nffix :
  string ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int ->
  ('kind, 'structure) typ

val rgx_nffix :
  string ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int ->
  ('kind, 'structure) typ

val zx_composedsum :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_compositum :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zpx_disc_val :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val zpx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zpx_monic_factor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpx_primedec :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zpx_reduced_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zpx_reduced_resultant_fast :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zpx_resultant_val :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long

val checkmodpr : ('kind, 'structure) typ -> unit

val compositum :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val compositum2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfdisc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val get_modpr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val indexpartial :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val modpr_genfq : ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_to_fq_init :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val nf_to_fq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfm_to_fqm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfv_to_fqv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfx_to_fqx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfx_to_monic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val nfbasis :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val nfcompositum :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfdiscfactors : ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfmaxord :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val nfmodpr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfmodprinit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfmodprinit0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfmodprlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfreducemodpr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val polcompositum0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealprimedec :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealprimedec_galois :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealprimedec_degrees :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealprimedec_kummer :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealprimedec_limit_f :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealprimedec_limit_norm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val poldiscfactors :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rnfbasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfdedekind :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfdet :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfdisc_factored :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rnfdiscf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfequation :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfequation0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfequation2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_pv_to_prv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_rnfeq :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_rnfeqsimple :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfequationall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rnfhnfbasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfisfree :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val rnflllgram :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfpolred :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfpseudobasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfsimplifybasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfsteinitz :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val factorial_lval : pari_ulong -> pari_ulong -> Signed.long

val zk_to_fq_init :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zk_to_fq :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxqv_to_fpm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zkmodprinit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealstar :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealstarprk :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgc_to_nfc :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_rgx_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_to_nfm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_to_nfx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_nfval : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val zc_nfvalrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val zc_prdvd : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val zm_zx_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_snf_gcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val algtobasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val basistoalg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ei_multable :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val get_nf_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  bb_field Ctypes.structure Ctypes_static.ptr

val famat_nfvalrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gpnfvalrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val idealfactorback :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int ->
  ('kind, 'structure) typ

val ideallist :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ideallist0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val gideallist :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ideallistarch :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ideallog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ideallogmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ideallog_units :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ideallog_units0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealprincipalunits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealstar0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealstarmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val indices_to_vec01 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val matalgtobasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val matbasistoalg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val multable :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_to_scalar_or_alg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_to_scalar_or_basis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_cxlog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfv_cxlog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfchecksigns :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val nfdiv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfdiveuc :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfembed :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfeltembed :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfeltembed_i :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfeltsign :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nffactorback :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfinv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfinvmodideal :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfissquare :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val nfispower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val nflogembed :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ

val nfm_det :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfm_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfm_ker :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfm_nfc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfmuli :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfnorm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfpolsturm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfpow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfpow_u :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val nfpowmodideal :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfsign :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfsign_arch :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfsign_from_logarch :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfsqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfsqri :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfsub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nftrace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val nfvalrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val polmod_nffix :
  string ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int ->
  ('kind, 'structure) typ

val polmod_nffix2 :
  string ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int ->
  ('kind, 'structure) typ

val pr_basis_perm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val pr_equal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val rnfalgtobasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfbasistoalg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfeltnorm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfelttrace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val set_sign_mod_divisor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val tablemul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val tablemul_ei :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val tablemul_ei_ej :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val tablemulvec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val tablesqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vec01_to_indices : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecsmall01_to_indices : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zk_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zk_multable :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zk_scalar_or_multable :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zkchinese :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zkchinese1 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zkchineseinit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zkmultable_capz : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zkmultable_inv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val fl_invgen :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val rm_round_maxrank : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_famat_limit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_cba : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_cba_extend :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_cba :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_ppgle :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_ppio :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_ppo :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famatv_factorback :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famatv_zv_factorback :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_z_gcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_div_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_idealfactor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_inv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val famat_inv_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_makecoprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val famat_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_mul_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_mulpow_shallow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val famat_mulpows_shallow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val famat_pow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_pow_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_pows_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val famat_reduce : ('kind, 'structure) typ -> ('kind, 'structure) typ
val famat_remove_trivial : ('kind, 'structure) typ -> ('kind, 'structure) typ
val famat_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_to_nf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val famat_to_nf_moddivisor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val famat_to_nf_modideal_coprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val famatsmall_reduce : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gpidealfactor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gpidealval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealhnf_z_factor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val idealhnf_z_factor_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val idealhnf_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealhnf_inv_z :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealhnf_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealadd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealaddmultoone :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealaddtoone :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealaddtoone0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealaddtoone_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealaddtoone_raw :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealappr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealappr0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealapprfact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealchinese :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealcoprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealcoprimefact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealdiv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealdiv0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealdivexact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealdivpowprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealdown :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealfactor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealfactor_limit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val idealfactor_partial :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealhnf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealhnf0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealhnf_principal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealhnf_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealhnf_two :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealintersect :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealinv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealismaximal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealispower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val idealmin :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealmul0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealmulpowprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealmulred :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealnumden :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealpow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealpow0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealpowred :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealpows :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val idealprod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealprodprime :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealprodval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val idealpseudomin :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealpseudomin_nonscalar :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealpseudominvec :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealpseudored :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealred0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealred_elt :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealredmodpower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val idealsqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealtwoelt :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealtwoelt0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealtwoelt2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val idealtyp :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val idealval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val isideal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long
val matreduce : ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfc_multable_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfc_nf_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nf_get_gtwist :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_get_gtwist1 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val nf_to_fp_coprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfdetint :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfdivmodpr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfhnf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfhnf0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfhnfmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfkermodpr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfmulmodpr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfpowmodpr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfreduce :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfsnf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfsnf0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfsolvemodpr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val prv_lcm_capz : ('kind, 'structure) typ -> ('kind, 'structure) typ
val prv_primes : ('kind, 'structure) typ -> ('kind, 'structure) typ

val pr_hnf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val pr_inv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val pr_inv_p : ('kind, 'structure) typ -> ('kind, 'structure) typ

val pr_uniformizer :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val sunits_makecoprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val to_famat :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val to_famat_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val u_ppo : pari_ulong -> pari_ulong -> pari_ulong

val vecdiv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecinv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecmul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecpow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecsqr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zkc_multable_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val eltreltoabs :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val eltabstorel :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val eltabstorel_lift :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_nfzk :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnf_build_nfabs :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rnf_zkabs : ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfeltup :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rnfcomplete : ('kind, 'structure) typ -> unit

val rnfeltabstorel :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfeltdown :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfeltdown0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfeltreltoabs :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfeltup :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfeltup0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfidealabstorel :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealdown :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealfactor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealhnf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rnfidealnormabs :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealnormrel :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealprimedec :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealreltoabs :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealreltoabs0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfidealtwoelement :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealup :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfidealup0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfinit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfinit0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val get_arith_zzm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val get_arith_z : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gen_ph_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_shanks_init :
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_shanks :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_shanks_sqrtn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_gener :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_ellgens :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_ellgroup :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_factored_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_order :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_select_order :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_plog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gen_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_pow_fold :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_pow_fold_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_pow_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_pow_init :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_pow_table :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  int ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_powu_fold :
  ('kind, 'structure) typ ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_powu_fold_i :
  ('kind, 'structure) typ ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_powu_i :
  ('kind, 'structure) typ ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_product :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val matdetmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val matimagemod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val matinvmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val matkermod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val matsolvemod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bernfrac : Signed.long -> ('kind, 'structure) typ
val bernpol : Signed.long -> Signed.long -> ('kind, 'structure) typ
val bernreal : Signed.long -> Signed.long -> ('kind, 'structure) typ
val bernvec : Signed.long -> ('kind, 'structure) typ
val constbern : Signed.long -> unit
val eulerfrac : Signed.long -> ('kind, 'structure) typ
val eulerpol : Signed.long -> Signed.long -> ('kind, 'structure) typ
val eulerreal : Signed.long -> Signed.long -> ('kind, 'structure) typ
val eulervec : Signed.long -> ('kind, 'structure) typ
val harmonic : pari_ulong -> ('kind, 'structure) typ
val harmonic0 : pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qr_init :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  int

val r_from_qr :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_babai :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_qr_init :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  int

val rgm_gram_schmidt :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val algdep : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val algdep0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val bestapprnf :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val forqfvec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  float ->
  Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val forqfvec0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val forqfvec1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val gaussred_from_qr :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lindep : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lindep_xadic : ('kind, 'structure) typ -> ('kind, 'structure) typ

val lindep_bit :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lindep_padic : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lindep0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val lindep2 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lindepfull_bit :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mathouseholder :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val matqr :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val minim :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val minim_raw :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val minim_zm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val minim2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfminim0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val qfperfection : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfrep0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val seralgdep :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val serdiffdep :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val vandermondeinverse :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val vandermondeinverseinit : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zncoppersmith :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxq_reverse :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vec_equiv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_polint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val vec_reduce :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgxq_reverse :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zc_union_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_indexsort : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_sort : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_sort_inplace : ('kind, 'structure) typ -> unit
val zv_sort_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_sort_uniq : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_sort_uniq_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_union_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val binomial : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val binomial0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val binomialuu : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val cmp_flx : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val cmp_rgx : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val cmp_nodata :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val cmp_prime_ideal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val cmp_prime_over_p : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val cmp_universal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val convol :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gen_cmp_rgx :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val polcyclo_eval :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val dirdiv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val dirmul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val eulerianpol : Signed.long -> Signed.long -> ('kind, 'structure) typ

val gprec_wensure :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gen_indexsort :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_indexsort_uniq :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_search :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  Signed.long

val gen_setminus :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  (('kind, 'structure) typ -> ('kind, 'structure) typ -> int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_sort :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_sort_inplace :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val gen_sort_shallow :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_sort_uniq :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val getstack : unit -> Signed.long
val gettime : unit -> Signed.long
val getabstime : unit -> Signed.long
val getwalltime : unit -> ('kind, 'structure) typ
val gprec : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gprec_wtrunc :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gprec_w : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val indexlexsort : ('kind, 'structure) typ -> ('kind, 'structure) typ
val indexsort : ('kind, 'structure) typ -> ('kind, 'structure) typ

val indexvecsort :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val laplace : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lexsort : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mathilbert : Signed.long -> ('kind, 'structure) typ

val matqpascal :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val merge_factor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val merge_sort_uniq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val modreverse : ('kind, 'structure) typ -> ('kind, 'structure) typ
val polhermite : Signed.long -> Signed.long -> ('kind, 'structure) typ

val polhermite_eval0 :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val polhermite_eval :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val pollaguerre :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val pollaguerre_eval :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val pollaguerre_eval0 :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val pollegendre : Signed.long -> Signed.long -> ('kind, 'structure) typ
val pollegendre_reduced : Signed.long -> Signed.long -> ('kind, 'structure) typ

val pollegendre_eval :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val pollegendre_eval0 :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val polint :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val polint_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val polintspec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val polchebyshev :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val polchebyshev_eval :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val polchebyshev1 : Signed.long -> Signed.long -> ('kind, 'structure) typ
val polchebyshev2 : Signed.long -> Signed.long -> ('kind, 'structure) typ
val polrecip : ('kind, 'structure) typ -> ('kind, 'structure) typ

val setbinop :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val setdelta :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val setintersect :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val setisset : ('kind, 'structure) typ -> Signed.long

val setminus :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val setunion :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val setunion_i :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val sort : ('kind, 'structure) typ -> ('kind, 'structure) typ

val sort_factor :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val stirling :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val stirling1 : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val stirling2 : pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val tablesearch :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  (('kind, 'structure) typ -> ('kind, 'structure) typ -> int)
  Ctypes_static.static_funptr ->
  Signed.long

val vecbinomial : Signed.long -> ('kind, 'structure) typ

val vecsearch :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val vecsort :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecsort0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zv_search : ('kind, 'structure) typ -> Signed.long -> Signed.long

val bits_to_int :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val bits_to_u : ('kind, 'structure) typ -> Signed.long -> pari_ulong
val binaire : ('kind, 'structure) typ -> ('kind, 'structure) typ

val binary_2k :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val binary_2k_nv :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val binary_zv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bittest : ('kind, 'structure) typ -> Signed.long -> Signed.long

val fromdigits_2k :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gbitand :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gbitneg : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gbitnegimply :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gbitor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gbittest : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gbitxor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hammingl : pari_ulong -> Signed.long
val hammingweight : ('kind, 'structure) typ -> Signed.long

val ibitand :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ibitnegimply :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ibitor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ibitxor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nv_fromdigits_2k :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val bnflogef :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnflog :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnflogdegree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfislocalpower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val rnfislocalcyclo : ('kind, 'structure) typ -> Signed.long

val bnfisunit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnfissunit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnfsunit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnfunits :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnfisunit0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val sunits_mod_units :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val buchquad :
  ('kind, 'structure) typ ->
  float ->
  float ->
  Signed.long ->
  ('kind, 'structure) typ

val quadclassno : ('kind, 'structure) typ -> ('kind, 'structure) typ
val quadclassnos : Signed.long -> Signed.long

val quadclassunit0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val buchall :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val buchall_param :
  ('kind, 'structure) typ ->
  float ->
  float ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val bnf_build_cheapfu : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_build_cycgen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_build_matalpha : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_build_units : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_compactfu : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_compactfu_mat : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_has_fu : ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnfinit0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnfisprincipal0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnfnewprec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val bnfnewprec_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val bnftestprimes : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val bnrnewprec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val bnrnewprec_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val isprincipalfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val isprincipalfact_or_fail :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val isprincipal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_cxlog_normalize :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val nfcyclotomicunits :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfsign_units :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int ->
  ('kind, 'structure) typ

val nfsign_tu :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfsign_fu :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val signunits : ('kind, 'structure) typ -> ('kind, 'structure) typ
val hermite_bound : Signed.long -> Signed.long -> ('kind, 'structure) typ

val bnr_subgroup_sanitize :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val bnr_char_sanitize :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val abc_to_bnr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int ->
  ('kind, 'structure) typ

val buchray :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val buchraymod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnrautmatrix :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnr_subgroup_check :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val bnrchar :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnrchar_primitive :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnrclassno :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnrclassno0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnrclassnolist :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnrchar_primitive_raw :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnrconductor_factored :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnrconductor_raw :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnrconductormod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnrconductor0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrconductor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrconductor_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrconductorofchar :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnrdisc0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrdisc :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrdisclist0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnrgaloismatrix :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnrgaloisapply :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnrinit0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrinitmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val bnrisconductor0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val bnrisconductor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val bnrisgalois :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val bnrisprincipalmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrisprincipal :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrmap :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnrsurjection :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnfnarrow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnfcertify : ('kind, 'structure) typ -> Signed.long
val bnfcertify0 : ('kind, 'structure) typ -> Signed.long -> Signed.long

val bnrcompositum :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val decodemodule :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val discrayabslist :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val discrayabslistarch :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val idealmoddivisor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val isprincipalray :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val isprincipalraygen :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nf_deg1_prime : ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfarchstar :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rnfconductor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfconductor0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfnormgroup :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val subgrouplist0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnfisnorm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfisnorm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfisnorminit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int ->
  ('kind, 'structure) typ

val coprimes_zv : pari_ulong -> ('kind, 'structure) typ
val char_check : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val charconj :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val charconj0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val chardiv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val chardiv0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val chareval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val chargalois :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val charker :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val charker0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val charmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val charmul0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val charorder :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val charorder0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val charpow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val charpow0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val char_denormalize :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val char_normalize :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val char_simplify :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val checkznstar_i : ('kind, 'structure) typ -> int
val cyc_normalize : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ncharvecexpo :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znchar : ('kind, 'structure) typ -> ('kind, 'structure) typ

val znchar_quad :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zncharcheck : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val zncharconductor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zncharconj :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znchardecompose :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val znchardiv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val znchareval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zncharinduce :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zncharisodd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val zncharker :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zncharmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zncharorder :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zncharpow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val znchartokronecker :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val znchartoprimitive :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znconrey_check : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val znconrey_normalized :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znconreychar :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znconreyfromchar_normalized :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znconreyconductor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val znconreyexp :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znconreyfromchar :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znconreylog :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znconreylog_normalize :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val znlog0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zv_cyc_minimal :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val zv_cyc_minimize :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val closure_deriv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val closure_derivn :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val localvars_find :
  ('kind, 'structure) typ ->
  entree Ctypes.structure Ctypes_static.ptr ->
  Signed.long

val localvars_read_str :
  string -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val snm_closure :
  entree Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val strtoclosure : string -> Signed.long -> ('kind, 'structure) typ
val strtofunction : string -> ('kind, 'structure) typ

val gconcat :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gconcat1 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val matconcat : ('kind, 'structure) typ -> ('kind, 'structure) typ

val shallowconcat :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val shallowconcat1 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val shallowmatconcat : ('kind, 'structure) typ -> ('kind, 'structure) typ

val vconcat :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val default0 : string -> string -> ('kind, 'structure) typ
val getrealprecision : unit -> Signed.long
val pari_is_default : string -> entree Ctypes.structure Ctypes_static.ptr
val sd_texstyle : string -> Signed.long -> ('kind, 'structure) typ
val sd_colors : string -> Signed.long -> ('kind, 'structure) typ
val sd_compatible : string -> Signed.long -> ('kind, 'structure) typ
val sd_datadir : string -> Signed.long -> ('kind, 'structure) typ
val sd_debug : string -> Signed.long -> ('kind, 'structure) typ
val sd_debugfiles : string -> Signed.long -> ('kind, 'structure) typ
val sd_debugmem : string -> Signed.long -> ('kind, 'structure) typ
val sd_factor_add_primes : string -> Signed.long -> ('kind, 'structure) typ
val sd_factor_proven : string -> Signed.long -> ('kind, 'structure) typ
val sd_format : string -> Signed.long -> ('kind, 'structure) typ
val sd_histsize : string -> Signed.long -> ('kind, 'structure) typ
val sd_log : string -> Signed.long -> ('kind, 'structure) typ
val sd_logfile : string -> Signed.long -> ('kind, 'structure) typ
val sd_nbthreads : string -> Signed.long -> ('kind, 'structure) typ
val sd_new_galois_format : string -> Signed.long -> ('kind, 'structure) typ
val sd_output : string -> Signed.long -> ('kind, 'structure) typ
val sd_parisize : string -> Signed.long -> ('kind, 'structure) typ
val sd_parisizemax : string -> Signed.long -> ('kind, 'structure) typ
val sd_path : string -> Signed.long -> ('kind, 'structure) typ
val sd_plothsizes : string -> Signed.long -> ('kind, 'structure) typ
val sd_prettyprinter : string -> Signed.long -> ('kind, 'structure) typ
val sd_primelimit : string -> Signed.long -> ('kind, 'structure) typ
val sd_realbitprecision : string -> Signed.long -> ('kind, 'structure) typ
val sd_realprecision : string -> Signed.long -> ('kind, 'structure) typ
val sd_secure : string -> Signed.long -> ('kind, 'structure) typ
val sd_seriesprecision : string -> Signed.long -> ('kind, 'structure) typ
val sd_simplify : string -> Signed.long -> ('kind, 'structure) typ
val sd_sopath : string -> int -> ('kind, 'structure) typ
val sd_strictargs : string -> Signed.long -> ('kind, 'structure) typ
val sd_strictmatch : string -> Signed.long -> ('kind, 'structure) typ

val sd_string :
  string ->
  Signed.long ->
  string ->
  string Ctypes_static.ptr ->
  ('kind, 'structure) typ

val sd_threadsize : string -> Signed.long -> ('kind, 'structure) typ
val sd_threadsizemax : string -> Signed.long -> ('kind, 'structure) typ

val sd_intarray :
  string ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  string ->
  ('kind, 'structure) typ

val sd_toggle :
  string ->
  Signed.long ->
  string ->
  int Ctypes_static.ptr ->
  ('kind, 'structure) typ

val sd_ulong :
  string ->
  Signed.long ->
  string ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  pari_ulong ->
  string Ctypes_static.ptr ->
  ('kind, 'structure) typ

val setdefault : string -> string -> Signed.long -> ('kind, 'structure) typ

val setrealprecision :
  Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val digits :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fromdigits :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fromdigitsu :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gen_digits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val gen_fromdigits :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ

val sumdigits : ('kind, 'structure) typ -> ('kind, 'structure) typ

val sumdigits0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val sumdigitsu : pari_ulong -> pari_ulong
val ecpp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ecpp0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ecppexport :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ecppisvalid : ('kind, 'structure) typ -> Signed.long
val isprimeecpp : ('kind, 'structure) typ -> Signed.long
val sd_breakloop : string -> Signed.long -> ('kind, 'structure) typ
val sd_echo : string -> Signed.long -> ('kind, 'structure) typ
val sd_graphcolormap : string -> Signed.long -> ('kind, 'structure) typ
val sd_graphcolors : string -> Signed.long -> ('kind, 'structure) typ
val sd_help : string -> Signed.long -> ('kind, 'structure) typ
val sd_histfile : string -> Signed.long -> ('kind, 'structure) typ
val sd_lines : string -> Signed.long -> ('kind, 'structure) typ
val sd_linewrap : string -> Signed.long -> ('kind, 'structure) typ
val sd_prompt : string -> Signed.long -> ('kind, 'structure) typ
val sd_prompt_cont : string -> Signed.long -> ('kind, 'structure) typ
val sd_psfile : string -> Signed.long -> ('kind, 'structure) typ
val sd_readline : string -> Signed.long -> ('kind, 'structure) typ
val sd_recover : string -> Signed.long -> ('kind, 'structure) typ
val sd_timer : string -> Signed.long -> ('kind, 'structure) typ
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
val gp_alarm : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val gp_input : unit -> ('kind, 'structure) typ
val gp_allocatemem : ('kind, 'structure) typ -> unit
val gp_handle_exception : Signed.long -> int
val gp_alarm_handler : int -> unit
val gp_sigint_fun : unit -> unit
val gp_help : string -> Signed.long -> unit
val gp_echo_and_log : string -> string -> unit
val print_fun_list : string Ctypes_static.ptr -> Signed.long -> unit
val strtime : Signed.long -> ('kind, 'structure) typ

val direuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val dirpowers :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val dirpowerssum :
  pari_ulong ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val dirpowerssumfun :
  pari_ulong ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val vecpowuu : Signed.long -> pari_ulong -> ('kind, 'structure) typ

val vecpowug :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellanalyticrank :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellanalyticrank_bitprec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellanal_globalred_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ellheegner : ('kind, 'structure) typ -> ('kind, 'structure) typ

val elll1 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val elll1_bitprec :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellconvertname : ('kind, 'structure) typ -> ('kind, 'structure) typ
val elldatagenerators : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellidentify : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellsearch : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellsearchcurve : ('kind, 'structure) typ -> ('kind, 'structure) typ

val forell :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  unit

val ellfromeqn : ('kind, 'structure) typ -> ('kind, 'structure) typ

val akell :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bilhell :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val checkell : ('kind, 'structure) typ -> unit
val checkell_fq : ('kind, 'structure) typ -> unit
val checkell_q : ('kind, 'structure) typ -> unit
val checkell_qp : ('kind, 'structure) typ -> unit
val checkellisog : ('kind, 'structure) typ -> unit
val checkellpt : ('kind, 'structure) typ -> unit
val checkell5 : ('kind, 'structure) typ -> unit

val cxredsl2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val cxredsl2_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ec_2divpol_evalx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ec_3divpol_evalx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ec_bmodel :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ec_f_evalx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ec_h_evalx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ec_dfdx_evalq :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ec_dfdy_evalq :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ec_dmfdy_evalq :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ec_half_deriv_2divpol :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ec_half_deriv_2divpol_evalx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ec_phi2 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ell_is_integral : ('kind, 'structure) typ -> int
val ellq_get_cm : ('kind, 'structure) typ -> Signed.long
val ellq_get_n : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellq_get_nfa :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val ellqp_tate_uniformization :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellqp_agm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellqp_u : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ellqp_u2 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ellqp_q : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ellqp_ab : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ellqp_l : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellqp_root :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellqtwist_bsdperiod :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellr_area :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellr_ab : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ellr_eta : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellr_omega :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellr_roots :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellan : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellanq_zv :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellanal_globalred :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ellap :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellap_cm_fast :
  ('kind, 'structure) typ -> pari_ulong -> Signed.long -> Signed.long

val ellbasechar : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellbsd : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellchangecurve :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellchangeinvert : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellchangepoint :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellchangepointinv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val elleisnum :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val elleta : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val elleulerf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellff_get_card : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellff_get_gens : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellff_get_group : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellff_get_o : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellff_get_p : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellff_get_m : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellff_get_d : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellfromj : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellgenerators : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellglobalred : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellgroup :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellgroup0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellheight0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellheight :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellheightmatrix :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellheightoo :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellintegralmodel :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ellintegralmodel_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val elliscm : ('kind, 'structure) typ -> Signed.long

val ellisoncurve :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellisotree : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellissupersingular :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val elljissupersingular : ('kind, 'structure) typ -> int

val elllseries :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val elllocalred :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val elllog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellminimaldisc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellminimalmodel :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ellminimaltwist : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellminimaltwist0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellminimaltwistcond : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellnf_vecarea :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellnf_veceta :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellnf_vecomega :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellneg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellorder :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellorder_q :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val ellordinate :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellpadicheight0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellpadicheightmatrix :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellperiods :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellrootno :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val ellrootno_global : ('kind, 'structure) typ -> Signed.long

val ellsaturation :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellsea : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellsigma :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellsupersingularj : ('kind, 'structure) typ -> ('kind, 'structure) typ
val elltamagawa : ('kind, 'structure) typ -> ('kind, 'structure) typ

val elltaniyama :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val elltatepairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val elltors : ('kind, 'structure) typ -> ('kind, 'structure) typ
val elltors0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val elltors_psylow :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val elltrace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val elltwist :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellwp :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellwp0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellwpseries :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellxn :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellzeta :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val oncurve : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val orderell :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val pointell :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val point_to_a4a6 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val point_to_a4a6_fl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zell :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qp_agm2_sequence :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qp_ascending_landen :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val qp_descending_landen :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val ellformaldifferential :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellformalexp :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellformallog :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellformalpoint :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellformalw :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellnonsingularmultiple :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellpadicl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellpadicbsd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellpadicfrobenius :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val ellpadicheight :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellpadiclog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellpadicregulator :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellpadics2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ell2cover :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ellrank :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ellrankinit :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val hyperell_locally_soluble :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val nf_hyperell_locally_soluble :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val nfhilbert :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val nfhilbert0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val ellisdivisible :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val ellisogenyapply :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellisogeny :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellisomat :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ellweilcurve :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val flxq_elldivpolmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val fp_ellcard_sea :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fq_ellcard_sea :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fq_elldivpolmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ellmodulareqn :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val externstr : string -> ('kind, 'structure) typ
val gp_filter : string -> string
val gpextern : string -> ('kind, 'structure) typ
val gpsystem : string -> Signed.long
val readstr : string -> ('kind, 'structure) typ
val gentogenstr_nospace : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gentogenstr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gentotexstr : ('kind, 'structure) typ -> string
val gentostr : ('kind, 'structure) typ -> string
val gentostr_raw : ('kind, 'structure) typ -> string
val gentostr_unquoted : ('kind, 'structure) typ -> string
val str : ('kind, 'structure) typ -> ('kind, 'structure) typ
val strexpand : ('kind, 'structure) typ -> ('kind, 'structure) typ
val strtex : ('kind, 'structure) typ -> ('kind, 'structure) typ
val brute : ('kind, 'structure) typ -> char -> Signed.long -> unit
val dbggen : ('kind, 'structure) typ -> Signed.long -> unit
val error0 : ('kind, 'structure) typ -> unit
val dbg_pari_heap : unit -> unit
val err_flush : unit -> unit
val err_printf : string -> unit
val gp_getenv : string -> ('kind, 'structure) typ
val gp_fileclose : Signed.long -> unit
val gp_fileextern : string -> Signed.long
val gp_fileflush : Signed.long -> unit
val gp_fileflush0 : ('kind, 'structure) typ -> unit
val gp_fileopen : string -> string -> Signed.long
val gp_fileread : Signed.long -> ('kind, 'structure) typ
val gp_filereadstr : Signed.long -> ('kind, 'structure) typ
val gp_filewrite : Signed.long -> string -> unit
val gp_filewrite1 : Signed.long -> string -> unit
val gp_read_file : string -> ('kind, 'structure) typ
val gp_read_str_multiline : string -> string -> ('kind, 'structure) typ
val gp_readvec_file : string -> ('kind, 'structure) typ
val gpinstall : string -> string -> string -> string -> unit
val gsprintf : string -> ('kind, 'structure) typ
val itostr : ('kind, 'structure) typ -> string
val matbrute : ('kind, 'structure) typ -> char -> Signed.long -> unit
val os_getenv : string -> string
val uordinal : pari_ulong -> string
val outmat : ('kind, 'structure) typ -> unit
val output : ('kind, 'structure) typ -> unit
val rgv_to_str : ('kind, 'structure) typ -> Signed.long -> string

val pari_add_hist :
  ('kind, 'structure) typ -> Signed.long -> Signed.long -> unit

val pari_ask_confirm : string -> unit
val pari_flush : unit -> unit
val pari_get_hist : Signed.long -> ('kind, 'structure) typ
val pari_get_histrtime : Signed.long -> Signed.long
val pari_get_histtime : Signed.long -> Signed.long
val pari_get_homedir : string -> string
val pari_histtime : Signed.long -> ('kind, 'structure) typ
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
val pari_sprint0 : string -> ('kind, 'structure) typ -> Signed.long -> string
val print : ('kind, 'structure) typ -> unit
val printp : ('kind, 'structure) typ -> unit
val print1 : ('kind, 'structure) typ -> unit
val printf0 : string -> ('kind, 'structure) typ -> unit
val printsep : string -> ('kind, 'structure) typ -> unit
val printsep1 : string -> ('kind, 'structure) typ -> unit
val printtex : ('kind, 'structure) typ -> unit
val stack_sprintf : string -> string
val str_init : pari_str Ctypes.structure Ctypes_static.ptr -> int -> unit
val str_printf : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val str_putc : pari_str Ctypes.structure Ctypes_static.ptr -> char -> unit
val str_puts : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val strftime_expand : string -> string -> Signed.long -> unit
val strprintf : string -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val term_color : Signed.long -> unit
val term_get_color : string -> Signed.long -> string
val texe : ('kind, 'structure) typ -> char -> Signed.long -> unit
val warning0 : ('kind, 'structure) typ -> unit
val write0 : string -> ('kind, 'structure) typ -> unit
val write1 : string -> ('kind, 'structure) typ -> unit
val writebin : string -> ('kind, 'structure) typ -> unit
val writetex : string -> ('kind, 'structure) typ -> unit
val bincopy_relink : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val bitprecision0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val bitprecision00 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val break0 : Signed.long -> ('kind, 'structure) typ

val call0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val closure_callgen0prec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val closure_callgen1 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val closure_callgen1prec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val closure_callgen2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val closure_callgenall :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val closure_callgenvec :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val closure_callgenvecdef :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val closure_callgenvecdefprec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val closure_callgenvecprec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val closure_callvoid1 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val closure_context : Signed.long -> Signed.long -> Signed.long
val closure_disassemble : ('kind, 'structure) typ -> unit
val closure_err : Signed.long -> unit

val closure_evalbrk :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val closure_evalgen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val closure_evalnobrk : ('kind, 'structure) typ -> ('kind, 'structure) typ
val closure_evalres : ('kind, 'structure) typ -> ('kind, 'structure) typ
val closure_evalvoid : ('kind, 'structure) typ -> unit
val closure_func_err : unit -> string

val closure_trapgen :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val copybin_unlink : ('kind, 'structure) typ -> ('kind, 'structure) typ
val getlocalprec : Signed.long -> Signed.long
val getlocalbitprec : Signed.long -> Signed.long
val get_lex : Signed.long -> ('kind, 'structure) typ
val get_localprec : unit -> Signed.long
val get_localbitprec : unit -> Signed.long

val gp_call :
  unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gp_callprec :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gp_call2 :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gp_callbool :
  unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long

val gp_callvoid :
  unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long

val gp_eval :
  unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gp_evalbool :
  unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long

val gp_evalprec :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gp_evalupto :
  unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gp_evalvoid :
  unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long

val localprec : ('kind, 'structure) typ -> unit
val localbitprec : ('kind, 'structure) typ -> unit
val loop_break : unit -> Signed.long
val next0 : Signed.long -> ('kind, 'structure) typ
val pareval : ('kind, 'structure) typ -> ('kind, 'structure) typ
val pari_self : unit -> ('kind, 'structure) typ

val parsum :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val parvector :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val pop_lex : Signed.long -> unit
val pop_localprec : unit -> unit

val precision0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val precision00 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val push_lex : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val push_localbitprec : Signed.long -> unit
val push_localprec : Signed.long -> unit
val return0 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val set_lex : Signed.long -> ('kind, 'structure) typ -> unit

val forcomposite_init :
  forcomposite_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val forcomposite_next :
  forcomposite_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val forprime_next :
  forprime_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val forprime_init :
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val forprimestep_init :
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
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
val prodprimes : unit -> ('kind, 'structure) typ

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

val ff_frobenius :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ff_z_z_muldiv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ff_q_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_z_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_add :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_charpoly : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_conjvec : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_ellcard : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_ellcard_sea :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ff_ellgroup :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ff_elllog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ff_ellmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ff_ellorder :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ff_elltwist : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_ellrandom : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_elltatepairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ff_ellweilpairing :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ff_equal0 : ('kind, 'structure) typ -> int
val ff_equal1 : ('kind, 'structure) typ -> int
val ff_equalm1 : ('kind, 'structure) typ -> int
val ff_f : ('kind, 'structure) typ -> Signed.long
val ff_gen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_inv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_issquare : ('kind, 'structure) typ -> Signed.long

val ff_issquareall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val ff_ispower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val ff_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ff_map :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_minpoly : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_mul2n : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ff_neg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_neg_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_norm : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_order :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_p_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_pow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_q : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_samefield : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val ff_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_sqrt : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_sqrtn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ff_sub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ff_to_f2xq : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_to_f2xq_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_to_flxq : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_to_flxq_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_to_fpxq : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_trace : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ff_var : ('kind, 'structure) typ -> Signed.long

val ffm_ffc_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffm_ffc_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffm_ffc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffm_deplin :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffm_det :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffm_gauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffm_image :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffm_indexrank :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffm_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffm_invimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffm_ker :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffm_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffm_rank : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val ffm_suppl :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffx_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffx_ddf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffx_degfact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffx_disc :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffx_extgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ffx_factor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffx_factor_squarefree :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffx_halfgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffx_halfgcd_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ffx_ispower :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val ffx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffx_preimage :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffx_preimagerel :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffx_rem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffx_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffx_roots :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffx_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffxq_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffxq_minpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffxq_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_to_ffx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fq_to_ff :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_ff_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffembed :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fffrobenius :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ffinvmap : ('kind, 'structure) typ -> ('kind, 'structure) typ

val fflog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val ffmap :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffmaprel :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffcompomap :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fforder :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffprimroot :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ffrandom : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rg_is_ff :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> int

val rgc_is_ffc :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> int

val rgm_is_ffm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> int

val p_to_ff : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val tp_to_ff :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val flx_factcyclo :
  pari_ulong -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val fpx_factcyclo :
  pari_ulong -> ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val factormodcyclo :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val checkgal : ('kind, 'structure) typ -> ('kind, 'structure) typ

val checkgroup :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val checkgroupelts : ('kind, 'structure) typ -> ('kind, 'structure) typ

val embed_disc :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val embed_roots :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val galois_group : ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoisconj :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoisconj0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val galoisconjclasses : ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoisexport :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val galoisfixedfield :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val galoisidentify : ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoisinit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoisisabelian :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val galoisisnormal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val galoispermtopol :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoissplittinginit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoissubgroups : ('kind, 'structure) typ -> ('kind, 'structure) typ

val galoissubfields :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val numberofconjugates : ('kind, 'structure) typ -> Signed.long -> Signed.long

val polgalois :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val galoisnbpol : Signed.long -> ('kind, 'structure) typ
val galoisgetgroup : Signed.long -> Signed.long -> ('kind, 'structure) typ
val galoisgetname : Signed.long -> Signed.long -> ('kind, 'structure) typ

val galoisgetpol :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val conj_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val conjvec : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val divrunextu :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val gadd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gaddsg : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val gconj : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gdiv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gdivgs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gdivgu : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val gdivgunextu :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val ginv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmul2n : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gmulsg : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val gmulug : pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val gsqr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gsub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gsubsg : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val mulcxi : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mulcxmi : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mulcxpowis :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val qdivii :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qdiviu : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val qdivis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ser_normalize : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gassoc_proto :
  (('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val map_proto_g :
  (('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val map_proto_lg :
  (('kind, 'structure) typ -> Signed.long) Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val map_proto_lgl :
  (('kind, 'structure) typ -> Signed.long -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val q_lval : ('kind, 'structure) typ -> pari_ulong -> Signed.long

val q_lvalrem :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val q_pval : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val q_pvalrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val rgx_val : ('kind, 'structure) typ -> Signed.long

val rgx_valrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val rgx_valrem_inexact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val rgxv_maxdegree : ('kind, 'structure) typ -> Signed.long
val zv_z_dvd : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val zv_pval : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val zv_pvalrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val zv_lval : ('kind, 'structure) typ -> pari_ulong -> Signed.long

val zv_lvalrem :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val zx_lvalrem :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val zx_pval : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val zx_pvalrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val z_lvalrem_stop :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong ->
  int Ctypes_static.ptr ->
  Signed.long

val cgetp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val cvstop2 : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val cvtop :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val cvtop2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val cx_approx_equal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val cx_approx0 : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val gabs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gaffect : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val gaffsg : Signed.long -> ('kind, 'structure) typ -> unit
val gcmp : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val gequal0 : ('kind, 'structure) typ -> int
val gequal1 : ('kind, 'structure) typ -> int
val gequalx : ('kind, 'structure) typ -> int
val gequalm1 : ('kind, 'structure) typ -> int
val gcmpsg : Signed.long -> ('kind, 'structure) typ -> int

val gcvtop :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gequal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> bool
val gequalsg : Signed.long -> ('kind, 'structure) typ -> int
val gexpo : ('kind, 'structure) typ -> Signed.long
val gexpo_safe : ('kind, 'structure) typ -> Signed.long
val gpexponent : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gpvaluation :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gvaluation :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val gidentical : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val gmax :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmaxgs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gmin :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmings : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gneg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gneg_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gsigne : ('kind, 'structure) typ -> int
val gtolist : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gtolong : ('kind, 'structure) typ -> Signed.long
val lexcmp : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val listinsert :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val listpop : ('kind, 'structure) typ -> Signed.long -> unit
val listpop0 : ('kind, 'structure) typ -> Signed.long -> unit

val listput :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val listput0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long -> unit

val listsort : ('kind, 'structure) typ -> Signed.long -> unit
val matsize : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mklist : unit -> ('kind, 'structure) typ
val mklist_typ : Signed.long -> ('kind, 'structure) typ
val mklistcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mkmap : unit -> ('kind, 'structure) typ
val normalizeser : ('kind, 'structure) typ -> ('kind, 'structure) typ
val normalizepol : ('kind, 'structure) typ -> ('kind, 'structure) typ

val normalizepol_approx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val padic_to_fl : ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val padic_to_fp :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val quadtofp : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val sizedigit : ('kind, 'structure) typ -> Signed.long
val u_lval : pari_ulong -> pari_ulong -> Signed.long

val u_lvalrem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val u_lvalrem_stop :
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  int Ctypes_static.ptr ->
  Signed.long

val u_pval : pari_ulong -> ('kind, 'structure) typ -> Signed.long

val u_pvalrem :
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong Ctypes_static.ptr ->
  Signed.long

val vecindexmax : ('kind, 'structure) typ -> Signed.long
val vecindexmin : ('kind, 'structure) typ -> Signed.long

val vecmax0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val vecmax : ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecmin0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val vecmin : ('kind, 'structure) typ -> ('kind, 'structure) typ
val z_lval : Signed.long -> pari_ulong -> Signed.long

val z_lvalrem :
  Signed.long -> pari_ulong -> Signed.long Ctypes_static.ptr -> Signed.long

val z_pval : Signed.long -> ('kind, 'structure) typ -> Signed.long

val z_pvalrem :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val zx_lval : ('kind, 'structure) typ -> Signed.long -> Signed.long
val hgmcyclo : ('kind, 'structure) typ -> ('kind, 'structure) typ
val hgmalpha : ('kind, 'structure) typ -> ('kind, 'structure) typ
val hgmgamma : ('kind, 'structure) typ -> ('kind, 'structure) typ

val hgminit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hgmparams : ('kind, 'structure) typ -> ('kind, 'structure) typ

val hgmeulerfactor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val hgmcoef :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val hgmcoefs :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val hgmtwist : ('kind, 'structure) typ -> ('kind, 'structure) typ
val hgmissymmetrical : ('kind, 'structure) typ -> Signed.long
val hgmbydegree : Signed.long -> ('kind, 'structure) typ

val lfunhgm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qp_zeta : ('kind, 'structure) typ -> ('kind, 'structure) typ

val lerchphi :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val lerchzeta :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zetahurwitz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rgx_to_ser :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_to_ser_inexact :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gtoser :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val gtoser_prec :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rfrac_to_ser :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rfrac_to_ser_i :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rfracrecip_to_ser_absolute :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rfracrecip :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val scalarser :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val sertoser : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val toser_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgv_to_ser :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ser0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val padic_to_q : ('kind, 'structure) typ -> ('kind, 'structure) typ
val padic_to_q_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qpv_to_qv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_rgv_mulrealsym :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_mulreal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_cxeval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgx_deflate_max :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgx_deflate_order : ('kind, 'structure) typ -> Signed.long
val rgx_degree : ('kind, 'structure) typ -> Signed.long -> Signed.long
val rgx_integ : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxy_cxevalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zx_deflate_order : ('kind, 'structure) typ -> Signed.long

val zx_deflate_max :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ceil_safe : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ceilr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val centerlift : ('kind, 'structure) typ -> ('kind, 'structure) typ

val centerlift0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val compo : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val deg1pol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val deg1pol_shallow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val deg2pol_shallow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val denom : ('kind, 'structure) typ -> ('kind, 'structure) typ
val denom_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val denominator :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val derivser : ('kind, 'structure) typ -> ('kind, 'structure) typ

val diffop :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val diffop0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val diviiround :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val floor_safe : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gceil : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gcvtoi :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gdeflate :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val gdivent :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gdiventgs :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gdiventsg :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gdiventres :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gdivmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gdivround :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gdvd : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val geq :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val geval : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gfloor : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gtrunc2n : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gfrac : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gge :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ggrando : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ggt :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gimag : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gisexactzero : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gle :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val glt :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmodgs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gmodsg : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmodulo :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmodulsg : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val gmodulss : Signed.long -> Signed.long -> ('kind, 'structure) typ

val gne :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gnot : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gpolvar : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gppadicprec :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gppoldegree :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gprecision : ('kind, 'structure) typ -> Signed.long

val gpserprec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val greal : ('kind, 'structure) typ -> ('kind, 'structure) typ

val grndtoi :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ground : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gshift : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gsubst :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gsubstpol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gsubstvec :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gtocol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gtocol0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gtocolrev : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gtocolrev0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gtopoly : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gtopolyrev :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gtovec : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gtovec0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gtovecrev : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gtovecrev0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gtovecsmall : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gtovecsmall0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gtrunc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gvar2 : ('kind, 'structure) typ -> Signed.long

val hqfeval :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val imag_i : ('kind, 'structure) typ -> ('kind, 'structure) typ
val integ : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val integser : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ser_inv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val iscomplex : ('kind, 'structure) typ -> int
val isexactzero : ('kind, 'structure) typ -> int
val isrationalzeroscalar : ('kind, 'structure) typ -> int
val isinexact : ('kind, 'structure) typ -> int
val isinexactreal : ('kind, 'structure) typ -> int

val isint :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> int

val isrationalzero : ('kind, 'structure) typ -> int
val issmall : ('kind, 'structure) typ -> Signed.long Ctypes_static.ptr -> int
val lift : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lift_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lift0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val liftall : ('kind, 'structure) typ -> ('kind, 'structure) typ
val liftall_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val liftint : ('kind, 'structure) typ -> ('kind, 'structure) typ
val liftint_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val liftpol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val liftpol_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mkcoln : Signed.long -> ('kind, 'structure) typ
val mkintn : Signed.long -> ('kind, 'structure) typ
val mkpoln : Signed.long -> ('kind, 'structure) typ
val mkvecn : Signed.long -> ('kind, 'structure) typ
val mkvecsmalln : Signed.long -> ('kind, 'structure) typ

val modrr_safe :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val modrr_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mulreal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val numer : ('kind, 'structure) typ -> ('kind, 'structure) typ
val numer_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val numerator :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val padicprec :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val padicprec_relative : ('kind, 'structure) typ -> Signed.long

val polcoef :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val polcoef_i :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val poldegree : ('kind, 'structure) typ -> Signed.long -> Signed.long
val pollead : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val precision : ('kind, 'structure) typ -> Signed.long

val qf_apply_rgm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qf_apply_zm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfb_apply_zm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfbil :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfeval :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfeval0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfevalb :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfnorm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val real_i : ('kind, 'structure) typ -> ('kind, 'structure) typ

val round0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val roundr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val roundr_safe : ('kind, 'structure) typ -> ('kind, 'structure) typ

val scalarpol :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val scalarpol_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ser_unscale :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val serprec : ('kind, 'structure) typ -> Signed.long -> Signed.long
val serreverse : ('kind, 'structure) typ -> ('kind, 'structure) typ
val simplify : ('kind, 'structure) typ -> ('kind, 'structure) typ
val simplify_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ

val tayl :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val trunc0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val uu32toi : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val uu32toineg : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val vars_sort_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vars_to_rgxv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val variables_vecsmall : ('kind, 'structure) typ -> ('kind, 'structure) typ
val variables_vec : ('kind, 'structure) typ -> ('kind, 'structure) typ

val genus2red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val genus2igusa :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gchar_conductor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gchar_identify :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gcharalgebraic :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gcharduallog :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gchareval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gchari_lfun :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gcharinit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gcharisalgebraic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val gcharlocal :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gcharlog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gcharnewprec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val is_gchar_group : ('kind, 'structure) typ -> int

val lfungchar :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecan_gchar :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val eulerf_gchar :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val group_ident :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val group_ident_trans :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

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
  ('kind, 'structure) typ

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
  (('kind, 'structure) typ -> ('kind, 'structure) typ -> int)
  Ctypes_static.static_funptr ->
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

val hash_keys :
  hashtable Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val hash_values :
  hashtable Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

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
val hash_gen : ('kind, 'structure) typ -> pari_ulong
val hash_zv : ('kind, 'structure) typ -> pari_ulong

val zx_hyperellred :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val hyperellcharpoly : ('kind, 'structure) typ -> ('kind, 'structure) typ

val hyperellchangecurve :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hyperelldisc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val hyperellisoncurve :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val hyperellminimaldisc :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hyperellminimalmodel :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val hyperellpadicfrobenius0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val hyperellpadicfrobenius :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val hyperellred :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val nfhyperellpadicfrobenius :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val hypergeom :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val airy : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_hnfall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ

val zm_hnf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zm_hnf_knapsack : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_hnfall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ

val zm_hnfall_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ

val zm_hnfcenter : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_hnflll :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int ->
  ('kind, 'structure) typ

val zv_extgcd : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_snfall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zv_snf_group :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zv_snf_rank_u : ('kind, 'structure) typ -> pari_ulong -> Signed.long
val zv_snf_trunc : ('kind, 'structure) typ -> unit

val zm_hnfmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_hnfmodall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zm_hnfmodall_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zm_hnfmodid :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_hnfmodprime :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_hnfperm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zm_snfclean :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val zm_snf : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_snf_group :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zm_snfall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zm_snfall_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ

val zv_snfclean : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zpm_echelon :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gsmith : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gsmithall : ('kind, 'structure) typ -> ('kind, 'structure) typ
val hnf : ('kind, 'structure) typ -> ('kind, 'structure) typ

val hnf_divscale :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val hnf_invscale :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hnf_solve :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hnf_invimage :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hnfall : ('kind, 'structure) typ -> ('kind, 'structure) typ
val hnfdivide : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val hnflll : ('kind, 'structure) typ -> ('kind, 'structure) typ

val hnfmerge_get_1 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hnfmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hnfmodid :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val hnfperm : ('kind, 'structure) typ -> ('kind, 'structure) typ

val matfrobenius :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mathnf0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val matsnf0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val smith : ('kind, 'structure) typ -> ('kind, 'structure) typ
val smithall : ('kind, 'structure) typ -> ('kind, 'structure) typ
val smithclean : ('kind, 'structure) typ -> ('kind, 'structure) typ
val snfrank : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val zlm_echelon :
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val zv_snf_rank : ('kind, 'structure) typ -> pari_ulong -> Signed.long

val z_ecm :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val z_factor : ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_factor_limit :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val z_factor_until :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_issmooth : ('kind, 'structure) typ -> pari_ulong -> Signed.long

val z_issmooth_fact :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val z_issquarefree : ('kind, 'structure) typ -> Signed.long

val z_pollardbrent :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val absz_factor : ('kind, 'structure) typ -> ('kind, 'structure) typ

val absz_factor_limit :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val absz_factor_limit_strict :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val coreu : pari_ulong -> pari_ulong
val coreu_fact : ('kind, 'structure) typ -> pari_ulong

val factorint :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val factoru : pari_ulong -> ('kind, 'structure) typ
val tridiv_boundu : pari_ulong -> pari_ulong
val ifac_isprime : ('kind, 'structure) typ -> int

val ifac_next :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  int

val ifac_read :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  int

val ifac_skip : ('kind, 'structure) typ -> unit
val ifac_start : ('kind, 'structure) typ -> int -> ('kind, 'structure) typ

val is_357_power :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  int

val is_pth_power :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  pari_ulong ->
  int

val ispowerful : ('kind, 'structure) typ -> Signed.long
val maxomegau : pari_ulong -> Signed.long
val maxomegaoddu : pari_ulong -> Signed.long
val moebius : ('kind, 'structure) typ -> Signed.long
val moebiusu : pari_ulong -> Signed.long
val moebiusu_fact : ('kind, 'structure) typ -> Signed.long
val nextprime : ('kind, 'structure) typ -> ('kind, 'structure) typ
val precprime : ('kind, 'structure) typ -> ('kind, 'structure) typ
val radicalu : pari_ulong -> pari_ulong
val tridiv_bound : ('kind, 'structure) typ -> pari_ulong

val uis_357_power :
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  int

val uis_357_powermod : pari_ulong -> pari_ulong Ctypes_static.ptr -> int
val unextprime : pari_ulong -> pari_ulong
val uprecprime : pari_ulong -> pari_ulong
val vecfactorsquarefreeu : pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val vecfactorsquarefreeu_coprime :
  pari_ulong -> pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecfactoru_i : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val vecfactoru : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val vecfactoroddu_i : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val vecfactoroddu : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val vecsquarefreeu : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val chk_gerepileupto : ('kind, 'structure) typ -> int

val copy_bin :
  ('kind, 'structure) typ -> genbin Ctypes.structure Ctypes_static.ptr

val copy_bin_canon :
  ('kind, 'structure) typ -> genbin Ctypes.structure Ctypes_static.ptr

val dbg_gerepile : pari_ulong -> unit
val dbg_gerepileupto : ('kind, 'structure) typ -> unit
val errname : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gclone : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gcloneref : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gclone_refc : ('kind, 'structure) typ -> unit
val gcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gcopy_avma :
  ('kind, 'structure) typ ->
  pari_ulong Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gcopy_lg : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gerepile :
  pari_ulong -> pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gerepileallsp : pari_ulong -> pari_ulong -> int -> unit

val gerepilecoeffssp :
  pari_ulong -> pari_ulong -> Signed.long Ctypes_static.ptr -> int -> unit

val gerepilemanysp :
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr Ctypes_static.ptr ->
  int ->
  unit

val getheap : unit -> ('kind, 'structure) typ
val gsizeword : ('kind, 'structure) typ -> Signed.long
val gsizebyte : ('kind, 'structure) typ -> Signed.long
val gunclone : ('kind, 'structure) typ -> unit
val gunclone_deep : ('kind, 'structure) typ -> unit
val listcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val listinit : ('kind, 'structure) typ -> ('kind, 'structure) typ
val msgtimer : string -> unit
val name_numerr : string -> Signed.long
val new_chunk_resize : int -> unit
val newblock : int -> ('kind, 'structure) typ
val numerr_name : Signed.long -> string

val obj_check :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val obj_checkbuild :
  ('kind, 'structure) typ ->
  Signed.long ->
  (('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val obj_checkbuild_padicprec :
  ('kind, 'structure) typ ->
  Signed.long ->
  (('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) typ

val obj_checkbuild_realprec :
  ('kind, 'structure) typ ->
  Signed.long ->
  (('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) typ

val obj_checkbuild_prec :
  ('kind, 'structure) typ ->
  Signed.long ->
  (('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (('kind, 'structure) typ -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) typ

val obj_free : ('kind, 'structure) typ -> unit
val obj_init : Signed.long -> Signed.long -> ('kind, 'structure) typ

val obj_insert :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val obj_insert_shallow :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val obj_reinit : ('kind, 'structure) typ -> ('kind, 'structure) typ
val pari_add_function : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_defaults_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_close : unit -> unit
val pari_close_opts : pari_ulong -> unit
val pari_compile_str : string -> ('kind, 'structure) typ
val pari_daemon : unit -> int
val pari_err : int -> unit
val pari_err_last : unit -> ('kind, 'structure) typ
val pari_err2str : ('kind, 'structure) typ -> string
val pari_init_opts : int -> pari_ulong -> pari_ulong -> unit
val pari_init : int -> pari_ulong -> unit
val pari_stackcheck_init : unit Ctypes_static.ptr -> unit
val pari_sighandler : int -> unit
val pari_sig_init : (int -> unit) Ctypes_static.static_funptr -> unit

val pari_thread_alloc :
  pari_thread Ctypes.structure Ctypes_static.ptr ->
  int ->
  ('kind, 'structure) typ ->
  unit

val pari_thread_close : unit -> unit
val pari_thread_free : pari_thread Ctypes.structure Ctypes_static.ptr -> unit
val pari_thread_init : unit -> unit

val pari_thread_start :
  pari_thread Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val pari_thread_valloc :
  pari_thread Ctypes.structure Ctypes_static.ptr ->
  int ->
  int ->
  ('kind, 'structure) typ ->
  unit

val pari_version : unit -> ('kind, 'structure) typ
val pari_warn : int -> unit
val paristack_newrsize : pari_ulong -> unit
val paristack_resize : pari_ulong -> unit
val paristack_setsize : int -> int -> unit
val parivstack_resize : pari_ulong -> unit
val parivstack_reset : unit -> unit
val setalldebug : Signed.long -> unit
val setdebug : string -> Signed.long -> ('kind, 'structure) typ
val shiftaddress : ('kind, 'structure) typ -> Signed.long -> unit
val shiftaddress_canon : ('kind, 'structure) typ -> Signed.long -> unit
val timer : unit -> Signed.long
val timer_delay : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long
val timer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val timer_printf :
  pari_timer Ctypes.structure Ctypes_static.ptr -> string -> unit

val timer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit
val timer2 : unit -> Signed.long

val trap0 :
  string ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val traverseheap :
  (('kind, 'structure) typ -> unit Ctypes_static.ptr -> unit)
  Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  unit

val walltimer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit

val walltimer_delay :
  pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val walltimer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val contfraceval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val contfracinit :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val intcirc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val intfuncinit :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val intnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val intnumgauss :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val intnumgaussinit : Signed.long -> Signed.long -> ('kind, 'structure) typ

val intnuminit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val intnumosc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val intnumromb :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val intnumromb_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val prodeulerrat :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val prodnumrat :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val quodif : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val sumeulerrat :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val sumnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumnumap :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumnumapinit :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val sumnuminit :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val sumnumlagrangeinit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumnumlagrange :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumnummonien :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumnummonieninit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumnumrat :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumnumsidi :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  float ->
  Signed.long ->
  ('kind, 'structure) typ

val z_isanypower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val z_ispow2 : ('kind, 'structure) typ -> Signed.long

val z_ispowerall :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val z_issquareall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val zn_ispower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val zn_issquare :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val zp_issquare :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val gisanypower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val gissquare : ('kind, 'structure) typ -> ('kind, 'structure) typ

val gissquareall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val ispolygonal :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val ispower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val isprimepower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val ispseudoprimepower :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val issquare : ('kind, 'structure) typ -> Signed.long

val issquareall :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long

val sqrtint : ('kind, 'structure) typ -> ('kind, 'structure) typ

val sqrtint0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val uisprimepower : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long
val uissquare : pari_ulong -> Signed.long
val uissquareall : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val ulogintall :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val padicfields0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val padicfields :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrclassfield :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rnfkummer :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val is_linit : ('kind, 'structure) typ -> Signed.long
val ldata_get_an : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ldata_get_dual : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ldata_get_gammavec : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ldata_get_degree : ('kind, 'structure) typ -> Signed.long
val ldata_get_k : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ldata_get_k1 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ldata_get_conductor : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ldata_get_rootno : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ldata_get_residue : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ldata_get_type : ('kind, 'structure) typ -> Signed.long
val ldata_isreal : ('kind, 'structure) typ -> Signed.long
val linit_get_type : ('kind, 'structure) typ -> Signed.long
val linit_get_ldata : ('kind, 'structure) typ -> ('kind, 'structure) typ
val linit_get_tech : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_domain : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_dom : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_factgammavec : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_step : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_pol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_residue : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_k2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_w2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_expot : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfun_get_bitprec : ('kind, 'structure) typ -> Signed.long

val lfun :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val lfun0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfuncheckfeq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long

val lfunconductor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfuncost :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfuncost0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfuncreate : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfundual : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lfuneuler :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunparams :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lfunan :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunhardy :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val lfuninit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfuninit0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfuninit_make :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val lfunlambda :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunlambda0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunmisc_to_ldata : ('kind, 'structure) typ -> ('kind, 'structure) typ

val lfunmisc_to_ldata_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ

val lfunmisc_to_ldata_shallow_i :
  ('kind, 'structure) typ -> ('kind, 'structure) typ

val lfunorderzero :
  ('kind, 'structure) typ -> Signed.long -> Signed.long -> Signed.long

val lfunprod_get_fact : ('kind, 'structure) typ -> ('kind, 'structure) typ

val lfunrootno :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lfunrootres :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lfunrtopoles : ('kind, 'structure) typ -> ('kind, 'structure) typ

val lfunshift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfuntwist :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val lfuntheta :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunthetacost0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long

val lfunthetacost :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long

val lfunthetainit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunthetacheckinit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunzeros :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val sdomain_isincl :
  float -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val theta_get_an : ('kind, 'structure) typ -> ('kind, 'structure) typ
val theta_get_k : ('kind, 'structure) typ -> ('kind, 'structure) typ
val theta_get_r : ('kind, 'structure) typ -> ('kind, 'structure) typ
val theta_get_bitprec : ('kind, 'structure) typ -> Signed.long
val theta_get_m : ('kind, 'structure) typ -> Signed.long
val theta_get_tdom : ('kind, 'structure) typ -> ('kind, 'structure) typ
val theta_get_isqrtn : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vgaeasytheta : ('kind, 'structure) typ -> int

val znchargauss :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val dirzetak :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ellmoddegree : ('kind, 'structure) typ -> ('kind, 'structure) typ
val eta_zxn : Signed.long -> Signed.long -> ('kind, 'structure) typ

val eta_product_zxn :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val etaquotype :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val galois_get_conj : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ldata_vecan :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ldata_newprec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lfunabelianrelinit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunartin :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val lfundiv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunellmfpeters :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lfunetaquo : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lfungenus2 : ('kind, 'structure) typ -> ('kind, 'structure) typ

val lfunmfspec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lfunmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val lfunqf : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val lfunsympow :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val lfunzetakinit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val qfiseven : ('kind, 'structure) typ -> Signed.long
val lfunquadneg : Signed.long -> Signed.long -> ('kind, 'structure) typ

val zm_lll_norms :
  ('kind, 'structure) typ ->
  float ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val kerint : ('kind, 'structure) typ -> ('kind, 'structure) typ

val lllfp :
  ('kind, 'structure) typ -> float -> Signed.long -> ('kind, 'structure) typ

val lllgen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllgram : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllgramgen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllgramint : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllgramkerim : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllgramkerimgen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllint : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllintpartial : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllintpartial_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllkerim : ('kind, 'structure) typ -> ('kind, 'structure) typ
val lllkerimgen : ('kind, 'structure) typ -> ('kind, 'structure) typ

val matkerint0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val qflll0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val qflllgram0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gtomap : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mapdelete : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val mapdomain : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mapdomain_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mapget :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mapisdefined :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val mapput :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val maptomat : ('kind, 'structure) typ -> ('kind, 'structure) typ
val maptomat_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val matpermanent : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zm_permanent : ('kind, 'structure) typ -> ('kind, 'structure) typ
val dbllemma526 : float -> float -> float -> float -> float
val dblcoro526 : float -> float -> float -> float

val gammamellininv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val gammamellininvasymp :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val gammamellininvinit :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val gammamellininvrt :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val member_a1 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_a2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_a3 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_a4 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_a6 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_area : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_b2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_b4 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_b6 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_b8 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_bid : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_bnf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_c4 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_c6 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_clgp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_codiff : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_cyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_diff : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_disc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_e : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_eta : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_f : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_fu : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_gen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_group : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_index : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_j : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_nf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_no : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_omega : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_orders : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_p : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_pol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_polabs : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_reg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_r1 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_r2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_roots : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_sign : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_t2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_tate : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_tu : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_zk : ('kind, 'structure) typ -> ('kind, 'structure) typ
val member_zkst : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_m : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_mindex : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_minv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_basis : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_dim : ('kind, 'structure) typ -> Signed.long
val mf_get_e : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_fields : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_newforms : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_space : ('kind, 'structure) typ -> Signed.long
val mf_get_s : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mfcusp_get_vmjd : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mfnew_get_vj : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qab_tracerel :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qabm_tracerel :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qabv_tracerel :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qab_trace_init :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val checkmf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val checkmf_i : ('kind, 'structure) typ -> int
val getcache : unit -> ('kind, 'structure) typ
val hclassno6u : pari_ulong -> pari_ulong
val hclassno6u_no_cache : pari_ulong -> pari_ulong

val lfunmf :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfdelta : unit -> ('kind, 'structure) typ
val mfeh : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mfek : Signed.long -> ('kind, 'structure) typ
val mftheta : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_chi : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_n : ('kind, 'structure) typ -> Signed.long
val mf_get_nk : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_field : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_gn : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_gk : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mf_get_k : ('kind, 'structure) typ -> Signed.long
val mf_get_r : ('kind, 'structure) typ -> Signed.long
val mf_get_type : ('kind, 'structure) typ -> Signed.long

val mfatkin :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfatkineigenvalues :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mfatkininit :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mfbasis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val mfbd : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mfbracket :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfcharorder : ('kind, 'structure) typ -> Signed.long
val mfcharmodulus : ('kind, 'structure) typ -> Signed.long
val mfcharpol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mfcoef : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mfcoefs :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mfconductor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val mfcosets : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfcuspdim :
  Signed.long -> Signed.long -> ('kind, 'structure) typ -> Signed.long

val mfcuspisregular :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val mfcusps : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfcuspval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfcuspwidth :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val mfderiv : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mfderive2 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mfdescribe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val mfdim : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mfdiv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfdiv_val :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfeigenbasis : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfeigensearch :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfeisenstein :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mfeisensteindim :
  Signed.long -> Signed.long -> ('kind, 'structure) typ -> Signed.long

val mfembed :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfembed0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfeval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mffields : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mffromell : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mffrometaquo :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mffromlfun :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mffromqf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mffulldim :
  Signed.long -> Signed.long -> ('kind, 'structure) typ -> Signed.long

val mfgaloisprojrep :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfgaloistype :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfhecke :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfheckemat :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfinit : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val mfiscm : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfiscuspidal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val mfisequal :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long

val mfisetaquo :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mfkohnenbasis : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mfkohnenbijection : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfkohneneigenbasis :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mflinear :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfmanin : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mfmatembed :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfmul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfnewdim :
  Signed.long -> Signed.long -> ('kind, 'structure) typ -> Signed.long

val mfolddim :
  Signed.long -> Signed.long -> ('kind, 'structure) typ -> Signed.long

val mfparams : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfperiodpol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mfperiodpolbasis : Signed.long -> Signed.long -> ('kind, 'structure) typ

val mfpetersson :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfpow : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mfsearch :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfshift : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mfshimura :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfslashexpansion :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ

val mfspace : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val mfsplit :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mfsturm : ('kind, 'structure) typ -> Signed.long
val mfsturmngk : Signed.long -> ('kind, 'structure) typ -> Signed.long
val mfsturmnk : Signed.long -> Signed.long -> Signed.long
val mfsturm_mf : ('kind, 'structure) typ -> Signed.long

val mfsymboleval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mfsymbol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mftaylor :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mftobasis :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mftobasises :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mftocol :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mftocoset :
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mftonew :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mftraceform :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mftwist :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfvecembed :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mfvectomat :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val fl_inv : pari_ulong -> pari_ulong -> pari_ulong
val fl_invsafe : pari_ulong -> pari_ulong -> pari_ulong

val fp_ratlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val zm2_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val abscmpii : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val abscmprr : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val absequalii : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val addii_sign :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val addir_sign :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val addmulii :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val addmulii_inplace :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val addrr_sign :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val addsi_sign :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val addsr : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val addui_sign :
  pari_ulong ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val addumului :
  pari_ulong -> pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val affir : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val affrr : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val cbezout :
  Signed.long ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val cgcd : Signed.long -> Signed.long -> Signed.long
val clcm : Signed.long -> Signed.long -> Signed.long
val cmpii : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val cmprr : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val dblexpo : float -> Signed.long
val dblmantissa : float -> pari_ulong
val dbltor : float -> ('kind, 'structure) typ

val divir :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val divis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val divis_rem :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val absdiviu_rem :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  ('kind, 'structure) typ

val diviuuexact :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val diviuexact :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val divri :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val divrr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val divrs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val divru : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val divsi : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val divsr : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val divur : pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val dvmdii :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val equalrr : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val floorr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val halfgcdii :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val int2n : Signed.long -> ('kind, 'structure) typ
val int2u : pari_ulong -> ('kind, 'structure) typ
val int2um1 : pari_ulong -> ('kind, 'structure) typ

val int_normalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val invmod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val invmod2bil : pari_ulong -> pari_ulong
val invr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mantissa_real :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val modiiz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val mulir :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mulrr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mulsi : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val mulsr : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val mulss : Signed.long -> Signed.long -> ('kind, 'structure) typ
val mului : pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val mulur : pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val muluu : pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val muluui :
  pari_ulong -> pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val pari_kernel_close : unit -> unit
val pari_kernel_init : unit -> unit
val pari_kernel_version : unit -> string
val remi2n : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val rtodbl : ('kind, 'structure) typ -> float
val shifti : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val sqri : ('kind, 'structure) typ -> ('kind, 'structure) typ
val sqrr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val sqrs : Signed.long -> ('kind, 'structure) typ
val sqrtr_abs : ('kind, 'structure) typ -> ('kind, 'structure) typ

val sqrtremi :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val sqru : pari_ulong -> ('kind, 'structure) typ
val subsr : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val truedvmdii :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val truedvmdis :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val truedvmdsi :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val trunc2nr : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mantissa2nr :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val truncr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ugcd : pari_ulong -> pari_ulong -> pari_ulong
val ulcm : pari_ulong -> pari_ulong -> pari_ulong
val umodiu : ('kind, 'structure) typ -> pari_ulong -> pari_ulong
val vals : pari_ulong -> Signed.long

val fpc_ratlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpm_ratlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_ratlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxqx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zxqx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nffactor :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nffactormod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfgcd_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val nfissquarefree : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val nfroots :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfroots_if_split :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfrootsof1 : ('kind, 'structure) typ -> ('kind, 'structure) typ

val polfnf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfabelianconjgen :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfisabelian :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val forpart :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val forpart_init :
  forpart_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val forpart_next :
  forpart_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val forpart_prev :
  forpart_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val numbpart : ('kind, 'structure) typ -> ('kind, 'structure) typ

val partitions :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val forperm :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  unit

val forperm_init :
  forperm_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  unit

val forperm_next :
  forperm_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val forallsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> Signed.long -> unit

val forksubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  unit

val forsubset_next :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val forsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  unit

val glambertw :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mplambertw :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mplambertx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mplambertx_logx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mplambertxlogx_x :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val z_to_perm :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val abelian_group : ('kind, 'structure) typ -> ('kind, 'structure) typ

val conjclasses_repr :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val cyc_pow : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val cyc_pow_perm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val cyclicgroup :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val dicyclicgroup :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val group_abelianhnf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val group_abeliansnf :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val group_domain : ('kind, 'structure) typ -> Signed.long

val group_elts :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val group_export :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val group_export_gap : ('kind, 'structure) typ -> ('kind, 'structure) typ
val group_export_magma : ('kind, 'structure) typ -> ('kind, 'structure) typ
val group_isa4s4 : ('kind, 'structure) typ -> Signed.long
val group_isabelian : ('kind, 'structure) typ -> Signed.long

val group_leftcoset :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val group_order : ('kind, 'structure) typ -> Signed.long

val group_perm_normalize :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val group_quotient :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val group_rightcoset :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val group_set :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val group_subgroup_is_faithful :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val group_subgroup_isnormal :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val group_subgroups : ('kind, 'structure) typ -> ('kind, 'structure) typ

val groupelts_solvablesubgroups :
  ('kind, 'structure) typ -> ('kind, 'structure) typ

val group_to_cc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val groupelts_abelian_group : ('kind, 'structure) typ -> ('kind, 'structure) typ
val groupelts_center : ('kind, 'structure) typ -> ('kind, 'structure) typ

val groupelts_conj_set :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val groupelts_conjclasses :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val groupelts_exponent : ('kind, 'structure) typ -> Signed.long

val groupelts_quotient :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val groupelts_set :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val groupelts_to_group : ('kind, 'structure) typ -> ('kind, 'structure) typ

val numtoperm :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val perm_commute : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val perm_cycles : ('kind, 'structure) typ -> ('kind, 'structure) typ
val perm_order : ('kind, 'structure) typ -> ('kind, 'structure) typ
val perm_orderu : ('kind, 'structure) typ -> pari_ulong

val perm_pow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val perm_powu : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val perm_sign : ('kind, 'structure) typ -> Signed.long
val perm_to_gap : ('kind, 'structure) typ -> ('kind, 'structure) typ
val perm_to_z : ('kind, 'structure) typ -> ('kind, 'structure) typ
val permcycles : ('kind, 'structure) typ -> ('kind, 'structure) typ
val permorder : ('kind, 'structure) typ -> ('kind, 'structure) typ
val permsign : ('kind, 'structure) typ -> Signed.long
val permtonum : ('kind, 'structure) typ -> ('kind, 'structure) typ

val quotient_group :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val quotient_groupelts : ('kind, 'structure) typ -> ('kind, 'structure) typ

val quotient_perm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val quotient_subgroup_lift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val subgroups_tableset :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val tableset_find_index :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val trivialgroup : unit -> ('kind, 'structure) typ

val vec_insert :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val vec_is1to1 : ('kind, 'structure) typ -> int
val vec_isconst : ('kind, 'structure) typ -> int

val vecperm_orbits :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vecsmall_duplicate : ('kind, 'structure) typ -> Signed.long
val vecsmall_duplicate_sorted : ('kind, 'structure) typ -> Signed.long
val vecsmall_indexsort : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecsmall_is1to1 : ('kind, 'structure) typ -> int
val vecsmall_isconst : ('kind, 'structure) typ -> int
val vecsmall_sort : ('kind, 'structure) typ -> unit
val vecsmall_uniq : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecsmall_uniq_sorted : ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecsmall_counting_indexsort :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vecsmall_counting_sort : ('kind, 'structure) typ -> Signed.long -> unit

val vecsmall_counting_uniq :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vecvecsmall_indexsort : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecvecsmall_max : ('kind, 'structure) typ -> Signed.long

val vecvecsmall_search :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val vecvecsmall_sort : ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecvecsmall_sort_inplace :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> unit

val vecvecsmall_sort_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecvecsmall_sort_uniq : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mt_broadcast : ('kind, 'structure) typ -> unit
val mt_nbthreads : unit -> Signed.long
val mt_queue_end : pari_mt Ctypes.structure Ctypes_static.ptr -> unit

val mt_queue_get :
  pari_mt Ctypes.structure Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val mt_queue_start :
  pari_mt Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ -> unit

val mt_queue_start_lim :
  pari_mt Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val mt_queue_submit :
  pari_mt Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ ->
  unit

val mt_sigint_block : unit -> unit
val mt_sigint_unblock : unit -> unit
val pari_mt_init : unit -> unit
val pari_mt_close : unit -> unit

val subcyclopclgp :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val subcycloiwasawa :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val subcyclohminus :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val color_to_rgb :
  ('kind, 'structure) typ ->
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
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val parplothexport :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val plotbox :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val plotclip : Signed.long -> unit

val plotcolor :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val plotcopy :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val plotcursor : Signed.long -> ('kind, 'structure) typ
val plotdraw : ('kind, 'structure) typ -> Signed.long -> unit

val plotexport :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val plothexport :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val plothraw :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val plothrawexport :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val plothsizes : Signed.long -> ('kind, 'structure) typ

val plotinit :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val plotkill : Signed.long -> unit

val plotline :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val plotlines :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val plotlinetype : Signed.long -> Signed.long -> unit

val plotmove :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val plotpoints :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val plotpointsize : Signed.long -> ('kind, 'structure) typ -> unit
val plotpointtype : Signed.long -> Signed.long -> unit

val plotrbox :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val plotrecth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val plotrecthraw :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val plotrline :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val plotrmove :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val plotrpoint :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val plotscale :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val plotstring : Signed.long -> string -> Signed.long -> unit
val psdraw : ('kind, 'structure) typ -> Signed.long -> unit

val psploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val psplothraw :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rect2ps :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  string

val rect2ps_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  int ->
  string

val rect2svg :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  string

val pariplot :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  unit

val zx_zp_root :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zp_appr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val cmp_padic : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val factorpadic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gdeuc :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val grem :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val padicappr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val poldivrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val polrootspadic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val flv_factorback :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> pari_ulong -> pari_ulong

val flxqv_factorback :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val fpv_factorback :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqv_factorback :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val q_content : ('kind, 'structure) typ -> ('kind, 'structure) typ
val q_content_safe : ('kind, 'structure) typ -> ('kind, 'structure) typ
val q_denom : ('kind, 'structure) typ -> ('kind, 'structure) typ
val q_denom_safe : ('kind, 'structure) typ -> ('kind, 'structure) typ

val q_div_to_int :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val q_gcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val q_mul_to_int :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val q_muli_to_int :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val q_primitive_part :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val q_primpart : ('kind, 'structure) typ -> ('kind, 'structure) typ

val q_remove_denom :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val q_factor : ('kind, 'structure) typ -> ('kind, 'structure) typ

val q_factor_limit :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val rg_type :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rgc_type :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rescale_to_int : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_type :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_type2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_rg_type :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_chinese_coprime :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgx_disc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_extgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgx_extgcd_simple :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgx_gcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_gcd_simple :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_halfgcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_halfgcd_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgx_rescale_to_int : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_resultant_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val rgx_sturmpart :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val rgx_sylvestermatrix :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_type :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type3 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type_decode :
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val rgx_type_is_composite : Signed.long -> int

val rgxq_charpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxq_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxq_minpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rgxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxq_ratlift :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val rgxq_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val z_content : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zx_content : ('kind, 'structure) typ -> ('kind, 'structure) typ

val centermod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val centermod_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val centermodii :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val content : ('kind, 'structure) typ -> ('kind, 'structure) typ

val content0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val deg1_from_roots :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val factor0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val factorback : ('kind, 'structure) typ -> ('kind, 'structure) typ

val factorback2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gbezout :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val gdivexact :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gen_factorback :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val ggcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ggcd0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ghalfgcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ginvmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val glcm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val glcm0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val newtonpoly :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfrootsq : ('kind, 'structure) typ -> ('kind, 'structure) typ
val poldisc0 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val polresultant0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val polsym : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val primitive_part :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val primpart : ('kind, 'structure) typ -> ('kind, 'structure) typ
val reduceddiscsmith : ('kind, 'structure) typ -> ('kind, 'structure) typ

val resultant2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val resultant :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rnfcharpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val roots_from_deg1 : ('kind, 'structure) typ -> ('kind, 'structure) typ

val roots_to_pol :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val roots_to_pol_r1 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val sturmpart :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long

val subresext :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val sylvestermatrix :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val trivial_fact : unit -> ('kind, 'structure) typ

val gcdext0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val polresultantext0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val polresultantext :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val prime_fact : ('kind, 'structure) typ -> ('kind, 'structure) typ
val row_q_primpart : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vec_q_primpart : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecprod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_lcm : ('kind, 'structure) typ -> ('kind, 'structure) typ

val flx_flxy_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxx_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val fpx_fpxy_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_translate :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxqx_normalize :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxv_fpc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxy_fpxq_evaly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fpxc_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxm_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_fp_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_halve :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_invsafe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_mulu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_neg :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_neg_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_sqrt :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_sqrtn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fq_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqc_fq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqc_fqv_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqc_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqc_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqv_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqv_roots_to_pol :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fqx_fq_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_fq_mul_to_monic :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_fq_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_translate :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_matrix_pow :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_powers :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxy_eval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxy_evalx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qx_disc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qx_gcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qx_resultant :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxq_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxq_intnorm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxq_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qxq_norm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qxq_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rg_is_fp :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> int

val rg_is_fpxq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val rg_to_fp :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rg_to_fpxq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgc_to_fpc :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgc_to_fqc :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgm_is_fpm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> int

val rgm_to_flm :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val rgm_to_fpm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgm_to_fqm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgv_is_fpv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> int

val rgv_to_flv :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val rgv_to_fpv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_is_fpx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> int

val rgx_to_fpx :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_is_fpxqx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val rgx_to_fpxqx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgx_to_fqx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val z_incremental_crt :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong ->
  int

val z_init_crt : pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val zm_incremental_crt :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong ->
  int

val zm_init_crt :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val zx_zxy_resultant :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_zxy_rnfequation :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zx_disc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_gcd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_gcd_all :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zx_incremental_crt :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong ->
  int

val zx_init_crt :
  ('kind, 'structure) typ ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ

val zx_is_squarefree : ('kind, 'structure) typ -> int
val zx_radical : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_resultant :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zxm_incremental_crt :
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  pari_ulong ->
  int

val zxm_init_crt :
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val zxq_minpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) typ

val zxq_charpoly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val characteristic : ('kind, 'structure) typ -> ('kind, 'structure) typ

val ffnbirred :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ffnbirred0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val ffsumnbirred :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val get_fq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  bb_field Ctypes.structure Ctypes_static.ptr

val init_flxq :
  pari_ulong -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val init_fq :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val nfx_disc :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfx_resultant :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val pol_x_powers : Signed.long -> Signed.long -> ('kind, 'structure) typ
val residual_characteristic : ('kind, 'structure) typ -> ('kind, 'structure) typ

val polclass :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val fp_modinv_to_j :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_polmodular_evalx :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  int ->
  ('kind, 'structure) typ

val check_modinv : Signed.long -> unit
val disc_best_modinv : Signed.long -> Signed.long
val modinv_height_factor : Signed.long -> Signed.long
val modinv_good_disc : Signed.long -> Signed.long -> int
val modinv_good_prime : Signed.long -> Signed.long -> int
val modinv_is_weber : Signed.long -> int
val modinv_is_double_eta : Signed.long -> int

val polmodular :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val polmodular_zm : Signed.long -> Signed.long -> ('kind, 'structure) typ

val polmodular_zxx :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val bpsw_isprime : ('kind, 'structure) typ -> Signed.long
val bpsw_psp : ('kind, 'structure) typ -> Signed.long
val addprimes : ('kind, 'structure) typ -> ('kind, 'structure) typ
val check_ecppcert : ('kind, 'structure) typ -> Signed.long
val gisprime : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gispseudoprime :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gprimepi_upper_bound : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gprimepi_lower_bound : ('kind, 'structure) typ -> ('kind, 'structure) typ
val isprime : ('kind, 'structure) typ -> Signed.long
val ispseudoprime : ('kind, 'structure) typ -> Signed.long -> Signed.long
val millerrabin : ('kind, 'structure) typ -> Signed.long -> Signed.long
val prime : Signed.long -> ('kind, 'structure) typ

val primecert :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val primecert0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val primecertexport :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val primecertisvalid : ('kind, 'structure) typ -> Signed.long
val primepi : ('kind, 'structure) typ -> ('kind, 'structure) typ
val primepi_upper_bound : float -> float
val primepi_lower_bound : float -> float
val primes : Signed.long -> ('kind, 'structure) typ

val primes_interval :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val primes_interval_zv : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val primes_upto_zv : pari_ulong -> ('kind, 'structure) typ
val primes0 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val primes_zv : Signed.long -> ('kind, 'structure) typ

val randomprime0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val removeprimes : ('kind, 'structure) typ -> ('kind, 'structure) typ
val uis2psp : pari_ulong -> int
val uispsp : pari_ulong -> pari_ulong -> int
val uislucaspsp : pari_ulong -> int
val uisprime : pari_ulong -> int
val uisprime_101 : pari_ulong -> int
val uisprime_661 : pari_ulong -> int
val uprime : Signed.long -> pari_ulong
val uprimepi : pari_ulong -> pari_ulong

val qfauto :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfauto0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfautoexport :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val qfisom :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfisom0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfisominit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfisominit0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qforbits :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfminimize : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfparam :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qfsolve : ('kind, 'structure) typ -> ('kind, 'structure) typ
val z_isfundamental : ('kind, 'structure) typ -> Signed.long
val classno : ('kind, 'structure) typ -> ('kind, 'structure) typ
val classno2 : ('kind, 'structure) typ -> ('kind, 'structure) typ

val hclassnof_fact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val hclassno : ('kind, 'structure) typ -> ('kind, 'structure) typ
val hclassno6 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val isfundamental : ('kind, 'structure) typ -> Signed.long
val qfb_equal1 : ('kind, 'structure) typ -> int

val qfbclassno0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val qfi_shanks :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qfi_log :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfi_order :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val quadclassnof :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val quadclassnof_fact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val quaddisc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val quadregulator :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val quadunit : ('kind, 'structure) typ -> ('kind, 'structure) typ

val quadunit0 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val quadunitindex :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val quadunitnorm : ('kind, 'structure) typ -> Signed.long
val sisfundamental : Signed.long -> Signed.long
val uhclassnof_fact : ('kind, 'structure) typ -> Signed.long -> Signed.long
val unegisfundamental : pari_ulong -> Signed.long
val unegquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val uposisfundamental : pari_ulong -> Signed.long
val uposquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val uquadclassnof_fact :
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong

val zn_quad_roots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val getrand : unit -> ('kind, 'structure) typ
val pari_rand : unit -> pari_ulong
val randomr : Signed.long -> ('kind, 'structure) typ
val random_f2x : Signed.long -> Signed.long -> ('kind, 'structure) typ
val random_fl : pari_ulong -> pari_ulong
val random_bits : Signed.long -> Signed.long
val random_zv : Signed.long -> ('kind, 'structure) typ
val setrand : ('kind, 'structure) typ -> unit

val ellratpoints :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val hyperellratpoints :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qx_complex_roots :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val fft :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fftinv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val cleanroots :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val fujiwara_bound : ('kind, 'structure) typ -> float
val fujiwara_bound_real : ('kind, 'structure) typ -> Signed.long -> float
val isrealappr : ('kind, 'structure) typ -> Signed.long -> int
val polgraeffe : ('kind, 'structure) typ -> ('kind, 'structure) typ

val polmod_to_embed :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val polrootsbound :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val roots : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val realroots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zx_graeffe : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zx_realroots_irred :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zx_sturm : ('kind, 'structure) typ -> Signed.long
val zx_sturm_irred : ('kind, 'structure) typ -> Signed.long

val zx_sturmpart :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val zx_uspensky :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val factor_aurifeuille :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val factor_aurifeuille_prime :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val galoissubcyclo :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val polsubcyclo :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val polsubcyclofast :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val znsubgroupgenerators :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val nfsubfields :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val nfsubfields0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val nfsubfieldscm :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val nfsubfieldsmax :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val nflist :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val nfresolvent :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val subgrouplist :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val forsubgroup :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val abmap_kernel : ('kind, 'structure) typ -> ('kind, 'structure) typ

val abmap_subgroup_image :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnrl1 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrrootnumber :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val bnrstark :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val cyc2elts : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qfbforms : ('kind, 'structure) typ -> ('kind, 'structure) typ

val quadhilbert :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val quadray :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val chartogenstr : char -> ('kind, 'structure) typ
val pari_strdup : string -> string
val pari_strndup : string -> Signed.long -> string
val stack_strcat : string -> string -> string
val stack_strdup : string -> string
val pari_strchr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val strjoin :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val strntogenstr : string -> Signed.long -> ('kind, 'structure) typ

val strsplit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val strtogenstr : string -> ('kind, 'structure) typ
val type_name : Signed.long -> string

val asympnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val asympnumraw :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val derivnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val derivnumk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val derivfun :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val derivfunk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val forvec_init :
  forvec_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  int

val forvec_next :
  forvec_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val laurentseries :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val limitnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val polzag : Signed.long -> Signed.long -> ('kind, 'structure) typ

val prodeuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val prodinf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val prodinf1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val solvestep :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val sumalt :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumalt2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumpos :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumpos2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val suminf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val suminf_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sumdivmultexpr :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zbrent :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val bnfisintnorm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnfisintnormabs :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val ideals_by_norm :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val thue :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val thueinit :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val pi2n : Signed.long -> Signed.long -> ('kind, 'structure) typ
val pii2 : Signed.long -> ('kind, 'structure) typ
val pii2n : Signed.long -> Signed.long -> ('kind, 'structure) typ
val qp_exp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qp_exp_prec : ('kind, 'structure) typ -> Signed.long
val qp_log : ('kind, 'structure) typ -> ('kind, 'structure) typ
val qp_sqrt : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qp_sqrtn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val zn_sqrt :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zp_teichmuller :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val agm :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val constcatalan : Signed.long -> ('kind, 'structure) typ
val consteuler : Signed.long -> ('kind, 'structure) typ
val constlog2 : Signed.long -> ('kind, 'structure) typ
val constpi : Signed.long -> ('kind, 'structure) typ
val cxexpm1 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val elle : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ellk : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val expir : ('kind, 'structure) typ -> ('kind, 'structure) typ
val exp1r_abs : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gcos : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gcotan : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gcotanh : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gexp : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gexpm1 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val glog : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val glog1p : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gpow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gpowers : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gpowers0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gpowgs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val grootsof1 : Signed.long -> Signed.long -> ('kind, 'structure) typ
val gsin : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gsinc : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gsincos :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  unit

val gsqrpowers :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gsqrt : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gsqrtn :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) typ

val gtan : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val logr_abs : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mpcos : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mpeuler : Signed.long -> ('kind, 'structure) typ
val mpcatalan : Signed.long -> ('kind, 'structure) typ

val mpsincosm1 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val mpexp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mpexpm1 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mplog : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mplog2 : Signed.long -> ('kind, 'structure) typ
val mppi : Signed.long -> ('kind, 'structure) typ
val mpsin : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpsincos :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val pow2pis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val powpis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val powcx :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val powcx_prec :
  Signed.long -> ('kind, 'structure) typ -> Signed.long -> Signed.long

val powersr : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val powiu : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val powrfrac :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val powrs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val powrshalf :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val powru : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val powruhalf : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val powgi :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rootsof1_cx :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rootsof1u_cx : pari_ulong -> Signed.long -> ('kind, 'structure) typ

val rootsof1q_cx :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val rootsof1powinit :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val rootsof1pow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val serchop : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val serchop_i :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val serchop0 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val sqrtnint : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val sqrtnr_abs :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val teich : ('kind, 'structure) typ -> ('kind, 'structure) typ
val teichmullerinit : Signed.long -> Signed.long -> ('kind, 'structure) typ

val teichmuller :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val trans_eval :
  string ->
  (('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val trans_evalgen :
  string ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val upowuu : pari_ulong -> pari_ulong -> pari_ulong
val upowers : pari_ulong -> Signed.long -> ('kind, 'structure) typ
val usqrtn : pari_ulong -> pari_ulong -> pari_ulong
val usqrt : pari_ulong -> pari_ulong
val qp_gamma : ('kind, 'structure) typ -> ('kind, 'structure) typ
val atanhuu : pari_ulong -> pari_ulong -> Signed.long -> ('kind, 'structure) typ

val atanhui :
  pari_ulong ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gacosh : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gacos : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val garg : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gasinh : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gasin : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gatan : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gatanh : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gcosh : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ggammah : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val ggamma : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ggamma1m1 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val glngamma : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gpsi : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gsinh : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gtanh : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val mpfactr : Signed.long -> Signed.long -> ('kind, 'structure) typ

val mpsinhcosh :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val psi1series :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val sumformal :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgv_is_arithprog :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  int

val besseljzero :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val besselyzero :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val constzeta : Signed.long -> Signed.long -> ('kind, 'structure) typ

val cxek :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val dblmodulus : ('kind, 'structure) typ -> float
val dilog : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val eint1 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val expipir : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val expipic : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val expixy :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val eta : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val eta0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val gerfc : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val gpolylog :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val gzeta : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val hbessel1 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val hbessel2 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val hyperu :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val ibessel :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val incgam :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val incgam0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val incgamc :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val jbessel :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val jbesselh :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val jell : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val kbessel :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mpeint1 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpveceint1 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val polylog0 :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val sumdedekind :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val sumdedekind_coprime :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val szeta : Signed.long -> Signed.long -> ('kind, 'structure) typ

val theta :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val thetanullk :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val trueeta : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val u_sumdedekind_coprime :
  Signed.long -> Signed.long -> ('kind, 'structure) typ

val upper_to_cx :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val veceint1 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val vecthetanullk :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val vecthetanullk_tau :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val veczeta :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val weber0 :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val weberf : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val weberf1 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val weberf2 : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val ybessel :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val sl2_inv_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qevproj_apply :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qevproj_apply_vecei :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val qevproj_down :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val qevproj_init : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_act_gl2q :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_act_zgl2q :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val checkms : ('kind, 'structure) typ -> unit
val checkmspadic : ('kind, 'structure) typ -> unit

val ellpadiclambdamu :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mfnumcusps : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mfnumcusps_fact : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mfnumcuspsu_fact : ('kind, 'structure) typ -> pari_ulong
val mfnumcuspsu : pari_ulong -> pari_ulong

val msfromcusp :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val msfromell :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val msfromhecke :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val msdim : ('kind, 'structure) typ -> Signed.long

val mseval2_ooq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val msgetlevel : ('kind, 'structure) typ -> Signed.long
val msgetsign : ('kind, 'structure) typ -> Signed.long
val msgetweight : ('kind, 'structure) typ -> Signed.long

val msatkinlehner :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mscuspidal :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mseisenstein : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mseval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mshecke :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val msinit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val msissymbol :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mslattice :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val msomseval :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mspadic_parse_chi :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  unit

val mspadic_unit_eigenvalue :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mspadicinit :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mspadicl :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mspadicmoments :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mspadicseries :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val mspathgens : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mspathlog :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val msnew : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mspetersson :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mspolygon :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val msstar :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val msqexpansion :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val mssplit :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val mstooms :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mscosets0 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mscosets :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ

val msfarey :
  ('kind, 'structure) typ ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) typ -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val msfarey0 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val checkfarey_i : ('kind, 'structure) typ -> int

val polylogmult :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val polylogmult_interpolate :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zetamult : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val zetamultdual : ('kind, 'structure) typ -> ('kind, 'structure) typ

val zetamult_interpolate :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val zetamultall :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val zetamultconvert :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

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
val abscmpiu : ('kind, 'structure) typ -> pari_ulong -> int
val abscmpui : pari_ulong -> ('kind, 'structure) typ -> int
val absequaliu : ('kind, 'structure) typ -> pari_ulong -> int
val absi : ('kind, 'structure) typ -> ('kind, 'structure) typ
val absi_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val absr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val absrnz_equal1 : ('kind, 'structure) typ -> int
val absrnz_equal2n : ('kind, 'structure) typ -> int

val addiiz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val addir :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val addirz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val addis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val addri :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val addriz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val addrr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val addrrz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val addrs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val addsi : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val addsiz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val addsrz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val addss : Signed.long -> Signed.long -> ('kind, 'structure) typ
val addssz : Signed.long -> Signed.long -> ('kind, 'structure) typ -> unit
val adduu : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val affii : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val affiz : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val affrr_fixlg : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val affsi : Signed.long -> ('kind, 'structure) typ -> unit
val affsr : Signed.long -> ('kind, 'structure) typ -> unit
val affsz : Signed.long -> ('kind, 'structure) typ -> unit
val affui : pari_ulong -> ('kind, 'structure) typ -> unit
val affur : pari_ulong -> ('kind, 'structure) typ -> unit
val cgetg_block : Signed.long -> Signed.long -> ('kind, 'structure) typ

val cgetg_copy :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val cgeti : Signed.long -> ('kind, 'structure) typ
val cgetineg : Signed.long -> ('kind, 'structure) typ
val cgetipos : Signed.long -> ('kind, 'structure) typ
val cgetr : Signed.long -> ('kind, 'structure) typ
val cgetr_block : Signed.long -> ('kind, 'structure) typ
val cmpir : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val cmpis : ('kind, 'structure) typ -> Signed.long -> int
val cmpiu : ('kind, 'structure) typ -> pari_ulong -> int
val cmpri : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val cmprs : ('kind, 'structure) typ -> Signed.long -> int
val cmpsi : Signed.long -> ('kind, 'structure) typ -> int
val cmpsr : Signed.long -> ('kind, 'structure) typ -> int
val cmpss : Signed.long -> Signed.long -> int
val cmpui : pari_ulong -> ('kind, 'structure) typ -> int
val cmpuu : pari_ulong -> pari_ulong -> int

val divii :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val diviiz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val divirz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val divisz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> unit

val divriz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val divrrz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val divrsz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> unit

val divsi_rem :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val divsiz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val divsrz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val divss : Signed.long -> Signed.long -> ('kind, 'structure) typ

val divss_rem :
  Signed.long ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) typ

val divssz : Signed.long -> Signed.long -> ('kind, 'structure) typ -> unit
val dvdii : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val dvdiiz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  int

val dvdis : ('kind, 'structure) typ -> Signed.long -> int

val dvdisz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> int

val dvdiu : ('kind, 'structure) typ -> pari_ulong -> int

val dvdiuz :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ -> int

val dvdsi : Signed.long -> ('kind, 'structure) typ -> int
val dvdui : pari_ulong -> ('kind, 'structure) typ -> int

val dvmdiiz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val dvmdis :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val dvmdisz :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val dvmdsbil : Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val dvmdsi :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val dvmdsiz :
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val dvmdss :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val dvmdssz :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val dvmdubil : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val equalis : ('kind, 'structure) typ -> Signed.long -> int
val equalsi : Signed.long -> ('kind, 'structure) typ -> int
val equalui : pari_ulong -> ('kind, 'structure) typ -> int
val equaliu : ('kind, 'structure) typ -> pari_ulong -> int
val absequalui : pari_ulong -> ('kind, 'structure) typ -> int
val ceildivuu : pari_ulong -> pari_ulong -> pari_ulong
val evalexpo : Signed.long -> Signed.long
val evallg : Signed.long -> Signed.long
val evalprecp : Signed.long -> Signed.long
val evalvalp : Signed.long -> Signed.long
val evalvalser : Signed.long -> Signed.long
val expi : ('kind, 'structure) typ -> Signed.long
val expu : pari_ulong -> Signed.long
val fixlg : ('kind, 'structure) typ -> Signed.long -> unit
val fractor : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gc_bool : pari_ulong -> int -> int
val gc_const : pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val gc_double : pari_ulong -> float -> float
val gc_int : pari_ulong -> int -> int
val gc_long : pari_ulong -> Signed.long -> Signed.long
val gc_stoi : pari_ulong -> Signed.long -> ('kind, 'structure) typ
val gc_ulong : pari_ulong -> pari_ulong -> pari_ulong
val gc_utoi : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val gc_utoipos : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val gc_null : pari_ulong -> ('kind, 'structure) typ
val icopy : ('kind, 'structure) typ -> ('kind, 'structure) typ

val icopyspec :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val int_bit : ('kind, 'structure) typ -> Signed.long -> pari_ulong
val itor : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val itos : ('kind, 'structure) typ -> Signed.long
val itos_or_0 : ('kind, 'structure) typ -> Signed.long
val itou : ('kind, 'structure) typ -> pari_ulong
val itou_or_0 : ('kind, 'structure) typ -> pari_ulong
val leafcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val maxdd : float -> float -> float
val maxss : Signed.long -> Signed.long -> Signed.long
val maxuu : pari_ulong -> pari_ulong -> Signed.long
val mindd : float -> float -> float
val minss : Signed.long -> Signed.long -> Signed.long
val minuu : pari_ulong -> pari_ulong -> Signed.long
val mod16 : ('kind, 'structure) typ -> Signed.long
val mod2 : ('kind, 'structure) typ -> Signed.long
val mod2bil : ('kind, 'structure) typ -> pari_ulong
val mod32 : ('kind, 'structure) typ -> Signed.long
val mod4 : ('kind, 'structure) typ -> Signed.long
val mod64 : ('kind, 'structure) typ -> Signed.long
val mod8 : ('kind, 'structure) typ -> Signed.long
val modis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val modisz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> unit

val modsi : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val modsiz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val modss : Signed.long -> Signed.long -> ('kind, 'structure) typ
val modssz : Signed.long -> Signed.long -> ('kind, 'structure) typ -> unit
val mpabs : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mpabs_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpadd :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpaddz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val mpaff : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val mpceil : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mpcmp : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val mpcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpdiv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpexpo : ('kind, 'structure) typ -> Signed.long
val mpfloor : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpmul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpmulz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val mpneg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mpodd : ('kind, 'structure) typ -> int
val mpround : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mpsqr : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpsub :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mpsubz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val mptrunc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val muliiz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val mulirz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val mulis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val muliu : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val mulri :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mulriz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val mulrrz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val mulrs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val mulru : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val mulsiz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val mulsrz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val mulssz : Signed.long -> Signed.long -> ('kind, 'structure) typ -> unit
val negr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val new_chunk : int -> ('kind, 'structure) typ
val rcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rdivii :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rdiviiz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val rdivis :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rdivsi :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val rdivss :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val real2n : Signed.long -> Signed.long -> ('kind, 'structure) typ
val real_m2n : Signed.long -> Signed.long -> ('kind, 'structure) typ
val real_0 : Signed.long -> ('kind, 'structure) typ
val real_0_bit : Signed.long -> ('kind, 'structure) typ
val real_1 : Signed.long -> ('kind, 'structure) typ
val real_1_bit : Signed.long -> ('kind, 'structure) typ
val real_m1 : Signed.long -> ('kind, 'structure) typ

val remii :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val remiiz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val remis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val remisz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> unit

val remlll_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val remsi : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val remsiz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val remss : Signed.long -> Signed.long -> ('kind, 'structure) typ
val remssz : Signed.long -> Signed.long -> ('kind, 'structure) typ -> unit
val rtor : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val sdivsi : Signed.long -> ('kind, 'structure) typ -> Signed.long

val sdivsi_rem :
  Signed.long ->
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val sdivss_rem :
  Signed.long -> Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val get_avma : unit -> pari_ulong
val set_avma : pari_ulong -> unit

val uabsdiviu_rem :
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong

val udivuu_rem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val umodi2n : ('kind, 'structure) typ -> Signed.long -> pari_ulong
val setabssign : ('kind, 'structure) typ -> unit

val shift_left :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  unit

val shift_right :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  unit

val shiftl : pari_ulong -> pari_ulong -> pari_ulong
val shiftlr : pari_ulong -> pari_ulong -> pari_ulong
val shiftr : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val shiftr_inplace : ('kind, 'structure) typ -> Signed.long -> unit
val smodis : ('kind, 'structure) typ -> Signed.long -> Signed.long
val smodss : Signed.long -> Signed.long -> Signed.long
val stackdummy : pari_ulong -> pari_ulong -> unit
val stack_malloc : int -> string
val stack_malloc_align : int -> Signed.long -> string
val stack_calloc : int -> string
val stack_calloc_align : int -> Signed.long -> string

val subiiz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val subir :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val subirz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val subis : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val subisz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> unit

val subri :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val subriz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val subrr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val subrrz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val subrs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val subrsz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> unit

val subsi : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val subsiz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val subsrz :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val subss : Signed.long -> Signed.long -> ('kind, 'structure) typ
val subssz : Signed.long -> Signed.long -> ('kind, 'structure) typ -> unit
val subuu : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val togglesign : ('kind, 'structure) typ -> unit
val togglesign_safe : ('kind, 'structure) typ Ctypes_static.ptr -> unit
val affectsign : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val affectsign_safe :
  ('kind, 'structure) typ -> ('kind, 'structure) typ Ctypes_static.ptr -> unit

val truedivii :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val truedivis :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val truedivsi :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val uabsdivui_rem :
  pari_ulong ->
  ('kind, 'structure) typ ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong

val umodsu : Signed.long -> pari_ulong -> pari_ulong
val umodui : pari_ulong -> ('kind, 'structure) typ -> pari_ulong
val ugcdiu : ('kind, 'structure) typ -> pari_ulong -> pari_ulong
val ugcdui : pari_ulong -> ('kind, 'structure) typ -> pari_ulong
val umuluu_le : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val umuluu_or_0 : pari_ulong -> pari_ulong -> pari_ulong
val utoi : pari_ulong -> ('kind, 'structure) typ
val utoineg : pari_ulong -> ('kind, 'structure) typ
val utoipos : pari_ulong -> ('kind, 'structure) typ
val utor : pari_ulong -> Signed.long -> ('kind, 'structure) typ
val uutoi : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val uutoineg : pari_ulong -> pari_ulong -> ('kind, 'structure) typ
val vali : ('kind, 'structure) typ -> Signed.long
val varncmp : Signed.long -> Signed.long -> int
val varnmax : Signed.long -> Signed.long -> Signed.long
val varnmin : Signed.long -> Signed.long -> Signed.long

val pari_err_component :
  string -> string -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val pari_err_dim : string -> unit

val pari_err_domain :
  string ->
  string ->
  string ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val pari_err_file : string -> string -> unit
val pari_err_filedesc : string -> Signed.long -> unit
val pari_err_flag : string -> unit
val pari_err_impl : string -> unit
val pari_err_inv : string -> ('kind, 'structure) typ -> unit
val pari_err_irredpol : string -> ('kind, 'structure) typ -> unit
val pari_err_maxprime : pari_ulong -> unit

val pari_err_modulus :
  string -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val pari_err_op :
  string -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val pari_err_overflow : string -> unit
val pari_err_package : string -> unit
val pari_err_prec : string -> unit
val pari_err_prime : string -> ('kind, 'structure) typ -> unit

val pari_err_priority :
  string -> ('kind, 'structure) typ -> string -> Signed.long -> unit

val pari_err_sqrtn : string -> ('kind, 'structure) typ -> unit
val pari_err_type : string -> ('kind, 'structure) typ -> unit

val pari_err_type2 :
  string -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val pari_err_var :
  string -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val pari_err_roots0 : string -> unit

val mkintmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkintmodu : pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val mkpolmod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkfrac :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkfracss : Signed.long -> Signed.long -> ('kind, 'structure) typ

val qtoss :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val sstoq : Signed.long -> Signed.long -> ('kind, 'structure) typ
val uutoq : pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val mkfraccopy :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkrfrac :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkrfraccopy :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gen_i : unit -> ('kind, 'structure) typ
val cgetc : Signed.long -> ('kind, 'structure) typ

val mkquad :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkvecsmall : Signed.long -> ('kind, 'structure) typ
val mkvecsmall2 : Signed.long -> Signed.long -> ('kind, 'structure) typ

val mkvecsmall3 :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val mkvecsmall4 :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mkvecsmall5 :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mkqfb :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkvec3 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkvec4 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkvec5 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkvecs : Signed.long -> ('kind, 'structure) typ
val mkvec2s : Signed.long -> Signed.long -> ('kind, 'structure) typ

val mkvec3s :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val mkvec4s :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mkveccopy : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkvec2copy :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkcol : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkcol2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkcol3 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkcol4 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkcol5 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkcol6 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkcols : Signed.long -> ('kind, 'structure) typ
val mkcol2s : Signed.long -> Signed.long -> ('kind, 'structure) typ

val mkcol3s :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) typ

val mkcol4s :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mkcolcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mkmat : ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkmat2 :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mkmat3 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkmat4 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkmat5 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mkmatcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val mkerr : Signed.long -> ('kind, 'structure) typ
val mkoo : unit -> ('kind, 'structure) typ
val mkmoo : unit -> ('kind, 'structure) typ
val inf_get_sign : ('kind, 'structure) typ -> Signed.long

val mkmat22s :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val mkmat22 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val pol_x : Signed.long -> ('kind, 'structure) typ
val pol_xn : Signed.long -> Signed.long -> ('kind, 'structure) typ
val pol_xnall : Signed.long -> Signed.long -> ('kind, 'structure) typ
val polxn_flx : Signed.long -> Signed.long -> ('kind, 'structure) typ
val pol_1 : Signed.long -> ('kind, 'structure) typ
val pol_0 : Signed.long -> ('kind, 'structure) typ

val const_vec :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val const_col :
  Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val const_vecsmall : Signed.long -> Signed.long -> ('kind, 'structure) typ

val zeropadic :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zeropadic_shallow :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zeroser : Signed.long -> Signed.long -> ('kind, 'structure) typ
val ser_isexactzero : ('kind, 'structure) typ -> int
val zeropol : Signed.long -> ('kind, 'structure) typ
val zerocol : Signed.long -> ('kind, 'structure) typ
val zerovec : Signed.long -> ('kind, 'structure) typ
val zeromat : Signed.long -> Signed.long -> ('kind, 'structure) typ
val zero_flx : Signed.long -> ('kind, 'structure) typ
val zero_flv : Signed.long -> ('kind, 'structure) typ
val zero_flm : Signed.long -> Signed.long -> ('kind, 'structure) typ
val zero_flm_copy : Signed.long -> Signed.long -> ('kind, 'structure) typ
val zero_f2v : Signed.long -> ('kind, 'structure) typ
val zero_f2m : Signed.long -> Signed.long -> ('kind, 'structure) typ
val zero_f2m_copy : Signed.long -> Signed.long -> ('kind, 'structure) typ
val zeromatcopy : Signed.long -> Signed.long -> ('kind, 'structure) typ
val col_ei : Signed.long -> Signed.long -> ('kind, 'structure) typ
val vec_ei : Signed.long -> Signed.long -> ('kind, 'structure) typ
val f2v_ei : Signed.long -> Signed.long -> ('kind, 'structure) typ
val vecsmall_ei : Signed.long -> Signed.long -> ('kind, 'structure) typ

val rg_col_ei :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val shallowcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vectrunc_init : Signed.long -> ('kind, 'structure) typ
val coltrunc_init : Signed.long -> ('kind, 'structure) typ
val lg_increase : ('kind, 'structure) typ -> unit
val vectrunc_append : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val vectrunc_append_batch :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val vecsmalltrunc_init : Signed.long -> ('kind, 'structure) typ
val vecsmalltrunc_append : ('kind, 'structure) typ -> Signed.long -> unit
val hash_str : string -> pari_ulong
val hash_str_len : string -> Signed.long -> pari_ulong

val vec_shorten :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vec_lengthen :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vec_append :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vec_prepend :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vec_setconst :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecsmall_shorten :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vecsmall_lengthen :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vec_to_vecsmall : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecsmall_to_vec : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecsmall_to_vec_inplace : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecsmall_to_col : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecsmall_lexcmp : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val vecsmall_prefixcmp :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val vecsmall_prepend :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vecsmall_append :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vecsmall_concat :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecsmall_coincidence :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val vecsmall_isin : ('kind, 'structure) typ -> Signed.long -> Signed.long

val vecsmall_pack :
  ('kind, 'structure) typ -> Signed.long -> Signed.long -> Signed.long

val vecsmall_indexmax : ('kind, 'structure) typ -> Signed.long
val vecsmall_max : ('kind, 'structure) typ -> Signed.long
val vecsmall_indexmin : ('kind, 'structure) typ -> Signed.long
val vecsmall_min : ('kind, 'structure) typ -> Signed.long
val zv_isscalar : ('kind, 'structure) typ -> int
val qv_isscalar : ('kind, 'structure) typ -> int
val rgv_isscalar : ('kind, 'structure) typ -> int
val rgx_isscalar : ('kind, 'structure) typ -> int

val rgx_equal_var :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val rgx_to_rgv :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_is_rational : ('kind, 'structure) typ -> int
val rgx_is_zx : ('kind, 'structure) typ -> int
val rgx_is_qx : ('kind, 'structure) typ -> int
val rgx_is_monomial : ('kind, 'structure) typ -> int
val rgv_is_zv : ('kind, 'structure) typ -> int
val rgv_is_qv : ('kind, 'structure) typ -> int

val rgv_isin_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long

val rgv_isin : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long

val vecslicepermute :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rowslicepermute :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rowslice :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val matslice :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val rowsplice :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val vecsplice :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_minor :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val row : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val flm_row : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val rowcopy : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val row_i :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val vecreverse : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecsmall_reverse : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecreverse_inplace : ('kind, 'structure) typ -> unit

val vecsmallpermute :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecpermute :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rowpermute :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val identity_zv : Signed.long -> ('kind, 'structure) typ
val identity_perm : Signed.long -> ('kind, 'structure) typ
val cyclic_perm : Signed.long -> Signed.long -> ('kind, 'structure) typ

val perm_mul :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val perm_sqr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val perm_inv : ('kind, 'structure) typ -> ('kind, 'structure) typ

val perm_conj :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val pari_free : unit Ctypes_static.ptr -> unit
val pari_malloc : int -> unit Ctypes_static.ptr
val pari_realloc : unit Ctypes_static.ptr -> int -> unit Ctypes_static.ptr
val pari_realloc_ip : unit Ctypes_static.ptr Ctypes_static.ptr -> int -> unit
val pari_calloc : int -> unit Ctypes_static.ptr
val cgetalloc : int -> Signed.long -> ('kind, 'structure) typ

val icopy_avma :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val leafcopy_avma :
  ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ

val gerepileuptoleaf :
  pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gerepileuptoint :
  pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gerepileupto :
  pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gerepilecopy :
  pari_ulong -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gunclonenull : ('kind, 'structure) typ -> unit
val gunclonenull_deep : ('kind, 'structure) typ -> unit

val gerepilemany :
  pari_ulong ->
  ('kind, 'structure) typ Ctypes_static.ptr Ctypes_static.ptr ->
  int ->
  unit

val gerepileall : pari_ulong -> int -> unit
val gc_all : pari_ulong -> int -> ('kind, 'structure) typ
val gerepilecoeffs : pari_ulong -> ('kind, 'structure) typ -> int -> unit

val bin_copy :
  genbin Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val genbinbase :
  genbin Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) typ

val cgiv : ('kind, 'structure) typ -> unit
val killblock : ('kind, 'structure) typ -> unit
val is_universal_constant : ('kind, 'structure) typ -> int

val cxcompotor :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val cxtofp : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val cxtoreal : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gtodouble : ('kind, 'structure) typ -> float
val gisdouble : ('kind, 'structure) typ -> float Ctypes_static.ptr -> int
val gtos : ('kind, 'structure) typ -> Signed.long
val gtou : ('kind, 'structure) typ -> pari_ulong
val absfrac : ('kind, 'structure) typ -> ('kind, 'structure) typ
val absfrac_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val q_abs : ('kind, 'structure) typ -> ('kind, 'structure) typ
val q_abs_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val r_abs_shallow : ('kind, 'structure) typ -> ('kind, 'structure) typ
val r_abs : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gtofp : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gtomp : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_gtofp :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgc_gtofp :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgv_gtofp :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_gtofp :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgc_gtomp :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_gtomp :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_fpnorml2 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgc_fpnorml2 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgm_fpnorml2 :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val affgr : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val affc_fixlg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val trunc_safe : ('kind, 'structure) typ -> ('kind, 'structure) typ
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
val bit_prec : ('kind, 'structure) typ -> Signed.long
val bit_accuracy : Signed.long -> Signed.long
val prec2ndec : Signed.long -> Signed.long
val nbits2ndec : Signed.long -> Signed.long
val precdbl : Signed.long -> Signed.long
val divsbil : Signed.long -> Signed.long
val remsbil : Signed.long -> Signed.long

val fp_red :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_neg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_halve :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_center :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_center_i :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_addmul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_mulu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_muls :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_inv :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_invsafe :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fp_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fp_divu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val flx_mulu :
  ('kind, 'structure) typ -> pari_ulong -> pari_ulong -> ('kind, 'structure) typ

val get_f2x_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val get_f2x_var : ('kind, 'structure) typ -> Signed.long
val get_f2x_degree : ('kind, 'structure) typ -> Signed.long
val get_f2xqx_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val get_f2xqx_var : ('kind, 'structure) typ -> Signed.long
val get_f2xqx_degree : ('kind, 'structure) typ -> Signed.long
val get_flx_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val get_flx_var : ('kind, 'structure) typ -> Signed.long
val get_flx_degree : ('kind, 'structure) typ -> Signed.long
val get_flxqx_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val get_flxqx_var : ('kind, 'structure) typ -> Signed.long
val get_flxqx_degree : ('kind, 'structure) typ -> Signed.long
val get_fpx_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val get_fpx_var : ('kind, 'structure) typ -> Signed.long
val get_fpx_degree : ('kind, 'structure) typ -> Signed.long
val get_fpxqx_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val get_fpxqx_var : ('kind, 'structure) typ -> Signed.long
val get_fpxqx_degree : ('kind, 'structure) typ -> Signed.long

val submulii :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val mulsubii :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val submuliu :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val addmuliu :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val submuliu_inplace :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val addmuliu_inplace :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val lincombii :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

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
val qfb_is_qfi : ('kind, 'structure) typ -> int
val cbrtr_abs : ('kind, 'structure) typ -> ('kind, 'structure) typ
val cbrtr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val sqrtnr : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val logint : ('kind, 'structure) typ -> ('kind, 'structure) typ -> Signed.long
val ulogint : pari_ulong -> pari_ulong -> pari_ulong
val ismpzero : ('kind, 'structure) typ -> int
val isintzero : ('kind, 'structure) typ -> int
val isint1 : ('kind, 'structure) typ -> int
val isintm1 : ('kind, 'structure) typ -> int
val equali1 : ('kind, 'structure) typ -> int
val equalim1 : ('kind, 'structure) typ -> int
val is_pm1 : ('kind, 'structure) typ -> int
val is_bigint : ('kind, 'structure) typ -> int
val odd : Signed.long -> int
val both_odd : Signed.long -> Signed.long -> int
val isonstack : ('kind, 'structure) typ -> int
val dbllog2r : ('kind, 'structure) typ -> float

val mul_content :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val inv_content : ('kind, 'structure) typ -> ('kind, 'structure) typ

val div_content :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val mul_denom :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val constant_coeff : ('kind, 'structure) typ -> ('kind, 'structure) typ
val leading_coeff : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flx_lead : ('kind, 'structure) typ -> pari_ulong
val flx_constant : ('kind, 'structure) typ -> pari_ulong
val degpol : ('kind, 'structure) typ -> Signed.long
val lgpol : ('kind, 'structure) typ -> Signed.long
val lgcols : ('kind, 'structure) typ -> Signed.long
val nbrows : ('kind, 'structure) typ -> Signed.long
val truecoef : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val zxq_sqr :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_coeff :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val rgx_renormalize : ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgx_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val rgxqx_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rgxqx_rem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpx_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val flx_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flx_div_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val f2x_div :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpv_fpc_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val pol0_flx : Signed.long -> ('kind, 'structure) typ
val pol1_flx : Signed.long -> ('kind, 'structure) typ
val polx_flx : Signed.long -> ('kind, 'structure) typ
val zero_zx : Signed.long -> ('kind, 'structure) typ
val polx_zx : Signed.long -> ('kind, 'structure) typ
val zx_shift : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val zero_f2x : Signed.long -> ('kind, 'structure) typ
val pol0_f2x : Signed.long -> ('kind, 'structure) typ
val pol1_f2x : Signed.long -> ('kind, 'structure) typ
val polx_f2x : Signed.long -> ('kind, 'structure) typ
val f2x_equal1 : ('kind, 'structure) typ -> int
val f2x_equal : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int
val f2x_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2v_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flv_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flx_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val vecsmall_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zx_equal1 : ('kind, 'structure) typ -> int
val zx_is_monic : ('kind, 'structure) typ -> int

val zx_renormalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val fpx_renormalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val fpxx_renormalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val fpxqx_renormalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2x_renormalize :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val f2xx_shift :
  ('kind, 'structure) typ ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) typ

val f2v_to_f2x :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val sturm : ('kind, 'structure) typ -> Signed.long
val gval : ('kind, 'structure) typ -> Signed.long -> Signed.long
val rgx_shift_inplace_init : Signed.long -> unit

val rgx_shift_inplace :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zc_to_zc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_to_zv : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zx_to_zv : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val zv_to_zx : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zm_to_zxv :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zero_zm : Signed.long -> Signed.long -> ('kind, 'structure) typ
val zero_zv : Signed.long -> ('kind, 'structure) typ
val zm_transpose : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zm_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zm_row : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val zc_hnfrem :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_hnfrem :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zm_lll :
  ('kind, 'structure) typ -> float -> Signed.long -> ('kind, 'structure) typ

val rgm_dimensions :
  ('kind, 'structure) typ ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val rgm_shallowcopy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f2m_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val f3m_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val flm_copy : ('kind, 'structure) typ -> ('kind, 'structure) typ
val zv_dvd : ('kind, 'structure) typ -> ('kind, 'structure) typ -> int

val zm_zv_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val zv_zv_mod :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecmodii :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val vecmoduu :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fq_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fq_to_fpxq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val rg_to_fq :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val gener_fq_local :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val random_fq :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val fpxqx_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val flxqx_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxqx_div_pre :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) typ

val f2xqx_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxy_fq_evaly :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ

val fqx_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_neg :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_fp_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_fq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_mulu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_powu :
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_halve :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_get_red :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_rem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_divrem :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fqx_div_by_x_x :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fqx_halfgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_gcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_extgcd :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ Ctypes_static.ptr ->
  ('kind, 'structure) typ

val fqx_normalize :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_deriv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_integ :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_factor :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_factor_squarefree :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_ddf :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_degfact :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_roots :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqx_to_mod :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_div :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_inv :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_invsafe :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_sqr :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxq_pow :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxn_expint :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxn_exp :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxn_inv :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxn_mul :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fqxn_sqr :
  ('kind, 'structure) typ ->
  Signed.long ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val fpxq_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val flxq_add :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val flxq_sub :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  pari_ulong ->
  ('kind, 'structure) typ

val f2x_coeff : ('kind, 'structure) typ -> Signed.long -> pari_ulong
val f2x_clear : ('kind, 'structure) typ -> Signed.long -> unit
val f2x_set : ('kind, 'structure) typ -> Signed.long -> unit
val f2x_flip : ('kind, 'structure) typ -> Signed.long -> unit
val f2v_coeff : ('kind, 'structure) typ -> Signed.long -> pari_ulong
val f2v_clear : ('kind, 'structure) typ -> Signed.long -> unit
val f2v_set : ('kind, 'structure) typ -> Signed.long -> unit
val f2v_flip : ('kind, 'structure) typ -> Signed.long -> unit

val f2m_coeff :
  ('kind, 'structure) typ -> Signed.long -> Signed.long -> pari_ulong

val f2m_clear : ('kind, 'structure) typ -> Signed.long -> Signed.long -> unit
val f2m_set : ('kind, 'structure) typ -> Signed.long -> Signed.long -> unit
val f2m_flip : ('kind, 'structure) typ -> Signed.long -> Signed.long -> unit

val f3m_coeff :
  ('kind, 'structure) typ -> Signed.long -> Signed.long -> pari_ulong

val f3m_set :
  ('kind, 'structure) typ -> Signed.long -> Signed.long -> pari_ulong -> unit

val matpascal : Signed.long -> ('kind, 'structure) typ
val z_issquare : ('kind, 'structure) typ -> Signed.long
val z_ispower : ('kind, 'structure) typ -> pari_ulong -> Signed.long
val sqrti : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gaddgs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gcmpgs : ('kind, 'structure) typ -> Signed.long -> int
val gequalgs : ('kind, 'structure) typ -> Signed.long -> int
val gmaxsg : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val gminsg : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ
val gmulgs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gmulgu : ('kind, 'structure) typ -> pari_ulong -> ('kind, 'structure) typ
val gsubgs : ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ
val gdivsg : Signed.long -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmax_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val gmin_shallow :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val cxnorm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val quadnorm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val quad_disc : ('kind, 'structure) typ -> ('kind, 'structure) typ

val qfb_disc3 :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ

val qfb_disc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val sqrfrac : ('kind, 'structure) typ -> ('kind, 'structure) typ
val normalize_frac : ('kind, 'structure) typ -> unit
val powis : Signed.long -> ('kind, 'structure) typ
val mpexpz : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val mplogz : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val mpcosz : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val mpsinz : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val gnegz : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val gabsz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> unit

val gaddz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val gsubz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val gmulz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val gdivz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val gdiventz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val gmodz :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val gmul2nz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> unit

val gshiftz :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ -> unit

val ell_get_a1 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_a2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_a3 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_a4 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_a6 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_b2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_b4 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_b6 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_b8 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_c4 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_c6 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ell_get_type : ('kind, 'structure) typ -> Signed.long
val ellff_get_field : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellff_get_a4a6 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellqp_get_zero : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellqp_get_prec : ('kind, 'structure) typ -> Signed.long
val ellqp_get_p : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellr_get_prec : ('kind, 'structure) typ -> Signed.long
val ellr_get_sign : ('kind, 'structure) typ -> Signed.long
val ellnf_get_nf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val ellnf_get_bnf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val checkell_i : ('kind, 'structure) typ -> int
val ell_is_inf : ('kind, 'structure) typ -> int
val modpr_get_pr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val modpr_get_p : ('kind, 'structure) typ -> ('kind, 'structure) typ
val modpr_get_t : ('kind, 'structure) typ -> ('kind, 'structure) typ
val pr_get_p : ('kind, 'structure) typ -> ('kind, 'structure) typ
val pr_get_gen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val pr_get_e : ('kind, 'structure) typ -> Signed.long
val pr_get_f : ('kind, 'structure) typ -> Signed.long
val pr_get_tau : ('kind, 'structure) typ -> ('kind, 'structure) typ
val pr_is_inert : ('kind, 'structure) typ -> int
val pr_norm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val upr_norm : ('kind, 'structure) typ -> pari_ulong
val nf_get_varn : ('kind, 'structure) typ -> Signed.long
val nf_get_pol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_degree : ('kind, 'structure) typ -> Signed.long
val nf_get_r1 : ('kind, 'structure) typ -> Signed.long
val nf_get_r2 : ('kind, 'structure) typ -> Signed.long
val nf_get_index : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_m : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_g : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_roundg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_tr : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_diff : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_ramified_primes : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_roots : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_zkprimpart : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_zkden : ('kind, 'structure) typ -> ('kind, 'structure) typ
val nf_get_invzk : ('kind, 'structure) typ -> ('kind, 'structure) typ
val cyc_get_expo : ('kind, 'structure) typ -> ('kind, 'structure) typ
val abgrp_get_no : ('kind, 'structure) typ -> ('kind, 'structure) typ
val abgrp_get_cyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val abgrp_get_gen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_nf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_clgp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_no : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_cyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_gen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_reg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_logfu : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_sunits : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_tuu : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnf_get_tun : ('kind, 'structure) typ -> Signed.long
val bnf_get_fu_nocheck : ('kind, 'structure) typ -> ('kind, 'structure) typ

val nfv_to_scalar_or_alg :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val bnf_get_fu : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnr_get_bnf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnr_get_bid : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnr_get_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnr_get_nf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnr_get_clgp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnr_get_no : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnr_get_cyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnr_get_gen_nocheck : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bnr_get_gen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val locs_get_cyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val locs_get_lsprk : ('kind, 'structure) typ -> ('kind, 'structure) typ
val locs_get_lgenfil : ('kind, 'structure) typ -> ('kind, 'structure) typ
val locs_get_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val locs_get_famod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val locs_get_m_infty : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_basis : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_bnf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_nf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_zm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_modp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_s : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_dldata : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_sfu : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_cyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_hnf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_u : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_ui : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_m0 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_u0 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_r1 : ('kind, 'structure) typ -> Signed.long
val gchar_get_r2 : ('kind, 'structure) typ -> Signed.long
val gchar_get_loccyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gchar_get_nc : ('kind, 'structure) typ -> Signed.long
val gchar_get_ns : ('kind, 'structure) typ -> Signed.long
val gchar_get_nm : ('kind, 'structure) typ -> Signed.long
val gchar_get_evalprec : ('kind, 'structure) typ -> Signed.long
val gchar_get_prec : ('kind, 'structure) typ -> Signed.long
val gchar_get_nfprec : ('kind, 'structure) typ -> Signed.long
val gchar_set_evalprec : ('kind, 'structure) typ -> Signed.long -> unit
val gchar_set_prec : ('kind, 'structure) typ -> Signed.long -> unit

val gchar_copy_precs :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val gchar_set_nfprec : ('kind, 'structure) typ -> Signed.long -> unit
val gchar_get_ntors : ('kind, 'structure) typ -> Signed.long
val gchar_get_nfree : ('kind, 'structure) typ -> Signed.long
val gchar_get_nalg : ('kind, 'structure) typ -> Signed.long
val gchar_set_basis : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val gchar_set_nf : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val gchar_set_ntors : ('kind, 'structure) typ -> Signed.long -> unit
val gchar_set_nfree : ('kind, 'structure) typ -> Signed.long -> unit
val gchar_set_nalg : ('kind, 'structure) typ -> Signed.long -> unit
val gchar_set_cyc : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val gchar_set_huui :
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  ('kind, 'structure) typ ->
  unit

val gchar_set_m0 : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val gchar_set_u0 : ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit
val bid_get_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_ideal : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_arch : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_grp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_fact : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_fact2 : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_sprk : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_sarch : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_archp : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_u : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_no : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_cyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_gen_nocheck : ('kind, 'structure) typ -> ('kind, 'structure) typ
val bid_get_gen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_n : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_fan : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_no : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_cyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_gen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_conreycyc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_conreygen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_ui : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_u : ('kind, 'structure) typ -> ('kind, 'structure) typ
val znstar_get_pe : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_pol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_p : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_e : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_mod : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_roots : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_invvdm : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_den : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_group : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_gen : ('kind, 'structure) typ -> ('kind, 'structure) typ
val gal_get_orders : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_degree : ('kind, 'structure) typ -> Signed.long
val rnf_get_nfdegree : ('kind, 'structure) typ -> Signed.long
val rnf_get_absdegree : ('kind, 'structure) typ -> Signed.long
val rnf_get_idealdisc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_k : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_alpha : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_nf : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_nfzk : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_polabs : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_pol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_disc : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_index : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_ramified_primes : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_varn : ('kind, 'structure) typ -> Signed.long
val rnf_get_nfpol : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_nfvarn : ('kind, 'structure) typ -> Signed.long
val rnf_get_zk : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_map : ('kind, 'structure) typ -> ('kind, 'structure) typ
val rnf_get_invzk : ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealred :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val idealchineseinit :
  ('kind, 'structure) typ -> ('kind, 'structure) typ -> ('kind, 'structure) typ

val closure_arity : ('kind, 'structure) typ -> Signed.long
val closure_is_variadic : ('kind, 'structure) typ -> Signed.long
val closure_codestr : ('kind, 'structure) typ -> string
val closure_get_code : ('kind, 'structure) typ -> ('kind, 'structure) typ
val closure_get_oper : ('kind, 'structure) typ -> ('kind, 'structure) typ
val closure_get_data : ('kind, 'structure) typ -> ('kind, 'structure) typ
val closure_get_dbg : ('kind, 'structure) typ -> ('kind, 'structure) typ
val closure_get_text : ('kind, 'structure) typ -> ('kind, 'structure) typ
val closure_get_frame : ('kind, 'structure) typ -> ('kind, 'structure) typ
val err_get_num : ('kind, 'structure) typ -> Signed.long

val err_get_compo :
  ('kind, 'structure) typ -> Signed.long -> ('kind, 'structure) typ

val pari_err_bug : string -> unit
val pari_err_constpol : string -> unit

val pari_err_coprime :
  string -> ('kind, 'structure) typ -> ('kind, 'structure) typ -> unit

val with_stack_clean :
  (unit -> ('kind, 'structure) typ) -> ('kind, 'structure) typ

val with_stack_clean_opt :
  (unit -> ('kind, 'structure) typ option) -> ('kind, 'structure) typ option

val with_stack_clean6 :
  ?av:pari_ulong ->
  (unit ->
  ('k1, 's1) typ
  * ('k2, 's2) typ
  * ('k3, 's3) typ
  * ('k4, 's4) typ
  * ('k5, 's5) typ
  * ('k6, 's6) typ) ->
  ('k1, 's1) typ
  * ('k2, 's2) typ
  * ('k3, 's3) typ
  * ('k4, 's4) typ
  * ('k5, 's5) typ
  * ('k6, 's6) typ

val gentobytes : ('kind, 'structure) typ -> bytes
