type pari_ulong = Unsigned.ULong.t

val pari_ulong : pari_ulong Ctypes.typ

type ('kind, 'structure) t

val t : ('kind, 'structure) t Ctypes.typ

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
type ring
type field

module rec Complex : sig
  type complex = private Complex
  type nonrec t = (complex, field) t

  val inv : t -> t
  val add : t -> t -> t
  val create : re:Real.t -> im:Real.t -> t
  val to_string : t -> string
end

and Real : sig
  type real = private Real
  type nonrec t = (real, field) t

  external inj_complex : t -> Complex.t = "%identity"
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
  type rational = private Rational
  type ('a, 'b) p := ('a, 'b) t
  type nonrec t = (rational, field) t
  type nonrec ring = (rational, ring) p

  external inj_ring : t -> ring = "%identity"
  external inj_real : t -> Real.t = "%identity"
  external inj_complex : t -> Complex.t = "%identity"
  val shift : t -> int -> t
end

and Integer : sig
  type integer = private Integer
  type nonrec t = (integer, ring) t

  external inj_rat : t -> Rational.t = "%identity"
  external inj_real : t -> Real.t = "%identity"
  external inj_complex : t -> Complex.t = "%identity"
  val of_int : int -> t
  val of_signed : Signed.long -> t
  val equal : t -> t -> bool
  val shift : t -> int -> t
  val sqrt : t -> Real.t
  val zero : t
  val mul : t -> t -> t
  val add : t -> t -> t
  val sub : t -> t -> t
  val neg : t -> t
  val pow : t -> t -> t
  val modulo : t -> t -> t
  val of_string : string -> t option
  val to_string : t -> string
  val random_prime : bits_amount:int -> t

  module Infix : sig
    val ( * ) : t -> t -> t
    val ( + ) : t -> t -> t
    val ( - ) : t -> t -> t
    val ( ~- ) : t -> t
    val ( mod ) : t -> t -> t
    val ( = ) : t -> t -> bool
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

module Set : sig
  type nonrec 'a t constraint 'a = ('b, 'c) t

  val length : 'a t -> Signed.Long.t
  val search : 'a t -> 'a -> Signed.Long.t -> Signed.Long.t
end

module Vector : sig
  type nonrec ('a, 'b) t
    constraint 'a = ('c, 'd) t constraint 'b = [< `COL | `ROW ]

  val length : ('a, 'b) t -> int
  val of_array : 'a array -> ('a, [ `ROW ]) t
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
  type nonrec 'a t constraint 'a = ('b, 'c) t

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

type finite_field = private Finite_field

module rec Polynomial : sig
  type polynomial
  type ('a, 'b) p := ('a, 'b) t
  type 'a t = (polynomial, ring) p constraint 'a = ('b, ring) p

  val to_string : 'a t -> string
  val mul : 'a t -> 'a t -> 'a t
  val div : 'a t -> 'a t -> 'a t
  val equal : 'a t -> 'a t -> bool
  val add : 'a t -> 'a t -> 'a t
  val sub : 'a t -> 'a t -> 'a t
  val neg : 'a t -> 'a t
  val eval : 'a t -> 'a -> 'a
  val degree : 'a t -> int
  val get_coeff : 'a t -> int -> 'a

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
      # (Polynomial.to_string q);;
      - : string = "x^3 - 111*x^2 + 6064*x - 189804"
      # let zero = Polynomial.create [| Integer.of_int 0 |];;
      val zero : Integer.t Polynomial.t = <abstr>
      # let qq = Polynomial.create [| q; q; zero; zero |];;
      val qq : Integer.t Polynomial.t Polynomial.t = <abstr>
      # (Polynomial.to_string qq);;
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
      val qmin : ('a, ring) t Polynomial.t = <abstr>
      # Polynomial.to_string qmin;;
      - : string = "x^3 - x^2 - 60*x - 364"
      # Number_field.(are_isomorphic (create q) (create qmin));
      - : bool = true
      ]} *)

  val ( .%[] ) : 'a t -> int -> 'a

  val roots_ff :
    (finite_field, ring) p t -> ((finite_field, field) p, [ `ROW ]) Vector.t
end

and Fp : sig
  type t = Integer.t

  val add : t -> t -> modulo:t -> t
  val pow : t -> exponent:t -> modulo:t -> t
end

and Finite_field : sig
  type ('a, 'b) p := ('a, 'b) t
  type t = (finite_field, field) p

  external inj_ring : t -> (finite_field, ring) p = "%identity"
  external inj_field : (finite_field, ring) p -> t = "%identity"
  val generator : order:Integer.t -> t
  val prime_field_element : Integer.t -> p:Integer.t -> t
  val finite_field_element : Integer.t array -> t -> t

  val create : p:int -> degree:int -> (finite_field, ring) p Polynomial.t
  (** [create p degree] returns a monic irreducible polynomial of the
      given [degree] over F_p[X]. *)

  val generator_from_irreducible_polynomial :
    (finite_field, ring) p Polynomial.t -> t

  val residue_class : t -> (finite_field, ring) p Polynomial.t
  val equal : t -> t -> bool
  val add : t -> t -> t
  val mul : t -> t -> t
  val pow : t -> Integer.t -> t
  val random : t -> t
  val zero : t -> t
  val one : t -> t

  val extend :
    (finite_field, field) p ->
    [< `Degree of int | `Quotient of (finite_field, ring) p Polynomial.t ] ->
    t
  (** extend the field {m K} of definition of {m a} by a root of the polynomial
     {m P\in K[X]} assumed to be irreducible over {m K}.  Return {m [r, m]} where {m r}
     is a root of {m P} in the extension field {m L} and {m m} is a map from {m K} to {m L},
     see [ffmap].
     If {m v} is given, the variable name is used to display the generator of {m L},
     else the name of the variable of {m P} is used.
     A generator of {m L} can be recovered using [b=ffgen(r)].
     The image of {m P} in {m L[X]} can be recovered using [PL=ffmap(m,P)]. *)

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

module Number_field : sig
  type number_field = private Number_field
  type structure
  type nonrec t = (number_field, field) t

  val create : Rational.ring Polynomial.t -> structure
  (** [create p] returns the number field Q(X)/(p) for a monic
      irreducible polynomial [p] over the field Q of the rationals. *)

  val are_isomorphic : structure -> structure -> bool
  (** [are_isomorphic a b] returns true if and only if number
      fields [a] and [b] are isomorphic. *)

  val sign : structure -> Signed.Long.t * Signed.Long.t
  val discriminant : structure -> Integer.t
  val z_basis : structure -> (t, [ `ROW ]) Vector.t
  val elt : Rational.t array -> t
  val add : structure -> t -> t -> t
  val mul : structure -> t -> t -> t
  val equal : t -> t -> bool
  val divrem : structure -> t -> t -> t * t
  val ideal_norm : structure -> t -> Integer.t

  val splitting :
    [< `F of structure | `P of Integer.t Polynomial.t ] ->
    Integer.t Polynomial.t
  (** [splitting (nf|p)] given the number field [nf = Q(x)/(p)]
      or polynomial [p], returns the polynomial over Q for the
      splitting field of [p], that is the smallest field over
      which [p] is totally split. *)

  module Infix : sig
    val ( = ) : t -> t -> bool
  end
end

module Elliptic_curve : sig
  type elliptic_curve
  type ('a, 'b) p := ('a, 'b) t
  type 'a structure constraint 'a = ('b, field) t
  type nonrec 'a t = (elliptic_curve, group) t constraint 'a = ('b, field) t

  val create :
    ?a1:'a ->
    ?a2:'a ->
    ?a3:'a ->
    ?a4:'a ->
    ?a6:'a ->
    ?dom:'a ->
    unit ->
    'a structure option
  (** [create ?a1 ?a2 ?a3 ?a4 ?a6] defines the curve
      {math Y^2 + a_1 XY + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6}
      Returns [None] if the input coefficients do not define an elliptic curve
      over the field from the coefficients. *)

  val get_a1 : 'a structure -> 'a
  val get_a2 : 'a structure -> 'a
  val get_a3 : 'a structure -> 'a
  val get_a4 : 'a structure -> 'a
  val get_a6 : 'a structure -> 'a
  val of_coordinates : x:'a -> y:'a -> 'a t

  val order : 'a structure -> Integer.t
  (** {@ocaml[
      # let g = Finite_field.generator ~order:(Integer.of_int 625) (* 5^4 *);;
      val g : Finite_field.t = <abstr>
      # let ell = Option.get (Elliptic_curve.create ~a4:(Finite_field.pow g (Integer.of_int 4)) ~a6:(Finite_field.pow g (Integer.of_int 6)) ());;
      val ell : Finite_field.t Elliptic_curve.structure = <abstr>
      # Integer.(to_string (Elliptic_curve.order ell));;
      - : string = "675"
      ]} *)

  val discriminant : 'a structure -> 'a
  val j_invariant : 'a structure -> 'a
  val random : 'a structure -> 'a t

  val weil_pairing_ff :
    Finite_field.t structure ->
    l:Integer.t ->
    p:Finite_field.t t ->
    q:Finite_field.t t ->
    Finite_field.t
  (** [weil_pairing_ff ell ~l ~p ~q] returns the Weil pairing of the two points
      of [l]-torsion [p] and [q] on the elliptic curve [ell].

      {@ocaml[
      # let l = Integer.of_int 3;;
      val l : Integer.t = <abstr>
      # let ord = Integer.of_int 103;;
      val ord : Integer.t = <abstr>
      # let ell = Option.get (Elliptic_curve.create ~a3:(Finite_field.prime_field_element (Integer.of_int 1) ~p:ord) ());;
      val ell : Finite_field.t Elliptic_curve.structure = <abstr>
      # let (p, q) = Elliptic_curve.(of_coordinates ~x:(Finite_field.prime_field_element (Integer.of_int 0) ~p:ord) ~y:(Finite_field.prime_field_element (Integer.of_int 0) ~p:ord), of_coordinates ~x:(Finite_field.prime_field_element (Integer.of_int 57) ~p:ord) ~y:(Finite_field.prime_field_element (Integer.of_int 46) ~p:ord));;
      val p : Finite_field.t Elliptic_curve.t = <abstr>
      val q : Finite_field.t Elliptic_curve.t = <abstr>
      # let scalar = (Elliptic_curve.weil_pairing_ff ell ~l ~p ~q);;
      val scalar : Finite_field.t = <abstr>
      # Finite_field.to_string scalar (* 56 mod 103 *);;
      - : string = "56"
      # Finite_field.(to_string (pow scalar l));;
      - : string = "1"
      ]}
      *)

  val l_division_polynomial :
    ('a, field) p structure -> l:Signed.Long.t -> ('a, ring) p Polynomial.t
  (** {@ocaml[
      # let g = Finite_field.generator ~order:(Integer.of_int 625) (* 5^4 *);;
      val g : Finite_field.t = <abstr>
      # let ell = Option.get (Elliptic_curve.create ~a6:(Finite_field.pow g (Integer.of_int 6)) ());;
      val ell : Finite_field.t Elliptic_curve.structure = <abstr>
      # let pdiv7 = (Elliptic_curve.l_division_polynomial ell ~l:(Signed.Long.of_int 7));;
      val pdiv7 : (finite_field, ring) t Polynomial.t = <abstr>
      # Polynomial.to_string pdiv7;;
      - : string =
      "2*x^24 + (3*x^3 + x + 2)*x^21 + (x^3 + x^2 + x + 2)*x^18 + (2*x^3 + 2*x^2 + 4*x)*x^15 + (2*x^3 + 4*x^2 + 4*x + 1)*x^12 + (3*x^3 + 4*x^2 + 1)*x^9 + (4*x^3 + x^2 + 4)*x^6 + (2*x^3 + 3*x^2 + 3*x + 4)*x^3 + (x^3 + x^2)"
      # Polynomial.degree pdiv7 = ((7 * 7 - 1) / 2);;
      - : bool = true
      ]} *)

  val to_string : 'a t -> string
  val add : 'a structure -> 'a t -> 'a t -> 'a t
  val sub : 'a structure -> 'a t -> 'a t -> 'a t
  val mul : 'a structure -> n:Integer.t -> p:'a t -> 'a t
  val equal : 'a t -> 'a t -> bool

  val generators_ff :
    Finite_field.t structure -> (Finite_field.t t, [ `ROW ]) Vector.t

  val zero : 'a structure -> 'a t
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
  (('kind, 'structure) t, forprime_t Ctypes.structure) Ctypes.field

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
  (('kind, 'structure) t, forprime_t Ctypes.structure) Ctypes.field

val forcomposite_t : forcomposite_t Ctypes.structure Ctypes.typ
val forcomposite_t_first : (int, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_b :
  (('kind, 'structure) t, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_n :
  (('kind, 'structure) t, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_p :
  (('kind, 'structure) t, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_T :
  (forprime_t Ctypes.structure, forcomposite_t Ctypes.structure) Ctypes.field

val forvec_t : forvec_t Ctypes.structure Ctypes.typ
val forvec_t_first : (Signed.long, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_a :
  ( ('kind, 'structure) t Ctypes_static.ptr,
    forvec_t Ctypes.structure )
  Ctypes.field

val forvec_t_m :
  ( ('kind, 'structure) t Ctypes_static.ptr,
    forvec_t Ctypes.structure )
  Ctypes.field

val forvec_t_M :
  ( ('kind, 'structure) t Ctypes_static.ptr,
    forvec_t Ctypes.structure )
  Ctypes.field

val forvec_t_n : (Signed.long, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_next :
  ( (forvec_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t)
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
  (('kind, 'structure) t, forpart_t Ctypes.structure) Ctypes.field

val forperm_t : forperm_t Ctypes.structure Ctypes.typ
val forperm_t_k : (Signed.long, forperm_t Ctypes.structure) Ctypes.field
val forperm_t_first : (Signed.long, forperm_t Ctypes.structure) Ctypes.field

val forperm_t_v :
  (('kind, 'structure) t, forperm_t Ctypes.structure) Ctypes.field

val forsubset_t : forsubset_t Ctypes.structure Ctypes.typ
val forsubset_t_n : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_k : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_all : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_first : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field

val forsubset_t_v :
  (('kind, 'structure) t, forsubset_t Ctypes.structure) Ctypes.field

val pari_plot : pari_plot Ctypes.structure Ctypes.typ

val pari_plot_draw :
  ( (pari_plot Ctypes.structure Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
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
val genbin_x : (('kind, 'structure) t, genbin Ctypes.structure) Ctypes.field
val genbin_base : (('kind, 'structure) t, genbin Ctypes.structure) Ctypes.field

val genbin_rebase :
  ( (('kind, 'structure) t -> Signed.long -> unit) Ctypes_static.static_funptr,
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
  (('kind, 'structure) t, pari_parsestate Ctypes.structure) Ctypes.field

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
  (('kind, 'structure) t, pari_global_state Ctypes.structure) Ctypes.field

val pari_global_state_seadata :
  (('kind, 'structure) t, pari_global_state Ctypes.structure) Ctypes.field

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
  (('kind, 'structure) t, pari_thread Ctypes.structure) Ctypes.field

val mt_state : mt_state Ctypes.structure Ctypes.typ

val mt_state_worker :
  (('kind, 'structure) t, mt_state Ctypes.structure) Ctypes.field

val mt_state_pending :
  (('kind, 'structure) t, mt_state Ctypes.structure) Ctypes.field

val mt_state_workid : (Signed.long, mt_state Ctypes.structure) Ctypes.field
val pari_mt : pari_mt Ctypes.structure Ctypes.typ

val pari_mt_mt :
  (mt_state Ctypes.structure, pari_mt Ctypes.structure) Ctypes.field

val pari_mt_get :
  ( (mt_state Ctypes.structure Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    pari_mt Ctypes.structure )
  Ctypes.field

val pari_mt_submit :
  ( (mt_state Ctypes.structure Ctypes_static.ptr ->
    Signed.long ->
    ('kind, 'structure) t ->
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
  (('kind, 'structure) t, parfor_iter Ctypes.structure) Ctypes.field

val parfor_iter_pt :
  (pari_mt Ctypes.structure, parfor_iter Ctypes.structure) Ctypes.field

val parfor_t : parfor_t Ctypes.structure Ctypes.typ
val parfor_t_a : (('kind, 'structure) t, parfor_t Ctypes.structure) Ctypes.field
val parfor_t_b : (('kind, 'structure) t, parfor_t Ctypes.structure) Ctypes.field

val parfor_t_iter :
  (parfor_iter Ctypes.structure, parfor_t Ctypes.structure) Ctypes.field

val parforeach_t : parforeach_t Ctypes.structure Ctypes.typ

val parforeach_t_x :
  (('kind, 'structure) t, parforeach_t Ctypes.structure) Ctypes.field

val parforeach_t_W :
  (('kind, 'structure) t, parforeach_t Ctypes.structure) Ctypes.field

val parforeach_t_i : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_l : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field

val parforeach_t_iter :
  (parfor_iter Ctypes.structure, parforeach_t Ctypes.structure) Ctypes.field

val parforprime_t : parforprime_t Ctypes.structure Ctypes.typ

val parforprime_t_v :
  (('kind, 'structure) t, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_forprime :
  (forprime_t Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_iter :
  (parfor_iter Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforvec_t : parforvec_t Ctypes.structure Ctypes.typ

val parforvec_t_v :
  (('kind, 'structure) t, parforvec_t Ctypes.structure) Ctypes.field

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
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dT :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_T0 :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_unscale :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dK :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_index :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_basis :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_r1 : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_basden :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dTP :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dTE :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dKP :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_dKE :
  (('kind, 'structure) t, nfmaxord_t Ctypes.structure) Ctypes.field

val nfmaxord_t_certify : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field
val qfr_data : qfr_data Ctypes.structure Ctypes.typ
val qfr_data_D : (('kind, 'structure) t, qfr_data Ctypes.structure) Ctypes.field

val qfr_data_sqrtD :
  (('kind, 'structure) t, qfr_data Ctypes.structure) Ctypes.field

val qfr_data_isqrtD :
  (('kind, 'structure) t, qfr_data Ctypes.structure) Ctypes.field

val fp_chk_fun : fp_chk_fun Ctypes.structure Ctypes.typ

val fp_chk_fun_f :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_init :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_post :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_data :
  (unit Ctypes_static.ptr, fp_chk_fun Ctypes.structure) Ctypes.field

val fp_chk_fun_skipfirst :
  (Signed.long, fp_chk_fun Ctypes.structure) Ctypes.field

val zlog_s : zlog_s Ctypes.structure Ctypes.typ
val zlog_s_bid : (('kind, 'structure) t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_P : (('kind, 'structure) t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_k : (('kind, 'structure) t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_sprk : (('kind, 'structure) t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_archp : (('kind, 'structure) t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_mod : (('kind, 'structure) t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_U : (('kind, 'structure) t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_hU : (Signed.long, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_no2 : (int, zlog_s Ctypes.structure) Ctypes.field
val bb_group : bb_group Ctypes.structure Ctypes.typ

val bb_group_mul :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_pow :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_rand :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_hash :
  ( (('kind, 'structure) t -> pari_ulong) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal :
  ( (('kind, 'structure) t -> ('kind, 'structure) t -> int)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal1 :
  ( (('kind, 'structure) t -> int) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_easylog :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_field : bb_field Ctypes.structure Ctypes.typ

val bb_field_red :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_add :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_mul :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_neg :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_inv :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_equal0 :
  ( (('kind, 'structure) t -> int) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_s :
  ( (unit Ctypes_static.ptr -> Signed.long -> ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_algebra : bb_algebra Ctypes.structure Ctypes.typ

val bb_algebra_red :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_add :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sub :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_mul :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sqr :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_one :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_zero :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_ring : bb_ring Ctypes.structure Ctypes.typ

val bb_ring_add :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_mul :
  ( (unit Ctypes_static.ptr ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t ->
    ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_sqr :
  ( (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
    Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val buchimag :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val buchreal :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zidealstar :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zidealstarinit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zidealstarinitgen :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val factmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mpbern : Signed.long -> Signed.long -> unit

val simplefactmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val listkill : ('kind, 'structure) t -> unit

val isprincipalforce :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val isprincipalgen :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val isprincipalgenforce :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2ms_ker : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f2ms_to_f2m : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f2c_to_zc : ('kind, 'structure) t -> ('kind, 'structure) t
val f2c_to_mod : ('kind, 'structure) t -> ('kind, 'structure) t

val f2m_f2c_gauss :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2m_f2c_invimage :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2m_f2c_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2m_deplin : ('kind, 'structure) t -> ('kind, 'structure) t
val f2m_det : ('kind, 'structure) t -> pari_ulong
val f2m_det_sp : ('kind, 'structure) t -> pari_ulong

val f2m_gauss :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2m_inv : ('kind, 'structure) t -> ('kind, 'structure) t

val f2m_invimage :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2m_ker : ('kind, 'structure) t -> ('kind, 'structure) t
val f2m_ker_sp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val f2m_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2m_powu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val f2m_rank : ('kind, 'structure) t -> Signed.long
val f2m_row : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val f2m_rowslice :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val f2m_to_f2ms : ('kind, 'structure) t -> ('kind, 'structure) t
val f2m_to_flm : ('kind, 'structure) t -> ('kind, 'structure) t
val f2m_to_zm : ('kind, 'structure) t -> ('kind, 'structure) t
val f2m_to_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val f2m_transpose : ('kind, 'structure) t -> ('kind, 'structure) t
val f2v_add_inplace : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val f2v_and_inplace : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val f2v_dotproduct :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong

val f2v_equal0 : ('kind, 'structure) t -> int
val f2v_hamming : ('kind, 'structure) t -> pari_ulong

val f2v_negimply_inplace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val f2v_or_inplace : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val f2v_slice :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val f2v_subset : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val f2v_to_flv : ('kind, 'structure) t -> ('kind, 'structure) t
val matid_f2m : Signed.long -> ('kind, 'structure) t

val f2x_f2xq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2x_f2xqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2x_frobenius : ('kind, 'structure) t -> ('kind, 'structure) t
val f2x_1_add : ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_deflate : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f2x_degfact : ('kind, 'structure) t -> ('kind, 'structure) t
val f2x_degree : ('kind, 'structure) t -> Signed.long
val f2x_deriv : ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val f2x_eval : ('kind, 'structure) t -> pari_ulong -> pari_ulong

val f2x_even_odd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val f2x_extgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val f2x_gcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_get_red : ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_halfgcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_issquare : ('kind, 'structure) t -> int
val f2x_matfrobenius : ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_recip : ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_rem :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_shift : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f2x_sqr : ('kind, 'structure) t -> ('kind, 'structure) t
val f2x_sqrt : ('kind, 'structure) t -> ('kind, 'structure) t
val f2x_to_f2v : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f2x_to_f2xx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f2x_to_flx : ('kind, 'structure) t -> ('kind, 'structure) t
val f2x_to_zx : ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_valrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val f2xc_to_flxc : ('kind, 'structure) t -> ('kind, 'structure) t
val f2xc_to_zxc : ('kind, 'structure) t -> ('kind, 'structure) t
val f2xv_to_f2m : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f2xv_to_flxv_inplace : ('kind, 'structure) t -> unit
val f2xv_to_zxv_inplace : ('kind, 'structure) t -> unit

val f2xx_f2x_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xx_f2x_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xx_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xx_deriv : ('kind, 'structure) t -> ('kind, 'structure) t

val f2xx_renormalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val f2xx_to_kronecker :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val f2xx_to_flxx : ('kind, 'structure) t -> ('kind, 'structure) t
val f2xx_to_zxx : ('kind, 'structure) t -> ('kind, 'structure) t

val f2xx_to_f2xc :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val f2xxv_to_f2xm :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val f2xxc_to_zxxc : ('kind, 'structure) t -> ('kind, 'structure) t

val f2xy_f2xq_evalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xy_f2xqv_evalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xy_degreex : ('kind, 'structure) t -> Signed.long

val f2xn_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val f2xn_inv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f2xn_red : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val f2xq_artin_schreier :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xq_autpow :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_conjvec :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xq_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xq_invsafe :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xq_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_matrix_pow :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_pow_init :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_pow_table :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xq_sqrt :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xq_sqrt_fast :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_sqrtn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val f2xq_trace : ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong

val f2xqx_f2xq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_f2xq_mul_to_monic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_f2xqxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_f2xqxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_disc :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val f2xqx_extgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val f2xqx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_get_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_halfgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_halfgcd_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val f2xqx_invbarrett :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_ispower :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val f2xqx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_normalize :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_rem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqx_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqxq_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqxq_invsafe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqxq_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqxq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqxq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqxq_autpow :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqxq_auttrace :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqxqv_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val flm_to_f2m : ('kind, 'structure) t -> ('kind, 'structure) t
val flv_to_f2v : ('kind, 'structure) t -> ('kind, 'structure) t
val flx_to_f2x : ('kind, 'structure) t -> ('kind, 'structure) t
val flxc_to_f2xc : ('kind, 'structure) t -> ('kind, 'structure) t
val flxx_to_f2xx : ('kind, 'structure) t -> ('kind, 'structure) t
val flxxc_to_f2xxc : ('kind, 'structure) t -> ('kind, 'structure) t

val kronecker_to_f2xqx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rg_to_f2xq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_to_f2m : ('kind, 'structure) t -> ('kind, 'structure) t
val rgv_to_f2v : ('kind, 'structure) t -> ('kind, 'structure) t
val rgx_to_f2x : ('kind, 'structure) t -> ('kind, 'structure) t
val z_to_f2x : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val zm_to_f2m : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_to_f2v : ('kind, 'structure) t -> ('kind, 'structure) t
val zx_to_f2x : ('kind, 'structure) t -> ('kind, 'structure) t
val zxx_to_f2xx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val const_f2v : Signed.long -> ('kind, 'structure) t

val gener_f2xq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val get_f2xq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  bb_field Ctypes.structure Ctypes_static.ptr

val monomial_f2x : Signed.long -> Signed.long -> ('kind, 'structure) t
val pol1_f2xx : Signed.long -> Signed.long -> ('kind, 'structure) t
val polx_f2xx : Signed.long -> Signed.long -> ('kind, 'structure) t

val random_f2xqx :
  Signed.long -> Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_teichmuller :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val f2xq_ellcard :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_ellgens :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xq_ellgroup :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val f2xq_elltwist :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val f2xqe_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_changepoint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_changepointinv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_dbl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_neg :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_tatepairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqe_weilpairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val get_f2xqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_f2xqe :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val random_f2xqe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f3c_to_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val f3c_to_zc : ('kind, 'structure) t -> ('kind, 'structure) t
val f3m_ker : ('kind, 'structure) t -> ('kind, 'structure) t
val f3m_ker_sp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val f3m_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f3m_row : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f3m_to_flm : ('kind, 'structure) t -> ('kind, 'structure) t
val f3m_to_zm : ('kind, 'structure) t -> ('kind, 'structure) t
val f3m_to_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val f3m_transpose : ('kind, 'structure) t -> ('kind, 'structure) t
val f3v_to_flv : ('kind, 'structure) t -> ('kind, 'structure) t
val f3v_coeff : ('kind, 'structure) t -> Signed.long -> pari_ulong
val f3v_clear : ('kind, 'structure) t -> Signed.long -> unit
val f3v_set : ('kind, 'structure) t -> Signed.long -> pari_ulong -> unit
val flm_to_f3m : ('kind, 'structure) t -> ('kind, 'structure) t
val flv_to_f3v : ('kind, 'structure) t -> ('kind, 'structure) t
val rgm_to_f3m : ('kind, 'structure) t -> ('kind, 'structure) t
val rgv_to_f3v : ('kind, 'structure) t -> ('kind, 'structure) t
val zm_to_f3m : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_to_f3v : ('kind, 'structure) t -> ('kind, 'structure) t
val zero_f3m_copy : Signed.long -> Signed.long -> ('kind, 'structure) t
val zero_f3v : Signed.long -> ('kind, 'structure) t
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
  ('kind, 'structure) t

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
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fle_dbl :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val fle_changepoint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val fle_changepointinv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val fle_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fle_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fle_mulu :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fle_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fle_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fle_tatepairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fle_to_flj : ('kind, 'structure) t -> ('kind, 'structure) t

val fle_weilpairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val flj_add_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flj_changepointinv_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flj_dbl_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flj_mulu_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flj_neg : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flj_to_fle : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flj_to_fle_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val fljv_factorback_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val random_fle : pari_ulong -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val random_fle_pre :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val random_flj_pre :
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flc_to_zc : ('kind, 'structure) t -> ('kind, 'structure) t
val flc_to_zc_inplace : ('kind, 'structure) t -> ('kind, 'structure) t

val flm_flc_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_flc_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_adjoint : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flm_deplin : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flm_det : ('kind, 'structure) t -> pari_ulong -> pari_ulong
val flm_det_sp : ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flm_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_intersect :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_intersect_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_inv : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flm_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_ker : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flm_ker_sp :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val flm_rank : ('kind, 'structure) t -> pari_ulong -> Signed.long
val flm_to_zm : ('kind, 'structure) t -> ('kind, 'structure) t
val flm_to_zm_inplace : ('kind, 'structure) t -> ('kind, 'structure) t
val flv_to_zv : ('kind, 'structure) t -> ('kind, 'structure) t
val fl_to_flx : pari_ulong -> Signed.long -> ('kind, 'structure) t
val fl2_equal1 : ('kind, 'structure) t -> int

val fl2_inv_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fl2_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fl2_norm_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fl2_pow_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fl2_sqr_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fl2_sqrt_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val fl2_sqrtn_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flm_to_flxv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val flm_to_flxx :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val flv_flm_polint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) t

val flv_inv : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flv_inv_inplace : ('kind, 'structure) t -> pari_ulong -> unit

val flv_inv_pre_inplace :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> unit

val flv_inv_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flv_invvandermonde :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flv_polint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) t

val flv_prod : ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flv_prod_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> pari_ulong

val flv_roots_to_pol :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val flv_to_flx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val flx_fl_add :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_fl_mul :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_fl_mul_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_fl_mul_to_monic :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_fl_sub :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_fl2_eval_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_flv_multieval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_flxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_flxq_eval_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_flxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_flxqv_eval_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_frobenius : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_frobenius_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_laplace : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_newton :
  ('kind, 'structure) t -> Signed.long -> pari_ulong -> ('kind, 'structure) t

val flx_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_blocks :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val flx_composedprod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_composedsum :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_convol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_deflate : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val flx_deriv : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flx_diff1 : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_digits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_div_by_x_x :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  ('kind, 'structure) t

val flx_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flx_divrem_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flx_double : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flx_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val flx_eval : ('kind, 'structure) t -> pari_ulong -> pari_ulong -> pari_ulong

val flx_eval_powers_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val flx_eval_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val flx_extgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flx_extgcd_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flx_extresultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong

val flx_extresultant_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong

val flx_fromnewton :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_gcd_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_get_red : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_get_red_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_halfgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_halfgcd_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flx_halfgcd_all_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flx_halfgcd_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_halve : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flx_inflate : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val flx_integ : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_invbarrett :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_invlaplace :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_is_squarefree : ('kind, 'structure) t -> pari_ulong -> int
val flx_is_smooth : ('kind, 'structure) t -> Signed.long -> pari_ulong -> int

val flx_is_smooth_pre :
  ('kind, 'structure) t -> Signed.long -> pari_ulong -> pari_ulong -> int

val flx_matfrobenius :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_matfrobenius_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_mod_xn1 :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_mod_xnm1 :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_neg : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_neg_inplace :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_normalize : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_powu :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_powu_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_recip : ('kind, 'structure) t -> ('kind, 'structure) t
val flx_red : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_rem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_rem_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_renormalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val flx_rescale :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_resultant :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flx_resultant_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val flx_shift : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val flx_splitting :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val flx_sqr : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_sqr_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_translate1 :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_translate1_basecase :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_to_flv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val flx_to_flxx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val flx_to_zx : ('kind, 'structure) t -> ('kind, 'structure) t
val flx_to_zx_inplace : ('kind, 'structure) t -> ('kind, 'structure) t
val flx_triple : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flx_val : ('kind, 'structure) t -> Signed.long

val flx_valrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val flxc_flxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxc_flxqv_eval_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxc_flxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxc_flxq_eval_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxc_eval_powers_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxc_neg : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxc_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxc_to_zxc : ('kind, 'structure) t -> ('kind, 'structure) t

val flxm_flx_add_shallow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxm_eval_powers_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxm_neg : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxm_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxm_to_flxxv :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val flxm_to_zxm : ('kind, 'structure) t -> ('kind, 'structure) t
val flxt_red : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxv_flc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxv_flv_multieval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxv_flx_fromdigits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxv_composedsum :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxv_prod : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flxv_red : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flxv_to_flm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val flxv_to_flxx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val flxv_to_zxv : ('kind, 'structure) t -> ('kind, 'structure) t
val flxv_to_zxv_inplace : ('kind, 'structure) t -> unit

val flxn_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) t

val flxn_div_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxn_exp :
  ('kind, 'structure) t -> Signed.long -> pari_ulong -> ('kind, 'structure) t

val flxn_expint :
  ('kind, 'structure) t -> Signed.long -> pari_ulong -> ('kind, 'structure) t

val flxn_inv :
  ('kind, 'structure) t -> Signed.long -> pari_ulong -> ('kind, 'structure) t

val flxn_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) t

val flxn_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxn_sqr :
  ('kind, 'structure) t -> Signed.long -> pari_ulong -> ('kind, 'structure) t

val flxn_sqr_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxn_red : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val flxq_autpow :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_autpow_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_autpowers :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_autsum :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_auttrace :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_auttrace_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_charpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_conjvec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_div_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_inv_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_invsafe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_invsafe_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_issquare :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> int

val flxq_is2npower :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  int

val flxq_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_lroot :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val flxq_lroot_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_lroot_fast :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val flxq_lroot_fast_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_matrix_pow :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_matrix_pow_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_minpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_minpoly_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_norm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flxq_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_pow_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_pow_init :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_pow_init_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_pow_table_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_pow_table :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_powu_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_powers_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_sqr_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_sqrt :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_sqrt_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_sqrtn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxq_trace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flxqc_flxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_flxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqv_dotproduct :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqv_dotproduct_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val rg_to_f2 : ('kind, 'structure) t -> pari_ulong
val rg_to_fl : ('kind, 'structure) t -> pari_ulong -> pari_ulong

val rg_to_flxq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val rgx_to_flx : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val rgxv_to_flxv : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val z_to_flx :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val zxv_to_flxv : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val zxt_to_flxt : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val gener_flxq :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val get_flxq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  pari_ulong ->
  bb_field Ctypes.structure Ctypes_static.ptr

val get_flxq_star :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  pari_ulong ->
  bb_group Ctypes.structure Ctypes_static.ptr

val monomial_flx :
  pari_ulong -> Signed.long -> Signed.long -> ('kind, 'structure) t

val random_flx :
  Signed.long -> Signed.long -> pari_ulong -> ('kind, 'structure) t

val zero_flxc : Signed.long -> Signed.long -> ('kind, 'structure) t

val zero_flxm :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val zlx_translate1 :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val zx_to_flx : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxx_fl_mul :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flxx_flx_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxx_flx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxx_flx_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxx_laplace : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxx_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxx_blocks :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val flxx_deriv : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flxx_double : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxx_invlaplace :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxx_neg : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxx_renormalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val flxx_shift :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val flxx_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxx_swap :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val flxx_to_flm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val flxx_to_flx : ('kind, 'structure) t -> ('kind, 'structure) t

val flxx_to_flxc :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val flxx_to_zxx : ('kind, 'structure) t -> ('kind, 'structure) t

val flxx_translate1 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val flxx_triple : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxxc_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxxc_to_zxxc : ('kind, 'structure) t -> ('kind, 'structure) t
val flxxm_to_zxxm : ('kind, 'structure) t -> ('kind, 'structure) t

val flxxv_to_flxm :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val flxxn_red : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val flxy_flx_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxy_flx_translate :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxy_flxqv_evalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxy_flxqv_evalx_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxy_flxq_evalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxy_flxq_evalx_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxy_evalx :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flxy_evalx_powers_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxy_evalx_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxyqq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqv_roots_to_pol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) t

val flxqxc_flxqxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxc_flxqxqv_eval_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxc_flxqxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_autpow :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_autpow_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_autsum :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_autsum_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_auttrace :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_auttrace_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_div_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_inv_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_invsafe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_invsafe_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_matrix_pow :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_minpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_minpoly_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_pow_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_powers_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_powu_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_sqr_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxv_prod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_flxqxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_flxqxqv_eval_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_flxqxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_flxqxq_eval_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_flxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_flxq_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_flxq_mul_to_monic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_flxq_mul_to_monic_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_newton :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_newton_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_composedsum :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_disc :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_div_by_x_x :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxqx_div_by_x_x_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxqx_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxqx_divrem_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxqx_dotproduct :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_extgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxqx_extgcd_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxqx_fromnewton :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_fromnewton_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_gcd_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_get_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_get_red_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_halfgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_halfgcd_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxqx_halfgcd_all_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxqx_halfgcd_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_invbarrett :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_invbarrett_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_normalize :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_normalize_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_powu_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_red_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_rem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_rem_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_resultant_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_safegcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_saferesultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_sqr_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxn_expint :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxn_expint_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxn_inv :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxn_inv_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxn_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxn_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxn_sqr :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxn_sqr_pre :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxy_degreex : ('kind, 'structure) t -> Signed.long

val flxy_eval_powers_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fly_to_flxy : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val kronecker_to_flxqx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val kronecker_to_flxqx_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val rgx_to_flxqx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val get_flxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val pol1_flxx : Signed.long -> Signed.long -> ('kind, 'structure) t
val polx_flxx : Signed.long -> Signed.long -> ('kind, 'structure) t

val random_flxqx :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val zlxx_translate1 :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val zxx_to_kronecker :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val flxq_ellcard :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_ellgens :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_ellgroup :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxq_elltwist :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val flxq_ellj :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_ellj_to_a4a6 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val flxqe_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_changepoint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_changepointinv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_dbl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_neg :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_tatepairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_weilpairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqe_weilpairing_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val zxx_to_flxx :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val zxxt_to_flxxt :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val zxxv_to_flxxv :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val get_flxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_flxqe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val random_flxqe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val polisclass : ('kind, 'structure) t -> Signed.long
val fl_elltrace : pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fl_elltrace_cm :
  Signed.long -> pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fp_ellcard :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_elldivpol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_ellgens :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_ellgroup :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fp_ellj :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_ellj_to_a4a6 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val fp_elljissupersingular :
  ('kind, 'structure) t -> ('kind, 'structure) t -> int

val fp_elltwist :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val fp_ffellcard :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_changepoint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_changepointinv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_dbl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_neg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpe_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_to_fpj : ('kind, 'structure) t -> ('kind, 'structure) t

val fpe_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpe_tatepairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpe_weilpairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpj_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpj_dbl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpj_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpj_neg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpj_to_fpe :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxq_ellcard :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_ellcard_supersingular :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_elldivpol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_ellgens :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_ellgroup :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpxq_ellj :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_elljissupersingular :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t -> int

val fpxq_elltwist :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val fpxqe_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_changepoint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_changepointinv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_dbl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_neg :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_tatepairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqe_weilpairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_elljissupersingular :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t -> int

val fq_ellcard_supersingular :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rge_to_fpe :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rge_to_fpxqe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val get_fpe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  bb_group Ctypes.structure Ctypes_static.ptr

val ellsupersingularj_fpxq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val elltrace_extension :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val random_fpe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val random_fpxqe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_issquare : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val fp_fpx_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_fpxq_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpv_fpm_polint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fpv_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpv_invvandermonde :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpv_polint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fpv_roots_to_pol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fpx_fp_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fp_add_shallow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fp_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fp_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fp_mul_to_monic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fp_mulspec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fpx_fp_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fp_sub_shallow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fpv_multieval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fpxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fpxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_fpxv_multirem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_frobenius :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_laplace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_newton :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_center_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_chinese_coprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_composedprod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_composedsum :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_convol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_deriv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_digits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_disc :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_div_by_x_x :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpx_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpx_divu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_dotproduct :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_extgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpx_extresultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpx_fromnewton :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_gcd_check :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_get_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_halve :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_halfgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_halfgcd_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpx_integ :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_invbarrett :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_invlaplace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_is_squarefree : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val fpx_matfrobenius :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_mulspec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val fpx_mulu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_neg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_normalize :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_rem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_rescale :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_valrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val fpxc_fpxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxc_fpxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxm_fpxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_autpow :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_autpowers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_autsum :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_auttrace :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_charpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_conjvec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_invsafe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_issquare :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t -> int

val fpxq_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_matrix_pow :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_minpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_norm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_sqrt :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_sqrtn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpxq_trace :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqc_to_mod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqm_autsum :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxt_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxv_fpx_fromdigits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxv_chinese :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpxv_composedsum :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxv_factorback :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fpxv_prod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxv_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxn_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxn_exp :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxn_expint :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxn_inv :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxn_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxn_sqr :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_issquare :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t -> int

val fq_ispower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fq_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqc_to_mod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_to_mod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqv_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val z_to_fpx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gener_fpxq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val gener_fpxq_local :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val get_fpxq_star :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpx :
  Signed.long -> Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2x_ddf : ('kind, 'structure) t -> ('kind, 'structure) t
val f2x_factor : ('kind, 'structure) t -> ('kind, 'structure) t
val f2x_factor_squarefree : ('kind, 'structure) t -> ('kind, 'structure) t
val f2x_is_irred : ('kind, 'structure) t -> int
val flx_ddf : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_ddf_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_is_irred : ('kind, 'structure) t -> pari_ulong -> int
val flx_is_totally_split : ('kind, 'structure) t -> pari_ulong -> int

val flx_ispower :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val flx_degfact : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flx_factor : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_factor_squarefree :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_factor_squarefree_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_nbfact : ('kind, 'structure) t -> pari_ulong -> Signed.long

val flx_nbfact_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> Signed.long

val flx_nbfact_frobenius :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> Signed.long

val flx_nbfact_frobenius_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  Signed.long

val flx_nbfact_by_degree :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_nbroots : ('kind, 'structure) t -> pari_ulong -> Signed.long
val flx_oneroot : ('kind, 'structure) t -> pari_ulong -> pari_ulong
val flx_oneroot_split : ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flx_oneroot_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> pari_ulong

val flx_oneroot_split_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> pari_ulong

val flx_roots : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_roots_pre :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flx_rootsff :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val fpx_ddf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_ddf_degree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fpx_degfact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_factor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_factor_squarefree :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_is_irred : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val fpx_is_totally_split : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val fpx_ispower :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val fpx_nbfact : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val fpx_nbfact_frobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fpx_nbroots : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val fpx_oneroot :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_oneroot_split :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_roots :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_roots_mult :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_rootsff :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_split_part :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_ddf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_degfact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_factor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_factor_squarefree :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqx_roots :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val flx_factorff_irred :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_ffintersect :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val flx_ffisom :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_ffisom_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_frobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_frobenius_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_ddf :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_ddf_degree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long

val flxqx_degfact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_factor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_factor_squarefree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_factor_squarefree_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_ispower :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val flxqx_is_squarefree :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> Signed.long

val flxqx_nbfact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> Signed.long

val flxqx_nbfact_frobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long

val flxqx_nbfact_by_degree :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_nbroots :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> Signed.long

val flxqx_roots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqxq_halffrobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val fpx_factorff :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_factorff_irred :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_ffintersect :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val fpx_ffisom :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_ffisom_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_frobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_ddf :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_ddf_degree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fpxqx_degfact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_factor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_factor_squarefree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_ispower :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val fpxqx_nbfact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fpxqx_nbfact_frobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fpxqx_nbroots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fpxqx_roots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_split_part :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_halffrobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_is_squarefree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fqx_ispower :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val fqx_nbfact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fqx_nbroots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val factorff :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val factormod0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val factormodddf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val factormodsqf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_parse_tp :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  int

val polrootsff :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val polrootsmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rootmod0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fpxqx_fpxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_fpxqxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_fpxqxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_digits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_disc :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_div_by_x_x :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpxqx_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpxqx_dotproduct :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_extgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpxqx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_get_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_halfgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_halfgcd_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpxqx_invbarrett :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_rem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_to_mod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_invsafe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_matrix_pow :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_minpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_autpow :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_autsum :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxq_auttrace :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxt_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxv_fpxqx_fromdigits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxv_prod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxv_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxn_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxn_exp :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxn_expint :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxn_inv :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxn_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqxn_sqr :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxx_fp_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxx_fpx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxx_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxx_deriv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxx_halve :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxx_integ :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxx_mulu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxx_neg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxx_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxx_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxy_fpxq_evalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxy_fpxqv_evalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxy_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxy_evalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxy_evaly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fpxyqq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxc_to_mod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxm_to_mod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val kronecker_to_fpxqx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val get_fpxqx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpxqx :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val flc_flv_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flc_to_mod : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flm_fl_add :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flm_fl_mul :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flm_fl_mul_inplace :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> unit

val flm_fl_mul_pre :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_fl_sub :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flm_flc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_flc_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_flc_mul_pre_flx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) t

val flm_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_center :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_mul_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_neg : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flm_powers :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flm_powu :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flm_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flm_to_mod : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flm_transpose : ('kind, 'structure) t -> ('kind, 'structure) t

val flv_fl_div :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flv_fl_div_inplace :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> unit

val flv_fl_mul :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flv_fl_mul_inplace :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> unit

val flv_fl_mul_part_inplace :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> Signed.long -> unit

val flv_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flv_add_inplace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> unit

val flv_center :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val flv_dotproduct :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flv_dotproduct_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val flv_neg : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flv_neg_inplace : ('kind, 'structure) t -> pari_ulong -> unit

val flv_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flv_sub_inplace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> unit

val flv_sum : ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flx_dotproduct :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flx_dotproduct_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val fp_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpc_fpv_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpc_fp_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpc_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpc_center_inplace :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val fpc_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpc_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpm_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_fp_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_fpc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_fpc_mul_fpx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fpm_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_center_inplace :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val fpm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpm_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpms_fpc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpms_fpcs_solve :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpms_fpcs_solve_safe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpms_leftkernel_elt :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpc_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpc_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpv_fpms_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpv_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpv_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpv_dotproduct :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpv_dotsquare :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpv_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpv_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpvv_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpx_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxc_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxm_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zabm_ker :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zabm_indexrank :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zabm_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zabm_inv_ratlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zabm_pseudoinv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zv_zms_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zpms_zpcs_solve :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gen_fpm_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gen_zpm_dixon_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gen_matid :
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val matid_flm : Signed.long -> ('kind, 'structure) t
val matid_f2xqm : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val matid_flxqm :
  Signed.long -> ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val random_flv : Signed.long -> pari_ulong -> ('kind, 'structure) t
val random_fpc : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val random_fpv : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val scalar_flm : Signed.long -> Signed.long -> ('kind, 'structure) t
val zcs_to_zc : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val zms_to_zm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zms_zc_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zmv_to_flmv : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flx_teichmuller :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val z2_sqrt : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zp_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zp_exp :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val zp_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zp_invlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zp_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val zp_sqrt :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zp_sqrtlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zp_sqrtnlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpm_invlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpx_frobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpx_zpxq_liftroot :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpx_zpxq_liftroot_ea :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val zpx_liftfact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpx_liftroot :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpx_liftroots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpx_roots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxq_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxq_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxq_invlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxq_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxq_sqrt :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxq_sqrtnlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxqm_prodfrobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxqx_digits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxqx_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zpxqx_liftfact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxqx_liftroot :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxqx_liftroot_vald :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxqx_liftroots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxqx_roots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpxqx_zpxqxq_liftroot :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zq_sqrtnlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zqx_zqxq_liftroot :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zqx_liftfact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zqx_liftroot :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zqx_roots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gen_zpm_dixon :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_zpm_newton :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_zpx_dixon :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_zpx_newton :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val polteichmuller :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val polhensellift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val quadratic_prec_mask : Signed.long -> pari_ulong
val qx_factor : ('kind, 'structure) t -> ('kind, 'structure) t
val zx_factor : ('kind, 'structure) t -> ('kind, 'structure) t
val zx_is_irred : ('kind, 'structure) t -> Signed.long

val zx_squff :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val polcyclofactors : ('kind, 'structure) t -> ('kind, 'structure) t
val poliscyclo : ('kind, 'structure) t -> Signed.long
val poliscycloprod : ('kind, 'structure) t -> Signed.long

val rg_rgc_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_rg_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_rg_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_rg_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_rg_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_rgm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_rgv_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_is_ei : ('kind, 'structure) t -> Signed.long
val rgc_neg : ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rg_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rg_add_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rg_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rg_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rg_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rg_sub_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rgc_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rgv_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_det_triangular : ('kind, 'structure) t -> ('kind, 'structure) t
val rgm_is_qm : ('kind, 'structure) t -> int
val rgm_is_zm : ('kind, 'structure) t -> int
val rgm_isdiagonal : ('kind, 'structure) t -> int
val rgm_isidentity : ('kind, 'structure) t -> int
val rgm_isscalar : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val rgm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_multosym :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_neg : ('kind, 'structure) t -> ('kind, 'structure) t
val rgm_powers : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgm_sqr : ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_sumcol : ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_transmul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_transmultosym :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgmrow_zc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgm_zc_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_zm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgmrow_rgc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgv_rgm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_rgc_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_rg_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_dotproduct :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_dotsquare : ('kind, 'structure) t -> ('kind, 'structure) t
val rgv_is_zmv : ('kind, 'structure) t -> int
val rgv_kill0 : ('kind, 'structure) t -> ('kind, 'structure) t
val rgv_neg : ('kind, 'structure) t -> ('kind, 'structure) t
val rgv_prod : ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_sum : ('kind, 'structure) t -> ('kind, 'structure) t
val rgv_sumpart : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgv_sumpart2 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rgv_zc_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_zm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rgm_eval :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rgmv_eval :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val isdiagonal : ('kind, 'structure) t -> int
val matid : Signed.long -> ('kind, 'structure) t
val scalarcol : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val scalarcol_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val scalarmat : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val scalarmat_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val scalarmat_s : Signed.long -> Signed.long -> ('kind, 'structure) t

val kronecker_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qx_zxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxq_charpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qxq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxq_to_mod_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxqc_to_mod_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxqm_to_mod_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxqv_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxqx_homogenous_evalpow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxqx_to_mod_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxqxv_to_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxv_qxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxy_qxq_evalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rg_rgx_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rg_get_0 : ('kind, 'structure) t -> ('kind, 'structure) t
val rg_get_1 : ('kind, 'structure) t -> ('kind, 'structure) t
val rg_to_rgc : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgm_to_rgxv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgm_to_rgxv_reverse :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgm_to_rgxx :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rgv_to_rgx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgv_to_rgm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgv_to_rgx_reverse :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_rgxq_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgx_rgxqv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgx_rgxn_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_rgxnv_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_rg_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rg_add_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rg_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rg_divexact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rg_eval_bk :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rg_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rg_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rgv_eval :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_addmulxn_shallow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_addmulxn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_addspec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_addspec_shallow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_affine :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgx_blocks :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rgx_deflate : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgx_deriv : ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_digits :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_div_by_x_x :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgx_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgx_divs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgx_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val rgx_even_odd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val rgx_homogenize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_homogenous_evalpow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgx_inflate : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_mul_i :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_mul_normalized :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_mul2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgx_mulxn : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_mulhigh_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_muls : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_mulspec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_neg : ('kind, 'structure) t -> ('kind, 'structure) t
val rgx_normalize : ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_pseudodivrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgx_pseudorem :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_recip : ('kind, 'structure) t -> ('kind, 'structure) t
val rgx_recip_i : ('kind, 'structure) t -> ('kind, 'structure) t
val rgx_recip_shallow : ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rem :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_renormalize_lg :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_rescale :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_rotate_shallow :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rgx_shift : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_shift_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_splitting :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_sqr : ('kind, 'structure) t -> ('kind, 'structure) t
val rgx_sqr_i : ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_sqrhigh_i :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_sqrspec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_to_rgc : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_translate :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_unscale :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxq_matrix_pow :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxq_norm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxq_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxq_trace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxqc_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxqm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxqm_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxqv_rgxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxqv_factorback :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxqv_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxqx_rgxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxqx_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgxqx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxqx_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxqx_pseudodivrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgxqx_pseudorem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxqx_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxqx_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxqx_translate :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxv_rgv_eval :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxv_prod : ('kind, 'structure) t -> ('kind, 'structure) t

val rgxv_rescale :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxv_to_rgm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgxv_unscale :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxx_to_rgm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgxy_degreex : ('kind, 'structure) t -> Signed.long
val rgxy_derivx : ('kind, 'structure) t -> ('kind, 'structure) t

val rgxy_swap :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rgxy_swapspec :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rgxn_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgxn_div_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgxn_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgxn_exp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgxn_expint : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgxn_inv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgxn_inv_i : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgxn_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgxn_powers :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rgxn_recip_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgxn_red_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgxn_reverse : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgxn_sqr : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgxn_sqrt : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgxnv_red_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgxn_powu :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val rgxn_powu_i :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val zx_translate :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_unscale2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zx_unscale :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_unscale_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_unscale_divpow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zx_z_unscale : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zxq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zxq_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zxqx_dvd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t -> int

val brent_kung_optpow : Signed.long -> Signed.long -> Signed.long -> Signed.long

val gen_bkeval :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  int ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_bkeval_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val get_rg_algebra : unit -> bb_algebra Ctypes.structure Ctypes_static.ptr
val rfrac_deflate_order : ('kind, 'structure) t -> Signed.long

val rfrac_deflate_max :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val rfrac_deflate :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zgc_g_mul_inplace : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val zgcs_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val g_zgc_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val g_zg_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zgc_g_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zgc_z_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zg_g_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zg_z_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zg_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zg_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zg_neg : ('kind, 'structure) t -> ('kind, 'structure) t
val zg_normalize : ('kind, 'structure) t -> ('kind, 'structure) t

val zg_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val flc_lincomb1_inplace :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  unit

val vecsmall_prod : ('kind, 'structure) t -> ('kind, 'structure) t

val qm_qc_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qm_det : ('kind, 'structure) t -> ('kind, 'structure) t
val qm_ker : ('kind, 'structure) t -> ('kind, 'structure) t

val qm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qm_sqr : ('kind, 'structure) t -> ('kind, 'structure) t
val rgm_check_zm : ('kind, 'structure) t -> string -> unit
val rgv_check_zv : ('kind, 'structure) t -> string -> unit

val z_zc_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zv_zc_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_q_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_z_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_z_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_z_divexact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_z_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_zv_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_divexactu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val zc_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_copy : ('kind, 'structure) t -> ('kind, 'structure) t

val zc_hnfremdiv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zc_is_ei : ('kind, 'structure) t -> Signed.long

val zc_lincomb :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zc_lincomb1_inplace :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val zc_lincomb1_inplace_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val zc_neg : ('kind, 'structure) t -> ('kind, 'structure) t

val zc_reducemodlll :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_reducemodmatrix :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_z_mul : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zm_q_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_z_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_z_divexact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_z_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_det_triangular : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_diag_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_divexactu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val zm_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val zm_equal0 : ('kind, 'structure) t -> int

val zm_hnfdivrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zm_ishnf : ('kind, 'structure) t -> int
val zm_isdiagonal : ('kind, 'structure) t -> int
val zm_isidentity : ('kind, 'structure) t -> int
val zm_isscalar : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val zm_max_lg : ('kind, 'structure) t -> Signed.long

val zm_mul_diag :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_multosym :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_neg : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_nm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_pow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_powu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val zm_reducemodlll :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_reducemodmatrix :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_sqr : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_supnorm : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_transmul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_transmultosym :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_togglesign : ('kind, 'structure) t -> unit

val zm_zm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zmrow_zc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zmrow_equal0 : ('kind, 'structure) t -> Signed.long -> int
val zv_abscmp : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val zv_cmp : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val zv_dotsquare : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_max_lg : ('kind, 'structure) t -> Signed.long
val zv_to_nv : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_togglesign : ('kind, 'structure) t -> unit
val gram_matrix : ('kind, 'structure) t -> ('kind, 'structure) t

val nm_z_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_to_flm : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val zm_to_zm : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_zc_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zmv_to_zmv : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_abs : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_content : ('kind, 'structure) t -> Signed.long

val zv_dotproduct :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val zv_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val zv_equal0 : ('kind, 'structure) t -> int
val zv_neg : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_neg_inplace : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_prod : ('kind, 'structure) t -> Signed.long
val zv_prod_z : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_sum : ('kind, 'structure) t -> Signed.long
val zv_sumpart : ('kind, 'structure) t -> Signed.long -> Signed.long
val zv_to_flv : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val zv_z_mul : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zv_zm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zvv_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val kronecker_to_zxqx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val kronecker_to_zxx :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val qx_zx_rem :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qx_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qx_sqr : ('kind, 'structure) t -> ('kind, 'structure) t

val qxqm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxqm_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxqx_qxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxqx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxqx_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxqx_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_check_qx : ('kind, 'structure) t -> string -> unit
val rgx_check_zx : ('kind, 'structure) t -> string -> unit
val rgx_check_zxx : ('kind, 'structure) t -> string -> unit

val z_zx_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_z_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_z_add_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_z_eval :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_z_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_z_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_affine :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zx_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val zx_deriv : ('kind, 'structure) t -> ('kind, 'structure) t

val zx_digits :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zxv_zx_fromdigits :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_div_by_x_1 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zx_divuexact : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val zx_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val zx_eval1 : ('kind, 'structure) t -> ('kind, 'structure) t
val zx_max_lg : ('kind, 'structure) t -> Signed.long
val zx_mod_xnm1 : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val zx_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_mulspec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val zx_mulu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val zx_neg : ('kind, 'structure) t -> ('kind, 'structure) t

val zx_rem :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_remi2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val zx_rescale2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zx_rescale :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_rescale_lt : ('kind, 'structure) t -> ('kind, 'structure) t
val zx_shifti : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val zx_sqr : ('kind, 'structure) t -> ('kind, 'structure) t
val zx_sqrspec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zx_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_val : ('kind, 'structure) t -> Signed.long

val zx_valrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val zxc_to_flxc :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val zxm_to_flxm :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val zxqm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zxqm_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zxqx_zxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zxqx_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zxqx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zxt_remi2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zxv_z_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zxv_dotproduct :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zxv_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val zxv_remi2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zxx_z_divexact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zxx_z_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zxx_z_add_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zxx_evalx0 : ('kind, 'structure) t -> ('kind, 'structure) t
val zxx_max_lg : ('kind, 'structure) t -> Signed.long

val zxx_mul_kronecker :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zxx_renormalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zxx_sqr_kronecker :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgxx_to_kronecker :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgxx_to_kronecker_spec :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val zxn_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zxn_sqr : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val scalar_zx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val scalar_zx_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zx_to_zx : ('kind, 'structure) t -> ('kind, 'structure) t

val zx_z_divexact :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val alg_centralproj :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val alg_complete :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val alg_csa_table :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val alg_cyclic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val alg_get_absdim : ('kind, 'structure) t -> Signed.long
val alg_get_abssplitting : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_aut : ('kind, 'structure) t -> ('kind, 'structure) t
val algaut : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_auts : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_b : ('kind, 'structure) t -> ('kind, 'structure) t
val algb : ('kind, 'structure) t -> ('kind, 'structure) t
val algcenter : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_center : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_char : ('kind, 'structure) t -> ('kind, 'structure) t
val algchar : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_degree : ('kind, 'structure) t -> Signed.long
val algdegree : ('kind, 'structure) t -> Signed.long
val alg_get_dim : ('kind, 'structure) t -> Signed.long
val algdim : ('kind, 'structure) t -> Signed.long -> Signed.long
val alg_get_hasse_f : ('kind, 'structure) t -> ('kind, 'structure) t
val alghassef : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_hasse_i : ('kind, 'structure) t -> ('kind, 'structure) t
val alghassei : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_invbasis : ('kind, 'structure) t -> ('kind, 'structure) t
val alginvbasis : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_multable : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_basis : ('kind, 'structure) t -> ('kind, 'structure) t
val algbasis : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_relmultable : ('kind, 'structure) t -> ('kind, 'structure) t
val algrelmultable : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_splitpol : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_splittingfield : ('kind, 'structure) t -> ('kind, 'structure) t
val algsplittingfield : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_splittingbasis : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_splittingbasisinv : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_splittingdata : ('kind, 'structure) t -> ('kind, 'structure) t
val algsplittingdata : ('kind, 'structure) t -> ('kind, 'structure) t
val alg_get_tracebasis : ('kind, 'structure) t -> ('kind, 'structure) t

val alg_hasse :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val alg_hilbert :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val alg_matrix :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val alg_model : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val alg_quotient :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val algradical : ('kind, 'structure) t -> ('kind, 'structure) t
val algsimpledec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val algsimpledec_ss :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val algsubalg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val alg_type : ('kind, 'structure) t -> Signed.long

val algadd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val algalgtobasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val algbasistoalg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val algcharpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val algdisc : ('kind, 'structure) t -> ('kind, 'structure) t

val algdivl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val algdivr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val alggroup :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val alggroupcenter :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val alghasse :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val alginit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val algindex : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val alginv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val algisassociative : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val algiscommutative : ('kind, 'structure) t -> int
val algisdivision : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val algisramified : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val algissemisimple : ('kind, 'structure) t -> int
val algissimple : ('kind, 'structure) t -> Signed.long -> int
val algissplit : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val algisdivl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val algisinv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val algmakeintegral :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val algmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val algmultable : ('kind, 'structure) t -> ('kind, 'structure) t
val alglat_get_primbasis : ('kind, 'structure) t -> ('kind, 'structure) t
val alglat_get_scalar : ('kind, 'structure) t -> ('kind, 'structure) t

val alglatadd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val alglatcontains :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val alglatelement :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val alglathnf :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val alglatindex :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val alglatinter :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val alglatmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val alglatlefttransporter :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val alglatrighttransporter :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val alglatsubset :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val algneg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val algnorm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val algpoleval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val algpow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val algprimesubalg : ('kind, 'structure) t -> ('kind, 'structure) t
val algramifiedplaces : ('kind, 'structure) t -> ('kind, 'structure) t

val algrandom :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val algsplit : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val algtomatrix :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val algsqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val algsub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val algtableinit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val algtensor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val algtrace :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val algtype : ('kind, 'structure) t -> Signed.long

val bnfgwgeneric :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val checkalg : ('kind, 'structure) t -> unit

val checkhasse :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val checklat : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val conjclasses_algcenter :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val galoischardet :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val galoischarpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val galoischartable : ('kind, 'structure) t -> ('kind, 'structure) t

val nfgrunwaldwang :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfgwkummer :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val f2ms_colelim : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val f2m_image : ('kind, 'structure) t -> ('kind, 'structure) t
val f2m_indexrank : ('kind, 'structure) t -> ('kind, 'structure) t
val f2m_suppl : ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqm_f2xqc_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqm_f2xqc_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqm_f2xqc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqm_deplin :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqm_det :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqm_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqm_ker :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqm_image :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqm_indexrank :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqm_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val f2xqm_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val f2xqm_rank : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val f2xqm_suppl :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val flm_image : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flm_indexrank : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flm_suppl : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val flxqm_flxqc_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_flxqc_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_flxqc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_deplin :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_det :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_ker :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_image :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_indexrank :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqm_rank :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> Signed.long

val flxqm_suppl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val fpm_fpc_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_fpc_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_deplin :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpm_det :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpm_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_image :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpm_indexrank :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpm_intersect :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_intersect_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpm_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_ker :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpm_rank : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val fpm_suppl :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fqm_fqc_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_fqc_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_fqc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_deplin :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_det :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_ker :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_image :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_indexrank :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqm_rank :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val fqm_suppl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qm_image_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val qm_image : ('kind, 'structure) t -> ('kind, 'structure) t

val qm_gauss :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qm_gauss_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qm_indexrank : ('kind, 'structure) t -> ('kind, 'structure) t
val qm_inv : ('kind, 'structure) t -> ('kind, 'structure) t
val qm_rank : ('kind, 'structure) t -> Signed.long

val rgm_fp_init :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgm_hadamard : ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rgc_invimage :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_diagonal : ('kind, 'structure) t -> ('kind, 'structure) t
val rgm_diagonal_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val rgm_inv : ('kind, 'structure) t -> ('kind, 'structure) t
val rgm_inv_upper : ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_invimage :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_solve :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_solve_realimag :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgms_structelim :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val zm_det : ('kind, 'structure) t -> ('kind, 'structure) t
val zm_detmult : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_gauss :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_ker : ('kind, 'structure) t -> ('kind, 'structure) t
val zm_imagecompl : ('kind, 'structure) t -> ('kind, 'structure) t
val zm_indeximage : ('kind, 'structure) t -> ('kind, 'structure) t
val zm_indexrank : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zm_inv_ratlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zm_pseudoinv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zm_rank : ('kind, 'structure) t -> Signed.long

val zlm_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val closemodinvertible :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val deplin : ('kind, 'structure) t -> ('kind, 'structure) t
val det : ('kind, 'structure) t -> ('kind, 'structure) t
val det0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val det2 : ('kind, 'structure) t -> ('kind, 'structure) t
val detint : ('kind, 'structure) t -> ('kind, 'structure) t
val eigen : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gauss :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gaussmodulo :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gaussmodulo2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gen_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_gauss_pivot :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_det :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_ker :
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_matcolinvimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_matcolmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_matinvimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_matmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val image : ('kind, 'structure) t -> ('kind, 'structure) t
val image2 : ('kind, 'structure) t -> ('kind, 'structure) t
val imagecompl : ('kind, 'structure) t -> ('kind, 'structure) t
val indexrank : ('kind, 'structure) t -> ('kind, 'structure) t

val inverseimage :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ker : ('kind, 'structure) t -> ('kind, 'structure) t

val mateigen :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val matimage0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val matker0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rank : ('kind, 'structure) t -> Signed.long

val reducemodinvertible :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val reducemodlll :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val split_realimag :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val suppl : ('kind, 'structure) t -> ('kind, 'structure) t
val flm_charpoly : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val flm_hess : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val fpm_charpoly :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpm_hess :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val frobeniusform :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val qm_minors_coprime :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qm_imz : ('kind, 'structure) t -> ('kind, 'structure) t

val qm_imz_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val qm_imz_hnf : ('kind, 'structure) t -> ('kind, 'structure) t

val qm_imz_hnfall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t

val qm_imq : ('kind, 'structure) t -> ('kind, 'structure) t

val qm_imq_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val qm_imq_hnf : ('kind, 'structure) t -> ('kind, 'structure) t

val qm_imq_hnfall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t

val qm_charpoly_zx : ('kind, 'structure) t -> ('kind, 'structure) t

val qm_charpoly_zx_bound :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zm_charpoly : ('kind, 'structure) t -> ('kind, 'structure) t
val adj : ('kind, 'structure) t -> ('kind, 'structure) t
val adjsafe : ('kind, 'structure) t -> ('kind, 'structure) t
val caract : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val caradj :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val carberkowitz : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val carhess : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val charpoly : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val charpoly0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val gnorm : ('kind, 'structure) t -> ('kind, 'structure) t
val gnorml1 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gnorml1_fake : ('kind, 'structure) t -> ('kind, 'structure) t

val gnormlp :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gnorml2 : ('kind, 'structure) t -> ('kind, 'structure) t
val gsupnorm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gsupnorm_aux :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  unit

val gtrace : ('kind, 'structure) t -> ('kind, 'structure) t
val hess : ('kind, 'structure) t -> ('kind, 'structure) t

val intersect :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val jacobi : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val matadjoint0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val matcompanion : ('kind, 'structure) t -> ('kind, 'structure) t

val matrixqz0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val minpoly : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val qfgaussred : ('kind, 'structure) t -> ('kind, 'structure) t
val qfgaussred_positive : ('kind, 'structure) t -> ('kind, 'structure) t
val qfsign : ('kind, 'structure) t -> ('kind, 'structure) t

val apply0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val diagonal : ('kind, 'structure) t -> ('kind, 'structure) t
val diagonal_shallow : ('kind, 'structure) t -> ('kind, 'structure) t

val extract0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fold0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val genapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val genfold :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val genindexselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val genselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gtomat : ('kind, 'structure) t -> ('kind, 'structure) t
val gtrans : ('kind, 'structure) t -> ('kind, 'structure) t

val matmuldiagonal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val matmultodiagonal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val matslice0 :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val parapply :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val parfor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parfor_init :
  parfor_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val parfor_next :
  parfor_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val parfor_stop : parfor_t Ctypes.structure Ctypes_static.ptr -> unit

val parforeach :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforeach_init :
  parforeach_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val parforeach_next :
  parforeach_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val parforeach_stop : parforeach_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprime_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val parforprime_next :
  parforprime_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val parforprime_stop : parforprime_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprimestep :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprimestep_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val parforvec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforvec_init :
  parforvec_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val parforvec_next :
  parforvec_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val parforvec_stop : parforvec_t Ctypes.structure Ctypes_static.ptr -> unit

val parselect :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val select0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val shallowextract :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val shallowmatextract :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val shallowtrans : ('kind, 'structure) t -> ('kind, 'structure) t

val vecapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val veccatapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val veccatselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val vecrange :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecrangess : Signed.long -> Signed.long -> ('kind, 'structure) t

val vecselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val vecselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val vecslice0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val vecsum : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_diagonal : ('kind, 'structure) t -> ('kind, 'structure) t
val addhelp : string -> string -> unit
val arity0 : ('kind, 'structure) t -> ('kind, 'structure) t
val alias0 : string -> string -> unit
val compile_str : string -> ('kind, 'structure) t
val delete_var : unit -> Signed.long
val fetch_user_var : string -> Signed.long
val fetch_var : unit -> Signed.long
val fetch_var_higher : unit -> Signed.long

val fetch_var_value :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val gp_embedded : string -> string
val gp_embedded_init : Signed.long -> Signed.long -> unit
val gp_read_str : string -> ('kind, 'structure) t
val gp_read_str_bitprec : string -> Signed.long -> ('kind, 'structure) t
val gp_read_str_prec : string -> Signed.long -> ('kind, 'structure) t

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
val readseq : string -> ('kind, 'structure) t

val safegel :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr

val safeel :
  ('kind, 'structure) t -> Signed.long -> Signed.long Ctypes_static.ptr

val safelistel :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr

val safegcoeff :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr

val strtoi : string -> ('kind, 'structure) t
val strtor : string -> Signed.long -> ('kind, 'structure) t
val varhigher : string -> Signed.long -> ('kind, 'structure) t
val varlower : string -> Signed.long -> ('kind, 'structure) t

val divisorslenstra :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val isprimeaprcl : ('kind, 'structure) t -> Signed.long

val qfb0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val check_quaddisc :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  string ->
  unit

val check_quaddisc_imag :
  ('kind, 'structure) t -> Signed.long Ctypes_static.ptr -> string -> unit

val check_quaddisc_real :
  ('kind, 'structure) t -> Signed.long Ctypes_static.ptr -> string -> unit

val cornacchia :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val cornacchia2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val cornacchia2_sqrt :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val nucomp :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nudupl :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nupow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val primeform :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val primeform_u : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val qfb_1 : ('kind, 'structure) t -> ('kind, 'structure) t

val qfbcomp :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfbcomp_i :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfbcompraw :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfbcompraw_i :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfbcornacchia :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfbpow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfbpow_i :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfbpowraw : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val qfbpows : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val qfbred : ('kind, 'structure) t -> ('kind, 'structure) t
val qfbred_i : ('kind, 'structure) t -> ('kind, 'structure) t

val qfbred0 :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfbredsl2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfbsolve :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qfbsqr : ('kind, 'structure) t -> ('kind, 'structure) t
val qfbsqr_i : ('kind, 'structure) t -> ('kind, 'structure) t

val qfisolvep :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfr3_comp :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val qfr3_compraw :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfr3_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val qfr3_red :
  ('kind, 'structure) t ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val qfr3_rho :
  ('kind, 'structure) t ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val qfr3_to_qfr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfr5_comp :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val qfr5_compraw :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfr5_dist :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qfr5_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val qfr5_red :
  ('kind, 'structure) t ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val qfr5_rho :
  ('kind, 'structure) t ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val qfr5_to_qfr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfr_data_init :
  ('kind, 'structure) t ->
  Signed.long ->
  qfr_data Ctypes.structure Ctypes_static.ptr ->
  unit

val qfr_to_qfr5 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val qfrsolvep :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val quadgen : ('kind, 'structure) t -> ('kind, 'structure) t
val quadgen0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val quadpoly : ('kind, 'structure) t -> ('kind, 'structure) t
val quadpoly_i : ('kind, 'structure) t -> ('kind, 'structure) t
val quadpoly0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
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
val fl_powers : pari_ulong -> Signed.long -> pari_ulong -> ('kind, 'structure) t

val fl_powers_pre :
  pari_ulong -> Signed.long -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

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

val fp_2gener : ('kind, 'structure) t -> ('kind, 'structure) t

val fp_2gener_i :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fp_factored_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_ispower :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t -> int

val fp_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_pow_init :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_pow_table :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_pows :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_sqrt :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fp_sqrt_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_sqrtn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fpv_prod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_zv_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_zv_mod_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val z_chinese :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val z_chinese_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val z_chinese_coprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val z_chinese_post :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val z_chinese_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val z_factor_listp :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_nv_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_nv_mod_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zv_allpnqn : ('kind, 'structure) t -> ('kind, 'structure) t

val zv_chinese :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zv_chinese_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zv_chinesetree :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zv_nv_mod_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zv_producttree : ('kind, 'structure) t -> ('kind, 'structure) t

val zx_nv_mod_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zxc_nv_mod_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zxm_nv_mod_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zxx_nv_mod_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zideallog :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bestappr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bestapprpade : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val chinese :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val chinese1 : ('kind, 'structure) t -> ('kind, 'structure) t
val chinese1_coprime_z : ('kind, 'structure) t -> ('kind, 'structure) t

val contfrac0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val contfracpnqn : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val fibo : Signed.long -> ('kind, 'structure) t
val gboundcf : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gcf : ('kind, 'structure) t -> ('kind, 'structure) t

val gcf2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val get_fp_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  bb_field Ctypes.structure Ctypes_static.ptr

val hilbert :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val hilbertii :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val istotient :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val krois : ('kind, 'structure) t -> Signed.long -> Signed.long
val kroiu : ('kind, 'structure) t -> pari_ulong -> Signed.long
val kronecker : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long
val krosi : Signed.long -> ('kind, 'structure) t -> Signed.long
val kross : Signed.long -> Signed.long -> Signed.long
val kroui : pari_ulong -> ('kind, 'structure) t -> Signed.long
val krouu : pari_ulong -> pari_ulong -> Signed.long

val lcmii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fp_invgen :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val logint0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val logintall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val mpfact : Signed.long -> ('kind, 'structure) t
val factorial_fl : Signed.long -> pari_ulong -> pari_ulong
val factorial_fp : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val muls_interval : Signed.long -> Signed.long -> ('kind, 'structure) t
val mulu_interval : pari_ulong -> pari_ulong -> ('kind, 'structure) t

val mulu_interval_step :
  pari_ulong -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val ncv_chinese_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ncv_chinese_center_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nmv_chinese_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val nmv_chinese_center_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nonsquare_fl : pari_ulong -> pari_ulong

val nxcv_chinese_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val nxcv_chinese_center_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nxmv_chinese_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val nxv_chinese_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val nxv_chinese_center_tree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zv_chinese_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val odd_prime_divisors : ('kind, 'structure) t -> ('kind, 'structure) t
val pgener_fl : pari_ulong -> pari_ulong
val pgener_fl_local : pari_ulong -> ('kind, 'structure) t -> pari_ulong
val pgener_fp : ('kind, 'structure) t -> ('kind, 'structure) t

val pgener_fp_local :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val pgener_zl : pari_ulong -> pari_ulong
val pgener_zp : ('kind, 'structure) t -> ('kind, 'structure) t
val pnqn : ('kind, 'structure) t -> ('kind, 'structure) t
val ramanujantau : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rootsof1_fl : pari_ulong -> pari_ulong -> pari_ulong

val rootsof1_fp :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rootsof1u_fp : pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t

val u_chinese_coprime :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val znlog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val znorder :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znprimroot : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgv_is_zvpos : ('kind, 'structure) t -> int
val rgv_is_zvnon0 : ('kind, 'structure) t -> int
val rgv_is_prv : ('kind, 'structure) t -> int
val z_issquarefree_fact : ('kind, 'structure) t -> Signed.long

val z_lsmoothen :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val z_smoothen :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val bigomega : ('kind, 'structure) t -> Signed.long
val bigomegau : pari_ulong -> Signed.long
val boundfact : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val check_arith_pos : ('kind, 'structure) t -> string -> ('kind, 'structure) t
val check_arith_non0 : ('kind, 'structure) t -> string -> ('kind, 'structure) t
val check_arith_all : ('kind, 'structure) t -> string -> ('kind, 'structure) t
val clean_z_factor : ('kind, 'structure) t -> ('kind, 'structure) t
val core : ('kind, 'structure) t -> ('kind, 'structure) t

val coredisc2_fact :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val coredisc2u_fact :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong

val corepartial : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val core0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val core2 : ('kind, 'structure) t -> ('kind, 'structure) t
val core2partial : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val coredisc : ('kind, 'structure) t -> ('kind, 'structure) t
val coredisc0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val coredisc2 : ('kind, 'structure) t -> ('kind, 'structure) t
val corediscs : Signed.long -> pari_ulong Ctypes_static.ptr -> Signed.long
val divisors : ('kind, 'structure) t -> ('kind, 'structure) t
val divisors_factored : ('kind, 'structure) t -> ('kind, 'structure) t
val divisors0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val divisorsu : pari_ulong -> ('kind, 'structure) t
val divisorsu_moebius : ('kind, 'structure) t -> ('kind, 'structure) t
val divisorsu_fact : ('kind, 'structure) t -> ('kind, 'structure) t
val divisorsu_fact_factored : ('kind, 'structure) t -> ('kind, 'structure) t
val eulerphi : ('kind, 'structure) t -> ('kind, 'structure) t
val eulerphiu : pari_ulong -> pari_ulong
val eulerphiu_fact : ('kind, 'structure) t -> pari_ulong
val factor_pn_1 : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val factor_pn_1_limit :
  ('kind, 'structure) t -> Signed.long -> pari_ulong -> ('kind, 'structure) t

val factoru_pow : pari_ulong -> ('kind, 'structure) t

val fuse_z_factor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val is_z_factor : ('kind, 'structure) t -> int
val is_z_factornon0 : ('kind, 'structure) t -> int
val is_z_factorpos : ('kind, 'structure) t -> int
val is_nf_factor : ('kind, 'structure) t -> int
val is_nf_extfactor : ('kind, 'structure) t -> int
val issquarefree : ('kind, 'structure) t -> Signed.long
val numdiv : ('kind, 'structure) t -> ('kind, 'structure) t
val numdivu : Signed.long -> Signed.long
val numdivu_fact : ('kind, 'structure) t -> Signed.long
val omega : ('kind, 'structure) t -> Signed.long
val omegau : pari_ulong -> Signed.long
val sumdiv : ('kind, 'structure) t -> ('kind, 'structure) t
val sumdivk : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val uissquarefree : pari_ulong -> Signed.long
val uissquarefree_fact : ('kind, 'structure) t -> Signed.long
val usumdiv_fact : ('kind, 'structure) t -> ('kind, 'structure) t
val usumdivk_fact : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val fpx_fpc_nfpoleval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val embed_t2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val embednorm_t2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val embed_norm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val check_zkmodule_i : ('kind, 'structure) t -> int
val check_zkmodule : ('kind, 'structure) t -> string -> unit
val checkbid : ('kind, 'structure) t -> unit
val checkbid_i : ('kind, 'structure) t -> ('kind, 'structure) t
val checkbnf : ('kind, 'structure) t -> ('kind, 'structure) t
val checkbnf_i : ('kind, 'structure) t -> ('kind, 'structure) t
val checkbnr : ('kind, 'structure) t -> unit
val checkbnr_i : ('kind, 'structure) t -> ('kind, 'structure) t
val checkabgrp : ('kind, 'structure) t -> unit
val checksqmat : ('kind, 'structure) t -> Signed.long -> unit
val checknf : ('kind, 'structure) t -> ('kind, 'structure) t
val checknf_i : ('kind, 'structure) t -> ('kind, 'structure) t

val checknfelt_mod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  string ->
  ('kind, 'structure) t

val checkprid : ('kind, 'structure) t -> unit
val checkprid_i : ('kind, 'structure) t -> int
val checkrnf : ('kind, 'structure) t -> unit
val checkrnf_i : ('kind, 'structure) t -> int

val factoredpolred :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val factoredpolred2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val galoisapply :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val get_bnf :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val get_bnfpol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val get_nf :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val get_nfpol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val get_prid : ('kind, 'structure) t -> ('kind, 'structure) t

val idealfrobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealfrobenius_aut :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealramfrobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealramfrobenius_aut :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealramgroups :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealramgroups_aut :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nf_get_allroots : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_prec : ('kind, 'structure) t -> Signed.long

val nfmaxord_to_nf :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfcertify : ('kind, 'structure) t -> ('kind, 'structure) t

val nfgaloismatrix :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfgaloismatrixapply :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfgaloispermtobasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfinit_basic :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t -> unit

val nfinit_complete :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val nfinit : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val nfinit0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val nfinitred : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val nfinitred2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val nfisincl :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfisincl0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfisisom :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfnewprec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val nfnewprec_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val nfpoleval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfsplitting :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfsplitting0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nftyp : ('kind, 'structure) t -> Signed.long
val polredord : ('kind, 'structure) t -> ('kind, 'structure) t
val polred : ('kind, 'structure) t -> ('kind, 'structure) t

val polred0 :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val polred2 : ('kind, 'structure) t -> ('kind, 'structure) t
val polredabs : ('kind, 'structure) t -> ('kind, 'structure) t
val polredabs0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val polredabs2 : ('kind, 'structure) t -> ('kind, 'structure) t
val polredabsall : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val polredbest : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val poltomonic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rnfpolredabs :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfpolredbest :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val smallpolred : ('kind, 'structure) t -> ('kind, 'structure) t
val smallpolred2 : ('kind, 'structure) t -> ('kind, 'structure) t
val tschirnhaus : ('kind, 'structure) t -> ('kind, 'structure) t

val zx_q_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_q_normalize :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zx_z_normalize :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zx_to_monic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zx_primitive_to_monic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zxx_q_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fq_to_nf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fqm_to_nfm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fqv_to_nfv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fqx_to_nfx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rg_nffix :
  string ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int ->
  ('kind, 'structure) t

val rgv_nffix :
  string ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int ->
  ('kind, 'structure) t

val rgx_nffix :
  string ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int ->
  ('kind, 'structure) t

val zx_composedsum :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_compositum :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val zpx_disc_val : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val zpx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zpx_monic_factor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpx_primedec :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zpx_reduced_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zpx_reduced_resultant_fast :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zpx_resultant_val :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long

val checkmodpr : ('kind, 'structure) t -> unit

val compositum :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val compositum2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfdisc : ('kind, 'structure) t -> ('kind, 'structure) t
val get_modpr : ('kind, 'structure) t -> ('kind, 'structure) t

val indexpartial :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val modpr_genfq : ('kind, 'structure) t -> ('kind, 'structure) t

val nf_to_fq_init :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val nf_to_fq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfm_to_fqm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfv_to_fqv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfx_to_fqx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfx_to_monic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val nfbasis :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val nfcompositum :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfdiscfactors : ('kind, 'structure) t -> ('kind, 'structure) t

val nfmaxord :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val nfmodpr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfmodprinit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfmodprinit0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfmodprlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfreducemodpr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val polcompositum0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealprimedec :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealprimedec_galois :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealprimedec_degrees :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealprimedec_kummer :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealprimedec_limit_f :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealprimedec_limit_norm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val poldiscfactors :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rnfbasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfdedekind :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfdet :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfdisc_factored :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rnfdiscf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfequation :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfequation0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfequation2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_pv_to_prv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_rnfeq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_rnfeqsimple :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfequationall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rnfhnfbasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfisfree : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val rnflllgram :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfpolred :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfpseudobasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfsimplifybasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfsteinitz :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val factorial_lval : pari_ulong -> pari_ulong -> Signed.long

val zk_to_fq_init :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zk_to_fq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxqv_to_fpm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zkmodprinit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealstar :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealstarprk :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rgc_to_nfc :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_rgx_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_to_nfm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_to_nfx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_nfval : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val zc_nfvalrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val zc_prdvd : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val zm_zx_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zv_snf_gcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val algtobasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val basistoalg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ei_multable : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val get_nf_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  bb_field Ctypes.structure Ctypes_static.ptr

val famat_nfvalrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val gpnfvalrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val idealfactorback :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int ->
  ('kind, 'structure) t

val ideallist : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ideallist0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val gideallist :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ideallistarch :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ideallog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ideallogmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ideallog_units :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ideallog_units0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealprincipalunits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealstar0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealstarmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val indices_to_vec01 :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val matalgtobasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val matbasistoalg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val multable :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_to_scalar_or_alg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_to_scalar_or_basis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_cxlog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfv_cxlog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfadd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfchecksigns :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t -> int

val nfdiv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfdiveuc :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfdivrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfembed :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfeltembed :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfeltembed_i :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfeltsign :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nffactorback :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfinv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfinvmodideal :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfissquare :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val nfispower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val nflogembed :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t

val nfm_det :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfm_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfm_ker :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfm_nfc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfmuli :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfnorm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfpolsturm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfpow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfpow_u :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val nfpowmodideal :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfsign :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfsign_arch :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfsign_from_logarch :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfsqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfsqri :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfsub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nftrace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val nfvalrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val polmod_nffix :
  string ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int ->
  ('kind, 'structure) t

val polmod_nffix2 :
  string ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int ->
  ('kind, 'structure) t

val pr_basis_perm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val pr_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val rnfalgtobasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfbasistoalg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfeltnorm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfelttrace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val set_sign_mod_divisor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val tablemul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val tablemul_ei :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val tablemul_ei_ej :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val tablemulvec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val tablesqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vec01_to_indices : ('kind, 'structure) t -> ('kind, 'structure) t
val vecsmall01_to_indices : ('kind, 'structure) t -> ('kind, 'structure) t

val zk_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zk_multable :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zk_scalar_or_multable :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zkchinese :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zkchinese1 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zkchineseinit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zkmultable_capz : ('kind, 'structure) t -> ('kind, 'structure) t
val zkmultable_inv : ('kind, 'structure) t -> ('kind, 'structure) t

val fl_invgen :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val rm_round_maxrank : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_famat_limit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zv_cba : ('kind, 'structure) t -> ('kind, 'structure) t

val zv_cba_extend :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_cba :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_ppgle :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_ppio :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_ppo :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famatv_factorback :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famatv_zv_factorback :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_z_gcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_div_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_idealfactor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_inv : ('kind, 'structure) t -> ('kind, 'structure) t
val famat_inv_shallow : ('kind, 'structure) t -> ('kind, 'structure) t

val famat_makecoprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val famat_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_mul_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_mulpow_shallow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val famat_mulpows_shallow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val famat_pow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_pow_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_pows_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val famat_reduce : ('kind, 'structure) t -> ('kind, 'structure) t
val famat_remove_trivial : ('kind, 'structure) t -> ('kind, 'structure) t
val famat_sqr : ('kind, 'structure) t -> ('kind, 'structure) t

val famat_to_nf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val famat_to_nf_moddivisor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val famat_to_nf_modideal_coprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val famatsmall_reduce : ('kind, 'structure) t -> ('kind, 'structure) t

val gpidealfactor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gpidealval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealhnf_z_factor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val idealhnf_z_factor_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val idealhnf_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealhnf_inv_z :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealhnf_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealadd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealaddmultoone :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealaddtoone :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealaddtoone0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealaddtoone_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealaddtoone_raw :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealappr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealappr0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealapprfact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealchinese :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealcoprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealcoprimefact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealdiv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealdiv0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealdivexact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealdivpowprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealdown :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealfactor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealfactor_limit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val idealfactor_partial :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealhnf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealhnf0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealhnf_principal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealhnf_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealhnf_two :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealintersect :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealinv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealismaximal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealispower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val idealmin :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealmul0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealmulpowprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealmulred :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealnorm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealnumden :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealpow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealpow0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealpowred :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealpows :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val idealprod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealprodprime :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealprodval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val idealpseudomin :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealpseudomin_nonscalar :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealpseudominvec :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealpseudored :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealred0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealred_elt :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealredmodpower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val idealsqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealtwoelt :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealtwoelt0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealtwoelt2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val idealtyp :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val idealval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val isideal : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long
val matreduce : ('kind, 'structure) t -> ('kind, 'structure) t

val nfc_multable_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfc_nf_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nf_get_gtwist :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_get_gtwist1 :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val nf_to_fp_coprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfdetint :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfdivmodpr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfhnf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfhnf0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfhnfmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfkermodpr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfmulmodpr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfpowmodpr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfreduce :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfsnf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfsnf0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfsolvemodpr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val prv_lcm_capz : ('kind, 'structure) t -> ('kind, 'structure) t
val prv_primes : ('kind, 'structure) t -> ('kind, 'structure) t

val pr_hnf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val pr_inv : ('kind, 'structure) t -> ('kind, 'structure) t
val pr_inv_p : ('kind, 'structure) t -> ('kind, 'structure) t

val pr_uniformizer :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val sunits_makecoprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val to_famat :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val to_famat_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val u_ppo : pari_ulong -> pari_ulong -> pari_ulong

val vecdiv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecinv : ('kind, 'structure) t -> ('kind, 'structure) t

val vecmul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecpow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecsqr : ('kind, 'structure) t -> ('kind, 'structure) t

val zkc_multable_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val eltreltoabs :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val eltabstorel :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val eltabstorel_lift :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_nfzk :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnf_build_nfabs :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rnf_zkabs : ('kind, 'structure) t -> ('kind, 'structure) t

val nfeltup :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rnfcomplete : ('kind, 'structure) t -> unit

val rnfeltabstorel :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfeltdown :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfeltdown0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfeltreltoabs :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfeltup :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfeltup0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfidealabstorel :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealdown :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealfactor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealhnf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rnfidealnormabs :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealnormrel :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealprimedec :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealreltoabs :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealreltoabs0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfidealtwoelement :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealup :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfidealup0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfinit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfinit0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val get_arith_zzm : ('kind, 'structure) t -> ('kind, 'structure) t
val get_arith_z : ('kind, 'structure) t -> ('kind, 'structure) t

val gen_ph_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_shanks_init :
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_shanks :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_shanks_sqrtn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_gener :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_ellgens :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_ellgroup :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_factored_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_order :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_select_order :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_plog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val gen_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_pow_fold :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_pow_fold_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_pow_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_pow_init :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_pow_table :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  int ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t) Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_powu_fold :
  ('kind, 'structure) t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_powu_fold_i :
  ('kind, 'structure) t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_powu_i :
  ('kind, 'structure) t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_product :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val matdetmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val matimagemod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val matinvmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val matkermod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val matsolvemod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bernfrac : Signed.long -> ('kind, 'structure) t
val bernpol : Signed.long -> Signed.long -> ('kind, 'structure) t
val bernreal : Signed.long -> Signed.long -> ('kind, 'structure) t
val bernvec : Signed.long -> ('kind, 'structure) t
val constbern : Signed.long -> unit
val eulerfrac : Signed.long -> ('kind, 'structure) t
val eulerpol : Signed.long -> Signed.long -> ('kind, 'structure) t
val eulerreal : Signed.long -> Signed.long -> ('kind, 'structure) t
val eulervec : Signed.long -> ('kind, 'structure) t
val harmonic : pari_ulong -> ('kind, 'structure) t
val harmonic0 : pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t

val qr_init :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  int

val r_from_qr : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgm_babai :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_qr_init :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  int

val rgm_gram_schmidt :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val algdep : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val algdep0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val bestapprnf :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val forqfvec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  float ->
  Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val forqfvec0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val forqfvec1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val gaussred_from_qr :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val lindep : ('kind, 'structure) t -> ('kind, 'structure) t
val lindep_xadic : ('kind, 'structure) t -> ('kind, 'structure) t
val lindep_bit : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val lindep_padic : ('kind, 'structure) t -> ('kind, 'structure) t
val lindep0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val lindep2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val lindepfull_bit :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mathouseholder :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val matqr :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val minim :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val minim_raw :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val minim_zm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val minim2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfminim0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val qfperfection : ('kind, 'structure) t -> ('kind, 'structure) t

val qfrep0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val seralgdep :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val serdiffdep :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val vandermondeinverse :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val vandermondeinverseinit : ('kind, 'structure) t -> ('kind, 'structure) t

val zncoppersmith :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxq_reverse :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vec_equiv : ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_polint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val vec_reduce :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgxq_reverse :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zc_union_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zv_indexsort : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_sort : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_sort_inplace : ('kind, 'structure) t -> unit
val zv_sort_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_sort_uniq : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_sort_uniq_shallow : ('kind, 'structure) t -> ('kind, 'structure) t

val zv_union_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val binomial : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val binomial0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val binomialuu : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val cmp_flx : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val cmp_rgx : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val cmp_nodata :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int

val cmp_prime_ideal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val cmp_prime_over_p : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val cmp_universal : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val convol :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gen_cmp_rgx :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int

val polcyclo : Signed.long -> Signed.long -> ('kind, 'structure) t

val polcyclo_eval :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val dirdiv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val dirmul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val eulerianpol : Signed.long -> Signed.long -> ('kind, 'structure) t

val gprec_wensure :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gen_indexsort :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_indexsort_uniq :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_search :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  Signed.long

val gen_setminus :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  (('kind, 'structure) t -> ('kind, 'structure) t -> int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_sort :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_sort_inplace :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val gen_sort_shallow :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_sort_uniq :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val getstack : unit -> Signed.long
val gettime : unit -> Signed.long
val getabstime : unit -> Signed.long
val getwalltime : unit -> ('kind, 'structure) t
val gprec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gprec_wtrunc : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gprec_w : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtoset : ('kind, 'structure) t -> ('kind, 'structure) t
val indexlexsort : ('kind, 'structure) t -> ('kind, 'structure) t
val indexsort : ('kind, 'structure) t -> ('kind, 'structure) t

val indexvecsort :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val laplace : ('kind, 'structure) t -> ('kind, 'structure) t
val lexsort : ('kind, 'structure) t -> ('kind, 'structure) t
val mathilbert : Signed.long -> ('kind, 'structure) t
val matqpascal : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val merge_factor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val merge_sort_uniq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val modreverse : ('kind, 'structure) t -> ('kind, 'structure) t
val polhermite : Signed.long -> Signed.long -> ('kind, 'structure) t

val polhermite_eval0 :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val polhermite_eval :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val pollaguerre :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val pollaguerre_eval :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val pollaguerre_eval0 :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val pollegendre : Signed.long -> Signed.long -> ('kind, 'structure) t
val pollegendre_reduced : Signed.long -> Signed.long -> ('kind, 'structure) t

val pollegendre_eval :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val pollegendre_eval0 :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val polint :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val polint_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val polintspec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val polchebyshev :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val polchebyshev_eval :
  Signed.long -> Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val polchebyshev1 : Signed.long -> Signed.long -> ('kind, 'structure) t
val polchebyshev2 : Signed.long -> Signed.long -> ('kind, 'structure) t
val polrecip : ('kind, 'structure) t -> ('kind, 'structure) t

val setbinop :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val setdelta :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val setintersect :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val setisset : ('kind, 'structure) t -> Signed.long

val setminus :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val setsearch :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long -> Signed.long

val setunion :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val setunion_i :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val sort : ('kind, 'structure) t -> ('kind, 'structure) t

val sort_factor :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val stirling :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val stirling1 : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val stirling2 : pari_ulong -> pari_ulong -> ('kind, 'structure) t

val tablesearch :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  (('kind, 'structure) t -> ('kind, 'structure) t -> int)
  Ctypes_static.static_funptr ->
  Signed.long

val vecbinomial : Signed.long -> ('kind, 'structure) t

val vecsearch :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val vecsort :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecsort0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zv_search : ('kind, 'structure) t -> Signed.long -> Signed.long
val bits_to_int : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val bits_to_u : ('kind, 'structure) t -> Signed.long -> pari_ulong
val binaire : ('kind, 'structure) t -> ('kind, 'structure) t
val binary_2k : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val binary_2k_nv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val binary_zv : ('kind, 'structure) t -> ('kind, 'structure) t
val bittest : ('kind, 'structure) t -> Signed.long -> Signed.long

val fromdigits_2k :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gbitand :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gbitneg : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gbitnegimply :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gbitor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gbittest : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gbitxor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hammingl : pari_ulong -> Signed.long
val hammingweight : ('kind, 'structure) t -> Signed.long

val ibitand :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ibitnegimply :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ibitor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ibitxor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nv_fromdigits_2k :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val bnflogef :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnflog :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnflogdegree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfislocalpower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val rnfislocalcyclo : ('kind, 'structure) t -> Signed.long

val bnfisunit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnfissunit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnfsunit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnfunits :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnfisunit0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val sunits_mod_units :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val buchquad :
  ('kind, 'structure) t ->
  float ->
  float ->
  Signed.long ->
  ('kind, 'structure) t

val quadclassno : ('kind, 'structure) t -> ('kind, 'structure) t
val quadclassnos : Signed.long -> Signed.long

val quadclassunit0 :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val buchall :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val buchall_param :
  ('kind, 'structure) t ->
  float ->
  float ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val bnf_build_cheapfu : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_build_cycgen : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_build_matalpha : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_build_units : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_compactfu : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_compactfu_mat : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_has_fu : ('kind, 'structure) t -> ('kind, 'structure) t

val bnfinit0 :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnfisprincipal0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnfnewprec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val bnfnewprec_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val bnftestprimes : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val bnrnewprec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val bnrnewprec_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val isprincipalfact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val isprincipalfact_or_fail :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val isprincipal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_cxlog_normalize :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val nfcyclotomicunits :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfsign_units :
  ('kind, 'structure) t -> ('kind, 'structure) t -> int -> ('kind, 'structure) t

val nfsign_tu :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfsign_fu :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val signunits : ('kind, 'structure) t -> ('kind, 'structure) t
val hermite_bound : Signed.long -> Signed.long -> ('kind, 'structure) t

val bnr_subgroup_sanitize :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val bnr_char_sanitize :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val abc_to_bnr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int ->
  ('kind, 'structure) t

val buchray :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val buchraymod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnrautmatrix :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnr_subgroup_check :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val bnrchar :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnrchar_primitive :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnrclassno :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnrclassno0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnrclassnolist :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnrchar_primitive_raw :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnrconductor_factored :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnrconductor_raw :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnrconductormod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnrconductor0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnrconductor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnrconductor_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnrconductorofchar :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnrdisc0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnrdisc :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnrdisclist0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnrgaloismatrix :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnrgaloisapply :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnrinit0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnrinitmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val bnrisconductor0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val bnrisconductor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val bnrisgalois :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val bnrisprincipalmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnrisprincipal :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnrmap :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnrsurjection :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnfnarrow : ('kind, 'structure) t -> ('kind, 'structure) t
val bnfcertify : ('kind, 'structure) t -> Signed.long
val bnfcertify0 : ('kind, 'structure) t -> Signed.long -> Signed.long

val bnrcompositum :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val decodemodule :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val discrayabslist :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val discrayabslistarch :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val idealmoddivisor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val isprincipalray :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val isprincipalraygen :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nf_deg1_prime : ('kind, 'structure) t -> ('kind, 'structure) t

val nfarchstar :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rnfconductor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfconductor0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfnormgroup :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val subgrouplist0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnfisnorm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfisnorm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rnfisnorminit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> int -> ('kind, 'structure) t

val coprimes_zv : pari_ulong -> ('kind, 'structure) t
val char_check : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val charconj :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val charconj0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val chardiv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val chardiv0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val chareval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val chargalois :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val charker :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val charker0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val charmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val charmul0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val charorder :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val charorder0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val charpow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val charpow0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val char_denormalize :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val char_normalize :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val char_simplify :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val checkznstar_i : ('kind, 'structure) t -> int
val cyc_normalize : ('kind, 'structure) t -> ('kind, 'structure) t

val ncharvecexpo :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znchar : ('kind, 'structure) t -> ('kind, 'structure) t

val znchar_quad :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zncharcheck : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val zncharconductor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zncharconj :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znchardecompose :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val znchardiv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val znchareval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zncharinduce :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zncharisodd : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val zncharker :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zncharmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zncharorder :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zncharpow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val znchartokronecker :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val znchartoprimitive :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znconrey_check : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val znconrey_normalized :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znconreychar :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znconreyfromchar_normalized :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znconreyconductor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val znconreyexp :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znconreyfromchar :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znconreylog :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znconreylog_normalize :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val znlog0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zv_cyc_minimal :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val zv_cyc_minimize :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val closure_deriv : ('kind, 'structure) t -> ('kind, 'structure) t

val closure_derivn :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val localvars_find :
  ('kind, 'structure) t ->
  entree Ctypes.structure Ctypes_static.ptr ->
  Signed.long

val localvars_read_str :
  string -> ('kind, 'structure) t -> ('kind, 'structure) t

val snm_closure :
  entree Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val strtoclosure : string -> Signed.long -> ('kind, 'structure) t
val strtofunction : string -> ('kind, 'structure) t

val gconcat :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gconcat1 : ('kind, 'structure) t -> ('kind, 'structure) t
val matconcat : ('kind, 'structure) t -> ('kind, 'structure) t

val shallowconcat :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val shallowconcat1 : ('kind, 'structure) t -> ('kind, 'structure) t
val shallowmatconcat : ('kind, 'structure) t -> ('kind, 'structure) t

val vconcat :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val default0 : string -> string -> ('kind, 'structure) t
val getrealprecision : unit -> Signed.long
val pari_is_default : string -> entree Ctypes.structure Ctypes_static.ptr
val sd_texstyle : string -> Signed.long -> ('kind, 'structure) t
val sd_colors : string -> Signed.long -> ('kind, 'structure) t
val sd_compatible : string -> Signed.long -> ('kind, 'structure) t
val sd_datadir : string -> Signed.long -> ('kind, 'structure) t
val sd_debug : string -> Signed.long -> ('kind, 'structure) t
val sd_debugfiles : string -> Signed.long -> ('kind, 'structure) t
val sd_debugmem : string -> Signed.long -> ('kind, 'structure) t
val sd_factor_add_primes : string -> Signed.long -> ('kind, 'structure) t
val sd_factor_proven : string -> Signed.long -> ('kind, 'structure) t
val sd_format : string -> Signed.long -> ('kind, 'structure) t
val sd_histsize : string -> Signed.long -> ('kind, 'structure) t
val sd_log : string -> Signed.long -> ('kind, 'structure) t
val sd_logfile : string -> Signed.long -> ('kind, 'structure) t
val sd_nbthreads : string -> Signed.long -> ('kind, 'structure) t
val sd_new_galois_format : string -> Signed.long -> ('kind, 'structure) t
val sd_output : string -> Signed.long -> ('kind, 'structure) t
val sd_parisize : string -> Signed.long -> ('kind, 'structure) t
val sd_parisizemax : string -> Signed.long -> ('kind, 'structure) t
val sd_path : string -> Signed.long -> ('kind, 'structure) t
val sd_plothsizes : string -> Signed.long -> ('kind, 'structure) t
val sd_prettyprinter : string -> Signed.long -> ('kind, 'structure) t
val sd_primelimit : string -> Signed.long -> ('kind, 'structure) t
val sd_realbitprecision : string -> Signed.long -> ('kind, 'structure) t
val sd_realprecision : string -> Signed.long -> ('kind, 'structure) t
val sd_secure : string -> Signed.long -> ('kind, 'structure) t
val sd_seriesprecision : string -> Signed.long -> ('kind, 'structure) t
val sd_simplify : string -> Signed.long -> ('kind, 'structure) t
val sd_sopath : string -> int -> ('kind, 'structure) t
val sd_strictargs : string -> Signed.long -> ('kind, 'structure) t
val sd_strictmatch : string -> Signed.long -> ('kind, 'structure) t

val sd_string :
  string ->
  Signed.long ->
  string ->
  string Ctypes_static.ptr ->
  ('kind, 'structure) t

val sd_threadsize : string -> Signed.long -> ('kind, 'structure) t
val sd_threadsizemax : string -> Signed.long -> ('kind, 'structure) t

val sd_intarray :
  string ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  string ->
  ('kind, 'structure) t

val sd_toggle :
  string ->
  Signed.long ->
  string ->
  int Ctypes_static.ptr ->
  ('kind, 'structure) t

val sd_ulong :
  string ->
  Signed.long ->
  string ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  pari_ulong ->
  string Ctypes_static.ptr ->
  ('kind, 'structure) t

val setdefault : string -> string -> Signed.long -> ('kind, 'structure) t

val setrealprecision :
  Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val digits :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fromdigits :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fromdigitsu :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gen_digits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val gen_fromdigits :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t

val sumdigits : ('kind, 'structure) t -> ('kind, 'structure) t

val sumdigits0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val sumdigitsu : pari_ulong -> pari_ulong
val ecpp : ('kind, 'structure) t -> ('kind, 'structure) t
val ecpp0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ecppexport : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ecppisvalid : ('kind, 'structure) t -> Signed.long
val isprimeecpp : ('kind, 'structure) t -> Signed.long
val sd_breakloop : string -> Signed.long -> ('kind, 'structure) t
val sd_echo : string -> Signed.long -> ('kind, 'structure) t
val sd_graphcolormap : string -> Signed.long -> ('kind, 'structure) t
val sd_graphcolors : string -> Signed.long -> ('kind, 'structure) t
val sd_help : string -> Signed.long -> ('kind, 'structure) t
val sd_histfile : string -> Signed.long -> ('kind, 'structure) t
val sd_lines : string -> Signed.long -> ('kind, 'structure) t
val sd_linewrap : string -> Signed.long -> ('kind, 'structure) t
val sd_prompt : string -> Signed.long -> ('kind, 'structure) t
val sd_prompt_cont : string -> Signed.long -> ('kind, 'structure) t
val sd_psfile : string -> Signed.long -> ('kind, 'structure) t
val sd_readline : string -> Signed.long -> ('kind, 'structure) t
val sd_recover : string -> Signed.long -> ('kind, 'structure) t
val sd_timer : string -> Signed.long -> ('kind, 'structure) t
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
val gp_alarm : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val gp_input : unit -> ('kind, 'structure) t
val gp_allocatemem : ('kind, 'structure) t -> unit
val gp_handle_exception : Signed.long -> int
val gp_alarm_handler : int -> unit
val gp_sigint_fun : unit -> unit
val gp_help : string -> Signed.long -> unit
val gp_echo_and_log : string -> string -> unit
val print_fun_list : string Ctypes_static.ptr -> Signed.long -> unit
val strtime : Signed.long -> ('kind, 'structure) t

val direuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val dirpowers :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val dirpowerssum :
  pari_ulong ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val dirpowerssumfun :
  pari_ulong ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> pari_ulong -> Signed.long -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val vecpowuu : Signed.long -> pari_ulong -> ('kind, 'structure) t

val vecpowug :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellanalyticrank :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellanalyticrank_bitprec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellanal_globalred_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ellheegner : ('kind, 'structure) t -> ('kind, 'structure) t

val elll1 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val elll1_bitprec :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellconvertname : ('kind, 'structure) t -> ('kind, 'structure) t
val elldatagenerators : ('kind, 'structure) t -> ('kind, 'structure) t
val ellidentify : ('kind, 'structure) t -> ('kind, 'structure) t
val ellsearch : ('kind, 'structure) t -> ('kind, 'structure) t
val ellsearchcurve : ('kind, 'structure) t -> ('kind, 'structure) t

val forell :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  unit

val ellfromeqn : ('kind, 'structure) t -> ('kind, 'structure) t

val akell :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bilhell :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val checkell : ('kind, 'structure) t -> unit
val checkell_fq : ('kind, 'structure) t -> unit
val checkell_q : ('kind, 'structure) t -> unit
val checkell_qp : ('kind, 'structure) t -> unit
val checkellisog : ('kind, 'structure) t -> unit
val checkellpt : ('kind, 'structure) t -> unit
val checkell5 : ('kind, 'structure) t -> unit

val cxredsl2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val cxredsl2_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ec_2divpol_evalx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ec_3divpol_evalx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ec_bmodel : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ec_f_evalx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ec_h_evalx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ec_dfdx_evalq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ec_dfdy_evalq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ec_dmfdy_evalq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ec_half_deriv_2divpol :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ec_half_deriv_2divpol_evalx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ec_phi2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ell_is_integral : ('kind, 'structure) t -> int
val ellq_get_cm : ('kind, 'structure) t -> Signed.long
val ellq_get_n : ('kind, 'structure) t -> ('kind, 'structure) t

val ellq_get_nfa :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val ellqp_tate_uniformization :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellqp_agm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellqp_u : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellqp_u2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellqp_q : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellqp_ab : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellqp_l : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellqp_root : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellqtwist_bsdperiod :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellr_area : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellr_ab : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellr_eta : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellr_omega : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellr_roots : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val elladd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellan : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellanq_zv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellanal_globalred :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ellap :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellap_cm_fast :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> Signed.long

val ellbasechar : ('kind, 'structure) t -> ('kind, 'structure) t
val ellbsd : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellcard :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellchangecurve :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellchangeinvert : ('kind, 'structure) t -> ('kind, 'structure) t

val ellchangepoint :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellchangepointinv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val elldivpol :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val elleisnum :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val elleta : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val elleulerf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellff_get_card : ('kind, 'structure) t -> ('kind, 'structure) t
val ellff_get_gens : ('kind, 'structure) t -> ('kind, 'structure) t
val ellff_get_group : ('kind, 'structure) t -> ('kind, 'structure) t
val ellff_get_o : ('kind, 'structure) t -> ('kind, 'structure) t
val ellff_get_p : ('kind, 'structure) t -> ('kind, 'structure) t
val ellff_get_m : ('kind, 'structure) t -> ('kind, 'structure) t
val ellff_get_d : ('kind, 'structure) t -> ('kind, 'structure) t
val ellfromj : ('kind, 'structure) t -> ('kind, 'structure) t
val ellgenerators : ('kind, 'structure) t -> ('kind, 'structure) t
val ellglobalred : ('kind, 'structure) t -> ('kind, 'structure) t

val ellgroup :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellgroup0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellheight0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellheight :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellheightmatrix :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellheightoo :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellinit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellintegralmodel :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ellintegralmodel_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val elliscm : ('kind, 'structure) t -> Signed.long

val ellisoncurve :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellisotree : ('kind, 'structure) t -> ('kind, 'structure) t
val ellissupersingular : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val elljissupersingular : ('kind, 'structure) t -> int

val elllseries :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val elllocalred :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val elllog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellminimaldisc : ('kind, 'structure) t -> ('kind, 'structure) t

val ellminimalmodel :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ellminimaltwist : ('kind, 'structure) t -> ('kind, 'structure) t

val ellminimaltwist0 :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellminimaltwistcond : ('kind, 'structure) t -> ('kind, 'structure) t

val ellmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellnf_vecarea :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellnf_veceta : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellnf_vecomega :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellneg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellorder :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellorder_q : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val ellordinate :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellpadicheight0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellpadicheightmatrix :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellperiods :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellrandom : ('kind, 'structure) t -> ('kind, 'structure) t
val ellrootno : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long
val ellrootno_global : ('kind, 'structure) t -> Signed.long

val ellsaturation :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val ellsea : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellsigma :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val ellsub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellsupersingularj : ('kind, 'structure) t -> ('kind, 'structure) t
val elltamagawa : ('kind, 'structure) t -> ('kind, 'structure) t
val elltaniyama : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val elltatepairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val elltors : ('kind, 'structure) t -> ('kind, 'structure) t
val elltors0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val elltors_psylow :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val elltrace :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val elltwist :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellweilpairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellwp :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellwp0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val ellwpseries :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellxn :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellzeta :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val oncurve : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val orderell :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val pointell :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val point_to_a4a6 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val point_to_a4a6_fl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  ('kind, 'structure) t

val zell :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qp_agm2_sequence :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qp_ascending_landen :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val qp_descending_landen :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val ellformaldifferential :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellformalexp :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellformallog :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellformalpoint :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellformalw :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellnonsingularmultiple :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellpadicl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellpadicbsd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellpadicfrobenius :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val ellpadicheight :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellpadiclog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellpadicregulator :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellpadics2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ell2cover : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ellrank :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ellrankinit : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val hyperell_locally_soluble :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val nf_hyperell_locally_soluble :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val nfhilbert :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val nfhilbert0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val ellisdivisible :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val ellisogenyapply :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellisogeny :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val ellisomat :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ellweilcurve :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val flxq_elldivpolmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val fp_ellcard_sea :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fq_ellcard_sea :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fq_elldivpolmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ellmodulareqn :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val externstr : string -> ('kind, 'structure) t
val gp_filter : string -> string
val gpextern : string -> ('kind, 'structure) t
val gpsystem : string -> Signed.long
val readstr : string -> ('kind, 'structure) t
val gentogenstr_nospace : ('kind, 'structure) t -> ('kind, 'structure) t
val gentogenstr : ('kind, 'structure) t -> ('kind, 'structure) t
val gentotexstr : ('kind, 'structure) t -> string
val gentostr : ('kind, 'structure) t -> string
val gentostr_raw : ('kind, 'structure) t -> string
val gentostr_unquoted : ('kind, 'structure) t -> string
val str : ('kind, 'structure) t -> ('kind, 'structure) t
val strexpand : ('kind, 'structure) t -> ('kind, 'structure) t
val strtex : ('kind, 'structure) t -> ('kind, 'structure) t
val brute : ('kind, 'structure) t -> char -> Signed.long -> unit
val dbggen : ('kind, 'structure) t -> Signed.long -> unit
val error0 : ('kind, 'structure) t -> unit
val dbg_pari_heap : unit -> unit
val err_flush : unit -> unit
val err_printf : string -> unit
val gp_getenv : string -> ('kind, 'structure) t
val gp_fileclose : Signed.long -> unit
val gp_fileextern : string -> Signed.long
val gp_fileflush : Signed.long -> unit
val gp_fileflush0 : ('kind, 'structure) t -> unit
val gp_fileopen : string -> string -> Signed.long
val gp_fileread : Signed.long -> ('kind, 'structure) t
val gp_filereadstr : Signed.long -> ('kind, 'structure) t
val gp_filewrite : Signed.long -> string -> unit
val gp_filewrite1 : Signed.long -> string -> unit
val gp_read_file : string -> ('kind, 'structure) t
val gp_read_str_multiline : string -> string -> ('kind, 'structure) t
val gp_readvec_file : string -> ('kind, 'structure) t
val gpinstall : string -> string -> string -> string -> unit
val gsprintf : string -> ('kind, 'structure) t
val itostr : ('kind, 'structure) t -> string
val matbrute : ('kind, 'structure) t -> char -> Signed.long -> unit
val os_getenv : string -> string
val uordinal : pari_ulong -> string
val outmat : ('kind, 'structure) t -> unit
val output : ('kind, 'structure) t -> unit
val rgv_to_str : ('kind, 'structure) t -> Signed.long -> string
val pari_add_hist : ('kind, 'structure) t -> Signed.long -> Signed.long -> unit
val pari_ask_confirm : string -> unit
val pari_flush : unit -> unit
val pari_get_hist : Signed.long -> ('kind, 'structure) t
val pari_get_histrtime : Signed.long -> Signed.long
val pari_get_histtime : Signed.long -> Signed.long
val pari_get_homedir : string -> string
val pari_histtime : Signed.long -> ('kind, 'structure) t
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
val pari_sprint0 : string -> ('kind, 'structure) t -> Signed.long -> string
val print : ('kind, 'structure) t -> unit
val printp : ('kind, 'structure) t -> unit
val print1 : ('kind, 'structure) t -> unit
val printf0 : string -> ('kind, 'structure) t -> unit
val printsep : string -> ('kind, 'structure) t -> unit
val printsep1 : string -> ('kind, 'structure) t -> unit
val printtex : ('kind, 'structure) t -> unit
val stack_sprintf : string -> string
val str_init : pari_str Ctypes.structure Ctypes_static.ptr -> int -> unit
val str_printf : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val str_putc : pari_str Ctypes.structure Ctypes_static.ptr -> char -> unit
val str_puts : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val strftime_expand : string -> string -> Signed.long -> unit
val strprintf : string -> ('kind, 'structure) t -> ('kind, 'structure) t
val term_color : Signed.long -> unit
val term_get_color : string -> Signed.long -> string
val texe : ('kind, 'structure) t -> char -> Signed.long -> unit
val warning0 : ('kind, 'structure) t -> unit
val write0 : string -> ('kind, 'structure) t -> unit
val write1 : string -> ('kind, 'structure) t -> unit
val writebin : string -> ('kind, 'structure) t -> unit
val writetex : string -> ('kind, 'structure) t -> unit
val bincopy_relink : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val bitprecision0 :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val bitprecision00 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val break0 : Signed.long -> ('kind, 'structure) t

val call0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val closure_callgen0prec :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val closure_callgen1 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val closure_callgen1prec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val closure_callgen2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val closure_callgenall :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val closure_callgenvec :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val closure_callgenvecdef :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val closure_callgenvecdefprec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val closure_callgenvecprec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val closure_callvoid1 : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val closure_context : Signed.long -> Signed.long -> Signed.long
val closure_disassemble : ('kind, 'structure) t -> unit
val closure_err : Signed.long -> unit

val closure_evalbrk :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val closure_evalgen : ('kind, 'structure) t -> ('kind, 'structure) t
val closure_evalnobrk : ('kind, 'structure) t -> ('kind, 'structure) t
val closure_evalres : ('kind, 'structure) t -> ('kind, 'structure) t
val closure_evalvoid : ('kind, 'structure) t -> unit
val closure_func_err : unit -> string

val closure_trapgen :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val copybin_unlink : ('kind, 'structure) t -> ('kind, 'structure) t
val getlocalprec : Signed.long -> Signed.long
val getlocalbitprec : Signed.long -> Signed.long
val get_lex : Signed.long -> ('kind, 'structure) t
val get_localprec : unit -> Signed.long
val get_localbitprec : unit -> Signed.long

val gp_call :
  unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t

val gp_callprec :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gp_call2 :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gp_callbool : unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long
val gp_callvoid : unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long

val gp_eval :
  unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t

val gp_evalbool : unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long

val gp_evalprec :
  unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gp_evalupto :
  unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t

val gp_evalvoid : unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long
val localprec : ('kind, 'structure) t -> unit
val localbitprec : ('kind, 'structure) t -> unit
val loop_break : unit -> Signed.long
val next0 : Signed.long -> ('kind, 'structure) t
val pareval : ('kind, 'structure) t -> ('kind, 'structure) t
val pari_self : unit -> ('kind, 'structure) t

val parsum :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val parvector : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val pop_lex : Signed.long -> unit
val pop_localprec : unit -> unit
val precision0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val precision00 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val push_lex : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val push_localbitprec : Signed.long -> unit
val push_localprec : Signed.long -> unit
val return0 : ('kind, 'structure) t -> ('kind, 'structure) t
val set_lex : Signed.long -> ('kind, 'structure) t -> unit

val forcomposite_init :
  forcomposite_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int

val forcomposite_next :
  forcomposite_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val forprime_next :
  forprime_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val forprime_init :
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  int

val forprimestep_init :
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
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
val prodprimes : unit -> ('kind, 'structure) t

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

val ff_1 : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_frobenius : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ff_z_z_muldiv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ff_q_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_z_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_z_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_add :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_charpoly : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_conjvec : ('kind, 'structure) t -> ('kind, 'structure) t

val ff_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_ellcard : ('kind, 'structure) t -> ('kind, 'structure) t

val ff_ellcard_sea :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ff_ellgens : ('kind, 'structure) t -> ('kind, 'structure) t

val ff_ellgroup :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ff_elllog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ff_ellmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ff_ellorder :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ff_elltwist : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_ellrandom : ('kind, 'structure) t -> ('kind, 'structure) t

val ff_elltatepairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ff_ellweilpairing :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ff_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val ff_equal0 : ('kind, 'structure) t -> int
val ff_equal1 : ('kind, 'structure) t -> int
val ff_equalm1 : ('kind, 'structure) t -> int
val ff_f : ('kind, 'structure) t -> Signed.long
val ff_gen : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_inv : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_issquare : ('kind, 'structure) t -> Signed.long

val ff_issquareall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val ff_ispower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val ff_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ff_map :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_minpoly : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_mod : ('kind, 'structure) t -> ('kind, 'structure) t

val ff_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_mul2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ff_neg : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_neg_i : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_norm : ('kind, 'structure) t -> ('kind, 'structure) t

val ff_order :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_p : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_p_i : ('kind, 'structure) t -> ('kind, 'structure) t

val ff_pow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_primroot :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ff_q : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_samefield : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val ff_sqr : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_sqrt : ('kind, 'structure) t -> ('kind, 'structure) t

val ff_sqrtn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ff_sub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ff_to_f2xq : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_to_f2xq_i : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_to_flxq : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_to_flxq_i : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_to_fpxq : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_to_fpxq_i : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_trace : ('kind, 'structure) t -> ('kind, 'structure) t
val ff_var : ('kind, 'structure) t -> Signed.long
val ff_zero : ('kind, 'structure) t -> ('kind, 'structure) t

val ffm_ffc_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffm_ffc_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffm_ffc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffm_deplin :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffm_det :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffm_gauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffm_image :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffm_indexrank :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffm_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffm_invimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffm_ker :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffm_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffm_rank : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val ffm_suppl :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffx_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffx_ddf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffx_degfact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffx_disc :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffx_extgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ffx_factor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffx_factor_squarefree :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffx_halfgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffx_halfgcd_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ffx_ispower :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val ffx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffx_preimage :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffx_preimagerel :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffx_rem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffx_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffx_roots :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffx_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffxq_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffxq_minpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffxq_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_to_ffx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fq_to_ff :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_ff_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffembed :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffextend :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fffrobenius : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ffgen : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ffinvmap : ('kind, 'structure) t -> ('kind, 'structure) t

val fflog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val ffmap :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffmaprel :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffcompomap :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fforder :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ffprimroot :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ffrandom : ('kind, 'structure) t -> ('kind, 'structure) t

val rg_is_ff :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> int

val rgc_is_ffc :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> int

val rgm_is_ffm :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> int

val p_to_ff : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val tp_to_ff :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val flx_factcyclo :
  pari_ulong -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val fpx_factcyclo :
  pari_ulong -> ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val factormodcyclo :
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val checkgal : ('kind, 'structure) t -> ('kind, 'structure) t

val checkgroup :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val checkgroupelts : ('kind, 'structure) t -> ('kind, 'structure) t

val embed_disc :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val embed_roots : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val galois_group : ('kind, 'structure) t -> ('kind, 'structure) t

val galoisconj :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val galoisconj0 :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val galoisconjclasses : ('kind, 'structure) t -> ('kind, 'structure) t
val galoisexport : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val galoisfixedfield :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val galoisidentify : ('kind, 'structure) t -> ('kind, 'structure) t

val galoisinit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val galoisisabelian :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val galoisisnormal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val galoispermtopol :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val galoissplittinginit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val galoissubgroups : ('kind, 'structure) t -> ('kind, 'structure) t

val galoissubfields :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val numberofconjugates : ('kind, 'structure) t -> Signed.long -> Signed.long
val polgalois : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val galoisnbpol : Signed.long -> ('kind, 'structure) t
val galoisgetgroup : Signed.long -> Signed.long -> ('kind, 'structure) t
val galoisgetname : Signed.long -> Signed.long -> ('kind, 'structure) t

val galoisgetpol :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val conj_i : ('kind, 'structure) t -> ('kind, 'structure) t
val conjvec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val divrunextu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val gadd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gaddsg : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val gconj : ('kind, 'structure) t -> ('kind, 'structure) t

val gdiv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gdivgs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gdivgu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val gdivgunextu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val ginv : ('kind, 'structure) t -> ('kind, 'structure) t

val gmul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gmul2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gmulsg : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val gmulug : pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t
val gsqr : ('kind, 'structure) t -> ('kind, 'structure) t

val gsub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gsubsg : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val mulcxi : ('kind, 'structure) t -> ('kind, 'structure) t
val mulcxmi : ('kind, 'structure) t -> ('kind, 'structure) t
val mulcxpowis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val qdivii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qdiviu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val qdivis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ser_normalize : ('kind, 'structure) t -> ('kind, 'structure) t

val gassoc_proto :
  (('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val map_proto_g :
  (('kind, 'structure) t -> ('kind, 'structure) t) Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val map_proto_lg :
  (('kind, 'structure) t -> Signed.long) Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val map_proto_lgl :
  (('kind, 'structure) t -> Signed.long -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val q_lval : ('kind, 'structure) t -> pari_ulong -> Signed.long

val q_lvalrem :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val q_pval : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val q_pvalrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val rgx_val : ('kind, 'structure) t -> Signed.long

val rgx_valrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val rgx_valrem_inexact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val rgxv_maxdegree : ('kind, 'structure) t -> Signed.long
val zv_z_dvd : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val zv_pval : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val zv_pvalrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val zv_lval : ('kind, 'structure) t -> pari_ulong -> Signed.long

val zv_lvalrem :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val zx_lvalrem :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val zx_pval : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val zx_pvalrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val z_lvalrem_stop :
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong ->
  int Ctypes_static.ptr ->
  Signed.long

val cgetp : ('kind, 'structure) t -> ('kind, 'structure) t
val cvstop2 : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val cvtop :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val cvtop2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val cx_approx_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val cx_approx0 : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val gabs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gaffect : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val gaffsg : Signed.long -> ('kind, 'structure) t -> unit
val gcmp : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val gequal0 : ('kind, 'structure) t -> int
val gequal1 : ('kind, 'structure) t -> int
val gequalx : ('kind, 'structure) t -> int
val gequalm1 : ('kind, 'structure) t -> int
val gcmpsg : Signed.long -> ('kind, 'structure) t -> int

val gcvtop :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gequal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val gequalsg : Signed.long -> ('kind, 'structure) t -> int
val gexpo : ('kind, 'structure) t -> Signed.long
val gexpo_safe : ('kind, 'structure) t -> Signed.long
val gpexponent : ('kind, 'structure) t -> ('kind, 'structure) t

val gpvaluation :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gvaluation : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long
val gidentical : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val glength : ('kind, 'structure) t -> Signed.long

val gmax :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gmaxgs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gmin :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gmings : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gneg : ('kind, 'structure) t -> ('kind, 'structure) t
val gneg_i : ('kind, 'structure) t -> ('kind, 'structure) t
val gsigne : ('kind, 'structure) t -> int
val gtolist : ('kind, 'structure) t -> ('kind, 'structure) t
val gtolong : ('kind, 'structure) t -> Signed.long
val lexcmp : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val listinsert :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val listpop : ('kind, 'structure) t -> Signed.long -> unit
val listpop0 : ('kind, 'structure) t -> Signed.long -> unit

val listput :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val listput0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long -> unit

val listsort : ('kind, 'structure) t -> Signed.long -> unit
val matsize : ('kind, 'structure) t -> ('kind, 'structure) t
val mklist : unit -> ('kind, 'structure) t
val mklist_typ : Signed.long -> ('kind, 'structure) t
val mklistcopy : ('kind, 'structure) t -> ('kind, 'structure) t
val mkmap : unit -> ('kind, 'structure) t
val normalizeser : ('kind, 'structure) t -> ('kind, 'structure) t
val normalizepol : ('kind, 'structure) t -> ('kind, 'structure) t

val normalizepol_approx :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val normalizepol_lg :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val padic_to_fl : ('kind, 'structure) t -> pari_ulong -> pari_ulong

val padic_to_fp :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val quadtofp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val sizedigit : ('kind, 'structure) t -> Signed.long
val u_lval : pari_ulong -> pari_ulong -> Signed.long

val u_lvalrem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val u_lvalrem_stop :
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  int Ctypes_static.ptr ->
  Signed.long

val u_pval : pari_ulong -> ('kind, 'structure) t -> Signed.long

val u_pvalrem :
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong Ctypes_static.ptr ->
  Signed.long

val vecindexmax : ('kind, 'structure) t -> Signed.long
val vecindexmin : ('kind, 'structure) t -> Signed.long

val vecmax0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val vecmax : ('kind, 'structure) t -> ('kind, 'structure) t

val vecmin0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val vecmin : ('kind, 'structure) t -> ('kind, 'structure) t
val z_lval : Signed.long -> pari_ulong -> Signed.long

val z_lvalrem :
  Signed.long -> pari_ulong -> Signed.long Ctypes_static.ptr -> Signed.long

val z_pval : Signed.long -> ('kind, 'structure) t -> Signed.long

val z_pvalrem :
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val zx_lval : ('kind, 'structure) t -> Signed.long -> Signed.long
val hgmcyclo : ('kind, 'structure) t -> ('kind, 'structure) t
val hgmalpha : ('kind, 'structure) t -> ('kind, 'structure) t
val hgmgamma : ('kind, 'structure) t -> ('kind, 'structure) t

val hgminit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hgmparams : ('kind, 'structure) t -> ('kind, 'structure) t

val hgmeulerfactor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val hgmcoef :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val hgmcoefs :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val hgmtwist : ('kind, 'structure) t -> ('kind, 'structure) t
val hgmissymmetrical : ('kind, 'structure) t -> Signed.long
val hgmbydegree : Signed.long -> ('kind, 'structure) t

val lfunhgm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qp_zeta : ('kind, 'structure) t -> ('kind, 'structure) t

val lerchphi :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val lerchzeta :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zetahurwitz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rgx_to_ser : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_to_ser_inexact :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gtoser :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val gtoser_prec :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rfrac_to_ser : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rfrac_to_ser_i :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rfracrecip_to_ser_absolute :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rfracrecip :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val scalarser :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val sertoser : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val toser_i : ('kind, 'structure) t -> ('kind, 'structure) t

val rgv_to_ser :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ser0 :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val padic_to_q : ('kind, 'structure) t -> ('kind, 'structure) t
val padic_to_q_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val qpv_to_qv : ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_rgv_mulrealsym :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_mulreal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_cxeval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgx_deflate_max :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgx_deflate_order : ('kind, 'structure) t -> Signed.long
val rgx_degree : ('kind, 'structure) t -> Signed.long -> Signed.long
val rgx_integ : ('kind, 'structure) t -> ('kind, 'structure) t

val rgxy_cxevalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zx_deflate_order : ('kind, 'structure) t -> Signed.long

val zx_deflate_max :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val ceil_safe : ('kind, 'structure) t -> ('kind, 'structure) t
val ceilr : ('kind, 'structure) t -> ('kind, 'structure) t
val centerlift : ('kind, 'structure) t -> ('kind, 'structure) t
val centerlift0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val compo : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val deg1pol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val deg1pol_shallow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val deg2pol_shallow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val degree : ('kind, 'structure) t -> Signed.long
val denom : ('kind, 'structure) t -> ('kind, 'structure) t
val denom_i : ('kind, 'structure) t -> ('kind, 'structure) t

val denominator :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val deriv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
(** [deriv] *)

val derivn :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val derivser : ('kind, 'structure) t -> ('kind, 'structure) t

val diffop :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val diffop0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val diviiround :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val floor_safe : ('kind, 'structure) t -> ('kind, 'structure) t
val gceil : ('kind, 'structure) t -> ('kind, 'structure) t

val gcvtoi :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val gdeflate :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val gdivent :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gdiventgs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gdiventsg : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val gdiventres :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gdivmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val gdivround :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gdvd : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val geq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val geval : ('kind, 'structure) t -> ('kind, 'structure) t
val gfloor : ('kind, 'structure) t -> ('kind, 'structure) t
val gtrunc2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gfrac : ('kind, 'structure) t -> ('kind, 'structure) t

val gge :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ggrando : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ggt :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gimag : ('kind, 'structure) t -> ('kind, 'structure) t
val gisexactzero : ('kind, 'structure) t -> ('kind, 'structure) t

val gle :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val glt :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gmodgs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gmodsg : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val gmodulo :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gmodulsg : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val gmodulss : Signed.long -> Signed.long -> ('kind, 'structure) t

val gne :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gnot : ('kind, 'structure) t -> ('kind, 'structure) t
val gpolvar : ('kind, 'structure) t -> ('kind, 'structure) t

val gppadicprec :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gppoldegree : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gprecision : ('kind, 'structure) t -> Signed.long
val gpserprec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val greal : ('kind, 'structure) t -> ('kind, 'structure) t

val grndtoi :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val ground : ('kind, 'structure) t -> ('kind, 'structure) t
val gshift : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gsubst :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gsubstpol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gsubstvec :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gtocol : ('kind, 'structure) t -> ('kind, 'structure) t
val gtocol0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtocolrev : ('kind, 'structure) t -> ('kind, 'structure) t
val gtocolrev0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtopoly : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtopolyrev : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtovec : ('kind, 'structure) t -> ('kind, 'structure) t
val gtovec0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtovecrev : ('kind, 'structure) t -> ('kind, 'structure) t
val gtovecrev0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtovecsmall : ('kind, 'structure) t -> ('kind, 'structure) t
val gtovecsmall0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtrunc : ('kind, 'structure) t -> ('kind, 'structure) t
val gvar : ('kind, 'structure) t -> Signed.long
val gvar2 : ('kind, 'structure) t -> Signed.long

val hqfeval :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val imag_i : ('kind, 'structure) t -> ('kind, 'structure) t
val integ : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val integser : ('kind, 'structure) t -> ('kind, 'structure) t
val ser_inv : ('kind, 'structure) t -> ('kind, 'structure) t
val iscomplex : ('kind, 'structure) t -> int
val isexactzero : ('kind, 'structure) t -> int
val isrationalzeroscalar : ('kind, 'structure) t -> int
val isinexact : ('kind, 'structure) t -> int
val isinexactreal : ('kind, 'structure) t -> int

val isint :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> int

val isrationalzero : ('kind, 'structure) t -> int
val issmall : ('kind, 'structure) t -> Signed.long Ctypes_static.ptr -> int
val lift : ('kind, 'structure) t -> ('kind, 'structure) t
val lift_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val lift0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val liftall : ('kind, 'structure) t -> ('kind, 'structure) t
val liftall_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val liftint : ('kind, 'structure) t -> ('kind, 'structure) t
val liftint_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val liftpol : ('kind, 'structure) t -> ('kind, 'structure) t
val liftpol_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val mkcoln : Signed.long -> ('kind, 'structure) t
val mkintn : Signed.long -> ('kind, 'structure) t
val mkpoln : Signed.long -> ('kind, 'structure) t
val mkvecn : Signed.long -> ('kind, 'structure) t
val mkvecsmalln : Signed.long -> ('kind, 'structure) t

val modrr_safe :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val modrr_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mulreal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val numer : ('kind, 'structure) t -> ('kind, 'structure) t
val numer_i : ('kind, 'structure) t -> ('kind, 'structure) t

val numerator :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val padicprec : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long
val padicprec_relative : ('kind, 'structure) t -> Signed.long

val polcoef :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val polcoef_i :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val poldegree : ('kind, 'structure) t -> Signed.long -> Signed.long

val poleval :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val pollead : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val precision : ('kind, 'structure) t -> Signed.long

val qf_apply_rgm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qf_apply_zm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfb_apply_zm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfbil :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfeval :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfeval0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfevalb :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfnorm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val real_i : ('kind, 'structure) t -> ('kind, 'structure) t

val round0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val roundr : ('kind, 'structure) t -> ('kind, 'structure) t
val roundr_safe : ('kind, 'structure) t -> ('kind, 'structure) t
val scalarpol : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val scalarpol_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ser_unscale :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val serprec : ('kind, 'structure) t -> Signed.long -> Signed.long
val serreverse : ('kind, 'structure) t -> ('kind, 'structure) t
val simplify : ('kind, 'structure) t -> ('kind, 'structure) t
val simplify_shallow : ('kind, 'structure) t -> ('kind, 'structure) t

val tayl :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val trunc0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val uu32toi : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val uu32toineg : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val vars_sort_inplace : ('kind, 'structure) t -> ('kind, 'structure) t
val vars_to_rgxv : ('kind, 'structure) t -> ('kind, 'structure) t
val variables_vecsmall : ('kind, 'structure) t -> ('kind, 'structure) t
val variables_vec : ('kind, 'structure) t -> ('kind, 'structure) t

val genus2red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val genus2igusa : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gchar_conductor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gchar_identify :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gcharalgebraic :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gcharduallog :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gchareval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gchari_lfun :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gcharinit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gcharisalgebraic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val gcharlocal :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val gcharlog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gcharnewprec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val is_gchar_group : ('kind, 'structure) t -> int

val lfungchar :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecan_gchar :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val eulerf_gchar :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val group_ident : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val group_ident_trans :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

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
  ('kind, 'structure) t

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
  (('kind, 'structure) t -> ('kind, 'structure) t -> int)
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
  hashtable Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val hash_values :
  hashtable Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

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
val hash_gen : ('kind, 'structure) t -> pari_ulong
val hash_zv : ('kind, 'structure) t -> pari_ulong

val zx_hyperellred :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val hyperellcharpoly : ('kind, 'structure) t -> ('kind, 'structure) t

val hyperellchangecurve :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hyperelldisc : ('kind, 'structure) t -> ('kind, 'structure) t
val hyperellisoncurve : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val hyperellminimaldisc :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hyperellminimalmodel :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val hyperellpadicfrobenius0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val hyperellpadicfrobenius :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val hyperellred :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val nfhyperellpadicfrobenius :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) t

val hypergeom :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val airy : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgm_hnfall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t

val zm_hnf : ('kind, 'structure) t -> ('kind, 'structure) t
val zm_hnf_knapsack : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_hnfall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t

val zm_hnfall_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t

val zm_hnfcenter : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_hnflll :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int ->
  ('kind, 'structure) t

val zv_extgcd : ('kind, 'structure) t -> ('kind, 'structure) t

val zv_snfall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zv_snf_group :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zv_snf_rank_u : ('kind, 'structure) t -> pari_ulong -> Signed.long
val zv_snf_trunc : ('kind, 'structure) t -> unit

val zm_hnfmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_hnfmodall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zm_hnfmodall_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zm_hnfmodid :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_hnfmodprime :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_hnfperm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zm_snfclean :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val zm_snf : ('kind, 'structure) t -> ('kind, 'structure) t

val zm_snf_group :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zm_snfall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zm_snfall_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t

val zv_snfclean : ('kind, 'structure) t -> ('kind, 'structure) t

val zpm_echelon :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gsmith : ('kind, 'structure) t -> ('kind, 'structure) t
val gsmithall : ('kind, 'structure) t -> ('kind, 'structure) t
val hnf : ('kind, 'structure) t -> ('kind, 'structure) t

val hnf_divscale :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val hnf_invscale :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hnf_solve :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hnf_invimage :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hnfall : ('kind, 'structure) t -> ('kind, 'structure) t
val hnfdivide : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val hnflll : ('kind, 'structure) t -> ('kind, 'structure) t

val hnfmerge_get_1 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hnfmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hnfmodid :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val hnfperm : ('kind, 'structure) t -> ('kind, 'structure) t

val matfrobenius :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mathnf0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val matsnf0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val smith : ('kind, 'structure) t -> ('kind, 'structure) t
val smithall : ('kind, 'structure) t -> ('kind, 'structure) t
val smithclean : ('kind, 'structure) t -> ('kind, 'structure) t
val snfrank : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val zlm_echelon :
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val zv_snf_rank : ('kind, 'structure) t -> pari_ulong -> Signed.long

val z_ecm :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) t

val z_factor : ('kind, 'structure) t -> ('kind, 'structure) t

val z_factor_limit :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val z_factor_until :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_issmooth : ('kind, 'structure) t -> pari_ulong -> Signed.long

val z_issmooth_fact :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val z_issquarefree : ('kind, 'structure) t -> Signed.long

val z_pollardbrent :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val absz_factor : ('kind, 'structure) t -> ('kind, 'structure) t

val absz_factor_limit :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val absz_factor_limit_strict :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val coreu : pari_ulong -> pari_ulong
val coreu_fact : ('kind, 'structure) t -> pari_ulong
val factorint : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val factoru : pari_ulong -> ('kind, 'structure) t
val tridiv_boundu : pari_ulong -> pari_ulong
val ifac_isprime : ('kind, 'structure) t -> int

val ifac_next :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  int

val ifac_read :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  int

val ifac_skip : ('kind, 'structure) t -> unit
val ifac_start : ('kind, 'structure) t -> int -> ('kind, 'structure) t

val is_357_power :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  int

val is_pth_power :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  pari_ulong ->
  int

val ispowerful : ('kind, 'structure) t -> Signed.long
val maxomegau : pari_ulong -> Signed.long
val maxomegaoddu : pari_ulong -> Signed.long
val moebius : ('kind, 'structure) t -> Signed.long
val moebiusu : pari_ulong -> Signed.long
val moebiusu_fact : ('kind, 'structure) t -> Signed.long
val nextprime : ('kind, 'structure) t -> ('kind, 'structure) t
val precprime : ('kind, 'structure) t -> ('kind, 'structure) t
val radicalu : pari_ulong -> pari_ulong
val tridiv_bound : ('kind, 'structure) t -> pari_ulong

val uis_357_power :
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  int

val uis_357_powermod : pari_ulong -> pari_ulong Ctypes_static.ptr -> int
val unextprime : pari_ulong -> pari_ulong
val uprecprime : pari_ulong -> pari_ulong
val vecfactorsquarefreeu : pari_ulong -> pari_ulong -> ('kind, 'structure) t

val vecfactorsquarefreeu_coprime :
  pari_ulong -> pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecfactoru_i : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val vecfactoru : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val vecfactoroddu_i : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val vecfactoroddu : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val vecsquarefreeu : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val chk_gerepileupto : ('kind, 'structure) t -> int

val copy_bin :
  ('kind, 'structure) t -> genbin Ctypes.structure Ctypes_static.ptr

val copy_bin_canon :
  ('kind, 'structure) t -> genbin Ctypes.structure Ctypes_static.ptr

val dbg_gerepile : pari_ulong -> unit
val dbg_gerepileupto : ('kind, 'structure) t -> unit
val errname : ('kind, 'structure) t -> ('kind, 'structure) t
val gclone : ('kind, 'structure) t -> ('kind, 'structure) t
val gcloneref : ('kind, 'structure) t -> ('kind, 'structure) t
val gclone_refc : ('kind, 'structure) t -> unit
val gcopy : ('kind, 'structure) t -> ('kind, 'structure) t

val gcopy_avma :
  ('kind, 'structure) t -> pari_ulong Ctypes_static.ptr -> ('kind, 'structure) t

val gcopy_lg : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gerepile :
  pari_ulong -> pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t

val gerepileallsp : pari_ulong -> pari_ulong -> int -> unit

val gerepilecoeffssp :
  pari_ulong -> pari_ulong -> Signed.long Ctypes_static.ptr -> int -> unit

val gerepilemanysp :
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr Ctypes_static.ptr ->
  int ->
  unit

val getheap : unit -> ('kind, 'structure) t
val gsizeword : ('kind, 'structure) t -> Signed.long
val gsizebyte : ('kind, 'structure) t -> Signed.long
val gunclone : ('kind, 'structure) t -> unit
val gunclone_deep : ('kind, 'structure) t -> unit
val listcopy : ('kind, 'structure) t -> ('kind, 'structure) t
val listinit : ('kind, 'structure) t -> ('kind, 'structure) t
val msgtimer : string -> unit
val name_numerr : string -> Signed.long
val new_chunk_resize : int -> unit
val newblock : int -> ('kind, 'structure) t
val numerr_name : Signed.long -> string
val obj_check : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val obj_checkbuild :
  ('kind, 'structure) t ->
  Signed.long ->
  (('kind, 'structure) t -> ('kind, 'structure) t) Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val obj_checkbuild_padicprec :
  ('kind, 'structure) t ->
  Signed.long ->
  (('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) t

val obj_checkbuild_realprec :
  ('kind, 'structure) t ->
  Signed.long ->
  (('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) t

val obj_checkbuild_prec :
  ('kind, 'structure) t ->
  Signed.long ->
  (('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (('kind, 'structure) t -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) t

val obj_free : ('kind, 'structure) t -> unit
val obj_init : Signed.long -> Signed.long -> ('kind, 'structure) t

val obj_insert :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val obj_insert_shallow :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val obj_reinit : ('kind, 'structure) t -> ('kind, 'structure) t
val pari_add_function : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_defaults_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_close : unit -> unit
val pari_close_opts : pari_ulong -> unit
val pari_compile_str : string -> ('kind, 'structure) t
val pari_daemon : unit -> int
val pari_err : int -> unit
val pari_err_last : unit -> ('kind, 'structure) t
val pari_err2str : ('kind, 'structure) t -> string
val pari_init_opts : int -> pari_ulong -> pari_ulong -> unit
val pari_init : int -> pari_ulong -> unit
val pari_stackcheck_init : unit Ctypes_static.ptr -> unit
val pari_sighandler : int -> unit
val pari_sig_init : (int -> unit) Ctypes_static.static_funptr -> unit

val pari_thread_alloc :
  pari_thread Ctypes.structure Ctypes_static.ptr ->
  int ->
  ('kind, 'structure) t ->
  unit

val pari_thread_close : unit -> unit
val pari_thread_free : pari_thread Ctypes.structure Ctypes_static.ptr -> unit
val pari_thread_init : unit -> unit

val pari_thread_start :
  pari_thread Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val pari_thread_valloc :
  pari_thread Ctypes.structure Ctypes_static.ptr ->
  int ->
  int ->
  ('kind, 'structure) t ->
  unit

val pari_version : unit -> ('kind, 'structure) t
val pari_warn : int -> unit
val paristack_newrsize : pari_ulong -> unit
val paristack_resize : pari_ulong -> unit
val paristack_setsize : int -> int -> unit
val parivstack_resize : pari_ulong -> unit
val parivstack_reset : unit -> unit
val setalldebug : Signed.long -> unit
val setdebug : string -> Signed.long -> ('kind, 'structure) t
val shiftaddress : ('kind, 'structure) t -> Signed.long -> unit
val shiftaddress_canon : ('kind, 'structure) t -> Signed.long -> unit
val timer : unit -> Signed.long
val timer_delay : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long
val timer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val timer_printf :
  pari_timer Ctypes.structure Ctypes_static.ptr -> string -> unit

val timer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit
val timer2 : unit -> Signed.long

val trap0 :
  string ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val traverseheap :
  (('kind, 'structure) t -> unit Ctypes_static.ptr -> unit)
  Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  unit

val walltimer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit

val walltimer_delay :
  pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val walltimer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val contfraceval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val contfracinit : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val intcirc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val intfuncinit :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val intnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val intnumgauss :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val intnumgaussinit : Signed.long -> Signed.long -> ('kind, 'structure) t

val intnuminit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val intnumosc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val intnumromb :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val intnumromb_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val prodeulerrat :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val prodnumrat :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val quodif : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val sumeulerrat :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val sumnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumnumap :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumnumapinit : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val sumnuminit : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val sumnumlagrangeinit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumnumlagrange :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumnummonien :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumnummonieninit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumnumrat :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumnumsidi :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  float ->
  Signed.long ->
  ('kind, 'structure) t

val z_isanypower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val z_ispow2 : ('kind, 'structure) t -> Signed.long

val z_ispowerall :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val z_issquareall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val zn_ispower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val zn_issquare : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long
val zp_issquare : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val gisanypower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val gissquare : ('kind, 'structure) t -> ('kind, 'structure) t

val gissquareall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val ispolygonal :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val ispower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val isprimepower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val ispseudoprimepower :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val issquare : ('kind, 'structure) t -> Signed.long

val issquareall :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long

val sqrtint : ('kind, 'structure) t -> ('kind, 'structure) t

val sqrtint0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val uisprimepower : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long
val uissquare : pari_ulong -> Signed.long
val uissquareall : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val ulogintall :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val padicfields0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val padicfields :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val bnrclassfield :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rnfkummer :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val is_linit : ('kind, 'structure) t -> Signed.long
val ldata_get_an : ('kind, 'structure) t -> ('kind, 'structure) t
val ldata_get_dual : ('kind, 'structure) t -> ('kind, 'structure) t
val ldata_get_gammavec : ('kind, 'structure) t -> ('kind, 'structure) t
val ldata_get_degree : ('kind, 'structure) t -> Signed.long
val ldata_get_k : ('kind, 'structure) t -> ('kind, 'structure) t
val ldata_get_k1 : ('kind, 'structure) t -> ('kind, 'structure) t
val ldata_get_conductor : ('kind, 'structure) t -> ('kind, 'structure) t
val ldata_get_rootno : ('kind, 'structure) t -> ('kind, 'structure) t
val ldata_get_residue : ('kind, 'structure) t -> ('kind, 'structure) t
val ldata_get_type : ('kind, 'structure) t -> Signed.long
val ldata_isreal : ('kind, 'structure) t -> Signed.long
val linit_get_type : ('kind, 'structure) t -> Signed.long
val linit_get_ldata : ('kind, 'structure) t -> ('kind, 'structure) t
val linit_get_tech : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_domain : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_dom : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_factgammavec : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_step : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_pol : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_residue : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_k2 : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_w2 : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_expot : ('kind, 'structure) t -> ('kind, 'structure) t
val lfun_get_bitprec : ('kind, 'structure) t -> Signed.long

val lfun :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val lfun0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfuncheckfeq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long -> Signed.long

val lfunconductor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfuncost :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfuncost0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfuncreate : ('kind, 'structure) t -> ('kind, 'structure) t
val lfundual : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val lfuneuler :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val lfunparams : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val lfunan :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val lfunhardy :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val lfuninit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfuninit0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfuninit_make :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val lfunlambda :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val lfunlambda0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfunmisc_to_ldata : ('kind, 'structure) t -> ('kind, 'structure) t
val lfunmisc_to_ldata_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val lfunmisc_to_ldata_shallow_i : ('kind, 'structure) t -> ('kind, 'structure) t

val lfunorderzero :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> Signed.long

val lfunprod_get_fact : ('kind, 'structure) t -> ('kind, 'structure) t
val lfunrootno : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val lfunrootres : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val lfunrtopoles : ('kind, 'structure) t -> ('kind, 'structure) t

val lfunshift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfuntwist :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val lfuntheta :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfunthetacost0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long

val lfunthetacost :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long

val lfunthetainit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfunthetacheckinit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfunzeros :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val sdomain_isincl :
  float -> ('kind, 'structure) t -> ('kind, 'structure) t -> int

val theta_get_an : ('kind, 'structure) t -> ('kind, 'structure) t
val theta_get_k : ('kind, 'structure) t -> ('kind, 'structure) t
val theta_get_r : ('kind, 'structure) t -> ('kind, 'structure) t
val theta_get_bitprec : ('kind, 'structure) t -> Signed.long
val theta_get_m : ('kind, 'structure) t -> Signed.long
val theta_get_tdom : ('kind, 'structure) t -> ('kind, 'structure) t
val theta_get_isqrtn : ('kind, 'structure) t -> ('kind, 'structure) t
val vgaeasytheta : ('kind, 'structure) t -> int

val znchargauss :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val dirzetak :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ellmoddegree : ('kind, 'structure) t -> ('kind, 'structure) t
val eta_zxn : Signed.long -> Signed.long -> ('kind, 'structure) t

val eta_product_zxn :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val etaquotype :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val galois_get_conj : ('kind, 'structure) t -> ('kind, 'structure) t

val ldata_vecan :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ldata_newprec :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val lfunabelianrelinit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfunartin :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val lfundiv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val lfunellmfpeters :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val lfunetaquo : ('kind, 'structure) t -> ('kind, 'structure) t
val lfungenus2 : ('kind, 'structure) t -> ('kind, 'structure) t
val lfunmfspec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val lfunmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val lfunqf : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val lfunsympow : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val lfunzetakinit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val qfiseven : ('kind, 'structure) t -> Signed.long
val lfunquadneg : Signed.long -> Signed.long -> ('kind, 'structure) t

val zm_lll_norms :
  ('kind, 'structure) t ->
  float ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val kerint : ('kind, 'structure) t -> ('kind, 'structure) t
val lll : ('kind, 'structure) t -> ('kind, 'structure) t

val lllfp :
  ('kind, 'structure) t -> float -> Signed.long -> ('kind, 'structure) t

val lllgen : ('kind, 'structure) t -> ('kind, 'structure) t
val lllgram : ('kind, 'structure) t -> ('kind, 'structure) t
val lllgramgen : ('kind, 'structure) t -> ('kind, 'structure) t
val lllgramint : ('kind, 'structure) t -> ('kind, 'structure) t
val lllgramkerim : ('kind, 'structure) t -> ('kind, 'structure) t
val lllgramkerimgen : ('kind, 'structure) t -> ('kind, 'structure) t
val lllint : ('kind, 'structure) t -> ('kind, 'structure) t
val lllintpartial : ('kind, 'structure) t -> ('kind, 'structure) t
val lllintpartial_inplace : ('kind, 'structure) t -> ('kind, 'structure) t
val lllkerim : ('kind, 'structure) t -> ('kind, 'structure) t
val lllkerimgen : ('kind, 'structure) t -> ('kind, 'structure) t
val matkerint0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val qflll0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val qflllgram0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtomap : ('kind, 'structure) t -> ('kind, 'structure) t
val mapdelete : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val mapdomain : ('kind, 'structure) t -> ('kind, 'structure) t
val mapdomain_shallow : ('kind, 'structure) t -> ('kind, 'structure) t

val mapget :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mapisdefined :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val mapput :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val maptomat : ('kind, 'structure) t -> ('kind, 'structure) t
val maptomat_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val matpermanent : ('kind, 'structure) t -> ('kind, 'structure) t
val zm_permanent : ('kind, 'structure) t -> ('kind, 'structure) t
val dbllemma526 : float -> float -> float -> float -> float
val dblcoro526 : float -> float -> float -> float

val gammamellininv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val gammamellininvasymp :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val gammamellininvinit :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val gammamellininvrt :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val member_a1 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_a2 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_a3 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_a4 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_a6 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_area : ('kind, 'structure) t -> ('kind, 'structure) t
val member_b2 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_b4 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_b6 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_b8 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_bid : ('kind, 'structure) t -> ('kind, 'structure) t
val member_bnf : ('kind, 'structure) t -> ('kind, 'structure) t
val member_c4 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_c6 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_clgp : ('kind, 'structure) t -> ('kind, 'structure) t
val member_codiff : ('kind, 'structure) t -> ('kind, 'structure) t
val member_cyc : ('kind, 'structure) t -> ('kind, 'structure) t
val member_diff : ('kind, 'structure) t -> ('kind, 'structure) t
val member_disc : ('kind, 'structure) t -> ('kind, 'structure) t
val member_e : ('kind, 'structure) t -> ('kind, 'structure) t
val member_eta : ('kind, 'structure) t -> ('kind, 'structure) t
val member_f : ('kind, 'structure) t -> ('kind, 'structure) t
val member_fu : ('kind, 'structure) t -> ('kind, 'structure) t
val member_gen : ('kind, 'structure) t -> ('kind, 'structure) t
val member_group : ('kind, 'structure) t -> ('kind, 'structure) t
val member_index : ('kind, 'structure) t -> ('kind, 'structure) t
val member_j : ('kind, 'structure) t -> ('kind, 'structure) t
val member_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val member_nf : ('kind, 'structure) t -> ('kind, 'structure) t
val member_no : ('kind, 'structure) t -> ('kind, 'structure) t
val member_omega : ('kind, 'structure) t -> ('kind, 'structure) t
val member_orders : ('kind, 'structure) t -> ('kind, 'structure) t
val member_p : ('kind, 'structure) t -> ('kind, 'structure) t
val member_pol : ('kind, 'structure) t -> ('kind, 'structure) t
val member_polabs : ('kind, 'structure) t -> ('kind, 'structure) t
val member_reg : ('kind, 'structure) t -> ('kind, 'structure) t
val member_r1 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_r2 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_roots : ('kind, 'structure) t -> ('kind, 'structure) t
val member_sign : ('kind, 'structure) t -> ('kind, 'structure) t
val member_t2 : ('kind, 'structure) t -> ('kind, 'structure) t
val member_tate : ('kind, 'structure) t -> ('kind, 'structure) t
val member_tu : ('kind, 'structure) t -> ('kind, 'structure) t
val member_zk : ('kind, 'structure) t -> ('kind, 'structure) t
val member_zkst : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_m : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_mindex : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_minv : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_basis : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_dim : ('kind, 'structure) t -> Signed.long
val mf_get_e : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_fields : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_newforms : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_space : ('kind, 'structure) t -> Signed.long
val mf_get_s : ('kind, 'structure) t -> ('kind, 'structure) t
val mfcusp_get_vmjd : ('kind, 'structure) t -> ('kind, 'structure) t
val mfnew_get_vj : ('kind, 'structure) t -> ('kind, 'structure) t

val qab_tracerel :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qabm_tracerel :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qabv_tracerel :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qab_trace_init :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val checkmf : ('kind, 'structure) t -> ('kind, 'structure) t
val checkmf_i : ('kind, 'structure) t -> int
val getcache : unit -> ('kind, 'structure) t
val hclassno6u : pari_ulong -> pari_ulong
val hclassno6u_no_cache : pari_ulong -> pari_ulong

val lfunmf :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfdelta : unit -> ('kind, 'structure) t
val mfeh : ('kind, 'structure) t -> ('kind, 'structure) t
val mfek : Signed.long -> ('kind, 'structure) t
val mftheta : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_chi : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_n : ('kind, 'structure) t -> Signed.long
val mf_get_nk : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_field : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_gn : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_gk : ('kind, 'structure) t -> ('kind, 'structure) t
val mf_get_k : ('kind, 'structure) t -> Signed.long
val mf_get_r : ('kind, 'structure) t -> Signed.long
val mf_get_type : ('kind, 'structure) t -> Signed.long

val mfatkin :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfatkineigenvalues :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mfatkininit :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mfbasis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mfbd : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mfbracket :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfcharorder : ('kind, 'structure) t -> Signed.long
val mfcharmodulus : ('kind, 'structure) t -> Signed.long
val mfcharpol : ('kind, 'structure) t -> ('kind, 'structure) t
val mfcoef : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mfcoefs :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mfconductor : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long
val mfcosets : ('kind, 'structure) t -> ('kind, 'structure) t

val mfcuspdim :
  Signed.long -> Signed.long -> ('kind, 'structure) t -> Signed.long

val mfcuspisregular :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val mfcusps : ('kind, 'structure) t -> ('kind, 'structure) t

val mfcuspval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfcuspwidth : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long
val mfderiv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mfderive2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mfdescribe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val mfdim : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mfdiv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfdiv_val :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfeigenbasis : ('kind, 'structure) t -> ('kind, 'structure) t

val mfeigensearch :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfeisenstein :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mfeisensteindim :
  Signed.long -> Signed.long -> ('kind, 'structure) t -> Signed.long

val mfembed :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfembed0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfeval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mffields : ('kind, 'structure) t -> ('kind, 'structure) t
val mffromell : ('kind, 'structure) t -> ('kind, 'structure) t
val mffrometaquo : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mffromlfun : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mffromqf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mffulldim :
  Signed.long -> Signed.long -> ('kind, 'structure) t -> Signed.long

val mfgaloisprojrep :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfgaloistype :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfhecke :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfheckemat :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfinit : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mfiscm : ('kind, 'structure) t -> ('kind, 'structure) t
val mfiscuspidal : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val mfisequal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long -> Signed.long

val mfisetaquo : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mfkohnenbasis : ('kind, 'structure) t -> ('kind, 'structure) t
val mfkohnenbijection : ('kind, 'structure) t -> ('kind, 'structure) t

val mfkohneneigenbasis :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mflinear :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfmanin : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mfmatembed :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfmul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfnewdim :
  Signed.long -> Signed.long -> ('kind, 'structure) t -> Signed.long

val mfolddim :
  Signed.long -> Signed.long -> ('kind, 'structure) t -> Signed.long

val mfparams : ('kind, 'structure) t -> ('kind, 'structure) t

val mfperiodpol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val mfperiodpolbasis : Signed.long -> Signed.long -> ('kind, 'structure) t

val mfpetersson :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfpow : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mfsearch :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfshift : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mfshimura :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfslashexpansion :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t

val mfspace : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val mfsplit :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mfsturm : ('kind, 'structure) t -> Signed.long
val mfsturmngk : Signed.long -> ('kind, 'structure) t -> Signed.long
val mfsturmnk : Signed.long -> Signed.long -> Signed.long
val mfsturm_mf : ('kind, 'structure) t -> Signed.long

val mfsymboleval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mfsymbol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mftaylor :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val mftobasis :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mftobasises :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mftocol :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mftocoset :
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mftonew :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mftraceform : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mftwist :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfvecembed :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mfvectomat :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val fl_inv : pari_ulong -> pari_ulong -> pari_ulong
val fl_invsafe : pari_ulong -> pari_ulong -> pari_ulong

val fp_ratlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val zm2_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val abscmpii : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val abscmprr : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val absequalii : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val addii_sign :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val addir_sign :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val addmulii :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val addmulii_inplace :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val addrr_sign :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val addsi_sign :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val addsr : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val addui_sign :
  pari_ulong -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val addumului :
  pari_ulong -> pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t

val affir : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val affrr : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val bezout :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val cbezout :
  Signed.long ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val cgcd : Signed.long -> Signed.long -> Signed.long
val clcm : Signed.long -> Signed.long -> Signed.long
val cmpii : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val cmprr : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val dblexpo : float -> Signed.long
val dblmantissa : float -> pari_ulong
val dbltor : float -> ('kind, 'structure) t

val diviiexact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val divir :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val divis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val divis_rem :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val absdiviu_rem :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  ('kind, 'structure) t

val diviuuexact :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val diviuexact : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val divri :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val divrr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val divrs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val divru : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val divsi : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val divsr : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val divur : pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t

val dvmdii :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val equalii : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val equalrr : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val floorr : ('kind, 'structure) t -> ('kind, 'structure) t

val gcdii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val halfgcdii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val int2n : Signed.long -> ('kind, 'structure) t
val int2u : pari_ulong -> ('kind, 'structure) t
val int2um1 : pari_ulong -> ('kind, 'structure) t

val int_normalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val invmod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val invmod2bil : pari_ulong -> pari_ulong
val invr : ('kind, 'structure) t -> ('kind, 'structure) t

val mantissa_real :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val modii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val modiiz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val mulii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mulir :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mulrr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mulsi : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val mulsr : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val mulss : Signed.long -> Signed.long -> ('kind, 'structure) t
val mului : pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t
val mulur : pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t
val muluu : pari_ulong -> pari_ulong -> ('kind, 'structure) t

val muluui :
  pari_ulong -> pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t

val pari_kernel_close : unit -> unit
val pari_kernel_init : unit -> unit
val pari_kernel_version : unit -> string
val remi2n : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rtodbl : ('kind, 'structure) t -> float
val shifti : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val sqri : ('kind, 'structure) t -> ('kind, 'structure) t
val sqrr : ('kind, 'structure) t -> ('kind, 'structure) t
val sqrs : Signed.long -> ('kind, 'structure) t
val sqrtr_abs : ('kind, 'structure) t -> ('kind, 'structure) t

val sqrtremi :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val sqru : pari_ulong -> ('kind, 'structure) t
val subsr : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val truedvmdii :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val truedvmdis :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val truedvmdsi :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val trunc2nr : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mantissa2nr : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val truncr : ('kind, 'structure) t -> ('kind, 'structure) t
val ugcd : pari_ulong -> pari_ulong -> pari_ulong
val ulcm : pari_ulong -> pari_ulong -> pari_ulong
val umodiu : ('kind, 'structure) t -> pari_ulong -> pari_ulong
val vals : pari_ulong -> Signed.long

val fpc_ratlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpm_ratlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_ratlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxqx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zxqx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nffactor :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nffactormod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfgcd_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val nfissquarefree : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val nfroots :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfroots_if_split :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfrootsof1 : ('kind, 'structure) t -> ('kind, 'structure) t

val polfnf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfabelianconjgen :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfisabelian : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val forpart :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val forpart_init :
  forpart_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val forpart_next :
  forpart_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val forpart_prev :
  forpart_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val numbpart : ('kind, 'structure) t -> ('kind, 'structure) t

val partitions :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val forperm :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  unit

val forperm_init :
  forperm_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t -> unit

val forperm_next :
  forperm_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val forallsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> Signed.long -> unit

val forksubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  unit

val forsubset_next :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val forsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  unit

val glambertw :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mplambertw : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mplambertx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mplambertx_logx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mplambertxlogx_x :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val z_to_perm : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val abelian_group : ('kind, 'structure) t -> ('kind, 'structure) t

val conjclasses_repr :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val cyc_pow : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val cyc_pow_perm : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val cyclicgroup : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val dicyclicgroup :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val group_abelianhnf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val group_abeliansnf :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val group_domain : ('kind, 'structure) t -> Signed.long
val group_elts : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val group_export : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val group_export_gap : ('kind, 'structure) t -> ('kind, 'structure) t
val group_export_magma : ('kind, 'structure) t -> ('kind, 'structure) t
val group_isa4s4 : ('kind, 'structure) t -> Signed.long
val group_isabelian : ('kind, 'structure) t -> Signed.long

val group_leftcoset :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val group_order : ('kind, 'structure) t -> Signed.long

val group_perm_normalize :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val group_quotient :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val group_rightcoset :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val group_set : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val group_subgroup_is_faithful :
  ('kind, 'structure) t -> ('kind, 'structure) t -> int

val group_subgroup_isnormal :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val group_subgroups : ('kind, 'structure) t -> ('kind, 'structure) t
val groupelts_solvablesubgroups : ('kind, 'structure) t -> ('kind, 'structure) t
val group_to_cc : ('kind, 'structure) t -> ('kind, 'structure) t
val groupelts_abelian_group : ('kind, 'structure) t -> ('kind, 'structure) t
val groupelts_center : ('kind, 'structure) t -> ('kind, 'structure) t

val groupelts_conj_set :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val groupelts_conjclasses :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val groupelts_exponent : ('kind, 'structure) t -> Signed.long

val groupelts_quotient :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val groupelts_set :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val groupelts_to_group : ('kind, 'structure) t -> ('kind, 'structure) t
val numtoperm : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val perm_commute : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val perm_cycles : ('kind, 'structure) t -> ('kind, 'structure) t
val perm_order : ('kind, 'structure) t -> ('kind, 'structure) t
val perm_orderu : ('kind, 'structure) t -> pari_ulong

val perm_pow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val perm_powu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val perm_sign : ('kind, 'structure) t -> Signed.long
val perm_to_gap : ('kind, 'structure) t -> ('kind, 'structure) t
val perm_to_z : ('kind, 'structure) t -> ('kind, 'structure) t
val permcycles : ('kind, 'structure) t -> ('kind, 'structure) t
val permorder : ('kind, 'structure) t -> ('kind, 'structure) t
val permsign : ('kind, 'structure) t -> Signed.long
val permtonum : ('kind, 'structure) t -> ('kind, 'structure) t

val quotient_group :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val quotient_groupelts : ('kind, 'structure) t -> ('kind, 'structure) t

val quotient_perm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val quotient_subgroup_lift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val subgroups_tableset :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val tableset_find_index :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val trivialgroup : unit -> ('kind, 'structure) t

val vec_insert :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val vec_is1to1 : ('kind, 'structure) t -> int
val vec_isconst : ('kind, 'structure) t -> int

val vecperm_orbits :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val vecsmall_duplicate : ('kind, 'structure) t -> Signed.long
val vecsmall_duplicate_sorted : ('kind, 'structure) t -> Signed.long
val vecsmall_indexsort : ('kind, 'structure) t -> ('kind, 'structure) t
val vecsmall_is1to1 : ('kind, 'structure) t -> int
val vecsmall_isconst : ('kind, 'structure) t -> int
val vecsmall_sort : ('kind, 'structure) t -> unit
val vecsmall_uniq : ('kind, 'structure) t -> ('kind, 'structure) t
val vecsmall_uniq_sorted : ('kind, 'structure) t -> ('kind, 'structure) t

val vecsmall_counting_indexsort :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val vecsmall_counting_sort : ('kind, 'structure) t -> Signed.long -> unit

val vecsmall_counting_uniq :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val vecvecsmall_indexsort : ('kind, 'structure) t -> ('kind, 'structure) t
val vecvecsmall_max : ('kind, 'structure) t -> Signed.long

val vecvecsmall_search :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val vecvecsmall_sort : ('kind, 'structure) t -> ('kind, 'structure) t

val vecvecsmall_sort_inplace :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> unit

val vecvecsmall_sort_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val vecvecsmall_sort_uniq : ('kind, 'structure) t -> ('kind, 'structure) t
val mt_broadcast : ('kind, 'structure) t -> unit
val mt_nbthreads : unit -> Signed.long
val mt_queue_end : pari_mt Ctypes.structure Ctypes_static.ptr -> unit

val mt_queue_get :
  pari_mt Ctypes.structure Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val mt_queue_start :
  pari_mt Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t -> unit

val mt_queue_start_lim :
  pari_mt Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val mt_queue_submit :
  pari_mt Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t ->
  unit

val mt_sigint_block : unit -> unit
val mt_sigint_unblock : unit -> unit
val pari_mt_init : unit -> unit
val pari_mt_close : unit -> unit

val subcyclopclgp :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val subcycloiwasawa :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val subcyclohminus :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val color_to_rgb :
  ('kind, 'structure) t ->
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
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val parplothexport :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val plotbox :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val plotclip : Signed.long -> unit
val plotcolor : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val plotcopy :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val plotcursor : Signed.long -> ('kind, 'structure) t
val plotdraw : ('kind, 'structure) t -> Signed.long -> unit

val plotexport :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val plothexport :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val plothraw :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val plothrawexport :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val plothsizes : Signed.long -> ('kind, 'structure) t

val plotinit :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val plotkill : Signed.long -> unit

val plotline :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val plotlines :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val plotlinetype : Signed.long -> Signed.long -> unit

val plotmove :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val plotpoints :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val plotpointsize : Signed.long -> ('kind, 'structure) t -> unit
val plotpointtype : Signed.long -> Signed.long -> unit

val plotrbox :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val plotrecth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val plotrecthraw :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val plotrline :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val plotrmove :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val plotrpoint :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val plotscale :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val plotstring : Signed.long -> string -> Signed.long -> unit
val psdraw : ('kind, 'structure) t -> Signed.long -> unit

val psploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val psplothraw :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rect2ps :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  string

val rect2ps_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  int ->
  string

val rect2svg :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  string

val pariplot :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  unit

val zx_zp_root :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zp_appr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val cmp_padic : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val factorpadic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gdeuc :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val grem :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val padicappr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val poldivrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val polrootspadic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val flv_factorback :
  ('kind, 'structure) t -> ('kind, 'structure) t -> pari_ulong -> pari_ulong

val flxqv_factorback :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val fpv_factorback :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqv_factorback :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val q_content : ('kind, 'structure) t -> ('kind, 'structure) t
val q_content_safe : ('kind, 'structure) t -> ('kind, 'structure) t
val q_denom : ('kind, 'structure) t -> ('kind, 'structure) t
val q_denom_safe : ('kind, 'structure) t -> ('kind, 'structure) t

val q_div_to_int :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val q_gcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val q_mul_to_int :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val q_muli_to_int :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val q_primitive_part :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val q_primpart : ('kind, 'structure) t -> ('kind, 'structure) t

val q_remove_denom :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val q_factor : ('kind, 'structure) t -> ('kind, 'structure) t

val q_factor_limit :
  ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val rg_type :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rgc_type :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rescale_to_int : ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_type :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_type2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_rg_type :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_chinese_coprime :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgx_disc : ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_extgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgx_extgcd_simple :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgx_gcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_gcd_simple :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_halfgcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_halfgcd_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgx_rescale_to_int : ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_resultant_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val rgx_sturmpart :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val rgx_sylvestermatrix :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_type :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type3 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type_decode :
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val rgx_type_is_composite : Signed.long -> int

val rgxq_charpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgxq_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxq_minpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rgxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxq_ratlift :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val rgxq_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val z_content : ('kind, 'structure) t -> ('kind, 'structure) t
val zx_content : ('kind, 'structure) t -> ('kind, 'structure) t

val centermod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val centermod_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val centermodii :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val content : ('kind, 'structure) t -> ('kind, 'structure) t

val content0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val deg1_from_roots :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val factor : ('kind, 'structure) t -> ('kind, 'structure) t

val factor0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val factorback : ('kind, 'structure) t -> ('kind, 'structure) t

val factorback2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gbezout :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val gdivexact :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gen_factorback :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t) Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val ggcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ggcd0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ghalfgcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ginvmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val glcm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val glcm0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val newtonpoly :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfrootsq : ('kind, 'structure) t -> ('kind, 'structure) t
val poldisc0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val polisirreducible : ('kind, 'structure) t -> Signed.long

val polresultant0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val polsym : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val primitive_part :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val primpart : ('kind, 'structure) t -> ('kind, 'structure) t
val reduceddiscsmith : ('kind, 'structure) t -> ('kind, 'structure) t

val resultant2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val resultant :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rnfcharpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val roots_from_deg1 : ('kind, 'structure) t -> ('kind, 'structure) t
val roots_to_pol : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val roots_to_pol_r1 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val sturmpart :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long

val subresext :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val sylvestermatrix :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val trivial_fact : unit -> ('kind, 'structure) t

val gcdext0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val polresultantext0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val polresultantext :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val prime_fact : ('kind, 'structure) t -> ('kind, 'structure) t
val row_q_primpart : ('kind, 'structure) t -> ('kind, 'structure) t
val vec_q_primpart : ('kind, 'structure) t -> ('kind, 'structure) t
val vecprod : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_lcm : ('kind, 'structure) t -> ('kind, 'structure) t

val flx_flxy_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxx_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) t

val fpx_fpxy_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_translate :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxqx_normalize :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxv_fpc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxy_fpxq_evaly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fpxc_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxm_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_fp_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_halve :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_invsafe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_mulu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_neg :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_neg_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_sqrt :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_sqrtn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fq_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqc_fq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqc_fqv_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqc_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqc_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqv_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqv_roots_to_pol :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fqx_fq_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_fq_mul_to_monic :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_fq_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_translate :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_matrix_pow :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_powers :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxy_eval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxy_evalx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qx_disc : ('kind, 'structure) t -> ('kind, 'structure) t

val qx_gcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qx_resultant :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxq_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxq_intnorm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxq_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qxq_norm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qxq_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rg_is_fp :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> int

val rg_is_fpxq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val rg_to_fp :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rg_to_fpxq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgc_to_fpc :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgc_to_fqc :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgm_is_fpm :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> int

val rgm_to_flm : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val rgm_to_fpm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgm_to_fqm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgv_is_fpv :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> int

val rgv_to_flv : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val rgv_to_fpv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_is_fpx :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> int

val rgx_to_fpx :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_is_fpxqx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val rgx_to_fpxqx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgx_to_fqx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val z_incremental_crt :
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong ->
  int

val z_init_crt : pari_ulong -> pari_ulong -> ('kind, 'structure) t

val zm_incremental_crt :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong ->
  int

val zm_init_crt : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val zx_zxy_resultant :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_zxy_rnfequation :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val zx_disc : ('kind, 'structure) t -> ('kind, 'structure) t

val zx_gcd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zx_gcd_all :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zx_incremental_crt :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong ->
  int

val zx_init_crt :
  ('kind, 'structure) t -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val zx_is_squarefree : ('kind, 'structure) t -> int
val zx_radical : ('kind, 'structure) t -> ('kind, 'structure) t

val zx_resultant :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zxm_incremental_crt :
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  pari_ulong ->
  int

val zxm_init_crt :
  ('kind, 'structure) t -> Signed.long -> pari_ulong -> ('kind, 'structure) t

val zxq_minpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  pari_ulong ->
  ('kind, 'structure) t

val zxq_charpoly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val characteristic : ('kind, 'structure) t -> ('kind, 'structure) t
val ffnbirred : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ffnbirred0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val ffsumnbirred : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val get_fq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  bb_field Ctypes.structure Ctypes_static.ptr

val init_flxq :
  pari_ulong -> Signed.long -> Signed.long -> ('kind, 'structure) t

val init_fq :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val nfx_disc :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val nfx_resultant :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val pol_x_powers : Signed.long -> Signed.long -> ('kind, 'structure) t
val residual_characteristic : ('kind, 'structure) t -> ('kind, 'structure) t

val polclass :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val fp_modinv_to_j :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_polmodular_evalx :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  int ->
  ('kind, 'structure) t

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
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val polmodular_zm : Signed.long -> Signed.long -> ('kind, 'structure) t

val polmodular_zxx :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val bpsw_isprime : ('kind, 'structure) t -> Signed.long
val bpsw_psp : ('kind, 'structure) t -> Signed.long
val addprimes : ('kind, 'structure) t -> ('kind, 'structure) t
val check_ecppcert : ('kind, 'structure) t -> Signed.long
val gisprime : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gispseudoprime :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gprimepi_upper_bound : ('kind, 'structure) t -> ('kind, 'structure) t
val gprimepi_lower_bound : ('kind, 'structure) t -> ('kind, 'structure) t
val isprime : ('kind, 'structure) t -> Signed.long
val ispseudoprime : ('kind, 'structure) t -> Signed.long -> Signed.long
val millerrabin : ('kind, 'structure) t -> Signed.long -> Signed.long
val prime : Signed.long -> ('kind, 'structure) t
val primecert : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val primecert0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val primecertexport :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val primecertisvalid : ('kind, 'structure) t -> Signed.long
val primepi : ('kind, 'structure) t -> ('kind, 'structure) t
val primepi_upper_bound : float -> float
val primepi_lower_bound : float -> float
val primes : Signed.long -> ('kind, 'structure) t

val primes_interval :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val primes_interval_zv : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val primes_upto_zv : pari_ulong -> ('kind, 'structure) t
val primes0 : ('kind, 'structure) t -> ('kind, 'structure) t
val primes_zv : Signed.long -> ('kind, 'structure) t
val randomprime : ('kind, 'structure) t -> ('kind, 'structure) t

val randomprime0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val removeprimes : ('kind, 'structure) t -> ('kind, 'structure) t
val uis2psp : pari_ulong -> int
val uispsp : pari_ulong -> pari_ulong -> int
val uislucaspsp : pari_ulong -> int
val uisprime : pari_ulong -> int
val uisprime_101 : pari_ulong -> int
val uisprime_661 : pari_ulong -> int
val uprime : Signed.long -> pari_ulong
val uprimepi : pari_ulong -> pari_ulong

val qfauto :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfauto0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfautoexport : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val qfisom :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfisom0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfisominit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfisominit0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qforbits :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qfminimize : ('kind, 'structure) t -> ('kind, 'structure) t

val qfparam :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qfsolve : ('kind, 'structure) t -> ('kind, 'structure) t
val z_isfundamental : ('kind, 'structure) t -> Signed.long
val classno : ('kind, 'structure) t -> ('kind, 'structure) t
val classno2 : ('kind, 'structure) t -> ('kind, 'structure) t

val hclassnof_fact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val hclassno : ('kind, 'structure) t -> ('kind, 'structure) t
val hclassno6 : ('kind, 'structure) t -> ('kind, 'structure) t
val isfundamental : ('kind, 'structure) t -> Signed.long
val qfb_equal1 : ('kind, 'structure) t -> int
val qfbclassno0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val qfi_shanks :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qfi_log :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfi_order :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val quadclassnof :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val quadclassnof_fact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val quaddisc : ('kind, 'structure) t -> ('kind, 'structure) t

val quadregulator :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val quadunit : ('kind, 'structure) t -> ('kind, 'structure) t
val quadunit0 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val quadunitindex :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val quadunitnorm : ('kind, 'structure) t -> Signed.long
val sisfundamental : Signed.long -> Signed.long
val uhclassnof_fact : ('kind, 'structure) t -> Signed.long -> Signed.long
val unegisfundamental : pari_ulong -> Signed.long
val unegquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val uposisfundamental : pari_ulong -> Signed.long
val uposquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val uquadclassnof_fact :
  pari_ulong ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong

val zn_quad_roots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val genrand : ('kind, 'structure) t -> ('kind, 'structure) t
val getrand : unit -> ('kind, 'structure) t
val pari_rand : unit -> pari_ulong
val randomi : ('kind, 'structure) t -> ('kind, 'structure) t
val randomr : Signed.long -> ('kind, 'structure) t
val random_f2x : Signed.long -> Signed.long -> ('kind, 'structure) t
val random_fl : pari_ulong -> pari_ulong
val random_bits : Signed.long -> Signed.long
val random_zv : Signed.long -> ('kind, 'structure) t
val setrand : ('kind, 'structure) t -> unit

val ellratpoints :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val hyperellratpoints :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qx_complex_roots :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val fft :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fftinv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val cleanroots : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val fujiwara_bound : ('kind, 'structure) t -> float
val fujiwara_bound_real : ('kind, 'structure) t -> Signed.long -> float
val isrealappr : ('kind, 'structure) t -> Signed.long -> int
val polgraeffe : ('kind, 'structure) t -> ('kind, 'structure) t

val polmod_to_embed :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val polrootsbound :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val roots : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val realroots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zx_graeffe : ('kind, 'structure) t -> ('kind, 'structure) t

val zx_realroots_irred :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zx_sturm : ('kind, 'structure) t -> Signed.long
val zx_sturm_irred : ('kind, 'structure) t -> Signed.long
val zx_sturmpart : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val zx_uspensky :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val factor_aurifeuille :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val factor_aurifeuille_prime :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val galoissubcyclo :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val polsubcyclo :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val polsubcyclofast :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val znsubgroupgenerators :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val nfsubfields : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val nfsubfields0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val nfsubfieldscm :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val nfsubfieldsmax :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val nflist :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val nfresolvent : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val subgrouplist :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val forsubgroup :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val abmap_kernel : ('kind, 'structure) t -> ('kind, 'structure) t

val abmap_subgroup_image :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnrl1 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val bnrrootnumber :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val bnrstark :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val cyc2elts : ('kind, 'structure) t -> ('kind, 'structure) t
val qfbforms : ('kind, 'structure) t -> ('kind, 'structure) t
val quadhilbert : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val quadray :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val chartogenstr : char -> ('kind, 'structure) t
val pari_strdup : string -> string
val pari_strndup : string -> Signed.long -> string
val stack_strcat : string -> string -> string
val stack_strdup : string -> string
val pari_strchr : ('kind, 'structure) t -> ('kind, 'structure) t

val strjoin :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val strntogenstr : string -> Signed.long -> ('kind, 'structure) t

val strsplit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val strtogenstr : string -> ('kind, 'structure) t
val type_name : Signed.long -> string
val type0 : ('kind, 'structure) t -> ('kind, 'structure) t

val asympnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val asympnumraw :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val derivnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val derivnumk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val derivfun :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val derivfunk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val forvec_init :
  forvec_t Ctypes.structure Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  int

val forvec_next :
  forvec_t Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val laurentseries :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val limitnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val polzag : Signed.long -> Signed.long -> ('kind, 'structure) t

val prodeuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val prodinf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val prodinf1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val solvestep :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val sumalt :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumalt2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumpos :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumpos2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val suminf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val suminf_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sumdivmultexpr :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zbrent :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val bnfisintnorm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnfisintnormabs :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val ideals_by_norm :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val thue :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val thueinit :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val pi2n : Signed.long -> Signed.long -> ('kind, 'structure) t
val pii2 : Signed.long -> ('kind, 'structure) t
val pii2n : Signed.long -> Signed.long -> ('kind, 'structure) t
val qp_exp : ('kind, 'structure) t -> ('kind, 'structure) t
val qp_exp_prec : ('kind, 'structure) t -> Signed.long
val qp_log : ('kind, 'structure) t -> ('kind, 'structure) t
val qp_sqrt : ('kind, 'structure) t -> ('kind, 'structure) t

val qp_sqrtn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val zn_sqrt :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zp_teichmuller :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val agm :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val constcatalan : Signed.long -> ('kind, 'structure) t
val consteuler : Signed.long -> ('kind, 'structure) t
val constlog2 : Signed.long -> ('kind, 'structure) t
val constpi : Signed.long -> ('kind, 'structure) t
val cxexpm1 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val elle : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ellk : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val expir : ('kind, 'structure) t -> ('kind, 'structure) t
val exp1r_abs : ('kind, 'structure) t -> ('kind, 'structure) t
val gcos : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gcotan : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gcotanh : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gexp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gexpm1 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val glog : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val glog1p : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gpow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val gpowers : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gpowers0 :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gpowgs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val grootsof1 : Signed.long -> Signed.long -> ('kind, 'structure) t
val gsin : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gsinc : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gsincos :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  unit

val gsqrpowers : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gsqrt : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gsqrtn :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  Signed.long ->
  ('kind, 'structure) t

val gtan : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val logr_abs : ('kind, 'structure) t -> ('kind, 'structure) t
val mpcos : ('kind, 'structure) t -> ('kind, 'structure) t
val mpeuler : Signed.long -> ('kind, 'structure) t
val mpcatalan : Signed.long -> ('kind, 'structure) t

val mpsincosm1 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val mpexp : ('kind, 'structure) t -> ('kind, 'structure) t
val mpexpm1 : ('kind, 'structure) t -> ('kind, 'structure) t
val mplog : ('kind, 'structure) t -> ('kind, 'structure) t
val mplog2 : Signed.long -> ('kind, 'structure) t
val mppi : Signed.long -> ('kind, 'structure) t
val mpsin : ('kind, 'structure) t -> ('kind, 'structure) t

val mpsincos :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val pow2pis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val powpis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val powcx :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val powcx_prec :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> Signed.long

val powersr : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val powiu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val powrfrac :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val powrs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val powrshalf : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val powru : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val powruhalf : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val powuu : pari_ulong -> pari_ulong -> ('kind, 'structure) t

val powgi :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rootsof1_cx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rootsof1u_cx : pari_ulong -> Signed.long -> ('kind, 'structure) t

val rootsof1q_cx :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rootsof1powinit :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rootsof1pow : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val serchop : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val serchop_i : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val serchop0 : ('kind, 'structure) t -> ('kind, 'structure) t
val sqrtnint : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val sqrtnr_abs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val teich : ('kind, 'structure) t -> ('kind, 'structure) t
val teichmullerinit : Signed.long -> Signed.long -> ('kind, 'structure) t

val teichmuller :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val trans_eval :
  string ->
  (('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val trans_evalgen :
  string ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val upowuu : pari_ulong -> pari_ulong -> pari_ulong
val upowers : pari_ulong -> Signed.long -> ('kind, 'structure) t
val usqrtn : pari_ulong -> pari_ulong -> pari_ulong
val usqrt : pari_ulong -> pari_ulong
val qp_gamma : ('kind, 'structure) t -> ('kind, 'structure) t
val atanhuu : pari_ulong -> pari_ulong -> Signed.long -> ('kind, 'structure) t

val atanhui :
  pari_ulong -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gacosh : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gacos : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val garg : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gasinh : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gasin : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gatan : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gatanh : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gcosh : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ggammah : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ggamma : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val ggamma1m1 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val glngamma : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gpsi : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gsinh : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtanh : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mpfactr : Signed.long -> Signed.long -> ('kind, 'structure) t

val mpsinhcosh :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val psi1series :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val sumformal : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgv_is_arithprog :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  int

val besseljzero :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val besselyzero :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val constzeta : Signed.long -> Signed.long -> ('kind, 'structure) t

val cxek :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val dblmodulus : ('kind, 'structure) t -> float
val dilog : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val eint1 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val expipir : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val expipic : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val expixy :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val eta : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val eta0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val gerfc : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gpolylog :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val gzeta : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val hbessel1 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val hbessel2 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val hyperu :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val ibessel :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val incgam :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val incgam0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val incgamc :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val jbessel :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val jbesselh :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val jell : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val kbessel :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mpeint1 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mpveceint1 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val polylog0 :
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val sumdedekind :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val sumdedekind_coprime :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val szeta : Signed.long -> Signed.long -> ('kind, 'structure) t

val theta :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val thetanullk :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val trueeta : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val u_sumdedekind_coprime : Signed.long -> Signed.long -> ('kind, 'structure) t

val upper_to_cx :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val veceint1 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val vecthetanullk :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val vecthetanullk_tau :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val veczeta :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val weber0 :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val weberf : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val weberf1 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val weberf2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val ybessel :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val sl2_inv_shallow : ('kind, 'structure) t -> ('kind, 'structure) t

val qevproj_apply :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qevproj_apply_vecei :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val qevproj_down :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val qevproj_init : ('kind, 'structure) t -> ('kind, 'structure) t
val rgx_act_gl2q : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgx_act_zgl2q :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val checkms : ('kind, 'structure) t -> unit
val checkmspadic : ('kind, 'structure) t -> unit

val ellpadiclambdamu :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val mfnumcusps : ('kind, 'structure) t -> ('kind, 'structure) t
val mfnumcusps_fact : ('kind, 'structure) t -> ('kind, 'structure) t
val mfnumcuspsu_fact : ('kind, 'structure) t -> pari_ulong
val mfnumcuspsu : pari_ulong -> pari_ulong

val msfromcusp :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val msfromell : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val msfromhecke :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val msdim : ('kind, 'structure) t -> Signed.long

val mseval2_ooq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val msgetlevel : ('kind, 'structure) t -> Signed.long
val msgetsign : ('kind, 'structure) t -> Signed.long
val msgetweight : ('kind, 'structure) t -> Signed.long

val msatkinlehner :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mscuspidal : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mseisenstein : ('kind, 'structure) t -> ('kind, 'structure) t

val mseval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mshecke :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val msinit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val msissymbol :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mslattice :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val msomseval :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mspadic_parse_chi :
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  unit

val mspadic_unit_eigenvalue :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mspadicinit :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val mspadicl :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mspadicmoments :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mspadicseries :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val mspathgens : ('kind, 'structure) t -> ('kind, 'structure) t

val mspathlog :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val msnew : ('kind, 'structure) t -> ('kind, 'structure) t

val mspetersson :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mspolygon : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val msstar :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val msqexpansion :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val mssplit :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val mstooms :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mscosets0 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mscosets :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t

val msfarey :
  ('kind, 'structure) t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> ('kind, 'structure) t -> Signed.long)
  Ctypes_static.static_funptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val msfarey0 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val checkfarey_i : ('kind, 'structure) t -> int

val polylogmult :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val polylogmult_interpolate :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zetamult : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val zetamultdual : ('kind, 'structure) t -> ('kind, 'structure) t

val zetamult_interpolate :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val zetamultall :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val zetamultconvert :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

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
val abscmpiu : ('kind, 'structure) t -> pari_ulong -> int
val abscmpui : pari_ulong -> ('kind, 'structure) t -> int
val absequaliu : ('kind, 'structure) t -> pari_ulong -> int
val absi : ('kind, 'structure) t -> ('kind, 'structure) t
val absi_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val absr : ('kind, 'structure) t -> ('kind, 'structure) t
val absrnz_equal1 : ('kind, 'structure) t -> int
val absrnz_equal2n : ('kind, 'structure) t -> int

val addii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val addiiz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val addir :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val addirz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val addis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val addri :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val addriz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val addrr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val addrrz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val addrs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val addsi : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val addsiz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val addsrz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val addss : Signed.long -> Signed.long -> ('kind, 'structure) t
val addssz : Signed.long -> Signed.long -> ('kind, 'structure) t -> unit
val adduu : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val affii : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val affiz : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val affrr_fixlg : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val affsi : Signed.long -> ('kind, 'structure) t -> unit
val affsr : Signed.long -> ('kind, 'structure) t -> unit
val affsz : Signed.long -> ('kind, 'structure) t -> unit
val affui : pari_ulong -> ('kind, 'structure) t -> unit
val affur : pari_ulong -> ('kind, 'structure) t -> unit
val cgetg : Signed.long -> Signed.long -> ('kind, 'structure) t
val cgetg_block : Signed.long -> Signed.long -> ('kind, 'structure) t

val cgetg_copy :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val cgeti : Signed.long -> ('kind, 'structure) t
val cgetineg : Signed.long -> ('kind, 'structure) t
val cgetipos : Signed.long -> ('kind, 'structure) t
val cgetr : Signed.long -> ('kind, 'structure) t
val cgetr_block : Signed.long -> ('kind, 'structure) t
val cmpir : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val cmpis : ('kind, 'structure) t -> Signed.long -> int
val cmpiu : ('kind, 'structure) t -> pari_ulong -> int
val cmpri : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val cmprs : ('kind, 'structure) t -> Signed.long -> int
val cmpsi : Signed.long -> ('kind, 'structure) t -> int
val cmpsr : Signed.long -> ('kind, 'structure) t -> int
val cmpss : Signed.long -> Signed.long -> int
val cmpui : pari_ulong -> ('kind, 'structure) t -> int
val cmpuu : pari_ulong -> pari_ulong -> int

val divii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val diviiz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val divirz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val divisz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> unit

val divriz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val divrrz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val divrsz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> unit

val divsi_rem :
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val divsiz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val divsrz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val divss : Signed.long -> Signed.long -> ('kind, 'structure) t

val divss_rem :
  Signed.long ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  ('kind, 'structure) t

val divssz : Signed.long -> Signed.long -> ('kind, 'structure) t -> unit
val dvdii : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val dvdiiz :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t -> int

val dvdis : ('kind, 'structure) t -> Signed.long -> int

val dvdisz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> int

val dvdiu : ('kind, 'structure) t -> pari_ulong -> int
val dvdiuz : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t -> int
val dvdsi : Signed.long -> ('kind, 'structure) t -> int
val dvdui : pari_ulong -> ('kind, 'structure) t -> int

val dvmdiiz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val dvmdis :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val dvmdisz :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val dvmdsbil : Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val dvmdsi :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val dvmdsiz :
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val dvmdss :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val dvmdssz :
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val dvmdubil : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val equalis : ('kind, 'structure) t -> Signed.long -> int
val equalsi : Signed.long -> ('kind, 'structure) t -> int
val equalui : pari_ulong -> ('kind, 'structure) t -> int
val equaliu : ('kind, 'structure) t -> pari_ulong -> int
val absequalui : pari_ulong -> ('kind, 'structure) t -> int
val ceildivuu : pari_ulong -> pari_ulong -> pari_ulong
val evalexpo : Signed.long -> Signed.long
val evallg : Signed.long -> Signed.long
val evalprecp : Signed.long -> Signed.long
val evalvalp : Signed.long -> Signed.long
val evalvalser : Signed.long -> Signed.long
val expi : ('kind, 'structure) t -> Signed.long
val expu : pari_ulong -> Signed.long
val fixlg : ('kind, 'structure) t -> Signed.long -> unit
val fractor : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gc_bool : pari_ulong -> int -> int
val gc_const : pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t
val gc_double : pari_ulong -> float -> float
val gc_int : pari_ulong -> int -> int
val gc_long : pari_ulong -> Signed.long -> Signed.long
val gc_stoi : pari_ulong -> Signed.long -> ('kind, 'structure) t
val gc_ulong : pari_ulong -> pari_ulong -> pari_ulong
val gc_utoi : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val gc_utoipos : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val gc_null : pari_ulong -> ('kind, 'structure) t
val icopy : ('kind, 'structure) t -> ('kind, 'structure) t
val icopyspec : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val int_bit : ('kind, 'structure) t -> Signed.long -> pari_ulong
val itor : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val itos : ('kind, 'structure) t -> Signed.long
val itos_or_0 : ('kind, 'structure) t -> Signed.long
val itou : ('kind, 'structure) t -> pari_ulong
val itou_or_0 : ('kind, 'structure) t -> pari_ulong
val leafcopy : ('kind, 'structure) t -> ('kind, 'structure) t
val maxdd : float -> float -> float
val maxss : Signed.long -> Signed.long -> Signed.long
val maxuu : pari_ulong -> pari_ulong -> Signed.long
val mindd : float -> float -> float
val minss : Signed.long -> Signed.long -> Signed.long
val minuu : pari_ulong -> pari_ulong -> Signed.long
val mod16 : ('kind, 'structure) t -> Signed.long
val mod2 : ('kind, 'structure) t -> Signed.long
val mod2bil : ('kind, 'structure) t -> pari_ulong
val mod32 : ('kind, 'structure) t -> Signed.long
val mod4 : ('kind, 'structure) t -> Signed.long
val mod64 : ('kind, 'structure) t -> Signed.long
val mod8 : ('kind, 'structure) t -> Signed.long
val modis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val modisz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> unit

val modsi : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val modsiz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val modss : Signed.long -> Signed.long -> ('kind, 'structure) t
val modssz : Signed.long -> Signed.long -> ('kind, 'structure) t -> unit
val mpabs : ('kind, 'structure) t -> ('kind, 'structure) t
val mpabs_shallow : ('kind, 'structure) t -> ('kind, 'structure) t

val mpadd :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mpaddz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val mpaff : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val mpceil : ('kind, 'structure) t -> ('kind, 'structure) t
val mpcmp : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val mpcopy : ('kind, 'structure) t -> ('kind, 'structure) t

val mpdiv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mpexpo : ('kind, 'structure) t -> Signed.long
val mpfloor : ('kind, 'structure) t -> ('kind, 'structure) t

val mpmul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mpmulz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val mpneg : ('kind, 'structure) t -> ('kind, 'structure) t
val mpodd : ('kind, 'structure) t -> int
val mpround : ('kind, 'structure) t -> ('kind, 'structure) t
val mpshift : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mpsqr : ('kind, 'structure) t -> ('kind, 'structure) t

val mpsub :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mpsubz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val mptrunc : ('kind, 'structure) t -> ('kind, 'structure) t

val muliiz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val mulirz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val mulis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val muliu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val mulri :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mulriz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val mulrrz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val mulrs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val mulru : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val mulsiz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val mulsrz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val mulssz : Signed.long -> Signed.long -> ('kind, 'structure) t -> unit
val negi : ('kind, 'structure) t -> ('kind, 'structure) t
val negr : ('kind, 'structure) t -> ('kind, 'structure) t
val new_chunk : int -> ('kind, 'structure) t
val rcopy : ('kind, 'structure) t -> ('kind, 'structure) t

val rdivii :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val rdiviiz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val rdivis :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val rdivsi :
  Signed.long -> ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rdivss : Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t
val real2n : Signed.long -> Signed.long -> ('kind, 'structure) t
val real_m2n : Signed.long -> Signed.long -> ('kind, 'structure) t
val real_0 : Signed.long -> ('kind, 'structure) t
val real_0_bit : Signed.long -> ('kind, 'structure) t
val real_1 : Signed.long -> ('kind, 'structure) t
val real_1_bit : Signed.long -> ('kind, 'structure) t
val real_m1 : Signed.long -> ('kind, 'structure) t

val remii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val remiiz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val remis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val remisz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> unit

val remlll_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val remsi : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val remsiz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val remss : Signed.long -> Signed.long -> ('kind, 'structure) t
val remssz : Signed.long -> Signed.long -> ('kind, 'structure) t -> unit
val rtor : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val sdivsi : Signed.long -> ('kind, 'structure) t -> Signed.long

val sdivsi_rem :
  Signed.long ->
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val sdivss_rem :
  Signed.long -> Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val get_avma : unit -> pari_ulong
val set_avma : pari_ulong -> unit

val uabsdiviu_rem :
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong

val udivuu_rem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val umodi2n : ('kind, 'structure) t -> Signed.long -> pari_ulong
val setabssign : ('kind, 'structure) t -> unit

val shift_left :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  unit

val shift_right :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  pari_ulong ->
  pari_ulong ->
  unit

val shiftl : pari_ulong -> pari_ulong -> pari_ulong
val shiftlr : pari_ulong -> pari_ulong -> pari_ulong
val shiftr : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val shiftr_inplace : ('kind, 'structure) t -> Signed.long -> unit
val smodis : ('kind, 'structure) t -> Signed.long -> Signed.long
val smodss : Signed.long -> Signed.long -> Signed.long
val stackdummy : pari_ulong -> pari_ulong -> unit
val stack_malloc : int -> string
val stack_malloc_align : int -> Signed.long -> string
val stack_calloc : int -> string
val stack_calloc_align : int -> Signed.long -> string

val subii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val subiiz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val subir :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val subirz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val subis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val subisz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> unit

val subri :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val subriz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val subrr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val subrrz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val subrs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val subrsz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> unit

val subsi : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val subsiz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val subsrz :
  Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val subss : Signed.long -> Signed.long -> ('kind, 'structure) t
val subssz : Signed.long -> Signed.long -> ('kind, 'structure) t -> unit
val subuu : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val togglesign : ('kind, 'structure) t -> unit
val togglesign_safe : ('kind, 'structure) t Ctypes_static.ptr -> unit
val affectsign : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val affectsign_safe :
  ('kind, 'structure) t -> ('kind, 'structure) t Ctypes_static.ptr -> unit

val truedivii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val truedivis : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val truedivsi : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val uabsdivui_rem :
  pari_ulong ->
  ('kind, 'structure) t ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong

val umodsu : Signed.long -> pari_ulong -> pari_ulong
val umodui : pari_ulong -> ('kind, 'structure) t -> pari_ulong
val ugcdiu : ('kind, 'structure) t -> pari_ulong -> pari_ulong
val ugcdui : pari_ulong -> ('kind, 'structure) t -> pari_ulong
val umuluu_le : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val umuluu_or_0 : pari_ulong -> pari_ulong -> pari_ulong
val utoi : pari_ulong -> ('kind, 'structure) t
val utoineg : pari_ulong -> ('kind, 'structure) t
val utoipos : pari_ulong -> ('kind, 'structure) t
val utor : pari_ulong -> Signed.long -> ('kind, 'structure) t
val uutoi : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val uutoineg : pari_ulong -> pari_ulong -> ('kind, 'structure) t
val vali : ('kind, 'structure) t -> Signed.long
val varncmp : Signed.long -> Signed.long -> int
val varnmax : Signed.long -> Signed.long -> Signed.long
val varnmin : Signed.long -> Signed.long -> Signed.long

val pari_err_component :
  string -> string -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val pari_err_dim : string -> unit

val pari_err_domain :
  string ->
  string ->
  string ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val pari_err_file : string -> string -> unit
val pari_err_filedesc : string -> Signed.long -> unit
val pari_err_flag : string -> unit
val pari_err_impl : string -> unit
val pari_err_inv : string -> ('kind, 'structure) t -> unit
val pari_err_irredpol : string -> ('kind, 'structure) t -> unit
val pari_err_maxprime : pari_ulong -> unit

val pari_err_modulus :
  string -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val pari_err_op :
  string -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val pari_err_overflow : string -> unit
val pari_err_package : string -> unit
val pari_err_prec : string -> unit
val pari_err_prime : string -> ('kind, 'structure) t -> unit

val pari_err_priority :
  string -> ('kind, 'structure) t -> string -> Signed.long -> unit

val pari_err_sqrtn : string -> ('kind, 'structure) t -> unit
val pari_err_type : string -> ('kind, 'structure) t -> unit

val pari_err_type2 :
  string -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val pari_err_var :
  string -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val pari_err_roots0 : string -> unit

val mkintmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkintmodu : pari_ulong -> pari_ulong -> ('kind, 'structure) t

val mkpolmod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkfrac :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkfracss : Signed.long -> Signed.long -> ('kind, 'structure) t

val qtoss :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val sstoq : Signed.long -> Signed.long -> ('kind, 'structure) t
val uutoq : pari_ulong -> pari_ulong -> ('kind, 'structure) t

val mkfraccopy :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkrfrac :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkrfraccopy :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkcomplex :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gen_i : unit -> ('kind, 'structure) t
val cgetc : Signed.long -> ('kind, 'structure) t

val mkquad :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkvecsmall : Signed.long -> ('kind, 'structure) t
val mkvecsmall2 : Signed.long -> Signed.long -> ('kind, 'structure) t

val mkvecsmall3 :
  Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mkvecsmall4 :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val mkvecsmall5 :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val mkqfb :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkvec : ('kind, 'structure) t -> ('kind, 'structure) t

val mkvec2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkvec3 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkvec4 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkvec5 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkvecs : Signed.long -> ('kind, 'structure) t
val mkvec2s : Signed.long -> Signed.long -> ('kind, 'structure) t
val mkvec3s : Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mkvec4s :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val mkveccopy : ('kind, 'structure) t -> ('kind, 'structure) t

val mkvec2copy :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkcol : ('kind, 'structure) t -> ('kind, 'structure) t

val mkcol2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkcol3 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkcol4 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkcol5 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkcol6 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkcols : Signed.long -> ('kind, 'structure) t
val mkcol2s : Signed.long -> Signed.long -> ('kind, 'structure) t
val mkcol3s : Signed.long -> Signed.long -> Signed.long -> ('kind, 'structure) t

val mkcol4s :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val mkcolcopy : ('kind, 'structure) t -> ('kind, 'structure) t
val mkmat : ('kind, 'structure) t -> ('kind, 'structure) t

val mkmat2 :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mkmat3 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkmat4 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkmat5 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mkmatcopy : ('kind, 'structure) t -> ('kind, 'structure) t
val mkerr : Signed.long -> ('kind, 'structure) t
val mkoo : unit -> ('kind, 'structure) t
val mkmoo : unit -> ('kind, 'structure) t
val inf_get_sign : ('kind, 'structure) t -> Signed.long

val mkmat22s :
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val mkmat22 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val pol_x : Signed.long -> ('kind, 'structure) t
val pol_xn : Signed.long -> Signed.long -> ('kind, 'structure) t
val pol_xnall : Signed.long -> Signed.long -> ('kind, 'structure) t
val polxn_flx : Signed.long -> Signed.long -> ('kind, 'structure) t
val pol_1 : Signed.long -> ('kind, 'structure) t
val pol_0 : Signed.long -> ('kind, 'structure) t
val const_vec : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val const_col : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val const_vecsmall : Signed.long -> Signed.long -> ('kind, 'structure) t
val zeropadic : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zeropadic_shallow :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zeroser : Signed.long -> Signed.long -> ('kind, 'structure) t
val ser_isexactzero : ('kind, 'structure) t -> int
val zeropol : Signed.long -> ('kind, 'structure) t
val zerocol : Signed.long -> ('kind, 'structure) t
val zerovec : Signed.long -> ('kind, 'structure) t
val zeromat : Signed.long -> Signed.long -> ('kind, 'structure) t
val zero_flx : Signed.long -> ('kind, 'structure) t
val zero_flv : Signed.long -> ('kind, 'structure) t
val zero_flm : Signed.long -> Signed.long -> ('kind, 'structure) t
val zero_flm_copy : Signed.long -> Signed.long -> ('kind, 'structure) t
val zero_f2v : Signed.long -> ('kind, 'structure) t
val zero_f2m : Signed.long -> Signed.long -> ('kind, 'structure) t
val zero_f2m_copy : Signed.long -> Signed.long -> ('kind, 'structure) t
val zeromatcopy : Signed.long -> Signed.long -> ('kind, 'structure) t
val zerovec_block : Signed.long -> ('kind, 'structure) t
val col_ei : Signed.long -> Signed.long -> ('kind, 'structure) t
val vec_ei : Signed.long -> Signed.long -> ('kind, 'structure) t
val f2v_ei : Signed.long -> Signed.long -> ('kind, 'structure) t
val vecsmall_ei : Signed.long -> Signed.long -> ('kind, 'structure) t

val rg_col_ei :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val shallowcopy : ('kind, 'structure) t -> ('kind, 'structure) t
val vectrunc_init : Signed.long -> ('kind, 'structure) t
val coltrunc_init : Signed.long -> ('kind, 'structure) t
val lg_increase : ('kind, 'structure) t -> unit
val vectrunc_append : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val vectrunc_append_batch :
  ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val vecsmalltrunc_init : Signed.long -> ('kind, 'structure) t
val vecsmalltrunc_append : ('kind, 'structure) t -> Signed.long -> unit
val hash_str : string -> pari_ulong
val hash_str_len : string -> Signed.long -> pari_ulong
val vec_shorten : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val vec_lengthen : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val vec_append :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vec_prepend :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vec_setconst :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecsmall_shorten :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val vecsmall_lengthen :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val vec_to_vecsmall : ('kind, 'structure) t -> ('kind, 'structure) t
val vecsmall_to_vec : ('kind, 'structure) t -> ('kind, 'structure) t
val vecsmall_to_vec_inplace : ('kind, 'structure) t -> ('kind, 'structure) t
val vecsmall_to_col : ('kind, 'structure) t -> ('kind, 'structure) t
val vecsmall_lexcmp : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val vecsmall_prefixcmp : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val vecsmall_prepend :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val vecsmall_append :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val vecsmall_concat :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecsmall_coincidence :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val vecsmall_isin : ('kind, 'structure) t -> Signed.long -> Signed.long

val vecsmall_pack :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> Signed.long

val vecsmall_indexmax : ('kind, 'structure) t -> Signed.long
val vecsmall_max : ('kind, 'structure) t -> Signed.long
val vecsmall_indexmin : ('kind, 'structure) t -> Signed.long
val vecsmall_min : ('kind, 'structure) t -> Signed.long
val zv_isscalar : ('kind, 'structure) t -> int
val qv_isscalar : ('kind, 'structure) t -> int
val rgv_isscalar : ('kind, 'structure) t -> int
val rgx_isscalar : ('kind, 'structure) t -> int

val rgx_equal_var :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val rgx_to_rgv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgx_is_rational : ('kind, 'structure) t -> int
val rgx_is_zx : ('kind, 'structure) t -> int
val rgx_is_qx : ('kind, 'structure) t -> int
val rgx_is_monomial : ('kind, 'structure) t -> int
val rgv_is_zv : ('kind, 'structure) t -> int
val rgv_is_qv : ('kind, 'structure) t -> int

val rgv_isin_i :
  ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long -> Signed.long

val rgv_isin : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long

val vecslice :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val vecslicepermute :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rowslicepermute :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rowslice :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val matslice :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val rowsplice : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val vecsplice : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val rgm_minor :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val row : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val flm_row : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rowcopy : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val row_i :
  ('kind, 'structure) t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  ('kind, 'structure) t

val vecreverse : ('kind, 'structure) t -> ('kind, 'structure) t
val vecsmall_reverse : ('kind, 'structure) t -> ('kind, 'structure) t
val vecreverse_inplace : ('kind, 'structure) t -> unit

val vecsmallpermute :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecpermute :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rowpermute :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val identity_zv : Signed.long -> ('kind, 'structure) t
val identity_perm : Signed.long -> ('kind, 'structure) t
val cyclic_perm : Signed.long -> Signed.long -> ('kind, 'structure) t

val perm_mul :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val perm_sqr : ('kind, 'structure) t -> ('kind, 'structure) t
val perm_inv : ('kind, 'structure) t -> ('kind, 'structure) t

val perm_conj :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val pari_free : unit Ctypes_static.ptr -> unit
val pari_malloc : int -> unit Ctypes_static.ptr
val pari_realloc : unit Ctypes_static.ptr -> int -> unit Ctypes_static.ptr
val pari_realloc_ip : unit Ctypes_static.ptr Ctypes_static.ptr -> int -> unit
val pari_calloc : int -> unit Ctypes_static.ptr
val cgetalloc : int -> Signed.long -> ('kind, 'structure) t
val icopy_avma : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val leafcopy_avma : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t

val gerepileuptoleaf :
  pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t

val gerepileuptoint :
  pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t

val gerepileupto : pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t
val gerepilecopy : pari_ulong -> ('kind, 'structure) t -> ('kind, 'structure) t
val gunclonenull : ('kind, 'structure) t -> unit
val gunclonenull_deep : ('kind, 'structure) t -> unit

val gerepilemany :
  pari_ulong ->
  ('kind, 'structure) t Ctypes_static.ptr Ctypes_static.ptr ->
  int ->
  unit

val gerepileall : pari_ulong -> int -> unit
val gc_all : pari_ulong -> int -> ('kind, 'structure) t
val gerepilecoeffs : pari_ulong -> ('kind, 'structure) t -> int -> unit

val bin_copy :
  genbin Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val genbinbase :
  genbin Ctypes.structure Ctypes_static.ptr -> ('kind, 'structure) t

val cgiv : ('kind, 'structure) t -> unit
val killblock : ('kind, 'structure) t -> unit
val is_universal_constant : ('kind, 'structure) t -> int
val cxcompotor : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val cxtofp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val cxtoreal : ('kind, 'structure) t -> ('kind, 'structure) t
val gtodouble : ('kind, 'structure) t -> float
val gisdouble : ('kind, 'structure) t -> float Ctypes_static.ptr -> int
val gtos : ('kind, 'structure) t -> Signed.long
val gtou : ('kind, 'structure) t -> pari_ulong
val absfrac : ('kind, 'structure) t -> ('kind, 'structure) t
val absfrac_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val q_abs : ('kind, 'structure) t -> ('kind, 'structure) t
val q_abs_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val r_abs_shallow : ('kind, 'structure) t -> ('kind, 'structure) t
val r_abs : ('kind, 'structure) t -> ('kind, 'structure) t
val gtofp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gtomp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgx_gtofp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgc_gtofp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgv_gtofp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgm_gtofp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgc_gtomp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgm_gtomp : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgx_fpnorml2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgc_fpnorml2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgm_fpnorml2 : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val affgr : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val affc_fixlg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val trunc_safe : ('kind, 'structure) t -> ('kind, 'structure) t
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
val bit_prec : ('kind, 'structure) t -> Signed.long
val bit_accuracy : Signed.long -> Signed.long
val prec2ndec : Signed.long -> Signed.long
val nbits2ndec : Signed.long -> Signed.long
val precdbl : Signed.long -> Signed.long
val divsbil : Signed.long -> Signed.long
val remsbil : Signed.long -> Signed.long

val fp_red :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fp_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_neg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fp_halve :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fp_center :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_center_i :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_addmul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fp_mulu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_muls :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_inv :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fp_invsafe :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fp_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fp_divu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val flx_mulu :
  ('kind, 'structure) t -> pari_ulong -> pari_ulong -> ('kind, 'structure) t

val get_f2x_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val get_f2x_var : ('kind, 'structure) t -> Signed.long
val get_f2x_degree : ('kind, 'structure) t -> Signed.long
val get_f2xqx_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val get_f2xqx_var : ('kind, 'structure) t -> Signed.long
val get_f2xqx_degree : ('kind, 'structure) t -> Signed.long
val get_flx_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val get_flx_var : ('kind, 'structure) t -> Signed.long
val get_flx_degree : ('kind, 'structure) t -> Signed.long
val get_flxqx_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val get_flxqx_var : ('kind, 'structure) t -> Signed.long
val get_flxqx_degree : ('kind, 'structure) t -> Signed.long
val get_fpx_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val get_fpx_var : ('kind, 'structure) t -> Signed.long
val get_fpx_degree : ('kind, 'structure) t -> Signed.long
val get_fpxqx_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val get_fpxqx_var : ('kind, 'structure) t -> Signed.long
val get_fpxqx_degree : ('kind, 'structure) t -> Signed.long

val submulii :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val mulsubii :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val submuliu :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val addmuliu :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val submuliu_inplace :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val addmuliu_inplace :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val lincombii :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

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
val qfb_is_qfi : ('kind, 'structure) t -> int
val sqrtr : ('kind, 'structure) t -> ('kind, 'structure) t
val cbrtr_abs : ('kind, 'structure) t -> ('kind, 'structure) t
val cbrtr : ('kind, 'structure) t -> ('kind, 'structure) t
val sqrtnr : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val logint : ('kind, 'structure) t -> ('kind, 'structure) t -> Signed.long
val ulogint : pari_ulong -> pari_ulong -> pari_ulong
val ismpzero : ('kind, 'structure) t -> int
val isintzero : ('kind, 'structure) t -> int
val isint1 : ('kind, 'structure) t -> int
val isintm1 : ('kind, 'structure) t -> int
val equali1 : ('kind, 'structure) t -> int
val equalim1 : ('kind, 'structure) t -> int
val is_pm1 : ('kind, 'structure) t -> int
val is_bigint : ('kind, 'structure) t -> int
val odd : Signed.long -> int
val both_odd : Signed.long -> Signed.long -> int
val isonstack : ('kind, 'structure) t -> int
val dbllog2r : ('kind, 'structure) t -> float

val mul_content :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val inv_content : ('kind, 'structure) t -> ('kind, 'structure) t

val div_content :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val mul_denom :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val constant_coeff : ('kind, 'structure) t -> ('kind, 'structure) t
val leading_coeff : ('kind, 'structure) t -> ('kind, 'structure) t
val flx_lead : ('kind, 'structure) t -> pari_ulong
val flx_constant : ('kind, 'structure) t -> pari_ulong
val degpol : ('kind, 'structure) t -> Signed.long
val lgpol : ('kind, 'structure) t -> Signed.long
val lgcols : ('kind, 'structure) t -> Signed.long
val nbrows : ('kind, 'structure) t -> Signed.long
val truecoef : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val zxq_sqr :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val rgx_coeff : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val rgx_renormalize : ('kind, 'structure) t -> ('kind, 'structure) t

val rgx_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val rgxqx_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rgxqx_rem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpx_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val flx_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flx_div_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val f2x_div :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpv_fpc_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val pol0_flx : Signed.long -> ('kind, 'structure) t
val pol1_flx : Signed.long -> ('kind, 'structure) t
val polx_flx : Signed.long -> ('kind, 'structure) t
val zero_zx : Signed.long -> ('kind, 'structure) t
val polx_zx : Signed.long -> ('kind, 'structure) t
val zx_shift : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val zero_f2x : Signed.long -> ('kind, 'structure) t
val pol0_f2x : Signed.long -> ('kind, 'structure) t
val pol1_f2x : Signed.long -> ('kind, 'structure) t
val polx_f2x : Signed.long -> ('kind, 'structure) t
val f2x_equal1 : ('kind, 'structure) t -> int
val f2x_equal : ('kind, 'structure) t -> ('kind, 'structure) t -> int
val f2x_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val f2v_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val flv_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val flx_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val vecsmall_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val flx_equal1 : ('kind, 'structure) t -> int
val zx_equal1 : ('kind, 'structure) t -> int
val zx_is_monic : ('kind, 'structure) t -> int

val zx_renormalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val fpx_renormalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val fpxx_renormalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val fpxqx_renormalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val f2x_renormalize :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val f2xx_shift :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> ('kind, 'structure) t

val f2v_to_f2x : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val sturm : ('kind, 'structure) t -> Signed.long
val gval : ('kind, 'structure) t -> Signed.long -> Signed.long
val rgx_shift_inplace_init : Signed.long -> unit

val rgx_shift_inplace :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zc_to_zc : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_to_zv : ('kind, 'structure) t -> ('kind, 'structure) t
val zx_to_zv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val zv_to_zx : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val zm_to_zxv : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val zero_zm : Signed.long -> Signed.long -> ('kind, 'structure) t
val zero_zv : Signed.long -> ('kind, 'structure) t
val zm_transpose : ('kind, 'structure) t -> ('kind, 'structure) t
val zm_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val zm_row : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val zc_hnfrem :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_hnfrem :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zm_lll :
  ('kind, 'structure) t -> float -> Signed.long -> ('kind, 'structure) t

val rgm_dimensions :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val rgm_shallowcopy : ('kind, 'structure) t -> ('kind, 'structure) t
val f2m_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val f3m_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val flm_copy : ('kind, 'structure) t -> ('kind, 'structure) t
val zv_dvd : ('kind, 'structure) t -> ('kind, 'structure) t -> int

val zm_zv_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val zv_zv_mod :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecmodii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val vecmoduu :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fq_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fq_to_fpxq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val rg_to_fq :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val gener_fq_local :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val random_fq :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val fpxqx_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val flxqx_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxqx_div_pre :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  pari_ulong ->
  ('kind, 'structure) t

val f2xqx_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxy_fq_evaly :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t

val fqx_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_neg :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_fp_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_fq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_mulu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_powu :
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_halve :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_get_red :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_rem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_divrem :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fqx_div_by_x_x :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fqx_halfgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_gcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_extgcd :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t Ctypes_static.ptr ->
  ('kind, 'structure) t

val fqx_normalize :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_deriv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_integ :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_factor :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_factor_squarefree :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_ddf :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_degfact :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_roots :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqx_to_mod :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_div :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_inv :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_invsafe :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_sqr :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxq_pow :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxn_expint :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxn_exp :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxn_inv :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxn_mul :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fqxn_sqr :
  ('kind, 'structure) t ->
  Signed.long ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val fpxq_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val flxq_add :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val flxq_sub :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  pari_ulong ->
  ('kind, 'structure) t

val f2x_coeff : ('kind, 'structure) t -> Signed.long -> pari_ulong
val f2x_clear : ('kind, 'structure) t -> Signed.long -> unit
val f2x_set : ('kind, 'structure) t -> Signed.long -> unit
val f2x_flip : ('kind, 'structure) t -> Signed.long -> unit
val f2v_coeff : ('kind, 'structure) t -> Signed.long -> pari_ulong
val f2v_clear : ('kind, 'structure) t -> Signed.long -> unit
val f2v_set : ('kind, 'structure) t -> Signed.long -> unit
val f2v_flip : ('kind, 'structure) t -> Signed.long -> unit

val f2m_coeff :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> pari_ulong

val f2m_clear : ('kind, 'structure) t -> Signed.long -> Signed.long -> unit
val f2m_set : ('kind, 'structure) t -> Signed.long -> Signed.long -> unit
val f2m_flip : ('kind, 'structure) t -> Signed.long -> Signed.long -> unit

val f3m_coeff :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> pari_ulong

val f3m_set :
  ('kind, 'structure) t -> Signed.long -> Signed.long -> pari_ulong -> unit

val matpascal : Signed.long -> ('kind, 'structure) t
val z_issquare : ('kind, 'structure) t -> Signed.long
val z_ispower : ('kind, 'structure) t -> pari_ulong -> Signed.long
val sqrti : ('kind, 'structure) t -> ('kind, 'structure) t
val gaddgs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gcmpgs : ('kind, 'structure) t -> Signed.long -> int
val gequalgs : ('kind, 'structure) t -> Signed.long -> int
val gmaxsg : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val gminsg : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t
val gmulgs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gmulgu : ('kind, 'structure) t -> pari_ulong -> ('kind, 'structure) t
val gsubgs : ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t
val gdivsg : Signed.long -> ('kind, 'structure) t -> ('kind, 'structure) t

val gmax_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val gmin_shallow :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val cxnorm : ('kind, 'structure) t -> ('kind, 'structure) t
val quadnorm : ('kind, 'structure) t -> ('kind, 'structure) t
val quad_disc : ('kind, 'structure) t -> ('kind, 'structure) t

val qfb_disc3 :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t

val qfb_disc : ('kind, 'structure) t -> ('kind, 'structure) t
val sqrfrac : ('kind, 'structure) t -> ('kind, 'structure) t
val normalize_frac : ('kind, 'structure) t -> unit

val powii :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val powis : Signed.long -> ('kind, 'structure) t
val mpexpz : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val mplogz : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val mpcosz : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val mpsinz : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val gnegz : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val gabsz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> unit

val gaddz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val gsubz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val gmulz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val gdivz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val gdiventz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val gmodz :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val gmul2nz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> unit

val gshiftz :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t -> unit

val ell_get_a1 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_a2 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_a3 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_a4 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_a6 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_b2 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_b4 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_b6 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_b8 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_c4 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_c6 : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_disc : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_j : ('kind, 'structure) t -> ('kind, 'structure) t
val ell_get_type : ('kind, 'structure) t -> Signed.long
val ellff_get_field : ('kind, 'structure) t -> ('kind, 'structure) t
val ellff_get_a4a6 : ('kind, 'structure) t -> ('kind, 'structure) t
val ellqp_get_zero : ('kind, 'structure) t -> ('kind, 'structure) t
val ellqp_get_prec : ('kind, 'structure) t -> Signed.long
val ellqp_get_p : ('kind, 'structure) t -> ('kind, 'structure) t
val ellr_get_prec : ('kind, 'structure) t -> Signed.long
val ellr_get_sign : ('kind, 'structure) t -> Signed.long
val ellnf_get_nf : ('kind, 'structure) t -> ('kind, 'structure) t
val ellnf_get_bnf : ('kind, 'structure) t -> ('kind, 'structure) t
val checkell_i : ('kind, 'structure) t -> int
val ell_is_inf : ('kind, 'structure) t -> int
val ellinf : unit -> ('kind, 'structure) t
val modpr_get_pr : ('kind, 'structure) t -> ('kind, 'structure) t
val modpr_get_p : ('kind, 'structure) t -> ('kind, 'structure) t
val modpr_get_t : ('kind, 'structure) t -> ('kind, 'structure) t
val pr_get_p : ('kind, 'structure) t -> ('kind, 'structure) t
val pr_get_gen : ('kind, 'structure) t -> ('kind, 'structure) t
val pr_get_e : ('kind, 'structure) t -> Signed.long
val pr_get_f : ('kind, 'structure) t -> Signed.long
val pr_get_tau : ('kind, 'structure) t -> ('kind, 'structure) t
val pr_is_inert : ('kind, 'structure) t -> int
val pr_norm : ('kind, 'structure) t -> ('kind, 'structure) t
val upr_norm : ('kind, 'structure) t -> pari_ulong
val nf_get_varn : ('kind, 'structure) t -> Signed.long
val nf_get_pol : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_degree : ('kind, 'structure) t -> Signed.long
val nf_get_r1 : ('kind, 'structure) t -> Signed.long
val nf_get_r2 : ('kind, 'structure) t -> Signed.long
val nf_get_disc : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_index : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_m : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_g : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_roundg : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_tr : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_diff : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_ramified_primes : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_roots : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_zk : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_zkprimpart : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_zkden : ('kind, 'structure) t -> ('kind, 'structure) t
val nf_get_invzk : ('kind, 'structure) t -> ('kind, 'structure) t

val nf_get_sign :
  ('kind, 'structure) t ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val cyc_get_expo : ('kind, 'structure) t -> ('kind, 'structure) t
val abgrp_get_no : ('kind, 'structure) t -> ('kind, 'structure) t
val abgrp_get_cyc : ('kind, 'structure) t -> ('kind, 'structure) t
val abgrp_get_gen : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_nf : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_clgp : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_no : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_cyc : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_gen : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_reg : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_logfu : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_sunits : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_tuu : ('kind, 'structure) t -> ('kind, 'structure) t
val bnf_get_tun : ('kind, 'structure) t -> Signed.long
val bnf_get_fu_nocheck : ('kind, 'structure) t -> ('kind, 'structure) t

val nfv_to_scalar_or_alg :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val bnf_get_fu : ('kind, 'structure) t -> ('kind, 'structure) t
val bnr_get_bnf : ('kind, 'structure) t -> ('kind, 'structure) t
val bnr_get_bid : ('kind, 'structure) t -> ('kind, 'structure) t
val bnr_get_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val bnr_get_nf : ('kind, 'structure) t -> ('kind, 'structure) t
val bnr_get_clgp : ('kind, 'structure) t -> ('kind, 'structure) t
val bnr_get_no : ('kind, 'structure) t -> ('kind, 'structure) t
val bnr_get_cyc : ('kind, 'structure) t -> ('kind, 'structure) t
val bnr_get_gen_nocheck : ('kind, 'structure) t -> ('kind, 'structure) t
val bnr_get_gen : ('kind, 'structure) t -> ('kind, 'structure) t
val locs_get_cyc : ('kind, 'structure) t -> ('kind, 'structure) t
val locs_get_lsprk : ('kind, 'structure) t -> ('kind, 'structure) t
val locs_get_lgenfil : ('kind, 'structure) t -> ('kind, 'structure) t
val locs_get_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val locs_get_famod : ('kind, 'structure) t -> ('kind, 'structure) t
val locs_get_m_infty : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_basis : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_bnf : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_nf : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_zm : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_modp : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_s : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_dldata : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_sfu : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_cyc : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_hnf : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_u : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_ui : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_m0 : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_u0 : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_r1 : ('kind, 'structure) t -> Signed.long
val gchar_get_r2 : ('kind, 'structure) t -> Signed.long
val gchar_get_loccyc : ('kind, 'structure) t -> ('kind, 'structure) t
val gchar_get_nc : ('kind, 'structure) t -> Signed.long
val gchar_get_ns : ('kind, 'structure) t -> Signed.long
val gchar_get_nm : ('kind, 'structure) t -> Signed.long
val gchar_get_evalprec : ('kind, 'structure) t -> Signed.long
val gchar_get_prec : ('kind, 'structure) t -> Signed.long
val gchar_get_nfprec : ('kind, 'structure) t -> Signed.long
val gchar_set_evalprec : ('kind, 'structure) t -> Signed.long -> unit
val gchar_set_prec : ('kind, 'structure) t -> Signed.long -> unit
val gchar_copy_precs : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val gchar_set_nfprec : ('kind, 'structure) t -> Signed.long -> unit
val gchar_get_ntors : ('kind, 'structure) t -> Signed.long
val gchar_get_nfree : ('kind, 'structure) t -> Signed.long
val gchar_get_nalg : ('kind, 'structure) t -> Signed.long
val gchar_set_basis : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val gchar_set_nf : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val gchar_set_ntors : ('kind, 'structure) t -> Signed.long -> unit
val gchar_set_nfree : ('kind, 'structure) t -> Signed.long -> unit
val gchar_set_nalg : ('kind, 'structure) t -> Signed.long -> unit
val gchar_set_cyc : ('kind, 'structure) t -> ('kind, 'structure) t -> unit

val gchar_set_huui :
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  ('kind, 'structure) t ->
  unit

val gchar_set_m0 : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val gchar_set_u0 : ('kind, 'structure) t -> ('kind, 'structure) t -> unit
val bid_get_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_ideal : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_arch : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_grp : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_fact : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_fact2 : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_sprk : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_sarch : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_archp : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_u : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_no : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_cyc : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_gen_nocheck : ('kind, 'structure) t -> ('kind, 'structure) t
val bid_get_gen : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_n : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_fan : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_no : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_cyc : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_gen : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_conreycyc : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_conreygen : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_ui : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_u : ('kind, 'structure) t -> ('kind, 'structure) t
val znstar_get_pe : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_pol : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_p : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_e : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_mod : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_roots : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_invvdm : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_den : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_group : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_gen : ('kind, 'structure) t -> ('kind, 'structure) t
val gal_get_orders : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_degree : ('kind, 'structure) t -> Signed.long
val rnf_get_nfdegree : ('kind, 'structure) t -> Signed.long
val rnf_get_absdegree : ('kind, 'structure) t -> Signed.long
val rnf_get_idealdisc : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_k : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_alpha : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_nf : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_nfzk : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_polabs : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_pol : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_disc : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_index : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_ramified_primes : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_varn : ('kind, 'structure) t -> Signed.long
val rnf_get_nfpol : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_nfvarn : ('kind, 'structure) t -> Signed.long
val rnf_get_zk : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_map : ('kind, 'structure) t -> ('kind, 'structure) t
val rnf_get_invzk : ('kind, 'structure) t -> ('kind, 'structure) t

val idealred :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val idealchineseinit :
  ('kind, 'structure) t -> ('kind, 'structure) t -> ('kind, 'structure) t

val closure_arity : ('kind, 'structure) t -> Signed.long
val closure_is_variadic : ('kind, 'structure) t -> Signed.long
val closure_codestr : ('kind, 'structure) t -> string
val closure_get_code : ('kind, 'structure) t -> ('kind, 'structure) t
val closure_get_oper : ('kind, 'structure) t -> ('kind, 'structure) t
val closure_get_data : ('kind, 'structure) t -> ('kind, 'structure) t
val closure_get_dbg : ('kind, 'structure) t -> ('kind, 'structure) t
val closure_get_text : ('kind, 'structure) t -> ('kind, 'structure) t
val closure_get_frame : ('kind, 'structure) t -> ('kind, 'structure) t
val err_get_num : ('kind, 'structure) t -> Signed.long

val err_get_compo :
  ('kind, 'structure) t -> Signed.long -> ('kind, 'structure) t

val pari_err_bug : string -> unit
val pari_err_constpol : string -> unit

val pari_err_coprime :
  string -> ('kind, 'structure) t -> ('kind, 'structure) t -> unit
