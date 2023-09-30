type pari_ulong = Unsigned.ULong.t

val pari_ulong : pari_ulong Ctypes.typ

type 'a t

val t : 'a t Ctypes.typ

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
type 'a ring
type 'a field

module rec Complex : sig
  type complex = private Complex
  type nonrec t = complex field t

  val inv : t -> t
  val add : t -> t -> t
  val create : re:Real.t -> im:Real.t -> t
  val to_string : t -> string
end

and Real : sig
  type real = private Real
  type nonrec t = real field t

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
  type 'a p := 'a t
  type nonrec t = rational field t
  type nonrec ring = rational ring p

  external inj_ring : t -> ring = "%identity"
  external inj_real : t -> Real.t = "%identity"
  external inj_complex : t -> Complex.t = "%identity"
  val shift : t -> int -> t
end

and Integer : sig
  type integer = private Integer
  type nonrec t = integer ring t

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

module Set : sig
  type nonrec 'a t constraint 'a = 'b t

  val length : 'a t -> Signed.Long.t
  val search : 'a t -> 'a -> Signed.Long.t -> Signed.Long.t
end

module Vector : sig
  type nonrec ('a, 'b) t constraint 'a = 'c t constraint 'b = [< `COL | `ROW ]

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
  type matrix = private Matrix
  type nonrec 'a t = matrix t constraint 'a = 'b t

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
  type polynomial
  type 'a p := 'a t
  type 'a t = polynomial ring p constraint 'a = 'b ring p

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
  val create : ('a * int) list -> 'a t
  val deriv : ?indeterminate:int -> 'a t -> 'a t
  val derivn : ?indeterminate:int -> 'a t -> int -> 'a t
  val cyclotomic : Signed.long -> Integer.t t
  val is_irreducible : 'a t -> bool

  val minimal : 'a t -> 'a t
  (** [minimal p] reduces [p] to be the minimal polynomial of the
      roots of [p] over the field of the rational numbers.
   
      {@ocaml[
      # let q = Polynomial.create
        [
          (Integer.of_int 1, 3);
          (Integer.of_int (-111), 2);
          (Integer.of_int 6064, 1);
          (Integer.of_int (-189804), 0);
        ];;
      val q : Integer.t Polynomial.t = <abstr>
      # Polynomial.to_string q;;
      - : string = "x^3 - 111*x^2 + 6064*x - 189804"
      # let qmin = Polynomial.minimal q;;
      val qmin : 'a ring t Polynomial.t = <abstr>
      # Polynomial.to_string qmin;;
      - : string = "x^3 - x^2 - 60*x - 364"
      # Number_field.(are_isomorphic (create q) (create qmin));
      - : bool = true
      ]} *)

  val ( .%[] ) : 'a t -> int -> 'a

  val roots_ff :
    _ Finite_field.finite_field ring p t ->
    (_ Finite_field.finite_field field p, [ `ROW ]) Vector.t
end

and Fp : sig
  type t = Integer.t

  val add : t -> t -> modulo:t -> t
  val pow : t -> exponent:t -> modulo:t -> t
end

and Finite_field : sig
  type 'a finite_field
  type 'a p := 'a t
  type 'a t = 'a finite_field field p
  type prime
  type prime_field = prime finite_field field p
  type nonrec 'a ring = 'a finite_field ring p

  external inj_ring : _ t -> _ ring = "%identity"
  external inj_field : _ ring -> _ t = "%identity"
  val generator : order:Integer.t -> _ t
  val prime_field_element : Integer.t -> p:Integer.t -> prime t
  val finite_field_element : Integer.t array -> 'a t -> 'a t

  val create : p:int -> degree:int -> _ ring Polynomial.t
  (** [create p degree] returns a monic irreducible polynomial of the
      given [degree] over F_p[X]. *)

  val generator_from_irreducible_polynomial : _ ring Polynomial.t -> _ t
  val residue_class : _ t -> Integer.t Polynomial.t
  val residue_class_prime : prime t -> Integer.t
  val equal : _ t -> _ t -> bool
  val add : _ t -> _ t -> _ t
  val mul : _ t -> _ t -> _ t
  val pow : _ t -> Integer.t -> _ t
  val random : _ t -> _ t
  val zero : _ t -> _ t

  val extend :
    'a field t ->
    [< `Degree of int | `Quotient of 'a ring Polynomial.t ] ->
    'a field t
  (** extend the field {m K} of definition of {m a} by a root of the polynomial
     {m P\in K[X]} assumed to be irreducible over {m K}.  Return {m [r, m]} where {m r}
     is a root of {m P} in the extension field {m L} and {m m} is a map from {m K} to {m L},
     see [ffmap].
     If {m v} is given, the variable name is used to display the generator of {m L},
     else the name of the variable of {m P} is used.
     A generator of {m L} can be recovered using [b=ffgen(r)].
     The image of {m P} in {m L[X]} can be recovered using [PL=ffmap(m,P)]. *)

  val fpxq_star :
    p:pari_ulong -> quotient:Fp.t Polynomial.t -> _ finite_field group_structure

  val to_string : _ t -> string

  module Infix : sig
    val ( ~- ) : _ t -> _ t
    val ( + ) : _ t -> _ t -> _ t
    val ( - ) : _ t -> _ t -> _ t
    val ( * ) : _ t -> _ t -> _ t
    val ( ^ ) : _ t -> Integer.t -> _ t
  end
end

module Number_field : sig
  type number_field
  type structure
  type nonrec t = number_field field t

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
  type 'a p := 'a t
  type 'a structure constraint 'a = 'b field t
  type nonrec 'a t = elliptic_curve group t constraint 'a = 'b field t

  val create :
    ?a1:'a ->
    ?a2:'a ->
    ?a3:'a ->
    ?a4:'a ->
    ?a6:'a ->
    unit ->
    'a structure option
  (** [create ?a1 ?a2 ?a3 ?a4 ?a6] defines the curve
        {%math: Y^2 + a_1 XY + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6%}.
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
    _ Finite_field.t structure ->
    l:Integer.t ->
    p:_ Finite_field.t t ->
    q:_ Finite_field.t t ->
    _ Finite_field.t
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
    'a field p structure -> l:Signed.Long.t -> 'a ring p Polynomial.t
  (** {@ocaml[
      # let g = Finite_field.generator ~order:(Integer.of_int 625) (* 5^4 *);;
      val g : Finite_field.t = <abstr>
      # let ell = Option.get (Elliptic_curve.create ~a6:(Finite_field.pow g (Integer.of_int 6)) ());;
      val ell : Finite_field.t Elliptic_curve.structure = <abstr>
      # let pdiv7 = (Elliptic_curve.l_division_polynomial ell ~l:(Signed.Long.of_int 7));;
      val pdiv7 : Finite_field.finite_field ring t Polynomial.t = <abstr>
      # Polynomial.to_string pdiv7;;
      - : string =
      "2*x^24 + (3*x^3 + x + 2)*x^21 + (x^3 + x^2 + x + 2)*x^18 + (2*x^3 + 2*x^2 + 4*x)*x^15 + (2*x^3 + 4*x^2 + 4*x + 1)*x^12 + (3*x^3 + 4*x^2 + 1)*x^9 + (4*x^3 + x^2 + 4)*x^6 + (2*x^3 + 3*x^2 + 3*x + 4)*x^3 + (x^3 + x^2)"
      # Polynomial.degree pdiv7 = Signed.Long.of_int ((7 * 7 - 1) / 2);;
      - : bool = true
      ]} *)

  val to_string : 'a t -> string
  val add : 'a structure -> 'a t -> 'a t -> 'a t
  val sub : 'a structure -> 'a t -> 'a t -> 'a t
  val mul : 'a structure -> n:Integer.t -> p:'a t -> 'a t
  val equal : 'a t -> 'a t -> bool

  val generators_ff :
    _ Finite_field.t structure -> (_ Finite_field.t t, [ `ROW ]) Vector.t

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
val forprime_t_bb : (_ t, forprime_t Ctypes.structure) Ctypes.field
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
val forprime_t_pp : (_ t, forprime_t Ctypes.structure) Ctypes.field
val forcomposite_t : forcomposite_t Ctypes.structure Ctypes.typ
val forcomposite_t_first : (int, forcomposite_t Ctypes.structure) Ctypes.field
val forcomposite_t_b : (_ t, forcomposite_t Ctypes.structure) Ctypes.field
val forcomposite_t_n : (_ t, forcomposite_t Ctypes.structure) Ctypes.field
val forcomposite_t_p : (_ t, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_T :
  (forprime_t Ctypes.structure, forcomposite_t Ctypes.structure) Ctypes.field

val forvec_t : forvec_t Ctypes.structure Ctypes.typ
val forvec_t_first : (Signed.long, forvec_t Ctypes.structure) Ctypes.field
val forvec_t_a : (_ t Ctypes_static.ptr, forvec_t Ctypes.structure) Ctypes.field
val forvec_t_m : (_ t Ctypes_static.ptr, forvec_t Ctypes.structure) Ctypes.field
val forvec_t_M : (_ t Ctypes_static.ptr, forvec_t Ctypes.structure) Ctypes.field
val forvec_t_n : (Signed.long, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_next :
  ( (forvec_t Ctypes.structure Ctypes_static.ptr -> _ t)
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
val forpart_t_v : (_ t, forpart_t Ctypes.structure) Ctypes.field
val forperm_t : forperm_t Ctypes.structure Ctypes.typ
val forperm_t_k : (Signed.long, forperm_t Ctypes.structure) Ctypes.field
val forperm_t_first : (Signed.long, forperm_t Ctypes.structure) Ctypes.field
val forperm_t_v : (_ t, forperm_t Ctypes.structure) Ctypes.field
val forsubset_t : forsubset_t Ctypes.structure Ctypes.typ
val forsubset_t_n : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_k : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_all : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_first : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_v : (_ t, forsubset_t Ctypes.structure) Ctypes.field
val pari_plot : pari_plot Ctypes.structure Ctypes.typ

val pari_plot_draw :
  ( (pari_plot Ctypes.structure Ctypes_static.ptr -> _ t -> _ t -> _ t -> unit)
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
val genbin_x : (_ t, genbin Ctypes.structure) Ctypes.field
val genbin_base : (_ t, genbin Ctypes.structure) Ctypes.field

val genbin_rebase :
  ( (_ t -> Signed.long -> unit) Ctypes_static.static_funptr,
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
  (_ t, pari_parsestate Ctypes.structure) Ctypes.field

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
  (_ t, pari_global_state Ctypes.structure) Ctypes.field

val pari_global_state_seadata :
  (_ t, pari_global_state Ctypes.structure) Ctypes.field

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

val pari_thread_data : (_ t, pari_thread Ctypes.structure) Ctypes.field
val mt_state : mt_state Ctypes.structure Ctypes.typ
val mt_state_worker : (_ t, mt_state Ctypes.structure) Ctypes.field
val mt_state_pending : (_ t, mt_state Ctypes.structure) Ctypes.field
val mt_state_workid : (Signed.long, mt_state Ctypes.structure) Ctypes.field
val pari_mt : pari_mt Ctypes.structure Ctypes.typ

val pari_mt_mt :
  (mt_state Ctypes.structure, pari_mt Ctypes.structure) Ctypes.field

val pari_mt_get :
  ( (mt_state Ctypes.structure Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    _ t)
    Ctypes_static.static_funptr,
    pari_mt Ctypes.structure )
  Ctypes.field

val pari_mt_submit :
  ( (mt_state Ctypes.structure Ctypes_static.ptr -> Signed.long -> _ t -> unit)
    Ctypes_static.static_funptr,
    pari_mt Ctypes.structure )
  Ctypes.field

val pari_mt_end :
  (unit Ctypes_static.static_funptr, pari_mt Ctypes.structure) Ctypes.field

val parfor_iter : parfor_iter Ctypes.structure Ctypes.typ

val parfor_iter_pending :
  (Signed.long, parfor_iter Ctypes.structure) Ctypes.field

val parfor_iter_worker : (_ t, parfor_iter Ctypes.structure) Ctypes.field

val parfor_iter_pt :
  (pari_mt Ctypes.structure, parfor_iter Ctypes.structure) Ctypes.field

val parfor_t : parfor_t Ctypes.structure Ctypes.typ
val parfor_t_a : (_ t, parfor_t Ctypes.structure) Ctypes.field
val parfor_t_b : (_ t, parfor_t Ctypes.structure) Ctypes.field

val parfor_t_iter :
  (parfor_iter Ctypes.structure, parfor_t Ctypes.structure) Ctypes.field

val parforeach_t : parforeach_t Ctypes.structure Ctypes.typ
val parforeach_t_x : (_ t, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_W : (_ t, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_i : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_l : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field

val parforeach_t_iter :
  (parfor_iter Ctypes.structure, parforeach_t Ctypes.structure) Ctypes.field

val parforprime_t : parforprime_t Ctypes.structure Ctypes.typ
val parforprime_t_v : (_ t, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_forprime :
  (forprime_t Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_iter :
  (parfor_iter Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforvec_t : parforvec_t Ctypes.structure Ctypes.typ
val parforvec_t_v : (_ t, parforvec_t Ctypes.structure) Ctypes.field

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
val nfmaxord_t_T : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dT : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_T0 : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_unscale : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dK : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_index : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_basis : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_r1 : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_basden : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dTP : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dTE : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dKP : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dKE : (_ t, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_certify : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field
val qfr_data : qfr_data Ctypes.structure Ctypes.typ
val qfr_data_D : (_ t, qfr_data Ctypes.structure) Ctypes.field
val qfr_data_sqrtD : (_ t, qfr_data Ctypes.structure) Ctypes.field
val qfr_data_isqrtD : (_ t, qfr_data Ctypes.structure) Ctypes.field
val fp_chk_fun : fp_chk_fun Ctypes.structure Ctypes.typ

val fp_chk_fun_f :
  ( (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_init :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr -> _ t -> _ t -> _ t)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_post :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr -> _ t -> _ t -> _ t)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_data :
  (unit Ctypes_static.ptr, fp_chk_fun Ctypes.structure) Ctypes.field

val fp_chk_fun_skipfirst :
  (Signed.long, fp_chk_fun Ctypes.structure) Ctypes.field

val zlog_s : zlog_s Ctypes.structure Ctypes.typ
val zlog_s_bid : (_ t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_P : (_ t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_k : (_ t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_sprk : (_ t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_archp : (_ t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_mod : (_ t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_U : (_ t, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_hU : (Signed.long, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_no2 : (int, zlog_s Ctypes.structure) Ctypes.field
val bb_group : bb_group Ctypes.structure Ctypes.typ

val bb_group_mul :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_pow :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_rand :
  ( (unit Ctypes_static.ptr -> _ t) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_hash :
  ( (_ t -> pari_ulong) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal :
  ( (_ t -> _ t -> int) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal1 :
  ( (_ t -> int) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_easylog :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t -> _ t)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_field : bb_field Ctypes.structure Ctypes.typ

val bb_field_red :
  ( (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_add :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_mul :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_neg :
  ( (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_inv :
  ( (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_equal0 :
  ( (_ t -> int) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_s :
  ( (unit Ctypes_static.ptr -> Signed.long -> _ t) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_algebra : bb_algebra Ctypes.structure Ctypes.typ

val bb_algebra_red :
  ( (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_add :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sub :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_mul :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sqr :
  ( (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_one :
  ( (unit Ctypes_static.ptr -> _ t) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_zero :
  ( (unit Ctypes_static.ptr -> _ t) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_ring : bb_ring Ctypes.structure Ctypes.typ

val bb_ring_add :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_mul :
  ( (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_sqr :
  ( (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val buchimag : _ t -> _ t -> _ t -> _ t -> _ t
val buchreal : _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zidealstar : _ t -> _ t -> _ t
val zidealstarinit : _ t -> _ t -> _ t
val zidealstarinitgen : _ t -> _ t -> _ t
val factmod : _ t -> _ t -> _ t
val mpbern : Signed.long -> Signed.long -> unit
val simplefactmod : _ t -> _ t -> _ t
val listkill : _ t -> unit
val isprincipalforce : _ t -> _ t -> _ t
val isprincipalgen : _ t -> _ t -> _ t
val isprincipalgenforce : _ t -> _ t -> _ t
val f2ms_ker : _ t -> Signed.long -> _ t
val f2ms_to_f2m : _ t -> Signed.long -> _ t
val f2c_to_zc : _ t -> _ t
val f2c_to_mod : _ t -> _ t
val f2m_f2c_gauss : _ t -> _ t -> _ t
val f2m_f2c_invimage : _ t -> _ t -> _ t
val f2m_f2c_mul : _ t -> _ t -> _ t
val f2m_deplin : _ t -> _ t
val f2m_det : _ t -> pari_ulong
val f2m_det_sp : _ t -> pari_ulong
val f2m_gauss : _ t -> _ t -> _ t
val f2m_inv : _ t -> _ t
val f2m_invimage : _ t -> _ t -> _ t
val f2m_ker : _ t -> _ t
val f2m_ker_sp : _ t -> Signed.long -> _ t
val f2m_mul : _ t -> _ t -> _ t
val f2m_powu : _ t -> pari_ulong -> _ t
val f2m_rank : _ t -> Signed.long
val f2m_row : _ t -> Signed.long -> _ t
val f2m_rowslice : _ t -> Signed.long -> Signed.long -> _ t
val f2m_to_f2ms : _ t -> _ t
val f2m_to_flm : _ t -> _ t
val f2m_to_zm : _ t -> _ t
val f2m_to_mod : _ t -> _ t
val f2m_transpose : _ t -> _ t
val f2v_add_inplace : _ t -> _ t -> unit
val f2v_and_inplace : _ t -> _ t -> unit
val f2v_dotproduct : _ t -> _ t -> pari_ulong
val f2v_equal0 : _ t -> int
val f2v_hamming : _ t -> pari_ulong
val f2v_negimply_inplace : _ t -> _ t -> unit
val f2v_or_inplace : _ t -> _ t -> unit
val f2v_slice : _ t -> Signed.long -> Signed.long -> _ t
val f2v_subset : _ t -> _ t -> int
val f2v_to_flv : _ t -> _ t
val matid_f2m : Signed.long -> _ t
val f2x_f2xq_eval : _ t -> _ t -> _ t -> _ t
val f2x_f2xqv_eval : _ t -> _ t -> _ t -> _ t
val f2x_frobenius : _ t -> _ t
val f2x_1_add : _ t -> _ t
val f2x_add : _ t -> _ t -> _ t
val f2x_deflate : _ t -> Signed.long -> _ t
val f2x_degfact : _ t -> _ t
val f2x_degree : _ t -> Signed.long
val f2x_deriv : _ t -> _ t
val f2x_divrem : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val f2x_eval : _ t -> pari_ulong -> pari_ulong
val f2x_even_odd : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit

val f2x_extgcd :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val f2x_gcd : _ t -> _ t -> _ t
val f2x_get_red : _ t -> _ t
val f2x_halfgcd : _ t -> _ t -> _ t
val f2x_issquare : _ t -> int
val f2x_matfrobenius : _ t -> _ t
val f2x_mul : _ t -> _ t -> _ t
val f2x_recip : _ t -> _ t
val f2x_rem : _ t -> _ t -> _ t
val f2x_shift : _ t -> Signed.long -> _ t
val f2x_sqr : _ t -> _ t
val f2x_sqrt : _ t -> _ t
val f2x_to_f2v : _ t -> Signed.long -> _ t
val f2x_to_f2xx : _ t -> Signed.long -> _ t
val f2x_to_flx : _ t -> _ t
val f2x_to_zx : _ t -> _ t
val f2x_valrem : _ t -> _ t Ctypes_static.ptr -> Signed.long
val f2xc_to_flxc : _ t -> _ t
val f2xc_to_zxc : _ t -> _ t
val f2xv_to_f2m : _ t -> Signed.long -> _ t
val f2xv_to_flxv_inplace : _ t -> unit
val f2xv_to_zxv_inplace : _ t -> unit
val f2xx_f2x_add : _ t -> _ t -> _ t
val f2xx_f2x_mul : _ t -> _ t -> _ t
val f2xx_add : _ t -> _ t -> _ t
val f2xx_deriv : _ t -> _ t
val f2xx_renormalize : _ t -> Signed.long -> _ t
val f2xx_to_kronecker : _ t -> Signed.long -> _ t
val f2xx_to_flxx : _ t -> _ t
val f2xx_to_zxx : _ t -> _ t
val f2xx_to_f2xc : _ t -> Signed.long -> Signed.long -> _ t
val f2xxv_to_f2xm : _ t -> Signed.long -> Signed.long -> _ t
val f2xxc_to_zxxc : _ t -> _ t
val f2xy_f2xq_evalx : _ t -> _ t -> _ t -> _ t
val f2xy_f2xqv_evalx : _ t -> _ t -> _ t -> _ t
val f2xy_degreex : _ t -> Signed.long
val f2xn_div : _ t -> _ t -> Signed.long -> _ t
val f2xn_inv : _ t -> Signed.long -> _ t
val f2xn_red : _ t -> Signed.long -> _ t
val f2xq_artin_schreier : _ t -> _ t -> _ t
val f2xq_autpow : _ t -> Signed.long -> _ t -> _ t
val f2xq_conjvec : _ t -> _ t -> _ t
val f2xq_div : _ t -> _ t -> _ t -> _ t
val f2xq_inv : _ t -> _ t -> _ t
val f2xq_invsafe : _ t -> _ t -> _ t
val f2xq_log : _ t -> _ t -> _ t -> _ t -> _ t
val f2xq_matrix_pow : _ t -> Signed.long -> Signed.long -> _ t -> _ t
val f2xq_mul : _ t -> _ t -> _ t -> _ t
val f2xq_order : _ t -> _ t -> _ t -> _ t
val f2xq_pow : _ t -> _ t -> _ t -> _ t
val f2xq_pow_init : _ t -> _ t -> Signed.long -> _ t -> _ t
val f2xq_pow_table : _ t -> _ t -> _ t -> _ t
val f2xq_powu : _ t -> pari_ulong -> _ t -> _ t
val f2xq_powers : _ t -> Signed.long -> _ t -> _ t
val f2xq_sqr : _ t -> _ t -> _ t
val f2xq_sqrt : _ t -> _ t -> _ t
val f2xq_sqrt_fast : _ t -> _ t -> _ t -> _ t
val f2xq_sqrtn : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val f2xq_trace : _ t -> _ t -> pari_ulong
val f2xqx_f2xq_mul : _ t -> _ t -> _ t -> _ t
val f2xqx_f2xq_mul_to_monic : _ t -> _ t -> _ t -> _ t
val f2xqx_f2xqxq_eval : _ t -> _ t -> _ t -> _ t -> _ t
val f2xqx_f2xqxqv_eval : _ t -> _ t -> _ t -> _ t -> _ t
val f2xqx_disc : _ t -> _ t -> _ t
val f2xqx_divrem : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t

val f2xqx_extgcd :
  _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val f2xqx_gcd : _ t -> _ t -> _ t -> _ t
val f2xqx_get_red : _ t -> _ t -> _ t
val f2xqx_halfgcd : _ t -> _ t -> _ t -> _ t

val f2xqx_halfgcd_all :
  _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val f2xqx_invbarrett : _ t -> _ t -> _ t

val f2xqx_ispower :
  _ t -> Signed.long -> _ t -> _ t Ctypes_static.ptr -> Signed.long

val f2xqx_mul : _ t -> _ t -> _ t -> _ t
val f2xqx_normalize : _ t -> _ t -> _ t
val f2xqx_powu : _ t -> pari_ulong -> _ t -> _ t
val f2xqx_red : _ t -> _ t -> _ t
val f2xqx_rem : _ t -> _ t -> _ t -> _ t
val f2xqx_resultant : _ t -> _ t -> _ t -> _ t
val f2xqx_sqr : _ t -> _ t -> _ t
val f2xqxq_inv : _ t -> _ t -> _ t -> _ t
val f2xqxq_invsafe : _ t -> _ t -> _ t -> _ t
val f2xqxq_mul : _ t -> _ t -> _ t -> _ t -> _ t
val f2xqxq_sqr : _ t -> _ t -> _ t -> _ t
val f2xqxq_pow : _ t -> _ t -> _ t -> _ t -> _ t
val f2xqxq_powers : _ t -> Signed.long -> _ t -> _ t -> _ t
val f2xqxq_autpow : _ t -> Signed.long -> _ t -> _ t -> _ t
val f2xqxq_auttrace : _ t -> Signed.long -> _ t -> _ t -> _ t
val f2xqxqv_red : _ t -> _ t -> _ t -> _ t
val flm_to_f2m : _ t -> _ t
val flv_to_f2v : _ t -> _ t
val flx_to_f2x : _ t -> _ t
val flxc_to_f2xc : _ t -> _ t
val flxx_to_f2xx : _ t -> _ t
val flxxc_to_f2xxc : _ t -> _ t
val kronecker_to_f2xqx : _ t -> _ t -> _ t
val rg_to_f2xq : _ t -> _ t -> _ t
val rgm_to_f2m : _ t -> _ t
val rgv_to_f2v : _ t -> _ t
val rgx_to_f2x : _ t -> _ t
val z_to_f2x : _ t -> Signed.long -> _ t
val zm_to_f2m : _ t -> _ t
val zv_to_f2v : _ t -> _ t
val zx_to_f2x : _ t -> _ t
val zxx_to_f2xx : _ t -> Signed.long -> _ t
val const_f2v : Signed.long -> _ t
val gener_f2xq : _ t -> _ t Ctypes_static.ptr -> _ t

val get_f2xq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  bb_field Ctypes.structure Ctypes_static.ptr

val monomial_f2x : Signed.long -> Signed.long -> _ t
val pol1_f2xx : Signed.long -> Signed.long -> _ t
val polx_f2xx : Signed.long -> Signed.long -> _ t
val random_f2xqx : Signed.long -> Signed.long -> _ t -> _ t
val f2x_teichmuller : _ t -> Signed.long -> _ t
val f2xq_ellcard : _ t -> _ t -> _ t -> _ t
val f2xq_ellgens : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val f2xq_ellgroup : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t

val f2xq_elltwist :
  _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit

val f2xqe_add : _ t -> _ t -> _ t -> _ t -> _ t
val f2xqe_changepoint : _ t -> _ t -> _ t -> _ t
val f2xqe_changepointinv : _ t -> _ t -> _ t -> _ t
val f2xqe_dbl : _ t -> _ t -> _ t -> _ t
val f2xqe_log : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val f2xqe_mul : _ t -> _ t -> _ t -> _ t -> _ t
val f2xqe_neg : _ t -> _ t -> _ t -> _ t
val f2xqe_order : _ t -> _ t -> _ t -> _ t -> _ t
val f2xqe_sub : _ t -> _ t -> _ t -> _ t -> _ t
val f2xqe_tatepairing : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val f2xqe_weilpairing : _ t -> _ t -> _ t -> _ t -> _ t -> _ t

val get_f2xqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  _ t ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_f2xqe : _ t -> _ t -> _ t
val random_f2xqe : _ t -> _ t -> _ t -> _ t
val f3c_to_mod : _ t -> _ t
val f3c_to_zc : _ t -> _ t
val f3m_ker : _ t -> _ t
val f3m_ker_sp : _ t -> Signed.long -> _ t
val f3m_mul : _ t -> _ t -> _ t
val f3m_row : _ t -> Signed.long -> _ t
val f3m_to_flm : _ t -> _ t
val f3m_to_zm : _ t -> _ t
val f3m_to_mod : _ t -> _ t
val f3m_transpose : _ t -> _ t
val f3v_to_flv : _ t -> _ t
val f3v_coeff : _ t -> Signed.long -> pari_ulong
val f3v_clear : _ t -> Signed.long -> unit
val f3v_set : _ t -> Signed.long -> pari_ulong -> unit
val flm_to_f3m : _ t -> _ t
val flv_to_f3v : _ t -> _ t
val rgm_to_f3m : _ t -> _ t
val rgv_to_f3v : _ t -> _ t
val zm_to_f3m : _ t -> _ t
val zv_to_f3v : _ t -> _ t
val zero_f3m_copy : Signed.long -> Signed.long -> _ t
val zero_f3v : Signed.long -> _ t
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
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> _ t

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

val fle_add : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val fle_dbl : _ t -> pari_ulong -> pari_ulong -> _ t
val fle_changepoint : _ t -> _ t -> pari_ulong -> _ t
val fle_changepointinv : _ t -> _ t -> pari_ulong -> _ t
val fle_log : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val fle_mul : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val fle_mulu : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val fle_order : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val fle_sub : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val fle_tatepairing :
  _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fle_to_flj : _ t -> _ t

val fle_weilpairing :
  _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val flj_add_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val flj_changepointinv_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flj_dbl_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t

val flj_mulu_pre :
  _ t -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> _ t

val flj_neg : _ t -> pari_ulong -> _ t
val flj_to_fle : _ t -> pari_ulong -> _ t
val flj_to_fle_pre : _ t -> pari_ulong -> pari_ulong -> _ t

val fljv_factorback_pre :
  _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t

val random_fle : pari_ulong -> pari_ulong -> pari_ulong -> _ t
val random_fle_pre : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val random_flj_pre : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val flc_to_zc : _ t -> _ t
val flc_to_zc_inplace : _ t -> _ t
val flm_flc_gauss : _ t -> _ t -> pari_ulong -> _ t
val flm_flc_invimage : _ t -> _ t -> pari_ulong -> _ t
val flm_adjoint : _ t -> pari_ulong -> _ t
val flm_deplin : _ t -> pari_ulong -> _ t
val flm_det : _ t -> pari_ulong -> pari_ulong
val flm_det_sp : _ t -> pari_ulong -> pari_ulong
val flm_gauss : _ t -> _ t -> pari_ulong -> _ t
val flm_intersect : _ t -> _ t -> pari_ulong -> _ t
val flm_intersect_i : _ t -> _ t -> pari_ulong -> _ t
val flm_inv : _ t -> pari_ulong -> _ t
val flm_invimage : _ t -> _ t -> pari_ulong -> _ t
val flm_ker : _ t -> pari_ulong -> _ t
val flm_ker_sp : _ t -> pari_ulong -> Signed.long -> _ t
val flm_rank : _ t -> pari_ulong -> Signed.long
val flm_to_zm : _ t -> _ t
val flm_to_zm_inplace : _ t -> _ t
val flv_to_zv : _ t -> _ t
val fl_to_flx : pari_ulong -> Signed.long -> _ t
val fl2_equal1 : _ t -> int
val fl2_inv_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val fl2_mul_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val fl2_norm_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val fl2_pow_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val fl2_sqr_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val fl2_sqrt_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t

val fl2_sqrtn_pre :
  _ t ->
  _ t ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t

val flm_to_flxv : _ t -> Signed.long -> _ t
val flm_to_flxx : _ t -> Signed.long -> Signed.long -> _ t
val flv_flm_polint : _ t -> _ t -> pari_ulong -> Signed.long -> _ t
val flv_inv : _ t -> pari_ulong -> _ t
val flv_inv_inplace : _ t -> pari_ulong -> unit
val flv_inv_pre_inplace : _ t -> pari_ulong -> pari_ulong -> unit
val flv_inv_pre : _ t -> pari_ulong -> pari_ulong -> _ t
val flv_invvandermonde : _ t -> pari_ulong -> pari_ulong -> _ t
val flv_polint : _ t -> _ t -> pari_ulong -> Signed.long -> _ t
val flv_prod : _ t -> pari_ulong -> pari_ulong
val flv_prod_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong
val flv_roots_to_pol : _ t -> pari_ulong -> Signed.long -> _ t
val flv_to_flx : _ t -> Signed.long -> _ t
val flx_fl_add : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_fl_mul : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_fl_mul_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val flx_fl_mul_to_monic : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_fl_sub : _ t -> pari_ulong -> pari_ulong -> _ t

val flx_fl2_eval_pre :
  _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t

val flx_flv_multieval : _ t -> _ t -> pari_ulong -> _ t
val flx_flxq_eval : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flx_flxq_eval_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flx_flxqv_eval : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flx_flxqv_eval_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flx_frobenius : _ t -> pari_ulong -> _ t
val flx_frobenius_pre : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_laplace : _ t -> pari_ulong -> _ t
val flx_newton : _ t -> Signed.long -> pari_ulong -> _ t
val flx_add : _ t -> _ t -> pari_ulong -> _ t
val flx_blocks : _ t -> Signed.long -> Signed.long -> _ t
val flx_composedprod : _ t -> _ t -> pari_ulong -> _ t
val flx_composedsum : _ t -> _ t -> pari_ulong -> _ t
val flx_convol : _ t -> _ t -> pari_ulong -> _ t
val flx_deflate : _ t -> Signed.long -> _ t
val flx_deriv : _ t -> pari_ulong -> _ t
val flx_diff1 : _ t -> pari_ulong -> _ t
val flx_digits : _ t -> _ t -> pari_ulong -> _ t

val flx_div_by_x_x :
  _ t -> pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> _ t

val flx_divrem : _ t -> _ t -> pari_ulong -> _ t Ctypes_static.ptr -> _ t

val flx_divrem_pre :
  _ t -> _ t -> pari_ulong -> pari_ulong -> _ t Ctypes_static.ptr -> _ t

val flx_double : _ t -> pari_ulong -> _ t
val flx_equal : _ t -> _ t -> int
val flx_eval : _ t -> pari_ulong -> pari_ulong -> pari_ulong
val flx_eval_powers_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong
val flx_eval_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val flx_extgcd :
  _ t ->
  _ t ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val flx_extgcd_pre :
  _ t ->
  _ t ->
  pari_ulong ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val flx_extresultant :
  _ t ->
  _ t ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  pari_ulong

val flx_extresultant_pre :
  _ t ->
  _ t ->
  pari_ulong ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  pari_ulong

val flx_fromnewton : _ t -> pari_ulong -> _ t
val flx_gcd : _ t -> _ t -> pari_ulong -> _ t
val flx_gcd_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flx_get_red : _ t -> pari_ulong -> _ t
val flx_get_red_pre : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_halfgcd : _ t -> _ t -> pari_ulong -> _ t

val flx_halfgcd_all :
  _ t ->
  _ t ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val flx_halfgcd_all_pre :
  _ t ->
  _ t ->
  pari_ulong ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val flx_halfgcd_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flx_halve : _ t -> pari_ulong -> _ t
val flx_inflate : _ t -> Signed.long -> _ t
val flx_integ : _ t -> pari_ulong -> _ t
val flx_invbarrett : _ t -> pari_ulong -> _ t
val flx_invlaplace : _ t -> pari_ulong -> _ t
val flx_is_squarefree : _ t -> pari_ulong -> int
val flx_is_smooth : _ t -> Signed.long -> pari_ulong -> int
val flx_is_smooth_pre : _ t -> Signed.long -> pari_ulong -> pari_ulong -> int
val flx_matfrobenius : _ t -> pari_ulong -> _ t
val flx_matfrobenius_pre : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_mod_xn1 : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_mod_xnm1 : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_mul : _ t -> _ t -> pari_ulong -> _ t
val flx_mul_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flx_neg : _ t -> pari_ulong -> _ t
val flx_neg_inplace : _ t -> pari_ulong -> _ t
val flx_normalize : _ t -> pari_ulong -> _ t
val flx_powu : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_powu_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val flx_recip : _ t -> _ t
val flx_red : _ t -> pari_ulong -> _ t
val flx_rem : _ t -> _ t -> pari_ulong -> _ t
val flx_rem_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flx_renormalize : _ t -> Signed.long -> _ t
val flx_rescale : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_resultant : _ t -> _ t -> pari_ulong -> pari_ulong
val flx_resultant_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong
val flx_shift : _ t -> Signed.long -> _ t
val flx_splitting : _ t -> Signed.long -> _ t
val flx_sqr : _ t -> pari_ulong -> _ t
val flx_sqr_pre : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_sub : _ t -> _ t -> pari_ulong -> _ t
val flx_translate1 : _ t -> pari_ulong -> _ t
val flx_translate1_basecase : _ t -> pari_ulong -> _ t
val flx_to_flv : _ t -> Signed.long -> _ t
val flx_to_flxx : _ t -> Signed.long -> _ t
val flx_to_zx : _ t -> _ t
val flx_to_zx_inplace : _ t -> _ t
val flx_triple : _ t -> pari_ulong -> _ t
val flx_val : _ t -> Signed.long
val flx_valrem : _ t -> _ t Ctypes_static.ptr -> Signed.long
val flxc_flxqv_eval : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxc_flxqv_eval_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxc_flxq_eval : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxc_flxq_eval_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxc_eval_powers_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxc_neg : _ t -> pari_ulong -> _ t
val flxc_sub : _ t -> _ t -> pari_ulong -> _ t
val flxc_to_zxc : _ t -> _ t
val flxm_flx_add_shallow : _ t -> _ t -> pari_ulong -> _ t
val flxm_eval_powers_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxm_neg : _ t -> pari_ulong -> _ t
val flxm_sub : _ t -> _ t -> pari_ulong -> _ t
val flxm_to_flxxv : _ t -> Signed.long -> _ t
val flxm_to_zxm : _ t -> _ t
val flxt_red : _ t -> pari_ulong -> _ t
val flxv_flc_mul : _ t -> _ t -> pari_ulong -> _ t
val flxv_flv_multieval : _ t -> _ t -> pari_ulong -> _ t
val flxv_flx_fromdigits : _ t -> _ t -> pari_ulong -> _ t
val flxv_composedsum : _ t -> pari_ulong -> _ t
val flxv_prod : _ t -> pari_ulong -> _ t
val flxv_red : _ t -> pari_ulong -> _ t
val flxv_to_flm : _ t -> Signed.long -> _ t
val flxv_to_flxx : _ t -> Signed.long -> _ t
val flxv_to_zxv : _ t -> _ t
val flxv_to_zxv_inplace : _ t -> unit
val flxn_div : _ t -> _ t -> Signed.long -> pari_ulong -> _ t
val flxn_div_pre : _ t -> _ t -> Signed.long -> pari_ulong -> pari_ulong -> _ t
val flxn_exp : _ t -> Signed.long -> pari_ulong -> _ t
val flxn_expint : _ t -> Signed.long -> pari_ulong -> _ t
val flxn_inv : _ t -> Signed.long -> pari_ulong -> _ t
val flxn_mul : _ t -> _ t -> Signed.long -> pari_ulong -> _ t
val flxn_mul_pre : _ t -> _ t -> Signed.long -> pari_ulong -> pari_ulong -> _ t
val flxn_sqr : _ t -> Signed.long -> pari_ulong -> _ t
val flxn_sqr_pre : _ t -> Signed.long -> pari_ulong -> pari_ulong -> _ t
val flxn_red : _ t -> Signed.long -> _ t
val flxq_autpow : _ t -> pari_ulong -> _ t -> pari_ulong -> _ t

val flxq_autpow_pre :
  _ t -> pari_ulong -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxq_autpowers : _ t -> pari_ulong -> _ t -> pari_ulong -> _ t
val flxq_autsum : _ t -> pari_ulong -> _ t -> pari_ulong -> _ t
val flxq_auttrace : _ t -> pari_ulong -> _ t -> pari_ulong -> _ t

val flxq_auttrace_pre :
  _ t -> pari_ulong -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxq_charpoly : _ t -> _ t -> pari_ulong -> _ t
val flxq_conjvec : _ t -> _ t -> pari_ulong -> _ t
val flxq_div : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxq_div_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_inv : _ t -> _ t -> pari_ulong -> _ t
val flxq_inv_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_invsafe : _ t -> _ t -> pari_ulong -> _ t
val flxq_invsafe_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_issquare : _ t -> _ t -> pari_ulong -> int
val flxq_is2npower : _ t -> Signed.long -> _ t -> pari_ulong -> int
val flxq_log : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxq_lroot : _ t -> _ t -> Signed.long -> _ t
val flxq_lroot_pre : _ t -> _ t -> Signed.long -> pari_ulong -> _ t
val flxq_lroot_fast : _ t -> _ t -> _ t -> Signed.long -> _ t
val flxq_lroot_fast_pre : _ t -> _ t -> _ t -> Signed.long -> pari_ulong -> _ t

val flxq_matrix_pow :
  _ t -> Signed.long -> Signed.long -> _ t -> pari_ulong -> _ t

val flxq_matrix_pow_pre :
  _ t -> Signed.long -> Signed.long -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxq_minpoly : _ t -> _ t -> pari_ulong -> _ t
val flxq_minpoly_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_mul : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxq_mul_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_norm : _ t -> _ t -> pari_ulong -> pari_ulong
val flxq_order : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxq_pow : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxq_pow_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_pow_init : _ t -> _ t -> Signed.long -> _ t -> pari_ulong -> _ t

val flxq_pow_init_pre :
  _ t -> _ t -> Signed.long -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxq_pow_table_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_pow_table : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxq_powu : _ t -> pari_ulong -> _ t -> pari_ulong -> _ t
val flxq_powu_pre : _ t -> pari_ulong -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_powers : _ t -> Signed.long -> _ t -> pari_ulong -> _ t

val flxq_powers_pre :
  _ t -> Signed.long -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxq_sqr : _ t -> _ t -> pari_ulong -> _ t
val flxq_sqr_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_sqrt : _ t -> _ t -> pari_ulong -> _ t
val flxq_sqrt_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxq_sqrtn : _ t -> _ t -> _ t -> pari_ulong -> _ t Ctypes_static.ptr -> _ t
val flxq_trace : _ t -> _ t -> pari_ulong -> pari_ulong
val flxqc_flxq_mul : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqm_flxq_mul : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqv_dotproduct : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqv_dotproduct_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val rg_to_f2 : _ t -> pari_ulong
val rg_to_fl : _ t -> pari_ulong -> pari_ulong
val rg_to_flxq : _ t -> _ t -> pari_ulong -> _ t
val rgx_to_flx : _ t -> pari_ulong -> _ t
val rgxv_to_flxv : _ t -> pari_ulong -> _ t
val z_to_flx : _ t -> pari_ulong -> Signed.long -> _ t
val zxv_to_flxv : _ t -> pari_ulong -> _ t
val zxt_to_flxt : _ t -> pari_ulong -> _ t
val gener_flxq : _ t -> pari_ulong -> _ t Ctypes_static.ptr -> _ t

val get_flxq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  pari_ulong ->
  bb_field Ctypes.structure Ctypes_static.ptr

val get_flxq_star :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  pari_ulong ->
  bb_group Ctypes.structure Ctypes_static.ptr

val monomial_flx : pari_ulong -> Signed.long -> Signed.long -> _ t
val random_flx : Signed.long -> Signed.long -> pari_ulong -> _ t
val zero_flxc : Signed.long -> Signed.long -> _ t
val zero_flxm : Signed.long -> Signed.long -> Signed.long -> _ t
val zlx_translate1 : _ t -> pari_ulong -> Signed.long -> _ t
val zx_to_flx : _ t -> pari_ulong -> _ t
val flxx_fl_mul : _ t -> pari_ulong -> pari_ulong -> _ t
val flxx_flx_add : _ t -> _ t -> pari_ulong -> _ t
val flxx_flx_mul : _ t -> _ t -> pari_ulong -> _ t
val flxx_flx_sub : _ t -> _ t -> pari_ulong -> _ t
val flxx_laplace : _ t -> pari_ulong -> _ t
val flxx_add : _ t -> _ t -> pari_ulong -> _ t
val flxx_blocks : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val flxx_deriv : _ t -> pari_ulong -> _ t
val flxx_double : _ t -> pari_ulong -> _ t
val flxx_invlaplace : _ t -> pari_ulong -> _ t
val flxx_neg : _ t -> pari_ulong -> _ t
val flxx_renormalize : _ t -> Signed.long -> _ t
val flxx_shift : _ t -> Signed.long -> Signed.long -> _ t
val flxx_sub : _ t -> _ t -> pari_ulong -> _ t
val flxx_swap : _ t -> Signed.long -> Signed.long -> _ t
val flxx_to_flm : _ t -> Signed.long -> _ t
val flxx_to_flx : _ t -> _ t
val flxx_to_flxc : _ t -> Signed.long -> Signed.long -> _ t
val flxx_to_zxx : _ t -> _ t
val flxx_translate1 : _ t -> Signed.long -> Signed.long -> _ t
val flxx_triple : _ t -> pari_ulong -> _ t
val flxxc_sub : _ t -> _ t -> pari_ulong -> _ t
val flxxc_to_zxxc : _ t -> _ t
val flxxm_to_zxxm : _ t -> _ t
val flxxv_to_flxm : _ t -> Signed.long -> Signed.long -> _ t
val flxxn_red : _ t -> Signed.long -> _ t
val flxy_flx_div : _ t -> _ t -> pari_ulong -> _ t
val flxy_flx_translate : _ t -> _ t -> pari_ulong -> _ t
val flxy_flxqv_evalx : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxy_flxqv_evalx_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxy_flxq_evalx : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxy_flxq_evalx_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxy_evalx : _ t -> pari_ulong -> pari_ulong -> _ t
val flxy_evalx_powers_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxy_evalx_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val flxyqq_pow : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqv_roots_to_pol : _ t -> _ t -> pari_ulong -> Signed.long -> _ t
val flxqxc_flxqxqv_eval : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t

val flxqxc_flxqxqv_eval_pre :
  _ t -> _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxc_flxqxq_eval : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqxq_autpow : _ t -> Signed.long -> _ t -> _ t -> pari_ulong -> _ t

val flxqxq_autpow_pre :
  _ t -> Signed.long -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxq_autsum : _ t -> Signed.long -> _ t -> _ t -> pari_ulong -> _ t

val flxqxq_autsum_pre :
  _ t -> Signed.long -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxq_auttrace : _ t -> pari_ulong -> _ t -> _ t -> pari_ulong -> _ t

val flxqxq_auttrace_pre :
  _ t -> pari_ulong -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxq_div : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqxq_div_pre : _ t -> _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqxq_inv : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqxq_inv_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqxq_invsafe : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqxq_invsafe_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxq_matrix_pow :
  _ t -> Signed.long -> Signed.long -> _ t -> _ t -> pari_ulong -> _ t

val flxqxq_minpoly : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqxq_minpoly_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqxq_mul : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqxq_mul_pre : _ t -> _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqxq_pow : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqxq_pow_pre : _ t -> _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqxq_powers : _ t -> Signed.long -> _ t -> _ t -> pari_ulong -> _ t

val flxqxq_powers_pre :
  _ t -> Signed.long -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxq_powu : _ t -> pari_ulong -> _ t -> _ t -> pari_ulong -> _ t

val flxqxq_powu_pre :
  _ t -> pari_ulong -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxq_sqr : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqxq_sqr_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqxv_prod : _ t -> _ t -> pari_ulong -> _ t
val flxqx_flxqxqv_eval : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t

val flxqx_flxqxqv_eval_pre :
  _ t -> _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqx_flxqxq_eval : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t

val flxqx_flxqxq_eval_pre :
  _ t -> _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqx_flxq_mul : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqx_flxq_mul_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_flxq_mul_to_monic : _ t -> _ t -> _ t -> pari_ulong -> _ t

val flxqx_flxq_mul_to_monic_pre :
  _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqx_newton : _ t -> Signed.long -> _ t -> pari_ulong -> _ t

val flxqx_newton_pre :
  _ t -> Signed.long -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqx_composedsum : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqx_disc : _ t -> _ t -> pari_ulong -> _ t

val flxqx_div_by_x_x :
  _ t -> _ t -> _ t -> pari_ulong -> _ t Ctypes_static.ptr -> _ t

val flxqx_div_by_x_x_pre :
  _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t Ctypes_static.ptr -> _ t

val flxqx_divrem :
  _ t -> _ t -> _ t -> pari_ulong -> _ t Ctypes_static.ptr -> _ t

val flxqx_divrem_pre :
  _ t -> _ t -> _ t -> pari_ulong -> Signed.long -> _ t Ctypes_static.ptr -> _ t

val flxqx_dotproduct : _ t -> _ t -> _ t -> pari_ulong -> _ t

val flxqx_extgcd :
  _ t ->
  _ t ->
  _ t ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val flxqx_extgcd_pre :
  _ t ->
  _ t ->
  _ t ->
  pari_ulong ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val flxqx_fromnewton : _ t -> _ t -> pari_ulong -> _ t
val flxqx_fromnewton_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_gcd : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqx_gcd_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_get_red : _ t -> _ t -> pari_ulong -> _ t
val flxqx_get_red_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_halfgcd : _ t -> _ t -> _ t -> pari_ulong -> _ t

val flxqx_halfgcd_all :
  _ t ->
  _ t ->
  _ t ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val flxqx_halfgcd_all_pre :
  _ t ->
  _ t ->
  _ t ->
  pari_ulong ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val flxqx_halfgcd_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_invbarrett : _ t -> _ t -> pari_ulong -> _ t
val flxqx_invbarrett_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_mul : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqx_mul_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_normalize : _ t -> _ t -> pari_ulong -> _ t
val flxqx_normalize_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_powu : _ t -> pari_ulong -> _ t -> pari_ulong -> _ t
val flxqx_powu_pre : _ t -> pari_ulong -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_red : _ t -> _ t -> pari_ulong -> _ t
val flxqx_red_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_rem : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqx_rem_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_resultant : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqx_resultant_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_safegcd : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqx_saferesultant : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqx_sqr : _ t -> _ t -> pari_ulong -> _ t
val flxqx_sqr_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqxn_expint : _ t -> Signed.long -> _ t -> pari_ulong -> _ t

val flxqxn_expint_pre :
  _ t -> Signed.long -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxn_inv : _ t -> Signed.long -> _ t -> pari_ulong -> _ t

val flxqxn_inv_pre :
  _ t -> Signed.long -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxn_mul : _ t -> _ t -> Signed.long -> _ t -> pari_ulong -> _ t

val flxqxn_mul_pre :
  _ t -> _ t -> Signed.long -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqxn_sqr : _ t -> Signed.long -> _ t -> pari_ulong -> _ t

val flxqxn_sqr_pre :
  _ t -> Signed.long -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxy_degreex : _ t -> Signed.long

val flxy_eval_powers_pre :
  _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong

val fly_to_flxy : _ t -> Signed.long -> _ t
val kronecker_to_flxqx : _ t -> _ t -> pari_ulong -> _ t
val kronecker_to_flxqx_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val rgx_to_flxqx : _ t -> _ t -> pari_ulong -> _ t

val get_flxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  pari_ulong ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val pol1_flxx : Signed.long -> Signed.long -> _ t
val polx_flxx : Signed.long -> Signed.long -> _ t
val random_flxqx : Signed.long -> Signed.long -> _ t -> pari_ulong -> _ t
val zlxx_translate1 : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val zxx_to_kronecker : _ t -> _ t -> _ t
val flxq_ellcard : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxq_ellgens : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t

val flxq_ellgroup :
  _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t Ctypes_static.ptr -> _ t

val flxq_elltwist :
  _ t ->
  _ t ->
  _ t ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  unit

val flxq_ellj : _ t -> _ t -> _ t -> pari_ulong -> _ t

val flxq_ellj_to_a4a6 :
  _ t ->
  _ t ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  unit

val flxqe_add : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqe_changepoint : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqe_changepointinv : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqe_dbl : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqe_log : _ t -> _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqe_mul : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqe_neg : _ t -> _ t -> pari_ulong -> _ t
val flxqe_order : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqe_sub : _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqe_tatepairing : _ t -> _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqe_weilpairing : _ t -> _ t -> _ t -> _ t -> _ t -> pari_ulong -> _ t

val flxqe_weilpairing_pre :
  _ t -> _ t -> _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val zxx_to_flxx : _ t -> pari_ulong -> Signed.long -> _ t
val zxxt_to_flxxt : _ t -> pari_ulong -> Signed.long -> _ t
val zxxv_to_flxxv : _ t -> pari_ulong -> Signed.long -> _ t

val get_flxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  _ t ->
  pari_ulong ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_flxqe : _ t -> _ t -> pari_ulong -> _ t
val random_flxqe : _ t -> _ t -> _ t -> pari_ulong -> _ t
val polisclass : _ t -> Signed.long
val fl_elltrace : pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fl_elltrace_cm :
  Signed.long -> pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fp_ellcard : _ t -> _ t -> _ t -> _ t
val fp_elldivpol : _ t -> _ t -> Signed.long -> _ t -> _ t
val fp_ellgens : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fp_ellgroup : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fp_ellj : _ t -> _ t -> _ t -> _ t

val fp_ellj_to_a4a6 :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit

val fp_elljissupersingular : _ t -> _ t -> int

val fp_elltwist :
  _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit

val fp_ffellcard : _ t -> _ t -> _ t -> Signed.long -> _ t -> _ t
val fpe_add : _ t -> _ t -> _ t -> _ t -> _ t
val fpe_changepoint : _ t -> _ t -> _ t -> _ t
val fpe_changepointinv : _ t -> _ t -> _ t -> _ t
val fpe_dbl : _ t -> _ t -> _ t -> _ t
val fpe_log : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpe_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fpe_neg : _ t -> _ t -> _ t
val fpe_order : _ t -> _ t -> _ t -> _ t -> _ t
val fpe_sub : _ t -> _ t -> _ t -> _ t -> _ t
val fpe_to_fpj : _ t -> _ t
val fpe_to_mod : _ t -> _ t -> _ t
val fpe_tatepairing : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpe_weilpairing : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpj_add : _ t -> _ t -> _ t -> _ t -> _ t
val fpj_dbl : _ t -> _ t -> _ t -> _ t
val fpj_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fpj_neg : _ t -> _ t -> _ t
val fpj_to_fpe : _ t -> _ t -> _ t
val fpxq_ellcard : _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_ellcard_supersingular : _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_elldivpol : _ t -> _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxq_ellgens : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t

val fpxq_ellgroup :
  _ t -> _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t

val fpxq_ellj : _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_elljissupersingular : _ t -> _ t -> _ t -> int

val fpxq_elltwist :
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  unit

val fpxqe_add : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqe_changepoint : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqe_changepointinv : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqe_dbl : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqe_log : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqe_mul : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqe_neg : _ t -> _ t -> _ t -> _ t
val fpxqe_order : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqe_sub : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqe_tatepairing : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqe_weilpairing : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fq_elljissupersingular : _ t -> _ t -> _ t -> int
val fq_ellcard_supersingular : _ t -> _ t -> _ t -> _ t -> _ t
val rge_to_fpe : _ t -> _ t -> _ t
val rge_to_fpxqe : _ t -> _ t -> _ t -> _ t

val get_fpe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  _ t ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  bb_group Ctypes.structure Ctypes_static.ptr

val ellsupersingularj_fpxq : _ t -> _ t -> _ t
val elltrace_extension : _ t -> Signed.long -> _ t -> _ t
val random_fpe : _ t -> _ t -> _ t -> _ t
val random_fpxqe : _ t -> _ t -> _ t -> _ t -> _ t
val fp_issquare : _ t -> _ t -> int
val fp_fpx_sub : _ t -> _ t -> _ t -> _ t
val fp_fpxq_log : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpv_fpm_polint : _ t -> _ t -> _ t -> Signed.long -> _ t
val fpv_inv : _ t -> _ t -> _ t
val fpv_invvandermonde : _ t -> _ t -> _ t -> _ t
val fpv_polint : _ t -> _ t -> _ t -> Signed.long -> _ t
val fpv_roots_to_pol : _ t -> _ t -> Signed.long -> _ t
val fpx_fp_add : _ t -> _ t -> _ t -> _ t
val fpx_fp_add_shallow : _ t -> _ t -> _ t -> _ t
val fpx_fp_div : _ t -> _ t -> _ t -> _ t
val fpx_fp_mul : _ t -> _ t -> _ t -> _ t
val fpx_fp_mul_to_monic : _ t -> _ t -> _ t -> _ t
val fpx_fp_mulspec : _ t -> _ t -> _ t -> Signed.long -> _ t
val fpx_fp_sub : _ t -> _ t -> _ t -> _ t
val fpx_fp_sub_shallow : _ t -> _ t -> _ t -> _ t
val fpx_fpv_multieval : _ t -> _ t -> _ t -> _ t
val fpx_fpxq_eval : _ t -> _ t -> _ t -> _ t -> _ t
val fpx_fpxqv_eval : _ t -> _ t -> _ t -> _ t -> _ t
val fpx_fpxv_multirem : _ t -> _ t -> _ t -> _ t
val fpx_frobenius : _ t -> _ t -> _ t
val fpx_laplace : _ t -> _ t -> _ t
val fpx_newton : _ t -> Signed.long -> _ t -> _ t
val fpx_add : _ t -> _ t -> _ t -> _ t
val fpx_center : _ t -> _ t -> _ t -> _ t
val fpx_center_i : _ t -> _ t -> _ t -> _ t
val fpx_chinese_coprime : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpx_composedprod : _ t -> _ t -> _ t -> _ t
val fpx_composedsum : _ t -> _ t -> _ t -> _ t
val fpx_convol : _ t -> _ t -> _ t -> _ t
val fpx_deriv : _ t -> _ t -> _ t
val fpx_digits : _ t -> _ t -> _ t -> _ t
val fpx_disc : _ t -> _ t -> _ t
val fpx_div_by_x_x : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fpx_divrem : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fpx_divu : _ t -> pari_ulong -> _ t -> _ t
val fpx_dotproduct : _ t -> _ t -> _ t -> _ t
val fpx_eval : _ t -> _ t -> _ t -> _ t

val fpx_extgcd :
  _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val fpx_extresultant :
  _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val fpx_fromnewton : _ t -> _ t -> _ t
val fpx_gcd : _ t -> _ t -> _ t -> _ t
val fpx_gcd_check : _ t -> _ t -> _ t -> _ t
val fpx_get_red : _ t -> _ t -> _ t
val fpx_halve : _ t -> _ t -> _ t
val fpx_halfgcd : _ t -> _ t -> _ t -> _ t

val fpx_halfgcd_all :
  _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val fpx_integ : _ t -> _ t -> _ t
val fpx_invbarrett : _ t -> _ t -> _ t
val fpx_invlaplace : _ t -> _ t -> _ t
val fpx_is_squarefree : _ t -> _ t -> int
val fpx_matfrobenius : _ t -> _ t -> _ t
val fpx_mul : _ t -> _ t -> _ t -> _ t
val fpx_mulspec : _ t -> _ t -> _ t -> Signed.long -> Signed.long -> _ t
val fpx_mulu : _ t -> pari_ulong -> _ t -> _ t
val fpx_neg : _ t -> _ t -> _ t
val fpx_normalize : _ t -> _ t -> _ t
val fpx_powu : _ t -> pari_ulong -> _ t -> _ t
val fpx_red : _ t -> _ t -> _ t
val fpx_rem : _ t -> _ t -> _ t -> _ t
val fpx_rescale : _ t -> _ t -> _ t -> _ t
val fpx_resultant : _ t -> _ t -> _ t -> _ t
val fpx_sqr : _ t -> _ t -> _ t
val fpx_sub : _ t -> _ t -> _ t -> _ t
val fpx_valrem : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val fpxc_fpxq_eval : _ t -> _ t -> _ t -> _ t -> _ t
val fpxc_fpxqv_eval : _ t -> _ t -> _ t -> _ t -> _ t
val fpxm_fpxqv_eval : _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_autpow : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fpxq_autpowers : _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxq_autsum : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fpxq_auttrace : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fpxq_charpoly : _ t -> _ t -> _ t -> _ t
val fpxq_conjvec : _ t -> _ t -> _ t -> _ t
val fpxq_div : _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_inv : _ t -> _ t -> _ t -> _ t
val fpxq_invsafe : _ t -> _ t -> _ t -> _ t
val fpxq_issquare : _ t -> _ t -> _ t -> int
val fpxq_log : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_matrix_pow : _ t -> Signed.long -> Signed.long -> _ t -> _ t -> _ t
val fpxq_minpoly : _ t -> _ t -> _ t -> _ t
val fpxq_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_norm : _ t -> _ t -> _ t -> _ t
val fpxq_order : _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_pow : _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_powu : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fpxq_powers : _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxq_red : _ t -> _ t -> _ t -> _ t
val fpxq_sqr : _ t -> _ t -> _ t -> _ t
val fpxq_sqrt : _ t -> _ t -> _ t -> _ t
val fpxq_sqrtn : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fpxq_trace : _ t -> _ t -> _ t -> _ t
val fpxqc_to_mod : _ t -> _ t -> _ t -> _ t
val fpxqm_autsum : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fpxt_red : _ t -> _ t -> _ t
val fpxv_fpx_fromdigits : _ t -> _ t -> _ t -> _ t
val fpxv_chinese : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fpxv_composedsum : _ t -> _ t -> _ t
val fpxv_factorback : _ t -> _ t -> _ t -> Signed.long -> _ t
val fpxv_prod : _ t -> _ t -> _ t
val fpxv_red : _ t -> _ t -> _ t
val fpxn_div : _ t -> _ t -> Signed.long -> _ t -> _ t
val fpxn_exp : _ t -> Signed.long -> _ t -> _ t
val fpxn_expint : _ t -> Signed.long -> _ t -> _ t
val fpxn_inv : _ t -> Signed.long -> _ t -> _ t
val fpxn_mul : _ t -> _ t -> Signed.long -> _ t -> _ t
val fpxn_sqr : _ t -> Signed.long -> _ t -> _ t
val fq_issquare : _ t -> _ t -> _ t -> int
val fq_ispower : _ t -> _ t -> _ t -> _ t -> Signed.long
val fq_log : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fqc_to_mod : _ t -> _ t -> _ t -> _ t
val fqm_to_mod : _ t -> _ t -> _ t -> _ t
val fqv_inv : _ t -> _ t -> _ t -> _ t
val z_to_fpx : _ t -> _ t -> Signed.long -> _ t
val gener_fpxq : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val gener_fpxq_local : _ t -> _ t -> _ t -> _ t

val get_fpxq_star :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpx : Signed.long -> Signed.long -> _ t -> _ t
val f2x_ddf : _ t -> _ t
val f2x_factor : _ t -> _ t
val f2x_factor_squarefree : _ t -> _ t
val f2x_is_irred : _ t -> int
val flx_ddf : _ t -> pari_ulong -> _ t
val flx_ddf_pre : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_is_irred : _ t -> pari_ulong -> int
val flx_is_totally_split : _ t -> pari_ulong -> int

val flx_ispower :
  _ t -> pari_ulong -> pari_ulong -> _ t Ctypes_static.ptr -> Signed.long

val flx_degfact : _ t -> pari_ulong -> _ t
val flx_factor : _ t -> pari_ulong -> _ t
val flx_factor_squarefree : _ t -> pari_ulong -> _ t
val flx_factor_squarefree_pre : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_nbfact : _ t -> pari_ulong -> Signed.long
val flx_nbfact_pre : _ t -> pari_ulong -> pari_ulong -> Signed.long
val flx_nbfact_frobenius : _ t -> _ t -> pari_ulong -> Signed.long

val flx_nbfact_frobenius_pre :
  _ t -> _ t -> pari_ulong -> pari_ulong -> Signed.long

val flx_nbfact_by_degree :
  _ t -> Signed.long Ctypes_static.ptr -> pari_ulong -> _ t

val flx_nbroots : _ t -> pari_ulong -> Signed.long
val flx_oneroot : _ t -> pari_ulong -> pari_ulong
val flx_oneroot_split : _ t -> pari_ulong -> pari_ulong
val flx_oneroot_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong
val flx_oneroot_split_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong
val flx_roots : _ t -> pari_ulong -> _ t
val flx_roots_pre : _ t -> pari_ulong -> pari_ulong -> _ t
val flx_rootsff : _ t -> _ t -> pari_ulong -> _ t
val fpx_ddf : _ t -> _ t -> _ t
val fpx_ddf_degree : _ t -> _ t -> _ t -> Signed.long
val fpx_degfact : _ t -> _ t -> _ t
val fpx_factor : _ t -> _ t -> _ t
val fpx_factor_squarefree : _ t -> _ t -> _ t
val fpx_is_irred : _ t -> _ t -> int
val fpx_is_totally_split : _ t -> _ t -> int

val fpx_ispower :
  _ t -> pari_ulong -> _ t -> _ t Ctypes_static.ptr -> Signed.long

val fpx_nbfact : _ t -> _ t -> Signed.long
val fpx_nbfact_frobenius : _ t -> _ t -> _ t -> Signed.long
val fpx_nbroots : _ t -> _ t -> Signed.long
val fpx_oneroot : _ t -> _ t -> _ t
val fpx_oneroot_split : _ t -> _ t -> _ t
val fpx_roots : _ t -> _ t -> _ t
val fpx_roots_mult : _ t -> Signed.long -> _ t -> _ t
val fpx_rootsff : _ t -> _ t -> _ t -> _ t
val fpx_split_part : _ t -> _ t -> _ t
val f2xqx_ddf : _ t -> _ t -> _ t
val f2xqx_degfact : _ t -> _ t -> _ t
val f2xqx_factor : _ t -> _ t -> _ t
val f2xqx_factor_squarefree : _ t -> _ t -> _ t
val f2xqx_roots : _ t -> _ t -> _ t
val flx_factorff_irred : _ t -> _ t -> pari_ulong -> _ t

val flx_ffintersect :
  _ t ->
  _ t ->
  Signed.long ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t ->
  _ t ->
  unit

val flx_ffisom : _ t -> _ t -> pari_ulong -> _ t
val flxq_ffisom_inv : _ t -> _ t -> pari_ulong -> _ t
val flxqx_frobenius : _ t -> _ t -> pari_ulong -> _ t
val flxqx_frobenius_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flxqx_ddf : _ t -> _ t -> pari_ulong -> _ t
val flxqx_ddf_degree : _ t -> _ t -> _ t -> pari_ulong -> Signed.long
val flxqx_degfact : _ t -> _ t -> pari_ulong -> _ t
val flxqx_factor : _ t -> _ t -> pari_ulong -> _ t
val flxqx_factor_squarefree : _ t -> _ t -> pari_ulong -> _ t
val flxqx_factor_squarefree_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flxqx_ispower :
  _ t -> pari_ulong -> _ t -> pari_ulong -> _ t Ctypes_static.ptr -> Signed.long

val flxqx_is_squarefree : _ t -> _ t -> pari_ulong -> Signed.long
val flxqx_nbfact : _ t -> _ t -> pari_ulong -> Signed.long
val flxqx_nbfact_frobenius : _ t -> _ t -> _ t -> pari_ulong -> Signed.long

val flxqx_nbfact_by_degree :
  _ t -> Signed.long Ctypes_static.ptr -> _ t -> pari_ulong -> _ t

val flxqx_nbroots : _ t -> _ t -> pari_ulong -> Signed.long
val flxqx_roots : _ t -> _ t -> pari_ulong -> _ t
val flxqxq_halffrobenius : _ t -> _ t -> _ t -> pari_ulong -> _ t
val fpx_factorff : _ t -> _ t -> _ t -> _ t
val fpx_factorff_irred : _ t -> _ t -> _ t -> _ t

val fpx_ffintersect :
  _ t ->
  _ t ->
  Signed.long ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t ->
  _ t ->
  unit

val fpx_ffisom : _ t -> _ t -> _ t -> _ t
val fpxq_ffisom_inv : _ t -> _ t -> _ t -> _ t
val fpxqx_frobenius : _ t -> _ t -> _ t -> _ t
val fpxqx_ddf : _ t -> _ t -> _ t -> _ t
val fpxqx_ddf_degree : _ t -> _ t -> _ t -> _ t -> Signed.long
val fpxqx_degfact : _ t -> _ t -> _ t -> _ t
val fpxqx_factor : _ t -> _ t -> _ t -> _ t
val fpxqx_factor_squarefree : _ t -> _ t -> _ t -> _ t

val fpxqx_ispower :
  _ t -> pari_ulong -> _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long

val fpxqx_nbfact : _ t -> _ t -> _ t -> Signed.long
val fpxqx_nbfact_frobenius : _ t -> _ t -> _ t -> _ t -> Signed.long
val fpxqx_nbroots : _ t -> _ t -> _ t -> Signed.long
val fpxqx_roots : _ t -> _ t -> _ t -> _ t
val fpxqx_split_part : _ t -> _ t -> _ t -> _ t
val fpxqxq_halffrobenius : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_is_squarefree : _ t -> _ t -> _ t -> Signed.long

val fqx_ispower :
  _ t -> pari_ulong -> _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long

val fqx_nbfact : _ t -> _ t -> _ t -> Signed.long
val fqx_nbroots : _ t -> _ t -> _ t -> Signed.long
val factorff : _ t -> _ t -> _ t -> _ t
val factormod0 : _ t -> _ t -> Signed.long -> _ t
val factormodddf : _ t -> _ t -> _ t
val factormodsqf : _ t -> _ t -> _ t

val ff_parse_tp :
  _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> Signed.long -> int

val polrootsff : _ t -> _ t -> _ t -> _ t
val polrootsmod : _ t -> _ t -> _ t
val rootmod0 : _ t -> _ t -> Signed.long -> _ t
val fpxqx_fpxq_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqx_fpxqxqv_eval : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqx_fpxqxq_eval : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqx_digits : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqx_disc : _ t -> _ t -> _ t -> _ t
val fpxqx_div_by_x_x : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fpxqx_divrem : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fpxqx_dotproduct : _ t -> _ t -> _ t -> _ t -> _ t

val fpxqx_extgcd :
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val fpxqx_gcd : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqx_get_red : _ t -> _ t -> _ t -> _ t
val fpxqx_halfgcd : _ t -> _ t -> _ t -> _ t -> _ t

val fpxqx_halfgcd_all :
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val fpxqx_invbarrett : _ t -> _ t -> _ t -> _ t
val fpxqx_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqx_powu : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fpxqx_red : _ t -> _ t -> _ t -> _ t
val fpxqx_rem : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqx_resultant : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqx_sqr : _ t -> _ t -> _ t -> _ t
val fpxqx_to_mod : _ t -> _ t -> _ t -> _ t
val fpxqxq_div : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqxq_inv : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqxq_invsafe : _ t -> _ t -> _ t -> _ t -> _ t

val fpxqxq_matrix_pow :
  _ t -> Signed.long -> Signed.long -> _ t -> _ t -> _ t -> _ t

val fpxqxq_minpoly : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqxq_mul : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqxq_pow : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpxqxq_powers : _ t -> Signed.long -> _ t -> _ t -> _ t -> _ t
val fpxqxq_sqr : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqxq_autpow : _ t -> Signed.long -> _ t -> _ t -> _ t -> _ t
val fpxqxq_autsum : _ t -> Signed.long -> _ t -> _ t -> _ t -> _ t
val fpxqxq_auttrace : _ t -> Signed.long -> _ t -> _ t -> _ t -> _ t
val fpxqxt_red : _ t -> _ t -> _ t -> _ t
val fpxqxv_fpxqx_fromdigits : _ t -> _ t -> _ t -> _ t -> _ t
val fpxqxv_prod : _ t -> _ t -> _ t -> _ t
val fpxqxv_red : _ t -> _ t -> _ t -> _ t
val fpxqxn_div : _ t -> _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxqxn_exp : _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxqxn_expint : _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxqxn_inv : _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxqxn_mul : _ t -> _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxqxn_sqr : _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxx_fp_mul : _ t -> _ t -> _ t -> _ t
val fpxx_fpx_mul : _ t -> _ t -> _ t -> _ t
val fpxx_add : _ t -> _ t -> _ t -> _ t
val fpxx_deriv : _ t -> _ t -> _ t
val fpxx_halve : _ t -> _ t -> _ t
val fpxx_integ : _ t -> _ t -> _ t
val fpxx_mulu : _ t -> pari_ulong -> _ t -> _ t
val fpxx_neg : _ t -> _ t -> _ t
val fpxx_red : _ t -> _ t -> _ t
val fpxx_sub : _ t -> _ t -> _ t -> _ t
val fpxy_fpxq_evalx : _ t -> _ t -> _ t -> _ t -> _ t
val fpxy_fpxqv_evalx : _ t -> _ t -> _ t -> _ t -> _ t
val fpxy_eval : _ t -> _ t -> _ t -> _ t -> _ t
val fpxy_evalx : _ t -> _ t -> _ t -> _ t
val fpxy_evaly : _ t -> _ t -> _ t -> Signed.long -> _ t
val fpxyqq_pow : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fqxc_to_mod : _ t -> _ t -> _ t -> _ t
val fqxm_to_mod : _ t -> _ t -> _ t -> _ t
val kronecker_to_fpxqx : _ t -> _ t -> _ t -> _ t

val get_fpxqx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  _ t ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpxqx : Signed.long -> Signed.long -> _ t -> _ t -> _ t
val flc_flv_mul : _ t -> _ t -> pari_ulong -> _ t
val flc_to_mod : _ t -> pari_ulong -> _ t
val flm_fl_add : _ t -> pari_ulong -> pari_ulong -> _ t
val flm_fl_mul : _ t -> pari_ulong -> pari_ulong -> _ t
val flm_fl_mul_inplace : _ t -> pari_ulong -> pari_ulong -> unit
val flm_fl_mul_pre : _ t -> pari_ulong -> pari_ulong -> pari_ulong -> _ t
val flm_fl_sub : _ t -> pari_ulong -> pari_ulong -> _ t
val flm_flc_mul : _ t -> _ t -> pari_ulong -> _ t
val flm_flc_mul_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t

val flm_flc_mul_pre_flx :
  _ t -> _ t -> pari_ulong -> pari_ulong -> Signed.long -> _ t

val flm_add : _ t -> _ t -> pari_ulong -> _ t
val flm_center : _ t -> pari_ulong -> pari_ulong -> _ t
val flm_mul : _ t -> _ t -> pari_ulong -> _ t
val flm_mul_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val flm_neg : _ t -> pari_ulong -> _ t
val flm_powers : _ t -> pari_ulong -> pari_ulong -> _ t
val flm_powu : _ t -> pari_ulong -> pari_ulong -> _ t
val flm_sub : _ t -> _ t -> pari_ulong -> _ t
val flm_to_mod : _ t -> pari_ulong -> _ t
val flm_transpose : _ t -> _ t
val flv_fl_div : _ t -> pari_ulong -> pari_ulong -> _ t
val flv_fl_div_inplace : _ t -> pari_ulong -> pari_ulong -> unit
val flv_fl_mul : _ t -> pari_ulong -> pari_ulong -> _ t
val flv_fl_mul_inplace : _ t -> pari_ulong -> pari_ulong -> unit

val flv_fl_mul_part_inplace :
  _ t -> pari_ulong -> pari_ulong -> Signed.long -> unit

val flv_add : _ t -> _ t -> pari_ulong -> _ t
val flv_add_inplace : _ t -> _ t -> pari_ulong -> unit
val flv_center : _ t -> pari_ulong -> pari_ulong -> _ t
val flv_dotproduct : _ t -> _ t -> pari_ulong -> pari_ulong
val flv_dotproduct_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong
val flv_neg : _ t -> pari_ulong -> _ t
val flv_neg_inplace : _ t -> pari_ulong -> unit
val flv_sub : _ t -> _ t -> pari_ulong -> _ t
val flv_sub_inplace : _ t -> _ t -> pari_ulong -> unit
val flv_sum : _ t -> pari_ulong -> pari_ulong
val flx_dotproduct : _ t -> _ t -> pari_ulong -> pari_ulong
val flx_dotproduct_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> pari_ulong
val fp_to_mod : _ t -> _ t -> _ t
val fpc_fpv_mul : _ t -> _ t -> _ t -> _ t
val fpc_fp_mul : _ t -> _ t -> _ t -> _ t
val fpc_center : _ t -> _ t -> _ t -> _ t
val fpc_center_inplace : _ t -> _ t -> _ t -> unit
val fpc_red : _ t -> _ t -> _ t
val fpc_to_mod : _ t -> _ t -> _ t
val fpm_add : _ t -> _ t -> _ t -> _ t
val fpm_fp_mul : _ t -> _ t -> _ t -> _ t
val fpm_fpc_mul : _ t -> _ t -> _ t -> _ t
val fpm_fpc_mul_fpx : _ t -> _ t -> _ t -> Signed.long -> _ t
val fpm_center : _ t -> _ t -> _ t -> _ t
val fpm_center_inplace : _ t -> _ t -> _ t -> unit
val fpm_mul : _ t -> _ t -> _ t -> _ t
val fpm_powu : _ t -> pari_ulong -> _ t -> _ t
val fpm_red : _ t -> _ t -> _ t
val fpm_sub : _ t -> _ t -> _ t -> _ t
val fpm_to_mod : _ t -> _ t -> _ t
val fpms_fpc_mul : _ t -> _ t -> _ t -> _ t
val fpms_fpcs_solve : _ t -> _ t -> Signed.long -> _ t -> _ t
val fpms_fpcs_solve_safe : _ t -> _ t -> Signed.long -> _ t -> _ t
val fpms_leftkernel_elt : _ t -> Signed.long -> _ t -> _ t
val fpc_add : _ t -> _ t -> _ t -> _ t
val fpc_sub : _ t -> _ t -> _ t -> _ t
val fpv_fpms_mul : _ t -> _ t -> _ t -> _ t
val fpv_add : _ t -> _ t -> _ t -> _ t
val fpv_sub : _ t -> _ t -> _ t -> _ t
val fpv_dotproduct : _ t -> _ t -> _ t -> _ t
val fpv_dotsquare : _ t -> _ t -> _ t
val fpv_red : _ t -> _ t -> _ t
val fpv_to_mod : _ t -> _ t -> _ t
val fpvv_to_mod : _ t -> _ t -> _ t
val fpx_to_mod : _ t -> _ t -> _ t
val fpxc_to_mod : _ t -> _ t -> _ t
val fpxm_to_mod : _ t -> _ t -> _ t
val zabm_ker : _ t -> _ t -> Signed.long -> _ t
val zabm_indexrank : _ t -> _ t -> Signed.long -> _ t
val zabm_inv : _ t -> _ t -> Signed.long -> _ t Ctypes_static.ptr -> _ t
val zabm_inv_ratlift : _ t -> _ t -> Signed.long -> _ t Ctypes_static.ptr -> _ t

val zabm_pseudoinv :
  _ t ->
  _ t ->
  Signed.long ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val zv_zms_mul : _ t -> _ t -> _ t
val zpms_zpcs_solve : _ t -> _ t -> Signed.long -> _ t -> Signed.long -> _ t

val gen_fpm_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  _ t

val gen_zpm_dixon_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val gen_matid :
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  _ t

val matid_flm : Signed.long -> _ t
val matid_f2xqm : Signed.long -> _ t -> _ t
val matid_flxqm : Signed.long -> _ t -> pari_ulong -> _ t
val random_flv : Signed.long -> pari_ulong -> _ t
val random_fpc : Signed.long -> _ t -> _ t
val random_fpv : Signed.long -> _ t -> _ t
val scalar_flm : Signed.long -> Signed.long -> _ t
val zcs_to_zc : _ t -> Signed.long -> _ t
val zms_to_zm : _ t -> Signed.long -> _ t
val zms_zc_mul : _ t -> _ t -> _ t
val zmv_to_flmv : _ t -> pari_ulong -> _ t
val flx_teichmuller : _ t -> pari_ulong -> Signed.long -> _ t
val z2_sqrt : _ t -> Signed.long -> _ t
val zp_div : _ t -> _ t -> _ t -> Signed.long -> _ t
val zp_exp : _ t -> _ t -> pari_ulong -> _ t
val zp_inv : _ t -> _ t -> Signed.long -> _ t
val zp_invlift : _ t -> _ t -> _ t -> Signed.long -> _ t
val zp_log : _ t -> _ t -> pari_ulong -> _ t
val zp_sqrt : _ t -> _ t -> Signed.long -> _ t
val zp_sqrtlift : _ t -> _ t -> _ t -> Signed.long -> _ t
val zp_sqrtnlift : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zpm_invlift : _ t -> _ t -> _ t -> Signed.long -> _ t
val zpx_frobenius : _ t -> _ t -> Signed.long -> _ t
val zpx_zpxq_liftroot : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t

val zpx_zpxq_liftroot_ea :
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val zpx_liftfact : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zpx_liftroot : _ t -> _ t -> _ t -> Signed.long -> _ t
val zpx_liftroots : _ t -> _ t -> _ t -> Signed.long -> _ t
val zpx_roots : _ t -> _ t -> Signed.long -> _ t
val zpxq_div : _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zpxq_inv : _ t -> _ t -> _ t -> Signed.long -> _ t
val zpxq_invlift : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zpxq_log : _ t -> _ t -> _ t -> Signed.long -> _ t
val zpxq_sqrt : _ t -> _ t -> _ t -> Signed.long -> _ t
val zpxq_sqrtnlift : _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zpxqm_prodfrobenius : _ t -> _ t -> _ t -> Signed.long -> _ t
val zpxqx_digits : _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t

val zpxqx_divrem :
  _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t Ctypes_static.ptr -> _ t

val zpxqx_liftfact : _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zpxqx_liftroot : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t

val zpxqx_liftroot_vald :
  _ t -> _ t -> Signed.long -> _ t -> _ t -> Signed.long -> _ t

val zpxqx_liftroots : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zpxqx_roots : _ t -> _ t -> _ t -> Signed.long -> _ t

val zpxqx_zpxqxq_liftroot :
  _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t

val zq_sqrtnlift : _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zqx_zqxq_liftroot : _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zqx_liftfact : _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zqx_liftroot : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val zqx_roots : _ t -> _ t -> _ t -> Signed.long -> _ t

val gen_zpm_dixon :
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t -> _ t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_zpm_newton :
  _ t ->
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t

val gen_zpx_dixon :
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t -> _ t)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_zpx_newton :
  _ t ->
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t

val polteichmuller : _ t -> pari_ulong -> Signed.long -> _ t
val polhensellift : _ t -> _ t -> _ t -> Signed.long -> _ t
val quadratic_prec_mask : Signed.long -> pari_ulong
val qx_factor : _ t -> _ t
val zx_factor : _ t -> _ t
val zx_is_irred : _ t -> Signed.long
val zx_squff : _ t -> _ t Ctypes_static.ptr -> _ t
val polcyclofactors : _ t -> _ t
val poliscyclo : _ t -> Signed.long
val poliscycloprod : _ t -> Signed.long
val rg_rgc_sub : _ t -> _ t -> _ t
val rgc_rg_add : _ t -> _ t -> _ t
val rgc_rg_div : _ t -> _ t -> _ t
val rgc_rg_mul : _ t -> _ t -> _ t
val rgc_rg_sub : _ t -> _ t -> _ t
val rgc_rgm_mul : _ t -> _ t -> _ t
val rgc_rgv_mul : _ t -> _ t -> _ t
val rgc_add : _ t -> _ t -> _ t
val rgc_is_ei : _ t -> Signed.long
val rgc_neg : _ t -> _ t
val rgc_sub : _ t -> _ t -> _ t
val rgm_rg_add : _ t -> _ t -> _ t
val rgm_rg_add_shallow : _ t -> _ t -> _ t
val rgm_rg_div : _ t -> _ t -> _ t
val rgm_rg_mul : _ t -> _ t -> _ t
val rgm_rg_sub : _ t -> _ t -> _ t
val rgm_rg_sub_shallow : _ t -> _ t -> _ t
val rgm_rgc_mul : _ t -> _ t -> _ t
val rgm_rgv_mul : _ t -> _ t -> _ t
val rgm_add : _ t -> _ t -> _ t
val rgm_det_triangular : _ t -> _ t
val rgm_is_qm : _ t -> int
val rgm_is_zm : _ t -> int
val rgm_isdiagonal : _ t -> int
val rgm_isidentity : _ t -> int
val rgm_isscalar : _ t -> _ t -> int
val rgm_mul : _ t -> _ t -> _ t
val rgm_multosym : _ t -> _ t -> _ t
val rgm_neg : _ t -> _ t
val rgm_powers : _ t -> Signed.long -> _ t
val rgm_sqr : _ t -> _ t
val rgm_sub : _ t -> _ t -> _ t
val rgm_sumcol : _ t -> _ t
val rgm_transmul : _ t -> _ t -> _ t
val rgm_transmultosym : _ t -> _ t -> _ t
val rgmrow_zc_mul : _ t -> _ t -> Signed.long -> _ t
val rgm_zc_mul : _ t -> _ t -> _ t
val rgm_zm_mul : _ t -> _ t -> _ t
val rgmrow_rgc_mul : _ t -> _ t -> Signed.long -> _ t
val rgv_rgm_mul : _ t -> _ t -> _ t
val rgv_rgc_mul : _ t -> _ t -> _ t
val rgv_rg_mul : _ t -> _ t -> _ t
val rgv_add : _ t -> _ t -> _ t
val rgv_dotproduct : _ t -> _ t -> _ t
val rgv_dotsquare : _ t -> _ t
val rgv_is_zmv : _ t -> int
val rgv_kill0 : _ t -> _ t
val rgv_neg : _ t -> _ t
val rgv_prod : _ t -> _ t
val rgv_sub : _ t -> _ t -> _ t
val rgv_sum : _ t -> _ t
val rgv_sumpart : _ t -> Signed.long -> _ t
val rgv_sumpart2 : _ t -> Signed.long -> Signed.long -> _ t
val rgv_zc_mul : _ t -> _ t -> _ t
val rgv_zm_mul : _ t -> _ t -> _ t
val rgx_rgm_eval : _ t -> _ t -> _ t
val rgx_rgmv_eval : _ t -> _ t -> _ t
val isdiagonal : _ t -> int
val matid : Signed.long -> _ t
val scalarcol : _ t -> Signed.long -> _ t
val scalarcol_shallow : _ t -> Signed.long -> _ t
val scalarmat : _ t -> Signed.long -> _ t
val scalarmat_shallow : _ t -> Signed.long -> _ t
val scalarmat_s : Signed.long -> Signed.long -> _ t
val kronecker_to_mod : _ t -> _ t -> _ t
val qx_zxqv_eval : _ t -> _ t -> _ t -> _ t
val qxq_charpoly : _ t -> _ t -> Signed.long -> _ t
val qxq_powers : _ t -> Signed.long -> _ t -> _ t
val qxq_to_mod_shallow : _ t -> _ t -> _ t
val qxqc_to_mod_shallow : _ t -> _ t -> _ t
val qxqm_to_mod_shallow : _ t -> _ t -> _ t
val qxqv_to_mod : _ t -> _ t -> _ t
val qxqx_homogenous_evalpow : _ t -> _ t -> _ t -> _ t -> _ t
val qxqx_to_mod_shallow : _ t -> _ t -> _ t
val qxqxv_to_mod : _ t -> _ t -> _ t
val qxv_qxq_eval : _ t -> _ t -> _ t -> _ t
val qxy_qxq_evalx : _ t -> _ t -> _ t -> _ t
val rg_rgx_sub : _ t -> _ t -> _ t
val rg_get_0 : _ t -> _ t
val rg_get_1 : _ t -> _ t
val rg_to_rgc : _ t -> Signed.long -> _ t
val rgm_to_rgxv : _ t -> Signed.long -> _ t
val rgm_to_rgxv_reverse : _ t -> Signed.long -> _ t
val rgm_to_rgxx : _ t -> Signed.long -> Signed.long -> _ t
val rgv_to_rgx : _ t -> Signed.long -> _ t
val rgv_to_rgm : _ t -> Signed.long -> _ t
val rgv_to_rgx_reverse : _ t -> Signed.long -> _ t
val rgx_rgxq_eval : _ t -> _ t -> _ t -> _ t
val rgx_rgxqv_eval : _ t -> _ t -> _ t -> _ t
val rgx_rgxn_eval : _ t -> _ t -> Signed.long -> _ t
val rgx_rgxnv_eval : _ t -> _ t -> Signed.long -> _ t
val rgx_rg_add : _ t -> _ t -> _ t
val rgx_rg_add_shallow : _ t -> _ t -> _ t
val rgx_rg_div : _ t -> _ t -> _ t
val rgx_rg_divexact : _ t -> _ t -> _ t
val rgx_rg_eval_bk : _ t -> _ t -> _ t
val rgx_rg_mul : _ t -> _ t -> _ t
val rgx_rg_sub : _ t -> _ t -> _ t
val rgx_rgv_eval : _ t -> _ t -> _ t
val rgx_add : _ t -> _ t -> _ t
val rgx_addmulxn_shallow : _ t -> _ t -> Signed.long -> _ t
val rgx_addmulxn : _ t -> _ t -> Signed.long -> _ t
val rgx_addspec : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val rgx_addspec_shallow : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val rgx_affine : _ t -> _ t -> _ t -> _ t
val rgx_blocks : _ t -> Signed.long -> Signed.long -> _ t
val rgx_deflate : _ t -> Signed.long -> _ t
val rgx_deriv : _ t -> _ t
val rgx_digits : _ t -> _ t -> _ t
val rgx_div_by_x_x : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val rgx_divrem : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val rgx_divs : _ t -> Signed.long -> _ t
val rgx_equal : _ t -> _ t -> Signed.long
val rgx_even_odd : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit
val rgx_homogenize : _ t -> Signed.long -> _ t
val rgx_homogenous_evalpow : _ t -> _ t -> _ t -> _ t
val rgx_inflate : _ t -> Signed.long -> _ t
val rgx_mul : _ t -> _ t -> _ t
val rgx_mul_i : _ t -> _ t -> _ t
val rgx_mul_normalized : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val rgx_mul2n : _ t -> Signed.long -> _ t
val rgx_mulxn : _ t -> Signed.long -> _ t
val rgx_mulhigh_i : _ t -> _ t -> Signed.long -> _ t
val rgx_muls : _ t -> Signed.long -> _ t
val rgx_mulspec : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val rgx_neg : _ t -> _ t
val rgx_normalize : _ t -> _ t
val rgx_pseudodivrem : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val rgx_pseudorem : _ t -> _ t -> _ t
val rgx_recip : _ t -> _ t
val rgx_recip_i : _ t -> _ t
val rgx_recip_shallow : _ t -> _ t
val rgx_rem : _ t -> _ t -> _ t
val rgx_renormalize_lg : _ t -> Signed.long -> _ t
val rgx_rescale : _ t -> _ t -> _ t
val rgx_rotate_shallow : _ t -> Signed.long -> Signed.long -> _ t
val rgx_shift : _ t -> Signed.long -> _ t
val rgx_shift_shallow : _ t -> Signed.long -> _ t
val rgx_splitting : _ t -> Signed.long -> _ t
val rgx_sqr : _ t -> _ t
val rgx_sqr_i : _ t -> _ t
val rgx_sqrhigh_i : _ t -> Signed.long -> _ t
val rgx_sqrspec : _ t -> Signed.long -> _ t
val rgx_sub : _ t -> _ t -> _ t
val rgx_to_rgc : _ t -> Signed.long -> _ t
val rgx_translate : _ t -> _ t -> _ t
val rgx_unscale : _ t -> _ t -> _ t
val rgxq_matrix_pow : _ t -> Signed.long -> Signed.long -> _ t -> _ t
val rgxq_norm : _ t -> _ t -> _ t
val rgxq_pow : _ t -> _ t -> _ t -> _ t
val rgxq_powers : _ t -> Signed.long -> _ t -> _ t
val rgxq_powu : _ t -> pari_ulong -> _ t -> _ t
val rgxq_trace : _ t -> _ t -> _ t
val rgxqc_red : _ t -> _ t -> _ t
val rgxqm_mul : _ t -> _ t -> _ t -> _ t
val rgxqm_red : _ t -> _ t -> _ t
val rgxqv_rgxq_mul : _ t -> _ t -> _ t -> _ t
val rgxqv_factorback : _ t -> _ t -> _ t -> _ t
val rgxqv_red : _ t -> _ t -> _ t
val rgxqx_rgxq_mul : _ t -> _ t -> _ t -> _ t
val rgxqx_divrem : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val rgxqx_mul : _ t -> _ t -> _ t -> _ t
val rgxqx_powers : _ t -> Signed.long -> _ t -> _ t
val rgxqx_pseudodivrem : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val rgxqx_pseudorem : _ t -> _ t -> _ t -> _ t
val rgxqx_red : _ t -> _ t -> _ t
val rgxqx_sqr : _ t -> _ t -> _ t
val rgxqx_translate : _ t -> _ t -> _ t -> _ t
val rgxv_rgv_eval : _ t -> _ t -> _ t
val rgxv_prod : _ t -> _ t
val rgxv_rescale : _ t -> _ t -> _ t
val rgxv_to_rgm : _ t -> Signed.long -> _ t
val rgxv_unscale : _ t -> _ t -> _ t
val rgxx_to_rgm : _ t -> Signed.long -> _ t
val rgxy_degreex : _ t -> Signed.long
val rgxy_derivx : _ t -> _ t
val rgxy_swap : _ t -> Signed.long -> Signed.long -> _ t
val rgxy_swapspec : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val rgxn_div : _ t -> _ t -> Signed.long -> _ t
val rgxn_div_i : _ t -> _ t -> Signed.long -> _ t
val rgxn_eval : _ t -> _ t -> Signed.long -> _ t
val rgxn_exp : _ t -> Signed.long -> _ t
val rgxn_expint : _ t -> Signed.long -> _ t
val rgxn_inv : _ t -> Signed.long -> _ t
val rgxn_inv_i : _ t -> Signed.long -> _ t
val rgxn_mul : _ t -> _ t -> Signed.long -> _ t
val rgxn_powers : _ t -> Signed.long -> Signed.long -> _ t
val rgxn_recip_shallow : _ t -> Signed.long -> _ t
val rgxn_red_shallow : _ t -> Signed.long -> _ t
val rgxn_reverse : _ t -> Signed.long -> _ t
val rgxn_sqr : _ t -> Signed.long -> _ t
val rgxn_sqrt : _ t -> Signed.long -> _ t
val rgxnv_red_shallow : _ t -> Signed.long -> _ t
val rgxn_powu : _ t -> pari_ulong -> Signed.long -> _ t
val rgxn_powu_i : _ t -> pari_ulong -> Signed.long -> _ t
val zx_translate : _ t -> _ t -> _ t
val zx_unscale2n : _ t -> Signed.long -> _ t
val zx_unscale : _ t -> _ t -> _ t
val zx_unscale_div : _ t -> _ t -> _ t
val zx_unscale_divpow : _ t -> _ t -> Signed.long -> _ t
val zx_z_unscale : _ t -> Signed.long -> _ t
val zxq_powers : _ t -> Signed.long -> _ t -> _ t
val zxq_powu : _ t -> pari_ulong -> _ t -> _ t
val zxqx_dvd : _ t -> _ t -> _ t -> int
val brent_kung_optpow : Signed.long -> Signed.long -> Signed.long -> Signed.long

val gen_bkeval :
  _ t ->
  Signed.long ->
  _ t ->
  int ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t -> _ t)
  Ctypes_static.static_funptr ->
  _ t

val gen_bkeval_powers :
  _ t ->
  Signed.long ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t -> _ t)
  Ctypes_static.static_funptr ->
  _ t

val get_rg_algebra : unit -> bb_algebra Ctypes.structure Ctypes_static.ptr
val rfrac_deflate_order : _ t -> Signed.long
val rfrac_deflate_max : _ t -> Signed.long Ctypes_static.ptr -> _ t
val rfrac_deflate : _ t -> Signed.long -> _ t
val zgc_g_mul_inplace : _ t -> _ t -> unit
val zgcs_add : _ t -> _ t -> _ t
val g_zgc_mul : _ t -> _ t -> _ t
val g_zg_mul : _ t -> _ t -> _ t
val zgc_g_mul : _ t -> _ t -> _ t
val zgc_z_mul : _ t -> _ t -> _ t
val zg_g_mul : _ t -> _ t -> _ t
val zg_z_mul : _ t -> _ t -> _ t
val zg_add : _ t -> _ t -> _ t
val zg_mul : _ t -> _ t -> _ t
val zg_neg : _ t -> _ t
val zg_normalize : _ t -> _ t
val zg_sub : _ t -> _ t -> _ t
val flc_lincomb1_inplace : _ t -> _ t -> pari_ulong -> pari_ulong -> unit
val vecsmall_prod : _ t -> _ t
val qm_qc_mul : _ t -> _ t -> _ t
val qm_det : _ t -> _ t
val qm_ker : _ t -> _ t
val qm_mul : _ t -> _ t -> _ t
val qm_sqr : _ t -> _ t
val rgm_check_zm : _ t -> string -> unit
val rgv_check_zv : _ t -> string -> unit
val z_zc_sub : _ t -> _ t -> _ t
val zv_zc_mul : _ t -> _ t -> _ t
val zc_q_mul : _ t -> _ t -> _ t
val zc_z_add : _ t -> _ t -> _ t
val zc_z_div : _ t -> _ t -> _ t
val zc_z_divexact : _ t -> _ t -> _ t
val zc_z_sub : _ t -> _ t -> _ t
val zc_zv_mul : _ t -> _ t -> _ t
val zc_divexactu : _ t -> pari_ulong -> _ t
val zc_add : _ t -> _ t -> _ t
val zc_copy : _ t -> _ t
val zc_hnfremdiv : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val zc_is_ei : _ t -> Signed.long
val zc_lincomb : _ t -> _ t -> _ t -> _ t -> _ t
val zc_lincomb1_inplace : _ t -> _ t -> _ t -> unit
val zc_lincomb1_inplace_i : _ t -> _ t -> _ t -> Signed.long -> unit
val zc_neg : _ t -> _ t
val zc_reducemodlll : _ t -> _ t -> _ t
val zc_reducemodmatrix : _ t -> _ t -> _ t
val zc_sub : _ t -> _ t -> _ t
val zc_z_mul : _ t -> Signed.long -> _ t
val zm_q_mul : _ t -> _ t -> _ t
val zm_z_div : _ t -> _ t -> _ t
val zm_z_divexact : _ t -> _ t -> _ t
val zm_z_mul : _ t -> _ t -> _ t
val zm_add : _ t -> _ t -> _ t
val zm_det_triangular : _ t -> _ t
val zm_diag_mul : _ t -> _ t -> _ t
val zm_divexactu : _ t -> pari_ulong -> _ t
val zm_equal : _ t -> _ t -> int
val zm_equal0 : _ t -> int
val zm_hnfdivrem : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val zm_ishnf : _ t -> int
val zm_isdiagonal : _ t -> int
val zm_isidentity : _ t -> int
val zm_isscalar : _ t -> _ t -> int
val zm_max_lg : _ t -> Signed.long
val zm_mul_diag : _ t -> _ t -> _ t
val zm_multosym : _ t -> _ t -> _ t
val zm_neg : _ t -> _ t
val zm_nm_mul : _ t -> _ t -> _ t
val zm_pow : _ t -> _ t -> _ t
val zm_powu : _ t -> pari_ulong -> _ t
val zm_reducemodlll : _ t -> _ t -> _ t
val zm_reducemodmatrix : _ t -> _ t -> _ t
val zm_sqr : _ t -> _ t
val zm_sub : _ t -> _ t -> _ t
val zm_supnorm : _ t -> _ t
val zm_transmul : _ t -> _ t -> _ t
val zm_transmultosym : _ t -> _ t -> _ t
val zm_togglesign : _ t -> unit
val zm_zm_mul : _ t -> _ t -> _ t
val zmrow_zc_mul : _ t -> _ t -> Signed.long -> _ t
val zmrow_equal0 : _ t -> Signed.long -> int
val zv_abscmp : _ t -> _ t -> int
val zv_cmp : _ t -> _ t -> int
val zv_dotsquare : _ t -> _ t
val zv_max_lg : _ t -> Signed.long
val zv_to_nv : _ t -> _ t
val zv_togglesign : _ t -> unit
val gram_matrix : _ t -> _ t
val nm_z_mul : _ t -> _ t -> _ t
val zm_mul : _ t -> _ t -> _ t
val zm_to_flm : _ t -> pari_ulong -> _ t
val zm_to_zm : _ t -> _ t
val zm_zc_mul : _ t -> _ t -> _ t
val zmv_to_zmv : _ t -> _ t
val zv_abs : _ t -> _ t
val zv_content : _ t -> Signed.long
val zv_dotproduct : _ t -> _ t -> Signed.long
val zv_equal : _ t -> _ t -> int
val zv_equal0 : _ t -> int
val zv_neg : _ t -> _ t
val zv_neg_inplace : _ t -> _ t
val zv_prod : _ t -> Signed.long
val zv_prod_z : _ t -> _ t
val zv_sum : _ t -> Signed.long
val zv_sumpart : _ t -> Signed.long -> Signed.long
val zv_to_flv : _ t -> pari_ulong -> _ t
val zv_z_mul : _ t -> Signed.long -> _ t
val zv_zm_mul : _ t -> _ t -> _ t
val zvv_equal : _ t -> _ t -> int
val kronecker_to_zxqx : _ t -> _ t -> _ t
val kronecker_to_zxx : _ t -> Signed.long -> Signed.long -> _ t
val qx_zx_rem : _ t -> _ t -> _ t
val qx_mul : _ t -> _ t -> _ t
val qx_sqr : _ t -> _ t
val qxqm_mul : _ t -> _ t -> _ t -> _ t
val qxqm_sqr : _ t -> _ t -> _ t
val qxqx_qxq_mul : _ t -> _ t -> _ t -> _ t
val qxqx_mul : _ t -> _ t -> _ t -> _ t
val qxqx_powers : _ t -> Signed.long -> _ t -> _ t
val qxqx_sqr : _ t -> _ t -> _ t
val rgx_check_qx : _ t -> string -> unit
val rgx_check_zx : _ t -> string -> unit
val rgx_check_zxx : _ t -> string -> unit
val z_zx_sub : _ t -> _ t -> _ t
val zx_z_add : _ t -> _ t -> _ t
val zx_z_add_shallow : _ t -> _ t -> _ t
val zx_z_eval : _ t -> _ t -> _ t
val zx_z_mul : _ t -> _ t -> _ t
val zx_z_sub : _ t -> _ t -> _ t
val zx_add : _ t -> _ t -> _ t
val zx_affine : _ t -> _ t -> _ t -> _ t
val zx_copy : _ t -> _ t
val zx_deriv : _ t -> _ t
val zx_digits : _ t -> _ t -> _ t
val zxv_zx_fromdigits : _ t -> _ t -> _ t
val zx_div_by_x_1 : _ t -> _ t Ctypes_static.ptr -> _ t
val zx_divuexact : _ t -> pari_ulong -> _ t
val zx_equal : _ t -> _ t -> int
val zx_eval1 : _ t -> _ t
val zx_max_lg : _ t -> Signed.long
val zx_mod_xnm1 : _ t -> pari_ulong -> _ t
val zx_mul : _ t -> _ t -> _ t
val zx_mulspec : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val zx_mulu : _ t -> pari_ulong -> _ t
val zx_neg : _ t -> _ t
val zx_rem : _ t -> _ t -> _ t
val zx_remi2n : _ t -> Signed.long -> _ t
val zx_rescale2n : _ t -> Signed.long -> _ t
val zx_rescale : _ t -> _ t -> _ t
val zx_rescale_lt : _ t -> _ t
val zx_shifti : _ t -> Signed.long -> _ t
val zx_sqr : _ t -> _ t
val zx_sqrspec : _ t -> Signed.long -> _ t
val zx_sub : _ t -> _ t -> _ t
val zx_val : _ t -> Signed.long
val zx_valrem : _ t -> _ t Ctypes_static.ptr -> Signed.long
val zxc_to_flxc : _ t -> pari_ulong -> Signed.long -> _ t
val zxm_to_flxm : _ t -> pari_ulong -> Signed.long -> _ t
val zxqm_mul : _ t -> _ t -> _ t -> _ t
val zxqm_sqr : _ t -> _ t -> _ t
val zxqx_zxq_mul : _ t -> _ t -> _ t -> _ t
val zxqx_sqr : _ t -> _ t -> _ t
val zxqx_mul : _ t -> _ t -> _ t -> _ t
val zxt_remi2n : _ t -> Signed.long -> _ t
val zxv_z_mul : _ t -> _ t -> _ t
val zxv_dotproduct : _ t -> _ t -> _ t
val zxv_equal : _ t -> _ t -> int
val zxv_remi2n : _ t -> Signed.long -> _ t
val zxx_z_divexact : _ t -> _ t -> _ t
val zxx_z_mul : _ t -> _ t -> _ t
val zxx_z_add_shallow : _ t -> _ t -> _ t
val zxx_evalx0 : _ t -> _ t
val zxx_max_lg : _ t -> Signed.long
val zxx_mul_kronecker : _ t -> _ t -> Signed.long -> _ t
val zxx_renormalize : _ t -> Signed.long -> _ t
val zxx_sqr_kronecker : _ t -> Signed.long -> _ t
val rgxx_to_kronecker : _ t -> Signed.long -> _ t
val rgxx_to_kronecker_spec : _ t -> Signed.long -> Signed.long -> _ t
val zxn_mul : _ t -> _ t -> Signed.long -> _ t
val zxn_sqr : _ t -> Signed.long -> _ t
val scalar_zx : _ t -> Signed.long -> _ t
val scalar_zx_shallow : _ t -> Signed.long -> _ t
val zx_to_zx : _ t -> _ t
val zx_z_divexact : _ t -> Signed.long -> _ t
val alg_centralproj : _ t -> _ t -> Signed.long -> _ t
val alg_complete : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val alg_csa_table : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val alg_cyclic : _ t -> _ t -> _ t -> Signed.long -> _ t
val alg_get_absdim : _ t -> Signed.long
val alg_get_abssplitting : _ t -> _ t
val alg_get_aut : _ t -> _ t
val algaut : _ t -> _ t
val alg_get_auts : _ t -> _ t
val alg_get_b : _ t -> _ t
val algb : _ t -> _ t
val algcenter : _ t -> _ t
val alg_get_center : _ t -> _ t
val alg_get_char : _ t -> _ t
val algchar : _ t -> _ t
val alg_get_degree : _ t -> Signed.long
val algdegree : _ t -> Signed.long
val alg_get_dim : _ t -> Signed.long
val algdim : _ t -> Signed.long -> Signed.long
val alg_get_hasse_f : _ t -> _ t
val alghassef : _ t -> _ t
val alg_get_hasse_i : _ t -> _ t
val alghassei : _ t -> _ t
val alg_get_invbasis : _ t -> _ t
val alginvbasis : _ t -> _ t
val alg_get_multable : _ t -> _ t
val alg_get_basis : _ t -> _ t
val algbasis : _ t -> _ t
val alg_get_relmultable : _ t -> _ t
val algrelmultable : _ t -> _ t
val alg_get_splitpol : _ t -> _ t
val alg_get_splittingfield : _ t -> _ t
val algsplittingfield : _ t -> _ t
val alg_get_splittingbasis : _ t -> _ t
val alg_get_splittingbasisinv : _ t -> _ t
val alg_get_splittingdata : _ t -> _ t
val algsplittingdata : _ t -> _ t
val alg_get_tracebasis : _ t -> _ t

val alg_hasse :
  _ t -> Signed.long -> _ t -> _ t -> Signed.long -> Signed.long -> _ t

val alg_hilbert : _ t -> _ t -> _ t -> Signed.long -> Signed.long -> _ t
val alg_matrix : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val alg_model : _ t -> _ t -> Signed.long
val alg_quotient : _ t -> _ t -> Signed.long -> _ t
val algradical : _ t -> _ t
val algsimpledec : _ t -> Signed.long -> _ t
val algsimpledec_ss : _ t -> Signed.long -> _ t
val algsubalg : _ t -> _ t -> _ t
val alg_type : _ t -> Signed.long
val algadd : _ t -> _ t -> _ t -> _ t
val algalgtobasis : _ t -> _ t -> _ t
val algbasistoalg : _ t -> _ t -> _ t
val algcharpoly : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val algdisc : _ t -> _ t
val algdivl : _ t -> _ t -> _ t -> _ t
val algdivr : _ t -> _ t -> _ t -> _ t
val alggroup : _ t -> _ t -> _ t
val alggroupcenter : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val alghasse : _ t -> _ t -> _ t
val alginit : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val algindex : _ t -> _ t -> Signed.long
val alginv : _ t -> _ t -> _ t
val algisassociative : _ t -> _ t -> int
val algiscommutative : _ t -> int
val algisdivision : _ t -> _ t -> int
val algisramified : _ t -> _ t -> int
val algissemisimple : _ t -> int
val algissimple : _ t -> Signed.long -> int
val algissplit : _ t -> _ t -> int
val algisdivl : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> int
val algisinv : _ t -> _ t -> _ t Ctypes_static.ptr -> int
val algmakeintegral : _ t -> Signed.long -> _ t
val algmul : _ t -> _ t -> _ t -> _ t
val algmultable : _ t -> _ t
val alglat_get_primbasis : _ t -> _ t
val alglat_get_scalar : _ t -> _ t
val alglatadd : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val alglatcontains : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> int
val alglatelement : _ t -> _ t -> _ t -> _ t
val alglathnf : _ t -> _ t -> _ t -> _ t
val alglatindex : _ t -> _ t -> _ t -> _ t
val alglatinter : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val alglatmul : _ t -> _ t -> _ t -> _ t
val alglatlefttransporter : _ t -> _ t -> _ t -> _ t
val alglatrighttransporter : _ t -> _ t -> _ t -> _ t
val alglatsubset : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> int
val algneg : _ t -> _ t -> _ t
val algnorm : _ t -> _ t -> Signed.long -> _ t
val algpoleval : _ t -> _ t -> _ t -> _ t
val algpow : _ t -> _ t -> _ t -> _ t
val algprimesubalg : _ t -> _ t
val algramifiedplaces : _ t -> _ t
val algrandom : _ t -> _ t -> _ t
val algsplit : _ t -> Signed.long -> _ t
val algtomatrix : _ t -> _ t -> Signed.long -> _ t
val algsqr : _ t -> _ t -> _ t
val algsub : _ t -> _ t -> _ t -> _ t
val algtableinit : _ t -> _ t -> _ t
val algtensor : _ t -> _ t -> Signed.long -> _ t
val algtrace : _ t -> _ t -> Signed.long -> _ t
val algtype : _ t -> Signed.long
val bnfgwgeneric : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val checkalg : _ t -> unit
val checkhasse : _ t -> _ t -> _ t -> Signed.long -> unit
val checklat : _ t -> _ t -> unit
val conjclasses_algcenter : _ t -> _ t -> _ t
val galoischardet : _ t -> _ t -> Signed.long -> _ t
val galoischarpoly : _ t -> _ t -> Signed.long -> _ t
val galoischartable : _ t -> _ t
val nfgrunwaldwang : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val nfgwkummer : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val f2ms_colelim : _ t -> Signed.long -> _ t
val f2m_image : _ t -> _ t
val f2m_indexrank : _ t -> _ t
val f2m_suppl : _ t -> _ t
val f2xqm_f2xqc_gauss : _ t -> _ t -> _ t -> _ t
val f2xqm_f2xqc_invimage : _ t -> _ t -> _ t -> _ t
val f2xqm_f2xqc_mul : _ t -> _ t -> _ t -> _ t
val f2xqm_deplin : _ t -> _ t -> _ t
val f2xqm_det : _ t -> _ t -> _ t
val f2xqm_gauss : _ t -> _ t -> _ t -> _ t
val f2xqm_ker : _ t -> _ t -> _ t
val f2xqm_image : _ t -> _ t -> _ t
val f2xqm_indexrank : _ t -> _ t -> _ t
val f2xqm_inv : _ t -> _ t -> _ t
val f2xqm_invimage : _ t -> _ t -> _ t -> _ t
val f2xqm_mul : _ t -> _ t -> _ t -> _ t
val f2xqm_rank : _ t -> _ t -> Signed.long
val f2xqm_suppl : _ t -> _ t -> _ t
val flm_image : _ t -> pari_ulong -> _ t
val flm_indexrank : _ t -> pari_ulong -> _ t
val flm_suppl : _ t -> pari_ulong -> _ t
val flxqm_flxqc_gauss : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqm_flxqc_invimage : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqm_flxqc_mul : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqm_deplin : _ t -> _ t -> pari_ulong -> _ t
val flxqm_det : _ t -> _ t -> pari_ulong -> _ t
val flxqm_gauss : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqm_ker : _ t -> _ t -> pari_ulong -> _ t
val flxqm_image : _ t -> _ t -> pari_ulong -> _ t
val flxqm_indexrank : _ t -> _ t -> pari_ulong -> _ t
val flxqm_inv : _ t -> _ t -> pari_ulong -> _ t
val flxqm_invimage : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqm_mul : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqm_rank : _ t -> _ t -> pari_ulong -> Signed.long
val flxqm_suppl : _ t -> _ t -> pari_ulong -> _ t
val fpm_fpc_gauss : _ t -> _ t -> _ t -> _ t
val fpm_fpc_invimage : _ t -> _ t -> _ t -> _ t
val fpm_deplin : _ t -> _ t -> _ t
val fpm_det : _ t -> _ t -> _ t
val fpm_gauss : _ t -> _ t -> _ t -> _ t
val fpm_image : _ t -> _ t -> _ t
val fpm_indexrank : _ t -> _ t -> _ t
val fpm_intersect : _ t -> _ t -> _ t -> _ t
val fpm_intersect_i : _ t -> _ t -> _ t -> _ t
val fpm_inv : _ t -> _ t -> _ t
val fpm_invimage : _ t -> _ t -> _ t -> _ t
val fpm_ker : _ t -> _ t -> _ t
val fpm_rank : _ t -> _ t -> Signed.long
val fpm_suppl : _ t -> _ t -> _ t
val fqm_fqc_gauss : _ t -> _ t -> _ t -> _ t -> _ t
val fqm_fqc_invimage : _ t -> _ t -> _ t -> _ t -> _ t
val fqm_fqc_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fqm_deplin : _ t -> _ t -> _ t -> _ t
val fqm_det : _ t -> _ t -> _ t -> _ t
val fqm_gauss : _ t -> _ t -> _ t -> _ t -> _ t
val fqm_ker : _ t -> _ t -> _ t -> _ t
val fqm_image : _ t -> _ t -> _ t -> _ t
val fqm_indexrank : _ t -> _ t -> _ t -> _ t
val fqm_inv : _ t -> _ t -> _ t -> _ t
val fqm_invimage : _ t -> _ t -> _ t -> _ t -> _ t
val fqm_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fqm_rank : _ t -> _ t -> _ t -> Signed.long
val fqm_suppl : _ t -> _ t -> _ t -> _ t
val qm_image_shallow : _ t -> _ t
val qm_image : _ t -> _ t
val qm_gauss : _ t -> _ t -> _ t
val qm_gauss_i : _ t -> _ t -> Signed.long -> _ t
val qm_indexrank : _ t -> _ t
val qm_inv : _ t -> _ t
val qm_rank : _ t -> Signed.long
val rgm_fp_init : _ t -> _ t -> pari_ulong Ctypes_static.ptr -> _ t
val rgm_hadamard : _ t -> _ t
val rgm_rgc_invimage : _ t -> _ t -> _ t
val rgm_diagonal : _ t -> _ t
val rgm_diagonal_shallow : _ t -> _ t
val rgm_inv : _ t -> _ t
val rgm_inv_upper : _ t -> _ t
val rgm_invimage : _ t -> _ t -> _ t
val rgm_solve : _ t -> _ t -> _ t
val rgm_solve_realimag : _ t -> _ t -> _ t

val rgms_structelim :
  _ t ->
  Signed.long ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  unit

val zm_det : _ t -> _ t
val zm_detmult : _ t -> _ t
val zm_gauss : _ t -> _ t -> _ t
val zm_ker : _ t -> _ t
val zm_imagecompl : _ t -> _ t
val zm_indeximage : _ t -> _ t
val zm_indexrank : _ t -> _ t
val zm_inv : _ t -> _ t Ctypes_static.ptr -> _ t
val zm_inv_ratlift : _ t -> _ t Ctypes_static.ptr -> _ t
val zm_pseudoinv : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t
val zm_rank : _ t -> Signed.long
val zlm_gauss : _ t -> _ t -> pari_ulong -> Signed.long -> _ t -> _ t
val closemodinvertible : _ t -> _ t -> _ t
val deplin : _ t -> _ t
val det : _ t -> _ t
val det0 : _ t -> Signed.long -> _ t
val det2 : _ t -> _ t
val detint : _ t -> _ t
val eigen : _ t -> Signed.long -> _ t
val gauss : _ t -> _ t -> _ t
val gaussmodulo : _ t -> _ t -> _ t -> _ t
val gaussmodulo2 : _ t -> _ t -> _ t -> _ t

val gen_gauss :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_gauss_pivot :
  _ t ->
  Signed.long Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_det :
  _ t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_ker :
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_matcolinvimage :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_matcolmul :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_matinvimage :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_matmul :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  _ t

val image : _ t -> _ t
val image2 : _ t -> _ t
val imagecompl : _ t -> _ t
val indexrank : _ t -> _ t
val inverseimage : _ t -> _ t -> _ t
val ker : _ t -> _ t
val mateigen : _ t -> Signed.long -> Signed.long -> _ t
val matimage0 : _ t -> Signed.long -> _ t
val matker0 : _ t -> Signed.long -> _ t
val rank : _ t -> Signed.long
val reducemodinvertible : _ t -> _ t -> _ t
val reducemodlll : _ t -> _ t -> _ t
val split_realimag : _ t -> Signed.long -> Signed.long -> _ t
val suppl : _ t -> _ t
val flm_charpoly : _ t -> pari_ulong -> _ t
val flm_hess : _ t -> pari_ulong -> _ t
val fpm_charpoly : _ t -> _ t -> _ t
val fpm_hess : _ t -> _ t -> _ t
val frobeniusform : _ t -> Signed.long -> _ t
val qm_minors_coprime : _ t -> _ t -> _ t
val qm_imz : _ t -> _ t

val qm_imz_all :
  _ t -> _ t Ctypes_static.ptr -> Signed.long -> Signed.long -> _ t

val qm_imz_hnf : _ t -> _ t
val qm_imz_hnfall : _ t -> _ t Ctypes_static.ptr -> Signed.long -> _ t
val qm_imq : _ t -> _ t

val qm_imq_all :
  _ t -> _ t Ctypes_static.ptr -> Signed.long -> Signed.long -> _ t

val qm_imq_hnf : _ t -> _ t
val qm_imq_hnfall : _ t -> _ t Ctypes_static.ptr -> Signed.long -> _ t
val qm_charpoly_zx : _ t -> _ t
val qm_charpoly_zx_bound : _ t -> Signed.long -> _ t
val zm_charpoly : _ t -> _ t
val adj : _ t -> _ t
val adjsafe : _ t -> _ t
val caract : _ t -> Signed.long -> _ t
val caradj : _ t -> Signed.long -> _ t Ctypes_static.ptr -> _ t
val carberkowitz : _ t -> Signed.long -> _ t
val carhess : _ t -> Signed.long -> _ t
val charpoly : _ t -> Signed.long -> _ t
val charpoly0 : _ t -> Signed.long -> Signed.long -> _ t
val gnorm : _ t -> _ t
val gnorml1 : _ t -> Signed.long -> _ t
val gnorml1_fake : _ t -> _ t
val gnormlp : _ t -> _ t -> Signed.long -> _ t
val gnorml2 : _ t -> _ t
val gsupnorm : _ t -> Signed.long -> _ t

val gsupnorm_aux :
  _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> Signed.long -> unit

val gtrace : _ t -> _ t
val hess : _ t -> _ t
val intersect : _ t -> _ t -> _ t
val jacobi : _ t -> Signed.long -> _ t
val matadjoint0 : _ t -> Signed.long -> _ t
val matcompanion : _ t -> _ t
val matrixqz0 : _ t -> _ t -> _ t
val minpoly : _ t -> Signed.long -> _ t
val qfgaussred : _ t -> _ t
val qfgaussred_positive : _ t -> _ t
val qfsign : _ t -> _ t
val apply0 : _ t -> _ t -> _ t
val diagonal : _ t -> _ t
val diagonal_shallow : _ t -> _ t
val extract0 : _ t -> _ t -> _ t -> _ t
val fold0 : _ t -> _ t -> _ t

val genapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val genfold :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val genindexselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val genselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val gtomat : _ t -> _ t
val gtrans : _ t -> _ t
val matmuldiagonal : _ t -> _ t -> _ t
val matmultodiagonal : _ t -> _ t -> _ t

val matslice0 :
  _ t -> Signed.long -> Signed.long -> Signed.long -> Signed.long -> _ t

val parapply : _ t -> _ t -> _ t

val parfor :
  _ t ->
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parfor_init :
  parfor_t Ctypes.structure Ctypes_static.ptr -> _ t -> _ t -> _ t -> unit

val parfor_next : parfor_t Ctypes.structure Ctypes_static.ptr -> _ t
val parfor_stop : parfor_t Ctypes.structure Ctypes_static.ptr -> unit

val parforeach :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforeach_init :
  parforeach_t Ctypes.structure Ctypes_static.ptr -> _ t -> _ t -> unit

val parforeach_next : parforeach_t Ctypes.structure Ctypes_static.ptr -> _ t
val parforeach_stop : parforeach_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprime :
  _ t ->
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprime_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr -> _ t -> _ t -> _ t -> unit

val parforprime_next : parforprime_t Ctypes.structure Ctypes_static.ptr -> _ t
val parforprime_stop : parforprime_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprimestep :
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprimestep_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr ->
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  unit

val parforvec :
  _ t ->
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforvec_init :
  parforvec_t Ctypes.structure Ctypes_static.ptr ->
  _ t ->
  _ t ->
  Signed.long ->
  unit

val parforvec_next : parforvec_t Ctypes.structure Ctypes_static.ptr -> _ t
val parforvec_stop : parforvec_t Ctypes.structure Ctypes_static.ptr -> unit
val parselect : _ t -> _ t -> Signed.long -> _ t
val select0 : _ t -> _ t -> Signed.long -> _ t
val shallowextract : _ t -> _ t -> _ t
val shallowmatextract : _ t -> _ t -> _ t -> _ t
val shallowtrans : _ t -> _ t

val vecapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val veccatapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val veccatselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val vecrange : _ t -> _ t -> _ t
val vecrangess : Signed.long -> Signed.long -> _ t

val vecselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val vecselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val vecslice0 : _ t -> Signed.long -> Signed.long -> _ t
val vecsum : _ t -> _ t
val zv_diagonal : _ t -> _ t
val addhelp : string -> string -> unit
val arity0 : _ t -> _ t
val alias0 : string -> string -> unit
val compile_str : string -> _ t
val delete_var : unit -> Signed.long
val fetch_user_var : string -> Signed.long
val fetch_var : unit -> Signed.long
val fetch_var_higher : unit -> Signed.long
val fetch_var_value : Signed.long -> _ t -> _ t
val gp_embedded : string -> string
val gp_embedded_init : Signed.long -> Signed.long -> unit
val gp_read_str : string -> _ t
val gp_read_str_bitprec : string -> Signed.long -> _ t
val gp_read_str_prec : string -> Signed.long -> _ t

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
val readseq : string -> _ t
val safegel : _ t -> Signed.long -> _ t Ctypes_static.ptr
val safeel : _ t -> Signed.long -> Signed.long Ctypes_static.ptr
val safelistel : _ t -> Signed.long -> _ t Ctypes_static.ptr
val safegcoeff : _ t -> Signed.long -> Signed.long -> _ t Ctypes_static.ptr
val strtoi : string -> _ t
val strtor : string -> Signed.long -> _ t
val varhigher : string -> Signed.long -> _ t
val varlower : string -> Signed.long -> _ t
val divisorslenstra : _ t -> _ t -> _ t -> _ t
val isprimeaprcl : _ t -> Signed.long
val qfb0 : _ t -> _ t -> _ t -> _ t

val check_quaddisc :
  _ t ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  string ->
  unit

val check_quaddisc_imag : _ t -> Signed.long Ctypes_static.ptr -> string -> unit
val check_quaddisc_real : _ t -> Signed.long Ctypes_static.ptr -> string -> unit

val cornacchia :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> Signed.long

val cornacchia2 :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> Signed.long

val cornacchia2_sqrt :
  _ t ->
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long

val nucomp : _ t -> _ t -> _ t -> _ t
val nudupl : _ t -> _ t -> _ t
val nupow : _ t -> _ t -> _ t -> _ t
val primeform : _ t -> _ t -> _ t
val primeform_u : _ t -> pari_ulong -> _ t
val qfb_1 : _ t -> _ t
val qfbcomp : _ t -> _ t -> _ t
val qfbcomp_i : _ t -> _ t -> _ t
val qfbcompraw : _ t -> _ t -> _ t
val qfbcompraw_i : _ t -> _ t -> _ t
val qfbcornacchia : _ t -> _ t -> _ t
val qfbpow : _ t -> _ t -> _ t
val qfbpow_i : _ t -> _ t -> _ t
val qfbpowraw : _ t -> Signed.long -> _ t
val qfbpows : _ t -> Signed.long -> _ t
val qfbred : _ t -> _ t
val qfbred_i : _ t -> _ t
val qfbred0 : _ t -> Signed.long -> _ t -> _ t -> _ t
val qfbredsl2 : _ t -> _ t -> _ t
val qfbsolve : _ t -> _ t -> Signed.long -> _ t
val qfbsqr : _ t -> _ t
val qfbsqr_i : _ t -> _ t
val qfisolvep : _ t -> _ t -> _ t
val qfr3_comp : _ t -> _ t -> qfr_data Ctypes.structure Ctypes_static.ptr -> _ t
val qfr3_compraw : _ t -> _ t -> _ t
val qfr3_pow : _ t -> _ t -> qfr_data Ctypes.structure Ctypes_static.ptr -> _ t
val qfr3_red : _ t -> qfr_data Ctypes.structure Ctypes_static.ptr -> _ t
val qfr3_rho : _ t -> qfr_data Ctypes.structure Ctypes_static.ptr -> _ t
val qfr3_to_qfr : _ t -> _ t -> _ t
val qfr5_comp : _ t -> _ t -> qfr_data Ctypes.structure Ctypes_static.ptr -> _ t
val qfr5_compraw : _ t -> _ t -> _ t
val qfr5_dist : _ t -> _ t -> Signed.long -> _ t
val qfr5_pow : _ t -> _ t -> qfr_data Ctypes.structure Ctypes_static.ptr -> _ t
val qfr5_red : _ t -> qfr_data Ctypes.structure Ctypes_static.ptr -> _ t
val qfr5_rho : _ t -> qfr_data Ctypes.structure Ctypes_static.ptr -> _ t
val qfr5_to_qfr : _ t -> _ t -> _ t -> _ t

val qfr_data_init :
  _ t -> Signed.long -> qfr_data Ctypes.structure Ctypes_static.ptr -> unit

val qfr_to_qfr5 : _ t -> Signed.long -> _ t
val qfrsolvep : _ t -> _ t -> _ t
val quadgen : _ t -> _ t
val quadgen0 : _ t -> Signed.long -> _ t
val quadpoly : _ t -> _ t
val quadpoly_i : _ t -> _ t
val quadpoly0 : _ t -> Signed.long -> _ t
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
val fl_powers : pari_ulong -> Signed.long -> pari_ulong -> _ t
val fl_powers_pre : pari_ulong -> Signed.long -> pari_ulong -> pari_ulong -> _ t
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

val fp_2gener : _ t -> _ t
val fp_2gener_i : _ t -> _ t -> _ t
val fp_factored_order : _ t -> _ t -> _ t -> _ t
val fp_ispower : _ t -> _ t -> _ t -> int
val fp_log : _ t -> _ t -> _ t -> _ t -> _ t
val fp_order : _ t -> _ t -> _ t -> _ t
val fp_pow : _ t -> _ t -> _ t -> _ t
val fp_pow_init : _ t -> _ t -> Signed.long -> _ t -> _ t
val fp_pow_table : _ t -> _ t -> _ t -> _ t
val fp_powers : _ t -> Signed.long -> _ t -> _ t
val fp_pows : _ t -> Signed.long -> _ t -> _ t
val fp_powu : _ t -> pari_ulong -> _ t -> _ t
val fp_sqrt : _ t -> _ t -> _ t
val fp_sqrt_i : _ t -> _ t -> _ t -> _ t
val fp_sqrtn : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fpv_prod : _ t -> _ t -> _ t
val z_zv_mod : _ t -> _ t -> _ t
val z_zv_mod_tree : _ t -> _ t -> _ t -> _ t
val z_chinese : _ t -> _ t -> _ t -> _ t -> _ t
val z_chinese_all : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val z_chinese_coprime : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val z_chinese_post : _ t -> _ t -> _ t -> _ t -> _ t -> _ t

val z_chinese_pre :
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  unit

val z_factor_listp : _ t -> _ t -> _ t
val z_nv_mod : _ t -> _ t -> _ t
val zm_nv_mod_tree : _ t -> _ t -> _ t -> _ t
val zv_allpnqn : _ t -> _ t
val zv_chinese : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val zv_chinese_tree : _ t -> _ t -> _ t -> _ t -> _ t
val zv_chinesetree : _ t -> _ t -> _ t
val zv_nv_mod_tree : _ t -> _ t -> _ t -> _ t
val zv_producttree : _ t -> _ t
val zx_nv_mod_tree : _ t -> _ t -> _ t -> _ t
val zxc_nv_mod_tree : _ t -> _ t -> _ t -> Signed.long -> _ t
val zxm_nv_mod_tree : _ t -> _ t -> _ t -> Signed.long -> _ t
val zxx_nv_mod_tree : _ t -> _ t -> _ t -> Signed.long -> _ t
val zideallog : _ t -> _ t -> _ t
val bestappr : _ t -> _ t -> _ t
val bestapprpade : _ t -> Signed.long -> _ t
val chinese : _ t -> _ t -> _ t
val chinese1 : _ t -> _ t
val chinese1_coprime_z : _ t -> _ t
val contfrac0 : _ t -> _ t -> Signed.long -> _ t
val contfracpnqn : _ t -> Signed.long -> _ t
val fibo : Signed.long -> _ t
val gboundcf : _ t -> Signed.long -> _ t
val gcf : _ t -> _ t
val gcf2 : _ t -> _ t -> _ t

val get_fp_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  bb_field Ctypes.structure Ctypes_static.ptr

val hilbert : _ t -> _ t -> _ t -> Signed.long
val hilbertii : _ t -> _ t -> _ t -> Signed.long
val istotient : _ t -> _ t Ctypes_static.ptr -> Signed.long
val krois : _ t -> Signed.long -> Signed.long
val kroiu : _ t -> pari_ulong -> Signed.long
val kronecker : _ t -> _ t -> Signed.long
val krosi : Signed.long -> _ t -> Signed.long
val kross : Signed.long -> Signed.long -> Signed.long
val kroui : pari_ulong -> _ t -> Signed.long
val krouu : pari_ulong -> pari_ulong -> Signed.long
val lcmii : _ t -> _ t -> _ t
val fp_invgen : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val logint0 : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val logintall : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val mpfact : Signed.long -> _ t
val factorial_fl : Signed.long -> pari_ulong -> pari_ulong
val factorial_fp : Signed.long -> _ t -> _ t
val muls_interval : Signed.long -> Signed.long -> _ t
val mulu_interval : pari_ulong -> pari_ulong -> _ t
val mulu_interval_step : pari_ulong -> pari_ulong -> pari_ulong -> _ t
val ncv_chinese_center : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val ncv_chinese_center_tree : _ t -> _ t -> _ t -> _ t -> _ t
val nmv_chinese_center : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val nmv_chinese_center_tree : _ t -> _ t -> _ t -> _ t -> _ t
val nonsquare_fl : pari_ulong -> pari_ulong
val nxcv_chinese_center : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val nxcv_chinese_center_tree : _ t -> _ t -> _ t -> _ t -> _ t
val nxmv_chinese_center : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val nxv_chinese_center : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val nxv_chinese_center_tree : _ t -> _ t -> _ t -> _ t -> _ t
val zv_chinese_center : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val odd_prime_divisors : _ t -> _ t
val pgener_fl : pari_ulong -> pari_ulong
val pgener_fl_local : pari_ulong -> _ t -> pari_ulong
val pgener_fp : _ t -> _ t
val pgener_fp_local : _ t -> _ t -> _ t
val pgener_zl : pari_ulong -> pari_ulong
val pgener_zp : _ t -> _ t
val pnqn : _ t -> _ t
val ramanujantau : _ t -> Signed.long -> _ t
val rootsof1_fl : pari_ulong -> pari_ulong -> pari_ulong
val rootsof1_fp : _ t -> _ t -> _ t
val rootsof1u_fp : pari_ulong -> _ t -> _ t

val u_chinese_coprime :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val znlog : _ t -> _ t -> _ t -> _ t
val znorder : _ t -> _ t -> _ t
val znprimroot : _ t -> _ t
val znstar : _ t -> _ t
val znstar0 : _ t -> Signed.long -> _ t
val rgv_is_zvpos : _ t -> int
val rgv_is_zvnon0 : _ t -> int
val rgv_is_prv : _ t -> int
val z_issquarefree_fact : _ t -> Signed.long

val z_lsmoothen :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val z_smoothen :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val bigomega : _ t -> Signed.long
val bigomegau : pari_ulong -> Signed.long
val boundfact : _ t -> pari_ulong -> _ t
val check_arith_pos : _ t -> string -> _ t
val check_arith_non0 : _ t -> string -> _ t
val check_arith_all : _ t -> string -> _ t
val clean_z_factor : _ t -> _ t
val core : _ t -> _ t

val coredisc2_fact :
  _ t -> Signed.long -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val coredisc2u_fact :
  _ t ->
  Signed.long ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  pari_ulong

val corepartial : _ t -> Signed.long -> _ t
val core0 : _ t -> Signed.long -> _ t
val core2 : _ t -> _ t
val core2partial : _ t -> Signed.long -> _ t
val coredisc : _ t -> _ t
val coredisc0 : _ t -> Signed.long -> _ t
val coredisc2 : _ t -> _ t
val corediscs : Signed.long -> pari_ulong Ctypes_static.ptr -> Signed.long
val divisors : _ t -> _ t
val divisors_factored : _ t -> _ t
val divisors0 : _ t -> Signed.long -> _ t
val divisorsu : pari_ulong -> _ t
val divisorsu_moebius : _ t -> _ t
val divisorsu_fact : _ t -> _ t
val divisorsu_fact_factored : _ t -> _ t
val eulerphi : _ t -> _ t
val eulerphiu : pari_ulong -> pari_ulong
val eulerphiu_fact : _ t -> pari_ulong
val factor_pn_1 : _ t -> pari_ulong -> _ t
val factor_pn_1_limit : _ t -> Signed.long -> pari_ulong -> _ t
val factoru_pow : pari_ulong -> _ t
val fuse_z_factor : _ t -> _ t -> _ t
val is_z_factor : _ t -> int
val is_z_factornon0 : _ t -> int
val is_z_factorpos : _ t -> int
val is_nf_factor : _ t -> int
val is_nf_extfactor : _ t -> int
val issquarefree : _ t -> Signed.long
val numdiv : _ t -> _ t
val numdivu : Signed.long -> Signed.long
val numdivu_fact : _ t -> Signed.long
val omega : _ t -> Signed.long
val omegau : pari_ulong -> Signed.long
val sumdiv : _ t -> _ t
val sumdivk : _ t -> Signed.long -> _ t
val uissquarefree : pari_ulong -> Signed.long
val uissquarefree_fact : _ t -> Signed.long
val usumdiv_fact : _ t -> _ t
val usumdivk_fact : _ t -> pari_ulong -> _ t
val fpx_fpc_nfpoleval : _ t -> _ t -> _ t -> _ t -> _ t
val embed_t2 : _ t -> Signed.long -> _ t
val embednorm_t2 : _ t -> Signed.long -> _ t
val embed_norm : _ t -> Signed.long -> _ t
val check_zkmodule_i : _ t -> int
val check_zkmodule : _ t -> string -> unit
val checkbid : _ t -> unit
val checkbid_i : _ t -> _ t
val checkbnf : _ t -> _ t
val checkbnf_i : _ t -> _ t
val checkbnr : _ t -> unit
val checkbnr_i : _ t -> _ t
val checkabgrp : _ t -> unit
val checksqmat : _ t -> Signed.long -> unit
val checknf : _ t -> _ t
val checknf_i : _ t -> _ t
val checknfelt_mod : _ t -> _ t -> string -> _ t
val checkprid : _ t -> unit
val checkprid_i : _ t -> int
val checkrnf : _ t -> unit
val checkrnf_i : _ t -> int
val factoredpolred : _ t -> _ t -> _ t
val factoredpolred2 : _ t -> _ t -> _ t
val galoisapply : _ t -> _ t -> _ t -> _ t
val get_bnf : _ t -> Signed.long Ctypes_static.ptr -> _ t
val get_bnfpol : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t
val get_nf : _ t -> Signed.long Ctypes_static.ptr -> _ t
val get_nfpol : _ t -> _ t Ctypes_static.ptr -> _ t
val get_prid : _ t -> _ t
val idealfrobenius : _ t -> _ t -> _ t -> _ t
val idealfrobenius_aut : _ t -> _ t -> _ t -> _ t -> _ t
val idealramfrobenius : _ t -> _ t -> _ t -> _ t -> _ t
val idealramfrobenius_aut : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val idealramgroups : _ t -> _ t -> _ t -> _ t
val idealramgroups_aut : _ t -> _ t -> _ t -> _ t -> _ t
val nf_get_allroots : _ t -> _ t
val nf_get_prec : _ t -> Signed.long

val nfmaxord_to_nf :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr -> _ t -> Signed.long -> _ t

val nfcertify : _ t -> _ t
val nfgaloismatrix : _ t -> _ t -> _ t
val nfgaloismatrixapply : _ t -> _ t -> _ t -> _ t
val nfgaloispermtobasis : _ t -> _ t -> _ t
val nfinit_basic : nfmaxord_t Ctypes.structure Ctypes_static.ptr -> _ t -> unit

val nfinit_complete :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  _ t

val nfinit : _ t -> Signed.long -> _ t
val nfinit0 : _ t -> Signed.long -> Signed.long -> _ t
val nfinitred : _ t -> Signed.long -> _ t
val nfinitred2 : _ t -> Signed.long -> _ t
val nfisincl : _ t -> _ t -> _ t
val nfisincl0 : _ t -> _ t -> Signed.long -> _ t
val nfisisom : _ t -> _ t -> _ t
val nfnewprec : _ t -> Signed.long -> _ t
val nfnewprec_shallow : _ t -> Signed.long -> _ t
val nfpoleval : _ t -> _ t -> _ t -> _ t
val nfsplitting : _ t -> _ t -> _ t
val nfsplitting0 : _ t -> _ t -> Signed.long -> _ t
val nftyp : _ t -> Signed.long
val polredord : _ t -> _ t
val polred : _ t -> _ t
val polred0 : _ t -> Signed.long -> _ t -> _ t
val polred2 : _ t -> _ t
val polredabs : _ t -> _ t
val polredabs0 : _ t -> Signed.long -> _ t
val polredabs2 : _ t -> _ t
val polredabsall : _ t -> Signed.long -> _ t
val polredbest : _ t -> Signed.long -> _ t
val poltomonic : _ t -> _ t Ctypes_static.ptr -> _ t
val rnfpolredabs : _ t -> _ t -> Signed.long -> _ t
val rnfpolredbest : _ t -> _ t -> Signed.long -> _ t
val smallpolred : _ t -> _ t
val smallpolred2 : _ t -> _ t
val tschirnhaus : _ t -> _ t
val zx_q_mul : _ t -> _ t -> _ t
val zx_q_normalize : _ t -> _ t Ctypes_static.ptr -> _ t
val zx_z_normalize : _ t -> _ t Ctypes_static.ptr -> _ t
val zx_to_monic : _ t -> _ t Ctypes_static.ptr -> _ t
val zx_primitive_to_monic : _ t -> _ t Ctypes_static.ptr -> _ t
val zxx_q_mul : _ t -> _ t -> _ t
val fq_to_nf : _ t -> _ t -> _ t
val fqm_to_nfm : _ t -> _ t -> _ t
val fqv_to_nfv : _ t -> _ t -> _ t
val fqx_to_nfx : _ t -> _ t -> _ t
val rg_nffix : string -> _ t -> _ t -> int -> _ t
val rgv_nffix : string -> _ t -> _ t -> int -> _ t
val rgx_nffix : string -> _ t -> _ t -> int -> _ t
val zx_composedsum : _ t -> _ t -> _ t
val zx_compositum : _ t -> _ t -> Signed.long Ctypes_static.ptr -> _ t
val zpx_disc_val : _ t -> _ t -> Signed.long
val zpx_gcd : _ t -> _ t -> _ t -> _ t -> _ t
val zpx_monic_factor : _ t -> _ t -> Signed.long -> _ t
val zpx_primedec : _ t -> _ t -> _ t
val zpx_reduced_resultant : _ t -> _ t -> _ t -> _ t -> _ t
val zpx_reduced_resultant_fast : _ t -> _ t -> _ t -> Signed.long -> _ t
val zpx_resultant_val : _ t -> _ t -> _ t -> Signed.long -> Signed.long
val checkmodpr : _ t -> unit
val compositum : _ t -> _ t -> _ t
val compositum2 : _ t -> _ t -> _ t
val nfdisc : _ t -> _ t
val get_modpr : _ t -> _ t
val indexpartial : _ t -> _ t -> _ t
val modpr_genfq : _ t -> _ t

val nf_to_fq_init :
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val nf_to_fq : _ t -> _ t -> _ t -> _ t
val nfm_to_fqm : _ t -> _ t -> _ t -> _ t
val nfv_to_fqv : _ t -> _ t -> _ t -> _ t
val nfx_to_fqx : _ t -> _ t -> _ t -> _ t
val nfx_to_monic : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val nfbasis : _ t -> _ t Ctypes_static.ptr -> _ t
val nfcompositum : _ t -> _ t -> _ t -> Signed.long -> _ t
val nfdiscfactors : _ t -> _ t

val nfmaxord :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr -> _ t -> Signed.long -> unit

val nfmodpr : _ t -> _ t -> _ t -> _ t
val nfmodprinit : _ t -> _ t -> _ t
val nfmodprinit0 : _ t -> _ t -> Signed.long -> _ t
val nfmodprlift : _ t -> _ t -> _ t -> _ t
val nfreducemodpr : _ t -> _ t -> _ t -> _ t
val polcompositum0 : _ t -> _ t -> Signed.long -> _ t
val idealprimedec : _ t -> _ t -> _ t
val idealprimedec_galois : _ t -> _ t -> _ t
val idealprimedec_degrees : _ t -> _ t -> _ t
val idealprimedec_kummer : _ t -> _ t -> Signed.long -> _ t -> _ t
val idealprimedec_limit_f : _ t -> _ t -> Signed.long -> _ t
val idealprimedec_limit_norm : _ t -> _ t -> _ t -> _ t
val poldiscfactors : _ t -> Signed.long -> _ t
val rnfbasis : _ t -> _ t -> _ t
val rnfdedekind : _ t -> _ t -> _ t -> Signed.long -> _ t
val rnfdet : _ t -> _ t -> _ t
val rnfdisc_factored : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val rnfdiscf : _ t -> _ t -> _ t
val rnfequation : _ t -> _ t -> _ t
val rnfequation0 : _ t -> _ t -> Signed.long -> _ t
val rnfequation2 : _ t -> _ t -> _ t
val nf_pv_to_prv : _ t -> _ t -> _ t
val nf_rnfeq : _ t -> _ t -> _ t
val nf_rnfeqsimple : _ t -> _ t -> _ t

val rnfequationall :
  _ t -> _ t -> Signed.long Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val rnfhnfbasis : _ t -> _ t -> _ t
val rnfisfree : _ t -> _ t -> Signed.long
val rnflllgram : _ t -> _ t -> _ t -> Signed.long -> _ t
val rnfpolred : _ t -> _ t -> Signed.long -> _ t
val rnfpseudobasis : _ t -> _ t -> _ t
val rnfsimplifybasis : _ t -> _ t -> _ t
val rnfsteinitz : _ t -> _ t -> _ t
val factorial_lval : pari_ulong -> pari_ulong -> Signed.long

val zk_to_fq_init :
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val zk_to_fq : _ t -> _ t -> _ t
val qxqv_to_fpm : _ t -> _ t -> _ t -> _ t
val zkmodprinit : _ t -> _ t -> _ t
val idealstar : _ t -> _ t -> Signed.long -> _ t
val idealstarprk : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val rgc_to_nfc : _ t -> _ t -> _ t
val rgm_rgx_mul : _ t -> _ t -> _ t
val rgm_to_nfm : _ t -> _ t -> _ t
val rgx_to_nfx : _ t -> _ t -> _ t
val zc_nfval : _ t -> _ t -> Signed.long
val zc_nfvalrem : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val zc_prdvd : _ t -> _ t -> int
val zm_zx_mul : _ t -> _ t -> _ t
val zv_snf_gcd : _ t -> _ t -> _ t
val algtobasis : _ t -> _ t -> _ t
val basistoalg : _ t -> _ t -> _ t
val ei_multable : _ t -> Signed.long -> _ t

val get_nf_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  bb_field Ctypes.structure Ctypes_static.ptr

val famat_nfvalrem : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val gpnfvalrem : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val idealfactorback : _ t -> _ t -> _ t -> int -> _ t
val ideallist : _ t -> Signed.long -> _ t
val ideallist0 : _ t -> Signed.long -> Signed.long -> _ t
val gideallist : _ t -> _ t -> Signed.long -> _ t
val ideallistarch : _ t -> _ t -> _ t -> _ t
val ideallog : _ t -> _ t -> _ t -> _ t
val ideallogmod : _ t -> _ t -> _ t -> _ t -> _ t
val ideallog_units : _ t -> _ t -> _ t
val ideallog_units0 : _ t -> _ t -> _ t -> _ t
val idealprincipalunits : _ t -> _ t -> Signed.long -> _ t
val idealstar0 : _ t -> _ t -> Signed.long -> _ t
val idealstarmod : _ t -> _ t -> Signed.long -> _ t -> _ t
val indices_to_vec01 : _ t -> Signed.long -> _ t
val matalgtobasis : _ t -> _ t -> _ t
val matbasistoalg : _ t -> _ t -> _ t
val multable : _ t -> _ t -> _ t
val nf_to_scalar_or_alg : _ t -> _ t -> _ t
val nf_to_scalar_or_basis : _ t -> _ t -> _ t
val nf_cxlog : _ t -> _ t -> Signed.long -> _ t
val nfv_cxlog : _ t -> _ t -> Signed.long -> _ t
val nfadd : _ t -> _ t -> _ t -> _ t
val nfchecksigns : _ t -> _ t -> _ t -> int
val nfdiv : _ t -> _ t -> _ t -> _ t
val nfdiveuc : _ t -> _ t -> _ t -> _ t
val nfdivrem : _ t -> _ t -> _ t -> _ t
val nfembed : _ t -> _ t -> Signed.long -> _ t
val nfeltembed : _ t -> _ t -> _ t -> Signed.long -> _ t
val nfeltembed_i : _ t Ctypes_static.ptr -> _ t -> _ t -> Signed.long -> _ t
val nfeltsign : _ t -> _ t -> _ t -> _ t
val nffactorback : _ t -> _ t -> _ t -> _ t
val nfinv : _ t -> _ t -> _ t
val nfinvmodideal : _ t -> _ t -> _ t -> _ t
val nfissquare : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long

val nfispower :
  _ t -> _ t -> Signed.long -> _ t Ctypes_static.ptr -> Signed.long

val nflogembed : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long -> _ t
val nfm_det : _ t -> _ t -> _ t
val nfm_inv : _ t -> _ t -> _ t
val nfm_ker : _ t -> _ t -> _ t
val nfm_mul : _ t -> _ t -> _ t -> _ t
val nfm_nfc_mul : _ t -> _ t -> _ t -> _ t
val nfmod : _ t -> _ t -> _ t -> _ t
val nfmul : _ t -> _ t -> _ t -> _ t
val nfmuli : _ t -> _ t -> _ t -> _ t
val nfnorm : _ t -> _ t -> _ t
val nfpolsturm : _ t -> _ t -> _ t -> _ t
val nfpow : _ t -> _ t -> _ t -> _ t
val nfpow_u : _ t -> _ t -> pari_ulong -> _ t
val nfpowmodideal : _ t -> _ t -> _ t -> _ t -> _ t
val nfsign : _ t -> _ t -> _ t
val nfsign_arch : _ t -> _ t -> _ t -> _ t
val nfsign_from_logarch : _ t -> _ t -> _ t -> _ t
val nfsqr : _ t -> _ t -> _ t
val nfsqri : _ t -> _ t -> _ t
val nfsub : _ t -> _ t -> _ t -> _ t
val nftrace : _ t -> _ t -> _ t
val nfval : _ t -> _ t -> _ t -> Signed.long
val nfvalrem : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val polmod_nffix : string -> _ t -> _ t -> int -> _ t
val polmod_nffix2 : string -> _ t -> _ t -> _ t -> int -> _ t
val pr_basis_perm : _ t -> _ t -> _ t
val pr_equal : _ t -> _ t -> int
val rnfalgtobasis : _ t -> _ t -> _ t
val rnfbasistoalg : _ t -> _ t -> _ t
val rnfeltnorm : _ t -> _ t -> _ t
val rnfelttrace : _ t -> _ t -> _ t
val set_sign_mod_divisor : _ t -> _ t -> _ t -> _ t -> _ t
val tablemul : _ t -> _ t -> _ t -> _ t
val tablemul_ei : _ t -> _ t -> Signed.long -> _ t
val tablemul_ei_ej : _ t -> Signed.long -> Signed.long -> _ t
val tablemulvec : _ t -> _ t -> _ t -> _ t
val tablesqr : _ t -> _ t -> _ t
val vec01_to_indices : _ t -> _ t
val vecsmall01_to_indices : _ t -> _ t
val zk_inv : _ t -> _ t -> _ t
val zk_multable : _ t -> _ t -> _ t
val zk_scalar_or_multable : _ t -> _ t -> _ t
val zkchinese : _ t -> _ t -> _ t -> _ t
val zkchinese1 : _ t -> _ t -> _ t
val zkchineseinit : _ t -> _ t -> _ t -> _ t -> _ t
val zkmultable_capz : _ t -> _ t
val zkmultable_inv : _ t -> _ t

val fl_invgen :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val rm_round_maxrank : _ t -> _ t
val zm_famat_limit : _ t -> _ t -> _ t
val zv_cba : _ t -> _ t
val zv_cba_extend : _ t -> _ t -> _ t
val z_cba : _ t -> _ t -> _ t
val z_ppgle : _ t -> _ t -> _ t
val z_ppio : _ t -> _ t -> _ t
val z_ppo : _ t -> _ t -> _ t
val famatv_factorback : _ t -> _ t -> _ t
val famatv_zv_factorback : _ t -> _ t -> _ t
val famat_z_gcd : _ t -> _ t -> _ t
val famat_div : _ t -> _ t -> _ t
val famat_div_shallow : _ t -> _ t -> _ t
val famat_idealfactor : _ t -> _ t -> _ t
val famat_inv : _ t -> _ t
val famat_inv_shallow : _ t -> _ t
val famat_makecoprime : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val famat_mul : _ t -> _ t -> _ t
val famat_mul_shallow : _ t -> _ t -> _ t
val famat_mulpow_shallow : _ t -> _ t -> _ t -> _ t
val famat_mulpows_shallow : _ t -> _ t -> Signed.long -> _ t
val famat_pow : _ t -> _ t -> _ t
val famat_pow_shallow : _ t -> _ t -> _ t
val famat_pows_shallow : _ t -> Signed.long -> _ t
val famat_reduce : _ t -> _ t
val famat_remove_trivial : _ t -> _ t
val famat_sqr : _ t -> _ t
val famat_to_nf : _ t -> _ t -> _ t
val famat_to_nf_moddivisor : _ t -> _ t -> _ t -> _ t -> _ t
val famat_to_nf_modideal_coprime : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val famatsmall_reduce : _ t -> _ t
val gpidealfactor : _ t -> _ t -> _ t -> _ t
val gpidealval : _ t -> _ t -> _ t -> _ t

val idealhnf_z_factor :
  _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val idealhnf_z_factor_i :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val idealhnf_inv : _ t -> _ t -> _ t
val idealhnf_inv_z : _ t -> _ t -> _ t
val idealhnf_mul : _ t -> _ t -> _ t -> _ t
val idealadd : _ t -> _ t -> _ t -> _ t
val idealaddmultoone : _ t -> _ t -> _ t
val idealaddtoone : _ t -> _ t -> _ t -> _ t
val idealaddtoone0 : _ t -> _ t -> _ t -> _ t
val idealaddtoone_i : _ t -> _ t -> _ t -> _ t
val idealaddtoone_raw : _ t -> _ t -> _ t -> _ t
val idealappr : _ t -> _ t -> _ t
val idealappr0 : _ t -> _ t -> Signed.long -> _ t
val idealapprfact : _ t -> _ t -> _ t
val idealchinese : _ t -> _ t -> _ t -> _ t
val idealcoprime : _ t -> _ t -> _ t -> _ t
val idealcoprimefact : _ t -> _ t -> _ t -> _ t
val idealdiv : _ t -> _ t -> _ t -> _ t
val idealdiv0 : _ t -> _ t -> _ t -> Signed.long -> _ t
val idealdivexact : _ t -> _ t -> _ t -> _ t
val idealdivpowprime : _ t -> _ t -> _ t -> _ t -> _ t
val idealdown : _ t -> _ t -> _ t
val idealfactor : _ t -> _ t -> _ t
val idealfactor_limit : _ t -> _ t -> pari_ulong -> _ t
val idealfactor_partial : _ t -> _ t -> _ t -> _ t
val idealhnf : _ t -> _ t -> _ t
val idealhnf0 : _ t -> _ t -> _ t -> _ t
val idealhnf_principal : _ t -> _ t -> _ t
val idealhnf_shallow : _ t -> _ t -> _ t
val idealhnf_two : _ t -> _ t -> _ t
val idealintersect : _ t -> _ t -> _ t -> _ t
val idealinv : _ t -> _ t -> _ t
val idealismaximal : _ t -> _ t -> _ t

val idealispower :
  _ t -> _ t -> Signed.long -> _ t Ctypes_static.ptr -> Signed.long

val idealmin : _ t -> _ t -> _ t -> _ t
val idealmul : _ t -> _ t -> _ t -> _ t
val idealmul0 : _ t -> _ t -> _ t -> Signed.long -> _ t
val idealmulpowprime : _ t -> _ t -> _ t -> _ t -> _ t
val idealmulred : _ t -> _ t -> _ t -> _ t
val idealnorm : _ t -> _ t -> _ t
val idealnumden : _ t -> _ t -> _ t
val idealpow : _ t -> _ t -> _ t -> _ t
val idealpow0 : _ t -> _ t -> _ t -> Signed.long -> _ t
val idealpowred : _ t -> _ t -> _ t -> _ t
val idealpows : _ t -> _ t -> Signed.long -> _ t
val idealprod : _ t -> _ t -> _ t
val idealprodprime : _ t -> _ t -> _ t
val idealprodval : _ t -> _ t -> _ t -> Signed.long
val idealpseudomin : _ t -> _ t -> _ t
val idealpseudomin_nonscalar : _ t -> _ t -> _ t
val idealpseudominvec : _ t -> _ t -> _ t
val idealpseudored : _ t -> _ t -> _ t
val idealred0 : _ t -> _ t -> _ t -> _ t
val idealred_elt : _ t -> _ t -> _ t
val idealredmodpower : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val idealsqr : _ t -> _ t -> _ t
val idealtwoelt : _ t -> _ t -> _ t
val idealtwoelt0 : _ t -> _ t -> _ t -> _ t
val idealtwoelt2 : _ t -> _ t -> _ t -> _ t
val idealtyp : _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> Signed.long
val idealval : _ t -> _ t -> _ t -> Signed.long
val isideal : _ t -> _ t -> Signed.long
val matreduce : _ t -> _ t
val nfc_multable_mul : _ t -> _ t -> _ t
val nfc_nf_mul : _ t -> _ t -> _ t -> _ t
val nf_get_gtwist : _ t -> _ t -> _ t
val nf_get_gtwist1 : _ t -> Signed.long -> _ t
val nf_to_fp_coprime : _ t -> _ t -> _ t -> _ t
val nfdetint : _ t -> _ t -> _ t
val nfdivmodpr : _ t -> _ t -> _ t -> _ t -> _ t
val nfhnf : _ t -> _ t -> _ t
val nfhnf0 : _ t -> _ t -> Signed.long -> _ t
val nfhnfmod : _ t -> _ t -> _ t -> _ t
val nfkermodpr : _ t -> _ t -> _ t -> _ t
val nfmulmodpr : _ t -> _ t -> _ t -> _ t -> _ t
val nfpowmodpr : _ t -> _ t -> _ t -> _ t -> _ t
val nfreduce : _ t -> _ t -> _ t -> _ t
val nfsnf : _ t -> _ t -> _ t
val nfsnf0 : _ t -> _ t -> Signed.long -> _ t
val nfsolvemodpr : _ t -> _ t -> _ t -> _ t -> _ t
val prv_lcm_capz : _ t -> _ t
val prv_primes : _ t -> _ t
val pr_hnf : _ t -> _ t -> _ t
val pr_inv : _ t -> _ t
val pr_inv_p : _ t -> _ t
val pr_uniformizer : _ t -> _ t -> _ t
val sunits_makecoprime : _ t -> _ t -> _ t -> _ t
val to_famat : _ t -> _ t -> _ t
val to_famat_shallow : _ t -> _ t -> _ t
val u_ppo : pari_ulong -> pari_ulong -> pari_ulong
val vecdiv : _ t -> _ t -> _ t
val vecinv : _ t -> _ t
val vecmul : _ t -> _ t -> _ t
val vecpow : _ t -> _ t -> _ t
val vecsqr : _ t -> _ t
val zkc_multable_mul : _ t -> _ t -> _ t
val eltreltoabs : _ t -> _ t -> _ t
val eltabstorel : _ t -> _ t -> _ t
val eltabstorel_lift : _ t -> _ t -> _ t
val nf_nfzk : _ t -> _ t -> _ t
val rnf_build_nfabs : _ t -> Signed.long -> _ t
val rnf_zkabs : _ t -> _ t
val nfeltup : _ t -> _ t -> _ t -> _ t
val rnfcomplete : _ t -> unit
val rnfeltabstorel : _ t -> _ t -> _ t
val rnfeltdown : _ t -> _ t -> _ t
val rnfeltdown0 : _ t -> _ t -> Signed.long -> _ t
val rnfeltreltoabs : _ t -> _ t -> _ t
val rnfeltup : _ t -> _ t -> _ t
val rnfeltup0 : _ t -> _ t -> Signed.long -> _ t
val rnfidealabstorel : _ t -> _ t -> _ t
val rnfidealdown : _ t -> _ t -> _ t
val rnfidealfactor : _ t -> _ t -> _ t
val rnfidealhnf : _ t -> _ t -> _ t
val rnfidealmul : _ t -> _ t -> _ t -> _ t
val rnfidealnormabs : _ t -> _ t -> _ t
val rnfidealnormrel : _ t -> _ t -> _ t
val rnfidealprimedec : _ t -> _ t -> _ t
val rnfidealreltoabs : _ t -> _ t -> _ t
val rnfidealreltoabs0 : _ t -> _ t -> Signed.long -> _ t
val rnfidealtwoelement : _ t -> _ t -> _ t
val rnfidealup : _ t -> _ t -> _ t
val rnfidealup0 : _ t -> _ t -> Signed.long -> _ t
val rnfinit : _ t -> _ t -> _ t
val rnfinit0 : _ t -> _ t -> Signed.long -> _ t
val get_arith_zzm : _ t -> _ t
val get_arith_z : _ t -> _ t

val gen_ph_log :
  _ t ->
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_shanks_init :
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_shanks :
  _ t ->
  _ t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_shanks_sqrtn :
  _ t ->
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_gener :
  _ t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_ellgens :
  _ t ->
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t -> _ t -> _ t)
  Ctypes_static.static_funptr ->
  _ t

val gen_ellgroup :
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t -> _ t -> _ t)
  Ctypes_static.static_funptr ->
  _ t

val gen_factored_order :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_order :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_select_order :
  _ t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_plog :
  _ t ->
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  _ t

val gen_pow :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_pow_fold :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_pow_fold_i :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_pow_i :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_pow_init :
  _ t ->
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_pow_table :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_powers :
  _ t ->
  Signed.long ->
  int ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_powu :
  _ t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_powu_fold :
  _ t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_powu_fold_i :
  _ t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_powu_i :
  _ t ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val gen_product :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t

val matdetmod : _ t -> _ t -> _ t
val matimagemod : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val matinvmod : _ t -> _ t -> _ t
val matkermod : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val matsolvemod : _ t -> _ t -> _ t -> Signed.long -> _ t
val bernfrac : Signed.long -> _ t
val bernpol : Signed.long -> Signed.long -> _ t
val bernreal : Signed.long -> Signed.long -> _ t
val bernvec : Signed.long -> _ t
val constbern : Signed.long -> unit
val eulerfrac : Signed.long -> _ t
val eulerpol : Signed.long -> Signed.long -> _ t
val eulerreal : Signed.long -> Signed.long -> _ t
val eulervec : Signed.long -> _ t
val harmonic : pari_ulong -> _ t
val harmonic0 : pari_ulong -> _ t -> _ t

val qr_init :
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long ->
  int

val r_from_qr : _ t -> Signed.long -> _ t
val rgm_babai : _ t -> _ t -> _ t

val rgm_qr_init :
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long ->
  int

val rgm_gram_schmidt : _ t -> _ t Ctypes_static.ptr -> _ t
val algdep : _ t -> Signed.long -> _ t
val algdep0 : _ t -> Signed.long -> Signed.long -> _ t
val bestapprnf : _ t -> _ t -> _ t -> Signed.long -> _ t

val forqfvec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> float -> Signed.long)
  Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  unit

val forqfvec0 : _ t -> _ t -> _ t -> unit

val forqfvec1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  unit

val gaussred_from_qr : _ t -> Signed.long -> _ t
val lindep : _ t -> _ t
val lindep_xadic : _ t -> _ t
val lindep_bit : _ t -> Signed.long -> _ t
val lindep_padic : _ t -> _ t
val lindep0 : _ t -> Signed.long -> _ t
val lindep2 : _ t -> Signed.long -> _ t
val lindepfull_bit : _ t -> Signed.long -> _ t
val mathouseholder : _ t -> _ t -> _ t
val matqr : _ t -> Signed.long -> Signed.long -> _ t
val minim : _ t -> _ t -> _ t -> _ t
val minim_raw : _ t -> _ t -> _ t -> _ t
val minim_zm : _ t -> _ t -> _ t -> _ t
val minim2 : _ t -> _ t -> _ t -> _ t
val qfminim0 : _ t -> _ t -> _ t -> Signed.long -> Signed.long -> _ t
val qfperfection : _ t -> _ t
val qfrep0 : _ t -> _ t -> Signed.long -> _ t
val seralgdep : _ t -> Signed.long -> Signed.long -> _ t
val serdiffdep : _ t -> Signed.long -> Signed.long -> _ t
val vandermondeinverse : _ t -> _ t -> _ t -> _ t -> _ t
val vandermondeinverseinit : _ t -> _ t
val zncoppersmith : _ t -> _ t -> _ t -> _ t -> _ t
val qxq_reverse : _ t -> _ t -> _ t
val vec_equiv : _ t -> _ t
val rgv_polint : _ t -> _ t -> Signed.long -> _ t
val vec_reduce : _ t -> _ t Ctypes_static.ptr -> _ t
val rgxq_reverse : _ t -> _ t -> _ t
val zc_union_shallow : _ t -> _ t -> _ t
val zv_indexsort : _ t -> _ t
val zv_sort : _ t -> _ t
val zv_sort_inplace : _ t -> unit
val zv_sort_shallow : _ t -> _ t
val zv_sort_uniq : _ t -> _ t
val zv_sort_uniq_shallow : _ t -> _ t
val zv_union_shallow : _ t -> _ t -> _ t
val binomial : _ t -> Signed.long -> _ t
val binomial0 : _ t -> _ t -> _ t
val binomialuu : pari_ulong -> pari_ulong -> _ t
val cmp_flx : _ t -> _ t -> int
val cmp_rgx : _ t -> _ t -> int
val cmp_nodata : unit Ctypes_static.ptr -> _ t -> _ t -> int
val cmp_prime_ideal : _ t -> _ t -> int
val cmp_prime_over_p : _ t -> _ t -> int
val cmp_universal : _ t -> _ t -> int
val convol : _ t -> _ t -> _ t
val gen_cmp_rgx : unit Ctypes_static.ptr -> _ t -> _ t -> int
val polcyclo : Signed.long -> Signed.long -> _ t
val polcyclo_eval : Signed.long -> _ t -> _ t
val dirdiv : _ t -> _ t -> _ t
val dirmul : _ t -> _ t -> _ t
val eulerianpol : Signed.long -> Signed.long -> _ t
val gprec_wensure : _ t -> Signed.long -> _ t

val gen_indexsort :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  _ t

val gen_indexsort_uniq :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  _ t

val gen_search :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  Signed.long

val gen_setminus :
  _ t -> _ t -> (_ t -> _ t -> int) Ctypes_static.static_funptr -> _ t

val gen_sort :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  _ t

val gen_sort_inplace :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  _ t Ctypes_static.ptr ->
  unit

val gen_sort_shallow :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  _ t

val gen_sort_uniq :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  _ t

val getstack : unit -> Signed.long
val gettime : unit -> Signed.long
val getabstime : unit -> Signed.long
val getwalltime : unit -> _ t
val gprec : _ t -> Signed.long -> _ t
val gprec_wtrunc : _ t -> Signed.long -> _ t
val gprec_w : _ t -> Signed.long -> _ t
val gtoset : _ t -> _ t
val indexlexsort : _ t -> _ t
val indexsort : _ t -> _ t
val indexvecsort : _ t -> _ t -> _ t
val laplace : _ t -> _ t
val lexsort : _ t -> _ t
val mathilbert : Signed.long -> _ t
val matqpascal : Signed.long -> _ t -> _ t

val merge_factor :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  _ t

val merge_sort_uniq :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  _ t

val modreverse : _ t -> _ t
val polhermite : Signed.long -> Signed.long -> _ t
val polhermite_eval0 : Signed.long -> _ t -> Signed.long -> _ t
val polhermite_eval : Signed.long -> _ t -> _ t
val pollaguerre : Signed.long -> _ t -> Signed.long -> _ t
val pollaguerre_eval : Signed.long -> _ t -> _ t -> _ t
val pollaguerre_eval0 : Signed.long -> _ t -> _ t -> Signed.long -> _ t
val pollegendre : Signed.long -> Signed.long -> _ t
val pollegendre_reduced : Signed.long -> Signed.long -> _ t
val pollegendre_eval : Signed.long -> _ t -> _ t
val pollegendre_eval0 : Signed.long -> _ t -> Signed.long -> _ t
val polint : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val polint_i : _ t -> _ t -> _ t -> Signed.long Ctypes_static.ptr -> _ t

val polintspec :
  _ t -> _ t -> _ t -> Signed.long -> Signed.long Ctypes_static.ptr -> _ t

val polchebyshev : Signed.long -> Signed.long -> Signed.long -> _ t
val polchebyshev_eval : Signed.long -> Signed.long -> _ t -> _ t
val polchebyshev1 : Signed.long -> Signed.long -> _ t
val polchebyshev2 : Signed.long -> Signed.long -> _ t
val polrecip : _ t -> _ t
val setbinop : _ t -> _ t -> _ t -> _ t
val setdelta : _ t -> _ t -> _ t
val setintersect : _ t -> _ t -> _ t
val setisset : _ t -> Signed.long
val setminus : _ t -> _ t -> _ t
val setsearch : _ t -> _ t -> Signed.long -> Signed.long
val setunion : _ t -> _ t -> _ t
val setunion_i : _ t -> _ t -> _ t
val sort : _ t -> _ t

val sort_factor :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> int) Ctypes_static.static_funptr ->
  _ t

val stirling : Signed.long -> Signed.long -> Signed.long -> _ t
val stirling1 : pari_ulong -> pari_ulong -> _ t
val stirling2 : pari_ulong -> pari_ulong -> _ t

val tablesearch :
  _ t -> _ t -> (_ t -> _ t -> int) Ctypes_static.static_funptr -> Signed.long

val vecbinomial : Signed.long -> _ t
val vecsearch : _ t -> _ t -> _ t -> Signed.long
val vecsort : _ t -> _ t -> _ t
val vecsort0 : _ t -> _ t -> Signed.long -> _ t
val zv_search : _ t -> Signed.long -> Signed.long
val bits_to_int : _ t -> Signed.long -> _ t
val bits_to_u : _ t -> Signed.long -> pari_ulong
val binaire : _ t -> _ t
val binary_2k : _ t -> Signed.long -> _ t
val binary_2k_nv : _ t -> Signed.long -> _ t
val binary_zv : _ t -> _ t
val bittest : _ t -> Signed.long -> Signed.long
val fromdigits_2k : _ t -> Signed.long -> _ t
val gbitand : _ t -> _ t -> _ t
val gbitneg : _ t -> Signed.long -> _ t
val gbitnegimply : _ t -> _ t -> _ t
val gbitor : _ t -> _ t -> _ t
val gbittest : _ t -> Signed.long -> _ t
val gbitxor : _ t -> _ t -> _ t
val hammingl : pari_ulong -> Signed.long
val hammingweight : _ t -> Signed.long
val ibitand : _ t -> _ t -> _ t
val ibitnegimply : _ t -> _ t -> _ t
val ibitor : _ t -> _ t -> _ t
val ibitxor : _ t -> _ t -> _ t
val nv_fromdigits_2k : _ t -> Signed.long -> _ t
val bnflogef : _ t -> _ t -> _ t
val bnflog : _ t -> _ t -> _ t
val bnflogdegree : _ t -> _ t -> _ t -> _ t
val nfislocalpower : _ t -> _ t -> _ t -> _ t -> Signed.long
val rnfislocalcyclo : _ t -> Signed.long
val bnfisunit : _ t -> _ t -> _ t
val bnfissunit : _ t -> _ t -> _ t -> _ t
val bnfsunit : _ t -> _ t -> Signed.long -> _ t
val bnfunits : _ t -> _ t -> _ t
val bnfisunit0 : _ t -> _ t -> _ t -> _ t
val sunits_mod_units : _ t -> _ t -> _ t
val buchquad : _ t -> float -> float -> Signed.long -> _ t
val quadclassno : _ t -> _ t
val quadclassnos : Signed.long -> Signed.long
val quadclassunit0 : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val buchall : _ t -> Signed.long -> Signed.long -> _ t

val buchall_param :
  _ t -> float -> float -> Signed.long -> Signed.long -> Signed.long -> _ t

val bnf_build_cheapfu : _ t -> _ t
val bnf_build_cycgen : _ t -> _ t
val bnf_build_matalpha : _ t -> _ t
val bnf_build_units : _ t -> _ t
val bnf_compactfu : _ t -> _ t
val bnf_compactfu_mat : _ t -> _ t
val bnf_has_fu : _ t -> _ t
val bnfinit0 : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val bnfisprincipal0 : _ t -> _ t -> Signed.long -> _ t
val bnfnewprec : _ t -> Signed.long -> _ t
val bnfnewprec_shallow : _ t -> Signed.long -> _ t
val bnftestprimes : _ t -> _ t -> unit
val bnrnewprec : _ t -> Signed.long -> _ t
val bnrnewprec_shallow : _ t -> Signed.long -> _ t
val isprincipalfact : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val isprincipalfact_or_fail : _ t -> _ t -> _ t -> _ t -> _ t
val isprincipal : _ t -> _ t -> _ t
val nf_cxlog_normalize : _ t -> _ t -> Signed.long -> _ t
val nfcyclotomicunits : _ t -> _ t -> _ t
val nfsign_units : _ t -> _ t -> int -> _ t
val nfsign_tu : _ t -> _ t -> _ t
val nfsign_fu : _ t -> _ t -> _ t
val signunits : _ t -> _ t
val hermite_bound : Signed.long -> Signed.long -> _ t

val bnr_subgroup_sanitize :
  _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit

val bnr_char_sanitize : _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit
val abc_to_bnr : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> int -> _ t
val buchray : _ t -> _ t -> Signed.long -> _ t
val buchraymod : _ t -> _ t -> Signed.long -> _ t -> _ t
val bnrautmatrix : _ t -> _ t -> _ t
val bnr_subgroup_check : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val bnrchar : _ t -> _ t -> _ t -> _ t
val bnrchar_primitive : _ t -> _ t -> _ t -> _ t
val bnrclassno : _ t -> _ t -> _ t
val bnrclassno0 : _ t -> _ t -> _ t -> _ t
val bnrclassnolist : _ t -> _ t -> _ t
val bnrchar_primitive_raw : _ t -> _ t -> _ t -> _ t
val bnrconductor_factored : _ t -> _ t -> _ t
val bnrconductor_raw : _ t -> _ t -> _ t
val bnrconductormod : _ t -> _ t -> _ t -> _ t
val bnrconductor0 : _ t -> _ t -> _ t -> Signed.long -> _ t
val bnrconductor : _ t -> _ t -> Signed.long -> _ t
val bnrconductor_i : _ t -> _ t -> Signed.long -> _ t
val bnrconductorofchar : _ t -> _ t -> _ t
val bnrdisc0 : _ t -> _ t -> _ t -> Signed.long -> _ t
val bnrdisc : _ t -> _ t -> Signed.long -> _ t
val bnrdisclist0 : _ t -> _ t -> _ t -> _ t
val bnrgaloismatrix : _ t -> _ t -> _ t
val bnrgaloisapply : _ t -> _ t -> _ t -> _ t
val bnrinit0 : _ t -> _ t -> Signed.long -> _ t
val bnrinitmod : _ t -> _ t -> Signed.long -> _ t -> _ t
val bnrisconductor0 : _ t -> _ t -> _ t -> Signed.long
val bnrisconductor : _ t -> _ t -> Signed.long
val bnrisgalois : _ t -> _ t -> _ t -> Signed.long
val bnrisprincipalmod : _ t -> _ t -> _ t -> Signed.long -> _ t
val bnrisprincipal : _ t -> _ t -> Signed.long -> _ t
val bnrmap : _ t -> _ t -> _ t
val bnrsurjection : _ t -> _ t -> _ t
val bnfnarrow : _ t -> _ t
val bnfcertify : _ t -> Signed.long
val bnfcertify0 : _ t -> Signed.long -> Signed.long
val bnrcompositum : _ t -> _ t -> _ t
val decodemodule : _ t -> _ t -> _ t
val discrayabslist : _ t -> _ t -> _ t
val discrayabslistarch : _ t -> _ t -> pari_ulong -> _ t
val idealmoddivisor : _ t -> _ t -> _ t
val isprincipalray : _ t -> _ t -> _ t
val isprincipalraygen : _ t -> _ t -> _ t
val nf_deg1_prime : _ t -> _ t
val nfarchstar : _ t -> _ t -> _ t -> _ t
val rnfconductor : _ t -> _ t -> _ t
val rnfconductor0 : _ t -> _ t -> Signed.long -> _ t
val rnfnormgroup : _ t -> _ t -> _ t
val subgrouplist0 : _ t -> _ t -> Signed.long -> _ t
val bnfisnorm : _ t -> _ t -> Signed.long -> _ t
val rnfisnorm : _ t -> _ t -> Signed.long -> _ t
val rnfisnorminit : _ t -> _ t -> int -> _ t
val coprimes_zv : pari_ulong -> _ t
val char_check : _ t -> _ t -> int
val charconj : _ t -> _ t -> _ t
val charconj0 : _ t -> _ t -> _ t
val chardiv : _ t -> _ t -> _ t -> _ t
val chardiv0 : _ t -> _ t -> _ t -> _ t
val chareval : _ t -> _ t -> _ t -> _ t -> _ t
val chargalois : _ t -> _ t -> _ t
val charker : _ t -> _ t -> _ t
val charker0 : _ t -> _ t -> _ t
val charmul : _ t -> _ t -> _ t -> _ t
val charmul0 : _ t -> _ t -> _ t -> _ t
val charorder : _ t -> _ t -> _ t
val charorder0 : _ t -> _ t -> _ t
val charpow : _ t -> _ t -> _ t -> _ t
val charpow0 : _ t -> _ t -> _ t -> _ t
val char_denormalize : _ t -> _ t -> _ t -> _ t
val char_normalize : _ t -> _ t -> _ t
val char_simplify : _ t -> _ t -> _ t
val checkznstar_i : _ t -> int
val cyc_normalize : _ t -> _ t
val ncharvecexpo : _ t -> _ t -> _ t
val znchar : _ t -> _ t
val znchar_quad : _ t -> _ t -> _ t
val zncharcheck : _ t -> _ t -> int
val zncharconductor : _ t -> _ t -> _ t
val zncharconj : _ t -> _ t -> _ t
val znchardecompose : _ t -> _ t -> _ t -> _ t
val znchardiv : _ t -> _ t -> _ t -> _ t
val znchareval : _ t -> _ t -> _ t -> _ t -> _ t
val zncharinduce : _ t -> _ t -> _ t -> _ t
val zncharisodd : _ t -> _ t -> Signed.long
val zncharker : _ t -> _ t -> _ t
val zncharmul : _ t -> _ t -> _ t -> _ t
val zncharorder : _ t -> _ t -> _ t
val zncharpow : _ t -> _ t -> _ t -> _ t
val znchartokronecker : _ t -> _ t -> Signed.long -> _ t
val znchartoprimitive : _ t -> _ t -> _ t
val znconrey_check : _ t -> _ t -> int
val znconrey_normalized : _ t -> _ t -> _ t
val znconreychar : _ t -> _ t -> _ t
val znconreyfromchar_normalized : _ t -> _ t -> _ t
val znconreyconductor : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val znconreyexp : _ t -> _ t -> _ t
val znconreyfromchar : _ t -> _ t -> _ t
val znconreylog : _ t -> _ t -> _ t
val znconreylog_normalize : _ t -> _ t -> _ t
val znlog0 : _ t -> _ t -> _ t -> _ t
val zv_cyc_minimal : _ t -> _ t -> _ t -> Signed.long
val zv_cyc_minimize : _ t -> _ t -> _ t -> Signed.long
val closure_deriv : _ t -> _ t
val closure_derivn : _ t -> Signed.long -> _ t

val localvars_find :
  _ t -> entree Ctypes.structure Ctypes_static.ptr -> Signed.long

val localvars_read_str : string -> _ t -> _ t
val snm_closure : entree Ctypes.structure Ctypes_static.ptr -> _ t -> _ t
val strtoclosure : string -> Signed.long -> _ t
val strtofunction : string -> _ t
val gconcat : _ t -> _ t -> _ t
val gconcat1 : _ t -> _ t
val matconcat : _ t -> _ t
val shallowconcat : _ t -> _ t -> _ t
val shallowconcat1 : _ t -> _ t
val shallowmatconcat : _ t -> _ t
val vconcat : _ t -> _ t -> _ t
val default0 : string -> string -> _ t
val getrealprecision : unit -> Signed.long
val pari_is_default : string -> entree Ctypes.structure Ctypes_static.ptr
val sd_texstyle : string -> Signed.long -> _ t
val sd_colors : string -> Signed.long -> _ t
val sd_compatible : string -> Signed.long -> _ t
val sd_datadir : string -> Signed.long -> _ t
val sd_debug : string -> Signed.long -> _ t
val sd_debugfiles : string -> Signed.long -> _ t
val sd_debugmem : string -> Signed.long -> _ t
val sd_factor_add_primes : string -> Signed.long -> _ t
val sd_factor_proven : string -> Signed.long -> _ t
val sd_format : string -> Signed.long -> _ t
val sd_histsize : string -> Signed.long -> _ t
val sd_log : string -> Signed.long -> _ t
val sd_logfile : string -> Signed.long -> _ t
val sd_nbthreads : string -> Signed.long -> _ t
val sd_new_galois_format : string -> Signed.long -> _ t
val sd_output : string -> Signed.long -> _ t
val sd_parisize : string -> Signed.long -> _ t
val sd_parisizemax : string -> Signed.long -> _ t
val sd_path : string -> Signed.long -> _ t
val sd_plothsizes : string -> Signed.long -> _ t
val sd_prettyprinter : string -> Signed.long -> _ t
val sd_primelimit : string -> Signed.long -> _ t
val sd_realbitprecision : string -> Signed.long -> _ t
val sd_realprecision : string -> Signed.long -> _ t
val sd_secure : string -> Signed.long -> _ t
val sd_seriesprecision : string -> Signed.long -> _ t
val sd_simplify : string -> Signed.long -> _ t
val sd_sopath : string -> int -> _ t
val sd_strictargs : string -> Signed.long -> _ t
val sd_strictmatch : string -> Signed.long -> _ t

val sd_string :
  string -> Signed.long -> string -> string Ctypes_static.ptr -> _ t

val sd_threadsize : string -> Signed.long -> _ t
val sd_threadsizemax : string -> Signed.long -> _ t

val sd_intarray :
  string -> Signed.long -> _ t Ctypes_static.ptr -> string -> _ t

val sd_toggle : string -> Signed.long -> string -> int Ctypes_static.ptr -> _ t

val sd_ulong :
  string ->
  Signed.long ->
  string ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  pari_ulong ->
  string Ctypes_static.ptr ->
  _ t

val setdefault : string -> string -> Signed.long -> _ t

val setrealprecision :
  Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val digits : _ t -> _ t -> _ t
val fromdigits : _ t -> _ t -> _ t
val fromdigitsu : _ t -> _ t -> _ t

val gen_digits :
  _ t ->
  _ t ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t)
  Ctypes_static.static_funptr ->
  _ t

val gen_fromdigits :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  _ t

val sumdigits : _ t -> _ t
val sumdigits0 : _ t -> _ t -> _ t
val sumdigitsu : pari_ulong -> pari_ulong
val ecpp : _ t -> _ t
val ecpp0 : _ t -> Signed.long -> _ t
val ecppexport : _ t -> Signed.long -> _ t
val ecppisvalid : _ t -> Signed.long
val isprimeecpp : _ t -> Signed.long
val sd_breakloop : string -> Signed.long -> _ t
val sd_echo : string -> Signed.long -> _ t
val sd_graphcolormap : string -> Signed.long -> _ t
val sd_graphcolors : string -> Signed.long -> _ t
val sd_help : string -> Signed.long -> _ t
val sd_histfile : string -> Signed.long -> _ t
val sd_lines : string -> Signed.long -> _ t
val sd_linewrap : string -> Signed.long -> _ t
val sd_prompt : string -> Signed.long -> _ t
val sd_prompt_cont : string -> Signed.long -> _ t
val sd_psfile : string -> Signed.long -> _ t
val sd_readline : string -> Signed.long -> _ t
val sd_recover : string -> Signed.long -> _ t
val sd_timer : string -> Signed.long -> _ t
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
val gp_alarm : Signed.long -> _ t -> _ t
val gp_input : unit -> _ t
val gp_allocatemem : _ t -> unit
val gp_handle_exception : Signed.long -> int
val gp_alarm_handler : int -> unit
val gp_sigint_fun : unit -> unit
val gp_help : string -> Signed.long -> unit
val gp_echo_and_log : string -> string -> unit
val print_fun_list : string Ctypes_static.ptr -> Signed.long -> unit
val strtime : Signed.long -> _ t

val direuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  _ t ->
  _ t

val dirpowers : Signed.long -> _ t -> Signed.long -> _ t
val dirpowerssum : pari_ulong -> _ t -> Signed.long -> Signed.long -> _ t

val dirpowerssumfun :
  pari_ulong ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> pari_ulong -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  _ t

val vecpowuu : Signed.long -> pari_ulong -> _ t
val vecpowug : Signed.long -> _ t -> Signed.long -> _ t
val ellanalyticrank : _ t -> _ t -> Signed.long -> _ t
val ellanalyticrank_bitprec : _ t -> _ t -> Signed.long -> _ t

val ellanal_globalred_all :
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val ellheegner : _ t -> _ t
val elll1 : _ t -> Signed.long -> Signed.long -> _ t
val elll1_bitprec : _ t -> Signed.long -> Signed.long -> _ t
val ellconvertname : _ t -> _ t
val elldatagenerators : _ t -> _ t
val ellidentify : _ t -> _ t
val ellsearch : _ t -> _ t
val ellsearchcurve : _ t -> _ t

val forell :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  unit

val ellfromeqn : _ t -> _ t
val akell : _ t -> _ t -> _ t
val bilhell : _ t -> _ t -> _ t -> Signed.long -> _ t
val checkell : _ t -> unit
val checkell_fq : _ t -> unit
val checkell_q : _ t -> unit
val checkell_qp : _ t -> unit
val checkellisog : _ t -> unit
val checkellpt : _ t -> unit
val checkell5 : _ t -> unit
val cxredsl2 : _ t -> _ t Ctypes_static.ptr -> _ t
val cxredsl2_i : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t
val ec_2divpol_evalx : _ t -> _ t -> _ t
val ec_3divpol_evalx : _ t -> _ t -> _ t
val ec_bmodel : _ t -> Signed.long -> _ t
val ec_f_evalx : _ t -> _ t -> _ t
val ec_h_evalx : _ t -> _ t -> _ t
val ec_dfdx_evalq : _ t -> _ t -> _ t
val ec_dfdy_evalq : _ t -> _ t -> _ t
val ec_dmfdy_evalq : _ t -> _ t -> _ t
val ec_half_deriv_2divpol : _ t -> Signed.long -> _ t
val ec_half_deriv_2divpol_evalx : _ t -> _ t -> _ t
val ec_phi2 : _ t -> Signed.long -> _ t
val ell_is_integral : _ t -> int
val ellq_get_cm : _ t -> Signed.long
val ellq_get_n : _ t -> _ t
val ellq_get_nfa : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit
val ellqp_tate_uniformization : _ t -> Signed.long -> _ t
val ellqp_agm : _ t -> Signed.long -> _ t
val ellqp_u : _ t -> Signed.long -> _ t
val ellqp_u2 : _ t -> Signed.long -> _ t
val ellqp_q : _ t -> Signed.long -> _ t
val ellqp_ab : _ t -> Signed.long -> _ t
val ellqp_l : _ t -> Signed.long -> _ t
val ellqp_root : _ t -> Signed.long -> _ t
val ellqtwist_bsdperiod : _ t -> Signed.long -> _ t
val ellr_area : _ t -> Signed.long -> _ t
val ellr_ab : _ t -> Signed.long -> _ t
val ellr_eta : _ t -> Signed.long -> _ t
val ellr_omega : _ t -> Signed.long -> _ t
val ellr_roots : _ t -> Signed.long -> _ t
val elladd : _ t -> _ t -> _ t -> _ t
val ellan : _ t -> Signed.long -> _ t
val ellanq_zv : _ t -> Signed.long -> _ t
val ellanal_globalred : _ t -> _ t Ctypes_static.ptr -> _ t
val ellap : _ t -> _ t -> _ t
val ellap_cm_fast : _ t -> pari_ulong -> Signed.long -> Signed.long
val ellbasechar : _ t -> _ t
val ellbsd : _ t -> Signed.long -> _ t
val ellcard : _ t -> _ t -> _ t
val ellchangecurve : _ t -> _ t -> _ t
val ellchangeinvert : _ t -> _ t
val ellchangepoint : _ t -> _ t -> _ t
val ellchangepointinv : _ t -> _ t -> _ t
val elldivpol : _ t -> Signed.long -> Signed.long -> _ t
val elleisnum : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val elleta : _ t -> Signed.long -> _ t
val elleulerf : _ t -> _ t -> _ t
val ellff_get_card : _ t -> _ t
val ellff_get_gens : _ t -> _ t
val ellff_get_group : _ t -> _ t
val ellff_get_o : _ t -> _ t
val ellff_get_p : _ t -> _ t
val ellff_get_m : _ t -> _ t
val ellff_get_d : _ t -> _ t
val ellfromj : _ t -> _ t
val ellgenerators : _ t -> _ t
val ellglobalred : _ t -> _ t
val ellgroup : _ t -> _ t -> _ t
val ellgroup0 : _ t -> _ t -> Signed.long -> _ t
val ellheight0 : _ t -> _ t -> _ t -> Signed.long -> _ t
val ellheight : _ t -> _ t -> Signed.long -> _ t
val ellheightmatrix : _ t -> _ t -> Signed.long -> _ t
val ellheightoo : _ t -> _ t -> Signed.long -> _ t
val ellinit : _ t -> _ t -> Signed.long -> _ t
val ellintegralmodel : _ t -> _ t Ctypes_static.ptr -> _ t
val ellintegralmodel_i : _ t -> _ t Ctypes_static.ptr -> _ t
val elliscm : _ t -> Signed.long
val ellisoncurve : _ t -> _ t -> _ t
val ellisotree : _ t -> _ t
val ellissupersingular : _ t -> _ t -> int
val elljissupersingular : _ t -> int
val elllseries : _ t -> _ t -> _ t -> Signed.long -> _ t
val elllocalred : _ t -> _ t -> _ t
val elllog : _ t -> _ t -> _ t -> _ t -> _ t
val ellminimaldisc : _ t -> _ t
val ellminimalmodel : _ t -> _ t Ctypes_static.ptr -> _ t
val ellminimaltwist : _ t -> _ t
val ellminimaltwist0 : _ t -> Signed.long -> _ t
val ellminimaltwistcond : _ t -> _ t
val ellmul : _ t -> _ t -> _ t -> _ t
val ellnf_vecarea : _ t -> Signed.long -> _ t
val ellnf_veceta : _ t -> Signed.long -> _ t
val ellnf_vecomega : _ t -> Signed.long -> _ t
val ellneg : _ t -> _ t -> _ t
val ellorder : _ t -> _ t -> _ t -> _ t
val ellorder_q : _ t -> _ t -> Signed.long
val ellordinate : _ t -> _ t -> Signed.long -> _ t
val ellpadicheight0 : _ t -> _ t -> Signed.long -> _ t -> _ t -> _ t
val ellpadicheightmatrix : _ t -> _ t -> Signed.long -> _ t -> _ t
val ellperiods : _ t -> Signed.long -> Signed.long -> _ t
val ellrandom : _ t -> _ t
val ellrootno : _ t -> _ t -> Signed.long
val ellrootno_global : _ t -> Signed.long
val ellsaturation : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val ellsea : _ t -> Signed.long -> _ t
val ellsigma : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val ellsub : _ t -> _ t -> _ t -> _ t
val ellsupersingularj : _ t -> _ t
val elltamagawa : _ t -> _ t
val elltaniyama : _ t -> Signed.long -> _ t
val elltatepairing : _ t -> _ t -> _ t -> _ t -> _ t
val elltors : _ t -> _ t
val elltors0 : _ t -> Signed.long -> _ t
val elltors_psylow : _ t -> pari_ulong -> _ t
val elltrace : _ t -> _ t -> _ t
val elltwist : _ t -> _ t -> _ t
val ellweilpairing : _ t -> _ t -> _ t -> _ t -> _ t
val ellwp : _ t -> _ t -> Signed.long -> _ t
val ellwp0 : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val ellwpseries : _ t -> Signed.long -> Signed.long -> _ t
val ellxn : _ t -> Signed.long -> Signed.long -> _ t
val ellzeta : _ t -> _ t -> Signed.long -> _ t
val oncurve : _ t -> _ t -> int
val orderell : _ t -> _ t -> _ t
val pointell : _ t -> _ t -> Signed.long -> _ t
val point_to_a4a6 : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t

val point_to_a4a6_fl :
  _ t -> _ t -> pari_ulong -> pari_ulong Ctypes_static.ptr -> _ t

val zell : _ t -> _ t -> Signed.long -> _ t
val qp_agm2_sequence : _ t -> _ t -> _ t

val qp_ascending_landen :
  _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit

val qp_descending_landen :
  _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit

val ellformaldifferential : _ t -> Signed.long -> Signed.long -> _ t
val ellformalexp : _ t -> Signed.long -> Signed.long -> _ t
val ellformallog : _ t -> Signed.long -> Signed.long -> _ t
val ellformalpoint : _ t -> Signed.long -> Signed.long -> _ t
val ellformalw : _ t -> Signed.long -> Signed.long -> _ t
val ellnonsingularmultiple : _ t -> _ t -> _ t
val ellpadicl : _ t -> _ t -> Signed.long -> _ t -> Signed.long -> _ t -> _ t
val ellpadicbsd : _ t -> _ t -> Signed.long -> _ t -> _ t
val ellpadicfrobenius : _ t -> pari_ulong -> Signed.long -> _ t
val ellpadicheight : _ t -> _ t -> Signed.long -> _ t -> _ t
val ellpadiclog : _ t -> _ t -> Signed.long -> _ t -> _ t
val ellpadicregulator : _ t -> _ t -> Signed.long -> _ t -> _ t
val ellpadics2 : _ t -> _ t -> Signed.long -> _ t
val ell2cover : _ t -> Signed.long -> _ t
val ellrank : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val ellrankinit : _ t -> Signed.long -> _ t
val hyperell_locally_soluble : _ t -> _ t -> Signed.long
val nf_hyperell_locally_soluble : _ t -> _ t -> _ t -> Signed.long
val nfhilbert : _ t -> _ t -> _ t -> Signed.long
val nfhilbert0 : _ t -> _ t -> _ t -> _ t -> Signed.long
val ellisdivisible : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val ellisogenyapply : _ t -> _ t -> _ t
val ellisogeny : _ t -> _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val ellisomat : _ t -> Signed.long -> Signed.long -> _ t
val ellweilcurve : _ t -> _ t Ctypes_static.ptr -> _ t

val flxq_elldivpolmod :
  _ t -> _ t -> Signed.long -> _ t -> _ t -> pari_ulong -> _ t

val fp_ellcard_sea : _ t -> _ t -> _ t -> Signed.long -> _ t
val fq_ellcard_sea : _ t -> _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val fq_elldivpolmod : _ t -> _ t -> Signed.long -> _ t -> _ t -> _ t -> _ t
val ellmodulareqn : Signed.long -> Signed.long -> Signed.long -> _ t
val externstr : string -> _ t
val gp_filter : string -> string
val gpextern : string -> _ t
val gpsystem : string -> Signed.long
val readstr : string -> _ t
val gentogenstr_nospace : _ t -> _ t
val gentogenstr : _ t -> _ t
val gentotexstr : _ t -> string
val gentostr : _ t -> string
val gentostr_raw : _ t -> string
val gentostr_unquoted : _ t -> string
val str : _ t -> _ t
val strexpand : _ t -> _ t
val strtex : _ t -> _ t
val brute : _ t -> char -> Signed.long -> unit
val dbggen : _ t -> Signed.long -> unit
val error0 : _ t -> unit
val dbg_pari_heap : unit -> unit
val err_flush : unit -> unit
val err_printf : string -> unit
val gp_getenv : string -> _ t
val gp_fileclose : Signed.long -> unit
val gp_fileextern : string -> Signed.long
val gp_fileflush : Signed.long -> unit
val gp_fileflush0 : _ t -> unit
val gp_fileopen : string -> string -> Signed.long
val gp_fileread : Signed.long -> _ t
val gp_filereadstr : Signed.long -> _ t
val gp_filewrite : Signed.long -> string -> unit
val gp_filewrite1 : Signed.long -> string -> unit
val gp_read_file : string -> _ t
val gp_read_str_multiline : string -> string -> _ t
val gp_readvec_file : string -> _ t
val gpinstall : string -> string -> string -> string -> unit
val gsprintf : string -> _ t
val itostr : _ t -> string
val matbrute : _ t -> char -> Signed.long -> unit
val os_getenv : string -> string
val uordinal : pari_ulong -> string
val outmat : _ t -> unit
val output : _ t -> unit
val rgv_to_str : _ t -> Signed.long -> string
val pari_add_hist : _ t -> Signed.long -> Signed.long -> unit
val pari_ask_confirm : string -> unit
val pari_flush : unit -> unit
val pari_get_hist : Signed.long -> _ t
val pari_get_histrtime : Signed.long -> Signed.long
val pari_get_histtime : Signed.long -> Signed.long
val pari_get_homedir : string -> string
val pari_histtime : Signed.long -> _ t
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
val pari_sprint0 : string -> _ t -> Signed.long -> string
val print : _ t -> unit
val printp : _ t -> unit
val print1 : _ t -> unit
val printf0 : string -> _ t -> unit
val printsep : string -> _ t -> unit
val printsep1 : string -> _ t -> unit
val printtex : _ t -> unit
val stack_sprintf : string -> string
val str_init : pari_str Ctypes.structure Ctypes_static.ptr -> int -> unit
val str_printf : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val str_putc : pari_str Ctypes.structure Ctypes_static.ptr -> char -> unit
val str_puts : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val strftime_expand : string -> string -> Signed.long -> unit
val strprintf : string -> _ t -> _ t
val term_color : Signed.long -> unit
val term_get_color : string -> Signed.long -> string
val texe : _ t -> char -> Signed.long -> unit
val warning0 : _ t -> unit
val write0 : string -> _ t -> unit
val write1 : string -> _ t -> unit
val writebin : string -> _ t -> unit
val writetex : string -> _ t -> unit
val bincopy_relink : _ t -> _ t -> unit
val bitprecision0 : _ t -> Signed.long -> _ t
val bitprecision00 : _ t -> _ t -> _ t
val break0 : Signed.long -> _ t
val call0 : _ t -> _ t -> _ t
val closure_callgen0prec : _ t -> Signed.long -> _ t
val closure_callgen1 : _ t -> _ t -> _ t
val closure_callgen1prec : _ t -> _ t -> Signed.long -> _ t
val closure_callgen2 : _ t -> _ t -> _ t -> _ t
val closure_callgenall : _ t -> Signed.long -> _ t
val closure_callgenvec : _ t -> _ t -> _ t
val closure_callgenvecdef : _ t -> _ t -> _ t -> _ t
val closure_callgenvecdefprec : _ t -> _ t -> _ t -> Signed.long -> _ t
val closure_callgenvecprec : _ t -> _ t -> Signed.long -> _ t
val closure_callvoid1 : _ t -> _ t -> unit
val closure_context : Signed.long -> Signed.long -> Signed.long
val closure_disassemble : _ t -> unit
val closure_err : Signed.long -> unit
val closure_evalbrk : _ t -> Signed.long Ctypes_static.ptr -> _ t
val closure_evalgen : _ t -> _ t
val closure_evalnobrk : _ t -> _ t
val closure_evalres : _ t -> _ t
val closure_evalvoid : _ t -> unit
val closure_func_err : unit -> string
val closure_trapgen : _ t -> Signed.long -> _ t
val copybin_unlink : _ t -> _ t
val getlocalprec : Signed.long -> Signed.long
val getlocalbitprec : Signed.long -> Signed.long
val get_lex : Signed.long -> _ t
val get_localprec : unit -> Signed.long
val get_localbitprec : unit -> Signed.long
val gp_call : unit Ctypes_static.ptr -> _ t -> _ t
val gp_callprec : unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t
val gp_call2 : unit Ctypes_static.ptr -> _ t -> _ t -> _ t
val gp_callbool : unit Ctypes_static.ptr -> _ t -> Signed.long
val gp_callvoid : unit Ctypes_static.ptr -> _ t -> Signed.long
val gp_eval : unit Ctypes_static.ptr -> _ t -> _ t
val gp_evalbool : unit Ctypes_static.ptr -> _ t -> Signed.long
val gp_evalprec : unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t
val gp_evalupto : unit Ctypes_static.ptr -> _ t -> _ t
val gp_evalvoid : unit Ctypes_static.ptr -> _ t -> Signed.long
val localprec : _ t -> unit
val localbitprec : _ t -> unit
val loop_break : unit -> Signed.long
val next0 : Signed.long -> _ t
val pareval : _ t -> _ t
val pari_self : unit -> _ t
val parsum : _ t -> _ t -> _ t -> _ t
val parvector : Signed.long -> _ t -> _ t
val pop_lex : Signed.long -> unit
val pop_localprec : unit -> unit
val precision0 : _ t -> Signed.long -> _ t
val precision00 : _ t -> _ t -> _ t
val push_lex : _ t -> _ t -> unit
val push_localbitprec : Signed.long -> unit
val push_localprec : Signed.long -> unit
val return0 : _ t -> _ t
val set_lex : Signed.long -> _ t -> unit

val forcomposite_init :
  forcomposite_t Ctypes.structure Ctypes_static.ptr -> _ t -> _ t -> int

val forcomposite_next : forcomposite_t Ctypes.structure Ctypes_static.ptr -> _ t
val forprime_next : forprime_t Ctypes.structure Ctypes_static.ptr -> _ t

val forprime_init :
  forprime_t Ctypes.structure Ctypes_static.ptr -> _ t -> _ t -> int

val forprimestep_init :
  forprime_t Ctypes.structure Ctypes_static.ptr -> _ t -> _ t -> _ t -> int

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
val prodprimes : unit -> _ t

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

val ff_1 : _ t -> _ t
val ff_frobenius : _ t -> Signed.long -> _ t
val ff_z_z_muldiv : _ t -> _ t -> _ t -> _ t
val ff_q_add : _ t -> _ t -> _ t
val ff_z_add : _ t -> _ t -> _ t
val ff_z_mul : _ t -> _ t -> _ t
val ff_add : _ t -> _ t -> _ t
val ff_charpoly : _ t -> _ t
val ff_conjvec : _ t -> _ t
val ff_div : _ t -> _ t -> _ t
val ff_ellcard : _ t -> _ t
val ff_ellcard_sea : _ t -> Signed.long -> _ t
val ff_ellgens : _ t -> _ t
val ff_ellgroup : _ t -> _ t Ctypes_static.ptr -> _ t
val ff_elllog : _ t -> _ t -> _ t -> _ t -> _ t
val ff_ellmul : _ t -> _ t -> _ t -> _ t
val ff_ellorder : _ t -> _ t -> _ t -> _ t
val ff_elltwist : _ t -> _ t
val ff_ellrandom : _ t -> _ t
val ff_elltatepairing : _ t -> _ t -> _ t -> _ t -> _ t
val ff_ellweilpairing : _ t -> _ t -> _ t -> _ t -> _ t
val ff_equal : _ t -> _ t -> int
val ff_equal0 : _ t -> int
val ff_equal1 : _ t -> int
val ff_equalm1 : _ t -> int
val ff_f : _ t -> Signed.long
val ff_gen : _ t -> _ t
val ff_inv : _ t -> _ t
val ff_issquare : _ t -> Signed.long
val ff_issquareall : _ t -> _ t Ctypes_static.ptr -> Signed.long
val ff_ispower : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val ff_log : _ t -> _ t -> _ t -> _ t
val ff_map : _ t -> _ t -> _ t
val ff_minpoly : _ t -> _ t
val ff_mod : _ t -> _ t
val ff_mul : _ t -> _ t -> _ t
val ff_mul2n : _ t -> Signed.long -> _ t
val ff_neg : _ t -> _ t
val ff_neg_i : _ t -> _ t
val ff_norm : _ t -> _ t
val ff_order : _ t -> _ t -> _ t
val ff_p : _ t -> _ t
val ff_p_i : _ t -> _ t
val ff_pow : _ t -> _ t -> _ t
val ff_primroot : _ t -> _ t Ctypes_static.ptr -> _ t
val ff_q : _ t -> _ t
val ff_samefield : _ t -> _ t -> int
val ff_sqr : _ t -> _ t
val ff_sqrt : _ t -> _ t
val ff_sqrtn : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val ff_sub : _ t -> _ t -> _ t
val ff_to_f2xq : _ t -> _ t
val ff_to_f2xq_i : _ t -> _ t
val ff_to_flxq : _ t -> _ t
val ff_to_flxq_i : _ t -> _ t
val ff_to_fpxq : _ t -> _ t
val ff_to_fpxq_i : _ t -> _ t
val ff_trace : _ t -> _ t
val ff_var : _ t -> Signed.long
val ff_zero : _ t -> _ t
val ffm_ffc_invimage : _ t -> _ t -> _ t -> _ t
val ffm_ffc_gauss : _ t -> _ t -> _ t -> _ t
val ffm_ffc_mul : _ t -> _ t -> _ t -> _ t
val ffm_deplin : _ t -> _ t -> _ t
val ffm_det : _ t -> _ t -> _ t
val ffm_gauss : _ t -> _ t -> _ t -> _ t
val ffm_image : _ t -> _ t -> _ t
val ffm_indexrank : _ t -> _ t -> _ t
val ffm_inv : _ t -> _ t -> _ t
val ffm_invimage : _ t -> _ t -> _ t -> _ t
val ffm_ker : _ t -> _ t -> _ t
val ffm_mul : _ t -> _ t -> _ t -> _ t
val ffm_rank : _ t -> _ t -> Signed.long
val ffm_suppl : _ t -> _ t -> _ t
val ffx_add : _ t -> _ t -> _ t -> _ t
val ffx_ddf : _ t -> _ t -> _ t
val ffx_degfact : _ t -> _ t -> _ t
val ffx_disc : _ t -> _ t -> _ t

val ffx_extgcd :
  _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val ffx_factor : _ t -> _ t -> _ t
val ffx_factor_squarefree : _ t -> _ t -> _ t
val ffx_gcd : _ t -> _ t -> _ t -> _ t
val ffx_halfgcd : _ t -> _ t -> _ t -> _ t

val ffx_halfgcd_all :
  _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val ffx_ispower :
  _ t -> Signed.long -> _ t -> _ t Ctypes_static.ptr -> Signed.long

val ffx_mul : _ t -> _ t -> _ t -> _ t
val ffx_preimage : _ t -> _ t -> _ t -> _ t
val ffx_preimagerel : _ t -> _ t -> _ t -> _ t
val ffx_rem : _ t -> _ t -> _ t -> _ t
val ffx_resultant : _ t -> _ t -> _ t -> _ t
val ffx_roots : _ t -> _ t -> _ t
val ffx_sqr : _ t -> _ t -> _ t
val ffxq_inv : _ t -> _ t -> _ t -> _ t
val ffxq_minpoly : _ t -> _ t -> _ t -> _ t
val ffxq_mul : _ t -> _ t -> _ t -> _ t -> _ t
val ffxq_sqr : _ t -> _ t -> _ t -> _ t
val fqx_to_ffx : _ t -> _ t -> _ t
val fq_to_ff : _ t -> _ t -> _ t
val z_ff_div : _ t -> _ t -> _ t
val ffembed : _ t -> _ t -> _ t
val ffextend : _ t -> _ t -> Signed.long -> _ t
val fffrobenius : _ t -> Signed.long -> _ t
val ffgen : _ t -> Signed.long -> _ t
val ffinvmap : _ t -> _ t
val fflog : _ t -> _ t -> _ t -> _ t
val ffmap : _ t -> _ t -> _ t
val ffmaprel : _ t -> _ t -> _ t
val ffcompomap : _ t -> _ t -> _ t
val fforder : _ t -> _ t -> _ t
val ffprimroot : _ t -> _ t Ctypes_static.ptr -> _ t
val ffrandom : _ t -> _ t
val rg_is_ff : _ t -> _ t Ctypes_static.ptr -> int
val rgc_is_ffc : _ t -> _ t Ctypes_static.ptr -> int
val rgm_is_ffm : _ t -> _ t Ctypes_static.ptr -> int
val p_to_ff : _ t -> Signed.long -> _ t
val tp_to_ff : _ t -> _ t -> _ t
val flx_factcyclo : pari_ulong -> pari_ulong -> pari_ulong -> _ t
val fpx_factcyclo : pari_ulong -> _ t -> pari_ulong -> _ t
val factormodcyclo : Signed.long -> _ t -> Signed.long -> Signed.long -> _ t
val checkgal : _ t -> _ t
val checkgroup : _ t -> _ t Ctypes_static.ptr -> _ t
val checkgroupelts : _ t -> _ t
val embed_disc : _ t -> Signed.long -> Signed.long -> _ t
val embed_roots : _ t -> Signed.long -> _ t
val galois_group : _ t -> _ t
val galoisconj : _ t -> _ t -> _ t
val galoisconj0 : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val galoisconjclasses : _ t -> _ t
val galoisexport : _ t -> Signed.long -> _ t
val galoisfixedfield : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val galoisidentify : _ t -> _ t
val galoisinit : _ t -> _ t -> _ t
val galoisisabelian : _ t -> Signed.long -> _ t
val galoisisnormal : _ t -> _ t -> Signed.long
val galoispermtopol : _ t -> _ t -> _ t
val galoissplittinginit : _ t -> _ t -> _ t
val galoissubgroups : _ t -> _ t
val galoissubfields : _ t -> Signed.long -> Signed.long -> _ t
val numberofconjugates : _ t -> Signed.long -> Signed.long
val polgalois : _ t -> Signed.long -> _ t
val galoisnbpol : Signed.long -> _ t
val galoisgetgroup : Signed.long -> Signed.long -> _ t
val galoisgetname : Signed.long -> Signed.long -> _ t
val galoisgetpol : Signed.long -> Signed.long -> Signed.long -> _ t
val conj_i : _ t -> _ t
val conjvec : _ t -> Signed.long -> _ t
val divrunextu : _ t -> pari_ulong -> _ t
val gadd : _ t -> _ t -> _ t
val gaddsg : Signed.long -> _ t -> _ t
val gconj : _ t -> _ t
val gdiv : _ t -> _ t -> _ t
val gdivgs : _ t -> Signed.long -> _ t
val gdivgu : _ t -> pari_ulong -> _ t
val gdivgunextu : _ t -> pari_ulong -> _ t
val ginv : 'a t -> 'a t
val gmul : 'a t -> 'a t -> 'a t
val gmul2n : _ t -> Signed.long -> _ t
val gmulsg : Signed.long -> _ t -> _ t
val gmulug : pari_ulong -> _ t -> _ t
val gsqr : 'a t -> 'a t
val gsub : 'a t -> 'a t -> 'a t
val gsubsg : Signed.long -> _ t -> _ t
val mulcxi : _ t -> _ t
val mulcxmi : _ t -> _ t
val mulcxpowis : _ t -> Signed.long -> _ t
val qdivii : _ t -> _ t -> _ t
val qdiviu : _ t -> pari_ulong -> _ t
val qdivis : _ t -> Signed.long -> _ t
val ser_normalize : _ t -> _ t

val gassoc_proto :
  (_ t -> _ t -> _ t) Ctypes_static.static_funptr -> _ t -> _ t -> _ t

val map_proto_g : (_ t -> _ t) Ctypes_static.static_funptr -> _ t -> _ t

val map_proto_lg :
  (_ t -> Signed.long) Ctypes_static.static_funptr -> _ t -> _ t

val map_proto_lgl :
  (_ t -> Signed.long -> Signed.long) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val q_lval : _ t -> pari_ulong -> Signed.long
val q_lvalrem : _ t -> pari_ulong -> _ t Ctypes_static.ptr -> Signed.long
val q_pval : _ t -> _ t -> Signed.long
val q_pvalrem : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val rgx_val : _ t -> Signed.long
val rgx_valrem : _ t -> _ t Ctypes_static.ptr -> Signed.long
val rgx_valrem_inexact : _ t -> _ t Ctypes_static.ptr -> Signed.long
val rgxv_maxdegree : _ t -> Signed.long
val zv_z_dvd : _ t -> _ t -> int
val zv_pval : _ t -> _ t -> Signed.long
val zv_pvalrem : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val zv_lval : _ t -> pari_ulong -> Signed.long
val zv_lvalrem : _ t -> pari_ulong -> _ t Ctypes_static.ptr -> Signed.long
val zx_lvalrem : _ t -> pari_ulong -> _ t Ctypes_static.ptr -> Signed.long
val zx_pval : _ t -> _ t -> Signed.long
val zx_pvalrem : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long

val z_lvalrem_stop :
  _ t Ctypes_static.ptr -> pari_ulong -> int Ctypes_static.ptr -> Signed.long

val cgetp : _ t -> _ t
val cvstop2 : Signed.long -> _ t -> _ t
val cvtop : _ t -> _ t -> Signed.long -> _ t
val cvtop2 : _ t -> _ t -> _ t
val cx_approx_equal : _ t -> _ t -> int
val cx_approx0 : _ t -> _ t -> int
val gabs : _ t -> Signed.long -> _ t
val gaffect : _ t -> _ t -> unit
val gaffsg : Signed.long -> _ t -> unit
val gcmp : _ t -> _ t -> int
val gequal0 : _ t -> int
val gequal1 : _ t -> int
val gequalx : _ t -> int
val gequalm1 : _ t -> int
val gcmpsg : Signed.long -> _ t -> int
val gcvtop : _ t -> _ t -> Signed.long -> _ t
val gequal : _ t -> _ t -> int
val gequalsg : Signed.long -> _ t -> int
val gexpo : _ t -> Signed.long
val gexpo_safe : _ t -> Signed.long
val gpexponent : _ t -> _ t
val gpvaluation : _ t -> _ t -> _ t
val gvaluation : _ t -> _ t -> Signed.long
val gidentical : _ t -> _ t -> int
val glength : _ t -> Signed.long
val gmax : _ t -> _ t -> _ t
val gmaxgs : _ t -> Signed.long -> _ t
val gmin : _ t -> _ t -> _ t
val gmings : _ t -> Signed.long -> _ t
val gneg : _ t -> _ t
val gneg_i : [ `int ] t -> [ `int ] t
val gsigne : _ t -> int
val gtolist : _ t -> _ t
val gtolong : _ t -> Signed.long
val lexcmp : _ t -> _ t -> int
val listinsert : _ t -> _ t -> Signed.long -> _ t
val listpop : _ t -> Signed.long -> unit
val listpop0 : _ t -> Signed.long -> unit
val listput : _ t -> _ t -> Signed.long -> _ t
val listput0 : _ t -> _ t -> Signed.long -> unit
val listsort : _ t -> Signed.long -> unit
val matsize : _ t -> _ t
val mklist : unit -> _ t
val mklist_typ : Signed.long -> _ t
val mklistcopy : _ t -> _ t
val mkmap : unit -> _ t
val normalizeser : _ t -> _ t
val normalizepol : _ t -> _ t
val normalizepol_approx : _ t -> Signed.long -> _ t
val normalizepol_lg : _ t -> Signed.long -> _ t
val padic_to_fl : _ t -> pari_ulong -> pari_ulong
val padic_to_fp : _ t -> _ t -> _ t
val quadtofp : _ t -> Signed.long -> _ t
val sizedigit : _ t -> Signed.long
val u_lval : pari_ulong -> pari_ulong -> Signed.long

val u_lvalrem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val u_lvalrem_stop :
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  int Ctypes_static.ptr ->
  Signed.long

val u_pval : pari_ulong -> _ t -> Signed.long
val u_pvalrem : pari_ulong -> _ t -> pari_ulong Ctypes_static.ptr -> Signed.long
val vecindexmax : _ t -> Signed.long
val vecindexmin : _ t -> Signed.long
val vecmax0 : _ t -> _ t Ctypes_static.ptr -> _ t
val vecmax : _ t -> _ t
val vecmin0 : _ t -> _ t Ctypes_static.ptr -> _ t
val vecmin : _ t -> _ t
val z_lval : Signed.long -> pari_ulong -> Signed.long

val z_lvalrem :
  Signed.long -> pari_ulong -> Signed.long Ctypes_static.ptr -> Signed.long

val z_pval : Signed.long -> _ t -> Signed.long

val z_pvalrem :
  Signed.long -> _ t -> Signed.long Ctypes_static.ptr -> Signed.long

val zx_lval : _ t -> Signed.long -> Signed.long
val hgmcyclo : _ t -> _ t
val hgmalpha : _ t -> _ t
val hgmgamma : _ t -> _ t
val hgminit : _ t -> _ t -> _ t
val hgmparams : _ t -> _ t
val hgmeulerfactor : _ t -> _ t -> Signed.long -> _ t Ctypes_static.ptr -> _ t
val hgmcoef : _ t -> _ t -> _ t -> _ t
val hgmcoefs : _ t -> _ t -> Signed.long -> _ t
val hgmtwist : _ t -> _ t
val hgmissymmetrical : _ t -> Signed.long
val hgmbydegree : Signed.long -> _ t
val lfunhgm : _ t -> _ t -> _ t -> Signed.long -> _ t
val qp_zeta : _ t -> _ t
val lerchphi : _ t -> _ t -> _ t -> Signed.long -> _ t
val lerchzeta : _ t -> _ t -> _ t -> Signed.long -> _ t
val zetahurwitz : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val rgx_to_ser : _ t -> Signed.long -> _ t
val rgx_to_ser_inexact : _ t -> Signed.long -> _ t
val gtoser : _ t -> Signed.long -> Signed.long -> _ t
val gtoser_prec : _ t -> Signed.long -> Signed.long -> _ t
val rfrac_to_ser : _ t -> Signed.long -> _ t
val rfrac_to_ser_i : _ t -> Signed.long -> _ t
val rfracrecip_to_ser_absolute : _ t -> Signed.long -> _ t
val rfracrecip : _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> Signed.long
val scalarser : _ t -> Signed.long -> Signed.long -> _ t
val sertoser : _ t -> Signed.long -> _ t
val toser_i : _ t -> _ t
val rgv_to_ser : _ t -> Signed.long -> Signed.long -> _ t
val ser0 : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val padic_to_q : _ t -> _ t
val padic_to_q_shallow : _ t -> _ t
val qpv_to_qv : _ t -> _ t
val rgc_rgv_mulrealsym : _ t -> _ t -> _ t
val rgm_mulreal : _ t -> _ t -> _ t
val rgx_cxeval : _ t -> _ t -> _ t -> _ t
val rgx_deflate_max : _ t -> Signed.long Ctypes_static.ptr -> _ t
val rgx_deflate_order : _ t -> Signed.long
val rgx_degree : _ t -> Signed.long -> Signed.long
val rgx_integ : _ t -> _ t
val rgxy_cxevalx : _ t -> _ t -> _ t -> _ t
val zx_deflate_order : _ t -> Signed.long
val zx_deflate_max : _ t -> Signed.long Ctypes_static.ptr -> _ t
val ceil_safe : _ t -> _ t
val ceilr : _ t -> _ t
val centerlift : _ t -> _ t
val centerlift0 : _ t -> Signed.long -> _ t
val compo : _ t -> Signed.long -> _ t
val deg1pol : _ t -> _ t -> Signed.long -> _ t
val deg1pol_shallow : _ t -> _ t -> Signed.long -> _ t
val deg2pol_shallow : _ t -> _ t -> _ t -> Signed.long -> _ t
val degree : _ t -> Signed.long
val denom : _ t -> _ t
val denom_i : _ t -> _ t
val denominator : _ t -> _ t -> _ t

val deriv : _ t -> Signed.long -> _ t
(** [deriv] *)

val derivn : _ t -> Signed.long -> Signed.long -> _ t
val derivser : _ t -> _ t
val diffop : _ t -> _ t -> _ t -> _ t
val diffop0 : _ t -> _ t -> _ t -> Signed.long -> _ t
val diviiround : _ t -> _ t -> _ t
val divrem : _ t -> _ t -> Signed.long -> _ t
val floor_safe : _ t -> _ t
val gceil : _ t -> _ t
val gcvtoi : _ t -> Signed.long Ctypes_static.ptr -> _ t
val gdeflate : _ t -> Signed.long -> Signed.long -> _ t
val gdivent : _ t -> _ t -> _ t
val gdiventgs : _ t -> Signed.long -> _ t
val gdiventsg : Signed.long -> _ t -> _ t
val gdiventres : _ t -> _ t -> _ t
val gdivmod : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val gdivround : _ t -> _ t -> _ t
val gdvd : _ t -> _ t -> int
val geq : _ t -> _ t -> _ t
val geval : _ t -> _ t
val gfloor : _ t -> _ t
val gtrunc2n : _ t -> Signed.long -> _ t
val gfrac : _ t -> _ t
val gge : _ t -> _ t -> _ t
val ggrando : _ t -> Signed.long -> _ t
val ggt : _ t -> _ t -> _ t
val gimag : _ t -> _ t
val gisexactzero : _ t -> _ t
val gle : _ t -> _ t -> _ t
val glt : _ t -> _ t -> _ t
val gmod : _ t -> _ t -> _ t
val gmodgs : _ t -> Signed.long -> _ t
val gmodsg : Signed.long -> _ t -> _ t
val gmodulo : _ t -> _ t -> _ t
val gmodulsg : Signed.long -> _ t -> _ t
val gmodulss : Signed.long -> Signed.long -> _ t
val gne : _ t -> _ t -> _ t
val gnot : _ t -> _ t
val gpolvar : _ t -> _ t
val gppadicprec : _ t -> _ t -> _ t
val gppoldegree : _ t -> Signed.long -> _ t
val gprecision : _ t -> Signed.long
val gpserprec : _ t -> Signed.long -> _ t
val greal : _ t -> _ t
val grndtoi : _ t -> Signed.long Ctypes_static.ptr -> _ t
val ground : _ t -> _ t
val gshift : _ t -> Signed.long -> _ t
val gsubst : _ t -> Signed.long -> _ t -> _ t
val gsubstpol : _ t -> _ t -> _ t -> _ t
val gsubstvec : _ t -> _ t -> _ t -> _ t
val gtocol : _ t -> _ t
val gtocol0 : _ t -> Signed.long -> _ t
val gtocolrev : _ t -> _ t
val gtocolrev0 : _ t -> Signed.long -> _ t
val gtopoly : _ t -> Signed.long -> _ t
val gtopolyrev : _ t -> Signed.long -> _ t
val gtovec : _ t -> _ t
val gtovec0 : _ t -> Signed.long -> _ t
val gtovecrev : _ t -> _ t
val gtovecrev0 : _ t -> Signed.long -> _ t
val gtovecsmall : _ t -> _ t
val gtovecsmall0 : _ t -> Signed.long -> _ t
val gtrunc : _ t -> _ t
val gvar : _ t -> Signed.long
val gvar2 : _ t -> Signed.long
val hqfeval : _ t -> _ t -> _ t
val imag_i : _ t -> _ t
val integ : _ t -> Signed.long -> _ t
val integser : _ t -> _ t
val ser_inv : _ t -> _ t
val iscomplex : _ t -> int
val isexactzero : _ t -> int
val isrationalzeroscalar : _ t -> int
val isinexact : _ t -> int
val isinexactreal : _ t -> int
val isint : _ t -> _ t Ctypes_static.ptr -> int
val isrationalzero : _ t -> int
val issmall : _ t -> Signed.long Ctypes_static.ptr -> int
val lift : _ t -> _ t
val lift_shallow : _ t -> _ t
val lift0 : _ t -> Signed.long -> _ t
val liftall : _ t -> _ t
val liftall_shallow : _ t -> _ t
val liftint : _ t -> _ t
val liftint_shallow : _ t -> _ t
val liftpol : _ t -> _ t
val liftpol_shallow : _ t -> _ t
val mkcoln : Signed.long -> _ t
val mkintn : Signed.long -> _ t
val mkpoln : Signed.long -> _ t
val mkvecn : Signed.long -> _ t
val mkvecsmalln : Signed.long -> _ t
val modrr_safe : _ t -> _ t -> _ t
val modrr_i : _ t -> _ t -> _ t -> _ t
val mulreal : _ t -> _ t -> _ t
val numer : _ t -> _ t
val numer_i : _ t -> _ t
val numerator : _ t -> _ t -> _ t
val padicprec : _ t -> _ t -> Signed.long
val padicprec_relative : _ t -> Signed.long
val polcoef : _ t -> Signed.long -> Signed.long -> _ t
val polcoef_i : _ t -> Signed.long -> Signed.long -> _ t
val poldegree : _ t -> Signed.long -> Signed.long
val poleval : _ t -> _ t -> _ t
val pollead : _ t -> Signed.long -> _ t
val precision : _ t -> Signed.long
val qf_apply_rgm : _ t -> _ t -> _ t
val qf_apply_zm : _ t -> _ t -> _ t
val qfb_apply_zm : _ t -> _ t -> _ t
val qfbil : _ t -> _ t -> _ t -> _ t
val qfeval : _ t -> _ t -> _ t
val qfeval0 : _ t -> _ t -> _ t -> _ t
val qfevalb : _ t -> _ t -> _ t -> _ t
val qfnorm : _ t -> _ t -> _ t
val real_i : _ t -> _ t
val round0 : _ t -> _ t Ctypes_static.ptr -> _ t
val roundr : _ t -> _ t
val roundr_safe : _ t -> _ t
val scalarpol : _ t -> Signed.long -> _ t
val scalarpol_shallow : _ t -> Signed.long -> _ t
val ser_unscale : _ t -> _ t -> _ t
val serprec : _ t -> Signed.long -> Signed.long
val serreverse : _ t -> _ t
val simplify : _ t -> _ t
val simplify_shallow : _ t -> _ t
val tayl : _ t -> Signed.long -> Signed.long -> _ t
val trunc0 : _ t -> _ t Ctypes_static.ptr -> _ t
val uu32toi : pari_ulong -> pari_ulong -> _ t
val uu32toineg : pari_ulong -> pari_ulong -> _ t
val vars_sort_inplace : _ t -> _ t
val vars_to_rgxv : _ t -> _ t
val variables_vecsmall : _ t -> _ t
val variables_vec : _ t -> _ t
val genus2red : _ t -> _ t -> _ t
val genus2igusa : _ t -> Signed.long -> _ t
val gchar_conductor : _ t -> _ t -> _ t
val gchar_identify : _ t -> _ t -> _ t -> Signed.long -> _ t
val gcharalgebraic : _ t -> _ t -> _ t
val gcharduallog : _ t -> _ t -> _ t
val gchareval : _ t -> _ t -> _ t -> Signed.long -> _ t
val gchari_lfun : _ t -> _ t -> _ t -> _ t
val gcharinit : _ t -> _ t -> Signed.long -> _ t
val gcharisalgebraic : _ t -> _ t -> _ t Ctypes_static.ptr -> int

val gcharlocal :
  _ t -> _ t -> _ t -> Signed.long -> _ t Ctypes_static.ptr -> _ t

val gcharlog : _ t -> _ t -> Signed.long -> _ t
val gcharnewprec : _ t -> Signed.long -> _ t
val is_gchar_group : _ t -> int
val lfungchar : _ t -> _ t -> _ t
val vecan_gchar : _ t -> Signed.long -> Signed.long -> _ t
val eulerf_gchar : _ t -> _ t -> Signed.long -> _ t
val group_ident : _ t -> _ t -> Signed.long
val group_ident_trans : _ t -> _ t -> Signed.long

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
  hashtable Ctypes.structure Ctypes_static.ptr -> unit Ctypes_static.ptr -> _ t

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
  (_ t -> _ t -> int) Ctypes_static.static_funptr ->
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

val hash_keys : hashtable Ctypes.structure Ctypes_static.ptr -> _ t
val hash_values : hashtable Ctypes.structure Ctypes_static.ptr -> _ t

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
val hash_gen : _ t -> pari_ulong
val hash_zv : _ t -> pari_ulong
val zx_hyperellred : _ t -> _ t Ctypes_static.ptr -> _ t
val hyperellcharpoly : _ t -> _ t
val hyperellchangecurve : _ t -> _ t -> _ t
val hyperelldisc : _ t -> _ t
val hyperellisoncurve : _ t -> _ t -> int
val hyperellminimaldisc : _ t -> _ t -> _ t
val hyperellminimalmodel : _ t -> _ t Ctypes_static.ptr -> _ t -> _ t
val hyperellpadicfrobenius0 : _ t -> _ t -> Signed.long -> _ t
val hyperellpadicfrobenius : _ t -> pari_ulong -> Signed.long -> _ t
val hyperellred : _ t -> _ t Ctypes_static.ptr -> _ t
val nfhyperellpadicfrobenius : _ t -> _ t -> pari_ulong -> Signed.long -> _ t
val hypergeom : _ t -> _ t -> _ t -> Signed.long -> _ t
val airy : _ t -> Signed.long -> _ t
val rgm_hnfall : _ t -> _ t Ctypes_static.ptr -> Signed.long -> _ t
val zm_hnf : _ t -> _ t
val zm_hnf_knapsack : _ t -> _ t
val zm_hnfall : _ t -> _ t Ctypes_static.ptr -> Signed.long -> _ t
val zm_hnfall_i : _ t -> _ t Ctypes_static.ptr -> Signed.long -> _ t
val zm_hnfcenter : _ t -> _ t
val zm_hnflll : _ t -> _ t Ctypes_static.ptr -> int -> _ t
val zv_extgcd : _ t -> _ t
val zv_snfall : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t
val zv_snf_group : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t
val zv_snf_rank_u : _ t -> pari_ulong -> Signed.long
val zv_snf_trunc : _ t -> unit
val zm_hnfmod : _ t -> _ t -> _ t
val zm_hnfmodall : _ t -> _ t -> Signed.long -> _ t
val zm_hnfmodall_i : _ t -> _ t -> Signed.long -> _ t
val zm_hnfmodid : _ t -> _ t -> _ t
val zm_hnfmodprime : _ t -> _ t -> _ t
val zm_hnfperm : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t
val zm_snfclean : _ t -> _ t -> _ t -> unit
val zm_snf : _ t -> _ t
val zm_snf_group : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t
val zm_snfall : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val zm_snfall_i :
  _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> Signed.long -> _ t

val zv_snfclean : _ t -> _ t
val zpm_echelon : _ t -> Signed.long -> _ t -> _ t -> _ t
val gsmith : _ t -> _ t
val gsmithall : _ t -> _ t
val hnf : _ t -> _ t
val hnf_divscale : _ t -> _ t -> _ t -> _ t
val hnf_invscale : _ t -> _ t -> _ t
val hnf_solve : _ t -> _ t -> _ t
val hnf_invimage : _ t -> _ t -> _ t
val hnfall : _ t -> _ t
val hnfdivide : _ t -> _ t -> int
val hnflll : _ t -> _ t
val hnfmerge_get_1 : _ t -> _ t -> _ t
val hnfmod : _ t -> _ t -> _ t
val hnfmodid : _ t -> _ t -> _ t
val hnfperm : _ t -> _ t
val matfrobenius : _ t -> Signed.long -> Signed.long -> _ t
val mathnf0 : _ t -> Signed.long -> _ t
val matsnf0 : _ t -> Signed.long -> _ t
val smith : _ t -> _ t
val smithall : _ t -> _ t
val smithclean : _ t -> _ t
val snfrank : _ t -> _ t -> Signed.long
val zlm_echelon : _ t -> Signed.long -> pari_ulong -> pari_ulong -> _ t
val zv_snf_rank : _ t -> pari_ulong -> Signed.long
val z_ecm : _ t -> Signed.long -> Signed.long -> pari_ulong -> _ t
val z_factor : _ t -> _ t
val z_factor_limit : _ t -> pari_ulong -> _ t
val z_factor_until : _ t -> _ t -> _ t
val z_issmooth : _ t -> pari_ulong -> Signed.long
val z_issmooth_fact : _ t -> pari_ulong -> _ t
val z_issquarefree : _ t -> Signed.long
val z_pollardbrent : _ t -> Signed.long -> Signed.long -> _ t
val absz_factor : _ t -> _ t
val absz_factor_limit : _ t -> pari_ulong -> _ t
val absz_factor_limit_strict : _ t -> pari_ulong -> _ t Ctypes_static.ptr -> _ t
val coreu : pari_ulong -> pari_ulong
val coreu_fact : _ t -> pari_ulong
val factorint : _ t -> Signed.long -> _ t
val factoru : pari_ulong -> _ t
val tridiv_boundu : pari_ulong -> pari_ulong
val ifac_isprime : _ t -> int

val ifac_next :
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  int

val ifac_read :
  _ t -> _ t Ctypes_static.ptr -> Signed.long Ctypes_static.ptr -> int

val ifac_skip : _ t -> unit
val ifac_start : _ t -> int -> _ t

val is_357_power :
  _ t -> _ t Ctypes_static.ptr -> pari_ulong Ctypes_static.ptr -> int

val is_pth_power :
  _ t ->
  _ t Ctypes_static.ptr ->
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  pari_ulong ->
  int

val ispowerful : _ t -> Signed.long
val maxomegau : pari_ulong -> Signed.long
val maxomegaoddu : pari_ulong -> Signed.long
val moebius : _ t -> Signed.long
val moebiusu : pari_ulong -> Signed.long
val moebiusu_fact : _ t -> Signed.long
val nextprime : _ t -> _ t
val precprime : _ t -> _ t
val radicalu : pari_ulong -> pari_ulong
val tridiv_bound : _ t -> pari_ulong

val uis_357_power :
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  int

val uis_357_powermod : pari_ulong -> pari_ulong Ctypes_static.ptr -> int
val unextprime : pari_ulong -> pari_ulong
val uprecprime : pari_ulong -> pari_ulong
val vecfactorsquarefreeu : pari_ulong -> pari_ulong -> _ t
val vecfactorsquarefreeu_coprime : pari_ulong -> pari_ulong -> _ t -> _ t
val vecfactoru_i : pari_ulong -> pari_ulong -> _ t
val vecfactoru : pari_ulong -> pari_ulong -> _ t
val vecfactoroddu_i : pari_ulong -> pari_ulong -> _ t
val vecfactoroddu : pari_ulong -> pari_ulong -> _ t
val vecsquarefreeu : pari_ulong -> pari_ulong -> _ t
val chk_gerepileupto : _ t -> int
val copy_bin : _ t -> genbin Ctypes.structure Ctypes_static.ptr
val copy_bin_canon : _ t -> genbin Ctypes.structure Ctypes_static.ptr
val dbg_gerepile : pari_ulong -> unit
val dbg_gerepileupto : _ t -> unit
val errname : _ t -> _ t
val gclone : _ t -> _ t
val gcloneref : _ t -> _ t
val gclone_refc : _ t -> unit
val gcopy : _ t -> _ t
val gcopy_avma : _ t -> pari_ulong Ctypes_static.ptr -> _ t
val gcopy_lg : _ t -> Signed.long -> _ t
val gerepile : pari_ulong -> pari_ulong -> _ t -> _ t
val gerepileallsp : pari_ulong -> pari_ulong -> int -> unit

val gerepilecoeffssp :
  pari_ulong -> pari_ulong -> Signed.long Ctypes_static.ptr -> int -> unit

val gerepilemanysp :
  pari_ulong ->
  pari_ulong ->
  _ t Ctypes_static.ptr Ctypes_static.ptr ->
  int ->
  unit

val getheap : unit -> _ t
val gsizeword : _ t -> Signed.long
val gsizebyte : _ t -> Signed.long
val gunclone : _ t -> unit
val gunclone_deep : _ t -> unit
val listcopy : _ t -> _ t
val listinit : _ t -> _ t
val msgtimer : string -> unit
val name_numerr : string -> Signed.long
val new_chunk_resize : int -> unit
val newblock : int -> _ t
val numerr_name : Signed.long -> string
val obj_check : _ t -> Signed.long -> _ t

val obj_checkbuild :
  _ t -> Signed.long -> (_ t -> _ t) Ctypes_static.static_funptr -> _ t

val obj_checkbuild_padicprec :
  _ t ->
  Signed.long ->
  (_ t -> Signed.long -> _ t) Ctypes_static.static_funptr ->
  Signed.long ->
  _ t

val obj_checkbuild_realprec :
  _ t ->
  Signed.long ->
  (_ t -> Signed.long -> _ t) Ctypes_static.static_funptr ->
  Signed.long ->
  _ t

val obj_checkbuild_prec :
  _ t ->
  Signed.long ->
  (_ t -> Signed.long -> _ t) Ctypes_static.static_funptr ->
  (_ t -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  _ t

val obj_free : _ t -> unit
val obj_init : Signed.long -> Signed.long -> _ t
val obj_insert : _ t -> Signed.long -> _ t -> _ t
val obj_insert_shallow : _ t -> Signed.long -> _ t -> _ t
val obj_reinit : _ t -> _ t
val pari_add_function : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_defaults_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_close : unit -> unit
val pari_close_opts : pari_ulong -> unit
val pari_compile_str : string -> _ t
val pari_daemon : unit -> int
val pari_err : int -> unit
val pari_err_last : unit -> _ t
val pari_err2str : _ t -> string
val pari_init_opts : int -> pari_ulong -> pari_ulong -> unit
val pari_init : int -> pari_ulong -> unit
val pari_stackcheck_init : unit Ctypes_static.ptr -> unit
val pari_sighandler : int -> unit
val pari_sig_init : (int -> unit) Ctypes_static.static_funptr -> unit

val pari_thread_alloc :
  pari_thread Ctypes.structure Ctypes_static.ptr -> int -> _ t -> unit

val pari_thread_close : unit -> unit
val pari_thread_free : pari_thread Ctypes.structure Ctypes_static.ptr -> unit
val pari_thread_init : unit -> unit
val pari_thread_start : pari_thread Ctypes.structure Ctypes_static.ptr -> _ t

val pari_thread_valloc :
  pari_thread Ctypes.structure Ctypes_static.ptr -> int -> int -> _ t -> unit

val pari_version : unit -> _ t
val pari_warn : int -> unit
val paristack_newrsize : pari_ulong -> unit
val paristack_resize : pari_ulong -> unit
val paristack_setsize : int -> int -> unit
val parivstack_resize : pari_ulong -> unit
val parivstack_reset : unit -> unit
val setalldebug : Signed.long -> unit
val setdebug : string -> Signed.long -> _ t
val shiftaddress : _ t -> Signed.long -> unit
val shiftaddress_canon : _ t -> Signed.long -> unit
val timer : unit -> Signed.long
val timer_delay : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long
val timer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val timer_printf :
  pari_timer Ctypes.structure Ctypes_static.ptr -> string -> unit

val timer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit
val timer2 : unit -> Signed.long
val trap0 : string -> _ t -> _ t -> _ t

val traverseheap :
  (_ t -> unit Ctypes_static.ptr -> unit) Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  unit

val walltimer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit

val walltimer_delay :
  pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val walltimer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long
val contfraceval : _ t -> _ t -> Signed.long -> _ t
val contfracinit : _ t -> Signed.long -> _ t

val intcirc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val intfuncinit :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  Signed.long ->
  _ t

val intnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val intnumgauss :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val intnumgaussinit : Signed.long -> Signed.long -> _ t
val intnuminit : _ t -> _ t -> Signed.long -> Signed.long -> _ t

val intnumosc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t ->
  Signed.long ->
  _ t

val intnumromb :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  Signed.long ->
  _ t

val intnumromb_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  Signed.long ->
  _ t

val prodeulerrat : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val prodnumrat : _ t -> Signed.long -> Signed.long -> _ t
val quodif : _ t -> Signed.long -> _ t
val sumeulerrat : _ t -> _ t -> Signed.long -> Signed.long -> _ t

val sumnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val sumnumap :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val sumnumapinit : _ t -> Signed.long -> _ t
val sumnuminit : _ t -> Signed.long -> _ t
val sumnumlagrangeinit : _ t -> _ t -> Signed.long -> _ t

val sumnumlagrange :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val sumnummonien :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val sumnummonieninit : _ t -> _ t -> _ t -> Signed.long -> _ t
val sumnumrat : _ t -> _ t -> Signed.long -> _ t

val sumnumsidi :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t ->
  float ->
  Signed.long ->
  _ t

val z_isanypower : _ t -> _ t Ctypes_static.ptr -> Signed.long
val z_ispow2 : _ t -> Signed.long
val z_ispowerall : _ t -> pari_ulong -> _ t Ctypes_static.ptr -> Signed.long
val z_issquareall : _ t -> _ t Ctypes_static.ptr -> Signed.long
val zn_ispower : _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val zn_issquare : _ t -> _ t -> Signed.long
val zp_issquare : _ t -> _ t -> Signed.long
val gisanypower : _ t -> _ t Ctypes_static.ptr -> Signed.long
val gissquare : _ t -> _ t
val gissquareall : _ t -> _ t Ctypes_static.ptr -> _ t
val ispolygonal : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val ispower : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long
val isprimepower : _ t -> _ t Ctypes_static.ptr -> Signed.long
val ispseudoprimepower : _ t -> _ t Ctypes_static.ptr -> Signed.long
val issquare : _ t -> Signed.long
val issquareall : _ t -> _ t Ctypes_static.ptr -> Signed.long
val sqrtint : _ t -> _ t
val sqrtint0 : _ t -> _ t Ctypes_static.ptr -> _ t
val uisprimepower : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long
val uissquare : pari_ulong -> Signed.long
val uissquareall : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val ulogintall :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val padicfields0 : _ t -> _ t -> Signed.long -> _ t
val padicfields : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val bnrclassfield : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val rnfkummer : _ t -> _ t -> Signed.long -> _ t
val is_linit : _ t -> Signed.long
val ldata_get_an : _ t -> _ t
val ldata_get_dual : _ t -> _ t
val ldata_get_gammavec : _ t -> _ t
val ldata_get_degree : _ t -> Signed.long
val ldata_get_k : _ t -> _ t
val ldata_get_k1 : _ t -> _ t
val ldata_get_conductor : _ t -> _ t
val ldata_get_rootno : _ t -> _ t
val ldata_get_residue : _ t -> _ t
val ldata_get_type : _ t -> Signed.long
val ldata_isreal : _ t -> Signed.long
val linit_get_type : _ t -> Signed.long
val linit_get_ldata : _ t -> _ t
val linit_get_tech : _ t -> _ t
val lfun_get_domain : _ t -> _ t
val lfun_get_dom : _ t -> _ t
val lfun_get_factgammavec : _ t -> _ t
val lfun_get_step : _ t -> _ t
val lfun_get_pol : _ t -> _ t
val lfun_get_residue : _ t -> _ t
val lfun_get_k2 : _ t -> _ t
val lfun_get_w2 : _ t -> _ t
val lfun_get_expot : _ t -> _ t
val lfun_get_bitprec : _ t -> Signed.long
val lfun : _ t -> _ t -> Signed.long -> _ t
val lfun0 : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfuncheckfeq : _ t -> _ t -> Signed.long -> Signed.long
val lfunconductor : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfuncost : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfuncost0 : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfuncreate : _ t -> _ t
val lfundual : _ t -> Signed.long -> _ t
val lfuneuler : _ t -> _ t -> Signed.long -> _ t
val lfunparams : _ t -> Signed.long -> _ t
val lfunan : _ t -> Signed.long -> Signed.long -> _ t
val lfunhardy : _ t -> _ t -> Signed.long -> _ t
val lfuninit : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfuninit0 : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfuninit_make : Signed.long -> _ t -> _ t -> _ t -> _ t
val lfunlambda : _ t -> _ t -> Signed.long -> _ t
val lfunlambda0 : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfunmisc_to_ldata : _ t -> _ t
val lfunmisc_to_ldata_shallow : _ t -> _ t
val lfunmisc_to_ldata_shallow_i : _ t -> _ t
val lfunorderzero : _ t -> Signed.long -> Signed.long -> Signed.long
val lfunprod_get_fact : _ t -> _ t
val lfunrootno : _ t -> Signed.long -> _ t
val lfunrootres : _ t -> Signed.long -> _ t
val lfunrtopoles : _ t -> _ t
val lfunshift : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfuntwist : _ t -> _ t -> Signed.long -> _ t
val lfuntheta : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfunthetacost0 : _ t -> _ t -> Signed.long -> Signed.long -> Signed.long
val lfunthetacost : _ t -> _ t -> Signed.long -> Signed.long -> Signed.long
val lfunthetainit : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfunthetacheckinit : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfunzeros : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val sdomain_isincl : float -> _ t -> _ t -> int
val theta_get_an : _ t -> _ t
val theta_get_k : _ t -> _ t
val theta_get_r : _ t -> _ t
val theta_get_bitprec : _ t -> Signed.long
val theta_get_m : _ t -> Signed.long
val theta_get_tdom : _ t -> _ t
val theta_get_isqrtn : _ t -> _ t
val vgaeasytheta : _ t -> int
val znchargauss : _ t -> _ t -> _ t -> Signed.long -> _ t
val dirzetak : _ t -> _ t -> _ t
val ellmoddegree : _ t -> _ t
val eta_zxn : Signed.long -> Signed.long -> _ t
val eta_product_zxn : _ t -> Signed.long -> _ t

val etaquotype :
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val galois_get_conj : _ t -> _ t
val ldata_vecan : _ t -> Signed.long -> Signed.long -> _ t
val ldata_newprec : _ t -> Signed.long -> _ t
val lfunabelianrelinit : _ t -> _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfunartin : _ t -> _ t -> _ t -> Signed.long -> Signed.long -> _ t
val lfundiv : _ t -> _ t -> Signed.long -> _ t
val lfunellmfpeters : _ t -> Signed.long -> _ t
val lfunetaquo : _ t -> _ t
val lfungenus2 : _ t -> _ t
val lfunmfspec : _ t -> Signed.long -> _ t
val lfunmul : _ t -> _ t -> Signed.long -> _ t
val lfunqf : _ t -> Signed.long -> _ t
val lfunsympow : _ t -> pari_ulong -> _ t
val lfunzetakinit : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val qfiseven : _ t -> Signed.long
val lfunquadneg : Signed.long -> Signed.long -> _ t
val zm_lll_norms : _ t -> float -> Signed.long -> _ t Ctypes_static.ptr -> _ t
val kerint : _ t -> _ t
val lll : _ t -> _ t
val lllfp : _ t -> float -> Signed.long -> _ t
val lllgen : _ t -> _ t
val lllgram : _ t -> _ t
val lllgramgen : _ t -> _ t
val lllgramint : _ t -> _ t
val lllgramkerim : _ t -> _ t
val lllgramkerimgen : _ t -> _ t
val lllint : _ t -> _ t
val lllintpartial : _ t -> _ t
val lllintpartial_inplace : _ t -> _ t
val lllkerim : _ t -> _ t
val lllkerimgen : _ t -> _ t
val matkerint0 : _ t -> Signed.long -> _ t
val qflll0 : _ t -> Signed.long -> _ t
val qflllgram0 : _ t -> Signed.long -> _ t
val gtomap : _ t -> _ t
val mapdelete : _ t -> _ t -> unit
val mapdomain : _ t -> _ t
val mapdomain_shallow : _ t -> _ t
val mapget : _ t -> _ t -> _ t
val mapisdefined : _ t -> _ t -> _ t Ctypes_static.ptr -> int
val mapput : _ t -> _ t -> _ t -> unit
val maptomat : _ t -> _ t
val maptomat_shallow : _ t -> _ t
val matpermanent : _ t -> _ t
val zm_permanent : _ t -> _ t
val dbllemma526 : float -> float -> float -> float -> float
val dblcoro526 : float -> float -> float -> float
val gammamellininv : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val gammamellininvasymp : _ t -> Signed.long -> Signed.long -> _ t
val gammamellininvinit : _ t -> Signed.long -> Signed.long -> _ t
val gammamellininvrt : _ t -> _ t -> Signed.long -> _ t
val member_a1 : _ t -> _ t
val member_a2 : _ t -> _ t
val member_a3 : _ t -> _ t
val member_a4 : _ t -> _ t
val member_a6 : _ t -> _ t
val member_area : _ t -> _ t
val member_b2 : _ t -> _ t
val member_b4 : _ t -> _ t
val member_b6 : _ t -> _ t
val member_b8 : _ t -> _ t
val member_bid : _ t -> _ t
val member_bnf : _ t -> _ t
val member_c4 : _ t -> _ t
val member_c6 : _ t -> _ t
val member_clgp : _ t -> _ t
val member_codiff : _ t -> _ t
val member_cyc : _ t -> _ t
val member_diff : _ t -> _ t
val member_disc : _ t -> _ t
val member_e : _ t -> _ t
val member_eta : _ t -> _ t
val member_f : _ t -> _ t
val member_fu : _ t -> _ t
val member_gen : _ t -> _ t
val member_group : _ t -> _ t
val member_index : _ t -> _ t
val member_j : _ t -> _ t
val member_mod : _ t -> _ t
val member_nf : _ t -> _ t
val member_no : _ t -> _ t
val member_omega : _ t -> _ t
val member_orders : _ t -> _ t
val member_p : _ t -> _ t
val member_pol : _ t -> _ t
val member_polabs : _ t -> _ t
val member_reg : _ t -> _ t
val member_r1 : _ t -> _ t
val member_r2 : _ t -> _ t
val member_roots : _ t -> _ t
val member_sign : _ t -> _ t
val member_t2 : _ t -> _ t
val member_tate : _ t -> _ t
val member_tu : _ t -> _ t
val member_zk : _ t -> _ t
val member_zkst : _ t -> _ t
val mf_get_m : _ t -> _ t
val mf_get_mindex : _ t -> _ t
val mf_get_minv : _ t -> _ t
val mf_get_basis : _ t -> _ t
val mf_get_dim : _ t -> Signed.long
val mf_get_e : _ t -> _ t
val mf_get_fields : _ t -> _ t
val mf_get_newforms : _ t -> _ t
val mf_get_space : _ t -> Signed.long
val mf_get_s : _ t -> _ t
val mfcusp_get_vmjd : _ t -> _ t
val mfnew_get_vj : _ t -> _ t
val qab_tracerel : _ t -> Signed.long -> _ t -> _ t
val qabm_tracerel : _ t -> Signed.long -> _ t -> _ t
val qabv_tracerel : _ t -> Signed.long -> _ t -> _ t
val qab_trace_init : Signed.long -> Signed.long -> _ t -> _ t -> _ t
val checkmf : _ t -> _ t
val checkmf_i : _ t -> int
val getcache : unit -> _ t
val hclassno6u : pari_ulong -> pari_ulong
val hclassno6u_no_cache : pari_ulong -> pari_ulong
val lfunmf : _ t -> _ t -> Signed.long -> _ t
val mfdelta : unit -> _ t
val mfeh : _ t -> _ t
val mfek : Signed.long -> _ t
val mftheta : _ t -> _ t
val mf_get_chi : _ t -> _ t
val mf_get_n : _ t -> Signed.long
val mf_get_nk : _ t -> _ t
val mf_get_field : _ t -> _ t
val mf_get_gn : _ t -> _ t
val mf_get_gk : _ t -> _ t
val mf_get_k : _ t -> Signed.long
val mf_get_r : _ t -> Signed.long
val mf_get_type : _ t -> Signed.long
val mfatkin : _ t -> _ t -> _ t
val mfatkineigenvalues : _ t -> Signed.long -> Signed.long -> _ t
val mfatkininit : _ t -> Signed.long -> Signed.long -> _ t
val mfbasis : _ t -> Signed.long -> _ t
val mfbd : _ t -> Signed.long -> _ t
val mfbracket : _ t -> _ t -> Signed.long -> _ t
val mfcharorder : _ t -> Signed.long
val mfcharmodulus : _ t -> Signed.long
val mfcharpol : _ t -> _ t
val mfcoef : _ t -> Signed.long -> _ t
val mfcoefs : _ t -> Signed.long -> Signed.long -> _ t
val mfconductor : _ t -> _ t -> Signed.long
val mfcosets : _ t -> _ t
val mfcuspdim : Signed.long -> Signed.long -> _ t -> Signed.long
val mfcuspisregular : _ t -> _ t -> Signed.long
val mfcusps : _ t -> _ t
val mfcuspval : _ t -> _ t -> _ t -> Signed.long -> _ t
val mfcuspwidth : _ t -> _ t -> Signed.long
val mfderiv : _ t -> Signed.long -> _ t
val mfderive2 : _ t -> Signed.long -> _ t
val mfdescribe : _ t -> _ t Ctypes_static.ptr -> _ t
val mfdim : _ t -> Signed.long -> _ t
val mfdiv : _ t -> _ t -> _ t
val mfdiv_val : _ t -> _ t -> Signed.long -> _ t
val mfeigenbasis : _ t -> _ t
val mfeigensearch : _ t -> _ t -> _ t
val mfeisenstein : Signed.long -> _ t -> _ t -> _ t
val mfeisensteindim : Signed.long -> Signed.long -> _ t -> Signed.long
val mfembed : _ t -> _ t -> _ t
val mfembed0 : _ t -> _ t -> Signed.long -> _ t
val mfeval : _ t -> _ t -> _ t -> Signed.long -> _ t
val mffields : _ t -> _ t
val mffromell : _ t -> _ t
val mffrometaquo : _ t -> Signed.long -> _ t
val mffromlfun : _ t -> Signed.long -> _ t
val mffromqf : _ t -> _ t -> _ t
val mffulldim : Signed.long -> Signed.long -> _ t -> Signed.long
val mfgaloisprojrep : _ t -> _ t -> Signed.long -> _ t
val mfgaloistype : _ t -> _ t -> _ t
val mfhecke : _ t -> _ t -> Signed.long -> _ t
val mfheckemat : _ t -> _ t -> _ t
val mfinit : _ t -> Signed.long -> _ t
val mfiscm : _ t -> _ t
val mfiscuspidal : _ t -> _ t -> Signed.long
val mfisequal : _ t -> _ t -> Signed.long -> Signed.long
val mfisetaquo : _ t -> Signed.long -> _ t
val mfkohnenbasis : _ t -> _ t
val mfkohnenbijection : _ t -> _ t
val mfkohneneigenbasis : _ t -> _ t -> _ t
val mflinear : _ t -> _ t -> _ t
val mfmanin : _ t -> Signed.long -> _ t
val mfmatembed : _ t -> _ t -> _ t
val mfmul : _ t -> _ t -> _ t
val mfnewdim : Signed.long -> Signed.long -> _ t -> Signed.long
val mfolddim : Signed.long -> Signed.long -> _ t -> Signed.long
val mfparams : _ t -> _ t
val mfperiodpol : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val mfperiodpolbasis : Signed.long -> Signed.long -> _ t
val mfpetersson : _ t -> _ t -> _ t
val mfpow : _ t -> Signed.long -> _ t
val mfsearch : _ t -> _ t -> Signed.long -> _ t
val mfshift : _ t -> Signed.long -> _ t
val mfshimura : _ t -> _ t -> Signed.long -> _ t

val mfslashexpansion :
  _ t ->
  _ t ->
  _ t ->
  Signed.long ->
  Signed.long ->
  _ t Ctypes_static.ptr ->
  Signed.long ->
  _ t

val mfspace : _ t -> _ t -> Signed.long
val mfsplit : _ t -> Signed.long -> Signed.long -> _ t
val mfsturm : _ t -> Signed.long
val mfsturmngk : Signed.long -> _ t -> Signed.long
val mfsturmnk : Signed.long -> Signed.long -> Signed.long
val mfsturm_mf : _ t -> Signed.long
val mfsymboleval : _ t -> _ t -> _ t -> Signed.long -> _ t
val mfsymbol : _ t -> _ t -> Signed.long -> _ t
val mftaylor : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val mftobasis : _ t -> _ t -> Signed.long -> _ t
val mftobasises : _ t -> _ t -> _ t
val mftocol : _ t -> Signed.long -> Signed.long -> _ t
val mftocoset : pari_ulong -> _ t -> _ t -> _ t
val mftonew : _ t -> _ t -> _ t
val mftraceform : _ t -> Signed.long -> _ t
val mftwist : _ t -> _ t -> _ t
val mfvecembed : _ t -> _ t -> _ t
val mfvectomat : _ t -> Signed.long -> Signed.long -> _ t
val fl_inv : pari_ulong -> pari_ulong -> pari_ulong
val fl_invsafe : pari_ulong -> pari_ulong -> pari_ulong

val fp_ratlift :
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  int

val zm2_mul : _ t -> _ t -> _ t
val abscmpii : _ t -> _ t -> int
val abscmprr : _ t -> _ t -> int
val absequalii : _ t -> _ t -> int
val addii_sign : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val addir_sign : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val addmulii : _ t -> _ t -> _ t -> _ t
val addmulii_inplace : _ t -> _ t -> _ t -> _ t
val addrr_sign : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val addsi_sign : Signed.long -> _ t -> Signed.long -> _ t
val addsr : Signed.long -> _ t -> _ t
val addui_sign : pari_ulong -> _ t -> Signed.long -> _ t
val addumului : pari_ulong -> pari_ulong -> _ t -> _ t
val affir : _ t -> _ t -> unit
val affrr : _ t -> _ t -> unit
val bezout : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val cbezout :
  Signed.long ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val cgcd : Signed.long -> Signed.long -> Signed.long
val clcm : Signed.long -> Signed.long -> Signed.long
val cmpii : _ t -> _ t -> int
val cmprr : _ t -> _ t -> int
val dblexpo : float -> Signed.long
val dblmantissa : float -> pari_ulong
val dbltor : float -> _ t
val diviiexact : _ t -> _ t -> _ t
val divir : _ t -> _ t -> _ t
val divis : _ t -> Signed.long -> _ t
val divis_rem : _ t -> Signed.long -> Signed.long Ctypes_static.ptr -> _ t
val absdiviu_rem : _ t -> pari_ulong -> pari_ulong Ctypes_static.ptr -> _ t
val diviuuexact : _ t -> pari_ulong -> pari_ulong -> _ t
val diviuexact : _ t -> pari_ulong -> _ t
val divri : _ t -> _ t -> _ t
val divrr : _ t -> _ t -> _ t
val divrs : _ t -> Signed.long -> _ t
val divru : _ t -> pari_ulong -> _ t
val divsi : Signed.long -> _ t -> _ t
val divsr : Signed.long -> _ t -> _ t
val divur : pari_ulong -> _ t -> _ t
val dvmdii : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val equalii : _ t -> _ t -> int
val equalrr : _ t -> _ t -> int
val floorr : _ t -> _ t
val gcdii : _ t -> _ t -> _ t
val halfgcdii : _ t -> _ t -> _ t
val int2n : Signed.long -> _ t
val int2u : pari_ulong -> _ t
val int2um1 : pari_ulong -> _ t
val int_normalize : _ t -> Signed.long -> _ t
val invmod : _ t -> _ t -> _ t Ctypes_static.ptr -> int
val invmod2bil : pari_ulong -> pari_ulong
val invr : _ t -> _ t
val mantissa_real : _ t -> Signed.long Ctypes_static.ptr -> _ t
val modii : _ t -> _ t -> _ t
val modiiz : _ t -> _ t -> _ t -> unit
val mulii : _ t -> _ t -> _ t
val mulir : _ t -> _ t -> _ t
val mulrr : _ t -> _ t -> _ t
val mulsi : Signed.long -> _ t -> _ t
val mulsr : Signed.long -> _ t -> _ t
val mulss : Signed.long -> Signed.long -> _ t
val mului : pari_ulong -> _ t -> _ t
val mulur : pari_ulong -> _ t -> _ t
val muluu : pari_ulong -> pari_ulong -> _ t
val muluui : pari_ulong -> pari_ulong -> _ t -> _ t
val pari_kernel_close : unit -> unit
val pari_kernel_init : unit -> unit
val pari_kernel_version : unit -> string
val remi2n : _ t -> Signed.long -> _ t
val rtodbl : _ t -> float
val shifti : _ t -> Signed.long -> _ t
val sqri : _ t -> _ t
val sqrr : _ t -> _ t
val sqrs : Signed.long -> _ t
val sqrtr_abs : _ t -> _ t
val sqrtremi : _ t -> _ t Ctypes_static.ptr -> _ t
val sqru : pari_ulong -> _ t
val subsr : Signed.long -> _ t -> _ t
val truedvmdii : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val truedvmdis : _ t -> Signed.long -> _ t Ctypes_static.ptr -> _ t
val truedvmdsi : Signed.long -> _ t -> _ t Ctypes_static.ptr -> _ t
val trunc2nr : _ t -> Signed.long -> _ t
val mantissa2nr : _ t -> Signed.long -> _ t
val truncr : _ t -> _ t
val ugcd : pari_ulong -> pari_ulong -> pari_ulong
val ulcm : pari_ulong -> pari_ulong -> pari_ulong
val umodiu : _ t -> pari_ulong -> pari_ulong
val vals : pari_ulong -> Signed.long
val fpc_ratlift : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpm_ratlift : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fpx_ratlift : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val qxqx_gcd : _ t -> _ t -> _ t -> _ t
val zxqx_gcd : _ t -> _ t -> _ t -> _ t
val nffactor : _ t -> _ t -> _ t
val nffactormod : _ t -> _ t -> _ t -> _ t
val nfgcd : _ t -> _ t -> _ t -> _ t -> _ t
val nfgcd_all : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val nfissquarefree : _ t -> _ t -> int
val nfroots : _ t -> _ t -> _ t
val nfroots_if_split : _ t Ctypes_static.ptr -> _ t -> _ t
val nfrootsof1 : _ t -> _ t
val polfnf : _ t -> _ t -> _ t
val rnfabelianconjgen : _ t -> _ t -> _ t
val rnfisabelian : _ t -> _ t -> Signed.long

val forpart :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  _ t ->
  _ t ->
  unit

val forpart_init :
  forpart_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  _ t ->
  _ t ->
  unit

val forpart_next : forpart_t Ctypes.structure Ctypes_static.ptr -> _ t
val forpart_prev : forpart_t Ctypes.structure Ctypes_static.ptr -> _ t
val numbpart : _ t -> _ t
val partitions : Signed.long -> _ t -> _ t -> _ t

val forperm :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  _ t ->
  unit

val forperm_init : forperm_t Ctypes.structure Ctypes_static.ptr -> _ t -> unit
val forperm_next : forperm_t Ctypes.structure Ctypes_static.ptr -> _ t

val forallsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> Signed.long -> unit

val forksubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  unit

val forsubset_next : forsubset_t Ctypes.structure Ctypes_static.ptr -> _ t

val forsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> _ t -> unit

val glambertw : _ t -> Signed.long -> Signed.long -> _ t
val mplambertw : _ t -> Signed.long -> _ t
val mplambertx : _ t -> Signed.long -> _ t
val mplambertx_logx : _ t -> _ t -> Signed.long -> _ t
val mplambertxlogx_x : _ t -> _ t -> Signed.long -> _ t
val z_to_perm : Signed.long -> _ t -> _ t
val abelian_group : _ t -> _ t
val conjclasses_repr : _ t -> Signed.long -> _ t
val cyc_pow : _ t -> Signed.long -> _ t
val cyc_pow_perm : _ t -> Signed.long -> _ t
val cyclicgroup : _ t -> Signed.long -> _ t
val dicyclicgroup : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val group_abelianhnf : _ t -> _ t -> _ t
val group_abeliansnf : _ t -> _ t -> _ t
val group_domain : _ t -> Signed.long
val group_elts : _ t -> Signed.long -> _ t
val group_export : _ t -> Signed.long -> _ t
val group_export_gap : _ t -> _ t
val group_export_magma : _ t -> _ t
val group_isa4s4 : _ t -> Signed.long
val group_isabelian : _ t -> Signed.long
val group_leftcoset : _ t -> _ t -> _ t
val group_order : _ t -> Signed.long
val group_perm_normalize : _ t -> _ t -> Signed.long
val group_quotient : _ t -> _ t -> _ t
val group_rightcoset : _ t -> _ t -> _ t
val group_set : _ t -> Signed.long -> _ t
val group_subgroup_is_faithful : _ t -> _ t -> int
val group_subgroup_isnormal : _ t -> _ t -> Signed.long
val group_subgroups : _ t -> _ t
val groupelts_solvablesubgroups : _ t -> _ t
val group_to_cc : _ t -> _ t
val groupelts_abelian_group : _ t -> _ t
val groupelts_center : _ t -> _ t
val groupelts_conj_set : _ t -> _ t -> _ t
val groupelts_conjclasses : _ t -> Signed.long Ctypes_static.ptr -> _ t
val groupelts_exponent : _ t -> Signed.long
val groupelts_quotient : _ t -> _ t -> _ t
val groupelts_set : _ t -> Signed.long -> _ t
val groupelts_to_group : _ t -> _ t
val numtoperm : Signed.long -> _ t -> _ t
val perm_commute : _ t -> _ t -> int
val perm_cycles : _ t -> _ t
val perm_order : _ t -> _ t
val perm_orderu : _ t -> pari_ulong
val perm_pow : _ t -> _ t -> _ t
val perm_powu : _ t -> pari_ulong -> _ t
val perm_sign : _ t -> Signed.long
val perm_to_gap : _ t -> _ t
val perm_to_z : _ t -> _ t
val permcycles : _ t -> _ t
val permorder : _ t -> _ t
val permsign : _ t -> Signed.long
val permtonum : _ t -> _ t
val quotient_group : _ t -> _ t -> _ t
val quotient_groupelts : _ t -> _ t
val quotient_perm : _ t -> _ t -> _ t
val quotient_subgroup_lift : _ t -> _ t -> _ t -> _ t
val subgroups_tableset : _ t -> Signed.long -> _ t
val tableset_find_index : _ t -> _ t -> Signed.long
val trivialgroup : unit -> _ t
val vec_insert : _ t -> Signed.long -> _ t -> _ t
val vec_is1to1 : _ t -> int
val vec_isconst : _ t -> int
val vecperm_orbits : _ t -> Signed.long -> _ t
val vecsmall_duplicate : _ t -> Signed.long
val vecsmall_duplicate_sorted : _ t -> Signed.long
val vecsmall_indexsort : _ t -> _ t
val vecsmall_is1to1 : _ t -> int
val vecsmall_isconst : _ t -> int
val vecsmall_sort : _ t -> unit
val vecsmall_uniq : _ t -> _ t
val vecsmall_uniq_sorted : _ t -> _ t
val vecsmall_counting_indexsort : _ t -> Signed.long -> _ t
val vecsmall_counting_sort : _ t -> Signed.long -> unit
val vecsmall_counting_uniq : _ t -> Signed.long -> _ t
val vecvecsmall_indexsort : _ t -> _ t
val vecvecsmall_max : _ t -> Signed.long
val vecvecsmall_search : _ t -> _ t -> Signed.long
val vecvecsmall_sort : _ t -> _ t
val vecvecsmall_sort_inplace : _ t -> _ t Ctypes_static.ptr -> unit
val vecvecsmall_sort_shallow : _ t -> _ t
val vecvecsmall_sort_uniq : _ t -> _ t
val mt_broadcast : _ t -> unit
val mt_nbthreads : unit -> Signed.long
val mt_queue_end : pari_mt Ctypes.structure Ctypes_static.ptr -> unit

val mt_queue_get :
  pari_mt Ctypes.structure Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  _ t

val mt_queue_start : pari_mt Ctypes.structure Ctypes_static.ptr -> _ t -> unit

val mt_queue_start_lim :
  pari_mt Ctypes.structure Ctypes_static.ptr -> _ t -> Signed.long -> unit

val mt_queue_submit :
  pari_mt Ctypes.structure Ctypes_static.ptr -> Signed.long -> _ t -> unit

val mt_sigint_block : unit -> unit
val mt_sigint_unblock : unit -> unit
val pari_mt_init : unit -> unit
val pari_mt_close : unit -> unit
val subcyclopclgp : _ t -> _ t -> Signed.long -> _ t
val subcycloiwasawa : _ t -> _ t -> Signed.long -> _ t
val subcyclohminus : _ t -> _ t -> _ t

val color_to_rgb :
  _ t ->
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
  _ t -> _ t -> _ t -> Signed.long -> Signed.long -> Signed.long -> _ t

val parplothexport :
  _ t -> _ t -> _ t -> _ t -> Signed.long -> Signed.long -> Signed.long -> _ t

val plotbox : Signed.long -> _ t -> _ t -> Signed.long -> unit
val plotclip : Signed.long -> unit
val plotcolor : Signed.long -> _ t -> _ t
val plotcopy : Signed.long -> Signed.long -> _ t -> _ t -> Signed.long -> unit
val plotcursor : Signed.long -> _ t
val plotdraw : _ t -> Signed.long -> unit
val plotexport : _ t -> _ t -> Signed.long -> _ t

val ploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  _ t

val plothexport :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  _ t

val plothraw : _ t -> _ t -> Signed.long -> _ t
val plothrawexport : _ t -> _ t -> _ t -> Signed.long -> _ t
val plothsizes : Signed.long -> _ t
val plotinit : Signed.long -> _ t -> _ t -> Signed.long -> unit
val plotkill : Signed.long -> unit
val plotline : Signed.long -> _ t -> _ t -> unit
val plotlines : Signed.long -> _ t -> _ t -> Signed.long -> unit
val plotlinetype : Signed.long -> Signed.long -> unit
val plotmove : Signed.long -> _ t -> _ t -> unit
val plotpoints : Signed.long -> _ t -> _ t -> unit
val plotpointsize : Signed.long -> _ t -> unit
val plotpointtype : Signed.long -> Signed.long -> unit
val plotrbox : Signed.long -> _ t -> _ t -> Signed.long -> unit

val plotrecth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  Signed.long ->
  _ t ->
  _ t ->
  pari_ulong ->
  Signed.long ->
  Signed.long ->
  _ t

val plotrecthraw : Signed.long -> _ t -> Signed.long -> _ t
val plotrline : Signed.long -> _ t -> _ t -> unit
val plotrmove : Signed.long -> _ t -> _ t -> unit
val plotrpoint : Signed.long -> _ t -> _ t -> unit
val plotscale : Signed.long -> _ t -> _ t -> _ t -> _ t -> unit
val plotstring : Signed.long -> string -> Signed.long -> unit
val psdraw : _ t -> Signed.long -> unit

val psploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  _ t

val psplothraw : _ t -> _ t -> Signed.long -> _ t

val rect2ps :
  _ t -> _ t -> _ t -> pari_plot Ctypes.structure Ctypes_static.ptr -> string

val rect2ps_i :
  _ t ->
  _ t ->
  _ t ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  int ->
  string

val rect2svg :
  _ t -> _ t -> _ t -> pari_plot Ctypes.structure Ctypes_static.ptr -> string

val pariplot :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  Signed.long ->
  unit

val zx_zp_root : _ t -> _ t -> _ t -> Signed.long -> _ t
val zp_appr : _ t -> _ t -> _ t
val cmp_padic : _ t -> _ t -> int
val factorpadic : _ t -> _ t -> Signed.long -> _ t
val gdeuc : _ t -> _ t -> _ t
val grem : _ t -> _ t -> _ t
val padicappr : _ t -> _ t -> _ t
val poldivrem : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val polrootspadic : _ t -> _ t -> Signed.long -> _ t
val flv_factorback : _ t -> _ t -> pari_ulong -> pari_ulong
val flxqv_factorback : _ t -> _ t -> _ t -> pari_ulong -> _ t
val fpv_factorback : _ t -> _ t -> _ t -> _ t
val fqv_factorback : _ t -> _ t -> _ t -> _ t -> _ t
val q_content : _ t -> _ t
val q_content_safe : _ t -> _ t
val q_denom : _ t -> _ t
val q_denom_safe : _ t -> _ t
val q_div_to_int : _ t -> _ t -> _ t
val q_gcd : _ t -> _ t -> _ t
val q_mul_to_int : _ t -> _ t -> _ t
val q_muli_to_int : _ t -> _ t -> _ t
val q_primitive_part : _ t -> _ t Ctypes_static.ptr -> _ t
val q_primpart : _ t -> _ t
val q_remove_denom : _ t -> _ t Ctypes_static.ptr -> _ t
val q_factor : _ t -> _ t
val q_factor_limit : _ t -> pari_ulong -> _ t

val rg_type :
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rgc_type :
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rescale_to_int : _ t -> _ t

val rgm_type :
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_type2 :
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type :
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type2 :
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_rg_type :
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_chinese_coprime : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val rgx_disc : _ t -> _ t

val rgx_extgcd :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val rgx_extgcd_simple :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val rgx_gcd : _ t -> _ t -> _ t
val rgx_gcd_simple : _ t -> _ t -> _ t
val rgx_halfgcd : _ t -> _ t -> _ t

val rgx_halfgcd_all :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val rgx_rescale_to_int : _ t -> _ t
val rgx_resultant_all : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val rgx_sturmpart : _ t -> _ t -> Signed.long
val rgx_sylvestermatrix : _ t -> _ t -> _ t

val rgx_type :
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type2 :
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type3 :
  _ t ->
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type_decode :
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val rgx_type_is_composite : Signed.long -> int
val rgxq_charpoly : _ t -> _ t -> Signed.long -> _ t
val rgxq_inv : _ t -> _ t -> _ t
val rgxq_minpoly : _ t -> _ t -> Signed.long -> _ t
val rgxq_mul : _ t -> _ t -> _ t -> _ t

val rgxq_ratlift :
  _ t ->
  _ t ->
  Signed.long ->
  Signed.long ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  int

val rgxq_sqr : _ t -> _ t -> _ t
val z_content : _ t -> _ t
val zx_content : _ t -> _ t
val centermod : _ t -> _ t -> _ t
val centermod_i : _ t -> _ t -> _ t -> _ t
val centermodii : _ t -> _ t -> _ t -> _ t
val content : _ t -> _ t
val content0 : _ t -> _ t -> _ t
val deg1_from_roots : _ t -> Signed.long -> _ t
val factor : _ t -> _ t
val factor0 : _ t -> _ t -> _ t
val factorback : _ t -> _ t
val factorback2 : _ t -> _ t -> _ t

val gbezout :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val gdivexact : _ t -> _ t -> _ t

val gen_factorback :
  _ t ->
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t -> _ t) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> _ t) Ctypes_static.static_funptr ->
  _ t

val ggcd : _ t -> _ t -> _ t
val ggcd0 : _ t -> _ t -> _ t
val ghalfgcd : _ t -> _ t -> _ t
val ginvmod : _ t -> _ t -> _ t
val glcm : _ t -> _ t -> _ t
val glcm0 : _ t -> _ t -> _ t
val newtonpoly : _ t -> _ t -> _ t
val nfrootsq : _ t -> _ t
val poldisc0 : _ t -> Signed.long -> _ t
val polisirreducible : _ t -> Signed.long
val polresultant0 : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val polsym : _ t -> Signed.long -> _ t
val primitive_part : _ t -> _ t Ctypes_static.ptr -> _ t
val primpart : _ t -> _ t
val reduceddiscsmith : _ t -> _ t
val resultant2 : _ t -> _ t -> _ t
val resultant : _ t -> _ t -> _ t
val rnfcharpoly : _ t -> _ t -> _ t -> Signed.long -> _ t
val roots_from_deg1 : _ t -> _ t
val roots_to_pol : _ t -> Signed.long -> _ t
val roots_to_pol_r1 : _ t -> Signed.long -> Signed.long -> _ t
val sturmpart : _ t -> _ t -> _ t -> Signed.long

val subresext :
  _ t -> _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> _ t

val sylvestermatrix : _ t -> _ t -> _ t
val trivial_fact : unit -> _ t
val gcdext0 : _ t -> _ t -> _ t
val polresultantext0 : _ t -> _ t -> Signed.long -> _ t
val polresultantext : _ t -> _ t -> _ t
val prime_fact : _ t -> _ t
val row_q_primpart : _ t -> _ t
val vec_q_primpart : _ t -> _ t
val vecprod : _ t -> _ t
val zv_lcm : _ t -> _ t
val flx_flxy_resultant : _ t -> _ t -> pari_ulong -> _ t
val flxx_resultant : _ t -> _ t -> pari_ulong -> Signed.long -> _ t
val fpx_fpxy_resultant : _ t -> _ t -> _ t -> _ t
val fpx_translate : _ t -> _ t -> _ t -> _ t
val fpxqx_normalize : _ t -> _ t -> _ t -> _ t
val fpxv_fpc_mul : _ t -> _ t -> _ t -> _ t
val fpxy_fpxq_evaly : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val fpxc_center : _ t -> _ t -> _ t -> _ t
val fpxm_center : _ t -> _ t -> _ t -> _ t
val fq_fp_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fq_add : _ t -> _ t -> _ t -> _ t -> _ t
val fq_div : _ t -> _ t -> _ t -> _ t -> _ t
val fq_halve : _ t -> _ t -> _ t -> _ t
val fq_inv : _ t -> _ t -> _ t -> _ t
val fq_invsafe : _ t -> _ t -> _ t -> _ t
val fq_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fq_mulu : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fq_neg : _ t -> _ t -> _ t -> _ t
val fq_neg_inv : _ t -> _ t -> _ t -> _ t
val fq_pow : _ t -> _ t -> _ t -> _ t -> _ t
val fq_powu : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fq_sqr : _ t -> _ t -> _ t -> _ t
val fq_sqrt : _ t -> _ t -> _ t -> _ t
val fq_sqrtn : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fq_sub : _ t -> _ t -> _ t -> _ t -> _ t
val fqc_fq_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fqc_fqv_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fqc_add : _ t -> _ t -> _ t -> _ t -> _ t
val fqc_sub : _ t -> _ t -> _ t -> _ t -> _ t
val fqv_red : _ t -> _ t -> _ t -> _ t
val fqv_roots_to_pol : _ t -> _ t -> _ t -> Signed.long -> _ t
val fqx_fq_add : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_fq_mul_to_monic : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_fq_sub : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_eval : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_translate : _ t -> _ t -> _ t -> _ t -> _ t

val fqxq_matrix_pow :
  _ t -> Signed.long -> Signed.long -> _ t -> _ t -> _ t -> _ t

val fqxq_powers : _ t -> Signed.long -> _ t -> _ t -> _ t -> _ t
val fqxy_eval : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fqxy_evalx : _ t -> _ t -> _ t -> _ t -> _ t
val qx_disc : _ t -> _ t
val qx_gcd : _ t -> _ t -> _ t
val qx_resultant : _ t -> _ t -> _ t
val qxq_div : _ t -> _ t -> _ t -> _ t
val qxq_intnorm : _ t -> _ t -> _ t
val qxq_inv : _ t -> _ t -> _ t
val qxq_mul : _ t -> _ t -> _ t -> _ t
val qxq_norm : _ t -> _ t -> _ t
val qxq_sqr : _ t -> _ t -> _ t
val rg_is_fp : _ t -> _ t Ctypes_static.ptr -> int
val rg_is_fpxq : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> int
val rg_to_fp : _ t -> _ t -> _ t
val rg_to_fpxq : _ t -> _ t -> _ t -> _ t
val rgc_to_fpc : _ t -> _ t -> _ t
val rgc_to_fqc : _ t -> _ t -> _ t -> _ t
val rgm_is_fpm : _ t -> _ t Ctypes_static.ptr -> int
val rgm_to_flm : _ t -> pari_ulong -> _ t
val rgm_to_fpm : _ t -> _ t -> _ t
val rgm_to_fqm : _ t -> _ t -> _ t -> _ t
val rgv_is_fpv : _ t -> _ t Ctypes_static.ptr -> int
val rgv_to_flv : _ t -> pari_ulong -> _ t
val rgv_to_fpv : _ t -> _ t -> _ t
val rgx_is_fpx : _ t -> _ t Ctypes_static.ptr -> int
val rgx_to_fpx : _ t -> _ t -> _ t
val rgx_is_fpxqx : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> int
val rgx_to_fpxqx : _ t -> _ t -> _ t -> _ t
val rgx_to_fqx : _ t -> _ t -> _ t -> _ t

val z_incremental_crt :
  _ t Ctypes_static.ptr ->
  pari_ulong ->
  _ t Ctypes_static.ptr ->
  pari_ulong ->
  int

val z_init_crt : pari_ulong -> pari_ulong -> _ t

val zm_incremental_crt :
  _ t Ctypes_static.ptr -> _ t -> _ t Ctypes_static.ptr -> pari_ulong -> int

val zm_init_crt : _ t -> pari_ulong -> _ t
val zx_zxy_resultant : _ t -> _ t -> _ t
val zx_zxy_rnfequation : _ t -> _ t -> Signed.long Ctypes_static.ptr -> _ t
val zx_disc : _ t -> _ t
val zx_gcd : _ t -> _ t -> _ t
val zx_gcd_all : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t

val zx_incremental_crt :
  _ t Ctypes_static.ptr -> _ t -> _ t Ctypes_static.ptr -> pari_ulong -> int

val zx_init_crt : _ t -> pari_ulong -> Signed.long -> _ t
val zx_is_squarefree : _ t -> int
val zx_radical : _ t -> _ t
val zx_resultant : _ t -> _ t -> _ t

val zxm_incremental_crt :
  _ t Ctypes_static.ptr -> _ t -> _ t Ctypes_static.ptr -> pari_ulong -> int

val zxm_init_crt : _ t -> Signed.long -> pari_ulong -> _ t
val zxq_minpoly : _ t -> _ t -> Signed.long -> pari_ulong -> _ t
val zxq_charpoly : _ t -> _ t -> Signed.long -> _ t
val characteristic : _ t -> _ t
val ffnbirred : _ t -> Signed.long -> _ t
val ffnbirred0 : _ t -> Signed.long -> Signed.long -> _ t
val ffsumnbirred : _ t -> Signed.long -> _ t

val get_fq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  _ t ->
  _ t ->
  bb_field Ctypes.structure Ctypes_static.ptr

val init_flxq : pari_ulong -> Signed.long -> Signed.long -> _ t
val init_fq : _ t -> Signed.long -> Signed.long -> _ t
val nfx_disc : _ t -> _ t -> _ t
val nfx_resultant : _ t -> _ t -> _ t -> _ t
val pol_x_powers : Signed.long -> Signed.long -> _ t
val residual_characteristic : _ t -> _ t
val polclass : _ t -> Signed.long -> Signed.long -> _ t
val fp_modinv_to_j : _ t -> Signed.long -> _ t -> _ t

val fp_polmodular_evalx :
  Signed.long -> Signed.long -> _ t -> _ t -> Signed.long -> int -> _ t

val check_modinv : Signed.long -> unit
val disc_best_modinv : Signed.long -> Signed.long
val modinv_height_factor : Signed.long -> Signed.long
val modinv_good_disc : Signed.long -> Signed.long -> int
val modinv_good_prime : Signed.long -> Signed.long -> int
val modinv_is_weber : Signed.long -> int
val modinv_is_double_eta : Signed.long -> int

val polmodular :
  Signed.long -> Signed.long -> _ t -> Signed.long -> Signed.long -> _ t

val polmodular_zm : Signed.long -> Signed.long -> _ t

val polmodular_zxx :
  Signed.long -> Signed.long -> Signed.long -> Signed.long -> _ t

val bpsw_isprime : _ t -> Signed.long
val bpsw_psp : _ t -> Signed.long
val addprimes : _ t -> _ t
val check_ecppcert : _ t -> Signed.long
val gisprime : _ t -> Signed.long -> _ t
val gispseudoprime : _ t -> Signed.long -> _ t
val gprimepi_upper_bound : _ t -> _ t
val gprimepi_lower_bound : _ t -> _ t
val isprime : _ t -> Signed.long
val ispseudoprime : _ t -> Signed.long -> Signed.long
val millerrabin : _ t -> Signed.long -> Signed.long
val prime : Signed.long -> _ t
val primecert : _ t -> Signed.long -> _ t
val primecert0 : _ t -> Signed.long -> Signed.long -> _ t
val primecertexport : _ t -> Signed.long -> _ t
val primecertisvalid : _ t -> Signed.long
val primepi : _ t -> _ t
val primepi_upper_bound : float -> float
val primepi_lower_bound : float -> float
val primes : Signed.long -> _ t
val primes_interval : _ t -> _ t -> _ t
val primes_interval_zv : pari_ulong -> pari_ulong -> _ t
val primes_upto_zv : pari_ulong -> _ t
val primes0 : _ t -> _ t
val primes_zv : Signed.long -> _ t
val randomprime : _ t -> _ t
val randomprime0 : _ t -> _ t -> _ t
val removeprimes : _ t -> _ t
val uis2psp : pari_ulong -> int
val uispsp : pari_ulong -> pari_ulong -> int
val uislucaspsp : pari_ulong -> int
val uisprime : pari_ulong -> int
val uisprime_101 : pari_ulong -> int
val uisprime_661 : pari_ulong -> int
val uprime : Signed.long -> pari_ulong
val uprimepi : pari_ulong -> pari_ulong
val qfauto : _ t -> _ t -> _ t
val qfauto0 : _ t -> _ t -> _ t
val qfautoexport : _ t -> Signed.long -> _ t
val qfisom : _ t -> _ t -> _ t -> _ t -> _ t
val qfisom0 : _ t -> _ t -> _ t -> _ t -> _ t
val qfisominit : _ t -> _ t -> _ t -> _ t
val qfisominit0 : _ t -> _ t -> _ t -> _ t
val qforbits : _ t -> _ t -> _ t
val qfminimize : _ t -> _ t
val qfparam : _ t -> _ t -> Signed.long -> _ t
val qfsolve : _ t -> _ t
val z_isfundamental : _ t -> Signed.long
val classno : _ t -> _ t
val classno2 : _ t -> _ t
val hclassnof_fact : _ t -> _ t -> _ t -> _ t
val hclassno : _ t -> _ t
val hclassno6 : _ t -> _ t
val isfundamental : _ t -> Signed.long
val qfb_equal1 : _ t -> int
val qfbclassno0 : _ t -> Signed.long -> _ t
val qfi_shanks : _ t -> _ t -> Signed.long -> _ t
val qfi_log : _ t -> _ t -> _ t -> _ t
val qfi_order : _ t -> _ t -> _ t
val quadclassnof : _ t -> _ t Ctypes_static.ptr -> _ t
val quadclassnof_fact : _ t -> _ t -> _ t -> _ t
val quaddisc : _ t -> _ t
val quadregulator : _ t -> Signed.long -> _ t
val quadunit : _ t -> _ t
val quadunit0 : _ t -> Signed.long -> _ t
val quadunitindex : _ t -> _ t -> _ t
val quadunitnorm : _ t -> Signed.long
val sisfundamental : Signed.long -> Signed.long
val uhclassnof_fact : _ t -> Signed.long -> Signed.long
val unegisfundamental : pari_ulong -> Signed.long
val unegquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val uposisfundamental : pari_ulong -> Signed.long
val uposquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val uquadclassnof_fact : pari_ulong -> Signed.long -> _ t -> _ t -> pari_ulong
val zn_quad_roots : _ t -> _ t -> _ t -> _ t
val genrand : _ t -> _ t
val getrand : unit -> _ t
val pari_rand : unit -> pari_ulong
val randomi : _ t -> _ t
val randomr : Signed.long -> _ t
val random_f2x : Signed.long -> Signed.long -> _ t
val random_fl : pari_ulong -> pari_ulong
val random_bits : Signed.long -> Signed.long
val random_zv : Signed.long -> _ t
val setrand : _ t -> unit
val ellratpoints : _ t -> _ t -> Signed.long -> _ t
val hyperellratpoints : _ t -> _ t -> Signed.long -> _ t
val qx_complex_roots : _ t -> Signed.long -> _ t
val fft : _ t -> _ t -> _ t
val fftinv : _ t -> _ t -> _ t
val cleanroots : _ t -> Signed.long -> _ t
val fujiwara_bound : _ t -> float
val fujiwara_bound_real : _ t -> Signed.long -> float
val isrealappr : _ t -> Signed.long -> int
val polgraeffe : _ t -> _ t
val polmod_to_embed : _ t -> Signed.long -> _ t
val polrootsbound : _ t -> _ t -> _ t
val roots : _ t -> Signed.long -> _ t
val realroots : _ t -> _ t -> Signed.long -> _ t
val zx_graeffe : _ t -> _ t
val zx_realroots_irred : _ t -> Signed.long -> _ t
val zx_sturm : _ t -> Signed.long
val zx_sturm_irred : _ t -> Signed.long
val zx_sturmpart : _ t -> _ t -> Signed.long
val zx_uspensky : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val factor_aurifeuille : _ t -> Signed.long -> _ t
val factor_aurifeuille_prime : _ t -> Signed.long -> _ t
val galoissubcyclo : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val polsubcyclo : Signed.long -> Signed.long -> Signed.long -> _ t
val polsubcyclofast : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val znsubgroupgenerators : _ t -> Signed.long -> _ t
val nfsubfields : _ t -> Signed.long -> _ t
val nfsubfields0 : _ t -> Signed.long -> Signed.long -> _ t
val nfsubfieldscm : _ t -> Signed.long -> _ t
val nfsubfieldsmax : _ t -> Signed.long -> _ t
val nflist : _ t -> _ t -> Signed.long -> _ t -> _ t
val nfresolvent : _ t -> Signed.long -> _ t
val subgrouplist : _ t -> _ t -> _ t

val forsubgroup :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  unit

val abmap_kernel : _ t -> _ t
val abmap_subgroup_image : _ t -> _ t -> _ t
val bnrl1 : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val bnrrootnumber : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val bnrstark : _ t -> _ t -> Signed.long -> _ t
val cyc2elts : _ t -> _ t
val qfbforms : _ t -> _ t
val quadhilbert : _ t -> Signed.long -> _ t
val quadray : _ t -> _ t -> Signed.long -> _ t
val chartogenstr : char -> _ t
val pari_strdup : string -> string
val pari_strndup : string -> Signed.long -> string
val stack_strcat : string -> string -> string
val stack_strdup : string -> string
val pari_strchr : _ t -> _ t
val strjoin : _ t -> _ t -> _ t
val strntogenstr : string -> Signed.long -> _ t
val strsplit : _ t -> _ t -> _ t
val strtogenstr : string -> _ t
val type_name : Signed.long -> string
val type0 : _ t -> _ t

val asympnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val asympnumraw :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  Signed.long ->
  _ t ->
  Signed.long ->
  _ t

val derivnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val derivnumk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val derivfun :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val derivfunk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val forvec_init :
  forvec_t Ctypes.structure Ctypes_static.ptr -> _ t -> Signed.long -> int

val forvec_next : forvec_t Ctypes.structure Ctypes_static.ptr -> _ t

val laurentseries :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  _ t

val limitnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val polzag : Signed.long -> Signed.long -> _ t

val prodeuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val prodinf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val prodinf1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val solvestep :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  _ t ->
  Signed.long ->
  Signed.long ->
  _ t

val sumalt :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val sumalt2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val sumpos :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val sumpos2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val suminf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val suminf_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val sumdivmultexpr :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t

val zbrent :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  _ t ->
  Signed.long ->
  _ t

val bnfisintnorm : _ t -> _ t -> _ t
val bnfisintnormabs : _ t -> _ t -> _ t
val ideals_by_norm : _ t -> _ t -> _ t
val thue : _ t -> _ t -> _ t -> _ t
val thueinit : _ t -> Signed.long -> Signed.long -> _ t
val pi2n : Signed.long -> Signed.long -> _ t
val pii2 : Signed.long -> _ t
val pii2n : Signed.long -> Signed.long -> _ t
val qp_exp : _ t -> _ t
val qp_exp_prec : _ t -> Signed.long
val qp_log : _ t -> _ t
val qp_sqrt : _ t -> _ t
val qp_sqrtn : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val zn_sqrt : _ t -> _ t -> _ t
val zp_teichmuller : _ t -> _ t -> Signed.long -> _ t -> _ t
val agm : _ t -> _ t -> Signed.long -> _ t
val constcatalan : Signed.long -> _ t
val consteuler : Signed.long -> _ t
val constlog2 : Signed.long -> _ t
val constpi : Signed.long -> _ t
val cxexpm1 : _ t -> Signed.long -> _ t
val elle : _ t -> Signed.long -> _ t
val ellk : _ t -> Signed.long -> _ t
val expir : _ t -> _ t
val exp1r_abs : _ t -> _ t
val gcos : _ t -> Signed.long -> _ t
val gcotan : _ t -> Signed.long -> _ t
val gcotanh : _ t -> Signed.long -> _ t
val gexp : _ t -> Signed.long -> _ t
val gexpm1 : _ t -> Signed.long -> _ t
val glog : _ t -> Signed.long -> _ t
val glog1p : _ t -> Signed.long -> _ t
val gpow : _ t -> _ t -> Signed.long -> _ t
val gpowers : _ t -> Signed.long -> _ t
val gpowers0 : _ t -> Signed.long -> _ t -> _ t
val gpowgs : _ t -> Signed.long -> _ t
val grootsof1 : Signed.long -> Signed.long -> _ t
val gsin : _ t -> Signed.long -> _ t
val gsinc : _ t -> Signed.long -> _ t

val gsincos :
  _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> Signed.long -> unit

val gsqrpowers : _ t -> Signed.long -> _ t
val gsqrt : _ t -> Signed.long -> _ t
val gsqrtn : _ t -> _ t -> _ t Ctypes_static.ptr -> Signed.long -> _ t
val gtan : _ t -> Signed.long -> _ t
val logr_abs : _ t -> _ t
val mpcos : _ t -> _ t
val mpeuler : Signed.long -> _ t
val mpcatalan : Signed.long -> _ t
val mpsincosm1 : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit
val mpexp : _ t -> _ t
val mpexpm1 : _ t -> _ t
val mplog : _ t -> _ t
val mplog2 : Signed.long -> _ t
val mppi : Signed.long -> _ t
val mpsin : _ t -> _ t
val mpsincos : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit
val pow2pis : _ t -> Signed.long -> _ t
val powpis : _ t -> Signed.long -> _ t
val powcx : _ t -> _ t -> _ t -> Signed.long -> _ t
val powcx_prec : Signed.long -> _ t -> Signed.long -> Signed.long
val powersr : _ t -> Signed.long -> _ t
val powiu : _ t -> pari_ulong -> _ t
val powrfrac : _ t -> Signed.long -> Signed.long -> _ t
val powrs : _ t -> Signed.long -> _ t
val powrshalf : _ t -> Signed.long -> _ t
val powru : _ t -> pari_ulong -> _ t
val powruhalf : _ t -> pari_ulong -> _ t
val powuu : pari_ulong -> pari_ulong -> _ t
val powgi : _ t -> _ t -> _ t
val rootsof1_cx : _ t -> Signed.long -> _ t
val rootsof1u_cx : pari_ulong -> Signed.long -> _ t
val rootsof1q_cx : Signed.long -> Signed.long -> Signed.long -> _ t
val rootsof1powinit : Signed.long -> Signed.long -> Signed.long -> _ t
val rootsof1pow : _ t -> Signed.long -> _ t
val serchop : _ t -> Signed.long -> _ t
val serchop_i : _ t -> Signed.long -> _ t
val serchop0 : _ t -> _ t
val sqrtnint : _ t -> Signed.long -> _ t
val sqrtnr_abs : _ t -> Signed.long -> _ t
val teich : _ t -> _ t
val teichmullerinit : Signed.long -> Signed.long -> _ t
val teichmuller : _ t -> _ t -> _ t

val trans_eval :
  string ->
  (_ t -> Signed.long -> _ t) Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val trans_evalgen :
  string ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long -> _ t)
  Ctypes_static.static_funptr ->
  _ t ->
  Signed.long ->
  _ t

val upowuu : pari_ulong -> pari_ulong -> pari_ulong
val upowers : pari_ulong -> Signed.long -> _ t
val usqrtn : pari_ulong -> pari_ulong -> pari_ulong
val usqrt : pari_ulong -> pari_ulong
val qp_gamma : _ t -> _ t
val atanhuu : pari_ulong -> pari_ulong -> Signed.long -> _ t
val atanhui : pari_ulong -> _ t -> Signed.long -> _ t
val gacosh : _ t -> Signed.long -> _ t
val gacos : _ t -> Signed.long -> _ t
val garg : _ t -> Signed.long -> _ t
val gasinh : _ t -> Signed.long -> _ t
val gasin : _ t -> Signed.long -> _ t
val gatan : _ t -> Signed.long -> _ t
val gatanh : _ t -> Signed.long -> _ t
val gcosh : _ t -> Signed.long -> _ t
val ggammah : _ t -> Signed.long -> _ t
val ggamma : _ t -> Signed.long -> _ t
val ggamma1m1 : _ t -> Signed.long -> _ t
val glngamma : _ t -> Signed.long -> _ t
val gpsi : _ t -> Signed.long -> _ t
val gsinh : _ t -> Signed.long -> _ t
val gtanh : _ t -> Signed.long -> _ t
val mpfactr : Signed.long -> Signed.long -> _ t
val mpsinhcosh : _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit
val psi1series : Signed.long -> Signed.long -> Signed.long -> _ t
val sumformal : _ t -> Signed.long -> _ t

val rgv_is_arithprog :
  _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> int

val besseljzero : _ t -> Signed.long -> Signed.long -> _ t
val besselyzero : _ t -> Signed.long -> Signed.long -> _ t
val constzeta : Signed.long -> Signed.long -> _ t
val cxek : _ t -> Signed.long -> Signed.long -> _ t
val dblmodulus : _ t -> float
val dilog : _ t -> Signed.long -> _ t
val eint1 : _ t -> Signed.long -> _ t
val expipir : _ t -> Signed.long -> _ t
val expipic : _ t -> Signed.long -> _ t
val expixy : _ t -> _ t -> Signed.long -> _ t
val eta : _ t -> Signed.long -> _ t
val eta0 : _ t -> Signed.long -> Signed.long -> _ t
val gerfc : _ t -> Signed.long -> _ t
val gpolylog : Signed.long -> _ t -> Signed.long -> _ t
val gzeta : _ t -> Signed.long -> _ t
val hbessel1 : _ t -> _ t -> Signed.long -> _ t
val hbessel2 : _ t -> _ t -> Signed.long -> _ t
val hyperu : _ t -> _ t -> _ t -> Signed.long -> _ t
val ibessel : _ t -> _ t -> Signed.long -> _ t
val incgam : _ t -> _ t -> Signed.long -> _ t
val incgam0 : _ t -> _ t -> _ t -> Signed.long -> _ t
val incgamc : _ t -> _ t -> Signed.long -> _ t
val jbessel : _ t -> _ t -> Signed.long -> _ t
val jbesselh : _ t -> _ t -> Signed.long -> _ t
val jell : _ t -> Signed.long -> _ t
val kbessel : _ t -> _ t -> Signed.long -> _ t
val mpeint1 : _ t -> _ t -> _ t
val mpveceint1 : _ t -> _ t -> Signed.long -> _ t
val polylog0 : Signed.long -> _ t -> Signed.long -> Signed.long -> _ t
val sumdedekind : _ t -> _ t -> _ t
val sumdedekind_coprime : _ t -> _ t -> _ t
val szeta : Signed.long -> Signed.long -> _ t
val theta : _ t -> _ t -> Signed.long -> _ t
val thetanullk : _ t -> Signed.long -> Signed.long -> _ t
val trueeta : _ t -> Signed.long -> _ t
val u_sumdedekind_coprime : Signed.long -> Signed.long -> _ t
val upper_to_cx : _ t -> Signed.long Ctypes_static.ptr -> _ t
val veceint1 : _ t -> _ t -> Signed.long -> _ t
val vecthetanullk : _ t -> Signed.long -> Signed.long -> _ t
val vecthetanullk_tau : _ t -> Signed.long -> Signed.long -> _ t
val veczeta : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val weber0 : _ t -> Signed.long -> Signed.long -> _ t
val weberf : _ t -> Signed.long -> _ t
val weberf1 : _ t -> Signed.long -> _ t
val weberf2 : _ t -> Signed.long -> _ t
val ybessel : _ t -> _ t -> Signed.long -> _ t
val sl2_inv_shallow : _ t -> _ t
val qevproj_apply : _ t -> _ t -> _ t
val qevproj_apply_vecei : _ t -> _ t -> Signed.long -> _ t
val qevproj_down : _ t -> _ t -> _ t
val qevproj_init : _ t -> _ t
val rgx_act_gl2q : _ t -> Signed.long -> _ t
val rgx_act_zgl2q : _ t -> Signed.long -> _ t
val checkms : _ t -> unit
val checkmspadic : _ t -> unit
val ellpadiclambdamu : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val mfnumcusps : _ t -> _ t
val mfnumcusps_fact : _ t -> _ t
val mfnumcuspsu_fact : _ t -> pari_ulong
val mfnumcuspsu : pari_ulong -> pari_ulong
val msfromcusp : _ t -> _ t -> _ t
val msfromell : _ t -> Signed.long -> _ t
val msfromhecke : _ t -> _ t -> _ t -> _ t
val msdim : _ t -> Signed.long
val mseval2_ooq : _ t -> _ t -> _ t -> _ t
val msgetlevel : _ t -> Signed.long
val msgetsign : _ t -> Signed.long
val msgetweight : _ t -> Signed.long
val msatkinlehner : _ t -> Signed.long -> _ t -> _ t
val mscuspidal : _ t -> Signed.long -> _ t
val mseisenstein : _ t -> _ t
val mseval : _ t -> _ t -> _ t -> _ t
val mshecke : _ t -> Signed.long -> _ t -> _ t
val msinit : _ t -> _ t -> Signed.long -> _ t
val msissymbol : _ t -> _ t -> _ t
val mslattice : _ t -> _ t -> _ t
val msomseval : _ t -> _ t -> _ t -> _ t

val mspadic_parse_chi :
  _ t -> _ t Ctypes_static.ptr -> _ t Ctypes_static.ptr -> unit

val mspadic_unit_eigenvalue : _ t -> Signed.long -> _ t -> Signed.long -> _ t
val mspadicinit : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val mspadicl : _ t -> _ t -> Signed.long -> _ t
val mspadicmoments : _ t -> _ t -> Signed.long -> _ t
val mspadicseries : _ t -> Signed.long -> _ t
val mspathgens : _ t -> _ t
val mspathlog : _ t -> _ t -> _ t
val msnew : _ t -> _ t
val mspetersson : _ t -> _ t -> _ t -> _ t
val mspolygon : _ t -> Signed.long -> _ t
val msstar : _ t -> _ t -> _ t
val msqexpansion : _ t -> _ t -> pari_ulong -> _ t
val mssplit : _ t -> _ t -> Signed.long -> _ t
val mstooms : _ t -> _ t -> _ t
val mscosets0 : _ t -> _ t -> _ t

val mscosets :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  _ t

val msfarey :
  _ t ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> _ t -> Signed.long) Ctypes_static.static_funptr ->
  _ t Ctypes_static.ptr ->
  _ t

val msfarey0 : _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val checkfarey_i : _ t -> int
val polylogmult : _ t -> _ t -> Signed.long -> _ t
val polylogmult_interpolate : _ t -> _ t -> _ t -> Signed.long -> _ t
val zetamult : _ t -> Signed.long -> _ t
val zetamultdual : _ t -> _ t
val zetamult_interpolate : _ t -> _ t -> Signed.long -> _ t
val zetamultall : Signed.long -> Signed.long -> Signed.long -> _ t
val zetamultconvert : _ t -> Signed.long -> _ t
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
val abscmpiu : _ t -> pari_ulong -> int
val abscmpui : pari_ulong -> _ t -> int
val absequaliu : _ t -> pari_ulong -> int
val absi : _ t -> _ t
val absi_shallow : _ t -> _ t
val absr : _ t -> _ t
val absrnz_equal1 : _ t -> int
val absrnz_equal2n : _ t -> int
val addii : _ t -> _ t -> _ t
val addiiz : _ t -> _ t -> _ t -> unit
val addir : _ t -> _ t -> _ t
val addirz : _ t -> _ t -> _ t -> unit
val addis : _ t -> Signed.long -> _ t
val addri : _ t -> _ t -> _ t
val addriz : _ t -> _ t -> _ t -> unit
val addrr : _ t -> _ t -> _ t
val addrrz : _ t -> _ t -> _ t -> unit
val addrs : _ t -> Signed.long -> _ t
val addsi : Signed.long -> _ t -> _ t
val addsiz : Signed.long -> _ t -> _ t -> unit
val addsrz : Signed.long -> _ t -> _ t -> unit
val addss : Signed.long -> Signed.long -> _ t
val addssz : Signed.long -> Signed.long -> _ t -> unit
val adduu : pari_ulong -> pari_ulong -> _ t
val affii : _ t -> _ t -> unit
val affiz : _ t -> _ t -> unit
val affrr_fixlg : _ t -> _ t -> unit
val affsi : Signed.long -> _ t -> unit
val affsr : Signed.long -> _ t -> unit
val affsz : Signed.long -> _ t -> unit
val affui : pari_ulong -> _ t -> unit
val affur : pari_ulong -> _ t -> unit
val cgetg : Signed.long -> Signed.long -> _ t
val cgetg_block : Signed.long -> Signed.long -> _ t
val cgetg_copy : _ t -> Signed.long Ctypes_static.ptr -> _ t
val cgeti : Signed.long -> _ t
val cgetineg : Signed.long -> _ t
val cgetipos : Signed.long -> _ t
val cgetr : Signed.long -> _ t
val cgetr_block : Signed.long -> _ t
val cmpir : _ t -> _ t -> int
val cmpis : _ t -> Signed.long -> int
val cmpiu : _ t -> pari_ulong -> int
val cmpri : _ t -> _ t -> int
val cmprs : _ t -> Signed.long -> int
val cmpsi : Signed.long -> _ t -> int
val cmpsr : Signed.long -> _ t -> int
val cmpss : Signed.long -> Signed.long -> int
val cmpui : pari_ulong -> _ t -> int
val cmpuu : pari_ulong -> pari_ulong -> int
val divii : _ t -> _ t -> _ t
val diviiz : _ t -> _ t -> _ t -> unit
val divirz : _ t -> _ t -> _ t -> unit
val divisz : _ t -> Signed.long -> _ t -> unit
val divriz : _ t -> _ t -> _ t -> unit
val divrrz : _ t -> _ t -> _ t -> unit
val divrsz : _ t -> Signed.long -> _ t -> unit
val divsi_rem : Signed.long -> _ t -> Signed.long Ctypes_static.ptr -> _ t
val divsiz : Signed.long -> _ t -> _ t -> unit
val divsrz : Signed.long -> _ t -> _ t -> unit
val divss : Signed.long -> Signed.long -> _ t

val divss_rem :
  Signed.long -> Signed.long -> Signed.long Ctypes_static.ptr -> _ t

val divssz : Signed.long -> Signed.long -> _ t -> unit
val dvdii : _ t -> _ t -> int
val dvdiiz : _ t -> _ t -> _ t -> int
val dvdis : _ t -> Signed.long -> int
val dvdisz : _ t -> Signed.long -> _ t -> int
val dvdiu : _ t -> pari_ulong -> int
val dvdiuz : _ t -> pari_ulong -> _ t -> int
val dvdsi : Signed.long -> _ t -> int
val dvdui : pari_ulong -> _ t -> int
val dvmdiiz : _ t -> _ t -> _ t -> _ t -> unit
val dvmdis : _ t -> Signed.long -> _ t Ctypes_static.ptr -> _ t
val dvmdisz : _ t -> Signed.long -> _ t -> _ t -> unit
val dvmdsbil : Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long
val dvmdsi : Signed.long -> _ t -> _ t Ctypes_static.ptr -> _ t
val dvmdsiz : Signed.long -> _ t -> _ t -> _ t -> unit
val dvmdss : Signed.long -> Signed.long -> _ t Ctypes_static.ptr -> _ t
val dvmdssz : Signed.long -> Signed.long -> _ t -> _ t -> unit
val dvmdubil : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val equalis : _ t -> Signed.long -> int
val equalsi : Signed.long -> _ t -> int
val equalui : pari_ulong -> _ t -> int
val equaliu : _ t -> pari_ulong -> int
val absequalui : pari_ulong -> _ t -> int
val ceildivuu : pari_ulong -> pari_ulong -> pari_ulong
val evalexpo : Signed.long -> Signed.long
val evallg : Signed.long -> Signed.long
val evalprecp : Signed.long -> Signed.long
val evalvalp : Signed.long -> Signed.long
val evalvalser : Signed.long -> Signed.long
val expi : _ t -> Signed.long
val expu : pari_ulong -> Signed.long
val fixlg : _ t -> Signed.long -> unit
val fractor : _ t -> Signed.long -> _ t
val gc_bool : pari_ulong -> int -> int
val gc_const : pari_ulong -> _ t -> _ t
val gc_double : pari_ulong -> float -> float
val gc_int : pari_ulong -> int -> int
val gc_long : pari_ulong -> Signed.long -> Signed.long
val gc_stoi : pari_ulong -> Signed.long -> _ t
val gc_ulong : pari_ulong -> pari_ulong -> pari_ulong
val gc_utoi : pari_ulong -> pari_ulong -> _ t
val gc_utoipos : pari_ulong -> pari_ulong -> _ t
val gc_null : pari_ulong -> _ t
val icopy : _ t -> _ t
val icopyspec : _ t -> Signed.long -> _ t
val int_bit : _ t -> Signed.long -> pari_ulong
val itor : _ t -> Signed.long -> _ t
val itos : _ t -> Signed.long
val itos_or_0 : _ t -> Signed.long
val itou : _ t -> pari_ulong
val itou_or_0 : _ t -> pari_ulong
val leafcopy : _ t -> _ t
val maxdd : float -> float -> float
val maxss : Signed.long -> Signed.long -> Signed.long
val maxuu : pari_ulong -> pari_ulong -> Signed.long
val mindd : float -> float -> float
val minss : Signed.long -> Signed.long -> Signed.long
val minuu : pari_ulong -> pari_ulong -> Signed.long
val mod16 : _ t -> Signed.long
val mod2 : _ t -> Signed.long
val mod2bil : _ t -> pari_ulong
val mod32 : _ t -> Signed.long
val mod4 : _ t -> Signed.long
val mod64 : _ t -> Signed.long
val mod8 : _ t -> Signed.long
val modis : _ t -> Signed.long -> _ t
val modisz : _ t -> Signed.long -> _ t -> unit
val modsi : Signed.long -> _ t -> _ t
val modsiz : Signed.long -> _ t -> _ t -> unit
val modss : Signed.long -> Signed.long -> _ t
val modssz : Signed.long -> Signed.long -> _ t -> unit
val mpabs : _ t -> _ t
val mpabs_shallow : _ t -> _ t
val mpadd : _ t -> _ t -> _ t
val mpaddz : _ t -> _ t -> _ t -> unit
val mpaff : _ t -> _ t -> unit
val mpceil : _ t -> _ t
val mpcmp : _ t -> _ t -> int
val mpcopy : _ t -> _ t
val mpdiv : _ t -> _ t -> _ t
val mpexpo : _ t -> Signed.long
val mpfloor : _ t -> _ t
val mpmul : _ t -> _ t -> _ t
val mpmulz : _ t -> _ t -> _ t -> unit
val mpneg : _ t -> _ t
val mpodd : _ t -> int
val mpround : _ t -> _ t
val mpshift : _ t -> Signed.long -> _ t
val mpsqr : _ t -> _ t
val mpsub : _ t -> _ t -> _ t
val mpsubz : _ t -> _ t -> _ t -> unit
val mptrunc : _ t -> _ t
val muliiz : _ t -> _ t -> _ t -> unit
val mulirz : _ t -> _ t -> _ t -> unit
val mulis : _ t -> Signed.long -> _ t
val muliu : _ t -> pari_ulong -> _ t
val mulri : _ t -> _ t -> _ t
val mulriz : _ t -> _ t -> _ t -> unit
val mulrrz : _ t -> _ t -> _ t -> unit
val mulrs : _ t -> Signed.long -> _ t
val mulru : _ t -> pari_ulong -> _ t
val mulsiz : Signed.long -> _ t -> _ t -> unit
val mulsrz : Signed.long -> _ t -> _ t -> unit
val mulssz : Signed.long -> Signed.long -> _ t -> unit
val negi : _ t -> _ t
val negr : _ t -> _ t
val new_chunk : int -> _ t
val rcopy : _ t -> _ t
val rdivii : _ t -> _ t -> Signed.long -> _ t
val rdiviiz : _ t -> _ t -> _ t -> unit
val rdivis : _ t -> Signed.long -> Signed.long -> _ t
val rdivsi : Signed.long -> _ t -> Signed.long -> _ t
val rdivss : Signed.long -> Signed.long -> Signed.long -> _ t
val real2n : Signed.long -> Signed.long -> _ t
val real_m2n : Signed.long -> Signed.long -> _ t
val real_0 : Signed.long -> _ t
val real_0_bit : Signed.long -> _ t
val real_1 : Signed.long -> _ t
val real_1_bit : Signed.long -> _ t
val real_m1 : Signed.long -> _ t
val remii : _ t -> _ t -> _ t
val remiiz : _ t -> _ t -> _ t -> unit
val remis : _ t -> Signed.long -> _ t
val remisz : _ t -> Signed.long -> _ t -> unit

val remlll_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val remsi : Signed.long -> _ t -> _ t
val remsiz : Signed.long -> _ t -> _ t -> unit
val remss : Signed.long -> Signed.long -> _ t
val remssz : Signed.long -> Signed.long -> _ t -> unit
val rtor : _ t -> Signed.long -> _ t
val sdivsi : Signed.long -> _ t -> Signed.long

val sdivsi_rem :
  Signed.long -> _ t -> Signed.long Ctypes_static.ptr -> Signed.long

val sdivss_rem :
  Signed.long -> Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val get_avma : unit -> pari_ulong
val set_avma : pari_ulong -> unit

val uabsdiviu_rem :
  _ t -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val udivuu_rem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val umodi2n : _ t -> Signed.long -> pari_ulong
val setabssign : _ t -> unit

val shift_left :
  _ t -> _ t -> Signed.long -> Signed.long -> pari_ulong -> pari_ulong -> unit

val shift_right :
  _ t -> _ t -> Signed.long -> Signed.long -> pari_ulong -> pari_ulong -> unit

val shiftl : pari_ulong -> pari_ulong -> pari_ulong
val shiftlr : pari_ulong -> pari_ulong -> pari_ulong
val shiftr : _ t -> Signed.long -> _ t
val shiftr_inplace : _ t -> Signed.long -> unit
val smodis : _ t -> Signed.long -> Signed.long
val smodss : Signed.long -> Signed.long -> Signed.long
val stackdummy : pari_ulong -> pari_ulong -> unit
val stack_malloc : int -> string
val stack_malloc_align : int -> Signed.long -> string
val stack_calloc : int -> string
val stack_calloc_align : int -> Signed.long -> string
val subii : _ t -> _ t -> _ t
val subiiz : _ t -> _ t -> _ t -> unit
val subir : _ t -> _ t -> _ t
val subirz : _ t -> _ t -> _ t -> unit
val subis : _ t -> Signed.long -> _ t
val subisz : _ t -> Signed.long -> _ t -> unit
val subri : _ t -> _ t -> _ t
val subriz : _ t -> _ t -> _ t -> unit
val subrr : _ t -> _ t -> _ t
val subrrz : _ t -> _ t -> _ t -> unit
val subrs : _ t -> Signed.long -> _ t
val subrsz : _ t -> Signed.long -> _ t -> unit
val subsi : Signed.long -> _ t -> _ t
val subsiz : Signed.long -> _ t -> _ t -> unit
val subsrz : Signed.long -> _ t -> _ t -> unit
val subss : Signed.long -> Signed.long -> _ t
val subssz : Signed.long -> Signed.long -> _ t -> unit
val subuu : pari_ulong -> pari_ulong -> _ t
val togglesign : _ t -> unit
val togglesign_safe : _ t Ctypes_static.ptr -> unit
val affectsign : _ t -> _ t -> unit
val affectsign_safe : _ t -> _ t Ctypes_static.ptr -> unit
val truedivii : _ t -> _ t -> _ t
val truedivis : _ t -> Signed.long -> _ t
val truedivsi : Signed.long -> _ t -> _ t

val uabsdivui_rem :
  pari_ulong -> _ t -> pari_ulong Ctypes_static.ptr -> pari_ulong

val umodsu : Signed.long -> pari_ulong -> pari_ulong
val umodui : pari_ulong -> _ t -> pari_ulong
val ugcdiu : _ t -> pari_ulong -> pari_ulong
val ugcdui : pari_ulong -> _ t -> pari_ulong
val umuluu_le : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val umuluu_or_0 : pari_ulong -> pari_ulong -> pari_ulong
val utoi : pari_ulong -> _ t
val utoineg : pari_ulong -> _ t
val utoipos : pari_ulong -> _ t
val utor : pari_ulong -> Signed.long -> _ t
val uutoi : pari_ulong -> pari_ulong -> _ t
val uutoineg : pari_ulong -> pari_ulong -> _ t
val vali : _ t -> Signed.long
val varncmp : Signed.long -> Signed.long -> int
val varnmax : Signed.long -> Signed.long -> Signed.long
val varnmin : Signed.long -> Signed.long -> Signed.long
val pari_err_component : string -> string -> _ t -> _ t -> unit
val pari_err_dim : string -> unit
val pari_err_domain : string -> string -> string -> _ t -> _ t -> unit
val pari_err_file : string -> string -> unit
val pari_err_filedesc : string -> Signed.long -> unit
val pari_err_flag : string -> unit
val pari_err_impl : string -> unit
val pari_err_inv : string -> _ t -> unit
val pari_err_irredpol : string -> _ t -> unit
val pari_err_maxprime : pari_ulong -> unit
val pari_err_modulus : string -> _ t -> _ t -> unit
val pari_err_op : string -> _ t -> _ t -> unit
val pari_err_overflow : string -> unit
val pari_err_package : string -> unit
val pari_err_prec : string -> unit
val pari_err_prime : string -> _ t -> unit
val pari_err_priority : string -> _ t -> string -> Signed.long -> unit
val pari_err_sqrtn : string -> _ t -> unit
val pari_err_type : string -> _ t -> unit
val pari_err_type2 : string -> _ t -> _ t -> unit
val pari_err_var : string -> _ t -> _ t -> unit
val pari_err_roots0 : string -> unit
val mkintmod : _ t -> _ t -> _ t
val mkintmodu : pari_ulong -> pari_ulong -> _ t
val mkpolmod : _ t -> _ t -> _ t
val mkfrac : _ t -> _ t -> _ t
val mkfracss : Signed.long -> Signed.long -> _ t

val qtoss :
  _ t -> Signed.long Ctypes_static.ptr -> Signed.long Ctypes_static.ptr -> unit

val sstoq : Signed.long -> Signed.long -> _ t
val uutoq : pari_ulong -> pari_ulong -> _ t
val mkfraccopy : _ t -> _ t -> _ t
val mkrfrac : _ t -> _ t -> _ t
val mkrfraccopy : _ t -> _ t -> _ t
val mkcomplex : _ t -> _ t -> _ t
val gen_i : unit -> _ t
val cgetc : Signed.long -> _ t
val mkquad : _ t -> _ t -> _ t -> _ t
val mkvecsmall : Signed.long -> _ t
val mkvecsmall2 : Signed.long -> Signed.long -> _ t
val mkvecsmall3 : Signed.long -> Signed.long -> Signed.long -> _ t

val mkvecsmall4 :
  Signed.long -> Signed.long -> Signed.long -> Signed.long -> _ t

val mkvecsmall5 :
  Signed.long -> Signed.long -> Signed.long -> Signed.long -> Signed.long -> _ t

val mkqfb : _ t -> _ t -> _ t -> _ t -> _ t
val mkvec : _ t -> _ t
val mkvec2 : _ t -> _ t -> _ t
val mkvec3 : _ t -> _ t -> _ t -> _ t
val mkvec4 : _ t -> _ t -> _ t -> _ t -> _ t
val mkvec5 : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val mkvecs : Signed.long -> _ t
val mkvec2s : Signed.long -> Signed.long -> _ t
val mkvec3s : Signed.long -> Signed.long -> Signed.long -> _ t
val mkvec4s : Signed.long -> Signed.long -> Signed.long -> Signed.long -> _ t
val mkveccopy : _ t -> _ t
val mkvec2copy : _ t -> _ t -> _ t
val mkcol : _ t -> _ t
val mkcol2 : _ t -> _ t -> _ t
val mkcol3 : _ t -> _ t -> _ t -> _ t
val mkcol4 : _ t -> _ t -> _ t -> _ t -> _ t
val mkcol5 : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val mkcol6 : _ t -> _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val mkcols : Signed.long -> _ t
val mkcol2s : Signed.long -> Signed.long -> _ t
val mkcol3s : Signed.long -> Signed.long -> Signed.long -> _ t
val mkcol4s : Signed.long -> Signed.long -> Signed.long -> Signed.long -> _ t
val mkcolcopy : _ t -> _ t
val mkmat : _ t -> _ t
val mkmat2 : _ t -> _ t -> _ t
val mkmat3 : _ t -> _ t -> _ t -> _ t
val mkmat4 : _ t -> _ t -> _ t -> _ t -> _ t
val mkmat5 : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val mkmatcopy : _ t -> _ t
val mkerr : Signed.long -> _ t
val mkoo : unit -> _ t
val mkmoo : unit -> _ t
val inf_get_sign : _ t -> Signed.long
val mkmat22s : Signed.long -> Signed.long -> Signed.long -> Signed.long -> _ t
val mkmat22 : _ t -> _ t -> _ t -> _ t -> _ t
val pol_x : Signed.long -> _ t
val pol_xn : Signed.long -> Signed.long -> _ t
val pol_xnall : Signed.long -> Signed.long -> _ t
val polxn_flx : Signed.long -> Signed.long -> _ t
val pol_1 : Signed.long -> _ t
val pol_0 : Signed.long -> _ t
val const_vec : Signed.long -> _ t -> _ t
val const_col : Signed.long -> _ t -> _ t
val const_vecsmall : Signed.long -> Signed.long -> _ t
val zeropadic : _ t -> Signed.long -> _ t
val zeropadic_shallow : _ t -> Signed.long -> _ t
val zeroser : Signed.long -> Signed.long -> _ t
val ser_isexactzero : _ t -> int
val zeropol : Signed.long -> _ t
val zerocol : Signed.long -> _ t
val zerovec : Signed.long -> _ t
val zeromat : Signed.long -> Signed.long -> _ t
val zero_flx : Signed.long -> _ t
val zero_flv : Signed.long -> _ t
val zero_flm : Signed.long -> Signed.long -> _ t
val zero_flm_copy : Signed.long -> Signed.long -> _ t
val zero_f2v : Signed.long -> _ t
val zero_f2m : Signed.long -> Signed.long -> _ t
val zero_f2m_copy : Signed.long -> Signed.long -> _ t
val zeromatcopy : Signed.long -> Signed.long -> _ t
val zerovec_block : Signed.long -> _ t
val col_ei : Signed.long -> Signed.long -> _ t
val vec_ei : Signed.long -> Signed.long -> _ t
val f2v_ei : Signed.long -> Signed.long -> _ t
val vecsmall_ei : Signed.long -> Signed.long -> _ t
val rg_col_ei : _ t -> Signed.long -> Signed.long -> _ t
val shallowcopy : _ t -> _ t
val vectrunc_init : Signed.long -> _ t
val coltrunc_init : Signed.long -> _ t
val lg_increase : _ t -> unit
val vectrunc_append : _ t -> _ t -> unit
val vectrunc_append_batch : _ t -> _ t -> unit
val vecsmalltrunc_init : Signed.long -> _ t
val vecsmalltrunc_append : _ t -> Signed.long -> unit
val hash_str : string -> pari_ulong
val hash_str_len : string -> Signed.long -> pari_ulong
val vec_shorten : _ t -> Signed.long -> _ t
val vec_lengthen : _ t -> Signed.long -> _ t
val vec_append : _ t -> _ t -> _ t
val vec_prepend : _ t -> _ t -> _ t
val vec_setconst : _ t -> _ t -> _ t
val vecsmall_shorten : _ t -> Signed.long -> _ t
val vecsmall_lengthen : _ t -> Signed.long -> _ t
val vec_to_vecsmall : _ t -> _ t
val vecsmall_to_vec : _ t -> _ t
val vecsmall_to_vec_inplace : _ t -> _ t
val vecsmall_to_col : _ t -> _ t
val vecsmall_lexcmp : _ t -> _ t -> int
val vecsmall_prefixcmp : _ t -> _ t -> int
val vecsmall_prepend : _ t -> Signed.long -> _ t
val vecsmall_append : _ t -> Signed.long -> _ t
val vecsmall_concat : _ t -> _ t -> _ t
val vecsmall_coincidence : _ t -> _ t -> Signed.long
val vecsmall_isin : _ t -> Signed.long -> Signed.long
val vecsmall_pack : _ t -> Signed.long -> Signed.long -> Signed.long
val vecsmall_indexmax : _ t -> Signed.long
val vecsmall_max : _ t -> Signed.long
val vecsmall_indexmin : _ t -> Signed.long
val vecsmall_min : _ t -> Signed.long
val zv_isscalar : _ t -> int
val qv_isscalar : _ t -> int
val rgv_isscalar : _ t -> int
val rgx_isscalar : _ t -> int
val rgx_equal_var : _ t -> _ t -> Signed.long
val rgx_to_rgv : _ t -> Signed.long -> _ t
val rgx_is_rational : _ t -> int
val rgx_is_zx : _ t -> int
val rgx_is_qx : _ t -> int
val rgx_is_monomial : _ t -> int
val rgv_is_zv : _ t -> int
val rgv_is_qv : _ t -> int
val rgv_isin_i : _ t -> _ t -> Signed.long -> Signed.long
val rgv_isin : _ t -> _ t -> Signed.long
val vecslice : _ t -> Signed.long -> Signed.long -> _ t
val vecslicepermute : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val rowslicepermute : _ t -> _ t -> Signed.long -> Signed.long -> _ t
val rowslice : _ t -> Signed.long -> Signed.long -> _ t

val matslice :
  _ t -> Signed.long -> Signed.long -> Signed.long -> Signed.long -> _ t

val rowsplice : _ t -> Signed.long -> _ t
val vecsplice : _ t -> Signed.long -> _ t
val rgm_minor : _ t -> Signed.long -> Signed.long -> _ t
val row : _ t -> Signed.long -> _ t
val flm_row : _ t -> Signed.long -> _ t
val rowcopy : _ t -> Signed.long -> _ t
val row_i : _ t -> Signed.long -> Signed.long -> Signed.long -> _ t
val vecreverse : _ t -> _ t
val vecsmall_reverse : _ t -> _ t
val vecreverse_inplace : _ t -> unit
val vecsmallpermute : _ t -> _ t -> _ t
val vecpermute : _ t -> _ t -> _ t
val rowpermute : _ t -> _ t -> _ t
val identity_zv : Signed.long -> _ t
val identity_perm : Signed.long -> _ t
val cyclic_perm : Signed.long -> Signed.long -> _ t
val perm_mul : _ t -> _ t -> _ t
val perm_sqr : _ t -> _ t
val perm_inv : _ t -> _ t
val perm_conj : _ t -> _ t -> _ t
val pari_free : unit Ctypes_static.ptr -> unit
val pari_malloc : int -> unit Ctypes_static.ptr
val pari_realloc : unit Ctypes_static.ptr -> int -> unit Ctypes_static.ptr
val pari_realloc_ip : unit Ctypes_static.ptr Ctypes_static.ptr -> int -> unit
val pari_calloc : int -> unit Ctypes_static.ptr
val cgetalloc : int -> Signed.long -> _ t
val icopy_avma : _ t -> pari_ulong -> _ t
val leafcopy_avma : _ t -> pari_ulong -> _ t
val gerepileuptoleaf : pari_ulong -> _ t -> _ t
val gerepileuptoint : pari_ulong -> _ t -> _ t
val gerepileupto : pari_ulong -> _ t -> _ t
val gerepilecopy : pari_ulong -> _ t -> _ t
val gunclonenull : _ t -> unit
val gunclonenull_deep : _ t -> unit

val gerepilemany :
  pari_ulong -> _ t Ctypes_static.ptr Ctypes_static.ptr -> int -> unit

val gerepileall : pari_ulong -> int -> unit
val gc_all : pari_ulong -> int -> _ t
val gerepilecoeffs : pari_ulong -> _ t -> int -> unit
val bin_copy : genbin Ctypes.structure Ctypes_static.ptr -> _ t
val genbinbase : genbin Ctypes.structure Ctypes_static.ptr -> _ t
val cgiv : _ t -> unit
val killblock : _ t -> unit
val is_universal_constant : _ t -> int
val cxcompotor : _ t -> Signed.long -> _ t
val cxtofp : _ t -> Signed.long -> _ t
val cxtoreal : _ t -> _ t
val gtodouble : _ t -> float
val gisdouble : _ t -> float Ctypes_static.ptr -> int
val gtos : _ t -> Signed.long
val gtou : _ t -> pari_ulong
val absfrac : _ t -> _ t
val absfrac_shallow : _ t -> _ t
val q_abs : _ t -> _ t
val q_abs_shallow : _ t -> _ t
val r_abs_shallow : _ t -> _ t
val r_abs : _ t -> _ t
val gtofp : _ t -> Signed.long -> _ t
val gtomp : _ t -> Signed.long -> _ t
val rgx_gtofp : _ t -> Signed.long -> _ t
val rgc_gtofp : _ t -> Signed.long -> _ t
val rgv_gtofp : _ t -> Signed.long -> _ t
val rgm_gtofp : _ t -> Signed.long -> _ t
val rgc_gtomp : _ t -> Signed.long -> _ t
val rgm_gtomp : _ t -> Signed.long -> _ t
val rgx_fpnorml2 : _ t -> Signed.long -> _ t
val rgc_fpnorml2 : _ t -> Signed.long -> _ t
val rgm_fpnorml2 : _ t -> Signed.long -> _ t
val affgr : _ t -> _ t -> unit
val affc_fixlg : _ t -> _ t -> _ t
val trunc_safe : _ t -> _ t
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
val bit_prec : _ t -> Signed.long
val bit_accuracy : Signed.long -> Signed.long
val prec2ndec : Signed.long -> Signed.long
val nbits2ndec : Signed.long -> Signed.long
val precdbl : Signed.long -> Signed.long
val divsbil : Signed.long -> Signed.long
val remsbil : Signed.long -> Signed.long
val fp_red : _ t -> _ t -> _ t
val fp_add : _ t -> _ t -> _ t -> _ t
val fp_sub : _ t -> _ t -> _ t -> _ t
val fp_neg : _ t -> _ t -> _ t
val fp_halve : _ t -> _ t -> _ t
val fp_center : _ t -> _ t -> _ t -> _ t
val fp_center_i : _ t -> _ t -> _ t -> _ t
val fp_addmul : _ t -> _ t -> _ t -> _ t -> _ t
val fp_mul : _ t -> _ t -> _ t -> _ t
val fp_sqr : _ t -> _ t -> _ t
val fp_mulu : _ t -> pari_ulong -> _ t -> _ t
val fp_muls : _ t -> Signed.long -> _ t -> _ t
val fp_inv : _ t -> _ t -> _ t
val fp_invsafe : _ t -> _ t -> _ t
val fp_div : _ t -> _ t -> _ t -> _ t
val fp_divu : _ t -> pari_ulong -> _ t -> _ t
val flx_mulu : _ t -> pari_ulong -> pari_ulong -> _ t
val get_f2x_mod : _ t -> _ t
val get_f2x_var : _ t -> Signed.long
val get_f2x_degree : _ t -> Signed.long
val get_f2xqx_mod : _ t -> _ t
val get_f2xqx_var : _ t -> Signed.long
val get_f2xqx_degree : _ t -> Signed.long
val get_flx_mod : _ t -> _ t
val get_flx_var : _ t -> Signed.long
val get_flx_degree : _ t -> Signed.long
val get_flxqx_mod : _ t -> _ t
val get_flxqx_var : _ t -> Signed.long
val get_flxqx_degree : _ t -> Signed.long
val get_fpx_mod : _ t -> _ t
val get_fpx_var : _ t -> Signed.long
val get_fpx_degree : _ t -> Signed.long
val get_fpxqx_mod : _ t -> _ t
val get_fpxqx_var : _ t -> Signed.long
val get_fpxqx_degree : _ t -> Signed.long
val submulii : _ t -> _ t -> _ t -> _ t
val mulsubii : _ t -> _ t -> _ t -> _ t
val submuliu : _ t -> _ t -> pari_ulong -> _ t
val addmuliu : _ t -> _ t -> pari_ulong -> _ t
val submuliu_inplace : _ t -> _ t -> pari_ulong -> _ t
val addmuliu_inplace : _ t -> _ t -> pari_ulong -> _ t
val lincombii : _ t -> _ t -> _ t -> _ t -> _ t
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
val qfb_is_qfi : _ t -> int
val sqrtr : _ t -> _ t
val cbrtr_abs : _ t -> _ t
val cbrtr : _ t -> _ t
val sqrtnr : _ t -> Signed.long -> _ t
val logint : _ t -> _ t -> Signed.long
val ulogint : pari_ulong -> pari_ulong -> pari_ulong
val ismpzero : _ t -> int
val isintzero : _ t -> int
val isint1 : _ t -> int
val isintm1 : _ t -> int
val equali1 : _ t -> int
val equalim1 : _ t -> int
val is_pm1 : _ t -> int
val is_bigint : _ t -> int
val odd : Signed.long -> int
val both_odd : Signed.long -> Signed.long -> int
val isonstack : _ t -> int
val dbllog2r : _ t -> float
val mul_content : _ t -> _ t -> _ t
val inv_content : _ t -> _ t
val div_content : _ t -> _ t -> _ t
val mul_denom : _ t -> _ t -> _ t
val constant_coeff : _ t -> _ t
val leading_coeff : _ t -> _ t
val flx_lead : _ t -> pari_ulong
val flx_constant : _ t -> pari_ulong
val degpol : _ t -> Signed.long
val lgpol : _ t -> Signed.long
val lgcols : _ t -> Signed.long
val nbrows : _ t -> Signed.long
val truecoef : _ t -> Signed.long -> _ t
val zxq_mul : _ t -> _ t -> _ t -> _ t
val zxq_sqr : _ t -> _ t -> _ t
val rgx_copy : _ t -> _ t
val rgx_coeff : _ t -> Signed.long -> _ t
val rgx_renormalize : _ t -> _ t
val rgx_div : _ t -> _ t -> _ t
val rgxqx_div : _ t -> _ t -> _ t -> _ t
val rgxqx_rem : _ t -> _ t -> _ t -> _ t
val fpx_div : _ t -> _ t -> _ t -> _ t
val flx_div : _ t -> _ t -> pari_ulong -> _ t
val flx_div_pre : _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val f2x_div : _ t -> _ t -> _ t
val fpv_fpc_mul : _ t -> _ t -> _ t -> _ t
val pol0_flx : Signed.long -> _ t
val pol1_flx : Signed.long -> _ t
val polx_flx : Signed.long -> _ t
val zero_zx : Signed.long -> _ t
val polx_zx : Signed.long -> _ t
val zx_shift : _ t -> Signed.long -> _ t
val zero_f2x : Signed.long -> _ t
val pol0_f2x : Signed.long -> _ t
val pol1_f2x : Signed.long -> _ t
val polx_f2x : Signed.long -> _ t
val f2x_equal1 : _ t -> int
val f2x_equal : _ t -> _ t -> int
val f2x_copy : _ t -> _ t
val f2v_copy : _ t -> _ t
val flv_copy : _ t -> _ t
val flx_copy : _ t -> _ t
val vecsmall_copy : _ t -> _ t
val flx_equal1 : _ t -> int
val zx_equal1 : _ t -> int
val zx_is_monic : _ t -> int
val zx_renormalize : _ t -> Signed.long -> _ t
val fpx_renormalize : _ t -> Signed.long -> _ t
val fpxx_renormalize : _ t -> Signed.long -> _ t
val fpxqx_renormalize : _ t -> Signed.long -> _ t
val f2x_renormalize : _ t -> Signed.long -> _ t
val f2xx_shift : _ t -> Signed.long -> Signed.long -> _ t
val f2v_to_f2x : _ t -> Signed.long -> _ t
val sturm : _ t -> Signed.long
val gval : _ t -> Signed.long -> Signed.long
val rgx_shift_inplace_init : Signed.long -> unit
val rgx_shift_inplace : _ t -> Signed.long -> _ t
val zc_to_zc : _ t -> _ t
val zv_to_zv : _ t -> _ t
val zx_to_zv : _ t -> Signed.long -> _ t
val zv_to_zx : _ t -> Signed.long -> _ t
val zm_to_zxv : _ t -> Signed.long -> _ t
val zero_zm : Signed.long -> Signed.long -> _ t
val zero_zv : Signed.long -> _ t
val zm_transpose : _ t -> _ t
val zm_copy : _ t -> _ t
val zv_copy : _ t -> _ t
val zm_row : _ t -> Signed.long -> _ t
val zc_hnfrem : _ t -> _ t -> _ t
val zm_hnfrem : _ t -> _ t -> _ t
val zm_lll : _ t -> float -> Signed.long -> _ t

val rgm_dimensions :
  _ t -> Signed.long Ctypes_static.ptr -> Signed.long Ctypes_static.ptr -> unit

val rgm_shallowcopy : _ t -> _ t
val f2m_copy : _ t -> _ t
val f3m_copy : _ t -> _ t
val flm_copy : _ t -> _ t
val zv_dvd : _ t -> _ t -> int
val zm_zv_mod : _ t -> _ t -> _ t
val zv_zv_mod : _ t -> _ t -> _ t
val vecmodii : _ t -> _ t -> _ t
val vecmoduu : _ t -> _ t -> _ t
val fq_red : _ t -> _ t -> _ t -> _ t
val fq_to_fpxq : _ t -> _ t -> _ t -> _ t
val rg_to_fq : _ t -> _ t -> _ t -> _ t
val gener_fq_local : _ t -> _ t -> _ t -> _ t
val random_fq : _ t -> _ t -> _ t
val fpxqx_div : _ t -> _ t -> _ t -> _ t -> _ t
val flxqx_div : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxqx_div_pre : _ t -> _ t -> _ t -> pari_ulong -> pari_ulong -> _ t
val f2xqx_div : _ t -> _ t -> _ t -> _ t
val fpxy_fq_evaly : _ t -> _ t -> _ t -> _ t -> Signed.long -> _ t
val fqx_red : _ t -> _ t -> _ t -> _ t
val fqx_add : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_neg : _ t -> _ t -> _ t -> _ t
val fqx_sub : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_fp_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_fq_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_mul : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_mulu : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fqx_sqr : _ t -> _ t -> _ t -> _ t
val fqx_powu : _ t -> pari_ulong -> _ t -> _ t -> _ t
val fqx_halve : _ t -> _ t -> _ t -> _ t
val fqx_div : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_get_red : _ t -> _ t -> _ t -> _ t
val fqx_rem : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_divrem : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fqx_div_by_x_x : _ t -> _ t -> _ t -> _ t -> _ t Ctypes_static.ptr -> _ t
val fqx_halfgcd : _ t -> _ t -> _ t -> _ t -> _ t
val fqx_gcd : _ t -> _ t -> _ t -> _ t -> _ t

val fqx_extgcd :
  _ t ->
  _ t ->
  _ t ->
  _ t ->
  _ t Ctypes_static.ptr ->
  _ t Ctypes_static.ptr ->
  _ t

val fqx_normalize : _ t -> _ t -> _ t -> _ t
val fqx_deriv : _ t -> _ t -> _ t -> _ t
val fqx_integ : _ t -> _ t -> _ t -> _ t
val fqx_factor : _ t -> _ t -> _ t -> _ t
val fqx_factor_squarefree : _ t -> _ t -> _ t -> _ t
val fqx_ddf : _ t -> _ t -> _ t -> _ t
val fqx_degfact : _ t -> _ t -> _ t -> _ t
val fqx_roots : _ t -> _ t -> _ t -> _ t
val fqx_to_mod : _ t -> _ t -> _ t -> _ t
val fqxq_add : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fqxq_sub : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fqxq_div : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fqxq_inv : _ t -> _ t -> _ t -> _ t -> _ t
val fqxq_invsafe : _ t -> _ t -> _ t -> _ t -> _ t
val fqxq_mul : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fqxq_sqr : _ t -> _ t -> _ t -> _ t -> _ t
val fqxq_pow : _ t -> _ t -> _ t -> _ t -> _ t -> _ t
val fqxn_expint : _ t -> Signed.long -> _ t -> _ t -> _ t
val fqxn_exp : _ t -> Signed.long -> _ t -> _ t -> _ t
val fqxn_inv : _ t -> Signed.long -> _ t -> _ t -> _ t
val fqxn_mul : _ t -> _ t -> Signed.long -> _ t -> _ t -> _ t
val fqxn_sqr : _ t -> Signed.long -> _ t -> _ t -> _ t
val fpxq_add : _ t -> _ t -> _ t -> _ t -> _ t
val fpxq_sub : _ t -> _ t -> _ t -> _ t -> _ t
val flxq_add : _ t -> _ t -> _ t -> pari_ulong -> _ t
val flxq_sub : _ t -> _ t -> _ t -> pari_ulong -> _ t
val f2x_coeff : _ t -> Signed.long -> pari_ulong
val f2x_clear : _ t -> Signed.long -> unit
val f2x_set : _ t -> Signed.long -> unit
val f2x_flip : _ t -> Signed.long -> unit
val f2v_coeff : _ t -> Signed.long -> pari_ulong
val f2v_clear : _ t -> Signed.long -> unit
val f2v_set : _ t -> Signed.long -> unit
val f2v_flip : _ t -> Signed.long -> unit
val f2m_coeff : _ t -> Signed.long -> Signed.long -> pari_ulong
val f2m_clear : _ t -> Signed.long -> Signed.long -> unit
val f2m_set : _ t -> Signed.long -> Signed.long -> unit
val f2m_flip : _ t -> Signed.long -> Signed.long -> unit
val f3m_coeff : _ t -> Signed.long -> Signed.long -> pari_ulong
val f3m_set : _ t -> Signed.long -> Signed.long -> pari_ulong -> unit
val matpascal : Signed.long -> _ t
val z_issquare : _ t -> Signed.long
val z_ispower : _ t -> pari_ulong -> Signed.long
val sqrti : _ t -> _ t
val gaddgs : _ t -> Signed.long -> _ t
val gcmpgs : _ t -> Signed.long -> int
val gequalgs : _ t -> Signed.long -> int
val gmaxsg : Signed.long -> _ t -> _ t
val gminsg : Signed.long -> _ t -> _ t
val gmulgs : _ t -> Signed.long -> _ t
val gmulgu : _ t -> pari_ulong -> _ t
val gsubgs : _ t -> Signed.long -> _ t
val gdivsg : Signed.long -> _ t -> _ t
val gmax_shallow : _ t -> _ t -> _ t
val gmin_shallow : _ t -> _ t -> _ t
val cxnorm : _ t -> _ t
val quadnorm : _ t -> _ t
val quad_disc : _ t -> _ t
val qfb_disc3 : _ t -> _ t -> _ t -> _ t
val qfb_disc : _ t -> _ t
val sqrfrac : _ t -> _ t
val normalize_frac : _ t -> unit
val powii : _ t -> _ t -> _ t
val powis : Signed.long -> _ t
val mpexpz : _ t -> _ t -> unit
val mplogz : _ t -> _ t -> unit
val mpcosz : _ t -> _ t -> unit
val mpsinz : _ t -> _ t -> unit
val gnegz : _ t -> _ t -> unit
val gabsz : _ t -> Signed.long -> _ t -> unit
val gaddz : _ t -> _ t -> _ t -> unit
val gsubz : _ t -> _ t -> _ t -> unit
val gmulz : _ t -> _ t -> _ t -> unit
val gdivz : _ t -> _ t -> _ t -> unit
val gdiventz : _ t -> _ t -> _ t -> unit
val gmodz : _ t -> _ t -> _ t -> unit
val gmul2nz : _ t -> Signed.long -> _ t -> unit
val gshiftz : _ t -> Signed.long -> _ t -> unit
val ell_get_a1 : _ t -> _ t
val ell_get_a2 : _ t -> _ t
val ell_get_a3 : _ t -> _ t
val ell_get_a4 : _ t -> _ t
val ell_get_a6 : _ t -> _ t
val ell_get_b2 : _ t -> _ t
val ell_get_b4 : _ t -> _ t
val ell_get_b6 : _ t -> _ t
val ell_get_b8 : _ t -> _ t
val ell_get_c4 : _ t -> _ t
val ell_get_c6 : _ t -> _ t
val ell_get_disc : _ t -> _ t
val ell_get_j : _ t -> _ t
val ell_get_type : _ t -> Signed.long
val ellff_get_field : _ t -> _ t
val ellff_get_a4a6 : _ t -> _ t
val ellqp_get_zero : _ t -> _ t
val ellqp_get_prec : _ t -> Signed.long
val ellqp_get_p : _ t -> _ t
val ellr_get_prec : _ t -> Signed.long
val ellr_get_sign : _ t -> Signed.long
val ellnf_get_nf : _ t -> _ t
val ellnf_get_bnf : _ t -> _ t
val checkell_i : _ t -> int
val ell_is_inf : _ t -> int
val ellinf : unit -> _ t
val modpr_get_pr : _ t -> _ t
val modpr_get_p : _ t -> _ t
val modpr_get_t : _ t -> _ t
val pr_get_p : _ t -> _ t
val pr_get_gen : _ t -> _ t
val pr_get_e : _ t -> Signed.long
val pr_get_f : _ t -> Signed.long
val pr_get_tau : _ t -> _ t
val pr_is_inert : _ t -> int
val pr_norm : _ t -> _ t
val upr_norm : _ t -> pari_ulong
val nf_get_varn : _ t -> Signed.long
val nf_get_pol : _ t -> _ t
val nf_get_degree : _ t -> Signed.long
val nf_get_r1 : _ t -> Signed.long
val nf_get_r2 : _ t -> Signed.long
val nf_get_disc : _ t -> _ t
val nf_get_index : _ t -> _ t
val nf_get_m : _ t -> _ t
val nf_get_g : _ t -> _ t
val nf_get_roundg : _ t -> _ t
val nf_get_tr : _ t -> _ t
val nf_get_diff : _ t -> _ t
val nf_get_ramified_primes : _ t -> _ t
val nf_get_roots : _ t -> _ t
val nf_get_zk : _ t -> _ t
val nf_get_zkprimpart : _ t -> _ t
val nf_get_zkden : _ t -> _ t
val nf_get_invzk : _ t -> _ t

val nf_get_sign :
  _ t -> Signed.long Ctypes_static.ptr -> Signed.long Ctypes_static.ptr -> unit

val cyc_get_expo : _ t -> _ t
val abgrp_get_no : _ t -> _ t
val abgrp_get_cyc : _ t -> _ t
val abgrp_get_gen : _ t -> _ t
val bnf_get_nf : _ t -> _ t
val bnf_get_clgp : _ t -> _ t
val bnf_get_no : _ t -> _ t
val bnf_get_cyc : _ t -> _ t
val bnf_get_gen : _ t -> _ t
val bnf_get_reg : _ t -> _ t
val bnf_get_logfu : _ t -> _ t
val bnf_get_sunits : _ t -> _ t
val bnf_get_tuu : _ t -> _ t
val bnf_get_tun : _ t -> Signed.long
val bnf_get_fu_nocheck : _ t -> _ t
val nfv_to_scalar_or_alg : _ t -> _ t -> _ t
val bnf_get_fu : _ t -> _ t
val bnr_get_bnf : _ t -> _ t
val bnr_get_bid : _ t -> _ t
val bnr_get_mod : _ t -> _ t
val bnr_get_nf : _ t -> _ t
val bnr_get_clgp : _ t -> _ t
val bnr_get_no : _ t -> _ t
val bnr_get_cyc : _ t -> _ t
val bnr_get_gen_nocheck : _ t -> _ t
val bnr_get_gen : _ t -> _ t
val locs_get_cyc : _ t -> _ t
val locs_get_lsprk : _ t -> _ t
val locs_get_lgenfil : _ t -> _ t
val locs_get_mod : _ t -> _ t
val locs_get_famod : _ t -> _ t
val locs_get_m_infty : _ t -> _ t
val gchar_get_basis : _ t -> _ t
val gchar_get_bnf : _ t -> _ t
val gchar_get_nf : _ t -> _ t
val gchar_get_zm : _ t -> _ t
val gchar_get_mod : _ t -> _ t
val gchar_get_modp : _ t -> _ t
val gchar_get_s : _ t -> _ t
val gchar_get_dldata : _ t -> _ t
val gchar_get_sfu : _ t -> _ t
val gchar_get_cyc : _ t -> _ t
val gchar_get_hnf : _ t -> _ t
val gchar_get_u : _ t -> _ t
val gchar_get_ui : _ t -> _ t
val gchar_get_m0 : _ t -> _ t
val gchar_get_u0 : _ t -> _ t
val gchar_get_r1 : _ t -> Signed.long
val gchar_get_r2 : _ t -> Signed.long
val gchar_get_loccyc : _ t -> _ t
val gchar_get_nc : _ t -> Signed.long
val gchar_get_ns : _ t -> Signed.long
val gchar_get_nm : _ t -> Signed.long
val gchar_get_evalprec : _ t -> Signed.long
val gchar_get_prec : _ t -> Signed.long
val gchar_get_nfprec : _ t -> Signed.long
val gchar_set_evalprec : _ t -> Signed.long -> unit
val gchar_set_prec : _ t -> Signed.long -> unit
val gchar_copy_precs : _ t -> _ t -> unit
val gchar_set_nfprec : _ t -> Signed.long -> unit
val gchar_get_ntors : _ t -> Signed.long
val gchar_get_nfree : _ t -> Signed.long
val gchar_get_nalg : _ t -> Signed.long
val gchar_set_basis : _ t -> _ t -> unit
val gchar_set_nf : _ t -> _ t -> unit
val gchar_set_ntors : _ t -> Signed.long -> unit
val gchar_set_nfree : _ t -> Signed.long -> unit
val gchar_set_nalg : _ t -> Signed.long -> unit
val gchar_set_cyc : _ t -> _ t -> unit
val gchar_set_huui : _ t -> _ t -> _ t -> _ t -> unit
val gchar_set_m0 : _ t -> _ t -> unit
val gchar_set_u0 : _ t -> _ t -> unit
val bid_get_mod : _ t -> _ t
val bid_get_ideal : _ t -> _ t
val bid_get_arch : _ t -> _ t
val bid_get_grp : _ t -> _ t
val bid_get_fact : _ t -> _ t
val bid_get_fact2 : _ t -> _ t
val bid_get_sprk : _ t -> _ t
val bid_get_sarch : _ t -> _ t
val bid_get_archp : _ t -> _ t
val bid_get_u : _ t -> _ t
val bid_get_no : _ t -> _ t
val bid_get_cyc : _ t -> _ t
val bid_get_gen_nocheck : _ t -> _ t
val bid_get_gen : _ t -> _ t
val znstar_get_n : _ t -> _ t
val znstar_get_fan : _ t -> _ t
val znstar_get_no : _ t -> _ t
val znstar_get_cyc : _ t -> _ t
val znstar_get_gen : _ t -> _ t
val znstar_get_conreycyc : _ t -> _ t
val znstar_get_conreygen : _ t -> _ t
val znstar_get_ui : _ t -> _ t
val znstar_get_u : _ t -> _ t
val znstar_get_pe : _ t -> _ t
val gal_get_pol : _ t -> _ t
val gal_get_p : _ t -> _ t
val gal_get_e : _ t -> _ t
val gal_get_mod : _ t -> _ t
val gal_get_roots : _ t -> _ t
val gal_get_invvdm : _ t -> _ t
val gal_get_den : _ t -> _ t
val gal_get_group : _ t -> _ t
val gal_get_gen : _ t -> _ t
val gal_get_orders : _ t -> _ t
val rnf_get_degree : _ t -> Signed.long
val rnf_get_nfdegree : _ t -> Signed.long
val rnf_get_absdegree : _ t -> Signed.long
val rnf_get_idealdisc : _ t -> _ t
val rnf_get_k : _ t -> _ t
val rnf_get_alpha : _ t -> _ t
val rnf_get_nf : _ t -> _ t
val rnf_get_nfzk : _ t -> _ t
val rnf_get_polabs : _ t -> _ t
val rnf_get_pol : _ t -> _ t
val rnf_get_disc : _ t -> _ t
val rnf_get_index : _ t -> _ t
val rnf_get_ramified_primes : _ t -> _ t
val rnf_get_varn : _ t -> Signed.long
val rnf_get_nfpol : _ t -> _ t
val rnf_get_nfvarn : _ t -> Signed.long
val rnf_get_zk : _ t -> _ t
val rnf_get_map : _ t -> _ t
val rnf_get_invzk : _ t -> _ t
val idealred : _ t -> _ t -> _ t
val idealchineseinit : _ t -> _ t -> _ t
val closure_arity : _ t -> Signed.long
val closure_is_variadic : _ t -> Signed.long
val closure_codestr : _ t -> string
val closure_get_code : _ t -> _ t
val closure_get_oper : _ t -> _ t
val closure_get_data : _ t -> _ t
val closure_get_dbg : _ t -> _ t
val closure_get_text : _ t -> _ t
val closure_get_frame : _ t -> _ t
val err_get_num : _ t -> Signed.long
val err_get_compo : _ t -> Signed.long -> _ t
val pari_err_bug : string -> unit
val pari_err_constpol : string -> unit
val pari_err_coprime : string -> _ t -> _ t -> unit
