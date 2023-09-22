type pari_ulong = Unsigned.ULong.t

val pari_ulong : pari_ulong Ctypes.typ

type gen = Signed.Long.t Ctypes.ptr

val gen : gen Ctypes.typ

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
val forprime_t_bb : (gen, forprime_t Ctypes.structure) Ctypes.field
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
val forprime_t_pp : (gen, forprime_t Ctypes.structure) Ctypes.field
val forcomposite_t : forcomposite_t Ctypes.structure Ctypes.typ
val forcomposite_t_first : (int, forcomposite_t Ctypes.structure) Ctypes.field
val forcomposite_t_b : (gen, forcomposite_t Ctypes.structure) Ctypes.field
val forcomposite_t_n : (gen, forcomposite_t Ctypes.structure) Ctypes.field
val forcomposite_t_p : (gen, forcomposite_t Ctypes.structure) Ctypes.field

val forcomposite_t_T :
  (forprime_t Ctypes.structure, forcomposite_t Ctypes.structure) Ctypes.field

val forvec_t : forvec_t Ctypes.structure Ctypes.typ
val forvec_t_first : (Signed.long, forvec_t Ctypes.structure) Ctypes.field
val forvec_t_a : (gen Ctypes_static.ptr, forvec_t Ctypes.structure) Ctypes.field
val forvec_t_m : (gen Ctypes_static.ptr, forvec_t Ctypes.structure) Ctypes.field
val forvec_t_M : (gen Ctypes_static.ptr, forvec_t Ctypes.structure) Ctypes.field
val forvec_t_n : (Signed.long, forvec_t Ctypes.structure) Ctypes.field

val forvec_t_next :
  ( (forvec_t Ctypes.structure Ctypes_static.ptr -> gen)
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
val forpart_t_v : (gen, forpart_t Ctypes.structure) Ctypes.field
val forperm_t : forperm_t Ctypes.structure Ctypes.typ
val forperm_t_k : (Signed.long, forperm_t Ctypes.structure) Ctypes.field
val forperm_t_first : (Signed.long, forperm_t Ctypes.structure) Ctypes.field
val forperm_t_v : (gen, forperm_t Ctypes.structure) Ctypes.field
val forsubset_t : forsubset_t Ctypes.structure Ctypes.typ
val forsubset_t_n : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_k : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_all : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_first : (Signed.long, forsubset_t Ctypes.structure) Ctypes.field
val forsubset_t_v : (gen, forsubset_t Ctypes.structure) Ctypes.field
val pari_plot : pari_plot Ctypes.structure Ctypes.typ

val pari_plot_draw :
  ( (pari_plot Ctypes.structure Ctypes_static.ptr -> gen -> gen -> gen -> unit)
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
val genbin_x : (gen, genbin Ctypes.structure) Ctypes.field
val genbin_base : (gen, genbin Ctypes.structure) Ctypes.field

val genbin_rebase :
  ( (gen -> Signed.long -> unit) Ctypes_static.static_funptr,
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
  (gen, pari_parsestate Ctypes.structure) Ctypes.field

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
  (gen, pari_global_state Ctypes.structure) Ctypes.field

val pari_global_state_seadata :
  (gen, pari_global_state Ctypes.structure) Ctypes.field

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

val pari_thread_data : (gen, pari_thread Ctypes.structure) Ctypes.field
val mt_state : mt_state Ctypes.structure Ctypes.typ
val mt_state_worker : (gen, mt_state Ctypes.structure) Ctypes.field
val mt_state_pending : (gen, mt_state Ctypes.structure) Ctypes.field
val mt_state_workid : (Signed.long, mt_state Ctypes.structure) Ctypes.field
val pari_mt : pari_mt Ctypes.structure Ctypes.typ

val pari_mt_mt :
  (mt_state Ctypes.structure, pari_mt Ctypes.structure) Ctypes.field

val pari_mt_get :
  ( (mt_state Ctypes.structure Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    Signed.long Ctypes_static.ptr ->
    gen)
    Ctypes_static.static_funptr,
    pari_mt Ctypes.structure )
  Ctypes.field

val pari_mt_submit :
  ( (mt_state Ctypes.structure Ctypes_static.ptr -> Signed.long -> gen -> unit)
    Ctypes_static.static_funptr,
    pari_mt Ctypes.structure )
  Ctypes.field

val pari_mt_end :
  (unit Ctypes_static.static_funptr, pari_mt Ctypes.structure) Ctypes.field

val parfor_iter : parfor_iter Ctypes.structure Ctypes.typ

val parfor_iter_pending :
  (Signed.long, parfor_iter Ctypes.structure) Ctypes.field

val parfor_iter_worker : (gen, parfor_iter Ctypes.structure) Ctypes.field

val parfor_iter_pt :
  (pari_mt Ctypes.structure, parfor_iter Ctypes.structure) Ctypes.field

val parfor_t : parfor_t Ctypes.structure Ctypes.typ
val parfor_t_a : (gen, parfor_t Ctypes.structure) Ctypes.field
val parfor_t_b : (gen, parfor_t Ctypes.structure) Ctypes.field

val parfor_t_iter :
  (parfor_iter Ctypes.structure, parfor_t Ctypes.structure) Ctypes.field

val parforeach_t : parforeach_t Ctypes.structure Ctypes.typ
val parforeach_t_x : (gen, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_W : (gen, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_i : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field
val parforeach_t_l : (Signed.long, parforeach_t Ctypes.structure) Ctypes.field

val parforeach_t_iter :
  (parfor_iter Ctypes.structure, parforeach_t Ctypes.structure) Ctypes.field

val parforprime_t : parforprime_t Ctypes.structure Ctypes.typ
val parforprime_t_v : (gen, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_forprime :
  (forprime_t Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforprime_t_iter :
  (parfor_iter Ctypes.structure, parforprime_t Ctypes.structure) Ctypes.field

val parforvec_t : parforvec_t Ctypes.structure Ctypes.typ
val parforvec_t_v : (gen, parforvec_t Ctypes.structure) Ctypes.field

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
val nfmaxord_t_T : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dT : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_T0 : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_unscale : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dK : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_index : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_basis : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_r1 : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_basden : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dTP : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dTE : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dKP : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_dKE : (gen, nfmaxord_t Ctypes.structure) Ctypes.field
val nfmaxord_t_certify : (Signed.long, nfmaxord_t Ctypes.structure) Ctypes.field
val qfr_data : qfr_data Ctypes.structure Ctypes.typ
val qfr_data_D : (gen, qfr_data Ctypes.structure) Ctypes.field
val qfr_data_sqrtD : (gen, qfr_data Ctypes.structure) Ctypes.field
val qfr_data_isqrtD : (gen, qfr_data Ctypes.structure) Ctypes.field
val fp_chk_fun : fp_chk_fun Ctypes.structure Ctypes.typ

val fp_chk_fun_f :
  ( (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_init :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr -> gen -> gen -> gen)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_f_post :
  ( (fp_chk_fun Ctypes.structure Ctypes_static.ptr -> gen -> gen -> gen)
    Ctypes_static.static_funptr,
    fp_chk_fun Ctypes.structure )
  Ctypes.field

val fp_chk_fun_data :
  (unit Ctypes_static.ptr, fp_chk_fun Ctypes.structure) Ctypes.field

val fp_chk_fun_skipfirst :
  (Signed.long, fp_chk_fun Ctypes.structure) Ctypes.field

val zlog_s : zlog_s Ctypes.structure Ctypes.typ
val zlog_s_bid : (gen, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_P : (gen, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_k : (gen, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_sprk : (gen, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_archp : (gen, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_mod : (gen, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_U : (gen, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_hU : (Signed.long, zlog_s Ctypes.structure) Ctypes.field
val zlog_s_no2 : (int, zlog_s Ctypes.structure) Ctypes.field
val bb_group : bb_group Ctypes.structure Ctypes.typ

val bb_group_mul :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_pow :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_rand :
  ( (unit Ctypes_static.ptr -> gen) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_hash :
  ( (gen -> pari_ulong) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal :
  ( (gen -> gen -> int) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_equal1 :
  ( (gen -> int) Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_group_easylog :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen -> gen)
    Ctypes_static.static_funptr,
    bb_group Ctypes.structure )
  Ctypes.field

val bb_field : bb_field Ctypes.structure Ctypes.typ

val bb_field_red :
  ( (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_add :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_mul :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_neg :
  ( (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_inv :
  ( (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_equal0 :
  ( (gen -> int) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_field_s :
  ( (unit Ctypes_static.ptr -> Signed.long -> gen) Ctypes_static.static_funptr,
    bb_field Ctypes.structure )
  Ctypes.field

val bb_algebra : bb_algebra Ctypes.structure Ctypes.typ

val bb_algebra_red :
  ( (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_add :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sub :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_mul :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_sqr :
  ( (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_one :
  ( (unit Ctypes_static.ptr -> gen) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_algebra_zero :
  ( (unit Ctypes_static.ptr -> gen) Ctypes_static.static_funptr,
    bb_algebra Ctypes.structure )
  Ctypes.field

val bb_ring : bb_ring Ctypes.structure Ctypes.typ

val bb_ring_add :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_mul :
  ( (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val bb_ring_sqr :
  ( (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr,
    bb_ring Ctypes.structure )
  Ctypes.field

val buchimag : gen -> gen -> gen -> gen -> gen
val buchreal : gen -> gen -> gen -> gen -> gen -> Signed.long -> gen
val zidealstar : gen -> gen -> gen
val zidealstarinit : gen -> gen -> gen
val zidealstarinitgen : gen -> gen -> gen
val factmod : gen -> gen -> gen
val mpbern : Signed.long -> Signed.long -> unit
val simplefactmod : gen -> gen -> gen
val listkill : gen -> unit
val isprincipalforce : gen -> gen -> gen
val isprincipalgen : gen -> gen -> gen
val isprincipalgenforce : gen -> gen -> gen
val f2ms_ker : gen -> Signed.long -> gen
val f2ms_to_f2m : gen -> Signed.long -> gen
val f2c_to_zc : gen -> gen
val f2c_to_mod : gen -> gen
val f2m_f2c_gauss : gen -> gen -> gen
val f2m_f2c_invimage : gen -> gen -> gen
val f2m_f2c_mul : gen -> gen -> gen
val f2m_deplin : gen -> gen
val f2m_det : gen -> pari_ulong
val f2m_det_sp : gen -> pari_ulong
val f2m_gauss : gen -> gen -> gen
val f2m_inv : gen -> gen
val f2m_invimage : gen -> gen -> gen
val f2m_ker : gen -> gen
val f2m_ker_sp : gen -> Signed.long -> gen
val f2m_mul : gen -> gen -> gen
val f2m_powu : gen -> pari_ulong -> gen
val f2m_rank : gen -> Signed.long
val f2m_row : gen -> Signed.long -> gen
val f2m_rowslice : gen -> Signed.long -> Signed.long -> gen
val f2m_to_f2ms : gen -> gen
val f2m_to_flm : gen -> gen
val f2m_to_zm : gen -> gen
val f2m_to_mod : gen -> gen
val f2m_transpose : gen -> gen
val f2v_add_inplace : gen -> gen -> unit
val f2v_and_inplace : gen -> gen -> unit
val f2v_dotproduct : gen -> gen -> pari_ulong
val f2v_equal0 : gen -> int
val f2v_hamming : gen -> pari_ulong
val f2v_negimply_inplace : gen -> gen -> unit
val f2v_or_inplace : gen -> gen -> unit
val f2v_slice : gen -> Signed.long -> Signed.long -> gen
val f2v_subset : gen -> gen -> int
val f2v_to_flv : gen -> gen
val matid_f2m : Signed.long -> gen
val f2x_f2xq_eval : gen -> gen -> gen -> gen
val f2x_f2xqv_eval : gen -> gen -> gen -> gen
val f2x_frobenius : gen -> gen
val f2x_1_add : gen -> gen
val f2x_add : gen -> gen -> gen
val f2x_deflate : gen -> Signed.long -> gen
val f2x_degfact : gen -> gen
val f2x_degree : gen -> Signed.long
val f2x_deriv : gen -> gen
val f2x_divrem : gen -> gen -> gen Ctypes_static.ptr -> gen
val f2x_eval : gen -> pari_ulong -> pari_ulong
val f2x_even_odd : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit

val f2x_extgcd :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val f2x_gcd : gen -> gen -> gen
val f2x_get_red : gen -> gen
val f2x_halfgcd : gen -> gen -> gen
val f2x_issquare : gen -> int
val f2x_matfrobenius : gen -> gen
val f2x_mul : gen -> gen -> gen
val f2x_recip : gen -> gen
val f2x_rem : gen -> gen -> gen
val f2x_shift : gen -> Signed.long -> gen
val f2x_sqr : gen -> gen
val f2x_sqrt : gen -> gen
val f2x_to_f2v : gen -> Signed.long -> gen
val f2x_to_f2xx : gen -> Signed.long -> gen
val f2x_to_flx : gen -> gen
val f2x_to_zx : gen -> gen
val f2x_valrem : gen -> gen Ctypes_static.ptr -> Signed.long
val f2xc_to_flxc : gen -> gen
val f2xc_to_zxc : gen -> gen
val f2xv_to_f2m : gen -> Signed.long -> gen
val f2xv_to_flxv_inplace : gen -> unit
val f2xv_to_zxv_inplace : gen -> unit
val f2xx_f2x_add : gen -> gen -> gen
val f2xx_f2x_mul : gen -> gen -> gen
val f2xx_add : gen -> gen -> gen
val f2xx_deriv : gen -> gen
val f2xx_renormalize : gen -> Signed.long -> gen
val f2xx_to_kronecker : gen -> Signed.long -> gen
val f2xx_to_flxx : gen -> gen
val f2xx_to_zxx : gen -> gen
val f2xx_to_f2xc : gen -> Signed.long -> Signed.long -> gen
val f2xxv_to_f2xm : gen -> Signed.long -> Signed.long -> gen
val f2xxc_to_zxxc : gen -> gen
val f2xy_f2xq_evalx : gen -> gen -> gen -> gen
val f2xy_f2xqv_evalx : gen -> gen -> gen -> gen
val f2xy_degreex : gen -> Signed.long
val f2xn_div : gen -> gen -> Signed.long -> gen
val f2xn_inv : gen -> Signed.long -> gen
val f2xn_red : gen -> Signed.long -> gen
val f2xq_artin_schreier : gen -> gen -> gen
val f2xq_autpow : gen -> Signed.long -> gen -> gen
val f2xq_conjvec : gen -> gen -> gen
val f2xq_div : gen -> gen -> gen -> gen
val f2xq_inv : gen -> gen -> gen
val f2xq_invsafe : gen -> gen -> gen
val f2xq_log : gen -> gen -> gen -> gen -> gen
val f2xq_matrix_pow : gen -> Signed.long -> Signed.long -> gen -> gen
val f2xq_mul : gen -> gen -> gen -> gen
val f2xq_order : gen -> gen -> gen -> gen
val f2xq_pow : gen -> gen -> gen -> gen
val f2xq_pow_init : gen -> gen -> Signed.long -> gen -> gen
val f2xq_pow_table : gen -> gen -> gen -> gen
val f2xq_powu : gen -> pari_ulong -> gen -> gen
val f2xq_powers : gen -> Signed.long -> gen -> gen
val f2xq_sqr : gen -> gen -> gen
val f2xq_sqrt : gen -> gen -> gen
val f2xq_sqrt_fast : gen -> gen -> gen -> gen
val f2xq_sqrtn : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val f2xq_trace : gen -> gen -> pari_ulong
val f2xqx_f2xq_mul : gen -> gen -> gen -> gen
val f2xqx_f2xq_mul_to_monic : gen -> gen -> gen -> gen
val f2xqx_f2xqxq_eval : gen -> gen -> gen -> gen -> gen
val f2xqx_f2xqxqv_eval : gen -> gen -> gen -> gen -> gen
val f2xqx_disc : gen -> gen -> gen
val f2xqx_divrem : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen

val f2xqx_extgcd :
  gen -> gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val f2xqx_gcd : gen -> gen -> gen -> gen
val f2xqx_get_red : gen -> gen -> gen
val f2xqx_halfgcd : gen -> gen -> gen -> gen

val f2xqx_halfgcd_all :
  gen -> gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val f2xqx_invbarrett : gen -> gen -> gen

val f2xqx_ispower :
  gen -> Signed.long -> gen -> gen Ctypes_static.ptr -> Signed.long

val f2xqx_mul : gen -> gen -> gen -> gen
val f2xqx_normalize : gen -> gen -> gen
val f2xqx_powu : gen -> pari_ulong -> gen -> gen
val f2xqx_red : gen -> gen -> gen
val f2xqx_rem : gen -> gen -> gen -> gen
val f2xqx_resultant : gen -> gen -> gen -> gen
val f2xqx_sqr : gen -> gen -> gen
val f2xqxq_inv : gen -> gen -> gen -> gen
val f2xqxq_invsafe : gen -> gen -> gen -> gen
val f2xqxq_mul : gen -> gen -> gen -> gen -> gen
val f2xqxq_sqr : gen -> gen -> gen -> gen
val f2xqxq_pow : gen -> gen -> gen -> gen -> gen
val f2xqxq_powers : gen -> Signed.long -> gen -> gen -> gen
val f2xqxq_autpow : gen -> Signed.long -> gen -> gen -> gen
val f2xqxq_auttrace : gen -> Signed.long -> gen -> gen -> gen
val f2xqxqv_red : gen -> gen -> gen -> gen
val flm_to_f2m : gen -> gen
val flv_to_f2v : gen -> gen
val flx_to_f2x : gen -> gen
val flxc_to_f2xc : gen -> gen
val flxx_to_f2xx : gen -> gen
val flxxc_to_f2xxc : gen -> gen
val kronecker_to_f2xqx : gen -> gen -> gen
val rg_to_f2xq : gen -> gen -> gen
val rgm_to_f2m : gen -> gen
val rgv_to_f2v : gen -> gen
val rgx_to_f2x : gen -> gen
val z_to_f2x : gen -> Signed.long -> gen
val zm_to_f2m : gen -> gen
val zv_to_f2v : gen -> gen
val zx_to_f2x : gen -> gen
val zxx_to_f2xx : gen -> Signed.long -> gen
val const_f2v : Signed.long -> gen
val gener_f2xq : gen -> gen Ctypes_static.ptr -> gen

val get_f2xq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  bb_field Ctypes.structure Ctypes_static.ptr

val monomial_f2x : Signed.long -> Signed.long -> gen
val pol1_f2xx : Signed.long -> Signed.long -> gen
val polx_f2xx : Signed.long -> Signed.long -> gen
val random_f2xqx : Signed.long -> Signed.long -> gen -> gen
val f2x_teichmuller : gen -> Signed.long -> gen
val f2xq_ellcard : gen -> gen -> gen -> gen
val f2xq_ellgens : gen -> gen -> gen -> gen -> gen -> gen -> gen
val f2xq_ellgroup : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen

val f2xq_elltwist :
  gen -> gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit

val f2xqe_add : gen -> gen -> gen -> gen -> gen
val f2xqe_changepoint : gen -> gen -> gen -> gen
val f2xqe_changepointinv : gen -> gen -> gen -> gen
val f2xqe_dbl : gen -> gen -> gen -> gen
val f2xqe_log : gen -> gen -> gen -> gen -> gen -> gen
val f2xqe_mul : gen -> gen -> gen -> gen -> gen
val f2xqe_neg : gen -> gen -> gen -> gen
val f2xqe_order : gen -> gen -> gen -> gen -> gen
val f2xqe_sub : gen -> gen -> gen -> gen -> gen
val f2xqe_tatepairing : gen -> gen -> gen -> gen -> gen -> gen
val f2xqe_weilpairing : gen -> gen -> gen -> gen -> gen -> gen

val get_f2xqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  gen ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_f2xqe : gen -> gen -> gen
val random_f2xqe : gen -> gen -> gen -> gen
val f3c_to_mod : gen -> gen
val f3c_to_zc : gen -> gen
val f3m_ker : gen -> gen
val f3m_ker_sp : gen -> Signed.long -> gen
val f3m_mul : gen -> gen -> gen
val f3m_row : gen -> Signed.long -> gen
val f3m_to_flm : gen -> gen
val f3m_to_zm : gen -> gen
val f3m_to_mod : gen -> gen
val f3m_transpose : gen -> gen
val f3v_to_flv : gen -> gen
val f3v_coeff : gen -> Signed.long -> pari_ulong
val f3v_clear : gen -> Signed.long -> unit
val f3v_set : gen -> Signed.long -> pari_ulong -> unit
val flm_to_f3m : gen -> gen
val flv_to_f3v : gen -> gen
val rgm_to_f3m : gen -> gen
val rgv_to_f3v : gen -> gen
val zm_to_f3m : gen -> gen
val zv_to_f3v : gen -> gen
val zero_f3m_copy : Signed.long -> Signed.long -> gen
val zero_f3v : Signed.long -> gen
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
  pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> gen

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

val fle_add : gen -> gen -> pari_ulong -> pari_ulong -> gen
val fle_dbl : gen -> pari_ulong -> pari_ulong -> gen
val fle_changepoint : gen -> gen -> pari_ulong -> gen
val fle_changepointinv : gen -> gen -> pari_ulong -> gen
val fle_log : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val fle_mul : gen -> gen -> pari_ulong -> pari_ulong -> gen
val fle_mulu : gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val fle_order : gen -> gen -> pari_ulong -> pari_ulong -> gen
val fle_sub : gen -> gen -> pari_ulong -> pari_ulong -> gen

val fle_tatepairing :
  gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val fle_to_flj : gen -> gen

val fle_weilpairing :
  gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val flj_add_pre : gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val flj_changepointinv_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flj_dbl_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen

val flj_mulu_pre :
  gen -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> gen

val flj_neg : gen -> pari_ulong -> gen
val flj_to_fle : gen -> pari_ulong -> gen
val flj_to_fle_pre : gen -> pari_ulong -> pari_ulong -> gen

val fljv_factorback_pre :
  gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen

val random_fle : pari_ulong -> pari_ulong -> pari_ulong -> gen
val random_fle_pre : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val random_flj_pre : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val flc_to_zc : gen -> gen
val flc_to_zc_inplace : gen -> gen
val flm_flc_gauss : gen -> gen -> pari_ulong -> gen
val flm_flc_invimage : gen -> gen -> pari_ulong -> gen
val flm_adjoint : gen -> pari_ulong -> gen
val flm_deplin : gen -> pari_ulong -> gen
val flm_det : gen -> pari_ulong -> pari_ulong
val flm_det_sp : gen -> pari_ulong -> pari_ulong
val flm_gauss : gen -> gen -> pari_ulong -> gen
val flm_intersect : gen -> gen -> pari_ulong -> gen
val flm_intersect_i : gen -> gen -> pari_ulong -> gen
val flm_inv : gen -> pari_ulong -> gen
val flm_invimage : gen -> gen -> pari_ulong -> gen
val flm_ker : gen -> pari_ulong -> gen
val flm_ker_sp : gen -> pari_ulong -> Signed.long -> gen
val flm_rank : gen -> pari_ulong -> Signed.long
val flm_to_zm : gen -> gen
val flm_to_zm_inplace : gen -> gen
val flv_to_zv : gen -> gen
val fl_to_flx : pari_ulong -> Signed.long -> gen
val fl2_equal1 : gen -> int
val fl2_inv_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val fl2_mul_pre : gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val fl2_norm_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val fl2_pow_pre : gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val fl2_sqr_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val fl2_sqrt_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen

val fl2_sqrtn_pre :
  gen ->
  gen ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen

val flm_to_flxv : gen -> Signed.long -> gen
val flm_to_flxx : gen -> Signed.long -> Signed.long -> gen
val flv_flm_polint : gen -> gen -> pari_ulong -> Signed.long -> gen
val flv_inv : gen -> pari_ulong -> gen
val flv_inv_inplace : gen -> pari_ulong -> unit
val flv_inv_pre_inplace : gen -> pari_ulong -> pari_ulong -> unit
val flv_inv_pre : gen -> pari_ulong -> pari_ulong -> gen
val flv_invvandermonde : gen -> pari_ulong -> pari_ulong -> gen
val flv_polint : gen -> gen -> pari_ulong -> Signed.long -> gen
val flv_prod : gen -> pari_ulong -> pari_ulong
val flv_prod_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong
val flv_roots_to_pol : gen -> pari_ulong -> Signed.long -> gen
val flv_to_flx : gen -> Signed.long -> gen
val flx_fl_add : gen -> pari_ulong -> pari_ulong -> gen
val flx_fl_mul : gen -> pari_ulong -> pari_ulong -> gen
val flx_fl_mul_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val flx_fl_mul_to_monic : gen -> pari_ulong -> pari_ulong -> gen
val flx_fl_sub : gen -> pari_ulong -> pari_ulong -> gen

val flx_fl2_eval_pre :
  gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen

val flx_flv_multieval : gen -> gen -> pari_ulong -> gen
val flx_flxq_eval : gen -> gen -> gen -> pari_ulong -> gen
val flx_flxq_eval_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flx_flxqv_eval : gen -> gen -> gen -> pari_ulong -> gen
val flx_flxqv_eval_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flx_frobenius : gen -> pari_ulong -> gen
val flx_frobenius_pre : gen -> pari_ulong -> pari_ulong -> gen
val flx_laplace : gen -> pari_ulong -> gen
val flx_newton : gen -> Signed.long -> pari_ulong -> gen
val flx_add : gen -> gen -> pari_ulong -> gen
val flx_blocks : gen -> Signed.long -> Signed.long -> gen
val flx_composedprod : gen -> gen -> pari_ulong -> gen
val flx_composedsum : gen -> gen -> pari_ulong -> gen
val flx_convol : gen -> gen -> pari_ulong -> gen
val flx_deflate : gen -> Signed.long -> gen
val flx_deriv : gen -> pari_ulong -> gen
val flx_diff1 : gen -> pari_ulong -> gen
val flx_digits : gen -> gen -> pari_ulong -> gen

val flx_div_by_x_x :
  gen -> pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> gen

val flx_divrem : gen -> gen -> pari_ulong -> gen Ctypes_static.ptr -> gen

val flx_divrem_pre :
  gen -> gen -> pari_ulong -> pari_ulong -> gen Ctypes_static.ptr -> gen

val flx_double : gen -> pari_ulong -> gen
val flx_equal : gen -> gen -> int
val flx_eval : gen -> pari_ulong -> pari_ulong -> pari_ulong
val flx_eval_powers_pre : gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong
val flx_eval_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong

val flx_extgcd :
  gen ->
  gen ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val flx_extgcd_pre :
  gen ->
  gen ->
  pari_ulong ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val flx_extresultant :
  gen ->
  gen ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  pari_ulong

val flx_extresultant_pre :
  gen ->
  gen ->
  pari_ulong ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  pari_ulong

val flx_fromnewton : gen -> pari_ulong -> gen
val flx_gcd : gen -> gen -> pari_ulong -> gen
val flx_gcd_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flx_get_red : gen -> pari_ulong -> gen
val flx_get_red_pre : gen -> pari_ulong -> pari_ulong -> gen
val flx_halfgcd : gen -> gen -> pari_ulong -> gen

val flx_halfgcd_all :
  gen ->
  gen ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val flx_halfgcd_all_pre :
  gen ->
  gen ->
  pari_ulong ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val flx_halfgcd_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flx_halve : gen -> pari_ulong -> gen
val flx_inflate : gen -> Signed.long -> gen
val flx_integ : gen -> pari_ulong -> gen
val flx_invbarrett : gen -> pari_ulong -> gen
val flx_invlaplace : gen -> pari_ulong -> gen
val flx_is_squarefree : gen -> pari_ulong -> int
val flx_is_smooth : gen -> Signed.long -> pari_ulong -> int
val flx_is_smooth_pre : gen -> Signed.long -> pari_ulong -> pari_ulong -> int
val flx_matfrobenius : gen -> pari_ulong -> gen
val flx_matfrobenius_pre : gen -> pari_ulong -> pari_ulong -> gen
val flx_mod_xn1 : gen -> pari_ulong -> pari_ulong -> gen
val flx_mod_xnm1 : gen -> pari_ulong -> pari_ulong -> gen
val flx_mul : gen -> gen -> pari_ulong -> gen
val flx_mul_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flx_neg : gen -> pari_ulong -> gen
val flx_neg_inplace : gen -> pari_ulong -> gen
val flx_normalize : gen -> pari_ulong -> gen
val flx_powu : gen -> pari_ulong -> pari_ulong -> gen
val flx_powu_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val flx_recip : gen -> gen
val flx_red : gen -> pari_ulong -> gen
val flx_rem : gen -> gen -> pari_ulong -> gen
val flx_rem_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flx_renormalize : gen -> Signed.long -> gen
val flx_rescale : gen -> pari_ulong -> pari_ulong -> gen
val flx_resultant : gen -> gen -> pari_ulong -> pari_ulong
val flx_resultant_pre : gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong
val flx_shift : gen -> Signed.long -> gen
val flx_splitting : gen -> Signed.long -> gen
val flx_sqr : gen -> pari_ulong -> gen
val flx_sqr_pre : gen -> pari_ulong -> pari_ulong -> gen
val flx_sub : gen -> gen -> pari_ulong -> gen
val flx_translate1 : gen -> pari_ulong -> gen
val flx_translate1_basecase : gen -> pari_ulong -> gen
val flx_to_flv : gen -> Signed.long -> gen
val flx_to_flxx : gen -> Signed.long -> gen
val flx_to_zx : gen -> gen
val flx_to_zx_inplace : gen -> gen
val flx_triple : gen -> pari_ulong -> gen
val flx_val : gen -> Signed.long
val flx_valrem : gen -> gen Ctypes_static.ptr -> Signed.long
val flxc_flxqv_eval : gen -> gen -> gen -> pari_ulong -> gen
val flxc_flxqv_eval_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxc_flxq_eval : gen -> gen -> gen -> pari_ulong -> gen
val flxc_flxq_eval_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxc_eval_powers_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxc_neg : gen -> pari_ulong -> gen
val flxc_sub : gen -> gen -> pari_ulong -> gen
val flxc_to_zxc : gen -> gen
val flxm_flx_add_shallow : gen -> gen -> pari_ulong -> gen
val flxm_eval_powers_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxm_neg : gen -> pari_ulong -> gen
val flxm_sub : gen -> gen -> pari_ulong -> gen
val flxm_to_flxxv : gen -> Signed.long -> gen
val flxm_to_zxm : gen -> gen
val flxt_red : gen -> pari_ulong -> gen
val flxv_flc_mul : gen -> gen -> pari_ulong -> gen
val flxv_flv_multieval : gen -> gen -> pari_ulong -> gen
val flxv_flx_fromdigits : gen -> gen -> pari_ulong -> gen
val flxv_composedsum : gen -> pari_ulong -> gen
val flxv_prod : gen -> pari_ulong -> gen
val flxv_red : gen -> pari_ulong -> gen
val flxv_to_flm : gen -> Signed.long -> gen
val flxv_to_flxx : gen -> Signed.long -> gen
val flxv_to_zxv : gen -> gen
val flxv_to_zxv_inplace : gen -> unit
val flxn_div : gen -> gen -> Signed.long -> pari_ulong -> gen
val flxn_div_pre : gen -> gen -> Signed.long -> pari_ulong -> pari_ulong -> gen
val flxn_exp : gen -> Signed.long -> pari_ulong -> gen
val flxn_expint : gen -> Signed.long -> pari_ulong -> gen
val flxn_inv : gen -> Signed.long -> pari_ulong -> gen
val flxn_mul : gen -> gen -> Signed.long -> pari_ulong -> gen
val flxn_mul_pre : gen -> gen -> Signed.long -> pari_ulong -> pari_ulong -> gen
val flxn_sqr : gen -> Signed.long -> pari_ulong -> gen
val flxn_sqr_pre : gen -> Signed.long -> pari_ulong -> pari_ulong -> gen
val flxn_red : gen -> Signed.long -> gen
val flxq_autpow : gen -> pari_ulong -> gen -> pari_ulong -> gen

val flxq_autpow_pre :
  gen -> pari_ulong -> gen -> pari_ulong -> pari_ulong -> gen

val flxq_autpowers : gen -> pari_ulong -> gen -> pari_ulong -> gen
val flxq_autsum : gen -> pari_ulong -> gen -> pari_ulong -> gen
val flxq_auttrace : gen -> pari_ulong -> gen -> pari_ulong -> gen

val flxq_auttrace_pre :
  gen -> pari_ulong -> gen -> pari_ulong -> pari_ulong -> gen

val flxq_charpoly : gen -> gen -> pari_ulong -> gen
val flxq_conjvec : gen -> gen -> pari_ulong -> gen
val flxq_div : gen -> gen -> gen -> pari_ulong -> gen
val flxq_div_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_inv : gen -> gen -> pari_ulong -> gen
val flxq_inv_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_invsafe : gen -> gen -> pari_ulong -> gen
val flxq_invsafe_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_issquare : gen -> gen -> pari_ulong -> int
val flxq_is2npower : gen -> Signed.long -> gen -> pari_ulong -> int
val flxq_log : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxq_lroot : gen -> gen -> Signed.long -> gen
val flxq_lroot_pre : gen -> gen -> Signed.long -> pari_ulong -> gen
val flxq_lroot_fast : gen -> gen -> gen -> Signed.long -> gen
val flxq_lroot_fast_pre : gen -> gen -> gen -> Signed.long -> pari_ulong -> gen

val flxq_matrix_pow :
  gen -> Signed.long -> Signed.long -> gen -> pari_ulong -> gen

val flxq_matrix_pow_pre :
  gen -> Signed.long -> Signed.long -> gen -> pari_ulong -> pari_ulong -> gen

val flxq_minpoly : gen -> gen -> pari_ulong -> gen
val flxq_minpoly_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_mul : gen -> gen -> gen -> pari_ulong -> gen
val flxq_mul_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_norm : gen -> gen -> pari_ulong -> pari_ulong
val flxq_order : gen -> gen -> gen -> pari_ulong -> gen
val flxq_pow : gen -> gen -> gen -> pari_ulong -> gen
val flxq_pow_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_pow_init : gen -> gen -> Signed.long -> gen -> pari_ulong -> gen

val flxq_pow_init_pre :
  gen -> gen -> Signed.long -> gen -> pari_ulong -> pari_ulong -> gen

val flxq_pow_table_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_pow_table : gen -> gen -> gen -> pari_ulong -> gen
val flxq_powu : gen -> pari_ulong -> gen -> pari_ulong -> gen
val flxq_powu_pre : gen -> pari_ulong -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_powers : gen -> Signed.long -> gen -> pari_ulong -> gen

val flxq_powers_pre :
  gen -> Signed.long -> gen -> pari_ulong -> pari_ulong -> gen

val flxq_sqr : gen -> gen -> pari_ulong -> gen
val flxq_sqr_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_sqrt : gen -> gen -> pari_ulong -> gen
val flxq_sqrt_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxq_sqrtn : gen -> gen -> gen -> pari_ulong -> gen Ctypes_static.ptr -> gen
val flxq_trace : gen -> gen -> pari_ulong -> pari_ulong
val flxqc_flxq_mul : gen -> gen -> gen -> pari_ulong -> gen
val flxqm_flxq_mul : gen -> gen -> gen -> pari_ulong -> gen
val flxqv_dotproduct : gen -> gen -> gen -> pari_ulong -> gen
val flxqv_dotproduct_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val rg_to_f2 : gen -> pari_ulong
val rg_to_fl : gen -> pari_ulong -> pari_ulong
val rg_to_flxq : gen -> gen -> pari_ulong -> gen
val rgx_to_flx : gen -> pari_ulong -> gen
val rgxv_to_flxv : gen -> pari_ulong -> gen
val z_to_flx : gen -> pari_ulong -> Signed.long -> gen
val zxv_to_flxv : gen -> pari_ulong -> gen
val zxt_to_flxt : gen -> pari_ulong -> gen
val gener_flxq : gen -> pari_ulong -> gen Ctypes_static.ptr -> gen

val get_flxq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  pari_ulong ->
  bb_field Ctypes.structure Ctypes_static.ptr

val get_flxq_star :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  pari_ulong ->
  bb_group Ctypes.structure Ctypes_static.ptr

val monomial_flx : pari_ulong -> Signed.long -> Signed.long -> gen
val random_flx : Signed.long -> Signed.long -> pari_ulong -> gen
val zero_flxc : Signed.long -> Signed.long -> gen
val zero_flxm : Signed.long -> Signed.long -> Signed.long -> gen
val zlx_translate1 : gen -> pari_ulong -> Signed.long -> gen
val zx_to_flx : gen -> pari_ulong -> gen
val flxx_fl_mul : gen -> pari_ulong -> pari_ulong -> gen
val flxx_flx_add : gen -> gen -> pari_ulong -> gen
val flxx_flx_mul : gen -> gen -> pari_ulong -> gen
val flxx_flx_sub : gen -> gen -> pari_ulong -> gen
val flxx_laplace : gen -> pari_ulong -> gen
val flxx_add : gen -> gen -> pari_ulong -> gen
val flxx_blocks : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val flxx_deriv : gen -> pari_ulong -> gen
val flxx_double : gen -> pari_ulong -> gen
val flxx_invlaplace : gen -> pari_ulong -> gen
val flxx_neg : gen -> pari_ulong -> gen
val flxx_renormalize : gen -> Signed.long -> gen
val flxx_shift : gen -> Signed.long -> Signed.long -> gen
val flxx_sub : gen -> gen -> pari_ulong -> gen
val flxx_swap : gen -> Signed.long -> Signed.long -> gen
val flxx_to_flm : gen -> Signed.long -> gen
val flxx_to_flx : gen -> gen
val flxx_to_flxc : gen -> Signed.long -> Signed.long -> gen
val flxx_to_zxx : gen -> gen
val flxx_translate1 : gen -> Signed.long -> Signed.long -> gen
val flxx_triple : gen -> pari_ulong -> gen
val flxxc_sub : gen -> gen -> pari_ulong -> gen
val flxxc_to_zxxc : gen -> gen
val flxxm_to_zxxm : gen -> gen
val flxxv_to_flxm : gen -> Signed.long -> Signed.long -> gen
val flxxn_red : gen -> Signed.long -> gen
val flxy_flx_div : gen -> gen -> pari_ulong -> gen
val flxy_flx_translate : gen -> gen -> pari_ulong -> gen
val flxy_flxqv_evalx : gen -> gen -> gen -> pari_ulong -> gen
val flxy_flxqv_evalx_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxy_flxq_evalx : gen -> gen -> gen -> pari_ulong -> gen
val flxy_flxq_evalx_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxy_evalx : gen -> pari_ulong -> pari_ulong -> gen
val flxy_evalx_powers_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxy_evalx_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val flxyqq_pow : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqv_roots_to_pol : gen -> gen -> pari_ulong -> Signed.long -> gen
val flxqxc_flxqxqv_eval : gen -> gen -> gen -> gen -> pari_ulong -> gen

val flxqxc_flxqxqv_eval_pre :
  gen -> gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxc_flxqxq_eval : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqxq_autpow : gen -> Signed.long -> gen -> gen -> pari_ulong -> gen

val flxqxq_autpow_pre :
  gen -> Signed.long -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxq_autsum : gen -> Signed.long -> gen -> gen -> pari_ulong -> gen

val flxqxq_autsum_pre :
  gen -> Signed.long -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxq_auttrace : gen -> pari_ulong -> gen -> gen -> pari_ulong -> gen

val flxqxq_auttrace_pre :
  gen -> pari_ulong -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxq_div : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqxq_div_pre : gen -> gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqxq_inv : gen -> gen -> gen -> pari_ulong -> gen
val flxqxq_inv_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqxq_invsafe : gen -> gen -> gen -> pari_ulong -> gen
val flxqxq_invsafe_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxq_matrix_pow :
  gen -> Signed.long -> Signed.long -> gen -> gen -> pari_ulong -> gen

val flxqxq_minpoly : gen -> gen -> gen -> pari_ulong -> gen
val flxqxq_minpoly_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqxq_mul : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqxq_mul_pre : gen -> gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqxq_pow : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqxq_pow_pre : gen -> gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqxq_powers : gen -> Signed.long -> gen -> gen -> pari_ulong -> gen

val flxqxq_powers_pre :
  gen -> Signed.long -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxq_powu : gen -> pari_ulong -> gen -> gen -> pari_ulong -> gen

val flxqxq_powu_pre :
  gen -> pari_ulong -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxq_sqr : gen -> gen -> gen -> pari_ulong -> gen
val flxqxq_sqr_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqxv_prod : gen -> gen -> pari_ulong -> gen
val flxqx_flxqxqv_eval : gen -> gen -> gen -> gen -> pari_ulong -> gen

val flxqx_flxqxqv_eval_pre :
  gen -> gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqx_flxqxq_eval : gen -> gen -> gen -> gen -> pari_ulong -> gen

val flxqx_flxqxq_eval_pre :
  gen -> gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqx_flxq_mul : gen -> gen -> gen -> pari_ulong -> gen
val flxqx_flxq_mul_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_flxq_mul_to_monic : gen -> gen -> gen -> pari_ulong -> gen

val flxqx_flxq_mul_to_monic_pre :
  gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqx_newton : gen -> Signed.long -> gen -> pari_ulong -> gen

val flxqx_newton_pre :
  gen -> Signed.long -> gen -> pari_ulong -> pari_ulong -> gen

val flxqx_composedsum : gen -> gen -> gen -> pari_ulong -> gen
val flxqx_disc : gen -> gen -> pari_ulong -> gen

val flxqx_div_by_x_x :
  gen -> gen -> gen -> pari_ulong -> gen Ctypes_static.ptr -> gen

val flxqx_div_by_x_x_pre :
  gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen Ctypes_static.ptr -> gen

val flxqx_divrem :
  gen -> gen -> gen -> pari_ulong -> gen Ctypes_static.ptr -> gen

val flxqx_divrem_pre :
  gen -> gen -> gen -> pari_ulong -> Signed.long -> gen Ctypes_static.ptr -> gen

val flxqx_dotproduct : gen -> gen -> gen -> pari_ulong -> gen

val flxqx_extgcd :
  gen ->
  gen ->
  gen ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val flxqx_extgcd_pre :
  gen ->
  gen ->
  gen ->
  pari_ulong ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val flxqx_fromnewton : gen -> gen -> pari_ulong -> gen
val flxqx_fromnewton_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_gcd : gen -> gen -> gen -> pari_ulong -> gen
val flxqx_gcd_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_get_red : gen -> gen -> pari_ulong -> gen
val flxqx_get_red_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_halfgcd : gen -> gen -> gen -> pari_ulong -> gen

val flxqx_halfgcd_all :
  gen ->
  gen ->
  gen ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val flxqx_halfgcd_all_pre :
  gen ->
  gen ->
  gen ->
  pari_ulong ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val flxqx_halfgcd_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_invbarrett : gen -> gen -> pari_ulong -> gen
val flxqx_invbarrett_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_mul : gen -> gen -> gen -> pari_ulong -> gen
val flxqx_mul_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_normalize : gen -> gen -> pari_ulong -> gen
val flxqx_normalize_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_powu : gen -> pari_ulong -> gen -> pari_ulong -> gen
val flxqx_powu_pre : gen -> pari_ulong -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_red : gen -> gen -> pari_ulong -> gen
val flxqx_red_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_rem : gen -> gen -> gen -> pari_ulong -> gen
val flxqx_rem_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_resultant : gen -> gen -> gen -> pari_ulong -> gen
val flxqx_resultant_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_safegcd : gen -> gen -> gen -> pari_ulong -> gen
val flxqx_saferesultant : gen -> gen -> gen -> pari_ulong -> gen
val flxqx_sqr : gen -> gen -> pari_ulong -> gen
val flxqx_sqr_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqxn_expint : gen -> Signed.long -> gen -> pari_ulong -> gen

val flxqxn_expint_pre :
  gen -> Signed.long -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxn_inv : gen -> Signed.long -> gen -> pari_ulong -> gen

val flxqxn_inv_pre :
  gen -> Signed.long -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxn_mul : gen -> gen -> Signed.long -> gen -> pari_ulong -> gen

val flxqxn_mul_pre :
  gen -> gen -> Signed.long -> gen -> pari_ulong -> pari_ulong -> gen

val flxqxn_sqr : gen -> Signed.long -> gen -> pari_ulong -> gen

val flxqxn_sqr_pre :
  gen -> Signed.long -> gen -> pari_ulong -> pari_ulong -> gen

val flxy_degreex : gen -> Signed.long

val flxy_eval_powers_pre :
  gen -> gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong

val fly_to_flxy : gen -> Signed.long -> gen
val kronecker_to_flxqx : gen -> gen -> pari_ulong -> gen
val kronecker_to_flxqx_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val rgx_to_flxqx : gen -> gen -> pari_ulong -> gen

val get_flxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  pari_ulong ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val pol1_flxx : Signed.long -> Signed.long -> gen
val polx_flxx : Signed.long -> Signed.long -> gen
val random_flxqx : Signed.long -> Signed.long -> gen -> pari_ulong -> gen
val zlxx_translate1 : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val zxx_to_kronecker : gen -> gen -> gen
val flxq_ellcard : gen -> gen -> gen -> pari_ulong -> gen
val flxq_ellgens : gen -> gen -> gen -> gen -> gen -> gen -> pari_ulong -> gen

val flxq_ellgroup :
  gen -> gen -> gen -> gen -> pari_ulong -> gen Ctypes_static.ptr -> gen

val flxq_elltwist :
  gen ->
  gen ->
  gen ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  unit

val flxq_ellj : gen -> gen -> gen -> pari_ulong -> gen

val flxq_ellj_to_a4a6 :
  gen ->
  gen ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  unit

val flxqe_add : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqe_changepoint : gen -> gen -> gen -> pari_ulong -> gen
val flxqe_changepointinv : gen -> gen -> gen -> pari_ulong -> gen
val flxqe_dbl : gen -> gen -> gen -> pari_ulong -> gen
val flxqe_log : gen -> gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqe_mul : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqe_neg : gen -> gen -> pari_ulong -> gen
val flxqe_order : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqe_sub : gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqe_tatepairing : gen -> gen -> gen -> gen -> gen -> pari_ulong -> gen
val flxqe_weilpairing : gen -> gen -> gen -> gen -> gen -> pari_ulong -> gen

val flxqe_weilpairing_pre :
  gen -> gen -> gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen

val zxx_to_flxx : gen -> pari_ulong -> Signed.long -> gen
val zxxt_to_flxxt : gen -> pari_ulong -> Signed.long -> gen
val zxxv_to_flxxv : gen -> pari_ulong -> Signed.long -> gen

val get_flxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  gen ->
  pari_ulong ->
  bb_group Ctypes.structure Ctypes_static.ptr

val rge_to_flxqe : gen -> gen -> pari_ulong -> gen
val random_flxqe : gen -> gen -> gen -> pari_ulong -> gen
val polisclass : gen -> Signed.long
val fl_elltrace : pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fl_elltrace_cm :
  Signed.long -> pari_ulong -> pari_ulong -> pari_ulong -> Signed.long

val fp_ellcard : gen -> gen -> gen -> gen
val fp_elldivpol : gen -> gen -> Signed.long -> gen -> gen
val fp_ellgens : gen -> gen -> gen -> gen -> gen -> gen -> gen
val fp_ellgroup : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fp_ellj : gen -> gen -> gen -> gen

val fp_ellj_to_a4a6 :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit

val fp_elljissupersingular : gen -> gen -> int

val fp_elltwist :
  gen -> gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit

val fp_ffellcard : gen -> gen -> gen -> Signed.long -> gen -> gen
val fpe_add : gen -> gen -> gen -> gen -> gen
val fpe_changepoint : gen -> gen -> gen -> gen
val fpe_changepointinv : gen -> gen -> gen -> gen
val fpe_dbl : gen -> gen -> gen -> gen
val fpe_log : gen -> gen -> gen -> gen -> gen -> gen
val fpe_mul : gen -> gen -> gen -> gen -> gen
val fpe_neg : gen -> gen -> gen
val fpe_order : gen -> gen -> gen -> gen -> gen
val fpe_sub : gen -> gen -> gen -> gen -> gen
val fpe_to_fpj : gen -> gen
val fpe_to_mod : gen -> gen -> gen
val fpe_tatepairing : gen -> gen -> gen -> gen -> gen -> gen
val fpe_weilpairing : gen -> gen -> gen -> gen -> gen -> gen
val fpj_add : gen -> gen -> gen -> gen -> gen
val fpj_dbl : gen -> gen -> gen -> gen
val fpj_mul : gen -> gen -> gen -> gen -> gen
val fpj_neg : gen -> gen -> gen
val fpj_to_fpe : gen -> gen -> gen
val fpxq_ellcard : gen -> gen -> gen -> gen -> gen
val fpxq_ellcard_supersingular : gen -> gen -> gen -> gen -> gen
val fpxq_elldivpol : gen -> gen -> Signed.long -> gen -> gen -> gen
val fpxq_ellgens : gen -> gen -> gen -> gen -> gen -> gen -> gen -> gen

val fpxq_ellgroup :
  gen -> gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen

val fpxq_ellj : gen -> gen -> gen -> gen -> gen
val fpxq_elljissupersingular : gen -> gen -> gen -> int

val fpxq_elltwist :
  gen ->
  gen ->
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  unit

val fpxqe_add : gen -> gen -> gen -> gen -> gen -> gen
val fpxqe_changepoint : gen -> gen -> gen -> gen -> gen
val fpxqe_changepointinv : gen -> gen -> gen -> gen -> gen
val fpxqe_dbl : gen -> gen -> gen -> gen -> gen
val fpxqe_log : gen -> gen -> gen -> gen -> gen -> gen -> gen
val fpxqe_mul : gen -> gen -> gen -> gen -> gen -> gen
val fpxqe_neg : gen -> gen -> gen -> gen
val fpxqe_order : gen -> gen -> gen -> gen -> gen -> gen
val fpxqe_sub : gen -> gen -> gen -> gen -> gen -> gen
val fpxqe_tatepairing : gen -> gen -> gen -> gen -> gen -> gen -> gen
val fpxqe_weilpairing : gen -> gen -> gen -> gen -> gen -> gen -> gen
val fq_elljissupersingular : gen -> gen -> gen -> int
val fq_ellcard_supersingular : gen -> gen -> gen -> gen -> gen
val rge_to_fpe : gen -> gen -> gen
val rge_to_fpxqe : gen -> gen -> gen -> gen

val get_fpe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  gen ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpxqe_group :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  gen ->
  gen ->
  bb_group Ctypes.structure Ctypes_static.ptr

val ellsupersingularj_fpxq : gen -> gen -> gen
val elltrace_extension : gen -> Signed.long -> gen -> gen
val random_fpe : gen -> gen -> gen -> gen
val random_fpxqe : gen -> gen -> gen -> gen -> gen
val fp_issquare : gen -> gen -> int
val fp_fpx_sub : gen -> gen -> gen -> gen
val fp_fpxq_log : gen -> gen -> gen -> gen -> gen -> gen
val fpv_fpm_polint : gen -> gen -> gen -> Signed.long -> gen
val fpv_inv : gen -> gen -> gen
val fpv_invvandermonde : gen -> gen -> gen -> gen
val fpv_polint : gen -> gen -> gen -> Signed.long -> gen
val fpv_roots_to_pol : gen -> gen -> Signed.long -> gen
val fpx_fp_add : gen -> gen -> gen -> gen
val fpx_fp_add_shallow : gen -> gen -> gen -> gen
val fpx_fp_div : gen -> gen -> gen -> gen
val fpx_fp_mul : gen -> gen -> gen -> gen
val fpx_fp_mul_to_monic : gen -> gen -> gen -> gen
val fpx_fp_mulspec : gen -> gen -> gen -> Signed.long -> gen
val fpx_fp_sub : gen -> gen -> gen -> gen
val fpx_fp_sub_shallow : gen -> gen -> gen -> gen
val fpx_fpv_multieval : gen -> gen -> gen -> gen
val fpx_fpxq_eval : gen -> gen -> gen -> gen -> gen
val fpx_fpxqv_eval : gen -> gen -> gen -> gen -> gen
val fpx_fpxv_multirem : gen -> gen -> gen -> gen
val fpx_frobenius : gen -> gen -> gen
val fpx_laplace : gen -> gen -> gen
val fpx_newton : gen -> Signed.long -> gen -> gen
val fpx_add : gen -> gen -> gen -> gen
val fpx_center : gen -> gen -> gen -> gen
val fpx_center_i : gen -> gen -> gen -> gen
val fpx_chinese_coprime : gen -> gen -> gen -> gen -> gen -> gen -> gen
val fpx_composedprod : gen -> gen -> gen -> gen
val fpx_composedsum : gen -> gen -> gen -> gen
val fpx_convol : gen -> gen -> gen -> gen
val fpx_deriv : gen -> gen -> gen
val fpx_digits : gen -> gen -> gen -> gen
val fpx_disc : gen -> gen -> gen
val fpx_div_by_x_x : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fpx_divrem : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fpx_divu : gen -> pari_ulong -> gen -> gen
val fpx_dotproduct : gen -> gen -> gen -> gen
val fpx_eval : gen -> gen -> gen -> gen

val fpx_extgcd :
  gen -> gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val fpx_extresultant :
  gen -> gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val fpx_fromnewton : gen -> gen -> gen
val fpx_gcd : gen -> gen -> gen -> gen
val fpx_gcd_check : gen -> gen -> gen -> gen
val fpx_get_red : gen -> gen -> gen
val fpx_halve : gen -> gen -> gen
val fpx_halfgcd : gen -> gen -> gen -> gen

val fpx_halfgcd_all :
  gen -> gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val fpx_integ : gen -> gen -> gen
val fpx_invbarrett : gen -> gen -> gen
val fpx_invlaplace : gen -> gen -> gen
val fpx_is_squarefree : gen -> gen -> int
val fpx_matfrobenius : gen -> gen -> gen
val fpx_mul : gen -> gen -> gen -> gen
val fpx_mulspec : gen -> gen -> gen -> Signed.long -> Signed.long -> gen
val fpx_mulu : gen -> pari_ulong -> gen -> gen
val fpx_neg : gen -> gen -> gen
val fpx_normalize : gen -> gen -> gen
val fpx_powu : gen -> pari_ulong -> gen -> gen
val fpx_red : gen -> gen -> gen
val fpx_rem : gen -> gen -> gen -> gen
val fpx_rescale : gen -> gen -> gen -> gen
val fpx_resultant : gen -> gen -> gen -> gen
val fpx_sqr : gen -> gen -> gen
val fpx_sub : gen -> gen -> gen -> gen
val fpx_valrem : gen -> gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val fpxc_fpxq_eval : gen -> gen -> gen -> gen -> gen
val fpxc_fpxqv_eval : gen -> gen -> gen -> gen -> gen
val fpxm_fpxqv_eval : gen -> gen -> gen -> gen -> gen
val fpxq_autpow : gen -> pari_ulong -> gen -> gen -> gen
val fpxq_autpowers : gen -> Signed.long -> gen -> gen -> gen
val fpxq_autsum : gen -> pari_ulong -> gen -> gen -> gen
val fpxq_auttrace : gen -> pari_ulong -> gen -> gen -> gen
val fpxq_charpoly : gen -> gen -> gen -> gen
val fpxq_conjvec : gen -> gen -> gen -> gen
val fpxq_div : gen -> gen -> gen -> gen -> gen
val fpxq_inv : gen -> gen -> gen -> gen
val fpxq_invsafe : gen -> gen -> gen -> gen
val fpxq_issquare : gen -> gen -> gen -> int
val fpxq_log : gen -> gen -> gen -> gen -> gen -> gen
val fpxq_matrix_pow : gen -> Signed.long -> Signed.long -> gen -> gen -> gen
val fpxq_minpoly : gen -> gen -> gen -> gen
val fpxq_mul : gen -> gen -> gen -> gen -> gen
val fpxq_norm : gen -> gen -> gen -> gen
val fpxq_order : gen -> gen -> gen -> gen -> gen
val fpxq_pow : gen -> gen -> gen -> gen -> gen
val fpxq_powu : gen -> pari_ulong -> gen -> gen -> gen
val fpxq_powers : gen -> Signed.long -> gen -> gen -> gen
val fpxq_red : gen -> gen -> gen -> gen
val fpxq_sqr : gen -> gen -> gen -> gen
val fpxq_sqrt : gen -> gen -> gen -> gen
val fpxq_sqrtn : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fpxq_trace : gen -> gen -> gen -> gen
val fpxqc_to_mod : gen -> gen -> gen -> gen
val fpxqm_autsum : gen -> pari_ulong -> gen -> gen -> gen
val fpxt_red : gen -> gen -> gen
val fpxv_fpx_fromdigits : gen -> gen -> gen -> gen
val fpxv_chinese : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fpxv_composedsum : gen -> gen -> gen
val fpxv_factorback : gen -> gen -> gen -> Signed.long -> gen
val fpxv_prod : gen -> gen -> gen
val fpxv_red : gen -> gen -> gen
val fpxn_div : gen -> gen -> Signed.long -> gen -> gen
val fpxn_exp : gen -> Signed.long -> gen -> gen
val fpxn_expint : gen -> Signed.long -> gen -> gen
val fpxn_inv : gen -> Signed.long -> gen -> gen
val fpxn_mul : gen -> gen -> Signed.long -> gen -> gen
val fpxn_sqr : gen -> Signed.long -> gen -> gen
val fq_issquare : gen -> gen -> gen -> int
val fq_ispower : gen -> gen -> gen -> gen -> Signed.long
val fq_log : gen -> gen -> gen -> gen -> gen -> gen
val fqc_to_mod : gen -> gen -> gen -> gen
val fqm_to_mod : gen -> gen -> gen -> gen
val fqv_inv : gen -> gen -> gen -> gen
val z_to_fpx : gen -> gen -> Signed.long -> gen
val gener_fpxq : gen -> gen -> gen Ctypes_static.ptr -> gen
val gener_fpxq_local : gen -> gen -> gen -> gen

val get_fpxq_star :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  bb_group Ctypes.structure Ctypes_static.ptr

val get_fpx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpx : Signed.long -> Signed.long -> gen -> gen
val f2x_ddf : gen -> gen
val f2x_factor : gen -> gen
val f2x_factor_squarefree : gen -> gen
val f2x_is_irred : gen -> int
val flx_ddf : gen -> pari_ulong -> gen
val flx_ddf_pre : gen -> pari_ulong -> pari_ulong -> gen
val flx_is_irred : gen -> pari_ulong -> int
val flx_is_totally_split : gen -> pari_ulong -> int

val flx_ispower :
  gen -> pari_ulong -> pari_ulong -> gen Ctypes_static.ptr -> Signed.long

val flx_degfact : gen -> pari_ulong -> gen
val flx_factor : gen -> pari_ulong -> gen
val flx_factor_squarefree : gen -> pari_ulong -> gen
val flx_factor_squarefree_pre : gen -> pari_ulong -> pari_ulong -> gen
val flx_nbfact : gen -> pari_ulong -> Signed.long
val flx_nbfact_pre : gen -> pari_ulong -> pari_ulong -> Signed.long
val flx_nbfact_frobenius : gen -> gen -> pari_ulong -> Signed.long

val flx_nbfact_frobenius_pre :
  gen -> gen -> pari_ulong -> pari_ulong -> Signed.long

val flx_nbfact_by_degree :
  gen -> Signed.long Ctypes_static.ptr -> pari_ulong -> gen

val flx_nbroots : gen -> pari_ulong -> Signed.long
val flx_oneroot : gen -> pari_ulong -> pari_ulong
val flx_oneroot_split : gen -> pari_ulong -> pari_ulong
val flx_oneroot_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong
val flx_oneroot_split_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong
val flx_roots : gen -> pari_ulong -> gen
val flx_roots_pre : gen -> pari_ulong -> pari_ulong -> gen
val flx_rootsff : gen -> gen -> pari_ulong -> gen
val fpx_ddf : gen -> gen -> gen
val fpx_ddf_degree : gen -> gen -> gen -> Signed.long
val fpx_degfact : gen -> gen -> gen
val fpx_factor : gen -> gen -> gen
val fpx_factor_squarefree : gen -> gen -> gen
val fpx_is_irred : gen -> gen -> int
val fpx_is_totally_split : gen -> gen -> int

val fpx_ispower :
  gen -> pari_ulong -> gen -> gen Ctypes_static.ptr -> Signed.long

val fpx_nbfact : gen -> gen -> Signed.long
val fpx_nbfact_frobenius : gen -> gen -> gen -> Signed.long
val fpx_nbroots : gen -> gen -> Signed.long
val fpx_oneroot : gen -> gen -> gen
val fpx_oneroot_split : gen -> gen -> gen
val fpx_roots : gen -> gen -> gen
val fpx_roots_mult : gen -> Signed.long -> gen -> gen
val fpx_rootsff : gen -> gen -> gen -> gen
val fpx_split_part : gen -> gen -> gen
val f2xqx_ddf : gen -> gen -> gen
val f2xqx_degfact : gen -> gen -> gen
val f2xqx_factor : gen -> gen -> gen
val f2xqx_factor_squarefree : gen -> gen -> gen
val f2xqx_roots : gen -> gen -> gen
val flx_factorff_irred : gen -> gen -> pari_ulong -> gen

val flx_ffintersect :
  gen ->
  gen ->
  Signed.long ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen ->
  gen ->
  unit

val flx_ffisom : gen -> gen -> pari_ulong -> gen
val flxq_ffisom_inv : gen -> gen -> pari_ulong -> gen
val flxqx_frobenius : gen -> gen -> pari_ulong -> gen
val flxqx_frobenius_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flxqx_ddf : gen -> gen -> pari_ulong -> gen
val flxqx_ddf_degree : gen -> gen -> gen -> pari_ulong -> Signed.long
val flxqx_degfact : gen -> gen -> pari_ulong -> gen
val flxqx_factor : gen -> gen -> pari_ulong -> gen
val flxqx_factor_squarefree : gen -> gen -> pari_ulong -> gen
val flxqx_factor_squarefree_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen

val flxqx_ispower :
  gen -> pari_ulong -> gen -> pari_ulong -> gen Ctypes_static.ptr -> Signed.long

val flxqx_is_squarefree : gen -> gen -> pari_ulong -> Signed.long
val flxqx_nbfact : gen -> gen -> pari_ulong -> Signed.long
val flxqx_nbfact_frobenius : gen -> gen -> gen -> pari_ulong -> Signed.long

val flxqx_nbfact_by_degree :
  gen -> Signed.long Ctypes_static.ptr -> gen -> pari_ulong -> gen

val flxqx_nbroots : gen -> gen -> pari_ulong -> Signed.long
val flxqx_roots : gen -> gen -> pari_ulong -> gen
val flxqxq_halffrobenius : gen -> gen -> gen -> pari_ulong -> gen
val fpx_factorff : gen -> gen -> gen -> gen
val fpx_factorff_irred : gen -> gen -> gen -> gen

val fpx_ffintersect :
  gen ->
  gen ->
  Signed.long ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen ->
  gen ->
  unit

val fpx_ffisom : gen -> gen -> gen -> gen
val fpxq_ffisom_inv : gen -> gen -> gen -> gen
val fpxqx_frobenius : gen -> gen -> gen -> gen
val fpxqx_ddf : gen -> gen -> gen -> gen
val fpxqx_ddf_degree : gen -> gen -> gen -> gen -> Signed.long
val fpxqx_degfact : gen -> gen -> gen -> gen
val fpxqx_factor : gen -> gen -> gen -> gen
val fpxqx_factor_squarefree : gen -> gen -> gen -> gen

val fpxqx_ispower :
  gen -> pari_ulong -> gen -> gen -> gen Ctypes_static.ptr -> Signed.long

val fpxqx_nbfact : gen -> gen -> gen -> Signed.long
val fpxqx_nbfact_frobenius : gen -> gen -> gen -> gen -> Signed.long
val fpxqx_nbroots : gen -> gen -> gen -> Signed.long
val fpxqx_roots : gen -> gen -> gen -> gen
val fpxqx_split_part : gen -> gen -> gen -> gen
val fpxqxq_halffrobenius : gen -> gen -> gen -> gen -> gen
val fqx_is_squarefree : gen -> gen -> gen -> Signed.long

val fqx_ispower :
  gen -> pari_ulong -> gen -> gen -> gen Ctypes_static.ptr -> Signed.long

val fqx_nbfact : gen -> gen -> gen -> Signed.long
val fqx_nbroots : gen -> gen -> gen -> Signed.long
val factorff : gen -> gen -> gen -> gen
val factormod0 : gen -> gen -> Signed.long -> gen
val factormodddf : gen -> gen -> gen
val factormodsqf : gen -> gen -> gen

val ff_parse_tp :
  gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> Signed.long -> int

val polrootsff : gen -> gen -> gen -> gen
val polrootsmod : gen -> gen -> gen
val rootmod0 : gen -> gen -> Signed.long -> gen
val fpxqx_fpxq_mul : gen -> gen -> gen -> gen -> gen
val fpxqx_fpxqxqv_eval : gen -> gen -> gen -> gen -> gen -> gen
val fpxqx_fpxqxq_eval : gen -> gen -> gen -> gen -> gen -> gen
val fpxqx_digits : gen -> gen -> gen -> gen -> gen
val fpxqx_disc : gen -> gen -> gen -> gen
val fpxqx_div_by_x_x : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fpxqx_divrem : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fpxqx_dotproduct : gen -> gen -> gen -> gen -> gen

val fpxqx_extgcd :
  gen ->
  gen ->
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val fpxqx_gcd : gen -> gen -> gen -> gen -> gen
val fpxqx_get_red : gen -> gen -> gen -> gen
val fpxqx_halfgcd : gen -> gen -> gen -> gen -> gen

val fpxqx_halfgcd_all :
  gen ->
  gen ->
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val fpxqx_invbarrett : gen -> gen -> gen -> gen
val fpxqx_mul : gen -> gen -> gen -> gen -> gen
val fpxqx_powu : gen -> pari_ulong -> gen -> gen -> gen
val fpxqx_red : gen -> gen -> gen -> gen
val fpxqx_rem : gen -> gen -> gen -> gen -> gen
val fpxqx_resultant : gen -> gen -> gen -> gen -> gen
val fpxqx_sqr : gen -> gen -> gen -> gen
val fpxqx_to_mod : gen -> gen -> gen -> gen
val fpxqxq_div : gen -> gen -> gen -> gen -> gen -> gen
val fpxqxq_inv : gen -> gen -> gen -> gen -> gen
val fpxqxq_invsafe : gen -> gen -> gen -> gen -> gen

val fpxqxq_matrix_pow :
  gen -> Signed.long -> Signed.long -> gen -> gen -> gen -> gen

val fpxqxq_minpoly : gen -> gen -> gen -> gen -> gen
val fpxqxq_mul : gen -> gen -> gen -> gen -> gen -> gen
val fpxqxq_pow : gen -> gen -> gen -> gen -> gen -> gen
val fpxqxq_powers : gen -> Signed.long -> gen -> gen -> gen -> gen
val fpxqxq_sqr : gen -> gen -> gen -> gen -> gen
val fpxqxq_autpow : gen -> Signed.long -> gen -> gen -> gen -> gen
val fpxqxq_autsum : gen -> Signed.long -> gen -> gen -> gen -> gen
val fpxqxq_auttrace : gen -> Signed.long -> gen -> gen -> gen -> gen
val fpxqxt_red : gen -> gen -> gen -> gen
val fpxqxv_fpxqx_fromdigits : gen -> gen -> gen -> gen -> gen
val fpxqxv_prod : gen -> gen -> gen -> gen
val fpxqxv_red : gen -> gen -> gen -> gen
val fpxqxn_div : gen -> gen -> Signed.long -> gen -> gen -> gen
val fpxqxn_exp : gen -> Signed.long -> gen -> gen -> gen
val fpxqxn_expint : gen -> Signed.long -> gen -> gen -> gen
val fpxqxn_inv : gen -> Signed.long -> gen -> gen -> gen
val fpxqxn_mul : gen -> gen -> Signed.long -> gen -> gen -> gen
val fpxqxn_sqr : gen -> Signed.long -> gen -> gen -> gen
val fpxx_fp_mul : gen -> gen -> gen -> gen
val fpxx_fpx_mul : gen -> gen -> gen -> gen
val fpxx_add : gen -> gen -> gen -> gen
val fpxx_deriv : gen -> gen -> gen
val fpxx_halve : gen -> gen -> gen
val fpxx_integ : gen -> gen -> gen
val fpxx_mulu : gen -> pari_ulong -> gen -> gen
val fpxx_neg : gen -> gen -> gen
val fpxx_red : gen -> gen -> gen
val fpxx_sub : gen -> gen -> gen -> gen
val fpxy_fpxq_evalx : gen -> gen -> gen -> gen -> gen
val fpxy_fpxqv_evalx : gen -> gen -> gen -> gen -> gen
val fpxy_eval : gen -> gen -> gen -> gen -> gen
val fpxy_evalx : gen -> gen -> gen -> gen
val fpxy_evaly : gen -> gen -> gen -> Signed.long -> gen
val fpxyqq_pow : gen -> gen -> gen -> gen -> gen -> gen
val fqxc_to_mod : gen -> gen -> gen -> gen
val fqxm_to_mod : gen -> gen -> gen -> gen
val kronecker_to_fpxqx : gen -> gen -> gen -> gen

val get_fpxqx_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  Signed.long ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val get_fpxqxq_algebra :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  gen ->
  bb_algebra Ctypes.structure Ctypes_static.ptr

val random_fpxqx : Signed.long -> Signed.long -> gen -> gen -> gen
val flc_flv_mul : gen -> gen -> pari_ulong -> gen
val flc_to_mod : gen -> pari_ulong -> gen
val flm_fl_add : gen -> pari_ulong -> pari_ulong -> gen
val flm_fl_mul : gen -> pari_ulong -> pari_ulong -> gen
val flm_fl_mul_inplace : gen -> pari_ulong -> pari_ulong -> unit
val flm_fl_mul_pre : gen -> pari_ulong -> pari_ulong -> pari_ulong -> gen
val flm_fl_sub : gen -> pari_ulong -> pari_ulong -> gen
val flm_flc_mul : gen -> gen -> pari_ulong -> gen
val flm_flc_mul_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen

val flm_flc_mul_pre_flx :
  gen -> gen -> pari_ulong -> pari_ulong -> Signed.long -> gen

val flm_add : gen -> gen -> pari_ulong -> gen
val flm_center : gen -> pari_ulong -> pari_ulong -> gen
val flm_mul : gen -> gen -> pari_ulong -> gen
val flm_mul_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val flm_neg : gen -> pari_ulong -> gen
val flm_powers : gen -> pari_ulong -> pari_ulong -> gen
val flm_powu : gen -> pari_ulong -> pari_ulong -> gen
val flm_sub : gen -> gen -> pari_ulong -> gen
val flm_to_mod : gen -> pari_ulong -> gen
val flm_transpose : gen -> gen
val flv_fl_div : gen -> pari_ulong -> pari_ulong -> gen
val flv_fl_div_inplace : gen -> pari_ulong -> pari_ulong -> unit
val flv_fl_mul : gen -> pari_ulong -> pari_ulong -> gen
val flv_fl_mul_inplace : gen -> pari_ulong -> pari_ulong -> unit

val flv_fl_mul_part_inplace :
  gen -> pari_ulong -> pari_ulong -> Signed.long -> unit

val flv_add : gen -> gen -> pari_ulong -> gen
val flv_add_inplace : gen -> gen -> pari_ulong -> unit
val flv_center : gen -> pari_ulong -> pari_ulong -> gen
val flv_dotproduct : gen -> gen -> pari_ulong -> pari_ulong
val flv_dotproduct_pre : gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong
val flv_neg : gen -> pari_ulong -> gen
val flv_neg_inplace : gen -> pari_ulong -> unit
val flv_sub : gen -> gen -> pari_ulong -> gen
val flv_sub_inplace : gen -> gen -> pari_ulong -> unit
val flv_sum : gen -> pari_ulong -> pari_ulong
val flx_dotproduct : gen -> gen -> pari_ulong -> pari_ulong
val flx_dotproduct_pre : gen -> gen -> pari_ulong -> pari_ulong -> pari_ulong
val fp_to_mod : gen -> gen -> gen
val fpc_fpv_mul : gen -> gen -> gen -> gen
val fpc_fp_mul : gen -> gen -> gen -> gen
val fpc_center : gen -> gen -> gen -> gen
val fpc_center_inplace : gen -> gen -> gen -> unit
val fpc_red : gen -> gen -> gen
val fpc_to_mod : gen -> gen -> gen
val fpm_add : gen -> gen -> gen -> gen
val fpm_fp_mul : gen -> gen -> gen -> gen
val fpm_fpc_mul : gen -> gen -> gen -> gen
val fpm_fpc_mul_fpx : gen -> gen -> gen -> Signed.long -> gen
val fpm_center : gen -> gen -> gen -> gen
val fpm_center_inplace : gen -> gen -> gen -> unit
val fpm_mul : gen -> gen -> gen -> gen
val fpm_powu : gen -> pari_ulong -> gen -> gen
val fpm_red : gen -> gen -> gen
val fpm_sub : gen -> gen -> gen -> gen
val fpm_to_mod : gen -> gen -> gen
val fpms_fpc_mul : gen -> gen -> gen -> gen
val fpms_fpcs_solve : gen -> gen -> Signed.long -> gen -> gen
val fpms_fpcs_solve_safe : gen -> gen -> Signed.long -> gen -> gen
val fpms_leftkernel_elt : gen -> Signed.long -> gen -> gen
val fpc_add : gen -> gen -> gen -> gen
val fpc_sub : gen -> gen -> gen -> gen
val fpv_fpms_mul : gen -> gen -> gen -> gen
val fpv_add : gen -> gen -> gen -> gen
val fpv_sub : gen -> gen -> gen -> gen
val fpv_dotproduct : gen -> gen -> gen -> gen
val fpv_dotsquare : gen -> gen -> gen
val fpv_red : gen -> gen -> gen
val fpv_to_mod : gen -> gen -> gen
val fpvv_to_mod : gen -> gen -> gen
val fpx_to_mod : gen -> gen -> gen
val fpxc_to_mod : gen -> gen -> gen
val fpxm_to_mod : gen -> gen -> gen
val zabm_ker : gen -> gen -> Signed.long -> gen
val zabm_indexrank : gen -> gen -> Signed.long -> gen
val zabm_inv : gen -> gen -> Signed.long -> gen Ctypes_static.ptr -> gen
val zabm_inv_ratlift : gen -> gen -> Signed.long -> gen Ctypes_static.ptr -> gen

val zabm_pseudoinv :
  gen ->
  gen ->
  Signed.long ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val zv_zms_mul : gen -> gen -> gen
val zpms_zpcs_solve : gen -> gen -> Signed.long -> gen -> Signed.long -> gen

val gen_fpm_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  gen

val gen_zpm_dixon_wiedemann :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen

val gen_matid :
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  gen

val matid_flm : Signed.long -> gen
val matid_f2xqm : Signed.long -> gen -> gen
val matid_flxqm : Signed.long -> gen -> pari_ulong -> gen
val random_flv : Signed.long -> pari_ulong -> gen
val random_fpc : Signed.long -> gen -> gen
val random_fpv : Signed.long -> gen -> gen
val scalar_flm : Signed.long -> Signed.long -> gen
val zcs_to_zc : gen -> Signed.long -> gen
val zms_to_zm : gen -> Signed.long -> gen
val zms_zc_mul : gen -> gen -> gen
val zmv_to_flmv : gen -> pari_ulong -> gen
val flx_teichmuller : gen -> pari_ulong -> Signed.long -> gen
val z2_sqrt : gen -> Signed.long -> gen
val zp_div : gen -> gen -> gen -> Signed.long -> gen
val zp_exp : gen -> gen -> pari_ulong -> gen
val zp_inv : gen -> gen -> Signed.long -> gen
val zp_invlift : gen -> gen -> gen -> Signed.long -> gen
val zp_log : gen -> gen -> pari_ulong -> gen
val zp_sqrt : gen -> gen -> Signed.long -> gen
val zp_sqrtlift : gen -> gen -> gen -> Signed.long -> gen
val zp_sqrtnlift : gen -> gen -> gen -> gen -> Signed.long -> gen
val zpm_invlift : gen -> gen -> gen -> Signed.long -> gen
val zpx_frobenius : gen -> gen -> Signed.long -> gen
val zpx_zpxq_liftroot : gen -> gen -> gen -> gen -> Signed.long -> gen

val zpx_zpxq_liftroot_ea :
  gen ->
  gen ->
  gen ->
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val zpx_liftfact : gen -> gen -> gen -> gen -> Signed.long -> gen
val zpx_liftroot : gen -> gen -> gen -> Signed.long -> gen
val zpx_liftroots : gen -> gen -> gen -> Signed.long -> gen
val zpx_roots : gen -> gen -> Signed.long -> gen
val zpxq_div : gen -> gen -> gen -> gen -> gen -> Signed.long -> gen
val zpxq_inv : gen -> gen -> gen -> Signed.long -> gen
val zpxq_invlift : gen -> gen -> gen -> gen -> Signed.long -> gen
val zpxq_log : gen -> gen -> gen -> Signed.long -> gen
val zpxq_sqrt : gen -> gen -> gen -> Signed.long -> gen
val zpxq_sqrtnlift : gen -> gen -> gen -> gen -> gen -> Signed.long -> gen
val zpxqm_prodfrobenius : gen -> gen -> gen -> Signed.long -> gen
val zpxqx_digits : gen -> gen -> gen -> gen -> gen -> Signed.long -> gen

val zpxqx_divrem :
  gen -> gen -> gen -> gen -> gen -> Signed.long -> gen Ctypes_static.ptr -> gen

val zpxqx_liftfact : gen -> gen -> gen -> gen -> gen -> Signed.long -> gen
val zpxqx_liftroot : gen -> gen -> gen -> gen -> Signed.long -> gen

val zpxqx_liftroot_vald :
  gen -> gen -> Signed.long -> gen -> gen -> Signed.long -> gen

val zpxqx_liftroots : gen -> gen -> gen -> gen -> Signed.long -> gen
val zpxqx_roots : gen -> gen -> gen -> Signed.long -> gen

val zpxqx_zpxqxq_liftroot :
  gen -> gen -> gen -> gen -> gen -> Signed.long -> gen

val zq_sqrtnlift : gen -> gen -> gen -> gen -> gen -> Signed.long -> gen
val zqx_zqxq_liftroot : gen -> gen -> gen -> gen -> gen -> Signed.long -> gen
val zqx_liftfact : gen -> gen -> gen -> gen -> gen -> Signed.long -> gen
val zqx_liftroot : gen -> gen -> gen -> gen -> Signed.long -> gen
val zqx_roots : gen -> gen -> gen -> Signed.long -> gen

val gen_zpm_dixon :
  gen ->
  gen ->
  gen ->
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen -> gen)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_zpm_newton :
  gen ->
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen

val gen_zpx_dixon :
  gen ->
  gen ->
  gen ->
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen -> gen)
  Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_zpx_newton :
  gen ->
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen

val polteichmuller : gen -> pari_ulong -> Signed.long -> gen
val polhensellift : gen -> gen -> gen -> Signed.long -> gen
val quadratic_prec_mask : Signed.long -> pari_ulong
val qx_factor : gen -> gen
val zx_factor : gen -> gen
val zx_is_irred : gen -> Signed.long
val zx_squff : gen -> gen Ctypes_static.ptr -> gen
val polcyclofactors : gen -> gen
val poliscyclo : gen -> Signed.long
val poliscycloprod : gen -> Signed.long
val rg_rgc_sub : gen -> gen -> gen
val rgc_rg_add : gen -> gen -> gen
val rgc_rg_div : gen -> gen -> gen
val rgc_rg_mul : gen -> gen -> gen
val rgc_rg_sub : gen -> gen -> gen
val rgc_rgm_mul : gen -> gen -> gen
val rgc_rgv_mul : gen -> gen -> gen
val rgc_add : gen -> gen -> gen
val rgc_is_ei : gen -> Signed.long
val rgc_neg : gen -> gen
val rgc_sub : gen -> gen -> gen
val rgm_rg_add : gen -> gen -> gen
val rgm_rg_add_shallow : gen -> gen -> gen
val rgm_rg_div : gen -> gen -> gen
val rgm_rg_mul : gen -> gen -> gen
val rgm_rg_sub : gen -> gen -> gen
val rgm_rg_sub_shallow : gen -> gen -> gen
val rgm_rgc_mul : gen -> gen -> gen
val rgm_rgv_mul : gen -> gen -> gen
val rgm_add : gen -> gen -> gen
val rgm_det_triangular : gen -> gen
val rgm_is_qm : gen -> int
val rgm_is_zm : gen -> int
val rgm_isdiagonal : gen -> int
val rgm_isidentity : gen -> int
val rgm_isscalar : gen -> gen -> int
val rgm_mul : gen -> gen -> gen
val rgm_multosym : gen -> gen -> gen
val rgm_neg : gen -> gen
val rgm_powers : gen -> Signed.long -> gen
val rgm_sqr : gen -> gen
val rgm_sub : gen -> gen -> gen
val rgm_sumcol : gen -> gen
val rgm_transmul : gen -> gen -> gen
val rgm_transmultosym : gen -> gen -> gen
val rgmrow_zc_mul : gen -> gen -> Signed.long -> gen
val rgm_zc_mul : gen -> gen -> gen
val rgm_zm_mul : gen -> gen -> gen
val rgmrow_rgc_mul : gen -> gen -> Signed.long -> gen
val rgv_rgm_mul : gen -> gen -> gen
val rgv_rgc_mul : gen -> gen -> gen
val rgv_rg_mul : gen -> gen -> gen
val rgv_add : gen -> gen -> gen
val rgv_dotproduct : gen -> gen -> gen
val rgv_dotsquare : gen -> gen
val rgv_is_zmv : gen -> int
val rgv_kill0 : gen -> gen
val rgv_neg : gen -> gen
val rgv_prod : gen -> gen
val rgv_sub : gen -> gen -> gen
val rgv_sum : gen -> gen
val rgv_sumpart : gen -> Signed.long -> gen
val rgv_sumpart2 : gen -> Signed.long -> Signed.long -> gen
val rgv_zc_mul : gen -> gen -> gen
val rgv_zm_mul : gen -> gen -> gen
val rgx_rgm_eval : gen -> gen -> gen
val rgx_rgmv_eval : gen -> gen -> gen
val isdiagonal : gen -> int
val matid : Signed.long -> gen
val scalarcol : gen -> Signed.long -> gen
val scalarcol_shallow : gen -> Signed.long -> gen
val scalarmat : gen -> Signed.long -> gen
val scalarmat_shallow : gen -> Signed.long -> gen
val scalarmat_s : Signed.long -> Signed.long -> gen
val kronecker_to_mod : gen -> gen -> gen
val qx_zxqv_eval : gen -> gen -> gen -> gen
val qxq_charpoly : gen -> gen -> Signed.long -> gen
val qxq_powers : gen -> Signed.long -> gen -> gen
val qxq_to_mod_shallow : gen -> gen -> gen
val qxqc_to_mod_shallow : gen -> gen -> gen
val qxqm_to_mod_shallow : gen -> gen -> gen
val qxqv_to_mod : gen -> gen -> gen
val qxqx_homogenous_evalpow : gen -> gen -> gen -> gen -> gen
val qxqx_to_mod_shallow : gen -> gen -> gen
val qxqxv_to_mod : gen -> gen -> gen
val qxv_qxq_eval : gen -> gen -> gen -> gen
val qxy_qxq_evalx : gen -> gen -> gen -> gen
val rg_rgx_sub : gen -> gen -> gen
val rg_get_0 : gen -> gen
val rg_get_1 : gen -> gen
val rg_to_rgc : gen -> Signed.long -> gen
val rgm_to_rgxv : gen -> Signed.long -> gen
val rgm_to_rgxv_reverse : gen -> Signed.long -> gen
val rgm_to_rgxx : gen -> Signed.long -> Signed.long -> gen
val rgv_to_rgx : gen -> Signed.long -> gen
val rgv_to_rgm : gen -> Signed.long -> gen
val rgv_to_rgx_reverse : gen -> Signed.long -> gen
val rgx_rgxq_eval : gen -> gen -> gen -> gen
val rgx_rgxqv_eval : gen -> gen -> gen -> gen
val rgx_rgxn_eval : gen -> gen -> Signed.long -> gen
val rgx_rgxnv_eval : gen -> gen -> Signed.long -> gen
val rgx_rg_add : gen -> gen -> gen
val rgx_rg_add_shallow : gen -> gen -> gen
val rgx_rg_div : gen -> gen -> gen
val rgx_rg_divexact : gen -> gen -> gen
val rgx_rg_eval_bk : gen -> gen -> gen
val rgx_rg_mul : gen -> gen -> gen
val rgx_rg_sub : gen -> gen -> gen
val rgx_rgv_eval : gen -> gen -> gen
val rgx_add : gen -> gen -> gen
val rgx_addmulxn_shallow : gen -> gen -> Signed.long -> gen
val rgx_addmulxn : gen -> gen -> Signed.long -> gen
val rgx_addspec : gen -> gen -> Signed.long -> Signed.long -> gen
val rgx_addspec_shallow : gen -> gen -> Signed.long -> Signed.long -> gen
val rgx_affine : gen -> gen -> gen -> gen
val rgx_blocks : gen -> Signed.long -> Signed.long -> gen
val rgx_deflate : gen -> Signed.long -> gen
val rgx_deriv : gen -> gen
val rgx_digits : gen -> gen -> gen
val rgx_div_by_x_x : gen -> gen -> gen Ctypes_static.ptr -> gen
val rgx_divrem : gen -> gen -> gen Ctypes_static.ptr -> gen
val rgx_divs : gen -> Signed.long -> gen
val rgx_equal : gen -> gen -> Signed.long
val rgx_even_odd : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit
val rgx_homogenize : gen -> Signed.long -> gen
val rgx_homogenous_evalpow : gen -> gen -> gen -> gen
val rgx_inflate : gen -> Signed.long -> gen
val rgx_mul : gen -> gen -> gen
val rgx_mul_i : gen -> gen -> gen
val rgx_mul_normalized : gen -> Signed.long -> gen -> Signed.long -> gen
val rgx_mul2n : gen -> Signed.long -> gen
val rgx_mulxn : gen -> Signed.long -> gen
val rgx_mulhigh_i : gen -> gen -> Signed.long -> gen
val rgx_muls : gen -> Signed.long -> gen
val rgx_mulspec : gen -> gen -> Signed.long -> Signed.long -> gen
val rgx_neg : gen -> gen
val rgx_normalize : gen -> gen
val rgx_pseudodivrem : gen -> gen -> gen Ctypes_static.ptr -> gen
val rgx_pseudorem : gen -> gen -> gen
val rgx_recip : gen -> gen
val rgx_recip_i : gen -> gen
val rgx_recip_shallow : gen -> gen
val rgx_rem : gen -> gen -> gen
val rgx_renormalize_lg : gen -> Signed.long -> gen
val rgx_rescale : gen -> gen -> gen
val rgx_rotate_shallow : gen -> Signed.long -> Signed.long -> gen
val rgx_shift : gen -> Signed.long -> gen
val rgx_shift_shallow : gen -> Signed.long -> gen
val rgx_splitting : gen -> Signed.long -> gen
val rgx_sqr : gen -> gen
val rgx_sqr_i : gen -> gen
val rgx_sqrhigh_i : gen -> Signed.long -> gen
val rgx_sqrspec : gen -> Signed.long -> gen
val rgx_sub : gen -> gen -> gen
val rgx_to_rgc : gen -> Signed.long -> gen
val rgx_translate : gen -> gen -> gen
val rgx_unscale : gen -> gen -> gen
val rgxq_matrix_pow : gen -> Signed.long -> Signed.long -> gen -> gen
val rgxq_norm : gen -> gen -> gen
val rgxq_pow : gen -> gen -> gen -> gen
val rgxq_powers : gen -> Signed.long -> gen -> gen
val rgxq_powu : gen -> pari_ulong -> gen -> gen
val rgxq_trace : gen -> gen -> gen
val rgxqc_red : gen -> gen -> gen
val rgxqm_mul : gen -> gen -> gen -> gen
val rgxqm_red : gen -> gen -> gen
val rgxqv_rgxq_mul : gen -> gen -> gen -> gen
val rgxqv_factorback : gen -> gen -> gen -> gen
val rgxqv_red : gen -> gen -> gen
val rgxqx_rgxq_mul : gen -> gen -> gen -> gen
val rgxqx_divrem : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val rgxqx_mul : gen -> gen -> gen -> gen
val rgxqx_powers : gen -> Signed.long -> gen -> gen
val rgxqx_pseudodivrem : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val rgxqx_pseudorem : gen -> gen -> gen -> gen
val rgxqx_red : gen -> gen -> gen
val rgxqx_sqr : gen -> gen -> gen
val rgxqx_translate : gen -> gen -> gen -> gen
val rgxv_rgv_eval : gen -> gen -> gen
val rgxv_prod : gen -> gen
val rgxv_rescale : gen -> gen -> gen
val rgxv_to_rgm : gen -> Signed.long -> gen
val rgxv_unscale : gen -> gen -> gen
val rgxx_to_rgm : gen -> Signed.long -> gen
val rgxy_degreex : gen -> Signed.long
val rgxy_derivx : gen -> gen
val rgxy_swap : gen -> Signed.long -> Signed.long -> gen
val rgxy_swapspec : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val rgxn_div : gen -> gen -> Signed.long -> gen
val rgxn_div_i : gen -> gen -> Signed.long -> gen
val rgxn_eval : gen -> gen -> Signed.long -> gen
val rgxn_exp : gen -> Signed.long -> gen
val rgxn_expint : gen -> Signed.long -> gen
val rgxn_inv : gen -> Signed.long -> gen
val rgxn_inv_i : gen -> Signed.long -> gen
val rgxn_mul : gen -> gen -> Signed.long -> gen
val rgxn_powers : gen -> Signed.long -> Signed.long -> gen
val rgxn_recip_shallow : gen -> Signed.long -> gen
val rgxn_red_shallow : gen -> Signed.long -> gen
val rgxn_reverse : gen -> Signed.long -> gen
val rgxn_sqr : gen -> Signed.long -> gen
val rgxn_sqrt : gen -> Signed.long -> gen
val rgxnv_red_shallow : gen -> Signed.long -> gen
val rgxn_powu : gen -> pari_ulong -> Signed.long -> gen
val rgxn_powu_i : gen -> pari_ulong -> Signed.long -> gen
val zx_translate : gen -> gen -> gen
val zx_unscale2n : gen -> Signed.long -> gen
val zx_unscale : gen -> gen -> gen
val zx_unscale_div : gen -> gen -> gen
val zx_unscale_divpow : gen -> gen -> Signed.long -> gen
val zx_z_unscale : gen -> Signed.long -> gen
val zxq_powers : gen -> Signed.long -> gen -> gen
val zxq_powu : gen -> pari_ulong -> gen -> gen
val zxqx_dvd : gen -> gen -> gen -> int
val brent_kung_optpow : Signed.long -> Signed.long -> Signed.long -> Signed.long

val gen_bkeval :
  gen ->
  Signed.long ->
  gen ->
  int ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen -> gen)
  Ctypes_static.static_funptr ->
  gen

val gen_bkeval_powers :
  gen ->
  Signed.long ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_algebra Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen -> gen)
  Ctypes_static.static_funptr ->
  gen

val get_rg_algebra : unit -> bb_algebra Ctypes.structure Ctypes_static.ptr
val rfrac_deflate_order : gen -> Signed.long
val rfrac_deflate_max : gen -> Signed.long Ctypes_static.ptr -> gen
val rfrac_deflate : gen -> Signed.long -> gen
val zgc_g_mul_inplace : gen -> gen -> unit
val zgcs_add : gen -> gen -> gen
val g_zgc_mul : gen -> gen -> gen
val g_zg_mul : gen -> gen -> gen
val zgc_g_mul : gen -> gen -> gen
val zgc_z_mul : gen -> gen -> gen
val zg_g_mul : gen -> gen -> gen
val zg_z_mul : gen -> gen -> gen
val zg_add : gen -> gen -> gen
val zg_mul : gen -> gen -> gen
val zg_neg : gen -> gen
val zg_normalize : gen -> gen
val zg_sub : gen -> gen -> gen
val flc_lincomb1_inplace : gen -> gen -> pari_ulong -> pari_ulong -> unit
val vecsmall_prod : gen -> gen
val qm_qc_mul : gen -> gen -> gen
val qm_det : gen -> gen
val qm_ker : gen -> gen
val qm_mul : gen -> gen -> gen
val qm_sqr : gen -> gen
val rgm_check_zm : gen -> string -> unit
val rgv_check_zv : gen -> string -> unit
val z_zc_sub : gen -> gen -> gen
val zv_zc_mul : gen -> gen -> gen
val zc_q_mul : gen -> gen -> gen
val zc_z_add : gen -> gen -> gen
val zc_z_div : gen -> gen -> gen
val zc_z_divexact : gen -> gen -> gen
val zc_z_sub : gen -> gen -> gen
val zc_zv_mul : gen -> gen -> gen
val zc_divexactu : gen -> pari_ulong -> gen
val zc_add : gen -> gen -> gen
val zc_copy : gen -> gen
val zc_hnfremdiv : gen -> gen -> gen Ctypes_static.ptr -> gen
val zc_is_ei : gen -> Signed.long
val zc_lincomb : gen -> gen -> gen -> gen -> gen
val zc_lincomb1_inplace : gen -> gen -> gen -> unit
val zc_lincomb1_inplace_i : gen -> gen -> gen -> Signed.long -> unit
val zc_neg : gen -> gen
val zc_reducemodlll : gen -> gen -> gen
val zc_reducemodmatrix : gen -> gen -> gen
val zc_sub : gen -> gen -> gen
val zc_z_mul : gen -> Signed.long -> gen
val zm_q_mul : gen -> gen -> gen
val zm_z_div : gen -> gen -> gen
val zm_z_divexact : gen -> gen -> gen
val zm_z_mul : gen -> gen -> gen
val zm_add : gen -> gen -> gen
val zm_det_triangular : gen -> gen
val zm_diag_mul : gen -> gen -> gen
val zm_divexactu : gen -> pari_ulong -> gen
val zm_equal : gen -> gen -> int
val zm_equal0 : gen -> int
val zm_hnfdivrem : gen -> gen -> gen Ctypes_static.ptr -> gen
val zm_ishnf : gen -> int
val zm_isdiagonal : gen -> int
val zm_isidentity : gen -> int
val zm_isscalar : gen -> gen -> int
val zm_max_lg : gen -> Signed.long
val zm_mul_diag : gen -> gen -> gen
val zm_multosym : gen -> gen -> gen
val zm_neg : gen -> gen
val zm_nm_mul : gen -> gen -> gen
val zm_pow : gen -> gen -> gen
val zm_powu : gen -> pari_ulong -> gen
val zm_reducemodlll : gen -> gen -> gen
val zm_reducemodmatrix : gen -> gen -> gen
val zm_sqr : gen -> gen
val zm_sub : gen -> gen -> gen
val zm_supnorm : gen -> gen
val zm_transmul : gen -> gen -> gen
val zm_transmultosym : gen -> gen -> gen
val zm_togglesign : gen -> unit
val zm_zm_mul : gen -> gen -> gen
val zmrow_zc_mul : gen -> gen -> Signed.long -> gen
val zmrow_equal0 : gen -> Signed.long -> int
val zv_abscmp : gen -> gen -> int
val zv_cmp : gen -> gen -> int
val zv_dotsquare : gen -> gen
val zv_max_lg : gen -> Signed.long
val zv_to_nv : gen -> gen
val zv_togglesign : gen -> unit
val gram_matrix : gen -> gen
val nm_z_mul : gen -> gen -> gen
val zm_mul : gen -> gen -> gen
val zm_to_flm : gen -> pari_ulong -> gen
val zm_to_zm : gen -> gen
val zm_zc_mul : gen -> gen -> gen
val zmv_to_zmv : gen -> gen
val zv_abs : gen -> gen
val zv_content : gen -> Signed.long
val zv_dotproduct : gen -> gen -> Signed.long
val zv_equal : gen -> gen -> int
val zv_equal0 : gen -> int
val zv_neg : gen -> gen
val zv_neg_inplace : gen -> gen
val zv_prod : gen -> Signed.long
val zv_prod_z : gen -> gen
val zv_sum : gen -> Signed.long
val zv_sumpart : gen -> Signed.long -> Signed.long
val zv_to_flv : gen -> pari_ulong -> gen
val zv_z_mul : gen -> Signed.long -> gen
val zv_zm_mul : gen -> gen -> gen
val zvv_equal : gen -> gen -> int
val kronecker_to_zxqx : gen -> gen -> gen
val kronecker_to_zxx : gen -> Signed.long -> Signed.long -> gen
val qx_zx_rem : gen -> gen -> gen
val qx_mul : gen -> gen -> gen
val qx_sqr : gen -> gen
val qxqm_mul : gen -> gen -> gen -> gen
val qxqm_sqr : gen -> gen -> gen
val qxqx_qxq_mul : gen -> gen -> gen -> gen
val qxqx_mul : gen -> gen -> gen -> gen
val qxqx_powers : gen -> Signed.long -> gen -> gen
val qxqx_sqr : gen -> gen -> gen
val rgx_check_qx : gen -> string -> unit
val rgx_check_zx : gen -> string -> unit
val rgx_check_zxx : gen -> string -> unit
val z_zx_sub : gen -> gen -> gen
val zx_z_add : gen -> gen -> gen
val zx_z_add_shallow : gen -> gen -> gen
val zx_z_eval : gen -> gen -> gen
val zx_z_mul : gen -> gen -> gen
val zx_z_sub : gen -> gen -> gen
val zx_add : gen -> gen -> gen
val zx_affine : gen -> gen -> gen -> gen
val zx_copy : gen -> gen
val zx_deriv : gen -> gen
val zx_digits : gen -> gen -> gen
val zxv_zx_fromdigits : gen -> gen -> gen
val zx_div_by_x_1 : gen -> gen Ctypes_static.ptr -> gen
val zx_divuexact : gen -> pari_ulong -> gen
val zx_equal : gen -> gen -> int
val zx_eval1 : gen -> gen
val zx_max_lg : gen -> Signed.long
val zx_mod_xnm1 : gen -> pari_ulong -> gen
val zx_mul : gen -> gen -> gen
val zx_mulspec : gen -> gen -> Signed.long -> Signed.long -> gen
val zx_mulu : gen -> pari_ulong -> gen
val zx_neg : gen -> gen
val zx_rem : gen -> gen -> gen
val zx_remi2n : gen -> Signed.long -> gen
val zx_rescale2n : gen -> Signed.long -> gen
val zx_rescale : gen -> gen -> gen
val zx_rescale_lt : gen -> gen
val zx_shifti : gen -> Signed.long -> gen
val zx_sqr : gen -> gen
val zx_sqrspec : gen -> Signed.long -> gen
val zx_sub : gen -> gen -> gen
val zx_val : gen -> Signed.long
val zx_valrem : gen -> gen Ctypes_static.ptr -> Signed.long
val zxc_to_flxc : gen -> pari_ulong -> Signed.long -> gen
val zxm_to_flxm : gen -> pari_ulong -> Signed.long -> gen
val zxqm_mul : gen -> gen -> gen -> gen
val zxqm_sqr : gen -> gen -> gen
val zxqx_zxq_mul : gen -> gen -> gen -> gen
val zxqx_sqr : gen -> gen -> gen
val zxqx_mul : gen -> gen -> gen -> gen
val zxt_remi2n : gen -> Signed.long -> gen
val zxv_z_mul : gen -> gen -> gen
val zxv_dotproduct : gen -> gen -> gen
val zxv_equal : gen -> gen -> int
val zxv_remi2n : gen -> Signed.long -> gen
val zxx_z_divexact : gen -> gen -> gen
val zxx_z_mul : gen -> gen -> gen
val zxx_z_add_shallow : gen -> gen -> gen
val zxx_evalx0 : gen -> gen
val zxx_max_lg : gen -> Signed.long
val zxx_mul_kronecker : gen -> gen -> Signed.long -> gen
val zxx_renormalize : gen -> Signed.long -> gen
val zxx_sqr_kronecker : gen -> Signed.long -> gen
val rgxx_to_kronecker : gen -> Signed.long -> gen
val rgxx_to_kronecker_spec : gen -> Signed.long -> Signed.long -> gen
val zxn_mul : gen -> gen -> Signed.long -> gen
val zxn_sqr : gen -> Signed.long -> gen
val scalar_zx : gen -> Signed.long -> gen
val scalar_zx_shallow : gen -> Signed.long -> gen
val zx_to_zx : gen -> gen
val zx_z_divexact : gen -> Signed.long -> gen
val alg_centralproj : gen -> gen -> Signed.long -> gen
val alg_complete : gen -> gen -> gen -> gen -> Signed.long -> gen
val alg_csa_table : gen -> gen -> Signed.long -> Signed.long -> gen
val alg_cyclic : gen -> gen -> gen -> Signed.long -> gen
val alg_get_absdim : gen -> Signed.long
val alg_get_abssplitting : gen -> gen
val alg_get_aut : gen -> gen
val algaut : gen -> gen
val alg_get_auts : gen -> gen
val alg_get_b : gen -> gen
val algb : gen -> gen
val algcenter : gen -> gen
val alg_get_center : gen -> gen
val alg_get_char : gen -> gen
val algchar : gen -> gen
val alg_get_degree : gen -> Signed.long
val algdegree : gen -> Signed.long
val alg_get_dim : gen -> Signed.long
val algdim : gen -> Signed.long -> Signed.long
val alg_get_hasse_f : gen -> gen
val alghassef : gen -> gen
val alg_get_hasse_i : gen -> gen
val alghassei : gen -> gen
val alg_get_invbasis : gen -> gen
val alginvbasis : gen -> gen
val alg_get_multable : gen -> gen
val alg_get_basis : gen -> gen
val algbasis : gen -> gen
val alg_get_relmultable : gen -> gen
val algrelmultable : gen -> gen
val alg_get_splitpol : gen -> gen
val alg_get_splittingfield : gen -> gen
val algsplittingfield : gen -> gen
val alg_get_splittingbasis : gen -> gen
val alg_get_splittingbasisinv : gen -> gen
val alg_get_splittingdata : gen -> gen
val algsplittingdata : gen -> gen
val alg_get_tracebasis : gen -> gen

val alg_hasse :
  gen -> Signed.long -> gen -> gen -> Signed.long -> Signed.long -> gen

val alg_hilbert : gen -> gen -> gen -> Signed.long -> Signed.long -> gen
val alg_matrix : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val alg_model : gen -> gen -> Signed.long
val alg_quotient : gen -> gen -> Signed.long -> gen
val algradical : gen -> gen
val algsimpledec : gen -> Signed.long -> gen
val algsimpledec_ss : gen -> Signed.long -> gen
val algsubalg : gen -> gen -> gen
val alg_type : gen -> Signed.long
val algadd : gen -> gen -> gen -> gen
val algalgtobasis : gen -> gen -> gen
val algbasistoalg : gen -> gen -> gen
val algcharpoly : gen -> gen -> Signed.long -> Signed.long -> gen
val algdisc : gen -> gen
val algdivl : gen -> gen -> gen -> gen
val algdivr : gen -> gen -> gen -> gen
val alggroup : gen -> gen -> gen
val alggroupcenter : gen -> gen -> gen Ctypes_static.ptr -> gen
val alghasse : gen -> gen -> gen
val alginit : gen -> gen -> Signed.long -> Signed.long -> gen
val algindex : gen -> gen -> Signed.long
val alginv : gen -> gen -> gen
val algisassociative : gen -> gen -> int
val algiscommutative : gen -> int
val algisdivision : gen -> gen -> int
val algisramified : gen -> gen -> int
val algissemisimple : gen -> int
val algissimple : gen -> Signed.long -> int
val algissplit : gen -> gen -> int
val algisdivl : gen -> gen -> gen -> gen Ctypes_static.ptr -> int
val algisinv : gen -> gen -> gen Ctypes_static.ptr -> int
val algmakeintegral : gen -> Signed.long -> gen
val algmul : gen -> gen -> gen -> gen
val algmultable : gen -> gen
val alglat_get_primbasis : gen -> gen
val alglat_get_scalar : gen -> gen
val alglatadd : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val alglatcontains : gen -> gen -> gen -> gen Ctypes_static.ptr -> int
val alglatelement : gen -> gen -> gen -> gen
val alglathnf : gen -> gen -> gen -> gen
val alglatindex : gen -> gen -> gen -> gen
val alglatinter : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val alglatmul : gen -> gen -> gen -> gen
val alglatlefttransporter : gen -> gen -> gen -> gen
val alglatrighttransporter : gen -> gen -> gen -> gen
val alglatsubset : gen -> gen -> gen -> gen Ctypes_static.ptr -> int
val algneg : gen -> gen -> gen
val algnorm : gen -> gen -> Signed.long -> gen
val algpoleval : gen -> gen -> gen -> gen
val algpow : gen -> gen -> gen -> gen
val algprimesubalg : gen -> gen
val algramifiedplaces : gen -> gen
val algrandom : gen -> gen -> gen
val algsplit : gen -> Signed.long -> gen
val algtomatrix : gen -> gen -> Signed.long -> gen
val algsqr : gen -> gen -> gen
val algsub : gen -> gen -> gen -> gen
val algtableinit : gen -> gen -> gen
val algtensor : gen -> gen -> Signed.long -> gen
val algtrace : gen -> gen -> Signed.long -> gen
val algtype : gen -> Signed.long
val bnfgwgeneric : gen -> gen -> gen -> gen -> Signed.long -> gen
val checkalg : gen -> unit
val checkhasse : gen -> gen -> gen -> Signed.long -> unit
val checklat : gen -> gen -> unit
val conjclasses_algcenter : gen -> gen -> gen
val galoischardet : gen -> gen -> Signed.long -> gen
val galoischarpoly : gen -> gen -> Signed.long -> gen
val galoischartable : gen -> gen
val nfgrunwaldwang : gen -> gen -> gen -> gen -> Signed.long -> gen
val nfgwkummer : gen -> gen -> gen -> gen -> Signed.long -> gen
val f2ms_colelim : gen -> Signed.long -> gen
val f2m_image : gen -> gen
val f2m_indexrank : gen -> gen
val f2m_suppl : gen -> gen
val f2xqm_f2xqc_gauss : gen -> gen -> gen -> gen
val f2xqm_f2xqc_invimage : gen -> gen -> gen -> gen
val f2xqm_f2xqc_mul : gen -> gen -> gen -> gen
val f2xqm_deplin : gen -> gen -> gen
val f2xqm_det : gen -> gen -> gen
val f2xqm_gauss : gen -> gen -> gen -> gen
val f2xqm_ker : gen -> gen -> gen
val f2xqm_image : gen -> gen -> gen
val f2xqm_indexrank : gen -> gen -> gen
val f2xqm_inv : gen -> gen -> gen
val f2xqm_invimage : gen -> gen -> gen -> gen
val f2xqm_mul : gen -> gen -> gen -> gen
val f2xqm_rank : gen -> gen -> Signed.long
val f2xqm_suppl : gen -> gen -> gen
val flm_image : gen -> pari_ulong -> gen
val flm_indexrank : gen -> pari_ulong -> gen
val flm_suppl : gen -> pari_ulong -> gen
val flxqm_flxqc_gauss : gen -> gen -> gen -> pari_ulong -> gen
val flxqm_flxqc_invimage : gen -> gen -> gen -> pari_ulong -> gen
val flxqm_flxqc_mul : gen -> gen -> gen -> pari_ulong -> gen
val flxqm_deplin : gen -> gen -> pari_ulong -> gen
val flxqm_det : gen -> gen -> pari_ulong -> gen
val flxqm_gauss : gen -> gen -> gen -> pari_ulong -> gen
val flxqm_ker : gen -> gen -> pari_ulong -> gen
val flxqm_image : gen -> gen -> pari_ulong -> gen
val flxqm_indexrank : gen -> gen -> pari_ulong -> gen
val flxqm_inv : gen -> gen -> pari_ulong -> gen
val flxqm_invimage : gen -> gen -> gen -> pari_ulong -> gen
val flxqm_mul : gen -> gen -> gen -> pari_ulong -> gen
val flxqm_rank : gen -> gen -> pari_ulong -> Signed.long
val flxqm_suppl : gen -> gen -> pari_ulong -> gen
val fpm_fpc_gauss : gen -> gen -> gen -> gen
val fpm_fpc_invimage : gen -> gen -> gen -> gen
val fpm_deplin : gen -> gen -> gen
val fpm_det : gen -> gen -> gen
val fpm_gauss : gen -> gen -> gen -> gen
val fpm_image : gen -> gen -> gen
val fpm_indexrank : gen -> gen -> gen
val fpm_intersect : gen -> gen -> gen -> gen
val fpm_intersect_i : gen -> gen -> gen -> gen
val fpm_inv : gen -> gen -> gen
val fpm_invimage : gen -> gen -> gen -> gen
val fpm_ker : gen -> gen -> gen
val fpm_rank : gen -> gen -> Signed.long
val fpm_suppl : gen -> gen -> gen
val fqm_fqc_gauss : gen -> gen -> gen -> gen -> gen
val fqm_fqc_invimage : gen -> gen -> gen -> gen -> gen
val fqm_fqc_mul : gen -> gen -> gen -> gen -> gen
val fqm_deplin : gen -> gen -> gen -> gen
val fqm_det : gen -> gen -> gen -> gen
val fqm_gauss : gen -> gen -> gen -> gen -> gen
val fqm_ker : gen -> gen -> gen -> gen
val fqm_image : gen -> gen -> gen -> gen
val fqm_indexrank : gen -> gen -> gen -> gen
val fqm_inv : gen -> gen -> gen -> gen
val fqm_invimage : gen -> gen -> gen -> gen -> gen
val fqm_mul : gen -> gen -> gen -> gen -> gen
val fqm_rank : gen -> gen -> gen -> Signed.long
val fqm_suppl : gen -> gen -> gen -> gen
val qm_image_shallow : gen -> gen
val qm_image : gen -> gen
val qm_gauss : gen -> gen -> gen
val qm_gauss_i : gen -> gen -> Signed.long -> gen
val qm_indexrank : gen -> gen
val qm_inv : gen -> gen
val qm_rank : gen -> Signed.long
val rgm_fp_init : gen -> gen -> pari_ulong Ctypes_static.ptr -> gen
val rgm_hadamard : gen -> gen
val rgm_rgc_invimage : gen -> gen -> gen
val rgm_diagonal : gen -> gen
val rgm_diagonal_shallow : gen -> gen
val rgm_inv : gen -> gen
val rgm_inv_upper : gen -> gen
val rgm_invimage : gen -> gen -> gen
val rgm_solve : gen -> gen -> gen
val rgm_solve_realimag : gen -> gen -> gen

val rgms_structelim :
  gen ->
  Signed.long ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  unit

val zm_det : gen -> gen
val zm_detmult : gen -> gen
val zm_gauss : gen -> gen -> gen
val zm_ker : gen -> gen
val zm_imagecompl : gen -> gen
val zm_indeximage : gen -> gen
val zm_indexrank : gen -> gen
val zm_inv : gen -> gen Ctypes_static.ptr -> gen
val zm_inv_ratlift : gen -> gen Ctypes_static.ptr -> gen
val zm_pseudoinv : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen
val zm_rank : gen -> Signed.long
val zlm_gauss : gen -> gen -> pari_ulong -> Signed.long -> gen -> gen
val closemodinvertible : gen -> gen -> gen
val deplin : gen -> gen
val det : gen -> gen
val det0 : gen -> Signed.long -> gen
val det2 : gen -> gen
val detint : gen -> gen
val eigen : gen -> Signed.long -> gen
val gauss : gen -> gen -> gen
val gaussmodulo : gen -> gen -> gen -> gen
val gaussmodulo2 : gen -> gen -> gen -> gen

val gen_gauss :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_gauss_pivot :
  gen ->
  Signed.long Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_det :
  gen ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_ker :
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_matcolinvimage :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_matcolmul :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_matinvimage :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_matmul :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_field Ctypes.structure Ctypes_static.ptr ->
  gen

val image : gen -> gen
val image2 : gen -> gen
val imagecompl : gen -> gen
val indexrank : gen -> gen
val inverseimage : gen -> gen -> gen
val ker : gen -> gen
val mateigen : gen -> Signed.long -> Signed.long -> gen
val matimage0 : gen -> Signed.long -> gen
val matker0 : gen -> Signed.long -> gen
val rank : gen -> Signed.long
val reducemodinvertible : gen -> gen -> gen
val reducemodlll : gen -> gen -> gen
val split_realimag : gen -> Signed.long -> Signed.long -> gen
val suppl : gen -> gen
val flm_charpoly : gen -> pari_ulong -> gen
val flm_hess : gen -> pari_ulong -> gen
val fpm_charpoly : gen -> gen -> gen
val fpm_hess : gen -> gen -> gen
val frobeniusform : gen -> Signed.long -> gen
val qm_minors_coprime : gen -> gen -> gen
val qm_imz : gen -> gen

val qm_imz_all :
  gen -> gen Ctypes_static.ptr -> Signed.long -> Signed.long -> gen

val qm_imz_hnf : gen -> gen
val qm_imz_hnfall : gen -> gen Ctypes_static.ptr -> Signed.long -> gen
val qm_imq : gen -> gen

val qm_imq_all :
  gen -> gen Ctypes_static.ptr -> Signed.long -> Signed.long -> gen

val qm_imq_hnf : gen -> gen
val qm_imq_hnfall : gen -> gen Ctypes_static.ptr -> Signed.long -> gen
val qm_charpoly_zx : gen -> gen
val qm_charpoly_zx_bound : gen -> Signed.long -> gen
val zm_charpoly : gen -> gen
val adj : gen -> gen
val adjsafe : gen -> gen
val caract : gen -> Signed.long -> gen
val caradj : gen -> Signed.long -> gen Ctypes_static.ptr -> gen
val carberkowitz : gen -> Signed.long -> gen
val carhess : gen -> Signed.long -> gen
val charpoly : gen -> Signed.long -> gen
val charpoly0 : gen -> Signed.long -> Signed.long -> gen
val gnorm : gen -> gen
val gnorml1 : gen -> Signed.long -> gen
val gnorml1_fake : gen -> gen
val gnormlp : gen -> gen -> Signed.long -> gen
val gnorml2 : gen -> gen
val gsupnorm : gen -> Signed.long -> gen

val gsupnorm_aux :
  gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> Signed.long -> unit

val gtrace : gen -> gen
val hess : gen -> gen
val intersect : gen -> gen -> gen
val jacobi : gen -> Signed.long -> gen
val matadjoint0 : gen -> Signed.long -> gen
val matcompanion : gen -> gen
val matrixqz0 : gen -> gen -> gen
val minpoly : gen -> Signed.long -> gen
val qfgaussred : gen -> gen
val qfgaussred_positive : gen -> gen
val qfsign : gen -> gen
val apply0 : gen -> gen -> gen
val diagonal : gen -> gen
val diagonal_shallow : gen -> gen
val extract0 : gen -> gen -> gen -> gen
val fold0 : gen -> gen -> gen

val genapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen

val genfold :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen

val genindexselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  gen ->
  gen

val genselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  gen ->
  gen

val gtomat : gen -> gen
val gtrans : gen -> gen
val matmuldiagonal : gen -> gen -> gen
val matmultodiagonal : gen -> gen -> gen

val matslice0 :
  gen -> Signed.long -> Signed.long -> Signed.long -> Signed.long -> gen

val parapply : gen -> gen -> gen

val parfor :
  gen ->
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parfor_init :
  parfor_t Ctypes.structure Ctypes_static.ptr -> gen -> gen -> gen -> unit

val parfor_next : parfor_t Ctypes.structure Ctypes_static.ptr -> gen
val parfor_stop : parfor_t Ctypes.structure Ctypes_static.ptr -> unit

val parforeach :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforeach_init :
  parforeach_t Ctypes.structure Ctypes_static.ptr -> gen -> gen -> unit

val parforeach_next : parforeach_t Ctypes.structure Ctypes_static.ptr -> gen
val parforeach_stop : parforeach_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprime :
  gen ->
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprime_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr -> gen -> gen -> gen -> unit

val parforprime_next : parforprime_t Ctypes.structure Ctypes_static.ptr -> gen
val parforprime_stop : parforprime_t Ctypes.structure Ctypes_static.ptr -> unit

val parforprimestep :
  gen ->
  gen ->
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforprimestep_init :
  parforprime_t Ctypes.structure Ctypes_static.ptr ->
  gen ->
  gen ->
  gen ->
  gen ->
  unit

val parforvec :
  gen ->
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> Signed.long)
  Ctypes_static.static_funptr ->
  unit

val parforvec_init :
  parforvec_t Ctypes.structure Ctypes_static.ptr ->
  gen ->
  gen ->
  Signed.long ->
  unit

val parforvec_next : parforvec_t Ctypes.structure Ctypes_static.ptr -> gen
val parforvec_stop : parforvec_t Ctypes.structure Ctypes_static.ptr -> unit
val parselect : gen -> gen -> Signed.long -> gen
val select0 : gen -> gen -> Signed.long -> gen
val shallowextract : gen -> gen -> gen
val shallowmatextract : gen -> gen -> gen -> gen
val shallowtrans : gen -> gen

val vecapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen

val veccatapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen

val veccatselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen

val vecrange : gen -> gen -> gen
val vecrangess : Signed.long -> Signed.long -> gen

val vecselapply :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen

val vecselect :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  gen ->
  gen

val vecslice0 : gen -> Signed.long -> Signed.long -> gen
val vecsum : gen -> gen
val zv_diagonal : gen -> gen
val addhelp : string -> string -> unit
val arity0 : gen -> gen
val alias0 : string -> string -> unit
val compile_str : string -> gen
val delete_var : unit -> Signed.long
val fetch_user_var : string -> Signed.long
val fetch_var : unit -> Signed.long
val fetch_var_higher : unit -> Signed.long
val fetch_var_value : Signed.long -> gen -> gen
val gp_embedded : string -> string
val gp_embedded_init : Signed.long -> Signed.long -> unit
val gp_read_str : string -> gen
val gp_read_str_bitprec : string -> Signed.long -> gen
val gp_read_str_prec : string -> Signed.long -> gen

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
val readseq : string -> gen
val safegel : gen -> Signed.long -> gen Ctypes_static.ptr
val safeel : gen -> Signed.long -> Signed.long Ctypes_static.ptr
val safelistel : gen -> Signed.long -> gen Ctypes_static.ptr
val safegcoeff : gen -> Signed.long -> Signed.long -> gen Ctypes_static.ptr
val strtoi : string -> gen
val strtor : string -> Signed.long -> gen
val varhigher : string -> Signed.long -> gen
val varlower : string -> Signed.long -> gen
val divisorslenstra : gen -> gen -> gen -> gen
val isprimeaprcl : gen -> Signed.long
val qfb0 : gen -> gen -> gen -> gen

val check_quaddisc :
  gen ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  string ->
  unit

val check_quaddisc_imag : gen -> Signed.long Ctypes_static.ptr -> string -> unit
val check_quaddisc_real : gen -> Signed.long Ctypes_static.ptr -> string -> unit

val cornacchia :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> Signed.long

val cornacchia2 :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> Signed.long

val cornacchia2_sqrt :
  gen ->
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long

val nucomp : gen -> gen -> gen -> gen
val nudupl : gen -> gen -> gen
val nupow : gen -> gen -> gen -> gen
val primeform : gen -> gen -> gen
val primeform_u : gen -> pari_ulong -> gen
val qfb_1 : gen -> gen
val qfbcomp : gen -> gen -> gen
val qfbcomp_i : gen -> gen -> gen
val qfbcompraw : gen -> gen -> gen
val qfbcompraw_i : gen -> gen -> gen
val qfbcornacchia : gen -> gen -> gen
val qfbpow : gen -> gen -> gen
val qfbpow_i : gen -> gen -> gen
val qfbpowraw : gen -> Signed.long -> gen
val qfbpows : gen -> Signed.long -> gen
val qfbred : gen -> gen
val qfbred_i : gen -> gen
val qfbred0 : gen -> Signed.long -> gen -> gen -> gen
val qfbredsl2 : gen -> gen -> gen
val qfbsolve : gen -> gen -> Signed.long -> gen
val qfbsqr : gen -> gen
val qfbsqr_i : gen -> gen
val qfisolvep : gen -> gen -> gen
val qfr3_comp : gen -> gen -> qfr_data Ctypes.structure Ctypes_static.ptr -> gen
val qfr3_compraw : gen -> gen -> gen
val qfr3_pow : gen -> gen -> qfr_data Ctypes.structure Ctypes_static.ptr -> gen
val qfr3_red : gen -> qfr_data Ctypes.structure Ctypes_static.ptr -> gen
val qfr3_rho : gen -> qfr_data Ctypes.structure Ctypes_static.ptr -> gen
val qfr3_to_qfr : gen -> gen -> gen
val qfr5_comp : gen -> gen -> qfr_data Ctypes.structure Ctypes_static.ptr -> gen
val qfr5_compraw : gen -> gen -> gen
val qfr5_dist : gen -> gen -> Signed.long -> gen
val qfr5_pow : gen -> gen -> qfr_data Ctypes.structure Ctypes_static.ptr -> gen
val qfr5_red : gen -> qfr_data Ctypes.structure Ctypes_static.ptr -> gen
val qfr5_rho : gen -> qfr_data Ctypes.structure Ctypes_static.ptr -> gen
val qfr5_to_qfr : gen -> gen -> gen -> gen

val qfr_data_init :
  gen -> Signed.long -> qfr_data Ctypes.structure Ctypes_static.ptr -> unit

val qfr_to_qfr5 : gen -> Signed.long -> gen
val qfrsolvep : gen -> gen -> gen
val quadgen : gen -> gen
val quadgen0 : gen -> Signed.long -> gen
val quadpoly : gen -> gen
val quadpoly_i : gen -> gen
val quadpoly0 : gen -> Signed.long -> gen
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
val fl_powers : pari_ulong -> Signed.long -> pari_ulong -> gen
val fl_powers_pre : pari_ulong -> Signed.long -> pari_ulong -> pari_ulong -> gen
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

val fp_2gener : gen -> gen
val fp_2gener_i : gen -> gen -> gen
val fp_factored_order : gen -> gen -> gen -> gen
val fp_ispower : gen -> gen -> gen -> int
val fp_log : gen -> gen -> gen -> gen -> gen
val fp_order : gen -> gen -> gen -> gen
val fp_pow : gen -> gen -> gen -> gen
val fp_pow_init : gen -> gen -> Signed.long -> gen -> gen
val fp_pow_table : gen -> gen -> gen -> gen
val fp_powers : gen -> Signed.long -> gen -> gen
val fp_pows : gen -> Signed.long -> gen -> gen
val fp_powu : gen -> pari_ulong -> gen -> gen
val fp_sqrt : gen -> gen -> gen
val fp_sqrt_i : gen -> gen -> gen -> gen
val fp_sqrtn : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fpv_prod : gen -> gen -> gen
val z_zv_mod : gen -> gen -> gen
val z_zv_mod_tree : gen -> gen -> gen -> gen
val z_chinese : gen -> gen -> gen -> gen -> gen
val z_chinese_all : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val z_chinese_coprime : gen -> gen -> gen -> gen -> gen -> gen
val z_chinese_post : gen -> gen -> gen -> gen -> gen -> gen

val z_chinese_pre :
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  unit

val z_factor_listp : gen -> gen -> gen
val z_nv_mod : gen -> gen -> gen
val zm_nv_mod_tree : gen -> gen -> gen -> gen
val zv_allpnqn : gen -> gen
val zv_chinese : gen -> gen -> gen Ctypes_static.ptr -> gen
val zv_chinese_tree : gen -> gen -> gen -> gen -> gen
val zv_chinesetree : gen -> gen -> gen
val zv_nv_mod_tree : gen -> gen -> gen -> gen
val zv_producttree : gen -> gen
val zx_nv_mod_tree : gen -> gen -> gen -> gen
val zxc_nv_mod_tree : gen -> gen -> gen -> Signed.long -> gen
val zxm_nv_mod_tree : gen -> gen -> gen -> Signed.long -> gen
val zxx_nv_mod_tree : gen -> gen -> gen -> Signed.long -> gen
val zideallog : gen -> gen -> gen
val bestappr : gen -> gen -> gen
val bestapprpade : gen -> Signed.long -> gen
val chinese : gen -> gen -> gen
val chinese1 : gen -> gen
val chinese1_coprime_z : gen -> gen
val contfrac0 : gen -> gen -> Signed.long -> gen
val contfracpnqn : gen -> Signed.long -> gen
val fibo : Signed.long -> gen
val gboundcf : gen -> Signed.long -> gen
val gcf : gen -> gen
val gcf2 : gen -> gen -> gen

val get_fp_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  bb_field Ctypes.structure Ctypes_static.ptr

val hilbert : gen -> gen -> gen -> Signed.long
val hilbertii : gen -> gen -> gen -> Signed.long
val istotient : gen -> gen Ctypes_static.ptr -> Signed.long
val krois : gen -> Signed.long -> Signed.long
val kroiu : gen -> pari_ulong -> Signed.long
val kronecker : gen -> gen -> Signed.long
val krosi : Signed.long -> gen -> Signed.long
val kross : Signed.long -> Signed.long -> Signed.long
val kroui : pari_ulong -> gen -> Signed.long
val krouu : pari_ulong -> pari_ulong -> Signed.long
val lcmii : gen -> gen -> gen
val fp_invgen : gen -> gen -> gen Ctypes_static.ptr -> gen
val logint0 : gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val logintall : gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val mpfact : Signed.long -> gen
val factorial_fl : Signed.long -> pari_ulong -> pari_ulong
val factorial_fp : Signed.long -> gen -> gen
val muls_interval : Signed.long -> Signed.long -> gen
val mulu_interval : pari_ulong -> pari_ulong -> gen
val mulu_interval_step : pari_ulong -> pari_ulong -> pari_ulong -> gen
val ncv_chinese_center : gen -> gen -> gen Ctypes_static.ptr -> gen
val ncv_chinese_center_tree : gen -> gen -> gen -> gen -> gen
val nmv_chinese_center : gen -> gen -> gen Ctypes_static.ptr -> gen
val nmv_chinese_center_tree : gen -> gen -> gen -> gen -> gen
val nonsquare_fl : pari_ulong -> pari_ulong
val nxcv_chinese_center : gen -> gen -> gen Ctypes_static.ptr -> gen
val nxcv_chinese_center_tree : gen -> gen -> gen -> gen -> gen
val nxmv_chinese_center : gen -> gen -> gen Ctypes_static.ptr -> gen
val nxv_chinese_center : gen -> gen -> gen Ctypes_static.ptr -> gen
val nxv_chinese_center_tree : gen -> gen -> gen -> gen -> gen
val zv_chinese_center : gen -> gen -> gen Ctypes_static.ptr -> gen
val odd_prime_divisors : gen -> gen
val pgener_fl : pari_ulong -> pari_ulong
val pgener_fl_local : pari_ulong -> gen -> pari_ulong
val pgener_fp : gen -> gen
val pgener_fp_local : gen -> gen -> gen
val pgener_zl : pari_ulong -> pari_ulong
val pgener_zp : gen -> gen
val pnqn : gen -> gen
val ramanujantau : gen -> Signed.long -> gen
val rootsof1_fl : pari_ulong -> pari_ulong -> pari_ulong
val rootsof1_fp : gen -> gen -> gen
val rootsof1u_fp : pari_ulong -> gen -> gen

val u_chinese_coprime :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val znlog : gen -> gen -> gen -> gen
val znorder : gen -> gen -> gen
val znprimroot : gen -> gen
val znstar : gen -> gen
val znstar0 : gen -> Signed.long -> gen
val rgv_is_zvpos : gen -> int
val rgv_is_zvnon0 : gen -> int
val rgv_is_prv : gen -> int
val z_issquarefree_fact : gen -> Signed.long

val z_lsmoothen :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val z_smoothen :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val bigomega : gen -> Signed.long
val bigomegau : pari_ulong -> Signed.long
val boundfact : gen -> pari_ulong -> gen
val check_arith_pos : gen -> string -> gen
val check_arith_non0 : gen -> string -> gen
val check_arith_all : gen -> string -> gen
val clean_z_factor : gen -> gen
val core : gen -> gen

val coredisc2_fact :
  gen -> Signed.long -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val coredisc2u_fact :
  gen ->
  Signed.long ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  pari_ulong

val corepartial : gen -> Signed.long -> gen
val core0 : gen -> Signed.long -> gen
val core2 : gen -> gen
val core2partial : gen -> Signed.long -> gen
val coredisc : gen -> gen
val coredisc0 : gen -> Signed.long -> gen
val coredisc2 : gen -> gen
val corediscs : Signed.long -> pari_ulong Ctypes_static.ptr -> Signed.long
val divisors : gen -> gen
val divisors_factored : gen -> gen
val divisors0 : gen -> Signed.long -> gen
val divisorsu : pari_ulong -> gen
val divisorsu_moebius : gen -> gen
val divisorsu_fact : gen -> gen
val divisorsu_fact_factored : gen -> gen
val eulerphi : gen -> gen
val eulerphiu : pari_ulong -> pari_ulong
val eulerphiu_fact : gen -> pari_ulong
val factor_pn_1 : gen -> pari_ulong -> gen
val factor_pn_1_limit : gen -> Signed.long -> pari_ulong -> gen
val factoru_pow : pari_ulong -> gen
val fuse_z_factor : gen -> gen -> gen
val is_z_factor : gen -> int
val is_z_factornon0 : gen -> int
val is_z_factorpos : gen -> int
val is_nf_factor : gen -> int
val is_nf_extfactor : gen -> int
val issquarefree : gen -> Signed.long
val numdiv : gen -> gen
val numdivu : Signed.long -> Signed.long
val numdivu_fact : gen -> Signed.long
val omega : gen -> Signed.long
val omegau : pari_ulong -> Signed.long
val sumdiv : gen -> gen
val sumdivk : gen -> Signed.long -> gen
val uissquarefree : pari_ulong -> Signed.long
val uissquarefree_fact : gen -> Signed.long
val usumdiv_fact : gen -> gen
val usumdivk_fact : gen -> pari_ulong -> gen
val fpx_fpc_nfpoleval : gen -> gen -> gen -> gen -> gen
val embed_t2 : gen -> Signed.long -> gen
val embednorm_t2 : gen -> Signed.long -> gen
val embed_norm : gen -> Signed.long -> gen
val check_zkmodule_i : gen -> int
val check_zkmodule : gen -> string -> unit
val checkbid : gen -> unit
val checkbid_i : gen -> gen
val checkbnf : gen -> gen
val checkbnf_i : gen -> gen
val checkbnr : gen -> unit
val checkbnr_i : gen -> gen
val checkabgrp : gen -> unit
val checksqmat : gen -> Signed.long -> unit
val checknf : gen -> gen
val checknf_i : gen -> gen
val checknfelt_mod : gen -> gen -> string -> gen
val checkprid : gen -> unit
val checkprid_i : gen -> int
val checkrnf : gen -> unit
val checkrnf_i : gen -> int
val factoredpolred : gen -> gen -> gen
val factoredpolred2 : gen -> gen -> gen
val galoisapply : gen -> gen -> gen -> gen
val get_bnf : gen -> Signed.long Ctypes_static.ptr -> gen
val get_bnfpol : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen
val get_nf : gen -> Signed.long Ctypes_static.ptr -> gen
val get_nfpol : gen -> gen Ctypes_static.ptr -> gen
val get_prid : gen -> gen
val idealfrobenius : gen -> gen -> gen -> gen
val idealfrobenius_aut : gen -> gen -> gen -> gen -> gen
val idealramfrobenius : gen -> gen -> gen -> gen -> gen
val idealramfrobenius_aut : gen -> gen -> gen -> gen -> gen -> gen
val idealramgroups : gen -> gen -> gen -> gen
val idealramgroups_aut : gen -> gen -> gen -> gen -> gen
val nf_get_allroots : gen -> gen
val nf_get_prec : gen -> Signed.long

val nfmaxord_to_nf :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr -> gen -> Signed.long -> gen

val nfcertify : gen -> gen
val nfgaloismatrix : gen -> gen -> gen
val nfgaloismatrixapply : gen -> gen -> gen -> gen
val nfgaloispermtobasis : gen -> gen -> gen
val nfinit_basic : nfmaxord_t Ctypes.structure Ctypes_static.ptr -> gen -> unit

val nfinit_complete :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  gen

val nfinit : gen -> Signed.long -> gen
val nfinit0 : gen -> Signed.long -> Signed.long -> gen
val nfinitred : gen -> Signed.long -> gen
val nfinitred2 : gen -> Signed.long -> gen
val nfisincl : gen -> gen -> gen
val nfisincl0 : gen -> gen -> Signed.long -> gen
val nfisisom : gen -> gen -> gen
val nfnewprec : gen -> Signed.long -> gen
val nfnewprec_shallow : gen -> Signed.long -> gen
val nfpoleval : gen -> gen -> gen -> gen
val nfsplitting : gen -> gen -> gen
val nfsplitting0 : gen -> gen -> Signed.long -> gen
val nftyp : gen -> Signed.long
val polredord : gen -> gen
val polred : gen -> gen
val polred0 : gen -> Signed.long -> gen -> gen
val polred2 : gen -> gen
val polredabs : gen -> gen
val polredabs0 : gen -> Signed.long -> gen
val polredabs2 : gen -> gen
val polredabsall : gen -> Signed.long -> gen
val polredbest : gen -> Signed.long -> gen
val poltomonic : gen -> gen Ctypes_static.ptr -> gen
val rnfpolredabs : gen -> gen -> Signed.long -> gen
val rnfpolredbest : gen -> gen -> Signed.long -> gen
val smallpolred : gen -> gen
val smallpolred2 : gen -> gen
val tschirnhaus : gen -> gen
val zx_q_mul : gen -> gen -> gen
val zx_q_normalize : gen -> gen Ctypes_static.ptr -> gen
val zx_z_normalize : gen -> gen Ctypes_static.ptr -> gen
val zx_to_monic : gen -> gen Ctypes_static.ptr -> gen
val zx_primitive_to_monic : gen -> gen Ctypes_static.ptr -> gen
val zxx_q_mul : gen -> gen -> gen
val fq_to_nf : gen -> gen -> gen
val fqm_to_nfm : gen -> gen -> gen
val fqv_to_nfv : gen -> gen -> gen
val fqx_to_nfx : gen -> gen -> gen
val rg_nffix : string -> gen -> gen -> int -> gen
val rgv_nffix : string -> gen -> gen -> int -> gen
val rgx_nffix : string -> gen -> gen -> int -> gen
val zx_composedsum : gen -> gen -> gen
val zx_compositum : gen -> gen -> Signed.long Ctypes_static.ptr -> gen
val zpx_disc_val : gen -> gen -> Signed.long
val zpx_gcd : gen -> gen -> gen -> gen -> gen
val zpx_monic_factor : gen -> gen -> Signed.long -> gen
val zpx_primedec : gen -> gen -> gen
val zpx_reduced_resultant : gen -> gen -> gen -> gen -> gen
val zpx_reduced_resultant_fast : gen -> gen -> gen -> Signed.long -> gen
val zpx_resultant_val : gen -> gen -> gen -> Signed.long -> Signed.long
val checkmodpr : gen -> unit
val compositum : gen -> gen -> gen
val compositum2 : gen -> gen -> gen
val nfdisc : gen -> gen
val get_modpr : gen -> gen
val indexpartial : gen -> gen -> gen
val modpr_genfq : gen -> gen

val nf_to_fq_init :
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val nf_to_fq : gen -> gen -> gen -> gen
val nfm_to_fqm : gen -> gen -> gen -> gen
val nfv_to_fqv : gen -> gen -> gen -> gen
val nfx_to_fqx : gen -> gen -> gen -> gen
val nfx_to_monic : gen -> gen -> gen Ctypes_static.ptr -> gen
val nfbasis : gen -> gen Ctypes_static.ptr -> gen
val nfcompositum : gen -> gen -> gen -> Signed.long -> gen
val nfdiscfactors : gen -> gen

val nfmaxord :
  nfmaxord_t Ctypes.structure Ctypes_static.ptr -> gen -> Signed.long -> unit

val nfmodpr : gen -> gen -> gen -> gen
val nfmodprinit : gen -> gen -> gen
val nfmodprinit0 : gen -> gen -> Signed.long -> gen
val nfmodprlift : gen -> gen -> gen -> gen
val nfreducemodpr : gen -> gen -> gen -> gen
val polcompositum0 : gen -> gen -> Signed.long -> gen
val idealprimedec : gen -> gen -> gen
val idealprimedec_galois : gen -> gen -> gen
val idealprimedec_degrees : gen -> gen -> gen
val idealprimedec_kummer : gen -> gen -> Signed.long -> gen -> gen
val idealprimedec_limit_f : gen -> gen -> Signed.long -> gen
val idealprimedec_limit_norm : gen -> gen -> gen -> gen
val poldiscfactors : gen -> Signed.long -> gen
val rnfbasis : gen -> gen -> gen
val rnfdedekind : gen -> gen -> gen -> Signed.long -> gen
val rnfdet : gen -> gen -> gen
val rnfdisc_factored : gen -> gen -> gen Ctypes_static.ptr -> gen
val rnfdiscf : gen -> gen -> gen
val rnfequation : gen -> gen -> gen
val rnfequation0 : gen -> gen -> Signed.long -> gen
val rnfequation2 : gen -> gen -> gen
val nf_pv_to_prv : gen -> gen -> gen
val nf_rnfeq : gen -> gen -> gen
val nf_rnfeqsimple : gen -> gen -> gen

val rnfequationall :
  gen -> gen -> Signed.long Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val rnfhnfbasis : gen -> gen -> gen
val rnfisfree : gen -> gen -> Signed.long
val rnflllgram : gen -> gen -> gen -> Signed.long -> gen
val rnfpolred : gen -> gen -> Signed.long -> gen
val rnfpseudobasis : gen -> gen -> gen
val rnfsimplifybasis : gen -> gen -> gen
val rnfsteinitz : gen -> gen -> gen
val factorial_lval : pari_ulong -> pari_ulong -> Signed.long

val zk_to_fq_init :
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val zk_to_fq : gen -> gen -> gen
val qxqv_to_fpm : gen -> gen -> gen -> gen
val zkmodprinit : gen -> gen -> gen
val idealstar : gen -> gen -> Signed.long -> gen
val idealstarprk : gen -> gen -> Signed.long -> Signed.long -> gen
val rgc_to_nfc : gen -> gen -> gen
val rgm_rgx_mul : gen -> gen -> gen
val rgm_to_nfm : gen -> gen -> gen
val rgx_to_nfx : gen -> gen -> gen
val zc_nfval : gen -> gen -> Signed.long
val zc_nfvalrem : gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val zc_prdvd : gen -> gen -> int
val zm_zx_mul : gen -> gen -> gen
val zv_snf_gcd : gen -> gen -> gen
val algtobasis : gen -> gen -> gen
val basistoalg : gen -> gen -> gen
val ei_multable : gen -> Signed.long -> gen

val get_nf_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  bb_field Ctypes.structure Ctypes_static.ptr

val famat_nfvalrem : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val gpnfvalrem : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val idealfactorback : gen -> gen -> gen -> int -> gen
val ideallist : gen -> Signed.long -> gen
val ideallist0 : gen -> Signed.long -> Signed.long -> gen
val gideallist : gen -> gen -> Signed.long -> gen
val ideallistarch : gen -> gen -> gen -> gen
val ideallog : gen -> gen -> gen -> gen
val ideallogmod : gen -> gen -> gen -> gen -> gen
val ideallog_units : gen -> gen -> gen
val ideallog_units0 : gen -> gen -> gen -> gen
val idealprincipalunits : gen -> gen -> Signed.long -> gen
val idealstar0 : gen -> gen -> Signed.long -> gen
val idealstarmod : gen -> gen -> Signed.long -> gen -> gen
val indices_to_vec01 : gen -> Signed.long -> gen
val matalgtobasis : gen -> gen -> gen
val matbasistoalg : gen -> gen -> gen
val multable : gen -> gen -> gen
val nf_to_scalar_or_alg : gen -> gen -> gen
val nf_to_scalar_or_basis : gen -> gen -> gen
val nf_cxlog : gen -> gen -> Signed.long -> gen
val nfv_cxlog : gen -> gen -> Signed.long -> gen
val nfadd : gen -> gen -> gen -> gen
val nfchecksigns : gen -> gen -> gen -> int
val nfdiv : gen -> gen -> gen -> gen
val nfdiveuc : gen -> gen -> gen -> gen
val nfdivrem : gen -> gen -> gen -> gen
val nfembed : gen -> gen -> Signed.long -> gen
val nfeltembed : gen -> gen -> gen -> Signed.long -> gen
val nfeltembed_i : gen Ctypes_static.ptr -> gen -> gen -> Signed.long -> gen
val nfeltsign : gen -> gen -> gen -> gen
val nffactorback : gen -> gen -> gen -> gen
val nfinv : gen -> gen -> gen
val nfinvmodideal : gen -> gen -> gen -> gen
val nfissquare : gen -> gen -> gen Ctypes_static.ptr -> Signed.long

val nfispower :
  gen -> gen -> Signed.long -> gen Ctypes_static.ptr -> Signed.long

val nflogembed : gen -> gen -> gen Ctypes_static.ptr -> Signed.long -> gen
val nfm_det : gen -> gen -> gen
val nfm_inv : gen -> gen -> gen
val nfm_ker : gen -> gen -> gen
val nfm_mul : gen -> gen -> gen -> gen
val nfm_nfc_mul : gen -> gen -> gen -> gen
val nfmod : gen -> gen -> gen -> gen
val nfmul : gen -> gen -> gen -> gen
val nfmuli : gen -> gen -> gen -> gen
val nfnorm : gen -> gen -> gen
val nfpolsturm : gen -> gen -> gen -> gen
val nfpow : gen -> gen -> gen -> gen
val nfpow_u : gen -> gen -> pari_ulong -> gen
val nfpowmodideal : gen -> gen -> gen -> gen -> gen
val nfsign : gen -> gen -> gen
val nfsign_arch : gen -> gen -> gen -> gen
val nfsign_from_logarch : gen -> gen -> gen -> gen
val nfsqr : gen -> gen -> gen
val nfsqri : gen -> gen -> gen
val nfsub : gen -> gen -> gen -> gen
val nftrace : gen -> gen -> gen
val nfval : gen -> gen -> gen -> Signed.long
val nfvalrem : gen -> gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val polmod_nffix : string -> gen -> gen -> int -> gen
val polmod_nffix2 : string -> gen -> gen -> gen -> int -> gen
val pr_basis_perm : gen -> gen -> gen
val pr_equal : gen -> gen -> int
val rnfalgtobasis : gen -> gen -> gen
val rnfbasistoalg : gen -> gen -> gen
val rnfeltnorm : gen -> gen -> gen
val rnfelttrace : gen -> gen -> gen
val set_sign_mod_divisor : gen -> gen -> gen -> gen -> gen
val tablemul : gen -> gen -> gen -> gen
val tablemul_ei : gen -> gen -> Signed.long -> gen
val tablemul_ei_ej : gen -> Signed.long -> Signed.long -> gen
val tablemulvec : gen -> gen -> gen -> gen
val tablesqr : gen -> gen -> gen
val vec01_to_indices : gen -> gen
val vecsmall01_to_indices : gen -> gen
val zk_inv : gen -> gen -> gen
val zk_multable : gen -> gen -> gen
val zk_scalar_or_multable : gen -> gen -> gen
val zkchinese : gen -> gen -> gen -> gen
val zkchinese1 : gen -> gen -> gen
val zkchineseinit : gen -> gen -> gen -> gen -> gen
val zkmultable_capz : gen -> gen
val zkmultable_inv : gen -> gen

val fl_invgen :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val rm_round_maxrank : gen -> gen
val zm_famat_limit : gen -> gen -> gen
val zv_cba : gen -> gen
val zv_cba_extend : gen -> gen -> gen
val z_cba : gen -> gen -> gen
val z_ppgle : gen -> gen -> gen
val z_ppio : gen -> gen -> gen
val z_ppo : gen -> gen -> gen
val famatv_factorback : gen -> gen -> gen
val famatv_zv_factorback : gen -> gen -> gen
val famat_z_gcd : gen -> gen -> gen
val famat_div : gen -> gen -> gen
val famat_div_shallow : gen -> gen -> gen
val famat_idealfactor : gen -> gen -> gen
val famat_inv : gen -> gen
val famat_inv_shallow : gen -> gen
val famat_makecoprime : gen -> gen -> gen -> gen -> gen -> gen -> gen
val famat_mul : gen -> gen -> gen
val famat_mul_shallow : gen -> gen -> gen
val famat_mulpow_shallow : gen -> gen -> gen -> gen
val famat_mulpows_shallow : gen -> gen -> Signed.long -> gen
val famat_pow : gen -> gen -> gen
val famat_pow_shallow : gen -> gen -> gen
val famat_pows_shallow : gen -> Signed.long -> gen
val famat_reduce : gen -> gen
val famat_remove_trivial : gen -> gen
val famat_sqr : gen -> gen
val famat_to_nf : gen -> gen -> gen
val famat_to_nf_moddivisor : gen -> gen -> gen -> gen -> gen
val famat_to_nf_modideal_coprime : gen -> gen -> gen -> gen -> gen -> gen
val famatsmall_reduce : gen -> gen
val gpidealfactor : gen -> gen -> gen -> gen
val gpidealval : gen -> gen -> gen -> gen

val idealhnf_z_factor :
  gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val idealhnf_z_factor_i :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val idealhnf_inv : gen -> gen -> gen
val idealhnf_inv_z : gen -> gen -> gen
val idealhnf_mul : gen -> gen -> gen -> gen
val idealadd : gen -> gen -> gen -> gen
val idealaddmultoone : gen -> gen -> gen
val idealaddtoone : gen -> gen -> gen -> gen
val idealaddtoone0 : gen -> gen -> gen -> gen
val idealaddtoone_i : gen -> gen -> gen -> gen
val idealaddtoone_raw : gen -> gen -> gen -> gen
val idealappr : gen -> gen -> gen
val idealappr0 : gen -> gen -> Signed.long -> gen
val idealapprfact : gen -> gen -> gen
val idealchinese : gen -> gen -> gen -> gen
val idealcoprime : gen -> gen -> gen -> gen
val idealcoprimefact : gen -> gen -> gen -> gen
val idealdiv : gen -> gen -> gen -> gen
val idealdiv0 : gen -> gen -> gen -> Signed.long -> gen
val idealdivexact : gen -> gen -> gen -> gen
val idealdivpowprime : gen -> gen -> gen -> gen -> gen
val idealdown : gen -> gen -> gen
val idealfactor : gen -> gen -> gen
val idealfactor_limit : gen -> gen -> pari_ulong -> gen
val idealfactor_partial : gen -> gen -> gen -> gen
val idealhnf : gen -> gen -> gen
val idealhnf0 : gen -> gen -> gen -> gen
val idealhnf_principal : gen -> gen -> gen
val idealhnf_shallow : gen -> gen -> gen
val idealhnf_two : gen -> gen -> gen
val idealintersect : gen -> gen -> gen -> gen
val idealinv : gen -> gen -> gen
val idealismaximal : gen -> gen -> gen

val idealispower :
  gen -> gen -> Signed.long -> gen Ctypes_static.ptr -> Signed.long

val idealmin : gen -> gen -> gen -> gen
val idealmul : gen -> gen -> gen -> gen
val idealmul0 : gen -> gen -> gen -> Signed.long -> gen
val idealmulpowprime : gen -> gen -> gen -> gen -> gen
val idealmulred : gen -> gen -> gen -> gen
val idealnorm : gen -> gen -> gen
val idealnumden : gen -> gen -> gen
val idealpow : gen -> gen -> gen -> gen
val idealpow0 : gen -> gen -> gen -> Signed.long -> gen
val idealpowred : gen -> gen -> gen -> gen
val idealpows : gen -> gen -> Signed.long -> gen
val idealprod : gen -> gen -> gen
val idealprodprime : gen -> gen -> gen
val idealprodval : gen -> gen -> gen -> Signed.long
val idealpseudomin : gen -> gen -> gen
val idealpseudomin_nonscalar : gen -> gen -> gen
val idealpseudominvec : gen -> gen -> gen
val idealpseudored : gen -> gen -> gen
val idealred0 : gen -> gen -> gen -> gen
val idealred_elt : gen -> gen -> gen
val idealredmodpower : gen -> gen -> pari_ulong -> pari_ulong -> gen
val idealsqr : gen -> gen -> gen
val idealtwoelt : gen -> gen -> gen
val idealtwoelt0 : gen -> gen -> gen -> gen
val idealtwoelt2 : gen -> gen -> gen -> gen
val idealtyp : gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> Signed.long
val idealval : gen -> gen -> gen -> Signed.long
val isideal : gen -> gen -> Signed.long
val matreduce : gen -> gen
val nfc_multable_mul : gen -> gen -> gen
val nfc_nf_mul : gen -> gen -> gen -> gen
val nf_get_gtwist : gen -> gen -> gen
val nf_get_gtwist1 : gen -> Signed.long -> gen
val nf_to_fp_coprime : gen -> gen -> gen -> gen
val nfdetint : gen -> gen -> gen
val nfdivmodpr : gen -> gen -> gen -> gen -> gen
val nfhnf : gen -> gen -> gen
val nfhnf0 : gen -> gen -> Signed.long -> gen
val nfhnfmod : gen -> gen -> gen -> gen
val nfkermodpr : gen -> gen -> gen -> gen
val nfmulmodpr : gen -> gen -> gen -> gen -> gen
val nfpowmodpr : gen -> gen -> gen -> gen -> gen
val nfreduce : gen -> gen -> gen -> gen
val nfsnf : gen -> gen -> gen
val nfsnf0 : gen -> gen -> Signed.long -> gen
val nfsolvemodpr : gen -> gen -> gen -> gen -> gen
val prv_lcm_capz : gen -> gen
val prv_primes : gen -> gen
val pr_hnf : gen -> gen -> gen
val pr_inv : gen -> gen
val pr_inv_p : gen -> gen
val pr_uniformizer : gen -> gen -> gen
val sunits_makecoprime : gen -> gen -> gen -> gen
val to_famat : gen -> gen -> gen
val to_famat_shallow : gen -> gen -> gen
val u_ppo : pari_ulong -> pari_ulong -> pari_ulong
val vecdiv : gen -> gen -> gen
val vecinv : gen -> gen
val vecmul : gen -> gen -> gen
val vecpow : gen -> gen -> gen
val vecsqr : gen -> gen
val zkc_multable_mul : gen -> gen -> gen
val eltreltoabs : gen -> gen -> gen
val eltabstorel : gen -> gen -> gen
val eltabstorel_lift : gen -> gen -> gen
val nf_nfzk : gen -> gen -> gen
val rnf_build_nfabs : gen -> Signed.long -> gen
val rnf_zkabs : gen -> gen
val nfeltup : gen -> gen -> gen -> gen
val rnfcomplete : gen -> unit
val rnfeltabstorel : gen -> gen -> gen
val rnfeltdown : gen -> gen -> gen
val rnfeltdown0 : gen -> gen -> Signed.long -> gen
val rnfeltreltoabs : gen -> gen -> gen
val rnfeltup : gen -> gen -> gen
val rnfeltup0 : gen -> gen -> Signed.long -> gen
val rnfidealabstorel : gen -> gen -> gen
val rnfidealdown : gen -> gen -> gen
val rnfidealfactor : gen -> gen -> gen
val rnfidealhnf : gen -> gen -> gen
val rnfidealmul : gen -> gen -> gen -> gen
val rnfidealnormabs : gen -> gen -> gen
val rnfidealnormrel : gen -> gen -> gen
val rnfidealprimedec : gen -> gen -> gen
val rnfidealreltoabs : gen -> gen -> gen
val rnfidealreltoabs0 : gen -> gen -> Signed.long -> gen
val rnfidealtwoelement : gen -> gen -> gen
val rnfidealup : gen -> gen -> gen
val rnfidealup0 : gen -> gen -> Signed.long -> gen
val rnfinit : gen -> gen -> gen
val rnfinit0 : gen -> gen -> Signed.long -> gen
val get_arith_zzm : gen -> gen
val get_arith_z : gen -> gen

val gen_ph_log :
  gen ->
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_shanks_init :
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_shanks :
  gen ->
  gen ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_shanks_sqrtn :
  gen ->
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_gener :
  gen ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_ellgens :
  gen ->
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen -> gen -> gen)
  Ctypes_static.static_funptr ->
  gen

val gen_ellgroup :
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen -> gen -> gen)
  Ctypes_static.static_funptr ->
  gen

val gen_factored_order :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_order :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_select_order :
  gen ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_plog :
  gen ->
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_group Ctypes.structure Ctypes_static.ptr ->
  gen

val gen_pow :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_pow_fold :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_pow_fold_i :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_pow_i :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_pow_init :
  gen ->
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_pow_table :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_powers :
  gen ->
  Signed.long ->
  int ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen) Ctypes_static.static_funptr ->
  gen

val gen_powu :
  gen ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_powu_fold :
  gen ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_powu_fold_i :
  gen ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_powu_i :
  gen ->
  pari_ulong ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val gen_product :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  gen

val matdetmod : gen -> gen -> gen
val matimagemod : gen -> gen -> gen Ctypes_static.ptr -> gen
val matinvmod : gen -> gen -> gen
val matkermod : gen -> gen -> gen Ctypes_static.ptr -> gen
val matsolvemod : gen -> gen -> gen -> Signed.long -> gen
val bernfrac : Signed.long -> gen
val bernpol : Signed.long -> Signed.long -> gen
val bernreal : Signed.long -> Signed.long -> gen
val bernvec : Signed.long -> gen
val constbern : Signed.long -> unit
val eulerfrac : Signed.long -> gen
val eulerpol : Signed.long -> Signed.long -> gen
val eulerreal : Signed.long -> Signed.long -> gen
val eulervec : Signed.long -> gen
val harmonic : pari_ulong -> gen
val harmonic0 : pari_ulong -> gen -> gen

val qr_init :
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long ->
  int

val r_from_qr : gen -> Signed.long -> gen
val rgm_babai : gen -> gen -> gen

val rgm_qr_init :
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long ->
  int

val rgm_gram_schmidt : gen -> gen Ctypes_static.ptr -> gen
val algdep : gen -> Signed.long -> gen
val algdep0 : gen -> Signed.long -> Signed.long -> gen
val bestapprnf : gen -> gen -> gen -> Signed.long -> gen

val forqfvec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> float -> Signed.long)
  Ctypes_static.static_funptr ->
  gen ->
  gen ->
  unit

val forqfvec0 : gen -> gen -> gen -> unit

val forqfvec1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  unit

val gaussred_from_qr : gen -> Signed.long -> gen
val lindep : gen -> gen
val lindep_xadic : gen -> gen
val lindep_bit : gen -> Signed.long -> gen
val lindep_padic : gen -> gen
val lindep0 : gen -> Signed.long -> gen
val lindep2 : gen -> Signed.long -> gen
val lindepfull_bit : gen -> Signed.long -> gen
val mathouseholder : gen -> gen -> gen
val matqr : gen -> Signed.long -> Signed.long -> gen
val minim : gen -> gen -> gen -> gen
val minim_raw : gen -> gen -> gen -> gen
val minim_zm : gen -> gen -> gen -> gen
val minim2 : gen -> gen -> gen -> gen
val qfminim0 : gen -> gen -> gen -> Signed.long -> Signed.long -> gen
val qfperfection : gen -> gen
val qfrep0 : gen -> gen -> Signed.long -> gen
val seralgdep : gen -> Signed.long -> Signed.long -> gen
val serdiffdep : gen -> Signed.long -> Signed.long -> gen
val vandermondeinverse : gen -> gen -> gen -> gen -> gen
val vandermondeinverseinit : gen -> gen
val zncoppersmith : gen -> gen -> gen -> gen -> gen
val qxq_reverse : gen -> gen -> gen
val vec_equiv : gen -> gen
val rgv_polint : gen -> gen -> Signed.long -> gen
val vec_reduce : gen -> gen Ctypes_static.ptr -> gen
val rgxq_reverse : gen -> gen -> gen
val zc_union_shallow : gen -> gen -> gen
val zv_indexsort : gen -> gen
val zv_sort : gen -> gen
val zv_sort_inplace : gen -> unit
val zv_sort_shallow : gen -> gen
val zv_sort_uniq : gen -> gen
val zv_sort_uniq_shallow : gen -> gen
val zv_union_shallow : gen -> gen -> gen
val binomial : gen -> Signed.long -> gen
val binomial0 : gen -> gen -> gen
val binomialuu : pari_ulong -> pari_ulong -> gen
val cmp_flx : gen -> gen -> int
val cmp_rgx : gen -> gen -> int
val cmp_nodata : unit Ctypes_static.ptr -> gen -> gen -> int
val cmp_prime_ideal : gen -> gen -> int
val cmp_prime_over_p : gen -> gen -> int
val cmp_universal : gen -> gen -> int
val convol : gen -> gen -> gen
val gen_cmp_rgx : unit Ctypes_static.ptr -> gen -> gen -> int
val polcyclo : Signed.long -> Signed.long -> gen
val polcyclo_eval : Signed.long -> gen -> gen
val dirdiv : gen -> gen -> gen
val dirmul : gen -> gen -> gen
val eulerianpol : Signed.long -> Signed.long -> gen
val gprec_wensure : gen -> Signed.long -> gen

val gen_indexsort :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  gen

val gen_indexsort_uniq :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  gen

val gen_search :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  Signed.long

val gen_setminus :
  gen -> gen -> (gen -> gen -> int) Ctypes_static.static_funptr -> gen

val gen_sort :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  gen

val gen_sort_inplace :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  gen Ctypes_static.ptr ->
  unit

val gen_sort_shallow :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  gen

val gen_sort_uniq :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  gen

val getstack : unit -> Signed.long
val gettime : unit -> Signed.long
val getabstime : unit -> Signed.long
val getwalltime : unit -> gen
val gprec : gen -> Signed.long -> gen
val gprec_wtrunc : gen -> Signed.long -> gen
val gprec_w : gen -> Signed.long -> gen
val gtoset : gen -> gen
val indexlexsort : gen -> gen
val indexsort : gen -> gen
val indexvecsort : gen -> gen -> gen
val laplace : gen -> gen
val lexsort : gen -> gen
val mathilbert : Signed.long -> gen
val matqpascal : Signed.long -> gen -> gen

val merge_factor :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  gen

val merge_sort_uniq :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  gen

val modreverse : gen -> gen
val polhermite : Signed.long -> Signed.long -> gen
val polhermite_eval0 : Signed.long -> gen -> Signed.long -> gen
val polhermite_eval : Signed.long -> gen -> gen
val pollaguerre : Signed.long -> gen -> Signed.long -> gen
val pollaguerre_eval : Signed.long -> gen -> gen -> gen
val pollaguerre_eval0 : Signed.long -> gen -> gen -> Signed.long -> gen
val pollegendre : Signed.long -> Signed.long -> gen
val pollegendre_reduced : Signed.long -> Signed.long -> gen
val pollegendre_eval : Signed.long -> gen -> gen
val pollegendre_eval0 : Signed.long -> gen -> Signed.long -> gen
val polint : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val polint_i : gen -> gen -> gen -> Signed.long Ctypes_static.ptr -> gen

val polintspec :
  gen -> gen -> gen -> Signed.long -> Signed.long Ctypes_static.ptr -> gen

val polchebyshev : Signed.long -> Signed.long -> Signed.long -> gen
val polchebyshev_eval : Signed.long -> Signed.long -> gen -> gen
val polchebyshev1 : Signed.long -> Signed.long -> gen
val polchebyshev2 : Signed.long -> Signed.long -> gen
val polrecip : gen -> gen
val setbinop : gen -> gen -> gen -> gen
val setdelta : gen -> gen -> gen
val setintersect : gen -> gen -> gen
val setisset : gen -> Signed.long
val setminus : gen -> gen -> gen
val setsearch : gen -> gen -> Signed.long -> Signed.long
val setunion : gen -> gen -> gen
val setunion_i : gen -> gen -> gen
val sort : gen -> gen

val sort_factor :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> int) Ctypes_static.static_funptr ->
  gen

val stirling : Signed.long -> Signed.long -> Signed.long -> gen
val stirling1 : pari_ulong -> pari_ulong -> gen
val stirling2 : pari_ulong -> pari_ulong -> gen

val tablesearch :
  gen -> gen -> (gen -> gen -> int) Ctypes_static.static_funptr -> Signed.long

val vecbinomial : Signed.long -> gen
val vecsearch : gen -> gen -> gen -> Signed.long
val vecsort : gen -> gen -> gen
val vecsort0 : gen -> gen -> Signed.long -> gen
val zv_search : gen -> Signed.long -> Signed.long
val bits_to_int : gen -> Signed.long -> gen
val bits_to_u : gen -> Signed.long -> pari_ulong
val binaire : gen -> gen
val binary_2k : gen -> Signed.long -> gen
val binary_2k_nv : gen -> Signed.long -> gen
val binary_zv : gen -> gen
val bittest : gen -> Signed.long -> Signed.long
val fromdigits_2k : gen -> Signed.long -> gen
val gbitand : gen -> gen -> gen
val gbitneg : gen -> Signed.long -> gen
val gbitnegimply : gen -> gen -> gen
val gbitor : gen -> gen -> gen
val gbittest : gen -> Signed.long -> gen
val gbitxor : gen -> gen -> gen
val hammingl : pari_ulong -> Signed.long
val hammingweight : gen -> Signed.long
val ibitand : gen -> gen -> gen
val ibitnegimply : gen -> gen -> gen
val ibitor : gen -> gen -> gen
val ibitxor : gen -> gen -> gen
val nv_fromdigits_2k : gen -> Signed.long -> gen
val bnflogef : gen -> gen -> gen
val bnflog : gen -> gen -> gen
val bnflogdegree : gen -> gen -> gen -> gen
val nfislocalpower : gen -> gen -> gen -> gen -> Signed.long
val rnfislocalcyclo : gen -> Signed.long
val bnfisunit : gen -> gen -> gen
val bnfissunit : gen -> gen -> gen -> gen
val bnfsunit : gen -> gen -> Signed.long -> gen
val bnfunits : gen -> gen -> gen
val bnfisunit0 : gen -> gen -> gen -> gen
val sunits_mod_units : gen -> gen -> gen
val buchquad : gen -> float -> float -> Signed.long -> gen
val quadclassno : gen -> gen
val quadclassnos : Signed.long -> Signed.long
val quadclassunit0 : gen -> Signed.long -> gen -> Signed.long -> gen
val buchall : gen -> Signed.long -> Signed.long -> gen

val buchall_param :
  gen -> float -> float -> Signed.long -> Signed.long -> Signed.long -> gen

val bnf_build_cheapfu : gen -> gen
val bnf_build_cycgen : gen -> gen
val bnf_build_matalpha : gen -> gen
val bnf_build_units : gen -> gen
val bnf_compactfu : gen -> gen
val bnf_compactfu_mat : gen -> gen
val bnf_has_fu : gen -> gen
val bnfinit0 : gen -> Signed.long -> gen -> Signed.long -> gen
val bnfisprincipal0 : gen -> gen -> Signed.long -> gen
val bnfnewprec : gen -> Signed.long -> gen
val bnfnewprec_shallow : gen -> Signed.long -> gen
val bnftestprimes : gen -> gen -> unit
val bnrnewprec : gen -> Signed.long -> gen
val bnrnewprec_shallow : gen -> Signed.long -> gen
val isprincipalfact : gen -> gen -> gen -> gen -> Signed.long -> gen
val isprincipalfact_or_fail : gen -> gen -> gen -> gen -> gen
val isprincipal : gen -> gen -> gen
val nf_cxlog_normalize : gen -> gen -> Signed.long -> gen
val nfcyclotomicunits : gen -> gen -> gen
val nfsign_units : gen -> gen -> int -> gen
val nfsign_tu : gen -> gen -> gen
val nfsign_fu : gen -> gen -> gen
val signunits : gen -> gen
val hermite_bound : Signed.long -> Signed.long -> gen

val bnr_subgroup_sanitize :
  gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit

val bnr_char_sanitize : gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit
val abc_to_bnr : gen -> gen -> gen -> gen Ctypes_static.ptr -> int -> gen
val buchray : gen -> gen -> Signed.long -> gen
val buchraymod : gen -> gen -> Signed.long -> gen -> gen
val bnrautmatrix : gen -> gen -> gen
val bnr_subgroup_check : gen -> gen -> gen Ctypes_static.ptr -> gen
val bnrchar : gen -> gen -> gen -> gen
val bnrchar_primitive : gen -> gen -> gen -> gen
val bnrclassno : gen -> gen -> gen
val bnrclassno0 : gen -> gen -> gen -> gen
val bnrclassnolist : gen -> gen -> gen
val bnrchar_primitive_raw : gen -> gen -> gen -> gen
val bnrconductor_factored : gen -> gen -> gen
val bnrconductor_raw : gen -> gen -> gen
val bnrconductormod : gen -> gen -> gen -> gen
val bnrconductor0 : gen -> gen -> gen -> Signed.long -> gen
val bnrconductor : gen -> gen -> Signed.long -> gen
val bnrconductor_i : gen -> gen -> Signed.long -> gen
val bnrconductorofchar : gen -> gen -> gen
val bnrdisc0 : gen -> gen -> gen -> Signed.long -> gen
val bnrdisc : gen -> gen -> Signed.long -> gen
val bnrdisclist0 : gen -> gen -> gen -> gen
val bnrgaloismatrix : gen -> gen -> gen
val bnrgaloisapply : gen -> gen -> gen -> gen
val bnrinit0 : gen -> gen -> Signed.long -> gen
val bnrinitmod : gen -> gen -> Signed.long -> gen -> gen
val bnrisconductor0 : gen -> gen -> gen -> Signed.long
val bnrisconductor : gen -> gen -> Signed.long
val bnrisgalois : gen -> gen -> gen -> Signed.long
val bnrisprincipalmod : gen -> gen -> gen -> Signed.long -> gen
val bnrisprincipal : gen -> gen -> Signed.long -> gen
val bnrmap : gen -> gen -> gen
val bnrsurjection : gen -> gen -> gen
val bnfnarrow : gen -> gen
val bnfcertify : gen -> Signed.long
val bnfcertify0 : gen -> Signed.long -> Signed.long
val bnrcompositum : gen -> gen -> gen
val decodemodule : gen -> gen -> gen
val discrayabslist : gen -> gen -> gen
val discrayabslistarch : gen -> gen -> pari_ulong -> gen
val idealmoddivisor : gen -> gen -> gen
val isprincipalray : gen -> gen -> gen
val isprincipalraygen : gen -> gen -> gen
val nf_deg1_prime : gen -> gen
val nfarchstar : gen -> gen -> gen -> gen
val rnfconductor : gen -> gen -> gen
val rnfconductor0 : gen -> gen -> Signed.long -> gen
val rnfnormgroup : gen -> gen -> gen
val subgrouplist0 : gen -> gen -> Signed.long -> gen
val bnfisnorm : gen -> gen -> Signed.long -> gen
val rnfisnorm : gen -> gen -> Signed.long -> gen
val rnfisnorminit : gen -> gen -> int -> gen
val coprimes_zv : pari_ulong -> gen
val char_check : gen -> gen -> int
val charconj : gen -> gen -> gen
val charconj0 : gen -> gen -> gen
val chardiv : gen -> gen -> gen -> gen
val chardiv0 : gen -> gen -> gen -> gen
val chareval : gen -> gen -> gen -> gen -> gen
val chargalois : gen -> gen -> gen
val charker : gen -> gen -> gen
val charker0 : gen -> gen -> gen
val charmul : gen -> gen -> gen -> gen
val charmul0 : gen -> gen -> gen -> gen
val charorder : gen -> gen -> gen
val charorder0 : gen -> gen -> gen
val charpow : gen -> gen -> gen -> gen
val charpow0 : gen -> gen -> gen -> gen
val char_denormalize : gen -> gen -> gen -> gen
val char_normalize : gen -> gen -> gen
val char_simplify : gen -> gen -> gen
val checkznstar_i : gen -> int
val cyc_normalize : gen -> gen
val ncharvecexpo : gen -> gen -> gen
val znchar : gen -> gen
val znchar_quad : gen -> gen -> gen
val zncharcheck : gen -> gen -> int
val zncharconductor : gen -> gen -> gen
val zncharconj : gen -> gen -> gen
val znchardecompose : gen -> gen -> gen -> gen
val znchardiv : gen -> gen -> gen -> gen
val znchareval : gen -> gen -> gen -> gen -> gen
val zncharinduce : gen -> gen -> gen -> gen
val zncharisodd : gen -> gen -> Signed.long
val zncharker : gen -> gen -> gen
val zncharmul : gen -> gen -> gen -> gen
val zncharorder : gen -> gen -> gen
val zncharpow : gen -> gen -> gen -> gen
val znchartokronecker : gen -> gen -> Signed.long -> gen
val znchartoprimitive : gen -> gen -> gen
val znconrey_check : gen -> gen -> int
val znconrey_normalized : gen -> gen -> gen
val znconreychar : gen -> gen -> gen
val znconreyfromchar_normalized : gen -> gen -> gen
val znconreyconductor : gen -> gen -> gen Ctypes_static.ptr -> gen
val znconreyexp : gen -> gen -> gen
val znconreyfromchar : gen -> gen -> gen
val znconreylog : gen -> gen -> gen
val znconreylog_normalize : gen -> gen -> gen
val znlog0 : gen -> gen -> gen -> gen
val zv_cyc_minimal : gen -> gen -> gen -> Signed.long
val zv_cyc_minimize : gen -> gen -> gen -> Signed.long
val closure_deriv : gen -> gen
val closure_derivn : gen -> Signed.long -> gen

val localvars_find :
  gen -> entree Ctypes.structure Ctypes_static.ptr -> Signed.long

val localvars_read_str : string -> gen -> gen
val snm_closure : entree Ctypes.structure Ctypes_static.ptr -> gen -> gen
val strtoclosure : string -> Signed.long -> gen
val strtofunction : string -> gen
val gconcat : gen -> gen -> gen
val gconcat1 : gen -> gen
val matconcat : gen -> gen
val shallowconcat : gen -> gen -> gen
val shallowconcat1 : gen -> gen
val shallowmatconcat : gen -> gen
val vconcat : gen -> gen -> gen
val default0 : string -> string -> gen
val getrealprecision : unit -> Signed.long
val pari_is_default : string -> entree Ctypes.structure Ctypes_static.ptr
val sd_texstyle : string -> Signed.long -> gen
val sd_colors : string -> Signed.long -> gen
val sd_compatible : string -> Signed.long -> gen
val sd_datadir : string -> Signed.long -> gen
val sd_debug : string -> Signed.long -> gen
val sd_debugfiles : string -> Signed.long -> gen
val sd_debugmem : string -> Signed.long -> gen
val sd_factor_add_primes : string -> Signed.long -> gen
val sd_factor_proven : string -> Signed.long -> gen
val sd_format : string -> Signed.long -> gen
val sd_histsize : string -> Signed.long -> gen
val sd_log : string -> Signed.long -> gen
val sd_logfile : string -> Signed.long -> gen
val sd_nbthreads : string -> Signed.long -> gen
val sd_new_galois_format : string -> Signed.long -> gen
val sd_output : string -> Signed.long -> gen
val sd_parisize : string -> Signed.long -> gen
val sd_parisizemax : string -> Signed.long -> gen
val sd_path : string -> Signed.long -> gen
val sd_plothsizes : string -> Signed.long -> gen
val sd_prettyprinter : string -> Signed.long -> gen
val sd_primelimit : string -> Signed.long -> gen
val sd_realbitprecision : string -> Signed.long -> gen
val sd_realprecision : string -> Signed.long -> gen
val sd_secure : string -> Signed.long -> gen
val sd_seriesprecision : string -> Signed.long -> gen
val sd_simplify : string -> Signed.long -> gen
val sd_sopath : string -> int -> gen
val sd_strictargs : string -> Signed.long -> gen
val sd_strictmatch : string -> Signed.long -> gen

val sd_string :
  string -> Signed.long -> string -> string Ctypes_static.ptr -> gen

val sd_threadsize : string -> Signed.long -> gen
val sd_threadsizemax : string -> Signed.long -> gen

val sd_intarray :
  string -> Signed.long -> gen Ctypes_static.ptr -> string -> gen

val sd_toggle : string -> Signed.long -> string -> int Ctypes_static.ptr -> gen

val sd_ulong :
  string ->
  Signed.long ->
  string ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  pari_ulong ->
  string Ctypes_static.ptr ->
  gen

val setdefault : string -> string -> Signed.long -> gen

val setrealprecision :
  Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val digits : gen -> gen -> gen
val fromdigits : gen -> gen -> gen
val fromdigitsu : gen -> gen -> gen

val gen_digits :
  gen ->
  gen ->
  Signed.long ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen Ctypes_static.ptr -> gen)
  Ctypes_static.static_funptr ->
  gen

val gen_fromdigits :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  bb_ring Ctypes.structure Ctypes_static.ptr ->
  gen

val sumdigits : gen -> gen
val sumdigits0 : gen -> gen -> gen
val sumdigitsu : pari_ulong -> pari_ulong
val ecpp : gen -> gen
val ecpp0 : gen -> Signed.long -> gen
val ecppexport : gen -> Signed.long -> gen
val ecppisvalid : gen -> Signed.long
val isprimeecpp : gen -> Signed.long
val sd_breakloop : string -> Signed.long -> gen
val sd_echo : string -> Signed.long -> gen
val sd_graphcolormap : string -> Signed.long -> gen
val sd_graphcolors : string -> Signed.long -> gen
val sd_help : string -> Signed.long -> gen
val sd_histfile : string -> Signed.long -> gen
val sd_lines : string -> Signed.long -> gen
val sd_linewrap : string -> Signed.long -> gen
val sd_prompt : string -> Signed.long -> gen
val sd_prompt_cont : string -> Signed.long -> gen
val sd_psfile : string -> Signed.long -> gen
val sd_readline : string -> Signed.long -> gen
val sd_recover : string -> Signed.long -> gen
val sd_timer : string -> Signed.long -> gen
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
val gp_alarm : Signed.long -> gen -> gen
val gp_input : unit -> gen
val gp_allocatemem : gen -> unit
val gp_handle_exception : Signed.long -> int
val gp_alarm_handler : int -> unit
val gp_sigint_fun : unit -> unit
val gp_help : string -> Signed.long -> unit
val gp_echo_and_log : string -> string -> unit
val print_fun_list : string Ctypes_static.ptr -> Signed.long -> unit
val strtime : Signed.long -> gen

val direuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  gen ->
  gen

val dirpowers : Signed.long -> gen -> Signed.long -> gen
val dirpowerssum : pari_ulong -> gen -> Signed.long -> Signed.long -> gen

val dirpowerssumfun :
  pari_ulong ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> pari_ulong -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  gen

val vecpowuu : Signed.long -> pari_ulong -> gen
val vecpowug : Signed.long -> gen -> Signed.long -> gen
val ellanalyticrank : gen -> gen -> Signed.long -> gen
val ellanalyticrank_bitprec : gen -> gen -> Signed.long -> gen

val ellanal_globalred_all :
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val ellheegner : gen -> gen
val elll1 : gen -> Signed.long -> Signed.long -> gen
val elll1_bitprec : gen -> Signed.long -> Signed.long -> gen
val ellconvertname : gen -> gen
val elldatagenerators : gen -> gen
val ellidentify : gen -> gen
val ellsearch : gen -> gen
val ellsearchcurve : gen -> gen

val forell :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  unit

val ellfromeqn : gen -> gen
val akell : gen -> gen -> gen
val bilhell : gen -> gen -> gen -> Signed.long -> gen
val checkell : gen -> unit
val checkell_fq : gen -> unit
val checkell_q : gen -> unit
val checkell_qp : gen -> unit
val checkellisog : gen -> unit
val checkellpt : gen -> unit
val checkell5 : gen -> unit
val cxredsl2 : gen -> gen Ctypes_static.ptr -> gen
val cxredsl2_i : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen
val ec_2divpol_evalx : gen -> gen -> gen
val ec_3divpol_evalx : gen -> gen -> gen
val ec_bmodel : gen -> Signed.long -> gen
val ec_f_evalx : gen -> gen -> gen
val ec_h_evalx : gen -> gen -> gen
val ec_dfdx_evalq : gen -> gen -> gen
val ec_dfdy_evalq : gen -> gen -> gen
val ec_dmfdy_evalq : gen -> gen -> gen
val ec_half_deriv_2divpol : gen -> Signed.long -> gen
val ec_half_deriv_2divpol_evalx : gen -> gen -> gen
val ec_phi2 : gen -> Signed.long -> gen
val ell_is_integral : gen -> int
val ellq_get_cm : gen -> Signed.long
val ellq_get_n : gen -> gen
val ellq_get_nfa : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit
val ellqp_tate_uniformization : gen -> Signed.long -> gen
val ellqp_agm : gen -> Signed.long -> gen
val ellqp_u : gen -> Signed.long -> gen
val ellqp_u2 : gen -> Signed.long -> gen
val ellqp_q : gen -> Signed.long -> gen
val ellqp_ab : gen -> Signed.long -> gen
val ellqp_l : gen -> Signed.long -> gen
val ellqp_root : gen -> Signed.long -> gen
val ellqtwist_bsdperiod : gen -> Signed.long -> gen
val ellr_area : gen -> Signed.long -> gen
val ellr_ab : gen -> Signed.long -> gen
val ellr_eta : gen -> Signed.long -> gen
val ellr_omega : gen -> Signed.long -> gen
val ellr_roots : gen -> Signed.long -> gen
val elladd : gen -> gen -> gen -> gen
val ellan : gen -> Signed.long -> gen
val ellanq_zv : gen -> Signed.long -> gen
val ellanal_globalred : gen -> gen Ctypes_static.ptr -> gen
val ellap : gen -> gen -> gen
val ellap_cm_fast : gen -> pari_ulong -> Signed.long -> Signed.long
val ellbasechar : gen -> gen
val ellbsd : gen -> Signed.long -> gen
val ellcard : gen -> gen -> gen
val ellchangecurve : gen -> gen -> gen
val ellchangeinvert : gen -> gen
val ellchangepoint : gen -> gen -> gen
val ellchangepointinv : gen -> gen -> gen
val elldivpol : gen -> Signed.long -> Signed.long -> gen
val elleisnum : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val elleta : gen -> Signed.long -> gen
val elleulerf : gen -> gen -> gen
val ellff_get_card : gen -> gen
val ellff_get_gens : gen -> gen
val ellff_get_group : gen -> gen
val ellff_get_o : gen -> gen
val ellff_get_p : gen -> gen
val ellff_get_m : gen -> gen
val ellff_get_d : gen -> gen
val ellfromj : gen -> gen
val ellgenerators : gen -> gen
val ellglobalred : gen -> gen
val ellgroup : gen -> gen -> gen
val ellgroup0 : gen -> gen -> Signed.long -> gen
val ellheight0 : gen -> gen -> gen -> Signed.long -> gen
val ellheight : gen -> gen -> Signed.long -> gen
val ellheightmatrix : gen -> gen -> Signed.long -> gen
val ellheightoo : gen -> gen -> Signed.long -> gen
val ellinit : gen -> gen -> Signed.long -> gen
val ellintegralmodel : gen -> gen Ctypes_static.ptr -> gen
val ellintegralmodel_i : gen -> gen Ctypes_static.ptr -> gen
val elliscm : gen -> Signed.long
val ellisoncurve : gen -> gen -> gen
val ellisotree : gen -> gen
val ellissupersingular : gen -> gen -> int
val elljissupersingular : gen -> int
val elllseries : gen -> gen -> gen -> Signed.long -> gen
val elllocalred : gen -> gen -> gen
val elllog : gen -> gen -> gen -> gen -> gen
val ellminimaldisc : gen -> gen
val ellminimalmodel : gen -> gen Ctypes_static.ptr -> gen
val ellminimaltwist : gen -> gen
val ellminimaltwist0 : gen -> Signed.long -> gen
val ellminimaltwistcond : gen -> gen
val ellmul : gen -> gen -> gen -> gen
val ellnf_vecarea : gen -> Signed.long -> gen
val ellnf_veceta : gen -> Signed.long -> gen
val ellnf_vecomega : gen -> Signed.long -> gen
val ellneg : gen -> gen -> gen
val ellorder : gen -> gen -> gen -> gen
val ellorder_q : gen -> gen -> Signed.long
val ellordinate : gen -> gen -> Signed.long -> gen
val ellpadicheight0 : gen -> gen -> Signed.long -> gen -> gen -> gen
val ellpadicheightmatrix : gen -> gen -> Signed.long -> gen -> gen
val ellperiods : gen -> Signed.long -> Signed.long -> gen
val ellrandom : gen -> gen
val ellrootno : gen -> gen -> Signed.long
val ellrootno_global : gen -> Signed.long
val ellsaturation : gen -> gen -> Signed.long -> Signed.long -> gen
val ellsea : gen -> Signed.long -> gen
val ellsigma : gen -> gen -> Signed.long -> Signed.long -> gen
val ellsub : gen -> gen -> gen -> gen
val ellsupersingularj : gen -> gen
val elltamagawa : gen -> gen
val elltaniyama : gen -> Signed.long -> gen
val elltatepairing : gen -> gen -> gen -> gen -> gen
val elltors : gen -> gen
val elltors0 : gen -> Signed.long -> gen
val elltors_psylow : gen -> pari_ulong -> gen
val elltrace : gen -> gen -> gen
val elltwist : gen -> gen -> gen
val ellweilpairing : gen -> gen -> gen -> gen -> gen
val ellwp : gen -> gen -> Signed.long -> gen
val ellwp0 : gen -> gen -> Signed.long -> Signed.long -> gen
val ellwpseries : gen -> Signed.long -> Signed.long -> gen
val ellxn : gen -> Signed.long -> Signed.long -> gen
val ellzeta : gen -> gen -> Signed.long -> gen
val oncurve : gen -> gen -> int
val orderell : gen -> gen -> gen
val pointell : gen -> gen -> Signed.long -> gen
val point_to_a4a6 : gen -> gen -> gen -> gen Ctypes_static.ptr -> gen

val point_to_a4a6_fl :
  gen -> gen -> pari_ulong -> pari_ulong Ctypes_static.ptr -> gen

val zell : gen -> gen -> Signed.long -> gen
val qp_agm2_sequence : gen -> gen -> gen

val qp_ascending_landen :
  gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit

val qp_descending_landen :
  gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit

val ellformaldifferential : gen -> Signed.long -> Signed.long -> gen
val ellformalexp : gen -> Signed.long -> Signed.long -> gen
val ellformallog : gen -> Signed.long -> Signed.long -> gen
val ellformalpoint : gen -> Signed.long -> Signed.long -> gen
val ellformalw : gen -> Signed.long -> Signed.long -> gen
val ellnonsingularmultiple : gen -> gen -> gen
val ellpadicl : gen -> gen -> Signed.long -> gen -> Signed.long -> gen -> gen
val ellpadicbsd : gen -> gen -> Signed.long -> gen -> gen
val ellpadicfrobenius : gen -> pari_ulong -> Signed.long -> gen
val ellpadicheight : gen -> gen -> Signed.long -> gen -> gen
val ellpadiclog : gen -> gen -> Signed.long -> gen -> gen
val ellpadicregulator : gen -> gen -> Signed.long -> gen -> gen
val ellpadics2 : gen -> gen -> Signed.long -> gen
val ell2cover : gen -> Signed.long -> gen
val ellrank : gen -> Signed.long -> gen -> Signed.long -> gen
val ellrankinit : gen -> Signed.long -> gen
val hyperell_locally_soluble : gen -> gen -> Signed.long
val nf_hyperell_locally_soluble : gen -> gen -> gen -> Signed.long
val nfhilbert : gen -> gen -> gen -> Signed.long
val nfhilbert0 : gen -> gen -> gen -> gen -> Signed.long
val ellisdivisible : gen -> gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val ellisogenyapply : gen -> gen -> gen
val ellisogeny : gen -> gen -> Signed.long -> Signed.long -> Signed.long -> gen
val ellisomat : gen -> Signed.long -> Signed.long -> gen
val ellweilcurve : gen -> gen Ctypes_static.ptr -> gen

val flxq_elldivpolmod :
  gen -> gen -> Signed.long -> gen -> gen -> pari_ulong -> gen

val fp_ellcard_sea : gen -> gen -> gen -> Signed.long -> gen
val fq_ellcard_sea : gen -> gen -> gen -> gen -> gen -> Signed.long -> gen
val fq_elldivpolmod : gen -> gen -> Signed.long -> gen -> gen -> gen -> gen
val ellmodulareqn : Signed.long -> Signed.long -> Signed.long -> gen
val externstr : string -> gen
val gp_filter : string -> string
val gpextern : string -> gen
val gpsystem : string -> Signed.long
val readstr : string -> gen
val gentogenstr_nospace : gen -> gen
val gentogenstr : gen -> gen
val gentotexstr : gen -> string
val gentostr : gen -> string
val gentostr_raw : gen -> string
val gentostr_unquoted : gen -> string
val str : gen -> gen
val strexpand : gen -> gen
val strtex : gen -> gen
val brute : gen -> char -> Signed.long -> unit
val dbggen : gen -> Signed.long -> unit
val error0 : gen -> unit
val dbg_pari_heap : unit -> unit
val err_flush : unit -> unit
val err_printf : string -> unit
val gp_getenv : string -> gen
val gp_fileclose : Signed.long -> unit
val gp_fileextern : string -> Signed.long
val gp_fileflush : Signed.long -> unit
val gp_fileflush0 : gen -> unit
val gp_fileopen : string -> string -> Signed.long
val gp_fileread : Signed.long -> gen
val gp_filereadstr : Signed.long -> gen
val gp_filewrite : Signed.long -> string -> unit
val gp_filewrite1 : Signed.long -> string -> unit
val gp_read_file : string -> gen
val gp_read_str_multiline : string -> string -> gen
val gp_readvec_file : string -> gen
val gpinstall : string -> string -> string -> string -> unit
val gsprintf : string -> gen
val itostr : gen -> string
val matbrute : gen -> char -> Signed.long -> unit
val os_getenv : string -> string
val uordinal : pari_ulong -> string
val outmat : gen -> unit
val output : gen -> unit
val rgv_to_str : gen -> Signed.long -> string
val pari_add_hist : gen -> Signed.long -> Signed.long -> unit
val pari_ask_confirm : string -> unit
val pari_flush : unit -> unit
val pari_get_hist : Signed.long -> gen
val pari_get_histrtime : Signed.long -> Signed.long
val pari_get_histtime : Signed.long -> Signed.long
val pari_get_homedir : string -> string
val pari_histtime : Signed.long -> gen
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
val pari_sprint0 : string -> gen -> Signed.long -> string
val print : gen -> unit
val printp : gen -> unit
val print1 : gen -> unit
val printf0 : string -> gen -> unit
val printsep : string -> gen -> unit
val printsep1 : string -> gen -> unit
val printtex : gen -> unit
val stack_sprintf : string -> string
val str_init : pari_str Ctypes.structure Ctypes_static.ptr -> int -> unit
val str_printf : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val str_putc : pari_str Ctypes.structure Ctypes_static.ptr -> char -> unit
val str_puts : pari_str Ctypes.structure Ctypes_static.ptr -> string -> unit
val strftime_expand : string -> string -> Signed.long -> unit
val strprintf : string -> gen -> gen
val term_color : Signed.long -> unit
val term_get_color : string -> Signed.long -> string
val texe : gen -> char -> Signed.long -> unit
val warning0 : gen -> unit
val write0 : string -> gen -> unit
val write1 : string -> gen -> unit
val writebin : string -> gen -> unit
val writetex : string -> gen -> unit
val bincopy_relink : gen -> gen -> unit
val bitprecision0 : gen -> Signed.long -> gen
val bitprecision00 : gen -> gen -> gen
val break0 : Signed.long -> gen
val call0 : gen -> gen -> gen
val closure_callgen0prec : gen -> Signed.long -> gen
val closure_callgen1 : gen -> gen -> gen
val closure_callgen1prec : gen -> gen -> Signed.long -> gen
val closure_callgen2 : gen -> gen -> gen -> gen
val closure_callgenall : gen -> Signed.long -> gen
val closure_callgenvec : gen -> gen -> gen
val closure_callgenvecdef : gen -> gen -> gen -> gen
val closure_callgenvecdefprec : gen -> gen -> gen -> Signed.long -> gen
val closure_callgenvecprec : gen -> gen -> Signed.long -> gen
val closure_callvoid1 : gen -> gen -> unit
val closure_context : Signed.long -> Signed.long -> Signed.long
val closure_disassemble : gen -> unit
val closure_err : Signed.long -> unit
val closure_evalbrk : gen -> Signed.long Ctypes_static.ptr -> gen
val closure_evalgen : gen -> gen
val closure_evalnobrk : gen -> gen
val closure_evalres : gen -> gen
val closure_evalvoid : gen -> unit
val closure_func_err : unit -> string
val closure_trapgen : gen -> Signed.long -> gen
val copybin_unlink : gen -> gen
val getlocalprec : Signed.long -> Signed.long
val getlocalbitprec : Signed.long -> Signed.long
val get_lex : Signed.long -> gen
val get_localprec : unit -> Signed.long
val get_localbitprec : unit -> Signed.long
val gp_call : unit Ctypes_static.ptr -> gen -> gen
val gp_callprec : unit Ctypes_static.ptr -> gen -> Signed.long -> gen
val gp_call2 : unit Ctypes_static.ptr -> gen -> gen -> gen
val gp_callbool : unit Ctypes_static.ptr -> gen -> Signed.long
val gp_callvoid : unit Ctypes_static.ptr -> gen -> Signed.long
val gp_eval : unit Ctypes_static.ptr -> gen -> gen
val gp_evalbool : unit Ctypes_static.ptr -> gen -> Signed.long
val gp_evalprec : unit Ctypes_static.ptr -> gen -> Signed.long -> gen
val gp_evalupto : unit Ctypes_static.ptr -> gen -> gen
val gp_evalvoid : unit Ctypes_static.ptr -> gen -> Signed.long
val localprec : gen -> unit
val localbitprec : gen -> unit
val loop_break : unit -> Signed.long
val next0 : Signed.long -> gen
val pareval : gen -> gen
val pari_self : unit -> gen
val parsum : gen -> gen -> gen -> gen
val parvector : Signed.long -> gen -> gen
val pop_lex : Signed.long -> unit
val pop_localprec : unit -> unit
val precision0 : gen -> Signed.long -> gen
val precision00 : gen -> gen -> gen
val push_lex : gen -> gen -> unit
val push_localbitprec : Signed.long -> unit
val push_localprec : Signed.long -> unit
val return0 : gen -> gen
val set_lex : Signed.long -> gen -> unit

val forcomposite_init :
  forcomposite_t Ctypes.structure Ctypes_static.ptr -> gen -> gen -> int

val forcomposite_next : forcomposite_t Ctypes.structure Ctypes_static.ptr -> gen
val forprime_next : forprime_t Ctypes.structure Ctypes_static.ptr -> gen

val forprime_init :
  forprime_t Ctypes.structure Ctypes_static.ptr -> gen -> gen -> int

val forprimestep_init :
  forprime_t Ctypes.structure Ctypes_static.ptr -> gen -> gen -> gen -> int

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
val prodprimes : unit -> gen

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

val ff_1 : gen -> gen
val ff_frobenius : gen -> Signed.long -> gen
val ff_z_z_muldiv : gen -> gen -> gen -> gen
val ff_q_add : gen -> gen -> gen
val ff_z_add : gen -> gen -> gen
val ff_z_mul : gen -> gen -> gen
val ff_add : gen -> gen -> gen
val ff_charpoly : gen -> gen
val ff_conjvec : gen -> gen
val ff_div : gen -> gen -> gen
val ff_ellcard : gen -> gen
val ff_ellcard_sea : gen -> Signed.long -> gen
val ff_ellgens : gen -> gen
val ff_ellgroup : gen -> gen Ctypes_static.ptr -> gen
val ff_elllog : gen -> gen -> gen -> gen -> gen
val ff_ellmul : gen -> gen -> gen -> gen
val ff_ellorder : gen -> gen -> gen -> gen
val ff_elltwist : gen -> gen
val ff_ellrandom : gen -> gen
val ff_elltatepairing : gen -> gen -> gen -> gen -> gen
val ff_ellweilpairing : gen -> gen -> gen -> gen -> gen
val ff_equal : gen -> gen -> int
val ff_equal0 : gen -> int
val ff_equal1 : gen -> int
val ff_equalm1 : gen -> int
val ff_f : gen -> Signed.long
val ff_gen : gen -> gen
val ff_inv : gen -> gen
val ff_issquare : gen -> Signed.long
val ff_issquareall : gen -> gen Ctypes_static.ptr -> Signed.long
val ff_ispower : gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val ff_log : gen -> gen -> gen -> gen
val ff_map : gen -> gen -> gen
val ff_minpoly : gen -> gen
val ff_mod : gen -> gen
val ff_mul : gen -> gen -> gen
val ff_mul2n : gen -> Signed.long -> gen
val ff_neg : gen -> gen
val ff_neg_i : gen -> gen
val ff_norm : gen -> gen
val ff_order : gen -> gen -> gen
val ff_p : gen -> gen
val ff_p_i : gen -> gen
val ff_pow : gen -> gen -> gen
val ff_primroot : gen -> gen Ctypes_static.ptr -> gen
val ff_q : gen -> gen
val ff_samefield : gen -> gen -> int
val ff_sqr : gen -> gen
val ff_sqrt : gen -> gen
val ff_sqrtn : gen -> gen -> gen Ctypes_static.ptr -> gen
val ff_sub : gen -> gen -> gen
val ff_to_f2xq : gen -> gen
val ff_to_f2xq_i : gen -> gen
val ff_to_flxq : gen -> gen
val ff_to_flxq_i : gen -> gen
val ff_to_fpxq : gen -> gen
val ff_to_fpxq_i : gen -> gen
val ff_trace : gen -> gen
val ff_var : gen -> Signed.long
val ff_zero : gen -> gen
val ffm_ffc_invimage : gen -> gen -> gen -> gen
val ffm_ffc_gauss : gen -> gen -> gen -> gen
val ffm_ffc_mul : gen -> gen -> gen -> gen
val ffm_deplin : gen -> gen -> gen
val ffm_det : gen -> gen -> gen
val ffm_gauss : gen -> gen -> gen -> gen
val ffm_image : gen -> gen -> gen
val ffm_indexrank : gen -> gen -> gen
val ffm_inv : gen -> gen -> gen
val ffm_invimage : gen -> gen -> gen -> gen
val ffm_ker : gen -> gen -> gen
val ffm_mul : gen -> gen -> gen -> gen
val ffm_rank : gen -> gen -> Signed.long
val ffm_suppl : gen -> gen -> gen
val ffx_add : gen -> gen -> gen -> gen
val ffx_ddf : gen -> gen -> gen
val ffx_degfact : gen -> gen -> gen
val ffx_disc : gen -> gen -> gen

val ffx_extgcd :
  gen -> gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val ffx_factor : gen -> gen -> gen
val ffx_factor_squarefree : gen -> gen -> gen
val ffx_gcd : gen -> gen -> gen -> gen
val ffx_halfgcd : gen -> gen -> gen -> gen

val ffx_halfgcd_all :
  gen -> gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val ffx_ispower :
  gen -> Signed.long -> gen -> gen Ctypes_static.ptr -> Signed.long

val ffx_mul : gen -> gen -> gen -> gen
val ffx_preimage : gen -> gen -> gen -> gen
val ffx_preimagerel : gen -> gen -> gen -> gen
val ffx_rem : gen -> gen -> gen -> gen
val ffx_resultant : gen -> gen -> gen -> gen
val ffx_roots : gen -> gen -> gen
val ffx_sqr : gen -> gen -> gen
val ffxq_inv : gen -> gen -> gen -> gen
val ffxq_minpoly : gen -> gen -> gen -> gen
val ffxq_mul : gen -> gen -> gen -> gen -> gen
val ffxq_sqr : gen -> gen -> gen -> gen
val fqx_to_ffx : gen -> gen -> gen
val fq_to_ff : gen -> gen -> gen
val z_ff_div : gen -> gen -> gen
val ffembed : gen -> gen -> gen
val ffextend : gen -> gen -> Signed.long -> gen
val fffrobenius : gen -> Signed.long -> gen
val ffgen : gen -> Signed.long -> gen
val ffinvmap : gen -> gen
val fflog : gen -> gen -> gen -> gen
val ffmap : gen -> gen -> gen
val ffmaprel : gen -> gen -> gen
val ffcompomap : gen -> gen -> gen
val fforder : gen -> gen -> gen
val ffprimroot : gen -> gen Ctypes_static.ptr -> gen
val ffrandom : gen -> gen
val rg_is_ff : gen -> gen Ctypes_static.ptr -> int
val rgc_is_ffc : gen -> gen Ctypes_static.ptr -> int
val rgm_is_ffm : gen -> gen Ctypes_static.ptr -> int
val p_to_ff : gen -> Signed.long -> gen
val tp_to_ff : gen -> gen -> gen
val flx_factcyclo : pari_ulong -> pari_ulong -> pari_ulong -> gen
val fpx_factcyclo : pari_ulong -> gen -> pari_ulong -> gen
val factormodcyclo : Signed.long -> gen -> Signed.long -> Signed.long -> gen
val checkgal : gen -> gen
val checkgroup : gen -> gen Ctypes_static.ptr -> gen
val checkgroupelts : gen -> gen
val embed_disc : gen -> Signed.long -> Signed.long -> gen
val embed_roots : gen -> Signed.long -> gen
val galois_group : gen -> gen
val galoisconj : gen -> gen -> gen
val galoisconj0 : gen -> Signed.long -> gen -> Signed.long -> gen
val galoisconjclasses : gen -> gen
val galoisexport : gen -> Signed.long -> gen
val galoisfixedfield : gen -> gen -> Signed.long -> Signed.long -> gen
val galoisidentify : gen -> gen
val galoisinit : gen -> gen -> gen
val galoisisabelian : gen -> Signed.long -> gen
val galoisisnormal : gen -> gen -> Signed.long
val galoispermtopol : gen -> gen -> gen
val galoissplittinginit : gen -> gen -> gen
val galoissubgroups : gen -> gen
val galoissubfields : gen -> Signed.long -> Signed.long -> gen
val numberofconjugates : gen -> Signed.long -> Signed.long
val polgalois : gen -> Signed.long -> gen
val galoisnbpol : Signed.long -> gen
val galoisgetgroup : Signed.long -> Signed.long -> gen
val galoisgetname : Signed.long -> Signed.long -> gen
val galoisgetpol : Signed.long -> Signed.long -> Signed.long -> gen
val conj_i : gen -> gen
val conjvec : gen -> Signed.long -> gen
val divrunextu : gen -> pari_ulong -> gen
val gadd : gen -> gen -> gen
val gaddsg : Signed.long -> gen -> gen
val gconj : gen -> gen
val gdiv : gen -> gen -> gen
val gdivgs : gen -> Signed.long -> gen
val gdivgu : gen -> pari_ulong -> gen
val gdivgunextu : gen -> pari_ulong -> gen
val ginv : gen -> gen
val gmul : gen -> gen -> gen
val gmul2n : gen -> Signed.long -> gen
val gmulsg : Signed.long -> gen -> gen
val gmulug : pari_ulong -> gen -> gen
val gsqr : gen -> gen
val gsub : gen -> gen -> gen
val gsubsg : Signed.long -> gen -> gen
val mulcxi : gen -> gen
val mulcxmi : gen -> gen
val mulcxpowis : gen -> Signed.long -> gen
val qdivii : gen -> gen -> gen
val qdiviu : gen -> pari_ulong -> gen
val qdivis : gen -> Signed.long -> gen
val ser_normalize : gen -> gen

val gassoc_proto :
  (gen -> gen -> gen) Ctypes_static.static_funptr -> gen -> gen -> gen

val map_proto_g : (gen -> gen) Ctypes_static.static_funptr -> gen -> gen

val map_proto_lg :
  (gen -> Signed.long) Ctypes_static.static_funptr -> gen -> gen

val map_proto_lgl :
  (gen -> Signed.long -> Signed.long) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val q_lval : gen -> pari_ulong -> Signed.long
val q_lvalrem : gen -> pari_ulong -> gen Ctypes_static.ptr -> Signed.long
val q_pval : gen -> gen -> Signed.long
val q_pvalrem : gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val rgx_val : gen -> Signed.long
val rgx_valrem : gen -> gen Ctypes_static.ptr -> Signed.long
val rgx_valrem_inexact : gen -> gen Ctypes_static.ptr -> Signed.long
val rgxv_maxdegree : gen -> Signed.long
val zv_z_dvd : gen -> gen -> int
val zv_pval : gen -> gen -> Signed.long
val zv_pvalrem : gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val zv_lval : gen -> pari_ulong -> Signed.long
val zv_lvalrem : gen -> pari_ulong -> gen Ctypes_static.ptr -> Signed.long
val zx_lvalrem : gen -> pari_ulong -> gen Ctypes_static.ptr -> Signed.long
val zx_pval : gen -> gen -> Signed.long
val zx_pvalrem : gen -> gen -> gen Ctypes_static.ptr -> Signed.long

val z_lvalrem_stop :
  gen Ctypes_static.ptr -> pari_ulong -> int Ctypes_static.ptr -> Signed.long

val cgetp : gen -> gen
val cvstop2 : Signed.long -> gen -> gen
val cvtop : gen -> gen -> Signed.long -> gen
val cvtop2 : gen -> gen -> gen
val cx_approx_equal : gen -> gen -> int
val cx_approx0 : gen -> gen -> int
val gabs : gen -> Signed.long -> gen
val gaffect : gen -> gen -> unit
val gaffsg : Signed.long -> gen -> unit
val gcmp : gen -> gen -> int
val gequal0 : gen -> int
val gequal1 : gen -> int
val gequalx : gen -> int
val gequalm1 : gen -> int
val gcmpsg : Signed.long -> gen -> int
val gcvtop : gen -> gen -> Signed.long -> gen
val gequal : gen -> gen -> int
val gequalsg : Signed.long -> gen -> int
val gexpo : gen -> Signed.long
val gexpo_safe : gen -> Signed.long
val gpexponent : gen -> gen
val gpvaluation : gen -> gen -> gen
val gvaluation : gen -> gen -> Signed.long
val gidentical : gen -> gen -> int
val glength : gen -> Signed.long
val gmax : gen -> gen -> gen
val gmaxgs : gen -> Signed.long -> gen
val gmin : gen -> gen -> gen
val gmings : gen -> Signed.long -> gen
val gneg : gen -> gen
val gneg_i : gen -> gen
val gsigne : gen -> int
val gtolist : gen -> gen
val gtolong : gen -> Signed.long
val lexcmp : gen -> gen -> int
val listinsert : gen -> gen -> Signed.long -> gen
val listpop : gen -> Signed.long -> unit
val listpop0 : gen -> Signed.long -> unit
val listput : gen -> gen -> Signed.long -> gen
val listput0 : gen -> gen -> Signed.long -> unit
val listsort : gen -> Signed.long -> unit
val matsize : gen -> gen
val mklist : unit -> gen
val mklist_typ : Signed.long -> gen
val mklistcopy : gen -> gen
val mkmap : unit -> gen
val normalizeser : gen -> gen
val normalizepol : gen -> gen
val normalizepol_approx : gen -> Signed.long -> gen
val normalizepol_lg : gen -> Signed.long -> gen
val padic_to_fl : gen -> pari_ulong -> pari_ulong
val padic_to_fp : gen -> gen -> gen
val quadtofp : gen -> Signed.long -> gen
val sizedigit : gen -> Signed.long
val u_lval : pari_ulong -> pari_ulong -> Signed.long

val u_lvalrem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val u_lvalrem_stop :
  pari_ulong Ctypes_static.ptr ->
  pari_ulong ->
  int Ctypes_static.ptr ->
  Signed.long

val u_pval : pari_ulong -> gen -> Signed.long
val u_pvalrem : pari_ulong -> gen -> pari_ulong Ctypes_static.ptr -> Signed.long
val vecindexmax : gen -> Signed.long
val vecindexmin : gen -> Signed.long
val vecmax0 : gen -> gen Ctypes_static.ptr -> gen
val vecmax : gen -> gen
val vecmin0 : gen -> gen Ctypes_static.ptr -> gen
val vecmin : gen -> gen
val z_lval : Signed.long -> pari_ulong -> Signed.long

val z_lvalrem :
  Signed.long -> pari_ulong -> Signed.long Ctypes_static.ptr -> Signed.long

val z_pval : Signed.long -> gen -> Signed.long

val z_pvalrem :
  Signed.long -> gen -> Signed.long Ctypes_static.ptr -> Signed.long

val zx_lval : gen -> Signed.long -> Signed.long
val hgmcyclo : gen -> gen
val hgmalpha : gen -> gen
val hgmgamma : gen -> gen
val hgminit : gen -> gen -> gen
val hgmparams : gen -> gen
val hgmeulerfactor : gen -> gen -> Signed.long -> gen Ctypes_static.ptr -> gen
val hgmcoef : gen -> gen -> gen -> gen
val hgmcoefs : gen -> gen -> Signed.long -> gen
val hgmtwist : gen -> gen
val hgmissymmetrical : gen -> Signed.long
val hgmbydegree : Signed.long -> gen
val lfunhgm : gen -> gen -> gen -> Signed.long -> gen
val qp_zeta : gen -> gen
val lerchphi : gen -> gen -> gen -> Signed.long -> gen
val lerchzeta : gen -> gen -> gen -> Signed.long -> gen
val zetahurwitz : gen -> gen -> Signed.long -> Signed.long -> gen
val rgx_to_ser : gen -> Signed.long -> gen
val rgx_to_ser_inexact : gen -> Signed.long -> gen
val gtoser : gen -> Signed.long -> Signed.long -> gen
val gtoser_prec : gen -> Signed.long -> Signed.long -> gen
val rfrac_to_ser : gen -> Signed.long -> gen
val rfrac_to_ser_i : gen -> Signed.long -> gen
val rfracrecip_to_ser_absolute : gen -> Signed.long -> gen
val rfracrecip : gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> Signed.long
val scalarser : gen -> Signed.long -> Signed.long -> gen
val sertoser : gen -> Signed.long -> gen
val toser_i : gen -> gen
val rgv_to_ser : gen -> Signed.long -> Signed.long -> gen
val ser0 : gen -> Signed.long -> gen -> Signed.long -> gen
val padic_to_q : gen -> gen
val padic_to_q_shallow : gen -> gen
val qpv_to_qv : gen -> gen
val rgc_rgv_mulrealsym : gen -> gen -> gen
val rgm_mulreal : gen -> gen -> gen
val rgx_cxeval : gen -> gen -> gen -> gen
val rgx_deflate_max : gen -> Signed.long Ctypes_static.ptr -> gen
val rgx_deflate_order : gen -> Signed.long
val rgx_degree : gen -> Signed.long -> Signed.long
val rgx_integ : gen -> gen
val rgxy_cxevalx : gen -> gen -> gen -> gen
val zx_deflate_order : gen -> Signed.long
val zx_deflate_max : gen -> Signed.long Ctypes_static.ptr -> gen
val ceil_safe : gen -> gen
val ceilr : gen -> gen
val centerlift : gen -> gen
val centerlift0 : gen -> Signed.long -> gen
val compo : gen -> Signed.long -> gen
val deg1pol : gen -> gen -> Signed.long -> gen
val deg1pol_shallow : gen -> gen -> Signed.long -> gen
val deg2pol_shallow : gen -> gen -> gen -> Signed.long -> gen
val degree : gen -> Signed.long
val denom : gen -> gen
val denom_i : gen -> gen
val denominator : gen -> gen -> gen

val deriv : gen -> Signed.long -> gen
(** [deriv] *)

val derivn : gen -> Signed.long -> Signed.long -> gen
val derivser : gen -> gen
val diffop : gen -> gen -> gen -> gen
val diffop0 : gen -> gen -> gen -> Signed.long -> gen
val diviiround : gen -> gen -> gen
val divrem : gen -> gen -> Signed.long -> gen
val floor_safe : gen -> gen
val gceil : gen -> gen
val gcvtoi : gen -> Signed.long Ctypes_static.ptr -> gen
val gdeflate : gen -> Signed.long -> Signed.long -> gen
val gdivent : gen -> gen -> gen
val gdiventgs : gen -> Signed.long -> gen
val gdiventsg : Signed.long -> gen -> gen
val gdiventres : gen -> gen -> gen
val gdivmod : gen -> gen -> gen Ctypes_static.ptr -> gen
val gdivround : gen -> gen -> gen
val gdvd : gen -> gen -> int
val geq : gen -> gen -> gen
val geval : gen -> gen
val gfloor : gen -> gen
val gtrunc2n : gen -> Signed.long -> gen
val gfrac : gen -> gen
val gge : gen -> gen -> gen
val ggrando : gen -> Signed.long -> gen
val ggt : gen -> gen -> gen
val gimag : gen -> gen
val gisexactzero : gen -> gen
val gle : gen -> gen -> gen
val glt : gen -> gen -> gen
val gmod : gen -> gen -> gen
val gmodgs : gen -> Signed.long -> gen
val gmodsg : Signed.long -> gen -> gen
val gmodulo : gen -> gen -> gen
val gmodulsg : Signed.long -> gen -> gen
val gmodulss : Signed.long -> Signed.long -> gen
val gne : gen -> gen -> gen
val gnot : gen -> gen
val gpolvar : gen -> gen
val gppadicprec : gen -> gen -> gen
val gppoldegree : gen -> Signed.long -> gen
val gprecision : gen -> Signed.long
val gpserprec : gen -> Signed.long -> gen
val greal : gen -> gen
val grndtoi : gen -> Signed.long Ctypes_static.ptr -> gen
val ground : gen -> gen
val gshift : gen -> Signed.long -> gen
val gsubst : gen -> Signed.long -> gen -> gen
val gsubstpol : gen -> gen -> gen -> gen
val gsubstvec : gen -> gen -> gen -> gen
val gtocol : gen -> gen
val gtocol0 : gen -> Signed.long -> gen
val gtocolrev : gen -> gen
val gtocolrev0 : gen -> Signed.long -> gen
val gtopoly : gen -> Signed.long -> gen
val gtopolyrev : gen -> Signed.long -> gen
val gtovec : gen -> gen
val gtovec0 : gen -> Signed.long -> gen
val gtovecrev : gen -> gen
val gtovecrev0 : gen -> Signed.long -> gen
val gtovecsmall : gen -> gen
val gtovecsmall0 : gen -> Signed.long -> gen
val gtrunc : gen -> gen
val gvar : gen -> Signed.long
val gvar2 : gen -> Signed.long
val hqfeval : gen -> gen -> gen
val imag_i : gen -> gen
val integ : gen -> Signed.long -> gen
val integser : gen -> gen
val ser_inv : gen -> gen
val iscomplex : gen -> int
val isexactzero : gen -> int
val isrationalzeroscalar : gen -> int
val isinexact : gen -> int
val isinexactreal : gen -> int
val isint : gen -> gen Ctypes_static.ptr -> int
val isrationalzero : gen -> int
val issmall : gen -> Signed.long Ctypes_static.ptr -> int
val lift : gen -> gen
val lift_shallow : gen -> gen
val lift0 : gen -> Signed.long -> gen
val liftall : gen -> gen
val liftall_shallow : gen -> gen
val liftint : gen -> gen
val liftint_shallow : gen -> gen
val liftpol : gen -> gen
val liftpol_shallow : gen -> gen
val mkcoln : Signed.long -> gen
val mkintn : Signed.long -> gen
val mkpoln : Signed.long -> gen
val mkvecn : Signed.long -> gen
val mkvecsmalln : Signed.long -> gen
val modrr_safe : gen -> gen -> gen
val modrr_i : gen -> gen -> gen -> gen
val mulreal : gen -> gen -> gen
val numer : gen -> gen
val numer_i : gen -> gen
val numerator : gen -> gen -> gen
val padicprec : gen -> gen -> Signed.long
val padicprec_relative : gen -> Signed.long
val polcoef : gen -> Signed.long -> Signed.long -> gen
val polcoef_i : gen -> Signed.long -> Signed.long -> gen
val poldegree : gen -> Signed.long -> Signed.long
val poleval : gen -> gen -> gen
val pollead : gen -> Signed.long -> gen
val precision : gen -> Signed.long
val qf_apply_rgm : gen -> gen -> gen
val qf_apply_zm : gen -> gen -> gen
val qfb_apply_zm : gen -> gen -> gen
val qfbil : gen -> gen -> gen -> gen
val qfeval : gen -> gen -> gen
val qfeval0 : gen -> gen -> gen -> gen
val qfevalb : gen -> gen -> gen -> gen
val qfnorm : gen -> gen -> gen
val real_i : gen -> gen
val round0 : gen -> gen Ctypes_static.ptr -> gen
val roundr : gen -> gen
val roundr_safe : gen -> gen
val scalarpol : gen -> Signed.long -> gen
val scalarpol_shallow : gen -> Signed.long -> gen
val ser_unscale : gen -> gen -> gen
val serprec : gen -> Signed.long -> Signed.long
val serreverse : gen -> gen
val simplify : gen -> gen
val simplify_shallow : gen -> gen
val tayl : gen -> Signed.long -> Signed.long -> gen
val trunc0 : gen -> gen Ctypes_static.ptr -> gen
val uu32toi : pari_ulong -> pari_ulong -> gen
val uu32toineg : pari_ulong -> pari_ulong -> gen
val vars_sort_inplace : gen -> gen
val vars_to_rgxv : gen -> gen
val variables_vecsmall : gen -> gen
val variables_vec : gen -> gen
val genus2red : gen -> gen -> gen
val genus2igusa : gen -> Signed.long -> gen
val gchar_conductor : gen -> gen -> gen
val gchar_identify : gen -> gen -> gen -> Signed.long -> gen
val gcharalgebraic : gen -> gen -> gen
val gcharduallog : gen -> gen -> gen
val gchareval : gen -> gen -> gen -> Signed.long -> gen
val gchari_lfun : gen -> gen -> gen -> gen
val gcharinit : gen -> gen -> Signed.long -> gen
val gcharisalgebraic : gen -> gen -> gen Ctypes_static.ptr -> int

val gcharlocal :
  gen -> gen -> gen -> Signed.long -> gen Ctypes_static.ptr -> gen

val gcharlog : gen -> gen -> Signed.long -> gen
val gcharnewprec : gen -> Signed.long -> gen
val is_gchar_group : gen -> int
val lfungchar : gen -> gen -> gen
val vecan_gchar : gen -> Signed.long -> Signed.long -> gen
val eulerf_gchar : gen -> gen -> Signed.long -> gen
val group_ident : gen -> gen -> Signed.long
val group_ident_trans : gen -> gen -> Signed.long

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
  hashtable Ctypes.structure Ctypes_static.ptr -> unit Ctypes_static.ptr -> gen

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
  (gen -> gen -> int) Ctypes_static.static_funptr ->
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

val hash_keys : hashtable Ctypes.structure Ctypes_static.ptr -> gen
val hash_values : hashtable Ctypes.structure Ctypes_static.ptr -> gen

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
val hash_gen : gen -> pari_ulong
val hash_zv : gen -> pari_ulong
val zx_hyperellred : gen -> gen Ctypes_static.ptr -> gen
val hyperellcharpoly : gen -> gen
val hyperellchangecurve : gen -> gen -> gen
val hyperelldisc : gen -> gen
val hyperellisoncurve : gen -> gen -> int
val hyperellminimaldisc : gen -> gen -> gen
val hyperellminimalmodel : gen -> gen Ctypes_static.ptr -> gen -> gen
val hyperellpadicfrobenius0 : gen -> gen -> Signed.long -> gen
val hyperellpadicfrobenius : gen -> pari_ulong -> Signed.long -> gen
val hyperellred : gen -> gen Ctypes_static.ptr -> gen
val nfhyperellpadicfrobenius : gen -> gen -> pari_ulong -> Signed.long -> gen
val hypergeom : gen -> gen -> gen -> Signed.long -> gen
val airy : gen -> Signed.long -> gen
val rgm_hnfall : gen -> gen Ctypes_static.ptr -> Signed.long -> gen
val zm_hnf : gen -> gen
val zm_hnf_knapsack : gen -> gen
val zm_hnfall : gen -> gen Ctypes_static.ptr -> Signed.long -> gen
val zm_hnfall_i : gen -> gen Ctypes_static.ptr -> Signed.long -> gen
val zm_hnfcenter : gen -> gen
val zm_hnflll : gen -> gen Ctypes_static.ptr -> int -> gen
val zv_extgcd : gen -> gen
val zv_snfall : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen
val zv_snf_group : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen
val zv_snf_rank_u : gen -> pari_ulong -> Signed.long
val zv_snf_trunc : gen -> unit
val zm_hnfmod : gen -> gen -> gen
val zm_hnfmodall : gen -> gen -> Signed.long -> gen
val zm_hnfmodall_i : gen -> gen -> Signed.long -> gen
val zm_hnfmodid : gen -> gen -> gen
val zm_hnfmodprime : gen -> gen -> gen
val zm_hnfperm : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen
val zm_snfclean : gen -> gen -> gen -> unit
val zm_snf : gen -> gen
val zm_snf_group : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen
val zm_snfall : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val zm_snfall_i :
  gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> Signed.long -> gen

val zv_snfclean : gen -> gen
val zpm_echelon : gen -> Signed.long -> gen -> gen -> gen
val gsmith : gen -> gen
val gsmithall : gen -> gen
val hnf : gen -> gen
val hnf_divscale : gen -> gen -> gen -> gen
val hnf_invscale : gen -> gen -> gen
val hnf_solve : gen -> gen -> gen
val hnf_invimage : gen -> gen -> gen
val hnfall : gen -> gen
val hnfdivide : gen -> gen -> int
val hnflll : gen -> gen
val hnfmerge_get_1 : gen -> gen -> gen
val hnfmod : gen -> gen -> gen
val hnfmodid : gen -> gen -> gen
val hnfperm : gen -> gen
val matfrobenius : gen -> Signed.long -> Signed.long -> gen
val mathnf0 : gen -> Signed.long -> gen
val matsnf0 : gen -> Signed.long -> gen
val smith : gen -> gen
val smithall : gen -> gen
val smithclean : gen -> gen
val snfrank : gen -> gen -> Signed.long
val zlm_echelon : gen -> Signed.long -> pari_ulong -> pari_ulong -> gen
val zv_snf_rank : gen -> pari_ulong -> Signed.long
val z_ecm : gen -> Signed.long -> Signed.long -> pari_ulong -> gen
val z_factor : gen -> gen
val z_factor_limit : gen -> pari_ulong -> gen
val z_factor_until : gen -> gen -> gen
val z_issmooth : gen -> pari_ulong -> Signed.long
val z_issmooth_fact : gen -> pari_ulong -> gen
val z_issquarefree : gen -> Signed.long
val z_pollardbrent : gen -> Signed.long -> Signed.long -> gen
val absz_factor : gen -> gen
val absz_factor_limit : gen -> pari_ulong -> gen
val absz_factor_limit_strict : gen -> pari_ulong -> gen Ctypes_static.ptr -> gen
val coreu : pari_ulong -> pari_ulong
val coreu_fact : gen -> pari_ulong
val factorint : gen -> Signed.long -> gen
val factoru : pari_ulong -> gen
val tridiv_boundu : pari_ulong -> pari_ulong
val ifac_isprime : gen -> int

val ifac_next :
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  int

val ifac_read :
  gen -> gen Ctypes_static.ptr -> Signed.long Ctypes_static.ptr -> int

val ifac_skip : gen -> unit
val ifac_start : gen -> int -> gen

val is_357_power :
  gen -> gen Ctypes_static.ptr -> pari_ulong Ctypes_static.ptr -> int

val is_pth_power :
  gen ->
  gen Ctypes_static.ptr ->
  forprime_t Ctypes.structure Ctypes_static.ptr ->
  pari_ulong ->
  int

val ispowerful : gen -> Signed.long
val maxomegau : pari_ulong -> Signed.long
val maxomegaoddu : pari_ulong -> Signed.long
val moebius : gen -> Signed.long
val moebiusu : pari_ulong -> Signed.long
val moebiusu_fact : gen -> Signed.long
val nextprime : gen -> gen
val precprime : gen -> gen
val radicalu : pari_ulong -> pari_ulong
val tridiv_bound : gen -> pari_ulong

val uis_357_power :
  pari_ulong ->
  pari_ulong Ctypes_static.ptr ->
  pari_ulong Ctypes_static.ptr ->
  int

val uis_357_powermod : pari_ulong -> pari_ulong Ctypes_static.ptr -> int
val unextprime : pari_ulong -> pari_ulong
val uprecprime : pari_ulong -> pari_ulong
val vecfactorsquarefreeu : pari_ulong -> pari_ulong -> gen
val vecfactorsquarefreeu_coprime : pari_ulong -> pari_ulong -> gen -> gen
val vecfactoru_i : pari_ulong -> pari_ulong -> gen
val vecfactoru : pari_ulong -> pari_ulong -> gen
val vecfactoroddu_i : pari_ulong -> pari_ulong -> gen
val vecfactoroddu : pari_ulong -> pari_ulong -> gen
val vecsquarefreeu : pari_ulong -> pari_ulong -> gen
val chk_gerepileupto : gen -> int
val copy_bin : gen -> genbin Ctypes.structure Ctypes_static.ptr
val copy_bin_canon : gen -> genbin Ctypes.structure Ctypes_static.ptr
val dbg_gerepile : pari_ulong -> unit
val dbg_gerepileupto : gen -> unit
val errname : gen -> gen
val gclone : gen -> gen
val gcloneref : gen -> gen
val gclone_refc : gen -> unit
val gcopy : gen -> gen
val gcopy_avma : gen -> pari_ulong Ctypes_static.ptr -> gen
val gcopy_lg : gen -> Signed.long -> gen
val gerepile : pari_ulong -> pari_ulong -> gen -> gen
val gerepileallsp : pari_ulong -> pari_ulong -> int -> unit

val gerepilecoeffssp :
  pari_ulong -> pari_ulong -> Signed.long Ctypes_static.ptr -> int -> unit

val gerepilemanysp :
  pari_ulong ->
  pari_ulong ->
  gen Ctypes_static.ptr Ctypes_static.ptr ->
  int ->
  unit

val getheap : unit -> gen
val gsizeword : gen -> Signed.long
val gsizebyte : gen -> Signed.long
val gunclone : gen -> unit
val gunclone_deep : gen -> unit
val listcopy : gen -> gen
val listinit : gen -> gen
val msgtimer : string -> unit
val name_numerr : string -> Signed.long
val new_chunk_resize : int -> unit
val newblock : int -> gen
val numerr_name : Signed.long -> string
val obj_check : gen -> Signed.long -> gen

val obj_checkbuild :
  gen -> Signed.long -> (gen -> gen) Ctypes_static.static_funptr -> gen

val obj_checkbuild_padicprec :
  gen ->
  Signed.long ->
  (gen -> Signed.long -> gen) Ctypes_static.static_funptr ->
  Signed.long ->
  gen

val obj_checkbuild_realprec :
  gen ->
  Signed.long ->
  (gen -> Signed.long -> gen) Ctypes_static.static_funptr ->
  Signed.long ->
  gen

val obj_checkbuild_prec :
  gen ->
  Signed.long ->
  (gen -> Signed.long -> gen) Ctypes_static.static_funptr ->
  (gen -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  gen

val obj_free : gen -> unit
val obj_init : Signed.long -> Signed.long -> gen
val obj_insert : gen -> Signed.long -> gen -> gen
val obj_insert_shallow : gen -> Signed.long -> gen -> gen
val obj_reinit : gen -> gen
val pari_add_function : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_add_defaults_module : entree Ctypes.structure Ctypes_static.ptr -> unit
val pari_close : unit -> unit
val pari_close_opts : pari_ulong -> unit
val pari_compile_str : string -> gen
val pari_daemon : unit -> int
val pari_err : int -> unit
val pari_err_last : unit -> gen
val pari_err2str : gen -> string
val pari_init_opts : int -> pari_ulong -> pari_ulong -> unit
val pari_init : int -> pari_ulong -> unit
val pari_stackcheck_init : unit Ctypes_static.ptr -> unit
val pari_sighandler : int -> unit
val pari_sig_init : (int -> unit) Ctypes_static.static_funptr -> unit

val pari_thread_alloc :
  pari_thread Ctypes.structure Ctypes_static.ptr -> int -> gen -> unit

val pari_thread_close : unit -> unit
val pari_thread_free : pari_thread Ctypes.structure Ctypes_static.ptr -> unit
val pari_thread_init : unit -> unit
val pari_thread_start : pari_thread Ctypes.structure Ctypes_static.ptr -> gen

val pari_thread_valloc :
  pari_thread Ctypes.structure Ctypes_static.ptr -> int -> int -> gen -> unit

val pari_version : unit -> gen
val pari_warn : int -> unit
val paristack_newrsize : pari_ulong -> unit
val paristack_resize : pari_ulong -> unit
val paristack_setsize : int -> int -> unit
val parivstack_resize : pari_ulong -> unit
val parivstack_reset : unit -> unit
val setalldebug : Signed.long -> unit
val setdebug : string -> Signed.long -> gen
val shiftaddress : gen -> Signed.long -> unit
val shiftaddress_canon : gen -> Signed.long -> unit
val timer : unit -> Signed.long
val timer_delay : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long
val timer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val timer_printf :
  pari_timer Ctypes.structure Ctypes_static.ptr -> string -> unit

val timer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit
val timer2 : unit -> Signed.long
val trap0 : string -> gen -> gen -> gen

val traverseheap :
  (gen -> unit Ctypes_static.ptr -> unit) Ctypes_static.static_funptr ->
  unit Ctypes_static.ptr ->
  unit

val walltimer_start : pari_timer Ctypes.structure Ctypes_static.ptr -> unit

val walltimer_delay :
  pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long

val walltimer_get : pari_timer Ctypes.structure Ctypes_static.ptr -> Signed.long
val contfraceval : gen -> gen -> Signed.long -> gen
val contfracinit : gen -> Signed.long -> gen

val intcirc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  gen ->
  Signed.long ->
  gen

val intfuncinit :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  Signed.long ->
  gen

val intnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  gen ->
  Signed.long ->
  gen

val intnumgauss :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  gen ->
  Signed.long ->
  gen

val intnumgaussinit : Signed.long -> Signed.long -> gen
val intnuminit : gen -> gen -> Signed.long -> Signed.long -> gen

val intnumosc :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen ->
  Signed.long ->
  gen

val intnumromb :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  Signed.long ->
  gen

val intnumromb_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  Signed.long ->
  gen

val prodeulerrat : gen -> gen -> Signed.long -> Signed.long -> gen
val prodnumrat : gen -> Signed.long -> Signed.long -> gen
val quodif : gen -> Signed.long -> gen
val sumeulerrat : gen -> gen -> Signed.long -> Signed.long -> gen

val sumnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen

val sumnumap :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen

val sumnumapinit : gen -> Signed.long -> gen
val sumnuminit : gen -> Signed.long -> gen
val sumnumlagrangeinit : gen -> gen -> Signed.long -> gen

val sumnumlagrange :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen

val sumnummonien :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen

val sumnummonieninit : gen -> gen -> gen -> Signed.long -> gen
val sumnumrat : gen -> gen -> Signed.long -> gen

val sumnumsidi :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen ->
  float ->
  Signed.long ->
  gen

val z_isanypower : gen -> gen Ctypes_static.ptr -> Signed.long
val z_ispow2 : gen -> Signed.long
val z_ispowerall : gen -> pari_ulong -> gen Ctypes_static.ptr -> Signed.long
val z_issquareall : gen -> gen Ctypes_static.ptr -> Signed.long
val zn_ispower : gen -> gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val zn_issquare : gen -> gen -> Signed.long
val zp_issquare : gen -> gen -> Signed.long
val gisanypower : gen -> gen Ctypes_static.ptr -> Signed.long
val gissquare : gen -> gen
val gissquareall : gen -> gen Ctypes_static.ptr -> gen
val ispolygonal : gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val ispower : gen -> gen -> gen Ctypes_static.ptr -> Signed.long
val isprimepower : gen -> gen Ctypes_static.ptr -> Signed.long
val ispseudoprimepower : gen -> gen Ctypes_static.ptr -> Signed.long
val issquare : gen -> Signed.long
val issquareall : gen -> gen Ctypes_static.ptr -> Signed.long
val sqrtint : gen -> gen
val sqrtint0 : gen -> gen Ctypes_static.ptr -> gen
val uisprimepower : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long
val uissquare : pari_ulong -> Signed.long
val uissquareall : pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val ulogintall :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> Signed.long

val padicfields0 : gen -> gen -> Signed.long -> gen
val padicfields : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val bnrclassfield : gen -> gen -> Signed.long -> Signed.long -> gen
val rnfkummer : gen -> gen -> Signed.long -> gen
val is_linit : gen -> Signed.long
val ldata_get_an : gen -> gen
val ldata_get_dual : gen -> gen
val ldata_get_gammavec : gen -> gen
val ldata_get_degree : gen -> Signed.long
val ldata_get_k : gen -> gen
val ldata_get_k1 : gen -> gen
val ldata_get_conductor : gen -> gen
val ldata_get_rootno : gen -> gen
val ldata_get_residue : gen -> gen
val ldata_get_type : gen -> Signed.long
val ldata_isreal : gen -> Signed.long
val linit_get_type : gen -> Signed.long
val linit_get_ldata : gen -> gen
val linit_get_tech : gen -> gen
val lfun_get_domain : gen -> gen
val lfun_get_dom : gen -> gen
val lfun_get_factgammavec : gen -> gen
val lfun_get_step : gen -> gen
val lfun_get_pol : gen -> gen
val lfun_get_residue : gen -> gen
val lfun_get_k2 : gen -> gen
val lfun_get_w2 : gen -> gen
val lfun_get_expot : gen -> gen
val lfun_get_bitprec : gen -> Signed.long
val lfun : gen -> gen -> Signed.long -> gen
val lfun0 : gen -> gen -> Signed.long -> Signed.long -> gen
val lfuncheckfeq : gen -> gen -> Signed.long -> Signed.long
val lfunconductor : gen -> gen -> Signed.long -> Signed.long -> gen
val lfuncost : gen -> gen -> Signed.long -> Signed.long -> gen
val lfuncost0 : gen -> gen -> Signed.long -> Signed.long -> gen
val lfuncreate : gen -> gen
val lfundual : gen -> Signed.long -> gen
val lfuneuler : gen -> gen -> Signed.long -> gen
val lfunparams : gen -> Signed.long -> gen
val lfunan : gen -> Signed.long -> Signed.long -> gen
val lfunhardy : gen -> gen -> Signed.long -> gen
val lfuninit : gen -> gen -> Signed.long -> Signed.long -> gen
val lfuninit0 : gen -> gen -> Signed.long -> Signed.long -> gen
val lfuninit_make : Signed.long -> gen -> gen -> gen -> gen
val lfunlambda : gen -> gen -> Signed.long -> gen
val lfunlambda0 : gen -> gen -> Signed.long -> Signed.long -> gen
val lfunmisc_to_ldata : gen -> gen
val lfunmisc_to_ldata_shallow : gen -> gen
val lfunmisc_to_ldata_shallow_i : gen -> gen
val lfunorderzero : gen -> Signed.long -> Signed.long -> Signed.long
val lfunprod_get_fact : gen -> gen
val lfunrootno : gen -> Signed.long -> gen
val lfunrootres : gen -> Signed.long -> gen
val lfunrtopoles : gen -> gen
val lfunshift : gen -> gen -> Signed.long -> Signed.long -> gen
val lfuntwist : gen -> gen -> Signed.long -> gen
val lfuntheta : gen -> gen -> Signed.long -> Signed.long -> gen
val lfunthetacost0 : gen -> gen -> Signed.long -> Signed.long -> Signed.long
val lfunthetacost : gen -> gen -> Signed.long -> Signed.long -> Signed.long
val lfunthetainit : gen -> gen -> Signed.long -> Signed.long -> gen
val lfunthetacheckinit : gen -> gen -> Signed.long -> Signed.long -> gen
val lfunzeros : gen -> gen -> Signed.long -> Signed.long -> gen
val sdomain_isincl : float -> gen -> gen -> int
val theta_get_an : gen -> gen
val theta_get_k : gen -> gen
val theta_get_r : gen -> gen
val theta_get_bitprec : gen -> Signed.long
val theta_get_m : gen -> Signed.long
val theta_get_tdom : gen -> gen
val theta_get_isqrtn : gen -> gen
val vgaeasytheta : gen -> int
val znchargauss : gen -> gen -> gen -> Signed.long -> gen
val dirzetak : gen -> gen -> gen
val ellmoddegree : gen -> gen
val eta_zxn : Signed.long -> Signed.long -> gen
val eta_product_zxn : gen -> Signed.long -> gen

val etaquotype :
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val galois_get_conj : gen -> gen
val ldata_vecan : gen -> Signed.long -> Signed.long -> gen
val ldata_newprec : gen -> Signed.long -> gen
val lfunabelianrelinit : gen -> gen -> gen -> Signed.long -> Signed.long -> gen
val lfunartin : gen -> gen -> gen -> Signed.long -> Signed.long -> gen
val lfundiv : gen -> gen -> Signed.long -> gen
val lfunellmfpeters : gen -> Signed.long -> gen
val lfunetaquo : gen -> gen
val lfungenus2 : gen -> gen
val lfunmfspec : gen -> Signed.long -> gen
val lfunmul : gen -> gen -> Signed.long -> gen
val lfunqf : gen -> Signed.long -> gen
val lfunsympow : gen -> pari_ulong -> gen
val lfunzetakinit : gen -> gen -> Signed.long -> Signed.long -> gen
val qfiseven : gen -> Signed.long
val lfunquadneg : Signed.long -> Signed.long -> gen
val zm_lll_norms : gen -> float -> Signed.long -> gen Ctypes_static.ptr -> gen
val kerint : gen -> gen
val lll : gen -> gen
val lllfp : gen -> float -> Signed.long -> gen
val lllgen : gen -> gen
val lllgram : gen -> gen
val lllgramgen : gen -> gen
val lllgramint : gen -> gen
val lllgramkerim : gen -> gen
val lllgramkerimgen : gen -> gen
val lllint : gen -> gen
val lllintpartial : gen -> gen
val lllintpartial_inplace : gen -> gen
val lllkerim : gen -> gen
val lllkerimgen : gen -> gen
val matkerint0 : gen -> Signed.long -> gen
val qflll0 : gen -> Signed.long -> gen
val qflllgram0 : gen -> Signed.long -> gen
val gtomap : gen -> gen
val mapdelete : gen -> gen -> unit
val mapdomain : gen -> gen
val mapdomain_shallow : gen -> gen
val mapget : gen -> gen -> gen
val mapisdefined : gen -> gen -> gen Ctypes_static.ptr -> int
val mapput : gen -> gen -> gen -> unit
val maptomat : gen -> gen
val maptomat_shallow : gen -> gen
val matpermanent : gen -> gen
val zm_permanent : gen -> gen
val dbllemma526 : float -> float -> float -> float -> float
val dblcoro526 : float -> float -> float -> float
val gammamellininv : gen -> gen -> Signed.long -> Signed.long -> gen
val gammamellininvasymp : gen -> Signed.long -> Signed.long -> gen
val gammamellininvinit : gen -> Signed.long -> Signed.long -> gen
val gammamellininvrt : gen -> gen -> Signed.long -> gen
val member_a1 : gen -> gen
val member_a2 : gen -> gen
val member_a3 : gen -> gen
val member_a4 : gen -> gen
val member_a6 : gen -> gen
val member_area : gen -> gen
val member_b2 : gen -> gen
val member_b4 : gen -> gen
val member_b6 : gen -> gen
val member_b8 : gen -> gen
val member_bid : gen -> gen
val member_bnf : gen -> gen
val member_c4 : gen -> gen
val member_c6 : gen -> gen
val member_clgp : gen -> gen
val member_codiff : gen -> gen
val member_cyc : gen -> gen
val member_diff : gen -> gen
val member_disc : gen -> gen
val member_e : gen -> gen
val member_eta : gen -> gen
val member_f : gen -> gen
val member_fu : gen -> gen
val member_gen : gen -> gen
val member_group : gen -> gen
val member_index : gen -> gen
val member_j : gen -> gen
val member_mod : gen -> gen
val member_nf : gen -> gen
val member_no : gen -> gen
val member_omega : gen -> gen
val member_orders : gen -> gen
val member_p : gen -> gen
val member_pol : gen -> gen
val member_polabs : gen -> gen
val member_reg : gen -> gen
val member_r1 : gen -> gen
val member_r2 : gen -> gen
val member_roots : gen -> gen
val member_sign : gen -> gen
val member_t2 : gen -> gen
val member_tate : gen -> gen
val member_tu : gen -> gen
val member_zk : gen -> gen
val member_zkst : gen -> gen
val mf_get_m : gen -> gen
val mf_get_mindex : gen -> gen
val mf_get_minv : gen -> gen
val mf_get_basis : gen -> gen
val mf_get_dim : gen -> Signed.long
val mf_get_e : gen -> gen
val mf_get_fields : gen -> gen
val mf_get_newforms : gen -> gen
val mf_get_space : gen -> Signed.long
val mf_get_s : gen -> gen
val mfcusp_get_vmjd : gen -> gen
val mfnew_get_vj : gen -> gen
val qab_tracerel : gen -> Signed.long -> gen -> gen
val qabm_tracerel : gen -> Signed.long -> gen -> gen
val qabv_tracerel : gen -> Signed.long -> gen -> gen
val qab_trace_init : Signed.long -> Signed.long -> gen -> gen -> gen
val checkmf : gen -> gen
val checkmf_i : gen -> int
val getcache : unit -> gen
val hclassno6u : pari_ulong -> pari_ulong
val hclassno6u_no_cache : pari_ulong -> pari_ulong
val lfunmf : gen -> gen -> Signed.long -> gen
val mfdelta : unit -> gen
val mfeh : gen -> gen
val mfek : Signed.long -> gen
val mftheta : gen -> gen
val mf_get_chi : gen -> gen
val mf_get_n : gen -> Signed.long
val mf_get_nk : gen -> gen
val mf_get_field : gen -> gen
val mf_get_gn : gen -> gen
val mf_get_gk : gen -> gen
val mf_get_k : gen -> Signed.long
val mf_get_r : gen -> Signed.long
val mf_get_type : gen -> Signed.long
val mfatkin : gen -> gen -> gen
val mfatkineigenvalues : gen -> Signed.long -> Signed.long -> gen
val mfatkininit : gen -> Signed.long -> Signed.long -> gen
val mfbasis : gen -> Signed.long -> gen
val mfbd : gen -> Signed.long -> gen
val mfbracket : gen -> gen -> Signed.long -> gen
val mfcharorder : gen -> Signed.long
val mfcharmodulus : gen -> Signed.long
val mfcharpol : gen -> gen
val mfcoef : gen -> Signed.long -> gen
val mfcoefs : gen -> Signed.long -> Signed.long -> gen
val mfconductor : gen -> gen -> Signed.long
val mfcosets : gen -> gen
val mfcuspdim : Signed.long -> Signed.long -> gen -> Signed.long
val mfcuspisregular : gen -> gen -> Signed.long
val mfcusps : gen -> gen
val mfcuspval : gen -> gen -> gen -> Signed.long -> gen
val mfcuspwidth : gen -> gen -> Signed.long
val mfderiv : gen -> Signed.long -> gen
val mfderive2 : gen -> Signed.long -> gen
val mfdescribe : gen -> gen Ctypes_static.ptr -> gen
val mfdim : gen -> Signed.long -> gen
val mfdiv : gen -> gen -> gen
val mfdiv_val : gen -> gen -> Signed.long -> gen
val mfeigenbasis : gen -> gen
val mfeigensearch : gen -> gen -> gen
val mfeisenstein : Signed.long -> gen -> gen -> gen
val mfeisensteindim : Signed.long -> Signed.long -> gen -> Signed.long
val mfembed : gen -> gen -> gen
val mfembed0 : gen -> gen -> Signed.long -> gen
val mfeval : gen -> gen -> gen -> Signed.long -> gen
val mffields : gen -> gen
val mffromell : gen -> gen
val mffrometaquo : gen -> Signed.long -> gen
val mffromlfun : gen -> Signed.long -> gen
val mffromqf : gen -> gen -> gen
val mffulldim : Signed.long -> Signed.long -> gen -> Signed.long
val mfgaloisprojrep : gen -> gen -> Signed.long -> gen
val mfgaloistype : gen -> gen -> gen
val mfhecke : gen -> gen -> Signed.long -> gen
val mfheckemat : gen -> gen -> gen
val mfinit : gen -> Signed.long -> gen
val mfiscm : gen -> gen
val mfiscuspidal : gen -> gen -> Signed.long
val mfisequal : gen -> gen -> Signed.long -> Signed.long
val mfisetaquo : gen -> Signed.long -> gen
val mfkohnenbasis : gen -> gen
val mfkohnenbijection : gen -> gen
val mfkohneneigenbasis : gen -> gen -> gen
val mflinear : gen -> gen -> gen
val mfmanin : gen -> Signed.long -> gen
val mfmatembed : gen -> gen -> gen
val mfmul : gen -> gen -> gen
val mfnewdim : Signed.long -> Signed.long -> gen -> Signed.long
val mfolddim : Signed.long -> Signed.long -> gen -> Signed.long
val mfparams : gen -> gen
val mfperiodpol : gen -> gen -> Signed.long -> Signed.long -> gen
val mfperiodpolbasis : Signed.long -> Signed.long -> gen
val mfpetersson : gen -> gen -> gen
val mfpow : gen -> Signed.long -> gen
val mfsearch : gen -> gen -> Signed.long -> gen
val mfshift : gen -> Signed.long -> gen
val mfshimura : gen -> gen -> Signed.long -> gen

val mfslashexpansion :
  gen ->
  gen ->
  gen ->
  Signed.long ->
  Signed.long ->
  gen Ctypes_static.ptr ->
  Signed.long ->
  gen

val mfspace : gen -> gen -> Signed.long
val mfsplit : gen -> Signed.long -> Signed.long -> gen
val mfsturm : gen -> Signed.long
val mfsturmngk : Signed.long -> gen -> Signed.long
val mfsturmnk : Signed.long -> Signed.long -> Signed.long
val mfsturm_mf : gen -> Signed.long
val mfsymboleval : gen -> gen -> gen -> Signed.long -> gen
val mfsymbol : gen -> gen -> Signed.long -> gen
val mftaylor : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val mftobasis : gen -> gen -> Signed.long -> gen
val mftobasises : gen -> gen -> gen
val mftocol : gen -> Signed.long -> Signed.long -> gen
val mftocoset : pari_ulong -> gen -> gen -> gen
val mftonew : gen -> gen -> gen
val mftraceform : gen -> Signed.long -> gen
val mftwist : gen -> gen -> gen
val mfvecembed : gen -> gen -> gen
val mfvectomat : gen -> Signed.long -> Signed.long -> gen
val fl_inv : pari_ulong -> pari_ulong -> pari_ulong
val fl_invsafe : pari_ulong -> pari_ulong -> pari_ulong

val fp_ratlift :
  gen ->
  gen ->
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  int

val zm2_mul : gen -> gen -> gen
val abscmpii : gen -> gen -> int
val abscmprr : gen -> gen -> int
val absequalii : gen -> gen -> int
val addii_sign : gen -> Signed.long -> gen -> Signed.long -> gen
val addir_sign : gen -> Signed.long -> gen -> Signed.long -> gen
val addmulii : gen -> gen -> gen -> gen
val addmulii_inplace : gen -> gen -> gen -> gen
val addrr_sign : gen -> Signed.long -> gen -> Signed.long -> gen
val addsi_sign : Signed.long -> gen -> Signed.long -> gen
val addsr : Signed.long -> gen -> gen
val addui_sign : pari_ulong -> gen -> Signed.long -> gen
val addumului : pari_ulong -> pari_ulong -> gen -> gen
val affir : gen -> gen -> unit
val affrr : gen -> gen -> unit
val bezout : gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val cbezout :
  Signed.long ->
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val cgcd : Signed.long -> Signed.long -> Signed.long
val clcm : Signed.long -> Signed.long -> Signed.long
val cmpii : gen -> gen -> int
val cmprr : gen -> gen -> int
val dblexpo : float -> Signed.long
val dblmantissa : float -> pari_ulong
val dbltor : float -> gen
val diviiexact : gen -> gen -> gen
val divir : gen -> gen -> gen
val divis : gen -> Signed.long -> gen
val divis_rem : gen -> Signed.long -> Signed.long Ctypes_static.ptr -> gen
val absdiviu_rem : gen -> pari_ulong -> pari_ulong Ctypes_static.ptr -> gen
val diviuuexact : gen -> pari_ulong -> pari_ulong -> gen
val diviuexact : gen -> pari_ulong -> gen
val divri : gen -> gen -> gen
val divrr : gen -> gen -> gen
val divrs : gen -> Signed.long -> gen
val divru : gen -> pari_ulong -> gen
val divsi : Signed.long -> gen -> gen
val divsr : Signed.long -> gen -> gen
val divur : pari_ulong -> gen -> gen
val dvmdii : gen -> gen -> gen Ctypes_static.ptr -> gen
val equalii : gen -> gen -> int
val equalrr : gen -> gen -> int
val floorr : gen -> gen
val gcdii : gen -> gen -> gen
val halfgcdii : gen -> gen -> gen
val int2n : Signed.long -> gen
val int2u : pari_ulong -> gen
val int2um1 : pari_ulong -> gen
val int_normalize : gen -> Signed.long -> gen
val invmod : gen -> gen -> gen Ctypes_static.ptr -> int
val invmod2bil : pari_ulong -> pari_ulong
val invr : gen -> gen
val mantissa_real : gen -> Signed.long Ctypes_static.ptr -> gen
val modii : gen -> gen -> gen
val modiiz : gen -> gen -> gen -> unit
val mulii : gen -> gen -> gen
val mulir : gen -> gen -> gen
val mulrr : gen -> gen -> gen
val mulsi : Signed.long -> gen -> gen
val mulsr : Signed.long -> gen -> gen
val mulss : Signed.long -> Signed.long -> gen
val mului : pari_ulong -> gen -> gen
val mulur : pari_ulong -> gen -> gen
val muluu : pari_ulong -> pari_ulong -> gen
val muluui : pari_ulong -> pari_ulong -> gen -> gen
val pari_kernel_close : unit -> unit
val pari_kernel_init : unit -> unit
val pari_kernel_version : unit -> string
val remi2n : gen -> Signed.long -> gen
val rtodbl : gen -> float
val shifti : gen -> Signed.long -> gen
val sqri : gen -> gen
val sqrr : gen -> gen
val sqrs : Signed.long -> gen
val sqrtr_abs : gen -> gen
val sqrtremi : gen -> gen Ctypes_static.ptr -> gen
val sqru : pari_ulong -> gen
val subsr : Signed.long -> gen -> gen
val truedvmdii : gen -> gen -> gen Ctypes_static.ptr -> gen
val truedvmdis : gen -> Signed.long -> gen Ctypes_static.ptr -> gen
val truedvmdsi : Signed.long -> gen -> gen Ctypes_static.ptr -> gen
val trunc2nr : gen -> Signed.long -> gen
val mantissa2nr : gen -> Signed.long -> gen
val truncr : gen -> gen
val ugcd : pari_ulong -> pari_ulong -> pari_ulong
val ulcm : pari_ulong -> pari_ulong -> pari_ulong
val umodiu : gen -> pari_ulong -> pari_ulong
val vals : pari_ulong -> Signed.long
val fpc_ratlift : gen -> gen -> gen -> gen -> gen -> gen
val fpm_ratlift : gen -> gen -> gen -> gen -> gen -> gen
val fpx_ratlift : gen -> gen -> gen -> gen -> gen -> gen
val qxqx_gcd : gen -> gen -> gen -> gen
val zxqx_gcd : gen -> gen -> gen -> gen
val nffactor : gen -> gen -> gen
val nffactormod : gen -> gen -> gen -> gen
val nfgcd : gen -> gen -> gen -> gen -> gen
val nfgcd_all : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val nfissquarefree : gen -> gen -> int
val nfroots : gen -> gen -> gen
val nfroots_if_split : gen Ctypes_static.ptr -> gen -> gen
val nfrootsof1 : gen -> gen
val polfnf : gen -> gen -> gen
val rnfabelianconjgen : gen -> gen -> gen
val rnfisabelian : gen -> gen -> Signed.long

val forpart :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  Signed.long ->
  gen ->
  gen ->
  unit

val forpart_init :
  forpart_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  gen ->
  gen ->
  unit

val forpart_next : forpart_t Ctypes.structure Ctypes_static.ptr -> gen
val forpart_prev : forpart_t Ctypes.structure Ctypes_static.ptr -> gen
val numbpart : gen -> gen
val partitions : Signed.long -> gen -> gen -> gen

val forperm :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  gen ->
  unit

val forperm_init : forperm_t Ctypes.structure Ctypes_static.ptr -> gen -> unit
val forperm_next : forperm_t Ctypes.structure Ctypes_static.ptr -> gen

val forallsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> Signed.long -> unit

val forksubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr ->
  Signed.long ->
  Signed.long ->
  unit

val forsubset_next : forsubset_t Ctypes.structure Ctypes_static.ptr -> gen

val forsubset_init :
  forsubset_t Ctypes.structure Ctypes_static.ptr -> gen -> unit

val glambertw : gen -> Signed.long -> Signed.long -> gen
val mplambertw : gen -> Signed.long -> gen
val mplambertx : gen -> Signed.long -> gen
val mplambertx_logx : gen -> gen -> Signed.long -> gen
val mplambertxlogx_x : gen -> gen -> Signed.long -> gen
val z_to_perm : Signed.long -> gen -> gen
val abelian_group : gen -> gen
val conjclasses_repr : gen -> Signed.long -> gen
val cyc_pow : gen -> Signed.long -> gen
val cyc_pow_perm : gen -> Signed.long -> gen
val cyclicgroup : gen -> Signed.long -> gen
val dicyclicgroup : gen -> gen -> Signed.long -> Signed.long -> gen
val group_abelianhnf : gen -> gen -> gen
val group_abeliansnf : gen -> gen -> gen
val group_domain : gen -> Signed.long
val group_elts : gen -> Signed.long -> gen
val group_export : gen -> Signed.long -> gen
val group_export_gap : gen -> gen
val group_export_magma : gen -> gen
val group_isa4s4 : gen -> Signed.long
val group_isabelian : gen -> Signed.long
val group_leftcoset : gen -> gen -> gen
val group_order : gen -> Signed.long
val group_perm_normalize : gen -> gen -> Signed.long
val group_quotient : gen -> gen -> gen
val group_rightcoset : gen -> gen -> gen
val group_set : gen -> Signed.long -> gen
val group_subgroup_is_faithful : gen -> gen -> int
val group_subgroup_isnormal : gen -> gen -> Signed.long
val group_subgroups : gen -> gen
val groupelts_solvablesubgroups : gen -> gen
val group_to_cc : gen -> gen
val groupelts_abelian_group : gen -> gen
val groupelts_center : gen -> gen
val groupelts_conj_set : gen -> gen -> gen
val groupelts_conjclasses : gen -> Signed.long Ctypes_static.ptr -> gen
val groupelts_exponent : gen -> Signed.long
val groupelts_quotient : gen -> gen -> gen
val groupelts_set : gen -> Signed.long -> gen
val groupelts_to_group : gen -> gen
val numtoperm : Signed.long -> gen -> gen
val perm_commute : gen -> gen -> int
val perm_cycles : gen -> gen
val perm_order : gen -> gen
val perm_orderu : gen -> pari_ulong
val perm_pow : gen -> gen -> gen
val perm_powu : gen -> pari_ulong -> gen
val perm_sign : gen -> Signed.long
val perm_to_gap : gen -> gen
val perm_to_z : gen -> gen
val permcycles : gen -> gen
val permorder : gen -> gen
val permsign : gen -> Signed.long
val permtonum : gen -> gen
val quotient_group : gen -> gen -> gen
val quotient_groupelts : gen -> gen
val quotient_perm : gen -> gen -> gen
val quotient_subgroup_lift : gen -> gen -> gen -> gen
val subgroups_tableset : gen -> Signed.long -> gen
val tableset_find_index : gen -> gen -> Signed.long
val trivialgroup : unit -> gen
val vec_insert : gen -> Signed.long -> gen -> gen
val vec_is1to1 : gen -> int
val vec_isconst : gen -> int
val vecperm_orbits : gen -> Signed.long -> gen
val vecsmall_duplicate : gen -> Signed.long
val vecsmall_duplicate_sorted : gen -> Signed.long
val vecsmall_indexsort : gen -> gen
val vecsmall_is1to1 : gen -> int
val vecsmall_isconst : gen -> int
val vecsmall_sort : gen -> unit
val vecsmall_uniq : gen -> gen
val vecsmall_uniq_sorted : gen -> gen
val vecsmall_counting_indexsort : gen -> Signed.long -> gen
val vecsmall_counting_sort : gen -> Signed.long -> unit
val vecsmall_counting_uniq : gen -> Signed.long -> gen
val vecvecsmall_indexsort : gen -> gen
val vecvecsmall_max : gen -> Signed.long
val vecvecsmall_search : gen -> gen -> Signed.long
val vecvecsmall_sort : gen -> gen
val vecvecsmall_sort_inplace : gen -> gen Ctypes_static.ptr -> unit
val vecvecsmall_sort_shallow : gen -> gen
val vecvecsmall_sort_uniq : gen -> gen
val mt_broadcast : gen -> unit
val mt_nbthreads : unit -> Signed.long
val mt_queue_end : pari_mt Ctypes.structure Ctypes_static.ptr -> unit

val mt_queue_get :
  pari_mt Ctypes.structure Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  gen

val mt_queue_start : pari_mt Ctypes.structure Ctypes_static.ptr -> gen -> unit

val mt_queue_start_lim :
  pari_mt Ctypes.structure Ctypes_static.ptr -> gen -> Signed.long -> unit

val mt_queue_submit :
  pari_mt Ctypes.structure Ctypes_static.ptr -> Signed.long -> gen -> unit

val mt_sigint_block : unit -> unit
val mt_sigint_unblock : unit -> unit
val pari_mt_init : unit -> unit
val pari_mt_close : unit -> unit
val subcyclopclgp : gen -> gen -> Signed.long -> gen
val subcycloiwasawa : gen -> gen -> Signed.long -> gen
val subcyclohminus : gen -> gen -> gen

val color_to_rgb :
  gen ->
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
  gen -> gen -> gen -> Signed.long -> Signed.long -> Signed.long -> gen

val parplothexport :
  gen -> gen -> gen -> gen -> Signed.long -> Signed.long -> Signed.long -> gen

val plotbox : Signed.long -> gen -> gen -> Signed.long -> unit
val plotclip : Signed.long -> unit
val plotcolor : Signed.long -> gen -> gen
val plotcopy : Signed.long -> Signed.long -> gen -> gen -> Signed.long -> unit
val plotcursor : Signed.long -> gen
val plotdraw : gen -> Signed.long -> unit
val plotexport : gen -> gen -> Signed.long -> gen

val ploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  gen

val plothexport :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  gen

val plothraw : gen -> gen -> Signed.long -> gen
val plothrawexport : gen -> gen -> gen -> Signed.long -> gen
val plothsizes : Signed.long -> gen
val plotinit : Signed.long -> gen -> gen -> Signed.long -> unit
val plotkill : Signed.long -> unit
val plotline : Signed.long -> gen -> gen -> unit
val plotlines : Signed.long -> gen -> gen -> Signed.long -> unit
val plotlinetype : Signed.long -> Signed.long -> unit
val plotmove : Signed.long -> gen -> gen -> unit
val plotpoints : Signed.long -> gen -> gen -> unit
val plotpointsize : Signed.long -> gen -> unit
val plotpointtype : Signed.long -> Signed.long -> unit
val plotrbox : Signed.long -> gen -> gen -> Signed.long -> unit

val plotrecth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  Signed.long ->
  gen ->
  gen ->
  pari_ulong ->
  Signed.long ->
  Signed.long ->
  gen

val plotrecthraw : Signed.long -> gen -> Signed.long -> gen
val plotrline : Signed.long -> gen -> gen -> unit
val plotrmove : Signed.long -> gen -> gen -> unit
val plotrpoint : Signed.long -> gen -> gen -> unit
val plotscale : Signed.long -> gen -> gen -> gen -> gen -> unit
val plotstring : Signed.long -> string -> Signed.long -> unit
val psdraw : gen -> Signed.long -> unit

val psploth :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  gen

val psplothraw : gen -> gen -> Signed.long -> gen

val rect2ps :
  gen -> gen -> gen -> pari_plot Ctypes.structure Ctypes_static.ptr -> string

val rect2ps_i :
  gen ->
  gen ->
  gen ->
  pari_plot Ctypes.structure Ctypes_static.ptr ->
  int ->
  string

val rect2svg :
  gen -> gen -> gen -> pari_plot Ctypes.structure Ctypes_static.ptr -> string

val pariplot :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  gen ->
  gen ->
  Signed.long ->
  unit

val zx_zp_root : gen -> gen -> gen -> Signed.long -> gen
val zp_appr : gen -> gen -> gen
val cmp_padic : gen -> gen -> int
val factorpadic : gen -> gen -> Signed.long -> gen
val gdeuc : gen -> gen -> gen
val grem : gen -> gen -> gen
val padicappr : gen -> gen -> gen
val poldivrem : gen -> gen -> gen Ctypes_static.ptr -> gen
val polrootspadic : gen -> gen -> Signed.long -> gen
val flv_factorback : gen -> gen -> pari_ulong -> pari_ulong
val flxqv_factorback : gen -> gen -> gen -> pari_ulong -> gen
val fpv_factorback : gen -> gen -> gen -> gen
val fqv_factorback : gen -> gen -> gen -> gen -> gen
val q_content : gen -> gen
val q_content_safe : gen -> gen
val q_denom : gen -> gen
val q_denom_safe : gen -> gen
val q_div_to_int : gen -> gen -> gen
val q_gcd : gen -> gen -> gen
val q_mul_to_int : gen -> gen -> gen
val q_muli_to_int : gen -> gen -> gen
val q_primitive_part : gen -> gen Ctypes_static.ptr -> gen
val q_primpart : gen -> gen
val q_remove_denom : gen -> gen Ctypes_static.ptr -> gen
val q_factor : gen -> gen
val q_factor_limit : gen -> pari_ulong -> gen

val rg_type :
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rgc_type :
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_rescale_to_int : gen -> gen

val rgm_type :
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgm_type2 :
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type :
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgv_type2 :
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_rg_type :
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_chinese_coprime : gen -> gen -> gen -> gen -> gen -> gen
val rgx_disc : gen -> gen

val rgx_extgcd :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val rgx_extgcd_simple :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val rgx_gcd : gen -> gen -> gen
val rgx_gcd_simple : gen -> gen -> gen
val rgx_halfgcd : gen -> gen -> gen

val rgx_halfgcd_all :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val rgx_rescale_to_int : gen -> gen
val rgx_resultant_all : gen -> gen -> gen Ctypes_static.ptr -> gen
val rgx_sturmpart : gen -> gen -> Signed.long
val rgx_sylvestermatrix : gen -> gen -> gen

val rgx_type :
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type2 :
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type3 :
  gen ->
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  Signed.long

val rgx_type_decode :
  Signed.long ->
  Signed.long Ctypes_static.ptr ->
  Signed.long Ctypes_static.ptr ->
  unit

val rgx_type_is_composite : Signed.long -> int
val rgxq_charpoly : gen -> gen -> Signed.long -> gen
val rgxq_inv : gen -> gen -> gen
val rgxq_minpoly : gen -> gen -> Signed.long -> gen
val rgxq_mul : gen -> gen -> gen -> gen

val rgxq_ratlift :
  gen ->
  gen ->
  Signed.long ->
  Signed.long ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  int

val rgxq_sqr : gen -> gen -> gen
val z_content : gen -> gen
val zx_content : gen -> gen
val centermod : gen -> gen -> gen
val centermod_i : gen -> gen -> gen -> gen
val centermodii : gen -> gen -> gen -> gen
val content : gen -> gen
val content0 : gen -> gen -> gen
val deg1_from_roots : gen -> Signed.long -> gen
val factor : gen -> gen
val factor0 : gen -> gen -> gen
val factorback : gen -> gen
val factorback2 : gen -> gen -> gen

val gbezout :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val gdivexact : gen -> gen -> gen

val gen_factorback :
  gen ->
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen -> gen -> gen) Ctypes_static.static_funptr ->
  (unit Ctypes_static.ptr -> gen) Ctypes_static.static_funptr ->
  gen

val ggcd : gen -> gen -> gen
val ggcd0 : gen -> gen -> gen
val ghalfgcd : gen -> gen -> gen
val ginvmod : gen -> gen -> gen
val glcm : gen -> gen -> gen
val glcm0 : gen -> gen -> gen
val newtonpoly : gen -> gen -> gen
val nfrootsq : gen -> gen
val poldisc0 : gen -> Signed.long -> gen
val polisirreducible : gen -> Signed.long
val polresultant0 : gen -> gen -> Signed.long -> Signed.long -> gen
val polsym : gen -> Signed.long -> gen
val primitive_part : gen -> gen Ctypes_static.ptr -> gen
val primpart : gen -> gen
val reduceddiscsmith : gen -> gen
val resultant2 : gen -> gen -> gen
val resultant : gen -> gen -> gen
val rnfcharpoly : gen -> gen -> gen -> Signed.long -> gen
val roots_from_deg1 : gen -> gen
val roots_to_pol : gen -> Signed.long -> gen
val roots_to_pol_r1 : gen -> Signed.long -> Signed.long -> gen
val sturmpart : gen -> gen -> gen -> Signed.long

val subresext :
  gen -> gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> gen

val sylvestermatrix : gen -> gen -> gen
val trivial_fact : unit -> gen
val gcdext0 : gen -> gen -> gen
val polresultantext0 : gen -> gen -> Signed.long -> gen
val polresultantext : gen -> gen -> gen
val prime_fact : gen -> gen
val row_q_primpart : gen -> gen
val vec_q_primpart : gen -> gen
val vecprod : gen -> gen
val zv_lcm : gen -> gen
val flx_flxy_resultant : gen -> gen -> pari_ulong -> gen
val flxx_resultant : gen -> gen -> pari_ulong -> Signed.long -> gen
val fpx_fpxy_resultant : gen -> gen -> gen -> gen
val fpx_translate : gen -> gen -> gen -> gen
val fpxqx_normalize : gen -> gen -> gen -> gen
val fpxv_fpc_mul : gen -> gen -> gen -> gen
val fpxy_fpxq_evaly : gen -> gen -> gen -> gen -> Signed.long -> gen
val fpxc_center : gen -> gen -> gen -> gen
val fpxm_center : gen -> gen -> gen -> gen
val fq_fp_mul : gen -> gen -> gen -> gen -> gen
val fq_add : gen -> gen -> gen -> gen -> gen
val fq_div : gen -> gen -> gen -> gen -> gen
val fq_halve : gen -> gen -> gen -> gen
val fq_inv : gen -> gen -> gen -> gen
val fq_invsafe : gen -> gen -> gen -> gen
val fq_mul : gen -> gen -> gen -> gen -> gen
val fq_mulu : gen -> pari_ulong -> gen -> gen -> gen
val fq_neg : gen -> gen -> gen -> gen
val fq_neg_inv : gen -> gen -> gen -> gen
val fq_pow : gen -> gen -> gen -> gen -> gen
val fq_powu : gen -> pari_ulong -> gen -> gen -> gen
val fq_sqr : gen -> gen -> gen -> gen
val fq_sqrt : gen -> gen -> gen -> gen
val fq_sqrtn : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fq_sub : gen -> gen -> gen -> gen -> gen
val fqc_fq_mul : gen -> gen -> gen -> gen -> gen
val fqc_fqv_mul : gen -> gen -> gen -> gen -> gen
val fqc_add : gen -> gen -> gen -> gen -> gen
val fqc_sub : gen -> gen -> gen -> gen -> gen
val fqv_red : gen -> gen -> gen -> gen
val fqv_roots_to_pol : gen -> gen -> gen -> Signed.long -> gen
val fqx_fq_add : gen -> gen -> gen -> gen -> gen
val fqx_fq_mul_to_monic : gen -> gen -> gen -> gen -> gen
val fqx_fq_sub : gen -> gen -> gen -> gen -> gen
val fqx_eval : gen -> gen -> gen -> gen -> gen
val fqx_translate : gen -> gen -> gen -> gen -> gen

val fqxq_matrix_pow :
  gen -> Signed.long -> Signed.long -> gen -> gen -> gen -> gen

val fqxq_powers : gen -> Signed.long -> gen -> gen -> gen -> gen
val fqxy_eval : gen -> gen -> gen -> gen -> gen -> gen
val fqxy_evalx : gen -> gen -> gen -> gen -> gen
val qx_disc : gen -> gen
val qx_gcd : gen -> gen -> gen
val qx_resultant : gen -> gen -> gen
val qxq_div : gen -> gen -> gen -> gen
val qxq_intnorm : gen -> gen -> gen
val qxq_inv : gen -> gen -> gen
val qxq_mul : gen -> gen -> gen -> gen
val qxq_norm : gen -> gen -> gen
val qxq_sqr : gen -> gen -> gen
val rg_is_fp : gen -> gen Ctypes_static.ptr -> int
val rg_is_fpxq : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> int
val rg_to_fp : gen -> gen -> gen
val rg_to_fpxq : gen -> gen -> gen -> gen
val rgc_to_fpc : gen -> gen -> gen
val rgc_to_fqc : gen -> gen -> gen -> gen
val rgm_is_fpm : gen -> gen Ctypes_static.ptr -> int
val rgm_to_flm : gen -> pari_ulong -> gen
val rgm_to_fpm : gen -> gen -> gen
val rgm_to_fqm : gen -> gen -> gen -> gen
val rgv_is_fpv : gen -> gen Ctypes_static.ptr -> int
val rgv_to_flv : gen -> pari_ulong -> gen
val rgv_to_fpv : gen -> gen -> gen
val rgx_is_fpx : gen -> gen Ctypes_static.ptr -> int
val rgx_to_fpx : gen -> gen -> gen
val rgx_is_fpxqx : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> int
val rgx_to_fpxqx : gen -> gen -> gen -> gen
val rgx_to_fqx : gen -> gen -> gen -> gen

val z_incremental_crt :
  gen Ctypes_static.ptr ->
  pari_ulong ->
  gen Ctypes_static.ptr ->
  pari_ulong ->
  int

val z_init_crt : pari_ulong -> pari_ulong -> gen

val zm_incremental_crt :
  gen Ctypes_static.ptr -> gen -> gen Ctypes_static.ptr -> pari_ulong -> int

val zm_init_crt : gen -> pari_ulong -> gen
val zx_zxy_resultant : gen -> gen -> gen
val zx_zxy_rnfequation : gen -> gen -> Signed.long Ctypes_static.ptr -> gen
val zx_disc : gen -> gen
val zx_gcd : gen -> gen -> gen
val zx_gcd_all : gen -> gen -> gen Ctypes_static.ptr -> gen

val zx_incremental_crt :
  gen Ctypes_static.ptr -> gen -> gen Ctypes_static.ptr -> pari_ulong -> int

val zx_init_crt : gen -> pari_ulong -> Signed.long -> gen
val zx_is_squarefree : gen -> int
val zx_radical : gen -> gen
val zx_resultant : gen -> gen -> gen

val zxm_incremental_crt :
  gen Ctypes_static.ptr -> gen -> gen Ctypes_static.ptr -> pari_ulong -> int

val zxm_init_crt : gen -> Signed.long -> pari_ulong -> gen
val zxq_minpoly : gen -> gen -> Signed.long -> pari_ulong -> gen
val zxq_charpoly : gen -> gen -> Signed.long -> gen
val characteristic : gen -> gen
val ffinit : gen -> Signed.long -> Signed.long -> gen
val ffnbirred : gen -> Signed.long -> gen
val ffnbirred0 : gen -> Signed.long -> Signed.long -> gen
val ffsumnbirred : gen -> Signed.long -> gen

val get_fq_field :
  unit Ctypes_static.ptr Ctypes_static.ptr ->
  gen ->
  gen ->
  bb_field Ctypes.structure Ctypes_static.ptr

val init_flxq : pari_ulong -> Signed.long -> Signed.long -> gen
val init_fq : gen -> Signed.long -> Signed.long -> gen
val nfx_disc : gen -> gen -> gen
val nfx_resultant : gen -> gen -> gen -> gen
val pol_x_powers : Signed.long -> Signed.long -> gen
val residual_characteristic : gen -> gen
val polclass : gen -> Signed.long -> Signed.long -> gen
val fp_modinv_to_j : gen -> Signed.long -> gen -> gen

val fp_polmodular_evalx :
  Signed.long -> Signed.long -> gen -> gen -> Signed.long -> int -> gen

val check_modinv : Signed.long -> unit
val disc_best_modinv : Signed.long -> Signed.long
val modinv_height_factor : Signed.long -> Signed.long
val modinv_good_disc : Signed.long -> Signed.long -> int
val modinv_good_prime : Signed.long -> Signed.long -> int
val modinv_is_weber : Signed.long -> int
val modinv_is_double_eta : Signed.long -> int

val polmodular :
  Signed.long -> Signed.long -> gen -> Signed.long -> Signed.long -> gen

val polmodular_zm : Signed.long -> Signed.long -> gen

val polmodular_zxx :
  Signed.long -> Signed.long -> Signed.long -> Signed.long -> gen

val bpsw_isprime : gen -> Signed.long
val bpsw_psp : gen -> Signed.long
val addprimes : gen -> gen
val check_ecppcert : gen -> Signed.long
val gisprime : gen -> Signed.long -> gen
val gispseudoprime : gen -> Signed.long -> gen
val gprimepi_upper_bound : gen -> gen
val gprimepi_lower_bound : gen -> gen
val isprime : gen -> Signed.long
val ispseudoprime : gen -> Signed.long -> Signed.long
val millerrabin : gen -> Signed.long -> Signed.long
val prime : Signed.long -> gen
val primecert : gen -> Signed.long -> gen
val primecert0 : gen -> Signed.long -> Signed.long -> gen
val primecertexport : gen -> Signed.long -> gen
val primecertisvalid : gen -> Signed.long
val primepi : gen -> gen
val primepi_upper_bound : float -> float
val primepi_lower_bound : float -> float
val primes : Signed.long -> gen
val primes_interval : gen -> gen -> gen
val primes_interval_zv : pari_ulong -> pari_ulong -> gen
val primes_upto_zv : pari_ulong -> gen
val primes0 : gen -> gen
val primes_zv : Signed.long -> gen
val randomprime : gen -> gen
val randomprime0 : gen -> gen -> gen
val removeprimes : gen -> gen
val uis2psp : pari_ulong -> int
val uispsp : pari_ulong -> pari_ulong -> int
val uislucaspsp : pari_ulong -> int
val uisprime : pari_ulong -> int
val uisprime_101 : pari_ulong -> int
val uisprime_661 : pari_ulong -> int
val uprime : Signed.long -> pari_ulong
val uprimepi : pari_ulong -> pari_ulong
val qfauto : gen -> gen -> gen
val qfauto0 : gen -> gen -> gen
val qfautoexport : gen -> Signed.long -> gen
val qfisom : gen -> gen -> gen -> gen -> gen
val qfisom0 : gen -> gen -> gen -> gen -> gen
val qfisominit : gen -> gen -> gen -> gen
val qfisominit0 : gen -> gen -> gen -> gen
val qforbits : gen -> gen -> gen
val qfminimize : gen -> gen
val qfparam : gen -> gen -> Signed.long -> gen
val qfsolve : gen -> gen
val z_isfundamental : gen -> Signed.long
val classno : gen -> gen
val classno2 : gen -> gen
val hclassnof_fact : gen -> gen -> gen -> gen
val hclassno : gen -> gen
val hclassno6 : gen -> gen
val isfundamental : gen -> Signed.long
val qfb_equal1 : gen -> int
val qfbclassno0 : gen -> Signed.long -> gen
val qfi_shanks : gen -> gen -> Signed.long -> gen
val qfi_log : gen -> gen -> gen -> gen
val qfi_order : gen -> gen -> gen
val quadclassnof : gen -> gen Ctypes_static.ptr -> gen
val quadclassnof_fact : gen -> gen -> gen -> gen
val quaddisc : gen -> gen
val quadregulator : gen -> Signed.long -> gen
val quadunit : gen -> gen
val quadunit0 : gen -> Signed.long -> gen
val quadunitindex : gen -> gen -> gen
val quadunitnorm : gen -> Signed.long
val sisfundamental : Signed.long -> Signed.long
val uhclassnof_fact : gen -> Signed.long -> Signed.long
val unegisfundamental : pari_ulong -> Signed.long
val unegquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val uposisfundamental : pari_ulong -> Signed.long
val uposquadclassnof : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val uquadclassnof_fact : pari_ulong -> Signed.long -> gen -> gen -> pari_ulong
val zn_quad_roots : gen -> gen -> gen -> gen
val genrand : gen -> gen
val getrand : unit -> gen
val pari_rand : unit -> pari_ulong
val randomi : gen -> gen
val randomr : Signed.long -> gen
val random_f2x : Signed.long -> Signed.long -> gen
val random_fl : pari_ulong -> pari_ulong
val random_bits : Signed.long -> Signed.long
val random_zv : Signed.long -> gen
val setrand : gen -> unit
val ellratpoints : gen -> gen -> Signed.long -> gen
val hyperellratpoints : gen -> gen -> Signed.long -> gen
val qx_complex_roots : gen -> Signed.long -> gen
val fft : gen -> gen -> gen
val fftinv : gen -> gen -> gen
val cleanroots : gen -> Signed.long -> gen
val fujiwara_bound : gen -> float
val fujiwara_bound_real : gen -> Signed.long -> float
val isrealappr : gen -> Signed.long -> int
val polgraeffe : gen -> gen
val polmod_to_embed : gen -> Signed.long -> gen
val polrootsbound : gen -> gen -> gen
val roots : gen -> Signed.long -> gen
val realroots : gen -> gen -> Signed.long -> gen
val zx_graeffe : gen -> gen
val zx_realroots_irred : gen -> Signed.long -> gen
val zx_sturm : gen -> Signed.long
val zx_sturm_irred : gen -> Signed.long
val zx_sturmpart : gen -> gen -> Signed.long
val zx_uspensky : gen -> gen -> Signed.long -> Signed.long -> gen
val factor_aurifeuille : gen -> Signed.long -> gen
val factor_aurifeuille_prime : gen -> Signed.long -> gen
val galoissubcyclo : gen -> gen -> Signed.long -> Signed.long -> gen
val polsubcyclo : Signed.long -> Signed.long -> Signed.long -> gen
val polsubcyclofast : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val znsubgroupgenerators : gen -> Signed.long -> gen
val nfsubfields : gen -> Signed.long -> gen
val nfsubfields0 : gen -> Signed.long -> Signed.long -> gen
val nfsubfieldscm : gen -> Signed.long -> gen
val nfsubfieldsmax : gen -> Signed.long -> gen
val nflist : gen -> gen -> Signed.long -> gen -> gen
val nfresolvent : gen -> Signed.long -> gen
val subgrouplist : gen -> gen -> gen

val forsubgroup :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  unit

val abmap_kernel : gen -> gen
val abmap_subgroup_image : gen -> gen -> gen
val bnrl1 : gen -> gen -> Signed.long -> Signed.long -> gen
val bnrrootnumber : gen -> gen -> Signed.long -> Signed.long -> gen
val bnrstark : gen -> gen -> Signed.long -> gen
val cyc2elts : gen -> gen
val qfbforms : gen -> gen
val quadhilbert : gen -> Signed.long -> gen
val quadray : gen -> gen -> Signed.long -> gen
val chartogenstr : char -> gen
val pari_strdup : string -> string
val pari_strndup : string -> Signed.long -> string
val stack_strcat : string -> string -> string
val stack_strdup : string -> string
val pari_strchr : gen -> gen
val strjoin : gen -> gen -> gen
val strntogenstr : string -> Signed.long -> gen
val strsplit : gen -> gen -> gen
val strtogenstr : string -> gen
val type_name : Signed.long -> string
val type0 : gen -> gen

val asympnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val asympnumraw :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  Signed.long ->
  gen ->
  Signed.long ->
  gen

val derivnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val derivnumk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen

val derivfun :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val derivfunk :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen

val forvec_init :
  forvec_t Ctypes.structure Ctypes_static.ptr -> gen -> Signed.long -> int

val forvec_next : forvec_t Ctypes.structure Ctypes_static.ptr -> gen

val laurentseries :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  Signed.long ->
  Signed.long ->
  Signed.long ->
  gen

val limitnum :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val polzag : Signed.long -> Signed.long -> gen

val prodeuler :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen

val prodinf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val prodinf1 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val solvestep :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  gen ->
  Signed.long ->
  Signed.long ->
  gen

val sumalt :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val sumalt2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val sumpos :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val sumpos2 :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val suminf :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val suminf_bitprec :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val sumdivmultexpr :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen

val zbrent :
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> gen) Ctypes_static.static_funptr ->
  gen ->
  gen ->
  Signed.long ->
  gen

val bnfisintnorm : gen -> gen -> gen
val bnfisintnormabs : gen -> gen -> gen
val ideals_by_norm : gen -> gen -> gen
val thue : gen -> gen -> gen -> gen
val thueinit : gen -> Signed.long -> Signed.long -> gen
val pi2n : Signed.long -> Signed.long -> gen
val pii2 : Signed.long -> gen
val pii2n : Signed.long -> Signed.long -> gen
val qp_exp : gen -> gen
val qp_exp_prec : gen -> Signed.long
val qp_log : gen -> gen
val qp_sqrt : gen -> gen
val qp_sqrtn : gen -> gen -> gen Ctypes_static.ptr -> gen
val zn_sqrt : gen -> gen -> gen
val zp_teichmuller : gen -> gen -> Signed.long -> gen -> gen
val agm : gen -> gen -> Signed.long -> gen
val constcatalan : Signed.long -> gen
val consteuler : Signed.long -> gen
val constlog2 : Signed.long -> gen
val constpi : Signed.long -> gen
val cxexpm1 : gen -> Signed.long -> gen
val elle : gen -> Signed.long -> gen
val ellk : gen -> Signed.long -> gen
val expir : gen -> gen
val exp1r_abs : gen -> gen
val gcos : gen -> Signed.long -> gen
val gcotan : gen -> Signed.long -> gen
val gcotanh : gen -> Signed.long -> gen
val gexp : gen -> Signed.long -> gen
val gexpm1 : gen -> Signed.long -> gen
val glog : gen -> Signed.long -> gen
val glog1p : gen -> Signed.long -> gen
val gpow : gen -> gen -> Signed.long -> gen
val gpowers : gen -> Signed.long -> gen
val gpowers0 : gen -> Signed.long -> gen -> gen
val gpowgs : gen -> Signed.long -> gen
val grootsof1 : Signed.long -> Signed.long -> gen
val gsin : gen -> Signed.long -> gen
val gsinc : gen -> Signed.long -> gen

val gsincos :
  gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> Signed.long -> unit

val gsqrpowers : gen -> Signed.long -> gen
val gsqrt : gen -> Signed.long -> gen
val gsqrtn : gen -> gen -> gen Ctypes_static.ptr -> Signed.long -> gen
val gtan : gen -> Signed.long -> gen
val logr_abs : gen -> gen
val mpcos : gen -> gen
val mpeuler : Signed.long -> gen
val mpcatalan : Signed.long -> gen
val mpsincosm1 : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit
val mpexp : gen -> gen
val mpexpm1 : gen -> gen
val mplog : gen -> gen
val mplog2 : Signed.long -> gen
val mppi : Signed.long -> gen
val mpsin : gen -> gen
val mpsincos : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit
val pow2pis : gen -> Signed.long -> gen
val powpis : gen -> Signed.long -> gen
val powcx : gen -> gen -> gen -> Signed.long -> gen
val powcx_prec : Signed.long -> gen -> Signed.long -> Signed.long
val powersr : gen -> Signed.long -> gen
val powiu : gen -> pari_ulong -> gen
val powrfrac : gen -> Signed.long -> Signed.long -> gen
val powrs : gen -> Signed.long -> gen
val powrshalf : gen -> Signed.long -> gen
val powru : gen -> pari_ulong -> gen
val powruhalf : gen -> pari_ulong -> gen
val powuu : pari_ulong -> pari_ulong -> gen
val powgi : gen -> gen -> gen
val rootsof1_cx : gen -> Signed.long -> gen
val rootsof1u_cx : pari_ulong -> Signed.long -> gen
val rootsof1q_cx : Signed.long -> Signed.long -> Signed.long -> gen
val rootsof1powinit : Signed.long -> Signed.long -> Signed.long -> gen
val rootsof1pow : gen -> Signed.long -> gen
val serchop : gen -> Signed.long -> gen
val serchop_i : gen -> Signed.long -> gen
val serchop0 : gen -> gen
val sqrtnint : gen -> Signed.long -> gen
val sqrtnr_abs : gen -> Signed.long -> gen
val teich : gen -> gen
val teichmullerinit : Signed.long -> Signed.long -> gen
val teichmuller : gen -> gen -> gen

val trans_eval :
  string ->
  (gen -> Signed.long -> gen) Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val trans_evalgen :
  string ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long -> gen)
  Ctypes_static.static_funptr ->
  gen ->
  Signed.long ->
  gen

val upowuu : pari_ulong -> pari_ulong -> pari_ulong
val upowers : pari_ulong -> Signed.long -> gen
val usqrtn : pari_ulong -> pari_ulong -> pari_ulong
val usqrt : pari_ulong -> pari_ulong
val qp_gamma : gen -> gen
val atanhuu : pari_ulong -> pari_ulong -> Signed.long -> gen
val atanhui : pari_ulong -> gen -> Signed.long -> gen
val gacosh : gen -> Signed.long -> gen
val gacos : gen -> Signed.long -> gen
val garg : gen -> Signed.long -> gen
val gasinh : gen -> Signed.long -> gen
val gasin : gen -> Signed.long -> gen
val gatan : gen -> Signed.long -> gen
val gatanh : gen -> Signed.long -> gen
val gcosh : gen -> Signed.long -> gen
val ggammah : gen -> Signed.long -> gen
val ggamma : gen -> Signed.long -> gen
val ggamma1m1 : gen -> Signed.long -> gen
val glngamma : gen -> Signed.long -> gen
val gpsi : gen -> Signed.long -> gen
val gsinh : gen -> Signed.long -> gen
val gtanh : gen -> Signed.long -> gen
val mpfactr : Signed.long -> Signed.long -> gen
val mpsinhcosh : gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit
val psi1series : Signed.long -> Signed.long -> Signed.long -> gen
val sumformal : gen -> Signed.long -> gen

val rgv_is_arithprog :
  gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> int

val besseljzero : gen -> Signed.long -> Signed.long -> gen
val besselyzero : gen -> Signed.long -> Signed.long -> gen
val constzeta : Signed.long -> Signed.long -> gen
val cxek : gen -> Signed.long -> Signed.long -> gen
val dblmodulus : gen -> float
val dilog : gen -> Signed.long -> gen
val eint1 : gen -> Signed.long -> gen
val expipir : gen -> Signed.long -> gen
val expipic : gen -> Signed.long -> gen
val expixy : gen -> gen -> Signed.long -> gen
val eta : gen -> Signed.long -> gen
val eta0 : gen -> Signed.long -> Signed.long -> gen
val gerfc : gen -> Signed.long -> gen
val gpolylog : Signed.long -> gen -> Signed.long -> gen
val gzeta : gen -> Signed.long -> gen
val hbessel1 : gen -> gen -> Signed.long -> gen
val hbessel2 : gen -> gen -> Signed.long -> gen
val hyperu : gen -> gen -> gen -> Signed.long -> gen
val ibessel : gen -> gen -> Signed.long -> gen
val incgam : gen -> gen -> Signed.long -> gen
val incgam0 : gen -> gen -> gen -> Signed.long -> gen
val incgamc : gen -> gen -> Signed.long -> gen
val jbessel : gen -> gen -> Signed.long -> gen
val jbesselh : gen -> gen -> Signed.long -> gen
val jell : gen -> Signed.long -> gen
val kbessel : gen -> gen -> Signed.long -> gen
val mpeint1 : gen -> gen -> gen
val mpveceint1 : gen -> gen -> Signed.long -> gen
val polylog0 : Signed.long -> gen -> Signed.long -> Signed.long -> gen
val sumdedekind : gen -> gen -> gen
val sumdedekind_coprime : gen -> gen -> gen
val szeta : Signed.long -> Signed.long -> gen
val theta : gen -> gen -> Signed.long -> gen
val thetanullk : gen -> Signed.long -> Signed.long -> gen
val trueeta : gen -> Signed.long -> gen
val u_sumdedekind_coprime : Signed.long -> Signed.long -> gen
val upper_to_cx : gen -> Signed.long Ctypes_static.ptr -> gen
val veceint1 : gen -> gen -> Signed.long -> gen
val vecthetanullk : gen -> Signed.long -> Signed.long -> gen
val vecthetanullk_tau : gen -> Signed.long -> Signed.long -> gen
val veczeta : gen -> gen -> Signed.long -> Signed.long -> gen
val weber0 : gen -> Signed.long -> Signed.long -> gen
val weberf : gen -> Signed.long -> gen
val weberf1 : gen -> Signed.long -> gen
val weberf2 : gen -> Signed.long -> gen
val ybessel : gen -> gen -> Signed.long -> gen
val sl2_inv_shallow : gen -> gen
val qevproj_apply : gen -> gen -> gen
val qevproj_apply_vecei : gen -> gen -> Signed.long -> gen
val qevproj_down : gen -> gen -> gen
val qevproj_init : gen -> gen
val rgx_act_gl2q : gen -> Signed.long -> gen
val rgx_act_zgl2q : gen -> Signed.long -> gen
val checkms : gen -> unit
val checkmspadic : gen -> unit
val ellpadiclambdamu : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val mfnumcusps : gen -> gen
val mfnumcusps_fact : gen -> gen
val mfnumcuspsu_fact : gen -> pari_ulong
val mfnumcuspsu : pari_ulong -> pari_ulong
val msfromcusp : gen -> gen -> gen
val msfromell : gen -> Signed.long -> gen
val msfromhecke : gen -> gen -> gen -> gen
val msdim : gen -> Signed.long
val mseval2_ooq : gen -> gen -> gen -> gen
val msgetlevel : gen -> Signed.long
val msgetsign : gen -> Signed.long
val msgetweight : gen -> Signed.long
val msatkinlehner : gen -> Signed.long -> gen -> gen
val mscuspidal : gen -> Signed.long -> gen
val mseisenstein : gen -> gen
val mseval : gen -> gen -> gen -> gen
val mshecke : gen -> Signed.long -> gen -> gen
val msinit : gen -> gen -> Signed.long -> gen
val msissymbol : gen -> gen -> gen
val mslattice : gen -> gen -> gen
val msomseval : gen -> gen -> gen -> gen

val mspadic_parse_chi :
  gen -> gen Ctypes_static.ptr -> gen Ctypes_static.ptr -> unit

val mspadic_unit_eigenvalue : gen -> Signed.long -> gen -> Signed.long -> gen
val mspadicinit : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val mspadicl : gen -> gen -> Signed.long -> gen
val mspadicmoments : gen -> gen -> Signed.long -> gen
val mspadicseries : gen -> Signed.long -> gen
val mspathgens : gen -> gen
val mspathlog : gen -> gen -> gen
val msnew : gen -> gen
val mspetersson : gen -> gen -> gen -> gen
val mspolygon : gen -> Signed.long -> gen
val msstar : gen -> gen -> gen
val msqexpansion : gen -> gen -> pari_ulong -> gen
val mssplit : gen -> gen -> Signed.long -> gen
val mstooms : gen -> gen -> gen
val mscosets0 : gen -> gen -> gen

val mscosets :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  gen

val msfarey :
  gen ->
  unit Ctypes_static.ptr ->
  (unit Ctypes_static.ptr -> gen -> Signed.long) Ctypes_static.static_funptr ->
  gen Ctypes_static.ptr ->
  gen

val msfarey0 : gen -> gen -> gen Ctypes_static.ptr -> gen
val checkfarey_i : gen -> int
val polylogmult : gen -> gen -> Signed.long -> gen
val polylogmult_interpolate : gen -> gen -> gen -> Signed.long -> gen
val zetamult : gen -> Signed.long -> gen
val zetamultdual : gen -> gen
val zetamult_interpolate : gen -> gen -> Signed.long -> gen
val zetamultall : Signed.long -> Signed.long -> Signed.long -> gen
val zetamultconvert : gen -> Signed.long -> gen
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
val abscmpiu : gen -> pari_ulong -> int
val abscmpui : pari_ulong -> gen -> int
val absequaliu : gen -> pari_ulong -> int
val absi : gen -> gen
val absi_shallow : gen -> gen
val absr : gen -> gen
val absrnz_equal1 : gen -> int
val absrnz_equal2n : gen -> int
val addii : gen -> gen -> gen
val addiiz : gen -> gen -> gen -> unit
val addir : gen -> gen -> gen
val addirz : gen -> gen -> gen -> unit
val addis : gen -> Signed.long -> gen
val addri : gen -> gen -> gen
val addriz : gen -> gen -> gen -> unit
val addrr : gen -> gen -> gen
val addrrz : gen -> gen -> gen -> unit
val addrs : gen -> Signed.long -> gen
val addsi : Signed.long -> gen -> gen
val addsiz : Signed.long -> gen -> gen -> unit
val addsrz : Signed.long -> gen -> gen -> unit
val addss : Signed.long -> Signed.long -> gen
val addssz : Signed.long -> Signed.long -> gen -> unit
val adduu : pari_ulong -> pari_ulong -> gen
val affii : gen -> gen -> unit
val affiz : gen -> gen -> unit
val affrr_fixlg : gen -> gen -> unit
val affsi : Signed.long -> gen -> unit
val affsr : Signed.long -> gen -> unit
val affsz : Signed.long -> gen -> unit
val affui : pari_ulong -> gen -> unit
val affur : pari_ulong -> gen -> unit
val cgetg : Signed.long -> Signed.long -> gen
val cgetg_block : Signed.long -> Signed.long -> gen
val cgetg_copy : gen -> Signed.long Ctypes_static.ptr -> gen
val cgeti : Signed.long -> gen
val cgetineg : Signed.long -> gen
val cgetipos : Signed.long -> gen
val cgetr : Signed.long -> gen
val cgetr_block : Signed.long -> gen
val cmpir : gen -> gen -> int
val cmpis : gen -> Signed.long -> int
val cmpiu : gen -> pari_ulong -> int
val cmpri : gen -> gen -> int
val cmprs : gen -> Signed.long -> int
val cmpsi : Signed.long -> gen -> int
val cmpsr : Signed.long -> gen -> int
val cmpss : Signed.long -> Signed.long -> int
val cmpui : pari_ulong -> gen -> int
val cmpuu : pari_ulong -> pari_ulong -> int
val divii : gen -> gen -> gen
val diviiz : gen -> gen -> gen -> unit
val divirz : gen -> gen -> gen -> unit
val divisz : gen -> Signed.long -> gen -> unit
val divriz : gen -> gen -> gen -> unit
val divrrz : gen -> gen -> gen -> unit
val divrsz : gen -> Signed.long -> gen -> unit
val divsi_rem : Signed.long -> gen -> Signed.long Ctypes_static.ptr -> gen
val divsiz : Signed.long -> gen -> gen -> unit
val divsrz : Signed.long -> gen -> gen -> unit
val divss : Signed.long -> Signed.long -> gen

val divss_rem :
  Signed.long -> Signed.long -> Signed.long Ctypes_static.ptr -> gen

val divssz : Signed.long -> Signed.long -> gen -> unit
val dvdii : gen -> gen -> int
val dvdiiz : gen -> gen -> gen -> int
val dvdis : gen -> Signed.long -> int
val dvdisz : gen -> Signed.long -> gen -> int
val dvdiu : gen -> pari_ulong -> int
val dvdiuz : gen -> pari_ulong -> gen -> int
val dvdsi : Signed.long -> gen -> int
val dvdui : pari_ulong -> gen -> int
val dvmdiiz : gen -> gen -> gen -> gen -> unit
val dvmdis : gen -> Signed.long -> gen Ctypes_static.ptr -> gen
val dvmdisz : gen -> Signed.long -> gen -> gen -> unit
val dvmdsbil : Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long
val dvmdsi : Signed.long -> gen -> gen Ctypes_static.ptr -> gen
val dvmdsiz : Signed.long -> gen -> gen -> gen -> unit
val dvmdss : Signed.long -> Signed.long -> gen Ctypes_static.ptr -> gen
val dvmdssz : Signed.long -> Signed.long -> gen -> gen -> unit
val dvmdubil : pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong
val equalis : gen -> Signed.long -> int
val equalsi : Signed.long -> gen -> int
val equalui : pari_ulong -> gen -> int
val equaliu : gen -> pari_ulong -> int
val absequalui : pari_ulong -> gen -> int
val ceildivuu : pari_ulong -> pari_ulong -> pari_ulong
val evalexpo : Signed.long -> Signed.long
val evallg : Signed.long -> Signed.long
val evalprecp : Signed.long -> Signed.long
val evalvalp : Signed.long -> Signed.long
val evalvalser : Signed.long -> Signed.long
val expi : gen -> Signed.long
val expu : pari_ulong -> Signed.long
val fixlg : gen -> Signed.long -> unit
val fractor : gen -> Signed.long -> gen
val gc_bool : pari_ulong -> int -> int
val gc_const : pari_ulong -> gen -> gen
val gc_double : pari_ulong -> float -> float
val gc_int : pari_ulong -> int -> int
val gc_long : pari_ulong -> Signed.long -> Signed.long
val gc_stoi : pari_ulong -> Signed.long -> gen
val gc_ulong : pari_ulong -> pari_ulong -> pari_ulong
val gc_utoi : pari_ulong -> pari_ulong -> gen
val gc_utoipos : pari_ulong -> pari_ulong -> gen
val gc_null : pari_ulong -> gen
val icopy : gen -> gen
val icopyspec : gen -> Signed.long -> gen
val int_bit : gen -> Signed.long -> pari_ulong
val itor : gen -> Signed.long -> gen
val itos : gen -> Signed.long
val itos_or_0 : gen -> Signed.long
val itou : gen -> pari_ulong
val itou_or_0 : gen -> pari_ulong
val leafcopy : gen -> gen
val maxdd : float -> float -> float
val maxss : Signed.long -> Signed.long -> Signed.long
val maxuu : pari_ulong -> pari_ulong -> Signed.long
val mindd : float -> float -> float
val minss : Signed.long -> Signed.long -> Signed.long
val minuu : pari_ulong -> pari_ulong -> Signed.long
val mod16 : gen -> Signed.long
val mod2 : gen -> Signed.long
val mod2bil : gen -> pari_ulong
val mod32 : gen -> Signed.long
val mod4 : gen -> Signed.long
val mod64 : gen -> Signed.long
val mod8 : gen -> Signed.long
val modis : gen -> Signed.long -> gen
val modisz : gen -> Signed.long -> gen -> unit
val modsi : Signed.long -> gen -> gen
val modsiz : Signed.long -> gen -> gen -> unit
val modss : Signed.long -> Signed.long -> gen
val modssz : Signed.long -> Signed.long -> gen -> unit
val mpabs : gen -> gen
val mpabs_shallow : gen -> gen
val mpadd : gen -> gen -> gen
val mpaddz : gen -> gen -> gen -> unit
val mpaff : gen -> gen -> unit
val mpceil : gen -> gen
val mpcmp : gen -> gen -> int
val mpcopy : gen -> gen
val mpdiv : gen -> gen -> gen
val mpexpo : gen -> Signed.long
val mpfloor : gen -> gen
val mpmul : gen -> gen -> gen
val mpmulz : gen -> gen -> gen -> unit
val mpneg : gen -> gen
val mpodd : gen -> int
val mpround : gen -> gen
val mpshift : gen -> Signed.long -> gen
val mpsqr : gen -> gen
val mpsub : gen -> gen -> gen
val mpsubz : gen -> gen -> gen -> unit
val mptrunc : gen -> gen
val muliiz : gen -> gen -> gen -> unit
val mulirz : gen -> gen -> gen -> unit
val mulis : gen -> Signed.long -> gen
val muliu : gen -> pari_ulong -> gen
val mulri : gen -> gen -> gen
val mulriz : gen -> gen -> gen -> unit
val mulrrz : gen -> gen -> gen -> unit
val mulrs : gen -> Signed.long -> gen
val mulru : gen -> pari_ulong -> gen
val mulsiz : Signed.long -> gen -> gen -> unit
val mulsrz : Signed.long -> gen -> gen -> unit
val mulssz : Signed.long -> Signed.long -> gen -> unit
val negi : gen -> gen
val negr : gen -> gen
val new_chunk : int -> gen
val rcopy : gen -> gen
val rdivii : gen -> gen -> Signed.long -> gen
val rdiviiz : gen -> gen -> gen -> unit
val rdivis : gen -> Signed.long -> Signed.long -> gen
val rdivsi : Signed.long -> gen -> Signed.long -> gen
val rdivss : Signed.long -> Signed.long -> Signed.long -> gen
val real2n : Signed.long -> Signed.long -> gen
val real_m2n : Signed.long -> Signed.long -> gen
val real_0 : Signed.long -> gen
val real_0_bit : Signed.long -> gen
val real_1 : Signed.long -> gen
val real_1_bit : Signed.long -> gen
val real_m1 : Signed.long -> gen
val remii : gen -> gen -> gen
val remiiz : gen -> gen -> gen -> unit
val remis : gen -> Signed.long -> gen
val remisz : gen -> Signed.long -> gen -> unit

val remlll_pre :
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong ->
  pari_ulong

val remsi : Signed.long -> gen -> gen
val remsiz : Signed.long -> gen -> gen -> unit
val remss : Signed.long -> Signed.long -> gen
val remssz : Signed.long -> Signed.long -> gen -> unit
val rtor : gen -> Signed.long -> gen
val sdivsi : Signed.long -> gen -> Signed.long

val sdivsi_rem :
  Signed.long -> gen -> Signed.long Ctypes_static.ptr -> Signed.long

val sdivss_rem :
  Signed.long -> Signed.long -> Signed.long Ctypes_static.ptr -> Signed.long

val get_avma : unit -> pari_ulong
val set_avma : pari_ulong -> unit

val uabsdiviu_rem :
  gen -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val udivuu_rem :
  pari_ulong -> pari_ulong -> pari_ulong Ctypes_static.ptr -> pari_ulong

val umodi2n : gen -> Signed.long -> pari_ulong
val setabssign : gen -> unit

val shift_left :
  gen -> gen -> Signed.long -> Signed.long -> pari_ulong -> pari_ulong -> unit

val shift_right :
  gen -> gen -> Signed.long -> Signed.long -> pari_ulong -> pari_ulong -> unit

val shiftl : pari_ulong -> pari_ulong -> pari_ulong
val shiftlr : pari_ulong -> pari_ulong -> pari_ulong
val shiftr : gen -> Signed.long -> gen
val shiftr_inplace : gen -> Signed.long -> unit
val smodis : gen -> Signed.long -> Signed.long
val smodss : Signed.long -> Signed.long -> Signed.long
val stackdummy : pari_ulong -> pari_ulong -> unit
val stack_malloc : int -> string
val stack_malloc_align : int -> Signed.long -> string
val stack_calloc : int -> string
val stack_calloc_align : int -> Signed.long -> string
val stoi : Signed.long -> gen
val stor : Signed.long -> Signed.long -> gen
val subii : gen -> gen -> gen
val subiiz : gen -> gen -> gen -> unit
val subir : gen -> gen -> gen
val subirz : gen -> gen -> gen -> unit
val subis : gen -> Signed.long -> gen
val subisz : gen -> Signed.long -> gen -> unit
val subri : gen -> gen -> gen
val subriz : gen -> gen -> gen -> unit
val subrr : gen -> gen -> gen
val subrrz : gen -> gen -> gen -> unit
val subrs : gen -> Signed.long -> gen
val subrsz : gen -> Signed.long -> gen -> unit
val subsi : Signed.long -> gen -> gen
val subsiz : Signed.long -> gen -> gen -> unit
val subsrz : Signed.long -> gen -> gen -> unit
val subss : Signed.long -> Signed.long -> gen
val subssz : Signed.long -> Signed.long -> gen -> unit
val subuu : pari_ulong -> pari_ulong -> gen
val togglesign : gen -> unit
val togglesign_safe : gen Ctypes_static.ptr -> unit
val affectsign : gen -> gen -> unit
val affectsign_safe : gen -> gen Ctypes_static.ptr -> unit
val truedivii : gen -> gen -> gen
val truedivis : gen -> Signed.long -> gen
val truedivsi : Signed.long -> gen -> gen

val uabsdivui_rem :
  pari_ulong -> gen -> pari_ulong Ctypes_static.ptr -> pari_ulong

val umodsu : Signed.long -> pari_ulong -> pari_ulong
val umodui : pari_ulong -> gen -> pari_ulong
val ugcdiu : gen -> pari_ulong -> pari_ulong
val ugcdui : pari_ulong -> gen -> pari_ulong
val umuluu_le : pari_ulong -> pari_ulong -> pari_ulong -> pari_ulong
val umuluu_or_0 : pari_ulong -> pari_ulong -> pari_ulong
val utoi : pari_ulong -> gen
val utoineg : pari_ulong -> gen
val utoipos : pari_ulong -> gen
val utor : pari_ulong -> Signed.long -> gen
val uutoi : pari_ulong -> pari_ulong -> gen
val uutoineg : pari_ulong -> pari_ulong -> gen
val vali : gen -> Signed.long
val varncmp : Signed.long -> Signed.long -> int
val varnmax : Signed.long -> Signed.long -> Signed.long
val varnmin : Signed.long -> Signed.long -> Signed.long
val pari_err_component : string -> string -> gen -> gen -> unit
val pari_err_dim : string -> unit
val pari_err_domain : string -> string -> string -> gen -> gen -> unit
val pari_err_file : string -> string -> unit
val pari_err_filedesc : string -> Signed.long -> unit
val pari_err_flag : string -> unit
val pari_err_impl : string -> unit
val pari_err_inv : string -> gen -> unit
val pari_err_irredpol : string -> gen -> unit
val pari_err_maxprime : pari_ulong -> unit
val pari_err_modulus : string -> gen -> gen -> unit
val pari_err_op : string -> gen -> gen -> unit
val pari_err_overflow : string -> unit
val pari_err_package : string -> unit
val pari_err_prec : string -> unit
val pari_err_prime : string -> gen -> unit
val pari_err_priority : string -> gen -> string -> Signed.long -> unit
val pari_err_sqrtn : string -> gen -> unit
val pari_err_type : string -> gen -> unit
val pari_err_type2 : string -> gen -> gen -> unit
val pari_err_var : string -> gen -> gen -> unit
val pari_err_roots0 : string -> unit
val mkintmod : gen -> gen -> gen
val mkintmodu : pari_ulong -> pari_ulong -> gen
val mkpolmod : gen -> gen -> gen
val mkfrac : gen -> gen -> gen
val mkfracss : Signed.long -> Signed.long -> gen

val qtoss :
  gen -> Signed.long Ctypes_static.ptr -> Signed.long Ctypes_static.ptr -> unit

val sstoq : Signed.long -> Signed.long -> gen
val uutoq : pari_ulong -> pari_ulong -> gen
val mkfraccopy : gen -> gen -> gen
val mkrfrac : gen -> gen -> gen
val mkrfraccopy : gen -> gen -> gen
val mkcomplex : gen -> gen -> gen
val gen_i : unit -> gen
val cgetc : Signed.long -> gen
val mkquad : gen -> gen -> gen -> gen
val mkvecsmall : Signed.long -> gen
val mkvecsmall2 : Signed.long -> Signed.long -> gen
val mkvecsmall3 : Signed.long -> Signed.long -> Signed.long -> gen

val mkvecsmall4 :
  Signed.long -> Signed.long -> Signed.long -> Signed.long -> gen

val mkvecsmall5 :
  Signed.long -> Signed.long -> Signed.long -> Signed.long -> Signed.long -> gen

val mkqfb : gen -> gen -> gen -> gen -> gen
val mkvec : gen -> gen
val mkvec2 : gen -> gen -> gen
val mkvec3 : gen -> gen -> gen -> gen
val mkvec4 : gen -> gen -> gen -> gen -> gen
val mkvec5 : gen -> gen -> gen -> gen -> gen -> gen
val mkvecs : Signed.long -> gen
val mkvec2s : Signed.long -> Signed.long -> gen
val mkvec3s : Signed.long -> Signed.long -> Signed.long -> gen
val mkvec4s : Signed.long -> Signed.long -> Signed.long -> Signed.long -> gen
val mkveccopy : gen -> gen
val mkvec2copy : gen -> gen -> gen
val mkcol : gen -> gen
val mkcol2 : gen -> gen -> gen
val mkcol3 : gen -> gen -> gen -> gen
val mkcol4 : gen -> gen -> gen -> gen -> gen
val mkcol5 : gen -> gen -> gen -> gen -> gen -> gen
val mkcol6 : gen -> gen -> gen -> gen -> gen -> gen -> gen
val mkcols : Signed.long -> gen
val mkcol2s : Signed.long -> Signed.long -> gen
val mkcol3s : Signed.long -> Signed.long -> Signed.long -> gen
val mkcol4s : Signed.long -> Signed.long -> Signed.long -> Signed.long -> gen
val mkcolcopy : gen -> gen
val mkmat : gen -> gen
val mkmat2 : gen -> gen -> gen
val mkmat3 : gen -> gen -> gen -> gen
val mkmat4 : gen -> gen -> gen -> gen -> gen
val mkmat5 : gen -> gen -> gen -> gen -> gen -> gen
val mkmatcopy : gen -> gen
val mkerr : Signed.long -> gen
val mkoo : unit -> gen
val mkmoo : unit -> gen
val inf_get_sign : gen -> Signed.long
val mkmat22s : Signed.long -> Signed.long -> Signed.long -> Signed.long -> gen
val mkmat22 : gen -> gen -> gen -> gen -> gen
val pol_x : Signed.long -> gen
val pol_xn : Signed.long -> Signed.long -> gen
val pol_xnall : Signed.long -> Signed.long -> gen
val polxn_flx : Signed.long -> Signed.long -> gen
val pol_1 : Signed.long -> gen
val pol_0 : Signed.long -> gen
val const_vec : Signed.long -> gen -> gen
val const_col : Signed.long -> gen -> gen
val const_vecsmall : Signed.long -> Signed.long -> gen
val zeropadic : gen -> Signed.long -> gen
val zeropadic_shallow : gen -> Signed.long -> gen
val zeroser : Signed.long -> Signed.long -> gen
val ser_isexactzero : gen -> int
val zeropol : Signed.long -> gen
val zerocol : Signed.long -> gen
val zerovec : Signed.long -> gen
val zeromat : Signed.long -> Signed.long -> gen
val zero_flx : Signed.long -> gen
val zero_flv : Signed.long -> gen
val zero_flm : Signed.long -> Signed.long -> gen
val zero_flm_copy : Signed.long -> Signed.long -> gen
val zero_f2v : Signed.long -> gen
val zero_f2m : Signed.long -> Signed.long -> gen
val zero_f2m_copy : Signed.long -> Signed.long -> gen
val zeromatcopy : Signed.long -> Signed.long -> gen
val zerovec_block : Signed.long -> gen
val col_ei : Signed.long -> Signed.long -> gen
val vec_ei : Signed.long -> Signed.long -> gen
val f2v_ei : Signed.long -> Signed.long -> gen
val vecsmall_ei : Signed.long -> Signed.long -> gen
val rg_col_ei : gen -> Signed.long -> Signed.long -> gen
val shallowcopy : gen -> gen
val vectrunc_init : Signed.long -> gen
val coltrunc_init : Signed.long -> gen
val lg_increase : gen -> unit
val vectrunc_append : gen -> gen -> unit
val vectrunc_append_batch : gen -> gen -> unit
val vecsmalltrunc_init : Signed.long -> gen
val vecsmalltrunc_append : gen -> Signed.long -> unit
val hash_str : string -> pari_ulong
val hash_str_len : string -> Signed.long -> pari_ulong
val vec_shorten : gen -> Signed.long -> gen
val vec_lengthen : gen -> Signed.long -> gen
val vec_append : gen -> gen -> gen
val vec_prepend : gen -> gen -> gen
val vec_setconst : gen -> gen -> gen
val vecsmall_shorten : gen -> Signed.long -> gen
val vecsmall_lengthen : gen -> Signed.long -> gen
val vec_to_vecsmall : gen -> gen
val vecsmall_to_vec : gen -> gen
val vecsmall_to_vec_inplace : gen -> gen
val vecsmall_to_col : gen -> gen
val vecsmall_lexcmp : gen -> gen -> int
val vecsmall_prefixcmp : gen -> gen -> int
val vecsmall_prepend : gen -> Signed.long -> gen
val vecsmall_append : gen -> Signed.long -> gen
val vecsmall_concat : gen -> gen -> gen
val vecsmall_coincidence : gen -> gen -> Signed.long
val vecsmall_isin : gen -> Signed.long -> Signed.long
val vecsmall_pack : gen -> Signed.long -> Signed.long -> Signed.long
val vecsmall_indexmax : gen -> Signed.long
val vecsmall_max : gen -> Signed.long
val vecsmall_indexmin : gen -> Signed.long
val vecsmall_min : gen -> Signed.long
val zv_isscalar : gen -> int
val qv_isscalar : gen -> int
val rgv_isscalar : gen -> int
val rgx_isscalar : gen -> int
val rgx_equal_var : gen -> gen -> Signed.long
val rgx_to_rgv : gen -> Signed.long -> gen
val rgx_is_rational : gen -> int
val rgx_is_zx : gen -> int
val rgx_is_qx : gen -> int
val rgx_is_monomial : gen -> int
val rgv_is_zv : gen -> int
val rgv_is_qv : gen -> int
val rgv_isin_i : gen -> gen -> Signed.long -> Signed.long
val rgv_isin : gen -> gen -> Signed.long
val vecslice : gen -> Signed.long -> Signed.long -> gen
val vecslicepermute : gen -> gen -> Signed.long -> Signed.long -> gen
val rowslicepermute : gen -> gen -> Signed.long -> Signed.long -> gen
val rowslice : gen -> Signed.long -> Signed.long -> gen

val matslice :
  gen -> Signed.long -> Signed.long -> Signed.long -> Signed.long -> gen

val rowsplice : gen -> Signed.long -> gen
val vecsplice : gen -> Signed.long -> gen
val rgm_minor : gen -> Signed.long -> Signed.long -> gen
val row : gen -> Signed.long -> gen
val flm_row : gen -> Signed.long -> gen
val rowcopy : gen -> Signed.long -> gen
val row_i : gen -> Signed.long -> Signed.long -> Signed.long -> gen
val vecreverse : gen -> gen
val vecsmall_reverse : gen -> gen
val vecreverse_inplace : gen -> unit
val vecsmallpermute : gen -> gen -> gen
val vecpermute : gen -> gen -> gen
val rowpermute : gen -> gen -> gen
val identity_zv : Signed.long -> gen
val identity_perm : Signed.long -> gen
val cyclic_perm : Signed.long -> Signed.long -> gen
val perm_mul : gen -> gen -> gen
val perm_sqr : gen -> gen
val perm_inv : gen -> gen
val perm_conj : gen -> gen -> gen
val pari_free : unit Ctypes_static.ptr -> unit
val pari_malloc : int -> unit Ctypes_static.ptr
val pari_realloc : unit Ctypes_static.ptr -> int -> unit Ctypes_static.ptr
val pari_realloc_ip : unit Ctypes_static.ptr Ctypes_static.ptr -> int -> unit
val pari_calloc : int -> unit Ctypes_static.ptr
val cgetalloc : int -> Signed.long -> gen
val icopy_avma : gen -> pari_ulong -> gen
val leafcopy_avma : gen -> pari_ulong -> gen
val gerepileuptoleaf : pari_ulong -> gen -> gen
val gerepileuptoint : pari_ulong -> gen -> gen
val gerepileupto : pari_ulong -> gen -> gen
val gerepilecopy : pari_ulong -> gen -> gen
val gunclonenull : gen -> unit
val gunclonenull_deep : gen -> unit

val gerepilemany :
  pari_ulong -> gen Ctypes_static.ptr Ctypes_static.ptr -> int -> unit

val gerepileall : pari_ulong -> int -> unit
val gc_all : pari_ulong -> int -> gen
val gerepilecoeffs : pari_ulong -> gen -> int -> unit
val bin_copy : genbin Ctypes.structure Ctypes_static.ptr -> gen
val genbinbase : genbin Ctypes.structure Ctypes_static.ptr -> gen
val cgiv : gen -> unit
val killblock : gen -> unit
val is_universal_constant : gen -> int
val cxcompotor : gen -> Signed.long -> gen
val cxtofp : gen -> Signed.long -> gen
val cxtoreal : gen -> gen
val gtodouble : gen -> float
val gisdouble : gen -> float Ctypes_static.ptr -> int
val gtos : gen -> Signed.long
val gtou : gen -> pari_ulong
val absfrac : gen -> gen
val absfrac_shallow : gen -> gen
val q_abs : gen -> gen
val q_abs_shallow : gen -> gen
val r_abs_shallow : gen -> gen
val r_abs : gen -> gen
val gtofp : gen -> Signed.long -> gen
val gtomp : gen -> Signed.long -> gen
val rgx_gtofp : gen -> Signed.long -> gen
val rgc_gtofp : gen -> Signed.long -> gen
val rgv_gtofp : gen -> Signed.long -> gen
val rgm_gtofp : gen -> Signed.long -> gen
val rgc_gtomp : gen -> Signed.long -> gen
val rgm_gtomp : gen -> Signed.long -> gen
val rgx_fpnorml2 : gen -> Signed.long -> gen
val rgc_fpnorml2 : gen -> Signed.long -> gen
val rgm_fpnorml2 : gen -> Signed.long -> gen
val affgr : gen -> gen -> unit
val affc_fixlg : gen -> gen -> gen
val trunc_safe : gen -> gen
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
val bit_prec : gen -> Signed.long
val bit_accuracy : Signed.long -> Signed.long
val prec2ndec : Signed.long -> Signed.long
val nbits2ndec : Signed.long -> Signed.long
val precdbl : Signed.long -> Signed.long
val divsbil : Signed.long -> Signed.long
val remsbil : Signed.long -> Signed.long
val fp_red : gen -> gen -> gen
val fp_add : gen -> gen -> gen -> gen
val fp_sub : gen -> gen -> gen -> gen
val fp_neg : gen -> gen -> gen
val fp_halve : gen -> gen -> gen
val fp_center : gen -> gen -> gen -> gen
val fp_center_i : gen -> gen -> gen -> gen
val fp_addmul : gen -> gen -> gen -> gen -> gen
val fp_mul : gen -> gen -> gen -> gen
val fp_sqr : gen -> gen -> gen
val fp_mulu : gen -> pari_ulong -> gen -> gen
val fp_muls : gen -> Signed.long -> gen -> gen
val fp_inv : gen -> gen -> gen
val fp_invsafe : gen -> gen -> gen
val fp_div : gen -> gen -> gen -> gen
val fp_divu : gen -> pari_ulong -> gen -> gen
val flx_mulu : gen -> pari_ulong -> pari_ulong -> gen
val get_f2x_mod : gen -> gen
val get_f2x_var : gen -> Signed.long
val get_f2x_degree : gen -> Signed.long
val get_f2xqx_mod : gen -> gen
val get_f2xqx_var : gen -> Signed.long
val get_f2xqx_degree : gen -> Signed.long
val get_flx_mod : gen -> gen
val get_flx_var : gen -> Signed.long
val get_flx_degree : gen -> Signed.long
val get_flxqx_mod : gen -> gen
val get_flxqx_var : gen -> Signed.long
val get_flxqx_degree : gen -> Signed.long
val get_fpx_mod : gen -> gen
val get_fpx_var : gen -> Signed.long
val get_fpx_degree : gen -> Signed.long
val get_fpxqx_mod : gen -> gen
val get_fpxqx_var : gen -> Signed.long
val get_fpxqx_degree : gen -> Signed.long
val submulii : gen -> gen -> gen -> gen
val mulsubii : gen -> gen -> gen -> gen
val submuliu : gen -> gen -> pari_ulong -> gen
val addmuliu : gen -> gen -> pari_ulong -> gen
val submuliu_inplace : gen -> gen -> pari_ulong -> gen
val addmuliu_inplace : gen -> gen -> pari_ulong -> gen
val lincombii : gen -> gen -> gen -> gen -> gen
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
val qfb_is_qfi : gen -> int
val sqrtr : gen -> gen
val cbrtr_abs : gen -> gen
val cbrtr : gen -> gen
val sqrtnr : gen -> Signed.long -> gen
val logint : gen -> gen -> Signed.long
val ulogint : pari_ulong -> pari_ulong -> pari_ulong
val ismpzero : gen -> int
val isintzero : gen -> int
val isint1 : gen -> int
val isintm1 : gen -> int
val equali1 : gen -> int
val equalim1 : gen -> int
val is_pm1 : gen -> int
val is_bigint : gen -> int
val odd : Signed.long -> int
val both_odd : Signed.long -> Signed.long -> int
val isonstack : gen -> int
val dbllog2r : gen -> float
val mul_content : gen -> gen -> gen
val inv_content : gen -> gen
val div_content : gen -> gen -> gen
val mul_denom : gen -> gen -> gen
val constant_coeff : gen -> gen
val leading_coeff : gen -> gen
val flx_lead : gen -> pari_ulong
val flx_constant : gen -> pari_ulong
val degpol : gen -> Signed.long
val lgpol : gen -> Signed.long
val lgcols : gen -> Signed.long
val nbrows : gen -> Signed.long
val truecoef : gen -> Signed.long -> gen
val zxq_mul : gen -> gen -> gen -> gen
val zxq_sqr : gen -> gen -> gen
val rgx_copy : gen -> gen
val rgx_coeff : gen -> Signed.long -> gen
val rgx_renormalize : gen -> gen
val rgx_div : gen -> gen -> gen
val rgxqx_div : gen -> gen -> gen -> gen
val rgxqx_rem : gen -> gen -> gen -> gen
val fpx_div : gen -> gen -> gen -> gen
val flx_div : gen -> gen -> pari_ulong -> gen
val flx_div_pre : gen -> gen -> pari_ulong -> pari_ulong -> gen
val f2x_div : gen -> gen -> gen
val fpv_fpc_mul : gen -> gen -> gen -> gen
val pol0_flx : Signed.long -> gen
val pol1_flx : Signed.long -> gen
val polx_flx : Signed.long -> gen
val zero_zx : Signed.long -> gen
val polx_zx : Signed.long -> gen
val zx_shift : gen -> Signed.long -> gen
val zero_f2x : Signed.long -> gen
val pol0_f2x : Signed.long -> gen
val pol1_f2x : Signed.long -> gen
val polx_f2x : Signed.long -> gen
val f2x_equal1 : gen -> int
val f2x_equal : gen -> gen -> int
val f2x_copy : gen -> gen
val f2v_copy : gen -> gen
val flv_copy : gen -> gen
val flx_copy : gen -> gen
val vecsmall_copy : gen -> gen
val flx_equal1 : gen -> int
val zx_equal1 : gen -> int
val zx_is_monic : gen -> int
val zx_renormalize : gen -> Signed.long -> gen
val fpx_renormalize : gen -> Signed.long -> gen
val fpxx_renormalize : gen -> Signed.long -> gen
val fpxqx_renormalize : gen -> Signed.long -> gen
val f2x_renormalize : gen -> Signed.long -> gen
val f2xx_shift : gen -> Signed.long -> Signed.long -> gen
val f2v_to_f2x : gen -> Signed.long -> gen
val sturm : gen -> Signed.long
val gval : gen -> Signed.long -> Signed.long
val rgx_shift_inplace_init : Signed.long -> unit
val rgx_shift_inplace : gen -> Signed.long -> gen
val zc_to_zc : gen -> gen
val zv_to_zv : gen -> gen
val zx_to_zv : gen -> Signed.long -> gen
val zv_to_zx : gen -> Signed.long -> gen
val zm_to_zxv : gen -> Signed.long -> gen
val zero_zm : Signed.long -> Signed.long -> gen
val zero_zv : Signed.long -> gen
val zm_transpose : gen -> gen
val zm_copy : gen -> gen
val zv_copy : gen -> gen
val zm_row : gen -> Signed.long -> gen
val zc_hnfrem : gen -> gen -> gen
val zm_hnfrem : gen -> gen -> gen
val zm_lll : gen -> float -> Signed.long -> gen

val rgm_dimensions :
  gen -> Signed.long Ctypes_static.ptr -> Signed.long Ctypes_static.ptr -> unit

val rgm_shallowcopy : gen -> gen
val f2m_copy : gen -> gen
val f3m_copy : gen -> gen
val flm_copy : gen -> gen
val zv_dvd : gen -> gen -> int
val zm_zv_mod : gen -> gen -> gen
val zv_zv_mod : gen -> gen -> gen
val vecmodii : gen -> gen -> gen
val vecmoduu : gen -> gen -> gen
val fq_red : gen -> gen -> gen -> gen
val fq_to_fpxq : gen -> gen -> gen -> gen
val rg_to_fq : gen -> gen -> gen -> gen
val gener_fq_local : gen -> gen -> gen -> gen
val random_fq : gen -> gen -> gen
val fpxqx_div : gen -> gen -> gen -> gen -> gen
val flxqx_div : gen -> gen -> gen -> pari_ulong -> gen
val flxqx_div_pre : gen -> gen -> gen -> pari_ulong -> pari_ulong -> gen
val f2xqx_div : gen -> gen -> gen -> gen
val fpxy_fq_evaly : gen -> gen -> gen -> gen -> Signed.long -> gen
val fqx_red : gen -> gen -> gen -> gen
val fqx_add : gen -> gen -> gen -> gen -> gen
val fqx_neg : gen -> gen -> gen -> gen
val fqx_sub : gen -> gen -> gen -> gen -> gen
val fqx_fp_mul : gen -> gen -> gen -> gen -> gen
val fqx_fq_mul : gen -> gen -> gen -> gen -> gen
val fqx_mul : gen -> gen -> gen -> gen -> gen
val fqx_mulu : gen -> pari_ulong -> gen -> gen -> gen
val fqx_sqr : gen -> gen -> gen -> gen
val fqx_powu : gen -> pari_ulong -> gen -> gen -> gen
val fqx_halve : gen -> gen -> gen -> gen
val fqx_div : gen -> gen -> gen -> gen -> gen
val fqx_get_red : gen -> gen -> gen -> gen
val fqx_rem : gen -> gen -> gen -> gen -> gen
val fqx_divrem : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fqx_div_by_x_x : gen -> gen -> gen -> gen -> gen Ctypes_static.ptr -> gen
val fqx_halfgcd : gen -> gen -> gen -> gen -> gen
val fqx_gcd : gen -> gen -> gen -> gen -> gen

val fqx_extgcd :
  gen ->
  gen ->
  gen ->
  gen ->
  gen Ctypes_static.ptr ->
  gen Ctypes_static.ptr ->
  gen

val fqx_normalize : gen -> gen -> gen -> gen
val fqx_deriv : gen -> gen -> gen -> gen
val fqx_integ : gen -> gen -> gen -> gen
val fqx_factor : gen -> gen -> gen -> gen
val fqx_factor_squarefree : gen -> gen -> gen -> gen
val fqx_ddf : gen -> gen -> gen -> gen
val fqx_degfact : gen -> gen -> gen -> gen
val fqx_roots : gen -> gen -> gen -> gen
val fqx_to_mod : gen -> gen -> gen -> gen
val fqxq_add : gen -> gen -> gen -> gen -> gen -> gen
val fqxq_sub : gen -> gen -> gen -> gen -> gen -> gen
val fqxq_div : gen -> gen -> gen -> gen -> gen -> gen
val fqxq_inv : gen -> gen -> gen -> gen -> gen
val fqxq_invsafe : gen -> gen -> gen -> gen -> gen
val fqxq_mul : gen -> gen -> gen -> gen -> gen -> gen
val fqxq_sqr : gen -> gen -> gen -> gen -> gen
val fqxq_pow : gen -> gen -> gen -> gen -> gen -> gen
val fqxn_expint : gen -> Signed.long -> gen -> gen -> gen
val fqxn_exp : gen -> Signed.long -> gen -> gen -> gen
val fqxn_inv : gen -> Signed.long -> gen -> gen -> gen
val fqxn_mul : gen -> gen -> Signed.long -> gen -> gen -> gen
val fqxn_sqr : gen -> Signed.long -> gen -> gen -> gen
val fpxq_add : gen -> gen -> gen -> gen -> gen
val fpxq_sub : gen -> gen -> gen -> gen -> gen
val flxq_add : gen -> gen -> gen -> pari_ulong -> gen
val flxq_sub : gen -> gen -> gen -> pari_ulong -> gen
val f2x_coeff : gen -> Signed.long -> pari_ulong
val f2x_clear : gen -> Signed.long -> unit
val f2x_set : gen -> Signed.long -> unit
val f2x_flip : gen -> Signed.long -> unit
val f2v_coeff : gen -> Signed.long -> pari_ulong
val f2v_clear : gen -> Signed.long -> unit
val f2v_set : gen -> Signed.long -> unit
val f2v_flip : gen -> Signed.long -> unit
val f2m_coeff : gen -> Signed.long -> Signed.long -> pari_ulong
val f2m_clear : gen -> Signed.long -> Signed.long -> unit
val f2m_set : gen -> Signed.long -> Signed.long -> unit
val f2m_flip : gen -> Signed.long -> Signed.long -> unit
val f3m_coeff : gen -> Signed.long -> Signed.long -> pari_ulong
val f3m_set : gen -> Signed.long -> Signed.long -> pari_ulong -> unit
val matpascal : Signed.long -> gen
val z_issquare : gen -> Signed.long
val z_ispower : gen -> pari_ulong -> Signed.long
val sqrti : gen -> gen
val gaddgs : gen -> Signed.long -> gen
val gcmpgs : gen -> Signed.long -> int
val gequalgs : gen -> Signed.long -> int
val gmaxsg : Signed.long -> gen -> gen
val gminsg : Signed.long -> gen -> gen
val gmulgs : gen -> Signed.long -> gen
val gmulgu : gen -> pari_ulong -> gen
val gsubgs : gen -> Signed.long -> gen
val gdivsg : Signed.long -> gen -> gen
val gmax_shallow : gen -> gen -> gen
val gmin_shallow : gen -> gen -> gen
val cxnorm : gen -> gen
val quadnorm : gen -> gen
val quad_disc : gen -> gen
val qfb_disc3 : gen -> gen -> gen -> gen
val qfb_disc : gen -> gen
val sqrfrac : gen -> gen
val normalize_frac : gen -> unit
val powii : gen -> gen -> gen
val powis : Signed.long -> gen
val mpexpz : gen -> gen -> unit
val mplogz : gen -> gen -> unit
val mpcosz : gen -> gen -> unit
val mpsinz : gen -> gen -> unit
val gnegz : gen -> gen -> unit
val gabsz : gen -> Signed.long -> gen -> unit
val gaddz : gen -> gen -> gen -> unit
val gsubz : gen -> gen -> gen -> unit
val gmulz : gen -> gen -> gen -> unit
val gdivz : gen -> gen -> gen -> unit
val gdiventz : gen -> gen -> gen -> unit
val gmodz : gen -> gen -> gen -> unit
val gmul2nz : gen -> Signed.long -> gen -> unit
val gshiftz : gen -> Signed.long -> gen -> unit
val ell_get_a1 : gen -> gen
val ell_get_a2 : gen -> gen
val ell_get_a3 : gen -> gen
val ell_get_a4 : gen -> gen
val ell_get_a6 : gen -> gen
val ell_get_b2 : gen -> gen
val ell_get_b4 : gen -> gen
val ell_get_b6 : gen -> gen
val ell_get_b8 : gen -> gen
val ell_get_c4 : gen -> gen
val ell_get_c6 : gen -> gen
val ell_get_disc : gen -> gen
val ell_get_j : gen -> gen
val ell_get_type : gen -> Signed.long
val ellff_get_field : gen -> gen
val ellff_get_a4a6 : gen -> gen
val ellqp_get_zero : gen -> gen
val ellqp_get_prec : gen -> Signed.long
val ellqp_get_p : gen -> gen
val ellr_get_prec : gen -> Signed.long
val ellr_get_sign : gen -> Signed.long
val ellnf_get_nf : gen -> gen
val ellnf_get_bnf : gen -> gen
val checkell_i : gen -> int
val ell_is_inf : gen -> int
val ellinf : unit -> gen
val modpr_get_pr : gen -> gen
val modpr_get_p : gen -> gen
val modpr_get_t : gen -> gen
val pr_get_p : gen -> gen
val pr_get_gen : gen -> gen
val pr_get_e : gen -> Signed.long
val pr_get_f : gen -> Signed.long
val pr_get_tau : gen -> gen
val pr_is_inert : gen -> int
val pr_norm : gen -> gen
val upr_norm : gen -> pari_ulong
val nf_get_varn : gen -> Signed.long
val nf_get_pol : gen -> gen
val nf_get_degree : gen -> Signed.long
val nf_get_r1 : gen -> Signed.long
val nf_get_r2 : gen -> Signed.long
val nf_get_disc : gen -> gen
val nf_get_index : gen -> gen
val nf_get_m : gen -> gen
val nf_get_g : gen -> gen
val nf_get_roundg : gen -> gen
val nf_get_tr : gen -> gen
val nf_get_diff : gen -> gen
val nf_get_ramified_primes : gen -> gen
val nf_get_roots : gen -> gen
val nf_get_zk : gen -> gen
val nf_get_zkprimpart : gen -> gen
val nf_get_zkden : gen -> gen
val nf_get_invzk : gen -> gen

val nf_get_sign :
  gen -> Signed.long Ctypes_static.ptr -> Signed.long Ctypes_static.ptr -> unit

val cyc_get_expo : gen -> gen
val abgrp_get_no : gen -> gen
val abgrp_get_cyc : gen -> gen
val abgrp_get_gen : gen -> gen
val bnf_get_nf : gen -> gen
val bnf_get_clgp : gen -> gen
val bnf_get_no : gen -> gen
val bnf_get_cyc : gen -> gen
val bnf_get_gen : gen -> gen
val bnf_get_reg : gen -> gen
val bnf_get_logfu : gen -> gen
val bnf_get_sunits : gen -> gen
val bnf_get_tuu : gen -> gen
val bnf_get_tun : gen -> Signed.long
val bnf_get_fu_nocheck : gen -> gen
val nfv_to_scalar_or_alg : gen -> gen -> gen
val bnf_get_fu : gen -> gen
val bnr_get_bnf : gen -> gen
val bnr_get_bid : gen -> gen
val bnr_get_mod : gen -> gen
val bnr_get_nf : gen -> gen
val bnr_get_clgp : gen -> gen
val bnr_get_no : gen -> gen
val bnr_get_cyc : gen -> gen
val bnr_get_gen_nocheck : gen -> gen
val bnr_get_gen : gen -> gen
val locs_get_cyc : gen -> gen
val locs_get_lsprk : gen -> gen
val locs_get_lgenfil : gen -> gen
val locs_get_mod : gen -> gen
val locs_get_famod : gen -> gen
val locs_get_m_infty : gen -> gen
val gchar_get_basis : gen -> gen
val gchar_get_bnf : gen -> gen
val gchar_get_nf : gen -> gen
val gchar_get_zm : gen -> gen
val gchar_get_mod : gen -> gen
val gchar_get_modp : gen -> gen
val gchar_get_s : gen -> gen
val gchar_get_dldata : gen -> gen
val gchar_get_sfu : gen -> gen
val gchar_get_cyc : gen -> gen
val gchar_get_hnf : gen -> gen
val gchar_get_u : gen -> gen
val gchar_get_ui : gen -> gen
val gchar_get_m0 : gen -> gen
val gchar_get_u0 : gen -> gen
val gchar_get_r1 : gen -> Signed.long
val gchar_get_r2 : gen -> Signed.long
val gchar_get_loccyc : gen -> gen
val gchar_get_nc : gen -> Signed.long
val gchar_get_ns : gen -> Signed.long
val gchar_get_nm : gen -> Signed.long
val gchar_get_evalprec : gen -> Signed.long
val gchar_get_prec : gen -> Signed.long
val gchar_get_nfprec : gen -> Signed.long
val gchar_set_evalprec : gen -> Signed.long -> unit
val gchar_set_prec : gen -> Signed.long -> unit
val gchar_copy_precs : gen -> gen -> unit
val gchar_set_nfprec : gen -> Signed.long -> unit
val gchar_get_ntors : gen -> Signed.long
val gchar_get_nfree : gen -> Signed.long
val gchar_get_nalg : gen -> Signed.long
val gchar_set_basis : gen -> gen -> unit
val gchar_set_nf : gen -> gen -> unit
val gchar_set_ntors : gen -> Signed.long -> unit
val gchar_set_nfree : gen -> Signed.long -> unit
val gchar_set_nalg : gen -> Signed.long -> unit
val gchar_set_cyc : gen -> gen -> unit
val gchar_set_huui : gen -> gen -> gen -> gen -> unit
val gchar_set_m0 : gen -> gen -> unit
val gchar_set_u0 : gen -> gen -> unit
val bid_get_mod : gen -> gen
val bid_get_ideal : gen -> gen
val bid_get_arch : gen -> gen
val bid_get_grp : gen -> gen
val bid_get_fact : gen -> gen
val bid_get_fact2 : gen -> gen
val bid_get_sprk : gen -> gen
val bid_get_sarch : gen -> gen
val bid_get_archp : gen -> gen
val bid_get_u : gen -> gen
val bid_get_no : gen -> gen
val bid_get_cyc : gen -> gen
val bid_get_gen_nocheck : gen -> gen
val bid_get_gen : gen -> gen
val znstar_get_n : gen -> gen
val znstar_get_fan : gen -> gen
val znstar_get_no : gen -> gen
val znstar_get_cyc : gen -> gen
val znstar_get_gen : gen -> gen
val znstar_get_conreycyc : gen -> gen
val znstar_get_conreygen : gen -> gen
val znstar_get_ui : gen -> gen
val znstar_get_u : gen -> gen
val znstar_get_pe : gen -> gen
val gal_get_pol : gen -> gen
val gal_get_p : gen -> gen
val gal_get_e : gen -> gen
val gal_get_mod : gen -> gen
val gal_get_roots : gen -> gen
val gal_get_invvdm : gen -> gen
val gal_get_den : gen -> gen
val gal_get_group : gen -> gen
val gal_get_gen : gen -> gen
val gal_get_orders : gen -> gen
val rnf_get_degree : gen -> Signed.long
val rnf_get_nfdegree : gen -> Signed.long
val rnf_get_absdegree : gen -> Signed.long
val rnf_get_idealdisc : gen -> gen
val rnf_get_k : gen -> gen
val rnf_get_alpha : gen -> gen
val rnf_get_nf : gen -> gen
val rnf_get_nfzk : gen -> gen
val rnf_get_polabs : gen -> gen
val rnf_get_pol : gen -> gen
val rnf_get_disc : gen -> gen
val rnf_get_index : gen -> gen
val rnf_get_ramified_primes : gen -> gen
val rnf_get_varn : gen -> Signed.long
val rnf_get_nfpol : gen -> gen
val rnf_get_nfvarn : gen -> Signed.long
val rnf_get_zk : gen -> gen
val rnf_get_map : gen -> gen
val rnf_get_invzk : gen -> gen
val idealred : gen -> gen -> gen
val idealchineseinit : gen -> gen -> gen
val closure_arity : gen -> Signed.long
val closure_is_variadic : gen -> Signed.long
val closure_codestr : gen -> string
val closure_get_code : gen -> gen
val closure_get_oper : gen -> gen
val closure_get_data : gen -> gen
val closure_get_dbg : gen -> gen
val closure_get_text : gen -> gen
val closure_get_frame : gen -> gen
val err_get_num : gen -> Signed.long
val err_get_compo : gen -> Signed.long -> gen
val pari_err_bug : string -> unit
val pari_err_constpol : string -> unit
val pari_err_coprime : string -> gen -> gen -> unit
