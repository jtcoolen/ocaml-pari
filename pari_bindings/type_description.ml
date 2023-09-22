open Ctypes

module Types (T : Ctypes.TYPE) = struct
  open T

  type pari_ulong = Unsigned.ULong.t

  let pari_ulong : pari_ulong typ = ulong

  type gen = Signed.Long.t ptr

  let gen : gen typ = ptr long

  type byteptr = Unsigned.uchar ptr

  let byteptr : byteptr typ = ptr uchar

  type pari_sp = pari_ulong

  let pari_sp : pari_sp typ = pari_ulong

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

  let logstyle_none = constant "logstyle_none" int64_t
  let logstyle_plain = constant "logstyle_plain" int64_t
  let logstyle_color = constant "logstyle_color" int64_t
  let logstyle_tex = constant "logstyle_TeX" int64_t

  let pari_logstyles =
    enum "pari_logstyles"
      [
        (LOGSTYLE_NONE, logstyle_none);
        (LOGSTYLE_PLAIN, logstyle_plain);
        (LOGSTYLE_COLOR, logstyle_color);
        (LOGSTYLE_TEX, logstyle_tex);
      ]

  let e_syntax = constant "e_SYNTAX" int64_t
  let e_bug = constant "e_BUG" int64_t
  let e_alarm = constant "e_ALARM" int64_t
  let e_file = constant "e_FILE" int64_t
  let e_misc = constant "e_MISC" int64_t
  let e_flag = constant "e_FLAG" int64_t
  let e_impl = constant "e_IMPL" int64_t
  let e_arch = constant "e_ARCH" int64_t
  let e_package = constant "e_PACKAGE" int64_t
  let e_notfunc = constant "e_NOTFUNC" int64_t
  let e_prec = constant "e_PREC" int64_t
  let e_type = constant "e_TYPE" int64_t
  let e_dim = constant "e_DIM" int64_t
  let e_var = constant "e_VAR" int64_t
  let e_priority = constant "e_PRIORITY" int64_t
  let e_user = constant "e_USER" int64_t
  let e_stack = constant "e_STACK" int64_t
  let e_stackthread = constant "e_STACKTHREAD" int64_t
  let e_overflow = constant "e_OVERFLOW" int64_t
  let e_domain = constant "e_DOMAIN" int64_t
  let e_component = constant "e_COMPONENT" int64_t
  let e_maxprime = constant "e_MAXPRIME" int64_t
  let e_constpol = constant "e_CONSTPOL" int64_t
  let e_irredpol = constant "e_IRREDPOL" int64_t
  let e_coprime = constant "e_COPRIME" int64_t
  let e_prime = constant "e_PRIME" int64_t
  let e_modulus = constant "e_MODULUS" int64_t
  let e_roots0 = constant "e_ROOTS0" int64_t
  let e_op = constant "e_OP" int64_t
  let e_type2 = constant "e_TYPE2" int64_t
  let e_inv = constant "e_INV" int64_t
  let e_mem = constant "e_MEM" int64_t
  let e_sqrtn = constant "e_SQRTN" int64_t
  let e_filedesc = constant "e_FILEDESC" int64_t
  let e_none = constant "e_NONE" int64_t

  let err_list =
    enum "err_list"
      [
        (E_SYNTAX, e_syntax);
        (E_BUG, e_bug);
        (E_ALARM, e_alarm);
        (E_FILE, e_file);
        (E_MISC, e_misc);
        (E_FLAG, e_flag);
        (E_IMPL, e_impl);
        (E_ARCH, e_arch);
        (E_PACKAGE, e_package);
        (E_NOTFUNC, e_notfunc);
        (E_PREC, e_prec);
        (E_TYPE, e_type);
        (E_DIM, e_dim);
        (E_VAR, e_var);
        (E_PRIORITY, e_priority);
        (E_USER, e_user);
        (E_STACK, e_stack);
        (E_STACKTHREAD, e_stackthread);
        (E_OVERFLOW, e_overflow);
        (E_DOMAIN, e_domain);
        (E_COMPONENT, e_component);
        (E_MAXPRIME, e_maxprime);
        (E_CONSTPOL, e_constpol);
        (E_IRREDPOL, e_irredpol);
        (E_COPRIME, e_coprime);
        (E_PRIME, e_prime);
        (E_MODULUS, e_modulus);
        (E_ROOTS0, e_roots0);
        (E_OP, e_op);
        (E_TYPE2, e_type2);
        (E_INV, e_inv);
        (E_MEM, e_mem);
        (E_SQRTN, e_sqrtn);
        (E_FILEDESC, e_filedesc);
        (E_NONE, e_none);
      ]

  let pari_timer : pari_timer structure typ =
    typedef (structure "pari_timer") "pari_timer"

  let pari_timer_s = field pari_timer "s" long
  let pari_timer_us = field pari_timer "us" long
  let () = seal pari_timer

  let pari_str : pari_str structure typ =
    typedef (structure "pari_str") "pari_str"

  let pari_str_string = field pari_str "string" string
  let pari_str_end = field pari_str "end" string
  let pari_str_cur = field pari_str "cur" string
  let pari_str_size = field pari_str "size" int
  let pari_str_use_stack = field pari_str "use_stack" int
  let () = seal pari_str
  let pari_sieve : pari_sieve structure typ = structure "pari_sieve"
  let pari_sieve_start = field pari_sieve "start" pari_ulong
  let pari_sieve_end = field pari_sieve "end" pari_ulong
  let pari_sieve_maxpos = field pari_sieve "maxpos" pari_ulong
  let pari_sieve_c = field pari_sieve "c" pari_ulong
  let pari_sieve_q = field pari_sieve "q" pari_ulong
  let pari_sieve_sieve = field pari_sieve "sieve" (ptr uchar)
  let () = seal pari_sieve

  let forprime_t : forprime_t structure typ =
    typedef (structure "forprime_t") "forprime_t"

  let forprime_t_strategy = field forprime_t "strategy" int
  let forprime_t_bb = field forprime_t "bb" gen
  let forprime_t_c = field forprime_t "c" pari_ulong
  let forprime_t_q = field forprime_t "q" pari_ulong
  let forprime_t_d = field forprime_t "d" byteptr
  let forprime_t_p = field forprime_t "p" pari_ulong
  let forprime_t_b = field forprime_t "b" pari_ulong
  let forprime_t_psieve = field forprime_t "psieve" (ptr pari_sieve)
  let forprime_t_sieve = field forprime_t "sieve" (ptr uchar)
  let forprime_t_isieve = field forprime_t "isieve" (ptr uchar)
  let forprime_t_cache = field forprime_t "cache" (array 9 pari_ulong)
  let forprime_t_chunk = field forprime_t "chunk" pari_ulong
  let forprime_t_a = field forprime_t "a" pari_ulong
  let forprime_t_end = field forprime_t "end" pari_ulong
  let forprime_t_sieveb = field forprime_t "sieveb" pari_ulong
  let forprime_t_pos = field forprime_t "pos" pari_ulong
  let forprime_t_maxpos = field forprime_t "maxpos" pari_ulong
  let forprime_t_pp = field forprime_t "pp" gen
  let () = seal forprime_t

  let forcomposite_t : forcomposite_t structure typ =
    typedef (structure "forcomposite_t") "forcomposite_t"

  let forcomposite_t_first = field forcomposite_t "first" int
  let forcomposite_t_b = field forcomposite_t "b" gen
  let forcomposite_t_n = field forcomposite_t "n" gen
  let forcomposite_t_p = field forcomposite_t "p" gen
  let forcomposite_t_T = field forcomposite_t "T" forprime_t
  let () = seal forcomposite_t

  let forvec_t : forvec_t structure typ =
    typedef (structure "forvec_t") "forvec_t"

  let forvec_t_first = field forvec_t "first" long
  let forvec_t_a = field forvec_t "a" (ptr gen)
  let forvec_t_m = field forvec_t "m" (ptr gen)
  let forvec_t_M = field forvec_t "M" (ptr gen)
  let forvec_t_n = field forvec_t "n" long

  let forvec_t_next =
    field forvec_t "next" (static_funptr T.(ptr forvec_t @-> returning gen))

  let () = seal forvec_t

  let forpart_t : forpart_t structure typ =
    typedef (structure "forpart_t") "forpart_t"

  let forpart_t_k = field forpart_t "k" long
  let forpart_t_amax = field forpart_t "amax" long
  let forpart_t_amin = field forpart_t "amin" long
  let forpart_t_nmin = field forpart_t "nmin" long
  let forpart_t_nmax = field forpart_t "nmax" long
  let forpart_t_strip = field forpart_t "strip" long
  let forpart_t_v = field forpart_t "v" gen
  let () = seal forpart_t

  let forperm_t : forperm_t structure typ =
    typedef (structure "forperm_t") "forperm_t"

  let forperm_t_k = field forperm_t "k" long
  let forperm_t_first = field forperm_t "first" long
  let forperm_t_v = field forperm_t "v" gen
  let () = seal forperm_t

  let forsubset_t : forsubset_t structure typ =
    typedef (structure "forsubset_t") "forsubset_t"

  let forsubset_t_n = field forsubset_t "n" long
  let forsubset_t_k = field forsubset_t "k" long
  let forsubset_t_all = field forsubset_t "all" long
  let forsubset_t_first = field forsubset_t "first" long
  let forsubset_t_v = field forsubset_t "v" gen
  let () = seal forsubset_t

  let pari_plot : pari_plot structure typ =
    typedef (structure "PARI_plot") "PARI_plot"

  let pari_plot_draw =
    field pari_plot "draw"
      (static_funptr
         T.(ptr pari_plot @-> gen @-> gen @-> gen @-> returning void))

  let pari_plot_width = field pari_plot "width" long
  let pari_plot_height = field pari_plot "height" long
  let pari_plot_hunit = field pari_plot "hunit" long
  let pari_plot_vunit = field pari_plot "vunit" long
  let pari_plot_fwidth = field pari_plot "fwidth" long
  let pari_plot_fheight = field pari_plot "fheight" long
  let pari_plot_dwidth = field pari_plot "dwidth" long
  let pari_plot_dheight = field pari_plot "dheight" long
  let () = seal pari_plot
  let genbin : genbin structure typ = typedef (structure "GENbin") "GENbin"
  let genbin_len = field genbin "len" int
  let genbin_x = field genbin "x" gen
  let genbin_base = field genbin "base" gen

  let genbin_rebase =
    field genbin "rebase" (static_funptr T.(gen @-> long @-> returning void))

  let () = seal genbin
  let pari_mainstack : pari_mainstack structure typ = structure "pari_mainstack"
  let pari_mainstack_top = field pari_mainstack "top" pari_sp
  let pari_mainstack_bot = field pari_mainstack "bot" pari_sp
  let pari_mainstack_vbot = field pari_mainstack "vbot" pari_sp
  let pari_mainstack_size = field pari_mainstack "size" int
  let pari_mainstack_rsize = field pari_mainstack "rsize" int
  let pari_mainstack_vsize = field pari_mainstack "vsize" int
  let pari_mainstack_memused = field pari_mainstack "memused" int
  let () = seal pari_mainstack
  let entree : entree structure typ = typedef (structure "entree") "entree"
  let entree_name = field entree "name" string
  let entree_valence = field entree "valence" pari_ulong
  let entree_value = field entree "value" (ptr void)
  let entree_menu = field entree "menu" long
  let entree_code = field entree "code" string
  let entree_help = field entree "help" string
  let entree_pvalue = field entree "pvalue" (ptr void)
  let entree_arity = field entree "arity" long
  let entree_hash = field entree "hash" pari_ulong
  let entree_next = field entree "next" (ptr entree)
  let () = seal entree

  let pari_parsestate : pari_parsestate structure typ =
    structure "pari_parsestate"

  let pari_parsestate_node = field pari_parsestate "node" long
  let pari_parsestate_once = field pari_parsestate "once" int
  let pari_parsestate_discarded = field pari_parsestate "discarded" long
  let pari_parsestate_lex_start = field pari_parsestate "lex_start" string
  let pari_parsestate_lasterror = field pari_parsestate "lasterror" gen
  let () = seal pari_parsestate

  let pari_compilestate : pari_compilestate structure typ =
    structure "pari_compilestate"

  let pari_compilestate_opcode = field pari_compilestate "opcode" long
  let pari_compilestate_operand = field pari_compilestate "operand" long
  let pari_compilestate_accesslex = field pari_compilestate "accesslex" long
  let pari_compilestate_data = field pari_compilestate "data" long
  let pari_compilestate_localvars = field pari_compilestate "localvars" long
  let pari_compilestate_frames = field pari_compilestate "frames" long
  let pari_compilestate_dbginfo = field pari_compilestate "dbginfo" long
  let pari_compilestate_offset = field pari_compilestate "offset" long
  let pari_compilestate_nblex = field pari_compilestate "nblex" long
  let pari_compilestate_dbgstart = field pari_compilestate "dbgstart" string
  let () = seal pari_compilestate
  let pari_mtstate : pari_mtstate structure typ = structure "pari_mtstate"
  let pari_mtstate_pending_threads = field pari_mtstate "pending_threads" long
  let pari_mtstate_is_thread = field pari_mtstate "is_thread" long
  let pari_mtstate_trace_level = field pari_mtstate "trace_level" long
  let () = seal pari_mtstate
  let pari_evalstate : pari_evalstate structure typ = structure "pari_evalstate"
  let pari_evalstate_avma = field pari_evalstate "avma" pari_sp
  let pari_evalstate_sp = field pari_evalstate "sp" long
  let pari_evalstate_rp = field pari_evalstate "rp" long
  let pari_evalstate_var = field pari_evalstate "var" long
  let pari_evalstate_lvars = field pari_evalstate "lvars" long
  let pari_evalstate_locks = field pari_evalstate "locks" long
  let pari_evalstate_prec = field pari_evalstate "prec" long
  let pari_evalstate_trace = field pari_evalstate "trace" long
  let pari_evalstate_mt = field pari_evalstate "mt" pari_mtstate
  let pari_evalstate_comp = field pari_evalstate "comp" pari_compilestate
  let () = seal pari_evalstate
  let pari_varstate : pari_varstate structure typ = structure "pari_varstate"
  let pari_varstate_nvar = field pari_varstate "nvar" long
  let pari_varstate_max_avail = field pari_varstate "max_avail" long
  let pari_varstate_min_priority = field pari_varstate "min_priority" long
  let pari_varstate_max_priority = field pari_varstate "max_priority" long
  let () = seal pari_varstate

  let pari_global_state : pari_global_state structure typ =
    structure "pari_global_state"

  let pari_global_state_debugvar = field pari_global_state "debugvar" long
  let pari_global_state_bitprec = field pari_global_state "bitprec" long
  let pari_global_state_primetab = field pari_global_state "primetab" gen
  let pari_global_state_seadata = field pari_global_state "seadata" gen

  let pari_global_state_varpriority =
    field pari_global_state "varpriority" (ptr long)

  let pari_global_state_varstate =
    field pari_global_state "varstate" pari_varstate

  let () = seal pari_global_state
  let pari_thread : pari_thread structure typ = structure "pari_thread"
  let pari_thread_st = field pari_thread "st" pari_mainstack
  let pari_thread_gs = field pari_thread "gs" pari_global_state
  let pari_thread_data = field pari_thread "data" gen
  let () = seal pari_thread
  let mt_state : mt_state structure typ = structure "mt_state"
  let mt_state_worker = field mt_state "worker" gen
  let mt_state_pending = field mt_state "pending" gen
  let mt_state_workid = field mt_state "workid" long
  let () = seal mt_state
  let pari_mt : pari_mt structure typ = structure "pari_mt"
  let pari_mt_mt = field pari_mt "mt" mt_state

  let pari_mt_get =
    field pari_mt "get"
      (static_funptr
         T.(ptr mt_state @-> ptr long @-> ptr long @-> returning gen))

  let pari_mt_submit =
    field pari_mt "submit"
      (static_funptr T.(ptr mt_state @-> long @-> gen @-> returning void))

  let pari_mt_end = field pari_mt "end" (static_funptr T.(returning void))
  let () = seal pari_mt
  let parfor_iter : parfor_iter structure typ = structure "parfor_iter"
  let parfor_iter_pending = field parfor_iter "pending" long
  let parfor_iter_worker = field parfor_iter "worker" gen
  let parfor_iter_pt = field parfor_iter "pt" pari_mt
  let () = seal parfor_iter

  let parfor_t : parfor_t structure typ =
    typedef (structure "parfor_t") "parfor_t"

  let parfor_t_a = field parfor_t "a" gen
  let parfor_t_b = field parfor_t "b" gen
  let parfor_t_iter = field parfor_t "iter" parfor_iter
  let () = seal parfor_t

  let parforeach_t : parforeach_t structure typ =
    typedef (structure "parforeach_t") "parforeach_t"

  let parforeach_t_x = field parforeach_t "x" gen
  let parforeach_t_W = field parforeach_t "W" gen
  let parforeach_t_i = field parforeach_t "i" long
  let parforeach_t_l = field parforeach_t "l" long
  let parforeach_t_iter = field parforeach_t "iter" parfor_iter
  let () = seal parforeach_t

  let parforprime_t : parforprime_t structure typ =
    typedef (structure "parforprime_t") "parforprime_t"

  let parforprime_t_v = field parforprime_t "v" gen
  let parforprime_t_forprime = field parforprime_t "forprime" forprime_t
  let parforprime_t_iter = field parforprime_t "iter" parfor_iter
  let () = seal parforprime_t

  let parforvec_t : parforvec_t structure typ =
    typedef (structure "parforvec_t") "parforvec_t"

  let parforvec_t_v = field parforvec_t "v" gen
  let parforvec_t_forvec = field parforvec_t "forvec" forvec_t
  let parforvec_t_iter = field parforvec_t "iter" parfor_iter
  let () = seal parforvec_t

  let hashentry : hashentry structure typ =
    typedef (structure "hashentry") "hashentry"

  let hashentry_key = field hashentry "key" (ptr void)
  let hashentry_val = field hashentry "val" (ptr void)
  let hashentry_hash = field hashentry "hash" pari_ulong
  let hashentry_next = field hashentry "next" (ptr hashentry)
  let () = seal hashentry

  let hashtable : hashtable structure typ =
    typedef (structure "hashtable") "hashtable"

  let hashtable_len = field hashtable "len" pari_ulong
  let hashtable_table = field hashtable "table" (ptr (ptr hashentry))
  let hashtable_nb = field hashtable "nb" pari_ulong
  let hashtable_maxnb = field hashtable "maxnb" pari_ulong
  let hashtable_pindex = field hashtable "pindex" pari_ulong

  let hashtable_hash =
    field hashtable "hash" (static_funptr T.(ptr void @-> returning pari_ulong))

  let hashtable_eq =
    field hashtable "eq"
      (static_funptr T.(ptr void @-> ptr void @-> returning int))

  let hashtable_use_stack = field hashtable "use_stack" int
  let () = seal hashtable
  let gp_path : gp_path structure typ = typedef (structure "gp_path") "gp_path"
  let gp_path_PATH = field gp_path "PATH" string
  let gp_path_dirs = field gp_path "dirs" (ptr string)
  let () = seal gp_path

  let pariout_t : pariout_t structure typ =
    typedef (structure "pariout_t") "pariout_t"

  let pariout_t_format = field pariout_t "format" char
  let pariout_t_sigd = field pariout_t "sigd" long
  let pariout_t_sp = field pariout_t "sp" int
  let pariout_t_prettyp = field pariout_t "prettyp" int
  let pariout_t_TeXstyle = field pariout_t "TeXstyle" int
  let () = seal pariout_t

  let nfmaxord_t : nfmaxord_t structure typ =
    typedef (structure "nfmaxord_t") "nfmaxord_t"

  let nfmaxord_t_T = field nfmaxord_t "T" gen
  let nfmaxord_t_dT = field nfmaxord_t "dT" gen
  let nfmaxord_t_T0 = field nfmaxord_t "T0" gen
  let nfmaxord_t_unscale = field nfmaxord_t "unscale" gen
  let nfmaxord_t_dK = field nfmaxord_t "dK" gen
  let nfmaxord_t_index = field nfmaxord_t "index" gen
  let nfmaxord_t_basis = field nfmaxord_t "basis" gen
  let nfmaxord_t_r1 = field nfmaxord_t "r1" long
  let nfmaxord_t_basden = field nfmaxord_t "basden" gen
  let nfmaxord_t_dTP = field nfmaxord_t "dTP" gen
  let nfmaxord_t_dTE = field nfmaxord_t "dTE" gen
  let nfmaxord_t_dKP = field nfmaxord_t "dKP" gen
  let nfmaxord_t_dKE = field nfmaxord_t "dKE" gen
  let nfmaxord_t_certify = field nfmaxord_t "certify" long
  let () = seal nfmaxord_t
  let qfr_data : qfr_data structure typ = structure "qfr_data"
  let qfr_data_D = field qfr_data "D" gen
  let qfr_data_sqrtD = field qfr_data "sqrtD" gen
  let qfr_data_isqrtD = field qfr_data "isqrtD" gen
  let () = seal qfr_data

  let fp_chk_fun : fp_chk_fun structure typ =
    typedef (structure "FP_chk_fun") "FP_chk_fun"

  let fp_chk_fun_f =
    field fp_chk_fun "f" (static_funptr T.(ptr void @-> gen @-> returning gen))

  let fp_chk_fun_f_init =
    field fp_chk_fun "f_init"
      (static_funptr T.(ptr fp_chk_fun @-> gen @-> gen @-> returning gen))

  let fp_chk_fun_f_post =
    field fp_chk_fun "f_post"
      (static_funptr T.(ptr fp_chk_fun @-> gen @-> gen @-> returning gen))

  let fp_chk_fun_data = field fp_chk_fun "data" (ptr void)
  let fp_chk_fun_skipfirst = field fp_chk_fun "skipfirst" long
  let () = seal fp_chk_fun
  let zlog_s : zlog_s structure typ = typedef (structure "zlog_S") "zlog_S"
  let zlog_s_bid = field zlog_s "bid" gen
  let zlog_s_P = field zlog_s "P" gen
  let zlog_s_k = field zlog_s "k" gen
  let zlog_s_sprk = field zlog_s "sprk" gen
  let zlog_s_archp = field zlog_s "archp" gen
  let zlog_s_mod = field zlog_s "mod" gen
  let zlog_s_U = field zlog_s "U" gen
  let zlog_s_hU = field zlog_s "hU" long
  let zlog_s_no2 = field zlog_s "no2" int
  let () = seal zlog_s
  let bb_group : bb_group structure typ = structure "bb_group"

  let bb_group_mul =
    field bb_group "mul"
      (static_funptr T.(ptr void @-> gen @-> gen @-> returning gen))

  let bb_group_pow =
    field bb_group "pow"
      (static_funptr T.(ptr void @-> gen @-> gen @-> returning gen))

  let bb_group_rand =
    field bb_group "rand" (static_funptr T.(ptr void @-> returning gen))

  let bb_group_hash =
    field bb_group "hash" (static_funptr T.(gen @-> returning pari_ulong))

  let bb_group_equal =
    field bb_group "equal" (static_funptr T.(gen @-> gen @-> returning int))

  let bb_group_equal1 =
    field bb_group "equal1" (static_funptr T.(gen @-> returning int))

  let bb_group_easylog =
    field bb_group "easylog"
      (static_funptr T.(ptr void @-> gen @-> gen @-> gen @-> returning gen))

  let () = seal bb_group
  let bb_field : bb_field structure typ = structure "bb_field"

  let bb_field_red =
    field bb_field "red" (static_funptr T.(ptr void @-> gen @-> returning gen))

  let bb_field_add =
    field bb_field "add"
      (static_funptr T.(ptr void @-> gen @-> gen @-> returning gen))

  let bb_field_mul =
    field bb_field "mul"
      (static_funptr T.(ptr void @-> gen @-> gen @-> returning gen))

  let bb_field_neg =
    field bb_field "neg" (static_funptr T.(ptr void @-> gen @-> returning gen))

  let bb_field_inv =
    field bb_field "inv" (static_funptr T.(ptr void @-> gen @-> returning gen))

  let bb_field_equal0 =
    field bb_field "equal0" (static_funptr T.(gen @-> returning int))

  let bb_field_s =
    field bb_field "s" (static_funptr T.(ptr void @-> long @-> returning gen))

  let () = seal bb_field
  let bb_algebra : bb_algebra structure typ = structure "bb_algebra"

  let bb_algebra_red =
    field bb_algebra "red"
      (static_funptr T.(ptr void @-> gen @-> returning gen))

  let bb_algebra_add =
    field bb_algebra "add"
      (static_funptr T.(ptr void @-> gen @-> gen @-> returning gen))

  let bb_algebra_sub =
    field bb_algebra "sub"
      (static_funptr T.(ptr void @-> gen @-> gen @-> returning gen))

  let bb_algebra_mul =
    field bb_algebra "mul"
      (static_funptr T.(ptr void @-> gen @-> gen @-> returning gen))

  let bb_algebra_sqr =
    field bb_algebra "sqr"
      (static_funptr T.(ptr void @-> gen @-> returning gen))

  let bb_algebra_one =
    field bb_algebra "one" (static_funptr T.(ptr void @-> returning gen))

  let bb_algebra_zero =
    field bb_algebra "zero" (static_funptr T.(ptr void @-> returning gen))

  let () = seal bb_algebra
  let bb_ring : bb_ring structure typ = structure "bb_ring"

  let bb_ring_add =
    field bb_ring "add"
      (static_funptr T.(ptr void @-> gen @-> gen @-> returning gen))

  let bb_ring_mul =
    field bb_ring "mul"
      (static_funptr T.(ptr void @-> gen @-> gen @-> returning gen))

  let bb_ring_sqr =
    field bb_ring "sqr" (static_funptr T.(ptr void @-> gen @-> returning gen))

  let () = seal bb_ring
end
