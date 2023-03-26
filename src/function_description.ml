open Ctypes
open Types_generated

module F0 (F : Ctypes.FOREIGN) = struct
  open F

  let buchimag =
    foreign "buchimag" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let buchreal =
    foreign "buchreal"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zidealstar = foreign "zidealstar" (gen @-> gen @-> returning gen)
  let zidealstarinit = foreign "zidealstarinit" (gen @-> gen @-> returning gen)

  let zidealstarinitgen =
    foreign "zidealstarinitgen" (gen @-> gen @-> returning gen)

  let factmod = foreign "factmod" (gen @-> gen @-> returning gen)
  let mpbern = foreign "mpbern" (long @-> long @-> returning void)
  let simplefactmod = foreign "simplefactmod" (gen @-> gen @-> returning gen)
  let listkill = foreign "listkill" (gen @-> returning void)

  let isprincipalforce =
    foreign "isprincipalforce" (gen @-> gen @-> returning gen)

  let isprincipalgen = foreign "isprincipalgen" (gen @-> gen @-> returning gen)

  let isprincipalgenforce =
    foreign "isprincipalgenforce" (gen @-> gen @-> returning gen)

  let f2ms_ker = foreign "F2Ms_ker" (gen @-> long @-> returning gen)
  let f2ms_to_f2m = foreign "F2Ms_to_F2m" (gen @-> long @-> returning gen)
  let f2c_to_zc = foreign "F2c_to_ZC" (gen @-> returning gen)
  let f2c_to_mod = foreign "F2c_to_mod" (gen @-> returning gen)
  let f2m_f2c_gauss = foreign "F2m_F2c_gauss" (gen @-> gen @-> returning gen)

  let f2m_f2c_invimage =
    foreign "F2m_F2c_invimage" (gen @-> gen @-> returning gen)

  let f2m_f2c_mul = foreign "F2m_F2c_mul" (gen @-> gen @-> returning gen)
  let f2m_deplin = foreign "F2m_deplin" (gen @-> returning gen)
  let f2m_det = foreign "F2m_det" (gen @-> returning pari_ulong)
  let f2m_det_sp = foreign "F2m_det_sp" (gen @-> returning pari_ulong)
  let f2m_gauss = foreign "F2m_gauss" (gen @-> gen @-> returning gen)
  let f2m_inv = foreign "F2m_inv" (gen @-> returning gen)
  let f2m_invimage = foreign "F2m_invimage" (gen @-> gen @-> returning gen)
  let f2m_ker = foreign "F2m_ker" (gen @-> returning gen)
  let f2m_ker_sp = foreign "F2m_ker_sp" (gen @-> long @-> returning gen)
  let f2m_mul = foreign "F2m_mul" (gen @-> gen @-> returning gen)
  let f2m_powu = foreign "F2m_powu" (gen @-> pari_ulong @-> returning gen)
  let f2m_rank = foreign "F2m_rank" (gen @-> returning long)
  let f2m_row = foreign "F2m_row" (gen @-> long @-> returning gen)

  let f2m_rowslice =
    foreign "F2m_rowslice" (gen @-> long @-> long @-> returning gen)

  let f2m_to_f2ms = foreign "F2m_to_F2Ms" (gen @-> returning gen)
  let f2m_to_flm = foreign "F2m_to_Flm" (gen @-> returning gen)
  let f2m_to_zm = foreign "F2m_to_ZM" (gen @-> returning gen)
  let f2m_to_mod = foreign "F2m_to_mod" (gen @-> returning gen)
  let f2m_transpose = foreign "F2m_transpose" (gen @-> returning gen)

  let f2v_add_inplace =
    foreign "F2v_add_inplace" (gen @-> gen @-> returning void)

  let f2v_and_inplace =
    foreign "F2v_and_inplace" (gen @-> gen @-> returning void)

  let f2v_dotproduct =
    foreign "F2v_dotproduct" (gen @-> gen @-> returning pari_ulong)

  let f2v_equal0 = foreign "F2v_equal0" (gen @-> returning int)
  let f2v_hamming = foreign "F2v_hamming" (gen @-> returning pari_ulong)

  let f2v_negimply_inplace =
    foreign "F2v_negimply_inplace" (gen @-> gen @-> returning void)

  let f2v_or_inplace = foreign "F2v_or_inplace" (gen @-> gen @-> returning void)
  let f2v_slice = foreign "F2v_slice" (gen @-> long @-> long @-> returning gen)
  let f2v_subset = foreign "F2v_subset" (gen @-> gen @-> returning int)
  let f2v_to_flv = foreign "F2v_to_Flv" (gen @-> returning gen)
  let matid_f2m = foreign "matid_F2m" (long @-> returning gen)

  let f2x_f2xq_eval =
    foreign "F2x_F2xq_eval" (gen @-> gen @-> gen @-> returning gen)

  let f2x_f2xqv_eval =
    foreign "F2x_F2xqV_eval" (gen @-> gen @-> gen @-> returning gen)

  let f2x_frobenius = foreign "F2x_Frobenius" (gen @-> returning gen)
  let f2x_1_add = foreign "F2x_1_add" (gen @-> returning gen)
  let f2x_add = foreign "F2x_add" (gen @-> gen @-> returning gen)
  let f2x_deflate = foreign "F2x_deflate" (gen @-> long @-> returning gen)
  let f2x_degfact = foreign "F2x_degfact" (gen @-> returning gen)
  let f2x_degree = foreign "F2x_degree" (gen @-> returning long)
  let f2x_deriv = foreign "F2x_deriv" (gen @-> returning gen)

  let f2x_divrem =
    foreign "F2x_divrem" (gen @-> gen @-> ptr gen @-> returning gen)

  let f2x_eval = foreign "F2x_eval" (gen @-> pari_ulong @-> returning pari_ulong)

  let f2x_even_odd =
    foreign "F2x_even_odd" (gen @-> ptr gen @-> ptr gen @-> returning void)

  let f2x_extgcd =
    foreign "F2x_extgcd" (gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let f2x_gcd = foreign "F2x_gcd" (gen @-> gen @-> returning gen)
  let f2x_get_red = foreign "F2x_get_red" (gen @-> returning gen)
  let f2x_halfgcd = foreign "F2x_halfgcd" (gen @-> gen @-> returning gen)
  let f2x_issquare = foreign "F2x_issquare" (gen @-> returning int)
  let f2x_matfrobenius = foreign "F2x_matFrobenius" (gen @-> returning gen)
  let f2x_mul = foreign "F2x_mul" (gen @-> gen @-> returning gen)
  let f2x_recip = foreign "F2x_recip" (gen @-> returning gen)
  let f2x_rem = foreign "F2x_rem" (gen @-> gen @-> returning gen)
  let f2x_shift = foreign "F2x_shift" (gen @-> long @-> returning gen)
  let f2x_sqr = foreign "F2x_sqr" (gen @-> returning gen)
  let f2x_sqrt = foreign "F2x_sqrt" (gen @-> returning gen)
  let f2x_to_f2v = foreign "F2x_to_F2v" (gen @-> long @-> returning gen)
  let f2x_to_f2xx = foreign "F2x_to_F2xX" (gen @-> long @-> returning gen)
  let f2x_to_flx = foreign "F2x_to_Flx" (gen @-> returning gen)
  let f2x_to_zx = foreign "F2x_to_ZX" (gen @-> returning gen)
  let f2x_valrem = foreign "F2x_valrem" (gen @-> ptr gen @-> returning long)
  let f2xc_to_flxc = foreign "F2xC_to_FlxC" (gen @-> returning gen)
  let f2xc_to_zxc = foreign "F2xC_to_ZXC" (gen @-> returning gen)
  let f2xv_to_f2m = foreign "F2xV_to_F2m" (gen @-> long @-> returning gen)

  let f2xv_to_flxv_inplace =
    foreign "F2xV_to_FlxV_inplace" (gen @-> returning void)

  let f2xv_to_zxv_inplace =
    foreign "F2xV_to_ZXV_inplace" (gen @-> returning void)

  let f2xx_f2x_add = foreign "F2xX_F2x_add" (gen @-> gen @-> returning gen)
  let f2xx_f2x_mul = foreign "F2xX_F2x_mul" (gen @-> gen @-> returning gen)
  let f2xx_add = foreign "F2xX_add" (gen @-> gen @-> returning gen)
  let f2xx_deriv = foreign "F2xX_deriv" (gen @-> returning gen)

  let f2xx_renormalize =
    foreign "F2xX_renormalize" (gen @-> long @-> returning gen)

  let f2xx_to_kronecker =
    foreign "F2xX_to_Kronecker" (gen @-> long @-> returning gen)

  let f2xx_to_flxx = foreign "F2xX_to_FlxX" (gen @-> returning gen)
  let f2xx_to_zxx = foreign "F2xX_to_ZXX" (gen @-> returning gen)

  let f2xx_to_f2xc =
    foreign "F2xX_to_F2xC" (gen @-> long @-> long @-> returning gen)

  let f2xxv_to_f2xm =
    foreign "F2xXV_to_F2xM" (gen @-> long @-> long @-> returning gen)

  let f2xxc_to_zxxc = foreign "F2xXC_to_ZXXC" (gen @-> returning gen)

  let f2xy_f2xq_evalx =
    foreign "F2xY_F2xq_evalx" (gen @-> gen @-> gen @-> returning gen)

  let f2xy_f2xqv_evalx =
    foreign "F2xY_F2xqV_evalx" (gen @-> gen @-> gen @-> returning gen)

  let f2xy_degreex = foreign "F2xY_degreex" (gen @-> returning long)
  let f2xn_div = foreign "F2xn_div" (gen @-> gen @-> long @-> returning gen)
  let f2xn_inv = foreign "F2xn_inv" (gen @-> long @-> returning gen)
  let f2xn_red = foreign "F2xn_red" (gen @-> long @-> returning gen)

  let f2xq_artin_schreier =
    foreign "F2xq_Artin_Schreier" (gen @-> gen @-> returning gen)
end

module F1 (F : Ctypes.FOREIGN) = struct
  open F

  let f2xq_autpow =
    foreign "F2xq_autpow" (gen @-> long @-> gen @-> returning gen)

  let f2xq_conjvec = foreign "F2xq_conjvec" (gen @-> gen @-> returning gen)
  let f2xq_div = foreign "F2xq_div" (gen @-> gen @-> gen @-> returning gen)
  let f2xq_inv = foreign "F2xq_inv" (gen @-> gen @-> returning gen)
  let f2xq_invsafe = foreign "F2xq_invsafe" (gen @-> gen @-> returning gen)

  let f2xq_log =
    foreign "F2xq_log" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xq_matrix_pow =
    foreign "F2xq_matrix_pow" (gen @-> long @-> long @-> gen @-> returning gen)

  let f2xq_mul = foreign "F2xq_mul" (gen @-> gen @-> gen @-> returning gen)
  let f2xq_order = foreign "F2xq_order" (gen @-> gen @-> gen @-> returning gen)
  let f2xq_pow = foreign "F2xq_pow" (gen @-> gen @-> gen @-> returning gen)

  let f2xq_pow_init =
    foreign "F2xq_pow_init" (gen @-> gen @-> long @-> gen @-> returning gen)

  let f2xq_pow_table =
    foreign "F2xq_pow_table" (gen @-> gen @-> gen @-> returning gen)

  let f2xq_powu =
    foreign "F2xq_powu" (gen @-> pari_ulong @-> gen @-> returning gen)

  let f2xq_powers =
    foreign "F2xq_powers" (gen @-> long @-> gen @-> returning gen)

  let f2xq_sqr = foreign "F2xq_sqr" (gen @-> gen @-> returning gen)
  let f2xq_sqrt = foreign "F2xq_sqrt" (gen @-> gen @-> returning gen)

  let f2xq_sqrt_fast =
    foreign "F2xq_sqrt_fast" (gen @-> gen @-> gen @-> returning gen)

  let f2xq_sqrtn =
    foreign "F2xq_sqrtn" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let f2xq_trace = foreign "F2xq_trace" (gen @-> gen @-> returning pari_ulong)

  let f2xqx_f2xq_mul =
    foreign "F2xqX_F2xq_mul" (gen @-> gen @-> gen @-> returning gen)

  let f2xqx_f2xq_mul_to_monic =
    foreign "F2xqX_F2xq_mul_to_monic" (gen @-> gen @-> gen @-> returning gen)

  let f2xqx_f2xqxq_eval =
    foreign "F2xqX_F2xqXQ_eval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqx_f2xqxqv_eval =
    foreign "F2xqX_F2xqXQV_eval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqx_disc = foreign "F2xqX_disc" (gen @-> gen @-> returning gen)

  let f2xqx_divrem =
    foreign "F2xqX_divrem" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let f2xqx_extgcd =
    foreign "F2xqX_extgcd"
      (gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let f2xqx_gcd = foreign "F2xqX_gcd" (gen @-> gen @-> gen @-> returning gen)
  let f2xqx_get_red = foreign "F2xqX_get_red" (gen @-> gen @-> returning gen)

  let f2xqx_halfgcd =
    foreign "F2xqX_halfgcd" (gen @-> gen @-> gen @-> returning gen)

  let f2xqx_invbarrett =
    foreign "F2xqX_invBarrett" (gen @-> gen @-> returning gen)

  let f2xqx_ispower =
    foreign "F2xqX_ispower" (gen @-> long @-> gen @-> ptr gen @-> returning long)

  let f2xqx_mul = foreign "F2xqX_mul" (gen @-> gen @-> gen @-> returning gen)
  let f2xqx_normalize = foreign "F2xqX_normalize" (gen @-> gen @-> returning gen)

  let f2xqx_powu =
    foreign "F2xqX_powu" (gen @-> pari_ulong @-> gen @-> returning gen)

  let f2xqx_red = foreign "F2xqX_red" (gen @-> gen @-> returning gen)
  let f2xqx_rem = foreign "F2xqX_rem" (gen @-> gen @-> gen @-> returning gen)

  let f2xqx_resultant =
    foreign "F2xqX_resultant" (gen @-> gen @-> gen @-> returning gen)

  let f2xqx_sqr = foreign "F2xqX_sqr" (gen @-> gen @-> returning gen)
  let f2xqxq_inv = foreign "F2xqXQ_inv" (gen @-> gen @-> gen @-> returning gen)

  let f2xqxq_invsafe =
    foreign "F2xqXQ_invsafe" (gen @-> gen @-> gen @-> returning gen)

  let f2xqxq_mul =
    foreign "F2xqXQ_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqxq_sqr = foreign "F2xqXQ_sqr" (gen @-> gen @-> gen @-> returning gen)

  let f2xqxq_pow =
    foreign "F2xqXQ_pow" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqxq_powers =
    foreign "F2xqXQ_powers" (gen @-> long @-> gen @-> gen @-> returning gen)

  let f2xqxq_autpow =
    foreign "F2xqXQ_autpow" (gen @-> long @-> gen @-> gen @-> returning gen)

  let f2xqxq_auttrace =
    foreign "F2xqXQ_auttrace" (gen @-> long @-> gen @-> gen @-> returning gen)

  let f2xqxqv_red = foreign "F2xqXQV_red" (gen @-> gen @-> gen @-> returning gen)
  let flm_to_f2m = foreign "Flm_to_F2m" (gen @-> returning gen)
  let flv_to_f2v = foreign "Flv_to_F2v" (gen @-> returning gen)
  let flx_to_f2x = foreign "Flx_to_F2x" (gen @-> returning gen)
  let flxc_to_f2xc = foreign "FlxC_to_F2xC" (gen @-> returning gen)
  let flxx_to_f2xx = foreign "FlxX_to_F2xX" (gen @-> returning gen)
  let flxxc_to_f2xxc = foreign "FlxXC_to_F2xXC" (gen @-> returning gen)

  let kronecker_to_f2xqx =
    foreign "Kronecker_to_F2xqX" (gen @-> gen @-> returning gen)

  let rg_to_f2xq = foreign "Rg_to_F2xq" (gen @-> gen @-> returning gen)
  let rgm_to_f2m = foreign "RgM_to_F2m" (gen @-> returning gen)
  let rgv_to_f2v = foreign "RgV_to_F2v" (gen @-> returning gen)
  let rgx_to_f2x = foreign "RgX_to_F2x" (gen @-> returning gen)
  let z_to_f2x = foreign "Z_to_F2x" (gen @-> long @-> returning gen)
  let zm_to_f2m = foreign "ZM_to_F2m" (gen @-> returning gen)
  let zv_to_f2v = foreign "ZV_to_F2v" (gen @-> returning gen)
  let zx_to_f2x = foreign "ZX_to_F2x" (gen @-> returning gen)
  let zxx_to_f2xx = foreign "ZXX_to_F2xX" (gen @-> long @-> returning gen)
  let const_f2v = foreign "const_F2v" (long @-> returning gen)
  let gener_f2xq = foreign "gener_F2xq" (gen @-> ptr gen @-> returning gen)

  let get_f2xq_field =
    foreign "get_F2xq_field"
      (ptr (ptr void) @-> gen @-> returning (ptr bb_field))

  let monomial_f2x = foreign "monomial_F2x" (long @-> long @-> returning gen)
  let pol1_f2xx = foreign "pol1_F2xX" (long @-> long @-> returning gen)
  let polx_f2xx = foreign "polx_F2xX" (long @-> long @-> returning gen)

  let random_f2xqx =
    foreign "random_F2xqX" (long @-> long @-> gen @-> returning gen)

  let f2x_teichmuller =
    foreign "F2x_Teichmuller" (gen @-> long @-> returning gen)

  let f2xq_ellcard =
    foreign "F2xq_ellcard" (gen @-> gen @-> gen @-> returning gen)

  let f2xq_ellgens =
    foreign "F2xq_ellgens"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xq_ellgroup =
    foreign "F2xq_ellgroup"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let f2xq_elltwist =
    foreign "F2xq_elltwist"
      (gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning void)

  let f2xqe_add =
    foreign "F2xqE_add" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqe_changepoint =
    foreign "F2xqE_changepoint" (gen @-> gen @-> gen @-> returning gen)

  let f2xqe_changepointinv =
    foreign "F2xqE_changepointinv" (gen @-> gen @-> gen @-> returning gen)

  let f2xqe_dbl = foreign "F2xqE_dbl" (gen @-> gen @-> gen @-> returning gen)

  let f2xqe_log =
    foreign "F2xqE_log" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqe_mul =
    foreign "F2xqE_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqe_neg = foreign "F2xqE_neg" (gen @-> gen @-> gen @-> returning gen)

  let f2xqe_order =
    foreign "F2xqE_order" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqe_sub =
    foreign "F2xqE_sub" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqe_tatepairing =
    foreign "F2xqE_tatepairing"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let f2xqe_weilpairing =
    foreign "F2xqE_weilpairing"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let get_f2xqe_group =
    foreign "get_F2xqE_group"
      (ptr (ptr void) @-> gen @-> gen @-> gen @-> returning (ptr bb_group))

  let rge_to_f2xqe = foreign "RgE_to_F2xqE" (gen @-> gen @-> returning gen)

  let random_f2xqe =
    foreign "random_F2xqE" (gen @-> gen @-> gen @-> returning gen)

  let f3c_to_mod = foreign "F3c_to_mod" (gen @-> returning gen)
  let f3c_to_zc = foreign "F3c_to_ZC" (gen @-> returning gen)
  let f3m_ker = foreign "F3m_ker" (gen @-> returning gen)
  let f3m_ker_sp = foreign "F3m_ker_sp" (gen @-> long @-> returning gen)
  let f3m_mul = foreign "F3m_mul" (gen @-> gen @-> returning gen)
  let f3m_row = foreign "F3m_row" (gen @-> long @-> returning gen)
  let f3m_to_flm = foreign "F3m_to_Flm" (gen @-> returning gen)
  let f3m_to_zm = foreign "F3m_to_ZM" (gen @-> returning gen)
  let f3m_to_mod = foreign "F3m_to_mod" (gen @-> returning gen)
  let f3m_transpose = foreign "F3m_transpose" (gen @-> returning gen)
  let f3v_to_flv = foreign "F3v_to_Flv" (gen @-> returning gen)
end

module F2 (F : Ctypes.FOREIGN) = struct
  open F

  let f3v_coeff = foreign "F3v_coeff" (gen @-> long @-> returning pari_ulong)
  let f3v_clear = foreign "F3v_clear" (gen @-> long @-> returning void)

  let f3v_set =
    foreign "F3v_set" (gen @-> long @-> pari_ulong @-> returning void)

  let flm_to_f3m = foreign "Flm_to_F3m" (gen @-> returning gen)
  let flv_to_f3v = foreign "Flv_to_F3v" (gen @-> returning gen)
  let rgm_to_f3m = foreign "RgM_to_F3m" (gen @-> returning gen)
  let rgv_to_f3v = foreign "RgV_to_F3v" (gen @-> returning gen)
  let zm_to_f3m = foreign "ZM_to_F3m" (gen @-> returning gen)
  let zv_to_f3v = foreign "ZV_to_F3v" (gen @-> returning gen)
  let zero_f3m_copy = foreign "zero_F3m_copy" (long @-> long @-> returning gen)
  let zero_f3v = foreign "zero_F3v" (long @-> returning gen)

  let fl_elldisc =
    foreign "Fl_elldisc"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_elldisc_pre =
    foreign "Fl_elldisc_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fl_ellj =
    foreign "Fl_ellj"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_ellj_to_a4a6 =
    foreign "Fl_ellj_to_a4a6"
      (pari_ulong @-> pari_ulong @-> ptr pari_ulong @-> ptr pari_ulong
     @-> returning void)

  let fl_ellptors =
    foreign "Fl_ellptors"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let fl_elltwist =
    foreign "Fl_elltwist"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> ptr pari_ulong
     @-> ptr pari_ulong @-> returning void)

  let fl_elltwist_disc =
    foreign "Fl_elltwist_disc"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> ptr pari_ulong @-> ptr pari_ulong @-> returning void)

  let fle_add =
    foreign "Fle_add"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fle_dbl =
    foreign "Fle_dbl" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fle_changepoint =
    foreign "Fle_changepoint" (gen @-> gen @-> pari_ulong @-> returning gen)

  let fle_changepointinv =
    foreign "Fle_changepointinv" (gen @-> gen @-> pari_ulong @-> returning gen)

  let fle_log =
    foreign "Fle_log"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fle_mul =
    foreign "Fle_mul"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fle_mulu =
    foreign "Fle_mulu"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fle_order =
    foreign "Fle_order"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fle_sub =
    foreign "Fle_sub"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fle_tatepairing =
    foreign "Fle_tatepairing"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fle_to_flj = foreign "Fle_to_Flj" (gen @-> returning gen)

  let fle_weilpairing =
    foreign "Fle_weilpairing"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let flj_add_pre =
    foreign "Flj_add_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flj_changepointinv_pre =
    foreign "Flj_changepointinv_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flj_dbl_pre =
    foreign "Flj_dbl_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flj_mulu_pre =
    foreign "Flj_mulu_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flj_neg = foreign "Flj_neg" (gen @-> pari_ulong @-> returning gen)
  let flj_to_fle = foreign "Flj_to_Fle" (gen @-> pari_ulong @-> returning gen)

  let flj_to_fle_pre =
    foreign "Flj_to_Fle_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fljv_factorback_pre =
    foreign "FljV_factorback_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let random_fle =
    foreign "random_Fle"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let random_fle_pre =
    foreign "random_Fle_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let random_flj_pre =
    foreign "random_Flj_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flc_to_zc = foreign "Flc_to_ZC" (gen @-> returning gen)
  let flc_to_zc_inplace = foreign "Flc_to_ZC_inplace" (gen @-> returning gen)

  let flm_flc_gauss =
    foreign "Flm_Flc_gauss" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flm_flc_invimage =
    foreign "Flm_Flc_invimage" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flm_adjoint = foreign "Flm_adjoint" (gen @-> pari_ulong @-> returning gen)
  let flm_deplin = foreign "Flm_deplin" (gen @-> pari_ulong @-> returning gen)
  let flm_det = foreign "Flm_det" (gen @-> pari_ulong @-> returning pari_ulong)

  let flm_det_sp =
    foreign "Flm_det_sp" (gen @-> pari_ulong @-> returning pari_ulong)

  let flm_gauss =
    foreign "Flm_gauss" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flm_intersect =
    foreign "Flm_intersect" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flm_intersect_i =
    foreign "Flm_intersect_i" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flm_inv = foreign "Flm_inv" (gen @-> pari_ulong @-> returning gen)

  let flm_invimage =
    foreign "Flm_invimage" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flm_ker = foreign "Flm_ker" (gen @-> pari_ulong @-> returning gen)

  let flm_ker_sp =
    foreign "Flm_ker_sp" (gen @-> pari_ulong @-> long @-> returning gen)

  let flm_rank = foreign "Flm_rank" (gen @-> pari_ulong @-> returning long)
  let flm_to_zm = foreign "Flm_to_ZM" (gen @-> returning gen)
  let flm_to_zm_inplace = foreign "Flm_to_ZM_inplace" (gen @-> returning gen)
  let flv_to_zv = foreign "Flv_to_ZV" (gen @-> returning gen)
  let fl_to_flx = foreign "Fl_to_Flx" (pari_ulong @-> long @-> returning gen)
  let fl2_equal1 = foreign "Fl2_equal1" (gen @-> returning int)

  let fl2_inv_pre =
    foreign "Fl2_inv_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fl2_mul_pre =
    foreign "Fl2_mul_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let fl2_norm_pre =
    foreign "Fl2_norm_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl2_pow_pre =
    foreign "Fl2_pow_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let fl2_sqr_pre =
    foreign "Fl2_sqr_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fl2_sqrtn_pre =
    foreign "Fl2_sqrtn_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> ptr gen
     @-> returning gen)

  let flm_to_flxv = foreign "Flm_to_FlxV" (gen @-> long @-> returning gen)

  let flm_to_flxx =
    foreign "Flm_to_FlxX" (gen @-> long @-> long @-> returning gen)

  let flv_flm_polint =
    foreign "Flv_Flm_polint"
      (gen @-> gen @-> pari_ulong @-> long @-> returning gen)

  let flv_inv = foreign "Flv_inv" (gen @-> pari_ulong @-> returning gen)

  let flv_inv_inplace =
    foreign "Flv_inv_inplace" (gen @-> pari_ulong @-> returning void)

  let flv_inv_pre_inplace =
    foreign "Flv_inv_pre_inplace"
      (gen @-> pari_ulong @-> pari_ulong @-> returning void)

  let flv_inv_pre =
    foreign "Flv_inv_pre" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flv_invvandermonde =
    foreign "Flv_invVandermonde"
      (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flv_polint =
    foreign "Flv_polint" (gen @-> gen @-> pari_ulong @-> long @-> returning gen)

  let flv_prod = foreign "Flv_prod" (gen @-> pari_ulong @-> returning pari_ulong)

  let flv_prod_pre =
    foreign "Flv_prod_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let flv_roots_to_pol =
    foreign "Flv_roots_to_pol" (gen @-> pari_ulong @-> long @-> returning gen)

  let flv_to_flx = foreign "Flv_to_Flx" (gen @-> long @-> returning gen)

  let flx_fl_add =
    foreign "Flx_Fl_add" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_fl_mul =
    foreign "Flx_Fl_mul" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_fl_mul_to_monic =
    foreign "Flx_Fl_mul_to_monic"
      (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_fl_sub =
    foreign "Flx_Fl_sub" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_fl2_eval_pre =
    foreign "Flx_Fl2_eval_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flx_flv_multieval =
    foreign "Flx_Flv_multieval" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_flxq_eval =
    foreign "Flx_Flxq_eval"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_flxq_eval_pre =
    foreign "Flx_Flxq_eval_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_flxqv_eval =
    foreign "Flx_FlxqV_eval"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_flxqv_eval_pre =
    foreign "Flx_FlxqV_eval_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_frobenius =
    foreign "Flx_Frobenius" (gen @-> pari_ulong @-> returning gen)

  let flx_frobenius_pre =
    foreign "Flx_Frobenius_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_laplace = foreign "Flx_Laplace" (gen @-> pari_ulong @-> returning gen)

  let flx_newton =
    foreign "Flx_Newton" (gen @-> long @-> pari_ulong @-> returning gen)

  let flx_add = foreign "Flx_add" (gen @-> gen @-> pari_ulong @-> returning gen)
  let flx_blocks = foreign "Flx_blocks" (gen @-> long @-> long @-> returning gen)
  let flx_deflate = foreign "Flx_deflate" (gen @-> long @-> returning gen)
  let flx_deriv = foreign "Flx_deriv" (gen @-> pari_ulong @-> returning gen)
  let flx_diff1 = foreign "Flx_diff1" (gen @-> pari_ulong @-> returning gen)
end

module F3 (F : Ctypes.FOREIGN) = struct
  open F

  let flx_digits =
    foreign "Flx_digits" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_div_by_x_x =
    foreign "Flx_div_by_X_x"
      (gen @-> pari_ulong @-> pari_ulong @-> ptr pari_ulong @-> returning gen)

  let flx_divrem =
    foreign "Flx_divrem"
      (gen @-> gen @-> pari_ulong @-> ptr gen @-> returning gen)

  let flx_divrem_pre =
    foreign "Flx_divrem_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> ptr gen @-> returning gen)

  let flx_double = foreign "Flx_double" (gen @-> pari_ulong @-> returning gen)
  let flx_equal = foreign "Flx_equal" (gen @-> gen @-> returning int)

  let flx_eval =
    foreign "Flx_eval"
      (gen @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let flx_eval_powers_pre =
    foreign "Flx_eval_powers_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let flx_eval_pre =
    foreign "Flx_eval_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let flx_extgcd =
    foreign "Flx_extgcd"
      (gen @-> gen @-> pari_ulong @-> ptr gen @-> ptr gen @-> returning gen)

  let flx_extgcd_pre =
    foreign "Flx_extgcd_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> ptr gen @-> ptr gen
     @-> returning gen)

  let flx_extresultant =
    foreign "Flx_extresultant"
      (gen @-> gen @-> pari_ulong @-> ptr gen @-> ptr gen
     @-> returning pari_ulong)

  let flx_fromnewton =
    foreign "Flx_fromNewton" (gen @-> pari_ulong @-> returning gen)

  let flx_gcd = foreign "Flx_gcd" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_gcd_pre =
    foreign "Flx_gcd_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_get_red = foreign "Flx_get_red" (gen @-> pari_ulong @-> returning gen)

  let flx_get_red_pre =
    foreign "Flx_get_red_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_halfgcd =
    foreign "Flx_halfgcd" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_halfgcd_pre =
    foreign "Flx_halfgcd_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_halve = foreign "Flx_halve" (gen @-> pari_ulong @-> returning gen)
  let flx_inflate = foreign "Flx_inflate" (gen @-> long @-> returning gen)
  let flx_integ = foreign "Flx_integ" (gen @-> pari_ulong @-> returning gen)

  let flx_invbarrett =
    foreign "Flx_invBarrett" (gen @-> pari_ulong @-> returning gen)

  let flx_invlaplace =
    foreign "Flx_invLaplace" (gen @-> pari_ulong @-> returning gen)

  let flx_is_squarefree =
    foreign "Flx_is_squarefree" (gen @-> pari_ulong @-> returning int)

  let flx_is_smooth =
    foreign "Flx_is_smooth" (gen @-> long @-> pari_ulong @-> returning int)

  let flx_is_smooth_pre =
    foreign "Flx_is_smooth_pre"
      (gen @-> long @-> pari_ulong @-> pari_ulong @-> returning int)

  let flx_matfrobenius =
    foreign "Flx_matFrobenius" (gen @-> pari_ulong @-> returning gen)

  let flx_matfrobenius_pre =
    foreign "Flx_matFrobenius_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_mod_xn1 =
    foreign "Flx_mod_Xn1" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_mod_xnm1 =
    foreign "Flx_mod_Xnm1" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_mul = foreign "Flx_mul" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_mul_pre =
    foreign "Flx_mul_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_neg = foreign "Flx_neg" (gen @-> pari_ulong @-> returning gen)

  let flx_neg_inplace =
    foreign "Flx_neg_inplace" (gen @-> pari_ulong @-> returning gen)

  let flx_normalize =
    foreign "Flx_normalize" (gen @-> pari_ulong @-> returning gen)

  let flx_powu =
    foreign "Flx_powu" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_powu_pre =
    foreign "Flx_powu_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_recip = foreign "Flx_recip" (gen @-> returning gen)
  let flx_red = foreign "Flx_red" (gen @-> pari_ulong @-> returning gen)
  let flx_rem = foreign "Flx_rem" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_rem_pre =
    foreign "Flx_rem_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_renormalize =
    foreign "Flx_renormalize" (gen @-> long @-> returning gen)

  let flx_rescale =
    foreign "Flx_rescale" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_resultant =
    foreign "Flx_resultant" (gen @-> gen @-> pari_ulong @-> returning pari_ulong)

  let flx_resultant_pre =
    foreign "Flx_resultant_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let flx_shift = foreign "Flx_shift" (gen @-> long @-> returning gen)
  let flx_splitting = foreign "Flx_splitting" (gen @-> long @-> returning gen)
  let flx_sqr = foreign "Flx_sqr" (gen @-> pari_ulong @-> returning gen)

  let flx_sqr_pre =
    foreign "Flx_sqr_pre" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_sub = foreign "Flx_sub" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_translate1 =
    foreign "Flx_translate1" (gen @-> pari_ulong @-> returning gen)

  let flx_translate1_basecase =
    foreign "Flx_translate1_basecase" (gen @-> pari_ulong @-> returning gen)

  let flx_to_flv = foreign "Flx_to_Flv" (gen @-> long @-> returning gen)
  let flx_to_flxx = foreign "Flx_to_FlxX" (gen @-> long @-> returning gen)
  let flx_to_zx = foreign "Flx_to_ZX" (gen @-> returning gen)
  let flx_to_zx_inplace = foreign "Flx_to_ZX_inplace" (gen @-> returning gen)
  let flx_triple = foreign "Flx_triple" (gen @-> pari_ulong @-> returning gen)
  let flx_val = foreign "Flx_val" (gen @-> returning long)
  let flx_valrem = foreign "Flx_valrem" (gen @-> ptr gen @-> returning long)

  let flxc_flxqv_eval =
    foreign "FlxC_FlxqV_eval"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxc_flxqv_eval_pre =
    foreign "FlxC_FlxqV_eval_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxc_flxq_eval =
    foreign "FlxC_Flxq_eval"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxc_flxq_eval_pre =
    foreign "FlxC_Flxq_eval_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxc_eval_powers_pre =
    foreign "FlxC_eval_powers_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxc_neg = foreign "FlxC_neg" (gen @-> pari_ulong @-> returning gen)

  let flxc_sub =
    foreign "FlxC_sub" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxc_to_zxc = foreign "FlxC_to_ZXC" (gen @-> returning gen)

  let flxm_flx_add_shallow =
    foreign "FlxM_Flx_add_shallow" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxm_eval_powers_pre =
    foreign "FlxM_eval_powers_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxm_neg = foreign "FlxM_neg" (gen @-> pari_ulong @-> returning gen)

  let flxm_sub =
    foreign "FlxM_sub" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxm_to_flxxv = foreign "FlxM_to_FlxXV" (gen @-> long @-> returning gen)
  let flxm_to_zxm = foreign "FlxM_to_ZXM" (gen @-> returning gen)
  let flxt_red = foreign "FlxT_red" (gen @-> pari_ulong @-> returning gen)

  let flxv_flc_mul =
    foreign "FlxV_Flc_mul" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxv_flv_multieval =
    foreign "FlxV_Flv_multieval" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxv_flx_fromdigits =
    foreign "FlxV_Flx_fromdigits" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxv_prod = foreign "FlxV_prod" (gen @-> pari_ulong @-> returning gen)
  let flxv_red = foreign "FlxV_red" (gen @-> pari_ulong @-> returning gen)
  let flxv_to_flm = foreign "FlxV_to_Flm" (gen @-> long @-> returning gen)
  let flxv_to_flxx = foreign "FlxV_to_FlxX" (gen @-> long @-> returning gen)
  let flxv_to_zxv = foreign "FlxV_to_ZXV" (gen @-> returning gen)

  let flxv_to_zxv_inplace =
    foreign "FlxV_to_ZXV_inplace" (gen @-> returning void)

  let flxn_div =
    foreign "Flxn_div" (gen @-> gen @-> long @-> pari_ulong @-> returning gen)

  let flxn_div_pre =
    foreign "Flxn_div_pre"
      (gen @-> gen @-> long @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxn_exp =
    foreign "Flxn_exp" (gen @-> long @-> pari_ulong @-> returning gen)

  let flxn_expint =
    foreign "Flxn_expint" (gen @-> long @-> pari_ulong @-> returning gen)

  let flxn_inv =
    foreign "Flxn_inv" (gen @-> long @-> pari_ulong @-> returning gen)

  let flxn_mul =
    foreign "Flxn_mul" (gen @-> gen @-> long @-> pari_ulong @-> returning gen)

  let flxn_mul_pre =
    foreign "Flxn_mul_pre"
      (gen @-> gen @-> long @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxn_sqr =
    foreign "Flxn_sqr" (gen @-> long @-> pari_ulong @-> returning gen)

  let flxn_sqr_pre =
    foreign "Flxn_sqr_pre"
      (gen @-> long @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxn_red = foreign "Flxn_red" (gen @-> long @-> returning gen)

  let flxq_autpow =
    foreign "Flxq_autpow"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> returning gen)

  let flxq_autpow_pre =
    foreign "Flxq_autpow_pre"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxq_autpowers =
    foreign "Flxq_autpowers"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> returning gen)

  let flxq_autsum =
    foreign "Flxq_autsum"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> returning gen)

  let flxq_auttrace =
    foreign "Flxq_auttrace"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> returning gen)

  let flxq_auttrace_pre =
    foreign "Flxq_auttrace_pre"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)
end

module F4 (F : Ctypes.FOREIGN) = struct
  open F

  let flxq_charpoly =
    foreign "Flxq_charpoly" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_conjvec =
    foreign "Flxq_conjvec" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_div =
    foreign "Flxq_div" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_div_pre =
    foreign "Flxq_div_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxq_inv =
    foreign "Flxq_inv" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_inv_pre =
    foreign "Flxq_inv_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxq_invsafe =
    foreign "Flxq_invsafe" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_invsafe_pre =
    foreign "Flxq_invsafe_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxq_issquare =
    foreign "Flxq_issquare" (gen @-> gen @-> pari_ulong @-> returning int)

  let flxq_is2npower =
    foreign "Flxq_is2npower"
      (gen @-> long @-> gen @-> pari_ulong @-> returning int)

  let flxq_log =
    foreign "Flxq_log"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_lroot = foreign "Flxq_lroot" (gen @-> gen @-> long @-> returning gen)

  let flxq_lroot_pre =
    foreign "Flxq_lroot_pre"
      (gen @-> gen @-> long @-> pari_ulong @-> returning gen)

  let flxq_lroot_fast =
    foreign "Flxq_lroot_fast" (gen @-> gen @-> gen @-> long @-> returning gen)

  let flxq_lroot_fast_pre =
    foreign "Flxq_lroot_fast_pre"
      (gen @-> gen @-> gen @-> long @-> pari_ulong @-> returning gen)

  let flxq_matrix_pow =
    foreign "Flxq_matrix_pow"
      (gen @-> long @-> long @-> gen @-> pari_ulong @-> returning gen)

  let flxq_matrix_pow_pre =
    foreign "Flxq_matrix_pow_pre"
      (gen @-> long @-> long @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxq_minpoly =
    foreign "Flxq_minpoly" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_minpoly_pre =
    foreign "Flxq_minpoly_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxq_mul =
    foreign "Flxq_mul" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_mul_pre =
    foreign "Flxq_mul_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxq_norm =
    foreign "Flxq_norm" (gen @-> gen @-> pari_ulong @-> returning pari_ulong)

  let flxq_order =
    foreign "Flxq_order" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_pow =
    foreign "Flxq_pow" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_pow_pre =
    foreign "Flxq_pow_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxq_pow_init =
    foreign "Flxq_pow_init"
      (gen @-> gen @-> long @-> gen @-> pari_ulong @-> returning gen)

  let flxq_pow_init_pre =
    foreign "Flxq_pow_init_pre"
      (gen @-> gen @-> long @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxq_pow_table_pre =
    foreign "Flxq_pow_table_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxq_pow_table =
    foreign "Flxq_pow_table"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_powu =
    foreign "Flxq_powu"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> returning gen)

  let flxq_powu_pre =
    foreign "Flxq_powu_pre"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxq_powers =
    foreign "Flxq_powers" (gen @-> long @-> gen @-> pari_ulong @-> returning gen)

  let flxq_powers_pre =
    foreign "Flxq_powers_pre"
      (gen @-> long @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxq_sqr =
    foreign "Flxq_sqr" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_sqr_pre =
    foreign "Flxq_sqr_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxq_sqrt =
    foreign "Flxq_sqrt" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_sqrtn =
    foreign "Flxq_sqrtn"
      (gen @-> gen @-> gen @-> pari_ulong @-> ptr gen @-> returning gen)

  let flxq_trace =
    foreign "Flxq_trace" (gen @-> gen @-> pari_ulong @-> returning pari_ulong)

  let flxqc_flxq_mul =
    foreign "FlxqC_Flxq_mul"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_flxq_mul =
    foreign "FlxqM_Flxq_mul"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqv_dotproduct =
    foreign "FlxqV_dotproduct"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqv_dotproduct_pre =
    foreign "FlxqV_dotproduct_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let rg_to_f2 = foreign "Rg_to_F2" (gen @-> returning pari_ulong)
  let rg_to_fl = foreign "Rg_to_Fl" (gen @-> pari_ulong @-> returning pari_ulong)

  let rg_to_flxq =
    foreign "Rg_to_Flxq" (gen @-> gen @-> pari_ulong @-> returning gen)

  let rgx_to_flx = foreign "RgX_to_Flx" (gen @-> pari_ulong @-> returning gen)

  let rgxv_to_flxv =
    foreign "RgXV_to_FlxV" (gen @-> pari_ulong @-> returning gen)

  let z_to_flx =
    foreign "Z_to_Flx" (gen @-> pari_ulong @-> long @-> returning gen)

  let zx_to_flx = foreign "ZX_to_Flx" (gen @-> pari_ulong @-> returning gen)
  let zxv_to_flxv = foreign "ZXV_to_FlxV" (gen @-> pari_ulong @-> returning gen)
  let zxt_to_flxt = foreign "ZXT_to_FlxT" (gen @-> pari_ulong @-> returning gen)

  let gener_flxq =
    foreign "gener_Flxq" (gen @-> pari_ulong @-> ptr gen @-> returning gen)

  let get_flxq_field =
    foreign "get_Flxq_field"
      (ptr (ptr void) @-> gen @-> pari_ulong @-> returning (ptr bb_field))

  let get_flxq_star =
    foreign "get_Flxq_star"
      (ptr (ptr void) @-> gen @-> pari_ulong @-> returning (ptr bb_group))

  let monomial_flx =
    foreign "monomial_Flx" (pari_ulong @-> long @-> long @-> returning gen)

  let random_flx =
    foreign "random_Flx" (long @-> long @-> pari_ulong @-> returning gen)

  let zero_flxc = foreign "zero_FlxC" (long @-> long @-> returning gen)
  let zero_flxm = foreign "zero_FlxM" (long @-> long @-> long @-> returning gen)

  let zlx_translate1 =
    foreign "zlx_translate1" (gen @-> pari_ulong @-> long @-> returning gen)

  let zx_to_flx = foreign "zx_to_Flx" (gen @-> pari_ulong @-> returning gen)
  let flxxc_to_zxxc = foreign "FlxXC_to_ZXXC" (gen @-> returning gen)
  let flxxm_to_zxxm = foreign "FlxXM_to_ZXXM" (gen @-> returning gen)
  let flxx_to_zxx = foreign "FlxX_to_ZXX" (gen @-> returning gen)
  let fly_to_flxy = foreign "Fly_to_FlxY" (gen @-> long @-> returning gen)
  let zxx_to_flxx = foreign "zxX_to_FlxX" (gen @-> pari_ulong @-> returning gen)

  let zxx_to_kronecker =
    foreign "zxX_to_Kronecker" (gen @-> gen @-> returning gen)

  let flxqx_gcd =
    foreign "FlxqX_gcd" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_gcd_pre =
    foreign "FlxqX_gcd_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_safegcd =
    foreign "FlxqX_safegcd"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxx_add =
    foreign "FlxX_add" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxx_sub =
    foreign "FlxX_sub" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_fromnewton =
    foreign "FlxqX_fromNewton" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_fromnewton_pre =
    foreign "FlxqX_fromNewton_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_flxq_mul =
    foreign "FlxqX_Flxq_mul"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_flxq_mul_to_monic =
    foreign "FlxqX_Flxq_mul_to_monic"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_flxq_mul_pre =
    foreign "FlxqX_Flxq_mul_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_flxq_mul_to_monic_pre =
    foreign "FlxqX_Flxq_mul_to_monic_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_flxqxqv_eval =
    foreign "FlxqX_FlxqXQV_eval"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_flxqxqv_eval_pre =
    foreign "FlxqX_FlxqXQV_eval_pre"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxy_flx_translate =
    foreign "FlxY_Flx_translate" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxy_flxqv_evalx =
    foreign "FlxY_FlxqV_evalx"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxy_flxqv_evalx_pre =
    foreign "FlxY_FlxqV_evalx_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxy_flxq_evalx =
    foreign "FlxY_Flxq_evalx"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxy_flxq_evalx_pre =
    foreign "FlxY_Flxq_evalx_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_newton =
    foreign "FlxqX_Newton"
      (gen @-> long @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_newton_pre =
    foreign "FlxqX_Newton_pre"
      (gen @-> long @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxx_blocks =
    foreign "FlxX_blocks" (gen @-> long @-> long @-> long @-> returning gen)

  let zlxx_translate1 =
    foreign "zlxX_translate1" (gen @-> long @-> long @-> long @-> returning gen)

  let flxx_translate1 =
    foreign "FlxX_translate1" (gen @-> long @-> long @-> returning gen)

  let flxqx_flxqxq_eval =
    foreign "FlxqX_FlxqXQ_eval"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_flxqxq_eval_pre =
    foreign "FlxqX_FlxqXQ_eval_pre"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxy_evalx =
    foreign "FlxY_evalx" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxy_evalx_pre =
    foreign "FlxY_evalx_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_get_red =
    foreign "FlxqX_get_red" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_get_red_pre =
    foreign "FlxqX_get_red_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_invbarrett =
    foreign "FlxqX_invBarrett" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_invbarrett_pre =
    foreign "FlxqX_invBarrett_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqxv_prod =
    foreign "FlxqXV_prod" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqv_roots_to_pol =
    foreign "FlxqV_roots_to_pol"
      (gen @-> gen @-> pari_ulong @-> long @-> returning gen)

  let flxqx_powu =
    foreign "FlxqX_powu"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> returning gen)
end

module F5 (F : Ctypes.FOREIGN) = struct
  open F

  let flxqx_powu_pre =
    foreign "FlxqX_powu_pre"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxqx_saferesultant =
    foreign "FlxqX_saferesultant"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_extgcd =
    foreign "FlxqX_extgcd"
      (gen @-> gen @-> gen @-> pari_ulong @-> ptr gen @-> ptr gen
     @-> returning gen)

  let flxqxn_mul =
    foreign "FlxqXn_mul"
      (gen @-> gen @-> long @-> gen @-> pari_ulong @-> returning gen)

  let flxqxn_mul_pre =
    foreign "FlxqXn_mul_pre"
      (gen @-> gen @-> long @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxxn_red = foreign "FlxXn_red" (gen @-> long @-> returning gen)

  let flxqxn_sqr =
    foreign "FlxqXn_sqr" (gen @-> long @-> gen @-> pari_ulong @-> returning gen)

  let flxqxn_sqr_pre =
    foreign "FlxqXn_sqr_pre"
      (gen @-> long @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxx_shift = foreign "FlxX_shift" (gen @-> long @-> long @-> returning gen)

  let flxqxq_autsum =
    foreign "FlxqXQ_autsum"
      (gen @-> long @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_autsum_pre =
    foreign "FlxqXQ_autsum_pre"
      (gen @-> long @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxy_degreex = foreign "FlxY_degreex" (gen @-> returning long)

  let get_flxqxq_algebra =
    foreign "get_FlxqXQ_algebra"
      (ptr (ptr void)
      @-> gen @-> gen @-> pari_ulong
      @-> returning (ptr bb_algebra))

  let random_flxqx =
    foreign "random_FlxqX"
      (long @-> long @-> gen @-> pari_ulong @-> returning gen)

  let flxx_to_flx = foreign "FlxX_to_Flx" (gen @-> returning gen)

  let flxqxn_inv =
    foreign "FlxqXn_inv" (gen @-> long @-> gen @-> pari_ulong @-> returning gen)

  let flxqxn_inv_pre =
    foreign "FlxqXn_inv_pre"
      (gen @-> long @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqxn_expint =
    foreign "FlxqXn_expint"
      (gen @-> long @-> gen @-> pari_ulong @-> returning gen)

  let flxqxn_expint_pre =
    foreign "FlxqXn_expint_pre"
      (gen @-> long @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxy_eval_powers_pre =
    foreign "FlxY_eval_powers_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let flxy_evalx_powers_pre =
    foreign "FlxY_evalx_powers_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxy_evalx_powers_pre =
    foreign "FlxY_evalx_powers_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxx_to_flm = foreign "FlxX_to_Flm" (gen @-> long @-> returning gen)

  let flxxv_to_flxm =
    foreign "FlxXV_to_FlxM" (gen @-> long @-> long @-> returning gen)

  let pol1_flxx = foreign "pol1_FlxX" (long @-> long @-> returning gen)
  let polx_flxx = foreign "polx_FlxX" (long @-> long @-> returning gen)

  let flxqxc_flxqxq_eval =
    foreign "FlxqXC_FlxqXQ_eval"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_inv =
    foreign "FlxqXQ_inv" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_invsafe =
    foreign "FlxqXQ_invsafe"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_minpoly =
    foreign "FlxqXQ_minpoly"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_minpoly_pre =
    foreign "FlxqXQ_minpoly_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqxq_sqr =
    foreign "FlxqXQ_sqr" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_divrem_pre =
    foreign "FlxqX_divrem_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> long @-> ptr gen @-> returning gen)

  let flxqxq_inv_pre =
    foreign "FlxqXQ_inv_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqxq_invsafe_pre =
    foreign "FlxqXQ_invsafe_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqxq_sqr_pre =
    foreign "FlxqXQ_sqr_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_rem_pre =
    foreign "FlxqX_rem_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_disc =
    foreign "FlxqX_disc" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_sqr =
    foreign "FlxqX_sqr" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_sqr_pre =
    foreign "FlxqX_sqr_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let rgx_to_flxqx =
    foreign "RgX_to_FlxqX" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxyqq_pow =
    foreign "FlxYqq_pow"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_pow =
    foreign "FlxqXQ_pow"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_pow_pre =
    foreign "FlxqXQ_pow_pre"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxqxc_flxqxqv_eval =
    foreign "FlxqXC_FlxqXQV_eval"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxc_flxqxqv_eval_pre =
    foreign "FlxqXC_FlxqXQV_eval_pre"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxqxq_div =
    foreign "FlxqXQ_div"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_mul =
    foreign "FlxqXQ_mul"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_div_pre =
    foreign "FlxqXQ_div_pre"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxqxq_mul_pre =
    foreign "FlxqXQ_mul_pre"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxqx_dotproduct =
    foreign "FlxqX_dotproduct"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_halfgcd =
    foreign "FlxqX_halfgcd"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_mul =
    foreign "FlxqX_mul" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_rem =
    foreign "FlxqX_rem" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_resultant =
    foreign "FlxqX_resultant"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_divrem =
    foreign "FlxqX_divrem"
      (gen @-> gen @-> gen @-> pari_ulong @-> ptr gen @-> returning gen)

  let flxqx_halfgcd_pre =
    foreign "FlxqX_halfgcd_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_mul_pre =
    foreign "FlxqX_mul_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_extgcd_pre =
    foreign "FlxqX_extgcd_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> ptr gen @-> ptr gen
     @-> returning gen)

  let flxxc_sub =
    foreign "FlxXC_sub" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxx_flx_mul =
    foreign "FlxX_Flx_mul" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxy_flx_div =
    foreign "FlxY_Flx_div" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxx_to_flxc =
    foreign "FlxX_to_FlxC" (gen @-> long @-> long @-> returning gen)

  let flxqxq_powers_pre =
    foreign "FlxqXQ_powers_pre"
      (gen @-> long @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxx_renormalize =
    foreign "FlxX_renormalize" (gen @-> long @-> returning gen)

  let flxqxq_autpow =
    foreign "FlxqXQ_autpow"
      (gen @-> long @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_autpow_pre =
    foreign "FlxqXQ_autpow_pre"
      (gen @-> long @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxqxq_powers =
    foreign "FlxqXQ_powers"
      (gen @-> long @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_matrix_pow =
    foreign "FlxqXQ_matrix_pow"
      (gen @-> long @-> long @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxx_swap = foreign "FlxX_swap" (gen @-> long @-> long @-> returning gen)

  let flxqxq_auttrace =
    foreign "FlxqXQ_auttrace"
      (gen @-> pari_ulong @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_auttrace_pre =
    foreign "FlxqXQ_auttrace_pre"
      (gen @-> pari_ulong @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxqxq_powu =
    foreign "FlxqXQ_powu"
      (gen @-> pari_ulong @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_powu_pre =
    foreign "FlxqXQ_powu_pre"
      (gen @-> pari_ulong @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let flxx_laplace =
    foreign "FlxX_Laplace" (gen @-> pari_ulong @-> returning gen)

  let flxx_double = foreign "FlxX_double" (gen @-> pari_ulong @-> returning gen)

  let flxx_invlaplace =
    foreign "FlxX_invLaplace" (gen @-> pari_ulong @-> returning gen)

  let flxx_neg = foreign "FlxX_neg" (gen @-> pari_ulong @-> returning gen)
  let flxx_triple = foreign "FlxX_triple" (gen @-> pari_ulong @-> returning gen)

  let flxx_fl_mul =
    foreign "FlxX_Fl_mul" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxx_flx_add =
    foreign "FlxX_Flx_add" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxx_flx_sub =
    foreign "FlxX_Flx_sub" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_normalize =
    foreign "FlxqX_normalize" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_red =
    foreign "FlxqX_red" (gen @-> gen @-> pari_ulong @-> returning gen)

  let kronecker_to_flxqx =
    foreign "Kronecker_to_FlxqX" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_normalize_pre =
    foreign "FlxqX_normalize_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_red_pre =
    foreign "FlxqX_red_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let kronecker_to_flxqx_pre =
    foreign "Kronecker_to_FlxqX_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxx_deriv = foreign "FlxX_deriv" (gen @-> pari_ulong @-> returning gen)

  let flxq_ellcard =
    foreign "Flxq_ellcard" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_ellgens =
    foreign "Flxq_ellgens"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> pari_ulong
     @-> returning gen)

  let flxq_ellgroup =
    foreign "Flxq_ellgroup"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> ptr gen @-> returning gen)

  let flxq_elltwist =
    foreign "Flxq_elltwist"
      (gen @-> gen @-> gen @-> pari_ulong @-> ptr gen @-> ptr gen
     @-> returning void)

  let flxq_ellj =
    foreign "Flxq_ellj" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_ellj_to_a4a6 =
    foreign "Flxq_ellj_to_a4a6"
      (gen @-> gen @-> pari_ulong @-> ptr gen @-> ptr gen @-> returning void)

  let flxqe_add =
    foreign "FlxqE_add"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_changepoint =
    foreign "FlxqE_changepoint"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_changepointinv =
    foreign "FlxqE_changepointinv"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_dbl =
    foreign "FlxqE_dbl" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_log =
    foreign "FlxqE_log"
      (gen @-> gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)
end

module F6 (F : Ctypes.FOREIGN) = struct
  open F

  let flxqe_mul =
    foreign "FlxqE_mul"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_neg =
    foreign "FlxqE_neg" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_order =
    foreign "FlxqE_order"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_sub =
    foreign "FlxqE_sub"
      (gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_tatepairing =
    foreign "FlxqE_tatepairing"
      (gen @-> gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_weilpairing =
    foreign "FlxqE_weilpairing"
      (gen @-> gen @-> gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqe_weilpairing_pre =
    foreign "FlxqE_weilpairing_pre"
      (gen @-> gen @-> gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong
     @-> returning gen)

  let zxx_to_flxx =
    foreign "ZXX_to_FlxX" (gen @-> pari_ulong @-> long @-> returning gen)

  let zxxt_to_flxxt =
    foreign "ZXXT_to_FlxXT" (gen @-> pari_ulong @-> long @-> returning gen)

  let zxxv_to_flxxv =
    foreign "ZXXV_to_FlxXV" (gen @-> pari_ulong @-> long @-> returning gen)

  let get_flxqe_group =
    foreign "get_FlxqE_group"
      (ptr (ptr void)
      @-> gen @-> gen @-> gen @-> pari_ulong
      @-> returning (ptr bb_group))

  let rge_to_flxqe =
    foreign "RgE_to_FlxqE" (gen @-> gen @-> pari_ulong @-> returning gen)

  let random_flxqe =
    foreign "random_FlxqE" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let fl_elltrace =
    foreign "Fl_elltrace"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning long)

  let fl_elltrace_cm =
    foreign "Fl_elltrace_CM"
      (long @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning long)

  let fp_ellcard = foreign "Fp_ellcard" (gen @-> gen @-> gen @-> returning gen)

  let fp_elldivpol =
    foreign "Fp_elldivpol" (gen @-> gen @-> long @-> gen @-> returning gen)

  let fp_ellgens =
    foreign "Fp_ellgens"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fp_ellgroup =
    foreign "Fp_ellgroup"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fp_ellj = foreign "Fp_ellj" (gen @-> gen @-> gen @-> returning gen)

  let fp_elljissupersingular =
    foreign "Fp_elljissupersingular" (gen @-> gen @-> returning int)

  let fp_elltwist =
    foreign "Fp_elltwist"
      (gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning void)

  let fp_ffellcard =
    foreign "Fp_ffellcard"
      (gen @-> gen @-> gen @-> long @-> gen @-> returning gen)

  let fpe_add = foreign "FpE_add" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpe_changepoint =
    foreign "FpE_changepoint" (gen @-> gen @-> gen @-> returning gen)

  let fpe_changepointinv =
    foreign "FpE_changepointinv" (gen @-> gen @-> gen @-> returning gen)

  let fpe_dbl = foreign "FpE_dbl" (gen @-> gen @-> gen @-> returning gen)

  let fpe_log =
    foreign "FpE_log" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpe_mul = foreign "FpE_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fpe_neg = foreign "FpE_neg" (gen @-> gen @-> returning gen)

  let fpe_order =
    foreign "FpE_order" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpe_sub = foreign "FpE_sub" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fpe_to_fpj = foreign "FpE_to_FpJ" (gen @-> returning gen)
  let fpe_to_mod = foreign "FpE_to_mod" (gen @-> gen @-> returning gen)

  let fpe_tatepairing =
    foreign "FpE_tatepairing"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpe_weilpairing =
    foreign "FpE_weilpairing"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpj_add = foreign "FpJ_add" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fpj_dbl = foreign "FpJ_dbl" (gen @-> gen @-> gen @-> returning gen)
  let fpj_mul = foreign "FpJ_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fpj_neg = foreign "FpJ_neg" (gen @-> gen @-> returning gen)
  let fpj_to_fpe = foreign "FpJ_to_FpE" (gen @-> gen @-> returning gen)

  let fpxq_ellcard =
    foreign "FpXQ_ellcard" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_elldivpol =
    foreign "FpXQ_elldivpol"
      (gen @-> gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxq_ellgens =
    foreign "FpXQ_ellgens"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_ellgroup =
    foreign "FpXQ_ellgroup"
      (gen @-> gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fpxq_ellj =
    foreign "FpXQ_ellj" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_elljissupersingular =
    foreign "FpXQ_elljissupersingular" (gen @-> gen @-> gen @-> returning int)

  let fpxq_elltwist =
    foreign "FpXQ_elltwist"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning void)

  let fpxqe_add =
    foreign "FpXQE_add" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqe_changepoint =
    foreign "FpXQE_changepoint" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqe_changepointinv =
    foreign "FpXQE_changepointinv"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqe_dbl =
    foreign "FpXQE_dbl" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqe_log =
    foreign "FpXQE_log"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqe_mul =
    foreign "FpXQE_mul" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqe_neg = foreign "FpXQE_neg" (gen @-> gen @-> gen @-> returning gen)

  let fpxqe_order =
    foreign "FpXQE_order" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqe_sub =
    foreign "FpXQE_sub" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqe_tatepairing =
    foreign "FpXQE_tatepairing"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqe_weilpairing =
    foreign "FpXQE_weilpairing"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let rge_to_fpe = foreign "RgE_to_FpE" (gen @-> gen @-> returning gen)

  let rge_to_fpxqe =
    foreign "RgE_to_FpXQE" (gen @-> gen @-> gen @-> returning gen)

  let get_fpe_group =
    foreign "get_FpE_group"
      (ptr (ptr void) @-> gen @-> gen @-> gen @-> returning (ptr bb_group))

  let get_fpxqe_group =
    foreign "get_FpXQE_group"
      (ptr (ptr void)
      @-> gen @-> gen @-> gen @-> gen
      @-> returning (ptr bb_group))

  let elltrace_extension =
    foreign "elltrace_extension" (gen @-> long @-> gen @-> returning gen)

  let random_fpe = foreign "random_FpE" (gen @-> gen @-> gen @-> returning gen)

  let random_fpxqe =
    foreign "random_FpXQE" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fp_issquare = foreign "Fp_issquare" (gen @-> gen @-> returning int)
  let fp_fpx_sub = foreign "Fp_FpX_sub" (gen @-> gen @-> gen @-> returning gen)

  let fp_fpxq_log =
    foreign "Fp_FpXQ_log" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpv_fpm_polint =
    foreign "FpV_FpM_polint" (gen @-> gen @-> gen @-> long @-> returning gen)

  let fpv_inv = foreign "FpV_inv" (gen @-> gen @-> returning gen)

  let fpv_invvandermonde =
    foreign "FpV_invVandermonde" (gen @-> gen @-> gen @-> returning gen)

  let fpv_polint =
    foreign "FpV_polint" (gen @-> gen @-> gen @-> long @-> returning gen)

  let fpv_roots_to_pol =
    foreign "FpV_roots_to_pol" (gen @-> gen @-> long @-> returning gen)

  let fpx_fp_add = foreign "FpX_Fp_add" (gen @-> gen @-> gen @-> returning gen)

  let fpx_fp_add_shallow =
    foreign "FpX_Fp_add_shallow" (gen @-> gen @-> gen @-> returning gen)

  let fpx_fp_div = foreign "FpX_Fp_div" (gen @-> gen @-> gen @-> returning gen)
  let fpx_fp_mul = foreign "FpX_Fp_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpx_fp_mul_to_monic =
    foreign "FpX_Fp_mul_to_monic" (gen @-> gen @-> gen @-> returning gen)

  let fpx_fp_mulspec =
    foreign "FpX_Fp_mulspec" (gen @-> gen @-> gen @-> long @-> returning gen)

  let fpx_fp_sub = foreign "FpX_Fp_sub" (gen @-> gen @-> gen @-> returning gen)

  let fpx_fp_sub_shallow =
    foreign "FpX_Fp_sub_shallow" (gen @-> gen @-> gen @-> returning gen)

  let fpx_fpv_multieval =
    foreign "FpX_FpV_multieval" (gen @-> gen @-> gen @-> returning gen)

  let fpx_fpxq_eval =
    foreign "FpX_FpXQ_eval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpx_fpxqv_eval =
    foreign "FpX_FpXQV_eval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpx_fpxv_multirem =
    foreign "FpX_FpXV_multirem" (gen @-> gen @-> gen @-> returning gen)

  let fpx_frobenius = foreign "FpX_Frobenius" (gen @-> gen @-> returning gen)
  let fpx_laplace = foreign "FpX_Laplace" (gen @-> gen @-> returning gen)
  let fpx_newton = foreign "FpX_Newton" (gen @-> long @-> gen @-> returning gen)
  let fpx_add = foreign "FpX_add" (gen @-> gen @-> gen @-> returning gen)
  let fpx_center = foreign "FpX_center" (gen @-> gen @-> gen @-> returning gen)

  let fpx_center_i =
    foreign "FpX_center_i" (gen @-> gen @-> gen @-> returning gen)

  let fpx_chinese_coprime =
    foreign "FpX_chinese_coprime"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpx_convol = foreign "FpX_convol" (gen @-> gen @-> gen @-> returning gen)
  let fpx_deriv = foreign "FpX_deriv" (gen @-> gen @-> returning gen)
  let fpx_digits = foreign "FpX_digits" (gen @-> gen @-> gen @-> returning gen)
  let fpx_disc = foreign "FpX_disc" (gen @-> gen @-> returning gen)

  let fpx_div_by_x_x =
    foreign "FpX_div_by_X_x" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fpx_divrem =
    foreign "FpX_divrem" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fpx_divu =
    foreign "FpX_divu" (gen @-> pari_ulong @-> gen @-> returning gen)
end

module F7 (F : Ctypes.FOREIGN) = struct
  open F

  let fpx_dotproduct =
    foreign "FpX_dotproduct" (gen @-> gen @-> gen @-> returning gen)

  let fpx_eval = foreign "FpX_eval" (gen @-> gen @-> gen @-> returning gen)

  let fpx_extgcd =
    foreign "FpX_extgcd"
      (gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let fpx_fromnewton = foreign "FpX_fromNewton" (gen @-> gen @-> returning gen)
  let fpx_gcd = foreign "FpX_gcd" (gen @-> gen @-> gen @-> returning gen)

  let fpx_gcd_check =
    foreign "FpX_gcd_check" (gen @-> gen @-> gen @-> returning gen)

  let fpx_get_red = foreign "FpX_get_red" (gen @-> gen @-> returning gen)
  let fpx_halve = foreign "FpX_halve" (gen @-> gen @-> returning gen)
  let fpx_halfgcd = foreign "FpX_halfgcd" (gen @-> gen @-> gen @-> returning gen)
  let fpx_integ = foreign "FpX_integ" (gen @-> gen @-> returning gen)
  let fpx_invbarrett = foreign "FpX_invBarrett" (gen @-> gen @-> returning gen)
  let fpx_invlaplace = foreign "FpX_invLaplace" (gen @-> gen @-> returning gen)

  let fpx_is_squarefree =
    foreign "FpX_is_squarefree" (gen @-> gen @-> returning int)

  let fpx_matfrobenius =
    foreign "FpX_matFrobenius" (gen @-> gen @-> returning gen)

  let fpx_mul = foreign "FpX_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpx_mulspec =
    foreign "FpX_mulspec"
      (gen @-> gen @-> gen @-> long @-> long @-> returning gen)

  let fpx_mulu =
    foreign "FpX_mulu" (gen @-> pari_ulong @-> gen @-> returning gen)

  let fpx_neg = foreign "FpX_neg" (gen @-> gen @-> returning gen)
  let fpx_normalize = foreign "FpX_normalize" (gen @-> gen @-> returning gen)

  let fpx_powu =
    foreign "FpX_powu" (gen @-> pari_ulong @-> gen @-> returning gen)

  let fpx_red = foreign "FpX_red" (gen @-> gen @-> returning gen)
  let fpx_rem = foreign "FpX_rem" (gen @-> gen @-> gen @-> returning gen)
  let fpx_rescale = foreign "FpX_rescale" (gen @-> gen @-> gen @-> returning gen)

  let fpx_resultant =
    foreign "FpX_resultant" (gen @-> gen @-> gen @-> returning gen)

  let fpx_sqr = foreign "FpX_sqr" (gen @-> gen @-> returning gen)
  let fpx_sub = foreign "FpX_sub" (gen @-> gen @-> gen @-> returning gen)

  let fpx_valrem =
    foreign "FpX_valrem" (gen @-> gen @-> gen @-> ptr gen @-> returning long)

  let fpxc_fpxq_eval =
    foreign "FpXC_FpXQ_eval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxc_fpxqv_eval =
    foreign "FpXC_FpXQV_eval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxm_fpxqv_eval =
    foreign "FpXM_FpXQV_eval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_autpow =
    foreign "FpXQ_autpow" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fpxq_autpowers =
    foreign "FpXQ_autpowers" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxq_autsum =
    foreign "FpXQ_autsum" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fpxq_auttrace =
    foreign "FpXQ_auttrace"
      (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fpxq_charpoly =
    foreign "FpXQ_charpoly" (gen @-> gen @-> gen @-> returning gen)

  let fpxq_conjvec =
    foreign "FpXQ_conjvec" (gen @-> gen @-> gen @-> returning gen)

  let fpxq_div =
    foreign "FpXQ_div" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_inv = foreign "FpXQ_inv" (gen @-> gen @-> gen @-> returning gen)

  let fpxq_invsafe =
    foreign "FpXQ_invsafe" (gen @-> gen @-> gen @-> returning gen)

  let fpxq_issquare =
    foreign "FpXQ_issquare" (gen @-> gen @-> gen @-> returning int)

  let fpxq_log =
    foreign "FpXQ_log" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_matrix_pow =
    foreign "FpXQ_matrix_pow"
      (gen @-> long @-> long @-> gen @-> gen @-> returning gen)

  let fpxq_minpoly =
    foreign "FpXQ_minpoly" (gen @-> gen @-> gen @-> returning gen)

  let fpxq_mul =
    foreign "FpXQ_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_norm = foreign "FpXQ_norm" (gen @-> gen @-> gen @-> returning gen)

  let fpxq_order =
    foreign "FpXQ_order" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_pow =
    foreign "FpXQ_pow" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_powu =
    foreign "FpXQ_powu" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fpxq_powers =
    foreign "FpXQ_powers" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxq_red = foreign "FpXQ_red" (gen @-> gen @-> gen @-> returning gen)
  let fpxq_sqr = foreign "FpXQ_sqr" (gen @-> gen @-> gen @-> returning gen)
  let fpxq_sqrt = foreign "FpXQ_sqrt" (gen @-> gen @-> gen @-> returning gen)

  let fpxq_sqrtn =
    foreign "FpXQ_sqrtn"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fpxq_trace = foreign "FpXQ_trace" (gen @-> gen @-> gen @-> returning gen)

  let fpxqc_to_mod =
    foreign "FpXQC_to_mod" (gen @-> gen @-> gen @-> returning gen)

  let fpxqm_autsum =
    foreign "FpXQM_autsum" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fpxt_red = foreign "FpXT_red" (gen @-> gen @-> returning gen)

  let fpxv_fpx_fromdigits =
    foreign "FpXV_FpX_fromdigits" (gen @-> gen @-> gen @-> returning gen)

  let fpxv_chinese =
    foreign "FpXV_chinese" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fpxv_factorback =
    foreign "FpXV_factorback" (gen @-> gen @-> gen @-> long @-> returning gen)

  let fpxv_prod = foreign "FpXV_prod" (gen @-> gen @-> returning gen)
  let fpxv_red = foreign "FpXV_red" (gen @-> gen @-> returning gen)

  let fpxn_div =
    foreign "FpXn_div" (gen @-> gen @-> long @-> gen @-> returning gen)

  let fpxn_exp = foreign "FpXn_exp" (gen @-> long @-> gen @-> returning gen)

  let fpxn_expint =
    foreign "FpXn_expint" (gen @-> long @-> gen @-> returning gen)

  let fpxn_inv = foreign "FpXn_inv" (gen @-> long @-> gen @-> returning gen)

  let fpxn_mul =
    foreign "FpXn_mul" (gen @-> gen @-> long @-> gen @-> returning gen)

  let fpxn_sqr = foreign "FpXn_sqr" (gen @-> long @-> gen @-> returning gen)
  let fq_issquare = foreign "Fq_issquare" (gen @-> gen @-> gen @-> returning int)

  let fq_ispower =
    foreign "Fq_ispower" (gen @-> gen @-> gen @-> gen @-> returning long)

  let fq_log =
    foreign "Fq_log" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqc_to_mod = foreign "FqC_to_mod" (gen @-> gen @-> gen @-> returning gen)
  let fqm_to_mod = foreign "FqM_to_mod" (gen @-> gen @-> gen @-> returning gen)
  let fqv_inv = foreign "FqV_inv" (gen @-> gen @-> gen @-> returning gen)
  let z_to_fpx = foreign "Z_to_FpX" (gen @-> gen @-> long @-> returning gen)

  let gener_fpxq =
    foreign "gener_FpXQ" (gen @-> gen @-> ptr gen @-> returning gen)

  let gener_fpxq_local =
    foreign "gener_FpXQ_local" (gen @-> gen @-> gen @-> returning gen)

  let get_fpxq_star =
    foreign "get_FpXQ_star"
      (ptr (ptr void) @-> gen @-> gen @-> returning (ptr bb_group))

  let get_fpx_algebra =
    foreign "get_FpX_algebra"
      (ptr (ptr void) @-> gen @-> long @-> returning (ptr bb_algebra))

  let get_fpxq_algebra =
    foreign "get_FpXQ_algebra"
      (ptr (ptr void) @-> gen @-> gen @-> returning (ptr bb_algebra))

  let random_fpx = foreign "random_FpX" (long @-> long @-> gen @-> returning gen)
  let f2x_ddf = foreign "F2x_ddf" (gen @-> returning gen)
  let f2x_factor = foreign "F2x_factor" (gen @-> returning gen)

  let f2x_factor_squarefree =
    foreign "F2x_factor_squarefree" (gen @-> returning gen)

  let f2x_is_irred = foreign "F2x_is_irred" (gen @-> returning int)
  let flx_ddf = foreign "Flx_ddf" (gen @-> pari_ulong @-> returning gen)

  let flx_ddf_pre =
    foreign "Flx_ddf_pre" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_is_irred =
    foreign "Flx_is_irred" (gen @-> pari_ulong @-> returning int)

  let flx_is_totally_split =
    foreign "Flx_is_totally_split" (gen @-> pari_ulong @-> returning int)

  let flx_ispower =
    foreign "Flx_ispower"
      (gen @-> pari_ulong @-> pari_ulong @-> ptr gen @-> returning long)

  let flx_degfact = foreign "Flx_degfact" (gen @-> pari_ulong @-> returning gen)
  let flx_factor = foreign "Flx_factor" (gen @-> pari_ulong @-> returning gen)

  let flx_factor_squarefree =
    foreign "Flx_factor_squarefree" (gen @-> pari_ulong @-> returning gen)

  let flx_factor_squarefree_pre =
    foreign "Flx_factor_squarefree_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_nbfact = foreign "Flx_nbfact" (gen @-> pari_ulong @-> returning long)

  let flx_nbfact_pre =
    foreign "Flx_nbfact_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> returning long)

  let flx_nbfact_frobenius =
    foreign "Flx_nbfact_Frobenius"
      (gen @-> gen @-> pari_ulong @-> returning long)

  let flx_nbfact_frobenius_pre =
    foreign "Flx_nbfact_Frobenius_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning long)

  let flx_nbfact_by_degree =
    foreign "Flx_nbfact_by_degree"
      (gen @-> ptr long @-> pari_ulong @-> returning gen)

  let flx_nbroots = foreign "Flx_nbroots" (gen @-> pari_ulong @-> returning long)
end

module F8 (F : Ctypes.FOREIGN) = struct
  open F

  let flx_oneroot =
    foreign "Flx_oneroot" (gen @-> pari_ulong @-> returning pari_ulong)

  let flx_oneroot_split =
    foreign "Flx_oneroot_split" (gen @-> pari_ulong @-> returning pari_ulong)

  let flx_oneroot_pre =
    foreign "Flx_oneroot_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let flx_oneroot_split_pre =
    foreign "Flx_oneroot_split_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let flx_roots = foreign "Flx_roots" (gen @-> pari_ulong @-> returning gen)

  let flx_roots_pre =
    foreign "Flx_roots_pre" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_rootsff =
    foreign "Flx_rootsff" (gen @-> gen @-> pari_ulong @-> returning gen)

  let fpx_ddf = foreign "FpX_ddf" (gen @-> gen @-> returning gen)

  let fpx_ddf_degree =
    foreign "FpX_ddf_degree" (gen @-> gen @-> gen @-> returning long)

  let fpx_degfact = foreign "FpX_degfact" (gen @-> gen @-> returning gen)
  let fpx_factor = foreign "FpX_factor" (gen @-> gen @-> returning gen)

  let fpx_factor_squarefree =
    foreign "FpX_factor_squarefree" (gen @-> gen @-> returning gen)

  let fpx_is_irred = foreign "FpX_is_irred" (gen @-> gen @-> returning int)

  let fpx_is_totally_split =
    foreign "FpX_is_totally_split" (gen @-> gen @-> returning int)

  let fpx_ispower =
    foreign "FpX_ispower"
      (gen @-> pari_ulong @-> gen @-> ptr gen @-> returning long)

  let fpx_nbfact = foreign "FpX_nbfact" (gen @-> gen @-> returning long)

  let fpx_nbfact_frobenius =
    foreign "FpX_nbfact_Frobenius" (gen @-> gen @-> gen @-> returning long)

  let fpx_nbroots = foreign "FpX_nbroots" (gen @-> gen @-> returning long)
  let fpx_oneroot = foreign "FpX_oneroot" (gen @-> gen @-> returning gen)

  let fpx_oneroot_split =
    foreign "FpX_oneroot_split" (gen @-> gen @-> returning gen)

  let fpx_roots = foreign "FpX_roots" (gen @-> gen @-> returning gen)

  let fpx_roots_mult =
    foreign "FpX_roots_mult" (gen @-> long @-> gen @-> returning gen)

  let fpx_rootsff = foreign "FpX_rootsff" (gen @-> gen @-> gen @-> returning gen)
  let fpx_split_part = foreign "FpX_split_part" (gen @-> gen @-> returning gen)
  let f2xqx_ddf = foreign "F2xqX_ddf" (gen @-> gen @-> returning gen)
  let f2xqx_degfact = foreign "F2xqX_degfact" (gen @-> gen @-> returning gen)
  let f2xqx_factor = foreign "F2xqX_factor" (gen @-> gen @-> returning gen)

  let f2xqx_factor_squarefree =
    foreign "F2xqX_factor_squarefree" (gen @-> gen @-> returning gen)

  let f2xqx_roots = foreign "F2xqX_roots" (gen @-> gen @-> returning gen)

  let flx_factorff_irred =
    foreign "Flx_factorff_irred" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_ffintersect =
    foreign "Flx_ffintersect"
      (gen @-> gen @-> long @-> pari_ulong @-> ptr gen @-> ptr gen @-> gen
     @-> gen @-> returning void)

  let flx_ffisom =
    foreign "Flx_ffisom" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_ffisom_inv =
    foreign "Flxq_ffisom_inv" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_frobenius =
    foreign "FlxqX_Frobenius" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_frobenius_pre =
    foreign "FlxqX_Frobenius_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_ddf =
    foreign "FlxqX_ddf" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_ddf_degree =
    foreign "FlxqX_ddf_degree"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning long)

  let flxqx_degfact =
    foreign "FlxqX_degfact" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_factor =
    foreign "FlxqX_factor" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_factor_squarefree =
    foreign "FlxqX_factor_squarefree"
      (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_factor_squarefree_pre =
    foreign "FlxqX_factor_squarefree_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flxqx_ispower =
    foreign "FlxqX_ispower"
      (gen @-> pari_ulong @-> gen @-> pari_ulong @-> ptr gen @-> returning long)

  let flxqx_is_squarefree =
    foreign "FlxqX_is_squarefree" (gen @-> gen @-> pari_ulong @-> returning long)

  let flxqx_nbfact =
    foreign "FlxqX_nbfact" (gen @-> gen @-> pari_ulong @-> returning long)

  let flxqx_nbfact_frobenius =
    foreign "FlxqX_nbfact_Frobenius"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning long)

  let flxqx_nbfact_by_degree =
    foreign "FlxqX_nbfact_by_degree"
      (gen @-> ptr long @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_nbroots =
    foreign "FlxqX_nbroots" (gen @-> gen @-> pari_ulong @-> returning long)

  let flxqx_roots =
    foreign "FlxqX_roots" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqxq_halffrobenius =
    foreign "FlxqXQ_halfFrobenius"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let fpx_factorff =
    foreign "FpX_factorff" (gen @-> gen @-> gen @-> returning gen)

  let fpx_factorff_irred =
    foreign "FpX_factorff_irred" (gen @-> gen @-> gen @-> returning gen)

  let fpx_ffintersect =
    foreign "FpX_ffintersect"
      (gen @-> gen @-> long @-> gen @-> ptr gen @-> ptr gen @-> gen @-> gen
     @-> returning void)

  let fpx_ffisom = foreign "FpX_ffisom" (gen @-> gen @-> gen @-> returning gen)

  let fpxq_ffisom_inv =
    foreign "FpXQ_ffisom_inv" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_frobenius =
    foreign "FpXQX_Frobenius" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_ddf = foreign "FpXQX_ddf" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_ddf_degree =
    foreign "FpXQX_ddf_degree" (gen @-> gen @-> gen @-> gen @-> returning long)

  let fpxqx_degfact =
    foreign "FpXQX_degfact" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_factor =
    foreign "FpXQX_factor" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_factor_squarefree =
    foreign "FpXQX_factor_squarefree" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_ispower =
    foreign "FpXQX_ispower"
      (gen @-> pari_ulong @-> gen @-> gen @-> ptr gen @-> returning long)

  let fpxqx_nbfact =
    foreign "FpXQX_nbfact" (gen @-> gen @-> gen @-> returning long)

  let fpxqx_nbfact_frobenius =
    foreign "FpXQX_nbfact_Frobenius"
      (gen @-> gen @-> gen @-> gen @-> returning long)

  let fpxqx_nbroots =
    foreign "FpXQX_nbroots" (gen @-> gen @-> gen @-> returning long)

  let fpxqx_roots = foreign "FpXQX_roots" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_split_part =
    foreign "FpXQX_split_part" (gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_halffrobenius =
    foreign "FpXQXQ_halfFrobenius"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_is_squarefree =
    foreign "FqX_is_squarefree" (gen @-> gen @-> gen @-> returning long)

  let fqx_ispower =
    foreign "FqX_ispower"
      (gen @-> pari_ulong @-> gen @-> gen @-> ptr gen @-> returning long)

  let fqx_nbfact = foreign "FqX_nbfact" (gen @-> gen @-> gen @-> returning long)

  let fqx_nbroots =
    foreign "FqX_nbroots" (gen @-> gen @-> gen @-> returning long)

  let factorff = foreign "factorff" (gen @-> gen @-> gen @-> returning gen)
  let factormod0 = foreign "factormod0" (gen @-> gen @-> long @-> returning gen)
  let factormodddf = foreign "factormodDDF" (gen @-> gen @-> returning gen)
  let factormodsqf = foreign "factormodSQF" (gen @-> gen @-> returning gen)

  let ff_parse_tp =
    foreign "ff_parse_Tp"
      (gen @-> ptr gen @-> ptr gen @-> long @-> returning int)

  let polrootsff = foreign "polrootsff" (gen @-> gen @-> gen @-> returning gen)
  let polrootsmod = foreign "polrootsmod" (gen @-> gen @-> returning gen)
  let rootmod0 = foreign "rootmod0" (gen @-> gen @-> long @-> returning gen)

  let fpxqx_fpxq_mul =
    foreign "FpXQX_FpXQ_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_fpxqxqv_eval =
    foreign "FpXQX_FpXQXQV_eval"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_fpxqxq_eval =
    foreign "FpXQX_FpXQXQ_eval"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_digits =
    foreign "FpXQX_digits" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_disc = foreign "FpXQX_disc" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_div_by_x_x =
    foreign "FpXQX_div_by_X_x"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fpxqx_divrem =
    foreign "FpXQX_divrem"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fpxqx_dotproduct =
    foreign "FpXQX_dotproduct" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_extgcd =
    foreign "FpXQX_extgcd"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let fpxqx_gcd =
    foreign "FpXQX_gcd" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_get_red =
    foreign "FpXQX_get_red" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_halfgcd =
    foreign "FpXQX_halfgcd" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_invbarrett =
    foreign "FpXQX_invBarrett" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_mul =
    foreign "FpXQX_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_powu =
    foreign "FpXQX_powu" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fpxqx_red = foreign "FpXQX_red" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_rem =
    foreign "FpXQX_rem" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_resultant =
    foreign "FpXQX_resultant" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqx_sqr = foreign "FpXQX_sqr" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_to_mod =
    foreign "FpXQX_to_mod" (gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_div =
    foreign "FpXQXQ_div" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)
end

module F9 (F : Ctypes.FOREIGN) = struct
  open F

  let fpxqxq_inv =
    foreign "FpXQXQ_inv" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_invsafe =
    foreign "FpXQXQ_invsafe" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_matrix_pow =
    foreign "FpXQXQ_matrix_pow"
      (gen @-> long @-> long @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_minpoly =
    foreign "FpXQXQ_minpoly" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_mul =
    foreign "FpXQXQ_mul" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_pow =
    foreign "FpXQXQ_pow" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_powers =
    foreign "FpXQXQ_powers"
      (gen @-> long @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_sqr =
    foreign "FpXQXQ_sqr" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_autpow =
    foreign "FpXQXQ_autpow"
      (gen @-> long @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_autsum =
    foreign "FpXQXQ_autsum"
      (gen @-> long @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxq_auttrace =
    foreign "FpXQXQ_auttrace"
      (gen @-> long @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxt_red = foreign "FpXQXT_red" (gen @-> gen @-> gen @-> returning gen)

  let fpxqxv_fpxqx_fromdigits =
    foreign "FpXQXV_FpXQX_fromdigits"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxqxv_prod = foreign "FpXQXV_prod" (gen @-> gen @-> gen @-> returning gen)
  let fpxqxv_red = foreign "FpXQXV_red" (gen @-> gen @-> gen @-> returning gen)

  let fpxqxn_div =
    foreign "FpXQXn_div" (gen @-> gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxqxn_exp =
    foreign "FpXQXn_exp" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxqxn_expint =
    foreign "FpXQXn_expint" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxqxn_inv =
    foreign "FpXQXn_inv" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxqxn_mul =
    foreign "FpXQXn_mul" (gen @-> gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxqxn_sqr =
    foreign "FpXQXn_sqr" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxx_fp_mul = foreign "FpXX_Fp_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpxx_fpx_mul =
    foreign "FpXX_FpX_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpxx_add = foreign "FpXX_add" (gen @-> gen @-> gen @-> returning gen)
  let fpxx_deriv = foreign "FpXX_deriv" (gen @-> gen @-> returning gen)
  let fpxx_halve = foreign "FpXX_halve" (gen @-> gen @-> returning gen)
  let fpxx_integ = foreign "FpXX_integ" (gen @-> gen @-> returning gen)

  let fpxx_mulu =
    foreign "FpXX_mulu" (gen @-> pari_ulong @-> gen @-> returning gen)

  let fpxx_neg = foreign "FpXX_neg" (gen @-> gen @-> returning gen)
  let fpxx_red = foreign "FpXX_red" (gen @-> gen @-> returning gen)
  let fpxx_sub = foreign "FpXX_sub" (gen @-> gen @-> gen @-> returning gen)

  let fpxy_fpxq_evalx =
    foreign "FpXY_FpXQ_evalx" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxy_fpxqv_evalx =
    foreign "FpXY_FpXQV_evalx" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxy_eval =
    foreign "FpXY_eval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxy_evalx = foreign "FpXY_evalx" (gen @-> gen @-> gen @-> returning gen)

  let fpxy_evaly =
    foreign "FpXY_evaly" (gen @-> gen @-> gen @-> long @-> returning gen)

  let fpxyqq_pow =
    foreign "FpXYQQ_pow" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxc_to_mod = foreign "FqXC_to_mod" (gen @-> gen @-> gen @-> returning gen)
  let fqxm_to_mod = foreign "FqXM_to_mod" (gen @-> gen @-> gen @-> returning gen)

  let kronecker_to_fpxqx =
    foreign "Kronecker_to_FpXQX" (gen @-> gen @-> gen @-> returning gen)

  let get_fpxqx_algebra =
    foreign "get_FpXQX_algebra"
      (ptr (ptr void) @-> gen @-> gen @-> long @-> returning (ptr bb_algebra))

  let get_fpxqxq_algebra =
    foreign "get_FpXQXQ_algebra"
      (ptr (ptr void) @-> gen @-> gen @-> gen @-> returning (ptr bb_algebra))

  let random_fpxqx =
    foreign "random_FpXQX" (long @-> long @-> gen @-> gen @-> returning gen)

  let flc_flv_mul =
    foreign "Flc_Flv_mul" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flc_to_mod = foreign "Flc_to_mod" (gen @-> pari_ulong @-> returning gen)

  let flm_fl_add =
    foreign "Flm_Fl_add" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flm_fl_mul =
    foreign "Flm_Fl_mul" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flm_fl_mul_inplace =
    foreign "Flm_Fl_mul_inplace"
      (gen @-> pari_ulong @-> pari_ulong @-> returning void)

  let flm_fl_mul_pre =
    foreign "Flm_Fl_mul_pre"
      (gen @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flm_fl_sub =
    foreign "Flm_Fl_sub" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flm_flc_mul =
    foreign "Flm_Flc_mul" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flm_flc_mul_pre =
    foreign "Flm_Flc_mul_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flm_flc_mul_pre_flx =
    foreign "Flm_Flc_mul_pre_Flx"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> long @-> returning gen)

  let flm_add = foreign "Flm_add" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flm_center =
    foreign "Flm_center" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flm_mul = foreign "Flm_mul" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flm_mul_pre =
    foreign "Flm_mul_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flm_neg = foreign "Flm_neg" (gen @-> pari_ulong @-> returning gen)

  let flm_powers =
    foreign "Flm_powers" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flm_powu =
    foreign "Flm_powu" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flm_sub = foreign "Flm_sub" (gen @-> gen @-> pari_ulong @-> returning gen)
  let flm_to_mod = foreign "Flm_to_mod" (gen @-> pari_ulong @-> returning gen)
  let flm_transpose = foreign "Flm_transpose" (gen @-> returning gen)

  let flv_fl_div =
    foreign "Flv_Fl_div" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flv_fl_div_inplace =
    foreign "Flv_Fl_div_inplace"
      (gen @-> pari_ulong @-> pari_ulong @-> returning void)

  let flv_fl_mul =
    foreign "Flv_Fl_mul" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flv_fl_mul_inplace =
    foreign "Flv_Fl_mul_inplace"
      (gen @-> pari_ulong @-> pari_ulong @-> returning void)

  let flv_fl_mul_part_inplace =
    foreign "Flv_Fl_mul_part_inplace"
      (gen @-> pari_ulong @-> pari_ulong @-> long @-> returning void)

  let flv_add = foreign "Flv_add" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flv_add_inplace =
    foreign "Flv_add_inplace" (gen @-> gen @-> pari_ulong @-> returning void)

  let flv_center =
    foreign "Flv_center" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flv_dotproduct =
    foreign "Flv_dotproduct"
      (gen @-> gen @-> pari_ulong @-> returning pari_ulong)

  let flv_dotproduct_pre =
    foreign "Flv_dotproduct_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let flv_neg = foreign "Flv_neg" (gen @-> pari_ulong @-> returning gen)

  let flv_neg_inplace =
    foreign "Flv_neg_inplace" (gen @-> pari_ulong @-> returning void)

  let flv_sub = foreign "Flv_sub" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flv_sub_inplace =
    foreign "Flv_sub_inplace" (gen @-> gen @-> pari_ulong @-> returning void)

  let flv_sum = foreign "Flv_sum" (gen @-> pari_ulong @-> returning pari_ulong)

  let flx_dotproduct =
    foreign "Flx_dotproduct"
      (gen @-> gen @-> pari_ulong @-> returning pari_ulong)

  let flx_dotproduct_pre =
    foreign "Flx_dotproduct_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fp_to_mod = foreign "Fp_to_mod" (gen @-> gen @-> returning gen)
  let fpc_fpv_mul = foreign "FpC_FpV_mul" (gen @-> gen @-> gen @-> returning gen)
  let fpc_fp_mul = foreign "FpC_Fp_mul" (gen @-> gen @-> gen @-> returning gen)
  let fpc_center = foreign "FpC_center" (gen @-> gen @-> gen @-> returning gen)

  let fpc_center_inplace =
    foreign "FpC_center_inplace" (gen @-> gen @-> gen @-> returning void)

  let fpc_red = foreign "FpC_red" (gen @-> gen @-> returning gen)
  let fpc_to_mod = foreign "FpC_to_mod" (gen @-> gen @-> returning gen)
  let fpm_add = foreign "FpM_add" (gen @-> gen @-> gen @-> returning gen)
  let fpm_fp_mul = foreign "FpM_Fp_mul" (gen @-> gen @-> gen @-> returning gen)
  let fpm_fpc_mul = foreign "FpM_FpC_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpm_fpc_mul_fpx =
    foreign "FpM_FpC_mul_FpX" (gen @-> gen @-> gen @-> long @-> returning gen)

  let fpm_center = foreign "FpM_center" (gen @-> gen @-> gen @-> returning gen)

  let fpm_center_inplace =
    foreign "FpM_center_inplace" (gen @-> gen @-> gen @-> returning void)

  let fpm_mul = foreign "FpM_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpm_powu =
    foreign "FpM_powu" (gen @-> pari_ulong @-> gen @-> returning gen)

  let fpm_red = foreign "FpM_red" (gen @-> gen @-> returning gen)
  let fpm_sub = foreign "FpM_sub" (gen @-> gen @-> gen @-> returning gen)
  let fpm_to_mod = foreign "FpM_to_mod" (gen @-> gen @-> returning gen)

  let fpms_fpc_mul =
    foreign "FpMs_FpC_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpms_fpcs_solve =
    foreign "FpMs_FpCs_solve" (gen @-> gen @-> long @-> gen @-> returning gen)
end

module F10 (F : Ctypes.FOREIGN) = struct
  open F

  let fpms_fpcs_solve_safe =
    foreign "FpMs_FpCs_solve_safe"
      (gen @-> gen @-> long @-> gen @-> returning gen)

  let fpms_leftkernel_elt =
    foreign "FpMs_leftkernel_elt" (gen @-> long @-> gen @-> returning gen)

  let fpc_add = foreign "FpC_add" (gen @-> gen @-> gen @-> returning gen)
  let fpc_sub = foreign "FpC_sub" (gen @-> gen @-> gen @-> returning gen)

  let fpv_fpms_mul =
    foreign "FpV_FpMs_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpv_add = foreign "FpV_add" (gen @-> gen @-> gen @-> returning gen)
  let fpv_sub = foreign "FpV_sub" (gen @-> gen @-> gen @-> returning gen)

  let fpv_dotproduct =
    foreign "FpV_dotproduct" (gen @-> gen @-> gen @-> returning gen)

  let fpv_dotsquare = foreign "FpV_dotsquare" (gen @-> gen @-> returning gen)
  let fpv_red = foreign "FpV_red" (gen @-> gen @-> returning gen)
  let fpv_to_mod = foreign "FpV_to_mod" (gen @-> gen @-> returning gen)
  let fpvv_to_mod = foreign "FpVV_to_mod" (gen @-> gen @-> returning gen)
  let fpx_to_mod = foreign "FpX_to_mod" (gen @-> gen @-> returning gen)
  let fpxc_to_mod = foreign "FpXC_to_mod" (gen @-> gen @-> returning gen)
  let fpxm_to_mod = foreign "FpXM_to_mod" (gen @-> gen @-> returning gen)
  let zabm_ker = foreign "ZabM_ker" (gen @-> gen @-> long @-> returning gen)

  let zabm_indexrank =
    foreign "ZabM_indexrank" (gen @-> gen @-> long @-> returning gen)

  let zabm_inv =
    foreign "ZabM_inv" (gen @-> gen @-> long @-> ptr gen @-> returning gen)

  let zabm_inv_ratlift =
    foreign "ZabM_inv_ratlift"
      (gen @-> gen @-> long @-> ptr gen @-> returning gen)

  let zabm_pseudoinv =
    foreign "ZabM_pseudoinv"
      (gen @-> gen @-> long @-> ptr gen @-> ptr gen @-> returning gen)

  let zv_zms_mul = foreign "ZV_zMs_mul" (gen @-> gen @-> returning gen)

  let zpms_zpcs_solve =
    foreign "ZpMs_ZpCs_solve"
      (gen @-> gen @-> long @-> gen @-> long @-> returning gen)

  let gen_fpm_wiedemann =
    foreign "gen_FpM_Wiedemann"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> returning gen)

  let gen_zpm_dixon_wiedemann =
    foreign "gen_ZpM_Dixon_Wiedemann"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> returning gen)

  let gen_matid =
    foreign "gen_matid" (long @-> ptr void @-> ptr bb_field @-> returning gen)

  let matid_flm = foreign "matid_Flm" (long @-> returning gen)
  let matid_f2xqm = foreign "matid_F2xqM" (long @-> gen @-> returning gen)

  let matid_flxqm =
    foreign "matid_FlxqM" (long @-> gen @-> pari_ulong @-> returning gen)

  let random_flv = foreign "random_Flv" (long @-> pari_ulong @-> returning gen)
  let random_fpc = foreign "random_FpC" (long @-> gen @-> returning gen)
  let random_fpv = foreign "random_FpV" (long @-> gen @-> returning gen)
  let scalar_flm = foreign "scalar_Flm" (long @-> long @-> returning gen)
  let zcs_to_zc = foreign "zCs_to_ZC" (gen @-> long @-> returning gen)
  let zms_to_zm = foreign "zMs_to_ZM" (gen @-> long @-> returning gen)
  let zms_zc_mul = foreign "zMs_ZC_mul" (gen @-> gen @-> returning gen)
  let zmv_to_flmv = foreign "ZMV_to_FlmV" (gen @-> pari_ulong @-> returning gen)

  let flx_teichmuller =
    foreign "Flx_Teichmuller" (gen @-> pari_ulong @-> long @-> returning gen)

  let z2_sqrt = foreign "Z2_sqrt" (gen @-> long @-> returning gen)
  let zp_div = foreign "Zp_div" (gen @-> gen @-> gen @-> long @-> returning gen)
  let zp_exp = foreign "Zp_exp" (gen @-> gen @-> pari_ulong @-> returning gen)
  let zp_inv = foreign "Zp_inv" (gen @-> gen @-> long @-> returning gen)

  let zp_invlift =
    foreign "Zp_invlift" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zp_log = foreign "Zp_log" (gen @-> gen @-> pari_ulong @-> returning gen)
  let zp_sqrt = foreign "Zp_sqrt" (gen @-> gen @-> long @-> returning gen)

  let zp_sqrtlift =
    foreign "Zp_sqrtlift" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zp_sqrtnlift =
    foreign "Zp_sqrtnlift"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpm_invlift =
    foreign "ZpM_invlift" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zpx_frobenius =
    foreign "ZpX_Frobenius" (gen @-> gen @-> long @-> returning gen)

  let zpx_zpxq_liftroot =
    foreign "ZpX_ZpXQ_liftroot"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpx_zpxq_liftroot_ea =
    foreign "ZpX_ZpXQ_liftroot_ea"
      (gen @-> gen @-> gen @-> gen @-> long @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let zpx_liftfact =
    foreign "ZpX_liftfact"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpx_liftroot =
    foreign "ZpX_liftroot" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zpx_liftroots =
    foreign "ZpX_liftroots" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zpx_roots = foreign "ZpX_roots" (gen @-> gen @-> long @-> returning gen)

  let zpxq_div =
    foreign "ZpXQ_div"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxq_inv =
    foreign "ZpXQ_inv" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxq_invlift =
    foreign "ZpXQ_invlift"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxq_log =
    foreign "ZpXQ_log" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxq_sqrt =
    foreign "ZpXQ_sqrt" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxq_sqrtnlift =
    foreign "ZpXQ_sqrtnlift"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxqm_prodfrobenius =
    foreign "ZpXQM_prodFrobenius"
      (gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxqx_digits =
    foreign "ZpXQX_digits"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxqx_divrem =
    foreign "ZpXQX_divrem"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> ptr gen
     @-> returning gen)

  let zpxqx_liftfact =
    foreign "ZpXQX_liftfact"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxqx_liftroot =
    foreign "ZpXQX_liftroot"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxqx_liftroot_vald =
    foreign "ZpXQX_liftroot_vald"
      (gen @-> gen @-> long @-> gen @-> gen @-> long @-> returning gen)

  let zpxqx_liftroots =
    foreign "ZpXQX_liftroots"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxqx_roots =
    foreign "ZpXQX_roots" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zpxqx_zpxqxq_liftroot =
    foreign "ZpXQX_ZpXQXQ_liftroot"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zq_sqrtnlift =
    foreign "Zq_sqrtnlift"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zqx_zqxq_liftroot =
    foreign "ZqX_ZqXQ_liftroot"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zqx_liftfact =
    foreign "ZqX_liftfact"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zqx_liftroot =
    foreign "ZqX_liftroot"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let zqx_roots =
    foreign "ZqX_roots" (gen @-> gen @-> gen @-> long @-> returning gen)

  let gen_zpm_dixon =
    foreign "gen_ZpM_Dixon"
      (gen @-> gen @-> gen @-> gen @-> long @-> ptr void
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> gen @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> returning gen)

  let gen_zpm_newton =
    foreign "gen_ZpM_Newton"
      (gen @-> gen @-> long @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> gen @-> gen @-> long @-> returning gen)
      @-> returning gen)

  let gen_zpx_dixon =
    foreign "gen_ZpX_Dixon"
      (gen @-> gen @-> gen @-> gen @-> long @-> ptr void
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> gen @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> returning gen)

  let gen_zpx_newton =
    foreign "gen_ZpX_Newton"
      (gen @-> gen @-> long @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> gen @-> gen @-> long @-> returning gen)
      @-> returning gen)

  let polteichmuller =
    foreign "polteichmuller" (gen @-> pari_ulong @-> long @-> returning gen)

  let polhensellift =
    foreign "polhensellift" (gen @-> gen @-> gen @-> long @-> returning gen)

  let quadratic_prec_mask =
    foreign "quadratic_prec_mask" (long @-> returning pari_ulong)

  let qx_factor = foreign "QX_factor" (gen @-> returning gen)
  let zx_factor = foreign "ZX_factor" (gen @-> returning gen)
  let zx_is_irred = foreign "ZX_is_irred" (gen @-> returning long)
  let zx_squff = foreign "ZX_squff" (gen @-> ptr gen @-> returning gen)
  let polcyclofactors = foreign "polcyclofactors" (gen @-> returning gen)
  let poliscyclo = foreign "poliscyclo" (gen @-> returning long)
  let poliscycloprod = foreign "poliscycloprod" (gen @-> returning long)
  let rg_rgc_sub = foreign "Rg_RgC_sub" (gen @-> gen @-> returning gen)
  let rgc_rg_add = foreign "RgC_Rg_add" (gen @-> gen @-> returning gen)
  let rgc_rg_div = foreign "RgC_Rg_div" (gen @-> gen @-> returning gen)
  let rgc_rg_mul = foreign "RgC_Rg_mul" (gen @-> gen @-> returning gen)
  let rgc_rg_sub = foreign "RgC_Rg_sub" (gen @-> gen @-> returning gen)
  let rgc_rgm_mul = foreign "RgC_RgM_mul" (gen @-> gen @-> returning gen)
  let rgc_rgv_mul = foreign "RgC_RgV_mul" (gen @-> gen @-> returning gen)
  let rgc_add = foreign "RgC_add" (gen @-> gen @-> returning gen)
  let rgc_is_ei = foreign "RgC_is_ei" (gen @-> returning long)
  let rgc_neg = foreign "RgC_neg" (gen @-> returning gen)
  let rgc_sub = foreign "RgC_sub" (gen @-> gen @-> returning gen)
  let rgm_rg_add = foreign "RgM_Rg_add" (gen @-> gen @-> returning gen)
end

module F11 (F : Ctypes.FOREIGN) = struct
  open F

  let rgm_rg_add_shallow =
    foreign "RgM_Rg_add_shallow" (gen @-> gen @-> returning gen)

  let rgm_rg_div = foreign "RgM_Rg_div" (gen @-> gen @-> returning gen)
  let rgm_rg_mul = foreign "RgM_Rg_mul" (gen @-> gen @-> returning gen)
  let rgm_rg_sub = foreign "RgM_Rg_sub" (gen @-> gen @-> returning gen)

  let rgm_rg_sub_shallow =
    foreign "RgM_Rg_sub_shallow" (gen @-> gen @-> returning gen)

  let rgm_rgc_mul = foreign "RgM_RgC_mul" (gen @-> gen @-> returning gen)
  let rgm_rgv_mul = foreign "RgM_RgV_mul" (gen @-> gen @-> returning gen)
  let rgm_zm_mul = foreign "RgM_ZM_mul" (gen @-> gen @-> returning gen)
  let rgm_add = foreign "RgM_add" (gen @-> gen @-> returning gen)
  let rgm_det_triangular = foreign "RgM_det_triangular" (gen @-> returning gen)
  let rgm_is_qm = foreign "RgM_is_QM" (gen @-> returning int)
  let rgm_is_zm = foreign "RgM_is_ZM" (gen @-> returning int)
  let rgm_isdiagonal = foreign "RgM_isdiagonal" (gen @-> returning int)
  let rgm_isidentity = foreign "RgM_isidentity" (gen @-> returning int)
  let rgm_isscalar = foreign "RgM_isscalar" (gen @-> gen @-> returning int)
  let rgm_mul = foreign "RgM_mul" (gen @-> gen @-> returning gen)
  let rgm_multosym = foreign "RgM_multosym" (gen @-> gen @-> returning gen)
  let rgm_neg = foreign "RgM_neg" (gen @-> returning gen)
  let rgm_powers = foreign "RgM_powers" (gen @-> long @-> returning gen)
  let rgm_sqr = foreign "RgM_sqr" (gen @-> returning gen)
  let rgm_sub = foreign "RgM_sub" (gen @-> gen @-> returning gen)
  let rgm_sumcol = foreign "RgM_sumcol" (gen @-> returning gen)
  let rgm_transmul = foreign "RgM_transmul" (gen @-> gen @-> returning gen)

  let rgm_transmultosym =
    foreign "RgM_transmultosym" (gen @-> gen @-> returning gen)

  let rgmrow_zc_mul =
    foreign "RgMrow_zc_mul" (gen @-> gen @-> long @-> returning gen)

  let rgm_zc_mul = foreign "RgM_zc_mul" (gen @-> gen @-> returning gen)
  let rgm_zm_mul = foreign "RgM_zm_mul" (gen @-> gen @-> returning gen)

  let rgmrow_rgc_mul =
    foreign "RgMrow_RgC_mul" (gen @-> gen @-> long @-> returning gen)

  let rgv_rgm_mul = foreign "RgV_RgM_mul" (gen @-> gen @-> returning gen)
  let rgv_rgc_mul = foreign "RgV_RgC_mul" (gen @-> gen @-> returning gen)
  let rgv_rg_mul = foreign "RgV_Rg_mul" (gen @-> gen @-> returning gen)
  let rgv_add = foreign "RgV_add" (gen @-> gen @-> returning gen)
  let rgv_dotproduct = foreign "RgV_dotproduct" (gen @-> gen @-> returning gen)
  let rgv_dotsquare = foreign "RgV_dotsquare" (gen @-> returning gen)
  let rgv_is_zmv = foreign "RgV_is_ZMV" (gen @-> returning int)
  let rgv_kill0 = foreign "RgV_kill0" (gen @-> returning gen)
  let rgv_neg = foreign "RgV_neg" (gen @-> returning gen)
  let rgv_prod = foreign "RgV_prod" (gen @-> returning gen)
  let rgv_sub = foreign "RgV_sub" (gen @-> gen @-> returning gen)
  let rgv_sum = foreign "RgV_sum" (gen @-> returning gen)
  let rgv_sumpart = foreign "RgV_sumpart" (gen @-> long @-> returning gen)

  let rgv_sumpart2 =
    foreign "RgV_sumpart2" (gen @-> long @-> long @-> returning gen)

  let rgv_zc_mul = foreign "RgV_zc_mul" (gen @-> gen @-> returning gen)
  let rgv_zm_mul = foreign "RgV_zm_mul" (gen @-> gen @-> returning gen)
  let rgx_rgm_eval = foreign "RgX_RgM_eval" (gen @-> gen @-> returning gen)
  let rgx_rgmv_eval = foreign "RgX_RgMV_eval" (gen @-> gen @-> returning gen)
  let isdiagonal = foreign "isdiagonal" (gen @-> returning int)
  let matid = foreign "matid" (long @-> returning gen)
  let scalarcol = foreign "scalarcol" (gen @-> long @-> returning gen)

  let scalarcol_shallow =
    foreign "scalarcol_shallow" (gen @-> long @-> returning gen)

  let scalarmat = foreign "scalarmat" (gen @-> long @-> returning gen)

  let scalarmat_shallow =
    foreign "scalarmat_shallow" (gen @-> long @-> returning gen)

  let scalarmat_s = foreign "scalarmat_s" (long @-> long @-> returning gen)

  let kronecker_to_mod =
    foreign "Kronecker_to_mod" (gen @-> gen @-> returning gen)

  let qx_zxqv_eval =
    foreign "QX_ZXQV_eval" (gen @-> gen @-> gen @-> returning gen)

  let qxq_charpoly =
    foreign "QXQ_charpoly" (gen @-> gen @-> long @-> returning gen)

  let qxq_powers = foreign "QXQ_powers" (gen @-> long @-> gen @-> returning gen)

  let qxq_to_mod_shallow =
    foreign "QXQ_to_mod_shallow" (gen @-> gen @-> returning gen)

  let qxqc_to_mod_shallow =
    foreign "QXQC_to_mod_shallow" (gen @-> gen @-> returning gen)

  let qxqm_to_mod_shallow =
    foreign "QXQM_to_mod_shallow" (gen @-> gen @-> returning gen)

  let qxqv_to_mod = foreign "QXQV_to_mod" (gen @-> gen @-> returning gen)

  let qxqx_homogenous_evalpow =
    foreign "QXQX_homogenous_evalpow"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let qxqx_to_mod_shallow =
    foreign "QXQX_to_mod_shallow" (gen @-> gen @-> returning gen)

  let qxqxv_to_mod = foreign "QXQXV_to_mod" (gen @-> gen @-> returning gen)

  let qxv_qxq_eval =
    foreign "QXV_QXQ_eval" (gen @-> gen @-> gen @-> returning gen)

  let qxy_qxq_evalx =
    foreign "QXY_QXQ_evalx" (gen @-> gen @-> gen @-> returning gen)

  let rg_rgx_sub = foreign "Rg_RgX_sub" (gen @-> gen @-> returning gen)
  let rg_get_0 = foreign "Rg_get_0" (gen @-> returning gen)
  let rg_get_1 = foreign "Rg_get_1" (gen @-> returning gen)
  let rg_to_rgc = foreign "Rg_to_RgC" (gen @-> long @-> returning gen)
  let rgm_to_rgxv = foreign "RgM_to_RgXV" (gen @-> long @-> returning gen)

  let rgm_to_rgxv_reverse =
    foreign "RgM_to_RgXV_reverse" (gen @-> long @-> returning gen)

  let rgm_to_rgxx =
    foreign "RgM_to_RgXX" (gen @-> long @-> long @-> returning gen)

  let rgv_to_rgx = foreign "RgV_to_RgX" (gen @-> long @-> returning gen)
  let rgv_to_rgm = foreign "RgV_to_RgM" (gen @-> long @-> returning gen)

  let rgv_to_rgx_reverse =
    foreign "RgV_to_RgX_reverse" (gen @-> long @-> returning gen)

  let rgx_rgxq_eval =
    foreign "RgX_RgXQ_eval" (gen @-> gen @-> gen @-> returning gen)

  let rgx_rgxqv_eval =
    foreign "RgX_RgXQV_eval" (gen @-> gen @-> gen @-> returning gen)

  let rgx_rgxn_eval =
    foreign "RgX_RgXn_eval" (gen @-> gen @-> long @-> returning gen)

  let rgx_rgxnv_eval =
    foreign "RgX_RgXnV_eval" (gen @-> gen @-> long @-> returning gen)

  let rgx_rg_add = foreign "RgX_Rg_add" (gen @-> gen @-> returning gen)

  let rgx_rg_add_shallow =
    foreign "RgX_Rg_add_shallow" (gen @-> gen @-> returning gen)

  let rgx_rg_div = foreign "RgX_Rg_div" (gen @-> gen @-> returning gen)
  let rgx_rg_divexact = foreign "RgX_Rg_divexact" (gen @-> gen @-> returning gen)
  let rgx_rg_eval_bk = foreign "RgX_Rg_eval_bk" (gen @-> gen @-> returning gen)
  let rgx_rg_mul = foreign "RgX_Rg_mul" (gen @-> gen @-> returning gen)
  let rgx_rg_sub = foreign "RgX_Rg_sub" (gen @-> gen @-> returning gen)
  let rgx_rgv_eval = foreign "RgX_RgV_eval" (gen @-> gen @-> returning gen)
  let rgx_add = foreign "RgX_add" (gen @-> gen @-> returning gen)

  let rgx_addmulxn_shallow =
    foreign "RgX_addmulXn_shallow" (gen @-> gen @-> long @-> returning gen)

  let rgx_addmulxn =
    foreign "RgX_addmulXn" (gen @-> gen @-> long @-> returning gen)

  let rgx_addspec =
    foreign "RgX_addspec" (gen @-> gen @-> long @-> long @-> returning gen)

  let rgx_addspec_shallow =
    foreign "RgX_addspec_shallow"
      (gen @-> gen @-> long @-> long @-> returning gen)

  let rgx_affine = foreign "RgX_affine" (gen @-> gen @-> gen @-> returning gen)
  let rgx_blocks = foreign "RgX_blocks" (gen @-> long @-> long @-> returning gen)
  let rgx_deflate = foreign "RgX_deflate" (gen @-> long @-> returning gen)
  let rgx_deriv = foreign "RgX_deriv" (gen @-> returning gen)
  let rgx_digits = foreign "RgX_digits" (gen @-> gen @-> returning gen)

  let rgx_div_by_x_x =
    foreign "RgX_div_by_X_x" (gen @-> gen @-> ptr gen @-> returning gen)

  let rgx_divrem =
    foreign "RgX_divrem" (gen @-> gen @-> ptr gen @-> returning gen)
end

module F12 (F : Ctypes.FOREIGN) = struct
  open F

  let rgx_divs = foreign "RgX_divs" (gen @-> long @-> returning gen)
  let rgx_equal = foreign "RgX_equal" (gen @-> gen @-> returning long)

  let rgx_even_odd =
    foreign "RgX_even_odd" (gen @-> ptr gen @-> ptr gen @-> returning void)

  let rgx_homogenize = foreign "RgX_homogenize" (gen @-> long @-> returning gen)

  let rgx_homogenous_evalpow =
    foreign "RgX_homogenous_evalpow" (gen @-> gen @-> gen @-> returning gen)

  let rgx_inflate = foreign "RgX_inflate" (gen @-> long @-> returning gen)
  let rgx_mul = foreign "RgX_mul" (gen @-> gen @-> returning gen)
  let rgx_mul_i = foreign "RgX_mul_i" (gen @-> gen @-> returning gen)

  let rgx_mul_normalized =
    foreign "RgX_mul_normalized"
      (gen @-> long @-> gen @-> long @-> returning gen)

  let rgx_mul2n = foreign "RgX_mul2n" (gen @-> long @-> returning gen)
  let rgx_mulxn = foreign "RgX_mulXn" (gen @-> long @-> returning gen)

  let rgx_mulhigh_i =
    foreign "RgX_mulhigh_i" (gen @-> gen @-> long @-> returning gen)

  let rgx_muls = foreign "RgX_muls" (gen @-> long @-> returning gen)

  let rgx_mulspec =
    foreign "RgX_mulspec" (gen @-> gen @-> long @-> long @-> returning gen)

  let rgx_neg = foreign "RgX_neg" (gen @-> returning gen)
  let rgx_normalize = foreign "RgX_normalize" (gen @-> returning gen)

  let rgx_pseudodivrem =
    foreign "RgX_pseudodivrem" (gen @-> gen @-> ptr gen @-> returning gen)

  let rgx_pseudorem = foreign "RgX_pseudorem" (gen @-> gen @-> returning gen)
  let rgx_recip = foreign "RgX_recip" (gen @-> returning gen)
  let rgx_recip_i = foreign "RgX_recip_i" (gen @-> returning gen)
  let rgx_recip_shallow = foreign "RgX_recip_shallow" (gen @-> returning gen)
  let rgx_rem = foreign "RgX_rem" (gen @-> gen @-> returning gen)

  let rgx_renormalize_lg =
    foreign "RgX_renormalize_lg" (gen @-> long @-> returning gen)

  let rgx_rescale = foreign "RgX_rescale" (gen @-> gen @-> returning gen)

  let rgx_rotate_shallow =
    foreign "RgX_rotate_shallow" (gen @-> long @-> long @-> returning gen)

  let rgx_shift = foreign "RgX_shift" (gen @-> long @-> returning gen)

  let rgx_shift_shallow =
    foreign "RgX_shift_shallow" (gen @-> long @-> returning gen)

  let rgx_splitting = foreign "RgX_splitting" (gen @-> long @-> returning gen)
  let rgx_sqr = foreign "RgX_sqr" (gen @-> returning gen)
  let rgx_sqr_i = foreign "RgX_sqr_i" (gen @-> returning gen)
  let rgx_sqrhigh_i = foreign "RgX_sqrhigh_i" (gen @-> long @-> returning gen)
  let rgx_sqrspec = foreign "RgX_sqrspec" (gen @-> long @-> returning gen)
  let rgx_sub = foreign "RgX_sub" (gen @-> gen @-> returning gen)
  let rgx_to_rgc = foreign "RgX_to_RgC" (gen @-> long @-> returning gen)
  let rgx_translate = foreign "RgX_translate" (gen @-> gen @-> returning gen)
  let rgx_unscale = foreign "RgX_unscale" (gen @-> gen @-> returning gen)

  let rgxq_matrix_pow =
    foreign "RgXQ_matrix_pow" (gen @-> long @-> long @-> gen @-> returning gen)

  let rgxq_norm = foreign "RgXQ_norm" (gen @-> gen @-> returning gen)
  let rgxq_pow = foreign "RgXQ_pow" (gen @-> gen @-> gen @-> returning gen)

  let rgxq_powers =
    foreign "RgXQ_powers" (gen @-> long @-> gen @-> returning gen)

  let rgxq_powu =
    foreign "RgXQ_powu" (gen @-> pari_ulong @-> gen @-> returning gen)

  let rgxq_trace = foreign "RgXQ_trace" (gen @-> gen @-> returning gen)
  let rgxqc_red = foreign "RgXQC_red" (gen @-> gen @-> returning gen)
  let rgxqm_mul = foreign "RgXQM_mul" (gen @-> gen @-> gen @-> returning gen)
  let rgxqm_red = foreign "RgXQM_red" (gen @-> gen @-> returning gen)

  let rgxqv_rgxq_mul =
    foreign "RgXQV_RgXQ_mul" (gen @-> gen @-> gen @-> returning gen)

  let rgxqv_factorback =
    foreign "RgXQV_factorback" (gen @-> gen @-> gen @-> returning gen)

  let rgxqv_red = foreign "RgXQV_red" (gen @-> gen @-> returning gen)

  let rgxqx_rgxq_mul =
    foreign "RgXQX_RgXQ_mul" (gen @-> gen @-> gen @-> returning gen)

  let rgxqx_divrem =
    foreign "RgXQX_divrem" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let rgxqx_mul = foreign "RgXQX_mul" (gen @-> gen @-> gen @-> returning gen)

  let rgxqx_powers =
    foreign "RgXQX_powers" (gen @-> long @-> gen @-> returning gen)

  let rgxqx_pseudodivrem =
    foreign "RgXQX_pseudodivrem"
      (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let rgxqx_pseudorem =
    foreign "RgXQX_pseudorem" (gen @-> gen @-> gen @-> returning gen)

  let rgxqx_red = foreign "RgXQX_red" (gen @-> gen @-> returning gen)
  let rgxqx_sqr = foreign "RgXQX_sqr" (gen @-> gen @-> returning gen)

  let rgxqx_translate =
    foreign "RgXQX_translate" (gen @-> gen @-> gen @-> returning gen)

  let rgxv_rgv_eval = foreign "RgXV_RgV_eval" (gen @-> gen @-> returning gen)
  let rgxv_prod = foreign "RgXV_prod" (gen @-> returning gen)
  let rgxv_rescale = foreign "RgXV_rescale" (gen @-> gen @-> returning gen)
  let rgxv_to_rgm = foreign "RgXV_to_RgM" (gen @-> long @-> returning gen)
  let rgxv_unscale = foreign "RgXV_unscale" (gen @-> gen @-> returning gen)
  let rgxx_to_rgm = foreign "RgXX_to_RgM" (gen @-> long @-> returning gen)
  let rgxy_degreex = foreign "RgXY_degreex" (gen @-> returning long)
  let rgxy_derivx = foreign "RgXY_derivx" (gen @-> returning gen)
  let rgxy_swap = foreign "RgXY_swap" (gen @-> long @-> long @-> returning gen)

  let rgxy_swapspec =
    foreign "RgXY_swapspec" (gen @-> long @-> long @-> long @-> returning gen)

  let rgxn_div = foreign "RgXn_div" (gen @-> gen @-> long @-> returning gen)
  let rgxn_div_i = foreign "RgXn_div_i" (gen @-> gen @-> long @-> returning gen)
  let rgxn_eval = foreign "RgXn_eval" (gen @-> gen @-> long @-> returning gen)
  let rgxn_exp = foreign "RgXn_exp" (gen @-> long @-> returning gen)
  let rgxn_expint = foreign "RgXn_expint" (gen @-> long @-> returning gen)
  let rgxn_inv = foreign "RgXn_inv" (gen @-> long @-> returning gen)
  let rgxn_inv_i = foreign "RgXn_inv_i" (gen @-> long @-> returning gen)
  let rgxn_mul = foreign "RgXn_mul" (gen @-> gen @-> long @-> returning gen)

  let rgxn_powers =
    foreign "RgXn_powers" (gen @-> long @-> long @-> returning gen)

  let rgxn_recip_shallow =
    foreign "RgXn_recip_shallow" (gen @-> long @-> returning gen)

  let rgxn_red_shallow =
    foreign "RgXn_red_shallow" (gen @-> long @-> returning gen)

  let rgxn_reverse = foreign "RgXn_reverse" (gen @-> long @-> returning gen)
  let rgxn_sqr = foreign "RgXn_sqr" (gen @-> long @-> returning gen)
  let rgxn_sqrt = foreign "RgXn_sqrt" (gen @-> long @-> returning gen)

  let rgxnv_red_shallow =
    foreign "RgXnV_red_shallow" (gen @-> long @-> returning gen)

  let rgxn_powu =
    foreign "RgXn_powu" (gen @-> pari_ulong @-> long @-> returning gen)

  let rgxn_powu_i =
    foreign "RgXn_powu_i" (gen @-> pari_ulong @-> long @-> returning gen)

  let zx_translate = foreign "ZX_translate" (gen @-> gen @-> returning gen)
  let zx_unscale2n = foreign "ZX_unscale2n" (gen @-> long @-> returning gen)
  let zx_unscale = foreign "ZX_unscale" (gen @-> gen @-> returning gen)
  let zx_unscale_div = foreign "ZX_unscale_div" (gen @-> gen @-> returning gen)

  let zx_unscale_divpow =
    foreign "ZX_unscale_divpow" (gen @-> gen @-> long @-> returning gen)

  let zx_z_unscale = foreign "ZX_z_unscale" (gen @-> long @-> returning gen)
  let zxq_powers = foreign "ZXQ_powers" (gen @-> long @-> gen @-> returning gen)

  let zxq_powu =
    foreign "ZXQ_powu" (gen @-> pari_ulong @-> gen @-> returning gen)

  let zxqx_dvd = foreign "ZXQX_dvd" (gen @-> gen @-> gen @-> returning int)

  let brent_kung_optpow =
    foreign "brent_kung_optpow" (long @-> long @-> long @-> returning long)

  let gen_bkeval =
    foreign "gen_bkeval"
      (gen @-> long @-> gen @-> int @-> ptr void @-> ptr bb_algebra
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> long @-> gen @-> returning gen)
      @-> returning gen)

  let gen_bkeval_powers =
    foreign "gen_bkeval_powers"
      (gen @-> long @-> gen @-> ptr void @-> ptr bb_algebra
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> long @-> gen @-> returning gen)
      @-> returning gen)

  let get_rg_algebra =
    foreign "get_Rg_algebra" (void @-> returning (ptr bb_algebra))

  let rfrac_deflate_order =
    foreign "rfrac_deflate_order" (gen @-> returning long)

  let rfrac_deflate_max =
    foreign "rfrac_deflate_max" (gen @-> ptr long @-> returning gen)

  let rfrac_deflate = foreign "rfrac_deflate" (gen @-> long @-> returning gen)
end

module F13 (F : Ctypes.FOREIGN) = struct
  open F

  let zgc_g_mul_inplace =
    foreign "ZGC_G_mul_inplace" (gen @-> gen @-> returning void)

  let zgcs_add = foreign "ZGCs_add" (gen @-> gen @-> returning gen)
  let g_zgc_mul = foreign "G_ZGC_mul" (gen @-> gen @-> returning gen)
  let g_zg_mul = foreign "G_ZG_mul" (gen @-> gen @-> returning gen)
  let zgc_g_mul = foreign "ZGC_G_mul" (gen @-> gen @-> returning gen)
  let zgc_z_mul = foreign "ZGC_Z_mul" (gen @-> gen @-> returning gen)
  let zg_g_mul = foreign "ZG_G_mul" (gen @-> gen @-> returning gen)
  let zg_z_mul = foreign "ZG_Z_mul" (gen @-> gen @-> returning gen)
  let zg_add = foreign "ZG_add" (gen @-> gen @-> returning gen)
  let zg_mul = foreign "ZG_mul" (gen @-> gen @-> returning gen)
  let zg_neg = foreign "ZG_neg" (gen @-> returning gen)
  let zg_normalize = foreign "ZG_normalize" (gen @-> returning gen)
  let zg_sub = foreign "ZG_sub" (gen @-> gen @-> returning gen)

  let flc_lincomb1_inplace =
    foreign "Flc_lincomb1_inplace"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning void)

  let vecsmall_prod = foreign "vecsmall_prod" (gen @-> returning gen)
  let qm_qc_mul = foreign "QM_QC_mul" (gen @-> gen @-> returning gen)
  let qm_det = foreign "QM_det" (gen @-> returning gen)
  let qm_ker = foreign "QM_ker" (gen @-> returning gen)
  let qm_mul = foreign "QM_mul" (gen @-> gen @-> returning gen)
  let qm_sqr = foreign "QM_sqr" (gen @-> returning gen)
  let rgm_check_zm = foreign "RgM_check_ZM" (gen @-> string @-> returning void)
  let rgv_check_zv = foreign "RgV_check_ZV" (gen @-> string @-> returning void)
  let z_zc_sub = foreign "Z_ZC_sub" (gen @-> gen @-> returning gen)
  let zv_zc_mul = foreign "ZV_zc_mul" (gen @-> gen @-> returning gen)
  let zc_q_mul = foreign "ZC_Q_mul" (gen @-> gen @-> returning gen)
  let zc_z_add = foreign "ZC_Z_add" (gen @-> gen @-> returning gen)
  let zc_z_div = foreign "ZC_Z_div" (gen @-> gen @-> returning gen)
  let zc_z_divexact = foreign "ZC_Z_divexact" (gen @-> gen @-> returning gen)
  let zc_z_mul = foreign "ZC_Z_mul" (gen @-> gen @-> returning gen)
  let zc_z_sub = foreign "ZC_Z_sub" (gen @-> gen @-> returning gen)
  let zc_zv_mul = foreign "ZC_ZV_mul" (gen @-> gen @-> returning gen)

  let zc_divexactu =
    foreign "ZC_divexactu" (gen @-> pari_ulong @-> returning gen)

  let zc_add = foreign "ZC_add" (gen @-> gen @-> returning gen)
  let zc_copy = foreign "ZC_copy" (gen @-> returning gen)

  let zc_hnfremdiv =
    foreign "ZC_hnfremdiv" (gen @-> gen @-> ptr gen @-> returning gen)

  let zc_is_ei = foreign "ZC_is_ei" (gen @-> returning long)

  let zc_lincomb =
    foreign "ZC_lincomb" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let zc_lincomb1_inplace =
    foreign "ZC_lincomb1_inplace" (gen @-> gen @-> gen @-> returning void)

  let zc_lincomb1_inplace_i =
    foreign "ZC_lincomb1_inplace_i"
      (gen @-> gen @-> gen @-> long @-> returning void)

  let zc_neg = foreign "ZC_neg" (gen @-> returning gen)
  let zc_reducemodlll = foreign "ZC_reducemodlll" (gen @-> gen @-> returning gen)

  let zc_reducemodmatrix =
    foreign "ZC_reducemodmatrix" (gen @-> gen @-> returning gen)

  let zc_sub = foreign "ZC_sub" (gen @-> gen @-> returning gen)
  let zc_z_mul = foreign "ZC_z_mul" (gen @-> long @-> returning gen)
  let zm_q_mul = foreign "ZM_Q_mul" (gen @-> gen @-> returning gen)
  let zm_zc_mul = foreign "ZM_ZC_mul" (gen @-> gen @-> returning gen)
  let zm_z_div = foreign "ZM_Z_div" (gen @-> gen @-> returning gen)
  let zm_z_divexact = foreign "ZM_Z_divexact" (gen @-> gen @-> returning gen)
  let zm_z_mul = foreign "ZM_Z_mul" (gen @-> gen @-> returning gen)
  let zm_add = foreign "ZM_add" (gen @-> gen @-> returning gen)
  let zm_copy = foreign "ZM_copy" (gen @-> returning gen)
  let zm_det_triangular = foreign "ZM_det_triangular" (gen @-> returning gen)
  let zm_diag_mul = foreign "ZM_diag_mul" (gen @-> gen @-> returning gen)

  let zm_divexactu =
    foreign "ZM_divexactu" (gen @-> pari_ulong @-> returning gen)

  let zm_equal = foreign "ZM_equal" (gen @-> gen @-> returning int)
  let zm_equal0 = foreign "ZM_equal0" (gen @-> returning int)

  let zm_hnfdivrem =
    foreign "ZM_hnfdivrem" (gen @-> gen @-> ptr gen @-> returning gen)

  let zm_ishnf = foreign "ZM_ishnf" (gen @-> returning int)
  let zm_isdiagonal = foreign "ZM_isdiagonal" (gen @-> returning int)
  let zm_isidentity = foreign "ZM_isidentity" (gen @-> returning int)
  let zm_isscalar = foreign "ZM_isscalar" (gen @-> gen @-> returning int)
  let zm_max_lg = foreign "ZM_max_lg" (gen @-> returning long)
  let zm_mul = foreign "ZM_mul" (gen @-> gen @-> returning gen)
  let zm_mul_diag = foreign "ZM_mul_diag" (gen @-> gen @-> returning gen)
  let zm_multosym = foreign "ZM_multosym" (gen @-> gen @-> returning gen)
  let zm_neg = foreign "ZM_neg" (gen @-> returning gen)
  let zm_nm_mul = foreign "ZM_nm_mul" (gen @-> gen @-> returning gen)
  let zm_pow = foreign "ZM_pow" (gen @-> gen @-> returning gen)
  let zm_powu = foreign "ZM_powu" (gen @-> pari_ulong @-> returning gen)
  let zm_reducemodlll = foreign "ZM_reducemodlll" (gen @-> gen @-> returning gen)

  let zm_reducemodmatrix =
    foreign "ZM_reducemodmatrix" (gen @-> gen @-> returning gen)

  let zm_sqr = foreign "ZM_sqr" (gen @-> returning gen)
  let zm_sub = foreign "ZM_sub" (gen @-> gen @-> returning gen)
  let zm_supnorm = foreign "ZM_supnorm" (gen @-> returning gen)
  let zm_to_flm = foreign "ZM_to_Flm" (gen @-> pari_ulong @-> returning gen)
  let zm_to_zm = foreign "ZM_to_zm" (gen @-> returning gen)
  let zm_transmul = foreign "ZM_transmul" (gen @-> gen @-> returning gen)

  let zm_transmultosym =
    foreign "ZM_transmultosym" (gen @-> gen @-> returning gen)

  let zmv_to_zmv = foreign "ZMV_to_zmV" (gen @-> returning gen)
  let zm_togglesign = foreign "ZM_togglesign" (gen @-> returning void)
  let zm_zc_mul = foreign "ZM_zc_mul" (gen @-> gen @-> returning gen)
  let zm_zm_mul = foreign "ZM_zm_mul" (gen @-> gen @-> returning gen)

  let zmrow_zc_mul =
    foreign "ZMrow_ZC_mul" (gen @-> gen @-> long @-> returning gen)

  let zmrow_equal0 = foreign "ZMrow_equal0" (gen @-> long @-> returning int)
  let zv_zm_mul = foreign "ZV_ZM_mul" (gen @-> gen @-> returning gen)
  let zv_abscmp = foreign "ZV_abscmp" (gen @-> gen @-> returning int)
  let zv_cmp = foreign "ZV_cmp" (gen @-> gen @-> returning int)
  let zv_content = foreign "ZV_content" (gen @-> returning gen)
  let zv_dotproduct = foreign "ZV_dotproduct" (gen @-> gen @-> returning gen)
  let zv_dotsquare = foreign "ZV_dotsquare" (gen @-> returning gen)
  let zv_equal = foreign "ZV_equal" (gen @-> gen @-> returning int)
  let zv_equal0 = foreign "ZV_equal0" (gen @-> returning int)
  let zv_max_lg = foreign "ZV_max_lg" (gen @-> returning long)
  let zv_neg_inplace = foreign "ZV_neg_inplace" (gen @-> returning void)
  let zv_prod = foreign "ZV_prod" (gen @-> returning gen)
  let zv_sum = foreign "ZV_sum" (gen @-> returning gen)
  let zv_to_flv = foreign "ZV_to_Flv" (gen @-> pari_ulong @-> returning gen)
  let zv_to_nv = foreign "ZV_to_nv" (gen @-> returning gen)
  let zv_togglesign = foreign "ZV_togglesign" (gen @-> returning void)
  let gram_matrix = foreign "gram_matrix" (gen @-> returning gen)
end

module F14 (F : Ctypes.FOREIGN) = struct
  open F

  let nm_z_mul = foreign "nm_Z_mul" (gen @-> gen @-> returning gen)
  let zm_mul = foreign "zm_mul" (gen @-> gen @-> returning gen)
  let zm_to_flm = foreign "zm_to_Flm" (gen @-> pari_ulong @-> returning gen)
  let zm_to_zm = foreign "zm_to_ZM" (gen @-> returning gen)
  let zm_zc_mul = foreign "zm_zc_mul" (gen @-> gen @-> returning gen)
  let zmv_to_zmv = foreign "zmV_to_ZMV" (gen @-> returning gen)
  let zv_abs = foreign "zv_abs" (gen @-> returning gen)
  let zv_content = foreign "zv_content" (gen @-> returning long)
  let zv_dotproduct = foreign "zv_dotproduct" (gen @-> gen @-> returning long)
  let zv_equal = foreign "zv_equal" (gen @-> gen @-> returning int)
  let zv_equal0 = foreign "zv_equal0" (gen @-> returning int)
  let zv_neg = foreign "zv_neg" (gen @-> returning gen)
  let zv_neg_inplace = foreign "zv_neg_inplace" (gen @-> returning gen)
  let zv_prod = foreign "zv_prod" (gen @-> returning long)
  let zv_prod_z = foreign "zv_prod_Z" (gen @-> returning gen)
  let zv_sum = foreign "zv_sum" (gen @-> returning long)
  let zv_sumpart = foreign "zv_sumpart" (gen @-> long @-> returning long)
  let zv_to_flv = foreign "zv_to_Flv" (gen @-> pari_ulong @-> returning gen)
  let zv_z_mul = foreign "zv_z_mul" (gen @-> long @-> returning gen)
  let zv_zm_mul = foreign "zv_ZM_mul" (gen @-> gen @-> returning gen)
  let zvv_equal = foreign "zvV_equal" (gen @-> gen @-> returning int)

  let kronecker_to_zxqx =
    foreign "Kronecker_to_ZXQX" (gen @-> gen @-> returning gen)

  let kronecker_to_zxx =
    foreign "Kronecker_to_ZXX" (gen @-> long @-> long @-> returning gen)

  let qx_zx_rem = foreign "QX_ZX_rem" (gen @-> gen @-> returning gen)
  let qx_mul = foreign "QX_mul" (gen @-> gen @-> returning gen)
  let qx_sqr = foreign "QX_sqr" (gen @-> returning gen)
  let qxqm_mul = foreign "QXQM_mul" (gen @-> gen @-> gen @-> returning gen)
  let qxqm_sqr = foreign "QXQM_sqr" (gen @-> gen @-> returning gen)

  let qxqx_qxq_mul =
    foreign "QXQX_QXQ_mul" (gen @-> gen @-> gen @-> returning gen)

  let qxqx_mul = foreign "QXQX_mul" (gen @-> gen @-> gen @-> returning gen)

  let qxqx_powers =
    foreign "QXQX_powers" (gen @-> long @-> gen @-> returning gen)

  let qxqx_sqr = foreign "QXQX_sqr" (gen @-> gen @-> returning gen)
  let rgx_check_qx = foreign "RgX_check_QX" (gen @-> string @-> returning void)
  let rgx_check_zx = foreign "RgX_check_ZX" (gen @-> string @-> returning void)
  let rgx_check_zxx = foreign "RgX_check_ZXX" (gen @-> string @-> returning void)
  let z_zx_sub = foreign "Z_ZX_sub" (gen @-> gen @-> returning gen)
  let zx_z_add = foreign "ZX_Z_add" (gen @-> gen @-> returning gen)

  let zx_z_add_shallow =
    foreign "ZX_Z_add_shallow" (gen @-> gen @-> returning gen)

  let zx_z_divexact = foreign "ZX_Z_divexact" (gen @-> gen @-> returning gen)
  let zx_z_eval = foreign "ZX_Z_eval" (gen @-> gen @-> returning gen)
  let zx_z_mul = foreign "ZX_Z_mul" (gen @-> gen @-> returning gen)
  let zx_z_sub = foreign "ZX_Z_sub" (gen @-> gen @-> returning gen)
  let zx_add = foreign "ZX_add" (gen @-> gen @-> returning gen)
  let zx_affine = foreign "ZX_affine" (gen @-> gen @-> gen @-> returning gen)
  let zx_copy = foreign "ZX_copy" (gen @-> returning gen)
  let zx_deriv = foreign "ZX_deriv" (gen @-> returning gen)
  let zx_digits = foreign "ZX_digits" (gen @-> gen @-> returning gen)

  let zxv_zx_fromdigits =
    foreign "ZXV_ZX_fromdigits" (gen @-> gen @-> returning gen)

  let zx_div_by_x_1 = foreign "ZX_div_by_X_1" (gen @-> ptr gen @-> returning gen)

  let zx_divuexact =
    foreign "ZX_divuexact" (gen @-> pari_ulong @-> returning gen)

  let zx_equal = foreign "ZX_equal" (gen @-> gen @-> returning int)
  let zx_eval1 = foreign "ZX_eval1" (gen @-> returning gen)
  let zx_max_lg = foreign "ZX_max_lg" (gen @-> returning long)
  let zx_mod_xnm1 = foreign "ZX_mod_Xnm1" (gen @-> pari_ulong @-> returning gen)
  let zx_mul = foreign "ZX_mul" (gen @-> gen @-> returning gen)

  let zx_mulspec =
    foreign "ZX_mulspec" (gen @-> gen @-> long @-> long @-> returning gen)

  let zx_mulu = foreign "ZX_mulu" (gen @-> pari_ulong @-> returning gen)
  let zx_neg = foreign "ZX_neg" (gen @-> returning gen)
  let zx_rem = foreign "ZX_rem" (gen @-> gen @-> returning gen)
  let zx_remi2n = foreign "ZX_remi2n" (gen @-> long @-> returning gen)
  let zx_rescale2n = foreign "ZX_rescale2n" (gen @-> long @-> returning gen)
  let zx_rescale = foreign "ZX_rescale" (gen @-> gen @-> returning gen)
  let zx_rescale_lt = foreign "ZX_rescale_lt" (gen @-> returning gen)
  let zx_shifti = foreign "ZX_shifti" (gen @-> long @-> returning gen)
  let zx_sqr = foreign "ZX_sqr" (gen @-> returning gen)
  let zx_sqrspec = foreign "ZX_sqrspec" (gen @-> long @-> returning gen)
  let zx_sub = foreign "ZX_sub" (gen @-> gen @-> returning gen)
  let zx_val = foreign "ZX_val" (gen @-> returning long)
  let zx_valrem = foreign "ZX_valrem" (gen @-> ptr gen @-> returning long)

  let zxc_to_flxc =
    foreign "ZXC_to_FlxC" (gen @-> pari_ulong @-> long @-> returning gen)

  let zxm_to_flxm =
    foreign "ZXM_to_FlxM" (gen @-> pari_ulong @-> long @-> returning gen)

  let zxqm_mul = foreign "ZXQM_mul" (gen @-> gen @-> gen @-> returning gen)
  let zxqm_sqr = foreign "ZXQM_sqr" (gen @-> gen @-> returning gen)

  let zxqx_zxq_mul =
    foreign "ZXQX_ZXQ_mul" (gen @-> gen @-> gen @-> returning gen)

  let zxqx_sqr = foreign "ZXQX_sqr" (gen @-> gen @-> returning gen)
  let zxqx_mul = foreign "ZXQX_mul" (gen @-> gen @-> gen @-> returning gen)
  let zxt_remi2n = foreign "ZXT_remi2n" (gen @-> long @-> returning gen)
  let zxv_z_mul = foreign "ZXV_Z_mul" (gen @-> gen @-> returning gen)
  let zxv_dotproduct = foreign "ZXV_dotproduct" (gen @-> gen @-> returning gen)
  let zxv_equal = foreign "ZXV_equal" (gen @-> gen @-> returning int)
  let zxv_remi2n = foreign "ZXV_remi2n" (gen @-> long @-> returning gen)
  let zxx_z_divexact = foreign "ZXX_Z_divexact" (gen @-> gen @-> returning gen)
  let zxx_z_mul = foreign "ZXX_Z_mul" (gen @-> gen @-> returning gen)

  let zxx_z_add_shallow =
    foreign "ZXX_Z_add_shallow" (gen @-> gen @-> returning gen)

  let zxx_evalx0 = foreign "ZXX_evalx0" (gen @-> returning gen)
  let zxx_max_lg = foreign "ZXX_max_lg" (gen @-> returning long)

  let zxx_mul_kronecker =
    foreign "ZXX_mul_Kronecker" (gen @-> gen @-> long @-> returning gen)

  let zxx_renormalize =
    foreign "ZXX_renormalize" (gen @-> long @-> returning gen)

  let zxx_sqr_kronecker =
    foreign "ZXX_sqr_Kronecker" (gen @-> long @-> returning gen)

  let rgxx_to_kronecker =
    foreign "RgXX_to_Kronecker" (gen @-> long @-> returning gen)

  let rgxx_to_kronecker_spec =
    foreign "RgXX_to_Kronecker_spec" (gen @-> long @-> long @-> returning gen)

  let zxn_mul = foreign "ZXn_mul" (gen @-> gen @-> long @-> returning gen)
  let zxn_sqr = foreign "ZXn_sqr" (gen @-> long @-> returning gen)
  let scalar_zx = foreign "scalar_ZX" (gen @-> long @-> returning gen)

  let scalar_zx_shallow =
    foreign "scalar_ZX_shallow" (gen @-> long @-> returning gen)

  let zx_to_zx = foreign "zx_to_ZX" (gen @-> returning gen)
  let zx_z_divexact = foreign "zx_z_divexact" (gen @-> long @-> returning gen)

  let alg_centralproj =
    foreign "alg_centralproj" (gen @-> gen @-> long @-> returning gen)

  let alg_complete =
    foreign "alg_complete"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let alg_csa_table =
    foreign "alg_csa_table" (gen @-> gen @-> long @-> long @-> returning gen)
end

module F15 (F : Ctypes.FOREIGN) = struct
  open F

  let alg_cyclic =
    foreign "alg_cyclic" (gen @-> gen @-> gen @-> long @-> returning gen)

  let alg_get_absdim = foreign "alg_get_absdim" (gen @-> returning long)

  let alg_get_abssplitting =
    foreign "alg_get_abssplitting" (gen @-> returning gen)

  let alg_get_aut = foreign "alg_get_aut" (gen @-> returning gen)
  let algaut = foreign "algaut" (gen @-> returning gen)
  let alg_get_auts = foreign "alg_get_auts" (gen @-> returning gen)
  let alg_get_b = foreign "alg_get_b" (gen @-> returning gen)
  let algb = foreign "algb" (gen @-> returning gen)
  let algcenter = foreign "algcenter" (gen @-> returning gen)
  let alg_get_center = foreign "alg_get_center" (gen @-> returning gen)
  let alg_get_char = foreign "alg_get_char" (gen @-> returning gen)
  let algchar = foreign "algchar" (gen @-> returning gen)
  let alg_get_degree = foreign "alg_get_degree" (gen @-> returning long)
  let algdegree = foreign "algdegree" (gen @-> returning long)
  let alg_get_dim = foreign "alg_get_dim" (gen @-> returning long)
  let algdim = foreign "algdim" (gen @-> long @-> returning long)
  let alg_get_hasse_f = foreign "alg_get_hasse_f" (gen @-> returning gen)
  let alghassef = foreign "alghassef" (gen @-> returning gen)
  let alg_get_hasse_i = foreign "alg_get_hasse_i" (gen @-> returning gen)
  let alghassei = foreign "alghassei" (gen @-> returning gen)
  let alg_get_invbasis = foreign "alg_get_invbasis" (gen @-> returning gen)
  let alginvbasis = foreign "alginvbasis" (gen @-> returning gen)
  let alg_get_multable = foreign "alg_get_multable" (gen @-> returning gen)
  let alg_get_basis = foreign "alg_get_basis" (gen @-> returning gen)
  let algbasis = foreign "algbasis" (gen @-> returning gen)
  let alg_get_relmultable = foreign "alg_get_relmultable" (gen @-> returning gen)
  let algrelmultable = foreign "algrelmultable" (gen @-> returning gen)
  let alg_get_splitpol = foreign "alg_get_splitpol" (gen @-> returning gen)

  let alg_get_splittingfield =
    foreign "alg_get_splittingfield" (gen @-> returning gen)

  let algsplittingfield = foreign "algsplittingfield" (gen @-> returning gen)

  let alg_get_splittingbasis =
    foreign "alg_get_splittingbasis" (gen @-> returning gen)

  let alg_get_splittingbasisinv =
    foreign "alg_get_splittingbasisinv" (gen @-> returning gen)

  let alg_get_splittingdata =
    foreign "alg_get_splittingdata" (gen @-> returning gen)

  let algsplittingdata = foreign "algsplittingdata" (gen @-> returning gen)
  let alg_get_tracebasis = foreign "alg_get_tracebasis" (gen @-> returning gen)

  let alg_hasse =
    foreign "alg_hasse"
      (gen @-> long @-> gen @-> gen @-> long @-> long @-> returning gen)

  let alg_hilbert =
    foreign "alg_hilbert"
      (gen @-> gen @-> gen @-> long @-> long @-> returning gen)

  let alg_matrix =
    foreign "alg_matrix" (gen @-> long @-> long @-> long @-> returning gen)

  let alg_model = foreign "alg_model" (gen @-> gen @-> returning long)

  let alg_quotient =
    foreign "alg_quotient" (gen @-> gen @-> long @-> returning gen)

  let algradical = foreign "algradical" (gen @-> returning gen)
  let algsimpledec = foreign "algsimpledec" (gen @-> long @-> returning gen)

  let algsimpledec_ss =
    foreign "algsimpledec_ss" (gen @-> long @-> returning gen)

  let algsubalg = foreign "algsubalg" (gen @-> gen @-> returning gen)
  let alg_type = foreign "alg_type" (gen @-> returning long)
  let algadd = foreign "algadd" (gen @-> gen @-> gen @-> returning gen)
  let algalgtobasis = foreign "algalgtobasis" (gen @-> gen @-> returning gen)
  let algbasistoalg = foreign "algbasistoalg" (gen @-> gen @-> returning gen)

  let algcharpoly =
    foreign "algcharpoly" (gen @-> gen @-> long @-> long @-> returning gen)

  let algdisc = foreign "algdisc" (gen @-> returning gen)
  let algdivl = foreign "algdivl" (gen @-> gen @-> gen @-> returning gen)
  let algdivr = foreign "algdivr" (gen @-> gen @-> gen @-> returning gen)
  let alggroup = foreign "alggroup" (gen @-> gen @-> returning gen)

  let alggroupcenter =
    foreign "alggroupcenter" (gen @-> gen @-> ptr gen @-> returning gen)

  let alghasse = foreign "alghasse" (gen @-> gen @-> returning gen)

  let alginit =
    foreign "alginit" (gen @-> gen @-> long @-> long @-> returning gen)

  let algindex = foreign "algindex" (gen @-> gen @-> returning long)
  let alginv = foreign "alginv" (gen @-> gen @-> returning gen)

  let algisassociative =
    foreign "algisassociative" (gen @-> gen @-> returning int)

  let algiscommutative = foreign "algiscommutative" (gen @-> returning int)
  let algisdivision = foreign "algisdivision" (gen @-> gen @-> returning int)
  let algisramified = foreign "algisramified" (gen @-> gen @-> returning int)
  let algissemisimple = foreign "algissemisimple" (gen @-> returning int)
  let algissimple = foreign "algissimple" (gen @-> long @-> returning int)
  let algissplit = foreign "algissplit" (gen @-> gen @-> returning int)

  let algisdivl =
    foreign "algisdivl" (gen @-> gen @-> gen @-> ptr gen @-> returning int)

  let algisinv = foreign "algisinv" (gen @-> gen @-> ptr gen @-> returning int)

  let algmakeintegral =
    foreign "algmakeintegral" (gen @-> long @-> returning gen)

  let algmul = foreign "algmul" (gen @-> gen @-> gen @-> returning gen)
  let algmultable = foreign "algmultable" (gen @-> returning gen)

  let alglat_get_primbasis =
    foreign "alglat_get_primbasis" (gen @-> returning gen)

  let alglat_get_scalar = foreign "alglat_get_scalar" (gen @-> returning gen)

  let alglatadd =
    foreign "alglatadd" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let alglatcontains =
    foreign "alglatcontains" (gen @-> gen @-> gen @-> ptr gen @-> returning int)

  let alglatelement =
    foreign "alglatelement" (gen @-> gen @-> gen @-> returning gen)

  let alglathnf = foreign "alglathnf" (gen @-> gen @-> gen @-> returning gen)
  let alglatindex = foreign "alglatindex" (gen @-> gen @-> gen @-> returning gen)

  let alglatinter =
    foreign "alglatinter" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let alglatmul = foreign "alglatmul" (gen @-> gen @-> gen @-> returning gen)

  let alglatlefttransporter =
    foreign "alglatlefttransporter" (gen @-> gen @-> gen @-> returning gen)

  let alglatrighttransporter =
    foreign "alglatrighttransporter" (gen @-> gen @-> gen @-> returning gen)

  let alglatsubset =
    foreign "alglatsubset" (gen @-> gen @-> gen @-> ptr gen @-> returning int)

  let algneg = foreign "algneg" (gen @-> gen @-> returning gen)
  let algnorm = foreign "algnorm" (gen @-> gen @-> long @-> returning gen)
  let algpoleval = foreign "algpoleval" (gen @-> gen @-> gen @-> returning gen)
  let algpow = foreign "algpow" (gen @-> gen @-> gen @-> returning gen)
  let algprimesubalg = foreign "algprimesubalg" (gen @-> returning gen)
  let algramifiedplaces = foreign "algramifiedplaces" (gen @-> returning gen)
  let algrandom = foreign "algrandom" (gen @-> gen @-> returning gen)
  let algsplit = foreign "algsplit" (gen @-> long @-> returning gen)

  let algtomatrix =
    foreign "algtomatrix" (gen @-> gen @-> long @-> returning gen)

  let algsqr = foreign "algsqr" (gen @-> gen @-> returning gen)
  let algsub = foreign "algsub" (gen @-> gen @-> gen @-> returning gen)
  let algtableinit = foreign "algtableinit" (gen @-> gen @-> returning gen)
  let algtensor = foreign "algtensor" (gen @-> gen @-> long @-> returning gen)
  let algtrace = foreign "algtrace" (gen @-> gen @-> long @-> returning gen)
  let algtype = foreign "algtype" (gen @-> returning long)

  let bnfgwgeneric =
    foreign "bnfgwgeneric"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let checkalg = foreign "checkalg" (gen @-> returning void)

  let checkhasse =
    foreign "checkhasse" (gen @-> gen @-> gen @-> long @-> returning void)
end

module F16 (F : Ctypes.FOREIGN) = struct
  open F

  let checklat = foreign "checklat" (gen @-> gen @-> returning void)

  let conjclasses_algcenter =
    foreign "conjclasses_algcenter" (gen @-> gen @-> returning gen)

  let galoischardet =
    foreign "galoischardet" (gen @-> gen @-> long @-> returning gen)

  let galoischarpoly =
    foreign "galoischarpoly" (gen @-> gen @-> long @-> returning gen)

  let galoischartable = foreign "galoischartable" (gen @-> returning gen)

  let nfgrunwaldwang =
    foreign "nfgrunwaldwang"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let nfgwkummer =
    foreign "nfgwkummer" (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let f2ms_colelim = foreign "F2Ms_colelim" (gen @-> long @-> returning gen)
  let f2m_image = foreign "F2m_image" (gen @-> returning gen)
  let f2m_indexrank = foreign "F2m_indexrank" (gen @-> returning gen)
  let f2m_suppl = foreign "F2m_suppl" (gen @-> returning gen)

  let f2xqm_f2xqc_gauss =
    foreign "F2xqM_F2xqC_gauss" (gen @-> gen @-> gen @-> returning gen)

  let f2xqm_f2xqc_invimage =
    foreign "F2xqM_F2xqC_invimage" (gen @-> gen @-> gen @-> returning gen)

  let f2xqm_f2xqc_mul =
    foreign "F2xqM_F2xqC_mul" (gen @-> gen @-> gen @-> returning gen)

  let f2xqm_deplin = foreign "F2xqM_deplin" (gen @-> gen @-> returning gen)
  let f2xqm_det = foreign "F2xqM_det" (gen @-> gen @-> returning gen)
  let f2xqm_gauss = foreign "F2xqM_gauss" (gen @-> gen @-> gen @-> returning gen)
  let f2xqm_ker = foreign "F2xqM_ker" (gen @-> gen @-> returning gen)
  let f2xqm_image = foreign "F2xqM_image" (gen @-> gen @-> returning gen)
  let f2xqm_indexrank = foreign "F2xqM_indexrank" (gen @-> gen @-> returning gen)
  let f2xqm_inv = foreign "F2xqM_inv" (gen @-> gen @-> returning gen)

  let f2xqm_invimage =
    foreign "F2xqM_invimage" (gen @-> gen @-> gen @-> returning gen)

  let f2xqm_mul = foreign "F2xqM_mul" (gen @-> gen @-> gen @-> returning gen)
  let f2xqm_rank = foreign "F2xqM_rank" (gen @-> gen @-> returning long)
  let f2xqm_suppl = foreign "F2xqM_suppl" (gen @-> gen @-> returning gen)
  let flm_image = foreign "Flm_image" (gen @-> pari_ulong @-> returning gen)

  let flm_indexrank =
    foreign "Flm_indexrank" (gen @-> pari_ulong @-> returning gen)

  let flm_suppl = foreign "Flm_suppl" (gen @-> pari_ulong @-> returning gen)

  let flxqm_flxqc_gauss =
    foreign "FlxqM_FlxqC_gauss"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_flxqc_invimage =
    foreign "FlxqM_FlxqC_invimage"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_flxqc_mul =
    foreign "FlxqM_FlxqC_mul"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_deplin =
    foreign "FlxqM_deplin" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_det =
    foreign "FlxqM_det" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_gauss =
    foreign "FlxqM_gauss" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_ker =
    foreign "FlxqM_ker" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_image =
    foreign "FlxqM_image" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_indexrank =
    foreign "FlxqM_indexrank" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_inv =
    foreign "FlxqM_inv" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_invimage =
    foreign "FlxqM_invimage"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_mul =
    foreign "FlxqM_mul" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqm_rank =
    foreign "FlxqM_rank" (gen @-> gen @-> pari_ulong @-> returning long)

  let flxqm_suppl =
    foreign "FlxqM_suppl" (gen @-> gen @-> pari_ulong @-> returning gen)

  let fpm_fpc_gauss =
    foreign "FpM_FpC_gauss" (gen @-> gen @-> gen @-> returning gen)

  let fpm_fpc_invimage =
    foreign "FpM_FpC_invimage" (gen @-> gen @-> gen @-> returning gen)

  let fpm_deplin = foreign "FpM_deplin" (gen @-> gen @-> returning gen)
  let fpm_det = foreign "FpM_det" (gen @-> gen @-> returning gen)
  let fpm_gauss = foreign "FpM_gauss" (gen @-> gen @-> gen @-> returning gen)
  let fpm_image = foreign "FpM_image" (gen @-> gen @-> returning gen)
  let fpm_indexrank = foreign "FpM_indexrank" (gen @-> gen @-> returning gen)

  let fpm_intersect =
    foreign "FpM_intersect" (gen @-> gen @-> gen @-> returning gen)

  let fpm_intersect_i =
    foreign "FpM_intersect_i" (gen @-> gen @-> gen @-> returning gen)

  let fpm_inv = foreign "FpM_inv" (gen @-> gen @-> returning gen)

  let fpm_invimage =
    foreign "FpM_invimage" (gen @-> gen @-> gen @-> returning gen)

  let fpm_ker = foreign "FpM_ker" (gen @-> gen @-> returning gen)
  let fpm_rank = foreign "FpM_rank" (gen @-> gen @-> returning long)
  let fpm_suppl = foreign "FpM_suppl" (gen @-> gen @-> returning gen)

  let fqm_fqc_gauss =
    foreign "FqM_FqC_gauss" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqm_fqc_invimage =
    foreign "FqM_FqC_invimage" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqm_fqc_mul =
    foreign "FqM_FqC_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqm_deplin = foreign "FqM_deplin" (gen @-> gen @-> gen @-> returning gen)
  let fqm_det = foreign "FqM_det" (gen @-> gen @-> gen @-> returning gen)

  let fqm_gauss =
    foreign "FqM_gauss" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqm_ker = foreign "FqM_ker" (gen @-> gen @-> gen @-> returning gen)
  let fqm_image = foreign "FqM_image" (gen @-> gen @-> gen @-> returning gen)

  let fqm_indexrank =
    foreign "FqM_indexrank" (gen @-> gen @-> gen @-> returning gen)

  let fqm_inv = foreign "FqM_inv" (gen @-> gen @-> gen @-> returning gen)

  let fqm_invimage =
    foreign "FqM_invimage" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqm_mul = foreign "FqM_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fqm_rank = foreign "FqM_rank" (gen @-> gen @-> gen @-> returning long)
  let fqm_suppl = foreign "FqM_suppl" (gen @-> gen @-> gen @-> returning gen)
  let qm_image_shallow = foreign "QM_image_shallow" (gen @-> returning gen)
  let qm_image = foreign "QM_image" (gen @-> returning gen)
  let qm_gauss = foreign "QM_gauss" (gen @-> gen @-> returning gen)
  let qm_gauss_i = foreign "QM_gauss_i" (gen @-> gen @-> long @-> returning gen)
  let qm_indexrank = foreign "QM_indexrank" (gen @-> returning gen)
  let qm_inv = foreign "QM_inv" (gen @-> returning gen)
  let qm_rank = foreign "QM_rank" (gen @-> returning long)

  let rgm_fp_init =
    foreign "RgM_Fp_init" (gen @-> gen @-> ptr pari_ulong @-> returning gen)

  let rgm_hadamard = foreign "RgM_Hadamard" (gen @-> returning gen)

  let rgm_rgc_invimage =
    foreign "RgM_RgC_invimage" (gen @-> gen @-> returning gen)

  let rgm_diagonal = foreign "RgM_diagonal" (gen @-> returning gen)

  let rgm_diagonal_shallow =
    foreign "RgM_diagonal_shallow" (gen @-> returning gen)

  let rgm_inv = foreign "RgM_inv" (gen @-> returning gen)
  let rgm_inv_upper = foreign "RgM_inv_upper" (gen @-> returning gen)
  let rgm_invimage = foreign "RgM_invimage" (gen @-> gen @-> returning gen)
  let rgm_solve = foreign "RgM_solve" (gen @-> gen @-> returning gen)

  let rgm_solve_realimag =
    foreign "RgM_solve_realimag" (gen @-> gen @-> returning gen)

  let rgms_structelim =
    foreign "RgMs_structelim"
      (gen @-> long @-> gen @-> ptr gen @-> ptr gen @-> returning void)

  let zm_det = foreign "ZM_det" (gen @-> returning gen)
  let zm_detmult = foreign "ZM_detmult" (gen @-> returning gen)
  let zm_gauss = foreign "ZM_gauss" (gen @-> gen @-> returning gen)
  let zm_ker = foreign "ZM_ker" (gen @-> returning gen)
  let zm_imagecompl = foreign "ZM_imagecompl" (gen @-> returning gen)
  let zm_indeximage = foreign "ZM_indeximage" (gen @-> returning gen)
  let zm_indexrank = foreign "ZM_indexrank" (gen @-> returning gen)
  let zm_inv = foreign "ZM_inv" (gen @-> ptr gen @-> returning gen)

  let zm_inv_ratlift =
    foreign "ZM_inv_ratlift" (gen @-> ptr gen @-> returning gen)

  let zm_pseudoinv =
    foreign "ZM_pseudoinv" (gen @-> ptr gen @-> ptr gen @-> returning gen)

  let zm_rank = foreign "ZM_rank" (gen @-> returning long)

  let zlm_gauss =
    foreign "ZlM_gauss"
      (gen @-> gen @-> pari_ulong @-> long @-> gen @-> returning gen)
end

module F17 (F : Ctypes.FOREIGN) = struct
  open F

  let closemodinvertible =
    foreign "closemodinvertible" (gen @-> gen @-> returning gen)

  let deplin = foreign "deplin" (gen @-> returning gen)
  let det = foreign "det" (gen @-> returning gen)
  let det0 = foreign "det0" (gen @-> long @-> returning gen)
  let det2 = foreign "det2" (gen @-> returning gen)
  let detint = foreign "detint" (gen @-> returning gen)
  let eigen = foreign "eigen" (gen @-> long @-> returning gen)
  let gauss = foreign "gauss" (gen @-> gen @-> returning gen)
  let gaussmodulo = foreign "gaussmodulo" (gen @-> gen @-> gen @-> returning gen)

  let gaussmodulo2 =
    foreign "gaussmodulo2" (gen @-> gen @-> gen @-> returning gen)

  let gen_gauss =
    foreign "gen_Gauss"
      (gen @-> gen @-> ptr void @-> ptr bb_field @-> returning gen)

  let gen_gauss_pivot =
    foreign "gen_Gauss_pivot"
      (gen @-> ptr long @-> ptr void @-> ptr bb_field @-> returning gen)

  let gen_det =
    foreign "gen_det" (gen @-> ptr void @-> ptr bb_field @-> returning gen)

  let gen_ker =
    foreign "gen_ker"
      (gen @-> long @-> ptr void @-> ptr bb_field @-> returning gen)

  let gen_matcolinvimage =
    foreign "gen_matcolinvimage"
      (gen @-> gen @-> ptr void @-> ptr bb_field @-> returning gen)

  let gen_matcolmul =
    foreign "gen_matcolmul"
      (gen @-> gen @-> ptr void @-> ptr bb_field @-> returning gen)

  let gen_matinvimage =
    foreign "gen_matinvimage"
      (gen @-> gen @-> ptr void @-> ptr bb_field @-> returning gen)

  let gen_matmul =
    foreign "gen_matmul"
      (gen @-> gen @-> ptr void @-> ptr bb_field @-> returning gen)

  let image = foreign "image" (gen @-> returning gen)
  let image2 = foreign "image2" (gen @-> returning gen)
  let imagecompl = foreign "imagecompl" (gen @-> returning gen)
  let indexrank = foreign "indexrank" (gen @-> returning gen)
  let inverseimage = foreign "inverseimage" (gen @-> gen @-> returning gen)
  let ker = foreign "ker" (gen @-> returning gen)
  let mateigen = foreign "mateigen" (gen @-> long @-> long @-> returning gen)
  let matimage0 = foreign "matimage0" (gen @-> long @-> returning gen)
  let matker0 = foreign "matker0" (gen @-> long @-> returning gen)
  let rank = foreign "rank" (gen @-> returning long)

  let reducemodinvertible =
    foreign "reducemodinvertible" (gen @-> gen @-> returning gen)

  let reducemodlll = foreign "reducemodlll" (gen @-> gen @-> returning gen)

  let split_realimag =
    foreign "split_realimag" (gen @-> long @-> long @-> returning gen)

  let suppl = foreign "suppl" (gen @-> returning gen)

  let flm_charpoly =
    foreign "Flm_charpoly" (gen @-> pari_ulong @-> returning gen)

  let flm_hess = foreign "Flm_hess" (gen @-> pari_ulong @-> returning gen)
  let fpm_charpoly = foreign "FpM_charpoly" (gen @-> gen @-> returning gen)
  let fpm_hess = foreign "FpM_hess" (gen @-> gen @-> returning gen)
  let frobeniusform = foreign "Frobeniusform" (gen @-> long @-> returning gen)

  let qm_minors_coprime =
    foreign "QM_minors_coprime" (gen @-> gen @-> returning gen)

  let qm_imz = foreign "QM_ImZ" (gen @-> returning gen)

  let qm_imz_all =
    foreign "QM_ImZ_all" (gen @-> ptr gen @-> long @-> long @-> returning gen)

  let qm_imz_hnf = foreign "QM_ImZ_hnf" (gen @-> returning gen)

  let qm_imz_hnfall =
    foreign "QM_ImZ_hnfall" (gen @-> ptr gen @-> long @-> returning gen)

  let qm_imq = foreign "QM_ImQ" (gen @-> returning gen)

  let qm_imq_all =
    foreign "QM_ImQ_all" (gen @-> ptr gen @-> long @-> long @-> returning gen)

  let qm_imq_hnf = foreign "QM_ImQ_hnf" (gen @-> returning gen)

  let qm_imq_hnfall =
    foreign "QM_ImQ_hnfall" (gen @-> ptr gen @-> long @-> returning gen)

  let qm_charpoly_zx = foreign "QM_charpoly_ZX" (gen @-> returning gen)

  let qm_charpoly_zx_bound =
    foreign "QM_charpoly_ZX_bound" (gen @-> long @-> returning gen)

  let zm_charpoly = foreign "ZM_charpoly" (gen @-> returning gen)
  let adj = foreign "adj" (gen @-> returning gen)
  let adjsafe = foreign "adjsafe" (gen @-> returning gen)
  let caract = foreign "caract" (gen @-> long @-> returning gen)
  let caradj = foreign "caradj" (gen @-> long @-> ptr gen @-> returning gen)
  let carberkowitz = foreign "carberkowitz" (gen @-> long @-> returning gen)
  let carhess = foreign "carhess" (gen @-> long @-> returning gen)
  let charpoly = foreign "charpoly" (gen @-> long @-> returning gen)
  let charpoly0 = foreign "charpoly0" (gen @-> long @-> long @-> returning gen)
  let gnorm = foreign "gnorm" (gen @-> returning gen)
  let gnorml1 = foreign "gnorml1" (gen @-> long @-> returning gen)
  let gnorml1_fake = foreign "gnorml1_fake" (gen @-> returning gen)
  let gnormlp = foreign "gnormlp" (gen @-> gen @-> long @-> returning gen)
  let gnorml2 = foreign "gnorml2" (gen @-> returning gen)
  let gsupnorm = foreign "gsupnorm" (gen @-> long @-> returning gen)

  let gsupnorm_aux =
    foreign "gsupnorm_aux"
      (gen @-> ptr gen @-> ptr gen @-> long @-> returning void)

  let gtrace = foreign "gtrace" (gen @-> returning gen)
  let hess = foreign "hess" (gen @-> returning gen)
  let intersect = foreign "intersect" (gen @-> gen @-> returning gen)
  let jacobi = foreign "jacobi" (gen @-> long @-> returning gen)
  let matadjoint0 = foreign "matadjoint0" (gen @-> long @-> returning gen)
  let matcompanion = foreign "matcompanion" (gen @-> returning gen)
  let matrixqz0 = foreign "matrixqz0" (gen @-> gen @-> returning gen)
  let minpoly = foreign "minpoly" (gen @-> long @-> returning gen)
  let qfgaussred = foreign "qfgaussred" (gen @-> returning gen)
  let qfgaussred_positive = foreign "qfgaussred_positive" (gen @-> returning gen)
  let qfsign = foreign "qfsign" (gen @-> returning gen)
  let apply0 = foreign "apply0" (gen @-> gen @-> returning gen)
  let diagonal = foreign "diagonal" (gen @-> returning gen)
  let diagonal_shallow = foreign "diagonal_shallow" (gen @-> returning gen)
  let extract0 = foreign "extract0" (gen @-> gen @-> gen @-> returning gen)
  let fold0 = foreign "fold0" (gen @-> gen @-> returning gen)

  let genapply =
    foreign "genapply"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> returning gen)

  let genfold =
    foreign "genfold"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> gen @-> returning gen)

  let genindexselect =
    foreign "genindexselect"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> gen @-> returning gen)

  let genselect =
    foreign "genselect"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> gen @-> returning gen)

  let gtomat = foreign "gtomat" (gen @-> returning gen)
  let gtrans = foreign "gtrans" (gen @-> returning gen)
  let matmuldiagonal = foreign "matmuldiagonal" (gen @-> gen @-> returning gen)

  let matmultodiagonal =
    foreign "matmultodiagonal" (gen @-> gen @-> returning gen)

  let matslice0 =
    foreign "matslice0"
      (gen @-> long @-> long @-> long @-> long @-> returning gen)

  let parapply = foreign "parapply" (gen @-> gen @-> returning gen)

  let parfor =
    foreign "parfor"
      (gen @-> gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning long)
      @-> returning void)

  let parfor_init =
    foreign "parfor_init"
      (ptr parfor_t @-> gen @-> gen @-> gen @-> returning void)

  let parfor_next = foreign "parfor_next" (ptr parfor_t @-> returning gen)
  let parfor_stop = foreign "parfor_stop" (ptr parfor_t @-> returning void)

  let parforeach =
    foreign "parforeach"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning long)
      @-> returning void)

  let parforeach_init =
    foreign "parforeach_init"
      (ptr parforeach_t @-> gen @-> gen @-> returning void)

  let parforeach_next =
    foreign "parforeach_next" (ptr parforeach_t @-> returning gen)

  let parforeach_stop =
    foreign "parforeach_stop" (ptr parforeach_t @-> returning void)

  let parforprime =
    foreign "parforprime"
      (gen @-> gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning long)
      @-> returning void)

  let parforprime_init =
    foreign "parforprime_init"
      (ptr parforprime_t @-> gen @-> gen @-> gen @-> returning void)
end

module F18 (F : Ctypes.FOREIGN) = struct
  open F

  let parforprime_next =
    foreign "parforprime_next" (ptr parforprime_t @-> returning gen)

  let parforprime_stop =
    foreign "parforprime_stop" (ptr parforprime_t @-> returning void)

  let parforprimestep =
    foreign "parforprimestep"
      (gen @-> gen @-> gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning long)
      @-> returning void)

  let parforprimestep_init =
    foreign "parforprimestep_init"
      (ptr parforprime_t @-> gen @-> gen @-> gen @-> gen @-> returning void)

  let parforvec =
    foreign "parforvec"
      (gen @-> gen @-> long @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning long)
      @-> returning void)

  let parforvec_init =
    foreign "parforvec_init"
      (ptr parforvec_t @-> gen @-> gen @-> long @-> returning void)

  let parforvec_next =
    foreign "parforvec_next" (ptr parforvec_t @-> returning gen)

  let parforvec_stop =
    foreign "parforvec_stop" (ptr parforvec_t @-> returning void)

  let parselect = foreign "parselect" (gen @-> gen @-> long @-> returning gen)
  let select0 = foreign "select0" (gen @-> gen @-> long @-> returning gen)
  let shallowextract = foreign "shallowextract" (gen @-> gen @-> returning gen)

  let shallowmatextract =
    foreign "shallowmatextract" (gen @-> gen @-> gen @-> returning gen)

  let shallowtrans = foreign "shallowtrans" (gen @-> returning gen)

  let vecapply =
    foreign "vecapply"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> returning gen)

  let veccatapply =
    foreign "veccatapply"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> returning gen)

  let veccatselapply =
    foreign "veccatselapply"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> returning gen)

  let vecrange = foreign "vecrange" (gen @-> gen @-> returning gen)
  let vecrangess = foreign "vecrangess" (long @-> long @-> returning gen)

  let vecselapply =
    foreign "vecselapply"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> returning gen)

  let vecselect =
    foreign "vecselect"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> gen @-> returning gen)

  let vecslice0 = foreign "vecslice0" (gen @-> long @-> long @-> returning gen)
  let vecsum = foreign "vecsum" (gen @-> returning gen)
  let zv_diagonal = foreign "zv_diagonal" (gen @-> returning gen)
  let addhelp = foreign "addhelp" (string @-> string @-> returning void)
  let arity0 = foreign "arity0" (gen @-> returning gen)
  let alias0 = foreign "alias0" (string @-> string @-> returning void)
  let compile_str = foreign "compile_str" (string @-> returning gen)
  let delete_var = foreign "delete_var" (void @-> returning long)
  let fetch_user_var = foreign "fetch_user_var" (string @-> returning long)
  let fetch_var = foreign "fetch_var" (void @-> returning long)
  let fetch_var_higher = foreign "fetch_var_higher" (void @-> returning long)

  let fetch_var_value =
    foreign "fetch_var_value" (long @-> gen @-> returning gen)

  let gp_embedded = foreign "gp_embedded" (string @-> returning string)

  let gp_embedded_init =
    foreign "gp_embedded_init" (long @-> long @-> returning void)

  let gp_read_str = foreign "gp_read_str" (string @-> returning gen)

  let gp_read_str_bitprec =
    foreign "gp_read_str_bitprec" (string @-> long @-> returning gen)

  let gp_read_str_prec =
    foreign "gp_read_str_prec" (string @-> long @-> returning gen)

  let install =
    foreign "install" (ptr void @-> string @-> string @-> returning (ptr entree))

  let is_entry = foreign "is_entry" (string @-> returning (ptr entree))
  let kill0 = foreign "kill0" (string @-> returning void)
  let pari_var_close = foreign "pari_var_close" (void @-> returning void)
  let pari_var_init = foreign "pari_var_init" (void @-> returning void)
  let pari_var_next = foreign "pari_var_next" (void @-> returning long)
  let pari_var_next_temp = foreign "pari_var_next_temp" (void @-> returning long)
  let pari_var_create = foreign "pari_var_create" (ptr entree @-> returning long)
  let name_var = foreign "name_var" (long @-> string @-> returning void)
  let readseq = foreign "readseq" (string @-> returning gen)
  let safegel = foreign "safegel" (gen @-> long @-> returning (ptr gen))
  let safeel = foreign "safeel" (gen @-> long @-> returning (ptr long))
  let safelistel = foreign "safelistel" (gen @-> long @-> returning (ptr gen))

  let safegcoeff =
    foreign "safegcoeff" (gen @-> long @-> long @-> returning (ptr gen))

  let strtoi = foreign "strtoi" (string @-> returning gen)
  let strtor = foreign "strtor" (string @-> long @-> returning gen)
  let varhigher = foreign "varhigher" (string @-> long @-> returning gen)
  let varlower = foreign "varlower" (string @-> long @-> returning gen)

  let divisorslenstra =
    foreign "divisorslenstra" (gen @-> gen @-> gen @-> returning gen)

  let isprimeaprcl = foreign "isprimeAPRCL" (gen @-> returning long)
  let qfb0 = foreign "Qfb0" (gen @-> gen @-> gen @-> returning gen)

  let check_quaddisc =
    foreign "check_quaddisc"
      (gen @-> ptr long @-> ptr long @-> string @-> returning void)

  let check_quaddisc_imag =
    foreign "check_quaddisc_imag"
      (gen @-> ptr long @-> string @-> returning void)

  let check_quaddisc_real =
    foreign "check_quaddisc_real"
      (gen @-> ptr long @-> string @-> returning void)

  let cornacchia =
    foreign "cornacchia" (gen @-> gen @-> ptr gen @-> ptr gen @-> returning long)

  let cornacchia2 =
    foreign "cornacchia2"
      (gen @-> gen @-> ptr gen @-> ptr gen @-> returning long)

  let cornacchia2_sqrt =
    foreign "cornacchia2_sqrt"
      (gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning long)

  let nucomp = foreign "nucomp" (gen @-> gen @-> gen @-> returning gen)
  let nudupl = foreign "nudupl" (gen @-> gen @-> returning gen)
  let nupow = foreign "nupow" (gen @-> gen @-> gen @-> returning gen)
  let primeform = foreign "primeform" (gen @-> gen @-> returning gen)
  let primeform_u = foreign "primeform_u" (gen @-> pari_ulong @-> returning gen)
  let qfb_1 = foreign "qfb_1" (gen @-> returning gen)
  let qfb_equal1 = foreign "qfb_equal1" (gen @-> returning int)
  let qfbcomp = foreign "qfbcomp" (gen @-> gen @-> returning gen)
  let qfbcomp_i = foreign "qfbcomp_i" (gen @-> gen @-> returning gen)
  let qfbcompraw = foreign "qfbcompraw" (gen @-> gen @-> returning gen)
  let qfbcompraw_i = foreign "qfbcompraw_i" (gen @-> gen @-> returning gen)
  let qfbcornacchia = foreign "qfbcornacchia" (gen @-> gen @-> returning gen)
  let qfbpow = foreign "qfbpow" (gen @-> gen @-> returning gen)
  let qfbpow_i = foreign "qfbpow_i" (gen @-> gen @-> returning gen)
  let qfbpowraw = foreign "qfbpowraw" (gen @-> long @-> returning gen)
  let qfbpows = foreign "qfbpows" (gen @-> long @-> returning gen)
  let qfbred = foreign "qfbred" (gen @-> returning gen)
  let qfbred_i = foreign "qfbred_i" (gen @-> returning gen)

  let qfbred0 =
    foreign "qfbred0" (gen @-> long @-> gen @-> gen @-> returning gen)

  let qfbredsl2 = foreign "qfbredsl2" (gen @-> gen @-> returning gen)
  let qfbsolve = foreign "qfbsolve" (gen @-> gen @-> long @-> returning gen)
  let qfbsqr = foreign "qfbsqr" (gen @-> returning gen)
  let qfbsqr_i = foreign "qfbsqr_i" (gen @-> returning gen)
  let qfi_shanks = foreign "qfi_Shanks" (gen @-> gen @-> long @-> returning gen)
  let qfi_log = foreign "qfi_log" (gen @-> gen @-> gen @-> returning gen)
  let qfi_order = foreign "qfi_order" (gen @-> gen @-> returning gen)
  let qfisolvep = foreign "qfisolvep" (gen @-> gen @-> returning gen)

  let qfr3_comp =
    foreign "qfr3_comp" (gen @-> gen @-> ptr qfr_data @-> returning gen)

  let qfr3_compraw = foreign "qfr3_compraw" (gen @-> gen @-> returning gen)

  let qfr3_pow =
    foreign "qfr3_pow" (gen @-> gen @-> ptr qfr_data @-> returning gen)

  let qfr3_red = foreign "qfr3_red" (gen @-> ptr qfr_data @-> returning gen)
  let qfr3_rho = foreign "qfr3_rho" (gen @-> ptr qfr_data @-> returning gen)
  let qfr3_to_qfr = foreign "qfr3_to_qfr" (gen @-> gen @-> returning gen)

  let qfr5_comp =
    foreign "qfr5_comp" (gen @-> gen @-> ptr qfr_data @-> returning gen)

  let qfr5_compraw = foreign "qfr5_compraw" (gen @-> gen @-> returning gen)
  let qfr5_dist = foreign "qfr5_dist" (gen @-> gen @-> long @-> returning gen)
end

module F19 (F : Ctypes.FOREIGN) = struct
  open F

  let qfr5_pow =
    foreign "qfr5_pow" (gen @-> gen @-> ptr qfr_data @-> returning gen)

  let qfr5_red = foreign "qfr5_red" (gen @-> ptr qfr_data @-> returning gen)
  let qfr5_rho = foreign "qfr5_rho" (gen @-> ptr qfr_data @-> returning gen)
  let qfr5_to_qfr = foreign "qfr5_to_qfr" (gen @-> gen @-> gen @-> returning gen)

  let qfr_data_init =
    foreign "qfr_data_init" (gen @-> long @-> ptr qfr_data @-> returning void)

  let qfr_to_qfr5 = foreign "qfr_to_qfr5" (gen @-> long @-> returning gen)
  let qfrsolvep = foreign "qfrsolvep" (gen @-> gen @-> returning gen)
  let quadgen = foreign "quadgen" (gen @-> returning gen)
  let quadgen0 = foreign "quadgen0" (gen @-> long @-> returning gen)
  let quadpoly = foreign "quadpoly" (gen @-> returning gen)
  let quadpoly_i = foreign "quadpoly_i" (gen @-> returning gen)
  let quadpoly0 = foreign "quadpoly0" (gen @-> long @-> returning gen)

  let fl_2gener_pre =
    foreign "Fl_2gener_pre" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_log =
    foreign "Fl_log"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fl_log_pre =
    foreign "Fl_log_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fl_order =
    foreign "Fl_order"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_powers =
    foreign "Fl_powers" (pari_ulong @-> long @-> pari_ulong @-> returning gen)

  let fl_powers_pre =
    foreign "Fl_powers_pre"
      (pari_ulong @-> long @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fl_powu =
    foreign "Fl_powu"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_powu_pre =
    foreign "Fl_powu_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fl_sqrt =
    foreign "Fl_sqrt" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_sqrt_pre =
    foreign "Fl_sqrt_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_sqrt_pre_i =
    foreign "Fl_sqrt_pre_i"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fl_sqrtl =
    foreign "Fl_sqrtl"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_sqrtl_pre =
    foreign "Fl_sqrtl_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fl_sqrtn =
    foreign "Fl_sqrtn"
      (pari_ulong @-> long @-> pari_ulong @-> ptr pari_ulong
     @-> returning pari_ulong)

  let fl_sqrtn_pre =
    foreign "Fl_sqrtn_pre"
      (pari_ulong @-> long @-> pari_ulong @-> pari_ulong @-> ptr pari_ulong
     @-> returning pari_ulong)

  let fp_2gener = foreign "Fp_2gener" (gen @-> returning gen)

  let fp_factored_order =
    foreign "Fp_factored_order" (gen @-> gen @-> gen @-> returning gen)

  let fp_ispower = foreign "Fp_ispower" (gen @-> gen @-> gen @-> returning int)
  let fp_log = foreign "Fp_log" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fp_order = foreign "Fp_order" (gen @-> gen @-> gen @-> returning gen)
  let fp_pow = foreign "Fp_pow" (gen @-> gen @-> gen @-> returning gen)

  let fp_pow_init =
    foreign "Fp_pow_init" (gen @-> gen @-> long @-> gen @-> returning gen)

  let fp_pow_table =
    foreign "Fp_pow_table" (gen @-> gen @-> gen @-> returning gen)

  let fp_powers = foreign "Fp_powers" (gen @-> long @-> gen @-> returning gen)
  let fp_pows = foreign "Fp_pows" (gen @-> long @-> gen @-> returning gen)
  let fp_powu = foreign "Fp_powu" (gen @-> pari_ulong @-> gen @-> returning gen)
  let fp_sqrt = foreign "Fp_sqrt" (gen @-> gen @-> returning gen)
  let fp_sqrt_i = foreign "Fp_sqrt_i" (gen @-> gen @-> gen @-> returning gen)

  let fp_sqrtn =
    foreign "Fp_sqrtn" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fpv_prod = foreign "FpV_prod" (gen @-> gen @-> returning gen)
  let z_zv_mod = foreign "Z_ZV_mod" (gen @-> gen @-> returning gen)

  let z_zv_mod_tree =
    foreign "Z_ZV_mod_tree" (gen @-> gen @-> gen @-> returning gen)

  let z_chinese =
    foreign "Z_chinese" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let z_chinese_all =
    foreign "Z_chinese_all"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let z_chinese_coprime =
    foreign "Z_chinese_coprime"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let z_chinese_post =
    foreign "Z_chinese_post"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let z_chinese_pre =
    foreign "Z_chinese_pre"
      (gen @-> gen @-> ptr gen @-> ptr gen @-> ptr gen @-> returning void)

  let z_factor_listp = foreign "Z_factor_listP" (gen @-> gen @-> returning gen)
  let z_nv_mod = foreign "Z_nv_mod" (gen @-> gen @-> returning gen)

  let zm_nv_mod_tree =
    foreign "ZM_nv_mod_tree" (gen @-> gen @-> gen @-> returning gen)

  let zv_allpnqn = foreign "ZV_allpnqn" (gen @-> returning gen)

  let zv_chinese =
    foreign "ZV_chinese" (gen @-> gen @-> ptr gen @-> returning gen)

  let zv_chinese_tree =
    foreign "ZV_chinese_tree" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let zv_chinesetree = foreign "ZV_chinesetree" (gen @-> gen @-> returning gen)

  let zv_nv_mod_tree =
    foreign "ZV_nv_mod_tree" (gen @-> gen @-> gen @-> returning gen)

  let zv_producttree = foreign "ZV_producttree" (gen @-> returning gen)

  let zx_nv_mod_tree =
    foreign "ZX_nv_mod_tree" (gen @-> gen @-> gen @-> returning gen)

  let zxc_nv_mod_tree =
    foreign "ZXC_nv_mod_tree" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zxm_nv_mod_tree =
    foreign "ZXM_nv_mod_tree" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zxx_nv_mod_tree =
    foreign "ZXX_nv_mod_tree" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zideallog = foreign "Zideallog" (gen @-> gen @-> returning gen)
  let bestappr = foreign "bestappr" (gen @-> gen @-> returning gen)
  let bestapprpade = foreign "bestapprPade" (gen @-> long @-> returning gen)
  let chinese = foreign "chinese" (gen @-> gen @-> returning gen)
  let chinese1 = foreign "chinese1" (gen @-> returning gen)
  let chinese1_coprime_z = foreign "chinese1_coprime_Z" (gen @-> returning gen)
  let contfrac0 = foreign "contfrac0" (gen @-> gen @-> long @-> returning gen)
  let contfracpnqn = foreign "contfracpnqn" (gen @-> long @-> returning gen)
  let fibo = foreign "fibo" (long @-> returning gen)
  let gboundcf = foreign "gboundcf" (gen @-> long @-> returning gen)
  let gcf = foreign "gcf" (gen @-> returning gen)
  let gcf2 = foreign "gcf2" (gen @-> gen @-> returning gen)

  let get_fp_field =
    foreign "get_Fp_field" (ptr (ptr void) @-> gen @-> returning (ptr bb_field))

  let hilbert = foreign "hilbert" (gen @-> gen @-> gen @-> returning long)
  let hilbertii = foreign "hilbertii" (gen @-> gen @-> gen @-> returning long)
  let istotient = foreign "istotient" (gen @-> ptr gen @-> returning long)
  let krois = foreign "krois" (gen @-> long @-> returning long)
  let kroiu = foreign "kroiu" (gen @-> pari_ulong @-> returning long)
  let kronecker = foreign "kronecker" (gen @-> gen @-> returning long)
  let krosi = foreign "krosi" (long @-> gen @-> returning long)
  let kross = foreign "kross" (long @-> long @-> returning long)
  let kroui = foreign "kroui" (pari_ulong @-> gen @-> returning long)
  let krouu = foreign "krouu" (pari_ulong @-> pari_ulong @-> returning long)
  let lcmii = foreign "lcmii" (gen @-> gen @-> returning gen)
  let fp_invgen = foreign "Fp_invgen" (gen @-> gen @-> ptr gen @-> returning gen)
  let logint0 = foreign "logint0" (gen @-> gen @-> ptr gen @-> returning long)

  let logintall =
    foreign "logintall" (gen @-> gen @-> ptr gen @-> returning long)

  let mpfact = foreign "mpfact" (long @-> returning gen)

  let factorial_fl =
    foreign "factorial_Fl" (long @-> pari_ulong @-> returning pari_ulong)

  let factorial_fp = foreign "factorial_Fp" (long @-> gen @-> returning gen)
  let muls_interval = foreign "muls_interval" (long @-> long @-> returning gen)

  let mulu_interval =
    foreign "mulu_interval" (pari_ulong @-> pari_ulong @-> returning gen)

  let mulu_interval_step =
    foreign "mulu_interval_step"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let ncv_chinese_center =
    foreign "ncV_chinese_center" (gen @-> gen @-> ptr gen @-> returning gen)

  let ncv_chinese_center_tree =
    foreign "ncV_chinese_center_tree"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let nmv_chinese_center =
    foreign "nmV_chinese_center" (gen @-> gen @-> ptr gen @-> returning gen)

  let nmv_chinese_center_tree =
    foreign "nmV_chinese_center_tree"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let nonsquare_fl = foreign "nonsquare_Fl" (pari_ulong @-> returning pari_ulong)
end

module F20 (F : Ctypes.FOREIGN) = struct
  open F

  let nxcv_chinese_center =
    foreign "nxCV_chinese_center" (gen @-> gen @-> ptr gen @-> returning gen)

  let nxcv_chinese_center_tree =
    foreign "nxCV_chinese_center_tree"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let nxmv_chinese_center =
    foreign "nxMV_chinese_center" (gen @-> gen @-> ptr gen @-> returning gen)

  let nxv_chinese_center =
    foreign "nxV_chinese_center" (gen @-> gen @-> ptr gen @-> returning gen)

  let nxv_chinese_center_tree =
    foreign "nxV_chinese_center_tree"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let zv_chinese_center =
    foreign "ZV_chinese_center" (gen @-> gen @-> ptr gen @-> returning gen)

  let odd_prime_divisors = foreign "odd_prime_divisors" (gen @-> returning gen)
  let pgener_fl = foreign "pgener_Fl" (pari_ulong @-> returning pari_ulong)

  let pgener_fl_local =
    foreign "pgener_Fl_local" (pari_ulong @-> gen @-> returning pari_ulong)

  let pgener_fp = foreign "pgener_Fp" (gen @-> returning gen)
  let pgener_fp_local = foreign "pgener_Fp_local" (gen @-> gen @-> returning gen)
  let pgener_zl = foreign "pgener_Zl" (pari_ulong @-> returning pari_ulong)
  let pgener_zp = foreign "pgener_Zp" (gen @-> returning gen)
  let pnqn = foreign "pnqn" (gen @-> returning gen)
  let ramanujantau = foreign "ramanujantau" (gen @-> long @-> returning gen)

  let rootsof1_fl =
    foreign "rootsof1_Fl" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let rootsof1_fp = foreign "rootsof1_Fp" (gen @-> gen @-> returning gen)

  let rootsof1u_fp =
    foreign "rootsof1u_Fp" (pari_ulong @-> gen @-> returning gen)

  let u_chinese_coprime =
    foreign "u_chinese_coprime"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let znlog = foreign "znlog" (gen @-> gen @-> gen @-> returning gen)
  let znorder = foreign "znorder" (gen @-> gen @-> returning gen)
  let znprimroot = foreign "znprimroot" (gen @-> returning gen)
  let znstar = foreign "znstar" (gen @-> returning gen)
  let znstar0 = foreign "znstar0" (gen @-> long @-> returning gen)
  let rgv_is_zvpos = foreign "RgV_is_ZVpos" (gen @-> returning int)
  let rgv_is_zvnon0 = foreign "RgV_is_ZVnon0" (gen @-> returning int)
  let rgv_is_prv = foreign "RgV_is_prV" (gen @-> returning int)

  let z_issquarefree_fact =
    foreign "Z_issquarefree_fact" (gen @-> returning long)

  let z_lsmoothen =
    foreign "Z_lsmoothen" (gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let z_smoothen =
    foreign "Z_smoothen" (gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let bigomega = foreign "bigomega" (gen @-> returning long)
  let bigomegau = foreign "bigomegau" (pari_ulong @-> returning long)
  let boundfact = foreign "boundfact" (gen @-> pari_ulong @-> returning gen)

  let check_arith_pos =
    foreign "check_arith_pos" (gen @-> string @-> returning gen)

  let check_arith_non0 =
    foreign "check_arith_non0" (gen @-> string @-> returning gen)

  let check_arith_all =
    foreign "check_arith_all" (gen @-> string @-> returning gen)

  let clean_z_factor = foreign "clean_Z_factor" (gen @-> returning gen)
  let core = foreign "core" (gen @-> returning gen)

  let coredisc2_fact =
    foreign "coredisc2_fact"
      (gen @-> long @-> ptr gen @-> ptr gen @-> returning gen)

  let coredisc2u_fact =
    foreign "coredisc2u_fact"
      (gen @-> long @-> ptr gen @-> ptr gen @-> returning pari_ulong)

  let corepartial = foreign "corepartial" (gen @-> long @-> returning gen)
  let core0 = foreign "core0" (gen @-> long @-> returning gen)
  let core2 = foreign "core2" (gen @-> returning gen)
  let core2partial = foreign "core2partial" (gen @-> long @-> returning gen)
  let coredisc = foreign "coredisc" (gen @-> returning gen)
  let coredisc0 = foreign "coredisc0" (gen @-> long @-> returning gen)
  let coredisc2 = foreign "coredisc2" (gen @-> returning gen)

  let corediscs =
    foreign "corediscs" (long @-> ptr pari_ulong @-> returning long)

  let divisors = foreign "divisors" (gen @-> returning gen)
  let divisors_factored = foreign "divisors_factored" (gen @-> returning gen)
  let divisors0 = foreign "divisors0" (gen @-> long @-> returning gen)
  let divisorsu = foreign "divisorsu" (pari_ulong @-> returning gen)
  let divisorsu_moebius = foreign "divisorsu_moebius" (gen @-> returning gen)
  let divisorsu_fact = foreign "divisorsu_fact" (gen @-> returning gen)

  let divisorsu_fact_factored =
    foreign "divisorsu_fact_factored" (gen @-> returning gen)

  let eulerphi = foreign "eulerphi" (gen @-> returning gen)
  let eulerphiu = foreign "eulerphiu" (pari_ulong @-> returning pari_ulong)
  let eulerphiu_fact = foreign "eulerphiu_fact" (gen @-> returning pari_ulong)
  let factor_pn_1 = foreign "factor_pn_1" (gen @-> pari_ulong @-> returning gen)

  let factor_pn_1_limit =
    foreign "factor_pn_1_limit" (gen @-> long @-> pari_ulong @-> returning gen)

  let factoru_pow = foreign "factoru_pow" (pari_ulong @-> returning gen)
  let fuse_z_factor = foreign "fuse_Z_factor" (gen @-> gen @-> returning gen)
  let is_z_factor = foreign "is_Z_factor" (gen @-> returning int)
  let is_z_factornon0 = foreign "is_Z_factornon0" (gen @-> returning int)
  let is_z_factorpos = foreign "is_Z_factorpos" (gen @-> returning int)
  let is_nf_factor = foreign "is_nf_factor" (gen @-> returning int)
  let is_nf_extfactor = foreign "is_nf_extfactor" (gen @-> returning int)
  let issquarefree = foreign "issquarefree" (gen @-> returning long)
  let numdiv = foreign "numdiv" (gen @-> returning gen)
  let numdivu = foreign "numdivu" (long @-> returning long)
  let numdivu_fact = foreign "numdivu_fact" (gen @-> returning long)
  let omega = foreign "omega" (gen @-> returning long)
  let omegau = foreign "omegau" (pari_ulong @-> returning long)
  let sumdiv = foreign "sumdiv" (gen @-> returning gen)
  let sumdivk = foreign "sumdivk" (gen @-> long @-> returning gen)
  let uissquarefree = foreign "uissquarefree" (pari_ulong @-> returning long)
  let uissquarefree_fact = foreign "uissquarefree_fact" (gen @-> returning long)
  let usumdiv_fact = foreign "usumdiv_fact" (gen @-> returning gen)

  let usumdivk_fact =
    foreign "usumdivk_fact" (gen @-> pari_ulong @-> returning gen)

  let fpx_fpc_nfpoleval =
    foreign "FpX_FpC_nfpoleval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let embed_t2 = foreign "embed_T2" (gen @-> long @-> returning gen)
  let embednorm_t2 = foreign "embednorm_T2" (gen @-> long @-> returning gen)
  let embed_norm = foreign "embed_norm" (gen @-> long @-> returning gen)
  let check_zkmodule_i = foreign "check_ZKmodule_i" (gen @-> returning int)

  let check_zkmodule =
    foreign "check_ZKmodule" (gen @-> string @-> returning void)

  let checkbid = foreign "checkbid" (gen @-> returning void)
  let checkbid_i = foreign "checkbid_i" (gen @-> returning gen)
  let checkbnf = foreign "checkbnf" (gen @-> returning gen)
  let checkbnf_i = foreign "checkbnf_i" (gen @-> returning gen)
  let checkbnr = foreign "checkbnr" (gen @-> returning void)
  let checkbnr_i = foreign "checkbnr_i" (gen @-> returning gen)
  let checkabgrp = foreign "checkabgrp" (gen @-> returning void)
  let checksqmat = foreign "checksqmat" (gen @-> long @-> returning void)
  let checknf = foreign "checknf" (gen @-> returning gen)
  let checknf_i = foreign "checknf_i" (gen @-> returning gen)

  let checknfelt_mod =
    foreign "checknfelt_mod" (gen @-> gen @-> string @-> returning gen)

  let checkprid = foreign "checkprid" (gen @-> returning void)
  let checkprid_i = foreign "checkprid_i" (gen @-> returning int)
  let checkrnf = foreign "checkrnf" (gen @-> returning void)
  let checkrnf_i = foreign "checkrnf_i" (gen @-> returning int)
end

module F21 (F : Ctypes.FOREIGN) = struct
  open F

  let factoredpolred = foreign "factoredpolred" (gen @-> gen @-> returning gen)
  let factoredpolred2 = foreign "factoredpolred2" (gen @-> gen @-> returning gen)
  let galoisapply = foreign "galoisapply" (gen @-> gen @-> gen @-> returning gen)
  let get_bnf = foreign "get_bnf" (gen @-> ptr long @-> returning gen)

  let get_bnfpol =
    foreign "get_bnfpol" (gen @-> ptr gen @-> ptr gen @-> returning gen)

  let get_nf = foreign "get_nf" (gen @-> ptr long @-> returning gen)
  let get_nfpol = foreign "get_nfpol" (gen @-> ptr gen @-> returning gen)
  let get_prid = foreign "get_prid" (gen @-> returning gen)

  let idealfrobenius =
    foreign "idealfrobenius" (gen @-> gen @-> gen @-> returning gen)

  let idealfrobenius_aut =
    foreign "idealfrobenius_aut" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let idealramfrobenius =
    foreign "idealramfrobenius" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let idealramfrobenius_aut =
    foreign "idealramfrobenius_aut"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let idealramgroups =
    foreign "idealramgroups" (gen @-> gen @-> gen @-> returning gen)

  let idealramgroups_aut =
    foreign "idealramgroups_aut" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let nf_get_allroots = foreign "nf_get_allroots" (gen @-> returning gen)
  let nf_get_prec = foreign "nf_get_prec" (gen @-> returning long)

  let nfmaxord_to_nf =
    foreign "nfmaxord_to_nf" (ptr nfmaxord_t @-> gen @-> long @-> returning gen)

  let nfcertify = foreign "nfcertify" (gen @-> returning gen)
  let nfgaloismatrix = foreign "nfgaloismatrix" (gen @-> gen @-> returning gen)

  let nfgaloismatrixapply =
    foreign "nfgaloismatrixapply" (gen @-> gen @-> gen @-> returning gen)

  let nfgaloispermtobasis =
    foreign "nfgaloispermtobasis" (gen @-> gen @-> returning gen)

  let nfinit_basic =
    foreign "nfinit_basic" (ptr nfmaxord_t @-> gen @-> returning void)

  let nfinit_complete =
    foreign "nfinit_complete"
      (ptr nfmaxord_t @-> long @-> long @-> returning gen)

  let nfinit = foreign "nfinit" (gen @-> long @-> returning gen)
  let nfinit0 = foreign "nfinit0" (gen @-> long @-> long @-> returning gen)
  let nfinitred = foreign "nfinitred" (gen @-> long @-> returning gen)
  let nfinitred2 = foreign "nfinitred2" (gen @-> long @-> returning gen)
  let nfisincl = foreign "nfisincl" (gen @-> gen @-> returning gen)
  let nfisincl0 = foreign "nfisincl0" (gen @-> gen @-> long @-> returning gen)
  let nfisisom = foreign "nfisisom" (gen @-> gen @-> returning gen)
  let nfnewprec = foreign "nfnewprec" (gen @-> long @-> returning gen)

  let nfnewprec_shallow =
    foreign "nfnewprec_shallow" (gen @-> long @-> returning gen)

  let nfpoleval = foreign "nfpoleval" (gen @-> gen @-> gen @-> returning gen)
  let nfsplitting = foreign "nfsplitting" (gen @-> gen @-> returning gen)

  let nfsplitting0 =
    foreign "nfsplitting0" (gen @-> gen @-> long @-> returning gen)

  let nftyp = foreign "nftyp" (gen @-> returning long)
  let polredord = foreign "polredord" (gen @-> returning gen)
  let polred = foreign "polred" (gen @-> returning gen)
  let polred0 = foreign "polred0" (gen @-> long @-> gen @-> returning gen)
  let polred2 = foreign "polred2" (gen @-> returning gen)
  let polredabs = foreign "polredabs" (gen @-> returning gen)
  let polredabs0 = foreign "polredabs0" (gen @-> long @-> returning gen)
  let polredabs2 = foreign "polredabs2" (gen @-> returning gen)
  let polredabsall = foreign "polredabsall" (gen @-> long @-> returning gen)
  let polredbest = foreign "polredbest" (gen @-> long @-> returning gen)
  let poltomonic = foreign "poltomonic" (gen @-> ptr gen @-> returning gen)

  let rnfpolredabs =
    foreign "rnfpolredabs" (gen @-> gen @-> long @-> returning gen)

  let rnfpolredbest =
    foreign "rnfpolredbest" (gen @-> gen @-> long @-> returning gen)

  let smallpolred = foreign "smallpolred" (gen @-> returning gen)
  let smallpolred2 = foreign "smallpolred2" (gen @-> returning gen)
  let tschirnhaus = foreign "tschirnhaus" (gen @-> returning gen)
  let zx_q_mul = foreign "ZX_Q_mul" (gen @-> gen @-> returning gen)

  let zx_q_normalize =
    foreign "ZX_Q_normalize" (gen @-> ptr gen @-> returning gen)

  let zx_z_normalize =
    foreign "ZX_Z_normalize" (gen @-> ptr gen @-> returning gen)

  let zx_to_monic = foreign "ZX_to_monic" (gen @-> ptr gen @-> returning gen)

  let zx_primitive_to_monic =
    foreign "ZX_primitive_to_monic" (gen @-> ptr gen @-> returning gen)

  let zxx_q_mul = foreign "ZXX_Q_mul" (gen @-> gen @-> returning gen)
  let fq_to_nf = foreign "Fq_to_nf" (gen @-> gen @-> returning gen)
  let fqm_to_nfm = foreign "FqM_to_nfM" (gen @-> gen @-> returning gen)
  let fqv_to_nfv = foreign "FqV_to_nfV" (gen @-> gen @-> returning gen)
  let fqx_to_nfx = foreign "FqX_to_nfX" (gen @-> gen @-> returning gen)

  let rg_nffix =
    foreign "Rg_nffix" (string @-> gen @-> gen @-> int @-> returning gen)

  let rgv_nffix =
    foreign "RgV_nffix" (string @-> gen @-> gen @-> int @-> returning gen)

  let rgx_nffix =
    foreign "RgX_nffix" (string @-> gen @-> gen @-> int @-> returning gen)

  let zx_compositum_disjoint =
    foreign "ZX_compositum_disjoint" (gen @-> gen @-> returning gen)

  let zx_compositum =
    foreign "ZX_compositum" (gen @-> gen @-> ptr long @-> returning gen)

  let zpx_disc_val = foreign "ZpX_disc_val" (gen @-> gen @-> returning long)
  let zpx_gcd = foreign "ZpX_gcd" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let zpx_monic_factor =
    foreign "ZpX_monic_factor" (gen @-> gen @-> long @-> returning gen)

  let zpx_primedec = foreign "ZpX_primedec" (gen @-> gen @-> returning gen)

  let zpx_reduced_resultant =
    foreign "ZpX_reduced_resultant"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let zpx_reduced_resultant_fast =
    foreign "ZpX_reduced_resultant_fast"
      (gen @-> gen @-> gen @-> long @-> returning gen)

  let zpx_resultant_val =
    foreign "ZpX_resultant_val" (gen @-> gen @-> gen @-> long @-> returning long)

  let checkmodpr = foreign "checkmodpr" (gen @-> returning void)
  let compositum = foreign "compositum" (gen @-> gen @-> returning gen)
  let compositum2 = foreign "compositum2" (gen @-> gen @-> returning gen)
  let nfdisc = foreign "nfdisc" (gen @-> returning gen)
  let get_modpr = foreign "get_modpr" (gen @-> returning gen)
  let indexpartial = foreign "indexpartial" (gen @-> gen @-> returning gen)
  let modpr_genfq = foreign "modpr_genFq" (gen @-> returning gen)

  let nf_to_fq_init =
    foreign "nf_to_Fq_init"
      (gen @-> ptr gen @-> ptr gen @-> ptr gen @-> returning gen)

  let nf_to_fq = foreign "nf_to_Fq" (gen @-> gen @-> gen @-> returning gen)
  let nfm_to_fqm = foreign "nfM_to_FqM" (gen @-> gen @-> gen @-> returning gen)
  let nfv_to_fqv = foreign "nfV_to_FqV" (gen @-> gen @-> gen @-> returning gen)
  let nfx_to_fqx = foreign "nfX_to_FqX" (gen @-> gen @-> gen @-> returning gen)

  let nfx_to_monic =
    foreign "nfX_to_monic" (gen @-> gen @-> ptr gen @-> returning gen)

  let nfbasis = foreign "nfbasis" (gen @-> ptr gen @-> returning gen)

  let nfcompositum =
    foreign "nfcompositum" (gen @-> gen @-> gen @-> long @-> returning gen)

  let nfdiscfactors = foreign "nfdiscfactors" (gen @-> returning gen)

  let nfmaxord =
    foreign "nfmaxord" (ptr nfmaxord_t @-> gen @-> long @-> returning void)

  let nfmodpr = foreign "nfmodpr" (gen @-> gen @-> gen @-> returning gen)
  let nfmodprinit = foreign "nfmodprinit" (gen @-> gen @-> returning gen)

  let nfmodprinit0 =
    foreign "nfmodprinit0" (gen @-> gen @-> long @-> returning gen)

  let nfmodprlift = foreign "nfmodprlift" (gen @-> gen @-> gen @-> returning gen)

  let nfreducemodpr =
    foreign "nfreducemodpr" (gen @-> gen @-> gen @-> returning gen)

  let polcompositum0 =
    foreign "polcompositum0" (gen @-> gen @-> long @-> returning gen)

  let idealprimedec = foreign "idealprimedec" (gen @-> gen @-> returning gen)

  let idealprimedec_galois =
    foreign "idealprimedec_galois" (gen @-> gen @-> returning gen)

  let idealprimedec_degrees =
    foreign "idealprimedec_degrees" (gen @-> gen @-> returning gen)

  let idealprimedec_kummer =
    foreign "idealprimedec_kummer"
      (gen @-> gen @-> long @-> gen @-> returning gen)
end

module F22 (F : Ctypes.FOREIGN) = struct
  open F

  let idealprimedec_limit_f =
    foreign "idealprimedec_limit_f" (gen @-> gen @-> long @-> returning gen)

  let idealprimedec_limit_norm =
    foreign "idealprimedec_limit_norm" (gen @-> gen @-> gen @-> returning gen)

  let poldiscfactors = foreign "poldiscfactors" (gen @-> long @-> returning gen)
  let rnfbasis = foreign "rnfbasis" (gen @-> gen @-> returning gen)

  let rnfdedekind =
    foreign "rnfdedekind" (gen @-> gen @-> gen @-> long @-> returning gen)

  let rnfdet = foreign "rnfdet" (gen @-> gen @-> returning gen)

  let rnfdisc_factored =
    foreign "rnfdisc_factored" (gen @-> gen @-> ptr gen @-> returning gen)

  let rnfdiscf = foreign "rnfdiscf" (gen @-> gen @-> returning gen)
  let rnfequation = foreign "rnfequation" (gen @-> gen @-> returning gen)

  let rnfequation0 =
    foreign "rnfequation0" (gen @-> gen @-> long @-> returning gen)

  let rnfequation2 = foreign "rnfequation2" (gen @-> gen @-> returning gen)
  let nf_pv_to_prv = foreign "nf_pV_to_prV" (gen @-> gen @-> returning gen)
  let nf_rnfeq = foreign "nf_rnfeq" (gen @-> gen @-> returning gen)
  let nf_rnfeqsimple = foreign "nf_rnfeqsimple" (gen @-> gen @-> returning gen)

  let rnfequationall =
    foreign "rnfequationall"
      (gen @-> gen @-> ptr long @-> ptr gen @-> returning gen)

  let rnfhnfbasis = foreign "rnfhnfbasis" (gen @-> gen @-> returning gen)
  let rnfisfree = foreign "rnfisfree" (gen @-> gen @-> returning long)

  let rnflllgram =
    foreign "rnflllgram" (gen @-> gen @-> gen @-> long @-> returning gen)

  let rnfpolred = foreign "rnfpolred" (gen @-> gen @-> long @-> returning gen)
  let rnfpseudobasis = foreign "rnfpseudobasis" (gen @-> gen @-> returning gen)

  let rnfsimplifybasis =
    foreign "rnfsimplifybasis" (gen @-> gen @-> returning gen)

  let rnfsteinitz = foreign "rnfsteinitz" (gen @-> gen @-> returning gen)

  let factorial_lval =
    foreign "factorial_lval" (pari_ulong @-> pari_ulong @-> returning long)

  let zk_to_fq_init =
    foreign "zk_to_Fq_init"
      (gen @-> ptr gen @-> ptr gen @-> ptr gen @-> returning gen)

  let zk_to_fq = foreign "zk_to_Fq" (gen @-> gen @-> returning gen)
  let qxqv_to_fpm = foreign "QXQV_to_FpM" (gen @-> gen @-> gen @-> returning gen)
  let zkmodprinit = foreign "zkmodprinit" (gen @-> gen @-> returning gen)
  let idealstar = foreign "Idealstar" (gen @-> gen @-> long @-> returning gen)

  let idealstarmod =
    foreign "Idealstarmod" (gen @-> gen @-> long @-> gen @-> returning gen)

  let idealstarprk =
    foreign "Idealstarprk" (gen @-> gen @-> long @-> long @-> returning gen)

  let rgc_to_nfc = foreign "RgC_to_nfC" (gen @-> gen @-> returning gen)
  let rgm_rgx_mul = foreign "RgM_RgX_mul" (gen @-> gen @-> returning gen)
  let rgm_to_nfm = foreign "RgM_to_nfM" (gen @-> gen @-> returning gen)
  let rgx_to_nfx = foreign "RgX_to_nfX" (gen @-> gen @-> returning gen)
  let zc_nfval = foreign "ZC_nfval" (gen @-> gen @-> returning long)

  let zc_nfvalrem =
    foreign "ZC_nfvalrem" (gen @-> gen @-> ptr gen @-> returning long)

  let zc_prdvd = foreign "ZC_prdvd" (gen @-> gen @-> returning int)
  let zm_zx_mul = foreign "ZM_ZX_mul" (gen @-> gen @-> returning gen)
  let zv_snf_gcd = foreign "ZV_snf_gcd" (gen @-> gen @-> returning gen)
  let algtobasis = foreign "algtobasis" (gen @-> gen @-> returning gen)
  let basistoalg = foreign "basistoalg" (gen @-> gen @-> returning gen)
  let ei_multable = foreign "ei_multable" (gen @-> long @-> returning gen)

  let get_nf_field =
    foreign "get_nf_field" (ptr (ptr void) @-> gen @-> returning (ptr bb_field))

  let famat_nfvalrem =
    foreign "famat_nfvalrem" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let gpnfvalrem =
    foreign "gpnfvalrem" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let idealfactorback =
    foreign "idealfactorback" (gen @-> gen @-> gen @-> int @-> returning gen)

  let ideallist = foreign "ideallist" (gen @-> long @-> returning gen)
  let ideallist0 = foreign "ideallist0" (gen @-> long @-> long @-> returning gen)
  let gideallist = foreign "gideallist" (gen @-> gen @-> long @-> returning gen)

  let ideallistarch =
    foreign "ideallistarch" (gen @-> gen @-> gen @-> returning gen)

  let ideallog = foreign "ideallog" (gen @-> gen @-> gen @-> returning gen)

  let ideallogmod =
    foreign "ideallogmod" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let ideallog_units = foreign "ideallog_units" (gen @-> gen @-> returning gen)

  let ideallog_units0 =
    foreign "ideallog_units0" (gen @-> gen @-> gen @-> returning gen)

  let idealprincipalunits =
    foreign "idealprincipalunits" (gen @-> gen @-> long @-> returning gen)

  let idealstar0 = foreign "idealstar0" (gen @-> gen @-> long @-> returning gen)

  let idealstarmod =
    foreign "idealstarmod" (gen @-> gen @-> long @-> gen @-> returning gen)

  let indices_to_vec01 =
    foreign "indices_to_vec01" (gen @-> long @-> returning gen)

  let matalgtobasis = foreign "matalgtobasis" (gen @-> gen @-> returning gen)
  let matbasistoalg = foreign "matbasistoalg" (gen @-> gen @-> returning gen)
  let multable = foreign "multable" (gen @-> gen @-> returning gen)

  let nf_to_scalar_or_alg =
    foreign "nf_to_scalar_or_alg" (gen @-> gen @-> returning gen)

  let nf_to_scalar_or_basis =
    foreign "nf_to_scalar_or_basis" (gen @-> gen @-> returning gen)

  let nf_cxlog = foreign "nf_cxlog" (gen @-> gen @-> long @-> returning gen)
  let nfv_cxlog = foreign "nfV_cxlog" (gen @-> gen @-> long @-> returning gen)
  let nfadd = foreign "nfadd" (gen @-> gen @-> gen @-> returning gen)

  let nfchecksigns =
    foreign "nfchecksigns" (gen @-> gen @-> gen @-> returning int)

  let nfdiv = foreign "nfdiv" (gen @-> gen @-> gen @-> returning gen)
  let nfdiveuc = foreign "nfdiveuc" (gen @-> gen @-> gen @-> returning gen)
  let nfdivrem = foreign "nfdivrem" (gen @-> gen @-> gen @-> returning gen)
  let nfembed = foreign "nfembed" (gen @-> gen @-> long @-> returning gen)

  let nfeltembed =
    foreign "nfeltembed" (gen @-> gen @-> gen @-> long @-> returning gen)

  let nfeltembed_i =
    foreign "nfeltembed_i" (ptr gen @-> gen @-> gen @-> long @-> returning gen)

  let nfeltsign = foreign "nfeltsign" (gen @-> gen @-> gen @-> returning gen)

  let nffactorback =
    foreign "nffactorback" (gen @-> gen @-> gen @-> returning gen)

  let nfinv = foreign "nfinv" (gen @-> gen @-> returning gen)

  let nfinvmodideal =
    foreign "nfinvmodideal" (gen @-> gen @-> gen @-> returning gen)

  let nfissquare =
    foreign "nfissquare" (gen @-> gen @-> ptr gen @-> returning long)

  let nfispower =
    foreign "nfispower" (gen @-> gen @-> long @-> ptr gen @-> returning long)

  let nflogembed =
    foreign "nflogembed" (gen @-> gen @-> ptr gen @-> long @-> returning gen)

  let nfm_det = foreign "nfM_det" (gen @-> gen @-> returning gen)
  let nfm_inv = foreign "nfM_inv" (gen @-> gen @-> returning gen)
  let nfm_ker = foreign "nfM_ker" (gen @-> gen @-> returning gen)
  let nfm_mul = foreign "nfM_mul" (gen @-> gen @-> gen @-> returning gen)
  let nfm_nfc_mul = foreign "nfM_nfC_mul" (gen @-> gen @-> gen @-> returning gen)
  let nfmod = foreign "nfmod" (gen @-> gen @-> gen @-> returning gen)
  let nfmul = foreign "nfmul" (gen @-> gen @-> gen @-> returning gen)
  let nfmuli = foreign "nfmuli" (gen @-> gen @-> gen @-> returning gen)
  let nfnorm = foreign "nfnorm" (gen @-> gen @-> returning gen)
  let nfpolsturm = foreign "nfpolsturm" (gen @-> gen @-> gen @-> returning gen)
  let nfpow = foreign "nfpow" (gen @-> gen @-> gen @-> returning gen)
  let nfpow_u = foreign "nfpow_u" (gen @-> gen @-> pari_ulong @-> returning gen)

  let nfpowmodideal =
    foreign "nfpowmodideal" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let nfsign = foreign "nfsign" (gen @-> gen @-> returning gen)
  let nfsign_arch = foreign "nfsign_arch" (gen @-> gen @-> gen @-> returning gen)

  let nfsign_from_logarch =
    foreign "nfsign_from_logarch" (gen @-> gen @-> gen @-> returning gen)

  let nfsqr = foreign "nfsqr" (gen @-> gen @-> returning gen)
  let nfsqri = foreign "nfsqri" (gen @-> gen @-> returning gen)
  let nfsub = foreign "nfsub" (gen @-> gen @-> gen @-> returning gen)
  let nftrace = foreign "nftrace" (gen @-> gen @-> returning gen)
end

module F23 (F : Ctypes.FOREIGN) = struct
  open F

  let nfval = foreign "nfval" (gen @-> gen @-> gen @-> returning long)

  let nfvalrem =
    foreign "nfvalrem" (gen @-> gen @-> gen @-> ptr gen @-> returning long)

  let polmod_nffix =
    foreign "polmod_nffix" (string @-> gen @-> gen @-> int @-> returning gen)

  let polmod_nffix2 =
    foreign "polmod_nffix2"
      (string @-> gen @-> gen @-> gen @-> int @-> returning gen)

  let pr_basis_perm = foreign "pr_basis_perm" (gen @-> gen @-> returning gen)
  let pr_equal = foreign "pr_equal" (gen @-> gen @-> returning int)
  let rnfalgtobasis = foreign "rnfalgtobasis" (gen @-> gen @-> returning gen)
  let rnfbasistoalg = foreign "rnfbasistoalg" (gen @-> gen @-> returning gen)
  let rnfeltnorm = foreign "rnfeltnorm" (gen @-> gen @-> returning gen)
  let rnfelttrace = foreign "rnfelttrace" (gen @-> gen @-> returning gen)

  let set_sign_mod_divisor =
    foreign "set_sign_mod_divisor"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let tablemul = foreign "tablemul" (gen @-> gen @-> gen @-> returning gen)

  let tablemul_ei =
    foreign "tablemul_ei" (gen @-> gen @-> long @-> returning gen)

  let tablemul_ei_ej =
    foreign "tablemul_ei_ej" (gen @-> long @-> long @-> returning gen)

  let tablemulvec = foreign "tablemulvec" (gen @-> gen @-> gen @-> returning gen)
  let tablesqr = foreign "tablesqr" (gen @-> gen @-> returning gen)
  let vec01_to_indices = foreign "vec01_to_indices" (gen @-> returning gen)

  let vecsmall01_to_indices =
    foreign "vecsmall01_to_indices" (gen @-> returning gen)

  let zk_inv = foreign "zk_inv" (gen @-> gen @-> returning gen)
  let zk_multable = foreign "zk_multable" (gen @-> gen @-> returning gen)

  let zk_scalar_or_multable =
    foreign "zk_scalar_or_multable" (gen @-> gen @-> returning gen)

  let zkchinese = foreign "zkchinese" (gen @-> gen @-> gen @-> returning gen)
  let zkchinese1 = foreign "zkchinese1" (gen @-> gen @-> returning gen)

  let zkchineseinit =
    foreign "zkchineseinit" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let zkmultable_capz = foreign "zkmultable_capZ" (gen @-> returning gen)
  let zkmultable_inv = foreign "zkmultable_inv" (gen @-> returning gen)

  let fl_invgen =
    foreign "Fl_invgen"
      (pari_ulong @-> pari_ulong @-> ptr pari_ulong @-> returning pari_ulong)

  let rm_round_maxrank = foreign "RM_round_maxrank" (gen @-> returning gen)
  let zm_famat_limit = foreign "ZM_famat_limit" (gen @-> gen @-> returning gen)
  let zv_cba = foreign "ZV_cba" (gen @-> returning gen)
  let zv_cba_extend = foreign "ZV_cba_extend" (gen @-> gen @-> returning gen)
  let z_cba = foreign "Z_cba" (gen @-> gen @-> returning gen)
  let z_ppgle = foreign "Z_ppgle" (gen @-> gen @-> returning gen)
  let z_ppio = foreign "Z_ppio" (gen @-> gen @-> returning gen)
  let z_ppo = foreign "Z_ppo" (gen @-> gen @-> returning gen)

  let famatv_factorback =
    foreign "famatV_factorback" (gen @-> gen @-> returning gen)

  let famatv_zv_factorback =
    foreign "famatV_zv_factorback" (gen @-> gen @-> returning gen)

  let famat_z_gcd = foreign "famat_Z_gcd" (gen @-> gen @-> returning gen)
  let famat_div = foreign "famat_div" (gen @-> gen @-> returning gen)

  let famat_div_shallow =
    foreign "famat_div_shallow" (gen @-> gen @-> returning gen)

  let famat_idealfactor =
    foreign "famat_idealfactor" (gen @-> gen @-> returning gen)

  let famat_inv = foreign "famat_inv" (gen @-> returning gen)
  let famat_inv_shallow = foreign "famat_inv_shallow" (gen @-> returning gen)

  let famat_makecoprime =
    foreign "famat_makecoprime"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let famat_mul = foreign "famat_mul" (gen @-> gen @-> returning gen)

  let famat_mul_shallow =
    foreign "famat_mul_shallow" (gen @-> gen @-> returning gen)

  let famat_mulpow_shallow =
    foreign "famat_mulpow_shallow" (gen @-> gen @-> gen @-> returning gen)

  let famat_mulpows_shallow =
    foreign "famat_mulpows_shallow" (gen @-> gen @-> long @-> returning gen)

  let famat_pow = foreign "famat_pow" (gen @-> gen @-> returning gen)

  let famat_pow_shallow =
    foreign "famat_pow_shallow" (gen @-> gen @-> returning gen)

  let famat_pows_shallow =
    foreign "famat_pows_shallow" (gen @-> long @-> returning gen)

  let famat_reduce = foreign "famat_reduce" (gen @-> returning gen)

  let famat_remove_trivial =
    foreign "famat_remove_trivial" (gen @-> returning gen)

  let famat_sqr = foreign "famat_sqr" (gen @-> returning gen)
  let famat_to_nf = foreign "famat_to_nf" (gen @-> gen @-> returning gen)

  let famat_to_nf_moddivisor =
    foreign "famat_to_nf_moddivisor"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let famat_to_nf_modideal_coprime =
    foreign "famat_to_nf_modideal_coprime"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let famatsmall_reduce = foreign "famatsmall_reduce" (gen @-> returning gen)

  let gpidealfactor =
    foreign "gpidealfactor" (gen @-> gen @-> gen @-> returning gen)

  let gpidealval = foreign "gpidealval" (gen @-> gen @-> gen @-> returning gen)

  let idealhnf_z_factor =
    foreign "idealHNF_Z_factor" (gen @-> ptr gen @-> ptr gen @-> returning gen)

  let idealhnf_z_factor_i =
    foreign "idealHNF_Z_factor_i"
      (gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let idealhnf_inv = foreign "idealHNF_inv" (gen @-> gen @-> returning gen)
  let idealhnf_inv_z = foreign "idealHNF_inv_Z" (gen @-> gen @-> returning gen)

  let idealhnf_mul =
    foreign "idealHNF_mul" (gen @-> gen @-> gen @-> returning gen)

  let idealadd = foreign "idealadd" (gen @-> gen @-> gen @-> returning gen)

  let idealaddmultoone =
    foreign "idealaddmultoone" (gen @-> gen @-> returning gen)

  let idealaddtoone =
    foreign "idealaddtoone" (gen @-> gen @-> gen @-> returning gen)

  let idealaddtoone0 =
    foreign "idealaddtoone0" (gen @-> gen @-> gen @-> returning gen)

  let idealaddtoone_i =
    foreign "idealaddtoone_i" (gen @-> gen @-> gen @-> returning gen)

  let idealaddtoone_raw =
    foreign "idealaddtoone_raw" (gen @-> gen @-> gen @-> returning gen)

  let idealappr = foreign "idealappr" (gen @-> gen @-> returning gen)
  let idealappr0 = foreign "idealappr0" (gen @-> gen @-> long @-> returning gen)
  let idealapprfact = foreign "idealapprfact" (gen @-> gen @-> returning gen)

  let idealchinese =
    foreign "idealchinese" (gen @-> gen @-> gen @-> returning gen)

  let idealcoprime =
    foreign "idealcoprime" (gen @-> gen @-> gen @-> returning gen)

  let idealcoprimefact =
    foreign "idealcoprimefact" (gen @-> gen @-> gen @-> returning gen)

  let idealdiv = foreign "idealdiv" (gen @-> gen @-> gen @-> returning gen)

  let idealdiv0 =
    foreign "idealdiv0" (gen @-> gen @-> gen @-> long @-> returning gen)

  let idealdivexact =
    foreign "idealdivexact" (gen @-> gen @-> gen @-> returning gen)

  let idealdivpowprime =
    foreign "idealdivpowprime" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let idealdown = foreign "idealdown" (gen @-> gen @-> returning gen)
  let idealfactor = foreign "idealfactor" (gen @-> gen @-> returning gen)

  let idealfactor_limit =
    foreign "idealfactor_limit" (gen @-> gen @-> pari_ulong @-> returning gen)

  let idealfactor_partial =
    foreign "idealfactor_partial" (gen @-> gen @-> gen @-> returning gen)

  let idealhnf = foreign "idealhnf" (gen @-> gen @-> returning gen)
  let idealhnf0 = foreign "idealhnf0" (gen @-> gen @-> gen @-> returning gen)

  let idealhnf_principal =
    foreign "idealhnf_principal" (gen @-> gen @-> returning gen)

  let idealhnf_shallow =
    foreign "idealhnf_shallow" (gen @-> gen @-> returning gen)

  let idealhnf_two = foreign "idealhnf_two" (gen @-> gen @-> returning gen)

  let idealintersect =
    foreign "idealintersect" (gen @-> gen @-> gen @-> returning gen)

  let idealinv = foreign "idealinv" (gen @-> gen @-> returning gen)
  let idealismaximal = foreign "idealismaximal" (gen @-> gen @-> returning gen)

  let idealispower =
    foreign "idealispower" (gen @-> gen @-> long @-> ptr gen @-> returning long)

  let idealmin = foreign "idealmin" (gen @-> gen @-> gen @-> returning gen)
  let idealmul = foreign "idealmul" (gen @-> gen @-> gen @-> returning gen)

  let idealmul0 =
    foreign "idealmul0" (gen @-> gen @-> gen @-> long @-> returning gen)

  let idealmulpowprime =
    foreign "idealmulpowprime" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let idealmulred = foreign "idealmulred" (gen @-> gen @-> gen @-> returning gen)
  let idealnorm = foreign "idealnorm" (gen @-> gen @-> returning gen)
end

module F24 (F : Ctypes.FOREIGN) = struct
  open F

  let idealnumden = foreign "idealnumden" (gen @-> gen @-> returning gen)
  let idealpow = foreign "idealpow" (gen @-> gen @-> gen @-> returning gen)

  let idealpow0 =
    foreign "idealpow0" (gen @-> gen @-> gen @-> long @-> returning gen)

  let idealpowred = foreign "idealpowred" (gen @-> gen @-> gen @-> returning gen)
  let idealpows = foreign "idealpows" (gen @-> gen @-> long @-> returning gen)
  let idealprod = foreign "idealprod" (gen @-> gen @-> returning gen)
  let idealprodprime = foreign "idealprodprime" (gen @-> gen @-> returning gen)

  let idealprodval =
    foreign "idealprodval" (gen @-> gen @-> gen @-> returning long)

  let idealpseudomin = foreign "idealpseudomin" (gen @-> gen @-> returning gen)

  let idealpseudomin_nonscalar =
    foreign "idealpseudomin_nonscalar" (gen @-> gen @-> returning gen)

  let idealpseudominvec =
    foreign "idealpseudominvec" (gen @-> gen @-> returning gen)

  let idealpseudored = foreign "idealpseudored" (gen @-> gen @-> returning gen)
  let idealred0 = foreign "idealred0" (gen @-> gen @-> gen @-> returning gen)
  let idealred_elt = foreign "idealred_elt" (gen @-> gen @-> returning gen)

  let idealredmodpower =
    foreign "idealredmodpower"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let idealsqr = foreign "idealsqr" (gen @-> gen @-> returning gen)
  let idealtwoelt = foreign "idealtwoelt" (gen @-> gen @-> returning gen)

  let idealtwoelt0 =
    foreign "idealtwoelt0" (gen @-> gen @-> gen @-> returning gen)

  let idealtwoelt2 =
    foreign "idealtwoelt2" (gen @-> gen @-> gen @-> returning gen)

  let idealtyp = foreign "idealtyp" (ptr gen @-> ptr gen @-> returning long)
  let idealval = foreign "idealval" (gen @-> gen @-> gen @-> returning long)
  let isideal = foreign "isideal" (gen @-> gen @-> returning long)
  let matreduce = foreign "matreduce" (gen @-> returning gen)

  let nfc_multable_mul =
    foreign "nfC_multable_mul" (gen @-> gen @-> returning gen)

  let nfc_nf_mul = foreign "nfC_nf_mul" (gen @-> gen @-> gen @-> returning gen)
  let nf_get_gtwist = foreign "nf_get_Gtwist" (gen @-> gen @-> returning gen)
  let nf_get_gtwist1 = foreign "nf_get_Gtwist1" (gen @-> long @-> returning gen)

  let nf_to_fp_coprime =
    foreign "nf_to_Fp_coprime" (gen @-> gen @-> gen @-> returning gen)

  let nfdetint = foreign "nfdetint" (gen @-> gen @-> returning gen)

  let nfdivmodpr =
    foreign "nfdivmodpr" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let nfhnf = foreign "nfhnf" (gen @-> gen @-> returning gen)
  let nfhnf0 = foreign "nfhnf0" (gen @-> gen @-> long @-> returning gen)
  let nfhnfmod = foreign "nfhnfmod" (gen @-> gen @-> gen @-> returning gen)
  let nfkermodpr = foreign "nfkermodpr" (gen @-> gen @-> gen @-> returning gen)

  let nfmulmodpr =
    foreign "nfmulmodpr" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let nfpowmodpr =
    foreign "nfpowmodpr" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let nfreduce = foreign "nfreduce" (gen @-> gen @-> gen @-> returning gen)
  let nfsnf = foreign "nfsnf" (gen @-> gen @-> returning gen)
  let nfsnf0 = foreign "nfsnf0" (gen @-> gen @-> long @-> returning gen)

  let nfsolvemodpr =
    foreign "nfsolvemodpr" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let prv_lcm_capz = foreign "prV_lcm_capZ" (gen @-> returning gen)
  let prv_primes = foreign "prV_primes" (gen @-> returning gen)
  let pr_hnf = foreign "pr_hnf" (gen @-> gen @-> returning gen)
  let pr_inv = foreign "pr_inv" (gen @-> returning gen)
  let pr_inv_p = foreign "pr_inv_p" (gen @-> returning gen)
  let pr_uniformizer = foreign "pr_uniformizer" (gen @-> gen @-> returning gen)

  let sunits_makecoprime =
    foreign "sunits_makecoprime" (gen @-> gen @-> gen @-> returning gen)

  let to_famat = foreign "to_famat" (gen @-> gen @-> returning gen)

  let to_famat_shallow =
    foreign "to_famat_shallow" (gen @-> gen @-> returning gen)

  let u_ppo =
    foreign "u_ppo" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let vecdiv = foreign "vecdiv" (gen @-> gen @-> returning gen)
  let vecinv = foreign "vecinv" (gen @-> returning gen)
  let vecmul = foreign "vecmul" (gen @-> gen @-> returning gen)
  let vecpow = foreign "vecpow" (gen @-> gen @-> returning gen)
  let vecsqr = foreign "vecsqr" (gen @-> returning gen)

  let zkc_multable_mul =
    foreign "zkC_multable_mul" (gen @-> gen @-> returning gen)

  let eltreltoabs = foreign "eltreltoabs" (gen @-> gen @-> returning gen)
  let eltabstorel = foreign "eltabstorel" (gen @-> gen @-> returning gen)

  let eltabstorel_lift =
    foreign "eltabstorel_lift" (gen @-> gen @-> returning gen)

  let nf_nfzk = foreign "nf_nfzk" (gen @-> gen @-> returning gen)

  let rnf_build_nfabs =
    foreign "rnf_build_nfabs" (gen @-> long @-> returning gen)

  let rnf_zkabs = foreign "rnf_zkabs" (gen @-> returning gen)
  let nfeltup = foreign "nfeltup" (gen @-> gen @-> gen @-> returning gen)
  let rnfcomplete = foreign "rnfcomplete" (gen @-> returning void)
  let rnfeltabstorel = foreign "rnfeltabstorel" (gen @-> gen @-> returning gen)
  let rnfeltdown = foreign "rnfeltdown" (gen @-> gen @-> returning gen)

  let rnfeltdown0 =
    foreign "rnfeltdown0" (gen @-> gen @-> long @-> returning gen)

  let rnfeltreltoabs = foreign "rnfeltreltoabs" (gen @-> gen @-> returning gen)
  let rnfeltup = foreign "rnfeltup" (gen @-> gen @-> returning gen)
  let rnfeltup0 = foreign "rnfeltup0" (gen @-> gen @-> long @-> returning gen)

  let rnfidealabstorel =
    foreign "rnfidealabstorel" (gen @-> gen @-> returning gen)

  let rnfidealdown = foreign "rnfidealdown" (gen @-> gen @-> returning gen)
  let rnfidealfactor = foreign "rnfidealfactor" (gen @-> gen @-> returning gen)
  let rnfidealhnf = foreign "rnfidealhnf" (gen @-> gen @-> returning gen)
  let rnfidealmul = foreign "rnfidealmul" (gen @-> gen @-> gen @-> returning gen)
  let rnfidealnormabs = foreign "rnfidealnormabs" (gen @-> gen @-> returning gen)
  let rnfidealnormrel = foreign "rnfidealnormrel" (gen @-> gen @-> returning gen)

  let rnfidealprimedec =
    foreign "rnfidealprimedec" (gen @-> gen @-> returning gen)

  let rnfidealreltoabs =
    foreign "rnfidealreltoabs" (gen @-> gen @-> returning gen)

  let rnfidealreltoabs0 =
    foreign "rnfidealreltoabs0" (gen @-> gen @-> long @-> returning gen)

  let rnfidealtwoelement =
    foreign "rnfidealtwoelement" (gen @-> gen @-> returning gen)

  let rnfidealup = foreign "rnfidealup" (gen @-> gen @-> returning gen)

  let rnfidealup0 =
    foreign "rnfidealup0" (gen @-> gen @-> long @-> returning gen)

  let rnfinit = foreign "rnfinit" (gen @-> gen @-> returning gen)
  let rnfinit0 = foreign "rnfinit0" (gen @-> gen @-> long @-> returning gen)
  let get_arith_zzm = foreign "get_arith_ZZM" (gen @-> returning gen)
  let get_arith_z = foreign "get_arith_Z" (gen @-> returning gen)

  let gen_ph_log =
    foreign "gen_PH_log"
      (gen @-> gen @-> gen @-> ptr void @-> ptr bb_group @-> returning gen)

  let gen_shanks_init =
    foreign "gen_Shanks_init"
      (gen @-> long @-> ptr void @-> ptr bb_group @-> returning gen)

  let gen_shanks =
    foreign "gen_Shanks"
      (gen @-> gen @-> pari_ulong @-> ptr void @-> ptr bb_group
     @-> returning gen)

  let gen_shanks_sqrtn =
    foreign "gen_Shanks_sqrtn"
      (gen @-> gen @-> gen @-> ptr gen @-> ptr void @-> ptr bb_group
     @-> returning gen)

  let gen_gener =
    foreign "gen_gener" (gen @-> ptr void @-> ptr bb_group @-> returning gen)

  let gen_ellgens =
    foreign "gen_ellgens"
      (gen @-> gen @-> gen @-> ptr void @-> ptr bb_group
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> gen @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let gen_ellgroup =
    foreign "gen_ellgroup"
      (gen @-> gen @-> ptr gen @-> ptr void @-> ptr bb_group
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> gen @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let gen_factored_order =
    foreign "gen_factored_order"
      (gen @-> gen @-> ptr void @-> ptr bb_group @-> returning gen)

  let gen_order =
    foreign "gen_order"
      (gen @-> gen @-> ptr void @-> ptr bb_group @-> returning gen)

  let gen_select_order =
    foreign "gen_select_order"
      (gen @-> ptr void @-> ptr bb_group @-> returning gen)

  let gen_plog =
    foreign "gen_plog"
      (gen @-> gen @-> gen @-> ptr void @-> ptr bb_group @-> returning gen)

  let gen_pow =
    foreign "gen_pow"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let gen_pow_fold =
    foreign "gen_pow_fold"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> returning gen)
end

module F25 (F : Ctypes.FOREIGN) = struct
  open F

  let gen_pow_fold_i =
    foreign "gen_pow_fold_i"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> returning gen)

  let gen_pow_i =
    foreign "gen_pow_i"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let gen_pow_init =
    foreign "gen_pow_init"
      (gen @-> gen @-> long @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let gen_pow_table =
    foreign "gen_pow_table"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let gen_powers =
    foreign "gen_powers"
      (gen @-> long @-> int @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> returning gen)
      @-> returning gen)

  let gen_powu =
    foreign "gen_powu"
      (gen @-> pari_ulong @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let gen_powu_fold =
    foreign "gen_powu_fold"
      (gen @-> pari_ulong @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> returning gen)

  let gen_powu_fold_i =
    foreign "gen_powu_fold_i"
      (gen @-> pari_ulong @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> returning gen)

  let gen_powu_i =
    foreign "gen_powu_i"
      (gen @-> pari_ulong @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let gen_product =
    foreign "gen_product"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> returning gen)

  let matdetmod = foreign "matdetmod" (gen @-> gen @-> returning gen)

  let matimagemod =
    foreign "matimagemod" (gen @-> gen @-> ptr gen @-> returning gen)

  let matinvmod = foreign "matinvmod" (gen @-> gen @-> returning gen)
  let matkermod = foreign "matkermod" (gen @-> gen @-> ptr gen @-> returning gen)

  let matsolvemod =
    foreign "matsolvemod" (gen @-> gen @-> gen @-> long @-> returning gen)

  let bernfrac = foreign "bernfrac" (long @-> returning gen)
  let bernpol = foreign "bernpol" (long @-> long @-> returning gen)
  let bernreal = foreign "bernreal" (long @-> long @-> returning gen)
  let bernvec = foreign "bernvec" (long @-> returning gen)
  let constbern = foreign "constbern" (long @-> returning void)
  let eulerfrac = foreign "eulerfrac" (long @-> returning gen)
  let eulerpol = foreign "eulerpol" (long @-> long @-> returning gen)
  let eulerreal = foreign "eulerreal" (long @-> long @-> returning gen)
  let eulervec = foreign "eulervec" (long @-> returning gen)
  let harmonic = foreign "harmonic" (pari_ulong @-> returning gen)
  let harmonic0 = foreign "harmonic0" (pari_ulong @-> gen @-> returning gen)

  let qr_init =
    foreign "QR_init"
      (gen @-> ptr gen @-> ptr gen @-> ptr gen @-> long @-> returning int)

  let r_from_qr = foreign "R_from_QR" (gen @-> long @-> returning gen)
  let rgm_babai = foreign "RgM_Babai" (gen @-> gen @-> returning gen)

  let rgm_qr_init =
    foreign "RgM_QR_init"
      (gen @-> ptr gen @-> ptr gen @-> ptr gen @-> long @-> returning int)

  let rgm_gram_schmidt =
    foreign "RgM_gram_schmidt" (gen @-> ptr gen @-> returning gen)

  let algdep = foreign "algdep" (gen @-> long @-> returning gen)
  let algdep0 = foreign "algdep0" (gen @-> long @-> long @-> returning gen)

  let bestapprnf =
    foreign "bestapprnf" (gen @-> gen @-> gen @-> long @-> returning gen)

  let forqfvec =
    foreign "forqfvec"
      (ptr void
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> gen @-> double @-> returning long)
      @-> gen @-> gen @-> returning void)

  let forqfvec0 = foreign "forqfvec0" (gen @-> gen @-> gen @-> returning void)

  let forqfvec1 =
    foreign "forqfvec1"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> gen @-> gen @-> returning void)

  let gaussred_from_qr =
    foreign "gaussred_from_QR" (gen @-> long @-> returning gen)

  let lindep = foreign "lindep" (gen @-> returning gen)
  let lindep_xadic = foreign "lindep_Xadic" (gen @-> returning gen)
  let lindep_bit = foreign "lindep_bit" (gen @-> long @-> returning gen)
  let lindep_padic = foreign "lindep_padic" (gen @-> returning gen)
  let lindep0 = foreign "lindep0" (gen @-> long @-> returning gen)
  let lindep2 = foreign "lindep2" (gen @-> long @-> returning gen)
  let lindepfull_bit = foreign "lindepfull_bit" (gen @-> long @-> returning gen)
  let mathouseholder = foreign "mathouseholder" (gen @-> gen @-> returning gen)
  let matqr = foreign "matqr" (gen @-> long @-> long @-> returning gen)
  let minim = foreign "minim" (gen @-> gen @-> gen @-> returning gen)
  let minim_raw = foreign "minim_raw" (gen @-> gen @-> gen @-> returning gen)
  let minim_zm = foreign "minim_zm" (gen @-> gen @-> gen @-> returning gen)
  let minim2 = foreign "minim2" (gen @-> gen @-> gen @-> returning gen)

  let qfminim0 =
    foreign "qfminim0" (gen @-> gen @-> gen @-> long @-> long @-> returning gen)

  let qfperfection = foreign "qfperfection" (gen @-> returning gen)
  let qfrep0 = foreign "qfrep0" (gen @-> gen @-> long @-> returning gen)
  let seralgdep = foreign "seralgdep" (gen @-> long @-> long @-> returning gen)
  let serdiffdep = foreign "serdiffdep" (gen @-> long @-> long @-> returning gen)

  let vandermondeinverse =
    foreign "vandermondeinverse" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let vandermondeinverseinit =
    foreign "vandermondeinverseinit" (gen @-> returning gen)

  let zncoppersmith =
    foreign "zncoppersmith" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let qxq_reverse = foreign "QXQ_reverse" (gen @-> gen @-> returning gen)
  let vec_equiv = foreign "vec_equiv" (gen @-> returning gen)
  let rgv_polint = foreign "RgV_polint" (gen @-> gen @-> long @-> returning gen)
  let vec_reduce = foreign "vec_reduce" (gen @-> ptr gen @-> returning gen)
  let rgxq_reverse = foreign "RgXQ_reverse" (gen @-> gen @-> returning gen)

  let zc_union_shallow =
    foreign "ZC_union_shallow" (gen @-> gen @-> returning gen)

  let zv_indexsort = foreign "ZV_indexsort" (gen @-> returning gen)
  let zv_search = foreign "ZV_search" (gen @-> gen @-> returning long)
  let zv_sort = foreign "ZV_sort" (gen @-> returning gen)
  let zv_sort_inplace = foreign "ZV_sort_inplace" (gen @-> returning void)
  let zv_sort_shallow = foreign "ZV_sort_shallow" (gen @-> returning gen)
  let zv_sort_uniq = foreign "ZV_sort_uniq" (gen @-> returning gen)

  let zv_sort_uniq_shallow =
    foreign "ZV_sort_uniq_shallow" (gen @-> returning gen)

  let zv_union_shallow =
    foreign "ZV_union_shallow" (gen @-> gen @-> returning gen)

  let binomial = foreign "binomial" (gen @-> long @-> returning gen)
  let binomial0 = foreign "binomial0" (gen @-> gen @-> returning gen)

  let binomialuu =
    foreign "binomialuu" (pari_ulong @-> pari_ulong @-> returning gen)

  let cmp_flx = foreign "cmp_Flx" (gen @-> gen @-> returning int)
  let cmp_rgx = foreign "cmp_RgX" (gen @-> gen @-> returning int)

  let cmp_nodata =
    foreign "cmp_nodata" (ptr void @-> gen @-> gen @-> returning int)

  let cmp_prime_ideal = foreign "cmp_prime_ideal" (gen @-> gen @-> returning int)

  let cmp_prime_over_p =
    foreign "cmp_prime_over_p" (gen @-> gen @-> returning int)

  let cmp_universal = foreign "cmp_universal" (gen @-> gen @-> returning int)
  let convol = foreign "convol" (gen @-> gen @-> returning gen)

  let gen_cmp_rgx =
    foreign "gen_cmp_RgX" (ptr void @-> gen @-> gen @-> returning int)

  let polcyclo = foreign "polcyclo" (long @-> long @-> returning gen)
  let polcyclo_eval = foreign "polcyclo_eval" (long @-> gen @-> returning gen)
  let dirdiv = foreign "dirdiv" (gen @-> gen @-> returning gen)
  let dirmul = foreign "dirmul" (gen @-> gen @-> returning gen)
  let eulerianpol = foreign "eulerianpol" (long @-> long @-> returning gen)
  let gprec_wensure = foreign "gprec_wensure" (gen @-> long @-> returning gen)

  let gen_indexsort =
    foreign "gen_indexsort"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> returning gen)

  let gen_indexsort_uniq =
    foreign "gen_indexsort_uniq"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> returning gen)

  let gen_search =
    foreign "gen_search"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> returning long)

  let gen_setminus =
    foreign "gen_setminus"
      (gen @-> gen
      @-> static_funptr Ctypes.(gen @-> gen @-> returning int)
      @-> returning gen)

  let gen_sort =
    foreign "gen_sort"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> returning gen)

  let gen_sort_inplace =
    foreign "gen_sort_inplace"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> ptr gen @-> returning void)

  let gen_sort_shallow =
    foreign "gen_sort_shallow"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> returning gen)

  let gen_sort_uniq =
    foreign "gen_sort_uniq"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> returning gen)

  let getstack = foreign "getstack" (void @-> returning long)
  let gettime = foreign "gettime" (void @-> returning long)
end

module F26 (F : Ctypes.FOREIGN) = struct
  open F

  let getabstime = foreign "getabstime" (void @-> returning long)
  let getwalltime = foreign "getwalltime" (void @-> returning gen)
  let gprec = foreign "gprec" (gen @-> long @-> returning gen)
  let gprec_wtrunc = foreign "gprec_wtrunc" (gen @-> long @-> returning gen)
  let gprec_w = foreign "gprec_w" (gen @-> long @-> returning gen)
  let gtoset = foreign "gtoset" (gen @-> returning gen)
  let indexlexsort = foreign "indexlexsort" (gen @-> returning gen)
  let indexsort = foreign "indexsort" (gen @-> returning gen)
  let indexvecsort = foreign "indexvecsort" (gen @-> gen @-> returning gen)
  let laplace = foreign "laplace" (gen @-> returning gen)
  let lexsort = foreign "lexsort" (gen @-> returning gen)
  let mathilbert = foreign "mathilbert" (long @-> returning gen)
  let matqpascal = foreign "matqpascal" (long @-> gen @-> returning gen)

  let merge_factor =
    foreign "merge_factor"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> returning gen)

  let merge_sort_uniq =
    foreign "merge_sort_uniq"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> returning gen)

  let modreverse = foreign "modreverse" (gen @-> returning gen)
  let polhermite = foreign "polhermite" (long @-> long @-> returning gen)

  let polhermite_eval0 =
    foreign "polhermite_eval0" (long @-> gen @-> long @-> returning gen)

  let polhermite_eval =
    foreign "polhermite_eval" (long @-> gen @-> returning gen)

  let pollaguerre =
    foreign "pollaguerre" (long @-> gen @-> long @-> returning gen)

  let pollaguerre_eval =
    foreign "pollaguerre_eval" (long @-> gen @-> gen @-> returning gen)

  let pollaguerre_eval0 =
    foreign "pollaguerre_eval0" (long @-> gen @-> gen @-> long @-> returning gen)

  let pollegendre = foreign "pollegendre" (long @-> long @-> returning gen)

  let pollegendre_reduced =
    foreign "pollegendre_reduced" (long @-> long @-> returning gen)

  let pollegendre_eval =
    foreign "pollegendre_eval" (long @-> gen @-> returning gen)

  let pollegendre_eval0 =
    foreign "pollegendre_eval0" (long @-> gen @-> long @-> returning gen)

  let polint =
    foreign "polint" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let polint_i =
    foreign "polint_i" (gen @-> gen @-> gen @-> ptr long @-> returning gen)

  let polintspec =
    foreign "polintspec"
      (gen @-> gen @-> gen @-> long @-> ptr long @-> returning gen)

  let polchebyshev =
    foreign "polchebyshev" (long @-> long @-> long @-> returning gen)

  let polchebyshev_eval =
    foreign "polchebyshev_eval" (long @-> long @-> gen @-> returning gen)

  let polchebyshev1 = foreign "polchebyshev1" (long @-> long @-> returning gen)
  let polchebyshev2 = foreign "polchebyshev2" (long @-> long @-> returning gen)
  let polrecip = foreign "polrecip" (gen @-> returning gen)
  let setbinop = foreign "setbinop" (gen @-> gen @-> gen @-> returning gen)
  let setdelta = foreign "setdelta" (gen @-> gen @-> returning gen)
  let setintersect = foreign "setintersect" (gen @-> gen @-> returning gen)
  let setisset = foreign "setisset" (gen @-> returning long)
  let setminus = foreign "setminus" (gen @-> gen @-> returning gen)
  let setsearch = foreign "setsearch" (gen @-> gen @-> long @-> returning long)
  let setunion = foreign "setunion" (gen @-> gen @-> returning gen)
  let setunion_i = foreign "setunion_i" (gen @-> gen @-> returning gen)
  let sort = foreign "sort" (gen @-> returning gen)

  let sort_factor =
    foreign "sort_factor"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning int)
      @-> returning gen)

  let stirling = foreign "stirling" (long @-> long @-> long @-> returning gen)

  let stirling1 =
    foreign "stirling1" (pari_ulong @-> pari_ulong @-> returning gen)

  let stirling2 =
    foreign "stirling2" (pari_ulong @-> pari_ulong @-> returning gen)

  let tablesearch =
    foreign "tablesearch"
      (gen @-> gen
      @-> static_funptr Ctypes.(gen @-> gen @-> returning int)
      @-> returning long)

  let vecbinomial = foreign "vecbinomial" (long @-> returning gen)
  let vecsearch = foreign "vecsearch" (gen @-> gen @-> gen @-> returning long)
  let vecsort = foreign "vecsort" (gen @-> gen @-> returning gen)
  let vecsort0 = foreign "vecsort0" (gen @-> gen @-> long @-> returning gen)
  let zv_search = foreign "zv_search" (gen @-> long @-> returning long)
  let bits_to_int = foreign "bits_to_int" (gen @-> long @-> returning gen)
  let bits_to_u = foreign "bits_to_u" (gen @-> long @-> returning pari_ulong)
  let binaire = foreign "binaire" (gen @-> returning gen)
  let binary_2k = foreign "binary_2k" (gen @-> long @-> returning gen)
  let binary_2k_nv = foreign "binary_2k_nv" (gen @-> long @-> returning gen)
  let binary_zv = foreign "binary_zv" (gen @-> returning gen)
  let bittest = foreign "bittest" (gen @-> long @-> returning long)
  let fromdigits_2k = foreign "fromdigits_2k" (gen @-> long @-> returning gen)
  let gbitand = foreign "gbitand" (gen @-> gen @-> returning gen)
  let gbitneg = foreign "gbitneg" (gen @-> long @-> returning gen)
  let gbitnegimply = foreign "gbitnegimply" (gen @-> gen @-> returning gen)
  let gbitor = foreign "gbitor" (gen @-> gen @-> returning gen)
  let gbittest = foreign "gbittest" (gen @-> long @-> returning gen)
  let gbitxor = foreign "gbitxor" (gen @-> gen @-> returning gen)
  let hammingl = foreign "hammingl" (pari_ulong @-> returning long)
  let hammingweight = foreign "hammingweight" (gen @-> returning long)
  let ibitand = foreign "ibitand" (gen @-> gen @-> returning gen)
  let ibitnegimply = foreign "ibitnegimply" (gen @-> gen @-> returning gen)
  let ibitor = foreign "ibitor" (gen @-> gen @-> returning gen)
  let ibitxor = foreign "ibitxor" (gen @-> gen @-> returning gen)

  let nv_fromdigits_2k =
    foreign "nv_fromdigits_2k" (gen @-> long @-> returning gen)

  let bnflogef = foreign "bnflogef" (gen @-> gen @-> returning gen)
  let bnflog = foreign "bnflog" (gen @-> gen @-> returning gen)

  let bnflogdegree =
    foreign "bnflogdegree" (gen @-> gen @-> gen @-> returning gen)

  let nfislocalpower =
    foreign "nfislocalpower" (gen @-> gen @-> gen @-> gen @-> returning long)

  let rnfislocalcyclo = foreign "rnfislocalcyclo" (gen @-> returning long)
  let bnfisunit = foreign "bnfisunit" (gen @-> gen @-> returning gen)
  let bnfissunit = foreign "bnfissunit" (gen @-> gen @-> gen @-> returning gen)
  let bnfsunit = foreign "bnfsunit" (gen @-> gen @-> long @-> returning gen)
  let bnfunits = foreign "bnfunits" (gen @-> gen @-> returning gen)
  let bnfisunit0 = foreign "bnfisunit0" (gen @-> gen @-> gen @-> returning gen)

  let sunits_mod_units =
    foreign "sunits_mod_units" (gen @-> gen @-> returning gen)

  let buchquad =
    foreign "Buchquad" (gen @-> double @-> double @-> long @-> returning gen)

  let quadclassno = foreign "quadclassno" (gen @-> returning gen)
  let quadclassnos = foreign "quadclassnos" (long @-> returning long)

  let quadclassunit0 =
    foreign "quadclassunit0" (gen @-> long @-> gen @-> long @-> returning gen)

  let buchall = foreign "Buchall" (gen @-> long @-> long @-> returning gen)

  let buchall_param =
    foreign "Buchall_param"
      (gen @-> double @-> double @-> long @-> long @-> long @-> returning gen)

  let bnf_build_cheapfu = foreign "bnf_build_cheapfu" (gen @-> returning gen)
  let bnf_build_cycgen = foreign "bnf_build_cycgen" (gen @-> returning gen)
  let bnf_build_matalpha = foreign "bnf_build_matalpha" (gen @-> returning gen)
  let bnf_build_units = foreign "bnf_build_units" (gen @-> returning gen)
  let bnf_compactfu = foreign "bnf_compactfu" (gen @-> returning gen)
  let bnf_compactfu_mat = foreign "bnf_compactfu_mat" (gen @-> returning gen)
  let bnf_has_fu = foreign "bnf_has_fu" (gen @-> returning gen)

  let bnfinit0 =
    foreign "bnfinit0" (gen @-> long @-> gen @-> long @-> returning gen)

  let bnfisprincipal0 =
    foreign "bnfisprincipal0" (gen @-> gen @-> long @-> returning gen)
end

module F27 (F : Ctypes.FOREIGN) = struct
  open F

  let bnfnewprec = foreign "bnfnewprec" (gen @-> long @-> returning gen)

  let bnfnewprec_shallow =
    foreign "bnfnewprec_shallow" (gen @-> long @-> returning gen)

  let bnftestprimes = foreign "bnftestprimes" (gen @-> gen @-> returning void)
  let bnrnewprec = foreign "bnrnewprec" (gen @-> long @-> returning gen)

  let bnrnewprec_shallow =
    foreign "bnrnewprec_shallow" (gen @-> long @-> returning gen)

  let isprincipalfact =
    foreign "isprincipalfact"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let isprincipalfact_or_fail =
    foreign "isprincipalfact_or_fail"
      (gen @-> gen @-> gen @-> gen @-> returning gen)

  let isprincipal = foreign "isprincipal" (gen @-> gen @-> returning gen)

  let nf_cxlog_normalize =
    foreign "nf_cxlog_normalize" (gen @-> gen @-> long @-> returning gen)

  let nfcyclotomicunits =
    foreign "nfcyclotomicunits" (gen @-> gen @-> returning gen)

  let nfsign_units =
    foreign "nfsign_units" (gen @-> gen @-> int @-> returning gen)

  let nfsign_tu = foreign "nfsign_tu" (gen @-> gen @-> returning gen)
  let nfsign_fu = foreign "nfsign_fu" (gen @-> gen @-> returning gen)
  let signunits = foreign "signunits" (gen @-> returning gen)
  let hermite_bound = foreign "Hermite_bound" (long @-> long @-> returning gen)

  let bnr_subgroup_sanitize =
    foreign "bnr_subgroup_sanitize" (ptr gen @-> ptr gen @-> returning void)

  let bnr_char_sanitize =
    foreign "bnr_char_sanitize" (ptr gen @-> ptr gen @-> returning void)

  let abc_to_bnr =
    foreign "ABC_to_bnr"
      (gen @-> gen @-> gen @-> ptr gen @-> int @-> returning gen)

  let buchray = foreign "Buchray" (gen @-> gen @-> long @-> returning gen)

  let buchraymod =
    foreign "Buchraymod" (gen @-> gen @-> long @-> gen @-> returning gen)

  let bnrautmatrix = foreign "bnrautmatrix" (gen @-> gen @-> returning gen)

  let bnr_subgroup_check =
    foreign "bnr_subgroup_check" (gen @-> gen @-> ptr gen @-> returning gen)

  let bnrchar = foreign "bnrchar" (gen @-> gen @-> gen @-> returning gen)

  let bnrchar_primitive =
    foreign "bnrchar_primitive" (gen @-> gen @-> gen @-> returning gen)

  let bnrclassno = foreign "bnrclassno" (gen @-> gen @-> returning gen)
  let bnrclassno0 = foreign "bnrclassno0" (gen @-> gen @-> gen @-> returning gen)
  let bnrclassnolist = foreign "bnrclassnolist" (gen @-> gen @-> returning gen)

  let bnrchar_primitive_raw =
    foreign "bnrchar_primitive_raw" (gen @-> gen @-> gen @-> returning gen)

  let bnrconductor_factored =
    foreign "bnrconductor_factored" (gen @-> gen @-> returning gen)

  let bnrconductor_raw =
    foreign "bnrconductor_raw" (gen @-> gen @-> returning gen)

  let bnrconductormod =
    foreign "bnrconductormod" (gen @-> gen @-> gen @-> returning gen)

  let bnrconductor0 =
    foreign "bnrconductor0" (gen @-> gen @-> gen @-> long @-> returning gen)

  let bnrconductor =
    foreign "bnrconductor" (gen @-> gen @-> long @-> returning gen)

  let bnrconductor_i =
    foreign "bnrconductor_i" (gen @-> gen @-> long @-> returning gen)

  let bnrconductorofchar =
    foreign "bnrconductorofchar" (gen @-> gen @-> returning gen)

  let bnrdisc0 =
    foreign "bnrdisc0" (gen @-> gen @-> gen @-> long @-> returning gen)

  let bnrdisc = foreign "bnrdisc" (gen @-> gen @-> long @-> returning gen)

  let bnrdisclist0 =
    foreign "bnrdisclist0" (gen @-> gen @-> gen @-> returning gen)

  let bnrgaloismatrix = foreign "bnrgaloismatrix" (gen @-> gen @-> returning gen)

  let bnrgaloisapply =
    foreign "bnrgaloisapply" (gen @-> gen @-> gen @-> returning gen)

  let bnrinit0 = foreign "bnrinit0" (gen @-> gen @-> long @-> returning gen)

  let bnrinitmod =
    foreign "bnrinitmod" (gen @-> gen @-> long @-> gen @-> returning gen)

  let bnrisconductor0 =
    foreign "bnrisconductor0" (gen @-> gen @-> gen @-> returning long)

  let bnrisconductor = foreign "bnrisconductor" (gen @-> gen @-> returning long)

  let bnrisgalois =
    foreign "bnrisgalois" (gen @-> gen @-> gen @-> returning long)

  let bnrisprincipalmod =
    foreign "bnrisprincipalmod" (gen @-> gen @-> gen @-> long @-> returning gen)

  let bnrisprincipal =
    foreign "bnrisprincipal" (gen @-> gen @-> long @-> returning gen)

  let bnrmap = foreign "bnrmap" (gen @-> gen @-> returning gen)
  let bnrsurjection = foreign "bnrsurjection" (gen @-> gen @-> returning gen)
  let bnfnarrow = foreign "bnfnarrow" (gen @-> returning gen)
  let bnfcertify = foreign "bnfcertify" (gen @-> returning long)
  let bnfcertify0 = foreign "bnfcertify0" (gen @-> long @-> returning long)
  let bnrcompositum = foreign "bnrcompositum" (gen @-> gen @-> returning gen)
  let decodemodule = foreign "decodemodule" (gen @-> gen @-> returning gen)
  let discrayabslist = foreign "discrayabslist" (gen @-> gen @-> returning gen)

  let discrayabslistarch =
    foreign "discrayabslistarch" (gen @-> gen @-> pari_ulong @-> returning gen)

  let idealmoddivisor = foreign "idealmoddivisor" (gen @-> gen @-> returning gen)
  let isprincipalray = foreign "isprincipalray" (gen @-> gen @-> returning gen)

  let isprincipalraygen =
    foreign "isprincipalraygen" (gen @-> gen @-> returning gen)

  let nf_deg1_prime = foreign "nf_deg1_prime" (gen @-> returning gen)
  let nfarchstar = foreign "nfarchstar" (gen @-> gen @-> gen @-> returning gen)
  let rnfconductor = foreign "rnfconductor" (gen @-> gen @-> returning gen)

  let rnfconductor0 =
    foreign "rnfconductor0" (gen @-> gen @-> long @-> returning gen)

  let rnfnormgroup = foreign "rnfnormgroup" (gen @-> gen @-> returning gen)

  let subgrouplist0 =
    foreign "subgrouplist0" (gen @-> gen @-> long @-> returning gen)

  let bnfisnorm = foreign "bnfisnorm" (gen @-> gen @-> long @-> returning gen)
  let rnfisnorm = foreign "rnfisnorm" (gen @-> gen @-> long @-> returning gen)

  let rnfisnorminit =
    foreign "rnfisnorminit" (gen @-> gen @-> int @-> returning gen)

  let coprimes_zv = foreign "coprimes_zv" (pari_ulong @-> returning gen)
  let char_check = foreign "char_check" (gen @-> gen @-> returning int)
  let charconj = foreign "charconj" (gen @-> gen @-> returning gen)
  let charconj0 = foreign "charconj0" (gen @-> gen @-> returning gen)
  let chardiv = foreign "chardiv" (gen @-> gen @-> gen @-> returning gen)
  let chardiv0 = foreign "chardiv0" (gen @-> gen @-> gen @-> returning gen)

  let chareval =
    foreign "chareval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let chargalois = foreign "chargalois" (gen @-> gen @-> returning gen)
  let charker = foreign "charker" (gen @-> gen @-> returning gen)
  let charker0 = foreign "charker0" (gen @-> gen @-> returning gen)
  let charmul = foreign "charmul" (gen @-> gen @-> gen @-> returning gen)
  let charmul0 = foreign "charmul0" (gen @-> gen @-> gen @-> returning gen)
  let charorder = foreign "charorder" (gen @-> gen @-> returning gen)
  let charorder0 = foreign "charorder0" (gen @-> gen @-> returning gen)
  let charpow = foreign "charpow" (gen @-> gen @-> gen @-> returning gen)
  let charpow0 = foreign "charpow0" (gen @-> gen @-> gen @-> returning gen)

  let char_denormalize =
    foreign "char_denormalize" (gen @-> gen @-> gen @-> returning gen)

  let char_normalize = foreign "char_normalize" (gen @-> gen @-> returning gen)
  let char_simplify = foreign "char_simplify" (gen @-> gen @-> returning gen)
  let checkznstar_i = foreign "checkznstar_i" (gen @-> returning int)
  let cyc_normalize = foreign "cyc_normalize" (gen @-> returning gen)
  let ncharvecexpo = foreign "ncharvecexpo" (gen @-> gen @-> returning gen)
  let znchar = foreign "znchar" (gen @-> returning gen)
  let znchar_quad = foreign "znchar_quad" (gen @-> gen @-> returning gen)
  let zncharcheck = foreign "zncharcheck" (gen @-> gen @-> returning int)
  let zncharconductor = foreign "zncharconductor" (gen @-> gen @-> returning gen)
  let zncharconj = foreign "zncharconj" (gen @-> gen @-> returning gen)

  let znchardecompose =
    foreign "znchardecompose" (gen @-> gen @-> gen @-> returning gen)

  let znchardiv = foreign "znchardiv" (gen @-> gen @-> gen @-> returning gen)

  let znchareval =
    foreign "znchareval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let zncharinduce =
    foreign "zncharinduce" (gen @-> gen @-> gen @-> returning gen)

  let zncharisodd = foreign "zncharisodd" (gen @-> gen @-> returning long)
end

module F28 (F : Ctypes.FOREIGN) = struct
  open F

  let zncharker = foreign "zncharker" (gen @-> gen @-> returning gen)
  let zncharmul = foreign "zncharmul" (gen @-> gen @-> gen @-> returning gen)
  let zncharorder = foreign "zncharorder" (gen @-> gen @-> returning gen)
  let zncharpow = foreign "zncharpow" (gen @-> gen @-> gen @-> returning gen)

  let znchartokronecker =
    foreign "znchartokronecker" (gen @-> gen @-> long @-> returning gen)

  let znchartoprimitive =
    foreign "znchartoprimitive" (gen @-> gen @-> returning gen)

  let znconrey_check = foreign "znconrey_check" (gen @-> gen @-> returning int)

  let znconrey_normalized =
    foreign "znconrey_normalized" (gen @-> gen @-> returning gen)

  let znconreychar = foreign "znconreychar" (gen @-> gen @-> returning gen)

  let znconreyfromchar_normalized =
    foreign "znconreyfromchar_normalized" (gen @-> gen @-> returning gen)

  let znconreyconductor =
    foreign "znconreyconductor" (gen @-> gen @-> ptr gen @-> returning gen)

  let znconreyexp = foreign "znconreyexp" (gen @-> gen @-> returning gen)

  let znconreyfromchar =
    foreign "znconreyfromchar" (gen @-> gen @-> returning gen)

  let znconreylog = foreign "znconreylog" (gen @-> gen @-> returning gen)

  let znconreylog_normalize =
    foreign "znconreylog_normalize" (gen @-> gen @-> returning gen)

  let znlog0 = foreign "znlog0" (gen @-> gen @-> gen @-> returning gen)

  let zv_cyc_minimal =
    foreign "zv_cyc_minimal" (gen @-> gen @-> gen @-> returning long)

  let zv_cyc_minimize =
    foreign "zv_cyc_minimize" (gen @-> gen @-> gen @-> returning long)

  let closure_deriv = foreign "closure_deriv" (gen @-> returning gen)
  let closure_derivn = foreign "closure_derivn" (gen @-> long @-> returning gen)

  let localvars_find =
    foreign "localvars_find" (gen @-> ptr entree @-> returning long)

  let localvars_read_str =
    foreign "localvars_read_str" (string @-> gen @-> returning gen)

  let snm_closure = foreign "snm_closure" (ptr entree @-> gen @-> returning gen)
  let strtoclosure = foreign "strtoclosure" (string @-> long @-> returning gen)
  let strtofunction = foreign "strtofunction" (string @-> returning gen)
  let gconcat = foreign "gconcat" (gen @-> gen @-> returning gen)
  let gconcat1 = foreign "gconcat1" (gen @-> returning gen)
  let matconcat = foreign "matconcat" (gen @-> returning gen)
  let shallowconcat = foreign "shallowconcat" (gen @-> gen @-> returning gen)
  let shallowconcat1 = foreign "shallowconcat1" (gen @-> returning gen)
  let shallowmatconcat = foreign "shallowmatconcat" (gen @-> returning gen)
  let vconcat = foreign "vconcat" (gen @-> gen @-> returning gen)
  let default0 = foreign "default0" (string @-> string @-> returning gen)
  let getrealprecision = foreign "getrealprecision" (void @-> returning long)

  let pari_is_default =
    foreign "pari_is_default" (string @-> returning (ptr entree))

  let sd_texstyle = foreign "sd_TeXstyle" (string @-> long @-> returning gen)
  let sd_colors = foreign "sd_colors" (string @-> long @-> returning gen)
  let sd_compatible = foreign "sd_compatible" (string @-> long @-> returning gen)
  let sd_datadir = foreign "sd_datadir" (string @-> long @-> returning gen)
  let sd_debug = foreign "sd_debug" (string @-> long @-> returning gen)
  let sd_debugfiles = foreign "sd_debugfiles" (string @-> long @-> returning gen)
  let sd_debugmem = foreign "sd_debugmem" (string @-> long @-> returning gen)

  let sd_factor_add_primes =
    foreign "sd_factor_add_primes" (string @-> long @-> returning gen)

  let sd_factor_proven =
    foreign "sd_factor_proven" (string @-> long @-> returning gen)

  let sd_format = foreign "sd_format" (string @-> long @-> returning gen)
  let sd_histsize = foreign "sd_histsize" (string @-> long @-> returning gen)
  let sd_log = foreign "sd_log" (string @-> long @-> returning gen)
  let sd_logfile = foreign "sd_logfile" (string @-> long @-> returning gen)
  let sd_nbthreads = foreign "sd_nbthreads" (string @-> long @-> returning gen)

  let sd_new_galois_format =
    foreign "sd_new_galois_format" (string @-> long @-> returning gen)

  let sd_output = foreign "sd_output" (string @-> long @-> returning gen)
  let sd_parisize = foreign "sd_parisize" (string @-> long @-> returning gen)

  let sd_parisizemax =
    foreign "sd_parisizemax" (string @-> long @-> returning gen)

  let sd_path = foreign "sd_path" (string @-> long @-> returning gen)
  let sd_plothsizes = foreign "sd_plothsizes" (string @-> long @-> returning gen)

  let sd_prettyprinter =
    foreign "sd_prettyprinter" (string @-> long @-> returning gen)

  let sd_primelimit = foreign "sd_primelimit" (string @-> long @-> returning gen)

  let sd_realbitprecision =
    foreign "sd_realbitprecision" (string @-> long @-> returning gen)

  let sd_realprecision =
    foreign "sd_realprecision" (string @-> long @-> returning gen)

  let sd_secure = foreign "sd_secure" (string @-> long @-> returning gen)

  let sd_seriesprecision =
    foreign "sd_seriesprecision" (string @-> long @-> returning gen)

  let sd_simplify = foreign "sd_simplify" (string @-> long @-> returning gen)
  let sd_sopath = foreign "sd_sopath" (string @-> int @-> returning gen)
  let sd_strictargs = foreign "sd_strictargs" (string @-> long @-> returning gen)

  let sd_strictmatch =
    foreign "sd_strictmatch" (string @-> long @-> returning gen)

  let sd_string =
    foreign "sd_string"
      (string @-> long @-> string @-> ptr string @-> returning gen)

  let sd_threadsize = foreign "sd_threadsize" (string @-> long @-> returning gen)

  let sd_threadsizemax =
    foreign "sd_threadsizemax" (string @-> long @-> returning gen)

  let sd_intarray =
    foreign "sd_intarray"
      (string @-> long @-> ptr gen @-> string @-> returning gen)

  let sd_toggle =
    foreign "sd_toggle"
      (string @-> long @-> string @-> ptr int @-> returning gen)

  let sd_ulong =
    foreign "sd_ulong"
      (string @-> long @-> string @-> ptr pari_ulong @-> pari_ulong
     @-> pari_ulong @-> ptr string @-> returning gen)

  let setdefault =
    foreign "setdefault" (string @-> string @-> long @-> returning gen)

  let setrealprecision =
    foreign "setrealprecision" (long @-> ptr long @-> returning long)

  let digits = foreign "digits" (gen @-> gen @-> returning gen)
  let fromdigits = foreign "fromdigits" (gen @-> gen @-> returning gen)
  let fromdigitsu = foreign "fromdigitsu" (gen @-> gen @-> returning gen)

  let gen_digits =
    foreign "gen_digits"
      (gen @-> gen @-> long @-> ptr void @-> ptr bb_ring
      @-> static_funptr
            Ctypes.(ptr void @-> gen @-> gen @-> ptr gen @-> returning gen)
      @-> returning gen)

  let gen_fromdigits =
    foreign "gen_fromdigits"
      (gen @-> gen @-> ptr void @-> ptr bb_ring @-> returning gen)

  let sumdigits = foreign "sumdigits" (gen @-> returning gen)
  let sumdigits0 = foreign "sumdigits0" (gen @-> gen @-> returning gen)
  let sumdigitsu = foreign "sumdigitsu" (pari_ulong @-> returning pari_ulong)
  let ecpp = foreign "ecpp" (gen @-> returning gen)
  let ecpp0 = foreign "ecpp0" (gen @-> long @-> returning gen)
  let ecppexport = foreign "ecppexport" (gen @-> long @-> returning gen)
  let ecppisvalid = foreign "ecppisvalid" (gen @-> returning long)
  let isprimeecpp = foreign "isprimeECPP" (gen @-> returning long)
  let sd_breakloop = foreign "sd_breakloop" (string @-> long @-> returning gen)
  let sd_echo = foreign "sd_echo" (string @-> long @-> returning gen)

  let sd_graphcolormap =
    foreign "sd_graphcolormap" (string @-> long @-> returning gen)

  let sd_graphcolors =
    foreign "sd_graphcolors" (string @-> long @-> returning gen)

  let sd_help = foreign "sd_help" (string @-> long @-> returning gen)
  let sd_histfile = foreign "sd_histfile" (string @-> long @-> returning gen)
  let sd_lines = foreign "sd_lines" (string @-> long @-> returning gen)
  let sd_linewrap = foreign "sd_linewrap" (string @-> long @-> returning gen)
  let sd_prompt = foreign "sd_prompt" (string @-> long @-> returning gen)

  let sd_prompt_cont =
    foreign "sd_prompt_cont" (string @-> long @-> returning gen)

  let sd_psfile = foreign "sd_psfile" (string @-> long @-> returning gen)
  let sd_readline = foreign "sd_readline" (string @-> long @-> returning gen)
  let sd_recover = foreign "sd_recover" (string @-> long @-> returning gen)
  let sd_timer = foreign "sd_timer" (string @-> long @-> returning gen)
end

module F29 (F : Ctypes.FOREIGN) = struct
  open F

  let pari_hit_return = foreign "pari_hit_return" (void @-> returning void)
  let gp_load_gprc = foreign "gp_load_gprc" (void @-> returning void)
  let gp_meta = foreign "gp_meta" (string @-> int @-> returning int)

  let gphelp_keyword_list =
    foreign "gphelp_keyword_list" (void @-> returning (ptr string))

  let pari_center = foreign "pari_center" (string @-> returning void)
  let pari_community = foreign "pari_community" (void @-> returning long)
  let pari_print_version = foreign "pari_print_version" (void @-> returning void)
  let gp_format_time = foreign "gp_format_time" (long @-> returning string)
  let gp_format_prompt = foreign "gp_format_prompt" (string @-> returning string)
  let pari_alarm = foreign "pari_alarm" (long @-> returning void)
  let gp_alarm = foreign "gp_alarm" (long @-> gen @-> returning gen)
  let gp_input = foreign "gp_input" (void @-> returning gen)
  let gp_allocatemem = foreign "gp_allocatemem" (gen @-> returning void)

  let gp_handle_exception =
    foreign "gp_handle_exception" (long @-> returning int)

  let gp_alarm_handler = foreign "gp_alarm_handler" (int @-> returning void)
  let gp_sigint_fun = foreign "gp_sigint_fun" (void @-> returning void)
  let gp_help = foreign "gp_help" (string @-> long @-> returning void)

  let gp_echo_and_log =
    foreign "gp_echo_and_log" (string @-> string @-> returning void)

  let print_fun_list =
    foreign "print_fun_list" (ptr string @-> long @-> returning void)

  let strtime = foreign "strtime" (long @-> returning gen)

  let direuler =
    foreign "direuler"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> gen @-> returning gen)

  let dirpowers = foreign "dirpowers" (long @-> gen @-> long @-> returning gen)

  let dirpowerssum =
    foreign "dirpowerssum"
      (pari_ulong @-> gen @-> long @-> long @-> returning gen)

  let dirpowerssumfun =
    foreign "dirpowerssumfun"
      (pari_ulong @-> gen @-> ptr void
      @-> static_funptr
            Ctypes.(ptr void @-> pari_ulong @-> long @-> returning gen)
      @-> long @-> long @-> returning gen)

  let vecpowuu = foreign "vecpowuu" (long @-> pari_ulong @-> returning gen)
  let vecpowug = foreign "vecpowug" (long @-> gen @-> long @-> returning gen)

  let ellanalyticrank =
    foreign "ellanalyticrank" (gen @-> gen @-> long @-> returning gen)

  let ellanalyticrank_bitprec =
    foreign "ellanalyticrank_bitprec" (gen @-> gen @-> long @-> returning gen)

  let ellanal_globalred_all =
    foreign "ellanal_globalred_all"
      (gen @-> ptr gen @-> ptr gen @-> ptr gen @-> returning gen)

  let ellheegner = foreign "ellheegner" (gen @-> returning gen)
  let elll1 = foreign "ellL1" (gen @-> long @-> long @-> returning gen)

  let elll1_bitprec =
    foreign "ellL1_bitprec" (gen @-> long @-> long @-> returning gen)

  let ellconvertname = foreign "ellconvertname" (gen @-> returning gen)
  let elldatagenerators = foreign "elldatagenerators" (gen @-> returning gen)
  let ellidentify = foreign "ellidentify" (gen @-> returning gen)
  let ellsearch = foreign "ellsearch" (gen @-> returning gen)
  let ellsearchcurve = foreign "ellsearchcurve" (gen @-> returning gen)

  let forell =
    foreign "forell"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> long @-> long @-> long @-> returning void)

  let ellfromeqn = foreign "ellfromeqn" (gen @-> returning gen)
  let akell = foreign "akell" (gen @-> gen @-> returning gen)

  let bilhell =
    foreign "bilhell" (gen @-> gen @-> gen @-> long @-> returning gen)

  let checkell = foreign "checkell" (gen @-> returning void)
  let checkell_fq = foreign "checkell_Fq" (gen @-> returning void)
  let checkell_q = foreign "checkell_Q" (gen @-> returning void)
  let checkell_qp = foreign "checkell_Qp" (gen @-> returning void)
  let checkellisog = foreign "checkellisog" (gen @-> returning void)
  let checkellpt = foreign "checkellpt" (gen @-> returning void)
  let checkell5 = foreign "checkell5" (gen @-> returning void)
  let cxredsl2 = foreign "cxredsl2" (gen @-> ptr gen @-> returning gen)

  let cxredsl2_i =
    foreign "cxredsl2_i" (gen @-> ptr gen @-> ptr gen @-> returning gen)

  let ec_2divpol_evalx =
    foreign "ec_2divpol_evalx" (gen @-> gen @-> returning gen)

  let ec_3divpol_evalx =
    foreign "ec_3divpol_evalx" (gen @-> gen @-> returning gen)

  let ec_bmodel = foreign "ec_bmodel" (gen @-> returning gen)
  let ec_f_evalx = foreign "ec_f_evalx" (gen @-> gen @-> returning gen)
  let ec_h_evalx = foreign "ec_h_evalx" (gen @-> gen @-> returning gen)
  let ec_dfdx_evalq = foreign "ec_dFdx_evalQ" (gen @-> gen @-> returning gen)
  let ec_dfdy_evalq = foreign "ec_dFdy_evalQ" (gen @-> gen @-> returning gen)
  let ec_dmfdy_evalq = foreign "ec_dmFdy_evalQ" (gen @-> gen @-> returning gen)

  let ec_half_deriv_2divpol_evalx =
    foreign "ec_half_deriv_2divpol_evalx" (gen @-> gen @-> returning gen)

  let ec_phi2 = foreign "ec_phi2" (gen @-> returning gen)
  let ell_is_integral = foreign "ell_is_integral" (gen @-> returning int)
  let ellq_get_cm = foreign "ellQ_get_CM" (gen @-> returning long)
  let ellq_get_n = foreign "ellQ_get_N" (gen @-> returning gen)

  let ellq_get_nfa =
    foreign "ellQ_get_Nfa" (gen @-> ptr gen @-> ptr gen @-> returning void)

  let ellqp_tate_uniformization =
    foreign "ellQp_Tate_uniformization" (gen @-> long @-> returning gen)

  let ellqp_agm = foreign "ellQp_AGM" (gen @-> long @-> returning gen)
  let ellqp_u = foreign "ellQp_u" (gen @-> long @-> returning gen)
  let ellqp_u2 = foreign "ellQp_u2" (gen @-> long @-> returning gen)
  let ellqp_q = foreign "ellQp_q" (gen @-> long @-> returning gen)
  let ellqp_ab = foreign "ellQp_ab" (gen @-> long @-> returning gen)
  let ellqp_l = foreign "ellQp_L" (gen @-> long @-> returning gen)
  let ellqp_root = foreign "ellQp_root" (gen @-> long @-> returning gen)

  let ellqtwist_bsdperiod =
    foreign "ellQtwist_bsdperiod" (gen @-> long @-> returning gen)

  let ellr_area = foreign "ellR_area" (gen @-> long @-> returning gen)
  let ellr_ab = foreign "ellR_ab" (gen @-> long @-> returning gen)
  let ellr_eta = foreign "ellR_eta" (gen @-> long @-> returning gen)
  let ellr_omega = foreign "ellR_omega" (gen @-> long @-> returning gen)
  let ellr_roots = foreign "ellR_roots" (gen @-> long @-> returning gen)
  let elladd = foreign "elladd" (gen @-> gen @-> gen @-> returning gen)
  let ellan = foreign "ellan" (gen @-> long @-> returning gen)
  let ellanq_zv = foreign "ellanQ_zv" (gen @-> long @-> returning gen)

  let ellanal_globalred =
    foreign "ellanal_globalred" (gen @-> ptr gen @-> returning gen)

  let ellap = foreign "ellap" (gen @-> gen @-> returning gen)

  let ellap_cm_fast =
    foreign "ellap_CM_fast" (gen @-> pari_ulong @-> long @-> returning long)

  let ellbasechar = foreign "ellbasechar" (gen @-> returning gen)
  let ellbsd = foreign "ellbsd" (gen @-> long @-> returning gen)
  let ellcard = foreign "ellcard" (gen @-> gen @-> returning gen)
  let ellchangecurve = foreign "ellchangecurve" (gen @-> gen @-> returning gen)
  let ellchangeinvert = foreign "ellchangeinvert" (gen @-> returning gen)
  let ellchangepoint = foreign "ellchangepoint" (gen @-> gen @-> returning gen)

  let ellchangepointinv =
    foreign "ellchangepointinv" (gen @-> gen @-> returning gen)

  let elldivpol = foreign "elldivpol" (gen @-> long @-> long @-> returning gen)

  let elleisnum =
    foreign "elleisnum" (gen @-> long @-> long @-> long @-> returning gen)

  let elleta = foreign "elleta" (gen @-> long @-> returning gen)
  let elleulerf = foreign "elleulerf" (gen @-> gen @-> returning gen)
  let ellff_get_card = foreign "ellff_get_card" (gen @-> returning gen)
  let ellff_get_gens = foreign "ellff_get_gens" (gen @-> returning gen)
  let ellff_get_group = foreign "ellff_get_group" (gen @-> returning gen)
  let ellff_get_o = foreign "ellff_get_o" (gen @-> returning gen)
  let ellff_get_p = foreign "ellff_get_p" (gen @-> returning gen)
end

module F30 (F : Ctypes.FOREIGN) = struct
  open F

  let ellff_get_m = foreign "ellff_get_m" (gen @-> returning gen)
  let ellff_get_d = foreign "ellff_get_D" (gen @-> returning gen)
  let ellfromj = foreign "ellfromj" (gen @-> returning gen)
  let ellgenerators = foreign "ellgenerators" (gen @-> returning gen)
  let ellglobalred = foreign "ellglobalred" (gen @-> returning gen)
  let ellgroup = foreign "ellgroup" (gen @-> gen @-> returning gen)
  let ellgroup0 = foreign "ellgroup0" (gen @-> gen @-> long @-> returning gen)

  let ellheight0 =
    foreign "ellheight0" (gen @-> gen @-> gen @-> long @-> returning gen)

  let ellheight = foreign "ellheight" (gen @-> gen @-> long @-> returning gen)

  let ellheightmatrix =
    foreign "ellheightmatrix" (gen @-> gen @-> long @-> returning gen)

  let ellheightoo =
    foreign "ellheightoo" (gen @-> gen @-> long @-> returning gen)

  let ellinit = foreign "ellinit" (gen @-> gen @-> long @-> returning gen)

  let ellintegralmodel =
    foreign "ellintegralmodel" (gen @-> ptr gen @-> returning gen)

  let ellintegralmodel_i =
    foreign "ellintegralmodel_i" (gen @-> ptr gen @-> returning gen)

  let ellisoncurve = foreign "ellisoncurve" (gen @-> gen @-> returning gen)
  let ellisotree = foreign "ellisotree" (gen @-> returning gen)

  let ellissupersingular =
    foreign "ellissupersingular" (gen @-> gen @-> returning int)

  let elljissupersingular = foreign "elljissupersingular" (gen @-> returning int)

  let elllseries =
    foreign "elllseries" (gen @-> gen @-> gen @-> long @-> returning gen)

  let elllocalred = foreign "elllocalred" (gen @-> gen @-> returning gen)
  let elllog = foreign "elllog" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let ellminimaldisc = foreign "ellminimaldisc" (gen @-> returning gen)

  let ellminimalmodel =
    foreign "ellminimalmodel" (gen @-> ptr gen @-> returning gen)

  let ellminimaltwist = foreign "ellminimaltwist" (gen @-> returning gen)

  let ellminimaltwist0 =
    foreign "ellminimaltwist0" (gen @-> long @-> returning gen)

  let ellminimaltwistcond = foreign "ellminimaltwistcond" (gen @-> returning gen)
  let ellmul = foreign "ellmul" (gen @-> gen @-> gen @-> returning gen)
  let ellnf_vecarea = foreign "ellnf_vecarea" (gen @-> long @-> returning gen)
  let ellnf_veceta = foreign "ellnf_veceta" (gen @-> long @-> returning gen)
  let ellnf_vecomega = foreign "ellnf_vecomega" (gen @-> long @-> returning gen)
  let ellneg = foreign "ellneg" (gen @-> gen @-> returning gen)
  let ellorder = foreign "ellorder" (gen @-> gen @-> gen @-> returning gen)
  let ellorder_q = foreign "ellorder_Q" (gen @-> gen @-> returning long)

  let ellordinate =
    foreign "ellordinate" (gen @-> gen @-> long @-> returning gen)

  let ellpadicheight0 =
    foreign "ellpadicheight0"
      (gen @-> gen @-> long @-> gen @-> gen @-> returning gen)

  let ellpadicheightmatrix =
    foreign "ellpadicheightmatrix"
      (gen @-> gen @-> long @-> gen @-> returning gen)

  let ellperiods = foreign "ellperiods" (gen @-> long @-> long @-> returning gen)
  let ellrandom = foreign "ellrandom" (gen @-> returning gen)
  let ellrootno = foreign "ellrootno" (gen @-> gen @-> returning long)
  let ellrootno_global = foreign "ellrootno_global" (gen @-> returning long)

  let ellsaturation =
    foreign "ellsaturation" (gen @-> gen @-> long @-> long @-> returning gen)

  let ellsea = foreign "ellsea" (gen @-> long @-> returning gen)

  let ellsigma =
    foreign "ellsigma" (gen @-> gen @-> long @-> long @-> returning gen)

  let ellsub = foreign "ellsub" (gen @-> gen @-> gen @-> returning gen)
  let elltamagawa = foreign "elltamagawa" (gen @-> returning gen)
  let elltaniyama = foreign "elltaniyama" (gen @-> long @-> returning gen)

  let elltatepairing =
    foreign "elltatepairing" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let elltors = foreign "elltors" (gen @-> returning gen)
  let elltors0 = foreign "elltors0" (gen @-> long @-> returning gen)

  let elltors_psylow =
    foreign "elltors_psylow" (gen @-> pari_ulong @-> returning gen)

  let elltrace = foreign "elltrace" (gen @-> gen @-> returning gen)
  let elltwist = foreign "elltwist" (gen @-> gen @-> returning gen)

  let ellweilpairing =
    foreign "ellweilpairing" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let ellwp = foreign "ellwp" (gen @-> gen @-> long @-> returning gen)
  let ellwp0 = foreign "ellwp0" (gen @-> gen @-> long @-> long @-> returning gen)

  let ellwpseries =
    foreign "ellwpseries" (gen @-> long @-> long @-> returning gen)

  let ellxn = foreign "ellxn" (gen @-> long @-> long @-> returning gen)
  let ellzeta = foreign "ellzeta" (gen @-> gen @-> long @-> returning gen)
  let oncurve = foreign "oncurve" (gen @-> gen @-> returning int)
  let orderell = foreign "orderell" (gen @-> gen @-> returning gen)
  let pointell = foreign "pointell" (gen @-> gen @-> long @-> returning gen)

  let point_to_a4a6 =
    foreign "point_to_a4a6" (gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let point_to_a4a6_fl =
    foreign "point_to_a4a6_Fl"
      (gen @-> gen @-> pari_ulong @-> ptr pari_ulong @-> returning gen)

  let zell = foreign "zell" (gen @-> gen @-> long @-> returning gen)

  let qp_agm2_sequence =
    foreign "Qp_agm2_sequence" (gen @-> gen @-> returning gen)

  let qp_ascending_landen =
    foreign "Qp_ascending_Landen"
      (gen @-> ptr gen @-> ptr gen @-> returning void)

  let qp_descending_landen =
    foreign "Qp_descending_Landen"
      (gen @-> ptr gen @-> ptr gen @-> returning void)

  let ellformaldifferential =
    foreign "ellformaldifferential" (gen @-> long @-> long @-> returning gen)

  let ellformalexp =
    foreign "ellformalexp" (gen @-> long @-> long @-> returning gen)

  let ellformallog =
    foreign "ellformallog" (gen @-> long @-> long @-> returning gen)

  let ellformalpoint =
    foreign "ellformalpoint" (gen @-> long @-> long @-> returning gen)

  let ellformalw = foreign "ellformalw" (gen @-> long @-> long @-> returning gen)

  let ellnonsingularmultiple =
    foreign "ellnonsingularmultiple" (gen @-> gen @-> returning gen)

  let ellpadicl =
    foreign "ellpadicL"
      (gen @-> gen @-> long @-> gen @-> long @-> gen @-> returning gen)

  let ellpadicbsd =
    foreign "ellpadicbsd" (gen @-> gen @-> long @-> gen @-> returning gen)

  let ellpadicfrobenius =
    foreign "ellpadicfrobenius" (gen @-> pari_ulong @-> long @-> returning gen)

  let ellpadicheight =
    foreign "ellpadicheight" (gen @-> gen @-> long @-> gen @-> returning gen)

  let ellpadiclog =
    foreign "ellpadiclog" (gen @-> gen @-> long @-> gen @-> returning gen)

  let ellpadicregulator =
    foreign "ellpadicregulator" (gen @-> gen @-> long @-> gen @-> returning gen)

  let ellpadics2 = foreign "ellpadics2" (gen @-> gen @-> long @-> returning gen)
  let ell2cover = foreign "ell2cover" (gen @-> long @-> returning gen)

  let ellrank =
    foreign "ellrank" (gen @-> long @-> gen @-> long @-> returning gen)

  let ellrankinit = foreign "ellrankinit" (gen @-> long @-> returning gen)

  let hyperell_locally_soluble =
    foreign "hyperell_locally_soluble" (gen @-> gen @-> returning long)

  let nf_hyperell_locally_soluble =
    foreign "nf_hyperell_locally_soluble"
      (gen @-> gen @-> gen @-> returning long)

  let nfhilbert = foreign "nfhilbert" (gen @-> gen @-> gen @-> returning long)

  let nfhilbert0 =
    foreign "nfhilbert0" (gen @-> gen @-> gen @-> gen @-> returning long)

  let ellisdivisible =
    foreign "ellisdivisible" (gen @-> gen @-> gen @-> ptr gen @-> returning long)

  let ellisogenyapply = foreign "ellisogenyapply" (gen @-> gen @-> returning gen)

  let ellisogeny =
    foreign "ellisogeny"
      (gen @-> gen @-> long @-> long @-> long @-> returning gen)

  let ellisomat = foreign "ellisomat" (gen @-> long @-> long @-> returning gen)
  let ellweilcurve = foreign "ellweilcurve" (gen @-> ptr gen @-> returning gen)

  let flxq_elldivpolmod =
    foreign "Flxq_elldivpolmod"
      (gen @-> gen @-> long @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let fp_ellcard_sea =
    foreign "Fp_ellcard_SEA" (gen @-> gen @-> gen @-> long @-> returning gen)

  let fq_ellcard_sea =
    foreign "Fq_ellcard_SEA"
      (gen @-> gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let fq_elldivpolmod =
    foreign "Fq_elldivpolmod"
      (gen @-> gen @-> long @-> gen @-> gen @-> gen @-> returning gen)

  let ellmodulareqn =
    foreign "ellmodulareqn" (long @-> long @-> long @-> returning gen)

  let externstr = foreign "externstr" (string @-> returning gen)
  let gp_filter = foreign "gp_filter" (string @-> returning string)
  let gpextern = foreign "gpextern" (string @-> returning gen)
end

module F31 (F : Ctypes.FOREIGN) = struct
  open F

  let gpsystem = foreign "gpsystem" (string @-> returning long)
  let readstr = foreign "readstr" (string @-> returning gen)
  let gentogenstr_nospace = foreign "GENtoGENstr_nospace" (gen @-> returning gen)
  let gentogenstr = foreign "GENtoGENstr" (gen @-> returning gen)
  let gentotexstr = foreign "GENtoTeXstr" (gen @-> returning string)
  let gentostr = foreign "GENtostr" (gen @-> returning string)
  let gentostr_raw = foreign "GENtostr_raw" (gen @-> returning string)
  let gentostr_unquoted = foreign "GENtostr_unquoted" (gen @-> returning string)
  let str = foreign "Str" (gen @-> returning gen)
  let strexpand = foreign "strexpand" (gen @-> returning gen)
  let strtex = foreign "strtex" (gen @-> returning gen)
  let brute = foreign "brute" (gen @-> char @-> long @-> returning void)
  let dbggen = foreign "dbgGEN" (gen @-> long @-> returning void)
  let error0 = foreign "error0" (gen @-> returning void)
  let dbg_pari_heap = foreign "dbg_pari_heap" (void @-> returning void)
  let err_flush = foreign "err_flush" (void @-> returning void)
  let err_printf = foreign "err_printf" (string @-> returning void)
  let gp_getenv = foreign "gp_getenv" (string @-> returning gen)
  let gp_fileclose = foreign "gp_fileclose" (long @-> returning void)
  let gp_fileextern = foreign "gp_fileextern" (string @-> returning long)
  let gp_fileflush = foreign "gp_fileflush" (long @-> returning void)
  let gp_fileflush0 = foreign "gp_fileflush0" (gen @-> returning void)
  let gp_fileopen = foreign "gp_fileopen" (string @-> string @-> returning long)
  let gp_fileread = foreign "gp_fileread" (long @-> returning gen)
  let gp_filereadstr = foreign "gp_filereadstr" (long @-> returning gen)
  let gp_filewrite = foreign "gp_filewrite" (long @-> string @-> returning void)

  let gp_filewrite1 =
    foreign "gp_filewrite1" (long @-> string @-> returning void)

  let gp_read_file = foreign "gp_read_file" (string @-> returning gen)

  let gp_read_str_multiline =
    foreign "gp_read_str_multiline" (string @-> string @-> returning gen)

  let gp_readvec_file = foreign "gp_readvec_file" (string @-> returning gen)

  let gpinstall =
    foreign "gpinstall"
      (string @-> string @-> string @-> string @-> returning void)

  let gsprintf = foreign "gsprintf" (string @-> returning gen)
  let itostr = foreign "itostr" (gen @-> returning string)
  let matbrute = foreign "matbrute" (gen @-> char @-> long @-> returning void)
  let os_getenv = foreign "os_getenv" (string @-> returning string)
  let uordinal = foreign "uordinal" (pari_ulong @-> returning string)
  let outmat = foreign "outmat" (gen @-> returning void)
  let output = foreign "output" (gen @-> returning void)
  let rgv_to_str = foreign "RgV_to_str" (gen @-> long @-> returning string)

  let pari_add_hist =
    foreign "pari_add_hist" (gen @-> long @-> long @-> returning void)

  let pari_ask_confirm = foreign "pari_ask_confirm" (string @-> returning void)
  let pari_flush = foreign "pari_flush" (void @-> returning void)

  let pari_fopen =
    foreign "pari_fopen" (string @-> string @-> returning (ptr parifile))

  let pari_fopen_or_fail =
    foreign "pari_fopen_or_fail" (string @-> string @-> returning (ptr parifile))

  let pari_fopengz = foreign "pari_fopengz" (string @-> returning (ptr parifile))
  let pari_get_hist = foreign "pari_get_hist" (long @-> returning gen)
  let pari_get_histrtime = foreign "pari_get_histrtime" (long @-> returning long)
  let pari_get_histtime = foreign "pari_get_histtime" (long @-> returning long)
  let pari_get_homedir = foreign "pari_get_homedir" (string @-> returning string)
  let pari_histtime = foreign "pari_histtime" (long @-> returning gen)
  let pari_is_dir = foreign "pari_is_dir" (string @-> returning int)
  let pari_is_file = foreign "pari_is_file" (string @-> returning int)

  let pari_last_was_newline =
    foreign "pari_last_was_newline" (void @-> returning int)

  let pari_set_last_newline =
    foreign "pari_set_last_newline" (int @-> returning void)

  let pari_nb_hist = foreign "pari_nb_hist" (void @-> returning pari_ulong)
  let pari_printf = foreign "pari_printf" (string @-> returning void)
  let pari_putc = foreign "pari_putc" (char @-> returning void)
  let pari_puts = foreign "pari_puts" (string @-> returning void)

  let pari_safefopen =
    foreign "pari_safefopen" (string @-> string @-> returning (ptr parifile))

  let pari_sprintf = foreign "pari_sprintf" (string @-> returning string)
  let pari_stdin_isatty = foreign "pari_stdin_isatty" (void @-> returning int)
  let pari_unique_dir = foreign "pari_unique_dir" (string @-> returning string)

  let pari_unique_filename =
    foreign "pari_unique_filename" (string @-> returning string)

  let pari_unique_filename_suffix =
    foreign "pari_unique_filename_suffix"
      (string @-> string @-> returning string)

  let pari_unlink = foreign "pari_unlink" (string @-> returning void)
  let path_expand = foreign "path_expand" (string @-> returning string)

  let pari_sprint0 =
    foreign "pari_sprint0" (string @-> gen @-> long @-> returning string)

  let print = foreign "print" (gen @-> returning void)
  let printp = foreign "printp" (gen @-> returning void)
  let print1 = foreign "print1" (gen @-> returning void)
  let printf0 = foreign "printf0" (string @-> gen @-> returning void)
  let printsep = foreign "printsep" (string @-> gen @-> returning void)
  let printsep1 = foreign "printsep1" (string @-> gen @-> returning void)
  let printtex = foreign "printtex" (gen @-> returning void)
  let stack_sprintf = foreign "stack_sprintf" (string @-> returning string)
  let str_init = foreign "str_init" (ptr pari_str @-> int @-> returning void)

  let str_printf =
    foreign "str_printf" (ptr pari_str @-> string @-> returning void)

  let str_putc = foreign "str_putc" (ptr pari_str @-> char @-> returning void)
  let str_puts = foreign "str_puts" (ptr pari_str @-> string @-> returning void)

  let strftime_expand =
    foreign "strftime_expand" (string @-> string @-> long @-> returning void)

  let strprintf = foreign "strprintf" (string @-> gen @-> returning gen)
  let term_color = foreign "term_color" (long @-> returning void)

  let term_get_color =
    foreign "term_get_color" (string @-> long @-> returning string)

  let texe = foreign "texe" (gen @-> char @-> long @-> returning void)
  let warning0 = foreign "warning0" (gen @-> returning void)
  let write0 = foreign "write0" (string @-> gen @-> returning void)
  let write1 = foreign "write1" (string @-> gen @-> returning void)
  let writebin = foreign "writebin" (string @-> gen @-> returning void)
  let writetex = foreign "writetex" (string @-> gen @-> returning void)
  let bincopy_relink = foreign "bincopy_relink" (gen @-> gen @-> returning void)
  let bitprecision0 = foreign "bitprecision0" (gen @-> long @-> returning gen)
  let bitprecision00 = foreign "bitprecision00" (gen @-> gen @-> returning gen)
  let break0 = foreign "break0" (long @-> returning gen)
  let call0 = foreign "call0" (gen @-> gen @-> returning gen)

  let closure_callgen0prec =
    foreign "closure_callgen0prec" (gen @-> long @-> returning gen)

  let closure_callgen1 =
    foreign "closure_callgen1" (gen @-> gen @-> returning gen)

  let closure_callgen1prec =
    foreign "closure_callgen1prec" (gen @-> gen @-> long @-> returning gen)

  let closure_callgen2 =
    foreign "closure_callgen2" (gen @-> gen @-> gen @-> returning gen)

  let closure_callgenall =
    foreign "closure_callgenall" (gen @-> long @-> returning gen)

  let closure_callgenvec =
    foreign "closure_callgenvec" (gen @-> gen @-> returning gen)
end

module F32 (F : Ctypes.FOREIGN) = struct
  open F

  let closure_callgenvecdef =
    foreign "closure_callgenvecdef" (gen @-> gen @-> gen @-> returning gen)

  let closure_callgenvecdefprec =
    foreign "closure_callgenvecdefprec"
      (gen @-> gen @-> gen @-> long @-> returning gen)

  let closure_callgenvecprec =
    foreign "closure_callgenvecprec" (gen @-> gen @-> long @-> returning gen)

  let closure_callvoid1 =
    foreign "closure_callvoid1" (gen @-> gen @-> returning void)

  let closure_context =
    foreign "closure_context" (long @-> long @-> returning long)

  let closure_disassemble =
    foreign "closure_disassemble" (gen @-> returning void)

  let closure_err = foreign "closure_err" (long @-> returning void)

  let closure_evalbrk =
    foreign "closure_evalbrk" (gen @-> ptr long @-> returning gen)

  let closure_evalgen = foreign "closure_evalgen" (gen @-> returning gen)
  let closure_evalnobrk = foreign "closure_evalnobrk" (gen @-> returning gen)
  let closure_evalres = foreign "closure_evalres" (gen @-> returning gen)
  let closure_evalvoid = foreign "closure_evalvoid" (gen @-> returning void)
  let closure_func_err = foreign "closure_func_err" (void @-> returning string)

  let closure_trapgen =
    foreign "closure_trapgen" (gen @-> long @-> returning gen)

  let copybin_unlink = foreign "copybin_unlink" (gen @-> returning gen)
  let getlocalprec = foreign "getlocalprec" (long @-> returning long)
  let getlocalbitprec = foreign "getlocalbitprec" (long @-> returning long)
  let get_lex = foreign "get_lex" (long @-> returning gen)
  let get_localprec = foreign "get_localprec" (void @-> returning long)
  let get_localbitprec = foreign "get_localbitprec" (void @-> returning long)
  let gp_call = foreign "gp_call" (ptr void @-> gen @-> returning gen)

  let gp_callprec =
    foreign "gp_callprec" (ptr void @-> gen @-> long @-> returning gen)

  let gp_call2 = foreign "gp_call2" (ptr void @-> gen @-> gen @-> returning gen)
  let gp_callbool = foreign "gp_callbool" (ptr void @-> gen @-> returning long)
  let gp_callvoid = foreign "gp_callvoid" (ptr void @-> gen @-> returning long)
  let gp_eval = foreign "gp_eval" (ptr void @-> gen @-> returning gen)
  let gp_evalbool = foreign "gp_evalbool" (ptr void @-> gen @-> returning long)

  let gp_evalprec =
    foreign "gp_evalprec" (ptr void @-> gen @-> long @-> returning gen)

  let gp_evalupto = foreign "gp_evalupto" (ptr void @-> gen @-> returning gen)
  let gp_evalvoid = foreign "gp_evalvoid" (ptr void @-> gen @-> returning long)
  let localprec = foreign "localprec" (gen @-> returning void)
  let localbitprec = foreign "localbitprec" (gen @-> returning void)
  let loop_break = foreign "loop_break" (void @-> returning long)
  let next0 = foreign "next0" (long @-> returning gen)
  let pareval = foreign "pareval" (gen @-> returning gen)
  let pari_self = foreign "pari_self" (void @-> returning gen)
  let parsum = foreign "parsum" (gen @-> gen @-> gen @-> returning gen)
  let parvector = foreign "parvector" (long @-> gen @-> returning gen)
  let pop_lex = foreign "pop_lex" (long @-> returning void)
  let pop_localprec = foreign "pop_localprec" (void @-> returning void)
  let precision0 = foreign "precision0" (gen @-> long @-> returning gen)
  let precision00 = foreign "precision00" (gen @-> gen @-> returning gen)
  let push_lex = foreign "push_lex" (gen @-> gen @-> returning void)
  let push_localbitprec = foreign "push_localbitprec" (long @-> returning void)
  let push_localprec = foreign "push_localprec" (long @-> returning void)
  let return0 = foreign "return0" (gen @-> returning gen)
  let set_lex = foreign "set_lex" (long @-> gen @-> returning void)

  let forcomposite_init =
    foreign "forcomposite_init"
      (ptr forcomposite_t @-> gen @-> gen @-> returning int)

  let forcomposite_next =
    foreign "forcomposite_next" (ptr forcomposite_t @-> returning gen)

  let forprime_next = foreign "forprime_next" (ptr forprime_t @-> returning gen)

  let forprime_init =
    foreign "forprime_init" (ptr forprime_t @-> gen @-> gen @-> returning int)

  let forprimestep_init =
    foreign "forprimestep_init"
      (ptr forprime_t @-> gen @-> gen @-> gen @-> returning int)

  let initprimes =
    foreign "initprimes"
      (pari_ulong @-> ptr long @-> ptr pari_ulong @-> returning byteptr)

  let initprimetable = foreign "initprimetable" (pari_ulong @-> returning void)

  let init_primepointer_geq =
    foreign "init_primepointer_geq"
      (pari_ulong @-> ptr byteptr @-> returning pari_ulong)

  let init_primepointer_gt =
    foreign "init_primepointer_gt"
      (pari_ulong @-> ptr byteptr @-> returning pari_ulong)

  let init_primepointer_leq =
    foreign "init_primepointer_leq"
      (pari_ulong @-> ptr byteptr @-> returning pari_ulong)

  let init_primepointer_lt =
    foreign "init_primepointer_lt"
      (pari_ulong @-> ptr byteptr @-> returning pari_ulong)

  let maxprime = foreign "maxprime" (void @-> returning pari_ulong)
  let maxprimen = foreign "maxprimeN" (void @-> returning pari_ulong)
  let maxprime_check = foreign "maxprime_check" (pari_ulong @-> returning void)
  let maxprimelim = foreign "maxprimelim" (void @-> returning pari_ulong)

  let pari_init_primes =
    foreign "pari_init_primes" (pari_ulong @-> returning void)

  let prodprimes = foreign "prodprimes" (void @-> returning gen)

  let u_forprime_next =
    foreign "u_forprime_next" (ptr forprime_t @-> returning pari_ulong)

  let u_forprime_init =
    foreign "u_forprime_init"
      (ptr forprime_t @-> pari_ulong @-> pari_ulong @-> returning int)

  let u_forprime_restrict =
    foreign "u_forprime_restrict"
      (ptr forprime_t @-> pari_ulong @-> returning void)

  let u_forprime_arith_init =
    foreign "u_forprime_arith_init"
      (ptr forprime_t @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> pari_ulong @-> returning int)

  let ff_1 = foreign "FF_1" (gen @-> returning gen)
  let ff_frobenius = foreign "FF_Frobenius" (gen @-> long @-> returning gen)

  let ff_z_z_muldiv =
    foreign "FF_Z_Z_muldiv" (gen @-> gen @-> gen @-> returning gen)

  let ff_q_add = foreign "FF_Q_add" (gen @-> gen @-> returning gen)
  let ff_z_add = foreign "FF_Z_add" (gen @-> gen @-> returning gen)
  let ff_z_mul = foreign "FF_Z_mul" (gen @-> gen @-> returning gen)
  let ff_add = foreign "FF_add" (gen @-> gen @-> returning gen)
  let ff_charpoly = foreign "FF_charpoly" (gen @-> returning gen)
  let ff_conjvec = foreign "FF_conjvec" (gen @-> returning gen)
  let ff_div = foreign "FF_div" (gen @-> gen @-> returning gen)
  let ff_ellcard = foreign "FF_ellcard" (gen @-> returning gen)
  let ff_ellcard_sea = foreign "FF_ellcard_SEA" (gen @-> long @-> returning gen)
  let ff_ellgens = foreign "FF_ellgens" (gen @-> returning gen)
  let ff_ellgroup = foreign "FF_ellgroup" (gen @-> ptr gen @-> returning gen)

  let ff_elllog =
    foreign "FF_elllog" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let ff_ellmul = foreign "FF_ellmul" (gen @-> gen @-> gen @-> returning gen)
  let ff_ellorder = foreign "FF_ellorder" (gen @-> gen @-> gen @-> returning gen)
  let ff_elltwist = foreign "FF_elltwist" (gen @-> returning gen)
  let ff_ellrandom = foreign "FF_ellrandom" (gen @-> returning gen)

  let ff_elltatepairing =
    foreign "FF_elltatepairing" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let ff_ellweilpairing =
    foreign "FF_ellweilpairing" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let ff_equal = foreign "FF_equal" (gen @-> gen @-> returning int)
  let ff_equal0 = foreign "FF_equal0" (gen @-> returning int)
  let ff_equal1 = foreign "FF_equal1" (gen @-> returning int)
  let ff_equalm1 = foreign "FF_equalm1" (gen @-> returning int)
  let ff_f = foreign "FF_f" (gen @-> returning long)
  let ff_gen = foreign "FF_gen" (gen @-> returning gen)
  let ff_inv = foreign "FF_inv" (gen @-> returning gen)
  let ff_issquare = foreign "FF_issquare" (gen @-> returning long)

  let ff_issquareall =
    foreign "FF_issquareall" (gen @-> ptr gen @-> returning long)

  let ff_ispower =
    foreign "FF_ispower" (gen @-> gen @-> ptr gen @-> returning long)

  let ff_log = foreign "FF_log" (gen @-> gen @-> gen @-> returning gen)
end

module F33 (F : Ctypes.FOREIGN) = struct
  open F

  let ff_map = foreign "FF_map" (gen @-> gen @-> returning gen)
  let ff_minpoly = foreign "FF_minpoly" (gen @-> returning gen)
  let ff_mod = foreign "FF_mod" (gen @-> returning gen)
  let ff_mul = foreign "FF_mul" (gen @-> gen @-> returning gen)
  let ff_mul2n = foreign "FF_mul2n" (gen @-> long @-> returning gen)
  let ff_neg = foreign "FF_neg" (gen @-> returning gen)
  let ff_neg_i = foreign "FF_neg_i" (gen @-> returning gen)
  let ff_norm = foreign "FF_norm" (gen @-> returning gen)
  let ff_order = foreign "FF_order" (gen @-> gen @-> returning gen)
  let ff_p = foreign "FF_p" (gen @-> returning gen)
  let ff_p_i = foreign "FF_p_i" (gen @-> returning gen)
  let ff_pow = foreign "FF_pow" (gen @-> gen @-> returning gen)
  let ff_primroot = foreign "FF_primroot" (gen @-> ptr gen @-> returning gen)
  let ff_q = foreign "FF_q" (gen @-> returning gen)
  let ff_samefield = foreign "FF_samefield" (gen @-> gen @-> returning int)
  let ff_sqr = foreign "FF_sqr" (gen @-> returning gen)
  let ff_sqrt = foreign "FF_sqrt" (gen @-> returning gen)
  let ff_sqrtn = foreign "FF_sqrtn" (gen @-> gen @-> ptr gen @-> returning gen)
  let ff_sub = foreign "FF_sub" (gen @-> gen @-> returning gen)
  let ff_to_f2xq = foreign "FF_to_F2xq" (gen @-> returning gen)
  let ff_to_f2xq_i = foreign "FF_to_F2xq_i" (gen @-> returning gen)
  let ff_to_flxq = foreign "FF_to_Flxq" (gen @-> returning gen)
  let ff_to_flxq_i = foreign "FF_to_Flxq_i" (gen @-> returning gen)
  let ff_to_fpxq = foreign "FF_to_FpXQ" (gen @-> returning gen)
  let ff_to_fpxq_i = foreign "FF_to_FpXQ_i" (gen @-> returning gen)
  let ff_trace = foreign "FF_trace" (gen @-> returning gen)
  let ff_var = foreign "FF_var" (gen @-> returning long)
  let ff_zero = foreign "FF_zero" (gen @-> returning gen)

  let ffm_ffc_invimage =
    foreign "FFM_FFC_invimage" (gen @-> gen @-> gen @-> returning gen)

  let ffm_ffc_gauss =
    foreign "FFM_FFC_gauss" (gen @-> gen @-> gen @-> returning gen)

  let ffm_ffc_mul = foreign "FFM_FFC_mul" (gen @-> gen @-> gen @-> returning gen)
  let ffm_deplin = foreign "FFM_deplin" (gen @-> gen @-> returning gen)
  let ffm_det = foreign "FFM_det" (gen @-> gen @-> returning gen)
  let ffm_gauss = foreign "FFM_gauss" (gen @-> gen @-> gen @-> returning gen)
  let ffm_image = foreign "FFM_image" (gen @-> gen @-> returning gen)
  let ffm_indexrank = foreign "FFM_indexrank" (gen @-> gen @-> returning gen)
  let ffm_inv = foreign "FFM_inv" (gen @-> gen @-> returning gen)

  let ffm_invimage =
    foreign "FFM_invimage" (gen @-> gen @-> gen @-> returning gen)

  let ffm_ker = foreign "FFM_ker" (gen @-> gen @-> returning gen)
  let ffm_mul = foreign "FFM_mul" (gen @-> gen @-> gen @-> returning gen)
  let ffm_rank = foreign "FFM_rank" (gen @-> gen @-> returning long)
  let ffm_suppl = foreign "FFM_suppl" (gen @-> gen @-> returning gen)
  let ffx_add = foreign "FFX_add" (gen @-> gen @-> gen @-> returning gen)
  let ffx_ddf = foreign "FFX_ddf" (gen @-> gen @-> returning gen)
  let ffx_degfact = foreign "FFX_degfact" (gen @-> gen @-> returning gen)
  let ffx_disc = foreign "FFX_disc" (gen @-> gen @-> returning gen)

  let ffx_extgcd =
    foreign "FFX_extgcd"
      (gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let ffx_factor = foreign "FFX_factor" (gen @-> gen @-> returning gen)

  let ffx_factor_squarefree =
    foreign "FFX_factor_squarefree" (gen @-> gen @-> returning gen)

  let ffx_gcd = foreign "FFX_gcd" (gen @-> gen @-> gen @-> returning gen)
  let ffx_halfgcd = foreign "FFX_halfgcd" (gen @-> gen @-> gen @-> returning gen)

  let ffx_ispower =
    foreign "FFX_ispower" (gen @-> long @-> gen @-> ptr gen @-> returning long)

  let ffx_mul = foreign "FFX_mul" (gen @-> gen @-> gen @-> returning gen)

  let ffx_preimage =
    foreign "FFX_preimage" (gen @-> gen @-> gen @-> returning gen)

  let ffx_preimagerel =
    foreign "FFX_preimagerel" (gen @-> gen @-> gen @-> returning gen)

  let ffx_rem = foreign "FFX_rem" (gen @-> gen @-> gen @-> returning gen)

  let ffx_resultant =
    foreign "FFX_resultant" (gen @-> gen @-> gen @-> returning gen)

  let ffx_roots = foreign "FFX_roots" (gen @-> gen @-> returning gen)
  let ffx_sqr = foreign "FFX_sqr" (gen @-> gen @-> returning gen)
  let ffxq_inv = foreign "FFXQ_inv" (gen @-> gen @-> gen @-> returning gen)

  let ffxq_minpoly =
    foreign "FFXQ_minpoly" (gen @-> gen @-> gen @-> returning gen)

  let ffxq_mul =
    foreign "FFXQ_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let ffxq_sqr = foreign "FFXQ_sqr" (gen @-> gen @-> gen @-> returning gen)
  let fqx_to_ffx = foreign "FqX_to_FFX" (gen @-> gen @-> returning gen)
  let fq_to_ff = foreign "Fq_to_FF" (gen @-> gen @-> returning gen)
  let z_ff_div = foreign "Z_FF_div" (gen @-> gen @-> returning gen)
  let ffembed = foreign "ffembed" (gen @-> gen @-> returning gen)
  let ffextend = foreign "ffextend" (gen @-> gen @-> long @-> returning gen)
  let fffrobenius = foreign "fffrobenius" (gen @-> long @-> returning gen)
  let ffgen = foreign "ffgen" (gen @-> long @-> returning gen)
  let ffinvmap = foreign "ffinvmap" (gen @-> returning gen)
  let fflog = foreign "fflog" (gen @-> gen @-> gen @-> returning gen)
  let ffmap = foreign "ffmap" (gen @-> gen @-> returning gen)
  let ffmaprel = foreign "ffmaprel" (gen @-> gen @-> returning gen)
  let ffcompomap = foreign "ffcompomap" (gen @-> gen @-> returning gen)
  let fforder = foreign "fforder" (gen @-> gen @-> returning gen)
  let ffprimroot = foreign "ffprimroot" (gen @-> ptr gen @-> returning gen)
  let ffrandom = foreign "ffrandom" (gen @-> returning gen)
  let rg_is_ff = foreign "Rg_is_FF" (gen @-> ptr gen @-> returning int)
  let rgc_is_ffc = foreign "RgC_is_FFC" (gen @-> ptr gen @-> returning int)
  let rgm_is_ffm = foreign "RgM_is_FFM" (gen @-> ptr gen @-> returning int)
  let p_to_ff = foreign "p_to_FF" (gen @-> long @-> returning gen)
  let tp_to_ff = foreign "Tp_to_FF" (gen @-> gen @-> returning gen)

  let flx_factcyclo =
    foreign "Flx_factcyclo"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fpx_factcyclo =
    foreign "FpX_factcyclo" (pari_ulong @-> gen @-> pari_ulong @-> returning gen)

  let factormodcyclo =
    foreign "factormodcyclo" (long @-> gen @-> long @-> long @-> returning gen)

  let checkgal = foreign "checkgal" (gen @-> returning gen)
  let checkgroup = foreign "checkgroup" (gen @-> ptr gen @-> returning gen)
  let checkgroupelts = foreign "checkgroupelts" (gen @-> returning gen)
  let embed_disc = foreign "embed_disc" (gen @-> long @-> long @-> returning gen)
  let embed_roots = foreign "embed_roots" (gen @-> long @-> returning gen)
  let galois_group = foreign "galois_group" (gen @-> returning gen)
  let galoisconj = foreign "galoisconj" (gen @-> gen @-> returning gen)

  let galoisconj0 =
    foreign "galoisconj0" (gen @-> long @-> gen @-> long @-> returning gen)

  let galoisconjclasses = foreign "galoisconjclasses" (gen @-> returning gen)
  let galoisexport = foreign "galoisexport" (gen @-> long @-> returning gen)

  let galoisfixedfield =
    foreign "galoisfixedfield" (gen @-> gen @-> long @-> long @-> returning gen)

  let galoisidentify = foreign "galoisidentify" (gen @-> returning gen)
  let galoisinit = foreign "galoisinit" (gen @-> gen @-> returning gen)

  let galoisisabelian =
    foreign "galoisisabelian" (gen @-> long @-> returning gen)
end

module F34 (F : Ctypes.FOREIGN) = struct
  open F

  let galoisisnormal = foreign "galoisisnormal" (gen @-> gen @-> returning long)
  let galoispermtopol = foreign "galoispermtopol" (gen @-> gen @-> returning gen)

  let galoissplittinginit =
    foreign "galoissplittinginit" (gen @-> gen @-> returning gen)

  let galoissubgroups = foreign "galoissubgroups" (gen @-> returning gen)

  let galoissubfields =
    foreign "galoissubfields" (gen @-> long @-> long @-> returning gen)

  let numberofconjugates =
    foreign "numberofconjugates" (gen @-> long @-> returning long)

  let polgalois = foreign "polgalois" (gen @-> long @-> returning gen)
  let galoisnbpol = foreign "galoisnbpol" (long @-> returning gen)
  let galoisgetgroup = foreign "galoisgetgroup" (long @-> long @-> returning gen)
  let galoisgetname = foreign "galoisgetname" (long @-> long @-> returning gen)

  let galoisgetpol =
    foreign "galoisgetpol" (long @-> long @-> long @-> returning gen)

  let conj_i = foreign "conj_i" (gen @-> returning gen)
  let conjvec = foreign "conjvec" (gen @-> long @-> returning gen)
  let divrunextu = foreign "divrunextu" (gen @-> pari_ulong @-> returning gen)
  let gadd = foreign "gadd" (gen @-> gen @-> returning gen)
  let gaddsg = foreign "gaddsg" (long @-> gen @-> returning gen)
  let gconj = foreign "gconj" (gen @-> returning gen)
  let gdiv = foreign "gdiv" (gen @-> gen @-> returning gen)
  let gdivgs = foreign "gdivgs" (gen @-> long @-> returning gen)
  let gdivgu = foreign "gdivgu" (gen @-> pari_ulong @-> returning gen)
  let gdivgunextu = foreign "gdivgunextu" (gen @-> pari_ulong @-> returning gen)
  let ginv = foreign "ginv" (gen @-> returning gen)
  let gmul = foreign "gmul" (gen @-> gen @-> returning gen)
  let gmul2n = foreign "gmul2n" (gen @-> long @-> returning gen)
  let gmulsg = foreign "gmulsg" (long @-> gen @-> returning gen)
  let gmulug = foreign "gmulug" (pari_ulong @-> gen @-> returning gen)
  let gsqr = foreign "gsqr" (gen @-> returning gen)
  let gsub = foreign "gsub" (gen @-> gen @-> returning gen)
  let gsubsg = foreign "gsubsg" (long @-> gen @-> returning gen)
  let mulcxi = foreign "mulcxI" (gen @-> returning gen)
  let mulcxmi = foreign "mulcxmI" (gen @-> returning gen)
  let mulcxpowis = foreign "mulcxpowIs" (gen @-> long @-> returning gen)
  let qdivii = foreign "Qdivii" (gen @-> gen @-> returning gen)
  let qdiviu = foreign "Qdiviu" (gen @-> pari_ulong @-> returning gen)
  let qdivis = foreign "Qdivis" (gen @-> long @-> returning gen)
  let ser_normalize = foreign "ser_normalize" (gen @-> returning gen)

  let gassoc_proto =
    foreign "gassoc_proto"
      (static_funptr Ctypes.(gen @-> gen @-> returning gen)
      @-> gen @-> gen @-> returning gen)

  let map_proto_g =
    foreign "map_proto_G"
      (static_funptr Ctypes.(gen @-> returning gen) @-> gen @-> returning gen)

  let map_proto_lg =
    foreign "map_proto_lG"
      (static_funptr Ctypes.(gen @-> returning long) @-> gen @-> returning gen)

  let map_proto_lgl =
    foreign "map_proto_lGL"
      (static_funptr Ctypes.(gen @-> long @-> returning long)
      @-> gen @-> long @-> returning gen)

  let q_lval = foreign "Q_lval" (gen @-> pari_ulong @-> returning long)

  let q_lvalrem =
    foreign "Q_lvalrem" (gen @-> pari_ulong @-> ptr gen @-> returning long)

  let q_pval = foreign "Q_pval" (gen @-> gen @-> returning long)

  let q_pvalrem =
    foreign "Q_pvalrem" (gen @-> gen @-> ptr gen @-> returning long)

  let rgx_val = foreign "RgX_val" (gen @-> returning long)
  let rgx_valrem = foreign "RgX_valrem" (gen @-> ptr gen @-> returning long)

  let rgx_valrem_inexact =
    foreign "RgX_valrem_inexact" (gen @-> ptr gen @-> returning long)

  let rgxv_maxdegree = foreign "RgXV_maxdegree" (gen @-> returning long)
  let zv_z_dvd = foreign "ZV_Z_dvd" (gen @-> gen @-> returning int)
  let zv_pval = foreign "ZV_pval" (gen @-> gen @-> returning long)

  let zv_pvalrem =
    foreign "ZV_pvalrem" (gen @-> gen @-> ptr gen @-> returning long)

  let zv_lval = foreign "ZV_lval" (gen @-> pari_ulong @-> returning long)

  let zv_lvalrem =
    foreign "ZV_lvalrem" (gen @-> pari_ulong @-> ptr gen @-> returning long)

  let zx_lvalrem =
    foreign "ZX_lvalrem" (gen @-> pari_ulong @-> ptr gen @-> returning long)

  let zx_lval = foreign "ZX_lval" (gen @-> pari_ulong @-> returning long)
  let zx_pval = foreign "ZX_pval" (gen @-> gen @-> returning long)

  let zx_pvalrem =
    foreign "ZX_pvalrem" (gen @-> gen @-> ptr gen @-> returning long)

  let z_lval = foreign "Z_lval" (gen @-> pari_ulong @-> returning long)

  let z_lvalrem =
    foreign "Z_lvalrem" (gen @-> pari_ulong @-> ptr gen @-> returning long)

  let z_lvalrem_stop =
    foreign "Z_lvalrem_stop"
      (ptr gen @-> pari_ulong @-> ptr int @-> returning long)

  let z_pval = foreign "Z_pval" (gen @-> gen @-> returning long)

  let z_pvalrem =
    foreign "Z_pvalrem" (gen @-> gen @-> ptr gen @-> returning long)

  let cgetp = foreign "cgetp" (gen @-> returning gen)
  let cvstop2 = foreign "cvstop2" (long @-> gen @-> returning gen)
  let cvtop = foreign "cvtop" (gen @-> gen @-> long @-> returning gen)
  let cvtop2 = foreign "cvtop2" (gen @-> gen @-> returning gen)
  let cx_approx_equal = foreign "cx_approx_equal" (gen @-> gen @-> returning int)
  let cx_approx0 = foreign "cx_approx0" (gen @-> gen @-> returning int)
  let gabs = foreign "gabs" (gen @-> long @-> returning gen)
  let gaffect = foreign "gaffect" (gen @-> gen @-> returning void)
  let gaffsg = foreign "gaffsg" (long @-> gen @-> returning void)
  let gcmp = foreign "gcmp" (gen @-> gen @-> returning int)
  let gequal0 = foreign "gequal0" (gen @-> returning int)
  let gequal1 = foreign "gequal1" (gen @-> returning int)
  let gequalx = foreign "gequalX" (gen @-> returning int)
  let gequalm1 = foreign "gequalm1" (gen @-> returning int)
  let gcmpsg = foreign "gcmpsg" (long @-> gen @-> returning int)
  let gcvtop = foreign "gcvtop" (gen @-> gen @-> long @-> returning gen)
  let gequal = foreign "gequal" (gen @-> gen @-> returning int)
  let gequalsg = foreign "gequalsg" (long @-> gen @-> returning int)
  let gexpo = foreign "gexpo" (gen @-> returning long)
  let gexpo_safe = foreign "gexpo_safe" (gen @-> returning long)
  let gpexponent = foreign "gpexponent" (gen @-> returning gen)
  let gpvaluation = foreign "gpvaluation" (gen @-> gen @-> returning gen)
  let gvaluation = foreign "gvaluation" (gen @-> gen @-> returning long)
  let gidentical = foreign "gidentical" (gen @-> gen @-> returning int)
  let glength = foreign "glength" (gen @-> returning long)
  let gmax = foreign "gmax" (gen @-> gen @-> returning gen)
  let gmaxgs = foreign "gmaxgs" (gen @-> long @-> returning gen)
  let gmin = foreign "gmin" (gen @-> gen @-> returning gen)
  let gmings = foreign "gmings" (gen @-> long @-> returning gen)
  let gneg = foreign "gneg" (gen @-> returning gen)
  let gneg_i = foreign "gneg_i" (gen @-> returning gen)
  let gsigne = foreign "gsigne" (gen @-> returning int)
  let gtolist = foreign "gtolist" (gen @-> returning gen)
  let gtolong = foreign "gtolong" (gen @-> returning long)
  let lexcmp = foreign "lexcmp" (gen @-> gen @-> returning int)
  let listinsert = foreign "listinsert" (gen @-> gen @-> long @-> returning gen)
  let listpop = foreign "listpop" (gen @-> long @-> returning void)
  let listpop0 = foreign "listpop0" (gen @-> long @-> returning void)
end

module F35 (F : Ctypes.FOREIGN) = struct
  open F

  let listput = foreign "listput" (gen @-> gen @-> long @-> returning gen)
  let listput0 = foreign "listput0" (gen @-> gen @-> long @-> returning gen)
  let listsort = foreign "listsort" (gen @-> long @-> returning void)
  let matsize = foreign "matsize" (gen @-> returning gen)
  let mklist = foreign "mklist" (void @-> returning gen)
  let mklist_typ = foreign "mklist_typ" (long @-> returning gen)
  let mklistcopy = foreign "mklistcopy" (gen @-> returning gen)
  let mkmap = foreign "mkmap" (void @-> returning gen)
  let normalizeser = foreign "normalizeser" (gen @-> returning gen)
  let normalizepol = foreign "normalizepol" (gen @-> returning gen)

  let normalizepol_approx =
    foreign "normalizepol_approx" (gen @-> long @-> returning gen)

  let normalizepol_lg =
    foreign "normalizepol_lg" (gen @-> long @-> returning gen)

  let padic_to_fl =
    foreign "padic_to_Fl" (gen @-> pari_ulong @-> returning pari_ulong)

  let padic_to_fp = foreign "padic_to_Fp" (gen @-> gen @-> returning gen)
  let quadtofp = foreign "quadtofp" (gen @-> long @-> returning gen)
  let sizedigit = foreign "sizedigit" (gen @-> returning long)
  let u_lval = foreign "u_lval" (pari_ulong @-> pari_ulong @-> returning long)

  let u_lvalrem =
    foreign "u_lvalrem"
      (pari_ulong @-> pari_ulong @-> ptr pari_ulong @-> returning long)

  let u_lvalrem_stop =
    foreign "u_lvalrem_stop"
      (ptr pari_ulong @-> pari_ulong @-> ptr int @-> returning long)

  let u_pval = foreign "u_pval" (pari_ulong @-> gen @-> returning long)

  let u_pvalrem =
    foreign "u_pvalrem"
      (pari_ulong @-> gen @-> ptr pari_ulong @-> returning long)

  let vecindexmax = foreign "vecindexmax" (gen @-> returning long)
  let vecindexmin = foreign "vecindexmin" (gen @-> returning long)
  let vecmax0 = foreign "vecmax0" (gen @-> ptr gen @-> returning gen)
  let vecmax = foreign "vecmax" (gen @-> returning gen)
  let vecmin0 = foreign "vecmin0" (gen @-> ptr gen @-> returning gen)
  let vecmin = foreign "vecmin" (gen @-> returning gen)
  let z_lval = foreign "z_lval" (long @-> pari_ulong @-> returning long)

  let z_lvalrem =
    foreign "z_lvalrem" (long @-> pari_ulong @-> ptr long @-> returning long)

  let z_pval = foreign "z_pval" (long @-> gen @-> returning long)

  let z_pvalrem =
    foreign "z_pvalrem" (long @-> gen @-> ptr long @-> returning long)

  let zx_lval = foreign "zx_lval" (gen @-> long @-> returning long)
  let hgmcyclo = foreign "hgmcyclo" (gen @-> returning gen)
  let hgmalpha = foreign "hgmalpha" (gen @-> returning gen)
  let hgmgamma = foreign "hgmgamma" (gen @-> returning gen)
  let hgminit = foreign "hgminit" (gen @-> gen @-> returning gen)
  let hgmparams = foreign "hgmparams" (gen @-> returning gen)

  let hgmeulerfactor =
    foreign "hgmeulerfactor" (gen @-> gen @-> long @-> ptr gen @-> returning gen)

  let hgmcoef = foreign "hgmcoef" (gen @-> gen @-> gen @-> returning gen)
  let hgmcoefs = foreign "hgmcoefs" (gen @-> gen @-> long @-> returning gen)
  let hgmtwist = foreign "hgmtwist" (gen @-> returning gen)
  let hgmissymmetrical = foreign "hgmissymmetrical" (gen @-> returning long)
  let hgmbydegree = foreign "hgmbydegree" (long @-> returning gen)

  let lfunhgm =
    foreign "lfunhgm" (gen @-> gen @-> gen @-> long @-> returning gen)

  let qp_zeta = foreign "Qp_zeta" (gen @-> returning gen)

  let lerchphi =
    foreign "lerchphi" (gen @-> gen @-> gen @-> long @-> returning gen)

  let lerchzeta =
    foreign "lerchzeta" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zetahurwitz =
    foreign "zetahurwitz" (gen @-> gen @-> long @-> long @-> returning gen)

  let rgx_to_ser = foreign "RgX_to_ser" (gen @-> long @-> returning gen)

  let rgx_to_ser_inexact =
    foreign "RgX_to_ser_inexact" (gen @-> long @-> returning gen)

  let gtoser = foreign "gtoser" (gen @-> long @-> long @-> returning gen)

  let gtoser_prec =
    foreign "gtoser_prec" (gen @-> long @-> long @-> returning gen)

  let rfrac_to_ser = foreign "rfrac_to_ser" (gen @-> long @-> returning gen)
  let rfrac_to_ser_i = foreign "rfrac_to_ser_i" (gen @-> long @-> returning gen)

  let rfracrecip_to_ser_absolute =
    foreign "rfracrecip_to_ser_absolute" (gen @-> long @-> returning gen)

  let rfracrecip = foreign "rfracrecip" (ptr gen @-> ptr gen @-> returning long)
  let scalarser = foreign "scalarser" (gen @-> long @-> long @-> returning gen)
  let sertoser = foreign "sertoser" (gen @-> long @-> returning gen)
  let toser_i = foreign "toser_i" (gen @-> returning gen)
  let rgv_to_ser = foreign "RgV_to_ser" (gen @-> long @-> long @-> returning gen)
  let ser0 = foreign "Ser0" (gen @-> long @-> gen @-> long @-> returning gen)
  let padic_to_q = foreign "padic_to_Q" (gen @-> returning gen)
  let padic_to_q_shallow = foreign "padic_to_Q_shallow" (gen @-> returning gen)
  let qpv_to_qv = foreign "QpV_to_QV" (gen @-> returning gen)

  let rgc_rgv_mulrealsym =
    foreign "RgC_RgV_mulrealsym" (gen @-> gen @-> returning gen)

  let rgm_mulreal = foreign "RgM_mulreal" (gen @-> gen @-> returning gen)
  let rgx_cxeval = foreign "RgX_cxeval" (gen @-> gen @-> gen @-> returning gen)

  let rgx_deflate_max =
    foreign "RgX_deflate_max" (gen @-> ptr long @-> returning gen)

  let rgx_deflate_order = foreign "RgX_deflate_order" (gen @-> returning long)
  let rgx_degree = foreign "RgX_degree" (gen @-> long @-> returning long)
  let rgx_integ = foreign "RgX_integ" (gen @-> returning gen)

  let rgxy_cxevalx =
    foreign "RgXY_cxevalx" (gen @-> gen @-> gen @-> returning gen)

  let zx_deflate_order = foreign "ZX_deflate_order" (gen @-> returning long)

  let zx_deflate_max =
    foreign "ZX_deflate_max" (gen @-> ptr long @-> returning gen)

  let ceil_safe = foreign "ceil_safe" (gen @-> returning gen)
  let ceilr = foreign "ceilr" (gen @-> returning gen)
  let centerlift = foreign "centerlift" (gen @-> returning gen)
  let centerlift0 = foreign "centerlift0" (gen @-> long @-> returning gen)
  let compo = foreign "compo" (gen @-> long @-> returning gen)
  let deg1pol = foreign "deg1pol" (gen @-> gen @-> long @-> returning gen)

  let deg1pol_shallow =
    foreign "deg1pol_shallow" (gen @-> gen @-> long @-> returning gen)

  let deg2pol_shallow =
    foreign "deg2pol_shallow" (gen @-> gen @-> gen @-> long @-> returning gen)

  let degree = foreign "degree" (gen @-> returning long)
  let denom = foreign "denom" (gen @-> returning gen)
  let denom_i = foreign "denom_i" (gen @-> returning gen)
  let denominator = foreign "denominator" (gen @-> gen @-> returning gen)
  let deriv = foreign "deriv" (gen @-> long @-> returning gen)
  let derivn = foreign "derivn" (gen @-> long @-> long @-> returning gen)
  let derivser = foreign "derivser" (gen @-> returning gen)
  let diffop = foreign "diffop" (gen @-> gen @-> gen @-> returning gen)

  let diffop0 =
    foreign "diffop0" (gen @-> gen @-> gen @-> long @-> returning gen)

  let diviiround = foreign "diviiround" (gen @-> gen @-> returning gen)
  let divrem = foreign "divrem" (gen @-> gen @-> long @-> returning gen)
  let floor_safe = foreign "floor_safe" (gen @-> returning gen)
  let gceil = foreign "gceil" (gen @-> returning gen)
  let gcvtoi = foreign "gcvtoi" (gen @-> ptr long @-> returning gen)
  let gdeflate = foreign "gdeflate" (gen @-> long @-> long @-> returning gen)
  let gdivent = foreign "gdivent" (gen @-> gen @-> returning gen)
  let gdiventgs = foreign "gdiventgs" (gen @-> long @-> returning gen)
  let gdiventsg = foreign "gdiventsg" (long @-> gen @-> returning gen)
end

module F36 (F : Ctypes.FOREIGN) = struct
  open F

  let gdiventres = foreign "gdiventres" (gen @-> gen @-> returning gen)
  let gdivmod = foreign "gdivmod" (gen @-> gen @-> ptr gen @-> returning gen)
  let gdivround = foreign "gdivround" (gen @-> gen @-> returning gen)
  let gdvd = foreign "gdvd" (gen @-> gen @-> returning int)
  let geq = foreign "geq" (gen @-> gen @-> returning gen)
  let geval = foreign "geval" (gen @-> returning gen)
  let gfloor = foreign "gfloor" (gen @-> returning gen)
  let gtrunc2n = foreign "gtrunc2n" (gen @-> long @-> returning gen)
  let gfrac = foreign "gfrac" (gen @-> returning gen)
  let gge = foreign "gge" (gen @-> gen @-> returning gen)
  let ggrando = foreign "ggrando" (gen @-> long @-> returning gen)
  let ggt = foreign "ggt" (gen @-> gen @-> returning gen)
  let gimag = foreign "gimag" (gen @-> returning gen)
  let gisexactzero = foreign "gisexactzero" (gen @-> returning gen)
  let gle = foreign "gle" (gen @-> gen @-> returning gen)
  let glt = foreign "glt" (gen @-> gen @-> returning gen)
  let gmod = foreign "gmod" (gen @-> gen @-> returning gen)
  let gmodgs = foreign "gmodgs" (gen @-> long @-> returning gen)
  let gmodsg = foreign "gmodsg" (long @-> gen @-> returning gen)
  let gmodulo = foreign "gmodulo" (gen @-> gen @-> returning gen)
  let gmodulsg = foreign "gmodulsg" (long @-> gen @-> returning gen)
  let gmodulss = foreign "gmodulss" (long @-> long @-> returning gen)
  let gne = foreign "gne" (gen @-> gen @-> returning gen)
  let gnot = foreign "gnot" (gen @-> returning gen)
  let gpolvar = foreign "gpolvar" (gen @-> returning gen)
  let gppadicprec = foreign "gppadicprec" (gen @-> gen @-> returning gen)
  let gppoldegree = foreign "gppoldegree" (gen @-> long @-> returning gen)
  let gprecision = foreign "gprecision" (gen @-> returning long)
  let gpserprec = foreign "gpserprec" (gen @-> long @-> returning gen)
  let greal = foreign "greal" (gen @-> returning gen)
  let grndtoi = foreign "grndtoi" (gen @-> ptr long @-> returning gen)
  let ground = foreign "ground" (gen @-> returning gen)
  let gshift = foreign "gshift" (gen @-> long @-> returning gen)
  let gsubst = foreign "gsubst" (gen @-> long @-> gen @-> returning gen)
  let gsubstpol = foreign "gsubstpol" (gen @-> gen @-> gen @-> returning gen)
  let gsubstvec = foreign "gsubstvec" (gen @-> gen @-> gen @-> returning gen)
  let gtocol = foreign "gtocol" (gen @-> returning gen)
  let gtocol0 = foreign "gtocol0" (gen @-> long @-> returning gen)
  let gtocolrev = foreign "gtocolrev" (gen @-> returning gen)
  let gtocolrev0 = foreign "gtocolrev0" (gen @-> long @-> returning gen)
  let gtopoly = foreign "gtopoly" (gen @-> long @-> returning gen)
  let gtopolyrev = foreign "gtopolyrev" (gen @-> long @-> returning gen)
  let gtovec = foreign "gtovec" (gen @-> returning gen)
  let gtovec0 = foreign "gtovec0" (gen @-> long @-> returning gen)
  let gtovecrev = foreign "gtovecrev" (gen @-> returning gen)
  let gtovecrev0 = foreign "gtovecrev0" (gen @-> long @-> returning gen)
  let gtovecsmall = foreign "gtovecsmall" (gen @-> returning gen)
  let gtovecsmall0 = foreign "gtovecsmall0" (gen @-> long @-> returning gen)
  let gtrunc = foreign "gtrunc" (gen @-> returning gen)
  let gvar = foreign "gvar" (gen @-> returning long)
  let gvar2 = foreign "gvar2" (gen @-> returning long)
  let hqfeval = foreign "hqfeval" (gen @-> gen @-> returning gen)
  let imag_i = foreign "imag_i" (gen @-> returning gen)
  let integ = foreign "integ" (gen @-> long @-> returning gen)
  let integser = foreign "integser" (gen @-> returning gen)
  let ser_inv = foreign "ser_inv" (gen @-> returning gen)
  let iscomplex = foreign "iscomplex" (gen @-> returning int)
  let isexactzero = foreign "isexactzero" (gen @-> returning int)

  let isrationalzeroscalar =
    foreign "isrationalzeroscalar" (gen @-> returning int)

  let isinexact = foreign "isinexact" (gen @-> returning int)
  let isinexactreal = foreign "isinexactreal" (gen @-> returning int)
  let isint = foreign "isint" (gen @-> ptr gen @-> returning int)
  let isrationalzero = foreign "isrationalzero" (gen @-> returning int)
  let issmall = foreign "issmall" (gen @-> ptr long @-> returning int)
  let lift = foreign "lift" (gen @-> returning gen)
  let lift_shallow = foreign "lift_shallow" (gen @-> returning gen)
  let lift0 = foreign "lift0" (gen @-> long @-> returning gen)
  let liftall = foreign "liftall" (gen @-> returning gen)
  let liftall_shallow = foreign "liftall_shallow" (gen @-> returning gen)
  let liftint = foreign "liftint" (gen @-> returning gen)
  let liftint_shallow = foreign "liftint_shallow" (gen @-> returning gen)
  let liftpol = foreign "liftpol" (gen @-> returning gen)
  let liftpol_shallow = foreign "liftpol_shallow" (gen @-> returning gen)
  let mkcoln = foreign "mkcoln" (long @-> returning gen)
  let mkintn = foreign "mkintn" (long @-> returning gen)
  let mkpoln = foreign "mkpoln" (long @-> returning gen)
  let mkvecn = foreign "mkvecn" (long @-> returning gen)
  let mkvecsmalln = foreign "mkvecsmalln" (long @-> returning gen)
  let modrr_safe = foreign "modRr_safe" (gen @-> gen @-> returning gen)
  let mulreal = foreign "mulreal" (gen @-> gen @-> returning gen)
  let numer = foreign "numer" (gen @-> returning gen)
  let numer_i = foreign "numer_i" (gen @-> returning gen)
  let numerator = foreign "numerator" (gen @-> gen @-> returning gen)
  let padicprec = foreign "padicprec" (gen @-> gen @-> returning long)
  let padicprec_relative = foreign "padicprec_relative" (gen @-> returning long)
  let polcoef = foreign "polcoef" (gen @-> long @-> long @-> returning gen)
  let polcoef_i = foreign "polcoef_i" (gen @-> long @-> long @-> returning gen)
  let poldegree = foreign "poldegree" (gen @-> long @-> returning long)
  let poleval = foreign "poleval" (gen @-> gen @-> returning gen)
  let pollead = foreign "pollead" (gen @-> long @-> returning gen)
  let precision = foreign "precision" (gen @-> returning long)
  let qf_apply_rgm = foreign "qf_apply_RgM" (gen @-> gen @-> returning gen)
  let qf_apply_zm = foreign "qf_apply_ZM" (gen @-> gen @-> returning gen)
  let qfb_apply_zm = foreign "qfb_apply_ZM" (gen @-> gen @-> returning gen)
  let qfbil = foreign "qfbil" (gen @-> gen @-> gen @-> returning gen)
  let qfeval = foreign "qfeval" (gen @-> gen @-> returning gen)
  let qfeval0 = foreign "qfeval0" (gen @-> gen @-> gen @-> returning gen)
  let qfevalb = foreign "qfevalb" (gen @-> gen @-> gen @-> returning gen)
  let qfnorm = foreign "qfnorm" (gen @-> gen @-> returning gen)
  let real_i = foreign "real_i" (gen @-> returning gen)
end

module F37 (F : Ctypes.FOREIGN) = struct
  open F

  let round0 = foreign "round0" (gen @-> ptr gen @-> returning gen)
  let roundr = foreign "roundr" (gen @-> returning gen)
  let roundr_safe = foreign "roundr_safe" (gen @-> returning gen)
  let scalarpol = foreign "scalarpol" (gen @-> long @-> returning gen)

  let scalarpol_shallow =
    foreign "scalarpol_shallow" (gen @-> long @-> returning gen)

  let ser_unscale = foreign "ser_unscale" (gen @-> gen @-> returning gen)
  let serprec = foreign "serprec" (gen @-> long @-> returning long)
  let serreverse = foreign "serreverse" (gen @-> returning gen)
  let simplify = foreign "simplify" (gen @-> returning gen)
  let simplify_shallow = foreign "simplify_shallow" (gen @-> returning gen)
  let tayl = foreign "tayl" (gen @-> long @-> long @-> returning gen)
  let trunc0 = foreign "trunc0" (gen @-> ptr gen @-> returning gen)
  let uu32toi = foreign "uu32toi" (pari_ulong @-> pari_ulong @-> returning gen)

  let uu32toineg =
    foreign "uu32toineg" (pari_ulong @-> pari_ulong @-> returning gen)

  let vars_sort_inplace = foreign "vars_sort_inplace" (gen @-> returning gen)
  let vars_to_rgxv = foreign "vars_to_RgXV" (gen @-> returning gen)
  let variables_vecsmall = foreign "variables_vecsmall" (gen @-> returning gen)
  let variables_vec = foreign "variables_vec" (gen @-> returning gen)
  let genus2red = foreign "genus2red" (gen @-> gen @-> returning gen)
  let genus2igusa = foreign "genus2igusa" (gen @-> long @-> returning gen)
  let gchar_conductor = foreign "gchar_conductor" (gen @-> gen @-> returning gen)

  let gchar_identify =
    foreign "gchar_identify" (gen @-> gen @-> gen @-> long @-> returning gen)

  let gcharalgebraic = foreign "gcharalgebraic" (gen @-> gen @-> returning gen)
  let gcharduallog = foreign "gcharduallog" (gen @-> gen @-> returning gen)

  let gchareval =
    foreign "gchareval" (gen @-> gen @-> gen @-> long @-> returning gen)

  let gchari_lfun = foreign "gchari_lfun" (gen @-> gen @-> gen @-> returning gen)
  let gcharinit = foreign "gcharinit" (gen @-> gen @-> long @-> returning gen)

  let gcharisalgebraic =
    foreign "gcharisalgebraic" (gen @-> gen @-> ptr gen @-> returning int)

  let gcharlocal =
    foreign "gcharlocal"
      (gen @-> gen @-> gen @-> long @-> ptr gen @-> returning gen)

  let gcharlog = foreign "gcharlog" (gen @-> gen @-> long @-> returning gen)
  let gcharnewprec = foreign "gcharnewprec" (gen @-> long @-> returning gen)
  let is_gchar_group = foreign "is_gchar_group" (gen @-> returning int)
  let lfungchar = foreign "lfungchar" (gen @-> gen @-> returning gen)

  let vecan_gchar =
    foreign "vecan_gchar" (gen @-> long @-> long @-> returning gen)

  let eulerf_gchar =
    foreign "eulerf_gchar" (gen @-> gen @-> long @-> returning gen)

  let group_ident = foreign "group_ident" (gen @-> gen @-> returning long)

  let group_ident_trans =
    foreign "group_ident_trans" (gen @-> gen @-> returning long)

  let hash_create_ulong =
    foreign "hash_create_ulong"
      (pari_ulong @-> long @-> returning (ptr hashtable))

  let hash_create_str =
    foreign "hash_create_str" (pari_ulong @-> long @-> returning (ptr hashtable))

  let hash_create =
    foreign "hash_create"
      (pari_ulong
      @-> static_funptr Ctypes.(ptr void @-> returning pari_ulong)
      @-> static_funptr Ctypes.(ptr void @-> ptr void @-> returning int)
      @-> int
      @-> returning (ptr hashtable))

  let hash_dbg = foreign "hash_dbg" (ptr hashtable @-> returning void)

  let hash_haskey_gen =
    foreign "hash_haskey_GEN" (ptr hashtable @-> ptr void @-> returning gen)

  let hash_haskey_long =
    foreign "hash_haskey_long"
      (ptr hashtable @-> ptr void @-> ptr long @-> returning int)

  let hash_init =
    foreign "hash_init"
      (ptr hashtable @-> pari_ulong
      @-> static_funptr Ctypes.(ptr void @-> returning pari_ulong)
      @-> static_funptr Ctypes.(ptr void @-> ptr void @-> returning int)
      @-> int @-> returning void)

  let hash_init_gen =
    foreign "hash_init_GEN"
      (ptr hashtable @-> pari_ulong
      @-> static_funptr Ctypes.(gen @-> gen @-> returning int)
      @-> int @-> returning void)

  let hash_init_ulong =
    foreign "hash_init_ulong"
      (ptr hashtable @-> pari_ulong @-> int @-> returning void)

  let hash_insert =
    foreign "hash_insert"
      (ptr hashtable @-> ptr void @-> ptr void @-> returning void)

  let hash_insert_long =
    foreign "hash_insert_long"
      (ptr hashtable @-> ptr void @-> long @-> returning void)

  let hash_insert2 =
    foreign "hash_insert2"
      (ptr hashtable @-> ptr void @-> ptr void @-> pari_ulong @-> returning void)

  let hash_keys = foreign "hash_keys" (ptr hashtable @-> returning gen)
  let hash_values = foreign "hash_values" (ptr hashtable @-> returning gen)

  let hash_search =
    foreign "hash_search"
      (ptr hashtable @-> ptr void @-> returning (ptr hashentry))

  let hash_search2 =
    foreign "hash_search2"
      (ptr hashtable @-> ptr void @-> pari_ulong @-> returning (ptr hashentry))

  let hash_select =
    foreign "hash_select"
      (ptr hashtable @-> ptr void @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> ptr hashentry @-> returning int)
      @-> returning (ptr hashentry))

  let hash_remove =
    foreign "hash_remove"
      (ptr hashtable @-> ptr void @-> returning (ptr hashentry))

  let hash_remove_select =
    foreign "hash_remove_select"
      (ptr hashtable @-> ptr void @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> ptr hashentry @-> returning int)
      @-> returning (ptr hashentry))

  let hash_destroy = foreign "hash_destroy" (ptr hashtable @-> returning void)
  let hash_gen = foreign "hash_GEN" (gen @-> returning pari_ulong)
  let hash_zv = foreign "hash_zv" (gen @-> returning pari_ulong)

  let zx_hyperellred =
    foreign "ZX_hyperellred" (gen @-> ptr gen @-> returning gen)

  let hyperellcharpoly = foreign "hyperellcharpoly" (gen @-> returning gen)

  let hyperellchangecurve =
    foreign "hyperellchangecurve" (gen @-> gen @-> returning gen)

  let hyperelldisc = foreign "hyperelldisc" (gen @-> returning gen)

  let hyperellisoncurve =
    foreign "hyperellisoncurve" (gen @-> gen @-> returning int)

  let hyperellminimaldisc =
    foreign "hyperellminimaldisc" (gen @-> gen @-> returning gen)

  let hyperellminimalmodel =
    foreign "hyperellminimalmodel" (gen @-> ptr gen @-> gen @-> returning gen)

  let hyperellpadicfrobenius0 =
    foreign "hyperellpadicfrobenius0" (gen @-> gen @-> long @-> returning gen)

  let hyperellpadicfrobenius =
    foreign "hyperellpadicfrobenius"
      (gen @-> pari_ulong @-> long @-> returning gen)

  let hyperellred = foreign "hyperellred" (gen @-> ptr gen @-> returning gen)

  let nfhyperellpadicfrobenius =
    foreign "nfhyperellpadicfrobenius"
      (gen @-> gen @-> pari_ulong @-> long @-> returning gen)

  let hypergeom =
    foreign "hypergeom" (gen @-> gen @-> gen @-> long @-> returning gen)

  let airy = foreign "airy" (gen @-> long @-> returning gen)

  let rgm_hnfall =
    foreign "RgM_hnfall" (gen @-> ptr gen @-> long @-> returning gen)

  let zm_hnf = foreign "ZM_hnf" (gen @-> returning gen)
  let zm_hnf_knapsack = foreign "ZM_hnf_knapsack" (gen @-> returning gen)

  let zm_hnfall =
    foreign "ZM_hnfall" (gen @-> ptr gen @-> long @-> returning gen)

  let zm_hnfall_i =
    foreign "ZM_hnfall_i" (gen @-> ptr gen @-> long @-> returning gen)

  let zm_hnfcenter = foreign "ZM_hnfcenter" (gen @-> returning gen)
  let zm_hnflll = foreign "ZM_hnflll" (gen @-> ptr gen @-> int @-> returning gen)
  let zv_extgcd = foreign "ZV_extgcd" (gen @-> returning gen)

  let zv_snfall =
    foreign "ZV_snfall" (gen @-> ptr gen @-> ptr gen @-> returning gen)

  let zv_snf_group =
    foreign "ZV_snf_group" (gen @-> ptr gen @-> ptr gen @-> returning gen)

  let zv_snf_rank = foreign "ZV_snf_rank" (gen @-> gen @-> returning long)

  let zv_snf_rank_u =
    foreign "ZV_snf_rank_u" (gen @-> pari_ulong @-> returning long)

  let zv_snf_trunc = foreign "ZV_snf_trunc" (gen @-> returning void)
  let zm_hnfmod = foreign "ZM_hnfmod" (gen @-> gen @-> returning gen)

  let zm_hnfmodall =
    foreign "ZM_hnfmodall" (gen @-> gen @-> long @-> returning gen)

  let zm_hnfmodall_i =
    foreign "ZM_hnfmodall_i" (gen @-> gen @-> long @-> returning gen)

  let zm_hnfmodid = foreign "ZM_hnfmodid" (gen @-> gen @-> returning gen)
  let zm_hnfmodprime = foreign "ZM_hnfmodprime" (gen @-> gen @-> returning gen)

  let zm_hnfperm =
    foreign "ZM_hnfperm" (gen @-> ptr gen @-> ptr gen @-> returning gen)

  let zm_snfclean =
    foreign "ZM_snfclean" (gen @-> gen @-> gen @-> returning void)

  let zm_snf = foreign "ZM_snf" (gen @-> returning gen)

  let zm_snf_group =
    foreign "ZM_snf_group" (gen @-> ptr gen @-> ptr gen @-> returning gen)

  let zm_snfall =
    foreign "ZM_snfall" (gen @-> ptr gen @-> ptr gen @-> returning gen)

  let zm_snfall_i =
    foreign "ZM_snfall_i"
      (gen @-> ptr gen @-> ptr gen @-> long @-> returning gen)

  let zv_snfclean = foreign "ZV_snfclean" (gen @-> returning gen)

  let zpm_echelon =
    foreign "ZpM_echelon" (gen @-> long @-> gen @-> gen @-> returning gen)

  let gsmith = foreign "gsmith" (gen @-> returning gen)
  let gsmithall = foreign "gsmithall" (gen @-> returning gen)
end

module F38 (F : Ctypes.FOREIGN) = struct
  open F

  let hnf = foreign "hnf" (gen @-> returning gen)

  let hnf_divscale =
    foreign "hnf_divscale" (gen @-> gen @-> gen @-> returning gen)

  let hnf_invscale = foreign "hnf_invscale" (gen @-> gen @-> returning gen)
  let hnf_solve = foreign "hnf_solve" (gen @-> gen @-> returning gen)
  let hnf_invimage = foreign "hnf_invimage" (gen @-> gen @-> returning gen)
  let hnfall = foreign "hnfall" (gen @-> returning gen)
  let hnfdivide = foreign "hnfdivide" (gen @-> gen @-> returning int)
  let hnflll = foreign "hnflll" (gen @-> returning gen)
  let hnfmerge_get_1 = foreign "hnfmerge_get_1" (gen @-> gen @-> returning gen)
  let hnfmod = foreign "hnfmod" (gen @-> gen @-> returning gen)
  let hnfmodid = foreign "hnfmodid" (gen @-> gen @-> returning gen)
  let hnfperm = foreign "hnfperm" (gen @-> returning gen)

  let matfrobenius =
    foreign "matfrobenius" (gen @-> long @-> long @-> returning gen)

  let mathnf0 = foreign "mathnf0" (gen @-> long @-> returning gen)
  let matsnf0 = foreign "matsnf0" (gen @-> long @-> returning gen)
  let smith = foreign "smith" (gen @-> returning gen)
  let smithall = foreign "smithall" (gen @-> returning gen)
  let smithclean = foreign "smithclean" (gen @-> returning gen)
  let snfrank = foreign "snfrank" (gen @-> gen @-> returning long)

  let zlm_echelon =
    foreign "zlm_echelon"
      (gen @-> long @-> pari_ulong @-> pari_ulong @-> returning gen)

  let zv_snf_rank = foreign "zv_snf_rank" (gen @-> pari_ulong @-> returning long)

  let z_ecm =
    foreign "Z_ECM" (gen @-> long @-> long @-> pari_ulong @-> returning gen)

  let z_factor = foreign "Z_factor" (gen @-> returning gen)

  let z_factor_limit =
    foreign "Z_factor_limit" (gen @-> pari_ulong @-> returning gen)

  let z_factor_until = foreign "Z_factor_until" (gen @-> gen @-> returning gen)
  let z_issmooth = foreign "Z_issmooth" (gen @-> pari_ulong @-> returning long)

  let z_issmooth_fact =
    foreign "Z_issmooth_fact" (gen @-> pari_ulong @-> returning gen)

  let z_issquarefree = foreign "Z_issquarefree" (gen @-> returning long)

  let z_pollardbrent =
    foreign "Z_pollardbrent" (gen @-> long @-> long @-> returning gen)

  let absz_factor = foreign "absZ_factor" (gen @-> returning gen)

  let absz_factor_limit =
    foreign "absZ_factor_limit" (gen @-> pari_ulong @-> returning gen)

  let absz_factor_limit_strict =
    foreign "absZ_factor_limit_strict"
      (gen @-> pari_ulong @-> ptr gen @-> returning gen)

  let coreu = foreign "coreu" (pari_ulong @-> returning pari_ulong)
  let coreu_fact = foreign "coreu_fact" (gen @-> returning pari_ulong)
  let factorint = foreign "factorint" (gen @-> long @-> returning gen)
  let factoru = foreign "factoru" (pari_ulong @-> returning gen)

  let tridiv_boundu =
    foreign "tridiv_boundu" (pari_ulong @-> returning pari_ulong)

  let ifac_isprime = foreign "ifac_isprime" (gen @-> returning int)

  let ifac_next =
    foreign "ifac_next" (ptr gen @-> ptr gen @-> ptr long @-> returning int)

  let ifac_read =
    foreign "ifac_read" (gen @-> ptr gen @-> ptr long @-> returning int)

  let ifac_skip = foreign "ifac_skip" (gen @-> returning void)
  let ifac_start = foreign "ifac_start" (gen @-> int @-> returning gen)

  let is_357_power =
    foreign "is_357_power" (gen @-> ptr gen @-> ptr pari_ulong @-> returning int)

  let is_pth_power =
    foreign "is_pth_power"
      (gen @-> ptr gen @-> ptr forprime_t @-> pari_ulong @-> returning int)

  let ispowerful = foreign "ispowerful" (gen @-> returning long)
  let maxomegau = foreign "maxomegau" (pari_ulong @-> returning long)
  let maxomegaoddu = foreign "maxomegaoddu" (pari_ulong @-> returning long)
  let moebius = foreign "moebius" (gen @-> returning long)
  let moebiusu = foreign "moebiusu" (pari_ulong @-> returning long)
  let moebiusu_fact = foreign "moebiusu_fact" (gen @-> returning long)
  let nextprime = foreign "nextprime" (gen @-> returning gen)
  let precprime = foreign "precprime" (gen @-> returning gen)
  let radicalu = foreign "radicalu" (pari_ulong @-> returning pari_ulong)
  let tridiv_bound = foreign "tridiv_bound" (gen @-> returning pari_ulong)

  let uis_357_power =
    foreign "uis_357_power"
      (pari_ulong @-> ptr pari_ulong @-> ptr pari_ulong @-> returning int)

  let uis_357_powermod =
    foreign "uis_357_powermod" (pari_ulong @-> ptr pari_ulong @-> returning int)

  let unextprime = foreign "unextprime" (pari_ulong @-> returning pari_ulong)
  let uprecprime = foreign "uprecprime" (pari_ulong @-> returning pari_ulong)

  let vecfactorsquarefreeu =
    foreign "vecfactorsquarefreeu" (pari_ulong @-> pari_ulong @-> returning gen)

  let vecfactorsquarefreeu_coprime =
    foreign "vecfactorsquarefreeu_coprime"
      (pari_ulong @-> pari_ulong @-> gen @-> returning gen)

  let vecfactoru_i =
    foreign "vecfactoru_i" (pari_ulong @-> pari_ulong @-> returning gen)

  let vecfactoru =
    foreign "vecfactoru" (pari_ulong @-> pari_ulong @-> returning gen)

  let vecfactoroddu_i =
    foreign "vecfactoroddu_i" (pari_ulong @-> pari_ulong @-> returning gen)

  let vecfactoroddu =
    foreign "vecfactoroddu" (pari_ulong @-> pari_ulong @-> returning gen)

  let vecsquarefreeu =
    foreign "vecsquarefreeu" (pari_ulong @-> pari_ulong @-> returning gen)

  let chk_gerepileupto = foreign "chk_gerepileupto" (gen @-> returning int)
  let copy_bin = foreign "copy_bin" (gen @-> returning (ptr genbin))
  let copy_bin_canon = foreign "copy_bin_canon" (gen @-> returning (ptr genbin))
  let dbg_gerepile = foreign "dbg_gerepile" (pari_sp @-> returning void)
  let dbg_gerepileupto = foreign "dbg_gerepileupto" (gen @-> returning void)
  let errname = foreign "errname" (gen @-> returning gen)
  let gclone = foreign "gclone" (gen @-> returning gen)
  let gcloneref = foreign "gcloneref" (gen @-> returning gen)
  let gclone_refc = foreign "gclone_refc" (gen @-> returning void)
  let gcopy = foreign "gcopy" (gen @-> returning gen)
  let gcopy_avma = foreign "gcopy_avma" (gen @-> ptr pari_sp @-> returning gen)
  let gcopy_lg = foreign "gcopy_lg" (gen @-> long @-> returning gen)

  let gerepile =
    foreign "gerepile" (pari_sp @-> pari_sp @-> gen @-> returning gen)

  let gerepileallsp =
    foreign "gerepileallsp" (pari_sp @-> pari_sp @-> int @-> returning void)

  let gerepilecoeffssp =
    foreign "gerepilecoeffssp"
      (pari_sp @-> pari_sp @-> ptr long @-> int @-> returning void)

  let gerepilemanysp =
    foreign "gerepilemanysp"
      (pari_sp @-> pari_sp @-> ptr (ptr gen) @-> int @-> returning void)

  let getheap = foreign "getheap" (void @-> returning gen)

  let gp_context_save =
    foreign "gp_context_save" (ptr gp_context @-> returning void)

  let gp_context_restore =
    foreign "gp_context_restore" (ptr gp_context @-> returning void)

  let gsizeword = foreign "gsizeword" (gen @-> returning long)
  let gsizebyte = foreign "gsizebyte" (gen @-> returning long)
  let gunclone = foreign "gunclone" (gen @-> returning void)
  let gunclone_deep = foreign "gunclone_deep" (gen @-> returning void)
  let listcopy = foreign "listcopy" (gen @-> returning gen)
  let listinit = foreign "listinit" (gen @-> returning gen)
  let msgtimer = foreign "msgtimer" (string @-> returning void)
  let name_numerr = foreign "name_numerr" (string @-> returning long)
  let new_chunk_resize = foreign "new_chunk_resize" (int @-> returning void)
  let newblock = foreign "newblock" (int @-> returning gen)
  let numerr_name = foreign "numerr_name" (long @-> returning string)
  let obj_check = foreign "obj_check" (gen @-> long @-> returning gen)

  let obj_checkbuild =
    foreign "obj_checkbuild"
      (gen @-> long
      @-> static_funptr Ctypes.(gen @-> returning gen)
      @-> returning gen)

  let obj_checkbuild_padicprec =
    foreign "obj_checkbuild_padicprec"
      (gen @-> long
      @-> static_funptr Ctypes.(gen @-> long @-> returning gen)
      @-> long @-> returning gen)

  let obj_checkbuild_realprec =
    foreign "obj_checkbuild_realprec"
      (gen @-> long
      @-> static_funptr Ctypes.(gen @-> long @-> returning gen)
      @-> long @-> returning gen)

  let obj_checkbuild_prec =
    foreign "obj_checkbuild_prec"
      (gen @-> long
      @-> static_funptr Ctypes.(gen @-> long @-> returning gen)
      @-> static_funptr Ctypes.(gen @-> returning long)
      @-> long @-> returning gen)
end

module F39 (F : Ctypes.FOREIGN) = struct
  open F

  let obj_free = foreign "obj_free" (gen @-> returning void)
  let obj_init = foreign "obj_init" (long @-> long @-> returning gen)
  let obj_insert = foreign "obj_insert" (gen @-> long @-> gen @-> returning gen)

  let obj_insert_shallow =
    foreign "obj_insert_shallow" (gen @-> long @-> gen @-> returning gen)

  let obj_reinit = foreign "obj_reinit" (gen @-> returning gen)

  let pari_add_function =
    foreign "pari_add_function" (ptr entree @-> returning void)

  let pari_add_module = foreign "pari_add_module" (ptr entree @-> returning void)

  let pari_add_defaults_module =
    foreign "pari_add_defaults_module" (ptr entree @-> returning void)

  let pari_close = foreign "pari_close" (void @-> returning void)
  let pari_close_opts = foreign "pari_close_opts" (pari_ulong @-> returning void)
  let pari_compile_str = foreign "pari_compile_str" (string @-> returning gen)
  let pari_daemon = foreign "pari_daemon" (void @-> returning int)
  let pari_err = foreign "pari_err" (int @-> returning void)
  let pari_err_last = foreign "pari_err_last" (void @-> returning gen)
  let pari_err2str = foreign "pari_err2str" (gen @-> returning string)

  let pari_init_opts =
    foreign "pari_init_opts"
      (int @-> pari_ulong @-> pari_ulong @-> returning void)

  let pari_init = foreign "pari_init" (int @-> pari_ulong @-> returning void)
  let pari_sighandler = foreign "pari_sighandler" (int @-> returning void)

  let pari_sig_init =
    foreign "pari_sig_init"
      (static_funptr Ctypes.(int @-> returning void) @-> returning void)

  let pari_thread_alloc =
    foreign "pari_thread_alloc"
      (ptr pari_thread @-> int @-> gen @-> returning void)

  let pari_thread_close = foreign "pari_thread_close" (void @-> returning void)

  let pari_thread_free =
    foreign "pari_thread_free" (ptr pari_thread @-> returning void)

  let pari_thread_init = foreign "pari_thread_init" (void @-> returning void)

  let pari_thread_start =
    foreign "pari_thread_start" (ptr pari_thread @-> returning gen)

  let pari_thread_valloc =
    foreign "pari_thread_valloc"
      (ptr pari_thread @-> int @-> int @-> gen @-> returning void)

  let pari_version = foreign "pari_version" (void @-> returning gen)
  let pari_warn = foreign "pari_warn" (int @-> returning void)

  let paristack_newrsize =
    foreign "paristack_newrsize" (pari_ulong @-> returning void)

  let paristack_resize =
    foreign "paristack_resize" (pari_ulong @-> returning void)

  let paristack_setsize =
    foreign "paristack_setsize" (int @-> int @-> returning void)

  let parivstack_resize =
    foreign "parivstack_resize" (pari_ulong @-> returning void)

  let parivstack_reset = foreign "parivstack_reset" (void @-> returning void)
  let setalldebug = foreign "setalldebug" (long @-> returning void)
  let setdebug = foreign "setdebug" (string @-> long @-> returning gen)
  let shiftaddress = foreign "shiftaddress" (gen @-> long @-> returning void)

  let shiftaddress_canon =
    foreign "shiftaddress_canon" (gen @-> long @-> returning void)

  let timer = foreign "timer" (void @-> returning long)
  let timer_delay = foreign "timer_delay" (ptr pari_timer @-> returning long)
  let timer_get = foreign "timer_get" (ptr pari_timer @-> returning long)

  let timer_printf =
    foreign "timer_printf" (ptr pari_timer @-> string @-> returning void)

  let timer_start = foreign "timer_start" (ptr pari_timer @-> returning void)
  let timer2 = foreign "timer2" (void @-> returning long)
  let trap0 = foreign "trap0" (string @-> gen @-> gen @-> returning gen)

  let traverseheap =
    foreign "traverseheap"
      (static_funptr Ctypes.(gen @-> ptr void @-> returning void)
      @-> ptr void @-> returning void)

  let walltimer_start =
    foreign "walltimer_start" (ptr pari_timer @-> returning void)

  let walltimer_delay =
    foreign "walltimer_delay" (ptr pari_timer @-> returning long)

  let walltimer_get = foreign "walltimer_get" (ptr pari_timer @-> returning long)

  let contfraceval =
    foreign "contfraceval" (gen @-> gen @-> long @-> returning gen)

  let contfracinit = foreign "contfracinit" (gen @-> long @-> returning gen)

  let intcirc =
    foreign "intcirc"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let intfuncinit =
    foreign "intfuncinit"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> long @-> returning gen)

  let intnum =
    foreign "intnum"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let intnumgauss =
    foreign "intnumgauss"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let intnumgaussinit =
    foreign "intnumgaussinit" (long @-> long @-> returning gen)

  let intnuminit =
    foreign "intnuminit" (gen @-> gen @-> long @-> long @-> returning gen)

  let intnumosc =
    foreign "intnumosc"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> gen @-> long @-> returning gen)

  let intnumromb =
    foreign "intnumromb"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> long @-> returning gen)

  let intnumromb_bitprec =
    foreign "intnumromb_bitprec"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> long @-> returning gen)

  let prodeulerrat =
    foreign "prodeulerrat" (gen @-> gen @-> long @-> long @-> returning gen)

  let prodnumrat = foreign "prodnumrat" (gen @-> long @-> long @-> returning gen)
  let quodif = foreign "quodif" (gen @-> long @-> returning gen)

  let sumeulerrat =
    foreign "sumeulerrat" (gen @-> gen @-> long @-> long @-> returning gen)

  let sumnum =
    foreign "sumnum"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> returning gen)

  let sumnumap =
    foreign "sumnumap"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> returning gen)

  let sumnumapinit = foreign "sumnumapinit" (gen @-> long @-> returning gen)
  let sumnuminit = foreign "sumnuminit" (gen @-> long @-> returning gen)

  let sumnumlagrangeinit =
    foreign "sumnumlagrangeinit" (gen @-> gen @-> long @-> returning gen)

  let sumnumlagrange =
    foreign "sumnumlagrange"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> gen @-> gen @-> long @-> returning gen)

  let sumnummonien =
    foreign "sumnummonien"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> returning gen)

  let sumnummonieninit =
    foreign "sumnummonieninit" (gen @-> gen @-> gen @-> long @-> returning gen)

  let sumnumrat = foreign "sumnumrat" (gen @-> gen @-> long @-> returning gen)

  let sumnumsidi =
    foreign "sumnumsidi"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> gen @-> double @-> long @-> returning gen)

  let z_isanypower = foreign "Z_isanypower" (gen @-> ptr gen @-> returning long)
  let z_ispow2 = foreign "Z_ispow2" (gen @-> returning long)

  let z_ispowerall =
    foreign "Z_ispowerall" (gen @-> pari_ulong @-> ptr gen @-> returning long)

  let z_issquareall =
    foreign "Z_issquareall" (gen @-> ptr gen @-> returning long)

  let zn_ispower =
    foreign "Zn_ispower" (gen @-> gen @-> gen @-> ptr gen @-> returning long)

  let zn_issquare = foreign "Zn_issquare" (gen @-> gen @-> returning long)
  let zp_issquare = foreign "Zp_issquare" (gen @-> gen @-> returning long)
  let gisanypower = foreign "gisanypower" (gen @-> ptr gen @-> returning long)
  let gissquare = foreign "gissquare" (gen @-> returning gen)
  let gissquareall = foreign "gissquareall" (gen @-> ptr gen @-> returning gen)

  let ispolygonal =
    foreign "ispolygonal" (gen @-> gen @-> ptr gen @-> returning long)

  let ispower = foreign "ispower" (gen @-> gen @-> ptr gen @-> returning long)
  let isprimepower = foreign "isprimepower" (gen @-> ptr gen @-> returning long)

  let ispseudoprimepower =
    foreign "ispseudoprimepower" (gen @-> ptr gen @-> returning long)

  let issquare = foreign "issquare" (gen @-> returning long)
  let issquareall = foreign "issquareall" (gen @-> ptr gen @-> returning long)
  let sqrtint = foreign "sqrtint" (gen @-> returning gen)
  let sqrtint0 = foreign "sqrtint0" (gen @-> ptr gen @-> returning gen)

  let uisprimepower =
    foreign "uisprimepower" (pari_ulong @-> ptr pari_ulong @-> returning long)

  let uissquare = foreign "uissquare" (pari_ulong @-> returning long)

  let uissquareall =
    foreign "uissquareall" (pari_ulong @-> ptr pari_ulong @-> returning long)

  let ulogintall =
    foreign "ulogintall"
      (pari_ulong @-> pari_ulong @-> ptr pari_ulong @-> returning long)

  let padicfields0 =
    foreign "padicfields0" (gen @-> gen @-> long @-> returning gen)

  let padicfields =
    foreign "padicfields" (gen @-> long @-> long @-> long @-> returning gen)

  let bnrclassfield =
    foreign "bnrclassfield" (gen @-> gen @-> long @-> long @-> returning gen)

  let rnfkummer = foreign "rnfkummer" (gen @-> gen @-> long @-> returning gen)
  let is_linit = foreign "is_linit" (gen @-> returning long)
  let ldata_get_an = foreign "ldata_get_an" (gen @-> returning gen)
end

module F40 (F : Ctypes.FOREIGN) = struct
  open F

  let ldata_get_dual = foreign "ldata_get_dual" (gen @-> returning gen)
  let ldata_get_gammavec = foreign "ldata_get_gammavec" (gen @-> returning gen)
  let ldata_get_degree = foreign "ldata_get_degree" (gen @-> returning long)
  let ldata_get_k = foreign "ldata_get_k" (gen @-> returning gen)
  let ldata_get_k1 = foreign "ldata_get_k1" (gen @-> returning gen)
  let ldata_get_conductor = foreign "ldata_get_conductor" (gen @-> returning gen)
  let ldata_get_rootno = foreign "ldata_get_rootno" (gen @-> returning gen)
  let ldata_get_residue = foreign "ldata_get_residue" (gen @-> returning gen)
  let ldata_get_type = foreign "ldata_get_type" (gen @-> returning long)
  let ldata_isreal = foreign "ldata_isreal" (gen @-> returning long)
  let linit_get_type = foreign "linit_get_type" (gen @-> returning long)
  let linit_get_ldata = foreign "linit_get_ldata" (gen @-> returning gen)
  let linit_get_tech = foreign "linit_get_tech" (gen @-> returning gen)
  let lfun_get_domain = foreign "lfun_get_domain" (gen @-> returning gen)
  let lfun_get_dom = foreign "lfun_get_dom" (gen @-> returning gen)
  let lfun_get_bitprec = foreign "lfun_get_bitprec" (gen @-> returning long)

  let lfun_get_factgammavec =
    foreign "lfun_get_factgammavec" (gen @-> returning gen)

  let lfun_get_step = foreign "lfun_get_step" (gen @-> returning gen)
  let lfun_get_pol = foreign "lfun_get_pol" (gen @-> returning gen)
  let lfun_get_residue = foreign "lfun_get_Residue" (gen @-> returning gen)
  let lfun_get_k2 = foreign "lfun_get_k2" (gen @-> returning gen)
  let lfun_get_w2 = foreign "lfun_get_w2" (gen @-> returning gen)
  let lfun_get_expot = foreign "lfun_get_expot" (gen @-> returning gen)
  let lfun_get_bitprec = foreign "lfun_get_bitprec" (gen @-> returning long)
  let lfun = foreign "lfun" (gen @-> gen @-> long @-> returning gen)
  let lfun0 = foreign "lfun0" (gen @-> gen @-> long @-> long @-> returning gen)
  let lfunan = foreign "lfunan" (gen @-> long @-> long @-> returning gen)

  let lfuncheckfeq =
    foreign "lfuncheckfeq" (gen @-> gen @-> long @-> returning long)

  let lfunconductor =
    foreign "lfunconductor" (gen @-> gen @-> long @-> long @-> returning gen)

  let lfuncost =
    foreign "lfuncost" (gen @-> gen @-> long @-> long @-> returning gen)

  let lfuncost0 =
    foreign "lfuncost0" (gen @-> gen @-> long @-> long @-> returning gen)

  let lfuncreate = foreign "lfuncreate" (gen @-> returning gen)
  let lfundual = foreign "lfundual" (gen @-> long @-> returning gen)
  let lfuneuler = foreign "lfuneuler" (gen @-> gen @-> long @-> returning gen)
  let lfunparams = foreign "lfunparams" (gen @-> long @-> returning gen)
  let lfunan = foreign "lfunan" (gen @-> long @-> long @-> returning gen)
  let lfunhardy = foreign "lfunhardy" (gen @-> gen @-> long @-> returning gen)

  let lfuninit =
    foreign "lfuninit" (gen @-> gen @-> long @-> long @-> returning gen)

  let lfuninit0 =
    foreign "lfuninit0" (gen @-> gen @-> long @-> long @-> returning gen)

  let lfuninit_make =
    foreign "lfuninit_make" (long @-> gen @-> gen @-> gen @-> returning gen)

  let lfunlambda = foreign "lfunlambda" (gen @-> gen @-> long @-> returning gen)

  let lfunlambda0 =
    foreign "lfunlambda0" (gen @-> gen @-> long @-> long @-> returning gen)

  let lfunmisc_to_ldata = foreign "lfunmisc_to_ldata" (gen @-> returning gen)

  let lfunmisc_to_ldata_shallow =
    foreign "lfunmisc_to_ldata_shallow" (gen @-> returning gen)

  let lfunmisc_to_ldata_shallow_i =
    foreign "lfunmisc_to_ldata_shallow_i" (gen @-> returning gen)

  let lfunorderzero =
    foreign "lfunorderzero" (gen @-> long @-> long @-> returning long)

  let lfunprod_get_fact = foreign "lfunprod_get_fact" (gen @-> returning gen)
  let lfunrootno = foreign "lfunrootno" (gen @-> long @-> returning gen)
  let lfunrootres = foreign "lfunrootres" (gen @-> long @-> returning gen)
  let lfunrtopoles = foreign "lfunrtopoles" (gen @-> returning gen)

  let lfunshift =
    foreign "lfunshift" (gen @-> gen @-> long @-> long @-> returning gen)

  let lfuntwist = foreign "lfuntwist" (gen @-> gen @-> long @-> returning gen)

  let lfuntheta =
    foreign "lfuntheta" (gen @-> gen @-> long @-> long @-> returning gen)

  let lfunthetacost0 =
    foreign "lfunthetacost0" (gen @-> gen @-> long @-> long @-> returning long)

  let lfunthetacost =
    foreign "lfunthetacost" (gen @-> gen @-> long @-> long @-> returning long)

  let lfunthetainit =
    foreign "lfunthetainit" (gen @-> gen @-> long @-> long @-> returning gen)

  let lfunthetacheckinit =
    foreign "lfunthetacheckinit"
      (gen @-> gen @-> long @-> long @-> returning gen)

  let lfunzeros =
    foreign "lfunzeros" (gen @-> gen @-> long @-> long @-> returning gen)

  let sdomain_isincl =
    foreign "sdomain_isincl" (double @-> gen @-> gen @-> returning int)

  let theta_get_an = foreign "theta_get_an" (gen @-> returning gen)
  let theta_get_k = foreign "theta_get_K" (gen @-> returning gen)
  let theta_get_r = foreign "theta_get_R" (gen @-> returning gen)
  let theta_get_bitprec = foreign "theta_get_bitprec" (gen @-> returning long)
  let theta_get_m = foreign "theta_get_m" (gen @-> returning long)
  let theta_get_tdom = foreign "theta_get_tdom" (gen @-> returning gen)
  let theta_get_isqrtn = foreign "theta_get_isqrtN" (gen @-> returning gen)
  let vgaeasytheta = foreign "Vgaeasytheta" (gen @-> returning int)

  let znchargauss =
    foreign "znchargauss" (gen @-> gen @-> gen @-> long @-> returning gen)

  let dirzetak = foreign "dirzetak" (gen @-> gen @-> returning gen)
  let ellmoddegree = foreign "ellmoddegree" (gen @-> returning gen)
  let eta_zxn = foreign "eta_ZXn" (long @-> long @-> returning gen)

  let eta_product_zxn =
    foreign "eta_product_ZXn" (gen @-> long @-> returning gen)

  let etaquotype =
    foreign "etaquotype"
      (ptr gen @-> ptr gen @-> ptr gen @-> ptr gen @-> ptr long @-> ptr long
     @-> ptr long @-> returning long)

  let galois_get_conj = foreign "galois_get_conj" (gen @-> returning gen)

  let ldata_vecan =
    foreign "ldata_vecan" (gen @-> long @-> long @-> returning gen)

  let ldata_newprec = foreign "ldata_newprec" (gen @-> long @-> returning gen)

  let lfunabelianrelinit =
    foreign "lfunabelianrelinit"
      (gen @-> gen @-> gen @-> long @-> long @-> returning gen)

  let lfunartin =
    foreign "lfunartin" (gen @-> gen @-> gen @-> long @-> long @-> returning gen)

  let lfundiv = foreign "lfundiv" (gen @-> gen @-> long @-> returning gen)

  let lfunellmfpeters =
    foreign "lfunellmfpeters" (gen @-> long @-> returning gen)

  let lfunetaquo = foreign "lfunetaquo" (gen @-> returning gen)
  let lfungenus2 = foreign "lfungenus2" (gen @-> returning gen)
  let lfunmfspec = foreign "lfunmfspec" (gen @-> long @-> returning gen)
  let lfunmul = foreign "lfunmul" (gen @-> gen @-> long @-> returning gen)
  let lfunqf = foreign "lfunqf" (gen @-> long @-> returning gen)
  let lfunsympow = foreign "lfunsympow" (gen @-> pari_ulong @-> returning gen)

  let lfunzetakinit =
    foreign "lfunzetakinit" (gen @-> gen @-> long @-> long @-> returning gen)

  let qfiseven = foreign "qfiseven" (gen @-> returning long)
  let lfunquadneg = foreign "lfunquadneg" (long @-> long @-> returning gen)

  let zm_lll_norms =
    foreign "ZM_lll_norms"
      (gen @-> double @-> long @-> ptr gen @-> returning gen)

  let kerint = foreign "kerint" (gen @-> returning gen)
  let lll = foreign "lll" (gen @-> returning gen)
  let lllfp = foreign "lllfp" (gen @-> double @-> long @-> returning gen)
  let lllgen = foreign "lllgen" (gen @-> returning gen)
  let lllgram = foreign "lllgram" (gen @-> returning gen)
  let lllgramgen = foreign "lllgramgen" (gen @-> returning gen)
  let lllgramint = foreign "lllgramint" (gen @-> returning gen)
  let lllgramkerim = foreign "lllgramkerim" (gen @-> returning gen)
  let lllgramkerimgen = foreign "lllgramkerimgen" (gen @-> returning gen)
  let lllint = foreign "lllint" (gen @-> returning gen)
end

module F41 (F : Ctypes.FOREIGN) = struct
  open F

  let lllintpartial = foreign "lllintpartial" (gen @-> returning gen)

  let lllintpartial_inplace =
    foreign "lllintpartial_inplace" (gen @-> returning gen)

  let lllkerim = foreign "lllkerim" (gen @-> returning gen)
  let lllkerimgen = foreign "lllkerimgen" (gen @-> returning gen)
  let matkerint0 = foreign "matkerint0" (gen @-> long @-> returning gen)
  let qflll0 = foreign "qflll0" (gen @-> long @-> returning gen)
  let qflllgram0 = foreign "qflllgram0" (gen @-> long @-> returning gen)
  let gtomap = foreign "gtomap" (gen @-> returning gen)
  let mapdelete = foreign "mapdelete" (gen @-> gen @-> returning void)
  let mapdomain = foreign "mapdomain" (gen @-> returning gen)
  let mapdomain_shallow = foreign "mapdomain_shallow" (gen @-> returning gen)
  let mapget = foreign "mapget" (gen @-> gen @-> returning gen)

  let mapisdefined =
    foreign "mapisdefined" (gen @-> gen @-> ptr gen @-> returning int)

  let mapput = foreign "mapput" (gen @-> gen @-> gen @-> returning void)
  let maptomat = foreign "maptomat" (gen @-> returning gen)
  let maptomat_shallow = foreign "maptomat_shallow" (gen @-> returning gen)
  let matpermanent = foreign "matpermanent" (gen @-> returning gen)
  let zm_permanent = foreign "zm_permanent" (gen @-> returning gen)
  let zm_permanent = foreign "ZM_permanent" (gen @-> returning gen)

  let dbllemma526 =
    foreign "dbllemma526"
      (double @-> double @-> double @-> double @-> returning double)

  let dblcoro526 =
    foreign "dblcoro526" (double @-> double @-> double @-> returning double)

  let gammamellininv =
    foreign "gammamellininv" (gen @-> gen @-> long @-> long @-> returning gen)

  let gammamellininvasymp =
    foreign "gammamellininvasymp" (gen @-> long @-> long @-> returning gen)

  let gammamellininvinit =
    foreign "gammamellininvinit" (gen @-> long @-> long @-> returning gen)

  let gammamellininvrt =
    foreign "gammamellininvrt" (gen @-> gen @-> long @-> returning gen)

  let member_a1 = foreign "member_a1" (gen @-> returning gen)
  let member_a2 = foreign "member_a2" (gen @-> returning gen)
  let member_a3 = foreign "member_a3" (gen @-> returning gen)
  let member_a4 = foreign "member_a4" (gen @-> returning gen)
  let member_a6 = foreign "member_a6" (gen @-> returning gen)
  let member_area = foreign "member_area" (gen @-> returning gen)
  let member_b2 = foreign "member_b2" (gen @-> returning gen)
  let member_b4 = foreign "member_b4" (gen @-> returning gen)
  let member_b6 = foreign "member_b6" (gen @-> returning gen)
  let member_b8 = foreign "member_b8" (gen @-> returning gen)
  let member_bid = foreign "member_bid" (gen @-> returning gen)
  let member_bnf = foreign "member_bnf" (gen @-> returning gen)
  let member_c4 = foreign "member_c4" (gen @-> returning gen)
  let member_c6 = foreign "member_c6" (gen @-> returning gen)
  let member_clgp = foreign "member_clgp" (gen @-> returning gen)
  let member_codiff = foreign "member_codiff" (gen @-> returning gen)
  let member_cyc = foreign "member_cyc" (gen @-> returning gen)
  let member_diff = foreign "member_diff" (gen @-> returning gen)
  let member_disc = foreign "member_disc" (gen @-> returning gen)
  let member_e = foreign "member_e" (gen @-> returning gen)
  let member_eta = foreign "member_eta" (gen @-> returning gen)
  let member_f = foreign "member_f" (gen @-> returning gen)
  let member_fu = foreign "member_fu" (gen @-> returning gen)
  let member_gen = foreign "member_gen" (gen @-> returning gen)
  let member_group = foreign "member_group" (gen @-> returning gen)
  let member_index = foreign "member_index" (gen @-> returning gen)
  let member_j = foreign "member_j" (gen @-> returning gen)
  let member_mod = foreign "member_mod" (gen @-> returning gen)
  let member_nf = foreign "member_nf" (gen @-> returning gen)
  let member_no = foreign "member_no" (gen @-> returning gen)
  let member_omega = foreign "member_omega" (gen @-> returning gen)
  let member_orders = foreign "member_orders" (gen @-> returning gen)
  let member_p = foreign "member_p" (gen @-> returning gen)
  let member_pol = foreign "member_pol" (gen @-> returning gen)
  let member_polabs = foreign "member_polabs" (gen @-> returning gen)
  let member_reg = foreign "member_reg" (gen @-> returning gen)
  let member_r1 = foreign "member_r1" (gen @-> returning gen)
  let member_r2 = foreign "member_r2" (gen @-> returning gen)
  let member_roots = foreign "member_roots" (gen @-> returning gen)
  let member_sign = foreign "member_sign" (gen @-> returning gen)
  let member_t2 = foreign "member_t2" (gen @-> returning gen)
  let member_tate = foreign "member_tate" (gen @-> returning gen)
  let member_tu = foreign "member_tu" (gen @-> returning gen)
  let member_zk = foreign "member_zk" (gen @-> returning gen)
  let member_zkst = foreign "member_zkst" (gen @-> returning gen)
  let mf_get_chi = foreign "MF_get_CHI" (gen @-> returning gen)
  let mf_get_m = foreign "MF_get_M" (gen @-> returning gen)
  let mf_get_mindex = foreign "MF_get_Mindex" (gen @-> returning gen)
  let mf_get_minv = foreign "MF_get_Minv" (gen @-> returning gen)
  let mf_get_n = foreign "MF_get_N" (gen @-> returning long)
  let mf_get_basis = foreign "MF_get_basis" (gen @-> returning gen)
  let mf_get_dim = foreign "MF_get_dim" (gen @-> returning long)
  let mf_get_e = foreign "MF_get_E" (gen @-> returning gen)
  let mf_get_fields = foreign "MF_get_fields" (gen @-> returning gen)
  let mf_get_gn = foreign "MF_get_gN" (gen @-> returning gen)
  let mf_get_gk = foreign "MF_get_gk" (gen @-> returning gen)
  let mf_get_k = foreign "MF_get_k" (gen @-> returning long)
  let mf_get_newforms = foreign "MF_get_newforms" (gen @-> returning gen)
  let mf_get_r = foreign "MF_get_r" (gen @-> returning long)
  let mf_get_space = foreign "MF_get_space" (gen @-> returning long)
  let mf_get_s = foreign "MF_get_S" (gen @-> returning gen)
  let mfcusp_get_vmjd = foreign "MFcusp_get_vMjd" (gen @-> returning gen)
  let mfnew_get_vj = foreign "MFnew_get_vj" (gen @-> returning gen)

  let qab_tracerel =
    foreign "Qab_tracerel" (gen @-> long @-> gen @-> returning gen)

  let qabm_tracerel =
    foreign "QabM_tracerel" (gen @-> long @-> gen @-> returning gen)

  let qabv_tracerel =
    foreign "QabV_tracerel" (gen @-> long @-> gen @-> returning gen)

  let qab_trace_init =
    foreign "Qab_trace_init" (long @-> long @-> gen @-> gen @-> returning gen)

  let checkmf = foreign "checkMF" (gen @-> returning gen)
  let checkmf_i = foreign "checkMF_i" (gen @-> returning gen)
  let checkmf_i = foreign "checkmf_i" (gen @-> returning int)
  let getcache = foreign "getcache" (void @-> returning gen)
  let hclassno6u = foreign "hclassno6u" (pari_ulong @-> returning pari_ulong)

  let hclassno6u_no_cache =
    foreign "hclassno6u_no_cache" (pari_ulong @-> returning pari_ulong)

  let lfunmf = foreign "lfunmf" (gen @-> gen @-> long @-> returning gen)
  let mfdelta = foreign "mfDelta" (void @-> returning gen)
end

module F42 (F : Ctypes.FOREIGN) = struct
  open F

  let mfeh = foreign "mfEH" (gen @-> returning gen)
  let mfek = foreign "mfEk" (long @-> returning gen)
  let mftheta = foreign "mfTheta" (gen @-> returning gen)
  let mf_get_chi = foreign "mf_get_CHI" (gen @-> returning gen)
  let mf_get_n = foreign "mf_get_N" (gen @-> returning long)
  let mf_get_nk = foreign "mf_get_NK" (gen @-> returning gen)
  let mf_get_field = foreign "mf_get_field" (gen @-> returning gen)
  let mf_get_gn = foreign "mf_get_gN" (gen @-> returning gen)
  let mf_get_gk = foreign "mf_get_gk" (gen @-> returning gen)
  let mf_get_k = foreign "mf_get_k" (gen @-> returning long)
  let mf_get_r = foreign "mf_get_r" (gen @-> returning long)
  let mf_get_type = foreign "mf_get_type" (gen @-> returning long)
  let mfatkin = foreign "mfatkin" (gen @-> gen @-> returning gen)

  let mfatkineigenvalues =
    foreign "mfatkineigenvalues" (gen @-> long @-> long @-> returning gen)

  let mfatkininit =
    foreign "mfatkininit" (gen @-> long @-> long @-> returning gen)

  let mfbasis = foreign "mfbasis" (gen @-> long @-> returning gen)
  let mfbd = foreign "mfbd" (gen @-> long @-> returning gen)
  let mfbracket = foreign "mfbracket" (gen @-> gen @-> long @-> returning gen)
  let mfcharorder = foreign "mfcharorder" (gen @-> returning long)
  let mfcharmodulus = foreign "mfcharmodulus" (gen @-> returning long)
  let mfcharpol = foreign "mfcharpol" (gen @-> returning gen)
  let mfcoef = foreign "mfcoef" (gen @-> long @-> returning gen)
  let mfcoefs = foreign "mfcoefs" (gen @-> long @-> long @-> returning gen)
  let mfconductor = foreign "mfconductor" (gen @-> gen @-> returning long)
  let mfcosets = foreign "mfcosets" (gen @-> returning gen)
  let mfcuspdim = foreign "mfcuspdim" (long @-> long @-> gen @-> returning long)

  let mfcuspisregular =
    foreign "mfcuspisregular" (gen @-> gen @-> returning long)

  let mfcusps = foreign "mfcusps" (gen @-> returning gen)

  let mfcuspval =
    foreign "mfcuspval" (gen @-> gen @-> gen @-> long @-> returning gen)

  let mfcuspwidth = foreign "mfcuspwidth" (gen @-> gen @-> returning long)
  let mfderiv = foreign "mfderiv" (gen @-> long @-> returning gen)
  let mfderive2 = foreign "mfderivE2" (gen @-> long @-> returning gen)
  let mfdescribe = foreign "mfdescribe" (gen @-> ptr gen @-> returning gen)
  let mfdim = foreign "mfdim" (gen @-> long @-> returning gen)
  let mfdiv = foreign "mfdiv" (gen @-> gen @-> returning gen)
  let mfdiv_val = foreign "mfdiv_val" (gen @-> gen @-> long @-> returning gen)
  let mfeigenbasis = foreign "mfeigenbasis" (gen @-> returning gen)
  let mfeigensearch = foreign "mfeigensearch" (gen @-> gen @-> returning gen)

  let mfeisenstein =
    foreign "mfeisenstein" (long @-> gen @-> gen @-> returning gen)

  let mfeisensteindim =
    foreign "mfeisensteindim" (long @-> long @-> gen @-> returning long)

  let mfembed = foreign "mfembed" (gen @-> gen @-> returning gen)
  let mfembed0 = foreign "mfembed0" (gen @-> gen @-> long @-> returning gen)
  let mfeval = foreign "mfeval" (gen @-> gen @-> gen @-> long @-> returning gen)
  let mffields = foreign "mffields" (gen @-> returning gen)
  let mffromell = foreign "mffromell" (gen @-> returning gen)
  let mffrometaquo = foreign "mffrometaquo" (gen @-> long @-> returning gen)
  let mffromlfun = foreign "mffromlfun" (gen @-> long @-> returning gen)
  let mffromqf = foreign "mffromqf" (gen @-> gen @-> returning gen)
  let mffulldim = foreign "mffulldim" (long @-> long @-> gen @-> returning long)

  let mfgaloisprojrep =
    foreign "mfgaloisprojrep" (gen @-> gen @-> long @-> returning gen)

  let mfgaloistype = foreign "mfgaloistype" (gen @-> gen @-> returning gen)
  let mfhecke = foreign "mfhecke" (gen @-> gen @-> long @-> returning gen)
  let mfheckemat = foreign "mfheckemat" (gen @-> gen @-> returning gen)
  let mfinit = foreign "mfinit" (gen @-> long @-> returning gen)
  let mfiscm = foreign "mfisCM" (gen @-> returning gen)
  let mfiscuspidal = foreign "mfiscuspidal" (gen @-> gen @-> returning long)
  let mfisequal = foreign "mfisequal" (gen @-> gen @-> long @-> returning long)
  let mfisetaquo = foreign "mfisetaquo" (gen @-> long @-> returning gen)
  let mfkohnenbasis = foreign "mfkohnenbasis" (gen @-> returning gen)
  let mfkohnenbijection = foreign "mfkohnenbijection" (gen @-> returning gen)

  let mfkohneneigenbasis =
    foreign "mfkohneneigenbasis" (gen @-> gen @-> returning gen)

  let mflinear = foreign "mflinear" (gen @-> gen @-> returning gen)
  let mfmanin = foreign "mfmanin" (gen @-> long @-> returning gen)
  let mfmatembed = foreign "mfmatembed" (gen @-> gen @-> returning gen)
  let mfmul = foreign "mfmul" (gen @-> gen @-> returning gen)
  let mfnewdim = foreign "mfnewdim" (long @-> long @-> gen @-> returning long)
  let mfolddim = foreign "mfolddim" (long @-> long @-> gen @-> returning long)
  let mfparams = foreign "mfparams" (gen @-> returning gen)

  let mfperiodpol =
    foreign "mfperiodpol" (gen @-> gen @-> long @-> long @-> returning gen)

  let mfperiodpolbasis =
    foreign "mfperiodpolbasis" (long @-> long @-> returning gen)

  let mfpetersson = foreign "mfpetersson" (gen @-> gen @-> returning gen)
  let mfpow = foreign "mfpow" (gen @-> long @-> returning gen)
  let mfsearch = foreign "mfsearch" (gen @-> gen @-> long @-> returning gen)
  let mfshift = foreign "mfshift" (gen @-> long @-> returning gen)
  let mfshimura = foreign "mfshimura" (gen @-> gen @-> long @-> returning gen)

  let mfslashexpansion =
    foreign "mfslashexpansion"
      (gen @-> gen @-> gen @-> long @-> long @-> ptr gen @-> long
     @-> returning gen)

  let mfspace = foreign "mfspace" (gen @-> gen @-> returning long)
  let mfsplit = foreign "mfsplit" (gen @-> long @-> long @-> returning gen)
  let mfsturm = foreign "mfsturm" (gen @-> returning long)
  let mfsturmngk = foreign "mfsturmNgk" (long @-> gen @-> returning long)
  let mfsturmnk = foreign "mfsturmNk" (long @-> long @-> returning long)
  let mfsturm_mf = foreign "mfsturm_mf" (gen @-> returning long)

  let mfsymboleval =
    foreign "mfsymboleval" (gen @-> gen @-> gen @-> long @-> returning gen)

  let mfsymbol = foreign "mfsymbol" (gen @-> gen @-> long @-> returning gen)

  let mftaylor =
    foreign "mftaylor" (gen @-> long @-> long @-> long @-> returning gen)

  let mftobasis = foreign "mftobasis" (gen @-> gen @-> long @-> returning gen)
  let mftobasises = foreign "mftobasisES" (gen @-> gen @-> returning gen)
  let mftocol = foreign "mftocol" (gen @-> long @-> long @-> returning gen)

  let mftocoset =
    foreign "mftocoset" (pari_ulong @-> gen @-> gen @-> returning gen)

  let mftonew = foreign "mftonew" (gen @-> gen @-> returning gen)
  let mftraceform = foreign "mftraceform" (gen @-> long @-> returning gen)
  let mftwist = foreign "mftwist" (gen @-> gen @-> returning gen)
  let mfvecembed = foreign "mfvecembed" (gen @-> gen @-> returning gen)
  let mfvectomat = foreign "mfvectomat" (gen @-> long @-> long @-> returning gen)

  let fl_inv =
    foreign "Fl_inv" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_invsafe =
    foreign "Fl_invsafe" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fp_ratlift =
    foreign "Fp_ratlift"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning int)

  let zm2_mul = foreign "ZM2_mul" (gen @-> gen @-> returning gen)
  let abscmpii = foreign "abscmpii" (gen @-> gen @-> returning int)
  let abscmprr = foreign "abscmprr" (gen @-> gen @-> returning int)
end

module F43 (F : Ctypes.FOREIGN) = struct
  open F

  let absequalii = foreign "absequalii" (gen @-> gen @-> returning int)

  let addii_sign =
    foreign "addii_sign" (gen @-> long @-> gen @-> long @-> returning gen)

  let addir_sign =
    foreign "addir_sign" (gen @-> long @-> gen @-> long @-> returning gen)

  let addmulii = foreign "addmulii" (gen @-> gen @-> gen @-> returning gen)

  let addmulii_inplace =
    foreign "addmulii_inplace" (gen @-> gen @-> gen @-> returning gen)

  let addrr_sign =
    foreign "addrr_sign" (gen @-> long @-> gen @-> long @-> returning gen)

  let addsi_sign = foreign "addsi_sign" (long @-> gen @-> long @-> returning gen)
  let addsr = foreign "addsr" (long @-> gen @-> returning gen)

  let addui_sign =
    foreign "addui_sign" (pari_ulong @-> gen @-> long @-> returning gen)

  let addumului =
    foreign "addumului" (pari_ulong @-> pari_ulong @-> gen @-> returning gen)

  let affir = foreign "affir" (gen @-> gen @-> returning void)
  let affrr = foreign "affrr" (gen @-> gen @-> returning void)

  let bezout =
    foreign "bezout" (gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let cbezout =
    foreign "cbezout"
      (long @-> long @-> ptr long @-> ptr long @-> returning long)

  let cgcd = foreign "cgcd" (long @-> long @-> returning long)
  let clcm = foreign "clcm" (long @-> long @-> returning long)
  let cmpii = foreign "cmpii" (gen @-> gen @-> returning int)
  let cmprr = foreign "cmprr" (gen @-> gen @-> returning int)
  let dblexpo = foreign "dblexpo" (double @-> returning long)
  let dblmantissa = foreign "dblmantissa" (double @-> returning pari_ulong)
  let dbltor = foreign "dbltor" (double @-> returning gen)
  let diviiexact = foreign "diviiexact" (gen @-> gen @-> returning gen)
  let divir = foreign "divir" (gen @-> gen @-> returning gen)
  let divis = foreign "divis" (gen @-> long @-> returning gen)

  let divis_rem =
    foreign "divis_rem" (gen @-> long @-> ptr long @-> returning gen)

  let absdiviu_rem =
    foreign "absdiviu_rem"
      (gen @-> pari_ulong @-> ptr pari_ulong @-> returning gen)

  let diviuuexact =
    foreign "diviuuexact" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let diviuexact = foreign "diviuexact" (gen @-> pari_ulong @-> returning gen)
  let divri = foreign "divri" (gen @-> gen @-> returning gen)
  let divrr = foreign "divrr" (gen @-> gen @-> returning gen)
  let divrs = foreign "divrs" (gen @-> long @-> returning gen)
  let divru = foreign "divru" (gen @-> pari_ulong @-> returning gen)
  let divsi = foreign "divsi" (long @-> gen @-> returning gen)
  let divsr = foreign "divsr" (long @-> gen @-> returning gen)
  let divur = foreign "divur" (pari_ulong @-> gen @-> returning gen)
  let dvmdii = foreign "dvmdii" (gen @-> gen @-> ptr gen @-> returning gen)
  let equalii = foreign "equalii" (gen @-> gen @-> returning int)
  let equalrr = foreign "equalrr" (gen @-> gen @-> returning int)
  let floorr = foreign "floorr" (gen @-> returning gen)
  let gcdii = foreign "gcdii" (gen @-> gen @-> returning gen)
  let halfgcdii = foreign "halfgcdii" (gen @-> gen @-> returning gen)
  let int2n = foreign "int2n" (long @-> returning gen)
  let int2u = foreign "int2u" (pari_ulong @-> returning gen)
  let int2um1 = foreign "int2um1" (pari_ulong @-> returning gen)
  let int_normalize = foreign "int_normalize" (gen @-> long @-> returning gen)
  let invmod = foreign "invmod" (gen @-> gen @-> ptr gen @-> returning int)
  let invmod2bil = foreign "invmod2BIL" (pari_ulong @-> returning pari_ulong)
  let invr = foreign "invr" (gen @-> returning gen)

  let mantissa_real =
    foreign "mantissa_real" (gen @-> ptr long @-> returning gen)

  let modii = foreign "modii" (gen @-> gen @-> returning gen)
  let modiiz = foreign "modiiz" (gen @-> gen @-> gen @-> returning void)
  let mulii = foreign "mulii" (gen @-> gen @-> returning gen)
  let mulir = foreign "mulir" (gen @-> gen @-> returning gen)
  let mulrr = foreign "mulrr" (gen @-> gen @-> returning gen)
  let mulsi = foreign "mulsi" (long @-> gen @-> returning gen)
  let mulsr = foreign "mulsr" (long @-> gen @-> returning gen)
  let mulss = foreign "mulss" (long @-> long @-> returning gen)
  let mului = foreign "mului" (pari_ulong @-> gen @-> returning gen)
  let mulur = foreign "mulur" (pari_ulong @-> gen @-> returning gen)
  let muluu = foreign "muluu" (pari_ulong @-> pari_ulong @-> returning gen)

  let muluui =
    foreign "muluui" (pari_ulong @-> pari_ulong @-> gen @-> returning gen)

  let pari_kernel_close = foreign "pari_kernel_close" (void @-> returning void)
  let pari_kernel_init = foreign "pari_kernel_init" (void @-> returning void)

  let pari_kernel_version =
    foreign "pari_kernel_version" (void @-> returning string)

  let remi2n = foreign "remi2n" (gen @-> long @-> returning gen)
  let rtodbl = foreign "rtodbl" (gen @-> returning double)
  let shifti = foreign "shifti" (gen @-> long @-> returning gen)
  let sqri = foreign "sqri" (gen @-> returning gen)
  let sqrr = foreign "sqrr" (gen @-> returning gen)
  let sqrs = foreign "sqrs" (long @-> returning gen)
  let sqrtr_abs = foreign "sqrtr_abs" (gen @-> returning gen)
  let sqrtremi = foreign "sqrtremi" (gen @-> ptr gen @-> returning gen)
  let sqru = foreign "sqru" (pari_ulong @-> returning gen)
  let subsr = foreign "subsr" (long @-> gen @-> returning gen)

  let truedvmdii =
    foreign "truedvmdii" (gen @-> gen @-> ptr gen @-> returning gen)

  let truedvmdis =
    foreign "truedvmdis" (gen @-> long @-> ptr gen @-> returning gen)

  let truedvmdsi =
    foreign "truedvmdsi" (long @-> gen @-> ptr gen @-> returning gen)

  let trunc2nr = foreign "trunc2nr" (gen @-> long @-> returning gen)
  let mantissa2nr = foreign "mantissa2nr" (gen @-> long @-> returning gen)
  let truncr = foreign "truncr" (gen @-> returning gen)
  let ugcd = foreign "ugcd" (pari_ulong @-> pari_ulong @-> returning pari_ulong)
  let ulcm = foreign "ulcm" (pari_ulong @-> pari_ulong @-> returning pari_ulong)
  let umodiu = foreign "umodiu" (gen @-> pari_ulong @-> returning pari_ulong)
  let vals = foreign "vals" (pari_ulong @-> returning long)

  let fpc_ratlift =
    foreign "FpC_ratlift" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpm_ratlift =
    foreign "FpM_ratlift" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpx_ratlift =
    foreign "FpX_ratlift" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let qxqx_gcd = foreign "QXQX_gcd" (gen @-> gen @-> gen @-> returning gen)
  let zxqx_gcd = foreign "ZXQX_gcd" (gen @-> gen @-> gen @-> returning gen)
  let nffactor = foreign "nffactor" (gen @-> gen @-> returning gen)
  let nffactormod = foreign "nffactormod" (gen @-> gen @-> gen @-> returning gen)
  let nfgcd = foreign "nfgcd" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let nfgcd_all =
    foreign "nfgcd_all"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let nfissquarefree = foreign "nfissquarefree" (gen @-> gen @-> returning int)
  let nfroots = foreign "nfroots" (gen @-> gen @-> returning gen)

  let nfroots_if_split =
    foreign "nfroots_if_split" (ptr gen @-> gen @-> returning gen)

  let nfrootsof1 = foreign "nfrootsof1" (gen @-> returning gen)
  let polfnf = foreign "polfnf" (gen @-> gen @-> returning gen)

  let rnfabelianconjgen =
    foreign "rnfabelianconjgen" (gen @-> gen @-> returning gen)

  let rnfisabelian = foreign "rnfisabelian" (gen @-> gen @-> returning long)
end

module F44 (F : Ctypes.FOREIGN) = struct
  open F

  let forpart =
    foreign "forpart"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> long @-> gen @-> gen @-> returning void)

  let forpart_init =
    foreign "forpart_init"
      (ptr forpart_t @-> long @-> gen @-> gen @-> returning void)

  let forpart_next = foreign "forpart_next" (ptr forpart_t @-> returning gen)
  let forpart_prev = foreign "forpart_prev" (ptr forpart_t @-> returning gen)
  let numbpart = foreign "numbpart" (gen @-> returning gen)
  let partitions = foreign "partitions" (long @-> gen @-> gen @-> returning gen)

  let forperm =
    foreign "forperm"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> gen @-> returning void)

  let forperm_init =
    foreign "forperm_init" (ptr forperm_t @-> gen @-> returning void)

  let forperm_next = foreign "forperm_next" (ptr forperm_t @-> returning gen)

  let forallsubset_init =
    foreign "forallsubset_init" (ptr forsubset_t @-> long @-> returning void)

  let forksubset_init =
    foreign "forksubset_init"
      (ptr forsubset_t @-> long @-> long @-> returning void)

  let forsubset_next =
    foreign "forsubset_next" (ptr forsubset_t @-> returning gen)

  let forsubset_init =
    foreign "forsubset_init" (ptr forsubset_t @-> gen @-> returning void)

  let glambertw = foreign "glambertW" (gen @-> long @-> long @-> returning gen)
  let mplambertw = foreign "mplambertW" (gen @-> long @-> returning gen)
  let mplambertx = foreign "mplambertX" (gen @-> long @-> returning gen)

  let mplambertx_logx =
    foreign "mplambertx_logx" (gen @-> gen @-> long @-> returning gen)

  let mplambertxlogx_x =
    foreign "mplambertxlogx_x" (gen @-> gen @-> long @-> returning gen)

  let z_to_perm = foreign "Z_to_perm" (long @-> gen @-> returning gen)
  let abelian_group = foreign "abelian_group" (gen @-> returning gen)

  let conjclasses_repr =
    foreign "conjclasses_repr" (gen @-> long @-> returning gen)

  let cyc_pow = foreign "cyc_pow" (gen @-> long @-> returning gen)
  let cyc_pow_perm = foreign "cyc_pow_perm" (gen @-> long @-> returning gen)
  let cyclicgroup = foreign "cyclicgroup" (gen @-> long @-> returning gen)

  let dicyclicgroup =
    foreign "dicyclicgroup" (gen @-> gen @-> long @-> long @-> returning gen)

  let group_abelianhnf =
    foreign "group_abelianHNF" (gen @-> gen @-> returning gen)

  let group_abeliansnf =
    foreign "group_abelianSNF" (gen @-> gen @-> returning gen)

  let group_domain = foreign "group_domain" (gen @-> returning long)
  let group_elts = foreign "group_elts" (gen @-> long @-> returning gen)
  let group_export = foreign "group_export" (gen @-> long @-> returning gen)
  let group_export_gap = foreign "group_export_GAP" (gen @-> returning gen)
  let group_export_magma = foreign "group_export_MAGMA" (gen @-> returning gen)
  let group_isa4s4 = foreign "group_isA4S4" (gen @-> returning long)
  let group_isabelian = foreign "group_isabelian" (gen @-> returning long)
  let group_leftcoset = foreign "group_leftcoset" (gen @-> gen @-> returning gen)
  let group_order = foreign "group_order" (gen @-> returning long)

  let group_perm_normalize =
    foreign "group_perm_normalize" (gen @-> gen @-> returning long)

  let group_quotient = foreign "group_quotient" (gen @-> gen @-> returning gen)

  let group_rightcoset =
    foreign "group_rightcoset" (gen @-> gen @-> returning gen)

  let group_set = foreign "group_set" (gen @-> long @-> returning gen)

  let group_subgroup_is_faithful =
    foreign "group_subgroup_is_faithful" (gen @-> gen @-> returning int)

  let group_subgroup_isnormal =
    foreign "group_subgroup_isnormal" (gen @-> gen @-> returning long)

  let group_subgroups = foreign "group_subgroups" (gen @-> returning gen)

  let groupelts_solvablesubgroups =
    foreign "groupelts_solvablesubgroups" (gen @-> returning gen)

  let group_to_cc = foreign "group_to_cc" (gen @-> returning gen)

  let groupelts_abelian_group =
    foreign "groupelts_abelian_group" (gen @-> returning gen)

  let groupelts_center = foreign "groupelts_center" (gen @-> returning gen)

  let groupelts_conj_set =
    foreign "groupelts_conj_set" (gen @-> gen @-> returning gen)

  let groupelts_conjclasses =
    foreign "groupelts_conjclasses" (gen @-> ptr long @-> returning gen)

  let groupelts_exponent = foreign "groupelts_exponent" (gen @-> returning long)

  let groupelts_quotient =
    foreign "groupelts_quotient" (gen @-> gen @-> returning gen)

  let groupelts_set = foreign "groupelts_set" (gen @-> long @-> returning gen)
  let groupelts_to_group = foreign "groupelts_to_group" (gen @-> returning gen)
  let numtoperm = foreign "numtoperm" (long @-> gen @-> returning gen)
  let perm_commute = foreign "perm_commute" (gen @-> gen @-> returning int)
  let perm_cycles = foreign "perm_cycles" (gen @-> returning gen)
  let perm_order = foreign "perm_order" (gen @-> returning gen)
  let perm_orderu = foreign "perm_orderu" (gen @-> returning pari_ulong)
  let perm_pow = foreign "perm_pow" (gen @-> gen @-> returning gen)
  let perm_powu = foreign "perm_powu" (gen @-> pari_ulong @-> returning gen)
  let perm_sign = foreign "perm_sign" (gen @-> returning long)
  let perm_to_gap = foreign "perm_to_GAP" (gen @-> returning gen)
  let perm_to_z = foreign "perm_to_Z" (gen @-> returning gen)
  let permcycles = foreign "permcycles" (gen @-> returning gen)
  let permorder = foreign "permorder" (gen @-> returning gen)
  let permsign = foreign "permsign" (gen @-> returning long)
  let permtonum = foreign "permtonum" (gen @-> returning gen)
  let quotient_group = foreign "quotient_group" (gen @-> gen @-> returning gen)
  let quotient_groupelts = foreign "quotient_groupelts" (gen @-> returning gen)
  let quotient_perm = foreign "quotient_perm" (gen @-> gen @-> returning gen)

  let quotient_subgroup_lift =
    foreign "quotient_subgroup_lift" (gen @-> gen @-> gen @-> returning gen)

  let subgroups_tableset =
    foreign "subgroups_tableset" (gen @-> long @-> returning gen)

  let tableset_find_index =
    foreign "tableset_find_index" (gen @-> gen @-> returning long)

  let trivialgroup = foreign "trivialgroup" (void @-> returning gen)
  let vec_insert = foreign "vec_insert" (gen @-> long @-> gen @-> returning gen)
  let vec_is1to1 = foreign "vec_is1to1" (gen @-> returning int)
  let vec_isconst = foreign "vec_isconst" (gen @-> returning int)
  let vecperm_orbits = foreign "vecperm_orbits" (gen @-> long @-> returning gen)
  let vecsmall_duplicate = foreign "vecsmall_duplicate" (gen @-> returning long)

  let vecsmall_duplicate_sorted =
    foreign "vecsmall_duplicate_sorted" (gen @-> returning long)

  let vecsmall_indexsort = foreign "vecsmall_indexsort" (gen @-> returning gen)
  let vecsmall_is1to1 = foreign "vecsmall_is1to1" (gen @-> returning int)
  let vecsmall_isconst = foreign "vecsmall_isconst" (gen @-> returning int)
  let vecsmall_sort = foreign "vecsmall_sort" (gen @-> returning void)
  let vecsmall_uniq = foreign "vecsmall_uniq" (gen @-> returning gen)

  let vecsmall_uniq_sorted =
    foreign "vecsmall_uniq_sorted" (gen @-> returning gen)

  let vecsmall_counting_indexsort =
    foreign "vecsmall_counting_indexsort" (gen @-> long @-> returning gen)

  let vecsmall_counting_sort =
    foreign "vecsmall_counting_sort" (gen @-> long @-> returning void)

  let vecsmall_counting_uniq =
    foreign "vecsmall_counting_uniq" (gen @-> long @-> returning gen)

  let vecvecsmall_indexsort =
    foreign "vecvecsmall_indexsort" (gen @-> returning gen)

  let vecvecsmall_max = foreign "vecvecsmall_max" (gen @-> returning long)

  let vecvecsmall_search =
    foreign "vecvecsmall_search" (gen @-> gen @-> returning long)

  let vecvecsmall_sort = foreign "vecvecsmall_sort" (gen @-> returning gen)

  let vecvecsmall_sort_inplace =
    foreign "vecvecsmall_sort_inplace" (gen @-> ptr gen @-> returning void)

  let vecvecsmall_sort_shallow =
    foreign "vecvecsmall_sort_shallow" (gen @-> returning gen)

  let vecvecsmall_sort_uniq =
    foreign "vecvecsmall_sort_uniq" (gen @-> returning gen)

  let mt_broadcast = foreign "mt_broadcast" (gen @-> returning void)
  let mt_nbthreads = foreign "mt_nbthreads" (void @-> returning long)
  let mt_queue_end = foreign "mt_queue_end" (ptr pari_mt @-> returning void)

  let mt_queue_get =
    foreign "mt_queue_get"
      (ptr pari_mt @-> ptr long @-> ptr long @-> returning gen)
end

module F45 (F : Ctypes.FOREIGN) = struct
  open F

  let mt_queue_start =
    foreign "mt_queue_start" (ptr pari_mt @-> gen @-> returning void)

  let mt_queue_start_lim =
    foreign "mt_queue_start_lim"
      (ptr pari_mt @-> gen @-> long @-> returning void)

  let mt_queue_submit =
    foreign "mt_queue_submit" (ptr pari_mt @-> long @-> gen @-> returning void)

  let mt_sigint_block = foreign "mt_sigint_block" (void @-> returning void)
  let mt_sigint_unblock = foreign "mt_sigint_unblock" (void @-> returning void)
  let pari_mt_init = foreign "pari_mt_init" (void @-> returning void)
  let pari_mt_close = foreign "pari_mt_close" (void @-> returning void)

  let subcyclopclgp =
    foreign "subcyclopclgp" (gen @-> gen @-> long @-> returning gen)

  let subcycloiwasawa =
    foreign "subcycloiwasawa" (gen @-> gen @-> long @-> returning gen)

  let subcyclohminus = foreign "subcyclohminus" (gen @-> gen @-> returning gen)

  let color_to_rgb =
    foreign "color_to_rgb"
      (gen @-> ptr int @-> ptr int @-> ptr int @-> returning void)

  let colorname_to_rgb =
    foreign "colorname_to_rgb"
      (string @-> ptr int @-> ptr int @-> ptr int @-> returning void)

  let long_to_rgb =
    foreign "long_to_rgb"
      (long @-> ptr int @-> ptr int @-> ptr int @-> returning void)

  let pari_plot_by_file =
    foreign "pari_plot_by_file" (string @-> string @-> string @-> returning void)

  let pari_set_plot_engine =
    foreign "pari_set_plot_engine"
      (static_funptr Ctypes.(ptr pari_plot @-> returning void)
      @-> returning void)

  let pari_kill_plot_engine =
    foreign "pari_kill_plot_engine" (void @-> returning void)

  let parploth =
    foreign "parploth"
      (gen @-> gen @-> gen @-> long @-> long @-> long @-> returning gen)

  let parplothexport =
    foreign "parplothexport"
      (gen @-> gen @-> gen @-> gen @-> long @-> long @-> long @-> returning gen)

  let plotbox =
    foreign "plotbox" (long @-> gen @-> gen @-> long @-> returning void)

  let plotclip = foreign "plotclip" (long @-> returning void)
  let plotcolor = foreign "plotcolor" (long @-> gen @-> returning gen)

  let plotcopy =
    foreign "plotcopy"
      (long @-> long @-> gen @-> gen @-> long @-> returning void)

  let plotcursor = foreign "plotcursor" (long @-> returning gen)
  let plotdraw = foreign "plotdraw" (gen @-> long @-> returning void)
  let plotexport = foreign "plotexport" (gen @-> gen @-> long @-> returning gen)

  let ploth =
    foreign "ploth"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> long @-> long @-> returning gen)

  let plothexport =
    foreign "plothexport"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> long @-> long @-> returning gen)

  let plothraw = foreign "plothraw" (gen @-> gen @-> long @-> returning gen)

  let plothrawexport =
    foreign "plothrawexport" (gen @-> gen @-> gen @-> long @-> returning gen)

  let plothsizes = foreign "plothsizes" (long @-> returning gen)

  let plotinit =
    foreign "plotinit" (long @-> gen @-> gen @-> long @-> returning void)

  let plotkill = foreign "plotkill" (long @-> returning void)
  let plotline = foreign "plotline" (long @-> gen @-> gen @-> returning void)

  let plotlines =
    foreign "plotlines" (long @-> gen @-> gen @-> long @-> returning void)

  let plotlinetype = foreign "plotlinetype" (long @-> long @-> returning void)
  let plotmove = foreign "plotmove" (long @-> gen @-> gen @-> returning void)
  let plotpoints = foreign "plotpoints" (long @-> gen @-> gen @-> returning void)
  let plotpointsize = foreign "plotpointsize" (long @-> gen @-> returning void)
  let plotpointtype = foreign "plotpointtype" (long @-> long @-> returning void)

  let plotrbox =
    foreign "plotrbox" (long @-> gen @-> gen @-> long @-> returning void)

  let plotrecth =
    foreign "plotrecth"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> long @-> gen @-> gen @-> pari_ulong @-> long @-> long
      @-> returning gen)

  let plotrecthraw =
    foreign "plotrecthraw" (long @-> gen @-> long @-> returning gen)

  let plotrline = foreign "plotrline" (long @-> gen @-> gen @-> returning void)
  let plotrmove = foreign "plotrmove" (long @-> gen @-> gen @-> returning void)
  let plotrpoint = foreign "plotrpoint" (long @-> gen @-> gen @-> returning void)

  let plotscale =
    foreign "plotscale" (long @-> gen @-> gen @-> gen @-> gen @-> returning void)

  let plotstring =
    foreign "plotstring" (long @-> string @-> long @-> returning void)

  let psdraw = foreign "psdraw" (gen @-> long @-> returning void)

  let psploth =
    foreign "psploth"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> long @-> long @-> returning gen)

  let psplothraw = foreign "psplothraw" (gen @-> gen @-> long @-> returning gen)

  let rect2ps =
    foreign "rect2ps"
      (gen @-> gen @-> gen @-> ptr pari_plot @-> returning string)

  let rect2ps_i =
    foreign "rect2ps_i"
      (gen @-> gen @-> gen @-> ptr pari_plot @-> int @-> returning string)

  let rect2svg =
    foreign "rect2svg"
      (gen @-> gen @-> gen @-> ptr pari_plot @-> returning string)

  let pariplot =
    foreign "pariplot"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> gen @-> gen @-> long @-> returning void)

  let zx_zp_root =
    foreign "ZX_Zp_root" (gen @-> gen @-> gen @-> long @-> returning gen)

  let zp_appr = foreign "Zp_appr" (gen @-> gen @-> returning gen)
  let cmp_padic = foreign "cmp_padic" (gen @-> gen @-> returning int)

  let factorpadic =
    foreign "factorpadic" (gen @-> gen @-> long @-> returning gen)

  let gdeuc = foreign "gdeuc" (gen @-> gen @-> returning gen)
  let grem = foreign "grem" (gen @-> gen @-> returning gen)
  let padicappr = foreign "padicappr" (gen @-> gen @-> returning gen)
  let poldivrem = foreign "poldivrem" (gen @-> gen @-> ptr gen @-> returning gen)

  let polrootspadic =
    foreign "polrootspadic" (gen @-> gen @-> long @-> returning gen)

  let flv_factorback =
    foreign "Flv_factorback"
      (gen @-> gen @-> pari_ulong @-> returning pari_ulong)

  let flxqv_factorback =
    foreign "FlxqV_factorback"
      (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let fpv_factorback =
    foreign "FpV_factorback" (gen @-> gen @-> gen @-> returning gen)

  let fqv_factorback =
    foreign "FqV_factorback" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let q_content = foreign "Q_content" (gen @-> returning gen)
  let q_content_safe = foreign "Q_content_safe" (gen @-> returning gen)
  let q_denom = foreign "Q_denom" (gen @-> returning gen)
  let q_denom_safe = foreign "Q_denom_safe" (gen @-> returning gen)
  let q_div_to_int = foreign "Q_div_to_int" (gen @-> gen @-> returning gen)
  let q_gcd = foreign "Q_gcd" (gen @-> gen @-> returning gen)
  let q_mul_to_int = foreign "Q_mul_to_int" (gen @-> gen @-> returning gen)
  let q_muli_to_int = foreign "Q_muli_to_int" (gen @-> gen @-> returning gen)

  let q_primitive_part =
    foreign "Q_primitive_part" (gen @-> ptr gen @-> returning gen)

  let q_primpart = foreign "Q_primpart" (gen @-> returning gen)

  let q_remove_denom =
    foreign "Q_remove_denom" (gen @-> ptr gen @-> returning gen)

  let q_factor = foreign "Q_factor" (gen @-> returning gen)

  let q_factor_limit =
    foreign "Q_factor_limit" (gen @-> pari_ulong @-> returning gen)

  let rg_type =
    foreign "Rg_type"
      (gen @-> ptr gen @-> ptr gen @-> ptr long @-> returning long)

  let rgm_rgc_type =
    foreign "RgM_RgC_type"
      (gen @-> gen @-> ptr gen @-> ptr gen @-> ptr long @-> returning long)

  let rgm_rescale_to_int = foreign "RgM_rescale_to_int" (gen @-> returning gen)

  let rgm_type =
    foreign "RgM_type"
      (gen @-> ptr gen @-> ptr gen @-> ptr long @-> returning long)

  let rgm_type2 =
    foreign "RgM_type2"
      (gen @-> gen @-> ptr gen @-> ptr gen @-> ptr long @-> returning long)

  let rgv_type =
    foreign "RgV_type"
      (gen @-> ptr gen @-> ptr gen @-> ptr long @-> returning long)

  let rgv_type2 =
    foreign "RgV_type2"
      (gen @-> gen @-> ptr gen @-> ptr gen @-> ptr long @-> returning long)

  let rgx_rg_type =
    foreign "RgX_Rg_type"
      (gen @-> gen @-> ptr gen @-> ptr gen @-> ptr long @-> returning long)

  let rgx_chinese_coprime =
    foreign "RgX_chinese_coprime"
      (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let rgx_disc = foreign "RgX_disc" (gen @-> returning gen)

  let rgx_extgcd =
    foreign "RgX_extgcd" (gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let rgx_extgcd_simple =
    foreign "RgX_extgcd_simple"
      (gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let rgx_gcd = foreign "RgX_gcd" (gen @-> gen @-> returning gen)
  let rgx_gcd_simple = foreign "RgX_gcd_simple" (gen @-> gen @-> returning gen)
  let rgx_halfgcd = foreign "RgX_halfgcd" (gen @-> gen @-> returning gen)
  let rgx_rescale_to_int = foreign "RgX_rescale_to_int" (gen @-> returning gen)

  let rgx_resultant_all =
    foreign "RgX_resultant_all" (gen @-> gen @-> ptr gen @-> returning gen)

  let rgx_sturmpart = foreign "RgX_sturmpart" (gen @-> gen @-> returning long)

  let rgx_sylvestermatrix =
    foreign "RgX_sylvestermatrix" (gen @-> gen @-> returning gen)

  let rgx_type =
    foreign "RgX_type"
      (gen @-> ptr gen @-> ptr gen @-> ptr long @-> returning long)
end

module F46 (F : Ctypes.FOREIGN) = struct
  open F

  let rgx_type2 =
    foreign "RgX_type2"
      (gen @-> gen @-> ptr gen @-> ptr gen @-> ptr long @-> returning long)

  let rgx_type3 =
    foreign "RgX_type3"
      (gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> ptr long
     @-> returning long)

  let rgx_type_decode =
    foreign "RgX_type_decode" (long @-> ptr long @-> ptr long @-> returning void)

  let rgx_type_is_composite =
    foreign "RgX_type_is_composite" (long @-> returning int)

  let rgxq_charpoly =
    foreign "RgXQ_charpoly" (gen @-> gen @-> long @-> returning gen)

  let rgxq_inv = foreign "RgXQ_inv" (gen @-> gen @-> returning gen)

  let rgxq_minpoly =
    foreign "RgXQ_minpoly" (gen @-> gen @-> long @-> returning gen)

  let rgxq_mul = foreign "RgXQ_mul" (gen @-> gen @-> gen @-> returning gen)

  let rgxq_ratlift =
    foreign "RgXQ_ratlift"
      (gen @-> gen @-> long @-> long @-> ptr gen @-> ptr gen @-> returning int)

  let rgxq_sqr = foreign "RgXQ_sqr" (gen @-> gen @-> returning gen)
  let z_content = foreign "Z_content" (gen @-> returning gen)
  let zx_content = foreign "ZX_content" (gen @-> returning gen)
  let centermod = foreign "centermod" (gen @-> gen @-> returning gen)
  let centermod_i = foreign "centermod_i" (gen @-> gen @-> gen @-> returning gen)
  let centermodii = foreign "centermodii" (gen @-> gen @-> gen @-> returning gen)
  let content = foreign "content" (gen @-> returning gen)
  let content0 = foreign "content0" (gen @-> gen @-> returning gen)

  let deg1_from_roots =
    foreign "deg1_from_roots" (gen @-> long @-> returning gen)

  let factor = foreign "factor" (gen @-> returning gen)
  let factor0 = foreign "factor0" (gen @-> gen @-> returning gen)
  let factorback = foreign "factorback" (gen @-> returning gen)
  let factorback2 = foreign "factorback2" (gen @-> gen @-> returning gen)

  let gbezout =
    foreign "gbezout" (gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let gdivexact = foreign "gdivexact" (gen @-> gen @-> returning gen)

  let gen_factorback =
    foreign "gen_factorback"
      (gen @-> gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> gen @-> gen @-> returning gen)
      @-> static_funptr Ctypes.(ptr void @-> returning gen)
      @-> returning gen)

  let ggcd = foreign "ggcd" (gen @-> gen @-> returning gen)
  let ggcd0 = foreign "ggcd0" (gen @-> gen @-> returning gen)
  let ghalfgcd = foreign "ghalfgcd" (gen @-> gen @-> returning gen)
  let ginvmod = foreign "ginvmod" (gen @-> gen @-> returning gen)
  let glcm = foreign "glcm" (gen @-> gen @-> returning gen)
  let glcm0 = foreign "glcm0" (gen @-> gen @-> returning gen)
  let newtonpoly = foreign "newtonpoly" (gen @-> gen @-> returning gen)
  let nfrootsq = foreign "nfrootsQ" (gen @-> returning gen)
  let poldisc0 = foreign "poldisc0" (gen @-> long @-> returning gen)
  let polisirreducible = foreign "polisirreducible" (gen @-> returning long)

  let polresultant0 =
    foreign "polresultant0" (gen @-> gen @-> long @-> long @-> returning gen)

  let polsym = foreign "polsym" (gen @-> long @-> returning gen)

  let primitive_part =
    foreign "primitive_part" (gen @-> ptr gen @-> returning gen)

  let primpart = foreign "primpart" (gen @-> returning gen)
  let reduceddiscsmith = foreign "reduceddiscsmith" (gen @-> returning gen)
  let resultant2 = foreign "resultant2" (gen @-> gen @-> returning gen)
  let resultant = foreign "resultant" (gen @-> gen @-> returning gen)

  let rnfcharpoly =
    foreign "rnfcharpoly" (gen @-> gen @-> gen @-> long @-> returning gen)

  let roots_from_deg1 = foreign "roots_from_deg1" (gen @-> returning gen)
  let roots_to_pol = foreign "roots_to_pol" (gen @-> long @-> returning gen)

  let roots_to_pol_r1 =
    foreign "roots_to_pol_r1" (gen @-> long @-> long @-> returning gen)

  let sturmpart = foreign "sturmpart" (gen @-> gen @-> gen @-> returning long)

  let subresext =
    foreign "subresext" (gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let sylvestermatrix = foreign "sylvestermatrix" (gen @-> gen @-> returning gen)
  let trivial_fact = foreign "trivial_fact" (void @-> returning gen)
  let gcdext0 = foreign "gcdext0" (gen @-> gen @-> returning gen)

  let polresultantext0 =
    foreign "polresultantext0" (gen @-> gen @-> long @-> returning gen)

  let polresultantext = foreign "polresultantext" (gen @-> gen @-> returning gen)
  let prime_fact = foreign "prime_fact" (gen @-> returning gen)
  let row_q_primpart = foreign "row_Q_primpart" (gen @-> returning gen)
  let vec_q_primpart = foreign "vec_Q_primpart" (gen @-> returning gen)
  let vecprod = foreign "vecprod" (gen @-> returning gen)
  let zv_lcm = foreign "ZV_lcm" (gen @-> returning gen)

  let flx_flxy_resultant =
    foreign "Flx_FlxY_resultant" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flxx_resultant =
    foreign "FlxX_resultant"
      (gen @-> gen @-> pari_ulong @-> long @-> returning gen)

  let fpx_fpxy_resultant =
    foreign "FpX_FpXY_resultant" (gen @-> gen @-> gen @-> returning gen)

  let fpx_translate =
    foreign "FpX_translate" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_normalize =
    foreign "FpXQX_normalize" (gen @-> gen @-> gen @-> returning gen)

  let fpxv_fpc_mul =
    foreign "FpXV_FpC_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpxy_fpxq_evaly =
    foreign "FpXY_FpXQ_evaly"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let fpxc_center = foreign "FpXC_center" (gen @-> gen @-> gen @-> returning gen)
  let fpxm_center = foreign "FpXM_center" (gen @-> gen @-> gen @-> returning gen)

  let fq_fp_mul =
    foreign "Fq_Fp_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fq_add = foreign "Fq_add" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fq_div = foreign "Fq_div" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fq_halve = foreign "Fq_halve" (gen @-> gen @-> gen @-> returning gen)
  let fq_inv = foreign "Fq_inv" (gen @-> gen @-> gen @-> returning gen)
  let fq_invsafe = foreign "Fq_invsafe" (gen @-> gen @-> gen @-> returning gen)
  let fq_mul = foreign "Fq_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fq_mulu =
    foreign "Fq_mulu" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fq_neg = foreign "Fq_neg" (gen @-> gen @-> gen @-> returning gen)
  let fq_neg_inv = foreign "Fq_neg_inv" (gen @-> gen @-> gen @-> returning gen)
  let fq_pow = foreign "Fq_pow" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fq_powu =
    foreign "Fq_powu" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fq_sqr = foreign "Fq_sqr" (gen @-> gen @-> gen @-> returning gen)
  let fq_sqrt = foreign "Fq_sqrt" (gen @-> gen @-> gen @-> returning gen)

  let fq_sqrtn =
    foreign "Fq_sqrtn"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fq_sub = foreign "Fq_sub" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqc_fq_mul =
    foreign "FqC_Fq_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqc_fqv_mul =
    foreign "FqC_FqV_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqc_add = foreign "FqC_add" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fqc_sub = foreign "FqC_sub" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fqv_red = foreign "FqV_red" (gen @-> gen @-> gen @-> returning gen)

  let fqv_roots_to_pol =
    foreign "FqV_roots_to_pol" (gen @-> gen @-> gen @-> long @-> returning gen)

  let fqx_fq_add =
    foreign "FqX_Fq_add" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_fq_mul_to_monic =
    foreign "FqX_Fq_mul_to_monic" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_fq_sub =
    foreign "FqX_Fq_sub" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_eval =
    foreign "FqX_eval" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_translate =
    foreign "FqX_translate" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_matrix_pow =
    foreign "FqXQ_matrix_pow"
      (gen @-> long @-> long @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_powers =
    foreign "FqXQ_powers"
      (gen @-> long @-> gen @-> gen @-> gen @-> returning gen)

  let fqxy_eval =
    foreign "FqXY_eval" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxy_evalx =
    foreign "FqXY_evalx" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let qx_disc = foreign "QX_disc" (gen @-> returning gen)
  let qx_gcd = foreign "QX_gcd" (gen @-> gen @-> returning gen)
end

module F47 (F : Ctypes.FOREIGN) = struct
  open F

  let qx_resultant = foreign "QX_resultant" (gen @-> gen @-> returning gen)
  let qxq_div = foreign "QXQ_div" (gen @-> gen @-> gen @-> returning gen)
  let qxq_intnorm = foreign "QXQ_intnorm" (gen @-> gen @-> returning gen)
  let qxq_inv = foreign "QXQ_inv" (gen @-> gen @-> returning gen)
  let qxq_mul = foreign "QXQ_mul" (gen @-> gen @-> gen @-> returning gen)
  let qxq_norm = foreign "QXQ_norm" (gen @-> gen @-> returning gen)
  let qxq_sqr = foreign "QXQ_sqr" (gen @-> gen @-> returning gen)
  let rg_is_fp = foreign "Rg_is_Fp" (gen @-> ptr gen @-> returning int)

  let rg_is_fpxq =
    foreign "Rg_is_FpXQ" (gen @-> ptr gen @-> ptr gen @-> returning int)

  let rg_to_fp = foreign "Rg_to_Fp" (gen @-> gen @-> returning gen)
  let rg_to_fpxq = foreign "Rg_to_FpXQ" (gen @-> gen @-> gen @-> returning gen)
  let rgc_to_fpc = foreign "RgC_to_FpC" (gen @-> gen @-> returning gen)
  let rgc_to_fqc = foreign "RgC_to_FqC" (gen @-> gen @-> gen @-> returning gen)
  let rgm_is_fpm = foreign "RgM_is_FpM" (gen @-> ptr gen @-> returning int)
  let rgm_to_flm = foreign "RgM_to_Flm" (gen @-> pari_ulong @-> returning gen)
  let rgm_to_fpm = foreign "RgM_to_FpM" (gen @-> gen @-> returning gen)
  let rgm_to_fqm = foreign "RgM_to_FqM" (gen @-> gen @-> gen @-> returning gen)
  let rgv_is_fpv = foreign "RgV_is_FpV" (gen @-> ptr gen @-> returning int)
  let rgv_to_flv = foreign "RgV_to_Flv" (gen @-> pari_ulong @-> returning gen)
  let rgv_to_fpv = foreign "RgV_to_FpV" (gen @-> gen @-> returning gen)
  let rgx_is_fpx = foreign "RgX_is_FpX" (gen @-> ptr gen @-> returning int)
  let rgx_to_fpx = foreign "RgX_to_FpX" (gen @-> gen @-> returning gen)

  let rgx_is_fpxqx =
    foreign "RgX_is_FpXQX" (gen @-> ptr gen @-> ptr gen @-> returning int)

  let rgx_to_fpxqx =
    foreign "RgX_to_FpXQX" (gen @-> gen @-> gen @-> returning gen)

  let rgx_to_fqx = foreign "RgX_to_FqX" (gen @-> gen @-> gen @-> returning gen)

  let z_incremental_crt =
    foreign "Z_incremental_CRT"
      (ptr gen @-> pari_ulong @-> ptr gen @-> pari_ulong @-> returning int)

  let z_init_crt =
    foreign "Z_init_CRT" (pari_ulong @-> pari_ulong @-> returning gen)

  let zm_incremental_crt =
    foreign "ZM_incremental_CRT"
      (ptr gen @-> gen @-> ptr gen @-> pari_ulong @-> returning int)

  let zm_init_crt = foreign "ZM_init_CRT" (gen @-> pari_ulong @-> returning gen)

  let zx_zxy_resultant =
    foreign "ZX_ZXY_resultant" (gen @-> gen @-> returning gen)

  let zx_zxy_rnfequation =
    foreign "ZX_ZXY_rnfequation" (gen @-> gen @-> ptr long @-> returning gen)

  let zx_disc = foreign "ZX_disc" (gen @-> returning gen)
  let zx_gcd = foreign "ZX_gcd" (gen @-> gen @-> returning gen)

  let zx_gcd_all =
    foreign "ZX_gcd_all" (gen @-> gen @-> ptr gen @-> returning gen)

  let zx_incremental_crt =
    foreign "ZX_incremental_CRT"
      (ptr gen @-> gen @-> ptr gen @-> pari_ulong @-> returning int)

  let zx_init_crt =
    foreign "ZX_init_CRT" (gen @-> pari_ulong @-> long @-> returning gen)

  let zx_is_squarefree = foreign "ZX_is_squarefree" (gen @-> returning int)
  let zx_radical = foreign "ZX_radical" (gen @-> returning gen)
  let zx_resultant = foreign "ZX_resultant" (gen @-> gen @-> returning gen)

  let zxm_incremental_crt =
    foreign "ZXM_incremental_CRT"
      (ptr gen @-> gen @-> ptr gen @-> pari_ulong @-> returning int)

  let zxm_init_crt =
    foreign "ZXM_init_CRT" (gen @-> long @-> pari_ulong @-> returning gen)

  let zxq_minpoly =
    foreign "ZXQ_minpoly" (gen @-> gen @-> long @-> pari_ulong @-> returning gen)

  let zxq_charpoly =
    foreign "ZXQ_charpoly" (gen @-> gen @-> long @-> returning gen)

  let characteristic = foreign "characteristic" (gen @-> returning gen)
  let ffinit = foreign "ffinit" (gen @-> long @-> long @-> returning gen)
  let ffnbirred = foreign "ffnbirred" (gen @-> long @-> returning gen)
  let ffnbirred0 = foreign "ffnbirred0" (gen @-> long @-> long @-> returning gen)
  let ffsumnbirred = foreign "ffsumnbirred" (gen @-> long @-> returning gen)

  let get_fq_field =
    foreign "get_Fq_field"
      (ptr (ptr void) @-> gen @-> gen @-> returning (ptr bb_field))

  let init_flxq =
    foreign "init_Flxq" (pari_ulong @-> long @-> long @-> returning gen)

  let init_fq = foreign "init_Fq" (gen @-> long @-> long @-> returning gen)
  let nfx_disc = foreign "nfX_disc" (gen @-> gen @-> returning gen)

  let nfx_resultant =
    foreign "nfX_resultant" (gen @-> gen @-> gen @-> returning gen)

  let pol_x_powers = foreign "pol_x_powers" (long @-> long @-> returning gen)

  let residual_characteristic =
    foreign "residual_characteristic" (gen @-> returning gen)

  let polclass = foreign "polclass" (gen @-> long @-> long @-> returning gen)

  let fp_modinv_to_j =
    foreign "Fp_modinv_to_j" (gen @-> long @-> gen @-> returning gen)

  let fp_polmodular_evalx =
    foreign "Fp_polmodular_evalx"
      (long @-> long @-> gen @-> gen @-> long @-> int @-> returning gen)

  let check_modinv = foreign "check_modinv" (long @-> returning void)
  let disc_best_modinv = foreign "disc_best_modinv" (long @-> returning long)

  let modinv_height_factor =
    foreign "modinv_height_factor" (long @-> returning long)

  let modinv_good_disc =
    foreign "modinv_good_disc" (long @-> long @-> returning int)

  let modinv_good_prime =
    foreign "modinv_good_prime" (long @-> long @-> returning int)

  let modinv_is_weber = foreign "modinv_is_Weber" (long @-> returning int)

  let modinv_is_double_eta =
    foreign "modinv_is_double_eta" (long @-> returning int)

  let polmodular =
    foreign "polmodular"
      (long @-> long @-> gen @-> long @-> long @-> returning gen)

  let polmodular_zm = foreign "polmodular_ZM" (long @-> long @-> returning gen)

  let polmodular_zxx =
    foreign "polmodular_ZXX" (long @-> long @-> long @-> long @-> returning gen)

  let bpsw_isprime = foreign "BPSW_isprime" (gen @-> returning long)
  let bpsw_psp = foreign "BPSW_psp" (gen @-> returning long)
  let addprimes = foreign "addprimes" (gen @-> returning gen)
  let check_ecppcert = foreign "check_ecppcert" (gen @-> returning long)
  let gisprime = foreign "gisprime" (gen @-> long @-> returning gen)
  let gispseudoprime = foreign "gispseudoprime" (gen @-> long @-> returning gen)

  let gprimepi_upper_bound =
    foreign "gprimepi_upper_bound" (gen @-> returning gen)

  let gprimepi_lower_bound =
    foreign "gprimepi_lower_bound" (gen @-> returning gen)

  let isprime = foreign "isprime" (gen @-> returning long)
  let ispseudoprime = foreign "ispseudoprime" (gen @-> long @-> returning long)
  let millerrabin = foreign "millerrabin" (gen @-> long @-> returning long)
  let prime = foreign "prime" (long @-> returning gen)
  let primecert = foreign "primecert" (gen @-> long @-> returning gen)
  let primecert0 = foreign "primecert0" (gen @-> long @-> long @-> returning gen)

  let primecertexport =
    foreign "primecertexport" (gen @-> long @-> returning gen)

  let primecertisvalid = foreign "primecertisvalid" (gen @-> returning long)
  let primepi = foreign "primepi" (gen @-> returning gen)

  let primepi_upper_bound =
    foreign "primepi_upper_bound" (double @-> returning double)

  let primepi_lower_bound =
    foreign "primepi_lower_bound" (double @-> returning double)

  let primes = foreign "primes" (long @-> returning gen)
  let primes_interval = foreign "primes_interval" (gen @-> gen @-> returning gen)

  let primes_interval_zv =
    foreign "primes_interval_zv" (pari_ulong @-> pari_ulong @-> returning gen)

  let primes_upto_zv = foreign "primes_upto_zv" (pari_ulong @-> returning gen)
  let primes0 = foreign "primes0" (gen @-> returning gen)
  let primes_zv = foreign "primes_zv" (long @-> returning gen)
  let randomprime = foreign "randomprime" (gen @-> returning gen)
  let randomprime0 = foreign "randomprime0" (gen @-> gen @-> returning gen)
  let removeprimes = foreign "removeprimes" (gen @-> returning gen)
  let uis2psp = foreign "uis2psp" (pari_ulong @-> returning int)
  let uispsp = foreign "uispsp" (pari_ulong @-> pari_ulong @-> returning int)
  let uislucaspsp = foreign "uislucaspsp" (pari_ulong @-> returning int)
  let uisprime = foreign "uisprime" (pari_ulong @-> returning int)
end

module F48 (F : Ctypes.FOREIGN) = struct
  open F

  let uisprime_101 = foreign "uisprime_101" (pari_ulong @-> returning int)
  let uisprime_661 = foreign "uisprime_661" (pari_ulong @-> returning int)
  let uprime = foreign "uprime" (long @-> returning pari_ulong)
  let uprimepi = foreign "uprimepi" (pari_ulong @-> returning pari_ulong)
  let qfauto = foreign "qfauto" (gen @-> gen @-> returning gen)
  let qfauto0 = foreign "qfauto0" (gen @-> gen @-> returning gen)
  let qfautoexport = foreign "qfautoexport" (gen @-> long @-> returning gen)
  let qfisom = foreign "qfisom" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let qfisom0 = foreign "qfisom0" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let qfisominit = foreign "qfisominit" (gen @-> gen @-> gen @-> returning gen)
  let qfisominit0 = foreign "qfisominit0" (gen @-> gen @-> gen @-> returning gen)
  let qforbits = foreign "qforbits" (gen @-> gen @-> returning gen)
  let qfminimize = foreign "qfminimize" (gen @-> returning gen)
  let qfparam = foreign "qfparam" (gen @-> gen @-> long @-> returning gen)
  let qfsolve = foreign "qfsolve" (gen @-> returning gen)
  let z_isfundamental = foreign "Z_isfundamental" (gen @-> returning long)
  let classno = foreign "classno" (gen @-> returning gen)
  let classno2 = foreign "classno2" (gen @-> returning gen)

  let hclassnof_fact =
    foreign "hclassnoF_fact" (gen @-> gen @-> gen @-> returning gen)

  let hclassno = foreign "hclassno" (gen @-> returning gen)
  let hclassno6 = foreign "hclassno6" (gen @-> returning gen)
  let isfundamental = foreign "isfundamental" (gen @-> returning long)
  let qfbclassno0 = foreign "qfbclassno0" (gen @-> long @-> returning gen)
  let quadclassnof = foreign "quadclassnoF" (gen @-> ptr gen @-> returning gen)

  let quadclassnof_fact =
    foreign "quadclassnoF_fact" (gen @-> gen @-> gen @-> returning gen)

  let quaddisc = foreign "quaddisc" (gen @-> returning gen)
  let quadregulator = foreign "quadregulator" (gen @-> long @-> returning gen)
  let quadunit = foreign "quadunit" (gen @-> returning gen)
  let quadunit0 = foreign "quadunit0" (gen @-> long @-> returning gen)
  let quadunitindex = foreign "quadunitindex" (gen @-> gen @-> returning gen)
  let quadunitnorm = foreign "quadunitnorm" (gen @-> returning long)
  let sisfundamental = foreign "sisfundamental" (long @-> returning long)

  let uhclassnof_fact =
    foreign "uhclassnoF_fact" (gen @-> long @-> returning long)

  let unegisfundamental =
    foreign "unegisfundamental" (pari_ulong @-> returning long)

  let unegquadclassnof =
    foreign "unegquadclassnoF"
      (pari_ulong @-> ptr pari_ulong @-> returning pari_ulong)

  let uposisfundamental =
    foreign "uposisfundamental" (pari_ulong @-> returning long)

  let uposquadclassnof =
    foreign "uposquadclassnoF"
      (pari_ulong @-> ptr pari_ulong @-> returning pari_ulong)

  let uquadclassnof_fact =
    foreign "uquadclassnoF_fact"
      (pari_ulong @-> long @-> gen @-> gen @-> returning pari_ulong)

  let zn_quad_roots =
    foreign "Zn_quad_roots" (gen @-> gen @-> gen @-> returning gen)

  let genrand = foreign "genrand" (gen @-> returning gen)
  let getrand = foreign "getrand" (void @-> returning gen)
  let pari_rand = foreign "pari_rand" (void @-> returning pari_ulong)
  let randomi = foreign "randomi" (gen @-> returning gen)
  let randomr = foreign "randomr" (long @-> returning gen)
  let random_f2x = foreign "random_F2x" (long @-> long @-> returning gen)
  let random_fl = foreign "random_Fl" (pari_ulong @-> returning pari_ulong)
  let random_bits = foreign "random_bits" (long @-> returning long)
  let random_zv = foreign "random_zv" (long @-> returning gen)
  let setrand = foreign "setrand" (gen @-> returning void)

  let ellratpoints =
    foreign "ellratpoints" (gen @-> gen @-> long @-> returning gen)

  let hyperellratpoints =
    foreign "hyperellratpoints" (gen @-> gen @-> long @-> returning gen)

  let qx_complex_roots =
    foreign "QX_complex_roots" (gen @-> long @-> returning gen)

  let fft = foreign "FFT" (gen @-> gen @-> returning gen)
  let fftinv = foreign "FFTinv" (gen @-> gen @-> returning gen)
  let cleanroots = foreign "cleanroots" (gen @-> long @-> returning gen)
  let fujiwara_bound = foreign "fujiwara_bound" (gen @-> returning double)

  let fujiwara_bound_real =
    foreign "fujiwara_bound_real" (gen @-> long @-> returning double)

  let isrealappr = foreign "isrealappr" (gen @-> long @-> returning int)
  let polgraeffe = foreign "polgraeffe" (gen @-> returning gen)

  let polmod_to_embed =
    foreign "polmod_to_embed" (gen @-> long @-> returning gen)

  let polrootsbound = foreign "polrootsbound" (gen @-> gen @-> returning gen)
  let roots = foreign "roots" (gen @-> long @-> returning gen)
  let realroots = foreign "realroots" (gen @-> gen @-> long @-> returning gen)
  let zx_graeffe = foreign "ZX_graeffe" (gen @-> returning gen)

  let zx_realroots_irred =
    foreign "ZX_realroots_irred" (gen @-> long @-> returning gen)

  let zx_sturm = foreign "ZX_sturm" (gen @-> returning long)
  let zx_sturm_irred = foreign "ZX_sturm_irred" (gen @-> returning long)
  let zx_sturmpart = foreign "ZX_sturmpart" (gen @-> gen @-> returning long)

  let zx_uspensky =
    foreign "ZX_Uspensky" (gen @-> gen @-> long @-> long @-> returning gen)

  let factor_aurifeuille =
    foreign "factor_Aurifeuille" (gen @-> long @-> returning gen)

  let factor_aurifeuille_prime =
    foreign "factor_Aurifeuille_prime" (gen @-> long @-> returning gen)

  let galoissubcyclo =
    foreign "galoissubcyclo" (gen @-> gen @-> long @-> long @-> returning gen)

  let polsubcyclo =
    foreign "polsubcyclo" (long @-> long @-> long @-> returning gen)

  let polsubcyclofast =
    foreign "polsubcyclofast" (gen @-> long @-> long @-> long @-> returning gen)

  let znsubgroupgenerators =
    foreign "znsubgroupgenerators" (gen @-> long @-> returning gen)

  let nfsubfields = foreign "nfsubfields" (gen @-> long @-> returning gen)

  let nfsubfields0 =
    foreign "nfsubfields0" (gen @-> long @-> long @-> returning gen)

  let nfsubfieldscm = foreign "nfsubfieldscm" (gen @-> long @-> returning gen)
  let nfsubfieldsmax = foreign "nfsubfieldsmax" (gen @-> long @-> returning gen)
  let nflist = foreign "nflist" (gen @-> gen @-> long @-> gen @-> returning gen)
  let nfresolvent = foreign "nfresolvent" (gen @-> long @-> returning gen)
  let subgrouplist = foreign "subgrouplist" (gen @-> gen @-> returning gen)

  let forsubgroup =
    foreign "forsubgroup"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> gen @-> gen @-> returning void)

  let abmap_kernel = foreign "abmap_kernel" (gen @-> returning gen)

  let abmap_subgroup_image =
    foreign "abmap_subgroup_image" (gen @-> gen @-> returning gen)

  let bnrl1 = foreign "bnrL1" (gen @-> gen @-> long @-> long @-> returning gen)

  let bnrrootnumber =
    foreign "bnrrootnumber" (gen @-> gen @-> long @-> long @-> returning gen)

  let bnrstark = foreign "bnrstark" (gen @-> gen @-> long @-> returning gen)
  let cyc2elts = foreign "cyc2elts" (gen @-> returning gen)
  let qfbforms = foreign "qfbforms" (gen @-> returning gen)
  let quadhilbert = foreign "quadhilbert" (gen @-> long @-> returning gen)
  let quadray = foreign "quadray" (gen @-> gen @-> long @-> returning gen)
  let chartogenstr = foreign "chartoGENstr" (char @-> returning gen)
  let pari_strdup = foreign "pari_strdup" (string @-> returning string)

  let pari_strndup =
    foreign "pari_strndup" (string @-> long @-> returning string)

  let stack_strcat =
    foreign "stack_strcat" (string @-> string @-> returning string)

  let stack_strdup = foreign "stack_strdup" (string @-> returning string)
  let pari_strchr = foreign "pari_strchr" (gen @-> returning gen)
  let strjoin = foreign "strjoin" (gen @-> gen @-> returning gen)
  let strntogenstr = foreign "strntoGENstr" (string @-> long @-> returning gen)
end

module F49 (F : Ctypes.FOREIGN) = struct
  open F

  let strsplit = foreign "strsplit" (gen @-> gen @-> returning gen)
  let strtogenstr = foreign "strtoGENstr" (string @-> returning gen)
  let type_name = foreign "type_name" (long @-> returning string)
  let type0 = foreign "type0" (gen @-> returning gen)

  let asympnum =
    foreign "asympnum"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let asympnumraw =
    foreign "asympnumraw"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> long @-> gen @-> long @-> returning gen)

  let derivnum =
    foreign "derivnum"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let derivnumk =
    foreign "derivnumk"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> gen @-> gen @-> long @-> returning gen)

  let derivfun =
    foreign "derivfun"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let derivfunk =
    foreign "derivfunk"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> gen @-> gen @-> long @-> returning gen)

  let forvec_init =
    foreign "forvec_init" (ptr forvec_t @-> gen @-> long @-> returning int)

  let forvec_next = foreign "forvec_next" (ptr forvec_t @-> returning gen)

  let laurentseries =
    foreign "laurentseries"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> long @-> long @-> long @-> returning gen)

  let limitnum =
    foreign "limitnum"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let polzag = foreign "polzag" (long @-> long @-> returning gen)

  let prodeuler =
    foreign "prodeuler"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> returning gen)

  let prodinf =
    foreign "prodinf"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let prodinf1 =
    foreign "prodinf1"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let solvestep =
    foreign "solvestep"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> gen @-> long @-> long @-> returning gen)

  let sumalt =
    foreign "sumalt"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let sumalt2 =
    foreign "sumalt2"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let sumpos =
    foreign "sumpos"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let sumpos2 =
    foreign "sumpos2"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let suminf =
    foreign "suminf"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let suminf_bitprec =
    foreign "suminf_bitprec"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let sumdivmultexpr =
    foreign "sumdivmultexpr"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> returning gen)

  let zbrent =
    foreign "zbrent"
      (ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning gen)
      @-> gen @-> gen @-> long @-> returning gen)

  let bnfisintnorm = foreign "bnfisintnorm" (gen @-> gen @-> returning gen)
  let bnfisintnormabs = foreign "bnfisintnormabs" (gen @-> gen @-> returning gen)
  let ideals_by_norm = foreign "ideals_by_norm" (gen @-> gen @-> returning gen)
  let thue = foreign "thue" (gen @-> gen @-> gen @-> returning gen)
  let thueinit = foreign "thueinit" (gen @-> long @-> long @-> returning gen)
  let pi2n = foreign "Pi2n" (long @-> long @-> returning gen)
  let pii2 = foreign "PiI2" (long @-> returning gen)
  let pii2n = foreign "PiI2n" (long @-> long @-> returning gen)
  let qp_exp = foreign "Qp_exp" (gen @-> returning gen)
  let qp_exp_prec = foreign "Qp_exp_prec" (gen @-> returning long)
  let qp_log = foreign "Qp_log" (gen @-> returning gen)
  let qp_sqrt = foreign "Qp_sqrt" (gen @-> returning gen)
  let qp_sqrtn = foreign "Qp_sqrtn" (gen @-> gen @-> ptr gen @-> returning gen)
  let zn_sqrt = foreign "Zn_sqrt" (gen @-> gen @-> returning gen)

  let zp_teichmuller =
    foreign "Zp_teichmuller" (gen @-> gen @-> long @-> gen @-> returning gen)

  let agm = foreign "agm" (gen @-> gen @-> long @-> returning gen)
  let constcatalan = foreign "constcatalan" (long @-> returning gen)
  let consteuler = foreign "consteuler" (long @-> returning gen)
  let constlog2 = foreign "constlog2" (long @-> returning gen)
  let constpi = foreign "constpi" (long @-> returning gen)
  let cxexpm1 = foreign "cxexpm1" (gen @-> long @-> returning gen)
  let elle = foreign "ellE" (gen @-> long @-> returning gen)
  let ellk = foreign "ellK" (gen @-> long @-> returning gen)
  let expir = foreign "expIr" (gen @-> returning gen)
  let exp1r_abs = foreign "exp1r_abs" (gen @-> returning gen)
  let gcos = foreign "gcos" (gen @-> long @-> returning gen)
  let gcotan = foreign "gcotan" (gen @-> long @-> returning gen)
  let gcotanh = foreign "gcotanh" (gen @-> long @-> returning gen)
  let gexp = foreign "gexp" (gen @-> long @-> returning gen)
  let gexpm1 = foreign "gexpm1" (gen @-> long @-> returning gen)
  let glog = foreign "glog" (gen @-> long @-> returning gen)
  let glog1p = foreign "glog1p" (gen @-> long @-> returning gen)
  let gpow = foreign "gpow" (gen @-> gen @-> long @-> returning gen)
  let gpowers = foreign "gpowers" (gen @-> long @-> returning gen)
  let gpowers0 = foreign "gpowers0" (gen @-> long @-> gen @-> returning gen)
  let gpowgs = foreign "gpowgs" (gen @-> long @-> returning gen)
  let grootsof1 = foreign "grootsof1" (long @-> long @-> returning gen)
  let gsin = foreign "gsin" (gen @-> long @-> returning gen)
  let gsinc = foreign "gsinc" (gen @-> long @-> returning gen)

  let gsincos =
    foreign "gsincos" (gen @-> ptr gen @-> ptr gen @-> long @-> returning void)

  let gsqrpowers = foreign "gsqrpowers" (gen @-> long @-> returning gen)
  let gsqrt = foreign "gsqrt" (gen @-> long @-> returning gen)

  let gsqrtn =
    foreign "gsqrtn" (gen @-> gen @-> ptr gen @-> long @-> returning gen)

  let gtan = foreign "gtan" (gen @-> long @-> returning gen)
  let logr_abs = foreign "logr_abs" (gen @-> returning gen)
  let mpcos = foreign "mpcos" (gen @-> returning gen)
  let mpeuler = foreign "mpeuler" (long @-> returning gen)
  let mpcatalan = foreign "mpcatalan" (long @-> returning gen)

  let mpsincosm1 =
    foreign "mpsincosm1" (gen @-> ptr gen @-> ptr gen @-> returning void)

  let mpexp = foreign "mpexp" (gen @-> returning gen)
  let mpexpm1 = foreign "mpexpm1" (gen @-> returning gen)
  let mplog = foreign "mplog" (gen @-> returning gen)
  let mplog2 = foreign "mplog2" (long @-> returning gen)
  let mppi = foreign "mppi" (long @-> returning gen)
  let mpsin = foreign "mpsin" (gen @-> returning gen)

  let mpsincos =
    foreign "mpsincos" (gen @-> ptr gen @-> ptr gen @-> returning void)

  let pow2pis = foreign "pow2Pis" (gen @-> long @-> returning gen)
  let powpis = foreign "powPis" (gen @-> long @-> returning gen)
  let powcx = foreign "powcx" (gen @-> gen @-> gen @-> long @-> returning gen)

  let powcx_prec =
    foreign "powcx_prec" (long @-> gen @-> long @-> returning long)

  let powersr = foreign "powersr" (gen @-> long @-> returning gen)
  let powis = foreign "powis" (gen @-> long @-> returning gen)
  let powiu = foreign "powiu" (gen @-> pari_ulong @-> returning gen)
  let powrfrac = foreign "powrfrac" (gen @-> long @-> long @-> returning gen)
  let powrs = foreign "powrs" (gen @-> long @-> returning gen)
  let powrshalf = foreign "powrshalf" (gen @-> long @-> returning gen)
  let powru = foreign "powru" (gen @-> pari_ulong @-> returning gen)
  let powruhalf = foreign "powruhalf" (gen @-> pari_ulong @-> returning gen)
  let powuu = foreign "powuu" (pari_ulong @-> pari_ulong @-> returning gen)
  let powgi = foreign "powgi" (gen @-> gen @-> returning gen)
  let rootsof1_cx = foreign "rootsof1_cx" (gen @-> long @-> returning gen)

  let rootsof1u_cx =
    foreign "rootsof1u_cx" (pari_ulong @-> long @-> returning gen)

  let rootsof1q_cx =
    foreign "rootsof1q_cx" (long @-> long @-> long @-> returning gen)
end

module F50 (F : Ctypes.FOREIGN) = struct
  open F

  let rootsof1powinit =
    foreign "rootsof1powinit" (long @-> long @-> long @-> returning gen)

  let rootsof1pow = foreign "rootsof1pow" (gen @-> long @-> returning gen)
  let serchop = foreign "serchop" (gen @-> long @-> returning gen)
  let serchop_i = foreign "serchop_i" (gen @-> long @-> returning gen)
  let serchop0 = foreign "serchop0" (gen @-> returning gen)
  let sqrtnint = foreign "sqrtnint" (gen @-> long @-> returning gen)
  let sqrtnr_abs = foreign "sqrtnr_abs" (gen @-> long @-> returning gen)
  let teich = foreign "teich" (gen @-> returning gen)

  let teichmullerinit =
    foreign "teichmullerinit" (long @-> long @-> returning gen)

  let teichmuller = foreign "teichmuller" (gen @-> gen @-> returning gen)

  let trans_eval =
    foreign "trans_eval"
      (string
      @-> static_funptr Ctypes.(gen @-> long @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let trans_evalgen =
    foreign "trans_evalgen"
      (string @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> long @-> returning gen)
      @-> gen @-> long @-> returning gen)

  let upowuu =
    foreign "upowuu" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let upowers = foreign "upowers" (pari_ulong @-> long @-> returning gen)

  let usqrtn =
    foreign "usqrtn" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let usqrt = foreign "usqrt" (pari_ulong @-> returning pari_ulong)
  let qp_gamma = foreign "Qp_gamma" (gen @-> returning gen)

  let atanhuu =
    foreign "atanhuu" (pari_ulong @-> pari_ulong @-> long @-> returning gen)

  let atanhui = foreign "atanhui" (pari_ulong @-> gen @-> long @-> returning gen)
  let gacosh = foreign "gacosh" (gen @-> long @-> returning gen)
  let gacos = foreign "gacos" (gen @-> long @-> returning gen)
  let garg = foreign "garg" (gen @-> long @-> returning gen)
  let gasinh = foreign "gasinh" (gen @-> long @-> returning gen)
  let gasin = foreign "gasin" (gen @-> long @-> returning gen)
  let gatan = foreign "gatan" (gen @-> long @-> returning gen)
  let gatanh = foreign "gatanh" (gen @-> long @-> returning gen)
  let gcosh = foreign "gcosh" (gen @-> long @-> returning gen)
  let ggammah = foreign "ggammah" (gen @-> long @-> returning gen)
  let ggamma = foreign "ggamma" (gen @-> long @-> returning gen)
  let ggamma1m1 = foreign "ggamma1m1" (gen @-> long @-> returning gen)
  let glngamma = foreign "glngamma" (gen @-> long @-> returning gen)
  let gpsi = foreign "gpsi" (gen @-> long @-> returning gen)
  let gsinh = foreign "gsinh" (gen @-> long @-> returning gen)
  let gtanh = foreign "gtanh" (gen @-> long @-> returning gen)
  let mpfactr = foreign "mpfactr" (long @-> long @-> returning gen)

  let mpsinhcosh =
    foreign "mpsinhcosh" (gen @-> ptr gen @-> ptr gen @-> returning void)

  let psi1series =
    foreign "psi1series" (long @-> long @-> long @-> returning gen)

  let sumformal = foreign "sumformal" (gen @-> long @-> returning gen)

  let rgv_is_arithprog =
    foreign "RgV_is_arithprog" (gen @-> ptr gen @-> ptr gen @-> returning int)

  let besseljzero =
    foreign "besseljzero" (gen @-> long @-> long @-> returning gen)

  let besselyzero =
    foreign "besselyzero" (gen @-> long @-> long @-> returning gen)

  let constzeta = foreign "constzeta" (long @-> long @-> returning gen)
  let cxek = foreign "cxEk" (gen @-> long @-> long @-> returning gen)
  let dblmodulus = foreign "dblmodulus" (gen @-> returning double)
  let dilog = foreign "dilog" (gen @-> long @-> returning gen)
  let eint1 = foreign "eint1" (gen @-> long @-> returning gen)
  let expipir = foreign "expIPiR" (gen @-> long @-> returning gen)
  let expipic = foreign "expIPiC" (gen @-> long @-> returning gen)
  let expixy = foreign "expIxy" (gen @-> gen @-> long @-> returning gen)
  let eta = foreign "eta" (gen @-> long @-> returning gen)
  let eta0 = foreign "eta0" (gen @-> long @-> long @-> returning gen)
  let gerfc = foreign "gerfc" (gen @-> long @-> returning gen)
  let gpolylog = foreign "gpolylog" (long @-> gen @-> long @-> returning gen)
  let gzeta = foreign "gzeta" (gen @-> long @-> returning gen)
  let hbessel1 = foreign "hbessel1" (gen @-> gen @-> long @-> returning gen)
  let hbessel2 = foreign "hbessel2" (gen @-> gen @-> long @-> returning gen)
  let hyperu = foreign "hyperu" (gen @-> gen @-> gen @-> long @-> returning gen)
  let ibessel = foreign "ibessel" (gen @-> gen @-> long @-> returning gen)
  let incgam = foreign "incgam" (gen @-> gen @-> long @-> returning gen)

  let incgam0 =
    foreign "incgam0" (gen @-> gen @-> gen @-> long @-> returning gen)

  let incgamc = foreign "incgamc" (gen @-> gen @-> long @-> returning gen)
  let jbessel = foreign "jbessel" (gen @-> gen @-> long @-> returning gen)
  let jbesselh = foreign "jbesselh" (gen @-> gen @-> long @-> returning gen)
  let jell = foreign "jell" (gen @-> long @-> returning gen)
  let kbessel = foreign "kbessel" (gen @-> gen @-> long @-> returning gen)
  let mpeint1 = foreign "mpeint1" (gen @-> gen @-> returning gen)
  let mpveceint1 = foreign "mpveceint1" (gen @-> gen @-> long @-> returning gen)

  let polylog0 =
    foreign "polylog0" (long @-> gen @-> long @-> long @-> returning gen)

  let sumdedekind = foreign "sumdedekind" (gen @-> gen @-> returning gen)

  let sumdedekind_coprime =
    foreign "sumdedekind_coprime" (gen @-> gen @-> returning gen)

  let szeta = foreign "szeta" (long @-> long @-> returning gen)
  let theta = foreign "theta" (gen @-> gen @-> long @-> returning gen)
  let thetanullk = foreign "thetanullk" (gen @-> long @-> long @-> returning gen)
  let trueeta = foreign "trueeta" (gen @-> long @-> returning gen)

  let u_sumdedekind_coprime =
    foreign "u_sumdedekind_coprime" (long @-> long @-> returning gen)

  let upper_to_cx = foreign "upper_to_cx" (gen @-> ptr long @-> returning gen)
  let veceint1 = foreign "veceint1" (gen @-> gen @-> long @-> returning gen)

  let vecthetanullk =
    foreign "vecthetanullk" (gen @-> long @-> long @-> returning gen)

  let vecthetanullk_tau =
    foreign "vecthetanullk_tau" (gen @-> long @-> long @-> returning gen)

  let veczeta =
    foreign "veczeta" (gen @-> gen @-> long @-> long @-> returning gen)

  let weber0 = foreign "weber0" (gen @-> long @-> long @-> returning gen)
  let weberf = foreign "weberf" (gen @-> long @-> returning gen)
  let weberf1 = foreign "weberf1" (gen @-> long @-> returning gen)
  let weberf2 = foreign "weberf2" (gen @-> long @-> returning gen)
  let ybessel = foreign "ybessel" (gen @-> gen @-> long @-> returning gen)
  let sl2_inv_shallow = foreign "SL2_inv_shallow" (gen @-> returning gen)
  let qevproj_apply = foreign "Qevproj_apply" (gen @-> gen @-> returning gen)

  let qevproj_apply_vecei =
    foreign "Qevproj_apply_vecei" (gen @-> gen @-> long @-> returning gen)

  let qevproj_down = foreign "Qevproj_down" (gen @-> gen @-> returning gen)
  let qevproj_init = foreign "Qevproj_init" (gen @-> returning gen)
  let rgx_act_gl2q = foreign "RgX_act_Gl2Q" (gen @-> long @-> returning gen)
  let rgx_act_zgl2q = foreign "RgX_act_ZGl2Q" (gen @-> long @-> returning gen)
  let checkms = foreign "checkms" (gen @-> returning void)
  let checkmspadic = foreign "checkmspadic" (gen @-> returning void)

  let ellpadiclambdamu =
    foreign "ellpadiclambdamu" (gen @-> long @-> long @-> long @-> returning gen)

  let mfnumcusps = foreign "mfnumcusps" (gen @-> returning gen)
  let mfnumcusps_fact = foreign "mfnumcusps_fact" (gen @-> returning gen)

  let mfnumcuspsu_fact =
    foreign "mfnumcuspsu_fact" (gen @-> returning pari_ulong)

  let mfnumcuspsu = foreign "mfnumcuspsu" (pari_ulong @-> returning pari_ulong)
  let msfromcusp = foreign "msfromcusp" (gen @-> gen @-> returning gen)
end

module F51 (F : Ctypes.FOREIGN) = struct
  open F

  let msfromell = foreign "msfromell" (gen @-> long @-> returning gen)
  let msfromhecke = foreign "msfromhecke" (gen @-> gen @-> gen @-> returning gen)
  let msdim = foreign "msdim" (gen @-> returning long)
  let mseval2_ooq = foreign "mseval2_ooQ" (gen @-> gen @-> gen @-> returning gen)
  let msgetlevel = foreign "msgetlevel" (gen @-> returning long)
  let msgetsign = foreign "msgetsign" (gen @-> returning long)
  let msgetweight = foreign "msgetweight" (gen @-> returning long)

  let msatkinlehner =
    foreign "msatkinlehner" (gen @-> long @-> gen @-> returning gen)

  let mscuspidal = foreign "mscuspidal" (gen @-> long @-> returning gen)
  let mseisenstein = foreign "mseisenstein" (gen @-> returning gen)
  let mseval = foreign "mseval" (gen @-> gen @-> gen @-> returning gen)
  let mshecke = foreign "mshecke" (gen @-> long @-> gen @-> returning gen)
  let msinit = foreign "msinit" (gen @-> gen @-> long @-> returning gen)
  let msissymbol = foreign "msissymbol" (gen @-> gen @-> returning gen)
  let mslattice = foreign "mslattice" (gen @-> gen @-> returning gen)
  let msomseval = foreign "msomseval" (gen @-> gen @-> gen @-> returning gen)

  let mspadic_parse_chi =
    foreign "mspadic_parse_chi" (gen @-> ptr gen @-> ptr gen @-> returning void)

  let mspadic_unit_eigenvalue =
    foreign "mspadic_unit_eigenvalue"
      (gen @-> long @-> gen @-> long @-> returning gen)

  let mspadicinit =
    foreign "mspadicinit" (gen @-> long @-> long @-> long @-> returning gen)

  let mspadicl = foreign "mspadicL" (gen @-> gen @-> long @-> returning gen)

  let mspadicmoments =
    foreign "mspadicmoments" (gen @-> gen @-> long @-> returning gen)

  let mspadicseries = foreign "mspadicseries" (gen @-> long @-> returning gen)
  let mspathgens = foreign "mspathgens" (gen @-> returning gen)
  let mspathlog = foreign "mspathlog" (gen @-> gen @-> returning gen)
  let msnew = foreign "msnew" (gen @-> returning gen)
  let mspetersson = foreign "mspetersson" (gen @-> gen @-> gen @-> returning gen)
  let mspolygon = foreign "mspolygon" (gen @-> long @-> returning gen)
  let msstar = foreign "msstar" (gen @-> gen @-> returning gen)

  let msqexpansion =
    foreign "msqexpansion" (gen @-> gen @-> pari_ulong @-> returning gen)

  let mssplit = foreign "mssplit" (gen @-> gen @-> long @-> returning gen)
  let mstooms = foreign "mstooms" (gen @-> gen @-> returning gen)
  let mscosets0 = foreign "mscosets0" (gen @-> gen @-> returning gen)

  let mscosets =
    foreign "mscosets"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> returning gen)

  let msfarey =
    foreign "msfarey"
      (gen @-> ptr void
      @-> static_funptr Ctypes.(ptr void @-> gen @-> returning long)
      @-> ptr gen @-> returning gen)

  let msfarey0 = foreign "msfarey0" (gen @-> gen @-> ptr gen @-> returning gen)
  let checkfarey_i = foreign "checkfarey_i" (gen @-> returning int)

  let polylogmult =
    foreign "polylogmult" (gen @-> gen @-> long @-> returning gen)

  let polylogmult_interpolate =
    foreign "polylogmult_interpolate"
      (gen @-> gen @-> gen @-> long @-> returning gen)

  let zetamult = foreign "zetamult" (gen @-> long @-> returning gen)
  let zetamultdual = foreign "zetamultdual" (gen @-> returning gen)

  let zetamult_interpolate =
    foreign "zetamult_interpolate" (gen @-> gen @-> long @-> returning gen)

  let zetamultall =
    foreign "zetamultall" (long @-> long @-> long @-> returning gen)

  let zetamultconvert =
    foreign "zetamultconvert" (gen @-> long @-> returning gen)

  let fl_add =
    foreign "Fl_add"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_addmul_pre =
    foreign "Fl_addmul_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fl_addmulmul_pre =
    foreign "Fl_addmulmul_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> pari_ulong @-> returning pari_ulong)

  let fl_center =
    foreign "Fl_center"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning long)

  let fl_div =
    foreign "Fl_div"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_double =
    foreign "Fl_double" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_ellj_pre =
    foreign "Fl_ellj_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fl_halve =
    foreign "Fl_halve" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_mul =
    foreign "Fl_mul"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_mul_pre =
    foreign "Fl_mul_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let fl_neg =
    foreign "Fl_neg" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_sqr =
    foreign "Fl_sqr" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_sqr_pre =
    foreign "Fl_sqr_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_sub =
    foreign "Fl_sub"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let fl_triple =
    foreign "Fl_triple" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let mod2 = foreign "Mod2" (gen @-> returning pari_ulong)
  let mod4 = foreign "Mod4" (gen @-> returning pari_ulong)
  let mod8 = foreign "Mod8" (gen @-> returning pari_ulong)
  let mod16 = foreign "Mod16" (gen @-> returning pari_ulong)
  let mod32 = foreign "Mod32" (gen @-> returning pari_ulong)
  let mod64 = foreign "Mod64" (gen @-> returning pari_ulong)
  let abscmpiu = foreign "abscmpiu" (gen @-> pari_ulong @-> returning int)
  let abscmpui = foreign "abscmpui" (pari_ulong @-> gen @-> returning int)
  let absequaliu = foreign "absequaliu" (gen @-> pari_ulong @-> returning int)
  let absi = foreign "absi" (gen @-> returning gen)
  let absi_shallow = foreign "absi_shallow" (gen @-> returning gen)
  let absr = foreign "absr" (gen @-> returning gen)
  let absrnz_equal1 = foreign "absrnz_equal1" (gen @-> returning int)
  let absrnz_equal2n = foreign "absrnz_equal2n" (gen @-> returning int)
  let addii = foreign "addii" (gen @-> gen @-> returning gen)
  let addiiz = foreign "addiiz" (gen @-> gen @-> gen @-> returning void)
  let addir = foreign "addir" (gen @-> gen @-> returning gen)
  let addirz = foreign "addirz" (gen @-> gen @-> gen @-> returning void)
  let addis = foreign "addis" (gen @-> long @-> returning gen)
  let addri = foreign "addri" (gen @-> gen @-> returning gen)
  let addriz = foreign "addriz" (gen @-> gen @-> gen @-> returning void)
  let addrr = foreign "addrr" (gen @-> gen @-> returning gen)
  let addrrz = foreign "addrrz" (gen @-> gen @-> gen @-> returning void)
  let addrs = foreign "addrs" (gen @-> long @-> returning gen)
  let addsi = foreign "addsi" (long @-> gen @-> returning gen)
  let addsiz = foreign "addsiz" (long @-> gen @-> gen @-> returning void)
  let addsrz = foreign "addsrz" (long @-> gen @-> gen @-> returning void)
  let addss = foreign "addss" (long @-> long @-> returning gen)
  let addssz = foreign "addssz" (long @-> long @-> gen @-> returning void)
  let adduu = foreign "adduu" (pari_ulong @-> pari_ulong @-> returning gen)
  let affgr = foreign "affgr" (gen @-> gen @-> returning void)
  let affii = foreign "affii" (gen @-> gen @-> returning void)
  let affiz = foreign "affiz" (gen @-> gen @-> returning void)
  let affrr_fixlg = foreign "affrr_fixlg" (gen @-> gen @-> returning void)
  let affsi = foreign "affsi" (long @-> gen @-> returning void)
  let affsr = foreign "affsr" (long @-> gen @-> returning void)
  let affsz = foreign "affsz" (long @-> gen @-> returning void)
  let affui = foreign "affui" (pari_ulong @-> gen @-> returning void)
  let affur = foreign "affur" (pari_ulong @-> gen @-> returning void)
  let cgetg = foreign "cgetg" (long @-> long @-> returning gen)
  let cgetg_block = foreign "cgetg_block" (long @-> long @-> returning gen)
  let cgetg_copy = foreign "cgetg_copy" (gen @-> ptr long @-> returning gen)
end

module F52 (F : Ctypes.FOREIGN) = struct
  open F

  let cgeti = foreign "cgeti" (long @-> returning gen)
  let cgetineg = foreign "cgetineg" (long @-> returning gen)
  let cgetipos = foreign "cgetipos" (long @-> returning gen)
  let cgetr = foreign "cgetr" (long @-> returning gen)
  let cgetr_block = foreign "cgetr_block" (long @-> returning gen)
  let cmpir = foreign "cmpir" (gen @-> gen @-> returning int)
  let cmpis = foreign "cmpis" (gen @-> long @-> returning int)
  let cmpiu = foreign "cmpiu" (gen @-> pari_ulong @-> returning int)
  let cmpri = foreign "cmpri" (gen @-> gen @-> returning int)
  let cmprs = foreign "cmprs" (gen @-> long @-> returning int)
  let cmpsi = foreign "cmpsi" (long @-> gen @-> returning int)
  let cmpsr = foreign "cmpsr" (long @-> gen @-> returning int)
  let cmpss = foreign "cmpss" (long @-> long @-> returning int)
  let cmpui = foreign "cmpui" (pari_ulong @-> gen @-> returning int)
  let cmpuu = foreign "cmpuu" (pari_ulong @-> pari_ulong @-> returning int)
  let cxtofp = foreign "cxtofp" (gen @-> long @-> returning gen)
  let divii = foreign "divii" (gen @-> gen @-> returning gen)
  let diviiz = foreign "diviiz" (gen @-> gen @-> gen @-> returning void)
  let divirz = foreign "divirz" (gen @-> gen @-> gen @-> returning void)
  let divisz = foreign "divisz" (gen @-> long @-> gen @-> returning void)
  let divriz = foreign "divriz" (gen @-> gen @-> gen @-> returning void)
  let divrrz = foreign "divrrz" (gen @-> gen @-> gen @-> returning void)
  let divrsz = foreign "divrsz" (gen @-> long @-> gen @-> returning void)

  let divsi_rem =
    foreign "divsi_rem" (long @-> gen @-> ptr long @-> returning gen)

  let divsiz = foreign "divsiz" (long @-> gen @-> gen @-> returning void)
  let divsrz = foreign "divsrz" (long @-> gen @-> gen @-> returning void)
  let divss = foreign "divss" (long @-> long @-> returning gen)

  let divss_rem =
    foreign "divss_rem" (long @-> long @-> ptr long @-> returning gen)

  let divssz = foreign "divssz" (long @-> long @-> gen @-> returning void)
  let dvdii = foreign "dvdii" (gen @-> gen @-> returning int)
  let dvdiiz = foreign "dvdiiz" (gen @-> gen @-> gen @-> returning int)
  let dvdis = foreign "dvdis" (gen @-> long @-> returning int)
  let dvdisz = foreign "dvdisz" (gen @-> long @-> gen @-> returning int)
  let dvdiu = foreign "dvdiu" (gen @-> pari_ulong @-> returning int)
  let dvdiuz = foreign "dvdiuz" (gen @-> pari_ulong @-> gen @-> returning int)
  let dvdsi = foreign "dvdsi" (long @-> gen @-> returning int)
  let dvdui = foreign "dvdui" (pari_ulong @-> gen @-> returning int)

  let dvmdiiz =
    foreign "dvmdiiz" (gen @-> gen @-> gen @-> gen @-> returning void)

  let dvmdis = foreign "dvmdis" (gen @-> long @-> ptr gen @-> returning gen)

  let dvmdisz =
    foreign "dvmdisz" (gen @-> long @-> gen @-> gen @-> returning void)

  let dvmdsbil = foreign "dvmdsBIL" (long @-> ptr long @-> returning long)
  let dvmdsi = foreign "dvmdsi" (long @-> gen @-> ptr gen @-> returning gen)

  let dvmdsiz =
    foreign "dvmdsiz" (long @-> gen @-> gen @-> gen @-> returning void)

  let dvmdss = foreign "dvmdss" (long @-> long @-> ptr gen @-> returning gen)

  let dvmdssz =
    foreign "dvmdssz" (long @-> long @-> gen @-> gen @-> returning void)

  let dvmdubil =
    foreign "dvmduBIL" (pari_ulong @-> ptr pari_ulong @-> returning pari_ulong)

  let equalis = foreign "equalis" (gen @-> long @-> returning int)
  let equalsi = foreign "equalsi" (long @-> gen @-> returning int)
  let equalui = foreign "equalui" (pari_ulong @-> gen @-> returning int)
  let equaliu = foreign "equaliu" (gen @-> pari_ulong @-> returning int)
  let absequalui = foreign "absequalui" (pari_ulong @-> gen @-> returning int)

  let ceildivuu =
    foreign "ceildivuu" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let evalexpo = foreign "evalexpo" (long @-> returning long)
  let evallg = foreign "evallg" (long @-> returning long)
  let evalprecp = foreign "evalprecp" (long @-> returning long)
  let evalvalp = foreign "evalvalp" (long @-> returning long)
  let evalvalser = foreign "evalvalser" (long @-> returning long)
  let expi = foreign "expi" (gen @-> returning long)
  let expu = foreign "expu" (pari_ulong @-> returning long)
  let fixlg = foreign "fixlg" (gen @-> long @-> returning void)
  let fractor = foreign "fractor" (gen @-> long @-> returning gen)
  let gc_bool = foreign "gc_bool" (pari_sp @-> int @-> returning int)
  let gc_const = foreign "gc_const" (pari_sp @-> gen @-> returning gen)
  let gc_double = foreign "gc_double" (pari_sp @-> double @-> returning double)
  let gc_int = foreign "gc_int" (pari_sp @-> int @-> returning int)
  let gc_long = foreign "gc_long" (pari_sp @-> long @-> returning long)

  let gc_ulong =
    foreign "gc_ulong" (pari_sp @-> pari_ulong @-> returning pari_ulong)

  let gc_null = foreign "gc_NULL" (pari_sp @-> returning gen)
  let icopy = foreign "icopy" (gen @-> returning gen)
  let icopyspec = foreign "icopyspec" (gen @-> long @-> returning gen)
  let icopy_avma = foreign "icopy_avma" (gen @-> pari_sp @-> returning gen)
  let int_bit = foreign "int_bit" (gen @-> long @-> returning pari_ulong)
  let itor = foreign "itor" (gen @-> long @-> returning gen)
  let itos = foreign "itos" (gen @-> returning long)
  let itos_or_0 = foreign "itos_or_0" (gen @-> returning long)
  let itou = foreign "itou" (gen @-> returning pari_ulong)
  let itou_or_0 = foreign "itou_or_0" (gen @-> returning pari_ulong)
  let leafcopy = foreign "leafcopy" (gen @-> returning gen)
  let leafcopy_avma = foreign "leafcopy_avma" (gen @-> pari_sp @-> returning gen)
  let maxdd = foreign "maxdd" (double @-> double @-> returning double)
  let maxss = foreign "maxss" (long @-> long @-> returning long)
  let maxuu = foreign "maxuu" (pari_ulong @-> pari_ulong @-> returning long)
  let mindd = foreign "mindd" (double @-> double @-> returning double)
  let minss = foreign "minss" (long @-> long @-> returning long)
  let minuu = foreign "minuu" (pari_ulong @-> pari_ulong @-> returning long)
  let mod16 = foreign "mod16" (gen @-> returning long)
  let mod2 = foreign "mod2" (gen @-> returning long)
  let mod2bil = foreign "mod2BIL" (gen @-> returning pari_ulong)
  let mod32 = foreign "mod32" (gen @-> returning long)
  let mod4 = foreign "mod4" (gen @-> returning long)
  let mod64 = foreign "mod64" (gen @-> returning long)
  let mod8 = foreign "mod8" (gen @-> returning long)
  let modis = foreign "modis" (gen @-> long @-> returning gen)
  let modisz = foreign "modisz" (gen @-> long @-> gen @-> returning void)
  let modsi = foreign "modsi" (long @-> gen @-> returning gen)
  let modsiz = foreign "modsiz" (long @-> gen @-> gen @-> returning void)
  let modss = foreign "modss" (long @-> long @-> returning gen)
  let modssz = foreign "modssz" (long @-> long @-> gen @-> returning void)
  let mpabs = foreign "mpabs" (gen @-> returning gen)
  let mpabs_shallow = foreign "mpabs_shallow" (gen @-> returning gen)
end

module F53 (F : Ctypes.FOREIGN) = struct
  open F

  let mpadd = foreign "mpadd" (gen @-> gen @-> returning gen)
  let mpaddz = foreign "mpaddz" (gen @-> gen @-> gen @-> returning void)
  let mpaff = foreign "mpaff" (gen @-> gen @-> returning void)
  let mpceil = foreign "mpceil" (gen @-> returning gen)
  let mpcmp = foreign "mpcmp" (gen @-> gen @-> returning int)
  let mpcopy = foreign "mpcopy" (gen @-> returning gen)
  let mpdiv = foreign "mpdiv" (gen @-> gen @-> returning gen)
  let mpexpo = foreign "mpexpo" (gen @-> returning long)
  let mpfloor = foreign "mpfloor" (gen @-> returning gen)
  let mpmul = foreign "mpmul" (gen @-> gen @-> returning gen)
  let mpmulz = foreign "mpmulz" (gen @-> gen @-> gen @-> returning void)
  let mpneg = foreign "mpneg" (gen @-> returning gen)
  let mpodd = foreign "mpodd" (gen @-> returning int)
  let mpround = foreign "mpround" (gen @-> returning gen)
  let mpshift = foreign "mpshift" (gen @-> long @-> returning gen)
  let mpsqr = foreign "mpsqr" (gen @-> returning gen)
  let mpsub = foreign "mpsub" (gen @-> gen @-> returning gen)
  let mpsubz = foreign "mpsubz" (gen @-> gen @-> gen @-> returning void)
  let mptrunc = foreign "mptrunc" (gen @-> returning gen)
  let muliiz = foreign "muliiz" (gen @-> gen @-> gen @-> returning void)
  let mulirz = foreign "mulirz" (gen @-> gen @-> gen @-> returning void)
  let mulis = foreign "mulis" (gen @-> long @-> returning gen)
  let muliu = foreign "muliu" (gen @-> pari_ulong @-> returning gen)
  let mulri = foreign "mulri" (gen @-> gen @-> returning gen)
  let mulriz = foreign "mulriz" (gen @-> gen @-> gen @-> returning void)
  let mulrrz = foreign "mulrrz" (gen @-> gen @-> gen @-> returning void)
  let mulrs = foreign "mulrs" (gen @-> long @-> returning gen)
  let mulru = foreign "mulru" (gen @-> pari_ulong @-> returning gen)
  let mulsiz = foreign "mulsiz" (long @-> gen @-> gen @-> returning void)
  let mulsrz = foreign "mulsrz" (long @-> gen @-> gen @-> returning void)
  let mulssz = foreign "mulssz" (long @-> long @-> gen @-> returning void)
  let negi = foreign "negi" (gen @-> returning gen)
  let negr = foreign "negr" (gen @-> returning gen)
  let new_chunk = foreign "new_chunk" (int @-> returning gen)
  let rcopy = foreign "rcopy" (gen @-> returning gen)
  let rdivii = foreign "rdivii" (gen @-> gen @-> long @-> returning gen)
  let rdiviiz = foreign "rdiviiz" (gen @-> gen @-> gen @-> returning void)
  let rdivis = foreign "rdivis" (gen @-> long @-> long @-> returning gen)
  let rdivsi = foreign "rdivsi" (long @-> gen @-> long @-> returning gen)
  let rdivss = foreign "rdivss" (long @-> long @-> long @-> returning gen)
  let real2n = foreign "real2n" (long @-> long @-> returning gen)
  let real_m2n = foreign "real_m2n" (long @-> long @-> returning gen)
  let real_0 = foreign "real_0" (long @-> returning gen)
  let real_0_bit = foreign "real_0_bit" (long @-> returning gen)
  let real_1 = foreign "real_1" (long @-> returning gen)
  let real_1_bit = foreign "real_1_bit" (long @-> returning gen)
  let real_m1 = foreign "real_m1" (long @-> returning gen)
  let remii = foreign "remii" (gen @-> gen @-> returning gen)
  let remiiz = foreign "remiiz" (gen @-> gen @-> gen @-> returning void)
  let remis = foreign "remis" (gen @-> long @-> returning gen)
  let remisz = foreign "remisz" (gen @-> long @-> gen @-> returning void)

  let remlll_pre =
    foreign "remlll_pre"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong @-> pari_ulong
     @-> returning pari_ulong)

  let remsi = foreign "remsi" (long @-> gen @-> returning gen)
  let remsiz = foreign "remsiz" (long @-> gen @-> gen @-> returning void)
  let remss = foreign "remss" (long @-> long @-> returning gen)
  let remssz = foreign "remssz" (long @-> long @-> gen @-> returning void)
  let rtor = foreign "rtor" (gen @-> long @-> returning gen)
  let sdivsi = foreign "sdivsi" (long @-> gen @-> returning long)

  let sdivsi_rem =
    foreign "sdivsi_rem" (long @-> gen @-> ptr long @-> returning long)

  let sdivss_rem =
    foreign "sdivss_rem" (long @-> long @-> ptr long @-> returning long)

  let get_avma = foreign "get_avma" (void @-> returning pari_ulong)
  let set_avma = foreign "set_avma" (pari_ulong @-> returning void)

  let uabsdiviu_rem =
    foreign "uabsdiviu_rem"
      (gen @-> pari_ulong @-> ptr pari_ulong @-> returning pari_ulong)

  let udivuu_rem =
    foreign "udivuu_rem"
      (pari_ulong @-> pari_ulong @-> ptr pari_ulong @-> returning pari_ulong)

  let umodi2n = foreign "umodi2n" (gen @-> long @-> returning pari_ulong)
  let setabssign = foreign "setabssign" (gen @-> returning void)

  let shift_left =
    foreign "shift_left"
      (gen @-> gen @-> long @-> long @-> pari_ulong @-> pari_ulong
     @-> returning void)

  let shift_right =
    foreign "shift_right"
      (gen @-> gen @-> long @-> long @-> pari_ulong @-> pari_ulong
     @-> returning void)

  let shiftl =
    foreign "shiftl" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let shiftlr =
    foreign "shiftlr" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let shiftr = foreign "shiftr" (gen @-> long @-> returning gen)
  let shiftr_inplace = foreign "shiftr_inplace" (gen @-> long @-> returning void)
  let smodis = foreign "smodis" (gen @-> long @-> returning long)
  let smodss = foreign "smodss" (long @-> long @-> returning long)
  let stackdummy = foreign "stackdummy" (pari_sp @-> pari_sp @-> returning void)
  let stack_malloc = foreign "stack_malloc" (int @-> returning string)

  let stack_malloc_align =
    foreign "stack_malloc_align" (int @-> long @-> returning string)

  let stack_calloc = foreign "stack_calloc" (int @-> returning string)

  let stack_calloc_align =
    foreign "stack_calloc_align" (int @-> long @-> returning string)

  let stoi = foreign "stoi" (long @-> returning gen)
  let stor = foreign "stor" (long @-> long @-> returning gen)
  let subii = foreign "subii" (gen @-> gen @-> returning gen)
  let subiiz = foreign "subiiz" (gen @-> gen @-> gen @-> returning void)
  let subir = foreign "subir" (gen @-> gen @-> returning gen)
  let subirz = foreign "subirz" (gen @-> gen @-> gen @-> returning void)
  let subis = foreign "subis" (gen @-> long @-> returning gen)
  let subisz = foreign "subisz" (gen @-> long @-> gen @-> returning void)
  let subri = foreign "subri" (gen @-> gen @-> returning gen)
  let subriz = foreign "subriz" (gen @-> gen @-> gen @-> returning void)
  let subrr = foreign "subrr" (gen @-> gen @-> returning gen)
  let subrrz = foreign "subrrz" (gen @-> gen @-> gen @-> returning void)
  let subrs = foreign "subrs" (gen @-> long @-> returning gen)
  let subrsz = foreign "subrsz" (gen @-> long @-> gen @-> returning void)
  let subsi = foreign "subsi" (long @-> gen @-> returning gen)
  let subsiz = foreign "subsiz" (long @-> gen @-> gen @-> returning void)
  let subsrz = foreign "subsrz" (long @-> gen @-> gen @-> returning void)
  let subss = foreign "subss" (long @-> long @-> returning gen)
  let subssz = foreign "subssz" (long @-> long @-> gen @-> returning void)
  let subuu = foreign "subuu" (pari_ulong @-> pari_ulong @-> returning gen)
  let togglesign = foreign "togglesign" (gen @-> returning void)
end

module F54 (F : Ctypes.FOREIGN) = struct
  open F

  let togglesign_safe = foreign "togglesign_safe" (ptr gen @-> returning void)
  let affectsign = foreign "affectsign" (gen @-> gen @-> returning void)

  let affectsign_safe =
    foreign "affectsign_safe" (gen @-> ptr gen @-> returning void)

  let truedivii = foreign "truedivii" (gen @-> gen @-> returning gen)
  let truedivis = foreign "truedivis" (gen @-> long @-> returning gen)
  let truedivsi = foreign "truedivsi" (long @-> gen @-> returning gen)

  let uabsdivui_rem =
    foreign "uabsdivui_rem"
      (pari_ulong @-> gen @-> ptr pari_ulong @-> returning pari_ulong)

  let umodsu = foreign "umodsu" (long @-> pari_ulong @-> returning pari_ulong)
  let umodui = foreign "umodui" (pari_ulong @-> gen @-> returning pari_ulong)
  let ugcdiu = foreign "ugcdiu" (gen @-> pari_ulong @-> returning pari_ulong)
  let ugcdui = foreign "ugcdui" (pari_ulong @-> gen @-> returning pari_ulong)

  let umuluu_le =
    foreign "umuluu_le"
      (pari_ulong @-> pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let umuluu_or_0 =
    foreign "umuluu_or_0" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let utoi = foreign "utoi" (pari_ulong @-> returning gen)
  let utoineg = foreign "utoineg" (pari_ulong @-> returning gen)
  let utoipos = foreign "utoipos" (pari_ulong @-> returning gen)
  let utor = foreign "utor" (pari_ulong @-> long @-> returning gen)
  let uutoi = foreign "uutoi" (pari_ulong @-> pari_ulong @-> returning gen)
  let uutoineg = foreign "uutoineg" (pari_ulong @-> pari_ulong @-> returning gen)
  let vali = foreign "vali" (gen @-> returning long)
  let varncmp = foreign "varncmp" (long @-> long @-> returning int)
  let varnmax = foreign "varnmax" (long @-> long @-> returning long)
  let varnmin = foreign "varnmin" (long @-> long @-> returning long)
  let zm_zv_mod = foreign "ZM_ZV_mod" (gen @-> gen @-> returning gen)
  let zv_zv_mod = foreign "ZV_ZV_mod" (gen @-> gen @-> returning gen)
  let abgrp_get_cyc = foreign "abgrp_get_cyc" (gen @-> returning gen)
  let abgrp_get_gen = foreign "abgrp_get_gen" (gen @-> returning gen)
  let abgrp_get_no = foreign "abgrp_get_no" (gen @-> returning gen)
  let bid_get_u = foreign "bid_get_U" (gen @-> returning gen)
  let bid_get_arch = foreign "bid_get_arch" (gen @-> returning gen)
  let bid_get_archp = foreign "bid_get_archp" (gen @-> returning gen)
  let bid_get_cyc = foreign "bid_get_cyc" (gen @-> returning gen)
  let bid_get_fact = foreign "bid_get_fact" (gen @-> returning gen)
  let bid_get_fact2 = foreign "bid_get_fact2" (gen @-> returning gen)
  let bid_get_gen = foreign "bid_get_gen" (gen @-> returning gen)
  let bid_get_gen_nocheck = foreign "bid_get_gen_nocheck" (gen @-> returning gen)
  let bid_get_grp = foreign "bid_get_grp" (gen @-> returning gen)
  let bid_get_ideal = foreign "bid_get_ideal" (gen @-> returning gen)
  let bid_get_mod = foreign "bid_get_mod" (gen @-> returning gen)
  let bid_get_no = foreign "bid_get_no" (gen @-> returning gen)
  let bid_get_sarch = foreign "bid_get_sarch" (gen @-> returning gen)
  let bid_get_sprk = foreign "bid_get_sprk" (gen @-> returning gen)
  let bnf_get_clgp = foreign "bnf_get_clgp" (gen @-> returning gen)
  let bnf_get_cyc = foreign "bnf_get_cyc" (gen @-> returning gen)
  let bnf_get_fu = foreign "bnf_get_fu" (gen @-> returning gen)
  let bnf_get_fu_nocheck = foreign "bnf_get_fu_nocheck" (gen @-> returning gen)
  let bnf_get_gen = foreign "bnf_get_gen" (gen @-> returning gen)
  let bnf_get_logfu = foreign "bnf_get_logfu" (gen @-> returning gen)
  let bnf_get_nf = foreign "bnf_get_nf" (gen @-> returning gen)
  let bnf_get_no = foreign "bnf_get_no" (gen @-> returning gen)
  let bnf_get_reg = foreign "bnf_get_reg" (gen @-> returning gen)
  let bnf_get_sunits = foreign "bnf_get_sunits" (gen @-> returning gen)
  let bnf_get_tuu = foreign "bnf_get_tuU" (gen @-> returning gen)
  let bnr_get_bid = foreign "bnr_get_bid" (gen @-> returning gen)
  let bnr_get_bnf = foreign "bnr_get_bnf" (gen @-> returning gen)
  let bnr_get_clgp = foreign "bnr_get_clgp" (gen @-> returning gen)
  let bnr_get_cyc = foreign "bnr_get_cyc" (gen @-> returning gen)
  let bnr_get_gen = foreign "bnr_get_gen" (gen @-> returning gen)
  let bnr_get_gen_nocheck = foreign "bnr_get_gen_nocheck" (gen @-> returning gen)
  let bnr_get_mod = foreign "bnr_get_mod" (gen @-> returning gen)
  let bnr_get_nf = foreign "bnr_get_nf" (gen @-> returning gen)
  let bnr_get_no = foreign "bnr_get_no" (gen @-> returning gen)
  let cyc_get_expo = foreign "cyc_get_expo" (gen @-> returning gen)
  let ellqp_get_p = foreign "ellQp_get_p" (gen @-> returning gen)
  let ellqp_get_zero = foreign "ellQp_get_zero" (gen @-> returning gen)
  let ell_get_a1 = foreign "ell_get_a1" (gen @-> returning gen)
  let ell_get_a2 = foreign "ell_get_a2" (gen @-> returning gen)
  let ell_get_a3 = foreign "ell_get_a3" (gen @-> returning gen)
  let ell_get_a4 = foreign "ell_get_a4" (gen @-> returning gen)
  let ell_get_a6 = foreign "ell_get_a6" (gen @-> returning gen)
  let ell_get_b2 = foreign "ell_get_b2" (gen @-> returning gen)
  let ell_get_b4 = foreign "ell_get_b4" (gen @-> returning gen)
  let ell_get_b6 = foreign "ell_get_b6" (gen @-> returning gen)
  let ell_get_b8 = foreign "ell_get_b8" (gen @-> returning gen)
  let ell_get_c4 = foreign "ell_get_c4" (gen @-> returning gen)
  let ell_get_c6 = foreign "ell_get_c6" (gen @-> returning gen)
  let ell_get_disc = foreign "ell_get_disc" (gen @-> returning gen)
  let ell_get_j = foreign "ell_get_j" (gen @-> returning gen)
  let ellff_get_a4a6 = foreign "ellff_get_a4a6" (gen @-> returning gen)
  let ellff_get_field = foreign "ellff_get_field" (gen @-> returning gen)
  let ellinf = foreign "ellinf" (void @-> returning gen)
  let ellnf_get_bnf = foreign "ellnf_get_bnf" (gen @-> returning gen)
  let ellnf_get_nf = foreign "ellnf_get_nf" (gen @-> returning gen)
  let gal_get_den = foreign "gal_get_den" (gen @-> returning gen)
  let gal_get_e = foreign "gal_get_e" (gen @-> returning gen)
  let gal_get_gen = foreign "gal_get_gen" (gen @-> returning gen)
  let gal_get_group = foreign "gal_get_group" (gen @-> returning gen)
  let gal_get_invvdm = foreign "gal_get_invvdm" (gen @-> returning gen)
  let gal_get_mod = foreign "gal_get_mod" (gen @-> returning gen)
  let gal_get_orders = foreign "gal_get_orders" (gen @-> returning gen)
  let gal_get_p = foreign "gal_get_p" (gen @-> returning gen)
  let gal_get_pol = foreign "gal_get_pol" (gen @-> returning gen)
  let gal_get_roots = foreign "gal_get_roots" (gen @-> returning gen)

  let idealchineseinit =
    foreign "idealchineseinit" (gen @-> gen @-> returning gen)

  let idealred = foreign "idealred" (gen @-> gen @-> returning gen)
  let modpr_get_t = foreign "modpr_get_T" (gen @-> returning gen)
  let modpr_get_p = foreign "modpr_get_p" (gen @-> returning gen)
  let modpr_get_pr = foreign "modpr_get_pr" (gen @-> returning gen)

  let nfv_to_scalar_or_alg =
    foreign "nfV_to_scalar_or_alg" (gen @-> gen @-> returning gen)

  let nf_get_g = foreign "nf_get_G" (gen @-> returning gen)
end

module F55 (F : Ctypes.FOREIGN) = struct
  open F

  let nf_get_m = foreign "nf_get_M" (gen @-> returning gen)
  let nf_get_tr = foreign "nf_get_Tr" (gen @-> returning gen)
  let nf_get_diff = foreign "nf_get_diff" (gen @-> returning gen)
  let nf_get_disc = foreign "nf_get_disc" (gen @-> returning gen)
  let nf_get_index = foreign "nf_get_index" (gen @-> returning gen)
  let nf_get_invzk = foreign "nf_get_invzk" (gen @-> returning gen)
  let nf_get_pol = foreign "nf_get_pol" (gen @-> returning gen)

  let nf_get_ramified_primes =
    foreign "nf_get_ramified_primes" (gen @-> returning gen)

  let nf_get_roots = foreign "nf_get_roots" (gen @-> returning gen)
  let nf_get_roundg = foreign "nf_get_roundG" (gen @-> returning gen)
  let nf_get_zk = foreign "nf_get_zk" (gen @-> returning gen)
  let nf_get_zkden = foreign "nf_get_zkden" (gen @-> returning gen)
  let nf_get_zkprimpart = foreign "nf_get_zkprimpart" (gen @-> returning gen)
  let pr_get_gen = foreign "pr_get_gen" (gen @-> returning gen)
  let pr_get_p = foreign "pr_get_p" (gen @-> returning gen)
  let pr_get_tau = foreign "pr_get_tau" (gen @-> returning gen)
  let pr_norm = foreign "pr_norm" (gen @-> returning gen)
  let rnf_get_alpha = foreign "rnf_get_alpha" (gen @-> returning gen)
  let rnf_get_disc = foreign "rnf_get_disc" (gen @-> returning gen)
  let rnf_get_idealdisc = foreign "rnf_get_idealdisc" (gen @-> returning gen)
  let rnf_get_index = foreign "rnf_get_index" (gen @-> returning gen)
  let rnf_get_invzk = foreign "rnf_get_invzk" (gen @-> returning gen)
  let rnf_get_k = foreign "rnf_get_k" (gen @-> returning gen)
  let rnf_get_map = foreign "rnf_get_map" (gen @-> returning gen)
  let rnf_get_nf = foreign "rnf_get_nf" (gen @-> returning gen)
  let rnf_get_nfpol = foreign "rnf_get_nfpol" (gen @-> returning gen)
  let rnf_get_nfzk = foreign "rnf_get_nfzk" (gen @-> returning gen)
  let rnf_get_pol = foreign "rnf_get_pol" (gen @-> returning gen)
  let rnf_get_polabs = foreign "rnf_get_polabs" (gen @-> returning gen)

  let rnf_get_ramified_primes =
    foreign "rnf_get_ramified_primes" (gen @-> returning gen)

  let rnf_get_zk = foreign "rnf_get_zk" (gen @-> returning gen)
  let vecmodii = foreign "vecmodii" (gen @-> gen @-> returning gen)
  let vecmoduu = foreign "vecmoduu" (gen @-> gen @-> returning gen)
  let znstar_get_n = foreign "znstar_get_N" (gen @-> returning gen)
  let znstar_get_u = foreign "znstar_get_U" (gen @-> returning gen)
  let znstar_get_ui = foreign "znstar_get_Ui" (gen @-> returning gen)

  let znstar_get_conreycyc =
    foreign "znstar_get_conreycyc" (gen @-> returning gen)

  let znstar_get_conreygen =
    foreign "znstar_get_conreygen" (gen @-> returning gen)

  let znstar_get_cyc = foreign "znstar_get_cyc" (gen @-> returning gen)
  let znstar_get_fan = foreign "znstar_get_faN" (gen @-> returning gen)
  let znstar_get_gen = foreign "znstar_get_gen" (gen @-> returning gen)
  let znstar_get_no = foreign "znstar_get_no" (gen @-> returning gen)
  let znstar_get_pe = foreign "znstar_get_pe" (gen @-> returning gen)
  let checkell_i = foreign "checkell_i" (gen @-> returning int)
  let ell_is_inf = foreign "ell_is_inf" (gen @-> returning int)
  let pr_is_inert = foreign "pr_is_inert" (gen @-> returning int)
  let bnf_get_tun = foreign "bnf_get_tuN" (gen @-> returning long)
  let ellqp_get_prec = foreign "ellQp_get_prec" (gen @-> returning long)
  let ellr_get_prec = foreign "ellR_get_prec" (gen @-> returning long)
  let ellr_get_sign = foreign "ellR_get_sign" (gen @-> returning long)
  let ell_get_type = foreign "ell_get_type" (gen @-> returning long)
  let logint = foreign "logint" (gen @-> gen @-> returning long)
  let nf_get_degree = foreign "nf_get_degree" (gen @-> returning long)
  let nf_get_r1 = foreign "nf_get_r1" (gen @-> returning long)
  let nf_get_r2 = foreign "nf_get_r2" (gen @-> returning long)
  let nf_get_varn = foreign "nf_get_varn" (gen @-> returning long)
  let pr_get_e = foreign "pr_get_e" (gen @-> returning long)
  let pr_get_f = foreign "pr_get_f" (gen @-> returning long)
  let rnf_get_absdegree = foreign "rnf_get_absdegree" (gen @-> returning long)
  let rnf_get_degree = foreign "rnf_get_degree" (gen @-> returning long)
  let rnf_get_nfdegree = foreign "rnf_get_nfdegree" (gen @-> returning long)
  let rnf_get_nfvarn = foreign "rnf_get_nfvarn" (gen @-> returning long)
  let rnf_get_varn = foreign "rnf_get_varn" (gen @-> returning long)
  let hash_str = foreign "hash_str" (string @-> returning pari_ulong)

  let hash_str_len =
    foreign "hash_str_len" (string @-> long @-> returning pari_ulong)

  let ulogint =
    foreign "ulogint" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let upr_norm = foreign "upr_norm" (gen @-> returning pari_ulong)

  let nf_get_sign =
    foreign "nf_get_sign" (gen @-> ptr long @-> ptr long @-> returning void)

  let closure_arity = foreign "closure_arity" (gen @-> returning long)
  let closure_codestr = foreign "closure_codestr" (gen @-> returning string)
  let closure_get_code = foreign "closure_get_code" (gen @-> returning gen)
  let closure_get_oper = foreign "closure_get_oper" (gen @-> returning gen)
  let closure_get_data = foreign "closure_get_data" (gen @-> returning gen)
  let closure_get_dbg = foreign "closure_get_dbg" (gen @-> returning gen)
  let closure_get_text = foreign "closure_get_text" (gen @-> returning gen)
  let closure_get_frame = foreign "closure_get_frame" (gen @-> returning gen)

  let closure_is_variadic =
    foreign "closure_is_variadic" (gen @-> returning long)

  let addmuliu =
    foreign "addmuliu" (gen @-> gen @-> pari_ulong @-> returning gen)

  let addmuliu_inplace =
    foreign "addmuliu_inplace" (gen @-> gen @-> pari_ulong @-> returning gen)

  let lincombii =
    foreign "lincombii" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let mulsubii = foreign "mulsubii" (gen @-> gen @-> gen @-> returning gen)
  let submulii = foreign "submulii" (gen @-> gen @-> gen @-> returning gen)

  let submuliu =
    foreign "submuliu" (gen @-> gen @-> pari_ulong @-> returning gen)

  let submuliu_inplace =
    foreign "submuliu_inplace" (gen @-> gen @-> pari_ulong @-> returning gen)

  let fpxq_add =
    foreign "FpXQ_add" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_sub =
    foreign "FpXQ_sub" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let flxq_add =
    foreign "Flxq_add" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_sub =
    foreign "Flxq_sub" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let fpxqx_div =
    foreign "FpXQX_div" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let flxqx_div =
    foreign "FlxqX_div" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_div_pre =
    foreign "FlxqX_div_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let f2xqx_div = foreign "F2xqX_div" (gen @-> gen @-> gen @-> returning gen)
  let rg_to_fq = foreign "Rg_to_Fq" (gen @-> gen @-> gen @-> returning gen)
  let fq_red = foreign "Fq_red" (gen @-> gen @-> gen @-> returning gen)
  let fq_to_fpxq = foreign "Fq_to_FpXQ" (gen @-> gen @-> gen @-> returning gen)

  let gener_fq_local =
    foreign "gener_Fq_local" (gen @-> gen @-> gen @-> returning gen)

  let fpxy_fq_evaly =
    foreign "FpXY_Fq_evaly"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let fqx_fp_mul =
    foreign "FqX_Fp_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_fq_mul =
    foreign "FqX_Fq_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_add = foreign "FqX_add" (gen @-> gen @-> gen @-> gen @-> returning gen)
end

module F56 (F : Ctypes.FOREIGN) = struct
  open F

  let fqx_ddf = foreign "FqX_ddf" (gen @-> gen @-> gen @-> returning gen)
  let fqx_degfact = foreign "FqX_degfact" (gen @-> gen @-> gen @-> returning gen)
  let fqx_deriv = foreign "FqX_deriv" (gen @-> gen @-> gen @-> returning gen)
  let fqx_div = foreign "FqX_div" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_div_by_x_x =
    foreign "FqX_div_by_X_x"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fqx_divrem =
    foreign "FqX_divrem"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fqx_extgcd =
    foreign "FqX_extgcd"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let fqx_factor = foreign "FqX_factor" (gen @-> gen @-> gen @-> returning gen)

  let fqx_factor_squarefree =
    foreign "FqX_factor_squarefree" (gen @-> gen @-> gen @-> returning gen)

  let fqx_gcd = foreign "FqX_gcd" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fqx_get_red = foreign "FqX_get_red" (gen @-> gen @-> gen @-> returning gen)

  let fqx_halfgcd =
    foreign "FqX_halfgcd" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_halve = foreign "FqX_halve" (gen @-> gen @-> gen @-> returning gen)
  let fqx_integ = foreign "FqX_integ" (gen @-> gen @-> gen @-> returning gen)
  let fqx_mul = foreign "FqX_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_mulu =
    foreign "FqX_mulu" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fqx_neg = foreign "FqX_neg" (gen @-> gen @-> gen @-> returning gen)

  let fqx_normalize =
    foreign "FqX_normalize" (gen @-> gen @-> gen @-> returning gen)

  let fqx_powu =
    foreign "FqX_powu" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fqx_red = foreign "FqX_red" (gen @-> gen @-> gen @-> returning gen)
  let fqx_rem = foreign "FqX_rem" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fqx_roots = foreign "FqX_roots" (gen @-> gen @-> gen @-> returning gen)
  let fqx_sqr = foreign "FqX_sqr" (gen @-> gen @-> gen @-> returning gen)
  let fqx_sub = foreign "FqX_sub" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fqx_to_mod = foreign "FqX_to_mod" (gen @-> gen @-> gen @-> returning gen)

  let fqxq_add =
    foreign "FqXQ_add" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_div =
    foreign "FqXQ_div" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_inv =
    foreign "FqXQ_inv" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_invsafe =
    foreign "FqXQ_invsafe" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_mul =
    foreign "FqXQ_mul" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_pow =
    foreign "FqXQ_pow" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_sqr =
    foreign "FqXQ_sqr" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_sub =
    foreign "FqXQ_sub" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxn_exp =
    foreign "FqXn_exp" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fqxn_expint =
    foreign "FqXn_expint" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fqxn_inv =
    foreign "FqXn_inv" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fqxn_mul =
    foreign "FqXn_mul" (gen @-> gen @-> long @-> gen @-> gen @-> returning gen)

  let fqxn_sqr =
    foreign "FqXn_sqr" (gen @-> long @-> gen @-> gen @-> returning gen)

  let get_f2x_degree = foreign "get_F2x_degree" (gen @-> returning long)
  let get_f2x_mod = foreign "get_F2x_mod" (gen @-> returning gen)
  let get_f2x_var = foreign "get_F2x_var" (gen @-> returning long)
  let get_f2xqx_degree = foreign "get_F2xqX_degree" (gen @-> returning long)
  let get_f2xqx_mod = foreign "get_F2xqX_mod" (gen @-> returning gen)
  let get_f2xqx_var = foreign "get_F2xqX_var" (gen @-> returning long)
  let get_flx_degree = foreign "get_Flx_degree" (gen @-> returning long)
  let get_flx_mod = foreign "get_Flx_mod" (gen @-> returning gen)
  let get_flx_var = foreign "get_Flx_var" (gen @-> returning long)
  let get_flxqx_degree = foreign "get_FlxqX_degree" (gen @-> returning long)
  let get_flxqx_mod = foreign "get_FlxqX_mod" (gen @-> returning gen)
  let get_flxqx_var = foreign "get_FlxqX_var" (gen @-> returning long)
  let get_fpx_degree = foreign "get_FpX_degree" (gen @-> returning long)
  let get_fpx_mod = foreign "get_FpX_mod" (gen @-> returning gen)
  let get_fpx_var = foreign "get_FpX_var" (gen @-> returning long)
  let get_fpxqx_degree = foreign "get_FpXQX_degree" (gen @-> returning long)
  let get_fpxqx_mod = foreign "get_FpXQX_mod" (gen @-> returning gen)
  let get_fpxqx_var = foreign "get_FpXQX_var" (gen @-> returning long)

  let f2m_coeff =
    foreign "F2m_coeff" (gen @-> long @-> long @-> returning pari_ulong)

  let f2m_clear = foreign "F2m_clear" (gen @-> long @-> long @-> returning void)
  let f2m_flip = foreign "F2m_flip" (gen @-> long @-> long @-> returning void)
  let f2m_set = foreign "F2m_set" (gen @-> long @-> long @-> returning void)
  let f2v_clear = foreign "F2v_clear" (gen @-> long @-> returning void)
  let f2v_coeff = foreign "F2v_coeff" (gen @-> long @-> returning pari_ulong)
  let f2v_flip = foreign "F2v_flip" (gen @-> long @-> returning void)
  let f2v_to_f2x = foreign "F2v_to_F2x" (gen @-> long @-> returning gen)
  let f2v_set = foreign "F2v_set" (gen @-> long @-> returning void)
  let f2x_clear = foreign "F2x_clear" (gen @-> long @-> returning void)
  let f2x_coeff = foreign "F2x_coeff" (gen @-> long @-> returning pari_ulong)
  let f2x_flip = foreign "F2x_flip" (gen @-> long @-> returning void)
  let f2x_set = foreign "F2x_set" (gen @-> long @-> returning void)
  let f2x_equal1 = foreign "F2x_equal1" (gen @-> returning int)
  let f2x_equal = foreign "F2x_equal" (gen @-> gen @-> returning int)
  let f2x_div = foreign "F2x_div" (gen @-> gen @-> returning gen)

  let f2x_renormalize =
    foreign "F2x_renormalize" (gen @-> long @-> returning gen)

  let f2m_copy = foreign "F2m_copy" (gen @-> returning gen)
  let f2v_copy = foreign "F2v_copy" (gen @-> returning gen)
  let f2x_copy = foreign "F2x_copy" (gen @-> returning gen)
  let f2v_ei = foreign "F2v_ei" (long @-> long @-> returning gen)

  let f3m_coeff =
    foreign "F3m_coeff" (gen @-> long @-> long @-> returning pari_ulong)

  let f3m_set =
    foreign "F3m_set" (gen @-> long @-> long @-> pari_ulong @-> returning void)

  let f3m_copy = foreign "F3m_copy" (gen @-> returning gen)
  let flm_copy = foreign "Flm_copy" (gen @-> returning gen)
  let flv_copy = foreign "Flv_copy" (gen @-> returning gen)
  let flx_equal1 = foreign "Flx_equal1" (gen @-> returning int)
  let flx_constant = foreign "Flx_constant" (gen @-> returning pari_ulong)
  let flx_copy = foreign "Flx_copy" (gen @-> returning gen)
  let flx_div = foreign "Flx_div" (gen @-> gen @-> pari_ulong @-> returning gen)

  let flx_div_pre =
    foreign "Flx_div_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let flx_lead = foreign "Flx_lead" (gen @-> returning pari_ulong)

  let flx_mulu =
    foreign "Flx_mulu" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let fp_divu = foreign "Fp_divu" (gen @-> pari_ulong @-> gen @-> returning gen)
  let fpv_fpc_mul = foreign "FpV_FpC_mul" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_renormalize =
    foreign "FpXQX_renormalize" (gen @-> long @-> returning gen)

  let fpxx_renormalize =
    foreign "FpXX_renormalize" (gen @-> long @-> returning gen)

  let fpx_div = foreign "FpX_div" (gen @-> gen @-> gen @-> returning gen)

  let fpx_renormalize =
    foreign "FpX_renormalize" (gen @-> long @-> returning gen)

  let fp_add = foreign "Fp_add" (gen @-> gen @-> gen @-> returning gen)

  let fp_addmul =
    foreign "Fp_addmul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fp_center = foreign "Fp_center" (gen @-> gen @-> gen @-> returning gen)
  let fp_center_i = foreign "Fp_center_i" (gen @-> gen @-> gen @-> returning gen)
  let fp_div = foreign "Fp_div" (gen @-> gen @-> gen @-> returning gen)
end

module F57 (F : Ctypes.FOREIGN) = struct
  open F

  let fp_halve = foreign "Fp_halve" (gen @-> gen @-> returning gen)
  let fp_inv = foreign "Fp_inv" (gen @-> gen @-> returning gen)
  let fp_invsafe = foreign "Fp_invsafe" (gen @-> gen @-> returning gen)
  let fp_mul = foreign "Fp_mul" (gen @-> gen @-> gen @-> returning gen)
  let fp_muls = foreign "Fp_muls" (gen @-> long @-> gen @-> returning gen)
  let fp_mulu = foreign "Fp_mulu" (gen @-> pari_ulong @-> gen @-> returning gen)
  let fp_neg = foreign "Fp_neg" (gen @-> gen @-> returning gen)
  let fp_red = foreign "Fp_red" (gen @-> gen @-> returning gen)
  let fp_sqr = foreign "Fp_sqr" (gen @-> gen @-> returning gen)
  let fp_sub = foreign "Fp_sub" (gen @-> gen @-> gen @-> returning gen)
  let genbinbase = foreign "GENbinbase" (ptr genbin @-> returning gen)
  let q_abs = foreign "Q_abs" (gen @-> returning gen)
  let q_abs_shallow = foreign "Q_abs_shallow" (gen @-> returning gen)
  let qv_isscalar = foreign "QV_isscalar" (gen @-> returning int)
  let qtoss = foreign "Qtoss" (gen @-> ptr long @-> ptr long @-> returning void)
  let r_abs = foreign "R_abs" (gen @-> returning gen)
  let r_abs_shallow = foreign "R_abs_shallow" (gen @-> returning gen)
  let rgc_fpnorml2 = foreign "RgC_fpnorml2" (gen @-> long @-> returning gen)
  let rgc_gtofp = foreign "RgC_gtofp" (gen @-> long @-> returning gen)
  let rgc_gtomp = foreign "RgC_gtomp" (gen @-> long @-> returning gen)

  let rgm_dimensions =
    foreign "RgM_dimensions" (gen @-> ptr long @-> ptr long @-> returning void)

  let rgm_fpnorml2 = foreign "RgM_fpnorml2" (gen @-> long @-> returning gen)
  let rgm_gtofp = foreign "RgM_gtofp" (gen @-> long @-> returning gen)
  let rgm_gtomp = foreign "RgM_gtomp" (gen @-> long @-> returning gen)
  let rgm_minor = foreign "RgM_minor" (gen @-> long @-> long @-> returning gen)
  let rgm_shallowcopy = foreign "RgM_shallowcopy" (gen @-> returning gen)
  let rgv_gtofp = foreign "RgV_gtofp" (gen @-> long @-> returning gen)
  let rgv_isscalar = foreign "RgV_isscalar" (gen @-> returning int)
  let rgv_isin = foreign "RgV_isin" (gen @-> gen @-> returning long)
  let rgv_isin_i = foreign "RgV_isin_i" (gen @-> gen @-> long @-> returning long)
  let rgv_is_zv = foreign "RgV_is_ZV" (gen @-> returning int)
  let rgv_is_qv = foreign "RgV_is_QV" (gen @-> returning int)
  let rgx_equal_var = foreign "RgX_equal_var" (gen @-> gen @-> returning long)
  let rgx_is_monomial = foreign "RgX_is_monomial" (gen @-> returning int)
  let rgx_is_rational = foreign "RgX_is_rational" (gen @-> returning int)
  let rgx_is_qx = foreign "RgX_is_QX" (gen @-> returning int)
  let rgx_is_zx = foreign "RgX_is_ZX" (gen @-> returning int)
  let rgx_isscalar = foreign "RgX_isscalar" (gen @-> returning int)

  let rgx_shift_inplace =
    foreign "RgX_shift_inplace" (gen @-> long @-> returning gen)

  let rgx_shift_inplace_init =
    foreign "RgX_shift_inplace_init" (long @-> returning void)

  let rgx_to_rgv = foreign "RgX_to_RgV" (gen @-> long @-> returning gen)
  let rgxqx_div = foreign "RgXQX_div" (gen @-> gen @-> gen @-> returning gen)
  let rgxqx_rem = foreign "RgXQX_rem" (gen @-> gen @-> gen @-> returning gen)
  let rgx_coeff = foreign "RgX_coeff" (gen @-> long @-> returning gen)
  let rgx_copy = foreign "RgX_copy" (gen @-> returning gen)
  let rgx_div = foreign "RgX_div" (gen @-> gen @-> returning gen)
  let rgx_fpnorml2 = foreign "RgX_fpnorml2" (gen @-> long @-> returning gen)
  let rgx_gtofp = foreign "RgX_gtofp" (gen @-> long @-> returning gen)
  let rgx_renormalize = foreign "RgX_renormalize" (gen @-> returning gen)
  let rg_col_ei = foreign "Rg_col_ei" (gen @-> long @-> long @-> returning gen)
  let zc_hnfrem = foreign "ZC_hnfrem" (gen @-> gen @-> returning gen)
  let zm_hnfrem = foreign "ZM_hnfrem" (gen @-> gen @-> returning gen)
  let zm_lll = foreign "ZM_lll" (gen @-> double @-> long @-> returning gen)
  let zv_dvd = foreign "ZV_dvd" (gen @-> gen @-> returning int)
  let zv_isscalar = foreign "ZV_isscalar" (gen @-> returning int)
  let zv_to_zv = foreign "ZV_to_zv" (gen @-> returning gen)
  let zx_equal1 = foreign "ZX_equal1" (gen @-> returning int)
  let zx_is_monic = foreign "ZX_is_monic" (gen @-> returning int)
  let zx_renormalize = foreign "ZX_renormalize" (gen @-> long @-> returning gen)
  let zxq_mul = foreign "ZXQ_mul" (gen @-> gen @-> gen @-> returning gen)
  let zxq_sqr = foreign "ZXQ_sqr" (gen @-> gen @-> returning gen)
  let z_ispower = foreign "Z_ispower" (gen @-> pari_ulong @-> returning long)
  let z_issquare = foreign "Z_issquare" (gen @-> returning long)
  let absfrac = foreign "absfrac" (gen @-> returning gen)
  let absfrac_shallow = foreign "absfrac_shallow" (gen @-> returning gen)
  let affc_fixlg = foreign "affc_fixlg" (gen @-> gen @-> returning gen)
  let bin_copy = foreign "bin_copy" (ptr genbin @-> returning gen)
  let bit_accuracy = foreign "bit_accuracy" (long @-> returning long)

  let bit_accuracy_mul =
    foreign "bit_accuracy_mul" (long @-> double @-> returning double)

  let bit_prec = foreign "bit_prec" (gen @-> returning long)
  let both_odd = foreign "both_odd" (long @-> long @-> returning int)
  let cbrtr = foreign "cbrtr" (gen @-> returning gen)
  let cbrtr_abs = foreign "cbrtr_abs" (gen @-> returning gen)
  let cgetc = foreign "cgetc" (long @-> returning gen)
  let cgetalloc = foreign "cgetalloc" (int @-> long @-> returning gen)
  let cgiv = foreign "cgiv" (gen @-> returning void)
  let col_ei = foreign "col_ei" (long @-> long @-> returning gen)
  let const_col = foreign "const_col" (long @-> gen @-> returning gen)
  let const_vec = foreign "const_vec" (long @-> gen @-> returning gen)
  let const_vecsmall = foreign "const_vecsmall" (long @-> long @-> returning gen)
  let constant_coeff = foreign "constant_coeff" (gen @-> returning gen)
  let cxcompotor = foreign "cxcompotor" (gen @-> long @-> returning gen)
  let cxnorm = foreign "cxnorm" (gen @-> returning gen)
  let cxtoreal = foreign "cxtoreal" (gen @-> returning gen)
  let cyclic_perm = foreign "cyclic_perm" (long @-> long @-> returning gen)
  let dbllog2r = foreign "dbllog2r" (gen @-> returning double)
  let degpol = foreign "degpol" (gen @-> returning long)
  let divsbil = foreign "divsBIL" (long @-> returning long)
  let gabsz = foreign "gabsz" (gen @-> long @-> gen @-> returning void)
  let gaddgs = foreign "gaddgs" (gen @-> long @-> returning gen)
  let gaddz = foreign "gaddz" (gen @-> gen @-> gen @-> returning void)
  let gc_all = foreign "gc_all" (pari_sp @-> int @-> returning gen)
  let gcmpgs = foreign "gcmpgs" (gen @-> long @-> returning int)
  let gdiventz = foreign "gdiventz" (gen @-> gen @-> gen @-> returning void)
  let gdivsg = foreign "gdivsg" (long @-> gen @-> returning gen)
  let gdivz = foreign "gdivz" (gen @-> gen @-> gen @-> returning void)
  let gen_i = foreign "gen_I" (void @-> returning gen)
  let gerepileall = foreign "gerepileall" (pari_sp @-> int @-> returning void)

  let gerepilecoeffs =
    foreign "gerepilecoeffs" (pari_sp @-> gen @-> int @-> returning void)

  let gerepilecopy = foreign "gerepilecopy" (pari_sp @-> gen @-> returning gen)
end

module F58 (F : Ctypes.FOREIGN) = struct
  open F

  let gerepilemany =
    foreign "gerepilemany" (pari_sp @-> ptr (ptr gen) @-> int @-> returning void)

  let gequalgs = foreign "gequalgs" (gen @-> long @-> returning int)
  let gerepileupto = foreign "gerepileupto" (pari_sp @-> gen @-> returning gen)

  let gerepileuptoint =
    foreign "gerepileuptoint" (pari_sp @-> gen @-> returning gen)

  let gerepileuptoleaf =
    foreign "gerepileuptoleaf" (pari_sp @-> gen @-> returning gen)

  let gisdouble = foreign "gisdouble" (gen @-> ptr double @-> returning int)
  let gmax_shallow = foreign "gmax_shallow" (gen @-> gen @-> returning gen)
  let gmaxsg = foreign "gmaxsg" (long @-> gen @-> returning gen)
  let gmin_shallow = foreign "gmin_shallow" (gen @-> gen @-> returning gen)
  let gminsg = foreign "gminsg" (long @-> gen @-> returning gen)
  let gmodz = foreign "gmodz" (gen @-> gen @-> gen @-> returning void)
  let gmul2nz = foreign "gmul2nz" (gen @-> long @-> gen @-> returning void)
  let gmulgs = foreign "gmulgs" (gen @-> long @-> returning gen)
  let gmulgu = foreign "gmulgu" (gen @-> pari_ulong @-> returning gen)
  let gmulz = foreign "gmulz" (gen @-> gen @-> gen @-> returning void)
  let gnegz = foreign "gnegz" (gen @-> gen @-> returning void)
  let gshiftz = foreign "gshiftz" (gen @-> long @-> gen @-> returning void)
  let gsubgs = foreign "gsubgs" (gen @-> long @-> returning gen)
  let gsubz = foreign "gsubz" (gen @-> gen @-> gen @-> returning void)
  let gtodouble = foreign "gtodouble" (gen @-> returning double)
  let gtofp = foreign "gtofp" (gen @-> long @-> returning gen)
  let gtomp = foreign "gtomp" (gen @-> long @-> returning gen)
  let gtos = foreign "gtos" (gen @-> returning long)
  let gtou = foreign "gtou" (gen @-> returning pari_ulong)
  let gunclonenull = foreign "guncloneNULL" (gen @-> returning void)
  let gunclonenull_deep = foreign "guncloneNULL_deep" (gen @-> returning void)
  let gval = foreign "gval" (gen @-> long @-> returning long)
  let identity_perm = foreign "identity_perm" (long @-> returning gen)
  let identity_zv = foreign "identity_zv" (long @-> returning gen)
  let identity_zv = foreign "identity_ZV" (long @-> returning gen)
  let equali1 = foreign "equali1" (gen @-> returning int)
  let equalim1 = foreign "equalim1" (gen @-> returning int)
  let inf_get_sign = foreign "inf_get_sign" (gen @-> returning long)
  let inv_content = foreign "inv_content" (gen @-> returning gen)
  let is_bigint = foreign "is_bigint" (gen @-> returning int)
  let is_const_t = foreign "is_const_t" (long @-> returning int)
  let is_extscalar_t = foreign "is_extscalar_t" (long @-> returning int)
  let is_intreal_t = foreign "is_intreal_t" (long @-> returning int)
  let is_matvec_t = foreign "is_matvec_t" (long @-> returning int)
  let is_noncalc_t = foreign "is_noncalc_t" (long @-> returning int)
  let is_pm1 = foreign "is_pm1" (gen @-> returning int)
  let is_qfb_t = foreign "is_qfb_t" (long @-> returning int)
  let is_rational_t = foreign "is_rational_t" (long @-> returning int)
  let is_real_t = foreign "is_real_t" (long @-> returning int)
  let is_recursive_t = foreign "is_recursive_t" (long @-> returning int)
  let is_scalar_t = foreign "is_scalar_t" (long @-> returning int)

  let is_universal_constant =
    foreign "is_universal_constant" (gen @-> returning int)

  let is_vec_t = foreign "is_vec_t" (long @-> returning int)
  let isint1 = foreign "isint1" (gen @-> returning int)
  let isintm1 = foreign "isintm1" (gen @-> returning int)
  let isintzero = foreign "isintzero" (gen @-> returning int)
  let ismpzero = foreign "ismpzero" (gen @-> returning int)
  let isonstack = foreign "isonstack" (gen @-> returning int)
  let killblock = foreign "killblock" (gen @-> returning void)
  let leading_coeff = foreign "leading_coeff" (gen @-> returning gen)
  let lg_increase = foreign "lg_increase" (gen @-> returning void)
  let lgcols = foreign "lgcols" (gen @-> returning long)
  let lgpol = foreign "lgpol" (gen @-> returning long)
  let div_content = foreign "div_content" (gen @-> gen @-> returning gen)
  let matpascal = foreign "matpascal" (long @-> returning gen)

  let matslice =
    foreign "matslice"
      (gen @-> long @-> long @-> long @-> long @-> returning gen)

  let mkcol = foreign "mkcol" (gen @-> returning gen)
  let mkcol2 = foreign "mkcol2" (gen @-> gen @-> returning gen)
  let mkcol2s = foreign "mkcol2s" (long @-> long @-> returning gen)
  let mkcol3 = foreign "mkcol3" (gen @-> gen @-> gen @-> returning gen)
  let mkcol3s = foreign "mkcol3s" (long @-> long @-> long @-> returning gen)
  let mkcol4 = foreign "mkcol4" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkcol4s =
    foreign "mkcol4s" (long @-> long @-> long @-> long @-> returning gen)

  let mkcol5 =
    foreign "mkcol5" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkcol6 =
    foreign "mkcol6"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkcolcopy = foreign "mkcolcopy" (gen @-> returning gen)
  let mkcols = foreign "mkcols" (long @-> returning gen)
  let mkcomplex = foreign "mkcomplex" (gen @-> gen @-> returning gen)
  let mkerr = foreign "mkerr" (long @-> returning gen)
  let mkmoo = foreign "mkmoo" (void @-> returning gen)
  let mkoo = foreign "mkoo" (void @-> returning gen)
  let mkfrac = foreign "mkfrac" (gen @-> gen @-> returning gen)
  let mkfracss = foreign "mkfracss" (long @-> long @-> returning gen)
  let mkfraccopy = foreign "mkfraccopy" (gen @-> gen @-> returning gen)
  let mkintmod = foreign "mkintmod" (gen @-> gen @-> returning gen)

  let mkintmodu =
    foreign "mkintmodu" (pari_ulong @-> pari_ulong @-> returning gen)

  let mkmat = foreign "mkmat" (gen @-> returning gen)
  let mkmat2 = foreign "mkmat2" (gen @-> gen @-> returning gen)
  let mkmat22 = foreign "mkmat22" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkmat22s =
    foreign "mkmat22s" (long @-> long @-> long @-> long @-> returning gen)

  let mkmat3 = foreign "mkmat3" (gen @-> gen @-> gen @-> returning gen)
  let mkmat4 = foreign "mkmat4" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkmat5 =
    foreign "mkmat5" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkmatcopy = foreign "mkmatcopy" (gen @-> returning gen)
  let mkpolmod = foreign "mkpolmod" (gen @-> gen @-> returning gen)
  let mkqfb = foreign "mkqfb" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let mkquad = foreign "mkquad" (gen @-> gen @-> gen @-> returning gen)
  let mkrfrac = foreign "mkrfrac" (gen @-> gen @-> returning gen)
  let mkrfraccopy = foreign "mkrfraccopy" (gen @-> gen @-> returning gen)
  let mkvec = foreign "mkvec" (gen @-> returning gen)
  let mkvec2 = foreign "mkvec2" (gen @-> gen @-> returning gen)
  let mkvec2copy = foreign "mkvec2copy" (gen @-> gen @-> returning gen)
  let mkvec2s = foreign "mkvec2s" (long @-> long @-> returning gen)
  let mkvec3 = foreign "mkvec3" (gen @-> gen @-> gen @-> returning gen)
  let mkvec3s = foreign "mkvec3s" (long @-> long @-> long @-> returning gen)
end

module F59 (F : Ctypes.FOREIGN) = struct
  open F

  let mkvec4 = foreign "mkvec4" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkvec4s =
    foreign "mkvec4s" (long @-> long @-> long @-> long @-> returning gen)

  let mkvec5 =
    foreign "mkvec5" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkveccopy = foreign "mkveccopy" (gen @-> returning gen)
  let mkvecs = foreign "mkvecs" (long @-> returning gen)
  let mkvecsmall = foreign "mkvecsmall" (long @-> returning gen)
  let mkvecsmall2 = foreign "mkvecsmall2" (long @-> long @-> returning gen)

  let mkvecsmall3 =
    foreign "mkvecsmall3" (long @-> long @-> long @-> returning gen)

  let mkvecsmall4 =
    foreign "mkvecsmall4" (long @-> long @-> long @-> long @-> returning gen)

  let mkvecsmall5 =
    foreign "mkvecsmall5"
      (long @-> long @-> long @-> long @-> long @-> returning gen)

  let mpcosz = foreign "mpcosz" (gen @-> gen @-> returning void)
  let mpexpz = foreign "mpexpz" (gen @-> gen @-> returning void)
  let mplogz = foreign "mplogz" (gen @-> gen @-> returning void)
  let mpsinz = foreign "mpsinz" (gen @-> gen @-> returning void)
  let mul_content = foreign "mul_content" (gen @-> gen @-> returning gen)
  let mul_denom = foreign "mul_denom" (gen @-> gen @-> returning gen)
  let nbits2nlong = foreign "nbits2nlong" (long @-> returning long)
  let nbits2extraprec = foreign "nbits2extraprec" (long @-> returning long)
  let nbits2ndec = foreign "nbits2ndec" (long @-> returning long)
  let nbits2prec = foreign "nbits2prec" (long @-> returning long)
  let nbits2lg = foreign "nbits2lg" (long @-> returning long)
  let nbrows = foreign "nbrows" (gen @-> returning long)
  let nchar2nlong = foreign "nchar2nlong" (long @-> returning long)
  let ndec2nbits = foreign "ndec2nbits" (long @-> returning long)
  let ndec2nlong = foreign "ndec2nlong" (long @-> returning long)
  let ndec2prec = foreign "ndec2prec" (long @-> returning long)
  let normalize_frac = foreign "normalize_frac" (gen @-> returning void)
  let odd = foreign "odd" (long @-> returning int)
  let pari_free = foreign "pari_free" (ptr void @-> returning void)
  let pari_calloc = foreign "pari_calloc" (int @-> returning (ptr void))
  let pari_malloc = foreign "pari_malloc" (int @-> returning (ptr void))

  let pari_realloc =
    foreign "pari_realloc" (ptr void @-> int @-> returning (ptr void))

  let pari_realloc_ip =
    foreign "pari_realloc_ip" (ptr (ptr void) @-> int @-> returning void)

  let perm_conj = foreign "perm_conj" (gen @-> gen @-> returning gen)
  let perm_inv = foreign "perm_inv" (gen @-> returning gen)
  let perm_mul = foreign "perm_mul" (gen @-> gen @-> returning gen)
  let perm_sqr = foreign "perm_sqr" (gen @-> returning gen)
  let pol_0 = foreign "pol_0" (long @-> returning gen)
  let pol_1 = foreign "pol_1" (long @-> returning gen)
  let pol_x = foreign "pol_x" (long @-> returning gen)
  let pol_xn = foreign "pol_xn" (long @-> long @-> returning gen)
  let pol_xnall = foreign "pol_xnall" (long @-> long @-> returning gen)
  let pol0_f2x = foreign "pol0_F2x" (long @-> returning gen)
  let pol0_flx = foreign "pol0_Flx" (long @-> returning gen)
  let pol1_f2x = foreign "pol1_F2x" (long @-> returning gen)
  let pol1_flx = foreign "pol1_Flx" (long @-> returning gen)
  let polx_flx = foreign "polx_Flx" (long @-> returning gen)
  let polx_f2x = foreign "polx_F2x" (long @-> returning gen)
  let polx_zx = foreign "polx_zx" (long @-> returning gen)
  let polxn_flx = foreign "polxn_Flx" (long @-> long @-> returning gen)
  let powii = foreign "powii" (gen @-> gen @-> returning gen)
  let powis = foreign "powIs" (long @-> returning gen)
  let prec2nbits = foreign "prec2nbits" (long @-> returning long)

  let prec2nbits_mul =
    foreign "prec2nbits_mul" (long @-> double @-> returning double)

  let prec2ndec = foreign "prec2ndec" (long @-> returning long)
  let precdbl = foreign "precdbl" (long @-> returning long)
  let quad_disc = foreign "quad_disc" (gen @-> returning gen)
  let qfb_disc = foreign "qfb_disc" (gen @-> returning gen)
  let qfb_disc3 = foreign "qfb_disc3" (gen @-> gen @-> gen @-> returning gen)
  let quadnorm = foreign "quadnorm" (gen @-> returning gen)
  let remsbil = foreign "remsBIL" (long @-> returning long)
  let row = foreign "row" (gen @-> long @-> returning gen)
  let flm_row = foreign "Flm_row" (gen @-> long @-> returning gen)
  let row_i = foreign "row_i" (gen @-> long @-> long @-> long @-> returning gen)
  let zm_row = foreign "zm_row" (gen @-> long @-> returning gen)
  let rowcopy = foreign "rowcopy" (gen @-> long @-> returning gen)
  let rowpermute = foreign "rowpermute" (gen @-> gen @-> returning gen)
  let rowslice = foreign "rowslice" (gen @-> long @-> long @-> returning gen)

  let rowslicepermute =
    foreign "rowslicepermute" (gen @-> gen @-> long @-> long @-> returning gen)

  let rowsplice = foreign "rowsplice" (gen @-> long @-> returning gen)
  let ser_isexactzero = foreign "ser_isexactzero" (gen @-> returning int)
  let shallowcopy = foreign "shallowcopy" (gen @-> returning gen)
  let sqrfrac = foreign "sqrfrac" (gen @-> returning gen)
  let sqrti = foreign "sqrti" (gen @-> returning gen)
  let sqrtnr = foreign "sqrtnr" (gen @-> long @-> returning gen)
  let sqrtr = foreign "sqrtr" (gen @-> returning gen)
  let sstoq = foreign "sstoQ" (long @-> long @-> returning gen)
  let uutoq = foreign "uutoQ" (pari_ulong @-> pari_ulong @-> returning gen)
  let sturm = foreign "sturm" (gen @-> returning long)
  let truecoef = foreign "truecoef" (gen @-> long @-> returning gen)
  let trunc_safe = foreign "trunc_safe" (gen @-> returning gen)
  let vec_ei = foreign "vec_ei" (long @-> long @-> returning gen)
  let vec_append = foreign "vec_append" (gen @-> gen @-> returning gen)
  let vec_lengthen = foreign "vec_lengthen" (gen @-> long @-> returning gen)
  let vec_prepend = foreign "vec_prepend" (gen @-> gen @-> returning gen)
  let vec_setconst = foreign "vec_setconst" (gen @-> gen @-> returning gen)
  let vec_shorten = foreign "vec_shorten" (gen @-> long @-> returning gen)
  let vec_to_vecsmall = foreign "vec_to_vecsmall" (gen @-> returning gen)
  let vecpermute = foreign "vecpermute" (gen @-> gen @-> returning gen)
  let vecreverse = foreign "vecreverse" (gen @-> returning gen)
  let vecreverse_inplace = foreign "vecreverse_inplace" (gen @-> returning void)
  let vecsmallpermute = foreign "vecsmallpermute" (gen @-> gen @-> returning gen)
  let vecslice = foreign "vecslice" (gen @-> long @-> long @-> returning gen)

  let vecslicepermute =
    foreign "vecslicepermute" (gen @-> gen @-> long @-> long @-> returning gen)

  let vecsplice = foreign "vecsplice" (gen @-> long @-> returning gen)

  let vecsmall_append =
    foreign "vecsmall_append" (gen @-> long @-> returning gen)

  let vecsmall_coincidence =
    foreign "vecsmall_coincidence" (gen @-> gen @-> returning long)

  let vecsmall_concat = foreign "vecsmall_concat" (gen @-> gen @-> returning gen)
  let vecsmall_copy = foreign "vecsmall_copy" (gen @-> returning gen)
  let vecsmall_ei = foreign "vecsmall_ei" (long @-> long @-> returning gen)
end

module F60 (F : Ctypes.FOREIGN) = struct
  open F

  let vecsmall_indexmax = foreign "vecsmall_indexmax" (gen @-> returning long)
  let vecsmall_indexmin = foreign "vecsmall_indexmin" (gen @-> returning long)
  let vecsmall_isin = foreign "vecsmall_isin" (gen @-> long @-> returning long)

  let vecsmall_lengthen =
    foreign "vecsmall_lengthen" (gen @-> long @-> returning gen)

  let vecsmall_lexcmp = foreign "vecsmall_lexcmp" (gen @-> gen @-> returning int)
  let vecsmall_max = foreign "vecsmall_max" (gen @-> returning long)
  let vecsmall_min = foreign "vecsmall_min" (gen @-> returning long)

  let vecsmall_pack =
    foreign "vecsmall_pack" (gen @-> long @-> long @-> returning long)

  let vecsmall_prefixcmp =
    foreign "vecsmall_prefixcmp" (gen @-> gen @-> returning int)

  let vecsmall_prepend =
    foreign "vecsmall_prepend" (gen @-> long @-> returning gen)

  let vecsmall_reverse = foreign "vecsmall_reverse" (gen @-> returning gen)

  let vecsmall_shorten =
    foreign "vecsmall_shorten" (gen @-> long @-> returning gen)

  let vecsmall_to_col = foreign "vecsmall_to_col" (gen @-> returning gen)
  let vecsmall_to_vec = foreign "vecsmall_to_vec" (gen @-> returning gen)

  let vecsmall_to_vec_inplace =
    foreign "vecsmall_to_vec_inplace" (gen @-> returning gen)

  let vecsmalltrunc_append =
    foreign "vecsmalltrunc_append" (gen @-> long @-> returning void)

  let vecsmalltrunc_init = foreign "vecsmalltrunc_init" (long @-> returning gen)

  let vectrunc_append =
    foreign "vectrunc_append" (gen @-> gen @-> returning void)

  let vectrunc_append_batch =
    foreign "vectrunc_append_batch" (gen @-> gen @-> returning void)

  let vectrunc_init = foreign "vectrunc_init" (long @-> returning gen)
  let coltrunc_init = foreign "coltrunc_init" (long @-> returning gen)
  let zc_to_zc = foreign "zc_to_ZC" (gen @-> returning gen)
  let zero_f2m = foreign "zero_F2m" (long @-> long @-> returning gen)
  let zero_f2m_copy = foreign "zero_F2m_copy" (long @-> long @-> returning gen)
  let zero_f2v = foreign "zero_F2v" (long @-> returning gen)
  let zero_f2x = foreign "zero_F2x" (long @-> returning gen)
  let zero_flm = foreign "zero_Flm" (long @-> long @-> returning gen)
  let zero_flm_copy = foreign "zero_Flm_copy" (long @-> long @-> returning gen)
  let zero_flv = foreign "zero_Flv" (long @-> returning gen)
  let zero_flx = foreign "zero_Flx" (long @-> returning gen)
  let zero_zm = foreign "zero_zm" (long @-> long @-> returning gen)
  let zero_zv = foreign "zero_zv" (long @-> returning gen)
  let zero_zx = foreign "zero_zx" (long @-> returning gen)
  let zerocol = foreign "zerocol" (long @-> returning gen)
  let zeromat = foreign "zeromat" (long @-> long @-> returning gen)
  let zeromatcopy = foreign "zeromatcopy" (long @-> long @-> returning gen)
  let zeropadic = foreign "zeropadic" (gen @-> long @-> returning gen)

  let zeropadic_shallow =
    foreign "zeropadic_shallow" (gen @-> long @-> returning gen)

  let zeropol = foreign "zeropol" (long @-> returning gen)
  let zeroser = foreign "zeroser" (long @-> long @-> returning gen)
  let zerovec = foreign "zerovec" (long @-> returning gen)
  let zerovec_block = foreign "zerovec_block" (long @-> returning gen)
  let zm_copy = foreign "zm_copy" (gen @-> returning gen)
  let zm_to_zxv = foreign "zm_to_zxV" (gen @-> long @-> returning gen)
  let zm_transpose = foreign "zm_transpose" (gen @-> returning gen)
  let zv_copy = foreign "zv_copy" (gen @-> returning gen)
  let zv_to_zv = foreign "zv_to_ZV" (gen @-> returning gen)
  let zv_to_zx = foreign "zv_to_zx" (gen @-> long @-> returning gen)
  let zx_renormalize = foreign "zx_renormalize" (gen @-> long @-> returning gen)
  let zx_shift = foreign "zx_shift" (gen @-> long @-> returning gen)
  let zx_to_zv = foreign "zx_to_zv" (gen @-> long @-> returning gen)
  let err_get_compo = foreign "err_get_compo" (gen @-> long @-> returning gen)
  let err_get_num = foreign "err_get_num" (gen @-> returning long)
  let pari_err_bug = foreign "pari_err_BUG" (string @-> returning void)

  let pari_err_component =
    foreign "pari_err_COMPONENT"
      (string @-> string @-> gen @-> gen @-> returning void)

  let pari_err_constpol = foreign "pari_err_CONSTPOL" (string @-> returning void)

  let pari_err_coprime =
    foreign "pari_err_COPRIME" (string @-> gen @-> gen @-> returning void)

  let pari_err_dim = foreign "pari_err_DIM" (string @-> returning void)

  let pari_err_domain =
    foreign "pari_err_DOMAIN"
      (string @-> string @-> string @-> gen @-> gen @-> returning void)

  let pari_err_file =
    foreign "pari_err_FILE" (string @-> string @-> returning void)

  let pari_err_filedesc =
    foreign "pari_err_FILEDESC" (string @-> long @-> returning void)

  let pari_err_flag = foreign "pari_err_FLAG" (string @-> returning void)
  let pari_err_impl = foreign "pari_err_IMPL" (string @-> returning void)
  let pari_err_inv = foreign "pari_err_INV" (string @-> gen @-> returning void)

  let pari_err_irredpol =
    foreign "pari_err_IRREDPOL" (string @-> gen @-> returning void)

  let pari_err_maxprime =
    foreign "pari_err_MAXPRIME" (pari_ulong @-> returning void)

  let pari_err_modulus =
    foreign "pari_err_MODULUS" (string @-> gen @-> gen @-> returning void)

  let pari_err_op =
    foreign "pari_err_OP" (string @-> gen @-> gen @-> returning void)

  let pari_err_overflow = foreign "pari_err_OVERFLOW" (string @-> returning void)
  let pari_err_package = foreign "pari_err_PACKAGE" (string @-> returning void)
  let pari_err_prec = foreign "pari_err_PREC" (string @-> returning void)

  let pari_err_prime =
    foreign "pari_err_PRIME" (string @-> gen @-> returning void)

  let pari_err_priority =
    foreign "pari_err_PRIORITY"
      (string @-> gen @-> string @-> long @-> returning void)

  let pari_err_sqrtn =
    foreign "pari_err_SQRTN" (string @-> gen @-> returning void)

  let pari_err_type = foreign "pari_err_TYPE" (string @-> gen @-> returning void)

  let pari_err_type2 =
    foreign "pari_err_TYPE2" (string @-> gen @-> gen @-> returning void)

  let pari_err_var =
    foreign "pari_err_VAR" (string @-> gen @-> gen @-> returning void)

  let pari_err_roots0 = foreign "pari_err_ROOTS0" (string @-> returning void)
  let mkintmod = foreign "mkintmod" (gen @-> gen @-> returning gen)

  let mkintmodu =
    foreign "mkintmodu" (pari_ulong @-> pari_ulong @-> returning gen)

  let mkpolmod = foreign "mkpolmod" (gen @-> gen @-> returning gen)
  let mkfrac = foreign "mkfrac" (gen @-> gen @-> returning gen)
  let mkfracss = foreign "mkfracss" (long @-> long @-> returning gen)
  let qtoss = foreign "Qtoss" (gen @-> ptr long @-> ptr long @-> returning void)
  let sstoq = foreign "sstoQ" (long @-> long @-> returning gen)
  let uutoq = foreign "uutoQ" (pari_ulong @-> pari_ulong @-> returning gen)
  let mkfraccopy = foreign "mkfraccopy" (gen @-> gen @-> returning gen)
  let mkrfrac = foreign "mkrfrac" (gen @-> gen @-> returning gen)
  let mkrfraccopy = foreign "mkrfraccopy" (gen @-> gen @-> returning gen)
  let mkcomplex = foreign "mkcomplex" (gen @-> gen @-> returning gen)
  let gen_i = foreign "gen_I" (void @-> returning gen)
  let cgetc = foreign "cgetc" (long @-> returning gen)
  let mkquad = foreign "mkquad" (gen @-> gen @-> gen @-> returning gen)
  let mkvecsmall = foreign "mkvecsmall" (long @-> returning gen)
  let mkvecsmall2 = foreign "mkvecsmall2" (long @-> long @-> returning gen)

  let mkvecsmall3 =
    foreign "mkvecsmall3" (long @-> long @-> long @-> returning gen)

  let mkvecsmall4 =
    foreign "mkvecsmall4" (long @-> long @-> long @-> long @-> returning gen)

  let mkvecsmall5 =
    foreign "mkvecsmall5"
      (long @-> long @-> long @-> long @-> long @-> returning gen)

  let mkqfb = foreign "mkqfb" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let mkvec = foreign "mkvec" (gen @-> returning gen)
end

module F61 (F : Ctypes.FOREIGN) = struct
  open F

  let mkvec2 = foreign "mkvec2" (gen @-> gen @-> returning gen)
  let mkvec3 = foreign "mkvec3" (gen @-> gen @-> gen @-> returning gen)
  let mkvec4 = foreign "mkvec4" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkvec5 =
    foreign "mkvec5" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkvecs = foreign "mkvecs" (long @-> returning gen)
  let mkvec2s = foreign "mkvec2s" (long @-> long @-> returning gen)
  let mkvec3s = foreign "mkvec3s" (long @-> long @-> long @-> returning gen)

  let mkvec4s =
    foreign "mkvec4s" (long @-> long @-> long @-> long @-> returning gen)

  let mkveccopy = foreign "mkveccopy" (gen @-> returning gen)
  let mkvec2copy = foreign "mkvec2copy" (gen @-> gen @-> returning gen)
  let mkcol = foreign "mkcol" (gen @-> returning gen)
  let mkcol2 = foreign "mkcol2" (gen @-> gen @-> returning gen)
  let mkcol3 = foreign "mkcol3" (gen @-> gen @-> gen @-> returning gen)
  let mkcol4 = foreign "mkcol4" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkcol5 =
    foreign "mkcol5" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkcol6 =
    foreign "mkcol6"
      (gen @-> gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkcols = foreign "mkcols" (long @-> returning gen)
  let mkcol2s = foreign "mkcol2s" (long @-> long @-> returning gen)
  let mkcol3s = foreign "mkcol3s" (long @-> long @-> long @-> returning gen)

  let mkcol4s =
    foreign "mkcol4s" (long @-> long @-> long @-> long @-> returning gen)

  let mkcolcopy = foreign "mkcolcopy" (gen @-> returning gen)
  let mkmat = foreign "mkmat" (gen @-> returning gen)
  let mkmat2 = foreign "mkmat2" (gen @-> gen @-> returning gen)
  let mkmat3 = foreign "mkmat3" (gen @-> gen @-> gen @-> returning gen)
  let mkmat4 = foreign "mkmat4" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkmat5 =
    foreign "mkmat5" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let mkmatcopy = foreign "mkmatcopy" (gen @-> returning gen)
  let mkerr = foreign "mkerr" (long @-> returning gen)
  let mkoo = foreign "mkoo" (void @-> returning gen)
  let mkmoo = foreign "mkmoo" (void @-> returning gen)
  let inf_get_sign = foreign "inf_get_sign" (gen @-> returning long)

  let mkmat22s =
    foreign "mkmat22s" (long @-> long @-> long @-> long @-> returning gen)

  let mkmat22 = foreign "mkmat22" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let pol_x = foreign "pol_x" (long @-> returning gen)
  let pol_xn = foreign "pol_xn" (long @-> long @-> returning gen)
  let pol_xnall = foreign "pol_xnall" (long @-> long @-> returning gen)
  let polxn_flx = foreign "polxn_Flx" (long @-> long @-> returning gen)
  let pol_1 = foreign "pol_1" (long @-> returning gen)
  let pol_0 = foreign "pol_0" (long @-> returning gen)
  let const_vec = foreign "const_vec" (long @-> gen @-> returning gen)
  let const_col = foreign "const_col" (long @-> gen @-> returning gen)
  let const_vecsmall = foreign "const_vecsmall" (long @-> long @-> returning gen)
  let zeropadic = foreign "zeropadic" (gen @-> long @-> returning gen)

  let zeropadic_shallow =
    foreign "zeropadic_shallow" (gen @-> long @-> returning gen)

  let zeroser = foreign "zeroser" (long @-> long @-> returning gen)
  let ser_isexactzero = foreign "ser_isexactzero" (gen @-> returning int)
  let zeropol = foreign "zeropol" (long @-> returning gen)
  let zerocol = foreign "zerocol" (long @-> returning gen)
  let zerovec = foreign "zerovec" (long @-> returning gen)
  let zeromat = foreign "zeromat" (long @-> long @-> returning gen)
  let zero_flx = foreign "zero_Flx" (long @-> returning gen)
  let zero_flv = foreign "zero_Flv" (long @-> returning gen)
  let zero_flm = foreign "zero_Flm" (long @-> long @-> returning gen)
  let zero_flm_copy = foreign "zero_Flm_copy" (long @-> long @-> returning gen)
  let zero_f2v = foreign "zero_F2v" (long @-> returning gen)
  let zero_f2m = foreign "zero_F2m" (long @-> long @-> returning gen)
  let zero_f2m_copy = foreign "zero_F2m_copy" (long @-> long @-> returning gen)
  let zeromatcopy = foreign "zeromatcopy" (long @-> long @-> returning gen)
  let zerovec_block = foreign "zerovec_block" (long @-> returning gen)
  let col_ei = foreign "col_ei" (long @-> long @-> returning gen)
  let vec_ei = foreign "vec_ei" (long @-> long @-> returning gen)
  let f2v_ei = foreign "F2v_ei" (long @-> long @-> returning gen)
  let vecsmall_ei = foreign "vecsmall_ei" (long @-> long @-> returning gen)
  let rg_col_ei = foreign "Rg_col_ei" (gen @-> long @-> long @-> returning gen)
  let shallowcopy = foreign "shallowcopy" (gen @-> returning gen)
  let vectrunc_init = foreign "vectrunc_init" (long @-> returning gen)
  let coltrunc_init = foreign "coltrunc_init" (long @-> returning gen)
  let lg_increase = foreign "lg_increase" (gen @-> returning void)

  let vectrunc_append =
    foreign "vectrunc_append" (gen @-> gen @-> returning void)

  let vectrunc_append_batch =
    foreign "vectrunc_append_batch" (gen @-> gen @-> returning void)

  let vecsmalltrunc_init = foreign "vecsmalltrunc_init" (long @-> returning gen)

  let vecsmalltrunc_append =
    foreign "vecsmalltrunc_append" (gen @-> long @-> returning void)

  let hash_str = foreign "hash_str" (string @-> returning pari_ulong)

  let hash_str_len =
    foreign "hash_str_len" (string @-> long @-> returning pari_ulong)

  let vec_shorten = foreign "vec_shorten" (gen @-> long @-> returning gen)
  let vec_lengthen = foreign "vec_lengthen" (gen @-> long @-> returning gen)
  let vec_append = foreign "vec_append" (gen @-> gen @-> returning gen)
  let vec_prepend = foreign "vec_prepend" (gen @-> gen @-> returning gen)
  let vec_setconst = foreign "vec_setconst" (gen @-> gen @-> returning gen)

  let vecsmall_shorten =
    foreign "vecsmall_shorten" (gen @-> long @-> returning gen)

  let vecsmall_lengthen =
    foreign "vecsmall_lengthen" (gen @-> long @-> returning gen)

  let vec_to_vecsmall = foreign "vec_to_vecsmall" (gen @-> returning gen)
  let vecsmall_to_vec = foreign "vecsmall_to_vec" (gen @-> returning gen)

  let vecsmall_to_vec_inplace =
    foreign "vecsmall_to_vec_inplace" (gen @-> returning gen)

  let vecsmall_to_col = foreign "vecsmall_to_col" (gen @-> returning gen)
  let vecsmall_lexcmp = foreign "vecsmall_lexcmp" (gen @-> gen @-> returning int)

  let vecsmall_prefixcmp =
    foreign "vecsmall_prefixcmp" (gen @-> gen @-> returning int)

  let vecsmall_prepend =
    foreign "vecsmall_prepend" (gen @-> long @-> returning gen)

  let vecsmall_append =
    foreign "vecsmall_append" (gen @-> long @-> returning gen)

  let vecsmall_concat = foreign "vecsmall_concat" (gen @-> gen @-> returning gen)

  let vecsmall_coincidence =
    foreign "vecsmall_coincidence" (gen @-> gen @-> returning long)

  let vecsmall_isin = foreign "vecsmall_isin" (gen @-> long @-> returning long)

  let vecsmall_pack =
    foreign "vecsmall_pack" (gen @-> long @-> long @-> returning long)

  let vecsmall_indexmax = foreign "vecsmall_indexmax" (gen @-> returning long)
  let vecsmall_max = foreign "vecsmall_max" (gen @-> returning long)
  let vecsmall_indexmin = foreign "vecsmall_indexmin" (gen @-> returning long)
  let vecsmall_min = foreign "vecsmall_min" (gen @-> returning long)
  let zv_isscalar = foreign "ZV_isscalar" (gen @-> returning int)
  let qv_isscalar = foreign "QV_isscalar" (gen @-> returning int)
  let rgv_isscalar = foreign "RgV_isscalar" (gen @-> returning int)
end

module F62 (F : Ctypes.FOREIGN) = struct
  open F

  let rgx_isscalar = foreign "RgX_isscalar" (gen @-> returning int)
  let rgx_equal_var = foreign "RgX_equal_var" (gen @-> gen @-> returning long)
  let rgx_to_rgv = foreign "RgX_to_RgV" (gen @-> long @-> returning gen)
  let rgx_is_rational = foreign "RgX_is_rational" (gen @-> returning int)
  let rgx_is_zx = foreign "RgX_is_ZX" (gen @-> returning int)
  let rgx_is_qx = foreign "RgX_is_QX" (gen @-> returning int)
  let rgx_is_monomial = foreign "RgX_is_monomial" (gen @-> returning int)
  let rgv_is_zv = foreign "RgV_is_ZV" (gen @-> returning int)
  let rgv_is_qv = foreign "RgV_is_QV" (gen @-> returning int)
  let rgv_isin_i = foreign "RgV_isin_i" (gen @-> gen @-> long @-> returning long)
  let rgv_isin = foreign "RgV_isin" (gen @-> gen @-> returning long)
  let vecslice = foreign "vecslice" (gen @-> long @-> long @-> returning gen)

  let vecslicepermute =
    foreign "vecslicepermute" (gen @-> gen @-> long @-> long @-> returning gen)

  let rowslicepermute =
    foreign "rowslicepermute" (gen @-> gen @-> long @-> long @-> returning gen)

  let rowslice = foreign "rowslice" (gen @-> long @-> long @-> returning gen)

  let matslice =
    foreign "matslice"
      (gen @-> long @-> long @-> long @-> long @-> returning gen)

  let rowsplice = foreign "rowsplice" (gen @-> long @-> returning gen)
  let vecsplice = foreign "vecsplice" (gen @-> long @-> returning gen)
  let rgm_minor = foreign "RgM_minor" (gen @-> long @-> long @-> returning gen)
  let row = foreign "row" (gen @-> long @-> returning gen)
  let flm_row = foreign "Flm_row" (gen @-> long @-> returning gen)
  let rowcopy = foreign "rowcopy" (gen @-> long @-> returning gen)
  let row_i = foreign "row_i" (gen @-> long @-> long @-> long @-> returning gen)
  let vecreverse = foreign "vecreverse" (gen @-> returning gen)
  let vecsmall_reverse = foreign "vecsmall_reverse" (gen @-> returning gen)
  let vecreverse_inplace = foreign "vecreverse_inplace" (gen @-> returning void)
  let vecsmallpermute = foreign "vecsmallpermute" (gen @-> gen @-> returning gen)
  let vecpermute = foreign "vecpermute" (gen @-> gen @-> returning gen)
  let rowpermute = foreign "rowpermute" (gen @-> gen @-> returning gen)
  let identity_zv = foreign "identity_zv" (long @-> returning gen)
  let identity_zv = foreign "identity_ZV" (long @-> returning gen)
  let identity_perm = foreign "identity_perm" (long @-> returning gen)
  let cyclic_perm = foreign "cyclic_perm" (long @-> long @-> returning gen)
  let perm_mul = foreign "perm_mul" (gen @-> gen @-> returning gen)
  let perm_sqr = foreign "perm_sqr" (gen @-> returning gen)
  let perm_inv = foreign "perm_inv" (gen @-> returning gen)
  let perm_conj = foreign "perm_conj" (gen @-> gen @-> returning gen)
  let pari_free = foreign "pari_free" (ptr void @-> returning void)
  let pari_malloc = foreign "pari_malloc" (int @-> returning (ptr void))

  let pari_realloc =
    foreign "pari_realloc" (ptr void @-> int @-> returning (ptr void))

  let pari_realloc_ip =
    foreign "pari_realloc_ip" (ptr (ptr void) @-> int @-> returning void)

  let pari_calloc = foreign "pari_calloc" (int @-> returning (ptr void))
  let cgetalloc = foreign "cgetalloc" (int @-> long @-> returning gen)
  let icopy_avma = foreign "icopy_avma" (gen @-> pari_sp @-> returning gen)
  let leafcopy_avma = foreign "leafcopy_avma" (gen @-> pari_sp @-> returning gen)

  let gerepileuptoleaf =
    foreign "gerepileuptoleaf" (pari_sp @-> gen @-> returning gen)

  let gerepileuptoint =
    foreign "gerepileuptoint" (pari_sp @-> gen @-> returning gen)

  let gerepileupto = foreign "gerepileupto" (pari_sp @-> gen @-> returning gen)
  let gerepilecopy = foreign "gerepilecopy" (pari_sp @-> gen @-> returning gen)
  let gunclonenull = foreign "guncloneNULL" (gen @-> returning void)
  let gunclonenull_deep = foreign "guncloneNULL_deep" (gen @-> returning void)

  let gerepilemany =
    foreign "gerepilemany" (pari_sp @-> ptr (ptr gen) @-> int @-> returning void)

  let gerepileall = foreign "gerepileall" (pari_sp @-> int @-> returning void)
  let gc_all = foreign "gc_all" (pari_sp @-> int @-> returning gen)

  let gerepilecoeffs =
    foreign "gerepilecoeffs" (pari_sp @-> gen @-> int @-> returning void)

  let bin_copy = foreign "bin_copy" (ptr genbin @-> returning gen)
  let genbinbase = foreign "GENbinbase" (ptr genbin @-> returning gen)
  let cgiv = foreign "cgiv" (gen @-> returning void)
  let killblock = foreign "killblock" (gen @-> returning void)

  let is_universal_constant =
    foreign "is_universal_constant" (gen @-> returning int)

  let cxcompotor = foreign "cxcompotor" (gen @-> long @-> returning gen)
  let cxtofp = foreign "cxtofp" (gen @-> long @-> returning gen)
  let cxtoreal = foreign "cxtoreal" (gen @-> returning gen)
  let gtodouble = foreign "gtodouble" (gen @-> returning double)
  let gisdouble = foreign "gisdouble" (gen @-> ptr double @-> returning int)
  let gtos = foreign "gtos" (gen @-> returning long)
  let gtou = foreign "gtou" (gen @-> returning pari_ulong)
  let absfrac = foreign "absfrac" (gen @-> returning gen)
  let absfrac_shallow = foreign "absfrac_shallow" (gen @-> returning gen)
  let q_abs = foreign "Q_abs" (gen @-> returning gen)
  let q_abs_shallow = foreign "Q_abs_shallow" (gen @-> returning gen)
  let r_abs_shallow = foreign "R_abs_shallow" (gen @-> returning gen)
  let r_abs = foreign "R_abs" (gen @-> returning gen)
  let gtofp = foreign "gtofp" (gen @-> long @-> returning gen)
  let gtomp = foreign "gtomp" (gen @-> long @-> returning gen)
  let rgx_gtofp = foreign "RgX_gtofp" (gen @-> long @-> returning gen)
  let rgc_gtofp = foreign "RgC_gtofp" (gen @-> long @-> returning gen)
  let rgv_gtofp = foreign "RgV_gtofp" (gen @-> long @-> returning gen)
  let rgm_gtofp = foreign "RgM_gtofp" (gen @-> long @-> returning gen)
  let rgc_gtomp = foreign "RgC_gtomp" (gen @-> long @-> returning gen)
  let rgm_gtomp = foreign "RgM_gtomp" (gen @-> long @-> returning gen)
  let rgx_fpnorml2 = foreign "RgX_fpnorml2" (gen @-> long @-> returning gen)
  let rgc_fpnorml2 = foreign "RgC_fpnorml2" (gen @-> long @-> returning gen)
  let rgm_fpnorml2 = foreign "RgM_fpnorml2" (gen @-> long @-> returning gen)
  let affgr = foreign "affgr" (gen @-> gen @-> returning void)
  let affc_fixlg = foreign "affc_fixlg" (gen @-> gen @-> returning gen)
  let trunc_safe = foreign "trunc_safe" (gen @-> returning gen)
  let ndec2nlong = foreign "ndec2nlong" (long @-> returning long)
  let ndec2prec = foreign "ndec2prec" (long @-> returning long)
  let ndec2nbits = foreign "ndec2nbits" (long @-> returning long)
  let nbits2nlong = foreign "nbits2nlong" (long @-> returning long)
  let nbits2extraprec = foreign "nbits2extraprec" (long @-> returning long)
  let nbits2prec = foreign "nbits2prec" (long @-> returning long)
  let nbits2lg = foreign "nbits2lg" (long @-> returning long)
  let nchar2nlong = foreign "nchar2nlong" (long @-> returning long)
  let prec2nbits = foreign "prec2nbits" (long @-> returning long)

  let bit_accuracy_mul =
    foreign "bit_accuracy_mul" (long @-> double @-> returning double)

  let prec2nbits_mul =
    foreign "prec2nbits_mul" (long @-> double @-> returning double)

  let bit_prec = foreign "bit_prec" (gen @-> returning long)
  let bit_accuracy = foreign "bit_accuracy" (long @-> returning long)
end

module F63 (F : Ctypes.FOREIGN) = struct
  open F

  let prec2ndec = foreign "prec2ndec" (long @-> returning long)
  let nbits2ndec = foreign "nbits2ndec" (long @-> returning long)
  let precdbl = foreign "precdbl" (long @-> returning long)
  let divsbil = foreign "divsBIL" (long @-> returning long)
  let remsbil = foreign "remsBIL" (long @-> returning long)
  let fp_red = foreign "Fp_red" (gen @-> gen @-> returning gen)
  let fp_add = foreign "Fp_add" (gen @-> gen @-> gen @-> returning gen)
  let fp_sub = foreign "Fp_sub" (gen @-> gen @-> gen @-> returning gen)
  let fp_neg = foreign "Fp_neg" (gen @-> gen @-> returning gen)
  let fp_halve = foreign "Fp_halve" (gen @-> gen @-> returning gen)
  let fp_center = foreign "Fp_center" (gen @-> gen @-> gen @-> returning gen)
  let fp_center_i = foreign "Fp_center_i" (gen @-> gen @-> gen @-> returning gen)

  let fp_addmul =
    foreign "Fp_addmul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fp_mul = foreign "Fp_mul" (gen @-> gen @-> gen @-> returning gen)
  let fp_sqr = foreign "Fp_sqr" (gen @-> gen @-> returning gen)
  let fp_mulu = foreign "Fp_mulu" (gen @-> pari_ulong @-> gen @-> returning gen)
  let fp_muls = foreign "Fp_muls" (gen @-> long @-> gen @-> returning gen)
  let fp_inv = foreign "Fp_inv" (gen @-> gen @-> returning gen)
  let fp_invsafe = foreign "Fp_invsafe" (gen @-> gen @-> returning gen)
  let fp_div = foreign "Fp_div" (gen @-> gen @-> gen @-> returning gen)
  let fp_divu = foreign "Fp_divu" (gen @-> pari_ulong @-> gen @-> returning gen)

  let flx_mulu =
    foreign "Flx_mulu" (gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let get_f2x_mod = foreign "get_F2x_mod" (gen @-> returning gen)
  let get_f2x_var = foreign "get_F2x_var" (gen @-> returning long)
  let get_f2x_degree = foreign "get_F2x_degree" (gen @-> returning long)
  let get_f2xqx_mod = foreign "get_F2xqX_mod" (gen @-> returning gen)
  let get_f2xqx_var = foreign "get_F2xqX_var" (gen @-> returning long)
  let get_f2xqx_degree = foreign "get_F2xqX_degree" (gen @-> returning long)
  let get_flx_mod = foreign "get_Flx_mod" (gen @-> returning gen)
  let get_flx_var = foreign "get_Flx_var" (gen @-> returning long)
  let get_flx_degree = foreign "get_Flx_degree" (gen @-> returning long)
  let get_flxqx_mod = foreign "get_FlxqX_mod" (gen @-> returning gen)
  let get_flxqx_var = foreign "get_FlxqX_var" (gen @-> returning long)
  let get_flxqx_degree = foreign "get_FlxqX_degree" (gen @-> returning long)
  let get_fpx_mod = foreign "get_FpX_mod" (gen @-> returning gen)
  let get_fpx_var = foreign "get_FpX_var" (gen @-> returning long)
  let get_fpx_degree = foreign "get_FpX_degree" (gen @-> returning long)
  let get_fpxqx_mod = foreign "get_FpXQX_mod" (gen @-> returning gen)
  let get_fpxqx_var = foreign "get_FpXQX_var" (gen @-> returning long)
  let get_fpxqx_degree = foreign "get_FpXQX_degree" (gen @-> returning long)
  let submulii = foreign "submulii" (gen @-> gen @-> gen @-> returning gen)
  let mulsubii = foreign "mulsubii" (gen @-> gen @-> gen @-> returning gen)

  let submuliu =
    foreign "submuliu" (gen @-> gen @-> pari_ulong @-> returning gen)

  let addmuliu =
    foreign "addmuliu" (gen @-> gen @-> pari_ulong @-> returning gen)

  let submuliu_inplace =
    foreign "submuliu_inplace" (gen @-> gen @-> pari_ulong @-> returning gen)

  let addmuliu_inplace =
    foreign "addmuliu_inplace" (gen @-> gen @-> pari_ulong @-> returning gen)

  let lincombii =
    foreign "lincombii" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let is_const_t = foreign "is_const_t" (long @-> returning int)
  let is_extscalar_t = foreign "is_extscalar_t" (long @-> returning int)
  let is_intreal_t = foreign "is_intreal_t" (long @-> returning int)
  let is_matvec_t = foreign "is_matvec_t" (long @-> returning int)
  let is_noncalc_t = foreign "is_noncalc_t" (long @-> returning int)
  let is_qfb_t = foreign "is_qfb_t" (long @-> returning int)
  let is_rational_t = foreign "is_rational_t" (long @-> returning int)
  let is_real_t = foreign "is_real_t" (long @-> returning int)
  let is_recursive_t = foreign "is_recursive_t" (long @-> returning int)
  let is_scalar_t = foreign "is_scalar_t" (long @-> returning int)
  let is_vec_t = foreign "is_vec_t" (long @-> returning int)
  let qfb_is_qfi = foreign "qfb_is_qfi" (gen @-> returning int)
  let sqrtr = foreign "sqrtr" (gen @-> returning gen)
  let cbrtr_abs = foreign "cbrtr_abs" (gen @-> returning gen)
  let cbrtr = foreign "cbrtr" (gen @-> returning gen)
  let sqrtnr = foreign "sqrtnr" (gen @-> long @-> returning gen)
  let logint = foreign "logint" (gen @-> gen @-> returning long)

  let ulogint =
    foreign "ulogint" (pari_ulong @-> pari_ulong @-> returning pari_ulong)

  let ismpzero = foreign "ismpzero" (gen @-> returning int)
  let isintzero = foreign "isintzero" (gen @-> returning int)
  let isint1 = foreign "isint1" (gen @-> returning int)
  let isintm1 = foreign "isintm1" (gen @-> returning int)
  let equali1 = foreign "equali1" (gen @-> returning int)
  let equalim1 = foreign "equalim1" (gen @-> returning int)
  let is_pm1 = foreign "is_pm1" (gen @-> returning int)
  let is_bigint = foreign "is_bigint" (gen @-> returning int)
  let odd = foreign "odd" (long @-> returning int)
  let both_odd = foreign "both_odd" (long @-> long @-> returning int)
  let isonstack = foreign "isonstack" (gen @-> returning int)
  let dbllog2r = foreign "dbllog2r" (gen @-> returning double)
  let mul_content = foreign "mul_content" (gen @-> gen @-> returning gen)
  let inv_content = foreign "inv_content" (gen @-> returning gen)
  let div_content = foreign "div_content" (gen @-> gen @-> returning gen)
  let mul_denom = foreign "mul_denom" (gen @-> gen @-> returning gen)
  let constant_coeff = foreign "constant_coeff" (gen @-> returning gen)
  let leading_coeff = foreign "leading_coeff" (gen @-> returning gen)
  let flx_lead = foreign "Flx_lead" (gen @-> returning pari_ulong)
  let flx_constant = foreign "Flx_constant" (gen @-> returning pari_ulong)
  let degpol = foreign "degpol" (gen @-> returning long)
  let lgpol = foreign "lgpol" (gen @-> returning long)
  let lgcols = foreign "lgcols" (gen @-> returning long)
  let nbrows = foreign "nbrows" (gen @-> returning long)
  let truecoef = foreign "truecoef" (gen @-> long @-> returning gen)
  let zxq_mul = foreign "ZXQ_mul" (gen @-> gen @-> gen @-> returning gen)
  let zxq_sqr = foreign "ZXQ_sqr" (gen @-> gen @-> returning gen)
  let rgx_copy = foreign "RgX_copy" (gen @-> returning gen)
  let rgx_coeff = foreign "RgX_coeff" (gen @-> long @-> returning gen)
  let rgx_renormalize = foreign "RgX_renormalize" (gen @-> returning gen)
  let rgx_div = foreign "RgX_div" (gen @-> gen @-> returning gen)
  let rgxqx_div = foreign "RgXQX_div" (gen @-> gen @-> gen @-> returning gen)
  let rgxqx_rem = foreign "RgXQX_rem" (gen @-> gen @-> gen @-> returning gen)
  let fpx_div = foreign "FpX_div" (gen @-> gen @-> gen @-> returning gen)
  let flx_div = foreign "Flx_div" (gen @-> gen @-> pari_ulong @-> returning gen)
end

module F64 (F : Ctypes.FOREIGN) = struct
  open F

  let flx_div_pre =
    foreign "Flx_div_pre"
      (gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let f2x_div = foreign "F2x_div" (gen @-> gen @-> returning gen)
  let fpv_fpc_mul = foreign "FpV_FpC_mul" (gen @-> gen @-> gen @-> returning gen)
  let pol0_flx = foreign "pol0_Flx" (long @-> returning gen)
  let pol1_flx = foreign "pol1_Flx" (long @-> returning gen)
  let polx_flx = foreign "polx_Flx" (long @-> returning gen)
  let zero_zx = foreign "zero_zx" (long @-> returning gen)
  let polx_zx = foreign "polx_zx" (long @-> returning gen)
  let zx_shift = foreign "zx_shift" (gen @-> long @-> returning gen)
  let zx_renormalize = foreign "zx_renormalize" (gen @-> long @-> returning gen)
  let zero_f2x = foreign "zero_F2x" (long @-> returning gen)
  let pol0_f2x = foreign "pol0_F2x" (long @-> returning gen)
  let pol1_f2x = foreign "pol1_F2x" (long @-> returning gen)
  let polx_f2x = foreign "polx_F2x" (long @-> returning gen)
  let f2x_equal1 = foreign "F2x_equal1" (gen @-> returning int)
  let f2x_equal = foreign "F2x_equal" (gen @-> gen @-> returning int)
  let f2x_copy = foreign "F2x_copy" (gen @-> returning gen)
  let f2v_copy = foreign "F2v_copy" (gen @-> returning gen)
  let flv_copy = foreign "Flv_copy" (gen @-> returning gen)
  let flx_copy = foreign "Flx_copy" (gen @-> returning gen)
  let vecsmall_copy = foreign "vecsmall_copy" (gen @-> returning gen)
  let flx_equal1 = foreign "Flx_equal1" (gen @-> returning int)
  let zx_equal1 = foreign "ZX_equal1" (gen @-> returning int)
  let zx_is_monic = foreign "ZX_is_monic" (gen @-> returning int)
  let zx_renormalize = foreign "ZX_renormalize" (gen @-> long @-> returning gen)

  let fpx_renormalize =
    foreign "FpX_renormalize" (gen @-> long @-> returning gen)

  let fpxx_renormalize =
    foreign "FpXX_renormalize" (gen @-> long @-> returning gen)

  let fpxqx_renormalize =
    foreign "FpXQX_renormalize" (gen @-> long @-> returning gen)

  let f2x_renormalize =
    foreign "F2x_renormalize" (gen @-> long @-> returning gen)

  let f2v_to_f2x = foreign "F2v_to_F2x" (gen @-> long @-> returning gen)
  let sturm = foreign "sturm" (gen @-> returning long)
  let gval = foreign "gval" (gen @-> long @-> returning long)

  let rgx_shift_inplace_init =
    foreign "RgX_shift_inplace_init" (long @-> returning void)

  let rgx_shift_inplace =
    foreign "RgX_shift_inplace" (gen @-> long @-> returning gen)

  let zv_to_zv = foreign "zv_to_ZV" (gen @-> returning gen)
  let zc_to_zc = foreign "zc_to_ZC" (gen @-> returning gen)
  let zv_to_zv = foreign "ZV_to_zv" (gen @-> returning gen)
  let zx_to_zv = foreign "zx_to_zv" (gen @-> long @-> returning gen)
  let zv_to_zx = foreign "zv_to_zx" (gen @-> long @-> returning gen)
  let zm_to_zxv = foreign "zm_to_zxV" (gen @-> long @-> returning gen)
  let zero_zm = foreign "zero_zm" (long @-> long @-> returning gen)
  let zero_zv = foreign "zero_zv" (long @-> returning gen)
  let zm_transpose = foreign "zm_transpose" (gen @-> returning gen)
  let zm_copy = foreign "zm_copy" (gen @-> returning gen)
  let zv_copy = foreign "zv_copy" (gen @-> returning gen)
  let zm_row = foreign "zm_row" (gen @-> long @-> returning gen)
  let zc_hnfrem = foreign "ZC_hnfrem" (gen @-> gen @-> returning gen)
  let zm_hnfrem = foreign "ZM_hnfrem" (gen @-> gen @-> returning gen)
  let zm_lll = foreign "ZM_lll" (gen @-> double @-> long @-> returning gen)

  let rgm_dimensions =
    foreign "RgM_dimensions" (gen @-> ptr long @-> ptr long @-> returning void)

  let rgm_shallowcopy = foreign "RgM_shallowcopy" (gen @-> returning gen)
  let f2m_copy = foreign "F2m_copy" (gen @-> returning gen)
  let f3m_copy = foreign "F3m_copy" (gen @-> returning gen)
  let flm_copy = foreign "Flm_copy" (gen @-> returning gen)
  let zv_dvd = foreign "ZV_dvd" (gen @-> gen @-> returning int)
  let zm_zv_mod = foreign "ZM_ZV_mod" (gen @-> gen @-> returning gen)
  let zv_zv_mod = foreign "ZV_ZV_mod" (gen @-> gen @-> returning gen)
  let vecmodii = foreign "vecmodii" (gen @-> gen @-> returning gen)
  let vecmoduu = foreign "vecmoduu" (gen @-> gen @-> returning gen)
  let fq_red = foreign "Fq_red" (gen @-> gen @-> gen @-> returning gen)
  let fq_to_fpxq = foreign "Fq_to_FpXQ" (gen @-> gen @-> gen @-> returning gen)
  let rg_to_fq = foreign "Rg_to_Fq" (gen @-> gen @-> gen @-> returning gen)

  let gener_fq_local =
    foreign "gener_Fq_local" (gen @-> gen @-> gen @-> returning gen)

  let fpxqx_div =
    foreign "FpXQX_div" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let flxqx_div =
    foreign "FlxqX_div" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxqx_div_pre =
    foreign "FlxqX_div_pre"
      (gen @-> gen @-> gen @-> pari_ulong @-> pari_ulong @-> returning gen)

  let f2xqx_div = foreign "F2xqX_div" (gen @-> gen @-> gen @-> returning gen)

  let fpxy_fq_evaly =
    foreign "FpXY_Fq_evaly"
      (gen @-> gen @-> gen @-> gen @-> long @-> returning gen)

  let fqx_red = foreign "FqX_red" (gen @-> gen @-> gen @-> returning gen)
  let fqx_add = foreign "FqX_add" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fqx_neg = foreign "FqX_neg" (gen @-> gen @-> gen @-> returning gen)
  let fqx_sub = foreign "FqX_sub" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_fp_mul =
    foreign "FqX_Fp_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_fq_mul =
    foreign "FqX_Fq_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_mul = foreign "FqX_mul" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_mulu =
    foreign "FqX_mulu" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fqx_sqr = foreign "FqX_sqr" (gen @-> gen @-> gen @-> returning gen)

  let fqx_powu =
    foreign "FqX_powu" (gen @-> pari_ulong @-> gen @-> gen @-> returning gen)

  let fqx_halve = foreign "FqX_halve" (gen @-> gen @-> gen @-> returning gen)
  let fqx_div = foreign "FqX_div" (gen @-> gen @-> gen @-> gen @-> returning gen)
  let fqx_get_red = foreign "FqX_get_red" (gen @-> gen @-> gen @-> returning gen)
  let fqx_rem = foreign "FqX_rem" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_divrem =
    foreign "FqX_divrem"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fqx_div_by_x_x =
    foreign "FqX_div_by_X_x"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> returning gen)

  let fqx_halfgcd =
    foreign "FqX_halfgcd" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_gcd = foreign "FqX_gcd" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqx_extgcd =
    foreign "FqX_extgcd"
      (gen @-> gen @-> gen @-> gen @-> ptr gen @-> ptr gen @-> returning gen)

  let fqx_normalize =
    foreign "FqX_normalize" (gen @-> gen @-> gen @-> returning gen)

  let fqx_deriv = foreign "FqX_deriv" (gen @-> gen @-> gen @-> returning gen)
  let fqx_integ = foreign "FqX_integ" (gen @-> gen @-> gen @-> returning gen)
  let fqx_factor = foreign "FqX_factor" (gen @-> gen @-> gen @-> returning gen)

  let fqx_factor_squarefree =
    foreign "FqX_factor_squarefree" (gen @-> gen @-> gen @-> returning gen)

  let fqx_ddf = foreign "FqX_ddf" (gen @-> gen @-> gen @-> returning gen)
  let fqx_degfact = foreign "FqX_degfact" (gen @-> gen @-> gen @-> returning gen)
  let fqx_roots = foreign "FqX_roots" (gen @-> gen @-> gen @-> returning gen)
  let fqx_to_mod = foreign "FqX_to_mod" (gen @-> gen @-> gen @-> returning gen)

  let fqxq_add =
    foreign "FqXQ_add" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_sub =
    foreign "FqXQ_sub" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_div =
    foreign "FqXQ_div" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_inv =
    foreign "FqXQ_inv" (gen @-> gen @-> gen @-> gen @-> returning gen)
end

module F65 (F : Ctypes.FOREIGN) = struct
  open F

  let fqxq_invsafe =
    foreign "FqXQ_invsafe" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_mul =
    foreign "FqXQ_mul" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_sqr =
    foreign "FqXQ_sqr" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxq_pow =
    foreign "FqXQ_pow" (gen @-> gen @-> gen @-> gen @-> gen @-> returning gen)

  let fqxn_expint =
    foreign "FqXn_expint" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fqxn_exp =
    foreign "FqXn_exp" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fqxn_inv =
    foreign "FqXn_inv" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fqxn_mul =
    foreign "FqXn_mul" (gen @-> gen @-> long @-> gen @-> gen @-> returning gen)

  let fqxn_sqr =
    foreign "FqXn_sqr" (gen @-> long @-> gen @-> gen @-> returning gen)

  let fpxq_add =
    foreign "FpXQ_add" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let fpxq_sub =
    foreign "FpXQ_sub" (gen @-> gen @-> gen @-> gen @-> returning gen)

  let flxq_add =
    foreign "Flxq_add" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let flxq_sub =
    foreign "Flxq_sub" (gen @-> gen @-> gen @-> pari_ulong @-> returning gen)

  let f2x_coeff = foreign "F2x_coeff" (gen @-> long @-> returning pari_ulong)
  let f2x_clear = foreign "F2x_clear" (gen @-> long @-> returning void)
  let f2x_set = foreign "F2x_set" (gen @-> long @-> returning void)
  let f2x_flip = foreign "F2x_flip" (gen @-> long @-> returning void)
  let f2v_coeff = foreign "F2v_coeff" (gen @-> long @-> returning pari_ulong)
  let f2v_clear = foreign "F2v_clear" (gen @-> long @-> returning void)
  let f2v_set = foreign "F2v_set" (gen @-> long @-> returning void)
  let f2v_flip = foreign "F2v_flip" (gen @-> long @-> returning void)

  let f2m_coeff =
    foreign "F2m_coeff" (gen @-> long @-> long @-> returning pari_ulong)

  let f2m_clear = foreign "F2m_clear" (gen @-> long @-> long @-> returning void)
  let f2m_set = foreign "F2m_set" (gen @-> long @-> long @-> returning void)
  let f2m_flip = foreign "F2m_flip" (gen @-> long @-> long @-> returning void)

  let f3m_coeff =
    foreign "F3m_coeff" (gen @-> long @-> long @-> returning pari_ulong)

  let f3m_set =
    foreign "F3m_set" (gen @-> long @-> long @-> pari_ulong @-> returning void)

  let matpascal = foreign "matpascal" (long @-> returning gen)
  let z_issquare = foreign "Z_issquare" (gen @-> returning long)
  let z_ispower = foreign "Z_ispower" (gen @-> pari_ulong @-> returning long)
  let sqrti = foreign "sqrti" (gen @-> returning gen)
  let gaddgs = foreign "gaddgs" (gen @-> long @-> returning gen)
  let gcmpgs = foreign "gcmpgs" (gen @-> long @-> returning int)
  let gequalgs = foreign "gequalgs" (gen @-> long @-> returning int)
  let gmaxsg = foreign "gmaxsg" (long @-> gen @-> returning gen)
  let gminsg = foreign "gminsg" (long @-> gen @-> returning gen)
  let gmulgs = foreign "gmulgs" (gen @-> long @-> returning gen)
  let gmulgu = foreign "gmulgu" (gen @-> pari_ulong @-> returning gen)
  let gsubgs = foreign "gsubgs" (gen @-> long @-> returning gen)
  let gdivsg = foreign "gdivsg" (long @-> gen @-> returning gen)
  let gmax_shallow = foreign "gmax_shallow" (gen @-> gen @-> returning gen)
  let gmin_shallow = foreign "gmin_shallow" (gen @-> gen @-> returning gen)
  let cxnorm = foreign "cxnorm" (gen @-> returning gen)
  let quadnorm = foreign "quadnorm" (gen @-> returning gen)
  let quad_disc = foreign "quad_disc" (gen @-> returning gen)
  let qfb_disc3 = foreign "qfb_disc3" (gen @-> gen @-> gen @-> returning gen)
  let qfb_disc = foreign "qfb_disc" (gen @-> returning gen)
  let sqrfrac = foreign "sqrfrac" (gen @-> returning gen)
  let normalize_frac = foreign "normalize_frac" (gen @-> returning void)
  let powii = foreign "powii" (gen @-> gen @-> returning gen)
  let powis = foreign "powIs" (long @-> returning gen)
  let mpexpz = foreign "mpexpz" (gen @-> gen @-> returning void)
  let mplogz = foreign "mplogz" (gen @-> gen @-> returning void)
  let mpcosz = foreign "mpcosz" (gen @-> gen @-> returning void)
  let mpsinz = foreign "mpsinz" (gen @-> gen @-> returning void)
  let gnegz = foreign "gnegz" (gen @-> gen @-> returning void)
  let gabsz = foreign "gabsz" (gen @-> long @-> gen @-> returning void)
  let gaddz = foreign "gaddz" (gen @-> gen @-> gen @-> returning void)
  let gsubz = foreign "gsubz" (gen @-> gen @-> gen @-> returning void)
  let gmulz = foreign "gmulz" (gen @-> gen @-> gen @-> returning void)
  let gdivz = foreign "gdivz" (gen @-> gen @-> gen @-> returning void)
  let gdiventz = foreign "gdiventz" (gen @-> gen @-> gen @-> returning void)
  let gmodz = foreign "gmodz" (gen @-> gen @-> gen @-> returning void)
  let gmul2nz = foreign "gmul2nz" (gen @-> long @-> gen @-> returning void)
  let gshiftz = foreign "gshiftz" (gen @-> long @-> gen @-> returning void)
  let ell_get_a1 = foreign "ell_get_a1" (gen @-> returning gen)
  let ell_get_a2 = foreign "ell_get_a2" (gen @-> returning gen)
  let ell_get_a3 = foreign "ell_get_a3" (gen @-> returning gen)
  let ell_get_a4 = foreign "ell_get_a4" (gen @-> returning gen)
  let ell_get_a6 = foreign "ell_get_a6" (gen @-> returning gen)
  let ell_get_b2 = foreign "ell_get_b2" (gen @-> returning gen)
  let ell_get_b4 = foreign "ell_get_b4" (gen @-> returning gen)
  let ell_get_b6 = foreign "ell_get_b6" (gen @-> returning gen)
  let ell_get_b8 = foreign "ell_get_b8" (gen @-> returning gen)
  let ell_get_c4 = foreign "ell_get_c4" (gen @-> returning gen)
  let ell_get_c6 = foreign "ell_get_c6" (gen @-> returning gen)
  let ell_get_disc = foreign "ell_get_disc" (gen @-> returning gen)
  let ell_get_j = foreign "ell_get_j" (gen @-> returning gen)
  let ell_get_type = foreign "ell_get_type" (gen @-> returning long)
  let ellff_get_field = foreign "ellff_get_field" (gen @-> returning gen)
  let ellff_get_a4a6 = foreign "ellff_get_a4a6" (gen @-> returning gen)
  let ellqp_get_zero = foreign "ellQp_get_zero" (gen @-> returning gen)
  let ellqp_get_prec = foreign "ellQp_get_prec" (gen @-> returning long)
  let ellqp_get_p = foreign "ellQp_get_p" (gen @-> returning gen)
  let ellr_get_prec = foreign "ellR_get_prec" (gen @-> returning long)
  let ellr_get_sign = foreign "ellR_get_sign" (gen @-> returning long)
  let ellnf_get_nf = foreign "ellnf_get_nf" (gen @-> returning gen)
  let ellnf_get_bnf = foreign "ellnf_get_bnf" (gen @-> returning gen)
  let checkell_i = foreign "checkell_i" (gen @-> returning int)
  let ell_is_inf = foreign "ell_is_inf" (gen @-> returning int)
  let ellinf = foreign "ellinf" (void @-> returning gen)
  let modpr_get_pr = foreign "modpr_get_pr" (gen @-> returning gen)
  let modpr_get_p = foreign "modpr_get_p" (gen @-> returning gen)
  let modpr_get_t = foreign "modpr_get_T" (gen @-> returning gen)
  let pr_get_p = foreign "pr_get_p" (gen @-> returning gen)
  let pr_get_gen = foreign "pr_get_gen" (gen @-> returning gen)
  let pr_get_e = foreign "pr_get_e" (gen @-> returning long)
  let pr_get_f = foreign "pr_get_f" (gen @-> returning long)
  let pr_get_tau = foreign "pr_get_tau" (gen @-> returning gen)
  let pr_is_inert = foreign "pr_is_inert" (gen @-> returning int)
end

module F66 (F : Ctypes.FOREIGN) = struct
  open F

  let pr_norm = foreign "pr_norm" (gen @-> returning gen)
  let upr_norm = foreign "upr_norm" (gen @-> returning pari_ulong)
  let nf_get_varn = foreign "nf_get_varn" (gen @-> returning long)
  let nf_get_pol = foreign "nf_get_pol" (gen @-> returning gen)
  let nf_get_degree = foreign "nf_get_degree" (gen @-> returning long)
  let nf_get_r1 = foreign "nf_get_r1" (gen @-> returning long)
  let nf_get_r2 = foreign "nf_get_r2" (gen @-> returning long)
  let nf_get_disc = foreign "nf_get_disc" (gen @-> returning gen)
  let nf_get_index = foreign "nf_get_index" (gen @-> returning gen)
  let nf_get_m = foreign "nf_get_M" (gen @-> returning gen)
  let nf_get_g = foreign "nf_get_G" (gen @-> returning gen)
  let nf_get_roundg = foreign "nf_get_roundG" (gen @-> returning gen)
  let nf_get_tr = foreign "nf_get_Tr" (gen @-> returning gen)
  let nf_get_diff = foreign "nf_get_diff" (gen @-> returning gen)

  let nf_get_ramified_primes =
    foreign "nf_get_ramified_primes" (gen @-> returning gen)

  let nf_get_roots = foreign "nf_get_roots" (gen @-> returning gen)
  let nf_get_zk = foreign "nf_get_zk" (gen @-> returning gen)
  let nf_get_zkprimpart = foreign "nf_get_zkprimpart" (gen @-> returning gen)
  let nf_get_zkden = foreign "nf_get_zkden" (gen @-> returning gen)
  let nf_get_invzk = foreign "nf_get_invzk" (gen @-> returning gen)

  let nf_get_sign =
    foreign "nf_get_sign" (gen @-> ptr long @-> ptr long @-> returning void)

  let cyc_get_expo = foreign "cyc_get_expo" (gen @-> returning gen)
  let abgrp_get_no = foreign "abgrp_get_no" (gen @-> returning gen)
  let abgrp_get_cyc = foreign "abgrp_get_cyc" (gen @-> returning gen)
  let abgrp_get_gen = foreign "abgrp_get_gen" (gen @-> returning gen)
  let bnf_get_nf = foreign "bnf_get_nf" (gen @-> returning gen)
  let bnf_get_clgp = foreign "bnf_get_clgp" (gen @-> returning gen)
  let bnf_get_no = foreign "bnf_get_no" (gen @-> returning gen)
  let bnf_get_cyc = foreign "bnf_get_cyc" (gen @-> returning gen)
  let bnf_get_gen = foreign "bnf_get_gen" (gen @-> returning gen)
  let bnf_get_reg = foreign "bnf_get_reg" (gen @-> returning gen)
  let bnf_get_logfu = foreign "bnf_get_logfu" (gen @-> returning gen)
  let bnf_get_sunits = foreign "bnf_get_sunits" (gen @-> returning gen)
  let bnf_get_tuu = foreign "bnf_get_tuU" (gen @-> returning gen)
  let bnf_get_tun = foreign "bnf_get_tuN" (gen @-> returning long)
  let bnf_get_fu_nocheck = foreign "bnf_get_fu_nocheck" (gen @-> returning gen)

  let nfv_to_scalar_or_alg =
    foreign "nfV_to_scalar_or_alg" (gen @-> gen @-> returning gen)

  let bnf_get_fu = foreign "bnf_get_fu" (gen @-> returning gen)
  let bnr_get_bnf = foreign "bnr_get_bnf" (gen @-> returning gen)
  let bnr_get_bid = foreign "bnr_get_bid" (gen @-> returning gen)
  let bnr_get_mod = foreign "bnr_get_mod" (gen @-> returning gen)
  let bnr_get_nf = foreign "bnr_get_nf" (gen @-> returning gen)
  let bnr_get_clgp = foreign "bnr_get_clgp" (gen @-> returning gen)
  let bnr_get_no = foreign "bnr_get_no" (gen @-> returning gen)
  let bnr_get_cyc = foreign "bnr_get_cyc" (gen @-> returning gen)
  let bnr_get_gen_nocheck = foreign "bnr_get_gen_nocheck" (gen @-> returning gen)
  let bnr_get_gen = foreign "bnr_get_gen" (gen @-> returning gen)
  let locs_get_cyc = foreign "locs_get_cyc" (gen @-> returning gen)
  let locs_get_lsprk = foreign "locs_get_Lsprk" (gen @-> returning gen)
  let locs_get_lgenfil = foreign "locs_get_Lgenfil" (gen @-> returning gen)
  let locs_get_mod = foreign "locs_get_mod" (gen @-> returning gen)
  let locs_get_famod = foreign "locs_get_famod" (gen @-> returning gen)
  let locs_get_m_infty = foreign "locs_get_m_infty" (gen @-> returning gen)
  let gchar_get_basis = foreign "gchar_get_basis" (gen @-> returning gen)
  let gchar_get_bnf = foreign "gchar_get_bnf" (gen @-> returning gen)
  let gchar_get_nf = foreign "gchar_get_nf" (gen @-> returning gen)
  let gchar_get_zm = foreign "gchar_get_zm" (gen @-> returning gen)
  let gchar_get_mod = foreign "gchar_get_mod" (gen @-> returning gen)
  let gchar_get_modp = foreign "gchar_get_modP" (gen @-> returning gen)
  let gchar_get_s = foreign "gchar_get_S" (gen @-> returning gen)
  let gchar_get_dldata = foreign "gchar_get_DLdata" (gen @-> returning gen)
  let gchar_get_sfu = foreign "gchar_get_sfu" (gen @-> returning gen)
  let gchar_get_cyc = foreign "gchar_get_cyc" (gen @-> returning gen)
  let gchar_get_hnf = foreign "gchar_get_hnf" (gen @-> returning gen)
  let gchar_get_u = foreign "gchar_get_U" (gen @-> returning gen)
  let gchar_get_ui = foreign "gchar_get_Ui" (gen @-> returning gen)
  let gchar_get_m0 = foreign "gchar_get_m0" (gen @-> returning gen)
  let gchar_get_u0 = foreign "gchar_get_u0" (gen @-> returning gen)
  let gchar_get_r1 = foreign "gchar_get_r1" (gen @-> returning long)
  let gchar_get_r2 = foreign "gchar_get_r2" (gen @-> returning long)
  let gchar_get_loccyc = foreign "gchar_get_loccyc" (gen @-> returning gen)
  let gchar_get_nc = foreign "gchar_get_nc" (gen @-> returning long)
  let gchar_get_ns = foreign "gchar_get_ns" (gen @-> returning long)
  let gchar_get_nm = foreign "gchar_get_nm" (gen @-> returning long)
  let gchar_get_evalprec = foreign "gchar_get_evalprec" (gen @-> returning long)
  let gchar_get_prec = foreign "gchar_get_prec" (gen @-> returning long)
  let gchar_get_nfprec = foreign "gchar_get_nfprec" (gen @-> returning long)

  let gchar_set_evalprec =
    foreign "gchar_set_evalprec" (gen @-> long @-> returning void)

  let gchar_set_prec = foreign "gchar_set_prec" (gen @-> long @-> returning void)

  let gchar_copy_precs =
    foreign "gchar_copy_precs" (gen @-> gen @-> returning void)

  let gchar_set_nfprec =
    foreign "gchar_set_nfprec" (gen @-> long @-> returning void)

  let gchar_get_ntors = foreign "gchar_get_ntors" (gen @-> returning long)
  let gchar_get_nfree = foreign "gchar_get_nfree" (gen @-> returning long)
  let gchar_get_nalg = foreign "gchar_get_nalg" (gen @-> returning long)

  let gchar_set_basis =
    foreign "gchar_set_basis" (gen @-> gen @-> returning void)

  let gchar_set_nf = foreign "gchar_set_nf" (gen @-> gen @-> returning void)

  let gchar_set_ntors =
    foreign "gchar_set_ntors" (gen @-> long @-> returning void)

  let gchar_set_nfree =
    foreign "gchar_set_nfree" (gen @-> long @-> returning void)

  let gchar_set_nalg = foreign "gchar_set_nalg" (gen @-> long @-> returning void)
  let gchar_set_cyc = foreign "gchar_set_cyc" (gen @-> gen @-> returning void)

  let gchar_set_huui =
    foreign "gchar_set_HUUi" (gen @-> gen @-> gen @-> gen @-> returning void)

  let gchar_set_m0 = foreign "gchar_set_m0" (gen @-> gen @-> returning void)
  let gchar_set_u0 = foreign "gchar_set_u0" (gen @-> gen @-> returning void)
  let bid_get_mod = foreign "bid_get_mod" (gen @-> returning gen)
  let bid_get_ideal = foreign "bid_get_ideal" (gen @-> returning gen)
  let bid_get_arch = foreign "bid_get_arch" (gen @-> returning gen)
  let bid_get_grp = foreign "bid_get_grp" (gen @-> returning gen)
  let bid_get_fact = foreign "bid_get_fact" (gen @-> returning gen)
  let bid_get_fact2 = foreign "bid_get_fact2" (gen @-> returning gen)
  let bid_get_sprk = foreign "bid_get_sprk" (gen @-> returning gen)
end

module F67 (F : Ctypes.FOREIGN) = struct
  open F

  let bid_get_sarch = foreign "bid_get_sarch" (gen @-> returning gen)
  let bid_get_archp = foreign "bid_get_archp" (gen @-> returning gen)
  let bid_get_u = foreign "bid_get_U" (gen @-> returning gen)
  let bid_get_no = foreign "bid_get_no" (gen @-> returning gen)
  let bid_get_cyc = foreign "bid_get_cyc" (gen @-> returning gen)
  let bid_get_gen_nocheck = foreign "bid_get_gen_nocheck" (gen @-> returning gen)
  let bid_get_gen = foreign "bid_get_gen" (gen @-> returning gen)
  let znstar_get_n = foreign "znstar_get_N" (gen @-> returning gen)
  let znstar_get_fan = foreign "znstar_get_faN" (gen @-> returning gen)
  let znstar_get_no = foreign "znstar_get_no" (gen @-> returning gen)
  let znstar_get_cyc = foreign "znstar_get_cyc" (gen @-> returning gen)
  let znstar_get_gen = foreign "znstar_get_gen" (gen @-> returning gen)

  let znstar_get_conreycyc =
    foreign "znstar_get_conreycyc" (gen @-> returning gen)

  let znstar_get_conreygen =
    foreign "znstar_get_conreygen" (gen @-> returning gen)

  let znstar_get_ui = foreign "znstar_get_Ui" (gen @-> returning gen)
  let znstar_get_u = foreign "znstar_get_U" (gen @-> returning gen)
  let znstar_get_pe = foreign "znstar_get_pe" (gen @-> returning gen)
  let gal_get_pol = foreign "gal_get_pol" (gen @-> returning gen)
  let gal_get_p = foreign "gal_get_p" (gen @-> returning gen)
  let gal_get_e = foreign "gal_get_e" (gen @-> returning gen)
  let gal_get_mod = foreign "gal_get_mod" (gen @-> returning gen)
  let gal_get_roots = foreign "gal_get_roots" (gen @-> returning gen)
  let gal_get_invvdm = foreign "gal_get_invvdm" (gen @-> returning gen)
  let gal_get_den = foreign "gal_get_den" (gen @-> returning gen)
  let gal_get_group = foreign "gal_get_group" (gen @-> returning gen)
  let gal_get_gen = foreign "gal_get_gen" (gen @-> returning gen)
  let gal_get_orders = foreign "gal_get_orders" (gen @-> returning gen)
  let rnf_get_degree = foreign "rnf_get_degree" (gen @-> returning long)
  let rnf_get_nfdegree = foreign "rnf_get_nfdegree" (gen @-> returning long)
  let rnf_get_absdegree = foreign "rnf_get_absdegree" (gen @-> returning long)
  let rnf_get_idealdisc = foreign "rnf_get_idealdisc" (gen @-> returning gen)
  let rnf_get_k = foreign "rnf_get_k" (gen @-> returning gen)
  let rnf_get_alpha = foreign "rnf_get_alpha" (gen @-> returning gen)
  let rnf_get_nf = foreign "rnf_get_nf" (gen @-> returning gen)
  let rnf_get_nfzk = foreign "rnf_get_nfzk" (gen @-> returning gen)
  let rnf_get_polabs = foreign "rnf_get_polabs" (gen @-> returning gen)
  let rnf_get_pol = foreign "rnf_get_pol" (gen @-> returning gen)
  let rnf_get_disc = foreign "rnf_get_disc" (gen @-> returning gen)
  let rnf_get_index = foreign "rnf_get_index" (gen @-> returning gen)

  let rnf_get_ramified_primes =
    foreign "rnf_get_ramified_primes" (gen @-> returning gen)

  let rnf_get_varn = foreign "rnf_get_varn" (gen @-> returning long)
  let rnf_get_nfpol = foreign "rnf_get_nfpol" (gen @-> returning gen)
  let rnf_get_nfvarn = foreign "rnf_get_nfvarn" (gen @-> returning long)
  let rnf_get_zk = foreign "rnf_get_zk" (gen @-> returning gen)
  let rnf_get_map = foreign "rnf_get_map" (gen @-> returning gen)
  let rnf_get_invzk = foreign "rnf_get_invzk" (gen @-> returning gen)
  let idealred = foreign "idealred" (gen @-> gen @-> returning gen)

  let idealchineseinit =
    foreign "idealchineseinit" (gen @-> gen @-> returning gen)

  let closure_arity = foreign "closure_arity" (gen @-> returning long)

  let closure_is_variadic =
    foreign "closure_is_variadic" (gen @-> returning long)

  let closure_codestr = foreign "closure_codestr" (gen @-> returning string)
  let closure_get_code = foreign "closure_get_code" (gen @-> returning gen)
  let closure_get_oper = foreign "closure_get_oper" (gen @-> returning gen)
  let closure_get_data = foreign "closure_get_data" (gen @-> returning gen)
  let closure_get_dbg = foreign "closure_get_dbg" (gen @-> returning gen)
  let closure_get_text = foreign "closure_get_text" (gen @-> returning gen)
  let closure_get_frame = foreign "closure_get_frame" (gen @-> returning gen)
  let err_get_num = foreign "err_get_num" (gen @-> returning long)
  let err_get_compo = foreign "err_get_compo" (gen @-> long @-> returning gen)
  let pari_err_bug = foreign "pari_err_BUG" (string @-> returning void)
  let pari_err_constpol = foreign "pari_err_CONSTPOL" (string @-> returning void)

  let pari_err_coprime =
    foreign "pari_err_COPRIME" (string @-> gen @-> gen @-> returning void)

  let pari_err_dim = foreign "pari_err_DIM" (string @-> returning void)

  let pari_err_file =
    foreign "pari_err_FILE" (string @-> string @-> returning void)

  let pari_err_filedesc =
    foreign "pari_err_FILEDESC" (string @-> long @-> returning void)

  let pari_err_flag = foreign "pari_err_FLAG" (string @-> returning void)
  let pari_err_impl = foreign "pari_err_IMPL" (string @-> returning void)
  let pari_err_inv = foreign "pari_err_INV" (string @-> gen @-> returning void)

  let pari_err_irredpol =
    foreign "pari_err_IRREDPOL" (string @-> gen @-> returning void)

  let pari_err_domain =
    foreign "pari_err_DOMAIN"
      (string @-> string @-> string @-> gen @-> gen @-> returning void)

  let pari_err_component =
    foreign "pari_err_COMPONENT"
      (string @-> string @-> gen @-> gen @-> returning void)

  let pari_err_maxprime =
    foreign "pari_err_MAXPRIME" (pari_ulong @-> returning void)

  let pari_err_op =
    foreign "pari_err_OP" (string @-> gen @-> gen @-> returning void)

  let pari_err_overflow = foreign "pari_err_OVERFLOW" (string @-> returning void)
  let pari_err_prec = foreign "pari_err_PREC" (string @-> returning void)
  let pari_err_package = foreign "pari_err_PACKAGE" (string @-> returning void)

  let pari_err_prime =
    foreign "pari_err_PRIME" (string @-> gen @-> returning void)

  let pari_err_modulus =
    foreign "pari_err_MODULUS" (string @-> gen @-> gen @-> returning void)

  let pari_err_roots0 = foreign "pari_err_ROOTS0" (string @-> returning void)

  let pari_err_sqrtn =
    foreign "pari_err_SQRTN" (string @-> gen @-> returning void)

  let pari_err_type = foreign "pari_err_TYPE" (string @-> gen @-> returning void)

  let pari_err_type2 =
    foreign "pari_err_TYPE2" (string @-> gen @-> gen @-> returning void)

  let pari_err_var =
    foreign "pari_err_VAR" (string @-> gen @-> gen @-> returning void)

  let pari_err_priority =
    foreign "pari_err_PRIORITY"
      (string @-> gen @-> string @-> long @-> returning void)
end

module Functions (F : Ctypes.FOREIGN) = struct
  include F0 (F)
  include F1 (F)
  include F2 (F)
  include F3 (F)
  include F4 (F)
  include F5 (F)
  include F6 (F)
  include F7 (F)
  include F8 (F)
  include F9 (F)
  include F10 (F)
  include F11 (F)
  include F12 (F)
  include F13 (F)
  include F14 (F)
  include F15 (F)
  include F16 (F)
  include F17 (F)
  include F18 (F)
  include F19 (F)
  include F20 (F)
  include F21 (F)
  include F22 (F)
  include F23 (F)
  include F24 (F)
  include F25 (F)
  include F26 (F)
  include F27 (F)
  include F28 (F)
  include F29 (F)
  include F30 (F)
  include F31 (F)
  include F32 (F)
  include F33 (F)
  include F34 (F)
  include F35 (F)
  include F36 (F)
  include F37 (F)
  include F38 (F)
  include F39 (F)
  include F40 (F)
  include F41 (F)
  include F42 (F)
  include F43 (F)
  include F44 (F)
  include F45 (F)
  include F46 (F)
  include F47 (F)
  include F48 (F)
  include F49 (F)
  include F50 (F)
  include F51 (F)
  include F52 (F)
  include F53 (F)
  include F54 (F)
  include F55 (F)
  include F56 (F)
  include F57 (F)
  include F58 (F)
  include F59 (F)
  include F60 (F)
  include F61 (F)
  include F62 (F)
  include F63 (F)
  include F64 (F)
  include F65 (F)
  include F66 (F)
end
