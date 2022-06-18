MTX_BEGIN

bool mtx_elem_pos(uint64_t &ln, uint64_t &col, uint64_t idx, uint64_t curr_col_cnt, uint64_t ln_from, uint64_t col_from, uint64_t ln_dilate, uint64_t col_dilate);
bool mtx_elem_pos(uint64_t &ln, uint64_t &col, uint64_t idx, uint64_t col_cnt);
uint64_t mtx_elem_pos(uint64_t ln, uint64_t col, uint64_t orgn_col_cnt, uint64_t ln_from, uint64_t col_from, uint64_t ln_dilate, uint64_t col_dilate);
uint64_t mtx_elem_pos(uint64_t ln, uint64_t col, uint64_t col_cnt);

__mtx_callback bool mtx_copy(__mtx &ans, __mtx &src, uint64_t elem_cnt);
__mtx_callback bool mtx_copy(__mtx &ans, __mtx &src, uint64_t ln_cnt, uint64_t col_cnt);
__mtx_callback bool mtx_move(__mtx &ans, __mtx &&src);

__mtx_callback bool mtx_fill(__mtx &src, __mtx_elem &&val, uint64_t elem_cnt);
__mtx_callback bool mtx_fill(__mtx &src, __mtx_elem &&val, uint64_t ln_cnt, uint64_t col_cnt);

__mtx_callback bool mtx_init_E(__mtx &src, uint64_t dms);
__mtx_callback bool mtx_init_rand(__mtx &src, uint64_t elem_cnt, __mtx_elem &&boundary_first = 0, __mtx_elem &&boundary_second = 0, __mtx_elem &&acc = 1e-5);
__mtx_callback bool mtx_print(__mtx &src, uint64_t ln_cnt, uint64_t col_cnt);
__mtx_callback bool mtx_equal(__mtx &first, __mtx &second, uint64_t elem_cnt);
__mtx_callback bool mtx_equal(__mtx &first, __mtx &second, uint64_t ln_cnt, uint64_t col_cnt);

uint64_t mtx_pad_cnt(uint64_t prev_pad, uint64_t rear_pad, uint64_t dir_cnt, uint64_t dir_distance);
__mtx_callback bool mtx_pad(__mtx &ans, uint64_t ans_col_cnt, __mtx &src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_t = 0, uint64_t col_r = 0, uint64_t ln_b = 0, uint64_t col_l = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0);
uint64_t mtx_crop_cnt(uint64_t prev_crop, uint64_t rear_crop, uint64_t dir_cnt, uint64_t dir_distance);
__mtx_callback bool mtx_crop(__mtx &ans, uint64_t ans_ln_cnt, uint64_t ans_col_cnt, __mtx &src, uint64_t ln_cnt, uint64_t ln_t = 0, uint64_t col_r = 0, uint64_t ln_b = 0, uint64_t col_l = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0);

__mtx_callback __mtx_elem mtx_extm_val(__mtx &src, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_cnt, uint64_t col_cnt, uint64_t *&ln_ls, uint64_t *&col_ls, uint64_t &pos_ls_len, bool max_flag = true, bool get_extm_pos = true, uint64_t ln_dilation = 0, uint64_t col_dilation = 0);
__mtx_callback __mtx_elem mtx_sum(__mtx &src, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_dilation = 0, uint64_t col_dilation = 0);

__mtx_callback bool mtx_add(__mtx &ans, __mtx &l_src, __mtx &r_src, uint64_t ln_cnt, uint64_t col_cnt, bool nega = false);
__mtx_callback bool mtx_mult(__mtx &ans, __mtx &l_src, __mtx &r_src, uint64_t l_ln_cnt, uint64_t l_col_cnt, uint64_t r_ln_cnt, uint64_t r_col_cnt);
__mtx_callback bool mtx_mult(__mtx &ans, __mtx &src, __mtx_elem &&val, uint64_t ln_cnt, uint64_t col_cnt);
__mtx_callback bool mtx_elem_cal_opt(__mtx &ans, __mtx &l_src, __mtx &r_src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t opt_idx);
__mtx_callback bool mtx_elem_cal_opt(__mtx &ans, __mtx &src, __mtx_elem &&para, uint64_t ln_cnt, uint64_t col_cnt, uint64_t opt_idx);
__mtx_callback bool mtx_broadcast_add(__mtx &ans, __mtx &src, uint64_t elem_cnt, __mtx_elem &&val, bool nega = false);

__mtx_callback bool mtx_abs(__mtx &src, uint64_t elem_cnt);
__mtx_callback bool mtx_transposition(__mtx &ans, __mtx &src, uint64_t ln_cnt, uint64_t col_cnt);
__mtx_callback bool mtx_swap_elem(__mtx &ans, __mtx &src, uint64_t l_pos, uint64_t r_pos, uint64_t ln_cnt, uint64_t col_cnt, bool is_ln = true);
__mtx_callback __mtx_elem mtx_det(__mtx &src, uint64_t dms);
__mtx_callback bool mtx_cofactor(__mtx &ans, __mtx &src, uint64_t ln, uint64_t col, uint64_t dms);
__mtx_callback __mtx_elem mtx_algebraic_cofactor(__mtx &src, uint64_t ln, uint64_t col, uint64_t dms);
__mtx_callback bool mtx_adjugate(__mtx &ans, __mtx &src, uint64_t dms);
__mtx_callback bool mtx_inverser(__mtx &ans, __mtx &src, uint64_t dms);
__mtx_callback __mtx_elem mtx_max_eigenvalue(__mtx &src, __mtx &w, uint64_t dms, __mtx_elem &&acc = 1e-5, __mtx_elem &&init_elem = 1);
__mtx_callback uint64_t mtx_rank(__mtx &src, uint64_t ln_cnt, uint64_t col_cnt);

__mtx_callback bool mtx_LU(__mtx &ans, __mtx &src, uint64_t dms);
__mtx_callback bool mtx_equation(__mtx &ans, __mtx &coefficient, __mtx &b, uint64_t dms);
__mtx_callback bool mtx_jacobi_iterate(__mtx &ans, __mtx &coefficient, __mtx &b, uint64_t dms, __mtx_elem &&acc = 1e-5, __mtx_elem &&init_elem = 1);

__mtx_callback bool mtx_rotate_rect(__mtx &ans, __mtx &src, int ln_cnt, int col_cnt, bool clock_wise = true);
__mtx_callback bool mtx_mirror_flip(__mtx &ans, __mtx &src, uint64_t ln_cnt, uint64_t col_cnt, bool symmetry_vertical = true);

MTX_END

#include "matrix_base.hpp"