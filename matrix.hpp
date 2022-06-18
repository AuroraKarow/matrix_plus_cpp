MTX_BEGIN

__mtx_callback_ inline int __mtx_t::mtx_strassen::ln_curr(int ln) { return ((ln + ln_begin) * col_cnt_orgn + col_begin); }
__mtx_callback_ inline void __mtx_t::mtx_strassen::alloc(int ln_cnt, int col_cnt)
{
    auto elem_cnt = ln_cnt * col_cnt;
    MTX_ALLOC(__mtx_elem, ptr, elem_cnt);
    mtx_fill(ptr, std::move(__mtx_elem(0)), elem_cnt);
    this->ln_cnt = ln_cnt; this->ln_cnt_orgn = ln_cnt;
    this->col_cnt = col_cnt; this->col_cnt_orgn = col_cnt;
}
__mtx_callback_ inline bool __mtx_t::mtx_strassen::convert(__mtx &ans, bool free)
{
    if(ans)
    {
        for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j) ans[mtx_elem_pos(i, j, col_cnt)] = ptr[ln_curr(i)+j];
        if(free) release();
        return true;
    }
    return false;
}
__mtx_callback_ inline __mtx __mtx_t::mtx_strassen::operator[](int ln) { return ptr + ln_curr(ln); }
__mtx_callback_ inline typename __mtx_t::mtx_strassen __mtx_t::mtx_strassen::operator+(mtx_strassen &right)
{
    mtx_strassen ans; ans.alloc(ln_cnt, col_cnt);
    for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j) ans.ptr[ans.ln_curr(i)+j] = ptr[ln_curr(i)+j] + right.ptr[right.ln_curr(i)+j];
    return ans;
}
__mtx_callback_ inline typename __mtx_t::mtx_strassen __mtx_t::mtx_strassen::operator-(mtx_strassen &right)
{
    mtx_strassen ans; ans.alloc(ln_cnt, col_cnt);
    for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j) ans.ptr[ans.ln_curr(i)+j] = ptr[ln_curr(i)+j] - right.ptr[right.ln_curr(i)+j];
    return ans;
}
__mtx_callback_ inline typename __mtx_t::mtx_strassen __mtx_t::mtx_strassen::operator*(mtx_strassen &right)
{
    mtx_strassen ans; ans.alloc(ln_cnt, right.col_cnt);
    for(int i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j)
    {
        auto coe = ptr[ln_curr(i)+j];
        for(auto k=0; k<ans.col_cnt; ++k) ans.ptr[ans.ln_curr(i)+k] += coe * right.ptr[right.ln_curr(j)+k];
    }
    return ans;
}
__mtx_callback_ inline bool __mtx_t::mtx_strassen::release()
{
    if(ln_cnt==ln_cnt_orgn && col_cnt==col_cnt_orgn && ptr)
    {
        MTX_RESET(ptr);
        return true;
    }
    else return false;
}
__mtx_callback_ inline void __mtx_t::mtx_strassen::print()
{
    for(auto i=0; i<ln_cnt; ++i)
    {
        for(auto j=0; j<col_cnt; ++j) std::cout << ptr[ln_curr(i)+j] << '\t';
        std::cout << std::endl;
    }
}
__mtx_callback_ inline void __mtx_t::mtx_strassen_child(mtx_strassen &ans, mtx_strassen &src, uint64_t ln_from, uint64_t col_from, uint64_t ln_cnt_child, uint64_t col_cnt_child)
{
    ans.ptr = src.ptr; ans.ln_cnt_orgn = src.ln_cnt_orgn; ans.col_cnt_orgn = src.col_cnt_orgn;
    ans.ln_cnt = ln_cnt_child;
    ans.col_cnt = col_cnt_child;
    ans.ln_begin = src.ln_begin + ln_from; ans.ln_end = src.ln_begin + ln_cnt_child - 1;
    ans.col_begin = src.col_begin + col_from; ans.ln_end = src.col_begin + col_cnt_child - 1;
}
__mtx_callback_ inline void __mtx_t::mtx_strassen_quartile(mtx_strassen &src, mtx_strassen &src00, mtx_strassen &src01, mtx_strassen &src10, mtx_strassen &src11)
{
    auto dms_sub = src.ln_cnt / 2;
    mtx_strassen_child(src00, src, 0, 0, dms_sub, dms_sub);
    mtx_strassen_child(src01, src, 0, dms_sub, dms_sub, dms_sub);
    mtx_strassen_child(src10, src, dms_sub, 0, dms_sub, dms_sub);
    mtx_strassen_child(src11, src, dms_sub, dms_sub, dms_sub, dms_sub);
}
__mtx_callback_ inline typename __mtx_t::mtx_strassen __mtx_t::mtx_strassen_mult(mtx_strassen &left, mtx_strassen &right, int recursive_gate)
{
    if(left.ln_cnt > recursive_gate)
    {
        mtx_strassen a00, a01, a10, a11, b00, b01, b10, b11;
        
        mtx_strassen_quartile(left, a00, a01, a10, a11);
        mtx_strassen_quartile(right, b00, b01, b10, b11);

        auto s0 = b01 - b11;
        auto s1 = a00 + a01;
        auto s2 = a10 + a11;
        auto s3 = b10 - b00;
        auto s4 = a00 + a11;
        auto s5 = b00 + b11;
        auto s6 = a01 - a11;
        auto s7 = b10 + b11;
        auto s8 = a00 - a10;
        auto s9 = b00 + b01;

        auto p0 = mtx_strassen_mult(a00, s0, recursive_gate);
        auto p1 = mtx_strassen_mult(s1, b11, recursive_gate);
        auto p2 = mtx_strassen_mult(s2, b00, recursive_gate);
        auto p3 = mtx_strassen_mult(a11, s3, recursive_gate);
        auto p4 = mtx_strassen_mult(s4, s5, recursive_gate);
        auto p5 = mtx_strassen_mult(s6, s7, recursive_gate);
        auto p6 = mtx_strassen_mult(s8, s9, recursive_gate);

        s0.release(); s1.release(); s2.release(); s3.release(); s4.release(); s5.release(); s6.release(); s7.release(); s8.release(); s9.release();

        mtx_strassen ans; ans.alloc(left.ln_cnt, right.col_cnt);
        auto dms_sub = ans.ln_cnt / 2;
        for(auto i=0; i<ans.ln_cnt; ++i) for(auto j=0; j<ans.col_cnt; ++j)
        {
            if(i<dms_sub)
            {
                if(j<dms_sub) ans.ptr[ans.ln_curr(i)+j] = p4[i][j] + p3[i][j] - p1[i][j] + p5[i][j];
                else ans.ptr[ans.ln_curr(i)+j] = p0[i][j-dms_sub] + p1[i][j-dms_sub];
            }
            else
            {
                if(j<dms_sub) ans.ptr[ans.ln_curr(i)+j] = p2[i-dms_sub][j] + p3[i-dms_sub][j];
                else ans.ptr[ans.ln_curr(i)+j] = p4[i-dms_sub][j-dms_sub] + p0[i-dms_sub][j-dms_sub] - p2[i-dms_sub][j-dms_sub] - p6[i-dms_sub][j-dms_sub];
            }
        }

        p0.release(); p1.release(); p2.release(); p3.release(); p4.release(); p5.release(); p6.release();
        return ans;
    }
    else return left * right;
}
__mtx_callback_ inline bool __mtx_t::mtx_strassen_mult(__mtx &ans, __mtx &l_src, __mtx &r_src, uint64_t l_ln_cnt, uint64_t l_col_cnt, uint64_t r_ln_cnt, uint64_t r_col_cnt, uint64_t recursive_gate)
{
    if(ans && l_src && r_src && l_ln_cnt && l_col_cnt && r_ln_cnt && r_col_cnt && recursive_gate)
    {
        auto r_p_l = l_col_cnt, b_p_l = l_ln_cnt, r_p_r = r_col_cnt, b_p_r = r_ln_cnt, mx_len = _BAGRT num_extreme({r_p_l, b_p_l, r_p_r, b_p_r});
        if(mx_len > recursive_gate)
        {
            auto pad_val = bagrt::num_pow_pad_cnt(mx_len, 2), pad_ans = mx_len + pad_val;
            r_p_l = pad_ans - r_p_l; b_p_l = pad_ans - b_p_l;
            r_p_r = pad_ans - r_p_r; b_p_r = pad_ans - b_p_r;
            auto l_ln_cnt_p = mtx_pad_cnt(0, b_p_l, l_ln_cnt, 0),
                l_col_cnt_p = mtx_pad_cnt(0, r_p_l, l_col_cnt, 0),
                r_ln_cnt_p = mtx_pad_cnt(0, b_p_r, r_ln_cnt, 0),
                r_col_cnt_p = mtx_pad_cnt(0, r_p_r, r_col_cnt, 0);
            mtx_strassen l_src_s, r_src_s;
            l_src_s.alloc(l_ln_cnt_p, l_col_cnt_p);
            r_src_s.alloc(r_ln_cnt_p, r_col_cnt_p);
            if((mtx_pad(l_src_s.ptr, l_col_cnt_p, l_src, l_ln_cnt, l_col_cnt, 0, r_p_l, b_p_l) && mtx_pad(r_src_s.ptr, r_col_cnt_p, r_src, r_ln_cnt, r_col_cnt, 0, r_p_r, b_p_r)))
            {
                l_src_s.ln_cnt = pad_ans; l_src_s.ln_cnt_orgn = pad_ans;
                l_src_s.col_cnt = pad_ans; l_src_s.col_cnt_orgn = pad_ans;
                r_src_s.ln_cnt = pad_ans; r_src_s.ln_cnt_orgn = pad_ans;
                r_src_s.col_cnt = pad_ans; r_src_s.col_cnt_orgn = pad_ans;
                auto ans_pad = mtx_strassen_mult(l_src_s, r_src_s, recursive_gate);
                mtx_crop(ans, l_ln_cnt, r_col_cnt, ans_pad.ptr, pad_ans, 0, r_p_r, b_p_l);
                ans_pad.release();
            }
            l_src_s.release(); r_src_s.release();
            return true;
        }
        else return mtx_mult(ans, l_src, r_src, l_ln_cnt, l_col_cnt, r_ln_cnt, r_col_cnt);
    }
    else return false;
}

__mtx_callback_ inline void __mtx_t::para_set(uint64_t ln_cnt, uint64_t col_cnt) { mtx_ln_cnt = ln_cnt; mtx_col_cnt = col_cnt; mtx_elem_cnt = ln_cnt * col_cnt; }
__mtx_callback_ template<typename _Ity> inline void __mtx_t::__value_copy(const matrix<_Ity> &src)
{
    auto p_src = const_cast<matrix<_Ity>*>(&src);
    if(mtx_elem_cnt != p_src->get_elem_cnt())
    {
        MTX_RESET(mtx_ptr);
        MTX_ALLOC(__mtx_elem, mtx_ptr, p_src->get_elem_cnt());
    }
    para_set(p_src->get_ln_cnt(), p_src->get_col_cnt());
    for(auto i=0; i<mtx_elem_cnt; ++i) mtx_ptr[i] = p_src->pos_idx(i);
    p_src = nullptr;
}
__mtx_callback_ template<typename _Ity> inline void __mtx_t::value_move(matrix<_Ity> &&src)
{
    if constexpr (std::is_same_v<_Ity, __mtx_elem>)
    {
        mtx_move(mtx_ptr, std::move(src.mtx_ptr));
        para_set(src.mtx_ln_cnt, src.mtx_col_cnt);
    }
    else __value_copy(src);
    src.reset();
}
__mtx_callback_ inline void __mtx_t::init_list_mtx(matrix &&_vect, bagrt::net_sequence<matrix> &vect_set, bagrt::net_sequence<mtx_pos> &left_top_pos, uint64_t curr_left, uint64_t curr_top, uint64_t &left_increase, uint64_t &top_increase)
{
    left_increase = _vect.mtx_col_cnt;
    top_increase = _vect.mtx_ln_cnt;
    vect_set.emplace_back(std::move(_vect));
    left_top_pos.emplace_back(mtx_pos(curr_top, curr_left));
}
__mtx_callback_ template<typename _Ity> inline void __mtx_t::init_list_mtx(std::initializer_list<std::initializer_list<_Ity>> &&_vect, bagrt::net_sequence<matrix> &vect_set, bagrt::net_sequence<mtx_pos> &left_top_pos, uint64_t curr_left, uint64_t curr_top, uint64_t &left_increase, uint64_t &top_increase)
{
    auto curr_sub_top = curr_top;
    for(auto ln_temp : _vect)
    {
        auto curr_sub_left = curr_left;
        auto curr_top_increase = 0ull;
        for(auto temp : ln_temp)
        {
            auto curr_left_increase = 0ull;
            init_list_mtx(std::move(temp), vect_set, left_top_pos, curr_sub_left, curr_sub_top, curr_left_increase, curr_top_increase);
            curr_sub_left += curr_left_increase;
        }
        curr_sub_top += curr_top_increase;
    }
}
__mtx_callback_ inline __mtx_t::col_pointer::col_pointer(matrix *_mtx_buf_ptr, uint64_t begin_idx, uint64_t _ptr_len) : _buf_ptr(_mtx_buf_ptr), _ptr_ln(begin_idx), _ptr_col_cnt(_ptr_len) {}
__mtx_callback_ inline __mtx_elem &__mtx_t::col_pointer::operator[](uint64_t col)
{
    assert(col < _ptr_col_cnt);
    return *(_buf_ptr->mtx_ptr + _ptr_ln * _ptr_col_cnt + col);
}
__mtx_callback_ inline __mtx_t::col_pointer::~col_pointer() { _buf_ptr = nullptr; _ptr_ln = 0; _ptr_col_cnt = 0; }
__mtx_callback_ template<typename _Ity> inline void __mtx_t::value_copy(matrix<_Ity> &src) { __value_copy(src); }
__mtx_callback_ inline __mtx_t::matrix(matrix &src) { value_copy(src); }
__mtx_callback_ inline __mtx_t::matrix(const matrix &src) { __value_copy(src); }
__mtx_callback_ inline __mtx_t::matrix(matrix &&src) { value_move(std::move(src)); }
__mtx_callback_ inline __mtx_t::matrix(uint64_t ln_cnt, uint64_t col_cnt, bool rand_mode, __mtx_elem &&rand_boundary_first, __mtx_elem &&rand_boundary_second, __mtx_elem &&rand_acc)
{
    if(ln_cnt && col_cnt)
    {
        para_set(ln_cnt, col_cnt);
        MTX_ALLOC(__mtx_elem, mtx_ptr, mtx_elem_cnt);
        if(rand_mode) mtx_init_rand(mtx_ptr, mtx_elem_cnt, std::move(rand_boundary_first), std::move(rand_boundary_second), std::move(rand_acc));
        else mtx_fill(mtx_ptr, std::move(__mtx_elem(0)), mtx_elem_cnt);
    }
}
__mtx_callback_ inline __mtx_t::matrix(__mtx &&ptr_src, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(ptr_src && ln_cnt && col_cnt)
    {
        para_set(ln_cnt, col_cnt);
        _BAGRT ptr_move(mtx_ptr, std::move(ptr_src)); 
    }
}
__mtx_callback_ inline __mtx_t::matrix(__mtx &ptr_src, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(ptr_src && ln_cnt && col_cnt)
    {
        para_set(ln_cnt, col_cnt);
        MTX_ALLOC(__mtx_elem, mtx_ptr, mtx_elem_cnt);
        _BAGRT ptr_copy(mtx_ptr, ptr_src, mtx_elem_cnt);   
    }
}
__mtx_callback_ template<typename _Ity> inline __mtx_t::matrix(_Ity &&atom)
{
    MTX_ALLOC(__mtx_elem, mtx_ptr, 1); mtx_ptr[0] = atom;
    para_set(1, 1);
}
__mtx_callback_ inline __mtx_t::matrix(bagrt::net_sequence<matrix> &_vect, bagrt::net_sequence<mtx_pos> &left_top_pos)
{
    if(_vect.length == left_top_pos.length)
    {
        auto ln_cnt = 0, col_cnt = 0;
        for(auto i=0; i<left_top_pos.length; ++i)
        {
            if(!left_top_pos[i].ln) col_cnt += _vect[i].mtx_col_cnt;
            if(!left_top_pos[i].col) ln_cnt += _vect[i].mtx_ln_cnt;
        }
        *this = matrix(ln_cnt, col_cnt);
        for(auto i=0; i<_vect.length; ++i) for(auto j=0; j<_vect[i].mtx_ln_cnt; ++j) for(auto k=0; k<_vect[i].mtx_col_cnt; ++k) (*this)[left_top_pos[i].ln+j][left_top_pos[i].col+k] = _vect[i][j][k];
    }
}
__mtx_callback_ template<typename _Ity> inline __mtx_t::matrix(std::initializer_list<std::initializer_list<_Ity>> _vect)
{
    if(_vect.size() && _vect.begin()->size())
    {
        bagrt::net_sequence<matrix> temp_seq;
        bagrt::net_sequence<mtx_pos> left_top_pos;
        auto top_temp = 0ull, left_temp = 0ull;
        init_list_mtx(std::move(_vect), temp_seq, left_top_pos, 0, 0, left_temp, top_temp);
        *this = matrix(temp_seq, left_top_pos);
    }
}
__mtx_callback_ inline uint64_t __mtx_t::get_ln_cnt() { return mtx_ln_cnt; }
__mtx_callback_ inline uint64_t __mtx_t::get_col_cnt() { return mtx_col_cnt; }
__mtx_callback_ inline uint64_t __mtx_t::get_elem_cnt() { return mtx_elem_cnt; }
__mtx_callback_ inline __mtx __mtx_t::get_mtx_ptr()
{
    MTX_INIT(__mtx_elem, ans, mtx_elem_cnt);
    if(! _BAGRT ptr_copy(ans, mtx_ptr, mtx_elem_cnt)) { MTX_RESET(ans); }
    return ans;
}
__mtx_callback_ inline bool __mtx_t::is_matrix() { return (mtx_ptr && mtx_ln_cnt && mtx_elem_cnt && mtx_col_cnt && (mtx_ln_cnt*mtx_col_cnt==mtx_elem_cnt)); }
__mtx_callback_ inline bool __mtx_t::fill_with(__mtx_elem &&val) { return mtx_fill(mtx_ptr, std::move(val), mtx_elem_cnt); }
__mtx_callback_ inline __mtx_elem __mtx_t::get_determinant() { assert(mtx_ln_cnt == mtx_col_cnt) return mtx_det(mtx_ptr, mtx_ln_cnt);}
__mtx_callback_ inline __mtx_t __mtx_t::get_inverser()
{
    assert(mtx_ln_cnt == mtx_col_cnt);
    auto det_ans = mtx_det(mtx_ptr, mtx_ln_cnt);
    matrix ans;
    if(det_ans != 0)
    {
        ans = matrix(mtx_ln_cnt, mtx_col_cnt);
        if(!mtx_inverser(ans.mtx_ptr, mtx_ptr, mtx_ln_cnt)) ans.reset();
    }
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::get_transposition()
{
    matrix ans(mtx_col_cnt, mtx_ln_cnt);
    if(!mtx_transposition(ans.mtx_ptr, mtx_ptr, mtx_ln_cnt, mtx_col_cnt)) ans.reset();
    return ans;
}
__mtx_callback_ inline __mtx_callback_t mtx_extm_data __mtx_t::extremum(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, bool max_flag, bool get_extm_pos, uint64_t ln_dilation, uint64_t col_dilation)
{
    mtx_extm_data ans;
    uint64_t *ln_ls = nullptr, *col_ls = nullptr, ls_cnt = 0;
    ans.ans = mtx_extm_val(mtx_ptr, from_ln, to_ln, from_col, to_col, mtx_ln_cnt, mtx_col_cnt, ln_ls, col_ls, ls_cnt, max_flag, get_extm_pos, ln_dilation, col_dilation);
    if(ls_cnt)
    {
        ans.ans_pos.init(ls_cnt);
        for(auto i=0; i<ls_cnt; ++i)
        {
            ans.ans_pos[i].ln = ln_ls[i];
            ans.ans_pos[i].col = col_ls[i];
        }
        MEM_RECYCLE(ln_ls); MEM_RECYCLE(col_ls);
    }
    return ans;
}
__mtx_callback_ inline __mtx_elem __mtx_t::get_atom() { assert(mtx_elem_cnt == 1); return mtx_ptr[0]; }
__mtx_callback_ inline __mtx_t __mtx_t::child(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilation, uint64_t col_dilation)
{
    auto child_ln_cnt = mtx_child_vec_dir_cnt(from_ln, to_ln, ln_dilation),
        child_col_cnt = mtx_child_vec_dir_cnt(from_col, to_col, col_dilation);
    matrix ans(child_ln_cnt, child_col_cnt);
    if(!mtx_child_vec(ans.mtx_ptr, child_ln_cnt, child_col_cnt, mtx_ptr, from_ln, to_ln, from_col, to_col, mtx_ln_cnt, mtx_col_cnt, ln_dilation, col_dilation)) ans.reset();
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::rotate_rect(bool clockwise)
{
    matrix ans(mtx_col_cnt, mtx_ln_cnt);
    if(!mtx_rotate_rect(ans.mtx_ptr, mtx_ptr, mtx_ln_cnt, mtx_col_cnt, clockwise)) ans.reset();
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::mirror_flip(bool symmetry_vertical)
{
    matrix ans(mtx_ln_cnt, mtx_col_cnt);
    if(!mtx_mirror_flip(ans.mtx_ptr, mtx_ptr, mtx_ln_cnt, mtx_col_cnt, symmetry_vertical)) ans.reset();
    return ans;
}
__mtx_callback_ inline bool __mtx_t::shape_valid(uint64_t ln_cnt, uint64_t col_cnt) { return (ln_cnt==mtx_ln_cnt && col_cnt==mtx_col_cnt); }
__mtx_callback_ inline bool __mtx_t::shape_valid(matrix &src) { return shape_valid(src.mtx_ln_cnt, src.mtx_col_cnt); }
__mtx_callback_ inline __mtx_t __mtx_t::reshape(uint64_t ln_cnt, uint64_t col_cnt)
{
    assert(ln_cnt*col_cnt == mtx_elem_cnt);
    auto ans = *this;
    ans.mtx_ln_cnt = ln_cnt;
    ans.mtx_col_cnt = col_cnt;
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::reshape(matrix &as_val) { assert(mtx_elem_cnt == as_val.mtx_elem_cnt); return reshape(as_val.mtx_ln_cnt, as_val.mtx_col_cnt); }
__mtx_callback_ inline __mtx_elem __mtx_t::elem_sum(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilation, uint64_t col_dilation) { return mtx_sum(mtx_ptr, from_ln, to_ln, from_col, to_col, mtx_ln_cnt, mtx_col_cnt, ln_dilation, col_dilation); }
__mtx_callback_ inline __mtx_elem __mtx_t::elem_sum() { return elem_sum(0, mtx_ln_cnt-1, 0, mtx_col_cnt-1); }
__mtx_callback_ inline __mtx_t __mtx_t::abs()
{
    auto ans = *this;
    if(!mtx_abs(ans.mtx_ptr, ans.mtx_elem_cnt)) ans.reset();
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::elem_cal_opt(matrix &r_val, uint64_t opt_idx)
{
    assert(shape_valid(r_val));
    matrix ans = *this;
    if(!mtx_elem_cal_opt(ans.mtx_ptr, mtx_ptr, r_val.mtx_ptr, mtx_ln_cnt, mtx_col_cnt, opt_idx)) ans.reset();
    return ans;        
}
__mtx_callback_ inline __mtx_t __mtx_t::elem_cal_opt(__mtx_elem &&para, uint64_t opt_idx)
{
    auto ans = *this;
    if(!mtx_elem_cal_opt(ans.mtx_ptr, mtx_ptr, std::move(para), mtx_ln_cnt, mtx_col_cnt, opt_idx)) ans.reset();
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::broadcast_add(__mtx_elem &&val, bool nega)
{
    auto ans = *this;
    if(!mtx_broadcast_add(ans.mtx_ptr, mtx_ptr, mtx_elem_cnt, std::move(val), nega)) ans.reset();
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::pad(uint64_t ln_t, uint64_t col_r, uint64_t ln_b, uint64_t col_l, uint64_t ln_dist, uint64_t col_dist)
{
    if(ln_t || col_r || ln_b || col_l || ln_dist || col_dist)
    {
        auto ln_cnt_p = mtx_pad_cnt(ln_t, ln_b, mtx_ln_cnt, ln_dist),
        col_cnt_p = mtx_pad_cnt(col_l, col_r, mtx_col_cnt, col_dist);
        matrix ans(ln_cnt_p, col_cnt_p);
        if(!mtx_pad(ans.mtx_ptr, col_cnt_p, mtx_ptr, mtx_ln_cnt, mtx_col_cnt, ln_t, col_r, ln_b, col_l, ln_dist, col_dist)) ans.reset();
        return ans;
    }
    else return *this;
}
__mtx_callback_ inline __mtx_t __mtx_t::crop(uint64_t ln_t, uint64_t col_r, uint64_t ln_b, uint64_t col_l, uint64_t ln_dist, uint64_t col_dist)
{
    if(ln_t || col_r || ln_b || col_l || ln_dist || col_dist)
    {
        auto ln_cnt_c = mtx_crop_cnt(ln_t, ln_b, mtx_ln_cnt, ln_dist),
            col_cnt_c = mtx_crop_cnt(col_l, col_r, mtx_col_cnt, col_dist);
        matrix ans(ln_cnt_c, col_cnt_c);
        if(!mtx_crop(ans.mtx_ptr, ln_cnt_c, col_cnt_c, mtx_ptr, mtx_col_cnt, ln_t, col_r, ln_b, col_l, ln_dist, col_dist)) ans.reset();
        return ans;
    }
    else return *this;
}
__mtx_callback_ inline __mtx_t __mtx_t::round_fit()
{
    matrix ans;
    if(is_matrix())
    {
        ans = *this;
        for(auto i=0; i<mtx_elem_cnt; ++i) ans.mtx_ptr[i] = _BAGRT __to_integer_(ans.mtx_ptr[i]);
    }
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::LU_decompose()
{
    auto ans = *this;
    if(!mtx_LU(ans.mtx_ptr, mtx_ptr, mtx_ln_cnt)) ans.reset();
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::swap_dir_elem(uint64_t l_idx, uint64_t r_idx, bool is_ln)
{
    MTX_INIT(__mtx_elem, ans_ptr, mtx_elem_cnt);
    if(!mtx_swap_elem(ans_ptr, mtx_ptr, l_idx, r_idx, mtx_ln_cnt, mtx_col_cnt, is_ln)) { MTX_RESET(ans_ptr); }
    return matrix(std::move(ans_ptr), mtx_ln_cnt, mtx_col_cnt);
}
__mtx_callback_ inline __mtx_t __mtx_t::get_adjugate()
{
    auto ans = *this;
    if(!mtx_adjugate(ans.mtx_ptr, mtx_ptr, mtx_ln_cnt)) ans.reset();
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::mult_strassen(matrix &r_src, uint64_t recursive_gate)
{
    matrix ans(mtx_ln_cnt, r_src.mtx_col_cnt);
    if(!mtx_strassen_mult(ans.mtx_ptr, mtx_ptr, r_src.mtx_ptr, mtx_ln_cnt, mtx_col_cnt, r_src.mtx_ln_cnt, r_src.mtx_col_cnt, recursive_gate)) ans.reset();
    return ans;
}
__mtx_callback_ inline uint64_t __mtx_t::get_rank() { return mtx_rank(mtx_ptr, mtx_ln_cnt, mtx_col_cnt); }
__mtx_callback_ inline matrix<_BAGRT decimal> __mtx_t::get_mtx_dec()
{
    static_assert(!std::is_same<__mtx_elem, _BAGRT decimal>::value, "Type \"bagrt::decimal\" could not be transformed.");
    matrix<_BAGRT decimal> ans(mtx_ln_cnt, mtx_col_cnt);
    for(auto i=0; i<mtx_elem_cnt; ++i) ans.pos_idx(i) = mtx_ptr[i];
    return ans;
}
__mtx_callback_ inline matrix<long double> __mtx_t::get_mtx_f()
{
    static_assert(!std::is_floating_point<__mtx_elem>::value, "Float point type could not be transformed.");
    matrix<long double> ans(mtx_ln_cnt, mtx_col_cnt);
    for(auto i=0; i<mtx_elem_cnt; ++i) ans.pos_idx(i) = __to_float_point_(mtx_ptr[i]);
    return ans;
}
__mtx_callback_ inline __mtx_elem &__mtx_t::pos_idx(uint64_t idx) { assert(idx < mtx_elem_cnt); return mtx_ptr[idx]; }
__mtx_callback_ inline __mtx_t __mtx_t::operator+(matrix &val)
{
    matrix ans;
    if(shape_valid(val))
    {
        ans = *this;
        if(!mtx_add(ans.mtx_ptr, mtx_ptr, val.mtx_ptr, mtx_ln_cnt, mtx_col_cnt)) ans.reset();
    }
    return ans;
}
__mtx_callback_ inline __mtx_t __mtx_t::operator-(matrix &val)
{
    matrix ans;
    if(shape_valid(val))
    {
        ans = *this;
        if(!mtx_add(ans.mtx_ptr, mtx_ptr, val.mtx_ptr, mtx_ln_cnt, mtx_col_cnt, true)) ans.reset();
    }
    return ans;
}
__mtx_callback_ inline void __mtx_t::operator+=(matrix &val) { *this = *this + val; }
__mtx_callback_ inline void __mtx_t::operator-=(matrix &val) { *this = *this - val; }
__mtx_callback_ inline __mtx_t __mtx_t::operator*(matrix &val)
{
    matrix ans(mtx_ln_cnt, val.mtx_col_cnt);
    if(!mtx_mult(ans.mtx_ptr, mtx_ptr, val.mtx_ptr, mtx_ln_cnt, mtx_col_cnt, val.mtx_ln_cnt, val.mtx_col_cnt)) ans.reset();
    return ans;
}
__mtx_callback_ inline void __mtx_t::operator*=(matrix &val) { *this = *this * val; }
__mtx_callback_ inline __mtx_t __mtx_t::operator*(__mtx_elem &&val)
{
    auto ans = *this;
    if(!mtx_mult(ans.mtx_ptr, mtx_ptr, std::move(val), mtx_ln_cnt, mtx_col_cnt)) ans.reset();
    return ans;
}
__mtx_callback_ inline void __mtx_t::operator*=(__mtx_elem &&val) { *this = *this * std::move(val); }
__mtx_callback_ template<typename _Ity> inline void __mtx_t::operator=(const matrix<_Ity> &src) { __value_copy(src); }
__mtx_callback_ template<typename _Ity> inline void __mtx_t::operator=(matrix<_Ity> &src) { value_copy(src); }
__mtx_callback_ template<typename _Ity> inline void __mtx_t::operator=(matrix<_Ity> &&src) { value_move(std::move(src)); }
__mtx_callback_ inline bool __mtx_t::operator==(matrix &val)
{
    if(shape_valid(val)) return mtx_equal(val.mtx_ptr, mtx_ptr, mtx_elem_cnt);
    else return false;
}
__mtx_callback_ inline bool __mtx_t::operator!=(matrix &val) { return !(*this == val); }
__mtx_callback_ inline typename __mtx_t::col_pointer __mtx_t::operator[](uint64_t ln) { assert(ln < mtx_ln_cnt); return col_pointer(this, ln, mtx_col_cnt); }
__mtx_callback_ inline __mtx_t __mtx_t::init_E(uint64_t dms)
{
    matrix ans(dms, dms);
    if(!mtx_init_E(ans.mtx_ptr, dms)) ans.reset();
    return ans;
}
__mtx_callback_ inline void __mtx_t::reset() { mtx_ln_cnt = 0; mtx_col_cnt = 0; mtx_elem_cnt = 0; MTX_RESET(mtx_ptr); }
__mtx_callback_ inline __mtx_t::~matrix() { reset();}

MTX_END