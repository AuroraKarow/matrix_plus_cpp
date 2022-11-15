MATRIX_BEGIN

matrix_declare struct net_strassen_matrix final {
    uint64_t ln_cnt       = 0,
             ln_cnt_orgn  = 0,
             col_cnt      = 0,
             col_cnt_orgn = 0,
             ln_begin     = 0,
             col_begin    = 0;
               
    matrix_ptr ptr = nullptr;

    uint64_t curr_ln(uint64_t ln) const { return (ln + ln_begin) * col_cnt_orgn + col_begin; }

    void init_para(uint64_t _ln_cnt, uint64_t _col_cnt) {
        ln_cnt       = _ln_cnt;
        col_cnt      = _col_cnt;
        ln_cnt_orgn  = _ln_cnt;
        col_cnt_orgn = _col_cnt;
    }

    void allocate(uint64_t _ln_cnt, uint64_t _col_cnt) {
        if (ln_cnt == ln_cnt_orgn && col_cnt == col_cnt_orgn) reset(ptr);
        ptr = init<matrix_elem_t>(_ln_cnt * _col_cnt);
        init_para(_ln_cnt, _col_cnt);
    }

    long double *operator[](uint64_t ln) const { return ptr + curr_ln(ln); }

    net_strassen_matrix operator+(const net_strassen_matrix &right) const {
        net_strassen_matrix ans;
        if (!(ln_cnt == right.ln_cnt && col_cnt == right.col_cnt)) return ans;
        ans.allocate(ln_cnt, col_cnt);
        for (auto i = 0ull; i < ln_cnt; ++i) for (auto j = 0ull; j < col_cnt; ++j) ans[i][j] = ptr[curr_ln(i) + j] + right.ptr[right.curr_ln(i) + j];
        return ans;
    }

    net_strassen_matrix operator-(const net_strassen_matrix &right) const {
        net_strassen_matrix ans;
        if (!(ln_cnt == right.ln_cnt && col_cnt == right.col_cnt)) return ans;
        ans.allocate(ln_cnt, col_cnt);
        for (auto i = 0ull; i < ln_cnt; ++i) for (auto j = 0ull; j < col_cnt; ++j) ans[i][j] = ptr[curr_ln(i) + j] - right.ptr[right.curr_ln(i) + j];
        return ans;
    }
    
    net_strassen_matrix operator*(const net_strassen_matrix &right) const {
        net_strassen_matrix ans;
        if (col_cnt != right.ln_cnt) return ans;
        ans.allocate(ln_cnt, right.col_cnt);
        for (auto i = 0ull; i < ln_cnt; ++i) for (auto j = 0ull; j < col_cnt; ++j) {
            auto coe = ptr[curr_ln(i) + j];
            for(auto k = 0ull; k < ans.col_cnt; ++k) ans[i][k] += coe * right.ptr[right.curr_ln(j) + k];
        }
        return ans;
    }

    void release() {
        if (ln_cnt == ln_cnt_orgn && col_cnt == col_cnt_orgn && ptr) reset(ptr);
        else ptr = nullptr;
    }

    void print() {
        for (auto i = 0ull; i < ln_cnt; ++i) {
            for (auto j = 0ull; j < col_cnt; ++j) std::cout << *(ptr + curr_ln(i) + j) << '\t';
            std::cout << '\n';
        }
    }
};

callback_matrix void strassen_child(net_strassen_matrix<matrix_elem_t> &ans, const net_strassen_matrix<matrix_elem_t> &src, uint64_t ln_from, uint64_t col_from, uint64_t ln_cnt_child, uint64_t col_cnt_child) {
    ans.ptr          = src.ptr; 
    ans.ln_cnt_orgn  = src.ln_cnt_orgn; 
    ans.col_cnt_orgn = src.col_cnt_orgn;
    ans.ln_cnt       = ln_cnt_child;
    ans.col_cnt      = col_cnt_child;
    ans.ln_begin     = src.ln_begin + ln_from;
    ans.col_begin    = src.col_begin + col_from;
}

callback_matrix void strassen_quartile(const net_strassen_matrix<matrix_elem_t> &src, net_strassen_matrix<matrix_elem_t> &src00, net_strassen_matrix<matrix_elem_t> &src01, net_strassen_matrix<matrix_elem_t> &src10, net_strassen_matrix<matrix_elem_t> &src11) {
    auto dim_cnt_sub = src.ln_cnt / 2;
    strassen_child(src00, src, 0, 0, dim_cnt_sub, dim_cnt_sub);
    strassen_child(src01, src, 0, dim_cnt_sub, dim_cnt_sub, dim_cnt_sub);
    strassen_child(src10, src, dim_cnt_sub, 0, dim_cnt_sub, dim_cnt_sub);
    strassen_child(src11, src, dim_cnt_sub, dim_cnt_sub, dim_cnt_sub, dim_cnt_sub);
}

callback_matrix net_strassen_matrix<matrix_elem_t> strassen_mult(const net_strassen_matrix<matrix_elem_t> &left, const net_strassen_matrix<matrix_elem_t> &right, uint64_t recursive_gate) {
    if (left.ln_cnt <= recursive_gate) return left * right;

    net_strassen_matrix<matrix_elem_t> a00, a01, a10, a11, b00, b01, b10, b11;
    
    strassen_quartile(left, a00, a01, a10, a11);
    strassen_quartile(right, b00, b01, b10, b11);

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

    auto p0 = strassen_mult(a00, s0, recursive_gate);
    auto p1 = strassen_mult(s1, b11, recursive_gate);
    auto p2 = strassen_mult(s2, b00, recursive_gate);
    auto p3 = strassen_mult(a11, s3, recursive_gate);
    auto p4 = strassen_mult(s4, s5, recursive_gate);
    auto p5 = strassen_mult(s6, s7, recursive_gate);
    auto p6 = strassen_mult(s8, s9, recursive_gate);

    net_strassen_matrix<matrix_elem_t> ans;
    ans.allocate(left.ln_cnt, right.col_cnt);
    auto dim_cnt_sub = ans.ln_cnt / 2;
    for (auto i = 0ull; i < ans.ln_cnt; ++i) for (auto j = 0ull; j < ans.col_cnt; ++j) if (i < dim_cnt_sub) {
        if (j < dim_cnt_sub) ans[i][j] = p4[i][j] + p3[i][j] - p1[i][j] + p5[i][j];
        else ans[i][j] = p0[i][j - dim_cnt_sub] + p1[i][j - dim_cnt_sub];
    } else {
        if (j < dim_cnt_sub) ans[i][j] = p2[i - dim_cnt_sub][j] + p3[i - dim_cnt_sub][j];
        else ans[i][j] = p4[i - dim_cnt_sub][j - dim_cnt_sub] + p0[i - dim_cnt_sub][j - dim_cnt_sub] - p2[i - dim_cnt_sub][j - dim_cnt_sub] - p6[i - dim_cnt_sub][j - dim_cnt_sub];
    }

    s0.release(); s1.release(); s2.release(); s3.release(); s4.release(); s5.release(); s6.release(); s7.release(); s8.release(); s9.release();

    p0.release(); p1.release(); p2.release(); p3.release(); p4.release(); p5.release(); p6.release();

    return ans;
}

callback_matrix matrix_ptr strassen_mult(const matrix_ptr left, uint64_t left_ln_cnt, uint64_t left_col_cnt, const matrix_ptr right, uint64_t right_col_cnt, uint64_t recursive_gate = 32) {
    auto r_p_l  = left_col_cnt,
         b_p_l  = left_ln_cnt,
         r_p_r  = right_col_cnt,
         b_p_r  = left_col_cnt,
         mx_len = num_extreme({r_p_l, b_p_l, r_p_r, b_p_r});
    if (mx_len <= recursive_gate) return mult(left, left_ln_cnt, left_col_cnt, right, right_col_cnt);
    auto pad_val = num_pad_pow(mx_len, 2),
         pad_ans = mx_len + pad_val;
         r_p_l   = pad_ans - r_p_l;
         b_p_l   = pad_ans - b_p_l;
         r_p_r   = pad_ans - r_p_r;
         b_p_r   = pad_ans - b_p_r;
    net_strassen_matrix<matrix_elem_t> left_stsn, right_stsn;
    left_stsn.init_para(pad_ans, pad_ans);
    right_stsn.init_para(pad_ans, pad_ans);
    left_stsn.ptr = pad(pad_ans, pad_ans, left, left_ln_cnt, left_col_cnt, 0, r_p_l, b_p_l),
    right_stsn.ptr = pad(pad_ans, pad_ans, right, left_col_cnt, right_col_cnt, 0, r_p_r, b_p_r);
    auto ans = strassen_mult(left_stsn, right_stsn, recursive_gate);
    auto ans_ptr = crop(left_ln_cnt, right_col_cnt, ans.ptr, pad_ans, pad_ans, 0, r_p_r, b_p_l);
    left_stsn.release();
    right_stsn.release();
    ans.release();
    return ans_ptr;
}

matrix_declare class net_matrix {
protected:
    void value_assign(const net_matrix &src) {
        ln_cnt   = src.ln_cnt;
        col_cnt  = src.col_cnt;
        elem_cnt = src.elem_cnt;
    }
    
    void value_copy(const net_matrix &src) {
        copy(ptr, elem_cnt, src.ptr, src.elem_cnt);
        value_assign(src);
    }
    template<typename s_arg> void value_copy(const net_matrix<s_arg> &src) {
        ptr_alter(ptr, elem_cnt, src.__elem_cnt__(), false);
        ln_cnt   = src.__ln_cnt__();
        col_cnt  = src.__col_cnt__();
        elem_cnt = src.__elem_cnt__();
        for (auto i = 0ull; i < elem_cnt; ++i) if constexpr (std::is_same_v<s_arg, net_decimal>) *(ptr + i) = src.index(i).number_format;
        else *(ptr + i) = src.index(i);
    }

    void value_move(net_matrix &&src) {
        move(ptr, std::move(src.ptr));
        value_assign(src);
        src.reset();
    }

    void init_list_mtx(net_matrix &&_vect, net_sequence<net_matrix> &vect_set, net_sequence<pos> &left_top_pos, uint64_t curr_left, uint64_t curr_top, uint64_t &left_increase, uint64_t &top_increase) const {
        left_increase = _vect.col_cnt;
        top_increase  = _vect.ln_cnt;
        vect_set.emplace_back(std::move(_vect));
        left_top_pos.emplace_back(pos{curr_top, curr_left});
    }

    void init_list_mtx(std::initializer_list<std::initializer_list<net_matrix>> _vect, net_sequence<net_matrix> &vect_set, net_sequence<pos> &left_top_pos, uint64_t curr_left, uint64_t curr_top, uint64_t &left_increase, uint64_t &top_increase) const {
        auto curr_sub_top = curr_top;
        for(auto ln_temp : _vect) {
            auto curr_sub_left     = curr_left,
                 curr_top_increase = 0ull;
            for(auto temp : ln_temp) {
                auto curr_left_increase = 0ull;
                init_list_mtx(std::move(temp), vect_set, left_top_pos, curr_sub_left, curr_sub_top, curr_left_increase, curr_top_increase);
                curr_sub_left += curr_left_increase;
            }
            curr_sub_top += curr_top_increase;
        }
    }

    struct line_data final {
    public:
        line_data(const net_matrix *buf_ptr = nullptr, uint64_t begin_idx = 0, uint64_t ptr_len = 0) :
            _buf_ptr_(buf_ptr),
            _ptr_ln_(begin_idx),
            _ptr_col_cnt(ptr_len) {}
        
        matrix_elem_t &operator[](uint64_t col) const {
            assert(col < _ptr_col_cnt);
            return *(_buf_ptr_->ptr + _ptr_ln_ * _ptr_col_cnt + col);
        }
        
        ~line_data() {
            _buf_ptr_    = nullptr;
            _ptr_ln_     = 0;
            _ptr_col_cnt = 0;
        }
    private:
        const net_matrix *_buf_ptr_ = nullptr;

        uint64_t _ptr_ln_     = 0,
                 _ptr_col_cnt = 0;
    };

public:
    net_matrix(uint64_t mtx_ln_cnt, uint64_t mtx_col_cnt, bool rand_init = false, const matrix_elem_t &fst_rng = 0, const matrix_elem_t &snd_rng = 0, uint64_t acc = 5) :
        ln_cnt(mtx_ln_cnt),
        col_cnt(mtx_col_cnt),
        elem_cnt(mtx_ln_cnt * mtx_col_cnt) {
            if (elem_cnt == 0) return;
            if (rand_init) ptr = init_rand<matrix_elem_t>(elem_cnt, fst_rng, snd_rng, acc);
            else ptr = init<matrix_elem_t>(elem_cnt);
        }
    template<typename matrix_elem_para, typename matrix_elem_para_v> net_matrix(const i_arg &src) :
        elem_cnt(1),
        ln_cnt(1),
        col_cnt(1),
        ptr(init<matrix_elem_t>(1)) {
        if constexpr (std::is_same_v<matrix_elem_para, net_decimal> && !std::is_same_v<matrix_elem_t, matrix_elem_para>) *ptr = src.number_format;
        else *ptr = src;
    }
    net_matrix(matrix_ptr &&src = nullptr, uint64_t mtx_ln_cnt = 0, uint64_t mtx_col_cnt = 0) :
        ln_cnt(mtx_ln_cnt),
        col_cnt(mtx_col_cnt),
        elem_cnt(mtx_ln_cnt * mtx_col_cnt) { ptr_move(ptr, std::move(src)); }
    net_matrix(const net_sequence<net_matrix> &_vect, const net_sequence<pos> &left_top_pos) {
        if (_vect.length != left_top_pos.length) return;
        uint64_t ln_cnt_temp  = 0,
                 col_cnt_temp = 0;
        for (auto i = 0ull; i < left_top_pos.length; ++i) {
            if(!left_top_pos[i].ln) col_cnt_temp += _vect[i].col_cnt;
            if(!left_top_pos[i].col) ln_cnt_temp += _vect[i].ln_cnt;
        }
        value_move(net_matrix(ln_cnt_temp, col_cnt_temp));
        for (auto i = 0ull; i < _vect.length; ++i) for (auto j = 0ull; j < _vect[i].ln_cnt; ++j) for(auto k = 0ull; k < _vect[i].col_cnt; ++k) (*this)[left_top_pos[i].ln + j][left_top_pos[i].col + k] = std::move(_vect[i][j][k]);
    }
    net_matrix(std::initializer_list<std::initializer_list<net_matrix>> _vect) {
        if(_vect.size() && _vect.begin()->size()) {
            net_sequence<net_matrix> temp_seq;
            net_sequence<pos> left_top_pos;
            uint64_t top_temp  = 0,
                     left_temp = 0;
            init_list_mtx(_vect, temp_seq, left_top_pos, 0, 0, left_temp, top_temp);
            value_move(net_matrix(temp_seq, left_top_pos));
        }
    }
    net_matrix(const net_matrix &src) { value_copy(src); }
    net_matrix(net_matrix &&src) { value_move(std::move(src)); }
    template<typename s_arg> net_matrix(const net_matrix<s_arg> &src) { value_copy(src); }

    bool is_matrix() const { return ptr && (ln_cnt * col_cnt == elem_cnt) && ln_cnt && col_cnt; }

    template<typename matrix_elem_para, typename matrix_elem_para_v> void fill_elem(const matrix_elem_para &src) {
        if constexpr (!std::is_same_v<matrix_elem_para, matrix_elem_t> && std::is_same_v<matrix_elem_para, net_decimal>) fill(ptr, elem_cnt, src.number_format);
        else fill<matrix_elem_t>(ptr, elem_cnt, src);
    }

    net_list<pos> extremum_position(bool get_max, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilate = 0, uint64_t col_dilate = 0) const {
        net_list<pos> ans;
        extremum(ans, get_max, ptr, ln_cnt, col_cnt, from_ln, to_ln, from_col, to_col, ln_dilate, col_dilate);
        return ans;
    }
    net_list<pos> extremum_position(bool get_max = true) const { return extremum_position(get_max, 0, ln_cnt - 1, 0, col_cnt - 1); }

    net_matrix child(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilate, uint64_t col_dilate) const {
        uint64_t ln_cnt_temp  = 0,
                 col_cnt_temp = 0;
        auto     ans_temp     = sub(ln_cnt_temp, col_cnt_temp, ptr, ln_cnt, col_cnt, from_ln, to_ln, from_col, to_col, ln_dilate, col_dilate);
        return net_matrix(std::move(ans_temp), ln_cnt_temp, col_cnt_temp);
    }
    net_matrix child(uint64_t from_ln, uint64_t c_ln_cnt, uint64_t from_col, uint64_t c_col_cnt) { return child(from_ln, from_ln + c_ln_cnt - 1, from_col, from_col + c_col_cnt - 1, 0, 0); }

    net_matrix rotate_rectangle(bool clockwise = true) const { return net_matrix(rotate_rect(ptr, ln_cnt, col_cnt, clockwise), col_cnt, ln_cnt); }

    net_matrix mirror_flipping(bool vertical_symmetry = true) const { return net_matrix(mirror_flip(ptr, ln_cnt, col_cnt, vertical_symmetry), ln_cnt, col_cnt); }

    bool shape_verify(uint64_t src_ln_cnt, uint64_t src_col_cnt) const { return ln_cnt == src_ln_cnt && col_cnt == src_col_cnt; }
    bool shape_verify(const net_matrix &src) const { return shape_verify(src.ln_cnt, src.col_cnt); }

    net_matrix reshape(uint64_t re_ln_cnt, uint64_t re_col_cnt) const {
        auto ans = *this;
        if (re_ln_cnt * re_col_cnt == elem_cnt) {
            ans.ln_cnt  = re_ln_cnt;
            ans.col_cnt = re_col_cnt;
        }
        return ans;
    }
    net_matrix reshape(const net_matrix &src) const { return reshape(src.ln_cnt, src.col_cnt); }

    net_matrix elem_swap(uint64_t fst_idx, uint64_t snd_idx, bool ln_swap = true) { return net_matrix(ln_col_swap(ptr, ln_cnt, col_cnt, fst_idx, snd_idx, ln_swap), ln_cnt, col_cnt); }

    matrix_elem_t elem_sum() {
        matrix_elem_t ans = 0;
        for (auto i = 0ull; i < elem_cnt; ++i) ans += (*(ptr + i));
        return ans;
    }
    matrix_elem_t elem_sum(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilate, uint64_t col_dilate) const { return child(from_ln, to_ln, from_col, to_col, ln_dilate, col_dilate).elem_sum(); }

    net_matrix elem_wise_opt(const net_matrix &val, uint64_t operation, bool elem_swap = false, long double epsilon = 1e-8) const {
        if (shape_verify(val)) return net_matrix(elem_operate(ptr, val.ptr, elem_cnt, operation, elem_swap, epsilon), ln_cnt, col_cnt);
        else return net_matrix();
    }
    template<typename matrix_elem_para, typename matrix_elem_para_v> net_matrix elem_wise_opt(const matrix_elem_para &para, uint64_t operation, bool para_fst = false, long double epsilon = 1e-8) const {
        if constexpr (std::is_same_v<matrix_elem_para, net_decimal> && !std::is_same_v<matrix_elem_t, net_decimal>) return net_matrix(elem_operate(ptr, elem_cnt, para.number_format, operation, para_fst, epsilon), ln_cnt, col_cnt);
        else return net_matrix(elem_operate<matrix_elem_t>(ptr, elem_cnt, para, operation, para_fst, epsilon), ln_cnt, col_cnt);
    }

    template<typename matrix_elem_para, typename matrix_elem_para_v> net_matrix broadcast_addition(const matrix_elem_para &para, bool subtract = false, bool para_fst = false) const {
        if constexpr (std::is_same_v<matrix_elem_para, net_decimal> && !std::is_same_v<matrix_elem_t, net_decimal>) return net_matrix(broadcast_add(ptr, elem_cnt, para.number_format, subtract, para_fst), ln_cnt, col_cnt);
        else return net_matrix(broadcast_add<matrix_elem_t>(ptr, elem_cnt, para, subtract, para_fst), ln_cnt, col_cnt);
    }

    net_matrix padding(uint64_t top = 0, uint64_t right = 0, uint64_t bottom = 0, uint64_t left = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0) const {
        uint64_t ln_cnt_temp  = 0,
                 col_cnt_temp = 0;
        auto     ans_temp     = pad(ln_cnt_temp, col_cnt_temp, ptr, ln_cnt, col_cnt, top, right, bottom, left, ln_dist, col_dist);
        return net_matrix(std::move(ans_temp), ln_cnt_temp, col_cnt_temp);
    }

    net_matrix cropping(uint64_t top = 0, uint64_t right = 0, uint64_t bottom = 0, uint64_t left = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0) {
                uint64_t ln_cnt_temp  = 0,
                 col_cnt_temp = 0;
        auto     ans_temp     = crop(ln_cnt_temp, col_cnt_temp, ptr, ln_cnt, col_cnt, top, right, bottom, left, ln_dist, col_dist);
        return net_matrix(std::move(ans_temp), ln_cnt_temp, col_cnt_temp);
    }

    matrix_elem_t &index(uint64_t idx) const {
        assert(idx < elem_cnt);
        return *(ptr + idx);
    }

    net_matrix<long double> float_point_format() const {
        if constexpr (std::is_same_v<net_decimal, matrix_elem_t>) {
            net_matrix<long double> ans(ln_cnt, col_cnt);
            for (auto i = 0ull; i < elem_cnt; ++i) ans.index(i) = (*(ptr + i)).number_format;
            return ans;
        } else return net_matrix(*this);
    }

    void reset() {
        recycle<matrix_elem_t>(ptr);
        ln_cnt   = 0;
        col_cnt  = 0;
        elem_cnt = 0;
    }

    ~net_matrix() { reset(); }

protected:
    uint64_t ln_cnt   = 0,
             col_cnt  = 0,
             elem_cnt = 0;
    
    matrix_ptr ptr = nullptr;

public:

    uint64_t __ln_cnt__() const { return ln_cnt; }

    uint64_t __col_cnt__() const { return col_cnt; }

    uint64_t __elem_cnt__() const { return elem_cnt; }
    
    matrix_elem_t __determinant__() const {
        assert(ln_cnt == col_cnt);
        return det(ptr, ln_cnt);
    }

    net_matrix __inverse_() const {
        assert(ln_cnt == col_cnt);
        return net_matrix(inverser(ptr, ln_cnt), ln_cnt, col_cnt);
    }

    net_matrix __tranpose__() const { return net_matrix(transposition(ptr, ln_cnt, col_cnt), col_cnt, ln_cnt); }

    matrix_elem_t __atom__() const {
        assert(ln_cnt == col_cnt && ln_cnt == 1);
        return *ptr;
    }

    net_matrix __abs__() const { return net_matrix(absolute(ptr, elem_cnt), ln_cnt, col_cnt); }

    net_matrix __LU_decompose__() const {
        assert(ln_cnt == col_cnt);
        return net_matrix(LU(ptr, ln_cnt), ln_cnt, col_cnt);
    }

    net_matrix __adjugation__() const {
        assert(ln_cnt == col_cnt);
        return net_matrix(adjugate(ptr, ln_cnt), ln_cnt, col_cnt);
    }

    uint64_t __ranking__() const { return rank(ptr, ln_cnt, col_cnt); }

    net_matrix<net_decimal> __decimal_format__() const {
        if constexpr (std::is_same_v<net_decimal, matrix_elem_t>) return net_matrix(*this);
        else {
            net_matrix<net_decimal> ans(ln_cnt, col_cnt);
            for (auto i = 0ull; i < elem_cnt; ++i) ans.index(i) = *(ptr + i);
            return ans;
        }
    }

    static net_matrix identity(uint64_t dim_cnt) { return net_matrix(init_identity<matrix_elem_t>(dim_cnt), dim_cnt, dim_cnt); }

    __declspec(property(get=__ln_cnt__))       uint64_t      line_count;
    __declspec(property(get=__col_cnt__))      uint64_t      column_count;
    __declspec(property(get=__elem_cnt__))     uint64_t      element_count;
    __declspec(property(get=is_matrix))        bool          verify;
    __declspec(property(get=__determinant__))  matrix_elem_t determinant;
    __declspec(property(get=__inverse_))       net_matrix    inverse;
    __declspec(property(get=__tranpose__))     net_matrix    transpose;
    __declspec(property(get=__atom__))         matrix_elem_t atom;
    __declspec(property(get=__abs__))          net_matrix    abs;
    __declspec(property(get=__LU_decompose__)) net_matrix    LU_decompse;
    __declspec(property(get=__adjugation__))   net_matrix    adjugation;
    __declspec(property(get=__ranking__))      uint64_t      ranking;

    __declspec(property(get=float_point_format)) net_matrix<long double> float_format;
    __declspec(property(get=__decimal_format__)) net_matrix<net_decimal> decimal_format;

    net_matrix operator+(const net_matrix &src) const {
        if (!shape_verify(src)) return net_matrix(*this);
        net_matrix ans;
        ans.ptr = add(ptr, src.ptr, elem_cnt);
        ans.value_assign(src);
        return ans;
    }
    void operator+=(const net_matrix &src) { *this = *this + src; }

    net_matrix operator-(const net_matrix &src) const {
        if (!shape_verify(src)) return net_matrix(*this);
        net_matrix ans;
        ans.ptr = add(ptr, src.ptr, elem_cnt, true);
        ans.value_assign(src);
        return ans;
    }
    void operator-=(const net_matrix &src) { *this = *this - src; }

    net_matrix operator*(const net_matrix &src) const {
        net_matrix ans;
        if (src.elem_cnt == 1) {
            ans.ptr      = mult(ptr, elem_cnt, *src.ptr);
            ans.col_cnt  = col_cnt;
            ans.elem_cnt = elem_cnt;
        } else {
            if (col_cnt != src.ln_cnt) return net_matrix(*this);
            ans.ptr      = mult(ptr, ln_cnt, col_cnt, src.ptr, src.col_cnt);
            ans.col_cnt  = src.col_cnt;
            ans.elem_cnt = ln_cnt * ans.col_cnt;
        }
        ans.ln_cnt = ln_cnt;
        return ans;
    }
    template<typename matrix_elem_para, typename matrix_elem_para_v> friend net_matrix operator*(const matrix_elem_para &para, const net_matrix &src) { return src * net_matrix(para); }
    void operator*=(const net_matrix &src) { *this = *this * src; }

    net_matrix &operator=(const net_matrix &src) {
        value_copy(src);
        return *this;
    }
    net_matrix &operator=(net_matrix &&src) {
        value_move(std::move(src));
        return *this;
    }
    template<typename s_arg> net_matrix &operator=(const net_matrix<s_arg> &src) {
        value_copy(src);
        return *this;
    }

    line_data operator[](uint64_t ln) const {
        assert(ln < ln_cnt);
        return line_data(this, ln, col_cnt);
    }

    // matrix_elem_t *operator[](uint64_t ln) const {
    //     assert(ln < ln_cnt);
    //     return ptr + ln * col_cnt;
    // }

    bool operator==(const net_matrix &src) const {
        if (shape_verify(src)) {
            for (auto i = 0ull; i < elem_cnt; ++i) if (*(ptr + i) != *(src.ptr + i)) return false;
            return true;
        } else return false;
    }

    bool operator!=(const net_matrix &src) const { return !(*this == src); }

    friend std::ostream &operator<<(std::ostream &out, const net_matrix &src) {
        if(src.is_matrix()) for(auto i = 0ull; i < src.elem_cnt; ++i) {
            out << *(src.ptr + i);
            if(!((i + 1) % src.col_cnt || (i + 1) == src.elem_cnt)) out << '\n';
            else out << '\t'; 
        }
        return out;
    }

};

MATRIX_END