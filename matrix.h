MTX_BEGIN

template<typename __mtx_elem, typename __mtx_elem_v> class matrix
{
public:
    struct mtx_strassen
    {
        int ln_curr(int ln);
        void alloc(int ln_cnt, int col_cnt);
        bool convert(__mtx &src, bool free = false);
        __mtx operator[](int ln);
        mtx_strassen operator+(mtx_strassen &right);
        mtx_strassen operator-(mtx_strassen &right);
        mtx_strassen operator*(mtx_strassen &right);
        bool release();
        void print();
        int ln_cnt = 0, ln_cnt_orgn = 0, col_cnt = 0, col_cnt_orgn = 0, ln_begin = 0, ln_end = 0, col_begin = 0, col_end = 0;
        __mtx ptr = nullptr;
    };
    struct mtx_pos
    {
        uint64_t ln = 0, col = 0;
        mtx_pos(uint64_t _ln = 0, uint64_t _col = 0) : ln(_ln), col(_col) {}
        friend std::ostream &operator<<(std::ostream &output, mtx_pos &src)
        {
            output << '(' << src.ln << ", " << src.col << ')';
            return output;
        }
    };
protected:
    void mtx_strassen_child(mtx_strassen &ans, mtx_strassen &src, uint64_t ln_from, uint64_t col_from, uint64_t ln_cnt_child, uint64_t col_cnt_child);
    void mtx_strassen_quartile(mtx_strassen &src, mtx_strassen &src00, mtx_strassen &src01, mtx_strassen &src10, mtx_strassen &src11);
    mtx_strassen mtx_strassen_mult(mtx_strassen &left, mtx_strassen &right, int recursive_gate);
    bool mtx_strassen_mult(__mtx &ans, __mtx &l_src, __mtx &r_src, uint64_t l_ln_cnt, uint64_t l_col_cnt, uint64_t r_ln_cnt, uint64_t r_col_cnt, uint64_t recursive_gate = 32);
    void para_set(uint64_t ln_cnt, uint64_t col_cnt);
    template<typename _Ity = __mtx_elem> void __value_copy(const matrix<_Ity> &src);
    void init_list_mtx(matrix &&_vect, bagrt::net_sequence<matrix> &vect_set, bagrt::net_sequence<mtx_pos> &left_top_pos, uint64_t curr_left, uint64_t curr_top, uint64_t &left_increase, uint64_t &top_increase);
    template<typename _Ity = __mtx_elem> void init_list_mtx(std::initializer_list<std::initializer_list<_Ity>> &&_vect, bagrt::net_sequence<matrix> &vect_set, bagrt::net_sequence<mtx_pos> &left_top_pos, uint64_t curr_left, uint64_t curr_top, uint64_t &left_increase, uint64_t &top_increase);
    struct mtx_extm_data { __mtx_elem ans = 0; bagrt::net_set<mtx_pos> ans_pos; };
    struct col_pointer
    {
    public:
        col_pointer(matrix *_mtx_buf_ptr = nullptr, uint64_t begin_idx = 0, uint64_t _ptr_len = 0);
        __mtx_elem &operator[](uint64_t col);
        ~col_pointer();
    private:
        matrix *_buf_ptr;
        uint64_t _ptr_ln, _ptr_col_cnt;
    };
public:
    matrix(matrix &src);
    matrix(const matrix &src);
    matrix(matrix &&src);
    template<typename _Ity = __mtx_elem> void value_move(matrix<_Ity> &&src);
    template<typename _Ity = __mtx_elem> void value_copy(matrix<_Ity> &src);
    matrix(uint64_t ln_cnt = 0, uint64_t col_cnt = 0, bool rand_mode = false, __mtx_elem &&rand_boundary_first = 0, __mtx_elem &&rand_boundary_second = 0, __mtx_elem &&rand_acc = 1e-8);
    matrix(__mtx &&ptr_src, uint64_t ln_cnt, uint64_t col_cnt);
    matrix(__mtx &ptr_src, uint64_t ln_cnt, uint64_t col_cnt);
    template<typename _Ity = __mtx_elem> matrix(_Ity &&atom);
    matrix(bagrt::net_sequence<matrix> &_vect, bagrt::net_sequence<mtx_pos> &left_top_pos);
    template<typename _Ity = __mtx_elem> matrix(std::initializer_list<std::initializer_list<_Ity>> _vect);
    uint64_t get_ln_cnt();
    uint64_t get_col_cnt();
    uint64_t get_elem_cnt();
    __mtx get_mtx_ptr();
    bool is_matrix();
    bool fill_with(__mtx_elem &&val);
    __mtx_elem get_determinant();
    matrix get_inverser();
    matrix get_transposition();
    mtx_extm_data extremum(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, bool max_flag = true, bool get_extm_pos = true, uint64_t ln_dilation = 0, uint64_t col_dilation = 0);
    __mtx_elem get_atom();
    matrix child(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilation = 0, uint64_t col_dilation = 0);
    matrix rotate_rect(bool clockwise = true);
    matrix mirror_flip(bool symmetry_vertical = true);
    bool shape_valid(uint64_t ln_cnt, uint64_t col_cnt);
    bool shape_valid(matrix &src);
    matrix reshape(uint64_t ln_cnt, uint64_t col_cnt);
    matrix reshape(matrix &as_val);
    __mtx_elem elem_sum(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilation = 0, uint64_t col_dilation = 0);
    __mtx_elem elem_sum();
    matrix abs();
    matrix elem_cal_opt(matrix &r_val, uint64_t opt_idx);
    matrix elem_cal_opt(__mtx_elem &&para, uint64_t opt_idx);
    matrix broadcast_add(__mtx_elem &&val, bool nega = false);
    matrix pad(uint64_t ln_t = 0, uint64_t col_r = 0, uint64_t ln_b = 0, uint64_t col_l = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0);
    matrix crop(uint64_t ln_t = 0, uint64_t col_r = 0, uint64_t ln_b = 0, uint64_t col_l = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0);
    matrix round_fit();
    matrix LU_decompose();
    matrix swap_dir_elem(uint64_t l_idx, uint64_t r_idx, bool is_ln = true);
    matrix get_adjugate();
    matrix mult_strassen(matrix &r_src, uint64_t recursive_gate = 32);
    uint64_t get_rank();
    matrix<_BAGRT decimal> get_mtx_dec();
    matrix<long double> get_mtx_f();
    __mtx_elem &pos_idx(uint64_t idx);
    matrix operator+(matrix &val);
    matrix operator-(matrix &val);
    void operator+=(matrix &val);
    void operator-=(matrix &val);
    matrix operator*(matrix &val);
    void operator*=(matrix &val);
    matrix operator*(__mtx_elem &&val);
    void operator*=(__mtx_elem &&val);
    template<typename _Ity = __mtx_elem> void operator=(const matrix<_Ity> &val);
    template<typename _Ity = __mtx_elem> void operator=(matrix<_Ity> &src);
    template<typename _Ity = __mtx_elem> void operator=(matrix<_Ity> &&src);
    bool operator==(matrix &val);
    bool operator!=(matrix &val);
    col_pointer operator[](uint64_t ln);
    void reset();
    ~matrix();
protected:
    uint64_t mtx_ln_cnt = 0, mtx_col_cnt = 0, mtx_elem_cnt = 0;
    __mtx mtx_ptr = nullptr;
public:
    __declspec(property(get=get_ln_cnt)) uint64_t line_count;
    __declspec(property(get=get_col_cnt)) uint64_t column_count;
    __declspec(property(get=get_elem_cnt)) uint64_t element_count;
    __declspec(property(get=abs)) matrix absolute;
    __declspec(property(get=round_fit)) matrix round_elements;
    __declspec(property(get=LU_decompose)) matrix LU;
    __declspec(property(get=get_adjugate)) matrix adjugate;
    __declspec(property(get=get_rank)) uint64_t rank;
    __declspec(property(get=get_determinant)) __mtx_elem determinant;
    __declspec(property(get=get_inverser)) matrix inverser;
    __declspec(property(get=get_transposition)) matrix transposition;
    __declspec(property(get=get_atom)) __mtx_elem atom;
    __declspec(property(get=get_mtx_ptr)) __mtx pointer;
    __declspec(property(get=get_mtx_dec)) matrix<_BAGRT decimal> decimal_matrix;
    __declspec(property(get=get_mtx_f)) matrix<long double> float_matrix;
    friend std::ostream &operator<<(std::ostream &output, matrix &src)
    {
        if(src.is_matrix()) for(auto i=0; i<src.mtx_elem_cnt; ++i)
        {
            output << src.mtx_ptr[i];
            if(!((i+1)%src.mtx_col_cnt || (i+1)==src.mtx_elem_cnt)) output << std::endl;
            else output << '\t'; 
        }
        return output;
    }
    friend matrix operator*(__mtx_elem &&val, matrix &src) { return src * std::move(val); }
    static matrix init_E(uint64_t dms);
};

MTX_END

#include "matrix.hpp"