MATRIX_BEGIN

std::default_random_engine  rand_e;

struct mtx_pos
{
    int ln = 0;
    int col = 0;
    friend std::ostream& operator<<(std::ostream &input, mtx_pos &val)
    {
        input << '(' << val.ln << ", " << val.col << ')';
        return input;
    }
};

struct mtx_info
{
    MATRIX mtx_val = nullptr;
    uint64_t ln_cnt = 0;
    uint64_t col_cnt = 0;
};

struct mtx_extm
{
    double val = 0;
    bagrt::net_list<mtx_pos> pos_list;
};

template<typename...args> void mtx_reset(args&...mtx_vals) { bagrt::reset_ptr(mtx_vals...); }

bool mtx_pos_valid(uint64_t ln, uint64_t col, uint64_t ln_cnt, uint64_t col_cnt)
{
    if (ln<ln_cnt && col<col_cnt) return true;
    else return false;
}

mtx_pos mtx_elem_pos(uint64_t idx, uint64_t col_cnt)
{
    mtx_pos pos;
    pos.ln = idx / col_cnt;
    pos.col = idx % col_cnt;
    return pos;
}

uint64_t mtx_elem_pos(uint64_t curr_ln, uint64_t curr_col, uint64_t col_cnt, uint64_t ln_from = 0, uint64_t ln_to = 0, uint64_t col_from = 0, uint64_t col_to = 0, uint64_t ln_dilate = 0, uint64_t col_dilate = 0)
{
    if(ln_from || ln_to || col_from || col_to || ln_dilate || col_dilate)
    {
        auto ln = ln_from + curr_ln * (1 + ln_dilate),
            col = col_from + curr_col * (1 + col_dilate);
        return ln * col_cnt + col;
    }
    return curr_ln * col_cnt + curr_col;
}

MATRIX mtx_init(uint64_t elem_cnt)
{
    auto mtx_ptr = std::make_unique<double[]>(elem_cnt);
    for(auto i=0; i<elem_cnt; ++i) mtx_ptr[i] = 0;
    return mtx_ptr;
}
MATRIX mtx_init(uint64_t ln_cnt, uint64_t col_cnt) {return mtx_init(ln_cnt * col_cnt);}

MATRIX mtx_init_E(uint64_t dms)
{
    auto elem_cnt = dms * dms;
    auto mtx_ptr = mtx_init(elem_cnt);
    for(auto i=0; i<elem_cnt; ++i) mtx_ptr[i*(dms+1)] = 1;
    return mtx_ptr;
}

bool mtx_equal(MATRIX &first, MATRIX &second, uint64_t elem_cnt)
{
    for(auto i=0; i<elem_cnt; ++i) if(first[i] != second[i])
    {
        auto dim = mtx_elem_pos(i, std::pow(elem_cnt, 0.5));
        std::cout << dim.ln << ", " << dim.col << std::endl;
        return false;
    }
    return true;
}

bool mtx_equal(MATRIX &first, MATRIX &second, uint64_t ln_cnt, uint64_t col_cnt) { return mtx_equal(first, second, ln_cnt*col_cnt); }

void mtx_abs(MATRIX &mtx_src, uint64_t elem_cnt) {if(mtx_src && elem_cnt) for(auto i=0; i<elem_cnt; ++i) mtx_src[i] = std::abs(mtx_src[i]);}

MATRIX mtx_init_rand(uint64_t elem_cnt, double boundry_first = 0, double boundry_second = 0, double acc = 1e-5)
{
    MATRIX mtx_val = mtx_init(elem_cnt);
    for (int i = 0; i < elem_cnt; ++i)
        if(boundry_first==boundry_second) mtx_val[i] = (((double)rand_e() / (double)rand_e._Max) - 0.5) * 2.0;
        else mtx_val[i] = bagrt::random_number(boundry_first, boundry_second, false, acc); 
    return mtx_val;
}

void mtx_child_vec(MATRIX &ans, uint64_t &child_ln_cnt, uint64_t &child_col_cnt, MATRIX &mtx_src, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_dilation = 0, uint64_t col_dilation = 0)
{
    if(mtx_src && col_cnt && ln_cnt &&
        from_ln>=0 && to_ln>=from_ln && ln_cnt>to_ln &&
        from_col>=0 && to_col>=from_col && col_cnt>to_col)
    {
        child_ln_cnt = (bagrt::num_cnt(from_ln, to_ln) + ln_dilation) / (ln_dilation + 1);
        child_col_cnt = (bagrt::num_cnt(from_col, to_col) + col_dilation) / (col_dilation + 1);
        if(!ans) ans = mtx_init(child_ln_cnt, child_col_cnt);
        for(auto i=0; i<child_ln_cnt; ++i)
            for(auto j=0; j<child_col_cnt; ++j)
            {
                auto curr_orgn_ln = from_ln + i * (1 + ln_dilation);
                auto curr_orgn_col = from_col + j * (1 + col_dilation);
                ans[mtx_elem_pos(i, j, child_col_cnt)] = mtx_src[mtx_elem_pos(curr_orgn_ln, curr_orgn_col, col_cnt)];
            }
    }
}

mtx_info mtx_child_vec(MATRIX &mtx_src, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_dilation = 0, uint64_t col_dilation = 0)
{
    mtx_info mtx_info;
    mtx_child_vec(mtx_info.mtx_val, mtx_info.ln_cnt, mtx_info.col_cnt, mtx_src, from_ln, to_ln, from_col, to_col, ln_cnt, col_cnt, ln_dilation, col_dilation);
    return mtx_info;
}

double mtx_atom(MATRIX &mtx_val) {return mtx_val[IDX_ZERO];}

bool mtx_refresh(MATRIX &mtx_val, uint64_t elem_cnt)
{
    if(mtx_val)
    {
        for(auto i=0; i<elem_cnt; ++i) if(mtx_val[i]) mtx_val[i] = 0;
        return true;
    }
    else return false;
}
bool mtx_refresh(MATRIX &mtx_val, uint64_t ln_cnt, uint64_t col_cnt) {return mtx_refresh(mtx_val, ln_cnt * col_cnt);}

void mtx_show(MATRIX &mtx_val, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(mtx_val)
    {
        auto elem_amt = ln_cnt * col_cnt;
        for (auto i = 0; i < elem_amt; ++i)
        {
            std::cout << mtx_val[i];
            if ((i + 1) % col_cnt) std::cout << '\t';
            else std::cout << std::endl;
        }
    }
}

double mtx_det(MATRIX &mtx_val, int dms)
{
    if(mtx_val)
    {
        if (dms == 1)
        return mtx_val[IDX_ZERO];
        auto ac = mtx_init((dms - 1) * (dms - 1));
        int mov = 0;
        double sum = 0.0;
        for (int arow = 0; arow < dms; arow++)
        {
            for (int brow = 0; brow < dms - 1; brow++)
            {
                mov = arow > brow ? 0 : 1;
                for (int j = 0; j < dms - 1; ++j)
                    ac[brow * (dms - 1) + j] = mtx_val[(brow + mov) * dms + j + 1];
            }
            int flag = (arow % 2 == 0 ? 1 : -1);
            sum += flag * mtx_val[arow * dms] * mtx_det(ac, dms - 1);
        }
        mtx_reset(ac);
        return sum;
    }
    else return NAN;    
}

mtx_info mtx_pad(MATRIX &mtx_val, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_t = 0, uint64_t col_r = 0, uint64_t ln_b = 0, uint64_t col_l = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0)
{
    mtx_info mtx_res;
    if(mtx_val && ln_cnt && col_cnt)
    {
        mtx_res.ln_cnt = ln_t + ln_b + ln_cnt + (ln_cnt - 1) * ln_dist;
        mtx_res.col_cnt = col_r + col_l + col_cnt + (col_cnt - 1) * col_dist;
        auto elem_cnt = mtx_res.ln_cnt * mtx_res.col_cnt;
        mtx_res.mtx_val = mtx_init(elem_cnt);
        for(auto i=0; i<ln_cnt; ++i)
            for(auto j=0; j<col_cnt; ++j)
            {
                auto curr_pad_no =  mtx_elem_pos(ln_t+i*(ln_dist+1), col_l+j*(col_dist+1), mtx_res.col_cnt),
                    curr_origin_no = mtx_elem_pos(i, j, col_cnt);
                mtx_res.mtx_val[curr_pad_no] = mtx_val[curr_origin_no];
            }
    }
    return mtx_res;
}

mtx_info mtx_crop(MATRIX &mtx_val, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_t = 0, uint64_t col_r = 0, uint64_t ln_b = 0, uint64_t col_l = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0)
{
    mtx_info mtx_res;
    if(mtx_val && ln_cnt && col_cnt)
    {
        mtx_res.ln_cnt = (ln_cnt - (ln_t + ln_b) + ln_dist) / (ln_dist + 1);
        mtx_res.col_cnt = (col_cnt - (col_r + col_l) + col_dist) / (col_dist + 1);
        auto elem_cnt = mtx_res.ln_cnt * mtx_res.col_cnt;
        mtx_res.mtx_val = mtx_init(elem_cnt);
        for(auto i=0; i<mtx_res.ln_cnt; ++i)
            for(auto j=0; j<mtx_res.col_cnt; ++j)
            {
                auto curr_res_no = mtx_elem_pos(i, j, mtx_res.col_cnt),
                    curr_origin_no = mtx_elem_pos(ln_t+i*(ln_dist+1), col_l+j*(col_dist+1), col_cnt);
                mtx_res.mtx_val[curr_res_no] = mtx_val[curr_origin_no];
            }
    }
    return mtx_res;
}

mtx_extm mtx_extm_val(MATRIX &mtx_val, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_dilation = 0, uint64_t col_dilation = 0, bool max_flag = true)
{
    mtx_extm mtx_info;
    if(mtx_val && col_cnt && ln_cnt &&
        from_ln>=0 && to_ln>=from_ln && ln_cnt>to_ln &&
        from_col>=0 && to_col>=from_col && col_cnt>to_col)
    {
        mtx_info.val = mtx_val[mtx_elem_pos(from_ln, from_col, col_cnt)];
        for(auto i=from_ln; i<=to_ln; ++(i+=ln_dilation))
            for(auto j=from_col; j<=to_col; ++(j+=col_dilation))
            {
                auto curr_no = mtx_elem_pos(i, j, col_cnt);
                mtx_pos curr_pos;
                curr_pos.ln = i;
                curr_pos.col = j;
                if((max_flag&&mtx_val[curr_no]>mtx_info.val) || (!max_flag&&mtx_val[curr_no]<mtx_info.val))
                {
                    mtx_info.val = mtx_val[curr_no];
                    mtx_info.pos_list.reset();
                    mtx_info.pos_list.emplace_back(curr_pos);
                }
                else if(mtx_val[curr_no] == mtx_info.val) mtx_info.pos_list.emplace_back(curr_pos);
                else continue;
            }
    }
    return mtx_info;
}

double mtx_sum(MATRIX &mtx_val, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_dilation = 0, uint64_t col_dilation = 0)
{
    auto sum = 0.0;
    if(mtx_val && col_cnt && ln_cnt &&
        from_ln>=0 && to_ln>=from_ln && ln_cnt>to_ln &&
        from_col>=0 && to_col>=from_col && col_cnt>to_col)
    {
        auto mtx_ch = mtx_child_vec(mtx_val, from_ln, to_ln, from_col, to_col, ln_cnt, col_cnt, ln_dilation, col_dilation);
        auto elem_cnt = mtx_ch.col_cnt * mtx_ch.ln_cnt;
        for(auto i=0; i<elem_cnt; ++i) sum += mtx_ch.mtx_val[i];
    }
    return sum;
}

MATRIX mtx_copy(MATRIX &mtx_src, uint64_t elem_cnt)
{
    auto mtx_val = mtx_init(elem_cnt);
    if(mtx_src) for(auto i=0; i<elem_cnt; ++i) mtx_val[i] = mtx_src[i];
    return mtx_val;
}
MATRIX mtx_copy(MATRIX &mtx_src, uint64_t ln_cnt, uint64_t col_cnt) {return mtx_copy(mtx_src, ln_cnt * col_cnt);}

MATRIX mtx_from_str(std::string &mtx_str, uint64_t &line, uint64_t &column)
{
    column = 0;
    line = 0;
    int cnt = 0;
    // get column amount
    while (mtx_str[cnt] != '\n' && mtx_str[cnt] != '\0')
    {
        if (mtx_str[cnt] == '.' || mtx_str[cnt] == '-' || mtx_str[cnt] == '+' ||
            (mtx_str[cnt] >= '0' && mtx_str[cnt] <= '9') || mtx_str[cnt] == '/')
            if (mtx_str[cnt - 1] == ' ' || mtx_str[cnt - 1] == '\t' || cnt - 1 < 0)
                column++;
        cnt++;
    }
    for (int i = 0; i < mtx_str.length(); ++i) if (mtx_str[i] == '\n') line++;
    return bagrt::extract_number(mtx_str);
}

MATRIX mtx_transposition(MATRIX &mtx_src, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(mtx_src && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        auto mtx_trans = mtx_init(elem_cnt);
        for (auto i = 0; i < ln_cnt; ++i)
            for (auto j = 0; j < col_cnt; ++j)
                mtx_trans[mtx_elem_pos(j, i, ln_cnt)] = mtx_src[mtx_elem_pos(i, j, col_cnt)];
        return mtx_trans;
    }
    else return MATRIX_NULL;
}

MATRIX mtx_swap_elem(MATRIX &mtx_val, uint64_t l_pos, uint64_t r_pos, uint64_t ln_cnt, uint64_t col_cnt, bool is_ln = true)
{
    auto ans = mtx_copy(mtx_val, ln_cnt, col_cnt);
    r_pos < l_pos ? std::swap(l_pos, r_pos) : 0;
    if(is_ln) for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j)
    {
        auto curr_ln = i;
        if(i == l_pos) curr_ln = r_pos;
        else if(i==r_pos) curr_ln = l_pos;
        else curr_ln = i;
        ans[mtx_elem_pos(i, j, col_cnt)] = mtx_val[mtx_elem_pos(curr_ln, j, col_cnt)];
    }
    else for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j)
    {
        auto curr_col = j;
        if(j == l_pos) curr_col = r_pos;
        else if(j==r_pos) curr_col = l_pos;
        else curr_col = j;
        ans[mtx_elem_pos(i, j, col_cnt)] = mtx_val[mtx_elem_pos(i, curr_col, col_cnt)];
    }
    return ans;
}

MATRIX mtx_add(MATRIX &l_mtx, MATRIX &r_mtx, uint64_t ln_cnt, uint64_t col_cnt, bool nega = false)
{
    if(l_mtx && r_mtx && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        auto mtx_val = mtx_init(elem_cnt);
        for(auto i=0; i<elem_cnt; ++i)
            if(nega) mtx_val[i] = l_mtx[i] - r_mtx[i];
            else mtx_val[i] = l_mtx[i] + r_mtx[i];
        return mtx_val;
    }
    else return MATRIX_NULL;
}

MATRIX mtx_mult(MATRIX &l_mtx, MATRIX &r_mtx, uint64_t l_ln_cnt, uint64_t l_col_cnt, uint64_t r_ln_cnt, uint64_t r_col_cnt)
{
    if(l_mtx && r_mtx && l_ln_cnt && l_col_cnt && r_ln_cnt && r_col_cnt && l_col_cnt==r_ln_cnt)
    {
        auto res_mtx = mtx_init(l_ln_cnt, r_col_cnt);
        for (int i=0; i<l_ln_cnt; ++i) for(auto j=0; j<l_col_cnt; ++j)
        {
            auto coe = l_mtx[mtx_elem_pos(i, j, l_col_cnt)];
            for(auto k=0; k<r_col_cnt; ++k)
            res_mtx[mtx_elem_pos(i, k, r_col_cnt)] += coe * r_mtx[mtx_elem_pos(j, k, r_col_cnt)];
        }
        return res_mtx;
    }
    else return MATRIX_NULL;
}
MATRIX mtx_mult(MATRIX &mtx_src, double val, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(mtx_src && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        auto res_mtx = mtx_copy(mtx_src, elem_cnt);
        for(auto i=0; i<elem_cnt; ++i) res_mtx[i] *= val;
        return res_mtx;
    }
    else return MATRIX_NULL;
}

MATRIX mtx_elem_cal_opt(MATRIX &l_mtx, MATRIX &r_mtx, uint64_t ln_cnt, uint64_t col_cnt, uint64_t opt_idx)
{
    if(l_mtx && r_mtx && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        MATRIX mtx_val = mtx_init(elem_cnt);
        for(auto i=0; i<elem_cnt; ++i)
            switch (opt_idx)
            {
            case MATRIX_ELEM_MULT:
                mtx_val[i] = l_mtx[i] * r_mtx[i];
                break;
            case MATRIX_ELEM_DIV:
                if(r_mtx[i])
                {
                    mtx_val[i] = l_mtx[i] / r_mtx[i];
                    break;
                }
                else return MATRIX_NULL;
            default: return MATRIX_NULL;
            }
        return mtx_val;
    }
    else return MATRIX_NULL;
}

MATRIX mtx_elem_cal_opt(MATRIX &mtx_val, double para, uint64_t ln_cnt, uint64_t col_cnt, uint64_t opt_idx)
{
    if(mtx_val && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        MATRIX mtx_res;
        if(para != 1)
        {
            mtx_res = mtx_init(ln_cnt, col_cnt);
            for(auto i=0; i<elem_cnt; ++i) switch (opt_idx)
            {
            case MATRIX_ELEM_POW:
                mtx_res[i] = std::pow(mtx_val[i], para);
                break;
            case MATRIX_ELEM_DIV:
                if(para)
                {
                    mtx_res[i] = mtx_val[i] / para;
                    break;
                }
            default: return MATRIX_NULL;
            }
        }
        else mtx_res = mtx_copy(mtx_val, ln_cnt, col_cnt);
        return mtx_res;
    }
    else return MATRIX_NULL;
}

MATRIX mtx_broadcast_add(MATRIX &mtx_val, uint64_t elem_cnt, double val)
{
    if(mtx_val && elem_cnt)
    {
        auto mtx_res = mtx_copy(mtx_val, elem_cnt);
        for(auto i=0; i<elem_cnt; ++i) mtx_res[i] += val;
        return mtx_res;
    }
    else return MATRIX_NULL;
}

MATRIX mtx_broadcast_subtract(MATRIX &mtx_val, uint64_t elem_cnt, double val, bool val_subtrahend = true)
{
    if(mtx_val && elem_cnt)
    {
        auto mtx_res = mtx_copy(mtx_val, elem_cnt);
        for(auto i=0; i<elem_cnt; ++i)
            if(val_subtrahend) mtx_res[i] -= val;
            else mtx_res[i] = val - mtx_res[i];
        return mtx_res;
    }
    else return MATRIX_NULL;
}

MATRIX mtx_LU(MATRIX &mtx_val, int dms)
{
    if(mtx_val && dms)
    {
        auto lu = mtx_copy(mtx_val, dms, dms);
        auto e = mtx_init_E(dms);
        auto l = mtx_copy(e, dms, dms);
        for (int i = 0; i < dms - 1; ++i)
        {
            auto tool = mtx_copy(e, dms, dms);
            for (int j = dms - 1; j > i; --j) tool[mtx_elem_pos(j, i, dms)] = (-1.0) * lu[mtx_elem_pos(j, i, dms)] / lu[mtx_elem_pos(i, i, dms)];
            for (int j = dms - 1; j > i; --j) l[mtx_elem_pos(j, i, dms)] = (-1.0) * tool[mtx_elem_pos(j, i, dms)];
            for (int j = dms - 1; j > i; --j) lu[mtx_elem_pos(j, i, dms)] = (-1.0) * tool[mtx_elem_pos(j, i, dms)];
            lu = mtx_mult(tool, lu, dms, dms, dms, dms);
            mtx_reset(tool);
        }
        // get the inversion matrix
        for (int i = 0; i < dms; ++i)
            for (int j = 0; j < i; ++j)
                lu[mtx_elem_pos(i, j, dms)] = l[mtx_elem_pos(i, j, dms)];
        mtx_reset(e, l);
        return lu;
    }
    else return MATRIX_NULL;
}

MATRIX mtx_equation(MATRIX &coefficient, MATRIX &b, int dms)
{
    if(coefficient && b && dms)
    {
        auto lu = mtx_LU(coefficient, dms);
        auto e = mtx_init_E(dms);
        auto l = mtx_copy(e, dms, dms);
        auto u = mtx_copy(e, dms, dms);
        for (int i = 0; i < dms; ++i)
            for (int j = 0; j < dms; ++j)
                if (i <= j) u[mtx_elem_pos(i, j, dms)] = lu[mtx_elem_pos(i, j, dms)];
                else l[mtx_elem_pos(i, j, dms)] = lu[mtx_elem_pos(i, j, dms)];
        auto y = mtx_init(dms), x = mtx_init(dms);
        for (int i = 0; i < dms; ++i) if (i)
        {
            double temp = 0.0;
            for (int j=0; j < i; ++j) temp += y[j] * l[mtx_elem_pos(i, j, dms)] * 1.0;
            y[i] = b[i] - temp;
        }
        else y[i] = b[i];
        for (int i=dms; i>0; --i) if (i - dms)
        {
            double temp = 0.0;
            for (int n=i-1; n<dms-1; ++n) temp += u[mtx_elem_pos(i-1, n+1, dms)] * x[n+1];
            x[i-1] = (y[i-1] - temp) / u[mtx_elem_pos(i-1, i-1, dms)];
        }
        else x[dms-1] = y[dms-1] / u[mtx_elem_pos(dms-1, dms-1, dms)];
        mtx_reset(lu, e, l, u, y);
        return x;
    }
    else return MATRIX_NULL;
}

MATRIX mtx_cofactor(MATRIX &mat_val, int ln, int col, int dms)
{
    auto ans_dms = dms - 1;
    auto ans = mtx_init(ans_dms, ans_dms);
    auto elem_cnt = dms * dms;
    for(auto i=0; i<elem_cnt; ++i)
    {
        auto curr_dim = mtx_elem_pos(i, dms);
        auto curr_ln = curr_dim.ln, curr_col = curr_dim.col;
        if(curr_ln!=ln && curr_col!=col)
        {
            auto curr_ans_ln = curr_ln, curr_ans_col = curr_col;
            if(curr_ln > ln) -- curr_ans_ln;
            if(curr_col > col) -- curr_ans_col;
            ans[mtx_elem_pos(curr_ans_ln, curr_ans_col, ans_dms)] = mat_val[i];
        }
    }
    return ans;
}

double mtx_algebraic_cofactor(MATRIX &mat_val, int ln, int col, int dms) { return mtx_det(mtx_cofactor(mat_val, ln, col, dms), dms-1); }

MATRIX mtx_adjugate(MATRIX &mat_val, int dms)
{
    auto elem_cnt = dms * dms;
    auto ans = mtx_init(elem_cnt);
    for(auto i=0; i<elem_cnt; ++i)
    {
        auto curr_dim = mtx_elem_pos(i, dms);
        auto coe = 1;
        if((curr_dim.ln+curr_dim.col) % 2) coe = -1;
        ans[mtx_elem_pos(curr_dim.col, curr_dim.ln, dms)] = coe * mtx_algebraic_cofactor(mat_val, curr_dim.ln, curr_dim.col, dms);
    }
    return ans;
}

MATRIX mtx_inverser(MATRIX &mat_val, int dms)
{
    if(mat_val && dms)
    {
        auto elem_cnt = dms * dms;
        auto ans = mtx_init(elem_cnt);
        auto val_det = mtx_det(mat_val, dms);
        auto val_adj = mtx_adjugate(mat_val, dms);
        for(auto i=0; i<elem_cnt; ++i) ans[i] = val_adj[i] / val_det;
        return ans;
    } return MATRIX_NULL;
}

double mtx_max_eigenvalue(MATRIX &mtx_val, int dms, MATRIX &w, double error = 1e-5, double init_elem = 1)
{
    if(mtx_val && dms)
    {
        double lambda = 0.0, lambda_temp = 0.0;
        auto x = mtx_init(dms), w = mtx_init(dms);
        for (int i=0; i<dms; ++i) x[i] = init_elem;
        do
        {
            lambda = lambda_temp;
            auto temp = mtx_init(dms);
            for (int i=0; i<dms; ++i) temp[i] = x[i];
            auto temp_abs =  mtx_copy(temp, dms);
            mtx_abs(temp_abs, dms);
            double max = mtx_extm_val(temp_abs, 0, dms-1, 0, 0, dms, 1).val;
            for (int i=0; i<dms; ++i) x[i] /= max;
            for (int i=0; i<dms; ++i) w[i] = x[i];
            double sum = 0.0;
            for (int i=0; i<dms; ++i) sum += w[i];
            for (int i=0; i<dms; ++i) w[i] /= sum;
            x = mtx_mult(mtx_val, x, dms, dms, dms, 1);
            for (int i=0; i<dms; ++i) temp[i] = x[i];
            lambda_temp = mtx_extm_val(temp, 0, dms-1, 0, 0, dms, 1).val;
            mtx_reset(temp, temp_abs);
        } while (std::abs(lambda_temp - lambda) > error);
        mtx_reset(x);
        return lambda;
    }
    else return NAN;
}

MATRIX mtx_jacobi_iterate(MATRIX &coefficient, MATRIX &b, int dms, double error = 1e-5, double init_elem = 1)
{
    auto x = mtx_init(dms), temp = mtx_init(dms);
    for (int i=0; i<dms; ++i) temp[i] = init_elem;
    bool flag = false;
    do
    {
        flag = false;
        for (int i=0; i<dms; ++i) x[i] = temp[i];
        for (int i=0; i<dms; ++i)
        {
            double sum = 0.0;
            for (int j=0; j<dms; ++j) if (i != j) sum += coefficient[mtx_elem_pos(i, j, dms)] * x[j];
            temp[i] = (b[i] - sum) / coefficient[mtx_elem_pos(i, i, dms)];
        }
        for (int i=0; i<dms; ++i) if (std::abs(x[i] - temp[i]) >= error)
        {
            flag = true;
            break;
        }
    } while (flag);
    mtx_reset(temp);
    return x;
}

MATRIX mtx_rotate_rect(MATRIX &mtx_val, int ln_cnt, int col_cnt, bool clock_wise = true)
{
    auto rot_rec_mtx = mtx_init(col_cnt, ln_cnt);
    for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j)
        if(clock_wise)  rot_rec_mtx[mtx_elem_pos(j, ln_cnt-i-1, ln_cnt)] = mtx_val[mtx_elem_pos(i, j, col_cnt)];
        else rot_rec_mtx[mtx_elem_pos(col_cnt-j-1, i, ln_cnt)] = mtx_val[mtx_elem_pos(i, j, col_cnt)];
    return rot_rec_mtx;
}

MATRIX mtx_mirror_flip(MATRIX &mtx_src, uint64_t ln_cnt, uint64_t col_cnt, bool symmetry_vertical = true)
{
    auto mtx_val = mtx_init(ln_cnt, col_cnt);
    for(auto i=0; i<ln_cnt; ++i)
        for(int j=0; j<col_cnt; ++j) if(symmetry_vertical)
        {
            auto src_pos = col_cnt-1-j;
            if(src_pos != j)mtx_val[mtx_elem_pos(i, j, col_cnt)] = mtx_src[mtx_elem_pos(i, src_pos, col_cnt)];
        }
        else
        {
            auto src_pos = ln_cnt-1-i;
            if(src_pos != i)mtx_val[mtx_elem_pos(i, j, col_cnt)] = mtx_src[mtx_elem_pos(src_pos, j, col_cnt)];
        }
    return mtx_val;
}

uint64_t mtx_rank(MATRIX &mtx_val, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(ln_cnt && col_cnt)
    {
        MATRIX ans;
        if(ln_cnt > col_cnt)
        {
            std::swap(ln_cnt, col_cnt);
            ans = mtx_transposition(mtx_val, ln_cnt, col_cnt);
        }
        else ans = mtx_copy(mtx_val, ln_cnt, col_cnt);
        auto rank_val = 0;
        for(auto i=0; i<col_cnt; ++i)
        {
            bool elim_flag = true;
            for(auto j=rank_val; j<ln_cnt; ++j)
            {
                auto elim_val = ans[mtx_elem_pos(j, i, col_cnt)];
                if(elim_val)
                {
                    if(elim_flag)
                    {
                        if(elim_val != 1) for(auto k=i; k<col_cnt; k++) ans[mtx_elem_pos(j, k, col_cnt)] /= elim_val;
                        if(j!=i) ans = mtx_swap_elem(ans, i, j, ln_cnt, col_cnt);
                        ++ rank_val;
                        elim_flag = false;
                    }
                    else for(auto k=i; k<col_cnt; ++k) ans[mtx_elem_pos(j, k, col_cnt)] -= ans[mtx_elem_pos(i, k, col_cnt)] * elim_val;
                }
            }
        }
        mtx_reset(ans);
        return rank_val;
    }
    else return NAN;
}

// Unsafe
struct mtx_strassen
{
    int ln_cnt = 0, ln_cnt_orgn = 0, col_cnt = 0, col_cnt_orgn = 0;
    int ln_begin = 0, ln_end = 0, col_begin = 0, col_end = 0;
    double *ptr = nullptr;
    int ln_curr(int ln) { return ((ln + ln_begin) * col_cnt_orgn + col_begin); }
    void alloc(int ln_cnt, int col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        ptr = new double[elem_cnt];
        memset(ptr,(double)0, elem_cnt*sizeof(double));
        this->ln_cnt = ln_cnt; this->ln_cnt_orgn = ln_cnt;
        this->col_cnt = col_cnt; this->col_cnt_orgn = col_cnt;
    }
    MATRIX convert(bool free = false)
    {
        auto ans = mtx_init(ln_cnt, col_cnt);
        for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j) ans[mtx_elem_pos(i, j, col_cnt)] = ptr[ln_curr(i)+j];
        if(free) release();
        return ans;
    }
    double *operator[](int ln) { return ptr + ln_curr(ln); }
    mtx_strassen operator+(mtx_strassen &right)
    {
        mtx_strassen ans; ans.alloc(ln_cnt, col_cnt);
        for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j) ans.ptr[ans.ln_curr(i)+j] = ptr[ln_curr(i)+j] + right.ptr[right.ln_curr(i)+j];
        return ans;
    }
    mtx_strassen operator-(mtx_strassen &right)
    {
        mtx_strassen ans; ans.alloc(ln_cnt, col_cnt);
        for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j) ans.ptr[ans.ln_curr(i)+j] = ptr[ln_curr(i)+j] - right.ptr[right.ln_curr(i)+j];
        return ans;
    }
    mtx_strassen operator*(mtx_strassen &right)
    {
        mtx_strassen ans; ans.alloc(ln_cnt, right.col_cnt);
        for(int i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j)
        {
            auto coe = ptr[ln_curr(i)+j];
            for(auto k=0; k<ans.col_cnt; ++k) ans.ptr[ans.ln_curr(i)+k] += coe * right.ptr[right.ln_curr(j)+k];
        }
        return ans;
    }
    void release() { if(ln_cnt==ln_cnt_orgn && col_cnt==col_cnt_orgn && ptr) { delete [] ptr; ptr = nullptr; } }
    void show()
    {
        for(auto i=0; i<ln_cnt; ++i)
        {
            for(auto j=0; j<col_cnt; ++j) std::cout << ptr[ln_curr(i)+j] << '\t';
            std::cout << std::endl;
        }
    }
};

// Unsafe
void mtx_strassen_child(mtx_strassen &ans, mtx_strassen &src, uint64_t ln_from, uint64_t col_from, uint64_t ln_cnt_child, uint64_t col_cnt_child)
{
    ans.ptr = src.ptr; ans.ln_cnt_orgn = src.ln_cnt_orgn; ans.col_cnt_orgn = src.col_cnt_orgn;
    ans.ln_cnt = ln_cnt_child;
    ans.col_cnt = col_cnt_child;
    ans.ln_begin = src.ln_begin + ln_from; ans.ln_end = src.ln_begin + ln_cnt_child - 1;
    ans.col_begin = src.col_begin + col_from; ans.ln_end = src.col_begin + col_cnt_child - 1;
}

// Unsafe
void mtx_strassen_quartile(mtx_strassen &src, mtx_strassen &src00, mtx_strassen &src01, mtx_strassen &src10, mtx_strassen &src11)
{
    auto dms_sub = src.ln_cnt / 2;
    mtx_strassen_child(src00, src, 0, 0, dms_sub, dms_sub);
    mtx_strassen_child(src01, src, 0, dms_sub, dms_sub, dms_sub);
    mtx_strassen_child(src10, src, dms_sub, 0, dms_sub, dms_sub);
    mtx_strassen_child(src11, src, dms_sub, dms_sub, dms_sub, dms_sub);
}

// Unsafe
mtx_strassen mtx_strassen_mult(mtx_strassen &left, mtx_strassen &right, int recursive_gate)
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

MATRIX mtx_strassen_mult(MATRIX &l_mtx, MATRIX &r_mtx, uint64_t l_ln_cnt, uint64_t l_col_cnt, uint64_t r_ln_cnt, uint64_t r_col_cnt, uint64_t recursive_gate = 32)
{
    auto r_p_l = l_col_cnt, b_p_l = l_ln_cnt, r_p_r = r_col_cnt, b_p_r = r_ln_cnt,
        mx_len = bagrt::num_extreme({r_p_l, b_p_l, r_p_r, b_p_r});
    if(mx_len > recursive_gate)
    {
        auto pad_val = bagrt::num_pow_pad_cnt(mx_len, 2),
            pad_ans = mx_len + pad_val;
        r_p_l = pad_ans - r_p_l; b_p_l = pad_ans - b_p_l;
        r_p_r = pad_ans - r_p_r; b_p_r = pad_ans - b_p_r;
        auto mtx_s_l = mtx_pad(l_mtx, l_ln_cnt, l_col_cnt, 0, r_p_l, b_p_l).mtx_val,
            mtx_s_r = mtx_pad(r_mtx, r_ln_cnt, r_col_cnt, 0, r_p_r, b_p_r).mtx_val;
        auto pad_elem_cnt = pad_ans * pad_ans;
        mtx_strassen left, right;
        left.ln_cnt = pad_ans; left.ln_cnt_orgn = pad_ans;
        left.col_cnt = pad_ans; left.col_cnt_orgn = pad_ans;
        right.ln_cnt = pad_ans; right.ln_cnt_orgn = pad_ans;
        right.col_cnt = pad_ans; right.col_cnt_orgn = pad_ans;
        left.ptr = mtx_s_l.get();
        right.ptr = mtx_s_r.get();
        auto ans_pad = mtx_strassen_mult(left, right, recursive_gate);
        auto ans = ans_pad.convert(true);
        return mtx_crop(ans, pad_ans, pad_ans, 0, r_p_r, b_p_l).mtx_val;
    }
    else return mtx_mult(l_mtx, r_mtx, l_ln_cnt, l_col_cnt, r_ln_cnt, r_col_cnt);
}

class matrix
{
protected:
    mtx_info info;
    uint64_t elem_cnt = 0;
    void _para_init(uint64_t ln_cnt, uint64_t col_cnt)
    {
        info.ln_cnt = ln_cnt;
        info.col_cnt = col_cnt;
        elem_cnt = ln_cnt * col_cnt;
    }
    void _init(uint64_t ln_cnt, uint64_t col_cnt)
    {
        _para_init(ln_cnt, col_cnt);
        bagrt::reset_ptr(info.mtx_val);
        info.mtx_val = mtx_init(ln_cnt, col_cnt);
    }
    void _para_reset()
    {
        elem_cnt = 0;
        info.col_cnt = 0;
        info.ln_cnt = 0;
    }
public:
    void reset()
    {
        _para_reset();
        bagrt::reset_ptr(info.mtx_val);
    }
    static matrix blank_matrix() { return matrix(); }
    bool is_matrix() { return (info.col_cnt && info.ln_cnt && info.mtx_val && elem_cnt && elem_cnt==info.ln_cnt*info.col_cnt); }
    MATRIX ptr() { return mtx_copy(info.mtx_val, info.ln_cnt, info.col_cnt); }
    matrix()
    {
        info.ln_cnt = 0;
        info.col_cnt = 0;
        info.mtx_val = nullptr;
    }
    matrix(uint64_t ln_cnt, uint64_t col_cnt, bool rand = false, double rand_boundry_first = 0, double rand_boundry_second = 0, double rand_acc = 1e-5)
    {
        if(ln_cnt && col_cnt)
        {
            reset();
            _para_init(ln_cnt, col_cnt);
            if(rand) info.mtx_val = mtx_init_rand(elem_cnt, rand_boundry_first, rand_boundry_second, rand_acc);
            else info.mtx_val = mtx_init(elem_cnt);
        }
    }
    matrix(MATRIX &&ptr_val, uint64_t ln_cnt, uint64_t col_cnt)
    {
        reset();
        _para_init(ln_cnt, col_cnt);
        info.mtx_val = std::move(ptr_val);
    }
    matrix(MATRIX &ptr_val, uint64_t ln_cnt, uint64_t col_cnt)
    {
        if(ln_cnt*col_cnt != elem_cnt)
        {
            reset();
            info.mtx_val = mtx_copy(ptr_val, ln_cnt*col_cnt);
        }
        else for(auto i=0; i<elem_cnt; ++i) info.mtx_val[i] = ptr_val[i];
        _para_init(ln_cnt, col_cnt);
    }
    matrix(double atom) 
    {  
        if(elem_cnt > 1) _init(1, 1);
        info.mtx_val[IDX_ZERO] = atom;
    }
    matrix(matrix &val) { value_copy(val); }
    matrix(matrix &&val) { value_move(std::move(val)); }
    matrix(const matrix &val)
    {
        _init(val.info.ln_cnt, val.info.col_cnt);
        for(auto i=0; i<elem_cnt; ++i) info.mtx_val[i] = val.info.mtx_val[i];
    }
    matrix(std::initializer_list<std::initializer_list<double>> _vect)
    {
        bagrt::net_list<double> elem_temp;
        auto _col_cnt = 0;
        bool asg_flag = true;
        for(auto ln_temp : _vect)
        {
            auto col_cnt_temp = 0;
            for(auto col_temp : ln_temp)
            {
                ++ col_cnt_temp;
                elem_temp.emplace_back(col_temp);
            }
            if(!_col_cnt) _col_cnt = col_cnt_temp;
            else if(_col_cnt != col_cnt_temp)
            {
                asg_flag = false;
                break;
            }
            else continue;
        }
        if(asg_flag)
        {
            auto _elem_cnt = elem_temp.size();
            if(_elem_cnt != elem_cnt || _col_cnt != info.col_cnt)  _init(_elem_cnt/_col_cnt, _col_cnt);
            for(auto i=0; i<elem_cnt; ++i) info.mtx_val[i] = elem_temp[i];
        }
    }
    matrix(std::initializer_list<std::initializer_list<matrix>> _vect)
    {
        int ln_cnt_hyper = _vect.size(), ln_cnt = 0, col_cnt_hyper = 0, col_cnt = 0;
        for(auto temp_ln : _vect)
        {
            ln_cnt += temp_ln.begin()->info.ln_cnt;
            for(auto temp_col : temp_ln)
            {
                if(col_cnt_hyper) break;
                else col_cnt += temp_col.info.col_cnt;
            }
            if(!col_cnt_hyper) col_cnt_hyper = temp_ln.size();
        }
        _init(ln_cnt, col_cnt);
        auto ln_iter = 0;
        for(auto temp_ln : _vect)
        {
            auto col_iter = 0;
            for(auto temp_col : temp_ln)
            {
                for(auto i=0; i<temp_col.elem_cnt; ++i)
                {
                    auto curr_child_dim = mtx_elem_pos(i, temp_col.info.col_cnt);
                    info.mtx_val[mtx_elem_pos(curr_child_dim.ln+ln_iter, curr_child_dim.col+col_iter, info.col_cnt)] = temp_col.info.mtx_val[i];
                }
                col_iter += temp_col.info.col_cnt;
                if(col_iter > info.col_cnt) return;
            }
            ln_iter += temp_ln.begin()->info.ln_cnt;
            if(ln_iter > info.ln_cnt) return;
        }
    }
    bool value_copy(matrix &val)
    {
        if(val.is_matrix())
        {
            if(val.elem_cnt != elem_cnt) _init(val.info.ln_cnt, val.info.col_cnt);
            else _para_init(val.info.ln_cnt, val.info.col_cnt);
            for(auto i=0; i<elem_cnt; ++i) info.mtx_val[i] = val.info.mtx_val[i];
            return true;
        }
        else return false;
    }
    bool value_move(matrix &&val)
    {
        if(val.is_matrix())
        {
            reset();
            _para_init(val.info.ln_cnt, val.info.col_cnt);
            info.mtx_val = std::move(val.info.mtx_val);
            val.reset();
            return true;
        }
        else return false;
    }
    void value_fill(double val) { for(auto i=0; i<elem_cnt; ++i) info.mtx_val[i] = val; }
    uint64_t get_ln_cnt() { return info.ln_cnt; }
    uint64_t get_col_cnt() { return info.col_cnt; }
    uint64_t get_elem_cnt() { return elem_cnt; }
    __declspec(property(get = get_ln_cnt)) uint64_t LN_CNT;
    __declspec(property(get = get_col_cnt)) uint64_t COL_CNT;
    __declspec(property(get = get_elem_cnt)) uint64_t ELEM_CNT;
    double determinant()
    {
        if(info.ln_cnt == info.col_cnt) return mtx_det(info.mtx_val, info.ln_cnt);
        else return NAN;
    }
    matrix inverser()
    {
        if(info.ln_cnt == info.col_cnt) return matrix(mtx_inverser(info.mtx_val, info.ln_cnt), info.ln_cnt, info.col_cnt);
        else return matrix();
    }
    matrix transposition()
    {
        if(is_matrix()) return matrix(mtx_transposition(info.mtx_val, info.ln_cnt, info.col_cnt), info.col_cnt, info.ln_cnt);
        else return matrix();
    }
    mtx_extm extremum(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilation = 0, uint64_t col_dilation = 0, bool max_flag = true) { return mtx_extm_val(info.mtx_val, from_ln, to_ln, from_col, to_col, info.ln_cnt, info.col_cnt, ln_dilation, col_dilation, max_flag); }
    double atom()
    {
        if(elem_cnt==1 && is_matrix()) return info.mtx_val[IDX_ZERO];
        else return NAN;
    }
    matrix child(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilation = 0, uint64_t col_dilation = 0) 
    {
        auto child_info = mtx_child_vec(info.mtx_val, from_ln, to_ln, from_col, to_col, info.ln_cnt, info.col_cnt, ln_dilation, col_dilation);
        return matrix(std::move(child_info.mtx_val), child_info.ln_cnt, child_info.col_cnt);
    }
    matrix rotate_rect(bool clockwise = true)
    {
        if(is_matrix()) return matrix(mtx_rotate_rect(info.mtx_val, info.ln_cnt, info.col_cnt, clockwise), info.col_cnt, info.ln_cnt);
        else return matrix();
    }
    matrix mirror_flip(bool is_vertical = true)
    {
        if(is_matrix()) return matrix(mtx_mirror_flip(info.mtx_val, info.ln_cnt, info.col_cnt, is_vertical), info.ln_cnt, info.col_cnt);
        else return matrix();
    }
    bool shape_valid(uint64_t ln_cnt, uint64_t col_cnt)
    {
        if(ln_cnt==info.ln_cnt && col_cnt==info.col_cnt) return true;
        else return false;
    }
    bool shape_valid(matrix &mtx_src) { return shape_valid(mtx_src.info.ln_cnt, mtx_src.info.col_cnt); }
    matrix reshape(uint64_t ln_cnt, uint64_t col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        if(this->elem_cnt==elem_cnt && is_matrix()) return matrix(info.mtx_val, ln_cnt, col_cnt);
        else return matrix();
    }
    matrix reshape(matrix &as_val) { return reshape(as_val.info.ln_cnt, as_val.info.col_cnt); }
    double elem_sum(uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilation = 0, uint64_t col_dilation = 0) { return mtx_sum(info.mtx_val, from_ln, to_ln, from_col, to_col, info.ln_cnt, info.col_cnt, ln_dilation, col_dilation); }
    double elem_sum() { return elem_sum(0, info.ln_cnt-1, 0, info.col_cnt-1); }
    matrix abs()
    {
        auto cpy_ptr = mtx_copy(info.mtx_val, elem_cnt);
        mtx_abs(cpy_ptr, elem_cnt);
        return matrix(std::move(cpy_ptr), info.ln_cnt, info.col_cnt);
    }
    matrix elem_cal_opt(matrix &r_val, uint64_t opt_idx)
    {
        if(shape_valid(r_val)) return matrix(mtx_elem_cal_opt(info.mtx_val, r_val.info.mtx_val, info.ln_cnt, info.col_cnt, opt_idx), info.ln_cnt, info.col_cnt);
        else return matrix();
    }
    matrix elem_cal_opt(double para, uint64_t opt_idx)
    {
        if(is_matrix()) return matrix(mtx_elem_cal_opt(info.mtx_val, para, info.ln_cnt, info.col_cnt, opt_idx), info.ln_cnt, info.col_cnt);
        else return matrix();
    }
    matrix broadcast_add(double val) { return matrix(mtx_broadcast_add(info.mtx_val, elem_cnt, val), info.ln_cnt, info.col_cnt); }
    matrix broadcast_subtract(double val, bool is_subtrahend = true) { return matrix(mtx_broadcast_subtract(info.mtx_val, elem_cnt, val, is_subtrahend), info.ln_cnt, info.col_cnt); }
    matrix pad(uint64_t ln_t = 0, uint64_t col_r = 0, uint64_t ln_b = 0, uint64_t col_l = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0)
    {
        if(ln_t || col_r || ln_b || col_l || ln_dist || col_dist)
        {
            auto pad_info = mtx_pad(info.mtx_val, info.ln_cnt, info.col_cnt, ln_t, col_r, ln_b, col_l, ln_dist, col_dist);
            return matrix(pad_info.mtx_val, pad_info.ln_cnt, pad_info.col_cnt);
        }
        else return *this;
    }
    matrix crop(uint64_t ln_t = 0, uint64_t col_r = 0, uint64_t ln_b = 0, uint64_t col_l = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0)
    {
        if(ln_t || col_r || ln_b || col_l || ln_dist || col_dist)
        {
            auto crop_info = mtx_crop(info.mtx_val, info.ln_cnt, info.col_cnt, ln_t, col_r, ln_b, col_l, ln_dist, col_dist);
            return matrix(crop_info.mtx_val, crop_info.ln_cnt, crop_info.col_cnt);
        }
        else return *this;  
    }
    matrix round_fit() 
    {
        if(is_matrix())
        {
            auto cpy_ptr = mtx_copy(info.mtx_val, elem_cnt);
            for(auto i=0; i<elem_cnt; ++i) cpy_ptr[i] = std::round(cpy_ptr[i]);
            return matrix(std::move(cpy_ptr), info.ln_cnt, info.col_cnt);
        }
        else return matrix();
    }
    matrix LU()
    {
        if(info.ln_cnt == info.col_cnt && is_matrix()) return matrix(mtx_LU(info.mtx_val, info.ln_cnt), info.ln_cnt, info.col_cnt);
        else return matrix();
    }
    matrix linear_eq(matrix &val_b, bool eq_idx = MATRIX_EQ_LU)
    {
        if(is_matrix() && info.ln_cnt == info.col_cnt && info.ln_cnt==val_b.info.ln_cnt && val_b.info.col_cnt==1)
            if(eq_idx) return matrix(mtx_equation(info.mtx_val, val_b.info.mtx_val, info.ln_cnt), info.ln_cnt, info.col_cnt);
            else return matrix(mtx_jacobi_iterate(info.mtx_val, val_b.info.mtx_val, info.ln_cnt), info.ln_cnt, info.col_cnt);
        else return matrix();
    }
    matrix swap_dir_elem(uint64_t l_idx, uint64_t r_idx, bool is_ln = true) { return matrix(mtx_swap_elem(info.mtx_val, l_idx, r_idx, info.ln_cnt, info.col_cnt, is_ln), info.ln_cnt, info.col_cnt); }
    matrix adjugate()
    {
        if(info.col_cnt == info.ln_cnt) return matrix(mtx_adjugate(info.mtx_val, info.ln_cnt), info.ln_cnt, info.col_cnt);
        else return blank_matrix();
    }
    matrix mult_strassen(matrix val) {}
    uint64_t rank() { return mtx_rank(info.mtx_val, info.ln_cnt, info.col_cnt); }
    double &pos_idx(uint64_t idx) { return info.mtx_val[idx]; }
    matrix operator+(matrix &val) { return matrix(mtx_add(info.mtx_val, val.info.mtx_val, info.ln_cnt, info.col_cnt), info.ln_cnt, info.col_cnt); }
    matrix operator-(matrix &val) { return matrix(mtx_add(info.mtx_val, val.info.mtx_val, info.ln_cnt, info.col_cnt, true), info.ln_cnt, info.col_cnt); }
    void operator+=(matrix &val) { *this = matrix(std::move(*this + val)); }
    void operator-=(matrix &val) { *this = matrix(std::move(*this - val)); }
    matrix operator*(matrix &val)
    {
        if(info.ln_cnt <= 1000) return matrix(mtx_mult(info.mtx_val, val.info.mtx_val, info.ln_cnt, info.col_cnt, val.info.ln_cnt, val.info.col_cnt), info.ln_cnt, val.info.col_cnt);
        else return matrix(mtx_strassen_mult(info.mtx_val, val.info.mtx_val, info.ln_cnt, info.col_cnt, val.info.ln_cnt, val.info.col_cnt), info.ln_cnt, val.info.col_cnt);
    }
    void operator*=(matrix &val) { *this = matrix(std::move(*this * val)); }
    matrix operator*(double val) { return matrix(mtx_mult(info.mtx_val, val, info.ln_cnt, info.col_cnt), info.ln_cnt, info.col_cnt); }
    void operator*=(double val) { *this = matrix(std::move(*this * val)); }
    friend matrix operator*(double val, matrix &r_val) { return r_val * val; }
    void operator=(matrix &val) { value_copy(val); }
    void operator=(matrix &&val) { value_move(std::move(val)); }
    void operator=(const matrix &val)
    {
        if(info.ln_cnt!=val.info.ln_cnt || info.col_cnt != val.info.col_cnt) _init(val.info.ln_cnt, val.info.col_cnt);
        for(auto i=0; i<elem_cnt; ++i) info.mtx_val[i] = val.info.mtx_val[i];
    }
    void operator=(std::initializer_list<std::initializer_list<double>> _vect) { *this = matrix(_vect); }
    bool operator==(matrix &val) { if(shape_valid(val)) return mtx_equal(info.mtx_val, val.info.mtx_val, elem_cnt); else return false; }
    bool operator!=(matrix &val) { return !(*this == val); }
    double *operator[](uint64_t ln) { return info.mtx_val.get() + ln * info.col_cnt; }
    friend std::ostream &operator<<(std::ostream &output, matrix &out_matrix)
    {
        if(out_matrix.is_matrix()) for(auto i=0; i<out_matrix.elem_cnt; ++i)
        {
            output << out_matrix.info.mtx_val[i];
            if(!((i+1)%out_matrix.info.col_cnt || (i+1)==out_matrix.elem_cnt)) output << std::endl;
            else output << '\t'; 
        }
        return output;
    }
    ~matrix() { reset(); }
};

MATRIX_END