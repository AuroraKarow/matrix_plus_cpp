MTX_BEGIN

bool mtx_elem_pos(uint64_t &ln, uint64_t &col, uint64_t idx, uint64_t curr_col_cnt, uint64_t ln_from, uint64_t col_from, uint64_t ln_dilate, uint64_t col_dilate)
{
    if(curr_col_cnt)
    {
        ln = 0; col = 0;
        if(ln_from || col_from  || ln_dilate || col_dilate)
        {
            ln = (idx / curr_col_cnt) * (ln_dilate + 1) + ln_from;
            col = (idx % curr_col_cnt) * (col_dilate + 1) + col_from;
        }
        else
        {
            ln = idx / curr_col_cnt;
            col = idx % curr_col_cnt;
        }
        return true;
    }
    else return false;
}
bool mtx_elem_pos(uint64_t &ln, uint64_t &col, uint64_t idx, uint64_t col_cnt) { return mtx_elem_pos(ln, col, idx, col_cnt, 0, 0, 0, 0); }
uint64_t mtx_elem_pos(uint64_t ln, uint64_t col, uint64_t orgn_col_cnt, uint64_t ln_from, uint64_t col_from, uint64_t ln_dilate, uint64_t col_dilate)
{
    auto curr_ln = 0, curr_col = 0;
    if(ln_from || col_from || ln_dilate || col_dilate)
    {
        curr_ln = ln_from + ln * (1 + ln_dilate),
        curr_col = col_from + col * (1 + col_dilate);
    }
    else
    {
        curr_ln = ln;
        curr_col = col;
    }
    return curr_ln * orgn_col_cnt + curr_col;
}
uint64_t mtx_elem_pos(uint64_t ln, uint64_t col, uint64_t col_cnt) { return mtx_elem_pos(ln, col, col_cnt, 0, 0, 0, 0); }

__mtx_callback bool mtx_copy(__mtx &ans, __mtx &src, uint64_t elem_cnt) { return _BAGRT ptr_copy(ans, src, elem_cnt); }
__mtx_callback bool mtx_copy(__mtx &ans, __mtx &src, uint64_t ln_cnt, uint64_t col_cnt) { return mtx_copy(ans, src, ln_cnt*col_cnt); }
__mtx_callback bool mtx_move(__mtx &ans, __mtx &&src) { return _BAGRT ptr_move(ans, std::move(src)); }

__mtx_callback bool mtx_fill(__mtx &src, __mtx_elem &&val, uint64_t elem_cnt)
{
    if(src)
    {
        std::fill_n(src, elem_cnt, val);
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_fill(__mtx &src, __mtx_elem &&val, uint64_t ln_cnt, uint64_t col_cnt) { return mtx_fill(src, std::move(val), ln_cnt*col_cnt); }

__mtx_callback bool mtx_init_E(__mtx &src, uint64_t dms)
{
    if(src && dms)
    {
        mtx_fill(src, __mtx_elem(0), dms*dms);
        for(auto i=0; i<dms; ++i) src[i*(dms+1)] = 1;
        return true;
    }
    return false;
}
__mtx_callback bool mtx_init_rand(__mtx &src, uint64_t elem_cnt, __mtx_elem &&boundary_first, __mtx_elem &&boundary_second, __mtx_elem &&acc)
{
    if(src)
    {
        for (int i = 0; i < elem_cnt; ++i)
            if(boundary_first == boundary_second) src[i] = (((double)_BAGRT rand_e() / (double)_BAGRT rand_e._Max) - 0.5) * 2.0;
            else src[i] = _BAGRT rand_num(boundary_first, boundary_second, false, acc);
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_print(__mtx &src, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(src)
    {
        auto elem_amt = ln_cnt * col_cnt;
        for (auto i = 0; i < elem_amt; ++i)
        {
            std::cout << src[i];
            if ((i + 1) % col_cnt) std::cout << '\t';
            else std::cout << std::endl;
        }
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_equal(__mtx &first, __mtx &second, uint64_t elem_cnt)
{
    for(auto i=0; i<elem_cnt; ++i) if(first[i] != second[i]) return false;
    return true;
}
__mtx_callback bool mtx_equal(__mtx &first, __mtx &second, uint64_t ln_cnt, uint64_t col_cnt) { return mtx_equal(first, second, ln_cnt*col_cnt); }

uint64_t mtx_pad_cnt(uint64_t prev_pad, uint64_t rear_pad, uint64_t dir_cnt, uint64_t dir_distance) { return prev_pad + rear_pad + dir_cnt + (dir_cnt - 1) * dir_distance; }
__mtx_callback bool mtx_pad(__mtx &ans, uint64_t ans_col_cnt, __mtx &src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_t, uint64_t col_r, uint64_t ln_b, uint64_t col_l, uint64_t ln_dist, uint64_t col_dist)
{
    if(src && ans && ans_col_cnt && ln_cnt && col_cnt)
    {
        for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j)
        {
            auto curr_pad_no = mtx_elem_pos(ln_t+i*(ln_dist+1), col_l+j*(col_dist+1), ans_col_cnt), curr_origin_no = mtx_elem_pos(i, j, col_cnt);
            ans[curr_pad_no] = src[curr_origin_no];
        }
        return true;
    }
    else return false;
}
uint64_t mtx_crop_cnt(uint64_t prev_crop, uint64_t rear_crop, uint64_t dir_cnt, uint64_t dir_distance) { return (dir_cnt - (prev_crop + rear_crop) + dir_distance) / (dir_distance + 1); }
__mtx_callback bool mtx_crop(__mtx &ans, uint64_t ans_ln_cnt, uint64_t ans_col_cnt, __mtx &src, uint64_t col_cnt, uint64_t ln_t, uint64_t col_r, uint64_t ln_b, uint64_t col_l, uint64_t ln_dist, uint64_t col_dist)
{
    if(src && ans && ans_col_cnt && ans_ln_cnt && col_cnt)
    {
        if(ln_t || col_r || ln_b || col_l || ln_dist || col_dist) for(auto i=0; i<ans_ln_cnt; ++i) for(auto j=0; j<ans_col_cnt; ++j)
        {
            auto curr_res_no = mtx_elem_pos(i, j, ans_col_cnt), curr_origin_no = mtx_elem_pos(ln_t+i*(ln_dist+1), col_l+j*(col_dist+1), col_cnt);
            ans[curr_res_no] = src[curr_origin_no];
        }
        else return mtx_copy(ans, src, ans_ln_cnt, ans_col_cnt);
        return true;
    }
    else return false;
}

__mtx_callback __mtx_elem mtx_extm_val(__mtx &src, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_cnt, uint64_t col_cnt, uint64_t *&ln_ls, uint64_t *&col_ls, uint64_t &pos_ls_len, bool max_flag, bool get_extm_pos, uint64_t ln_dilation, uint64_t col_dilation)
{
    MEM_RECYCLE(ln_ls); MEM_RECYCLE(col_ls);
    if(src && col_cnt && ln_cnt &&
        from_ln>=0 && to_ln>=from_ln && ln_cnt>to_ln &&
        from_col>=0 && to_col>=from_col && col_cnt>to_col)
    {
        pos_ls_len = 0;
        auto ans = src[mtx_elem_pos(from_ln, from_col, col_cnt)];
        for(auto i=from_ln; i<=to_ln; ++(i+=ln_dilation)) for(auto j=from_col; j<=to_col; ++(j+=col_dilation))
        {
            auto curr_no = mtx_elem_pos(i, j, col_cnt);
            auto curr_ln = i, curr_col = j;
            if((max_flag&&src[curr_no]>ans) || (!max_flag&&src[curr_no]<ans))
            {
                ans = src[curr_no];
                if(get_extm_pos)
                {
                    MEM_RECYCLE(ln_ls); MEM_RECYCLE(col_ls); pos_ls_len = 0;
                    _BAGRT ptr_insert(ln_ls, std::move(curr_ln), pos_ls_len, pos_ls_len);
                    _BAGRT ptr_insert(col_ls, std::move(curr_col), pos_ls_len, pos_ls_len);
                    ++ pos_ls_len;
                }
            }
            else if(src[curr_no]==ans && get_extm_pos)
            {
                _BAGRT ptr_insert(ln_ls, std::move(curr_ln), pos_ls_len, pos_ls_len);
                _BAGRT ptr_insert(col_ls, std::move(curr_col), pos_ls_len, pos_ls_len);
                ++ pos_ls_len;
            }
            else continue;
        }
        return ans;
    }
    else return 0;
}
__mtx_callback __mtx_elem mtx_sum(__mtx &src, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_cnt, uint64_t col_cnt, uint64_t ln_dilation, uint64_t col_dilation)
{
    auto sum = 0.0;
    if(src && col_cnt && ln_cnt &&
        from_ln>=0 && to_ln>=from_ln && ln_cnt>to_ln &&
        from_col>=0 && to_col>=from_col && col_cnt>to_col)
    {
        auto child_ln_cnt = mtx_child_vec_dir_cnt(from_ln, to_ln, ln_dilation),
            child_col_cnt = mtx_child_vec_dir_cnt(from_col, to_col, col_dilation);
        auto temp_elem_cnt = child_ln_cnt * child_col_cnt;
        MTX_INIT(__mtx_elem, temp, temp_elem_cnt);
        if(mtx_child_vec(temp, child_ln_cnt, child_col_cnt, src, from_ln, to_ln, from_col, to_col, ln_cnt, col_cnt, ln_dilation, col_dilation)) for(auto i=0; i<temp_elem_cnt; ++i) sum += temp[i];
        MTX_RESET(temp);        
    }
    return sum;
}

__mtx_callback bool mtx_add(__mtx &ans, __mtx &l_src, __mtx &r_src, uint64_t ln_cnt, uint64_t col_cnt, bool nega)
{
    if(ans && l_src && r_src && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        for(auto i=0; i<elem_cnt; ++i)
            if(nega) ans[i] = l_src[i] - r_src[i];
            else ans[i] = l_src[i] + r_src[i];
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_mult(__mtx &ans, __mtx &l_src, __mtx &r_src, uint64_t l_ln_cnt, uint64_t l_col_cnt, uint64_t r_ln_cnt, uint64_t r_col_cnt)
{
    if(ans && l_src && r_src && l_ln_cnt && l_col_cnt && r_ln_cnt && r_col_cnt && l_col_cnt==r_ln_cnt)
    {
        mtx_fill(ans, __mtx_elem(0), l_ln_cnt*r_col_cnt);
        for (int i=0; i<l_ln_cnt; ++i) for(auto j=0; j<l_col_cnt; ++j)
        {
            auto coe = l_src[mtx_elem_pos(i, j, l_col_cnt)];
            for(auto k=0; k<r_col_cnt; ++k)
            ans[mtx_elem_pos(i, k, r_col_cnt)] += coe * r_src[mtx_elem_pos(j, k, r_col_cnt)];
        }
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_mult(__mtx &ans, __mtx &src, __mtx_elem &&val, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(ans && src && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        for(auto i=0; i<elem_cnt; ++i) ans[i] = src[i] * val;
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_elem_cal_opt(__mtx &ans, __mtx &l_src, __mtx &r_src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t opt_idx)
{
    if(ans && l_src && r_src && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        for(auto i=0; i<elem_cnt; ++i)
            switch (opt_idx)
            {
            case MTX_ELEM_MULT:
                ans[i] = l_src[i] * r_src[i];
                break;
            case MTX_ELEM_DIV:
                if(r_src[i])
                {
                    ans[i] = l_src[i] / r_src[i];
                    break;
                }
                else return false;
            default: return false;
            }
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_elem_cal_opt(__mtx &ans, __mtx &src, __mtx_elem &&para, uint64_t ln_cnt, uint64_t col_cnt, uint64_t opt_idx)
{
    if(ans && src && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        if(para != 1)
        {
            for(auto i=0; i<elem_cnt; ++i) switch (opt_idx)
            {
            case MTX_ELEM_POW:
                ans[i] = _BAGRT __power_(src[i], para);
                break;
            case MTX_ELEM_DIV:
                if(para != 0)
                {
                    ans[i] = src[i] / para;
                    break;
                }
                else return false;
            default: return false;
            }
        }
        else mtx_copy(ans, src, ln_cnt, col_cnt);
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_broadcast_add(__mtx &ans, __mtx &src, uint64_t elem_cnt, __mtx_elem &&val, bool nega)
{
    if(ans && src && elem_cnt)
    {
        for(auto i=0; i<elem_cnt; ++i)
            if(nega) ans[i] -= val;
            else ans[i] += val;
        return true;
    }
    else return false;
}

__mtx_callback bool mtx_abs(__mtx &src, uint64_t elem_cnt)
{
    if(src)
    {
        for(auto i=0; i<elem_cnt; ++i) if(src[i] < 0) src[i] *= (-1);
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_transposition(__mtx &ans, __mtx &src, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(ans && src && ln_cnt && col_cnt)
    {
        auto elem_cnt = ln_cnt * col_cnt;
        for (auto i = 0; i < ln_cnt; ++i) for (auto j = 0; j < col_cnt; ++j)
            ans[mtx_elem_pos(j, i, ln_cnt)] = src[mtx_elem_pos(i, j, col_cnt)];
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_swap_elem(__mtx &ans, __mtx &src, uint64_t l_pos, uint64_t r_pos, uint64_t ln_cnt, uint64_t col_cnt, bool is_ln)
{
    if(ans && src && ln_cnt && col_cnt)
    {
        mtx_copy(ans, src, ln_cnt, col_cnt);
        if(is_ln) for(auto i=0; i<col_cnt; ++i) std::swap(ans[mtx_elem_pos(l_pos, i, col_cnt)], ans[mtx_elem_pos(r_pos, i, col_cnt)]);
        else for(auto i=0; i<ln_cnt; ++i) std::swap(ans[mtx_elem_pos(i, l_pos, col_cnt)], ans[mtx_elem_pos(i, r_pos, col_cnt)]);
        return true;
    }
    else return false;        
}
__mtx_callback __mtx_elem mtx_det(__mtx &src, uint64_t dms)
{
    if(src && dms)
    {
        if (dms == 1)
        return src[IDX_ZERO];
        MTX_INIT(__mtx_elem, ac, (dms-1)*(dms-1));
        int mov = 0;
        __mtx_elem sum = 0.0;
        for (int arow = 0; arow < dms; arow++)
        {
            for (int brow = 0; brow < dms - 1; brow++)
            {
                mov = arow > brow ? 0 : 1;
                for (int j = 0; j < dms - 1; ++j)
                    ac[brow * (dms - 1) + j] = src[(brow + mov) * dms + j + 1];
            }
            int flag = (arow % 2 == 0 ? 1 : -1);
            sum += flag * src[arow * dms] * mtx_det(ac, dms - 1);
        }
        MTX_RESET(ac);
        return sum;
    }
    else return NAN;
}
__mtx_callback bool mtx_cofactor(__mtx &ans, __mtx &src, uint64_t ln, uint64_t col, uint64_t dms)
{
    if(ans && src && dms)
    {
        auto ans_dms = dms - 1;
        auto elem_cnt = dms * dms;
        for(auto i=0; i<elem_cnt; ++i)
        {
            auto curr_ln = 0ull, curr_col = 0ull;
            mtx_elem_pos(curr_ln, curr_col, i, dms);
            if(curr_ln!=ln && curr_col!=col)
            {
                auto curr_ans_ln = curr_ln, curr_ans_col = curr_col;
                if(curr_ln > ln) -- curr_ans_ln;
                if(curr_col > col) -- curr_ans_col;
                ans[mtx_elem_pos(curr_ans_ln, curr_ans_col, ans_dms)] = src[i];
            }
        }
        return true;
    }
    else return false;
}
__mtx_callback __mtx_elem mtx_algebraic_cofactor(__mtx &src, uint64_t ln, uint64_t col, uint64_t dms)
{
    MTX_INIT(__mtx_elem, cofactor, (dms-1)*(dms-1));
    __mtx_elem ans = 0.0;
    if(mtx_cofactor(cofactor, src, ln, col, dms)) ans = mtx_det(cofactor, dms-1);
    MTX_RESET(cofactor);
    return ans;
}
__mtx_callback bool mtx_adjugate(__mtx &ans, __mtx &src, uint64_t dms)
{
    if(ans && src && dms)
    {
        auto elem_cnt = dms * dms;
        for(auto i=0; i<elem_cnt; ++i)
        {
            auto curr_ln = 0ull, curr_col = 0ull;
            if(!mtx_elem_pos(curr_ln, curr_col, i, dms)) return false;
            auto coe = 1;
            if((curr_ln+curr_col) % 2) coe = -1;
            ans[mtx_elem_pos(curr_col, curr_ln, dms)] = coe * mtx_algebraic_cofactor(src, curr_ln, curr_col, dms);
        }
        return true;
    }
    return true;
}
__mtx_callback bool mtx_inverser(__mtx &ans, __mtx &src, uint64_t dms)
{
    if(ans && src && dms)
    {
        auto elem_cnt = dms * dms;
        auto val_det = mtx_det(src, dms);
        MTX_INIT(__mtx_elem, val_adj, elem_cnt);
        mtx_adjugate(val_adj, src, dms);
        for(auto i=0; i<elem_cnt; ++i) ans[i] = val_adj[i] / val_det;
        MTX_RESET(val_adj);
        return true;
    }
    else return false;
}
__mtx_callback __mtx_elem mtx_max_eigenvalue(__mtx &src, __mtx &w, uint64_t dms, __mtx_elem &&acc, __mtx_elem &&init_elem)
{
    if(src && dms && w)
    {
        __mtx_elem lambda = 0.0, lambda_temp = 0.0;
        MTX_INIT(__mtx_elem, x, dms);
        uint64_t *ptr_temp_true = nullptr, *ptr_temp_false = nullptr, ptr_temp_len = 0;
        mtx_fill(x, std::move(init_elem), dms);
        do
        {
            lambda = lambda_temp;
            MTX_INIT(__mtx_elem, temp, dms); MTX_INIT(__mtx_elem, temp_abs, dms);
            mtx_copy(temp, x, dms);
            mtx_copy(temp_abs, temp, dms);
            mtx_abs(temp_abs, dms);
            __mtx_elem max = mtx_extm_val(temp_abs, 0, dms-1, 0, 0, dms, 1, ptr_temp_true, ptr_temp_false, ptr_temp_len, true, false);
            for (int i=0; i<dms; ++i) x[i] /= max;
            for (int i=0; i<dms; ++i) w[i] = x[i];
            __mtx_elem sum = 0.0;
            for (int i=0; i<dms; ++i) sum += w[i];
            for (int i=0; i<dms; ++i) w[i] /= sum;
            mtx_mult(x, src, temp, dms, dms, dms, 1);
            mtx_copy(temp, x, dms);
            lambda_temp = mtx_extm_val(temp, 0, dms-1, 0, 0, dms, 1, ptr_temp_true, ptr_temp_false, ptr_temp_len, true, false);
            MTX_RESET(temp); MTX_RESET(temp_abs);
        } while (std::abs(lambda_temp - lambda) > acc);
        MTX_RESET(x);
        return lambda;
    }
    else return NAN;
}
__mtx_callback uint64_t mtx_rank(__mtx &src, uint64_t ln_cnt, uint64_t col_cnt)
{
    if(ln_cnt && col_cnt && src)
    {
        MTX_PTR(__mtx_elem, ans);
        if(ln_cnt > col_cnt)
        {
            std::swap(ln_cnt, col_cnt);
            MTX_ALLOC(__mtx_elem, ans, col_cnt*ln_cnt);
            mtx_transposition(ans, src, ln_cnt, col_cnt);
        }
        else
        {
            MTX_ALLOC(__mtx_elem, ans, ln_cnt*col_cnt);
            mtx_copy(ans, src, ln_cnt, col_cnt);
        }
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
                        if(j!=i)
                        {
                            MTX_INIT(__mtx_elem, temp, ln_cnt*col_cnt);
                            mtx_copy(temp, ans, ln_cnt, col_cnt);
                            mtx_swap_elem(ans, temp, i, j, ln_cnt, col_cnt);
                            MTX_RESET(temp);
                        }
                        ++ rank_val;
                        elim_flag = false;
                    }
                    else for(auto k=i; k<col_cnt; ++k) ans[mtx_elem_pos(j, k, col_cnt)] -= ans[mtx_elem_pos(i, k, col_cnt)] * elim_val;
                }
            }
        }
        MTX_RESET(ans);
        return rank_val;
    }
    else return NAN;
}

__mtx_callback bool mtx_LU(__mtx &ans, __mtx &src, uint64_t dms)
{
    if(ans && src && dms)
    {
        auto elem_cnt = dms * dms;
        mtx_copy(ans, src, elem_cnt);
        MTX_INIT(__mtx_elem, e, elem_cnt);
        MTX_INIT(__mtx_elem, l, elem_cnt);
        mtx_init_E(l, dms);
        MTX_INIT(__mtx_elem, temp, elem_cnt);
        for (int i = 0; i < dms - 1; ++i)
        {
            mtx_init_E(e, dms);
            for (int j = dms - 1; j > i; --j) e[mtx_elem_pos(j, i, dms)] = (-1.0) * ans[mtx_elem_pos(j, i, dms)] / ans[mtx_elem_pos(i, i, dms)];
            for (int j = dms - 1; j > i; --j) l[mtx_elem_pos(j, i, dms)] = (-1.0) * e[mtx_elem_pos(j, i, dms)];
            for (int j = dms - 1; j > i; --j) ans[mtx_elem_pos(j, i, dms)] = (-1.0) * e[mtx_elem_pos(j, i, dms)];
            mtx_copy(temp, ans, elem_cnt);
            mtx_mult(ans, e, temp, dms, dms, dms, dms);
        }
        // get the inversion matrix
        for (int i = 0; i < dms; ++i) for (int j = 0; j < i; ++j) ans[mtx_elem_pos(i, j, dms)] = l[mtx_elem_pos(i, j, dms)];
        MTX_RESET(e); MTX_RESET(l); MTX_RESET(temp);
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_equation(__mtx &ans, __mtx &coefficient, __mtx &b, uint64_t dms)
{
    if(ans && coefficient && b && dms)
    {
        auto elem_cnt = dms * dms;
        MTX_INIT(__mtx_elem, lu, elem_cnt);
        mtx_LU(lu, coefficient, dms);
        MTX_INIT(__mtx_elem, e, elem_cnt); MTX_INIT(__mtx_elem, l, elem_cnt); MTX_INIT(__mtx_elem, u, elem_cnt);
        mtx_init_E(e, dms); mtx_init_E(l, dms); mtx_init_E(u, dms);
        MTX_INIT(__mtx_elem, y, dms);
        for (int i = 0; i < dms; ++i) for (int j = 0; j < dms; ++j)
            if (i <= j) u[mtx_elem_pos(i, j, dms)] = lu[mtx_elem_pos(i, j, dms)];
            else l[mtx_elem_pos(i, j, dms)] = lu[mtx_elem_pos(i, j, dms)];
        for (int i = 0; i < dms; ++i) if (i)
        {
            __mtx_elem temp = 0.0;
            for (int j=0; j < i; ++j) temp += y[j] * l[mtx_elem_pos(i, j, dms)] * 1.0;
            y[i] = b[i] - temp;
        }
        else y[i] = b[i];
        for (int i=dms; i>0; --i) if (i - dms)
        {
            __mtx_elem temp = 0.0;
            for (int n=i-1; n<dms-1; ++n) temp += u[mtx_elem_pos(i-1, n+1, dms)] * ans[n+1];
            ans[i-1] = (y[i-1] - temp) / u[mtx_elem_pos(i-1, i-1, dms)];
        }
        else ans[dms-1] = y[dms-1] / u[mtx_elem_pos(dms-1, dms-1, dms)];
        MTX_RESET(lu); MTX_RESET(e); MTX_RESET(l); MTX_RESET(u); MTX_RESET(y);
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_jacobi_iterate(__mtx &ans, __mtx &coefficient, __mtx &b, uint64_t dms, __mtx_elem &&acc, __mtx_elem &&init_elem)
{
    if(ans && coefficient && b && dms && init_elem>0)
    {
        MTX_INIT(__mtx_elem, temp, dms);
        mtx_fill(temp, std::move(init_elem), dms);
        bool flag = false;
        do
        {
            flag = false;
            for (int i=0; i<dms; ++i) ans[i] = temp[i];
            for (int i=0; i<dms; ++i)
            {
                __mtx_elem sum = 0.0;
                for (int j=0; j<dms; ++j) if (i != j) sum += coefficient[mtx_elem_pos(i, j, dms)] * ans[j];
                temp[i] = (b[i] - sum) / coefficient[mtx_elem_pos(i, i, dms)];
            }
            for (int i=0; i<dms; ++i) if (std::abs(ans[i] - temp[i]) >= acc)
            {
                flag = true;
                break;
            }
        } while (flag);
        MTX_RESET(temp);
        return true;
    }
    else return false;
}

__mtx_callback bool mtx_rotate_rect(__mtx &ans, __mtx &src, int ln_cnt, int col_cnt, bool clock_wise)
{
    if(ans && src && ln_cnt && col_cnt)
    {
        for(auto i=0; i<ln_cnt; ++i) for(auto j=0; j<col_cnt; ++j)
            if(clock_wise)  ans[mtx_elem_pos(j, ln_cnt-i-1, ln_cnt)] = src[mtx_elem_pos(i, j, col_cnt)];
            else ans[mtx_elem_pos(col_cnt-j-1, i, ln_cnt)] = src[mtx_elem_pos(i, j, col_cnt)];
        return true;
    }
    else return false;
}
__mtx_callback bool mtx_mirror_flip(__mtx &ans, __mtx &src, uint64_t ln_cnt, uint64_t col_cnt, bool symmetry_vertical)
{
    if(ans && src && ln_cnt && col_cnt)
    {
        
        for(auto i=0; i<ln_cnt; ++i) for(int j=0; j<col_cnt; ++j) if(symmetry_vertical)
        {
            auto src_pos = col_cnt-1-j;
            ans[mtx_elem_pos(i, j, col_cnt)] = src[mtx_elem_pos(i, src_pos, col_cnt)];
        }
        else
        {
            auto src_pos = ln_cnt-1-i;
            ans[mtx_elem_pos(i, j, col_cnt)] = src[mtx_elem_pos(src_pos, j, col_cnt)];
        }
        return true;
    }
    else return false;
}

MTX_END