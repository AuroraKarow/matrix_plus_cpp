MATRIX_BEGIN

struct pos { uint64_t ln = 0, col = 0; };

pos elem_pos(uint64_t idx, uint64_t curr_col_cnt, uint64_t ln_from, uint64_t col_from, uint64_t ln_dilate, uint64_t col_dilate) {
    pos ans{};
    if (curr_col_cnt == 0) return ans;
    if (ln_from || col_from  || ln_dilate || col_dilate) {
        ans.ln  = (idx / curr_col_cnt) * (ln_dilate + 1) + ln_from;
        ans.col = (idx % curr_col_cnt) * (col_dilate + 1) + col_from;
    } else {
        ans.ln  = idx / curr_col_cnt;
        ans.col = idx % curr_col_cnt;
    }
    return ans;
}
pos elem_pos(uint64_t idx, uint64_t col_cnt) { return elem_pos(idx, col_cnt, 0, 0, 0, 0); }
uint64_t elem_pos(uint64_t ln, uint64_t col, uint64_t orgn_col_cnt, uint64_t ln_from, uint64_t col_from, uint64_t ln_dilate, uint64_t col_dilate) {
    auto curr_ln  = 0,
         curr_col = 0;
    if (ln_from || col_from || ln_dilate || col_dilate) {
        curr_ln  = ln_from + ln * (1 + ln_dilate),
        curr_col = col_from + col * (1 + col_dilate);
    } else {
        curr_ln  = ln;
        curr_col = col;
    }
    return curr_ln * orgn_col_cnt + curr_col;
}
uint64_t elem_pos(uint64_t ln, uint64_t col, uint64_t col_cnt) { return elem_pos(ln, col, col_cnt, 0, 0, 0, 0); }

callback_matrix matrix_ptr init(uint64_t elem_cnt) {
    if (elem_cnt) {
        auto ans = new matrix_elem_t[elem_cnt];
        std::fill_n(ans, elem_cnt, 0);
        return ans;
    }
    else return nullptr;
}
callback_matrix_n auto init(std::initializer_list<std::initializer_list<matrix_elem_t>> src, uint64_t &ln_cnt, uint64_t &col_cnt) {
    ln_cnt   = src.size(),
    col_cnt  = src.begin()->size();
    auto cnt = 0ull;
    if constexpr (std::is_integral_v<matrix_elem_t>) {
        auto ans = init<long double>(ln_cnt * col_cnt);
        for (auto ln_temp : src) for (auto temp : ln_temp) *(ans + cnt++) = std::move(temp);
        return ans;
    } else {
        auto ans = init<matrix_elem_t>(ln_cnt * col_cnt);
        for (auto ln_temp : src) for (auto temp : ln_temp) *(ans + cnt++) = std::move(temp);
        return ans;
    }
}

callback_matrix matrix_ptr init_identity(uint64_t dim_cnt) {
    if (dim_cnt == 0) return nullptr;
    auto ans = init<matrix_elem_t>(dim_cnt * dim_cnt);
    for(auto i = 0ull; i < dim_cnt; ++i) *(ans + i * (dim_cnt + 1)) = 1;
    return ans;
}

callback_matrix matrix_ptr init_rand(uint64_t elem_cnt, const matrix_elem_t &fst_rng = 0, const matrix_elem_t &snd_rng = 0, uint64_t acc = 8) {
    if (elem_cnt == 0) return nullptr;
    auto ans = init<matrix_elem_t>(elem_cnt);
    for (auto i = 0ull; i < elem_cnt; ++i)
        if constexpr (std::is_same_v<matrix_elem_t, net_decimal>) *(ans + i) = num_rand(fst_rng.number_format, snd_rng.number_format, acc);
        else *(ans + i) = num_rand(fst_rng, snd_rng, acc);
    return ans;
}

callback_matrices void recycle(matrices_ptr &...val) { ptr_reset(val...); }

callback_matrix matrix_ptr copy(const matrix_ptr src, uint64_t elem_cnt) {
    if (src && elem_cnt) return ptr_copy(src, elem_cnt);
    else return nullptr;
}
callback_matrix bool copy(matrix_ptr &dest, uint64_t dest_elem_cnt, const matrix_ptr src, uint64_t elem_cnt) {
    if (dest_elem_cnt != elem_cnt) ptr_alter(dest, dest_elem_cnt, elem_cnt, false);
    return ptr_copy(dest, src, elem_cnt);
}

callback_matrix void move(matrix_ptr &dest, matrix_ptr &&src) { ptr_move(dest, std::move(src)); }

callback_matrix void fill(matrix_ptr &dest, uint64_t elem_cnt, const matrix_elem_t &src = 0) { std::fill_n(dest, elem_cnt, src); }

callback_matrix void print(matrix_ptr &src, uint64_t ln_cnt, uint64_t col_cnt) {
    auto elem_cnt = ln_cnt * col_cnt;
    if (elem_cnt == 0) std::cout << "[Null]" << std::endl;
    for (auto i = 0ull; i < elem_cnt; ++i) {
        std::cout << *(src + i);
        if ((i + 1) % col_cnt) std::cout << '\t';
        else std::cout << '\n';
    }
}

callback_matrix bool compare(const matrix_ptr fst, const matrix_ptr snd, uint64_t elem_cnt) {
    if (elem_cnt && fst && snd) {
        for (auto i = 0ull; i < elem_cnt; ++i) if (*(fst + i) != *(snd + i)) return false;
        return true; 
    } else return false;
}

callback_matrix matrix_ptr absolute(const matrix_ptr src, uint64_t elem_cnt) { 
    matrix_ptr ans = nullptr;
    if (src && elem_cnt) {
        ans = init<matrix_elem_t>(elem_cnt);
        for (auto i = 0ull; i < elem_cnt; ++i)
            if (*(src + i) < 0) *(ans + i) = (-1) * (*(src + i));
            else *(ans + i) = *(src + i);
    }
    return ans;
}

uint64_t pad_res_dir_cnt(uint64_t prev_pad, uint64_t rear_pad, uint64_t dir_cnt, uint64_t dir_dist) { return prev_pad + rear_pad + dir_cnt + (dir_cnt - 1) * dir_dist; }

uint64_t crop_res_dir_cnt(uint64_t prev_crop, uint64_t rear_crop, uint64_t dir_cnt, uint64_t dir_dist) { return (dir_cnt - (prev_crop + rear_crop) + dir_dist) / (dir_dist + 1); }

callback_matrix matrix_ptr pad(uint64_t &pad_ln_cnt, uint64_t &pad_col_cnt, const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t top = 0, uint64_t right = 0, uint64_t bottom = 0, uint64_t left = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0) {
    pad_ln_cnt  = pad_res_dir_cnt(top, bottom, ln_cnt, ln_dist);
    pad_col_cnt = pad_res_dir_cnt(left, right, col_cnt, col_dist);
    if (src && ln_cnt && col_cnt && !(pad_ln_cnt == ln_cnt && pad_col_cnt == col_cnt)) {
        auto ans = init<matrix_elem_t>(pad_ln_cnt * pad_col_cnt);
        for (auto i = 0ull; i < ln_cnt; ++i) for (auto j = 0ull; j < col_cnt; ++j) {
            auto curr_origin_no = elem_pos(i, j, col_cnt),
                 curr_pad_no    = elem_pos(top + i * (ln_dist + 1), left + j * (col_dist + 1), pad_col_cnt);
            *(ans + curr_pad_no) = *(src + curr_origin_no);
        }
        return ans;
    } else return copy(src, ln_cnt * col_cnt);
}

callback_matrix matrix_ptr crop(uint64_t &crop_ln_cnt, uint64_t &crop_col_cnt, const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t top = 0, uint64_t right = 0, uint64_t bottom = 0, uint64_t left = 0, uint64_t ln_dist = 0, uint64_t col_dist = 0) {
    crop_ln_cnt  = crop_res_dir_cnt(top, bottom, ln_cnt, ln_dist);
    crop_col_cnt = crop_res_dir_cnt(left, right, col_cnt, col_dist);
    if (src && ln_cnt && col_cnt && !(crop_ln_cnt == ln_cnt && crop_col_cnt == col_cnt)) {
        auto ans = init<matrix_elem_t>(crop_ln_cnt * crop_col_cnt);
        for (auto i = 0ull; i < crop_ln_cnt; ++i) for (auto j = 0ull; j < crop_col_cnt; ++j) {
            auto curr_res_no     = elem_pos(i, j, crop_col_cnt),
                 curr_origin_no  = elem_pos(top + i * (ln_dist + 1), left + j * (col_dist + 1), col_cnt);
            *(ans + curr_res_no) = *(src + curr_origin_no);
        }
        return ans;
    } else return copy(src, ln_cnt * col_cnt);
}

callback_matrix matrix_elem_t extremum(net_set<pos> &elem_pos_set, bool get_max, const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilate = 0, uint64_t col_dilate = 0) {
    matrix_elem_t ans {};
    if(src && col_cnt && ln_cnt && from_ln >= 0 && to_ln >= from_ln && ln_cnt > to_ln && from_col >= 0 && to_col >= from_col && col_cnt > to_col) {
        net_sequence<pos> elem_pos_seq;
        ans = *(src + elem_pos(from_ln, from_col, col_cnt));
        for (auto i = from_ln; i <= to_ln; ++(i += ln_dilate)) for (auto j = from_col; j <= to_col; ++(j += col_dilate)) {
            auto curr_no  = elem_pos(i, j, col_cnt);
            pos  curr_pos {i, j};
            if ((get_max && src[curr_no] > ans) || (!get_max && src[curr_no] < ans)) {
                ans = *(src + curr_no);
                elem_pos_seq.clear();
                elem_pos_seq.emplace_back(curr_pos);
            } else if (*(src + curr_no) == ans) elem_pos_seq.emplace_back(curr_pos);
        }
        elem_pos_seq.shrink();
        elem_pos_set = std::move(elem_pos_seq);
    }
    return ans;
}

callback_matrix matrix_ptr sub(uint64_t &sub_ln_cnt, uint64_t &sub_col_cnt, const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilate = 0, uint64_t col_dilate = 0) {
    if (src && col_cnt && ln_cnt && from_ln >= 0 && to_ln >= from_ln && ln_cnt > to_ln && from_col >= 0 && to_col >= from_col && col_cnt > to_col) {
        sub_ln_cnt  = (to_ln - from_ln) / (1 + ln_dilate) + 1;
        sub_col_cnt = (to_col - from_col) / (1 + col_dilate) + 1;
        if (sub_ln_cnt == ln_cnt && sub_col_cnt == col_cnt) return copy(src, ln_cnt * col_cnt);
        auto ans = init<matrix_elem_t>(sub_ln_cnt * sub_col_cnt);
        auto cnt = 0ull;
        for (auto i = from_ln; i <= to_ln; ++(i += ln_dilate)) for (auto j = from_col; j <= to_col; ++(j += col_dilate)) *(ans + cnt++) = *(src + elem_pos(i, j, col_cnt));
        return ans;
    }
    return nullptr;
}

callback_matrix matrix_elem_t sum(const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t from_ln, uint64_t to_ln, uint64_t from_col, uint64_t to_col, uint64_t ln_dilate = 0, uint64_t col_dilate = 0) {
    auto sub_src = sub(ln_cnt, col_cnt, src, ln_cnt, col_cnt, from_ln, to_ln, from_col, to_col, ln_dilate, col_dilate);
    if (sub_src == nullptr) return 0;
    matrix_elem_t ans      = 0;
    auto          elem_cnt = ln_cnt * col_cnt;
    for (auto i = 0ull; i < elem_cnt; ++i) ans += *(sub_src + i);
    recycle(sub_src);
    return ans;
}

callback_matrix matrix_ptr add(const matrix_ptr fst, const matrix_ptr snd, uint64_t elem_cnt, bool subtract = false) {
    if (elem_cnt && fst && snd) {
        auto ans = init<matrix_elem_t>(elem_cnt);
        for (auto i = 0ull; i < elem_cnt; ++i)
            if (subtract) *(ans + i) = *(fst + i) - *(snd + i);
            else *(ans + i) = *(fst + i) + *(snd + i);
        return ans;
    }
    return nullptr;
}

/* ans_ln_cnt  = fst_ln_cnt
 * ans_col_cnt = snd_col_cnt
 */
callback_matrix matrix_ptr mult(const matrix_ptr fst, uint64_t fst_ln_cnt, uint64_t fst_col_cnt, const matrix_ptr snd, uint64_t snd_col_cnt) {
    if (fst && snd && fst_ln_cnt && fst_col_cnt && snd_col_cnt) {
        auto ans = init<matrix_elem_t>(fst_ln_cnt * snd_col_cnt);
        for (auto i = 0ull; i < fst_ln_cnt; ++i) for (auto j = 0ull; j < fst_col_cnt; ++j) {
            auto coe = *(fst + elem_pos(i, j, fst_col_cnt));
            for (auto k = 0ull; k < snd_col_cnt; ++k) *(ans + elem_pos(i, k, snd_col_cnt)) += coe * (*(snd + elem_pos(j, k, snd_col_cnt)));
        }
        return ans;
    }
    return nullptr;
}
callback_matrix matrix_ptr mult(const matrix_ptr val, uint64_t elem_cnt, const matrix_elem_t &coe) {
    if (coe == 1) return copy(val, elem_cnt);
    auto ans = init<matrix_elem_t>(elem_cnt);
    if (coe == 0) return ans;
    for (auto i = 0ull; i < elem_cnt; ++i) *(ans + i) = *(val + i) * coe;
    return ans;
}

callback_matrix matrix_ptr elem_operate(const matrix_ptr val, uint64_t elem_cnt, const matrix_elem_t &para, uint64_t operation, bool para_fst = false, long double epsilon = 1e-8) {
    matrix_ptr ans = nullptr;
    if (val && elem_cnt) {
        ans = init<matrix_elem_t>(elem_cnt);
        for (auto i = 0ull; i < elem_cnt; ++i) {
            auto elem      = *(val + i),
                 curr_val  = para;
            if (para_fst) std::swap(elem, curr_val);
            switch (operation) {
            case MATRIX_ELEM_DIV:
                if (curr_val == 0) curr_val = epsilon;
                *(ans + i) = elem / curr_val;
                break;
            case MATRIX_ELEM_POW:
                *(ans + i) = std::pow(elem, curr_val);
                break;
            default: recycle(ans); break;
            }
            if (ans == nullptr) break;
        }
    }
    return ans;
}
callback_matrix matrix_ptr elem_operate(const matrix_ptr fst, const matrix_ptr snd, uint64_t elem_cnt, uint64_t operation, bool elem_swap = false, long double epsilon = 1e-8) {
    matrix_ptr ans = nullptr;
    if (fst && snd && elem_cnt) {
        ans = init<matrix_elem_t>(elem_cnt);
        for (auto i = 0ull; i < elem_cnt; ++i) {
            switch (operation) {
            case MATRIX_ELEM_DIV: {
                auto fst_elem = *(fst + i),
                     snd_elem = *(snd + i);
                if (elem_swap) std::swap(fst_elem, snd_elem);
                if (snd_elem == 0) snd_elem = epsilon;
                *(ans + i) = fst_elem / snd_elem;
                break;
            }
            case MATRIX_ELEM_MULT: *(ans + i) = *(fst + i) * (*(snd + i)); break;
            default: recycle(ans); break;
            }
            if (ans == nullptr) break;
        }
    }
    return ans;
}

callback_matrix matrix_ptr broadcast_add(const matrix_ptr val, uint64_t elem_cnt, const matrix_elem_t &para, bool subtract = false, bool para_fst = false) {
    matrix_ptr ans = nullptr;
    if (val && elem_cnt) {
        ans = init<matrix_elem_t>(elem_cnt);
        for (auto i = 0ull; i < elem_cnt; ++i)
            if (subtract) {
                auto subtrahend = *(val + i),
                     minuend    = para;
                if (para_fst) std::swap(subtrahend, minuend);
                *(ans + i) = subtrahend - minuend;
            }
            else *(ans + i) = *(val + i) + para;
    }
    return ans;
}

callback_matrix matrix_ptr transposition(const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt) {
    matrix_ptr ans = nullptr;
    if (src && ln_cnt && col_cnt) {
        ans = init<matrix_elem_t>(ln_cnt * col_cnt);
        for (auto i = 0ull; i < ln_cnt; ++i) for (auto j = 0ull; j < col_cnt; ++j) *(ans + elem_pos(j, i, ln_cnt)) = *(src + elem_pos(i, j, col_cnt));
    }
    return ans;
}

callback_matrix matrix_ptr ln_col_swap(const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt, uint64_t fst_ln_col, uint64_t snd_ln_col, bool ln_dir = true) {
    matrix_ptr ans = nullptr;
    if (snd_ln_col < fst_ln_col) std::swap(snd_ln_col, fst_ln_col);
    if(src && ln_cnt && col_cnt && ((ln_dir && snd_ln_col < ln_cnt) || (!ln_dir && snd_ln_col < col_cnt))) {
        ans = copy(src, ln_cnt * col_cnt);
        if (fst_ln_col == snd_ln_col) return ans;
        if (ln_dir) for (auto i = 0ull; i < col_cnt; ++i) std::swap(*(ans + elem_pos(fst_ln_col, i, col_cnt)), *(ans + elem_pos(snd_ln_col, i, col_cnt)));
        else for (auto i = 0ull; i < ln_cnt; ++i) std::swap(*(ans + elem_pos(i, fst_ln_col, col_cnt)), *(ans + elem_pos(i, snd_ln_col, col_cnt)));
    }
    return ans;
}

callback_matrix uint64_t rank(const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt) {
    if (ln_cnt && col_cnt && src) {
        auto elem_cnt = ln_cnt * col_cnt;
        auto ans      = new long double[elem_cnt];
        for (auto i = 0ull; i < elem_cnt; ++i) {
            if constexpr (std::is_same_v<matrix_elem_t, net_decimal>) *(ans + i) = (*(src + i)).number_format;
            else *(ans + i) = (*(src + i));
        }
        if (ln_cnt > col_cnt) {
            ans = transposition(ans, ln_cnt, col_cnt);
            std::swap(ln_cnt, col_cnt);
        }
        auto rank_val = 0ull;
        for (auto i = 0ull; i < col_cnt; ++i) {
            bool elim_flag = true;
            for(auto j = rank_val; j < ln_cnt; ++j) {
                auto elim_val = *(ans + elem_pos(j, i, col_cnt));
                if (elim_val) {
                    if (elim_flag) {
                        if (elim_val != 1) for (auto k = i; k < col_cnt; ++k) *(ans + elem_pos(j, k, col_cnt)) /= elim_val;
                        if (j != i) ans = ln_col_swap(ans, ln_cnt, col_cnt, i, j);
                        ++ rank_val;
                        elim_flag = false;
                    }
                    else for(auto k = i; k < col_cnt; ++k) *(ans + elem_pos(j, k, col_cnt)) -= *(ans + elem_pos(i, k, col_cnt)) * elim_val;
                }
            }
        }
        recycle(ans);
        return rank_val;
    }
    else return NAN;
}

callback_matrix matrix_elem_t det(const matrix_ptr src, uint64_t dim_cnt)
{
    if(src && dim_cnt) {
        if (dim_cnt == 1) return *src;
        auto ac = init<matrix_elem_t>((dim_cnt - 1) * (dim_cnt - 1));
        auto mov = 0;
        matrix_elem_t sum = 0;
        for (auto arow = 0ull; arow < dim_cnt; ++arow) {
            for (auto brow = 0ull; brow < dim_cnt - 1; ++brow) {
                mov = arow > brow ? 0 : 1;
                for (auto i = 0ull; i < dim_cnt - 1; ++i)
                    *(ac + brow * (dim_cnt - 1) + i) = *(src + (brow + mov) * dim_cnt + i + 1);
            }
            matrix_elem_t flag = (arow % 2 == 0 ? 1 : -1);
            sum += flag * (*(src + arow * dim_cnt)) * det(ac, dim_cnt - 1);
        }
        recycle(ac);
        return sum;
    }
    else return NAN;
}

// ans_dim_cnt = dim_cnt - 1
callback_matrix matrix_ptr cofactor(const matrix_ptr src, uint64_t ln, uint64_t col, uint64_t dim_cnt) {
    if(src && dim_cnt && ln < dim_cnt && col < dim_cnt) {
        auto elem_cnt     = dim_cnt * dim_cnt,
             ans_dim_cnt  = dim_cnt - 1,
             ans_elem_cnt = ans_dim_cnt * ans_dim_cnt;
        auto ans = init<matrix_elem_t>(ans_elem_cnt);
        for(auto i = 0ull; i < elem_cnt; ++i) {
            auto curr_pos = elem_pos(i, dim_cnt);
            if (curr_pos.ln != ln && curr_pos.col != col) {
                auto curr_ans_ln  = curr_pos.ln,
                     curr_ans_col = curr_pos.col;
                if (curr_ans_ln > ln) -- curr_ans_ln;
                if (curr_ans_col > col) -- curr_ans_col;
                *(ans + elem_pos(curr_ans_ln, curr_ans_col, ans_dim_cnt)) = *(src + i);
            }
        }
        return ans;
    }
    return nullptr;
}

callback_matrix matrix_elem_t cofactor_algebra(const matrix_ptr src, uint64_t ln, uint64_t col, uint64_t dim_cnt) {
    auto cof = cofactor(src, ln, col, dim_cnt);
    if (cof == nullptr) return 0;
    auto ans = det(cof, dim_cnt - 1);
    recycle(cof);
    return ans;
}

callback_matrix matrix_ptr adjugate(const matrix_ptr src, uint64_t dim_cnt) {
    matrix_ptr ans = nullptr;
    if (src && dim_cnt) {
        auto elem_cnt = dim_cnt * dim_cnt;
        ans = init<matrix_elem_t>(elem_cnt);
        for (auto i = 0ull; i < elem_cnt; ++i) {
            auto curr_pos = elem_pos(i, dim_cnt);
            auto coe = 1;
            if ((curr_pos.ln + curr_pos.col) % 2) coe = -1;
            *(ans + elem_pos(curr_pos.col, curr_pos.ln, dim_cnt)) = coe * cofactor_algebra(src, curr_pos.ln, curr_pos.col, dim_cnt);
        }
    }
    return ans;
}

callback_matrix matrix_ptr inverser(const matrix_ptr src, uint64_t dim_cnt) {
    auto det_val = det(src, dim_cnt);
    if (src && dim_cnt && det_val != 0) {
        auto elem_cnt = dim_cnt * dim_cnt;
        auto ajg_val = adjugate(src, dim_cnt);
        for (auto i = 0ull; i < elem_cnt; ++i) *(ajg_val + i) /= det_val;
        return ajg_val;
    } else return nullptr;
}

callback_matrix auto max_eigen(const matrix_ptr src, matrix_ptr &w, uint64_t dim_cnt, const matrix_elem_t &init_elem = 1, long double acc = 1e-8) {
    if (src && dim_cnt) {
        matrix_elem_t lambda      = 0,
                      lambda_temp = 0,
                      acc_dist    = 0;
        w      = init<matrix_elem_t>(dim_cnt);
        auto x = init<matrix_elem_t>(dim_cnt);
        fill(x, dim_cnt, init_elem);
        do {
            net_set<pos> temp_max_pos;
            lambda        = lambda_temp;
            auto temp     = copy(x, dim_cnt);
            auto temp_abs = absolute(temp, dim_cnt);
            auto max      = extremum(temp_max_pos, true, temp_abs, dim_cnt, 1, 0, dim_cnt - 1, 0, 0);
            for (auto i = 0ull; i < dim_cnt; ++i) *(x + i) /= max;
            for (auto i = 0ull; i < dim_cnt; ++i) *(w + i)  = *(x + i);
            matrix_elem_t sum = 0;
            for (auto i = 0ull; i < dim_cnt; ++i) sum      += *(w + i);
            for (auto i = 0ull; i < dim_cnt; ++i) *(w + i) /= sum;
            move(x, mult(src, dim_cnt, dim_cnt, x, 1));
            temp = copy(x, dim_cnt);
            lambda_temp = extremum(temp_max_pos, true, temp, dim_cnt, 1, 0, dim_cnt - 1, 0, 0);
            recycle(temp, temp_abs);
            if constexpr (std::is_same_v<matrix_elem_t, net_decimal>) acc_dist = (lambda - lambda_temp).absolute;
            else acc_dist = std::abs(lambda - lambda_temp);
        } while (acc_dist > acc);
        recycle(x);
        return lambda;
    } else return nullptr;
}

callback_matrix matrix_ptr LU(const matrix_ptr src, uint64_t dim_cnt) {
    if (src && dim_cnt) {
        auto elem_cnt = dim_cnt * dim_cnt;
        auto ans      = copy(src, elem_cnt);
        auto l        = init_identity<matrix_elem_t>(dim_cnt);
        matrix_ptr e  = nullptr;
        for (auto i = 0ull; i < dim_cnt - 1; ++i) {
            move(e, init_identity<matrix_elem_t>(dim_cnt));
            for (auto j = dim_cnt - 1; j > i; --j) *(e + elem_pos(j, i, dim_cnt)) = (-1) * (*(ans + elem_pos(j, i, dim_cnt))) / (*(ans + elem_pos(i, i, dim_cnt)));
            for (auto j = dim_cnt - 1; j > i; --j) *(l + elem_pos(j, i, dim_cnt)) = (-1) * (*(e + elem_pos(j, i, dim_cnt)));
            for (auto j = dim_cnt - 1; j > i; --j) *(ans + elem_pos(j, i, dim_cnt)) = (-1) * (*(e + elem_pos(j, i, dim_cnt)));
            move(ans, mult(e, dim_cnt, dim_cnt, ans, dim_cnt));
        }
        for (auto i = 0ull; i < dim_cnt; ++i) for (auto j = 0ull; j < i; ++j) *(ans + elem_pos(i, j, dim_cnt)) = *(l + elem_pos(i, j, dim_cnt));
        recycle(e, l);
        return ans;
    } else return nullptr;
}

callback_matrix matrix_ptr linear_equation(const matrix_ptr coefficient, const matrix_ptr b, uint64_t dim_cnt) {
    if (b && coefficient && dim_cnt) {
        auto elem_cnt = dim_cnt * dim_cnt;
        auto lu       = LU(coefficient, dim_cnt),
             e        = init_identity<matrix_elem_t>(dim_cnt),
             l        = copy(e, elem_cnt),
             u        = copy(e, elem_cnt),
             y        = init<matrix_elem_t>(dim_cnt),
             ans      = init<matrix_elem_t>(dim_cnt);
        for (auto i = 0ull; i < dim_cnt; ++i) for (auto j = 0ull; j < dim_cnt; ++j)
            if (i <= j) *(u + elem_pos(i, j, dim_cnt)) = *(lu + elem_pos(i, j, dim_cnt));
            else *(l + elem_pos(i, j, dim_cnt)) = *(lu + elem_pos(i, j, dim_cnt));
        for (auto i = 0ull; i < dim_cnt; ++i) if (i) {
            matrix_elem_t temp = 0;
            for (auto j = 0ull; j < i; ++j) temp += *(y + j) * (*(l + elem_pos(i, j, dim_cnt)));
            *(y + i) = *(b + i) - temp;
        } else *(y + i) = *(b + i);
        for (auto i = dim_cnt; i > 0; --i) if (i - dim_cnt) {
            matrix_elem_t temp = 0;
            for (auto n= i - 1; n < dim_cnt - 1; ++n) temp += *(u + elem_pos(i - 1, n + 1, dim_cnt)) * (*(ans  + n + 1));
            *(ans + i - 1) = (*(y + i - 1) - temp) / (*(u + elem_pos(i - 1, i - 1, dim_cnt)));
        } else *(ans + dim_cnt - 1) = *(y + dim_cnt - 1) / (*(u + elem_pos(dim_cnt - 1, dim_cnt - 1, dim_cnt)));
        recycle(lu, e, l, u, y);
        return ans;
    } else return nullptr;
}

callback_matrix matrix_ptr jacobi_iter(const matrix_ptr coefficient, const matrix_ptr b, uint64_t dim_cnt, const matrix_elem_t &init_elem = 1, long double acc = 1e-8) {
    if (coefficient && b && dim_cnt && init_elem > 0) {
        auto temp = init<matrix_elem_t>(dim_cnt),
             ans  = init<matrix_elem_t>(dim_cnt);
        fill(temp, dim_cnt, init_elem);
        bool flag = false;
        do {
            ans = copy(temp, dim_cnt);
            for (auto i = 0ull; i < dim_cnt; ++i) {
                matrix_elem_t sum = 0;
                for (auto j = 0ull; j < dim_cnt; ++j) if (i != j) sum += *(coefficient + elem_pos(i, j, dim_cnt)) * (*(ans + j));
                *(temp + i) = (*(b + i) - sum) / (*(coefficient + elem_pos(i, i, dim_cnt)));
            }
            for (auto i = 0ull; i < dim_cnt; ++i) {
                auto acc_dist = *(ans + i) - (*(temp + i));
                if constexpr (std::is_same_v<matrix_elem_t, net_decimal>) acc_dist = acc_dist.absolute;
                else acc_dist = std::abs(acc_dist);
                if (acc_dist >= acc) {
                    flag = true;
                    break;
                }
            }
            flag = false;
        } while (flag);
        ptr_reset(temp);
    } else return nullptr;
}

callback_matrix matrix_ptr rotate_rect(const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt, bool clockwise = true) {
    if (src && ln_cnt && col_cnt) {
        auto ans = init<matrix_elem_t>(ln_cnt * col_cnt);
        for (auto i = 0ull; i < ln_cnt; ++i) for (auto j = 0ull; j < col_cnt; ++j)
            if(clockwise) *(ans + elem_pos(j, ln_cnt - i - 1, ln_cnt)) = *(src + elem_pos(i, j, col_cnt));
            else *(ans + elem_pos(col_cnt - j - 1, i, ln_cnt)) = *(src + elem_pos(i, j, col_cnt));
        return ans;
    } else return nullptr;
}

callback_matrix matrix_ptr mirror_flip(const matrix_ptr src, uint64_t ln_cnt, uint64_t col_cnt, bool symmetry_vertical = true) {
    if (src && ln_cnt && col_cnt) {
        auto ans = init<matrix_elem_t>(ln_cnt * col_cnt);
        for(auto i = 0ull; i < ln_cnt; ++i) for (auto j = 0ull; j < col_cnt; ++j) {
            auto          src_pos   = col_cnt;
            matrix_elem_t curr_elem = 0;
            if(symmetry_vertical) {
                src_pos   = col_cnt - 1 - j;
                curr_elem = *(src + elem_pos(i, src_pos, col_cnt));
            } else {
                src_pos   = ln_cnt - 1 - i;
                curr_elem = *(src + elem_pos(src_pos, j, col_cnt));
            }
            *(ans + elem_pos(i, j, col_cnt)) = std::move(curr_elem);
        }
        return ans;
    } else return nullptr;
}

MATRIX_END