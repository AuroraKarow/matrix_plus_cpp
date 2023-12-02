NEUNET_BEGIN

struct net_decimal_data final {
    /* 1|1010|0100|0131|0800|1920|0001.0001|0400|0021|7080|0400|005
     * {1,1920,800,131,100,1010,1}     {1,400,21,7080,400,50}
     */
    net_set<uint64_t> it, ft;

    void value_copy(const net_decimal_data &src) {
        if (this == &src) return;
        it = src.it;
        ft = src.ft;
    }

    void value_move(net_decimal_data &&src) {
        if (this == &src) return;
        it = std::move(src.it);
        ft = std::move(src.ft);
    }

    net_decimal_data() {}
    net_decimal_data(const net_decimal_data &src) { value_copy(src); }
    net_decimal_data(net_decimal_data &&src) { value_move(std::move(src)); }

    void reset() {
        it.reset();
        ft.reset();
    }

    ~net_decimal_data() { reset(); }

    net_decimal_data &operator=(const net_decimal_data &src) {
        value_copy(src);
        return *this;
    }
    net_decimal_data &operator=(net_decimal_data &&src) {
        value_move(std::move(src));
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const net_decimal_data &src) {
        if (src.it.length) for (auto i = src.it.length; i; --i) {
            auto pow_seg = NEUNET_DEC_SEG_MAX / 10,
                 seg_tmp = src.it[i - 1];
            if (i < src.it.length) while (pow_seg > seg_tmp) {
                pow_seg /= 10;
                os << 0;
            }
            if (pow_seg) os << seg_tmp;
        } else os << 0;
        if (!src.ft.length) return os;
        os << '.';
        for (auto i = 0ull; i < src.ft.length; ++i) {
            auto pow_seg = NEUNET_DEC_SEG_MAX / 10,
                 seg_tmp = src.ft[i];
            while (pow_seg > seg_tmp) {
                pow_seg /= 10;
                os << 0;
            }
            if ((i + 1) == src.ft.length) while (!(seg_tmp % 10)) seg_tmp /= 10;
            if (seg_tmp) os << seg_tmp;
        }
        return os;
    }
};

uint64_t dec_dig_cnt(uint64_t seg_last_idx) { return seg_last_idx * NEUNET_DEC_DIG_MAX; }
uint64_t dec_dig_cnt(const net_decimal_data &src, bool is_it) {
    if ((is_it && !src.it.length) || !(is_it || src.ft.length)) return 0;
    uint64_t ans = 0,   
             tmp = 0;
    if(is_it){
        tmp = src.it.length - 1;
        ans = dec_dig_cnt(tmp) + std::log10(src.it[tmp]) + 1;
        return ans;
    }
    ans = dec_dig_cnt(src.ft.length);
    tmp = src.ft[src.ft.length - 1];
    while (!(tmp % 10)) {
        tmp /= 10;
        --ans;
    }
    return ans;
}

bool dec_is_int(const net_decimal_data &src) { return !src.ft.length; }

bool dec_is_zero(const net_decimal_data &src) { return !(src.it.length || src.ft.length); }

bool dec_is_one(const net_decimal_data &src) { return !src.ft.length && src.it.length == 1 && src.it[0] == 1; }

net_decimal_data dec_init(bool &sgn, long double src) {
    sgn = src < 0;
    if (!src) return {};
    src = std::abs(src);
    net_decimal_data ans;
    auto it_seg = (uint64_t)src;
    auto ft_seg = src - it_seg;
    if (it_seg > 0) {
        ans.it.init(1);
        ans.it[0] = it_seg;
    }
    if (!ft_seg) return ans;
    ans.ft.init(1);
    uint64_t dig_cnt = 0;
    while (dig_cnt++ < NEUNET_DEC_DIG_MAX) {
        if (ft_seg) ft_seg *= 10;
        ans.ft[0] *= 10;
        ans.ft[0] += ft_seg;
        if (ft_seg) ft_seg -= (int)ft_seg;
    }
    return ans;
}
net_decimal_data dec_init(bool &sgn, const std::string &src) {
    // verify
    auto str_len  = (uint64_t)src.length(),
         dot_idx  = str_len,
         bgn_idx  = 0ull,
         end_idx  = dot_idx - 1;
    // symbol
    auto dot_flag = false;
    net_decimal_data ans;
    for (auto i = 0ull; i < str_len; ++i) if (src[i] < '0' || src[i] > '9') {
        if (src[i] == '.' && !dot_flag) {
            dot_flag = true;
            dot_idx  = i;
            continue;
        } 
        if (src[i] == '-' && !i) {
            sgn = true;
            ++bgn_idx;
            continue;
        }
        if (src[i] == '+' && !i) {
            ++bgn_idx;
            continue;
        }
        return ans;
    }
    // zero
    while (end_idx > dot_idx && src[end_idx] == '0') --end_idx;
    while (bgn_idx < dot_idx && src[bgn_idx] == '0') ++bgn_idx;
    if (end_idx == bgn_idx && dot_idx == bgn_idx) return ans;
    uint64_t seg_cnt = 0,
             seg_tmp = 0,
             tmp_cnt = 0;
    // integer
    if (bgn_idx < dot_idx) {
        tmp_cnt = dot_idx - bgn_idx;
        seg_cnt = tmp_cnt / NEUNET_DEC_DIG_MAX;
        if (tmp_cnt % NEUNET_DEC_DIG_MAX) ++seg_cnt;
        ans.it.init(seg_cnt, false);
        seg_cnt = 0;
        tmp_cnt = 1;
        for (auto i = dot_idx; i > bgn_idx; --i) {
            seg_tmp += (src[i - 1] - '0') * tmp_cnt;
            tmp_cnt *= 10;
            if (tmp_cnt == NEUNET_DEC_SEG_MAX) {
                ans.it[seg_cnt++] = seg_tmp;
                tmp_cnt           = 1;
                seg_tmp           = 0;
            }
        }
        if (seg_tmp) {
            ans.it[seg_cnt] = seg_tmp;
            seg_tmp     = 0;
        }
        seg_cnt = 0;
    }
    // float
    if (end_idx > dot_idx) {
        tmp_cnt   = end_idx - dot_idx;
        seg_cnt   = tmp_cnt / NEUNET_DEC_DIG_MAX;
        if (tmp_cnt % NEUNET_DEC_DIG_MAX) ++seg_cnt;
        ans.ft.init(seg_cnt, false);
        seg_cnt = 0;
        tmp_cnt = 0;
        for (auto i = dot_idx + 1; i <= end_idx; ++i) {
            seg_tmp *= 10;
            seg_tmp += src[i] - '0';
            ++tmp_cnt;
            if (tmp_cnt == NEUNET_DEC_DIG_MAX) {
                ans.ft[seg_cnt++] = seg_tmp;
                tmp_cnt           = 0;
                seg_tmp           = 0;
            }
        }
        if (seg_tmp) {
            while (tmp_cnt++ < NEUNET_DEC_DIG_MAX) seg_tmp *= 10;
            ans.ft[seg_cnt] = seg_tmp;
        }
    }
    return ans;
}

int64_t dec2i(bool sgn, const net_decimal_data &src) {
    if (dec_is_zero(src)) return 0;
    int64_t ans = src.it[0];
    if (sgn) ans = 0 - ans;
    return ans;
}

long double dec2f(bool sgn, const net_decimal_data &src) {
    if (dec_is_zero(src)) return 0;
    long double ans = 0;
    if (src.ft.length) {
        ans = src.ft[0];
        while (ans > 1) ans /= 10;
    }
    if (src.it.length) return ans += src.it[0];
    if (sgn) ans = 0 - ans;
    return ans;
}

net_decimal_data dec_float_part(net_decimal_data &int_part, const net_decimal_data &src) {
    int_part.it = src.it;
    int_part.ft.reset();
    net_decimal_data ans;
    ans.ft = src.ft;
    return ans;
}

/* unsigned decimal digit part
 * NEUNET_DEC_CMP_EQL first is equal to second
 * NEUNET_DEC_CMP_LES first is less than second
 * NEUNET_DEC_CMP_GTR first is greater than second
 */
int dec_comp(const net_decimal_data &fst, const net_decimal_data &snd) {
    if (fst.it.length > snd.it.length) return NEUNET_DEC_CMP_GTR;
    if (fst.it.length < snd.it.length) return NEUNET_DEC_CMP_LES;
    for (auto i = fst.it.length; i; --i) {
        auto idx = i - 1;
        if (fst.it[idx] > snd.it[idx]) return NEUNET_DEC_CMP_GTR;
        if (fst.it[idx] < snd.it[idx]) return NEUNET_DEC_CMP_LES;
    }
    for (auto i = 0; i < fst.ft.length; ++i) {
        if (i == snd.ft.length) return NEUNET_DEC_CMP_GTR;
        if (fst.ft[i] > snd.ft[i]) return NEUNET_DEC_CMP_GTR;
        if (fst.ft[i] < snd.ft[i]) return NEUNET_DEC_CMP_LES;
    }
    if (fst.ft.length == snd.ft.length) return NEUNET_DEC_CMP_EQL;
    else return NEUNET_DEC_CMP_LES;
}
int dec_comp(bool sgn, int comp_res) {
    if (comp_res == NEUNET_DEC_CMP_GTR && sgn) return NEUNET_DEC_CMP_LES;
    if (comp_res == NEUNET_DEC_CMP_LES && sgn) return NEUNET_DEC_CMP_GTR;
    return comp_res;
}

uint64_t dec_sub(bool &carry, uint64_t minu, uint64_t subt) {
    auto cay = carry;
    carry    = subt > minu;
    if (carry) return NEUNET_DEC_SEG_MAX - subt + minu - cay;
    auto ans = minu - subt;
    if (cay) return dec_sub(carry, ans, 1);
    return ans;
}

// unsigned number digit segment
uint64_t dec_add(bool &carry, uint64_t fst, uint64_t snd) {
    auto dif = NEUNET_DEC_SEG_MAX - fst;
    auto cay = carry;
    carry    = dif <= snd;
    if (carry) return snd - dif + cay;
    auto ans = fst + snd;
    if (cay) return dec_add(carry, ans, 1);
    return ans;
}
uint64_t dec_add(bool &carry, uint64_t fst, uint64_t snd, bool subt) { return subt ? dec_sub(carry, fst, snd) : dec_add(carry, fst, snd); }
/* unsigned segment number
 * minuhend (first) segment should be greater than subtrahend (second) segment for true value of parameter subtract
 */
net_decimal_data dec_add(const net_decimal_data &fst, const net_decimal_data &snd, bool subt = false) {
    auto it_len  = std::max(fst.it.length, snd.it.length),
         ft_len  = std::max(fst.ft.length, snd.ft.length),
         ans_idx = it_len + ft_len + 1;
    net_set<uint64_t> ans_tmp(ans_idx);
    net_decimal_data ans;
    auto carry = false;
    for (auto i = ft_len; i; --i) {
        auto idx_tmp = i - 1;
        auto seg_tmp = dec_add(carry,
                               fst.ft.length > idx_tmp ? fst.ft[idx_tmp] : 0,
                               snd.ft.length > idx_tmp ? snd.ft[idx_tmp] : 0,
                               subt);
        if (seg_tmp || ans_idx < ans_tmp.length) ans_tmp[--ans_idx] = seg_tmp;
    }
    if (ans_idx < ans_tmp.length) ans.ft = ans_tmp.sub_set(ans_idx, ans_tmp.length - 1);
    ans_idx = 0;
    for (auto i = 0ull; i < it_len; ++i) ans_tmp[ans_idx++] = dec_add(carry,
                                                                      fst.it.length > i ? fst.it[i] : 0,
                                                                      snd.it.length > i ? snd.it[i] : 0,
                                                                      subt);
    if (carry) ans_tmp[ans_idx++] = carry;
    else while (ans_idx && !ans_tmp[ans_idx - 1]) --ans_idx;
    if (ans_idx) ans.it = ans_tmp.sub_set(0, ans_idx - 1);
    return ans;
}

void dec_mul_coe(uint64_t coe[]) {
    coe[0] %= NEUNET_DEC_MUL_POW;
    coe[1] %= NEUNET_DEC_MUL_SQR;
    coe[1] /= NEUNET_DEC_MUL_POW;
    coe[2] /= NEUNET_DEC_MUL_SQR;
}

uint64_t dec_mul_val(const net_decimal_data &src, uint64_t idx) {
    if (idx < src.ft.length) return src.ft[src.ft.length - idx - 1];
    else return src.it[idx - src.ft.length];
}

void dec_mul_carry(uint64_t &carry, bool &ca_add) { if (ca_add) {
    ca_add = false;
    ++carry;
} }

uint64_t dec_mul(uint64_t &carry, uint64_t fst, uint64_t snd) {
    carry = 0;
    if (!(snd && fst)) return 0;
    if (snd < NEUNET_DEC_SEG_MAX / fst) return fst * snd;
    /*
    9999999|999999|999999 ^ 2 =
    9999999999999999998|0000000000000000001
                               999998000001
                        1999996000002
                2099997|6000003
          1999997800000|2
    99999980000001
    */
    uint64_t fst_coe[3] = {fst, fst, fst},
             snd_coe[3] = {snd, snd, snd},
             ans_coe[5] = {0};
    dec_mul_coe(fst_coe);
    dec_mul_coe(snd_coe);
    for (auto i = 0; i < 3; ++i) for (auto j = 0; j < 3; ++j) ans_coe[i + j] += fst_coe[i] * snd_coe[j];
    auto ca_add = false;
    auto ans_1  = ans_coe[1] * NEUNET_DEC_MUL_POW,
         ans_2  = ans_coe[2] % NEUNET_DEC_MUL_END * NEUNET_DEC_MUL_SQR,
         ans_3  = ans_coe[3] % 10 * NEUNET_DEC_MUL_CUB,
         ans    = dec_add(ca_add, ans_coe[0], ans_1);
    dec_mul_carry(carry, ca_add);
    ans = dec_add(ca_add, ans, ans_2);
    dec_mul_carry(carry, ca_add);
    ans = dec_add(ca_add, ans, ans_3);
    dec_mul_carry(carry, ca_add);
    ans_coe[2] /= NEUNET_DEC_MUL_END;
    ans_coe[3] /= 10;
    carry      += ans_coe[2] + ans_coe[3] + ans_coe[4] * (NEUNET_DEC_MUL_POW / 10);
    return ans;
}
void dec_mul(net_set<uint64_t> &ans_coe, uint64_t fst, uint64_t snd, uint64_t fst_idx, uint64_t snd_idx) {
    uint64_t coe_idx = fst_idx + snd_idx,
             carry   = 0;
    auto     ca_add  = false;
    ans_coe[coe_idx] = dec_add(ca_add, ans_coe[coe_idx], dec_mul(carry, fst, snd));
    carry           += ca_add;
    ca_add           = false;
    ++coe_idx;
    ans_coe[coe_idx] = dec_add(ca_add, ans_coe[coe_idx], carry);
    while (ca_add) {
        ++coe_idx;
        ans_coe[coe_idx] = dec_add(ca_add, ans_coe[coe_idx], 0);
    }
}
net_decimal_data dec_mul(const net_decimal_data &fst, const net_decimal_data &snd) {
    auto fst_len = fst.it.length + fst.ft.length,
         snd_len = snd.it.length + snd.ft.length,
         ft_len  = fst.ft.length + snd.ft.length;
    net_set<uint64_t> ans_coe(fst_len * snd_len * 2);
    for (auto i = 0ull; i < fst_len; ++i) {
        auto fst_val = dec_mul_val(fst, i);
        for (auto j = 0ull; j < snd_len; ++j) {
            auto snd_val = dec_mul_val(snd, j);
            dec_mul(ans_coe, fst_val, snd_val, i, j);
        }
    }
    net_decimal_data ans;
    auto idx = 0ull;
    while (idx < ft_len && !ans_coe[idx]) ++idx;
    if (idx < ft_len) {
        ans.ft = ans_coe.sub_set(idx, ft_len - 1);
        ans.ft.reverse();
    }
    idx = ans_coe.length;
    while (ft_len < idx && !ans_coe[idx - 1]) --idx;
    if (ft_len < idx) ans.it = ans_coe.sub_set(ft_len, idx - 1);
    return ans;
}

uint64_t dec_div(uint64_t divd_0, uint64_t divd_1, uint64_t divr) {
    /*
    4316845796773021554|0000000000000000000
                      / 9632515584713416829
                      = 4481535232212608781
    precision lost -> 4316845796773021554 / 0.9632515584713416829 + 0
                      = 4481535232212609024
                      * 9632515584713416829
    4316845796773021787|9395928199718864896
    or return 4481535232212609024
    */
    auto     divr_d = divr / 1e19;
    uint64_t ans_0  = divd_0 / divr_d + divd_1 / divr * 1.,
             divd_2 = 0,
             divd_3 = dec_mul(divd_2, divr, ans_0);
    if (divd_0 == divd_2 && divd_1 == divd_3) return ans_0;
    /*
    compare = 4316845796773021554|0000000000000000000 > 4316845796773021787|9395928199718864896
    true  4316845796773021554|0000000000000000000 - 4316845796773021787|9395928199718864896
    false 4316845796773021787|9395928199718864896 - 4316845796773021554|0000000000000000000
    = 233|9395928199718864896
    */
    auto cmp = divd_0 > divd_2 || (divd_0 == divd_2 && divd_1 > divd_3),
         tmp = false;
    if (cmp) {
        divd_3 = dec_sub(tmp, divd_1, divd_3);
        divd_2 = divd_0 - divd_2;
    } else {
        divd_3  = dec_sub(tmp, divd_3, divd_1);
        divd_2 -= divd_0;
    }
    if (tmp) --divd_2;
    /*
    242 = 233 / 0.9632515584713416829 + 9395928199718864896 / 9632515584713416829
    compare ->
    true  4481535232212609024 + 242
    truncating for rounding down
    -> false 4481535232212609024 - (242 + 1)
    Rounding up
    return 4481535232212608781
    */
    uint64_t ans_1 = divd_2 / divr_d + divd_3 * 1. / divr;
    if (cmp) return ans_0 + ans_1;
    return ans_0 - ans_1;
}
uint64_t dec_div_seg_cnt(const net_decimal_data &divd, const net_decimal_data &divr) {
    auto divd_seg_cnt = divd.it.length + divd.ft.length,
         divr_seg_cnt = divr.it.length + divr.ft.length;
    return (std::max)(divd.it.length + divd.ft.length, divr.it.length + divr.ft.length);
}

net_set<uint64_t> dec_div_seg(const net_decimal_data &src, uint64_t len) {
    net_set<uint64_t> ans(len);
    auto idx = 0ull;
    for (auto i = src.it.length; i; --i) ans[idx++] = src.it[i - 1];
    for (auto i = 0ull; i < src.ft.length; ++i) ans[idx++] = src.ft[i];
    return ans;
}

void dec_div_seg(const net_decimal_data &divd, net_set<uint64_t> &divd_seg, uint64_t &divd_it_len, const net_decimal_data &divr, net_set<uint64_t> &divr_seg, uint64_t &divr_it_len) {
    if (divr.it.length || divr.ft[0]) {
        auto seg_cnt = dec_div_seg_cnt(divd, divr);
        divd_seg     = dec_div_seg(divd, seg_cnt);
        divr_seg     = dec_div_seg(divr, seg_cnt);
        divd_it_len  = divd.it.length;
        divr_it_len  = divr.it.length;
        return;
    }
    divr_it_len  = 0;
    auto seg_cnt = 0ull;
    while (!divr.ft[seg_cnt]) ++seg_cnt;
    divd_it_len = divd.it.length + seg_cnt;
    auto divd_seg_cnt = divd.it.length + divd.ft.length;
    if (divd.ft.length < seg_cnt) divd_seg_cnt += seg_cnt - divd.ft.length;
    auto divr_seg_cnt = divr.ft.length - seg_cnt,
         coe_seg_len  = (std::max)(divr_seg_cnt, divd_seg_cnt);
    divd_seg = dec_div_seg(divd, coe_seg_len);
    divr_seg.init(coe_seg_len);
    for (auto i = seg_cnt; i < divr.ft.length; ++i) divr_seg[i - seg_cnt] = divr.ft[i];
}

int dec_div_comp(const net_set<uint64_t> &fst, const net_set<uint64_t> &snd, uint64_t fst_high = 0, uint64_t snd_high = 0) {
    if (fst_high > snd_high) return NEUNET_DEC_CMP_GTR;
    if (fst_high < snd_high) return NEUNET_DEC_CMP_LES;
    for (auto i = 0ull; i < fst.length; ++i) {
        if (fst[i] > snd[i]) return NEUNET_DEC_CMP_GTR;
        if (fst[i] < snd[i]) return NEUNET_DEC_CMP_LES;
    }
    return NEUNET_DEC_CMP_EQL;
}

// divd_high = 0, divd_sgn = false, ans_sgn = false, ans_seg_idx = 0, end = false
uint64_t dec_div_coe(net_set<uint64_t> &divd_seg, const net_set<uint64_t> &divr_seg, uint64_t &divd_high, uint64_t &ans_seg_idx, bool &divd_sgn, bool &ans_sgn, bool &div_end) {
    ans_sgn = divd_sgn;
    if (!divd_high && divd_seg[0] < divr_seg[0]) {
        divd_high = divd_seg[0];
        for (auto i = 0ull; i < divd_seg.length - 1; ++i) divd_seg[i] = divd_seg[i + 1];
        divd_seg[divd_seg.length - 1] = 0;
        ++ans_seg_idx;
        return 0;
    }
    auto ans = dec_div(divd_high, divd_seg[0], divr_seg[0]);
    net_set<uint64_t> prod_divr_ans(divr_seg.length);
    auto carry = 0ull;
    for (auto i = divr_seg.length; i; --i) {
        auto idx           = i - 1,
             mul_carry     = carry;
        auto add_carry     = false;
        prod_divr_ans[idx] = dec_add(add_carry, dec_mul(carry, ans, divr_seg[idx]), mul_carry);
        if (add_carry) ++carry;
    }
    auto cmp = dec_div_comp(prod_divr_ans, divd_seg, carry, divd_high);
    if (cmp == NEUNET_DEC_CMP_EQL) {
        div_end = true;
        return ans;
    }
    auto prod_gtr  = cmp == NEUNET_DEC_CMP_GTR;
    auto sub_carry = false;
    if (prod_gtr) divd_sgn = !divd_sgn;
    for (auto i = divd_seg.length; i; --i) {
        auto idx = i - 1;
        if (prod_gtr) divd_seg[idx] = dec_sub(sub_carry, prod_divr_ans[idx], divd_seg[idx]);
        else divd_seg[idx] = dec_sub(sub_carry, divd_seg[idx], prod_divr_ans[idx]);
    }
    if (sub_carry) {
        if (prod_gtr) --carry;
        else --divd_high;
    }
    if (carry == divd_high) {
        divd_high = 0;
        return ans;
    }
    if (prod_gtr) divd_high = carry - divd_high;
    else divd_high -= carry;
    return ans;
}

net_set<uint64_t> dec_div(const net_decimal_data &divd, const net_decimal_data &divr, uint64_t prec, uint64_t &ans_itsz, bool &ans_sgn, bool &div_end) {
    auto divd_high = 0ull,
         ans_idx   = 0ull,
         divd_itsz = 0ull,
         divr_itsz = 0ull,
         ans_segsz = 1 + prec / NEUNET_DEC_DIG_MAX;
    auto divd_sgn  = false;
    auto add_carry = 0;
    net_set<uint64_t> divd_seg, divr_seg;
    dec_div_seg(divd, divd_seg, divd_itsz, divr, divr_seg, divr_itsz);
    if (prec % NEUNET_DEC_DIG_MAX) ++ans_segsz;
    ans_itsz = 1;
    ans_sgn  = false;
    div_end  = false;
    if (divd_itsz > divr_itsz) {
        ans_itsz  += divd_itsz - divr_itsz;
        ans_segsz += ans_itsz;
        --ans_segsz;
    }
    if (divd_itsz < divr_itsz) ans_idx = divr_itsz - divd_itsz;
    net_set<uint64_t> ans_seg(ans_segsz);
    neunet_dec_loop {
        auto ans = dec_div_coe(divd_seg, divr_seg, divd_high, ans_idx, divd_sgn, ans_sgn, div_end);
        if (ans_idx == ans_seg.length) break;
        if (!ans) continue;
        auto carry = false;
        if (ans_sgn) ans_seg[ans_idx] = dec_sub(carry, ans_seg[ans_idx], ans);
        else ans_seg[ans_idx] = dec_add(carry, ans_seg[ans_idx], ans);
        for (auto j = ans_idx; j && carry; --j) {
            auto idx = j - 1;
            if (ans_sgn) ans_seg[idx] = dec_sub(carry, ans_seg[idx], 0);
            else ans_seg[idx] = dec_add(carry, ans_seg[idx], 0);
        }
        if (carry) {
            if (ans_sgn) --add_carry;
            else ++add_carry;
        }
        if (div_end) break;
    }
    return ans_seg;
}
net_decimal_data dec_div_ans(const net_set<uint64_t> &ans_seg, uint64_t ans_itsz) {
    auto ans_idx = 0ull;
    net_decimal_data ans;
    for (ans_idx = 0; ans_idx < ans_itsz; ++ans_idx) if (ans_seg[ans_idx]) break;
    if (ans_idx < ans_itsz) {
        ans.it.init(ans_itsz - ans_idx);
        for (auto i = 0ull; i < ans.it.length; ++i) ans.it[i] = ans_seg[ans_itsz - 1 - i];
    }
    for (ans_idx = ans_seg.length; ans_idx > ans_itsz; --ans_idx) if (ans_seg[ans_idx - 1]) break;
    if (ans_itsz < ans_idx) {
        ans.ft.init(ans_idx - ans_itsz);
        for (auto i = 0ull; i < ans.ft.length; ++i) ans.ft[i] = ans_seg[ans_itsz + i];
    }
    return ans;
}
net_decimal_data dec_div(const net_decimal_data &divd, const net_decimal_data &divr, uint64_t prec) {
    auto ans_itsz = 0ull;
    auto ans_sgn  = false,
         div_end  = false;
    auto ans_seg  = dec_div(divd, divr, prec, ans_itsz, ans_sgn, div_end);
    return dec_div_ans(ans_seg, ans_itsz);
}

// return quotient
net_decimal_data dec_rem(net_decimal_data &divd_rem, const net_decimal_data &divr) {
    auto cmp = dec_comp(divd_rem, divr);
    if (cmp == NEUNET_DEC_CMP_EQL) {
        divd_rem.it.reset();
        auto sgn = false;
        return dec_init(sgn, 1);
    }
    if (cmp == NEUNET_DEC_CMP_LES) return {};
    auto ans_itsz = 0ull,
         ans_idx  = 0ull;
    auto ans_sgn  = false,
         div_end  = false;
    auto ans_seg  = dec_div(divd_rem, divr, 0, ans_itsz, ans_sgn, div_end);
    auto ans_quot = dec_div_ans(ans_seg, ans_itsz);
    if (div_end) {
        divd_rem.reset();
        return ans_quot;
    }
    if (ans_sgn && !dec_is_one(ans_quot)) ans_quot = dec_add(ans_quot, dec_init(div_end, 1), true);
    divd_rem = dec_add(divd_rem, dec_mul(ans_quot, divr), true);
    return ans_quot;
}

bool dec_gcd(net_decimal_data &fst, net_decimal_data &snd) {
    auto comp = dec_comp(fst, snd);
    if (comp == NEUNET_DEC_CMP_EQL) return false;
    if (comp == NEUNET_DEC_CMP_GTR) {
        dec_rem(fst, snd);
        // std::cout << fst << '\n' << snd << '\n' << std::endl;
        if (dec_is_zero(fst)) return true;
    }
    if (comp == NEUNET_DEC_CMP_LES) {
        dec_rem(snd, fst);
        // std::cout << fst << '\n' << snd << '\n' << std::endl;
        if (dec_is_zero(snd)) return false;
    }
    return dec_gcd(fst, snd);
}

void dec_e10_rsh(uint64_t &src, uint64_t rem_pow, uint64_t seg_pow, uint64_t &rem_dig) { if (seg_pow) {
    auto rem_tmp = src % rem_pow;
    src         /= rem_pow;
    src         += rem_dig * seg_pow;
    rem_dig      = rem_tmp;
} }

net_decimal_data dec_e10_rsh(const net_decimal_data &src, uint64_t dig_cnt) {
    auto seg_dif = dig_cnt / NEUNET_DEC_DIG_MAX,
         rem_dig = dig_cnt % NEUNET_DEC_DIG_MAX,
         rem_pow = 1ull,
         seg_pow = 0ull;
    if (rem_dig) {
        auto tmp = rem_dig;
        while (tmp--) rem_pow *= 10;
        seg_pow = NEUNET_DEC_SEG_MAX / rem_pow;
    }
    auto ans_ft_idx = 0;
    net_decimal_data ans;
    ans.ft.init(src.ft.length + seg_dif + (rem_dig ? 1 : 0));
    rem_dig = 0;
    if (src.it.length > seg_dif) {
        ans.it.init(src.it.length - seg_dif);
        auto src_it_idx = src.it.length;
        for (auto i = ans.it.length; i; --i) {
            auto idx     = i - 1;
            ans.it[idx]  = src.it[--src_it_idx];
            dec_e10_rsh(ans.it[idx], rem_pow, seg_pow, rem_dig);
        }
        seg_dif = ans.it.length - 1;
        if (!ans.it[seg_dif]) {
            if (seg_dif--) ans.it = ans.it.sub_set(0, seg_dif);
            else ans.it.reset();
        }
        for (auto i = src_it_idx; i; --i) {
            ans.ft[ans_ft_idx] = src.it[i - 1];
            dec_e10_rsh(ans.ft[ans_ft_idx++], rem_pow, seg_pow, rem_dig);
        }
    } else {
        ans_ft_idx = seg_dif - src.it.length;
        if (src.it.length) for (auto i = src.it.length; i; --i) {
            ans.ft[ans_ft_idx] = src.it[i - 1];
            dec_e10_rsh(ans.ft[ans_ft_idx++], rem_pow, seg_pow, rem_dig);
        }
    }
    for (auto i = 0ull; i < src.ft.length; ++i) {
        ans.ft[ans_ft_idx] = src.ft[i];
        dec_e10_rsh(ans.ft[ans_ft_idx++], rem_pow, seg_pow, rem_dig);
    }
    if (rem_dig) {
        dec_e10_rsh(ans.ft[ans_ft_idx++], rem_pow, seg_pow, rem_dig);
        return ans;
    }
    while (ans_ft_idx && !ans.ft[ans_ft_idx - 1]) --ans_ft_idx;
    if (ans_ft_idx--) ans.ft = ans.ft.sub_set(0, ans_ft_idx);
    else ans.ft.reset();
    return ans;
}

void dec_e10_lsh(uint64_t &src, uint64_t rem_pow, uint64_t seg_pow, uint64_t &rem_dig) { if (rem_pow) {
    auto rem_tmp = src / rem_pow;
    src         %= rem_pow;
    src         *= seg_pow;
    src         += rem_dig;
    rem_dig      = rem_tmp;
} }

net_decimal_data dec_e10_lsh(const net_decimal_data &src, uint64_t dig_cnt) {
    auto seg_dif = dig_cnt / NEUNET_DEC_DIG_MAX,
         rem_dig = dig_cnt % NEUNET_DEC_DIG_MAX,
         rem_pow = 0ull,
         seg_pow = 1ull;
    if (rem_dig) {
        auto tmp = rem_dig;
        while (tmp--) seg_pow *= 10;
        rem_pow = NEUNET_DEC_SEG_MAX / seg_pow;
    }
    auto ans_it_idx = 0ull;
    net_decimal_data ans;
    ans.it.init(src.it.length + seg_dif + (rem_dig ? 1 : 0));
    rem_dig = 0;
    if (src.ft.length > seg_dif) {
        ans.ft.init(src.ft.length - seg_dif);
        auto src_ft_idx = src.ft.length;
        for (auto i = ans.ft.length; i; --i) {
            auto idx    = i - 1;
            ans.ft[idx] = src.ft[--src_ft_idx];
            dec_e10_lsh(ans.ft[idx], rem_pow, seg_pow, rem_dig);
        }
        seg_dif = ans.ft.length - 1;
        if (!ans.ft[seg_dif]) {
            if (seg_dif--) ans.ft = ans.ft.sub_set(0, seg_dif);
            else ans.ft.reset();
        }
        while (src_ft_idx) {
            ans.it[ans_it_idx] = src.ft[--src_ft_idx];
            dec_e10_lsh(ans.it[ans_it_idx++], rem_pow, seg_pow, rem_dig);
        }
    } else {
        ans_it_idx = seg_dif - src.ft.length;
        if (src.ft.length) for (auto i = src.ft.length; i; --i) {
            ans.it[ans_it_idx] = src.ft[i - 1];
            dec_e10_lsh(ans.it[ans_it_idx++], rem_pow, seg_pow, rem_dig);
        }
    }
    for (auto i = 0ull; i < src.it.length; ++i) {
        ans.it[ans_it_idx] = src.it[i];
        dec_e10_lsh(ans.it[ans_it_idx++], rem_pow, seg_pow, rem_dig);
    }
    if (rem_dig) {
        dec_e10_lsh(ans.it[ans_it_idx++], rem_pow, seg_pow, rem_dig);
        return ans;
    }
    while (ans_it_idx && !ans.it[ans_it_idx - 1]) --ans_it_idx;
    if (ans_it_idx--) ans.it = ans.it.sub_set(0, ans_it_idx);
    else ans.it.reset();
    return ans;
}

// 5 -> 100000
net_decimal_data dec_e10(uint64_t e10) {
    net_decimal_data ans;
    if (!e10) {
        ans.it.init(1);
        ans.it[0] = 1;
        return ans;
    }
    auto seg_cnt = e10 / NEUNET_DEC_DIG_MAX + 1,
         dig_cnt = e10 % NEUNET_DEC_DIG_MAX;
    ans.it.init(seg_cnt);
    if (!dig_cnt) {
        ans.it[seg_cnt - 1] = 1;
        return ans;
    }
    seg_cnt = 1;
    for (auto i = 0; i < dig_cnt; ++i) seg_cnt *= 10;
    ans.it[ans.it.length - 1] = seg_cnt;
    return ans;
}
/* 100000000000000000000
 * total == true -> 1 * 100000000000000000000
 * total == false -> 10 * 10000000000000000000
 */
void dec_e10(net_decimal_data &e10, net_decimal_data &src, bool total = false) {
    e10.reset();
    e10.it = {1};
    if (!src.ft.length) goto it_part;
    src.ft.reverse();
    if (src.it.length || src.ft[src.ft.length - 1]) src.it = src.ft.unit(src.it);
    else for (auto i = src.ft.length; i; --i) if (src.ft[i - 1]) {
        src.it = src.ft.sub_set(0, i - 1);
        break;
    }
    e10.reset();
    e10.ft.init(src.ft.length);
    e10.ft[src.ft.length - 1] = 1;
    src.ft.reset();
    it_part: if (!src.it.length || src.it[0]) goto total_part;
    for (auto i = 0ull; i < src.it.length; ++i) if (src.it[i]) {
        src.it = src.it.sub_set(i, src.it.length - 1);
        if (i < e10.ft.length) {
            e10.ft.init(e10.ft.length - i, false);
            e10.ft[e10.ft.length - 1] = NEUNET_DEC_SEG_MAX / 10;
        } else {
            e10.it.init(i - e10.ft.length + 1, false);
            e10.ft.reset();
            e10.it[e10.it.length - 1] = 1;
        }
        break;
    }
    total_part: if (!total || src.it[0] % 10) return;
    auto seg_tmp = src.it[0],
         dig_cnt = 0ull;
    do {
        seg_tmp /= 10;
        ++dig_cnt;
    } while (!(seg_tmp % 10));
    src = dec_e10_rsh(src, dig_cnt);
    e10 = dec_e10_lsh(e10, dig_cnt);
}

void dec_truncate(net_decimal_data &src, uint64_t prec) {
    auto dig_cnt = dec_dig_cnt(src, false);
    auto sgn_tmp = false;
    if (prec >= dig_cnt) return;
    net_decimal_data carry;
    if (prec) {
        auto carry_seg = NEUNET_DEC_SEG_MAX / 10,
             seg_cnt   = 1ull;
        for (auto i = 0ull; i < prec; ++i) if (carry_seg == 1) {
            carry_seg = NEUNET_DEC_SEG_MAX / 10;
            ++seg_cnt;
        } else carry_seg /= 10;
        carry.ft.init(seg_cnt, false);
        carry.ft[seg_cnt - 1] = carry_seg;
    } else carry = dec_init(sgn_tmp, 1);
    prec += 2;
    uint64_t tgt_idx = prec / NEUNET_DEC_DIG_MAX,
             seg_len = prec % NEUNET_DEC_DIG_MAX,
             dig_idx = NEUNET_DEC_DIG_MAX;
    if (seg_len) ++tgt_idx;
    else seg_len = dig_idx;
    src.ft.init(tgt_idx--);
    auto curr_seg = src.ft[tgt_idx],
         seg_pow  = 1ull,
         last_dig = 0ull;
    // get last 2 digits after accuracy digits
    while (dig_idx-- > seg_len) {
        last_dig         = seg_pow * (curr_seg % 10);
        curr_seg        /= 10;
        seg_pow         *= 10;
        src.ft[tgt_idx] -= last_dig;
    }
    last_dig = curr_seg % 10;
    curr_seg /= 10;
    // last digit of the 2;
    if (last_dig >= 5) {
        // carry
        src = dec_add(src, carry);
        ++curr_seg;
    }
    last_dig *= seg_pow;
    seg_pow  *= 10;
    if (dig_idx) {
        // not the end of current segment
        src.ft[tgt_idx] -= last_dig;
        last_dig         = curr_seg % 10;
        if (last_dig < 5) src.ft[tgt_idx] -= last_dig * seg_pow;
        while (tgt_idx && !src.ft[tgt_idx]) --tgt_idx;
        if (tgt_idx) src.ft = src.ft.sub_set(0, tgt_idx);
        else src.ft.reset();
    } else {
        // current segment ends
        src.ft.init(tgt_idx--);
        last_dig = src.ft[tgt_idx] % 10;
        if (last_dig < 5) src.ft[tgt_idx] -= last_dig;
    }
}

NEUNET_END