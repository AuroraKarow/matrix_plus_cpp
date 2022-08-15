NEUNET_BEGIN

/* polynomial */

net_set<uint8_t> dec_polynomial_add(bool &ans_sgn, const net_set<uint8_t> &coe_a, bool sgn_a, const net_set<uint8_t> &coe_b, bool sgn_b) {
    net_set<uint8_t> ans;
    if (coe_a.length != coe_b.length) return ans;
    auto sgn_same = sgn_a == sgn_b;
    auto ans_len  = sgn_same ? coe_a.length + 1 : coe_a.length;
    auto conj_a   = sgn_a ? sgn_same ? 1 : -1 : 1,
         conj_b   = sgn_b ? sgn_same ? 1 : -1 : 1,
         carry    = 0;
         ans_sgn  = false;
    ans.init(ans_len);
    for (auto i = 0ull; i < coe_a.length; ++i) {
             carry   += conj_a * coe_a[i] + conj_b * coe_b[i];
        auto curr_bit = carry % 10;
             carry   /= 10;
        if(curr_bit < 0)
        {
            curr_bit = 10 + curr_bit;
            --carry;
        }
        ans[i] = curr_bit;
    }
    if (carry < 0) {
        for (auto i = 0ull; i < coe_a.length; ++i) {
            ans[i] = 10 - ans[i];
            if (i) --ans[i];
        }
        if (ans.length > coe_a.length) ans[coe_a.length] = (-1) * carry - 1;
        ans_sgn = true;
    } else {
        if(ans.length > coe_a.length) ans[coe_a.length] = carry;
        ans_sgn = sgn_a && sgn_same;
    }
    return ans;
}

net_set<uint64_t> dec_polynomial_mult(const net_set<uint8_t> &coe_a, const net_set<uint8_t> &coe_b)
{
    net_set<uint64_t> ans(coe_a.length + coe_b.length - 1);
    for (auto i = 0ull; i < coe_a.length; ++i) for (auto j = 0ull; j < coe_b.length; ++j)
    {
        auto curr_a = coe_a[i],
             curr_b = coe_b[j];
        ans[i + j] += curr_a * curr_b;
    }
    return ans;
}

/* fourier transform */

bool dec_fft(std::complex<long double> *&src, uint64_t len, bool inverse = false) {
    auto flag = false;
    if(len == 1) return true;
    if(len % 2) return false;
    auto src_sub_left  = ptr_init<std::complex<long double>>(len >> 1),
         src_sub_right = ptr_init<std::complex<long double>>(len >> 1);
    for(auto i = 0ull; i < len; i += 2) {
        *(src_sub_left + (i >> 1))  = *(src + i);
        *(src_sub_right + (i >> 1)) = *(src + i + 1);
    }
    if (dec_fft(src_sub_left, len >> 1, inverse) && dec_fft(src_sub_right, len >> 1, inverse)) {
        // dft
        // idft
        auto conj = 1;
        if(inverse) conj = -1;
        // omega
        std::complex<long double> omega(1, 0),
                                  omega_unit(std::cos(2 * NEUNET_PI / len), conj * std::sin(2 * NEUNET_PI / len));
        len >>= 1;
        for(auto i = 0ull; i < len; ++i) {
            *(src + i)       = *(src_sub_left + i) + omega * (*(src_sub_right + i));
            *(src + i + len) = *(i + src_sub_left) - omega * (*(src_sub_right + i));
            omega           *= omega_unit;
        }
        flag = true;
    }
    ptr_reset(src_sub_left, src_sub_right);
    return flag;
}

bool dec_dft(std::complex<long double> *&c, std::complex<long double> *&a, std::complex<long double> *&b, uint64_t len) {
    if (c && a && b && len && dec_fft(a, len) && dec_fft(b, len)) {
        for (auto i = 0ull; i < len; ++i)
            if(c + i) *(c + i) = *(a + i) * (*(b + i));
            else return false;
        return true;
    } else return false;
}

bool dec_idft(uint64_t *&d, std::complex<long double> *&c, uint64_t len) {
    if(c && d && len && dec_fft(c, len, true)) {
        for(auto i = 0ull; i < len; ++i)
            if(d + i) *(d + i) = (uint64_t)((c + i)->real() / len + 0.5);
            else return false;
        return true;
    } else return false;
}

/* number theoretical transform */

uint64_t dec_euler_pow(uint64_t base, uint64_t times) {
    uint64_t ans = 1;
    while(times)
    {
        if(times & 1) ans = ans * base % NEUNET_EULER_MOD;
        base = base * base % NEUNET_EULER_MOD;
        times >>= 1;
    }
    return ans;
}

bool dec_ntt(uint64_t *&src, uint64_t len, bool inverse = false) {
    if (num_pad_pow(len, 2) || len < 4 && src == nullptr) return false;
    auto idx_bit_cnt = num_bit_cnt(len - 1);
    auto temp = ptr_init<uint64_t>(len);
    for (auto i = 0ull; i < len; ++i) *(temp + i) = *(src + num_bit_reverse(i, idx_bit_cnt));
    ptr_move(src, std::move(temp));
    for (auto i = 1ull; i < len; i <<= 1) {
        auto g_root = dec_euler_pow(NEUNET_EULER_MOD_G, (NEUNET_EULER_MOD - 1) / (i << 1));
        if (inverse) g_root = dec_euler_pow(g_root, NEUNET_EULER_MOD - 2);
        for (auto j = 0ull; j < len; j += (i << 1)) {
            auto g_unit = 1ull;
            for (auto k = 0ull; k < i; ++k) {
                auto g_unit_front       = *(src + j + k),
                     g_unit_rear        = g_unit * (*(src + i + j + k)) % NEUNET_EULER_MOD;
                *(src + j + k)          = (g_unit_front + g_unit_rear) % NEUNET_EULER_MOD;
                *(src + i + j + k)      = (g_unit_front - g_unit_rear + NEUNET_EULER_MOD) % NEUNET_EULER_MOD;
                g_unit                  = g_unit * g_root % NEUNET_EULER_MOD;
            }
        }
    }
    if (inverse) {
        auto inverser = dec_euler_pow(len, NEUNET_EULER_MOD - 2);
        for (auto i = 0ull; i < len; ++i) *(src + i) = *(src + i) * inverser % NEUNET_EULER_MOD;
    }
    return true;
}

bool dec_fnt(uint64_t *&coe_c, uint64_t *&coe_a, uint64_t *&coe_b, uint64_t len)
{
    if (coe_c && coe_a && coe_b && len && dec_ntt(coe_a, len) && dec_ntt(coe_b, len)) {
        for(auto i = 0ull; i < len; ++i) *(coe_c + i) = *(coe_a + i) * (*(coe_b + i)) % NEUNET_EULER_MOD;
        return true;
    } else return false;
}

bool dec_ifnt(uint64_t *&coe_c, uint64_t len) { return (coe_c && len && dec_ntt(coe_c, len, true)); }

/* high precision decimal */

class net_decimal {
protected:
    void value_copy(const net_decimal &src) {
        val[0] = src.val[0];
        val[1] = src.val[1];
        sgn    = src.sgn;
    }

    void value_move(net_decimal &&src) {
        val[0] = std::move(src.val[0]);
        val[1] = std::move(src.val[1]);
        sgn    = src.sgn;
        src.reset();
    }
    
    net_sequence<uint8_t> dec_coe(bool dec_part, bool from_low) const {
        if (val[dec_part].length == 0) return net_sequence<uint8_t>();
        net_sequence<uint8_t> ans;
        for (auto i = 0ull; i < val[dec_part].length; ++i) {
            auto curr_seg = val[dec_part][i],
                 curr_dig = 0ull;
            while (curr_seg) {
                ans.emplace_back(curr_seg % 10);
                curr_seg /= 10;
                ++curr_dig;
            }
            if (i + 1 != val[dec_part].length) while (curr_dig != NEUNET_DEC_DIG_MAX) {
                ans.emplace_back(0);
                ++curr_dig;
            }
        }
        if ((from_low && dec_part == ft) || (!from_low && dec_part == it)) ans.reverse();
        return ans;
    }
    void dec_coe(net_set<uint8_t> &&coe_src, bool dec_part, bool from_low) {
        val[dec_part].reset();
        if (coe_src.length == 0) return;
        net_sequence<uint64_t> ans;
        if ((from_low && dec_part == ft) || (!from_low && dec_part == it)) coe_src.reverse();
        auto end_pt = coe_src.length;
        // calibrate
        while (end_pt && coe_src[end_pt - 1] == 0) --end_pt;
        if (end_pt == 0) return;
        // digit
        uint64_t curr_seg   = 0,
                 curr_times = 1;
        for (auto i = 0; i < end_pt; ++i) {
            if (curr_times == NEUNET_DEC_SEG_MAX) {
                ans.emplace_back(curr_seg);
                curr_times = 1;
                curr_seg   = 0;
            }
            curr_seg   += curr_times * coe_src[i];
            curr_times *= 10;
        }
        ans.emplace_back(curr_seg);
        ans.shrink();
        val[dec_part] = std::move(ans);
        coe_src.reset();
    }

    net_sequence<bool> dec_bin(bool from_low) const {
        net_sequence<bool> ans;
        if (val[ft].length) return ans;
        auto temp     = *this;
             temp.sgn = false;
        while (temp > 1) {
            ans.emplace_back((bool)(temp % 2).to_num());
            temp = (temp / 2).dec_num_part(it);
        }
        ans.emplace_back(temp.to_num());
        if (!from_low) ans.reverse();
        return ans;
    }
    void dec_bin(net_set<bool> &&bin_set, bool dec_sgn, bool from_low) {
        reset();
        if (bin_set.length == 0) return;
        if (!from_low) bin_set.reverse();
        net_decimal bin_cnt = 1;
        for (auto i = 0ull; i < bin_set.length; ++i) {
            *this += (int)bin_set[i] * bin_cnt;
            bin_cnt *= 2;
        }
        sgn = dec_sgn;
        bin_set.reset();
    }

    static void dec_bin_1st_com(net_set<bool> &bin_set) { for (auto i = 0ull; i < bin_set.length; ++i) bin_set[i] = !bin_set[i]; }

    static void dec_bin_2nd_com(net_set<bool> &bin_set, bool from_low) {
        dec_bin_1st_com(bin_set);
        for (auto i = 0ull; i < bin_set.length; ++i) {
            auto curr_idx = bin_set.length - i - 1;
            if (from_low) curr_idx = i;
            if (bin_set[curr_idx]) {
                bin_set[curr_idx] = false;
                if (i + 1 == bin_set.length) bin_set.insert(0, true);
            } else {
                bin_set[curr_idx] = true;
                break;
            }
        }
    }

    static void dec_bin_2nd_com_inverse(net_set<bool> &bin_set, bool from_low) {
        for (auto i = 0ull; i < bin_set.length; ++i) {
            auto curr_idx = bin_set.length - i - 1;
            if (from_low) curr_idx = i;
            if (bin_set[curr_idx]) {
                bin_set[curr_idx] = false;
                break;
            } else {
                bin_set[curr_idx] = true;
                if (i + 1 == bin_set.length) {
                    from_low ? bin_set.emplace_back(false) : bin_set.insert(0, false);
                    break;
                }
            }
        }
        dec_bin_1st_com(bin_set);
    }

    net_set<bool> dec_bin_or_and_xor(bool &ans_sgn, const net_decimal &src, uint8_t boolean_type) const {
        auto bin     = dec_bin(true),
             src_bin = src.dec_bin(true);
        if (sgn) dec_bin_2nd_com(bin, true);
        if (src.sgn) dec_bin_2nd_com(src_bin, true);
        net_sequence<bool> ans_bin;
        auto bit = bin.length > src_bin.length ? bin.length : src_bin.length;
        for (auto i = 0ull; i < bit; ++i) {
            auto bit_curr = false,
                 src_bit_curr = false;
            if (sgn) bit_curr = true;
            if (src.sgn) src_bit_curr = true;
            if (i < bin.length) bit_curr = bin[i];
            if (i < src_bin.length) src_bit_curr = src_bin[i];
            ans_bin.emplace_back(num_booleam(bit_curr, src_bit_curr, boolean_type));
        }
        ans_sgn = (sgn && src.sgn && boolean_type == 1) || ((sgn || src.sgn) && boolean_type == 0) || (src.sgn != sgn && boolean_type == 2);
        if (ans_sgn) dec_bin_2nd_com_inverse(ans_bin, true);
        return ans_bin;
    }
    
    static void dec_acc() {
        if (net_dec_con_acc == calculate_accuracy) return;
        net_dec_con_acc = calculate_accuracy;
        auto acc = calculate_accuracy;
        std::string temp = "0.";
        while (--acc) temp.push_back('0');
        temp.push_back('1');
        net_dec_con = temp;
    }

protected:
    /* <real>                      <- write ->                   <real>
     * 1|101|001|000|131|080|019|200|001.000|104|000|021|708|004|000|05
     * write -> (i*10 +) <real>   -> (i*10 +)        <real>
     * {1,200,19,80,131,0,1,101,1}, {0,401,0,120,807,400,0,50}
     * read -> (%10 /=10 low)     -> (%10 /=10 high)
     */
    net_set<uint64_t> val[2];

    bool sgn = false;
    
    static const bool it = false,
                      ft = true;

    static uint64_t     net_dec_pi_acc,
                        net_dec_con_acc;
    static net_decimal  net_dec_pi_cnt,
                        net_dec_pi_frc,
                        net_dec_pi,
                        net_dec_con;

public:
    net_decimal() {}
    callback_dec_arg net_decimal(arg src) {
        if (src == 0) return;
        if (src < 0) {
            sgn  = true;
            src *= (-1);
        }
        uint64_t it_part = src;
        src -= it_part;
        if (it_part) {
            net_sequence<uint8_t> temp;
            while (it_part) {
                temp.emplace_back(it_part % 10);
                it_part /= 10;
            }
            dec_coe(std::move(temp), it, true);
        }
        if (src) {
            net_sequence<uint8_t> temp;
            for (auto i = 0ull; i < calculate_accuracy; ++i) {
                src *= 10;
                uint8_t curr_dig = src;
                temp.emplace_back(curr_dig);
                src -= curr_dig;
            }
            dec_coe(std::move(temp), ft, false);
        }
    }
    net_decimal(const ch_str src) { *this = net_decimal(std::string(src)); }
    net_decimal(const std::string &src) {
        std::string temp = src;
        // calibrate
        while (temp[0] == '0') temp.erase(0, 1);
        while (temp[temp.length() - 1] == '0') temp.erase(temp.length() -  1, temp.length());
        if (temp.length() == 1 && temp[0] == '.') return;
        // symbol check
        auto dot_idx = temp.length();
        auto dot_flag = false;
        if (temp[temp.length() - 1] == '.') {
            temp.erase(temp.length() - 1, temp.length());
            dot_flag = true;
        }
        if (temp[0] == '.') {
            dot_idx  = 0;
            temp[0]  = '0';
            dot_flag = true;
        }
        if (temp[0] < '0' || temp[0] > '9') {
            if (temp[0] == '-') sgn = true;
            else if(temp[0] == '+') sgn = false;
            else return;
            temp.erase(0, 1);
        }
        for (auto i = 0ull; i < temp.length(); ++i) {
            if ((temp[i] < '0' || temp[i] > '9') && dot_flag) return;
            if(temp[i] == '.') {
                if(dot_flag) return;
                dot_flag = true;
                dot_idx = i;
            }
        }
        // integer
        if(dot_idx > 0) {
            net_sequence<uint8_t> it_temp;
            for (auto i = dot_idx; i > 0; --i) it_temp.emplace_back(temp[i-1] - '0');
            dec_coe(std::move(it_temp), it, true);
        }
        // float
        if(dot_idx < temp.length()) {
            net_sequence<uint8_t> ft_temp;
            for(auto i = dot_idx + 1; i < temp.length(); ++i) ft_temp.emplace_back(temp[i] - '0');
            dec_coe(std::move(ft_temp), ft, false);
        }
    }
    net_decimal(net_decimal &&src) { value_move(std::move(src)); }
    net_decimal(const net_decimal &src) { value_copy(src); }

    long double to_num() const {
        auto        coe_it = dec_coe(it, false),
                    coe_ft = dec_coe(ft, true);
        uint64_t    ans_it = 0;
        long double ans    = 0;
        for (auto i = 0ull; i < coe_it.length; ++i) {
            ans_it *= 10;
            ans_it += coe_it[i];
        }
        for (auto i = 0ull; i < coe_ft.length; ++i) {
            ans += coe_ft[i];
            ans /= 10;
        }
        ans += ans_it;
        if (sgn) ans *= (-1);
        return ans;
    }

    ch_str to_str() const {
        if (*this == 0) {
            auto zero_str = ptr_init<char>(2);
            *(zero_str)     = '0';
            *(zero_str + 1) = '\0';
            return zero_str;
        }
        auto   coe_it = dec_coe(it, false),
               coe_ft = dec_coe(ft, false);
        ch_str ans    = nullptr;
        auto   cnt    = 0ull;
        if (coe_ft.length && coe_it.length) {
            ans = ptr_init<char>(coe_it.length + coe_ft.length + (sgn ? 3 : 2));
            if (sgn) *(ans + cnt++) = '-';
            for (auto i = 0ull; i < coe_it.length; ++i) *(ans + cnt++) = coe_it[i] + '0';
            *(ans + cnt++) = '.';
            for (auto i = 0ull; i < coe_ft.length; ++i) *(ans + cnt++) = coe_ft[i] + '0';
        }
        else if (coe_ft.length) {
            ans = ptr_init<char>(coe_ft.length + (sgn ? 4 : 3));
            if (sgn) *(ans + cnt++) = '-';
            *(ans + cnt++) = '0';
            *(ans + cnt++) = '.';
            for (auto i = 0ull; i < coe_ft.length; ++i) *(ans + cnt++) = coe_ft[i] + '0';
        }
        else {
            ans = ptr_init<char>(coe_it.length + (sgn ? 2 : 1));
            if (sgn) *(ans + cnt++) = '-';
            for (auto i = 0ull; i < coe_it.length; ++i) *(ans + cnt++) = coe_it[i] + '0';
        }
        *(ans + cnt) = '\0';
        return ans;
    }

    net_decimal dec_num_part(bool dec_part) const {
        net_decimal ans;
        ans.val[dec_part] = val[dec_part];
        ans.sgn           = sgn;
        return ans;
    }

    uint64_t dec_dig_cnt(bool dec_part) const {
        if (val[dec_part].length == 0) return 0;
        auto idx = val[dec_part].length - 1,
             ans = idx * NEUNET_DEC_DIG_MAX,
             seg = val[dec_part][idx],
             cnt = 0ull;
        while (seg) {
            ++cnt;
            seg /= 10;
        }
        return cnt + ans;
    }

    void dec_truncate(uint64_t dig_cnt) {
        if (dig_cnt >= dec_dig_cnt(ft)) return;
        if (dig_cnt == 0) {
            *this += 0.5;
            val[ft].reset();
        }
        std::string round = "0.";
        for (auto i = 0ull; i < dig_cnt; ++i) round.push_back('0');
        round.push_back('5');
        if (sgn) *this -= round;
        else *this += round;
        auto coe = dec_coe(ft, false);
        coe.cut(dig_cnt);
        dec_coe(std::move(coe), ft, false);
        if (val[0].length == 0 && val[1].length == 0) sgn = false;
    }

    net_decimal dec_recip() const {
        if (*this == 0) return 0;
        if (*this == 1) return 1;
        dec_acc();
        net_decimal curr = 0,
                    temp = *this;
        temp.sgn = false;
        std::string next_temp = "0.";
        auto it_len_temp      = dec_dig_cnt(it);
        if(it_len_temp) while(-- it_len_temp) next_temp.push_back('0');
        next_temp.push_back('1');
        net_decimal next = next_temp;
        do {
            curr = next;
            next = curr * (2 - temp * curr);
            if (!high_precision_mode) next.dec_truncate(calculate_accuracy);
        } while ((next - curr).__abs__() > net_dec_con);
        if (sgn) next.sgn = true;
        if (high_precision_mode) next.dec_truncate(calculate_accuracy);
        return next;
    }

    void reset() {
        val[0].reset();
        val[1].reset();
        sgn = false;
    }

    ~net_decimal() { reset(); }

private:
    static net_decimal dec_sin_prototype(const net_decimal &src) {
        dec_acc();
        net_decimal ans      = 0,
                    cnt      = 3,
                    coe      = -1,
                    frac_div = 6,
                    curr_ans = src,
                    frac_sgl = src * src,
                    pow_num  = frac_sgl * src;
        do {
            ans       = curr_ans;
            curr_ans += coe * pow_num * frac_div.dec_recip();
            pow_num  *= frac_sgl;
            coe.sgn   = !coe.sgn;
            frac_div *= (cnt + 1) * (cnt + 2);
            cnt      += 2;
            if (!high_precision_mode) curr_ans.dec_truncate(calculate_accuracy);
            if (!high_precision_mode) pow_num.dec_truncate(calculate_accuracy);
        } while ((curr_ans - ans).__abs__() > net_dec_con);
        if (high_precision_mode) curr_ans.dec_truncate(calculate_accuracy);
        return curr_ans;
    }

    static net_decimal dec_cos_prototype(const net_decimal &src) {
        dec_acc();
        net_decimal ans      = 0,
                    coe      = -1,
                    cnt      = 2,
                    curr_ans = 1,
                    frac_div = 2,
                    pow_num  = src * src,
                    frac_sgl = pow_num;
        do {
            ans       = curr_ans;
            curr_ans += coe * pow_num * frac_div.dec_recip();
            coe.sgn   = !coe.sgn;
            pow_num  *= frac_sgl;
            frac_div *= (cnt + 1) * (cnt + 2);
            cnt      += 2;
            if (!high_precision_mode) curr_ans.dec_truncate(calculate_accuracy);
            if (!high_precision_mode) pow_num.dec_truncate(calculate_accuracy);
        } while ((ans - curr_ans).__abs__() > net_dec_con);
        if (high_precision_mode) curr_ans.dec_truncate(calculate_accuracy);
        return curr_ans;
    }

    static net_decimal dec_ln_prototype(const net_decimal &src) {
        dec_acc();
        net_decimal ans      = 0,
                    cnt      = 1,
                    _src     = (src - 1) * (src + 1).dec_recip(),
                    curr_ans = _src,
                    pow_num  = _src;
        do
        {
            ans      = curr_ans;
            cnt     += 2;
            pow_num *= _src * _src;
            curr_ans = ans + pow_num * cnt.dec_recip();
            if (!high_precision_mode) curr_ans.dec_truncate(calculate_accuracy);
            if (!high_precision_mode) pow_num.dec_truncate(calculate_accuracy);
        } while ((curr_ans - ans).__abs__() > net_dec_con);
        curr_ans *= 2;
        if (high_precision_mode) curr_ans.dec_truncate(calculate_accuracy);
        return curr_ans;
    }

    static void dec_pi() {
        if (calculate_accuracy <= net_dec_pi_acc) return;
        dec_acc();
        net_dec_pi      = net_dec_pi * net_decimal(2).dec_recip();
        net_dec_pi_acc  = calculate_accuracy;
        net_decimal ans = 0;
        do
        {
            ++net_dec_pi_cnt;
            ans = net_dec_pi;
            net_dec_pi_frc *= net_dec_pi_cnt * (net_dec_pi_cnt * 2 + 1).dec_recip();
            net_dec_pi     += net_dec_pi_frc;
        } while ((ans - net_dec_pi).__abs__() > net_dec_con);
        net_dec_pi *= 2;
    }

    static net_decimal euler_period(const net_decimal &times, const net_decimal &k) { return (2 * k + 1) * net_dec_pi * times; }

    static net_decimal euler_formula(const net_decimal &times) {
        dec_pi();
        net_decimal cnt    = 0,
                    ans    = 0,
                    limit  = (3 * times + 2) * net_dec_pi,
                    curr_p = 0;
        do {
            curr_p   = euler_period(times, cnt);
            auto ans = dec_cos(curr_p),
                 img = dec_sin(curr_p);
            img.dec_truncate(10);
            if (!(img.val[0].length || img.val[1].length)) return ans;
            ++ cnt;
        } while (curr_p.__abs__() < limit.__abs__());
        ans = 0;
        return ans;
    }

public:
    friend net_decimal operator+(const net_decimal &src) { return src; }
    net_decimal operator+(const net_decimal &src) const {
        if (src == 0) return net_decimal(*this);
        if (*this == 0) return net_decimal(src);
        auto coe_ft     = dec_coe(ft, false),
             coe_it     = dec_coe(it, true),
             src_coe_ft = src.dec_coe(ft, false),
             src_coe_it = src.dec_coe(it, true);
        for(auto i = coe_ft.length; i < src_coe_ft.length; ++i) coe_ft.emplace_back(0);
        for(auto i = src_coe_ft.length; i < coe_ft.length; ++i) src_coe_ft.emplace_back(0);
        coe_ft.reverse();
        src_coe_ft.reverse();
        for(auto i = coe_it.length; i < src_coe_it.length; ++i) coe_it.emplace_back(0);
        for(auto i = src_coe_it.length; i < coe_it.length; ++i) src_coe_it.emplace_back(0);
        net_decimal ans;
        auto ans_coe = dec_polynomial_add(ans.sgn, coe_ft.unit(coe_it), sgn, src_coe_ft.unit(src_coe_it), src.sgn);
        auto begin_pt = 0ull;
        while (src_coe_ft.length - begin_pt > calculate_accuracy) ++begin_pt;
        ans.dec_coe(ans_coe.sub_set(begin_pt, src_coe_ft.length - 1), ft, true);
        ans.dec_coe(ans_coe.sub_set(src_coe_ft.length, ans_coe.length - 1), it, true);
        return ans;
    }
    callback_dec_arg friend net_decimal operator+(arg val, const net_decimal &src) { return src + val; }
    void operator+=(const net_decimal &src) { *this = *this + src; }
    callback_dec_arg friend void operator+=(arg &val, const net_decimal &src) { val = (val + src).to_num(); }
    net_decimal &operator++() {
        *this += 1;
        return *this;
    }
    net_decimal operator++(int) {
        auto temp = *this;
        ++(*this);
        return temp;
    }

    friend net_decimal operator-(const net_decimal &src) {
        auto temp = src;
        temp.sgn  = !temp.sgn;
        return temp;
    }
    net_decimal operator-(const net_decimal &src) const {
        auto temp = src;
        temp.sgn  = !temp.sgn;
        return *this + temp;
    }
    callback_dec_arg friend net_decimal operator-(arg val, const net_decimal &src) { return net_decimal(val) - src; }
    void operator-=(const net_decimal &src) { *this = *this - src; }
    callback_dec_arg friend void operator-=(arg &val, const net_decimal &src) { val = (val - src).to_num(); }
    net_decimal &operator--() {
        *this -= 1;
        return *this;
    }
    net_decimal operator--(int) {
        auto temp = *this;
        --(*this);
        return temp;
    }

    net_decimal operator*(const net_decimal &src) const {
        if (src == 0 || *this == 0) return 0;
        if (*this == 1) return net_decimal(src);
        if (src == 1) return net_decimal(*this);
        net_decimal ans;
        if (src == -1) {
            ans = *this;
            ans.sgn = !ans.sgn;
            return ans;
        }
        if (*this == -1) {
            ans = src;
            ans.sgn = !ans.sgn;
            return ans;
        }
        auto coe_a = dec_coe(it, false).unit(dec_coe(ft, false)),
             coe_b = src.dec_coe(it, false).unit(src.dec_coe(ft, false));
        coe_a.reverse();
        coe_b.reverse();
        auto coe_ans = dec_polynomial_mult(coe_a, coe_b);
        net_sequence<uint8_t> ft_temp,
                              it_temp;
        auto cnt        = 0ull,
             carry      = 0ull,
             ft_bit_cnt = dec_dig_cnt(ft) + src.dec_dig_cnt(ft);
        while (ft_bit_cnt || carry || cnt < coe_ans.length) {
            if (cnt < coe_ans.length) carry += coe_ans[cnt++];
            if (ft_bit_cnt) {
                ft_temp.emplace_back(carry % 10);
                --ft_bit_cnt;
            } else it_temp.emplace_back(carry % 10);
            carry /= 10;
        }
        if (sgn != src.sgn) ans.sgn = true;
        ans.dec_coe(std::move(it_temp), it, true);
        ans.dec_coe(std::move(ft_temp), ft, true);
        if (!high_precision_mode) ans.dec_truncate(calculate_accuracy);
        return ans;
    }
    callback_dec_arg friend net_decimal operator*(arg val, const net_decimal &src) { return src * val; }
    void operator*=(const net_decimal &src) { *this = *this * src; }
    callback_dec_arg friend void operator*=(arg &val, const net_decimal &src) { val = (val * src).to_num(); }

    net_decimal operator/(const net_decimal &src) const { return *this * src.dec_recip(); }
    callback_dec_arg friend net_decimal operator/(arg val, const net_decimal &src) { return net_decimal(val) / src; }
    void operator/=(const net_decimal &src) { *this = *this / src; }
    callback_dec_arg friend void operator/=(arg &val, const net_decimal &src) { val = (val / src).to_num(); }

    net_decimal operator%(const net_decimal &src) const {
        auto multiple = (*this * src.dec_recip()).dec_num_part(it);
        if (modulo_mode && multiple.sgn) --multiple;
        return *this - src * multiple;
    }
    callback_dec_arg friend net_decimal operator%(arg val, const net_decimal &src) { return net_decimal(val) % src; }
    void operator%=(const net_decimal &src) { *this = *this % src; }
    callback_dec_arg friend void operator%=(arg &val, const net_decimal &src) { val = (val % src).to_num(); }

    net_decimal operator~() {
        net_decimal ans = 0;
        if (*this == 0) return ans;
        auto bin_set = dec_bin(true);
        if (sgn) {
            dec_bin_2nd_com(bin_set, true);
            dec_bin_1st_com(bin_set);
            ans.dec_bin(std::move(bin_set), false, true);
        } else {
            dec_bin_1st_com(bin_set);
            dec_bin_2nd_com_inverse(bin_set, true);
            ans.dec_bin(std::move(bin_set), true, true);
        }
        return ans;
    }

    net_decimal operator<<(const net_decimal &src) const {
        auto bit = src.dec_num_part(it);
        auto bin = dec_bin(false);
        for (net_decimal i = 0; i < bit; ++i) bin.emplace_back(0);
        net_decimal ans;
        ans.dec_bin(std::move(bin), sgn, false);
        return ans;
    }
    callback_dec_arg friend net_decimal operator<<(arg val, const net_decimal &src) { return net_decimal(val) << src; }
    void operator<<=(const net_decimal &src) { *this = *this << src; }
    callback_dec_arg friend void operator<<=(arg &val, const net_decimal &src) { val = (val << src).to_num(); }

    net_decimal operator>>(const net_decimal &src) const {
        auto bit = src.dec_num_part(it);
        auto bin = dec_bin(false);
        if (bit >= bin.length) return 0;
        for (net_decimal i = 0; i < bit; ++i) bin.erase(bin.length - 1, false);
        net_decimal ans;
        ans.dec_bin(std::move(bin), sgn, false);
        return ans;
    }
    callback_dec_arg friend net_decimal operator>>(arg val, const net_decimal &src) { return net_decimal(val) >> src; }
    void operator>>=(const net_decimal &src) { *this = *this >> src; }
    callback_dec_arg friend void operator>>=(arg &val, const net_decimal &src) { val = (val >> src).to_num(); }

    net_decimal operator&(const net_decimal &src) const {
        auto ans_sgn = false;
        auto ans_bin = dec_bin_or_and_xor(ans_sgn, src, 1);
        net_decimal ans;
        ans.dec_bin(std::move(ans_bin), ans_sgn, true);
        return ans;
    }
    callback_dec_arg friend net_decimal operator&(arg val, const net_decimal &src) { return net_decimal(val) & src; }
    void operator&=(const net_decimal &src) { *this = *this & src; }
    callback_dec_arg friend void operator&=(arg &val, const net_decimal &src) { val = (val & src).to_num(); }

    net_decimal operator|(const net_decimal &src) const {
        auto ans_sgn = false;
        auto ans_bin = dec_bin_or_and_xor(ans_sgn, src, 0);
        net_decimal ans;
        ans.dec_bin(std::move(ans_bin), ans_sgn, true);
        return ans;
    }
    callback_dec_arg friend net_decimal operator|(arg val, const net_decimal &src) { return net_decimal(val) | src; }
    void operator|=(const net_decimal &src) { *this = *this | src; }
    callback_dec_arg friend void operator|=(arg &val, const net_decimal &src) { val = (val | src).to_num(); }

    net_decimal operator^(const net_decimal &src) const {
        auto ans_sgn = false;
        auto ans_bin = dec_bin_or_and_xor(ans_sgn, src, 2);
        net_decimal ans;
        ans.dec_bin(std::move(ans_bin), ans_sgn, true);
        return ans;
    }
    callback_dec_arg friend net_decimal operator^(arg val, const net_decimal &src) { return net_decimal(val) ^ src; }
    void operator^=(const net_decimal &src) { *this = *this ^ src; }
    callback_dec_arg friend void operator^=(arg &val, const net_decimal &src) { val = (val ^ src).to_num(); }

    bool operator<(const net_decimal &src) const {
        if (!sgn && src.sgn) return false;
        if (sgn && !src.sgn) return true;
        // integer part
        if (val[it].length < src.val[it].length) {
            if (sgn) return false;
            else return true;
        }
        if (val[it].length > src.val[it].length) {
            if (sgn) return true;
            else return false;
        }
        auto coe     = dec_coe(it, false),
             src_coe = src.dec_coe(it, false);
        auto dig_cnt = coe.length;
        if (dig_cnt > src_coe.length) dig_cnt = src_coe.length;
        for (auto i = 0ull; i < dig_cnt; ++i) {
            if (coe[i] > src_coe[i]) {
                if (sgn) return true;
                else return false;
            }
            if (coe[i] < src_coe[i]) {
                if (sgn) return false;
                else return true;
            }
        }
        // float part
        coe     = dec_coe(ft, false),
        src_coe = src.dec_coe(ft, false);
        dig_cnt = coe.length;
        if (dig_cnt > src_coe.length) dig_cnt = src_coe.length;
        for (auto i = 0ull; i < dig_cnt; ++i) {
            if (coe[i] > src_coe[i]) {
                if (sgn) return true;
                else return false;
            }
            if (coe[i] < src_coe[i]) {
                if (sgn) return false;
                else return true;
            }
        }
        if ((dig_cnt < src_coe.length && !sgn) || (dig_cnt > src_coe.length && sgn)) return true;
        return false;
    }
    callback_dec_arg friend bool operator<(arg val, const net_decimal &src) { return net_decimal(val) < src; }

    bool operator<=(const net_decimal &src) const { return *this == src || *this < src; }
    callback_dec_arg friend bool operator<=(arg val, const net_decimal &src) { return net_decimal(val) <= src; }

    bool operator>(const net_decimal &src) const { return !((*this < src) || (*this == src)); }
    callback_dec_arg friend bool operator>(arg val, const net_decimal &src) { return net_decimal(val) > src; }

    bool operator>=(const net_decimal &src) const { return !(*this < src); }
    callback_dec_arg friend bool operator>=(arg val, const net_decimal &src) { return net_decimal(val) >= src; }

    bool operator==(const net_decimal &src) const { return src.val[0] == val[0] && src.val[1] == val[1] && sgn == src.sgn; }
    callback_dec_arg friend bool operator==(arg val, const net_decimal &src) { return net_decimal(val) == src; }

    bool operator!=(const net_decimal &src) const { return !(*this == src); }
    callback_dec_arg friend bool operator!=(arg val, const net_decimal &src) { return net_decimal(val) != src; }

    net_decimal &operator=(const net_decimal &src) {
        value_copy(src);
        return *this;
    }
    net_decimal &operator=(net_decimal &&src) {
        value_move(std::move(src));
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &out, const net_decimal &src) {
        auto str_temp = src.to_str();
        out << str_temp;
        ptr_reset(str_temp);
        return out;
    }

    explicit operator long double () const { return to_num(); }

    std::string __to_str__() const {
        auto temp = to_str();
        std::string ans = temp;
        ptr_reset(temp);
        return ans;
    }

    net_decimal __int_part__() const { return dec_num_part(it); }

    net_decimal __ftp_part__() const { return dec_num_part(ft); }

    uint64_t __int_dig_cnt__() const { return dec_dig_cnt(it); }

    uint64_t __ftp_dig_cnt__() const { return dec_dig_cnt(ft); }

    net_decimal __abs__() const {
        auto ans = *this;
        ans.sgn  = false;
        return ans;
    }

    static net_decimal dec_sin(const net_decimal &src) {
        dec_pi();
        auto temp     = src;
        auto period   = 2 * net_dec_pi;
             temp.sgn = false;
        while((temp - period) > 0) temp -= period;
        if(temp > net_dec_pi) {
            temp    -= net_dec_pi;
            auto ans = dec_sin_prototype(temp);
            if(!src.sgn) ans.sgn = !ans.sgn;
            return ans;
        } else {
            auto ans = dec_sin_prototype(temp);
            if(src.sgn) ans.sgn = !ans.sgn;
            return ans;
        }
    }

    static net_decimal dec_cos(const net_decimal &src) {
        dec_pi();
        auto period = 2 * net_dec_pi;
        auto temp   = src;
        temp.sgn    = false;
        while((temp - period) > 0) temp -= period;
        if(temp > net_dec_pi)
        {
            temp    -= net_dec_pi;
            auto ans = dec_cos_prototype(temp);
            ans.sgn  = !ans.sgn;
            return ans;
        }
        else return dec_cos_prototype(temp);
    }

    static net_decimal dec_ln(const net_decimal &src) {
        if (src <= 0 || src == 1) return 0;
        if (src > 10) {
            net_decimal times = 0,
                        temp  = src;
            while (temp > 10) {
                ++ times;
                temp /= 10;
            }
            // ln(c) = ln(a * 10 ^ b) = ln(a) + bln10
            return dec_ln_prototype(temp) + times * dec_ln_prototype(10);
        }
        if (src > 1) return dec_ln_prototype(src);
        net_decimal times = 0,
                    temp  = src;
        while (temp < 1)
        {
            ++ times;
            temp *= 10;
        }
        // ln(c) = ln(a * 10 ^ -b) = ln(a) - bln10
        return dec_ln_prototype(temp) - times * dec_ln_prototype(10);
    }

    static net_decimal dec_exp(const net_decimal &src) {
        dec_acc();
        net_decimal ans      = 1,
                    cnt      = 1,
                    pow_num  = src,
                    curr_ans = 1,
                    fact_num = 1;
        do {
            ++ cnt;
            ans       = curr_ans;
            curr_ans  = ans + pow_num * fact_num.dec_recip();
            fact_num *= cnt;
            pow_num  *= src;
            if (!high_precision_mode) curr_ans.dec_truncate(calculate_accuracy);
            if (!high_precision_mode) pow_num.dec_truncate(calculate_accuracy);
        }
        while ((curr_ans - ans).__abs__() > net_dec_con);
        if (high_precision_mode) curr_ans.dec_truncate(calculate_accuracy);
        return ans;
    }

    static net_decimal dec_pow(const net_decimal &base, const net_decimal &times) {
        if (base == 0) return 0;
        if (base == 1 || times == 0) return 1;
        if (times == 1) return net_decimal(base);
        if (times == -1) return base.dec_recip();
        if (times.val[ft].length) if (base > 0) return dec_exp(times * dec_ln(base));
        else {
            auto temp      = base,
                 euler_ans = euler_formula(times);
            temp.sgn       = false;
            return dec_exp(times * dec_ln(temp)) * euler_ans;
        } else {
            auto ans        = base,
                 times_temp = times;
            times_temp.sgn  = false;
            for (auto i = 1ull; i < times_temp; ++i) ans *= base;
            if (times.sgn) return ans.dec_recip();
            else return ans;
        }
    }

    // Float digit remain of calculation accuracy truncating
    static uint64_t calculate_accuracy;

    // This may slow down the calculating speed
    static bool high_precision_mode;

    // Switch on the mode for modulo operation from off status for remain mode
    static bool modulo_mode;

    __declspec(property(get = __ftp_dig_cnt__,
                        put = dec_truncate))    uint64_t    float_digit_count;
    __declspec(property(get = __int_dig_cnt__)) uint64_t    integer_digit_count;
    __declspec(property(get = to_num))          long double number_format;
    __declspec(property(get = __to_str__))      std::string string_format;
    __declspec(property(get = __int_part__))    net_decimal integer_part;
    __declspec(property(get = __ftp_part__))    net_decimal float_part;
    __declspec(property(get = __abs__))         net_decimal absolute;
    __declspec(property(get = dec_recip))       net_decimal reciprocal;
};

bool        net_decimal::modulo_mode         = false,
            net_decimal::high_precision_mode = false;
uint64_t    net_decimal::calculate_accuracy  = 32,
            net_decimal::net_dec_pi_acc      = 0,
            net_decimal::net_dec_con_acc     = 0;
net_decimal net_decimal::net_dec_pi_cnt      = 0,
            net_decimal::net_dec_pi_frc      = 1,
            net_decimal::net_dec_pi          = 2,
            net_decimal::net_dec_con         = 0;

NEUNET_END

_STD_BEGIN

neunet::net_decimal pow(const neunet::net_decimal &_Xx, const neunet::net_decimal &_Yx) { return neunet::net_decimal::dec_pow(_Xx, _Yx); }

neunet::net_decimal exp(const neunet::net_decimal &_Xx) { return neunet::net_decimal::dec_exp(_Xx); }

_STD_END