NEUNET_BEGIN

class net_decimal {
public:
    // Change the operation mode between modulus and remainder operation - '%'
    bool modulus_mode = false;

    // Default precision of infinite decimal
    inline static uint64_t default_infinite_precision = 64;
    
    // Binary bit count for bit operation. If this value is 0, bit count would be changed automatically.
    uint64_t binary_bit_set = 0;

    // Float point digit count for infinite decimal
    uint64_t infinite_precision = default_infinite_precision;

    __declspec(property(get = abs))      net_decimal absolute;
    __declspec(property(get = to_float)) long double float_point_format;
    __declspec(property(get = to_int))   int64_t     integer_format;

protected:
    void value_assign(const net_decimal &src) {
        sgn                = src.sgn;
        red                = src.red;
        modulus_mode       = src.modulus_mode;
        binary_bit_set     = src.binary_bit_set;
        infinite_precision = src.infinite_precision;
    }

    void value_copy(const net_decimal &src) {
        value_assign(src);
        num = src.num;
        den = src.den;
    }

    void value_move(net_decimal &&src) {
        value_assign(src);
        num = std::move(src.num);
        den = std::move(src.den);
    }

    template <uint64_t opt_idx> static net_decimal binary_operator(const net_decimal &fst, const net_decimal &snd) {
        net_decimal ans, fst_red, snd_red;
        if (!dec_is_zero(fst.den)) {
            fst_red = fst;
            fst_red.reduct();
            if (!dec_is_zero(fst_red.den)) return ans;
            return binary_operator<opt_idx>(fst_red, snd);
        }
        if (!dec_is_zero(snd.den)) {
            snd_red = fst;
            snd_red.reduct();
            if (!dec_is_zero(snd_red.den)) return ans;
            return binary_operator<opt_idx>(fst, snd_red);
        }
        if (fst.num.ft.length || snd.num.ft.length) return ans;
        auto bit_set = 0ull;
        if (fst.binary_bit_set && snd.binary_bit_set) bit_set = (std::max)(fst.binary_bit_set, snd.binary_bit_set);
        if constexpr (opt_idx == NEUNET_DEC_BIN_OR) ans.num = dec_bit_or(fst.num, snd.num, fst.sgn, snd.sgn, bit_set);
        else if constexpr (opt_idx == NEUNET_DEC_BIN_AND) ans.num = dec_bit_and(fst.num, snd.num, fst.sgn, snd.sgn, bit_set);
        else if constexpr (opt_idx == NEUNET_DEC_BIN_XOR) ans.num = dec_bit_xor(fst.num, snd.num, fst.sgn, snd.sgn, bit_set);
        return ans;
    }

public:
    static net_decimal pi(uint64_t places = default_infinite_precision) {
        net_decimal ans;
        ans.num = dec_series_pi::value(places);
        return ans;
    }

    net_decimal() { num = dec_init(sgn, 0); }
    net_decimal(const net_decimal &src) { value_copy(src); }
    net_decimal(net_decimal &&src) { value_move(std::move(src)); }
    callback_arg net_decimal(const arg &src) {
        static_assert(neunet_dec_init_expr);
        num = dec_init(sgn, src);
    }

    void reduct() {
        if (red) return;
        red = true;
        dec_frac_red(num, den);
    }

    int64_t to_int() { return dec_frac2i(sgn, num, den); }

    long double to_float() { return dec_frac2f(sgn, num, den); }

    net_decimal int_part(net_decimal &float_part) const {
        float_part.num = num;
        float_part.den = den;
        net_decimal ans;
        ans.num = dec_frac_int_part(float_part.num, float_part.den);
        return ans;
    }

    net_decimal abs() const {
        if (!sgn) return *this;
        auto ans = *this;
        ans.sgn  = false;
        return ans;
    }

    net_decimal ln() const {
        if (dec_is_zero(num) || dec_frac_is_one(num, den) || sgn) return {};
        net_decimal ans;
        ans.num = dec_series_ln(ans.sgn, num, infinite_precision);
        if (dec_is_zero(den)) return ans;
        auto sgn = false;
        auto tmp = dec_series_ln(sgn, den, infinite_precision);
        ans.num  = dec_sub(ans.sgn, ans.num, ans.sgn, ans.den, sgn);
        ans.den.reset();
        return ans;
    }

    net_decimal exp() const {
        if (dec_is_zero(num)) return {1};
        net_decimal ans;
        ans.num = dec_series_exp(ans.sgn, num, den, sgn, infinite_precision);
        return ans;
    }

    net_decimal sin() const {
        net_decimal ans;
        auto num_term = num,
             den_term = den,
             fra_form = dec_init(ans.sgn, 1);
        ans.num = dec_series_sin_cos(ans.sgn, num, den, sgn, infinite_precision, num_term, den_term, fra_form);
        return ans;
    }

    net_decimal cos() const {
        net_decimal ans;
        auto num_term = dec_init(ans.sgn, 1),
             den_term = dec_init(ans.sgn, 0),
             fra_form = den_term;
        ans.num = dec_series_sin_cos(ans.sgn, num, den, sgn, infinite_precision, num_term, den_term, fra_form);
        return ans;
    }

    net_decimal pow(const net_decimal &times) const {
        if (dec_is_zero(num)) return {};
        if (dec_is_one(num)){
            if (!sgn) return {1};
            if (dec_is_zero(times.den)) {
                auto tmp = 0ull;
                if (times.num.ft.length) tmp = times.num.ft[times.num.ft.length - 1];
                else tmp = times.num.it[0];
                if (tmp % 2) return {-1};
                return {};
            }
            auto num_tmp = num,
                 den_tmp = den;
            dec_frac_red(num_tmp, den_tmp);
            if (dec_is_zero(num_tmp)) {
                net_decimal times_tmp;
                times_tmp.num = std::move(num_tmp);
                times_tmp.den = std::move(den_tmp);
                return pow(times_tmp);
            }
            auto tmp = 0ull;
            if (den_tmp.ft.length) tmp = den_tmp.ft[den_tmp.ft.length - 1];
            else tmp = den_tmp.it[0];
            if (!(tmp % 2)) return {};
            // den = -1
            if (num_tmp.ft.length) tmp = num_tmp.ft[num_tmp.ft.length - 1];
            else tmp = num_tmp.it[0];
            if (tmp % 2) return {1};
            return {-1};
        }
        if (dec_frac_is_zero(times.num, times.den)) return {1};
        net_decimal ans;
        if (dec_frac_is_one(times.num, times.den)) {
            if (!times.sgn) return *this;
            ans.num = den;
            ans.den = num;
            if (dec_is_zero(den)) ans.num = dec_init(ans.sgn, 1);
            return ans;
        }
        ans.num = dec_series_pow(ans.sgn, num, den, sgn, times.num, times.den, times.sgn, infinite_precision);
        return ans;
    }

    void reset() {
        sgn = false;
        red = true;
        num.reset();
        den.reset();
        modulus_mode       = false;
        binary_bit_set     = 0;
        infinite_precision = infinite_precision;
    }

    ~net_decimal() { reset(); }

protected:
    // numerator / denominator
    net_decimal_data num, den;

    bool sgn = false,
         red = true;

public:
    callback_dec_arg explicit operator arg() const {
        if (dec_is_zero(den)) {
            if constexpr (std::is_unsigned_v<neunet_dec_type(arg)>) {
                if (sgn) std::abort();
                if (num.it.length) return num.it[0];
                return 0;
            }
            if constexpr (std::is_integral_v<neunet_dec_type(arg)>) return dec2i(sgn, num);
            return dec2f(sgn, num);
        }
        net_decimal tmp;
        tmp.num = dec_div(num, den, infinite_precision);
        return arg(tmp);
    }

    net_decimal &operator=(const net_decimal &src) {
        value_copy(src);
        return *this;
    }
    net_decimal &operator=(net_decimal &&src) {
        value_move(std::move(src));
        return *this;
    }

    friend bool operator==(const net_decimal &fst, const net_decimal &snd) { return dec_frac_comp(fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn) == NEUNET_DEC_CMP_EQL; }
    friend std::strong_ordering operator<=>(const net_decimal &fst, const net_decimal &snd) {
        auto cmp = dec_frac_comp(fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn);
        if (cmp == NEUNET_DEC_CMP_LES) return std::strong_ordering::less;
        if (cmp == NEUNET_DEC_CMP_GTR) return std::strong_ordering::greater;
        return std::strong_ordering::equal;
    }

    net_decimal operator+() { return *this; }
    friend net_decimal operator+(const net_decimal &fst, const net_decimal &snd) {
        net_decimal ans;
        ans.sgn = dec_frac_add(ans.num, ans.den, fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn);
        ans.red = dec_frac1(ans.num, ans.den) || dec_frac0(ans.num, ans.den) || dec_is_zero(ans.den);
        return ans;
    }
    net_decimal &operator+=(const net_decimal &snd) {
        *this = *this + snd;
        return *this;
    }
    callback_dec_arg friend void operator+=(arg &fst, const net_decimal &snd) { fst += snd.to_float(); }
    net_decimal &operator++() {
        *this += 1;
        return *this;
    }
    net_decimal operator++(int) {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    net_decimal operator-() { return 0 - *this; }
    friend net_decimal operator-(const net_decimal &fst, const net_decimal &snd) {
        net_decimal ans;
        ans.sgn = dec_frac_sub(ans.num, ans.den, fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn);
        ans.red = dec_frac1(ans.num, ans.den) || dec_frac0(ans.num, ans.den) || dec_is_zero(ans.den);
        return ans;
    }
    net_decimal &operator-=(const net_decimal &snd) {
        *this = *this - snd;
        return *this;
    }
    callback_dec_arg friend arg operator-=(arg &fst, const net_decimal &snd) { return fst -= snd.to_float(); }
    net_decimal &operator--() {
        *this -= 1;
        return *this;
    }
    net_decimal operator--(int) {
        auto tmp = *this;
        --(*this);
        return tmp;
    }

    friend net_decimal operator*(const net_decimal &fst, const net_decimal &snd) {
        net_decimal ans;
        ans.sgn = dec_frac_mul(ans.num, ans.den, fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn);
        ans.red = dec_frac1(ans.num, ans.den) || dec_frac0(ans.num, ans.den) || dec_is_zero(ans.den);
        return ans;
    }
    net_decimal &operator*=(const net_decimal &snd) {
        *this = *this * snd;
        return *this;
    }
    callback_dec_arg friend arg operator*=(arg &fst, const net_decimal &snd) { return fst *= snd.to_float(); }

    friend net_decimal operator/(const net_decimal &divd, const net_decimal &divr) {
        net_decimal ans;
        ans.sgn = dec_frac_div(ans.num, ans.den, divd.num, divd.den, divd.sgn, divr.num, divr.den, divr.sgn);
        ans.red = dec_frac1(ans.num, ans.den) || dec_frac0(ans.num, ans.den) || dec_is_zero(ans.den);
        return ans;
    }
    net_decimal &operator/=(const net_decimal &divr) {
        *this = *this / divr;
        return *this;
    }
    callback_dec_arg friend arg operator/=(arg &divd, net_decimal &divr) {
        if (dec_is_zero(divr.den)) return divd /= dec2f(divr.sgn, divr.num);
        return divd /= divr.to_float();
    }

    friend net_decimal operator%(const net_decimal &divd, const net_decimal &divr) {
        net_decimal ans;
        dec_frac_rem(ans.num, ans.den, divd.num, divd.den, divr.num, divr.den);
        if (divr.modulus_mode && divd.sgn != divr.sgn) {
            if (divr.sgn) ans += divr;
            else ans = divr - ans;
            return ans;
        }
        ans.sgn = divd.sgn;
        return ans;
    }
    net_decimal &operator%=(const net_decimal &divr) {
        *this = *this % divr;
        return *this;
    }
    callback_dec_arg friend arg operator%=(arg &divd, const net_decimal &divr) {
        auto divr_int = divr.to_int();
        divd         %= divr_int;
        if (divr.modulus_mode && divr.sgn) divd += divr_int;
        if (divd < 0) divd *= -1;
        return divd;
    }

    friend net_decimal operator<<(const net_decimal &src, const net_decimal &bit) {
        auto ans = src;
        ans <<= bit;
        return ans;
    }
    net_decimal &operator<<=(const net_decimal &bit) {
        if (dec_is_zero(bit.den)) {
            if (bit.num.ft.length) return *this;
            reduct();
            if (dec_is_zero(den)) dec_bit_lsh(num, bit.num, binary_bit_set);
            return *this;
        }
        net_decimal bit_cnt;
        if (dec_frac_bit_verify(bit_cnt.num, bit.num, bit.den)) return *this <<= bit_cnt;
        return *this;
    }
    callback_dec_arg friend arg operator<<=(arg &src, const net_decimal &bit) {
        net_decimal_data bit_cnt;
        if (!dec_frac_bit_verify(bit_cnt, bit.num, bit.den)) return src;
        if (dec_is_zero(den) && bit_cnt.it.length) src <<= bit_cnt.it[0];
        return src;
    }

    friend net_decimal operator>>(const net_decimal &src, const net_decimal &bit) {
        auto ans = src;
        ans >>= bit;
        return ans;
    }
    net_decimal &operator>>=(const net_decimal &bit) {
        if (dec_is_zero(bit.den)) {
            if (bit.num.ft.length) return *this;
            reduct();
            if (dec_is_zero(den)) dec_bit_rsh(num, bit.num, binary_bit_set);
            return *this;
        }
        net_decimal bit_cnt;
        if (dec_frac_bit_verify(bit_cnt.num, bit.num, bit.den)) return *this >>= bit_cnt;
        return *this;
    }
    callback_dec_arg friend arg operator>>=(arg &src, const net_decimal &bit) {
        net_decimal_data bit_cnt;
        if (!dec_frac_bit_verify(bit_cnt, bit.num, bit.den)) return src;
        if (dec_is_zero(den) && bit_cnt.it.length) src >>= bit_cnt.it[0];
        return src;
    }

    friend net_decimal operator&(const net_decimal &fst, const net_decimal &snd) { return binary_operator<NEUNET_DEC_BIN_AND>(fst, snd); }
    net_decimal &operator&=(const net_decimal &snd) {
        *this = *this & snd;
        return *this;
    }
    callback_dec_arg friend arg operator&=(arg &fst, const net_decimal &snd) {
        if (dec_is_zero(snd.den)) {
            if (snd.num.ft.length) return fst;
            if (snd.num.it.length) return fst &= snd.num.it[0];
            return 0;
        }
        auto tmp = snd.num,
             opt = dec_rem(tmp, snd.den);
        if (dec_is_zero(tmp) && !tmp.ft.length) return fst &= opt.it[0];
        return fst;
    }

    friend net_decimal operator|(const net_decimal &fst, const net_decimal &snd) { return binary_operator<NEUNET_DEC_BIN_OR>(fst, snd); }
    net_decimal &operator|=(const net_decimal &snd) {
        *this = *this | snd;
        return *this;
    }
    callback_dec_arg friend arg operator|=(arg &fst, const net_decimal &snd) {
        if (dec_is_zero(snd.den)) {
            if (snd.num.ft.length) return fst;
            if (snd.num.it.length) return fst |= snd.num.it[0];
            return fst;
        }
        auto tmp = snd.num,
             opt = dec_rem(tmp, snd.den);
        if (dec_is_zero(tmp) && !tmp.ft.length) {
            if (opt.it.length) return fst |= opt.it[0];
            return fst;
        }
        return fst;
    }

    friend net_decimal operator^(const net_decimal &fst, const net_decimal &snd) { return binary_operator<NEUNET_DEC_BIN_XOR>(fst, snd); }
    net_decimal &operator^=(const net_decimal &snd) {
        *this = *this ^ snd;
        return *this;
    }
    callback_dec_arg friend arg operator^=(arg &fst, const net_decimal &snd) {
        if (dec_is_zero(snd.den)) {
            if (snd.num.ft.length) return fst;
            if (snd.num.it.length) return fst ^= snd.num.it[0];
            return fst ^= 0;
        }
        auto tmp = snd.num,
             opt = dec_rem(tmp, snd.den);
        if (dec_is_zero(tmp) && !tmp.ft.length) {
            if (opt.it.length) return fst ^= opt.it[0];
            return fst ^= 0;
        }
        return fst;
    }

    net_decimal operator~() {
        if (dec_is_zero(den)) {
            if (num.ft.length) return *this;
            dec_bit_not(num, false, binary_bit_set, sgn);
            return *this;
        }
        auto rem = num,
             opt = dec_rem(rem, den);
        if (dec_is_zero(rem)) {
            num = std::move(opt);
            dec_bit_not(num, false, binary_bit_set, sgn);
        }
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const net_decimal &src) {
        if (src.sgn) os << '-';
        if (dec_is_zero(src.den)) return os << src.num;
        return os << dec_div(src.num, src.den, src.infinite_precision);
    }
};

net_decimal operator""_d(const char *src, size_t len) { return net_decimal(std::string(src)); }
net_decimal operator""_d(long double src) { return net_decimal(src); }
net_decimal operator""_d(uint64_t src) { return net_decimal(src); }

NEUNET_END

_STD_BEGIN

neunet::net_decimal log(const neunet::net_decimal &_Xx) { return _Xx.ln(); }

neunet::net_decimal abs(const neunet::net_decimal &_Xx) { return _Xx.abs(); }

neunet::net_decimal exp(const neunet::net_decimal &_Xx) { return _Xx.exp(); }

neunet::net_decimal sin(const neunet::net_decimal &_Xx) { return _Xx.sin(); }

neunet::net_decimal cos(const neunet::net_decimal &_Xx) { return _Xx.cos(); }

neunet::net_decimal pow(const neunet::net_decimal &_Xx, const neunet::net_decimal &_Yx) { return _Xx.pow(_Yx); }

_STD_END