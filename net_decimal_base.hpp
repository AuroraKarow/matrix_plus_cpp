NEUNET_BEGIN

int dec_comp(const net_decimal_data &fst, bool fst_sgn, const net_decimal_data &snd, bool snd_sgn) {
    if (fst_sgn == snd_sgn) return dec_comp(fst_sgn, dec_comp(fst, snd));
    if (fst_sgn) return NEUNET_DEC_CMP_LES;
    return NEUNET_DEC_CMP_GTR;
}

net_decimal_data dec_add(bool &ans_sgn, const net_decimal_data &fst, bool fst_sgn, const net_decimal_data &snd, bool snd_sgn) {
    if (dec_is_zero(fst)) {
        ans_sgn = snd_sgn;
        return snd;
    }
    if (dec_is_zero(snd)) {
        ans_sgn = fst_sgn;
        return fst;
    }
    ans_sgn = false;
    if (fst_sgn != snd_sgn) {
        auto cmp_res = dec_comp(snd, fst);
        if (cmp_res == NEUNET_DEC_CMP_GTR) {
            if (snd_sgn) ans_sgn = true;
            return dec_add(snd, fst, true);
        }
        if (cmp_res == NEUNET_DEC_CMP_LES) {
            if (fst_sgn) ans_sgn = true;
            return dec_add(fst, snd, true);
        }
        return {};
    }
    if (fst_sgn) ans_sgn = true;
    return dec_add(fst, snd);
}

net_decimal_data dec_sub(bool &ans_sgn, const net_decimal_data &minu, bool minu_sgn, const net_decimal_data &subt, bool subt_sgn) {
    if (dec_is_zero(subt)) {
        ans_sgn = minu_sgn;
        return minu;
    }
    if (dec_is_zero(minu)) {
        ans_sgn = !minu_sgn;
        return subt;
    }
    ans_sgn = false;
    return dec_add(ans_sgn, minu, minu_sgn, subt, !subt_sgn);
}

net_decimal_data dec_mul(bool &ans_sgn, const net_decimal_data &fst, bool fst_sgn, const net_decimal_data &snd, bool snd_sgn) {
    ans_sgn = fst_sgn != snd_sgn;
    if (dec_is_zero(fst) || dec_is_zero(snd)) return {};
    if (dec_is_one(fst)) return snd;
    if (dec_is_one(snd)) return fst;
    return dec_mul(fst, snd);
}

net_decimal_data dec_div(bool &ans_sgn, const net_decimal_data &divd, bool divd_sgn, const net_decimal_data &divr, bool divr_sgn, uint64_t prec) {
    ans_sgn = divd_sgn != divr_sgn;
    if (dec_is_zero(divd)) return {};
    if (dec_is_zero(divr)) return divd;
    return dec_div(divd, divr, prec);
}

// kp
net_decimal_data dec_bit_k0(bool sgn = false) { return dec_init(sgn, 0); }
// k, k1
net_decimal_data dec_bit_k1(bool sgn = false) { return dec_init(sgn, 1); }

/** 
 * @brief Left shift
 * @param src   [IO]    source
 * @param bit   [In]    shift count
 * @param sta   [In]    bit stability
 */
void dec_bit_lsh(net_decimal_data &src, uint64_t bit = 1, bool sta = false){
    if (src.ft.length || !src.it.length) return;
    if (!(src.it.length - 1) && log(src.it[0]) + bit < 64) { src.it[0] <<= bit; return; }
    uint64_t b = bit % NEUNET_DEC_DIG_MAX,
             B = bit / NEUNET_DEC_DIG_MAX,
             l = src.it.length,
             T = 0,
             t;
    if (!sta && (B || src.it[l - 1] >= (NEUNET_DEC_SEG_MAX >> b))) src.it.init(l + (bit / 63) + 1);
    if (src.it[l - 1] >= (NEUNET_DEC_SEG_MAX >> b)) l++;
    // std::cout << "S" << std::endl;
    if (b){
        uint64_t bas = NEUNET_DEC_SEG_MAX >> b;
        for(uint64_t m = 0; m < l; m++) {
            t = src.it[m] / bas;
            src.it[m] %= bas;
            src.it[m] <<= b;
            src.it[m] += T;
            T = t;
        }
    }
    // std::cout << "A" << std::endl;
    uint64_t Bas = NEUNET_DEC_SEG_MAX >> NEUNET_DEC_DIG_MAX;
    while (B--){
        T = 0;
        for(uint64_t m = 0; m < src.it.length; m++) {
            t = src.it[m] / Bas;
            src.it[m] %= Bas;
            src.it[m] <<= NEUNET_DEC_DIG_MAX;
            src.it[m] += T;
            T = t;
        }
    }
    l = src.it.length;
    // std::cout << "B" << std::endl;
    if (!sta) {
        while(!src.it[--l] && l > 0) continue;
        if (!src.it[l]){
            src.it.reset();
            return;
        }
        src.it = src.it.sub_set(0, l);
    }
}
void dec_bit_lsh(net_decimal_data &src, net_decimal_data bit, bool sta = false){
    if (src.ft.length || !src.it.length || bit.ft.length || !bit.it.length) return;
    if (!(src.it.length - 1) && !(bit.it.length - 1) && log(src.it[0]) + bit.it[0] < 64) { src.it[0] <<= bit.it[0]; return; }
    bool sgn = false;
    if (!(bit.it.length - 1)){
        dec_bit_lsh(src, bit.it[0]);
        return;
    }
    net_decimal_data b = dec_init(sgn, NEUNET_DEC_SEG_MAX),
                     B,
                     z = dec_bit_k0(),
                     o = dec_bit_k1();
    B = dec_div(src, b, 0);
    B.ft.reset();
    uint64_t l = src.it[0],
             T = 0,
             t;
    for (auto i = z; dec_comp(i, B) == NEUNET_DEC_CMP_LES; i = dec_add(i, o)){ dec_bit_lsh(src, NEUNET_DEC_SEG_MAX); }
    dec_bit_lsh(src, bit.it[0]);
}

/**
 * @brief Right shift
 * @param src   [IO]    source
 * @param bit   [In]    shift count
 * @param sta   [In]    bit stability
 */
void dec_bit_rsh(net_decimal_data &src, uint64_t bit = 1, bool sta = false){
    if (src.ft.length || !src.it.length) return;
    if (!(src.it.length - 1)) { src.it[0] >>= bit; return; }
    uint64_t b = bit % NEUNET_DEC_DIG_MAX,
             B = bit / NEUNET_DEC_DIG_MAX,
             l = src.it.length,
             T = 0,
             t;
    // std::cout << "S" << std::endl;
    if (b){
        uint64_t bas = NEUNET_DEC_SEG_MAX >> b,
                 div = 1 << b;
        for(uint64_t m = l; m > 0; m--) {
            t = (src.it[m - 1] % div) * bas;
            src.it[m - 1] >>= b;
            src.it[m - 1] += T;
            T = t;
        }
    }
    // std::cout << "A" << std::endl;
    uint64_t Bas = NEUNET_DEC_SEG_MAX >> NEUNET_DEC_DIG_MAX,
             Div = 1 << NEUNET_DEC_DIG_MAX;
    while (B--){
        T = 0;
        for(uint64_t m = l; m > 0; m--) {
            t = (src.it[m - 1] % Div) * Bas;
            src.it[m - 1] >>= NEUNET_DEC_DIG_MAX;
            src.it[m - 1] += T;
            T = t;
        }
    }
    // std::cout << "B" << std::endl;
    if (!sta) {
        while(!src.it[--l] && l > 0) continue;
        if (!src.it[l]){
            src.it.reset();
            return;
        }
        src.it = src.it.sub_set(0, l);
    }
}
void dec_bit_rsh(net_decimal_data &src, net_decimal_data bit, bool sta = false){
    if (src.ft.length || !src.it.length || bit.ft.length || !bit.it.length) return;
    if (!(src.it.length - 1)){
        if (!(bit.it.length - 1)) {
            src.it[0] >>= bit.it[0];
            return; 
        }
        else {
            src.it.reset();
            return; 
        }
    }
    bool sgn = true;
    if (!(bit.it.length - 1)){
        dec_bit_rsh(src, bit.it[0]);
        return;
    }
    net_decimal_data b = dec_init(sgn, NEUNET_DEC_SEG_MAX),
             B,
             z = dec_bit_k0(),
             o = dec_bit_k1();
    B = dec_div(src, b, 0);
    B.ft.reset();
    uint64_t l = src.it[0],
             T = 0,
             t;
    for (auto i = z; dec_comp(i, B) == NEUNET_DEC_CMP_LES; i = dec_add(i, o)){ dec_bit_rsh(src, NEUNET_DEC_SEG_MAX); }
    dec_bit_rsh(src, bit.it[0]);
}

void dec_bit_lsh1(net_decimal_data &src, bool sta = false){
    if (src.ft.length || !src.it.length) return;
    if (!sta && src.it[src.it.length - 1] >= NEUNET_DEC_BIT_BAS) src.it.init(src.it.length + 1);
    int t = 0;
    for(uint64_t m = 0; m < src.it.length; m++) if(src.it[m] >= NEUNET_DEC_BIT_BAS) {
        src.it[m] -= NEUNET_DEC_BIT_BAS;
        src.it[m] <<= 1;
        src.it[m] += t;
        t = 1;
    } else {
        src.it[m] <<= 1;
        src.it[m] += t;
        t = 0;
    }
}

void dec_bit_rsh1(net_decimal_data &src, bool sta = false){
    if (src.ft.length || !src.it.length) return;
    //if (!(src.it.length - 1)) src.it[0] >>= 1; return;
    uint64_t t = 0;
    for(uint64_t m = src.it.length; m > 0; m--) if(src.it[m - 1] % 2 == 1){ 
        src.it[m - 1] >>= 1;
        src.it[m - 1] += t;
        t = NEUNET_DEC_BIT_BAS;
    } else{
        src.it[m - 1] >>= 1;
        src.it[m - 1] += t;
        t = 0;
    }
    if (!sta && !src.it[src.it.length - 1]) src.it.length - 1 ? src.it[src.it.length - 2] < NEUNET_DEC_BIT_TOP ? src.it.init(src.it.length - 1, true) : src.it.init(src.it.length, true) : src.it.init(src.it.length - 1, true);
}

uint64_t dec_bit_cnt(net_decimal_data &src){
    uint64_t b = 0;
    auto K = src,
         k = dec_bit_k0();
    while(dec_comp(k, K) != NEUNET_DEC_CMP_EQL){ 
        dec_bit_rsh1(K);
        b++;
    }
    return b;
}

/**
 * @brief Binary not
 * @param src       [IO]    source
 * @param comple    [In]    one's complement
 * @param bit       [In]    lower bit count, minimum is the bit count of source value
 * @param sgn       [In]    Sign of source
 */
void dec_bit_not(net_decimal_data &src, bool comple = false, uint64_t bit = 0, bool sgn = false) {
    auto k0 = dec_bit_k0(),
         k1 = dec_bit_k1(),
         k  = k0,
         K  = k1;
    auto b  = 0;
    auto k_sgn = false;
    if (dec_comp(src, k0) == NEUNET_DEC_CMP_EQL && !bit) {src = std::move(k1); return;}
    if (comple && !sgn) return;
    while(dec_comp(src, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (src.it.length || !(src.it.length - 1) && src.it[0] != 0) { if (!(src.it[0] % 2)) k = dec_add(k_sgn, k, false, k1, false); }
        else if (!src.it.length) k = dec_add(k_sgn, k, false, k1, false);
        dec_bit_lsh1(k1);
        dec_bit_rsh1(src);
        b++;
    }
    if (comple) {
        k = dec_add(k_sgn, k, false, K, false);
        if (sgn) { if (dec_comp(k, k1) != NEUNET_DEC_CMP_LES) { dec_bit_lsh1(k1); } k = dec_add(k_sgn, k, false, k1, false); }
    }
    else if (!sgn) k = dec_add(k_sgn, k, false, k1, false);
    /*
    while(dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
        if(!(D1.it[0] % 2)){ k = dec_add(sgn, k, sgn, kp, sgn); }
        dec_bit_lsh_one(kp);
        dec_bit_rsh_one(D1);
    }
    */
    src = std::move(k);
}

/**
 * @brief Binary and
 * @param fst   [In]    first source
 * @param snd   [In]    second source
 * @param sgn1  [In]    first source sign
 * @param sgn2  [In]    second source sign
 * @param bit   [In]    lower bit count, minimum is the bit count of sources value
 */
net_decimal_data dec_bit_and(const net_decimal_data &fst, const net_decimal_data &snd, bool sgn1 = false, bool sgn2 = false, uint64_t bit = 0) {
    auto D  = fst,
         d  = snd,
         k0 = dec_bit_k0(),
         k1 = dec_bit_k1(),
         k  = k0;
    auto b  = 0;
    auto k_sgn = false;
    if (!bit) {
        auto a = dec_bit_cnt(D),
             b = dec_bit_cnt(d);
        bit = std::max(a, b);
    }
    if (sgn1) dec_bit_not(D, true, bit, true);
    if (sgn2) dec_bit_not(d, true, bit, true);
    while(dec_comp(d, k0) != NEUNET_DEC_CMP_EQL || dec_comp(D, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (!(!D.it.length || !d.it.length)) if (D.it[0] % 2 && d.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false);
        dec_bit_lsh1(k1);
        dec_bit_rsh1(D);
        dec_bit_rsh1(d);
        b++;
    }
    /*
    while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
        if(D1.it[0] % 2 && d1.it[0] % 2){ k = dec_add(sgn, k, sgn, kp, sgn); }
        dec_bit_lsh_one(kp);
        dec_bit_rsh_one(D1);
        dec_bit_rsh_one(d1);
    }
    */
    return k;
}

/**
 * @brief Binary or
 * @param fst   [In]    first source
 * @param snd   [In]    second source
 * @param sgn1  [In]    first source sign
 * @param sgn2  [In]    second source sign
 * @param bit   [In]    lower bit count, minimum is the bit count of sources value
 */
net_decimal_data dec_bit_or(const net_decimal_data &fst, const net_decimal_data &snd, bool sgn1 = false, bool sgn2 = false, uint64_t bit = 0) {
    auto D  = fst,
         d  = snd,
         k0 = dec_bit_k0(),
         k1 = dec_bit_k1(),
         k  = k0;
    auto b  = 0;
    auto k_sgn = false;
    if (!bit) {
        auto a = dec_bit_cnt(D),
             b = dec_bit_cnt(d);
        bit = a > b ? a : b;
    }
    if (sgn1) dec_bit_not(D, true, bit, true);
    if (sgn2) dec_bit_not(d, true, bit, true);
    while(dec_comp(d, k0) != NEUNET_DEC_CMP_EQL || dec_comp(D, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (!(!D.it.length || !d.it.length)) { if((D.it[0] % 2) || (d.it[0] % 2)) k = dec_add(k_sgn, k, true, k1, true); }
        else if (!D.it.length) { if (d.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false); }
        else if (!d.it.length) { if (D.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false); }
        dec_bit_lsh1(k1);
        dec_bit_rsh1(D);
        dec_bit_rsh1(d);
        b++;
    }
    /*
    while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
        if(!(!(D1.it[0] % 2) && !(d1.it[0] % 2))){ k = dec_add(sgn, k, sgn, kp, sgn); }
        dec_bit_lsh_one(kp);
        dec_bit_rsh_one(D1);
        dec_bit_rsh_one(d1);
    }
    */
    return k;
}

/**
 * @brief Binary xor
 * @param fst   [In]    first source
 * @param snd   [In]    second source
 * @param sgn1  [In]    first source sign
 * @param sgn2  [In]    second source sign
 * @param bit   [In]    lower bit count, minimum is the bit count of sources value
 */
net_decimal_data dec_bit_xor(const net_decimal_data &fst, const net_decimal_data &snd, bool sgn1 = false, bool sgn2 = false, uint64_t bit = 0) {
    auto k0 = dec_bit_k0(),
         k1 = dec_bit_k1(),
         k  = k0,
         D  = fst,
         d  = snd;
    auto b  = 0;
    auto k_sgn = false;
    if (!bit) {
        auto a = dec_bit_cnt(D),
             b = dec_bit_cnt(d);
        bit = a > b ? a : b;
    }
    if (sgn1) dec_bit_not(D, true, bit, true);
    if (sgn2) dec_bit_not(d, true, bit, true);
    while(dec_comp(d, k0) != NEUNET_DEC_CMP_EQL || dec_comp(D, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (!(!D.it.length || !d.it.length)) { if((D.it[0] % 2 && !(d.it[0] % 2)) || (!(D.it[0] % 2) && d.it[0] % 2)) k = dec_add(k_sgn, k, true, k1, true); }
        else if (!D.it.length) { if (d.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false); }
        else if (!d.it.length) { if (D.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false); }
        dec_bit_lsh1(k1);
        dec_bit_rsh1(D);
        dec_bit_rsh1(d);
        b++;
    }
    /*
    while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
        if((D1.it[0] % 2 && !(d1.it[0] % 2)) || (!(D1.it[0] % 2) && d1.it[0] % 2)){ k = dec_add(sgn, k, sgn, kp, sgn); }
        dec_bit_lsh_one(kp);
        dec_bit_rsh_one(D1);
        dec_bit_rsh_one(d1);
    }
    */
    return k;
}

bool dec_frac_is_zero(const net_decimal_data &num, const net_decimal_data &den) { return dec_is_zero(den) && dec_is_zero(num); }

bool dec_frac_is_one(const net_decimal_data &num, const net_decimal_data &den) { return dec_is_zero(den) && dec_is_one(num); }

long double dec_frac2f(bool sgn, net_decimal_data &num, net_decimal_data &den) { return dec2f(sgn, dec_div(num, den, 2)); }

int64_t dec_frac2i(bool sgn, net_decimal_data &num, net_decimal_data &den) {
    if (dec_is_zero(den)) return dec2i(sgn, num);
    return dec_frac2f(sgn, num, den);
}

int dec_frac_comp(const net_decimal_data &fst_num, const net_decimal_data &fst_den, bool fst_sgn, const net_decimal_data &snd_num, const net_decimal_data &snd_den, bool snd_sgn) {
    if (dec_comp(fst_den, snd_den) == NEUNET_DEC_CMP_EQL) return dec_comp(fst_num, fst_sgn, snd_num, snd_sgn);
    auto fst_is_frac = !dec_is_zero(fst_den),
         snd_is_frac = !dec_is_zero(snd_den);
    if (fst_is_frac && snd_is_frac) return dec_comp(dec_mul(fst_num, snd_den), fst_sgn, dec_mul(snd_num, fst_den), snd_sgn);
    if (fst_is_frac) return dec_comp(fst_num, fst_sgn, dec_mul(snd_num, fst_den), snd_sgn);
    return dec_comp(dec_mul(fst_num, snd_den), fst_sgn, snd_num, snd_sgn);
}

bool dec_frac_add(net_decimal_data &ans_num, net_decimal_data &ans_den, const net_decimal_data &fst_num, const net_decimal_data &fst_den, bool fst_sgn, const net_decimal_data &snd_num, const net_decimal_data &snd_den, bool snd_sgn) {
    if (dec_frac_is_zero(fst_num, fst_den)) {
        ans_num = snd_num;
        ans_den = snd_den;
        return snd_sgn;
    }
    if (dec_frac_is_zero(snd_num, snd_den)) {
        ans_num = fst_num;
        ans_den = fst_den;
        return snd_sgn;
    }
    auto ans_sgn = false;
    if (dec_comp(fst_den, snd_den) == NEUNET_DEC_CMP_EQL) {
        ans_num = dec_add(ans_sgn, fst_num, fst_sgn, snd_num, snd_sgn);
        ans_den = fst_den;
        return ans_sgn;
    }
    auto fst_is_frac = !dec_is_zero(fst_den),
         snd_is_frac = !dec_is_zero(snd_den);
    if (fst_is_frac && snd_is_frac) {
        ans_num = dec_add(ans_sgn, dec_mul(fst_num, snd_den), fst_sgn, dec_mul(snd_num, fst_den), snd_sgn);
        ans_den = dec_mul(fst_den, snd_den);
        return ans_sgn;
    }
    if (fst_is_frac) {
        ans_num = dec_add(ans_sgn, fst_num, fst_sgn, dec_mul(snd_num, fst_den), snd_sgn);
        ans_den = fst_den;
        return ans_sgn;
    }
    ans_num = dec_add(ans_sgn, dec_mul(fst_num, snd_den), fst_sgn, snd_num, snd_sgn);
    ans_den = snd_den;
    return ans_sgn;
}

bool dec_frac_sub(net_decimal_data &ans_num, net_decimal_data &ans_den, const net_decimal_data &minu_num, const net_decimal_data &minu_den, bool minu_sgn, const net_decimal_data &subt_num, const net_decimal_data &subt_den, bool subt_sgn) {
    if (dec_frac_is_zero(subt_num, subt_den)) {
        ans_num = minu_num;
        ans_den = minu_den;
        return minu_sgn;
    }
    if (dec_frac_is_zero(minu_num, minu_den)) {
        ans_num = subt_num;
        ans_den = subt_den;
        return !subt_sgn;
    }
    return dec_frac_add(ans_num, ans_den, minu_num, minu_den, minu_sgn, subt_num, subt_den, !subt_sgn);
}

bool dec_frac_mul(net_decimal_data &ans_num, net_decimal_data &ans_den, const net_decimal_data &fst_num, const net_decimal_data &fst_den, bool fst_sgn, const net_decimal_data &snd_num, const net_decimal_data &snd_den, bool snd_sgn) {
    if (dec_frac_is_zero(fst_num, fst_den) || dec_frac_is_zero(snd_num, snd_den)) {
        ans_num.reset();
        ans_den.reset();
        return false;
    }
    auto ans_sgn = fst_sgn != snd_sgn;
    if (dec_frac_is_one(snd_num, snd_den)) {
        ans_num = fst_num;
        ans_den = fst_den;
        return ans_sgn;
    }
    if (dec_frac_is_one(fst_num, fst_den)) {
        ans_num = snd_num;
        ans_den = snd_den;
        return ans_sgn;
    }
    auto fst_is_frac = !dec_is_zero(fst_den),
         snd_is_frac = !dec_is_zero(snd_den);
    ans_num          = dec_mul(fst_num, snd_num);
    if (!(fst_is_frac || snd_is_frac)) {
        ans_den.reset();
        return ans_sgn;
    }
    if (fst_is_frac && snd_is_frac) {
        ans_den = dec_mul(fst_den, snd_den);
        return ans_sgn;
    }
    if (snd_is_frac) ans_den = snd_den;
    else ans_den = fst_den;
    return ans_sgn;
}

bool dec_frac_div(net_decimal_data &ans_num, net_decimal_data &ans_den, const net_decimal_data &divd_num, const net_decimal_data &divd_den, bool divd_sgn, const net_decimal_data &divr_num, const net_decimal_data &divr_den, bool divr_sgn) {
    assert(!dec_is_zero(divr_num));
    if (dec_frac_is_zero(divd_num, divd_num)) {
        ans_num.reset();
        ans_den.reset();
        return false;
    }
    if (dec_frac_is_one(divr_num, divr_den)) {
        ans_num = divd_num;
        ans_den = divd_den;
        return divd_sgn;
    }
    if (dec_frac_is_one(divd_num, divd_den)) {
        ans_num = divr_den;
        ans_den = divr_num;
        if (dec_is_zero(ans_num)) ans_num.it = {1};
        return divr_sgn;
    }
    if (dec_is_zero(divd_den) && dec_is_zero(divr_den)) {
        ans_num = divd_num;
        ans_den = divr_num;
        return divd_sgn != divr_sgn;
    }
    if (dec_is_zero(divr_den)) {
        auto sgn = false;
        return dec_frac_mul(ans_num, ans_den, divd_num, divd_den, divd_sgn, dec_init(sgn, 1), divr_num, divr_sgn);
    }
    return dec_frac_mul(ans_num, ans_den, divd_num, divd_den, divd_sgn, divr_den, divr_num, divr_sgn);
}

// return quotient, unsigned
net_decimal_data dec_frac_rem(net_decimal_data &ans_num, net_decimal_data &ans_den, const net_decimal_data &divd_num, const net_decimal_data &divd_den, const net_decimal_data &divr_num, const net_decimal_data &divr_den) {
    dec_frac_div(ans_num, ans_den, divd_num, divd_den, false, divr_num, divr_den, false);
    auto ans = dec_rem(ans_num, ans_den);
    ans_den.reset();
    if (dec_is_zero(divd_den) && dec_is_zero(divr_den)) return ans;
    if (!(dec_is_zero(divd_den) || dec_is_zero(divr_den))) {
        ans_num = dec_div(ans_num, dec_mul(divd_den, divr_den), 0);
        return ans;
    }
    if (dec_is_zero(divd_den)) {
        ans_num = dec_div(ans_num, divr_den, 0);
        return ans;
    }
    ans_num = dec_div(ans_num, divd_den, 0);
    return ans;
}

bool dec_frac0(net_decimal_data &num, net_decimal_data &den) { if (dec_is_zero(num)) {
    den.reset();
    return true;
} return false; }

bool dec_frac1(net_decimal_data &num, net_decimal_data &den) { if (!dec_is_zero(num) && dec_comp(num, den) == NEUNET_DEC_CMP_EQL) {
    den.reset();
    num.reset();
    num.it = {1};
    return true;
} return false; }

void dec_frac_norm(net_decimal_data &num, net_decimal_data &den, bool total = false) {
    if (!(num.ft.length || den.ft.length)) return;
    auto lsh_cnt = 0ull;
    if (total) lsh_cnt = (std::max)(dec_dig_cnt(num, false), dec_dig_cnt(den, false));
    else lsh_cnt = (std::max)(den.ft.length, num.ft.length) * NEUNET_DEC_DIG_MAX;
    num = dec_e10_lsh(num, lsh_cnt);
    den = dec_e10_lsh(den, lsh_cnt);
}

void dec_frac_red(net_decimal_data &num, net_decimal_data &den) {
    if (dec_is_zero(den)) return;
    // dec_frac_norm(num, den, true);
    auto fst = num,
         snd = den;
    auto gcd = dec_gcd(fst, snd);
    if (gcd) {
        num = dec_div(num, snd, 0);
        den = dec_div(den, snd, 0);
    } else {
        num = dec_div(num, fst, 0);
        den = dec_div(den, fst, 0);
    }
    if (dec_is_one(den)) den.reset();
}

net_decimal_data dec_frac_int_part(net_decimal_data &num, net_decimal_data &den) {
    if (!dec_is_zero(den)) return dec_rem(num, den);
    net_decimal_data ans;
    num = dec_float_part(ans, num);
    return ans;
}

bool dec_frac_bit_verify(net_decimal_data &ans, const net_decimal_data &num, const net_decimal_data &den) {
    if (dec_is_zero(den) && num.ft.length) return false;
    auto tmp = num;
    ans      = dec_rem(tmp, den);
    return dec_is_zero(tmp);
}

bool dec_series_check(const net_decimal_data &divd, const net_decimal_data &divr, uint64_t prec) { return dec_comp(dec_e10_lsh(divd, prec), divr) == NEUNET_DEC_CMP_GTR; }

// NOT thread safety
struct dec_series_prec {
private: inline static net_set<net_decimal_data> prec {};

public: static net_decimal_data value(uint64_t precision) {
    if (!precision) return dec_e10(0);
    --precision;
    if (precision >= prec.length) {
        auto len = prec.length;
        if (!len) len = 128;
        while (len < precision) len <<= 1;
        prec.init(len);
    }
    if (dec_is_zero(prec[precision])) prec[precision] = dec_e10(precision + 1);
    return prec[precision];
}
};

// NOT thread safety
struct dec_series_pi {
private:
    inline static bool on  = false,
                       sgn = false;

    inline static std::atomic_uint64_t prec = 0;

    inline static net_decimal_data p, q, a, // +2
                                   b,       // +4
                                   c, d,    // +8
                                   o, s, dec_2, dec_4, dec_8;
    
    static void init() {
        if (on) return;
        p = dec_init(sgn, 0);
        q = dec_init(sgn, 1);
        a = dec_init(sgn, 0.25);
        b = dec_init(sgn, 2);
        c = dec_init(sgn, 5);
        d = dec_init(sgn, 6);
        o = dec_init(sgn, 0.0625);
        s = q;
        dec_2 = dec_init(sgn, 2);
        dec_4 = dec_init(sgn, 4);
        dec_8 = dec_init(sgn, 8);
        on = true;
    }

    static void run(uint64_t precision) { while (prec < precision) {
        init();
        auto u = dec_init(sgn, 1),
             v = a;
        u = dec_add(dec_mul(b, u), v, true);
        v = dec_mul(v, b);
        u = dec_add(dec_mul(c, u), v, true);
        v = dec_mul(v, c);
        u = dec_add(dec_mul(d, u), v, true);
        v = dec_mul(v, d);
        u = dec_mul(u, s);
        if (dec_comp(q, v) == NEUNET_DEC_CMP_EQL) p = dec_add(p, u);
        else {
            p = dec_add(dec_mul(v, p), dec_mul(u, q));
            q = dec_mul(v, q);
        }
        s = dec_mul(s, o);
        a = dec_add(a, dec_2);
        b = dec_add(b, dec_4);
        c = dec_add(c, dec_8);
        d = dec_add(d, dec_8);
        ++prec;
    } }

public: static net_decimal_data value(uint64_t precision) {
    run(precision);
    return dec_div(p, q, precision);
}
};

// NOT thread safety
struct dec_series_ln4 {
private:
    inline static bool sgn = false;

    inline static uint64_t prec = 0;

    inline static net_decimal_data num_form = dec_init(sgn, "0.36"),
                                   den_form = dec_init(sgn, "2"),
                                   num_term = dec_init(sgn, "0.6"),
                                   den_term = dec_init(sgn, "1"),
                                   
                                   den = den_term,
                                   num,
                                   ans;

public: static net_decimal_data value(uint64_t precision) {
    if (precision <= prec) return ans;
    // auto iter_cnt = 0ull;
    do {
        if (dec_comp(den, den_term) == NEUNET_DEC_CMP_EQL) num = dec_add(num, num_term);
        else {
            num = dec_add(dec_mul(den_term, num), dec_mul(num_term, den));
            den = dec_mul(den, den_term);
        }
        // ++iter_cnt;
        num_term = dec_mul(num_term, num_form);
        den_term = dec_add(den_term, den_form);
    } while (dec_series_check(num_term, den_term, prec));
    prec = precision;
    ans  = dec_mul(den_form, dec_div(num, den, prec));
    // std::cout << "next ln4 - " << iter_cnt << std::endl;
    return ans;
}
};

net_decimal_data dec_series_ln_2(bool &ans_sgn, const net_decimal_data &src, uint64_t prec) {
    // auto iter_cnt = 0ull;
    auto form_sgn = false,
         term_sgn = false,
         coef_sgn = false;

    auto den_form = dec_init(ans_sgn, "1"),
         num_form = dec_sub(form_sgn, src, false, den_form, false),
         den_term = den_form,
         num_term = num_form,
         
         num = dec_init(ans_sgn, 0),
         den = den_form;
    term_sgn = form_sgn;

    do {
        if (dec_comp(den, den_term) == NEUNET_DEC_CMP_EQL) num = dec_add(ans_sgn, num, ans_sgn, num_term, term_sgn);
        else {
            auto sgn = false;
            auto fst = dec_mul(sgn, num, false, den_term, false),
                 snd = dec_mul(sgn, num_term, false, den, false);
            num      = dec_add(ans_sgn, fst, ans_sgn, snd, term_sgn != coef_sgn);
            den      = dec_mul(den, den_term);
        }
        num_term = dec_mul(term_sgn, num_term, term_sgn, num_form, form_sgn);
        den_term = dec_add(den_term, den_form);
        coef_sgn = !coef_sgn;
        // ++iter_cnt;
    } while (dec_series_check(num_term, den_term, prec));
    // std::cout << "next ln2 - " << iter_cnt << std::endl;
    return dec_div(ans_sgn, num, ans_sgn, den, false, prec);
}

net_decimal_data dec_series_ln(bool &ans_sgn, const net_decimal_data &src, uint64_t prec) {
    if (dec_is_zero(src)) return {};
    auto one = dec_init(ans_sgn, 1);
    auto cmp = dec_comp(src, one);
    if (cmp == NEUNET_DEC_CMP_EQL) return {};
    if (cmp == NEUNET_DEC_CMP_LES) return dec_series_ln_2(ans_sgn, src, prec);
    auto ans = src,
         cnt = dec_init(ans_sgn, 0),
         ofs = dec_init(ans_sgn, 0.25);
    while (dec_comp(ans, one) == NEUNET_DEC_CMP_GTR) {
        ans = dec_mul(ans, ofs);
        cnt = dec_add(cnt, one);
    }
    ans = dec_series_ln_2(ans_sgn, ans, prec);
    return dec_add(ans_sgn, ans, ans_sgn, dec_mul(cnt, dec_series_ln4::value(prec)), false);
}


net_decimal_data dec_series_exp(bool &ans_sgn, const net_decimal_data &src_num, const net_decimal_data &src_den, bool src_sgn, uint64_t prec) {
    auto num_form = src_num,
         den_form = src_den;
    dec_frac_red(num_form, den_form);
    
    // auto two = dec_init(ans_sgn, 2),
    //      hlf = dec_init(ans_sgn, 0.5),
    //      ten = dec_init(ans_sgn, 10);
    // auto cnt = 0;
    // while (dec_comp(num_form, ten) == NEUNET_DEC_CMP_GTR) {
    //     num_form = dec_mul(num_form, hlf);
    //     ++cnt;
    // }
    // prec         *= cnt;

    auto num_term = num_form,
         den_term = den_form,
         fra_term = dec_init(ans_sgn, 1),
         frac_cnt = fra_term,
         one_stat = fra_term;
    auto term_sgn = src_sgn;
    if (dec_is_zero(src_den)) {
        den_form = one_stat;
        den_term = one_stat;
    }
    auto num = one_stat,
         den = one_stat,
         tmp = dec_init(ans_sgn, 0);
    do {
        tmp = dec_mul(den_term, fra_term);
        if (dec_comp(den, tmp) == NEUNET_DEC_CMP_EQL) num = dec_add(ans_sgn, num, ans_sgn, num_term, term_sgn);
        else {
            num = dec_add(ans_sgn, dec_mul(num, tmp), ans_sgn, dec_mul(num_term, den), term_sgn);
            den = dec_mul(den, tmp);
        }
        frac_cnt = dec_add(frac_cnt, one_stat);
        fra_term = dec_mul(fra_term, frac_cnt);
        num_term = dec_mul(term_sgn, num_term, term_sgn, num_form, src_sgn);
        den_term = dec_mul(den_term, den_form);
    } while (dec_series_check(num_term, tmp, prec));
    auto ans = dec_div(ans_sgn, num, ans_sgn, den, false, prec);
    
    // for (auto i = 0; i < cnt; ++i) ans = dec_mul(ans, ans);
    // dec_truncate(ans, 32);

    return ans;
}

/* General trigonometry series
 * num_term sine - src_num, cosine - 1
 * den_term sine - src_den, cosine - 0
 * fra_form sine - 1      , cosine - 0
 */
net_decimal_data dec_series_sin_cos(bool &ans_sgn, const  net_decimal_data &src_num, const net_decimal_data &src_den, bool src_sgn, uint64_t prec, net_decimal_data &num_term, net_decimal_data &den_term, net_decimal_data &fra_form) {
    auto sers_sgn = false;
    auto one_form = dec_init(ans_sgn, 1);
    if (dec_is_zero(den_term)) den_term = one_form;
    auto num_form = dec_mul(src_num, src_num),
         den_form = dec_mul(src_den, src_den),
         fra_term = one_form,
         num      = num_term,
         den      = den_term,
         tmp      = dec_init(ans_sgn, 0);
    auto prec_seg = prec / NEUNET_DEC_DIG_MAX,
         prec_rem = prec % NEUNET_DEC_DIG_MAX;
    if (dec_is_zero(den_form)) den_form = one_form;
    ans_sgn = src_sgn;
    do {
        sers_sgn = !sers_sgn;
        num_term = dec_mul(num_term, num_form);
        den_term = dec_mul(den_term, den_form);
        for (auto i = 0; i < 2; ++i) {
            fra_form = dec_add(fra_form, one_form);
            fra_term = dec_mul(fra_term, fra_form);
        }
        tmp = dec_mul(fra_term, den_term);
        num = dec_add(ans_sgn, dec_mul(num, tmp), ans_sgn, dec_mul(num_term, den), sers_sgn != src_sgn);
        den = dec_mul(den, tmp);
    } while (dec_series_check(num_term, tmp, prec));
    return dec_div(ans_sgn, num, ans_sgn, den, false, prec);
}

/* (2k + 1) a = {1, 2}
 * as sin x = 0
 * then x = kπ (k ∈ Z)
 * period = {1, 2} ⊆ Z
 */
bool dec_series_pow_check(const net_decimal_data &times_num, const net_decimal_data &times_den, bool times_sgn, const net_decimal_data &period) {
    auto k_num = times_den,
         k_den = times_num;
    if (dec_is_zero(k_num)) k_num = period;
    else k_num = dec_mul(k_num, period);
    k_num = dec_sub(times_sgn, k_num, times_sgn, k_den, false);
    k_den = dec_mul(k_den, dec_init(times_sgn, 2));
    dec_frac_red(k_num, k_den);
    return dec_is_zero(k_den);
}

// x ^ a, a = times, x = base
net_decimal_data dec_series_pow(bool &ans_sgn, const net_decimal_data &base_num, const net_decimal_data &base_den, bool base_sgn, const net_decimal_data &times_num, const net_decimal_data &times_den, bool times_sgn, uint64_t prec) {
    if (dec_is_zero(base_num)) return {}; // base = 0
    if (base_sgn) {
        // base < 0
        auto prev_sgn = false,
             perd_sgn = false;
        auto prev_ans = dec_series_pow(prev_sgn, base_num, base_den, false, times_num, times_den, times_sgn, prec),
             val_term = dec_init(perd_sgn, 1);
        if (dec_series_pow_check(times_num, times_den, times_sgn, val_term)) goto cplx_pow;
        val_term = dec_init(perd_sgn, 2);
        if (!dec_series_pow_check(times_num, times_den, times_sgn, val_term)) goto null_ret;
        cplx_pow: {
            auto num = dec_init(ans_sgn, 1),
                 den = dec_init(ans_sgn, 0),
                 fra = dec_init(ans_sgn, 0),
                 tmp = dec_mul(val_term, dec_series_pi::value(prec)),
                 ans = dec_mul(ans_sgn, dec_series_sin_cos(perd_sgn, tmp, dec_init(perd_sgn, 0), perd_sgn, prec, num, den, fra), perd_sgn, prev_ans, prev_sgn);
            dec_truncate(ans, prec); // truncate the float point part for specified precision
            return ans;
        } null_ret: return {};
    }
    // base > 0
    auto num_time = times_num,
         den_time = times_den,
         int_time = dec_frac_int_part(num_time, den_time),
         one_term = dec_init(ans_sgn, 1),
         prev_num = dec_init(ans_sgn, 0),
         prev_den = prev_num;
    auto prev_sgn = false,
         rear_sgn = false;
    // base * (int_time + num_time / dem_time) = (base ^ int_time) e ^ (num_time / den_time)ln(base)
    if (!dec_is_zero(int_time)) {
        // base ^ int_time
        prev_num = base_num;
        prev_den = base_den;
        // TODO: optimize -> integral times could be a power of 2
        for (auto i = dec_init(ans_sgn, 0); dec_comp(i, int_time) == NEUNET_DEC_CMP_LES; i = dec_add(i, one_term)) {
            prev_num = dec_mul(prev_num, base_num);
            if (!dec_is_zero(base_den)) prev_den = dec_mul(prev_den, base_den);
        }
        prev_sgn = int_time.it[0] % 2;
    }
    // ln(base_num) - ln(base_den)
    num_time = dec_mul(num_time, dec_add(ans_sgn, dec_series_ln(rear_sgn, base_num, prec), false, dec_series_ln(rear_sgn, base_den, prec), true));
    rear_sgn = rear_sgn != times_sgn;
    auto ans = dec_series_exp(rear_sgn, num_time, den_time, rear_sgn, prec);
    if (!dec_is_zero(int_time)) ans = dec_mul(ans_sgn, dec_div(prev_num, prev_den, prec), prev_sgn != base_sgn, ans, rear_sgn);
    dec_truncate(ans, prec); // truncate the float point part for specified precision
    return ans;
}

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
        if(curr_bit < 0) {
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
                                  omega_unit(std::cos(2 * NEUNET_FFT_PI / len), conj * std::sin(2 * NEUNET_FFT_PI / len));
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
        if(times & 1) ans = ans * base % NEUNET_NTT_EULER_MOD;
        base = base * base % NEUNET_NTT_EULER_MOD;
        times >>= 1;
    }
    return ans;
}

bool dec_ntt(uint64_t *&src, uint64_t len, bool inverse = false) {
    if (num_pad_pow(len, 2) || len < 4 && src == nullptr) return false;
    auto idx_bit_cnt = num_bit_cnt(len - 1);
    auto temp = ptr_init<uint64_t>(len);
    for (auto i = 0ull; i < len; ++i) *(temp + i) = *(src + num_bit_inverse(i, idx_bit_cnt));
    ptr_move(src, std::move(temp));
    for (auto i = 1ull; i < len; i <<= 1) {
        auto g_root = dec_euler_pow(NEUNET_NTT_EULER_MOD_G, (NEUNET_NTT_EULER_MOD - 1) / (i << 1));
        if (inverse) g_root = dec_euler_pow(g_root, NEUNET_NTT_EULER_MOD - 2);
        for (auto j = 0ull; j < len; j += (i << 1)) {
            auto g_unit = 1ull;
            for (auto k = 0ull; k < i; ++k) {
                auto g_unit_front       = *(src + j + k),
                     g_unit_rear        = g_unit * (*(src + i + j + k)) % NEUNET_NTT_EULER_MOD;
                *(src + j + k)          = (g_unit_front + g_unit_rear) % NEUNET_NTT_EULER_MOD;
                *(src + i + j + k)      = (g_unit_front - g_unit_rear + NEUNET_NTT_EULER_MOD) % NEUNET_NTT_EULER_MOD;
                g_unit                  = g_unit * g_root % NEUNET_NTT_EULER_MOD;
            }
        }
    }
    if (inverse) {
        auto inverser = dec_euler_pow(len, NEUNET_NTT_EULER_MOD - 2);
        for (auto i = 0ull; i < len; ++i) *(src + i) = *(src + i) * inverser % NEUNET_NTT_EULER_MOD;
    }
    return true;
}

bool dec_fnt(uint64_t *&coe_c, uint64_t *&coe_a, uint64_t *&coe_b, uint64_t len)
{
    if (coe_c && coe_a && coe_b && len && dec_ntt(coe_a, len) && dec_ntt(coe_b, len)) {
        for(auto i = 0ull; i < len; ++i) *(coe_c + i) = *(coe_a + i) * (*(coe_b + i)) % NEUNET_NTT_EULER_MOD;
        return true;
    } else return false;
}

bool dec_ifnt(uint64_t *&coe_c, uint64_t len) { return (coe_c && len && dec_ntt(coe_c, len, true)); }

NEUNET_END