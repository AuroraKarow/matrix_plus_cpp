NEUNET_BEGIN

/* Pointer */

callback_arg arg *ptr_init(uint64_t len) {
    if (len == 0) return nullptr;
    auto ans = new arg [len];
    if constexpr (std::is_copy_assignable_v<arg>) std::fill_n(ans, len, arg{});
    else for (auto i = 0ull; i < len; ++i) *(ans + i) = arg{};
    return ans;
}

callback_arg arg *ptr_init(uint64_t &len, std::initializer_list<arg> init_list) {
         len = init_list.size();
    auto ans = ptr_init<arg>(len);
    if (ans) {
        auto cnt = 0ull;
        for (auto temp : init_list) *(ans + cnt++) = std::move(temp);
    }
    return ans;
}

callback_args void ptr_reset(arg *&first, args *&...others) {
    while (first) {
        delete [] first;
        first = nullptr;
    }
    if constexpr (sizeof...(others) > 0) ptr_reset(others...);
    return;
}

callback_arg bool ptr_elem_equal(const arg *fst_src, uint64_t fst_len, const arg *snd_src, uint64_t snd_len) {
    if(snd_len == fst_len) {
        for (auto i = 0ull; i < fst_len; ++i) if (*(fst_src + i) != *(snd_src + i)) return false;
        return true;
    }
    return false;
}

callback_arg void ptr_print(const arg *src, uint64_t len, bool enter = true) { 
    for(auto i = 0ull; i < len; ++i)  std::cout << *(src + i) << '\n';
    if (enter) std::cout << std::endl;
}

callback_arg void ptr_move(arg *&dest, arg *&&src) {
    ptr_reset(dest);
    dest = std::move(src);
    src  = nullptr;
}

callback_arg bool ptr_shuffle(arg *&src, uint64_t len) {
    if(len && src) {
        std::srand((unsigned)std::time(0));
        for (auto i = len; i > 0; --i) std::swap(*(src + i - 1), *(src + std::rand() % i));
        return true;
    } else return false;
}

callback_arg void ptr_reverse(arg *&src, uint64_t len) { if (len > 1) for (auto i = 0ull, j = len - 1; i < j; ++i, --j) std::swap(*(src + i), *(src + j)); }

callback_arg bool ptr_copy(arg *&dest, const arg *src, uint64_t src_len) {
    if(dest && src && src_len) {
        if constexpr (std::is_copy_assignable_v<arg> && std::is_copy_constructible_v<arg>) std::copy(src, src + src_len, dest);
        else for (auto i = 0ull; i < src_len; ++i) *(dest + i) = *(src + i);
        return true;
    } else return false;
}
callback_arg arg *ptr_copy(const arg *src, uint64_t len) {
    auto ans = ptr_init<arg>(len);
    ptr_copy(ans, src, len);
    return ans;
}

callback_arg void ptr_alter(arg *&src, uint64_t src_len, uint64_t alter_len, bool remain = true) {
    if (src_len == alter_len) {
        if (!remain) std::fill_n(src, src_len, arg());
        return;
    }
    if (alter_len == 0) {
        ptr_reset(src);
        return;
    }
    if (src && src_len && remain) {
        auto temp = ptr_init<arg>(alter_len);
        auto remain_len = src_len > alter_len ? alter_len : src_len;
        for(auto i = 0ull; i < remain_len; ++i) *(temp + i) = std::move(*(src + i));
        ptr_move(src, std::move(temp));
    } else {
        ptr_reset(src);
        src = ptr_init<arg>(alter_len);
    }
}

callback_arg bool ptr_insert(arg *&dest, uint64_t elem_cnt, arg &&src, uint64_t tgt_idx, bool append = true) {
    if (tgt_idx > elem_cnt) return false;
    if (append) ptr_alter(dest, elem_cnt, elem_cnt + 1);
    for(auto i = elem_cnt; i > tgt_idx; --i) *(dest + i) = std::move(*(dest + i - 1));
    *(dest + tgt_idx) = std::move(src);
    return true;
}

callback_arg arg ptr_erase(arg *&src, uint64_t elem_cnt, uint64_t tgt_idx, bool shrink = true) {
    if (tgt_idx < elem_cnt) {
        auto ret = std::move(*(src + tgt_idx));
        for (auto i = tgt_idx; i < elem_cnt; ++i) *(src + i) = std::move(*(src + i + 1));
        if(shrink) ptr_alter(src, elem_cnt, elem_cnt - 1);
        return ret;
    } else return arg();
}

callback_arg bool ptr_sort(arg *&seq_val, uint64_t begin, uint64_t end, bool asc = true) {
    if (end == begin) return true;
    else if (seq_val + begin && seq_val + end) {
        auto pivot = begin,
             slide = end;
        while (slide != pivot)
            if (pivot<slide) {
                if ((asc  && seq_val[slide] < seq_val[pivot]) ||
                    (!asc && seq_val[slide] > seq_val[pivot])) {
                    std::swap(seq_val[slide], seq_val[pivot]);
                    std::swap(slide, pivot);
                    ++slide;
                }
                else --slide;
            } else {
                if ((asc  && seq_val[slide] > seq_val[pivot]) ||
                    (!asc && seq_val[slide] < seq_val[pivot])) {
                    std::swap(seq_val[slide], seq_val[pivot]);
                    std::swap(slide, pivot);
                    --slide;
                }
                else ++slide;
            }
        auto begin_flag = true, 
             end_flag   = true;
        if (pivot != begin) begin_flag = ptr_sort(seq_val, begin, pivot - 1, asc);
        if (pivot != end) end_flag = ptr_sort(seq_val, pivot + 1, end, asc);
        return (begin_flag && end_flag);
    } else return false;
}

callback_arg bool ptr_dup_remove(uint64_t &ans_len, arg *&src, uint64_t len) {
    if (len == 0 || src == nullptr) return false;
    if (len == 1) return true;
    ans_len = 0;
    ptr_sort(src, 0, len - 1);
    for (auto i = 1ull; i < len; ++i) if (*(src + i) != *(src + ans_len)) {
        ++ans_len;
        *(src + ans_len) = *(src + i);
    }
    ++ ans_len;
    if (ans_len != len) ptr_alter(src, len, ans_len);
    return true;
}

callback_arg arg *ptr_sub(uint64_t &sub_len, const arg *src, uint64_t src_len, uint64_t fst_rng, uint64_t snd_rng) {
    if (fst_rng > snd_rng) std::swap(fst_rng, snd_rng);
    if (snd_rng >= src_len) {
        sub_len = 0;
        return nullptr;
    }
    auto rng = snd_rng - fst_rng + 1;
    if(rng == src_len) {
        sub_len = src_len;
        return ptr_copy(src, src_len);
    } else {
        sub_len = rng;
        auto ptr_temp = src + fst_rng;
        return ptr_copy(ptr_temp, rng);
    }
}
callback_arg arg *ptr_sub(uint64_t &sub_len, const arg *src, uint64_t src_len, uint64_t *idx_arr, uint64_t arr_len) {
    if (src == nullptr || src_len == 0 || !ptr_dup_remove(arr_len, idx_arr, arr_len)) return nullptr;
    auto ans     = ptr_init<arg>(arr_len);
         sub_len = arr_len;
    for (auto i = 0ull; i < arr_len; ++i)
        if (*(idx_arr + i) < src_len) *(ans + i) = *(src + (*(idx_arr + i)));
        else {
            sub_len = 0;
            ptr_reset(ans);
            break;
        }
    return ans;
}

callback_arg arg *ptr_concat(const arg *src_fst, uint64_t len_fst, const arg *src_snd, uint64_t len_snd) {
    if (src_fst && len_fst && src_snd && len_snd) {
        auto ans = ptr_init<arg>(len_fst + len_snd);
        ptr_copy(ans, src_fst, len_fst);
        auto ans_snd = ans + len_fst;
        ptr_copy(ans_snd, src_snd, len_snd);
        ans_snd = nullptr;
        return ans;
    } else return nullptr;
}

callback_arg arg *ptr_union(uint64_t &ans_len, const arg *src_fst, uint64_t len_fst, const arg *src_snd, uint64_t len_snd) {
    auto temp_snd     = ptr_copy(src_snd, len_snd);
    auto temp_snd_len = len_snd;
    for (auto i = 0ull; i < len_fst; ++i) for (auto j = 0ull; j < temp_snd_len; ++j) if (*(temp_snd + j) == *(src_fst + i)) {
        ptr_erase(temp_snd, temp_snd_len, j, false);
        --temp_snd_len;
        break;
    }
    auto ans = ptr_concat(src_fst, len_fst, temp_snd, temp_snd_len);
    ptr_reset(temp_snd);
    ans_len = len_fst + temp_snd_len;
    return ans;
}

callback_arg arg *ptr_intersect(uint64_t &ans_len, const arg *src_fst, uint64_t len_fst, const arg *src_snd, uint64_t len_snd) {
    auto temp_snd     = ptr_copy(src_snd, len_snd);
    auto temp_snd_len = len_snd,
         ans_len_temp = len_fst > len_snd ? len_fst : len_snd;
    auto ans          = ptr_init<arg>(ans_len_temp);
         ans_len      = 0;
    for (auto i = 0ull; i < len_fst; ++i) for (auto j = 0ull; j < temp_snd_len; ++j) if (*(temp_snd + j) == *(src_fst + i)) {
        ptr_insert(ans, ans_len, ptr_erase(temp_snd, temp_snd_len, j, false), ans_len++, false);
        --temp_snd_len;
        break;
    }
    ptr_reset(temp_snd);
    if (ans_len != ans_len_temp) ptr_alter(ans, ans_len_temp, ans_len);
    return ans;
}

/** [Arrays common elements index]
 * @brief Find common element in two arrays
 * @param comm_cnt  [Out]   Common elements count
 * @param axis      [In]    Axis array, at rest array for comparing.
 * @param axis_len  [In]    Axis array length.
 * @param src       [In]    Source array, find the index of common element comparing with the axis array.
 * @param src_len   [In]    Source array length
 * @return Boolean array. Length of this array equals to the source array. Each element value of this array which is the index of common element in axis array would be true, otherwise false.
 */
callback_arg bool *ptr_com_elem_idx(uint64_t &comm_cnt, const arg *axis, uint64_t axis_len, const arg *src, uint64_t src_len) {
    auto ans = ptr_init<bool>(src_len);
    comm_cnt = 0;
    for (auto i = 0ull; i < axis_len; ++i) if (src_len != comm_cnt) for (auto j = 0ull; j < src_len; ++j) if (*(src + j) == *(axis + i) && !(*(ans + j))) {
        *(ans + j) = true;
        ++comm_cnt;
        break;
    }
    return ans;
}

callback_arg bool ptr_cut(arg *&src, uint64_t src_len, uint64_t tgt_idx, bool successor = true) {
    if(src_len && tgt_idx < src_len && src) {
        if(successor) ptr_alter(src, src_len, tgt_idx);
        else {
            auto ans_len  = src_len - tgt_idx - 1;
            auto ans_addr = src + tgt_idx + 1,
                 ans      = ptr_copy(ans_addr, ans_len);
            ptr_move(src, std::move(ans));
        }
        return true;
    } else return false;
}

callback_arg ul_ptr ptr_find(uint64_t &ans_len, const arg *src, uint64_t src_len, const arg &tgt, uint64_t fst_rng, uint64_t snd_rng) {
    if (src && src_len && fst_rng < src_len && snd_rng < src_len) {
        auto ans     = ptr_init<uint64_t>(src_len);
             ans_len = 0;
        if (fst_rng > snd_rng) std::swap(fst_rng, snd_rng);
        for (auto i = fst_rng; i <= snd_rng; ++i) if(*(src + i) == tgt) *(ans + ans_len++) = i;
        if (ans_len != src_len) ptr_alter(ans, src_len, ans_len);
        return ans;
    } else return nullptr;
}
callback_arg ul_ptr ptr_find(uint64_t &ans_len, const arg *src, uint64_t src_len, const arg &tgt) { return ptr_find(ans_len, src, src_len, tgt, 0, src_len - 1); }
callback_arg ul_ptr ptr_find(uint64_t &ans_len, const arg *src, uint64_t src_len, const arg &tgt, uint64_t *idx_arr, uint64_t arr_len) {
    auto sub_temp = ptr_sub(arr_len, src, src_len, idx_arr, arr_len);
    auto ans      = ptr_find(ans_len, sub_temp, arr_len, tgt);
    for(auto i = 0ull; i < ans_len; ++i) *(ans + i) = *(idx_arr + (*(ans + i)));
    return ans;
}

double *ptr_narr_float(const long double *src, uint64_t len) {
    if (!(len && src)) return nullptr;
    auto ans = ptr_init<double>(len);
    if constexpr (sizeof(double) == sizeof(long double)) std::copy(src, src + len, ans);
    else for (auto i = 0ull; i < len; ++i) *(ans + i) = (*(src + i));
    return ans;
}

/* Number */

bool num_booleam(bool fst, bool snd, uint8_t boolean_type) {
    switch (boolean_type)
    {
    case 0: return fst || snd;
    case 1: return fst && snd;
    case 2: return fst != snd;
    default: return false;
    }
}

ul_ptr num_primes(uint64_t &len, uint64_t upper) {
    auto ans = ptr_init<uint64_t>(upper + 1);
    auto pmf = ptr_init<bool>(upper + 1);
    len = 0;
    for (auto i = 2ull; i <= upper; ++i) {
        if (!(*(pmf + i))) *(ans + len++) = i;
        for (auto j = 0ull; j < len; ++j) {
            if (i * (*(ans + j)) > upper) break;
            *(pmf + (i * (*(ans + j)))) = true;
            if (!(i % (*(ans + j)))) break;
        }
    }
    ptr_alter(ans, upper, len);
    ptr_reset(pmf);
    return ans;
}

ul_ptr num_primes_factor(uint64_t &len, uint64_t val) {
    auto prime_len  = 0ull;
    auto prime_temp = num_primes(prime_len, val),
         ans        = ptr_init<uint64_t>(prime_len);
         len        = 0;
    for (auto i = 0ull; i < prime_len; ++i)  while (val % (*(prime_temp + i)) == 0) {
        *(ans + len++) = *(prime_temp + i);
        val           /= (*(prime_temp + i));
    }
    if(len != prime_len) ptr_alter(ans, prime_len, len);
    return ans;
}

uint64_t num_gcd_lcm(uint64_t v_fst, uint64_t v_snd, bool gcd = true) {
    if (v_fst > v_snd) std::swap(v_fst, v_snd);
    if ((v_fst == v_snd) || (v_snd % v_fst == 0)) return v_fst;
    auto pm_fact_fst_len = 0ull,
         pm_fact_snd_len = 0ull,
         unit_len        = 0ull,
         ans             = 1ull;
    auto pm_fact_fst     = num_primes_factor(pm_fact_fst_len, v_fst),
         pm_fact_snd     = num_primes_factor(pm_fact_snd_len, v_snd),
         unit_ptr        = (uint64_t *)nullptr;
    if (gcd) unit_ptr = ptr_intersect(unit_len, pm_fact_fst, pm_fact_fst_len, pm_fact_snd, pm_fact_snd_len);
    else unit_ptr = ptr_union(unit_len, pm_fact_fst, pm_fact_fst_len, pm_fact_snd, pm_fact_snd_len);
    for (auto i = 0ull; i < unit_len; ++i) ans *= (*(unit_ptr + i));
    ptr_reset(unit_ptr, pm_fact_fst, pm_fact_snd);
    return ans;
}

uint64_t num_cnt(uint64_t first, uint64_t second, uint64_t dilate = 0) {
    if (first > second) std::swap(first, second);
    auto ans = second - first;
    return (ans / (dilate + 1) + 1);
}

long double num_rate(long double numerator, long double denominator) {
    assert(denominator);
    return numerator / denominator;
}

callback_arg arg num_extreme(std::initializer_list<arg> init_num, bool max = true) {
    auto ext_val = *init_num.begin();
    for (auto temp : init_num)
    {
        if (max && temp > ext_val) ext_val = temp;
        if (!max && temp < ext_val) ext_val = temp;
    }
    return ext_val;
}

uint64_t num_pad_pow(uint64_t val, uint64_t base, uint64_t min_size = 1, uint64_t pad = 0, uint64_t fold = 0) {
    if (val <= min_size) return pad;
    else if (val % base)
    {
        auto pad_curr = 0ull;
        while (val % base) {
            ++val;
            ++pad_curr;
        }
        pad_curr *= (fold + 1);
        return num_pad_pow(val / base, base, min_size, pad + pad_curr, (fold + 1) * (base - 1) + fold);
    } else return num_pad_pow(val / base, base, min_size, pad, (fold + 1) * (base - 1) + fold);
}

uint64_t num_unsign(uint64_t val) { return val * (-1) < val ? 0 : val; }

uint32_t num_swap_endian(uint32_t val) {
	val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF);
	return (val << 16) | (val >> 16);
}

long long num_bit_reverse(long long src, uint8_t bit_cnt = 3) {
    long long ans = 0;
    while (bit_cnt)
    {
        ans <<= 1;
        ans  += src & 1;
        src >>= 1;
        --bit_cnt;
    }
    return ans;    
}

uint64_t num_bit_cnt(long long src) {
    auto ans = 0ull;
    while(src)
    {
        ++ ans;
        src >>= 1;
    }
    return ans;    
}

long double num_rand(long double fst_rng = 0, long double snd_rng = 0, uint64_t acc = 8) {
    if (fst_rng == snd_rng) return (((long double)lib_rand_e() / (long double)lib_rand_e._Max) - .5l) * 2.0l;
    else {
        // random seed
        auto curr_time = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now()).time_since_epoch().count();
        curr_time /= 100;
        // interval
        if (fst_rng > snd_rng) std::swap(fst_rng, snd_rng);
        long double ans   = 0;
        auto        times = 1ull;
        for (auto i = 0ull; i < acc || ans < snd_rng; ++i) {
            ans       *= 10;
            ans       += curr_time % 10;
            curr_time /= 10;
            times     *= 10;
        }
        // rectify
        auto rng = snd_rng - fst_rng;
        ans /= times / rng;
        ans += fst_rng;
        return ans;
    }
}

callback_arg arg *num_rand(uint64_t amt, arg fst_rng = 0, arg snd_rng = 0, bool order = true, uint64_t acc = 8) {
    if (amt == 0) return nullptr;
    auto ans = ptr_init<arg>(amt);
    auto cnt = 0;
    while (cnt < amt) {
        arg temp = num_rand(fst_rng, snd_rng, acc);
        for (auto i = 0ull; i < cnt; ++i) if (temp == *(ans + i)) continue;
        *(ans + cnt++) = temp;
    }
    if (order) ptr_sort(ans, 0, amt - 1);
    return ans;
}

/* character array */

void str_arr_reset(ch_str *&src, uint64_t arr_len) {
    for(auto i = 0ull; i < arr_len; ++i) ptr_reset(*(src + i));
    ptr_reset(src);
}

ch_str str_stream_in(uint64_t buffer_len = 1e5) {
    std::fflush(stdin);
    std::cout << "Please end submit with double enter." << std::endl;
    std::cout << std::endl;
    auto cmtr  = ptr_init<char>(buffer_len);
    auto i     = 0ull;
    char temp  = 0, 
         enter = 0; 
    while (i < buffer_len) {
        temp = std::getchar();
        if (temp == '\n') {
            enter = std::getchar();
            if(enter == '\n') break;
            else {
                cmtr[i]   = '\n';
                cmtr[++i] = enter;
            }
        } else {
            if (i == buffer_len) {
                ptr_alter(cmtr, buffer_len, buffer_len + buffer_len);
                buffer_len += buffer_len;
            }
            cmtr[i] = temp;
        }
        ++i;
    }
    cmtr[i] = '\0';
    if (i != buffer_len) ptr_alter(cmtr, buffer_len, i);
    return cmtr;
}

ch_str *str_split(uint64_t &len, const ch_str src, const char syb) {
    auto cnt     = 0ull,
         src_len = std::strlen(src);
    auto ans     = ptr_init<ch_str>(src_len);
         len     = 0;
    while (*(src + cnt) != '\0') {
        auto cnt_temp = 0ull;
        *(ans + len) = ptr_init<char>(src_len);
        while (*(src + cnt) != syb && *(src + cnt) != '\0') *(*(ans + len) + cnt_temp++) = *(src + cnt++);
        *(*(ans + len) + cnt_temp) = '\0';
        if (cnt_temp != src_len) ptr_alter(*(ans + len), src_len, ++cnt_temp);
        ++len;
        if (*(src + cnt) != '\0') ++cnt;
    }
    return ans;
}

ch_str str_charset_exchange(const wch_str src)
{
    auto len     = std::wcslen(src),
         buf_len = 0ull;
    auto ans     = ptr_init<char>(len + 1);
    wcstombs_s(&buf_len, ans, len + 1, src, len);
    *(ans + len) = '\0';
    return ans;
}
wch_str str_charset_exchange(const ch_str src)
{
    // setlocale(LC_ALL, "zh_CN.UTF-8");
    auto len     = std::strlen(src),
         buf_len = 0ull;
    auto ans     = ptr_init<wchar_t>(len + 1);
    mbstowcs_s(&buf_len, ans, len + 1, src, len);
    *(ans + len) = L'\0';
    return ans;
}

ch_str str_cat(const ch_str fst, const ch_str snd) {
    auto fst_len = std::strlen(fst),
         snd_len = std::strlen(snd),
         ans_len = fst_len + snd_len;
    auto ans     = ptr_init<char>(ans_len + 1);
    auto p_tool = ans;
    for (auto i = 0ull; i < fst_len; ++i) *(ans + i) = *(fst + i);
    for (auto i = 0ull; i < snd_len; ++i) *(ans + i + fst_len) = *(snd + i);
    ans[ans_len] = '\0';
    return ans;    
}
template <typename ... arg> ch_str str_cat(const ch_str fst, const ch_str snd, arg &&... src) {
    auto curr_ans = str_cat(fst, snd);
    auto ans      = str_cat(curr_ans, src ...);
    ptr_reset(curr_ans);
    return ans;
}

/* pointer */
template <typename arg> struct net_ptr_base {
    arg      *ptr_base = nullptr;
    uint64_t len       = 0;
    static net_ptr_base init(arg *&&src, uint64_t src_len) {
        net_ptr_base ans;
        ptr_move(ans.ptr_base, std::move(src));
        ans.len = src_len;
        return ans;
    }
    void init(uint64_t alloc_size) {
        ptr_alter(ptr_base, len, alloc_size, false);
        len = alloc_size;
    }
    void reset() {
        len = 0;
        ptr_reset(ptr_base);
    }
};

/* iterator */
template <typename arg, typename inst_t> struct net_iterator_base {
public:
    net_iterator_base(const inst_t *ptr_src = nullptr) : ptr(ptr_src) {}

    virtual bool operator==(const net_iterator_base &val) const { return ptr == val.ptr;}

    virtual arg operator*() const = 0;

    virtual net_iterator_base &operator++() { return *this; }

    virtual net_iterator_base &operator--() { return *this; }

    virtual ~net_iterator_base() { ptr = nullptr; }

protected: const inst_t *ptr = nullptr;
};

/* Hash key */

template <typename arg> uint64_t hash_in_built(const arg &src) {
    if constexpr (std::is_same_v<arg, std::string>) return std::hash<std::string>{}(src);
    else if constexpr (std::is_integral_v<arg>) return src;
    else if constexpr (std::is_floating_point_v<arg>) return std::hash<arg>{}(src);
    else return 0;
}

long long hash_detect(long long &threshold, bool &sgn) {
    if (threshold) if (sgn) { 
        auto temp = (-1) * threshold * threshold;
              sgn = false;
        ++ threshold;
        return temp;
    } else {
        sgn = true;
        return threshold * threshold;
    } else {
        auto temp = threshold;
             sgn  = false;
        ++ threshold;
        return temp;
    }
}
long long hash_next_key(long long hash_key, long long detect_v, long long curr_mem_len) { return (hash_key + detect_v) % curr_mem_len; }

NEUNET_END