BAGRT_BEGIN

/* Declaration */

/* Number operation */
uint64_t num_cnt(uint64_t first, uint64_t second, uint64_t dilate = 0);
template<typename arg> arg num_extreme(std::initializer_list<arg> init_num, bool max = true, std::function<bool(arg, arg)> bigger_comp = [](arg _first, arg _second) { return _first > _second; });
uint64_t num_pow_pad_cnt(uint64_t val, uint64_t base, uint64_t min_size = 1, uint64_t pad = 0, uint64_t fold = 0);

/* Net sequence */
template<typename _Ty> class net_sequence
{
protected:
    bool __value_copy(const net_sequence &src);
public:
    net_sequence(uint64_t _size = 0, uint64_t _alloc_size = 128);
    net_sequence(std::initializer_list<_Ty> _init_list);
    net_sequence(net_sequence &src);
    net_sequence(const net_sequence &src);
    net_sequence(net_sequence &&src);
    bool value_copy(net_sequence &src);
    bool value_move(net_sequence &&src);
    void init(uint64_t _size = 1, uint64_t _alloc_size = 128);
    bool ptr_array(_Ty *&&ptr_array, uint64_t array_len);
    bool extend(uint64_t ex_alloc_size, bool extend_size = false);
    uint64_t size();
    uint64_t mem_size();
    template<typename ... args> bool insert(uint64_t idx, args &&...paras);
    template<typename ... args> bool emplace_back(args &&...paras);
    bool push_back(_Ty val);
    _Ty erase(uint64_t idx);
    net_sequence sub_sequence(uint64_t range_first, uint64_t range_second);
    net_sequence sub_sequence(net_sequence<uint64_t> &idx_set);
    bool sort(bool asc = true, std::function<bool(_Ty&, _Ty&)> _func = [](_Ty &_first, _Ty &_second) { return _first > _second; });
    net_sequence unit(net_sequence &src);
    net_sequence unit_union(net_sequence &src);
    net_sequence unit_intersect(net_sequence &src);
    net_sequence<uint64_t> find(_Ty &&target, uint64_t range_first = 0, uint64_t range_second = 0);
    net_sequence<uint64_t> find(_Ty &target, uint64_t range_first = 0, uint64_t range_second = 0);
    _Ty sum(std::function<_Ty(_Ty&, _Ty&)> add_func = [](_Ty &first, _Ty &second) { return first + second; });
    void fill_with(_Ty &&src);
    void memory_set(_Ty &&src);
    bool alter_batch(_Ty &&target, _Ty &&src);
    bool alter_batch(_Ty &target, _Ty &src);
    bool shuffle();
    bool shrink();
    bool reverse();
    bool cut(uint64_t idx, bool sub_seq = true);
    _Ty &operator[](uint64_t idx);
    bool operator==(net_sequence &val);
    bool operator!=(net_sequence &val);
    void operator=(net_sequence &src);
    void operator=(const net_sequence &src);
    void operator=(net_sequence &&src);
    void clear();
    void reset();
    ~net_sequence();
protected:
    _Ty *ptr = nullptr;
    uint64_t len = 0;
    uint64_t mem_len = 0;
public:
    __declspec (property(get=size)) uint64_t length;
    __declspec (property(get=mem_size)) uint64_t memory_length;
    friend std::ostream& operator<<(std::ostream &output, net_sequence &val)
    {
        output << "[Length/Memory][" << val.size() << '/' << val.mem_size() << ']' << std::endl;
        for(auto i=0; i<val.len; ++i)
        {
            output << '[' << i << "][" << std::endl << (int)val.ptr[i] << std::endl << ']';
            if(i+1 < val.len) output << std::endl;
        }
        return output;
    }
};

/* Net set */
template<typename _Ty> class net_set : public net_sequence<_Ty>
{
public:
    net_set(uint64_t _size = 0);
    net_set(std::initializer_list<_Ty> _init_list);
    net_set(net_set &src);
    net_set(const net_set &src);
    net_set(net_set &&src);
    void init(uint64_t _size = 1);
    bool operator==(net_set &val);
    bool operator!=(net_set &val);
    void operator=(net_set &src);
    void operator=(const net_set &src);
    void operator=(net_set &&src);
    ~net_set();
};

/* Net list */
template<typename _Ty> class net_list
{
protected:
struct net_node
{
    _Ty elem;
    net_node *prev = nullptr;
    net_node *next = nullptr;
    ~net_node();
};
struct net_node_elem
{
public:
    net_node_elem(net_node *_curr_node);
    net_node_elem(net_node_elem &src);
    net_node_elem(const net_node_elem &src);
    bool is_null();
    void swap_node(net_node_elem &src);
    const _Ty _node_elem();
    bool _node_elem(_Ty &&src);
    bool _node_elem(_Ty &src);
    net_node_elem _prev_node();
    net_node_elem _next_node();
    void operator=(net_node_elem &src);
    void operator=(const net_node_elem &src);
    bool operator==(net_node_elem &src);
    bool operator!=(net_node_elem &src);
    ~net_node_elem();
protected:
    net_node *curr_node = nullptr;
public:
    __declspec(property(get=_node_elem, put=_node_elem)) _Ty elem;
    __declspec(property(get=_prev_node)) net_node_elem prev;
    __declspec(property(get=_next_node)) net_node_elem next;
};
protected:
    void __value_copy(const net_list &src);
public:
    net_list();
    net_list(std::initializer_list<_Ty> _init_list);
    net_list(net_list &src);
    net_list(const net_list &src);
    net_list(net_list &&src);
    void value_copy(net_list &src);
    void value_move(net_list &&src);
    uint64_t size();
    template<typename ... args> bool insert(uint64_t idx, args &&...paras);
    template<typename ... args> bool emplace_back(args &&...paras);
    bool push_back(_Ty val);
    _Ty erase(uint64_t idx);
    net_list sub_list(uint64_t range_first, uint64_t range_second);
    net_list sub_list(net_sequence<uint64_t> &idx_set);
    net_list unit(net_list &src);
    net_list unit_union(net_list &src);
    net_list unit_intersect(net_list &src);
    net_sequence<uint64_t> find(_Ty &&target, uint64_t range_first = 0, uint64_t range_second = 0);
    net_sequence<uint64_t> find(_Ty &target, uint64_t range_first = 0, uint64_t range_second = 0);
    _Ty sum(std::function<_Ty(_Ty&, _Ty&)> add_func = [](_Ty &first, _Ty &second){ return first + second; });
    bool alter_batch(_Ty &&target, _Ty &&src);
    bool alter_batch(_Ty &target, _Ty &src);
    net_node_elem head_node();
    net_node_elem tail_node();
    net_sequence<_Ty> convert_sequence();
    bool cut(uint64_t idx);
    void reverse();
    _Ty &operator[](uint64_t idx);
    bool operator==(net_list &val);
    bool operator!=(net_list &val);
    void operator=(net_list &src);
    void operator=(const net_list &src);
    void operator=(net_list &&src);
    void reset();
    ~net_list();
protected:
    net_node *head = nullptr, *tail = head;
    uint64_t len = 0;
public:
    __declspec (property(get=size)) uint64_t length;
    __declspec (property(get=head_node)) net_node_elem begin;
    __declspec (property(get=tail_node)) net_node_elem end;
    friend std::ostream& operator<<(std::ostream &output, net_list &val)
    {
        output << "[Length][" << val.len << ']' << std::endl;
        auto p_tool = val.head;
        for(auto i=0; i<val.len; ++i)
        {
            output << '[' << i << "][" << std::endl << p_tool->elem << std::endl << ']';
            if(i+1 < val.len) output << std::endl;
            p_tool = p_tool->next;
        }
        return output;
    }
};

/* Net map */
template<typename _Kty = uint64_t, typename _Ty = double> class net_map
{
protected:
struct net_kv
{
    _Kty key; _Ty value;
    net_kv();
    net_kv(net_kv &src);
    net_kv(const net_kv &src);
    net_kv(net_kv &&src);
    void operator=(net_kv &src);
    void operator=(const net_kv &src);
    void operator=(net_kv &&src);
    bool operator==(net_kv &src);
    bool operator!=(net_kv &src);
    friend std::ostream &operator<<(std::ostream &output, net_kv &src)
    {
        output << "[K][" << std::endl;
        output << src.key << std::endl;
        output << "] -> [V][" << std::endl;
        output << src.value << std::endl << ']';
        return output;
    }
};
    long long detective(long long &threshold, bool &sgn);
    long long next_key(long long hash_key, long long d, long long curr_mem_len);
    uint64_t _find_idx(_Kty &&_key, bool first_lex);
    void _rehash();
    void __value_copy(const net_map &src);
public:
    net_map(std::function<uint64_t(_Kty&&)> hash_key_func = [](_Kty&& src) { return __hash_in_build(std::move(src)); }, uint64_t _alloc_size = 128, uint64_t _rehash_load = 4, double _hash_factor = 0.75);
    net_map(std::initializer_list<_Kty> k_init_list, std::initializer_list<_Ty> v_init_list, std::function<uint64_t(_Kty&&)> hash_key_func = [](_Kty&& src) { return __hash_in_build(std::move(src)); }, uint64_t _alloc_size = 128, uint64_t _rehash_load = 4, double _hash_factor = 0.75);
    net_map(net_map &src);
    net_map(const net_map &src);
    net_map(net_map &&src);
    bool set_hash_func(std::function<uint64_t(_Kty&&)> &&hash_key_func);
    void value_copy(net_map &src);
    void value_move(net_map &&src);
    uint64_t size();
    uint64_t mem_size();
    void rehash(bool enforce = false);
    long long find_idx(_Kty &&key);
    long long find_idx(_Kty &key);
    net_sequence<_Kty> find_key(_Ty &&value);
    net_sequence<_Kty> find_key(_Ty &value);
    bool insert(_Kty &&key, _Ty &&value);
    bool insert(_Kty &key, _Ty &value);
    net_kv erase(_Kty &&key);
    net_kv erase(_Kty &key);
    net_kv index(uint64_t idx);
    bool occupy(uint64_t idx);
    _Ty &operator[](_Kty &&key);
    _Ty &operator[](_Kty &key);
    void operator=(net_map &src);
    void operator=(const net_map &src);
    void operator=(net_map &&src);
    bool operator==(net_map &val);
    bool operator!=(net_map &val);
    void reset();
    ~net_map();
protected:
    double hash_factor = 0.75;
    uint64_t len[2], rehash_load = 4, base_load = rehash_load;
    bool backup = true;
    std::function<uint64_t(_Kty&&)> hash_func;
    net_sequence<net_kv> kv_data[2];
    net_sequence<bool> kv_occupy[2];
public:
    __declspec (property(get=size)) uint64_t length;
    __declspec (property(get=mem_size)) uint64_t memory_length;
    friend std::ostream &operator<<(std::ostream &output, net_map &src)
    {
        output << "[Length/Backup][" << src.len[!src.backup] << '/' << src.kv_data[src.backup].memory_length << ']' << std::endl;
        if(src.len[!src.backup])
        {
            auto out_cnt = 0;
            for(auto i=0; i<src.kv_data[!src.backup].length; ++i) if(src.kv_occupy[!src.backup][i])
            {
                output << src.kv_data[!src.backup][i];
                ++ out_cnt;
                if(out_cnt == src.len[!src.backup]) break;
                else output << std::endl;
            }
        }
        return output;
    }
};

template<typename _Kty = uint64_t> struct clock_timer final
{
private:
struct _dur
{
    long begin  = 0;
    long end    = 0;
    long dur    = 0; 
};
public:
    clock_timer(std::function<uint64_t(_Kty&&)> _hash_func = hash_integer);
    clock_timer(clock_timer &src);
    clock_timer(const clock_timer &src);
    clock_timer(clock_timer &&src);
    bool clock_begin(_Kty &&clock_log_id);
    bool clock_begin(_Kty &clock_log_id);
    bool clock_end(_Kty &&clock_log_id);
    bool clock_end(_Kty &clock_log_id);
    long duration(_Kty &&clock_log_id);
    long duration(_Kty &clock_log_id);
    void operator=(clock_timer &src);
    void operator=(const clock_timer &src);
    void operator=(clock_timer &&src);
    void reset();
    ~clock_timer();
private:
    net_map<_Kty, _dur> clock_log;
};

template<typename _Ty> class net_memory
{
protected:
struct net_ptr
{
public:
    net_ptr(const net_ptr &src);
    net_ptr(net_ptr &src);
    net_ptr(net_ptr &&src);
    net_ptr(net_memory *src, uint64_t _begin_addr, uint64_t _len);
    uint64_t size();
    void operator=(const net_ptr &src);
    void operator=(net_ptr &src);
    void operator=(net_ptr &&src);
    _Ty &operator[](uint64_t addr);
    ~net_ptr();
private:
    net_memory *ptr;
    uint64_t len = 0, begin_addr = 0;
public:
    __declspec(property(get=size)) uint64_t length;
};
struct mem_block
{
    bool occupy = false;
    uint64_t begin_addr = 0, end_addr = 0, len = 0;
    int prev_id = -1, next_id = -1;
};
public:
    net_memory(uint64_t _alloc_size = IDX_MAX);
    net_memory(net_memory &src);
    net_memory(const net_memory &src);
    net_memory(net_memory &&src);
    uint64_t size();
    uint64_t mem_size();
    uint64_t max_mem_size();
    void print_block_info(int id, bool detail = false);
    void print_mem_info(bool detail = false);
    int alloc_mem(uint64_t _alloc_size = 1, _Ty *&&src = nullptr);
    bool free_mem(int &id, bool remain = true);
    void shrink(bool shrink_blk = false);
    bool exist(int id);
    void operator=(net_memory &src);
    void operator=(const net_memory &src);
    void operator=(net_memory &&src);
    net_ptr operator[](int id);
    void reset();
    ~net_memory();
protected:
    net_sequence<_Ty> mem_val;
    net_sequence<mem_block> idx_seq;
    uint64_t head_id = 0, rear_id = 0, len = 0;
public:
    __declspec(property(get=size)) uint64_t length;
    __declspec(property(get=mem_size)) uint64_t memory_length;
    __declspec(property(get=max_mem_size)) uint64_t max_memory_length;
protected:
    void value_assign(net_memory &src);
    uint64_t mem_block_mem_len(mem_block &&src);
};

/* Basic algorithm */
/* Integer sort */
net_set<uint64_t> integer_radix_sort(net_set<uint64_t> &src);
/* Prime */
net_set<uint64_t> integer_primes(uint64_t upper);
net_set<uint64_t> integer_primes_fact(uint64_t val);
/* Fraction */
uint64_t integer_greatest_common_divisor(uint64_t l_val, uint64_t r_val);
uint64_t integer_least_common_multiple(uint64_t l_val, uint64_t r_val);
net_sequence<uint8_t> integer_polynomial_add(bool &ans_sgn, net_sequence<uint8_t> &coe_a, net_sequence<uint8_t> &coe_b, bool sgn_a, bool sgn_b);
net_sequence<uint64_t> integer_polynomial_mult(net_sequence<uint8_t> &coe_a, net_sequence<uint8_t> &coe_b);
// Fast fourier transform
bool integer_fft(std::complex<long double> *&src, uint64_t len, bool inverse = false);
bool integer_dft(std::complex<long double> *&c, std::complex<long double> *&a, std::complex<long double> *&b, uint64_t len);
bool integer_idft(uint64_t *&d, std::complex<long double> *&c, uint64_t len);
// Number theoretical transform
bool integer_ntt(uint64_t *&src, uint64_t len, bool inverse = false);
bool integer_fnt(uint64_t *&coe_c, uint64_t *&coe_a, uint64_t *&coe_b, uint64_t len);
bool integer_ifnt(uint64_t *&coe_c, uint64_t len);

/* Accuracy */
double acc_round(double _val, double acc = 1e-8);
/* Pseudo random number */
double rand_num(double boundary_first, double boundary_second, bool to_sleep = false, double acc = 1e-8);
net_set<uint64_t> rand_idx(uint64_t seq_size, uint64_t amt);
/* Character setting */
std::string charset_input_paragraph(uint64_t buffer_length = 2000);
bool charset_extract_number(std::string &num_str, double *&num_arr, uint64_t &len);
std::string charset_exchange(std::wstring &str_src);
std::wstring charset_exchange(std::string &str_src);

struct decimal final
{
private:
    void __value_copy(const decimal &src);
    decimal euler_itr(decimal &times, decimal& k);
    net_sequence<std::complex<decimal>> euler_eqt(decimal &times);
    net_sequence<uint8_t> dec_coe(bool it_or_ft, bool from_low = true);
    void dec_coe(bool it_or_ft, net_sequence<uint8_t> &src, bool from_low = true);
    static decimal acc_convert(uint64_t con_acc);
public:
    void value_copy(decimal &src);
    void value_move(decimal &&src);
    decimal(decimal &src);
    decimal(const decimal &src);
    decimal(decimal &&src);
    decimal(long double src = 0);
    decimal(std::string src);
    uint64_t it_len();
    uint64_t ft_len();
    bool zero();
    std::string to_string();
    long double to_float();
    int64_t to_integer();
    void _acc(uint64_t con_acc);
    uint64_t _acc();
    decimal abs();
    decimal operator+(decimal &&src);
    decimal operator+(decimal &src);
    decimal operator+(long double src);
    void operator+=(decimal &&src);
    void operator+=(decimal &src);
    void operator+=(long double src);
    decimal &operator++();
    decimal operator++(int);
    decimal operator-(decimal &&src);
    decimal operator-(decimal &src);
    decimal operator-(long double src);
    void operator-=(decimal &&src);
    void operator-=(decimal &src);
    void operator-=(long double src);
    decimal &operator--();
    decimal operator--(int);
    decimal operator*(decimal &&src);
    decimal operator*(decimal &src);
    decimal operator*(long double src);
    void operator*=(decimal &&src);
    void operator*=(decimal &src);
    void operator*=(long double src);
    decimal operator/(decimal &&src);
    decimal operator/(decimal &src);
    decimal operator/(long double src);
    void operator/=(decimal &&src);
    void operator/=(decimal &src);
    void operator/=(long double src);
    decimal operator^(decimal &&src);
    decimal operator^(decimal &src);
    decimal operator^(long double src);
    int64_t operator<<(int64_t bit);
    decimal operator<<(decimal &&bit);
    decimal operator<<(decimal &bit);
    void operator<<=(int64_t bit);
    void operator<<=(decimal &&bit);
    void operator<<=(decimal &bit);
    int64_t operator>>(int64_t bit);
    decimal operator>>(decimal &&bit);
    decimal operator>>(decimal &bit);
    void operator>>=(int64_t bit);
    void operator>>=(decimal &&bit);
    void operator>>=(decimal &bit);
    int64_t operator|(int64_t bit);
    decimal operator|(decimal &&bit);
    decimal operator|(decimal &bit);
    void operator|=(int64_t bit);
    void operator|=(decimal &&bit);
    void operator|=(decimal &bit);
    int64_t operator&(int64_t src);
    decimal operator&(decimal &&src);
    decimal operator&(decimal &src);
    void operator&=(int64_t src);
    void operator&=(decimal &&src);
    void operator&=(decimal &src);
    int64_t operator~();
    int64_t operator%(int64_t src);
    decimal operator%(decimal &&bit);
    decimal operator%(decimal &bit);
    void operator%=(int64_t bit);
    void operator%=(decimal &&bit);
    void operator%=(decimal &bit);
    void operator=(decimal &src);
    void operator=(const decimal &src);
    void operator=(decimal &&src);
    void operator=(long double src);
    void operator=(std::string src);
    bool operator>(decimal &&src);
    bool operator>(decimal &src);
    bool operator>(long double src);
    bool operator>=(decimal &&src);
    bool operator>=(decimal &src);
    bool operator>=(long double src);
    bool operator<(decimal &&src);
    bool operator<(decimal &src);
    bool operator<(long double src);
    bool operator<=(decimal &&src);
    bool operator<=(decimal &src);
    bool operator<=(long double src);
    bool operator==(decimal &&src);
    bool operator==(decimal &src);
    bool operator==(long double src);
    bool operator!=(decimal &&src);
    bool operator!=(decimal &src);
    bool operator!=(long double src);
    decimal reciprocal();
    decimal power(decimal &times);
    void reset();
    ~decimal();
private:
    /* <real>                      <- write ->                   <real>
     * 1|101|001|000|131|080|019|200|001.000|104|000|021|708|004|000|05
     * write -> (i*10 +) <real>   -> (i*10 +)        <real>
     * {1,200,19,80,131,1,101,1}, {0,401,0,120,807,400,0,50}
     * read -> (%10 /=10 low)     -> (%10 /=10 high)
     */
    net_sequence<uint64_t> val[2];
    const bool it = false, ft = true;
    bool sgn = false;
public:
    static uint64_t calculate_digit;
    __declspec(property(get=_acc, put=_acc)) uint64_t accuracy;
    __declspec(property(get=it_len)) uint64_t integer_length;
    __declspec(property(get=ft_len)) uint64_t float_length;
    __declspec(property(get=abs)) decimal absolute;
    friend decimal operator+(long double src, decimal &val) { return val + src; }
    friend void operator+=(long double &src, decimal &val) { src = (val + src).to_float(); }
    friend decimal operator-(long double src, decimal &val) { return decimal(src) - val; }
    friend void operator-=(long double &src, decimal &val) { src = (src - val).to_float(); }
    friend decimal operator*(long double src, decimal &val) { return val * src; }
    friend void operator*=(long double &src, decimal &val) { src = (val * src).to_float(); }
    friend decimal operator/(long double src, decimal &val) { return decimal(src) / val; }
    friend void operator/=(long double &src, decimal &val) { src = (src / val).to_float(); }
    friend int64_t operator<<(int64_t src, decimal &val) { return src << val.to_integer(); }
    friend void operator<<=(int64_t &src, decimal &val) { src <<= val.to_integer(); }
    friend int64_t operator>>(int64_t src, decimal &val) { return src >> val.to_integer(); }
    friend void operator>>=(int64_t &src, decimal &val) { src >>= val.to_integer(); }
    friend int64_t operator|(int64_t src, decimal &val) { return src | val.to_integer(); }
    friend void operator|=(int64_t &src, decimal &val) { src |= val.to_integer(); }
    friend int64_t operator&(int64_t src, decimal &val) { return src & val.to_integer(); }
    friend void operator&=(int64_t &src, decimal &val) { src &= val.to_integer(); }
    friend int64_t operator%(int64_t src, decimal &val) { return src % val.to_integer(); }
    friend void operator%=(int64_t &src, decimal &val) { src %= val.to_integer(); }
    friend bool operator>(long double src, decimal &val) { return val > decimal(src); }
    friend bool operator>=(long double src, decimal &val) { return val >= decimal(src); }
    friend bool operator<(long double src, decimal &val) { return val < decimal(src); }
    friend bool operator<=(long double src, decimal &val) { return val <= decimal(src); }
    friend bool operator==(long double src, decimal &val) { return val == decimal(src); }
    friend bool operator!=(long double src, decimal &val) { return val != decimal(src); }
    friend std::ostream &operator<<(std::ostream &out, decimal &src) { out << src.to_string(); return out; }
    static decimal sin(decimal &src);
    static decimal cos(decimal &src);
    static decimal ln(decimal &src);
    static decimal exp(decimal &src);
};
uint64_t decimal::calculate_digit = 32;

BAGRT_END

#include "bagrt.hpp"