NEUNET_BEGIN

/* net_set */
template <typename arg> class net_set {
public:
    /* net_set iterator */
    struct iterator final : net_iterator_base<arg, net_set<arg>> {
    public:
        iterator(const net_set *src = nullptr, uint64_t idx = 0) : net_iterator_base<arg, net_set<arg>>(src),
            curr_idx(idx) {}

        virtual bool operator==(const iterator &val) const { return curr_idx == val.curr_idx && net_iterator_base<arg, net_set<arg>>::operator==(val); }

        virtual bool operator!=(const iterator &val) const { return !(*this == val); }
    
        virtual arg operator*() const {
            if (this->ptr) return (*(this->ptr))[curr_idx];
            else return arg();
        }
    
        virtual iterator &operator++() {
            if (this->ptr) {
                ++curr_idx;
                if (curr_idx == this->ptr->len) {
                    this->ptr = nullptr;
                    curr_idx  = 0;
                }
            }
            return *this;
        }
    
        virtual iterator operator++(int) { auto temp = *this; ++*this; return temp; }
    
        virtual iterator &operator--() {
            if (curr_idx == 0) this->ptr = nullptr;
            else --curr_idx;
            return *this;
        }
    
        virtual iterator operator--(int) { auto temp = *this; --*this; return temp; }
    
        ~iterator() { curr_idx = 0; }
    private: uint64_t curr_idx = 0;
    };

protected:
    void value_copy(const net_set &src) {
        ptr_alter(ptr, len, src.len, false);
        ptr_copy(ptr, src.ptr, src.len);
        len = src.len;
    }
    
    void value_move(net_set &&src) {
        len = src.len;
        ptr_move(ptr, std::move(src.ptr));
        src.reset();
    }

public:
    net_set(uint64_t alloc_size = 0) : ptr(ptr_init<arg>(alloc_size)), len(alloc_size) {}
    net_set(std::initializer_list<arg> init_list) :
        ptr(ptr_init(len, init_list)),
        len(init_list.size()) {}
    net_set(net_ptr_base<arg> &&src) :
        len(src.len) { ptr_move(ptr, std::move(src.ptr_base)); src.len = 0; }
    net_set(const net_set &src) { value_copy(src); }
    net_set(net_set &&src) { value_move(std::move(src)); }
    
    void init(uint64_t alloc_size, bool remain = true) {
        if (alloc_size == 0) {
            reset();
            return;
        }
        auto pre_len = len;
        len = alloc_size;
        ptr_alter(ptr, pre_len, len, remain);
    }
    
    void ptr_array(net_ptr_base<arg> &&src) {
        *this        = net_set(std::move(src));
        src.ptr_base = nullptr;
        src.len      = 0;
    }
    net_ptr_base<arg> ptr_array() const {
        net_ptr_base<arg> ans;
        ans.ptr_base = ptr_copy(ptr, len);
        ans.len      = len;
        return ans;
    }
    
    uint64_t size() const { return len; }
    
    template<typename ... args> bool insert(uint64_t idx, args &&...paras) {
        if (idx > len) return false;
        return ptr_insert(ptr, len++, std::move(arg(std::forward<args>(paras)...)), idx);
    }
    
    template<typename ... args> bool emplace_back(args &&...paras) { return insert(len, std::forward<args>(paras)...); }
    
    bool push_back(arg val) { return insert(len, std::move(val)); }
    
    arg erase(uint64_t idx) { if (idx < len) return ptr_erase(ptr, len--, idx); else return arg(); }

    net_set sub_set(uint64_t fst_rng, uint64_t snd_rng) const {
        net_set ans;
        ans.ptr = ptr_sub(ans.len, ptr, len, fst_rng, snd_rng);
        return ans;
    }
    net_set sub_set(const net_set<uint64_t> &idx_set) const {
        net_set ans;
        auto idx_ptr     = idx_set.ptr_array();
             ans.ptr     = ptr_sub(ans.len, ptr, len, idx_ptr.ptr_base, idx_ptr.len);
             idx_ptr.len = 0;
        ptr_reset(idx_ptr.ptr_base);
        return ans;
    }

    bool sort(bool asc = true) { return ptr_sort(ptr, 0, len - 1, asc); }

    net_set unit(const net_set &src) const {
        if (len && src.len) return net_set(net_ptr_base<arg>::init(ptr_concat(ptr, len, src.ptr, src.len), len + src.len));
        else if (src.len) return net_set(src);
        else return net_set(*this);
    }

    net_set unit_union(const net_set &src) const {
        if (len && src.len) {
            net_set ans;
            ans.ptr = ptr_union(ans.len, ptr, len, src.ptr, src.len);
            return ans;
        } else if (src.len) return net_set(src);
        else return net_set(*this);
    }

    net_set unit_intersect(const net_set &src) const {
        if (!(len && src.len)) return net_set();
        net_set ans;
        ans.ptr = ptr_intersect(ans.len, ptr, len, src.ptr, src.len);
        return ans;
    }

    net_set<uint64_t> find(const arg &tgt, uint64_t fst_rng, uint64_t snd_rng) const {
        if (fst_rng == snd_rng) {
            if (*(ptr + fst_rng) == tgt) return {fst_rng};
            else return net_set<uint64_t>();
        }
        net_set<uint64_t> ans;
        net_ptr_base<uint64_t> ptr_temp;
        ptr_temp.ptr_base = ptr_find(ptr_temp.len, ptr, len, tgt, fst_rng, snd_rng);
        ans.ptr_array(std::move(ptr_temp));
        return ans;
    }
    net_set<uint64_t> find(const arg &tgt) const { return find(tgt, 0, len - 1); }
    net_set<uint64_t> find(const arg &tgt, const net_set<uint64_t> &idx_set) const {
        net_set<uint64_t> ans;
        net_ptr_base<uint64_t> ptr_temp;
        auto idx_ptr           = idx_set.ptr_array();
             ptr_temp.ptr_base = ptr_find(ptr_temp.len, ptr, len, tgt, idx_ptr.ptr_base, idx_ptr.len);
        idx_ptr.reset();
        ans.ptr_array(std::move(ptr_temp));
        return ans;
    }

    arg sigma() const {
        arg ans {};
        if (len) ans = *ptr;
        else return ans;
        for (auto i = 1ull; i < len; ++i) ans += (*(ptr + i));
        return ans;
    }
    arg sigma(uint64_t fst_rng, uint64_t snd_rng) const { return sub_set(fst_rng, snd_rng).sigma(); }
    arg sigma(const net_set<uint64_t> &idx_set) const { return sub_set(idx_set).sigma(); }

    arg pi() const {
        arg ans {};
        if (len) ans = *ptr;
        else return ans;
        for (auto i = 1ull; i < len; ++i) ans *= (*(ptr + i));
        return ans;
    }
    arg pi(uint64_t fst_rng, uint64_t snd_rng) const { return sub_set(fst_rng, snd_rng).pi(); }
    arg pi(const net_set<uint64_t> &idx_set) const { return sub_set(idx_set).pi(); }

    void fill_with(const arg &src) { if (len) std::fill_n(ptr, len, src); }

    void supersede(const arg &tgt, const arg &src, const net_set<uint64_t> &idx_set) {
        auto tgt_idx = find(tgt, idx_set);
        for (auto i = 0ull; i < tgt_idx.size(); ++i) *(ptr + tgt_idx[i]) = src;
    }
    void supersede(const arg &tgt, const arg &src, uint64_t fst_rng, uint64_t snd_rng) { return supersede(tgt, src, find(tgt, fst_rng, snd_rng)); }
    void supersede(const arg &tgt, const arg &src) { return supersede(tgt, src, 0, len - 1); }

    bool shuffle() { return ptr_shuffle(ptr, len); }

    void reverse() { ptr_reverse(ptr, len); }

    bool cut(uint64_t idx, bool successor = true) {
        auto prev_len = len;
        if (successor) len = idx;
        else len -= idx + 1;
        return ptr_cut(ptr, prev_len, idx, successor);
    }
    
    iterator begin() const {
        if (len) return iterator(this, 0);
        else return end();
    }
    
    iterator end() const { return iterator(nullptr, 0); }
    
    void reset() { len = 0; ptr_reset(ptr); }
    
    arg &operator[](uint64_t idx) const {
        assert(idx < len);
        return *(ptr + idx);
    }

    bool operator==(const net_set &src) const { return ptr_elem_equal(ptr, len, src.ptr, src.len); }

    bool operator!=(const net_set &src) const { return !(*this == src); }
    
    net_set &operator=(const net_set &src) { value_copy(src); return *this; }
    net_set &operator=(net_set &&src) { value_move(std::move(src)); return *this; }

    template <typename seq> explicit operator seq () const {
        seq ans;
        ans.init(len);
        for (auto i = 0ull; i < len; ++i) ans[i] = *(ptr + i);
        return ans;
    }
    
    virtual ~net_set() { reset(); }

protected:
    arg *ptr     = nullptr;
    uint64_t len = 0;

public:
    __declspec(property(get = ptr_array, put = ptr_array)) net_ptr_base<arg> pointer;
    __declspec(property(get = size))                       uint64_t          length;
    __declspec(property(get = sigma))                      arg               sum;
    __declspec(property(get = pi))                         arg               product;
    
    friend std::ostream &operator<<(std::ostream &out, const net_set &src) {
        out << "[Length " << src.len << "]\n";
        for (auto i = 0ull; i < src.len; ++i) {
            out << '[' << i << "][\n";
            out << *(src.ptr + i) << "\n]";
            if (i + 1 != src.len) out << '\n';
        }
        return out;
    }
};

NEUNET_END