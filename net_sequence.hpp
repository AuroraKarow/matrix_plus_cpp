NEUNET_BEGIN

template <typename arg> class net_sequence : public net_set<arg> {
public:
    net_sequence(uint64_t init_size = 0, uint64_t alloc_size = 128) : net_set<arg>(alloc_size),
        mem_len(alloc_size) { if (init_size != this->len) this->len = init_size; }
    net_sequence(std::initializer_list<arg> init_list) : net_set<arg>(init_list),
        mem_len(init_list.size()) {}
    net_sequence(net_ptr_base<arg> &&src) : net_set<arg>(std::move(src)) { mem_len = this->len; }
    net_sequence(const net_sequence &src) : net_set<arg>(src),
        mem_len(src.mem_len) { if(mem_len != this->len) ptr_alter(this->ptr, this->len, mem_len); }
    net_sequence(net_sequence &&src) : net_set<arg>(std::move(src)),
        mem_len(src.mem_len) { src.reset(); }
    
    void init(uint64_t init_size, uint64_t alloc_size = 128, bool remain = true) {
        while (init_size > alloc_size) alloc_size <<= 1;
        net_set<arg>::init(alloc_size, remain);
        mem_len   = alloc_size;
        this->len = init_size;
    }

    uint64_t mem_size() const { return mem_len; }

    void ptr_array(net_ptr_base<arg> &&src) {
        net_set<arg>::ptr_array(std::move(src));
        mem_len = this->len;
    }
    net_ptr_base<arg> ptr_array() const { return net_set<arg>::ptr_array(); }

    net_sequence unit(const net_sequence &src) const {
        net_sequence ans(0, 0);
        if (src.len && this->len) {
            ans.init(src.len + this->len);
            for (auto i = 0ull; i < this->len; ++i) *(ans.ptr + i) = *(this->ptr + i);
            for (auto i = 0ull; i < src.len; ++i) *(ans.ptr + i + this->len) = *(src.ptr + i);
        } else if (src.len) ans = src;
        else if (this->len) ans = *this;
        return ans;
    }

    net_sequence unit_union(const net_sequence &src) const {
        net_sequence ans(0, 0);
        if (src.len && this->len) {
            auto cnt  = 0ull;
            auto temp = ptr_com_elem_idx(cnt, this->ptr, this->len, src.ptr, src.len);
            ans.init(this->len + src.len - cnt);
            cnt = 0;
            ptr_copy(ans.ptr, this->ptr, this->len);
            for (auto i = 0ull; i < src.len; ++i) if (!*(temp + i)) *(ans.ptr + this->len + cnt++) = *(src.ptr + i);
            ptr_reset(temp);
        } else if (src.len) ans = src;
        else if (this->len) ans = *this;
        return ans;
    }

    net_sequence unit_intersect(const net_sequence &src) const {
        if (src.len && this->len) {
            auto cnt  = 0ull;
            auto temp = ptr_com_elem_idx(cnt, this->ptr, this->len, src.ptr, src.len);
            net_sequence ans(cnt);
            cnt = 0;
            for (auto i = 0ull; i < src.len; ++i) if (*(temp + i)) *(ans.ptr + cnt++) = *(src.ptr + i);
            ptr_reset(temp);
            return ans;
        } else return net_sequence();
    }

    template<typename ... args> bool insert(uint64_t idx, args &&...paras) {
        if (idx > this->len) return false;
        if (mem_len == 0) init(this->len);
        if (mem_len == this->len) {
            mem_len += mem_len;
            init(this->len, mem_len);
        }
        return ptr_insert(this->ptr, this->len++, std::move(arg(std::forward<args>(paras)...)), idx, false);
    }
    
    template<typename ... args> bool emplace_back(args &&...paras) { return insert(this->len, std::forward<args>(paras)...); }
    
    bool push_back(arg val) { return insert(this->len, std::move(val)); }
    
    arg erase(uint64_t idx, bool mem_shrink = true, long double rate = 0.1) {
        if (idx < this->len) {
            if (mem_shrink && mem_len * rate > this->len) {
                // Release blank
                mem_len >>= 1;
                while (mem_len * rate > this->len) mem_len >>= 1;
                init(this->len, mem_len);
            }
            return ptr_erase(this->ptr, this->len--, idx, false);
        } else return arg();
    }

    bool cut(uint64_t idx, bool successor = true) {
        if (idx >= this->len) return false;
        if ((idx == 0 && successor) || (idx == this->len - 1 && !successor)) {
            clear();
            return true;
        }
        if (successor) this->len = idx;
        else {
            this->len -= idx + 1;
            for (auto i = 0ull; i < this->len; ++i) *(this->ptr + i) = std::move(*(this->ptr + i + idx + 1));
        }
        return true;
    }

    void shrink() { if (this->len != mem_len) init(this->len, this->len); }

    void clear() { this->len = 0; }

    void reset() { net_set<arg>::reset(); mem_len = 0; }

    net_sequence &operator=(const net_sequence &src) {
        net_set<arg>::operator=(src);
        mem_len = src.mem_len;
        if (mem_len != src.len) ptr_alter(this->ptr, this->len, mem_len);
        return *this;
    }
    net_sequence &operator=(net_sequence &&src) {
        net_set<arg>::operator=(std::move(src));
        mem_len = src.mem_len;
        return *this;
    }
    
    ~net_sequence() { mem_len = 0; }

protected:
    uint64_t mem_len = 0;

public:
    __declspec(property(get = ptr_array, put = ptr_array)) net_ptr_base<arg> pointer;
    __declspec(property(get = mem_size))                   uint64_t          memory_length;

    friend std::ostream &operator<<(std::ostream &out, const net_sequence &src) {
        out << "[Length/Memory " << src.len << '/' << src.mem_len << "]\n";
        for (auto i = 0ull; i < src.len; ++i) {
            out << '[' << i << "][\n";
            out << *(src.ptr + i) << "\n]";
            if (i + 1 != src.len) out << '\n';
        }
        return out;
    }
};

NEUNET_END