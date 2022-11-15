NEUNET_BEGIN

/* node */
template <typename arg> struct net_node {
    arg      elem  {};
    net_node *prev = nullptr,
             *next = nullptr;
    ~net_node() {
        prev = nullptr;
        while(next)
        {
            delete next;
            next = nullptr;
        }
    }
};

template <typename arg> class net_list {
public:
    struct iterator final : net_iterator_base<arg, net_node<arg>> {
    public:
        iterator(const net_node<arg> *node = nullptr) : net_iterator_base<arg, net_node<arg>>(node) {}

        virtual bool operator==(const iterator &val) const { return net_iterator_base<arg, net_node<arg>>::operator==(val); }

        virtual bool operator!=(const iterator &val) const { return !(*this == val); }
    
        virtual arg operator*() const {
            if (this->ptr) return this->ptr->elem;
            else return arg();
        }
    
        virtual iterator &operator++() {
            if (this->ptr) this->ptr = this->ptr->next;
            return *this;
        }
    
        virtual iterator operator++(int) { auto temp = *this; ++*this; return temp; }
    
        virtual iterator &operator--() {
            if (this->ptr) this->ptr = this->ptr->prev;
            return *this;
        }
    
        virtual iterator operator--(int) { auto temp = *this; --*this; return temp; }
    };

protected:
    void value_copy(const net_list &src) {
        if (src.head == nullptr) { reset(); return; }
        len = src.len;
        if(!head) {
            head       = new net_node<arg>;
            head->next = nullptr;
            head->prev = nullptr;
            tail       = head;
        }
        head->elem  = src.head->elem;
        auto p_tool = head, src_tool = src.head;
        while (src_tool->next) {
            if (!p_tool->next) {
                p_tool->next       = new net_node<arg>;
                p_tool->next->prev = p_tool;
            }
            p_tool->next->elem = src_tool->next->elem;
            p_tool             = p_tool->next;
            src_tool           = src_tool->next;
        }
        tail       = p_tool;
        p_tool     = p_tool->next;
        tail->next = nullptr;
        delete p_tool; p_tool = nullptr;
    }

    void value_move(net_list &&src) {
        if(len) delete head;
        head     = std::move(src.head);
        src.head = nullptr;
        tail     = src.tail;
        src.tail = nullptr;
        len      = src.len;
        src.len  = 0;
    }

public:
    net_list() {}
    net_list(std::initializer_list<arg> init_list) { for (auto temp : init_list) emplace_back(std::move(temp)); }
    net_list(const net_list &src) { value_copy(src); }
    net_list(net_list &&src) { value_move(std::move(src)); }

    uint64_t size() const { return len; }

    template<typename ... args> bool insert(uint64_t idx, args &&...paras) {
        if (idx > len) return false;
        auto i_node  = new net_node<arg>;
        i_node->elem = arg(std::forward<args>(paras)...);
        i_node->next = nullptr;
        i_node->prev = nullptr;
        if (len) if (idx) if (idx == len) {
            tail->next   = i_node;
            i_node->prev = tail;
            tail         = i_node;
        } else {
            auto p_tool = head;
            for(auto i = 0ull; i < idx; ++i) p_tool = p_tool->next;
            i_node->next       = p_tool;
            p_tool->prev->next = i_node;
            i_node->prev       = p_tool->prev;
            p_tool->prev       = i_node;
        } else {
            i_node->next = head;
            head->prev   = i_node;
            head         = i_node;
        } else {
            head = i_node;
            tail = head;
        }
        ++len;
        return true;
    }

    template<typename ... args> bool emplace_back(args &&...paras) { return insert(len, std::forward<args>(paras)...); }

    bool push_back(arg val) { return insert(len, std::move(val)); }

    arg erase(uint64_t idx) {
        arg temp{};
        if (idx >= len) return temp;
        net_node<arg> *p_tool = nullptr;
        if (idx) if (idx + 1 == len) {
            temp   = std::move(tail->elem);
            p_tool = tail;
            tail   = tail->prev;
            if (tail) tail->next = nullptr;
            else head = tail;
        } else {
            p_tool = head;
            for(auto i = 0ull; i < idx; ++i) p_tool = p_tool->next;
            temp = std::move(p_tool->elem);
            p_tool->prev->next = p_tool->next;
            p_tool->next->prev = p_tool->prev;
        } else {
            temp   = std::move(head->elem);
            p_tool = head;
            head   = head->next;
            if (head) head->prev = nullptr;
            else tail = head;
        }
        p_tool->prev = nullptr;
        p_tool->next = nullptr;
        delete p_tool; p_tool = nullptr;
        --len;
        return temp;
    }

    net_list sub_list(uint64_t fst_rng, uint64_t snd_rng) const {
        if (fst_rng > snd_rng) std::swap(snd_rng, fst_rng);
        if (snd_rng >= len) return net_list();
        net_list ans;
        auto p_tool = head;
        for (auto i = 0ull; i < fst_rng; ++i) p_tool = p_tool->next;
        for (auto i = fst_rng; i <= snd_rng; ++i) {
            ans.emplace_back(p_tool->elem);
            p_tool = p_tool->next;
        }
        return ans;
    }
    net_list sub_list(const net_set<uint64_t> &idx_set) const {
        auto ptr_temp = idx_set.pointer;
        ptr_dup_remove(ptr_temp.len, ptr_temp.ptr_base, ptr_temp.len);
        net_list ans;
        if (*(ptr_temp.ptr_base + ptr_temp.len - 1) < len) {
            auto p_tool = head;
            auto cnt    = 0ull;
            for (auto i = 0ull; i < ptr_temp.len; ++i) {
                while (cnt != *(ptr_temp.ptr_base + i)) {
                    p_tool = p_tool->next;
                    ++cnt;
                }
                ans.emplace_back(p_tool->elem);
            }
        }
        ptr_temp.reset();
        return ans;
    }

    net_list unit(const net_list &src) const {
        if (src.len && len) {
            auto ans  = *this;
            auto temp = src.head;
            for (auto i = 0ull; i < src.len; ++i)
            {
                ans.emplace_back(temp->elem);
                temp = temp->next;
            }
            return ans;
        } else if (len) return net_list(*this);
        else return net_list(src);
    }

    net_list unit_union(const net_list &src) const {
        if (src.len && len) {
            auto ans  = *this,
                 temp = src;
            auto tool = head;
            for (auto i = 0ull; i < len; ++i) {
                auto tool_temp = temp.head;
                for (auto j = 0ull; j < temp.size(); ++j) 
                    if (tool_temp->elem == tool->elem) {
                        temp.erase(j);
                        break;
                    } else tool_temp = tool_temp->next;
                tool = tool->next;
            }
            return ans.unit(temp);
        } else if (len) return net_list(*this);
        else return net_list(src);
    }

    net_list unit_intersect(const net_list &src) const {
        if (!(len && src.len)) return net_list();
        net_list ans;
        auto temp = src;
        auto tool = head;
        for (auto i = 0ull; i < len; ++i) {
            auto tool_temp = temp.head;
            for (auto j = 0ull; j < temp.size(); ++j)
                if (tool_temp->elem == tool->elem) {
                    ans.emplace_back(temp.erase(j));
                    break;
                } else tool_temp = tool_temp->next;
            tool = tool->next;
        }
        return ans;
    }

    net_set<uint64_t> find(const arg &tgt, uint64_t fst_rng, uint64_t snd_rng) const {
        if (fst_rng > snd_rng) std::swap(fst_rng, snd_rng);
        if (snd_rng > len) return net_set<uint64_t>();
        net_ptr_base<uint64_t> ans_ptr;
        ans_ptr.init(snd_rng - fst_rng + 1);
        auto temp = head;
        auto cnt  = 0ull;
        for (auto i = 0ull; i < fst_rng; ++i) temp = temp->next;
        for (auto i = fst_rng; i <= snd_rng; ++i) {
            if (temp->elem == tgt) *(ans_ptr.ptr_base + cnt++) = i;
            temp = temp->next;
        }
        if (cnt != ans_ptr.len) ptr_alter(ans_ptr.ptr_base, ans_ptr.len, cnt);
        ans_ptr.len = cnt;
        net_set<uint64_t> ans;
        ans.pointer = std::move(ans_ptr);
        return ans;
    }
    net_set<uint64_t> find(const arg &tgt) const { return find(tgt, 0, len - 1); }
    net_set<uint64_t> find(const arg &tgt, const net_set<uint64_t> &idx_set) const {
        auto idx_ptr = idx_set.pointer;
        ptr_dup_remove(idx_ptr.len, idx_ptr.ptr_base, idx_ptr.len);
        if (*(idx_ptr.ptr_base + idx_ptr.len - 1) >= len) return net_set<uint64_t>();
        net_ptr_base<uint64_t> ans_ptr;
        ans_ptr.init(idx_ptr.len);
        ans_ptr.len = 0;
        auto cnt    = 0ull;
        auto temp   = head;
        for (auto i = 0ull; i < idx_ptr.len; ++i) {
            while (cnt != *(idx_ptr.ptr_base + i)) {
                ++cnt;
                temp = temp->next;
            }
            if (temp->elem == tgt) *(ans_ptr.ptr_base + ans_ptr.len++) = cnt;
        }
        if (idx_ptr.len != ans_ptr.len) ptr_alter(ans_ptr.ptr_base, idx_ptr.len, cnt);
        net_set<uint64_t> ans(std::move(ans_ptr));
        idx_ptr.reset();
        return ans;
    }

    bool cut(uint64_t idx, bool successor = true) {
        if (idx >= len) return false;
        if ((idx == 0 && successor) || (idx == len - 1 && !successor)) reset();
        else if (successor) {
            auto tool = tail;
            for (auto i = len - 1; i >= idx; --i) tool = tool->prev;
            delete tool->next;
            tool->next = nullptr;
            tail       = tool;
            len        = idx;
        } else {
            auto tool = head;
            for (auto i = 0; i <= idx; ++i) head = head->next;
            head->prev->next = nullptr;
            head->prev       = nullptr;
            delete tool;
            tool =  nullptr;
            len  -= idx + 1;
        }
        return true;
    }

    void reverse() {
        if (len < 2) return;
        net_node<arg> *p_tool      = tail,
                      *p_tool_prev = nullptr;
        while (p_tool)
        {
            p_tool->next = p_tool->prev;
            p_tool->prev = p_tool_prev;
            p_tool_prev  = p_tool;
            p_tool       = p_tool->next;
        }
        head = tail;
        tail = p_tool_prev;
    }

    iterator begin() const {
        if (len) return iterator(head);
        else return end();
    }

    iterator end() const { return iterator(nullptr); }

    arg sigma() {
        arg ans {};
        if (len) ans = head->elem;
        else return ans;
        for (auto tool = head->next; tool; tool = tool->next) ans += tool->elem;
        return ans;
    }
    arg sigma(uint64_t fst_rng, uint64_t snd_rng) { return sub_list(fst_rng, snd_rng).sigma(); }
    arg sigma(const net_set<uint64_t> &idx_set) { return sub_list(idx_set).sigma(); }

    arg pi() {
        arg ans {};
        if (len) ans = head->elem;
        else return ans;
        for (auto tool = head->next; tool; tool = tool->next) ans *= tool->elem;
        return ans;
    }
    arg pi(uint64_t fst_rng, uint64_t snd_rng) { return sub_list(fst_rng, snd_rng).pi(); }
    arg pi(const net_set<uint64_t> &idx_set) { return sub_list(idx_set).pi(); }

    net_set<arg> set_output() {
        net_set<arg> ans(len);
        auto cnt = 0ull;
        for(auto temp : *this) ans[cnt++] = temp;
        return ans;
    }

    void set_input(const net_set<arg> &src) {
        if (src.length == 0) return;
        if (head == nullptr) head = new net_node<arg>;
        net_node<arg> *temp = nullptr;
        for (auto i = 0; i < src.length; ++i) {
            if (temp == nullptr) temp = head;
            else temp = temp->next;
            temp->elem = src[i];
        }
        tail = temp;
        len  = src.length;
    }

    bool sort(bool asc = true) {
        auto temp = set_output();
        auto flag = temp.sort(asc);
        set_input(temp);
        return flag;
    }

    bool shuffle() {
        if (len < 2) return false;
        auto temp = set_output();
        auto flag = temp.shuffle();
        set_input(temp);
        return flag;
    }

    void supersede(const arg &tgt, const arg &src, net_set<uint64_t> idx_set) {
        auto tgt_idx = find(tgt, idx_set);
        auto cnt     = 0ull;
        auto temp    = head;
        for (auto i = 0ull; i < len; ++i) {
            if (tgt_idx.length == cnt) break;
            if (tgt_idx[cnt] == i) {
                if(tgt == temp->elem) temp->elem = src;
                ++cnt;
            }
            temp = temp->next;
        }
    }
    void supersede(const arg &tgt, const arg &src, uint64_t fst_rng, uint64_t snd_rng) { return supersede(tgt, src, find(tgt, fst_rng, snd_rng)); }
    void supersede(const arg &tgt, const arg &src) { return supersede(tgt, src, 0, len - 1); }

    void reset() {
        delete head;
        head = nullptr; 
        tail = head;
        len  = 0;
    }

    arg &operator[](uint64_t idx) const {
        assert(idx < len);
        if (idx) if ((idx + 1) == len) return tail->elem;
        else {
            auto tool = head;
            for (auto i = 0ull; i < idx; ++i) tool = tool->next;
            return tool->elem;
        } else return head->elem;
    }

    bool operator==(const net_list &val) const {
        if (val.len != len) return false;
        auto temp     = head,
             val_temp = val.head;
        for (auto i = 0ull; i < len; ++i)
            if (temp->elem == val_temp->elem) {
                temp     = temp->next;
                val_temp = val_temp->next;
            } else return false;
        return true;
    }

    bool operator!=(const net_list &val) const { return !(*this == val); }

    net_list &operator=(const net_list &src) { value_copy(src); return *this; }
    net_list &operator=(net_list &&src) { value_move(std::move(src)); return *this; }

    ~net_list() { reset(); }

protected:
    net_node<arg> *head = nullptr,
                  *tail = head;
    uint64_t      len   = 0;
public:
    __declspec(property(get = size))  uint64_t length;
    __declspec(property(get = sigma)) arg      sum;
    __declspec(property(get = pi))    arg      product;

    friend std::ostream& operator<<(std::ostream &out, const net_list &src) {
        out << "[Length " << src.len << "]\n";
        auto p_tool = src.head;
        for (auto i = 0ull; i < src.len; ++i) {
            out << '[' << i << "][\n" << p_tool->elem << "\n]";
            if (i + 1 < src.len) out << '\n';
            p_tool = p_tool->next;
        }
        return out;
    }
};

NEUNET_END