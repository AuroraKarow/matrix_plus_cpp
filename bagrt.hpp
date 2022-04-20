BASEALGO_BEGIN

uint64_t num_cnt(uint64_t first, uint64_t second, uint64_t dilate = 0)
{
    first > second ? std::swap(first, second) : 0;
    auto ans = second - first;
    if(ans % (dilate+1)) return 0;
    else return (ans / (dilate+1) + 1);
}

double num_rate(double numerator, double denominator)
{
    if(denominator) return numerator / denominator;
    else return 0;
}

template<typename _T> _T num_extreme(std::initializer_list<_T> init_num, bool max = true, std::function<bool(_T, _T)>bigger_comp = [](_T _first, _T _second) { return _first > _second; })
{
    auto ext_val = *init_num.begin();
    for(auto temp : init_num)
    {
        if(max && !bigger_comp(ext_val, temp)) ext_val = temp;
        if(!max && bigger_comp(ext_val, temp)) ext_val = temp;
    }
    return ext_val;
}

uint64_t num_pow_pad_cnt(uint64_t val, uint64_t base, uint64_t min_size = 1, uint64_t pad = 0, uint64_t fold = 0)
{
    if(val <= min_size) return pad;
    else if(val % base)
    {
        auto pad_curr = 0;
        while (val % base) { ++ val; ++ pad_curr; }
        pad_curr *= (fold + 1);
        return num_pow_pad_cnt(val/base, base, min_size, pad+pad_curr, (fold+1)*(base-1)+fold);
    }
    else return num_pow_pad_cnt(val/base, base, min_size, pad, (fold+1)*(base-1)+fold);
}

uint64_t num_unsign(uint64_t val) { return val*(-1)<val ? 0 : val; }

template<typename _T> void quick_sort(std::unique_ptr<_T[]> &seq_val, uint64_t begin, uint64_t end, bool asc = true, std::function<bool(_T&, _T&)>func_comp = [](_T &_first, _T &_second){return _first > _second;})
{
    if(end == begin) return;
    else
    {
        auto pivot = begin, slide = end;
        while(slide != pivot)
            if(pivot<slide)
            {
                if(!func_comp(seq_val[slide], seq_val[pivot]))
                {
                    std::swap(seq_val[slide], seq_val[pivot]);
                    std::swap(slide, pivot);
                    slide ++;
                }
                else -- slide;
            }
            else
            {
                if(func_comp(seq_val[slide], seq_val[pivot]))
                {
                    std::swap(seq_val[slide], seq_val[pivot]);
                    std::swap(slide, pivot);
                    -- slide;
                }
                else ++ slide;
            }
        if(pivot != begin) quick_sort(seq_val, begin, pivot-1, asc, func_comp);
        if(pivot != end) quick_sort(seq_val, pivot+1, end, asc, func_comp);
    }
}

void reset_ptr() {}
template<typename T> void reset_ptr(std::unique_ptr<T> &val)
{
    if(val)
    {
        val.reset();
        val.release();
    }
}
template<typename arg, typename...args> void reset_ptr(arg &&first, args &&...others)
{
    reset_ptr(first);
    reset_ptr(others...);
}

template<typename _Ty> class net_queue
{
protected:
    std::unique_ptr<_Ty[]> _ptr;
    uint64_t len = 0;
    void realloc_inc_ptr(uint64_t _size)
    {
        if(_size > len)
        {
            auto p_tool = std::make_unique<_Ty[]>(_size);
            for(auto i=0; i<len; ++i) p_tool[i] = std::move(_ptr[i]);
            reset();
            len = _size;
            _ptr = std::move(p_tool);
        }
    }
    template <typename ... Args> net_queue<_Ty> realloc_dec_ptr(Args &&...idx)
    {
        int idx_arr[] = {idx ...};
        auto cut_len = sizeof(idx_arr)/sizeof(uint64_t);
        auto uniq_idx = std::make_unique<uint64_t[]>(cut_len);
        for(auto i=0; i<cut_len; ++i) if(idx_arr[i] < len) uniq_idx[i] = idx_arr[i];
        else
        {
            reset_ptr(uniq_idx);
            return net_queue<_Ty>::blank_queue();
        }
        net_queue<_Ty> out_temp(cut_len);
        quick_sort(uniq_idx, 0, cut_len - 1);
        auto p_tool = std::make_unique<_Ty[]>(len - cut_len);
        for(auto i=0,re_loc=0,idx_cnt=0; i<len; ++i)
            if(i==uniq_idx[idx_cnt])
            {
                out_temp[idx_cnt++] = _ptr[i];
                continue;
            }
            else p_tool[re_loc++] = _ptr[i];
        reset_ptr(_ptr, uniq_idx);
        _ptr = std::move(p_tool);
        len -= cut_len;
        return out_temp;
    }
    _Ty temp;
public:
    net_queue(uint64_t _size = 0) {realloc_inc_ptr(_size);}
    net_queue(std::initializer_list<_Ty> _init_list)
    {
        init(_init_list.size());
        auto cnt_temp = 0;
        for(auto temp : _init_list) _ptr[cnt_temp++] = std::move(temp);
    }
    void value_copy(net_queue &cpy_val)
    {
        if(cpy_val.len)
        {
            if(len != cpy_val.len)
            {
                reset();
                len = cpy_val.len;
                _ptr = std::make_unique<_Ty []>(len);
            }
            for(auto i=0; i<len; ++i) _ptr[i] = cpy_val._ptr[i];
        }
    }
    void value_move(net_queue &&mov_val)
    {
        reset();
        _ptr = std::move(mov_val._ptr);
        len = mov_val.len;
        mov_val.reset();
    }
    net_queue(net_queue &cpy_val) { value_copy(cpy_val); }
    net_queue(net_queue &&mov_val) { value_move(std::move(mov_val)); }
    bool init(uint64_t _size = 1)
    {
        if(_size > 0)
        {
            reset();
            _ptr = std::make_unique<_Ty[]>(_size);
            len = _size;
            return true;
        }
        else return false;
    }
    static net_queue blank_queue() { return net_queue(); }
    uint64_t size() { return len; }
    template<typename...Args> bool insert(uint64_t idx, Args &&...args)
    {
        if(idx > len) return false;
        else
        {
            realloc_inc_ptr(len + 1);
            for(auto i=len-1; i>idx; ++i) _ptr[i] = _ptr[i-1];
            _ptr[idx] = std::move(_Ty(std::forward<Args>(args)...));
            return true;
        }
    }
    template<typename...Args> bool emplace_back(Args &&...args) { return insert(len, std::forward<Args>(args)...); }
    bool push_back(_Ty val) { return insert(len, val); }
    template<typename... Args> net_queue<_Ty> erase(Args &&...idx) { return realloc_dec_ptr(idx...); }
    net_queue sub_queue(uint64_t idx_first, uint64_t idx_second)
    {
        net_queue sub_out;
        if(idx_first<len && idx_second<len)
        {
            idx_second<idx_first ? std::swap(idx_first, idx_second) : NULL;
            auto mem_size = num_cnt(idx_first, idx_second);
            if(mem_size == len) sub_out = *this;
            else
            {
                sub_out.init(mem_size);
                for(auto i=0; i<mem_size; ++i) sub_out._ptr[i] = _ptr[i+idx_first];
            }
        }
        return sub_out;
    }
    net_queue sub_queue(net_queue<uint64_t> &idx_seq)
    {
        if(idx_seq.size())
        {
            net_queue sub_out(idx_seq.size());
            if(idx_seq.size() == len) sub_out = *this;
            else for(auto i=0; i<sub_out.size(); ++i) sub_out._ptr[i] = _ptr[idx_seq[i]];
            return sub_out;
        }
        else return blank_queue();
    }
    void sort(bool asc = true, std::function<bool(_Ty&, _Ty&)> _func = [](_Ty &_first, _Ty &_second){return _first > _second;}) { quick_sort(_ptr, 0, len - 1, asc, _func); }
    net_queue unit(net_queue &val)
    {

        net_queue u_set(val.len + len);
        for(auto i=0; i<u_set.len; ++i)
        {
            if(i<len) u_set._ptr[i] = _ptr[i];
            else u_set._ptr[i] = val._ptr[i-len];
        }
        return u_set;
    }
    net_queue unit_union(net_queue &val)
    {
        auto src_val = val;
        for(auto i=0; i<len; ++i) for(auto j=0; j<src_val.len; ++j)
            if(_ptr[i] == src_val._ptr[j])
            {
                src_val.erase(j);
                break;
            }
        return unit(src_val);
    }
    net_queue unit_intersect(net_queue &val)
    {
        net_queue itrs;
        auto src_val = val;
        for(auto i=0; i<len; ++i) for(auto j=0; j<src_val.len; ++j)
            if(_ptr[i] == src_val._ptr[j])
            {
                itrs.push_back(src_val._ptr[j]);
                src_val.erase(j);
                break;
            }
        return itrs;
    }
    net_queue<uint64_t> find(_Ty &&target, uint64_t range_first = 0, uint64_t range_second = 0)
    {
        net_queue<uint64_t> idx_set;
        if(range_first<len && range_second<len && range_second!=range_first)
        {
            if(range_first < range_second) std::swap(range_second, range_first);
            for(auto i=range_first; i<=range_second; ++i)
                if(_ptr[i] == target) idx_set.push_back(i);
        }
        return idx_set;
    }
    net_queue<uint64_t> find(_Ty &target, uint64_t range_first = 0, uint64_t range_second = 0) {return find(std::move(target), range_first, range_second);}
    _Ty sum(std::function<_Ty(_Ty&, _Ty&)> add_func = [](_Ty &first, _Ty &second){return first + second;})
    {
        auto rtn_val = _ptr[IDX_ZERO];
        for(auto i=1; i<len; ++i) rtn_val = std::move(add_func(rtn_val, _ptr[i]));
        return rtn_val;
    }
    bool shuffle()
    {
        if(len)
        {
            std::srand((unsigned)std::time(NULL));
            for(auto i=len; i>0; --i) std::swap(_ptr[i-1], _ptr[std::rand()%i]);
            return true;
        }
        else return false;
    }
    _Ty &operator[](uint64_t idx)
    {
        if(idx < len) return _ptr[idx];
        else return temp;
    }
    bool operator==(net_queue &val)
    {
        if(len == val.len)
        {
            for(auto i=0; i<len; ++i) if(val._ptr[i] == _ptr[i]) continue;
            else return false;
            return true;
        }
        else return false;
    }
    bool operator!=(net_queue &val) { return !(val == *this); }
    void operator=(net_queue &val) { value_copy(val); }
    void operator=(net_queue &&val) { value_move(std::move(val)); }
    friend std::ostream& operator<<(std::ostream &output, net_queue &val)
    {
        for(auto i=0; i<val.len; ++i)
        {
            output << '[' << i << "][" << std::endl << val._ptr[i] << std::endl << ']';
            if(i+1 < val.len) output << std::endl;
        }
        return output;
    }
    void reset()
    {
        len = 0;
        reset_ptr(_ptr);
    }
    ~net_queue() { reset(); }
};

template<typename _Ty> class net_sequence
{
protected:
    std::unique_ptr<_Ty[]> p_val;
    uint64_t mem_len = 0, len = 0;
    _Ty temp_elem;
    void value_assign(net_sequence &src) { mem_len = src.mem_len; len = src.len; }
public:
    void value_copy(net_sequence &src)
    {
        if(src.len > mem_len)
        {
            value_assign(src);
            reset_ptr(p_val);
            p_val = std::make_unique<_Ty[]>(mem_len);
        }
        else len = src.len;
        for(auto i=0; i<len; ++i) p_val[i] = src.p_val[i];
    }
    void value_move(net_sequence &&src)
    {
        value_assign(src);
        reset_ptr(p_val);
        p_val = std::move(src.p_val);
        src.reset();
    }
    net_sequence(uint64_t _size = IDX_ZERO, uint64_t alloc_size = IDX_MAX) { init(_size, alloc_size);}
    net_sequence(net_sequence &src) { value_copy(src); }
    net_sequence(net_sequence &&src) { value_move(std::move(src)); }
    net_sequence(std::initializer_list<_Ty> _init_list)
    {
        init(_init_list.size());
        auto cnt_temp = 0;
        for(auto temp : _init_list) p_val[cnt_temp++] = std::move(temp);
    }
    net_sequence blank_sequence() { return net_sequence(IDX_ZERO); }
    uint64_t size() { return len; }
    uint64_t mem_size() { return mem_len; }
    void realloc(uint64_t _mem_size = IDX_MAX, bool append_size = false)
    {
        auto p_tool = std::make_unique<_Ty[]>(_mem_size);
        mem_len = _mem_size;
        if(len > mem_len) len = mem_len;
        for(auto i=0; i<len; ++i) p_tool[i] = std::move(p_val[i]);
        reset_ptr(p_val);
        p_val = std::move(p_tool);
        if(append_size) len = _mem_size;
    }
    bool init(uint64_t _size = 1, uint64_t _alloc_size = IDX_MAX)
    {
        if(_alloc_size < _size) return false;
        else
        {
            if(mem_len != _alloc_size)
            {
                mem_len = _alloc_size;
                reset_ptr(p_val);
                p_val = std::make_unique<_Ty[]>(mem_len);
            }
            else for(auto i=0; i<len; ++i) p_val[i] = _Ty();
            len = _size;
            return true;
        }
    }
    _Ty erase(uint64_t idx)
    {
        if(idx < len)
        {
            temp_elem = std::move(p_val[idx]);
            for(auto i=idx; i<len; ++i) p_val[i] = std::move(p_val[i+1]);
            -- len;
        }
        return temp_elem;
    }
    template<typename ... Args> bool insert(uint64_t idx, Args &&...args)
    {
        if(idx > len) return false;
        else
        {
            if(len == mem_len) realloc(mem_len+IDX_MAX);
            for(auto i=len; i>idx; --i) p_val[i] = std::move(p_val[i-1]);
            p_val[idx] = std::move(_Ty(std::forward<Args>(args)...));
            ++ len;
            return true;
        }
    }
    template<typename ... Args> bool emplace_back(Args &&...args) { return insert(len, std::forward<Args>(args)...); }
    bool push_back(_Ty val) { return insert(len, std::move(val)); }
    net_sequence sub_sequence(uint64_t first_idx, uint64_t second_idx)
    {
        if(first_idx<len && second_idx<len)
        {
            net_sequence temp;
            if(first_idx == second_idx)
            {
                temp.init(IDX_SGL);
                temp.p_val[second_idx] = p_val[second_idx];
                return temp;
            }
            else
            {
                if(first_idx > second_idx) std::swap(first_idx, second_idx);
                temp.init(num_cnt(first_idx, second_idx));
                for(auto i=0; i<temp.size(); ++i) temp.p_val[i] = p_val[i+first_idx];
                return temp;
            }
        }
        else return blank_sequence();
    }
    void sort(bool asc = true, std::function<bool(_Ty&, _Ty&)> _func = [](_Ty &_first, _Ty &_second){return _first > _second;}) { quick_sort(p_val, 0, len - 1, asc, _func); }
    net_sequence unit(net_sequence &val)
    {
        net_sequence u_seq(val.len + len);
        u_seq.len = u_seq.mem_len;
        for(auto i=0; i<u_seq.len; ++i)
        {
            if(i<len) u_seq.p_val[i] = p_val[i];
            else u_seq.p_val[i] = val.p_val[i-len];
        }
        return u_seq;
    }
    net_sequence unit_union(net_sequence &val)
    {
        auto src_val = val;
        for(auto i=0; i<len; ++i) for(auto j=0; j<src_val.len; ++j)
            if(p_val[i] == src_val.p_val[j])
            {
                src_val.erase(j);
                break;
            }
        return unit(src_val);
    }
    net_sequence unit_intersect(net_sequence &val)
    {
        net_sequence itrs;
        auto src_val = val;
        for(auto i=0; i<len; ++i) for(auto j=0; j<src_val.len; ++j)
            if(p_val[i] == src_val.p_val[j])
            {
                itrs.push_back(src_val.p_val[j]);
                src_val.erase(j);
                break;
            }
        return itrs;
    }
    net_sequence<uint64_t> find(_Ty &&target, uint64_t range_first, uint64_t range_second)
    {
        net_sequence<uint64_t> idx_set;
        if(range_first<len && range_second<len && range_second!=range_first)
        {
            if(range_first > range_second) std::swap(range_second, range_first);
            for(auto i=range_first; i<=range_second; ++i)
                if(p_val[i] == target) idx_set.push_back(i);
        }
        return idx_set;
    }
    net_sequence<uint64_t> find(_Ty &&target) { return find(std::move(target), 0, len-1); }
    net_sequence<uint64_t> find(_Ty &target, uint64_t range_first = 0, uint64_t range_second = 0) {return find(std::move(target), range_first, range_second);}
    net_sequence<uint64_t> find(_Ty &target) { return find(target, 0, len-1); }
    _Ty sum(std::function<_Ty(_Ty&, _Ty&)> add_func = [](_Ty &first, _Ty &second){return first + second;})
    {
        auto rtn_val = p_val[IDX_ZERO];
        for(auto i=1; i<len; ++i) rtn_val = std::move(add_func(rtn_val, p_val[i]));
        return rtn_val;
    }
    bool shuffle()
    {
        if(len)
        {
            std::srand((unsigned)std::time(NULL));
            for(auto i=len; i>0; --i) std::swap(p_val[i-1], p_val[std::rand()%i]);
            return true;
        }
        else return false;
    }
    std::unique_ptr<_Ty[]> &get() { return p_val; }
    _Ty &operator[](uint64_t idx)
    {
        if(idx < len) return p_val[idx];
        else return temp_elem;
    }
    bool operator==(net_sequence &val)
    {
        if(len == val.len)
        {
            for(auto i=0; i<len; ++i) if(val.p_val[i] == p_val[i]) continue;
            else return false;
            return true;
        }
        else return false;
    }
    bool operator!=(net_sequence &val) { return !(val == *this); }
    void operator=(net_sequence &val) { value_copy(val); }
    void operator=(net_sequence &&val) { value_move(std::move(val)); }
    friend std::ostream& operator<<(std::ostream &output, net_sequence &val)
    {
        for(auto i=0; i<val.len; ++i)
        {
            output << '[' << i << "][" << std::endl << val.p_val[i] << std::endl << ']';
            if(i+1 < val.len) output << std::endl;
        }
        return output;
    }
    void shrink()
    {
        if(len)
        {
            if(len != mem_len)
            {
                mem_len = len;
                auto p_tool = std::make_unique<_Ty[]>(len);
                for(auto i=0; i<len; ++i) p_tool[i] = std::move(p_val[i]);
                reset_ptr(p_val);
                p_val = std::move(p_tool);
            }
        }
        else reset();
    }
    void reset()
    {
        len = 0;
        mem_len = 0;
        reset_ptr(p_val);
    }
    ~net_sequence() { reset(); }
};

template<typename _Ty> class net_list
{
protected:
    struct node {node *prev = nullptr; std::unique_ptr<node> next = nullptr /* node *next */; _Ty elem;};
    std::unique_ptr<node> head;
    node *tail = nullptr, *itr = nullptr;
    uint64_t len = 0, itr_idx = 0;
    _Ty temp;
    std::unique_ptr<node> create_node() { return std::make_unique<node>(); }
    node *idx_node(uint64_t idx)
    {
        if(idx)
            if(idx+1 == len) itr = tail;
            else
            {
                auto tml_cnt = len - 1 - idx,
                    fnt_cnt = idx;
                uint64_t itr_mov = std::abs((int)itr_idx - (int)idx);
                if(itr && itr_mov<tml_cnt && itr_mov<fnt_cnt) while(itr_mov)
                {
                    if(itr_idx < idx) itr = (itr->next).get();
                    else itr = itr -> prev;
                    -- itr_mov;
                }
                else if(fnt_cnt < itr_mov)
                {
                    itr = head.get();
                    while(fnt_cnt)
                    {
                        itr = (itr->next).get();
                        -- fnt_cnt;
                    }
                }
                else
                {
                    itr = tail;
                    while(tml_cnt)
                    {
                        itr = itr -> prev;
                        -- tml_cnt;
                    }
                }
            }
        else itr = head.get();
        itr_idx = idx;
        return itr;
    }
    void clear()
    {
        while(head)
        {
            auto p_tool = std::move(head.get()->next);
            head.get() -> next = nullptr;
            head.get() -> prev = nullptr;
            reset_ptr(head);
            head = std::move(p_tool);
        }
        tail = nullptr;
        itr = nullptr;
    }
public:
    uint64_t size() { return len; }
    void reset()
    {
        clear();
        len = 0;
    }
    net_list() {}
    bool empty() {return !len;}
    net_list(std::initializer_list<_Ty> src) { for(auto temp : src) emplace_back(std::move(temp)); }
    void value_copy(net_list &src)
    {
        len = src.len;
        if(len)
        {
            clear();
            head = create_node();
            auto src_ptr = src.head.get(), temp_ptr = head.get();
            temp_ptr -> elem = src_ptr -> elem;
            tail = temp_ptr;
            for(auto i=1; i<len; ++i)
            {
                (temp_ptr -> next) = create_node();
                (temp_ptr -> next).get() -> prev = temp_ptr;
                (temp_ptr -> next).get() -> elem = src_ptr -> next -> elem;
                temp_ptr = (temp_ptr -> next).get();
                tail = temp_ptr;
                src_ptr = (src_ptr -> next).get();
            }
        }
    }
    void value_move(net_list &&src)
    {
        len = src.len;
        clear();
        head = std::move(src.head);
        tail = src.tail;
        itr = src.itr;
        itr_idx = src.itr_idx;
        src.reset();
    }
    net_list(net_list &src) { value_copy(src); }
    net_list(net_list &&src) { value_move(std::move(src)); }
    template<typename...args>bool insert(uint64_t idx, args &&...src)
    {
        if(idx > len) return false;
        else
        {
            auto temp_ptr = create_node();
            auto temp_node = temp_ptr.get();
            temp_node -> elem = std::move(_Ty(std::forward<args>(src)...));
            if(idx)
                if(idx == len)
                {
                    if(tail)
                    {
                        tail -> next = std::move(temp_ptr);
                        temp_node -> prev = tail; 
                        tail = temp_node;
                    }
                    else
                    {
                        head = std::move(temp_ptr);
                        tail = head.get();
                    }
                    itr = tail;
                }
                else
                {
                    auto tgt_node = idx_node(idx);
                    temp_node -> prev = tgt_node -> prev;
                    tgt_node -> prev = temp_node;
                    temp_node -> next = std::move(temp_node -> prev -> next);
                    temp_node -> prev -> next = std::move(temp_ptr);
                    itr = temp_node;
                }
            else
            {
                temp_node -> next = std::move(head);
                head = std::move(temp_ptr);
                if(temp_node -> next) (temp_node->next).get() -> prev = temp_node;
                else tail = temp_node;
                itr = head.get();
            }
            ++ len;
            itr_idx = idx;
            return true;
        }
    }
    template<typename...Args> bool emplace_back(Args &&...args) {return insert(len, args...);}
    _Ty &erase(uint64_t idx)
    {
        if(idx < len)
        {
            auto tgt_node = idx_node(idx);
            temp = std::move(tgt_node -> elem);
            -- len;
            if(tgt_node -> prev)
            {
                tgt_node = tgt_node -> prev;
                auto other_ptr = std::move((tgt_node->next).get()->next);
                (tgt_node->next).get() -> prev = nullptr;
                (tgt_node->next).get() -> next = nullptr;
                reset_ptr(tgt_node->next);
                tgt_node -> next = std::move(other_ptr);
                if(tgt_node->next) (tgt_node->next).get()->prev = tgt_node;
                else tail = tgt_node;
                -- itr_idx;
                itr = tgt_node;
            }
            else if(tgt_node -> next)
            {
                auto other_ptr = std::move(head.get()->next);
                other_ptr.get() -> prev = nullptr;
                head.get() -> next = nullptr;
                reset_ptr(head);
                head = std::move(other_ptr);
                itr = head.get();
                itr_idx = 0;
            }
            else reset();
        }
        return temp;
    }
    net_list unit(net_list &src)
    {
        net_list tool_cpy = *this;
        if(src.len)
        {
            net_list src_cpy = src;
            tool_cpy.len += src_cpy.len;
            src_cpy.head.get() -> prev = tool_cpy.tail;
            tool_cpy.tail -> next = std::move(src_cpy.head);
            tool_cpy.tail = src_cpy.tail;
        }
        return tool_cpy;
    }
    net_list unit_union(net_list &src)
    {
        net_list src_cpy = src, tool_cpy = *this;
        for(auto i=0; i<len; ++i) for(auto j=0; j<src_cpy.len; ++j) if(src_cpy[j] == tool_cpy[i])
        {
            src_cpy.erase(j);
            break;
        }
        return tool_cpy.unit(src_cpy);
    }
    net_list unit_intersect(net_list &src)
    {
        net_list src_cpy = src, ls_temp;
        for(auto i=0; i<len; ++i) for(auto j=0; j<src_cpy.len; ++j) if(src_cpy[j] == (*this)[i])
        {
            ls_temp.insert(ls_temp.len, src_cpy.erase(j));
            break;
        }
        return ls_temp;
    }
    _Ty &operator[](uint64_t idx)
    {
        if(idx < len) return idx_node(idx) -> elem;
        else return temp;
    }
    void operator=(net_list &src) { value_copy(src); }
    void operator=(net_list &&src) { value_move(std::move(src)); }
    bool operator==(net_list &src)
    {
        if(len == src.len)
        {
            auto tool_ptr = head.get(),
                src_ptr = src.head.get();
            for(auto i=0; i<len; ++i)
                if(tool_ptr->elem != src_ptr->elem) return false;
                else
                {
                    tool_ptr = (tool_ptr -> next).get();
                    src_ptr = (src_ptr -> next).get();
                }
            return true;
        }
        else return false;
    }
    bool operator!=(net_list &src) { return !(*this == src); }
    friend std::ostream &operator<<(std::ostream &output, net_list &src)
    {
        auto tool_ptr = src.head.get();
        for(auto i=0; i<src.len; ++i)
        {
            output << '[' << i << "][" << std::endl << tool_ptr -> elem << std::endl << ']';
            if(i+1 != src.len) output << std::endl;
            tool_ptr = (tool_ptr -> next).get();
        }
        return output;
    }
    ~net_list() { reset(); }
};

template<typename _Ty> class net_link
{
protected:
    struct node { node *prev = nullptr; std::unique_ptr<node> next = nullptr; _Ty elem; };
    std::unique_ptr<node> head = nullptr;
    node *tail = nullptr;
    uint64_t len = 0;
    _Ty temp;
    std::unique_ptr<node> create_node() { return std::make_unique<node>(); }
    void clear()
    {
        while(head)
        {
            auto p_tool = std::move(head.get()->next);
            head.get() -> next = nullptr;
            head.get() -> prev = nullptr;
            reset_ptr(head);
            head = std::move(p_tool);
        }
        tail = nullptr;
    }
    node *idx_node(uint64_t idx)
    {
        auto cnt = 0;
        auto p_tool = head.get();
        for(auto i=0; i<idx; ++i) p_tool = (p_tool->next).get();
        return p_tool;
    }
public:
    uint64_t size() { return len; }
    void reset() { len = 0; clear(); }
    net_link() {}
    net_link(std::initializer_list<_Ty> src) { for(auto temp : src) emplace_back(std::move(temp)); }
    bool empty() { return !size(); }
    void value_copy(net_link &src)
    {
        len = src.len;
        if(len)
        {
            clear();
            head = create_node();
            auto src_ptr = src.head.get(), temp_ptr = head.get();
            temp_ptr -> elem = src_ptr -> elem;
            tail = temp_ptr;
            for(auto i=1; i<len; ++i)
            {
                (temp_ptr -> next) = create_node();
                (temp_ptr -> next).get() -> prev = temp_ptr;
                (temp_ptr -> next).get() -> elem = src_ptr -> next -> elem;
                temp_ptr = (temp_ptr -> next).get();
                tail = temp_ptr;
                src_ptr = (src_ptr -> next).get();
            }
        }
    }
    void value_move(net_link &&src)
    {
        len = src.len;
        clear();
        head = std::move(src.head);
        tail = src.tail;
        src.reset();
    }
    net_link(net_link &src) { value_copy(src); }
    net_link(net_link &&src) { value_move(std::move(src)); }
    template<typename...args>bool insert(uint64_t idx, args &&...src)
    {
        if(idx > len) return false;
        else
        {
            auto temp_ptr = create_node();
            auto temp_node = temp_ptr.get();
            temp_node -> elem = std::move(_Ty(std::forward<args>(src)...));
            if(idx)
                if(idx == len)
                {
                    if(tail)
                    {
                        tail -> next = std::move(temp_ptr);
                        temp_node -> prev = tail; 
                        tail = temp_node;
                    }
                    else
                    {
                        head = std::move(temp_ptr);
                        tail = head.get();
                    }
                }
                else
                {
                    auto tgt_node = idx_node(idx);
                    temp_node -> prev = tgt_node -> prev;
                    tgt_node -> prev = temp_node;
                    temp_node -> next = std::move(temp_node -> prev -> next);
                    temp_node -> prev -> next = std::move(temp_ptr);
                }
            else
            {
                temp_node -> next = std::move(head);
                head = std::move(temp_ptr);
                if(temp_node -> next) (temp_node->next).get() -> prev = temp_node;
                else tail = temp_node;
            }
            ++ len;
            return true;
        }
    }
    template<typename...Args> bool emplace_back(Args &&...args) {return insert(len, std::forward<Args>(args)...);}
    _Ty erase(uint64_t idx)
    {
        if(idx < len)
        {
            auto tgt_node = idx_node(idx);
            temp = std::move(tgt_node -> elem);
            -- len;
            if(tgt_node -> prev)
            {
                tgt_node = tgt_node -> prev;
                auto other_ptr = std::move((tgt_node->next).get()->next);
                (tgt_node->next).get() -> prev = nullptr;
                (tgt_node->next).get() -> next = nullptr;
                reset_ptr(tgt_node->next);
                tgt_node -> next = std::move(other_ptr);
                if(tgt_node->next) (tgt_node->next).get()->prev = tgt_node;
                else tail = tgt_node;
            }
            else if(tgt_node -> next)
            {
                auto other_ptr = std::move(head.get()->next);
                other_ptr.get() -> prev = nullptr;
                head.get() -> next = nullptr;
                reset_ptr(head);
                head = std::move(other_ptr);
            }
            else reset();
        }
        return temp;
    }
    net_link unit(net_link &src)
    {
        net_link tool_cpy = *this;
        if(src.len)
        {
            net_link src_cpy = src;
            tool_cpy.len += src_cpy.len;
            src_cpy.head.get() -> prev = tool_cpy.tail;
            tool_cpy.tail -> next = std::move(src_cpy.head);
            tool_cpy.tail = src_cpy.tail;
        }
        return tool_cpy;
    }
    net_link unit_union(net_link &src)
    {
        net_link src_cpy = src, tool_cpy = *this;
        for(auto i=0; i<len; ++i) for(auto j=0; j<src_cpy.len; ++j) if(src_cpy[j] == tool_cpy[i])
        {
            src_cpy.erase(j);
            break;
        }
        return tool_cpy.unit(src_cpy);
    }
    net_link unit_intersect(net_link &src)
    {
        net_link src_cpy = src, ls_temp;
        for(auto i=0; i<len; ++i) for(auto j=0; j<src_cpy.len; ++j) if(src_cpy[j] == (*this)[i])
        {
            ls_temp.insert(ls_temp.len, src_cpy.erase(j));
            break;
        }
        return ls_temp;
    }
    _Ty &operator[](uint64_t idx)
    {
        if(idx < len) return idx_node(idx) -> elem;
        else return temp;
    }
    void operator=(net_link &src) { value_copy(src); }
    void operator=(net_link &&src) { value_move(std::move(src)); }
    bool operator==(net_link &src)
    {
        if(len == src.len)
        {
            auto tool_ptr = head.get(),
                src_ptr = src.head.get();
            for(auto i=0; i<len; ++i)
                if(tool_ptr->elem != src_ptr->elem) return false;
                else
                {
                    tool_ptr = (tool_ptr -> next).get();
                    src_ptr = (src_ptr -> next).get();
                }
            return true;
        }
        else return false;
    }
    bool operator!=(net_link &src) { return !(*this == src); }
    friend std::ostream &operator<<(std::ostream &output, net_link &src)
    {
        auto tool_ptr = src.head.get();
        for(auto i=0; i<src.len; ++i)
        {
            output << '[' << i << "][" << std::endl << tool_ptr -> elem << std::endl << ']';
            if(i+1 != src.len) output << std::endl;
            tool_ptr = (tool_ptr -> next).get();
        }
        return output;
    }
    ~net_link() { reset(); }
};

template<typename _K, typename _V>class net_map
{
public:
    struct kv
    {
        _K key; _V value;
        void operator=(kv &val)
        {
            key = val.key;
            value = val.value;
        }
        void operator==(kv &val) {return ((key==val.key) && (value==val.value));}
        bool operator!=(kv &val) {return !(*this == val);}
        friend std::ostream &operator<<(std::ostream &output, kv &out_val)
        {
            output << "[Key]" << std::endl << out_val.key << std::endl << "[Value]" << std::endl << out_val.value;
            return output;
        }
    };
protected:
    net_sequence<kv> val;
    _V v_temp;
    kv kv_temp;
public:
    net_map(uint64_t alloc_size = IDX_MAX) : val(IDX_ZERO, alloc_size) {}
    net_map(net_map &src) { val.value_copy(src.val); }
    net_map(net_map &&src) { val.value_move(std::move(src.val)); }
    void reset() { val.reset(); }
    uint64_t size() { return val.size(); }
    uint64_t mem_size() { return val.mem_size(); }
    int find_idx(_K &&key)
    {
        for(auto i=0; i<size(); ++i) if(val[i].key == key) return i;
        return -1;
    }
    int find_idx(_K &key) { return find_idx(std::move(key)); }
    net_list<_K> find_key(_V &&value)
    {
        net_list<_K> key_set;
        for(auto i=0; i<size(); ++i) if(val[i].value == value) key_set.emplace_back(val[i].key);
        return key_set;
    }
    net_list<_K> find_key(_V &value) { return find_key(std::move(value)); }
    bool insert(_K &&key, _V &&value)
    {
        if(this->find_idx(key)>=0) return false;
        else
        {
            kv in_temp;
            in_temp.key = key;
            in_temp.value = value;
            return val.emplace_back(std::move(in_temp));
        }
    }
    bool insert(_K &key, _V &value) { return insert(std::move(key), std::move(value)); }
    kv erase(_K &&key)
    {
        for(auto i=0; i<size(); ++i) if(val[i].key == key) kv_temp = val.erase(i);
        return kv_temp;
    }
    kv erase(_K &key) { return erase(std::move(key)); }
    kv &index(uint64_t idx)
    {
        if(idx < size()) return val[idx];
        else return kv_temp;
    }
    _V &operator[](_K &&key) 
    {
        auto tar_idx = find_idx(key);
        if(tar_idx < 0) return v_temp;
        else return val[tar_idx].value;
    }
    _V &operator[](_K &key) { return this->operator[](std::move(key)); }
    void operator=(net_map &src) { val.value_copy(src.val); }
    void operator=(net_map &&src) { val.value_move(std::move(src.val)); }
    bool operator==(net_map &val) { return val == val.val; }
    bool operator!=(net_map &val) { return val != val.val; }
    friend std::ostream &operator<<(std::ostream &output, net_map &out_val)
    {
        for(auto i=0; i<out_val.val.size(); ++i)
        {
            output << '[' << i << ']' << out_val.val[i];
            if(i + 1 < out_val.val.size()) output << std::endl;
        }
        return output;
    }
    ~net_map() { reset(); }
};

template<typename k_type = uint64_t> struct clock_timer
{
private:
    struct _dur
    {
        long begin;
        long end;
        long dur; 
    };
    net_map<k_type, _dur> clock_log;
public:
    clock_timer(clock_timer &src) : clock_log(src) {}
    clock_timer(clock_timer &&src) : clock_log(std::move(src)) {}
    void operator=(clock_timer &src) { clock_log = src.clock_log }
    void operator=(clock_timer &&src) { clock_log = std::move(src.clock_log); }

    clock_timer(uint64_t buf_len = IDX_MAX) : clock_log(IDX_MAX) {}
    bool clock_begin(k_type &&clock_log_id)
    {
        auto begin_point = clock();
        auto id_temp = clock_log.find_idx(clock_log_id);
        if(id_temp < 0)
        {
            _dur dur_temp;
            dur_temp.begin = begin_point;
            return clock_log.insert(clock_log_id, dur_temp);
        }
        else
        {
            clock_log.index(id_temp).value.begin = begin_point;
            return true;
        }
    }
    bool clock_begin(k_type &clock_log_id) { return clock_begin(std::move(clock_log_id)); }
    bool clock_end(k_type &&clock_log_id)
    {
        auto end_point = clock();
        auto log_idx = clock_log.find_idx(clock_log_id);
        if(log_idx < 0) return false;
        else
        {
            clock_log.index(log_idx).value.end = end_point;
            clock_log.index(log_idx).value.dur = (clock_log.index(log_idx).value.end - clock_log.index(log_idx).value.begin) * 1000 / CLOCKS_PER_SEC;
        }
        return true;
    }
    bool clock_end(k_type &clock_log_id) { return clock_end(std::move(clock_log_id)); }
    long duration(k_type &&clock_log_id)
    {
        auto log_idx = clock_log.find_idx(clock_log_id);
        if(log_idx < 0) return -1;
        else return clock_log.index(log_idx).value.dur;
    }
    long duration(k_type &clock_log_id) { return duration(std::move(clock_log_id)); }
    void reset() { clock_log.reset(); }
    ~clock_timer() { reset(); }
};

template<typename _Ty> class memory_sequence
{
public:
    struct _pointer
    {
    private:
        _Ty *p_ret; uint64_t len = 0; 
    public:
        _pointer() {}
        _pointer(_pointer &src) { *this = src; }
        _pointer(_pointer &&src) { *this = std::move(src); }
        void operator=(_pointer &src) { p_ret = src.p_ret; len = src.len; }
        void operator=(_pointer &&src) { p_ret = src.p_ret; src.p_ret = nullptr; len = src.len; src.len = 0;}
        _pointer(std::unique_ptr<_Ty[]> &ptr, uint64_t begin_addr, uint64_t _len) : p_ret(ptr.get()+begin_addr), len(_len) {}
        _Ty &operator[](int addr)
        {
            if(addr>=0 && addr<len) return p_ret[addr];
            else throw std::logic_error("Out of boundary!");
        }
        uint64_t size() { return len; }
        std::unique_ptr<_Ty[]> copy_ptr()
        {
            if(len)
            {
                auto p_tool = std::make_unique<_Ty[]>(len);
                for(auto i=0; i<len; ++i) p_tool[i] = p_ret[i];
                return p_tool;
            }
            else return nullptr;
        }
        ~_pointer() { len = 0; p_ret = nullptr; }
    };
protected:
    struct mem_block
    {
        bool occupy = false;
        uint64_t begin_addr = 0, end_addr = 0, size = 0;
        int prev_id = -1, next_id = -1;
        uint64_t length() { return bagrt::num_cnt(begin_addr, end_addr); }
    };
    net_sequence<_Ty> mem_val;
    net_sequence<mem_block> idx_seq;
    uint64_t head_id = 0, rear_id = 0, len = 0;
    void value_assign(memory_sequence &src) { len = src.len; head_id = src.head_id; rear_id = src.rear_id; }
public:
    memory_sequence(uint64_t mem_size = IDX_MAX) : mem_val(mem_size, mem_size) {}
    memory_sequence(memory_sequence &src) { *this = src; }
    memory_sequence(memory_sequence &&src) { *this = std::move(src.idx_seq); }
    void operator=(memory_sequence &src) { mem_val = src.mem_val; idx_seq = src.idx_seq; value_assign(src); }
    void operator=(memory_sequence &&src) { mem_val = std::move(src.mem_val); idx_seq = std::move(src.idx_seq); value_assign(src); }
    uint64_t max_size() { return mem_val.mem_size(); }
    uint64_t size() { return len; }
    uint64_t mem_length() { return (idx_seq[rear_id].end_addr + 1); }
    void print_block_info(int id, bool detail = false)
    {
        if(id>=0 && id<idx_seq.size())
        {
            std::cout << "[ID " << id << "][Size ";
            if(idx_seq[id].occupy) std::cout << idx_seq[id].size;
            else std::cout << 0;
            std::cout << '/' << idx_seq[id].length() << ']';
            if(detail) std::cout << "[Address " << idx_seq[id].begin_addr << '-' << idx_seq[id].end_addr << "][Previous " << idx_seq[id].prev_id << " Next " << idx_seq[id].next_id << ']';
        }
        else std::cerr << "[Illegal ID]" << std::endl;
    }
    void print_mem_info(bool detail = false)
    {
        int tool_idx = head_id;
        std::cout << "[Memory Structure]";
        while(tool_idx >= 0)
        {
            if(tool_idx != head_id) std::cout << "->";
            print_block_info(tool_idx, detail);
            tool_idx = idx_seq[tool_idx].next_id;
        }
        std::cout << std::endl;
        std::cout << "[Memory Size][Block " << size() << '/' << mem_length() << "][Max Size " << max_size() << ']' << std::endl;
    }
    int alloc_mem(uint64_t _size = 1, std::unique_ptr<_Ty[]> &&src = nullptr)
    {
        _size = num_unsign(_size);
        auto tar_idx = -1;
        if(_size)
        {
            for(auto i=0; i<idx_seq.size(); ++i) if(!idx_seq[i].occupy && _size<=idx_seq[i].length())
            {
                tar_idx = i;
                idx_seq[i].occupy = true;
                idx_seq[i].size = _size;
                break;
            }
            if(tar_idx < 0)
            {
                mem_block blk_temp;
                blk_temp.occupy = true;
                blk_temp.begin_addr = 0;
                if(idx_seq.size()) blk_temp.begin_addr = idx_seq[rear_id].end_addr + 1;
                blk_temp.end_addr = blk_temp.begin_addr + _size - 1;
                blk_temp.size = _size;
                tar_idx = idx_seq.size();
                if(tar_idx)
                {
                    blk_temp.prev_id = rear_id;
                    idx_seq[rear_id].next_id = tar_idx;
                }
                rear_id = tar_idx;
                if(mem_val.mem_size() <= blk_temp.end_addr)
                {
                    auto realloc_size = mem_val.mem_size();
                    while((realloc_size+=realloc_size) <= blk_temp.end_addr);
                    mem_val.realloc(realloc_size, true);
                }
                idx_seq.emplace_back(std::move(blk_temp));
            }
            if(src)
            {
                for(auto i=0; i<_size; ++i) mem_val[idx_seq[tar_idx].begin_addr+i] = std::move(src[i]);
                reset_ptr(src);
            }
            len += _size;
        }
        return tar_idx;
    }
    bool free_mem(int &id)
    {
        if(id>=0 && id<idx_seq.size()) if(idx_seq[id].occupy) { idx_seq[id].occupy = false; len -= idx_seq[id].size; id = -1; return true; }
        return false;
    }
    void re_arrange()
    {
        if(head_id!=rear_id)
        {
            int next_id_temp = -1, tool_idx = head_id; len = 0;
            while(tool_idx >= 0)
            {
                next_id_temp = idx_seq[tool_idx].next_id;
                if(idx_seq[tool_idx].occupy)
                {
                    len += idx_seq[tool_idx].size;
                    auto begin_addr_temp = 0;
                    if(idx_seq[tool_idx].prev_id >= 0) begin_addr_temp = idx_seq[idx_seq[tool_idx].prev_id].end_addr + 1;
                    if(begin_addr_temp != idx_seq[tool_idx].begin_addr)
                    {
                        for(auto i=0; i<idx_seq[tool_idx].size; ++i) mem_val[begin_addr_temp+i] = std::move(mem_val[idx_seq[tool_idx].begin_addr+i]);
                        idx_seq[tool_idx].end_addr -= idx_seq[tool_idx].begin_addr - begin_addr_temp;
                        idx_seq[tool_idx].begin_addr = begin_addr_temp;
                    }
                    if(!idx_seq[tool_idx].begin_addr) head_id = tool_idx;
                }
                else if(idx_seq[tool_idx].size)
                {
                    idx_seq[tool_idx].size = 0;
                    idx_seq[idx_seq[tool_idx].prev_id].next_id = idx_seq[tool_idx].next_id;
                    idx_seq[idx_seq[tool_idx].next_id].prev_id = idx_seq[tool_idx].prev_id;
                    idx_seq[tool_idx].prev_id = rear_id;
                    idx_seq[rear_id].next_id = tool_idx;
                    rear_id = tool_idx;
                    idx_seq[tool_idx].next_id = -1;
                }
                else
                {
                    auto blk_size = idx_seq[tool_idx].length();
                    idx_seq[tool_idx].begin_addr = idx_seq[idx_seq[tool_idx].prev_id].end_addr + 1;
                    idx_seq[tool_idx].end_addr = idx_seq[tool_idx].begin_addr + blk_size - 1;
                }
                tool_idx = next_id_temp;
            }
        }
    }
    bool exist(int id) { return (id>=0 && id<idx_seq.size() && idx_seq[id].occupy); }
    void reset() { mem_val.reset(); idx_seq.reset(); }
    _pointer operator[](int id)
    {
        if(id>=0 && id<idx_seq.size() && idx_seq[id].occupy) return _pointer(mem_val.get(), idx_seq[id].begin_addr, idx_seq[id].size);
        else return _pointer(mem_val.get(), 0, 0);
    }
    ~memory_sequence() { reset(); }
};

bool acc_valid(double acc)
{
    while(acc < 1) acc *= 10;
    return acc == 1;
}

double get_rand_base(long long seed, double acc = 1e-5)
{
    if(acc_valid(acc))
    {
        uint64_t expon = 10;
        while(expon * acc < 1) expon *= 10;
        auto base = seed % expon;
        return base * acc;
    }
    else return 0;
}

double random_number(double boundry_first, double boundry_second, bool to_sleep = false, double acc = 1e-5)
{
    if(boundry_first == boundry_second) return boundry_second;
    else if(acc_valid(acc))
    {
        if(to_sleep) _sleep(1);
        long long raw_pt = std::chrono::system_clock::now().time_since_epoch().count();
        auto seed_pt = get_rand_base(raw_pt, acc);
        return (boundry_first<boundry_second) ? (seed_pt*(boundry_second-boundry_first)+boundry_first) : ((seed_pt*(boundry_first-boundry_second))+boundry_second);
    }
    else return 0;
}

net_queue<uint64_t> random_index(uint64_t seq_size, uint64_t amt)
{
    if(seq_size && amt)
    {
        net_queue<uint64_t> idx_seq(seq_size);
        for(auto i=0; i<seq_size; ++i) idx_seq[i] = i;
        idx_seq.shuffle();
        net_queue<uint64_t> idx_seq_shuffle(seq_size);
        for(auto i=0; i<seq_size; ++i) idx_seq_shuffle[i] = i;
        idx_seq_shuffle.shuffle();
        net_queue<uint64_t> ans(amt);
        for(auto i=0; i<ans.size(); ++i) ans[i] = idx_seq[idx_seq_shuffle[i]];
        ans.sort();
        idx_seq.reset();
        idx_seq_shuffle.reset();
        return ans;
    }
    else return net_queue<uint64_t>::blank_queue();
}

uint32_t swap_endian(uint32_t val)
{
	val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF);
	return (val << 16) | (val >> 16);
}

std::pair<uint64_t, uint64_t> sci_num(double dec_part)
{
    dec_part = std::abs(dec_part);
    if(dec_part > 1) return std::make_pair(0, 0);
    else
    {
        uint64_t sum = 0, expon = 1, i_temp = 0;
        auto cnt = 0;
        do
        {
            auto d_temp = (float)dec_part;
            i_temp = (int)(d_temp * 10);
            dec_part *= 10;
            dec_part -= i_temp;
            expon *= 10;
            sum *= 10;
            sum += i_temp;
            ++ cnt;
        } while (i_temp != 0);
        return std::make_pair(sum/10, expon/10);
    }
}

std::string input_paragraph(uint64_t buffer_length = 2000)
{
    std::fflush(stdin);
    std::cout << "Please end submit with double enter." << std::endl;
    std::cout << std::endl; 
    auto cmtr = std::make_unique<char[]>(buffer_length);
    int i = 0;
    char temp, enter; 
    while(i < buffer_length){
        temp = std::getchar();
        if(temp == '\n'){
            enter = std::getchar();
            if(enter == '\n') break;
            else {
                cmtr[i] = '\n';
                cmtr[++i] = enter;
            }
        }else cmtr[i] = temp; i ++;
    }cmtr[i] = '\0';
    auto len = i + 1;
    std::string init_str = "";
    for(i=0; i<len; ++i) init_str += cmtr[i];
    std::fflush(stdin);
    return init_str;
}

std::unique_ptr<double[]> extract_number(std::string &num_str)
{
    int i = 0, cnt = 0, len = num_str.length(), amt = 0;
    // get the amount of the number from string
    while(i < len){
        if(num_str[i]=='.' || num_str[i]=='-' || num_str[i]=='+' || (num_str[i]>='0' && num_str[i]<='9'))
            if(i == 0 || num_str[i-1]==' '|| num_str[i-1]=='\t' || num_str[i-1]=='\n') cnt ++;
        i ++;
    } std::unique_ptr<double[]>num_seq = std::make_unique<double[]>(cnt);
    amt = cnt;
    cnt = 0;
    for(int j=0; j<amt; ++j){
        num_seq[j] = 0.0;
    } i = 0;
    while(i < len){
        double temp = 0.0, div = 0.0, dec = 0.0;
        bool neg = false, di = false;
        /*
        When iterate to a new number, go into the loop
        */
        if(num_str[i]=='.' || num_str[i]=='-' || num_str[i]=='+' || (num_str[i]>='0' && num_str[i]<='9')){
            // looping when meet a punctuated symbol
            bool de = false;
            int de_num = 0;
            while(num_str[i]!=' ' && num_str[i]!='\t' && num_str[i]!='\n' && num_str[i]!='\0'){
                if(num_str[i] == '-'){
                    // when meet a negative symbol sign convert to true
                    neg = true;
                }if(num_str[i]>='0' && num_str[i]<='9'){
                    // if meet decimal sign is true
                    if(de){
                        dec += (num_str[i]-'0')*1.0/std::pow(10, de_num);
                        de_num ++;
                    }else{
                        // if the fractional sign is true
                        if(di) div = div*10.0+(num_str[i]-'0')*1.0;
                        else temp = temp*10.0+(num_str[i]-'0')*1.0;
                    }
                }if(num_str[i] == '.'){
                    // this one is a decimal symbol
                    de = true;
                    de_num ++;
                }if(num_str[i] == '/'){
                    // this one is a fractional symbol
                    di = true;
                    temp += dec;
                    dec = 0.0;
                    de_num = 0;
                    de = false;
                }i ++;
            }if(di){
                div += dec;
                if(neg) div *= -1.0;
                temp /= div;
            }else{
                temp += dec;
                if(neg)temp *= -1.0;
            }
        } num_seq[cnt ++] = temp;
        i ++;
    }
    return num_seq;
}

std::string charset_exchange(std::wstring &str_src){
    auto nLen = str_src.length();
    char *psBuf = new char[nLen + 1];
    wcstombs(psBuf, str_src.c_str(), nLen);
    psBuf[nLen] = '\0';
    std::string strDest(psBuf);
    delete []psBuf;
    psBuf = nullptr;
    if(strDest.length() == nLen) return strDest;
    else return "";
}

std::wstring charset_exchange(std::string str_src){
    // setlocale(LC_ALL, "zh_CN.UTF-8");
    auto nLen = str_src.length();
    wchar_t *pwsBuf = new wchar_t[nLen + 1];
    mbstowcs(pwsBuf, str_src.c_str(), nLen);
    pwsBuf[nLen] = L'\0';
    std::wstring wstrDest(pwsBuf);
    delete []pwsBuf;
    pwsBuf = nullptr;
    if(wstrDest.length() == nLen) return wstrDest;
    else return L"";
}

net_queue<uint64_t> primes(uint64_t upper)
{
    net_queue<uint64_t> prm_set;
    prm_set.emplace_back(2);
    if(upper > 2) for(auto i=3; i<=upper; ++i)
    {
        bool p_flag = true;
        for(auto j=0; j<prm_set.size(); ++j) if(i%prm_set[j] == 0)
        {
            p_flag = false;
            break;
        }
        if(p_flag) prm_set.emplace_back(i);  
    }
    return prm_set;
}

net_queue<uint64_t> primes_fact(uint64_t val)
{
    auto prm_set = primes(val);
    auto itr_cnt = 0;
    net_queue<uint64_t> prm_fact;
    while(val > 1)
    {
        if(val%prm_set[itr_cnt] == 0)
        {
            val /= prm_set[itr_cnt];
            prm_fact.push_back(prm_set[itr_cnt]);
        }
        else ++ itr_cnt;
    }
    return prm_fact;
}

uint64_t greatest_common_divisor(uint64_t l_val, uint64_t r_val)
{
    auto itr_prm = primes_fact(l_val).unit_intersect(primes_fact(r_val));
    auto res = 1;
    for(auto i=0; i<itr_prm.size(); ++i) res *= itr_prm[i];
    return res;
}

uint64_t least_common_multiple(uint64_t l_val, uint64_t r_val)
{
    auto itr_prm = primes_fact(l_val).unit_union(primes_fact(r_val));
    auto res = 1;
    for(auto i=0; i<itr_prm.size(); ++i) res *= itr_prm[i];
    return res;
}

// initializer_list to net_queue
template<typename T> net_queue<T> initilaize_net_queue(std::initializer_list<T> &src)
{
    net_queue<T> cpy(src.size());
    auto cnt = 0;
    for(auto elem : src) cpy[cnt ++] = elem;
    return cpy;
}

template<typename T> net_queue<T> initilaize_net_list(std::initializer_list<T> &src)
{
    net_list<T> cpy;
    for(auto elem : src) cpy.emplace_back(elem);
    return cpy;
}

BASEALGO_END