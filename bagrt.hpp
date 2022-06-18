BAGRT_BEGIN

/* Definition */

/* Number operation */
uint64_t num_cnt(uint64_t first, uint64_t second, uint64_t dilate)
{
    first > second ? std::swap(first, second) : 0;
    auto ans = second - first;
    return (ans / (dilate+1) + 1);
}
uint64_t num_idx_reverse(uint64_t idx, uint64_t len)
{
    if(idx < len) return len - idx - 1;
    else return 0;
}
double num_rate(double numerator, double denominator)
{
    if(denominator) return numerator / denominator;
    else return 0;
}
template<typename arg> arg num_extreme(std::initializer_list<arg> init_num, bool max, std::function<bool(arg, arg)> bigger_comp)
{
    auto ext_val = *init_num.begin();
    for(auto temp : init_num)
    {
        if(max && !bigger_comp(ext_val, temp)) ext_val = temp;
        if(!max && bigger_comp(ext_val, temp)) ext_val = temp;
    }
    return ext_val;
}
uint64_t num_pow_pad_cnt(uint64_t val, uint64_t base, uint64_t min_size, uint64_t pad, uint64_t fold)
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
uint32_t num_swap_endian(uint32_t val)
{
	val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF);
	return (val << 16) | (val >> 16);
}
std::pair<uint64_t, uint64_t> num_sci(double dec_part)
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

/* Pointer */
template<typename arg> bool ptr_copy(arg *&dst, arg *src, uint64_t len)
{
    if(len && src && dst)
    {
        for(auto i=0; i<len; ++i) *(dst + i) = *(src + i);
        // std::memmove(dst, src, sizeof(arg)*len);
        return true;
    }
    else return false;
}
template<typename arg> bool ptr_move(arg *&dst, arg *&&src)
{
    if(src)
    {
        if(dst) { MEM_RECYCLE(dst); }
        dst = std::move(src);
        src = nullptr;
        return true;
    }
    else return false;
}
template<typename arg> bool ptr_insert(arg *&ptr, arg &&src, uint64_t tgt_idx, uint64_t len, bool is_full = true)
{
    if(tgt_idx > len) return false;
    else
    {
        // if(is_full)
        // {
        //     auto ptr_tool = new arg[len+1];
        //     if(tgt_idx)
        //         if(tgt_idx == len) std::memmove(ptr_tool, ptr, sizeof(arg)*(len));
        //         else
        //         {
        //             std::memmove(ptr_tool, ptr, sizeof(arg)*(tgt_idx));
        //             std::memmove(ptr_tool+tgt_idx+1, ptr+tgt_idx, sizeof(arg)*(len-tgt_idx));
        //         }
        //     else std::memmove(ptr_tool+1, ptr, sizeof(arg)*len);
        //     ptr_tool[tgt_idx] = std::move(src);
        //     ptr_move(ptr, std::move(ptr_tool));
        // }
        // else
        // {
        //     if(tgt_idx != len) std::memmove(ptr+tgt_idx+1, ptr+tgt_idx, sizeof(arg)*(len-tgt_idx));
        //     ptr[tgt_idx] = std::move(src);
        // }
        if(is_full)
        {
            auto ptr_tool = new arg[len+1];
            for(auto i=0; i<len; ++i)
                if(i < tgt_idx) *(ptr_tool + i) = std::move(*(ptr+i));
                else *(ptr_tool + i + 1) = std::move(*(ptr+i));
            ptr_move(ptr, std::move(ptr_tool));
        }
        else if(tgt_idx != len) for(auto i=len; i>tgt_idx; --i) *(ptr + i) = std::move(*(ptr+i-1));
        ptr[tgt_idx] = std::move(src);
        return true;
    }
}
template<typename arg> bool ptr_erase(arg *&ptr, uint64_t tgt_idx, uint64_t len, bool is_full = true)
{
    if(tgt_idx<len && ptr)
    {
        // if(is_full)
        // {
        //     auto ptr_tool = new arg[len-1];
        //     if(tgt_idx)
        //         if(tgt_idx == len) std::memmove(ptr_tool, ptr, sizeof(arg)*(len-1));
        //         else
        //         {
        //             std::memmove(ptr_tool, ptr, sizeof(arg)*tgt_idx);
        //             std::memmove(ptr_tool+tgt_idx, ptr+tgt_idx+1, sizeof(arg)*(len-tgt_idx-1));
        //         }
        //     else std::memmove(ptr_tool, ptr+1, sizeof(arg)*(len-1));
        //     ptr_move(ptr, std::move(ptr_tool));
        // }
        // else if(tgt_idx+1 < len) std::memmove(ptr+tgt_idx, ptr+tgt_idx+1, sizeof(arg)*(len-tgt_idx-1));
        if(is_full)
        {
            auto ptr_tool = new arg[len-1];
            for(auto i=0; i<len-1; ++i)
                if(i < tgt_idx) *(ptr_tool + i) = std::move(*(ptr+i));
                else *(ptr_tool + i) = std::move(*(ptr+i+1));
            ptr_move(ptr, std::move(ptr_tool));
        }
        else if(tgt_idx+1 < len) for(auto i=tgt_idx; i<len; ++i) *(ptr + i) = std::move(*(ptr+i+1));
        return true;
    }
    else return false;
}
template<typename arg> bool ptr_concat(arg *&ans, arg *src_first, arg *src_second, uint64_t len_first, uint64_t len_second)
{
    if(ans && src_first && src_second)
    {
        if(len_first) ptr_copy(ans, src_first, len_first);
        auto ans_rear = ans + len_first;
        if(len_second) ptr_copy(ans_rear, src_second, len_second);
        return true;
    }
    else return false;
}
template<typename arg> bool ptr_cut(arg *&ptr, uint64_t tgt_idx, uint64_t len, bool sub_seq = true, bool is_full = true)
{
    if(ptr && tgt_idx<len && len)
    {
        if(is_full)
            if((!tgt_idx&&sub_seq) || (tgt_idx+1==len&&!sub_seq)) { MEM_RECYCLE(ptr); }
            else
            {
                MEM_PTR(arg, ptr_tool);
                if(sub_seq)
                {
                    MEM_ALLOC(arg, ptr_tool, tgt_idx);
                    ptr_copy(ptr_tool, ptr, tgt_idx);
                }
                else
                {
                    auto sub_len = len - tgt_idx;
                    MEM_ALLOC(arg, ptr_tool, sub_len);
                    ptr_copy(ptr_tool, ptr+tgt_idx+1, sub_len-1);
                }
                ptr_move(ptr, std::move(ptr_tool));
            }
        else if(!sub_seq && tgt_idx+1!=len) for(auto i=tgt_idx+1; i<len; ++i) *(ptr + i - tgt_idx - 1) = std::move(*(ptr+i));
        return true;
    }
    return false;
}
template<typename arg> bool ptr_sort(arg *&seq_val, uint64_t begin, uint64_t end, bool asc, std::function<bool(arg&, arg&)> func_comp)
{
    if(end == begin) return true;
    else if(seq_val+begin && seq_val+end)
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
        auto begin_flag = true, end_flag = true;
        if(pivot != begin) begin_flag = ptr_sort(seq_val, begin, pivot-1, asc, func_comp);
        if(pivot != end) end_flag = ptr_sort(seq_val, pivot+1, end, asc, func_comp);
        return (begin_flag && end_flag);
    }
    else return false;
}
template<typename arg> bool ptr_shrink(arg *&src, uint64_t len)
{
    auto p_tool = new arg[len];
    return (ptr_copy(p_tool, src, len) && ptr_move(src, std::move(p_tool)));
}
template<typename arg> bool ptr_extend(arg *&src, uint64_t len, uint64_t ex_len)
{
    if(ex_len == len) return true;
    else if(ex_len > len)
    {
        auto p_tool = new arg[ex_len];
        return (ptr_copy(p_tool, src, len) && ptr_move(src, std::move(p_tool)));
    }
    else return false;
}

/* Hash function */
uint64_t hash_string(std::string &&src)
{
    auto sum = 0;
    for(auto i=0; i<src.length(); ++i) sum += (i+1)*src[i];
    return sum;
}
uint64_t hash_integer(uint64_t &&src) { return src; }
uint64_t hash_float(double &&src)
{
    if(src)
    {
        double a = src;
        uint64_t ans = *(uint64_t*)(&a);
        return ans;
    }
    else return 0;
}
uint64_t __hash_in_build(std::string &&src) { return hash_string(std::move(src)); }
uint64_t __hash_in_build(uint64_t &&src) { return hash_integer(std::move(src)); }
uint64_t __hash_in_build(double &&src) { return hash_float(std::move(src)); }

/* Net sequence */
template<typename _Ty> inline bool net_sequence<_Ty>::__value_copy(const net_sequence &src)
{
    if(mem_len < src.len)
    {
        MEM_RECYCLE(ptr);
        init(src.len, src.mem_len);
    }
    else len = src.len;
    if(len) return ptr_copy(ptr, src.ptr, len);
    else return true;
}
template<typename _Ty> inline net_sequence<_Ty>::net_sequence(uint64_t _size, uint64_t _alloc_size)
{ if(_size && _alloc_size) init(_size, _alloc_size); }
template<typename _Ty> inline net_sequence<_Ty>::net_sequence(std::initializer_list<_Ty> _init_list) : len(_init_list.size()), mem_len(_init_list.size())
{
    MEM_ALLOC(_Ty, ptr, len);
    auto cnt_temp = 0;
    for(auto temp : _init_list) ptr[cnt_temp++] = std::move(temp);
}
template<typename _Ty> inline net_sequence<_Ty>::net_sequence(const net_sequence &src) { __value_copy(src); }
template<typename _Ty> inline net_sequence<_Ty>::net_sequence(net_sequence &src) { value_copy(src); }
template<typename _Ty> inline net_sequence<_Ty>::net_sequence(net_sequence &&src) { value_move(std::move(src)); }
template<typename _Ty> inline bool net_sequence<_Ty>::value_copy(net_sequence &src) { return __value_copy(src); }
template<typename _Ty> inline bool net_sequence<_Ty>::value_move(net_sequence &&src)
{
    len = src.len;
    mem_len = src.mem_len;
    src.len = 0;
    src.mem_len = 0;
    return ptr_move(ptr, std::move(src.ptr));
}
template<typename _Ty> inline void net_sequence<_Ty>::init(uint64_t _size, uint64_t _alloc_size)
{
    while(_size > _alloc_size) _alloc_size += _alloc_size;
    if(mem_len != _alloc_size)
    {
        mem_len = _alloc_size;
        MEM_RECYCLE(ptr);
        MEM_ALLOC(_Ty, ptr, _alloc_size);
    }
    len = _size;
}
template<typename _Ty> inline bool net_sequence<_Ty>::ptr_array(_Ty *&&ptr_array, uint64_t array_len)
{
    if(ptr_array && array_len)
    {
        ptr_move(ptr, std::move(ptr_array));
        len = array_len;
        mem_len = array_len;
        return true;
    }
    else return false;
}
template<typename _Ty> inline bool net_sequence<_Ty>::extend(uint64_t ex_alloc_size, bool extend_size)
{
    if(ptr_extend(ptr, mem_len, ex_alloc_size))
    {
        mem_len = ex_alloc_size;
        if(extend_size) len = ex_alloc_size;
        return true;
    }
    else return false;
}
template<typename _Ty> inline uint64_t net_sequence<_Ty>::size() { return len; }
template<typename _Ty> inline uint64_t net_sequence<_Ty>::mem_size() { return mem_len; }
template<typename _Ty> template<typename ... args> inline bool net_sequence<_Ty>::insert(uint64_t idx, args &&...paras)
{
    if(idx > len) return false;
    else
    {
        if(!mem_len) init(0);
        else if(mem_len == len)
        {
            auto curr_mem_len = mem_len;
            mem_len += mem_len;
            ptr_extend(ptr, curr_mem_len, mem_len);
        }
        ptr_insert(ptr, std::move(_Ty(std::forward<args>(paras)...)), idx, len, false);
        ++ len;
    }
    return true;
}
template<typename _Ty> template<typename ... args> inline bool net_sequence<_Ty>::emplace_back(args &&...paras) { return insert(len, std::forward<args>(paras)...); }
template<typename _Ty> inline bool net_sequence<_Ty>::push_back(_Ty val) { return insert(len, std::move(val)); }
template<typename _Ty> inline _Ty net_sequence<_Ty>::erase(uint64_t idx)
{
    if(idx < len)
    {
        auto temp = std::move(ptr[idx]);
        ptr_erase(ptr, idx, len, false);
        -- len;
        return temp;
    }
    else return _Ty();
}
template<typename _Ty> inline net_sequence<_Ty> net_sequence<_Ty>::sub_sequence(uint64_t range_first, uint64_t range_second)
{
    net_sequence<_Ty> ans(0, 0);
    if(range_first<len && range_second<len)
    {
        if(range_first > range_second) std::swap(range_first, range_second);
        auto elem_cnt = num_cnt(range_first, range_second);
        ans.init(elem_cnt, elem_cnt);
        if(elem_cnt == 1) ans.ptr[0] = ptr[range_first];
        else ptr_copy(ans.ptr, ptr+range_first, elem_cnt);
    }
    return ans;
}
template<typename _Ty> inline net_sequence<_Ty> net_sequence<_Ty>::sub_sequence(net_sequence<uint64_t> &idx_set)
{
    net_sequence<_Ty> ans(0, 0);
    if(idx_set.size())
    {
        ans.init(idx_set.size(), idx_set.size());
        if(len == idx_set.size()) ptr_copy(ans.ptr, ptr, len);
        else for(auto i=0; i<idx_set.size(); ++i) ans.ptr[i] = ptr[idx_set[i]];
    }
    return ans;
}
template<typename _Ty> inline bool net_sequence<_Ty>::sort(bool asc, std::function<bool(_Ty&, _Ty&)> _func) { return ptr_sort(ptr, 0, len-1, asc, _func); }
template<typename _Ty> inline net_sequence<_Ty> net_sequence<_Ty>::unit(net_sequence &src)
{
    if(len && src.len)
    {
        net_sequence<_Ty> ans(len+src.len, len+src.len);
        ptr_concat(ans.ptr, ptr, src.ptr, len, src.len);
        return ans;
    }
    else if(len) return net_sequence<_Ty>(*this);
    else if(src.len) return net_sequence<_Ty>(src);
    else return net_sequence<_Ty>();
}
template<typename _Ty> inline net_sequence<_Ty> net_sequence<_Ty>::unit_union(net_sequence &src)
{
    if(len && src.len)
    {
        net_sequence<_Ty> ans(0, 0);
        if(src.len)
        {
            _Ty *p_tool = nullptr;
            auto len_tool = src.len;
            MEM_ALLOC(_Ty, p_tool, len_tool);
            ptr_copy(p_tool, src.ptr, len_tool);
            for(auto i=0; i<len; ++i) for(auto j=0; j<len_tool; ++j) if(ptr[i] == p_tool[j])
            {
                ptr_erase(p_tool, j, len_tool--);
                break;
            }
            ans.init(len_tool+len, len_tool+len);
            ptr_concat(ans.ptr, ptr, p_tool, len, len_tool);
            MEM_RECYCLE(p_tool);
        }
        else if(len) ans = src;
        return ans;
    }
    else if(len) return net_sequence<_Ty>(*this);
    else if(src.len) return net_sequence<_Ty>(src);
    else return net_sequence<_Ty>();
}
template<typename _Ty> inline net_sequence<_Ty> net_sequence<_Ty>::unit_intersect(net_sequence &src)
{
    if(len && src.len)
    {
        net_sequence<_Ty> ans(0, 0);
        auto len_tool = src.len;
        _Ty *p_tool = nullptr;
        MEM_ALLOC(_Ty, p_tool, len_tool);
        ptr_copy(p_tool, src.ptr, len_tool);
        for(auto i=0; i<len; ++i) for(auto j=0; j<len_tool; ++j) if(ptr[i] == p_tool[j])
        {
            ptr_insert(ans.ptr, std::move(p_tool[j]), ans.mem_len, ans.mem_len++);
            ptr_erase(p_tool, j, len_tool--);
            break;
        }
        MEM_RECYCLE(p_tool);
        ans.len = ans.mem_len;
        return ans;
    }
    else return net_sequence<_Ty>();
}
template<typename _Ty> inline net_sequence<uint64_t> net_sequence<_Ty>::find(_Ty &&target, uint64_t range_first, uint64_t range_second)
{
    net_sequence<uint64_t> ans(0, 0);
    if(range_first<len && range_second<len)
    {
        if(range_first > range_second) std::swap(range_first, range_second);
        for(auto i=range_first; i<=range_second; ++i) if(ptr[i] == target) ans.emplace_back(i);
    }
    return ans;
}
template<typename _Ty> inline net_sequence<uint64_t> net_sequence<_Ty>::find(_Ty &target, uint64_t range_first, uint64_t range_second) { return find(std::move(target), range_first, range_second); }
template<typename _Ty> inline _Ty net_sequence<_Ty>::sum(std::function<_Ty(_Ty&, _Ty&)> add_func)
{
    if(len)
    {
        auto ans = ptr[0];
        for(auto i=1; i<len; ++i) ans = add_func(ans, ptr[i]);
        return ans;
    }
    else return _Ty();
}
template<typename _Ty> inline void net_sequence<_Ty>::fill_with(_Ty &&src) { for(auto i=0; i<len; ++i) ptr[i] = src; }
template<typename _Ty> inline void net_sequence<_Ty>::memory_set(_Ty &&src) { std::fill_n(ptr, mem_len, src); }
template<typename _Ty> inline bool net_sequence<_Ty>::alter_batch(_Ty &&target, _Ty &&src)
{
    auto tgt_idx = find(target, 0, len-1);
    if(tgt_idx.size())
    {
        for(auto i=0; i<tgt_idx.size(); ++i) ptr[tgt_idx[i]] = src;
        return true;
    }
    else return false;
}
template<typename _Ty> inline bool net_sequence<_Ty>::alter_batch(_Ty &target, _Ty &src) { return alter_batch(std::move(target), std::move(src)); }
template<typename _Ty> inline bool net_sequence<_Ty>::shuffle()
{
    if(len)
    {
        std::srand((unsigned)std::time(NULL));
        for(auto i=len; i>0; --i) std::swap(ptr[i-1], ptr[std::rand()%i]);
        return true;
    }
    else return false;
}
template<typename _Ty> inline bool net_sequence<_Ty>::shrink()
{
    if(len)
        if(len != mem_len)
        {
            mem_len = len;
            return ptr_shrink(ptr, len);
        }
        else return true;
    else
    {
        reset();
        return true;
    }
}
template<typename _Ty> inline bool net_sequence<_Ty>::reverse()
{
    if(len)
    {
        for(int i=0,j=len-1; i<j; ++i,--j)
        {
            auto temp = std::move(ptr[i]);
            ptr[i] = std::move(ptr[j]);
            ptr[j] = std::move(temp);
        }
        return true;
    }
    else return false;
}
template<typename _Ty> inline bool net_sequence<_Ty>::cut(uint64_t idx, bool sub_seq)
{
    if(idx < len)
        if(sub_seq)
        {
            len = idx;
            return true;
        }
        else
        {
            auto flag =  ptr_cut(ptr, idx, len, sub_seq, false);
            len -= idx + 1;
            return flag;
        }
    else return false;
}
template<typename _Ty> inline _Ty &net_sequence<_Ty>::operator[](uint64_t idx)
{
    assert(idx < len);
    return ptr[idx];
}
template<typename _Ty> inline bool net_sequence<_Ty>::operator==(net_sequence &val)
{
    if(len == val)
    {
        for(auto i=0; i<len; ++i) if(val.ptr[i] != ptr[i]) return false;
        return true;
    }
    else return false;
}
template<typename _Ty> inline bool net_sequence<_Ty>::operator!=(net_sequence &val) { return !(*this == val); }
template<typename _Ty> inline void net_sequence<_Ty>::operator=(net_sequence &src) { value_copy(src); }
template<typename _Ty> inline void net_sequence<_Ty>::operator=(const net_sequence &src) { __value_copy(src); }
template<typename _Ty> inline void net_sequence<_Ty>::operator=(net_sequence &&src) { value_move(std::move(src)); }
template<typename _Ty> inline void net_sequence<_Ty>::clear() { len = 0; }
template<typename _Ty> inline void net_sequence<_Ty>::reset() { clear(); mem_len = 0; MEM_RECYCLE(ptr); }
template<typename _Ty> inline net_sequence<_Ty>::~net_sequence() { reset(); }

/* Net set */
template<typename _Ty> inline net_set<_Ty>::net_set(uint64_t _size) : net_sequence<_Ty>(_size, _size) {}
template<typename _Ty> inline net_set<_Ty>::net_set(std::initializer_list<_Ty> _init_list) : net_sequence<_Ty>(_init_list) {}
template<typename _Ty> inline net_set<_Ty>::net_set(const net_set &src) : net_sequence<_Ty>(src) {}
template<typename _Ty> inline net_set<_Ty>::net_set(net_set &src) : net_sequence<_Ty>(src) {}
template<typename _Ty> inline net_set<_Ty>::net_set(net_set &&src) : net_sequence<_Ty>(std::move(src)) {}
template<typename _Ty> inline void net_set<_Ty>::init(uint64_t _size) { net_sequence<_Ty>::init(_size, _size); }
template<typename _Ty> inline bool net_set<_Ty>::operator==(net_set<_Ty> &val) { net_sequence<_Ty>::operator==(val); }
template<typename _Ty> inline bool net_set<_Ty>::operator!=(net_set<_Ty> &val) { net_sequence<_Ty>::operator!=(val); }
template<typename _Ty> inline void net_set<_Ty>::operator=(net_set<_Ty> &src) { net_sequence<_Ty>::operator=(src); }
template<typename _Ty> inline void net_set<_Ty>::operator=(const net_set<_Ty> &src) { net_sequence<_Ty>::operator=(src); }
template<typename _Ty> inline void net_set<_Ty>::operator=(net_set<_Ty> &&src) { net_sequence<_Ty>::operator=(std::move(src)); }
template<typename _Ty> inline net_set<_Ty>::~net_set() {}

/* Net list */
/* node */
template<typename _Ty> inline net_list<_Ty>::net_node::~net_node()
{
    prev = nullptr;
    while(next)
    {
        delete next;
        next = nullptr;
    }
}
template<typename _Ty> inline net_list<_Ty>::net_node_elem::net_node_elem(net_node *_curr_node) : curr_node(_curr_node) {}
template<typename _Ty> inline net_list<_Ty>::net_node_elem::net_node_elem(net_node_elem &src) : curr_node(src.curr_node) {}
template<typename _Ty> inline net_list<_Ty>::net_node_elem::net_node_elem(const net_node_elem &src) : curr_node(src.curr_node) {}
template<typename _Ty> inline bool net_list<_Ty>::net_node_elem::is_null() { return (curr_node == nullptr); }
template<typename _Ty> inline void net_list<_Ty>::net_node_elem::swap_node(net_node_elem &src)
{
    auto p_tool = src.curr_node;
    src.curr_node = curr_node;
    curr_node = p_tool;
}
template<typename _Ty> inline const _Ty net_list<_Ty>::net_node_elem::_node_elem()
{
    if(curr_node)
    {
        const _Ty temp = curr_node->elem;
        return temp;
    }
    else
    {
        const _Ty temp = _Ty();
        return temp;
    }
}
template<typename _Ty> inline bool net_list<_Ty>::net_node_elem::_node_elem(_Ty &&src)
{
    if(curr_node)
    {
        curr_node->elem = std::move(src);
        return true;
    }
    else return false;
}
template<typename _Ty> inline bool net_list<_Ty>::net_node_elem::_node_elem(_Ty &src)
{
    auto temp = src;
    return _node_elem(std::move(temp));
}
template<typename _Ty> inline typename net_list<_Ty>::net_node_elem net_list<_Ty>::net_node_elem::_prev_node() { return net_list<_Ty>::net_node_elem(curr_node->prev); }
template<typename _Ty> inline typename net_list<_Ty>::net_node_elem net_list<_Ty>::net_node_elem::_next_node() { return net_list<_Ty>::net_node_elem(curr_node->next); }
template<typename _Ty> inline void net_list<_Ty>::net_node_elem::operator=(net_node_elem &src) { curr_node = src.curr_node; }
template<typename _Ty> inline void net_list<_Ty>::net_node_elem::operator=(const net_node_elem &src) { curr_node = src.curr_node; }
template<typename _Ty> inline bool net_list<_Ty>::net_node_elem::operator==(net_node_elem &src) { return (curr_node == src.curr_node); }
template<typename _Ty> inline bool net_list<_Ty>::net_node_elem::operator!=(net_node_elem &src) { return (curr_node != src.curr_node); }
template<typename _Ty> inline net_list<_Ty>::net_node_elem::~net_node_elem() { curr_node = nullptr; }
/* list */
template<typename _Ty> inline void net_list<_Ty>::__value_copy(const net_list &src)
{
    if(src.head)
    {
        if(!head)
        {
            head = new net_node;
            head->next = nullptr;
            head->prev = nullptr;
            tail = head;
        }
        head->elem = src.head->elem;
        auto p_tool = head, src_tool = src.head;
        while (src_tool->next)
        {
            if(!p_tool->next)
            {
                p_tool->next = new net_node;
                p_tool->next->prev = p_tool;
            }
            p_tool->next->elem = src_tool->next->elem;
            p_tool = p_tool->next;
            src_tool = src_tool->next;
        }
        tail = p_tool;
        p_tool = p_tool->next;
        tail->next = nullptr;
        delete p_tool;
        p_tool = nullptr;
        len = src.len;
    }
    else reset();
}
template<typename _Ty> inline net_list<_Ty>::net_list() {}
template<typename _Ty> inline net_list<_Ty>::net_list(std::initializer_list<_Ty> _init_list) { for(auto temp : _init_list) emplace_back(std::move(temp)); }
template<typename _Ty> inline net_list<_Ty>::net_list(net_list &src) { value_copy(src); }
template<typename _Ty> inline net_list<_Ty>::net_list(const net_list &src) { __value_copy(src); }
template<typename _Ty> inline net_list<_Ty>::net_list(net_list &&src) { value_move(std::move(src)); }
template<typename _Ty> inline void net_list<_Ty>::value_copy(net_list &src) { __value_copy(src); }
template<typename _Ty> inline void net_list<_Ty>::value_move(net_list &&src)
{
    if(len) delete head;
    head = std::move(src.head); src.head = nullptr;
    tail = src.tail; src.tail = nullptr;
    len = src.len;
    src.len = 0;
}
template<typename _Ty> inline uint64_t net_list<_Ty>::size() { return len; }
template<typename _Ty> template<typename ... args> inline bool net_list<_Ty>::insert(uint64_t idx, args &&...paras)
{
    if(idx > len) return false;
    else
    {
        auto i_node = new net_node;
        i_node->elem = std::move(_Ty(std::forward<args>(paras)...));
        i_node->next = nullptr;
        i_node->prev = nullptr;
        if(len)
            if(idx)
                if(idx == len)
                {
                    tail->next = i_node;
                    i_node->prev = tail;
                    tail = i_node;
                }
                else
                {
                    auto p_tool = head;
                    for(auto i=0; i<idx; ++i) p_tool = p_tool->next;
                    i_node->next = p_tool;
                    p_tool->prev->next = i_node;
                    i_node->prev = p_tool->prev;
                    p_tool->prev = i_node;
                }
            else
            {
                i_node->next = head;
                head->prev = i_node;
                head = i_node;
            }
        else 
        {
            head = i_node;
            tail = head;
        }
        ++ len;
        return true;
    }
}
template<typename _Ty> template<typename ... args> inline bool net_list<_Ty>::emplace_back(args &&...paras) { return insert(len, std::forward<args>(paras)...); }
template<typename _Ty> inline bool net_list<_Ty>::push_back(_Ty val) { return insert(len, std::move(val)); }
template<typename _Ty> inline _Ty net_list<_Ty>::erase(uint64_t idx)
{
    _Ty temp;
    if(idx < len)
    {
        net_node *p_tool = nullptr;
        if(idx)
            if(idx+1 == len)
            {
                temp = std::move(tail->elem);
                p_tool = tail;
                tail = tail->prev;
                if(tail) tail->next = nullptr;
                else head = tail;
            }
            else
            {
                p_tool = head;
                for(auto i=0; i<idx; ++i) p_tool = p_tool->next;
                temp = std::move(p_tool->elem);
                p_tool->prev->next = p_tool->next;
                p_tool->next->prev = p_tool->prev;
            }
        else
        {
            temp = std::move(head->elem);
            p_tool = head;
            head = head->next;
            if(head) head->prev = nullptr;
            else tail = head;
        }
        p_tool->prev = nullptr;
        p_tool->next = nullptr;
        delete p_tool; p_tool = nullptr;
        -- len;
    }
    return temp;
}
template<typename _Ty> inline net_list<_Ty> net_list<_Ty>::sub_list(uint64_t range_first, uint64_t range_second)
{
    net_list<_Ty> ans;
    if(range_first<len && range_second<len)
    {
        auto p_tool = head;
        if(range_first > range_second) std::swap(range_first, range_second);
        for(auto i=0; i<range_first; ++i) p_tool = p_tool->next;
        for(auto i=range_first; i<=range_second; ++i)
        {
            ans.emplace_back(p_tool->elem);
            p_tool = p_tool->next;
        }
    }
    return ans;
}
template<typename _Ty> inline net_list<_Ty> net_list<_Ty>::sub_list(net_sequence<uint64_t> &idx_set)
{
    net_list<_Ty> ans;
    if(idx_set.size())
    {
        for(auto i=0; i<idx_set.size(); ++i) if(idx_set[i] >= len) return ans;
        idx_set.sort();
        net_node *p_tool = head;
        auto idx_cnt = 0;
        for(auto i=0; i<idx_set.size(); ++i)
        {
            while(idx_cnt++ != idx_set[i]) p_tool = p_tool->next;
            ans.emplace_back(p_tool->elem);
            p_tool = p_tool->next;
        }
    }
    return ans;
}
template<typename _Ty> inline net_list<_Ty> net_list<_Ty>::unit(net_list &src)
{
    if(len && src.len)
    {
        auto ans = *this;
        if(src.len)
        {
            auto p_tool = head, src_tool = src.head;
            for(auto i=0; i<src.len; ++i)
            {
                ans.emplace_back(src_tool->elem);
                p_tool = p_tool->next;
                src_tool = src_tool->next;
            }
        }
        return ans;
    }
    else if(len) return *this;
    else if(src.len) return src;
    else return net_list<_Ty>();
}
template<typename _Ty> inline net_list<_Ty> net_list<_Ty>::unit_union(net_list &src)
{
    if(len && src)
    {
        auto ans = *this, src_temp = src;
        auto p_tool = ans.head, src_tool = src_temp.head;
        for(auto i=0; i<len; ++i)
        {
            for(auto j=0; j<src_temp.len; ++j)
                if(src_tool->elem == p_tool->elem)
                {
                    src_temp.erase(j);
                    break;
                }
                else src_tool = src_tool->next;
            p_tool = p_tool->next;
            src_tool = src_temp.head;
        }
        return ans.unit(src_temp);
    }
    else if(len) return net_sequence<_Ty>(*this);
    else if(src.len) return net_sequence<_Ty>(src);
    else return net_list<_Ty>();
}
template<typename _Ty> inline net_list<_Ty> net_list<_Ty>::unit_intersect(net_list &src)
{
    if(len && src)
    {
        auto src_temp = src, ans;
        auto p_tool = head, src_tool = src_temp.head;
        for(auto i=0; i<len; ++i)
        {
            for(auto j=0; j<src_temp.len; ++j)
                if(src_tool->elem == p_tool->elem)
                {
                    ans.emplace_back(std::move(src_temp.erase(j)));
                    break;
                }
                else src_tool = src_tool->next;
            p_tool = p_tool->next;
            src_tool = src_temp.head;
        }
        return ans;
    }
    else return net_list<_Ty>();
}
template<typename _Ty> inline net_sequence<uint64_t> net_list<_Ty>::find(_Ty &&target, uint64_t range_first, uint64_t range_second)
{
    net_list<uint64_t> ans_ls;
    if(range_first<len && range_second<len)
    {
        auto p_tool = head;
        if(range_first > range_second) std::swap(range_first, range_second);
        for(auto i=0; i<range_first; ++i) p_tool = p_tool->next;
        for(auto i=range_first; i<=range_second; ++i)
        {
            if(p_tool->elem == target) ans_ls.emplace_back(i);
            p_tool = p_tool->next;
        }
    }
    return ans_ls.convert_sequence();
}
template<typename _Ty> inline net_sequence<uint64_t> net_list<_Ty>::find(_Ty &target, uint64_t range_first, uint64_t range_second) { return find(std::move(target), range_first, range_second); }
template<typename _Ty> inline _Ty net_list<_Ty>::sum(std::function<_Ty(_Ty&, _Ty&)> add_func)
{
    _Ty ans;
    auto p_tool = head;
    for(auto i=0; i<len; ++i)
    {
        ans = add_func(ans, p_tool->elem);
        p_tool = p_tool->next;
    }
    return ans;
}
template<typename _Ty> inline bool net_list<_Ty>::alter_batch(_Ty &&target, _Ty &&src)
{
    auto p_tool = head;
    auto alter_flag = false;
    for(auto i=0; i<len; ++i)
    {
        if(p_tool->elem == target)
        {
            p_tool->elem = src;
            alter_flag = true;
        }
        p_tool = p_tool->next;
    }
    return alter_flag;
}
template<typename _Ty> inline bool net_list<_Ty>::alter_batch(_Ty &target, _Ty &src) { alter_batch(std::move(target), std::move(src)); }
template<typename _Ty> inline typename net_list<_Ty>::net_node_elem net_list<_Ty>::head_node() { return net_node_elem(head); }
template<typename _Ty> inline typename net_list<_Ty>::net_node_elem net_list<_Ty>::tail_node() { return net_node_elem(tail); }
template<typename _Ty> inline net_sequence<_Ty> net_list<_Ty>::convert_sequence()
{
    net_sequence<_Ty> ans(0, 0);
    if(len)
    {
        ans.init(len, len);
        auto p_tool = head;
        for(auto i=0; i<len; ++i)
        {
            ans[i] = p_tool->elem;
            p_tool = p_tool->next;
        }
    }
    return ans;
}
template<typename _Ty> inline bool net_list<_Ty>::cut(uint64_t idx)
{
    if(idx < len)
    {
        if(idx)
        {
            auto p_tool = head;
            for(auto i=1; i<idx; ++i) p_tool = p_tool->next;
            delete p_tool->next;
            p_tool->next = nullptr;
            tail = p_tool;
            len = idx;
        }
        else reset();
        return true;
    }
    else return false;
}
template<typename _Ty> inline void net_list<_Ty>::reverse()
{
    net_node *p_tool = tail,
            *p_tool_prev = nullptr;
    while(p_tool)
    {
        p_tool->next = p_tool->prev;
        p_tool->prev = p_tool_prev;
        p_tool_prev = p_tool;
        p_tool = p_tool->next;
    }
    head = tail;
    tail = p_tool_prev;
}
template<typename _Ty> inline _Ty &net_list<_Ty>::operator[](uint64_t idx)
{
    assert(idx < len);
    auto p_tool = head;
    for(auto i=0; i<idx; ++i) p_tool = p_tool->next;
    return p_tool->elem;
}
template<typename _Ty> inline bool net_list<_Ty>::operator==(net_list &val)
{
    if(val.len == len)
    {
        auto p_tool = head, src_tool = val.head;
        for(auto i=0; i<len; ++i) if(p_tool->elem != src_tool->elem) return false;
        return true;
    }
    else return false;
}
template<typename _Ty> inline bool net_list<_Ty>::operator!=(net_list &val) { return !(*this == src); }
template<typename _Ty> inline void net_list<_Ty>::operator=(net_list &src) { value_copy(src); }
template<typename _Ty> inline void net_list<_Ty>::operator=(const net_list &src) { __value_copy(src); }
template<typename _Ty> inline void net_list<_Ty>::operator=(net_list &&src) { value_move(std::move(src)); }
template<typename _Ty> inline void net_list<_Ty>::reset() { delete head; len = 0; head = nullptr; tail = head;}
template<typename _Ty> inline net_list<_Ty>::~net_list() { reset(); }

/* Net map */
/* kv */
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::net_kv::net_kv() {}
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::net_kv::net_kv(const net_kv &src) { *this = src; }
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::net_kv::net_kv(net_kv &src) { *this = src; }
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::net_kv::net_kv(net_kv &&src) { *this = std::move(src); }
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::net_kv::operator=(const net_kv &src) { key = src.key; value = src.value; }
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::net_kv::operator=(net_kv &src) { key = src.key; value = src.value; }
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::net_kv::operator=(net_kv &&src) { key = std::move(src.key); value = std::move(src.value); }
template<typename _Kty, typename _Ty> inline bool net_map<_Kty, _Ty>::net_kv::operator==(net_kv &src) { return (key==src.key && value==src.value); }
template<typename _Kty, typename _Ty> inline bool net_map<_Kty, _Ty>::net_kv::operator!=(net_kv &src) { return !(*this == src); }
/* map */
template<typename _Kty, typename _Ty> inline long long net_map<_Kty, _Ty>::detective(long long &threshold, bool &sgn)
{
    if(threshold)
        if(sgn) { sgn = false; auto temp = (-1) * threshold * threshold; ++ threshold; return temp; }
        else { sgn = true; return threshold * threshold; }
    else
    {
        auto temp = threshold;
        ++ threshold;
        sgn = false;
        return temp;
    }
}
template<typename _Kty, typename _Ty> inline long long net_map<_Kty, _Ty>::next_key(long long hash_key, long long d, long long curr_mem_len) { return (hash_key + d) % curr_mem_len; }
template<typename _Kty, typename _Ty> inline uint64_t net_map<_Kty, _Ty>::_find_idx(_Kty &&_key, bool lex_idx)
{
    long long threshold = 0, idx = hash_func(std::move(_key));
    auto sgn = false;
    while (true)
    {
        idx = next_key(idx, detective(threshold, sgn), kv_data[lex_idx].size());
        if(idx>=0 && (!kv_occupy[lex_idx][idx]||(kv_occupy[lex_idx][idx]&&kv_data[lex_idx][idx].key==_key))) break;
    }
    return idx;
}
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::_rehash()
{
    int data_size = rehash_load, data_cnt = 0;
    if(len[!backup] < data_size) data_size = len[!backup];
    for(auto i=0; i<kv_occupy[!backup].size(); ++i) if(kv_occupy[!backup][i])
    {
        auto tgt_idx = _find_idx(std::move(kv_data[!backup][i].key), backup);
        if(!kv_occupy[backup][tgt_idx])
        {
            kv_data[backup][tgt_idx].key = kv_data[!backup][i].key;
            kv_data[backup][tgt_idx].value = kv_data[!backup][i].value;
            kv_occupy[backup][tgt_idx] = true;
            ++ len[backup];
            ++ data_cnt;
        }
        if(data_cnt == data_size) break;
    }
    if(len[backup] == len[!backup])
    {
        auto mem_len = kv_data[backup].mem_size();
        mem_len += mem_len;
        kv_data[!backup].reset();
        kv_data[!backup].init(mem_len, mem_len);
        kv_occupy[!backup].reset();
        kv_occupy[!backup].init(mem_len, mem_len);
        kv_occupy[!backup].memory_set(false);
        len[!backup] = 0;
        backup = !backup;
        rehash_load = base_load;
    }
}
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::__value_copy(const net_map &src)
{
    hash_factor = src.hash_factor;
    rehash_load = src.rehash_load;
    base_load = src.base_load;
    len[0] = src.len[0]; len[1] = src.len[1];
    backup = src.backup;
    hash_func = src.hash_func;
    kv_data[0] = src.kv_data[0]; kv_data[1] = src.kv_data[1];
    kv_occupy[0] = src.kv_occupy[0]; kv_occupy[1] = src.kv_occupy[1];
}
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::net_map(std::function<uint64_t(_Kty&&)> hash_key_func, uint64_t _alloc_size, uint64_t _rehash_load, double _hash_factor) : hash_func(hash_key_func), hash_factor(_hash_factor), rehash_load(_rehash_load), base_load(_rehash_load)
{
    uint64_t base_mem_len = _alloc_size/_hash_factor;
    kv_data[!backup].init(base_mem_len, base_mem_len);
    kv_data[backup].init(base_mem_len+base_mem_len, base_mem_len+base_mem_len);
    kv_occupy[!backup].init(base_mem_len, base_mem_len);
    kv_occupy[!backup].fill_with(false);
    kv_occupy[backup].init(base_mem_len+base_mem_len, base_mem_len+base_mem_len);
    kv_occupy[backup].fill_with(false);
    len[0] = 0; len[1] = 0;
}
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::net_map(std::initializer_list<_Kty> k_init_list, std::initializer_list<_Ty> v_init_list, std::function<uint64_t(_Kty&&)> hash_key_func, uint64_t _alloc_size, uint64_t _rehash_load, double _hash_factor)
{
    assert(k_init_list.size() >= v_init_list.size());
    *this = net_map(hash_key_func, _alloc_size, _rehash_load, _hash_factor);
    auto k_init_head = k_init_list.begin();
    auto v_init_head = v_init_list.begin();
    while(v_init_head != v_init_list.end())
    {
        auto k_temp = *k_init_head;
        auto v_temp = *v_init_head;
        assert(insert(std::move(k_temp), std::move(v_temp)));
        k_init_head ++;
        v_init_head ++;
    }
    while(k_init_head != k_init_list.end())
    {
        auto k_temp = *k_init_head;
        assert(insert(std::move(k_temp), _Ty()));
        k_init_head ++;
    }
}
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::net_map(net_map &src) { value_copy(src); }
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::net_map(const net_map &src) { __value_copy(src); }
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::net_map(net_map &&src) { value_move(std::move(src)); }
template<typename _Kty, typename _Ty> inline bool net_map<_Kty, _Ty>::set_hash_func(std::function<uint64_t(_Kty&&)> &&hash_key_func)
{
    if(len[!backup] || len[backup]) return false;
    else hash_func = std::move(hash_key_func);
    return true;
}
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::value_copy(net_map &src) { __value_copy(src); }
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::value_move(net_map &&src)
{
    hash_factor = src.hash_factor;
    rehash_load = src.rehash_load;
    base_load = src.base_load;
    len[0] = src.len[0]; len[1] = src.len[1];
    backup = src.backup;
    hash_func = std::move(src.hash_func);
    kv_data[0] = std::move(src.kv_data[0]); kv_data[1] = std::move(src.kv_data[1]);
    kv_occupy[0] = std::move(src.kv_occupy[0]); kv_occupy[1] = std::move(src.kv_occupy[1]);
}
template<typename _Kty, typename _Ty> inline uint64_t net_map<_Kty, _Ty>::size() { return len[!backup]; }
template<typename _Kty, typename _Ty> inline uint64_t net_map<_Kty, _Ty>::mem_size(){ return kv_data[!backup].mem_size(); }
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::rehash(bool enforce)
{
    if(len[!backup]>kv_data[!backup].size()*hash_factor || enforce)
    {
        _rehash();
        rehash_load += rehash_load;
    }
}
template<typename _Kty, typename _Ty> inline long long net_map<_Kty, _Ty>::find_idx(_Kty &&key)
{
    auto tgt_idx = _find_idx(std::move(key), !backup);
    if(kv_occupy[!backup][tgt_idx]) return tgt_idx;
    else return -1;
}
template<typename _Kty, typename _Ty> inline long long net_map<_Kty, _Ty>::find_idx(_Kty &key) { return find_idx(std::move(key)); }
template<typename _Kty, typename _Ty> inline net_sequence<_Kty> net_map<_Kty, _Ty>::find_key(_Ty &&value)
{
    net_sequence<_Kty> ans;
    for(auto i=0; i<kv_occupy[!backup].size(); ++i) if(kv_occupy[!backup][i] && kv_data[!backup][i].value == value) ans.emplace_back(kv_data[!backup][i].key);
    return ans;
}
template<typename _Kty, typename _Ty> inline net_sequence<_Kty> net_map<_Kty, _Ty>::find_key(_Ty &value) { return find_key(std::move(value)); }
template<typename _Kty, typename _Ty> inline bool net_map<_Kty, _Ty>::insert(_Kty &&key, _Ty &&value)
{
    rehash();
    auto tgt_idx = _find_idx(std::move(key), !backup);
    if(kv_occupy[!backup][tgt_idx]) return false;
    else
    {
        kv_occupy[!backup][tgt_idx] = true;
        kv_data[!backup][tgt_idx].key = std::move(key);
        kv_data[!backup][tgt_idx].value = std::move(value);
        ++ len[!backup];
        return true;
    }
}
template<typename _Kty, typename _Ty> inline bool net_map<_Kty, _Ty>::insert(_Kty &key, _Ty &value)
{
    auto _key = key;
    auto _value = value;
    return insert(std::move(_key), std::move(_value));
}
template<typename _Kty, typename _Ty> inline typename net_map<_Kty, _Ty>::net_kv net_map<_Kty, _Ty>::erase(_Kty &&key)
{
    rehash();
    net_kv temp;
    auto tgt_idx = _find_idx(std::move(key), !backup);
    if(kv_occupy[!backup][tgt_idx])
    {
        temp.key = std::move(kv_data[!backup][tgt_idx].key);
        temp.value = std::move(kv_data[!backup][tgt_idx].value);
        kv_occupy[!backup][tgt_idx] = false;
        -- len[!backup];
    }
    return temp;
}
template<typename _Kty, typename _Ty> inline typename net_map<_Kty, _Ty>::net_kv net_map<_Kty, _Ty>::erase(_Kty &key) { return erase(std::move(key)); }
template<typename _Kty, typename _Ty> inline typename net_map<_Kty, _Ty>::net_kv net_map<_Kty, _Ty>::index(uint64_t idx)
{
    rehash();
    if(idx < kv_data[!backup].size()) return kv_data[!backup][idx];
    else return net_kv();
}
template<typename _Kty, typename _Ty> inline bool net_map<_Kty, _Ty>::occupy(uint64_t idx) { return kv_occupy[!backup][idx]; }
template<typename _Kty, typename _Ty> inline _Ty &net_map<_Kty, _Ty>::operator[](_Kty &&key)
{
    rehash();
    auto tgt_idx = _find_idx(std::move(key), !backup);
    assert(kv_occupy[!backup][tgt_idx]);
    return kv_data[!backup][tgt_idx].value;
}
template<typename _Kty, typename _Ty> inline _Ty &net_map<_Kty, _Ty>::operator[](_Kty &key) { return operator[](std::move(key)); }
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::operator=(net_map &src) { value_copy(src); }
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::operator=(const net_map &src) { __value_copy(src); }
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::operator=(net_map &&src) { value_move(std::move(src)); }
template<typename _Kty, typename _Ty> inline bool net_map<_Kty, _Ty>::operator==(net_map &val)
{
    if(val.len[!backup]==len[!backup])
    {
        for(auto i=0; i<kv_occupy[!backup].size(); ++i) if((val.kv_occupy[!backup][i]!=kv_occupy[!backup][i]) || (kv_occupy[!backup][i]&&kv_data[!backup][i]!=val.kv_data[!backup][i])) return false;
        return true;
    }
    else return false;
}
template<typename _Kty, typename _Ty> inline bool net_map<_Kty, _Ty>::operator!=(net_map &val) { return !(*this == val); }
template<typename _Kty, typename _Ty> inline void net_map<_Kty, _Ty>::reset()
{
    hash_factor = 0;
    backup = false;
    len[0] = 0; len[1] = 0;
    rehash_load = 0; base_load = 0;
    kv_data[0].reset(); kv_data[1].reset();
    kv_occupy[0].reset(); kv_occupy[1].reset();
}
template<typename _Kty, typename _Ty> inline net_map<_Kty, _Ty>::~net_map() { reset(); }

/* Time log */
template<typename _Kty> inline clock_timer<_Kty>::clock_timer(std::function<uint64_t(_Kty&&)> _hash_func) : clock_log(_hash_func) {}
template<typename _Kty> inline clock_timer<_Kty>::clock_timer(clock_timer &src) { *this = src; }
template<typename _Kty> inline clock_timer<_Kty>::clock_timer(const clock_timer &src) { *this = src; }
template<typename _Kty> inline clock_timer<_Kty>::clock_timer(clock_timer &&src) { *this = std::move(src); }
template<typename _Kty> inline bool clock_timer<_Kty>::clock_begin(_Kty &&clock_log_id)
{
    auto begin_pt = clock();
    _dur dur_temp;
    dur_temp.begin = begin_pt;
    auto id_temp = clock_log.find_idx(clock_log_id);
    if(id_temp < 0) return clock_log.insert(std::move(clock_log_id), std::move(dur_temp));
    else
    {
        clock_log[clock_log_id] = dur_temp;
        return true;
    }
}
template<typename _Kty> inline bool clock_timer<_Kty>::clock_begin(_Kty &clock_log_id) { return clock_begin(std::move(clock_log_id)); }
template<typename _Kty> inline bool clock_timer<_Kty>::clock_end(_Kty &&clock_log_id)
{
    auto end_point = clock();
    auto log_idx = clock_log.find_idx(clock_log_id);
    if(log_idx < 0) return false;
    else
    {
        clock_log[clock_log_id].end = end_point;
        clock_log[clock_log_id].dur = (clock_log.index(log_idx).value.end - clock_log.index(log_idx).value.begin) * 1000 / CLOCKS_PER_SEC;
        return true;
    }
}
template<typename _Kty> inline bool clock_timer<_Kty>::clock_end(_Kty &clock_log_id) { return clock_end(std::move(clock_log_id)); }
template<typename _Kty> inline long clock_timer<_Kty>::duration(_Kty &&clock_log_id)
{
    auto log_idx = clock_log.find_idx(clock_log_id);
    if(log_idx < 0) return -1;
    else return clock_log[clock_log_id].dur;
}
template<typename _Kty> inline long clock_timer<_Kty>::duration(_Kty &clock_log_id) { return duration(std::move(clock_log_id)); }
template<typename _Kty> inline void clock_timer<_Kty>::operator=(const clock_timer &src) { clock_log = src.clock_log; }
template<typename _Kty> inline void clock_timer<_Kty>::operator=(clock_timer &src) { clock_log = src.clock_log; }
template<typename _Kty> inline void clock_timer<_Kty>::operator=(clock_timer &&src) { clock_log = std::move(src.clock_log); }
template<typename _Kty> inline void clock_timer<_Kty>::reset() { clock_log.reset(); }
template<typename _Kty> inline clock_timer<_Kty>::~clock_timer() { reset(); }

/* Net memory */
/* Pointer */
template<typename _Ty> inline net_memory<_Ty>::net_ptr::net_ptr(net_memory *_ptr, uint64_t _begin_addr, uint64_t _len) : ptr(_ptr), len(_len), begin_addr(_begin_addr) {}
template<typename _Ty> inline net_memory<_Ty>::net_ptr::net_ptr(const net_ptr &src) { *this = src; }
template<typename _Ty> inline net_memory<_Ty>::net_ptr::net_ptr(net_ptr &src) { *this = src; }
template<typename _Ty> inline net_memory<_Ty>::net_ptr::net_ptr(net_ptr &&src) { *this = std::move(src); }
template<typename _Ty> inline uint64_t net_memory<_Ty>::net_ptr::size() { return len; }
template<typename _Ty> inline void net_memory<_Ty>::net_ptr::operator=(net_ptr &src) { ptr = src.ptr; len = src.len; begin_addr = src.begin_addr; }
template<typename _Ty> inline void net_memory<_Ty>::net_ptr::operator=(const net_ptr &src) { ptr = src.ptr; len = src.len; begin_addr = src.begin_addr; }
template<typename _Ty> inline void net_memory<_Ty>::net_ptr::operator=(net_ptr &&src)
{
    len = src.len; begin_addr = src.begin_addr; ptr = src.ptr;
    src.ptr = nullptr; src.len = 0; src.begin_addr = 0;
}
template<typename _Ty> inline _Ty &net_memory<_Ty>::net_ptr::operator[](uint64_t addr)
{
    assert(addr < len);
    return ptr->mem_val[begin_addr+addr];
}
template<typename _Ty> inline net_memory<_Ty>::net_ptr::~net_ptr() { ptr = nullptr; len = 0; begin_addr = 0; }
/* Memory */
template<typename _Ty> inline net_memory<_Ty>::net_memory(uint64_t _alloc_size) : mem_val(_alloc_size, _alloc_size) {}
template<typename _Ty> inline net_memory<_Ty>::net_memory(const net_memory &src) { *this = src; }
template<typename _Ty> inline net_memory<_Ty>::net_memory(net_memory &src) { *this = src; }
template<typename _Ty> inline net_memory<_Ty>::net_memory(net_memory &&src) { *this = std::move(src); }
template<typename _Ty> inline uint64_t net_memory<_Ty>::size() { return len; }
template<typename _Ty> inline uint64_t net_memory<_Ty>::mem_size() { return (idx_seq.length ? (idx_seq[rear_id].end_addr + 1) : 0); }
template<typename _Ty> inline uint64_t net_memory<_Ty>::max_mem_size() { return mem_val.mem_size(); }
template<typename _Ty> inline void net_memory<_Ty>::print_block_info(int id, bool detail)
{
    if(id>=0 && id<idx_seq.size())
    {
        std::cout << "[ID " << id << "][Size ";
        if(idx_seq[id].occupy) std::cout << idx_seq[id].len;
        else std::cout << 0;
        std::cout << '/' << mem_block_mem_len(std::move(idx_seq[id])) << ']';
        if(detail) std::cout << "[Address " << idx_seq[id].begin_addr << '-' << idx_seq[id].end_addr << "][Previous " << idx_seq[id].prev_id << " Next " << idx_seq[id].next_id << ']';
    }
    else std::cerr << "[Illegal ID]" << std::endl;
}
template<typename _Ty> inline void net_memory<_Ty>::print_mem_info(bool detail)
{
    int tool_idx = head_id;
    std::cout << "[Memory Structure]" << std::endl;
    if(idx_seq.length) while(tool_idx>=0)
    {
        if(tool_idx != head_id)
        {
            std::cout << " -> ";
            if(detail) std::cout << std::endl;
        }
        print_block_info(tool_idx, detail);
        tool_idx = idx_seq[tool_idx].next_id;
    }
    else std::cout << "[Blank]";
    std::cout << std::endl;
    std::cout << "[Memory Size]" << std::endl;
    std::cout << "[Block " << size() << '/' << mem_size() << "][Max Size " << max_mem_size() << ']' << std::endl;
}
template<typename _Ty> inline int net_memory<_Ty>::alloc_mem(uint64_t _alloc_size, _Ty *&&src)
{
    _alloc_size = num_unsign(_alloc_size);
    auto tar_idx = -1;
    if(_alloc_size)
    {
        for(auto i=0; i<idx_seq.size(); ++i) if(!idx_seq[i].occupy && _alloc_size<=mem_block_mem_len(std::move(idx_seq[i])))
        {
            tar_idx = i;
            idx_seq[i].occupy = true;
            idx_seq[i].len = _alloc_size;
            break;
        }
        if(tar_idx < 0)
        {
            mem_block blk_temp;
            blk_temp.occupy = true;
            blk_temp.begin_addr = 0;
            if(idx_seq.size()) blk_temp.begin_addr = idx_seq[rear_id].end_addr + 1;
            blk_temp.end_addr = blk_temp.begin_addr + _alloc_size - 1;
            blk_temp.len = _alloc_size;
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
                mem_val.extend(realloc_size, true);
            }
            idx_seq.emplace_back(std::move(blk_temp));
        }
        if(src)
        {
            for(auto i=0; i<_alloc_size; ++i) mem_val[idx_seq[tar_idx].begin_addr+i] = std::move(src[i]);
            MEM_RECYCLE(src);
        }
        len += _alloc_size;
    }
    return tar_idx;
}
template<typename _Ty> inline bool net_memory<_Ty>::free_mem(int &id, bool remain)
{
    if(id>=0 && id<idx_seq.size()) if(idx_seq[id].occupy)
    {
        len -= idx_seq[id].len;
        if(len || remain)
        {
            idx_seq[id].occupy = false;
            id = -1; return true;
        }
        else
        {
            *this = net_memory();
            return true;
        }
    }
    return false;
}
template<typename _Ty> inline void net_memory<_Ty>::shrink(bool shrink_blk)
{
    if(len)
    {
        int next_id_temp = -1, tool_idx = head_id; len = 0;
        while(tool_idx >= 0)
        {
            next_id_temp = idx_seq[tool_idx].next_id;
            if(idx_seq[tool_idx].occupy)
            {
                len += idx_seq[tool_idx].len;
                auto begin_addr_temp = 0;
                if(idx_seq[tool_idx].prev_id >= 0) begin_addr_temp = idx_seq[idx_seq[tool_idx].prev_id].end_addr + 1;
                if(begin_addr_temp != idx_seq[tool_idx].begin_addr)
                {
                    for(auto i=0; i<idx_seq[tool_idx].len; ++i) mem_val[begin_addr_temp+i] = std::move(mem_val[idx_seq[tool_idx].begin_addr+i]);
                    idx_seq[tool_idx].end_addr -= idx_seq[tool_idx].begin_addr - begin_addr_temp;
                    idx_seq[tool_idx].begin_addr = begin_addr_temp;
                }
                if(shrink_blk)
                {
                    auto seg_len = idx_seq[tool_idx].end_addr - idx_seq[tool_idx].begin_addr + 1 - idx_seq[tool_idx].len;
                    if(seg_len) idx_seq[tool_idx].end_addr -= seg_len;
                }
                if(!idx_seq[tool_idx].begin_addr) head_id = tool_idx;
            }
            else if(idx_seq[tool_idx].len)
            {
                idx_seq[tool_idx].len = 0;
                if(idx_seq[tool_idx].prev_id >= 0) idx_seq[idx_seq[tool_idx].prev_id].next_id = idx_seq[tool_idx].next_id;
                if(idx_seq[tool_idx].next_id >= 0) idx_seq[idx_seq[tool_idx].next_id].prev_id = idx_seq[tool_idx].prev_id;
                idx_seq[tool_idx].prev_id = rear_id;
                idx_seq[rear_id].next_id = tool_idx;
                rear_id = tool_idx;
                idx_seq[tool_idx].next_id = -1;
            }
            else
            {
                auto blk_size = mem_block_mem_len(std::move(idx_seq[tool_idx]));
                if(idx_seq[tool_idx].prev_id < 0) idx_seq[tool_idx].begin_addr = 0;
                else idx_seq[tool_idx].begin_addr = idx_seq[idx_seq[tool_idx].prev_id].end_addr + 1;
                idx_seq[tool_idx].end_addr = idx_seq[tool_idx].begin_addr + blk_size - 1;
            }
            tool_idx = next_id_temp;
        }
    }
    else *this = net_memory();
}
template<typename _Ty> inline bool net_memory<_Ty>::exist(int id) { return (id>=0 && id<idx_seq.size() && idx_seq[id].occupy); }
template<typename _Ty> inline void net_memory<_Ty>::operator=(net_memory &src) { value_assign(src); mem_val = src.mem_val; idx_seq = src.idx_seq; }
template<typename _Ty> inline void net_memory<_Ty>::operator=(const net_memory &src) { value_assign(src); mem_val = src.mem_val; idx_seq = src.idx_seq; }
template<typename _Ty> inline void net_memory<_Ty>::operator=(net_memory &&src) { value_assign(src); mem_val = std::move(src.mem_val); idx_seq = std::move(src.idx_seq); }
template<typename _Ty> inline typename net_memory<_Ty>::net_ptr net_memory<_Ty>::operator[](int id)
{
    if(exist(id)) return net_ptr(this, idx_seq[id].begin_addr, idx_seq[id].len);
    else return net_ptr(nullptr, 0, 0);
}
template<typename _Ty> inline void net_memory<_Ty>::reset() { head_id = 0; rear_id = 0; len = 0; mem_val.reset(); idx_seq.reset(); }
template<typename _Ty> inline net_memory<_Ty>::~net_memory() { reset(); }
template<typename _Ty> inline void net_memory<_Ty>::value_assign(net_memory &src) { head_id = src.head_id; rear_id = src.rear_id; len = src.len; }
template<typename _Ty> inline uint64_t net_memory<_Ty>::mem_block_mem_len(mem_block &&src) { return num_cnt(src.begin_addr, src.end_addr); }

/* Basic algorithm */
uint64_t integer_digit_cnt(uint64_t src)
{
    auto cnt = 0;
    while(src /= 10) ++ cnt;
    return cnt;
}
/* Integer operation */ 
net_set<uint64_t> integer_radix_sort(net_set<uint64_t> &src)
{
    net_set<net_sequence<uint64_t>> tool(10);
    auto ans = src;
    uint64_t base = 1;
    while(true)
    {
        for(auto i=0; i<src.size(); ++i) tool[(ans[i]/base)%10].push_back(ans[i]);
        if(tool[0].size() == src.size()) break;
        auto cnt = 0;
        for(auto i=0; i<tool.size(); ++i)
        {
            for(auto j=0; j<tool[i].size(); ++j) ans[cnt++] = tool[i][j];
            tool[i].clear();
        }
        base *= 10;
    }
    return ans;
}
int integer_bit_reverse(int src, int bit_cnt = 32)
{
    auto ans = 0;
    while(bit_cnt)
    {
        ans <<= 1;
        ans += src & 1;
        src >>= 1;
        -- bit_cnt;
    }
    return ans;
}
int integer_bit_cnt(int src)
{
    auto ans = 0;
    while(src)
    {
        ++ ans;
        src >>= 1;
    }
    return ans;
}
/* Prime */
net_set<uint64_t> integer_primes(uint64_t upper)
{
    net_set<uint64_t> prm_set;
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
net_set<uint64_t> integer_primes_fact(uint64_t val)
{
    auto prm_set = integer_primes(val);
    auto itr_cnt = 0;
    net_set<uint64_t> prm_fact;
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
/* Fraction */
uint64_t integer_greatest_common_divisor(uint64_t l_val, uint64_t r_val)
{
    auto itr_prm = integer_primes_fact(l_val).unit_intersect(integer_primes_fact(r_val));
    auto res = 1;
    for(auto i=0; i<itr_prm.size(); ++i) res *= itr_prm[i];
    return res;
}
uint64_t integer_least_common_multiple(uint64_t l_val, uint64_t r_val)
{
    auto itr_prm = integer_primes_fact(l_val).unit_union(integer_primes_fact(r_val));
    auto res = 1;
    for(auto i=0; i<itr_prm.size(); ++i) res *= itr_prm[i];
    return res;
}
/* Polynormial calculation */
net_sequence<uint8_t> integer_polynomial_add(bool &ans_sgn, net_sequence<uint8_t> &coe_a, net_sequence<uint8_t> &coe_b, bool sgn_a, bool sgn_b)
{
    if(coe_b.length == coe_a.length)
    {
        auto sgn_same = sgn_a == sgn_b;
        auto ans_len = sgn_same ? coe_a.length + 1 : coe_a.length;
        net_sequence<uint8_t> ans(ans_len);
        auto conj_a = sgn_a ? sgn_same ? 1 : -1 : 1,
            conj_b = sgn_b ? sgn_same ? 1 : -1 : 1,
            carry = 0;
        ans_sgn = false;
        for(auto i=0; i<coe_a.length; ++i)
        {
            carry += conj_a * coe_a[i] + conj_b * coe_b[i];
            auto curr_bit = carry % 10;
            carry /= 10;
            if(curr_bit < 0)
            {
                curr_bit = 10 + curr_bit;
                -- carry;
            }
            ans[i] = curr_bit;
        }
        if(carry < 0)
        {
            for(auto i=0; i<coe_a.length; ++i)
            {
                ans[i] = 10 - ans[i];
                if(i) -- ans[i];
            }
            if(ans.length > coe_a.length) ans[coe_a.length] = (-1) * carry - 1;
            ans_sgn = true;
        }
        else
        {
            if(ans.length > coe_a.length)  ans[coe_a.length] = carry;
            ans_sgn = sgn_a && sgn_same;
        }
        return ans;
    }
    else return net_sequence<uint8_t>();
}
net_sequence<uint64_t> integer_polynomial_mult(net_sequence<uint8_t> &coe_a, net_sequence<uint8_t> &coe_b)
{
    net_sequence<uint64_t> ans(coe_a.length + coe_b.length - 1);
    ans.memory_set(0);
    for(auto i=0; i<coe_a.length; ++i) for(auto j=0; j<coe_b.length; ++j)
    {
        uint64_t curr_a = coe_a[i], curr_b = coe_b[j];
        ans[i+j] += curr_a * curr_b;
    }
    return ans;
}
/* FFT */
bool integer_fft(std::complex<long double> *&src, uint64_t len, bool inverse)
{
    if(len == 1) return true;
    else if(len%2) return false;
    // Sub source
    MEM_INIT(std::complex<long double>, src_sub_left, len>>1);
    MEM_INIT(std::complex<long double>, src_sub_right, len>>1);
    for(auto i=0; i<len; i+=2)
    {
        src_sub_left[i>>1] = src[i];
        src_sub_right[i>>1] = src[i+1];
    }
    auto flag = integer_fft(src_sub_left, len>>1, inverse) && integer_fft(src_sub_right, len>>1, inverse);
    if(flag)
    {
        // dft
        // idft
        auto conj = 1;
        if(inverse) conj = -1;
        // omega
        std::complex<long double> omega_unit(std::cos(2*pi/len), conj*std::sin(2*pi/len)), omega(1, 0);
        len >>= 1;
        for(auto i=0; i<len; ++i)
        {
            src[i] = src_sub_left[i] + omega * src_sub_right[i];
            src[i+len] = src_sub_left[i] - omega * src_sub_right[i];
            omega *= omega_unit;
        }
    }
    MEM_RECYCLE(src_sub_left);
    MEM_RECYCLE(src_sub_right);
    return flag;
}
bool integer_dft(std::complex<long double> *&c, std::complex<long double> *&a, std::complex<long double> *&b, uint64_t len)
{
    if(c && a && b && len && integer_fft(a, len) && integer_fft(b, len))
    {
        for(auto i=0; i<len; ++i)
        {
            if(c+i) c[i] = a[i] * b[i];
            else return false;
        }
        return true;
    }
    else return false;
}
bool integer_idft(uint64_t *&d, std::complex<long double> *&c, uint64_t len)
{
    if(c && d && len && integer_fft(c, len, true))
    {
        for(auto i=0; i<len; ++i)
        {
            if(d+i) d[i] = (uint64_t)(c[i].real() / len + 0.5);
            else return false;
        }
        return true;
    }
    else return false;
}
/* NTT */
uint64_t integer_euler_pow(uint64_t base, uint64_t times)
{
    uint64_t ans = 1;
    while(times)
    {
        if(times & 1) ans = ans * base % EULER_MOD;
        base = base * base % EULER_MOD;
        times >>= 1;
    }
    return ans;
}
bool integer_ntt(uint64_t *&src, uint64_t len, bool inverse)
{
    if(!num_pow_pad_cnt(len, 2) && len>=4 && src)
    {
        auto idx_bit_cnt = integer_bit_cnt(len-1);
        MEM_INIT(uint64_t, temp, len);
        for(auto i=0; i<len; ++i) temp[i] = src[integer_bit_reverse(i, idx_bit_cnt)];
        ptr_move(src, std::move(temp));
        for(auto i=1; i<len; i<<=1)
        {
            auto g_root = integer_euler_pow(EULER_MOD_G, (EULER_MOD-1)/(i<<1));
            if(inverse) g_root = integer_euler_pow(g_root, EULER_MOD-2);
            for(auto j=0; j<len; j+=(i<<1))
            {
                auto g_unit = 1ull;
                for(auto k=0; k<i; ++k)
                {
                    auto g_unit_front = src[j+k], g_unit_rear = g_unit * src[i+j+k] % EULER_MOD;
                    src[j+k] = (g_unit_front + g_unit_rear) % EULER_MOD;
                    src[i+j+k] = (g_unit_front - g_unit_rear + EULER_MOD) % EULER_MOD;
                    g_unit = g_unit * g_root % EULER_MOD;
                }
            }
        }
        if(inverse)
        {
            auto inverser = integer_euler_pow(len, EULER_MOD-2);
            for(auto i=0; i<len; ++i) src[i] = src[i] * inverser % EULER_MOD;
        }
        return true;
    }
    else return false;
}
bool integer_fnt(uint64_t *&coe_c, uint64_t *&coe_a, uint64_t *&coe_b, uint64_t len)
{
    if(coe_c && coe_a && coe_b && len && integer_ntt(coe_a, len) && integer_ntt(coe_b, len))
    {
        for(auto i=0; i<len; ++i) coe_c[i] = coe_a[i] * coe_b[i] % EULER_MOD;
        return true;
    }
    else return false;
}
bool integer_ifnt(uint64_t *&coe_c, uint64_t len) { return (coe_c && len && integer_ntt(coe_c, len, true)); }

/* Accuracy */
bool acc_valid(double acc) { while(acc < 1) acc *= 10; return acc == 1; }
double acc_norm(double acc)
{
    if(acc<1 && acc>0)
    {
        auto _acc = acc;
        auto bit_cnt = 1;
        while(_acc > 0)
        {
            bit_cnt *= 10;
            _acc *= 10;
            _acc -= (int)_acc;
        }
        acc = 1.0 / bit_cnt;
    }
    return acc;
}
double acc_round(double _val, double acc)
{
    acc = acc_norm(acc);
    auto _itgr = (int)_val;
    auto _dcml = (_val - _itgr) / acc;
    _dcml -= (int)_dcml;
    _dcml *= acc;
    _val -= _dcml;
    if(_dcml < 0.5*acc) return _val;
    else return _val + acc;
}

/* Pseudo random number */
double rand_base(long long seed, double acc)
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
double rand_num(double boundary_first, double boundary_second, bool to_sleep, double acc)
{
    if(boundary_first == boundary_second) return boundary_second;
    else if(acc_valid(acc))
    {
        if(to_sleep) _sleep(1);
        long long raw_pt = std::chrono::system_clock::now().time_since_epoch().count();
        auto seed_pt = rand_base(raw_pt, acc);
        return (boundary_first<boundary_second) ? (seed_pt*(boundary_second-boundary_first)+boundary_first) : ((seed_pt*(boundary_first-boundary_second))+boundary_second);
    }
    else return 0;
}
decimal rand_num(decimal &boundary_first, decimal &boundary_second, bool to_sleep, decimal &acc) { return rand_num(boundary_first.to_float(), boundary_second.to_float(), to_sleep, acc.to_float()); }
net_set<uint64_t> rand_idx(uint64_t seq_size, uint64_t amt)
{
    if(seq_size && amt)
    {
        net_set<uint64_t> idx_seq(seq_size);
        for(auto i=0; i<seq_size; ++i) idx_seq[i] = i;
        idx_seq.shuffle();
        net_set<uint64_t> idx_seq_shuffle(seq_size);
        for(auto i=0; i<seq_size; ++i) idx_seq_shuffle[i] = i;
        idx_seq_shuffle.shuffle();
        net_set<uint64_t> ans(amt);
        for(auto i=0; i<ans.size(); ++i) ans[i] = idx_seq[idx_seq_shuffle[i]];
        ans.sort();
        idx_seq.reset();
        idx_seq_shuffle.reset();
        return ans;
    }
    else return net_set<uint64_t>();
}

/* Character setting */
std::string charset_input_paragraph(uint64_t buffer_length)
{
    std::fflush(stdin);
    std::cout << "Please end submit with double enter." << std::endl;
    std::cout << std::endl; 
    char *cmtr = nullptr;
    MEM_ALLOC(char, cmtr, buffer_length);
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
    MEM_RECYCLE(cmtr);
    return init_str;
}
bool charset_extract_number(std::string &num_str, double *&num_arr, uint64_t &arr_len)
{
    int i = 0, cnt = 0, len = num_str.length(), amt = 0;
    // get the amount of the number from string
    while(i < len){
        if(num_str[i]=='.' || num_str[i]=='-' || num_str[i]=='+' || (num_str[i]>='0' && num_str[i]<='9'))
            if(i == 0 || num_str[i-1]==' '|| num_str[i-1]=='\t' || num_str[i-1]=='\n') cnt ++;
        i ++;
    }
    if(!cnt) return false;
    if(num_arr) { MEM_RECYCLE(num_arr); }
    MEM_ALLOC(double, num_arr, cnt);
    arr_len = cnt;
    amt = cnt;
    cnt = 0;
    for(int j=0; j<amt; ++j){
        num_arr[j] = 0.0;
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
        } num_arr[cnt ++] = temp;
        i ++;
    }
    return true;
}
std::string charset_exchange(std::wstring &str_src)
{
    auto nLen = str_src.length();
    char *psBuf = nullptr;
    MEM_ALLOC(char, psBuf, nLen+1);
    wcstombs(psBuf, str_src.c_str(), nLen);
    psBuf[nLen] = '\0';
    std::string strDest(psBuf);
    MEM_RECYCLE(psBuf);
    if(strDest.length() == nLen) return strDest;
    else return "";
}
std::wstring charset_exchange(std::string &str_src)
{
    // setlocale(LC_ALL, "zh_CN.UTF-8");
    auto nLen = str_src.length();
    wchar_t *pwsBuf = nullptr;
    MEM_ALLOC(wchar_t, pwsBuf, nLen+1);
    mbstowcs(pwsBuf, str_src.c_str(), nLen);
    pwsBuf[nLen] = L'\0';
    std::wstring wstrDest(pwsBuf);
    MEM_RECYCLE(pwsBuf);
    if(wstrDest.length() == nLen) return wstrDest;
    else return L"";
}

/* Decimal */
inline net_sequence<uint8_t> decimal::dec_coe(bool it_or_ft, bool from_low)
{
    net_sequence<uint8_t> ans;
    for(auto i=0ull; i<val[it_or_ft].length; ++i)
    {
        auto curr_seg = val[it_or_ft][i];
        for(auto j=0; j<DECIMAL_DIG_MAX; ++j)
        {
            if(curr_seg)
            {
                ans.emplace_back(curr_seg % 10);
                curr_seg /= 10;
            }
            else if(i+1 == val[it_or_ft].length) break;
            else ans.emplace_back(0);
        }
    }
    if((it_or_ft==it&&!from_low) || (it_or_ft==ft&&from_low)) ans.reverse();
    return ans;
}
inline void decimal::dec_coe(bool it_or_ft, net_sequence<uint8_t> &src, bool from_low)
{
    val[it_or_ft].reset();
    // Calibrate
    if((it_or_ft==ft&&from_low) || (it_or_ft==it&&!from_low)) while(src.length && !src[0]) src.erase(0);
    if((it_or_ft==ft&&!from_low) || (it_or_ft==it&&from_low)) while(src.length && !src[src.length-1]) src.erase(src.length-1);
    if(!src.length) return;
    auto seg_times = DECIMAL_SEG_MAX;
    for(auto i=it_or_ft ? (from_low?src.length:0) : (from_low?0:src.length);
        it_or_ft ? (from_low?i>0:i<src.length) : (from_low?i<src.length:i>0);
        it_or_ft ? (from_low?--i:++i) : (from_low?++i:--i))
    {
        auto curr_idx = it_or_ft ? (from_low?(i-1):i) : (from_low?i:(i-1));
        if(seg_times == DECIMAL_SEG_MAX)
        {
            val[it_or_ft].emplace_back(0);
            seg_times = 1;
        }
        val[it_or_ft][val[it_or_ft].length-1] += seg_times * src[curr_idx];
        seg_times *= 10;
    }
}
inline void decimal::__value_copy(const decimal &src)
{
    sgn = src.sgn;
    val[it] = src.val[it];
    val[ft] = src.val[ft];
    val[ft].shrink(); val[it].shrink();
}
inline void decimal::value_copy(decimal &src) { __value_copy(src); }
inline void decimal::value_move(decimal &&src)
{
    sgn = src.sgn;
    val[it] = std::move(src.val[it]);
    val[ft] = std::move(src.val[ft]);
    val[ft].shrink(); val[it].shrink();
    src.reset();
}
inline decimal::decimal(const decimal &src) { __value_copy(src); }
inline decimal::decimal(decimal &src) { value_copy(src); }
inline decimal::decimal(decimal &&src) { value_move(std::move(src)); }
inline decimal::decimal(long double src) : sgn(src<0)
{
    if(src)
    {
        if(sgn) src *= (-1);
        auto it_side = (uint64_t)src;
        auto ft_side = src - it_side;
        auto ft_dig_cnt = calculate_digit;
        if(it_side)
        {
            net_sequence<uint8_t> it_temp;
            while(it_side)
            {
                it_temp.emplace_back(it_side%10);
                it_side /= 10;
            }
            dec_coe(it, it_temp);
        }
        if(ft_dig_cnt && ft_side)
        {
            net_sequence<uint8_t> ft_temp;
            while(ft_dig_cnt)
            {
                -- ft_dig_cnt;
                ft_side *= 10;
                ft_temp.emplace_back((uint8_t)ft_side);
                ft_side -= (uint8_t)ft_side;
                if(!ft_dig_cnt)
                {
                   ft_side *= 10;
                   if(ft_side >= 5) ++ ft_temp[ft_temp.length-1];
                }
            }
            dec_coe(ft, ft_temp, false);
        }
        val[ft].shrink(); val[it].shrink(); 
    }
}
inline decimal::decimal(std::string src)
{
    if(src.length())
    {
        // Calibrate
        while(src[0] == '0') src.erase(0, 1);
        while(src[src.length()-1] == '0') src.erase(src.length()-1, src.length());
        // Symbol check
        auto dot_idx = src.length();
        auto dot_flag = false;
        if(src[src.length()-1] == '.')
        {
            src.erase(src.length()-1, src.length());
            dot_flag = true;
        }
        if(src[0] == '.')
        {
            dot_idx = 0;
            src[0] = '0';
            dot_flag = true;
        }
        if(src[0]<'0' || src[0]>'9')
        {
            if(src[0] == '-') sgn = true;
            else if(src[0] == '+') sgn = false;
            else return;
            src.erase(0, 1);
        }
        for(auto i=0ull; i<src.length(); ++i)
        {
            if((src[i]<'0' || src[i]>'9') && dot_flag) return;
            if(src[i] == '.')
            {
                if(dot_flag) return;
                dot_flag = true;
                dot_idx = i;
            }
        }
        // integer
        if(dot_idx > 0)
        {
            net_sequence<uint8_t> it_temp;
            for(auto i=dot_idx; i>0; --i) it_temp.emplace_back(src[i-1]-'0');
            dec_coe(it, it_temp);
        }
        // float
        if(dot_idx < src.length())
        {
            net_sequence<uint8_t> ft_temp;
            for(auto i=dot_idx+1; i<src.length(); ++i) ft_temp.emplace_back(src[i]-'0');
            dec_coe(ft, ft_temp, false);
        }
        val[ft].shrink(); val[it].shrink(); 
    }
}
inline uint64_t decimal::it_len() { return dec_coe(it).length; }
inline uint64_t decimal::ft_len() { return dec_coe(ft).length; }
inline bool decimal::zero()
{
    if(val[it].length || val[ft].length) return false;
    else return true;
}
inline std::string decimal::to_string()
{
    auto it_coe = dec_coe(it, false), ft_coe = dec_coe(ft, false);
    std::string ans = "";
    if(sgn) ans.push_back('-');
    if(it_coe.length) for(auto i=0ull; i<it_coe.length; ++i) ans.push_back(it_coe[i]+'0');
    else ans.push_back('0');
    if(ft_coe.length) 
    {
        ans.push_back('.');
        for(auto i=0ull; i<ft_coe.length; ++i) ans.push_back(ft_coe[i]+'0');
    }
    return ans;
}
inline long double decimal::to_float()
{
    long double ans = to_integer(), ans_ft = 0;
    auto ft_coe = dec_coe(ft);
    for(auto i=0ull; i<ft_coe.length; ++i)
    {
        ans_ft += ft_coe[i];
        ans_ft /= 10;
    }
    ans += ans_ft;
    return ans;
}
inline int64_t decimal::to_integer()
{
    auto ans = 0ull;
    auto it_coe = dec_coe(it, false);
    for(auto i=0ull; i<it_coe.length; ++i)
    {
        ans *= 10;
        ans += it_coe[i];
    }
    return ans;
}
inline void decimal::_acc(uint64_t con_acc)
{
    if(!con_acc) val[ft].reset();
    else
    {
        auto coe_ft = dec_coe(ft, false);
        if(coe_ft.length > con_acc)
        {
            bool carry = false;
            if(coe_ft[con_acc] >= 5) carry = true;
            coe_ft.cut(con_acc);
            -- con_acc;
            while(carry)
            {
                ++ coe_ft[con_acc];
                if(coe_ft[con_acc] == 10)
                {
                    carry = true;
                    coe_ft.erase(con_acc);
                }
                else carry = false;
                if(con_acc) -- con_acc;
                else break;
            }
            dec_coe(ft, coe_ft, false);
            if(carry)
            {
                auto coe_it = dec_coe(it);
                for(auto i=0ull; i<coe_it.length&&carry; ++i)
                {
                    ++ coe_it[i];
                    if(coe_it[i] == 10)
                    {
                        carry = true;
                        coe_it[i] = 0;
                    }
                    else carry = false;
                }
                if(carry) coe_it.emplace_back(1);
                dec_coe(it, coe_it);
            }
        }
    }
}
inline uint64_t decimal::_acc() { return ft_len(); }
inline decimal decimal::abs()
{
    auto ans = *this;
    ans.sgn = false;
    return ans;
}
inline decimal decimal::operator+(decimal &&src)
{
    if(!(zero() || src.zero()))
    {
        auto l_coe_ft = dec_coe(ft, false), l_coe_it = dec_coe(it),
            r_coe_ft = src.dec_coe(ft, false), r_coe_it = src.dec_coe(it);
        if(l_coe_ft.length < r_coe_ft.length) for(auto i=l_coe_ft.length; i<r_coe_ft.length; ++i) l_coe_ft.emplace_back(0);
        if(l_coe_ft.length > r_coe_ft.length) for(auto i=r_coe_ft.length; i<l_coe_ft.length; ++i) r_coe_ft.emplace_back(0);
        l_coe_ft.reverse(); r_coe_ft.reverse();
        if(l_coe_it.length < r_coe_it.length) for(auto i=l_coe_it.length; i<r_coe_it.length; ++i) l_coe_it.emplace_back(0);
        if(l_coe_it.length > r_coe_it.length) for(auto i=r_coe_it.length; i<l_coe_it.length; ++i) r_coe_it.emplace_back(0);
        decimal ans;
        auto ans_coe = integer_polynomial_add(ans.sgn, l_coe_ft.unit(l_coe_it), r_coe_ft.unit(r_coe_it), sgn, src.sgn);
        ans.dec_coe(ft, ans_coe.sub_sequence(0, r_coe_ft.length-1));
        ans.dec_coe(it, ans_coe.sub_sequence(r_coe_ft.length, ans_coe.length-1));
        return ans;
    }
    else if(src.zero()) return decimal(*this);
    else if(zero()) return decimal(src);
    else return decimal();
}
inline decimal decimal::operator+(decimal &src) { return *this + std::move(src); }
inline decimal decimal::operator+(long double src) { return *this + decimal(src); }
inline void decimal::operator+=(decimal &&src) { *this = *this + std::move(src); }
inline void decimal::operator+=(decimal &src) { *this = *this + src; }
inline void decimal::operator+=(long double src) { *this = *this + src; }
inline decimal &decimal::operator++()
{
    *this += 1;
    return *this;
}
inline decimal decimal::operator++(int)
{
    auto temp = *this;
    *this += 1;
    return temp;
}
inline decimal decimal::operator-(decimal &&src)
{
    auto temp = src;
    temp.sgn = !temp.sgn;
    return *this + temp;
}
inline decimal decimal::operator-(decimal &src) { return *this - std::move(src); }
inline decimal decimal::operator-(long double src) { return *this - decimal(src); }
inline void decimal::operator-=(decimal &&src) { *this = *this - std::move(src); }
inline void decimal::operator-=(decimal &src) { *this = *this - src; }
inline void decimal::operator-=(long double src) { *this = *this - src; }
inline decimal &decimal::operator--()
{
    *this -= 1;
    return *this;
}
inline decimal decimal::operator--(int)
{
    auto temp = *this;
    *this -= 1;
    return temp;
}
inline decimal decimal::operator*(decimal &&src)
{
    if(!(zero() || src.zero()))
    {
        auto l_coe_ft = dec_coe(ft, false), l_coe_it = dec_coe(it, false),
            r_coe_ft = src.dec_coe(ft, false), r_coe_it = src.dec_coe(it, false);
        uint64_t coe_a_len = (l_coe_it.length ? l_coe_it.length : 1) + l_coe_ft.length,
                coe_b_len = (r_coe_it.length ? r_coe_it.length : 1) + r_coe_ft.length,
                coe_c_len_base = coe_a_len + coe_b_len - 1,
                coe_c_len_pad = num_pow_pad_cnt(coe_c_len_base, 2);
        net_sequence<uint64_t> coe_ans(0, 0);
        // ntt
        if(coe_c_len_pad>coe_c_len_base/2 && coe_c_len_base>=4)
        {
            auto coe_c_len = coe_c_len_base + coe_c_len_pad;
            auto cnt = 0ull;
            MEM_INIT(uint64_t, coe_a, coe_c_len);
            MEM_INIT(uint64_t, coe_b, coe_c_len);
            MEM_INIT(uint64_t, coe_c, coe_c_len);
            std::fill_n(coe_a, coe_c_len, 0);
            std::fill_n(coe_b, coe_c_len, 0);
            std::fill_n(coe_c, coe_c_len, 0);
            if(l_coe_ft.length) for(auto i=l_coe_ft.length; i; --i) coe_a[cnt++] = l_coe_ft[i-1];
            if(l_coe_it.length) for(auto i=l_coe_it.length; i; --i) coe_a[cnt++] = l_coe_it[i-1];
            cnt = 0;
            if(r_coe_ft.length) for(auto i=r_coe_ft.length; i; --i) coe_b[cnt++] = r_coe_ft[i-1];
            if(r_coe_it.length) for(auto i=r_coe_it.length; i; --i) coe_b[cnt++] = r_coe_it[i-1];
            if(integer_fnt(coe_c, coe_a, coe_b, coe_c_len) && integer_ifnt(coe_c, coe_c_len))
            coe_ans.ptr_array(std::move(coe_c), coe_c_len);
            MEM_RECYCLE(coe_a);
            MEM_RECYCLE(coe_b);
        }
        else
        {
            // normal
            auto coe_a = l_coe_it.unit(l_coe_ft), coe_b = r_coe_it.unit(r_coe_ft);
            coe_a.reverse();
            coe_b.reverse();
            coe_ans = integer_polynomial_mult(coe_a, coe_b);
        }
        net_sequence<uint8_t> ft_temp, it_temp;
        auto carry = 0ull;
        long long ft_bit_cnt = l_coe_ft.length + r_coe_ft.length;
        for(auto i=0ull; i<coe_ans.length; ++i)
        {
            carry += coe_ans[i];
            if(ft_bit_cnt-- > 0) ft_temp.insert(0, carry%10);
            else it_temp.insert(0, carry%10);
            carry /= 10;
        }
        while(ft_bit_cnt -- > 0)
        {
            ft_temp.insert(0, carry%10);
            carry /= 10;
        }
        while (carry)
        {
            it_temp.insert(0, carry%10);
            carry /= 10;
        }
        decimal ans;
        if(sgn != src.sgn) ans.sgn = true;
        ans.dec_coe(it, it_temp, false);
        ans.dec_coe(ft, ft_temp, false);
        return ans;
    }
    else return decimal();    
}
inline decimal decimal::operator*(decimal &src) { return *this * std::move(src); }
inline decimal decimal::operator*(long double src) { return *this * decimal(src); }
inline void decimal::operator*=(decimal &&src) { *this = *this * std::move(src); }
inline void decimal::operator*=(decimal &src) { *this = *this * src; }
inline void decimal::operator*=(long double src) { *this = *this * src; }
inline decimal decimal::operator/(decimal &&src) { return *this * src.reciprocal(); }
inline decimal decimal::operator/(decimal &src) { return *this / std::move(src); }
inline decimal decimal::operator/(long double src) { return *this / decimal(src); }
inline void decimal::operator/=(decimal &&src) { *this = *this / std::move(src); }
inline void decimal::operator/=(decimal &src) { *this = *this / src; }
inline void decimal::operator/=(long double src) { *this = *this / src; }
inline int64_t decimal::operator<<(int64_t bit) { return to_integer() << bit; }
inline decimal decimal::operator<<(decimal &&bit) { return *this << bit.to_integer(); }
inline decimal decimal::operator<<(decimal &bit) { return *this << std::move(bit); }
inline void decimal::operator<<=(int64_t bit) { *this = to_integer() << bit; }
inline void decimal::operator<<=(decimal &&bit) { *this <<= bit.to_integer(); }
inline void decimal::operator<<=(decimal &bit) { *this <<= std::move(bit); }
inline int64_t decimal::operator>>(int64_t bit) { return to_integer() >> bit; }
inline decimal decimal::operator>>(decimal &&bit) { return *this >> bit.to_integer(); }
inline decimal decimal::operator>>(decimal &bit) { return *this >> std::move(bit); }
inline void decimal::operator>>=(int64_t bit) { *this = to_integer() >> bit; }
inline void decimal::operator>>=(decimal &&bit) { *this >>= bit.to_integer(); }
inline void decimal::operator>>=(decimal &bit) { *this >>= std::move(bit); }
inline int64_t decimal::operator|(int64_t bit) { return to_integer() | bit; }
inline decimal decimal::operator|(decimal &&bit) { return *this | bit.to_integer(); }
inline decimal decimal::operator|(decimal &bit) { return *this | std::move(bit); }
inline void decimal::operator|=(int64_t bit) { *this = to_integer() | bit; }
inline void decimal::operator|=(decimal &&bit) { *this |= bit.to_integer(); }
inline void decimal::operator|=(decimal &bit) { *this |= std::move(bit); }
inline int64_t decimal::operator&(int64_t bit) { return to_integer() & bit; }
inline decimal decimal::operator&(decimal &&bit) { return *this & bit.to_integer(); }
inline decimal decimal::operator&(decimal &bit) { return *this & std::move(bit); }
inline void decimal::operator&=(int64_t bit) { *this = to_integer() & bit; }
inline void decimal::operator&=(decimal &&bit) { *this &= bit.to_integer(); }
inline void decimal::operator&=(decimal &bit) { *this &= std::move(bit); }
inline int64_t decimal::operator~() { return ~to_integer(); }
inline int64_t decimal::operator%(int64_t src) { assert(src); return to_integer() % src; }
inline decimal decimal::operator%(decimal &&src) { return *this % src.to_integer(); }
inline decimal decimal::operator%(decimal &src) { return *this % std::move(src); }
inline void decimal::operator%=(int64_t src) { *this %= to_integer() % src;  }
inline void decimal::operator%=(decimal &&src) { *this %= src.to_integer(); }
inline void decimal::operator%=(decimal &src) { *this %= std::move(src); }
inline void decimal::operator=(decimal &src) { value_copy(src); }
inline void decimal::operator=(const decimal &src) { __value_copy(src); }
inline void decimal::operator=(decimal &&src) { value_move(std::move(src)); }
inline void decimal::operator=(long double src) { *this = decimal(src); }
inline void decimal::operator=(std::string src) { *this = decimal(src); }
inline bool decimal::operator>(decimal &&src)
{
    if(src.sgn == sgn)
    {
        auto l_coe = dec_coe(it, false), r_coe = src.dec_coe(it, false);
        if(l_coe.length > r_coe.length)
        {
            if(sgn) return false;
            else return true;
        }
        else if(l_coe.length < r_coe.length)
        {
            if(sgn) return true;
            else return false;
        }
        else for(auto i=0ull; i<l_coe.length; ++i) if(l_coe[i] > r_coe[i])
        {
            if(sgn) return false;
            else return true;
        }
        else if(l_coe[i] < r_coe[i])
        {
            if(sgn) return true;
            else return false;
        }
        l_coe = dec_coe(ft, false);
        r_coe = src.dec_coe(ft, false);
        auto ft_dig_cnt = l_coe.length;
        if(r_coe.length < ft_dig_cnt) ft_dig_cnt = r_coe.length;
        for(auto i=0ull; i<ft_dig_cnt; ++i) if(l_coe[i] > r_coe[i])
        {
            if(sgn) return false;
            else return true;
        }
        else if(l_coe[i] < r_coe[i])
        {
            if(sgn) return true;
            else return false;
        }
        return l_coe.length > ft_dig_cnt;
    }
    else if(sgn) return false;
    else return true;
}
inline bool decimal::operator>(decimal &src) { return *this > std::move(src); }
inline bool decimal::operator>(long double src) { return *this > decimal(src); }
inline bool decimal::operator>=(decimal &&src) { return ((*this>std::move(src)) || (*this==std::move(src))); }
inline bool decimal::operator>=(decimal &src) { return *this >= std::move(src); }
inline bool decimal::operator>=(long double src) { return *this >= decimal(src); }
inline bool decimal::operator<(decimal &&src) { return !((*this>src) || (*this==src)); }
inline bool decimal::operator<(decimal &src) { return *this < std::move(src); }
inline bool decimal::operator<(long double src) { return *this < decimal(src); }
inline bool decimal::operator<=(decimal &&src) { return ((*this<std::move(src)) || (*this==std::move(src))); }
inline bool decimal::operator<=(decimal &src) { return *this<=std::move(src); }
inline bool decimal::operator<=(long double src) { return *this <= decimal(src); }
inline bool decimal::operator==(decimal &&src)
{
    if(val[ft].length==src.val[ft].length && val[it].length==src.val[it].length && sgn==src.sgn)
    {
        for(auto i=0ull; i<val[ft].length; ++i) if(val[ft][i] != src.val[ft][i]) return false;
        for(auto i=0ull; i<val[it].length; ++i) if(val[it][i] != src.val[it][i]) return false;
        return true;
    }
    else return false;
}
inline bool decimal::operator==(decimal &src) { return *this == std::move(src); }
inline bool decimal::operator==(long double src) { return *this == decimal(src); }
inline bool decimal::operator!=(decimal &&src) { return !(*this == std::move(src)); }
inline bool decimal::operator!=(decimal &src) { return *this != std::move(src); }
inline bool decimal::operator!=(long double src) { return *this != decimal(src); }
inline void decimal::reset() { sgn = false; val[it].reset(); val[ft].reset(); }
inline decimal::~decimal() { reset(); }
inline decimal decimal::acc_convert(uint64_t con_acc)
{
    decimal ans;
    if(con_acc)
    {
        std::string temp_str = "0.";
        while(--con_acc) temp_str.push_back('0');
        temp_str.push_back('1');
        ans = temp_str;
    }
    return ans;
}
inline decimal decimal::reciprocal()
{
    if(zero()) return 0;
    else
    {
        decimal curr = 0, deriv = 2, convergency = acc_convert(calculate_digit), temp = *this;
        temp.sgn = false;
        std::string next_str = "0.";
        auto temp_it_len = it_len();
        if(temp_it_len) while(-- temp_it_len) next_str.push_back('0');
        next_str.push_back('1');
        decimal next = next_str;
        do
        {
            curr = next;
            next = curr * (deriv - temp * curr);
            next._acc(calculate_digit);
        }
        while ((next-curr).abs() > convergency);
        if(sgn) curr.sgn = true;
        return curr;
    }
}
inline decimal decimal::sin(decimal &src)
{
    decimal ans = src, pow_num = src, curr_ans = src, frct_num = 1, coe = -1, cnt = 1,convergency = acc_convert(calculate_digit);
    do
    {
        ans = curr_ans;
        pow_num *= src * src;
        pow_num._acc(calculate_digit);
        frct_num *= (cnt + 1) * (cnt + 2);
        curr_ans = ans + coe * pow_num * frct_num.reciprocal();
        curr_ans._acc(calculate_digit);
        coe.sgn = !coe.sgn;
        cnt += 2;
    } while ((curr_ans-ans).abs() > convergency);
    return ans;
}
inline decimal decimal::cos(decimal &src)
{
    decimal ans = 1, pow_num = src * src, curr_ans = 1, frct_num = 1, coe = -1, cnt = 0, convergency = acc_convert(calculate_digit);
    do
    {
        ans = curr_ans;
        frct_num *= (cnt + 1) * (cnt + 2);
        curr_ans = ans + coe * pow_num * frct_num.reciprocal();
        curr_ans._acc(calculate_digit);
        pow_num *= src * src;
        pow_num._acc(calculate_digit);
        cnt += 2;
        coe.sgn = !coe.sgn;
    }
    while ((curr_ans-ans).abs() > convergency);
    return ans;
}
inline decimal decimal::ln(decimal &src)
{
    if(src > 0)
    {
        src = (src - 1) * (src + 1).reciprocal();
        auto ans = src, curr_ans = src, pow_num = src, convergency = acc_convert(calculate_digit), cnt = decimal(1);
        do
        {
            ans = curr_ans;
            cnt += 2;
            pow_num *= src * src;
            pow_num._acc(calculate_digit);
            curr_ans = ans + pow_num * cnt.reciprocal();
            curr_ans._acc(calculate_digit);
        } while ((curr_ans-ans).abs() > convergency);
        return 2 * ans;
    }
    else return 0;
}
inline decimal decimal::exp(decimal &src)
{
    decimal ans = 1, pow_num = src, curr_ans = 1, fact_num = 1, cnt = 1, convergency = acc_convert(calculate_digit);
    do
    {
        ans = curr_ans;
        curr_ans = ans + pow_num * fact_num.reciprocal();
        curr_ans._acc(calculate_digit);
        ++ cnt;
        fact_num *= cnt;
        pow_num *= src;
        pow_num._acc(calculate_digit);
    }
    while ((curr_ans-ans).abs() > convergency);
    return ans;
}
inline decimal decimal::euler_itr(decimal &times, decimal &k) { return (2 * k + 1) * decimal(pi) * times; }
inline net_sequence<std::complex<decimal>> decimal::euler_eqt(decimal &times)
{
    net_sequence<std::complex<decimal>> ans;
    auto convergency = acc_convert(calculate_digit);
    decimal cnt = 0;
    auto flag = true;
    auto limit = (3 * times + 2) * decimal(pi);
    decimal curr_itr = 0;
    do
    {
        curr_itr = euler_itr(times, cnt);
        auto real = cos(curr_itr), img = sin(curr_itr);
        std::complex<decimal> curr_ans(real, img);
        ans.emplace_back(std::move(curr_ans));
        ++ cnt;
    } while (curr_itr.abs() < limit.abs());
    return ans;
}
inline decimal decimal::power(decimal &times)
{
    if(zero()) return 0;
    else
    {
        auto temp = *this;
        auto sub_acc = calculate_digit / 2;
        if(times.val[ft].length)
            if(*this > 0)
            {
                auto ans = decimal::exp(times*decimal::ln(temp));
                ans._acc(sub_acc);
                return ans;
            }
            else
            {
                auto euler_ans = euler_eqt(times);
                if(euler_ans.length) for(auto i=0; i<euler_ans.length; ++i)
                {
                    auto curr_i = euler_ans[i].imag();
                    curr_i._acc(8);
                    if(curr_i == 0)
                    {
                        temp.sgn = !temp.sgn;
                        auto curr_r = euler_ans[i].real();
                        auto ans =  decimal::exp(times*decimal::ln(temp))*euler_ans[i].real();
                        ans._acc(sub_acc);
                        return ans;
                    }
                }
                return 0;
            }
        else if(val[ft].length)
        {
            auto temp = times.to_integer();
            decimal ans = *this;
            for(auto i=1; i<temp; ++i) ans *= (*this);
            if(times.sgn) ans = ans.reciprocal();
            return ans;
        }
        else return integer_euler_pow(to_integer(), times.to_integer());
    }
}
inline decimal decimal::operator^(decimal &&src) { return power(src); }
inline decimal decimal::operator^(decimal &src) { return *this ^ std::move(src); }
inline decimal decimal::operator^(long double src) { return (*this) ^ std::move(decimal(src)); }
inline decimal operator^(long double src, decimal &val) { return decimal(src) ^ std::move(val); }

decimal __power_(decimal &base, decimal &times) { return base.power(times); }
long double __power_(long double base, long double times) { return std::pow(base, times); }

int64_t __to_integer_(decimal &src) { src.accuracy = 0; return src.to_integer(); }
int64_t __to_integer_(long double src) { src < 0? src -= 0.5 : src += 0.5; return (int64_t)src; }

long double __to_float_point_(decimal src) { return src.to_float(); }
long double __to_float_point_(int64_t src) { return src; }

BAGRT_END