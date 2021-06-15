template<typename T>
struct Array {
    // using Type = typename std::conditional<std::is_const<T>::value, const T, T>::type;
    // using T_Nonconst = typename std::remove_const<T>::type;
    T *data;
    int size;
    int capacity;

    T& operator[](size_t idx) {
        internal_check_limits(idx);
        return data[idx];
    }

    const T& operator[](size_t idx) const {
        internal_check_limits(idx);
        return data[idx];
    }

    void add(const T& elem) {
        assert(size < capacity);
        data[size] = elem;
        ++size;
    }

    void internal_check_limits(size_t idx) const {
        assert(idx <= 0x7FFFFFFF);
        assert((int)idx < size);
    }

    Array() {}
    Array(T *arr, int n, int cap = 0x7FFFFFFF) : data(arr), size(n), capacity((cap == 0x7FFFFFFF) ? n : cap) {}
    Array(const Array<T>& arr) = delete;
};

#define For(Arr) \
    for (auto [it_index, it] = std::tuple{0, Arr.data[0]}; it_index < Arr.size; ++it_index, it = Arr.data[it_index])

#define FOR(It, Min, Max) assert(Min <= Max); for (int It = Min; It != Max; ++It)

// #define For(item, arr) \
//     for (auto [it_index, item] = std::tuple{0, arr.data[0]}; it_index < arr.size; ++it_index, item = arr.data[it_index])
