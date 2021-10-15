#pragma once
#include<memory>

template<class T>
class PimplPtr {
public:
    using ElementType = typename std::unique_ptr<T>::element_type;

    PimplPtr() : p_(std::make_unique<T>()) {
        static_assert(sizeof(T) > 0, "Probably, you forgot to declare constructor explicitly");
    }
    explicit PimplPtr(std::unique_ptr<T>&& p) : p_(std::move(p)) {}

    PimplPtr(PimplPtr&&) noexcept = default;
    PimplPtr& operator =(PimplPtr&&) noexcept = default;

    ~PimplPtr() {
        static_assert(sizeof(T) > 0, "Probably, you forgot to declare destructor explicitly");
    }

    const ElementType* get() const noexcept { return p_.get(); }
    const ElementType* operator->() const noexcept { return get(); }
    const ElementType& operator*() const noexcept { return *get(); }
    explicit operator const ElementType* () const noexcept { return get(); }

    ElementType* get() noexcept { return p_.get(); }
    ElementType* operator->() noexcept { return get(); }
    ElementType& operator*() noexcept { return *get(); }
    explicit operator ElementType* () noexcept { return get(); }

private:
    std::unique_ptr<T> p_;
};