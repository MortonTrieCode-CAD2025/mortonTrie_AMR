// MortonTrieMap<T>: STL-compatible wrapper around MortonTrieGeneric<T>.
// Routes find/at/operator[]/count/insert/exists through trie_ directly.
// Uses lazy Morton-order bitmap cache for begin()/end() iteration only.

#include "header/mortonTrie.hpp"
#include "header/Grid_Class.h"   // Node, Node_IB, Infor_near_solid, Node*
#include "header/user.h"         // Coord<T>, D_uvw = Coord<D_Phy_Velocity>
#include "header/Lat_Manager.hpp" // Cell
#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>
#include <vector>

// Lazy Morton-order cache for begin()/end() iteration

template <typename Value>
void MortonTrieMap<Value>::update_cache() const
{
    if (cache_valid_) return;

    cached_data_.clear();

    using Trie = MortonTrieGeneric<Value>;
    constexpr uint32_t bg_bits_start = refine_level * DIM;

    // Walk existence bitmap in Morton order → sorted output
    for (uint64_t bg_morton = 0; bg_morton <= Trie::MAX_BG_MORTON; ++bg_morton)
    {
        for (uint32_t ref_index = 0; ref_index < Trie::BITMAP_BITS_PER_CELL; ++ref_index)
        {
            D_morton code = D_morton((bg_morton << bg_bits_start) | ref_index);
            if (!trie_.exists(code, default_level_))
                continue;
            auto c = trie_.cursor(code, default_level_);
            if (!c.valid())
                continue;
            cached_data_.emplace_back(code, trie_.deref(c));
        }
    }

    cache_valid_ = true;
}

template <typename Value>
typename MortonTrieMap<Value>::iterator MortonTrieMap<Value>::begin() const
{
    update_cache();
    if (cached_data_.empty())
        return end();  // avoid deref nullptr on empty map
    auto& mutable_cache =
        const_cast<std::vector<std::pair<D_morton, Value>>&>(cached_data_);
    return iterator(const_cast<MortonTrieMap<Value>*>(this),
                    mutable_cache.begin(), mutable_cache.end());
}

template <typename Value>
typename MortonTrieMap<Value>::iterator MortonTrieMap<Value>::end() const
{
    return iterator(const_cast<MortonTrieMap<Value>*>(this), nullptr);
}

template <typename Value>
bool MortonTrieMap<Value>::empty() const
{
    update_cache();
    return cached_data_.empty();
}

template <typename Value>
typename MortonTrieMap<Value>::size_type MortonTrieMap<Value>::size() const
{
    update_cache();
    return cached_data_.size();
}

template <typename Value>
void MortonTrieMap<Value>::clear()
{
    cached_data_.clear();
    invalidate_cache();
}

template <typename Value>
std::pair<typename MortonTrieMap<Value>::iterator, bool>
MortonTrieMap<Value>::insert(const value_type& value)
{
    const bool already = trie_.exists(value.first, default_level_);
    if (already)
        return std::make_pair(find(value.first), false);
    trie_.insert(value.first, default_level_, value.second);
    invalidate_cache();
    return std::make_pair(find(value.first), true);
}

template <typename Value>
typename MortonTrieMap<Value>::size_type
MortonTrieMap<Value>::erase(const D_morton& key)
{
    if (!trie_.exists(key, default_level_))
        return 0;
    trie_.erase(key, default_level_);
    invalidate_cache();
    return 1;
}

template <typename Value>
template <typename... Args>
std::pair<typename MortonTrieMap<Value>::iterator, bool>
MortonTrieMap<Value>::emplace(Args&&... args)
{
    value_type value(std::forward<Args>(args)...);
    return insert(value);
}

template <typename Value>
typename MortonTrieMap<Value>::iterator
MortonTrieMap<Value>::find(const D_morton& key) const
{
    if (!trie_.exists(key, default_level_))
        return end();
    auto c = trie_.cursor(key, default_level_);
    if (!c.valid())
        return end();
    return iterator(const_cast<MortonTrieMap<Value>*>(this), c, key);
}

template <typename Value>
typename MortonTrieMap<Value>::size_type
MortonTrieMap<Value>::count(const D_morton& key) const
{
    return trie_.exists(key, default_level_) ? 1 : 0;
}

template <typename Value>
Value& MortonTrieMap<Value>::operator[](const D_morton& key)
{
    const bool existed = trie_.exists(key, default_level_);
    Value& slot = trie_.touch(key, default_level_);
    if (!existed)
        invalidate_cache();
    return slot;
}

template <typename Value>
Value& MortonTrieMap<Value>::at(const D_morton& key)
{
    if (!trie_.exists(key, default_level_))
        throw std::out_of_range("MortonTrieMap::at: key not found");
    return trie_.at(key, default_level_);
}

template <typename Value>
const Value& MortonTrieMap<Value>::at(const D_morton& key) const
{
    if (!trie_.exists(key, default_level_))
        throw std::out_of_range("MortonTrieMap::at: key not found");
    return trie_.at(key, default_level_);
}

template <typename Value>
template <typename InputIt>
void MortonTrieMap<Value>::insert(InputIt first, InputIt last)
{
    for (auto it = first; it != last; ++it)
        insert(*it);
}

// Explicit template instantiations
template class MortonTrieMap<int>;
template class MortonTrieMap<float>;
template class MortonTrieMap<double>;
template class MortonTrieMap<Node>;
template class MortonTrieMap<Node_IB>;
template class MortonTrieMap<Coord<double>>;            // D_uvw
template class MortonTrieMap<std::array<int, 2>>;       // D_mapint2

// Additional instantiations required by the LBM build (see link errors).
template class MortonTrieMap<bool>;                                   // visited flags
template class MortonTrieMap<Node*>;                                  // D_mapNodePtr
template class MortonTrieMap<Infor_near_solid>;                       // D_map_Infor_near_solid
template class MortonTrieMap<Cell>;                                   // D_mapLat
template class MortonTrieMap<std::array<D_real, C_Q-1>>;              // Lat_Manager::dis
template class MortonTrieMap<std::vector<std::array<D_int, 2>>>;      // triFaceID_in_srf

// Member-template insert(InputIt, InputIt) explicit instantiations
template void MortonTrieMap<int>::insert<MortonTrieIterator<int>>(
    MortonTrieIterator<int>, MortonTrieIterator<int>);
template void MortonTrieMap<double>::insert<MortonTrieIterator<double>>(
    MortonTrieIterator<double>, MortonTrieIterator<double>);
template void MortonTrieMap<Coord<double>>::insert<MortonTrieIterator<Coord<double>>>(
    MortonTrieIterator<Coord<double>>, MortonTrieIterator<Coord<double>>);
