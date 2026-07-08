#pragma once

// MortonTrieGeneric<T>: per-instance, templated Morton trie used by MortonTrieMap<T>
// in the C_MAP_TYPE=3 (USE_MORTONTRIE) path. Replaces the legacy singleton MortonTrie
// (which is double-only and shared across all map instances) for the LBM hot path.
//
// Layout:
//   * Internal nodes: TrieNodeT<T>::children[BRANCH] allocated lazily
//   * Leaf nodes:    TrieNodeT<T>::data[BRANCH] (one T per cell in the leaf's BRANCH bucket)
//   * One existence bitmap per instance for O(1) exists() (per-instance state, ~2.8 MB)
//   * One search_path_stack for cursor / search_halfway path-amortized lookups
//
// The legacy MortonTrie singleton (mortonTrie.hpp / mortonTrie.cpp) is preserved
// untouched so the existing test/ablation programs keep building.

#include "settings.hpp"
#include "General.h"
#include "Morton_Assist.h"
#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

template <typename T>
struct TrieNodeT
{
    TrieNodeT* children = nullptr;
    T* data = nullptr;
#ifdef ABPattern
    T* data2 = nullptr;
#endif
    uint8_t level = 0;
};

template <typename T>
class MortonTrieGeneric
{
public:
    using value_type = T;

    struct Cursor
    {
        TrieNodeT<T>* leaf = nullptr;
        uint8_t idx_in_leaf = 0;
        bool valid() const { return leaf != nullptr; }
    };

    MortonTrieGeneric()
        : depth_(static_cast<uint8_t>(std::log2(max_of_three<uint>(Nx, Ny, Nz)) + 1 + refine_level)),
          root_(new TrieNodeT<T>()),
          search_path_stack_(depth_, nullptr),
          stack_code_(0),
          bitmap_storage_(MAX_BG_MORTON + 1)
    {
    }

    ~MortonTrieGeneric()
    {
        free_subtree(root_);
        delete root_;
    }

    // Deep copy: the trie owns raw tree pointers and a per-instance bitmap.
    // search_path_stack_ is reset (lazily rebuilt on next cursor_halfway); the
    // bitmap is value-copied. Required so MortonTrieMap can be a drop-in for
    // std::unordered_map (which is copyable) under C_MAP_TYPE=3.
    MortonTrieGeneric(const MortonTrieGeneric& other)
        : depth_(other.depth_),
          root_(new TrieNodeT<T>()),
          search_path_stack_(other.depth_, nullptr),
          stack_code_(0),
          bitmap_storage_(other.bitmap_storage_)
    {
        clone_node(root_, other.root_);
    }

    MortonTrieGeneric& operator=(const MortonTrieGeneric& other)
    {
        if (this != &other)
        {
            free_subtree(root_);
            delete root_;
            root_ = new TrieNodeT<T>();
            clone_node(root_, other.root_);
            depth_ = other.depth_;
            search_path_stack_.assign(other.depth_, nullptr);
            stack_code_ = 0;
            bitmap_storage_ = other.bitmap_storage_;
        }
        return *this;
    }

    // O(1) whole-trie swap: exchange root pointer + bitmap + path stack. Used by
    // MortonTrieMap::swap to back the STL "swap-with-temp-to-clear" idiom (e.g.
    // Lat_Voxelize.cpp's `temp.swap(lat_overlap_C2F.at(l))`). All members are
    // exchanged together so the stack_code_ <-> search_path_stack_ invariant
    // stays consistent on both sides after the swap.
    void swap(MortonTrieGeneric& other) noexcept
    {
        using std::swap;
        swap(depth_, other.depth_);
        swap(root_, other.root_);
        search_path_stack_.swap(other.search_path_stack_);
        swap(stack_code_, other.stack_code_);
        bitmap_storage_.swap(other.bitmap_storage_);
    }

    // Soft erase: clear the existence bit (so exists/find/count return false)
    // and reset the data slot to T{} so any stale cursor read returns a
    // default value. Trie structure is kept (no path collapse); memory is
    // reclaimed only when the whole map is destroyed. Acceptable for AMR where
    // these maps are short-lived.
    void erase(const D_morton& morton, uint8_t level)
    {
        if (!exists(morton, level))
            return;
        const uint64_t morton_val = morton.to_ullong();
        constexpr uint32_t bg_bits_start = refine_level * DIM;
        constexpr uint64_t bg_mask = (1ULL << (BG_DEPTH * DIM)) - 1;
        const uint64_t bg_morton = (morton_val >> bg_bits_start) & bg_mask;
        constexpr uint64_t ref_mask = (1ULL << (refine_level * DIM)) - 1;
        const uint32_t ref_index = static_cast<uint32_t>(morton_val & ref_mask);
        if (__builtin_expect(bg_morton <= MAX_BG_MORTON && ref_index < BITMAP_BITS_PER_CELL, 1))
            bitmap_storage_[bg_morton][ref_index] = false;
        Cursor c = build_path_to_leaf(morton, level, /*create=*/false);
        if (c.valid())
            c.leaf->data[c.idx_in_leaf] = T{};
    }

    uint8_t depth() const { return depth_; }
    TrieNodeT<T>* root() const { return root_; }
    size_t bitmap_size() const { return bitmap_storage_.size(); }

    // Bitmap-based O(1) existence query
    static constexpr uint32_t BITMAP_BITS_PER_CELL = 1u << (refine_level * DIM);
    static constexpr uint32_t MAX_BG_MORTON = 5446;
    static constexpr uint32_t BG_DEPTH = 6;
    using CellBitmap = std::bitset<BITMAP_BITS_PER_CELL>;

    inline bool exists(const D_morton& morton, uint8_t /*level*/) const
    {
        const uint64_t morton_val = morton.to_ullong();
        constexpr uint32_t bg_bits_start = refine_level * DIM;
        constexpr uint64_t bg_mask = (1ULL << (BG_DEPTH * DIM)) - 1;
        const uint64_t bg_morton = (morton_val >> bg_bits_start) & bg_mask;
        constexpr uint64_t ref_mask = (1ULL << (refine_level * DIM)) - 1;
        const uint32_t ref_index = static_cast<uint32_t>(morton_val & ref_mask);
        if (__builtin_expect(bg_morton > MAX_BG_MORTON || ref_index >= BITMAP_BITS_PER_CELL, 0))
            return false;
        return bitmap_storage_[bg_morton][ref_index];
    }

    // Insert / at / update
    void insert(const D_morton& morton, uint8_t level, const T& value)
    {
        Cursor c = build_path_to_leaf(morton, level, /*create=*/true);
        c.leaf->data[c.idx_in_leaf] = value;
        set_existence_bit(morton);
    }

    // Returns reference to T at (morton, level). Inserts default-constructed T
    // if absent; matches std::unordered_map::operator[] semantics.
    T& touch(const D_morton& morton, uint8_t level)
    {
        Cursor c = build_path_to_leaf(morton, level, /*create=*/true);
        if (!exists(morton, level))
        {
            c.leaf->data[c.idx_in_leaf] = T{};
            set_existence_bit(morton);
        }
        return c.leaf->data[c.idx_in_leaf];
    }

    // Returns reference, throws if not present (std::unordered_map::at semantics).
    T& at(const D_morton& morton, uint8_t level)
    {
        Cursor c = build_path_to_leaf(morton, level, /*create=*/false);
        if (!c.valid())
            throw std::out_of_range("MortonTrieGeneric::at: key not found");
        return c.leaf->data[c.idx_in_leaf];
    }

    const T& at(const D_morton& morton, uint8_t level) const
    {
        Cursor c = const_cast<MortonTrieGeneric*>(this)->build_path_to_leaf(morton, level, /*create=*/false);
        if (!c.valid())
            throw std::out_of_range("MortonTrieGeneric::at: key not found");
        return c.leaf->data[c.idx_in_leaf];
    }

    // Cursor API: amortized walks for batched access
    // Walk from root once, return handle. Caller can then read/write via
    // c.leaf->data[c.idx_in_leaf] directly (zero further trie cost).
    Cursor cursor(const D_morton& morton, uint8_t level)
    {
        return build_path_to_leaf(morton, level, /*create=*/false);
    }

    // Path-amortized lookup: reuse the path stack from the last call, descend
    // only from the common-ancestor of `stack_code_` and `morton`. Mirrors the
    // legacy MortonTrie::search_halfway algorithm but generalised to template T.
    Cursor cursor_halfway(const D_morton& morton, uint8_t level)
    {
        const uint8_t search_depth = depth_ - refine_level + level;
        const D_morton xor_diff = stack_code_ ^ morton;
        const uint64_t xor_ull = xor_diff.to_ullong();

        // Fully reuse path: same morton as last time.
        if (xor_ull == 0 && search_path_stack_[search_depth - 2] != nullptr)
        {
            Cursor c;
            c.leaf = search_path_stack_[search_depth - 2];
            c.idx_in_leaf = EXTRACT_MORTON_CONVERT_INDEX(morton, DIM);
            return c;
        }

        // Find divergence depth from the high (root) end of the morton code.
        int start_search_depth =
            (xor_ull == 0)
                ? -1
                : (static_cast<int>(__builtin_clzll(xor_ull)) - (64 - DIM * static_cast<int>(depth_))) / DIM - 1;

        // Cold-start or stale-stack safety: if the path-stack entry we'd resume
        // from is null (no prior cursor() call on this instance, or after clear),
        // fall back to a full descent from root.
        if (start_search_depth >= 0 && search_path_stack_[start_search_depth] == nullptr)
            return build_path_to_leaf(morton, level, /*create=*/false);

        TrieNodeT<T>* node = (start_search_depth < 0) ? root_ : search_path_stack_[start_search_depth];

        for (int i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth)
        {
            const uint8_t start_pos = (depth_ - i_depth) * DIM;
            const uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
            if (!node->children)
                return Cursor{};
            node = &node->children[index];
            search_path_stack_[i_depth] = node;
        }

        if (!node->data)
            return Cursor{};

        stack_code_ = morton;
        Cursor c;
        c.leaf = node;
        c.idx_in_leaf = EXTRACT_MORTON_CONVERT_INDEX(morton, DIM);
        return c;
    }

    // Walk to a neighbour of the current node, reusing the path stack. The
    // neighbour morton code is computed via Morton_Assist::find_neighbor.
    Cursor cursor_neighbor(const D_morton& code,
                           uint8_t level,
                           int x_shift, int y_shift, int z_shift)
    {
        const D_morton ncode = Morton_Assist::find_neighbor(code, level, x_shift, y_shift, z_shift);
        return cursor_halfway(ncode, level);
    }

    T& deref(const Cursor& c)
    {
        return c.leaf->data[c.idx_in_leaf];
    }

    const T& deref(const Cursor& c) const
    {
        return c.leaf->data[c.idx_in_leaf];
    }

#ifdef ABPattern
    T& deref2(const Cursor& c)
    {
        return c.leaf->data2[c.idx_in_leaf];
    }
    const T& deref2(const Cursor& c) const
    {
        return c.leaf->data2[c.idx_in_leaf];
    }
#endif

private:
    uint8_t depth_;
    TrieNodeT<T>* root_;
    std::vector<TrieNodeT<T>*> search_path_stack_;
    D_morton stack_code_;
    std::vector<CellBitmap> bitmap_storage_;

    inline void set_existence_bit(const D_morton& morton)
    {
        const uint64_t morton_val = morton.to_ullong();
        constexpr uint32_t bg_bits_start = refine_level * DIM;
        constexpr uint64_t bg_mask = (1ULL << (BG_DEPTH * DIM)) - 1;
        const uint64_t bg_morton = (morton_val >> bg_bits_start) & bg_mask;
        constexpr uint64_t ref_mask = (1ULL << (refine_level * DIM)) - 1;
        const uint32_t ref_index = static_cast<uint32_t>(morton_val & ref_mask);
        if (__builtin_expect(bg_morton <= MAX_BG_MORTON && ref_index < BITMAP_BITS_PER_CELL, 1))
            bitmap_storage_[bg_morton][ref_index] = true;
    }

    // Build (or follow, if !create) the trie path down to the leaf containing
    // the cell for `morton`. Returns a Cursor pointing at leaf + idx_in_leaf.
    Cursor build_path_to_leaf(const D_morton& morton, uint8_t level, bool create)
    {
        TrieNodeT<T>* node = root_;
        const uint8_t search_depth = depth_ - refine_level + level;

        for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
        {
            const uint8_t start_pos = (depth_ - i_depth) * DIM;
            const uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
            if (!node->children)
            {
                if (!create)
                    return Cursor{};
                node->children = new TrieNodeT<T>[BRANCH];
                node->children[index].level = i_depth;
            }
            node = &node->children[index];
            search_path_stack_[i_depth] = node;
        }

        if (!node->data)
        {
            if (!create)
                return Cursor{};
            node->data = new T[BRANCH]();
#ifdef ABPattern
            node->data2 = new T[BRANCH]();
#endif
        }

        stack_code_ = morton;
        Cursor c;
        c.leaf = node;
        c.idx_in_leaf = EXTRACT_MORTON_CONVERT_INDEX(morton, DIM);
        return c;
    }

    void free_subtree(TrieNodeT<T>* node)
    {
        if (!node) return;
        if (node->children)
        {
            for (uint8_t i = 0; i < BRANCH; ++i)
                free_subtree(&node->children[i]);
            delete[] node->children;
            node->children = nullptr;
        }
        if (node->data)
        {
            delete[] node->data;
            node->data = nullptr;
        }
#ifdef ABPattern
        if (node->data2)
        {
            delete[] node->data2;
            node->data2 = nullptr;
        }
#endif
    }

    // Deep-copy a single TrieNodeT (and its subtree). dst must be a freshly
    // default-constructed node (children/data null); we allocate child arrays
    // and copy element-by-element so levels and ABPattern data2 stay in sync.
    static void clone_node(TrieNodeT<T>* dst, const TrieNodeT<T>* src)
    {
        dst->level = src->level;
        dst->children = nullptr;
        dst->data = nullptr;
#ifdef ABPattern
        dst->data2 = nullptr;
#endif
        if (src->children)
        {
            dst->children = new TrieNodeT<T>[BRANCH];
            for (uint8_t i = 0; i < BRANCH; ++i)
                clone_node(&dst->children[i], &src->children[i]);
        }
        if (src->data)
        {
            dst->data = new T[BRANCH];
            for (uint8_t i = 0; i < BRANCH; ++i)
                dst->data[i] = src->data[i];
        }
#ifdef ABPattern
        if (src->data2)
        {
            dst->data2 = new T[BRANCH];
            for (uint8_t i = 0; i < BRANCH; ++i)
                dst->data2[i] = src->data2[i];
        }
#endif
    }
};
