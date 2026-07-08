#pragma once

// MortonTrieDDF: dedicated Morton trie for the LBM distribution function
// arrays f[C_Q] / fcol[C_Q]. All C_Q directions of a single field share one
// trie skeleton + bitmap, and each leaf stores `double[BRANCH][C_Q]` so that
// the 27 directions of one cell are contiguous in memory.
//
// Why this layout (vs. C_Q independent tries):
//   * Collision hot path (Collision_Cumulant.cpp:37-40) reads all 27 directions
//     of a single cell into a local ddf[27] buffer. With one shared trie + fat
//     leaf, this is 1 trie walk + 27 contiguous double loads (~4 cache lines).
//     With 27 independent tries it would be 27 cold trie walks.
//   * Stream write of the current cell's 27 directions is the same: 1 walk +
//     27 contiguous stores.
//   * Stream read of 27 different neighbours benefits from cursor_halfway
//     amortisation, which reuses the path stack from the previous neighbour
//     based on the XOR-diff common ancestor.
//   * Bitmap cost is paid once per field (~2.8 MB) instead of C_Q times
//     (~150 MB).
//
// Use via MortonTrieDDFSlot views: each f[i_q] is a thin (parent, slot_idx)
// wrapper modelling the std::unordered_map<D_morton, double> interface, so
// existing LBM call sites of the form `f[i_q].at(code)` keep working.

#include "settings.hpp"
#include "General.h"
#include "Morton_Assist.h"
#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

class MortonTrieDDF
{
public:
    static constexpr uint32_t Q_DIM = C_Q;          // distribution-function directions
    static constexpr uint32_t BRANCH_LOCAL = BRANCH; // 8 in 3D

    // Leaf payload: 27 doubles per cell × 8 cells in the leaf bucket.
    // Layout is [idx_in_leaf][direction] so the 27 directions of a single
    // cell are contiguous (~216 bytes, ~4 cache lines).
    struct LeafData
    {
        double f[BRANCH_LOCAL][Q_DIM] = {};
#ifdef ABPattern
        double f2[BRANCH_LOCAL][Q_DIM] = {};
#endif
    };

    struct DDFTrieNode
    {
        DDFTrieNode* children = nullptr;
        LeafData* data = nullptr; // only set on leaf nodes
        uint8_t level = 0;
    };

    struct Cursor
    {
        DDFTrieNode* leaf = nullptr;
        uint8_t idx_in_leaf = 0;
        bool valid() const { return leaf != nullptr; }
    };

    MortonTrieDDF()
        : depth_(static_cast<uint8_t>(std::log2(max_of_three<uint>(Nx, Ny, Nz)) + 1 + refine_level)),
          root_(new DDFTrieNode()),
          search_path_stack_(depth_, nullptr),
          stack_code_(0),
          bitmap_storage_(MAX_BG_MORTON + 1)
    {
    }

    ~MortonTrieDDF()
    {
        free_subtree(root_);
        delete root_;
    }

    MortonTrieDDF(const MortonTrieDDF&) = delete;
    MortonTrieDDF& operator=(const MortonTrieDDF&) = delete;

    uint8_t depth() const { return depth_; }

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

    // Core API: cursor + per-slot access
    Cursor cursor(const D_morton& morton, uint8_t level)
    {
        return build_path_to_leaf(morton, level, /*create=*/false);
    }

    Cursor cursor_create(const D_morton& morton, uint8_t level)
    {
        return build_path_to_leaf(morton, level, /*create=*/true);
    }

    // Amortised neighbour cursor (Stream hot path).
    Cursor cursor_neighbor(const D_morton& code, uint8_t level,
                           int x_shift, int y_shift, int z_shift)
    {
        const D_morton ncode = Morton_Assist::find_neighbor(code, level, x_shift, y_shift, z_shift);
        return cursor_halfway(ncode, level);
    }

    Cursor cursor_halfway(const D_morton& morton, uint8_t level)
    {
        const uint8_t search_depth = depth_ - refine_level + level;
        const D_morton xor_diff = stack_code_ ^ morton;
        const uint64_t xor_ull = xor_diff.to_ullong();

        if (xor_ull == 0 && search_path_stack_[search_depth - 2] != nullptr)
        {
            Cursor c;
            c.leaf = search_path_stack_[search_depth - 2];
            c.idx_in_leaf = EXTRACT_MORTON_CONVERT_INDEX(morton, DIM);
            return c;
        }

        int start_search_depth =
            (xor_ull == 0)
                ? -1
                : (static_cast<int>(__builtin_clzll(xor_ull)) - (64 - DIM * static_cast<int>(depth_))) / DIM - 1;

        if (start_search_depth >= 0 && search_path_stack_[start_search_depth] == nullptr)
            return build_path_to_leaf(morton, level, /*create=*/false);

        DDFTrieNode* node = (start_search_depth < 0) ? root_ : search_path_stack_[start_search_depth];

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

    // Single-direction access via cursor.
    inline double& deref_slot(const Cursor& c, uint8_t slot)
    {
        return c.leaf->data->f[c.idx_in_leaf][slot];
    }

    inline const double& deref_slot(const Cursor& c, uint8_t slot) const
    {
        return c.leaf->data->f[c.idx_in_leaf][slot];
    }

#ifdef ABPattern
    inline double& deref_slot2(const Cursor& c, uint8_t slot)
    {
        return c.leaf->data->f2[c.idx_in_leaf][slot];
    }
    inline const double& deref_slot2(const Cursor& c, uint8_t slot) const
    {
        return c.leaf->data->f2[c.idx_in_leaf][slot];
    }
#endif

    // Convenience: one trie walk, single slot.
    double& at_slot(const D_morton& morton, uint8_t level, uint8_t slot)
    {
        Cursor c = cursor(morton, level);
        if (!c.valid())
            throw std::out_of_range("MortonTrieDDF::at_slot: key not found");
        return deref_slot(c, slot);
    }

    const double& at_slot(const D_morton& morton, uint8_t level, uint8_t slot) const
    {
        Cursor c = const_cast<MortonTrieDDF*>(this)->build_path_to_leaf(morton, level, false);
        if (!c.valid())
            throw std::out_of_range("MortonTrieDDF::at_slot: key not found");
        return c.leaf->data->f[c.idx_in_leaf][slot];
    }

    // Returns reference for slot, inserting (and zero-initialising the entire
    // 27-wide cell) if absent. Matches unordered_map::operator[] semantics for
    // each individual slot view.
    double& touch_slot(const D_morton& morton, uint8_t level, uint8_t slot)
    {
        Cursor c = build_path_to_leaf(morton, level, /*create=*/true);
        if (!exists(morton, level))
        {
            // First touch on this cell: zero the 27-wide row and mark present.
            std::memset(&c.leaf->data->f[c.idx_in_leaf][0], 0, sizeof(double) * Q_DIM);
#ifdef ABPattern
            std::memset(&c.leaf->data->f2[c.idx_in_leaf][0], 0, sizeof(double) * Q_DIM);
#endif
            set_existence_bit(morton);
        }
        return c.leaf->data->f[c.idx_in_leaf][slot];
    }

    // Bulk read/write of all C_Q directions for one cell (Collision hot path).
    void load_all(const D_morton& morton, uint8_t level, double* out)
    {
        Cursor c = cursor(morton, level);
        if (!c.valid())
            throw std::out_of_range("MortonTrieDDF::load_all: key not found");
        std::memcpy(out, &c.leaf->data->f[c.idx_in_leaf][0], sizeof(double) * Q_DIM);
    }

    void store_all(const D_morton& morton, uint8_t level, const double* in)
    {
        Cursor c = build_path_to_leaf(morton, level, /*create=*/true);
        if (!exists(morton, level))
            set_existence_bit(morton);
        std::memcpy(&c.leaf->data->f[c.idx_in_leaf][0], in, sizeof(double) * Q_DIM);
    }

private:
    uint8_t depth_;
    DDFTrieNode* root_;
    std::vector<DDFTrieNode*> search_path_stack_;
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

    Cursor build_path_to_leaf(const D_morton& morton, uint8_t level, bool create)
    {
        DDFTrieNode* node = root_;
        const uint8_t search_depth = depth_ - refine_level + level;

        for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
        {
            const uint8_t start_pos = (depth_ - i_depth) * DIM;
            const uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
            if (!node->children)
            {
                if (!create)
                    return Cursor{};
                node->children = new DDFTrieNode[BRANCH]();
                node->children[index].level = i_depth;
            }
            node = &node->children[index];
            search_path_stack_[i_depth] = node;
        }

        if (!node->data)
        {
            if (!create)
                return Cursor{};
            node->data = new LeafData();
        }

        stack_code_ = morton;
        Cursor c;
        c.leaf = node;
        c.idx_in_leaf = EXTRACT_MORTON_CONVERT_INDEX(morton, DIM);
        return c;
    }

    void free_subtree(DDFTrieNode* node)
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
            delete node->data;
            node->data = nullptr;
        }
    }
};

// MortonTrieDDFSlot: single-direction view over shared MortonTrieDDF,
// drop-in replacement for std::unordered_map<D_morton, double>.

class MortonTrieDDFSlot
{
public:
    using key_type = D_morton;
    using mapped_type = double;
    using value_type = std::pair<const D_morton, double>;
    using size_type = std::size_t;

    MortonTrieDDFSlot() = default;
    MortonTrieDDFSlot(MortonTrieDDF* parent, uint8_t slot)
        : parent_(parent), slot_(slot) {}

    void bind(MortonTrieDDF* parent, uint8_t slot)
    {
        parent_ = parent;
        slot_ = slot;
    }

    // Hot-path APIs (Stream / Collision inner loops)
    inline double& at(const D_morton& key)
    {
        return parent_->at_slot(key, refine_level, slot_);
    }

    inline const double& at(const D_morton& key) const
    {
        return parent_->at_slot(key, refine_level, slot_);
    }

    inline double& operator[](const D_morton& key)
    {
        return parent_->touch_slot(key, refine_level, slot_);
    }

    inline bool exists(const D_morton& key) const
    {
        return parent_->exists(key, refine_level);
    }

    inline bool exists(const D_morton& key, uint8_t level) const
    {
        return parent_->exists(key, level);
    }

    inline size_type count(const D_morton& key) const
    {
        return exists(key) ? 1 : 0;
    }

    // Init-time APIs
    // Note: existence is tracked per-cell (one bitmap shared by all C_Q slots
    // of the field), not per-slot. The natural LBM init pattern initialises
    // all C_Q slots for each cell, so this is benign. The first slot insert for
    // a cell zero-initialises the row and sets the bitmap; subsequent slot
    // inserts for the same cell still write their value (we cannot distinguish
    // "slot was explicitly set" from "slot is zero because cell was just
    // created").
    std::pair<bool, bool> insert(const std::pair<D_morton, double>& kv)
    {
        const bool cell_existed = parent_->exists(kv.first, refine_level);
        // touch_slot ensures the leaf exists and zero-inits the 27-wide row on
        // first touch of the cell.
        parent_->touch_slot(kv.first, refine_level, slot_) = kv.second;
        return {true, !cell_existed};
    }

    // Range insert: matches the std::unordered_map signature used by some
    // legacy / DEPRECATED LBM init paths. Loops over the iterator range.
    template <typename InputIt>
    void insert(InputIt first, InputIt last)
    {
        for (auto it = first; it != last; ++it)
            insert(std::make_pair(it->first, it->second));
    }

    // Shared-trie access
    MortonTrieDDF* parent() const { return parent_; }
    uint8_t slot() const { return slot_; }

private:
    MortonTrieDDF* parent_ = nullptr;
    uint8_t slot_ = 0;
};
