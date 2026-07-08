#pragma once

#include "settings.hpp"
#include "General.h"
#include "mortonTrieGeneric.hpp"  // per-instance, templated trie used by MortonTrieMap<T>
#include <array>
#include <vector>
#include <utility>
#include <iterator>
#include <cmath>
#include <bitset>
#include <unordered_map>

class TrieNode {
public:
    TrieNode* children; //in Z-order: 000, 001, 010, 011, 100, 101, 110, 111
    uint8_t level;
    double* data;
    #ifdef ABPattern
    double* data2;
    #endif

    TrieNode() : children(nullptr), level(0), data(nullptr) {
        #ifdef ABPattern
        data2 = nullptr;
        #endif
    }
};

using NgbrList = std::array<TrieNode*, STENCIL_SIZE-1>;

class MortonTrie {
public:

    void insert(const D_morton morton, const uint8_t level, const double value);
    double search(const D_morton morton, const uint8_t level) const; //Top-down search
    double search_halfway(const D_morton morton, const uint8_t level); //Bottom-up search
    double search_halfway_noUpdateStack(const D_morton morton, const uint8_t level) const; //Bottom-up search without updating the search path stack
    void init_searchStack(D_morton morton, uint8_t level);
    void update(const D_morton morton, const uint8_t level, const double value);
    void update_halfway(const D_morton morton, const uint8_t level, const double value);
    void update_searchStack(const D_morton morton, const uint8_t level);

    // Bitmap-based O(1) existence query using Morton-ordered array
    // OPTIMIZED: Inline implementation for maximum performance
    inline bool exists(const D_morton& morton, const uint8_t level) const {
        // Single bitset to uint64_t conversion (avoid redundant conversions)
        const uint64_t morton_val = morton.to_ullong();

        // Compile-time constants for bit manipulation
        constexpr uint32_t bg_bits_start = refine_level * DIM;
        constexpr uint64_t bg_mask = (1ULL << (BG_DEPTH * DIM)) - 1;

        // Extract background Morton code from high-order bits
        const uint64_t bg_morton = (morton_val >> bg_bits_start) & bg_mask;

        // Extract refinement index from low-order bits (inline, no function call)
        const uint32_t bits_to_extract = level * DIM;
        const uint64_t ref_mask = (1ULL << bits_to_extract) - 1;
        const uint32_t ref_index = static_cast<uint32_t>(morton_val & ref_mask);

        // Bounds check with branch prediction hint (expect valid indices)
        if (__builtin_expect(bg_morton > MAX_BG_MORTON || ref_index >= BITMAP_BITS_PER_CELL, 0)) {
            return false;
        }

        // Direct bitmap access - single memory lookup
        return existence_bitmaps[bg_morton][ref_index];
    }

    // Tree traversal-based O(depth) existence query (for ablation study baseline)
    bool exists_tree_traversal(const D_morton& morton, const uint8_t level) const;

    // Memory measurement utilities
    size_t get_bitmap_memory_usage() const;
    size_t get_tree_memory_usage() const;

    uint8_t get_depth() const { return depth; }
    TrieNode* get_root() const { return root; }
    
    static MortonTrie& getInstance() {
        static MortonTrie instance;
        return instance;
    };


private:

    TrieNode* root;
    uint8_t depth;

    // Bitmap-based existence query data structures
    static constexpr uint32_t BITMAP_BITS_PER_CELL = 1 << (refine_level * DIM); // 2^(L*DIM) bits per background cell
    static constexpr uint32_t BACKGROUND_GRID_SIZE = Nx * Ny * Nz; // Total background grid cells
    using CellBitmap = std::bitset<BITMAP_BITS_PER_CELL>;

    // Morton-ordered bitmap array for O(1) existence queries:
    // Bitmaps are indexed directly by background Morton code, enabling single-access lookup
    // Maximum background Morton code for grid (Nx-1, Ny-1, Nz-1) = (20, 9, 5) is 5446
    static constexpr uint32_t MAX_BG_MORTON = 5446; // Pre-computed for 21×10×6 grid
    static constexpr uint32_t BG_DEPTH = 6; // Background depth: depth - refine_level

    std::vector<CellBitmap> existence_bitmaps; // Indexed by background Morton code

    MortonTrie() {
        depth = static_cast<uint8_t>(log2(max_of_three(Nx, Ny, Nz)) + 1 + refine_level);
        // not containing the single root level , +1 for std::ceil(log2(max_of_three(Nx, Ny, Nz)))
        root = new TrieNode();
        search_path_stack = new std::vector<TrieNode*>(depth);

        // Initialize Morton-ordered bitmap array
        existence_bitmaps.resize(MAX_BG_MORTON + 1);
    }

    ~MortonTrie() = default;

    // Helper methods for bitmap operations
    uint32_t extract_refinement_index(const D_morton& morton, const uint8_t level) const;

    // OPTIMIZED: Inline implementation for maximum performance
    inline void set_existence_bit(const D_morton& morton, const uint8_t level) {
        // Single bitset to uint64_t conversion
        const uint64_t morton_val = morton.to_ullong();

        // Compile-time constants for bit manipulation
        constexpr uint32_t bg_bits_start = refine_level * DIM;
        constexpr uint64_t bg_mask = (1ULL << (BG_DEPTH * DIM)) - 1;

        // Extract background Morton code
        const uint64_t bg_morton = (morton_val >> bg_bits_start) & bg_mask;

        // Extract refinement index (inline, no function call)
        const uint32_t bits_to_extract = level * DIM;
        const uint64_t ref_mask = (1ULL << bits_to_extract) - 1;
        const uint32_t ref_index = static_cast<uint32_t>(morton_val & ref_mask);

        // Bounds check with branch prediction hint (expect valid indices)
        if (__builtin_expect(bg_morton <= MAX_BG_MORTON && ref_index < BITMAP_BITS_PER_CELL, 1)) {
            existence_bitmaps[bg_morton][ref_index] = true;
        }
    }

    std::vector<TrieNode*> *search_path_stack;
    D_morton stack_code = 0;
};

// Forward declaration for MortonTrieMap
template<typename Value>
class MortonTrieMap;

// Iterator class for MortonTrieMap.
//
// The iterator can be in one of three modes:
//   * END     — past-the-end sentinel. Returned by end(). Does NOT trigger
//               cache rebuild, so `m.find(k) == m.end()` is O(1).
//   * CACHE   — backed by cached_data_ vector iterators. Returned by begin()
//               for forward iteration. Triggers update_cache() once at begin().
//   * CURSOR  — backed by a MortonTrieGeneric::Cursor + key. Returned by
//               find() on a hit. No cache work at all; the value is fetched
//               lazily from the trie on first dereference.
//
// This three-mode design is what makes find()/end() cheap: previously both
// called update_cache() (a 22M-bitmap scan), which dominated the LBM hot path
// because Grid_Manager / Lat_Voxelize interleave find() + insert() at high
// frequency and each insert invalidated the cache.
template<typename Value>
class MortonTrieIterator {
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = std::pair<D_morton, Value>;  // Remove const from key for simplicity
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

private:
    enum class Mode : uint8_t { END, CACHE, CURSOR };
    using Cursor = typename MortonTrieGeneric<Value>::Cursor;

    MortonTrieMap<Value>* map_ = nullptr;
    Mode mode_ = Mode::END;

    // CACHE mode: half-open range into cached_data_.
    typename std::vector<std::pair<D_morton, Value>>::iterator iter_;
    typename std::vector<std::pair<D_morton, Value>>::iterator iter_end_;

    // CURSOR mode: trie leaf handle + key. mat_ is materialised on first
    // dereference; subsequent dereferences of the same iterator reuse it.
    Cursor cursor_;
    D_morton key_;
    mutable bool mat_valid_ = false;
    mutable std::pair<D_morton, Value> mat_;

    void ensure_mat() const {
        if (!mat_valid_) {
            mat_.first = key_;
            mat_.second = map_->deref(cursor_);
            mat_valid_ = true;
        }
    }

public:
    MortonTrieIterator() = default;

    // END sentinel constructor.
    MortonTrieIterator(MortonTrieMap<Value>* map, std::nullptr_t)
        : map_(map), mode_(Mode::END) {}

    // CACHE iterator constructor (used by begin()).
    MortonTrieIterator(MortonTrieMap<Value>* map,
                       typename std::vector<std::pair<D_morton, Value>>::iterator iter,
                       typename std::vector<std::pair<D_morton, Value>>::iterator iter_end)
        : map_(map), mode_(Mode::CACHE), iter_(iter), iter_end_(iter_end) {}

    // CACHE iterator constructor (const-map variant, mirrors original API).
    MortonTrieIterator(const MortonTrieMap<Value>* map,
                       typename std::vector<std::pair<D_morton, Value>>::const_iterator iter,
                       typename std::vector<std::pair<D_morton, Value>>::const_iterator iter_end)
        : map_(const_cast<MortonTrieMap<Value>*>(map)),
          mode_(Mode::CACHE),
          iter_(const_cast<typename std::vector<std::pair<D_morton, Value>>::iterator>(iter)),
          iter_end_(const_cast<typename std::vector<std::pair<D_morton, Value>>::iterator>(iter_end)) {}

    // CURSOR iterator constructor (used by find()).
    MortonTrieIterator(MortonTrieMap<Value>* map, Cursor c, const D_morton& key)
        : map_(map), mode_(Mode::CURSOR), cursor_(c), key_(key) {}

    reference operator*() const {
        if (mode_ == Mode::CACHE)
            return *iter_;
        ensure_mat();
        return mat_;
    }
    pointer operator->() const {
        if (mode_ == Mode::CACHE)
            return &(*iter_);
        ensure_mat();
        return &mat_;
    }

    MortonTrieIterator& operator++() {
        if (mode_ == Mode::CACHE) {
            ++iter_;
            if (iter_ == iter_end_)
                mode_ = Mode::END;
        } else if (mode_ == Mode::CURSOR) {
            // find() returns a single-element iterator; advancing moves to end.
            mode_ = Mode::END;
        }
        return *this;
    }
    MortonTrieIterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }

    bool operator==(const MortonTrieIterator& other) const {
        // A CACHE iterator that has reached iter_end_ is logically END, so treat
        // it as equal to an END sentinel. This covers the empty-cache case
        // (begin()==end()==nullptr) even if a caller somehow builds a CACHE
        // iterator from an empty range.
        const bool this_end = (mode_ == Mode::END) ||
                              (mode_ == Mode::CACHE && iter_ == iter_end_);
        const bool other_end = (other.mode_ == Mode::END) ||
                               (other.mode_ == Mode::CACHE && other.iter_ == other.iter_end_);
        if (this_end && other_end)
            return true;
        if (this_end != other_end)
            return false;
        if (mode_ == Mode::CACHE)
            return iter_ == other.iter_;
        return key_ == other.key_;  // CURSOR
    }
    bool operator!=(const MortonTrieIterator& other) const { return !(*this == other); }
};

/**
 * @brief STL-compatible wrapper for MortonTrie that provides map-like interface
 * This class allows MortonTrie to be used as a drop-in replacement for std::unordered_map and std::map
 */
template<typename Value>
class MortonTrieMap {
public:
    using key_type = D_morton;
    using mapped_type = Value;
    using value_type = std::pair<const D_morton, Value>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using iterator = MortonTrieIterator<Value>;
    using const_iterator = MortonTrieIterator<Value>;

private:
    // Per-instance trie storage (replaces the legacy MortonTrie singleton). Each
    // MortonTrieMap<Value> owns its own trie so different value types do not
    // collide on a shared singleton.
    mutable MortonTrieGeneric<Value> trie_;
    uint8_t default_level_;
    // Lazily-rebuilt Morton-ordered snapshot used only for iteration (begin/end).
    // find/at/operator[]/exists/insert/count all go through trie_ directly.
    mutable std::vector<std::pair<D_morton, Value>> cached_data_;
    mutable bool cache_valid_;

    void invalidate_cache() { cache_valid_ = false; }
    void update_cache() const;

public:
    MortonTrieMap(uint8_t level = refine_level)
        : trie_(), default_level_(level), cache_valid_(false) {}

    // O(1) swap. Backs the STL idiom `temp.swap(map.at(k))` used in
    // Lat_Voxelize.cpp to clear a level's overlap map by trading guts with an
    // empty temp. The trie_, iteration cache, and cache_valid_ flag are all
    // exchanged together so neither side sees a stale cache.
    void swap(MortonTrieMap& other) noexcept
    {
        using std::swap;
        trie_.swap(other.trie_);
        swap(default_level_, other.default_level_);
        cached_data_.swap(other.cached_data_);
        swap(cache_valid_, other.cache_valid_);
    }

    // STL container interface
    iterator begin() const;
    iterator end() const;
    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const { return end(); }

    bool empty() const;
    size_type size() const;
    void clear();

    // Map interface
    // Single insert overload: value_type is std::pair<const D_morton, Value>,
    // which binds cleanly to both brace-init {k, v} and std::make_pair(k, v)
    // (the latter converts via the implicit const-key conversion). The earlier
    // second overload `insert(const std::pair<D_morton, Value>&)` made
    // map.insert({k, v}) ambiguous under C_MAP_TYPE=3.
    std::pair<iterator, bool> insert(const value_type& value);
    size_type erase(const D_morton& key);
    template<typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args);

    iterator find(const D_morton& key) const;
    size_type count(const D_morton& key) const;

    Value& operator[](const D_morton& key);
    Value& at(const D_morton& key);
    const Value& at(const D_morton& key) const;

    // Bitmap-based O(1) existence query
    bool exists(const D_morton& key) const {
        return trie_.exists(key, default_level_);
    }

    bool exists(const D_morton& key, uint8_t level) const {
        return trie_.exists(key, level);
    }

    // ---- Trie-native hot-path APIs (P1-2 neighbour direct query) ------------
    // Returns a cursor (leaf handle + index) so callers can amortise repeated
    // accesses on the same cell across e.g. 27 DDF directions.
    using Cursor = typename MortonTrieGeneric<Value>::Cursor;

    Cursor cursor(const D_morton& key) {
        return trie_.cursor(key, default_level_);
    }

    Value& deref(const Cursor& c) { return trie_.deref(c); }
    const Value& deref(const Cursor& c) const { return trie_.deref(c); }

    // Stream-pattern neighbour lookup: given the current node code and a
    // direction (ex, ey, ez each in {-1, 0, +1}), descend to the neighbour
    // leaf reusing the path stack from the last access (search_halfway).
    Value& at_neighbor(const D_morton& code, int ex, int ey, int ez) {
        Cursor c = trie_.cursor_neighbor(code, default_level_, ex, ey, ez);
        if (!c.valid())
            throw std::out_of_range("MortonTrieMap::at_neighbor: neighbour not found");
        return trie_.deref(c);
    }

    // Range insertion for compatibility with existing code
    template<typename InputIt>
    void insert(InputIt first, InputIt last);
};