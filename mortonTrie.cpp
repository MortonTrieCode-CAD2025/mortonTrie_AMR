#include "header/mortonTrie.hpp"
#include "header/Morton_Assist.h"
#include <iostream>
#include <functional>
#include <unordered_map>

#define LOOP_METHOD_FINDFIRST1

void MortonTrie::insert(const D_morton morton, const uint8_t level, const double value) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            node->children = new TrieNode[BRANCH];
            node->children[index].level = i_depth;
        }
        
        node = &node->children[index];
    }
    // last level
        uint8_t start_pos = DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->data) {
            node->data = new double[BRANCH];
            #ifdef ABPattern
            node->data2 = new double[BRANCH];
            #endif
            // node->children[index].level = search_depth - 1;
        }
        node->data[index] = value;

    // Update existence bitmap
    set_existence_bit(morton, level);
}

double MortonTrie::search(const D_morton morton, const uint8_t level) const {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[search] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
    }

    return node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
}

void MortonTrie::init_searchStack(const D_morton morton,const uint8_t level) {
    TrieNode* node = root; search_path_stack->clear(); search_path_stack->resize(depth);

    uint8_t search_depth = depth - refine_level + level;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        node = &node->children[EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos)];
        search_path_stack->at(i_depth) = node;
    }

    stack_code = morton;
}

double MortonTrie::search_halfway(const D_morton morton, const uint8_t level) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;

    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        return search_path_stack->at(search_depth - 2)->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
    }

    node = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[Bottom-Up search] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
        search_path_stack->at(i_depth) = node;
    }

    stack_code = morton;

    return node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
    
}


double MortonTrie::search_halfway_noUpdateStack(const D_morton morton, const uint8_t level) const {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;

    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        return search_path_stack->at(search_depth - 2)->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
    }

    node = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[search_bottom_up] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
    }

    return node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
}

void MortonTrie::update_searchStack(const D_morton morton, const uint8_t level) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;
    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        return;
    }

    node = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        node =&node->children[EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos)];
        search_path_stack->at(i_depth) = node;
    }

    stack_code = morton;
}

void MortonTrie::update(const D_morton morton, const uint8_t level, const double value) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[Update] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
    }

    #ifdef AAPattern
    node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
    #endif
    #ifdef ABPattern
    node->data2[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
    #endif

    // Update existence bitmap (in case this is a new cell)
    set_existence_bit(morton, level);
}

void MortonTrie::update_halfway(const D_morton morton, const uint8_t level, const double value) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;

    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        #ifdef AAPattern
        search_path_stack->at(search_depth - 2)->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
        #endif
        #ifdef ABPattern
        search_path_stack->at(search_depth - 2)->data2[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
        #endif
    }

    node = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[Bottom-Up search] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
        search_path_stack->at(i_depth) = node;
    }

    stack_code = morton;

    #ifdef AAPattern
    node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
    #endif
    #ifdef ABPattern
    node->data2[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
    #endif

}

// Bitmap-based O(1) existence query implementation moved to header file for inlining

// Background index extraction is no longer needed with Morton-ordered bitmaps
// Morton code is used directly as array index

uint32_t MortonTrie::extract_refinement_index(const D_morton& morton, const uint8_t level) const {
    // Extract refinement field from lower bits
    // The refinement field is the lower (level * DIM) bits of the Morton code
    uint32_t bits_to_extract = level * DIM;

    // Convert bitset to uint64_t and mask to get lower bits
    uint64_t morton_val = morton.to_ullong();
    uint64_t mask = (1ULL << bits_to_extract) - 1;

    return static_cast<uint32_t>(morton_val & mask);
}

// set_existence_bit() moved to header file for inlining

// Morton-ordered bitmap implementation uses direct indexing

// Tree traversal-based O(depth) existence query (for ablation study)
// This method checks existence by traversing the trie tree structure
// Time complexity: O(depth), where depth = search_depth - 1
bool MortonTrie::exists_tree_traversal(const D_morton& morton, const uint8_t level) const {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;

    // Traverse the tree to search_depth - 1
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);

        // If children don't exist at this level, the node doesn't exist
        if (!node->children) {
            return false;
        }

        node = &node->children[index];
    }

    // Check if data exists at the final level
    if (!node->data) {
        return false;
    }

    // Check if the specific index has been initialized
    // Note: This is a simplified check. In practice, you might need additional
    // bookkeeping to track which indices in the data array are actually in use
    uint8_t final_index = EXTRACT_MORTON_CONVERT_INDEX(morton, DIM);

    // For this implementation, we assume if data array exists, the element exists
    // This matches the behavior of the insert/search methods
    return true;
}

// Morton-ordered bitmap array requires no additional initialization

// Get bitmap memory usage in bytes
size_t MortonTrie::get_bitmap_memory_usage() const {
    // Each CellBitmap is a std::bitset<BITMAP_BITS_PER_CELL>
    // Morton-ordered array: (MAX_BG_MORTON + 1) positions, includes empty slots
    // Memory = (MAX_BG_MORTON + 1) * sizeof(CellBitmap)
    size_t bytes = existence_bitmaps.size() * sizeof(CellBitmap);

    // Add vector overhead
    bytes += sizeof(std::vector<CellBitmap>);

    return bytes;
}

// Morton-ordered bitmap memory usage

// Get tree structure memory usage in bytes (approximate)
size_t MortonTrie::get_tree_memory_usage() const {
    size_t total_memory = sizeof(MortonTrie);

    // Count root node
    total_memory += sizeof(TrieNode);

    // Recursively count allocated nodes only
    std::function<void(TrieNode*)> count_nodes = [&](TrieNode* node) {
        if (!node) return;

        // Count children array if allocated
        if (node->children) {
            total_memory += sizeof(TrieNode) * BRANCH;
            // Recursively count each child
            for (uint8_t i = 0; i < BRANCH; ++i) {
                count_nodes(&node->children[i]);
            }
        }

        // Count data array if allocated
        if (node->data) {
            total_memory += sizeof(double) * BRANCH;
            #ifdef ABPattern
            total_memory += sizeof(double) * BRANCH;
            #endif
        }
    };

    count_nodes(root);

    // Add search path stack
    total_memory += sizeof(std::vector<TrieNode*>) + sizeof(TrieNode*) * depth;

    return total_memory;
}

