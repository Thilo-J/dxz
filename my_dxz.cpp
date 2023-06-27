#include <iostream>
extern "C"{
#include "quad_linked_list.h"
}
extern "C" {
#include "sparse_matrix.h"
}
#include <map>
#include <tuple>
#include <list>
#include <fstream>
#include <string>
#include <sstream>
#include <windows.h>
#include "my_dxz.h"
#include <vector>
#include <unordered_map>
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

int node_counter = 1;

struct VectorHasher {
    int operator()(const std::vector<int>& V) const {
        int hash = V.size();
        for (auto& i : V) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

struct Node {
    int value;
    const Node* left;
    const Node* right;
};

const Node* top = new Node({ -1, NULL, NULL });
const Node* bottom = new Node({ -2, NULL, NULL });

struct dxz_manager {
    const Node* root = bottom;
    std::unordered_map<std::vector<int>, const Node*, VectorHasher> cache;
};

const Node* search(list sparse_matrix, int k, int max, dxz_manager* manager) {
    list col, row, next;

    // Base case
    if (get_right(sparse_matrix) == sparse_matrix) {
        return top;
    }

    //Check cache
    std::vector<int> cache_key;
    list key_col = sparse_matrix;
    while ((key_col = get_right(key_col)) != sparse_matrix) {
        cache_key.push_back((int)key_col);
    }
    auto iter = manager->cache.find(cache_key);
    if (iter != manager->cache.end()) {
        return iter->second;
    }

    // Main algorithm:
    col = choose_column_with_min_data(sparse_matrix, max);
    if (get_data(col)->data == 0) return bottom;

    cover_column(col);
    const Node* x = bottom;
    for (row = col; (row = get_down(row)) != col; ) {
        int row_index = get_data(row)->data;  // save the row number
        for (next = row; (next = get_right(next)) != row; )
            cover_column(get_data(next)->list_data);
        
        const Node* y = search(sparse_matrix, k+1, max, manager);        
        for (next = row; (next = get_left(next)) != row; )
            uncover_column(get_data(next)->list_data);
        if (y->value != -2) {
            //x = get_node(row_index, x, y, manager);
            x = new Node({
                row_index,
                x,
                y
            });
            node_counter++;
            manager->root = x;            
        }
    }
    uncover_column(col);

    manager->cache[cache_key] = x; // Save in cache
    return x;
    
}

const Node* dxz_get_exact_covers(int rows, int cols, int matrix[]) {
    dxz_manager manager;
    list sparse_matrix;
    sparse_matrix = create_sparse(rows, cols, matrix);
    free(matrix);
    search(sparse_matrix, 0, rows, &manager);
    return manager.root;
}

void preorder_traversal(const Node* node, const Node* prev, bool add, std::vector<int>* solution, std::list<std::vector<int>>* solutions)
{
    if (add)
        solution->push_back(prev->value);
    if (node->value == -2)
        return;
    if (node->value == -1)
        return solutions->push_back(*solution);
    preorder_traversal(node->left, node, false, solution, solutions);
    preorder_traversal(node->right, node, true, solution, solutions);
    solution->pop_back();
}

std::list<std::vector<int>> get_solutions(const Node* root) {
    std::list<std::vector<int>> solutions;
    std::vector<int> solution;
    if (root == bottom) return solutions;
    preorder_traversal(root->left, root, false, &solution, &solutions);
    preorder_traversal(root->right, root, true, &solution, &solutions);
    return solutions;
}


std::list<std::vector<int>> solve_matrix(int rows, int cols, int matrix[])
{
    const Node* node = dxz_get_exact_covers(rows, cols, matrix);
    return get_solutions(node);
}

PYBIND11_MODULE(dxz, module_handle)
{
    module_handle.doc() = "C++ dxz solver";
    module_handle.def("dxz_solve", &solve_matrix);
}