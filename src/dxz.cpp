extern "C"{
#include "quad_linked_list.h"
}
extern "C" {
#include "sparse_matrix.h"
}
#include <map>
#include <tuple>
#include <list>
#include <sstream>
#include <vector>
#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

int node_counter = 1;

struct VectorHasher {
    int operator()(const std::vector<size_t>& V) const {
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

const Node* top = new Node({-1, NULL, NULL});
const Node* bottom = new Node({-2, NULL, NULL});

const Node* get_node(int value, const Node* left, const Node* right)
{
    if (right->value == bottom->value)
        return left;
    
    return  new Node({
            value,
            left,
            right
    });
    
}


struct dxz_manager {
    const Node* root = bottom;
    std::unordered_map<std::vector<size_t>, const Node*, VectorHasher> cache;

};

const Node* search(list sparse_matrix, int k, int max, dxz_manager* manager) {
    list col, row, next;

    // Base case
    if (get_right(sparse_matrix) == sparse_matrix) {
        return top;
    }

    //Check cache
    std::vector<size_t> cache_key;
    list key_col = sparse_matrix;
    while ((key_col = get_right(key_col)) != sparse_matrix) {
        cache_key.push_back((size_t)key_col);
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
        if (y->value != bottom->value) {
            x = get_node(row_index, x, y);
            node_counter++;
            manager->root = x;            
        }
    }
    uncover_column(col);

    manager->cache[cache_key] = x; // Save in cache
    return x;
    
}


void preorder_traversal(const Node* node, std::vector<int>* solution, std::list<std::vector<int>>* solutions, int max_solutions)
{
    if(max_solutions != 0 && solutions->size() >= max_solutions)
        return;
    if(node->value == bottom->value)
        return;
    if(node->value == top->value)
        return solutions->push_back(*solution);
    preorder_traversal(node->left, solution, solutions, max_solutions);
    solution->push_back(node->value);
    preorder_traversal(node->right, solution, solutions, max_solutions);
    solution->pop_back();
}


std::list<std::vector<int>> zdd_to_list(const Node* root, int max_solutions)
{
    std::list<std::vector<int>> solutions;
    std::vector<int> solution;
    preorder_traversal(root, &solution, &solutions, max_solutions);
    return solutions;
}

std::list<std::vector<int>>  dxz_get_exact_covers(int rows, int cols, std::list<int> matrix, int max_solutions) {
    int arr[matrix.size()];
    std::copy(matrix.begin(), matrix.end(), arr);
    dxz_manager manager;
    list sparse_matrix;
    sparse_matrix = create_sparse(rows, cols, arr);
    search(sparse_matrix, 0, rows, &manager);
    std::list<std::vector<int>>  solutions = zdd_to_list(manager.root, max_solutions);
    return solutions;
}

int node_count(const Node* node)
{
    if (node->value == bottom->value) return 0; 
    if (node->value == top->value)  return 1;
    return node_count(node->left) + node_count(node->right);
}


int dxz_get_number_of_solutions(int rows, int cols, std::list<int> matrix) {
    int arr[matrix.size()];
    std::copy(matrix.begin(), matrix.end(), arr);
    dxz_manager manager;
    list sparse_matrix;
    sparse_matrix = create_sparse(rows, cols, arr);
    search(sparse_matrix, 0, rows, &manager);
    return node_count(manager.root);
}

PYBIND11_MODULE(dxz, module_handle)
{
    module_handle.doc() = "Exact cover solver that returns all solutions";
    module_handle.def("dxz_solve", &dxz_get_exact_covers, 
    "Returns all solutions to a given exact cover problem.\n\
    Paramaters:\n\
        arg0: number of rows\n\
        arg1: number of columns\n\
        arg2: matrix in a single list\n\
        arg3: max number of solutions to be returned (0 for all)"
    );
    module_handle.def("dxz_count", &dxz_get_number_of_solutions, 
    "Returns the number of solutions of a given exact cover problem.\n\
    Paramaters:\n\
        arg0: number of rows\n\
        arg1: number of columns\n\
        arg2: matrix in a single list"
    );
}

