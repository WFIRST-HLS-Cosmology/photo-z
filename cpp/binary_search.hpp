/*! \file
 *
 * Contains a binary search class, designed to quickly search a sorted array.
 *
 * The BinaryTree class searches a sorted array of type T objects and returns the index of the value
 * that most closely matches the value for which the user is searching. If no value in the array is
 * within a specified threshold, the search fails catastrophically (by terminating the program).
 *
 * \todo generalize so that types other than float, double, long double can be used as value type
 */

#ifndef BINARY_SEARCH_HPP
#define BINARY_SEARCH_HPP

#ifdef HAS_FOLLY
    #include <folly/FBVector.h>
    #define VECTOR folly::fbvector
#else
    #include <vector>
    #define VECTOR std::vector
#endif

#include <iostream>
#include <cmath>

/*! A single node in the BinaryTree class.
 *
 */
template<typename T>
class TreeNode
{
public:

    /*!
     * \brief TreeNode creates a node.
     * \param value is the value stored at the node.
     * \param index is the index associated with the value stored at the node.
     */
    TreeNode(const T value, const uint index)
    : value_(value), index_(index)
    {}

    /*!
    * \brief Next identifies the next node in tree search. The goal is to eventually find a node
    * whose value is approximately equal to the target.
    * \param target is the value for which we are searching.
    * \return returns a pointer to the next node in the search process. If the current node is the
    * node of interest (the approximate match), then nullptr is returned.
    */
    TreeNode* Next(const T& target) const
    {
        const T diff = target - value_;

        if (diff * diff < 2.5e-7) // The match is sufficient (value_ is close enough to target)
        {
            return nullptr;
        }
        else if (target > value_) // target is larger than value_ of this node
        {
             return greater_;
        }
        return less_; // target is smaller than the value_ of this node
    }

    /*!
     * \brief Insert inserts a new node into tree (i.e., it sets the less_ or greater_ pointer).
     * \param child_node is the node that will be insterted.
     * \return Returns nullptr if the node was successfully insterted as a child of the current node.
     * Returns the address of this node if the new node is approximately equal to this node.
     * Otherwise, it returns a pointer to the next candidate parent node (the user needs to step
     * forward and call Insert on the result unlit either nullptr or this* is reached).
     */
    TreeNode<T>* Insert(TreeNode<T>* child_node)
    {
        const T diff = child_node->value() - value_;

        if (diff * diff < 2.5e-7) /// \todo generalize this so that the constant is an input parameter.
        {
            return this;
        }
        else if (child_node->value() > value_)
        {
            if (greater_ == nullptr)
            {
                greater_ = child_node;
                return nullptr; // success!
            }
            else { return greater_; }
        }

        // at this point, child_node.value() < value_

        if (less_ == nullptr)
        {
            less_ = child_node;
            return nullptr; // success!
        }
        else { return less_; }

        return nullptr;
    }

    /*!
     * IsTerminal returns true if this is a terminal node; false otherwise
     * */
    const bool IsTerminal() const { return (less_ == nullptr && greater_ == nullptr); }

    const uint index() const { return index_; }

    const T value() const { return value_; }

private:

    T      value_;
    uint   index_;
    TreeNode<T>* less_ = nullptr;
    TreeNode<T>* greater_ = nullptr;
};



/*!  Implements a binary search tree.
 * */
template<typename T>
class BinaryTree
{
public:

    /*!
     * \brief BinaryTree constructs a search tree from the data provided in data_vector.
     * \param data_vector is a vector of objects that will be repeatedly searched (It makes little
     * sense to use a search tree if the search will only happen once or twice.) This vector should
     * sorted by value before the tree is constructed, for best results.
     */
    BinaryTree(VECTOR<T>& data_vector) : data_vector_(data_vector)
    {
        insert_nodes( optimal_node_order(data_vector_.size()) );
    }

    /*!
     * \brief index_of finds the index in data_vector of the closest match to the value parameter.
     * \param value is the value whose index we are searching for.
     * \return Returns the index data_vector of the closest match to the value parameter.
     */
    uint index_of(const T& value)
    {
        TreeNode<T>* next_node = &nodes_[0]; // the root node

        while ( next_node->Next(value) != nullptr ) { next_node = next_node->Next(value); }

        /// \todo generalize this so that the constant is an input parameter.
        if (std::abs(next_node->value() - value) > 0.0005)
        {
            std::cerr << "Error: Value is not in the tree.\n";
            exit(2);
        }

        return next_node->index();
    }

    private:

    /*!
     * \brief insert_nodes inserts the nodes into the tree.
     * \param node_order is the order in which the nodes will be added. This should be computed by
     * the optimal_node_order method.
     */
    void insert_nodes(VECTOR<uint> node_order)
    {
        // make sure that the vector has enough storage to hold the entire tree to begin with:

        nodes_.reserve(node_order.size());

        // add the root node
        const uint root_idx = node_order[0];
        const T root_val = data_vector_[root_idx];
        nodes_.emplace_back(root_val, root_idx);

        TreeNode<T>& root_node = nodes_[0];

        // add all of the other nodes, discarding duplicates.

        for (uint i = 1; i < node_order.size(); ++i)
        {
            const uint idx = node_order[i];
            const T value = data_vector_[idx];
            nodes_.emplace_back(value, idx);

            TreeNode<T>* result = root_node.Insert(&nodes_.back());

            while (result != nullptr) // if result == nullptr, success! Otherwise, walk the tree
            {
                if (result == result->Insert(&nodes_.back())) // if node with this value already exists
                { /// \todo remove the duplicate element from the list of nodes (don't pop_back() )
                    break;
                }
                result = result->Insert(&nodes_.back());
            }
        }
    }

    /*!
     * \brief optimal_node_order computes the approximate optimal order in which to add the nodes so
     * that the tree is as shallow as possible.
     * \param size_of_data_vector is simply the length of the data_vector being searched.
     * \return
     */
    VECTOR<uint> optimal_node_order(uint size_of_data_vector)
    {
        VECTOR<uint> ordering;

        // for very small inputs, we just go in linear order. The optimal ordering function fails
        // when the input is too small.
        if (size_of_data_vector < 4)
        {
            for (uint i = 0; i < size_of_data_vector; ++i)
            {
                ordering.push_back(i);
            }

            return ordering;
        }

        // compute the depth of the tree
        // depth = log2(size + 1) - 1, rounded up to the nearest integer, or...

        const int depth = floor(log2(size_of_data_vector + 1));

        float portion = size_of_data_vector / 2;

        ordering.push_back(size_of_data_vector / 2);

        for (int level = 0; level < depth; ++level)
        {
            portion *= 0.5; // N/4, N/8, N/16, ..., N/2^(level + 1)

            for (uint i = 1, v = (uint) floor(i * portion); v < size_of_data_vector; i += 2)
            {
                if (ordering.size() > 0 && v != ordering.back()) { ordering.push_back(v); }

                v = (uint) floor(i * portion);
            }
        }

        return ordering;
    }

    VECTOR<TreeNode<T>> nodes_;
    VECTOR<T>& data_vector_;
};


#endif // BINARY_SEARCH_HPP
