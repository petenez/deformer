/* 
 * File:   Node.h
 * Author: Petri Hirvonen
 * 8 April 2022
 */

#ifndef NODE_H
#define NODE_H

#include <memory>
#include <vector>
#include <map>
#include <set>

#include "Vector.h"
 
// Nodes in a square grid connected to their neighbors with springs.
// The nodes also form a complete quadtree for more efficient, multiscale solutions.
class Node final : public std::enable_shared_from_this<Node> {
    
private:
  
  // Selection of nodes (i, j, 0) indicated by their linear indices n = size*j + i
  typedef std::set<int> selection;

  /*
   * Types of nodes:
   * FIXED - nodes whose positions are fixed and cannot be relaxed
   * ACTIVE - nodes whose positions are relaxed
   * SAMPLED - nodes whose positions are up- or downsampled
   * TOP - nodes above the topmost sampled nodes in the quadtree
   * BOTTOM - nodes below the bottommost sampled nodes in the quadtree
   */
  enum class Types {
    FIXED,
    ACTIVE,
    SAMPLED,
    TOP,
    BOTTOM
  };

  // Map for printing nodes' types
  std::map<enum Types, std::string> Type_labels = {
    {Types::FIXED, "FIXED"},
    {Types::ACTIVE, "ACTIVE"},
    {Types::SAMPLED, "SAMPLED"},
    {Types::TOP, "TOP"},
    {Types::BOTTOM, "BOTTOM"}
  };
    
  // Type of a node
  enum Types m_type = Types::TOP;
  // Quadtree i index of a node
  const int m_i;
  // Quadtree j index of a node
  const int m_j;
  // Quadtree k index of a node
  const int m_k;
  // Child index of a node (m_cth child of its parent)
  const int m_c;
  // Equilibrium distance from its neighbors of a node.
  // A value of zero indicates that there is no equilibrium distance.
  double m_scale {1.0};
  // Stiffness of the springs to its neighbors of a node
  double m_stiffness {1.0};
  // Position of a node
  Vector m_position;
  // Update in position of a node
  Vector m_update;
  // Parent of a node
  std::weak_ptr<Node> m_parent;
  // Children of a node.
  // Leaf nodes have zero children, all others have four.
  std::vector<std::shared_ptr<Node>> m_children;
  // Neighbors of a node.
  // All nodes have four neighbors, for edge nodes some of them are nullptr.
  std::vector<std::weak_ptr<Node>> m_neighbors;
  
  /*
   * Populates the quadtree with nodes.
   */
  void populate();
  
  /*
   * Checks if node (p_i, p_j, p_k) is a descendant of the current node.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_k - quadtree k index
   * Returns:
   *   A bool indicating whether the node is a descendant of the current node.
   */
  bool is_descendant(const int p_i, const int p_j, const int p_k) const;
  
  /*
   * Returns the root node of the quadtree.
   * Returns:
   *   A shared pointer to the root node of the quadtree.
   */
  std::shared_ptr<Node> get_root();
  
  /*
   * Returns the node (p_i, p_j, p_k).
   * Note:
   *   This function is limited to searching the current node's descendants.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_k - quadtree k index
   * Returns:
   *   A shared pointer to the node indicated by the quadtree indices.
   *   A nullptr if it is not found.
   */
  std::shared_ptr<Node> get_node_recursive(const int p_i, const int p_j, const int p_k);
  
  /*
   * Returns the node (p_i, p_j, p_k).
   * Note:
   *   This function is not limited to searching the current node's descendants.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_k - quadtree k index
   * Returns:
   *   A shared pointer to the node indicated by the quadtree indices.
   *   A nullptr if it is not found.
   */
  std::shared_ptr<Node> get_node(const int p_i, const int p_j, const int p_k);
  
  /*
   * Sets the nodes' neighbors in a complete quadtree.
   */
  void set_neighbors();
  
  /*
   * Updates the active nodes from around the current node in the quadtree.
   * Note:
   *   See the documentation for more details.
   */
  void set_active();
  
  /*
   * Sets the bottom nodes in the quadtree.
   * Note:
   *   See the documentation for more details.
   */
  void set_bottom();
  
  /*
   * Sets the sampled nodes in the quadtree.
   * Note:
   *   See the documentation for more details.
   */
  void set_sampled();
  
  /*
   * Checks if the current node is at the edge of the selection.
   * Note:
   *   The current node must be part of the selection.
   * Parameters:
   *   p_nodes - a selection of nodes
   * Returns:
   *   A bool indicating whether the current node is at the edge of the selection.
   */
  bool is_edge(const selection& p_nodes);
  
  /*
   * Sets the equilibrium length scale of node (p_i, p_j, 0).
   * Note:
   *   Also updates the equilibrium length scales of the ancestor nodes.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_scale - length scale
   */
  void set_scale(const int p_i, const int p_j, const double p_scale);
  
  /*
   * Sets the stiffness of node (p_i, p_j, 0).
   * Note:
   *   Also updates the stiffnesses of the ancestor nodes.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_stiffness - stiffness
   */
  void set_stiffness(const int p_i, const int p_j, const double p_stiffness);
  
  /*
   * Applies a translation to node (p_i, p_j, 0).
   * Note:
   *   Also updates the positions of the ancestor nodes.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_translation - translation
   */
  void apply_translation(const int p_i, const int p_j, const Vector& p_translation);
  
  /*
   * Applies a rotation to node (p_i, p_j, 0).
   * Note:
   *   Also updates the positions of the ancestor nodes.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_center - center of rotation
   *   p_angle - rotation angle
   */
  void apply_rotation(
    const int p_i,
    const int p_j,
    const Vector& p_center,
    const double p_angle
  );
  
  /*
   * Applies scaling to node (p_i, p_j, 0).
   * Note:
   *   Also updates the positions of the ancestor nodes.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_center - center of scaling
   *   p_factor - scaling factor
   */
  void apply_scaling(
    const int p_i,
    const int p_j,
    const Vector& p_center,
    const double p_factor
  );
  
  /*
   * Returns the upsampled nodes.
   * Note:
   *   See the documentation for more details.
   * Returns:
   *   A vector containing shared pointers to the upsampled nodes.
   */
  std::vector<std::shared_ptr<Node>> get_upsampled();
  
  /*
   * Returns the active nodes.
   * Note:
   *   See the documentation for more details.
   * Returns:
   *   A vector containing shared pointers to the active nodes.
   */
  std::vector<std::shared_ptr<Node>> get_active();
  
  /*
   * Updates the downsampled nodes.
   */
  void downsample();
  
  /*
   * Updates the upsampled nodes.
   * Note:
   *   See the documentation for more details.
   */
  void upsample(const std::vector<std::shared_ptr<Node>>& p_nodes);
  
  /*
   * Updates the leaf nodes (i, j, 0) of the quadtree.
   */
  void sample_top();
  
  /*
   * Performs a single relaxation step of active nodes.
   * Note:
   *   See the documentation for more details.
   * Parameters:
   *   p_nodes - the active nodes to be relaxed
   *   p_step_size - relaxation step size
   *   p_eliminate_drift - whether to eliminate linear and angular drift of the nodes
   *   p_step_limit - limit for step size (see the documentation for more details)
   */
  void relax(
    const std::vector<std::shared_ptr<Node>>& p_nodes,
    const double p_step_size,
    const bool p_eliminate_drift,
    const double p_step_limit
  );
  
  // Copy constructor
  // Note: Copying is not allowed.
  Node(const Node& p_node) = delete;
  
public:
    
  /*
   * Constructor
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_k - quadtree k index
   *   p_c - node child index
   *   p_parent - parent of the node
   */
  Node(
    const int p_i,
    const int p_j,
    const int p_k,
    const int p_c,
    const std::shared_ptr<Node>& p_parent
  );
  
  /*
   * Creates a complete quadtree of nodes with a p_size-by-p_size grid of leaf nodes.
   * Note:
   *   This must be called by the root node of the quadtree only.
   *   Parameter p_size must be a power of and greater than two.
   * Parameters:
   *   p_size - side length of leaf node grid
   * Returns:
   *   A shared pointer to the root node of the quadtree
   *   A null pointer if called by a node other than the root node or if p_size is invalid
   */
  static std::shared_ptr<Node> create_nodes(const int p_size);
  
  /*
   * Returns the side length of the leaf node grid.
   * Note:
   *   This must be called by the root node only.
   * Returns:
   *   The side length of the leaf node grid
   */
  int get_size();
  
  /*
   * Returns the position of node (p_i, p_j, 0).
   * Note:
   *   This must be called by the root node only.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   * Returns:
   *   The position of node (p_i, p_j, 0)
   */
  Vector get_position(const int p_i, const int p_j);
  
  /*
   * Sets the position of node (p_i, p_j, 0).
   * Note:
   *   This must be called by the root node only.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   *   p_position - a position
   */
  void set_position(const int p_i, const int p_j, const Vector& p_position);
  
  /*
   * Sets the node (p_i, p_j, 0) fixed.
   * Note:
   *   This must be called by the root node only.
   *   Updates the active nodes.
   * Parameters:
   *   p_i - quadtree i index
   *   p_j - quadtree j index
   */
  void set_fixed(const int p_i, const int p_j);
  
  /*
   * Sets the nodes in the selection fixed.
   * Note:
   *   This must be called by the root node only.
   *   Updates the active nodes.
   * Parameters:
   *   p_nodes - a selection of nodes
   *   p_edges_only - a bool indicating whether to fix the edge nodes of the selection only
   */
  void set_fixed(const selection& p_nodes, const bool p_edges_only = false);
  
  /*
   * Sets the scale of the nodes in the selection.
   * Note:
   *   This must be called by the root node only.
   *   Updates the active nodes.
   * Parameters:
   *   p_nodes - a selection of nodes
   *   p_scale - a length scale
   */
  void set_scale(const selection& p_nodes, const double p_scale);
  
  /*
   * Sets the stiffness of the nodes in the selection.
   * Note:
   *   This must be called by the root node only.
   *   Updates the active nodes.
   * Parameters:
   *   p_nodes - a selection of nodes
   *   p_stiffness - a stiffness
   */
  void set_stiffness(const selection& p_nodes, const double p_stiffness);
  
  /*
   * Applies a translation to the nodes in the selection.
   * Note:
   *   This must be called by the root node only.
   *   Updates the active nodes.
   * Parameters:
   *   p_nodes - a selection of nodes
   *   p_translation - a translation
   *   p_edges only - a bool indicating whether to translate the edge nodes of the selection only
   */
  void apply_translation(
    const selection& p_nodes,
    const Vector& p_translation,
    const bool p_edges_only = false
  );
  
  /*
   * Applies a rotation to the nodes in the selection.
   * Note:
   *   This must be called by the root node only.
   *   Updates the active nodes.
   * Parameters:
   *   p_nodes - a selection of nodes
   *   p_center - the axis of rotation
   *   p_angle - a rotation angle
   *   p_edges only - a bool indicating whether to rotate the edge nodes of the selection only
   */
  void apply_rotation(
    const selection& p_nodes,
    const Vector& p_center,
    const double p_angle,
    const bool p_edges_only = false
  );
  
  /*
   * Applies a scaling to the nodes in the selection.
   * Note:
   *   This must be called by the root node only.
   *   Updates the active nodes.
   * Parameters:
   *   p_nodes - a selection of nodes
   *   p_center - the center of scaling
   *   p_factor - a scaling factor
   *   p_edges only - a bool indicating whether to scale the edge nodes of the selection only
   */
  void apply_scaling(
    const selection& p_nodes,
    const Vector& p_center,
    const double p_factor,
    const bool p_edges_only = false
  );
  
  /*
   * Relaxes the system of nodes.
   * Note:
   *   This must be called by the root node only.
   * Parameters:
   *   p_num_steps - the number of relaxation steps to be performed
   *   p_step_size - relaxation step size
   *   p_eliminate_drift - a bool indicating whether to eliminate linear and angular drift
   *   p_step_limit - relaxation step limit  (see the documentation for more details)
   */
  void relax(
    const int p_num_steps,
    const double p_step_size,
    const bool p_eliminate_drift = false,
    const double p_step_limit = 0.0
  );
  
  /*
   * Writes the nodes' details into a text file.
   * Note:
   *   This must be called by the root node only.
   * Parameters:
   *   p_filename - a filename
   */
  void write_nodes(const std::string& p_filename);
  
  /*
   * Prints a visualization of the quadtree with the nodes' types.
   * Note:
   *   This must be called by the root node only.
   */
  void print_types();
};

#endif /* NODE_H */

