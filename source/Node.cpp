/* 
 * File:   Node.cpp
 * Author: Petri Hirvonen
 * 8 April 2022
 */

#include <iostream>
#include <fstream>
// remember to build release version with -DNDEBUG --> disables assertions
#include <assert.h>
#include <omp.h>
#include <set>
#include <algorithm>

#include "Node.h"

namespace {
  // for looping over a pixel's neighbors
  // mth neighbor of pixel (i, j) has indices i + di[m] and j + dj[m]
  const int di[] { 1, 0, -1, 0 };
  const int dj[] { 0, 1, 0, -1 };
  // same but including diagonal neighbors also
  const int di_[] { 1, 1, 0, -1, -1, -1, 0, 1 };
  const int dj_[] { 0, 1, 1, 1, 0, -1, -1, -1 };
  // child indices (A node's children are not created in their logical order.)
  const int c[] {2, 3, 1, 0};
}

// *****

// private:

// Populates the quadtree with nodes.
void Node::populate() {
  if(m_k > 0) {
    // node (i, j, m_k > 0) has 2-by-2 children
    for(int v {0}; v < 2; ++v) {
      for(int u {0}; u < 2; ++u) {
        const std::shared_ptr<Node>& child {
          std::make_shared<Node>(
            2*m_i + u,
            2*m_j + v,
            m_k - 1,
            // child index according to logical order
            c[m_children.size()],
            shared_from_this()
          )
        };
        m_children.push_back(child);
        child->populate();
      }
    }
  }
}

// Checks if node (p_i, p_j, p_k) is a descendant of the current node.
bool Node::is_descendant(const int p_i, const int p_j, const int p_k) const {
  assert(("Error: Node::get: Negative node index k!", p_k >= 0));
  if(p_k < m_k) {
    const int dk {m_k - p_k};
    // extrema of i and j at quadtree level p_k of current node's descendants
    const int imin {m_i << dk};
    const int imax {(m_i + 1) << dk};
    const int jmin {m_j << dk};
    const int jmax {(m_j + 1) << dk};
    if(p_i >= imin && p_i < imax && p_j >= jmin && p_j < jmax) {
      return true;
    }
  }
  return false;
}

// Returns the root node of the quadtree.
std::shared_ptr<Node> Node::get_root() {
  std::shared_ptr<Node> ancestor;
  for(
    ancestor = shared_from_this();
    // parent exists?
    ancestor->m_parent.lock();
    ancestor = ancestor->m_parent.lock()
  );
  return ancestor;
}

// Returns the node (p_i, p_j, p_k).
// Note: The node can only be found if it is a descendant of the current node.
std::shared_ptr<Node> Node::get_node_recursive(
  const int p_i,
  const int p_j,
  const int p_k
) {
  if(p_i == m_i && p_j == m_j && p_k == m_k) {
    return shared_from_this();
  }
  if(is_descendant(p_i, p_j, p_k)) {
    std::shared_ptr<Node> node;
    for(const std::shared_ptr<Node>& child : m_children) {
      node = child->get_node_recursive(p_i, p_j, p_k);
      if(node) {
        return node;
      }
    }
    // this should be reached only if the quadtree is corrupted somehow
    assert(("Error: Node::get_node_recursive: Failed to get the node!", node));
  }
  return nullptr;
}

// Returns the node (p_i, p_j, p_k).
std::shared_ptr<Node> Node::get_node(const int p_i, const int p_j, const int p_k) {
  return get_root()->get_node_recursive(p_i, p_j, p_k);
}

// Sets the nodes' neighbors in a complete quadtree.
void Node::set_neighbors() {
  for(int m {0}; m < 4; ++m) {
    m_neighbors.push_back(get_node(m_i + di[m], m_j + dj[m], m_k));
  }
  for(const std::shared_ptr<Node>& child : m_children) {
    child->set_neighbors();
  }
}

// Updates the active nodes from around the current node in the quadtree.
// Note: See the documentation for more details.
void Node::set_active() {
  bool is_modified {false};
  if(m_type == Types::TOP) {
    m_type = Types::ACTIVE;
    is_modified = true;
  }
  // nearest and diagonal neighbors
  for(int n {0}; n < 8; ++n) {
    const std::shared_ptr<Node>& neighbor {get_node(m_i + di_[n], m_j + dj_[n], m_k)};
    if(neighbor) {
      if(neighbor->m_type == Types::TOP) {
        neighbor->m_type = Types::ACTIVE;
        is_modified = true;
      }
      const std::shared_ptr<Node>& parent {neighbor->m_parent.lock()};
      if(parent) {
        for(const std::shared_ptr<Node>& sibling : parent->m_children) {
          if(sibling->m_type == Types::TOP) {
            sibling->m_type = Types::ACTIVE;
            is_modified = true;
          }
        }
        if(is_modified) {
          parent->set_active();
        }
      }
    }
  }
}

// Sets the bottom nodes in the quadtree.
// Note: See the documentation for more details.
void Node::set_bottom() {
  if(m_type == Types::ACTIVE) {
    for(const std::shared_ptr<Node>& child : m_children) {
      if(child->m_type == Types::ACTIVE) {
        m_type = Types::BOTTOM;
        child->set_bottom();
      }
    }
  }
}

// Sets the sampled nodes in the quadtree.
// Note: See the documentation for more details.
void Node::set_sampled() {
  if(m_type == Types::ACTIVE) {
    for(const std::weak_ptr<Node>& w_neighbor : m_neighbors) {
      const std::shared_ptr<Node>& neighbor {w_neighbor.lock()};
      if(
        neighbor
        && (neighbor->m_type == Types::TOP || neighbor->m_type == Types::BOTTOM)
      ) {
        neighbor->m_type = Types::SAMPLED;
      }
    }
  }
  else {
    for(const std::shared_ptr<Node>& child : m_children) {
      child->set_sampled();
    }
  }
}

// Checks if the current node is at the edge of the selection.
bool Node::is_edge(const selection& p_nodes) {
  assert((
    "Error: Node::is_edge: This must be called by leaf nodes only!",
    m_k == 0
  ));
  const int size {1 << get_root()->m_k};
  assert((
    "Error: Node::is_edge: This must be called by nodes in the selection only!",
    p_nodes.count(size*m_j + m_i) == 1
  ));
  // if any of the current node's existing neighbors are not part of the selection 
  // --> the node is at the edge of the selection
  return std::any_of(
    m_neighbors.begin(), m_neighbors.end(),
    [&](const std::weak_ptr<Node>& w_neighbor) {
      const std::shared_ptr<Node>& neighbor {w_neighbor.lock()};
      return neighbor && p_nodes.count(size*neighbor->m_j + neighbor->m_i) == 0;
    }
  );
}

// Sets the equilibrium length scale of node (p_i, p_j, 0).
// Note: Also updates the equilibrium length scales of the ancestor nodes.
void Node::set_scale(const int p_i, const int p_j, const double p_scale) {
  if(m_k == 0 && m_i == p_i && m_j == p_j) {
    m_scale = p_scale;
  }
  else if(is_descendant(p_i, p_j, 0)) {
    double sum {0.0};
    for(const std::shared_ptr<Node>& child : m_children) {
      child->set_scale(p_i, p_j, p_scale);
      sum += child->m_scale;
    }
    m_scale = 0.25*sum;
  }
}

// Sets the stiffness of node (p_i, p_j, 0).
// Note: Also updates the stiffnesses of the ancestor nodes.
void Node::set_stiffness(const int p_i, const int p_j, const double p_stiffness) {
  if(m_k == 0 && m_i == p_i && m_j == p_j) {
    m_stiffness = p_stiffness;
  }
  else if(is_descendant(p_i, p_j, 0)) {
    double sum {0.0};
    for(const std::shared_ptr<Node>& child : m_children) {
      child->set_stiffness(p_i, p_j, p_stiffness);
      sum += child->m_stiffness;
    }
    m_stiffness = 0.25*sum;
  }
}

// Applies a translation to node (p_i, p_j, 0).
// Note: Also updates the positions of the ancestor nodes.
void Node::apply_translation(const int p_i, const int p_j, const Vector& p_translation) {
  if(m_k == 0 && m_i == p_i && m_j == p_j) {
    m_position += p_translation;
  }
  else if(is_descendant(p_i, p_j, 0)) {
    Vector sum;
    for(const std::shared_ptr<Node>& child : m_children) {
      child->apply_translation(p_i, p_j, p_translation);
      sum += child->m_position;
    }
    m_position = 0.25*sum;
  }
}

// Applies a rotation to node (p_i, p_j, 0).
// Note: Also updates the positions of the ancestor nodes.
void Node::apply_rotation(
  const int p_i,
  const int p_j,
  const Vector& p_center,
  const double p_angle
) {
  if(m_k == 0 && m_i == p_i && m_j == p_j) {
    m_position -= p_center;
    m_position = m_position.rotate(p_angle);
    m_position += p_center;
  }
  else if(is_descendant(p_i, p_j, 0)) {
    Vector sum;
    for(const std::shared_ptr<Node>& child : m_children) {
      child->apply_rotation(p_i, p_j, p_center, p_angle);
      sum += child->m_position;
    }
    m_position = 0.25*sum;
  }
}

// Applies scaling to node (p_i, p_j, 0).
// Note: Also updates the positions of the ancestor nodes.
void Node::apply_scaling(
  const int p_i,
  const int p_j,
  const Vector& p_center,
  const double p_factor
) {
  if(m_k == 0 && m_i == p_i && m_j == p_j) {
    m_position -= p_center;
    m_position *= p_factor;
    m_position += p_center;
  }
  else if(is_descendant(p_i, p_j, 0)) {
    Vector sum;
    for(const std::shared_ptr<Node>& child : m_children) {
      child->apply_scaling(p_i, p_j, p_center, p_factor);
      sum += child->m_position;
   }
    m_position = 0.25*sum;
  }
}

// Returns the upsampled nodes.
// Note: See the documentation for more details.
std::vector<std::shared_ptr<Node>> Node::get_upsampled() {
  std::vector<std::shared_ptr<Node>> upsampled;
  if(m_parent.lock() && m_parent.lock()->m_type == Types::ACTIVE) {
    if(m_type == Types::SAMPLED) {
      upsampled.push_back(shared_from_this());
    }
    return upsampled;
  }
  for(const std::shared_ptr<Node>& child : m_children) {
    const std::vector<std::shared_ptr<Node>> upsampled_ {child->get_upsampled()};
    upsampled.insert(upsampled.end(), upsampled_.begin(), upsampled_.end());
  }
  return upsampled;
}

// Returns the active nodes.
// Note: See the documentation for more details.
std::vector<std::shared_ptr<Node>> Node::get_active() {
  std::vector<std::shared_ptr<Node>> active;
  if(m_type == Types::ACTIVE) {
    active.push_back(shared_from_this());
    return active;
  }
  for(const std::shared_ptr<Node>& child : m_children) {
    const std::vector<std::shared_ptr<Node>> active_ {child->get_active()};
    active.insert(active.end(), active_.begin(), active_.end());
  }
  return active;
}

// Updates the downsampled nodes.
void Node::downsample() {
  if(m_type == Types::SAMPLED) {
    m_position = Vector();
    for(const std::shared_ptr<Node>& child : m_children) {
      child->downsample();
      m_position += child->m_position;
    }
    m_position *= 0.25;
  }
  else if(m_type == Types::BOTTOM) {
    for(const std::shared_ptr<Node>& child : m_children) {
      child->downsample();
    }
  }
}

// Updates the upsampled nodes.
// Note: See the documentation for more details.
void Node::upsample(const std::vector<std::shared_ptr<Node>>& p_nodes) {
  assert((
    "Error: Node::upsample: This must be called by the root node only!",
    shared_from_this() == get_root()
  ));
  #pragma omp parallel for
  for(int n = 0; n < p_nodes.size(); ++n) {
    const std::shared_ptr<Node>& node {p_nodes[n]};
    const std::shared_ptr<Node>& parent {node->m_parent.lock()};
    Vector offset;
    double sign {1.0};
    std::shared_ptr<Node> neighbor {parent->m_neighbors[node->m_c].lock()};
    // edge node without neighbor?
    if(!neighbor) {
      // use the opposite neighbor
      neighbor = parent->m_neighbors[(node->m_c + 2)%4].lock();
      sign = -1.0;
    }
    offset += sign*(neighbor->m_position - parent->m_position);
    sign = 1.0;
    neighbor = parent->m_neighbors[(node->m_c + 1)%4].lock();
    if(!neighbor) {
      neighbor = parent->m_neighbors[(node->m_c + 3)%4].lock();
      sign = -1.0;
    }
    offset += sign*(neighbor->m_position - parent->m_position);
    node->m_update = parent->m_position + 0.25*offset;
  }
  #pragma omp parallel for
  for(int n = 0; n < p_nodes.size(); ++n) {
    p_nodes[n]->m_position = p_nodes[n]->m_update;
  }
}

// Updates the leaf nodes (i, j, 0) of the quadtree.
void Node::sample_top() {
  assert((
    "Error: Node::sample_top: This must be called by the root node only!",
    shared_from_this() == get_root()
  ));
  // start almost at the bottom of the quadtree
  int size {4};
  for(int k {m_k - 2}; k >= 0; --k) {
    #pragma omp parallel for
    for(int j = 0; j < size; ++j) {
      for(int i {0}; i < size; ++i) {
        const std::shared_ptr<Node>& node {get_node(i, j, k)};
        if(node->m_type == Types::TOP) {
          const std::shared_ptr<Node>& parent {node->m_parent.lock()};
          Vector offset;
          double sign {1.0};
          std::shared_ptr<Node> neighbor {parent->m_neighbors[node->m_c].lock()};
          // edge node without neighbor?
          if(!neighbor) {
            // use the opposite neighbor
            neighbor = parent->m_neighbors[(node->m_c + 2)%4].lock();
            sign = -1.0;
          }
          offset += sign*(neighbor->m_position - parent->m_position);
          sign = 1.0;
          neighbor = parent->m_neighbors[(node->m_c + 1)%4].lock();
          if(!neighbor) {
            neighbor = parent->m_neighbors[(node->m_c + 3)%4].lock();
            sign = -1.0;
          }
          offset += sign*(neighbor->m_position - parent->m_position);
          node->m_position = parent->m_position + 0.25*offset;
        }
      }
    }
    size *= 2;
  }
}

// Performs a single relaxation step of active nodes.
// Note: See the documentation for more details.
void Node::relax(
  const std::vector<std::shared_ptr<Node>>& p_nodes,
  const double p_step_size,
  const bool p_eliminate_drift,
  const double p_step_limit
) {
  assert((
    "Error: Node::relax: This must be called by the root node only!",
    shared_from_this() == get_root()
  ));
  
  // reset nodes' steps
  #pragma omp parallel for
  for(int n = 0; n < p_nodes.size(); ++n) {
    p_nodes[n]->m_update = Vector();
  }
  
  // ***
  
  // update steps
  #pragma omp parallel for
  // loop over active nodes
  for(int n = 0; n < p_nodes.size(); ++n) {
    const std::shared_ptr<Node>& node {p_nodes[n]};
    const double sqrt2 {std::sqrt(2.0)};
    // loop over nearest neighbors
    for(int m {0}; m < 4; ++m) {
      const std::shared_ptr<Node>& neighbor_m {node->m_neighbors[m].lock()};
      // node isn't an edge node without neighbor?
      if(neighbor_m) {
        const Vector diff_m {neighbor_m->m_position - node->m_position};
        const double length_m {diff_m.magnitude()};
        if(length_m > 0.0) {
          const double lsqrt2 {length_m*sqrt2};
          const std::vector<double> other_lengths {lsqrt2, 2.0*length_m, lsqrt2};
          const Vector direction {diff_m.normalize()};
          const double stiffness {0.5*(node->m_stiffness + neighbor_m->m_stiffness)};
          // the current node and its neighbor have equilibrium length scales?
          if(node->m_scale > 0.0 && neighbor_m->m_scale > 0.0) {
            const double length0 {0.5*(node->m_scale + neighbor_m->m_scale)*(1 << node->m_k)};
            const double magnitude {stiffness*(length_m - length0)/length0};
            node->m_update += magnitude*direction;
          }
          // loop over the three other nearest neighbors
          for(int l {0}; l < 3; ++l) {
            const std::shared_ptr<Node>& neighbor_l {node->m_neighbors[(m + l + 1)%4].lock()};
            // neighbor_m isn't an edge node without neighbor?
            if(neighbor_l) {
              const Vector diff_l {neighbor_l->m_position - neighbor_m->m_position};
              const double length_l {diff_l.magnitude()};
              if(length_l > 0.0) {
                const double magnitude {stiffness*(other_lengths[l] - length_l)/length_l};
                node->m_update += magnitude*direction;
              }
            }
          }
        }
      }
    }
  }
  
  // ***
  
  // eliminate linear and angular drift
  if(p_eliminate_drift) {
    double mass {0.0};
    Vector center;
    #pragma omp parallel
    {
      double mass_local {0.0};
      Vector center_local;
      #pragma omp for nowait
      for(int n = 0; n < p_nodes.size(); ++n) {
        const std::shared_ptr<Node>& node {p_nodes[n]};
        const int size {1 << node->m_k};
        const int size2 {size*size};
        mass_local += size2;
        center_local += size2*(node->m_position + p_step_size*node->m_update);
      }
      #pragma omp critical
      {
        mass += mass_local;
        center += center_local;
      }
    }
    const double _step_size {1.0/p_step_size};
    center /= mass;
    const Vector position_root {m_position};
    const Vector offset = {center - position_root};
    double angle {0.0};
    #pragma omp parallel
    {
      double angle_local {0.0};
      #pragma omp for nowait
      for(int n = 0; n < p_nodes.size(); ++n) {
        const std::shared_ptr<Node>& node {p_nodes[n]};
        node->m_update -= offset*_step_size;
        const Vector diff {node->m_position - position_root};
        const double magnitude {diff.magnitude()};
        angle_local += Vector::cross_product(diff, node->m_update)/(magnitude*magnitude);
      }
      #pragma omp critical
      {
        angle += angle_local;
      }
    }
    angle /= p_nodes.size();
    #pragma omp parallel for
    for(int n = 0; n < p_nodes.size(); ++n) {
      p_nodes[n]->m_update -= Vector::cross_product(
        angle,
        p_nodes[n]->m_position - position_root
      );
    }
  }
  
  // ***
  
  // limit step size
  if(p_step_limit > 0.0) {
    double mean {0.0};
    #pragma omp parallel
    {
      double mean_local {0.0};
      #pragma omp for nowait
      for(int n = 0; n < p_nodes.size(); ++n) {
        const std::shared_ptr<Node>& node {p_nodes[n]};
        mean += node->m_update.magnitude();
      }
      #pragma omp critical
      {
        mean += mean_local;
      }
    }
    mean /= p_nodes.size();
    const double coef {p_step_limit*mean};
    const double _coef {1.0/coef};
    #pragma omp parallel for
    for(int n = 0; n < p_nodes.size(); ++n) {
      const std::shared_ptr<Node>& node {p_nodes[n]};
      const double magnitude {node->m_update.magnitude()};
      if(magnitude > 0.0) {
        node->m_update = coef/magnitude*std::tanh(magnitude*_coef)*node->m_update;
      }
    }
  }
  
  // ***
  
  // update positions
  #pragma omp parallel for
  for(int n = 0; n < p_nodes.size(); ++n) {
    p_nodes[n]->m_position += p_step_size*p_nodes[n]->m_update;
  }
}
  
// *****

// public:

// Constructor
Node::Node(
  const int p_i,
  const int p_j,
  const int p_k,
  const int p_c,
  const std::shared_ptr<Node>& parent
) : m_i {p_i}, m_j {p_j}, m_k {p_k}, m_c {p_c}, m_parent {parent} {
  if(
    m_k < 0
    || m_i < 0
    || m_j < 0
    || m_c < 0
    || m_c >= 4
    || m_parent.lock() ?
      !m_parent.lock()->get_root()->is_descendant(m_i, m_j, m_k)
      : false
  ) {
    throw std::invalid_argument("Error: Node::Node: Invalid node indices!");
  }
  const int offset {1 << p_k};
  m_position = Vector(offset*p_i + 0.5*(offset - 1), offset*p_j + 0.5*(offset - 1));
}

/* Creates a complete quadtree of nodes with a p_size-by-p_size grid of leaf nodes.
 * Note:
 *   This must be called by the root node of the quadtree only.
 *   Parameter p_size must be a power of and greater than two.
 */
std::shared_ptr<Node> Node::create_nodes(const int p_size) {
  if(p_size < 4 || p_size & (p_size - 1) != 0) {
    std::cout << "Warning: Node::create_nodes: The size must be a power of and greater "
      "than 2! Returned a nullptr!" << std::endl;
    return nullptr;
  }
  if(p_size > 1 << 15) {
    std::cout << "Warning: Node::create_nodes: The linear grid size must not be greater "
      "than 2^15! Returned a nullptr" << std::endl;
    return nullptr;
  }
  int k;
  for(k = 0; p_size >> k > 1; ++k);
  const std::shared_ptr<Node>& root {std::make_shared<Node>(0, 0, k, 0, nullptr)};
  root->populate();
  root->set_neighbors();
  return root;
}

// Returns the side length of the leaf node grid.
// Note: This must be called by the root node only.
int Node::get_size() {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::get_size: This must be called by the root node only! "
      "Returned -1!" << std::endl;
    return -1;
  }
  return 1 << m_k;
}

// Returns the position of node (p_i, p_j, 0).
// Note: This must be called by the root node only.
Vector Node::get_position(const int p_i, const int p_j) {
  if(shared_from_this() != get_root()) {
    std::cout << "Error: Node::get_position: This must be called by the root node only! "
      "Returned a zero vector!" << std::endl;
    return Vector();
  }
  if(!is_descendant(p_i, p_j, 0)) {
    std::cout << "Warning: Node::get_position: Invalid node indices! Returned a zero "
      "vector!" << std::endl;
  }
  return get_node(p_i, p_j, 0)->m_position;
}

// Sets the position of node (p_i, p_j, 0).
// Note: This must be called by the root node only.
void Node::set_position(const int p_i, const int p_j, const Vector& p_position) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::set_position: This must be called by the root node only! "
      "No changes were made!" << std::endl;
    return;
  }
  if(!is_descendant(p_i, p_j, 0)) {
    std::cout << "Warning: Node::set_position: Invalid node indices! No changes were made!"
      << std::endl;
    return;
  }
  const std::shared_ptr<Node>& node {get_node(p_i, p_j, 0)};
  if(node) {
    node->m_position = p_position;
    node->set_active();
  }
}

/* Sets the node (p_i, p_j, 0) fixed.
 * Note:
 *   This must be called by the root node only.
 *   Updates the active nodes.
*/
void Node::set_fixed(const int p_i, const int p_j) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::set_fixed: This must be called by the root node only! "
      "No changes were made!" << std::endl;
    return;
  }
  if(!is_descendant(p_i, p_j, 0)) {
    std::cout << "Warning: Node::set_fixed: Invalid node indices! No changes were made!"
      << std::endl;
    return;
  }
  const std::shared_ptr<Node>& node {get_node(p_i, p_j, 0)};
  if(node) {
    node->m_type = Types::FIXED;
    node->set_active();
  }
}

/*
 * Sets the nodes in the selection fixed.
 * Note:
 *   This must be called by the root node only.
 *   Updates the active nodes.
 */
void Node::set_fixed(const selection& p_nodes, const bool p_edges_only) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::set_fixed: This must be called by the root node only! "
      "No changes were made!" << std::endl;
    return;
  }
  const int size {1 << m_k};
  if(p_nodes.size() == 0 || *p_nodes.begin() < 0 || *p_nodes.rbegin() >= size*size) {
    std::cout << "Warning: Node::set_fixed: Invalid selection! No changes were made!"
      << std::endl;
    return;
  }
  for(const int n : p_nodes) {
    const int i {n%size};
    const int j {n/size};
    const std::shared_ptr<Node>& node {get_node(i, j, 0)};
    const bool is_edge {node->is_edge(p_nodes)};
    if(!p_edges_only || is_edge) {
      node->set_active();
      node->m_type = Types::FIXED;
    }
  }
}

/*
 * Sets the scale of the nodes in the selection.
 * Note:
 *   This must be called by the root node only.
 *   Updates the active nodes.
 */
void Node::set_scale(const selection& p_nodes, const double p_scale) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::set_scale: This must be called by the root node only! "
      "No changes were made!" << std::endl;
    return;
  }
  const int size {1 << m_k};
  if(p_nodes.size() == 0 || *p_nodes.begin() < 0 || *p_nodes.rbegin() >= size*size) {
    std::cout << "Warning: Node::set_scale: Invalid selection! No changes were made!"
      << std::endl;
    return;
  }
  if(p_scale < 0.0) {
    std::cout << "Warning: Node::set_scale: The length scale must not be negative! "
      "No changes were made!" << std::endl;
    return;
  }
  for(const int n : p_nodes) {
    const int i {n%size};
    const int j {n/size};
    const std::shared_ptr<Node>& node {get_node(i, j, 0)};
    set_scale(i, j, p_scale);
    if(node->is_edge(p_nodes)) {
      node->set_active();
    }
  }
}

/*
 * Sets the stiffness of the nodes in the selection.
 * Note:
 *   This must be called by the root node only.
 *   Updates the active nodes.
 */
void Node::set_stiffness(const selection& p_nodes, const double p_stiffness) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::set_stiffness: This must be called by the root node only!"
      " No changes were made!" << std::endl;
    return;
  }
  const int size {1 << m_k};
  if(p_nodes.size() == 0 || *p_nodes.begin() < 0 || *p_nodes.rbegin() >= size*size) {
    std::cout << "Warning: Node::set_stiffness: Invalid selection! No changes were made!" 
      << std::endl;
    return;
  }
  if(p_stiffness <= 0.0) {
    std::cout << "Warning: Node::set_stiffness: The stiffness must be greater than zero! "
      "No changes were made!" << std::endl;
    return;
  }
  for(const int n : p_nodes) {
    const int i {n%size};
    const int j {n/size};
    const std::shared_ptr<Node>& node {get_node(i, j, 0)};
    set_stiffness(i, j, p_stiffness);
    if(node->is_edge(p_nodes)) {
      node->set_active();
    }
  }
}

/*
 * Applies a translation to the nodes in the selection.
 * Note:
 *   This must be called by the root node only.
 *   Updates the active nodes.
 */
void Node::apply_translation(
  const selection& p_nodes,
  const Vector& p_translation,
  const bool p_edges_only
) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::translate: This must be called by the root node only! "
      "No changes were made!" << std::endl;
    return;
  }
  const int size {1 << m_k};
  if(p_nodes.size() == 0 || *p_nodes.begin() < 0 || *p_nodes.rbegin() >= size*size) {
    std::cout << "Warning: Node::translate: Invalid selection! No changes were made!"
      << std::endl;
    return;
  }
  for(const int n : p_nodes) {
    const int i {n%size};
    const int j {n/size};
    const std::shared_ptr<Node>& node {get_node(i, j, 0)};
    const bool is_edge {node->is_edge(p_nodes)};
    if(!p_edges_only || is_edge) {
      apply_translation(i, j, p_translation);
    }
    if(is_edge) {
      node->set_active();
    }
  }
}

/*
 * Applies a rotation to the nodes in the selection.
 * Note:
 *   This must be called by the root node only.
 *   Updates the active nodes.
 */
void Node::apply_rotation(
  const selection& p_nodes,
  const Vector& p_axis,
  const double p_angle,
  const bool p_edges_only
) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::rotate: This must be called by the root node only! "
      "No changes were made!" << std::endl;
    return;
  }
  const int size = {1 << m_k};
  if(p_nodes.size() == 0 || *p_nodes.begin() < 0 || *p_nodes.rbegin() >= size*size) {
    std::cout << "Warning: Node::rotate: Invalid selection! No changes were made!"
      << std::endl;
    return;
  }
  for(const int n : p_nodes) {
    const int i {n%size};
    const int j {n/size};
    const std::shared_ptr<Node>& node {get_node(i, j, 0)};
    const bool is_edge {node->is_edge(p_nodes)};
    if(!p_edges_only || is_edge) {
      apply_rotation(i, j, p_axis, p_angle);
    }
    if(is_edge) {
      node->set_active();
    }
  }
}

/*
 * Applies a scaling to the nodes in the selection.
 * Note:
 *   This must be called by the root node only.
 *   Updates the active nodes.
 */
void Node::apply_scaling(
  const selection& p_nodes,
  const Vector& p_center,
  const double p_factor,
  const bool p_edges_only
) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::scale: This must be called by the root node only! "
      "No changes were made!" << std::endl;
    return;
  }
  const int size {1 << m_k};
  if(p_nodes.size() == 0 || *p_nodes.begin() < 0 || *p_nodes.rbegin() >= size*size) {
    std::cout << "Warning: Node::scale: Invalid selection! No changes were made!"
      << std::endl;
    return;
  }
  if(p_factor <= 0.0) {
    std::cout << "Warning: Node::scale: The scaling factor must be greater than zero! "
      "No changes were made!" << std::endl;
    return;
  }
  for(const int n : p_nodes) {
    const int i {n%size};
    const int j {n/size};
    const std::shared_ptr<Node>& node {get_node(i, j, 0)};
    const bool is_edge {node->is_edge(p_nodes)};
    if(!p_edges_only || is_edge) {
      apply_scaling(i, j, p_center, p_factor);
    }
    if(is_edge) {
      node->set_active();
    }
  }
}

// Relaxes the system of nodes.
// Note: This must be called by the root node only.
void Node::relax(
  const int p_num_updates,
  const double p_step_size,
  const bool p_eliminate_drift,
  const double p_step_limit
) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::relax: This must be called by the root node only! "
      "No changes were made!" << std::endl;
    return;
  }
  if(p_num_updates <= 0) {
    std::cout << "Warning: Node::relax: The number of steps must be greater than zero! "
      "No changes were made!" << std::endl;
    return;
  }
  if(p_step_size <= 0.0) {
    std::cout << "Warning: Node::relax: The step size must be greater than zero! "
      "No changes were made!" << std::endl;
    return;
  }
  if(p_step_limit < 0.0) {
    std::cout << "Warning: Node::relax: The step limit must not be negative! "
      "No changes were made!" << std::endl;
    return;
  }
  set_bottom();
  set_sampled();
  const std::vector<std::shared_ptr<Node>> active {get_active()};
  const std::vector<std::shared_ptr<Node>> upsampled {get_upsampled()};
  for(int n {0}; n < p_num_updates; ++n) {
    downsample();
    upsample(upsampled);
    relax(active, p_step_size, p_eliminate_drift, p_step_limit);
  }
  sample_top();
}

// Writes the nodes' details into a text file.
// Note: This must be called by the root node only.
void Node::write_nodes(const std::string& p_filename) {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::write_nodes: This must be called by the root node only!"
      << std::endl;
    return;
  }
  std::ofstream file(p_filename);
  if(!file.is_open()) {
    std::cout << "Warning: Node::write_nodes: Failed to open the output file! "
      "No output was generated!" << std::endl;
    return;
  }
  int size {1 << m_k};
  file << "width height" << std::endl;
  file << size << " " << size << std::endl;
  file << std::endl;
  file << 
    "m_i m_j m_k m_c m_type m_scale m_stiffness m_position.x m_position.y "
    "m_parent.m_i m_parent.m_j neighbor0.m_i neighbor0.m_j neighbor1.m_i neighbor1.m_j "
    "neighbor2.m_i neighbor2.m_j neighbor3.m_i neighbor3.m_j child0.m_i child0.m_j "
    "child1.m_i child1.m_j child2.m_i child2.m_j child3.m_i child3.m_j"
  << std::endl;
  for(int k {m_k}; k >= 0; --k) {
    size = 1 << (m_k - k);
    for(int j {0}; j < size; ++j) {
      for(int i {0}; i < size; ++i) {
        const std::shared_ptr<Node>& node {get_node(i, j, k)};
        file
          << node->m_i << " "
          << node->m_j << " "
          << node->m_k << " "
          << node->m_c << " ";
        file << Type_labels[node->m_type] << " ";
        file << node->m_scale << " " << node->m_stiffness << " ";
        file << node->m_position.x() << " " << node->m_position.y() << " ";
        // parent
        const std::shared_ptr<Node>& parent {node->m_parent.lock()};
        if(parent) {
          file << parent->m_i << " " << parent->m_j << " ";
        }
        else {
          file << "-1 -1 ";
        }
        // neighbors
        for(const std::weak_ptr<Node>& w_neighbor : node->m_neighbors) {
          const std::shared_ptr<Node>& neighbor {w_neighbor.lock()};
          if(neighbor) {
            file << neighbor->m_i << " " << neighbor->m_j << " ";
          }
          else {
            file << "-1 -1 ";
          }
        }
        // children
        if(k == 0) {
          file << "-1 -1 -1 -1 -1 -1 -1 -1 ";
        }
        else {
          for(const std::shared_ptr<Node>& child : node->m_children) {
            file << child->m_i << " " << child->m_j << " ";
          }
        }
        file << std::endl;
      }
    }
  }
  file.close();
}

// Prints a visualization of the quadtree with the nodes' types.
// Note: This must be called by the root node only.
void Node::print_types() {
  if(shared_from_this() != get_root()) {
    std::cout << "Warning: Node::print_types: This must be called by the root node only!" 
      << std::endl;
    return;
  }
  set_bottom();
  set_sampled();
  for(int k {0}; k <= get_root()->m_k; ++k) {
    const int size {1 << get_root()->m_k - k};
    for(int j {size - 1}; j >= 0; --j) {
      for(int i {0}; i < size; ++i) {
        const std::shared_ptr<Node>& node {get_node(i, j, k)};
        if(node) {
          if(node->m_type == Types::FIXED) {
            std::cout << "F ";
          }
          else if(node->m_type == Types::ACTIVE) {
            std::cout << "A ";
          }
          else if(node->m_type == Types::SAMPLED) {
            std::cout << "S ";
          }
          else if(node->m_type == Types::TOP) {
            std::cout << "^ ";
          }
          else if(node->m_type == Types::BOTTOM) {
            std::cout << ". ";
          }
        }
        else {
          std::cout << "@ ";
        }
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}