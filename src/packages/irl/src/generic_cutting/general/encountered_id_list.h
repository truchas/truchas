// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_GENERAL_ENCOUNTERED_ID_LIST_H_
#define SRC_GENERIC_CUTTING_GENERAL_ENCOUNTERED_ID_LIST_H_

#include <functional>
#include <unordered_map>
#include <utility>

#include "src/data_structures/short_alloc.h"
#include "src/parameters/defined_types.h"

namespace IRL {

class EncounteredIdList {
  static constexpr std::size_t initial_capacity = 50;
  using StackAllocator = short_alloc::short_alloc<
      std::pair<const UnsignedIndex_t, UnsignedIndex_t>,
      static_cast<std::size_t>(initial_capacity) *
          sizeof(std::pair<const UnsignedIndex_t, UnsignedIndex_t>)>;

 public:

  EncounteredIdList(void)
      : valid_size_m(0), memory_pool_m(), encountered_list_m(initial_capacity, std::hash<UnsignedIndex_t>(), std::equal_to<UnsignedIndex_t>(), memory_pool_m) {
    encountered_list_m.reserve(initial_capacity);
  }

  void addId(const UnsignedIndex_t a_id) {
    encountered_list_m[a_id] = this->size();
    ++valid_size_m;
  }

  UnsignedIndex_t size(void) const { return valid_size_m; }

  void resize(const UnsignedIndex_t a_size) { valid_size_m = a_size; }

  bool isIdPresent(const UnsignedIndex_t a_id) {
    return encountered_list_m.find(a_id) != encountered_list_m.end()
               ? encountered_list_m.at(a_id) < valid_size_m
               : false;
  }

  ~EncounteredIdList(void) = default;

 private:
  UnsignedIndex_t valid_size_m;
  typename StackAllocator::arena_type memory_pool_m;
  std::unordered_map<UnsignedIndex_t, UnsignedIndex_t,
                     std::hash<UnsignedIndex_t>, std::equal_to<UnsignedIndex_t>,
                     StackAllocator>
      encountered_list_m;
};

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_GENERAL_ENCOUNTERED_ID_LIST_H_
