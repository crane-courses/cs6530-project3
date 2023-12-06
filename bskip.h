/*
 * ============================================================================
 *
 *       Filename:  btree.h
 *
 *         Author:  Helen Xu, hjxu@lbl.gov
 *   Organization:  Lawrence Berkeley Laboratory
 *
 * ============================================================================
 */

#ifndef _BSKIP_H_
#define _BSKIP_H_

#include <string.h>
#include <random>
#include <iostream>
#include <stack>
#include <utility>
#include <vector>
#include <limits>
#include <algorithm> 
#include <cmath>

#include "tbassert.h"
// #include "array.h"

#define WEIGHTED 0
#define MAX_HEIGHT 8
#define BINARY_SEARCH 1

// leaf node is the default
#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS> class BSkipNode {
#else
template <typename K, size_t MAX_KEYS> class BSkipNode {
#endif
public:
  uint32_t num_elts;

#if WEIGHTED
  BSkipNode<K, V, MAX_KEYS>* next;
#else
  BSkipNode<K, MAX_KEYS>* next;
#endif

	K keys[MAX_KEYS];

#if WEIGHTED
  V vals[MAX_KEYS];
#endif

#if STATS
	uint64_t num_comparisons;
#endif

  // initialize the skiplist node with default container size
  BSkipNode() {
    num_elts = 0;
    next = NULL;

#if STATS
		num_comparisons = 0;
#endif
  }

	// TODO: this could have an optimized version depending on the container type
  inline K get_header() {
    return keys[0];
  }

	uint64_t sum_keys() {
		uint64_t result = 0;
		for(uint32_t i = 0; i < num_elts; i++) {
			result += keys[i];
		}
		return result;
	}

	void print_keys() {
		for(int i = 0; i < num_elts; i++) {
			printf("\t\tkey[%d] = %lu\n", i, keys[i]);
		}
		printf("\n");
	}

	// return rank of key largest key at most k
	uint32_t find_key(K k) {
#if BINARY_SEARCH
		return find_index_binary(k);
#else
#if DEBUG_PRINT
		print_keys();
#endif
		return find_index_linear(k);
#endif
	}

	// add key elt at rank, shifting everything down if necessary
	void insert_key_at_rank(uint32_t rank, K elt) {
		// keep track of num_elts here
#if DEBUG
		assert(num_elts + 1 <= MAX_KEYS);
#endif

		// shift everything over by 1
		memmove(keys + rank + 1, keys + rank, (num_elts - rank) * sizeof(K));
		// set it
		keys[rank] = elt;		
	}	

	inline K get_key_at_rank(uint32_t rank) {
		return keys[rank];
	}

	inline void set_key_at_rank(uint32_t rank, K k) {
#if DEBUG
		assert(rank < num_elts);
#endif
		keys[rank] = k;
	}

#if WEIGHTED
	// add key elt at rank, shifting everything down if necessary
	void insert_val_at_rank(uint32_t rank, V elt) {
		// keep track of num_elts here
		assert(num_elts + 1 < MAX_KEYS);

		// shift everything over by 1
		memmove(vals + rank + 1, vals + rank, (num_elts - rank) * sizeof(V));
		// set it
		vals[rank] = elt;		
	}
#endif

// move keys starting from starting_rank into dest node starting from dest_rank
#if WEIGHTED
	int split_keys(BSkipNode<K, V, MAX_KEYS>* dest, uint32_t starting_rank, uint32_t dest_rank = 1) {
#else
	int split_keys(BSkipNode<K, MAX_KEYS>* dest, uint32_t starting_rank, uint32_t dest_rank = 1) {
#endif
		uint32_t num_elts_to_move = num_elts - starting_rank;
#if DEBUG
		assert(starting_rank <= num_elts);
		assert(num_elts_to_move <= num_elts);
		assert(num_elts_to_move < MAX_KEYS);
#endif

		memmove(dest->keys + dest_rank, keys + starting_rank, num_elts_to_move * sizeof(K));

		num_elts = starting_rank;
		dest->num_elts += num_elts_to_move;

		return num_elts_to_move;
	}

private:
  uint32_t find_index_linear(K k) {
		uint32_t i;
#if DEBUG
		for(i = 1; i < num_elts; i++) {
			assert(keys[i] > keys[i-1]); 
		}
		assert(k >= keys[0]);
#endif
		for (i = 0; i < num_elts; i++) {
#if STATS
			num_comparisons++;
#endif
			if (keys[i] < k)
				continue;
			else if (k == keys[i]) {
				return i;
			} else
				break;
		}
		return i - 1;
	}

	// TODO: add binary search
  uint32_t find_index_binary(K k) {
		uint32_t left = 0;
		uint32_t right = num_elts - 1;
		while (left <= right) {
			int mid = left + (right - left) / 2;
			if (keys[mid] == k) {
				left = mid;
				break;
			} else if (keys[mid] < k) {
				left = mid + 1;
			} else {
				if (mid == 0) {
					break;
				}
				right = mid - 1;
			}
		}
		if (left == num_elts || keys[left] > k) {
#if DEBUG
			assert(left > 0);
#endif
			left--;
		}

#if DEBUG
	 	tbassert(left < num_elts, "left = %u, num elts %u\n", left, num_elts);
		tbassert(left == find_index_linear(k), "k %lu, binary = %u, linear = %u\n", k, left, find_index_linear(k));
#endif
		return left;
	}
};

// extra case for internal node
#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS> class BSkipNodeInternal: public BSkipNode<K, V, MAX_KEYS> {
#else
template <typename K, size_t MAX_KEYS> class BSkipNodeInternal : public BSkipNode<K, MAX_KEYS> {
#endif
public:
#if WEIGHTED
  BSkipNode<K, V, MAX_KEYS>* children[MAX_KEYS];
#else
  BSkipNode<K, MAX_KEYS>* children[MAX_KEYS];
#endif

	BSkipNodeInternal() { }

	// add key elt at rank, shifting everything down if necessary
#if WEIGHTED
	void insert_child_at_rank(uint32_t rank, BSkipNode<K, V, MAX_KEYS>* elt) {
#else
	void insert_child_at_rank(uint32_t rank, BSkipNode<K, MAX_KEYS>* elt) {
#endif
		// shift everything over by 1
		memmove(children + rank + 1, children + rank, (this->num_elts - rank) * sizeof(K));
		// set it
		children[rank] = elt;		
	}

#if WEIGHTED
	BSkipNode<K, V, MAX_KEYS>* get_child_at_rank(uint32_t rank) {
#else
	BSkipNode<K, MAX_KEYS>* get_child_at_rank(uint32_t rank) {
#endif
#if DEBUG
		tbassert(this->num_elts >= rank, "num elts %d, asked for rank %d\n", this->num_elts, rank);
#endif
		return children[rank];
	}

#if WEIGHTED
	void set_child_at_rank(uint32_t rank, BSkipNode<K, V, MAX_KEYS>* elt) {
#else
	void set_child_at_rank(uint32_t rank, BSkipNode<K, MAX_KEYS>* elt) {
#endif

#if DEBUG
		tbassert(this->num_elts >= rank, "num elts %d, asked for rank %d\n", this->num_elts, rank);
#endif
		children[rank] = elt;
	}

#if WEIGHTED
	void move_children(BSkipNodeInternal<K, V, MAX_KEYS>* dest, uint32_t starting_rank, uint32_t num_elts_to_move, uint32_t dest_rank = 1) {
#else
	void move_children(BSkipNodeInternal<K, MAX_KEYS>* dest, uint32_t starting_rank, uint32_t num_elts_to_move, uint32_t dest_rank = 1) {
#endif

#if DEBUG
		assert(num_elts_to_move < MAX_KEYS);
#endif 
#if WEIGHTED
		memmove(dest->children + dest_rank, children + starting_rank, num_elts_to_move * sizeof(BSkipNode<K, V, MAX_KEYS>*));
#else
		memmove(dest->children + dest_rank, children + starting_rank, num_elts_to_move * sizeof(BSkipNode<K, MAX_KEYS>*));
#endif
	}
};

#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS> class BSkip {
#else
template <typename K, size_t MAX_KEYS> class BSkip {
#endif

public:
  int curr_max_height; // max height of current sl
  double promotion_probability; // promotion probability
  int boost;
  bool has_zero = false; // using 0 as first sentinel, so just keep it separate

	uint64_t elts_in_sl;

  std::random_device random_device;
  std::mt19937 random_engine{random_device()};
  // headers for each node
#if WEIGHTED
  BSkipNode<K, V, MAX_KEYS>* headers[MAX_HEIGHT];
#else
  BSkipNode<K, MAX_KEYS>* headers[MAX_HEIGHT];
#endif

#if WEIGHTED
  bool insert(K k, V v);
  BSkipNode<K, V, MAX_KEYS> *find(K k);
#else
  bool insert(K k, uint32_t idx);
  BSkipNode<K, MAX_KEYS> *find(K k);
#endif
	uint32_t last_promoted = 0;

	uint64_t sum();
	uint64_t psum();
	void get_size_stats();
	void get_avg_comparisons();
  int flip_coins();
  int node_size;

	// testing scaling down variance
	K s;

  BSkip(double p, uint32_t _boost = 1) {
	  /*
		elts_in_sl = 0;

		// flip a coin every s elements added
		s = p / 2;

	  node_size = p;
    // https://cpppatterns.com/patterns/flip-a-biased-coin.html
    promotion_probability = (double)s/(double)p;
		printf("s = %lu, promotion prob %f\n", s, promotion_probability);
    */
	promotion_probability = (double)(1.0)/(double)p;
	boost = _boost;
    curr_max_height = 0;
	elts_in_sl = 0;
    
    // initialize bskip with array of headers
    for(int i = 0; i < MAX_HEIGHT; i++) {
			if (i > 0) {
#if WEIGHTED
				headers[i] = new BSkipNodeInternal<K, V, MAX_KEYS>();
				auto end_sentinel = new BSkipNodeInternal<K, V, MAX_KEYS>();
#else
				headers[i] = new BSkipNodeInternal<K, MAX_KEYS>();
				auto end_sentinel = new BSkipNodeInternal<K, MAX_KEYS>();
#endif

				end_sentinel->num_elts++;
				end_sentinel->set_key_at_rank(0, std::numeric_limits<K>::max());
				headers[i]->next = end_sentinel;

				// TODO: do the end sentinels need a down pointer?
			} else {
#if WEIGHTED
				headers[i] = new BSkipNode<K, V, MAX_KEYS>();
				auto end_sentinel = new BSkipNode<K, V, MAX_KEYS>();
#else
				headers[i] = new BSkipNode<K, MAX_KEYS>();
				auto end_sentinel = new BSkipNode<K, MAX_KEYS>();
#endif
				end_sentinel->num_elts++;
				end_sentinel->set_key_at_rank(0, std::numeric_limits<K>::max());
				headers[i]->next = end_sentinel;
			}
		}

    // set min elt at start of header
    for(int i = MAX_HEIGHT - 1; i >= 0; i--) {
			headers[i]->num_elts++;
			headers[i]->set_key_at_rank(0, std::numeric_limits<K>::min());
			if (i > 0) {
#if WEIGHTED
				((BSkipNodeInternal<K, V, MAX_KEYS>*)headers[i])->children[0] = headers[i-1];
#else
				((BSkipNodeInternal<K, MAX_KEYS>*)headers[i])->set_child_at_rank(0, headers[i-1]);
#endif
			}
#if WEIGHTED
			headers[i]->vals[0] = std::numeric_limits<V>::min();
#endif
    }

#if DEBUG
		for(int i = 0; i < MAX_HEIGHT; i++) {
			assert(headers[i]->num_elts == 1);
			assert(headers[i]->next->num_elts == 1);
		}
#endif
  }

	private:
		void sum_helper(std::vector<uint64_t>& sums, BSkipNode<K, MAX_KEYS>* node, int level, K max);
};


#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS>
#else
template <typename K, size_t MAX_KEYS>
#endif
#if WEIGHTED
int BSkip<K, V, MAX_KEYS>::flip_coins() {
#else
int BSkip<K, MAX_KEYS>::flip_coins() {
#endif
	int result = 0;
	if (elts_in_sl % s == 0) {

#if DEBUG_PRINT
		printf("flip coin at elt %lu\n", elts_in_sl);
#endif
		// std::uniform_real_distribution<double> distribution(0.0, 1.0);
		std::negative_binomial_distribution<int> distribution(boost, 1 - promotion_probability);
		result = distribution(random_engine);
		int remainder = result % boost;
		result = result / boost;
		// std::uniform_real_distribution<double> distr(0.0, 1.0);
		// if (distr(random_engine) <= (double)remainder / (double)(boost)) {
		// 	result += 1;
		// }
		// while(distribution(random_engine) <= promotion_probability && result < MAX_HEIGHT - 1) { 
			// result++; 
		// }	
	}
	elts_in_sl++;
	return std::min(result, MAX_HEIGHT - 1);
	/*
  int result = 0;
  std::bernoulli_distribution coin_distribution{promotion_probability};
  while(coin_distribution(random_engine) && result < MAX_HEIGHT) { result++; }
	if (result > MAX_HEIGHT - 1) { result = MAX_HEIGHT - 1; }
  return result;
	*/
}

#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS>
#else
template <typename K, size_t MAX_KEYS>
#endif
#if WEIGHTED
bool BSkip<K, V, MAX_KEYS>::insert(K k, V v) {
#else
bool BSkip<K, MAX_KEYS>::insert(K k, uint32_t idx) {
#endif
	// special case for 0 since we use it as the front sentinel
	if (k == 0ULL) { 
		has_zero = true; 
		return true; 
	}
	
  // flip coins to determine your promotion level
  int level_to_promote = flip_coins();

#if DEBUG_PRINT
	if (level_to_promote > 0) {
		if (last_promoted > 0) {
			printf("dist to last promoted = %u\n", idx - last_promoted);
		}
		last_promoted = idx;
	}
#endif

#if DEBUG
  assert(level_to_promote < MAX_HEIGHT); 
#endif
#if DEBUG_PRINT
  printf("\n\ninserting key %lu, promoted to level %d\n", k, level_to_promote); 
#endif
  // update max level if it has been raised
  if (level_to_promote > curr_max_height) {
    curr_max_height = level_to_promote;
  }
  // start inserting from the top level that exists
  auto curr_node = headers[curr_max_height];

  int rank_above; // save rank of key in parent to fixup pointers if splitting down

  // save node from above that we will insert parent into after split
#if WEIGHTED
  BSkipNodeInternal<K, V, MAX_KEYS>* node_above = NULL;
#else
  BSkipNodeInternal<K, MAX_KEYS>* node_above = NULL;
#endif
  for(int level = curr_max_height; level >= 0; level--) {
#if DEBUG_PRINT
	  printf("\nloop start: curr node head %lu at level %d\n", curr_node->get_header(), level);
#endif
#if WEIGHTED
    BSkipNode<K, V, MAX_KEYS>* next_node = curr_node->next;
#else
    BSkipNode<K, MAX_KEYS>* next_node = curr_node->next;
#endif

#if DEBUG_PRINT
		for (int i = 0; i < curr_node->num_elts; i++) {
			if (level > 0) {
				if (((BSkipNodeInternal<K, MAX_KEYS>*)curr_node)->get_child_at_rank(i) == NULL) {
					curr_node->print_keys();
					printf("missing child at rank %u\n", i);
					assert(false);
				}
			}
		}
#endif

    // find the node to insert the key in in this level
    while(next_node->get_header() < k) {
#if DEBUG_PRINT
			printf("\tcurrent header %lu, next header %lu\n", curr_node->get_header(), next_node->get_header());
#endif
#if DEBUG
			assert(curr_node->get_header() < next_node->get_header());
#endif
      curr_node = next_node;
      next_node = curr_node->next;
#if DEBUG_PRINT
			assert(curr_node->get_header() != 0);
			for (int i = 0; i < curr_node->num_elts; i++) {
				if (level > 0) {
					if (((BSkipNodeInternal<K, MAX_KEYS>*)curr_node)->get_child_at_rank(i) == NULL) {
						curr_node->print_keys();
						printf("missing child at rank %u\n", i);
						assert(false);
					}
				}
			}
#endif
		}

		// now we are at the correct node - look for the key
		uint32_t rank = curr_node->find_key(k);

#if DEBUG_PRINT
		printf("\trank found %u\n", rank);
#endif
		// if the key was found
		if(curr_node->get_key_at_rank(rank) == k) {
			printf("found elt %lu at rank %u on level %d\n", k, rank, level);
#if WEIGHTED
			// if the key is there, update values all the way down (in the weighted case)
			curr_node->set_val_at_rank(rank, v);
			level--;
			if (level > 0) {
				curr_node = curr_node->children->get_elt_at_rank(rank);
				level--;
				while(level >= 0) {
					curr_node->set_val_at_rank(0, v);
					if (level > 0) {
					  curr_node = curr_node->children->get_elt_at_rank(0);
					}
					level--;
				}
			} 
# endif
			return true;
		} else { // otherwise, this key was not found at this level
#if DEBUG_PRINT
		  printf("\tnot found at level %u\n", level);
#endif
			if (level_to_promote < level) {
				// case 1: do not promote to this level.
#if DEBUG_PRINT
				printf("\t\tcurr node num elts %u\n", curr_node->num_elts);
				printf("\t\tcurr node elt %lu, level %d\n", curr_node->get_key_at_rank(rank), level);
#endif

				// drop down a level (if you are here, you are at an internal node)
#if WEIGHTED
				curr_node = ((BSkipNodeInternal<K, V, MAX_KEYS>*)curr_node)->get_child_at_rank(rank);
#else
				curr_node = ((BSkipNodeInternal<K, MAX_KEYS>*)curr_node)->get_child_at_rank(rank);
#endif

#if DEBUG
				assert(curr_node != NULL);
#endif
				continue;
			} else if (level_to_promote == level) {
				// Case 2: insert but not split due to promotion
				// split if overfull
				if (curr_node->num_elts + 1 > MAX_KEYS) {

#if DEBUG_PRINT
					printf("splitting because overfull\n");
#endif
#if WEIGHTED
					BSkipNode<K, V, MAX_KEYS>* new_node;
#else
					BSkipNode<K, MAX_KEYS>* new_node;
#endif
					if (level > 0) {
#if WEIGHTED
						new_node = new BSkipNodeInternal<K, V, MAX_KEYS>();
#else
						new_node = new BSkipNodeInternal<K, MAX_KEYS>();
#endif
					} else {
#if WEIGHTED
						new_node = new BSkipNode<K, V, MAX_KEYS>();
#else
						new_node = new BSkipNode<K, MAX_KEYS>();
#endif
					}
					// fixup next pointers
					new_node->next = curr_node->next;
					curr_node->next = new_node;

					// do the split
					int half_keys = curr_node->num_elts / 2;

					// move second half of keys into new node
					// returns the number of elements that were moved
					// updates the number of elts in each node
					uint32_t elts_moved = curr_node->split_keys(new_node, half_keys, 0);

					// move children if necessary
					if (level > 0) {
#if DEBUG_PRINT
						printf("moving children\n");
#endif
#if WEIGHTED
						((BSkipNodeInternal<K, V, MAX_KEYS>*)curr_node)->move_children(((BSkipNodeInternal<K, V, C>*)new_node), half_keys, elts_moved, 0);
#else
						((BSkipNodeInternal<K, MAX_KEYS>*)curr_node)->move_children(((BSkipNodeInternal<K, MAX_KEYS>*)new_node), half_keys, elts_moved, 0);
#endif
					}

#if DEBUG_PRINT
				printf("split moved %u nodes\n", elts_moved);
				printf("curr node keys:\n");
				curr_node->print_keys();
				printf("new node keys\n");
				new_node->print_keys();
#endif
					// new elt goes into first node
					if (rank + 1 <= curr_node->num_elts) {
#if DEBUG
						assert(k < new_node->get_header());
						assert(k > curr_node->get_header());
#endif
						curr_node->insert_key_at_rank(rank + 1, k);
						curr_node->num_elts++;

#if DEBUG_PRINT
						printf("\tinserting key in first node at rank %u\n", rank + 1);
						curr_node->print_keys();
#endif
						// if you are an internal node, make space for the new element's child
						if (level > 0) {
	#if WEIGHTED
							((BSkipNodeInternal<K, V, MAX_KEYS>*)curr_node)->insert_child_at_rank(rank + 1, NULL);
	#else
							((BSkipNodeInternal<K, MAX_KEYS>*)curr_node)->insert_child_at_rank(rank + 1, NULL);
	#endif
	#if WEIGHTED
							BSkipNodeInternal<K, V, MAX_KEYS>* curr_node_cast = (BSkipNodeInternal<K, V, MAX_KEYS>*)curr_node;
	#else
							BSkipNodeInternal<K, MAX_KEYS>* curr_node_cast = (BSkipNodeInternal<K, MAX_KEYS>*)curr_node;
	#endif
							node_above = curr_node_cast;
							rank_above = rank + 1;
							curr_node = curr_node_cast->get_child_at_rank(rank);
						}
					} else { // insert it into the new node
#if DEBUG
						assert(k > new_node->get_header());
#endif
						int rank_to_insert = rank - curr_node->num_elts;

#if DEBUG_PRINT
						printf("rank %u, prev node num elts %u, rank to insert %u\n", rank, curr_node->num_elts, rank_to_insert);
#endif
						new_node->insert_key_at_rank(rank_to_insert + 1, k);
						new_node->num_elts++;

#if DEBUG_PRINT
						printf("insert key into new node at rank %d\n", rank_to_insert + 1);
						new_node->print_keys();
#endif

						if (level > 0) {
						// make space for the child pointer
#if WEIGHTED
							((BSkipNodeInternal<K, V, MAX_KEYS>*)new_node)->insert_child_at_rank(rank_to_insert + 1, NULL);
#else
							((BSkipNodeInternal<K, MAX_KEYS>*)new_node)->insert_child_at_rank(rank_to_insert + 1, NULL);
#endif

#if WEIGHTED
							BSkipNodeInternal<K, V, MAX_KEYS>* new_node_cast = (BSkipNodeInternal<K, V, MAX_KEYS>*)new_node;
#else
							BSkipNodeInternal<K, MAX_KEYS>* new_node_cast = (BSkipNodeInternal<K, MAX_KEYS>*)new_node;
#endif
							node_above = new_node_cast;
							rank_above = rank_to_insert + 1;
							curr_node = new_node_cast->get_child_at_rank(rank_to_insert);

#if DEBUG
							assert(curr_node);
							assert(curr_node->get_header() == new_node_cast->get_child_at_rank(rank_to_insert)->get_header());
#endif
						}
					}
				} else {
					// add key (and if needed, value) to this node 
					curr_node->insert_key_at_rank(rank + 1, k);

					// update num elts
					curr_node->num_elts++;
	#if DEBUG_PRINT
					printf("\n\tinserted %lu into level %u at rank %u\n", k, level, rank + 1);
					curr_node->print_keys();
	#endif

	#if WEIGHTED
					curr_node->vals->insert_elt_at_rank(rank + 1, v);
	#endif
					// make space for the child pointer
					if (level > 0) {
	#if WEIGHTED
						BSkipNodeInternal<K, V, MAX_KEYS>* curr_node_cast = (BSkipNodeInternal<K, V, MAX_KEYS>*)curr_node;
	#else
						BSkipNodeInternal<K, MAX_KEYS>* curr_node_cast = (BSkipNodeInternal<K, MAX_KEYS>*)curr_node;
	#endif
	#if WEIGHTED
						curr_node_cast->insert_child_at_rank(rank+1, NULL);
	#else
						curr_node_cast->insert_child_at_rank(rank+1, NULL);
	#endif

						rank_above = rank + 1; // save current rank in case we need it for a lower split
						node_above = curr_node_cast; // save current node for pointers in lower split
	#if WEIGHTED
						curr_node = curr_node_cast->get_child_at_rank(rank);
	#else
						curr_node = curr_node_cast->get_child_at_rank(rank);
	#endif
						node_above = curr_node_cast;

	#if DEBUG_PRINT
						printf("\tmake space at rank %u on level %d, node head %lu\n", rank + 1, level, curr_node->get_header());
	#endif
	#if DEBUG
						assert(curr_node != NULL);
	#endif
					}
				}
			} else { // case 3: do a split
#if DEBUG_PRINT
			  printf("split at level %d with elt %lu\n", level, k);
#endif
#if WEIGHTED
				BSkipNode<K, V, MAX_KEYS>* new_node;
#else
				BSkipNode<K, MAX_KEYS>* new_node;
#endif
				if (level > 0) {
#if WEIGHTED
					new_node = new BSkipNodeInternal<K, V, MAX_KEYS>();
#else
					new_node = new BSkipNodeInternal<K, MAX_KEYS>();
#endif
				} else {
#if WEIGHTED
			  	new_node = new BSkipNode<K, V, MAX_KEYS>();
#else
					new_node = new BSkipNode<K, MAX_KEYS>();
#endif
				}
#if DEBUG
				assert(new_node);
#endif
				new_node->insert_key_at_rank(0, k);
				new_node->num_elts++;

				// move all keys in the curr node starting from rank + 1 into the new node
				// returns the number of elements that were moved
				// updates the number of elts in each node
				uint32_t elts_moved = curr_node->split_keys(new_node, rank + 1);

#if DEBUG_PRINT
				printf("split moved %u nodes\n", elts_moved);
				printf("curr node keys:\n");
				curr_node->print_keys();
				printf("new node keys\n");
				new_node->print_keys();
#endif

# if WEIGHTED
				new_node->insert_val_at_rank(0, v);

				// TODO: add split method for vals (this is currently nonexistent)
				// move all vals in the curr node starting from rank + 1 into new node
			  curr_node->vals->split_into(new_node->vals, rank + 1);
#endif

				// if you are an internal node, make space for the child of this head
				if (level > 0) {
#if WEIGHTED
					((BSkipNodeInternal<K, V, MAX_KEYS>*)new_node)->insert_child_at_rank(0, NULL);
#else
					((BSkipNodeInternal<K, MAX_KEYS>*)new_node)->insert_child_at_rank(0, NULL);
#endif
				
#if WEIGHTED
					((BSkipNodeInternal<K, V, MAX_KEYS>*)curr_node)->move_children(((BSkipNodeInternal<K, V, C>*)new_node), rank + 1, elts_moved);
#else
					((BSkipNodeInternal<K, MAX_KEYS>*)curr_node)->move_children(((BSkipNodeInternal<K, MAX_KEYS>*)new_node), rank + 1, elts_moved);
#endif
				}

				// fixup next pointers
#if DEBUG
				assert(curr_node->get_header() < curr_node->next->get_header());
				assert(new_node->get_header() < curr_node->next->get_header());
				assert(curr_node->get_header() < new_node->get_header());
#endif
#if DEBUG_PRINT
				printf("curr header %lu, next header %lu, new header %lu\n", curr_node->get_header(), curr_node->next->get_header(), new_node->get_header());
#endif
				new_node->next = curr_node->next;
				curr_node->next = new_node;

				// fixup child pointer above
#if DEBUG_PRINT
				printf("setting above pointer at rank %d with header %lu\n", rank_above, node_above->get_header());
#endif
#if DEBUG
				assert(node_above);
				assert(node_above->get_key_at_rank(rank_above) == k);
#endif

				node_above->set_child_at_rank(rank_above, new_node);

				// if we are at an internal level, move on and save info for lower splits
				if (level > 0) {
#if WEIGHTED
					curr_node = ((BSkipNodeInternal<K, V, MAX_KEYS>*)curr_node)->get_child_at_rank(rank);
#else
					curr_node = ((BSkipNodeInternal<K, MAX_KEYS>*)curr_node)->get_child_at_rank(rank);
#endif

#if DEBUG_PRINT
					printf("level above 0, rank %u\n", rank);
#endif
#if DEBUG
					assert(curr_node);
#endif
					rank_above = 0;
#if WEIGHTED
					node_above = (BSkipNodeInternal<K, V, MAX_KEYS>*)new_node;
#else
					node_above = (BSkipNodeInternal<K, MAX_KEYS>*)new_node;
#endif
					assert(node_above);
			  }
			}
		}
  }
	return true;
}


#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS>
#else
template <typename K, size_t MAX_KEYS>
#endif
#if WEIGHTED
  BSkipNode<K, V, MAX_KEYS>* BSkip<K, V, MAX_KEYS>::find(K k) {
#else
  BSkipNode<K, MAX_KEYS>* BSkip<K, MAX_KEYS>::find(K k) {
#endif
	// TODO: deal with if if you look for 0. should this return true for the unweighted case and the val in the weighted case? why is the node the thing getting returned?
#if DEBUG_PRINT
	printf("searching for %lu\n", k);
#endif
	// start search from the top node
	auto curr_node = headers[curr_max_height];
	
	for(int level = curr_max_height; level >= 0; level--) {
		auto next_node = curr_node->next;
#if DEBUG_PRINT
		curr_node->print_keys();
#endif

		// move forward until we find the right node that contains the key range
		while(next_node->get_header() <= k) {
#if DEBUG_PRINT
			curr_node->print_keys();
#endif
			curr_node = next_node;
			next_node = curr_node->next;
		}

		// look for the largest element that is at most the search key
		uint32_t rank = curr_node->find_key(k);

		// if it is found, return the node (returns the topmost node the key is found in
		if(curr_node->get_key_at_rank(rank) == k) {
#if DEBUG_PRINT
			printf("found key at rank %u\n", rank);
			curr_node->print_keys();
#endif
			return curr_node;
		}

		// if not found, drop down a level
		if (level > 0) {
#if WEIGHTED
			curr_node = ((BSkipNodeInternal<K, V, MAX_KEYS>*)curr_node)->get_child_at_rank(rank);
#else
			curr_node = ((BSkipNodeInternal<K, MAX_KEYS>*)curr_node)->get_child_at_rank(rank);
#endif
		}
	}
	return NULL;
}

// sums all keys
#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS>
#else
template <typename K, size_t MAX_KEYS>
#endif
#if WEIGHTED
  uint64_t BSkip<K, V, MAX_KEYS>::sum() {
#else
  uint64_t BSkip<K, MAX_KEYS>::sum() {
#endif
  auto curr_node = headers[0];
	uint64_t result = 0;
	while(curr_node->get_header() < std::numeric_limits<K>::max()) {
#if DEBUG_PRINT
		printf("node header %lu, add sum %lu\n", curr_node->get_header(), curr_node->sum_keys());
#endif
		result += curr_node->sum_keys();
		curr_node = curr_node->next;
	}
	return result;
}


// sums all keys
#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS>
#else
template <typename K, size_t MAX_KEYS>
#endif
#if WEIGHTED
  uint64_t BSkip<K, V, MAX_KEYS>::sum_helper() {
#else
  void BSkip<K, MAX_KEYS>::sum_helper(std::vector<uint64_t>& sums, BSkipNode<K, MAX_KEYS>* node, int level, K local_max) {
#endif
	assert(node);
	assert(node->next);
	if(level == 0) {
		if (node->next->get_header() < local_max) {
			// cilk_spawn sum_helper(sums, node->next, level, local_max); 
			sum_helper(sums, node->next, level, local_max); 
		}
#if DEBUG_PRINT
		printf("\tnode header %lu, add sum %lu\n", node->get_header(), node->sum_keys());
#endif
		// sums[ParallelTools::getWorkerNum() * 8] += node->sum_keys();
		sums[0] += node->sum_keys();
	} else {
		BSkipNodeInternal<K, MAX_KEYS>* curr_node_cast = (BSkipNodeInternal<K, MAX_KEYS>*)node;
		if (curr_node_cast->next->get_header() < local_max) {
			sum_helper(sums, node->next, level, local_max);
		}
#if DEBUG_PRINT
		printf("start key %lu, end key %lu\n", node->get_header(), node->get_key_at_rank(node->num_elts - 1));
#endif
		for(int i = 0; i < node->num_elts; i++) {
			if (i < node->num_elts - 1) {
				sum_helper(sums, curr_node_cast->get_child_at_rank(i), level - 1, node->get_key_at_rank(i+1));
			} else {
				sum_helper(sums, curr_node_cast->get_child_at_rank(i), level - 1, node->next->get_header());
			}
		}
	}
}


// sums all keys
#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS>
#else
template <typename K, size_t MAX_KEYS>
#endif
#if WEIGHTED
  uint64_t BSkip<K, V, MAX_KEYS>::psum() {
#else
  uint64_t BSkip<K, MAX_KEYS>::psum() {
#endif
	// std::vector<uint64_t> partial_sums(ParallelTools::getWorkers() * 8);
	std::vector<uint64_t> partial_sums(1);
	//  auto curr_node = headers[curr_max_height];
	int start_level = 2;
	auto top_node = headers[start_level];
	while(top_node->get_header() < std::numeric_limits<K>::max()) {
		BSkipNodeInternal<K, MAX_KEYS>* top_node_cast = (BSkipNodeInternal<K, MAX_KEYS>*)top_node;	
		// iterate over the 2nd level
		// cilk_for(int i = 0; i < top_node->num_elts; i++) {
		for(uint64_t i = 0; i < top_node->num_elts; i++) {
			auto curr_node = top_node_cast->get_child_at_rank(i);
			K end;
			if (i < top_node->num_elts - 1) {
				end = top_node->get_key_at_rank(i+1);
			} else {
				end = top_node->next->get_header();
			}

			while(curr_node->get_header() < end) {
				BSkipNodeInternal<K, MAX_KEYS>* curr_node_cast = (BSkipNodeInternal<K, MAX_KEYS>*)curr_node;
				// cilk_for(int j = 0; j < curr_node->num_elts; j++) {
					for(uint64_t j = 0; j < curr_node->num_elts; j++) {
					auto curr_leaf = curr_node_cast->get_child_at_rank(j);
					K leaf_end = 0;
					if ( j < curr_node->num_elts - 1 ) {
						leaf_end = curr_node->get_key_at_rank(j+1);
					} else {
						leaf_end = curr_node->next->get_header();
					}

					while(curr_leaf->get_header() < leaf_end) {
						// partial_sums[ParallelTools::getWorkerNum() * 8] += curr_leaf->sum_keys();
						partial_sums[0] += curr_leaf->sum_keys();
						curr_leaf = curr_leaf->next;
					}
				}

				curr_node = curr_node->next;
			}
		}
		top_node = top_node->next;
	}
		/*
		if (i < curr_node->num_elts - 1) {
			sum_helper(partial_sums, curr_node_cast->get_child_at_rank(i), 0, curr_node->get_key_at_rank(i+1));
		} else {
			sum_helper(partial_sums, curr_node_cast->get_child_at_rank(i), 0, curr_node->next->get_header());
		}
		*/
	
	// sum_helper(partial_sums, curr_node, start_level, std::numeric_limits<K>::max());

	// add up results
	uint64_t result = 0;
	// for(int i = 0; i < ParallelTools::getWorkers(); i++) {
	for(int i = 0; i < 1; i++) {
		result += partial_sums[i * 8];
	}
	return result;
}

#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS>
#else
template <typename K, size_t MAX_KEYS>
#endif
#if WEIGHTED
  void BSkip<K, V, MAX_KEYS>::get_size_stats() {
#else
  void BSkip<K, MAX_KEYS>::get_size_stats() {
#endif

	// we want the density per level
	uint64_t elts_per_level[MAX_HEIGHT];
	uint64_t nodes_per_level[MAX_HEIGHT];

	for(int level = MAX_HEIGHT - 1; level >= 0; level--) {
		// empty level if above the current max promoted elt
		if (level > curr_max_height) { 
			elts_per_level[level] = 2;
			nodes_per_level[level] = 2;
			continue;
		}

		// count nodes and actual elts in them
  	auto curr_node = headers[level];
		nodes_per_level[level] = 0;
		elts_per_level[level] = 0;
		while(curr_node->get_header() < std::numeric_limits<K>::max()) {
			nodes_per_level[level]++;
			elts_per_level[level] += curr_node->num_elts;
			curr_node = curr_node->next;
		}

		// add in last sentinel
		nodes_per_level[level]++;
		elts_per_level[level]++;
	}

	uint64_t total_elts = 0;
	uint64_t total_nodes = 0;
	double density_per_level[MAX_HEIGHT];
	uint64_t num_internal_nodes = 0;
	for(int level = MAX_HEIGHT - 1; level >= 0; level--) {
		// count up for total size and density
		total_elts += elts_per_level[level];
		total_nodes += nodes_per_level[level];

		if (level > 0) { num_internal_nodes += nodes_per_level[level]; }

		// get density at this level
		double density = (double)elts_per_level[level] / (double)(nodes_per_level[level] * MAX_KEYS);
		printf("level %d, elts %lu, nodes %lu, total slots %lu, density = %f\n", level, elts_per_level[level], nodes_per_level[level], nodes_per_level[level] * MAX_KEYS, density);
		
		density_per_level[level] = density;
	}

	// printf("size of internal %lu, size of leaf %lu\n", sizeof(BSkipNodeInternal<K, MAX_KEYS>), sizeof(BSkipNode<K, MAX_KEYS>));
	uint64_t internal_size = num_internal_nodes * sizeof(BSkipNodeInternal<K, MAX_KEYS>);
	uint64_t size_in_bytes = internal_size +  nodes_per_level[0] * sizeof(BSkipNode<K, MAX_KEYS>);
	double overhead = (double)internal_size / (double) size_in_bytes;
	double overall_density = (double)total_elts/(double)(total_nodes * MAX_KEYS);


	double leaf_avg = elts_per_level[0] / nodes_per_level[0];
	// min, max, var of leaves
  double var_numerator = 0;
	uint32_t min_leaf = std::numeric_limits<uint32_t>::max();
	uint32_t max_leaf = std::numeric_limits<uint32_t>::min();
	auto curr_node = headers[0];
	while(curr_node->get_header() < std::numeric_limits<K>::max()) {
		if (curr_node->num_elts < min_leaf) { min_leaf = curr_node->num_elts; }
		if (curr_node->num_elts > max_leaf) { max_leaf = curr_node->num_elts; }
		double x = (double)curr_node->num_elts - leaf_avg;
		var_numerator += (x * x);
		curr_node = curr_node->next;
	}
	double var = var_numerator / (double)(nodes_per_level[0]); 
	double stdev = sqrt(var);

	printf("avg %f, min %u, max %u, var %f, stddev %f\n", leaf_avg, min_leaf, max_leaf, var, stdev);
  FILE* file = fopen("bskip_sizes.csv", "a+");
  fprintf(file, "%lu,%lu,%d,%lu,%lu,%f,%f,%lu,%f,%lu,%f,%lu,%f,%lu,%f,%lu,%f,%f,%u,%u,%f,%f\n", 
	elts_per_level[0], MAX_KEYS, curr_max_height, internal_size, size_in_bytes, overhead, overall_density,
	nodes_per_level[4], density_per_level[4], nodes_per_level[3], density_per_level[3], nodes_per_level[2], density_per_level[2],
	nodes_per_level[1], density_per_level[1], nodes_per_level[0], density_per_level[0], leaf_avg, min_leaf, max_leaf, var, stdev);
}


#if WEIGHTED
template <typename K, typename V, size_t MAX_KEYS>
#else
template <typename K, size_t MAX_KEYS>
#endif
#if WEIGHTED
  void BSkip<K, V, MAX_KEYS>::get_avg_comparisons() {
#else
  void BSkip<K, MAX_KEYS>::get_avg_comparisons() {
#endif
		uint64_t total_comparisons = 0;
		uint64_t num_nodes = 0;
		for(int level = curr_max_height; level >= 0; level--) {
			auto curr_node = headers[level];
			while(curr_node->get_header() < std::numeric_limits<K>::max()) {
				total_comparisons += curr_node->num_comparisons;
				num_nodes++;
				curr_node = curr_node->next;
			}
		}
		double avg_comparisons = (double) total_comparisons / (double) num_nodes;
		printf("avg comparisons = %f\n", avg_comparisons);
 }

#endif
