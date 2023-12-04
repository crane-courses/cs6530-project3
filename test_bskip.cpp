// test driver for bskip

#define DEBUG 0
#define DEBUG_PRINT 0
#define STATS 0

// #include "ParallelTools/reducer.h"
#include <assert.h>
#include <random>
#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <sys/time.h>
// #include <ParallelTools/parallel.h>

// #include "array.h"
//#include "parallel.h"
//#include "cilk/cilk_api.h"

#include "bskip.h"

#if CILK != 1
#define cilk_for for
#endif

static long get_usecs() {
    struct timeval st;
    gettimeofday(&st,NULL);
    return st.tv_sec*1000000 + st.tv_usec;
}

template <class T>
std::vector<T> create_random_data(size_t n, size_t max_val,
                                  std::seed_seq &seed) {

  std::mt19937_64 eng(seed); // a source of random data

  std::uniform_int_distribution<T> dist(0, max_val);
  std::vector<T> v(n);

  generate(begin(v), end(v), bind(dist, eng));
  return v;
}

template <class T, size_t MAX_KEYS> void test_ordered_insert(uint64_t max_size, int p) {
  if (max_size > std::numeric_limits<T>::max()) {
    max_size = std::numeric_limits<T>::max();
  }
  uint64_t start, end;
  
#if WEIGHTED
  BSkip<T, T, MAX_KEYS> s(p);
#else
  BSkip<T, MAX_KEYS> s(p);
#endif
  
  start = get_usecs();
  for (uint32_t i = 1; i < max_size; i++) {
    s.insert(i);
  }
  end = get_usecs();
  printf("\ninsertion,\t %lu,", end - start);

  start = get_usecs();
  for (uint32_t i = 1; i < max_size; i++) {
    auto node = s.find(i);
    if (node == NULL) {
      printf("\ncouldn't find key %lu in skip at index %u\n", i, i);
      exit(0);
    }
  }
  end = get_usecs();
  printf("\nfind all,\t %lu,", end - start);

  start = get_usecs();
  uint64_t sum = s.sum();
  end = get_usecs();

  printf("\nsum_time, %lu, sum_total, %lu\n", end - start, sum);
	printf("\n");
}

template <class T, size_t MAX_KEYS>
void test_unordered_insert(uint64_t max_size, std::seed_seq &seed, double p, int boost) {
  if (max_size > std::numeric_limits<T>::max()) {
    max_size = std::numeric_limits<T>::max();
  }
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
      //create_random_data<T>(max_size, 100, seed);

  // std::set<T> inserted_data;

  uint64_t start, end;

  // save insertion, find, iter sum, naive sum times

#if WEIGHTED
  BSkip<T, T, MAX_KEYS> s(p, boost);
#else
  BSkip<T, MAX_KEYS> s(p, boost);
#endif
	printf("finished creating data\n");
  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    s.insert(data[i], i);
    // inserted_data.insert(data[i]);
  }
  end = get_usecs();
	uint64_t insert_time = end - start;
  printf("\ninsertion,\t %lu,", end - start);

	// serial find
  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    auto node = s.find(data[i]);
		
    if (node == nullptr) {
      printf("\ncouldn't find key %lu in skip at index %u\n", data[i], i);
      exit(0);
    }
  }
  end = get_usecs();
	uint64_t find_time = end - start;
  printf("\nfind all,\t %lu\n", end - start);

	// parallel find
  std::seed_seq seed2{1};

	// generate n / 10 random elts
  std::vector<T> data_to_search =
      create_random_data<T>(max_size / 10, std::numeric_limits<T>::max(), seed2);

	// pick n/10 from the input
	for(uint32_t i = 0; i < max_size; i+=10) {
		if (i < max_size) { data_to_search.push_back(data[i]); }
	}

	// shuffle them
  std::mt19937_64 g(seed); // a source of random data
	std::shuffle(data_to_search.begin(), data_to_search.end(), g);

	// std::vector<T> partial_sums(ParallelTools::getWorkers() * 8);
	std::vector<T> partial_sums(1);		//disable parallelism for now

  // serial longest and average
	T result = 0;
	uint64_t serial_total_time = 0;
	uint64_t serial_longest_find = 0;
  for (uint32_t i = 0; i < data_to_search.size(); i++) {
		start = get_usecs();
    auto node = s.find(data_to_search[i]);
  	end = get_usecs();
		result += !(node == nullptr);
		serial_total_time += end - start;
		if (end - start > serial_longest_find) {
			serial_longest_find = end - start;
		}
  }
	double serial_average_find = (double)serial_total_time / (double) data_to_search.size();

  printf("\nserial find (B = %lu),\tlongest %lu,\taverage %f,\tnum found %lu\n", MAX_KEYS, serial_longest_find, serial_average_find, result);

  // span
  start = get_usecs();
  cilk_for (uint32_t i = 0; i < data_to_search.size(); i++) {
    auto node = s.find(data_to_search[i]);
		// partial_sums[ParallelTools::getWorkerNum() * 8] += !(node == nullptr);
		partial_sums[0] += !(node == nullptr);
  }
  end = get_usecs();

	uint64_t parallel_longest_find = end - start;
	// sum up results
	result = 0;
	// for(int i = 0; i < ParallelTools::getWorkers(); i++) {
	for(int i = 0; i < 1; i++) {
		// result += partial_sums[i * 8];
		result += partial_sums[i];
	}

  printf("\nparallel longest find (B = %lu),\t %lu,\tnum found %lu\n", MAX_KEYS, parallel_longest_find, result);

	// average
	std::fill(partial_sums.begin(), partial_sums.end(), 0);
	// std::vector<uint64_t> partial_times(ParallelTools::getWorkers() * 8);
	std::vector<uint64_t> partial_times(1);

   cilk_for (uint32_t i = 0; i < data_to_search.size(); i++) {
		start = get_usecs();
    auto node = s.find(data_to_search[i]);
		end = get_usecs();
		// partial_times[ParallelTools::getWorkerNum() * 8] += end - start;
		// partial_sums[ParallelTools::getWorkerNum() * 8] += !(node == nullptr);
		partial_times[0] += end - start;
		partial_sums[0] += !(node == nullptr);
  }

	// sum up results
	result = 0;
	// for(int i = 0; i < ParallelTools::getWorkers(); i++) {
	for(int i = 0; i < 1; i++) {
		// result += partial_sums[i * 8];
		result += partial_sums[i];
	}

	uint64_t total_time = 0;
	// for(int i = 0; i < ParallelTools::getWorkers(); i++) {
	for(int i = 0; i < 1; i++) {
		total_time += partial_times[i];
	}
	double parallel_average_find = (double)total_time / (double) data_to_search.size();

	printf("parallel avg find (B = %lu), \t %f, \tnum found %lu\n", MAX_KEYS, parallel_average_find, result);

  start = get_usecs();
  uint64_t sum = s.sum();
  end = get_usecs();

  printf("\nsum_time, %lu, sum_total, %lu\n", end - start, sum);
	uint64_t sum_time = end - start;

  start = get_usecs();
  uint64_t psum = s.psum();
  end = get_usecs();

  printf("psum_time, %lu, psum_total, %lu\n", end - start, psum);
	uint64_t psum_time = end - start;

	uint64_t correct_sum = 0;
	for(int i = 0; i < max_size; i++) {
		correct_sum += data[i];
	}

	tbassert(correct_sum == sum, "got sum %lu, should be %lu\n", sum, correct_sum);
	tbassert(correct_sum == psum, "got psum %lu, should be %lu\n", psum, correct_sum);

	FILE* file = fopen("bskip_times.csv", "a+");
	fprintf(file, "%d,%lu,%lu,%lu,%lu,%lu,%lu,%f,%lu,%f\n", p, max_size, insert_time, find_time, sum_time, psum_time, serial_longest_find, serial_average_find, parallel_longest_find, parallel_average_find);

	s.get_size_stats();
#if STATS
	// s.get_avg_comparisons();
	printf("\n");
#endif
}

int main(int argc, char** argv) {
  std::seed_seq seed{0};
  // int p = atoi(argv[1]);
	int n = atoi(argv[1]);
	// uint32_t s = atoi(argv[2]);
	// int n = 100000000; // 100M
  // __cilkrts_set_param("nworkers","16");
//   printf("num workers %d\n", ParallelTools::getWorkers());
  printf("num workers %d\n", 1);
	printf("TESTING LOG BOOSTING AND STANDARD: sqrt\n");
	// 2**5 -> 2**13
	constexpr uint32_t p_32 = (1 << 5);
	printf("\nexp size %u\n", p_32);
	// printf("\nscale: %u\n", s);
	// test_unordered_insert<uint64_t, p_32 / 2>(n, seed, p_32);

	double s = 1;

	// test_unordered_insert<uint64_t, p_32>(n, seed, p_32 * s, 7);
	// test_unordered_insert<uint64_t, p_32>(n, seed, p_32 * s, 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 7);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 2);
	test_unordered_insert<uint64_t, p_32>(n, seed, std::sqrt(p_32 * s), 7);
	// test_unordered_insert<uint64_t, p_32>(n, seed, std::sqrt(p_32 * s), 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, std::cbrt(p_32 * s), 7);
	// test_unordered_insert<uint64_t, p_32>(n, seed, std::cbrt(p_32 * s), 2);


	// test_unordered_insert<uint64_t, p_32>(n, seed, p_32 * s, 1);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 1);
	test_unordered_insert<uint64_t, p_32>(n, seed, std::sqrt(p_32 * s), 1);

	printf("\nexp size 64\n");
	constexpr uint32_t p_64 = 1 << 6;
	// test_unordered_insert<uint64_t, p_64>(n, seed, p_64 * s, 9);
	// test_unordered_insert<uint64_t, p_64>(n, seed, p_64 * s, 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 9);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 2);
	test_unordered_insert<uint64_t, p_64>(n, seed, std::sqrt(p_64 * s), 9);
	// test_unordered_insert<uint64_t, p_64>(n, seed, std::sqrt(p_64 * s), 2);
	// test_unordered_insert<uint64_t, p_64>(n, seed, std::cbrt(p_64 * s), 9);
	// test_unordered_insert<uint64_t, p_64>(n, seed, std::cbrt(p_64 * s), 2);


	// test_unordered_insert<uint64_t, p_64>(n, seed, p_64 * s, 1);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 1);
	test_unordered_insert<uint64_t, p_64>(n, seed, std::sqrt(p_64 * s), 1);


	printf("\nexp size 128\n");
	constexpr uint32_t p_128 = (1 << 7);
	// test_unordered_insert<uint64_t, p_128>(n, seed, p_128 * s, 10);
	// test_unordered_insert<uint64_t, p_128>(n, seed, p_128 * s, 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 10);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 2);
	test_unordered_insert<uint64_t, p_128>(n, seed, std::sqrt(p_128 * s), 10);
	// test_unordered_insert<uint64_t, p_128>(n, seed, std::sqrt(p_128 * s), 2);
	// test_unordered_insert<uint64_t, p_128>(n, seed, std::cbrt(p_128 * s), 10);
	// test_unordered_insert<uint64_t, p_128>(n, seed, std::cbrt(p_128 * s), 2);


	// test_unordered_insert<uint64_t, p_128>(n, seed, p_128 * s, 1);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 1);
	test_unordered_insert<uint64_t, p_128>(n, seed, std::sqrt(p_128 * s), 1);

	printf("\nexp size 256\n");
	constexpr uint32_t p_256 = 1 << 8;
	// test_unordered_insert<uint64_t, p_256>(n, seed, p_256 * s, 12);
	// test_unordered_insert<uint64_t, p_256>(n, seed, p_256 * s, 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 12);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 2);
	test_unordered_insert<uint64_t, p_256>(n, seed, std::sqrt(p_256 * s), 12);
	// test_unordered_insert<uint64_t, p_256>(n, seed, std::sqrt(p_256 * s), 2);
	// test_unordered_insert<uint64_t, p_256>(n, seed, std::cbrt(p_256 * s), 12);
	// test_unordered_insert<uint64_t, p_256>(n, seed, std::cbrt(p_256 * s), 2);


	// test_unordered_insert<uint64_t, p_256>(n, seed, p_256 * s, 1);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 1);
	test_unordered_insert<uint64_t, p_256>(n, seed, std::sqrt(p_256 * s), 1);

	printf("\nexp size 512\n");
	constexpr uint32_t p_512 = 1 << 9;
	// test_unordered_insert<uint64_t, p_512>(n, seed, p_512 * s, 13);
	// test_unordered_insert<uint64_t, p_512>(n, seed, p_512 * s, 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 13);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 2);
	test_unordered_insert<uint64_t, p_512>(n, seed, std::sqrt(p_512 * s), 13);
	// test_unordered_insert<uint64_t, p_512>(n, seed, std::sqrt(p_512 * s), 2);
	// test_unordered_insert<uint64_t, p_512>(n, seed, std::cbrt(p_512 * s), 13);
	// test_unordered_insert<uint64_t, p_512>(n, seed, std::cbrt(p_512 * s), 2);


	// test_unordered_insert<uint64_t, p_512>(n, seed, p_512 * s, 1);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 1);
	test_unordered_insert<uint64_t, p_512>(n, seed, std::sqrt(p_512 * s), 1);

	printf("\nexp size 1024\n");
	constexpr uint32_t p_1024 = 1 << 10;
	// test_unordered_insert<uint64_t, p_1024>(n, seed, p_1024 * s, 14);
	// test_unordered_insert<uint64_t, p_1024>(n, seed, p_1024 * s, 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 14);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 2);
	test_unordered_insert<uint64_t, p_1024>(n, seed, std::sqrt(p_1024 * s), 14);
	// test_unordered_insert<uint64_t, p_1024>(n, seed, std::sqrt(p_1024 * s), 2);
	// test_unordered_insert<uint64_t, p_1024>(n, seed, std::cbrt(p_1024 * s), 14);
	// test_unordered_insert<uint64_t, p_1024>(n, seed, std::cbrt(p_1024 * s), 2);


	// test_unordered_insert<uint64_t, p_1024>(n, seed, p_1024 * s, 1);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 1);
	test_unordered_insert<uint64_t, p_1024>(n, seed, std::sqrt(p_1024 * s), 1);


	printf("\nexp size 2048\n");
	constexpr uint32_t p_2048 = 1 << 11;
	// test_unordered_insert<uint64_t, p_2048>(n, seed, p_2048 * s, 16);
	// test_unordered_insert<uint64_t, p_2048>(n, seed, p_2048 * s, 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 15);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 2);
	test_unordered_insert<uint64_t, p_2048>(n, seed, std::sqrt(p_2048 * s), 16);
	// test_unordered_insert<uint64_t, p_2048>(n, seed, std::sqrt(p_2048 * s), 2);
	// test_unordered_insert<uint64_t, p_2048>(n, seed, std::cbrt(p_2048 * s), 16);
	// test_unordered_insert<uint64_t, p_2048>(n, seed, std::cbrt(p_2048 * s), 2);


	// test_unordered_insert<uint64_t, p_2048>(n, seed, p_2048 * s, 1);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 1);
	test_unordered_insert<uint64_t, p_2048>(n, seed, std::sqrt(p_2048 * s), 1);


	printf("\nexp size 4096\n");
	constexpr uint32_t p_4096 = 1 << 12;
	// test_unordered_insert<uint64_t, p_4096 >(n, seed, p_4096 * s, 17);
	// test_unordered_insert<uint64_t, p_4096 >(n, seed, p_4096 * s, 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 17);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 2);
	test_unordered_insert<uint64_t, p_4096 >(n, seed, std::sqrt(p_4096 * s), 17);
	// test_unordered_insert<uint64_t, p_4096 >(n, seed, std::sqrt(p_4096 * s), 2);
	// test_unordered_insert<uint64_t, p_4096 >(n, seed, std::cbrt(p_4096 * s), 17);
	// test_unordered_insert<uint64_t, p_4096 >(n, seed, std::cbrt(p_4096 * s), 2);


	// test_unordered_insert<uint64_t, p_4096 >(n, seed, p_4096 * s, 1);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 1);
	test_unordered_insert<uint64_t, p_4096 >(n, seed, std::sqrt(p_4096 * s), 1);


	printf("\nexp size 8192\n");
	constexpr uint32_t p_8192 = 1 << 13;
	// test_unordered_insert<uint64_t, p_8192>(n, seed, p_8192 * s, 19);
	// test_unordered_insert<uint64_t, p_8192>(n, seed, p_8192 * s, 2);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 19);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 2);
	test_unordered_insert<uint64_t, p_8192>(n, seed, std::sqrt(p_8192 * s), 19);
	// test_unordered_insert<uint64_t, p_8192>(n, seed, std::sqrt(p_8192 * s), 2);
	// test_unordered_insert<uint64_t, p_8192>(n, seed, std::cbrt(p_8192 * s), 19);
	// test_unordered_insert<uint64_t, p_8192>(n, seed, std::cbrt(p_8192 * s), 2);


	// test_unordered_insert<uint64_t, p_8192>(n, seed, p_8192 * s, 1);
	// test_unordered_insert<uint64_t, p_32>(n, seed, 2, 1);
	test_unordered_insert<uint64_t, p_8192>(n, seed, std::sqrt(p_8192 * s), 1);


	/*
	constexpr int MAX_KEYS = 4096;
	// printf("------- ORDERED INSERT --------\n");
  // test_ordered_insert<uint64_t>(n, p);
  printf("------- UNORDERED INSERT --------\n");
  test_unordered_insert<uint64_t, MAX_KEYS>(n, seed, p);
	*/
}

