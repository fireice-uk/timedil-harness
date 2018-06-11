/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "thdq.hpp"
#include <thread>
#include <atomic>
#include <chrono>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <random>
#include <iomanip>
#include <algorithm>

struct block_found
{
	uint64_t diff;
	uint64_t timestamp;
	bool honest;
};

constexpr size_t   TIME_DILATION_MULT = 100;
constexpr uint64_t ATTACK_START_BLOCK = 0;

thdq<block_found> blk_q;
std::atomic<uint64_t> block_diff;
uint64_t base_walltime;
std::chrono::time_point<std::chrono::steady_clock> base_timestamp;

inline int64_t elapsed_time()
{
	using namespace std::chrono;
	double elapsed_time = duration_cast<milliseconds>(steady_clock::now() - base_timestamp).count();
	elapsed_time /= 1000.0 / TIME_DILATION_MULT;
	return elapsed_time;
}

inline int64_t dilated_time()
{
	return base_walltime + elapsed_time();
}

inline uint64_t hash(std::mt19937_64 &gen)
{
	std::uniform_int_distribution<uint64_t> dis;
	return dis(gen);
}

inline int64_t get_walltime()
{
	using namespace std::chrono;
	return duration_cast<seconds>(system_clock::now().time_since_epoch()).count();
}

// Each mining thread does 1000H/S (dilated)
void honest_miner()
{
	std::random_device rd;
	std::mt19937_64 gen(rd());

	while(true)
	{
		for(size_t i=0; i < TIME_DILATION_MULT; i++)
		{
			uint64_t diff = 0xFFFFFFFFFFFFFFFFULL / hash(gen);
			if(diff > block_diff)
			{
				block_found blk;
				blk.diff = diff;
				blk.timestamp = dilated_time();
				blk.honest = true;
				blk_q.push(blk);
				break;
			}
		}

		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	}
}

void attack_miner()
{
	std::random_device rd;
	std::mt19937_64 gen(rd());

	for(size_t i=0; i < TIME_DILATION_MULT; i++)
	{
		uint64_t diff = 0xFFFFFFFFFFFFFFFFULL / hash(gen);
		if(diff > block_diff)
		{
			block_found blk;
			blk.diff = diff;
			blk.timestamp = base_walltime;
			blk.honest = false;
			blk_q.push(blk);
			break;
		}
	}

	std::this_thread::sleep_for(std::chrono::milliseconds(1));
}

#define DIFFICULTY_TARGET                               240  // seconds
#define DIFFICULTY_WINDOW                               720  // blocks
#define DIFFICULTY_LAG                                  15   // !!!
#define DIFFICULTY_CUT                                  60   // timestamps to cut after sorting
#define DIFFICULTY_BLOCKS_COUNT                         DIFFICULTY_WINDOW + DIFFICULTY_LAG

#define DIFFICULTY_WINDOW_V2							              17
#define DIFFICULTY_CUT_V2                               6
#define DIFFICULTY_BLOCKS_COUNT_V2                      DIFFICULTY_WINDOW_V2 + DIFFICULTY_CUT_V2*2
#define MAX_AVERAGE_TIMESPAN          (uint64_t) DIFFICULTY_TARGET*6   // 24 minutes
#define MIN_AVERAGE_TIMESPAN          (uint64_t) DIFFICULTY_TARGET/24  // 10s


namespace epee
{
namespace misc_utils
{
  template<class type_vec_type>
  type_vec_type median(std::vector<type_vec_type> &v)
  {
    if(v.empty())
      return 0;
    if(v.size() == 1)
      return v[0];

    size_t n = (v.size()) / 2;
    std::sort(v.begin(), v.end());
    //nth_element(v.begin(), v.begin()+n-1, v.end());
    if(v.size()%2)
    {//1, 3, 5...
      return v[n];
    }else 
    {//2, 4, 6...
      return (v[n-1] + v[n])/2;
    }
  }
}
}

static inline void mul(uint64_t a, uint64_t b, uint64_t &low, uint64_t &high) {
	unsigned __int128 r = (unsigned __int128)a * (unsigned __int128)b;
	high = r >> 64;
	low = (uint64_t)r;
}

typedef uint64_t difficulty_type;

difficulty_type difficulty_sumo (std::vector<std::uint64_t> timestamps, std::vector<difficulty_type> cumulative_difficulties, size_t target_seconds) {
    if (timestamps.size() > DIFFICULTY_BLOCKS_COUNT_V2)
    {
      timestamps.resize(DIFFICULTY_BLOCKS_COUNT_V2);
      cumulative_difficulties.resize(DIFFICULTY_BLOCKS_COUNT_V2);
    }

    size_t length = timestamps.size();
    if (length <= 1) {
      return 1;
    }

    std::sort(timestamps.begin(), timestamps.end());
    size_t cut_begin, cut_end;
    static_assert(2 * DIFFICULTY_CUT_V2 <= DIFFICULTY_BLOCKS_COUNT_V2 - 2, "Cut length is too large");
    if (length <= DIFFICULTY_BLOCKS_COUNT_V2 - 2 * DIFFICULTY_CUT_V2) {
      cut_begin = 0;
      cut_end = length;
    }
    else {
      cut_begin = (length - (DIFFICULTY_BLOCKS_COUNT_V2 - 2 * DIFFICULTY_CUT_V2) + 1) / 2;
      cut_end = cut_begin + (DIFFICULTY_BLOCKS_COUNT_V2 - 2 * DIFFICULTY_CUT_V2);
    }
    
    uint64_t total_timespan = timestamps[cut_end - 1] - timestamps[cut_begin];
    if (total_timespan == 0) {
      total_timespan = 1;
    }

    uint64_t timespan_median = 0;
    if (cut_begin > 0 && length >= cut_begin * 2 + 3){
      std::vector<std::uint64_t> time_spans;
      for (size_t i = length - cut_begin * 2 - 3; i < length - 1; i++){
        uint64_t time_span = timestamps[i + 1] - timestamps[i];
        if (time_span == 0) {
          time_span = 1;
        }
        time_spans.push_back(time_span);

        std::cout << "Timespan " << i << ": " << (time_span / 60) / 60 << ":" << (time_span > 3600 ? (time_span % 3600) / 60 : time_span / 60) << ":" << time_span % 60 << " (" << time_span << ")\n";
      }
      timespan_median = epee::misc_utils::median(time_spans);
    }

    uint64_t timespan_length = length - cut_begin * 2 - 1;
    std::cout << "Timespan Median: " << timespan_median << ", Timespan Average: " << total_timespan / timespan_length << "\n";

    uint64_t total_timespan_median = timespan_median > 0 ? timespan_median * timespan_length : total_timespan * 7 / 10;
    uint64_t adjusted_total_timespan = (total_timespan * 8 + total_timespan_median * 3) / 10; //  0.8A + 0.3M (the median of a poisson distribution is 70% of the mean, so 0.25A = 0.25/0.7 = 0.285M)
    if (adjusted_total_timespan > MAX_AVERAGE_TIMESPAN * timespan_length){
      adjusted_total_timespan = MAX_AVERAGE_TIMESPAN * timespan_length;
    }
    if (adjusted_total_timespan < MIN_AVERAGE_TIMESPAN * timespan_length){
      adjusted_total_timespan = MIN_AVERAGE_TIMESPAN * timespan_length;
    }
    
    difficulty_type total_work = cumulative_difficulties[cut_end - 1] - cumulative_difficulties[cut_begin];

    uint64_t low, high;
    mul(total_work, target_seconds, low, high);
    if (high != 0) {
      return 0;
    }

    uint64_t next_diff = (low + adjusted_total_timespan - 1) / adjusted_total_timespan;
    if (next_diff < 1) next_diff = 1;
    std::cout << "Total timespan: " << total_timespan << ", Adjusted total timespan: " << adjusted_total_timespan << ", Total work: " << total_work << ", Next diff: " << next_diff << ", Hashrate (H/s): " << next_diff / target_seconds << "\n";

    return next_diff;
  }

//Const diff assuming a single miner
difficulty_type difficulty_const(std::vector<std::uint64_t> timestamps, std::vector<difficulty_type> cumulative_difficulties, size_t target_seconds) 
{
	return target_seconds * 1000;
}

int main(int argc, char **argv) 
{
	base_walltime = get_walltime();
	base_timestamp = std::chrono::steady_clock::now();

	block_diff = 1;

	std::thread thd_1(honest_miner);
	std::thread athd_1;

	std::vector<uint64_t> timestamps;
	std::vector<uint64_t> cum_diffs;
	uint64_t diff_sum = 0;
	uint64_t block = 0;

	while(true)
	{
		char wtime[64];
		char diltime[64];
		char tbuf[128];
		block_found blk = blk_q.pop();
		diff_sum += blk.diff;
		block++;
		timestamps.emplace_back(blk.timestamp);
		cum_diffs.emplace_back(diff_sum);

		int64_t time_dil = dilated_time();
		strftime(diltime, sizeof(diltime), "%X", gmtime(&time_dil));
		int64_t time_wall = get_walltime();
		strftime(wtime, sizeof(wtime), "%X", gmtime(&time_wall));
		snprintf(tbuf, sizeof(tbuf), "[ %s | %s ] : ", wtime, diltime);

		std::cout << tbuf << "block (" << block << ") found by " << (blk.honest ? "honest" : "attacker") << " diff: " << 
			std::setw(8) << std::setfill(' ') << blk.diff << " timestamp: " << blk.timestamp << "\n\n";

		//block_diff = next_difficulty_v2(timestamps, cum_diffs, DIFFICULTY_TARGET);
		block_diff = difficulty_const(timestamps, cum_diffs, DIFFICULTY_TARGET);
		std::cout << "\n" << tbuf << "next diff is " << block_diff << "\n" << "We found " << block << " blocks with average window of " << elapsed_time() / block << "s\n";

		if(block == ATTACK_START_BLOCK)
		{
			std::cout << "!!! Attack miner starting!" << "\n";
			athd_1 = std::thread(attack_miner);
		}
	}
}
