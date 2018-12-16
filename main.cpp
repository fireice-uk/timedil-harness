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
 *
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
#include <assert.h>

#include "misc.hpp"
#include "functors.hpp"

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

	while(true)
	{
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


// LWMA-3 difficulty algorithm
// Copyright (c) 2017-2018 Zawy, MIT License
// https://github.com/zawy12/difficulty-algorithms/issues/3
// See commented version for required config file changes. Fix your FTL and MTP.

// difficulty_type should be uint64_t
difficulty_type next_difficulty_v3_1(std::vector<uint64_t> timestamps,
    std::vector<difficulty_type> cumulative_difficulties) {

    uint64_t  T = 120;
    uint64_t  N = 60; // N=45, 60, and 90 for T=600, 120, 60.
    uint64_t  L(0), sum_3_ST(0), next_D, prev_D;
    int64_t ST, previous_timestamp;

    // If it's a new coin, do startup code.
    // Increase difficulty_guess if it needs to be much higher, but guess lower than lowest guess.
    uint64_t difficulty_guess = 100;
    if (timestamps.size() <= 10 ) {   return difficulty_guess;   }
    if ( timestamps.size() < N +1 ) { N = timestamps.size()-1;  }

    // If hashrate/difficulty ratio after a fork is < 1/3 prior ratio, hardcode D for N+1 blocks after fork.
    // difficulty_guess = 100; //  Dev may change.  Guess low.
    // if (height <= UPGRADE_HEIGHT + N+1 ) { return difficulty_guess;  }

    previous_timestamp = static_cast<int64_t>(timestamps[0]);
    for ( uint64_t i = 1; i <= N; i++) {
	   ST = static_cast<int64_t>(timestamps[i]) - previous_timestamp;
       ST = std::max(1l, std::min(ST, static_cast<int64_t>(6*T)));
	   previous_timestamp += ST;

       L +=  ST * i ;
	   //std::cout << "ST : " << ST << std::endl;
       // delete the following line if you do not want the "jump rule"
       if ( i > N-3 ) { sum_3_ST += ST; }
    }

    std::cout << "avg L " <<  L/((N*N*1+N*1)/2) << std::endl;
    next_D = ((cumulative_difficulties[N] - cumulative_difficulties[0])*T*(N+1)*99)/(100*2*L);
    prev_D = cumulative_difficulties[N] - cumulative_difficulties[N-1];
    next_D = std::max((prev_D*67)/100, std::min(next_D, (prev_D*150)/100));

    // delete the following line if you do not want the "jump rule"
    if ( sum_3_ST < (8*T)/10) {  next_D = std::max(next_D,(prev_D*108)/100); }

    return next_D;
}

template <typename T>
inline T clamp(T lo, T v, T hi)
{
	return v < lo ? lo : v > hi ? hi : v;
}

constexpr uint64_t T = 240;
constexpr uint64_t N = 45;

difficulty_type next_difficulty_v4(std::vector<uint64_t> timestamps, std::vector<difficulty_type> cumulative_difficulties)
{
	if(timestamps.size() != N + 1 || cumulative_difficulties.size() != N + 1)
		abort();

	uint64_t L = 0;
	uint64_t prev_t = timestamps[0];
	for(uint64_t i = 1; i <= N; i++)
	{
		uint64_t this_t = std::max(timestamps[i], prev_t);
		L += std::min(this_t - prev_t, 6*T) * i * i;
		prev_t = this_t;
	}

	// Let's take CD as a sum of N difficulties. Sum of weights is (n*(n+1)*(2n+1))/6 (SUM)
	// L is a sigma(timeperiods * weights)
	// D = CD*T*SUM / NL
	// D = CD*T*N*(N+1)*(2N+1) / 6NL
	// D = CD*T*(N+1)*(2N+1) / 6L
	// TSUM = T*(N+1)*(2N+1) / 6 (const)
	// D = CD*TSUM / L

	// By a happy accident most time units are a multiple of 6 so we can prepare a TSUM without loosing accuracy
	constexpr uint64_t TSUM = (T * (N+1) * (2*N+1)) / 6;
	uint64_t next_D = ((cumulative_difficulties[N] - cumulative_difficulties[0]) * TSUM) / L;

	// 99/100 adds a small bias towards decreasing diff, unlike zawy we do it in a separate step to avoid an overflow at 6GH/s
	next_D = (next_D * 99ull) / 100ull;

	// Sanity limits
	uint64_t prev_D = cumulative_difficulties[N] - cumulative_difficulties[N-1];
    next_D = std::max((prev_D*67)/100, std::min(next_D, (prev_D*150)/100));

	return next_D;
}

difficulty_type next_difficulty_v3(std::vector<uint64_t> timestamps, std::vector<difficulty_type> cumulative_difficulties)
{
	if(timestamps.size() != N + 1 || cumulative_difficulties.size() != N + 1)
		abort();

	int64_t FTL = DIFFICULTY_TARGET * 3;
	int64_t L = 0;
	for(int64_t i = 1; i <= N; i++)
		L += clamp(-FTL, int64_t(timestamps[i]) - int64_t(timestamps[i - 1]), int64_t(6 * T)) * i;

	constexpr int64_t clamp_increase = (T * N * (N + 1) * 99) / int64_t(100.0 * 2.0 * 2.5);
	constexpr int64_t clamp_decrease = (T * N * (N + 1) * 99) / int64_t(100.0 * 2.0 * 0.2);
	L = clamp(clamp_increase, L, clamp_decrease); // This guarantees positive L

	// Commetary by fireice
	// Let's take CD as a sum of N difficulties. Sum of weights is (n*(n+1)*(2n+1))/6
	// L is a sigma(timeperiods * weights)
	// Therefore D = CD*0.5*(2N^2+2N)*T / NL
	// Therefore D = 0.5*CD*(2N+2)*T / L
	// Therefore D = CD*(2N+2)*T / 2L
	uint64_t next_D = ((cumulative_difficulties[N] - cumulative_difficulties[0]) * T * (N+1)) / (2 * L);

	// 99/100 adds a small bias towards decreasing diff, unlike zawy we do it in a separate step to avoid an overflow at 6GH/s
	next_D = (next_D * 99ull) / 100ull;
	return next_D;
}

difficulty_type next_difficulty_v4_interp6(std::vector<uint64_t> timestamps, std::vector<difficulty_type> cumulative_difficulties)
{
	if(timestamps.size() != N + 1 || cumulative_difficulties.size() != N + 1)
		abort();

	// highest timestamp must higher than the lowest
	timestamps[N] = std::max(timestamps[N], timestamps[0] + 1);
	// the newest timestamp is not allwed to be older than the previous
	uint64_t t_last = std::max(timestamps[N], timestamps[N - 1]);
	// the newest timestamp can only be 5 target times in the future
	timestamps[N] = std::min(t_last, timestamps[N - 1] + 5 * T);

	uint64_t lastValid = timestamps[0];
	uint64_t maxValid = timestamps[N];
	for(size_t i = 1; i < N; i++)
	{
		/* 
		 * Mask timestamp if it is smaller or equal to last valid timestamp
		 * or if it is larger or equal to largest timestamp
		 */
		if(timestamps[i] <= lastValid || timestamps[i] >= maxValid)
			timestamps[i] = 0;
		else
			lastValid = timestamps[i];
	}

	for(size_t i = 2; i < N; i++)
	{
		/* 
		 * Mask timestamp if the next timestamp is masked
		 */
		if(timestamps[i] == 0)
			timestamps[i-1] = 0;
	}

	// Now replace zeros with number of masked timestamps before this one (inclusive)
	uint64_t mctr = 0;
	for(size_t i = 1; i < N; i++)
	{
		if(timestamps[i] == 0)
			timestamps[i] = ++mctr;
		else
			mctr = 0;
	}

	// interpolate timestamps of masked times
	for(uint64_t i = N-1; i > 0; i--)
	{
		if(timestamps[i] <= N)
		{
			// denominator -- NOT THE SAME AS [i+1]
			uint64_t den = timestamps[i] + 1; 
			// numerator
			uint64_t num = timestamps[i];
			uint64_t delta = timestamps[i+1] - timestamps[i-num];
			timestamps[i] = timestamps[i-num] + (delta * num) / den;
		}
	}

	uint64_t L = 0;
	for(uint64_t i = 1; i <= N; i++)
	{
		assert(timestamps[i] > timestamps[i - 1]);
		assert(timestamps[i] != 0);
		L += std::min(timestamps[i] - timestamps[i - 1], 5 * T) * i * i;
	}
	L = std::max(L, 1ul);

	// Let's take CD as a sum of N difficulties. Sum of weights is (n*(n+1)*(2n+1))/6 (SUM)
	// L is a sigma(timeperiods * weights)
	// D = CD*T*SUM / NL
	// D = CD*T*N*(N+1)*(2N+1) / 6NL
	// D = CD*T*(N+1)*(2N+1) / 6L
	// TSUM = T*(N+1)*(2N+1) / 6 (const)
	// D = CD*TSUM / L

	// By a happy accident most time units are a multiple of 6 so we can prepare a TSUM without loosing accuracy
	constexpr uint64_t TSUM = (T * (N + 1) * (2 * N + 1)) / 6;
	uint64_t next_D = ((cumulative_difficulties[N] - cumulative_difficulties[0]) * TSUM) / L;

	// Sanity limits
	uint64_t prev_D = cumulative_difficulties[N] - cumulative_difficulties[N - 1];
	next_D = std::max((prev_D * 67) / 100, std::min(next_D, (prev_D * 150) / 100));

	return next_D;
}

difficulty_type next_difficulty_v4_interp2(std::vector<uint64_t> timestamps, std::vector<difficulty_type> cumulative_difficulties)
{
        if(timestamps.size() != N + 1 || cumulative_difficulties.size() != N + 1)
                abort();

		// the newest timestamp is not allwed to be older than the previous
		uint64_t t_last = std::max(timestamps[N], timestamps[N - 1] - 0 * T);
		// the newest timestamp can only be 3 target times in the future
        timestamps[N] = std::min(t_last, timestamps[N - 1] + 2 * T);

		// mask all invalid timestamps
        std::vector<bool> mask(timestamps.size());

        uint64_t lastValid=timestamps[0];
	    uint64_t maxValid=timestamps[N];
        for(uint64_t i = 1; i < N; i++)
        {
			/* check for invalid time stamps
			 *
			 * invalid is each timestamp which is older than the previous
			 * or newer than the latest block timestamp
			 *
			 * If a timestamp is invalid the current and the last timestamp is masked as invalid
			 */
			if(timestamps[i] < lastValid || timestamps[i] > maxValid )
			{
				// do not invalidate the first timestamp
				if(i-1 != 0)
				{
					mask[i-1] = true;
				}
				// do not invalidate the last timestamp
				if(i != N)
				{
					mask[i] = true;
				}
			}
			else
				lastValid=timestamps[i];
        }


		// interpolate timestamps of masked times
        for(uint64_t i = 1; i < N; i++)
        {
			if(mask[i])
			{
				uint64_t x = i + 1;
				for(; x < N; ++x)
				{
					if(!mask[x])
					{
						break;
					}
				}
				std::cerr<< "Interp" << timestamps[i] << " to " <<timestamps[i - 1] + (timestamps[x] - timestamps[i - 1]) / (x - i + 1) << std::endl;
				timestamps[i] = timestamps[i - 1] + (timestamps[x] - timestamps[i - 1]) / (x - i + 1);
			}
        }

        uint64_t L = 0;

        for(uint64_t i = 1; i <= N; i++)
        {
			L+= (timestamps[i] - timestamps[i -1]) * i * i;
        }

        // Let's take CD as a sum of N difficulties. Sum of weights is (n*(n+1)*(2n+1))/6 (SUM)
        // L is a sigma(timeperiods * weights)
        // D = CD*T*SUM / NL
        // D = CD*T*N*(N+1)*(2N+1) / 6NL
        // D = CD*T*(N+1)*(2N+1) / 6L
        // TSUM = T*(N+1)*(2N+1) / 6 (const)
        // D = CD*TSUM / L

        // By a happy accident most time units are a multiple of 6 so we can prepare a TSUM without loosing accuracy
        constexpr uint64_t TSUM = (T * (N+1) * (2*N+1)) / 6;
        uint64_t next_D = ((cumulative_difficulties[N] - cumulative_difficulties[0]) * TSUM) / L;

        // Sanity limits
        uint64_t prev_D = cumulative_difficulties[N] - cumulative_difficulties[N-1];
    	next_D = std::max((prev_D*80)/100, std::min(next_D, (prev_D*120)/100));

        return next_D;
}

int main(int argc, char **argv)
{
    FunctorConnector::getInstance().addRules(argv[2]);

	std::vector<HashPowerFunction*> hashPowerFunctors = {
        new functors::hrSet{},
        new functors::hrMul
    };

    std::vector<TimestampFunction*> timestampFunctors = {
        new functors::addScaled{},
        new functors::set{},
		new functors::add{}
    };

	std::vector<uint64_t> timestamps;
	std::vector<uint64_t> cum_diffs;

#if 1
	uint64_t timestamp = 1500000000;
	uint64_t diff = 0;
	for(size_t i=0; i <= N; i++)
	{
		diff += 1000000 * T;
		timestamp += T;
		timestamps.emplace_back(timestamp);

		cum_diffs.emplace_back(diff);
	}

    size_t max_rounds = std::stoull(argv[1]);

    uint64_t hr = 1000000;
	uint64_t firstHr = hr;
	uint64_t firsDifficulty = 1;

	for(size_t i = 1; i < max_rounds; i++)
	{
		uint64_t d;
		// select diff algorithm
		if(argc == 4 && std::string(argv[3]) == "v3")
		{
			if(i == 1)
				std::cerr<<"algorithm v3"<<std::endl;
            d = next_difficulty_v3(timestamps, cum_diffs);
		}
		else if(argc == 4 && std::string(argv[3]) == "v3_1")
		{
			if(i == 1)
				std::cerr<<"algorithm v3_1"<<std::endl;
			d = next_difficulty_v3_1(timestamps, cum_diffs);
		}
		else if(argc == 4 && std::string(argv[3]) == "v4-interp6")
		{
			if(i == 1)
				std::cerr<<"algorithm v4-interp6"<<std::endl;
			d = next_difficulty_v4_interp6(timestamps, cum_diffs);
		}
		else if(argc == 4 && std::string(argv[3]) == "v4-interp2")
		{
			if(i == 1)
				std::cerr<<"algorithm v4-interp2"<<std::endl;
			d = next_difficulty_v4_interp2(timestamps, cum_diffs);
		}
		else
		{
			if(i == 1)
				std::cerr<<"algorithm v4"<<std::endl;
			d = next_difficulty_v4(timestamps, cum_diffs);
		}
		diff += d;
		if(i == 1)
			firsDifficulty = d;

		// adjust hash rate via functor chain
        FunctorConnector::getInstance().hashPower(hashPowerFunctors, i, hr);
		uint64_t solve_time =  d / hr;

        timestamp += solve_time;
        uint64_t fakeTimestamp = timestamp;

		// manipulate block time stamp via functor chain
        FunctorConnector::getInstance().timeStamp(timestampFunctors, i, fakeTimestamp, solve_time);

        std::cout << i <<
			" timestamp " << timestamp << " fake_timestamp " << fakeTimestamp <<
			" diff " << d << " normalized_diff "<< (double)d/(double)firsDifficulty <<
			" solve_time " << solve_time << " hash_rate "<<hr<<" normalized_hash_rate "<<(double)hr/(double)firstHr<<std::endl;

		cum_diffs.emplace_back(diff);
	    timestamps.emplace_back(fakeTimestamp);

		cum_diffs.erase(cum_diffs.begin());
		timestamps.erase(timestamps.begin());
	}

#else

	base_walltime = get_walltime();
	base_timestamp = std::chrono::steady_clock::now();

	block_diff = 1;

	std::thread thd_1(honest_miner);
	std::thread athd_1;

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

		// we need N + 1 blcoks for v4 difficulty algorithm
		if(block > N + 1)
		{
			cum_diffs.erase(cum_diffs.begin());
			timestamps.erase(timestamps.begin());
		}

		int64_t time_dil = dilated_time();
		strftime(diltime, sizeof(diltime), "%X", gmtime(&time_dil));
		int64_t time_wall = get_walltime();
		strftime(wtime, sizeof(wtime), "%X", gmtime(&time_wall));
		snprintf(tbuf, sizeof(tbuf), "[ %s | %s ] : ", wtime, diltime);

		std::cout << tbuf << "block (" << block << ") found by " << (blk.honest ? "honest" : "attacker") << " diff: " <<
			std::setw(8) << std::setfill(' ') << blk.diff << " timestamp: " << blk.timestamp << "\n\n";

		if(block <= N + 1)
			block_diff = difficulty_const(timestamps, cum_diffs, DIFFICULTY_TARGET);
		else
		{
			// select diff algorithm
			if(argc == 4 && std::string(argv[3]) == "v3")
			{
				block_diff = next_difficulty_v3(timestamps, cum_diffs);
			}
			else if(argc == 4 && std::string(argv[3]) == "v3_1")
			{
				block_diff = next_difficulty_v3_1(timestamps, cum_diffs);
			}
			else
			{
				block_diff = next_difficulty_v4(timestamps, cum_diffs);
			}
		}

		std::cout << "\n" << tbuf << "next diff is " << block_diff << "\n" << "We found " << block << " blocks with average window of " << elapsed_time() / block << "s\n";

		if(block == ATTACK_START_BLOCK)
		{
			std::cout << "!!! Attack miner starting!" << "\n";
			athd_1 = std::thread(attack_miner);
		}
	}
#endif

}
