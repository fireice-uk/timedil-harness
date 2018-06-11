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

constexpr size_t   TIME_DILATION_MULT = 1000;
constexpr uint64_t ATTACK_START_BLOCK = 150;

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
	//uint64_t stamp = dilated_time();

	while(true)
	{
		for(size_t i=0; i < TIME_DILATION_MULT; i++)
		{
			uint64_t diff = 0xFFFFFFFFFFFFFFFFULL / hash(gen);
			if(diff > block_diff)
			{
				block_found blk;
				blk.diff = diff;
				blk.timestamp = dilated_time() + 3600*2;
				blk.honest = false;
				blk_q.push(blk);
				break;
			}
		}
	}

	std::this_thread::sleep_for(std::chrono::milliseconds(1));
}


// LWMA-2 difficulty algorithm (commented version)
// Copyright (c) 2017-2018 Zawy
// https://github.com/zawy12/difficulty-algorithms/issues/3
// Bitcoin clones must lower their FTL. 
// Cryptonote et al coins must make the following changes:
#define DIFFICULTY_TARGET 240
#define BLOCKCHAIN_TIMESTAMP_CHECK_WINDOW    11
#define CRYPTONOTE_BLOCK_FUTURE_TIME_LIMIT        3 * DIFFICULTY_TARGET 
#define DIFFICULTY_WINDOW                      (60 + 1)
typedef uint64_t difficulty_type;

// Do not sort timestamps.  CN coins must deploy the Jagerman MTP Patch. See:
// https://github.com/loki-project/loki/pull/26   or
// https://github.com/wownero/wownero/commit/1f6760533fcec0d84a6bd68369e0ea26716b01e7

// New coins:  "return 100;" can change the 100 based on expected hashrate.
// D should not be < 50 for the life of the coin because this uses integer math,

difficulty_type next_difficulty(std::vector<std::uint64_t> timestamps, std::vector<difficulty_type> cumulative_difficulties)
{
	int64_t T    = DIFFICULTY_TARGET; // target solvetime seconds
	int64_t N   = DIFFICULTY_WINDOW - 1; //  N=45, 60, and 90 for T=600, 120, 60.
	int64_t FTL = CRYPTONOTE_BLOCK_FUTURE_TIME_LIMIT; // < 3xT
	int64_t L(0), ST, sum_3_ST(0), next_D, prev_D, SMA;

	// If expecting a 10x decrease or 1000x increase in D after a fork, seriously consider: 
	// if ( height >= fork_height && height <= fork_height+N )  { return difficulty_guess; }

	// TS and CD vectors must be size N+1 after startup, and element N is most recent block.

	// If coin is starting, this will be activated.
	uint64_t initial_difficulty_guess = 240000; // Dev needs to select this. Guess low.
	if (timestamps.size() <= static_cast<uint64_t>(N) )  {  
		return initial_difficulty_guess;  
	}

	// N is most recently solved block. i must be signed
	for ( int64_t i = 1; i <= N; i++) {  
		// +/- FTL limits are bad timestamp protection.  6xT limits drop in D to reduce oscillations.
		ST = std::max(-FTL, std::min( (int64_t)(timestamps[i]) - (int64_t)(timestamps[i-1]), 6*T));
		L +=  ST * i ; // Give more weight to most recent blocks.
		// Do these inside loop to capture -FTL and +6*T limitations.
		if ( i > N-3 ) { sum_3_ST += ST; }      
	}
	// Calculate next_D = avgD * T / LWMA(STs) using integer math
	std::cout << "diff sum: " << (cumulative_difficulties[N] - cumulative_difficulties[0]) << " L " << L << " sizes " << timestamps.size() << " " << cumulative_difficulties.size() << "\n";
	next_D = ((cumulative_difficulties[N] - cumulative_difficulties[0])*T*(N+1)*99)/(100*2*L);

	std::cout << "next_D " << next_D << "\n";

	// begin LWMA-2 changes from LWMA
	prev_D = cumulative_difficulties[N] - cumulative_difficulties[N-1];
	SMA = ( cumulative_difficulties[N] - cumulative_difficulties[N-N/2] ) *4 *T / 
		( 3*(N-N/2)*T + (int64_t)(timestamps[N]) - (int64_t)(timestamps[N-N/2]) ); 

	std::cout << "SMA " << SMA << "\n";
	if ( sum_3_ST < (8*T)/10 && (prev_D*109)/100 < (12*SMA)/10) {  
		next_D = (prev_D*109)/100; 
	}
	std::cout << "next_D " << next_D << "\n";
	// end LWMA-2 section
	return static_cast<uint64_t>(next_D);
	// next_Target = sumTargets*L*2/0.998/T/(N+1)/N/N; // To show the difference.
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

	block_diff = 240000;

	std::thread thd_1(honest_miner);
	std::thread athd_1, athd_2;

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
		diff_sum += block_diff;//blk.diff;
		block++;
		timestamps.emplace_back(blk.timestamp);
		cum_diffs.emplace_back(diff_sum);
		
		if(timestamps.size() > DIFFICULTY_WINDOW)
			timestamps.erase(timestamps.begin());
		if(cum_diffs.size() > DIFFICULTY_WINDOW)
			cum_diffs.erase(cum_diffs.begin());

		int64_t time_dil = dilated_time();
		strftime(diltime, sizeof(diltime), "%X", gmtime(&time_dil));
		int64_t time_wall = get_walltime();
		strftime(wtime, sizeof(wtime), "%X", gmtime(&time_wall));
		snprintf(tbuf, sizeof(tbuf), "[ %s | %s ] : ", wtime, diltime);

		std::cout << tbuf << "block (" << block << ") found by " << (blk.honest ? "honest" : "attacker") << " diff: " << 
			std::setw(8) << std::setfill(' ') << blk.diff << " timestamp: " << blk.timestamp << "\n\n";

		//block_diff = next_difficulty_v2(timestamps, cum_diffs, DIFFICULTY_TARGET);
		block_diff = next_difficulty(timestamps, cum_diffs);
		std::cout << "\n" << tbuf << "next diff is " << block_diff << "\n" << "We found " << block << " blocks with average window of " << elapsed_time() / block << "s\n";

		if(block == ATTACK_START_BLOCK)
		{
			std::cout << "!!! Attack miner starting!" << "\n";
			athd_1 = std::thread(attack_miner);
			athd_2 = std::thread(attack_miner);
		}
	}
}
