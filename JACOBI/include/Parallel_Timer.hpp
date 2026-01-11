/**
 * 
 * Use this as
 * {
 *      Parallel_Timer t("what we are timing")
 *      ... code to be timed
 * }
 * 
 * at the end of main put
 * Parallel_Timer::gather_timing_data(MPI_COMM_WORLD, 0);
 * 
 */

#ifndef PARALLEL_TIMER_H
#define PARALLEL_TIMER_H

#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <map>
#include <stdexcept>


//#include "timer_tools.hpp"


struct TimerData{
    int calls{0};
    size_t time{0}; // Time in milliseconds
};

std::map<std::string, TimerData> timing_table;

class Parallel_Timer{
public:
    //using time_units = std::chrono::milliseconds;
    using time_units = std::chrono::microseconds;
    
    Parallel_Timer(const std::string& name0) : name(name0){
        if (timing_table.find(name) == timing_table.end()) {
            timing_table[name] = TimerData();
        }
        
        // Increment the call count for this timing label
        timing_table[name].calls++;

        start_time = std::chrono::steady_clock::now();
    }

    ~Parallel_Timer(){
        auto end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<time_units>(end_time - start_time).count();
        
        // Update the cumulative time for this label
        timing_table[name].time += duration; 
    }

    static void gather_timing_data(MPI_Comm comm, int root);
    static void print_timing_results();
    
    //auxiliar functions 
    template <typename T>
    static T max_element(const std::vector<T>& v);
    template <typename T>
    static T min_element(const std::vector<T>& v);
    template <typename T>
    static T accumulate_elements(const std::vector<T>& v, T base_value);
    template <typename T>
    static double average(const std::vector<T>& v);
    

private:
    const std::string name;
    std::chrono::time_point<std::chrono::steady_clock> start_time;

};

/**
 * 
 * We use MPI_Type_create_struct to create a custom MPI datatype that matches the memory layout of TimerData.
 * Each process sends the lengths of its function names.
 * We gather these lengths on the root process to know how to reconstruct the names.
 * Concatenate all function names into a single std::vector<char> per process.
 * Use MPI_Gatherv to gather these character arrays on the root process.
 * The root process uses the gathered lengths and character arrays to reconstruct the function names.
 * 
 * We use MPI_Gatherv with the custom MPI_TIMERDATA type to gather timing data from all processes.
 * Once the root process has all the data, it can organize it by function name and calculate the statistics 
 * using max_element, min_element, and average functions from math_utils.hpp
 * 
 * 
 * 
 */
void Parallel_Timer::gather_timing_data(MPI_Comm comm, int root) {
    int world_rank, world_size;
    MPI_Comm_rank(comm, &world_rank);
    MPI_Comm_size(comm, &world_size);

    // Define MPI Datatype for TimerData
    MPI_Datatype MPI_TIMERDATA;
    int lengths[2] = {1, 1};
    MPI_Aint displacements[2];
    MPI_Datatype types[2] = {MPI_INT, MPI_UNSIGNED_LONG};

    // Calculate displacements
    displacements[0] = offsetof(TimerData, calls);
    displacements[1] = offsetof(TimerData, time);

    MPI_Type_create_struct(2, lengths, displacements, types, &MPI_TIMERDATA);
    MPI_Type_commit(&MPI_TIMERDATA);

    // Prepare local data
    std::vector<std::string> names;
    std::vector<TimerData> times;

    // Copy timing data to vectors
    for (const auto& entry : timing_table) {

        names.push_back(entry.first);
        times.push_back(entry.second);
    }

    // Prepare data for function names
    int local_num_functions = names.size();
    std::vector<int> name_lengths(local_num_functions);
    std::vector<char> name_chars;

    for (int i = 0; i < local_num_functions; ++i) {
        name_lengths[i] = names[i].size();
        name_chars.insert(name_chars.end(), names[i].begin(), names[i].end());
    }

    // 1. Send number of functions each process has
    std::vector<int> all_num_functions(world_size);
    MPI_Gather(&local_num_functions, 1, MPI_INT, all_num_functions.data(), 1, MPI_INT, root, comm);

    // 2. Send lengths of function names
    std::vector<int> all_name_lengths;
    std::vector<int> name_lengths_displacements;
    if (world_rank == root) {
        int total_name_lengths = accumulate_elements(all_num_functions, 0);
        all_name_lengths.resize(total_name_lengths);
        name_lengths_displacements.resize(world_size);
        int disp = 0;
        for (int i = 0; i < world_size; ++i) {
            name_lengths_displacements[i] = disp;
            disp += all_num_functions[i];
        }
    }

    MPI_Gatherv(name_lengths.data(), local_num_functions, MPI_INT,
                all_name_lengths.data(), all_num_functions.data(), name_lengths_displacements.data(), MPI_INT, root, comm);
    
    // 3. Send the actual names concatenated into a char vector
    int local_name_chars_size = name_chars.size();
    std::vector<int> all_name_chars_sizes(world_size);
    MPI_Gather(&local_name_chars_size, 1, MPI_INT, all_name_chars_sizes.data(), 1, MPI_INT, root, comm);

    std::vector<int> name_chars_displacements;
    std::vector<char> all_name_chars;
    if (world_rank == root) {
        int total_name_chars_size = accumulate_elements(all_name_chars_sizes, 0);
        all_name_chars.resize(total_name_chars_size);
        name_chars_displacements.resize(world_size);
        int disp = 0;
        for (int i = 0; i < world_size; ++i) {
            name_chars_displacements[i] = disp;
            disp += all_name_chars_sizes[i];
        }
    }

    MPI_Gatherv(name_chars.data(), local_name_chars_size, MPI_CHAR,
                all_name_chars.data(), all_name_chars_sizes.data(), name_chars_displacements.data(), MPI_CHAR, root, comm);


    // 4. Send TimerData using MPI_Type
    // First, gather counts
    std::vector<int> timerdata_displacements;
    std::vector<TimerData> all_times;
    if (world_rank == root) {
        int total_times = accumulate_elements(all_num_functions, 0);
        all_times.resize(total_times);
        timerdata_displacements.resize(world_size);
        int disp = 0;
        for (int i = 0; i < world_size; ++i) {
            timerdata_displacements[i] = disp;
            disp += all_num_functions[i];
        }
    }

    MPI_Gatherv(times.data(), local_num_functions, MPI_TIMERDATA,
                all_times.data(), all_num_functions.data(), timerdata_displacements.data(), MPI_TIMERDATA,
                root, comm);

    // Now, the root process reconstructs the names and times
    if (world_rank == root) {
        std::vector<std::string> all_names;
        size_t idx = 0; // Index into all_name_chars
        for (size_t i = 0; i < all_name_lengths.size(); ++i) {
            int name_len = all_name_lengths[i];
            std::string name(all_name_chars.begin() + idx, all_name_chars.begin() + idx + name_len);
            all_names.push_back(name);
            idx += name_len;
        }

        // Now, we can organize the timing data
        std::map<std::string, std::vector<size_t>> timing_stats;
        for (size_t i = 0; i < all_names.size(); ++i) {
            timing_stats[all_names[i]].push_back(all_times[i].time);
        }

        // Calculate and print the max, min, and average for each function
        std::cout << "Timing Statistics:\n";
        for (const auto& entry : timing_stats) {
            const std::string& function_name = entry.first;
            const std::vector<size_t>& times_vec = entry.second;

            size_t max_time = max_element(times_vec);
            size_t min_time = min_element(times_vec);
            double avg_time = average(times_vec);

            std::cout <<function_name<<" Avg: " <<avg_time << " micro_s"<<std::endl;
            /*
            std::cout << "Function: " << function_name
                      << " | Max: " << max_time << " ms"
                      << " | Min: " << min_time << " ms"
                      << " | Avg: " << avg_time << " ms" << std::endl;
            */
        }
    }

    // Clean up the MPI datatype
    MPI_Type_free(&MPI_TIMERDATA);
}


void Parallel_Timer::print_timing_results(){
    std::cout << "Timing results:\n";
    for (const auto &entry : timing_table) {
        std::cout << entry.first << " -> Total time: " << entry.second.time << " ms, Calls: " << entry.second.calls << std::endl;
    }
}

////////////////////////// obtain min, max average and more ////////////////

template <typename T>
T Parallel_Timer::max_element(const std::vector<T>& v) {
   if (v.empty()) {
     throw std::invalid_argument("Vector is empty");
   }
   T max = v[0];
   for (size_t i = 1; i < v.size(); i++) {
     if (v[i] > max) {
       max = v[i];
     }
   }
  return max;
 }
 
 template <typename T>
 T Parallel_Timer::min_element(const std::vector<T>& v) {
   if (v.empty()) {
     throw std::invalid_argument("Vector is empty");
   }
   T min = v[0];
   for (size_t i = 1; i < v.size(); i++) {
     if (v[i] < min) {
         min = v[i];
     }
   }
   return min;
 }
 
 template <typename T>
 T Parallel_Timer::accumulate_elements(const std::vector<T>& v, T base_value){
   if (v.empty()) {
     throw std::invalid_argument("Vector is empty");
   }
   T sum = base_value;
   for (size_t i = 0; i < v.size(); i++) {
     sum += v[i];
   }
   return sum;
 }
 
 template <typename T>
 double Parallel_Timer::average(const std::vector<T>& v) {
   return accumulate_elements<T>(v, 0) / static_cast<double>(v.size());
 }

////////////////////////////////////////////////////////////////////////////

#endif
