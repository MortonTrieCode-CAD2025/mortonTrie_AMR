/**
* @file
* @brief Main function.
* @note .
*/
#include "Obj_Manager.h"
#include <iomanip>
#include <memory>
#include <typeinfo>
#include <fstream>
#if(C_MAP_TYPE == 2)
#include <algorithm>
#include <numeric>
#include <list>
#include <forward_list>
#endif

Grid_Manager* Grid_Manager::pointer_me;
Lat_Manager* Lat_Manager::pointer_me;
IO_Manager* IO_Manager::pointer_me;
Solid_Manager* Solid_Manager::pointer_me;

// declarations
Grid_Manager gr_manager;
Solid_Manager solid_manager;
IO_Manager io_manager;
Lat_Manager lat_manager;

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#elif defined(__unix__) || defined(__linux__) || defined(__APPLE__)
#include <sys/resource.h>
#include <unistd.h>
#endif

#if (C_MAP_TYPE == 1)
// Memory measurement utility class
class MapMemoryMeasurer {
public:
	// Method 1: Estimate memory usage using bucket information
	template <typename K, typename V, typename Hash = std::hash<K>,
				typename Pred = std::equal_to<K>,
				typename Alloc = std::allocator<std::pair<const K, V>>>
	static size_t estimateMemoryUsage(const std::unordered_map<K, V, Hash, Pred, Alloc>& map) {
		// Basic container size
		size_t memoryUsage = sizeof(std::unordered_map<K, V, Hash, Pred, Alloc>);

		// Add bucket array size
		size_t bucketCount = map.bucket_count();
		memoryUsage += bucketCount * sizeof(void*);

		// Add node size
		size_t nodeSize = sizeof(std::pair<const K, V>) + sizeof(void*); // Node contains pair and pointer
		memoryUsage += map.size() * nodeSize;

		// For string keys, consider the memory of the strings themselves
		if constexpr(std::is_same_v<K, std::string>) {
			for (const auto& pair : map) {
				// String capacity (usually larger than size) plus pointer and size_t (except for SSO optimization)
				memoryUsage += pair.first.capacity() + 1 + sizeof(size_t) + sizeof(void*);
			}
		}

		// For string values, also consider string memory
		if constexpr(std::is_same_v<V, std::string>) {
			for (const auto& pair : map) {
				memoryUsage += pair.second.capacity() + 1 + sizeof(size_t) + sizeof(void*);
			}
		}

		return memoryUsage;
	}

	// Method 2: More detailed memory analysis
	template <typename K, typename V, typename Hash = std::hash<K>,
				typename Pred = std::equal_to<K>,
				typename Alloc = std::allocator<std::pair<const K, V>>>
	static void analyzeMapMemory(const std::unordered_map<K, V, Hash, Pred, Alloc>& map) {
		// Get basic information
		size_t size = map.size();
		size_t bucketCount = map.bucket_count();
		float loadFactor = size / static_cast<float>(bucketCount > 0 ? bucketCount : 1);

		// Calculate estimated memory
		size_t estimatedMemory = estimateMemoryUsage(map);

		// Output analysis results
		std::cout << "===== unordered_map Memory Analysis =====\n";
		std::cout << "Type: unordered_map<"
					<< typeid(K).name() << ", "
					<< typeid(V).name() << ">\n";
		std::cout << "Element count: " << size << "\n";
		std::cout << "Bucket count: " << bucketCount << "\n";
		std::cout << "Load factor: " << loadFactor << "\n";
		std::cout << "Max load factor: " << map.max_load_factor() << "\n";
		std::cout << "Estimated memory usage: " << estimatedMemory << " bytes ("
					<< (estimatedMemory / 1024.0) << " KB)\n";
		std::cout << "Average memory per element: "
					<< (size > 0 ? estimatedMemory / static_cast<double>(size) : 0)
					<< " bytes\n";

		// Analyze bucket distribution
		size_t emptyBuckets = 0;
		size_t maxBucketSize = 0;
		size_t totalChainLength = 0;

		for (size_t i = 0; i < bucketCount; ++i) {
			size_t bucketSize = map.bucket_size(i);
			if (bucketSize == 0) {
				emptyBuckets++;
			}
			if (bucketSize > maxBucketSize) {
				maxBucketSize = bucketSize;
			}
			totalChainLength += bucketSize;
		}

		std::cout << "Empty buckets: " << emptyBuckets << " ("
					<< (bucketCount > 0 ? emptyBuckets * 100.0 / bucketCount : 0) << "%)\n";
		std::cout << "Max chain length: " << maxBucketSize << "\n";
		std::cout << "Average chain length: " << (bucketCount > 0 ? totalChainLength / static_cast<double>(bucketCount) : 0) << "\n";
		std::cout << "=========================================\n";
	}

	// Method 3: Measure memory difference before and after copying (requires OS support)
	template <typename MapType>
	static size_t measureByDifference(const MapType& map) {
		size_t beforeMem = getCurrentRSS();

		// Create a copy of the map, which will allocate new memory
		MapType mapCopy = map;

		size_t afterMem = getCurrentRSS();

		// Return the difference, which should approximate the actual memory usage of the map
		return afterMem - beforeMem;
	}

private:
	// Get current process memory usage (Resident Set Size)
	static size_t getCurrentRSS() {
#if defined(_WIN32)
		// Windows implementation
		PROCESS_MEMORY_COUNTERS pmc;
		if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
			return pmc.WorkingSetSize;
		}
		return 0;
#elif defined(__unix__) || defined(__linux__) || defined(__APPLE__)
		// Unix/Linux/Mac implementation
		FILE* file = fopen("/proc/self/statm", "r");
		if (file) {
			long rss = 0;
			if (fscanf(file, "%*s %ld", &rss) == 1) {
				fclose(file);
				return rss * sysconf(_SC_PAGESIZE);
			}
			fclose(file);
		}

		// If the above method fails, try using getrusage
		struct rusage usage;
		if (getrusage(RUSAGE_SELF, &usage) == 0) {
			return usage.ru_maxrss * 1024;
		}
		return 0;
#else
		return 0; // Unsupported platform
#endif
	}
};
#elif(C_MAP_TYPE == 2)
class MapMemoryMeasurer {
	public:
		// Measure std::map memory usage
		template <typename K, typename V, typename Compare = std::less<K>,
				  typename Alloc = std::allocator<std::pair<const K, V>>>
		static size_t estimateMemoryUsage(const std::map<K, V, Compare, Alloc>& map) {
			// Basic container size
			size_t memoryUsage = sizeof(std::map<K, V, Compare, Alloc>);

			// std::map typically uses red-black tree implementation, each node contains pair, color flag and three pointers (parent, left, right)
			size_t nodeOverhead = 3 * sizeof(void*) + sizeof(char); // Three pointers and color flag
			size_t nodeSize = sizeof(std::pair<const K, V>) + nodeOverhead;

			// Add memory for all nodes
			memoryUsage += map.size() * nodeSize;

			// Handle string type keys
			if constexpr(std::is_same_v<K, std::string>) {
				for (const auto& pair : map) {
					memoryUsage += estimateStringMemory(pair.first);
				}
			}

			// Handle string type values
			if constexpr(std::is_same_v<V, std::string>) {
				for (const auto& pair : map) {
					memoryUsage += estimateStringMemory(pair.second);
				}
			}

			// Handle nested containers
			if constexpr(is_container_v<V>) {
				for (const auto& pair : map) {
					memoryUsage += estimateContainerMemory(pair.second);
				}
			}

			return memoryUsage;
		}
		template <typename K, typename V, typename Compare = std::less<K>,
              typename Alloc = std::allocator<std::pair<const K, V>>>
    static void analyzeMapMemory(const std::map<K, V, Compare, Alloc>& map) {
        // Get basic information
        size_t size = map.size();

        // Calculate estimated memory
        size_t estimatedMemory = estimateMemoryUsage(map);

        // Output analysis results
        std::cout << "===== std::map Memory Analysis =====\n";
        std::cout << "Type: map<"
                  << typeid(K).name() << ", "
                  << typeid(V).name() << ">\n";
        std::cout << "Element count: " << size << "\n";
        std::cout << "Estimated memory usage: " << estimatedMemory << " bytes ("
                  << (estimatedMemory / 1024.0) << " KB)\n";
        std::cout << "Average memory per element: "
                  << (size > 0 ? estimatedMemory / static_cast<double>(size) : 0)
                  << " bytes\n";

        // Red-black tree height estimation (log2(n))
        if (size > 0) {
            double estimatedHeight = std::log2(size);
            std::cout << "Estimated tree height: ~" << estimatedHeight << "\n";
        }

        std::cout << "=====================================\n";
    }
	// Measure memory difference before and after copying
    template <typename MapType>
    static size_t measureByDifference(const MapType& map) {
        // Force GC
        std::malloc(1);
        std::free(std::malloc(1));

        size_t beforeMem = getCurrentRSS();

        // Create a copy of the map, which will allocate new memory
        MapType mapCopy = map;

        // Ensure all memory is actually allocated
        volatile char dummy = 0;
        for (const auto& pair : mapCopy) {
            dummy += reinterpret_cast<const char*>(&pair)[0];
        }

        // Force GC
        std::malloc(1);
        std::free(std::malloc(1));

        size_t afterMem = getCurrentRSS();

        // Return the difference, which should approximate the actual memory usage of the map
        return afterMem - beforeMem;
    }
	private:
    // Trait to detect if a type is a container
    template <typename T, typename = void>
    struct is_container : std::false_type {};

    template <typename T>
    struct is_container<T,
        std::void_t<
            typename T::value_type,
            typename T::size_type,
            typename T::iterator,
            decltype(std::declval<T>().size()),
            decltype(std::declval<T>().begin()),
            decltype(std::declval<T>().end())
        >> : std::true_type {};

    template <typename T>
    static constexpr bool is_container_v = is_container<T>::value;

    // Estimate std::string memory usage
    static size_t estimateStringMemory(const std::string& str) {
        // Check if SSO (Small String Optimization) is used
        // This threshold varies by implementation, typically between 15-23 bytes
        constexpr size_t ssoThreshold = 15; // Assume SSO threshold is 15

        if (str.size() <= ssoThreshold) {
            return 0; // String is contained within the string object itself
        } else {
            // Non-SSO strings require additional heap allocation
            // Calculate actual allocated capacity (usually power of 2 or close value)
            return str.capacity() + 1; // +1 for null terminator
        }
    }

    // Estimate container memory usage
    template <typename Container>
    static size_t estimateContainerMemory(const Container& container) {
        using ValueType = typename Container::value_type;

        size_t memoryUsage = 0;

        // For vector and similar containers, consider capacity
        if constexpr(has_capacity<Container>::value) {
            memoryUsage += container.capacity() * sizeof(ValueType);
        }

		// For list and similar containers, consider node overhead
		else if constexpr(std::is_same_v<Container, std::list<ValueType>> ||
			std::is_same_v<Container, std::forward_list<ValueType>>) {
		  // Each node contains value and at least one pointer
		  memoryUsage += container.size() * (sizeof(ValueType) + sizeof(void*));
	  }
	  // For set/map and other tree-based containers
	  else if constexpr(std::is_same_v<Container, std::set<ValueType>> ||
						std::is_same_v<Container, std::multiset<ValueType>> ||
						std::is_same_v<Container, std::map<typename Container::key_type, typename Container::mapped_type>>) {
		  // Red-black tree nodes typically contain value, color flag and three pointers
		  memoryUsage += container.size() * (sizeof(ValueType) + sizeof(char) + 3 * sizeof(void*));
	  }
	  // Default case
	  else {
		  memoryUsage += container.size() * sizeof(ValueType);
	  }

	  // Recursively handle strings or nested containers within the container
	  if constexpr(is_string_or_container_v<ValueType>) {
		  for (const auto& item : container) {
			  memoryUsage += estimateNestedMemory(item);
		  }
	  }

	  return memoryUsage;
  }

  // Handle nested strings or containers
  template <typename T>
  static size_t estimateNestedMemory(const T& value) {
	  if constexpr(std::is_same_v<T, std::string>) {
		  return estimateStringMemory(value);
	  }
	  else if constexpr(is_container_v<T>) {
		  return estimateContainerMemory(value);
	  }
	  else {
		  return 0;
	  }
  }

  // Detect if type is string or container
  template <typename T>
  static constexpr bool is_string_or_container_v =
	  std::is_same_v<T, std::string> || is_container_v<T>;

  // Detect if container has capacity method
  template <typename T, typename = void>
  struct has_capacity : std::false_type {};

  template <typename T>
  struct has_capacity<T, std::void_t<decltype(std::declval<T>().capacity())>>
	  : std::true_type {};

  // Get current process memory usage (Resident Set Size)
  static size_t getCurrentRSS() {
#ifdef _WIN32
	  // Windows implementation
	  PROCESS_MEMORY_COUNTERS pmc;
	  if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
		  return pmc.WorkingSetSize;
	  }
	  return 0;
#elif defined(__unix__) || defined(__linux__) || defined(__APPLE__)
	  // Unix/Linux/Mac implementation
	  FILE* file = fopen("/proc/self/statm", "r");
	  if (file) {
		  long rss = 0;
		  if (fscanf(file, "%*s %ld", &rss) == 1) {
			  fclose(file);
			  return rss * sysconf(_SC_PAGESIZE);
		  }
		  fclose(file);
	  }

	  // If the above method fails, try using getrusage
	  struct rusage usage;
	  if (getrusage(RUSAGE_SELF, &usage) == 0) {
		  return usage.ru_maxrss * 1024;
	  }
	  return 0;
#else
	  return 0; // Unsupported platform
#endif
  }

};
#endif

void Obj_Manager::initial()
{

	// Initialize solid geometry
	std::vector<Ini_Shape> ini_shape(1); // number of solids

	// ini_shape.at(0).shape_type = geofile;
	ini_shape.at(0).shape_type = geofile_stl;
	// ini_shape.at(0).bool_moving = true;
	ini_shape.at(0).bool_moving = false;
	// Set solid position: STL local origin (0,0,0) will be placed at
	// C_solid_origin (defined per-namespace in Constants.h). Edit there to
	// manually control where the solid sits inside the domain.
	ini_shape.at(0).length.push_back(C_solid_origin[0]);  // x origin
	ini_shape.at(0).length.push_back(C_solid_origin[1]);  // y origin
	ini_shape.at(0).length.push_back(C_solid_origin[2]);  // z origin

	// ini_shape.at(0).shape_type = circle;
	// ini_shape.at(0).bool_moving = true;
	// ini_shape.at(0).numb_nodes = 2000;
	// ini_shape.at(0).length.push_back(C_xb / 2);
	// ini_shape.at(0).length.push_back(C_yb / 2);
	// ini_shape.at(0).length.push_back(1);

	// ini_shape.at(1).shape_type = geofile_stl;
	// // ini_shape.at(0).bool_moving = true;
	// ini_shape.at(1).bool_moving = false;
	// ini_shape.at(1).length.push_back(C_xb*2 / 3);
    // ini_shape.at(1).length.push_back(C_yb*2 / 3);
	// ini_shape.at(1).length.push_back(C_zb*2 / 3);

	//ini_shape.at(1).shape_type = line_fillx;
	//ini_shape.at(1).x0 = 20.;
	//ini_shape.at(1).y0 = 20.;
	//ini_shape.at(1).numb_nodes = 200;
	//ini_shape.at(1).length.push_back(-1.);
	//ini_shape.at(1).length.push_back(20.);

	Solid_Manager::pointer_me = &solid_manager;
	solid_manager.initial(ini_shape);
	// solid_manager.output_pointCloud(ini_shape);

	std::vector <Ini_Shape>().swap(ini_shape);

	// Generate and initialize mesh
	Grid_Manager::pointer_me = &gr_manager;
	Timer tmmr;
	double st = tmmr.elapsed();
	gr_manager.initial();
	double ed = tmmr.elapsed();

	std::cout << "Grid Size: " << gr_manager.gr_inner.grid.size() << std::endl;
	std::cout << "Grid space: " << gr_manager.gr_inner.get_dx() << std::endl;
	std::cout << "Generation time: " << ed - st << " s." << std::endl;

#if (C_MAP_TYPE == 1) || (C_MAP_TYPE == 2)
	MapMemoryMeasurer::analyzeMapMemory(gr_manager.gr_inner.grid);

	size_t diffMemory = MapMemoryMeasurer::measureByDifference(gr_manager.gr_inner.grid);
	std::cout << "Memory measured by difference method: " << diffMemory << " bytes ("
	<< (diffMemory / 1024.0) << " KB)\n\n";
#endif


	auto now = std::chrono::system_clock::now();
	auto now_time_t = std::chrono::system_clock::to_time_t(now);
	std::tm local_tm{};
    #ifdef _WIN32
    localtime_s(&local_tm, &now_time_t);
    #else
    localtime_r(&now_time_t, &local_tm);  // POSIX standard
    #endif
	std::stringstream ss;
    ss << OUTPUT_NAME << "_dx" << gr_manager.gr_inner.get_dx() << "_"
       << std::put_time(&local_tm, "%Y%m%d_%H%M%S") << ".txt";



	std::ofstream amr_nodes_file;
	amr_nodes_file.open(outputPath(ss.str()), std::ios::out);
	for (auto iter =  gr_manager.gr_inner.grid.begin(); iter != gr_manager.gr_inner.grid.end(); ++iter) {
		amr_nodes_file << iter->first << "\n";
	}
	amr_nodes_file.close();
	

	// Generate lattices from mesh (node = vertex)
	Lat_Manager::pointer_me = &lat_manager;
	lat_manager.voxelize();
	// lat_manager.initial(); // old version
}

void Obj_Manager::time_marching_management()
{
	D_real sum_t = 0.;
	std::array<D_mapint, C_max_level + 1> map_add_nodes, map_remove_nodes;
	while (sum_t < 3 * C_dx)
	{
		time_marching(sum_t, map_add_nodes, map_remove_nodes);
		sum_t += C_dx;
	}

	
}

void Obj_Manager::time_marching(D_real sum_t, std::array<D_mapint, C_max_level + 1>  &map_add_nodes, std::array<D_mapint, C_max_level + 1>  &map_remove_nodes)
{
	D_real dt = C_dx / static_cast<D_real> (two_power_n(C_max_level));
	std::array<unsigned int, C_max_level + 1> accumulate_t{};
	std::array<unsigned int, C_max_level + 1> flag_time_step{}; // record number of time step at ilevel

	for (std::vector<unsigned int>::iterator iter = run_order.begin(); iter != run_order.end(); ++iter)
	{
		Timer tmr;
		double t0 = tmr.elapsed();

		unsigned int ilevel = *iter;
		++flag_time_step[ilevel];

		D_real current_t;

		accumulate_t[ilevel] += two_power_n(C_max_level - ilevel);
		current_t = sum_t + dt * static_cast<D_real> (accumulate_t[ilevel]);

		if (ilevel == C_max_level)
		{
#if (C_SOLID_BOUNDARY == 2)
			// update shape information
			unsigned int numb_solids = Solid_Manager::pointer_me->numb_solids;
			for (unsigned int ishape = 0; ishape < numb_solids; ++ishape)
			{
				if (Solid_Manager::pointer_me->shape_solids.at(ishape).bool_moving)
				{
					Solid_Manager::pointer_me->renew(ishape, current_t);				
					Grid_Manager::pointer_me->update_nodes_near_solid(ishape, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel));
					Grid_Manager::pointer_me->update_map_node_IB(ilevel, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel));

					if (Solid_Manager::pointer_me->shape_solids.at(ishape).bool_enclosed)
					{
						Grid_Manager::pointer_me->update_ghost_node(ilevel, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel));
					}
				}

			}
			
			if ((flag_time_step[ilevel] % 2) == 0)
			{
				Grid_Manager::pointer_me->call_update_nodes(ilevel, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel), map_add_nodes.at(ilevel - 1), map_remove_nodes.at(ilevel - 1));
				map_add_nodes[ilevel].clear();
				map_remove_nodes[ilevel].clear();
			}

#endif			

		}
		else if (ilevel == 0)
		{

		}
		else
		{
#if (C_SOLID_BOUNDARY == 2)
			if ((flag_time_step[ilevel] % 2) == 0)
			{
				Grid_Manager::pointer_me->call_update_nodes(ilevel, map_add_nodes[ilevel], map_remove_nodes[ilevel], map_add_nodes[ilevel - 1], map_remove_nodes[ilevel - 1]);
				map_add_nodes[ilevel].clear();
				map_remove_nodes[ilevel].clear();
			}
#endif	
		}

		//double t1 = tmr.elapsed();
		//double t2 = tmr.elapsed();
		//std::cout << ilevel<< ": " << t1 - t0 << ", " << t2 - t1 << std::endl;
	}
	
}

void Obj_Manager::output()
{
	// time
	auto now = std::chrono::system_clock::now();
	auto now_time_t = std::chrono::system_clock::to_time_t(now);
	std::tm local_tm{};
    #ifdef _WIN32
    localtime_s(&local_tm, &now_time_t);
    #else
    localtime_r(&now_time_t, &local_tm);  // POSIX standard
    #endif
	std::stringstream ss;
    ss << OUTPUT_NAME << "_"
       << std::put_time(&local_tm, "%Y%m%d_%H%M%S");

	// Write mesh using Tecplot format
	IO_Manager::pointer_me = &io_manager;
	io_manager.outfile = ss.str();

	std::vector<unsigned int> out_vlevel = {C_max_level};
	// std::vector<unsigned int> out_vlevel = {3};
	// std::vector<unsigned int> out_vlevel;
	// for (unsigned ilevel = 0; ilevel < C_max_level + 1; ++ilevel)
	// {
	// 	out_vlevel.push_back(ilevel);
	// }
	io_manager.vlevel = out_vlevel;

	io_manager.writeMesh();  // Write mesh in Tecplot format
}