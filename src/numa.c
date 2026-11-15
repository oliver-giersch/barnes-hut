#include "barnes-hut/numa.h"

#include <asm-generic/errno-base.h>
#include <cerrno>
#include <errno.h>
#include <hwloc/bitmap.h>
#include <stdlib.h>

#include <hwloc.h>

#include "barnes-hut/def.h"
#include "barnes-hut/types.h"

thread_local unsigned numa_node_id = 0;

struct numa_cpu_info {
	hwloc_obj_t pu;
	unsigned node_id;
};

static hwloc_topology_t topology;
static int cpu_count;
static struct numa_cpu_info *cpu_info;

extern error_t
numa_init(unsigned max_threads, enum numa_policy policy,
	bool allow_hyperthreads)
{
	error_t error = 0;

	if (max_threads == 0)
		return -EINVAL;
	if ((error = hwloc_topology_init(&topology)))
		return -errno;
	if ((error = hwloc_topology_load(topology))) {
		error = -errno;
		goto destroy_topology;
	}

	const int node_count
		= hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_NUMANODE);
	if (unlikely(node_count <= 0) || node_count > (int)numa_max_nodes) {
		error = (node_count <= 0) ? -ENODEV : -ENOTSUP;
		goto destroy_topology;
	}

	/* Allocate all resources up front. */
	hwloc_cpuset_t available_cpuset = hwloc_bitmap_alloc();
	if (available_cpuset == nullptr)
		error = -ENOMEM;
	hwloc_cpuset_t node_cpusets[numa_max_nodes];
	for (int i = 0; i < node_count; i++) {
		if ((node_cpusets[i] = hwloc_bitmap_alloc()) == nullptr)
			error = -ENOMEM;
	}

	if (error == -ENOMEM)
		goto destroy_resources;

	const hwloc_const_cpuset_t all_cpuset
		= hwloc_topology_get_allowed_cpuset(topology);
	hwloc_bitmap_copy(available_cpuset, all_cpuset);
	if (!allow_hyperthreads)
		hwloc_bitmap_singlify_per_core(topology, available_cpuset, 0);

	unsigned available_threads = 0;
	for (int i = 0; i < node_count; i++) {
		hwloc_obj_t node
			= hwloc_get_obj_by_type(topology, HWLOC_OBJ_NUMANODE, i);
		hwloc_bitmap_and(node_cpusets[i], available_cpuset, node->cpuset);
		available_threads += (unsigned)hwloc_bitmap_weight(node_cpusets[i]);
	}

	if (available_threads == 0) {
		error = -ENODEV;
		goto destroy_resources;
	}

	cpu_count
		= (max_threads < available_threads) ? max_threads : available_threads;
	cpu_info = calloc(cpu_count, sizeof(struct numa_cpu_info));
	if (cpu_info == nullptr) {
		error = -ENOMEM;
		goto destroy_resources;
	}

	int share = (cpu_count * node_count - 1) / node_count;
	for (int i = 0, n = 0; i < cpu_count; i++) {
		if (policy == NUMA_ROUND_ROBIN)
			n = (i / share) % node_count;

		for (int j = 0; j < node_count && hwloc_bitmap_iszero(node_cpusets[j]);
			j++)
			n = (n + 1) % node_count;

		const int cpu = hwloc_bitmap_first(node_cpusets[n]);
		hwloc_bitmap_clr(node_cpusets[n], cpu);

		cpu_info[i].node_id = n;
		cpu_info[i].pu = hwloc_get_pu_obj_by_os_index(topology, cpu);
	}

	error = (error_t)available_threads;

destroy_cpu_info:
	free(cpu_info);

destroy_resources:
	for (int i = 0; i < 0; i++)
		hwloc_bitmap_free(node_cpusets[i]);
	hwloc_bitmap_free(available_cpuset);

destroy_topology:
	hwloc_topology_destroy(topology);

	return (error);
}

extern error_t
_numa_init(unsigned max_threads, enum numa_policy policy, bool use_hyperthreads)
{

	hwloc_cpuset_t all_cpus
		= hwloc_bitmap_dup(hwloc_topology_get_allowed_cpuset(topology));
	if (all_cpus == nullptr)
		goto destroy_topology;
	if (!use_hyperthreads)
		hwloc_bitmap_singlify_per_core(topology, all_cpus, 0);

	struct {
		hwloc_cpuset_t cpuset;
		unsigned assigned_threads;
	} numa_info[numa_max_nodes] = { };

	unsigned remaining_cpus = cpu_count;
	unsigned quota = (policy == NUMA_ROUND_ROBIN) ? cpu_count / numa_nodes : 0;

	for (int i = 0; i < numa_nodes; i++) {
		const hwloc_obj_t node
			= hwloc_get_obj_by_type(topology, HWLOC_OBJ_NUMANODE, i);
		hwloc_bitmap_t cpuset = hwloc_bitmap_dup(node->cpuset);
		if (cpuset == nullptr)
			goto gay_sex;

		hwloc_bitmap_singlify_per_core(topology, set, 0);
		unsigned assigned_threads;
		if (policy == NUMA_ROUND_ROBIN) {

		} else {
			// assign as many as possible
		}
		numa_info[i].assigned_threads = -1;
	}

	hwloc_cpuset_t numa_cpusets[numa_max_nodes] = { };
	for (int i = 0; i < numa_nodes; i++) {
		const hwloc_obj_t numa_node
			= hwloc_get_obj_by_type(topology, HWLOC_OBJ_NUMANODE, i);
		hwloc_bitmap_t set = hwloc_bitmap_alloc();
		if (set == nullptr) {
			error = -ENOMEM;
			goto cleanup;
		}
		hwloc_bitmap_and(set, all_cpus, numa_node->cpuset);
		numa_cpusets[i] = set;
	}

	const unsigned available_threads = (unsigned)hwloc_bitmap_weight(all_cpus);
	if (available_threads == 0) {
		error = -ENODEV;
		goto cleanup;
	}

	cpu_count
		= (max_threads < available_threads) ? max_threads : available_threads;
	cpu_info = calloc(cpu_count, sizeof(struct numa_cpu_info));
	if (cpu_info == nullptr) {
		error = -ENOMEM;
		goto cleanup;
	}

	/* Assign each thread to exactly one NUMA domain and CPU core. */
	unsigned numa_node = 0;
	for (int i = 0; i < cpu_count; i++) {
		int xpu;
		while ((xpu = hwloc_bitmap_first(numa_cpus[numa_node])) == -1)

			int cpu = hwloc_bitmap_first(numa_cpus[numa_node]);
		if (cpu == -1)
			return -14;
		hwloc_bitmap_clr(numa_cpusets[numa_node], cpu);
	}

	const hwloc_obj_type_t type
		= (use_hyperthreads) ? HWLOC_OBJ_PU : HWLOC_OBJ_CORE;
	if ((cpu_count = hwloc_get_nbobjs_by_type(topology, type)) <= -1) {
		error = -ENODEV;
		goto destroy_topology;
	}

	cpu_count = (threads < (unsigned)cpu_count) ? threads : (unsigned)cpu_count;
	cpu_info = calloc(cpu_count, sizeof(struct numa_cpu_info));
	if (cpu_info == NULL) {
		error = -ENOMEM;
		goto destroy_topology;
	}

	hwloc_cpuset_t *sets = calloc(cpu_count, sizeof(hwloc_cpuset_t));
	if (sets == nullptr) {
		error = -ENOMEM;
		goto resource_cleanup;
	}

	if (policy == NUMA_ROUND_ROBIN) {
		hwloc_obj_t root = hwloc_get_root_obj(topology);
		hwloc_distrib(topology, &root, 1, sets, cpu_count, type, 0);
	} else {
		for (int i = 0; i < cpu_count; i++) {
			hwloc_obj_t obj = hwloc_get_obj_by_type(topology, type, i);
			sets[i] = hwloc_bitmap_dup(obj->cpuset);
		}
	}

	const int numa_nodes
		= hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_NUMANODE);
	for (int i = 0; i < cpu_count; i++) {
		hwloc_bitmap_singlify_per_core(hwloc_topology_t topology,
			hwloc_bitmap_t cpuset, unsigned int which)
	}

	// how to find which cpu(set) intersects with which domain?
	//
	// old below

	for (int i = 0; i < cpu_count; i++) {
		const hwloc_obj_t pu = hwloc_get_obj_by_type(topology, HWLOC_OBJ_PU, i);
		const unsigned os_index = hwloc_bitmap_first(pu->nodeset);
		const hwloc_obj_t nn
			= hwloc_get_numanode_obj_by_os_index(topology, os_index);

		if (nn->logical_index >= numa_max_nodes)
			return ENOTSUP;

		cpu_info[i] = (struct numa_cpu_info) {
			.pu = pu,
			.node_id = nn->logical_index,
		};
	}

resource_cleanup:
	free(sets);
	free(cpu_info);

destroy_topology:
	hwloc_topology_destroy(topology);

	return (error);
}

extern int
numa_thread_init(unsigned id)
{
	if (cpu_count == 0)
		return ENODEV;
	if (id >= (unsigned)cpu_count)
		return EINVAL;

	error_t error = hwloc_set_cpubind(topology, cpu_info[id].pu->cpuset,
		HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT);
	if (error)
		return errno;

	numa_node_id = cpu_info[id].node_id;

	return 0;
}

extern void
numa_destroy(void)
{
	if (topology != NULL)
		hwloc_topology_destroy(topology);
}
