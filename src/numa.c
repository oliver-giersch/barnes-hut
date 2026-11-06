#include "barnes-hut/numa.h"

#include <errno.h>
#include <stdlib.h>

#include <hwloc.h>

#include "barnes-hut/types.h"

thread_local unsigned numa_node_id = 0;

struct numa_cpu_info {
	hwloc_obj_t	pu;
	unsigned	node_id;
};

static hwloc_topology_t topology;
static int cpu_count;
static struct numa_cpu_info *cpu_info;

extern int
numa_init(void)
{
	error_t error;

	if ((error = hwloc_topology_init(&topology)))
		return errno;
	if ((error = hwloc_topology_load(topology)))
		return errno;

	if ((cpu_count = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU)) == -1)
		return errno;

	if (cpu_count == 0)
		return ENODEV;

	cpu_info = malloc(cpu_count * sizeof(struct numa_cpu_info));
	if (cpu_info == NULL)
		return ENOMEM;

	for (int i = 0; i < cpu_count; i++) {
		const hwloc_obj_t pu = hwloc_get_obj_by_type(topology, HWLOC_OBJ_PU, i);
		const unsigned os_index = hwloc_bitmap_first(pu->nodeset);
		const hwloc_obj_t nn = hwloc_get_numanode_obj_by_os_index(
			topology, os_index);

		if (nn->logical_index >= numa_max_nodes)
			return ENOTSUP;

		cpu_info[i] = (struct numa_cpu_info){
			.pu = pu,
			.node_id = nn->logical_index,
		};
	}

	return (0);
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
