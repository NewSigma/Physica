/*
 * Copyright 2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <hwloc.h>
#include <hwloc/cpukinds.h>
#include <cassert>
#include <climits>
#include <limits>
#include <print>

namespace {
    hwloc_bitmap_t getPerformanceCore(hwloc_topology_t topo) {
        const int kinds = hwloc_cpukinds_get_nr(topo, 0);
        assert(kinds >= 0);
        const auto* support = hwloc_topology_get_support(topo);
        const int effKnown = support != nullptr && support->discovery != nullptr && support->discovery->cpukind_efficiency;

        hwloc_bitmap_t result = hwloc_bitmap_alloc();
        assert(result != nullptr && "hwloc_bitmap_alloc");
        bool isHybrid = kinds > 0;
        if (isHybrid && effKnown) {
            int max = std::numeric_limits<int>::min();
            int cpu = 0;
            for (int k = 0; k < kinds; ++k) {
                int eff = 0;
                hwloc_cpukinds_get_info(topo, k, nullptr, &eff, nullptr, nullptr, 0);
                if (eff > max) {
                    max = eff;
                    cpu = k;
                }
            }
            hwloc_cpukinds_get_info(topo, cpu, result, nullptr, nullptr, nullptr, 0);
        }
        else {
            hwloc_obj_t pu = hwloc_get_obj_by_type(topo, HWLOC_OBJ_PU, 0);
            assert(pu != nullptr && "no PU in topology");
            hwloc_bitmap_set(result, pu->os_index);
        }
        return result;
    }

    hwloc_uint64_t getCacheSize(hwloc_topology_t topo, hwloc_obj_type_t type, int cpu) {
        const int depth = hwloc_get_type_depth(topo, type);
        assert(depth >= 0 && "cache type not found in topology");
        const int nb = hwloc_get_nbobjs_by_depth(topo, depth);
        hwloc_uint64_t found = 0;
        for (int i = 0; i < nb; ++i) {
            hwloc_obj_t obj = hwloc_get_obj_by_depth(topo, depth, i);
            if (obj == nullptr || obj->cpuset == nullptr)
                continue;
            if (!hwloc_bitmap_isset(obj->cpuset, cpu))
                continue;
            assert((found == 0 || found == obj->attr->cache.size) && "inconsistent cache size across objects covering the selected cpu");
            found = obj->attr->cache.size;
        }
        assert(found != 0 && "no cache object covers the selected cpu");
        return found;
    }

    unsigned getLineSize(hwloc_topology_t topo, int cpu) {
        const int depth = hwloc_get_type_depth(topo, HWLOC_OBJ_L1CACHE);
        assert(depth >= 0 && "L1D cache type not found in topology");
        const int nb = hwloc_get_nbobjs_by_depth(topo, depth);
        unsigned line = 0;
        for (int i = 0; i < nb; ++i) {
            hwloc_obj_t obj = hwloc_get_obj_by_depth(topo, depth, i);
            if (obj != nullptr && obj->cpuset != nullptr && hwloc_bitmap_isset(obj->cpuset, cpu)) {
                line = obj->attr->cache.linesize;
                break;
            }
        }
        assert(line != 0 && "cache line size unknown");
        return line;
    }
}

int main() {
    hwloc_topology_t topo{};
    hwloc_topology_init(&topo);
    hwloc_topology_load(topo);

    auto* coreset = getPerformanceCore(topo);
    const int cpu = hwloc_bitmap_first(coreset);
    assert(cpu >= 0 && "selected cpu kind has no cpu");

    const hwloc_uint64_t l1d = getCacheSize(topo, HWLOC_OBJ_L1CACHE, cpu);
    const hwloc_uint64_t l2 = getCacheSize(topo, HWLOC_OBJ_L2CACHE, cpu);
    const hwloc_uint64_t l3 = getCacheSize(topo, HWLOC_OBJ_L3CACHE, cpu);
    const unsigned line = getLineSize(topo, cpu);
    hwloc_bitmap_free(coreset);
    hwloc_topology_destroy(topo);
    std::println("CacheProbe L1D={} L2={} L3={} line={}", l1d, l2, l3, line);
    return 0;
}
