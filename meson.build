project('barnes-hut', 'c',
    version: '0.0.0',
    meson_version: '>=1.1',
    default_options: ['c_std=c23', 'warning_level=3']
)

if get_option('double_precision')
  add_global_arguments('-DOPTION_DOUBLE_PRECISION=1', language: 'c')
endif

numa_aware = get_option('numa_aware')
if numa_aware.enabled()
  add_project_arguments('-DOPTION_NUMA_AWARE=1', language: 'c')
endif

if get_option('print_progress')
  add_global_arguments('-DPRINT_PROGRESS=1', language: 'c')
endif

include_dir = include_directories('include')
subdir('src')
