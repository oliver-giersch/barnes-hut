project('barnes-hut', 'c',
    version: '0.0.0',
    meson_version: '>=1.1',
    default_options: ['c_std=c23', 'warning_level=3']
)

add_project_arguments('-D_DEFAULT_SOURCE', language: 'c')

if get_option('double_precision')
  add_project_arguments('-DOPTION_DOUBLE_PRECISION=1', language: 'c')
endif

if get_option('large_simulation')
  add_project_arguments('-DOPTION_LARGE_SIMULATION=1', language: 'c')
endif

numa_aware = get_option('numa_aware')
if numa_aware.enabled()
  add_project_arguments('-DOPTION_NUMA_AWARE=1', language: 'c')
endif

if get_option('print_progress')
  add_project_arguments('-DOPTION_PRINT_PROGRESS=1', language: 'c')
endif

render = get_option('render')
if render.enabled()
  add_project_arguments('-DOPTION_RENDER=1', language: 'c')
endif

include_dir = include_directories('include')
subdir('src')
