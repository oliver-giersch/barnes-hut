src = [
  'options.c',
  'physics.c',
  'thread.c'
]

hwloc = dependency('hwloc', required: numa_aware)
if numa_aware.enabled()
  src += 'numa.c'
endif

os = build_machine.system()
c_args = []
link_args = []

executable('barnes_hut',
  src,
  include_directories: include_dir,
  c_args: c_args,
  link_args: link_args,
  dependencies: [hwloc]
)
