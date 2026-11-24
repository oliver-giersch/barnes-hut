src = ['main.c', 'options.c', 'physics.c', 'thread.c']

hwloc = dependency('hwloc', required: numa_aware)
if numa_aware.enabled()
  src += 'numa.c'
endif

sdl3 = dependency('SDL3', required: render)
gl = dependency('GL', required: render)
if render.enabled()
  src += 'render.c'
endif

os = build_machine.system()
c_args = []
link_args = []

executable(
  'barnes_hut',
  src,
  include_directories: include_dir,
  c_args: c_args,
  link_args: link_args,
  dependencies: [gl, hwloc, sdl3],
)
