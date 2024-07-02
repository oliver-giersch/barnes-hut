# Barnes-Hut

## Build

### Dependencies

- GL
- GLU
- SDL2

#### Fedora/RHEL

```console
$ sudo dnf install mesa-libGL-devel mesa-libGLU-devel SDL2-devel
```

#### Debian/Ubuntu

```console
$ sudo apt-get install libgl1-mesa-dev libglu1-mesa libsdl2-dev
```

### Compile

```console
$ make BUILD=optimize-lto RENDER=1
```