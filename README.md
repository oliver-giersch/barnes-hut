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
$ sudo apt-get install libgl1-mesa-dev libglu1-mesa-dev libsdl2-dev
```

### Compile

For default settings *with* 3D rendering:

```console
$ make RENDER=1
```

For a debug build:

```console
$ make BUILD=debug
```
