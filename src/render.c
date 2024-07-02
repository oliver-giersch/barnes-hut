#include <stdio.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "barnes-hut/common.h"
#include "barnes-hut/options.h"
#include "barnes-hut/phys.h"

static const unsigned width	 = 1280;
static const unsigned height = 960;

static SDL_Window *window	  = NULL;
static SDL_GLContext *context = NULL;

static void render_axes(float radius);
static void render_point(const struct vec3 *v, float radius);

int
render_init(void)
{
	SDL_SetHint(SDL_HINT_NO_SIGNAL_HANDLERS, "1");
	SDL_SetHint(SDL_HINT_VIDEODRIVER, "wayland");
	SDL_SetHint(SDL_HINT_VIDEO_WAYLAND_PREFER_LIBDECOR, "1");
	if (SDL_Init(SDL_INIT_VIDEO) < 0)
		goto error;

	window = SDL_CreateWindow("barnes-hut", SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED, width, height,
		SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
	if (window == NULL)
		goto deinit_sdl;

	context = SDL_GL_CreateContext(window);
	if (context == NULL)
		goto deinit_window;

	glViewport(0, 0, width, height);

	glEnable(GL_DEPTH_TEST);
	glPointSize(1.25);

	render_axes(options.radius);

	SDL_GL_SwapWindow(window);
	return 0;

deinit_window:
	SDL_DestroyWindow(window);
deinit_sdl:
	SDL_Quit();
error:
	return BHE_RENDER_ERROR;
}

void
render_deinit(void)
{
	SDL_GL_DeleteContext(context);
	SDL_DestroyWindow(window);
	SDL_Quit();
}

bool
render_scene(const struct accel_particle particles[], float radius)
{
	SDL_Event event;
	while (SDL_PollEvent(&event))
		if (event.type == SDL_QUIT)
			return true;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	render_axes(radius);

	glBegin(GL_POINTS);
	for (size_t p = 0; p < options.particles; p++)
		render_point(&particles[p].part.pos, radius);
	glEnd();

	SDL_GL_SwapWindow(window);

	return false;
}

static void
render_axes(float radius)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(-radius, radius, -radius, radius, 0.1, 3.0 * radius);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(radius, radius, radius, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

	glBegin(GL_LINES);

	// x-axis (red)
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(-radius, 0.0, 0.0);
	glVertex3f(radius, 0.0, 0.0);
	// y-axis (green)
	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, -radius, 0.0);
	glVertex3f(0.0, radius, 0.0);
	// z-axis (blue)
	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, -radius);
	glVertex3f(0.0, 0.0, radius);

	glEnd();
}

/*static inline float
sq(float x)
{
	return x * x;
}

static inline float
vec3_dist_sq(const struct vec3 *v, const struct vec3 *u)
{
	return sq(v->x - u->x) + sq(v->y - u->y) + sq(v->z - u->z);
}

static inline float
vec3_dist(const struct vec3 *v, const struct vec3 *u)
{
	return sqrtf(vec3_dist_sq(v, u));
}*/

// looks cool
/*static void
render_point(const struct vec3 *v, float radius)
{
	// static const float sqrt_twelve = 3.464101615;

	// const struct vec3 cam = { radius, radius, radius };
	// const float dist	  = vec3_dist(v, &cam);
	const struct vec3 zz = { 0.0, 0.0, 0.0 };
	const float dist	 = vec3_dist(v, &zz);
	const float blue	 = (dist / radius);
	// extern long int random(void);
	// const float red	  = (float)random() / (float)RAND_MAX;
	// const float green = (float)random() / (float)RAND_MAX;
	// blue			  = (float)random() / (float)RAND_MAX;
	// fprintf(stderr, "%f,%f,%f\n", red, green, blue);
	glColor3f(0.0, blue, blue);

	glVertex3f(v->x, v->y, v->z);
}*/

static void
render_point(const struct vec3 *v, float radius)
{
	const float d	  = 2 * radius;
	const float red	  = (v->x + radius) / d;
	const float green = (v->y + radius) / d;
	const float blue  = (v->z + radius) / d;

	glColor3f(red, green, blue);
	glVertex3f(v->x, v->y, v->z);
}