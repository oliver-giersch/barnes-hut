#include <stdio.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "barnes-hut/common.h"
#include "barnes-hut/options.h"
#include "barnes-hut/phys.h"

static const unsigned width	 = 1024;
static const unsigned height = 768;

static SDL_Window *window	  = NULL;
static SDL_GLContext *context = NULL;

static void render_axes(float radius);

int
render_init(void)
{
	SDL_SetHint(SDL_HINT_NO_SIGNAL_HANDLERS, "1");
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

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-500.0, 500.0, -500.0, 500.0, 0.1, 1000.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glEnable(GL_DEPTH_TEST);
	glPointSize(1.0);

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

void
render_scene(const struct accel_particle particles[], float radius)
{
	// TODO: check for close event
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	gluLookAt(1.25 * radius, 1.25 * radius, 1.25 * radius, 0.0, 0.0, 0.0, 0.0,
		1.0, 0.0);

	render_axes(radius);

	glBegin(GL_POINTS);
	for (size_t p = 0; p < options.particles; p++) {
		const struct vec3 *pos = &particles[p].part.pos;
		glVertex3f(pos->x, pos->y, pos->z);
	}

	glEnd();
	SDL_GL_SwapWindow(window);

	// preparation: mmapped SHM file, byte buffer:
	// glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE,
	// pixelData);
}

static void
render_axes(float radius)
{
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