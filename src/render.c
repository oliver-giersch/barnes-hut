#include "barnes-hut/render.h"

static void render_event_loop(SDL_Window *window,
	const struct particle particles[static particle_count]);

static void
render_event_loop(SDL_Window *window,
	const struct particle particles[static particle_count])
{
	SDL_Event event;
	while (true) {
		if (unlikely(!SDL_WaitEvent(&event)))
			return; // error
		else if (unlikely(event.type == SDL_EVENT_QUIT)
			return;
		else if (event.type != SDL_EVENT_USER)
			continue;

		if ((error = render_draw_scene(window)
			return; // error
		// try to SDL_GL_SwapWindow(window)
		// store render state idle
		// wake notify
	}
}
