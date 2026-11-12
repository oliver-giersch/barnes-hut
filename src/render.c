#include "barnes-hut/render.h"

#include <glew/glew.h>

#include "barnes-hut/types.h"

typedef long gl_error_t;

static constexpr int width = 1280;
static constexpr int height = 960;

static constexpr char fragment_shader[] = {
	#embed "shader/fragment.glsl"
	, '\0'
};

static constexpr char vertex_shader[] = {
	#embed "shader/vertex.glsl"
	, '\0'
};

static struct gl_state {
	GLuint program;
	/* Program arguments */
	GLuint view, radius, axis;
	GLuint particle_vbo, particle_vao;
	GLuint axes_vbo, axes_vao
} gl_state;

static error_t render_init_state(void);
static error_t render_compile_shader(const char name[]);
static void render_event_loop(SDL_Window *window,
	const struct particle particles[static particle_count]);

extern error_t
render_main(void)
{
	
}

static error_t
render_init_state(void)
{
	if (!SDL_Init(SDL_INIT_VIDEO)
		return -1;

	SDL_Window *window = SDL_CreateWindow("barnes-hut", width, height,
		SDL_WINDOW_OPENGL);
	if (window == nullptr)
		return -1;

	SDL_GL_Context *context = SDL_GL_CreateContext(window);
	if (context == nullptr)
		goto destroy_window;

	if (glewContextInit != GLEW_OK)
		goto destroy_context;

destroy_context:
	SDL_GL_DestroyContext(context);
destroy_window:
	SDL_DestroyWindow(window);
destroy_sdl:
	SDL_Quit();
}

static gl_error_t
render_gl_state_init(void)
{
	gl_error_t error;

	error = render_link_program();
	if (error <= 0)
		return error;

	const GLuint program = (GLuint)error;

	GLint view = glGetUniformLocation(program, "view");
	if (view == -1)
		return -1;

	GLint radius = glGetUniformLocation(program, "radius");
	if (view == -1)
		return -1;

	GLint axis = glGetUniformLocation(program, "axis");
	if (axis == -1)
		return -1;

	GLuint particle_vao, particle_vbo;
	glGenVertexArrays(1, &particle_vao);
	glBindVertexArray(particle_vao);
	glGenBuffers(1, &particle_vbo);

	GLuint axes_vao, axes_vbo;
	glBindVertexArray(&axes_vao);
	glBindVertexArray(axes_vao);
	glGenBuffers(1, &axes_vbo);

	glPointSize(0.5);
	glEnable(GL_DEPTH_TEST);

	gl_state = (struct gl_state){
		.program = program,
		.view = view,
		.radius = radius,
		.axis = axis,
		.particle_vao = particle_vao,
		.particle_vbo = particle_vbo,
		.axes_vao = axes_vao,
		.axes_vbo = axes_vbo,
	};

	return 0;
}

static void
render_gl_state_destroy(void)
{
	glDeleteBuffers(gl_state.axes_vbo);
	glDeleteVertexArrays(gl_state.axes_vao);
	glDeleteBuffers(gl_state.particle_vbo);
	glDeleteVertexArrays(gl_state.particle_vao);
	glDeleteProgram(gl_state.program);
}

static gl_error_t
render_link_program(void)
{
	gl_error_t error;

	error = render_compile_shader(GL_FRAGMENT_SHADER, fragment_shader);
	if (unlikely(error < 0))
		return error;
	const GLuint fragment_shader = (GLuint)error;

	error = render_compile_shader(GL_VERTEX_SHADER, vertex_shader);
	if (unlikely(error < 0))
		return error;
	const GLuint vertex_shader = (GLuint)error;

	const GLuint program = glCreateProgram();
	if (program == 0)
		return -1;

	glAttachShader(program, fragment_shader);
	glAttachShader(program, vertex_shader);
	glLinkProgram(program);

	GLint ok;
	glGetProgramiv(program, GL_LINK_STATUS, &ok);
	if (ok == GL_FALSE) {
	
	}

	return (gl_error_t)program;
}

static gl_error_t
render_compile_shader(GLenum shader_type, const char source[])
{
	GLuint shader = glCreateShader(shader_type);
	if (shader == 0)
		return -1;

	glShaderSource(shader, 1, source, nullptr); //fixme: errno?
	glCompileShader(shader);

	GLint ok = 0;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
	if (ok == GL_FALSE) {

		goto delete_shader;
	}

delete_shader:
	
}

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
