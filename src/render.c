#include <GL/glu.h>

#include "barnes-hut/options.h"
#include "barnes-hut/phys.h"

static void
render_axes(float radius)
{
	glBegin(GL_LINES);

	// x-axis (red)
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(1.0f, 0.0f, 0.0f);
	// y-axis (green)
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 1.0f, 0.0f);
	// z-axis (blue)
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 1.0f);

	glEnd();
}

void
render_scene(struct accel_particle particles[], float radius)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	gluLookAt(radius, radius, -radius, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

	glBegin(GL_POINTS);
	for (size_t p = 0; p < options.particles; p++)
		glVertex3f(particles[p].x, particles[p].y, particles[p].z);
	glEnd();

	// preparation: mmapped SHM file, byte buffer:
	// glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE,
	// pixelData);
}