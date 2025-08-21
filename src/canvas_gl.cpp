#include <dymp/canvas_gl.h>
#include <dymp/util.h>
#include <dymp/model.h>

#include <GL/glut.h>

using namespace std;
using namespace Eigen;

namespace dymp{;
namespace render{;

//-------------------------------------------------------------------------------------------------
// CanvasGL

CanvasGL::CanvasGL(){

}
void CanvasGL::Line(const vec3_t& p0, const vec3_t& p1){
	glLineWidth(lineWidth);
	glBegin(GL_LINES);
	glColor4fv(&lineColor[0].rgba[0]);                                     glVertex3d(p0.x(), p0.y(), p0.z());
	glColor4fv(gradation ? &lineColor[1].rgba[0] : &lineColor[0].rgba[0]); glVertex3d(p1.x(), p1.y(), p1.z());
	glEnd();
}

void CanvasGL::BeginPath(){
	glLineWidth(lineWidth);
	glBegin(GL_LINE_STRIP);
	glColor4fv(&lineColor[0].rgba[0]);
}

void CanvasGL::EndPath(){
	glEnd();
}

void CanvasGL::MoveTo(const vec3_t& p){
	glVertex3d(p.x(), p.y(), p.z());
}

void CanvasGL::LineTo(const vec3_t& p){
	glVertex3d(p.x(), p.y(), p.z());
}

void CanvasGL::Point(const vec3_t& p){
	glPointSize(pointSize);
	glBegin(GL_POINTS);
	glColor4fv(&pointColor.rgba[0]);
	glVertex3d(p.x(), p.y(), p.z());
	glEnd();
}

inline void Rect(const vec2_t& rmin, const vec2_t& rmax){
	glVertex3d(rmin.x(), rmin.y(), 0.0);
	glVertex3d(rmax.x(), rmin.y(), 0.0);
	glVertex3d(rmax.x(), rmax.y(), 0.0);
	glVertex3d(rmin.x(), rmax.y(), 0.0);
}

void CanvasGL::Rectangle(const vec2_t& rmin, const vec2_t& rmax){
	if(drawFill){
		glBegin(GL_TRIANGLE_FAN);
		glColor4fv(&fillColor.rgba[0]);
		Rect(rmin, rmax);
		glEnd();
	}
	if(drawLine){
		glLineWidth(lineWidth);
		glBegin(GL_LINE_LOOP);
		glColor4fv(&lineColor[0].rgba[0]);
		Rect(rmin, rmax);
		glEnd();
	}
}

inline void Circ(const vec2_t& c, float r, const Rotation2D<real_t>& R, size_t ndiv){
	vec2_t p(r, 0.0);
	for(size_t i = 0; i <= ndiv; i++){
		glVertex3d(c.x() + p.x(), c.y() + p.y(), 0.0);
		p = R * p;
	}
}

void CanvasGL::Circle(const vec2_t& center, real_t r){
	Rotation2D<real_t> R(2.0*_pi/(real_t)nCircleDiv);
	if(drawFill){
		glBegin(GL_TRIANGLE_FAN);
		glColor4fv(&fillColor.rgba[0]);
		glVertex3f(center.x(), center.y(), 0.0);
		Circ(center, r, R, nCircleDiv);
		glEnd();
	}
	if(drawLine){
		glLineWidth(lineWidth);
		glBegin(GL_LINE_STRIP);
		glColor4fv(&lineColor[0].rgba[0]);
		Circ(center, r, R, nCircleDiv);
		glEnd();
	}
}

void CanvasGL::Box(const vec3_t& rmin, const vec3_t& rmax){
	if(drawFill){
		glColor4fv(&fillColor.rgba[0]);
		glBegin(GL_TRIANGLE_FAN);
		glVertex3d(rmin[0], rmin[1], rmin[2]);
		glVertex3d(rmin[0], rmax[1], rmin[2]);
		glVertex3d(rmin[0], rmax[1], rmax[2]);
		glVertex3d(rmin[0], rmin[1], rmax[2]);
		glEnd();

		glBegin(GL_TRIANGLE_FAN);
		glVertex3d(rmax[0], rmin[1], rmin[2]);
		glVertex3d(rmax[0], rmin[1], rmax[2]);
		glVertex3d(rmax[0], rmax[1], rmax[2]);
		glVertex3d(rmax[0], rmax[1], rmin[2]);
		glEnd();

		glBegin(GL_TRIANGLE_FAN);
		glVertex3d(rmin[0], rmin[1], rmin[2]);
		glVertex3d(rmin[0], rmin[1], rmax[2]);
		glVertex3d(rmax[0], rmin[1], rmax[2]);
		glVertex3d(rmax[0], rmin[1], rmin[2]);
		glEnd();

		glBegin(GL_TRIANGLE_FAN);
		glVertex3d(rmin[0], rmax[1], rmin[2]);
		glVertex3d(rmax[0], rmax[1], rmin[2]);
		glVertex3d(rmax[0], rmax[1], rmax[2]);
		glVertex3d(rmin[0], rmax[1], rmax[2]);
		glEnd();

		glBegin(GL_TRIANGLE_FAN);
		glVertex3d(rmin[0], rmin[1], rmin[2]);
		glVertex3d(rmax[0], rmin[1], rmin[2]);
		glVertex3d(rmax[0], rmax[1], rmin[2]);
		glVertex3d(rmin[0], rmax[1], rmin[2]);
		glEnd();

		glBegin(GL_TRIANGLE_FAN);
		glVertex3d(rmin[0], rmin[1], rmax[2]);
		glVertex3d(rmin[0], rmax[1], rmax[2]);
		glVertex3d(rmax[0], rmax[1], rmax[2]);
		glVertex3d(rmax[0], rmin[1], rmax[2]);
		glEnd();
	}
	if(drawLine){
		glLineWidth(lineWidth);
		glColor4fv(&lineColor[0].rgba[0]);

		glBegin(GL_LINE_LOOP);
		glVertex3d(rmin[0], rmin[1], rmin[2]);
		glVertex3d(rmin[0], rmax[1], rmin[2]);
		glVertex3d(rmin[0], rmax[1], rmax[2]);
		glVertex3d(rmin[0], rmin[1], rmax[2]);
		glEnd();

		glBegin(GL_LINE_LOOP);
		glVertex3d(rmax[0], rmin[1], rmin[2]);
		glVertex3d(rmax[0], rmax[1], rmin[2]);
		glVertex3d(rmax[0], rmax[1], rmax[2]);
		glVertex3d(rmax[0], rmin[1], rmax[2]);
		glEnd();

		glBegin(GL_LINES);
		glVertex3d(rmin[0], rmin[1], rmin[2]);
		glVertex3d(rmax[0], rmin[1], rmin[2]);
		glVertex3d(rmin[0], rmax[1], rmin[2]);
		glVertex3d(rmax[0], rmax[1], rmin[2]);
		glVertex3d(rmin[0], rmax[1], rmax[2]);
		glVertex3d(rmax[0], rmax[1], rmax[2]);
		glVertex3d(rmin[0], rmin[1], rmax[2]);
		glVertex3d(rmax[0], rmin[1], rmax[2]);
		glEnd();
	}
}

// precomputed array of cos and sin
const int angle_res = 12;
const double cos_array[] = {1.0, 0.86, 0.5,  0.0, -0.5 , -0.86, -1.0, -0.86, -0.5 ,  0.0,  0.5 ,  0.86, 1.0};
const double sin_array[] = {0.0, 0.5,  0.86, 1.0,  0.86,  0.5 ,  0.0, -0.5 , -0.86, -1.0, -0.86, -0.5 , 0.0};

void CanvasGL::Sphere(const vec3_t& center, real_t r){
	glLineWidth(lineWidth);
	glColor4fv(&lineColor[0].rgba[0]);

	/// three orthogonal circles
	glBegin(GL_LINE_LOOP);
	for(int i = 0; i < angle_res; i++){
		glVertex3d(r * cos_array[i], r * sin_array[i], 0.0);
	}
	glEnd();
	glBegin(GL_LINE_LOOP);
	for(int i = 0; i < angle_res; i++){
		glVertex3d(0.0, r * cos_array[i], r * sin_array[i]);
	}
	glEnd();
	glBegin(GL_LINE_LOOP);
	for(int i = 0; i < angle_res; i++){
		glVertex3d(r * sin_array[i], 0.0, r * cos_array[i]);
	}
	glEnd();
	
}

void CanvasGL::Cylinder(real_t r, real_t l){
	glLineWidth(lineWidth);
	glColor4fv(&lineColor[0].rgba[0]);

	real_t l_half = 0.5 * l;
	// caps
	glBegin(GL_LINE_LOOP);
	for(int i = 0; i < angle_res; i++){
		glVertex3d(r * cos_array[i], r * sin_array[i], l_half);
	}
	glEnd();
	
	glBegin(GL_LINE_LOOP);
	for(int i = 0; i < angle_res; i++){
		glVertex3d(r * cos_array[i], r * sin_array[i], -l_half);
	}
	glEnd();
	
	// sides
	glBegin(GL_LINES);
	glVertex3d( r  ,  0.0, l_half); glVertex3f( r  ,  0.0, -l_half);
	glVertex3d(-r  ,  0.0, l_half); glVertex3f(-r  ,  0.0, -l_half);
	glVertex3d( 0.0,  r  , l_half); glVertex3f( 0.0,  r  , -l_half);
	glVertex3d( 0.0, -r  , l_half); glVertex3f( 0.0, -r  , -l_half);
	glEnd();
	
}

}
}
