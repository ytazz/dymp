#pragma once

#include <dymp/canvas.h>

namespace dymp{;
namespace render{;

class CanvasGL : public Canvas{
public:
	virtual void Point          (const vec3_t& p);
	virtual void Line           (const vec3_t& p0, const vec3_t& p1);
	virtual void BeginPath      ();
	virtual void EndPath        ();
	virtual void MoveTo         (const vec3_t& p);
	virtual void LineTo         (const vec3_t& p);
	virtual void Rectangle      (const vec2_t& rmin, const vec2_t& rmax);
	virtual void Circle         (const vec2_t& center, real_t r);
	virtual void Box            (const vec3_t& rmin, const vec3_t& rmax);
	virtual void Sphere         (const vec3_t& center, real_t r);
	virtual void Cylinder       (real_t r, real_t l);
	
    CanvasGL();
};

}
}
