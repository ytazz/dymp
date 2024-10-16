#pragma once

#include <dymp/types.h>

namespace dymp{;

class Model;

namespace render{;

/// things to be drawn
struct Item{
	enum{
		GlobalAxis,         ///< global coordinate axis
		CentroidPos,
		CentroidEnd,
		CentroidEndTraj,
		CentroidFace,
		CentroidForce,
        CentroidTorso,
		WholebodyLink,
		WholebodyLinkTraj,
	};
};

class Canvas;

/**
	Configuration for Visualization
 */
class Config{
public:
	/** override this function to customize material settings
		@param render	rendering context
		@param id		identifier of an object
		@return			return true if the object specified by id is to be drawn, return false otherwise.
	 */
	virtual bool Set(Canvas* canvas, int attr, Model* model);

	virtual float Scale(int attr, Model* model);

	Config();
};

class Color{
public:
    static bool FromName(std::string name, Eigen::Vector4f& c);
	
public:
	std::string      name;
	Eigen::Vector4f  rgba;

public:
	void Init();

	Color(std::string n = "white"){
		name = n;
		Init();
	}
};

class Canvas{
public:
	Color   pointColor;
	Color   lineColor[2];		///< line color
	Color   fillColor;		    ///< fill color
	
	float   pointSize;			///< point size
	float   lineWidth;			///< line width
	size_t	nCircleDiv;			///< number of corners for drawing circle

	bool	drawLine;			///< draw outline?
	bool	drawFill;			///< fill inside shape?
	bool	transparent;		///< enable transparent?
	bool	gradation;			///< line gradation?

	void SetPointColor(const std::string& c, float alpha = 1.0f);
	void SetLineColor (const std::string& c, int idx = 0, float alpha = 1.0f);
	void SetFillColor (const std::string& c, float alpha = 1.0f);
	void SetPointSize (float ps);
	void SetLineWidth (float lw);
	
	virtual void Point          (const vec3_t& p){}
	virtual void Line           (const vec3_t& p0, const vec3_t& p1){}
	virtual void BeginPath      (){}
	virtual void EndPath        (){}
	virtual void MoveTo         (const vec3_t& p){}
	virtual void LineTo         (const vec3_t& p){}
	virtual void Rectangle      (const vec2_t& rmin, const vec2_t& rmax){}
	virtual void Circle         (const vec2_t& center, float r){}
	virtual void Box            (const vec3_t& rmin, const vec3_t& rmax){}
	virtual void Sphere         (const vec3_t& center, float r){}
	virtual void Cylinder       (float r, float l){}
	virtual void BeginLayer     (std::string name, bool shown){}
	virtual void EndLayer       (){}

	Canvas();
};

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
