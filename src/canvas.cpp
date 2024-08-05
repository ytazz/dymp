#include <canvas.h>
#include <util.h>
#include <model.h>

#include <GL/glut.h>

using namespace std;
using namespace Eigen;

namespace dimp3{;
namespace render{;

Config::Config(){
}

bool Config::Set(Canvas* canvas, int attr, Model* model){
	string c;
	float  lw = 1.0f;
	float  ps = 1.0f;
	
	     if(attr == Item::CentroidPos     ){ c = "black"  ; lw = 3.0f; }
	else if(attr == Item::CentroidEnd     ){ c = "black"  ; lw = 0.5f; }
	else if(attr == Item::CentroidEndTraj ){ c = "black"  ; lw = 1.5f; }
	else if(attr == Item::CentroidFace    ){ c = "green"  ; lw = 1.0f; }
	else if(attr == Item::CentroidForce   ){ c = "magenta"; lw = 1.0f; }
    else if(attr == Item::CentroidTorso   ){ c = "black"  ; lw = 1.0f; }

	canvas->SetLineColor(c, 0);
	canvas->SetLineColor(c, 1);
	canvas->SetFillColor(c);
	canvas->SetPointColor(c);
	canvas->SetPointSize(ps);
	canvas->SetLineWidth(lw);
	
	return true;
}

float Config::Scale(int attr, Model* model){
	return 1.0f;
}

struct ColorInfo{
	const char* name;
	Eigen::Vector4f rgba;
};

bool ColorFromName(string name, Eigen::Vector4f& c){
	#undef  RGB
	#define RGB(x, y, z) Eigen::Vector4f((float)x/255.0f, (float)y/255.0f, (float)z/255.0f, 1.0f)
	static const ColorInfo colorInfo[] = {
		{"indianred"           , RGB(205,  92,  92)},
		{"lightcoral"          , RGB(240, 128, 128)},
		{"salmon"              , RGB(250, 128, 114)},
		{"darksalmon"          , RGB(233, 150, 122)},
		{"lightsalmon"         , RGB(255, 160, 122)},
		{"red"                 , RGB(255,   0,   0)},
		{"crimson"             , RGB(220,  20,  60)},
		{"firebrick"           , RGB(178,  34,  34)},
		{"darkred"             , RGB(139,   0,   0)},
		{"pink"                , RGB(255, 192, 203)},
		{"lightpink"           , RGB(255, 182, 193)},
		{"hotpink"             , RGB(255, 105, 180)},
		{"deeppink"            , RGB(255,  20, 147)},
		{"mudiumvioletred"     , RGB(199,  21, 133)},
		{"palevioletred"       , RGB(219, 112, 147)},
		{"coral"               , RGB(255, 127,  80)},
		{"tomato"              , RGB(255,  99,  71)},
		{"orangered"           , RGB(255,  69,   0)},
		{"darkorange"          , RGB(255, 140,   0)},
		{"orange"              , RGB(255, 165,   0)},
		{"gold"                , RGB(255, 215,   0)},
		{"yellow"              , RGB(255, 255,   0)},
		{"lightyellow"         , RGB(255, 255, 224)},
		{"lemonchiffon"        , RGB(255, 250, 205)},
		{"lightgoldenrodyellow", RGB(250, 250, 210)},
		{"papayawhip"          , RGB(255, 239, 213)},
		{"moccasin"            , RGB(255, 228, 181)},
		{"peachpuff"           , RGB(255, 218, 185)},
		{"palegoldenrod"       , RGB(238, 232, 170)},
		{"khaki"               , RGB(240, 230, 140)},
		{"darkkhaki"           , RGB(189, 183, 107)},
		{"lavender"            , RGB(230, 230, 250)},
		{"thistle"             , RGB(216, 191, 216)},
		{"plum"                , RGB(221, 160, 221)},
		{"violet"              , RGB(238, 130, 238)},
		{"orchild"             , RGB(218, 112, 214)},
		{"fuchsia"             , RGB(255,   0, 255)},
		{"magenta"             , RGB(255,   0, 255)},
		{"mediumorchild"       , RGB(186,  85, 211)},
		{"mediumpurple"        , RGB(147, 112, 219)},
		{"blueviolet"          , RGB(138,  43, 226)},
		{"darkviolet"          , RGB(148,   0, 211)},
		{"darkorchild"         , RGB(153,  50, 204)},
		{"darkmagenta"         , RGB(139,   0, 139)},
		{"purple"              , RGB(128,   0, 128)},
		{"indigo"              , RGB( 75,   0, 130)},
		{"darkslateblue"       , RGB( 72,  61, 139)},
		{"slateblue"           , RGB(106,  90, 205)},
		{"mediumslateblue"     , RGB(123, 104, 238)},
		{"greenyellow"         , RGB(173, 255,  47)},
		{"chartreuse"          , RGB(127, 255,   0)},
		{"lawngreen"           , RGB(124, 252,   0)},
		{"lime"                , RGB(  0, 252,   0)},
		{"limegreen"           , RGB( 50, 205,  50)},
		{"palegreen"           , RGB(152, 251, 152)},
		{"lightgreen"          , RGB(144, 238, 144)},
		{"mediumspringgreen"   , RGB(  0, 250, 154)},
		{"springgreen"         , RGB(  0, 255, 127)},
		{"mediumseagreen"      , RGB( 60, 179, 113)},
		{"seagreen"            , RGB( 46, 139,  87)},
		{"forestgreen"         , RGB( 34, 139,  34)},
		{"green"               , RGB(  0, 128,   0)},
		{"darkgreen"           , RGB(  0, 100,   0)},
		{"yellowgreen"         , RGB(154, 205,  50)},
		{"olivedrab"           , RGB(107, 142,  35)},
		{"olive"               , RGB(128, 128,   0)},
		{"darkolivegreen"      , RGB( 85, 107,  47)},
		{"mediumaquamarine"    , RGB(102, 205, 170)},
		{"darkseagreen"        , RGB(143, 188, 143)},
		{"lightseagreen"       , RGB( 32, 178, 170)},
		{"darkcyan"            , RGB(  0, 139, 139)},
		{"teal"                , RGB(  0, 128, 128)},
		{"aqua"                , RGB(  0, 255, 255)},
		{"cyan"                , RGB(  0, 255, 255)},
		{"lightcyan"           , RGB(224, 255, 255)},
		{"paleturquoise"       , RGB(175, 238, 238)},
		{"aquamarine"          , RGB(127, 255, 212)},
		{"turquoise"           , RGB( 64, 224, 208)},
		{"mediumturquoise"     , RGB( 72, 209, 204)},
		{"darkturquoise"       , RGB(  0, 206, 209)},
		{"cadetblue"           , RGB( 95, 158, 160)},
		{"steelblue"           , RGB( 70, 130, 180)},
		{"lightsteelblue"      , RGB(176, 196, 222)},
		{"powderblue"          , RGB(176, 224, 230)},
		{"lightblue"           , RGB(173, 216, 230)},
		{"skyblue"             , RGB(135, 206, 235)},
		{"lightskyblue"        , RGB(135, 206, 250)},
		{"deepskyblue"         , RGB(  0, 191, 255)},
		{"dodgerblue"          , RGB( 30, 144, 255)},
		{"cornflowerblue"      , RGB(100, 149, 237)},
		{"royalblue"           , RGB( 65, 105, 225)},
		{"blue"                , RGB(  0,   0, 255)},
		{"mediumblue"          , RGB(  0,   0, 205)},
		{"darkblue"            , RGB(  0,   0, 139)},
		{"navy"                , RGB(  0,   0, 128)},
		{"midnightblue"        , RGB( 25,  25, 112)},
		{"cornsilk"            , RGB(255, 248, 220)},
		{"blanchedalmond"      , RGB(255, 235, 205)},
		{"bisque"              , RGB(255, 228, 196)},
		{"navajowhite"         , RGB(255, 222, 173)},
		{"wheat"               , RGB(245, 222, 179)},
		{"burlywood"           , RGB(222, 184, 135)},
		{"tan"                 , RGB(210, 180, 140)},
		{"rosybrown"           , RGB(188, 143, 143)},
		{"sandybrown"          , RGB(244, 164,  96)},
		{"goldenrod"           , RGB(218, 165,  32)},
		{"darkgoldenrod"       , RGB(184, 134,  11)},
		{"peru"                , RGB(205, 133,  63)},
		{"chocolate"           , RGB(210, 105,  30)},
		{"saddlebrown"         , RGB(139,  69,  19)},
		{"sienna"              , RGB(160,  82,  45)},
		{"brown"               , RGB(165,  42,  42)},
		{"maroon"              , RGB(128,   0,   0)},
		{"white"               , RGB(255, 255, 255)},
		{"snow"                , RGB(255, 250, 250)},
		{"honeydew"            , RGB(240, 255, 240)},
		{"mintcream"           , RGB(245, 255, 250)},
		{"azure"               , RGB(240, 255, 255)},
		{"aliceblue"           , RGB(240, 248, 255)},
		{"ghostwhite"          , RGB(248, 248, 255)},
		{"whitesmoke"          , RGB(245, 245, 245)},
		{"seashell"            , RGB(255, 245, 238)},
		{"beige"               , RGB(245, 245, 220)},
		{"oldlace"             , RGB(253, 245, 230)},
		{"floralwhite"         , RGB(255, 250, 240)},
		{"ivory"               , RGB(255, 255, 240)},
		{"antiquewhite"        , RGB(250, 235, 215)},
		{"linen"               , RGB(250, 240, 230)},
		{"lavenderblush"       , RGB(255, 240, 245)},
		{"mistyrose"           , RGB(255, 228, 225)},
		{"gainsboro"           , RGB(220, 220, 220)},
		{"lightgray"           , RGB(211, 211, 211)},
		{"silver"              , RGB(192, 192, 192)},
		{"darkgray"            , RGB(169, 169, 169)},
		{"gray"                , RGB(128, 128, 128)},
		{"dimgray"             , RGB(105, 105, 105)},
		{"lightslategray"      , RGB(119, 136, 153)},
		{"slategray"           , RGB(112, 128, 144)},
		{"darkslategray"       , RGB( 47,  79,  79)},
		{"black"               , RGB(  0,   0,   0)}
	};
	#undef RGB
	
    if(name == "none"){
		c = Eigen::Vector4f();
		return true;
	}

	const int n = sizeof(colorInfo)/sizeof(ColorInfo);
	for(int i = 0; i < n; i++){
		if(name == colorInfo[i].name){
			c = colorInfo[i].rgba;
			return true;
		}
	}
	return false;
}

void Color::Init(){
	ColorFromName(name, rgba);
}

//-------------------------------------------------------------------------------------------------
// Canvas

Canvas::Canvas(){
	SetPointColor("black");
	SetLineColor ("black", 0);
	SetLineColor ("black", 1);
	SetFillColor ("white");
	SetPointSize (3.0f);
	SetLineWidth (1.0f);
	nCircleDiv	= 18;
	drawFill	= false;
	drawLine	= true;
	gradation	= false;
}
void Canvas::SetPointColor(const string& c, float alpha){
	if(pointColor.name != c){
		pointColor.name = c;
		pointColor.Init();
	}
	pointColor.rgba[3] = alpha;
}
void Canvas::SetLineColor(const string& c, int idx, float alpha){
	if(lineColor[idx].name != c){
		lineColor[idx].name = c;
		lineColor[idx].Init();
	}
	lineColor[idx].rgba[3] = alpha;
}
void Canvas::SetFillColor(const string& c, float alpha){
	if(fillColor.name != c){
		fillColor.name = c;
		fillColor.Init();
	}
	fillColor.rgba[3] = alpha;
}
void Canvas::SetPointSize(float ps){
	pointSize = ps;
}
void Canvas::SetLineWidth(float lw){
	lineWidth = lw;
}

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
