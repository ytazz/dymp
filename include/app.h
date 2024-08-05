#pragma once

#include <types.h>
#include <world.h>

#include <GLFW/glfw3.h>

namespace dymp{;

struct ButtonState{
    enum{
        Left   = 0x01 << 0,
        Middle = 0x01 << 1,
        Right  = 0x01 << 2,
        Shift  = 0x01 << 3,
        Alt    = 0x01 << 4,
        Ctrl   = 0x01 << 5
    };
};

class Camera{
public:
    float front, back;
    float aspect, fov;
    float width, height;
    float latitude, longitude, distance;
    Eigen::Vector3f  target;

    int    rotMask, zoomMask, trnMask;
    float  rotGain, zoomGain, trnGain;
    Eigen::Vector2f   lonRange, latRange, distRange;

    Eigen::Vector3f   pos;
    Eigen::Matrix3f   ori;
    Eigen::Matrix4f   affProj;
    Eigen::Transform<float, 3, Eigen::Isometry> affView, affViewInv;

public:
    void OnMouseMove(int button, int dx, int dy);
    void CalcTransform();

    Camera();
};

class App{
public:
    static App*  instance;

    static void KeyboardCallback   (GLFWwindow* window, int key, int scancode, int act, int mods);
    static void MouseButtonCallback(GLFWwindow* window, int button, int act, int mods);
    static void MouseMoveCallback  (GLFWwindow* window, double xpos, double ypos);
    static void ScrollCallback     (GLFWwindow* window, double xoffset, double yoffset);


public:
    GLFWwindow*  window;
    std::unique_ptr<World>   world;
    render::Config*          conf;
	std::unique_ptr<render::CanvasGL>  canvasGL;
	
    int    viewport_width, viewport_height;
    bool   button_left, button_middle, button_right;
    bool   mod_shift, mod_alt, mod_ctrl;
    double lastx, lasty;   //< previous mouse position

    Camera  camera;

public:
    virtual void BuildScene();
    virtual void OnStep();
    virtual void OnKeyboard(GLFWwindow* window, int key, int scancode, int act, int mods);
    virtual void OnMouseButton(GLFWwindow* window, int button, int act, int mods);
    virtual void OnMouseMove(GLFWwindow* window, double xpos, double ypos);
    virtual void OnScroll(GLFWwindow* window, double xoffset, double yoffset);

    bool Init();
    void Loop();
    void Cleanup();

    App();
    virtual ~App();
};

}