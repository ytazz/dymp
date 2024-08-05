#include <app.h>
#include <world.h>
#include <solver.h>
#include <model.h>
#include <variable.h>
#include <constraint.h>
#include <timer.h>

namespace dymp{;

Camera::Camera(){
    front = 0.01f;
    back  = 100.0f;
    aspect = 1.0f;
    fov    = 1.0f;
    width  = 1.0f;
    height = 1.0f;
    latitude  = 0.0f;
    longitude = 0.0f;
    distance  = 1.0f;
    
    target = Vector3f::Zero();

    rotMask  = ButtonState::Left;
    zoomMask = ButtonState::Right;
    trnMask  = ButtonState::Left | ButtonState::Alt;
    
    rotGain  = 0.01f;
    zoomGain = 0.01f;
    trnGain  = 1.0f;
    
    lonRange = Vector2f(deg_to_rad(-180.0), deg_to_rad(180.0));
    latRange = Vector2f(deg_to_rad(-90.0), deg_to_rad(90.0));
    distRange = Vector2f(0.01f, 100.0f);

    CalcTransform();
}

void Camera::CalcTransform(){
	affProj = Matrix4f::Zero();
	affProj(0,0) =  2.0f/width;
	affProj(1,1) =  2.0f/height;
	affProj(2,2) = -2.0f/(back - front);
	affProj(3,3) =  1.0f;
	affProj(0,3) =  0.0f;
	affProj(1,3) =  0.0f;
    affProj(2,3) = -(back + front) / (back - front);

    pos = distance * Eigen::Vector3f(
			cos(latitude) * sin(longitude),
			sin(latitude),
			cos(latitude) * cos(longitude));
    ori = Matrix3f(AngleAxisf(longitude, Vector3f::UnitY())*AngleAxisf(-latitude, Vector3f::UnitX()));

    affView.translation() = pos;
    affView.linear() = ori;

    affViewInv = affView.inverse();
}

void Camera::OnMouseMove(int button, int dx, int dy){
    if(button == rotMask){
		longitude -= (float)dx * rotGain;
		latitude  += (float)dy * rotGain;
		longitude = min(max(lonRange[0], longitude), lonRange[1]);
		latitude  = min(max(latRange[0], latitude ), latRange[1]);
	}
	if(button == zoomMask){
		distance *= (float)exp((double)dy * zoomGain);
	    distance = min(max(distRange[0], distance), distRange[1]);
	}
	if(button == trnMask){
		target += ori * Vector3f(-(float)dx * distance * trnGain, (float)dy * distance * trnGain, 0.0f);
	}
    CalcTransform();
}

App* App::instance = 0;

void App::KeyboardCallback   (GLFWwindow* window, int key, int scancode, int act, int mods){
    instance->OnKeyboard(window, key, scancode, act, mods);
}

void App::MouseButtonCallback(GLFWwindow* window, int button, int act, int mods){
    instance->OnMouseButton(window, button, act, mods);
}

void App::MouseMoveCallback  (GLFWwindow* window, double xpos, double ypos){
    instance->OnMouseMove(window, xpos, ypos);
}

void App::ScrollCallback     (GLFWwindow* window, double xoffset, double yoffset){
    instance->OnScroll(window, xoffset, yoffset);
}

App::App(){
    instance = this;

    button_left   = false;
    button_middle = false;
    button_right  = false;
    mod_shift     = false;
    mod_alt       = false;
    mod_ctrl      = false;

    lastx  = 0;
    lasty  = 0;

    running = false;
}

App::~App(){

}

void App::BuildScene(){

}

void App::OnStep(){

}

// keyboard callback
void App::OnKeyboard(GLFWwindow* window, int key, int scancode, int act, int mods) {
    /*
    // backspace: reset simulation
    if (act==GLFW_PRESS && key==GLFW_KEY_BACKSPACE) {

    }
    */
    if (act==GLFW_PRESS && key==GLFW_KEY_SPACE) {
        running = !running;
    }
}

// mouse button callback
void App::OnMouseButton(GLFWwindow* window, int button, int act, int mods) {
    // update button state
    button_left   = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT  ) == GLFW_PRESS);
    button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS);
    button_right  = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT ) == GLFW_PRESS);

    // update mouse position
    glfwGetCursorPos(window, &lastx, &lasty);
}

// mouse move callback
void App::OnMouseMove(GLFWwindow* window, double xpos, double ypos) {
    // no buttons down: nothing to do
    if (!button_left && !button_middle && !button_right) {
        return;
    }

    // compute mouse displacement, save
    double dx = xpos - lastx;
    double dy = ypos - lasty;
    lastx = xpos;
    lasty = ypos;

    // get current window size
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // get shift key state
    bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT ) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);

    int button = 
        (button_left   ? ButtonState::Left   : 0) | 
        (button_middle ? ButtonState::Middle : 0) | 
        (button_right  ? ButtonState::Right  : 0) | 
        (mod_shift     ? ButtonState::Shift  : 0) | 
        (mod_alt       ? ButtonState::Alt    : 0) | 
        (mod_ctrl      ? ButtonState::Ctrl   : 0);

    camera.OnMouseMove(button, dx, dy);

    /*
    // determine action based on mouse button
    mjtMouse action;
    if (button_right) {
        action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
    } 
    else if (button_left) {
        action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
    } 
    else {
        action = mjMOUSE_ZOOM;
    }
    */
}

// scroll callback
void App::OnScroll(GLFWwindow* window, double xoffset, double yoffset) {

}

bool App::Init(){
    if (!glfwInit()) {
        printf("Could not initialize GLFW\n");
        return false;
    }

    window = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // install GLFW mouse and keyboard callbacks
    glfwSetKeyCallback        (window, &App::KeyboardCallback   );
    glfwSetMouseButtonCallback(window, &App::MouseButtonCallback);
    glfwSetCursorPosCallback  (window, &App::MouseMoveCallback  );
    glfwSetScrollCallback     (window, &App::ScrollCallback     );

    world    = unique_ptr<World>(new World());
    conf     = 0;
	canvasGL = unique_ptr<render::CanvasGL>(new render::CanvasGL());
	
    BuildScene();

    return true;
}

void App::Loop(){
    // run main loop, target real-time simulation and 60 fps rendering
    while (!glfwWindowShouldClose(window)) {
        // get framebuffer viewport
        glfwGetFramebufferSize(window, &viewport_width, &viewport_height);

        if(running){
            world->Step();
            OnStep();
        }

        glClearColor(1.0f, 0.8f, 1.0f, 1.0f);
		glClearDepth(1.0); 
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	    
	    camera.fov    = 0.3f;
        camera.aspect = (float)viewport_height/(float)viewport_width;
	    camera.front  =   0.01f;
	    camera.back   = 100.0f;
	    camera.width  = 2.0f*camera.fov*camera.distance;
	    camera.height = camera.width*camera.aspect;
	    camera.CalcTransform();

        glMatrixMode(GL_PROJECTION);
	    glLoadMatrixf(camera.affProj.data());
	    glMatrixMode(GL_MODELVIEW);
        glLoadMatrixf(camera.affViewInv.data());
	    glRotatef(-90.0, 1.0f, 0.0f, 0.0f);
        
        world->Draw(canvasGL.get(), conf);

        // swap OpenGL buffers (blocking call due to v-sync)
        glfwSwapBuffers(window);

        // process pending GUI events, call GLFW callbacks
        glfwPollEvents();

        Timer::Sleep(10);
    }
}

void App::Cleanup(){
    // terminate GLFW (crashes with Linux NVidia drivers)
    glfwTerminate();
}

}