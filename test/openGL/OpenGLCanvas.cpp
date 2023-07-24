#include <vector>

#include "OpenGLCanvas.h"
#include "OpenGLRenderer.h"
#include "../meshdata.h"
#include "../argparser.h"
#include "../camera.h"
#include "../mesh.h"
#include "../raytrace/raytracer.h"

// ========================================================
// static variables of OpenGLCanvas class

ArgParser* OpenGLCanvas::args = NULL;
MeshData* OpenGLCanvas::mesh_data = NULL;
OpenGLRenderer* OpenGLCanvas::renderer = NULL;
GLFWwindow* OpenGLCanvas::window = NULL;

// mouse position
int OpenGLCanvas::mouseX = 0;
int OpenGLCanvas::mouseY = 0;
// which mouse button
bool OpenGLCanvas::keys[1024] = { false };
bool OpenGLCanvas::leftMousePressed = false;
bool OpenGLCanvas::middleMousePressed = false;
bool OpenGLCanvas::rightMousePressed = false;
// current state of modifier keys
bool OpenGLCanvas::shiftKeyPressed = false;
bool OpenGLCanvas::controlKeyPressed = false;
bool OpenGLCanvas::altKeyPressed = false;
bool OpenGLCanvas::superKeyPressed = false;

// ========================================================
// Initialize all appropriate OpenGL variables, set
// callback functions, and start the main event loop.
// This function will not return but can be terminated
// by calling 'exit(0)'
// ========================================================

void OpenGLCanvas::initialize(ArgParser *_args, MeshData *_mesh_data, OpenGLRenderer *_renderer) {
    args = _args;
    mesh_data = _mesh_data;
    renderer = _renderer;

    glfwSetErrorCallback(error_callback);

    // Initialize GLFW
    if( !glfwInit() ) {
        std::cerr << "ERROR: Failed to initialize GLFW" << std::endl;
        exit(1);
    }

    // We will ask it to specifically open an OpenGL 3.2 context
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create a GLFW window
    window = glfwCreateWindow(mesh_data->width,mesh_data->height, "ACG HW1 MESHES", NULL, NULL);
    if (!window) {
        std::cerr << "ERROR: Failed to open GLFW window" << std::endl;
        glfwTerminate();
        exit(1);
    }
    glfwMakeContextCurrent(window);
    HandleGLError("in glcanvas first");

    // Initialize GLEW
    glewExperimental = true; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        std::cerr << "ERROR: Failed to initialize GLEW" << std::endl;
        glfwTerminate();
        exit(1);
    }

    // there seems to be a "GL_INVALID_ENUM" error in glewInit that is a
    // know issue, but can safely be ignored
    HandleGLError("after glewInit()",true);

    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "OpenGL Version: " << (char*)glGetString(GL_VERSION) << '\n';
    std::cout << "-------------------------------------------------------" << std::endl;

    // Initialize callback functions
    glfwSetCursorPosCallback(OpenGLCanvas::window,OpenGLCanvas::mousemotionCB);
    glfwSetMouseButtonCallback(OpenGLCanvas::window,OpenGLCanvas::mousebuttonCB);
    glfwSetKeyCallback(OpenGLCanvas::window,OpenGLCanvas::keyboardCB);

    GLOBAL_args->camera->glPlaceCamera();

    HandleGLError("finished glcanvas initialize");
}

// ========================================================
// Callback function for mouse click or release
// ========================================================

void OpenGLCanvas::mousebuttonCB(GLFWwindow* /*window*/, int which_button, int action, int /*mods*/) {
    // store the current state of the mouse buttons
    if (which_button == GLFW_MOUSE_BUTTON_1) {
        if (action == GLFW_PRESS) {
            leftMousePressed = true;
        } else {
            assert (action == GLFW_RELEASE);
            leftMousePressed = false;
        }
    } else if (which_button == GLFW_MOUSE_BUTTON_2) {
        if (action == GLFW_PRESS) {
            rightMousePressed = true;
        } else {
            assert (action == GLFW_RELEASE);
            rightMousePressed = false;
        }
    } else if (which_button == GLFW_MOUSE_BUTTON_3) {
        if (action == GLFW_PRESS) {
            middleMousePressed = true;
        } else {
            assert (action == GLFW_RELEASE);
            middleMousePressed = false;
        }
    }
}	

// ========================================================
// Callback function for mouse drag
// ========================================================

void OpenGLCanvas::mousemotionCB(GLFWwindow* /*window*/, double x, double y) {
    // camera controls that work well for a 3 button mouse
    if (!shiftKeyPressed && !controlKeyPressed && !altKeyPressed) {
        if (leftMousePressed) {
            GLOBAL_args->camera->rotateCamera(mouseX-x,mouseY-y);
        } else if (middleMousePressed)  {
            GLOBAL_args->camera->truckCamera(mouseX-x, y-mouseY);
        } else if (rightMousePressed) {
            GLOBAL_args->camera->dollyCamera(mouseY-y);
        }
    }
    if (leftMousePressed || middleMousePressed || rightMousePressed) {
        if (shiftKeyPressed) {
            GLOBAL_args->camera->zoomCamera(mouseY-y);
        }
        // allow reasonable control for a non-3 button mouse
        if (controlKeyPressed) {
            GLOBAL_args->camera->truckCamera(mouseX-x, y-mouseY);    
        }
        if (altKeyPressed) {
            GLOBAL_args->camera->dollyCamera(y-mouseY);    
        }
    }
    mouseX = x;
    mouseY = y;
    mesh_data->raytracing_animation = false;
}

// ========================================================
// Callback function for keyboard events
// ========================================================


// NOTE: These functions are also called by the Mac Metal Objective-C
// code, so we need this extern to allow C code to call C++ functions
// (without function name mangling confusion).

extern "C" {
void Load();
void RayTreeActivate();
void RayTreeDeactivate();
void PhotonMappingTracePhotons();
void RaytracerClear();
void PhotonMappingClear();
void PackMesh();
void Step();
void Animate_Fluid();
}

void OpenGLCanvas::keyboardCB(GLFWwindow* /*window*/, int key, int /*scancode*/, int action, int mods) {
    // store the modifier keys
    shiftKeyPressed = (GLFW_MOD_SHIFT & mods);
    controlKeyPressed = (GLFW_MOD_CONTROL & mods);
    altKeyPressed = (GLFW_MOD_ALT & mods);
    superKeyPressed = (GLFW_MOD_SUPER & mods);
    // non modifier key actions
    if (key == GLFW_KEY_ESCAPE || key == 'q' || key == 'Q') {
        glfwSetWindowShouldClose(OpenGLCanvas::window, GL_TRUE);
        // force quit
        exit(0);
    }
    else if (keys['L']) {
        if (GLOBAL_args->camera->lock) {
            GLOBAL_args->camera->lock = !GLOBAL_args->camera->lock;
            std::cout << "unlock focus point" << std::endl;
        }
        else {
            GLOBAL_args->camera->lock = !GLOBAL_args->camera->lock;
            GLOBAL_args->camera->focus = GLOBAL_args->camera->point_of_interest;
            std::cout << "lock focus point" << std::endl;
        }
    }else if (keys['P']) {
        if (GLOBAL_args->camera->first_person) {
            GLOBAL_args->camera->first_person = !GLOBAL_args->camera->first_person;
            GLOBAL_args->camera->focus = GLOBAL_args->mesh->getPersonVertex(5)->get();
            std::cout << "start first person" << std::endl;
        }else {
            GLOBAL_args->camera->first_person = !GLOBAL_args->camera->first_person;
            GLOBAL_args->camera->focus = GLOBAL_args->camera->point_of_interest;
            std::cout << "stop first person" << std::endl;
        }
    }
    if (action == GLFW_PRESS) {
       keys[key] = true;
    }else if (action == GLFW_RELEASE)
       keys[key] = false;
}

void OpenGLCanvas::move_detect() {
    Vec3f pos = GLOBAL_args->camera->getDirection();
    double temp = pos.z();
    Vec3f dir;
    dir = pos;
    dir.setz(pos.x());
    dir.setx(-temp);
    if (keys['W']) {
        if (!GLOBAL_args->camera->lock)
            GLOBAL_args->camera->moveCamera("forward");
        else{
            if (keys['S'])
                return;
            GLOBAL_args->mesh->move_person(pos);
            if (keys['D']) {
                GLOBAL_args->mesh->move_person(dir);
            }
            else if (keys['A']) {
                GLOBAL_args->mesh->move_person(-dir);
            }
        }
        PackMesh();
    }else if (keys['S']) {
        if (!GLOBAL_args->camera->lock)
            GLOBAL_args->camera->moveCamera("back");
        else {
            if (keys['W'])
                return;
            GLOBAL_args->mesh->move_person(-pos);
            if (keys['D']) {
                GLOBAL_args->mesh->move_person(dir);
            }else if (keys['A']) {
                GLOBAL_args->mesh->move_person(-dir);
            }
        }
        PackMesh();
    }else if (keys['D']) {
        if (!GLOBAL_args->camera->lock)
            GLOBAL_args->camera->moveCamera("right");
        else {
            if (keys['A'])
                return;
            GLOBAL_args->mesh->move_person(dir);
        }
        PackMesh();
    }else if (keys['A']) {
        if (!GLOBAL_args->camera->lock)
            GLOBAL_args->camera->moveCamera("left");
        else {
            if (keys['D'])
                return;
            GLOBAL_args->mesh->move_person(-dir);
            PackMesh();
        }
    }
}

void OpenGLCanvas::gravity() {
    //GLOBAL_args->mesh->move_person("down");
    
}
// ========================================================
// Load the vertex & fragment shaders
// ========================================================

GLuint LoadShaders(const std::string &vertex_file_path,const std::string &fragment_file_path){

    // Create the shaders
    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

    // Read the Vertex Shader code from the file
    std::string VertexShaderCode;
    std::ifstream VertexShaderStream(vertex_file_path.c_str(), std::ios::in);
    if (VertexShaderStream.is_open()){
        std::string Line = "";
        while(getline(VertexShaderStream, Line))
            VertexShaderCode += "\n" + Line;
        VertexShaderStream.close();
    } else {
        std::cerr << "ERROR: cannot open " << vertex_file_path << std::endl;
        exit(0);
    }
    // Read the Fragment Shader code from the file
    std::string FragmentShaderCode;
    std::ifstream FragmentShaderStream(fragment_file_path.c_str(), std::ios::in);
    if(FragmentShaderStream.is_open()){
        std::string Line = "";
        while(getline(FragmentShaderStream, Line)) {
            FragmentShaderCode += "\n" + Line;
        }
        FragmentShaderStream.close();
    } else {
        std::cerr << "ERROR: cannot open " << fragment_file_path << std::endl;
        exit(0);
    }

    GLint Result = GL_FALSE;
    int InfoLogLength;

    // Compile Vertex Shader
    std::cout << "Compiling shader : " << vertex_file_path << std::endl;
    char const * VertexSourcePointer = VertexShaderCode.c_str();
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
    glCompileShader(VertexShaderID);
    // Check Vertex Shader
    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
        std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
        glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
        //std::cerr << "ERROR: " << VertexShaderErrorMessage[0] << std::endl;
    }

    // Compile Fragment Shader
    std::cout << "Compiling shader : " << fragment_file_path << std::endl;
    char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
    glCompileShader(FragmentShaderID);
    // Check Fragment Shader
    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
        std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
        glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
        //std::cerr << "ERROR: " << FragmentShaderErrorMessage[0] << std::endl;
    }

    // Link the program
    std::cout << "Linking program" << std::endl;
    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, VertexShaderID);
    glAttachShader(ProgramID, FragmentShaderID);
    glLinkProgram(ProgramID);
    // Check the program
    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
        std::vector<char> ProgramErrorMessage(InfoLogLength+1);
        glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
        //std::cerr << "ERROR: " << ProgramErrorMessage[0] << std::endl;
    }

    glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);

    return ProgramID;
}

// ========================================================
// Functions related to error handling
// ========================================================

void OpenGLCanvas::error_callback(int /*error*/, const char* description) {
    std::cerr << "ERROR CALLBACK: " << description << std::endl;
}

std::string WhichGLError(GLenum &error) {
    switch (error) {
        case GL_NO_ERROR:
            return "NO_ERROR";
        case GL_INVALID_ENUM:
            return "GL_INVALID_ENUM";
        case GL_INVALID_VALUE:
            return "GL_INVALID_VALUE";
        case GL_INVALID_OPERATION:
            return "GL_INVALID_OPERATION";
        case GL_INVALID_FRAMEBUFFER_OPERATION:
            return "GL_INVALID_FRAMEBUFFER_OPERATION";
        case GL_OUT_OF_MEMORY:
            return "GL_OUT_OF_MEMORY";
        case GL_STACK_UNDERFLOW:
            return "GL_STACK_UNDERFLOW";
        case GL_STACK_OVERFLOW:
            return "GL_STACK_OVERFLOW";
        default:
            return "OTHER GL ERROR";
    }
}

int HandleGLError(const std::string &message, bool ignore) {
    GLenum error;
    int i = 0;
    while ((error = glGetError()) != GL_NO_ERROR) {
        if (!ignore) {
            if (message != "") {
                std::cout << "[" << message << "] ";
            }
            std::cout << "GL ERROR(" << i << ") " << WhichGLError(error) << std::endl;
        }
        i++;
    }
    if (i == 0) return 1;
    return 0;
}



// ====================================================================
// Construct the ViewMatrix & ProjectionMatrix for GL Rendering
// ====================================================================

void convert(Matrix &b, const glm::mat4 &a);

void OrthographicCamera::glPlaceCamera() {

    if (OpenGLCanvas::window == NULL) return;

    glfwGetWindowSize(OpenGLCanvas::window, &GLOBAL_args->mesh_data->width, &GLOBAL_args->mesh_data->height);
    float aspect = GLOBAL_args->mesh_data->width / (float)GLOBAL_args->mesh_data->height;
    float w;
    float h;
    // handle non square windows
    if (aspect < 1.0) {
    w = size / 2.0;
        h = w / aspect;
    } else {
        h = size / 2.0;
        w = h * aspect;
    }
    
    glm::mat4 projmat = glm::ortho<float>(-w,w,-h,h, 0.1f, 100.0f) ;
    glm::vec3 campos(camera_position.x(),camera_position.y(), camera_position.z());
    glm::vec3 poi(point_of_interest.x(),point_of_interest.y(),point_of_interest.z());
    glm::vec3 screenup(getScreenUp().x(),getScreenUp().y(),getScreenUp().z());
    glm::mat4 viewmat = glm::lookAt(campos,poi,screenup);

    convert(ProjectionMatrix,projmat);
    convert(ViewMatrix,viewmat);
}

void PerspectiveCamera::glPlaceCamera() {

    if (OpenGLCanvas::window == NULL) return;

    glfwGetWindowSize(OpenGLCanvas::window, &GLOBAL_args->mesh_data->width, &GLOBAL_args->mesh_data->height);
    float aspect = GLOBAL_args->mesh_data->width / (float)GLOBAL_args->mesh_data->height;
    // must convert angle from degrees to radians
 
    Vec3f point;

    if(lock){
        point = point_of_interest;
    }
    else{
        point = focus; 
    }

    //std::cout<<camera_position.x()<<","<<camera_position.y()<<","<<camera_position.z()<<std::endl;


    glm::vec3 campos(camera_position.x(),camera_position.y(), camera_position.z());

    glm::vec3 poi(point.x(),point.y(),point.z());
    glm::vec3 screenup(getScreenUp().x(),getScreenUp().y(),getScreenUp().z());

    glm::mat4 projmat = glm::perspective<float>(glm::radians(angle), aspect, 0.1f, 1000.0f);
    convert (ProjectionMatrix,projmat);
    glm::mat4 viewmat =  glm::lookAt(campos,poi,screenup);
    convert (ViewMatrix,viewmat);
}

// ========================================================
// ========================================================
